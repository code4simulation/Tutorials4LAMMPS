#!/usr/bin/env C:/Users/user/Downloads/dev_w_antigravity/pf_clustering/TI_lammps/venv_TI/Scripts/python.exe
"""
TI Workflow: Solid-phase free energy via hybrid/scaled Frenkel-Ladd + Reversible Scaling.

Implements the methodology in Sum.md §3 + §4 (hybrid/scaled calphy approach).

Compliant with Sum.md milestones:
  M1 - LAMMPS build (EXTRA-PAIR) + sandbox test:   hybrid/scaled verified with Fe.eam
  M2 - Solid-phase F(T):
        MSD-based spring constant per species     [Sum.md §5.2 milestone 2-1, 3.2]
        COM momentum fix                          [Sum.md §3.1]
        Hysteresis (eV/atom) check + plot         [Sum.md §5.2 milestone 2-4, checklist]
  M3 - Extension to liquid / multi-phase:           not yet (future work)

Usage:
    python run_ti.py                          # run full pipeline
    python run_ti.py --phase bcc              # single phase
    python run_ti.py --step fl --phase fcc    # only FL for FCC
    python run_ti.py --step analysis          # only analysis + plots

Config is in the SystemConfig dataclass below.  To change the system, edit that block.
"""

from __future__ import annotations
import os, sys, subprocess, argparse, textwrap, logging, yaml
from typing import List, Dict, Optional
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(), logging.FileHandler("ti_workflow.log", mode="w")],
)
log = logging.getLogger("TI")

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants  (from shared module)
# ─────────────────────────────────────────────────────────────────────────────
from ti_common import KB, HBAR_SI, KB_SI, H_SI, EV_TO_J, AMU_KG
from ti_common import run_lammps as _run_lammps_common


# ─────────────────────────────────────────────────────────────────────────────
# System Configuration
# ─────────────────────────────────────────────────────────────────────────────
class PotentialConfig:
    def __init__(self, style, coeff):
        self.style = style
        self.coeff = coeff

class PhaseConfig:
    def __init__(self, name, structure_dict, supercell=None, T0=None, T_min=None, T_max=None):
        self.name = name
        self.structure = structure_dict or {}
        self.supercell = supercell or [1, 1, 1]
        self.T0 = T0
        self.T_min = T_min
        self.T_max = T_max
        self.n_atoms = None
        self.vol_per_atom = None

class SystemConfig:
    def __init__(self, path="config.yaml"):
        with open(path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
            
        sys_cfg = cfg.get('system', {})
        self.lmp_exe = sys_cfg.get('lmp_exe', "lmp")
        self.python_exe = sys_cfg.get('python_exe', sys.executable)
        self.workdir = sys_cfg.get('workdir', ".")
        
        p = cfg.get('parameters', {})
        self.T0 = float(p.get('T0', 1000.0))
        self.T_min = float(p.get('T_min', 1000.0))
        self.T_max = float(p.get('T_max', 1000.0))
        self.pressure = float(p.get('pressure', 0.0))
        self.N_eq = int(p.get('N_eq', 20000))
        self.N_fl = int(p.get('N_fl', 50000))
        self.N_rs = int(p.get('N_rs', 80000))
        self.hysteresis_warn = float(p.get('hysteresis_warn', 0.05))
        self.rs_scheduling = p.get('rs_scheduling', 'linear-T')
        self.msd_per_element = p.get('msd_per_element', True)
        self.elements = cfg.get('elements', [])
        if not self.elements and 'element' in cfg:
            self.elements = [cfg['element']]
        self.element = self.elements[0] if self.elements else "Si"
        self.phases = []
        for ph in cfg.get('phases', []):
            self.phases.append(PhaseConfig(
                ph['name'], 
                ph.get('structure', {}),
                supercell=ph.get('supercell', [1, 1, 1]),
                T0=float(ph['T0']) if 'T0' in ph else None,
                T_min=float(ph['T_min']) if 'T_min' in ph else None,
                T_max=float(ph['T_max']) if 'T_max' in ph else None,
            ))
            
        self.potentials = []
        for pot in cfg.get('potentials', []):
            self.potentials.append(PotentialConfig(pot['style'], pot['coeff']))
            
        # Species information: { type_id (1-based): {"symbol": str, "mass": float} }
        self.species = {} 

    def phase(self, name: str) -> PhaseConfig:
        for p in self.phases:
            if p.name == name:
                return p
        raise KeyError(f"Unknown phase '{name}'")

    def data_file(self, phase_name: str) -> str:
        return os.path.join(self.workdir, f"ti_{phase_name}.data")

    def conf_eq_file(self, phase_name: str) -> str:
        return os.path.join(self.workdir, f"conf_eq_{phase_name}.data")
        
    def get_pair_commands(self, indent=12) -> str:
        spaces = " " * indent
        if len(self.potentials) > 1:
            styles = " ".join([p.style for p in self.potentials])
            out = f"pair_style hybrid/scaled v_lambda hybrid/overlay {styles}\n"
        else:
            style = self.potentials[0].style
            out = f"pair_style hybrid/scaled v_lambda {style}\n"
            
        for p in self.potentials:
            out += f"{spaces}pair_coeff {p.coeff}\n"
        return out.strip()

    def get_mass_commands(self) -> str:
        return "\n".join(f"mass {tid} {info['mass']}" for tid, info in self.species.items())


# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS runner  (thin wrapper around ti_common)
# ─────────────────────────────────────────────────────────────────────────────
def run_lmp(cfg: SystemConfig, script_text: str, script_name: str):
    _run_lammps_common(cfg.lmp_exe, script_text, script_name,
                       workdir=cfg.workdir, logger=log)


# ─────────────────────────────────────────────────────────────────────────────
# Step 0: Prepare structures with ASE
# ─────────────────────────────────────────────────────────────────────────────
def prepare_structures(cfg: SystemConfig):
    from ase.io import read, write
    from ase.data import atomic_masses, atomic_numbers
    os.makedirs(cfg.workdir, exist_ok=True)

    for ph in cfg.phases:
        out = cfg.data_file(ph.name)
        if os.path.exists(out):
            log.info(f"[SKIP] {out} already exists, loading metadata only")
            atoms = read(out, format="lammps-data", style="atomic")
            ph.n_atoms = len(atoms)
            ph.vol_per_atom = atoms.get_volume() / len(atoms)
            
            # Count atoms per type
            types = atoms.get_array('type') if atoms.has('type') else np.ones(len(atoms), dtype=int)
            unique_types = np.unique(types)
            ph.n_per_type = {int(t): int(np.sum(types == t)) for t in unique_types}
            
            from ase.data import atomic_masses, atomic_numbers
            for i, tid in enumerate(unique_types):
                element = cfg.elements[i] if i < len(cfg.elements) else "Si"
                mass = float(atomic_masses[atomic_numbers[element]])
                cfg.species[int(tid)] = {"symbol": element, "mass": mass}
        else:
            s_type = ph.structure.get('type', 'bulk')
            if s_type == 'file':
                file_path = ph.structure.get('file_path')
                if not file_path or not os.path.exists(file_path):
                    raise FileNotFoundError(f"Structure path missing or not found: {file_path}")
                log.info(f"[PREPARE] Reading {file_path} for {ph.name.upper()} phase")
                base_atoms = read(file_path)
            elif s_type == 'bulk':
                from ase.build import bulk
                params = ph.structure.get('bulk_params', {})
                crystal = params.get('crystal', ph.name)
                a_ref = params.get('a_ref')
                if not a_ref:
                    raise ValueError(f"Missing a_ref in bulk_params for {ph.name}")
                log.info(f"[PREPARE] Generating {crystal} bulk for {ph.name.upper()} phase")
                base_atoms = bulk(cfg.elements[0], crystal, a=a_ref, cubic=True)
            else:
                raise ValueError(f"Unknown structure type: {s_type}")

            # Apply supercell expansion
            atoms = base_atoms.repeat(ph.supercell)
            ph.n_atoms = len(atoms)
            ph.vol_per_atom = atoms.get_volume() / float(len(atoms))
            
            # Map types
            symbols = np.array(atoms.get_chemical_symbols())
            unique_symbols, first_indices = np.unique(symbols, return_index=True)
            order = np.argsort(first_indices)
            unique_symbols = unique_symbols[order]
            
            type_arr = np.zeros(len(atoms), dtype=int)
            mass_arr = np.zeros(len(atoms))
            for i, sym in enumerate(unique_symbols):
                tid = i + 1
                type_arr[symbols == sym] = tid
                mass_val = float(atomic_masses[atomic_numbers[sym]])
                mass_arr[symbols == sym] = mass_val
                cfg.species[tid] = {"symbol": sym, "mass": mass_val}
            
            atoms.set_array('type', type_arr)
            atoms.set_masses(mass_arr)
            ph.n_per_type = {int(t): int(np.sum(type_arr == t)) for t in np.unique(type_arr)}
            
            write(out, atoms, format="lammps-data", atom_style="atomic")

        # Fallback element-level data for logging
        if not hasattr(cfg, 'element') or cfg.element is None:
            if 1 in cfg.species:
                cfg.element = cfg.species[1]["symbol"]
                cfg.mass_amu = cfg.species[1]["mass"]

        log.info(f"  -> {ph.name.upper()} phase prepared. Species counts: {ph.n_per_type}")



# ─────────────────────────────────────────────────────────────────────────────
# Step 1: NPT equilibration + MSD -> spring constants
#         Saves: msd_{phase}.dat, msd_{phase}_plot.png, spring_constants.txt
# ─────────────────────────────────────────────────────────────────────────────
K_OPT: Dict[str, Dict[int, float]] = {}   # {phase_name: {type_id: k_value}}


def _plot_msd(cfg: SystemConfig, ph: PhaseConfig, msd_path: str, k_map: Dict[int, float], avg_msd_map: Dict[int, float]):
    """Plot MSD time-series for all species."""
    data = np.loadtxt(msd_path, skiprows=2)  # LAMMPS ave/time: 2 comment lines
    steps = data[:, 0]
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # We assume col 1 is type 1, col 2 is type 2... based on how we wrote fix ave/time
    for tid, k_val in k_map.items():
        msd_vals = data[:, tid] # column index for tid is tid in our fix ave/time setup
        avg_v = avg_msd_map[tid]
        label = cfg.species[tid]["symbol"]
        line, = ax.plot(steps, msd_vals, linewidth=0.8, label=f"MSD ({label})")
        ax.axhline(avg_v, color=line.get_color(), ls="--", linewidth=1.2,
                   label=f"mean {label} = {avg_v:.4f}")

    ax.set_xlabel("Timestep")
    ax.set_ylabel("MSD (Å²)")
    title = f"{ph.name.upper()} | T={cfg.T0}K | k_avg={np.mean(list(k_map.values())):.3f}"
    ax.set_title(title)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    out = os.path.join(cfg.workdir, f"msd_{ph.name}_plot.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info(f"  MSD plot saved -> {out}")


def step1_spring_constants(cfg: SystemConfig):
    """Compute spring constants from MSD, or load from cache."""
    cache = os.path.join(cfg.workdir, "spring_constants.txt")
    if os.path.exists(cache):
        with open(cache) as f:
            for line in f:
                if "=" not in line:
                    continue
                # Format: k_{phase}_type{tid} = {val}
                key, val = line.split("=")
                key = key.strip()
                if not key.startswith("k_"): continue
                
                parts = key[2:].split("_type")
                if len(parts) == 2:
                    ph_name, tid_str = parts
                    tid = int(tid_str)
                    if ph_name not in K_OPT: K_OPT[ph_name] = {}
                    K_OPT[ph_name][tid] = float(val)
                else:
                    # Legacy support for single element: k_{phase} = {val}
                    ph_name = parts[0]
                    if ph_name not in K_OPT: K_OPT[ph_name] = {}
                    K_OPT[ph_name][1] = float(val)

        # Check if all phases and types are loaded
        all_loaded = True
        for ph in cfg.phases:
            if ph.name not in K_OPT:
                all_loaded = False; break
            for tid in cfg.species:
                if tid not in K_OPT[ph.name]:
                    all_loaded = False; break
        
        if all_loaded:
            log.info(f"[LOAD] spring constants from {cache}")
            return

    for ph in cfg.phases:
        msd_file = f"msd_{ph.name}.dat"
        
        # Build groups and computes for each species
        group_cmds = ""
        compute_cmds = ""
        fix_vars = ""
        for tid in cfg.species:
            group_cmds += f"group type{tid} type {tid}\n            "
            compute_cmds += f"compute msd{tid} type{tid} msd com yes\n            "
            fix_vars += f"c_msd{tid}[4] "
            
        T0 = ph.T0 if ph.T0 is not None else cfg.T0
        script = textwrap.dedent(f"""
            units metal
            atom_style atomic
            boundary p p p
            read_data {cfg.data_file(ph.name)}
            {cfg.get_mass_commands()}
            variable lambda equal 1.0
            {cfg.get_pair_commands(12)}
            velocity all create {T0} {np.random.randint(1,99999)} mom yes rot yes
            fix 1 all npt temp {T0} {T0} 0.5 iso {cfg.pressure} {cfg.pressure} 5.0
            thermo 1000
            thermo_style custom step temp pe press vol
            run {cfg.N_eq}
            unfix 1
            reset_timestep 0
            {group_cmds}
            {compute_cmds}
            fix 2 all nvt temp {T0} {T0} 0.5
            fix 3 all momentum 100 linear 1 1 1
            fix msd_out all ave/time 1 1000 1000 {fix_vars} file {msd_file}
            run {cfg.N_eq}
        """).strip()

        run_lmp(cfg, script, f"in.msd_{ph.name}")

        raw = np.loadtxt(os.path.join(cfg.workdir, msd_file), comments="#")
        # Column 0 is step, Columns 1..N are MSDs for types 1..N
        
        K_OPT[ph.name] = {}
        avg_msd_map = {}
        for i, tid in enumerate(cfg.species):
            col_idx = i + 1
            msd_arr = raw[:, col_idx]
            half = len(msd_arr) // 2
            avg_v = float(np.mean(msd_arr[half:]))
            avg_msd_map[tid] = avg_v
            
            k = (3 * KB * T0) / avg_v
            K_OPT[ph.name][tid] = k
            
            log.info(f"[MSD] {ph.name} Type {tid} ({cfg.species[tid]['symbol']}): <dr2> = {avg_v:.5f} A^2 -> k = {k:.4f}")

        _plot_msd(cfg, ph, os.path.join(cfg.workdir, msd_file), K_OPT[ph.name], avg_msd_map)

    with open(cache, "w") as f:
        # Use first phase T0 as a representative header if needed, but per-type is accurate
        f.write(f"# Spring constants generated with reference T0\n")
        for ph_name, k_map in K_OPT.items():
            for tid, k_val in k_map.items():
                f.write(f"k_{ph_name}_type{tid} = {k_val}\n")
    log.info(f"Spring constants saved -> {cache}")



# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Frenkel-Ladd TI  (lambda 1→0 fwd, 0→1 bwd)
# ─────────────────────────────────────────────────────────────────────────────
def step2_frenkel_ladd(cfg: SystemConfig, ph: PhaseConfig):
    k_map = K_OPT[ph.name]
    out_fwd = f"fl_{ph.name}_fwd.dat"
    out_bwd = f"fl_{ph.name}_bwd.dat"

    # Define the tethering logic
    # If msd_per_element is True, we use a variable for k based on atom type.
    # Otherwise, we use the average k.
    if cfg.msd_per_element and len(cfg.species) > 1:
        # variable k_atom atom (type==1)*k1 + (type==2)*k2 + ...
        k_logic_parts = [f"(type=={tid})*{val}" for tid, val in k_map.items()]
        k_logic = " + ".join(k_logic_parts)
        tether_k = f"variable k_atom atom ({k_logic})"
        k_ref_for_script = "v_k_atom"
    else:
        # Single k for all
        k_avg = np.mean(list(k_map.values()))
        tether_k = f"variable k_fixed equal {k_avg}"
        k_ref_for_script = "v_k_fixed"

    T0 = ph.T0 if ph.T0 is not None else cfg.T0
    # 2a: Equilibrate at lambda=1 (pure potential), save conf
    eq_script = textwrap.dedent(f"""
        units metal
        atom_style atomic
        boundary p p p
        read_data {cfg.data_file(ph.name)}
        {cfg.get_mass_commands()}
        variable lambda equal 1.0
        {cfg.get_pair_commands(8)}
        velocity all create {T0} {np.random.randint(1,99999)} mom yes rot yes
        fix nvt1 all nvt temp {T0} {T0} 0.5
        fix mom  all momentum 100 linear 1 1 1
        thermo 1000
        thermo_style custom step temp pe ke
        run {cfg.N_eq}
        write_data conf_eq_{ph.name}.data nocoeff
    """).strip()
    run_lmp(cfg, eq_script, f"in.eq_{ph.name}")

    # 2b: Forward + Backward switching
    for direction, lam_s, lam_e, out in [
            ("fwd", 1.0, 0.0, out_fwd),
            ("bwd", 0.0, 1.0, out_bwd)]:
        script = textwrap.dedent(f"""
            units metal
            atom_style atomic
            boundary p p p
            read_data conf_eq_{ph.name}.data
            {cfg.get_mass_commands()}
            
            # lambda: scales EAM potential (1=full EAM, 0=Einstein crystal)
            variable lambda  equal ramp({lam_s},{lam_e})
            variable inv_lam equal 1.0-v_lambda
            
            # tethering spring k*(1-λ)
            {tether_k}
            variable k_sp    atom {k_ref_for_script}*v_inv_lam

            {cfg.get_pair_commands(12)}

            fix tether all spring/self v_k_sp
            fix nvt1   all nvt temp {T0} {T0} 0.5
            fix mom    all momentum 100 linear 1 1 1

            compute pe_eam all pe pair
            variable teth_e equal f_tether

            thermo 1000
            thermo_style custom step temp v_lambda c_pe_eam v_teth_e

            variable s equal step
            variable l equal v_lambda
            variable e equal c_pe_eam
            variable t equal v_teth_e

            fix print_out all print 100 "$s $l $e $t" &
                title "step lambda pe_eam pe_tether" &
                file {out} screen no

            run {cfg.N_fl}
        """).strip()
        run_lmp(cfg, script, f"in.fl_{ph.name}_{direction}")
        log.info(f"[FL] {ph.name} {direction} -> {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Reversible Scaling  (lambda = T0/T sweep)
# ─────────────────────────────────────────────────────────────────────────────
def step3_reversible_scaling(cfg: SystemConfig, ph: PhaseConfig):
    for direction, T_target in [("up", ph.T_min), ("down", ph.T_max)]:
        T0 = ph.T0 if ph.T0 is not None else cfg.T0
        if T_target == T0:
            log.info(f"[RS] Skipping '{direction}' sweep because T_target == T0 ({T_target}K)")
            continue

        out_fwd = f"rs_{ph.name}_{direction}_fwd.dat"
        out_bwd = f"rs_{ph.name}_{direction}_bwd.dat"

        if cfg.rs_scheduling == "linear-T":
            lam_logic_fwd = textwrap.dedent(f"""
                variable f_fwd    equal (elapsed/v_t_total)
                variable lambda   equal 1.0/(1.0 + v_f_fwd*(v_T_target/v_T0 - 1.0))
            """).strip()
            lam_logic_bwd = textwrap.dedent(f"""
                variable f_bwd    equal (elapsed/v_t_total)
                variable lambda   equal 1.0/(1.0 + (1.0 - v_f_bwd)*(v_T_target/v_T0 - 1.0))
            """).strip()
        else:
            # linear-lambda
            lam_e = cfg.T0 / T_target
            lam_logic_fwd = f"variable lambda equal ramp(1.0, {lam_e})"
            lam_logic_bwd = f"variable lambda equal ramp({lam_e}, 1.0)"

        script = textwrap.dedent(f"""
            units metal
            atom_style atomic
            boundary p p p
            read_data conf_eq_{ph.name}.data
            {cfg.get_mass_commands()}
            
            variable T_target equal {T_target}
            variable T0       equal {T0}
            variable t_total  equal {cfg.N_rs}

            # ── 1. Forward Sweep (T0 -> T_target) ──
            {lam_logic_fwd}
            variable T_eff    equal {T0}/v_lambda

            {cfg.get_pair_commands(12)}

            fix nvt1 all nvt temp {T0} {T0} 0.5
            fix mom  all momentum 100 linear 1 1 1

            compute pe_pair all pe pair

            variable s  equal step
            variable la equal v_lambda
            variable te equal v_T_eff
            variable ep equal c_pe_pair

            fix print_fwd all print 200 "${{s}} ${{la}} ${{te}} ${{ep}}" &
                title "step lambda T_eff pe_eam" file {out_fwd} screen no

            thermo 2000
            thermo_style custom step temp v_lambda v_T_eff c_pe_pair

            run {cfg.N_rs}
            unfix print_fwd

            # ── 2. Equilibrate at T_target ──
            run {cfg.N_eq}

            # ── 3. Backward Sweep (T_target -> T0) ──
            {lam_logic_bwd}
            fix print_bwd all print 200 "${{s}} ${{la}} ${{te}} ${{ep}}" &
                title "step lambda T_eff pe_eam" file {out_bwd} screen no
                
            run {cfg.N_rs}
            unfix print_bwd
        """).strip()
        run_lmp(cfg, script, f"in.rs_{ph.name}_{direction}")
        log.info(f"[RS] {ph.name} {direction} bidir ({cfg.T0}K <-> {T_target}K) -> fwd/bwd.dat")


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Analysis  (F_Ein, F_CM, FL work, RS curve, Tc)
# ─────────────────────────────────────────────────────────────────────────────

# Global variable to track current phase context for Nj-per-type lookups
ph_current_context: Optional[PhaseConfig] = None

def _f_einstein(T: float, k_map: Dict[int, float], Nat: int, species: Dict[int, dict]) -> float:
    """
    Helmholtz free energy of Einstein crystal in eV.
    F_Ein = sum_j [ N_j * (3*kB*T*ln(hbar*omega_j/kB*T)) ]
    where omega_j = sqrt(k_j/m_j)
    """
    if T <= 0: return 0.0
    f_total = 0.0
    for tid, k in k_map.items():
        mass_amu = species[tid]["mass"]
        omega = np.sqrt(k / mass_amu) * 9.82269e13
        val = 3.0 * KB * T * np.log(HBAR_SI * omega / (KB_SI * T))
        # nj logic
        nj = Nat if len(k_map) == 1 else (Nat / len(k_map))
        if ph_current_context and hasattr(ph_current_context, 'n_per_type') and tid in ph_current_context.n_per_type:
            nj = ph_current_context.n_per_type[tid]
        f_total += nj * val
    return f_total


def _f_cm(T: float, V: float, Nat: int, k_map: Dict[int, float]) -> float:
    """
    Center-of-mass correction (eV)  [FreeEnergyLAMMPS integrate.py L49].
    F_CM = kB*T * ln( (N/V) * (2*pi*kB*T / (N*k_avg))^1.5 )
    Uses V/N as the accessible CM volume (periodic crystal, see Sum.md §3.1).
    """
    if T <= 0:
        return 0.0
    k_avg = float(np.mean(list(k_map.values())))
    kBT = KB * T  # eV
    arg = (Nat / V) * (2.0 * np.pi * kBT / (Nat * k_avg)) ** 1.5
    return kBT * np.log(arg)

def _integrate_fl(cfg: SystemConfig, ph: PhaseConfig):
    """
    Frenkel-Ladd integration [FreeEnergyLAMMPS Eq.12]:
      W = (I_forward - I_backward) / 2
    """
    global ph_current_context
    ph_current_context = ph
    
    I = {}
    for direction, sign in (("fwd", -1), ("bwd", +1)):
        fname = os.path.join(cfg.workdir, f"fl_{ph.name}_{direction}.dat")
        data  = np.loadtxt(fname, skiprows=1)
        lam   = data[:, 1]
        pe_e  = data[:, 2]   # lambda * U_EAM
        pe_t  = data[:, 3]   # (1-lambda) * U_Ein
        lam_safe = np.clip(lam, 1e-9, 1 - 1e-9)
        U_EAM = pe_e / lam_safe
        U_Ein = pe_t / (1 - lam_safe)
        dU    = U_EAM - U_Ein            # dH/dlambda
        # sign: fwd runs 1->0, trapz gives negative, negate to get 0->1 convention
        I[direction] = sign * trapezoid(dU, lam)  # both now: eV, 0->1 direction

    W_avg = (I["fwd"] + I["bwd"]) / 2.0
    Nat   = ph.n_atoms
    hys_per_atom = abs(I["fwd"] - I["bwd"]) / Nat

    log.info(f"[FL-{ph.name}] Integrals: I_fwd={I['fwd']:.4f} eV, I_bwd={I['bwd']:.4f} eV")
    log.info(f"[FL-{ph.name}] Decision: Average Work W = {W_avg:.4f} eV ({W_avg/Nat:.5f} eV/atom)")
    
    if hys_per_atom > cfg.hysteresis_warn:
        log.warning(f"[FL-{ph.name}] High hysteresis: {hys_per_atom:.5f} eV/atom > limit {cfg.hysteresis_warn}")
    else:
        log.info(f"[FL-{ph.name}] Hysteresis check passed: {hys_per_atom:.5f} eV/atom")
    return I["fwd"], I["bwd"], W_avg



def _integrate_rs_bidirectional(fwd_file: str, bwd_file: str, N: int):
    from scipy.integrate import cumulative_trapezoid
    
    # FWD: lambda goes 1.0 -> lam_e
    data_f = np.loadtxt(fwd_file, skiprows=1)
    lam_f = data_f[:, 1]
    T_eff_f = data_f[:, 2]
    U_f = data_f[:, 3] / np.clip(lam_f, 1e-10, None) / N
    
    # BWD: lambda goes lam_e -> 1.0
    data_b = np.loadtxt(bwd_file, skiprows=1)
    lam_b = data_b[:, 1]
    T_eff_b = data_b[:, 2]
    U_b = data_b[:, 3] / np.clip(lam_b, 1e-10, None) / N
    
    # Cumulative integral W(lam) = int_1^lam U dlam'
    # For FWD: Since lam_f starts at 1.0, integration with lam yields int_1^lam directly
    W_f = cumulative_trapezoid(U_f, lam_f, initial=0)
    
    # For BWD: Starts at lam_e, ends at 1.0. W_raw(lam) = int_{lam_e}^lam U dlam'
    W_b_raw = cumulative_trapezoid(U_b, lam_b, initial=0)
    # The total integral from lam_e to 1.0 is the last element
    I_total = W_b_raw[-1]
    # We want W_b(lam) = int_1^lam U dlam' = int_{lam_e}^lam U dlam' - int_{lam_e}^1 U dlam'
    W_b = W_b_raw - I_total
    
    # Average W_b and W_f onto the same lambda grid.
    # lam_f is decreasing (e.g. 1.0 -> 0.1), lam_b is increasing (e.g. 0.1 -> 1.0).
    # numpy interp requires the x-coordinates to be increasing.
    W_b_interp = np.interp(lam_f[::-1], lam_b, W_b)[::-1]
    
    W_avg = 0.5 * (W_f + W_b_interp)
    return lam_f, T_eff_f, W_f, W_b_interp, W_avg


def _plot_rs_stability(cfg: SystemConfig, ph: PhaseConfig, T: np.ndarray, F_fwd: np.ndarray, F_bwd: np.ndarray, F_avg: np.ndarray):
    """Plot Fwd vs Bwd vs Average for a single phase to check RS convergence."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(T, F_fwd, "g--", lw=1.0, label="Forward (T0 -> T_target)")
    ax.plot(T, F_bwd, "b--", lw=1.0, label="Backward (T_target -> T0)")
    ax.plot(T, F_avg, "r-", lw=1.5, label="Average (Hysteresis-corrected)")
    
    # Shade the hysteresis gap
    ax.fill_between(T, F_fwd, F_bwd, color="gray", alpha=0.1, label="Hysteresis Loop")
    
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("F (eV/atom)")
    ax.set_title(f"RS Stability: {cfg.element} {ph.name.upper()} | t_s = {cfg.N_rs}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(cfg.workdir, f"rs_convergence_{ph.name}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info(f"[PLOT] RS Convergence check saved -> {out}")


def _integrate_rs(cfg: SystemConfig, ph: PhaseConfig):
    """
    RS reference formula [FreeEnergyLAMMPS Eq.21-22]:
      1) Unscale energy: U_true = pe_recorded / lambda
      2) Cumulative integral W_avg from bidirectional sweeps.
      3) F(T) = F0/lambda + 1.5*kB*T*ln(lambda) + W_avg  [per atom]
    """
    lam_full, T_full, F_avg_full = [], [], []
    F_fwd_full, F_bwd_full = [], []

    # Process "up" direction (T0 -> T_min)
    if cfg.T_min != cfg.T0:
        fwd = os.path.join(cfg.workdir, f"rs_{ph.name}_up_fwd.dat")
        bwd = os.path.join(cfg.workdir, f"rs_{ph.name}_up_bwd.dat")
        log.info(f"  [RS-Analysis] Processing {ph.name} UP bidirectional pass")
        lam_up, T_up, W_f, W_b, W_a = _integrate_rs_bidirectional(fwd, bwd, ph.n_atoms)
        
        # Calculate F_T0_pa needed for F(T) assembly
        # We'll pull this from the first call in step4_analysis? No, let's keep it abstract.
        # But wait, we need F_T0_pa to calculate F_fwd/F_bwd. 
        # For now let's just use the logic to return these W components.
        
        lam_full.extend(lam_up[::-1])
        T_full.extend(T_up[::-1])
        # We will collect W and convert to F later in step4_analysis.
        return lam_up, T_up, W_f, W_b, W_a

    # Process "down" direction (T0 -> T_max)
    if cfg.T_max != cfg.T0:
        fwd = os.path.join(cfg.workdir, f"rs_{ph.name}_down_fwd.dat")
        bwd = os.path.join(cfg.workdir, f"rs_{ph.name}_down_bwd.dat")
        log.info(f"  [RS-Analysis] Processing {ph.name} DOWN bidirectional pass")
        return _integrate_rs_bidirectional(fwd, bwd, ph.n_atoms)

    return np.array([1.0]), np.array([cfg.T0]), np.array([0.0]), np.array([0.0]), np.array([0.0])



def _plot_hysteresis(cfg: SystemConfig, results: dict):
    """
    Bar chart of |hysteresis| in eV/atom for each phase.
    """
    phases  = list(results.keys())
    hys_vals = [results[ph]["hys_per_atom"] for ph in phases]
    colors  = ["steelblue" if v < cfg.hysteresis_warn else "tomato" for v in hys_vals]

    fig, ax = plt.subplots(figsize=(5, 4))
    bars    = ax.bar(phases, hys_vals, color=colors, edgecolor="k", linewidth=0.8)
    ax.axhline(cfg.hysteresis_warn, color="red", ls="--", linewidth=1,
               label=f"threshold = {cfg.hysteresis_warn} eV/atom")
    for bar, val in zip(bars, hys_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                f"{val:.4f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("|W_fwd − W_bwd| / N_atoms  (eV/atom)")
    ax.set_title(f"{cfg.element} FL Hysteresis Check  (T0={cfg.T0} K)")
    ax.set_ylim(0, max(hys_vals) * 1.3 + cfg.hysteresis_warn * 0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(cfg.workdir, "hysteresis_check.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info(f"[PLOT] Hysteresis bar chart -> {out}")


def step4_analysis(cfg: SystemConfig, active_phases: List[PhaseConfig]):
    results = {}

    for ph in active_phases:
        Nat = ph.n_atoms
        k_map = K_OPT[ph.name]
        V_tot = ph.vol_per_atom * Nat           # Å³ (approximate lattice)

        # Context for Nj tracking
        global ph_current_context
        ph_current_context = ph
        
        F_ein = _f_einstein(cfg.T0, k_map, Nat, cfg.species)
        F_cm  = _f_cm(cfg.T0, V_tot, Nat, k_map)

        I_fwd, I_bwd, W_avg = _integrate_fl(cfg, ph)
        hys_per_atom = abs(I_fwd - I_bwd) / Nat

        # F_solid(T0) = F_Ein + F_CM + W  [eV, whole box]
        F_T0  = F_ein + F_cm + W_avg
        F_T0_pa = F_T0 / Nat

        log.info(f"\n[{ph.name.upper()}] N={Nat}")
        log.info(f"  F_Ein = {F_ein:.4f} eV ({F_ein/Nat:.5f} eV/atom)")
        log.info(f"  F_CM  = {F_cm:.5f} eV ({F_cm/Nat:.6f} eV/atom)")
        log.info(f"  W_FL  = {W_avg:.4f} eV ({W_avg/Nat:.5f} eV/atom)")
        log.info(f"  F(T0) = {F_T0:.4f} eV ({F_T0_pa:.5f} eV/atom)")

        # RS: F(T) = F0/lambda + 1.5*kB*T*ln(lambda) + W_rs/lambda  [per atom]
        # [FreeEnergyLAMMPS RS integrate.py Line 30, 34]
        lam_rs, T_rs, W_f, W_b, W_a = _integrate_rs(cfg, ph)
        
        # Assemble components
        F_fwd = F_T0_pa / lam_rs + 1.5 * KB * T_rs * np.log(lam_rs) + W_f / lam_rs
        F_bwd = F_T0_pa / lam_rs + 1.5 * KB * T_rs * np.log(lam_rs) + W_b / lam_rs
        F_avg = F_T0_pa / lam_rs + 1.5 * KB * T_rs * np.log(lam_rs) + W_a / lam_rs

        # Convergence Plot (visualize loop)
        _plot_rs_stability(cfg, ph, T_rs, F_fwd, F_bwd, F_avg)

        results[ph.name] = {
            "T": T_rs, "F_per_atom": F_avg,
            "hys_per_atom": hys_per_atom,
            "F_T0_pa": F_T0_pa,
            "F_fwd": F_fwd, "F_bwd": F_bwd
        }

    # ── Hysteresis plot ──────────────────────────────────────────────────────
    _plot_hysteresis(cfg, results)

    # ── Free energy crossing ─────────────────────────────────────────────────
    names     = [ph.name for ph in cfg.phases if ph.name in results]
    if len(names) == 0:
        return None
        
    for n in names:
        if n == 'solid':
            np.savez("solid_ft.npz", T=results[n]["T"], F=results[n]["F_per_atom"])
            
    T_min_all = max(results[n]["T"].min() for n in names)
    T_max_all = min(results[n]["T"].max() for n in names)
    T_grid    = np.linspace(T_min_all, T_max_all, 600)

    F_interp = {n: np.interp(T_grid, results[n]["T"], results[n]["F_per_atom"])
                for n in names}

    Tc = None
    if len(names) == 2:
        diff = F_interp[names[0]] - F_interp[names[1]]
        crosses = np.where(np.diff(np.sign(diff)))[0]
        if crosses.size:
            i  = crosses[0]
            dT = T_grid[1] - T_grid[0]
            Tc = T_grid[i] + (-diff[i]) / (diff[i+1] - diff[i]) * dT
            log.info(f"\n[RESULT] {names[0].upper()}-{names[1].upper()} "
                     f"crossing Tc = {Tc:.1f} K")
        else:
            log.warning("[RESULT] No crossing found in temperature range.")

    # ── Free energy plot ─────────────────────────────────────────────────────
    colors = plt.cm.tab10.colors
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, n in enumerate(names):
        ax.plot(results[n]["T"], results[n]["F_per_atom"],
                color=colors[i], lw=1.5, label=n.upper())
    if Tc:
        ax.axvline(Tc, color="k", ls="--", lw=1.2,
                   label=f"Tc ≈ {Tc:.0f} K")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("F (eV/atom)")
    ax.set_title(f"{cfg.element}  Free Energy: {' vs '.join(n.upper() for n in names)}")
    ax.legend()
    fig.tight_layout()
    out_fe = os.path.join(cfg.workdir, f"{cfg.element}_free_energy.png")
    fig.savefig(out_fe, dpi=150)
    plt.close(fig)
    log.info(f"[PLOT] Free energy -> {out_fe}")

    return Tc


def _plot_rs_rate_convergence(cfg: SystemConfig, ph: PhaseConfig, records: list):
    """Plot G(T_target) vs Total Steps (ts) matching Fig 2 in the user's reference."""
    import matplotlib.pyplot as plt
    records = sorted(records, key=lambda x: x["steps"])
    steps = [r["steps"] / 1000.0 for r in records]  # in 10^3 steps
    f_fwd = [r["fwd"] for r in records]
    f_bwd = [r["bwd"] for r in records]
    f_avg = [r["avg"] for r in records]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(steps, f_fwd, "go-", ms=6, lw=1.0, label="Forward (T0 -> T_target)")
    ax.plot(steps, f_bwd, "bo-", ms=6, lw=1.0, label="Backward (T_target -> T0)")
    ax.plot(steps, f_avg, "ro-", ms=7, lw=1.5, label="Average")
    
    # reference line if N_rs is high
    ax.axhline(f_avg[-1], color="k", ls="--", lw=0.8, alpha=0.5, label=f"Final Avg: {f_avg[-1]:.4f}")

    ax.set_xlabel("t_s (10^3 MD steps)")
    ax.set_ylabel(f"Free Energy at {cfg.T_max} K (eV/atom)")
    ax.set_title(f"RS Convergence Analysis: {cfg.element} {ph.name.upper()}")
    ax.legend(fontsize=9)
    ax.grid(True, ls=":", alpha=0.6)
    
    fig.tight_layout()
    out = os.path.join(cfg.workdir, f"convergence_{ph.name}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info(f"[PLOT] RS Convergence Rate plot -> {out}")


def step5_convergence(cfg: SystemConfig, active_phases: List[PhaseConfig]):
    """Automatically run RS multiple times with different step counts to check rate convergence."""
    # Define a range of step counts
    base = cfg.N_rs
    # Use a set of fractions of the base steps, or absolute counts
    test_steps = [20000, 40000, 60000, 80000, 100000, 150000]
    # Filter out steps that are too similar or redundant if user set a small N_rs
    test_steps = sorted(list(set([s for s in test_steps if s <= base] + [base])))

    for ph in active_phases:
        log.info(f"\n[CONVERGENCE] Starting RS Sweep for {ph.name.upper()}...")
        records = []
        
        # Need F(T0) from FL step (Step 2 must have been run)
        # We RE-USE the spring constants and FL work already calculated.
        k_map = K_OPT[ph.name]
        Nat = ph.n_atoms
        V_tot = ph.vol_per_atom * Nat
        # Fetching W_avg from step 2
        try:
            _, _, W_avg = _integrate_fl(cfg, ph)
        except Exception as e:
            log.error(f"[CONVERGENCE] Failed to get FL work for {ph.name}. Ensure --step fl or all was run.")
            continue

        F_T0_pa = (_f_einstein(cfg.T0, k_map, Nat, cfg.species) + 
                   _f_cm(cfg.T0, V_tot, Nat, k_map) + 
                   W_avg) / Nat

        # Backup original N_rs
        orig_nrs = cfg.N_rs
        
        try:
            for s in test_steps:
                cfg.N_rs = s
                log.info(f"  -> Testing t_s = {s} steps...")
                step3_reversible_scaling(cfg, ph)
                
                # Analyze this pass
                lam_rs, T_rs, W_f, W_b, W_a = _integrate_rs(cfg, ph)
                
                # Extract values at the endpoint (T_max)
                idx_end = -1
                t_actual = T_rs[idx_end]
                f_fwd = F_T0_pa / lam_rs[idx_end] + 1.5 * KB * t_actual * np.log(lam_rs[idx_end]) + W_f[idx_end] / lam_rs[idx_end]
                f_bwd = F_T0_pa / lam_rs[idx_end] + 1.5 * KB * t_actual * np.log(lam_rs[idx_end]) + W_b[idx_end] / lam_rs[idx_end]
                f_avg = F_T0_pa / lam_rs[idx_end] + 1.5 * KB * t_actual * np.log(lam_rs[idx_end]) + W_a[idx_end] / lam_rs[idx_end]
                
                records.append({"steps": s, "fwd": f_fwd, "bwd": f_bwd, "avg": f_avg})
                log.info(f"     G({t_actual:.0f}K): Fwd={f_fwd:.5f}, Bwd={f_bwd:.5f}, Avg={f_avg:.5f}")

            # General Plotting
            _plot_rs_rate_convergence(cfg, ph, records)
            
        finally:
            cfg.N_rs = orig_nrs


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    # ── CLI ──────────────────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        description="Solid-phase TI: hybrid/scaled Frenkel-Ladd + Reversible Scaling")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--phase", default="all",
                        help="Phase name or 'all' (default all)")
    parser.add_argument("--step", default="all",
                        choices=["all", "structures", "msd", "fl", "rs", "analysis", "convergence"],
                        help="Pipeline step to run (default all)")
    args = parser.parse_args()

    # ── Config ───────────────────────────────────────────────────────────────
    cfg = SystemConfig(args.config)    # Parse from YAML

    active_phases = (cfg.phases if args.phase == "all"
                     else [cfg.phase(args.phase)])

    # Ensure structure and mass data are globally available before TI logic
    prepare_structures(cfg)

    log.info("=" * 60)
    log.info(f"System : {cfg.element if hasattr(cfg, 'element') else 'Multi-Species'} | Phases: {[p.name for p in cfg.phases]}")
    log.info(f"N_eq={cfg.N_eq}  N_fl={cfg.N_fl}  N_rs={cfg.N_rs}")
    for ph in active_phases:
        T0 = ph.T0 if ph.T0 is not None else cfg.T0
        Tmin = ph.T_min if ph.T_min is not None else cfg.T_min
        Tmax = ph.T_max if ph.T_max is not None else cfg.T_max
        log.info(f"Phase {ph.name}: T0={T0} K, Range=[{Tmin}, {Tmax}] K")
    if cfg.msd_per_element:
        log.info(f"MSD Mode: Per-Element Spring Constants")
    log.info("=" * 60)

    # ── Execute ───────────────────────────────────────────────────────────────
    # We define a helper to ensure spring constants are loaded/run
    def ensure_spring_constants():
        step1_spring_constants(cfg)

    if args.step in ("all", "structures"):
        pass  # structures prepared automatically above

    if args.step == "msd":
        ensure_spring_constants()
        return

    if args.step in ("all", "fl"):
        ensure_spring_constants()
        for ph in active_phases:
            step2_frenkel_ladd(cfg, ph)

    if args.step in ("all", "rs"):
        ensure_spring_constants()
        for ph in active_phases:
            step3_reversible_scaling(cfg, ph)

    if args.step in ("all", "analysis"):
        ensure_spring_constants()
        step4_analysis(cfg, active_phases)

    if args.step == "convergence":
        ensure_spring_constants()
        step5_convergence(cfg, active_phases)



if __name__ == "__main__":
    main()

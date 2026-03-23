"""
ti_common.py  –  Shared utilities for TI / QHA workflows.

Consolidates physical constants, LAMMPS runner, and force/energy
calculators that were previously duplicated across run_msd.py,
run_qha.py, and run_ti.py.
"""

from __future__ import annotations
import os, sys, subprocess, logging
from typing import Optional
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Physical Constants  (CODATA 2018)
# ─────────────────────────────────────────────────────────────────────────────
KB      = 8.617333e-5          # eV/K   (Boltzmann)
HBAR_SI = 1.054571817e-34      # J·s    (reduced Planck)
KB_SI   = 1.380649e-23         # J/K    (Boltzmann, SI)
H_SI    = 6.62607015e-34       # J·s    (Planck)
EV_TO_J = 1.60218e-19          # J/eV
AMU_KG  = 1.66053906660e-27    # kg/amu

# Phonopy unit-conversion factor  (LAMMPS metal → THz)
# Nu [THz] = sqrt(k[eV/Å²] / m[amu]) * 15.633302
PHONOPY_FACTOR_LAMMPS = 15.633302

# kJ/mol → eV/particle
KJ_MOL_TO_EV = 0.01036427


# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS Runner
# ─────────────────────────────────────────────────────────────────────────────
_log = logging.getLogger("ti_common")


def run_lammps(
    lmp_exe: str,
    script_text: str,
    script_name: str,
    workdir: str = ".",
    log_file: Optional[str] = None,
    capture: bool = False,
    logger: Optional[logging.Logger] = None,
) -> subprocess.CompletedProcess:
    """
    Write *script_text* to *script_name* inside *workdir* and execute LAMMPS.

    Parameters
    ----------
    lmp_exe : str
        Path to the LAMMPS executable.
    script_text : str
        Full LAMMPS input script content.
    script_name : str
        Filename for the input script (e.g. ``in.msd_bcc``).
    workdir : str
        Working directory for LAMMPS execution.
    log_file : str or None
        If given, ``-log <log_file>`` is appended to the LAMMPS command.
    capture : bool
        If True, capture stdout/stderr (silent mode). Default False.
    logger : Logger or None
        Logger instance for info/error messages.

    Returns
    -------
    subprocess.CompletedProcess
    """
    lg = logger or _log
    path = os.path.join(workdir, script_name)
    with open(path, "w") as f:
        f.write(script_text)

    cmd = [lmp_exe, "-in", script_name]
    if log_file:
        cmd += ["-log", log_file]

    lg.info(f"Running LAMMPS: {script_name}")
    result = subprocess.run(cmd, cwd=workdir, capture_output=capture)
    if result.returncode != 0:
        lg.error(f"LAMMPS failed: {script_name}")
        sys.exit(1)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Static Force / Energy calculators  (used by QHA)
# ─────────────────────────────────────────────────────────────────────────────
def get_forces_lammps(
    lmp_exe: str,
    atoms,
    workdir: str,
    pot_style: str,
    pot_coeff: str,
    name: str = "tmp",
) -> np.ndarray:
    """
    Run a single-point LAMMPS calculation and return per-atom forces [N, 3].

    The *atoms* object must be an ASE ``Atoms`` instance.
    """
    from ase.io import write

    data_path = os.path.join(workdir, f"{name}.data")
    write(data_path, atoms, format="lammps-data", atom_style="atomic")

    if isinstance(pot_style, list):
        styles = " ".join(pot_style)
        pair_cmds = f"pair_style hybrid/overlay {styles}\n"
        for coeff in pot_coeff:
            pair_cmds += f"pair_coeff {coeff}\n"
    else:
        pair_cmds = f"pair_style {pot_style}\npair_coeff {pot_coeff}"

    script = f"""\
units metal
atom_style atomic
boundary p p p
read_data {name}.data

{pair_cmds}

thermo 1
thermo_style custom step pe press
dump 1 all custom 1 {name}.dump id fx fy fz
run 0
"""
    run_lammps(lmp_exe, script, f"in.{name}", workdir,
               log_file=f"log.{name}", capture=True)

    forces = []
    dump_file = os.path.join(workdir, f"{name}.dump")
    with open(dump_file, "r") as f:
        capture_flag = False
        for line in f:
            if "ITEM: ATOMS id fx fy fz" in line:
                capture_flag = True
                continue
            if capture_flag:
                vals = line.split()
                if len(vals) == 4:
                    forces.append([int(vals[0]), float(vals[1]),
                                   float(vals[2]), float(vals[3])])

    forces.sort(key=lambda x: x[0])
    return np.array([f[1:] for f in forces])


def get_energy_lammps(
    lmp_exe: str,
    atoms,
    workdir: str,
    pot_style: str,
    pot_coeff: str,
    name: str = "energy",
) -> float:
    """
    Run a single-point LAMMPS calculation and return the total potential energy.
    """
    from ase.io import write

    data_path = os.path.join(workdir, f"{name}.data")
    write(data_path, atoms, format="lammps-data", atom_style="atomic")

    if isinstance(pot_style, list):
        styles = " ".join(pot_style)
        pair_cmds = f"pair_style hybrid/overlay {styles}\n"
        for coeff in pot_coeff:
            pair_cmds += f"pair_coeff {coeff}\n"
    else:
        pair_cmds = f"pair_style {pot_style}\npair_coeff {pot_coeff}"

    script = f"""\
units metal
atom_style atomic
read_data {name}.data
{pair_cmds}
thermo 1
run 0
"""
    run_lammps(lmp_exe, script, f"in.{name}", workdir,
               log_file=f"log.{name}", capture=True)

    pe = 0.0
    log_path = os.path.join(workdir, f"log.{name}")
    with open(log_path, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "PotEng" in line:
                for j in range(i + 1, len(lines)):
                    parts = lines[j].split()
                    if len(parts) > 2 and parts[0].isdigit():
                        pe = float(parts[2])
                        break
                break
    return pe

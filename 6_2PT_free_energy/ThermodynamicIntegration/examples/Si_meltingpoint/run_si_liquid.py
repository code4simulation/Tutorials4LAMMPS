#!/usr/bin/env python3
import os, sys, argparse, textwrap, logging, yaml
import numpy as np
from scipy.integrate import trapezoid, cumulative_trapezoid
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("SI-LIQUID")

sys.path.append(os.path.abspath("."))
sys.path.append(os.path.abspath("../Fe")) # For ti_common
from ti_common import run_lammps, KB, HBAR_SI, KB_SI, EV_TO_J, AMU_KG
import uf_math

def load_cfg(path="config_si.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)

def generate_table(sigma, epsilon=0.05, r_outer=8.0, n_points=5000):
    r_arr = np.linspace(0.01, r_outer, n_points)
    with open("uf.table", "w") as f:
        f.write("# UF table\n\nUF_POT\n")
        f.write(f"N {n_points} R 0.01 {r_outer}\n\n")
        for i, r in enumerate(r_arr):
            v = max(1 - np.exp(-(r/sigma)**2), 1e-15)
            ene = -epsilon * np.log(v)
            frc = epsilon * (np.exp(-(r/sigma)**2) / v) * (2*r / (sigma**2))
            f.write(f"{i+1} {r:.6f} {ene:.6f} {frc:.6f}\n")

def liquid_equil(cfg_yaml, T0, N_eq):
    log.info(f"=== Equilibrating Liquid at {T0}K ===")
    
    script = textwrap.dedent(f"""
        units metal
        boundary p p p
        atom_style atomic
        
        lattice diamond 5.43
        region box block 0 4 0 4 0 4
        create_box 1 box
        create_atoms 1 box
        mass 1 28.0855
        
        pair_style sw
        pair_coeff * * Si.sw Si
        
        # Melt at 4000K, extremely fast
        velocity all create 4000.0 12345 rot yes mom yes
        fix 1 all npt temp 4000.0 4000.0 0.5 iso 0.0 0.0 5.0
        run 10000
        unfix 1
        
        # Cool down to T0
        fix 2 all npt temp {T0} {T0} 0.5 iso 0.0 0.0 5.0
        thermo 1000
        run {N_eq}
        unfix 2
        
        # NVT to sample Density
        fix 3 all nvt temp {T0} {T0} 0.5
        variable d equal density
        fix def all ave/time 1 {max(1, N_eq//4)} {N_eq} v_d file density.dat
        run {N_eq}
        
        write_data conf_liquid.data nocoeff
    """)
    run_lammps("lmp", script, "in.liq_eq", workdir=".", logger=log)
    
    # Extract density
    data = np.atleast_2d(np.loadtxt("density.dat", skiprows=2))
    vol_density = float(data[-1, 1]) # atoms / A^3
    log.info(f"Computed Volume Density: {vol_density:.5f} atoms/A^3")
    return vol_density

def liquid_fl(cfg_yaml, rho, T0, N_fl):
    
    sig_opt = (2.0 / (np.pi**1.5 * rho))**(1/3)
    epsilon = 0.05
    log.info(f"Optimized UFM sigma for target x=1.0: {sig_opt:.4f} A")
    
    generate_table(sigma=sig_opt, epsilon=epsilon)
    
    for d, ls, le in [("fwd", 1, 0), ("bwd", 0, 1)]:
        script = textwrap.dedent(f"""
            units metal
            boundary p p p
            atom_style atomic
            read_data conf_liquid.data
            
            variable lambda equal ramp({ls},{le})
            variable inv_lam equal 1.0-v_lambda
            
            pair_style hybrid/scaled v_lambda sw v_inv_lam table linear 10000
            pair_coeff * * sw Si.sw Si
            pair_coeff 1 1 table uf.table UF_POT 8.0
            
            fix 1 all nvt temp {T0} {T0} 0.5
            
            variable lam equal v_lambda
            variable st equal step
            dump 1 all atom 100 dump.liq_fl_{d}.lammpstrj
            fix p1 all print 100 "${{st}} ${{lam}}" file lam_{d}.txt screen no
            run {N_fl}
            undump 1
        """)
        run_lammps("lmp", script, f"in.liq_fl_{d}", workdir=".", logger=log)
        
        # Rerun to extract pure SW and pure UF energies
        for pot_type, p_style, p_coeff in [
            ("sw", "sw", "Si.sw Si"),
            ("uf", "table linear 10000", "uf.table UF_POT 8.0")
        ]:
            rerun_script = textwrap.dedent(f"""
                units metal
                boundary p p p
                atom_style atomic
                read_data conf_liquid.data
                
                pair_style {p_style}
                pair_coeff * * {p_coeff}
                
                variable e equal pe
                fix p2 all print 100 "${{e}}" file ene_{pot_type}_{d}.txt screen no
                rerun dump.liq_fl_{d}.lammpstrj dump x y z box no
            """)
            run_lammps("lmp", rerun_script, f"in.rerun_{pot_type}_{d}", workdir=".", logger=log)
        
    # Assemble integral
    I_vals = {}
    for d in ["fwd", "bwd"]:
        lam = np.loadtxt(f"lam_{d}.txt", skiprows=1)[:, 1]
        u_sw = np.loadtxt(f"ene_sw_{d}.txt", skiprows=1)
        u_uf = np.loadtxt(f"ene_uf_{d}.txt", skiprows=1)
        
        dU_dlam = u_sw - u_uf
        sign = -1 if d == "fwd" else 1
        I_vals[d] = sign * trapezoid(dU_dlam, lam)
        log.info(f"Liquid FL {d} Integral: {I_vals[d]:.5f} eV (Total)")
        
    W_avg = (I_vals["fwd"] + I_vals["bwd"]) / 2.0
    return sig_opt, W_avg

def liquid_rs(cfg_yaml, T0, T_min, T_max, N_rs):
    
    for direction, T_target in [("up", T_min), ("down", T_max)]:
        if T_target == T0: continue
        
        lam_logic_fwd = textwrap.dedent(f"""
            variable f_fwd    equal (elapsed/v_t_total)
            variable lambda   equal 1.0/(1.0 + v_f_fwd*(v_T_target/v_T0 - 1.0))
        """).strip()
        lam_logic_bwd = textwrap.dedent(f"""
            variable f_bwd    equal (elapsed/v_t_total)
            variable lambda   equal 1.0/(1.0 + (1.0 - v_f_bwd)*(v_T_target/v_T0 - 1.0))
        """).strip()

        script = textwrap.dedent(f"""
            units metal
            boundary p p p
            atom_style atomic
            read_data conf_liquid.data
            
            variable T_target equal {T_target}
            variable T0       equal {T0}
            variable t_total  equal {N_rs}
            
            {lam_logic_fwd}
            variable T_eff    equal {T0}/v_lambda
            
            pair_style sw
            pair_coeff * * Si.sw Si
            
            fix 1 all nvt temp {T0} {T0} 0.5
            compute pe_pair all pe pair
            
            variable s  equal step
            variable la equal v_lambda
            variable te equal v_T_eff
            variable ep equal c_pe_pair
            
            fix print_fwd all print 200 "${{s}} ${{la}} ${{te}} ${{ep}}" file rs_liq_{direction}_fwd.dat screen no
            run {N_rs}
            unfix print_fwd
            
            # Equilibrate at T_target (effective)
            run 10000
            
            {lam_logic_bwd}
            fix print_bwd all print 200 "${{s}} ${{la}} ${{te}} ${{ep}}" file rs_liq_{direction}_bwd.dat screen no
            run {N_rs}
            unfix print_bwd
        """)
        run_lammps("lmp", script, f"in.liq_rs_{direction}", workdir=".", logger=log)

def integrate_rs_bidirectional(fwd_file, bwd_file, N):
    data_f = np.loadtxt(fwd_file, skiprows=1)
    lam_f = data_f[:, 1]
    T_eff_f = data_f[:, 2]
    U_f = data_f[:, 3] / np.clip(lam_f, 1e-10, None) / N
    
    data_b = np.loadtxt(bwd_file, skiprows=1)
    lam_b = data_b[:, 1]
    T_eff_b = data_b[:, 2]
    U_b = data_b[:, 3] / np.clip(lam_b, 1e-10, None) / N
    
    W_f = cumulative_trapezoid(U_f, lam_f, initial=0)
    W_b_raw = cumulative_trapezoid(U_b, lam_b, initial=0)
    W_b = W_b_raw - W_b_raw[-1]
    W_b_interp = np.interp(lam_f[::-1], lam_b, W_b)[::-1]
    W_avg = 0.5 * (W_f + W_b_interp)
    return lam_f, T_eff_f, W_f, W_b_interp, W_avg

def process_liquid(cfg_yaml, N_atoms):
    # Find liquid phase params
    liq_ph = next((p for p in cfg_yaml.get('phases', []) if p['name'] == 'liquid'), {})
    params = cfg_yaml.get('parameters', {})
    
    T0    = liq_ph.get('T0', params.get('T0', 1500.0))
    T_min = liq_ph.get('T_min', params.get('T_min', 1000.0))
    T_max = liq_ph.get('T_max', params.get('T_max', 2000.0))
    N_eq  = params.get('N_eq', 10000)
    N_fl  = params.get('N_fl', 40000)
    N_rs  = params.get('N_rs', 60000)
    
    rho = liquid_equil(cfg_yaml, T0, N_eq)
    sig_opt, W_fl = liquid_fl(cfg_yaml, rho, T0, N_fl)
    
    F_ideal_pa = uf_math.get_ideal_gas_fe(T0, rho, 28.0855) 
    F_res_pa = uf_math.get_uhlenbeck_ford_fe(T0, rho, sig_opt)
    
    F_T0_pa = F_ideal_pa + F_res_pa + (W_fl / N_atoms)
    log.info(f"=== Liquid F(T0) Analysis ===")
    log.info(f" Analytic F_ideal = {F_ideal_pa:.6f} eV/atom")
    log.info(f" Analytic F_res_UF = {F_res_pa:.6f} eV/atom")
    log.info(f" W_FL = {W_fl / N_atoms:.6f} eV/atom")
    log.info(f" Total F_liquid({T0}K) = {F_T0_pa:.6f} eV/atom")
    
    # Run RS
    liquid_rs(cfg_yaml, T0, T_min, T_max, N_rs)
    
    lam_full, T_full, F_avg_full = [], [], []
    for direction in ["up", "down"]:
        fwd = f"rs_liq_{direction}_fwd.dat"
        bwd = f"rs_liq_{direction}_bwd.dat"
        if not os.path.exists(fwd): continue
        
        lam_rs, T_rs, W_f, W_b, W_a = integrate_rs_bidirectional(fwd, bwd, N_atoms)
        
        if direction == "up":
            lam_full.extend(lam_rs[::-1])
            T_full.extend(T_rs[::-1])
            F_avg_full.extend( (F_T0_pa / lam_rs + 1.5 * KB * T_rs * np.log(lam_rs) + W_a / lam_rs)[::-1] )
        else:
            lam_full.extend(lam_rs)
            T_full.extend(T_rs)
            F_avg_full.extend(F_T0_pa / lam_rs + 1.5 * KB * T_rs * np.log(lam_rs) + W_a / lam_rs)
            
    np.savez("liquid_ft.npz", T=T_full, F=F_avg_full)
    
    plt.plot(T_full, F_avg_full, 'ro-', label='Liquid')
    plt.xlabel("T (K)"); plt.ylabel("Free Energy (eV/atom)")
    plt.legend(); plt.savefig("liquid_fe.png")
    log.info("Saved liquid_fe.png and liquid_ft.npz")

if __name__ == "__main__":
    cfg = load_cfg()
    # N is defined by the 4 4 4 supercell of diamond -> 512 atoms
    process_liquid(cfg, 512)

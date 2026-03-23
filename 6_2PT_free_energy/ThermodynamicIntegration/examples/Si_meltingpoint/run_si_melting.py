#!/usr/bin/env python3
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=" * 60)
    print(" 1) Running Solid Phase Thermodynamic Integration")
    print("=" * 60)
    res_solid = subprocess.run([sys.executable, "run_si_solid.py", "--config", "config_si.yaml", "--phase", "solid"])
    if res_solid.returncode != 0:
        print("Error in Solid TI pipeline.")
        sys.exit(1)

    print("=" * 60)
    print(" 2) Running Liquid Phase Thermodynamic Integration (UFM Reference)")
    print("=" * 60)
    res_liquid = subprocess.run([sys.executable, "run_si_liquid.py"])
    if res_liquid.returncode != 0:
        print("Error in Liquid TI pipeline.")
        sys.exit(1)

    print("=" * 60)
    print(" 3) Analyzing Melting Point Intersection")
    print("=" * 60)
    
    if not os.path.exists("solid_ft.npz") or not os.path.exists("liquid_ft.npz"):
        print("Missing output files (solid_ft.npz or liquid_ft.npz).")
        sys.exit(1)
        
    sol = np.load("solid_ft.npz")
    liq = np.load("liquid_ft.npz")
    
    T_s, F_s = sol['T'], sol['F']
    T_l, F_l = liq['T'], liq['F']
    
    # Quadratic Fit for robust extrapolation
    p_s = np.polyfit(T_s, F_s, 2)
    p_l = np.polyfit(T_l, F_l, 2)
    
    # Common grid for intersection search and plotting
    T_grid = np.linspace(1000, 2200, 2400)
    F_s_fit = np.polyval(p_s, T_grid)
    F_l_fit = np.polyval(p_l, T_grid)
    
    # Analytical intersection of quadratic fits
    # (p_s[0]-p_l[0])T^2 + (p_s[1]-p_l[1])T + (p_s[2]-p_l[2]) = 0
    diff_coeffs = p_s - p_l
    roots = np.roots(diff_coeffs)
    valid_roots = [r.real for r in roots if np.isreal(r) and 1000 <= r.real <= 2200]
    
    Tc = valid_roots[0] if valid_roots else None
    if Tc:
        print(f"\n>>> EXTRAPOLATED MELTING POINT (Tc) = {Tc:.1f} K <<<\n")
    else:
        print("\n>>> WARNING: No melting point intersection found in 1000K-2200K. <<<\n")

    # Plot
    plt.figure(figsize=(8, 6))
    # Original Data
    plt.plot(T_s, F_s, 'bo', ms=4, alpha=0.3, label='Solid Data')
    plt.plot(T_l, F_l, 'ro', ms=4, alpha=0.3, label='Liquid Data')
    
    # Fitted/Extrapolated Lines
    # Use Solid line where T < Tc and Dashed where T > Tc (vice versa for liquid)
    mask_s = T_grid <= (Tc if Tc else 2200)
    mask_l = T_grid >= (Tc if Tc else 1000)
    
    plt.plot(T_grid, F_s_fit, 'b-', lw=2, label='Solid Fit/Extrap')
    plt.plot(T_grid, F_l_fit, 'r-', lw=2, label='Liquid Fit/Extrap')
    
    if Tc is not None:
        Fc = np.polyval(p_s, Tc)
        plt.plot(Tc, Fc, 'k*', ms=12, label=f'Melting Point ≈ {Tc:.1f} K')
        plt.axvline(Tc, color='k', ls='--', alpha=0.3)

    plt.xlabel("Temperature (K)")
    plt.ylabel("Free Energy (eV/atom)")
    plt.title("Silicon Melting Point Extrapolation (Stillinger-Weber)")
    plt.legend()
    plt.grid(ls=':', alpha=0.7)
    plt.xlim(1000, 2200)
    plt.tight_layout()
    plt.savefig("si_melting_point.png", dpi=200)
    print("Saved plot to 'si_melting_point.png'")

if __name__ == "__main__":
    main()

"""
Thermochemistry & Free Energy Module
Author: Antigravity (Scientific Intelligence & Validation Lead)
Date: 2026-04-16
Calculates Gibbs free energies for molecules (gas phase), substrates (solids), 
and adsorbates, mimicking the capabilities of VASPKIT (e.g., module 511, 512).
Physics/Algorithm:
------------------
1. Substrate/Adsorbate: Uses the Harmonic approximation.
   G(T) = E_DFT + ZPE + U_vib(T) - T * S_vib(T)
   (PV term is assumed to be ~0 for condensed phases).
2. Molecule: Uses the Ideal Gas approximation.
   G(T, P) = E_DFT + ZPE + ΔH_trans + ΔH_rot + ΔH_vib - T*(S_trans + S_rot + S_vib)
Note on Imaginary Frequencies:
------------------------------
Real systems or stable intermediates should only have positive frequencies. 
If evaluating a Transition State (TS), the imaginary mode (representing the reaction 
coordinate) must be dropped before calculating the partition function.
"""
from typing import List, Tuple, Dict
import numpy as np
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.units import invcm
from ase import Atoms
def clean_frequencies(frequencies_cm_inv: List[float], is_transition_state: bool = False) -> np.ndarray:
    """
    Cleans frequency list by removing imaginary frequencies (negative or complex).
    Interprets the structural stability based on the number of imaginary modes.
    
    Args:
        frequencies_cm_inv: List of frequencies in cm^-1.
        is_transition_state: Flag indicating if the user expects a transition state.
        
    Returns:
        np.ndarray: Array of vibrational energies in eV (imaginary modes removed).
    """
    valid_freqs = []
    imag_freqs = []
    
    for f in frequencies_cm_inv:
        if hasattr(f, "imag") and f.imag != 0:
            imag_freqs.append(f)
        elif getattr(f, "real", f) < 0:
            imag_freqs.append(f)
        else:
            valid_freqs.append(getattr(f, "real", f))
            
    n_imag = len(imag_freqs)
    
    print("\n" + "="*40)
    print(f"{'Structural Stability Assessment':^40}")
    print("="*40)
    
    if n_imag == 0:
        print("Status: Stable Local Minimum (Ground State)")
        if is_transition_state:
            print("[Warning] Expected a Transition State, but found 0 imaginary modes.")
    elif n_imag == 1:
        print("Status: 1st-Order Saddle Point (Transition State)")
        if not is_transition_state:
            print("[Warning] Structure is a TS, but `is_transition_state=False`.")
            print("          Thermodynamic properties correspond to an activated state.")
    else:
        print(f"Status: {n_imag}th-Order Saddle Point (Highly Unstable)")
        print("[CRITICAL WARNING] Multiple imaginary frequencies detected!")
        print("This structure is NOT a valid local minimum nor a standard transition state.")
        print("The calculated free energy will likely be physically meaningless.")
        
    print("="*40)
        
    # Convert cm^-1 to eV (1 cm^-1 ~ 1.23984e-4 eV)
    energies_ev = np.array(valid_freqs) * invcm
    return energies_ev
def calc_substrate_free_energy(
    electronic_energy: float,
    frequencies_cm_inv: List[float],
    temperature: float = 298.15,
    is_transition_state: bool = False
) -> float:
    """
    Calculates the Free Energy of a solid substrate using Harmonic approximation.
    """
    vib_energies = clean_frequencies(frequencies_cm_inv, is_transition_state)
    thermo = HarmonicThermo(vib_energies=vib_energies, potentialenergy=electronic_energy)
    
    # get_helmholtz_energy returns A = E_pot + ZPE + U_vib - T*S_vib (approx G for solids)
    free_energy = thermo.get_helmholtz_energy(temperature=temperature)
    
    print("\n" + "-"*40)
    print(f"{'Substrate / Solid Free Energy':^40}")
    print("-" * 40)
    print(f"Temperature         : {temperature:10.2f} K")
    print(f"Electronic Energy   : {electronic_energy:10.4f} eV")
    zpe = np.sum(vib_energies) / 2.0
    print(f"ZPE                 : {zpe:10.4f} eV")
    print(f"Entropy (S_vib)     : {thermo.get_entropy(temperature):10.7f} eV/K")
    print(f"Gibbs Free Energy   : {free_energy:10.4f} eV")
    print("-" * 40)
    
    return free_energy
def calc_adsorbate_free_energy(
    electronic_energy: float,
    frequencies_cm_inv: List[float],
    temperature: float = 298.15,
    is_transition_state: bool = False
) -> float:
    """
    Calculates the Free Energy of an adsorbate. 
    In standard HTST models, an adsorbate's translational and rotational 
    degrees of freedom are converted into frustrated vibrations.
    Therefore, the Harmonic approximation is mathematically identical.
    """
    print("\n[Note] For Adsorbates, translational/rotational modes are treated as frustrated vibrations.")
    return calc_substrate_free_energy(electronic_energy, frequencies_cm_inv, temperature, is_transition_state)
def calc_molecule_free_energy(
    electronic_energy: float,
    frequencies_cm_inv: List[float],
    atoms: Atoms,
    geometry: str,
    symmetrynumber: int = 1,
    spin: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 101325.0,  # 1 atm in Pa
    is_transition_state: bool = False
) -> float:
    """
    Calculates Free Energy of a molecule in the Ideal Gas phase.
    
    Args:
        electronic_energy: E_DFT (eV)
        frequencies_cm_inv: List of frequencies (cm^-1)
        atoms: ASE Atoms object of the molecule (needed for mass/inertia)
        geometry: "monatomic", "linear", or "nonlinear"
        symmetrynumber: Rotational symmetry number (e.g. 2 for H2, 12 for CH4)
        spin: Total electronic spin (0 for singlet, 0.5 for doublet, etc.)
        temperature: Temperature in K
        pressure: Pressure in Pa
        is_transition_state: Flag indicating if the molecule is a TS
    """
    vib_energies = clean_frequencies(frequencies_cm_inv, is_transition_state)
    
    thermo = IdealGasThermo(
        vib_energies=vib_energies,
        potentialenergy=electronic_energy,
        atoms=atoms,
        geometry=geometry,
        symmetrynumber=symmetrynumber,
        spin=spin
    )
    
    free_energy = thermo.get_gibbs_energy(temperature=temperature, pressure=pressure)
    
    print("\n" + "-"*40)
    print(f"{'Molecule (Ideal Gas) Free Energy':^40}")
    print("-" * 40)
    print(f"Temperature         : {temperature:10.2f} K")
    print(f"Pressure            : {pressure:10.1f} Pa")
    print(f"Electronic Energy   : {electronic_energy:10.4f} eV")
    zpe = np.sum(vib_energies) / 2.0
    print(f"ZPE                 : {zpe:10.4f} eV")
    print(f"Entropy (S_tot)     : {thermo.get_entropy(temperature, pressure):10.7f} eV/K")
    print(f"Gibbs Free Energy   : {free_energy:10.4f} eV")
    print("-" * 40)
    
    return free_energy
if __name__ == "__main__":
    from ase import Atoms
    
    # Mock Validation Data (E.g. from VASP OUTCAR or PHVA analysis)
    E_dft = -25.500 # eV
    mock_freqs = [3000.0, 1500.0, 1000.0, -150.0] # 1 imaginary mode included
    
    print("="*40)
    print(f"{'THERMOCHEMISTRY VALIDATION SUITE':^40}")
    print("="*40)
    
    # 1. Substrate / Solid validation
    calc_substrate_free_energy(
        electronic_energy=E_dft,
        frequencies_cm_inv=mock_freqs,
        temperature=300.0,
        is_transition_state=True # We expect a TS since there is 1 imaginary mode
    )
    
    # 2. Molecule Gas Validation (e.g. H2O)
    h2o = Atoms('H2O', positions=[(0, 0, 0), (0, 0.76, 0.58), (0, -0.76, 0.58)])
    calc_molecule_free_energy(
        electronic_energy=E_dft,
        frequencies_cm_inv=mock_freqs,
        atoms=h2o,
        geometry='nonlinear',
        symmetrynumber=2, # C2v symmetry
        temperature=300.0
    )

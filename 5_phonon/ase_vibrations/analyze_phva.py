"""
PHVA (Partial Hessian Vibrational Analysis) Script
Author: Antigravity (Scientific Intelligence & Validation Lead)
Date: 2026-04-16
This script performs localized vibrational analysis on a specified subset of atoms
using the ASE-LAMMPS interface. It is primarily used for calculating prefactors
in Harmonic Transition State Theory (HTST).
Physics/Algorithm:
------------------
The vibrational frequencies ν are calculated from the eigenvalues λ of the
mass-weighted Hessian matrix H̃. In PHVA, only a subset of atoms (S) is
displaced, effectively assuming the rest of the system is frozen (infinite mass).
H̃_ij = (1 / sqrt(m_i * m_j)) * (∂²V / ∂R_i ∂R_j)
λ_k = (2πν_k)²
Units:
------
- Energy: eV
- Length: Å
- Mass: amu
- Frequency: cm⁻¹ (standard ASE output)
"""
import os
from typing import List, Dict, Any
from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
from ase.vibrations import Vibrations
def setup_lammps_calculator(
    executable: str,
    potential_settings: Dict[str, Any]
) -> LAMMPS:
    """
    Configures and returns an ASE LAMMPS calculator.
    
    Args:
        executable: Path to the LAMMPS binary (e.g., 'lmp_serial').
        potential_settings: Dictionary containing arbitrary ASE LAMMPS kwargs 
                            (e.g., 'pair_style', 'pair_coeff', 'files', 'parameters').        
    Returns:
        LAMMPS: Initialized ASE calculator.
    """
    # Explicitly link the LAMMPS binary via the 'executable' argument.
    # Alternatively, one could set os.environ['ASE_LAMMPSRUN_COMMAND'] = f"{executable} -in %s"
    
    # Use dictionary unpacking to allow for flexible potential configurations
    # such as 'hybrid/overlay' and explicitly provided parameters/files.
    calc = LAMMPS(
        executable=executable,
        **potential_settings
    )
    return calc
def run_phva(
    atoms_path: str,
    active_indices: List[int],
    lammps_exe: str,
    potential_dict: Dict[str, Any],
    delta: float = 0.01
):
    """
    Executes the Partial Hessian Vibrational Analysis.
    
    Args:
        atoms_path: Path to the structure file (e.g., geometry.in).
        active_indices: List of 0-based atomic indices for the partial Hessian.
        lammps_exe: Path to LAMMPS executable.
        potential_dict: Dictionary with LAMMPS potential parameters.
        delta: Displacement distance for finite difference (Å).
    """
    # 1. Load Geometry
    # units: Å, amu
    atoms = read(atoms_path)
    
    # 2. Attach Calculator
    calc = setup_lammps_calculator(lammps_exe, potential_dict)
    atoms.set_calculator(calc)
    
    # 3. Setup Vibrations (Partial Hessian)
    # The 'indices' argument restricts displacements to the specified atoms.
    # vibration_name is used for naming the .pkl log files.
    vib_name = "phva_results"
    if os.path.exists(vib_name):
        import shutil
        shutil.rmtree(vib_name) # Clean previous runs
        
    vib = Vibrations(atoms, indices=active_indices, name=vib_name, delta=delta)
    
    # 4. Run Calculation
    # Performs 3 * len(active_indices) * 2 force calls.
    print(f"Starting PHVA for {len(active_indices)} active atoms...")
    vib.run()
    
    # 5. Extract and Report Frequencies
    # units: cm⁻¹
    frequencies = vib.get_frequencies()
    
    print("\n" + "="*40)
    print(f"{'Vibrational Frequency Analysis (PHVA)':^40}")
    print("="*40)
    print(f"Active Indices: {active_indices}")
    print(f"Total Modes:   {len(frequencies)}")
    print("-"*40)
    
    for i, f in enumerate(frequencies):
        if isinstance(f, complex) or f < 0:
            # ASE represent imaginary frequencies as complex or negative depending on version
            # Here we print as 'imaginary' for TS identification.
            val = f.imag if hasattr(f, "imag") else abs(f)
            print(f"Mode {i:3d}: {val:10.4f} cm⁻¹ (Imaginary/TS)")
        else:
            print(f"Mode {i:3d}: {f:10.4f} cm⁻¹")
    
    print("="*40)
    
    # Cleanup (Optional)
    # vib.clean() # Uncomment to remove .pkl files after execution
    
    return frequencies
def run_fhva(
    atoms_path: str,
    lammps_exe: str,
    potential_dict: Dict[str, Any],
    delta: float = 0.01
) -> List[float]:
    """
    Executes the Full Hessian Vibrational Analysis (FHVA).
    
    Args:
        atoms_path: Path to the structure file (e.g., geometry.in).
        lammps_exe: Path to LAMMPS executable.
        potential_dict: Dictionary with LAMMPS potential parameters.
        delta: Displacement distance for finite difference (Å).
    """
    atoms = read(atoms_path)
    calc = setup_lammps_calculator(lammps_exe, potential_dict)
    atoms.set_calculator(calc)
    
    vib_name = "fhva_results"
    if os.path.exists(vib_name):
        import shutil
        shutil.rmtree(vib_name) # Clean previous runs
        
    vib = Vibrations(atoms, name=vib_name, delta=delta)
    
    print(f"Starting FHVA for all {len(atoms)} atoms...")
    vib.run()
    
    frequencies = vib.get_frequencies()
    
    print("\n" + "="*40)
    print(f"{'Vibrational Frequency Analysis (FHVA)':^40}")
    print("="*40)
    print(f"Total Modes:   {len(frequencies)}")
    print("-"*40)
    
    for i, f in enumerate(frequencies[:10]):  # Show max 10 to prevent terminal flood
        if isinstance(f, complex) or f < 0:
            val = f.imag if hasattr(f, "imag") else abs(f)
            print(f"Mode {i:3d}: {val:10.4f} cm⁻¹ (Imaginary/TS)")
        else:
            print(f"Mode {i:3d}: {f:10.4f} cm⁻¹")
            
    if len(frequencies) > 10:
         print(f"... and {len(frequencies) - 10} more modes.")
         
    print("="*40)
    return frequencies
if __name__ == "__main__":
    # --- Input Parameters ---
    STRUCTURE_FILE = "geometry.in"
    LAMMPS_EXECUTABLE = "/path/to/lmp_serial"  # User must update this path
    
    # Example: Hybrid Overlay (Tersoff + ZBL)
    # ZBL parameters: <inner_cutoff> <outer_cutoff> and Zi Zj (14.0 for Silicon)
    POTENTIAL_DICT = {
        "pair_style": "hybrid/overlay tersoff zbl 4.0 4.8",
        "pair_coeff": [
            "* * tersoff Si.tersoff Si",
            "* * zbl 14.0 14.0"
        ],
        "files": ["Si.tersoff"],
        "parameters": {
            "mass": ["1 28.0855"]
        }
    }
    
    # Provided indices: [100, 101, 102, 103, 104, 105]
    ACTIVE_ATOM_INDICES = [100, 101, 102, 103, 104, 105]
    DISPLACEMENT_DELTA = 0.01  # units: Å
    
    # check if file exists before running
    if os.path.exists(STRUCTURE_FILE):
        import numpy as np
        
        print("\n[Validation Suite] Verifying PHVA vs FHVA")
        atoms = read(STRUCTURE_FILE)
        all_indices = list(range(len(atoms))) # Select all atoms for equivalence test
        
        # 1. Run FHVA
        freqs_fhva = run_fhva(STRUCTURE_FILE, LAMMPS_EXECUTABLE, POTENTIAL_DICT, delta=DISPLACEMENT_DELTA)
        
        # 2. Run PHVA on all indices
        freqs_phva = run_phva(STRUCTURE_FILE, all_indices, LAMMPS_EXECUTABLE, POTENTIAL_DICT, delta=DISPLACEMENT_DELTA)
        
        # 3. Direct compare
        mag_fhva = np.abs(freqs_fhva)
        mag_phva = np.abs(freqs_phva)
        
        max_diff = np.max(np.abs(mag_fhva - mag_phva))
        
        print("\n" + "="*40)
        print(f"{'Validation Report: FHVA vs PHVA':^40}")
        print("="*40)
        if max_diff < 1e-4:
            print("Result: SUCCESS! Identical frequencies.")
        else:
            print("Result: WARNING! Discrepancy observed.")
        print(f"Max absolute difference: {max_diff:.2e} cm⁻¹")
        print("="*40)
        
    else:
        print(f"Error: {STRUCTURE_FILE} not found. Please ensure the geometry file exists.")

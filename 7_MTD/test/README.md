# Metadynamics Research Environment (MTD)

This folder contains a dedicated Python virtual environment and resources for Metadynamics simulations using **LAMMPS** and **PLUMED**.

## 1. Environment Activation
To use the dedicated environment, activate it using PowerShell or Command Prompt:

**PowerShell:**
```powershell
.\venv\Scripts\Activate.ps1
```

**Command Prompt:**
```cmd
venv\Scripts\activate.bat
```

## 2. Installed Packages
- **Scientific:** `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`
- **Metadynamics:** `plumed` (Python interface)
- **Interactive:** `jupyterlab`, `notebook`

## 3. Simulation Engine Setup (LAMMPS + PLUMED)
Since you are on Windows, there are two ways to run the simulations:

### Option A: WSL2 (Windows Subsystem for Linux) - **Recommended**
LAMMPS and PLUMED are natively supported and more stable on Linux.
1. Install WSL2: `wsl --install`
2. Install Ubuntu from the Microsoft Store.
3. Follow the [PLUMED Installation Guide](https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html) inside Ubuntu.

### Option B: Native Windows
If you already have LAMMPS installed on Windows:
1. Ensure the LAMMPS executable (e.g., `lmp.exe`) is in your system `PATH`.
2. To use PLUMED with LAMMPS, you must have a version of LAMMPS that was compiled with the `PKG-PLUMED` package.
3. You can verify your LAMMPS packages by running:
   ```cmd
   lmp -h
   ```
   Look for `plumed` in the "Installed packages" section.

## 4. Getting Started
Open the included Jupyter notebook to start with the theory and basic analysis:
```bash
jupyter notebook Getting_Started_Metadynamics.ipynb
```

---
> [!NOTE]
> This environment focuses on **Metadynamics Theory and Practice**. MD Analysis libraries (like MDAnalysis or pyscal) have been excluded per your request, but can be added later if needed.

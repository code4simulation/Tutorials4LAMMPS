import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.integrate import simps
from sklearn.linear_model import LinearRegression

def extract_temperature(filename):
    match = re.search(r'MSD_(\d+).txt', filename)
    if match:
        return int(match.group(1))
    else:
        return None

def calculate_diffusion_msd(filename):
    data = np.loadtxt(filename, skiprows=2)
    steps = data[:, 0]
    msd = data[:, 1]
    
    # Fit a linear model to the MSD data to find the slope
    model = LinearRegression().fit(steps.reshape(-1, 1), msd)
    slope = model.coef_[0]
    
    # Calculate the diffusion coefficient (D)
    D_msd = slope / 6  # for 3D diffusion
    
    return D_msd, steps, msd

def calculate_diffusion_vacf(filename):
    data = np.loadtxt(filename, skiprows=2)
    steps = data[:, 0]
    vacf = data[:, 3]
    
    # Integrate VACF to find the diffusion coefficient
    D_vacf = simps(vacf, steps) / 3  # for 3D diffusion
    
    return D_vacf, steps, vacf

# List of MSD files
msd_files = glob.glob("MSD_*.txt")

# Sort files by extracted temperature
msd_files.sort(key=extract_temperature)

# Dictionaries to store results
msd_results = {}
vacf_results = {}
temperatures = []

# Process each MSD file
fig, ax = plt.subplots(constrained_layout=True)
for msd_file in msd_files:
    temperature = extract_temperature(msd_file)
    temperatures.append(temperature)
    
    # Calculate diffusion coefficients using MSD and VACF
    D_msd, steps_msd, msd = calculate_diffusion_msd(msd_file)
    D_vacf, steps_vacf, vacf = calculate_diffusion_vacf(msd_file)
    
    msd_results[temperature] = D_msd
    vacf_results[temperature] = D_vacf
    
    # Plot MSD vs TimeStep
    ax[0].plot(steps_msd, msd, label=f'T={temperature}K')
    ax[0].set_xlabel('Simulation Step')
    ax[0].set_ylabel('Mean Squared Displacement (MSD)')
    ax[0].set_title('MSD vs Simulation Step')
    
    # Plot VACF vs TimeStep
    ax[1].plot(steps_vacf, vacf, label=f'T={temperature}K')
    ax[1].set_xlabel('Simulation Step')
    ax[1].set_ylabel('Velocity Auto-Correlation Function (VACF)')
    ax[1].set_title('VACF vs Simulation Step')

ax[0].legend()
ax[1].legend()

# Plot diffusion coefficients comparison
temperatures.sort()
D_msd_values = [msd_results[temp] for temp in temperatures]
D_vacf_values = [vacf_results[temp] for temp in temperatures]

ax[2].plot(temperatures, D_msd_values, 'o-', label='MSD Method')
ax[2].plot(temperatures, D_vacf_values, 's-', label='VACF Method')
ax[2].set_xlabel('Temperature (K)')
ax[2].set_ylabel('Diffusion Coefficient (D)')
ax[2].set_title('Diffusion Coefficient vs Temperature')
ax[2]t.legend()

# Print diffusion coefficients for each temperature
print("Temperature (K) | Diffusion Coefficient (MSD) | Diffusion Coefficient (VACF)")
print("----------------------------------------------------------------------------")
for temp in temperatures:
    print(f"{temp:15} | {msd_results[temp]:27.4E} | {vacf_results[temp]:27.4f}")

plt.show()

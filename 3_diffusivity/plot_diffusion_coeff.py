import numpy as np
import scipy.constants as CONST
import matplotlib.pyplot as plt
import glob
import re

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
    coeffs = np.polyfit(steps, msd, 1)
    model = np.poly1d(coeffs)

    # Calculate the diffusion coefficient (D)
    D_msd = coeffs[0] / 6  # for 3D diffusion

    return D_msd, steps, msd

# List of MSD files
msd_files = glob.glob("MSD_*.txt")

# Sort files by extracted temperature
msd_files.sort(key=extract_temperature)

# Dictionaries to store results
msd_results = {}
temperatures = []

# Process each MSD file
fig, ax = plt.subplots(1,3, figsize=(15,5), constrained_layout=True)
for msd_file in msd_files:
    temperature = extract_temperature(msd_file)
    temperatures.append(temperature)

    # Calculate diffusion coefficients using MSD
    D_msd, steps_msd, msd = calculate_diffusion_msd(msd_file)

    msd_results[temperature] = D_msd

    # Plot MSD vs TimeStep
    ax[0].plot(steps_msd, msd, label=f'T={temperature}K')
    ax[0].set_xlabel('Simulation Step')
    ax[0].set_ylabel('Mean Squared Displacement (MSD)')
    ax[0].set_title('MSD vs Simulation Step')

# Plot diffusion coefficients
temperatures.sort()
temperatures = np.array(temperatures)
D_msd_values = [msd_results[temp] for temp in temperatures]

ax[1].plot(temperatures, D_msd_values, 'o-', label='MSD Method')
ax[1].set_xlabel('Temperature (K)')
ax[1].set_ylabel('Diffusion Coefficient (D)')
ax[1].set_title('Diffusion Coefficient vs Temperature')

logD = np.log(D_msd_values)
coeffs = np.polyfit(1/temperatures, logD, 1)
print(f"Ea: {-coeffs[0]*CONST.k/CONST.e:.2f} (eV)")
print(f"D0: {np.exp(coeffs[1]):.2E} (cm2/s)")
fit_line = np.poly1d(coeffs)
ax[2].plot(1/temperatures, logD, 'o-', label='MSD Method')
ax[2].plot(1/temperatures, fit_line(1/temperatures), 'r--', label='Fit line')
ax[2].set_xlabel(r'1/T (K$^{-1}$)')
ax[2].set_ylabel('log(D)')
ax[2].set_title('Diffusion Coefficient vs Temperature')

for i in range(3):
    ax[i].legend()
plt.show()

# Print diffusion coefficients for each temperature
print(f"{'Temperature (K)':^15s} | {'Diffusion Coefficient (MSD)':^27s} |")
print("-----------------------------------------------")
for temp in temperatures:
    print(f"{temp:^15d} | {msd_results[temp]:^27.4E}")

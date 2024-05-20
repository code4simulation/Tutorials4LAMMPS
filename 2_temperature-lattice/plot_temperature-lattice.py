import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# initial setting
temperatures = range(10, 101, 10)
lattice_parameters = []
fig, ax = plt.subplots(1,2, constrained_layout=True)
ax[0].set_xlabel('Step')
ax[0].set_ylabel(r'Lattice Parameter ($\mathrm{\AA}$)')

# Post-processing
for temp in temperatures:
    filename = f"2_{temp}.dat"
    data = np.loadtxt(filename, skiprows=1)
    
    # Average
    lattice_param_avg = np.mean(data[-100:, 4])
    lattice_parameters.append(lattice_param_avg)
    
    # Print the average lattice paramter at temp
    print(f"Temperature: {temp}K, Average Lattice Parameter: {lattice_param_avg:.4f} ansgtrom")

    # Each temperature
    ax[0].plot(data[:,0], data[:,4], label=f'{temp} K')

# Plot
ax[1].plot(temperatures, lattice_parameters, marker='o')
ax[1].set_xlabel('Temperature (K)')
ax[1].set_ylabel(r'Lattice Parameter ($\mathrm{\AA}$)')
plt.show()

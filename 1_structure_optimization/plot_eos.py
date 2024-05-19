import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# Birch-Murnaghan equation of state
def birch_murnaghan(V, E0, V0, B0, B0_prime):
    term1 = (V0 / V) ** (2 / 3)
    term2 = (term1 - 1)
    return E0 + (9 * V0 * B0 / 16) * ((term2 ** 3) * B0_prime + (term2 ** 2) * (6 - 4 * term1))

# Load data
data = np.loadtxt('eos.dat')
volumes = data[:, 2]
energies = data[:, 1]

# Set an initial values
E0_initial = min(energies)
V0_initial = volumes[np.argmin(energies)]
B0_initial = 1.0
B0_prime_initial = 4.0

# Fitting
params_initial = [E0_initial, V0_initial, B0_initial, B0_prime_initial]
params_opt, params_cov = scipy.optimize.curve_fit(birch_murnaghan, volumes, energies, p0=params_initial)

# Print outputs
E0, V0, B0, B0_prime = params_opt
print(f"Minimum Energy (E0): {E0}")
print(f"Equilibrium Volume (V0): {V0}")
print(f"Bulk Modulus (B0): {B0}")
print(f"Bulk Modulus Derivative (B0_prime): {B0_prime}")

# Plotting
V_fit = np.linspace(min(volumes), max(volumes), 100)
E_fit = birch_murnaghan(V_fit, *params_opt)

plt.plot(volumes, energies, 'o', label='Data')
plt.plot(V_fit, E_fit, '-', label='Fit')
plt.xlabel('Volume (Å^3)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.show()

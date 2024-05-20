import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

file_path = sys.argv[1]
data = pd.read_csv(file_path, delim_whitespace=True, header=None,
                   names=['step', 'temperature', 'energy', 'atom energy','lattice'],
                   skiprows=1)
print(data.head())

def analyze_convergence(data):
    summary = data.describe()
    print(summary)
    return summary

def plot_convergence(data):
    plt.figure(figsize=(20, 5))
    plt.subplot(1, 4, 1)
    plt.plot(data['step'], data['temperature'], label='Temperature')
    plt.xlabel('Step')
    plt.ylabel('Temperature')
    plt.title('Temperature vs Step')
    plt.legend()

    plt.subplot(1, 4, 2)
    plt.plot(data['step'], data['energy'], label='Energy', color='orange')
    plt.xlabel('Step')
    plt.ylabel('Energy')
    plt.title('Energy vs Step')
    plt.legend()

    plt.subplot(1, 4, 3)
    plt.plot(data['step'], data['atom energy'], label='Per-atom Energy', color='orange')
    plt.xlabel('Step')
    plt.ylabel('Per-atom Energy')
    plt.title('Energy vs Step')
    plt.legend()

    plt.subplot(1, 4, 4)
    plt.plot(data['step'], data['lattice'], label='Lattice', color='green')
    plt.xlabel('Step')
    plt.ylabel('Lattice')
    plt.title('Lattice vs Step')
    plt.legend()

    plt.tight_layout()
    plt.show()

# Analysis
summary = analyze_convergence(data)

# Plot
plot_convergence(data)

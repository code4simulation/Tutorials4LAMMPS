import sys, os
import itertools
from ase.io import read
from ase.data import atomic_numbers, atomic_masses

input_elem = sys.argv[2:]
elements = input_elem.copy()
elements.sort(key=lambda x: atomic_numbers[x])

print('# Define interatomic potential.\n')

if sys.argv[1]=='serial':
  # 7-net (serial)
  print('pair_style       e3gnn\n')
  print('pair_coeff       * * ${POT} ' + f'{" ".join(elements)}\n')
elif sys.argv[1]=='parallel':
  # 7-net (parallel)
  print('pair_style       e3gnn/parallel\n')
  print('pair_coeff       * * ${NMPL} ${POT} ' + f'{" ".join(elements)}\n')

for idx, elem in enumerate(input_elem):
    print(f"mass\t{idx+1} {atomic_masses[atomic_numbers[elem]]:.4f}")

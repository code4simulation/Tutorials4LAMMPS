from ase.data import atomic_masses, atomic_numbers
from scipy.constants import N_A, calorie, eV
from ase.io import read, write

def kcal_mol_to_eV(val):
  return val * 1000 * calorie / N_A / eV

init_str = 'Ar_FCC_init.lammps'
rlx_str = 'Ar_FCC_relaxed.lammps'

mass = atomic_masses[atomic_numbers['Ar']]
# J. Chem. Phys. 107, 4618 (1997)
#rcut = 10 # Angstrom
#When volume is larger than a specific value, potential energy shitfs up abruptly.
rcut = 20 # Angstrom
eps = 0.238122 #kcal/mol
eps_eV = kcal_mol_to_eV(eps)
sigma = 3.405 # Angstrom

# Minimization with box/relax
script=f'''# LAMMPS input script for Argon with Lennard-Jones potential
units       metal
atom_style  atomic

# Create a box and atoms
lattice         fcc {5.30}
region          box block 0 1 0 1 0 1
create_box      1 box
create_atoms    1 box

mass 1 {mass:.4f}

write_data      {init_str}

# Define the Lennard-Jones potential
pair_style      lj/cut {rcut}
pair_coeff      1 1 {eps_eV} {sigma} {rcut}

# Define thermodynamic output
thermo          10
thermo_style    custom step temp pe etotal press vol

fix 1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0
min_style cg
minimize 1.0e-4 1.0e-6 100 1000
write_data      {rlx_str}
'''
with open('script_min.in','w') as o:
    o.write(script)


# Equation of states
script=f'''# LAMMPS input script for Argon with Lennard-Jones potential
units       metal
atom_style  atomic

read_data   {rlx_str}

mass 1 {mass:.4f}

# Define the Lennard-Jones potential
pair_style      lj/cut {rcut}
pair_coeff      1 1 {eps_eV} {sigma} {rcut}

# Define thermodynamic output
thermo          10
thermo_style    custom step temp pe etotal press vol

#
variable vmin equal 0.97
variable vmax equal 1.03
variable vstep equal 0.005
variable vloop equal (${{vmax}}-${{vmin}})/${{vstep}}+1

variable i loop ${{vloop}}
label loop

if "${{i}} > 1" then &
  "variable scale equal (${{vmin}}+(${{i}}-1)*${{vstep}})/(${{vmin}}+(${{i}}-2)*${{vstep}})" &
else &
  "variable scale equal ${{vmin}}+(${{i}}-1)*${{vstep}}"

variable ratio equal ${{vmin}}+(${{i}}-1)*${{vstep}}
change_box all x scale ${{scale}} y scale ${{scale}} z scale ${{scale}} remap

minimize 1.0e-4 1.0e-6 100 1000

variable E equal pe
variable V equal vol
print "${{ratio}} ${{E}} ${{V}}" append eos.dat

next i
jump SELF loop

# 종료
write_data final.data
'''
with open('script_eos.in','w') as o:
    o.write(script)

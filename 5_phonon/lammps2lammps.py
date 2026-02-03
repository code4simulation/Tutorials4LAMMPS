import sys
from ase.io import read
from ase.io.lammpsdata import write_lammps_data

lammps = read(sys.argv[1], style='atomic',format='lammps-data')
lammps = lammps[lammps.get_array('id')-1]

#lammps = lammps.repeat([3,3,3])
#lammps.center(vacuum=10,axis=2)

write_lammps_data(f'{sys.argv[2]}/unit.lammps', lammps, force_skew=True)

atom_type_labels='Atom Type Labels\n\n'
for i, e in enumerate(sys.argv[3:]):
    atom_type_labels += f'{i+1:>6d} {e:>4s}\n'

file_path = f'{sys.argv[2]}/unit.lammps'
with open(file_path, 'r') as f:
    lines = f.readlines()

insert_index = -1
for i, line in enumerate(lines):
    if "Atoms # atomic" in line:
        insert_index = i
        break

if insert_index != -1:
    lines.insert(insert_index, '\n')
    lines.insert(insert_index, atom_type_labels)

    with open(file_path, 'w') as f:
        f.writelines(lines)

from pyscf import gto, scf
import numpy as np
from pyscf.symm import rotation_mat


mat = []
with open('ExMol5.zmat', 'r') as file:
    z_mat = file.read()

atom_list = []
int_coordinates = np.zeros((23,9))
j = 0
for line in z_mat.splitlines():

    line = line.strip()
    rawd = line.split()
    for i in range(len(rawd)):

        if i == 0:
            atom_list.append(rawd[0])
        else:
            int_coordinates[j, i-1] =rawd[i]
    j = j + 1


print(int_coordinates)
print(atom_list)

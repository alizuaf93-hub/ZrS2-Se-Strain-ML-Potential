import os
from ase.io import read, write
from ase.build import make_supercell
import numpy as np

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
relaxed_path = os.path.join(base_dir, "02_relaxation", "ZrS2_relaxed.cif")
output_dir = os.path.join(base_dir, "03_supercell")

os.makedirs(output_dir, exist_ok=True)

# Read relaxed structure
atoms = read(relaxed_path)

# 2x2x1 transformation matrix
P = np.array([[2,0,0],
              [0,2,0],
              [0,0,1]])

supercell = make_supercell(atoms, P)

print("Supercell created")
print("Number of atoms:", len(supercell))

# Save
write(os.path.join(output_dir, "ZrS2_2x2x1.cif"), supercell)

print("Saved in 03_supercell folder")
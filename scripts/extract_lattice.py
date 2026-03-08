import os
from ase.io import read

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
relaxed_path = os.path.join(base_dir, "02_relaxation", "ZrS2_relaxed.cif")

atoms = read(relaxed_path)

cell = atoms.get_cell()
a = cell.lengths()[0]
b = cell.lengths()[1]
c = cell.lengths()[2]

print("\nRelaxed Lattice Parameters:")
print("a =", round(a, 4), "Angstrom")
print("b =", round(b, 4), "Angstrom")
print("c =", round(c, 4), "Angstrom")

print("\nNumber of atoms:", len(atoms))
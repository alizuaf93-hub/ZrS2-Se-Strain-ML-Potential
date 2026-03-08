import os
from ase.io import read, write
from ase.optimize import BFGS
from mattersim.forcefield import MatterSimCalculator

# -----------------------------
# Paths
# -----------------------------

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
structure_path = os.path.join(base_dir, "01_structure", "ZrS2_primitive.cif")
output_dir = os.path.join(base_dir, "02_relaxation")

os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# Read structure
# -----------------------------

atoms = read(structure_path)

print("Initial number of atoms:", len(atoms))
print("Initial cell:")
print(atoms.get_cell())

# -----------------------------
# Attach MatterSim calculator
# -----------------------------

atoms.calc = MatterSimCalculator()

# -----------------------------
# Geometry optimization
# -----------------------------

optimizer = BFGS(atoms, trajectory=os.path.join(output_dir, "relax.traj"))
optimizer.run(fmax=0.02)

# -----------------------------
# Final results
# -----------------------------

final_energy = atoms.get_potential_energy()

print("\nRelaxation completed")
print("Final energy (eV):", final_energy)

# Save relaxed structure
relaxed_path = os.path.join(output_dir, "ZrS2_relaxed.cif")
write(relaxed_path, atoms)

print("Relaxed structure saved at:", relaxed_path)
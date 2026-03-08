import os
from ase.io import read, write
from ase.optimize import BFGS
from mattersim.forcefield import MatterSimCalculator

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
supercell_path = os.path.join(base_dir, "03_supercell", "ZrS2_2x2x1.cif")
output_dir = os.path.join(base_dir, "03_supercell")

# Read structure
atoms = read(supercell_path)

print("Initial number of atoms:", len(atoms))

# Attach calculator
atoms.calc = MatterSimCalculator()

# Relaxation
optimizer = BFGS(atoms, trajectory=os.path.join(output_dir, "supercell_relax.traj"))
optimizer.run(fmax=0.02)

# Final energy
energy = atoms.get_potential_energy()

print("\nSupercell relaxation completed")
print("Final energy (eV):", energy)

# Save relaxed structure
write(os.path.join(output_dir, "ZrS2_2x2x1_relaxed.cif"), atoms)

print("Relaxed supercell saved")
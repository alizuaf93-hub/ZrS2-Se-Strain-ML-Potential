import os
from ase.io import read, write
from ase.optimize import BFGS
from mattersim.forcefield import MatterSimCalculator

# -----------------------------
# Paths
# -----------------------------

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
supercell_path = os.path.join(base_dir, "03_supercell", "ZrS2_2x2x1_relaxed.cif")
doping_dir = os.path.join(base_dir, "05_doping", "S_to_Se")

os.makedirs(doping_dir, exist_ok=True)

# -----------------------------
# Read supercell
# -----------------------------

atoms = read(supercell_path)

print("Initial composition:", atoms.get_chemical_formula())
print("Total atoms:", len(atoms))

# -----------------------------
# Replace first S with Se
# -----------------------------

for i, atom in enumerate(atoms):
    if atom.symbol == "S":
        atoms[i].symbol = "Se"
        print("Replaced atom index:", i)
        break

print("New composition:", atoms.get_chemical_formula())

# -----------------------------
# Attach calculator
# -----------------------------

atoms.calc = MatterSimCalculator()

# -----------------------------
# Relaxation
# -----------------------------

optimizer = BFGS(atoms, trajectory=os.path.join(doping_dir, "relax.traj"))
optimizer.run(fmax=0.02)

energy = atoms.get_potential_energy()
energy_per_atom = energy / len(atoms)

print("\nDoped structure relaxed")
print("Total Energy (eV):", energy)
print("Energy per atom (eV):", energy_per_atom)

# Save structure
write(os.path.join(doping_dir, "ZrS2_S_to_Se_relaxed.cif"), atoms)

print("Relaxed doped structure saved.")
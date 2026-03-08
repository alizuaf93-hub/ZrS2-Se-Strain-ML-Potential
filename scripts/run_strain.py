import os
import numpy as np
import pandas as pd
from ase.io import read, write
from ase.optimize import BFGS
from mattersim.forcefield import MatterSimCalculator

# -----------------------------
# Base paths
# -----------------------------

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
supercell_path = os.path.join(base_dir, "03_supercell", "ZrS2_2x2x1_relaxed.cif")
strain_base_dir = os.path.join(base_dir, "04_strain_study")
results_dir = os.path.join(base_dir, "06_results")

os.makedirs(results_dir, exist_ok=True)

# Strain values
strain_values = {
    "strain_m4": -4,
    "strain_m2": -2,
    "strain_0": 0,
    "strain_p2": 2,
    "strain_p4": 4
}

results = []

# -----------------------------
# Loop over strain cases
# -----------------------------

for folder, strain_percent in strain_values.items():

    print("\nRunning case:", folder)

    strain = strain_percent / 100.0

    atoms = read(supercell_path)
    cell = atoms.get_cell()

    new_cell = cell.copy()
    new_cell[0] *= (1 + strain)
    new_cell[1] *= (1 + strain)

    atoms.set_cell(new_cell, scale_atoms=True)

    atoms.calc = MatterSimCalculator()

    case_dir = os.path.join(strain_base_dir, folder)

    optimizer = BFGS(atoms, trajectory=os.path.join(case_dir, "relax.traj"))
    optimizer.run(fmax=0.02)

    energy = atoms.get_potential_energy()
    energy_per_atom = energy / len(atoms)

    print("Final energy (eV):", energy)

    write(os.path.join(case_dir, "ZrS2_strained_relaxed.cif"), atoms)

    results.append([strain_percent, energy, energy_per_atom])

# -----------------------------
# Save CSV
# -----------------------------

df = pd.DataFrame(results, columns=["Strain (%)", "Total Energy (eV)", "Energy per Atom (eV)"])
df = df.sort_values("Strain (%)")

csv_path = os.path.join(results_dir, "strain_energy_results.csv")
df.to_csv(csv_path, index=False)

print("\nAll strain cases completed.")
print("Results saved at:", csv_path)
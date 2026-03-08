import os
import pandas as pd
from ase.io import read

# Base paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
strain_base_dir = os.path.join(base_dir, "04_strain_study")
output_dir = os.path.join(base_dir, "06_results")

os.makedirs(output_dir, exist_ok=True)

# Strain mapping
strain_values = {
    "strain_m4": -4,
    "strain_m2": -2,
    "strain_0": 0,
    "strain_p2": 2,
    "strain_p4": 4
}

data = []

for folder, strain in strain_values.items():
    structure_path = os.path.join(strain_base_dir, folder, "ZrS2_strained_relaxed.cif")
    atoms = read(structure_path)

    atoms.calc = None  # No recalculation
    energy = None

    # Energy manually insert karenge (since already printed)
    # For now leave placeholder

    data.append([strain, len(atoms)])

df = pd.DataFrame(data, columns=["Strain (%)", "Number of atoms"])

df.to_csv(os.path.join(output_dir, "strain_summary.csv"), index=False)

print("Basic strain summary saved.")
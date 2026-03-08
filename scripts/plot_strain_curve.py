import os
import pandas as pd
import matplotlib.pyplot as plt

# Base paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, "06_results")

csv_path = os.path.join(results_dir, "strain_energy_results.csv")

# Read data
df = pd.read_csv(csv_path)

# Sort just in case
df = df.sort_values("Strain (%)")

strain = df["Strain (%)"]
energy_per_atom = df["Energy per Atom (eV)"]

# Plot
plt.figure()
plt.plot(strain, energy_per_atom, marker='o')
plt.xlabel("Strain (%)")
plt.ylabel("Energy per Atom (eV)")
plt.title("Energy vs Biaxial Strain for ZrS2")
plt.grid(True)

# Save figure
fig_path = os.path.join(results_dir, "energy_vs_strain.png")
plt.savefig(fig_path, dpi=300)

plt.show()

print("Figure saved at:", fig_path)
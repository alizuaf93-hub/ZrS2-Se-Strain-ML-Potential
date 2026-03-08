import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 12,
    "font.family": "Arial",
    "axes.linewidth": 1.2
})

# Constants
mu_S = -1.1706586
mu_Se = -1.008476
mu_diff = mu_S - mu_Se

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, "06_results")

undoped_csv = os.path.join(results_dir, "strain_energy_results.csv")
doped_csv = os.path.join(results_dir, "doped_strain_energy_results.csv")

df_u = pd.read_csv(undoped_csv).sort_values("Strain (%)")
df_d = pd.read_csv(doped_csv).sort_values("Strain (%)")

strain = df_u["Strain (%)"].values
E_u = df_u["Total Energy (eV)"].values
E_d = df_d["Total Energy (eV)"].values

# Formation energy per supercell
E_form = E_d - E_u + mu_diff

# Per atom normalization (optional)
E_form_atom = E_form / 12.0

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12,5))

# Panel A – Total formation energy
axes[0].plot(strain, E_form, 'o-', color='darkred')
axes[0].set_title("A) Formation Energy vs Strain (Supercell)")
axes[0].set_xlabel("Strain (%)")
axes[0].set_ylabel("Formation Energy (eV)")
axes[0].grid(False)

# Panel B – Per atom
axes[1].plot(strain, E_form_atom, 's--', color='black')
axes[1].set_title("B) Formation Energy per Atom")
axes[1].set_xlabel("Strain (%)")
axes[1].set_ylabel("Formation Energy (eV/atom)")
axes[1].grid(False)

plt.tight_layout()

fig_path = os.path.join(results_dir, "Figure5_Formation_vs_Strain.png")
plt.savefig(fig_path, dpi=900, bbox_inches='tight')
plt.show()

print("Figure 5 saved at:", fig_path)
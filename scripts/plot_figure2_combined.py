import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 12,
    "font.family": "Arial",
    "axes.linewidth": 1.2
})

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, "06_results")

undoped_csv = os.path.join(results_dir, "strain_energy_results.csv")
doped_csv = os.path.join(results_dir, "doped_strain_energy_results.csv")

df_u = pd.read_csv(undoped_csv).sort_values("Strain (%)")
df_d = pd.read_csv(doped_csv).sort_values("Strain (%)")

strain_u = df_u["Strain (%)"].values
energy_u = df_u["Energy per Atom (eV)"].values
rel_u = energy_u - np.min(energy_u)

strain_d = df_d["Strain (%)"].values
energy_d = df_d["Energy per Atom (eV)"].values
rel_d = energy_d - np.min(energy_d)

# Quadratic fit
fit_u = np.polyfit(strain_u, energy_u, 2)
fit_d = np.polyfit(strain_d, energy_d, 2)

x_smooth = np.linspace(-4, 4, 300)
y_fit_u = np.polyval(fit_u, x_smooth)
y_fit_d = np.polyval(fit_d, x_smooth)

# Plot
fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

# ---------- Panel A: Relative Energy ----------
axes[0].plot(strain_u, rel_u, 'o-', color='black', label='Undoped')
axes[0].plot(strain_d, rel_d, 's--', color='blue', label='Se-Doped')
axes[0].set_ylabel("Relative Energy (eV)")
axes[0].set_title("A) Relative Energy vs Strain")
axes[0].legend(frameon=False)

# ---------- Panel B: Absolute Energy ----------
axes[1].plot(strain_u, energy_u, 'o', color='black')
axes[1].plot(x_smooth, y_fit_u, '-', color='black', linewidth=1.8, label='Undoped Fit')

axes[1].plot(strain_d, energy_d, 's', color='blue')
axes[1].plot(x_smooth, y_fit_d, '--', color='blue', linewidth=1.8, label='Se-Doped Fit')

axes[1].set_ylabel("Energy per Atom (eV)")
axes[1].set_xlabel("Strain (%)")
axes[1].set_title("B) Absolute Energy vs Strain")
axes[1].legend(frameon=False)

for ax in axes:
    ax.grid(False)

plt.tight_layout()

fig_path = os.path.join(results_dir, "Figure2_Combined_Final.png")
plt.savefig(fig_path, dpi=900, bbox_inches='tight')
plt.show()

print("Combined Figure 2 saved at:", fig_path)
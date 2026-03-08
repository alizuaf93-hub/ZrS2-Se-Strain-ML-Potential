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

strain = df_u["Strain (%)"].values
energy_u = df_u["Energy per Atom (eV)"].values
energy_d = df_d["Energy per Atom (eV)"].values

# Quadratic fits
fit_u = np.polyfit(strain, energy_u, 2)
fit_d = np.polyfit(strain, energy_d, 2)

x_smooth = np.linspace(-4, 4, 400)
y_fit_u = np.polyval(fit_u, x_smooth)
y_fit_d = np.polyval(fit_d, x_smooth)

curvature_u = 2 * fit_u[0]
curvature_d = 2 * fit_d[0]

fig, axes = plt.subplots(1, 2, figsize=(12,5))

# Panel A: Direct overlay zoomed
axes[0].plot(x_smooth, y_fit_u, color='black', linewidth=2, label='Undoped')
axes[0].plot(x_smooth, y_fit_d, color='blue', linestyle='--', linewidth=2, label='Se-Doped')
axes[0].set_xlim(-3,3)
axes[0].set_title("A) Zoomed Energy Curvature Comparison")
axes[0].set_xlabel("Strain (%)")
axes[0].set_ylabel("Energy per Atom (eV)")
axes[0].legend(frameon=False)
axes[0].grid(False)

# Panel B: Second derivative comparison
axes[1].plot(x_smooth, np.full_like(x_smooth, curvature_u),
             color='black', linewidth=2, label='Undoped')
axes[1].plot(x_smooth, np.full_like(x_smooth, curvature_d),
             color='blue', linestyle='--', linewidth=2, label='Se-Doped')
axes[1].set_title("B) Second Derivative (Elastic Indicator)")
axes[1].set_xlabel("Strain (%)")
axes[1].set_ylabel("d²E/dε²")
axes[1].legend(frameon=False)
axes[1].grid(False)

plt.tight_layout()

fig_path = os.path.join(results_dir, "Figure4_Upgraded.png")
plt.savefig(fig_path, dpi=900, bbox_inches='tight')
plt.show()

print("Upgraded Figure 4 saved at:", fig_path)
print("Second derivative Undoped:", curvature_u)
print("Second derivative Doped:", curvature_d)
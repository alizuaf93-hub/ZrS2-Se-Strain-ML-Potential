import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams.update({
    "font.size": 12,
    "font.family": "Arial",
    "axes.linewidth": 1.2
})

# -------- Paths --------
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, "06_results")

undoped_csv = os.path.join(results_dir, "strain_energy_results.csv")
doped_csv = os.path.join(results_dir, "doped_strain_energy_results.csv")

df_u = pd.read_csv(undoped_csv).sort_values("Strain (%)")
df_d = pd.read_csv(doped_csv).sort_values("Strain (%)")

strain_u = df_u["Strain (%)"].values
energy_u = df_u["Energy per Atom (eV)"].values

strain_d = df_d["Strain (%)"].values
energy_d = df_d["Energy per Atom (eV)"].values

# Quadratic fitting
fit_u = np.polyfit(strain_u, energy_u, 2)
fit_d = np.polyfit(strain_d, energy_d, 2)

x_smooth = np.linspace(-4, 4, 400)
y_fit_u = np.polyval(fit_u, x_smooth)
y_fit_d = np.polyval(fit_d, x_smooth)

# Mesh grid
Y_range = np.linspace(min(min(y_fit_u), min(y_fit_d)),
                      max(max(y_fit_u), max(y_fit_d)), 400)

X_u, Y_u = np.meshgrid(x_smooth, Y_range)
Z_u = np.abs(Y_u - np.polyval(fit_u, X_u))

X_d, Y_d = np.meshgrid(x_smooth, Y_range)
Z_d = np.abs(Y_d - np.polyval(fit_d, X_d))

# Common normalization
norm = Normalize(vmin=0, vmax=max(Z_u.max(), Z_d.max()))

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

cont1 = axes[0].contourf(X_u, Y_u, Z_u, levels=60, cmap='plasma', norm=norm)
axes[0].plot(x_smooth, y_fit_u, color='white', linewidth=2.5)
axes[0].set_title("A) Undoped Energy Landscape")
axes[0].set_xlabel("Strain (%)")
axes[0].set_ylabel("Energy per Atom (eV)")

cont2 = axes[1].contourf(X_d, Y_d, Z_d, levels=60, cmap='plasma', norm=norm)
axes[1].plot(x_smooth, y_fit_d, color='white', linewidth=2.5)
axes[1].set_title("B) Se-Substituted Energy Landscape")
axes[1].set_xlabel("Strain (%)")

for ax in axes:
    ax.grid(False)

# Single shared colorbar
cbar = fig.colorbar(cont1, ax=axes.ravel().tolist(), shrink=0.85)
cbar.set_label("Energy Density (eV)")

plt.tight_layout()

fig_path = os.path.join(results_dir, "Figure3_Final.png")
plt.savefig(fig_path, dpi=900, bbox_inches='tight')
plt.show()

print("Final Figure 3 saved at:", fig_path)
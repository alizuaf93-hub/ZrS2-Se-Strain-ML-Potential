import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.neighborlist import neighbor_list
from scipy.stats import gaussian_kde

plt.rcParams.update({
    "font.size": 12,
    "font.family": "Arial",
    "axes.linewidth": 1.2
})

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

undoped_path = os.path.join(base_dir, "03_supercell", "ZrS2_2x2x1_relaxed.cif")
doped_path = os.path.join(base_dir, "05_doping", "S_to_Se", "ZrS2_S_to_Se_relaxed.cif")

atoms_u = read(undoped_path)
atoms_d = read(doped_path)

cutoff = 3.0

# ---- Extract Bonds ----
i_u, j_u, d_u = neighbor_list("ijd", atoms_u, cutoff)
bond_u = [d for i,j,d in zip(i_u,j_u,d_u)
          if atoms_u[i].symbol=="Zr" and atoms_u[j].symbol=="S"]

i_d, j_d, d_d = neighbor_list("ijd", atoms_d, cutoff)

bond_zr_s = []
bond_zr_se = []

for i,j,d in zip(i_d,j_d,d_d):
    if atoms_d[i].symbol=="Zr" and atoms_d[j].symbol=="S":
        bond_zr_s.append(d)
    if atoms_d[i].symbol=="Zr" and atoms_d[j].symbol=="Se":
        bond_zr_se.append(d)

bond_u = np.array(bond_u)
bond_zr_s = np.array(bond_zr_s)
bond_zr_se = np.array(bond_zr_se)

# ---- KDE ----
x_range = np.linspace(2.55, 2.75, 400)

kde_u = gaussian_kde(bond_u)
kde_s = gaussian_kde(bond_zr_s)
kde_se = gaussian_kde(bond_zr_se)

# ---- Plot ----
fig, axes = plt.subplots(1,2, figsize=(12,5))

# Panel A – Smooth Density
axes[0].plot(x_range, kde_u(x_range), color='black', linewidth=2, label="Undoped Zr–S")
axes[0].plot(x_range, kde_s(x_range), color='blue', linewidth=2, label="Doped Zr–S")
axes[0].plot(x_range, kde_se(x_range), color='red', linewidth=2, label="Zr–Se")

axes[0].axvline(np.mean(bond_u), color='black', linestyle='--', alpha=0.7)
axes[0].axvline(np.mean(bond_zr_s), color='blue', linestyle='--', alpha=0.7)
axes[0].axvline(np.mean(bond_zr_se), color='red', linestyle='--', alpha=0.7)

axes[0].set_title("A) Bond Length Density Distribution")
axes[0].set_xlabel("Bond Length (Å)")
axes[0].set_ylabel("Density")
axes[0].legend(frameon=False)
axes[0].grid(False)

# Panel B – Clean Boxplot
box = axes[1].boxplot([bond_u, bond_zr_s, bond_zr_se],
                      labels=["Undoped Zr–S","Doped Zr–S","Zr–Se"],
                      patch_artist=True)

colors = ['black','blue','red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.5)

axes[1].set_title("B) Statistical Bond Comparison")
axes[1].set_ylabel("Bond Length (Å)")
axes[1].grid(False)

plt.tight_layout()

results_dir = os.path.join(base_dir, "06_results")
fig_path = os.path.join(results_dir, "Figure6_Bond_Distribution_Polished.png")
plt.savefig(fig_path, dpi=900, bbox_inches='tight')
plt.show()

print("Polished Figure saved at:", fig_path)
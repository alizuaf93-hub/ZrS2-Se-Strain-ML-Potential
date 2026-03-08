import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

plt.rcParams.update({
    "font.size": 13,
    "font.family": "Arial",
    "axes.linewidth": 1.4
})

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
undoped_path = os.path.join(base_dir, "03_supercell", "ZrS2_2x2x1_relaxed.cif")
doped_path = os.path.join(base_dir, "05_doping", "S_to_Se", "ZrS2_S_to_Se_relaxed.cif")

atoms_u = read(undoped_path)
atoms_d = read(doped_path)

# Reference S in undoped
s_index = [i for i,a in enumerate(atoms_u) if a.symbol=="S"][0]

# Se in doped
se_index = [i for i,a in enumerate(atoms_d) if a.symbol=="Se"][0]

# Distances from S (undoped)
dist_u = []
for i,a in enumerate(atoms_u):
    if a.symbol == "Zr":
        dist_u.append(atoms_u.get_distance(s_index, i, mic=True))

# Distances from Se (doped)
dist_d = []
for i,a in enumerate(atoms_d):
    if a.symbol == "Zr":
        dist_d.append(atoms_d.get_distance(se_index, i, mic=True))

dist_u = np.sort(np.array(dist_u))
dist_d = np.sort(np.array(dist_d))

min_len = min(len(dist_u), len(dist_d))
dist_u = dist_u[:min_len]
dist_d = dist_d[:min_len]

delta_d = dist_d - dist_u

neighbor_order = np.arange(1, min_len+1)

# ---- Plot ----
fig, axes = plt.subplots(1,2, figsize=(13,5))

# Panel A – Distance comparison
axes[0].plot(neighbor_order, dist_u,
             marker='o', markersize=7,
             linewidth=2.5, color='black',
             label='Undoped (S reference)')

axes[0].plot(neighbor_order, dist_d,
             marker='s', markersize=7,
             linewidth=2.5, linestyle='--',
             color='#1f77b4',
             label='Se-doped (Se reference)')

axes[0].set_title("A) Shell-Matched Distance Comparison")
axes[0].set_xlabel("Coordination Shell")
axes[0].set_ylabel("Distance (Å)")
axes[0].set_xticks(neighbor_order)
axes[0].legend(frameon=False)
axes[0].grid(False)

# Panel B – Distortion decay
axes[1].plot(neighbor_order, delta_d,
             marker='D', markersize=7,
             linewidth=2.5,
             color='#8B0000')

axes[1].axhline(0, color='black', linewidth=1.2)

axes[1].set_title("B) True Local Distortion Decay")
axes[1].set_xlabel("Coordination Shell")
axes[1].set_ylabel("Δd (Å)")
axes[1].set_xticks(neighbor_order)
axes[1].set_ylim(0, 0.15)
axes[1].grid(False)

plt.tight_layout()

results_dir = os.path.join(base_dir, "06_results")
fig_path = os.path.join(results_dir, "Figure9_Distortion_Decay_Ultra.png")
plt.savefig(fig_path, dpi=1200, bbox_inches='tight')
plt.show()

print("Ultra-polished distortion figure saved at:", fig_path)
print("Distortion values:", delta_d)
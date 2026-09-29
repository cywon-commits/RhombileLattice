"""Figure for the corrected finding: the fully-frustrated triangular
lattice's D3 transition is a SHARP kink at D3=3.0 (not the broad
D3~2-3.5 crossover originally reported), confirmed by heavy SA matching
the exact two-line envelope E*(D3)=min(144*D3, 432) at every checked
point. Shows the original (light-scan, under-searched) values alongside
the corrected heavy-search values against the exact trivial baseline.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

N = 432
N_BONDS = 1296

# Original light scan (triangular_lattice_D3_sweep.py, 5 random + 3 seeded
# restarts): systematically above the true optimum, badly so at D3=2.75.
D3_orig = np.array([1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5])
n3_orig = np.array([1/3, 0.3241, 0.3171, 0.3079, 0.322, 0.0949, 0.0])
nv_orig = np.array([0.0, 0.0093, 0.0162, 0.0255, 0.116, 0.2384, 0.3333])
E_orig = n3_orig * N * D3_orig + nv_orig * N_BONDS

# Heavy re-verification (triangular_lattice_sharp_transition_check.py +
# the D3=2.75 standalone reverify): exactly matches the trivial two-line
# envelope everywhere checked.
D3_heavy = np.array([2.0, 2.25, 2.5, 2.75, 3.0, 3.25])
n3_heavy = np.array([1/3, 1/3, 1/3, 1/3, 0.0856, 0.0])
nv_heavy = np.array([0.0, 0.0, 0.0, 0.0, 0.2477, 1/3])
E_heavy = n3_heavy * N * D3_heavy + nv_heavy * N_BONDS

D3_line = np.linspace(0, 4.0, 400)
E_trivial = np.minimum(N / 3 * D3_line, N_BONDS / 3)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

ax = axes[0]
ax.plot(D3_line, E_trivial, color="#3b4ba8", lw=2, zorder=1,
        label="exact baseline: min(144·D3, 432)")
ax.scatter(D3_orig, E_orig, marker="x", s=70, color="#c5791f", zorder=3,
           label="original light scan (under-searched)")
ax.scatter(D3_heavy, E_heavy, marker="o", s=60, facecolor="none",
           edgecolor="#3d7a4f", linewidth=2, zorder=4,
           label="heavy re-verification (matches exactly)")
ax.axvline(3.0, color="gray", ls=":", lw=1, zorder=0)
ax.set_xlabel("D3")
ax.set_ylabel("Ground-state energy E")
ax.set_title("Energy: sharp kink at D3=3.0, not a broad crossover")
ax.legend(fontsize=8, loc="lower right")
ax.set_xlim(1.5, 4.0)
ax.set_ylim(200, 570)

ax2 = axes[1]
n3_line = np.where(D3_line < 3.0, 1/3, 0.0)
ax2.plot(D3_line, n3_line, color="#3b4ba8", lw=2, zorder=1,
         label="exact baseline")
ax2.scatter(D3_orig, n3_orig, marker="x", s=70, color="#c5791f", zorder=3,
            label="original light scan")
ax2.scatter(D3_heavy, n3_heavy, marker="o", s=60, facecolor="none",
            edgecolor="#3d7a4f", linewidth=2, zorder=4,
            label="heavy re-verification")
ax2.axvline(3.0, color="gray", ls=":", lw=1, zorder=0)
ax2.set_xlabel("D3")
ax2.set_ylabel("n_state3 / N")
ax2.set_title("State-3 fraction: 1/3 exactly, then drops to 0 at D3=3")
ax2.legend(fontsize=8, loc="upper right")
ax2.set_xlim(1.5, 4.0)

fig.tight_layout()
fig.savefig("triangular_lattice_sharp_kink_figure.png", dpi=150)
print("saved triangular_lattice_sharp_kink_figure.png")

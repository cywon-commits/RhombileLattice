"""Combined D2 (-1.2 to +1.6) vs n_state3 plot from the fine, rigorous
re-verification sweeps, replacing the original coarse-grid picture that
suggested smooth gradual ramps on both sides. Also plots E(D2) for
context. Pulls in the original coarse sweep's own D2>=1.1 tail (already
deep in the saturated regime, consistent with the fine sweep's D2=1.0
point) since the fine sweeps stopped at +1.0/-1.2.
"""
import pickle

import matplotlib.pyplot as plt

with open("D2_sweep_fine_negative_results.pkl", "rb") as f:
    neg = pickle.load(f)
with open("D2_sweep_fine_positive_results.pkl", "rb") as f:
    pos1 = pickle.load(f)
with open("D2_sweep_fine_positive_part2_results.pkl", "rb") as f:
    pos2 = pickle.load(f)
with open("D2_sweep_positive_results.pkl", "rb") as f:
    pos_coarse = pickle.load(f)

combined = {}
combined.update(neg)
combined.update(pos1)
combined.update(pos2)
for D2, row in pos_coarse.items():
    if D2 > max(pos2.keys()):
        combined[D2] = row

D2_vals = sorted(combined.keys())
n3_vals = [combined[d]["n3"] for d in D2_vals]
nv_vals = [combined[d]["n_violated"] for d in D2_vals]
E_vals = [combined[d]["E"] for d in D2_vals]

fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
ax_n3, ax_e = axes

ax_n3.plot(D2_vals, n3_vals, "o-", color="tab:orange", label="n_state3")
ax_n3.plot(D2_vals, nv_vals, "s--", color="tab:purple", alpha=0.7, label="n_violated")
ax_n3.axvline(0, color="gray", lw=0.7)
ax_n3.axvline(1.0, color="tab:red", lw=0.7, ls=":", label="D2=+D3")
ax_n3.axvline(-1.0, color="tab:red", lw=0.7, ls=":")
ax_n3.set_ylabel("count")
ax_n3.set_title("D1=0, D3=1 fixed, Case2 hexagon-detour: n_state3 vs D2 "
                 "(fine, rigorous re-verification)")
ax_n3.legend()

ax_e.plot(D2_vals, E_vals, "o-", color="tab:blue")
ax_e.axvline(0, color="gray", lw=0.7)
ax_e.axvline(1.0, color="tab:red", lw=0.7, ls=":")
ax_e.axvline(-1.0, color="tab:red", lw=0.7, ls=":")
ax_e.set_xlabel("D2")
ax_e.set_ylabel("ground-state energy")

plt.tight_layout()
fig.savefig("D2_sweep_fine_combined.png", dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved D2_sweep_fine_combined.png")
print(f"D2 range covered: {min(D2_vals)} to {max(D2_vals)}, {len(D2_vals)} points")

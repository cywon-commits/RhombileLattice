"""Figure for (b): exact E_inf over loop-updated dimer coverings with the
defects held fixed (results/B1_pair_h*.log, results/B2_p*.log)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path, col):
    return np.array([int(ln.split()[col]) for ln in open(path) if ln[0].isdigit()])


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax = axes[0]
for h, color in [(3, "#3d7a4f"), (6, "#3b4ba8"), (10, "#a83b3b")]:
    path = f"results/B1_pair_h{h}.log"
    d = next(ln for ln in open(path) if ln.startswith("# pair")).split("=")[-1].strip()
    e = load(path, 2)
    ax.plot(e, "o-", ms=3.5, lw=1, color=color,
            label=f"d_tri={d}: E_inf " + "/".join(map(str, sorted(set(e)))))
ax.set_xlabel("covering index (closed-loop worm updates)")
ax.set_ylabel("exact E_inf (violated bonds)")
ax.set_title("one defect pair, 41 coverings each")
ax.legend(fontsize=8)

ax = axes[1]
sa = {"0.04": 110, "0.08": 163, "0.12": 227}
for p, color in [("0.04", "#3d7a4f"), ("0.08", "#3b4ba8"), ("0.12", "#a83b3b")]:
    path = f"results/B2_p{p}.log"
    if not os.path.exists(path) or not any(ln[0].isdigit() for ln in open(path)):
        continue
    e = load(path, 2)
    lo = load(path, 1)
    ax.plot(e / e[0], "o-", ms=3.5, lw=1, color=color,
            label=f"p={p}: E_inf {e.min()}-{e.max()}  (SA <= {sa[p]})")
ax.axhline(1, color="k", lw=0.6, ls=":")
ax.set_xlabel("covering index (defects frozen)")
ax.set_ylabel("E_inf / E_inf(covering 0)")
ax.set_title("many defects (worm-built, 16x16): exact MILP")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("b_path_independence_figure.png", dpi=150)
print("saved b_path_independence_figure.png")

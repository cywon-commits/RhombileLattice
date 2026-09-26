"""Time series and distributions of |psi|^2 and E from annealed_hysteresis.py:
telegraph-like switching between an ordered and a disordered state, with a
bimodal distribution, signals two-phase coexistence (first-order); a single
broad distribution with both starts converging signals slow but continuous
ordering."""
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

cases = [(48, 0.515), (48, 0.525), (48, 0.535), (64, 0.515), (64, 0.525)]
fig, ax = plt.subplots(2, len(cases), figsize=(4 * len(cases), 7))
for j, (L, T) in enumerate(cases):
    allp, allE = [], []
    for st, c in (("ordered", "#3b4ba8"), ("random", "#a83b3b")):
        f = glob.glob(f"results/fss/hyst_L{L}_T{T:g}_{st}_s1.npz")
        if not f:
            continue
        d = np.load(f[0])
        t = np.arange(len(d["p2"])) * int(d["every"])
        ax[0, j].plot(t / 1000, d["p2"], lw=0.4, color=c, label=f"{st} start")
        allp.append(d["p2"][len(d["p2"]) // 10:])
        allE.append(d["E"][len(d["E"]) // 10:] / int(d["N_tri"]))
    p = np.concatenate(allp)
    ax[1, j].hist(p, bins=60, density=True, color="#6b6558")
    ax[0, j].set_title(f"L={L}, T={T}", fontsize=10)
    ax[0, j].set_xlabel("sweeps / 1000")
    ax[0, j].set_ylabel("|psi|^2")
    ax[0, j].legend(fontsize=7)
    ax[1, j].set_xlabel("|psi|^2 (both starts, first 10% dropped)")
    ax[1, j].set_ylabel("density")
fig.tight_layout()
fig.savefig("annealed_hysteresis_figure.png", dpi=130)
print("saved annealed_hysteresis_figure.png")

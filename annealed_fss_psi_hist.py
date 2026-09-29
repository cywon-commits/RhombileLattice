"""Order-parameter distribution at the pseudo-critical temperature (C peak)
for the largest sizes: areal density P(|psi|)/|psi| in the complex psi
plane. A single ring = continuous-transition-like; a central (disordered)
plateau coexisting with an ordered ring = incipient phase coexistence."""
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TPK = {24: 0.5406, 32: 0.5340, 48: 0.5284}          # C-peak temperatures (annealed_fss_analysis.py)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for (L, t), c in zip(TPK.items(), ["#a83b3b", "#6b4ba8", "#201d18"]):
    xs, pts = [], []
    for f in glob.glob(f"results/fss/fss_L{L}_d21_s*.npz"):
        d = np.load(f)
        k = np.argmin(abs(1 / d["betas"] - t))
        if abs(1 / d["betas"][k] - t) < 0.003:
            xs.append(np.hypot(d["re"][:, k], d["im"][:, k]))
            pts.append((d["re"][:, k], d["im"][:, k]))
    x = np.concatenate(xs)
    h, e = np.histogram(x, bins=25, range=(0, 0.7), density=True)
    r = 0.5 * (e[1:] + e[:-1])
    ax[0].plot(r, h / r / (h / r)[:6].max(), "o-", ms=3, color=c, label=f"L={L}, T={t}")
    if L == 48:
        re, im = np.concatenate([p[0] for p in pts]), np.concatenate([p[1] for p in pts])
        ax[1].hist2d(re, im, bins=60, range=[[-0.45, 0.45], [-0.45, 0.45]], cmap="Greys")
        ax[1].set_aspect("equal")
        ax[1].set_title(f"L=48, T={t}: psi in the complex plane", fontsize=10)
ax[0].set_xlabel("|psi|")
ax[0].set_ylabel("areal density P(|psi|)/|psi| (normalised)")
ax[0].legend(fontsize=8)
ax[1].set_xlabel("Re psi")
ax[1].set_ylabel("Im psi")
fig.tight_layout()
fig.savefig("annealed_fss_psi_hist.png", dpi=150)
print("saved annealed_fss_psi_hist.png")

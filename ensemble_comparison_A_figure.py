"""Fits and figure for ensemble_comparison_A.py results (results/A*.log).

A1 dimer: slope of F(r) against ln r over an intermediate window (small r
excluded for lattice-scale effects, r beyond ~half the torus diameter
excluded for finite-size saturation). Prediction 0.5.

A2 potts: F(L) against L for both seeds; compares a pure-linear fit with
a pure-log fit (same number of parameters) over L >= 6.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    rows = [ln.split() for ln in open(path) if ln[0].isdigit()]
    return np.array([[float(x) for x in r] for r in rows])


def fit_line(x, y):
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    return coef, float(np.sqrt(np.mean(resid ** 2)))


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

ax = axes[0]
for path, nx, color, rmax in [("results/A1_dimer_NX24.log", 24, "#c5791f", 23),
                              ("results/A1_dimer_NX32.log", 32, "#3b4ba8", 31)]:
    d = load(path)
    r, F = d[:, 0], d[:, 3]
    m = (r >= 5) & (r <= rmax)
    (a, k), rms = fit_line(np.log(r[m]), F[m])
    print(f"A1 NX={nx}: F = {a:.3f} + {k:.3f} ln r  (r in [5,{rmax}], rms {rms:.4f})")
    ax.plot(np.log(r), F - a, "o", ms=4, color=color, label=f"dimer, {nx}x{nx} torus (fit slope {k:.2f})")
lr = np.linspace(0, 4.2, 50)
ax.plot(lr, 0.5 * lr, "k--", lw=1.2, label="prediction 0.5 ln r")
ax.set_xlabel("ln r   (r = hop distance between monomers)")
ax.set_ylabel("F(r) = -ln[P(r)/g(r)]  (shifted)")
ax.set_title("A1: bonds fluctuate (dimer ensemble) - logarithmic")
ax.legend(fontsize=8, loc="upper left")

ax = axes[1]
for path, color, lab in [("results/A2_potts_s1.log", "#3d7a4f", "seed 1"),
                         ("results/A2_potts_s2.log", "#a83b3b", "seed 2")]:
    d = load(path)
    L, F = d[:, 0], d[:, 2]
    m = L >= 6
    (a, s), rms_lin = fit_line(L[m], F[m])
    (b, k), rms_log = fit_line(np.log(L[m]), F[m])
    print(f"A2 {lab}: linear F = {a:.3f} + {s:.4f} L (rms {rms_lin:.3f});  "
          f"log F = {b:.3f} + {k:.3f} ln L (rms {rms_log:.3f})")
    ax.plot(L, F, "o", ms=4, color=color, label=f"Potts, {lab} (slope {s:.3f}/hop)")
    ax.plot(L[m], a + s * L[m], "-", color=color, lw=1)
ax.set_xlabel("L = string length (hops)")
ax.set_ylabel("F(L) = -S(L) + const")
ax.set_title("A2: bonds fixed, D3=0 Potts spins - linear")
ax.legend(fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig("ensemble_comparison_A_figure.png", dpi=150)
print("saved ensemble_comparison_A_figure.png")

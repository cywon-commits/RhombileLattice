"""Analysis + figure for annealed_thermo.py (results/thermo/*.log).

Entropy per triangle by thermodynamic integration from beta=0:
  s(beta) = lnZ(0)/N_tri - int_0^beta e dbeta' + beta e(beta)
  lnZ(0)/N_tri = (N_site/N_tri) ln 3 + lnM/N_tri = 0.5 ln 3 + lnM/N_tri
with e(beta=0) exact: (1.5 - d_inf)/3 + D3/6 + mu (1 - 2 d_inf), d_inf =
<n_dimer>/N_tri of uniform matchings (matchings log). The beta grid is the
high-T supplement (L=12, size effects negligible there) followed by each
L's own scan. Consistency check: s(T->0) should approach the honeycomb
dimer (lozenge-tiling) entropy 0.16153 per triangle.

Usage: python3 annealed_thermo_analysis.py <D3> <mu> L1 L2 ...
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S_TILING = 0.323065947 / 2          # honeycomb dimer entropy per vertex (= per triangle)


def load(path):
    rows = [ln.split() for ln in open(path) if ln[0].isdigit()]
    return np.array([[float(x) for x in r] for r in rows])


def matchings_info(path="results/thermo/matchings_L16.log"):
    txt = open(path).read()
    lnM = float(re.search(r"lnM/N_tri = ([\d.]+)", txt).group(1))
    d = load(path)
    return lnM, d[-1, 1], d[-1, 2]            # lnM, dimers/tri and monomer frac at z=1


def entropy(D3, mu, main, high, lnM, d_inf):
    e0 = (1.5 - d_inf) / 3 + D3 / 6 + mu * (1 - 2 * d_inf)
    hb = high[:, 1][high[:, 1] < main[0, 1] - 1e-9]
    he = high[:, 2][high[:, 1] < main[0, 1] - 1e-9]
    beta = np.concatenate([[0.0], hb, main[:, 1]])
    e = np.concatenate([[e0], he, main[:, 2]])
    order = np.argsort(beta)
    beta, e = beta[order], e[order]
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(beta))])
    s = 0.5 * np.log(3) + lnM - integ + beta * e
    return beta, e, s


def main():
    D3, mu = float(sys.argv[1]), float(sys.argv[2])
    Ls = [int(x) for x in sys.argv[3:]]
    tag = f"D{D3:g}_mu{mu:g}"
    lnM, d_inf, mon_inf = matchings_info()
    high = load(f"results/thermo/highT_L12_{tag}.log")
    print(f"lnM/N_tri={lnM:.4f}  s(inf)={0.5 * np.log(3) + lnM:.4f}  n_mon(inf)={mon_inf:.4f}")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    colors = ["#c5791f", "#3d7a4f", "#3b4ba8", "#a83b3b", "#6b6558"]
    print(" L   T_peak   C_max/N_tri   s(T_min)   n_mon(T_min)")
    for L, col in zip(Ls, colors):
        d = load(f"results/thermo/scan_L{L}_{tag}.log")
        T, C, nm = d[:, 0], d[:, 3], d[:, 4]
        beta, e, s = entropy(D3, mu, d, high, lnM, d_inf)
        i = np.argmax(C)
        print(f"{L:3d}  {T[i]:.3f}   {C[i]:.4f}      {s[-1]:.4f}    {nm[-1]:.4f}")
        lab = f"L={L} ({6 * L * L} triangles)"
        ax[0, 0].plot(T, d[:, 2], "o-", ms=2.5, lw=1, color=col, label=lab)
        ax[0, 1].plot(T, C, "o-", ms=2.5, lw=1, color=col, label=lab)
        ax[1, 0].plot(T, nm, "o-", ms=2.5, lw=1, color=col, label=lab)
        m = beta > 0
        ax[1, 1].plot(1 / beta[m], s[m], "-", lw=1.4, color=col, label=lab)
    ax[0, 0].set_ylabel("E / N_tri")
    ax[0, 1].set_ylabel("C / N_tri")
    ax[1, 0].set_ylabel("frustrated-triangle (monomer) fraction")
    ax[1, 0].axhline(mon_inf, color="k", ls=":", lw=1, label=f"T=inf limit {mon_inf:.3f}")
    ax[1, 1].set_ylabel("entropy S / N_tri")
    ax[1, 1].axhline(S_TILING, color="k", ls="--", lw=1, label=f"lozenge tilings {S_TILING:.4f}")
    ax[1, 1].axhline(0.5 * np.log(3) + lnM, color="k", ls=":", lw=1, label="T=inf: ln3/2 + lnM")
    ax[1, 1].set_xscale("log")
    for a in ax.flat:
        a.set_xlabel("T")
        a.legend(fontsize=7)
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.set_xlim(0, 4)
    fig.suptitle(f"Annealed bonds + Potts spins, D3={D3:g}, mu={mu:g}")
    fig.tight_layout()
    out = f"annealed_thermo_{tag}_figure.png"
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()

"""D1=D2 versus 0=D1<D2 in the annealed-bond model (results/thermo/*).

D2=0: T=0 manifold = all lozenge tilings (s(0) = 0.1615 per triangle).
D2>0: T=0 state = the unique rhombile tiling (s(0) = 0), so the extra
0.1615 per triangle has to be released at finite T.

Usage: python3 annealed_thermo_compare.py <L> [<L2> ...]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from annealed_thermo_analysis import load, matchings_info, entropy, S_TILING

R = "results/thermo"
RUNS = [  # (D2, scan file suffix, high-T file, colour)
    (0.0, "_lowT", f"{R}/highT_L12_D10_mu0.log", "#3b4ba8"),
    (0.5, "_d20.5", f"{R}/highT_L12_D10_mu0_d20.5.log", "#c5791f"),
    (1.0, "_d21", f"{R}/highT_L12_D10_mu0_d21.log", "#a83b3b"),
]


def main():
    Ls = [int(x) for x in sys.argv[1:]]
    lnM, d_inf, mon_inf = matchings_info()
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.2))
    styles = ["-", "--", ":"]
    print(" D2   L   T_peak(s)            C_peaks            s(Tmin)  hub(Tmin)  n_mon(Tmin)  int C/T dT")
    for D2, suf, hi, col in RUNS:
        for L, ls in zip(Ls, styles):
            try:
                d = load(f"{R}/scan_L{L}_D10_mu0{suf}.log")
                high = load(hi)
            except (FileNotFoundError, IndexError):
                continue
            T, C, nm = d[:, 0], d[:, 3], d[:, 4]
            beta, e, s = entropy(10.0, 0.0, d, high, lnM, d_inf, D2=D2)
            # local maxima of C
            pk = [i for i in range(1, len(C) - 1) if C[i] >= C[i - 1] and C[i] >= C[i + 1] and C[i] > 0.1]
            dS = np.trapezoid(C[::-1] / T[::-1], T[::-1])
            print(f"{D2:4.1f} {L:3d}  {', '.join(f'{T[i]:.3f}' for i in pk):20s} "
                  f"{', '.join(f'{C[i]:.3f}' for i in pk):18s} {s[-1]:.4f}   {d[-1, 9]:.4f}     {nm[-1]:.4f}      {dS:.4f}")
            lab = f"D2={D2:g}, L={L}"
            ax[0, 0].plot(T, C, ls, color=col, lw=1.3, label=lab)
            ax[0, 1].plot(T, d[:, 2], ls, color=col, lw=1.3, label=lab)
            m = beta > 0
            ax[0, 2].plot(1 / beta[m], s[m], ls, color=col, lw=1.3, label=lab)
            ax[1, 0].plot(T, nm, ls, color=col, lw=1.3, label=lab)
            ax[1, 1].plot(T, d[:, 9], ls, color=col, lw=1.3, label=lab)
            ax[1, 2].plot(T, d[:, 8], ls, color=col, lw=1.3, label=lab)
    ax[0, 0].set_ylabel("C / N_tri")
    ax[0, 1].set_ylabel("E / N_tri")
    ax[0, 2].set_ylabel("S / N_tri")
    ax[0, 2].axhline(S_TILING, color="k", ls="--", lw=0.8, label="lozenge tilings 0.1615")
    ax[0, 2].axhline(0, color="k", lw=0.5)
    ax[0, 2].set_xscale("log")
    ax[1, 0].set_ylabel("frustrated-triangle fraction")
    ax[1, 1].set_ylabel("hub-type site fraction (all 6 bonds active)")
    ax[1, 1].axhline(1 / 3, color="k", ls=":", lw=0.8, label="rhombile tiling 1/3")
    ax[1, 2].set_ylabel("state-2 fraction")
    ax[1, 2].axhline(1 / 3, color="k", ls=":", lw=0.8)
    for a in ax.flat:
        a.set_xlabel("T")
        a.legend(fontsize=6.5)
    for a in (ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1], ax[1, 2]):
        a.set_xlim(0, 2.5)
    fig.suptitle("Annealed bonds + Potts spins (D3=10, mu=0): D1=D2 versus 0=D1<D2")
    fig.tight_layout()
    fig.savefig("annealed_thermo_D2_compare_figure.png", dpi=150)
    print("saved annealed_thermo_D2_compare_figure.png")


if __name__ == "__main__":
    main()

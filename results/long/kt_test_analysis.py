"""Potts vs KT test from the Binder cumulant of the spin order parameter psi_s.

At a 3-state Potts point the Binder curves cross with slope max dU4/dbeta ~ L^(1/nu), 1/nu = 1.2.
At a KT-type (locking) transition the curves merge instead of crossing and the slope grows only
logarithmically, so the apparent 1/nu drifts towards 0 with L. Compared at mu=0 and mu=10 with the
same pipeline (exact-level WHAM over the single-T long runs, spin order parameter).

Usage: python3 results/long/kt_test_analysis.py
"""
import sys

import numpy as np

sys.path.insert(0, ".")
import annealed_fss_analysis as fa  # noqa: E402
import annealed_long_analysis as la  # noqa: E402


def analyse(mu, Ls, order):
    la.ORDER, la.MU, la.DIR = order, mu, "results/long"
    la._cache.clear()
    res = {}
    for L in Ls:
        runs = la.load_blocks(1.0, L)
        if not runs:
            continue
        b = np.array(sorted({r["betas"][0] for r in runs}))
        g = np.linspace(b.min(), b.max(), 400)
        res[L] = (g, fa.curves(runs, g), 1 / b)
    Ls = sorted(res)
    print(f"mu={mu} order={order}: sizes {Ls}")
    for L in Ls:
        g, cv, Ts = res[L]
        i = np.argmax(cv["dU"])
        print(f"  L={L:3d}  T range {Ts.min():.3f}-{Ts.max():.3f} ({len(Ts)} T)  max dU/dbeta={cv['dU'][i]:7.2f} at T={1 / g[i]:.4f}"
              f"  U(max)={cv['U'].max():.3f}")
    for a, b in zip(Ls[:-1], Ls[1:]):
        Tx = fa.crossing(res[a][0], res[a][1]["U"], res[b][0], res[b][1]["U"])
        sa, sb = np.max(res[a][1]["dU"]), np.max(res[b][1]["dU"])
        print(f"  L={a}/{b}: Binder crossing T={Tx:.4f}   apparent 1/nu = {np.log(sb / sa) / np.log(b / a):.2f}")
    lg = np.log(Ls)
    s = np.array([np.max(res[L][1]["dU"]) for L in Ls])
    print(f"  fit over all L: 1/nu = {np.polyfit(lg, np.log(s), 1)[0]:.2f}   (3-state Potts 1.2; KT -> 0)")


def direct(mu, Ls, Tmin, Tmax):
    """No reweighting (T spacing too coarse for WHAM overlap at mu=10): U4 per (L, T) from the raw
    series, with time-block errors; slope from finite differences between neighbouring T."""
    import glob
    import re
    print(f"mu={mu} spin order, direct per-T values (no reweighting), T in [{Tmin}, {Tmax}]")
    tab = {}
    for f in glob.glob(f"results/long/long_L*_d21_mu{mu:g}_T*_s*.npz"):
        m = re.search(r"L(\d+)_d21_mu[\d.]+_T([\d.]+)_s(\d+)\.npz$", f)
        L, T = int(m[1]), float(m[2])
        if L not in Ls or not Tmin <= T <= Tmax:
            continue
        z = np.load(f)
        k = 100000 // int(z["thin"])
        p2 = z["sre"][k:, 0].astype(float) ** 2 + z["sim"][k:, 0].astype(float) ** 2
        n = len(p2) // 8 * 8
        bl = p2[:n].reshape(8, -1)
        Ub = 1 - (bl ** 2).mean(1) / (2 * bl.mean(1) ** 2)
        tab.setdefault(L, {})[T] = (1 - (p2 ** 2).mean() / (2 * p2.mean() ** 2), Ub.std(ddof=1) / np.sqrt(8))
    Ts = sorted({T for L in tab for T in tab[L]})
    print("  T      " + "  ".join(f"L={L:<9d}" for L in sorted(tab)))
    for T in Ts:
        print(f"  {T:.3f}  " + "  ".join(f"{tab[L][T][0]:.3f}({1e3 * tab[L][T][1]:3.0f})" if T in tab[L] else " " * 11
                                       for L in sorted(tab)))
    return tab


if __name__ == "__main__":
    analyse(0.0, [32, 48, 64], "hub")
    direct(10.0, [24, 32, 48, 64, 96], 1.0, 1.2)

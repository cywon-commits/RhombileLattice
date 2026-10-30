"""Part IX task 2b: T_c(D2) at small D2, annealed model (pure Ising sweep, mu=0)
vs the fixed-bond TAFM, from Binder crossings of the SPIN order parameter psi_s.

For every (D2, model, L, T) the U4 of psi_s is measured directly (seeds pooled,
jackknife over 8 time blocks); the annealed runs are results/long/*_ising.npz
(sre/sim), the TAFM runs results/tafm (re/im already is the spin order parameter).
Crossings of successive sizes are interpolated on the common T grid (see
crossing()); histogram reweighting is not used because at L>=48 the T spacing
exceeds the energy-distribution width.
Also checks ordered vs random starts (L=32, files *_random.npz) run by run.

Usage (repo root): [DISCARD=20000] python3 results/smallD2/tc_small_d2.py
Writes results/smallD2/summary.txt, results/smallD2/u4_<model>_d2<D2>.txt,
results/smallD2/tc_vs_D2.png
"""
import glob
import os
import re

import numpy as np

HC = 0.266                     # TAFM ground-state KT field (Bloete-Nightingale)
D2S = [0.5, 0.3, 0.2, 0.1]
DISCARD = os.environ.get("DISCARD", "20000")
NB = 8
OUT = "results/smallD2"
KNOWN = {("tafm", 1.0): (0.4713, 0.0010), ("ann", 1.0): (0.528, 0.002), ("tafm", 3.0): (0.6718, 0.0007)}


def u4_points(model, D2):
    """{L: {T: (U4, err)}}: all seeds of one (L, T) pooled; error from a jackknife over
    8 time blocks per run (blocks of all seeds dropped together)."""
    pat = (f"results/long/long_L*_d2{D2:g}_T*_s*_ising.npz" if model == "ann"
           else f"results/tafm/long_L*_d2{D2:g}_T*_s*.npz")
    runs = {}
    for f in glob.glob(pat):
        if "_random" in f:
            continue
        d = np.load(f)
        k = int(DISCARD) // int(d["thin"])
        re_, im_ = (d["sre"], d["sim"]) if "sre" in d.files else (d["re"], d["im"])
        p2 = re_[k:, 0].astype(float) ** 2 + im_[k:, 0].astype(float) ** 2
        n = len(p2) // NB * NB
        runs.setdefault((int(d["L"]), round(float(d["T"]), 4)), []).append(p2[:n].reshape(NB, -1))
    U = lambda x: 1 - np.mean(x ** 2) / (2 * np.mean(x) ** 2)
    out = {}
    for (L, T), bl in runs.items():
        full = U(np.concatenate([b.ravel() for b in bl]))
        jk = np.array([U(np.concatenate([np.delete(b, i, 0).ravel() for b in bl])) for i in range(NB)])
        out.setdefault(L, {})[T] = (full, float(np.sqrt((NB - 1) / NB * np.sum((jk - jk.mean()) ** 2))), len(bl))
    return out


def crossing(a, b, nmc=2000, seed=0):
    """Crossing of U4(L_b) and U4(L_a), L_b > L_a, on their common T grid: the highest T
    at which d = U4(L_b) - U4(L_a) changes from >= 0 (larger L more ordered) to < 0,
    linearly interpolated. Median and 68% range over Gaussian resamplings of the points."""
    Ts = sorted(set(a) & set(b))
    if len(Ts) < 3:
        return np.nan, np.nan, np.nan
    d0 = np.array([b[T][0] - a[T][0] for T in Ts])
    e = np.array([np.hypot(b[T][1], a[T][1]) for T in Ts])
    rng = np.random.default_rng(seed)
    Tarr = np.array(Ts)
    res = []
    for _ in range(nmc):
        d = d0 + e * rng.standard_normal(len(Ts))
        idx = [i for i in range(len(Ts) - 1) if d[i] >= 0 > d[i + 1]]
        if not idx:
            res.append(Tarr[0] if d[0] < 0 else np.nan)
            continue
        i = idx[-1]
        res.append(Tarr[i] + d[i] * (Tarr[i + 1] - Tarr[i]) / (d[i] - d[i + 1]))
    res = np.array(res, float)
    ok = res[np.isfinite(res)]
    if len(ok) < nmc / 2:
        return np.nan, np.nan, np.nan
    lo, med, hi = np.percentile(ok, [16, 50, 84])
    return med, (hi - lo) / 2, np.mean(ok == Tarr[0])


def crossings(model, D2):
    pts = u4_points(model, D2)
    Ls = sorted(L for L in pts if len(pts[L]) >= 3)
    pairs = []
    for La, Lb in zip(Ls[:-1], Ls[1:]):
        t, e, edge = crossing(pts[La], pts[Lb])
        pairs.append((La, Lb, t, e, edge))
    with open(f"{OUT}/u4_{model}_d2{D2:g}.txt", "w") as fh:
        fh.write(f"# U4(psi_s) per (L, T), {model}, D2={D2:g}, discard {DISCARD}; value(err)[#seeds]\n")
        for T in sorted({T for L in Ls for T in pts[L]}):
            fh.write(f"T={T:<7g} " + "  ".join(f"L{L}:{pts[L][T][0]:+.4f}({pts[L][T][1]:.4f})[{pts[L][T][2]}]"
                                                 if T in pts[L] else f"L{L}:{'--':>17s}" for L in Ls) + "\n")
    return Ls, pairs


def best(pairs):
    """T_c estimate: crossing of the largest pair with a finite value."""
    ok = [p for p in pairs if np.isfinite(p[2])]
    return (ok[-1][2], ok[-1][3], f"{ok[-1][0]}/{ok[-1][1]}") if ok else (np.nan, np.nan, "-")


def start_check():
    """Ordered vs random start, run by run (same L, T, D2, model)."""
    rows = []
    for fr in sorted(glob.glob("results/long/long_L*_d2*_T*_s1_random_ising.npz")
                     + glob.glob("results/tafm/long_L*_d2*_T*_s1_random.npz")):
        fo = fr.replace("_random", "")
        if not os.path.exists(fo) or "_d21_" in fr:
            continue
        out = []
        for f in (fo, fr):
            d = np.load(f)
            k = int(DISCARD) // int(d["thin"])
            re_, im_ = (d["sre"], d["sim"]) if "sre" in d.files else (d["re"], d["im"])
            p2 = re_[k:, 0].astype(float) ** 2 + im_[k:, 0].astype(float) ** 2
            E = d["E"][k:, 0].astype(float) / int(d["N_tri"])
            nb = 8
            n = len(p2) // nb * nb
            U = lambda x: 1 - np.mean(x ** 2) / (2 * np.mean(x) ** 2)
            bl = lambda f_, x: np.array([f_(np.delete(x[:n], np.s_[b * n // nb:(b + 1) * n // nb])) for b in range(nb)])
            jk = lambda v: np.sqrt((nb - 1) / nb * np.sum((v - v.mean()) ** 2))
            out.append([(U(p2), jk(bl(U, p2))), (p2.mean(), jk(bl(np.mean, p2))), (E.mean(), jk(bl(np.mean, E)))])
        name = os.path.basename(fo)
        model = "ann" if "_ising" in name else "tafm"
        L = int(re.search(r"_L(\d+)_", name).group(1))
        D2 = float(re.search(r"_d2([\d.]+)_", name).group(1))
        T = float(re.search(r"_T([\d.]+)_", name).group(1))
        rows.append((model, D2, L, T, out))
    return sorted(rows)


def main():
    lines = []
    w = lines.append
    w("# Part IX 2b: T_c(D2) at small D2 -- annealed (pure Ising sweep, mu=0) vs fixed-bond TAFM")
    w(f"# Binder crossings of the spin order parameter psi_s (direct per-T U4, interpolated); jackknife over 8 time blocks; discard {DISCARD} sweeps")
    w(f"# H = D2/(2 T_c) (reduced field; TAFM ground-state KT value H_c = {HC}); asymptote T_c -> D2/(2 H_c) = {1 / (2 * HC):.3f} D2")
    w("")
    res = {}
    for D2 in D2S:
        for model in ("tafm", "ann"):
            Ls, pairs = crossings(model, D2)
            res[(model, D2)] = (Ls, pairs)
            w(f"D2={D2:g} {model:4s} sizes {Ls}: " + ", ".join(
                f"L={a}/{b}: {t:.3f}({e:.3f})" + (" [at grid edge]" if edge > 0.3 else "") for a, b, t, e, edge in pairs))
    w("")
    hdr = f"{'D2':>4} | {'T_c TAFM':>15} {'pair':>6} | {'T_c ann':>15} {'pair':>6} | {'ann/TAFM':>14} | {'H TAFM':>6} {'H ann':>6} | {'D2/(2Hc)':>8}"
    w(hdr)
    w("-" * len(hdr))
    table = []
    for D2 in D2S:
        (tt, et, pt), (ta, ea, pa) = best(res[("tafm", D2)][1]), best(res[("ann", D2)][1])
        r = ta / tt
        er = r * np.hypot(ea / ta, et / tt) if np.isfinite(r) else np.nan
        table.append((D2, tt, et, ta, ea))
        w(f"{D2:4g} | {tt:.4f}({et:.4f}) {pt:>6} | {ta:.4f}({ea:.4f}) {pa:>6} | {r:.3f}({er:.3f}) | "
          f"{D2 / (2 * tt):6.3f} {D2 / (2 * ta):6.3f} | {D2 / (2 * HC):8.3f}")
    for (m, D2), (t, e) in sorted(KNOWN.items(), key=lambda x: x[0][1]):
        w(f"# reference: D2={D2:g} {m}: T_c = {t}({e})  H = {D2 / (2 * t):.3f}")
    w("")
    w("# ordered vs random start (same L, T): U4, psi_s^2, E/N_tri  [z = (ordered - random)/combined error]")
    for model, D2, L, T, (o, r_) in start_check():
        zs = [(a[0] - b[0]) / np.hypot(a[1], b[1]) if np.hypot(a[1], b[1]) > 0 else 0 for a, b in zip(o, r_)]
        w(f"  {model:4s} D2={D2:<4g} L={L} T={T:<6g} U4 {o[0][0]:.3f}/{r_[0][0]:.3f} (z={zs[0]:+.1f})  "
          f"psi2 {o[1][0]:.4f}/{r_[1][0]:.4f} (z={zs[1]:+.1f})  E {o[2][0]:.5f}/{r_[2][0]:.5f} (z={zs[2]:+.1f})")
    text = "\n".join(lines) + "\n"
    open(f"{OUT}/summary.txt", "w").write(text)
    print(text)
    plot(table)


def plot(table):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array(table, float)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    a = ax[0]
    a.errorbar(t[:, 0], t[:, 1], t[:, 2], fmt="o", color="#3b4ba8", label="fixed-bond TAFM")
    a.errorbar(t[:, 0], t[:, 3], t[:, 4], fmt="s", color="#c5791f", label="annealed (Ising sweep)")
    a.errorbar([1.0, 3.0], [0.4713, 0.6718], [0.001, 0.0007], fmt="o", mfc="none", color="#3b4ba8", label="TAFM, earlier runs")
    a.errorbar([1.0], [0.528], [0.002], fmt="s", mfc="none", color="#c5791f", label="annealed D3=10, earlier runs")
    x = np.linspace(0, 1.05, 50)
    a.plot(x, x / (2 * HC), "k--", lw=1, label=r"$D_2/(2H_c)$, $H_c=0.266$")
    a.set_xlim(0, 1.05)
    a.set_ylim(0, 0.8)
    a.set_xlabel(r"$D_2$")
    a.set_ylabel(r"$T_c$")
    a.legend(fontsize=7)
    a = ax[1]
    a.errorbar(t[:, 0], t[:, 3] / t[:, 1], t[:, 3] / t[:, 1] * np.hypot(t[:, 4] / t[:, 3], t[:, 2] / t[:, 1]),
               fmt="o", color="k", label="this work")
    a.plot([1.0], [0.528 / 0.4713], "o", mfc="none", color="k", label="D2=1, earlier runs")
    a.axhline(1, color="gray", lw=0.8)
    a.set_xlabel(r"$D_2$")
    a.set_ylabel(r"$T_c^{\rm ann}/T_c^{\rm TAFM}$")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{OUT}/tc_vs_D2.png", dpi=150)


if __name__ == "__main__":
    main()

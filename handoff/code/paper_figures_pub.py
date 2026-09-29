"""Publication figures for the Part IX paper, drawn from the raw time series.

Writes figures/pub/fig{2..6}_*.png/.pdf and, for every plotted curve, the numbers
behind it as handoff/data/figdata/*.csv, so the figures can be redrawn without the
raw .npz files.

  fig2_binder   Binder U4(T): annealed (D2=1, mu=0; hub order parameter) and the
                fixed-bond TAFM at the same field h/J=1 (minority-spin order parameter)
  fig3_xi_def   D1=D2 (h=0): correlation length and spin-defect density, annealed vs TAFM
  fig4_heat     heat capacity near T_c at L=64: annealed C_mu and C_n, TAFM
  fig5_tc_mu    T_c versus mu
  fig6_eta      effective eta (psi_s^2 ~ L^-eta, L=32-96) versus T/T_c for mu=1..10,
                and <cos 3theta> versus L at mu=10

Usage: python3 paper_figures_pub.py
"""
import csv
import glob
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import annealed_fss_analysis as fa
import annealed_long_analysis as la

OUT = "figures/pub"
DAT = "handoff/data/figdata"
# validated palette (dataviz reference instance): categorical slots 1-4 for system sizes,
# a 5-step blue ordinal ramp for mu
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
ORD = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
MARK = ["o", "s", "^", "D", "v"]
LCOL = {32: CAT[0], 48: CAT[1], 64: CAT[2], 24: CAT[3]}   # colour follows the system size in every panel

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "legend.frameon": False, "savefig.dpi": 200, "figure.facecolor": "white",
})


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_csv(name, header, rows):
    with open(f"{DAT}/{name}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def reweighted(dirname, Ls, order, mu=0.0, npts=300, Tmin=0.0):
    """U4, C, C_n on a T grid per L from the single-T long runs (exact-level WHAM)."""
    la.DIR, la.ORDER, la.MU = dirname, order, mu
    la._cache.clear()
    out = {}
    for L in Ls:
        runs = [r for r in la.load_blocks(1.0, L) if 1 / r["betas"][0] >= Tmin]   # drop coarse low-T points: no histogram overlap
        b = np.array([r["betas"][0] for r in runs])
        g = np.linspace(b.min(), b.max(), npts)
        cv = fa.curves(runs, g)
        out[L] = (1 / g, cv, sorted(set(1 / b)))
    return out


def fig2():
    ann = reweighted("results/long", [32, 48, 64], "hub")
    taf = reweighted("results/tafm", [24, 32, 48, 64], "hub")   # tafm files store the spin psi in re/im
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    rows = []
    for a, res, tc, title, tag in ((ax[0], ann, 0.528, "annealed bonds (D2=1, μ=0)", "annealed"),
                                   (ax[1], taf, 0.4713, "fixed-bond TAFM (h/J=1)", "tafm")):
        for L, (T, cv, Ts) in sorted(res.items()):
            c = LCOL[L]
            a.plot(T, cv["U"], color=c, label=f"L={L}")
            rows += [(tag, L, f"{t:.5f}", f"{u:.5f}") for t, u in zip(T, cv["U"])]
        a.axvline(tc, color=INK2, lw=1, ls=":")
        a.text(tc, 0.26, f"T_c={tc} ", color=INK2, fontsize=8, ha="right")
        a.set_title(title, fontsize=10, color=INK)
        a.set_xlabel("T")
        a.legend(fontsize=8, loc="lower left")
    ax[0].set_ylabel("Binder cumulant U4")
    ax[0].set_ylim(0, 0.52)
    save(fig, "fig2_binder")
    write_csv("fig2_binder", ["model", "L", "T", "U4"], rows)


def xi_pt(f):
    d = np.load(f)
    L = int(d["L"])
    ph = d["phi"]
    S = np.abs(ph) ** 2 * int(d["N_site"])
    s0, s1 = S[:, :, 0].mean(0), S[:, :, 1:].mean((0, 2))
    return 1 / d["betas"], np.sqrt(np.clip(s0 / s1 - 1, 0, None)) / (2 * np.sin(np.pi / L))


def fig3():
    xs = [xi_pt(f) for f in sorted(glob.glob("results/d1d2/pt_L32_s*.npz"))]
    T_a = xs[0][0]
    xa = np.array([x[1] for x in xs])
    xa_m, xa_e = xa.mean(0), xa.std(0, ddof=1) / np.sqrt(len(xa))
    t = np.loadtxt("results/tafm/xi_L32_s1.txt")
    d = np.loadtxt("results/tafm/defects_L24_s1.txt")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    a = ax[0]
    a.errorbar(1 / T_a, xa_m, yerr=xa_e, fmt="o-", ms=4, color=CAT[0], label="annealed (4 seeds)")
    a.plot(1 / t[:, 0], t[:, 3], "s-", ms=4, color=CAT[1], label="fixed-bond TAFM (1 seed)")
    a.axhline(16, color=INK2, lw=1, ls=":")
    a.text(1.0, 16, " L/2", color=INK2, fontsize=8, va="bottom")
    a.set_yscale("log")
    a.set_xlabel("1/T")
    a.set_ylabel("ξ (L=32, zero field)")
    a.legend(fontsize=8)
    a = ax[1]
    a.plot(1 / d[:, 0], d[:, 1], "o-", ms=4, color=CAT[0], label="spin defects, annealed")
    a.plot(1 / d[:, 0], d[:, 2], "s-", ms=4, color=CAT[1], label="spin defects, TAFM")
    a.plot(1 / d[:, 0], d[:, 3], "^--", ms=4, lw=1.5, color=CAT[2], label="frustrated triangles, annealed")
    a.set_yscale("log")
    a.set_xlabel("1/T")
    a.set_ylabel("density per triangle (L=24)")
    a.legend(fontsize=8)
    save(fig, "fig3_xi_defects")
    write_csv("fig3a_xi", ["model", "T", "xi", "err"],
              [("annealed", f"{a_:.4f}", f"{b:.4f}", f"{c:.4f}") for a_, b, c in zip(T_a, xa_m, xa_e)] +
              [("tafm", f"{r[0]:.4f}", f"{r[3]:.4f}", "") for r in t])
    write_csv("fig3b_defects", ["T", "n_def_annealed", "n_def_tafm", "n_mon_annealed"],
              [[f"{x:.4g}" for x in r] for r in d])


def fig4():
    ann = reweighted("results/long", [64], "hub", Tmin=0.51)
    taf = reweighted("results/tafm", [64], "hub")
    T, cv, _ = ann[64]
    Tt, ct, _ = taf[64]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.plot(T / 0.528, cv["C"], color=CAT[0], label="annealed C_μ")
    ax.plot(T / 0.528, cv["Cn"], color=CAT[0], ls="--", label="annealed C_n (fixed n)")
    ax.plot(Tt / 0.4713, ct["C"], color=CAT[1], label="fixed-bond TAFM")
    ax.set_xlim(0.97, 1.035)
    ax.set_xlabel("T / T_c")
    ax.set_ylabel("C / N_tri  (L=64)")
    ax.legend(fontsize=8)
    save(fig, "fig4_heat")
    write_csv("fig4_heat", ["model", "T", "T_over_Tc", "C", "C_n"],
              [("annealed", f"{a:.5f}", f"{a / 0.528:.5f}", f"{b:.5f}", f"{c:.5f}") for a, b, c in zip(T, cv["C"], cv["Cn"])] +
              [("tafm", f"{a:.5f}", f"{a / 0.4713:.5f}", f"{b:.5f}", "") for a, b in zip(Tt, ct["C"])])


TC = {0: (0.528, 0.002), 0.5: (0.620, 0.007), 1: (0.696, 0.004), 2: (0.80, 0.03), 3: (0.95, 0.03),
      5: (1.05, 0.03), 10: (1.141, 0.003)}


def fig5():
    mu = np.array(sorted(TC))
    tc = np.array([TC[m][0] for m in mu])
    e = np.array([TC[m][1] for m in mu])
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.errorbar(mu, tc, yerr=e, fmt="o-", ms=6, color=CAT[0], capsize=3, label="annealed, D2=1")
    ax.axhline(0.4713, color=CAT[1], lw=1.5, ls="--")
    ax.text(10, 0.4713, "fixed-bond TAFM, h/J=1", color=INK2, fontsize=8, ha="right", va="bottom")
    ax.set_xlabel("μ (frustrated-triangle cost)")
    ax.set_ylabel("T_c")
    ax.set_ylim(0.4, 1.25)
    save(fig, "fig5_tc_mu")
    write_csv("fig5_tc_mu", ["mu", "T_c", "err", "note"],
              [(m, TC[m][0], TC[m][1], "approximate" if TC[m][1] >= 0.02 else "") for m in mu])


def spin_series():
    """psi_s^2 per (mu, T, L) pooled over seeds with tau-based errors."""
    d = defaultdict(lambda: defaultdict(list))
    for f in glob.glob("results/long/long_L*_d21*_T*_s*.npz"):
        m = re.search(r"L(\d+)_d21(?:_mu([\d.]+))?_T([\d.]+)_s(\d+)\.npz$", f)
        if not m:
            continue
        L, mu, T = int(m[1]), float(m[2] or 0), float(m[3])
        z = np.load(f)
        if "sre" not in z.files or L < 32:
            continue
        k = 100000 // int(z["thin"])
        p = z["sre"][k:, 0].astype(float) ** 2 + z["sim"][k:, 0].astype(float) ** 2
        th = np.angle(z["sre"][k:, 0] + 1j * z["sim"][k:, 0])
        d[(mu, T)][L].append((p.mean(), np.sqrt(p.var() * 2 * la.tau_int(p) / len(p)), np.cos(3 * th).mean()))
    return d


def eta_fit(dd):
    Ls = np.array(sorted(dd))
    ps, ss = [], []
    for L in Ls:
        a = np.array(dd[L])[:, :2]
        w = 1 / a[:, 1] ** 2
        p = (a[:, 0] * w).sum() / w.sum()
        e = np.sqrt(1 / w.sum())
        if len(a) > 1:
            e *= np.sqrt(max(1, ((a[:, 0] - p) ** 2 * w).sum() / (len(a) - 1)))
        ps.append(p)
        ss.append(e / p)
    y, s = np.log(ps), np.array(ss)
    X = np.vstack([np.ones(len(Ls)), np.log(Ls)]).T
    W = np.diag(1 / s ** 2)
    C = np.linalg.inv(X.T @ W @ X)
    b = C @ X.T @ W @ y
    chi = ((y - X @ b) ** 2 / s ** 2).sum() / max(1, len(Ls) - 2)
    return -b[1], np.sqrt(C[1, 1] * max(1, chi)), Ls


def fig6():
    d = spin_series()
    mus = [1, 2, 3, 5, 10]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    rows = []
    a = ax[0]
    for mu, c, mk in zip(mus, ORD, MARK):
        pts = []
        for (m, T), dd in d.items():
            if m != mu or len(dd) < 3:
                continue
            eta, err, Ls = eta_fit(dd)
            if eta > 1.0:
                continue
            pts.append((T / TC[mu][0], eta, err, T, Ls.max()))
        pts.sort()
        if not pts:
            continue
        x, y, e = (np.array([p[i] for p in pts]) for i in range(3))
        a.errorbar(x, y, yerr=e, fmt=mk + "-", ms=5, lw=1.5, color=c, capsize=2, label=f"μ={mu}")
        rows += [(mu, f"{p[3]:.3f}", f"{p[0]:.4f}", f"{p[1]:.4f}", f"{p[2]:.4f}", p[4]) for p in pts]
    for v, lab in ((4 / 15, "4/15 (3-state Potts)"), (4 / 9, "4/9 (Z3 locking)")):
        a.axhline(v, color=INK2, lw=1, ls=":")
        a.text(0.605, v, lab, color=INK2, fontsize=8, va="bottom")
    a.set_xlabel("T / T_c")
    a.set_ylabel("effective η  (ψ_s² ∝ L^−η, L=32–96)")
    a.set_xlim(0.6, 1.02)
    a.set_ylim(-0.1, 0.9)
    a.legend(fontsize=8, loc="upper left")
    write_csv("fig6a_eta", ["mu", "T", "T_over_Tc", "eta_eff", "err", "L_max"], rows)
    a = ax[1]
    rows = []
    for T, c, mk in zip((0.7, 1.04, 1.08, 1.14), CAT, MARK):
        dd = d.get((10.0, T))
        if not dd:
            continue
        Ls = sorted(dd)
        y = [np.mean([r[2] for r in dd[L]]) for L in Ls]
        a.plot(Ls, y, mk + "-", ms=6, color=c, label=f"T/T_c={T / TC[10][0]:.2f}")
        rows += [(T, L, f"{v:.4f}") for L, v in zip(Ls, y)]
    a.set_xscale("log")
    a.set_xticks([32, 48, 64, 96])
    a.set_xticklabels(["32", "48", "64", "96"])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xlabel("L")
    a.set_ylabel("⟨cos 3θ⟩  (μ=10)")
    a.set_ylim(0, 1.05)
    a.legend(fontsize=8)
    write_csv("fig6b_cos3", ["T", "L", "cos3theta"], rows)
    save(fig, "fig6_eta")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(DAT, exist_ok=True)
    for f in (fig2, fig3, fig4, fig5, fig6):
        f()
        print("done", f.__name__, flush=True)

"""Publication figures (vector PDF) for main.tex.   Run:  cd paper && python3 make_figs.py

Writes figs/fig_*.pdf and rewrites the generated timing table in main.tex (between the
'% BEGIN generated timing table' / '% END generated timing table' markers) from
../hub_runs/timing_cert.jsonl.  Re-runnable: the timing benchmark appends lines, and a rerun
picks up whatever lines exist.

Data sources (all in d3_complexity/):
  report/penrose_states.json   Penrose face-adjacency ground states (radius-6 patch, 140 rhombi)
  penrose_runs/r*_s*.json      exact energy curves of the ten Penrose patches (segments)
  hub_runs/timing_cert.jsonl   certification timing (fam, n, fr, healed, t_match, t_2sat)
  lattices3.py, heal_check.py, ising_tjoin.py, potts_exact.py   lattices and solvers

Style: state 1 = blue circle, state 2 = orange square, state 3 = yellow triangle / hatched
tile; frustrated (same-state) bonds thick dark red with a white casing.  Small graphs carry
the state digit inside each site, so every figure reads in greyscale.
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np
import matplotlib
import matplotlib.ticker

matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection, LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "report"))
OUT = os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.1,
    "hatch.linewidth": 0.45,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})
COL1 = 1.0 * 3.4
COL2 = 7.0

# Okabe-Ito based; grey levels differ (dark / medium / light) so states survive greyscale.
SC = {1: "#0072B2", 2: "#E69F00", 3: "#F0E442"}
TXT = {1: "white", 2: "black", 3: "black"}
MARK = {1: "o", 2: "s", 3: "^"}
FR = "#A50F15"          # frustrated bond
BOND = "0.45"
ODD = "0.86"            # odd faces
INK = "0.15"


def panel_label(ax, s, x=0.0, y=1.0, **kw):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, va="top", ha="left", **kw)


def blank(ax):
    ax.set_aspect("equal")
    ax.axis("off")


def draw_graph(ax, pos, edges, col, r=0.07, fs=6.5, ring=(), lw=0.8, digits=True):
    """Small graph: bonds (frustrated = thick red), sites as circles with state digits."""
    for u, v in edges:
        (x1, y1), (x2, y2) = pos[u], pos[v]
        cu, cv = col.get(u), col.get(v)
        if cu is not None and cu == cv:
            ax.plot([x1, x2], [y1, y2], color="white", lw=3.6, solid_capstyle="butt", zorder=1)
            ax.plot([x1, x2], [y1, y2], color=FR, lw=2.2, solid_capstyle="butt", zorder=1.1)
        else:
            ax.plot([x1, x2], [y1, y2], color=BOND, lw=lw, zorder=1)
    for v, (x, y) in pos.items():
        c = col.get(v)
        fc = SC[c] if c else "white"
        ax.add_patch(Circle((x, y), r, fc=fc, ec=INK, lw=0.6, zorder=3))
        if digits and c:
            ax.text(x, y - 0.004, str(c), ha="center", va="center", fontsize=fs, color=TXT[c],
                    zorder=4, fontweight="bold")
        if v in ring:
            ax.add_patch(Circle((x, y), r * 1.55, fc="none", ec=INK, lw=0.7, ls=(0, (2, 1.2)), zorder=3))


def draw_sites(ax, P, s, ms=3.0, big=1.35, lw=0.35):
    P = np.asarray(P, float)
    s = np.asarray(s)
    for c in (1, 2, 3):
        m = s == c
        if m.any():
            size = (ms * (big if c == 3 else 1.0)) ** 2
            ax.scatter(P[m, 0], P[m, 1], s=size, marker=MARK[c], c=SC[c], edgecolors=INK,
                       linewidths=lw, zorder=3)


def draw_bonds(ax, P, edges, s=None, lw=0.5, frlw=1.6):
    P = np.asarray(P, float)
    norm, fr = [], []
    for u, v in edges:
        (fr if (s is not None and s[u] == s[v]) else norm).append([P[u], P[v]])
    ax.add_collection(LineCollection(norm, colors=BOND, linewidths=lw, zorder=1))
    if fr:
        ax.add_collection(LineCollection(fr, colors="white", linewidths=frlw + 1.2, zorder=1.5,
                                         capstyle="butt"))
        ax.add_collection(LineCollection(fr, colors=FR, linewidths=frlw, zorder=1.6, capstyle="butt"))


def state_legend(fig_or_ax, loc, ncol=4, frus=True, odd=False, bbox=None, marker=True, **kw):
    h = []
    for c in (1, 2, 3):
        if marker:
            h.append(Line2D([], [], ls="", marker=MARK[c], mfc=SC[c], mec=INK, mew=0.5,
                            ms=5.5 if c == 3 else 4.5, label=f"state {c}"))
        else:
            h.append(Polygon([[0, 0]], fc=SC[c], ec=INK, lw=0.4, hatch="////" if c == 3 else None,
                             label=f"state {c}"))
    if frus:
        h.append(Line2D([], [], color=FR, lw=2.0, label="frustrated bond"))
    if odd:
        h.append(Polygon([[0, 0]], fc=ODD, ec="none", label="odd face"))
    return fig_or_ax.legend(handles=h, loc=loc, ncol=ncol, frameon=False, bbox_to_anchor=bbox,
                            handlelength=1.4, handletextpad=0.4, columnspacing=1.0, borderaxespad=0.0, **kw)


def star_pos(cx, cy, R, d, start=90):
    pos = {"c": (cx, cy)}
    for k in range(d):
        a = math.radians(start - 360 * k / d)
        pos[k] = (cx + R * math.cos(a), cy + R * math.sin(a))
    return pos


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
                                 lw=0.7, color=INK, zorder=2))


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


# ================================================================== 1. recolouring
def fig_recolour():
    W, H = COL1, 2.40
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    blank(ax)
    r = 0.075
    # (a) triangle
    ax.text(0.02, H - 0.02, "(a)", fontsize=9, va="top")
    for cx, cols, t1, t2 in ((0.95, (1, 2, 1), "two states", "$E=1$"),
                             (2.35, (1, 2, 3), "third state", "$E=D_3$")):
        s = 0.27
        pos = {0: (cx, 2.28), 1: (cx - s, 2.28 - 1.732 * s), 2: (cx + s, 2.28 - 1.732 * s)}
        draw_graph(ax, pos, [(0, 1), (1, 2), (0, 2)], dict(enumerate(cols)), r=r)
        ax.text(cx + 0.42, 2.03, t1 + "\n" + t2, ha="left", va="center", fontsize=8, linespacing=1.3)
    # (b) local recolouring, d = 4
    ax.plot([0.05, W - 0.05], [1.66, 1.66], color="0.75", lw=0.5)
    ax.text(0.02, 1.62, "(b)", fontsize=9, va="top")
    R = 0.27
    for y, nb, lab in ((1.22, [1, 2, 1, 2], "$n_1{:}n_2=2{:}2$\n$\\Delta E=2-D_3$"),
                       (0.40, [1, 2, 2, 2], "$n_1{:}n_2=1{:}3$\n$\\Delta E=1-D_3$")):
        for cx, cc in ((0.62, 3), (1.72, 1)):
            pos = star_pos(cx, y, R, 4, start=45)
            col = {"c": cc}; col.update(dict(enumerate(nb)))
            draw_graph(ax, pos, [("c", k) for k in range(4)], col, r=r)
        arrow(ax, 0.98, y, 1.36, y)
        ax.text(1.17, y + 0.06, "$3\\to1$", ha="center", va="bottom", fontsize=7)
        ax.text(2.20, y, lab, ha="left", va="center", fontsize=8, linespacing=1.35)
    save(fig, "fig_recolour.pdf")


# ================================================================== 2. wheels
def wheel_pos(cx, cy, R, D):
    pos = {"h": (cx, cy)}
    for k in range(D):
        a = math.radians(90 - 360 * k / D)
        pos[k] = (cx + R * math.cos(a), cy + R * math.sin(a))
    E = [("h", k) for k in range(D)] + [(k, (k + 1) % D) for k in range(D)]
    return pos, E


def fig_wheel():
    fig = plt.figure(figsize=(COL1, 3.05))
    axg = fig.add_axes([0, 0.60, 1, 0.40])
    Wd, Hd = COL1, 0.40 * 3.05
    axg.set_xlim(0, Wd); axg.set_ylim(0, Hd); blank(axg)
    R = 0.30
    specs = [(0.43, 4, 1, [2, 1, 2, 1], "$E=2$", "$D_3>2$"),
             (1.23, 4, 3, [2, 1, 2, 1], "$E=D_3$", "$D_3<2$"),
             (2.17, 5, 1, [2, 1, 2, 1, 2], "$E=3$", "$D_3>2$"),
             (2.97, 5, 3, [2, 1, 2, 1, 2], "$E=D_3+1$", "$D_3<2$")]
    cy = 0.66
    for cx, D, hub, rim, e, reg in specs:
        pos, E = wheel_pos(cx, cy, R, D)
        col = {"h": hub}; col.update(dict(enumerate(rim)))
        draw_graph(axg, pos, E, col, r=0.068, fs=6)
        axg.text(cx, cy - R - 0.14, e, ha="center", va="center", fontsize=8)
        axg.text(cx, cy - R - 0.30, "g.s. for " + reg, ha="center", va="center", fontsize=7, color="0.3")
    axg.text(0.83, Hd - 0.02, "$W_4$", ha="center", va="top", fontsize=8.5)
    axg.text(2.57, Hd - 0.02, "$W_5$", ha="center", va="top", fontsize=8.5)
    axg.plot([1.70, 1.70], [0.1, Hd - 0.05], color="0.75", lw=0.5)
    panel_label(axg, "(a)", x=0.005, y=1.0)

    ax = fig.add_axes([0.15, 0.10, 0.82, 0.44])
    x = np.linspace(0, 3.5, 200)
    for (two, off), ls, name in (((2, 0), "-", "W_4"), ((3, 1), "--", "W_5")):
        ax.plot(x, np.minimum(two, x + off), color="k", ls=ls, lw=1.3, label=f"${name}$", zorder=3)
        ax.plot([0, 2], [two, two], color="0.6", ls=":", lw=0.8)
        ax.plot([2, 3.5], [2 + off, 3.5 + off], color="0.6", ls=":", lw=0.8)
    ax.axvline(2, color=FR, lw=0.7, ls=(0, (4, 2)))
    ax.text(2.06, 0.25, "$\\lfloor\\Delta/2\\rfloor=2$", color=FR, fontsize=7.5, ha="left")
    ax.text(0.62, 0.30, "hub in state 3", fontsize=7, color="0.25")
    ax.text(2.55, 3.12, "two states", fontsize=7, color="0.25")
    ax.set_xlim(0, 3.5); ax.set_ylim(0, 4)
    ax.set_xlabel("$D_3$"); ax.set_ylabel("$E_0$")
    ax.set_xticks([0, 1, 2, 3]); ax.set_yticks([0, 1, 2, 3, 4])
    ax.legend(loc="lower right", frameon=False, handlelength=2.2)
    panel_label(ax, "(b)", x=-0.17, y=1.08)
    save(fig, "fig_wheel.pdf")


# ================================================================== 3. triangular lattice
def tri_patch(colfun, nx=8, ny=6):
    pos, col, E = {}, {}, []
    for j in range(ny):
        for i in range(nx):
            pos[(i, j)] = (i + 0.5 * (j % 2), -0.8660254 * j)
            col[(i, j)] = colfun(i, j)
    for j in range(ny):
        for i in range(nx):
            nb = [(i + 1, j)] + ([(i, j + 1), (i - 1, j + 1)] if j % 2 == 0 else [(i, j + 1), (i + 1, j + 1)])
            E += [((i, j), q) for q in nb if q in pos]
    keys = list(pos)
    idx = {k: t for t, k in enumerate(keys)}
    return np.array([pos[k] for k in keys]), [(idx[a], idx[b]) for a, b in E], [col[k] for k in keys]


def fig_triangular():
    def three(i, j):
        q = i - (j - (j & 1)) // 2
        return 1 + ((q - j) % 3)

    def stripes(i, j):
        return 1 if j % 2 == 0 else 2

    fig = plt.figure(figsize=(COL1, 3.0))
    for k, (f, lab, e) in enumerate(((three, "(a)", "$E/N=D_3/3$"), (stripes, "(b)", "$E/N=1$"))):
        ax = fig.add_axes([0.02 + 0.5 * k, 0.53, 0.46, 0.40])
        P, E, s = tri_patch(f)
        draw_bonds(ax, P, E, s, lw=0.55, frlw=1.5)
        draw_sites(ax, P, s, ms=4.2)
        blank(ax)
        ax.set_xlim(P[:, 0].min() - 0.35, P[:, 0].max() + 0.35)
        ax.set_ylim(P[:, 1].min() - 0.35, P[:, 1].max() + 0.35)
        panel_label(ax, lab, x=-0.03, y=1.10)
        ax.text(0.5, -0.04, e, transform=ax.transAxes, ha="center", va="top", fontsize=8)
    state_legend(fig, "upper center", ncol=4, bbox=(0.5, 1.0))
    ax = fig.add_axes([0.15, 0.09, 0.82, 0.31])
    x = np.linspace(0, 4.5, 300)
    ax.plot(x, np.minimum(x / 3, 1), color="k", lw=1.4, zorder=3, label="$E_0/N$")
    ax.plot(x, x / 3, color=SC[1], ls="--", lw=0.9, label="three-sublattice, $D_3/3$")
    ax.plot(x, np.ones_like(x), color=SC[2], ls="-.", lw=0.9, label="Wannier, $1$")
    ax.axvline(3, color=FR, lw=0.7, ls=(0, (4, 2)))
    ax.text(3.07, 0.12, "$D_3^*=3$", color=FR, fontsize=7.5)
    ax.set_xlim(0, 4.5); ax.set_ylim(0, 1.4)
    ax.set_xlabel("$D_3$"); ax.set_ylabel("$E/N$")
    ax.set_yticks([0, 0.5, 1]); ax.set_xticks([0, 1, 2, 3, 4])
    ax.legend(loc="upper left", frameon=False, handlelength=2.0, borderaxespad=0.2)
    panel_label(ax, "(c)", x=-0.17, y=1.12)
    save(fig, "fig_triangular.pdf")


# ================================================================== 4. Penrose face-adjacency
def load_penrose_curves():
    curves = {}
    for p in sorted(glob.glob(os.path.join(ROOT, "penrose_runs", "r*_s*.json"))):
        try:
            txt = open(p).read()
            d = json.loads(txt[txt.index("{"):])      # some files start with HiGHS log lines
        except (ValueError, OSError):
            continue                    # run not finished
        key = os.path.basename(p)[:-5]
        curves[key] = (d["n"], [(g["mono"], g["n3"]) for g in d["segments"]])
    return curves


def envelope(lines, x):
    return np.min([m + k * x for m, k in lines], axis=0)


def fig_penrose():
    d = json.load(open(os.path.join(ROOT, "report", "penrose_states.json")))
    polys = [np.array(p) for p in d["polys"]]
    fig = plt.figure(figsize=(COL2, 2.25))
    w = 0.215
    for k, key in enumerate(("0.5", "1.5", "2.5")):
        ax = fig.add_axes([0.005 + k * (w + 0.01), 0.02, w, 0.76])
        st = d["states"][key]
        col = st["col"]
        for c in (1, 2, 3):
            pc = [p for p, cc in zip(polys, col) if cc == c]
            if pc:
                ax.add_collection(PolyCollection(pc, facecolors=SC[c], edgecolors="0.3", linewidths=0.3,
                                                 hatch="//////" if c == 3 else None, zorder=1))
        seg = [e[2] for e in d["edges"] if col[e[0]] == col[e[1]]]
        if seg:
            ax.add_collection(LineCollection(seg, colors="white", linewidths=2.6, capstyle="round", zorder=2))
            ax.add_collection(LineCollection(seg, colors=FR, linewidths=1.5, capstyle="round", zorder=2.1))
        allp = np.vstack(polys)
        ax.set_xlim(allp[:, 0].min() - 0.1, allp[:, 0].max() + 0.1)
        ax.set_ylim(allp[:, 1].min() - 0.1, allp[:, 1].max() + 0.1)
        ax.set_aspect("equal", adjustable="datalim"); ax.axis("off")
        panel_label(ax, "(" + "abc"[k] + ")", x=0.0, y=1.12)
        ax.text(0.13, 1.12, f"$D_3={key}$:  $n_3={st['n3']}$, $m={st['mono']}$", transform=ax.transAxes,
                ha="left", va="top", fontsize=7.5)
    state_legend(fig, "upper left", ncol=4, marker=False, bbox=(0.005, 1.0))

    ax = fig.add_axes([0.735, 0.17, 0.255, 0.63])
    curves = load_penrose_curves()
    x = np.linspace(0, 3, 601)
    for key, (n, lines) in curves.items():
        if key != "r6_s0":
            ax.plot(x, envelope(lines, x) / n, color="0.72", lw=0.6, zorder=1)
    n, lines = curves["r6_s0"]
    ax.plot(x, envelope(lines, x) / n, color="k", lw=1.3, zorder=3)
    kinks = []
    for (m1, k1), (m2, k2) in zip(lines, lines[1:]):
        xk = (m2 - m1) / (k1 - k2)
        kinks.append((xk, (m1 + k1 * xk) / n))
    ax.plot(*zip(*kinks), ls="", marker="o", ms=3.2, mfc="white", mec="k", mew=0.8, zorder=4)
    for xs in (0.5, 1.5, 2.5):
        ax.axvline(xs, color="0.55", lw=0.5, ls=(0, (1, 1.5)), zorder=0)
    ax.text(0.52, 0.02, "(a)", fontsize=7, color="0.35", ha="left", va="bottom")
    ax.text(1.52, 0.02, "(b)", fontsize=7, color="0.35", ha="left", va="bottom")
    ax.text(2.52, 0.02, "(c)", fontsize=7, color="0.35", ha="left", va="bottom")
    ax.axvline(2, color=FR, lw=0.7, ls=(0, (4, 2)), zorder=0)
    ax.set_xlim(0, 3); ax.set_ylim(0, 0.6)
    ax.set_xticks([0, 1, 2, 3]); ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("$D_3$"); ax.set_ylabel("$E_0/N$", labelpad=2)
    panel_label(ax, "(d)", x=-0.30, y=1.13)
    save(fig, "fig_penrose.pdf")
    return curves


# ================================================================== 5. healing
def fig_healing():
    from ising_tjoin import ising_ground_state_tjoin
    from heal_check import two_sat_heal
    from lattices3 import truncated_penrose
    from potts_exact import solve

    # (a) triangular prism, an Ising ground state whose frustrated bonds are joined by spokes
    E = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)]
    col = [1, 1, 2, 2, 2, 1]
    fr, _ = ising_ground_state_tjoin(6, E)
    assert fr == sum(col[u] == col[v] for u, v in E) == 2
    chosen, M = two_sat_heal(6, E, col)
    healed = list(col)
    for v in chosen:
        healed[v] = 3
    assert all(healed[u] != healed[v] for u, v in E)

    fig = plt.figure(figsize=(COL1, 2.85))
    axg = fig.add_axes([0, 0.60, 1, 0.40])
    Wd, Hd = COL1, 0.40 * 2.85
    axg.set_xlim(0, Wd); axg.set_ylim(0, Hd); blank(axg)
    R1, R2, cy = 0.36, 0.14, 0.66
    titles = ("Ising g.s., $\\mathrm{fr}=2$", "one endpoint each", "healed, $E=2D_3$")
    for k, (cx, cc, ring) in enumerate(((0.52, col, ()), (1.70, col, set(chosen)), (2.88, healed, ()))):
        pos = {}
        for i in range(3):
            a = math.radians(90 + 120 * i)
            pos[i] = (cx + R1 * math.cos(a), cy + 0.03 + R1 * math.sin(a))
            pos[i + 3] = (cx + R2 * math.cos(a), cy + 0.03 + R2 * math.sin(a))
        draw_graph(axg, pos, E, dict(enumerate(cc)), r=0.065, fs=6, ring=ring)
        axg.text(cx, 0.17, titles[k], ha="center", va="center", fontsize=7.5)
        if k < 2:
            arrow(axg, cx + 0.40, cy - 0.05, cx + 0.78, cy - 0.05)
    panel_label(axg, "(a)", x=0.005, y=1.0)

    # (b) truncated Penrose patch: certified E0 = min(D3,1) fr, checked against exact MILP
    n, Et, pos, info = truncated_penrose(4, seed=0)
    frt, colt = ising_ground_state_tjoin(n, Et)
    ch, _ = two_sat_heal(n, Et, colt)
    assert ch is not None
    xs = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]
    ex = []
    for D in xs:
        mono, n3, _ = solve(n, Et, D)
        ex.append(mono + D * n3)
    ex = np.array(ex)
    assert np.allclose(ex, np.minimum(xs, 1) * frt)
    ax = fig.add_axes([0.15, 0.11, 0.82, 0.42])
    x = np.linspace(0, 2, 200)
    ax.plot(x, np.minimum(x, 1) * frt / n, color="k", lw=1.3, label="$\\min(D_3,1)\\,\\mathrm{fr}/N$", zorder=2)
    ax.plot([1, 2], [frt / n, 2 * frt / n], color="0.55", ls=":", lw=0.8)
    ax.plot([0, 1], [frt / n, frt / n], color="0.55", ls=":", lw=0.8)
    ax.plot(xs, ex / n, ls="", marker="o", ms=3.8, mfc="white", mec="k", mew=0.8, zorder=3,
            label="exact (MILP)")
    ax.set_xlim(0, 2); ax.set_ylim(0, 0.3)
    ax.set_xticks([0, 0.5, 1, 1.5, 2]); ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_xlabel("$D_3$"); ax.set_ylabel("$E_0/N$")
    ax.legend(loc="lower right", frameon=False, handlelength=1.8)
    ax.text(0.04, 0.265, f"truncated Penrose, $N={n}$, $\\mathrm{{fr}}={frt}$", fontsize=7.5)
    panel_label(ax, "(b)", x=-0.17, y=1.10)
    save(fig, "fig_healing.pdf")
    return n, frt


# ================================================================== 6. truncation + gallery
def fig_truncation():
    from heal_check import healed_state
    from lattices3 import truncated_penrose, random_cubic
    from deg3_gallery import geometric_faces, star_lattice

    fig = plt.figure(figsize=(COL2, 2.2))
    # (a) schematic
    ax = fig.add_axes([0.0, 0.10, 0.20, 0.80])
    ax.set_xlim(0, 1.40); ax.set_ylim(0, 1.76); blank(ax)
    for y0, d in ((1.30, 4), (0.45, 5)):
        R = 0.30
        for cx, trunc in ((0.33, False), (1.07, True)):
            ang = [math.radians(90 - 360 * k / d + (45 if d == 4 else 0)) for k in range(d)]
            out = [(cx + R * math.cos(a), y0 + R * math.sin(a)) for a in ang]
            if not trunc:
                for p in out:
                    ax.plot([cx, p[0]], [y0, p[1]], color=INK, lw=0.9)
                ax.add_patch(Circle((cx, y0), 0.04, fc="k", ec="none", zorder=3))
                ax.text(cx + 0.02, y0 - R - 0.06, f"$d={d}$", ha="center", va="top", fontsize=7.5)
            else:
                inn = [(cx + 0.11 * math.cos(a), y0 + 0.11 * math.sin(a)) for a in ang]
                ax.add_patch(Polygon(inn, closed=True, fc=ODD if d % 2 else "white", ec=INK, lw=0.9, zorder=1))
                for p, q in zip(inn, out):
                    ax.plot([p[0], q[0]], [p[1], q[1]], color=INK, lw=0.9)
                for p in inn:
                    ax.add_patch(Circle(p, 0.028, fc="white", ec=INK, lw=0.7, zorder=3))
                ax.text(cx, y0 - R - 0.06, "even face" if d % 2 == 0 else "odd face", ha="center",
                        va="top", fontsize=7.5)
        arrow(ax, 0.62, y0, 0.84, y0)
    panel_label(ax, "(a)", x=0.0, y=1.08)

    lat = [("(b)", "truncated Penrose", lambda: truncated_penrose(3.2, seed=0)),
           ("(c)", "Voronoi foam", lambda: random_cubic(45, 5)),
           ("(d)", "star lattice", lambda: star_lattice(5))]
    stats = {}
    for k, (lab, name, fn) in enumerate(lat):
        ax = fig.add_axes([0.215 + k * 0.265, 0.10, 0.25, 0.80])
        n, E, pos, info = fn()
        P = np.asarray(pos, float)
        fr, s = healed_state(n, E)
        faces = geometric_faces(n, E, [tuple(p) for p in P])
        odd = [P[f] for f in faces if len(f) % 2 == 1]
        ax.add_collection(PolyCollection(odd, facecolors=ODD, edgecolors="none", zorder=0))
        draw_bonds(ax, P, E, s, lw=0.55)
        sc = 2.6 if n > 150 else 3.2
        draw_sites(ax, P, s, ms=sc, big=1.45, lw=0.3)
        pad = 0.03 * np.ptp(P, axis=0).max()
        ax.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad)
        ax.set_ylim(P[:, 1].min() - pad, P[:, 1].max() + pad)
        ax.set_aspect("equal", adjustable="datalim"); ax.axis("off")
        n3 = sum(1 for x in s if x == 3)
        stats[name] = (n, fr, n3)
        panel_label(ax, lab, x=0.0, y=1.08)
        ax.text(1.0, 1.08, f"{name}\n$N={n}$, $\\mathrm{{fr}}={fr}$", transform=ax.transAxes, ha="right",
                va="top", fontsize=7.5, linespacing=1.2)
    state_legend(fig, "lower center", ncol=5, frus=False, odd=True, bbox=(0.6, -0.04))
    save(fig, "fig_truncation.pdf")
    return stats


# ================================================================== 7. timing + table
FAMNAME = {"truncated_penrose": "truncated Penrose", "voronoi_foam": "Voronoi foam",
           "defect_honeycomb": "defect honeycomb"}
FAMSTYLE = {"truncated_penrose": ("o", SC[1]), "voronoi_foam": ("s", "#D55E00"),
            "defect_honeycomb": ("D", "#009E73")}


def load_timing():
    p = os.path.join(ROOT, "hub_runs", "timing_cert.jsonl")
    rows = []
    for line in open(p):
        line = line.strip()
        if line:
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass                    # a partially written last line
    return rows


def fig_timing(rows):
    fig = plt.figure(figsize=(COL1, 2.35))
    ax = fig.add_axes([0.15, 0.17, 0.82, 0.80])
    fam_h = []
    for fam in FAMNAME:
        R = sorted((r for r in rows if r["fam"] == fam), key=lambda r: r["n"])
        if not R:
            continue
        mk, c = FAMSTYLE[fam]
        n = [r["n"] for r in R]
        ax.plot(n, [r["t_match"] for r in R], marker=mk, color=c, ms=3.8, mec=c, lw=1.0)
        R2 = [r for r in R if r["t_2sat"] > 0]          # 0.0 = below the 1 ms resolution
        ax.plot([r["n"] for r in R2], [r["t_2sat"] for r in R2], marker=mk, color=c, ms=3.8,
                mfc="white", mec=c, lw=0.9, ls="--")
        fam_h.append(Line2D([], [], marker=mk, color=c, ms=3.8, lw=1.0, label=FAMNAME[fam]))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("number of sites $N$"); ax.set_ylabel("time (s)")
    ax.set_ylim(3e-4, 1e3)
    ax.set_xticks([200, 500, 1000, 2000, 5000])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    l1 = ax.legend(handles=fam_h, loc="upper left", frameon=False, handlelength=2.0, borderaxespad=0.3)
    ax.add_artist(l1)
    step_h = [Line2D([], [], color="0.3", marker="o", ms=3.8, lw=1.0, label="matching"),
              Line2D([], [], color="0.3", marker="o", ms=3.8, mfc="white", ls="--", lw=0.9, label="2-SAT")]
    ax.legend(handles=step_h, loc="center right", frameon=False, handlelength=2.4, borderaxespad=0.3)
    save(fig, "fig_timing.pdf")


def timing_table_tex(rows):
    fams = [f for f in FAMNAME if any(r["fam"] == f for r in rows)]
    L = [r"\begin{table}",
         r"\caption{Certification algorithm on coordination-three lattices: number of sites $N$,",
         r"Ising frustration $\fr$, wall-clock times of the matching step (Ising ground state as a",
         r"minimum $T$-join) and of the 2-SAT step (choice of endpoints), and whether the first Ising",
         r"ground state was healed, which certifies $E_0(D_3)=\min(D_3,1)\fr$. Python/NetworkX implementation.}",
         r"\label{tab:timing}",
         r"\begin{ruledtabular}",
         r"\begin{tabular}{lrrrrc}",
         r"lattice & $N$ & $\fr$ & matching (s) & 2-SAT (ms) & healed\\",
         r"\hline"]
    for f in fams:
        for r in sorted((r for r in rows if r["fam"] == f), key=lambda r: r["n"]):
            t2 = f"{1000 * r['t_2sat']:.0f}" if r["t_2sat"] > 0 else "$<1$"
            ok = "yes" if r["healed"] else "no"
            L.append(f"{FAMNAME[f]} & {r['n']} & {r['fr']} & {r['t_match']:.2f} & {t2} & {ok}\\\\")
    L += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    return "\n".join(L)


def update_main_tex(table):
    p = os.path.join(HERE, "main.tex")
    src = open(p).read()
    a, b = "% BEGIN generated timing table (make_figs.py)", "% END generated timing table"
    if a not in src or b not in src:
        print("timing-table markers not found in main.tex; table not written")
        return
    pat = re.compile(re.escape(a) + r".*?" + re.escape(b), re.S)
    new = pat.sub(lambda m: a + "\n" + table + "\n" + b, src)
    if new != src:
        open(p, "w").write(new)
        print("updated timing table in main.tex")


if __name__ == "__main__":
    which = set(sys.argv[1:])
    run = lambda k: not which or k in which
    if run("recolour"):
        fig_recolour()
    if run("wheel"):
        fig_wheel()
    if run("triangular"):
        fig_triangular()
    if run("penrose"):
        fig_penrose()
    if run("healing"):
        print("healing", fig_healing())
    if run("truncation"):
        print("truncation", fig_truncation())
    if run("timing"):
        rows = load_timing()
        fig_timing(rows)
        update_main_tex(timing_table_tex(rows))
        print("timing rows:", len(rows))

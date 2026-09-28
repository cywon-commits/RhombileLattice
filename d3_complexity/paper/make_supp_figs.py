"""Figure for Supplementary Note S1 (run from paper/):  python3 make_supp_figs.py

fig_supp_minor.pdf: (a) the reduction graph G(phi*) of Ref. [Johnson2025], Thm 11, for
phi* = (x or y)(not x or y)(x or not y), built by verify_delta3_reduction.build; the shaded
sets are contracted (path P -> p, clause gadget j -> C_j), unused inputs u_j deleted.
(b) the resulting minor: the incidence graph of phi* plus an apex p adjacent to every clause,
which is K_{3,3} with parts {C1, C2, C3} and {p, x, y}.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from verify_delta3_reduction import build
import networkx as nx

plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 8, "pdf.fonttype": 42,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
INK, BOND = "0.15", "0.45"
GADGET, PATH, LIT = "#E69F00", "#0072B2", "#F0E442"

PHI = [[1, 2], [-1, 2], [1, -2]]


def fig_minor():
    n, E, idx = build(2, PHI)
    inv = {v: k for k, v in idx.items()}
    m = len(PHI)
    pos = {}
    X = [0.9, 2.5, 4.1]                                   # gadget centres
    lit_x = {(0, 1): 0.9, (0, 0): 1.9, (1, 1): 3.1, (1, 0): 4.1}
    for key, x in lit_x.items():
        pos[("x",) + key] = (x, 3.05)
    tmpl = {"a": (-0.35, 2.35), "b": (0.35, 2.35), "d": (0.0, 1.95),
            "e": (0.0, 1.55), "f": (0.42, 1.2), "c": (-0.3, 1.1)}
    for j in range(m):
        for t, (dx, y) in tmpl.items():
            pos[("g", j, t)] = (X[j] + dx, y)
        pos[("u", j)] = (X[j] + 0.75, 1.55)
    for k in range(2 * m):
        pos[("P", k)] = (0.55 + 0.72 * k, 0.35)
    P = {v: pos[inv[v]] for v in range(n)}
    G = nx.Graph(); G.add_edges_from(E)
    _, K = nx.check_planarity(G, counterexample=True)
    Kset = {tuple(sorted(e)) for e in K.edges()}

    fig = plt.figure(figsize=(7.0, 2.35))
    ax = fig.add_axes([0.0, 0.0, 0.62, 1.0]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(0.1, 5.2); ax.set_ylim(0.0, 3.45)
    for j in range(m):   # contraction sets
        ax.add_patch(FancyBboxPatch((X[j] - 0.55, 0.92), 1.1, 1.62, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=GADGET, alpha=0.18, ec="none", zorder=0))
        ax.text(X[j] - 0.47, 1.75, f"$Q_{j+1}$", fontsize=8, color="#8a5a00", va="center")
    ax.add_patch(FancyBboxPatch((0.35, 0.18), 3.95, 0.34, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=PATH, alpha=0.15, ec="none", zorder=0))
    ax.text(4.42, 0.35, "$P\\to p$", fontsize=8, va="center", color="#00507d")
    for u, v in E:
        (x1, y1), (x2, y2) = P[u], P[v]
        inK = tuple(sorted((u, v))) in Kset
        ax.plot([x1, x2], [y1, y2], color="#A50F15" if inK else BOND, lw=1.5 if inK else 0.6, zorder=1)
    for v in range(n):
        key = inv[v]; x, y = P[v]
        fc = LIT if key[0] == "x" else ("white" if key[0] in ("u",) else "0.85")
        ax.add_patch(Circle((x, y), 0.07, fc=fc, ec=INK, lw=0.6, zorder=3))
    for (i, s), x in lit_x.items():
        name = "xy"[i] if s == 1 else "\\bar{" + "xy"[i] + "}"
        ax.text(x, 3.22, f"${name}$", ha="center", va="bottom", fontsize=8)
    for j in range(m):
        ax.text(X[j] + 0.75, 1.72, f"$u_{j+1}$", ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.text(0.12, 3.42, "(a)", fontsize=9, va="top")

    # (b) K_{3,3}
    ax2 = fig.add_axes([0.66, 0.08, 0.33, 0.84]); ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_xlim(-0.3, 2.9); ax2.set_ylim(-0.35, 2.45)
    top = {"C_1": (0.2, 1.9), "C_2": (1.3, 1.9), "C_3": (2.4, 1.9)}
    bot = {"p": (0.2, 0.3), "x": (1.3, 0.3), "y": (2.4, 0.3)}
    for a in top.values():
        for b in bot.values():
            ax2.plot([a[0], b[0]], [a[1], b[1]], color=INK, lw=0.9, zorder=1)
    for lab, (x, y) in top.items():
        ax2.add_patch(Circle((x, y), 0.14, fc=GADGET, ec=INK, lw=0.6, zorder=3))
        ax2.text(x, y + 0.22, f"${lab}$", ha="center", va="bottom", fontsize=8)
    for lab, (x, y) in bot.items():
        ax2.add_patch(Circle((x, y), 0.14, fc=PATH if lab == "p" else LIT, ec=INK, lw=0.6, zorder=3))
        ax2.text(x, y - 0.22, f"${lab}$", ha="center", va="top", fontsize=8)
    ax2.text(-0.28, 2.43, "(b)", fontsize=9, va="top")
    fig.savefig(os.path.join(HERE, "figs", "fig_supp_minor.pdf"))
    print("wrote figs/fig_supp_minor.pdf; Kuratowski subgraph:", K.number_of_nodes(), "vertices")


if __name__ == "__main__":
    fig_minor()

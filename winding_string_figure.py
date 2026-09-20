"""Figure for the direct-vs-winding flux-sector experiment
(winding_string_experiment2.py): same two defect sites on the same 30x8
torus, connected either directly (never touches the periodic seam) or by
deliberately routing through it (winds once). Shows the two lattices side
by side plus a small energy-vs-D3 summary panel using the exact (mincut/
best-SA) values already computed."""
import pickle

import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import RhombileLattice, apply_string_defect, frustrated_triangles
from dual_string_demo import draw_lattice

NX, NY = 30, 8
Y_CROSS = 2.1
P_START = np.array([2.0, Y_CROSS])
P_END = np.array([27.0, Y_CROSS])


def build(wrap):
    lat = RhombileLattice(NX, NY)
    flipped, p_end_used = apply_string_defect(lat, P_START, P_END, wrap=wrap)
    on = [b for b in flipped if b["J"] != 0.0]
    off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    return lat, on, off, ft, p_end_used


def main():
    with open("winding_string_results2.pkl", "rb") as f:
        results = pickle.load(f)

    lat_d, on_d, off_d, ft_d, p_end_d = build(wrap=False)
    lat_w, on_w, off_w, ft_w, p_end_w = build(wrap=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    ax_d, ax_w, ax_bar, ax_text = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    draw_lattice(ax_d, lat_d, highlight_on=on_d, highlight_off=off_d, frustrated=ft_d,
                 box=(NX, NY), targets=[P_START, P_END], off_lw=0.6,
                 title=f"DIRECT: straight through the middle\nL=25, never touches the periodic seam")
    ax_d.set_xlim(-1, NX)

    draw_lattice(ax_w, lat_w, highlight_on=on_w, highlight_off=off_w, frustrated=ft_w,
                 box=(NX, NY), targets=[P_START, P_END], off_lw=0.6,
                 title=f"WINDING: same 2 defect sites, routed through the seam instead\n"
                       f"L=5, wraps once (shown near the box's right/left edge)")
    ax_w.set_xlim(-1, NX)

    D3_vals = sorted(results.keys())
    e_d = [results[d]["E_direct"] for d in D3_vals]
    e_w = [results[d]["E_winding"] for d in D3_vals]
    src_d = [results[d]["src_direct"] for d in D3_vals]
    src_w = [results[d]["src_winding"] for d in D3_vals]
    x = np.arange(len(D3_vals))
    w = 0.35
    ax_bar.bar(x - w / 2, e_d, width=w, label="direct (L=25)", color="tab:orange")
    ax_bar.bar(x + w / 2, e_w, width=w, label="winding (L=5, wraps once)", color="tab:blue")
    for xi, e, s in zip(x - w / 2, e_d, src_d):
        ax_bar.text(xi, e + 1, s, ha="center", fontsize=7, rotation=90, va="bottom")
    for xi, e, s in zip(x + w / 2, e_w, src_w):
        ax_bar.text(xi, e + 1, s, ha="center", fontsize=7, rotation=90, va="bottom")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([f"D3={d}" for d in D3_vals])
    ax_bar.set_ylabel("best energy found\n(mincut-seed + gentle re-anneal + 10 random restarts)")
    ax_bar.set_title("Same 2 defect sites, 2 topologically different\nconnecting paths: energy never converges")
    ax_bar.legend(fontsize=8)

    ax_text.axis("off")
    text = (
        "Setup: two fixed lattice sites at x=2 and x=27 on a 30-wide periodic (torus) lattice.\n"
        "  DIRECT distance (through the middle) = 25 unit cells.\n"
        "  WINDING distance (through the periodic seam) = 30-25 = 5 unit cells.\n\n"
        "Both constructions isolate exactly 2 real topological defects at the same 2 sites\n"
        "(verified: frustrated_triangles == 2 for both). Only WINDING's flipped-bond set\n"
        "includes an actual wrap-flagged bond (verified directly).\n\n"
        "Result: extensive SA (exact mincut seed + gentle re-anneal + 10 random restarts,\n"
        "4000 sweeps each) never closes the gap. DIRECT stays near its own L=25 scaling\n"
        "(25 -> 48, a small local correction below the naive hub-fixed value of 50),\n"
        "WINDING stays near its own L=5 scaling (5 -> 10). Local Monte Carlo dynamics\n"
        "cannot move between the two: exactly the flip-connectivity-only-within-a-flux-\n"
        "sector picture from Thurston's theorem for dimer coverings on a torus.\n\n"
        "Caveat: this is numerical (SA never finding a bridge), not a rigorous inescapability\n"
        "proof -- but it is the expected signature if the sectors are genuinely disconnected\n"
        "under any local (single-site) update, matching the theory."
    )
    ax_text.text(0.02, 0.98, text, va="top", ha="left", fontsize=9.5, family="monospace",
                 transform=ax_text.transAxes)

    plt.suptitle("Same two defects, two flux sectors: winding beats direct and SA can't bridge them",
                  y=1.01, fontsize=13)
    plt.tight_layout()
    fig.savefig("winding_string_comparison.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_string_comparison.png")


if __name__ == "__main__":
    main()

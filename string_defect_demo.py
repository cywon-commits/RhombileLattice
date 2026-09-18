"""String-defect experiments on the rhombile lattice: straight-line J
dilution/reconnection, the resulting topologically-frustrated triangles,
SA-based defect reconnection, and the D3 (state-2 anisotropy) threshold.

Convention: H = sum_<i,j> (-J_ij) delta(s_i, s_j) + sum_i D[s_i], AF for the
lattice's default J < 0. See rhombile_lattice.py docstrings for details.
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, frustrated_triangles, matched_bonds,
    simulated_annealing, total_energy,
)

NX, NY = 10, 10


def draw_lattice(ax, lattice, highlight_bonds=None, frustrated=None,
                  states=None, title=""):
    """Shared plotting helper: black = ordinary active bonds, red = bonds
    in `highlight_bonds`, gold stars = triangles in `frustrated`, and (if
    given) site fill color = spin state (white/black/gray for 0/1/2)."""
    highlight_ids = set(id(b) for b in highlight_bonds) if highlight_bonds else set()
    for b in lattice.bonds:
        if b["J"] != 0.0:
            is_hl = id(b) in highlight_ids
            color = "red" if is_hl else "black"
            lw = 2.2 if is_hl else 0.8
            ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                     color=color, linewidth=lw, zorder=2 if is_hl else 1)

    pos, sub_of = lattice.site_positions()
    if states is None:
        colors = {"r1": "tab:red", "r2": "tab:blue", "r3": "tab:green"}
        for sub, c in colors.items():
            mask = sub_of == sub
            ax.scatter(pos[mask, 0], pos[mask, 1], s=22, color=c, zorder=3)
    else:
        state_color = {0: "white", 1: "black", 2: "tab:orange"}
        for k, c in state_color.items():
            mask = states == k
            ax.scatter(pos[mask, 0], pos[mask, 1], s=35, color=c,
                       edgecolor="gray", linewidth=0.5, zorder=3)

    if frustrated:
        centroids = np.array([f["centroid"] for f in frustrated])
        ax.scatter(centroids[:, 0], centroids[:, 1], marker="*", s=260,
                   color="gold", edgecolor="black", linewidth=0.8, zorder=4)

    ax.set_aspect("equal")
    ax.set_xlim(-1, NX - 1)
    ax.set_ylim(-1, NY - 1)
    ax.set_title(title)


def single_string_demo(y_cross=2.1, x0=0.5, x1=7.5):
    """One horizontal string: shows the 'connecting' (0->-1) vs 'breaking'
    (-1->0) bond counts, and the 2 topologically-frustrated triangles that
    end up exactly at the string's two ends."""
    lat = RhombileLattice(NX, NY)
    flipped, _ = apply_string_defect(lat, np.array([x0, y_cross]), np.array([x1, y_cross]), wrap=False)
    on = [b for b in flipped if b["J"] != 0.0]
    off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    print(f"[single string, y={y_cross}] connecting(0->-1)={len(on)}  "
          f"breaking(-1->0)={len(off)}  frustrated triangles={len(ft)}")

    fig, ax = plt.subplots(figsize=(9, 8))
    draw_lattice(ax, lat, highlight_bonds=on, frustrated=ft,
                 title=f"single string y={y_cross}: {len(ft)} topological defects")
    plt.tight_layout()
    fig.savefig("string_defect_single.png", dpi=150)
    plt.close(fig)
    return lat


def two_parallel_strings_demo(y1=2.1, y2=3.0, x0=0.5, x1=7.5, n_trials=8, n_sweeps=1200):
    """Two parallel strings: draws 4 topological defects (2 pairs stacked
    vertically), then runs 2-state-forced SA to show the true ground state
    reconnects the *near* (vertical) pairs instead of following each
    string's own long horizontal path -- a minimum-weight-matching effect."""
    lat = RhombileLattice(NX, NY)
    f1, _ = apply_string_defect(lat, np.array([x0, y1]), np.array([x1, y1]), wrap=False)
    f2, _ = apply_string_defect(lat, np.array([x0, y2]), np.array([x1, y2]), wrap=False)
    on = [b for b in f1 + f2 if b["J"] != 0.0]
    ft = frustrated_triangles(lat)
    print(f"[two strings, y={y1},{y2}] frustrated triangles={len(ft)}")

    D = (0.0, 0.0, 1e6)  # force 2-state (Ising) sector
    rng = np.random.default_rng(0)
    best_states, best_E = None, np.inf
    for _ in range(n_trials):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps,
                                                 T_start=5.0, T_end=1e-4, record_energy=True)
        if energies[-1] < best_E:
            best_E, best_states = energies[-1], states
    print(f"[two strings] SA ground energy = {best_E} "
          f"(naive independent-strings estimate: ~2x single-string E)")

    mb = matched_bonds(lat, best_states)
    fig, ax = plt.subplots(figsize=(9, 8))
    draw_lattice(ax, lat, highlight_bonds=on, frustrated=ft, states=best_states,
                 title=f"two strings, SA ground state E={best_E}\nmagenta-eligible bonds shown below are the real cost")
    for b in mb:
        ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]], color="magenta", linewidth=3.0, zorder=5)
    plt.tight_layout()
    fig.savefig("string_defect_two_reconnected.png", dpi=150)
    plt.close(fig)
    return lat, best_E


def single_triangle_D3_threshold():
    """Exact brute-force threshold for when the 3rd Potts state becomes
    favorable on one isolated AF triangle (3 sites, pairwise |J|=1):
    D_c = |J| = 1. Below it, using state 2 on one corner (making its two
    bonds delta=0 for free) beats the best achievable 2-state assignment."""
    def triangle_energy(states, D):
        e = sum(D[s] for s in states)
        for a, b in [(states[0], states[1]), (states[1], states[2]), (states[0], states[2])]:
            if a == b:
                e += 1.0
        return e

    print("\n[single-triangle D3 threshold]")
    for D3 in [2.0, 1.5, 1.001, 1.0, 0.999, 0.5, 0.0]:
        D = (0.0, 0.0, D3)
        configs = list(product(range(3), repeat=3))
        best = min(triangle_energy(s, D) for s in configs)
        uses_state2 = any(2 in s for s in configs if triangle_energy(s, D) == best)
        print(f"  D3={D3:6.3f}  min_E={best:.3f}  uses_state2={uses_state2}")


def d3_sweep_strings(D3_values=(1e6, 1.0, 0.5), n_trials=10, n_sweeps=1200):
    """Full 3-state Potts SA on the single- and double-string lattices,
    scanning D3 down through and below the isolated-triangle threshold."""
    strings_single = [((0.5, 2.1), (7.5, 2.1))]
    strings_double = [((0.5, 2.1), (7.5, 2.1)), ((0.5, 3.0), (7.5, 3.0))]

    def run(D, strings, rng):
        lat = RhombileLattice(NX, NY)
        for p0, p1 in strings:
            apply_string_defect(lat, np.array(p0), np.array(p1), wrap=False)
        best_E, best_states = np.inf, None
        for _ in range(n_trials):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps,
                                                     T_start=4.0, T_end=1e-3, record_energy=True)
            if energies[-1] < best_E:
                best_E, best_states = energies[-1], states
        return best_E, float(np.mean(best_states == 2))

    print("\n[D3 sweep: single- vs double-string energy and state-2 usage]")
    rng = np.random.default_rng(42)
    for D3 in D3_values:
        D = (0.0, 0.0, D3)
        E1, f1 = run(D, strings_single, rng)
        E2, f2 = run(D, strings_double, rng)
        print(f"  D3={D3:8.3g}:  single E={E1:6.2f} (state2 frac={f1:.3f})   "
              f"double E={E2:6.2f} (state2 frac={f2:.3f})")


if __name__ == "__main__":
    single_string_demo()
    two_parallel_strings_demo()
    single_triangle_D3_threshold()
    d3_sweep_strings()

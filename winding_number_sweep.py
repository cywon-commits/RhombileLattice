"""Generalizes winding_string_experiment2.py from a single "direct (n=0) vs
winds-once (n=-1)" comparison to several winding numbers, to check whether
the sector-separation seen there is a one-off or a genuine pattern: for the
SAME two physical defect sites (x=2 and x=27 on a 30-wide torus, y=2.1),
build a straight string that goes around the periodic x-direction exactly
n times before reaching the target (n=0 is the plain direct path; n=-1 is
the earlier "winding" path; n=+1, -2, +2, ... wind further).

Since minimum_image_endpoint always folds to the *shortest* representative
(so plain apply_string_defect(wrap=True) can only ever give n=0 or the one
minimal wind), each n>{-1,0} here is built by hand: target = the true site
plus n full box-widths, fed through the same shift-tested bond-crossing
logic apply_string_defect uses internally, just with a wider shift range
and without the automatic minimum-imaging step.

Expected pattern if the earlier result generalizes: each winding number is
its own topological sector (not reachable from another via local SA moves)
whose floor energy tracks its own string length L(n) = |25 + 30n|, so
energy vs n should look like several straight-ish branches with a minimum
at n=-1 (L=5), not a single curve that all collapse onto.

SUPERSEDED / INCOMPLETE: n=-2 fails its own sanity check (0 defects
instead of 2). Root cause, found by hand-tracing a 2-row staircase
variant afterwards: "the string" is not an object intrinsic to a
covering, only its symmetric difference against a chosen reference
background is -- so trying to dial in a specific winding number by
construction procedure (how many times a hand-drawn line loops) is
fragile: a >1-box-width single-row line starts re-toggling the same
finite set of ~NX bonds it already touched, XOR-canceling itself, and a
naive multi-row "staircase" isn't safe either (apply_string_defect's
wrap=True tests the full 2D (kx,ky) shift grid, not just x). The
principled fix is to measure flux directly off the FINAL bond/state
configuration (a gauge-invariant quantity, like a height-function
winding or a whole-graph 2-coloring role-swap count along a fixed
non-contractible loop) rather than trust that a construction procedure
achieved the intended winding number. See flux_measure.py.
"""
import pickle

import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, A1, segments_intersect, frustrated_triangles,
    simulated_annealing, total_energy, full_state_via_mincut,
)
from dual_string_demo import draw_lattice

NX, NY = 30, 8
Y_CROSS = 2.1
P_START = np.array([2.0, Y_CROSS])
P_END_TRUE = np.array([27.0, Y_CROSS])

N_WINDS = [-2, -1, 0, 1, 2]
D3_GRID = [1.0, 2.0, 3.0]
RESULT_PATH = "winding_number_sweep_results.pkl"


def build_winding(n):
    """n full box-widths added to the true target before drawing a plain
    straight segment; toggles every real bond crossed by ANY periodic
    translate of that segment, exactly like apply_string_defect's own
    wrap=True path but skipping its minimum-imaging step so a specific,
    possibly non-minimal, winding number can be forced."""
    lat = RhombileLattice(NX, NY)
    p_end = P_END_TRUE + n * NX * A1
    shift_range = range(-(abs(n) + 2), 3)
    shifts = [kx * NX * A1 for kx in shift_range]  # only x wraps matter (string stays at y=Y_CROSS)
    flipped = []
    for shift in shifts:
        s_start, s_end = P_START - shift, p_end - shift
        for b in lat.bonds:
            if segments_intersect(s_start, s_end, b["p1"], b["p2"]):
                b["J"] = 0.0 if b["J"] != 0.0 else -1.0
                flipped.append(b)
    L = abs(p_end[0] - P_START[0])
    return lat, flipped, p_end, L


def best_energy(lat, D3, seed, n_random=8, n_sweeps_random=3500,
                 n_seeded_gentle=5, n_sweeps_gentle=3500):
    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(seed)
    candidates = []
    try:
        seed_states, _ = full_state_via_mincut(lat, D3)
        candidates.append((total_energy(lat, seed_states, D), seed_states.copy()))
        for _ in range(n_seeded_gentle):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_gentle,
                                                     T_start=0.5, T_end=1e-6,
                                                     states=seed_states, record_energy=True)
            candidates.append((min(energies), states))
    except ValueError:
        pass
    for _ in range(n_random):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_random,
                                                 T_start=5.0, T_end=1e-6, record_energy=True)
        candidates.append((min(energies), states))
    best_e, best_states = min(candidates, key=lambda c: c[0])
    return best_e, best_states


def main():
    print(f"{'n':>3}  {'L':>5}  {'#flipped':>9}  {'#defects':>9}")
    builds = {}
    for n in N_WINDS:
        lat, flipped, p_end, L = build_winding(n)
        ft = frustrated_triangles(lat)
        print(f"{n:>3}  {L:>5.0f}  {len(flipped):>9d}  {len(ft):>9d}")
        assert len(ft) == 2, f"n={n}: expected 2 defects, got {len(ft)}"
        builds[n] = {"lat": lat, "flipped": flipped, "p_end": p_end, "L": L}

    results = {"L": {n: builds[n]["L"] for n in N_WINDS}, "D3": {}}
    print()
    header = "  D3  " + "  ".join(f"n={n:>3}" for n in N_WINDS)
    print(header)
    for D3 in D3_GRID:
        row_e = {}
        for n in N_WINDS:
            e, states = best_energy(builds[n]["lat"], D3, seed=1000 + int(D3 * 10) + n)
            row_e[n] = e
        print(f"{D3:>5.1f}  " + "  ".join(f"{row_e[n]:>6.1f}" for n in N_WINDS))
        results["D3"][D3] = row_e
        with open(RESULT_PATH, "wb") as f:
            pickle.dump(results, f)

    # figure: energy vs L, one line per D3, plus lattice panels for the
    # extreme cases (n=-2 and n=+2) to show the multi-wrap string visually
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    ax_curve = axes[0, 0]
    Ls = [builds[n]["L"] for n in N_WINDS]
    order = np.argsort(Ls)
    for D3 in D3_GRID:
        es = [results["D3"][D3][n] for n in N_WINDS]
        Ls_sorted = [Ls[i] for i in order]
        es_sorted = [es[i] for i in order]
        ax_curve.plot(Ls_sorted, es_sorted, "o-", label=f"D3={D3}")
        for n, L, e in zip(N_WINDS, Ls, es):
            ax_curve.annotate(f"n={n}", (L, e), textcoords="offset points",
                               xytext=(0, 6), fontsize=8, ha="center")
    ax_curve.set_xlabel("string length L(n) = |25 + 30n|")
    ax_curve.set_ylabel("best energy found")
    ax_curve.set_title("Energy vs winding number: each n sits on its own\n"
                        "length-dependent branch, none collapse together")
    ax_curve.legend()

    ax_text = axes[0, 1]
    ax_text.axis("off")
    lines = ["n    L    D3=1   D3=2   D3=3"]
    for n in N_WINDS:
        lines.append(f"{n:>2}  {builds[n]['L']:>4.0f}   " +
                      "  ".join(f"{results['D3'][D3][n]:>5.1f}" for D3 in D3_GRID))
    ax_text.text(0.05, 0.9, "\n".join(lines), family="monospace", fontsize=11,
                 va="top", transform=ax_text.transAxes)
    ax_text.set_title("Same two physical defect sites, different winding numbers")

    for ax, n in zip([axes[1, 0], axes[1, 1]], [-2, 2]):
        lat, flipped, p_end, L = builds[n]["lat"], builds[n]["flipped"], builds[n]["p_end"], builds[n]["L"]
        on = [b for b in flipped if b["J"] != 0.0]
        off = [b for b in flipped if b["J"] == 0.0]
        ft = frustrated_triangles(lat)
        draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                     box=(NX, NY), targets=[P_START, P_END_TRUE], off_lw=0.4,
                     title=f"n={n} (winds {abs(n)}x), L={L:.0f}: same 2 defect sites, drawn\n"
                           f"across {abs(n)} extra trip(s) around the periodic x-direction")
        ax.set_xlim(-1, NX)

    plt.tight_layout()
    fig.savefig("winding_number_sweep.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_number_sweep.png")


if __name__ == "__main__":
    main()

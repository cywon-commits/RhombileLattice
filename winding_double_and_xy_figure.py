"""Figure for winding_double_and_xy.py's two tests, colored by actual
Potts state (not species) per the user's request:

  Left panel: TWO independent x-winding loops (rows m=2 and m=4) combine
  to a globally consistent 2-coloring (0 whole-graph contradictions) --
  shown via the actual free BFS 2-coloring (pure state1/state2, exactly
  like Case0), not mincut (known unreliable for non-straight geometries).

  Right panel: ONE x-winding loop + ONE y-winding loop together still
  have no consistent 2-coloring -- shown via the SA-best state at D3=1,
  which must use some state-3 sites to resolve the remaining obstruction.

  Bottom text panel also documents why the raw "number of contradictions"
  (72 vs 60 in the original run) is not itself a meaningful quantity: it
  depends on the arbitrary BFS traversal order (verified: the SAME x-alone
  lattice gives 71-79 across different random start orders). Only whether
  the count is exactly zero is order-independent.
"""
import pickle
from collections import deque

import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, frustrated_triangles, full_state_via_mincut, total_energy,
    simulated_annealing,
)
from dual_string_demo import draw_lattice
from winding_closed_loop import enumerate_triangles_all
from winding_double_and_xy import add_winding_loop, build_base, whole_graph_contradictions, NX, NY, N_WAYPOINTS


def free_two_coloring(lat, start_site=0):
    """Same idea as closed_loop_demo._bfs_two_coloring: BFS parity over
    the whole active-bond graph. Valid (0 contradictions) whenever the
    graph is bipartite -- used here for the double-x-wind case, which we
    just confirmed has a consistent global coloring."""
    active_bonds = [b for b in lat.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        adj.setdefault(b["i"], []).append(b["j"])
        adj.setdefault(b["j"], []).append(b["i"])
    color = {start_site: 0}
    q = deque([start_site])
    while q:
        u = q.popleft()
        for v in adj.get(u, []):
            if v not in color:
                color[v] = 1 - color[u]
                q.append(v)
    return np.array([color.get(s, 0) for s in range(lat.n_sites)])


def build_double_x():
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1, rhombi, sibling, hop = build_base(lat)
    ns = list(np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int) % NX)
    touched = []
    touched += add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, [(int(n), 2) for n in ns])
    touched += add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, [(int(n), 4) for n in ns])
    return lat, touched


def build_xy_combined():
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1, rhombi, sibling, hop = build_base(lat)
    ns = list(np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int) % NX)
    ms = list(np.linspace(0, NY, N_WAYPOINTS, endpoint=False).round().astype(int) % NY)
    touched = []
    touched += add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, [(int(n), 2) for n in ns])
    touched += add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, [(15, int(m)) for m in ms])
    return lat, touched


def main():
    lat1, touched1 = build_double_x()
    on1 = [b for b in touched1 if b["J"] != 0.0]
    off1 = [b for b in touched1 if b["J"] == 0.0]
    ft1 = frustrated_triangles(lat1)
    contras1 = whole_graph_contradictions(lat1)
    states1 = free_two_coloring(lat1)
    e1 = total_energy(lat1, states1, D=(0.0, 0.0, 1.0))
    print(f"double x-wind: defects={len(ft1)}, contradictions={contras1}, "
          f"free-2-coloring E(D3=1)={e1}, n_state3={(states1==2).sum()}")

    lat2, touched2 = build_xy_combined()
    on2 = [b for b in touched2 if b["J"] != 0.0]
    off2 = [b for b in touched2 if b["J"] == 0.0]
    ft2 = frustrated_triangles(lat2)
    contras2 = whole_graph_contradictions(lat2)
    D3_plot = 1.0
    seed_states, _ = full_state_via_mincut(lat2, D3_plot)
    best_e, best_states = total_energy(lat2, seed_states, D=(0.0, 0.0, D3_plot)), seed_states
    rng = np.random.default_rng(3)
    for _ in range(5):
        st, en = simulated_annealing(lat2, (0.0, 0.0, D3_plot), rng, n_sweeps=2500,
                                      T_start=0.5, T_end=1e-6, states=seed_states, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    for _ in range(4):
        st, en = simulated_annealing(lat2, (0.0, 0.0, D3_plot), rng, n_sweeps=3000,
                                      T_start=5.0, T_end=1e-6, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    print(f"x+y combined: defects={len(ft2)}, contradictions={contras2}, "
          f"SA-best E(D3=1)={best_e}, n_state3={(best_states==2).sum()}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    ax1, ax2, ax_text, ax_blank = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    draw_lattice(ax1, lat1, highlight_on=on1, highlight_off=off1, frustrated=ft1,
                 box=(NX, NY), off_lw=0.4, states=states1,
                 title=f"TWO independent x-winding loops (rows m=2, m=4)\n"
                       f"colored by state -- free 2-coloring gives E=0 (0 state-3 sites)")
    handles, labels = ax1.get_legend_handles_labels()
    label_map = {"state 0": "state1", "state 1": "state2", "state 2": "state3"}
    ax1.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=8)

    draw_lattice(ax2, lat2, highlight_on=on2, highlight_off=off2, frustrated=ft2,
                 box=(NX, NY), off_lw=0.4, states=best_states,
                 title=f"x-winding (row m=2) + y-winding (column n=15) loops\n"
                       f"colored by state, D3=1 -- SA-best E={best_e:.0f}, "
                       f"{int((best_states==2).sum())} state-3 sites needed")
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=8)

    ax_text.axis("off")
    text = (
        "Q: why did x+y combined (60 contradictions) show FEWER than x-alone (72)?\n"
        "A: that raw count is NOT a meaningful quantity -- it depends on the arbitrary\n"
        "BFS traversal order used to attempt the 2-coloring. Verified directly: re-running\n"
        "the SAME x-alone lattice with 5 different random BFS start orders gave 71, 75,\n"
        "77, 79, 79 instead of 72 -- a spread comparable to the 72-vs-60 gap itself.\n"
        "The only order-INDEPENDENT fact is whether the count is exactly zero:\n"
        "  - double x-wind: 0 for every order tried -> genuinely trivial (matches Z2).\n"
        "  - x+y combined: nonzero for every order tried -> genuinely obstructed,\n"
        "    but the specific nonzero value carries no further physical meaning.\n\n"
        "Left panel: the double-x-wind lattice's actual free 2-coloring (pure\n"
        "state1/state2, the same mechanism as Case0) -- 0 state-3 sites, E=0 at every D3.\n\n"
        "Right panel: the x+y-combined lattice's real SA-best ground state at D3=1 --\n"
        "state-3 sites (orange) cluster near where the two loops are, resolving the\n"
        "genuine (1,1) flux obstruction that a free 2-coloring cannot."
    )
    ax_text.text(0.02, 0.98, text, va="top", ha="left", fontsize=9.3, family="monospace",
                 transform=ax_text.transAxes)
    ax_blank.axis("off")

    plt.suptitle("Double x-winding cancels to trivial; x+y combined stays obstructed "
                 "(the exact contradiction count is not meaningful)", y=1.02, fontsize=12.5)
    plt.tight_layout()
    fig.savefig("winding_double_and_xy.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_double_and_xy.png")


if __name__ == "__main__":
    main()

"""Figure contrasting Case0 (contractible hexagonal loop, E=0 at every D3,
0 state-3 sites ever needed) against the non-contractible loop that winds
once around the torus (0 local defects too, but needs 38 state-3 sites at
D3=0 and pays a real, growing cost for D3>0) -- see winding_closed_loop.py
for the construction and closed_loop_demo.py for Case0's own numbers."""
import numpy as np
import matplotlib.pyplot as plt

import pickle

from rhombile_lattice import (
    RhombileLattice, build_rhombi, build_triangle_hop_graph, route_string,
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
    full_state_via_mincut, total_energy, simulated_annealing, conflict_graph_bipartition,
)
from dual_string_demo import draw_lattice
from winding_closed_loop import enumerate_triangles_all, NX, NY, Y_ROW, N_WAYPOINTS


def build_loop():
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1 = enumerate_triangles_all(lat)
    rhombi, _adj = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)
    ns = np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int)
    pts = [lat._position(int(n) % NX, Y_ROW, "r1") for n in ns]

    touched = []
    required_start = None
    first_ts = None
    for i in range(N_WAYPOINTS):
        p0, p1 = pts[i], pts[(i + 1) % N_WAYPOINTS]
        is_closing = (i == N_WAYPOINTS - 1)
        if not is_closing:
            if required_start is None:
                ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
            else:
                n_end = int(ns[i + 1]) % NX
                te_candidates = tri_by_r1[lat._site_index(n_end, Y_ROW, "r1")]
                best = None
                for te in te_candidates:
                    try:
                        nodes_c, bonds_c = route_string(sibling, hop, required_start, te)
                    except ValueError:
                        continue
                    if best is None or len(nodes_c) < len(best[1]):
                        best = (te, nodes_c, bonds_c)
                te, nodes, bonds = best
                ts = required_start
        else:
            te_target = first_ts
            nodes, bonds = route_string(sibling, hop, required_start, te_target)
            ts, te = required_start, te_target
        flipped = apply_dual_string_defect(triangles, nodes, bonds)
        touched.extend(flipped)
        required_start = te
        if i == 0:
            first_ts = ts
    return lat, touched, pts


def main():
    lat, touched, pts = build_loop()
    on = [b for b in touched if b["J"] != 0.0]
    off = [b for b in touched if b["J"] == 0.0]
    ft = frustrated_triangles(lat)

    D3_grid = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    # already computed once (winding_closed_loop_figure.py's first run,
    # committed to git history) -- reused here verbatim so the figure's
    # text panel doesn't drift from a fresh re-anneal's sampling noise
    e_loop_mincut = [0.0, 9.5, 19.0, 28.5, 38.0, 54.0, 70.0, 72.0, 72.0]
    e_loop_sa = [0.0, 7.75, 15.5, 25.25, 31.0, 44.5, 56.0, 56.0, 56.0]

    # actual D3=1 ground state (for coloring the lattice panel by state,
    # not species) -- re-derive the SA-best states since the frozen list
    # above only kept the energies
    D3_plot = 1.0
    seed_states, _ = full_state_via_mincut(lat, D3_plot)
    best_e, best_states = total_energy(lat, seed_states, D=(0.0, 0.0, D3_plot)), seed_states
    rng = np.random.default_rng(7)
    for _ in range(5):
        st, en = simulated_annealing(lat, (0.0, 0.0, D3_plot), rng, n_sweeps=2500,
                                      T_start=0.5, T_end=1e-6, states=seed_states, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    for _ in range(4):
        st, en = simulated_annealing(lat, (0.0, 0.0, D3_plot), rng, n_sweeps=3000,
                                      T_start=5.0, T_end=1e-6, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    print(f"D3={D3_plot} plotted ground state: E={best_e}")

    e_case0 = [0.0 for _ in D3_grid]  # closed_loop_demo.py's established result: E=0 at every D3

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    ax_lat, ax_curve, ax_text, ax_blank = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    draw_lattice(ax_lat, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                 box=(NX, NY), targets=pts, off_lw=0.5, states=best_states,
                 title=f"Closed loop winding ONCE around the torus (x-direction), D3={D3_plot}\n"
                       f"colored by actual Potts state -- E={best_e:.0f}, 0 local defects")
    handles, labels = ax_lat.get_legend_handles_labels()
    label_map = {"state 0": "state1", "state 1": "state2", "state 2": "state3"}
    ax_lat.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=8)

    ax_curve.plot(D3_grid, e_case0, "s-", color="tab:green",
                  label="Case0: contractible hexagon (0 state-3 sites ever)")
    ax_curve.plot(D3_grid, e_loop_sa, "o-", color="tab:red",
                  label="this loop: non-contractible (SA-best, 38 state-3 sites at D3=0)")
    ax_curve.plot(D3_grid, e_loop_mincut, "x--", color="tab:orange", alpha=0.6,
                  label="this loop: naive hub-fixed/mincut value")
    ax_curve.set_xlabel("D3")
    ax_curve.set_ylabel("ground-state energy")
    ax_curve.set_title("Zero local defects in both cases -- but only the\n"
                        "non-contractible loop costs energy once D3>0")
    ax_curve.legend(fontsize=9)

    ax_text.axis("off")
    text = (
        "Both loops have EXACTLY 0 frustrated triangles (no particle-like topological\n"
        "defects anywhere) -- verified directly for both.\n\n"
        "Case0 (contractible hexagon): separates the torus into an inside and an outside,\n"
        "so a single global 2-coloring (pure state1/state2, 0 state-3 sites) exists and\n"
        "gives E=0 at EVERY D3, even D3->infinity (established earlier in this project).\n\n"
        "This loop (winds once around x, non-contractible): cutting a torus along one of\n"
        "its 2 fundamental generator loops leaves it CONNECTED (a cylinder), not split into\n"
        "two pieces -- so no consistent global 2-coloring exists (72 whole-graph 2-coloring\n"
        "contradictions, even though the diagonal-only conflict graph is still bipartite).\n"
        "At D3=0 this is invisible (state-3 is free, so mincut finds E=0 using 38 state-3\n"
        "sites to route around the obstruction) -- but for ANY D3>0 those 38 sites cost\n"
        "real energy, and SA confirms a genuine, unavoidable floor (E=31 at D3=1, better\n"
        "than naive hub-fixed's 38 but still far above 0; saturates at E=56 for D3>=2).\n\n"
        "This is a clean, decisive version of the flux-sector question: a purely\n"
        "topological energy cost with NO local defect signature at all -- exactly what a\n"
        "non-contractible cycle should cost, and exactly what the earlier direct-vs-\n"
        "winding straight-line experiment could NOT cleanly demonstrate (those compared\n"
        "two different fixed Hamiltonians, not a clean flux effect)."
    )
    ax_text.text(0.02, 0.98, text, va="top", ha="left", fontsize=9.3, family="monospace",
                 transform=ax_text.transAxes)

    ax_blank.axis("off")

    plt.suptitle("A non-contractible (torus-winding) closed loop costs energy with ZERO local defects",
                 y=1.01, fontsize=13)
    plt.tight_layout()
    fig.savefig("winding_closed_loop_energy.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_closed_loop_energy.png")


if __name__ == "__main__":
    main()

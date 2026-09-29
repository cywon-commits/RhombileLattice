"""Two follow-ups to dual_string_demo / density_demo:

1. Minimum-weight matching, extended from 4 to 6 defects: pairing 6
   target points via the true minimum-weight matching (brute force over
   all 15 pairings, weighted by triangle-hop graph distance) gives a
   clean, predictable result -- non-crossing strings, lengths add up
   exactly. Any other pairing not only costs more, its strings can cross
   and partially cancel, so neither the total length nor even the defect
   count comes out as naively expected.

2. Node (state-2) vs. link (violated-bond) excitations across the D3
   sweep, split out of the same exact mincut solutions used in
   dual_string_demo/density_demo: a single isolated defect pair switches
   from all-node to all-link in one sharp jump at D3=2; a system with
   several defects steps through the switch gradually, one branch point
   at a time.
"""
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, frustrated_triangles,
    route_string_between_points, apply_dual_string_defect, conflict_graph_bipartition,
    state2_and_violations_via_mincut,
)
from dual_string_demo import draw_lattice, build_dual
from density_demo import build_dense_grid

NX, NY = 12, 12


def _all_perfect_matchings(points):
    if not points:
        yield []
        return
    first = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for m in _all_perfect_matchings(rest):
            yield [(first, points[i])] + m


def matching_6defects_demo(pts=((1.0, 1.0), (9.5, 1.5), (2.0, 9.0),
                                 (9.0, 9.5), (5.5, 1.2), (5.0, 9.3))):
    lat0, tri0, rho0, sib0, hop0 = build_dual(NX, NY)
    n = len(pts)
    dist = {}
    for a, b in combinations(range(n), 2):
        _, _, nodes, _ = route_string_between_points(rho0, sib0, hop0, pts[a], pts[b])
        dist[(a, b)] = dist[(b, a)] = len(nodes)

    best_cost, best_m, worst_cost, worst_m = None, None, None, None
    for m in _all_perfect_matchings(list(range(n))):
        cost = sum(dist[(a, b)] for a, b in m)
        if best_cost is None or cost < best_cost:
            best_cost, best_m = cost, m
        if worst_cost is None or cost > worst_cost:
            worst_cost, worst_m = cost, m

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
    for ax, matching, label in zip(axes, [best_m, worst_m],
                                    ["BEST (min-weight) matching", "WORST matching (crossing paths)"]):
        lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)
        on, off = [], []
        for a, b in matching:
            _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, pts[a], pts[b])
            flipped = apply_dual_string_defect(triangles, nodes, bonds)
            on.extend(x for x in flipped if x["J"] != 0.0)
            off.extend(x for x in flipped if x["J"] == 0.0)
        ft = frustrated_triangles(lat)
        edges, _, _, _ = conflict_graph_bipartition(lat)
        draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft, box=(NX, NY),
                     targets=list(pts), title=f"{label}\nmatching={matching}, L={len(edges)}, defects={len(ft)}")
    plt.suptitle("Same 6 target defect points: optimal matching stays clean;\n"
                  "a bad matching's strings cross and partially cancel", y=1.03, fontsize=12.5)
    plt.tight_layout()
    fig.savefig("matching_6defects.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"best matching={best_m} L={best_cost}; worst matching={worst_m} predicted_L={worst_cost}")


def node_vs_link_demo():
    D3_grid = np.linspace(0, 3.6, 300)

    lat1 = RhombileLattice(10, 10)
    apply_string_defect(lat1, np.array([0.5, 2.1]), np.array([7.5, 2.1]), wrap=False)
    edges1, _, _, _ = conflict_graph_bipartition(lat1)
    L1 = len(edges1)
    n2_1, nv_1 = zip(*[state2_and_violations_via_mincut(lat1, d) for d in D3_grid])

    lat2, _ = build_dense_grid(16, 16, 4.0)
    edges2, _, _, _ = conflict_graph_bipartition(lat2)
    L2 = len(edges2)
    n2_2, nv_2 = zip(*[state2_and_violations_via_mincut(lat2, d) for d in D3_grid])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
    panels = [(axes[0], n2_1, nv_1, L1, f"single isolated string (L={L1})\nnode<->link switch is a sharp jump"),
              (axes[1], n2_2, nv_2, L2, f"sparse grid, 16 strings (L={L2})\nswitch is a gradual staircase")]
    for ax, n2, nv, L, title in panels:
        ax.plot(D3_grid, np.array(n2) / L, color="#e8792c", lw=2.3, label="node-type (state-2 fraction)")
        ax.plot(D3_grid, np.array(nv) / L, color="#3b4ba8", lw=2.3, label="link-type (violated-bond fraction)")
        ax.set_xlabel("D3"); ax.set_title(title, fontsize=10.5)
        ax.legend(fontsize=8.5, loc="center right")
    axes[0].set_ylabel("fraction of L")
    plt.suptitle("Node vs. link excitations across the D3 sweep", y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("node_vs_link.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    matching_6defects_demo()
    node_vs_link_demo()

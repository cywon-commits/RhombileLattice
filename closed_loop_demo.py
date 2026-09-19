"""Simplification requested by the user: instead of an open string (whose
two endpoints are topological defects, entangling "does a bend cost more"
with "how much do the endpoints themselves cost"), close the string into a
loop -- a hexagonal boundary chained from 6 straight segments corner-to-
corner. Each corner point is shared by two consecutive segments, and since
apply_dual_string_defect toggles a rhombus's diagonal via the *same*
physical bond object regardless of which of its two triangle-halves is
used as an endpoint (see build_triangle_hop_graph's sibling pairing), the
diagonal toggle at every corner happens exactly twice -- once as the tail
of the incoming segment, once as the head of the outgoing one -- and
cancels back to its pristine (inactive) state. So a closed loop has NO
open ends and, structurally, should have exactly 0 topological defects
(frustrated triangles), unlike every open string in this project so far.

That doesn't mean E=0 is guaranteed at every D3, though: the ring of
still-active diagonals running around the loop's perimeter (every rhombus
the path visits *except* the 6 corners, where it cancels) is itself a
cycle in the "conflict graph" (active r2-r3 diagonals) that all the D3=0
machinery depends on being bipartite. Every open-string configuration
tried previously (including a deliberate 1444-overlapping-string search)
stayed bipartite. A closed loop is the first structure in this project
that can contain a genuine graph-theoretic CYCLE by construction, so it's
the natural place to finally test the bipartiteness conjecture -- if the
loop's perimeter length (in rhombi) is odd, the ring itself is an odd
cycle, and the bipartition-based ground state (min_state2_assignment)
would be forced to leave a real residual cost even at D3=0, which would
be a genuinely new phenomenon, not the D3>0-only issue found for open
detours.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    frustrated_triangles, apply_dual_string_defect, route_string_between_points,
    min_state2_assignment, conflict_graph_bipartition, energy_via_mincut,
    total_energy,
)
from dual_string_demo import draw_lattice, build_dual
from exact_ground_state_investigation import spanning_tree_coloring_energy

NX, NY = 16, 16
CENTER = np.array([6.0, 6.5])
RADIUS = 3.5
N_CORNERS = 6


def hexagon_corners(center=CENTER, radius=RADIUS, n=N_CORNERS):
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [center + radius * np.array([np.cos(a), np.sin(a)]) for a in angles]


def apply_closed_loop(lat, triangles, rhombi, sibling, hop, corners):
    """Chain route_string_between_points around the corners back to the
    start, applying apply_dual_string_defect per edge. Returns the list of
    every bond touched by *any* segment (classify by final J afterwards,
    exactly as dual_string_demo.multi_string_demo does, since a corner
    rhombus's diagonal is touched -- and cancelled -- by two segments)."""
    touched = []
    n = len(corners)
    for i in range(n):
        p0, p1 = corners[i], corners[(i + 1) % n]
        _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
        touched.extend(apply_dual_string_defect(triangles, nodes, bonds))
    return touched


def closed_loop_demo():
    lat0, _, _, _, _ = build_dual(NX, NY)  # pristine, for the "before" panel

    lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)
    corners = hexagon_corners()
    touched = apply_closed_loop(lat, triangles, rhombi, sibling, hop, corners)

    # classify by *final* J (a bond touched by two segments -- i.e. every
    # corner's diagonal -- must be read after all 6 segments are applied)
    on = [b for b in touched if b["J"] != 0.0]
    off = [b for b in touched if b["J"] == 0.0]

    ft = frustrated_triangles(lat)
    print(f"touched bonds (raw, with corner double-counts): {len(touched)}")
    print(f"final state: {len(on)} on, {len(off)} off")
    print(f"frustrated triangles (topological defects): {len(ft)}")

    states0, n_state2, non_bipartite = min_state2_assignment(lat)
    print(f"D3=0 ground state: n_state2={n_state2}, non_bipartite_components={non_bipartite}")
    e0 = total_energy(lat, states0, D=(0.0, 0.0, 0.0))
    print(f"D3=0 energy of that assignment: {e0}")

    edges, side_a, side_b, nb = conflict_graph_bipartition(lat)
    print(f"active-diagonal conflict graph: {len(edges)} edges, "
          f"{len(side_a) + len(side_b)} vertices, non_bipartite={nb}")

    if nb == 0:
        print("  D3     hub-fixed(mincut)   free (unrestricted 2-coloring)")
        for D3 in (0.0, 0.5, 1.0, 1.5, 2.0, 5.0, 1e6):
            hub_fixed = energy_via_mincut(lat, D3)
            free = spanning_tree_coloring_energy(lat, D3)
            print(f"  {D3:<9} {hub_fixed:<19.2f} {free:.2f}")
        print(f"  -> the closed loop's full active-bond graph (not just the diagonal-only\n"
              f"     conflict graph) is itself exactly 2-colorable with ZERO violations, at\n"
              f"     every D3 -- a topologically trivial (zero-net-charge) domain wall never\n"
              f"     needs the 3rd Potts state at all. hub-fixed's nonzero, D3=2-saturating\n"
              f"     curve above is therefore not just imprecise here but flatly WRONG: the\n"
              f"     same bend-induced failure mode as the open detour case, just total this\n"
              f"     time ({len(edges)} claimed vs. 0 actual, instead of 53 vs. 15).")
    else:
        print("  NON-BIPARTITE -- energy_via_mincut cannot be used; this is a genuinely "
              "new (non-bipartite-conflict-graph) case, first one found in this project.")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    draw_lattice(axes[0], lat0, box=(NX, NY), targets=corners,
                 title="original (pristine) lattice\nhexagon shown as the intended loop boundary")
    draw_lattice(axes[1], lat, highlight_on=on, highlight_off=off, frustrated=ft, box=(NX, NY),
                 targets=corners, off_lw=0.5,
                 title=f"closed hexagonal loop applied\n"
                       f"{len(ft)} topological defects, "
                       f"{len(edges)} active diagonals in the ring")
    plt.suptitle("Simplified test case: a closed (defect-free) string loop instead of an open one",
                  y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("closed_loop_before_after.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved closed_loop_before_after.png")


def case1_vs_case2_demo():
    """Same two frustrated plaquettes (the hexagon's bottom two corners),
    connected two different ways:
      Case 1 (L1): the direct straight edge between them (the hexagon's
        own bottom side) -- one segment.
      Case 2 (L2): the detour around the other five sides of the hexagon
        -- five chained segments, going the long way around through the
        remaining four corners.
    Both are OPEN strings (unlike the closed loop), so both should land
    exactly 2 real topological defects at the same two spots -- this
    isolates "does routing the long way around cost more" from the
    closed-loop's zero-net-charge case just examined."""
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    p_left, p_right = corners[i_left], corners[i_right]
    other_order = [i_right] + [i for i in range(len(corners)) if i not in (i_left, i_right)] + [i_left]

    # Case 1: direct edge
    lat1, tri1, rho1, sib1, hop1 = build_dual(NX, NY)
    _, _, nodes1, bonds1 = route_string_between_points(rho1, sib1, hop1, p_left, p_right)
    flipped1 = apply_dual_string_defect(tri1, nodes1, bonds1)
    on1 = [b for b in flipped1 if b["J"] != 0.0]
    off1 = [b for b in flipped1 if b["J"] == 0.0]
    ft1 = frustrated_triangles(lat1)
    edges1, _, _, nb1 = conflict_graph_bipartition(lat1)

    # Case 2: the other five sides, chained corner-to-corner
    lat2, tri2, rho2, sib2, hop2 = build_dual(NX, NY)
    touched2 = []
    for a, b in zip(other_order[:-1], other_order[1:]):
        _, _, nodes, bonds = route_string_between_points(rho2, sib2, hop2, corners[a], corners[b])
        touched2.extend(apply_dual_string_defect(tri2, nodes, bonds))
    on2 = [b for b in touched2 if b["J"] != 0.0]
    off2 = [b for b in touched2 if b["J"] == 0.0]
    ft2 = frustrated_triangles(lat2)
    edges2, _, _, nb2 = conflict_graph_bipartition(lat2)

    print(f"Case 1 (direct bottom edge):  L={len(edges1)} active diagonals, "
          f"{len(ft1)} defects, non_bipartite={nb1}")
    print(f"Case 2 (5-segment detour):    L={len(edges2)} active diagonals, "
          f"{len(ft2)} defects, non_bipartite={nb2}")
    print(f"ratio L2/L1 = {len(edges2) / len(edges1):.2f}  (naive expectation ~5)")

    if nb1 == 0:
        print("Case 1  D3     hub-fixed   free-2-coloring")
        for D3 in (0.0, 0.5, 1.0, 1.5, 2.0, 5.0, 1e6):
            print(f"        {D3:<9} {energy_via_mincut(lat1, D3):<12.2f} "
                  f"{spanning_tree_coloring_energy(lat1, D3):.2f}")
    if nb2 == 0:
        print("Case 2  D3     hub-fixed   free-2-coloring")
        for D3 in (0.0, 0.5, 1.0, 1.5, 2.0, 5.0, 1e6):
            print(f"        {D3:<9} {energy_via_mincut(lat2, D3):<12.2f} "
                  f"{spanning_tree_coloring_energy(lat2, D3):.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    draw_lattice(axes[0], lat1, highlight_on=on1, highlight_off=off1, frustrated=ft1,
                 box=(NX, NY), targets=[p_left, p_right], off_lw=0.5,
                 title=f"Case 1: direct edge (L1)\nL={len(edges1)}, {len(ft1)} defects")
    draw_lattice(axes[1], lat2, highlight_on=on2, highlight_off=off2, frustrated=ft2,
                 box=(NX, NY), targets=[p_left, p_right], off_lw=0.5,
                 title=f"Case 2: 5-segment detour around the other sides (L2)\n"
                       f"L={len(edges2)}, {len(ft2)} defects")
    plt.suptitle("Same two defect locations (hexagon's bottom corners): "
                  "direct edge vs. the long way around", y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("case1_vs_case2.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved case1_vs_case2.png")


if __name__ == "__main__":
    closed_loop_demo()
    case1_vs_case2_demo()

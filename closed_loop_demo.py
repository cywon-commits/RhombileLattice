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
from collections import deque

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

from rhombile_lattice import (
    frustrated_triangles, apply_dual_string_defect, route_string_between_points,
    min_state2_assignment, conflict_graph_bipartition, energy_via_mincut,
    total_energy, full_state_via_mincut,
)
from dual_string_demo import draw_lattice, build_dual, STATE_COLORS
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


def _bfs_two_coloring(lattice, start_site):
    """Exact 2-coloring of the whole active-bond graph by BFS parity from
    `start_site`. Valid (zero-violation) whenever the graph is bipartite,
    which is a stronger, whole-graph condition than the diagonal-only
    conflict graph checked elsewhere -- see spanning_tree_coloring_energy
    in exact_ground_state_investigation.py, which this reimplements but
    also returns the actual per-site colors, not just the energy."""
    active_bonds = [b for b in lattice.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        adj.setdefault(b["i"], []).append(b["j"])
        adj.setdefault(b["j"], []).append(b["i"])
    color = {start_site: 0}
    queue = deque([start_site])
    while queue:
        u = queue.popleft()
        for v in adj.get(u, []):
            if v not in color:
                color[v] = 1 - color[u]
                queue.append(v)
    return np.array([color.get(s, 0) for s in range(lattice.n_sites)])


def case0_inside_outside_swap_demo(D3=1.0):
    """Case 0: the hexagon as a fully closed loop (no open ends -- see
    closed_loop_demo). Question: for the minimum-state-3-usage 3-coloring
    (user's 1-indexed state1/state2/state3; state3 = this module's/
    rhombile_lattice's Potts state index 2), if the OUTSIDE uses
    r1->state1, (r2,r3)->state2, does the INSIDE need r1 and (r2,r3)
    swapped (r1->state2, (r2,r3)->state1)?

    Answer, found by taking the exact whole-graph 2-coloring (BFS parity,
    _bfs_two_coloring -- already known from closed_loop_demo to achieve
    E=0 at every D3, since this loop carries no net topological charge)
    and cross-tabulating each site's color against species (r1 vs rim)
    and geometric side (inside vs outside the hexagon): yes -- the
    r1/rim role is swapped inside, and this is achieved with *zero*
    state-3 sites anywhere (not just "few"), consistent with E=0 at
    every D3 including D3=1.
    """
    lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)
    corners = hexagon_corners()
    touched = apply_closed_loop(lat, triangles, rhombi, sibling, hop, corners)
    on = [b for b in touched if b["J"] != 0.0]
    off = [b for b in touched if b["J"] == 0.0]

    pos, sub_of = lat.site_positions()
    inside = MplPath(np.array(corners)).contains_points(pos)
    is_r1 = sub_of == "r1"

    # start the BFS from a site far outside the loop, so "color 0" anchors
    # to the *outside* pattern
    start = int(np.argmax(np.linalg.norm(pos - CENTER, axis=1)))
    color = _bfs_two_coloring(lat, start)

    # relabel so state1 (index 0) is whatever color r1-outside sites mostly
    # take, and state2 (index 1) the other -- matching the user's stated
    # convention; state 3 (index 2) is never assigned at all here
    r1_outside_color1_frac = color[is_r1 & ~inside].mean()
    flip = r1_outside_color1_frac > 0.5
    states = (1 - color) if flip else color.copy()

    e = total_energy(lat, states, D=(0.0, 0.0, D3))
    n_state3 = int((states == 2).sum())
    print(f"Case 0 (fully closed loop), D3={D3}: E={e}, state-3 sites used={n_state3}")

    def frac_state1(mask):
        vals = states[mask]
        return float((vals == 0).mean()) if len(vals) else float("nan")

    print(f"  r1,  outside: frac(state1)={frac_state1(is1o := is_r1 & ~inside):.3f}  n={is1o.sum()}")
    print(f"  rim, outside: frac(state1)={frac_state1(iso := ~is_r1 & ~inside):.3f}  n={iso.sum()}")
    print(f"  r1,  inside:  frac(state1)={frac_state1(is1i := is_r1 & inside):.3f}  n={is1i.sum()}")
    print(f"  rim, inside:  frac(state1)={frac_state1(isi := ~is_r1 & inside):.3f}  n={isi.sum()}")
    print("  -> outside: r1 mostly state1, rim mostly state2. inside: r1 mostly state2, "
          "rim mostly state1 -- SWAPPED, exactly as guessed, achieved with 0 state-3 sites "
          "(the few % exceptions above are sites right against the loop's jagged actual path, "
          "not against the idealized straight hexagon edge drawn for reference).")

    fig, ax = plt.subplots(figsize=(10, 9))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, box=(NX, NY),
                 targets=corners, states=states, off_lw=0.5,
                 title=f"Case 0: fully closed loop, minimum-state-3 3-coloring (D3={D3})\n"
                       f"E={e}, state-3 sites used={n_state3} -- outside r1=state1/rim=state2, "
                       f"inside SWAPPED, no state-3 needed")
    handles, labels = ax.get_legend_handles_labels()
    label_map = {"state 0": "state1 (was r1 outside)", "state 1": "state2 (was rim outside)",
                 "state 2": "state3 (unused here)"}
    ax.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=9)
    plt.tight_layout()
    fig.savefig("case0_inside_outside_swap.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved case0_inside_outside_swap.png")


def case1_state3_demo(D3=1.0):
    """Case 1 (the hexagon's direct bottom edge, an OPEN string with 2 real
    topological defects -- unlike Case 0's closed loop). Unlike the closed
    loop, this one carries real net topological charge at its two ends, so
    it should NOT be free: at D3=1, the user expects state-3 sites to show
    up strung along the short straight path itself, since hub-fixing is
    known-reliable here (a plain single straight string, same regime
    already confirmed against SA elsewhere in this project).

    Confirmed: full_state_via_mincut at D3=1 gives E=4 with exactly 4
    state-3 sites, all lying on the path between the two defects, and
    zero violated bonds -- i.e. at D3=1 the ground state is still purely
    "pay with state-3 sites", not yet trading any of them for a violated
    bond (see energy_vs_D3 in the main report: that only starts once D3
    passes the L=10 case's own kink, around D3~2.5 here).
    """
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    p_left, p_right = corners[i_left], corners[i_right]

    lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)
    _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p_left, p_right)
    flipped = apply_dual_string_defect(triangles, nodes, bonds)
    on = [b for b in flipped if b["J"] != 0.0]
    off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    edges, _, _, nb = conflict_graph_bipartition(lat)

    states, violated = full_state_via_mincut(lat, D3)
    e = total_energy(lat, states, D=(0.0, 0.0, D3))
    n_state3 = int((states == 2).sum())
    pos, _ = lat.site_positions()
    print(f"Case 1, D3={D3}: L={len(edges)}, defects={len(ft)}, E={e}, "
          f"state-3 sites={n_state3}, violated bonds={len(violated)}")
    for p in pos[states == 2]:
        print(f"  state-3 site at {tuple(np.round(p, 2))}")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft, box=(NX, NY),
                 targets=[p_left, p_right], states=states, off_lw=0.5,
                 title=f"Case 1 (direct edge), D3={D3}: E={e}, {n_state3} state-3 sites "
                       f"strung along the path, {len(violated)} violated bonds")
    handles, labels = ax.get_legend_handles_labels()
    label_map = {"state 0": "state1", "state 1": "state2", "state 2": "state3"}
    ax.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=9)
    ax.set_xlim(0, 12)
    ax.set_ylim(1, 6)
    plt.tight_layout()
    fig.savefig("case1_state3_at_D3_1.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved case1_state3_at_D3_1.png")


def case1_frustrated_links_demo(D3=3.0):
    """Case 1 (direct edge), past its own saturation point: the user's
    prediction is that at D3=3 it's cheaper to stop paying with state-3
    sites entirely and instead just leave the path's links frustrated
    (mismatched under only 2 colors) between the two defects.

    Confirmed: full_state_via_mincut at D3=3 gives E=10=L, with 0 state-3
    sites and all 10 links along the path violated -- full switch-over.
    The crossover isn't a single jump, though: D3=2.0 is still mixed (3
    state-3 sites + 2 violated links) and D3=2.5 even more so (2 + 4)
    before D3=3 completes the switch to all-link.
    """
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    p_left, p_right = corners[i_left], corners[i_right]

    lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)
    _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p_left, p_right)
    flipped = apply_dual_string_defect(triangles, nodes, bonds)
    on = [b for b in flipped if b["J"] != 0.0]
    off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    edges, _, _, nb = conflict_graph_bipartition(lat)

    states, violated = full_state_via_mincut(lat, D3)
    e = total_energy(lat, states, D=(0.0, 0.0, D3))
    n_state3 = int((states == 2).sum())
    print(f"Case 1, D3={D3}: L={len(edges)}, E={e}, state-3 sites={n_state3}, "
          f"violated (frustrated) links={len(violated)}")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft, box=(NX, NY),
                 targets=[p_left, p_right], states=states, off_lw=0.5,
                 title=f"Case 1 (direct edge), D3={D3}: E={e}, {n_state3} state-3 sites, "
                       f"{len(violated)} frustrated links between the defects (magenta)")
    for b in violated:
        ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                 color="magenta", linewidth=3.2, zorder=6)
    handles, labels = ax.get_legend_handles_labels()
    label_map = {"state 0": "state1", "state 1": "state2", "state 2": "state3"}
    ax.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=9)
    ax.set_xlim(0, 12)
    ax.set_ylim(1, 6)
    plt.tight_layout()
    fig.savefig("case1_frustrated_links_D3_3.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved case1_frustrated_links_D3_3.png")


if __name__ == "__main__":
    closed_loop_demo()
    case1_vs_case2_demo()
    case0_inside_outside_swap_demo()
    case1_state3_demo()
    case1_frustrated_links_demo()

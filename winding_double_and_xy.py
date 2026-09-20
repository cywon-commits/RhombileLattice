"""Two cheap follow-up tests to winding_closed_loop.py, prompted by the
user's own predictions:

  (1) Winding the SAME direction (x) TWICE: predicted to be TRIVIAL
      (resolvable to E=0 at every D3), since the whole-graph 2-coloring
      consistency check is a Z2 (mod-2) invariant -- crossing the domain
      wall's own "seam" an even number of times flips the antiphase role
      an even number of times, returning to consistency. Tested here as
      two INDEPENDENT single-x-wind loops (each individually identical in
      construction to winding_closed_loop.py's own loop) placed on two
      different rows of the SAME lattice -- if the invariant is additive
      mod 2 over independent cycles, their combination should cancel.

  (2) Winding x ONCE and y ONCE simultaneously: predicted to be a
      DIFFERENT, still-nontrivial sector (Z2 x Z2 has 4 elements: (0,0)
      trivial, (1,0), (0,1), (1,1) -- all distinct), NOT reducible to
      trivial by "cutting a crossing", contrary to the initial guess that
      it must resolve like two independent trivial loops. Tested here as
      one x-winding loop (row m=2) plus one y-winding loop (column n=15)
      placed on the SAME lattice.

Both reuse the exact corner-cancelling waypoint-chaining method validated
in winding_closed_loop.py (each independent loop verified to have 0 local
defects on its own), extended to also route in the y (m) direction using
the same unfiltered triangle enumeration.
"""
import numpy as np
from collections import deque

from rhombile_lattice import (
    RhombileLattice, build_rhombi, build_triangle_hop_graph, route_string,
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
    full_state_via_mincut, total_energy, simulated_annealing, conflict_graph_bipartition,
)
from winding_closed_loop import enumerate_triangles_all

NX, NY = 30, 8
N_WAYPOINTS = 6


def add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, waypoints):
    """Chains route_string/route_string_between_points around `waypoints`
    (a list of (n, m) unit-cell coordinates, already reduced mod NX/NY)
    back to the first one, pinning each segment's start half to the
    previous segment's end half so corners cancel exactly as in Case0 /
    winding_closed_loop.py. Returns the list of touched bonds."""
    pts = [lat._position(n, m, "r1") for n, m in waypoints]
    n = len(waypoints)
    touched = []
    required_start = None
    first_ts = None
    for i in range(n):
        p0, p1 = pts[i], pts[(i + 1) % n]
        is_closing = (i == n - 1)
        if not is_closing:
            if required_start is None:
                ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
            else:
                n_end, m_end = waypoints[i + 1]
                te_candidates = tri_by_r1[lat._site_index(n_end, m_end, "r1")]
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
    return touched


def whole_graph_contradictions(lat):
    active_bonds = [b for b in lat.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        adj.setdefault(b["i"], []).append(b["j"])
        adj.setdefault(b["j"], []).append(b["i"])
    color = {}
    contradiction_pairs = set()
    for start_site in adj:
        if start_site in color:
            continue
        color[start_site] = 0
        q = deque([start_site])
        while q:
            u = q.popleft()
            for v in adj.get(u, []):
                if v not in color:
                    color[v] = 1 - color[u]
                    q.append(v)
                elif color[v] == color[u]:
                    contradiction_pairs.add((min(u, v), max(u, v)))
    return len(contradiction_pairs)


def report(lat, label):
    ft = frustrated_triangles(lat)
    edges, side_a, side_b, nb = conflict_graph_bipartition(lat)
    contras = whole_graph_contradictions(lat)
    print(f"[{label}] local defects={len(ft)}, conflict-graph non_bipartite={nb}, "
          f"whole-graph contradictions={contras}")
    for D3 in (0.0, 1.0, 2.0):
        try:
            states, violated = full_state_via_mincut(lat, D3)
            e = total_energy(lat, states, D=(0.0, 0.0, D3))
            print(f"    D3={D3}: mincut E={e}, n3={(states==2).sum()}, n_violated={len(violated)}")
        except ValueError as exc:
            print(f"    D3={D3}: mincut failed ({exc})")
    return ft, contras


def build_base(lat):
    triangles, tri_by_r1 = enumerate_triangles_all(lat)
    rhombi, _adj = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)
    return triangles, tri_by_r1, rhombi, sibling, hop


def test_double_x_wind():
    print("\n=== TEST 1: two independent x-winding loops (rows m=2 and m=4) ===")
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1, rhombi, sibling, hop = build_base(lat)
    ns = list(np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int) % NX)

    wp1 = [(int(n), 2) for n in ns]
    add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, wp1)
    report(lat, "after loop 1 alone (row m=2)")

    wp2 = [(int(n), 4) for n in ns]
    add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, wp2)
    ft, contras = report(lat, "after BOTH loops (rows m=2 and m=4)")
    return ft, contras


def test_xy_combined():
    print("\n=== TEST 2: one x-winding loop (row m=2) + one y-winding loop (column n=15) ===")
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1, rhombi, sibling, hop = build_base(lat)
    ns = list(np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int) % NX)
    wp_x = [(int(n), 2) for n in ns]
    add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, wp_x)
    report(lat, "after x-winding loop alone")

    ms = list(np.linspace(0, NY, N_WAYPOINTS, endpoint=False).round().astype(int) % NY)
    wp_y = [(15, int(m)) for m in ms]
    add_winding_loop(lat, triangles, tri_by_r1, rhombi, sibling, hop, wp_y)
    ft, contras = report(lat, "after BOTH x-winding and y-winding loops")
    return ft, contras


if __name__ == "__main__":
    ft1, c1 = test_double_x_wind()
    ft2, c2 = test_xy_combined()
    print("\n=== SUMMARY ===")
    print(f"double x-wind:  local defects={len(ft1)}, whole-graph contradictions={c1} "
          f"-> {'TRIVIAL (matches Z2 prediction)' if c1 == 0 else 'STILL OBSTRUCTED (prediction wrong)'}")
    print(f"x + y combined: local defects={len(ft2)}, whole-graph contradictions={c2} "
          f"-> {'TRIVIAL' if c2 == 0 else 'STILL OBSTRUCTED (matches Z2xZ2 prediction: (1,1) != (0,0))'}")

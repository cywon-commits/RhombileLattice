"""Applies height_function.py's Burgers-vector diagnostic to the 4 cases
that matter for the user's question:
  1. pristine lattice -- baseline (already checked in height_function.py
     itself: perfectly consistent, zero holonomy everywhere).
  2. Case0: a closed, CONTRACTIBLE hexagonal loop with 0 local defects
     (known E=0 at every D3) -- prediction: still zero holonomy
     everywhere (a step/terrace loop, no dislocation).
  3. the non-contractible loop that winds once around the torus (0 local
     defects, but known nonzero energy at D3>0 from a whole-graph
     2-coloring obstruction invisible to local defect counting) --
     prediction: some cycle now shows NONZERO holonomy -- the first
     direct, quantitative Burgers-vector reading of that obstruction.
  4. a plain open string between 2 real point defects -- prediction: the
     2 endpoints are a genuine dislocation DIPOLE, i.e. carry exactly
     opposite charge (+3 / -3), consistent with the matching principle
     (any path between the same 2 endpoints is a different choice of
     Volterra cut but must leave the same 2 core charges).
"""
from collections import deque

import numpy as np

from rhombile_lattice import (
    RhombileLattice, build_rhombi, build_triangle_hop_graph, route_string,
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
)
from dual_string_demo import build_dual
from closed_loop_demo import hexagon_corners, apply_closed_loop
from height_function import (
    enumerate_triangles_all, build_increments, compute_height, triangle_charge,
)


def report(label, lat):
    triangles = enumerate_triangles_all(lat)
    inc, build_mismatches = build_increments(lat, triangles)
    h, unreached, checks = compute_height(lat, inc)
    nonzero = [(i, j, v) for i, j, v in checks if abs(v) > 1e-9]
    ft = frustrated_triangles(lat)
    charges = [round(triangle_charge(t), 3) for t in triangles if abs(triangle_charge(t)) > 1e-9]
    print(f"\n=== {label} ===")
    print(f"  local (frustrated-triangle) defects: {len(ft)}, per-defect charges: {charges}")
    print(f"  build-consistency mismatches: {len(build_mismatches)} (should be 0)")
    print(f"  unreached sites: {len(unreached)} (should be 0)")
    print(f"  non-tree bonds checked: {len(checks)}, nonzero holonomy: {len(nonzero)}")
    if nonzero:
        vals = sorted(set(round(v, 3) for _, _, v in nonzero))
        print(f"  DISTINCT nonzero holonomy values: {vals}")
    return h, checks, ft, charges


def case0_hexagon():
    lat, triangles, rhombi, sibling, hop = build_dual(16, 16)
    corners = hexagon_corners()
    apply_closed_loop(lat, triangles, rhombi, sibling, hop, corners)
    return lat


def winding_loop(nx=30, ny=8, y_row=2, n_waypoints=6):
    lat = RhombileLattice(nx, ny)
    triangles, tri_by_r1 = enumerate_triangles_all_with_index(lat)
    rhombi, _adj = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)

    ns = np.linspace(0, nx, n_waypoints, endpoint=False).round().astype(int)
    pts = [lat._position(int(n) % nx, y_row, "r1") for n in ns]

    required_start = None
    first_ts = None
    for i in range(n_waypoints):
        p0, p1 = pts[i], pts[(i + 1) % n_waypoints]
        is_closing = (i == n_waypoints - 1)
        if not is_closing:
            if required_start is None:
                ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
            else:
                n_end = int(ns[i + 1]) % nx
                te_candidates = tri_by_r1[lat._site_index(n_end, y_row, "r1")]
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
            n_start = int(ns[i]) % nx
            te_target = first_ts
            ts_candidates = (tri_by_r1[lat._site_index(n_start, y_row, "r1")]
                             if required_start is None else [required_start])
            best = None
            for ts in ts_candidates:
                try:
                    nodes_c, bonds_c = route_string(sibling, hop, ts, te_target)
                except ValueError:
                    continue
                if best is None or len(nodes_c) < len(best[1]):
                    best = (ts, nodes_c, bonds_c)
            ts, nodes, bonds = best
            te = te_target
        apply_dual_string_defect(triangles, nodes, bonds)
        required_start = te
        if i == 0:
            first_ts = ts
    return lat


def enumerate_triangles_all_with_index(lattice):
    """Same as height_function.enumerate_triangles_all, but also returns
    tri_by_r1 (needed for winding_loop's waypoint pinning, same as
    winding_closed_loop.py's own version)."""
    bond_by_pair = {frozenset((b["i"], b["j"])): b for b in lattice.bonds}
    triangles = []
    tri_by_r1 = {}
    for m in range(lattice.ny):
        for n in range(lattice.nx):
            r1_idx = lattice._site_index(n, m, "r1")
            r1_pos = lattice._position(n, m, "r1")
            nbrs = []
            for other_idx, b in lattice.neighbor_table()[r1_idx]:
                other_pos = b["p2"] if b["i"] == r1_idx else b["p1"]
                ang = np.arctan2(*(other_pos - r1_pos)[::-1])
                nbrs.append((ang, other_idx, other_pos, b))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, pos1, b1 = nbrs[i]
                _, idx2, pos2, b2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None or closing["type"] != "r2r3":
                    continue
                edges = {b1["type"]: b1, b2["type"]: b2}
                if "r1r2" not in edges or "r1r3" not in edges:
                    continue
                ti = len(triangles)
                centroid = (r1_pos + pos1 + pos2) / 3
                triangles.append({"sites": (r1_idx, idx1, idx2), "r1r2": edges["r1r2"],
                                   "r1r3": edges["r1r3"], "r2r3": closing, "centroid": centroid})
                tri_by_r1.setdefault(r1_idx, []).append(ti)
    return triangles, tri_by_r1


def open_string_pair():
    lat, triangles, rhombi, sibling, hop = build_dual(16, 16)
    p_left = np.array([4.0, 5.0])
    p_right = np.array([10.0, 5.0])
    _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p_left, p_right)
    apply_dual_string_defect(triangles, nodes, bonds)
    return lat


if __name__ == "__main__":
    print("### 1. pristine baseline ###")
    lat1 = RhombileLattice(16, 16)
    report("pristine 16x16", lat1)

    print("\n### 2. Case0: closed contractible hexagon loop ###")
    lat2 = case0_hexagon()
    report("Case0 hexagon loop", lat2)

    print("\n### 3. non-contractible winding loop ###")
    lat3 = winding_loop()
    report("winding loop (30x8, winds once in x)", lat3)

    print("\n### 4. open string between 2 real point defects ###")
    lat4 = open_string_pair()
    h4, checks4, ft4, charges4 = report("open string, 2 endpoints", lat4)
    print(f"  endpoint charges sum to: {sum(charges4):.3f} (expect 0: opposite-sign dipole)")

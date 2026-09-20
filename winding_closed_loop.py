"""Directly tests the user's sharper question: the two "winding" tests so
far (straight lines) aren't really testable for winding in a meaningful
sense (an open string's homotopy class relative to its 2 endpoints isn't
gauge-invariant -- see the conversation). Case0's closed hexagonal loop
*was* well-defined and reached E=0 at every D3 because a CONTRACTIBLE loop
separates the torus into two disjoint regions (inside/outside), giving a
consistent global 2-coloring (role swap) once you cross it.

A loop that winds once around the torus (a non-contractible "meridian"),
even with zero local (triangle) defects by construction, does NOT
separate the torus into two disjoint pieces -- cutting a torus along one
of its two fundamental generator loops leaves it connected (a cylinder),
unlike cutting along a contractible loop (which leaves two disjoint
pieces). So the whole-graph 2-coloring trick that made Case0's E=0 result
work may simply not be available here, even at 0 local defects -- this
script builds that loop and checks directly.

Requires routing THROUGH the periodic seam, which enumerate_triangles
structurally excludes (see flux_winding_audit.py) -- so this reimplements
an unfiltered version of it. build_rhombi/build_triangle_hop_graph need no
changes: both use only id(bond) matching, never centroids, so they are
already correct for wrap-touching triangles; only position-based lookups
(nearest_rhombus/nearest_triangle) are unreliable there, so this uses
site-index-based lookups instead of point-based ones wherever the route
must actually cross the seam.
"""
import numpy as np
from collections import deque

from rhombile_lattice import (
    RhombileLattice, build_rhombi, build_triangle_hop_graph, route_string,
    apply_dual_string_defect, frustrated_triangles, conflict_graph_bipartition,
    total_energy, nearest_rhombus, route_string_between_points,
)
from dual_string_demo import draw_lattice

NX, NY = 30, 8
Y_ROW = 2  # the m-row (unit-cell index) the loop runs along
N_WAYPOINTS = 6  # matches Case0's hexagon segment count


def enumerate_triangles_all(lattice):
    """enumerate_triangles without the wrap-exclusion filter -- see that
    function's own docstring for why centroids of wrap-touching triangles
    can be positioned inconsistently; this script never trusts centroids
    for triangles it knows touch the seam."""
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
                centroid = (r1_pos + pos1 + pos2) / 3
                ti = len(triangles)
                triangles.append({
                    "sites": (r1_idx, idx1, idx2), "r1r2": edges["r1r2"],
                    "r1r3": edges["r1r3"], "r2r3": closing, "centroid": centroid,
                    "n": n, "m": m,
                })
                tri_by_r1.setdefault(r1_idx, []).append(ti)
    return triangles, tri_by_r1


def main():
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1 = enumerate_triangles_all(lat)
    rhombi, _adj = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)
    print(f"{len(triangles)} triangles total (unfiltered), {len(rhombi)} rhombi")

    # waypoints: N_WAYPOINTS points evenly spaced in n around the SAME
    # m-row, closing the last one back to the first through the seam
    ns = np.linspace(0, NX, N_WAYPOINTS, endpoint=False).round().astype(int)
    pts = [lat._position(int(n) % NX, Y_ROW, "r1") for n in ns]
    print("waypoint n-values:", list(ns))

    # segment i's END triangle-half must be EXACTLY segment i+1's START
    # triangle-half for the shared waypoint to cancel (same physical half,
    # not just "some triangle at that site") -- so track and reuse the
    # exact index returned by route_string_between_points instead of
    # letting each segment re-pick independently (which is what
    # apply_closed_loop/Case0 does implicitly for its non-wrapping
    # segments, fine there since none of them needed a fallback lookup).
    touched = []
    required_start = None
    for i in range(N_WAYPOINTS):
        p0, p1 = pts[i], pts[(i + 1) % N_WAYPOINTS]
        is_closing = (i == N_WAYPOINTS - 1)
        if not is_closing:
            if required_start is None:
                ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
            else:
                # start half is pinned; still search both possible end halves
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
            # closing segment crosses the seam: end half is pinned to
            # match segment 0's own start half exactly (captured below)
            n_start = int(ns[i]) % NX
            te_target = first_ts
            ts_candidates = (tri_by_r1[lat._site_index(n_start, Y_ROW, "r1")]
                             if required_start is None else [required_start])
            best = None
            for ts in ts_candidates:
                try:
                    nodes_c, bonds_c = route_string(sibling, hop, ts, te_target)
                except ValueError:
                    continue
                if best is None or len(nodes_c) < len(best[1]):
                    best = (ts, nodes_c, bonds_c)
            if best is None:
                raise RuntimeError("no closing route found across the seam")
            ts, nodes, bonds = best
            te = te_target
        flipped = apply_dual_string_defect(triangles, nodes, bonds)
        touched.extend(flipped)
        required_start = te
        if i == 0:
            first_ts = ts
        ft_now = frustrated_triangles(lat)
        print(f"segment {i} ({'CLOSING/seam' if is_closing else 'normal'}): "
              f"ts={ts}, te={te}, {len(nodes)} triangle-halves, {len(bonds)} outer bonds crossed, "
              f"running local defects={len(ft_now)}")

    on = [b for b in touched if b["J"] != 0.0]
    off = [b for b in touched if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    print(f"\nFINAL local (triangle) defects: {len(ft)}")

    edges, side_a, side_b, nb = conflict_graph_bipartition(lat)
    print(f"conflict graph (active diagonals): {len(edges)} edges, "
          f"{len(side_a) + len(side_b)} vertices, non_bipartite_components={nb}")

    # whole active-bond-graph consistency check (the same test that
    # verified Case0's E=0 antiphase swap), watching for contradictions
    active_bonds = [b for b in lat.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        adj.setdefault(b["i"], []).append(b["j"])
        adj.setdefault(b["j"], []).append(b["i"])
    color = {}
    contradiction_sites = set()
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
                    contradiction_sites.add((min(u, v), max(u, v)))
    print(f"whole-graph 2-coloring contradictions (distinct bond pairs): {len(contradiction_sites)}")

    if len(ft) == 0 and nb == 0 and len(contradiction_sites) == 0:
        print("\n-> a valid global E=0 assignment exists at every D3, just like Case0.")
    elif len(ft) == 0:
        print("\n-> ZERO local triangle defects, but the whole-graph 2-coloring still fails: "
              "a genuine flux obstruction with no local (triangle-level) signature at all -- "
              "this is exactly the non-contractible-loop effect the user predicted.")
    else:
        print(f"\n-> construction did not fully cancel: {len(ft)} residual local defects remain "
              f"(construction needs adjusting).")

    fig_title = (f"Closed loop winding once around x (n={list(ns)}), row m={Y_ROW}\n"
                 f"{len(ft)} local defects, non_bipartite_components={nb}, "
                 f"{len(contradiction_sites)} 2-coloring contradictions")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 6))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                 box=(NX, NY), targets=pts, off_lw=0.5, title=fig_title)
    plt.tight_layout()
    fig.savefig("winding_closed_loop.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_closed_loop.png")


if __name__ == "__main__":
    main()

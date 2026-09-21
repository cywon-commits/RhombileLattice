"""Quantifies the "phason"/slope coordinate the user asked about: for each
non-defect triangle, exactly one of its 3 edges is inactive (matched) --
label it by TYPE (r1r2 / r1r3 / r2r3). The one-hot vector of this label is
already a perfectly well-defined POINTWISE (per-triangle) discrete field,
no spatial averaging required -- averaging it over a region gives the
region's (p12, p13, p23) "slope"/phason coordinate, but a string defect
shows up directly, unaveraged, as a literal LINE of anomalous one-hot
labels against the uniform (0,0,1) background, exactly as predicted.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import RhombileLattice, route_string_between_points, apply_dual_string_defect
from dual_string_demo import build_dual
from height_function import enumerate_triangles_all

TYPE_NAMES = ["r1r2", "r1r3", "r2r3"]
TYPE_COLORS = {"r1r2": [1.0, 0.15, 0.15], "r1r3": [0.15, 0.75, 0.15], "r2r3": [0.15, 0.35, 1.0]}


def triangle_inactive_type(tri):
    """Index (0=r1r2, 1=r1r3, 2=r2r3) of the ONE inactive edge, or None if
    the triangle is a defect (0 or >=2 inactive edges -- frustrated
    all-active, or the "over-matched" all-inactive anomaly found for
    Case0's corners)."""
    edges = [tri["r1r2"], tri["r1r3"], tri["r2r3"]]
    inactive = [i for i, e in enumerate(edges) if e["J"] == 0.0]
    return inactive[0] if len(inactive) == 1 else None


def region_slope(triangles):
    """(p_r1r2, p_r1r3, p_r2r3) averaged over every non-defect triangle
    given -- the region-averaged phason/slope coordinate."""
    counts = np.zeros(3)
    n = 0
    for t in triangles:
        k = triangle_inactive_type(t)
        if k is not None:
            counts[k] += 1
            n += 1
    return counts / n if n else counts, n


def pointwise_along_path(triangles, node_indices, label=""):
    """No averaging at all: the literal per-triangle type sequence along a
    string's own routed path, directly showing it as a discrete line of
    anomalous labels against the (assumed) r2r3 background."""
    seq = [TYPE_NAMES[triangle_inactive_type(triangles[t])]
           if triangle_inactive_type(triangles[t]) is not None else "DEFECT"
           for t in node_indices]
    print(f"{label} path type-sequence ({len(seq)} triangles): {seq}")
    return seq


def main():
    print("=== pristine: region slope should be exactly (0,0,1) ===")
    lat0 = RhombileLattice(20, 20)
    tri0 = enumerate_triangles_all(lat0)
    slope0, n0 = region_slope(tri0)
    print(f"  slope (p_r1r2,p_r1r3,p_r2r3) = {np.round(slope0,4)}  over {n0} triangles")

    print("\n=== straight open string: pointwise type sequence along its own path ===")
    lat, triangles, rhombi, sibling, hop = build_dual(20, 20)
    p_left = np.array([4.0, 8.0])
    p_right = np.array([14.0, 8.0])
    ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p_left, p_right)
    apply_dual_string_defect(triangles, nodes, bonds)
    # nodes are indices into build_dual's OWN (filtered) triangle list, and the
    # bonds it mutated are the same physical objects, so read the type straight
    # off that list -- no re-enumeration, no index realignment needed here.
    pointwise_along_path(triangles, nodes, label="open string")
    # for the WHOLE-lattice slope, a fresh unfiltered enumeration (different,
    # unrelated indexing) is fine since we only aggregate over all of it.
    tri_full = enumerate_triangles_all(lat)
    slope1, n1 = region_slope(tri_full)
    print(f"  WHOLE-LATTICE slope after inserting the string = {np.round(slope1,4)} over {n1} "
          f"triangles (barely shifted from background -- the string is a 1D line, "
          f"a vanishing fraction of the whole 2D area, so it barely moves a bulk average)")

    print("\n=== visualization: pristine vs string, colored by inactive-edge type ===")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, lat_, tri_, title in [(axes[0], lat0, tri0, "pristine (uniform r2r3 background)"),
                                    (axes[1], lat, tri_full, "one straight open string")]:
        for t in tri_:
            k = triangle_inactive_type(t)
            color = TYPE_COLORS[TYPE_NAMES[k]] if k is not None else "black"
            r1, a, b = t["sites"]
            pos_r1 = lat_._position(*_unindex(lat_, r1))
            # centroid via the 3 sites' actual positions (recompute directly, since
            # enumerate_triangles_all's stored "centroid" key was stripped in
            # height_function.py's version)
            c = _triangle_centroid(lat_, t)
            ax.scatter(*c, s=6, color=color)
        ax.set_aspect("equal")
        ax.set_title(title)
    plt.tight_layout()
    plt.savefig("phason_slope_map.png", dpi=130)
    print("saved phason_slope_map.png")


def _unindex(lattice, site_idx):
    sub_i = site_idx % 3
    cell = site_idx // 3
    m = cell // lattice.nx
    n = cell % lattice.nx
    subs = ["r1", "r2", "r3"]
    return n, m, subs[sub_i]


def _triangle_centroid(lattice, tri):
    r1, a, b = tri["sites"]
    p_r1 = lattice._position(*_unindex(lattice, r1))
    # for a and b, use the bond's own stored (unwrapped) p1/p2 relative to r1's bond
    def other_pos(edge, anchor):
        return edge["p2"] if edge["i"] == anchor else edge["p1"]
    p_a = other_pos(tri["r1r2"] if {tri["r1r2"]["i"], tri["r1r2"]["j"]} == {r1, a} else tri["r1r3"], r1)
    p_b = other_pos(tri["r1r3"] if {tri["r1r3"]["i"], tri["r1r3"]["j"]} == {r1, b} else tri["r1r2"], r1)
    return (p_r1 + p_a + p_b) / 3


if __name__ == "__main__":
    main()

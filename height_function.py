"""Discrete height-function / Burgers-vector formalism for the rhombile
lattice, testing the user's question: since strings/loops here rearrange
which lattice EDGES carry an active AF bond (not the site positions
themselves), is there a genuine "dislocation" description via a Burgers
vector, analogous to the standard height-function representation of
dimer coverings / lozenge tilings?

Key correspondence (see conversation for the full derivation):
  - Each elementary (r1,r2,r3) triangle always has, away from any defect,
    exactly 1 inactive edge and 2 active edges. This is exactly a perfect
    matching of the triangle-ADJACENCY graph (a honeycomb lattice: each
    triangle has 3 neighbors, one via each edge), where "matched via edge
    e" <=> "edge e is inactive" (a monomer/dimer covering). The default
    (no string) config matches every triangle to its rhombus sibling via
    the diagonal (r2r3) edge; a string defect is a sequence of MONOMER-
    DIMER FLIPS along an augmenting path, exactly the standard operation
    of moving a monomer pair along a string in dimer-model language.
  - This bipartite matching structure gives a standard discrete height
    function h on the SITES (not the triangles): color each triangle by
    its chirality (does its first CCW spoke edge from r1 go to the r2 or
    the r3 corner?); for a triangle with CCW corner order (c0,c1,c2),
    define increment(c_k -> c_{k+1}) = sign(color) * (+1 if that edge is
    active else -2). Any two triangles sharing an edge have OPPOSITE
    color (adjacent triangles in a triangular lattice always alternate
    orientation), which makes this single-valued on every ordinary
    (non-frustrated) triangle (2*(+1) + 1*(-2) = 0 automatically, exactly
    because every triangle has exactly 1 inactive edge). A genuinely
    frustrated triangle (all 3 edges active -- a monomer / point defect)
    breaks this: the sum around it is sign(color)*3, a topological
    "charge" -- exactly a Burgers vector of a single dislocation core.

  - A closed, CONTRACTIBLE domain wall (Case0) with no local defects
    should show h consistent everywhere (zero holonomy on every cycle):
    a closed step/terrace loop, not a dislocation.
  - A non-contractible loop that winds around the torus, even with zero
    LOCAL defects, can still show a nonzero holonomy around the torus's
    own generator cycle -- exactly a screw dislocation with Burgers
    vector along the height axis, threaded through the periodic
    direction. This is the rigorous version of the Z2xZ2 flux sector
    already established (Kasteleyn/Thurston/Kenyon), now phrased in
    explicit Burgers-vector language.
  - Two real point defects connected by an open string should show equal
    and opposite charge (+3/-3) at the two endpoints: a genuine
    dislocation DIPOLE, with the string itself playing the role of the
    (unphysical, gauge/path-independent) Volterra cut.
"""
from collections import deque

import numpy as np

from rhombile_lattice import RhombileLattice


def enumerate_triangles_all(lattice):
    """enumerate_triangles without the interior-only wrap-exclusion filter
    (see rhombile_lattice.enumerate_triangles's own docstring): every bond
    -- including ones that cross the periodic seam -- has a perfectly
    well-defined *local* geometry (its own p1/p2), so per-triangle CCW
    order / color / increments are fine even when the triangle touches
    the seam. Only cross-triangle CENTROID averaging is unreliable there
    (irrelevant for a height function, which only uses site indices)."""
    bond_by_pair = {frozenset((b["i"], b["j"])): b for b in lattice.bonds}
    triangles = []
    for m in range(lattice.ny):
        for n in range(lattice.nx):
            r1_idx = lattice._site_index(n, m, "r1")
            r1_pos = lattice._position(n, m, "r1")
            nbrs = []
            for other_idx, b in lattice.neighbor_table()[r1_idx]:
                other_pos = b["p2"] if b["i"] == r1_idx else b["p1"]
                ang = np.arctan2(*(other_pos - r1_pos)[::-1])
                nbrs.append((ang, other_idx, b))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, b1 = nbrs[i]
                _, idx2, b2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None or closing["type"] != "r2r3":
                    continue
                edges = {b1["type"]: b1, b2["type"]: b2}
                if "r1r2" not in edges or "r1r3" not in edges:
                    continue
                triangles.append({
                    "sites": (r1_idx, idx1, idx2),
                    "r1r2": edges["r1r2"], "r1r3": edges["r1r3"], "r2r3": closing,
                })
    return triangles


def color_of(tri):
    """'A' if the CCW-first spoke from r1 (sites[0]->sites[1]) is the
    r1r2 edge, else 'B' (it's r1r3). Adjacent triangles (sharing any
    edge) always have opposite color -- the triangular lattice's
    inherent up/down bipartiteness -- which is what makes the height
    function below single-valued on every ordinary triangle."""
    c0, c1, _ = tri["sites"]
    return "A" if frozenset((tri["r1r2"]["i"], tri["r1r2"]["j"])) == frozenset((c0, c1)) else "B"


def f_of(bond):
    return 1.0 if bond["J"] != 0.0 else -2.0


def triangle_charge(tri):
    """Sum of increments going around this triangle's own 3 edges, in
    its own CCW order -- 0 for an ordinary triangle, sign*3 for a
    frustrated (all-active) one. Independent of any global height
    computation; a pure local consistency fact of the construction."""
    sign = 1.0 if color_of(tri) == "A" else -1.0
    return sign * (f_of(tri["r1r2"]) + f_of(tri["r1r3"]) + f_of(tri["r2r3"]))


def build_increments(lattice, triangles):
    """Directed increment dict {(i, j): value} with value = -dict[(j, i)],
    built from every triangle's 3 edges, checking that any bond shared by
    2 triangles gets the same value from both (should always hold given
    the color argument above; mismatches are returned for diagnosis)."""
    inc = {}
    mismatches = []
    for tri in triangles:
        c0, c1, c2 = tri["sites"]
        sign = 1.0 if color_of(tri) == "A" else -1.0
        for u, v, bond in ((c0, c1, tri["r1r2"] if color_of(tri) == "A" else tri["r1r3"]),
                           (c1, c2, tri["r2r3"]),
                           (c2, c0, tri["r1r3"] if color_of(tri) == "A" else tri["r1r2"])):
            val = sign * f_of(bond)
            if (u, v) in inc:
                if abs(inc[(u, v)] - val) > 1e-9:
                    mismatches.append((u, v, inc[(u, v)], val))
            else:
                inc[(u, v)] = val
                inc[(v, u)] = -val
    return inc, mismatches


def compute_height(lattice, increments, ref=0):
    """BFS spanning tree over every site reachable via a bond that has a
    recorded increment (site index already canonical mod nx,ny, so a
    wrap bond just looks like an ordinary edge here). Returns
    (h, tree_edges, nontree_mismatches): h[ref]=0; tree_edges is the set
    of (i,j) used to build h; nontree_mismatches lists every OTHER bond
    (i,j) with recorded increment where h[j]-h[i] != increments[(i,j)] --
    each such mismatch is the holonomy (Burgers charge) of the minimal
    cycle formed by the tree path i->...->ref->...->j plus this edge."""
    adj = {}
    for (i, j) in increments:
        adj.setdefault(i, []).append(j)

    h = {ref: 0.0}
    tree_edges = set()
    queue = deque([ref])
    while queue:
        u = queue.popleft()
        for v in adj.get(u, []):
            if v not in h:
                h[v] = h[u] + increments[(u, v)]
                tree_edges.add((u, v))
                tree_edges.add((v, u))
                queue.append(v)

    unreached = [s for s in range(lattice.n_sites) if s not in h]
    checks = []
    for (i, j), val in increments.items():
        if (i, j) in tree_edges:
            continue
        if i not in h or j not in h:
            continue
        actual = h[j] - h[i]
        checks.append((i, j, actual - val))
    return h, unreached, checks


def summarize(label, lattice):
    triangles = enumerate_triangles_all(lattice)
    inc, build_mismatches = build_increments(lattice, triangles)
    h, unreached, holonomies = compute_height(lattice, inc)
    ft_charges = sorted(set(round(triangle_charge(t), 6) for t in triangles) - {0.0})
    n_defect_tri = sum(1 for t in triangles if abs(triangle_charge(t)) > 1e-9)
    nonzero_holo = [(i, j, v) for i, j, v in holonomies if abs(v) > 1e-9]
    print(f"\n=== {label} ===")
    print(f"  {len(triangles)} triangles, {len(build_mismatches)} build-consistency "
          f"mismatches (should be 0)")
    print(f"  {n_defect_tri} frustrated triangles, per-triangle charges present: {ft_charges}")
    print(f"  {len(unreached)} unreached sites (should be 0)")
    print(f"  {len(nonzero_holo)} / {len(holonomies)} non-tree bonds show nonzero holonomy")
    if nonzero_holo:
        vals = [v for _, _, v in nonzero_holo]
        print(f"    holonomy values seen: {sorted(set(round(v, 3) for v in vals))}")
    return {
        "triangles": triangles, "increments": inc, "h": h,
        "build_mismatches": build_mismatches, "unreached": unreached,
        "holonomies": holonomies, "n_defect_tri": n_defect_tri,
    }


if __name__ == "__main__":
    NX, NY = 12, 12
    lat = RhombileLattice(NX, NY)
    summarize("pristine (no string defect)", lat)

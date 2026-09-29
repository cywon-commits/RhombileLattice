"""Rhombile lattice (PBC) construction, string-defect insertion, visualization.

Lattice: 3-site basis (r1, r2, r3) per unit cell, bond offsets and geometry
follow rhombile_topo_report.md / rhombile_classical_proposal.md.
"""
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow
import matplotlib.pyplot as plt

A1 = np.array([1.0, 0.0])
A2 = np.array([-0.5, np.sqrt(3) / 2])
A_MAT = np.column_stack([A1, A2])
A_MAT_INV = np.linalg.inv(A_MAT)

BASIS = {
    "r1": np.array([0.0, 0.0]),
    "r2": np.array([0.5, np.sqrt(3) / 6]),
    "r3": np.array([0.0, np.sqrt(3) / 3]),
}
SUB_INDEX = {"r1": 0, "r2": 1, "r3": 2}

# (site_a, site_b) -> offsets (dn, dm) from a(n, m) to b(n+dn, m+dm)
BOND_OFFSETS = {
    ("r1", "r2"): [(0, 0), (-1, 0), (-1, -1)],
    ("r1", "r3"): [(0, 0), (0, -1), (-1, -1)],
    ("r2", "r3"): [(0, 0), (0, -1), (1, 0)],
}

# range of periodic box translations tried when resolving a wrapping string;
# +-2 comfortably covers the minimum-image displacement for any nx, ny >= 3
WRAP_SHIFT_RANGE = range(-2, 3)


class RhombileLattice:
    def __init__(self, nx, ny, j12=-1.0, j13=-1.0, j23=0.0):
        self.nx, self.ny = nx, ny
        self.j0 = {"r1r2": j12, "r1r3": j13, "r2r3": j23}
        self.n_sites = nx * ny * 3
        self.bonds = []  # dicts: i, j, type, J, p1, p2 (unwrapped positions)
        self._neighbors = None  # lazy cache, see neighbor_table()
        self._build()

    def _site_index(self, n, m, sub):
        return ((m % self.ny) * self.nx + (n % self.nx)) * 3 + SUB_INDEX[sub]

    def _position(self, n, m, sub):
        return n * A1 + m * A2 + BASIS[sub]

    def _build(self):
        for m in range(self.ny):
            for n in range(self.nx):
                for (s1, s2), offsets in BOND_OFFSETS.items():
                    i = self._site_index(n, m, s1)
                    p1 = self._position(n, m, s1)
                    for dn, dm in offsets:
                        j = self._site_index(n + dn, m + dm, s2)
                        p2 = self._position(n + dn, m + dm, s2)
                        wraps = not (0 <= n + dn < self.nx and 0 <= m + dm < self.ny)
                        self.bonds.append({
                            "i": i, "j": j, "type": s1 + s2, "wraps": wraps,
                            "J": self.j0[s1 + s2], "p1": p1, "p2": p2,
                        })

    def site_positions(self):
        pos = np.zeros((self.n_sites, 2))
        sub_of = np.empty(self.n_sites, dtype=object)
        for m in range(self.ny):
            for n in range(self.nx):
                for sub in BASIS:
                    idx = self._site_index(n, m, sub)
                    pos[idx] = self._position(n, m, sub)
                    sub_of[idx] = sub
        return pos, sub_of

    def neighbor_table(self):
        """list[site] -> [(other_site, bond_dict), ...], each bond referenced
        (not copied) so a later string-defect J flip is seen automatically."""
        if self._neighbors is None:
            nbrs = [[] for _ in range(self.n_sites)]
            for b in self.bonds:
                nbrs[b["i"]].append((b["j"], b))
                nbrs[b["j"]].append((b["i"], b))
            self._neighbors = nbrs
        return self._neighbors


def random_states(lattice, rng):
    """Random initial Potts state (0, 1, or 2) for every site."""
    return rng.integers(0, 3, size=lattice.n_sites)


def total_energy(lattice, states, D):
    """Anisotropic 3-state Potts Hamiltonian, generalizing H = -J S_i.S_j to
    Potts variables via delta(s_i, s_j) in place of the dot product:

        H = sum_<i,j> (-J_ij) * delta(s_i, s_j)  +  sum_i D[s_i]

    J_ij < 0 (the lattice default) therefore means antiferromagnetic:
    matching neighbors cost +|J_ij|, mismatched neighbors cost 0.
    D = (D_1, D_2, D_3) is the single-site anisotropy energy of each state.
    """
    e = float(np.asarray(D, dtype=float)[states].sum())
    for b in lattice.bonds:
        if states[b["i"]] == states[b["j"]]:
            e -= b["J"]
    return e


def state_energies(lattice, states, site, D):
    """Energy of `site` for each candidate state (0, 1, 2), holding every
    other site fixed -- the effective field used to drive Metropolis/SA
    updates. The cost of flipping `site` to state k is

        delta_E = state_energies(...)[k] - state_energies(...)[states[site]]
    """
    e = np.asarray(D, dtype=float).copy()
    for j, b in lattice.neighbor_table()[site]:
        e[states[j]] -= b["J"]
    return e


def frustrated_triangles(lattice):
    """Enumerate every elementary (r1, r2, r3) triangle and return the ones
    that are topologically frustrated: all 3 edges active (nonzero J). Such
    a triangle is an odd all-AF cycle -- no spin assignment can satisfy
    every bond on it, regardless of what the rest of the lattice does.

    Returns a list of dicts {"r1", "a", "b", "centroid"} (site indices of
    the triangle's corners and its plotting centroid).
    """
    bond_by_pair = {frozenset((b["i"], b["j"])): b for b in lattice.bonds}
    result = []
    for m in range(lattice.ny):
        for n in range(lattice.nx):
            r1_idx = lattice._site_index(n, m, "r1")
            r1_pos = lattice._position(n, m, "r1")
            nbrs = []
            for other_idx, b in lattice.neighbor_table()[r1_idx]:
                other_pos = b["p2"] if b["i"] == r1_idx else b["p1"]
                ang = np.arctan2(*(other_pos - r1_pos)[::-1])
                nbrs.append((ang, other_idx, other_pos, b["J"]))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, pos1, J1 = nbrs[i]
                _, idx2, pos2, J2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None:
                    continue
                if J1 != 0 and J2 != 0 and closing["J"] != 0:
                    centroid = (r1_pos + pos1 + pos2) / 3
                    result.append({"r1": r1_idx, "a": idx1, "b": idx2, "centroid": centroid})
    return result


def matched_bonds(lattice, states):
    """Bonds whose two endpoints share the same state: the ones actually
    contributing +|J| to total_energy (the real physical cost), as opposed
    to bonds that were merely toggled when a string defect was drawn."""
    return [b for b in lattice.bonds if b["J"] != 0.0 and states[b["i"]] == states[b["j"]]]


def enumerate_triangles(lattice):
    """Every elementary (r1, r2, r3) triangle whose 3 edges all stay
    inside the box (regardless of J), keyed by edge type. Each triangle
    has exactly one r1-r2 edge, one r1-r3 edge (the two default-active
    "outer" edges of a rhombus), and one r2-r3 edge (the default-inactive
    "diagonal" of that rhombus).

    Triangles that touch a periodic-wraparound bond are left out: with
    PBC, a bond crossing the box edge is a real short bond on the torus,
    but its stored p1/p2 uses one endpoint's *unwrapped* position (see
    RhombileLattice._build), so such a triangle's naive centroid can land
    a full box-width away from where it physically is. That's harmless
    for frustrated_triangles (which never averages positions across
    triangles), but build_rhombi/build_triangle_hop_graph do, so string
    routing is restricted to the box's interior, away from that seam.

    Returns a list of dicts {"sites", "r1r2", "r1r3", "r2r3", "centroid"};
    the 3 edge entries are bond dicts shared by reference with
    lattice.bonds, so toggling their "J" is visible everywhere else.
    """
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
                nbrs.append((ang, other_idx, other_pos, b))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, pos1, b1 = nbrs[i]
                _, idx2, pos2, b2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None or closing["type"] != "r2r3":
                    continue
                edges = {b1["type"]: b1, b2["type"]: b2, "r2r3": closing}
                if "r1r2" not in edges or "r1r3" not in edges:
                    continue
                if edges["r1r2"]["wraps"] or edges["r1r3"]["wraps"] or closing["wraps"]:
                    continue
                centroid = (r1_pos + pos1 + pos2) / 3
                triangles.append({
                    "sites": (r1_idx, idx1, idx2), "r1r2": edges["r1r2"],
                    "r1r3": edges["r1r3"], "r2r3": edges["r2r3"], "centroid": centroid,
                })
    return triangles


def build_rhombi(triangles):
    """Pair triangles that share a diagonal (r2-r3 edge) into rhombi, and
    connect rhombi that share an outer edge (r1-r2 or r1-r3) into an
    adjacency graph -- the dual graph a "string" is routed through.

    Returns (rhombi, adjacency): rhombi[i] = {"halves": (tri_a, tri_b),
    "diagonal": bond, "centroid"}; adjacency[i] = [(j, shared_outer_bond), ...].
    """
    by_diagonal = {}
    for ti, tri in enumerate(triangles):
        by_diagonal.setdefault(id(tri["r2r3"]), []).append(ti)

    rhombi = []
    tri_to_rhombus = {}
    for tri_ids in by_diagonal.values():
        if len(tri_ids) != 2:
            continue  # only possible with a too-small non-periodic patch
        rid = len(rhombi)
        centroid = np.mean([triangles[t]["centroid"] for t in tri_ids], axis=0)
        rhombi.append({
            "halves": tuple(tri_ids), "diagonal": triangles[tri_ids[0]]["r2r3"],
            "centroid": centroid,
        })
        for t in tri_ids:
            tri_to_rhombus[t] = rid

    by_outer = {}
    for ti, tri in enumerate(triangles):
        if ti not in tri_to_rhombus:
            continue
        for key in ("r1r2", "r1r3"):
            by_outer.setdefault(id(tri[key]), []).append((ti, tri[key]))

    adjacency = [[] for _ in rhombi]
    for entries in by_outer.values():
        if len(entries) != 2:
            continue
        (t1, b1), (t2, _) = entries
        r1_id, r2_id = tri_to_rhombus[t1], tri_to_rhombus[t2]
        if r1_id != r2_id:
            adjacency[r1_id].append((r2_id, b1))
            adjacency[r2_id].append((r1_id, b1))
    return rhombi, adjacency


def nearest_rhombus(rhombi, point):
    """Index of the rhombus whose centroid is closest to `point`."""
    point = np.asarray(point, dtype=float)
    dists = [np.linalg.norm(r["centroid"] - point) for r in rhombi]
    return int(np.argmin(dists))


def nearest_triangle(triangles, point):
    """Index of the triangle (rhombus half) whose centroid is closest to
    `point` -- the fine-grained counterpart of nearest_rhombus, needed
    because a string's routing must pick a specific half, not just a
    rhombus (see build_triangle_hop_graph)."""
    point = np.asarray(point, dtype=float)
    dists = [np.linalg.norm(t["centroid"] - point) for t in triangles]
    return int(np.argmin(dists))


def build_triangle_hop_graph(triangles):
    """The graph a string is actually routed through, one level finer than
    build_rhombi's rhombus graph: nodes are triangle halves, and crossing
    a rhombus's diagonal to its sibling half is forced, not a free choice.

    Naively routing on rhombi (pick any of a rhombus's up to 4 outer
    edges to enter/exit by) can enter and leave through the *same* half,
    leaving the other, untouched half with a bare diagonal flip and no
    compensating outer-edge change -- a spurious frustrated triangle in
    the middle of the path. Forcing every hop to land on the sibling half
    (this function) rules that out structurally: every interior half ends
    up with exactly one outer edge toggled plus the (shared) diagonal
    toggle, i.e. still 2-active/1-inactive, just with a different edge
    inactive -- the domain-wall picture -- so frustration is confined to
    the two ends of the path regardless of its shape.

    Returns (sibling, hop): sibling[t] = t's diagonal partner; hop[t] =
    [(next_far_half, outer_bond, entered_near_half), ...], where
    next_far_half is already the sibling of whatever near half the outer
    bond leads to (the forced cross is baked in).
    """
    sibling = {}
    diag_groups = {}
    for ti, tri in enumerate(triangles):
        diag_groups.setdefault(id(tri["r2r3"]), []).append(ti)
    for tis in diag_groups.values():
        if len(tis) == 2:
            a, b = tis
            sibling[a], sibling[b] = b, a

    outer_groups = {}
    outer_bond_of = {}
    for ti, tri in enumerate(triangles):
        for key in ("r1r2", "r1r3"):
            bond = tri[key]
            outer_groups.setdefault(id(bond), []).append(ti)
            outer_bond_of[id(bond)] = bond

    hop = [[] for _ in triangles]
    for bond_id, tis in outer_groups.items():
        if len(tis) != 2:
            continue
        ta, tb = tis
        if ta not in sibling or tb not in sibling:
            continue
        bond = outer_bond_of[bond_id]
        hop[ta].append((sibling[tb], bond, tb))
        hop[tb].append((sibling[ta], bond, ta))
    return sibling, hop


def route_string(sibling, hop, t_start, t_end):
    """Shortest path from t_start to t_end through the triangle hop graph.

    t_start and t_end are left as the two path *endpoints* (the ones that
    end up frustrated); the walk itself starts at sibling[t_start] (the
    half that is immediately ready to hop onward) and must arrive
    exactly at t_end. Returns (nodes, bonds): nodes = every rhombus-half
    visited (whose diagonal gets toggled once each, including the
    endpoints), bonds = the outer edges crossed between them.
    """
    start_node = sibling[t_start]
    if start_node == t_end:
        return [t_end], []
    prev, prev_bond = {start_node: None}, {}
    queue = deque([start_node])
    while queue:
        cur = queue.popleft()
        if cur == t_end:
            break
        for nxt, bond, _near_half in hop[cur]:
            if nxt not in prev:
                prev[nxt] = cur
                prev_bond[nxt] = bond
                queue.append(nxt)
    if t_end not in prev:
        raise ValueError("no path between the two triangles (graph disconnected?)")
    nodes, bonds, node = [t_end], [], t_end
    while prev[node] is not None:
        bonds.append(prev_bond[node])
        node = prev[node]
        nodes.append(node)
    nodes.reverse()
    bonds.reverse()
    return nodes, bonds


def route_string_between_points(rhombi, sibling, hop, p_start, p_end):
    """Convenience wrapper around route_string for two arbitrary points.

    Each rhombus's 2 triangle halves have opposite up/down parity, and
    route_string can only connect same-parity endpoints (crossing a
    diagonal then an outer edge always returns to the starting parity),
    so only 2 of the 4 (half-of-start, half-of-end) combinations are
    actually reachable. This tries all 4 and keeps the shortest that
    works. Returns (t_start, t_end, nodes, bonds).
    """
    rid_s = nearest_rhombus(rhombi, p_start)
    rid_e = nearest_rhombus(rhombi, p_end)
    best = None
    for ts in rhombi[rid_s]["halves"]:
        for te in rhombi[rid_e]["halves"]:
            try:
                nodes, bonds = route_string(sibling, hop, ts, te)
            except ValueError:
                continue
            if best is None or len(nodes) < len(best[2]):
                best = (ts, te, nodes, bonds)
    if best is None:
        raise ValueError("no valid string between these two rhombi")
    return best


def _toggle_bond(b):
    b["J"] = 0.0 if b["J"] != 0.0 else -1.0


def apply_dual_string_defect(triangles, nodes, bonds):
    """Toggle ON the diagonal of every rhombus-half in `nodes` (one
    toggle per rhombus, since a diagonal is shared by both its halves)
    and toggle OFF every outer edge in `bonds` connecting consecutive
    halves. See route_string / build_triangle_hop_graph for why this
    confines the resulting frustration to exactly nodes[0] and nodes[-1].
    """
    flipped = []
    for t in nodes:
        d = triangles[t]["r2r3"]
        _toggle_bond(d)
        flipped.append(d)
    for b in bonds:
        _toggle_bond(b)
        flipped.append(b)
    return flipped


def min_state2_assignment(lattice):
    """The actual D3=0 ground state for a diluted/reconnected lattice:
    r1 sites fixed at state 0 (always mismatches every r1-rim bond, active
    or not); every rim (r2/r3) site defaults to state 1. The only bonds
    that can still cost anything are *active* r2-r3 (diagonal) edges,
    since both endpoints default to 1. Recoloring one endpoint of each to
    the otherwise-unused state 2 clears it for free -- as long as the two
    recolored endpoints of any two such edges never coincide by being
    forced onto opposite requirements, which is exactly the 2-coloring
    (bipartition) of the graph formed by those diagonal edges.

    A naive minimum *vertex cover* of that graph is NOT the same thing --
    it only guarantees every edge has >=1 endpoint in the cover, not that
    the two endpoints differ, so 2 cover vertices sharing an edge would
    silently reintroduce a cost. Bipartition (proper 2-coloring, taking
    the smaller color class per component as state 2) is both correct and
    minimal, since each connected component's smaller side is the fewest
    sites that can possibly cover its edges 1-for-1.

    Returns (states, n_state2, n_non_bipartite_components). A nonzero
    last value means that component contains an odd cycle of active
    diagonals -- no zero-cost 2-coloring exists there, and the true D3=0
    ground state (found separately, e.g. by SA) has residual energy.
    """
    edges = [(b["i"], b["j"]) for b in lattice.bonds if b["type"] == "r2r3" and b["J"] != 0.0]
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    pos, sub_of = lattice.site_positions()
    states = np.where(sub_of == "r1", 0, 1).astype(int)

    color, non_bipartite = {}, 0
    for start in adj:
        if start in color:
            continue
        color[start] = 0
        comp, queue = [start], deque([start])
        ok = True
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    comp.append(v)
                    queue.append(v)
                elif color[v] == color[u]:
                    ok = False
        if not ok:
            non_bipartite += 1
        side0 = [v for v in comp if color[v] == 0]
        side1 = [v for v in comp if color[v] == 1]
        for v in (side0 if len(side0) <= len(side1) else side1):
            states[v] = 2
    return states, int(np.sum(states == 2)), non_bipartite


def conflict_graph_bipartition(lattice):
    """2-color the graph of active r2-r3 (diagonal) bonds -- the same
    graph min_state2_assignment works with, but returning the raw edges
    plus a single global 2-coloring (one side per color, pooled across
    every connected component) rather than picking the smaller side per
    component. That per-component choice is what makes D3=0 minimal;
    a global 2-coloring is what energy_via_mincut needs instead, since
    it lets the flow solver pick the right combination itself for
    whatever D3 is asked for.

    Returns (edges, side_a, side_b, non_bipartite) -- same meaning as
    in min_state2_assignment.
    """
    edges = [(b["i"], b["j"]) for b in lattice.bonds if b["type"] == "r2r3" and b["J"] != 0.0]
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    color, non_bipartite = {}, 0
    for start in adj:
        if start in color:
            continue
        color[start] = 0
        queue, ok = deque([start]), True
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    queue.append(v)
                elif color[v] == color[u]:
                    ok = False
        if not ok:
            non_bipartite += 1
    side_a = [v for v, c in color.items() if c == 0]
    side_b = [v for v, c in color.items() if c == 1]
    return edges, side_a, side_b, non_bipartite


def energy_via_mincut(lattice, D3, scale=10000):
    """D3=0 always gives the true exact energy (0), which is provably
    optimal regardless of anything below. For D3>0, treat the result as
    an UPPER BOUND, not a confirmed exact minimum: this function fixes
    every r1 site to Potts state 0 and only lets rim sites choose
    between 1 and 2, and a later investigation
    (exact_ground_state_investigation.py) found configurations where
    letting r1 vary too gives a strictly lower energy at large D3 (a
    bent/detour string: 15 found vs. 53 claimed here). An attempted fix
    via Toulouse/Barahona planar-matching theory produced smaller
    numbers but could not be verified with a working reconstruction,
    and simulated annealing sided with THIS function's numbers at one
    clearly-discriminating test point. Net effect: this is the best
    verified upper bound at D3>0, but neither this nor the "fix" has a
    fully confirmed lower bound to match it -- see that file's docstring
    before trusting D3>0 numbers from this function in anything meant
    to be exact.

    Exact ground-state energy at a given D3, for a lattice whose
    active-diagonal conflict graph is bipartite (checked; raises if not,
    since the reduction below assumes it).

    Minimizing E(x) = D3*sum(x_v) + sum_edges[x_u == x_v] over x in
    {0,1}^rim (x_v=1 meaning "site v is state 2") is an antiferromagnetic
    binary MRF -- NP-hard in general, but exactly solvable in polynomial
    time here because the conflict graph is bipartite: relabeling one
    side (y_v = x_v on side A, y_v = 1-x_v on side B) turns it into a
    *ferromagnetic* (submodular) energy, which reduces to a standard
    minimum s-t cut: source->v capacity D3 for v in A, v->sink capacity
    D3 for v in B, and capacity 1 both ways on every original edge
    (see rhombile_lattice module notes / the research doc for the
    algebra). The reduction is exact -- min-cut value equals E(D3)
    with no leftover additive constant -- and was checked against a
    brute-force search over all state-2 subsets on a small case.

    This is what makes it possible to sweep D3 on conflict graphs with
    hundreds of vertices (e.g. a dense grid of many strings), where the
    2^n brute force used for the first single-string check is hopeless.
    """
    edges, side_a, side_b, non_bipartite = conflict_graph_bipartition(lattice)
    if non_bipartite:
        raise ValueError(
            f"conflict graph has {non_bipartite} non-bipartite component(s); "
            "the min-cut reduction here assumes bipartite (no odd cycles)."
        )
    idx_a = {v: i + 1 for i, v in enumerate(side_a)}
    idx_b = {v: i + 1 + len(side_a) for i, v in enumerate(side_b)}
    n = len(side_a) + len(side_b) + 2
    source, sink = 0, n - 1
    rows, cols, caps = [], [], []
    cap_d3 = int(round(D3 * scale))
    for v in side_a:
        rows.append(source); cols.append(idx_a[v]); caps.append(cap_d3)
    for v in side_b:
        rows.append(idx_b[v]); cols.append(sink); caps.append(cap_d3)
    for a, b in edges:
        u, w = (a, b) if a in idx_a else (b, a)
        rows.append(idx_a[u]); cols.append(idx_b[w]); caps.append(scale)
        rows.append(idx_b[w]); cols.append(idx_a[u]); caps.append(scale)
    capacity = csr_matrix((caps, (rows, cols)), shape=(n, n))
    return maximum_flow(capacity, source, sink).flow_value / scale


def state2_and_violations_via_mincut(lattice, D3, scale=10000):
    """See energy_via_mincut's docstring first: for D3>0 this shares its
    r1-fixed restriction and the same open reliability question (exact
    at D3=0 only; an upper bound, not a confirmed exact value, above
    that). exact_ground_state_investigation.py has the details.

    Same exact optimum as energy_via_mincut, but split into its two
    physical pieces instead of just the total: n_state2 (node-type
    excitations -- sites that took the 3rd Potts state) and n_violated
    (link-type excitations -- active diagonals whose endpoints still
    match). E(D3) = D3*n_state2 + n_violated exactly.

    Recovers the actual minimum cut (not just its value) by BFS-ing the
    residual capacity graph from the source; the reachable set is one
    side of the optimal cut.

    Returns (n_state2, n_violated).
    """
    edges, side_a, side_b, non_bipartite = conflict_graph_bipartition(lattice)
    if non_bipartite:
        raise ValueError(
            f"conflict graph has {non_bipartite} non-bipartite component(s); "
            "the min-cut reduction here assumes bipartite (no odd cycles)."
        )
    idx_a = {v: i + 1 for i, v in enumerate(side_a)}
    idx_b = {v: i + 1 + len(side_a) for i, v in enumerate(side_b)}
    n = len(side_a) + len(side_b) + 2
    source, sink = 0, n - 1
    rows, cols, caps = [], [], []
    cap_d3 = int(round(D3 * scale))
    for v in side_a:
        rows.append(source); cols.append(idx_a[v]); caps.append(cap_d3)
    for v in side_b:
        rows.append(idx_b[v]); cols.append(sink); caps.append(cap_d3)
    for a, b in edges:
        u, w = (a, b) if a in idx_a else (b, a)
        rows.append(idx_a[u]); cols.append(idx_b[w]); caps.append(scale)
        rows.append(idx_b[w]); cols.append(idx_a[u]); caps.append(scale)
    capacity = csr_matrix((caps, (rows, cols)), shape=(n, n))
    result = maximum_flow(capacity, source, sink)
    residual = (capacity - result.flow).tocsr()

    reachable = {source}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        row = residual.getrow(u)
        for v, c in zip(row.indices, row.data):
            if c > 1e-9 and v not in reachable:
                reachable.add(v)
                queue.append(v)

    n_state2 = (sum(1 for v in side_a if idx_a[v] not in reachable)
                + sum(1 for v in side_b if idx_b[v] in reachable))
    n_violated = 0
    for a, b in edges:
        u, w = (a, b) if a in idx_a else (b, a)
        x_u = idx_a[u] not in reachable
        x_w = idx_b[w] in reachable
        if x_u == x_w:
            n_violated += 1
    return n_state2, n_violated


def full_state_via_mincut(lattice, D3, scale=10000):
    """See energy_via_mincut's docstring first: for D3>0 this shares its
    r1-fixed restriction and the same open reliability question. Exact
    at D3=0 only.

    Same exact optimum as state2_and_violations_via_mincut, but
    returning the full per-site state assignment and the actual list of
    violated (still-matching, active r2-r3) bonds, rather than just their
    counts -- what you need to draw the ground state, not just quote its
    energy.

    Returns (states, violated_bonds).
    """
    edges, side_a, side_b, non_bipartite = conflict_graph_bipartition(lattice)
    if non_bipartite:
        raise ValueError(
            f"conflict graph has {non_bipartite} non-bipartite component(s); "
            "the min-cut reduction here assumes bipartite (no odd cycles)."
        )
    idx_a = {v: i + 1 for i, v in enumerate(side_a)}
    idx_b = {v: i + 1 + len(side_a) for i, v in enumerate(side_b)}
    n = len(side_a) + len(side_b) + 2
    source, sink = 0, n - 1
    rows, cols, caps = [], [], []
    cap_d3 = int(round(D3 * scale))
    for v in side_a:
        rows.append(source); cols.append(idx_a[v]); caps.append(cap_d3)
    for v in side_b:
        rows.append(idx_b[v]); cols.append(sink); caps.append(cap_d3)
    for a, b in edges:
        u, w = (a, b) if a in idx_a else (b, a)
        rows.append(idx_a[u]); cols.append(idx_b[w]); caps.append(scale)
        rows.append(idx_b[w]); cols.append(idx_a[u]); caps.append(scale)
    capacity = csr_matrix((caps, (rows, cols)), shape=(n, n))
    result = maximum_flow(capacity, source, sink)
    residual = (capacity - result.flow).tocsr()

    reachable = {source}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        row = residual.getrow(u)
        for v, c in zip(row.indices, row.data):
            if c > 1e-9 and v not in reachable:
                reachable.add(v)
                queue.append(v)

    pos, sub_of = lattice.site_positions()
    states = np.where(sub_of == "r1", 0, 1).astype(int)
    for v in side_a:
        if idx_a[v] not in reachable:
            states[v] = 2
    for v in side_b:
        if idx_b[v] in reachable:
            states[v] = 2

    violated_bonds = [b for b in lattice.bonds
                       if b["type"] == "r2r3" and b["J"] != 0.0 and states[b["i"]] == states[b["j"]]]
    return states, violated_bonds


def heat_bath_sweep(lattice, states, D, T, rng):
    """One sweep = resample every site's state from the Gibbs distribution
    p(k) ~ exp(-state_energies(...)[k] / T), in random order, in place."""
    for site in rng.permutation(lattice.n_sites):
        e = state_energies(lattice, states, site, D)
        w = np.exp(-(e - e.min()) / T)
        states[site] = rng.choice(3, p=w / w.sum())
    return states


def simulated_annealing(lattice, D, rng, n_sweeps=400, T_start=3.0, T_end=0.01,
                         states=None, record_energy=False):
    """Anneal from T_start down to T_end over n_sweeps heat-bath sweeps
    (geometric schedule). Returns the final states, and the energy trace
    (per sweep) if record_energy=True."""
    if states is None:
        states = random_states(lattice, rng)
    else:
        states = states.copy()
    energies = [] if record_energy else None
    for T in np.geomspace(T_start, T_end, n_sweeps):
        heat_bath_sweep(lattice, states, D, T, rng)
        if record_energy:
            energies.append(total_energy(lattice, states, D))
    return states, energies


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def segments_intersect(p1, p2, p3, p4):
    d1, d2 = _cross(p3, p4, p1), _cross(p3, p4, p2)
    d3, d4 = _cross(p1, p2, p3), _cross(p1, p2, p4)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def minimum_image_endpoint(lattice, p_start, p_end):
    """Periodic image of p_end closest to p_start, i.e. the endpoint of the
    shortest string between the two points on the torus."""
    frac_start = A_MAT_INV @ p_start
    frac_end = A_MAT_INV @ p_end
    du, dv = frac_end - frac_start
    du -= lattice.nx * np.round(du / lattice.nx)
    dv -= lattice.ny * np.round(dv / lattice.ny)
    return p_start + du * A1 + dv * A2


def apply_string_defect(lattice, p_start, p_end, wrap=True):
    """Toggle J (-1 <-> 0) on every bond crossed by the string p_start -> p_end.

    With wrap=True, p_end is first replaced by its minimum-image copy (the
    shortest string on the torus), and the crossing test is repeated over a
    small set of periodic box translations, so a string that runs off one
    edge of the simulation box and reappears on the opposite edge is handled
    correctly. With wrap=False the string is treated as a plain segment in
    the plane (no PBC wrap-around).
    """
    if wrap:
        p_end = minimum_image_endpoint(lattice, p_start, p_end)
        shifts = [kx * lattice.nx * A1 + ky * lattice.ny * A2
                  for kx in WRAP_SHIFT_RANGE for ky in WRAP_SHIFT_RANGE]
    else:
        shifts = [np.zeros(2)]

    flipped = []
    for shift in shifts:
        s_start, s_end = p_start - shift, p_end - shift
        for b in lattice.bonds:
            if segments_intersect(s_start, s_end, b["p1"], b["p2"]):
                b["J"] = 0.0 if b["J"] != 0.0 else -1.0
                flipped.append(b)
    return flipped, p_end


def clip_string_for_display(lattice, p_start, p_end_image, min_length=0.05):
    """Split the (possibly out-of-box) string p_start -> p_end_image at box
    boundaries and translate each piece back into the fundamental box, so a
    wrapping string can be drawn as separate segments re-entering the plot.

    When the string grazes a corner of the box, u and v can cross their
    periodic boundaries at two slightly different parameter values, briefly
    assigning a sliver of the path to a diagonally-adjacent cell copy; that
    sliver is real (the straight line genuinely passes through it) but only
    spans a fraction of a bond length, so pieces shorter than min_length are
    dropped here to avoid a stray dash far from the rest of the string."""
    frac0 = A_MAT_INV @ p_start
    frac1 = A_MAT_INV @ p_end_image
    du, dv = frac1 - frac0
    ts = {0.0, 1.0}
    if du != 0:
        lo, hi = sorted([frac0[0], frac1[0]])
        for k in range(int(np.floor(lo / lattice.nx)), int(np.ceil(hi / lattice.nx)) + 1):
            t = (k * lattice.nx - frac0[0]) / du
            if 0 < t < 1:
                ts.add(t)
    if dv != 0:
        lo, hi = sorted([frac0[1], frac1[1]])
        for k in range(int(np.floor(lo / lattice.ny)), int(np.ceil(hi / lattice.ny)) + 1):
            t = (k * lattice.ny - frac0[1]) / dv
            if 0 < t < 1:
                ts.add(t)
    ts = sorted(ts)

    segments = []
    for t0, t1 in zip(ts[:-1], ts[1:]):
        tm = 0.5 * (t0 + t1)
        fm = frac0 + tm * np.array([du, dv])
        kx, ky = int(np.floor(fm[0] / lattice.nx)), int(np.floor(fm[1] / lattice.ny))
        shift = kx * lattice.nx * A1 + ky * lattice.ny * A2
        q0 = p_start + t0 * (p_end_image - p_start) - shift
        q1 = p_start + t1 * (p_end_image - p_start) - shift
        if np.linalg.norm(q1 - q0) >= min_length:
            segments.append((q0, q1))
    return segments


def random_site_pair(lattice, rng):
    """Two random points = two random lattice sites."""
    pos, _ = lattice.site_positions()
    i, j = rng.choice(lattice.n_sites, size=2, replace=False)
    return pos[i], pos[j]


def random_point_pair(lattice, rng):
    """Two random points = two arbitrary continuous points inside the box."""
    frac = rng.uniform(0, 1, size=(2, 2)) * [lattice.nx, lattice.ny]
    return frac[0, 0] * A1 + frac[0, 1] * A2, frac[1, 0] * A1 + frac[1, 1] * A2


def plot_lattice(lattice, string_endpoints=None, wrap=True, ax=None, title=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    for b in lattice.bonds:
        if b["J"] != 0.0:
            ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                     color="black", linewidth=1.2, zorder=1)
    pos, sub_of = lattice.site_positions()
    colors = {"r1": "tab:red", "r2": "tab:blue", "r3": "tab:green"}
    for sub, c in colors.items():
        mask = sub_of == sub
        ax.scatter(pos[mask, 0], pos[mask, 1], s=25, color=c, label=sub, zorder=2)
    if string_endpoints is not None:
        p_start, p_end = string_endpoints
        p_end_image = minimum_image_endpoint(lattice, p_start, p_end) if wrap else p_end
        for k, (q0, q1) in enumerate(clip_string_for_display(lattice, p_start, p_end_image)):
            ax.plot([q0[0], q1[0]], [q0[1], q1[1]], "--", color="orange",
                     linewidth=2, zorder=3, label="string" if k == 0 else None)
        ax.scatter(*p_start, marker="x", color="orange", s=100, zorder=4)
        ax.scatter(*p_end, marker="x", color="orange", s=100, zorder=4)
    ax.set_aspect("equal")
    ax.legend(loc="upper right")
    if title:
        ax.set_title(title)
    return ax


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    nx, ny = 8, 8

    # 1. site-based string endpoints, PBC wrap allowed
    lat = RhombileLattice(nx, ny)
    p_start, p_end = random_site_pair(lat, rng)
    flipped, _ = apply_string_defect(lat, p_start, p_end, wrap=True)
    print(f"[site pair] {p_start} -> {p_end}: flipped {len(flipped)} bonds")
    ax = plot_lattice(lat, string_endpoints=(p_start, p_end),
                       title=f"Rhombile {nx}x{ny} (PBC): site-pair string")
    plt.savefig("rhombile_site_string.png", dpi=150)
    plt.close()

    # 2. continuous random-point string endpoints
    lat2 = RhombileLattice(nx, ny)
    p_start2, p_end2 = random_point_pair(lat2, rng)
    flipped2, _ = apply_string_defect(lat2, p_start2, p_end2, wrap=True)
    print(f"[continuous pair] {p_start2} -> {p_end2}: flipped {len(flipped2)} bonds")
    ax = plot_lattice(lat2, string_endpoints=(p_start2, p_end2),
                       title=f"Rhombile {nx}x{ny} (PBC): continuous-point string")
    plt.savefig("rhombile_continuous_string.png", dpi=150)
    plt.close()

    # 3. deliberately wrapping string (endpoints near opposite edges)
    lat3 = RhombileLattice(nx, ny)
    p_a = lat3._position(1, ny // 2, "r1")
    p_b = lat3._position(nx - 2, ny // 2, "r1")
    flipped3, p_b_image = apply_string_defect(lat3, p_a, p_b, wrap=True)
    print(f"[wrap demo] {p_a} -> {p_b} (image {p_b_image}): flipped {len(flipped3)} bonds")
    ax = plot_lattice(lat3, string_endpoints=(p_a, p_b),
                       title=f"Rhombile {nx}x{ny} (PBC): wrap-around string")
    plt.savefig("rhombile_wrap_string.png", dpi=150)
    plt.close()

    print("saved rhombile_site_string.png, rhombile_continuous_string.png, rhombile_wrap_string.png")

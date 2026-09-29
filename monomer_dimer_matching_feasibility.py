"""Direct combinatorial re-test of the density question, per the user's
own better idea: instead of building up isolated monomers one string at
a time (which stalled at p~4-5% due to string-PATH crossing artifacts,
not the actual combinatorial limit -- see isolated_monomer_density_
feasibility.py), fix a target monomer SET directly (an independent set,
balanced across the two bipartite classes) and ask the pure graph
question -- does the induced remainder graph admit a perfect matching?
-- via max-flow, the same tool this project already uses for its D3=0
ground-state construction (energy_via_mincut / full_state_via_mincut).

This is exactly the "mutilated chessboard" question the user pointed
at: the triangle-adjacency graph is bipartite (up/down triangles), so
removing an UNBALANCED monomer set makes a perfect matching of the rest
provably impossible (same parity argument as the classic two-missing-
same-color-squares puzzle) -- balance is necessary but not sufficient;
this checks the real (Hall's-theorem) condition directly via max-flow,
with no string-routing/crossing issues at all.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

from rhombile_lattice import RhombileLattice

NX, NY = 16, 16
DENSITIES = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
SEED = 3


def enumerate_all_triangles(lattice):
    """Every elementary (r1,r2,r3) triangle, INCLUDING ones that touch a
    periodic-wraparound bond (unlike enumerate_triangles, which drops
    those for centroid/plotting reasons this script doesn't need) --
    the genuine torus triangle count, always even and evenly split
    between the two bipartite classes, matching frustrated_triangles'
    own neighbor_table-based logic."""
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
                edges = {b1["type"]: b1, b2["type"]: b2, "r2r3": closing}
                if "r1r2" not in edges or "r1r3" not in edges:
                    continue
                triangles.append({"sites": (r1_idx, idx1, idx2),
                                   "r1r2": edges["r1r2"], "r1r3": edges["r1r3"], "r2r3": closing})
    return triangles


def build_triangle_adjacency(triangles):
    """adj[t] = set of triangle indices sharing an edge (any of the 3
    edge types) with t -- the full honeycomb triangle-adjacency graph,
    not just the forced diagonal-sibling relation."""
    n = len(triangles)
    adj = [set() for _ in range(n)]
    for key in ("r1r2", "r1r3", "r2r3"):
        groups = {}
        for ti, tri in enumerate(triangles):
            groups.setdefault(id(tri[key]), []).append(ti)
        for tis in groups.values():
            if len(tis) == 2:
                a, b = tis
                adj[a].add(b)
                adj[b].add(a)
    return adj


def bipartition(adj):
    """2-color the (bipartite, by construction) triangle-adjacency graph."""
    n = len(adj)
    color = [-1] * n
    for start in range(n):
        if color[start] != -1:
            continue
        color[start] = 0
        stack = [start]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1 - color[u]
                    stack.append(v)
                elif color[v] == color[u]:
                    raise ValueError("triangle-adjacency graph is not bipartite (unexpected)")
    return color


def greedy_balanced_independent_set(adj, color, target_k, rng):
    """Greedily grow a random independent set, alternating which
    bipartite class we draw from to keep it balanced, stopping early if
    we run out of eligible candidates on either side."""
    n = len(adj)
    order_by_class = {0: list(np.where(np.array(color) == 0)[0]), 1: list(np.where(np.array(color) == 1)[0])}
    for c in (0, 1):
        rng.shuffle(order_by_class[c])
    blocked = [False] * n
    chosen = {0: [], 1: []}
    idx = {0: 0, 1: 0}
    turn = 0
    stalled_class = set()
    while len(chosen[0]) + len(chosen[1]) < target_k and len(stalled_class) < 2:
        c = turn % 2
        turn += 1
        placed = False
        while idx[c] < len(order_by_class[c]):
            cand = order_by_class[c][idx[c]]
            idx[c] += 1
            if blocked[cand]:
                continue
            chosen[c].append(cand)
            blocked[cand] = True
            for v in adj[cand]:
                blocked[v] = True
            placed = True
            break
        if not placed:
            stalled_class.add(c)
    return chosen[0] + chosen[1], len(chosen[0]), len(chosen[1])


def perfect_matching_exists(adj, color, monomers, return_matching=False):
    """Max-flow check: does the graph induced on (all triangles minus
    monomers) admit a perfect matching? With return_matching=True, also
    extracts the actual matched pairs from the flow's saturated edges."""
    mono = set(monomers)
    left = [i for i, c in enumerate(color) if c == 0 and i not in mono]
    right = [i for i, c in enumerate(color) if c == 1 and i not in mono]
    left_pos = {v: i for i, v in enumerate(left)}
    right_pos = {v: i for i, v in enumerate(right)}
    n_left, n_right = len(left), len(right)
    if n_left != n_right:
        return (False, n_left, n_right, []) if return_matching else (False, n_left, n_right)
    n_nodes = 2 + n_left + n_right
    source, sink = 0, n_nodes - 1
    rows, cols, caps = [], [], []
    edge_pairs = []
    for v in left:
        rows.append(source); cols.append(1 + left_pos[v]); caps.append(1)
    for v in right:
        rows.append(1 + n_left + right_pos[v]); cols.append(sink); caps.append(1)
    for v in left:
        for u in adj[v]:
            if u in right_pos:
                rows.append(1 + left_pos[v]); cols.append(1 + n_left + right_pos[u]); caps.append(1)
                edge_pairs.append((v, u, 1 + left_pos[v], 1 + n_left + right_pos[u]))
    capacity = csr_matrix((caps, (rows, cols)), shape=(n_nodes, n_nodes))
    result = maximum_flow(capacity, source, sink)
    ok = result.flow_value == n_left
    if not return_matching:
        return ok, n_left, n_right
    matching = []
    if ok:
        flow = result.flow.tocsr()
        for v, u, iv, iu in edge_pairs:
            if flow[iv, iu] > 0:
                matching.append((v, u))
    return ok, n_left, n_right, matching


def main():
    lat = RhombileLattice(NX, NY)
    triangles = enumerate_all_triangles(lat)
    n_tri = len(triangles)
    adj = build_triangle_adjacency(triangles)
    color = bipartition(adj)
    n0, n1 = color.count(0), color.count(1)
    print(f"N_triangles={n_tri}  bipartite classes: {n0} / {n1}", flush=True)

    n_trials = 25
    rng = np.random.default_rng(SEED)
    for p in DENSITIES:
        target_k = int(round(p * n_tri))
        target_k -= target_k % 2  # keep it even (equal split target)
        successes = 0
        achieved_p = None
        for _ in range(n_trials):
            monomers, k0, k1 = greedy_balanced_independent_set(adj, color, target_k, rng)
            achieved_p = len(monomers) / n_tri
            ok, n_left, n_right = perfect_matching_exists(adj, color, monomers)
            successes += int(ok)
        print(f"target p={p:.2f}  achieved p={achieved_p:.4f}  monomers={target_k}  "
              f"success_rate={successes}/{n_trials}", flush=True)


if __name__ == "__main__":
    main()

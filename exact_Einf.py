"""Exact D3=infinity (pure 2-state AF Ising) ground-state energy for an
arbitrary bond pattern of the rhombile torus, with homology handling and
a correctness certificate.

Faces of the active-bond graph = triangles glued across inactive (J=0)
bonds; a face is frustrated (in T) iff it has an odd number of active
boundary bonds. Any spin configuration's violated-bond set V is a T-join
on the dual, so |V| >= min T-join (a lower bound, also on the torus).

On the torus V must additionally lie in the right Z2xZ2 homology class.
Shortest dual paths are therefore computed per class, by a 0-1 BFS on
(triangle, x-parity, y-parity): crossing an active bond costs 1,
crossing an inactive one (staying inside a face) costs 0, and a step
whose two hubs sit on opposite sides of the box seam flips the parity.

Candidates: the min-weight perfect matching with each pair on its
cheapest class (the lower bound), then, if that fails the certificate,
single-pair class changes in increasing cost order. The certificate
2-colours the lattice with "same" across V and "different" elsewhere;
success means E = |V| is achieved. Returns (lower, upper, certified).
For two defects only one class can certify, so the certified value is exact
even when it exceeds the unconstrained lower bound.
"""
from collections import deque
import itertools

import networkx as nx
import numpy as np

from monomer_dimer_matching_feasibility import enumerate_all_triangles

CLASSES = [(0, 0), (1, 0), (0, 1), (1, 1)]


class _Graph:
    def __init__(self, lat):
        self.lat = lat
        self.tri = enumerate_all_triangles(lat)
        n = len(self.tri)
        cell = {}
        for m in range(lat.ny):
            for k in range(lat.nx):
                cell[lat._site_index(k, m, "r1")] = (k, m)
        hub = [cell[t["sites"][0]] for t in self.tri]
        self.adj = [[] for _ in range(n)]
        by_bond = {}
        for ti, t in enumerate(self.tri):
            for key in ("r1r2", "r1r3", "r2r3"):
                by_bond.setdefault(id(t[key]), (t[key], []))[1].append(ti)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for bond, tis in by_bond.values():
            if len(tis) != 2:
                continue
            a, b = tis
            fx = int(abs(hub[a][0] - hub[b][0]) > 1)
            fy = int(abs(hub[a][1] - hub[b][1]) > 1)
            cost = 0 if bond["J"] == 0.0 else 1
            self.adj[a].append((b, bond, cost, fx, fy))
            self.adj[b].append((a, bond, cost, fx, fy))
            if cost == 0:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
        face = [find(t) for t in range(n)]
        deg = {}
        for bond, tis in by_bond.values():
            if bond["J"] != 0.0:
                for ti in tis:
                    deg[face[ti]] = deg.get(face[ti], 0) + 1
        self.T = sorted(f for f, d in deg.items() if d % 2 == 1)  # face reps = triangle ids

    def paths_from(self, s):
        n = len(self.tri)
        dist = np.full((n, 2, 2), 1 << 30, dtype=np.int64)
        prev = {}
        dist[s, 0, 0] = 0
        dq = deque([(s, 0, 0)])
        while dq:
            t, px, py = dq.popleft()
            d = dist[t, px, py]
            for u, bond, c, fx, fy in self.adj[t]:
                qx, qy = px ^ fx, py ^ fy
                if d + c < dist[u, qx, qy]:
                    dist[u, qx, qy] = d + c
                    prev[(u, qx, qy)] = ((t, px, py), bond, c)
                    (dq.appendleft if c == 0 else dq.append)((u, qx, qy))
        return dist, prev

    @staticmethod
    def path_bonds(prev, target_state):
        out, st = [], target_state
        while st in prev:
            st, bond, c = prev[st]
            if c:
                out.append(bond)
        return out


def _certify(lat, vset):
    color = [-1] * lat.n_sites
    nbr = lat.neighbor_table()
    for s0 in range(lat.n_sites):
        if color[s0] != -1:
            continue
        color[s0] = 0
        q = deque([s0])
        while q:
            x = q.popleft()
            for y, b in nbr[x]:
                if b["J"] == 0.0:
                    continue
                want = color[x] if id(b) in vset else 1 - color[x]
                if color[y] == -1:
                    color[y] = want
                    q.append(y)
                elif color[y] != want:
                    return False
    return True


def _xor(bond_lists):
    acc = {}
    for bl in bond_lists:
        for b in bl:
            if id(b) in acc:
                del acc[id(b)]
            else:
                acc[id(b)] = b
    return acc


def exact_two_state_energy(lat, max_fixups=400):
    G = _Graph(lat)
    T = G.T
    if not T:
        return 0, 0, _certify(lat, set())
    info = {s: G.paths_from(s) for s in T}
    cost = {}
    for s, t in itertools.combinations(T, 2):
        dist = info[s][0]
        cost[(s, t)] = {c: int(dist[t, c[0], c[1]]) for c in CLASSES}

    def bonds_of(s, t, c):
        return G.path_bonds(info[s][1], (t, c[0], c[1]))

    g = nx.Graph()
    for (s, t), cc in cost.items():
        g.add_edge(s, t, weight=min(cc.values()))
    matching = [tuple(sorted(e)) for e in nx.min_weight_matching(g)]
    best_c = {p: min(CLASSES, key=lambda c: cost[p][c]) for p in matching}
    lower = sum(cost[p][best_c[p]] for p in matching)

    def trial(assign):
        v = _xor(bonds_of(s, t, assign[(s, t)]) for (s, t) in matching)
        return len(v), _certify(lat, set(v))

    e, ok = trial(best_c)
    if ok:
        return lower, e, True
    options = sorted(((cost[p][c] - cost[p][best_c[p]], p, c)
                      for p in matching for c in CLASSES if c != best_c[p]))
    for extra, p, c in options[:max_fixups]:
        assign = dict(best_c)
        assign[p] = c
        e, ok = trial(assign)
        if ok:
            return lower, e, True
    return lower, None, False


def exact_two_state_energy_milp(lat, time_limit=600, knn=None, targets=None, bounds=False):
    """Exact E_inf for many defects: min-weight T-join with the global
    homology constraint, as a MILP (scipy/HiGHS). Variables x[pair, class]
    (each defect covered once); the XOR of the chosen classes must equal a
    target class. Any T-join in a class decomposes into paths pairing T plus
    cycles, and a cycle can be absorbed into one pair's walk class, so the
    MILP optimum for the certifying target is exact. All four targets are
    solved; exactly one of them certifies. Returns (E, per_target) where
    per_target maps class -> (milp optimum, realized |V|, certified).

    For many defects: knn=k keeps only pairs among each defect's k cheapest
    partners (the result is then a certified UPPER bound), targets restricts
    the classes solved, and bounds=True makes per_target entries
    (value, |V|, certified, dual_bound) with HiGHS's rigorous dual bound
    (a LOWER bound even when the time limit stops the search)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    G = _Graph(lat)
    T = G.T
    if not T:
        return 0, {}
    idx = {s: i for i, s in enumerate(T)}
    info = {s: G.paths_from(s) for s in T}
    keep = None
    if knn is not None:
        keep = set()
        for s in T:
            d = info[s][0]
            near = sorted((int(d[t].min()), t) for t in T if t != s)[:knn]
            keep.update(frozenset((s, t)) for _, t in near)
    var = []
    for s, t in itertools.combinations(T, 2):
        if keep is not None and frozenset((s, t)) not in keep:
            continue
        dist = info[s][0]
        for c in CLASSES:
            var.append((s, t, c, int(dist[t, c[0], c[1]])))
    nv = len(var)
    cost = np.array([v[3] for v in var] + [0, 0], dtype=float)  # + kx, ky
    A = lil_matrix((len(T) + 2, nv + 2))
    for j, (s, t, c, _) in enumerate(var):
        A[idx[s], j] = 1
        A[idx[t], j] = 1
        A[len(T), j] = c[0]
        A[len(T) + 1, j] = c[1]
    A[len(T), nv] = -2
    A[len(T) + 1, nv + 1] = -2
    A = A.tocsr()
    ub = np.concatenate([np.ones(nv), [len(T), len(T)]])
    out = {}
    for tgt in (targets or CLASSES):
        lo = np.concatenate([np.ones(len(T)), tgt])
        res = milp(cost, constraints=LinearConstraint(A, lo, lo),
                   integrality=np.ones(nv + 2), bounds=Bounds(0, ub),
                   options={"time_limit": time_limit})
        dual = getattr(res, "mip_dual_bound", None)
        if res.x is None:
            out[tgt] = (None, None, False) + ((dual,) if bounds else ())
            continue
        chosen = [var[j] for j in range(nv) if res.x[j] > 0.5]
        v = _xor(G.path_bonds(info[s][1], (t, c[0], c[1])) for s, t, c, _ in chosen)
        out[tgt] = (int(round(res.fun)), len(v), _certify(lat, set(v))) + ((dual,) if bounds else ())
    good = [r for r in out.values() if r[2]]
    return (min(r[0] for r in good) if good else None), out


def _validate():
    from rhombile_lattice import (RhombileLattice, route_string_between_points,
                                  apply_dual_string_defect, apply_string_defect)
    from dual_string_demo import build_dual
    from closed_loop_demo import hexagon_corners, NX, NY
    from case2_sa_check import build_case2

    print("pristine:", exact_two_state_energy(RhombileLattice(NX, NY)))
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    lat1, tri1, rho1, sib1, hop1 = build_dual(NX, NY)
    _, _, nodes, bonds = route_string_between_points(rho1, sib1, hop1, corners[i_left], corners[i_right])
    apply_dual_string_defect(tri1, nodes, bonds)
    print("Case 1:", exact_two_state_energy(lat1))
    lat2, *_ = build_case2()
    print("Case 2:", exact_two_state_energy(lat2))
    for dist in (5.0, 20.0):
        lat = RhombileLattice(30, 8)
        apply_string_defect(lat, np.array([2.0, 2.1]), np.array([2.0 + dist, 2.1]), wrap=False)
        print(f"Part II straight string dist={dist} (L={int(2 * dist)}):", exact_two_state_energy(lat))


if __name__ == "__main__":
    _validate()

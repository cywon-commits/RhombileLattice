"""Search for a planar, max-degree-3 CROSSOVER gadget for the D3 in (0,1) reduction.

Setting (see theorem_delta3.md).  The global lower bound is  cost >= D3 * (#packing triangles),
with equality iff the colouring is PERFECT: proper, colour 3 on exactly one vertex of every
packing triangle and nowhere else.  So a crossover gadget only has to be right on perfect
colourings:

  X is a planar graph, max degree 3, all triangles vertex-disjoint (they join the packing),
  four degree-1 terminals a, b, a2, b2 on the outer face in this cyclic order (never colour 3).
  Allowed boundary set  B(X) = {(a,b,a2,b2) in {1,2}^4 : some perfect colouring of X extends it}.
  Wanted:  B(X) = {a2 = a, b2 = b}  (4 patterns).

Abstract model used for speed: edges carry a parity p (p = 1: a real edge, colours differ;
p = 0: an edge subdivided once, colours equal).  A perfect colouring = a choice of one vertex
per triangle (removed, colour 3) such that the remaining parity system is consistent.
"""
import itertools
import random
import sys

import networkx as nx

A, B_, A2, B2 = 0, 1, 2, 3          # terminal ids
TARGET = {(x, y, x, y) for x in (0, 1) for y in (0, 1)}


class UF:
    def __init__(self, n):
        self.p = list(range(n)); self.par = [0] * n

    def find(self, x):
        if self.p[x] == x:
            return x, 0
        r, q = self.find(self.p[x])
        self.p[x] = r; self.par[x] ^= q
        return r, self.par[x]

    def union(self, x, y, d):          # colour(x) xor colour(y) == d
        (rx, px), (ry, py) = self.find(x), self.find(y)
        if rx == ry:
            return (px ^ py) == d
        self.p[rx] = ry; self.par[rx] = px ^ py ^ d
        return True


def boundary_set(n, edges, triangles, stop_if_bad=True):
    """edges: list of (u, v, parity). triangles: list of vertex triples. Returns B(X) as a set."""
    Bset = set()
    for choice in itertools.product(range(3), repeat=len(triangles)):
        removed = {t[c] for t, c in zip(triangles, choice)}
        uf = UF(n); ok = True
        for u, v, p in edges:
            if u in removed or v in removed:
                continue
            if not uf.union(u, v, p):
                ok = False; break
        if not ok:
            continue
        roots = [uf.find(t) for t in (A, B_, A2, B2)]
        for pat in itertools.product((0, 1), repeat=4):
            good = True
            seen = {}
            for (r, q), x in zip(roots, pat):
                if r in seen and seen[r] != (x ^ q):
                    good = False; break
                seen[r] = x ^ q
            if good:
                Bset.add(pat)
                if stop_if_bad and pat not in TARGET:
                    return Bset
    return Bset


def score(Bset):
    return len(Bset - TARGET) * 2 + len(TARGET - Bset)


def planar_with_order(n, edges):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from((u, v) for u, v, _ in edges)
    # force terminals onto one face in cyclic order a, b, a2, b2
    ring = [(A, B_), (B_, A2), (A2, B2), (B2, A)]
    G.add_edges_from(ring)
    G.add_edges_from((n, t) for t in (A, B_, A2, B2))
    return nx.check_planarity(G)[0]


def random_gadget(k, m, rng, p_even=0.35, tries=200):
    """4 terminals, k triangles, m plain vertices; random subcubic completion."""
    n = 4 + 3 * k + m
    tri = [tuple(range(4 + 3 * i, 7 + 3 * i)) for i in range(k)]
    plain = list(range(4 + 3 * k, n))
    cap = {t: 1 for t in range(4)}
    cap.update({v: 1 for t in tri for v in t})
    cap.update({v: 3 for v in plain})
    for _ in range(tries):
        edges = [(t[i], t[j], 1) for t in tri for i, j in ((0, 1), (1, 2), (0, 2))]
        free = dict(cap)
        adj = {v: set() for v in range(n)}
        for u, v, _ in edges:
            adj[u].add(v); adj[v].add(u)
        slots = [v for v in range(n) for _ in range(free[v])]
        rng.shuffle(slots)
        for _ in range(4 * n):
            cand = [v for v in range(n) if free[v] > 0]
            if len(cand) < 2:
                break
            u, v = rng.sample(cand, 2)
            if v in adj[u] or (u < 4 and v < 4):
                continue
            p = 0 if rng.random() < p_even else 1
            if p == 1 and adj[u] & adj[v]:
                continue                        # would create a new triangle
            edges.append((u, v, p)); adj[u].add(v); adj[v].add(u)
            free[u] -= 1; free[v] -= 1
        if any(free[t] for t in range(4)):
            continue
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((u, v) for u, v, _ in edges)
        if not nx.is_connected(G) or not planar_with_order(n, edges):
            continue
        return n, edges, tri
    return None


def mutate(n, edges, tri, rng, p_even=0.35):
    """remove one non-triangle edge and try to re-add up to two random edges."""
    tri_edges = {frozenset((t[i], t[j])) for t in tri for i, j in ((0, 1), (1, 2), (0, 2))}
    E = list(edges)
    loose = [e for e in E if frozenset(e[:2]) not in tri_edges]
    if not loose:
        return None
    E.remove(rng.choice(loose))
    cap = lambda v: 1 if v < 4 or v < 4 + 3 * len(tri) else 3
    for _ in range(rng.randint(1, 2)):
        deg = {v: 0 for v in range(n)}
        adj = {v: set() for v in range(n)}
        for u, v, _ in E:
            deg[u] += 1; deg[v] += 1; adj[u].add(v); adj[v].add(u)
        ext = lambda v: deg[v] - (2 if 4 <= v < 4 + 3 * len(tri) else 0)
        cand = [v for v in range(n) if ext(v) < cap(v)]
        if len(cand) < 2:
            break
        u, v = rng.sample(cand, 2)
        if v in adj[u] or (u < 4 and v < 4):
            continue
        p = 0 if rng.random() < p_even else 1
        if p == 1 and adj[u] & adj[v]:
            continue
        E.append((u, v, p))
    deg = {v: 0 for v in range(n)}
    for u, v, _ in E:
        deg[u] += 1; deg[v] += 1
    if any(deg[t] != 1 for t in range(4)):
        return None
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((u, v) for u, v, _ in E)
    if not nx.is_connected(G) or not planar_with_order(n, E):
        return None
    return E


def search(k, m, seed, restarts=50, steps=400):
    rng = random.Random(seed)
    best_overall = None
    for r in range(restarts):
        g = random_gadget(k, m, rng)
        if g is None:
            continue
        n, E, tri = g
        cur = score(boundary_set(n, E, tri, stop_if_bad=False))
        for s in range(steps):
            if cur == 0:
                return n, E, tri
            E2 = mutate(n, E, tri, rng)
            if E2 is None:
                continue
            sc = score(boundary_set(n, E2, tri, stop_if_bad=False))
            if sc <= cur:
                E, cur = E2, sc
        if best_overall is None or cur < best_overall[0]:
            best_overall = (cur, n, E, tri)
    return best_overall


if __name__ == "__main__":
    k, m, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    res = search(k, m, seed)
    print(k, m, seed, res if res and len(res) == 3 else ("best score", res[0] if res else None))

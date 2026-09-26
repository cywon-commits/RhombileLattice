"""Planar lattices with every interior site of degree 3 (non-bipartite through odd faces).

Each builder returns (n, edges, pos, info): n sites, edge list, 2D positions, and a dict with
'odd_faces' (number of odd interior faces) and 'interior' (sites of degree 3).

- delaunay_dual(points): sites = Delaunay triangles, bonds = shared triangle sides.  Faces are
  the Delaunay vertices; a face is odd iff that point has odd Delaunay degree.
    * random_cubic(N)       : random points (a Voronoi / 2D-foam-like network)
    * defect_honeycomb(L,p) : triangular-lattice points with a fraction p removed (5/7 defects)
- truncated(G, pos): replace every vertex of degree d by a d-cycle.  Odd faces = odd-degree
    vertices of G.  truncated_penrose(radius) applies it to a Penrose P3 patch.
"""
import math

import numpy as np
from scipy.spatial import Delaunay


def delaunay_dual(pts):
    tri = Delaunay(pts)
    T = tri.simplices
    n = len(T)
    edges = set()
    for i, nb in enumerate(tri.neighbors):
        for j in nb:
            if j >= 0:
                edges.add((min(i, int(j)), max(i, int(j))))
    pos = pts[T].mean(axis=1)
    # interior Delaunay vertices = faces of the dual
    hull = set(np.unique(tri.convex_hull))
    deg = {}
    for a, b in edges:
        pass
    vdeg = np.zeros(len(pts), int)
    nbrs = [set() for _ in pts]
    for t in T:
        for x in t:
            for y in t:
                if x != y:
                    nbrs[x].add(y)
    odd = sum(1 for v in range(len(pts)) if v not in hull and len(nbrs[v]) % 2 == 1)
    d = np.zeros(n, int)
    for a, b in edges:
        d[a] += 1; d[b] += 1
    return n, sorted(edges), pos, {"odd_faces": odd, "interior": [i for i in range(n) if d[i] == 3]}


def random_cubic(N, seed=0):
    rng = np.random.default_rng(seed)
    return delaunay_dual(rng.random((N, 2)))


def defect_honeycomb(L, p, seed=0, jitter=0.02):
    rng = np.random.default_rng(seed)
    pts = []
    for j in range(L):
        for i in range(L):
            if rng.random() >= p:
                pts.append((i + 0.5 * (j % 2), j * math.sqrt(3) / 2))
    pts = np.array(pts) + rng.normal(0, jitter, (len(pts), 2))
    return delaunay_dual(pts)


def truncated(Vpos, Eorig):
    """Vpos: {v: (x,y)}, Eorig: list of (u,v).  Returns the truncation as (n, edges, pos, info)."""
    adj = {v: [] for v in Vpos}
    for u, v in Eorig:
        adj[u].append(v); adj[v].append(u)
    idx = {}
    pos = []
    for v, nb in adj.items():
        x0, y0 = Vpos[v]
        nb.sort(key=lambda w: math.atan2(Vpos[w][1] - y0, Vpos[w][0] - x0))
        for w in nb:
            idx[(v, w)] = len(pos)
            x1, y1 = Vpos[w]
            pos.append((x0 + 0.3 * (x1 - x0), y0 + 0.3 * (y1 - y0)))
    edges = set()
    for v, nb in adj.items():
        d = len(nb)
        if d >= 2:
            for k in range(d):
                a, b = idx[(v, nb[k])], idx[(v, nb[(k + 1) % d])]
                if a != b and (d > 2 or k == 0):
                    edges.add((min(a, b), max(a, b)))
    for u, v in Eorig:
        a, b = idx[(u, v)], idx[(v, u)]
        edges.add((min(a, b), max(a, b)))
    n = len(pos)
    deg = np.zeros(n, int)
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    return n, sorted(edges), np.array(pos), {"interior": [i for i in range(n) if deg[i] == 3]}


def truncated_penrose(radius, seed=0):
    from penrose import tiling_graphs, default_gamma, E as EV
    R, vdeg, eo, dual = tiling_graphs(radius, default_gamma(seed))
    Vpos = {}
    Eorig = []
    for e in eo:
        u, v = tuple(e)
        for w in (u, v):
            Vpos[w] = tuple(np.array(w) @ EV)
        Eorig.append((u, v))
    # interior original vertices: all incident edges shared by two rhombi
    interior_v = {v for v in Vpos if all(len(o) == 2 for e, o in eo.items() if v in e)}
    n, edges, pos, info = truncated(Vpos, Eorig)
    info["odd_faces"] = sum(1 for v in interior_v if vdeg[v] % 2 == 1)
    return n, edges, pos, info

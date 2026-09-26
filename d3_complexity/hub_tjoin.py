"""Polynomial algorithm for the hub relaxation on Delaunay-dual (degree-3) lattices.

Lattice G: sites = Delaunay triangles, bonds = shared sides (lattices3.delaunay_dual).
Its planar dual is the Delaunay graph itself (faces of G = Delaunay points).  A 2-colouring's
same-colour bonds form a T-join in the Delaunay graph, T = odd faces.  A colour-3 site t
merges its three corner faces: model it as a hub node h_t joined to the three corners with
weight D3/2 (passing through the hub = paying D3 once).  Hull points are merged into one
outer node.  Minimum T-join = shortest paths between T-nodes + minimum-weight perfect matching.

The result equals the hub relaxation  min( mono among colours 1,2 + D3 * n3 ).
"""
import itertools

import networkx as nx
import numpy as np
from scipy.spatial import Delaunay


def hub_relaxation_value(pts, D3):
    tri = Delaunay(pts)
    hull = set(int(x) for x in np.unique(tri.convex_hull))
    OUT = "out"
    node = lambda p: OUT if p in hull else int(p)
    H = nx.Graph()

    def add(a, b, w):
        if a == b:
            return
        if H.has_edge(a, b):
            H[a][b]["weight"] = min(H[a][b]["weight"], w)
        else:
            H.add_edge(a, b, weight=w)

    deg = {}
    for t_idx, t in enumerate(tri.simplices):
        for a, b in itertools.combinations(t, 2):
            add(node(int(a)), node(int(b)), 1.0)
        h = ("hub", t_idx)
        for a in t:
            add(h, node(int(a)), D3 / 2)
    # odd faces: interior points of odd Delaunay degree
    nbrs = {}
    for t in tri.simplices:
        for a in t:
            for b in t:
                if a != b:
                    nbrs.setdefault(int(a), set()).add(int(b))
    T = [p for p in nbrs if p not in hull and len(nbrs[p]) % 2 == 1]
    if len(T) % 2 == 1:
        T.append(OUT)
    # minimum T-join: metric closure on T + min-weight perfect matching
    dist = {s: nx.single_source_dijkstra_path_length(H, s) for s in T}
    K = nx.Graph()
    for a, b in itertools.combinations(T, 2):
        K.add_edge(a, b, weight=-dist[a][b])
    M = nx.max_weight_matching(K, maxcardinality=True)
    return sum(dist[a][b] for a, b in M), len(T)

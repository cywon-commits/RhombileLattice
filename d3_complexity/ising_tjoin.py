"""Exact planar Ising ground state (antiferromagnetic, unit bonds) by minimum T-join in the dual.

Faces come from a planar embedding (networkx); T = odd faces (the outer face included, as it is
a real cycle of the patch).  Minimum T-join = shortest paths between T-faces + minimum-weight
perfect matching; the dual edges used are the frustrated bonds.  Returns (fr, colouring 1/2).
"""
import itertools

import networkx as nx


def faces_of(G):
    ok, emb = nx.check_planarity(G)
    assert ok
    seen = set()
    faces = []
    for u, v in emb.edges():
        if (u, v) in seen:
            continue
        faces.append(emb.traverse_face(u, v, mark_half_edges=seen))
    return faces


def ising_ground_state_tjoin(n, edges):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    faces = faces_of(G)
    side = {}
    for fi, f in enumerate(faces):
        for k in range(len(f)):
            a, b = f[k], f[(k + 1) % len(f)]
            side.setdefault(frozenset((a, b)), []).append(fi)
    between = {}
    Ds = nx.Graph()
    Ds.add_nodes_from(range(len(faces)))
    for e, fs in side.items():
        if len(fs) == 2 and fs[0] != fs[1]:
            Ds.add_edge(fs[0], fs[1])
            between.setdefault(frozenset(fs), e)
    T = [fi for fi, f in enumerate(faces) if len(f) % 2 == 1]
    paths = {s: nx.single_source_shortest_path(Ds, s) for s in T}
    K = nx.Graph()
    for a, b in itertools.combinations(T, 2):
        K.add_edge(a, b, weight=-(len(paths[a][b]) - 1))
    M = nx.max_weight_matching(K, maxcardinality=True)
    J = set()
    for a, b in M:
        p = paths[a][b]
        for x, y in zip(p, p[1:]):
            J ^= {between[frozenset((x, y))]}
    col = [0] * n
    for s in range(n):
        if col[s]:
            continue
        col[s] = 1
        stack = [s]
        while stack:
            u = stack.pop()
            for w in G[u]:
                want = col[u] if frozenset((u, w)) in J else 3 - col[u]
                if col[w] == 0:
                    col[w] = want
                    stack.append(w)
                else:
                    assert col[w] == want, "inconsistent T-join colouring"
    fr = sum(col[u] == col[v] for u, v in edges)
    return fr, col

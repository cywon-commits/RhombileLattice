"""Exhaustive check of Conjecture C.5 on all connected graphs with max degree <= 3 up to N vertices.

Generation: every connected graph has a non-cut vertex v, and G - v is connected and subcubic,
so all connected subcubic graphs on n+1 vertices arise from those on n vertices by adding a
vertex joined to 1-3 vertices of degree <= 2.  Isomorphism classes via WL hash + nx.is_isomorphic.
Check: some maximum cut admits an independent transversal of its uncut matching (ioct = fr).
"""
import itertools
import sys
from collections import defaultdict

import networkx as nx

from heal_moves import all_max_cuts, sat_for_cut


def extend(graphs):
    out = defaultdict(list)
    for G in graphs:
        n = G.number_of_nodes()
        free = [v for v in G if G.degree(v) <= 2]
        for k in (1, 2, 3):
            for nb in itertools.combinations(free, k):
                H = G.copy(); H.add_node(n)
                H.add_edges_from((n, v) for v in nb)
                h = nx.weisfeiler_lehman_graph_hash(H, iterations=4)
                if any(nx.is_isomorphic(H, K) for K in out[h]):
                    continue
                out[h].append(H)
    return [G for L in out.values() for G in L]


def healable(G):
    n, E = G.number_of_nodes(), list(G.edges())
    if n == 4 and len(E) == 6:
        return None            # K4
    fr, cuts = all_max_cuts(n, E)
    return any(sat_for_cut(n, E, t) for t in cuts)


if __name__ == "__main__":
    N = int(sys.argv[1])
    level = [nx.path_graph(2)]
    for n in range(3, N + 1):
        level = extend(level)
        bad = [G for G in level if healable(G) is False]
        print(f"n={n}: connected subcubic graphs {len(level)}, non-healable (excluding K4): {len(bad)}", flush=True)
        for G in bad[:3]:
            print("   ", sorted(G.edges()))

"""Exploration toward Conjecture C.5: which local moves always remove conflicts?

State: a maximum cut t (2-colouring with fr(G) same-colour edges; they form a matching M) and a
transversal S (one endpoint of every edge of M).  Conflicts kappa = #edges inside S.

Question: from EVERY (max cut, transversal) pair, can conflicts be driven to 0 by
  (a) endpoint swaps (same cut, other endpoint of one M-edge), and
  (b) moving to another maximum cut (any other max cut, with any transversal)?
Trivially (b) makes the state space "all max cuts x all transversals"; the conjecture is that
some state has kappa = 0.  Here we measure how far plain endpoint swaps go on each max cut,
i.e. for each max cut: is its 2-SAT satisfiable?  And we record the graphs/cuts where it is not.
"""
import itertools

import networkx as nx


def all_max_cuts(n, E):
    best, cuts = None, []
    for bits in range(1 << (n - 1)):          # vertex n-1 fixed to side 0
        t = [(bits >> v) & 1 for v in range(n - 1)] + [0]
        m = sum(t[u] == t[v] for u, v in E)
        if best is None or m < best:
            best, cuts = m, [t]
        elif m == best:
            cuts.append(t)
    return best, cuts


def sat_for_cut(n, E, t):
    adj = [set() for _ in range(n)]
    for u, v in E:
        adj[u].add(v); adj[v].add(u)
    M = [(u, v) for u, v in E if t[u] == t[v]]
    for choice in itertools.product((0, 1), repeat=len(M)):
        S = {e[c] for e, c in zip(M, choice)}
        if all(not (adj[a] & S) for a in S):
            return True
    return False


def analyse(G):
    G = nx.convert_node_labels_to_integers(G)
    n, E = G.number_of_nodes(), list(G.edges())
    fr, cuts = all_max_cuts(n, E)
    good = [t for t in cuts if sat_for_cut(n, E, t)]
    return fr, len(cuts), len(good)

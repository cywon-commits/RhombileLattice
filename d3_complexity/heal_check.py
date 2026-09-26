"""Healability certificate for degree-3 lattices (theorem_C.md, C.2 (iii)).

Given one Ising ground state (from any exact solver; here the MILP with colour 3 priced out),
its frustrated bonds form a matching.  Choosing one endpoint per frustrated bond so that the
chosen sites are pairwise non-adjacent is a 2-SAT problem.  If it is satisfiable, the Potts
ground state for every 0 < D3 < 1 is: Ising ground state with the chosen sites set to state 3,
and E(D3) = D3 * fr(G) exactly.  (If it fails, another Ising ground state may still work.)
"""
import networkx as nx

from potts_exact import solve
from ising_tjoin import ising_ground_state_tjoin


def ising_ground_state(n, edges):
    m, k, col = solve(n, edges, 10.0)       # D3 = 10 >= floor(3/2): colour 3 never used
    assert k == 0
    return m, col


def two_sat_heal(n, edges, col):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v); adj[v].add(u)
    M = [(u, v) for u, v in edges if col[u] == col[v]]
    ends = [v for e in M for v in e]
    assert len(ends) == len(set(ends)), "frustrated bonds do not form a matching"
    # literal (i, 0): choose M[i][0]; (i, 1): choose M[i][1]
    owner = {}
    for i, (a, b) in enumerate(M):
        owner[a] = (i, 0); owner[b] = (i, 1)
    G = nx.DiGraph()
    for i in range(len(M)):
        G.add_node((i, 0)); G.add_node((i, 1))
    neg = lambda L: (L[0], 1 - L[1])
    for p, Lp in owner.items():
        for q in adj[p]:
            if q in owner and owner[q][0] != Lp[0]:
                Lq = owner[q]
                # not (Lp and Lq):  Lp -> not Lq,  Lq -> not Lp
                G.add_edge(Lp, neg(Lq)); G.add_edge(Lq, neg(Lp))
    comp = {}
    for ci, C in enumerate(nx.strongly_connected_components(G)):
        for x in C:
            comp[x] = ci
    if any(comp[(i, 0)] == comp[(i, 1)] for i in range(len(M))):
        return None, M
    # standard assignment: condensation topological order
    C = nx.condensation(G, scc=list(nx.strongly_connected_components(G)))
    order = {c: k for k, c in enumerate(nx.topological_sort(C))}
    cmap = C.graph["mapping"]
    chosen = []
    for i, (a, b) in enumerate(M):
        pick0 = order[cmap[(i, 0)]] > order[cmap[(i, 1)]]
        chosen.append(a if pick0 else b)
    return chosen, M


def healed_state(n, edges, exact_milp=False):
    fr, col = ising_ground_state(n, edges) if exact_milp else ising_ground_state_tjoin(n, edges)
    chosen, M = two_sat_heal(n, edges, col)
    if chosen is None:
        return fr, None
    s = list(col)
    for v in chosen:
        s[v] = 3
    assert all(s[u] != s[v] for u, v in edges), "healed state is not proper"
    return fr, s

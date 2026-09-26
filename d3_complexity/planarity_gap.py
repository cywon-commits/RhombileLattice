"""Counterexample to the planarity claim for G(phi) in Johnson et al., arXiv:2211.12887v5, Thm 11.

phi below is a 2/3-CNF with every literal in <= 2 clauses and a planar variable-clause
incidence graph (even with each variable split into an x -- not-x edge).  In G(phi),
contract the path P to one vertex and each clause gadget {a,b,d,e,f,c} to its output c
(delete unused third inputs).  The result contains the literal-split incidence graph
plus one apex adjacent to every clause.  That minor is non-planar, so G(phi) is
non-planar for EVERY clause order and literal order.
"""
import networkx as nx

phi = [[5, 2], [4, -5], [1, -6], [6, 3], [3, -4], [6, 4, -2], [-2, 1], [-1, 5], [-6, -1]]
G = nx.Graph()
for v in range(1, 7):
    G.add_edge(("L", v), ("L", -v))
for j, cl in enumerate(phi):
    assert len({abs(l) for l in cl}) == len(cl)
    for l in cl:
        G.add_edge(("L", l), ("c", j))
occ = {}
for cl in phi:
    for l in cl:
        occ[l] = occ.get(l, 0) + 1
assert max(occ.values()) <= 2
H = G.copy()
for j in range(len(phi)):
    H.add_edge(("c", j), "apex")
print("literal-split incidence graph planar:", nx.check_planarity(G)[0])
print("minor of G(phi) (apex joined to all clauses) planar:", nx.check_planarity(H)[0])

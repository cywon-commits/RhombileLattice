"""Planar-CSP test for perfect colourings at Delta = 3 (see perfect_colouring_obstruction.md).

A perfect colouring (proper, colour 3 exactly once per packing triangle, nowhere else) of a
max-degree-3 graph whose triangles are vertex-disjoint is a Boolean CSP on the non-triangle
vertices: edges give the relation "!=", and every TRIANGLE CLUSTER (triangles joined by
triangle-triangle edges) gives a relation R_C on its outside neighbours.  R_C is always
self-complementary (colour 1<->2 symmetry).  By Dvorak-Kupec / Kazda-Kolmogorov-Rolinek
(arXiv:1602.03124, Thms 15-16) the planar CSP is polynomial iff every dR is an even
Delta-matroid, and NP-hard as soon as one dR is not (dT = cyclic differences of T, taken in
the cyclic order of the constraint's face).

This script enumerates small clusters with all terminals on one face, computes R_C and dR in
the face's cyclic order, and tests the even Delta-matroid exchange axiom.
"""
import itertools
import sys

import networkx as nx


def is_even_delta_matroid(M):
    M = set(M)
    if not M:
        return False
    if len({sum(x) % 2 for x in M}) != 1:
        return False
    for f in M:
        for g in M:
            diff = [i for i in range(len(f)) if f[i] != g[i]]
            for v in diff:
                ok = False
                for u in diff:
                    if u == v:
                        continue
                    h = list(f); h[u] ^= 1; h[v] ^= 1
                    if tuple(h) in M:
                        ok = True; break
                if not ok:
                    return False
    return True


def d_of(R):
    return {tuple((t[i] + t[(i + 1) % len(t)]) % 2 for i in range(len(t))) for t in R}


PERMS = list(itertools.permutations((1, 2, 3)))


def cluster_relation(j, links):
    """j triangles; triangle i has ports (i,0),(i,1),(i,2).  links: list of port pairs joined by
    an edge.  Every free port gets a pendant terminal.  Returns (terminal ports, relation over
    colours in {0,1} of the terminals, graph)."""
    used = {p for l in links for p in l}
    terms = [(i, k) for i in range(j) for k in range(3) if (i, k) not in used]
    R = set()
    for state in itertools.product(PERMS, repeat=j):
        col = {(i, k): state[i][k] for i in range(j) for k in range(3)}
        if any(col[a] == col[b] for a, b in links):
            continue
        # terminal t adjacent to port p: t in {1,2}, t != col[p]
        choices = []
        for p in terms:
            choices.append([c for c in (1, 2) if c != col[p]])
        for tc in itertools.product(*choices):
            R.add(tuple(c - 1 for c in tc))
    return terms, R


def outer_order(j, links, terms):
    """Cyclic order of terminals if they can all lie on one face; else None."""
    G = nx.Graph()
    for i in range(j):
        for a, b in ((0, 1), (1, 2), (0, 2)):
            G.add_edge((i, a), (i, b))
    for a, b in links:
        G.add_edge(a, b)
    for p in terms:
        G.add_edge(p, ("t", p))
    H = G.copy()
    for p in terms:
        H.add_edge("apex", ("t", p))
    ok, emb = nx.check_planarity(H)
    if not ok:
        return None
    order = [nb[1] for nb in emb.neighbors_cw_order("apex")]
    return order


def enumerate_links(j, max_links):
    ports = [(i, k) for i in range(j) for k in range(3)]
    pairs = [(a, b) for a, b in itertools.combinations(ports, 2) if a[0] != b[0]]
    for m in range(j - 1, max_links + 1):
        for L in itertools.combinations(pairs, m):
            used = [p for l in L for p in l]
            if len(used) != len(set(used)):
                continue
            G = nx.Graph(); G.add_nodes_from(range(j)); G.add_edges_from((a[0], b[0]) for a, b in L)
            if not nx.is_connected(G):
                continue
            # no extra triangles: two triangles joined by two links at ports sharing... (fine: 4-cycles)
            yield list(L)


if __name__ == "__main__":
    jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    seen = 0; bad = []
    for j in range(1, jmax + 1):
        for L in enumerate_links(j, min(3 * j // 2, j + 2)):
            terms, R = cluster_relation(j, L)
            if not terms:
                continue
            order = outer_order(j, L, terms)
            if order is None:
                continue
            idx = [terms.index(p) for p in order]
            Ro = {tuple(t[i] for i in idx) for t in R}
            assert Ro == {tuple(1 - x for x in t) for t in Ro}   # self-complementary
            seen += 1
            if not is_even_delta_matroid(d_of(Ro)):
                bad.append((j, L, len(terms), sorted(Ro)))
                print("NOT even Delta-matroid:", j, L, "terminals", len(terms), "|R| =", len(Ro))
                sys.stdout.flush()
        print(f"j={j}: clusters tested so far {seen}, failures {len(bad)}")

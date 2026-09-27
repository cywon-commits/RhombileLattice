"""Choi-Nakajima-Rim (1989, Thm 2) planar Delta=4 construction, re-used for the D3 problem.

Construction (SRC TR 87-203, pp. 9-12).  Variable v with n occurrences: a cycle
b^{k1} b^{k2} b^{k3} b^{k4} (k = 1..n) of length 4n; triangle tops a^{kl} over the bottom edge
(b^{kl}, b^{k,l+1}) (indices cyclic).  b^{k1}, b^{k3} <-> v ; b^{k2}, b^{k4} <-> not v.
Occurrence number k of v in clause c: positive -> tops a^{k2}, a^{k3} (they share b^{k3});
negative -> tops a^{k1}, a^{k2} (they share b^{k2}... see note).  Clause c: vertices w1..w4, edge
w1w4, and literal l = 1,2,3 joins w_l to the first top and w_{l+1} to the second top, so the
clause cycle w1 ->(lit 1)-> w2 ->(lit 2)-> w3 ->(lit 3)-> w4 -> w1 has length 13.

Note on the TR: its text pairs "not v" with tops a^{k1}, a^{k2}.  With tops over (b^{kl}, b^{k,l+1})
those share b^{k2} (a "not v" vertex), and a^{k2}, a^{k3} share b^{k3} (a "v" vertex), matching
the TR's example clause cycle.  Deleting the literal's own b breaks the cycle iff the literal is
true when true literals' b-vertices are deleted.

Claim tested here: for every fixed 0 < D3 < 2,
    OPT_D3(G) <= 6m * D3   <=>   phi satisfiable     (exactly-3-literal clauses, 28m vertices)
Proof idea (to be written): the 4n triangles of a variable component are edge-disjoint and each
b lies in two of them; charging a state-3 b by D3/2 to each of its triangles, every triangle costs
>= min(D3/2, 1) = D3/2, with equality iff the b's alternate around the cycle and nothing else pays.
"""
import itertools
import random
import sys

import networkx as nx

from potts_exact import solve


def build(nv, clauses):
    occ = {v: [] for v in range(1, nv + 1)}
    for j, cl in enumerate(clauses):
        for pos, lit in enumerate(cl):
            occ[abs(lit)].append((j, pos, lit > 0))
    idx = {}
    E = set()

    def node(x):
        if x not in idx:
            idx[x] = len(idx)
        return idx[x]

    def edge(x, y):
        a, b = node(x), node(y)
        E.add((min(a, b), max(a, b)))

    tops = {}
    for v, L in occ.items():
        n = len(L)
        if n == 0:
            continue
        cyc = [("b", v, k, l) for k in range(n) for l in range(4)]
        for i in range(4 * n):
            edge(cyc[i], cyc[(i + 1) % (4 * n)])
            top = ("a", v, i)
            edge(top, cyc[i]); edge(top, cyc[(i + 1) % (4 * n)])
        for k, (j, pos, positive) in enumerate(L):
            # positive: tops over (b^{k2},b^{k3}) and (b^{k3},b^{k4}) -> indices 4k+1, 4k+2 (share b^{k3})
            # negative: tops over (b^{k1},b^{k2}) and (b^{k2},b^{k3}) -> indices 4k+0, 4k+1 (share b^{k2})
            t1, t2 = (4 * k + 1, 4 * k + 2) if positive else (4 * k, 4 * k + 1)
            tops[(j, pos)] = (("a", v, t1), ("a", v, t2))
    for j, cl in enumerate(clauses):
        w = [("w", j, i) for i in range(4)]
        edge(w[0], w[3])
        for pos in range(3):
            x, y = tops[(j, pos)]
            edge(w[pos], x); edge(w[pos + 1], y)
    return len(idx), sorted(E), idx


def sat(nv, clauses):
    return any(all(any((l > 0) == bool(a[abs(l) - 1]) for l in cl) for cl in clauses)
               for a in itertools.product((0, 1), repeat=nv))


if __name__ == "__main__":
    seed, count = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    bad = 0; stats = {True: 0, False: 0}
    for t in range(count):
        nv = rng.randint(3, 4); m = rng.randint(2, 4)
        clauses = [[v * rng.choice((1, -1)) for v in rng.sample(range(1, nv + 1), 3)] for _ in range(m)]
        n, E, idx = build(nv, clauses)
        G = nx.Graph(E)
        assert max(dict(G.degree()).values()) <= 4
        s = sat(nv, clauses); stats[s] += 1
        for D3 in (0.3, 1.0, 1.5, 1.9):
            mono, n3, col = solve(n, E, D3)
            ok = (mono + D3 * n3 <= 6 * m * D3 + 1e-7) == s
            if not ok:
                bad += 1
                print("MISMATCH", clauses, D3, mono, n3, 6 * m * D3, s, flush=True)
        print(t, "m", m, "n", n, "sat", s, flush=True)
    print("sat/unsat", stats, "mismatches", bad)

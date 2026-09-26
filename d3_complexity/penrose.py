"""Penrose P3 rhombus patches via de Bruijn's pentagrid, and their face-adjacency (dual) graph.

Vertices of the tiling are integer 5-tuples K (exact identification); a rhombus is
created at every intersection of grid lines (j, r) and (k, s).  gamma sums to 0
(the Penrose condition) and is chosen generically, so no three lines meet.
"""
import itertools
import math
import numpy as np

E = np.array([[math.cos(2 * math.pi * j / 5), math.sin(2 * math.pi * j / 5)] for j in range(5)])


def default_gamma(seed=0):
    rng = np.random.default_rng(seed)
    g = rng.uniform(-0.5, 0.5, 5)
    return g - g.mean()


def rhombi(radius, gamma=None):
    """Return list of rhombi (K0, K1, K2, K3 in cyclic order, (j, k), centre) within `radius`."""
    g = default_gamma() if gamma is None else np.asarray(gamma, float)
    rmax = int(math.ceil(radius)) + 2
    out = []
    for j, k in itertools.combinations(range(5), 2):
        A = np.array([E[j], E[k]])
        for r in range(-rmax, rmax + 1):
            for s in range(-rmax, rmax + 1):
                z = np.linalg.solve(A, [r - g[j], s - g[k]])
                K = [math.ceil(E[i] @ z + g[i]) for i in range(5)]
                corners = []
                for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    Kc = list(K)
                    Kc[j], Kc[k] = r + a, s + b
                    corners.append(tuple(Kc))
                centre = sum(np.array(c) @ E for c in corners) / 4
                if np.hypot(*centre) <= radius:
                    out.append((tuple(corners), (j, k), centre))
    return out


def tiling_graphs(radius, gamma=None):
    """Return (rhombus list, tiling vertex degree dict, dual adjacency list)."""
    R = rhombi(radius, gamma)
    edge_owner = {}
    vdeg = {}
    for idx, (cs, _, _) in enumerate(R):
        for a in range(4):
            e = frozenset((cs[a], cs[(a + 1) % 4]))
            edge_owner.setdefault(e, []).append(idx)
    for e in edge_owner:
        for v in e:
            vdeg[v] = vdeg.get(v, 0) + 1
    dual = [set() for _ in R]
    for e, owners in edge_owner.items():
        assert len(owners) <= 2, "edge shared by more than two rhombi: not a tiling"
        if len(owners) == 2:
            a, b = owners
            dual[a].add(b)
            dual[b].add(a)
    return R, vdeg, edge_owner, dual


def is_thick(jk):
    """Thick rhombus (72 deg) iff |j-k| mod 5 in {1,4}; thin (36 deg) otherwise."""
    d = (jk[1] - jk[0]) % 5
    return d in (1, 4)

"""Monomer-sliding ("worm") dynamics on the triangle dimer covering.

State: every triangle is either matched to exactly one neighbor (their
shared bond is the one inactive edge, J=0) or is a monomer (all 3 edges
active, i.e. a frustrated triangle). The pristine lattice (j23=0) is the
rhombus tiling: every triangle matched to its r2r3 sibling.

split(t): turn on t's matched bond -> t and its partner become an
adjacent monomer pair.

move(t, u): monomer t, neighbor u (itself matched to w). Pair t with u
through their shared edge and free w -- the monomer jumps t -> w. Every
other triangle keeps exactly one inactive edge, so the frustrated set is
always exactly the monomer set, whatever path the walkers take (crossing
their own or each other's trails included), because each move is
defined on the CURRENT covering rather than on the original tiling.
"""
from collections import deque

import numpy as np

from rhombile_lattice import RhombileLattice, frustrated_triangles
from monomer_dimer_matching_feasibility import enumerate_all_triangles


class DimerState:
    def __init__(self, nx, ny):
        self.lat = RhombileLattice(nx, ny)
        self.tri = enumerate_all_triangles(self.lat)
        n = len(self.tri)
        self.n = n
        self.nbrs = [[] for _ in range(n)]
        groups = {}
        for ti, t in enumerate(self.tri):
            for key in ("r1r2", "r1r3", "r2r3"):
                groups.setdefault(id(t[key]), []).append((ti, t[key]))
        for entries in groups.values():
            assert len(entries) == 2, "expected every bond shared by exactly 2 triangles"
            (a, b), (c, _) = entries
            self.nbrs[a].append((c, b))
            self.nbrs[c].append((a, b))
        self.partner = [None] * n
        for ti in range(n):
            for u, b in self.nbrs[ti]:
                if b["J"] == 0.0:
                    self.partner[ti] = u
        assert all(p is not None for p in self.partner)
        self.monomers = set()
        self.site_to_tris = {}
        for ti, t in enumerate(self.tri):
            for s in t["sites"]:
                self.site_to_tris.setdefault(s, []).append(ti)

    def bond_between(self, a, b):
        for u, bond in self.nbrs[a]:
            if u == b:
                return bond
        raise ValueError("not adjacent")

    def split(self, t):
        u = self.partner[t]
        self.bond_between(t, u)["J"] = -1.0
        self.partner[t] = self.partner[u] = None
        self.monomers.update((t, u))
        return t, u

    def move(self, t, u):
        w = self.partner[u]
        self.bond_between(t, u)["J"] = 0.0
        self.bond_between(u, w)["J"] = -1.0
        self.partner[t], self.partner[u], self.partner[w] = u, t, None
        self.monomers.discard(t)
        self.monomers.add(w)
        return w

    def site_neighbors(self, t):
        """Triangles sharing at least one corner site with t (t excluded)."""
        out = set()
        for s in self.tri[t]["sites"]:
            out.update(self.site_to_tris[s])
        out.discard(t)
        return out

    def hop_distance(self, a, b):
        seen = {a: 0}
        q = deque([a])
        while q:
            x = q.popleft()
            if x == b:
                return seen[x]
            for y, _ in self.nbrs[x]:
                if y not in seen:
                    seen[y] = seen[x] + 1
                    q.append(y)
        return None


def sanity_single_pair(nx=16, ny=16, n_steps=3000, seed=0):
    rng = np.random.default_rng(seed)
    st = DimerState(nx, ny)
    t0 = int(rng.integers(st.n))
    a, b = st.split(t0)
    walkers = {a: None, b: None}  # current monomer -> last pivot (forbidden next neighbor)
    print(f"split triangle {t0}: monomers {a},{b}  hop distance {st.hop_distance(a, b)}", flush=True)
    for step in range(1, n_steps + 1):
        t = list(walkers)[rng.integers(2)]
        forbidden = walkers[t]
        cands = [u for u, _ in st.nbrs[t] if u not in st.monomers and u != forbidden]
        if not cands:
            cands = [u for u, _ in st.nbrs[t] if u not in st.monomers]
        u = cands[rng.integers(len(cands))]
        w = st.move(t, u)
        del walkers[t]
        walkers[w] = u
        if step % 300 == 0:
            ft = frustrated_triangles(st.lat)
            ms = list(walkers)
            print(f"step {step:5d}  frustrated={len(ft)}  "
                  f"hop_distance={st.hop_distance(ms[0], ms[1])}", flush=True)


if __name__ == "__main__":
    sanity_single_pair()

"""(a) Two-ensemble comparison of the effective interaction between a
defect (monomer) pair, same observable F(r) = -ln[P(r)/g(r)] in both.

A1 "dimer": bonds fluctuate, spins ignored. Uniform ensemble of dimer
coverings of the (honeycomb) triangle graph with exactly two monomers,
sampled by symmetric worm moves (pick a monomer, pick one of its 3
neighbours uniformly, always accept unless it is the other monomer) --
detailed balance holds, so the stationary measure is uniform and
P(r) = g(r) * W(r), with g(r) the number of (class-0, class-1) triangle
pairs at hop distance r. Prediction (critical dimer model, Fisher-
Stephenson / Ciucu): W(r) ~ r^(-1/2), i.e. F(r) = 0.5 ln r + const.

A2 "potts": bonds fixed to the Part V straight-string family, D3=0 Potts
spins fluctuate. Defect B moves only along a fixed straight track (worm
moves by one hop, +-1, range-limited), accepted only if the newly active
bond is satisfied by the current coloring; spins updated by zero-energy
heat bath plus Kempe-chain (Wang-Swendsen-Kotecky) swaps for ergodicity.
Stationary measure is uniform over valid (L, coloring) pairs, so
P(L) ~ N_colorings(L) and F(L) = -ln P(L) = -S(L) + const -- the actual
entropy, not the windowed state-3 proxy of entropic_force_T0_*.py.
Prediction (Kotecky-Salas-Sokal long-range order -> smooth height
phase): linear F(L) (entropic string tension).

Usage:
  python3 ensemble_comparison_A.py dimer <NX> <n_steps> <seed>
  python3 ensemble_comparison_A.py potts <NX> <NY> <Lmax> <n_steps> <seed> [bias]
"""
import sys
from collections import deque

import numpy as np

from rhombile_lattice import (
    enumerate_triangles, build_rhombi, build_triangle_hop_graph,
    route_string_between_points, frustrated_triangles, total_energy,
)
from worm_monomer_walk import DimerState


# ---------------------------------------------------------------- A1 dimer

def all_pairs_hop(st):
    n = st.n
    dist = np.full((n, n), -1, dtype=np.int16)
    for s in range(n):
        d = dist[s]
        d[s] = 0
        q = deque([s])
        while q:
            x = q.popleft()
            for y, _ in st.nbrs[x]:
                if d[y] < 0:
                    d[y] = d[x] + 1
                    q.append(y)
    return dist


def run_dimer(nx, n_steps, seed):
    rng = np.random.default_rng(seed)
    st = DimerState(nx, nx)
    dist = all_pairs_hop(st)
    a, b = st.split(int(rng.integers(st.n)))
    # class of each triangle = parity of hop distance from triangle 0
    cls = dist[0] % 2
    ca, cb = cls[a], cls[b]
    g = np.bincount(dist[np.ix_(np.where(cls == ca)[0], np.where(cls == cb)[0])].ravel())
    hist = np.zeros(len(g), dtype=np.int64)
    mons = [a, b]
    burn = n_steps // 10
    for step in range(n_steps):
        i = int(rng.integers(2))
        t = mons[i]
        u = st.nbrs[t][int(rng.integers(3))][0]
        if u != mons[1 - i]:
            mons[i] = st.move(t, u)
        if step >= burn:
            hist[dist[mons[0], mons[1]]] += 1
    assert len(frustrated_triangles(st.lat)) == 2
    return g, hist


# ---------------------------------------------------------------- A2 potts

class PottsTrack:
    def __init__(self, nx, ny, lmax, seed, y_cross=None, bias=0.0):
        self.rng = np.random.default_rng(seed)
        # multicanonical weight w(L) = exp(bias*L); undone in the output as
        # F(L) = -ln count(L) + bias*L. bias ~ ln 2 flattens the histogram.
        self.bias = bias
        self.st = DimerState(nx, ny)
        lat = self.st.lat
        y = y_cross if y_cross is not None else ny * np.sqrt(3) / 4
        tri_int = enumerate_triangles(lat)
        rho, _ = build_rhombi(tri_int)
        sib, hop = build_triangle_hop_graph(tri_int)
        # a long straight track through the middle row, away from the seam
        p0, p1 = np.array([3.0, y]), np.array([3.0 + lmax + 2.0, y])
        ts, te, nodes, bonds = route_string_between_points(rho, sib, hop, p0, p1)
        key = {frozenset(t["sites"]): i for i, t in enumerate(self.st.tri)}
        self.anchor = key[frozenset(tri_int[ts]["sites"])]
        self.track = [key[frozenset(tri_int[t]["sites"])] for t in nodes]
        assert self.st.partner[self.anchor] == self.track[0]
        self.st.split(self.anchor)      # monomers: anchor (fixed) and track[0]
        self.L = 0                      # B sits at track[L]
        self.lmax = min(lmax, len(self.track) - 1)
        self.pivot = []                 # pivots used, for retraction
        self.nbr = lat.neighbor_table()
        pos, sub_of = lat.site_positions()
        self.states = np.where(sub_of == "r1", 0, 1).astype(int)
        self.states = self._valid_start()

    def _valid_start(self):
        # D3=0 ground state for the initial (length-0) string: greedy fix-up
        s = self.states.copy()
        for _ in range(50):
            bad = [b for b in self.st.lat.bonds if b["J"] != 0.0 and s[b["i"]] == s[b["j"]]]
            if not bad:
                return s
            for b in bad:
                for v in (b["i"], b["j"]):
                    used = {s[j] for j, bb in self.nbr[v] if bb["J"] != 0.0}
                    free = [c for c in range(3) if c not in used]
                    if free:
                        s[v] = free[0]
                        break
        raise RuntimeError("could not build a zero-energy start")

    def _ok(self, bond):
        return self.states[bond["i"]] != self.states[bond["j"]]

    def track_move(self):
        st = self.st
        if self.rng.random() < 0.5:           # extend
            if self.L >= self.lmax:
                return
            t, w = self.track[self.L], self.track[self.L + 1]
            u = st.partner[w]
            f = st.bond_between(u, w)          # becomes active
            if not self._ok(f):
                return
            st.move(t, u)
            self.pivot.append(u)
            self.L += 1
        else:                                  # retract
            if self.L == 0:
                return
            if self.bias > 0 and self.rng.random() >= np.exp(-self.bias):
                return
            t = self.track[self.L]
            u = self.pivot[-1]
            w = st.partner[u]
            f = st.bond_between(u, w)          # becomes active again
            if not self._ok(f):
                return
            st.move(t, u)
            self.pivot.pop()
            self.L -= 1

    def spin_heat_bath(self):
        v = int(self.rng.integers(len(self.states)))
        used = {self.states[j] for j, b in self.nbr[v] if b["J"] != 0.0}
        free = [c for c in range(3) if c not in used]
        self.states[v] = free[int(self.rng.integers(len(free)))]

    def kempe(self):
        s = self.states
        v = int(self.rng.integers(len(s)))
        a = s[v]
        b = (a + 1 + int(self.rng.integers(2))) % 3
        comp = [v]
        seen = {v}
        q = deque([v])
        while q:
            x = q.popleft()
            for y, bond in self.nbr[x]:
                if bond["J"] != 0.0 and y not in seen and s[y] in (a, b):
                    seen.add(y)
                    comp.append(y)
                    q.append(y)
        for x in comp:
            s[x] = b if s[x] == a else a

    def _front_sites(self):
        """Sites near each possible front position (triangles within 3 hops of
        track[L] and their neighbour sites). Chosen by L alone, so extra
        heat-bath updates there keep detailed balance (they never change L)."""
        cache = []
        for L in range(self.lmax + 1):
            tris = self.track[max(0, L - 3): L + 4] + [self.anchor]
            sites = set()
            for t in tris:
                for v in self.st.tri[t]["sites"]:
                    sites.add(v)
                    sites.update(j for j, _ in self.nbr[v])
            cache.append(np.array(sorted(sites)))
        return cache

    def _heat_bath_at(self, v):
        used = {self.states[j] for j, b in self.nbr[v] if b["J"] != 0.0}
        free = [c for c in range(3) if c not in used]
        self.states[v] = free[int(self.rng.integers(len(free)))]

    def run(self, n_steps, n_local=20, n_global=10, kempe_every=10):
        front = self._front_sites()
        hist = np.zeros(self.lmax + 1, dtype=np.int64)
        burn = n_steps // 10
        for step in range(n_steps):
            loc = front[self.L]
            for v in loc[self.rng.integers(len(loc), size=n_local)]:
                self._heat_bath_at(int(v))
            for _ in range(n_global):
                self.spin_heat_bath()
            if step % kempe_every == 0:
                self.kempe()
            self.track_move()
            if step >= burn:
                hist[self.L] += 1
        e = total_energy(self.st.lat, self.states, (0.0, 0.0, 0.0))
        assert abs(e) < 1e-9, f"left the E=0 manifold: E={e}"
        assert len(frustrated_triangles(self.st.lat)) == 2
        return hist


def main():
    mode = sys.argv[1]
    if mode == "dimer":
        nx, n_steps, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        g, hist = run_dimer(nx, n_steps, seed)
        print(f"# dimer NX={nx} steps={n_steps} seed={seed}")
        print("# r  g(r)  count  F=-ln(count/g)")
        for r in range(len(g)):
            if g[r] > 0 and hist[r] > 0:
                print(f"{r} {g[r]} {hist[r]} {-np.log(hist[r] / g[r]):.4f}")
    elif mode == "potts":
        nx, ny, lmax, n_steps, seed = map(int, sys.argv[2:7])
        bias = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
        pt = PottsTrack(nx, ny, lmax, seed, bias=bias)
        hist = pt.run(n_steps)
        print(f"# potts NX={nx} NY={ny} Lmax={pt.lmax} steps={n_steps} seed={seed} bias={bias}")
        print("# L  count  F=-ln(count)+bias*L")
        for L, c in enumerate(hist):
            if c > 0:
                print(f"{L} {c} {-np.log(c) + bias * L:.4f}")


if __name__ == "__main__":
    main()

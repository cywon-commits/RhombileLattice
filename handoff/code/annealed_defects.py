"""Mechanism test for Part IX: density of spin-defect triangles (all three
Ising spins equal) at D2 = 0, annealed model vs fixed-bond TAFM.

Energy units: an equal-spin active bond costs 1. In the TAFM a pair of
defect triangles costs 2 (flip a spin with two equal neighbours). In the
annealed model each triangle can inactivate only one bond (its dimer), so a
defect triangle keeps >= 2 uncovered equal bonds: a pair costs >= 4. If this
is the operative excitation, n_def ~ exp(-2/T) (annealed) vs exp(-1/T)
(TAFM), and xi ~ n_def^(-1/2) ~ exp(2/T) vs exp(1/T).
Also records the frustrated-triangle (monomer) density of the annealed model.

Usage: python3 annealed_defects.py <L> <T1,T2,...> <n_eq> <n_meas> <seed>
Writes results/tafm/defects_L<L>_s<seed>.txt
"""
import os
import sys

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, sweep, total_energy, _seed
from annealed_fss import count_monomers
from annealed_tafm import neighbours
from worm_monomer_walk import DimerState

D3 = 10.0


@njit(cache=True)
def n_defects(s, tri_sites):
    n = 0
    for t in range(tri_sites.shape[0]):
        a, b, c = tri_sites[t, 0], tri_sites[t, 1], tri_sites[t, 2]
        if s[a] == s[b] and s[b] == s[c]:
            n += 1
    return n


@njit(cache=True)
def tafm_sweeps(s, nbr, beta, n):
    ns = s.shape[0]
    for _ in range(n):
        for _ in range(ns):
            v = np.random.randint(ns)
            old = s[v]
            same = 0
            for q in range(6):
                if s[nbr[v, q]] == old:
                    same += 1
            dE = (6 - same) - same
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                s[v] = 1 - old


def main():
    L = int(sys.argv[1])
    Ts = sorted([float(x) for x in sys.argv[2].split(",")], reverse=True)
    n_eq, n_meas, seed = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    A = build_arrays(L)
    nbr = neighbours(A)
    tri = np.array([t["sites"] for t in DimerState(L, L).tri], np.int64)   # same order as build_arrays
    assert all(len(set(r)) == 3 for r in tri)
    _seed(seed)
    rng = np.random.default_rng(seed)
    s_a = rng.integers(2, size=A["ns"]).astype(np.int64)
    s_t = s_a.copy()
    dim, mate = A["dimer"].copy(), A["mate"].copy()
    E = total_energy(s_a, dim, A["bi"], A["bj"], 0.0, D3, 0.0, mate)
    nt = A["nt"]
    lines = ["# T  n_def/N_tri(annealed)  n_def/N_tri(TAFM)  n_mon/N_tri(annealed)"]
    for T in Ts:
        b = 1.0 / T
        for _ in range(n_eq):
            E = sweep(s_a, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], b, 0.0, D3, 0.0, E)
        tafm_sweeps(s_t, nbr, b, n_eq)
        da = dt = dm = 0.0
        for _ in range(n_meas // 10):
            for _ in range(10):
                E = sweep(s_a, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], b, 0.0, D3, 0.0, E)
            tafm_sweeps(s_t, nbr, b, 10)
            da += n_defects(s_a, tri)
            dt += n_defects(s_t, tri)
            dm += count_monomers(mate)
        k = n_meas // 10
        lines.append(f"{T:.4f} {da / k / nt:.3e} {dt / k / nt:.3e} {dm / k / nt:.3e}")
        print(lines[-1], flush=True)
    os.makedirs("results/tafm", exist_ok=True)
    open(f"results/tafm/defects_L{L}_s{seed}.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

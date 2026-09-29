"""Equilibration / hysteresis test for the Z3 crystallisation (0=D1<D2, D3=10,
mu=0): single-temperature Metropolis runs (no tempering) from two starts,
  ordered : the rhombile crystal ground state (hubs in state 2, rims in 1)
  random  : pristine tiling with random spins, quenched to T directly
recording E/N_tri, |psi|^2 and the frustrated fraction every `every` sweeps.
If both starts converge to the same averages the long-time state is the
equilibrium one (slow dynamics only); if each stays in its own basin for
the whole run, the transition has a free-energy barrier (metastability).

Usage: python3 annealed_hysteresis.py <L> <T> <start: ordered|random> <n_sweeps> <every> <seed> [D2]
Writes results/fss/hyst_L<L>_T<T>_<start>_s<seed>.npz
"""
import os
import sys

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, sweep, total_energy, _seed
from annealed_fss import psi_components, count_monomers

D3, MU = 10.0, 0.0


@njit(cache=True)
def run(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, sub, beta, D2, E, n_sweeps, every, oE, oP, oN):
    rec = 0
    for it in range(n_sweeps):
        E = sweep(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, MU, E)
        if it % every == 0:
            re, im = psi_components(dimer, site_b, sub)
            oE[rec] = E
            oP[rec] = re * re + im * im
            oN[rec] = count_monomers(mate)
            rec += 1
    return E


def main():
    L, T, start = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    n_sweeps, every, seed = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    D2 = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    A = build_arrays(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    dimer, mate = A["dimer"].copy(), A["mate"].copy()          # pristine rhombile tiling, hubs on r1
    if start == "ordered":
        s = np.where(A["sub"] == 0, 1, 0).astype(np.int64)     # hubs -> state 2 (index 1), rims -> state 1
    else:
        s = rng.integers(2, size=A["ns"]).astype(np.int64)
    E = total_energy(s, dimer, A["bi"], A["bj"], D2, D3, MU, mate)
    n = n_sweeps // every + 1
    oE, oP, oN = np.zeros(n), np.zeros(n), np.zeros(n)
    run(s, dimer, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], A["sub"],
        1.0 / T, D2, E, n_sweeps, every, oE, oP, oN)
    os.makedirs("results/fss", exist_ok=True)
    fn = f"results/fss/hyst_L{L}_T{T:g}_{start}_s{seed}.npz"
    np.savez_compressed(fn, L=L, T=T, start=start, every=every, N_tri=A["nt"], E=oE, p2=oP, nm=oN)
    q = len(oP) // 4
    print(f"L={L} T={T} {start}: <psi2> quarters " + " ".join(f"{oP[i*q:(i+1)*q].mean():.4f}" for i in range(4)) +
          f"  E/Ntri last half {oE[len(oE)//2:].mean() / A['nt']:.5f}  -> {fn}", flush=True)


if __name__ == "__main__":
    main()

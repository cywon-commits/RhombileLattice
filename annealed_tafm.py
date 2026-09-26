"""Fixed-bond reference for Part IX: the triangular Ising antiferromagnet in a
field, in the same units as the annealed model.

All 9L^2 bonds of the triangular lattice are active (no dimers, no
frustrated-triangle freedom); spins s in {0, 1} (the D3 state is dropped):
  H = sum_{all bonds} delta(s_i, s_j) + D2 * n_1
i.e. Ising J = 1/2 and field h = D2/2, h/J = D2. For 0 < D2 < 6 the ground
state is the sqrt3 x sqrt3 state with the minority (state 1) on one of the
three triangular sublattices -- the same state as the rhombile crystal of
the annealed model (hubs = minority sites). The comparison isolates what
annealing the bonds does to T_c and to the transition.

Order parameter (the natural analogue of the annealed psi): m_k = fraction of
minority (state-1) sites on sublattice k, psi = sum_k m_k w^k. Records are in
the annealed_long.py format (E, re, im, nm; nm = n_1) so the same WHAM
analysis applies (annealed_long_analysis.py with PREFIX=results/tafm/long_).

Usage: python3 annealed_tafm.py <L> <T> <n_total> <thin> <chunk> <seed> [D2]
Writes results/tafm/long_L<L>_d2<D2>_T<T>_s<seed>.npz (checkpointed).
"""
import os
import sys
import time

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, _seed


@njit(cache=True)
def energy(s, nbr, D2):
    e = 0.0
    for i in range(s.shape[0]):
        for q in range(6):
            j = nbr[i, q]
            if j > i and s[i] == s[j]:
                e += 1.0
        if s[i] == 1:
            e += D2
    return e


@njit(cache=True)
def run_chunk(s, nbr, sub, beta, D2, E, n, thin, it0, oE, ore, oim, onm, rec0):
    ns = s.shape[0]
    cnt = np.zeros(3)
    for i in range(ns):
        cnt[sub[i]] += 1
    rec = rec0
    # precomputed acceptance factors: dE in {-6..6} + D2 * {-1, 0, 1}
    for it in range(it0, it0 + n):
        for _ in range(ns):
            v = np.random.randint(ns)
            old = s[v]
            new = 1 - old
            same_old = 0
            for q in range(6):
                if s[nbr[v, q]] == old:
                    same_old += 1
            dE = (6 - same_old) - same_old + D2 * (new - old)
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                s[v] = new
                E += dE
        if it % thin == 0:
            m = np.zeros(3)
            n1 = 0
            for i in range(ns):
                if s[i] == 1:
                    m[sub[i]] += 1
                    n1 += 1
            for k in range(3):
                m[k] /= cnt[k]
            oE[rec] = E
            ore[rec] = m[0] - 0.5 * (m[1] + m[2])
            oim[rec] = 0.8660254037844386 * (m[1] - m[2])
            onm[rec] = n1
            rec += 1
    return E, rec


def neighbours(A):
    ns = A["ns"]
    nbr = np.full((ns, 6), -1, np.int64)
    for i in range(ns):
        for q in range(6):
            k = A["site_b"][i, q]
            nbr[i, q] = A["bj"][k] if A["bi"][k] == i else A["bi"][k]
    return nbr


def main():
    L, T, n_total, thin = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    chunk, seed = int(sys.argv[5]), int(sys.argv[6])
    D2 = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    tag = f"L{L}_d2{D2:g}_T{T:g}_s{seed}"
    os.makedirs("results/tafm", exist_ok=True)
    ck, fn = f"results/tafm/ckpt_{tag}.npz", f"results/tafm/long_{tag}.npz"
    if os.path.exists(fn):
        print("done already:", fn)
        return
    A = build_arrays(L)
    nbr = neighbours(A)
    nrec = (n_total + thin - 1) // thin
    oE, ore, oim, onm = (np.zeros(nrec) for _ in range(4))
    if os.path.exists(ck):
        d = np.load(ck)
        s, it, rec = d["s"].copy(), int(d["it"]), int(d["rec"])
        for a, k in ((oE, "E"), (ore, "re"), (oim, "im"), (onm, "nm")):
            a[:rec] = d[k][:rec]
        _seed(seed * 1000003 + it // chunk)
    else:
        s = np.where(A["sub"] == 0, 1, 0).astype(np.int64)      # ordered start: minority on sublattice 0
        it, rec = 0, 0
        _seed(seed)
    E = energy(s, nbr, D2)
    t0 = time.time()
    while it < n_total:
        n = min(chunk, n_total - it)
        E, rec = run_chunk(s, nbr, A["sub"], 1.0 / T, D2, E, n, thin, it, oE, ore, oim, onm, rec)
        it += n
        assert abs(energy(s, nbr, D2) - E) < 1e-6
        tmp = ck + ".tmp.npz"
        np.savez(tmp, s=s, it=it, rec=rec, E=oE[:rec], re=ore[:rec], im=oim[:rec], nm=onm[:rec])
        os.replace(tmp, ck)
    f32 = lambda a: a[:rec, None].astype(np.float32)
    np.savez_compressed(fn, betas=np.array([1.0 / T]), T=T, L=L, N_site=A["ns"], N_tri=A["nt"], D2=D2, mu=0.0,
                        thin=thin, start="ordered", E=f32(oE), re=f32(ore), im=f32(oim), nm=f32(onm))
    os.remove(ck)
    p2 = ore[:rec] ** 2 + oim[:rec] ** 2
    print(f"{tag}: {time.time() - t0:.0f}s <psi2>={p2[rec // 10:].mean():.4f} "
          f"E/Ntri={oE[rec // 10:rec].mean() / A['nt']:.5f} -> {fn}", flush=True)


if __name__ == "__main__":
    main()

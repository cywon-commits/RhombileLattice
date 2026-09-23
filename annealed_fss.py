"""Finite-size scaling of the Z3 crystallisation transition of the annealed
model (0=D1<D2, D3=10, mu=0): parallel tempering over a narrow temperature
window around T_c, time series of (E, psi) saved for multi-histogram
reweighting (annealed_fss_analysis.py).

psi = sum_k m_k w^k, m_k = hub-type fraction on triangular sublattice k
(the Z3 order parameter of the rhombile crystal: which sublattice hosts
the hubs). One process = one (L, seed); replicas at K temperatures, one
sweep each per step, then adjacent swap attempts (even/odd alternating).

Usage:
  python3 annealed_fss.py <L> <D2> <Tmin> <Tmax> <K> <n_eq> <n_meas> <thin> <seed>
Writes results/fss/fss_L<L>_d2<D2>_s<seed>.npz
"""
import os
import sys
import time

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, sweep, total_energy, _seed

D3, MU = 10.0, 0.0


@njit(cache=True)
def psi_components(dimer, site_b, sub):
    m = np.zeros(3)
    cnt = np.zeros(3)
    for i in range(sub.shape[0]):
        cnt[sub[i]] += 1
        full = True
        for q in range(6):
            if dimer[site_b[i, q]]:
                full = False
        if full:
            m[sub[i]] += 1
    for k in range(3):
        m[k] /= cnt[k]
    return m[0] - 0.5 * (m[1] + m[2]), 0.8660254037844386 * (m[1] - m[2])


@njit(cache=True)
def count_monomers(mate):
    n = 0
    for t in range(mate.shape[0]):
        if mate[t] < 0:
            n += 1
    return n


@njit(cache=True)
def run_pt(S, DIM, MATE, E, betas, bi, bj, bta, btb, site_b, hexb, sub, D2,
           n_eq, n_meas, thin, out_E, out_re, out_im, out_nm, acc_swap, rep_at):
    K = betas.shape[0]
    rec = 0
    for it in range(n_eq + n_meas):
        for k in range(K):
            r = rep_at[k]
            E[r] = sweep(S[r], DIM[r], MATE[r], bi, bj, bta, btb, site_b, hexb,
                         betas[k], D2, D3, MU, E[r])
        start = it % 2
        for k in range(start, K - 1, 2):
            a, b = rep_at[k], rep_at[k + 1]
            x = (betas[k] - betas[k + 1]) * (E[a] - E[b])
            if x >= 0 or np.random.random() < np.exp(x):
                rep_at[k], rep_at[k + 1] = b, a
                if it >= n_eq:
                    acc_swap[k] += 1
        if it >= n_eq and (it - n_eq) % thin == 0:
            for k in range(K):
                r = rep_at[k]
                out_E[rec, k] = E[r]
                re, im = psi_components(DIM[r], site_b, sub)
                out_re[rec, k] = re
                out_im[rec, k] = im
                out_nm[rec, k] = count_monomers(MATE[r])
            rec += 1


def main():
    L, D2, tmin, tmax, K = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    n_eq, n_meas, thin, seed = int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9])
    A = build_arrays(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    betas = np.linspace(1 / tmax, 1 / tmin, K)
    S = rng.integers(3, size=(K, A["ns"])).astype(np.int64)
    DIM = np.tile(A["dimer"], (K, 1))
    MATE = np.tile(A["mate"], (K, 1))
    E = np.array([total_energy(S[r], DIM[r], A["bi"], A["bj"], D2, D3, MU, MATE[r]) for r in range(K)])
    nrec = n_meas // thin
    out = {k: np.zeros((nrec, K), dtype=np.float32) for k in ("E", "re", "im", "nm")}
    acc = np.zeros(K - 1)
    rep_at = np.arange(K)
    t0 = time.time()
    # pre-anneal every replica from T=1 down to its own T (avoids quench traps)
    for r in range(K):
        for b in np.linspace(1.0, betas[r], 20):
            for _ in range(max(1, n_eq // 40)):
                E[r] = sweep(S[r], DIM[r], MATE[r], A["bi"], A["bj"], A["bta"], A["btb"],
                             A["site_b"], A["hexb"], b, D2, D3, MU, E[r])
    run_pt(S, DIM, MATE, E, betas, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"],
           A["sub"], D2, n_eq, n_meas, thin, out["E"], out["re"], out["im"], out["nm"], acc, rep_at)
    for r in range(K):
        e = total_energy(S[r], DIM[r], A["bi"], A["bj"], D2, D3, MU, MATE[r])
        assert abs(e - E[r]) < 1e-6
    os.makedirs("results/fss", exist_ok=True)
    fn = f"results/fss/fss_L{L}_d2{D2:g}_s{seed}.npz"
    np.savez_compressed(fn, betas=betas, L=L, N_site=A["ns"], N_tri=A["nt"], D2=D2,
                        swap_acc=acc / (n_meas / 2), **out)
    print(f"L={L} D2={D2} seed={seed} K={K} T=[{tmin},{tmax}] n_meas={n_meas} "
          f"swap acc min/mean={acc.min() / (n_meas / 2):.2f}/{acc.mean() / (n_meas / 2):.2f} "
          f"time={time.time() - t0:.0f}s -> {fn}", flush=True)


if __name__ == "__main__":
    main()

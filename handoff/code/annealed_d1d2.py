"""D1=D2 annealed model (D2=0, D3=10, mu=0): parallel tempering over a wide
temperature window, several seeds, time series saved for reweighting.

Besides E and the frustrated-triangle number n_mon, each record carries the
spin order parameter at the sqrt3 x sqrt3 wavevector K. At D2=0 the spins are
effectively Ising (state 2 costs D3=10); in a T=0 configuration every
triangle has exactly one equal-spin bond and it is the dimer, so the T=0
ensemble is the ground-state ensemble of the triangular Ising
antiferromagnet (TAFM; Wannier entropy 0.3231/site = 0.16153/triangle),
critical with <s s> ~ r^(-1/2) cos(K.r). We record
  phi(q) = (1/N) sum_i s_i exp(i q.r_i),   s_i = +-1,
at q = K and at q = K + dq1, K + dq2 (the two smallest torus steps), which
gives the second-moment correlation length
  xi = sqrt(S(K)/S(K+dq) - 1) / (2 sin(|dq|/2)),
and the Binder cumulant U = 1 - <|phi|^4>/(2<|phi|^2>^2) of the Z6 order
parameter phi(K). A finite-T transition would show as a crossing of U or of
xi/L; the TAFM-like alternative is xi finite at all T>0, growing
exponentially as T -> 0.

Usage: python3 annealed_d1d2.py <L> <Tmin> <Tmax> <K> <n_eq> <n_meas> <thin> <seed>
Writes results/d1d2/pt_L<L>_s<seed>.npz
"""
import os
import sys
import time

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, sweep, total_energy, _seed
from annealed_fss import count_monomers
from rhombile_lattice import A_MAT_INV

D2, D3, MU = 0.0, 10.0, 0.0


def phases(L):
    """exp(i q.r_i) at q = K, K+dq1, K+dq2 (r in units of A1, A2; K = 2pi(1,0))."""
    from rhombile_lattice import RhombileLattice
    pos, _ = RhombileLattice(L, L).site_positions()
    frac = (A_MAT_INV @ pos.T).T
    qs = 2 * np.pi * np.array([[1.0, 0.0], [1.0 + 1.0 / L, 0.0], [1.0, 1.0 / L]])
    return np.exp(1j * frac @ qs.T)          # (N_site, 3)


@njit(cache=True)
def phi3(s, ph):
    out = np.zeros(3, np.complex128)
    for i in range(s.shape[0]):
        x = 2.0 * s[i] - 1.0 if s[i] < 2 else 0.0
        for q in range(3):
            out[q] += x * ph[i, q]
    return out / s.shape[0]


@njit(cache=True)
def run_pt(S, DIM, MATE, E, betas, bi, bj, bta, btb, site_b, hexb, ph,
           n_eq, n_meas, thin, out_E, out_nm, out_phi, acc_swap, rep_at):
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
                out_nm[rec, k] = count_monomers(MATE[r])
                out_phi[rec, k] = phi3(S[r], ph)
            rec += 1


def main():
    L, tmin, tmax, K = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    n_eq, n_meas, thin, seed = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    A = build_arrays(L)
    ph = phases(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    betas = np.exp(np.linspace(np.log(1 / tmax), np.log(1 / tmin), K))
    S = rng.integers(2, size=(K, A["ns"])).astype(np.int64)
    DIM = np.tile(A["dimer"], (K, 1))
    MATE = np.tile(A["mate"], (K, 1))
    E = np.array([total_energy(S[r], DIM[r], A["bi"], A["bj"], D2, D3, MU, MATE[r]) for r in range(K)])
    nrec = n_meas // thin
    oE = np.zeros((nrec, K), np.float32)
    onm = np.zeros((nrec, K), np.float32)
    ophi = np.zeros((nrec, K, 3), np.complex64)
    acc = np.zeros(K - 1)
    rep_at = np.arange(K)
    t0 = time.time()
    run_pt(S, DIM, MATE, E, betas, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], ph,
           n_eq, n_meas, thin, oE, onm, ophi, acc, rep_at)
    for r in range(K):
        e = total_energy(S[r], DIM[r], A["bi"], A["bj"], D2, D3, MU, MATE[r])
        assert abs(e - E[r]) < 1e-6
    os.makedirs("results/d1d2", exist_ok=True)
    fn = f"results/d1d2/pt_L{L}_s{seed}.npz"
    np.savez_compressed(fn, betas=betas, L=L, N_site=A["ns"], N_tri=A["nt"], D2=D2, mu=MU,
                        swap_acc=acc / (n_meas / 2), E=oE, nm=onm, phi=ophi)
    print(f"L={L} seed={seed} K={K} T=[{tmin},{tmax}] n_meas={n_meas} "
          f"swap acc min/mean={acc.min() / (n_meas / 2):.2f}/{acc.mean() / (n_meas / 2):.2f} "
          f"time={time.time() - t0:.0f}s -> {fn}", flush=True)


if __name__ == "__main__":
    main()

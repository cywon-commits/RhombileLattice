"""Correlation length of the fixed-bond triangular Ising antiferromagnet in
zero field, measured exactly as for the D1=D2 annealed model
(annealed_d1d2.py): phi(q) = (1/N) sum_i s_i exp(i q.r_i) at q = K, K+dq1,
K+dq2 on the same lattice and torus, xi = sqrt(S(K)/S(K+dq) - 1)/(2 sin(pi/L)).
Energy units as in the annealed model: an equal-spin bond costs 1 (Ising J=1/2).

Usage: python3 annealed_tafm_xi.py <L> <T1,T2,...> <n_eq> <n_meas> <seed>
Prints T, E/N_tri, S(K), xi, U(phi) per temperature; writes results/tafm/xi_L<L>_s<seed>.txt
"""
import os
import sys

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, _seed
from annealed_tafm import neighbours, energy
from annealed_d1d2 import phases, phi3


@njit(cache=True)
def run(s, nbr, beta, n_eq, n_meas, ph, out):
    ns = s.shape[0]
    for it in range(n_eq + n_meas):
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
        if it >= n_eq:
            out[it - n_eq] = phi3(s, ph)


def main():
    L = int(sys.argv[1])
    Ts = [float(x) for x in sys.argv[2].split(",")]
    n_eq, n_meas, seed = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    A = build_arrays(L)
    nbr = neighbours(A)
    ph = phases(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    s = rng.integers(2, size=A["ns"]).astype(np.int64)
    N = A["ns"]
    lines = ["# T  E/N_tri  S(K)  xi  U"]
    for T in sorted(Ts, reverse=True):              # anneal from high T
        out = np.zeros((n_meas, 3), np.complex128)
        run(s, nbr, 1.0 / T, n_eq, n_meas, ph, out)
        a = np.abs(out) ** 2 * N
        s0, s1 = a[:, 0].mean(), a[:, 1:].mean()
        xi = np.sqrt(max(s0 / s1 - 1, 0)) / (2 * np.sin(np.pi / L))
        m2 = np.abs(out[:, 0]) ** 2
        U = 1 - (m2 ** 2).mean() / (2 * m2.mean() ** 2)
        lines.append(f"{T:.4f} {energy(s, nbr, 0.0) / A['nt']:.5f} {s0:.3f} {xi:.3f} {U:.4f}")
        print(lines[-1], flush=True)
    os.makedirs("results/tafm", exist_ok=True)
    open(f"results/tafm/xi_L{L}_s{seed}.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

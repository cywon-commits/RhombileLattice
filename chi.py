"""Uniform susceptibility and reduced-field response of the annealed-bond model
(pure Ising sweep, ising_sweep.sweep_ising) and of the fixed-bond TAFM.

  chi = -(1/N_site) d<n1>/dD2 = var(n1) / (N_site T),   m = <n1>/N_site

Modes (one result file per run in results/susceptibility/raw/, skipped if it
exists, so a job queue can be relaunched after a container restart):

  T  <model> <L> <neq> <nmeas> <seed> <mu> <D2> <start> T1 T2 ...
     model = annealed | tafm; start = gs (random exact ground state: random
     lozenge tiling from hexagon flips + its zero-energy colouring) or ord
     (rhombile crystal, minority on sublattice 0 = the exact ground state for
     D2 > 0). Measures n1 every 2 sweeps; records T*chi, m, and the excitation
     densities (annealed: frustrated triangles n_mon/N_tri; TAFM: defect
     triangles with three equal spins / N_tri).
  gs <L> <nsweep> <seed> <start> H1 H2 ...
     Ground-state manifold with reduced field H = D2/(2T): only zero-energy TAFM
     flips (spin with exactly 3 equal neighbours), accepted with
     min(1, exp(-2 H dn1)). T*chi = var(n1)/N, m, and the uniform structure
     factor S(q) = |sum_i (s_i - <s>) e^{iq.r_i}|^2 / N at q = 2 pi k/L
     (k = 0..KMAX) along the two reciprocal directions, averaged. Local flips
     conserve the tiling tilt (flux sector); both starts are in tilt 0.

Started from chi/chi.py of the chat session (annealed_tafm_paper.zip).
"""
import os
import sys
import time

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, total_energy, _seed
from ising_sweep import sweep_ising

OUT = "results/susceptibility/raw"
KMAX = 6


@njit(cache=True)
def gs_sweeps(s, nbr, H, n, every, phase_c, phase_s):
    """Zero-energy TAFM flips with field H; returns n1 every sweep and S(q) sums every `every`."""
    N = s.shape[0]
    out = np.empty(n)
    nk = phase_c.shape[0]
    sq = np.zeros(nk)
    nsq = 0
    acc_up = np.exp(-2.0 * H)
    for it in range(n):
        for _ in range(N):
            v = np.random.randint(N)
            eq = 0
            for q in range(6):
                eq += s[nbr[v, q]] == s[v]
            if eq == 3:
                if s[v] == 1 or np.random.random() < acc_up:     # 0 -> 1 raises n1 by 1
                    s[v] = 1 - s[v]
        n1 = 0
        for i in range(N):
            n1 += s[i]
        out[it] = n1
        if it % every == 0:
            mean = n1 / N
            for k in range(nk):
                re = 0.0
                im = 0.0
                for i in range(N):
                    x = s[i] - mean
                    re += x * phase_c[k, i]
                    im += x * phase_s[k, i]
                sq[k] += (re * re + im * im) / N
            nsq += 1
    return out, sq / max(nsq, 1)


@njit(cache=True)
def tafm_run(s, nbr, tri, beta, D2, n, every, out_n1, out_def):
    N = s.shape[0]
    k = 0
    for it in range(n):
        for _ in range(N):
            v = np.random.randint(N)
            eq = 0
            for q in range(6):
                eq += s[nbr[v, q]] == s[v]
            dE = 6 - 2 * eq + D2 * (1 - 2 * s[v])        # equal bonds after - before, + D2 dn1
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                s[v] = 1 - s[v]
        if out_n1.shape[0] > 0 and it % every == 0 and k < out_n1.shape[0]:
            n1 = 0
            for i in range(N):
                n1 += s[i]
            nd = 0
            for t in range(tri.shape[0]):
                if s[tri[t, 0]] == s[tri[t, 1]] and s[tri[t, 1]] == s[tri[t, 2]]:
                    nd += 1
            out_n1[k] = n1
            out_def[k] = nd
            k += 1


@njit(cache=True)
def annealed_run(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E, n, every, out_n1, out_nm):
    k = 0
    for it in range(n):
        E = sweep_ising(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E)
        if out_n1.shape[0] > 0 and it % every == 0 and k < out_n1.shape[0]:
            n1 = 0
            for i in range(s.shape[0]):
                n1 += s[i]
            nm = 0
            for t in range(mate.shape[0]):
                nm += mate[t] < 0
            out_n1[k] = n1
            out_nm[k] = nm
            k += 1
    return E


def setup(L):
    A = build_arrays(L)
    nbr = np.empty((A["ns"], 6), np.int64)
    for v in range(A["ns"]):
        for q in range(6):
            k = A["site_b"][v, q]
            nbr[v, q] = A["bj"][k] if A["bi"][k] == v else A["bi"][k]
    # site triples of the triangles (for TAFM defect triangles)
    tri = np.empty((A["nt"], 3), np.int64)
    for t in range(A["nt"]):
        sites = set()
        for k in A["tri_b"][t]:
            sites |= {A["bi"][k], A["bj"][k]}
        tri[t] = sorted(sites)
    return A, nbr, tri


def phases(L):
    """cos/sin(q.r) for q = 2 pi k/L along the two reciprocal directions (k = 0..KMAX)."""
    from rhombile_lattice import BASIS, SUB_INDEX, A1, A2
    M = np.array([A1, A2]).T
    frac = {SUB_INDEX[k]: np.linalg.solve(M, v) for k, v in BASIS.items()}
    ns = 3 * L * L
    idx = np.arange(ns)
    cell, sub = idx // 3, idx % 3
    n, m = cell % L, cell // L
    f = np.array([frac[x] for x in sub])
    u = n + f[:, 0]
    w = m + f[:, 1]
    cs, sn = [], []
    for coord in (u, w):
        for k in range(KMAX + 1):
            ph = 2 * np.pi * k * coord / L
            cs.append(np.cos(ph))
            sn.append(np.sin(ph))
    return np.array(cs), np.array(sn)


def exact_gs(A, seed, nflip=None):
    """Random lozenge tiling (hexagon flips from the rhombile tiling) + its zero-energy colouring."""
    rng = np.random.default_rng(seed)
    dim, mate = A["dimer"].copy(), A["mate"].copy()
    hexb = A["hexb"]
    for v in rng.integers(A["ns"], size=nflip or 200 * A["ns"]):
        d = dim[hexb[v]]
        if (d[0] and d[2] and d[4] and not (d[1] or d[3] or d[5])) or \
                (d[1] and d[3] and d[5] and not (d[0] or d[2] or d[4])):
            for k in hexb[v]:
                dim[k] = not dim[k]
                if dim[k]:
                    mate[A["bta"][k]] = k
                    mate[A["btb"][k]] = k
    s = -np.ones(A["ns"], np.int64)
    s[0] = int(rng.integers(2))
    stack = [0]
    while stack:
        v = stack.pop()
        for k in A["site_b"][v]:
            if dim[k]:
                continue
            j = A["bj"][k] if A["bi"][k] == v else A["bi"][k]
            if s[j] < 0:
                s[j] = 1 - s[v]
                stack.append(j)
    return s, dim, mate


def start_state(A, start, seed):
    if start == "gs":
        return exact_gs(A, seed)
    return np.where(A["sub"] == 0, 1, 0).astype(np.int64), A["dimer"].copy(), A["mate"].copy()


def stats(x, nb=10):
    """mean and variance with block-jackknife errors."""
    x = np.asarray(x, float)
    n = len(x) // nb * nb
    x = x[:n]
    bl = x.reshape(nb, -1)
    jm = np.array([np.delete(bl, b, 0).mean() for b in range(nb)])
    jv = np.array([np.delete(bl, b, 0).var() for b in range(nb)])
    jk = lambda v: float(np.sqrt((nb - 1) / nb * np.sum((v - v.mean()) ** 2)))
    return x.mean(), jk(jm), x.var(), jk(jv)


def main():
    mode = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    if mode == "gs":
        L, n, seed, start = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        A, nbr, _ = setup(L)
        N = A["ns"]
        pc, ps = phases(L)
        for H in [float(x) for x in sys.argv[6:]]:
            fn = f"{OUT}/gs_L{L}_H{H:g}_s{seed}_{start}.txt"
            if os.path.exists(fn):
                continue
            t0 = time.time()
            np.random.seed(seed)
            s, _, _ = start_state(A, start, seed)
            n1, sq = gs_sweeps(s, nbr, H, n, 10, pc, ps)
            n1 = n1[n // 5:]
            m, em, v, ev = stats(n1 / N)
            sqa = 0.5 * (sq[:KMAX + 1] + sq[KMAX + 1:])
            with open(fn, "w") as f:
                f.write(f"# gs L={L} N={N} H={H} seed={seed} start={start} sweeps={n} (first 20% dropped)\n")
                f.write(f"m {m:.7f} {em:.7f}\nTchi {N * v:.6f} {N * ev:.6f}\n")
                f.write("Sq " + " ".join(f"{x:.6f}" for x in sqa) + "\n")
                f.write(f"time {time.time() - t0:.0f}\n")
            print(f"gs L={L} H={H:g} s={seed} {start}: m={m:.5f}({em:.5f}) Tchi={N * v:.4f}({N * ev:.4f}) "
                  f"S(q_k)={np.round(sqa, 4)} [{time.time() - t0:.0f}s]", flush=True)
        return
    model, L, neq, nmeas, seed = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    mu, D2, start = float(sys.argv[7]), float(sys.argv[8]), sys.argv[9]
    A, nbr, tri = setup(L)
    N, Nt = A["ns"], A["nt"]
    for T in [float(x) for x in sys.argv[10:]]:
        tag = f"{model}" + (f"_mu{mu:g}" if model == "annealed" else "") + f"_d2{D2:g}_L{L}_T{T:g}_s{seed}_{start}"
        fn = f"{OUT}/T_{tag}.txt"
        if os.path.exists(fn):
            continue
        t0 = time.time()
        rs = seed * 100003 + int(round(1000 * T))
        s, dim, mate = start_state(A, start, rs)
        nrec = nmeas // 2
        o1, o2 = np.zeros(nrec), np.zeros(nrec)
        if model == "tafm":
            np.random.seed(rs)
            tafm_run(s, nbr, tri, 1 / T, D2, neq, 2, np.zeros(0), np.zeros(0))
            tafm_run(s, nbr, tri, 1 / T, D2, nmeas, 2, o1, o2)
        else:
            _seed(rs)
            E = total_energy(s, dim, A["bi"], A["bj"], D2, 0.0, mu, mate)
            args = (s, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], 1 / T, D2, mu)
            E = annealed_run(*args, E, neq, 2, np.zeros(0), np.zeros(0))
            E = annealed_run(*args, E, nmeas, 2, o1, o2)
            assert abs(E - total_energy(s, dim, A["bi"], A["bj"], D2, 0.0, mu, mate)) < 1e-6
        m, em, v, ev = stats(o1 / N)
        x, ex, _, _ = stats(o2 / Nt)
        with open(fn, "w") as f:
            f.write(f"# {model} mu={mu} D2={D2} L={L} T={T} seed={seed} start={start} neq={neq} nmeas={nmeas}\n")
            f.write(f"m {m:.7f} {em:.7f}\nTchi {N * v:.6f} {N * ev:.6f}\n")
            f.write(f"{'nmon' if model == 'annealed' else 'defect'} {x:.7f} {ex:.7f}\n")
            f.write(f"time {time.time() - t0:.0f}\n")
        print(f"{tag}: m={m:.5f} Tchi={N * v:.4f}({N * ev:.4f}) exc={x:.5f} [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()

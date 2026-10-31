"""Uniform susceptibility chi = beta*var(n1)/N_site at D2=0 for the annealed-bond model and the
fixed-bond TAFM, and the ground-state-manifold constant chi0 = lim T*chi.
Usage: chi.py gs L nsweep | chi.py T model L neq nmeas T1 T2 ...  (model = annealed | tafm)"""
import sys, time
import numpy as np
from numba import njit
sys.path.insert(0, "/home/claude")
from annealed_thermo import build_arrays, total_energy, _seed
from ising_sweep import sweep_ising


@njit(cache=True)
def gs_sweeps(s, nbr, n, seed):
    """Zero-temperature TAFM dynamics: flip only spins with exactly 3 equal neighbours (dE=0).
    Uniform sampling of the TAFM ground-state manifold within a flux sector."""
    np.random.seed(seed)
    N = s.shape[0]; out = np.empty(n)
    for it in range(n):
        for _ in range(N):
            v = np.random.randint(N); eq = 0
            for q in range(6):
                eq += s[nbr[v, q]] == s[v]
            if eq == 3:
                s[v] = 1 - s[v]
        out[it] = s.sum()
    return out


@njit(cache=True)
def tafm_sweeps(s, nbr, beta, n, nmeas_every, seed):
    np.random.seed(seed)
    N = s.shape[0]; out = np.empty(n // nmeas_every)
    k = 0
    for it in range(n):
        for _ in range(N):
            v = np.random.randint(N); eq = 0
            for q in range(6):
                eq += s[nbr[v, q]] == s[v]
            dE = 6 - 2 * eq                      # equal bonds after flip minus before
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                s[v] = 1 - s[v]
        if it % nmeas_every == 0 and k < out.shape[0]:
            out[k] = s.sum(); k += 1
    return out


def setup(L):
    A = build_arrays(L)
    nbr = np.empty((A["ns"], 6), np.int64)
    for v in range(A["ns"]):
        for q in range(6):
            k = A["site_b"][v, q]
            nbr[v, q] = A["bj"][k] if A["bi"][k] == v else A["bi"][k]
    return A, nbr


def exact_gs(A, seed, nflip=None):
    """Random lozenge tiling (hexagon flips from the rhombile tiling) + its zero-energy colouring."""
    rng = np.random.default_rng(seed); dim, mate = A["dimer"].copy(), A["mate"].copy(); hexb = A["hexb"]
    for v in rng.integers(A["ns"], size=nflip or 200 * A["ns"]):
        d = dim[hexb[v]]
        if (d[0] and d[2] and d[4] and not (d[1] or d[3] or d[5])) or (d[1] and d[3] and d[5] and not (d[0] or d[2] or d[4])):
            for k in hexb[v]:
                dim[k] = not dim[k]
                if dim[k]:
                    mate[A["bta"][k]] = k; mate[A["btb"][k]] = k
    s = -np.ones(A["ns"], np.int64); s[0] = int(rng.integers(2)); stack = [0]
    while stack:
        v = stack.pop()
        for k in A["site_b"][v]:
            if dim[k]:
                continue
            j = A["bj"][k] if A["bi"][k] == v else A["bi"][k]
            if s[j] < 0:
                s[j] = 1 - s[v]; stack.append(j)
    return s, dim, mate


def chi_from(n1, N):
    m = n1 / N - 0.5                             # <m>=0 by spin-flip symmetry at D2=0
    blocks = np.array_split(m, 10)
    vals = [N * (b ** 2).mean() for b in blocks]
    return N * (m ** 2).mean(), np.std(vals) / np.sqrt(len(vals))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "gs":
        L, n = int(sys.argv[2]), int(sys.argv[3])
        A, nbr = setup(L); N = A["ns"]
        s = np.where(A["sub"] == 0, 1, 0).astype(np.int64)
        n1 = gs_sweeps(s, nbr, n, 1)
        tchi, err = chi_from(n1[n // 5:], N)
        print(f"GS manifold L={L} N={N}: T*chi -> chi0 = {tchi:.4f} ({err:.4f})", flush=True)
    else:
        model, L, neq, nmeas = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        A, nbr = setup(L); N = A["ns"]
        for T in [float(x) for x in sys.argv[6:]]:
            t0 = time.time(); seed = int(1000 * T) + 7
            s, dim, mate = exact_gs(A, seed)
            if model == "tafm":
                tafm_sweeps(s, nbr, 1 / T, neq, neq + 1, seed)
                n1 = tafm_sweeps(s, nbr, 1 / T, nmeas, 2, seed + 1)
                nmon = np.nan
            else:
                _seed(seed); np.random.seed(seed)
                E = total_energy(s, dim, A["bi"], A["bj"], 0.0, 0.0, 0.0, mate)
                for _ in range(neq):
                    E = sweep_ising(s, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], 1 / T, 0.0, 0.0, E)
                n1 = []; nm = []
                for k in range(nmeas):
                    E = sweep_ising(s, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], 1 / T, 0.0, 0.0, E)
                    if k % 2 == 0:
                        n1.append(s.sum()); nm.append((mate < 0).mean())
                n1 = np.array(n1, float); nmon = np.mean(nm)
            tchi, err = chi_from(n1, N)
            print(f"{model} L={L} T={T:.3f} Tchi={tchi:.4f} ({err:.4f}) nmon={nmon:.4f} [{time.time()-t0:.0f}s]", flush=True)

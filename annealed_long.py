"""Long single-temperature runs for the Z3 crystallisation (0=D1<D2, D3=10),
with checkpoint/resume, for equilibrated finite-size scaling at L >= 48.

One process = one (L, T, seed). Starts from the ordered rhombile crystal
(hubs on sublattice r1 in state index 1, rims in 0), runs Metropolis sweeps
(annealed_thermo.sweep) and records (E, Re psi, Im psi, n_mon) every `thin`
sweeps. Every `chunk` sweeps the state and the time series so far are
written to results/long/ckpt_*.npz, so a killed container resumes where it
stopped. When n_total sweeps are done the series goes to
results/long/long_L<L>_d2<D2>[_mu<mu>]_T<T>_s<seed>.npz, in the format of
annealed_fss.py (betas of length 1, arrays of shape (n_rec, 1)), so
annealed_fss_analysis.py can reweight single-T runs after discarding the
first n_discard sweeps (see annealed_long_analysis.py).

Usage: python3 annealed_long.py <L> <T> <n_total> <thin> <chunk> <seed> [D2] [mu] [start]
"""
import os
import sys
import time

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, sweep, total_energy, _seed
from annealed_fss import psi_components, count_monomers

D3 = 10.0


@njit(cache=True)
def run_chunk(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, sub, beta, D2, mu, E,
              n, thin, it0, oE, ore, oim, onm, rec0):
    rec = rec0
    for it in range(it0, it0 + n):
        E = sweep(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, mu, E)
        if it % thin == 0:
            re, im = psi_components(dimer, site_b, sub)
            oE[rec] = E
            ore[rec] = re
            oim[rec] = im
            onm[rec] = count_monomers(mate)
            rec += 1
    return E, rec


def main():
    L, T, n_total, thin = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    chunk, seed = int(sys.argv[5]), int(sys.argv[6])
    D2 = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    mu = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0
    start = sys.argv[9] if len(sys.argv) > 9 else "ordered"
    tag = f"L{L}_d2{D2:g}" + (f"_mu{mu:g}" if mu else "") + f"_T{T:g}_s{seed}" + ("" if start == "ordered" else f"_{start}")
    os.makedirs("results/long", exist_ok=True)
    ck, fn = f"results/long/ckpt_{tag}.npz", f"results/long/long_{tag}.npz"
    if os.path.exists(fn):
        print("done already:", fn)
        return
    A = build_arrays(L)
    nrec = (n_total + thin - 1) // thin
    if os.path.exists(ck):
        d = np.load(ck)
        s, dimer, mate = d["s"].copy(), d["dimer"].copy(), d["mate"].copy()
        it, rec = int(d["it"]), int(d["rec"])
        oE, ore, oim, onm = (np.zeros(nrec) for _ in range(4))
        for a, k in ((oE, "E"), (ore, "re"), (oim, "im"), (onm, "nm")):
            a[:rec] = d[k][:rec]
        _seed(seed * 1000003 + it // chunk)          # fresh stream per resumed chunk
    else:
        rng = np.random.default_rng(seed)
        dimer, mate = A["dimer"].copy(), A["mate"].copy()
        if start == "ordered":
            s = np.where(A["sub"] == 0, 1, 0).astype(np.int64)
        else:
            s = rng.integers(2, size=A["ns"]).astype(np.int64)
        it, rec = 0, 0
        oE, ore, oim, onm = (np.zeros(nrec) for _ in range(4))
        _seed(seed)
    E = total_energy(s, dimer, A["bi"], A["bj"], D2, D3, mu, mate)
    t0 = time.time()
    while it < n_total:
        n = min(chunk, n_total - it)
        E, rec = run_chunk(s, dimer, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"],
                           A["sub"], 1.0 / T, D2, mu, E, n, thin, it, oE, ore, oim, onm, rec)
        it += n
        e = total_energy(s, dimer, A["bi"], A["bj"], D2, D3, mu, mate)
        assert abs(e - E) < 1e-6, (e, E)
        tmp = ck + ".tmp.npz"
        np.savez(tmp, s=s, dimer=dimer, mate=mate, it=it, rec=rec,
                 E=oE[:rec], re=ore[:rec], im=oim[:rec], nm=onm[:rec])
        os.replace(tmp, ck)
        p2 = ore[:rec] ** 2 + oim[:rec] ** 2
        print(f"{tag}: {it}/{n_total} sweeps, {time.time() - t0:.0f}s, last-chunk <psi2>="
              f"{p2[-(n // thin):].mean():.4f} E/Ntri={oE[rec - 1] / A['nt']:.5f}", flush=True)
    f32 = lambda a: a[:rec, None].astype(np.float32)
    np.savez_compressed(fn, betas=np.array([1.0 / T]), T=T, L=L, N_site=A["ns"], N_tri=A["nt"], D2=D2, mu=mu,
                        thin=thin, start=start, E=f32(oE), re=f32(ore), im=f32(oim), nm=f32(onm))
    os.remove(ck)
    print("wrote", fn, flush=True)


if __name__ == "__main__":
    main()

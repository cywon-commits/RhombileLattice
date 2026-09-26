"""FSS of the Z3 crystallisation from long single-T runs (annealed_long.py).

For each (L, T) run: drop the first n_discard sweeps, report the
integrated autocorrelation time of |psi|^2 and E, the number of tau's in
the kept series, and a drift test (difference of the two halves in units
of its error). Then all T of one L are combined by the exact-level
multi-histogram reweighting of annealed_fss_analysis.py; errors come from a
jackknife over NB time blocks (block b of every run is dropped together),
which replaces the seed-group jackknife of the PT analysis.

Usage: [MU=<mu>] [TC=<T_c>] [DISCARD=100000] [NB=8] [CKPT=1] \
       python3 annealed_long_analysis.py <D2> L1 L2 ...
CKPT=1 also reads unfinished checkpoints (monitoring only).
"""
import glob
import os
import sys

import numpy as np

import annealed_fss_analysis as fa

DISCARD = int(os.environ.get("DISCARD", 100000))
NB = int(os.environ.get("NB", 8))
USE_CKPT = os.environ.get("CKPT") == "1"
MU = float(os.environ["MU"]) if os.environ.get("MU") else 0.0


def tau_int(x):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac /= ac[0]
    tau = 0.5
    for k in range(1, n):
        tau += ac[k]
        if k >= 6 * tau:
            break
    return tau


def files(D2, L):
    tag = f"L{L}_d2{D2:g}" + (f"_mu{MU:g}" if MU else "") + "_T"
    fs = sorted(glob.glob(f"results/long/long_{tag}*_s*.npz"))
    if USE_CKPT:
        done = {f.replace("long_", "ckpt_") for f in fs}
        fs += [f for f in sorted(glob.glob(f"results/long/ckpt_{tag}*_s*.npz")) if f not in done]
    return [f for f in fs if "_random" not in f and not f.endswith(".tmp.npz")]


def read(f):
    d = np.load(f)
    if "betas" in d:
        thin = int(d["thin"])
        out = dict(beta=float(d["betas"][0]), E=d["E"][:, 0].astype(float), re=d["re"][:, 0].astype(float),
                   im=d["im"][:, 0].astype(float), nm=d["nm"][:, 0].astype(float), thin=thin)
    else:                                   # checkpoint: thin from filename is unknown -> assume 10
        T = float(f.split("_T")[1].split("_s")[0])
        out = dict(beta=1 / T, E=d["E"].astype(float), re=d["re"].astype(float), im=d["im"].astype(float),
                   nm=d["nm"].astype(float), thin=10)
    k = DISCARD // out["thin"]
    for key in ("E", "re", "im", "nm"):
        out[key] = out[key][k:]
    return out


_cache = {}


def load_blocks(D2, L):
    """Runs in annealed_fss_analysis format, one entry per (T, time block);
    group = block index, so the jackknife drops block b of every T."""
    if (D2, L) in _cache:
        return _cache[(D2, L)]
    L_ = int(L)
    N, Nt = 3 * L_ * L_, 6 * L_ * L_
    runs = []
    for f in files(D2, L_):
        r = read(f)
        n = len(r["E"]) // NB * NB
        if n < NB * 10:
            continue
        for b in range(NB):
            sl = slice(b * n // NB, (b + 1) * n // NB)
            runs.append(dict(group=b, betas=np.array([r["beta"]]), E=r["E"][sl, None], re=r["re"][sl, None],
                             im=r["im"][sl, None], nm=r["nm"][sl, None], N=N, Nt=Nt))
    _cache[(D2, L)] = runs
    return runs


def diagnostics(D2, Ls):
    print(f"# per-run diagnostics (discard {DISCARD} sweeps): T, kept sweeps, tau(psi2), tau(E) [sweeps], "
          f"n_tau, <psi2>, drift(psi2) in sigma, <E>/N_tri")
    for L in Ls:
        for f in files(D2, L):
            r = read(f)
            p2 = r["re"] ** 2 + r["im"] ** 2
            if len(p2) < 100:
                continue
            tp, te = tau_int(p2) * r["thin"], tau_int(r["E"]) * r["thin"]
            h = len(p2) // 2
            a, b = p2[:h], p2[h:]
            ta = tp / r["thin"]
            err = np.sqrt((a.var() + b.var()) * 2 * ta / h)
            drift = (b.mean() - a.mean()) / err if err > 0 else 0
            kept = len(p2) * r["thin"]
            print(f"  L={L:3d} T={1 / r['beta']:.4f} kept={kept:8d} tau_psi={tp:7.0f} tau_E={te:6.0f} "
                  f"n_tau={kept / tp:6.1f} <psi2>={p2.mean():.4f} drift={drift:+5.1f} "
                  f"E={r['E'].mean() / (6 * L * L):.5f}{'  (ckpt)' if 'ckpt_' in f else ''}")


if __name__ == "__main__":
    D2 = float(sys.argv[1])
    Ls = [int(x) for x in sys.argv[2:]]
    diagnostics(D2, Ls)
    fa.load_L = load_blocks
    fa.main()

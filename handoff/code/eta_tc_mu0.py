"""eta at the mu=0 critical point: hub <|psi|^2>(L) at T = 0.522-0.528, L=32/48/64 by WHAM over the
fine-T long runs, L=96 from its single-T run at that T (CKPT=1 also reads unfinished L=96 checkpoints).
3-state Potts: eta = 4/15 = 0.267; Z3 locking: 4/9 = 0.444.
Usage: [CKPT=1] python3 results/long/eta_tc_mu0.py"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, ".")
import annealed_fss_analysis as fa  # noqa: E402
import annealed_long_analysis as la  # noqa: E402

Ts = [0.522, 0.524, 0.526, 0.528]
la.ORDER, la.MU, la.DIR = "hub", 0.0, "results/long"
m2 = {}
for L in (32, 48, 64):
    runs = [r for r in la.load_blocks(1.0, L) if 1 / r["betas"][0] >= 0.51]
    rw = fa.Reweighter(runs)
    # block jackknife over the NB time blocks
    groups = sorted({r["group"] for r in runs})
    full = [rw.moments(1 / T)["m2"] for T in Ts]
    jk = []
    for g in groups:
        rj = fa.Reweighter([r for r in runs if r["group"] != g])
        jk.append([rj.moments(1 / T)["m2"] for T in Ts])
    jk = np.array(jk)
    err = np.sqrt((len(jk) - 1) / len(jk) * ((jk - jk.mean(0)) ** 2).sum(0))
    m2[L] = (np.array(full), err)
pat = "ckpt" if os.environ.get("CKPT") == "1" else "long"
v, e = [], []
for T in Ts:
    fs = glob.glob(f"results/long/long_L96_d21_T{T:g}_s1.npz") or glob.glob(f"results/long/{pat}_L96_d21_T{T:g}_s1.npz")
    d = np.load(fs[0])
    re_, im_ = (d["re"], d["im"]) if d["re"].ndim == 1 else (d["re"][:, 0], d["im"][:, 0])
    p = re_[10000:].astype(float) ** 2 + im_[10000:].astype(float) ** 2
    v.append(p.mean())
    e.append(np.sqrt(p.var() * 2 * la.tau_int(p) / len(p)))
    print(f"L=96 T={T}: {os.path.basename(fs[0])}, {len(p) * 10} sweeps kept")
m2[96] = (np.array(v), np.array(e))
Ls = np.array(sorted(m2))
print("T      " + "  ".join(f"L={L:<14d}" for L in Ls) + "  eta(32-96)      eta(48-96)")
for i, T in enumerate(Ts):
    y = np.array([m2[L][0][i] for L in Ls])
    s = np.array([m2[L][1][i] for L in Ls]) / y
    out = []
    for sel in (Ls >= 32, Ls >= 48):
        X = np.vstack([np.ones(sel.sum()), np.log(Ls[sel])]).T
        W = np.diag(1 / s[sel] ** 2)
        C = np.linalg.inv(X.T @ W @ X)
        b = C @ X.T @ W @ np.log(y[sel])
        chi = ((np.log(y[sel]) - X @ b) ** 2 / s[sel] ** 2).sum() / max(1, sel.sum() - 2)
        out.append(f"{-b[1]:.3f}±{np.sqrt(C[1, 1] * max(1, chi)):.3f}")
    print(f"{T:.3f}  " + "  ".join(f"{a:.5f}({1e5 * b:4.0f})" for a, b in zip(y, y * s)) + "  " + "   ".join(out))

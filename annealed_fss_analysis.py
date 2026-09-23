"""Multi-histogram (Ferrenberg-Swendsen / WHAM) finite-size-scaling analysis
of annealed_fss.py output (results/fss/fss_L*_d2<D2>_s*.npz).

Per L, all seeds x all PT temperatures are combined into one reweighting
estimate as a continuous function of beta, giving:
  U4(T)        = 1 - <|psi|^4> / (2 <|psi|^2>^2)
  chi'(T)      = beta N (<|psi|^2> - <|psi|>^2)          (connected)
  dU4/dbeta    (numerical derivative of the reweighted U4)
  C(T)         = beta^2 (<E^2> - <E>^2) / N_tri
Fits (log-log in L): max dU/dbeta ~ L^(1/nu), max chi' ~ L^(gamma/nu),
<|psi|^2>(T_c) ~ L^(-2 beta/nu), C_max ~ L^(alpha/nu); T_c from Binder
crossings of successive sizes. Errors: leave-one-seed-out jackknife.

Usage: python3 annealed_fss_analysis.py <D2> L1 L2 ...
"""
import glob
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp

CANDIDATES = {  # 1/nu, gamma/nu, 2beta/nu, alpha/nu
    "2D 3-state Potts": (1.2, 26 / 15, 4 / 15, 0.4),
    "2D Ising": (1.0, 1.75, 0.25, 0.0),
    "first order": (2.0, 2.0, 0.0, 2.0),
}


def load_L(D2, L):
    files = sorted(glob.glob(f"results/fss/fss_L{L}_d2{D2:g}_s*.npz"))
    runs = []
    for f in files:
        d = np.load(f)
        runs.append(dict(betas=d["betas"], E=d["E"].astype(float), re=d["re"].astype(float),
                         im=d["im"].astype(float), N=int(d["N_site"]), Nt=int(d["N_tri"])))
    return runs


class Reweighter:
    """WHAM over samples: sample i from ensemble k (beta_k); all ensembles of
    all runs enter with their own sample counts."""

    def __init__(self, runs):
        E, P2, beta_of, n = [], [], [], []
        self.N, self.Nt = runs[0]["N"], runs[0]["Nt"]
        for r in runs:
            for k, b in enumerate(r["betas"]):
                E.append(r["E"][:, k])
                P2.append(r["re"][:, k] ** 2 + r["im"][:, k] ** 2)
                beta_of.append(b)
                n.append(r["E"].shape[0])
        self.E = np.concatenate(E)
        self.p2 = np.concatenate(P2)
        self.b = np.array(beta_of)
        self.n = np.array(n, float)
        self.E0 = self.E.mean()
        e = self.E - self.E0
        f = np.zeros(len(self.b))
        # log denominator for each sample: log sum_k n_k exp(f_k - b_k e_i)
        for _ in range(2000):
            logden = logsumexp(np.log(self.n)[:, None] + f[:, None] - np.outer(self.b, e), axis=0)
            fn = -logsumexp(-np.outer(self.b, e) - logden[None, :], axis=1)
            fn -= fn[0]
            if np.max(np.abs(fn - f)) < 1e-9:
                f = fn
                break
            f = fn
        self.logden = logsumexp(np.log(self.n)[:, None] + f[:, None] - np.outer(self.b, e), axis=0)
        self.e = e

    def moments(self, beta):
        lw = -beta * self.e - self.logden
        w = np.exp(lw - lw.max())
        w /= w.sum()
        p2, E = self.p2, self.E
        m2, m4, m1 = w @ p2, w @ p2 ** 2, w @ np.sqrt(p2)
        e1, e2 = w @ E, w @ E ** 2
        return dict(U=1 - m4 / (2 * m2 ** 2), chi=beta * self.N * (m2 - m1 ** 2), m2=m2,
                    C=beta ** 2 * (e2 - e1 ** 2) / self.Nt)


def curves(runs, bgrid):
    rw = Reweighter(runs)
    out = {k: np.array([rw.moments(b)[k] for b in bgrid]) for k in ("U", "chi", "m2", "C")}
    out["dU"] = np.gradient(out["U"], bgrid)
    out["N"] = rw.N
    return out


def observables(D2, Ls, bgrid, drop=None):
    """Curves per L; drop = seed index to leave out (jackknife)."""
    res = {}
    for L in Ls:
        runs = load_L(D2, L)
        if drop is not None:
            runs = [r for i, r in enumerate(runs) if i != drop]
        bmin = max(r["betas"].min() for r in runs)
        bmax = min(r["betas"].max() for r in runs)
        g = bgrid[(bgrid >= bmin) & (bgrid <= bmax)]
        res[L] = (g, curves(runs, g))
    return res


def crossing(g1, U1, g2, U2):
    g = np.intersect1d(np.round(g1, 10), np.round(g2, 10))
    a = np.interp(g, g1, U1)
    b = np.interp(g, g2, U2)
    d = b - a
    idx = np.where(np.sign(d[1:]) != np.sign(d[:-1]))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[np.argmin(np.abs(g[idx] - g.mean()))]      # crossing nearest window centre
    return 1 / (g[i] - d[i] * (g[i + 1] - g[i]) / (d[i + 1] - d[i]))


def fit_exponents(res, Ls, Tc):
    Ls = np.array(Ls, float)
    lg = np.log(Ls)
    dUmax = np.array([np.max(res[L][1]["dU"]) for L in Ls.astype(int)])
    chimax = np.array([np.max(res[L][1]["chi"]) for L in Ls.astype(int)])
    Cmax = np.array([np.max(res[L][1]["C"]) for L in Ls.astype(int)])
    m2c = np.array([np.interp(1 / Tc, res[L][0], res[L][1]["m2"]) for L in Ls.astype(int)])
    chiTc = np.array([np.interp(1 / Tc, res[L][0], res[L][1]["m2"]) * res[L][1]["N"] / Tc for L in Ls.astype(int)])
    dUTc = np.array([np.interp(1 / Tc, res[L][0], res[L][1]["dU"]) for L in Ls.astype(int)])
    Umin = np.array([np.min(res[L][1]["U"]) for L in Ls.astype(int)])
    return dict(inv_nu=np.polyfit(lg, np.log(dUmax), 1)[0],
                inv_nu_Tc=np.polyfit(lg, np.log(dUTc), 1)[0],
                gamma_nu_Tc=np.polyfit(lg, np.log(chiTc), 1)[0],
                chiTc=chiTc, dUTc=dUTc, Umin=Umin,
                gamma_nu=np.polyfit(lg, np.log(chimax), 1)[0],
                two_beta_nu=-np.polyfit(lg, np.log(m2c), 1)[0],
                alpha_nu=np.polyfit(lg, np.log(Cmax), 1)[0],
                Cmax=Cmax, dUmax=dUmax, chimax=chimax, m2c=m2c)


def main():
    D2 = float(sys.argv[1])
    Ls = [int(x) for x in sys.argv[2:]]
    nseed = min(len(load_L(D2, L)) for L in Ls)
    allb = np.concatenate([r["betas"] for L in Ls for r in load_L(D2, L)])
    bgrid = np.linspace(allb.min(), allb.max(), 400)

    def estimate(drop=None):
        res = observables(D2, Ls, bgrid, drop)
        tcs = [crossing(res[a][0], res[a][1]["U"], res[b][0], res[b][1]["U"])
               for a, b in zip(Ls[:-1], Ls[1:])]
        Tc = tcs[-1] if np.isfinite(tcs[-1]) else np.nanmean(tcs)
        return res, tcs, Tc, fit_exponents(res, Ls, Tc)

    res, tcs, Tc, ex = estimate()
    jk = [estimate(i) for i in range(nseed)] if nseed > 1 else []

    def err(get):
        if not jk:
            return float("nan")
        v = np.array([get(j) for j in jk], float)
        return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))

    print(f"D2={D2}  sizes {Ls}  seeds per size {nseed}")
    for (a, b), t, i in zip(zip(Ls[:-1], Ls[1:]), tcs, range(len(tcs))):
        print(f"  Binder crossing L={a}/{b}: T = {t:.4f} +- {err(lambda j: j[1][i]):.4f}")
    print(f"  T_c (largest pair) = {Tc:.4f} +- {err(lambda j: j[2]):.4f}")
    names = [("inv_nu", "1/nu (max dU)"), ("inv_nu_Tc", "1/nu (dU at T_c)"), ("gamma_nu", "gamma/nu (max chi')"),
             ("gamma_nu_Tc", "gamma/nu (chi at T_c)"), ("two_beta_nu", "2beta/nu"), ("alpha_nu", "alpha/nu")]
    for key, lab in names:
        print(f"  {lab:22s} = {ex[key]:.3f} +- {err(lambda j: j[3][key]):.3f}")
    print("  candidates (1/nu, gamma/nu, 2beta/nu, alpha/nu):")
    for k, v in CANDIDATES.items():
        print(f"    {k:18s} {v[0]:.3f} {v[1]:.3f} {v[2]:.3f} {v[3]:.3f}")
    print("  per L: max dU/dbeta, dU/dbeta(T_c), max chi', chi(T_c), <|psi|^2>(T_c), C_max, min U4")
    for i, L in enumerate(Ls):
        print(f"    L={L:3d}  {ex['dUmax'][i]:8.3f}  {ex['dUTc'][i]:8.3f}  {ex['chimax'][i]:9.3f}  {ex['chiTc'][i]:9.2f}"
              f"  {ex['m2c'][i]:.5f}  {ex['Cmax'][i]:.4f}  {ex['Umin'][i]:+.3f}")

    cols = ["#c5791f", "#3d7a4f", "#3b4ba8", "#a83b3b", "#6b4ba8"]
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for L, c in zip(Ls, cols):
        g, cv = res[L]
        T = 1 / g
        ax[0, 0].plot(T, cv["U"], color=c, label=f"L={L}")
        ax[0, 1].plot(T, cv["chi"] / L ** (26 / 15), color=c, label=f"L={L}")
        ax[0, 2].plot(T, cv["C"], color=c, label=f"L={L}")
        for r in load_L(D2, L)[:1]:
            p2 = r["re"] ** 2 + r["im"] ** 2
            ax[0, 0].plot(1 / r["betas"], 1 - np.mean(p2 ** 2, 0) / (2 * np.mean(p2, 0) ** 2), "o", ms=2.5, color=c)
    ax[0, 0].axvline(Tc, color="k", ls=":", lw=0.8)
    ax[0, 0].set_ylabel("Binder U4 (lines: reweighted, dots: raw, seed 1)")
    ax[0, 1].set_ylabel("chi' / L^(26/15)")
    ax[0, 2].set_ylabel("C / N_tri")
    for a in ax[0]:
        a.set_xlabel("T")
        a.legend(fontsize=7)
    Larr = np.array(Ls, float)
    panels = [(ax[1, 0], ex["dUmax"], "max dU4/dbeta", "inv_nu", 0), (ax[1, 1], ex["chimax"], "max chi'", "gamma_nu", 1),
              (ax[1, 2], ex["m2c"], "<|psi|^2> at T_c", "two_beta_nu", 2)]
    for a, y, lab, key, ci in panels:
        a.loglog(Larr, y, "o", color="k", label=f"data: slope {'-' if key == 'two_beta_nu' else ''}{ex[key]:.3f}")
        for (nm, v), c in zip(CANDIDATES.items(), ["#3b4ba8", "#c5791f", "#a83b3b"]):
            s = -v[ci] if key == "two_beta_nu" else v[ci]
            a.loglog(Larr, y[-1] * (Larr / Larr[-1]) ** s, "-", color=c, lw=1, label=f"{nm} ({s:+.3f})")
        a.set_xlabel("L")
        a.set_ylabel(lab)
        a.legend(fontsize=7)
    fig.suptitle(f"Z3 crystallisation of the annealed model, 0=D1<D2={D2:g}: multi-histogram FSS (T_c = {Tc:.4f})")
    fig.tight_layout()
    out = f"annealed_fss_d2{D2:g}_figure.png"
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()

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

Usage: [TC=<T_c>] [MU=<mu>] python3 annealed_fss_analysis.py <D2> L1 L2 ...
(TC fixes T_c, e.g. from pseudo_critical(); default: largest-pair Binder crossing)
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


MU_TAG = ""


def load_L(D2, L):
    files = sorted(glob.glob(f"results/fss/fss_L{L}_d2{D2:g}{MU_TAG}_s*.npz"))
    runs = []
    for f in files:
        d = np.load(f)
        seed = int(f.rsplit("_s", 1)[1].split(".")[0])
        runs.append(dict(group=seed % 10, betas=d["betas"], E=d["E"].astype(float), re=d["re"].astype(float),
                         im=d["im"].astype(float), nm=d["nm"].astype(float),
                         N=int(d["N_site"]), Nt=int(d["N_tri"])))
    return runs


class Reweighter:
    """Multi-histogram (WHAM) over all ensembles of all runs. Energies are
    binned exactly on their unique values (integers here, since D2, D3 and
    the bond costs are integers), and per-bin sums of the psi moments are
    pooled over all ensembles, so reweighting costs O(#energy levels)."""

    def __init__(self, runs):
        self.N, self.Nt = runs[0]["N"], runs[0]["Nt"]
        E_all, P2_all, NM_all, beta_of, n = [], [], [], [], []
        for r in runs:
            for k, b in enumerate(r["betas"]):
                E_all.append(np.round(r["E"][:, k], 6))
                P2_all.append(r["re"][:, k] ** 2 + r["im"][:, k] ** 2)
                NM_all.append(r["nm"][:, k])
                beta_of.append(b)
                n.append(r["E"].shape[0])
        E_all, P2_all, NM_all = np.concatenate(E_all), np.concatenate(P2_all), np.concatenate(NM_all)
        self.lev, inv = np.unique(E_all, return_inverse=True)
        nl = len(self.lev)
        self.H = np.bincount(inv, minlength=nl).astype(float)
        self.S1 = np.bincount(inv, weights=np.sqrt(P2_all), minlength=nl)
        self.S2 = np.bincount(inv, weights=P2_all, minlength=nl)
        self.S4 = np.bincount(inv, weights=P2_all ** 2, minlength=nl)
        self.Sn = np.bincount(inv, weights=NM_all, minlength=nl)
        self.Snn = np.bincount(inv, weights=NM_all ** 2, minlength=nl)
        self.b = np.array(beta_of)
        self.n = np.array(n, float)
        e = self.lev - self.lev.mean()
        self.e = e
        f = np.zeros(len(self.b))
        logH = np.log(self.H)
        for _ in range(20000):
            logden = logsumexp(np.log(self.n)[:, None] + f[:, None] - np.outer(self.b, e), axis=0)
            fn = -logsumexp(logH[None, :] - np.outer(self.b, e) - logden[None, :], axis=1)
            fn -= fn[0]
            if np.max(np.abs(fn - f)) < 1e-10:
                f = fn
                break
            f = fn
        self.logg = logH - logsumexp(np.log(self.n)[:, None] + f[:, None] - np.outer(self.b, e), axis=0)

    def moments(self, beta):
        lw = self.logg - beta * self.e
        w = np.exp(lw - lw.max())
        Z = w.sum()
        m1 = (w * self.S1 / self.H).sum() / Z
        m2 = (w * self.S2 / self.H).sum() / Z
        m4 = (w * self.S4 / self.H).sum() / Z
        e1 = (w * self.lev).sum() / Z
        e2 = (w * self.lev ** 2).sum() / Z
        n1 = (w * self.Sn / self.H).sum() / Z
        n2 = (w * self.Snn / self.H).sum() / Z
        en = (w * self.lev * self.Sn / self.H).sum() / Z
        varE, varn, cov = e2 - e1 ** 2, n2 - n1 ** 2, en - e1 * n1
        return dict(U=1 - m4 / (2 * m2 ** 2), chi=beta * self.N * (m2 - m1 ** 2), m2=m2,
                    C=beta ** 2 * varE / self.Nt,
                    # canonical (fixed monomer number) heat capacity, Fisher's constrained ensemble:
                    # C_n = C_mu - beta^2 cov(E,n)^2 / var(n)
                    Cn=beta ** 2 * (varE - cov ** 2 / varn) / self.Nt,
                    rho=cov / np.sqrt(varE * varn), nmon=n1 / self.Nt)


def curves(runs, bgrid):
    rw = Reweighter(runs)
    out = {k: np.array([rw.moments(b)[k] for b in bgrid]) for k in ("U", "chi", "m2", "C", "Cn", "rho", "nmon")}
    out["dU"] = np.gradient(out["U"], bgrid)
    out["N"] = rw.N
    return out


def observables(D2, Ls, bgrid, drop=None):
    """Curves per L; drop = seed index to leave out (jackknife)."""
    res = {}
    for L in Ls:
        runs = load_L(D2, L)
        if drop is not None:
            runs = [r for r in runs if r["group"] != drop]
        bmin = min(r["betas"].min() for r in runs)
        bmax = max(r["betas"].max() for r in runs)
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


TC_FIXED = None


def main():
    global TC_FIXED
    import os
    TC_FIXED = float(os.environ["TC"]) if os.environ.get("TC") else None
    global MU_TAG
    MU_TAG = f"_mu{float(os.environ['MU']):g}" if os.environ.get("MU") else ""
    D2 = float(sys.argv[1])
    Ls = [int(x) for x in sys.argv[2:]]
    groups = sorted(set.intersection(*[{r["group"] for r in load_L(D2, L)} for L in Ls]))
    nseed = len(groups)
    allb = np.concatenate([r["betas"] for L in Ls for r in load_L(D2, L)])
    bgrid = np.linspace(allb.min(), allb.max(), 400)

    def estimate(drop=None):
        res = observables(D2, Ls, bgrid, drop)
        tcs = [crossing(res[a][0], res[a][1]["U"], res[b][0], res[b][1]["U"])
               for a, b in zip(Ls[:-1], Ls[1:])]
        Tc = TC_FIXED if TC_FIXED else (tcs[-1] if np.isfinite(tcs[-1]) else np.nanmean(tcs))
        return res, tcs, Tc, fit_exponents(res, Ls, Tc)

    res, tcs, Tc, ex = estimate()
    jk = [estimate(g) for g in groups] if nseed > 1 else []

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

    # first-order test: reweighted energy distribution at each L's C-peak temperature
    print("  energy distribution at the C peak: T_peak, #maxima, barrier ln(Pmax/Pmin_between)")
    pe = {}
    for L in Ls:
        runs = load_L(D2, L)
        rw = Reweighter(runs)
        g, cv = res[L]
        bpk = g[np.argmax(cv["C"])]
        lp = rw.logg - bpk * rw.e
        P = np.exp(lp - lp.max())
        width = max(1.0, 0.02 * np.sqrt(rw.Nt))
        grid = np.arange(rw.lev.min(), rw.lev.max() + 1)
        dens = np.interp(grid, rw.lev, P)
        ker = np.exp(-0.5 * (np.arange(-4 * width, 4 * width + 1) / width) ** 2)
        sm = np.convolve(dens, ker / ker.sum(), mode="same")
        mx = [i for i in range(1, len(sm) - 1) if sm[i] > sm[i - 1] and sm[i] >= sm[i + 1] and sm[i] > 0.05 * sm.max()]
        bar = np.log(max(sm[mx[0]], sm[mx[-1]]) / sm[mx[0]:mx[-1] + 1].min()) if len(mx) > 1 else 0.0
        print(f"    L={L:3d}  T={1 / bpk:.4f}  maxima={len(mx)}  barrier={bar:.3f}")
        pe[L] = (grid / rw.Nt, sm / sm.max(), 1 / bpk)
    cols = ["#c5791f", "#3d7a4f", "#3b4ba8", "#a83b3b", "#6b4ba8", "#201d18"]
    fig, ax = plt.subplots(2, 4, figsize=(18, 8))
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
    for L, c in zip(Ls, cols):
        x, y, tp = pe[L]
        ax[0, 3].plot(x, y, color=c, label=f"L={L}, T={tp:.4f}")
    ax[0, 3].set_xlabel("E / N_tri")
    ax[0, 3].set_ylabel("P(E) at the C peak (smoothed, max=1)")
    ax[0, 3].legend(fontsize=7)
    a = ax[1, 3]
    a.loglog(Larr, ex["Cmax"], "o", color="k", label=f"C_max: slope {ex['alpha_nu']:.3f}")
    a.set_xlabel("L")
    a.set_ylabel("C_max / N_tri")
    a.legend(fontsize=7)
    fig.suptitle(f"Z3 crystallisation of the annealed model, 0=D1<D2={D2:g}: multi-histogram FSS (T_c = {Tc:.4f})")
    fig.tight_layout()
    out = os.environ.get("FIG", f"annealed_fss_d2{D2:g}_figure.png")
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()


def pseudo_critical(D2, Ls, U_level=0.25):
    """Pseudo-critical temperatures per L (C max, chi' max, dU/dbeta max,
    U4 = U_level on the ordered-side rise) and a joint extrapolation
    T_L = T_c + a_q L^(-1/nu) with common T_c and 1/nu."""
    from scipy.optimize import least_squares
    allb = np.concatenate([r["betas"] for L in Ls for r in load_L(D2, L)])
    bgrid = np.linspace(allb.min(), allb.max(), 800)
    res = observables(D2, Ls, bgrid)
    T = {}
    for L in Ls:
        g, cv = res[L]
        Tg = 1 / g
        up = np.where(cv["U"] >= U_level)[0]
        tU = np.nan
        if len(up) and up.min() > 0:
            i = up.min()                      # first beta (from the high-T side) with U above the level
            tU = 1 / np.interp(U_level, [cv["U"][i - 1], cv["U"][i]], [g[i - 1], g[i]])
        edge = lambda k: np.nan if k in (0, len(g) - 1) else Tg[k]
        T[L] = dict(C=edge(np.argmax(cv["C"])), chi=edge(np.argmax(cv["chi"])),
                    dU=edge(np.argmax(cv["dU"])), U=tU)
    keys = ["C", "chi", "dU", "U"]
    rows = [(L, q, T[L][q]) for L in Ls for q in keys if np.isfinite(T[L][q])]

    def resid(p):
        tc, x = p[0], p[1]
        a = dict(zip(keys, p[2:]))
        return [t - (tc + a[q] * L ** (-x)) for L, q, t in rows]

    fit = least_squares(resid, x0=[0.52, 1.0, 0.3, 0.3, 0.3, 0.3])
    return T, fit.x[0], fit.x[1], rows


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "pc":
    pass

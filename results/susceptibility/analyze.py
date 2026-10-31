"""Analysis of the susceptibility runs (chi.py -> results/susceptibility/raw/).

Usage (repo root): python3 results/susceptibility/analyze.py
Writes results/susceptibility/summary.txt, chi_vs_T.pdf, m_vs_H.pdf.
Errors: for seed-averaged values, the larger of the seed scatter (std/sqrt(n))
and the combined block-jackknife error (as required by the task).
"""
import glob
import os
import re

import numpy as np

D = "results/susceptibility"
CW_THETA = {"tafm": 3.0, "annealed": 2.4}


def read(fn):
    out = {}
    for line in open(fn):
        if line.startswith("#"):
            continue
        k, *v = line.split()
        out[k] = [float(x) for x in v]
    return out


def runs_T():
    rows = []
    for fn in glob.glob(f"{D}/raw/T_*.txt"):
        b = os.path.basename(fn)[2:-4]
        m = re.match(r"(tafm|annealed)(?:_mu([\d.]+))?_d2([\d.e-]+)_L(\d+)_T([\d.]+)_s(\d+)_(\w+)$", b)
        model, mu, D2, L, T, s, st = m.groups()
        r = read(fn)
        rows.append(dict(model=model, mu=float(mu or 0), D2=float(D2), L=int(L), T=float(T), seed=int(s), start=st,
                         m=r["m"], tchi=r["Tchi"], exc=r.get("nmon", r.get("defect"))))
    return rows


def runs_gs():
    rows = []
    for fn in glob.glob(f"{D}/raw/gs_*.txt"):
        m = re.match(r"gs_L(\d+)_H([\d.]+)_s(\d+)_(\w+)\.txt", os.path.basename(fn))
        L, H, s, st = m.groups()
        r = read(fn)
        rows.append(dict(L=int(L), H=float(H), seed=int(s), start=st, m=r["m"], tchi=r["Tchi"], sq=np.array(r["Sq"])))
    return rows


def combine(vals, errs):
    vals, errs = np.asarray(vals, float), np.asarray(errs, float)
    n = len(vals)
    jk = np.sqrt(np.sum(errs ** 2)) / n
    sc = vals.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
    return vals.mean(), max(jk, sc), n


def group(rows, keys, field):
    g = {}
    for r in rows:
        g.setdefault(tuple(r[k] for k in keys), []).append(r)
    return {k: combine([r[field][0] for r in v], [r[field][1] for r in v]) for k, v in g.items()}, g


def main():
    rows = runs_T()
    gsr = runs_gs()
    lines = []
    w = lines.append
    w("# Uniform susceptibility chi = var(n1)/(N_site T) and reduced-field scaling")
    w("# annealed: pure Ising sweep (sweep_ising); TAFM: fixed bonds; all runs start from an exact ground state")
    w("# errors: max(seed scatter, block jackknife); 'n' = number of seeds")
    w("")

    # ---------------- 1. chi(T) ----------------
    r1 = [r for r in rows if r["D2"] == 0 and r["start"] == "gs"]
    tchi, _ = group(r1, ("model", "mu", "L", "T"), "tchi")
    exc, _ = group(r1, ("model", "mu", "L", "T"), "exc")
    w("## 1. T*chi(T) at D2=0")
    sets = sorted({k[:3] for k in tchi})
    Ts = sorted({k[3] for k in tchi})
    hdr = f"{'T':>6} | " + " | ".join(f"{m[:4]} mu={mu:g} L={L:<3d}      " for m, mu, L in sets)
    w(hdr)
    w("-" * len(hdr))
    for T in Ts:
        cells = []
        for st in sets:
            v = tchi.get(st + (T,))
            cells.append(f"{v[0]:.4f}({v[1]:.4f})[{v[2]}]" if v else f"{'--':>19s}")
        w(f"{T:6g} | " + " | ".join(cells))
    w("")
    w("# size check L=24 vs L=48 (z = difference / combined error):")
    for model in ("tafm", "annealed"):
        zs = []
        for T in Ts:
            a, b = tchi.get((model, 0.0, 24, T)), tchi.get((model, 0.0, 48, T))
            if a and b:
                zs.append((T, (b[0] - a[0]) / np.hypot(a[1], b[1])))
        if zs:
            w(f"#   {model:8s}: " + " ".join(f"{T:g}:{z:+.1f}" for T, z in zs) +
              f"   max|z| = {max(abs(z) for _, z in zs):.1f}")
    w("")

    # high-T Curie-Weiss: T*chi = (T/4)/(T+theta) -> theta_eff = T/(4 Tchi) - T
    w("# high-T Curie-Weiss check: theta_eff(T) = T/(4 T*chi) - T  (expected TAFM 3, annealed mu=0 ~2.4)")
    w("# theta(T->inf): linear extrapolation of theta_eff in 1/T over T = 3, 5, 10")
    theta = {}
    for st in sets:
        pts = []
        for T in (3.0, 5.0, 10.0):
            v = tchi.get(st + (T,))
            if v:
                th = T / (4 * v[0]) - T
                eth = T / (4 * v[0] ** 2) * v[1]
                pts.append((T, th, eth))
        if len(pts) >= 2:
            x = np.array([1 / p[0] for p in pts])
            y = np.array([p[1] for p in pts])
            e = np.array([p[2] for p in pts])
            A = np.vstack([np.ones_like(x), x]).T
            W = np.diag(1 / e ** 2)
            cov = np.linalg.inv(A.T @ W @ A)
            p = cov @ A.T @ W @ y
            theta[st] = (p[0], np.sqrt(cov[0, 0]))
            w(f"#   {st[0]:8s} mu={st[1]:<3g} L={st[2]:<3d}: " +
              " ".join(f"T={T:g}: {th:.2f}({eth:.2f})" for T, th, eth in pts) +
              f"  -> theta(inf) = {p[0]:.2f}({np.sqrt(cov[0, 0]):.2f})")
    w("")

    # 4. TAFM vs annealed dip
    w("## 4. Low-T T*chi and excitation densities (L=24 and 48, mu=0)")
    w("#   annealed exc = frustrated triangles / N_tri; TAFM exc = triangles with three equal spins / N_tri")
    for L in (24, 48):
        for T in Ts:
            if T > 0.8:
                continue
            a, b = tchi.get(("tafm", 0.0, L, T)), tchi.get(("annealed", 0.0, L, T))
            ea, eb = exc.get(("tafm", 0.0, L, T)), exc.get(("annealed", 0.0, L, T))
            if a and b:
                w(f"#   L={L} T={T:<5g} TAFM {a[0]:.4f}({a[1]:.4f}) exc {ea[0]:.2e} | annealed {b[0]:.4f}({b[1]:.4f}) "
                  f"exc {eb[0]:.2e} | diff ann-TAFM = {b[0] - a[0]:+.4f} ({(b[0] - a[0]) / np.hypot(a[1], b[1]):+.1f} sigma)")
    w("")

    # ---------------- 2. chi0 ----------------
    w("## 2. chi0 = T*chi on the ground-state manifold (H=0, zero-energy flips, tilt-0 sector)")
    g0 = [r for r in gsr if r["H"] == 0 and r["start"] == "gs"]
    c0, gg = group(g0, ("L",), "tchi")
    chi0 = {}
    for (L,), v in sorted(c0.items()):
        chi0[L] = v
        w(f"#   L={L:3d}: chi0 = {v[0]:.5f}({v[1]:.5f})  [{v[2]} seeds]")
    if len(chi0) >= 2:
        Ls = np.array(sorted(chi0))
        y = np.array([chi0[L][0] for L in Ls])
        e = np.array([chi0[L][1] for L in Ls])
        sl = np.polyfit(np.log(Ls), y, 1, w=1 / e)
        w(f"#   fit chi0 = a + b ln L: b = {sl[0]:+.5f}  (log divergence would give b > 0 well outside errors)")
    w("# uniform structure factor S(q_k), q_k = 2 pi k / L (both reciprocal directions averaged), seeds averaged;")
    w("# S(q->0) should approach chi0; height-model picture: S(q) - chi0 ~ q^2 (up to logs)")
    for L in sorted({r['L'] for r in g0}):
        sq = np.mean([r["sq"] for r in g0 if r["L"] == L], axis=0)
        esq = np.std([r["sq"] for r in g0 if r["L"] == L], axis=0, ddof=1) / np.sqrt(sum(r["L"] == L for r in g0))
        q = 2 * np.pi * np.arange(len(sq)) / L
        w(f"#   L={L:3d}: " + " ".join(f"k={k}:{sq[k]:.4f}({esq[k]:.4f})" for k in range(1, len(sq))))
        if L in chi0:
            w(f"#          (S(q)-chi0)/q^2: " + " ".join(f"{(sq[k] - chi0[L][0]) / q[k] ** 2:+.4f}" for k in range(1, len(sq))))
    w("")

    # ---------------- 3. m(H) ----------------
    w("## 3. Reduced-field scaling: m = <n1>/N_site vs H = D2/(2T)")
    mgs, _ = group([r for r in gsr], ("L", "H", "start"), "m")
    Hs = sorted({k[1] for k in mgs})
    w("# m_GS(H): ground-state manifold with field H; ord / gs = crystal / random-tiling start")
    for H in Hs:
        w(f"#   H={H:<5g} " + "  ".join(
            f"L={L} {st}: {mgs[(L, H, st)][0]:.5f}({mgs[(L, H, st)][1]:.5f})"
            for L in (24, 48) for st in ("ord", "gs") if (L, H, st) in mgs))
    r3 = [r for r in rows if r["start"] == "ord" and r["mu"] == 0]
    for r in r3:
        r["H"] = round(r["D2"] / (2 * r["T"]), 4)
    m3, _ = group(r3, ("model", "L", "T", "H"), "m")
    w("# m(D2, T) of the two models at H = D2/(2T) (ordered start, 2 seeds), deviation from m_GS(H) (same L, ord start)")
    for L in (24, 48):
        for T in (0.2, 0.3, 0.4, 0.5):
            for H in sorted({k[3] for k in m3 if k[1] == L and k[2] == T}):
                ref = mgs.get((L, H, "ord"))
                cells = []
                for model in ("tafm", "annealed"):
                    v = m3.get((model, L, T, H))
                    if v:
                        d = f" d={v[0] - ref[0]:+.5f}" if ref else ""
                        cells.append(f"{model:8s} {v[0]:.5f}({v[1]:.5f}){d}")
                w(f"#   L={L} T={T} H={H:<5g} " + " | ".join(cells))
    w("# closeness to the reference curve per T (L=48): rms over H of (m - m_GS), TAFM vs annealed")
    closeness = {}
    for T in (0.2, 0.3, 0.4, 0.5):
        out = []
        for model in ("tafm", "annealed"):
            d = [m3[(model, 48, T, H)][0] - mgs[(48, H, "ord")][0] for H in Hs
                 if (model, 48, T, H) in m3 and (48, H, "ord") in mgs]
            if d:
                closeness[(model, T)] = np.sqrt(np.mean(np.square(d)))
                out.append(f"{model} {closeness[(model, T)]:.5f}")
        w(f"#   T={T}: " + "  ".join(out))
    chk = [r for r in rows if r["start"] == "gs" and r["D2"] > 0]
    if chk:
        w("# start check (L=48, T=0.3): random-ground-state start vs ordered start")
        for r in chk:
            o = [x for x in r3 if x["model"] == r["model"] and x["L"] == r["L"] and x["T"] == r["T"]
                 and x["D2"] == r["D2"] and x["seed"] == r["seed"]]
            if o:
                o = o[0]
                w(f"#   {r['model']:8s} D2={r['D2']:g}: gs {r['m'][0]:.5f}({r['m'][1]:.5f}) ord {o['m'][0]:.5f}({o['m'][1]:.5f}) "
                  f"z={(r['m'][0] - o['m'][0]) / np.hypot(r['m'][1], o['m'][1]):+.1f}")
    text = "\n".join(lines) + "\n"
    open(f"{D}/summary.txt", "w").write(text)
    print(text)
    plots(tchi, sets, chi0, theta, mgs, m3)


def plots(tchi, sets, chi0, theta, mgs, m3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {("tafm", 0.0): "#3b4ba8", ("annealed", 0.0): "#c5791f", ("annealed", 2.0): "#3d7a4f",
           ("annealed", 10.0): "#a83b3b"}
    mk = {24: "o", 48: "s"}
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    c0 = np.mean([v[0] for v in chi0.values()]) if chi0 else 0.041
    Tg = np.logspace(np.log10(0.15), 1.05, 100)
    for st in sets:
        Ts = sorted(k[3] for k in tchi if k[:3] == st)
        y = np.array([tchi[st + (T,)][0] for T in Ts])
        e = np.array([tchi[st + (T,)][1] for T in Ts])
        lab = f"{'TAFM' if st[0] == 'tafm' else 'annealed'}{'' if st[0] == 'tafm' else f' mu={st[1]:g}'} L={st[2]}"
        c = col.get(st[:2], "k")
        ax[0].errorbar(Ts, y / np.array(Ts), e / np.array(Ts), fmt=mk[st[2]], ms=3.5, mfc="none" if st[2] == 24 else c,
                       color=c, label=lab, lw=0.8)
        ax[1].errorbar(Ts, y, e, fmt=mk[st[2]] + "-", ms=3.5, mfc="none" if st[2] == 24 else c, color=c, label=lab, lw=0.8)
    ax[0].plot(Tg, c0 / Tg, "k:", lw=1, label=rf"$\chi_0/T$, $\chi_0={c0:.4f}$")
    for m, th in CW_THETA.items():
        ax[0].plot(Tg, 0.25 / (Tg + th), "--", color=col[(m, 0.0)], lw=0.8, label=rf"$(1/4)/(T+{th:g})$")
        ax[1].plot(Tg[Tg > 1], (Tg[Tg > 1] / 4) / (Tg[Tg > 1] + th), "--", color=col[(m, 0.0)], lw=0.8)
    ax[1].axhline(c0, color="k", ls=":", lw=1)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("$T$")
    ax[0].set_ylabel(r"$\chi$")
    ax[0].set_title("(a)", loc="left")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("$T$")
    ax[1].set_ylabel(r"$T\chi$")
    ax[1].set_title("(b)", loc="left")
    ax[0].legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(f"{D}/chi_vs_T.pdf")
    fig.savefig(f"{D}/chi_vs_T.png", dpi=130)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    tc = {0.2: "#2b2b8f", 0.3: "#3d7a4f", 0.4: "#c5791f", 0.5: "#a83b3b"}
    for a, model in zip(ax, ("tafm", "annealed")):
        for L, ls in ((48, "-"), (24, "--")):
            H = sorted(k[1] for k in mgs if k[0] == L and k[2] == "ord")
            if H:
                a.plot(H, [mgs[(L, h, "ord")][0] for h in H], "k" + ls, lw=1, label=f"$m_{{GS}}(H)$, L={L}")
        for T, c in tc.items():
            H = sorted(k[3] for k in m3 if k[0] == model and k[1] == 48 and k[2] == T)
            if H:
                a.errorbar(H, [m3[(model, 48, T, h)][0] for h in H], [m3[(model, 48, T, h)][1] for h in H],
                           fmt="o", ms=3.5, color=c, label=f"T={T} (L=48)")
        a.axvline(0.266, color="gray", lw=0.7, ls=":")
        a.set_xlabel(r"$H = D_2/(2T)$")
        a.set_title("fixed-bond TAFM" if model == "tafm" else "annealed (Ising sweep, $\\mu=0$)")
        a.legend(fontsize=7)
    ax[0].set_ylabel(r"$m = \langle n_1\rangle / N_{\rm site}$")
    fig.tight_layout()
    fig.savefig(f"{D}/m_vs_H.pdf")
    fig.savefig(f"{D}/m_vs_H.png", dpi=130)


if __name__ == "__main__":
    main()

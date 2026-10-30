"""Part IX task 2b: T_c(D2) at small D2, annealed model (pure Ising sweep, mu=0)
vs the fixed-bond TAFM, from Binder crossings of the SPIN order parameter psi_s.

For every (D2, model) the long single-T runs of all sizes are combined per L by
multi-histogram reweighting (annealed_long_analysis.py: ISING=1 ORDER=spin for
the annealed runs in results/long, DIR=results/tafm for the TAFM, whose re/im
already are the spin order parameter). Errors: jackknife over NB time blocks.
Also checks ordered vs random starts (L=32, files *_random.npz) run by run.

Usage (repo root): [DISCARD=20000] python3 results/smallD2/tc_small_d2.py
Writes results/smallD2/summary.txt, results/smallD2/crossings_<model>_d2<D2>.txt,
results/smallD2/tc_vs_D2.png
"""
import glob
import os
import re
import subprocess

import numpy as np

HC = 0.266                     # TAFM ground-state KT field (Bloete-Nightingale)
D2S = [0.5, 0.3, 0.2, 0.1]
DISCARD = os.environ.get("DISCARD", "20000")
OUT = "results/smallD2"
KNOWN = {("tafm", 1.0): (0.4713, 0.0010), ("ann", 1.0): (0.528, 0.002), ("tafm", 3.0): (0.6718, 0.0007)}


def sizes(model, D2):
    pat = (f"results/long/long_L*_d2{D2:g}_T*_s1_ising.npz" if model == "ann"
           else f"results/tafm/long_L*_d2{D2:g}_T*_s1.npz")
    Ls = sorted({int(re.search(r"long_L(\d+)_", f).group(1)) for f in glob.glob(pat)})
    return [L for L in Ls if len(glob.glob(pat.replace("L*", f"L{L}"))) >= 3]


def crossings(model, D2):
    Ls = sizes(model, D2)
    if len(Ls) < 2:
        return Ls, [], None
    env = dict(os.environ, DISCARD=DISCARD, FIG=f"{OUT}/fss_{model}_d2{D2:g}.png")
    if model == "ann":
        env.update(ISING="1", ORDER="spin", DIR="results/long")
    else:
        env.update(DIR="results/tafm")
        env.pop("ISING", None)
    txt = subprocess.run(["python3", "annealed_long_analysis.py", f"{D2:g}"] + [str(L) for L in Ls],
                         env=env, capture_output=True, text=True).stdout
    open(f"{OUT}/crossings_{model}_d2{D2:g}.txt", "w").write(txt)
    pairs = [(int(a), int(b), float(t), float(e)) for a, b, t, e in
             re.findall(r"Binder crossing L=(\d+)/(\d+): T = (\S+) \+- (\S+)", txt)]
    return Ls, pairs, txt


def best(pairs):
    """T_c estimate: crossing of the largest pair with a finite value."""
    ok = [p for p in pairs if np.isfinite(p[2])]
    return (ok[-1][2], ok[-1][3], f"{ok[-1][0]}/{ok[-1][1]}") if ok else (np.nan, np.nan, "-")


def start_check():
    """Ordered vs random start, run by run (same L, T, D2, model)."""
    rows = []
    for fr in sorted(glob.glob("results/long/long_L*_d2*_T*_s1_random_ising.npz")
                     + glob.glob("results/tafm/long_L*_d2*_T*_s1_random.npz")):
        fo = fr.replace("_random", "")
        if not os.path.exists(fo) or "_d21_" in fr:
            continue
        out = []
        for f in (fo, fr):
            d = np.load(f)
            k = int(DISCARD) // int(d["thin"])
            re_, im_ = (d["sre"], d["sim"]) if "sre" in d.files else (d["re"], d["im"])
            p2 = re_[k:, 0].astype(float) ** 2 + im_[k:, 0].astype(float) ** 2
            E = d["E"][k:, 0].astype(float) / int(d["N_tri"])
            nb = 8
            n = len(p2) // nb * nb
            U = lambda x: 1 - np.mean(x ** 2) / (2 * np.mean(x) ** 2)
            bl = lambda f_, x: np.array([f_(np.delete(x[:n], np.s_[b * n // nb:(b + 1) * n // nb])) for b in range(nb)])
            jk = lambda v: np.sqrt((nb - 1) / nb * np.sum((v - v.mean()) ** 2))
            out.append([(U(p2), jk(bl(U, p2))), (p2.mean(), jk(bl(np.mean, p2))), (E.mean(), jk(bl(np.mean, E)))])
        name = os.path.basename(fo)
        model = "ann" if "_ising" in name else "tafm"
        L = int(re.search(r"_L(\d+)_", name).group(1))
        D2 = float(re.search(r"_d2([\d.]+)_", name).group(1))
        T = float(re.search(r"_T([\d.]+)_", name).group(1))
        rows.append((model, D2, L, T, out))
    return sorted(rows)


def main():
    lines = []
    w = lines.append
    w("# Part IX 2b: T_c(D2) at small D2 -- annealed (pure Ising sweep, mu=0) vs fixed-bond TAFM")
    w(f"# Binder crossings of the spin order parameter psi_s; WHAM per L; jackknife over 8 time blocks; discard {DISCARD} sweeps")
    w(f"# H = D2/(2 T_c) (reduced field; TAFM ground-state KT value H_c = {HC}); asymptote T_c -> D2/(2 H_c) = {1 / (2 * HC):.3f} D2")
    w("")
    res = {}
    for D2 in D2S:
        for model in ("tafm", "ann"):
            Ls, pairs, _ = crossings(model, D2)
            res[(model, D2)] = (Ls, pairs)
            w(f"D2={D2:g} {model:4s} sizes {Ls}: " + ", ".join(f"L={a}/{b}: {t:.4f}({e:.4f})" for a, b, t, e in pairs))
    w("")
    hdr = f"{'D2':>4} | {'T_c TAFM':>15} {'pair':>6} | {'T_c ann':>15} {'pair':>6} | {'ann/TAFM':>14} | {'H TAFM':>6} {'H ann':>6} | {'D2/(2Hc)':>8}"
    w(hdr)
    w("-" * len(hdr))
    table = []
    for D2 in D2S:
        (tt, et, pt), (ta, ea, pa) = best(res[("tafm", D2)][1]), best(res[("ann", D2)][1])
        r = ta / tt
        er = r * np.hypot(ea / ta, et / tt) if np.isfinite(r) else np.nan
        table.append((D2, tt, et, ta, ea))
        w(f"{D2:4g} | {tt:.4f}({et:.4f}) {pt:>6} | {ta:.4f}({ea:.4f}) {pa:>6} | {r:.3f}({er:.3f}) | "
          f"{D2 / (2 * tt):6.3f} {D2 / (2 * ta):6.3f} | {D2 / (2 * HC):8.3f}")
    for (m, D2), (t, e) in sorted(KNOWN.items(), key=lambda x: x[0][1]):
        w(f"# reference: D2={D2:g} {m}: T_c = {t}({e})  H = {D2 / (2 * t):.3f}")
    w("")
    w("# ordered vs random start (same L, T): U4, psi_s^2, E/N_tri  [z = (ordered - random)/combined error]")
    for model, D2, L, T, (o, r_) in start_check():
        zs = [(a[0] - b[0]) / np.hypot(a[1], b[1]) if np.hypot(a[1], b[1]) > 0 else 0 for a, b in zip(o, r_)]
        w(f"  {model:4s} D2={D2:<4g} L={L} T={T:<6g} U4 {o[0][0]:.3f}/{r_[0][0]:.3f} (z={zs[0]:+.1f})  "
          f"psi2 {o[1][0]:.4f}/{r_[1][0]:.4f} (z={zs[1]:+.1f})  E {o[2][0]:.5f}/{r_[2][0]:.5f} (z={zs[2]:+.1f})")
    text = "\n".join(lines) + "\n"
    open(f"{OUT}/summary.txt", "w").write(text)
    print(text)
    plot(table)


def plot(table):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array(table, float)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    a = ax[0]
    a.errorbar(t[:, 0], t[:, 1], t[:, 2], fmt="o", color="#3b4ba8", label="fixed-bond TAFM")
    a.errorbar(t[:, 0], t[:, 3], t[:, 4], fmt="s", color="#c5791f", label="annealed (Ising sweep)")
    a.errorbar([1.0, 3.0], [0.4713, 0.6718], [0.001, 0.0007], fmt="o", mfc="none", color="#3b4ba8", label="TAFM, earlier runs")
    a.errorbar([1.0], [0.528], [0.002], fmt="s", mfc="none", color="#c5791f", label="annealed D3=10, earlier runs")
    x = np.linspace(0, 1.05, 50)
    a.plot(x, x / (2 * HC), "k--", lw=1, label=r"$D_2/(2H_c)$, $H_c=0.266$")
    a.set_xlim(0, 1.05)
    a.set_ylim(0, 0.8)
    a.set_xlabel(r"$D_2$")
    a.set_ylabel(r"$T_c$")
    a.legend(fontsize=7)
    a = ax[1]
    a.errorbar(t[:, 0], t[:, 3] / t[:, 1], t[:, 3] / t[:, 1] * np.hypot(t[:, 4] / t[:, 3], t[:, 2] / t[:, 1]),
               fmt="o", color="k", label="this work")
    a.plot([1.0], [0.528 / 0.4713], "o", mfc="none", color="k", label="D2=1, earlier runs")
    a.axhline(1, color="gray", lw=0.8)
    a.set_xlabel(r"$D_2$")
    a.set_ylabel(r"$T_c^{\rm ann}/T_c^{\rm TAFM}$")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{OUT}/tc_vs_D2.png", dpi=150)


if __name__ == "__main__":
    main()

"""Ising-vs-(D3=10 Potts) validation of the Part IX annealed-bond results.

For every (mu, T, L) with both a pure-Ising run (annealed_long.py with
ISING=1, files *_ising.npz) and the existing D3=10 reference run (same name
without the suffix), compare

  U4      = 1 - <|psi_s|^4> / (2 <|psi_s|^2>^2)   (spin order parameter)
  psi_s^2 = <|psi_s|^2>
  E/N_tri = <E>/N_tri

after discarding the first DISCARD sweeps, with errors from a jackknife over
NB equal time blocks (the convention of annealed_long_analysis.py).

E: the D3=10 energy contains D3 * n3 (third-state sites), which the Ising
model does not have. Reference runs older than this check did not record n3,
so the third-state density per site is taken from short D3=10 runs (L=32,
seed 77, same mu and T; n3 is a local quantity) unless the reference file
itself carries n3. Both raw E and E - D3*n3 are reported.

Pass criterion: seed-to-seed scatter in these runs is up to ~2x the
autocorrelation error (HANDOFF result 14), so the combined error is inflated
by SEED_FACTOR=2 and a quantity passes if |difference| < 2 * inflated error.

Usage: python3 results/ising_check/ising_check.py  (from the repo root)
Writes results/ising_check/summary.txt
"""
import glob
import os

import numpy as np

DISCARD = int(os.environ.get("DISCARD", 100000))
NB = int(os.environ.get("NB", 8))
D3 = 10.0
SEED_FACTOR = 2.0
OUT = "results/ising_check/summary.txt"


def blocks_stat(arrs, f):
    """f(list of block-trimmed arrays) -> value; jackknife over NB time blocks."""
    n = len(arrs[0])
    edges = np.linspace(0, n, NB + 1).astype(int)
    full = f(arrs)
    jk = []
    for b in range(NB):
        keep = np.r_[0:edges[b], edges[b + 1]:n]
        jk.append(f([a[keep] for a in arrs]))
    jk = np.array(jk)
    return full, float(np.sqrt((NB - 1) / NB * np.sum((jk - jk.mean()) ** 2)))


def measure(fn):
    d = np.load(fn)
    thin = int(d["thin"])
    k0 = DISCARD // thin
    E = d["E"][k0:, 0].astype(float) / int(d["N_tri"])
    p2 = (d["sre"][k0:, 0].astype(float) ** 2 + d["sim"][k0:, 0].astype(float) ** 2)
    out = {
        "U4": blocks_stat([p2], lambda a: 1 - np.mean(a[0] ** 2) / (2 * np.mean(a[0]) ** 2)),
        "p2": blocks_stat([p2], lambda a: np.mean(a[0])),
        "E": blocks_stat([E], lambda a: np.mean(a[0])),
        "sweeps": thin * d["E"].shape[0],
        "Nsite": int(d["N_site"]), "Ntri": int(d["N_tri"]),
    }
    if "n3" in d.files:
        n3 = d["n3"][k0:, 0].astype(float) / int(d["N_site"])
        out["n3"] = blocks_stat([n3], lambda a: np.mean(a[0]))
    return out


def n3_density(mu, T):
    """Third-state density per site of the D3=10 model at (mu, T), from short L=32 runs."""
    tag = f"_mu{mu:g}" if mu else ""
    f = f"results/long/long_L32_d21{tag}_T{T:g}_s77.npz"
    if not os.path.exists(f):
        return None
    d = np.load(f)
    k0 = 50000 // int(d["thin"])
    n3 = d["n3"][k0:, 0].astype(float) / int(d["N_site"])
    return blocks_stat([n3], lambda a: np.mean(a[0]))


def main():
    rows, skipped = [], []
    for fi in sorted(glob.glob("results/long/long_L*_d21*_s1_ising.npz")):
        fr = fi.replace("_ising.npz", ".npz")
        if not os.path.exists(fr):
            continue
        name = os.path.basename(fr)
        L = int(name.split("_")[1][1:])
        mu = float(name.split("_mu")[1].split("_")[0]) if "_mu" in name else 0.0
        T = float(name.split("_T")[1].split("_")[0])
        if not all("sre" in np.load(f).files for f in (fi, fr)):
            skipped.append(name)                     # reference predates the spin order parameter
            continue
        rows.append((mu, T, L, measure(fi), measure(fr)))
    rows.sort()
    lines = []
    w = lines.append
    w("# Ising (sigma in {0,1}) vs D3=10 Potts implementation of the annealed-bond model, D2=1")
    w(f"# discard {DISCARD} sweeps, jackknife over {NB} time blocks; pass = |diff| < 2 * {SEED_FACTOR:g} * combined error")
    w("# E_D3 includes D3*n3; E_D3-D3n3 subtracts it (n3 from the reference file or from L=32 seed-77 D3=10 runs)")
    w("")
    hdr = (f"{'mu':>4} {'T':>6} {'L':>3} | {'U4 Ising':>15} {'U4 D3=10':>15} {'z':>5} | "
           f"{'psi_s^2 Ising':>17} {'psi_s^2 D3=10':>17} {'z':>5} | {'E/Ntri Ising':>19} {'E/Ntri D3=10':>19} "
           f"{'z':>5} | {'E-D3n3 D3=10':>19} {'z':>5} | {'n3/site':>9} | pass")
    w(hdr)
    w("-" * len(hdr))
    n_pass = n_tot = 0
    for mu, T, L, a, b in rows:
        n3 = b.get("n3") or n3_density(mu, T)
        zs = []
        for key in ("U4", "p2", "E"):
            (va, ea), (vb, eb) = a[key], b[key]
            zs.append((va - vb) / np.hypot(ea, eb))
        (vb, eb) = b["E"]
        if n3 is not None:
            ec = vb - D3 * n3[0] * b["Nsite"] / b["Ntri"]
            ecerr = np.hypot(eb, D3 * n3[1] * b["Nsite"] / b["Ntri"])
            zc = (a["E"][0] - ec) / np.hypot(a["E"][1], ecerr)
            ec_s, zc_s, n3_s = f"{ec:.6f}({ecerr:.6f})", f"{zc:+5.1f}", f"{n3[0]:.2e}"
            ok_e = abs(zc) < 2 * SEED_FACTOR
        else:
            ec_s, zc_s, n3_s, ok_e = "n/a", "  n/a", "n/a", abs(zs[2]) < 2 * SEED_FACTOR
        passed = abs(zs[0]) < 2 * SEED_FACTOR and abs(zs[1]) < 2 * SEED_FACTOR and ok_e
        n_pass += passed
        n_tot += 1
        w(f"{mu:4g} {T:6g} {L:3d} | {a['U4'][0]:.4f}({a['U4'][1]:.4f}) {b['U4'][0]:.4f}({b['U4'][1]:.4f}) {zs[0]:+5.1f} | "
          f"{a['p2'][0]:.5f}({a['p2'][1]:.5f}) {b['p2'][0]:.5f}({b['p2'][1]:.5f}) {zs[1]:+5.1f} | "
          f"{a['E'][0]:.6f}({a['E'][1]:.6f}) {b['E'][0]:.6f}({b['E'][1]:.6f}) {zs[2]:+5.1f} | "
          f"{ec_s:>19} {zc_s} | {n3_s:>9} | {'PASS' if passed else 'FAIL'}")
    w("")
    w(f"# {n_pass}/{n_tot} (mu, T, L) points pass (U4, psi_s^2 and E-D3*n3 all within 2 x {SEED_FACTOR:g} x combined error)")
    w("# z = (Ising - D3=10) / combined jackknife error (not inflated); sweeps per run as in the files")
    w("")
    if skipped:
        w("# skipped (reference run has no spin order parameter): " + ", ".join(skipped))
    w("# third-state density per site of the D3=10 model (short L=32 runs, seed 77):")
    for mu, T in sorted({(r[0], r[1]) for r in rows} | {(10.0, 1.04), (10.0, 1.08)}):
        n3 = n3_density(mu, T)
        if n3:
            w(f"#   mu={mu:g} T={T:g}: n3/site = {n3[0]:.2e} +- {n3[1]:.1e}")
    text = "\n".join(lines) + "\n"
    open(OUT, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()

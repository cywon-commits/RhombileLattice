"""Entropy of the pure-Ising annealed-bond model by thermodynamic integration
(HANDOFF task 3): the same sequential-cooling protocol as annealed_thermo.scan
(D2=0, mu=0, L=12; main grid T=0.06-6, 80 points uniform in beta; high-T
supplement T=6-300, 24 points), but with ising_sweep.sweep_ising, and the
beta=0 anchor of the Ising model:

  ln Z(0)/N_tri = (N_site/N_tri) ln 2 + ln M/N_tri = 0.5 ln 2 + 0.582
  e(beta=0)     = (1.5 - d_inf)/2 + D2/4 + mu (1 - 2 d_inf)

(d_inf = dimers per triangle of uniform matchings, from
results/thermo/matchings_L16.log). The D3=10 integration started from
0.5 ln 3 + 0.582 instead; the third state freezes out above T ~ 2, so S(T)
below that, and S(T=0.06) = 0.150, should be unchanged.

Usage (repo root): PYTHONPATH=. python3 results/ising_check/ising_entropy.py
Writes results/ising_check/entropy_ising_L12.txt
"""
import re

import numpy as np
from numba import njit

from annealed_thermo import build_arrays, total_energy, observe, _seed
from ising_sweep import sweep_ising

S_TILING = 0.323065947 / 2


@njit(cache=True)
def run_T(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E, n_eq, n_meas):
    for _ in range(n_eq):
        E = sweep_ising(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E)
    acc = np.zeros(3)
    for _ in range(n_meas):
        E = sweep_ising(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E)
        nm, n3, n2, hub = observe(s, mate, dimer, site_b)
        acc[0] += E
        acc[1] += E * E
        acc[2] += nm
    return E, acc / n_meas


def scan(A, seed, tmin, tmax, nT, n_eq, n_meas, D2=0.0, mu=0.0):
    _seed(seed)
    rng = np.random.default_rng(seed)
    s = rng.integers(2, size=A["ns"]).astype(np.int64)
    dimer, mate = A["dimer"].copy(), A["mate"].copy()
    E = total_energy(s, dimer, A["bi"], A["bj"], D2, 10.0, mu, mate)
    out = []
    for beta in np.linspace(1.0 / tmax, 1.0 / tmin, nT):
        E, a = run_T(s, dimer, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"],
                     beta, D2, mu, E, n_eq, n_meas)
        assert abs(total_energy(s, dimer, A["bi"], A["bj"], D2, 10.0, mu, mate) - E) < 1e-6
        nt = A["nt"]
        out.append((1 / beta, beta, a[0] / nt, beta ** 2 * (a[1] - a[0] ** 2) / nt, a[2] / nt))
    return np.array(out)


def main():
    txt = open("results/thermo/matchings_L16.log").read()
    lnM = float(re.search(r"lnM/N_tri = ([\d.]+)", txt).group(1))
    d_inf = [float(l.split()[1]) for l in txt.splitlines() if l[:1].isdigit()][-1]
    A = build_arrays(12)
    main_ = scan(A, 2, 0.06, 6.0, 80, 3000, 20000)
    high = scan(A, 5, 6.0, 300.0, 24, 1000, 8000)
    e0 = (1.5 - d_inf) / 2
    hb = high[high[:, 1] < main_[0, 1] - 1e-9]
    beta = np.concatenate([[0.0], hb[:, 1], main_[:, 1]])
    e = np.concatenate([[e0], hb[:, 2], main_[:, 2]])
    o = np.argsort(beta)
    beta, e = beta[o], e[o]
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(beta))])
    s = 0.5 * np.log(2) + lnM - integ + beta * e
    # D3=10 reference curve (same protocol, from results/thermo)
    ref = [l.split() for l in open("results/thermo/scan_L12_D10_mu0_lowT.log") if l[:1].isdigit()]
    ref = {round(float(r[0]), 4): float(r[2]) for r in ref}
    lines = ["# Pure-Ising annealed model, D2=0, mu=0, L=12: entropy by thermodynamic integration",
             f"# anchor s(beta=0) = 0.5 ln2 + lnM = {0.5 * np.log(2) + lnM:.4f} (lnM/N_tri = {lnM:.4f}); e(beta=0) = {e0:.5f}",
             f"# lozenge-tiling (TAFM ground-state) entropy per triangle: {S_TILING:.5f}",
             "# T  beta  E/Ntri(Ising)  E/Ntri(D3=10 same T)  C/Ntri  n_mon/Ntri  S/Ntri(Ising)"]
    sm = s[-len(main_):]
    for row, sv in zip(main_, sm):
        r = ref.get(round(row[0], 4))
        lines.append(f"{row[0]:.4f} {row[1]:.4f} {row[2]:.6f} {r if r is not None else float('nan'):.6f} "
                     f"{row[3]:.5f} {row[4]:.5f} {sv:.5f}")
    lines.append(f"# S(T={main_[-1, 0]:.3f}) = {sm[-1]:.4f}   (D3=10 integration: 0.150; exact T=0 value {S_TILING:.4f})")
    text = "\n".join(lines) + "\n"
    open("results/ising_check/entropy_ising_L12.txt", "w").write(text)
    print("\n".join(lines[:4] + lines[-6:]))


if __name__ == "__main__":
    main()

"""(c) D_c(p) on worm-built, site-isolated monomer configurations.

For each D3: SA ground state, energy E, total number of state-3 sites
n3 (at T=0, dE/dD3 = n3, so kinks in E(D3) are exactly where n3 jumps),
and the fraction of monomers with a state-3 corner, split into
site-isolated vs. site-sharing monomers.

Usage: python3 worm_Dc_vs_density.py <p> [seed]
"""
import sys

import numpy as np

from rhombile_lattice import simulated_annealing, frustrated_triangles
from worm_density_construction import build, report

D3_GRID = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.3]


def sa_best(lat, D, seed, n_random=3, n_sweeps=2000):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_random):
        states, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps,
                                          T_start=4.0, T_end=1e-6, record_energy=True)
        if best is None or en[-1] < best[0]:
            best = (en[-1], states)
    return best


def main():
    p = float(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 11
    rng = np.random.default_rng(seed)
    st, site_nb = build(p, rng, n_sweeps=600)
    mons = sorted(st.monomers)
    n_ft = len(frustrated_triangles(st.lat))
    edge_pairs, site_pairs, site_iso = report(st, site_nb)
    iso = [m for m in mons if not (site_nb[m] & st.monomers)]
    shared = [m for m in mons if site_nb[m] & st.monomers]
    print(f"p={len(mons)/st.n:.4f}  monomers={len(mons)}  frustrated_check={n_ft}  "
          f"edge_pairs={edge_pairs}  site_pairs={site_pairs}  isolated={len(iso)}  "
          f"sharing={len(shared)}", flush=True)

    for D3 in D3_GRID:
        e, states = sa_best(st.lat, (0.0, 0.0, D3), seed=int(D3 * 1000) + seed)
        n3 = int((states == 2).sum())

        def resolved(ms):
            return sum(1 for m in ms if any(states[s] == 2 for s in st.tri[m]["sites"]))

        print(f"D3={D3:.2f}  E={e:.1f}  n_state3_sites={n3}  "
              f"resolved_isolated={resolved(iso)}/{len(iso)}  "
              f"resolved_sharing={resolved(shared)}/{len(shared)}", flush=True)


if __name__ == "__main__":
    main()

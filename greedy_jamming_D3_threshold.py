"""Same D3-threshold measurement as isolated_monomer_D3_threshold.py, but
for a configuration built via greedy_random_matching_density.py's method
(random greedy maximal matching from a fully-frustrated start) instead
of the independent-set-first + max-flow method. That method's natural
jamming density (~0.12, always succeeds by construction, no retries
needed) gives a second, independent (p, D_c) data point -- lets us see
whether D_c(p) is still climbing from ~1.9 (at p~0.08) toward 3, or has
already saturated near 2.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, simulated_annealing, frustrated_triangles
from monomer_dimer_matching_feasibility import enumerate_all_triangles
from greedy_random_matching_density import build_edge_list, greedy_maximal_matching

NX, NY = 16, 16
D3_GRID = [1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.3]
SEED = 5


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
    lat = RhombileLattice(NX, NY, j12=-1.0, j13=-1.0, j23=-1.0)
    triangles = enumerate_all_triangles(lat)
    n_tri = len(triangles)
    edges = build_edge_list(triangles)

    rng = np.random.default_rng(SEED)
    monomers, chosen_bonds = greedy_maximal_matching(edges, n_tri, rng)
    for b in chosen_bonds:
        b["J"] = 0.0
    ft = frustrated_triangles(lat)
    print(f"N_triangles={n_tri}  monomers={len(monomers)}  "
          f"density={len(monomers)/n_tri:.4f}  frustrated_triangles_check={len(ft)}", flush=True)
    assert len(ft) == len(monomers), "construction mismatch"

    mono_sites = [triangles[m]["sites"] for m in monomers]
    all_mono_sites = [s for tri_sites in mono_sites for s in tri_sites]
    site_use_count = {}
    for s in all_mono_sites:
        site_use_count[s] = site_use_count.get(s, 0) + 1
    shared_sites = {s: c for s, c in site_use_count.items() if c > 1}
    print(f"Sites shared between >=2 monomer triangles: {len(shared_sites)} "
          f"(out of {len(monomers)} monomers, 3 corners each)", flush=True)

    for D3 in D3_GRID:
        D = (0.0, 0.0, D3)
        e, states = sa_best(lat, D, seed=int(D3 * 1000) + 1)
        n_resolved = sum(1 for tri_sites in mono_sites if any(states[s] == 2 for s in tri_sites))
        print(f"D3={D3:.2f}  E={e:.1f}  monomers_resolved_via_state3={n_resolved}/{len(monomers)}",
              flush=True)


if __name__ == "__main__":
    main()

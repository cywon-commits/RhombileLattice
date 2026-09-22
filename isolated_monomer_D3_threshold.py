"""Builds an actual lattice realizing a successful p~0.08 isolated-monomer
configuration (monomer_dimer_matching_feasibility.py's existence check,
now with the real matching extracted and materialized as bond states),
then sweeps D3 to find the empirical threshold -- testing the prediction
from Part III's greedy-heuristic section: a TRULY isolated single defect
(no neighboring conflict to domino-share with) should cross over at
D3*=1 (cost D3 to resolve via state3, vs cost 1 to leave one bond
violated), distinctly different from both the string-interior domino
value D3*=2 and the fully-frustrated bulk value D3*=3.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, simulated_annealing, frustrated_triangles
from monomer_dimer_matching_feasibility import (
    enumerate_all_triangles, build_triangle_adjacency, bipartition,
    greedy_balanced_independent_set, perfect_matching_exists,
)

NX, NY = 16, 16
TARGET_P = 0.08
D3_GRID = [0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4]
SEED = 3


def find_successful_config(adj, color, n_tri, target_p, rng, max_tries=200):
    target_k = int(round(target_p * n_tri))
    target_k -= target_k % 2
    for _ in range(max_tries):
        monomers, k0, k1 = greedy_balanced_independent_set(adj, color, target_k, rng)
        ok, n_left, n_right, matching = perfect_matching_exists(adj, color, monomers, return_matching=True)
        if ok:
            return monomers, matching
    raise RuntimeError("no successful configuration found in max_tries")


def shared_bond(tri_a, tri_b):
    for key in ("r1r2", "r1r3", "r2r3"):
        for key2 in ("r1r2", "r1r3", "r2r3"):
            if id(tri_a[key]) == id(tri_b[key2]):
                return tri_a[key]
    raise ValueError("no shared bond between matched triangles")


def build_realized_lattice(monomers, matching):
    # Fresh lattice -> fresh bond objects; enumerate_all_triangles is
    # deterministic in NX,NY alone, so the SAME triangle index still
    # refers to the same geometric triangle here as it did for the
    # lattice `monomers`/`matching` were computed on -- but the bond
    # dicts themselves must come from THIS lattice, not the old one.
    lat = RhombileLattice(NX, NY, j12=-1.0, j13=-1.0, j23=-1.0)  # every edge active by default
    triangles = enumerate_all_triangles(lat)
    for v, u in matching:
        b = shared_bond(triangles[v], triangles[u])
        b["J"] = 0.0  # this pair's shared edge is their matched (inactive) dimer bond
    ft = frustrated_triangles(lat)
    ft_site_sets = [frozenset((t["r1"], t["a"], t["b"])) for t in ft]
    mono_site_sets = [frozenset(triangles[m]["sites"]) for m in monomers]
    matched_ok = set(ft_site_sets) == set(mono_site_sets)
    return lat, triangles, len(ft), matched_ok


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
    lat0 = RhombileLattice(NX, NY)
    triangles = enumerate_all_triangles(lat0)
    n_tri = len(triangles)
    adj = build_triangle_adjacency(triangles)
    color = bipartition(adj)

    rng = np.random.default_rng(SEED)
    monomers, matching = find_successful_config(adj, color, n_tri, TARGET_P, rng)
    print(f"Found successful config: {len(monomers)} monomers "
          f"(density {len(monomers)/n_tri:.4f}), {len(matching)} matched pairs", flush=True)

    lat, triangles, n_ft, matched_ok = build_realized_lattice(monomers, matching)
    print(f"Realized lattice: frustrated_triangles count={n_ft} "
          f"(expected {len(monomers)}), sites match exactly: {matched_ok}", flush=True)
    if not matched_ok or n_ft != len(monomers):
        print("WARNING: construction mismatch, aborting D3 sweep.", flush=True)
        return

    mono_sites = [triangles[m]["sites"] for m in monomers]
    all_mono_hub_sites = [s for tri_sites in mono_sites for s in tri_sites]
    site_use_count = {}
    for s in all_mono_hub_sites:
        site_use_count[s] = site_use_count.get(s, 0) + 1
    shared_sites = {s: c for s, c in site_use_count.items() if c > 1}
    print(f"Sites shared between >=2 monomer triangles (not necessarily edge-adjacent): "
          f"{len(shared_sites)} (out of {len(monomers)} monomers, 3 corners each)", flush=True)

    for D3 in D3_GRID:
        D = (0.0, 0.0, D3)
        e, states = sa_best(lat, D, seed=int(D3 * 1000) + 1)
        n_resolved = sum(1 for tri_sites in mono_sites if any(states[s] == 2 for s in tri_sites))
        print(f"D3={D3:.2f}  E={e:.1f}  monomers_resolved_via_state3={n_resolved}/{len(monomers)}",
              flush=True)


if __name__ == "__main__":
    main()

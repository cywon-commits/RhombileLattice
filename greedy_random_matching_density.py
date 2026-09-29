"""The complementary construction the user proposed: instead of starting
EMPTY (independent set of monomers + max-flow matching of the rest,
monomer_dimer_matching_feasibility.py), start FULLY FROSTRATED (every
triangle a monomer) and randomly "cut" links to form rhombi one at a
time -- exactly the classical random GREEDY MAXIMAL MATCHING process:
shuffle all edges, and for each edge in turn, cut it (forming a dimer)
only if BOTH endpoint triangles are still untouched; once cut, both
triangles (and hence their OTHER edges) are permanently protected from
being cut again, so no triangle is ever claimed by two different dimers
-- exactly the "protect the neighboring links" rule the user described.

This is a genuinely different stochastic process from the independent-
set-based one (greedy/local, no global optimization), so its natural
JAMMING density (the leftover monomer fraction once no more edges can
be cut) is an independent empirical answer to the same question, not
guaranteed to match the ~7-8% ceiling found via careful independent-set
construction -- random greedy matching is known in general to be
suboptimal (leaves more unmatched vertices than a maximum matching
would), so this is expected to land at an EQUAL OR HIGHER monomer
density than that ceiling.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, frustrated_triangles
from monomer_dimer_matching_feasibility import enumerate_all_triangles

NX, NY = 16, 16
N_TRIALS = 10
SEED = 5


def build_edge_list(triangles):
    """Every (t1, t2, bond) triple, one per shared edge -- same grouping
    as build_triangle_adjacency but keeping the actual bond object and
    without deduplicating into an adjacency set."""
    edges = []
    for key in ("r1r2", "r1r3", "r2r3"):
        groups = {}
        for ti, tri in enumerate(triangles):
            groups.setdefault(id(tri[key]), []).append((ti, tri[key]))
        for entries in groups.values():
            if len(entries) == 2:
                (t1, b1), (t2, _) = entries
                edges.append((t1, t2, b1))
    return edges


def greedy_maximal_matching(edges, n_tri, rng):
    order = list(range(len(edges)))
    rng.shuffle(order)
    matched = [False] * n_tri
    chosen_bonds = []
    for ei in order:
        t1, t2, bond = edges[ei]
        if not matched[t1] and not matched[t2]:
            matched[t1] = matched[t2] = True
            chosen_bonds.append(bond)
    monomers = [t for t in range(n_tri) if not matched[t]]
    return monomers, chosen_bonds


def main():
    lat0 = RhombileLattice(NX, NY)
    triangles0 = enumerate_all_triangles(lat0)
    n_tri = len(triangles0)
    edges0 = build_edge_list(triangles0)
    print(f"N_triangles={n_tri}  N_edges={len(edges0)}", flush=True)

    rng = np.random.default_rng(SEED)
    densities = []
    for trial in range(N_TRIALS):
        lat = RhombileLattice(NX, NY, j12=-1.0, j13=-1.0, j23=-1.0)
        triangles = enumerate_all_triangles(lat)
        edges = build_edge_list(triangles)
        monomers, chosen_bonds = greedy_maximal_matching(edges, n_tri, rng)
        for b in chosen_bonds:
            b["J"] = 0.0
        ft = frustrated_triangles(lat)
        density = len(monomers) / n_tri
        densities.append(density)
        matches_expected = len(ft) == len(monomers)
        print(f"trial {trial}: monomers={len(monomers)}  density={density:.4f}  "
              f"frustrated_triangles_check={len(ft)} (expected {len(monomers)}, "
              f"match={matches_expected})", flush=True)

    print(f"\nMean jamming density over {N_TRIALS} trials: "
          f"{np.mean(densities):.4f} +- {np.std(densities, ddof=1):.4f}", flush=True)


if __name__ == "__main__":
    main()

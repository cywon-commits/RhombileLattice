"""Decisive follow-up to triangular_lattice_D3_2p75_reverify.py.

That heavy re-verification found D3=2.75's true ground state is EXACTLY
the pure proper-3-coloring (E=396.0=144*2.75, zero violated bonds) --
much better than the original light scan's reported 'bump' (E implied
higher, n3=0.322/nv=0.116). This raises a bigger question: is the whole
'broad D3~2-3.5 crossover' picture in the artifact itself just an SA
under-search artifact, with the TRUE transition being a single SHARP
level-crossing?

Exact counting gives two trivially-known candidate ground states:
  (a) proper 3-coloring: n_state3=N/3=144, n_violated=0 -> E_a(D3)=144*D3
  (b) Wannier's pure 2-coloring: n_state3=0, n_violated=N_bonds/3=432
      (exact) -> E_b=432, independent of D3
These two lines cross EXACTLY at D3*=432/144=3.0 -- notably equal to
N_bonds/N = avg_degree/2 = 3 for this coordination-6 lattice (every site
gained coordination 6 once j23 activated every triangle), a DIFFERENT
mechanism from the dilute limit's D3=2 (which comes from 1D domino
economy-of-scale, not a bulk coordination-number count).

If v(k) (violated bonds when k state-3 sites are kept, out of 144) is
EXACTLY LINEAR in k (each removed state-3 site costs a flat +3 violated
bonds, no economy of scale from removing more), the true ground-state
energy is E*(D3)=min(144*D3, 432) at EVERY D3 -- a single sharp kink at
D3=3, no broadening at all, and every D3<3 point's true optimum is
simply the pure coloring untouched.

This script heavy-reverifies D3=2.0, 2.25, 2.5, 3.0, 3.25 (the whole
region the original scan reported as 'gradually departing') against this
exact two-line baseline, seeding SA from BOTH reference states plus many
random restarts, to see whether every point matches the trivial baseline
exactly (confirming a sharp, not broad, transition) or whether any point
beats it (confirming genuine intermediate physics after all).
"""
import numpy as np

from rhombile_lattice import RhombileLattice, simulated_annealing

NX, NY = 12, 12
D3_CHECK = [2.0, 2.25, 2.5, 3.0, 3.25]
N_RANDOM = 20
N_SWEEPS_RANDOM = 5000
N_SEEDED_PER_REF = 10
N_SWEEPS_SEEDED = 4000


def get_pure_coloring_state(lat, seed=1):
    """D3=0's own energy landscape is exactly flat (any valid coloring
    costs nothing), so light SA trivially lands on a proper 3-coloring."""
    rng = np.random.default_rng(seed)
    states, en = simulated_annealing(lat, (0.0, 0.0, 0.0), rng, n_sweeps=1500,
                                      T_start=3.0, T_end=1e-6, record_energy=True)
    assert en[-1] == 0.0, f"expected E=0 pure coloring, got {en[-1]}"
    n3 = int((states == 2).sum())
    assert n3 == lat.n_sites // 3, f"expected exactly N/3 state3 sites, got {n3}"
    return states


def get_wannier_state(lat, seed=2, n_random=15, n_sweeps=4000):
    """Deep in the D3>>3 regime, SA should readily fall into one of
    Wannier's exponentially many degenerate ground states."""
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_random):
        states, en = simulated_annealing(lat, (0.0, 0.0, 6.0), rng, n_sweeps=n_sweeps,
                                          T_start=6.0, T_end=1e-6, record_energy=True)
        if best is None or en[-1] < best[0]:
            best = (en[-1], states)
    n_bonds = len(lat.bonds)
    print(f"Wannier reference: E={best[0]:.1f} (expect {n_bonds/3:.1f}), "
          f"n_state3={(best[1] == 2).sum()} (expect 0)", flush=True)
    return best[1]


def heavy_best(lat, D, seed, seed_states_list=(), n_random=N_RANDOM,
               n_sweeps_random=N_SWEEPS_RANDOM, n_seeded=N_SEEDED_PER_REF,
               n_sweeps_seeded=N_SWEEPS_SEEDED):
    rng = np.random.default_rng(seed)
    candidates = []
    for seed_states in seed_states_list:
        for _ in range(n_seeded):
            states, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_seeded,
                                              T_start=0.8, T_end=1e-6,
                                              states=seed_states.copy(), record_energy=True)
            candidates.append((en[-1], states))
    for _ in range(n_random):
        states, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_random,
                                          T_start=6.0, T_end=1e-6, record_energy=True)
        candidates.append((en[-1], states))
    return min(candidates, key=lambda c: c[0])


def main():
    lat = RhombileLattice(NX, NY, j12=-1.0, j13=-1.0, j23=-1.0)
    N = lat.n_sites
    n_bonds = len(lat.bonds)
    print(f"N={N}  n_bonds={n_bonds}  analytic crossing D3*={n_bonds/N:.4f}\n", flush=True)

    pure_state = get_pure_coloring_state(lat)
    wannier_state = get_wannier_state(lat)
    print(flush=True)

    for D3 in D3_CHECK:
        D = (0.0, 0.0, D3)
        e_trivial = min(N / 3 * D3, n_bonds / 3)
        e, states = heavy_best(lat, D, seed=int(D3 * 1000) + 7,
                                seed_states_list=[pure_state, wannier_state])
        n3 = int((states == 2).sum())
        active = [b for b in lat.bonds if b["J"] != 0.0]
        nv = sum(1 for b in active if states[b["i"]] == states[b["j"]])
        match = "MATCHES trivial baseline" if abs(e - e_trivial) < 1e-6 else \
                (f"BEATS baseline by {e_trivial - e:.2f}" if e < e_trivial - 1e-6
                 else f"WORSE than baseline by {e - e_trivial:.2f} (search failed)")
        print(f"D3={D3:.2f}  E_heavy={e:.1f}  E_trivial_baseline={e_trivial:.1f}  "
              f"n_state3_frac={n3/N:.4f}  n_violated_frac={nv/n_bonds:.4f}  -> {match}",
              flush=True)


if __name__ == "__main__":
    main()

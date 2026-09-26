"""Heavier re-verification of the D3=2.75 anomaly flagged in
triangular_lattice_D3_sweep.py: that scan's chained SA (3 seeded restarts
from the previous D3 point + 5 random restarts, 1500-2500 sweeps) found
n_state3_frac=0.322 / n_violated_frac=0.116 at D3=2.75, a non-monotonic
BUMP relative to its neighbors (D3=2.5: 0.308/0.026; D3=3.0: 0.095/0.238)
-- suspicious because both n_state3 and n_violated should move smoothly
between two points that are each individually well-converged (D3=2.5 sits
on the flat pre-crossover plateau, D3=3.0 is well past it).

This project already caught exactly this failure mode once before (the D2
staircase in the shared artifact's Part VI): a coarse sweep's SA got stuck
in a locally-optimal-looking state that a much heavier search (more random
restarts, more sweeps, AND seeding from neighboring points in both
directions) beats.

Method: build the D3=2.5 and D3=3.0 ground states independently (heavy
search, not reused from the original light chained scan), then attack
D3=2.75 with many random restarts, many sweeps, AND seeded restarts from
both neighbors -- whichever gives the true minimum should not depend on
which neighbor it started from.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, simulated_annealing

NX, NY = 12, 12
N_RANDOM = 25
N_SWEEPS_RANDOM = 5000
N_SEEDED_PER_NEIGHBOR = 10
N_SWEEPS_SEEDED = 4000


def heavy_best(lat, D, seed, seed_states_list=(), n_random=N_RANDOM,
               n_sweeps_random=N_SWEEPS_RANDOM, n_seeded=N_SEEDED_PER_NEIGHBOR,
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


def summarize(lat, D3, e, states):
    N = lat.n_sites
    n3 = int((states == 2).sum())
    active = [b for b in lat.bonds if b["J"] != 0.0]
    nv = sum(1 for b in active if states[b["i"]] == states[b["j"]])
    print(f"D3={D3:.2f}  E={e:.1f}  n_state3_frac={n3/N:.4f}  "
          f"n_violated_frac={nv/len(lat.bonds):.4f}", flush=True)
    return n3 / N, nv / len(lat.bonds)


def main():
    lat = RhombileLattice(NX, NY, j12=-1.0, j13=-1.0, j23=-1.0)

    print("Re-establishing neighbor ground states independently (heavy search):", flush=True)
    e25, s25 = heavy_best(lat, (0.0, 0.0, 2.5), seed=2501, n_random=10, n_sweeps_random=3000)
    summarize(lat, 2.5, e25, s25)
    e30, s30 = heavy_best(lat, (0.0, 0.0, 3.0), seed=3001, n_random=10, n_sweeps_random=3000)
    summarize(lat, 3.0, e30, s30)

    print("\nAttacking D3=2.75 with heavy search (random + both-neighbor-seeded):", flush=True)
    e275, s275 = heavy_best(lat, (0.0, 0.0, 2.75), seed=2751,
                             seed_states_list=[s25, s30])
    summarize(lat, 2.75, e275, s275)

    print("\nComparison against the original light scan's flagged values:")
    print("  original: D3=2.5 n3=0.308/nv=0.026, D3=2.75 n3=0.322/nv=0.116 (BUMP), D3=3.0 n3=0.095/nv=0.238")
    print(f"  heavy:    D3=2.5 n3={int((s25==2).sum())/lat.n_sites:.4f}, "
          f"D3=2.75 n3={int((s275==2).sum())/lat.n_sites:.4f}, "
          f"D3=3.0 n3={int((s30==2).sum())/lat.n_sites:.4f}")


if __name__ == "__main__":
    main()

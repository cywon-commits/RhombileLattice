"""Answers a question raised in conversation: the E(D3)=min(D3*L/2,L) kink
at exactly D3=2 comes from the DILUTE limit (isolated string defects,
where 'domino' resolves 2 conflicts per D3-cost site, isolated defects
resolve 1) -- does the same D3=2 threshold survive in the opposite,
FULLY-FRUSTRATED limit, where every triangle is frustrated?

That limit is literally the classical triangular-lattice 3-state Potts
AF (D3=0 endpoint) interpolating to Wannier's triangular Ising AF
(D3->infinity), NOT this project's own diced-lattice D3=0 sector (which
Kotecky-Salas-Sokal proved has genuine long-range order -- the standard
triangular lattice's own 3-coloring problem is instead T=0 CRITICAL,
Nightingale-Schick/Baxter). No new machinery needed: RhombileLattice's
j23 parameter defaults to 0 (diagonal inactive); setting it to -1 turns
EVERY unit triangle frustrated for free.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, simulated_annealing

NX, NY = 12, 12
D3_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0]


def sa_best(lat, D, seed, seed_states=None, n_seeded=3, n_sweeps_seeded=1500,
            n_random=5, n_sweeps_random=2500):
    rng = np.random.default_rng(seed)
    candidates = []
    if seed_states is not None:
        for _ in range(n_seeded):
            states, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_seeded,
                                              T_start=0.6, T_end=1e-6,
                                              states=seed_states, record_energy=True)
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

    current = None
    for D3 in D3_GRID:
        D = (0.0, 0.0, D3)
        e, states = sa_best(lat, D, seed=int(D3 * 1000) + 1, seed_states=current)
        current = states.copy()
        n3 = int((states == 2).sum())
        active = [b for b in lat.bonds if b["J"] != 0.0]
        nv = sum(1 for b in active if states[b["i"]] == states[b["j"]])
        print(f"D3={D3:.2f}  E={e:.1f}  n_state3_frac={n3/N:.4f}  "
              f"n_violated_frac={nv/n_bonds:.4f}", flush=True)


if __name__ == "__main__":
    main()

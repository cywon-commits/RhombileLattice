"""Gradual D2 sweep for Case2, D1=0, D3=1 fixed, D2 going from 0 down
to -1.2 in small steps (per the user's request: increase |D2| slowly,
since D2 magnitude approaching/exceeding D3=1 is expected to be a
different regime -- see D2_sweep_case2.py's boundary-erosion
hypothesis). Uses continuation (each D2 warm-started from the previous
one's best state) so the run is efficient and the sequence of pictures
shows a smooth evolution, with periodic fresh restarts to avoid
getting stuck.

Saves the full state array per D2 point (pickled) so a separate
plotting pass can make one figure per D2 (or a grid of them) without
re-running any SA.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing, total_energy
from case2_sa_check import build_case2

D3 = 1.0
D2_GRID = [0.0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -1.1, -1.2]
RESULT_PATH = "D2_sweep_gradual_results.pkl"


def sweep(D2_grid=D2_GRID, seed=11, n_warm_sweeps=1500, n_fresh_restarts=2, n_fresh_sweeps=2500):
    with open("case2_sa_result.pkl", "rb") as f:
        current_states = pickle.load(f)["best_states"]

    lat, touched, corners, (i_left, i_right) = build_case2()
    rng = np.random.default_rng(seed)

    results = {}
    for D2 in D2_grid:
        D = (0.0, D2, D3)
        candidates = []
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_warm_sweeps,
                                                T_start=0.8, T_end=1e-6,
                                                states=current_states, record_energy=True)
        candidates.append((energies[-1], states))
        for _ in range(n_fresh_restarts):
            states_r, energies_r = simulated_annealing(lat, D, rng, n_sweeps=n_fresh_sweeps,
                                                         T_start=5.0, T_end=1e-6, record_energy=True)
            candidates.append((energies_r[-1], states_r))
        best_e, best_states = min(candidates, key=lambda c: c[0])
        current_states = best_states.copy()

        n3 = int((best_states == 2).sum())
        active = [b for b in lat.bonds if b["J"] != 0.0]
        violated = [b for b in active if best_states[b["i"]] == best_states[b["j"]]]
        print(f"D2={D2:+.2f}  E={best_e:.3f}  n_state3={n3}  n_violated={len(violated)}", flush=True)
        results[D2] = {"E": best_e, "states": best_states.copy(), "n3": n3, "n_violated": len(violated)}

        with open(RESULT_PATH, "wb") as f:
            pickle.dump(results, f)

    return results


if __name__ == "__main__":
    sweep()

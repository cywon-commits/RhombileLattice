"""Fine-grained D2 sweep in [0, 0.4] to resolve a question raised after
the fact: the coarse grid in D2_sweep_positive.py (0.0, 0.2, 0.4, ...)
showed n_state3 = 3, 11, 15 -- looks like a gradual ramp, but with only 3
points in this range it's impossible to tell whether this is a genuinely
smooth/continuous crossover or a "devil's staircase" of many small
discrete jumps (individual hub sites flipping state2<->state3 at their
own, slightly different, effective thresholds) that only LOOKS smooth at
coarse resolution.

Physical expectation (stated before running): since each hub site's own
choice between state2 (cost D2) and state3 (cost D3=1, if locally useful)
is a discrete combinatorial comparison, no single site's own preference
can change continuously -- so if n_state3 does increase smoothly with D2,
it must be because many sites' individual (slightly different) thresholds
are densely packed across this range, not because any single site's
contribution is itself continuous. A fine grid should reveal this as a
staircase of small, mostly single-site jumps rather than a smooth ramp.

Same continuation-method SA as D2_sweep_positive.py, just at 0.025
resolution across [0, 0.4] instead of 0.2.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing
from case2_sa_check import build_case2

D3 = 1.0
D2_GRID = [round(0.025 * k, 3) for k in range(17)]  # 0.0, 0.025, ..., 0.4
RESULT_PATH = "D2_sweep_fine_positive_results.pkl"


def sweep(D2_grid=D2_GRID, seed=17, n_warm_sweeps=1800, n_fresh_restarts=3, n_fresh_sweeps=2800):
    with open("case2_sa_result.pkl", "rb") as f:
        current_states = pickle.load(f)["best_states"]

    lat, touched, corners, (i_left, i_right) = build_case2()
    rng = np.random.default_rng(seed)

    results = {}
    for D2 in D2_grid:
        D = (0.0, D2, D3)
        candidates = []
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_warm_sweeps,
                                                T_start=0.6, T_end=1e-6,
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
        print(f"D2={D2:+.3f}  E={best_e:.3f}  n_state3={n3}  n_violated={len(violated)}", flush=True)
        results[D2] = {"E": best_e, "states": best_states.copy(), "n3": n3, "n_violated": len(violated)}

        with open(RESULT_PATH, "wb") as f:
            pickle.dump(results, f)

    return results


if __name__ == "__main__":
    sweep()

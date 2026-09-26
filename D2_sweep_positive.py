"""D2>0 counterpart to D2_sweep_gradual.py. Prediction made before running
(shared with the user): this should NOT mirror the D2<0 three-stage
collapse. For D2<0, D_rim_default=min(D1,D2)=D2 keeps decreasing, which
is an entropic/opportunity-cost trade-off against a FIXED +1 bond-
violation cost -- hence the gradual, three-stage erosion seen there.

For D2>0, D_rim_default=min(D1,D2)=D1=0 stays FIXED (D1 is always
cheaper), so rim's role and the local defect patch's threshold
(D3 < 2 + D_rim_default = 2) never move. Instead, the relevant
comparison is hub's OWN choice of which state to use for its "pay
extra" role: state2 (cost D2) vs state3 (cost D3=1). This is a simple
threshold comparison, not an opportunity-cost-vs-entropy trade-off, so
the predicted transition is SHARP and located exactly at D2=D3=1 --
and much larger in scale (potentially all ~256 hub sites currently
paying D2, not just the ~33-bond loop), not a gradual 3-stage erosion.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing
from case2_sa_check import build_case2

D3 = 1.0
D2_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]
RESULT_PATH = "D2_sweep_positive_results.pkl"


def sweep(D2_grid=D2_GRID, seed=13, n_warm_sweeps=1500, n_fresh_restarts=2, n_fresh_sweeps=2500):
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

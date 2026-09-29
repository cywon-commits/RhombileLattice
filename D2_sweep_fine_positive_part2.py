"""Continues D2_sweep_fine_positive.py's re-verification from D2=0.4 up
to D2=1.0 (step 0.05), since that script's own D2=0.4 endpoint already
beat the original D2_sweep_positive.py's result there (E=118.8 vs 119.6),
suggesting the originally-reported "gradual precursor" (n_state3 growing
3->17 through D2=0.9) may have been substantially an SA-optimization
artifact rather than the true ground-state trend. Uses MORE SA effort
than either earlier script (more restarts, more sweeps) specifically
because we now know this problem needs it to avoid getting stuck.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing
from case2_sa_check import build_case2

D3 = 1.0
D2_GRID = [round(0.4 + 0.05 * k, 3) for k in range(13)]  # 0.4, 0.45, ..., 1.0
RESULT_PATH = "D2_sweep_fine_positive_part2_results.pkl"


def sweep(D2_grid=D2_GRID, seed=29, n_warm_sweeps=2200, n_fresh_restarts=5, n_fresh_sweeps=3200):
    with open("D2_sweep_fine_positive_results.pkl", "rb") as f:
        fine1 = pickle.load(f)
    current_states = fine1[0.4]["states"]

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

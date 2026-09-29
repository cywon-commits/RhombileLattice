"""Re-verifies the D2<0 branch (originally D2_sweep_gradual.py, coarse
step 0.1) at fine resolution (step 0.05) with the SAME extra SA rigor
that revealed the D2>0 branch's "gradual precursor" was largely an
under-optimized-SA artifact (the true D2>0 trend is a long flat plateau
then a sharp jump, not a smooth ramp). Original D2<0 claim: stable
through D2=-0.5, a "boundary erosion" precursor growing 3->16 from
D2~-0.6 to -0.9, then full collapse (0 state-3, all 33 bonds violated)
by D2<=-1.1. Checking whether that precursor stage is real or is ALSO
mostly a flat plateau with a sharp jump once given enough SA effort.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing
from case2_sa_check import build_case2

D3 = 1.0
D2_GRID = [round(-0.05 * k, 3) for k in range(25)]  # 0.0, -0.05, ..., -1.2
RESULT_PATH = "D2_sweep_fine_negative_results.pkl"


def sweep(D2_grid=D2_GRID, seed=41, n_warm_sweeps=2200, n_fresh_restarts=5, n_fresh_sweeps=3200):
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

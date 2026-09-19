"""Full D3 sweep for Case 2 (5-segment detour), using the same
SA-based approach that found case2_sa_check.py's D3=1 result (E=3,
Case0-style bulk swap + a local 3-site patch on the SHORT direct path
between the two defects) -- neither hub-fixed (energy_via_mincut) nor
the naive free-2-coloring are trustworthy for a bent/branching path
like this one, so this sweeps by simulated annealing instead.

Uses a continuation/adiabatic strategy rather than independent restarts
at every D3: starts at D3=0 (seeded from the free-coloring, which is
already known-good there since state-3 is free) and walks D3 upward in
small steps, warm-starting each step from the previous step's best
state (physically sensible -- the optimum should move smoothly with
D3 outside of a genuine phase-transition-like jump) plus a periodic
fresh random restart to catch a different basin. Results are appended
to case2_d3_sweep_results.csv after every D3 point, and the driver
script (run_and_commit.sh-equivalent logic is left to the caller) is
expected to git-commit that file incrementally rather than waiting for
the whole sweep -- see the __main__ block.
"""
import csv
import os

import numpy as np

from rhombile_lattice import total_energy, simulated_annealing
from closed_loop_demo import CENTER, _bfs_two_coloring
from case2_sa_check import build_case2

CSV_PATH = "case2_d3_sweep_results.csv"
D3_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
           2.75, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]


def sweep(n_refine_sweeps=1200, n_fresh_restarts=2, n_fresh_sweeps=2500, seed=0):
    lat, touched, corners, (i_left, i_right) = build_case2()
    pos, _ = lat.site_positions()
    start_site = int(np.argmax(np.linalg.norm(pos - CENTER, axis=1)))
    free_coloring = _bfs_two_coloring(lat, start_site)

    rng = np.random.default_rng(seed)
    current_states = free_coloring.copy()

    write_header = not os.path.exists(CSV_PATH)
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["D3", "best_E", "n_state3", "n_violated"])

        for D3 in D3_GRID:
            D = (0.0, 0.0, D3)
            candidates = []

            # warm start: refine from the previous D3's solution
            states, energies = simulated_annealing(
                lat, D, rng, n_sweeps=n_refine_sweeps, T_start=0.8, T_end=1e-5,
                states=current_states, record_energy=True)
            candidates.append((energies[-1], states))

            # a couple of fresh restarts to catch a different basin
            for _ in range(n_fresh_restarts):
                states_r, energies_r = simulated_annealing(
                    lat, D, rng, n_sweeps=n_fresh_sweeps, T_start=5.0, T_end=1e-5,
                    record_energy=True)
                candidates.append((energies_r[-1], states_r))

            best_e, best_states = min(candidates, key=lambda c: c[0])
            current_states = best_states.copy()
            n_state3 = int((best_states == 2).sum())
            active = [b for b in lat.bonds if b["J"] != 0.0]
            n_violated = sum(1 for b in active if best_states[b["i"]] == best_states[b["j"]])

            print(f"D3={D3:<6} best_E={best_e:<7} n_state3={n_state3:<4} "
                  f"n_violated={n_violated}", flush=True)
            writer.writerow([D3, best_e, n_state3, n_violated])
            f.flush()


if __name__ == "__main__":
    sweep()

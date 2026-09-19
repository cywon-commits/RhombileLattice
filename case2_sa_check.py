"""SA cross-check for Case 2 (5-segment detour, hexagon's other 5 sides)
at D3=1: does the true low-energy state beat the naive free-2-coloring
(E=17, 0 state-3) by paying with a few state-3 sites instead, closer to
Case 1's own D3=1 answer (E=4) for the identical two defect points?

Early result (partial run, see git history/commit message for the run
this file accompanies): seeding SA from the free-2-coloring and cooling
further at D3=1 repeatedly converged to E=3.0 -- BETTER than Case 1's
E=4 for the same two defects, and much better than the naive E=17.
That's notable on its own: it means the long way around can, after
state reorganization, end up net *cheaper* than the direct edge at this
D3, presumably because part of the interior swap (free, needed anyway
per Case 0) also helps pay down some of what would otherwise be a pure
defect-pair cost in Case 1. Not yet proven optimal, but consistent
across multiple independent SA restarts from the same seed.
"""
import pickle

import numpy as np

from rhombile_lattice import (
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
    total_energy, simulated_annealing,
)
from dual_string_demo import build_dual
from closed_loop_demo import hexagon_corners, NX, NY, CENTER, _bfs_two_coloring

RESULT_PATH = "case2_sa_result.pkl"


def build_case2():
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    other_order = ([i_right] + [i for i in range(len(corners)) if i not in (i_left, i_right)]
                   + [i_left])
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    touched = []
    for a, b in zip(other_order[:-1], other_order[1:]):
        _, _, nodes, bonds = route_string_between_points(rho, sib, hop, corners[a], corners[b])
        touched.extend(apply_dual_string_defect(tri, nodes, bonds))
    return lat, touched, corners, (i_left, i_right)


def run(n_seeded=6, n_random=6, n_sweeps_seeded=2500, n_sweeps_random=3000, D3=1.0, seed=0):
    lat, touched, corners, (i_left, i_right) = build_case2()
    pos, sub_of = lat.site_positions()
    start = int(np.argmax(np.linalg.norm(pos - CENTER, axis=1)))
    free_coloring = _bfs_two_coloring(lat, start)
    e_free = total_energy(lat, free_coloring, D=(0.0, 0.0, D3))
    print(f"free-coloring baseline energy at D3={D3}: {e_free}", flush=True)

    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(seed)
    best_E, best_states = np.inf, None

    for trial in range(n_seeded):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_seeded,
                                                 T_start=2.0, T_end=1e-5,
                                                 states=free_coloring, record_energy=True)
        print(f"  seeded-from-free trial {trial}: final E={energies[-1]}", flush=True)
        if energies[-1] < best_E:
            best_E, best_states = energies[-1], states.copy()

    for trial in range(n_random):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_random,
                                                 T_start=5.0, T_end=1e-5, record_energy=True)
        print(f"  random-start trial {trial}: final E={energies[-1]}", flush=True)
        if energies[-1] < best_E:
            best_E, best_states = energies[-1], states.copy()

    print(f"BEST energy found: {best_E}", flush=True)
    n_state3 = int((best_states == 2).sum())
    print(f"n_state3 in best: {n_state3}", flush=True)
    with open(RESULT_PATH, "wb") as f:
        pickle.dump({"best_E": best_E, "best_states": best_states, "e_free": e_free}, f)
    print(f"saved {RESULT_PATH}", flush=True)
    return best_E, best_states


if __name__ == "__main__":
    run()

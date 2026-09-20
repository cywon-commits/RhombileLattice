"""Follow-up to winding_string_experiment.py after noticing its SA run was
under-tuned: it reported E_direct=28 even though full_state_via_mincut
already gives E=25 for the identical (fixed) bond pattern -- meaning the
first run's "best of several SA trials" never even matched its own best
available seed. That's a convergence bug, not a physics result, so it's
fixed here by (1) always including the raw mincut-seed energy itself as a
candidate, (2) a gentler seeded anneal (lower T_start, more sweeps) so it
doesn't get kicked out of a good seed and fail to return, (3) many more
random restarts to genuinely stress-test whether a cheaper state exists.

This also answers the sharper question the first run left open: is
full_state_via_mincut's E=25 for the direct (non-winding, L=25) construction
truly the GLOBAL optimum, or -- like Case0's closed loop and Case2's bent
detour, both of which beat their own naive hub-fixed prediction via a
large-scale r1<->rim relabeling trick -- can a big chunk of the lattice
relabel itself for free and undercut it? hub-fixed/mincut has repeatedly
been *exact* for straight, unbent single strings elsewhere in this project
(Case1, distance_scaling_demo), but "straight yet wrapping through the
periodic seam" was never actually tested against that assumption before.
"""
import pickle

import numpy as np

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, frustrated_triangles,
    simulated_annealing, total_energy, full_state_via_mincut,
)

NX, NY = 30, 8
Y_CROSS = 2.1
P_START = np.array([2.0, Y_CROSS])
P_END = np.array([27.0, Y_CROSS])

D3_GRID = [1.0, 2.0, 3.0]
RESULT_PATH = "winding_string_results2.pkl"


def build(wrap):
    lat = RhombileLattice(NX, NY)
    flipped, p_end_used = apply_string_defect(lat, P_START, P_END, wrap=wrap)
    return lat, flipped, p_end_used


def best_energy(lat, D3, seed, n_random=10, n_sweeps_random=4000,
                 n_seeded_gentle=6, n_sweeps_gentle=4000):
    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(seed)
    candidates = []

    try:
        seed_states, _ = full_state_via_mincut(lat, D3)
        candidates.append((total_energy(lat, seed_states, D), seed_states.copy(), "mincut-seed"))
        for _ in range(n_seeded_gentle):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_gentle,
                                                     T_start=0.5, T_end=1e-6,
                                                     states=seed_states, record_energy=True)
            candidates.append((min(energies), states, "seeded-gentle"))
    except ValueError as exc:
        print(f"  [note] mincut seed unavailable at D3={D3}: {exc}")

    for _ in range(n_random):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_random,
                                                 T_start=5.0, T_end=1e-6, record_energy=True)
        candidates.append((min(energies), states, "random-restart"))

    best_e, best_states, best_src = min(candidates, key=lambda c: c[0])
    return best_e, best_states, best_src


def main():
    lat_d, flipped_d, p_end_d = build(wrap=False)
    lat_w, flipped_w, p_end_w = build(wrap=True)
    print(f"direct: {len(flipped_d)} bonds flipped, endpoint {tuple(np.round(p_end_d,2))}, "
          f"defects={len(frustrated_triangles(lat_d))}")
    print(f"winding: {len(flipped_w)} bonds flipped, endpoint {tuple(np.round(p_end_w,2))}, "
          f"defects={len(frustrated_triangles(lat_w))}\n")

    results = {}
    print(f"{'D3':>5}  {'E_direct':>10}  {'src':>14}  {'E_winding':>10}  {'src':>14}")
    for D3 in D3_GRID:
        e_d, states_d, src_d = best_energy(lat_d, D3, seed=101 + int(D3 * 10))
        e_w, states_w, src_w = best_energy(lat_w, D3, seed=202 + int(D3 * 10))
        print(f"{D3:>5.1f}  {e_d:>10.3f}  {src_d:>14}  {e_w:>10.3f}  {src_w:>14}")
        results[D3] = {"E_direct": e_d, "states_direct": states_d, "src_direct": src_d,
                        "E_winding": e_w, "states_winding": states_w, "src_winding": src_w}
        with open(RESULT_PATH, "wb") as f:
            pickle.dump(results, f)
    print(f"\nsaved {RESULT_PATH}")


if __name__ == "__main__":
    main()

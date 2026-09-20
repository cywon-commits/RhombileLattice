"""Part (b): does a genuinely WINDING connecting path (crosses the torus's
own periodic seam) between two fixed defect sites reach the same ground
state as a NON-winding path between the identical two sites, once SA is
allowed to fully reorganize both -- or does it get stuck in a distinct,
higher-energy flux sector, as the Thurston/dimer-flux literature found for
periodic (as opposed to simply-connected) geometries?

Setup, on the same NX=30, NY=8, Y_CROSS=2.1 lattice already validated in
distance_scaling_demo.py (a single straight apply_string_defect string
reliably confines frustration to exactly its 2 endpoints):

  p_start = (2.0, Y_CROSS)     p_end = (2.0 + 25.0, Y_CROSS) = (27.0, Y_CROSS)

  DIRECT construction (control):  apply_string_defect(..., wrap=False)
    -- treated as a plain segment in the plane; goes straight across the
    middle of the box, length 25, never touches a wrap bond (winding = 0
    by construction -- this is exactly the same regime as every earlier
    test in the project).

  WINDING construction:           apply_string_defect(..., wrap=True)
    -- minimum_image_endpoint sees that the *other* way around (length
    30 - 25 = 5) is shorter than the direct 25, so it automatically
    reroutes the string through the box's periodic seam. Both endpoints
    (p_start and p_end) are the same two physical lattice sites as the
    direct case; only the connecting path's homotopy class differs.

Both constructions are built on independent fresh lattices (no shared
mutable state) and checked for (1) exactly 2 frustrated triangles at
construction time and (2) whether any flipped bond has wraps=True, before
handing them to the same continuation-method SA used for the D2/D3 sweeps
elsewhere in this project.
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
P_END = np.array([27.0, Y_CROSS])  # direct distance 25; wrap distance 30-25=5

D3_GRID = [1.0, 2.0, 3.0]
RESULT_PATH = "winding_string_results.pkl"


def build(wrap):
    lat = RhombileLattice(NX, NY)
    flipped, p_end_used = apply_string_defect(lat, P_START, P_END, wrap=wrap)
    uses_wrap = any(b["wraps"] for b in flipped)
    ft = frustrated_triangles(lat)
    return lat, flipped, p_end_used, uses_wrap, ft


def run_sa(lat, D3, seed, n_seeded=4, n_random=4, n_sweeps_seeded=2000, n_sweeps_random=2500,
           seed_states=None):
    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(seed)
    best_E, best_states = np.inf, None
    if seed_states is not None:
        for _ in range(n_seeded):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_seeded,
                                                     T_start=2.0, T_end=1e-5,
                                                     states=seed_states, record_energy=True)
            if energies[-1] < best_E:
                best_E, best_states = energies[-1], states.copy()
    for _ in range(n_random):
        states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps_random,
                                                 T_start=5.0, T_end=1e-5, record_energy=True)
        if energies[-1] < best_E:
            best_E, best_states = energies[-1], states.copy()
    return best_E, best_states


def main():
    print(f"Direct vs winding string between the SAME two sites {tuple(P_START)} -> "
          f"{tuple(P_END)} on a {NX}x{NY} torus (direct distance 25, wrap distance 5).\n")

    lat_d, flipped_d, p_end_d, wrap_d, ft_d = build(wrap=False)
    print(f"DIRECT   construction: {len(flipped_d)} bonds flipped, endpoint used={tuple(np.round(p_end_d,2))}, "
          f"uses_periodic_wrap={wrap_d}, frustrated triangles={len(ft_d)}")

    lat_w, flipped_w, p_end_w, wrap_w, ft_w = build(wrap=True)
    print(f"WINDING  construction: {len(flipped_w)} bonds flipped, endpoint used={tuple(np.round(p_end_w,2))}, "
          f"uses_periodic_wrap={wrap_w}, frustrated triangles={len(ft_w)}")

    assert len(ft_d) == 2, f"direct construction should have exactly 2 defects, got {len(ft_d)}"
    assert len(ft_w) == 2, f"winding construction should have exactly 2 defects, got {len(ft_w)}"
    assert not wrap_d, "direct construction should never touch a wrap bond"
    assert wrap_w, "winding construction should use the periodic seam"
    print("\n(sanity checks passed: both constructions isolate exactly 2 real defects; "
          "only the winding one actually crosses the seam)\n")

    results = {"setup": {"NX": NX, "NY": NY, "p_start": P_START, "p_end": P_END,
                          "flipped_direct": len(flipped_d), "flipped_winding": len(flipped_w)},
               "D3": {}}

    def seed_or_none(lat, D3, tag):
        try:
            states, _ = full_state_via_mincut(lat, D3)
            return states
        except ValueError as exc:
            print(f"  [note] mincut seed unavailable for {tag} at D3={D3}: {exc}")
            return None

    print(f"{'D3':>5}  {'E_direct':>10}  {'n3_direct':>10}  {'nv_direct':>10}   "
          f"{'E_winding':>10}  {'n3_winding':>11}  {'nv_winding':>11}   match?")
    for D3 in D3_GRID:
        seed_states_d = seed_or_none(lat_d, D3, "direct")
        e_d, states_d = run_sa(lat_d, D3, seed=1, seed_states=seed_states_d)
        active = [b for b in lat_d.bonds if b["J"] != 0.0]
        nv_d = sum(1 for b in active if states_d[b["i"]] == states_d[b["j"]])
        n3_d = int((states_d == 2).sum())

        seed_states_w = seed_or_none(lat_w, D3, "winding")
        e_w, states_w = run_sa(lat_w, D3, seed=2, seed_states=seed_states_w)
        active_w = [b for b in lat_w.bonds if b["J"] != 0.0]
        nv_w = sum(1 for b in active_w if states_w[b["i"]] == states_w[b["j"]])
        n3_w = int((states_w == 2).sum())

        match = "YES" if abs(e_d - e_w) < 1e-6 else "NO"
        print(f"{D3:>5.1f}  {e_d:>10.3f}  {n3_d:>10d}  {nv_d:>10d}   "
              f"{e_w:>10.3f}  {n3_w:>11d}  {nv_w:>11d}   {match}")

        results["D3"][D3] = {
            "E_direct": e_d, "n3_direct": n3_d, "nv_direct": nv_d, "states_direct": states_d,
            "E_winding": e_w, "n3_winding": n3_w, "nv_winding": nv_w, "states_winding": states_w,
        }
        with open(RESULT_PATH, "wb") as f:
            pickle.dump(results, f)

    print(f"\nsaved {RESULT_PATH}")


if __name__ == "__main__":
    main()

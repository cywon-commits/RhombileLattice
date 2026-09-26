"""Properly-controlled redesign of D2_threshold_vs_area.py: the two real
defects (p_left, p_right) are FIXED (same direct separation for every
run), so Case1 (the direct edge between them) is IDENTICAL across all
sizes/shapes. Only Case2's detour changes: instead of a hexagon, it is a
2-segment path p_right -> p_top -> p_left through one adjustable apex
point. Enclosed area = 0.5 * base * height exactly (base fixed), so
sweeping the apex's HEIGHT cleanly varies area alone. Sweeping the
apex's HORIZONTAL offset at FIXED height keeps the area exactly fixed
while changing the triangle's shape (a skewed vs. symmetric triangle) --
a clean curvature/shape test with no area confound.
"""
import pickle

import numpy as np

from rhombile_lattice import (
    route_string_between_points, apply_dual_string_defect, simulated_annealing,
    frustrated_triangles,
)
from dual_string_demo import build_dual

D3 = 1.0
NX, NY = 20, 24
P_LEFT = np.array([6.25, 3.0])
P_RIGHT = np.array([9.75, 3.0])
BASE = np.linalg.norm(P_RIGHT - P_LEFT)


def build_triangle_case(p_top):
    """Case1: direct edge P_LEFT-P_RIGHT. Case2: P_RIGHT -> p_top -> P_LEFT."""
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    touched = []
    for a, b in [(P_RIGHT, p_top), (p_top, P_LEFT)]:
        _, _, nodes, bonds = route_string_between_points(rho, sib, hop, a, b)
        touched.extend(apply_dual_string_defect(tri, nodes, bonds))
    return lat, touched


def triangle_area(p_top):
    return 0.5 * BASE * abs(p_top[1] - P_LEFT[1])


def sa_best(lat, D, seed, seed_states=None, n_seeded=4, n_sweeps_seeded=2200,
            n_random=5, n_sweeps_random=3000):
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
                                          T_start=5.0, T_end=1e-6, record_energy=True)
        candidates.append((en[-1], states))
    return min(candidates, key=lambda c: c[0])


def scan(label, p_top, D2_grid, seed_base):
    area = triangle_area(p_top)
    lat, touched = build_triangle_case(p_top)
    ft = frustrated_triangles(lat)
    print(f"\n=== {label}: p_top={tuple(np.round(p_top,2))}, area={area:.2f}, "
          f"{len(touched)} bonds touched, {len(ft)} defects ===")
    assert len(ft) == 2, f"expected 2 defects, got {len(ft)}"

    current_states = None
    results = {"area": area, "p_top": p_top, "D2": {}}
    for D2 in D2_grid:
        D = (0.0, D2, D3)
        best_e, best_states = sa_best(lat, D, seed=seed_base + int(D2 * 1000),
                                       seed_states=current_states)
        current_states = best_states.copy()
        n3 = int((best_states == 2).sum())
        active = [b for b in lat.bonds if b["J"] != 0.0]
        nv = sum(1 for b in active if best_states[b["i"]] == best_states[b["j"]])
        print(f"  D2={D2:+.3f}  E={best_e:.3f}  n_state3={n3}  n_violated={nv}", flush=True)
        results["D2"][D2] = {"E": best_e, "n3": n3, "n_violated": nv}
    with open(f"D2_threshold_v2_{label}.pkl", "wb") as f:
        pickle.dump(results, f)
    return results


if __name__ == "__main__":
    D2_grid_default = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    mid_x = (P_LEFT[0] + P_RIGHT[0]) / 2

    # AREA sweep: symmetric apex (offset=0), varying height
    for h, label in [(2.0, "h2"), (4.0, "h4"), (6.0, "h6"), (9.0, "h9"), (13.0, "h13")]:
        p_top = np.array([mid_x, P_LEFT[1] + h])
        scan(f"area_{label}", p_top, D2_grid_default, seed_base=1000 + int(h * 100))

    # CURVATURE sweep: fixed height (area), varying horizontal skew
    h_fixed = 6.0
    for dx, label in [(0.0, "skew0"), (2.0, "skew2"), (4.0, "skew4")]:
        p_top = np.array([mid_x + dx, P_LEFT[1] + h_fixed])
        scan(f"curv_{label}", p_top, D2_grid_default, seed_base=5000 + int(dx * 100))

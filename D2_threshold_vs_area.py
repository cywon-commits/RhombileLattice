"""Tests the area-dependence prediction for the direct<->detour (Case1
vs Case2) D2 threshold, discussed but deferred earlier: if the swap's
cost is D2*(enclosed area) while the local defect-patch benefit is a
fixed O(1) amount (this model's antiphase boundary carries NO perimeter
tension, unlike standard Ising nucleation), the critical D2* where the
detour stops being worth it should scale roughly as 1/Area -- i.e. a
bigger loop should collapse at a SMALLER D2. A same-area but differently
-shaped (elongated) loop is also tested to check the predicted
curvature-independence (no perimeter term to reward compactness here).

Reuses the exact Case2 construction method (route_string_between_points
+ apply_dual_string_defect chained around hexagon_corners), generalized
to take an explicit (nx, ny, center, radius, n_corners) instead of
closed_loop_demo's module-level constants, so different sizes/shapes can
share one code path.
"""
import pickle

import numpy as np

from rhombile_lattice import (
    route_string_between_points, apply_dual_string_defect, simulated_annealing,
    total_energy,
)
from dual_string_demo import build_dual

D3 = 1.0


def polygon_corners(center, rx, ry, n):
    """Regular n-gon by default (rx=ry); an ellipse-like elongated shape
    if rx != ry, for the curvature-independence check."""
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [center + np.array([rx * np.cos(a), ry * np.sin(a)]) for a in angles]


def polygon_area(corners):
    """Shoelace formula."""
    c = np.array(corners)
    x, y = c[:, 0], c[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def build_case(nx, ny, corners):
    """Case1 (direct edge between the two 'bottom' corners) and Case2
    (the other 5 sides, chained) on the same lattice, same convention as
    case2_sa_check.build_case2."""
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    p_left, p_right = corners[i_left], corners[i_right]
    other_order = ([i_right] + [i for i in range(len(corners)) if i not in (i_left, i_right)]
                   + [i_left])

    lat, tri, rho, sib, hop = build_dual(nx, ny)
    touched = []
    for a, b in zip(other_order[:-1], other_order[1:]):
        _, _, nodes, bonds = route_string_between_points(rho, sib, hop, corners[a], corners[b])
        touched.extend(apply_dual_string_defect(tri, nodes, bonds))
    return lat, touched, (p_left, p_right)


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


def scan_size(label, nx, ny, center, rx, ry, D2_grid, seed_base):
    corners = polygon_corners(np.array(center), rx, ry, 6)
    area = polygon_area(corners)
    lat, touched, endpoints = build_case(nx, ny, corners)
    print(f"\n=== {label}: nx={nx},ny={ny}, rx={rx},ry={ry}, area={area:.1f} ===")

    current_states = None
    results = {"area": area, "D2": {}}
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
    with open(f"D2_threshold_{label}.pkl", "wb") as f:
        pickle.dump(results, f)
    return results


if __name__ == "__main__":
    D2_grid_default = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # small hexagon, same 16x16 box as the original (radius=3.5) case
    scan_size("small_r2.0", nx=16, ny=16, center=(6.0, 6.5), rx=2.0, ry=2.0,
              D2_grid=D2_grid_default, seed_base=1000)

    # original size, re-verified for a direct area/threshold data point
    scan_size("medium_r3.5", nx=16, ny=16, center=(6.0, 6.5), rx=3.5, ry=3.5,
              D2_grid=D2_grid_default, seed_base=2000)

    # large hexagon, needs a bigger box
    scan_size("large_r5.5", nx=22, ny=22, center=(9.0, 9.5), rx=5.5, ry=5.5,
              D2_grid=D2_grid_default, seed_base=3000)

    # same area as medium (r=3.5, area ~31.8), but elongated (curvature check)
    # regular hexagon area = (3*sqrt(3)/2)*r^2; for an ellipse-like hexagon
    # with rx*ry = r^2, area scales similarly -- pick rx=5.0 => ry=r^2/rx=2.45
    scan_size("elongated_matched_area", nx=16, ny=16, center=(6.0, 6.5), rx=5.0, ry=2.45,
              D2_grid=D2_grid_default, seed_base=4000)

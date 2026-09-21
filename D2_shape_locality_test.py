"""Tests whether a single detour loop's D2-driven state3 transition is a
genuine global all-or-nothing swap (as it appeared for convex, uniform-
curvature hexagon/triangle loops) or can be spatially graded when the loop
is non-convex (a wide "bulge" next to a narrow "pinch").

Same fixed real-defect pair (P_LEFT, P_RIGHT) as D2_threshold_vs_area_v2.py
so Case1 is identical; Case2 here is a 3-segment path
P_RIGHT -> A (wide bulge) -> B (reflex vertex / pinch) -> P_LEFT.
Convexity is checked explicitly via consecutive edge-vector cross products
(a sign change means a reflex vertex, i.e. non-convex).

For each D2 we record not just n_state3 but the SET of state3 site
positions, so onset order can be checked against distance-to-A (bulge)
vs distance-to-B (pinch).
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
A = np.array([11.5, 9.0])   # wide bulge
B = np.array([8.0, 3.3])    # pinch (close to the P_LEFT-P_RIGHT base line)

WAYPOINTS = [P_RIGHT, A, B, P_LEFT]


def polygon_area_and_convexity(pts):
    pts = np.array(pts)
    n = len(pts)
    area2 = 0.0
    signs = []
    for k in range(n):
        p0, p1, p2 = pts[k - 1], pts[k], pts[(k + 1) % n]
        e1, e2 = p1 - p0, p2 - p1
        cross = e1[0] * e2[1] - e1[1] * e2[0]
        signs.append(np.sign(cross))
        area2 += p0[0] * p1[1] - p1[0] * p0[1]
    is_convex = len(set(s for s in signs if s != 0)) <= 1
    return abs(area2) / 2.0, is_convex, signs


def build_case():
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    touched = []
    for a, b in zip(WAYPOINTS[:-1], WAYPOINTS[1:]):
        _, _, nodes, bonds = route_string_between_points(rho, sib, hop, a, b)
        touched.extend(apply_dual_string_defect(tri, nodes, bonds))
    return lat, touched


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


def main():
    poly_order = [P_LEFT, P_RIGHT, A, B]
    area, is_convex, signs = polygon_area_and_convexity(poly_order)
    print(f"polygon area={area:.3f}, convex={is_convex}, edge-cross signs={signs}")
    assert not is_convex, "expected a reflex vertex (non-convex) for this test"

    lat, touched = build_case()
    ft = frustrated_triangles(lat)
    print(f"{len(touched)} bonds touched, {len(ft)} defects")
    assert len(ft) == 2, f"expected 2 defects, got {len(ft)}"

    pos, sub_of = lat.site_positions()

    D2_grid = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
               0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

    current_states = None
    results = {"area": area, "D2": {}, "pos": pos}
    for D2 in D2_grid:
        D = (0.0, D2, D3)
        best_e, best_states = sa_best(lat, D, seed=7000 + int(D2 * 1000),
                                       seed_states=current_states)
        current_states = best_states.copy()
        n3 = int((best_states == 2).sum())
        state3_sites = np.where(best_states == 2)[0]
        d_to_A = np.linalg.norm(pos[state3_sites] - A, axis=1) if n3 else np.array([])
        d_to_B = np.linalg.norm(pos[state3_sites] - B, axis=1) if n3 else np.array([])
        mean_dA = float(d_to_A.mean()) if n3 else float("nan")
        mean_dB = float(d_to_B.mean()) if n3 else float("nan")
        print(f"  D2={D2:+.3f}  E={best_e:.3f}  n_state3={n3}  "
              f"mean_dist_to_bulgeA={mean_dA:.2f}  mean_dist_to_pinchB={mean_dB:.2f}",
              flush=True)
        results["D2"][D2] = {"E": best_e, "n3": n3, "state3_sites": state3_sites.copy(),
                              "mean_dA": mean_dA, "mean_dB": mean_dB}

    with open("D2_shape_locality_test.pkl", "wb") as f:
        pickle.dump(results, f)
    print("Saved D2_shape_locality_test.pkl")


if __name__ == "__main__":
    main()

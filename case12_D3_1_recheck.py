"""Heavy SA re-check of Case 1 vs Case 2 at D3=1 (after exact_Einf.py showed
Case 1's D3=inf value is 6, not the 10 reported earlier). 12 random restarts
x 4000 sweeps each, plus the exact D3=inf 2-colouring energy as the D3-free
upper bound."""
import numpy as np

from rhombile_lattice import route_string_between_points, apply_dual_string_defect, simulated_annealing
from dual_string_demo import build_dual
from closed_loop_demo import hexagon_corners, NX, NY
from case2_sa_check import build_case2
from exact_Einf import exact_two_state_energy


def case1():
    corners = hexagon_corners()
    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    _, _, nodes, bonds = route_string_between_points(rho, sib, hop, corners[i_left], corners[i_right])
    apply_dual_string_defect(tri, nodes, bonds)
    return lat


for name, lat in [("Case 1", case1()), ("Case 2", build_case2()[0])]:
    rng = np.random.default_rng(7)
    es = []
    for _ in range(12):
        states, en = simulated_annealing(lat, (0.0, 0.0, 1.0), rng, n_sweeps=4000,
                                         T_start=4.0, T_end=1e-6, record_energy=True)
        es.append(en[-1])
    print(f"{name}: D3=1 SA best {min(es):.1f}  (all: {sorted(es)})  exact E_inf {exact_two_state_energy(lat)[1]}", flush=True)

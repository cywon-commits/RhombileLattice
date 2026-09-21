"""Tests a hypothesis raised in conversation: given a genuine ground state
at some D3<2 with state3 sites present, does greedily converting each
state3 site to whichever of state1/state2 has FEWER neighbors among its
own active-bond neighbors (breaking ties arbitrarily) -- thereby
eliminating state3 usage entirely -- automatically reproduce the true
D3=infinity (forced 2-state Ising) ground state?

Locally, for a single site holding its neighbors fixed, "convert to the
minority state" is trivially the cost-minimizing choice (converting to
state k only creates new conflicts with neighbors already in state k).
But this can fail globally for two reasons: (1) multiple state3 sites
that are mutually adjacent get converted independently, using each
other's PRE-conversion states, so they can end up matching each other
post-conversion even though a joint choice could have avoided it; (2) a
D3<2 ground state's own "background" pattern (e.g. whether the bulk uses
Case0's antiphase swap) was optimized assuming state3 remains available,
and might not be the right scaffold for the true D3=infinity optimum,
which could need a different background choice entirely.

Tested on two cases:
  1. A single, well-isolated straight string (no background-swap freedom,
     D1=D2=0, state3 sites well-separated) -- expect the rule to work.
  2. Case 2's bent 5-segment detour (known from case2_sa_check.py to have
     a real background antiphase-swap in its own D3=1 ground state) --
     the interesting test of whether background-swap freedom breaks it.
"""
import numpy as np

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, total_energy, simulated_annealing,
    full_state_via_mincut, frustrated_triangles,
)
from case2_sa_check import build_case2


def greedy_remove_state3(lat, states):
    states = states.copy()
    state3_sites = np.where(states == 2)[0]
    nbrs = lat.neighbor_table()
    for site in state3_sites:
        cnt = [0, 0, 0]
        for j, b in nbrs[site]:
            if b["J"] != 0.0 and states[j] != 2:
                cnt[states[j]] += 1
        states[site] = 0 if cnt[0] <= cnt[1] else 1
    return states


def sa_best(lat, D, rng, n_restarts=10, n_sweeps=3000, seed_states=None):
    best = None
    if seed_states is not None:
        s, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps, T_start=0.6,
                                     T_end=1e-6, states=seed_states, record_energy=True)
        best = (en[-1], s)
    for _ in range(n_restarts):
        s, en = simulated_annealing(lat, D, rng, n_sweeps=n_sweeps, T_start=5.0,
                                     T_end=1e-6, record_energy=True)
        if best is None or en[-1] < best[0]:
            best = (en[-1], s)
    return best


def test_simple_string():
    print("=== 1. simple isolated straight string ===")
    NX, NY = 20, 20
    lat = RhombileLattice(NX, NY)
    p_left, p_right = np.array([4.0, 10.0]), np.array([16.0, 10.0])
    apply_string_defect(lat, p_left, p_right, wrap=False)

    D3 = 1.9
    D = (0.0, 0.0, D3)
    states, violated = full_state_via_mincut(lat, D3)
    e_before = total_energy(lat, states, D)
    n3_before = int((states == 2).sum())
    print(f"D3={D3}: exact ground state E={e_before}, n_state3={n3_before}")

    states_fixed = greedy_remove_state3(lat, states)
    e_after = total_energy(lat, states_fixed, D)
    print(f"after greedy conversion: E={e_after} (state3 remaining: {(states_fixed == 2).sum()})")

    rng = np.random.default_rng(0)
    _, states_inf = sa_best(lat, (0.0, 0.0, 1000.0), rng, n_restarts=8)
    e_inf = total_energy(lat, states_inf, (0.0, 0.0, 1.0))
    print(f"SA-found true D3=inf ground state energy: {e_inf}")
    print(f"MATCH: {abs(e_after - e_inf) < 1e-9}\n")


def test_bent_case2():
    print("=== 2. Case 2 bent detour (has background antiphase-swap freedom) ===")
    lat, touched, corners, (i_left, i_right) = build_case2()
    ft = frustrated_triangles(lat)
    print(f"{len(ft)} real defects")

    D3 = 1.0
    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(0)
    e_before, states = sa_best(lat, D, rng, n_restarts=10)
    n3_before = int((states == 2).sum())
    print(f"D3={D3}: SA-best ground state E={e_before}, n_state3={n3_before}")

    states_fixed = greedy_remove_state3(lat, states)
    e_after = total_energy(lat, states_fixed, D)
    print(f"after greedy conversion: E={e_after} (state3 remaining: {(states_fixed == 2).sum()})")

    rng2 = np.random.default_rng(1)
    e_inf, _ = sa_best(lat, (0.0, 0.0, 1000.0), rng2, n_restarts=10)
    print(f"SA-found true D3=inf ground state energy: {e_inf}")
    print(f"MATCH: {abs(e_after - e_inf) < 1e-9}\n")


if __name__ == "__main__":
    test_simple_string()
    test_bent_case2()

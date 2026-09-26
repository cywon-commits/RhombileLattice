"""Generalizes hub_rim_order_transition.py to take D3 as a parameter, to
test the prediction that raising D3 (making the state-3 "safety valve"
more expensive) pushes the empirical hub/rim order-disorder T_c up
toward the pure-2-state dice-lattice theory value (1.203), since D3=1's
own T_c (~0.79 from the fine rescan) sits well below it specifically
because of substantial thermal state-3 population that the theory
doesn't account for. D3=3 should show much less state-3 population and
a T_c closer to 1.203.
"""
import csv
import os
import sys

import numpy as np

from rhombile_lattice import RhombileLattice, heat_bath_sweep
from hub_rim_order_transition import NX, NY, pristine_ground_state, order_parameter


def sample_at_T(lat, D, is_r1, T, seed, init_states, n_equil=2000, n_samples=2400, thin=4, block=60):
    rng = np.random.default_rng(seed)
    states = init_states.copy()
    for _ in range(n_equil):
        heat_bath_sweep(lat, states, D, T, rng)
    ms, n2_fracs = [], []
    for _ in range(n_samples):
        for _ in range(thin):
            heat_bath_sweep(lat, states, D, T, rng)
        m, n2f = order_parameter(states, is_r1)
        ms.append(abs(m))
        n2_fracs.append(n2f)
    ms = np.array(ms)
    n_blocks = len(ms) // block
    block_means = ms[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    return block_means, np.mean(n2_fracs)


def scan(D3, T_grid, seed=131):
    D = (0.0, 0.0, D3)
    csv_path = f"hub_rim_order_transition_D3_{D3:.1f}_results.csv"
    lat = RhombileLattice(NX, NY)
    pos, sub_of = lat.site_positions()
    is_r1 = sub_of == "r1"
    init_states = pristine_ground_state(lat)

    write_header = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["T", "abs_m_mean", "abs_m_se", "n2_frac_hub"])
        for T in T_grid:
            bm, n2_frac = sample_at_T(lat, D, is_r1, T, seed, init_states)
            m_mean, m_se = bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))
            print(f"D3={D3} T={T:<6} <|m|>={m_mean:.4f}+-{m_se:.4f}  hub_state3_frac={n2_frac:.4f}", flush=True)
            writer.writerow([T, m_mean, m_se, n2_frac])
            f.flush()


if __name__ == "__main__":
    D3 = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
    T_grid = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0]
    scan(D3, T_grid)

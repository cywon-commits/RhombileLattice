"""Fine re-scan of hub_rim_order_transition.py in T=[0.5, 1.0] (step
0.05) to pin down O's actual order-disorder crossover more precisely
than the original coarse grid (0.6, 0.8 straddled a drop from 0.99 to
0.48). More samples than the coarse scan since this is the region that
actually matters and is presumably the noisiest (near-critical
fluctuations).
"""
import csv
import os

import numpy as np

from rhombile_lattice import RhombileLattice, heat_bath_sweep
from hub_rim_order_transition import (
    NX, NY, D, pristine_ground_state, order_parameter,
)

CSV_PATH = "hub_rim_order_transition_fine_results.csv"


def sample_at_T(lat, is_r1, T, seed, init_states, n_equil=2000, n_samples=2400, thin=4, block=60):
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


def scan(T_grid, seed=73):
    lat = RhombileLattice(NX, NY)
    pos, sub_of = lat.site_positions()
    is_r1 = sub_of == "r1"
    init_states = pristine_ground_state(lat)

    write_header = not os.path.exists(CSV_PATH) or os.path.getsize(CSV_PATH) == 0
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["T", "abs_m_mean", "abs_m_se", "n2_frac_hub"])
        for T in T_grid:
            bm, n2_frac = sample_at_T(lat, is_r1, T, seed, init_states)
            m_mean, m_se = bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))
            print(f"T={T:<6} <|m|>={m_mean:.4f}+-{m_se:.4f}  hub_state3_frac={n2_frac:.4f}", flush=True)
            writer.writerow([T, m_mean, m_se, n2_frac])
            f.flush()


if __name__ == "__main__":
    scan([round(0.5 + 0.05 * k, 3) for k in range(11)])  # 0.5, 0.55, ..., 1.0

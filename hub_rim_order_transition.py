"""Locates the order-disorder transition temperature for the PRISTINE
(defect-free) lattice O, as a calibration baseline before redoing the
finite-T A-vs-B (short vs long defect-pair separation) comparison from
finite_T_scan.py more carefully.

At D1=D2=0 (used throughout this project's D3-sweeps), the pristine
lattice's hub (r1) sublattice spontaneously picks ONE of state1/state2 to
be the majority at T=0 (a Z2-like symmetry breaking: r1 all-state1/
rim all-state2, or the reverse -- both cost zero either way, an
accidental degeneracy from D1=D2, already central to Part IV/VI of this
project). Order parameter here: m = frac(r1 in state1) - frac(r1 in
state2) (near +-1 at T=0, ignoring the O(1) negligible chance of any r1
site using state3 in the bulk away from any string). Track <|m|> vs T,
seeded from the T=0 ground state at every T (required for equilibration,
per finite_T_scan.py's own hard-won lesson), to locate where it drops
from near its T=0 value toward 0 -- the crossover/transition temperature.

Uses the same lattice size (NX=30, NY=8) and D3=1 as finite_T_scan.py's
A/B for direct comparability, plus reference-free (no subtraction needed
here, this is the reference itself) batch-means error bars.
"""
import csv
import os

import numpy as np

from rhombile_lattice import RhombileLattice, heat_bath_sweep

NX, NY = 30, 8
D3 = 1.0
D = (0.0, 0.0, D3)
CSV_PATH = "hub_rim_order_transition_results.csv"


def pristine_ground_state(lat):
    pos, sub_of = lat.site_positions()
    return np.where(sub_of == "r1", 0, 1).astype(int)


def order_parameter(states, is_r1):
    hub = states[is_r1]
    n0 = (hub == 0).sum()
    n1 = (hub == 1).sum()
    n2 = (hub == 2).sum()
    denom = n0 + n1
    m = (n0 - n1) / denom if denom > 0 else 0.0
    return m, n2 / len(hub)


def sample_at_T(lat, is_r1, T, seed, init_states, n_equil=1200, n_samples=1200, thin=3, block=40):
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


def scan(T_grid, seed=71):
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
    scan([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0])

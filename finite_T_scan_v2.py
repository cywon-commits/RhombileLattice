"""More careful re-do of finite_T_scan.py (flagged earlier as needing
more rigor): same A (short) vs B (long, 4x separation) vs C (pristine)
reference-subtraction setup, but with ~2.7x more equilibration/sampling
sweeps and a T grid concentrated where it now matters -- we know from
hub_rim_order_transition.py that O's own hub/rim order collapses
somewhere around T~0.6-0.9, which is presumably why the original run's
T=0.8 point was unusably noisy (dB/dA = 3.81+-2.87). Denser grid through
that region, still going up to T=5 for the high-T confinement-ratio
tail.
"""
import csv
import os

import numpy as np

from rhombile_lattice import (
    RhombileLattice, heat_bath_sweep, total_energy, full_state_via_mincut,
)
from distance_scaling_demo import build, DIST_A, DIST_B, NX, NY

D3 = 1.0
D = (0.0, 0.0, D3)
CSV_PATH = "finite_T_scan_v2_results.csv"


def ground_state(lat):
    states, _ = full_state_via_mincut(lat, D3)
    return states


def pristine_ground_state(lat):
    pos, sub_of = lat.site_positions()
    return np.where(sub_of == "r1", 0, 1).astype(int)


def sample_at_T(lat, T, seed, init_states, n_equil=2000, n_samples=2400, thin=4, block=60):
    rng = np.random.default_rng(seed)
    states = init_states.copy()
    for _ in range(n_equil):
        heat_bath_sweep(lat, states, D, T, rng)
    energies = []
    for _ in range(n_samples):
        for _ in range(thin):
            heat_bath_sweep(lat, states, D, T, rng)
        energies.append(total_energy(lat, states, D))
    energies = np.array(energies)
    n_blocks = len(energies) // block
    block_means = energies[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    return block_means, energies


def scan(T_grid, seed=97):
    lat_A, lat_B = build(DIST_A), build(DIST_B)
    lat_C = RhombileLattice(NX, NY)
    init_A, init_B, init_C = ground_state(lat_A), ground_state(lat_B), pristine_ground_state(lat_C)

    write_header = not os.path.exists(CSV_PATH) or os.path.getsize(CSV_PATH) == 0
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["T", "dA_mean", "dA_se", "dB_mean", "dB_se", "ratio", "ratio_se"])
        for T in T_grid:
            bm_A, _ = sample_at_T(lat_A, T, seed, init_A)
            bm_B, _ = sample_at_T(lat_B, T, seed, init_B)
            bm_C, _ = sample_at_T(lat_C, T, seed, init_C)
            dA = bm_A - bm_C
            dB = bm_B - bm_C
            dA_mean, dA_se = dA.mean(), dA.std(ddof=1) / np.sqrt(len(dA))
            dB_mean, dB_se = dB.mean(), dB.std(ddof=1) / np.sqrt(len(dB))
            ratio = dB_mean / dA_mean if abs(dA_mean) > 1e-9 else float("nan")
            # error propagation for a ratio of two correlated-ish means (treat as independent,
            # a conservative/slightly-overestimated error bar given dA,dB share the same C draws)
            ratio_se = (abs(ratio) * np.sqrt((dA_se / dA_mean) ** 2 + (dB_se / dB_mean) ** 2)
                        if abs(dA_mean) > 1e-9 and abs(dB_mean) > 1e-9 else float("nan"))
            print(f"T={T:<6} dA={dA_mean:.3f}+-{dA_se:.3f}  dB={dB_mean:.3f}+-{dB_se:.3f}  "
                  f"ratio={ratio:.3f}+-{ratio_se:.3f}", flush=True)
            writer.writerow([T, dA_mean, dA_se, dB_mean, dB_se, ratio, ratio_se])
            f.flush()


if __name__ == "__main__":
    T_grid = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0]
    scan(T_grid)

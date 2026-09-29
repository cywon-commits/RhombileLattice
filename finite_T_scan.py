"""Finite-T companion to distance_scaling_demo.py: fix D3=1, raise T,
and track how the SHORT (Case A) vs LONG (Case B) defect-pair's extra
energy (relative to a same-size pristine reference C, no string at
all) behaves. Goal: does the T=0 confinement (E_B/E_A = 4 exactly,
distance_scaling_demo.py) persist, shrink, or something else as T
rises -- i.e. is there a finite-T entropic effect of the kind Moore et
al. describe (see research_notes_next_directions.md, section 14.1)?

Two things needed to get a trustworthy signal, both because the bulk
lattice's own thermal energy is extensive (~N) and can easily swamp
the O(L) defect-pair signal if not controlled:

1. Reference subtraction: sample a same-size (NX=30, NY=8, N=720)
   PRISTINE lattice C (no string) at the same T, and report
   dX(T) = <E_X>_T - <E_C>_T for X in {A, B}, not raw <E_X>_T.
2. Common random numbers: A, B, C are sampled with independently-
   seeded but IDENTICALLY-SEEDED rng streams (same seed => same
   sequence of site-visitation permutations and Gibbs draws each
   sweep, since heat_bath_sweep consumes a fixed-shape draw regardless
   of the local weights) -- this correlates their bulk thermal
   fluctuations and should reduce the variance of dX relative to fully
   independent sampling.

Also important, found the hard way in an earlier quick prototype:
equilibration MUST start from the known T=0 ground state (from
full_state_via_mincut for A/B, the trivial hub/rim checkerboard for
C), not a fully random state -- a random start at low T did not
equilibrate in a few hundred sweeps and gave energies far above the
true thermal average (a burn-in artifact, not a physical result).

Error bars use batch means (split the post-equilibration samples into
blocks and treat block means as the i.i.d. quantity), which is at
least aware of the autocorrelation problem the user flagged, if not a
full integrated-autocorrelation-time estimate.
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
CSV_PATH = "finite_T_scan_results.csv"


def ground_state(lat):
    states, _ = full_state_via_mincut(lat, D3)
    return states


def pristine_ground_state(lat):
    pos, sub_of = lat.site_positions()
    return np.where(sub_of == "r1", 0, 1).astype(int)


def sample_at_T(lat, T, seed, init_states, n_equil=800, n_samples=900, thin=3, block=30):
    """Returns (block_means, all_samples) after equilibrating from
    init_states. block_means are used for the batch-means error bar."""
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


def scan(T_grid, seed=42):
    lat_A, lat_B = build(DIST_A), build(DIST_B)
    lat_C = RhombileLattice(NX, NY)
    init_A, init_B, init_C = ground_state(lat_A), ground_state(lat_B), pristine_ground_state(lat_C)

    write_header = not os.path.exists(CSV_PATH) or os.path.getsize(CSV_PATH) == 0
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["T", "dA_mean", "dA_se", "dB_mean", "dB_se", "ratio"])
        for T in T_grid:
            bm_A, _ = sample_at_T(lat_A, T, seed, init_A)
            bm_B, _ = sample_at_T(lat_B, T, seed, init_B)
            bm_C, _ = sample_at_T(lat_C, T, seed, init_C)
            # batch means are paired sample-by-sample (same block index -> same
            # rng draws consumed up to that point), so subtract block-by-block
            dA = bm_A - bm_C
            dB = bm_B - bm_C
            dA_mean, dA_se = dA.mean(), dA.std(ddof=1) / np.sqrt(len(dA))
            dB_mean, dB_se = dB.mean(), dB.std(ddof=1) / np.sqrt(len(dB))
            ratio = dB_mean / dA_mean if abs(dA_mean) > 1e-9 else float("nan")
            print(f"T={T:<6} dA={dA_mean:.3f}+-{dA_se:.3f}  dB={dB_mean:.3f}+-{dB_se:.3f}  "
                  f"ratio={ratio:.3f}", flush=True)
            writer.writerow([T, dA_mean, dA_se, dB_mean, dB_se, ratio])
            f.flush()


if __name__ == "__main__":
    scan([0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0])

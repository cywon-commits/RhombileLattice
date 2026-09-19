"""Does the T=0 state-3-vs-violated-link transition (exactly D3=2 for a
plain straight string, section 6/12 of the report) stay at D3=2 once
T>0, or shift due to an entropy difference between the two ways of
paying for a defect pair?

Quick combinatorial argument (shared with the user before running
this): for a plain straight string with L conflict edges, using k
non-overlapping state-3 "domino" placements always leaves exactly
L-2k violated links, so E(k) = D3*k + (L-2k) = (D3-2)*k + L is EXACTLY
LINEAR in k. The number of ways to place k such dominoes, Omega(k), is
a combinatorial function of k alone (independent of D3, T) that is
zero at the two extremes (k=L/2 and k=0, each realized only 1 way) and
larger in between. Minimizing F(k) = E(k) - T*log(Omega(k)) gives
(D3-2) = T * d(log Omega)/dk -- so the point where the optimal k sits
at Omega's own maximum (d(log Omega)/dk = 0) is EXACTLY D3=2 for any
T, in this idealized picture. Prediction: the transition should
BROADEN with T but its center might not shift. This rests on an
idealized (hub-fixed, no other entropy channels) picture though, so
it's a hypothesis to test, not a confident prediction -- real thermal
fluctuations elsewhere in the lattice could break the symmetry.

Method: fix T=0.5, sweep D3 finely (extra resolution near 2) for
Case B (L=40, the long string from distance_scaling_demo.py). At each
D3, seed from the EXACT T=0 ground state for that D3 (full_state_via_
mincut) to avoid the burn-in problem found in the first finite-T
prototype, equilibrate, then sample <n_state3> and <n_violated> with
block-mean error bars.

IMPORTANT FIX (caught by the user before trusting any numbers):
n_state3 must be counted only among the string's own local
conflict-graph sites (the endpoints of its L active diagonals), NOT
over the whole lattice. state-3 isn't special to the string -- any
bulk site can thermally flip into it too, and with ~700 bulk sites
versus the string's ~40 relevant ones, that background swamped the
signal in a first attempt (D3=0.5 gave n_state3=142, when the string
alone can hold at most L/2=20 -- a 7x contamination). n_violated
doesn't have this problem since r2-r3 bonds are off everywhere in the
pristine lattice except where this string turned them on.

Two crossover definitions are reported: (a) where <n_state3> crosses
half its own D3->0 plateau (L/4), and (b) where 2*<n_state3> (edges
resolved) crosses <n_violated> (edges left unresolved) directly --
(b) doesn't need to assume what the plateau value is, so it's the
more robust one; both reduce exactly to D3=2 in the T->0 limit.
"""
import csv
import os

import numpy as np

from rhombile_lattice import heat_bath_sweep, full_state_via_mincut, conflict_graph_bipartition
from distance_scaling_demo import build, DIST_B

T_FIXED = 0.5
CSV_PATH = "D3_transition_finite_T_results.csv"
D3_GRID = [0.5, 0.75, 1.0, 1.25, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
           2.1, 2.2, 2.3, 2.4, 2.5, 2.75, 3.0, 3.5]


def _string_sites(lat):
    """Sites that are endpoints of the string's own active r2-r3
    diagonals -- the only places n_state3 should be counted, to avoid
    bulk thermal contamination (see module docstring)."""
    edges, _, _, _ = conflict_graph_bipartition(lat)
    sites = set()
    for i, j in edges:
        sites.add(i); sites.add(j)
    return np.array(sorted(sites))


def n_state3_and_violated(lat, states, string_sites):
    n_state3 = int((states[string_sites] == 2).sum())
    n_violated = sum(1 for b in lat.bonds
                      if b["type"] == "r2r3" and b["J"] != 0.0 and states[b["i"]] == states[b["j"]])
    return n_state3, n_violated


def sample_at_D3(lat, D3, T, seed, string_sites, n_equil=800, n_samples=900, thin=3, block=30):
    D = (0.0, 0.0, D3)
    init_states, _ = full_state_via_mincut(lat, D3)
    rng = np.random.default_rng(seed)
    states = init_states.copy()
    for _ in range(n_equil):
        heat_bath_sweep(lat, states, D, T, rng)
    n3_samples, nv_samples = [], []
    for _ in range(n_samples):
        for _ in range(thin):
            heat_bath_sweep(lat, states, D, T, rng)
        n3, nv = n_state3_and_violated(lat, states, string_sites)
        n3_samples.append(n3)
        nv_samples.append(nv)
    n3_samples, nv_samples = np.array(n3_samples), np.array(nv_samples)
    n_blocks = len(n3_samples) // block
    n3_blocks = n3_samples[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    nv_blocks = nv_samples[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    return n3_blocks, nv_blocks


def scan(D3_grid=D3_GRID, T=T_FIXED, seed=7):
    lat = build(DIST_B)
    string_sites = _string_sites(lat)
    print(f"string_sites: {len(string_sites)} local sites (L={len(string_sites) - 1} edges expected)")
    write_header = not os.path.exists(CSV_PATH) or os.path.getsize(CSV_PATH) == 0
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["D3", "n3_mean", "n3_se", "nv_mean", "nv_se"])
        for D3 in D3_grid:
            n3_blocks, nv_blocks = sample_at_D3(lat, D3, T, seed, string_sites)
            n3_mean, n3_se = n3_blocks.mean(), n3_blocks.std(ddof=1) / np.sqrt(len(n3_blocks))
            nv_mean, nv_se = nv_blocks.mean(), nv_blocks.std(ddof=1) / np.sqrt(len(nv_blocks))
            print(f"D3={D3:<6} n3={n3_mean:.2f}+-{n3_se:.2f}  nv={nv_mean:.2f}+-{nv_se:.2f}", flush=True)
            writer.writerow([D3, n3_mean, n3_se, nv_mean, nv_se])
            f.flush()


if __name__ == "__main__":
    scan()

"""Finite-size scan for entropic_force_T0_fixedwin.py: Moore et al.
themselves (Table I of the paper) checked their own A(L) coefficient
across 4 lattice sizes (32,64,128,256) and extrapolated A(L)=A+CL^-alpha
as L->infinity, explicitly because their smaller sizes showed "considerable
finite-size effects" (their words) -- exactly the kind of check our own
single-box-size (NX=40) run has not had yet. Our own max separation there
(L=20) was already half the box width, so the apparent saturation between
L=16 and L=20 could easily be a boundary artifact rather than real physics.

This reruns the same fixed-window methodology on a much bigger torus
(NX=80, NY=30) at the SAME L values (12,16,20) to check whether they
still agree with the NX=40 result (a direct box-size-robustness check),
plus extends to L=28,36 (safely small relative to the bigger box) to see
whether growth resumes, keeps saturating, or does something else once
boundary effects are removed.
"""
import pickle

import numpy as np

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, total_energy, state_energies,
    full_state_via_mincut,
)

D3 = 0.0
NX, NY = 80, 30
Y_CROSS = 15
N_SWEEPS_EQUIL = 1200
N_SWEEPS_MEASURE = 5000
SEPARATIONS = [12, 16, 20, 28, 36]


def zero_temperature_sweep(lattice, states, D, rng):
    """One sweep of Moore et al.'s own T=0 dynamics: visit every site in
    random order, propose a uniformly random NEW state, accept iff it
    does not increase that site's own local energy (ties always
    accepted, so this is a random walk WITHIN the ground-state manifold
    once E=0 is reached, exactly matching their zero-temperature
    protocol -- not a designed-in-advance heat-bath schedule)."""
    n_accepted = 0
    for site in rng.permutation(lattice.n_sites):
        e = state_energies(lattice, states, site, D)
        cur = states[site]
        candidates = [k for k in range(3) if k != cur]
        new = candidates[rng.integers(len(candidates))]
        if e[new] <= e[cur] + 1e-12:
            states[site] = new
            n_accepted += 1
    return n_accepted


def pristine_ground_state(lattice):
    pos, sub_of = lattice.site_positions()
    return np.where(sub_of == "r1", 0, 1).astype(int)


def run_one_separation(L, seed):
    """Reference-subtracted, exactly as elsewhere in this project
    (finite_T_scan*.py): the raw n3 count is dominated by a huge,
    string-INDEPENDENT bulk gas of state-3 sites (D3=0 makes the 3rd
    state free everywhere, not just near the string), so the string's
    OWN contribution has to be isolated by subtracting a same-size,
    same-seed pristine run -- common random numbers, run sweep-by-sweep
    in parallel on both lattices."""
    lat = RhombileLattice(NX, NY)
    p_left = np.array([NX / 2 - L / 2, Y_CROSS])
    p_right = np.array([NX / 2 + L / 2, Y_CROSS])
    apply_string_defect(lat, p_left, p_right, wrap=False)
    states, violated = full_state_via_mincut(lat, D3=0.0)
    e0 = total_energy(lat, states, D=(0.0, 0.0, D3))
    assert abs(e0) < 1e-9, f"expected E=0 at D3=0, got {e0} (L={L})"

    lat_ref = RhombileLattice(NX, NY)
    states_ref = pristine_ground_state(lat_ref)
    assert abs(total_energy(lat_ref, states_ref, (0.0, 0.0, D3))) < 1e-9

    # FIXED window (independent of L, unlike entropic_force_T0.py's own
    # L-scaled window): the same absolute box for every separation, large
    # enough to contain the longest string (L=20) plus margin, centered on
    # the lattice midpoint -- isolates the L-dependence of the SIGNAL from
    # the L-dependence of the window's own bulk-background contribution.
    pos, sub_of = lat.site_positions()
    margin = 4.0
    half_box = max(SEPARATIONS) / 2 + margin
    mid_x = NX / 2
    window = (np.abs(pos[:, 1] - Y_CROSS) <= margin) & \
             (pos[:, 0] >= mid_x - half_box) & (pos[:, 0] <= mid_x + half_box)

    rng = np.random.default_rng(seed)
    rng_ref = np.random.default_rng(seed)
    for _ in range(N_SWEEPS_EQUIL):
        zero_temperature_sweep(lat, states, (0.0, 0.0, D3), rng)
        zero_temperature_sweep(lat_ref, states_ref, (0.0, 0.0, D3), rng_ref)

    d_trace, n3_trace, n3ref_trace = [], [], []
    for _ in range(N_SWEEPS_MEASURE):
        zero_temperature_sweep(lat, states, (0.0, 0.0, D3), rng)
        zero_temperature_sweep(lat_ref, states_ref, (0.0, 0.0, D3), rng_ref)
        n3 = int((states[window] == 2).sum())
        n3ref = int((states_ref[window] == 2).sum())
        n3_trace.append(n3)
        n3ref_trace.append(n3ref)
        d_trace.append(n3 - n3ref)

    d_trace = np.array(d_trace, dtype=float)
    block = 100
    n_blocks = len(d_trace) // block
    block_means = d_trace[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    d_sem = float(block_means.std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else float("nan")
    return {
        "L": L,
        "d_mean": float(d_trace.mean()), "d_std": float(d_trace.std()), "d_sem": d_sem,
        "d_trace": d_trace, "block_means": block_means,
        "n3_mean": float(np.mean(n3_trace)), "n3ref_mean": float(np.mean(n3ref_trace)),
    }


def main():
    results = {}
    for L in SEPARATIONS:
        seed_vals = []
        for seed_base in (1000, 7777):
            r = run_one_separation(L, seed=seed_base + L)
            print(f"L={L:3d} seed={seed_base+L}  d_mean={r['d_mean']:+.3f} +- {r['d_sem']:.3f} (SEM)  "
                  f"(n3={r['n3_mean']:.1f}, n3ref={r['n3ref_mean']:.1f})", flush=True)
            seed_vals.append(r)
        results[L] = seed_vals
        means = [x["d_mean"] for x in seed_vals]
        print(f"  -> L={L}: seed-avg d_mean = {np.mean(means):.3f} "
              f"(spread {np.std(means, ddof=1):.3f})", flush=True)
    with open("entropic_force_T0_sizescan.pkl", "wb") as f:
        pickle.dump(results, f)
    print("saved entropic_force_T0_sizescan.pkl")


if __name__ == "__main__":
    main()

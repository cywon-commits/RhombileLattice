"""Faithful(-as-possible) replication of Moore, Nordahl, Minar & Shalizi
(PRE 60, 5344, 1999), Sec. IV-V, on THIS lattice, at D3=0 (our own
zero-bare-energy sector -- the honest place to look for their T=0
ground-state-ENTROPY-driven Coulomb force, since it has nothing to do
with raising temperature past any order-disorder transition).

Structural caveat, made explicit rather than glossed over: in their model
a "defect" is a property of the CURRENT SPIN CONFIGURATION (a frustrated
bond that moves as spins flip -- genuine vortex dynamics). In ours, a
frustrated triangle is a property of the FIXED BOND PATTERN (which edges
are active), so it structurally cannot move under single-site Potts
flips alone; only the escape/state-3 assignment satisfying it can. So
instead of tracking a moving defect's drift, we track the natural analog
available here: the escape-site "cloud" that resolves the 2 fixed
defects' conflict graph, under the SAME T=0 dynamics they use (single
random-site update, accepted iff it does not increase the energy),
and ask whether its size/fluctuations depend on the defect separation r
-- a direct proxy for whether the local ground-state ENTROPY depends on r,
which is the actual quantity their Gaussian-interface argument predicts
(via log r), not an energy effect (energy is exactly 0 for every r here).
"""
import pickle

import numpy as np

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, total_energy, state_energies,
    full_state_via_mincut,
)

D3 = 0.0
NX, NY = 40, 20
Y_CROSS = 10
N_SWEEPS_EQUIL = 2000
N_SWEEPS_MEASURE = 8000
SEPARATIONS = [4, 8, 12, 16, 20]


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

    # windowed count: a fixed perpendicular margin around the string's own
    # row, restricted horizontally to the string's own span + a fixed
    # margin -- this excludes the huge, string-INDEPENDENT bulk gas of
    # state-3 sites elsewhere in the lattice, which otherwise swamps the
    # signal with irrelevant variance (checked empirically: raw whole-
    # lattice n3 fluctuates by +-30 with essentially no L-dependence
    # visible against that noise floor).
    pos, sub_of = lat.site_positions()
    margin = 4.0
    window = (np.abs(pos[:, 1] - Y_CROSS) <= margin) & \
             (pos[:, 0] >= p_left[0] - margin) & (pos[:, 0] <= p_right[0] + margin)

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
        r = run_one_separation(L, seed=1000 + L)
        print(f"L={L:3d}  d_mean={r['d_mean']:+.3f} +- {r['d_sem']:.3f} (SEM)  "
              f"d_std={r['d_std']:.3f}  (n3={r['n3_mean']:.1f}, n3ref={r['n3ref_mean']:.1f})",
              flush=True)
        results[L] = r
    with open("entropic_force_T0.pkl", "wb") as f:
        pickle.dump(results, f)
    print("saved entropic_force_T0.pkl")


if __name__ == "__main__":
    main()

"""Finite-size scaling for O's own hub/rim order-disorder transition
(pristine lattice, D3=1), still flagged open in the artifact's checklist.
Everything reported so far (T_c~0.79 at D3=1) came from a single lattice
size (NX=30, NY=8) -- on any finite lattice this is technically a smooth
crossover, not a true singularity, so before trusting it as a genuine
thermodynamic-limit phase transition we need to see it SHARPEN with size
and locate where it actually converges.

Two standard diagnostics, both computed from the same raw per-sample
order-parameter trace (not just its mean, unlike hub_rim_order_transition*.py):
  - <|m|>(T) at several sizes: a genuine transition should sharpen (steeper
    drop) as the system grows, not just look like a fixed-width crossover.
  - Binder cumulant U_L(T) = 1 - <m^4>/(3<m^2>^2): curves for different L
    should cross at approximately the same T, which converges to T_c as
    L->infinity (the standard, size-independent way to locate a critical
    point without needing to know its universality class in advance).

Uses SQUARE boxes (NX=NY=L) rather than the original elongated 30x8, since
clean finite-size scaling theory assumes an isotropic linear size; this
is a different (cleaner) comparison than a literal reproduction of the
earlier 0.79 number, which used a different aspect ratio.
"""
import pickle

import numpy as np

from rhombile_lattice import RhombileLattice, heat_bath_sweep

D3 = 1.0
D = (0.0, 0.0, D3)
SIZES = [6, 8, 10, 14]
T_GRID = [round(0.5 + 0.05 * k, 3) for k in range(11)]  # 0.5 .. 1.0


def pristine_ground_state(lat):
    pos, sub_of = lat.site_positions()
    return np.where(sub_of == "r1", 0, 1).astype(int)


def signed_order_parameter(states, is_r1):
    hub = states[is_r1]
    n0 = (hub == 0).sum()
    n1 = (hub == 1).sum()
    denom = n0 + n1
    return (n0 - n1) / denom if denom > 0 else 0.0


def sample_at_T(lat, is_r1, T, seed, init_states, n_equil=800, n_samples=1000, thin=3, block=50):
    rng = np.random.default_rng(seed)
    states = init_states.copy()
    for _ in range(n_equil):
        heat_bath_sweep(lat, states, D, T, rng)
    ms = np.empty(n_samples)
    for k in range(n_samples):
        for _ in range(thin):
            heat_bath_sweep(lat, states, D, T, rng)
        ms[k] = signed_order_parameter(states, is_r1)
    n_blocks = len(ms) // block
    abs_block_means = np.abs(ms[:n_blocks * block]).reshape(n_blocks, block).mean(axis=1)
    m2_block_means = (ms[:n_blocks * block] ** 2).reshape(n_blocks, block).mean(axis=1)
    m4_block_means = (ms[:n_blocks * block] ** 4).reshape(n_blocks, block).mean(axis=1)
    return {
        "abs_m": float(abs_block_means.mean()),
        "abs_m_se": float(abs_block_means.std(ddof=1) / np.sqrt(n_blocks)),
        "m2": float(m2_block_means.mean()),
        "m4": float(m4_block_means.mean()),
    }


def binder(m2, m4):
    return 1.0 - m4 / (3.0 * m2 ** 2)


def scan():
    results = {}
    for L in SIZES:
        lat = RhombileLattice(L, L)
        pos, sub_of = lat.site_positions()
        is_r1 = sub_of == "r1"
        init_states = pristine_ground_state(lat)
        per_T = {}
        for T in T_GRID:
            r = sample_at_T(lat, is_r1, T, seed=1000 + L, init_states=init_states)
            U = binder(r["m2"], r["m4"])
            print(f"L={L:3d} T={T:<5} <|m|>={r['abs_m']:.4f}+-{r['abs_m_se']:.4f}  "
                  f"Binder U={U:.4f}", flush=True)
            per_T[T] = {**r, "U": U}
            results[L] = per_T
            with open("hub_rim_Tc_finite_size.pkl", "wb") as f:
                pickle.dump(results, f)
    print("saved hub_rim_Tc_finite_size.pkl")


if __name__ == "__main__":
    scan()

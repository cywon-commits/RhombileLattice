"""Calibration step suggested in conversation: rather than guessing D3
values for the finite-size T_c scan, first fix T=1.5 (safely above every
T_c seen so far, so we're solidly in the disordered phase where thermal
state-3 excitation is most visible) and sweep D3 on a single moderate
lattice, measuring the hub 3rd-state thermal fraction directly, to find
the D3 where it's actually suppressed -- rather than assuming D3=3 was
already "large enough" (it visibly wasn't: T_c(D3=3, L->inf)~1.05-1.11,
still well below 1.2027, with observable curvature in T_c(D3) so far).

RESULT (L=10, T=1.5): hub_state3_frac decays roughly exponentially with
D3: 0.314(D3=1) -> 0.186(3) -> 0.119(4) -> 0.040(6) -> 0.011(8) ->
0.0032(10) -> 0.0001(15) -> 0.0000(20). D3=10-12 gives sub-0.5% thermal
state3 population -- a reasonable "sufficiently suppressed" calibration
point for the follow-up finite-size scan.
"""
import numpy as np

from rhombile_lattice import RhombileLattice, heat_bath_sweep
from hub_rim_order_transition import pristine_ground_state, order_parameter

T_FIXED = 1.5
L = 10
D3_GRID = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0]


def sample(lat, is_r1, D3, seed, init_states, n_equil=1000, n_samples=1500, thin=3, block=50):
    D = (0.0, 0.0, D3)
    rng = np.random.default_rng(seed)
    states = init_states.copy()
    for _ in range(n_equil):
        heat_bath_sweep(lat, states, D, T_FIXED, rng)
    n2_fracs = np.empty(n_samples)
    for k in range(n_samples):
        for _ in range(thin):
            heat_bath_sweep(lat, states, D, T_FIXED, rng)
        _, n2f = order_parameter(states, is_r1)
        n2_fracs[k] = n2f
    n_blocks = len(n2_fracs) // block
    block_means = n2_fracs[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    return float(block_means.mean()), float(block_means.std(ddof=1) / np.sqrt(n_blocks))


def main():
    lat = RhombileLattice(L, L)
    pos, sub_of = lat.site_positions()
    is_r1 = sub_of == "r1"
    init_states = pristine_ground_state(lat)

    for D3 in D3_GRID:
        mean, se = sample(lat, is_r1, D3, seed=int(D3 * 100) + 1, init_states=init_states)
        print(f"D3={D3:5.1f}  hub_state3_frac={mean:.4f} +- {se:.4f}", flush=True)


if __name__ == "__main__":
    main()

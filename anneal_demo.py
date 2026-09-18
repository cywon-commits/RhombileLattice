"""Sanity check for the anisotropic 3-state Potts SA loop on a clean
(defect-free) rhombile lattice: D = (D1, D2, D3) = (0, 1, 10), default J.
"""
import numpy as np

from rhombile_lattice import (
    RhombileLattice, SUB_INDEX, simulated_annealing, total_energy,
)

D = (0.0, 1.0, 10.0)
nx, ny = 8, 8
lat = RhombileLattice(nx, ny)
N = nx * ny  # sites per sublattice

sub_of = np.array([SUB_INDEX[s] for s in ["r1", "r2", "r3"]])  # 0,1,2
sub_labels = np.empty(lat.n_sites, dtype=int)
for m in range(ny):
    for n in range(nx):
        for sub, k in SUB_INDEX.items():
            sub_labels[lat._site_index(n, m, sub)] = k

# reference configurations, evaluated with the real total_energy()
uniform_state1 = np.zeros(lat.n_sites, dtype=int)
guess_a = np.where(sub_labels == 0, 0, 1)  # r1 -> state1(0), r2,r3 -> state2(1)
guess_b = np.where(sub_labels == 0, 1, 0)  # r1 -> state2(1), r2,r3 -> state1(0)
print(f"E[all state1]                 = {total_energy(lat, uniform_state1, D):.1f}")
print(f"E[r1=state1, r2=r3=state2]    = {total_energy(lat, guess_a, D):.1f}")
print(f"E[r1=state2, r2=r3=state1]    = {total_energy(lat, guess_b, D):.1f}")

rng = np.random.default_rng(0)
best_states, best_E = None, np.inf
for trial in range(5):
    states, energies = simulated_annealing(
        lat, D, rng, n_sweeps=300, T_start=4.0, T_end=0.005, record_energy=True
    )
    E = energies[-1]
    print(f"trial {trial}: final E = {E:.1f}")
    if E < best_E:
        best_E, best_states = E, states

print(f"\nbest SA energy = {best_E:.1f}")
for label, k in SUB_INDEX.items():
    counts = np.bincount(best_states[sub_labels == k], minlength=3)
    print(f"  {label}: state counts (1,2,3) = {counts}")

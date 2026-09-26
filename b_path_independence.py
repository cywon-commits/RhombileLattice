"""(b) Is the defect energy path-independent? Exact D3=infinity energies
(exact_Einf.py) for fixed defect positions on many different dimer
coverings.

Coverings are changed by closed-loop worm updates that never touch the
real defects: split a random rhombus, let one of the two new monomers do
a symmetric random walk (never onto a real defect) until it steps onto
its own partner, then re-pair them. Each update XORs a random closed
loop into the covering -- contractible or winding -- and leaves every
real defect exactly where it was. If the energy depended only on defect
positions (and flux sector), E_inf would be the same for all coverings
in a sector.

Usage:
  python3 b_path_independence.py pair <hops> <n_loops> <seed>
  python3 b_path_independence.py many <p> <n_loops> <seed> [knn]
"""
import sys

import numpy as np

from rhombile_lattice import frustrated_triangles
from exact_Einf import exact_two_state_energy, exact_two_state_energy_milp
from ensemble_comparison_A import PottsTrack, all_pairs_hop

NX = NY = 16


def loop_update(st, rng, frozen, max_steps=200000):
    """XOR one random closed loop into the covering; frozen = real defects."""
    while True:
        x = int(rng.integers(st.n))
        if st.partner[x] is not None and x not in frozen and st.partner[x] not in frozen:
            break
    a, b = st.split(x)
    walker, other = (a, b) if rng.random() < 0.5 else (b, a)
    for _ in range(max_steps):
        u = st.nbrs[walker][int(rng.integers(3))][0]
        if u == other:
            st.bond_between(walker, other)["J"] = 0.0
            st.partner[walker], st.partner[other] = other, walker
            st.monomers.discard(walker)
            st.monomers.discard(other)
            return True
        if u in frozen or u in st.monomers:
            continue
        walker = st.move(walker, u)
    raise RuntimeError("loop did not close")


def sector_signature(st):
    """Z2 x Z2 winding of the covering relative to the pristine tiling:
    parity of inactive bonds across each box seam (a bond that wraps at a
    corner crosses both seams and counts for both)."""
    from rhombile_lattice import A_MAT_INV
    wx = wy = 0
    for b in st.lat.bonds:
        if b.get("wraps") and b["J"] == 0.0:
            fu, fv = A_MAT_INV @ (b["p2"] - b["p1"])
            wx ^= int(abs(fu) > st.lat.nx / 2)
            wy ^= int(abs(fv) > st.lat.ny / 2)
    return wx, wy


def run_pair(hops, n_loops, seed):
    pt = PottsTrack(NX, NY, lmax=hops, seed=seed)
    st = pt.st
    for _ in range(hops):                   # walk B straight out along the track
        t, w = pt.track[pt.L], pt.track[pt.L + 1]
        st.move(t, st.partner[w])
        pt.L += 1
    A, B = pt.anchor, pt.track[pt.L]
    assert st.monomers == {A, B}
    dist = all_pairs_hop(st)
    rng = np.random.default_rng(seed + 100)
    print(f"# pair: hops={hops}  triangle-graph distance d_tri(A,B)={dist[A, B]}")
    print("# loop  lower  E_inf  certified  sector")
    lo, e, ok = exact_two_state_energy(st.lat)
    print(f"0 {lo} {e} {ok} {sector_signature(st)}", flush=True)
    for k in range(1, n_loops + 1):
        loop_update(st, rng, frozen={A, B})
        lo, e, ok = exact_two_state_energy(st.lat)
        print(f"{k} {lo} {e} {ok} {sector_signature(st)}", flush=True)
    assert len(frustrated_triangles(st.lat)) == 2


def run_many(p, n_loops, seed, knn=None):
    from worm_density_construction import build
    rng = np.random.default_rng(seed)
    st, _ = build(p, rng, n_sweeps=600)
    frozen = set(st.monomers)
    m = len(frozen)
    print(f"# many: p={m / st.n:.4f}  monomers={m}")
    note = "exact, MILP" if knn is None else f"MILP over {knn} nearest partners: certified upper bound"
    print(f"# loop  lower(no homology)  E_inf({note})  sector")
    rng2 = np.random.default_rng(seed + 100)
    for k in range(n_loops + 1):
        if k:
            loop_update(st, rng2, frozen=frozen)
        lo, _, _ = exact_two_state_energy(st.lat, max_fixups=0)
        e, _ = exact_two_state_energy_milp(st.lat, knn=knn)
        print(f"{k} {lo} {e} {sector_signature(st)}", flush=True)
    assert st.monomers == frozen


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "pair":
        run_pair(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    else:
        run_many(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]),
                 int(sys.argv[5]) if len(sys.argv) > 5 else None)

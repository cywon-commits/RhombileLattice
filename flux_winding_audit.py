"""Part (a) of the flux/winding-sector check requested after noticing that
dimer/lozenge-tiling models on a TORUS have winding-number sectors that are
NOT connected by local flips (Thurston's theorem only guarantees full flip-
connectivity for simply-connected regions). Question: did any construction
used so far in this project (Case0/1/2, the D2/D3/finite-T sweeps, the
distance-scaling demo) ever actually cross the lattice's own periodic seam?

Two independent lines of evidence:
  1. Structural: enumerate_triangles() (used by build_dual -> build_rhombi
     -> build_triangle_hop_graph, i.e. everything built on top of
     route_string_between_points/apply_dual_string_defect: Case0, Case1,
     Case2, and every D2/D3/finite-T sweep that warm-starts from
     case2_sa_check.build_case2()) explicitly EXCLUDES any triangle that
     touches a wrap bond, so *no* path built through that dual graph can
     ever use one -- this is true by construction, checked here by
     directly counting how many candidate triangles got excluded.
  2. Empirical: for the specific bond sets these demos actually toggled
     (Case0's closed loop, Case1's direct edge, Case2's 5-segment detour),
     check directly whether any flipped bond has wraps=True.

distance_scaling_demo.py uses the older straight-line apply_string_defect
with wrap=False, which is even more directly incapable of touching a wrap
bond (no shift loop is ever run), so it isn't re-checked here.
"""
import numpy as np

from rhombile_lattice import (
    RhombileLattice, enumerate_triangles, route_string_between_points,
    apply_dual_string_defect,
)
from dual_string_demo import build_dual
from closed_loop_demo import hexagon_corners, NX, NY, apply_closed_loop
from case2_sa_check import build_case2


def uses_periodic_wrap(flipped_bonds):
    return any(b["wraps"] for b in flipped_bonds)


def count_all_vs_enumerated_triangles(nx, ny):
    """How many candidate (r1,r2,r3) triangles enumerate_triangles considers
    in total (i.e. before its own wrap-exclusion check) vs. how many survive
    it -- reimplements the same r1-centered walk without the final `wraps`
    filter, purely to count what gets thrown away."""
    lat = RhombileLattice(nx, ny)
    bond_by_pair = {frozenset((b["i"], b["j"])): b for b in lat.bonds}
    total = 0
    for m in range(lat.ny):
        for n in range(lat.nx):
            r1_idx = lat._site_index(n, m, "r1")
            r1_pos = lat._position(n, m, "r1")
            nbrs = []
            for other_idx, b in lat.neighbor_table()[r1_idx]:
                other_pos = b["p2"] if b["i"] == r1_idx else b["p1"]
                ang = np.arctan2(*(other_pos - r1_pos)[::-1])
                nbrs.append((ang, other_idx, b))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, b1 = nbrs[i]
                _, idx2, b2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None or closing["type"] != "r2r3":
                    continue
                edges = {b1["type"]: b1, b2["type"]: b2}
                if "r1r2" not in edges or "r1r3" not in edges:
                    continue
                total += 1
    kept = len(enumerate_triangles(lat))
    return total, kept


def main():
    total, kept = count_all_vs_enumerated_triangles(NX, NY)
    print(f"[structural] {NX}x{NY} lattice: {total} elementary triangles total, "
          f"{kept} kept by enumerate_triangles ({total - kept} excluded for touching "
          f"a wrap bond).")
    print("  -> every dual-graph string/loop route (Case0, Case1, Case2, and every "
          "D2/D3/finite-T sweep built on build_case2) is confined to those "
          f"{kept} triangles and can therefore never use a wrap bond, by construction.\n")

    print("[empirical] re-building each construction and checking its own flipped bonds:")

    lat0, tri0, rho0, sib0, hop0 = build_dual(NX, NY)
    corners = hexagon_corners()
    touched0 = apply_closed_loop(lat0, tri0, rho0, sib0, hop0, corners)
    print(f"  Case0 (closed hexagonal loop): {len(touched0)} bonds touched, "
          f"uses_periodic_wrap={uses_periodic_wrap(touched0)}")

    bottom = sorted(range(len(corners)), key=lambda i: corners[i][1])[:2]
    i_left, i_right = sorted(bottom, key=lambda i: corners[i][0])
    p_left, p_right = corners[i_left], corners[i_right]
    lat1, tri1, rho1, sib1, hop1 = build_dual(NX, NY)
    _, _, nodes1, bonds1 = route_string_between_points(rho1, sib1, hop1, p_left, p_right)
    flipped1 = apply_dual_string_defect(tri1, nodes1, bonds1)
    print(f"  Case1 (direct edge): {len(flipped1)} bonds touched, "
          f"uses_periodic_wrap={uses_periodic_wrap(flipped1)}")

    lat2, touched2, corners2, _ = build_case2()
    print(f"  Case2 (5-segment detour): {len(touched2)} bonds touched, "
          f"uses_periodic_wrap={uses_periodic_wrap(touched2)}")

    print("\nConclusion: every matching-principle comparison made so far (Case1 vs "
          "Case2, and everything downstream of it) took place entirely within the "
          "trivial (non-winding) flux sector. This does not contradict any earlier "
          "result -- it just means the winding-sector question was never actually "
          "exercised by any test run before now.")


if __name__ == "__main__":
    main()

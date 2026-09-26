"""Pure construction feasibility test, no SA: can we place isolated
(non-adjacent) monomer defects at a target density p via the project's
own augmenting-path machinery, with everything else automatically
re-matched into rhombi/dimers?

Motivation (from conversation): tuning the density of frustrated
triangles continuously between the dilute limit (D3*=2, isolated
defects) and the fully-frustrated limit (D3*=3, every triangle a
monomer) requires GENUINELY ISOLATED monomers at intermediate density,
not domino-pairs -- randomly activating a random subset of rhombus
diagonals always creates monomers in adjacent PAIRS (one j23 bond is
shared by both triangle-halves of its own rhombus), which is a
different (and less interesting) sweep.

First attempt used the naive straight-line apply_string_defect with
random endpoints and found it frequently produces 0 or 1 (not 2)
frustrated triangles -- turns out that method is only reliable along
specially-chosen lines (distance_scaling_demo.py's own comment: "a
different height silently produced a non-canonical...conflict graph").
The project's actual robust tool for arbitrary-endpoint 2-defect
placement is the dual-graph router (route_string_between_points +
apply_dual_string_defect, via build_triangle_hop_graph), used throughout
Part IV/VI for exactly this reason ("without spurious mid-path
defects"). This version uses that instead.

Open question this answers: is there a maximum achievable density,
since isolated monomers must form an independent set AND the remaining
graph must still admit a perfect matching (Hall's theorem)? p=0.5 is a
natural test point since it's exactly the density of one whole bipartite
class (all "up" or all "down" triangles) on this honeycomb-like
triangle-adjacency graph -- at that extreme, no perfect matching of the
remainder can possibly exist, so construction is expected to fail well
before reaching it.
"""
import numpy as np

from rhombile_lattice import (
    RhombileLattice, enumerate_triangles, build_rhombi, build_triangle_hop_graph,
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
)

NX, NY = 10, 10
TARGET_DENSITIES = [0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
MAX_ATTEMPTS_PER_PAIR = 30
SEED = 7

FAIL_COUNTS = {"routing_error": 0, "count_mismatch": 0, "not_two_new": 0, "adjacency": 0}


def triangle_key_set(tri):
    return frozenset((tri["r1"], tri["a"], tri["b"]))


def adjacent(tri1, tri2):
    """Share any corner site -- catches both true rhombus siblings and
    same-hub neighbors, i.e. anything close enough to plausibly
    domino-share a single state-3 site."""
    return len(triangle_key_set(tri1) & triangle_key_set(tri2)) > 0


def random_point(nx, ny, rng):
    from rhombile_lattice import A1, A2
    fx, fy = rng.uniform(0, nx), rng.uniform(0, ny)
    return fx * A1 + fy * A2


def try_place_one_pair(lat, tri, rho, sib, hop, rng, existing):
    for _ in range(MAX_ATTEMPTS_PER_PAIR):
        a = random_point(NX, NY, rng)
        b = random_point(NX, NY, rng)
        try:
            _, _, nodes, bonds = route_string_between_points(rho, sib, hop, a, b)
        except ValueError:
            FAIL_COUNTS["routing_error"] += 1
            continue
        touched = apply_dual_string_defect(tri, nodes, bonds)
        ft = frustrated_triangles(lat)
        if len(ft) != len(existing) + 2:
            FAIL_COUNTS["count_mismatch"] += 1
            for b_ in touched:
                b_["J"] = 0.0 if b_["J"] != 0.0 else -1.0
            continue
        new_tris = [t for t in ft if not any(triangle_key_set(t) == triangle_key_set(e) for e in existing)]
        if len(new_tris) != 2:
            FAIL_COUNTS["not_two_new"] += 1
            for b_ in touched:
                b_["J"] = 0.0 if b_["J"] != 0.0 else -1.0
            continue
        bad = any(adjacent(new_tris[0], e) or adjacent(new_tris[1], e) for e in existing) or \
              adjacent(new_tris[0], new_tris[1])
        if bad:
            FAIL_COUNTS["adjacency"] += 1
            for b_ in touched:
                b_["J"] = 0.0 if b_["J"] != 0.0 else -1.0
            continue
        return True, new_tris
    return False, None


def main():
    lat = RhombileLattice(NX, NY)
    tri = enumerate_triangles(lat)
    rho, _ = build_rhombi(tri)
    sib, hop = build_triangle_hop_graph(tri)
    n_tri = len(tri)
    print(f"N_sites={lat.n_sites}  N_triangles={n_tri}", flush=True)

    rng = np.random.default_rng(SEED)
    existing = []
    for target_p in TARGET_DENSITIES:
        target_k = int(round(target_p * n_tri))
        target_pairs = target_k // 2
        stalled = False
        while len(existing) < target_pairs * 2:
            ok, new_tris = try_place_one_pair(lat, tri, rho, sib, hop, rng, existing)
            if not ok:
                stalled = True
                break
            existing.extend(new_tris)
        achieved_p = len(existing) / n_tri
        print(f"target p={target_p:.2f}  achieved p={achieved_p:.4f}  "
              f"n_monomers={len(existing)}  stalled={stalled}", flush=True)
        if stalled:
            print(f"  -> stopped: could not place another isolated pair after "
                  f"{MAX_ATTEMPTS_PER_PAIR} attempts at density {achieved_p:.4f}", flush=True)
            break

    print(f"\nFinal: {len(existing)} isolated monomers placed "
          f"(density {len(existing)/n_tri:.4f}) before construction stalled.", flush=True)
    print(f"Failure breakdown: {FAIL_COUNTS}", flush=True)


if __name__ == "__main__":
    main()

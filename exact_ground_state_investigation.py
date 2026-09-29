"""Investigation into the exact D3-dependent ground state energy, beyond
D3=0 -- STATUS: UNRESOLVED. Kept here (not polished into a demo) so the
next session can pick up exactly where this one stopped, including the
dead ends.

## What's solid

- D3=0: E=0 exactly, for any string/matching configuration tried so far.
  Verified by explicit construction (min_state2_assignment) AND
  independently by SA. Provably optimal since 0 is the energy floor.
- All the purely structural findings (face sizes, defect counts, L,
  matching-driven cancellation) don't depend on spin optimization at all
  and remain correct regardless of anything below.

## What's NOT settled: the D3>0 curve shape

Timeline of this investigation, in order, so the dead ends are visible:

1. `energy_via_mincut` / `state2_and_violations_via_mincut` /
   `full_state_via_mincut` (in rhombile_lattice.py) fix r1 sites to
   Potts state 0 always, letting only rim sites choose between 1 and 2,
   then solve *that restricted problem* exactly via min-cut. This gives
   E(D3) = min(D3*L/2, L) for a simple string (kink at D3=2).

2. Checked whether letting r1 ALSO vary ever helps. For a plain
   straight single string (L=14), an exhaustive 3^15 search over the
   sites touching the string confirmed the r1-fixed restriction is
   *not* limiting anything -- BUT that exhaustive search itself only
   covered sites near the string, not the whole lattice, which turned
   out to matter (see next point).

3. For the bent detour configuration, an unrestricted BFS 2-coloring
   over the WHOLE active-bond graph found E=15 at large D3 -- much
   better than energy_via_mincut's claimed 53. This showed the
   r1-fixed restriction *does* matter once the configuration isn't a
   simple straight path.

4. Attempted a principled fix via the classical Toulouse (1977) /
   Bieche-Maynard-Rammal-Uhry (1980) / Barahona (1982) theory: a
   uniform-AF planar Ising ground state's violated-bond count equals a
   minimum-weight matching of the "frustrated plaquettes" (odd faces)
   in the planar dual graph; for exactly 2 defects this is just the
   shortest path between them in the dual. Implemented
   `trace_faces_with_edges` + dual BFS distance. This gave much SMALLER
   numbers (7 instead of 14 for the simple string, 12/13 instead of
   16/53 for shortest-path/detour, 6/6 instead of 29/9 for
   natural/rematch pairing) -- and for D3=0 and small-D3 comparisons
   this looked consistent.

5. **Tried to verify step 4 by actually reconstructing the coloring
   that achieves that distance, and failed twice**:
   - Removing the K matched-path bonds and 2-coloring the rest via BFS
     left 37 bonds still violated (should be 0 if the theory was
     applied correctly here).
   - A spanning-tree BFS 2-coloring of the whole active-bond graph gave
     21 violated bonds -- worse than even the original r1-fixed
     estimate (14), let alone the claimed 7.
   Neither reconstruction attempt is a rigorous disproof of the
   distance value (7 could still be a valid lower bound that neither
   attempt achieved), but it means step 4's numbers are NOT verified
   constructions, just a distance calculation whose achievability is
   unconfirmed.

6. Cross-checked against simulated annealing (single-site heat bath,
   30 restarts, 3000 sweeps, T from 5.0 to 1e-5) for the simple string
   at D3 = 0.5, 1.0, 1.5, 2.0 -- ALL FOUR POINTS NOW COMPLETE:
     D3=0.5  SA=3.500   original=3.500   corrected=3.500   (tie)
     D3=1.0  SA=7.000   original=7.000   corrected=7.000   (tie)
     D3=1.5  SA=10.500  original=10.500  corrected=7.000   -> ORIGINAL
     D3=2.0  SA=14.000  original=14.000  corrected=7.000   -> ORIGINAL
   D3=1.5 and 2.0 are the only two points where the formulas actually
   differ (2.0 is also the original formula's saturation kink), and SA
   sided decisively with the ORIGINAL (hub-fixed) formula at BOTH,
   matching it to 3 decimal places every time. Combined with point 5's
   two failed reconstructions of the "corrected" formula, this is no
   longer just "weighs against" -- for the plain straight-string case,
   the original hub-fixed formula should now be treated as CONFIRMED
   and the Toulouse-distance "correction" as DISPROVEN (not just
   unverified: it was actually tested against ground truth and lost).

## Where this leaves things

Resolved for the single straight string: `energy_via_mincut`'s formula
E(D3)=min(D3*L/2, L) is correct there, confirmed independently by SA at
every tested point including the kink. Still open for anything bent or
branched: the detour counterexample (point 3 above) proves the same
hub-fixing restriction is NOT always exact once the path isn't a
straight line (15 achievable vs. 53 claimed), and no correct general
method has been found yet -- the Toulouse/Barahona attempt that would
have been the natural fix is now doubly discredited (two failed
reconstructions, AND disproof by SA on the one case where it made a
distinguishable prediction). What's missing is a reconstruction
algorithm for the bent/branched case that's actually correct (the
Toulouse/Barahona theory itself is standard and not in doubt; the bug
is somewhere in translating "shortest path in the dual" into "which
vertices get which spin", not in the theorem), or a smarter search
(cluster-move SA / parallel tempering, or genuinely exhaustive search
on a small-enough full lattice, not just "involved" sites) that can
find the true optimum for a bent configuration.

## Concrete next steps for whoever (human or Claude) picks this up

1. Get one FULLY exhaustive ground truth on a lattice small enough that
   literally every site (not just ones "involved" near the string) can
   be brute-forced -- e.g. a short BENT string (the straight case is
   now closed, see above) on the smallest lattice size RhombileLattice
   will build correctly (try nx=ny=5 or 6) and 3^n over the WHOLE
   lattice if n is small enough, or a smarter DP if not. This is the
   one thing neither side of the investigation actually has for the
   bent/branched case: full-lattice ground truth at an intermediate D3.
2. If reconstructing the Toulouse/Barahona coloring: the standard
   reference construction is worth reading directly from Barahona
   (1982) or Bieche et al. (1980) rather than re-deriving from memory
   -- re-derivation is exactly where this investigation went wrong,
   and it's now disproven (not just unverified) for the straight-string
   case, which is a strong signal the same bug affects the bent case.
3. If falling back to SA: use cluster-type moves (e.g. Wolff-style,
   or explicitly try flipping whole geometric regions at once) since
   single-site heat bath demonstrably cannot find the better solution
   in the one case (detour, unrestricted BFS beat it) where we know
   for certain a better solution exists.

The functions below are kept as working (if not fully verified) code,
not wired into any demo script -- the plain-straight-string case is now
resolved (energy_via_mincut is confirmed correct there, see point 6
above); treat everything involving a bend, branch, or multiple strings
in this file's outputs, and in dual_string_demo/density_demo/
matching_and_excitations_demo/detour_demo's D3>0 sections, as
provisional pending step 1 above.
"""
import numpy as np
from collections import deque
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

from rhombile_lattice import RhombileLattice, apply_string_defect, total_energy


def trace_faces_with_edges(lattice):
    """Planar faces of the active-bond graph, each as its list of bond
    dicts in traversal order. Used to find odd (frustrated) faces."""
    active_bonds = [b for b in lattice.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        i, j = b["i"], b["j"]
        pi, pj = b["p1"], b["p2"]
        adj.setdefault(i, []).append((j, np.arctan2(*(pj - pi)[::-1]), b))
        adj.setdefault(j, []).append((i, np.arctan2(*(pi - pj)[::-1]), b))
    for k in adj:
        adj[k].sort(key=lambda x: x[1])

    def next_in_face(u, v):
        nbrs = adj[v]
        idx_u = next(i for i, (n, _, _) in enumerate(nbrs) if n == u)
        return nbrs[(idx_u - 1) % len(nbrs)]

    visited = set()
    faces = []
    for b in active_bonds:
        for (u, v) in [(b["i"], b["j"]), (b["j"], b["i"])]:
            if (u, v) in visited:
                continue
            face_edges = [b]
            visited.add((u, v))
            cur_u, cur_v = u, v
            steps = 0
            while True:
                nxt, _, bond2 = next_in_face(cur_u, cur_v)
                if (cur_v, nxt) in visited:
                    break
                visited.add((cur_v, nxt))
                face_edges.append(bond2)
                cur_u, cur_v = cur_v, nxt
                steps += 1
                if cur_u == u and cur_v == v:
                    break
                if steps > 400:
                    break
            faces.append(face_edges)
    return faces


def odd_face_dual_distance(lattice):
    """UNVERIFIED beyond its distance value -- see module docstring
    point 5. Returns the dual-graph shortest-path distance between the
    (assumed exactly 2) odd faces. Do not treat this as a confirmed
    achievable energy without an accompanying working reconstruction."""
    faces = trace_faces_with_edges(lattice)
    bond_to_faces = {}
    for fi, face in enumerate(faces):
        for b in face:
            bond_to_faces.setdefault(id(b), []).append(fi)
    dual_adj = [[] for _ in faces]
    for bid, flist in bond_to_faces.items():
        if len(flist) == 2:
            f1, f2 = flist
            dual_adj[f1].append(f2)
            dual_adj[f2].append(f1)
    sizes = [len(f) for f in faces]
    outer = int(np.argmax(sizes))
    odd_faces = [fi for fi, face in enumerate(faces) if len(face) % 2 == 1 and fi != outer]
    if len(odd_faces) != 2:
        raise ValueError(f"expected 2 odd faces, got {len(odd_faces)}")
    s, t = odd_faces
    dist = {s: 0}
    q = deque([s])
    while q:
        u = q.popleft()
        for v in dual_adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist[t]


def spanning_tree_coloring_energy(lattice, D3):
    """A cheap (not optimal) 2-coloring baseline, for comparison: BFS
    over the whole active-bond graph from an arbitrary site. Useful as
    a sanity floor, not as a claimed ground state."""
    active_bonds = [b for b in lattice.bonds if b["J"] != 0.0]
    adj = {}
    for b in active_bonds:
        adj.setdefault(b["i"], []).append(b["j"])
        adj.setdefault(b["j"], []).append(b["i"])
    color = {}
    start = next(iter(adj))
    color[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in color:
                color[v] = 1 - color[u]
                q.append(v)
    states = np.array([color.get(s, 0) for s in range(lattice.n_sites)])
    return total_energy(lattice, states, D=(0.0, 0.0, D3))


if __name__ == "__main__":
    lat = RhombileLattice(10, 10)
    apply_string_defect(lat, np.array([0.5, 2.1]), np.array([7.5, 2.1]), wrap=False)
    print("dual-distance (UNVERIFIED as achievable):", odd_face_dual_distance(lat))
    print("spanning-tree baseline at D3=1e6:", spanning_tree_coloring_energy(lat, 1e6))
    from rhombile_lattice import energy_via_mincut
    print("original r1-fixed energy_via_mincut at D3=1e6:", energy_via_mincut(lat, 1e6))

# Next research directions (captured mid-session, not yet implemented)

## 1. Finite-temperature entropic force / possible confinement-deconfinement crossover

So far every D3>0 result in this project (energy_via_mincut, all the SA
cross-checks, Case0/1/2) is a **T=0 ground-state** calculation. The
original literature motivation (Moore, Nordahl, Minář, Shalizi 1999)
was about a **distance-independent entropic force** between defects in
a *different* AF Potts model -- a finite-T free-energy statement, not a
ground-state one. We never actually tested whether that shows up here.

Reasoning for why it plausibly could:
- Ground-state energy of a defect pair: E(L) ~ (D3-dependent slope) * L,
  L = matching distance (established repeatedly: L/2 result, Case1/2
  matching principle, etc.) -- an extensive "string tension" epsilon(D3).
- The number of distinct domain-wall paths of length L between two
  fixed points also grows extensively with L (self-avoiding-walk-like
  branching in the triangle-hop graph), giving an entropy S(L) ~ L * log(mu).
- Free energy F(L,T) = E(L) - T*S(L) = [epsilon - T*s] * L. The
  effective tension sigma(T) = epsilon - T*s could shrink with T and
  cross zero at some T* -- beyond which defect pairs would effectively
  stop caring about separation (entropy wins), matching the Moore et
  al. picture.
- Direct analogy: BKT vortex-antivortex unbinding in the 2D XY model,
  confinement/deconfinement in gauge theories (Wilson loop area law ->
  perimeter law), polymer looping entropy.

**Not yet computed.** Would need genuine finite-T equilibrium sampling
(not annealing to T=0) of F(L,T) for several L at several T, extracting
sigma(T) = d F/dL, and checking whether/where it changes sign. This is
a different kind of calculation than anything done so far in this
project (everything so far anneals to T=0).

## 2. Breaking the D1=D2 symmetry (currently D=(0,0,D3) always)

Everything so far uses D1=D2=0, D3>0. The pristine lattice's active
bonds (r1-r2, r1-r3; r2-r3 off) form a strict BIPARTITE graph between
{r1} (1/3 of sites) and {r2,r3} (2/3 of sites), so *any* assignment
where each site differs from its neighbors gives zero coupling energy
-- the only thing that matters then is the D-term, which is purely
extensive per site. Consequence: the ground state puts the cheaper
state on the MAJORITY sublattice (rim, 2/3) and the costlier one on
the MINORITY sublattice (r1, 1/3), since paying a per-site cost on
fewer sites is always better. With D1=D2=0 this cost is zero either
way, which is an accidental (D1=D2) degeneracy -- and *that* degeneracy
is exactly why Case 0 (the closed loop) could invert the r1<->rim
checkerboard for free: both choices cost nothing.

**If D2 > D1** (breaking that degeneracy, D3 free to be anything):
- The uniform ground state becomes UNIQUE: r1 -> the D2-costed state,
  rim -> the D1 (cheaper) state, always (paying D2 on the 1/3-minority
  sublattice minimizes total cost).
- A Case-0-style closed loop that inverts this assignment inside the
  loop no longer costs zero: every rim site inside now pays D2, i.e. an
  **area-law (bulk) cost proportional to the enclosed region**, not
  just a perimeter-law cost from the domain wall itself.
- This turns the whole "closed loop" problem into an **Ising-droplet
  nucleation problem**: perimeter tension sigma*L vs. enclosed-area
  bulk cost Delta_f * A, with the usual competition (critical droplet
  radius R* ~ sigma / Delta_f, Wulff-construction-shaped droplets,
  etc.) -- a substantially richer structure than anything explored so
  far, and a natural way to make (D1, D2, D3) into a genuine 2-parameter
  (differences from a common baseline) phase diagram instead of a
  single D3 sweep.

**Not yet computed at all.** Everything in rhombile_lattice.py currently
hardcodes D=(0,0,D3) at every call site; would need D1/D2 threaded
through as free parameters and the whole D3=0 zero-energy story
(min_state2_assignment, the bipartition machinery) re-derived for
general (D1,D2) since it currently silently assumes D1=D2=0.

## 3. Does the direct<->detour transition point depend on the detour's enclosed area/curvature?

Raised by the user while re-verifying the D2>0 sweep (see below): Part VI's
D1=0, D3=1, D2-sweep threshold sits at |D2|~1 for the specific hexagon
used throughout this project (radius 3.5, ~720-site lattice). If the
droplet-nucleation framing (perimeter tension sigma*L vs. enclosed-area
bulk cost Delta_f*A, critical radius R* ~ sigma/Delta_f) is the right
picture, a hexagon of DIFFERENT size or a detour of different SHAPE
(different curvature, not just area) should in general shift where the
"cheaper to swap the whole enclosed region" transition happens, since a
bigger enclosed area amplifies the D2-driven bulk term relative to the
fixed local-defect (D3) cost near the two real endpoints.

Concretely worth testing later (not urgent, user explicitly said this
can wait): build 2-3 different hexagon sizes (or non-hexagonal detour
shapes with different curvature at fixed enclosed area) and re-run the
same D1=0, D3=1, D2-sweep methodology on each, checking whether the
direct-vs-detour crossover (and the D2~D3 collapse threshold from Part
VI) shifts with area/curvature or stays pinned at the same |D2|.

## Correction in progress (as of the latest session): the D2>0 "gradual precursor" may have been an SA artifact

D2_sweep_positive.py's original coarse grid (step 0.2) reported n_state3
growing gradually 3->11->15->17 as D2 went from 0 to 0.9, which made it
into the artifact's Part VI. A finer re-sweep (D2_sweep_fine_positive.py,
step 0.025, MORE SA effort) found n_state3 stays flat at 3 across the
ENTIRE D2=0.00-0.35 range, only stepping to 4 by D2=0.4 -- and critically,
this run's D2=0.4 energy (118.8) is LOWER than the original coarse run's
own D2=0.4 result (119.6, n3=15), meaning the original run was stuck in a
worse local optimum, not at the true ground state. The true D2>0 trend
looks much more like a long flat plateau with small discrete steps than
a smooth gradual ramp. D2_sweep_fine_positive_part2.py is re-verifying
0.4-1.0 with the same (extra) rigor; the D2<0 branch (D2_sweep_gradual.py's
original three-stage-collapse claim, jump at -0.6, full collapse by -1.1)
needs the same re-verification treatment next -- flagged by the user,
not yet started. The artifact's Part VI text will need correcting once
both branches are re-verified.

## Context: is this publication-worthy?

Discussed with the user: the D3=0 endpoint, the zero-energy closed-loop
theorem, and the matching-principle confirmation into the non-bipartite
D3>0 regime are solid, but two things are missing for a real paper: (a)
a proven/exact algorithm for the general D3>0 ground state (currently
only SA + hand-picked geometries), and (b) an actual finite-T phase
transition story. Both open directions above are exactly aimed at (b)
and would also enrich (a)'s scope quite a bit. Not urgent to implement
immediately, but should not be lost -- hence this file.

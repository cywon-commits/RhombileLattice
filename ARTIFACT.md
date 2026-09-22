# Shared research write-up

The narrative, presentation-quality write-up of this project lives in a
Claude artifact (not in this repo), kept up to date as the research
progresses:

https://claude.ai/artifact/9jxHSSeFMssZzDUT7Ph4eS

As of this commit it is titled **"Confinement, Matching & Symmetry
Breaking on a Tunable Potts–Ising Lattice"** (v6 restructure) and is
organized thematically rather than chronologically:

1. Why this is worth doing (novelty/significance framing)
2. The model
3. Part I — the exact D3=0 limit
4. Part II — a single string at D3>0
5. Part III — the hub-fixed trap (methodological lesson); now includes
   a tested greedy state-3-removal heuristic (convert each D3<2 ground
   state's 3rd-state site to whichever of state0/state1 has fewer
   active-bond neighbors, to guess the D3=infinity ground state): exact
   when there's no background-swap freedom (a simple string), but fails
   on Case 2's bent detour (greedy E=9 vs true optimum E=7) — pinned
   down to a reversed global antiphase-swap-direction mismatch, not a
   mutual-adjacency conflict, sharpening what a correct general
   algorithm still needs to get right
6. Part IV — closed loops resolve it: the matching principle survives
   into the non-bipartite D3>0 regime (Case 0/1/2)
7. Part V — finite/zero temperature, now split into two distinct,
   correctly-separated mechanisms: (a) Moore et al.'s (1999) actual
   claim — a T=0 ground-state-ENTROPY force within D3=0's own
   zero-bare-energy manifold — directly replicated here for the first
   time (reference-subtracted, windowed, finite-size-checked NX=40 vs.
   NX=80), confirming a real, monotonically-growing, non-saturating
   signal out to defect separation L=36 (an earlier apparent plateau at
   L=16-20 turned out to be a small-box boundary artifact); (b) a
   separate finite-T confinement crossover traced to O's own hub/rim
   order-disorder transition (exact dice-lattice Ising T_c=1.2027
   theory, empirical T_c=0.79 at D3=1 rising to ~1.07 at D3=3 as thermal
   3rd-state population shrinks), explained by the Kadanoff-Ceva
   disorder-operator / Z2 gauge theory duality (this project's "strings"
   are disorder lines), plus an exact T->infinity endpoint (dA=dB=1/3,
   a pure bond-counting fact) and the original D3=2 threshold
   broadening result. Two capstones added: (i) genuine finite-size
   scaling (Binder cumulant + order-parameter sharpening on square
   L=6-14 boxes) now confirms this is a real thermodynamic-limit
   transition and extrapolates T_c(L->infinity) to 0.72 (D3=1),
   ~1.05-1.11 (D3=3), and 1.206 at a calibrated D3=10 — converging
   essentially exactly onto the theoretical 1.2027, with a T=1.5
   calibration sweep showing D3~10-12 (not the initially-guessed D3=3)
   is needed to suppress thermal state-3 occupation below 0.5%; (ii) a
   direct test of whether the D3=2 kink (a dilute-limit combinatorial
   fact) survives in a fully-frustrated lattice (every triangle active,
   the literal triangular-lattice Potts AF) — first reported as a broad
   D3~2-3.5 crossover, then CORRECTED by heavy re-verification: the
   true transition is sharp (not broad) and sits at an exactly
   analytically predicted D3=3.0 (the level crossing between the
   proper-3-coloring line 144*D3 and Wannier's fixed N_bonds/3=432,
   exact once every site's coordination is 6 with j23 active) — the
   original light scan's whole [2.0,3.0] region was under-searched,
   catastrophically so at D3=2.75 (E=533 reported vs. the true 396);
   heavy SA now matches the exact two-line baseline at every point
   checked. Also refined: the T=0 entropic-force functional-form fit
   (5 seeds/point, up from 2) now favors a power law with exponent
   ~0.99 (essentially linear in L) over Moore et al.'s own claimed
   log(r) form (reduced chi^2 0.46 vs 1.16)
8. Part VI — breaking D1=D2 (droplet-nucleation-style symmetry
   breaking); D2 sweeps re-verified as a sharp few-stage staircase on
   both branches (the original coarse sweeps' "gradual precursor" was
   substantially an SA-optimization artifact), a delocalization
   finding for the D2>=1 bulk collapse, and — now resolved — the
   direct<->detour threshold's area/curvature dependence: area
   dependence confirmed (bigger enclosed area -> lower threshold,
   bigger swap), but curvature is NOT independent as originally
   predicted (skewing the detour at fixed area suppresses the swap
   entirely, traced to the detour path's own D3-proportional length
   cost, a three-term not two-term competition)
9. Part VII — flux sectors on the torus: a non-contractible closed loop
   with zero local (frustrated-triangle) defects still costs real energy
   for any D3>0, invisible to both frustrated-triangle counting and
   conflict-graph bipartiteness (ties into Thurston's flip-connectivity
   theorem and Kasteleyn/Fisher–Wu's four-Pfaffian toroidal dimer count).
   Follow-up added: the flux group is confirmed Z2×Z2 (winding the same
   generator twice cancels; (2,1)≡(0,1), (3,1)≡(1,1)), and the matching
   principle itself survives within a fixed nontrivial sector (two very
   different constructions of the same (1,1) sector converge under SA) —
   plus a literature check (Kenyon et al.) on how flux-sector probability
   scales with torus size under random dimerization, and a discrete
   height-function/Burgers-vector formalism making the "dislocation"
   language literal (the winding loop is a genuine screw dislocation,
   B_x=+3 on every row; an open string's two endpoints are a genuine
   dislocation dipole, +3/-3).
10. Experimental and real-material connections (artificial spin ice
    dislocations, NiO domain-wall pinning, delafossite/kagome endpoint
    realizations — with an explicit caveat on what's analogy vs. what's
    the same physics)
11. What's next
12. References

Shared with the whole organization (set via the artifact's own Share
menu, not something this repo controls).

The artifact's version history retains the earlier, chronologically-
ordered notes (every dead end and correction as it happened), which
this v2 restructure intentionally compresses into a forward-looking
summary. `research_notes_next_directions.md` in this repo has more
detail on the two research directions (finite-T, D1≠D2) than fits in
the artifact's "What's next" section.

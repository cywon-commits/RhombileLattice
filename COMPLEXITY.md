# Handoff: computational complexity of the D3-interpolated ground-state problem

This is a self-contained starting point for a side project, separate from
Part IX. The deliverable is a **proof**, not a simulation. Numerics serve
only as a cross-check.

## The problem

Given a graph G = (V, E) and D3 ≥ 0, minimise over 3-colourings
s : V → {1, 2, 3}:

```
cost(s) = #{ (i,j) in E : s_i = s_j }  +  D3 · #{ v : s_v = 3 }
```

This is the ground-state energy of the antiferromagnetic 3-state Potts model
with single-ion anisotropy D = (0, 0, D3) on the "active-bond" graph of the
RhombileLattice project. Colour 3 is **not** vertex deletion: two adjacent
colour-3 vertices form a violated edge and cost 1. (An earlier message
framed it as odd-cycle transversal; that was imprecise.)

Main interest: **planar G of maximum degree Δ.** The project's lattices are
subgraphs of the triangular lattice (planar, Δ ≤ 6).

## Known endpoints

- **D3 → ∞.** Colour 3 is never used, so the problem is the planar Ising
  (max-cut) ground state, polynomial via T-join/matching on the dual
  (Hadlock 1975; Bieche–Maynard–Rammal–Uhry 1980; Barahona 1982).
- **D3 → 0.** "Cost 0?" asks whether G is 3-colourable, which is
  NP-complete for planar graphs even with Δ = 4 (Garey–Johnson–Stockmeyer
  1976). Note this uses D3 → 0 (or D3·n < 1), not a fixed D3.
- **General (non-planar) graphs.** Every fixed D3 is covered by the
  valued-CSP dichotomy (Thapper–Živný; Kolmogorov–Krokhin–Rolínek); max-cut
  alone is already NP-hard there. Planarity is outside that framework, so
  the planar case is the open one.
- **Barahona (1982)** also showed that the planar Ising spin glass **in a
  magnetic field** is NP-hard. Since D3 acts like a field on the third
  state, this is a likely source of hardness gadgets.

## Lemma (proved; please write it up formally and double-check)

**Statement.** If D3 ≥ ⌊Δ/2⌋, some optimal colouring uses no colour 3, so
the problem reduces to planar max-cut and is polynomial for planar G.

**Proof.** Take an optimal s with a colour-3 vertex v. Let v have a
colour-3 neighbours and n₁ + n₂ = deg(v) − a neighbours in colours 1 and 2.
Recolour v with the colour c ∈ {1, 2} that is rarer among its neighbours.
The cost changes by

min(n₁, n₂) − a − D3 ≤ ⌊(Δ − a)/2⌋ − a − D3 ≤ ⌊Δ/2⌋ − D3 ≤ 0.

No new colour-3 vertex is created, so repeating this removes every colour-3
vertex without increasing the cost.

**Per-vertex version.** A vertex of degree d never needs colour 3 once
D3 ≥ ⌊d/2⌋.

## Evidence from the project that the bound is sharp / relevant

- **Fully frustrated triangular lattice** (all bonds active, Δ = 6, bound 3).
  The exact kink sits at D3* = 3. It is the level crossing between the proper
  3-colouring (E = (N/3)·D3) and Wannier's 2-colouring (E = N_bonds/3 = N),
  verified by heavy SA; see `triangular_lattice_sharp_transition_check.py`
  and `triangular_lattice_sharp_kink_figure.png`. So the lemma is tight
  there.
- **Dilute defects on the rhombile lattice.**
  - Frustrated triangles are paired by strings; D_c ≈ 2 at densities
    p = 0.04–0.12 (`worm_Dc_envelope.py`, `Dc_vs_density_figure.png`).
  - The rim sites on strings have degree 4–5, so the per-vertex bound is 2,
    which is consistent.
  - The empirical "absorption rule" (D_c = 1 for a length-1 string, 2 for an
    interior domino, 3 at junctions/saturation, 2L/(L+1) for odd length L)
    is always ≤ ⌊deg/2⌋.
- **Exact D3 = ∞ solver** for these lattices: `exact_Einf.py` (T-join +
  Z2×Z2 torus homology, certified by an explicit 2-colouring; MILP for many
  defects). It can serve as the polynomial-side reference implementation.

## Conjecture (the research question)

For planar graphs of maximum degree Δ, the problem is polynomial for
D3 ≥ ⌊Δ/2⌋, and NP-hard for every fixed D3 < ⌊Δ/2⌋.

- If true, it is a clean dichotomy with a physical knob.
- If false, i.e. a polynomial window exists below ⌊Δ/2⌋, that is the more
  surprising result.
- Natural first case: **Δ = 4, D3 ∈ (0, 2)**. Also worth checking are
  Δ = 3 (planar cubic, bound 1; 3-colourability is easy by Brooks'
  theorem, so the low-D3 side may be polynomial too) and Δ = 6.

## Suggested plan

1. **Literature check.** Look for:
   - planar valued CSPs and planar Max-3-Cut / Max-k-Cut
   - Potts ground-state complexity on planar graphs
   - "independent odd cycle transversal" (e.g. results on planar and
     fullerene graphs)
   - complexity of planar Ising with fields (Barahona 1982)
   - list-colouring / weighted-colouring variants with a penalised colour class

   Decide whether the conjecture is already known.
2. **Write up the lemma formally**, with tightness examples (the triangular
   lattice, and small explicit planar graphs where colour 3 is needed for
   every D3 < ⌊Δ/2⌋).
3. **Attempt NP-hardness for fixed D3 < ⌊Δ/2⌋.** Start with Δ = 4. Candidate
   sources: planar 3-colouring with Δ = 4, planar Ising in a field, or planar
   vertex cover / independent set. Design gadgets that force or penalise
   colour 3 locally.
4. **Numerical sanity checks.** Exact ILP or brute force on small planar
   graphs:
   - Map where colour 3 appears in optima as a function of D3.
   - Search for instances just below ⌊Δ/2⌋ that need colour 3.
   - Test any gadget exhaustively.

   These checks do not prove complexity.

## Risks

- The result may already be known, or may collapse to the "obvious"
  threshold ⌊Δ/2⌋ with a routine hardness proof. That still gives a modest
  theorem.
- Gadget proofs are error-prone. Verify every gadget exhaustively by
  computer.
- **Scope.** This is essentially independent of the rest of the repository's
  physics (the project's lattices are 3-colourable, so hard instances never
  arise there). Beyond the tightness examples above, do not pull in Parts
  I–IX.

## Practical notes

- **Branch.** Work on your own branch `claude/d3-complexity` (create it from
  `claude/modest-feynman-dtkpqn`). Two other sessions push to the parent
  branch; do not push there.
- **SA helpers.** `rhombile_lattice.py` has `simulated_annealing`,
  `total_energy` and lattice builders. `networkx` and `scipy` (MILP via
  HiGHS) are available. numpy/scipy/matplotlib were reinstalled with
  `pip install --user`.
- **Never `pkill -f <pattern>`** when the pattern appears in your own command
  line; kill by PID instead.
- **Containers restart often.** Commit and push frequently.

# Second adversarial review of Theorem C.5 (healing theorem), `paper/healing_proof.tex`

Date: 2026-09-27. Scope: `paper/healing_proof.tex` (primary), `theorem_C.md` §C.5b (Korean notes, Lemma 4 and the
rotated-path exclusion), and the first review `reviews/c5b_review.md`. I made no changes to repository files other
than adding the three `reviews/c5_review2_*.py` scripts and this report.

## Verdict: **SOUND.** I found no fatal error and no gap that needs a new idea.

I re-derived every step independently. Each inference is valid as written, or valid after a one-line justification
that the text leaves implicit (issues 3–9 below). None of them changes the structure of the argument. The rotated-path
case analysis in the chain lemma is exhaustive. The finiteness conclusion re-establishes every hypothesis it needs at
each step.

My computational tests of the intermediate lemmas found nothing. This covers relaxed, non-vacuous versions of the
lemmas (see §3), about 34k + 101k hypothesis-satisfying chain states. They also show that the chain lemma uses
maximality in an essential way, in exactly the places the proof says it does (§3, "instructive non-violations").

The conflict with Johnson et al. (Algorithmica 87 (2025), Thm 11) therefore remains. The logic of this proof holds up
under scrutiny. What still needs outside confirmation is the independent-OCT reduction in that paper; see the existing
`theorem_delta3.md` note on its planarity gap. I recommend telling the authors and getting a second external reader
before submission.

## 1. Step-by-step check

**(F1).** Re-inserting a vertex with at most 3 present neighbours on the minority side creates at most 1 monochromatic
edge. Therefore `fr(G) <= |X| + fr(G-X)` for every X, independent or not. A tight S with `fr(G-S)=0` is then an
independent OCT of size `fr(G)`, and every OCT has at least `fr(G)` vertices. Both `fr` and `ioct` are additive over
components. The empty set is tight and tight sets have `|S| <= fr`, so a maximum tight set exists. ✓

**Lemma aug.** ✓ This needs only that `v` is an endpoint of a monochromatic edge in *some* optimal colouring of `G-S`.

**Lemma stuck (blocked structure), question (d).**
- The chosen optimal colouring is **arbitrary**. The lemma needs only that a frustrated component contains a
  monochromatic edge, which is true in every optimal colouring.
- Every recolouring in the walk keeps the colouring optimal for `G-T`. It changes only K, and an optimal colouring of
  `G-T` is optimal on each component. Lemma aug is therefore available at every step.
- The walk argument is correct. Suppose `v_{i+1} = v_l` with `l < i`. Since `v_i` is adjacent to `v_l`, we get
  `v_i ∈ {v_{l-1}, v_{l+1}}`. This forces either `i = l-1 < l`, which is impossible, or `v_{i+1} = v_{i-1}`, which
  would be backtracking. Backtracking cannot happen because `deg_K v_i = 2` in a simple graph.
- `m >= 1`, so the cycle has at least 3 vertices. `u` has degree 2 in K because it is also an endpoint of the
  monochromatic edge.
- Each cycle vertex is an endpoint of a monochromatic edge in some optimal colouring. So it has at least one
  T-neighbour, and with Δ ≤ 3 it has degree exactly 3 and exactly one T-neighbour.
- The lemma uses only "no augmentation is possible", which is weaker than maximality. ✓

**Lemma exchange, question (a).**
- *Tightness of T'.* `T'` is independent because `y`'s only T-neighbour is `t`. Also `fr(G-T-y) = k-1` because
  `fr(Y) = 1` and `fr(Y-y) = 0`. The vertex `t` has at most 2 neighbours in `G-T'` (T is independent), so re-inserting
  it adds at most 1. Hence `fr(G-T') <= k`, (F1) gives equality, and `|T'| = |T|` makes T' maximum. This step does not
  use maximality of T; I tested it as claim (A) in §3. ✓
- *t adjacent to another odd-cycle component Y2.* The component of `G-T'` containing Y2 is frustrated and has a
  vertex of degree 3, which contradicts Lemma stuck for T'. ✓
- *`fr(G-T') = (k-1) + fr(D)`.* This uses Lemma stuck for T: all k frustrated components of `G-T` are odd cycles with
  `fr = 1`. The other k−1 of them are untouched and are not in D.
  - Degenerate case "D − t bipartite": then `fr(D) = 0` and `fr(G-T') = k-1`, which contradicts (F1). This is covered.
  - Degenerate case "t has neighbours in two different components of G−T−y", for example one in `Y-y` and one in a
    path. Then `D - t` would be disconnected, but `D - t` is a path. This is covered implicitly.
  - Degenerate cases `deg t <= 2` and `|Q| = 0`: excluded, because Lemma stuck for T' forces `deg t = 3` and
    `t1 != t2`. `|Q| = 2` gives a triangle and is allowed.
- *The K4 case.*
  - Apply the same argument to `y+`, with `T'' = T - t + y+`. By the first paragraph T'' is also maximum tight.
  - The component W'' for `y+` contains `y ∈ Y - y+`, hence `W'' = Y - y+`.
  - Its endpoints are `{y, y++}`, and these must equal `t`'s other neighbours `{y, y-}`. So `|Y| = 3`.
  - `Y ∪ {t}` is then a 3-regular K4, hence a whole component. Since G is connected, `G = K4`.
  - The configuration "t adjacent to y and to exactly one non-neighbour of y on Y" is impossible: W would have to be
    `Y - y`, and then t's neighbours in G−T' would be its endpoints `y±`. ✓
- *Otherwise.* W is a component of `G-T-y` other than `Y-y`, hence a component of `G-T` other than Y, and it is a
  path. (a), (b) and (c) follow, including "y is t's only neighbour on Y": all of t's neighbours in `G-T'` lie on D. ✓
- *Disconnected G* is handled by additivity. *Multiple S-neighbours* of a cycle vertex are excluded by Lemma stuck.
  *Degree-2 vertices* can never lie on a frustrated component, be a closing vertex t, or lie on Y. They can appear
  only in untouched components, which the argument never uses. ✓

**Chain lemma, question (b).**
- *Component list (i)–(iii) of `G - T_j`.* The text states it without proof. It is correct and follows from three
  facts:
  - `N(w_i) = {p_i} ∪ ends(R_i)`, by exchange (a),(b) at stage i−1.
  - `N(p_i) ⊆ R_{i-1} ∪ {w_{i-1}, w_i}`, because `Y_{i-1}` is an odd-cycle component, `p_i` has degree 3, and its unique
    `T_{i-1}`-neighbour is `w_i`.
  - The R_i are pairwise disjoint, by (I2).
  - Details:
    - A vertex `r ∈ R_i` other than `p_{i+1}` is adjacent to no `w_m` except `w_i`. If `r ∈ N(w_m)`, then either
      `R_m = R_i` or `p_m ∈ R_i`, which means `m = i+1` and `r = p_{i+1}`.
    - `Z - p_1` is not adjacent to `w_1`, because `N(w_1) = {p_1} ∪ ends(R_1)`.
    - Type-(i) components meet no `N(w_i)`.
  - (I3) is not actually used in the proof of (I2). It is a by-product.
- *`w_{j+1} ∈ S \ W_j`.* ✓ The case i = 1 uses `R_0 = Z` and "R_j ≠ Z" from (I2).
- *Rotated-path exclusion; this answers the open doubt at `theorem_C.md` l.227.*
  - `N(p_i) ⊆ R_{i-1} ∪ {w_{i-1}, w_i}` covers every case: interior p_i, endpoint p_i (which gets `w_{i-1}`), and
    i = 1 (`N(p_1) ⊆ Z ∪ {w_1}`).
  - The three candidates for a non-new `R_{j+1}` are exhaustive, because `R_{j+1}` is a path component of `G-T_j`
    and `Y_j` is a cycle:
    - `R_{j+1} = Z - p_1`: excluded by distinctness of S-neighbours on Z. This needs `|Z| >= 3`, so `z+ != z-`.
    - `R_{j+1} = ρ_i` with `p_{i+1}` interior to R_i: excluded because `S \ W_j ⊆ T_i` and by distinctness at stage i.
    - `R_{j+1} = ρ_i` with `p_{i+1}` an endpoint of R_i: excluded because S is independent. This includes
      `|R_i| = 2`, where ρ_i is the edge `{w_i, other vertex}`.
  - The case i = j cannot occur, because at stage j the cycle `Y_j` is not a path.
  - I tried to build `R_j = Z` or `R_j = R_i`. Each attempt forces two vertices of some `Y_i` to share a
    `T_i`-neighbour. Distinctness at stage i is a consequence of maximality of `T_i` (via the K4 argument), so this
    cannot be arranged. §3 shows concretely that when maximality is dropped, exactly this coincidence occurs.
- *(I3)* ✓

**Finiteness, question (c).**
- The inductive hypotheses at step j+1 are exactly: `T_j` is maximum tight, `Y_j` is an odd-cycle component of
  `G - T_j` (by (I1), from exchange (c)), `R_j` is non-empty, and `G != K4`.
- Distinctness and (a) at the earlier stages i ≤ j come from applying the exchange lemma to `(T_i, Y_i)`, which is
  legitimate by (I1) at stage i.
- No "blocked structure" beyond Lemma stuck for `T_{j+1}` is needed, and that is automatic because `T_{j+1}` is
  maximum tight.
- Base case: Z exists because `k >= 1` and Lemma stuck applies to S. So an infinite chain has pairwise distinct
  `R_1, R_2, ...` in a finite graph. ✓

**Question (e).** Δ ≤ 3 is used in (F1) (at most 1 new monochromatic edge per re-insertion), in Lemma stuck (a T-neighbour
forces `deg_K <= 2`, and "exactly one T-neighbour, degree 3"), and in the exchange (t has at most 2 remaining
neighbours). 3-regularity is never assumed. Where degree 3 is needed, it is *derived*. ✓

## 2. Issues

1. **(none fatal, none major).**
2. **Minor – K4 case, Lemma exchange l.73–74.** "Applying the same argument to y+" hides two steps:
   - `T - t + y+` is again maximum tight (first paragraph of the proof).
   - The component W for `y+` contains `y`, so it equals `Y - y+`.

   *Fix:* add these two clauses.
3. **Minor – Lemma exchange l.70.** "W = D − t is a component of G − T − y" needs the reason: D − t is connected,
   being a path, and it is the union of the components of G − T − y adjacent to t.
4. **Minor – Chain lemma l.97–99.** The list (i)–(iii) of components of `G-T_j` is asserted without proof. *Fix:*
   - Move the neighbourhood sentence (l.99) before the list.
   - Add: "a vertex of `R_i` other than `p_{i+1}` has no neighbour in `W_j` other than `w_i`, since
     `r ∈ N(w_m)` forces `R_m = R_i` or `p_m ∈ R_i`; type-(i) components meet no `N(w_i)`."
   - Note that the list uses disjointness of the R_i, which is (I2) for indices ≤ j.
5. **Minor – Chain lemma, rotated-path case l.108–109.** The contradiction is with the "distinct T-neighbours"
   clause, equivalently with (a) for `y`. It helps to say explicitly that ρ_i's endpoints are the two
   `Y_i`-neighbours of `p_{i+1}`, and that `i = j` does not occur because `Y_j` is a cycle.
6. **Minor – Chain definition l.82.** "(so p_{j+1} ≠ w_j)" is meaningless for j = 0. Uniqueness of the
   `T_j`-neighbour `w_{j+1}` is used before (I1) is stated; say "by (I1) for j and Lemma stuck".
   (I1) should also be stated for j = 0 (S maximum tight, Z an odd-cycle component, by Lemma stuck), since the
   conclusion needs it for the first step.
7. **Presentation – (I3).** It is not used in the proof of (I2). Either drop it, or keep it and say it is a
   by-product. Its second clause is only true because p_{j+1} is not adjacent to R_{j+1}. The proof says this only
   implicitly, through `N(p_{j+1}) ⊆ R_j ∪ {w_j, w_{j+1}}`.
8. **Presentation – Lemma stuck l.44–46.** Replace "since all visited vertices have both neighbours accounted for"
   with the explicit first-repetition argument from §1. The first review asked for this (R2), and the .tex version is
   still terse.
9. **Presentation – Energy paragraph.** The case `0 < D3 <= 1` depends on Theorem sandwich, which is not in this file.
   The D3 ≥ 1 bound applies (F1) to the *state-3 set*, which need not be independent. (F1) allows that, so it is fine,
   but it is worth saying.
10. **Remark (not an error).** The proof is non-constructive. It does not give a polynomial algorithm to *find* the
    independent OCT, only its size. For planar G, the claim "OPT computable in polynomial time" is about the value,
    which is what the paper says. Keep that wording precise.

## 3. Computational checks (all scripts in `reviews/`)

**Why relaxed tests are needed.** Lemmas stuck, exchange and chain all assume "maximum tight set with k ≥ 1". If the
theorem is true, that assumption is vacuous. My exhaustive runs confirm this: through n = 11, every tight set that is
non-augmentable with k ≥ 1 is in K4. So I tested relaxed versions whose proofs go through verbatim.

**`c5_review2_lemmas.py 11`.** All 8,093 connected graphs with max degree ≤ 3 on 3–11 vertices (`nauty-geng -c -D3`),
all independent sets. Runtime about 3 min.

| Check | Cases tested | Violations |
|---|---|---|
| (F1), all independent T | all | 0 |
| (A) exchange keeps tightness, without maximality | 693 swap cases on 75,079 tight sets | 0 |
| Lemma aug, all tight sets | all | 0 |
| Lemma stuck under the weaker hypothesis "tight + non-augmentable" | 4 cases, all K4 | 0 |
| Non-augmentable tight k ≥ 1 in a non-K4 graph | — | 0 |

- **(B)/(C) relaxed hypothesis ("good").** T is *blocked*: every frustrated component is an odd cycle whose vertices
  have degree 3 and exactly one T-neighbour. In addition, every single exchange at every vertex of the current cycle
  gives a blocked set with the same fr. These are precisely the sets on which the proof invokes Lemma stuck.
  - 13,029 blocked (S, k ≥ 1) configurations were tested.
  - `good` held in 1 non-K4 state and 4 K4 states. There were 3 chain steps with the hypotheses satisfied, and 0
    violations.
- **(E) structural chain core.** Hypotheses H_0..H_j are the conclusions of stuck and exchange at each stage, stated
  locally. Checked conclusions: `w ∈ S\W`, `R_{j+1}` is a new component of G−S, `Y_{j+1}` is an odd-cycle component,
  the component list (i)–(iii), and (I3).
  - In this range H_0 holds only 3 times (n = 11), with 0 violations.

**Instructive non-violations: maximality is essential, and it enters exactly where the proof says.** When I first
required only "T and T' blocked with equal fr", I got apparent counterexamples. Both come from dropping the hypothesis
for the *other* exchanges, which the proof really does use.

- **Exchange (a) fails (n = 7).**
  - Edges: `[(0,2),(0,4),(0,5),(1,3),(1,5),(1,6),(2,4),(2,6),(3,5),(3,6)]`, S = {3,4}.
  - Z = 0-2-6-1-5 is a 5-cycle. Vertex 3 is adjacent to three vertices of Z.
  - T' = {1,4} is again blocked with fr 1. But `T - 3 + 5` is not blocked, so the proof's "apply the same argument
    to y+" step is exactly what fails.
- **R_2 = Z − p_1 (n = 8).**
  - Edges: `[(0,3),(0,5),(0,6),(1,4),(1,6),(1,7),(2,5),(3,5),(3,7),(4,6),(4,7)]`, S = {2,3,6}.
  - Z = {1,4,7}, p_1 = 7, w_1 = 3, R_1 = {0,5}, p_2 = 0, w_2 = 6, R_2 = {1,4} = Z − p_1.
  - This is possible only because vertices 1 and 4 share the S-neighbour 6. That violates the distinctness which
    the proof derives from maximality.

Neither is a counterexample to the paper; both violate hypotheses that do hold for maximum tight sets. They confirm
that the "distinct T-neighbours" consequence carries the whole rotated-path exclusion.

**`c5_review2_structural.py`.** Random "tree-of-closers" graphs, 25–90 vertices. Every vertex that needs an
S-neighbour gets a fresh closer of a fresh even path, or with some probability reuses an existing S-vertex with free
degree, adversarially. Random extra edges are added. Every chain is run under H_0..H_j.

| Seed | Graphs | H-satisfying states (by depth 0/1/2/3/4/5) | Chain steps checked | Violations |
|---|---|---|---|---|
| 1 | 20,000 | 14,239 / 16,252 / 3,065 / 264 / 14 / 2 | 114,357 | 0 |
| 2 | 60,000 | 42,769 / 48,545 / 8,990 / 811 / 60 / 1 | 341,975 | 0 |

Chains reached depth 6, so the Z−p_1 and both rotated-path cases were exercised many times, and none ever fired.

**`c5_review2_milp.py`.** Exact MILP (HiGHS) comparison of ioct and fr. Families:
- random cubic graphs, n = 30–80;
- random cubic graphs with gadget substitutions: diamond (K4−e) edge gadgets, prism-minus-rung, K4 with two edges
  subdivided, triangle subdivision, truncation;
- planar duals of random Delaunay triangulations, plain and gadgetized (planarity preserved);
- triangle-rich graphs built from C3/C5/diamond/K4-with-subdivided-edge blocks, plus a planarity-filtered version.

Results:
- 1,331 graphs from seeds 11/12/13 (900 s each), plus 55 from a trial run with seed 7: **1,386 graphs in total**.
- n = 16–100, mean about 55.
- 699 + 30 were planar.
- 2 MILP timeouts, which were skipped.
- **0 mismatches.**

**Counterexample: none found.**

## 4. Relation to the first review

The first review verified the lemmas logically and tested the final theorem and (F1)/aug. This review adds:
- a line-level check of the degenerate cases of the exchange lemma;
- a proof that the (i)–(iii) component list is complete, which the text only asserts;
- confirmation that `N(p_i) ⊆ W ∪ R_{i−1}` covers all cases (the open doubt at theorem_C.md l.227);
- non-vacuous machine tests of the combinatorial core of the chain lemma;
- two explicit examples showing which uses of maximality are indispensable.

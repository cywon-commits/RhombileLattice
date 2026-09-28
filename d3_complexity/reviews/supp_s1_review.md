# Referee report: Supplemental Note S1 (planarity of G_phi in Johnson et al., Thm 11)

Reviewer stance: skeptical, as for a note that criticises published work. Date 2026-09-28.
Files read: `paper/supplementary.tex` (+ PDF), `supp_planarity.py`, `verify_delta3_reduction.py`,
`planarity_gap.py`, `supp_runs/planarity_evidence.json`, `paper/make_supp_figs.py` (header),
`theorem_delta3.md`, `lit_review.md`, the relevant `main.tex` paragraphs (abstract, l. 474-485, l. 559-571).

## Verdict

The core argument holds. Lemma S2 is correct under (F1)-(F3). Every number in Table S1 that can
be checked was reproduced independently, except the |K| column, which depends on the implementation.
The S1.6 argument also holds and does not depend on which planar 3-SAT variant is meant.
I found **no fatal issue**. There are **three major issues**, all fixable in a few lines:
1. A false general inequality in S1.1.
2. An imprecise co-nested statement in S1.6. As worded, it describes an NP-complete class.
3. The "strictest variant" row phi_3 lies in a class that Tovey showed is always satisfiable.

The rest are minor points, and the tone is mostly already modest. I also found new instances
that close the "your variant is not ours" escape routes (Sec. B below). They should be added.

---

## A. Independent recomputation (my own code, scratchpad only; no project file modified)

Script: my own `G_phi` builder, written from the gadget edge list in `theorem_delta3.md`
(c1-a, c2-b, c3-f, a-b, a-d, b-d, d-e, e-f, f-c, e-c, c-P(2j), u_j-P(2j+1) for 2-clauses).
It is isomorphic to `verify_delta3_reduction.build` on all 7 formulas tested. networkx 3.6.1, brute-force SAT.

| inst | sizes | lit max | var occ | both signs | I plain | I split | H plain | H split | G: V,E (Delta) | G planar | #sat assignments |
|---|---|---|---|---|---|---|---|---|---|---|---|
| phi*  | 2 | 2 | 3 | yes | planar | planar | NON | NON | 31,43 (3) | no | 1 |
| phi3  | 3 | 2 | 3 | yes | planar | planar | NON | NON | 40,55 (3) | no | 8 |
| phi23 | 2,3 | 2 | 3 | yes | planar | planar | NON | NON | 41,57 (3) | no | 2 |
| phi0  | 2,3 | 2 | 2,3,4 | **no** (x3 only positive) | planar | planar | NON | NON | 92,130 (3) | no | 3 |

- All table entries agree, apart from the |K| column (issue 10).
- K_{3,3} of phi*: parts {C1,C2,C3}, {p,x,y}. Confirmed.
- phi23 path claim: parts {C1,C3,C4}, {p,x1,x3}, and C4-x1 realised as C4 x2 C2 x1. Confirmed. C4=(-2,3) contains x2 and x3. C2=(1,2) contains x1 and x2. x2 and C2 are not used elsewhere.
- phi3 by hand, which the text says "can be read off" but does not give: parts {C2,C3,C4} and {p,x3,x4}, with C2-x4 realised as the path C2 x2 C1 x4. Suggest adding this.
- "H split without the x-xbar edge" (literal vertices not merged) is **planar for all four instances**. So (F1) is load-bearing: see issue 5.
- Random family (re-ran their generator, seed 1, N=300, in scratch): 131 formulas, 130 non-planar H **both** plain and split (no disagreement), 130 non-planar G. Confirmed.
- `verify_delta3_reduction.py` re-run: 70 sat, 10 unsat, 0 mismatches. Confirmed.
- **Variable cycle (Lichtenstein).** Every table instance also has I(phi) + a cycle through the variables planar:
  - phi3: order 1,2,3,4
  - phi23: order 1,2,3
  - phi0: order 1,2,4,3,5,6
  - phi*: only 2 variables, so the cycle is degenerate (a single edge).
  - The note does not say this. It should (issue 4).

## B. New instances found (all verified: planar I, non-planar H plain, non-planar G, max degree 3, MILP equivalence OK)

1. **phi_B2 (exactly 3 distinct variables per clause, every literal exactly twice, every variable exactly 4 times, 2+2; I + variable cycle planar with order 1,2,4,5,3,6).** Satisfiable (16 assignments). G: 76 V, 109 E.

   `(-5,-4,2)(-6,-3,4)(3,5,1)(-6,3,1)(4,6,-2)(-1,-5,2)(5,-4,-3)(-2,-1,6)`

   This is the strictest variant compatible with a subcubic literal vertex that is not trivially satisfiable.
2. **phi_u (unsatisfiable; clause sizes 2 and 3; lit <= 2; every variable exactly 3 times with both signs; I + variable cycle planar with order 1..5; the literal-split incidence graph and the split variable cycle are also planar).** G: 72 V, 101 E.

   `(-4,-2)(4,5)(3,1)(-2,4,-5)(-5,-3)(2,-1)(1,-3)`

   MILP: OPT(0.3)=4.5 > 2m*0.3=4.2, and OPT(0.9)=13.5 > 12.6.
3. phi_2u (unsatisfiable pure 2-CNF, 4 variables, variable cycle planar). Less useful because it is 2-SAT:

   `(1,-3)(4,-2)(-2,1)(-3,4)(-4,-1)(3,2)`

No unsatisfiable exact-3/literal-exactly-twice instance with non-planar H turned up for up to 12
variables in a random search. Small unsatisfiable (3,B2) formulas are rare, so this is expected.
It is not needed: S1.6 covers every variant anyway.

---

## C. Issues

### 1. [major] S1.1: "Always ioct(G) >= fr(G)" is false for general graphs
- The wheel W4 (= K_{1,2,2}) has ioct = 1 (the hub) but fr = 2. The main text itself uses the fact that a degree-4 site relieves two bonds.
- The inequality holds for maximum degree <= 3: give each s in S the minority colour of its <= 3 neighbours, so it causes at most 1 monochromatic edge.
- A hostile reader will spot this at once.
- **Fix:** "For graphs of maximum degree at most three, ioct(G) >= fr(G) (colour G-S properly and give each vertex of S the minority colour of its neighbours)."

### 2. [major] S1.6: the co-nested class is described imprecisely, and as written it is the wrong (hard) class
- "Adding a cycle through the clause vertices ... keeps the drawing planar. Satisfiability of such formulas ('co-nested') is decidable in polynomial time." A reader will read "such formulas" as "incidence graph plus a clause cycle is planar".
- But, to my recollection, Planar 3-SAT with a clause cycle is NP-complete. Pilz 2019's abstract says its result "complements previous results demanding cycles only through either the variables or clauses". Linked Planar 3-SAT is hard too. Please check this against Pilz 2019.
- What H-planarity actually gives is stronger: the cycle can be added **so that all variable vertices lie on one side of it** (the cycle bounds a face containing no variable). To my understanding, this is what "co-nested" means in Kratochvíl-Křivánek: the dual of Knuth's nested SAT. It is also equivalent to H(phi) being planar.
- **Fix:** state it as: "H(phi) is planar iff I(phi) has a planar drawing in which a cycle through all clause vertices can be added with all variable vertices on the same side. These are the co-nested formulas of [KK93] (for m >= 3)."
- Keep the \verify, and extend it to "confirm that KK93's co-nested = clause cycle with all variables on one side, and that the result covers CNF with clauses of size 2 and 3."
- **Also fix the justification of "adding a cycle keeps the drawing planar".** As written it is only asserted, and a face boundary need not be a simple cycle (clause vertices may be cut vertices and recur on the boundary walk).
- Clean proof: take a planar drawing of H. Replace p by a small circle around p, meeting the spokes p-C_j in the rotation order at p. Take the m crossing points as new vertices and contract each remaining spoke segment. The result is I(phi) plus the cycle C_{pi(1)}...C_{pi(m)}, drawn planar, with the disc containing no variable.
- This avoids face-boundary subtleties entirely. For m <= 2 the claim is trivial.

### 3. [major] S1.5 / Table S1: phi_3 is in a class that is always satisfiable, and "strictest combination we are aware of" is misleading
- Tovey (Discrete Appl. Math. 8 (1984) 85) proved that every CNF with exactly 3 distinct variables per clause and every variable in at most 3 clauses is satisfiable (by Hall's theorem).
- So "exactly three literals per clause + exactly three occurrences per variable + both signs" is not a planar 3-SAT *variant* used for hardness. It is a trivial class, and an author could fairly say the table's "strictest" row is irrelevant.
- With the literal <= 2 constraint needed for a subcubic literal vertex, the only non-trivial exact-3 variants have some variables with 4 occurrences.
- **Fix:** add phi_B2 (Sec. B.1) as the "exactly 3 literals, each literal exactly twice" row. Keep phi_3 only with a sentence such as "(this class is always satisfiable [Tovey]; we include it only because it is an instance of every variant that bounds occurrences by three)", or drop it.
- Replace "the strictest combination we are aware of" with an explicit list of which instance satisfies which condition. \verify which variant Johnson et al. use, and whether that variant is proved NP-complete with clause size exactly 3.

### 4. [major-for-persuasiveness] Coverage of variants: say what the instances satisfy, and add an unsatisfiable one
- The likely author response is "our Planar 3-SAT variant excludes such instances." Three points answer it:
  - (a) S1.6 already answers this for *every* variant. Say so explicitly in S1.5: "The table is illustrative; by S1.6 no NP-complete variant can avoid such instances, unless P = NP."
  - (b) Add a column or sentence: "I(phi) plus a cycle through the variables is planar (Lichtenstein's original planar 3-SAT; equivalently a rectilinear drawing)". This holds for phi3, phi23, phi0, phi_B2 and phi_u. For phi* the cycle is degenerate, so point to phi23 as the smallest non-degenerate example.
  - (c) The caption says "all instances are satisfiable", which invites the objection "only easy instances". Satisfiability is irrelevant to whether planarity holds for all inputs, but add phi_u (unsatisfiable, sizes 2-3, lit <= 2, var = 3 with both signs, variable cycle planar) to close this.
- Also note that phi0 has a pure variable (x3 occurs only positively, twice). The "var <= 4" column hides this; it matters for "both signs" variants.
- **Fix:** add phi_B2 and phi_u rows plus a "variable cycle" column. Drop or rephrase "all instances are satisfiable".

### 5. [minor] Lemma S2: state all the hypotheses actually used; (F1) is essential
The proof also uses disjointness:
- P, the gadgets Q_j, and the literal vertices are pairwise disjoint.
- Distinct gadgets are disjoint.
- u_j is not part of Q_j or P.

State this as part of F1-F3, or as "(F0) these sets are pairwise disjoint". The \verify on F1-F3 should include it.

(F1) is load-bearing. Without merging x and xbar, the apex graph is planar for all four table instances (checked). If the original variable gadget is not a single edge, the lemma still holds provided the vertices representing x_i induce a connected subgraph disjoint from P and from the gadgets. Say so: it makes the lemma robust to small misreadings of Fig. 3.

(F3) is also load-bearing for the all-3-clause instances (phi3, phi_B2). The verify item should check explicitly that **every** clause gadget, not only the 2-literal ones, has an edge to P (the c_j-P(2j) edge in `build`).

"Deleting any further edges and parallel copies leaves exactly H" should read "... and any remaining vertices ...". The original may contain vertices our reconstruction lacks.

Parallel edges and 2-literal clauses are handled correctly. Contraction never merges distinct branch sets, because the sets are disjoint. The gadget {a,b,d,e,f,c} is connected (triangle abd - de - triangle efc). u_j is adjacent to P in `build`, so it is not "private" in the sense of degree 1. This is harmless because it is deleted, but F2 and the Fig. S1 caption ("unused third inputs") should not suggest that it is isolated.

### 6. [minor] S1.6: logical placement of "unless P = NP", and scope of "any reduction"
- "Hence any reduction that attaches one connected reference structure to every clause gadget produces planar graphs only on a polynomially solvable class of formulas, unless P = NP." The containment in the co-nested class is unconditional. Only the consequence "so it cannot be NP-hard" needs P != NP.
- "Any reduction" is too broad without the hypotheses. It needs a connected reference structure disjoint from the gadgets, and connected, disjoint variable parts adjacent to their clause gadgets.
- **Fix:** "For every construction with properties (F0)-(F3), the formulas phi for which G_phi is planar are co-nested, and co-nested satisfiability is in P [KK93]. Hence, unless P = NP, no such construction yields planar graphs on all instances of an NP-complete variant of SAT."
- Also replace "it must fail on hard formulas" (vague) with that sentence.

### 7. [minor] Remark S4 is circular and slightly combative
- It invokes Theorem S1 to dismiss Linked Planar 3-SAT. That only restates S1.1. The more useful, independent observation is this: in Linked Planar 3-SAT, variables lie on both sides of the clause cycle, so H(phi) is generically non-planar and Lemma S2 still applies to a single-reference-path design.
- **Fix:** either say that, or delete the remark. Add \verify for "Linked Planar 3-SAT NP-complete (Pilz 2019, Thm 10?)". The working notes cite Thm 10; that citation is unverified in the tex.

### 8. [minor] Reconstruction vs. original: numbers are ours
- |V|, |E|, |K| and Fig. S1(a) are properties of `build()`, which is our reconstruction. For example, the u_j-P(2j+1) edge and the attachment positions on P may differ from Fig. 3 of the original.
- **Fix:** in the table caption and S1.5, write "G_phi as reconstructed in `build` (the non-planarity conclusion does not depend on the reconstruction, by Proposition S3)".
- S1.2 item 1: "We checked it by exact integer programming" should read "We checked it, for our reconstruction and via the Potts formulation of the main text (OPT_{D3} <= 2m D3), by exact integer programming ...". The MILP does not test IOCT directly.

### 9. [minor] Statements that also need a \verify (beyond the existing red markers)
- S1.2 item 1: "a subcubic graph G_phi with 2m vertex-disjoint triangles such that G_phi has an IOCT of size 2m iff phi is satisfiable". This attributes specific content to their proof before the S1.3 verify appears. Move or extend the verify.
- Remark S4: Pilz's NP-completeness of Linked Planar 3-SAT (theorem number).
- The new Tovey citation (issue 3).
- The comment inside the S1.3 verify, "per our comparison the arXiv v5 text and the journal text of this proof coincide". `lit_review.md` and `theorem_delta3.md` claim the journal text was compared (2026-09-27), but the brief says the journal PDF is not available. Resolve this. If no comparison was actually made, delete the claim before submission.
- S1.2 item 2: add "without K_4 components" for disconnected cubic graphs (currently "other than K_4").

### 10. [minor] |K| column is implementation-dependent
The Kuratowski subgraph returned by networkx depends on vertex insertion order. With my isomorphic but differently labelled graph, networkx 3.6.1 returned 27, 27, 26, 40 vertices, against the table's 26, 25, 33, 55. The number carries no information.
- **Fix:** drop |K|, or say "a K_{3,3} subdivision (found by networkx 3.6; its size is not canonical)". Keep "all are K_{3,3} subdivisions". That holds: each has 6 branch vertices of degree 3, consistent with H being bipartite-plus-apex.

### 11. [minor] S1.7 item 3 numbers are not traceable to the main text
- "1,386 further graphs of 16-100 vertices (main text, Appendix)": 1,386 does not appear in `main.tex`. It comes from `reviews/c5b_review2.md` (seeds 7/11/12/13).
- "251 graphs of 14-60 vertices": the main text says "up to 60".
- **Fix:** either add the 1,386 run to the main-text Appendix with its script, or cite it correctly. Harmonise the vertex ranges.

### 12. [minor] Code/text mismatch on H
- The text defines I(phi) and H(phi) with one vertex per *variable*. `supp_planarity.py` tests the literal-split graphs (x-xbar edge present).
- I recomputed both: they agree on every instance and on all 131 random formulas. Split planarity of I implies plain planarity, but split non-planarity of H does not by itself imply plain non-planarity.
- **Fix:** test the plain H in the script, as the text states, or define H in the text as the split version. The split version is also a minor of G_phi (skip the x-xbar contraction). Mention in the table caption that both versions were checked.

### 13. [minor] "The smallest instance" needs a one-line justification
- H is bipartite (clauses vs. variables plus p), so a non-planar H needs a K_{3,3} minor, hence at least 6 vertices: m + n >= 5.
- m = 2 or n = 1 gives a planar H. So (m, n) = (3, 2) with every clause containing both variables is forced, and phi* is minimal up to renaming and signs.
- **Fix:** add "(smallest in m + n; with fewer clauses or variables H(phi) is planar)".

### 14. [minor] Pre-empt the symmetric objection to the paper's own Delta = 4 reduction
A referee will ask whether Theorem hard4's G_phi has the same problem. `theorem_delta3.md`/`theorem_D.md` note that it has no global reference path. Add one sentence to S1.6 or the main text: "The reduction of Theorem [hard4] has no reference structure adjacent to all clause gadgets, so Lemma S2 does not apply to it."

### 15. [wording] Tone and precision
- S1 intro: "where the proof of that statement does not carry over to the planar case". The proof *is* for the planar case. Better: "and the single step of its proof (planarity of the constructed graph) that we believe does not hold in general."
- S1.2 item 1 / S1.7 item 1 / main.tex l. 478: "are correct" and "the reduction given there is correct as a reduction to general subcubic graphs". We reproduced a reconstruction. Better: "we reproduced the equivalence and have no objection to it; it establishes NP-completeness on general subcubic graphs."
- S1.6 title, "Why the planarity step cannot be repaired within this design": acceptable only with the hypotheses of issue 6. Suggest "Why a different drawing or planar 3-SAT variant does not help".
- S1.6: "the planarity claim fails" should be "G_phi is not planar for these formulas".
- S1.5: "together with the instance from our earlier analysis" refers to nothing the reader has seen. Better: "and a larger random instance phi0".
- Numbering clash: Note S1, Theorem S1, eq. (S1), Table S1 and Fig. S1 all appear. Label the equation or drop its number.
- Fig. S1(a) labels the shaded boxes "C1, C2, C3" while the caption calls them Q1, Q2, Q3. Use Q_j in panel (a) and C_j in panel (b).
- Courtesy: consider contacting the corresponding author of Johnson et al. before submission. A referee may ask, and it reads as good faith.

### 16. S1.1 (checked; no issue beyond issue 1)
- The yes-iff condition "no K_4 component and fr(G) <= k" is correct:
  - Theorem S1 covers disconnected graphs and vertices of degree <= 2, since it is stated for maximum degree <= 3.
  - A K_4 component admits no IOCT: removing <= 1 vertex leaves a triangle.
  - fr is additive over components and is computable on planar graphs without a given embedding.
- Optionally add that the *existence* version ("does G have any IOCT?") is equally trivial: yes iff there is no K_4 component. This covers any problem definition in the original.
- The P = NP implication correctly rests on Theorem S1. It might say "(given the proof of Theorem S1 in the main text)".

---

## D. What I could not check
- All items tied to the original text (F1-F3, the planar 3-SAT variant used, the two-literal clause handling, Fig. 3).
- The exact definition in KK93 and Pilz 2019.

These remain \verify items, extended as described in issues 2, 5 and 9.

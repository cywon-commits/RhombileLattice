# Code and data for "Coordination three heals, coordination four is hard"

This folder holds the code, raw logs and manuscript for a study of the ground-state
complexity of the antiferromagnetic three-state Potts model with a costly third state,

    E(s) = #{bonds with equal states} + D3 * #{sites in state 3},   s : V -> {1, 2, 3},

on planar lattices, as a function of the maximum coordination Δ and of D3.

- Manuscript: `paper/main.tex` (REVTeX 4.2, Phys. Rev. E format); see `paper/README.md`.
- Proof notes: `theorem_A.md` (threshold, local structure), `theorem_C.md` (coordination three,
  healing theorem), `theorem_D.md` (Δ = 4 hardness), `lit_review.md` (literature checks).
- Referee passes: `reviews/`.
- Working summary (Korean): `HANDOFF.md`.

## Setup

```
pip install -r requirements.txt
```

Every script runs from this folder (`d3_complexity/`) unless noted; scripts in `reviews/` import
from this folder, so run them as `PYTHONPATH=. python3 reviews/<script>.py ...`. All exact optima use
`scipy.optimize.milp` (HiGHS); no commercial solver is needed. To build the paper you need
TeX Live with REVTeX (`texlive-publishers` on Debian/Ubuntu); then `cd paper && latexmk -pdf main.tex`.

## Where each result in the paper comes from

| Paper item | Script(s) | Command | Stored output |
|---|---|---|---|
| Lemma (local recolouring), Theorem (threshold), Proposition (wheels) | `verify_theorem_A.py` | `python3 verify_theorem_A.py` | printed |
| Exact optimum and envelope OPT(D3) | `potts_exact.py` | library: `solve(n, edges, D3)`, `curve(...)` | — |
| Healing theorem, exhaustive check (all connected Δ ≤ 3 graphs, n ≤ 12) | `exhaustive_heal.py` | `python3 exhaustive_heal.py 12` | `hub_runs/exhaustive_12.log` |
| Healing theorem, intermediate lemmas (n ≤ 9) | `reviews/tight_exhaustive.py` | `PYTHONPATH=. python3 reviews/tight_exhaustive.py 9` | `reviews/c5b_review.md` |
| Healing theorem, MILP stress test (14–60 sites) | `reviews/milp_stress.py` | `PYTHONPATH=. python3 reviews/milp_stress.py SEED N1 N2` | `reviews/c5b_review.md` |
| Healing theorem, second referee pass (lemmas n ≤ 11; chain test; MILP 16–100 sites) | `reviews/c5_review2_*.py` | `..._lemmas.py 11`, `..._structural.py SEED NGRAPHS`, `..._milp.py SEED SECONDS` | `reviews/c5b_review2.md` |
| Greedy augmentation never gets stuck (12–20 sites) | `greedy_heal.py` | `python3 greedy_heal.py SEED COUNT NMIN NMAX` | `hub_runs/greedy_12_20.log` |
| Stronger version ("every max cut heals") is false | `heal_core.py`, `heal_moves.py` | — | `hub_runs/strong_version_counterexample.json` |
| Certification algorithm (planar T-join + 2-SAT) | `ising_tjoin.py`, `heal_check.py` | library | — |
| Timing table `tab:timing`, Fig. `fig:timing` | `timing_cert.py` | `python3 timing_cert.py > hub_runs/timing_cert.jsonl` | `hub_runs/timing_cert.jsonl` |
| Coordination-three lattices (truncated Penrose, Voronoi foam, defect honeycomb) | `lattices3.py` | library | — |
| Theorem D (Δ = 4, 0 < D3 < 2): satisfiable instances reach 6m·D3 | `cnr_reduction.py` | `python3 cnr_reduction.py SEED COUNT` | printed |
| Theorem D: unsatisfiable instance, OPT − 6m·D3 = min(D3, 2 − D3) | `cnr_unsat_test.py` | `python3 cnr_unsat_test.py` | `penrose_runs/cnr_unsat.log` |
| Theorem D: planarity of the construction (307 planar 3-SAT instances) | `cnr_planarity_check.py` | `python3 cnr_planarity_check.py 300 1` | `reviews/cnr_planarity.md` |
| Gap in the planarity argument of Johnson et al. (2025), Thm 11 | `planarity_gap.py`, `verify_delta3_reduction.py` | `python3 planarity_gap.py` | `theorem_delta3.md` |
| Penrose face-adjacency graphs (Δ = 4): envelopes and kinks | `penrose.py`, `penrose_scan.py`, `penrose_boundary.py` | `python3 penrose_scan.py RADIUS [SEED]` | `penrose_runs/SUMMARY.md` |
| All figures (`paper/figs/*.pdf`) and the generated timing table | `paper/make_figs.py` | `cd paper && python3 make_figs.py` | `paper/figs/` |
| Explainer web page | `report/gen_svgs.py`, `report/build.py` | `cd report && python3 gen_svgs.py && python3 build.py` | `report/d3_explainer.html` |

Running times: the exhaustive check at n = 12 and the largest timing runs (about 8,000 sites)
take minutes to hours on a laptop; everything else finishes in seconds to minutes.

## Exploratory scripts (not used in the paper)

`analyze_gjs_or.py`, `crossover_search.py`, `gadget.py`, `test_gadget.py`
(gadget search for the withdrawn Δ = 3 hardness attempt, see `theorem_delta3.md`);
`cluster_dmatroid.py` (`perfect_colouring_obstruction.md`); `gap_hunt.py`, `hub_test.py`,
`hub_tjoin.py` (hub relaxation); `heal_search.py`, `stuck_test.py`, `penrose_inspect.py`.

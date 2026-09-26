# Part IX handoff: annealed-bond thermodynamics (candidate standalone paper)

Working title: *Thermally created frustration on a rhombile lattice:
residual tiling entropy and a Z3 crystallisation transition.*

This file is the self-contained starting point for a session that works only
on Part IX. The rest of the repository (Parts I–VIII: quenched defects,
matching principle, exact D3=∞ solver, flux sectors) is **out of scope** for
this paper; do not pull it in.

Shared write-up: https://claude.ai/artifact/9jxHSSeFMssZzDUT7Ph4eS, §11
"Part IX" (v7.5). Read it with the Artifact tool's `read` action before
republishing, and update that same URL.

## Model

A configuration is a matching of the triangle-adjacency graph (a honeycomb:
bipartite, 3-regular) plus a 3-state Potts spin per site of the triangular
lattice.
- A matched pair of triangles shares one inactive (J=0) bond, which forms a rhombus.
- An unmatched triangle (a monomer) has all 3 bonds active. It is a frustrated triangle.
- The number of frustrated triangles is **not** constrained.

```
H = sum_{active bonds} delta(s_i, s_j) + D2 * n_2 + D3 * n_3 + mu * n_mon
```

Standard setting: D1 = 0, D3 = 10 (Ising-like), mu = 0. We compare D2 = 0
(D1 = D2) with D2 > 0 (0 = D1 < D2). The system is an L×L torus with
N_site = 3L², N_tri = 6L², and N_bond = 9L².

## Solid results (keep)

1. **T=0, D1=D2 (exact).** Every lozenge tiling in a compatible flux sector admits a
   zero-energy 2-colouring, so the ground state is the whole random-tiling
   manifold: S(0)/N_tri = 0.16153 (honeycomb dimer entropy).
2. **T=0, 0=D1<D2 (exact).** The active graph is bipartite with Σ_A deg = Σ_B deg = 2N
   and deg ≤ 6. Hence the minority (state-2) sublattice has at least N/3 sites, with
   equality iff every minority site has all 6 bonds active. That is exactly the
   rhombile tiling (3 hub-sublattice choices). So the ground state is unique,
   S(0) = 0 and E/N_tri = D2/6, valid for 0 < D2 < 6.
3. **T=∞ (exact/MC).** The frustrated fraction is 0.399. ln M/N_tri = 0.582
   (all matchings, via fugacity integration), so s(∞) = ½ln3 + 0.582 = 1.131.
4. **Thermodynamic integration closes.** For D2=0, S(T=0.06) = 0.150 against the exact
   0.1615 (combined integration and sampling error). For D2>0, S → −0.01 ≈ 0, and the
   extra 0.16/triangle of tiling entropy is released (∫C/T = 1.08 vs 0.91).
5. **D1=D2 heat capacity.**
   - The C peak (T≈0.55–0.6, C/N_tri≈0.62–0.67) does not grow for L=8–24, so it is
     a crossover (sequential cooling, single seed; redo with PT and several seeds).
   - The fixed lattice's Ising transition (T_c=1.2027 for the pristine dice lattice) is
     washed out, because about 31% of triangles are frustrated by T≈1.
6. **D2>0 has genuine Z3 long-range order.**
   - Order parameter: psi = Σ_k m_k ω^k, where m_k is the hub-type (all-6-bonds-active)
     fraction on triangular sublattice k.
   - In the ordered phase <|psi|²>N ∝ N, and the Binder cumulant
     U4 = 1 − <|psi|⁴>/(2<|psi|²>²) → 1/2.
   - T_c(D2=1) ≈ 0.52–0.53 from equilibrated single-T runs at L=48/64. D2=0.5 orders
     at a lower T.
7. **C decomposition.**
   - About half of C_mu is the frustrated-number fluctuation term, with corr(E,n) = 0.71
     at every size.
   - The canonical heat capacity is C_n = C_mu − β²cov(E,n)²/var(n).
   - mu>0 raises T_c (≈0.625 at mu=0.5, ≈0.77 at mu=1), cuts the frustrated fraction
     at T_c (0.17 → 0.08 → 0.05), and moves the frustrated-triangle C bump above the
     ordering.
8. **The psi-plane triangle is geometric, not physical.** The distribution of psi fills a
   triangle because 0 ≤ m_k ≤ 1. This is not a coexistence signal.

## Withdrawn (do not reuse)

All parallel-tempering data for L ≥ 48 were **not equilibrated**: psi² drifts
between the two halves of every run, and tau is 2,500–5,000 sweeps. This
invalidates the following:
- T_c = 0.518
- the "Potts-like" exponents measured at T_c (γ/ν = 1.72, 2β/ν = 0.28)
- the C_max "saturation" (for L ≤ 32, C_max = 0.893, 0.925, 0.965, 0.995 grows steadily)
- the Fisher-renormalisation discussion built on that saturation

Equilibration is fine for L ≤ 24 and marginal at L = 32.

## Hysteresis test (latest)

We ran 200k sweeps at a single T from ordered and from random starts
(`annealed_hysteresis.py`, `results/fss/hyst_*`).
- **L=64:** both starts converge (T=0.515 gives psi² ≈ 0.11; T=0.525 is mostly
  ordered, and the random start needs about 60k sweeps to get there).
- **L=48:** the system switches between ordered (psi²≈0.1) and disordered (≈0.01)
  every 10⁴–10⁵ sweeps. The two phases have about equal weight near T≈0.525.
- **Energy gap:** the ordered-vs-disordered energy gap ΔE/N_tri is 0.006 at L=48 and
  0.004 at L=64, against a within-phase spread of 0.003–0.004.

Conclusion: the dynamics are slow, but no phase is trapped. The transition is
continuous or at most very weakly first order. **The universality class is open.**

## Next steps (in order)

1. **Equilibrate L ≥ 48.** Either
   - (a) run single-T or PT for about 10⁶ sweeps per temperature, starting from
     the ordered state and discarding at least 100k sweeps; or
   - (b) build a faster algorithm: cluster/loop moves for the tiling, a
     Swendsen–Wang-type spin update, or multicanonical sampling in |psi|.

   Then redo the FSS with `annealed_fss_analysis.py` (exact-level WHAM,
   seed-group jackknife) and compare against 3-state Potts (1/ν=1.2,
   γ/ν=1.733, 2β/ν=0.267, α/ν=0.4) and against weak first order.
2. **Redo D1=D2 with PT and several seeds.** Also ask whether a confinement
   transition with a weak thermal signature hides at low T.
3. **mu > 0 at L = 32–48** (equilibrated) to measure α/ν away from the
   frustrated-triangle background.
4. **Literature positioning:** interacting classical dimers (Alet et al.),
   dilute/annealed Potts (Nienhuis et al.; tricritical Potts), lozenge-tiling
   roughening, and Fisher renormalisation.

## Code

| file | role |
|---|---|
| `annealed_thermo.py` | numba MC: spin Metropolis, dimer toggle, worm slide, hexagon (cube) flip; `scan` (sequential cooling) and `matchings` (ln M) modes |
| `annealed_thermo_analysis.py` | entropy by thermodynamic integration from β=0 |
| `annealed_thermo_compare.py`, `annealed_thermo_z3_figure.py` | D2 comparison and Z3 figures |
| `annealed_fss.py` | parallel tempering; args `L D2 Tmin Tmax K n_eq n_meas thin seed [mu]` → `results/fss/fss_L*_d2*[_mu*]_s*.npz` |
| `annealed_fss_analysis.py` | exact-level multi-histogram reweighting; C_mu, C_n, U4, chi', pseudo-critical fits; env `TC`, `MU` |
| `annealed_fss_psi_hist.py` | psi distribution at the C peak |
| `annealed_hysteresis.py`, `_figure.py` | ordered vs random starts, single T |
| `worm_monomer_walk.py` (DimerState), `rhombile_lattice.py`, `monomer_dimer_matching_feasibility.py` | lattice and matching infrastructure used by the above |

Data are in `results/thermo/` (sequential scans, logs) and `results/fss/`
(PT npz and hysteresis npz; job lists in `jobs_*.txt`).

## Practical notes

- **numba** is installed with `pip install --user numba`, and the functions use
  `cache=True`. Launching several runs at once right after a code change races
  on the cache. Run one tiny job first to compile.
- **Never `pkill -f <pattern>`** when the pattern also appears in your own shell
  command, because it kills the calling shell. Kill by PID instead.
- **Containers restart often.** Launch long jobs with `nohup ... | tee log` and
  write output files per run. Commit and push results as soon as they land.
- **Rough cost per sweep:** 0.24 ms at L=16, ~1 ms at L=32, ~2.1 ms at L=48,
  ~4.5 ms at L=64, ~8.5 ms at L=96. The machine has 4 cores.
- **Energies are integer-valued** when D2 and mu are integers (or multiples of
  0.5 for mu=0.5). The reweighter bins on exact levels for that reason.

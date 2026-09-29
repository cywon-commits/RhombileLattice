# Part IX: literature positioning (working notes)

## The key identification: the model is an annealed-bond TAFM

With D3 = 10 the third Potts state is frozen out, so the spins are Ising (s = 0, 1).
Energy units: an equal-spin active bond costs 1, i.e. Ising J = 1/2 in
H = J Σ σσ − h Σ σ (σ = ±1), and the D2 cost on state 1 is a uniform field h = D2/2,
so **h/J = D2**.

- **D1 = D2 (h = 0).** In a zero-energy configuration every triangle has exactly
  one equal-spin bond, and that bond is its dimer. The two spins at the ends of the
  dimer are opposite corners of the rhombus's 4-cycle, so they are equal. The T=0
  ensemble is therefore exactly the **ground-state ensemble of the triangular Ising
  antiferromagnet (TAFM)**. The dimer is the frustrated bond: the standard
  TAFM ↔ honeycomb-dimer ↔ lozenge-tiling map.
  - Check: Wannier's entropy 0.323066/site = 0.161533/triangle, which is our S(0)/N_tri.
  - Consequence: T=0 is a critical point with Coulomb-gas (Gaussian height) exponents.
    Spin correlations go as r^(−1/2) cos(K·r) (Stephenson), and the order parameter is
    φ = Σ σ ω^sub at the √3×√3 wavevector K.
  - What we expect at T>0 is the TAFM behaviour: the correlation length is finite,
    growing like exp(const/T), with no finite-T transition. Our results must test
    this (Step 2).
- **0 = D1 < D2 (h > 0).** The T=0 state is the rhombile tiling. Hubs are the minority
  spins with all 6 bonds active; the rims carry the equal-spin (dimer) bonds. This is
  exactly the **√3×√3 ground state of the TAFM in a field**, 2 up : 1 down.
  - The existence range 0 < D2 < 6 (Part IX result 2) coincides with the TAFM range
    0 < h < 6J of the √3×√3 phase.
  - In the fixed-bond TAFM this phase melts through a **3-state Potts** transition
    (β = 1/9, ν = 5/6, η = 4/15, c = 4/5), and T_c(h) → 0 as h → 0 (a Kosterlitz–Thouless
    endpoint at T=0).
  - Our Z3 crystallisation is the same symmetry breaking. Differences from the TAFM:
    the bonds are annealed, frustrated triangles are an extra fluctuating species, and
    about half of C_mu comes from their number fluctuations. The generic expectation
    is therefore 3-state Potts, possibly with Fisher-renormalised thermal exponents
    if the constraint matters. This expectation is what Step 1 has to confirm or reject.
- What is new relative to the TAFM: at a given spin configuration, the weight is not
  e^{−βJ n_eq} but e^{−β n_eq} × Z_MD(σ). Here Z_MD is the monomer-dimer partition
  function on the honeycomb with bond activity e^{β} on equal-spin bonds and 1
  otherwise. The frustrated triangles (monomers) are thermally created
  frustration. They are confined at T=0 by strings of unsatisfied bonds (energy ∝
  separation).

## References to cite (checked)

TAFM and height / Coulomb-gas picture
- G. H. Wannier, Phys. Rev. 79, 357 (1950): TAFM ground-state entropy.
- J. Stephenson, J. Math. Phys. 11, 420 (1970): T=0 spin correlations ~ r^(−1/2).
- H. W. J. Blöte and H. J. Hilhorst, J. Phys. A 15, L631 (1982): TAFM ↔ SOS / height mapping.
- B. Nienhuis, H. J. Hilhorst and H. W. J. Blöte, J. Phys. A 17, 3559 (1984): triangular SOS models.
- H. W. J. Blöte and M. P. Nightingale, Phys. Rev. B 47, 15046 (1993): critical behaviour of the ground state, c = 1, Coulomb gas.
- H. Yin, N. Gross (?) and B. Chakraborty, Phys. Rev. E 61, 6426 (2000), cond-mat/9911129:
  zero-T field theory, and a finite-T defect Coulomb gas with scaling variable t = e^(−2/T).
  **Check the author list before citing.**

TAFM in a field (3-state Potts line)
- arXiv:2306.09046 (2023): critical line of the TAFM in a field from a C3-symmetric
  corner transfer matrix. Says the √3×√3 phase for 0<h<6J melts via 3-state Potts, and
  that near h=6J the transition is hard-hexagon-like. **Author list and journal ref to be added.**
- The Blöte–Nightingale / Qian et al. transfer-matrix work on the same line. **Exact refs to be looked up.**

Lattice deformation coupled to a TAFM (closest in spirit: annealed "bond" degrees of freedom)
- Y. Shokef, A. Souslov and T. C. Lubensky, PNAS 108, 11804 (2011): the TAFM on an
  elastic lattice; phonon entropy lifts the degeneracy (order by disorder).

Interacting dimers and monomers
- F. Alet, J. L. Jacobsen, G. Misguich, V. Pasquier, F. Mila and M. Troyer, Phys. Rev. Lett. 94,
  235702 (2005): interacting square-lattice dimers, KT transition to a columnar phase.
- The Kasteleyn transition, and its rounding by a finite monomer density. This is
  relevant because our monomer number is unconstrained.
- Phys. Rev. Research 5, 043061 (2023): Ashkin–Teller/multicritical behaviour in a
  classical monomer-dimer model.

Annealed dilution and constraints
- B. Nienhuis, A. N. Berker, E. K. Riedel and M. Schick, Phys. Rev. Lett. 43, 737 (1979):
  dilute Potts model, critical and tricritical branches.
- X. Qian, Y. Deng and H. W. J. Blöte, Phys. Rev. E 72, 056132 (2005): dilute Potts model in 2D.
- M. E. Fisher, Phys. Rev. 176, 257 (1968): renormalised exponents under an annealed constraint.
  This is relevant to C_mu vs C_n (Part IX result 7).

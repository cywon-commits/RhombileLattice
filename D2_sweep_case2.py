"""D1=0, D3=1 fixed for Case2 (5-segment detour); vary D2 near 0 (+-0.1)
and see how the local "state3 line" (found in case2_sa_check.py: 3
sites at (5,3.46),(6,3.46),(7,3.46), the short direct path between the
two defects) responds.

Physical reasoning worked out with the user before running: the
pristine lattice's active bonds are a strict bipartite graph between
hub (r1, 1/3 of sites) and rim (r2/r3, 2/3 of sites), so with
D1=D2=0 any assignment where neighbors differ costs nothing -- a
Z2-like degeneracy between "hub=state1/rim=state2" and its swap. This
degeneracy is exactly why Case 0's closed loop, and Case 2's bulk
region here, can swap between the two patterns for free. Breaking the
degeneracy (D2 != 0) makes ONE of the two uniform patterns strictly
cheaper (whichever puts the D2-costly state on the smaller sublattice):
  - D2 > 0: hub->state2 (pays D2 on the 240-site minority),
            rim->state1 (free) is the cheap uniform choice.
  - D2 < 0: rim->state2 (gets a D2<0 bonus on the 480-site majority),
            hub->state1 (free) is the cheap uniform choice.
The D2=0 SA solution happened to put "outside" (658 sites, the larger
region) in the hub=state1/rim=state2 pattern and "inside" (110 sites)
in the opposite -- an arbitrary, symmetry-broken choice, since both
cost 0 at D2=0. Prediction: D2<0 should be a free lunch (that existing
choice already matches the newly-cheap pattern on the large majority-
rim outside region, so energy should drop by roughly n_rim_outside*D2
with NO change to the local structure). D2>0 should be the interesting
case: the existing choice is now the EXPENSIVE one on the large outside
region, so the true optimum needs to re-label a big coherent region --
exactly the kind of global, correlated move that single-site heat-bath
SA struggles with. Confirmed: D2=-0.1 reproduced the local 3-site
structure exactly with E=-44.6 (matches the free-bonus prediction);
D2=+0.1 gave E=31.9. Region breakdown of the SA solution: r1-outside
0% state1 (i.e. ~100% state2, costly, 219 sites), rim-outside 99%
state1 (cheap, 439 sites), r1-inside 89% state1 (cheap, 37 sites),
rim-inside 12% state1 (i.e. ~88% state2, costly, 73 sites) -- i.e. SA
did NOT get stuck; it found a FULL relabeling of both inside and
outside (compare to the D2=0 baseline: r1-outside was ~99% state1,
rim-outside ~2%, i.e. the opposite convention). Total costly sites
~219+0.88*73~283, predicted energy 0.1*283+3(local patch)=31.3,
matching the actual 31.9 well. This IS the correct global optimum
among the two topologically-valid choices (the other -- keep the old
D2=0 labeling -- would cost 0.1*(439+37)+3=50.6, confirmed worse by
direct calculation). So single-site SA handled this global relabeling
fine here, likely because the flip doesn't need to happen atomically:
a domain of the cheaper pattern can nucleate and grow one site at a
time without a large collective energy barrier.

FOLLOW-UP QUESTION FROM THE USER (not yet tested): as D2 becomes more
negative than -0.1, will the enforced inside/outside topological
boundary itself start to "leak" -- i.e. will actual bond violations
appear along the loop's own path (not just relabeling), trading a
fixed +1 violation cost for freeing a site to also collect the |D2|
bonus? And will the same mechanism eventually erode the local defect's
own state-3 patch, since D_rim_default = min(D1, D2) = D2 for D2<0
lowers the local state-3-vs-violated threshold (originally D3 < edges
resolved + D_rim_default) as D2 becomes more negative -- at D2 around
-(2-D3)=-1 for a 2-edge/"domino" site, the local patch itself should
flip over to violated links too. Both phenomena are the same
underlying trade-off (opportunity cost of a forced non-bonus site vs.
a fixed per-bond violation cost) and should have thresholds in the
same rough range (order 1-3, matching this project's other established
thresholds: 1 for isolated/single-edge sites, 2 for bulk "domino"
sites, 3 for degree-3 branch points) -- exact values depend on local
bond connectivity and haven't been computed. Proposed next test: sweep
D2 to larger magnitude (-0.5, -1.0, -1.5, -2.0) and watch for (a) the
local state-3 count changing, (b) new violated bonds appearing along
the loop away from the local patch, (c) eventual convergence to a
single rim=state2 domain that ignores the defects entirely.
"""
import pickle

import numpy as np

from rhombile_lattice import simulated_annealing, total_energy
from case2_sa_check import build_case2

D3 = 1.0
SEED_PATH = "case2_sa_result.pkl"


def run(D2_values=(-0.1, 0.0, 0.1), seed=3, n_warm=4, n_warm_sweeps=2000,
        n_fresh=3, n_fresh_sweeps=3000):
    with open(SEED_PATH, "rb") as f:
        seed_states = pickle.load(f)["best_states"]

    lat, touched, corners, (i_left, i_right) = build_case2()
    pos, sub_of = lat.site_positions()

    results = {}
    for D2 in D2_values:
        D = (0.0, D2, D3)
        rng = np.random.default_rng(seed)
        best_E, best_states = np.inf, None
        for _ in range(n_warm):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_warm_sweeps,
                                                     T_start=1.0, T_end=1e-6,
                                                     states=seed_states, record_energy=True)
            if energies[-1] < best_E:
                best_E, best_states = energies[-1], states.copy()
        for _ in range(n_fresh):
            states, energies = simulated_annealing(lat, D, rng, n_sweeps=n_fresh_sweeps,
                                                     T_start=5.0, T_end=1e-6, record_energy=True)
            if energies[-1] < best_E:
                best_E, best_states = energies[-1], states.copy()

        n3 = int((best_states == 2).sum())
        active = [b for b in lat.bonds if b["J"] != 0.0]
        nv = sum(1 for b in active if best_states[b["i"]] == best_states[b["j"]])
        print(f"D2={D2:+.2f}: E={best_E:.3f}  n_state3={n3}  n_violated={nv}")
        for p in pos[best_states == 2]:
            print("   state3 at", np.round(p, 2))
        results[D2] = (best_E, best_states, n3, nv)
    return results


if __name__ == "__main__":
    run()

"""Adds 3 more independent seeds per separation to entropic_force_T0_sizescan.pkl
(which only had 2), to pin down the T=0 entropic force's functional form --
still open in the shared artifact ("needs more L values and seeds to fit
reliably"). More seeds matter more than more L values here: the existing
5-point d(L) curve (L=12,16,20,28,36) already spans a healthy dynamic range,
but with only 2 seeds each the seed-to-seed spread (e.g. L=12: 1.48 vs 3.53,
L=20: 3.77 vs 7.31) is comparable to the signal itself, so no functional
form can be distinguished yet. Reuses run_one_separation unchanged (same
lattice, same fixed window, same T=0 dynamics) -- only adds seeds.
"""
import pickle

import numpy as np

from entropic_force_T0_sizescan import run_one_separation, SEPARATIONS

NEW_SEED_BASES = (2000, 3000, 4004)


def main():
    with open("entropic_force_T0_sizescan.pkl", "rb") as f:
        results = pickle.load(f)

    for L in SEPARATIONS:
        for seed_base in NEW_SEED_BASES:
            r = run_one_separation(L, seed=seed_base + L)
            print(f"L={L:3d} seed={seed_base+L}  d_mean={r['d_mean']:+.3f} +- {r['d_sem']:.3f} (SEM)  "
                  f"(n3={r['n3_mean']:.1f}, n3ref={r['n3ref_mean']:.1f})", flush=True)
            results[L].append(r)
        means = [x["d_mean"] for x in results[L]]
        print(f"  -> L={L}: {len(means)}-seed avg d_mean = {np.mean(means):.3f} "
              f"+- {np.std(means, ddof=1)/np.sqrt(len(means)):.3f} (seed-to-seed SEM)", flush=True)
        with open("entropic_force_T0_sizescan.pkl", "wb") as f:
            pickle.dump(results, f)

    print("\nFinal combined table:")
    for L in SEPARATIONS:
        means = [x["d_mean"] for x in results[L]]
        print(f"L={L:3d}  n_seeds={len(means)}  d_mean={np.mean(means):.3f} "
              f"+- {np.std(means, ddof=1)/np.sqrt(len(means)):.3f}")
    print("saved (merged) entropic_force_T0_sizescan.pkl")


if __name__ == "__main__":
    main()

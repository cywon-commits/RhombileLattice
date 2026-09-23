"""(c) Controlled-density isolated monomers via worm dynamics.

Split K random rhombi (2K adjacent monomers, density p = 2K/N_tri), then
let every monomer slide (worm moves, never into another monomer, no
immediate backtracking) with a Metropolis rule on the number of
site-sharing monomer pairs, annealed to T=0. Every intermediate state is
a valid covering by construction, so whatever density reaches zero
conflicts is a constructive witness: p_max >= that density (for both the
weaker edge-isolation and the strict site-isolation criterion -- the
latter also excludes the hub-sharing that pulled the earlier p~0.08 and
p~0.12 configurations' D_c up to ~2).
"""
import numpy as np

from rhombile_lattice import frustrated_triangles
from worm_monomer_walk import DimerState

NX, NY = 16, 16
DENSITIES = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.24]
N_SWEEPS = 300
T_START = 0.6
SEED = 11


def relax(st, rng, n_sweeps=N_SWEEPS, t_start=T_START):
    site_nb = [st.site_neighbors(t) for t in range(st.n)]
    last_pivot = {m: None for m in st.monomers}

    def conflicts(x, exclude=None):
        c = len(site_nb[x] & st.monomers)
        if exclude is not None and exclude in site_nb[x]:
            c -= 1
        return c

    n_anneal = int(0.7 * n_sweeps)

    def reseed_stuck_pairs():
        """Edge-adjacent monomer pairs trapped in a local minimum: re-pair
        them into a rhombus and split a fresh, conflict-free one elsewhere
        -- keeps the monomer count fixed while escaping the trap."""
        mons = st.monomers
        stuck = {(min(t, u), max(t, u)) for t in mons for u, _ in st.nbrs[t] if u in mons}
        for t, u in stuck:
            if t not in st.monomers or u not in st.monomers:
                continue
            st.bond_between(t, u)["J"] = 0.0
            st.partner[t], st.partner[u] = u, t
            st.monomers.discard(t)
            st.monomers.discard(u)
            last_pivot.pop(t, None)
            last_pivot.pop(u, None)
            for _ in range(200):
                x = int(rng.integers(st.n))
                y = st.partner[x]
                if y is None:
                    continue
                if (site_nb[x] | site_nb[y]) & st.monomers:
                    continue
                a, b = st.split(x)
                last_pivot[a] = last_pivot[b] = None
                break
            else:
                st.split(t)  # no clean spot found: put the pair back
                last_pivot[t] = last_pivot[u] = None

    for sweep in range(n_sweeps):
        T = t_start * max(0.0, 1 - sweep / max(1, n_anneal))
        if T == 0 and sweep % 20 == 0 and sweep < n_sweeps - 60:
            reseed_stuck_pairs()
        order = list(st.monomers)
        rng.shuffle(order)
        for t in order:
            if t not in st.monomers:
                continue
            cands = [u for u, _ in st.nbrs[t] if u not in st.monomers and u != last_pivot.get(t)]
            if not cands:
                continue
            u = cands[rng.integers(len(cands))]
            w = st.partner[u]
            delta = conflicts(w, exclude=t) - conflicts(t)
            if delta <= 0 or (T > 0 and rng.random() < np.exp(-delta / T)):
                st.move(t, u)
                last_pivot.pop(t, None)
                last_pivot[w] = u
        total = sum(conflicts(m) for m in st.monomers) // 2
        if total == 0 and T == 0:
            break
    return site_nb


def report(st, site_nb):
    mons = st.monomers
    edge_pairs = sum(1 for m in mons for u, _ in st.nbrs[m] if u in mons) // 2
    site_pairs = sum(len(site_nb[m] & mons) for m in mons) // 2
    site_isolated = sum(1 for m in mons if not (site_nb[m] & mons))
    return edge_pairs, site_pairs, site_isolated


def build(p, rng, n_sweeps=N_SWEEPS):
    st = DimerState(NX, NY)
    k = int(round(p * st.n / 2))
    pool = list(range(st.n))
    rng.shuffle(pool)
    placed = 0
    for t in pool:
        if placed == k:
            break
        if t in st.monomers:
            continue
        st.split(t)
        placed += 1
    site_nb = relax(st, rng, n_sweeps=n_sweeps)
    return st, site_nb


def main():
    import sys
    densities = [float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else DENSITIES
    n_sweeps = int(sys.argv[2]) if len(sys.argv) > 2 else N_SWEEPS
    rng = np.random.default_rng(SEED)
    for p in densities:
        st, site_nb = build(p, rng, n_sweeps)
        n_ft = len(frustrated_triangles(st.lat))
        edge_pairs, site_pairs, site_iso = report(st, site_nb)
        m = len(st.monomers)
        print(f"p={m/st.n:.4f}  monomers={m}  frustrated_check={n_ft}  "
              f"edge-adjacent pairs={edge_pairs}  site-sharing pairs={site_pairs}  "
              f"site-isolated monomers={site_iso}/{m}  (sweeps={n_sweeps})", flush=True)


if __name__ == "__main__":
    main()

"""Thermodynamics of the annealed-bond model: bonds AND spins both fluctuate.

Configuration = (matching of the triangle graph, Potts spins). A matched
pair of triangles shares one inactive (J=0) bond; a bond of the triangular
lattice is active iff it is not a dimer. Unmatched triangles are monomers
(= frustrated triangles). No constraint on the monomer number; an optional
chemical potential mu costs mu per monomer (default 0: the only cost of a
monomer is what its strings do to the spins).

  H = sum_{active b} delta(s_i, s_j) + D2 * n_2 + D3 * n_3 + mu * n_mon

D2>0 (with D1=0) breaks the state1/state2 symmetry: at T=0 the minority
sublattice (state2) must be as small as possible, which forces every
state2 site to have all 6 bonds active -> the unique rhombile tiling
(3 hub-sublattice choices), S(0)=0, for 0<D2<6.

T=0: every perfect matching (lozenge tiling) in a flux sector compatible
with a 2-colouring has E=0 at any D3>0, so the ground state is the whole
random-tiling manifold, S(0)/N_tri = honeycomb dimer entropy 0.16153...
T=inf: spins uniform over 3 states, bonds uniform over all matchings.

Moves (all symmetric proposals + Metropolis):
  spin   : site -> one of its 2 other states
  toggle : bond b between triangles a,c: dimer -> removed (2 monomers), or
           both monomers -> dimer added
  slide  : worm move, monomer t next to u (matched to w): t-u matched, w freed
  hexflip: the 6 bonds around a site alternate dimer/active -> swap the two
           halves (lozenge-tiling cube flip), optionally with a new random
           state for the centre site

Entropy by thermodynamic integration in beta from beta=0:
  ln Z(beta) = ln Z(0) - int_0^beta E dbeta',  S = ln Z + beta E,
  ln Z(0) = N_site ln 3 + ln M,  ln M = int_0^1 <n_dimer>_z / z dz (dimer
fugacity z, spins absent; run 'matchings').

Usage:
  python3 annealed_thermo.py scan <L> <D3> <mu> <n_eq> <n_meas> <seed> [Tmin Tmax nT [D2]]
  python3 annealed_thermo.py matchings <L> <n_sweeps> <seed>
"""
import sys

import numpy as np
from numba import njit

from worm_monomer_walk import DimerState


def build_arrays(L):
    st = DimerState(L, L)
    lat = st.lat
    bid, bi, bj, bta, btb, dimer = {}, [], [], [], [], []
    for ti in range(st.n):
        for u, b in st.nbrs[ti]:
            if id(b) in bid:
                continue
            bid[id(b)] = len(bi)
            bi.append(b["i"]); bj.append(b["j"]); bta.append(ti); btb.append(u)
            dimer.append(b["J"] == 0.0)
    nb, nt, ns = len(bi), st.n, lat.n_sites
    tri_b = np.full((nt, 3), -1, np.int64)
    for ti in range(nt):
        for k, (u, b) in enumerate(st.nbrs[ti]):
            tri_b[ti, k] = bid[id(b)]
    site_b = np.full((ns, 6), -1, np.int64)
    deg = np.zeros(ns, np.int64)
    for k in range(nb):
        for s in (bi[k], bj[k]):
            site_b[s, deg[s]] = k
            deg[s] += 1
    assert (deg == 6).all()
    # cyclic order of the 6 bonds around each site: consecutive bonds share
    # a triangle (the hexagon of the honeycomb around that site)
    bta_a, btb_a = np.array(bta), np.array(btb)
    hexb = np.full((ns, 6), -1, np.int64)
    for s in range(ns):
        bs = list(site_b[s])
        tri_of = {b: (bta_a[b], btb_a[b]) for b in bs}
        order = [bs[0]]
        prev_t = tri_of[bs[0]][0]
        cur_t = tri_of[bs[0]][1]
        while len(order) < 6:
            nxt = [b for b in bs if b != order[-1] and cur_t in tri_of[b]]
            assert len(nxt) == 1
            b = nxt[0]
            order.append(b)
            ta, tb = tri_of[b]
            prev_t, cur_t = cur_t, (tb if ta == cur_t else ta)
        assert cur_t == tri_of[bs[0]][0]
        hexb[s] = order
    mate = np.full(nt, -1, np.int64)
    for k in range(nb):
        if dimer[k]:
            mate[bta[k]] = k
            mate[btb[k]] = k
    return dict(bi=np.array(bi), bj=np.array(bj), bta=bta_a, btb=btb_a,
                dimer=np.array(dimer), mate=mate, tri_b=tri_b, site_b=site_b,
                hexb=hexb, ns=ns, nt=nt, nb=nb)


@njit(cache=True)
def _seed(s):
    np.random.seed(s)


@njit(cache=True)
def total_energy(s, dimer, bi, bj, D2, D3, mu, mate):
    e = 0.0
    for k in range(bi.shape[0]):
        if not dimer[k] and s[bi[k]] == s[bj[k]]:
            e += 1.0
    for i in range(s.shape[0]):
        if s[i] == 2:
            e += D3
        elif s[i] == 1:
            e += D2
    for t in range(mate.shape[0]):
        if mate[t] < 0:
            e += mu
    return e


@njit(cache=True)
def sweep(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, mu, E):
    ns, nb = s.shape[0], bi.shape[0]
    # spins
    for _ in range(ns):
        v = np.random.randint(ns)
        old = s[v]
        new = (old + 1 + np.random.randint(2)) % 3
        dE = 0.0
        for q in range(6):
            k = site_b[v, q]
            if dimer[k]:
                continue
            j = bj[k] if bi[k] == v else bi[k]
            dE += int(s[j] == new) - int(s[j] == old)
        dE += D3 * (int(new == 2) - int(old == 2)) + D2 * (int(new == 1) - int(old == 1))
        if dE <= 0 or np.random.random() < np.exp(-beta * dE):
            s[v] = new
            E += dE
    # bond toggles and slides
    for _ in range(nb):
        k = np.random.randint(nb)
        a, c = bta[k], btb[k]
        sat = 1.0 if s[bi[k]] == s[bj[k]] else 0.0
        if np.random.random() < 0.5:                       # toggle
            if dimer[k]:
                dE = sat + 2 * mu
                if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                    dimer[k] = False
                    mate[a] = -1
                    mate[c] = -1
                    E += dE
            elif mate[a] < 0 and mate[c] < 0:
                dE = -sat - 2 * mu
                if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                    dimer[k] = True
                    mate[a] = k
                    mate[c] = k
                    E += dE
        else:                                              # slide
            if np.random.random() < 0.5:
                t, u = a, c
            else:
                t, u = c, a
            if mate[t] >= 0 or mate[u] < 0:
                continue
            k2 = mate[u]
            w = btb[k2] if bta[k2] == u else bta[k2]
            sat2 = 1.0 if s[bi[k2]] == s[bj[k2]] else 0.0
            dE = sat2 - sat
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                dimer[k] = True
                dimer[k2] = False
                mate[t] = k
                mate[u] = k
                mate[w] = -1
                E += dE
    # hexagon (cube) flips
    for _ in range(ns):
        v = np.random.randint(ns)
        par = -1
        if dimer[hexb[v, 0]] and dimer[hexb[v, 2]] and dimer[hexb[v, 4]] and not (
                dimer[hexb[v, 1]] or dimer[hexb[v, 3]] or dimer[hexb[v, 5]]):
            par = 0
        elif dimer[hexb[v, 1]] and dimer[hexb[v, 3]] and dimer[hexb[v, 5]] and not (
                dimer[hexb[v, 0]] or dimer[hexb[v, 2]] or dimer[hexb[v, 4]]):
            par = 1
        if par < 0:
            continue
        old = s[v]
        new = old if np.random.random() < 0.5 else np.random.randint(3)
        dE = D3 * (int(new == 2) - int(old == 2)) + D2 * (int(new == 1) - int(old == 1))
        for q in range(6):
            k = hexb[v, q]
            j = bj[k] if bi[k] == v else bi[k]
            act_old = not dimer[k]
            act_new = not act_old
            dE += int(act_new and s[j] == new) - int(act_old and s[j] == old)
        if dE <= 0 or np.random.random() < np.exp(-beta * dE):
            for q in range(6):
                k = hexb[v, q]
                dimer[k] = not dimer[k]
                if dimer[k]:
                    mate[bta[k]] = k
                    mate[btb[k]] = k
            s[v] = new
            E += dE
    return E


@njit(cache=True)
def observe(s, mate, dimer, site_b):
    nm = 0
    for t in range(mate.shape[0]):
        if mate[t] < 0:
            nm += 1
    n3 = 0
    for i in range(s.shape[0]):
        if s[i] == 2:
            n3 += 1
    n2 = 0
    hub = 0
    for i in range(s.shape[0]):
        if s[i] == 1:
            n2 += 1
        full = True
        for q in range(6):
            if dimer[site_b[i, q]]:
                full = False
        if full:
            hub += 1
    return nm, n3, n2, hub


@njit(cache=True)
def run_T(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, mu, E, n_eq, n_meas):
    for _ in range(n_eq):
        E = sweep(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, mu, E)
    acc = np.zeros(8)
    for _ in range(n_meas):
        E = sweep(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, D3, mu, E)
        nm, n3, n2, hub = observe(s, mate, dimer, site_b)
        acc[0] += E
        acc[1] += E * E
        acc[2] += nm
        acc[3] += nm * nm
        acc[4] += n3
        acc[5] += E * nm
        acc[6] += n2
        acc[7] += hub
    return E, acc / n_meas


@njit(cache=True)
def matching_sweeps(dimer, mate, bta, btb, z, n_eq, n_meas):
    nb = bta.shape[0]
    nd = 0
    for k in range(nb):
        nd += int(dimer[k])
    tot = 0.0
    for it in range(n_eq + n_meas):
        for _ in range(nb):
            k = np.random.randint(nb)
            a, c = bta[k], btb[k]
            if dimer[k]:
                if np.random.random() < 1.0 / z:
                    dimer[k] = False
                    mate[a] = -1
                    mate[c] = -1
                    nd -= 1
            elif mate[a] < 0 and mate[c] < 0 and np.random.random() < z:
                dimer[k] = True
                mate[a] = k
                mate[c] = k
                nd += 1
        if it >= n_eq:
            tot += nd
    return tot / n_meas


def scan(L, D3, mu, n_eq, n_meas, seed, tmin=0.25, tmax=6.0, nT=48, D2=0.0):
    A = build_arrays(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    s = rng.integers(3, size=A["ns"]).astype(np.int64)
    dimer, mate = A["dimer"].copy(), A["mate"].copy()
    E = total_energy(s, dimer, A["bi"], A["bj"], D2, D3, mu, mate)
    # grid uniform in beta, plus beta=0 point handled analytically
    betas = np.linspace(1.0 / tmax, 1.0 / tmin, nT)
    print(f"# annealed L={L} N_site={A['ns']} N_tri={A['nt']} N_bond={A['nb']} D2={D2} D3={D3} mu={mu} "
          f"n_eq={n_eq} n_meas={n_meas} seed={seed}")
    print("# T  beta  E/Ntri  C/Ntri  n_mon(frac)  chi_mon  n3/Nsite  cov(E,nmon)/Ntri  n2/Nsite  hub_frac")
    for beta in betas:
        E, a = run_T(s, dimer, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"],
                     beta, D2, D3, mu, E, n_eq, n_meas)
        e_chk = total_energy(s, dimer, A["bi"], A["bj"], D2, D3, mu, mate)
        assert abs(e_chk - E) < 1e-6, (e_chk, E)
        nt = A["nt"]
        C = beta ** 2 * (a[1] - a[0] ** 2) / nt
        chi = (a[3] - a[2] ** 2) / nt
        cov = (a[5] - a[0] * a[2]) / nt
        print(f"{1 / beta:.4f} {beta:.4f} {a[0] / nt:.6f} {C:.5f} {a[2] / nt:.6f} {chi:.5f} "
              f"{a[4] / A['ns']:.6f} {cov:.5f} {a[6] / A['ns']:.6f} {a[7] / A['ns']:.6f}", flush=True)


def matchings(L, n_sweeps, seed):
    A = build_arrays(L)
    _seed(seed)
    dimer, mate = A["dimer"].copy(), A["mate"].copy()
    zs = np.concatenate([np.linspace(0.005, 0.1, 10), np.linspace(0.12, 1.0, 45)])
    nd = []
    print(f"# matchings L={L} N_tri={A['nt']} N_bond={A['nb']} sweeps={n_sweeps} seed={seed}")
    print("# z  <n_dimer>/N_tri  monomer_frac")
    for z in zs:
        m = matching_sweeps(dimer, mate, A["bta"], A["btb"], z, n_sweeps // 3, n_sweeps)
        nd.append(m / A["nt"])
        print(f"{z:.4f} {nd[-1]:.6f} {1 - 2 * nd[-1]:.6f}", flush=True)
    x = np.concatenate([[0.0], zs])
    f = np.concatenate([[A["nb"] / A["nt"]], np.array(nd) / zs])
    print(f"# lnM/N_tri = {np.trapezoid(f, x):.6f}")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "scan":
        L, D3, mu = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        n_eq, n_meas, seed = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
        extra = [float(x) for x in sys.argv[8:10]] + ([int(sys.argv[10])] if len(sys.argv) > 10 else [])
        D2 = float(sys.argv[11]) if len(sys.argv) > 11 else 0.0
        scan(L, D3, mu, n_eq, n_meas, seed, *extra, D2=D2)
    else:
        matchings(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))

"""Pure Ising (sigma in {0, 1}) Metropolis sweep for the annealed-bond model,

  H = sum_{active b} delta(sigma_i, sigma_j) + D2 * n_1 + mu * n_mon,

i.e. the model of annealed_thermo.sweep with the third Potts state removed
outright instead of frozen by D3 = 10. Same move set and the same per-sweep
move counts, so runs are directly comparable with the D3 = 10 implementation:

  spin   : single-site flip sigma -> 1 - sigma (symmetric)
  toggle : dimer removed (two monomers) / two adjacent monomers paired
  slide  : worm move, monomer t next to u (matched to w): t-u matched, w freed
  hexflip: the 6 bonds around a site alternate dimer/active -> swap the two
           halves; the centre spin is flipped together with probability 1/2
           (symmetric proposal; the D3 = 10 version instead draws a new
           3-state value with probability 1/2)

Spins are stored as int64 0/1 so the arrays, total_energy() and the order
parameters of annealed_thermo / annealed_fss / annealed_long work unchanged
(index 1 is the D2-penalised minority state, index 2 never occurs).
"""
import numpy as np
from numba import njit


@njit(cache=True)
def sweep_ising(s, dimer, mate, bi, bj, bta, btb, site_b, hexb, beta, D2, mu, E):
    ns, nb = s.shape[0], bi.shape[0]
    # spins: single-site flips
    for _ in range(ns):
        v = np.random.randint(ns)
        old = s[v]
        new = 1 - old
        dE = 0.0
        for q in range(6):
            k = site_b[v, q]
            if dimer[k]:
                continue
            j = bj[k] if bi[k] == v else bi[k]
            dE += int(s[j] == new) - int(s[j] == old)
        dE += D2 * (new - old)
        if dE <= 0 or np.random.random() < np.exp(-beta * dE):
            s[v] = new
            E += dE
    # bond toggles and slides (identical to annealed_thermo.sweep)
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
    # hexagon (cube) flips, centre spin flipped with probability 1/2
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
        new = old if np.random.random() < 0.5 else 1 - old
        dE = D2 * (new - old)
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

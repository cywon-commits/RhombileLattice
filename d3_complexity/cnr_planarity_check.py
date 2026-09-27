"""Planarity check for the Choi-Nakajima-Rim (1989) Delta=4 construction used in Theorem D.

Run from d3_complexity/:   python cnr_planarity_check.py [n_random] [seed]

G_phi is built as in cnr_reduction.build, but with explicit
  * occ_order[v]  : the cyclic order of the occurrence blocks k = 0..n_v-1 around the cycle of v
                    (a list of (clause j, position pos) pairs),
  * lit_order[j]  : the order in which clause j's literals are placed on w1->w2->w3->w4.
Signs fix which pair of tops is used (positive: tops 4k+1, 4k+2 sharing b^{k3};
negative: tops 4k, 4k+1 sharing b^{k2}), exactly as in cnr_reduction.build.

Tests
  A. Planar 3-SAT instances (each clause has exactly 3 distinct variables, incidence graph planar;
     part of them also with a Lichtenstein variable cycle).  Orders are read off a planar
     embedding of the incidence graph (rotation system) with the four conventions
     (variable cw/ccw) x (clause cw/ccw).
  B. For comparison: random orders (how often is G_phi planar without the embedding?).
  C. Small instances: brute force over all occurrence orders and literal orders.
  D. MILP energy check OPT == 6 m D3 <=> satisfiable on a few planar instances.
"""
import itertools
import math
import random
import sys

import networkx as nx
import numpy as np
from scipy.spatial import Delaunay

from potts_exact import solve


# ---------------------------------------------------------------- construction
def build_ordered(nv, clauses, occ_order=None, lit_order=None):
    """G_phi with explicit occurrence / literal orders.  Returns an nx.Graph (labelled nodes)."""
    occ = {v: [] for v in range(1, nv + 1)}
    for j, cl in enumerate(clauses):
        for pos, lit in enumerate(cl):
            occ[abs(lit)].append((j, pos))
    if occ_order is None:
        occ_order = occ
    if lit_order is None:
        lit_order = {j: [0, 1, 2] for j in range(len(clauses))}
    G = nx.Graph()
    tops = {}
    for v, L in occ_order.items():
        n = len(L)
        if n == 0:
            continue
        assert sorted(L) == sorted(occ[v])
        cyc = [("b", v, i) for i in range(4 * n)]
        for i in range(4 * n):
            G.add_edge(cyc[i], cyc[(i + 1) % (4 * n)])
            G.add_edge(("a", v, i), cyc[i])
            G.add_edge(("a", v, i), cyc[(i + 1) % (4 * n)])
        for k, (j, pos) in enumerate(L):
            positive = clauses[j][pos] > 0
            t1, t2 = (4 * k + 1, 4 * k + 2) if positive else (4 * k, 4 * k + 1)
            tops[(j, pos)] = (("a", v, t1), ("a", v, t2))
    for j in range(len(clauses)):
        w = [("w", j, i) for i in range(4)]
        G.add_edge(w[0], w[3])
        for slot, pos in enumerate(lit_order[j]):
            x, y = tops[(j, pos)]
            G.add_edge(w[slot], x)
            G.add_edge(w[slot + 1], y)
    return G


def incidence(nv, clauses, var_cycle=None):
    H = nx.Graph()
    for j, cl in enumerate(clauses):
        for lit in cl:
            H.add_edge(("c", j), ("v", abs(lit)))
    if var_cycle:
        for a, b in zip(var_cycle, var_cycle[1:] + var_cycle[:1]):
            H.add_edge(("v", a), ("v", b))
    return H


def orders_from_embedding(nv, clauses, emb, var_dir, cl_dir):
    """var_dir / cl_dir = +1 (clockwise rotation) or -1 (counter-clockwise)."""
    occ_order, lit_order = {}, {}
    for v in range(1, nv + 1):
        if ("v", v) not in emb:
            occ_order[v] = []
            continue
        rot = [u for u in emb.neighbors_cw_order(("v", v)) if u[0] == "c"]
        if var_dir < 0:
            rot = rot[::-1]
        occ_order[v] = [(u[1], [abs(l) for l in clauses[u[1]]].index(v)) for u in rot]
    for j, cl in enumerate(clauses):
        rot = [u for u in emb.neighbors_cw_order(("c", j))]
        if cl_dir < 0:
            rot = rot[::-1]
        lit_order[j] = [[abs(l) for l in cl].index(u[1]) for u in rot]
    return occ_order, lit_order


def planar_ok(G):
    return nx.check_planarity(G)[0] and max(d for _, d in G.degree()) <= 4


# ---------------------------------------------------------------- instance generators
def gen_delaunay(rng, npts):
    """Variables = random points; clauses = subset of Delaunay triangles (one clause node per face)."""
    pts = np.array([[rng.random(), rng.random()] for _ in range(npts)])
    tri = Delaunay(pts).simplices
    faces = [f for f in tri if rng.random() < 0.6] or [tri[0]]
    clauses = [[(int(x) + 1) * rng.choice((1, -1)) for x in rng.sample(list(f), 3)] for f in faces]
    used = sorted({abs(l) for c in clauses for l in c})
    ren = {v: i + 1 for i, v in enumerate(used)}
    clauses = [[ren[abs(l)] * (1 if l > 0 else -1) for l in c] for c in clauses]
    return len(used), clauses, None


def gen_rejection(rng, nv, m, lichtenstein):
    """Random exactly-3 clauses, kept if incidence graph (+ variable cycle) is planar."""
    for _ in range(1000):
        clauses = [[v * rng.choice((1, -1)) for v in rng.sample(range(1, nv + 1), 3)] for _ in range(m)]
        used = sorted({abs(l) for c in clauses for l in c})
        if len(used) < nv:
            continue
        cyc = None
        if lichtenstein:
            cyc = list(range(1, nv + 1)); rng.shuffle(cyc)
        if nx.check_planarity(incidence(nv, clauses, cyc))[0]:
            return nv, clauses, cyc
    return None


HAND = [
    # (nv, clauses)
    (3, [[1, 2, 3]]),
    (3, [[1, 2, 3], [-1, -2, -3]]),
    (4, [[1, 2, 3], [-1, 2, -4], [1, -2, -3], [-1, -2, 4]]),
    (4, [[1, 2, 3], [-1, 2, 4], [1, -3, -4], [-2, 3, 4]]),            # K4-like faces
    (4, [[1, 2, 3], [-1, -2, 4], [2, -3, -4], [1, 3, -4]]),
    (5, [[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 5, -2], [-2, -3, -4]]),  # wheel around x1
    # planar UNSAT: (x1 v x2),(x1 v -x2),(-x1 v x2),(-x1 v -x2), each split with a fresh variable
    (6, [[1, 2, 3], [1, 2, -3], [1, -2, 4], [1, -2, -4], [-1, 2, 5], [-1, 2, -5], [-1, -2, 6], [-1, -2, -6]]),
]


def sat(nv, clauses):
    return any(all(any((l > 0) == bool(a[abs(l) - 1]) for l in cl) for cl in clauses)
               for a in itertools.product((0, 1), repeat=nv))


def random_orders(rng, nv, clauses):
    occ = {v: [] for v in range(1, nv + 1)}
    for j, cl in enumerate(clauses):
        for pos, lit in enumerate(cl):
            occ[abs(lit)].append((j, pos))
    for L in occ.values():
        rng.shuffle(L)
    lo = {j: rng.sample([0, 1, 2], 3) for j in range(len(clauses))}
    return occ, lo


def brute_force(nv, clauses):
    """Fraction of all (occurrence cyclic order, literal order) choices that give a planar G_phi.
    Occurrence orders are taken up to rotation (fix first element); literal orders all 3! ."""
    occ = {v: [] for v in range(1, nv + 1)}
    for j, cl in enumerate(clauses):
        for pos, lit in enumerate(cl):
            occ[abs(lit)].append((j, pos))
    vs = [v for v in occ if occ[v]]
    per_var = [[[occ[v][0]] + list(p) for p in itertools.permutations(occ[v][1:])] for v in vs]
    per_cl = list(itertools.permutations([0, 1, 2]))
    tot = good = 0
    for oo in itertools.product(*per_var):
        occ_order = dict(zip(vs, oo))
        for lo in itertools.product(per_cl, repeat=len(clauses)):
            tot += 1
            good += planar_ok(build_ordered(nv, clauses, occ_order, dict(enumerate(lo))))
    return good, tot


# ---------------------------------------------------------------- main
def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    inst = [(nv, cl, None, "hand") for nv, cl in HAND]
    while len(inst) < len(HAND) + N:
        r = rng.random()
        if r < 0.4:
            nv, cl, cyc = gen_delaunay(rng, rng.randint(4, 30)); kind = "delaunay"
        else:
            lich = r < 0.7
            g = gen_rejection(rng, rng.randint(3, 9), rng.randint(2, 10), lich)
            if g is None:
                continue
            nv, cl, cyc = g; kind = "lichtenstein" if lich else "rejection"
        inst.append((nv, cl, cyc, kind))

    conv = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    conv_ok = {c: 0 for c in conv}
    any_emb_ok = 0
    rand_ok = rand_tot = 0
    kinds = {}
    maxdeg_bad = 0
    sizes = []
    for nv, cl, cyc, kind in inst:
        kinds[kind] = kinds.get(kind, 0) + 1
        H = incidence(nv, cl, cyc)
        ok, emb = nx.check_planarity(H)
        assert ok
        res = {}
        for c in conv:
            oo, lo = orders_from_embedding(nv, cl, emb, *c)
            G = build_ordered(nv, cl, oo, lo)
            if max(d for _, d in G.degree()) > 4:
                maxdeg_bad += 1
            assert G.number_of_nodes() == 28 * len(cl) and G.number_of_edges() == 43 * len(cl)
            res[c] = planar_ok(G)
            conv_ok[c] += res[c]
        sizes.append(len(cl))
        any_emb_ok += any(res.values())
        for _ in range(5):
            oo, lo = random_orders(rng, nv, cl)
            rand_tot += 1
            rand_ok += planar_ok(build_ordered(nv, cl, oo, lo))

    print(f"A. instances: {len(inst)}  kinds={kinds}  m in [{min(sizes)}, {max(sizes)}]")
    print(f"   max degree > 4 occurrences: {maxdeg_bad}")
    for c in conv:
        print(f"   embedding orders, variable {'cw' if c[0] > 0 else 'ccw'} / clause "
              f"{'cw' if c[1] > 0 else 'ccw'}: planar {conv_ok[c]}/{len(inst)}")
    print(f"   planar for at least one convention: {any_emb_ok}/{len(inst)}")
    print(f"B. random orders: planar {rand_ok}/{rand_tot}")

    # C. brute force on small instances
    print("C. brute force over all orders (small instances):")
    small = [(nv, cl) for nv, cl, _, _ in inst
             if len(cl) <= 4 and np.prod([1] + [max(1, math.factorial(sum(abs(l) == v for c in cl for l in c) - 1))
                                          for v in range(1, nv + 1)]) * 6 ** len(cl) <= 20000]
    seen = set()
    shown = 0
    for nv, cl in small:
        key = str(cl)
        if key in seen:
            continue
        seen.add(key)
        g, t = brute_force(nv, cl)
        print(f"   nv={nv} m={len(cl)} {cl}: planar {g}/{t}")
        shown += 1
        if shown >= 8:
            break

    # C'. smallest obstruction found for the same-sense convention
    obs = (4, [[-4, -1, 3], [1, -2, -4], [-1, -3, -4]])
    _, emb = nx.check_planarity(incidence(*obs))
    same = planar_ok(build_ordered(*obs, *orders_from_embedding(*obs, emb, 1, 1)))
    opp = planar_ok(build_ordered(*obs, *orders_from_embedding(*obs, emb, 1, -1)))
    g, t = brute_force(*obs)
    print(f"   obstruction {obs[1]}: same-sense planar={same}, opposite-sense planar={opp}, "
          f"all orders {g}/{t}")

    # D. MILP energy check
    print("D. MILP: OPT == 6 m D3  <=>  satisfiable  (on embedding-ordered planar G_phi)")
    milp_inst = [(nv, cl, cyc) for nv, cl, cyc, _ in inst if len(cl) <= 5][:10]
    milp_inst.append((HAND[-1][0], HAND[-1][1], None))  # planar unsat, m = 8
    bad = 0
    for nv, cl, cyc in milp_inst:
        H = incidence(nv, cl, cyc)
        _, emb = nx.check_planarity(H)
        oo, lo = orders_from_embedding(nv, cl, emb, 1, -1)
        G = build_ordered(nv, cl, oo, lo)
        assert planar_ok(G)
        lab = {x: i for i, x in enumerate(G.nodes())}
        E = sorted((min(lab[u], lab[v]), max(lab[u], lab[v])) for u, v in G.edges())
        s = sat(nv, cl)
        row = []
        for D3 in (0.5, 1.5):
            mono, n3, _ = solve(len(lab), E, D3)
            opt = mono + D3 * n3
            eq = abs(opt - 6 * len(cl) * D3) < 1e-7
            bad += eq != s
            row.append(f"D3={D3}: OPT={opt:g} vs {6 * len(cl) * D3:g}")
        print(f"   m={len(cl)} sat={s}  " + "; ".join(row))
    print(f"   mismatches: {bad}")


if __name__ == "__main__":
    main()

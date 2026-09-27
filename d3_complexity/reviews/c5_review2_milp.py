"""Second referee pass on Theorem C.5: targeted random search for counterexamples to ioct = fr on larger
subcubic graphs (30-80 vertices), planar and non-planar, biased toward small odd cycles and near-K4 gadgets.
Exact MILPs (scipy/HiGHS).  Usage: python3 c5_review2_milp.py SEED SECONDS
"""
import sys, time, random, collections
import numpy as np, networkx as nx
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from scipy.spatial import Delaunay

seed, BUDGET = int(sys.argv[1]), float(sys.argv[2])
rng = random.Random(seed); np.random.seed(seed)

def fr_milp(n, E):
    N = n + len(E); c = np.zeros(N); c[n:] = 1
    A = lil_matrix((2 * len(E), N)); lo = []; hi = []
    for k, (u, w) in enumerate(E):
        A[2*k, u] = 1; A[2*k, w] = 1; A[2*k, n+k] = 1; lo.append(1); hi.append(np.inf)
        A[2*k+1, u] = 1; A[2*k+1, w] = 1; A[2*k+1, n+k] = -1; lo.append(-np.inf); hi.append(1)
    cons = [LinearConstraint(A.tocsr(), lo, hi)]
    r = milp(c, constraints=cons, integrality=np.ones(N), bounds=Bounds(0, 1), options={"time_limit": 120})
    return round(r.fun) if r.status == 0 else None

def ioct_milp(n, E):
    N = 2 * n; c = np.zeros(N); c[n:] = 1
    A = lil_matrix((3 * len(E), N)); lo = []; hi = []
    for k, (u, w) in enumerate(E):
        A[3*k, n+u] = 1; A[3*k, n+w] = 1; lo.append(-np.inf); hi.append(1)
        A[3*k+1, u] = 1; A[3*k+1, w] = 1; A[3*k+1, n+u] = 1; A[3*k+1, n+w] = 1; lo.append(1); hi.append(np.inf)
        A[3*k+2, u] = 1; A[3*k+2, w] = 1; A[3*k+2, n+u] = -1; A[3*k+2, n+w] = -1; lo.append(-np.inf); hi.append(1)
    r = milp(c, constraints=[LinearConstraint(A.tocsr(), lo, hi)], integrality=np.ones(N), bounds=Bounds(0, 1),
             options={"time_limit": 120})
    return round(r.fun) if r.status == 0 else None

# ---------------- generators ----------------
def planar_cubicish(m):
    pts = np.random.rand(m, 2); tri = Delaunay(pts)
    G = nx.Graph(); G.add_nodes_from(range(len(tri.simplices)))
    for i, nb in enumerate(tri.neighbors):
        for j in nb:
            if j >= 0: G.add_edge(i, j)
    return G  # planar dual, max degree 3 (boundary faces degree < 3)

def random_cubic(n): return nx.random_regular_graph(3, n, seed=rng.randrange(10**9))

def subst_edge(G, u, v, kind):
    """replace edge uv by a gadget with two ports (planarity preserved)."""
    G.remove_edge(u, v); b = max(G.nodes) + 1
    if kind == "diamond":   # K4 minus edge {a,d}; ports a,d
        a, x, y, d = b, b+1, b+2, b+3
        G.add_edges_from([(a, x), (a, y), (x, y), (x, d), (y, d), (u, a), (d, v)])
    elif kind == "prism":   # triangular prism minus one rung; ports = rung ends
        a1, a2, a3, b1, b2, b3 = range(b, b+6)
        G.add_edges_from([(a1, a2), (a2, a3), (a1, a3), (b1, b2), (b2, b3), (b1, b3), (a2, b2), (a3, b3), (u, a1), (b1, v)])
    elif kind == "K4sub":   # K4 on {p,q,r,s} with edge pq subdivided by x and rs subdivided by z; ports x,z
        p, q, r, s, x, z = range(b, b+6)
        G.add_edges_from([(p, x), (x, q), (p, r), (p, s), (q, r), (q, s), (r, z), (z, s), (u, x), (z, v)])
    elif kind == "tri":     # subdivide by a triangle with one pendant-free vertex (degree 2)
        a, c, d = b, b+1, b+2
        G.add_edges_from([(a, c), (c, d), (a, d), (u, a), (d, v)])
    elif kind == "sub":
        G.add_edges_from([(u, b), (b, v)])

def truncate(G, v):
    nb = list(G[v]);
    if len(nb) != 3: return
    b = max(G.nodes) + 1; T = [b, b+1, b+2]
    G.remove_node(v); G.add_edges_from([(T[0], T[1]), (T[1], T[2]), (T[0], T[2])])
    for t, x in zip(T, nb): G.add_edge(t, x)

def triangle_rich(m):
    """disjoint small odd cycles (3,3,5) + near-K4s, joined by random edges respecting degree 3."""
    G = nx.Graph(); b = 0
    while b < m:
        kind = rng.choice(["C3", "C3", "C5", "diamond", "K4sub1"])
        if kind == "C3": G.add_edges_from([(b, b+1), (b+1, b+2), (b, b+2)]); b += 3
        elif kind == "C5": G.add_edges_from([(b+i, b+(i+1) % 5) for i in range(5)]); b += 5
        elif kind == "diamond": G.add_edges_from([(b, b+1), (b, b+2), (b+1, b+2), (b+1, b+3), (b+2, b+3)]); b += 4
        else:  # K4 with one subdivided edge
            G.add_edges_from([(b, b+4), (b+4, b+1), (b, b+2), (b, b+3), (b+1, b+2), (b+1, b+3), (b+2, b+3)]); b += 5
    for _ in range(4 * m):
        a, c = rng.sample(list(G.nodes), 2)
        if G.degree(a) < 3 and G.degree(c) < 3 and not G.has_edge(a, c): G.add_edge(a, c)
    return G

def gadgetize(G, frac_edges, frac_vertices, kinds):
    G = nx.convert_node_labels_to_integers(G)
    for v in rng.sample(list(G.nodes), int(frac_vertices * G.number_of_nodes())):
        if v in G and G.degree(v) == 3: truncate(G, v)
    for e in rng.sample(list(G.edges), int(frac_edges * G.number_of_edges())):
        if G.has_edge(*e): subst_edge(G, e[0], e[1], rng.choice(kinds))
    return G

def has_K4_comp(G):
    return any(len(C) == 4 and G.subgraph(C).number_of_edges() == 6 for C in nx.connected_components(G))

def sample():
    fam = rng.choice(["cubic", "cubic_gadget", "planar", "planar_gadget", "tririch", "tririch_planar_filter"])
    if fam == "cubic": G = random_cubic(rng.randrange(30, 81, 2))
    elif fam == "cubic_gadget":
        G = gadgetize(random_cubic(rng.randrange(10, 31, 2)), rng.uniform(0.1, 0.5), rng.uniform(0, 0.5),
                      ["diamond", "prism", "K4sub", "tri", "sub"])
    elif fam == "planar": G = planar_cubicish(rng.randint(20, 45))
    elif fam == "planar_gadget":
        G = gadgetize(planar_cubicish(rng.randint(10, 25)), rng.uniform(0.1, 0.5), rng.uniform(0, 0.6),
                      ["diamond", "prism", "K4sub", "tri"])
    elif fam == "tririch": G = triangle_rich(rng.randint(30, 70))
    else:
        G = triangle_rich(rng.randint(30, 60))
        while not nx.check_planarity(G)[0]:
            G.remove_edge(*rng.choice(list(G.edges)))
    G = nx.convert_node_labels_to_integers(G)
    return fam, G

t0 = time.time(); cnt = collections.Counter(); planar = 0; mism = 0; sizes = []
while time.time() - t0 < BUDGET:
    fam, G = sample()
    if G.number_of_nodes() > 100 or G.number_of_edges() == 0: continue
    assert max(dict(G.degree()).values()) <= 3
    if has_K4_comp(G): cnt["skipped_K4comp"] += 1; continue
    E = list(G.edges()); n = G.number_of_nodes()
    f = fr_milp(n, E); i = ioct_milp(n, E)
    if f is None or i is None: cnt["timeout"] += 1; continue
    cnt[fam] += 1; sizes.append(n); planar += nx.check_planarity(G)[0]
    if f != i:
        mism += 1; print("MISMATCH", fam, n, "fr", f, "ioct", i, E, flush=True)
    if i < f: print("IMPOSSIBLE ioct<fr", E)
print("families", dict(cnt), "planar", planar, "total", sum(v for k, v in cnt.items() if k not in ("timeout", "skipped_K4comp")),
      "n-range", (min(sizes), max(sizes)) if sizes else None, "mean n", np.mean(sizes) if sizes else None, "MISMATCHES", mism)

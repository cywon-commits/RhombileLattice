"""Brute-force checks of Theorem A (d3_complexity/theorem_A.md). Sanity checks, not a proof.

1. Prop A.6: wheels W_D (D = 3..9): every optimum uses colour 3 iff D3 < floor(D/2).
2. Octahedron: same with threshold 2.
3. Random planar graphs (Delaunay minus random edges, max degree <= 5): in EVERY optimum,
   each colour-3 vertex satisfies min(n1,n2) >= a + D3 (A.2c), the A.5 consequences hold,
   and OPT = fr(G) whenever D3 >= floor(Delta/2) (A.3).
"""
import itertools, random
import numpy as np
from scipy.spatial import Delaunay

def all_optima(n, E, D3):
    best, opts = None, []
    for s in itertools.product((1, 2, 3), repeat=n):
        if s[0] == 2:           # 1<->2 symmetry
            continue
        c = sum(s[u] == s[v] for u, v in E) + D3 * s.count(3)
        if best is None or c < best - 1e-9:
            best, opts = c, [s]
        elif abs(c - best) < 1e-9:
            opts.append(s)
    return best, opts

def fr(n, E):
    return min(sum(s[u] == s[v] for u, v in E) for s in itertools.product((1, 2), repeat=n))

def wheel(D):
    return D + 1, [(0, i) for i in range(1, D + 1)] + [(i, i % D + 1) for i in range(1, D + 1)]

def check_threshold(name, n, E, k):
    bad = 0
    for D3 in (0.1, 0.5, k - 1 + 0.5, k - 0.01, k + 0.01, k + 0.5):
        if D3 < 0:
            continue
        _, opts = all_optima(n, E, D3)
        allc3 = all(3 in s for s in opts); anyc3 = any(3 in s for s in opts)
        ok = allc3 if D3 < k else (not anyc3)
        bad += not ok
    print(f"{name:12s} threshold {k}: {'ok' if not bad else 'FAIL'}")
    return bad

bad = 0
for D in range(3, 10):
    n, E = wheel(D); bad += check_threshold(f"W_{D}", n, E, D // 2)
octa = [(u, v) for u, v in itertools.combinations(range(6), 2) if abs(u - v) != 3]
bad += check_threshold("octahedron", 6, octa, 2)

random.seed(5); np.random.seed(5); cnt = viol = 0
while cnt < 200:
    n = random.randint(6, 10)
    pts = np.random.rand(n, 2)
    E = set()
    for tri in Delaunay(pts).simplices:
        for a, b in itertools.combinations(sorted(map(int, tri)), 2):
            E.add((a, b))
    E = [e for e in E if random.random() > 0.25]
    deg = [0] * n
    for u, v in E:
        deg[u] += 1; deg[v] += 1
    Dm = max(deg)
    if Dm > 5 or Dm < 3:
        continue
    D3 = random.choice([0.2, 0.7, 1.0, 1.3, 1.8, 2.0, 2.4, 3.1])
    best, opts = all_optima(n, E, D3)
    adj = [[] for _ in range(n)]
    for u, v in E:
        adj[u].append(v); adj[v].append(u)
    for s in opts:
        for v in range(n):
            if s[v] != 3:
                continue
            nb = [s[w] for w in adj[v]]; a = nb.count(3); n1, n2 = nb.count(1), nb.count(2)
            if min(n1, n2) < a + D3 - 1e-9: viol += 1
            if D3 > (Dm - 1) // 2 - 1 and a > 0: viol += 1          # A.5(i)
            if Dm <= 4 and 1 < D3 < 2 and not (len(nb) == 4 and n1 == n2 == 2): viol += 1  # A.5(iii)
    if D3 >= Dm // 2 and abs(best - fr(n, E)) > 1e-9: viol += 1        # A.3
    cnt += 1
print(f"random planar graphs: {cnt}, violations: {viol}")
print("ALL OK" if bad == 0 and viol == 0 else "FAILURES")

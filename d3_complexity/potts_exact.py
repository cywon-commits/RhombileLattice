"""Exact solver for  min_s  #mono edges + D3 * #colour-3  (s: V -> {1,2,3}) via MILP (HiGHS),
plus the exact piecewise-linear curve OPT(D3) on an interval (Eisner-Severance style).
"""
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


def solve(n, edges, D3, fix=None, time_limit=None):
    """Return (mono, n3, colours) of an optimal colouring. colours[v] in {1,2,3}.
    `fix` maps vertex -> set of forbidden colours (optional)."""
    ne = len(edges)
    N = 3 * n + ne
    c = np.zeros(N)
    c[[3 * v + 2 for v in range(n)]] = D3
    c[3 * n:] = 1.0
    A = lil_matrix((n + 3 * ne, N))
    lo = np.empty(n + 3 * ne)
    hi = np.empty(n + 3 * ne)
    for v in range(n):
        A[v, 3 * v:3 * v + 3] = 1
        lo[v] = hi[v] = 1
    r = n
    for k, (u, w) in enumerate(edges):
        for col in range(3):
            A[r, 3 * u + col] = 1
            A[r, 3 * w + col] = 1
            A[r, 3 * n + k] = -1
            lo[r], hi[r] = -np.inf, 1
            r += 1
    ub = np.ones(N)
    ub[3 * 0 + 1] = 0  # 1<->2 symmetry: vertex 0 never colour 2
    for v, cols in (fix or {}).items():
        for col in cols:
            ub[3 * v + col - 1] = 0
    opts = {} if time_limit is None else {"time_limit": time_limit}
    res = milp(c, constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=np.ones(N),
               bounds=Bounds(0, ub), options=opts)
    if res.status != 0:
        raise RuntimeError(res.message)
    x = np.round(res.x).astype(int)
    col = np.array([1 + int(np.argmax(x[3 * v:3 * v + 3])) for v in range(n)])
    mono = sum(col[u] == col[w] for u, w in edges)
    n3 = int((col == 3).sum())
    return int(mono), n3, col


def curve(n, edges, a, b, eps=1e-9):
    """Exact lower envelope OPT(D3) = min_s (mono_s + D3*n3_s) on [a, b].
    Returns sorted list of (D3_start, mono, n3, colours) segments."""
    cache = {}

    def at(d):
        if d not in cache:
            cache[d] = solve(n, edges, d)
        return cache[d]

    lines = {}

    def rec(lo, hi):
        ml, nl, cl = at(lo)
        mh, nh, ch = at(hi)
        lines[(ml, nl)] = cl
        lines[(mh, nh)] = ch
        if (ml, nl) == (mh, nh) or nl == nh:
            return
        d = (mh - ml) / (nl - nh)
        md, nd, cd = at(d)
        if md + d * nd >= ml + d * nl - eps:
            lines[(md, nd)] = cd
            return
        rec(lo, d)
        rec(d, hi)

    rec(a, b)
    # lower envelope of collected lines on [a,b]
    segs = []
    d = a
    while d < b - eps:
        best = min(lines, key=lambda L: (L[0] + d * L[1], L[1]))  # at d, prefer fewer colour-3 to the right
        # next breakpoint: smallest D > d where another line with fewer n3 becomes better
        nxt = b
        for L in lines:
            if L[1] < best[1]:
                x = (L[0] - best[0]) / (best[1] - L[1])
                if d + eps < x < nxt:
                    nxt = x
        segs.append((d, best[0], best[1], lines[best]))
        d = nxt
    return segs

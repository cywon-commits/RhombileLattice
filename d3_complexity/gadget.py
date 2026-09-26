"""Gadget verifier for  cost(s) = #mono edges + D3 * #colour-3 vertices.

A gadget is a small graph with an ordered list of terminal vertices.  For every colouring b
of the terminals it computes, exactly, the minimum internal cost as a function of D3:

    f_b(D3) = min over colourings of the internal vertices of (mono + D3 * n3),

where mono counts every gadget edge (internal-internal, terminal-internal and
terminal-terminal) and n3 counts colour-3 INTERNAL vertices (terminal D3 costs belong to
whoever owns the terminal; pass count_terminal_d3=True to include them).

f_b is the lower envelope of the lines mono + D3*n3, so it is stored exactly as its
Pareto frontier {n3: min mono}.  Enumeration is exhaustive (3^k internal colourings,
vectorised with numpy); keep k <= ~13.
"""
import itertools
from fractions import Fraction

import numpy as np

COLOURS = (1, 2, 3)


class Gadget:
    def __init__(self, n, edges, terminals, names=None, count_terminal_d3=False):
        self.n = n
        self.edges = [tuple(e) for e in edges]
        self.terminals = list(terminals)
        self.internal = [v for v in range(n) if v not in set(self.terminals)]
        self.names = names or {v: str(v) for v in range(n)}
        self.count_terminal_d3 = count_terminal_d3
        k = len(self.internal)
        if k > 14:
            raise ValueError(f"{k} internal vertices: too many for exhaustive enumeration")
        pos = {v: i for i, v in enumerate(self.internal)}
        tpos = {v: i for i, v in enumerate(self.terminals)}
        # all internal colourings, shape (3^k, k)
        self.C = (np.array(list(itertools.product(COLOURS, repeat=k)), dtype=np.int8).reshape(-1, k)
                  if k else np.zeros((1, 0), dtype=np.int8))
        self.n3 = (self.C == 3).sum(axis=1).astype(np.int32)
        ii = [(pos[u], pos[v]) for u, v in self.edges if u in pos and v in pos]
        self.base = np.zeros(len(self.C), dtype=np.int32)
        for a, b in ii:
            self.base += (self.C[:, a] == self.C[:, b])
        self.ti = [(tpos[u], pos[v]) if u in tpos else (tpos[v], pos[u])
                   for u, v in self.edges if (u in tpos) != (v in tpos)]
        self.tt = [(tpos[u], tpos[v]) for u, v in self.edges if u in tpos and v in tpos]

    # ---- core ------------------------------------------------------------------------
    def _mono(self, b):
        m = self.base.copy()
        for t, i in self.ti:
            m += (self.C[:, i] == b[t])
        m += sum(b[s] == b[t] for s, t in self.tt)
        return m

    def frontier(self, b):
        """Pareto frontier of (n3, mono) for terminal colouring b: dict n3 -> min mono."""
        m = self._mono(b)
        extra = sum(1 for c in b if c == 3) if self.count_terminal_d3 else 0
        out = {}
        for n3 in np.unique(self.n3):
            out[int(n3) + extra] = int(m[self.n3 == n3].min())
        # keep only lines that can be optimal for some D3 >= 0
        pts = sorted(out.items())
        keep, best = {}, None
        for n3, mono in pts:                  # increasing n3: must strictly decrease mono
            if best is None or mono < best:
                keep[n3] = mono
                best = mono
        return keep

    def table(self):
        """frontier for every terminal colouring (up to the global 1<->2 swap)."""
        res = {}
        for b in itertools.product(COLOURS, repeat=len(self.terminals)):
            sw = tuple({1: 2, 2: 1, 3: 3}[c] for c in b)
            if sw in res:
                continue
            res[b] = self.frontier(b)
        return res

    def optima(self, b, D3):
        """All optimal internal colourings for terminal colouring b at a given D3 (rows of C)."""
        m = self._mono(b) + D3 * self.n3
        best = m.min()
        return self.C[np.isclose(m, best)], float(best)

    def choice_report(self, b, D3):
        """Which internal vertices are colour 3 in all / some / no optimal colourings."""
        opt, val = self.optima(b, D3)
        always = [self.names[v] for j, v in enumerate(self.internal) if (opt[:, j] == 3).all()]
        some = [self.names[v] for j, v in enumerate(self.internal)
                if (opt[:, j] == 3).any() and not (opt[:, j] == 3).all()]
        n3_values = sorted(set(int(x) for x in (opt == 3).sum(axis=1)))
        return dict(value=val, n_optima=len(opt), always3=always, sometimes3=some, n3=n3_values)


# ---- piecewise-linear helpers (exact, on an interval) ---------------------------------
def value(front, D3):
    return min(m + D3 * k for k, m in front.items())


def breakpoints(front, lo, hi):
    """D3 values in (lo, hi) where the envelope of `front` changes slope."""
    lines = sorted(front.items())
    xs = set()
    for (k1, m1), (k2, m2) in itertools.combinations(lines, 2):
        if k1 != k2:
            x = Fraction(m1 - m2, k2 - k1)
            if lo < x < hi and abs(value(front, float(x)) - (m1 + float(x) * k1)) < 1e-9:
                xs.add(x)
    return sorted(xs)


def min_gap(front_a, front_b, lo, hi, open_ends=True):
    """min over D3 in [lo,hi] of f_a(D3) - f_b(D3); both piecewise linear, so the minimum is at
    an endpoint or a breakpoint of either envelope.  Returns (gap, argmin).  With open_ends the
    endpoints are approached from inside (useful for open intervals such as (0,1))."""
    eps = 1e-9 if open_ends else 0.0
    xs = [lo + eps, hi - eps] + [float(x) for x in breakpoints(front_a, lo, hi) + breakpoints(front_b, lo, hi)]
    return min((value(front_a, x) - value(front_b, x), x) for x in xs)


def describe(front):
    return " / ".join(f"{m}+{k}·D3" if k else f"{m}" for k, m in sorted(front.items()))

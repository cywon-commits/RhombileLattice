"""Unsatisfiable (and one satisfiable) exactly-3-literal instance for the CNR reduction test."""
import itertools, sys, time
from cnr_reduction import build, sat
from potts_exact import solve
allc = [[s1 * 1, s2 * 2, s3 * 3] for s1, s2, s3 in itertools.product((1, -1), repeat=3)]
cases = [("unsat: all 8 clauses on 3 vars", allc), ("sat: 7 clauses", allc[:7])]
for name, cl in cases:
    n, E, idx = build(3, cl); m = len(cl)
    for D3 in (0.5, 1.5, 1.9):
        t = time.time(); mono, n3, col = solve(n, E, D3)
        val = mono + D3 * n3
        print(name, "n", n, "D3", D3, "OPT", round(val, 4), "6mD3", round(6 * m * D3, 4),
              "OPT<=6mD3:", val <= 6 * m * D3 + 1e-7, "sat:", sat(3, cl), f"({time.time()-t:.0f}s)", flush=True)

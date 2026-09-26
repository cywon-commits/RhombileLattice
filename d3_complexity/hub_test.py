"""Is the hub relaxation exact on degree-3 planar lattices?

true  : min #mono + D3 * #colour-3                      (the real problem)
relax : same, but colour-3/colour-3 edges cost nothing   (the "hub" model: planar T-join with
        a hub of cost D3 inside every dual triangle; solvable by matching)
relax <= true always.  If the relaxed optimum has no adjacent colour-3 pair it is also a true
optimum, and the problem is solved by the polynomial hub model on that instance.
"""
import json
import sys
import time

from lattices3 import random_cubic, defect_honeycomb, truncated_penrose
from potts_exact import solve

FAMILIES = {
    "random": lambda s: random_cubic(60, s),
    "defect": lambda s: defect_honeycomb(10, 0.08, s),
    "Tpenrose": lambda s: truncated_penrose(4, s),
}


def run(fam, seeds, d3s):
    rows = []
    for s in seeds:
        n, E, pos, info = FAMILIES[fam](s)
        for D3 in d3s:
            t = time.time()
            m1, k1, c1 = solve(n, E, D3)
            m2, k2, c2 = solve(n, E, D3, count33=False)
            adj33 = sum(1 for u, v in E if c2[u] == 3 and c2[v] == 3)
            true_of_relax = m2 + adj33 + D3 * k2
            row = dict(fam=fam, seed=s, n=n, odd=info.get("odd_faces"), D3=D3,
                       true=m1 + D3 * k1, true_mono=m1, true_n3=k1,
                       relax=m2 + D3 * k2, relax_n3=k2, relax_adj33=adj33,
                       true_of_relax=true_of_relax, sec=round(time.time() - t, 1))
            row["gap"] = round(row["true"] - row["relax"], 6)
            rows.append(row)
            print(json.dumps(row), flush=True)
    return rows


if __name__ == "__main__":
    fam = sys.argv[1]
    seeds = range(int(sys.argv[2]))
    d3s = [0.1, 0.3, 0.5, 0.7, 0.9]
    run(fam, seeds, d3s)

"""Timing of the certification algorithm (matching + 2-SAT) on coordination-three lattices.
For each lattice: Ising ground state by planar T-join, then 2-SAT healing; if healed, the
state has energy D3*fr and is certified optimal (sandwich bound).  Prints one JSON per line."""
import json, sys, time
from lattices3 import truncated_penrose, random_cubic, defect_honeycomb
from ising_tjoin import ising_ground_state_tjoin
from heal_check import two_sat_heal

def run(fam, arg, build):
    n, edges, pos, info = build()
    t0 = time.time(); fr, col = ising_ground_state_tjoin(n, edges); t1 = time.time()
    chosen, M = two_sat_heal(n, edges, col); t2 = time.time()
    ok = chosen is not None
    if ok:
        s = list(col)
        for v in chosen: s[v] = 3
        ok = all(s[u] != s[v] for u, v in edges)
    print(json.dumps(dict(fam=fam, arg=arg, n=n, fr=fr, healed=ok,
                          t_match=round(t1 - t0, 2), t_2sat=round(t2 - t1, 3))), flush=True)

if __name__ == "__main__":
    for r in (4, 6, 8, 10, 12, 14):
        run("truncated_penrose", r, lambda r=r: truncated_penrose(r, seed=0))
    for N in (250, 500, 1000, 2000, 4000):
        run("voronoi_foam", N, lambda N=N: random_cubic(N, seed=1))
    for L in (10, 20, 30, 45):
        run("defect_honeycomb", L, lambda L=L: defect_honeycomb(L, 0.05, seed=2))

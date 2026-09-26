"""Search small degree-3 planar lattices for instances where the hub relaxation is NOT exact."""
import json, sys
import numpy as np
from lattices3 import delaunay_dual
from potts_exact import solve

def hunt(seed0, count, Nmin, Nmax, d3s=(0.2, 0.5, 0.8, 0.95)):
    rng = np.random.default_rng(seed0)
    found = []
    for it in range(count):
        N = int(rng.integers(Nmin, Nmax + 1))
        pts = rng.random((N, 2))
        n, E, pos, info = delaunay_dual(pts)
        for D3 in d3s:
            m1, k1, c1 = solve(n, E, D3)
            m2, k2, c2 = solve(n, E, D3, count33=False)
            gap = (m1 + D3 * k1) - (m2 + D3 * k2)
            if gap > 1e-6:
                rec = dict(it=it, N=N, D3=D3, gap=gap, true=(m1, k1), relax=(m2, k2),
                           pts=pts.tolist(), relax_col=[int(x) for x in c2], true_col=[int(x) for x in c1])
                found.append(rec)
                print(json.dumps({k: rec[k] for k in ("it", "N", "D3", "gap", "true", "relax")}), flush=True)
    return found

if __name__ == "__main__":
    s, c = int(sys.argv[1]), int(sys.argv[2])
    f = hunt(s, c, int(sys.argv[3]), int(sys.argv[4]))
    json.dump(f, open(f"hub_runs/gaps_{s}.json", "w"))
    print("done", len(f), "gap cases")

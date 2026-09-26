"""Ground states of the Penrose-dual Potts problem for the explainer figures (radius-6 patch)."""
import json, sys
sys.path.insert(0, "..")
import numpy as np
from penrose import tiling_graphs, E as EV
from potts_exact import solve

R, vdeg, eo, dual = tiling_graphs(6)
n = len(R)
edges = sorted({(min(a, b), max(a, b)) for a in range(n) for b in dual[a]})
out = {"polys": [], "states": {}}
for cs, jk, c in R:
    out["polys"].append([list(map(float, np.array(k) @ EV)) for k in cs])
# shared edge endpoints for each dual edge
shared = {}
for e, owners in eo.items():
    if len(owners) == 2:
        a, b = sorted(owners)
        shared[(a, b)] = [list(map(float, np.array(v) @ EV)) for v in e]
out["edges"] = [[a, b, shared[(a, b)]] for a, b in edges]
out["deg"] = [len(d) for d in dual]
for D3 in (0.5, 1.5, 2.5):
    mono, n3, col = solve(n, edges, D3)
    out["states"][str(D3)] = {"mono": mono, "n3": n3, "col": [int(x) for x in col]}
    print(D3, mono, n3, flush=True)
json.dump(out, open("penrose_states.json", "w"))

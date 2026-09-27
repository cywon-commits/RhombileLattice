"""Are the small kinks of the Penrose-dual energy curve (D3 = 1, 1.25) boundary effects?

For several patch radii and phason offsets (gamma seeds) compute the exact envelope on
[D3_lo, D3_hi], and at every breakpoint locate the rhombi whose state-3 status changes.
Report their graph distance to the patch edge (rhombi with fewer than 4 neighbours).
A boundary effect shows up as changes confined to small distance; a bulk effect as changes
spread through the interior with a count that grows with the area.
"""
import json
import sys
import time
from collections import deque

from penrose import tiling_graphs, default_gamma
from potts_exact import curve


def distance_to_edge(n, dual):
    edge = [v for v in range(n) if len(dual[v]) < 4]
    dist = [None] * n
    q = deque(edge)
    for v in edge:
        dist[v] = 0
    while q:
        u = q.popleft()
        for w in dual[u]:
            if dist[w] is None:
                dist[w] = dist[u] + 1
                q.append(w)
    return dist


def run(radius, seed, lo=0.9, hi=2.1):
    R, vdeg, eo, dual = tiling_graphs(radius, default_gamma(seed))
    n = len(R)
    edges = sorted({(min(a, b), max(a, b)) for a in range(n) for b in dual[a]})
    t = time.time()
    segs = curve(n, edges, lo, hi)
    dist = distance_to_edge(n, dual)
    out = dict(radius=radius, seed=seed, n=n, edge_sites=sum(1 for d in dist if d == 0),
               max_depth=max(dist), sec=round(time.time() - t, 1), segments=[], changes=[])
    for d, mono, n3, col in segs:
        out["segments"].append(dict(D3=round(float(d), 4), mono=mono, n3=n3))
    for (d1, m1, k1, c1), (d2, m2, k2, c2) in zip(segs, segs[1:]):
        changed = [v for v in range(n) if (c1[v] == 3) != (c2[v] == 3)]
        out["changes"].append(dict(at=round(float(d2), 4), n_changed=len(changed),
                                   depths=sorted(dist[v] for v in changed)))
    return out


if __name__ == "__main__":
    radius, seed = float(sys.argv[1]), int(sys.argv[2])
    res = run(radius, seed)
    print(json.dumps(res), flush=True)

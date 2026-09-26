"""Scan OPT(D3) on the face-adjacency (dual) graph of Penrose P3 patches."""
import sys, time, collections
from penrose import tiling_graphs, is_thick
from potts_exact import curve

def main(radius, seed=0):
    from penrose import default_gamma
    R, vdeg, eo, dual = tiling_graphs(radius, default_gamma(seed))
    n = len(R)
    edges = sorted({(min(a, b), max(a, b)) for a in range(n) for b in dual[a]})
    t = time.time()
    segs = curve(n, edges, 0.0, 3.0)
    print(f"radius={radius} seed={seed}: {n} rhombi, {len(edges)} dual edges, "
          f"max dual deg {max(len(d) for d in dual)}, {time.time()-t:.1f}s")
    for d, mono, n3, col in segs:
        print(f"  D3 >= {d:.4f}: mono={mono:4d}  n3={n3:4d}")
    return R, vdeg, eo, dual, edges, segs

if __name__ == "__main__":
    main(float(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 0)

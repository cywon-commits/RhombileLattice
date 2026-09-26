import sys, time, collections
import numpy as np
from penrose import tiling_graphs, is_thick, default_gamma
from potts_exact import solve

radius=float(sys.argv[1]); D3=float(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 0
R,vdeg,eo,dual=tiling_graphs(radius, default_gamma(seed))
n=len(R); edges=sorted({(min(a,b),max(a,b)) for a in range(n) for b in dual[a]})
t=time.time(); mono,n3,col=solve(n,edges,D3); opt=mono+D3*n3
print(f"D3={D3}: mono={mono} n3={n3} OPT={opt:.4f} ({time.time()-t:.1f}s)")
S=[v for v in range(n) if col[v]==3]
indep=all(col[w]!=3 for v in S for w in dual[v])
nbr=collections.Counter(tuple(sorted(col[w] for w in dual[v])) for v in S)
print("independent:",indep," neighbour colours of colour-3:",dict(nbr))
print("deg of colour-3:",collections.Counter(len(dual[v]) for v in S))
def corner_degs(v): return tuple(sorted(vdeg[c] for c in R[v][0]))
print("colour-3 rhombus (thick?, corner degrees):",collections.Counter((is_thick(R[v][1]),corner_degs(v)) for v in S))
deg4=[v for v in range(n) if len(dual[v])==4]
print("all deg-4 rhombi (thick?, corner degrees):",collections.Counter((is_thick(R[v][1]),corner_degs(v)) for v in deg4))
# forcedness: for each colour-3 rhombus, forbid colour 3 there and re-solve
if "--force" in sys.argv:
    forced=0; alt=0
    for v in S:
        m2,k2,_=solve(n,edges,D3,fix={v:{3}})
        if m2+D3*k2 > opt+1e-7: forced+=1
        else: alt+=1
    print(f"colour-3 rhombi forced: {forced}, replaceable: {alt}")

"""Exhaustive check on all connected subcubic graphs up to N vertices:
 (F1) phi(X) >= fr for all X;  maximum tight independent set has k=0 (non-K4);
 every NON-AUGMENTABLE tight independent set (Lemma 1 fails for all v) has k=0 unless G=K4;
 Lemma 1 itself (augmentation keeps tightness) on every tight set."""
import itertools, sys, networkx as nx
sys.path.insert(0,'/home/user/RhombileLattice/d3_complexity')
from exhaustive_heal import extend

def fr_sub(n, E, keep):
    # brute-force min mono edges on induced subgraph on 'keep' (list); returns (fr, list of optimal colourings as dicts)
    keep=list(keep); idx={v:i for i,v in enumerate(keep)}
    Es=[(idx[a],idx[b]) for a,b in E if a in idx and b in idx]
    m=len(keep)
    if m==0: return 0,[{}]
    best=None; cuts=[]
    for mask in range(1<<(m-1)):
        f=sum(1 for a,b in Es if ((mask>>a)&1)==((mask>>b)&1))
        if best is None or f<best: best=f; cuts=[mask]
        elif f==best: cuts.append(mask)
    return best,[{v:(c>>idx[v])&1 for v in keep} for c in cuts]

N=int(sys.argv[1])
level=[nx.path_graph(2)]; tot=0; viol=0
for n in range(3,N+1):
    level=extend(level)
    for G in level:
        V=list(G); E=list(G.edges()); tot+=1
        isK4=(n==4 and len(E)==6)
        fr,_=fr_sub(n,E,V)
        # F1 on all X
        cache={}
        for r in range(n+1):
            for X in itertools.combinations(V,r):
                k,_=fr_sub(n,E,[v for v in V if v not in X])
                cache[X]=k
                if r+k<fr: viol+=1; print("F1 violated",E,X)
        maxt=-1; 
        for X,k in cache.items():
            if any(G.has_edge(a,b) for a,b in itertools.combinations(X,2)): continue
            if len(X)+k!=fr: continue
            maxt=max(maxt,len(X))
            if k==0: continue
            Xs=set(X); rest=[v for v in V if v not in Xs]
            _,cuts=fr_sub(n,E,rest)
            aug=False
            for c in cuts:
                for a,b in E:
                    if a in Xs or b in Xs: continue
                    if c[a]==c[b]:
                        for v in (a,b):
                            if not any(w in Xs for w in G[v]):
                                aug=True
                                Y=tuple(sorted(X+(v,)))
                                # Lemma 1: S+v tight
                                if len(Y)+cache[Y]!=fr: viol+=1; print("Lemma1 violated",E,X,v)
            if not aug and not isK4: viol+=1; print("stuck non-K4",E,X)
        # max tight has k=0
        for X,k in cache.items():
            if len(X)==maxt and len(X)+k==fr and not any(G.has_edge(a,b) for a,b in itertools.combinations(X,2)):
                if k!=0 and not isK4: viol+=1; print("max tight k>0",E,X)
    print(f"n={n}: {len(level)} graphs, cumulative violations {viol}",flush=True)
print("total graphs",tot,"violations",viol)

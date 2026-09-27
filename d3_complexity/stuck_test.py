import itertools, random, networkx as nx, sys
def maxcuts(G):
    nodes=list(G); idx={v:i for i,v in enumerate(nodes)}; E=[(idx[a],idx[b]) for a,b in G.edges()]
    n=len(nodes); best=None; cuts=[]
    if n==0: return 0,[{}]
    for m in range(1<<(n-1)):
        f=sum(1 for a,b in E if ((m>>a)&1)==((m>>b)&1))
        if best is None or f<best: best=f; cuts=[m]
        elif f==best: cuts.append(m)
    return best,[{v:(m>>idx[v])&1 for v in nodes} for m in cuts]
def frv(G): return maxcuts(G)[0]
rng=random.Random(int(sys.argv[1])); checked=0; bad=0; stuckk=0
for trial in range(int(sys.argv[2])):
    n=rng.randint(6,12)
    G=nx.Graph(); G.add_nodes_from(range(n))
    if trial==0: G=nx.complete_graph(4); n=4
    for t in range(rng.randint(0,n//3)):
        a,b,c=rng.sample(range(n),3)
        if all(G.degree(x)<=1 for x in (a,b,c)): G.add_edges_from([(a,b),(b,c),(a,c)])
    for _ in range(3*n):
        a,b=rng.sample(range(n),2)
        if G.degree(a)<3 and G.degree(b)<3: G.add_edge(a,b)
    if not nx.is_connected(G): continue
    if trial==0: print('K4 run')
    fr=frv(G)
    for r in range(0,fr+1):
        for S in itertools.combinations(range(n),r):
            if any(G.has_edge(a,b) for a,b in itertools.combinations(S,2)): continue
            H=G.copy(); H.remove_nodes_from(S); k,cuts=maxcuts(H)
            assert k+r>=fr
            if k+r!=fr or k==0: continue
            Sset=set(S); aug=False
            for c in cuts:
                for a,b in H.edges():
                    if c[a]==c[b]:
                        for v in (a,b):
                            if not any(w in Sset for w in G[v]): aug=True
            if aug: continue
            stuckk+=1
            # check conclusion: frustrated comps are odd cycles, each vertex deg3 in G with exactly one S-nbr
            for comp in nx.connected_components(H):
                K=H.subgraph(comp)
                if frv(K)==0: continue
                ok = all(K.degree(v)==2 for v in K) and len(K)%2==1 and all(G.degree(v)==3 and sum(w in Sset for w in G[v])==1 for v in K)
                checked+=1
                if not ok: bad+=1; print("VIOLATION",sorted(G.edges()),S)
print("stuck tight sets with k>=1:",stuckk,"frustrated comps checked:",checked,"violations:",bad)

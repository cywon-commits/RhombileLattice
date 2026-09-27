"""MILP stress test of ioct = fr on larger subcubic graphs (beyond the n<=12 exhaustive range):
random cubic / subcubic graphs (n=14..60, triangle-rich variants) and the JMOPSvL G(phi) graphs."""
import random, sys, numpy as np, networkx as nx
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
sys.path.insert(0,'/home/user/RhombileLattice/d3_complexity')
from verify_delta3_reduction import build, randphi, sat

def fr_milp(n,E):
    N=n+len(E); c=np.zeros(N); c[n:]=1
    A=lil_matrix((2*len(E),N)); lo=[];hi=[]
    for k,(u,w) in enumerate(E):
        A[2*k,u]=1;A[2*k,w]=1;A[2*k,n+k]=1; lo.append(1);hi.append(np.inf)
        A[2*k+1,u]=1;A[2*k+1,w]=1;A[2*k+1,n+k]=-1; lo.append(-np.inf);hi.append(1)
    r=milp(c,constraints=LinearConstraint(A.tocsr(),lo,hi),integrality=np.ones(N),bounds=Bounds(0,1))
    return round(r.fun)
def ioct_milp(n,E):
    # vars: c_0..c_{n-1}, s_0..s_{n-1}
    N=2*n; c=np.zeros(N); c[n:]=1
    A=lil_matrix((3*len(E),N)); lo=[];hi=[]
    for k,(u,w) in enumerate(E):
        A[3*k,n+u]=1;A[3*k,n+w]=1; lo.append(-np.inf);hi.append(1)
        A[3*k+1,u]=1;A[3*k+1,w]=1;A[3*k+1,n+u]=1;A[3*k+1,n+w]=1; lo.append(1);hi.append(np.inf)
        A[3*k+2,u]=1;A[3*k+2,w]=1;A[3*k+2,n+u]=-1;A[3*k+2,n+w]=-1; lo.append(-np.inf);hi.append(1)
    r=milp(c,constraints=LinearConstraint(A.tocsr(),lo,hi),integrality=np.ones(N),bounds=Bounds(0,1))
    return None if r.status!=0 else round(r.fun)
def has_K4_comp(G):
    return any(len(C)==4 and G.subgraph(C).number_of_edges()==6 for C in nx.connected_components(G))
rng=random.Random(int(sys.argv[1])); random.seed(int(sys.argv[1]))
cnt=0; bad=0
def check(G,tag):
    global cnt,bad
    G=nx.convert_node_labels_to_integers(G); E=list(G.edges()); n=G.number_of_nodes()
    assert max(dict(G.degree()).values())<=3
    if has_K4_comp(G): return
    f=fr_milp(n,E); i=ioct_milp(n,E); cnt+=1
    if f!=i: bad+=1; print("MISMATCH",tag,n,f,i,E)
for t in range(int(sys.argv[2])):
    n=rng.choice(range(14,61,2))
    check(nx.random_regular_graph(3,n,seed=rng.randint(0,10**9)),"cubic")
    # triangle-rich subcubic
    m=rng.randint(14,50); G=nx.Graph(); G.add_nodes_from(range(m))
    for _ in range(m//3):
        a,b,c=rng.sample(range(m),3)
        if all(G.degree(x)<=1 for x in (a,b,c)): G.add_edges_from([(a,b),(b,c),(a,c)])
    for _ in range(3*m):
        a,b=rng.sample(range(m),2)
        if G.degree(a)<3 and G.degree(b)<3: G.add_edge(a,b)
    check(G,"tri")
    # replace vertices of a cubic graph by triangles (truncation) and subdivide a few edges
    H=nx.random_regular_graph(3,rng.choice(range(8,21,2)),seed=rng.randint(0,10**9))
    L=nx.line_graph(H)  # 4-regular; skip; use truncation instead
    T=nx.Graph()
    for v in H:
        nb=list(H[v])
        for i in range(3):
            for j in range(i+1,3): T.add_edge((v,nb[i]),(v,nb[j]))
    for u,v in H.edges(): T.add_edge((u,v),(v,u))
    for e in rng.sample(list(T.edges()),rng.randint(0,5)):
        if T.has_edge(*e): T.remove_edge(*e); x=('s',e); T.add_edge(e[0],x); T.add_edge(x,e[1])
    check(T,"trunc")
nsat=nuns=0
for t in range(int(sys.argv[3])):
    nv=rng.randint(3,7); m=rng.randint(3,min(10,2*nv))
    try: phi=randphi(nv,m)
    except ValueError: continue
    n,E,_=build(nv,phi); G=nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
    s=sat(nv,phi); nsat+=s; nuns+=(not s)
    check(G,"Gphi")
print("graphs checked",cnt,"mismatches",bad,"Gphi sat/unsat",nsat,nuns)

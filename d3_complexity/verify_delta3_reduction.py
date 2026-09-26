"""Exact cross-check of Theorem B1 (d3_complexity/theorem_delta3.md).

Builds the planar subcubic graph G(phi) of Johnson et al., arXiv:2211.12887,
Thm 11 (Planar 3-SAT, each literal in <= 2 clauses; 2m vertex-disjoint
triangles), solves min #mono edges + D3 * #colour-3 exactly by MILP (HiGHS),
and checks  OPT_D3(G) <= 2m*D3  <=>  phi satisfiable  for fixed D3 in (0,1).
Planarity is irrelevant to the equivalence, so random (non-planar) formulas
are used; max degree <= 3 is asserted.  This is a sanity check, not a proof.
"""
# Verify: for the JMOPSvL (2211.12887, Thm 11) reduction graph G(phi),
# OPT_D3(G) <= 2m*D3  <=>  phi satisfiable, for fixed D3 in (0,1).
import itertools, random, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

def build(nv, clauses):
    V=[]; E=set(); idx={}
    def node(name):
        if name not in idx: idx[name]=len(idx)
        return idx[name]
    def edge(a,b): E.add((min(node(a),node(b)),max(node(a),node(b))))
    for i in range(nv): edge(('x',i,1),('x',i,0))
    m=len(clauses)
    for k in range(2*m-1): edge(('P',k),('P',k+1))
    for j,cl in enumerate(clauses):
        ins=[('x',abs(l)-1,1 if l>0 else 0) for l in cl]
        if len(ins)==2: ins.append(('u',j))
        a,b,d,e,f,c=[('g',j,t) for t in 'abdefc']
        for p,q in [(ins[0],a),(ins[1],b),(ins[2],f),(a,b),(a,d),(b,d),(d,e),(e,f),(f,c),(e,c)]: edge(p,q)
        edge(c,('P',2*j))            # j-th odd vertex (1-indexed odd) = false
        if len(cl)==2: edge(ins[2],('P',2*j+1))   # j-th even vertex = true
    n=len(idx); return n,sorted(E),idx

def opt(n,E,D3):
    nx_=3*n; ne=len(E); N=nx_+ne
    c=np.zeros(N); c[[3*v+2 for v in range(n)]]=D3; c[nx_:]=1
    rows=n+3*ne; A=lil_matrix((rows,N)); lo=np.zeros(rows); hi=np.zeros(rows)
    for v in range(n):
        A[v,3*v:3*v+3]=1; lo[v]=hi[v]=1
    r=n
    for k,(u,w) in enumerate(E):
        for col in range(3):
            A[r,3*u+col]=1; A[r,3*w+col]=1; A[r,nx_+k]=-1; lo[r]=-np.inf; hi[r]=1; r+=1
    res=milp(c,constraints=LinearConstraint(A.tocsr(),lo,hi),integrality=np.ones(N),bounds=Bounds(0,1))
    return res.fun

def sat(nv,clauses):
    return any(all(any((l>0)==bool(a[abs(l)-1]) for l in cl) for cl in clauses) for a in itertools.product((0,1),repeat=nv))

def randphi(nv,m):
    for _ in range(2000):
        occ={}; cls=[]
        ok=True
        for j in range(m):
            k=min(nv,random.choice([2,2,3]))
            vs=random.sample(range(1,nv+1),k)
            cl=[v*random.choice([1,-1]) for v in vs]
            for l in cl: occ[l]=occ.get(l,0)+1
            cls.append(cl)
        if max(occ.values())<=2: return cls
    raise ValueError


if __name__ == "__main__":
    random.seed(11); cases=[]; tries=0
    while len(cases) < 80 and tries < 50000:
        tries += 1; nv = random.randint(2, 5); m = random.randint(2, min(7, 2*nv))
        try: phi = randphi(nv, m)
        except ValueError: continue
        cases.append((nv, phi))
    nsat = nunsat = bad = 0; mingap = float("inf")
    for nv, phi in cases:
        n, E, _ = build(nv, phi); m = len(phi); s = sat(nv, phi)
        assert max(sum(1 for e in E if v in e) for v in range(n)) <= 3
        nsat += s; nunsat += (not s)
        for D3 in (0.05, 0.3, 0.6, 0.9, 0.99):
            o = opt(n, E, D3)
            if (o <= 2*m*D3 + 1e-7) != s:
                bad += 1; print("MISMATCH", phi, D3, o)
            if not s: mingap = min(mingap, o - 2*m*D3)
    print(f"formulas: {nsat} sat, {nunsat} unsat; mismatches: {bad}; min unsat gap: {mingap:.4f}")

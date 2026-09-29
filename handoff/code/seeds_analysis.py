import numpy as np,glob,re
from collections import defaultdict
from annealed_long_analysis import tau_int
byk=defaultdict(lambda: defaultdict(list))
for f in sorted(glob.glob('results/long/long_L*_d21_mu[35]_T*_s*.npz')):
    m=re.search(r'L(\d+)_d21_mu(\d+)_T([\d.]+)_s(\d+)',f); L,mu,T=int(m[1]),int(m[2]),float(m[3])
    if T not in (0.78,0.83,0.88,0.86,0.93,1.0) or L<32: continue
    d=np.load(f); k=100000//int(d['thin'])
    p2=d['sre'][k:,0].astype(float)**2+d['sim'][k:,0].astype(float)**2
    err=np.sqrt(p2.var()*2*tau_int(p2)/len(p2))
    byk[(mu,T)][L].append((p2.mean(),err))
for (mu,T),dd in sorted(byk.items()):
    print(f"mu={mu} T={T}"); prev=None
    for L in sorted(dd):
        a=np.array(dd[L]); w=1/a[:,1]**2; p=(a[:,0]*w).sum()/w.sum(); e=np.sqrt(1/w.sum())
        chi=((a[:,0]-p)**2*w).sum()/(len(a)-1) if len(a)>1 else np.nan
        e2=e*np.sqrt(max(1,chi)) if len(a)>1 else e
        s=" ".join(f"{x:.4f}" for x in a[:,0])
        if prev:
            eta=np.log(prev[1]/p)/np.log(L/prev[0]); de=np.hypot(prev[2]/prev[1],e2/p)/np.log(L/prev[0])
            es=f"eta={eta:+.3f}±{de:.3f}"
        else: es=""
        print(f"  L={L:3d} seeds[{s}] mean={p:.4f}±{e2:.4f} chi2/dof={chi:.1f} {es}")
        prev=(L,p,e2)
print("# weighted fit psi2 ~ L^-eta over all L>=32 (errors inflated by sqrt(chi2/dof) of seed scatter)")
for (mu,T),dd in sorted(byk.items()):
    x=[];y=[];s=[]
    for L in sorted(dd):
        a=np.array(dd[L]); w=1/a[:,1]**2; p=(a[:,0]*w).sum()/w.sum(); e=np.sqrt(1/w.sum())
        if len(a)>1: e*=np.sqrt(max(1,((a[:,0]-p)**2*w).sum()/(len(a)-1)))
        x.append(np.log(L)); y.append(np.log(p)); s.append(e/p)
    x,y,s=map(np.array,(x,y,s)); W=np.diag(1/s**2); X=np.vstack([np.ones_like(x),x]).T
    C=np.linalg.inv(X.T@W@X); b=C@X.T@W@y; chi=((y-X@b)**2/s**2).sum()/max(1,len(x)-2)
    print(f"mu={mu} T={T}: eta={-b[1]:.3f}±{np.sqrt(C[1,1]*max(1,chi)):.3f} chi2/dof={chi:.2f} (L={len(x)} sizes)")

import re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":8.5,"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
BLUE,ORG,INK2,GRID="#2a78d6","#eb6834","#52514e","#e4e3df"
def load(f,model):
    d=[]
    for l in open(f):
        m=re.match(rf"{model} L=24 T=([\d.]+) Tchi=([\d.]+) \(([\d.]+)\)",l)
        if m: d.append([float(x) for x in m.groups()])
    return np.array(d)
a=load("chi_b.log","annealed"); t=load("chi_a.log","tafm"); chi0=0.041
fig,(p,q)=plt.subplots(1,2,figsize=(7.0,2.8))
TT=np.logspace(np.log10(0.18),1.05,200)
for d,c,mk,lb in [(a,BLUE,"o","annealed bonds"),(t,ORG,"s","fixed-bond TAFM")]:
    p.errorbar(d[:,0],d[:,1]/d[:,0],yerr=d[:,2]/d[:,0],fmt=mk,ms=3.5,color=c,label=lb,zorder=3)
    q.errorbar(d[:,0],d[:,1],yerr=d[:,2],fmt=mk+"-",ms=3.5,lw=1,color=c,label=lb,zorder=3)
p.plot(TT,chi0/TT,color=INK2,ls="--",lw=1); p.text(0.21,chi0/0.21*1.25,r"$\chi_0/T$ (ground-state manifold)",fontsize=7,color=INK2)
p.plot(TT[TT>1.4],0.25/(TT[TT>1.4]+3),color=ORG,ls=":",lw=1)
p.plot(TT[TT>1.4],0.25/(TT[TT>1.4]+2.4),color=BLUE,ls=":",lw=1)
p.text(6.0,0.25/(6+2.4)*1.45,"Curie–Weiss",fontsize=7,color=INK2,ha="center")
p.set_xscale("log"); p.set_yscale("log"); p.set_xlabel("T"); p.set_ylabel(r"$\chi=-\partial\langle n_1\rangle/N\partial D_2$"); p.legend(loc="lower left",frameon=False); p.grid(color=GRID,lw=0.5)
q.axhline(chi0,color=INK2,ls="--",lw=0.8); q.text(4.0,chi0+0.006,r"$\chi_0=0.041$",fontsize=7,color=INK2)
q.axhline(0.25,color=INK2,ls=":",lw=0.8); q.text(0.2,0.238,"free spins 1/4",fontsize=7,color=INK2,va="top")
q.set_xscale("log"); q.set_xlabel("T"); q.set_ylabel(r"$T\chi$"); q.set_ylim(0,0.27); q.legend(loc="center left",frameon=False); q.grid(color=GRID,lw=0.5)
for ax,l in [(p,"(a)"),(q,"(b)")]: ax.text(-0.16,1.02,l,transform=ax.transAxes,fontsize=10,fontweight="bold")
fig.tight_layout(w_pad=2); fig.savefig("/mnt/user-data/outputs/chi/chi_vs_T.png",dpi=200)

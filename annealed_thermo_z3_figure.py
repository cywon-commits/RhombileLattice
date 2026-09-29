"""Z3 crystallisation of the annealed model for 0=D1<D2 (results/thermo/z3_*).

psi = sum_k m_k w^k, m_k = hub-type fraction on triangular sublattice k.
Disordered: <|psi|^2> ~ 1/N, so <|psi|^2> N is size-independent.
Ordered (long-range order): <|psi|^2> finite, <|psi|^2> N grows ~ N.
Binder U4 = 1 - <|psi|^4>/(2<|psi|^2>^2): 0 for a Gaussian 2-component
order parameter, 1/2 for long-range order."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(p):
    return np.array([[float(x) for x in l.split()] for l in open(p) if l[0].isdigit()])


cols = {8: "#c5791f", 12: "#3d7a4f", 16: "#3b4ba8", 24: "#a83b3b"}
fig, ax = plt.subplots(2, 3, figsize=(14, 7.6))
for r, (d2, xl) in enumerate([("1", (0.3, 1.0)), ("0.5", (0.08, 0.6))]):
    for L, c in cols.items():
        d = load(f"results/thermo/z3_L{L}_d2{d2}.log")
        T, N = d[:, 0], 3 * L * L
        ax[r, 0].semilogy(T, d[:, 10] * N, "o-", ms=3, lw=1, color=c, label=f"L={L}")
        ax[r, 1].plot(T, d[:, 11], "o-", ms=3, lw=1, color=c, label=f"L={L}")
        ax[r, 2].plot(T, d[:, 3], "o-", ms=3, lw=1, color=c, label=f"L={L}")
    ax[r, 0].set_ylabel("<|psi|^2> x N_site")
    ax[r, 1].set_ylabel("Binder U4")
    ax[r, 1].axhline(0.5, color="k", ls=":", lw=0.8)
    ax[r, 1].axhline(0, color="k", lw=0.5)
    ax[r, 2].set_ylabel("C / N_tri")
    for a in ax[r]:
        a.set_xlim(*xl)
        a.set_xlabel("T")
        a.set_title(f"0=D1<D2={d2}", fontsize=10)
        a.legend(fontsize=7)
fig.suptitle("Annealed model, D3=10: Z3 order parameter of the rhombile crystal (which sublattice hosts the hubs)")
fig.tight_layout()
fig.savefig("annealed_thermo_z3_figure.png", dpi=150)
print("saved annealed_thermo_z3_figure.png")

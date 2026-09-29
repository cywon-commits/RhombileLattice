"""Illustrative figures for the Part IX paper draft.

fig_snapshots.png : Monte Carlo configurations (L=12): (a) rhombile ground state,
    (b) D1=D2 at T=0.40 (random tiling, TAFM-like spins), (c) D2=1 at T=0.50
    (sqrt3 x sqrt3 crystal), (d) D2=1 at T=0.60 (disordered). Lozenges = pairs of
    matched triangles (the shared inactive bond is not drawn); red triangles =
    frustrated (unmatched); dots = spins (filled = state 1, open = state 0).
fig_mechanism.png : why a spin defect costs twice as much with annealed bonds.

Usage: python3 paper_figures.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from annealed_thermo import build_arrays, sweep, total_energy, _seed
from rhombile_lattice import RhombileLattice, A_MAT, A_MAT_INV
from worm_monomer_walk import DimerState

L = 12
D3 = 10.0


def unwrap_tri(pts, box):
    """Bring the 3 corners of a triangle next to the first one (minimum image)."""
    f = (A_MAT_INV @ pts.T).T
    d = f[1:] - f[0]
    d -= np.round(d / box) * box
    f = np.vstack([f[0], f[0] + d])
    return (A_MAT @ f.T).T


def simulate(D2, T, n, seed, mu=0.0):
    A = build_arrays(L)
    _seed(seed)
    rng = np.random.default_rng(seed)
    s = rng.integers(2, size=A["ns"]).astype(np.int64)
    dim, mate = A["dimer"].copy(), A["mate"].copy()
    if D2 > 0:
        s = np.where(A["sub"] == 0, 1, 0).astype(np.int64)
    E = total_energy(s, dim, A["bi"], A["bj"], D2, D3, mu, mate)
    for Tk in np.linspace(max(T, 1.0), T, 10):
        for _ in range(n // 10):
            E = sweep(s, dim, mate, A["bi"], A["bj"], A["bta"], A["btb"], A["site_b"], A["hexb"], 1 / Tk, D2, D3, mu, E)
    return A, s, dim, mate


def draw(ax, A, s, dim, mate, title):
    st = DimerState(L, L)
    pos, _ = st.lat.site_positions()
    tri_sites = [t["sites"] for t in st.tri]
    # draw each matched pair as one lozenge, frustrated triangles in red
    done = set()
    for t, sites in enumerate(tri_sites):
        pts = unwrap_tri(pos[list(sites)], L)
        c = pts.mean(0)
        f = A_MAT_INV @ c
        if not (0 <= f[0] < L and 0 <= f[1] < L):
            continue
        if mate[t] < 0:
            ax.add_patch(Polygon(pts, closed=True, fc="#d9544d", ec="#444", lw=0.6, alpha=0.85))
        else:
            k = mate[t]
            u = A["btb"][k] if A["bta"][k] == t else A["bta"][k]
            if (min(t, u), max(t, u)) in done:
                continue
            done.add((min(t, u), max(t, u)))
            q = unwrap_tri(pos[list(tri_sites[u])], L)
            # shift u's triangle next to t's
            fu = A_MAT_INV @ q.mean(0) - A_MAT_INV @ c
            q = q - (A_MAT @ (np.round(fu / L) * L))
            allp = np.vstack([pts, q])
            # the lozenge = convex hull of the 4 distinct corners
            uniq = []
            for p in allp:
                if not any(np.hypot(*(p - v)) < 1e-6 for v in uniq):
                    uniq.append(p)
            uniq = np.array(uniq)
            cen = uniq.mean(0)
            order = np.argsort(np.arctan2(uniq[:, 1] - cen[1], uniq[:, 0] - cen[0]))
            ax.add_patch(Polygon(uniq[order], closed=True, fc="#eef2f7", ec="#444", lw=0.6))
    fr = (A_MAT_INV @ pos.T).T
    inside = (fr[:, 0] < L) & (fr[:, 1] < L)
    one = (s == 1) & inside
    zero = (s == 0) & inside
    ax.scatter(pos[one, 0], pos[one, 1], s=9, c="#1d2a44", zorder=3)
    ax.scatter(pos[zero, 0], pos[zero, 1], s=9, facecolors="white", edgecolors="#1d2a44", lw=0.6, zorder=3)
    corners = np.array([[0, 0], [L, 0], [L, L], [0, L]]) @ A_MAT.T
    ax.set_xlim(corners[:, 0].min() - 0.5, corners[:, 0].max() + 0.5)
    ax.set_ylim(-0.5, corners[:, 1].max() + 1)
    ax.set_aspect("equal")
    ax.axis("off")
    nm = (mate < 0).mean()
    ax.set_title(f"{title}\nfrustrated fraction {nm:.3f}", fontsize=10)


def snapshots():
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
    ax = axs.ravel()
    A = build_arrays(L)
    s0 = np.where(A["sub"] == 0, 1, 0).astype(np.int64)
    draw(ax[0], A, s0, A["dimer"], A["mate"], "(a) T=0, 0<D2<6: rhombile crystal")
    draw(ax[1], *simulate(0.0, 0.40, 20000, 3), "(b) D1=D2, T=0.40: random tiling")
    draw(ax[2], *simulate(1.0, 0.50, 20000, 4), "(c) D2=1, T=0.50 < T_c: sqrt3 crystal")
    draw(ax[3], *simulate(1.0, 0.60, 20000, 5), "(d) D2=1, T=0.60 > T_c: disordered")
    fig.tight_layout()
    fig.savefig("figures/fig_snapshots.png", dpi=150)


def mechanism():
    """Two minimal pictures: the same spin defect in the fixed-bond TAFM and in
    the annealed model. A triangle whose three spins are equal (defect) shown
    with its bonds; equal-spin active bonds (cost 1 each) in red, the dimer
    (inactive) dashed."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    h = np.sqrt(3) / 2
    tri = np.array([[0, 0], [1, 0], [0.5, h]])

    def base(a, title):
        a.set_aspect("equal")
        a.axis("off")
        a.set_title(title, fontsize=11)

    # (a) TAFM ground-state triangle: one equal bond, cost 1 per triangle-bond (shared), the 'minimal' frustration
    a = ax[0]
    base(a, "(a) TAFM: every triangle has one equal bond\n(fixed bonds: that bond always costs)")
    for (i, j), eq in (((0, 1), True), ((1, 2), False), ((2, 0), False)):
        a.plot(*tri[[i, j]].T, color="#d9544d" if eq else "#444", lw=3 if eq else 1.5)
    for p, v in zip(tri, (1, 1, 0)):
        a.scatter(*p, s=260, c="#1d2a44" if v else "white", edgecolors="#1d2a44", zorder=3)
    a.text(0.5, -0.22, "defect (all three equal): +2 per pair", ha="center", fontsize=10)

    # (b) annealed: the equal bond is the dimer -> free; a defect triangle keeps 2 exposed equal bonds
    a = ax[1]
    base(a, "(b) Annealed: the dimer hides the equal bond\n(ground state costs 0)")
    for (i, j), kind in (((0, 1), "dimer"), ((1, 2), "ok"), ((2, 0), "ok")):
        if kind == "dimer":
            a.plot(*tri[[i, j]].T, color="#888", lw=2, ls="--")
        else:
            a.plot(*tri[[i, j]].T, color="#444", lw=1.5)
    for p, v in zip(tri, (1, 1, 0)):
        a.scatter(*p, s=260, c="#1d2a44" if v else "white", edgecolors="#1d2a44", zorder=3)
    a.text(0.5, -0.22, "dashed = dimer (inactive bond)", ha="center", fontsize=10)

    a = ax[2]
    base(a, "(c) Annealed defect: only one bond can be hidden\n-> 2 exposed equal bonds, pair cost >= 4")
    for (i, j), kind in (((0, 1), "dimer"), ((1, 2), "bad"), ((2, 0), "bad")):
        if kind == "dimer":
            a.plot(*tri[[i, j]].T, color="#888", lw=2, ls="--")
        else:
            a.plot(*tri[[i, j]].T, color="#d9544d", lw=3)
    for p in tri:
        a.scatter(*p, s=260, c="#1d2a44", edgecolors="#1d2a44", zorder=3)
    a.text(0.5, -0.22, "red = equal-spin active bond (cost 1)", ha="center", fontsize=10)
    for a in ax:
        a.set_xlim(-0.3, 1.3)
        a.set_ylim(-0.35, 1.2)
    fig.tight_layout()
    fig.subplots_adjust(top=0.80)
    fig.savefig("figures/fig_mechanism.png", dpi=150)


if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)
    snapshots()
    mechanism()
    print("wrote figures/fig_snapshots.png, figures/fig_mechanism.png")

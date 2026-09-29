"""D_c (last kink of E(D3)) versus frustrated-triangle density p, from every
construction used: max-flow isolated monomers, greedy jamming, worm-built
configurations (envelope kinks), and the fully frustrated lattice (exact)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pts = [  # (p, D_c, label, marker, colour)
    (0.0404, 2.06, "worm-built (envelope kink)", "o", "#3b4ba8"),
    (0.0794, 1.89, None, "o", "#3b4ba8"),
    (0.1198, 2.26, None, "o", "#3b4ba8"),
    (0.08, 1.9, "max-flow isolated monomers", "s", "#3d7a4f"),
    (0.12, 2.0, "greedy jamming (p=0.120)", "^", "#c5791f"),
    (1.0, 3.0, "fully frustrated (exact level crossing)", "*", "#a83b3b"),
]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for p, d, lab, m, c in pts:
    ax.plot(p, d, m, ms=11 if m == "*" else 8, color=c, label=lab)
ax.axhline(2, color="k", ls="--", lw=1, label="dilute string/domino value 2")
ax.axhline(3, color="k", ls=":", lw=1, label="junction / saturated value 3")
ax.axvspan(0.125, 1.0, color="#e4ddcd", alpha=0.5, lw=0)
ax.text(0.4, 2.5, "not reached by our isolated-monomer\nconstructions (all are lower bounds)", ha="center", fontsize=8, color="#6b6558")
ax.set_xscale("log")
ax.set_xlim(0.03, 1.3)
ax.set_ylim(1.5, 3.3)
ax.set_xlabel("frustrated-triangle density p (fraction of triangles)")
ax.set_ylabel("D_c (last kink of E(D3))")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig("Dc_vs_density_figure.png", dpi=150)
print("saved Dc_vs_density_figure.png")

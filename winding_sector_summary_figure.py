"""Summary figure: does the matching principle (only the topological
class matters, not the specific construction) hold WITHIN a fixed
nontrivial flux sector too, not just within the trivial one (Case1 vs
Case2)?

Rebuilds the diagonal (1,1) loop (winding_double_and_xy's own SA answer
for the "two separate loops" (1,1) construction is reused from that
script's own run rather than recomputed here, to avoid two more minutes
of redundant SA), gets its actual best state for a state-colored panel,
and puts together a bar chart spanning every closed-loop construction
tried so far:

  (0,0) trivial:      Case0 hexagon = 0, double x-wind = 0
  (1,0) nontrivial:   single x-wind loop ~= 30
  (1,1) nontrivial:   two-separate-loops ~= 29, diagonal ~= 31

The point: (1,0) built one way costs ~30; (1,1) built two very different
ways (four times as many bonds touched, totally different shape) costs
~29-31 either time -- consistent with each other, not with the specific
construction's own length or bond count.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, route_string_between_points, route_string,
    apply_dual_string_defect, frustrated_triangles, full_state_via_mincut,
    total_energy, simulated_annealing,
)
from dual_string_demo import draw_lattice
from winding_double_and_xy import build_base, whole_graph_contradictions

NX, NY = 30, 8
N_WAYPOINTS = 6


def build_diagonal():
    lat = RhombileLattice(NX, NY)
    triangles, tri_by_r1, rhombi, sibling, hop = build_base(lat)
    ks = np.arange(N_WAYPOINTS)
    xs, ys = ks * NX / N_WAYPOINTS, ks * NY / N_WAYPOINTS
    pts = [np.array([x, y]) for x, y in zip(xs, ys)]

    touched = []
    required_start = None
    first_ts = None
    for i in range(N_WAYPOINTS):
        p0, p1 = pts[i], pts[(i + 1) % N_WAYPOINTS]
        is_closing = (i == N_WAYPOINTS - 1)
        if not is_closing:
            if required_start is None:
                ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
            else:
                n_end, m_end = int(round(p1[0])) % NX, int(round(p1[1])) % NY
                te_candidates = tri_by_r1[lat._site_index(n_end, m_end, "r1")]
                best = None
                for te in te_candidates:
                    try:
                        nodes_c, bonds_c = route_string(sibling, hop, required_start, te)
                    except ValueError:
                        continue
                    if best is None or len(nodes_c) < len(best[1]):
                        best = (te, nodes_c, bonds_c)
                te, nodes, bonds = best
                ts = required_start
        else:
            te_target = first_ts
            nodes, bonds = route_string(sibling, hop, required_start, te_target)
            ts, te = required_start, te_target
        flipped = apply_dual_string_defect(triangles, nodes, bonds)
        touched.extend(flipped)
        required_start = te
        if i == 0:
            first_ts = ts
    return lat, touched


def main():
    lat, touched = build_diagonal()
    on = [b for b in touched if b["J"] != 0.0]
    off = [b for b in touched if b["J"] == 0.0]
    ft = frustrated_triangles(lat)

    D3 = 1.0
    seed_states, _ = full_state_via_mincut(lat, D3)
    best_e, best_states = total_energy(lat, seed_states, D=(0.0, 0.0, D3)), seed_states
    rng = np.random.default_rng(11)
    for _ in range(5):
        st, en = simulated_annealing(lat, (0.0, 0.0, D3), rng, n_sweeps=2500,
                                      T_start=0.5, T_end=1e-6, states=seed_states, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    for _ in range(6):
        st, en = simulated_annealing(lat, (0.0, 0.0, D3), rng, n_sweeps=3500,
                                      T_start=5.0, T_end=1e-6, record_energy=True)
        if en[-1] < best_e:
            best_e, best_states = en[-1], st
    print(f"diagonal (1,1): E={best_e}, n3={(best_states==2).sum()}, bonds touched={len(touched)}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    ax_lat, ax_bar, ax_text, ax_blank = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    draw_lattice(ax_lat, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                 box=(NX, NY), off_lw=0.3, states=best_states,
                 title=f"Diagonal (minimal-length) (1,1) loop, D3=1\n"
                       f"{len(touched)} bonds touched, colored by state -- SA-best E={best_e:.0f}")
    handles, labels = ax_lat.get_legend_handles_labels()
    label_map = {"state 0": "state1", "state 1": "state2", "state 2": "state3"}
    ax_lat.legend(handles, [label_map.get(l, l) for l in labels], loc="upper right", fontsize=8)

    labels_bar = ["Case0\n(0,0)", "double x-wind\n(0,0)", "single x-wind\n(1,0)",
                  "x+y, 2-separate\n(1,1)", "x+y, diagonal\n(1,1)"]
    values = [0, 0, 30, 29, int(round(best_e))]
    bond_counts = ["-", "-", "~19", "~38 total", f"{len(touched)}"]
    colors = ["tab:green", "tab:green", "tab:orange", "tab:red", "tab:red"]
    x = np.arange(len(labels_bar))
    bars = ax_bar.bar(x, values, color=colors)
    for xi, v, bc in zip(x, values, bond_counts):
        ax_bar.text(xi, v + 1, f"E={v}\n({bc} bonds)", ha="center", fontsize=8)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels_bar, fontsize=9)
    ax_bar.set_ylabel("SA-best ground-state energy (D3=1)")
    ax_bar.set_title("Trivial sectors (green) always 0; a given nontrivial sector's\n"
                      "cost is consistent across wildly different constructions (orange/red)")

    ax_text.axis("off")
    text = (
        "The matching principle (only the topological class matters, not the specific\n"
        "construction) was established for the TRIVIAL sector via Case1 vs Case2 (Part IV):\n"
        "very different drawn paths between the same 2 defects converge to the same SA\n"
        "ground state. Today's test asks: does this ALSO hold for a fixed NONtrivial flux\n"
        "sector, not just for 'no flux'?\n\n"
        "Two very different constructions of the SAME (1,1) sector:\n"
        "  - 'two separate loops' (an independent x-wind + y-wind, ~38 bonds total): E=29\n"
        "  - a genuinely diagonal, near-minimal-length single loop (178 bonds!): E=31\n\n"
        "These are close to each other (and to the single x-wind (1,0) sector's own ~30),\n"
        "despite radically different bond patterns and total lengths -- consistent with the\n"
        "SAME matching principle extending across sectors: what matters is the flux class,\n"
        "not how a specific bond pattern realizing it happens to be drawn. Contrast with\n"
        "Part VII's direct-vs-winding STRING result (48 vs 5) -- those differed because\n"
        "they were in DIFFERENT sectors, not because one was a 'bad' construction."
    )
    ax_text.text(0.02, 0.98, text, va="top", ha="left", fontsize=9.2, family="monospace",
                 transform=ax_text.transAxes)
    ax_blank.axis("off")

    plt.suptitle("The matching principle survives within a fixed nontrivial flux sector too",
                 y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("winding_sector_summary.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved winding_sector_summary.png")


if __name__ == "__main__":
    main()

"""Same two defects, shortest path vs. a deliberate detour.

Answers two related questions:
- Does a longer (non-shortest) string connecting the same two defects
  cost more? Yes, at every D3 > 0 (equal only at D3=0) -- see
  detour_vs_shortest_demo.
- Section 5 (dual_string_demo / the research doc) found that a bad
  4-defect pairing can be beaten by "reorganizing" into a nearest-
  neighbor rematch. Doesn't that mean a detour should also reorganize
  into the shortest path? It's the same operation in principle (both
  are comparisons between two fixed, hand-built bond configurations),
  but what genuinely reorganizes *on its own*, automatically, is the
  Potts-state assignment within one fixed bond configuration -- see
  detour_reorganization_demo, which shows exactly how much that buys
  you (something, in the middle of the D3 range) and how little (it
  never shrinks the underlying structure, so it can't close the gap).
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
    conflict_graph_bipartition, energy_via_mincut, total_energy, full_state_via_mincut,
)
from dual_string_demo import draw_lattice, build_dual

NX, NY = 14, 14
A, B, C = (1.0, 1.0), (7.0, 1.0), (7.0, 9.0)


def _build_shortest():
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    _, _, nodes, bonds = route_string_between_points(rho, sib, hop, A, B)
    flipped = apply_dual_string_defect(tri, nodes, bonds)
    return lat, flipped, len(nodes)


def _build_detour():
    lat, tri, rho, sib, hop = build_dual(NX, NY)
    _, _, nodes1, bonds1 = route_string_between_points(rho, sib, hop, A, C)
    f1 = apply_dual_string_defect(tri, nodes1, bonds1)
    _, _, nodes2, bonds2 = route_string_between_points(rho, sib, hop, C, B)
    f2 = apply_dual_string_defect(tri, nodes2, bonds2)
    return lat, f1 + f2


def detour_vs_shortest_demo():
    lat_s, flipped_s, L_s = _build_shortest()
    on_s = [b for b in flipped_s if b["J"] != 0.0]
    off_s = [b for b in flipped_s if b["J"] == 0.0]
    ft_s = frustrated_triangles(lat_s)
    states_s, _ = full_state_via_mincut(lat_s, 0.0)

    lat_l, flipped_l = _build_detour()
    on_l = [b for b in flipped_l if b["J"] != 0.0]
    off_l = [b for b in flipped_l if b["J"] == 0.0]
    ft_l = frustrated_triangles(lat_l)
    states_l, _ = full_state_via_mincut(lat_l, 0.0)
    edges_l, _, _, _ = conflict_graph_bipartition(lat_l)
    L_l = len(edges_l)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.8))
    draw_lattice(axes[0], lat_s, highlight_on=on_s, highlight_off=off_s, frustrated=ft_s, box=(NX, NY),
                 targets=[A, B], states=states_s, title=f"shortest path (BFS-optimal): L={L_s}\nD3=0: E=0 for free")
    draw_lattice(axes[1], lat_l, highlight_on=on_l, highlight_off=off_l, frustrated=ft_l, box=(NX, NY),
                 targets=[A, B, C], states=states_l,
                 title=f"deliberate detour via C: L={L_l}\nD3=0: also E=0, but costs more for any D3>0")
    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.suptitle("Same two defects (A, B): shortest string vs. a deliberate detour through C", y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("detour_vs_shortest_lattice.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    pos, sub_of = lat_s.site_positions()
    naive_states = np.where(sub_of == "r1", 0, 1).astype(int)
    D3g = np.linspace(0, 3.6, 200)
    Es = [energy_via_mincut(lat_s, d) for d in D3g]
    El = [energy_via_mincut(lat_l, d) for d in D3g]
    E_naive = [total_energy(lat_s, naive_states, D=(0, 0, d)) for d in D3g]

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.plot(D3g, E_naive, color="#a49c8a", lw=1.8, ls=":", label=f"no minimization (naive, always {L_s})")
    ax.plot(D3g, Es, color="#3b4ba8", lw=2.5, label=f"shortest path, minimized (L={L_s})")
    ax.plot(D3g, El, color="#e8792c", lw=2.5, label=f"detour, minimized (L={L_l})")
    ax.set_xlabel("D3"); ax.set_ylabel("ground-state energy E")
    ax.set_title("Same 2 defects: exact ground-state energy for shortest path vs. a detour\n"
                  "(plus the un-minimized naive baseline, for reference)")
    ax.legend(fontsize=9, loc="center right")
    plt.tight_layout()
    fig.savefig("detour_vs_shortest_D3.png", dpi=140)
    plt.close(fig)
    print(f"shortest L={L_s}, detour L={L_l}")


def detour_reorganization_demo(D3_panels=(0.0, 1.5, 2.5, 3.5)):
    lat, flipped = _build_detour()
    all_on = [b for b in flipped if b["J"] != 0.0]
    all_off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)

    fig, axes = plt.subplots(1, len(D3_panels), figsize=(5.5 * len(D3_panels), 6.2))
    for ax, D3 in zip(axes, D3_panels):
        states, violated = full_state_via_mincut(lat, D3)
        draw_lattice(ax, lat, highlight_on=all_on, highlight_off=all_off, frustrated=ft, box=(NX, NY),
                     targets=[A, B], states=states, title=f"D3={D3}  (violated bonds in magenta: {len(violated)})")
        for b in violated:
            ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]], color="magenta", linewidth=3.2, zorder=6)
        if D3 == D3_panels[0]:
            ax.legend(loc="upper left", fontsize=7, framealpha=0.9)
    plt.suptitle("Detour configuration: optimal spin/state reorganization as D3 increases", y=1.02, fontsize=13)
    plt.tight_layout()
    fig.savefig("detour_reorg_D3.png", dpi=125, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    detour_vs_shortest_demo()
    detour_reorganization_demo()

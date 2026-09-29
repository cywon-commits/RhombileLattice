"""Step 1 of testing the finite-T entropic-force idea (research_notes_
next_directions.md, section 14.1): before any finite-temperature
calculation, first nail down the exact T=0 ground-state energy as a
function of defect separation L, for two very different separations,
using the straight-string regime that's already solidly confirmed
(section 6/12 of the report -- hub-fixed is exact here, verified
against SA at multiple D3 points).

Case A: short separation (~5 units). Case B: long separation (~20
units, 4x farther). Both are plain horizontal strings at the same
y_cross used throughout this project's straight-string demos (2.1 --
an arbitrary but already-validated height, not a special/degenerate
one; a different height silently produced a non-canonical, ~1
diagonal-per-unit conflict graph in an earlier attempt, so this
specific height is not interchangeable without re-checking).

Result: L_B/L_A = 4.0 exactly, and E_B(D3)/E_A(D3) = 4.0 exactly at
every D3 tested (0 through past both curves' own saturation) -- not
just a small-D3 approximation. This is the clean baseline the entropy
calculation (S(L) from counting paths, then F(L,T) = E(L) - T*S(L))
will need to compare against.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, frustrated_triangles,
    conflict_graph_bipartition, energy_via_mincut, full_state_via_mincut,
)
from dual_string_demo import draw_lattice

NX, NY = 30, 8
Y_CROSS = 2.1
DIST_A, DIST_B = 5.0, 20.0


def build(distance):
    lat = RhombileLattice(NX, NY)
    apply_string_defect(lat, np.array([2.0, Y_CROSS]), np.array([2.0 + distance, Y_CROSS]), wrap=False)
    return lat


def short_vs_long_demo():
    lat_A, lat_B = build(DIST_A), build(DIST_B)
    ft_A, ft_B = frustrated_triangles(lat_A), frustrated_triangles(lat_B)
    edges_A, _, _, nb_A = conflict_graph_bipartition(lat_A)
    edges_B, _, _, nb_B = conflict_graph_bipartition(lat_B)
    L_A, L_B = len(edges_A), len(edges_B)
    print(f"Case A: distance={DIST_A}, L={L_A}, defects={len(ft_A)}, non_bipartite={nb_A}")
    print(f"Case B: distance={DIST_B}, L={L_B}, defects={len(ft_B)}, non_bipartite={nb_B}")
    print(f"L ratio: {L_B / L_A}")

    print("D3      E_A     E_B     ratio")
    for D3 in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        e_A, e_B = energy_via_mincut(lat_A, D3), energy_via_mincut(lat_B, D3)
        ratio = e_B / e_A if e_A > 0 else float("inf")
        print(f"{D3:<7} {e_A:<7} {e_B:<7} {ratio:.3f}")

    states_A, _ = full_state_via_mincut(lat_A, 0.0)
    states_B, _ = full_state_via_mincut(lat_B, 0.0)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7))
    draw_lattice(axes[0], lat_A, frustrated=ft_A, box=(NX, NY), states=states_A, off_lw=0.5,
                 title=f"Case A: short distance ~{DIST_A:.0f} units (L={L_A})")
    axes[0].set_xlim(-1, 12); axes[0].set_ylim(0, 5)
    draw_lattice(axes[1], lat_B, frustrated=ft_B, box=(NX, NY), states=states_B, off_lw=0.5,
                 title=f"Case B: long distance ~{DIST_B:.0f} units (L={L_B})")
    axes[1].set_xlim(-1, 27); axes[1].set_ylim(0, 5)
    plt.tight_layout()
    fig.savefig("short_vs_long_lattice.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    D3_grid = np.linspace(0, 4, 200)
    E_A = [energy_via_mincut(lat_A, d) for d in D3_grid]
    E_B = [energy_via_mincut(lat_B, d) for d in D3_grid]

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.plot(D3_grid, E_A, color="#3b4ba8", lw=2.5, label=f"Case A (L={L_A}, distance~{DIST_A:.0f})")
    ax.plot(D3_grid, E_B, color="#8f2d56", lw=2.5, label=f"Case B (L={L_B}, distance~{DIST_B:.0f})")
    ax.set_xlabel("D3"); ax.set_ylabel("ground-state energy E")
    ax.set_title("Short vs. long defect separation: E(D3) scales EXACTLY 4x at every D3\n"
                  "(both linear regime and after saturation)")
    ax.legend(fontsize=10)
    plt.tight_layout()
    fig.savefig("short_vs_long_D3.png", dpi=140)
    plt.close(fig)
    print("saved short_vs_long_lattice.png, short_vs_long_D3.png")


if __name__ == "__main__":
    short_vs_long_demo()

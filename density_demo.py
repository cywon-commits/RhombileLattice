"""Defect-density test: tile a lattice with many short strings and see
whether the D3=2 threshold found for a single isolated string (see
dual_string_demo docstrings) survives once there are many defects.

Short answer (see rhombile_lattice.energy_via_mincut for the exact
method): the threshold moves to D3=3 and the sharp kink smooths into a
multi-segment curve, but this is driven by the *branching* (degree-3
vertices) in each short string's own conflict subgraph, not by density
as such -- a sparse grid of the same short strings shows almost the
same curve as a dense one.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, enumerate_triangles, build_rhombi, build_triangle_hop_graph,
    route_string_between_points, apply_dual_string_defect, frustrated_triangles,
    conflict_graph_bipartition, energy_via_mincut,
)
from dual_string_demo import draw_lattice

NX, NY = 16, 16


def build_dense_grid(nx, ny, spacing, dx=0.6):
    """Tile the interior with short (minimal, 2-defect) strings on a
    regular grid, `spacing` apart. Returns (lattice, anchors)."""
    lat = RhombileLattice(nx, ny)
    triangles = enumerate_triangles(lat)
    rhombi, _ = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)
    anchors, y = [], 1.5
    while y < ny - 1.5:
        x = 1.5
        while x < nx - 1.5:
            anchors.append((x, y))
            x += spacing
        y += spacing
    for (x, y) in anchors:
        _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, (x, y), (x + dx, y))
        apply_dual_string_defect(triangles, nodes, bonds)
    return lat, anchors


def density_comparison_demo():
    """Sparse vs. dense grid: defect count falls well short of the naive
    2-per-string count at high density, because neighboring strings'
    diagonal flips overlap and cancel in pairs."""
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.3))
    for ax, spacing in zip(axes, [4.0, 1.0]):
        lat, anchors = build_dense_grid(NX, NY, spacing)
        # bonds actually flipped: everything that ended up wrong or right of default
        on = [b for b in lat.bonds if b["type"] == "r2r3" and b["J"] != 0.0]
        off = [b for b in lat.bonds if b["type"] in ("r1r2", "r1r3") and b["J"] == 0.0]
        ft = frustrated_triangles(lat)
        print(f"spacing={spacing}: strings={len(anchors)} defects={len(ft)} "
              f"(naive expectation {2 * len(anchors)})")
        draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft, box=(NX, NY),
                     title=f"spacing={spacing}: {len(anchors)} short strings -> {len(ft)} actual defects\n"
                           f"(naive expectation: {2 * len(anchors)})")
    plt.suptitle("Dense-packing test: many short strings on a regular grid -- "
                  "at high density, defects annihilate in pairs", y=1.0, fontsize=13)
    plt.tight_layout()
    fig.savefig("density_comparison.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def density_D3_sweep_demo():
    """Exact E(D3)/L for an isolated string vs. sparse/dense grids of
    short strings, via energy_via_mincut (exact, polynomial-time)."""
    D3_grid = np.linspace(0, 3.2, 200)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.plot(D3_grid, np.minimum(D3_grid * 7, 14) / 14, color="#a49c8a", lw=2, ls=":",
            label="single isolated string (L=14): sharp kink at D3=2")
    for spacing, color in [(4.0, "#3b4ba8"), (1.0, "#e8792c")]:
        lat, anchors = build_dense_grid(NX, NY, spacing)
        edges, side_a, side_b, non_bip = conflict_graph_bipartition(lat)
        L = len(edges)
        assert non_bip == 0, "expected a bipartite conflict graph"
        Es = [energy_via_mincut(lat, d) for d in D3_grid]
        ax.plot(D3_grid, np.array(Es) / L, color=color, lw=2.3,
                label=f"{'sparse' if spacing > 2 else 'dense'} grid "
                      f"({len(anchors)} strings, L={L})")
    ax.axvline(2.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.axvline(3.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(2.02, 0.04, "D3=2", fontsize=8, color="gray")
    ax.text(3.02, 0.04, "D3=3", fontsize=8, color="gray")
    ax.set_xlabel("D3"); ax.set_ylabel("E(D3) / L  (fraction of full-Ising energy)")
    ax.set_title("More strings -> branchier conflict graph -> threshold shifts from D3=2 to D3=3,\n"
                  "and the sharp kink smooths into a multi-segment curve")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig("density_D3_sweep.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    density_comparison_demo()
    density_D3_sweep_demo()

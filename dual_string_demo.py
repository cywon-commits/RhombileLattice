"""Systematic string construction via the triangle-hop dual graph: instead
of a straight geometric line (rhombile_lattice.apply_string_defect), a
string is routed as a shortest path through the lattice's rhombi, so its
two endpoints -- and therefore its two resulting defects -- can be placed
at *any* two locations, not just ones a straight line happens to connect.

See rhombile_lattice.build_triangle_hop_graph / route_string for the
derivation: crossing a rhombus's diagonal always lands on its sibling
half, and the next hop must leave via an outer edge, so every interior
half keeps exactly 2 active edges (just a different one inactive) and
only the path's two ends end up with all 3 edges active -- frustrated --
regardless of how much the path bends.
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import (
    RhombileLattice, apply_string_defect, frustrated_triangles,
    enumerate_triangles, build_rhombi, build_triangle_hop_graph,
    route_string_between_points, apply_dual_string_defect,
)

NX, NY = 10, 10


STATE_COLORS = {0: "#8f2d56", 1: "#2a6f97", 2: "#e8792c"}  # r1-fixed, rim-default, escaped


def draw_lattice(ax, lattice, highlight_on=None, highlight_off=None, frustrated=None,
                  targets=None, title="", states=None, box=None):
    """Same convention as string_defect_demo.draw_lattice (black = active
    bonds, gold stars = frustrated triangle centroids), plus: red solid =
    bonds in `highlight_on` (rhombus diagonals just turned on), orange
    dashed = bonds in `highlight_off` (outer edges just turned off --
    drawn explicitly since they're no longer J!=0 and would otherwise
    vanish, breaking the string's visual continuity). `targets`, if
    given, marks the originally requested (possibly off-lattice) points
    with orange x's. `states`, if given (array of 0/1/2 per site), colors
    every site by its Potts *state* instead of by sublattice (r1/r2/r3) --
    use this to show an actual ground-state assignment (e.g. the minimum
    state-2 set that zeroes out a string's diagonal conflicts). `box`, if
    given as (nx, ny), sets the axis limits for a lattice of that size
    instead of the module-level NX, NY (needed for any lattice built at a
    different size than the rest of this file's demos)."""
    on_ids = set(id(b) for b in highlight_on) if highlight_on else set()
    for b in lattice.bonds:
        if b["J"] != 0.0:
            is_hl = id(b) in on_ids
            color = "red" if is_hl else "black"
            lw = 2.2 if is_hl else 0.8
            ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                     color=color, linewidth=lw, zorder=2 if is_hl else 1)
    for b in (highlight_off or []):
        ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                 color="darkorange", linewidth=1.8, linestyle="--", zorder=2)

    pos, sub_of = lattice.site_positions()
    if states is None:
        colors = {"r1": "tab:red", "r2": "tab:blue", "r3": "tab:green"}
        for sub, c in colors.items():
            mask = sub_of == sub
            ax.scatter(pos[mask, 0], pos[mask, 1], s=18, color=c, zorder=3)
    else:
        states = np.asarray(states)
        for k, c in STATE_COLORS.items():
            mask = states == k
            ax.scatter(pos[mask, 0], pos[mask, 1], s=18, color=c, zorder=3,
                       label=f"state {k}")

    if frustrated:
        centroids = np.array([f["centroid"] for f in frustrated])
        ax.scatter(centroids[:, 0], centroids[:, 1], marker="*", s=260,
                   color="gold", edgecolor="black", linewidth=0.8, zorder=5)

    if targets:
        pts = np.array(targets)
        ax.scatter(pts[:, 0], pts[:, 1], marker="x", s=90, color="orange",
                   linewidth=2.2, zorder=4)

    box_nx, box_ny = box if box is not None else (NX, NY)
    ax.set_aspect("equal")
    ax.set_xlim(-1, box_nx - 1)
    ax.set_ylim(-1, box_ny - 1)
    ax.set_title(title)


def build_dual(nx=NX, ny=NY):
    lat = RhombileLattice(nx, ny)
    triangles = enumerate_triangles(lat)
    rhombi, _ = build_rhombi(triangles)
    sibling, hop = build_triangle_hop_graph(triangles)
    return lat, triangles, rhombi, sibling, hop


def straight_line_vs_dual_demo(y_cross=2.1, x0=0.5, x1=7.5):
    """Sanity check: the new dual-graph router, given the same two
    endpoints as the original straight-line string, should reproduce the
    same qualitative result -- 2 frustrated triangles at the same spots."""
    lat_line = RhombileLattice(NX, NY)
    flipped_line, _ = apply_string_defect(
        lat_line, np.array([x0, y_cross]), np.array([x1, y_cross]), wrap=False)
    on_line = [b for b in flipped_line if b["J"] != 0.0]
    off_line = [b for b in flipped_line if b["J"] == 0.0]
    ft_line = frustrated_triangles(lat_line)

    lat_dual, triangles, rhombi, sibling, hop = build_dual()
    _, _, nodes, bonds = route_string_between_points(
        rhombi, sibling, hop, (x0, y_cross), (x1, y_cross))
    flipped_dual = apply_dual_string_defect(triangles, nodes, bonds)
    on_dual = [b for b in flipped_dual if b["J"] != 0.0]
    off_dual = [b for b in flipped_dual if b["J"] == 0.0]
    ft_dual = frustrated_triangles(lat_dual)

    print(f"[straight-line]  flipped={len(flipped_line)}  frustrated={len(ft_line)}")
    print(f"[dual-graph]     flipped={len(flipped_dual)}  frustrated={len(ft_dual)}  "
          f"rhombi visited={len(nodes)}")

    fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
    draw_lattice(axes[0], lat_line, highlight_on=on_line, highlight_off=off_line,
                 frustrated=ft_line, targets=[(x0, y_cross), (x1, y_cross)],
                 title=f"straight-line method\n{len(ft_line)} defects")
    draw_lattice(axes[1], lat_dual, highlight_on=on_dual, highlight_off=off_dual,
                 frustrated=ft_dual, targets=[(x0, y_cross), (x1, y_cross)],
                 title=f"dual-graph method (same endpoints)\n{len(ft_dual)} defects, "
                       f"{len(nodes)} rhombi visited")
    plt.tight_layout()
    fig.savefig("dual_string_straight_comparison.png", dpi=150)
    plt.close(fig)


def bent_string_demo(p_start=(1.0, 1.0), p_end=(7.0, 6.5)):
    """The actual point of the exercise: two arbitrary, non-collinear
    points, connected by a routed (bent) string, with defects landing
    exactly there instead of wherever a straight line happens to fall."""
    lat, triangles, rhombi, sibling, hop = build_dual()
    ts, te, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p_start, p_end)
    flipped = apply_dual_string_defect(triangles, nodes, bonds)
    on = [b for b in flipped if b["J"] != 0.0]
    off = [b for b in flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    print(f"[bent string {p_start} -> {p_end}] rhombi visited={len(nodes)} "
          f"flipped={len(flipped)} frustrated={len(ft)} "
          f"defect centroids={[tuple(np.round(f['centroid'], 2)) for f in ft]}")

    fig, ax = plt.subplots(figsize=(9, 8))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                 targets=[p_start, p_end],
                 title=f"arbitrarily-placed defect pair via a routed (bent) string\n"
                       f"{len(nodes)} rhombi visited, {len(ft)} defects "
                       "(orange x = requested points, dashed = broken outer edges)")
    plt.tight_layout()
    fig.savefig("dual_string_bent.png", dpi=150)
    plt.close(fig)


def multi_string_demo(pairs=(((1.0, 1.0), (3.0, 2.0)),
                              ((7.0, 1.0), (8.3, 3.0)),
                              ((1.5, 7.5), (4.0, 8.3)))):
    """Several independently-routed strings placing several arbitrary
    defect pairs at once. Strings that end up sharing a bond partially
    cancel/reconnect where they cross -- a real interaction, not a bug,
    the same phenomenon as string_defect_demo.two_parallel_strings_demo."""
    lat, triangles, rhombi, sibling, hop = build_dual()
    all_flipped, targets = [], []
    for p0, p1 in pairs:
        _, _, nodes, bonds = route_string_between_points(rhombi, sibling, hop, p0, p1)
        all_flipped.extend(apply_dual_string_defect(triangles, nodes, bonds))
        targets.extend([p0, p1])
    # classify by *final* J: a bond crossed by 2 strings gets toggled twice,
    # so its status right after either individual call can be stale.
    all_on = [b for b in all_flipped if b["J"] != 0.0]
    all_off = [b for b in all_flipped if b["J"] == 0.0]
    ft = frustrated_triangles(lat)
    print(f"[{len(pairs)} independent strings] requested defects={2 * len(pairs)} "
          f"actual frustrated={len(ft)} "
          f"(fewer means some pairs' paths crossed and partially cancelled)")

    fig, ax = plt.subplots(figsize=(9, 8))
    draw_lattice(ax, lat, highlight_on=all_on, highlight_off=all_off, frustrated=ft,
                 targets=targets,
                 title=f"{len(pairs)} arbitrarily-placed string defects at once\n"
                       f"requested {2 * len(pairs)} defect points -> {len(ft)} actual defects")
    plt.tight_layout()
    fig.savefig("dual_string_multi.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    straight_line_vs_dual_demo()
    bent_string_demo()
    multi_string_demo()

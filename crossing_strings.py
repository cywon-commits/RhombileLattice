"""Tests the user's question: can two OPEN strings be built so their drawn
paths actually cross each other? Prediction: apply_dual_string_defect just
XORs bond states, so unless the two strings happen to toggle a literally
SHARED bond at the crossing, they act as two fully independent defect
pairs -- 4 frustrated triangles total (string A's 2 endpoints + string B's
2 endpoints), with the crossing point itself looking perfectly ordinary
(no special cancellation or extra defect there), unlike a CLOSED loop's
corner (which cancels by construction because the same triangle-half is
deliberately shared between consecutive segments).
"""
import numpy as np
import matplotlib.pyplot as plt

from rhombile_lattice import route_string_between_points, apply_dual_string_defect, frustrated_triangles
from dual_string_demo import build_dual, draw_lattice

NX, NY = 20, 20


def main():
    lat, triangles, rhombi, sibling, hop = build_dual(NX, NY)

    # String A: roughly horizontal, string B: roughly vertical, crossing
    # near the box center.
    a_left = np.array([3.0, 10.0])
    a_right = np.array([16.0, 9.0])
    b_top = np.array([9.5, 16.0])
    b_bottom = np.array([9.5, 3.0])

    _, _, nodes_a, bonds_a = route_string_between_points(rhombi, sibling, hop, a_left, a_right)
    touched_a = apply_dual_string_defect(triangles, nodes_a, bonds_a)

    _, _, nodes_b, bonds_b = route_string_between_points(rhombi, sibling, hop, b_top, b_bottom)
    touched_b = apply_dual_string_defect(triangles, nodes_b, bonds_b)

    shared = set(id(b) for b in touched_a) & set(id(b) for b in touched_b)
    print(f"string A touched {len(touched_a)} bonds, string B touched {len(touched_b)} bonds")
    print(f"bonds touched by BOTH strings (would cancel via double-toggle): {len(shared)}")

    ft = frustrated_triangles(lat)
    print(f"frustrated triangles: {len(ft)} (predicted 4: 2 endpoints per string, "
          f"independent since no shared bond)")
    for f in ft:
        print(f"  defect at centroid {np.round(f['centroid'], 2)}")

    on = [b for b in (touched_a + touched_b) if b["J"] != 0.0]
    off = [b for b in (touched_a + touched_b) if b["J"] == 0.0]
    fig, ax = plt.subplots(figsize=(8, 8))
    draw_lattice(ax, lat, highlight_on=on, highlight_off=off, frustrated=ft,
                 targets=[a_left, a_right, b_top, b_bottom], box=(NX, NY),
                 title=f"Two crossing open strings: {len(ft)} frustrated triangles")
    plt.tight_layout()
    plt.savefig("crossing_strings.png", dpi=140, bbox_inches="tight")
    print("saved crossing_strings.png")


if __name__ == "__main__":
    main()

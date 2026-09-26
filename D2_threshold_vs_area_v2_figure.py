"""Analysis + figure for the corrected (non-confounded) D2-threshold vs
area/curvature sweep (D2_threshold_vs_area_v2.py). Loads whichever
D2_threshold_v2_{label}.pkl files are present (robust to a still-running
sweep), extracts a threshold D2* per config (first D2 where n_state3 jumps
out of the low ~4-6 plateau into a higher plateau, if any such jump exists
strictly below D2=1.0=D3), and plots:
  (1) n_state3 vs D2 for all AREA configs (h2,h4,h6,h9,h13)
  (2) n_state3 vs D2 for all CURVATURE configs (skew0,skew2,skew4 @ area=10.5)
  (3) extracted D2* vs area (log-log, with a 1/Area reference line)
"""
import glob
import pickle

import matplotlib.pyplot as plt
import numpy as np

AREA_LABELS = [("area_h2", 3.5), ("area_h4", 7.0), ("area_h6", 10.5),
               ("area_h9", 15.75), ("area_h13", 22.75)]
CURV_LABELS = [("curv_skew0", 10.5), ("curv_skew2", 10.5), ("curv_skew4", 10.5)]


def load(label):
    try:
        with open(f"D2_threshold_v2_{label}.pkl", "rb") as f:
            return pickle.load(f)
    except FileNotFoundError:
        return None


def extract_threshold(results):
    """First D2 (strictly < 1.0) where n3 jumps to a new plateau above the
    initial low-D2 baseline plateau, if any. Returns None if flat throughout
    (0, 1.0)."""
    if results is None:
        return None
    D2s = sorted(results["D2"].keys())
    D2s = [d for d in D2s if d < 0.999]
    n3s = [results["D2"][d]["n3"] for d in D2s]
    if len(n3s) < 2:
        return None
    baseline = n3s[1] if len(n3s) > 1 else n3s[0]  # skip D2=0 (often degenerate)
    for d, n in zip(D2s[1:], n3s[1:]):
        if n > baseline + 3:
            return d
    return None


def main():
    area_data = {lab: load(lab) for lab, _ in AREA_LABELS}
    curv_data = {lab: load(lab) for lab, _ in CURV_LABELS}

    present_area = [(lab, a) for lab, a in AREA_LABELS if area_data[lab] is not None]
    present_curv = [(lab, a) for lab, a in CURV_LABELS if curv_data[lab] is not None]
    print(f"Loaded {len(present_area)}/{len(AREA_LABELS)} area configs, "
          f"{len(present_curv)}/{len(CURV_LABELS)} curvature configs")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    for lab, area in present_area:
        r = area_data[lab]
        D2s = sorted(r["D2"].keys())
        n3s = [r["D2"][d]["n3"] for d in D2s]
        ax.plot(D2s, n3s, "o-", label=f"{lab} (A={area:.2f})")
    ax.set_xlabel("D2")
    ax.set_ylabel("n_state3")
    ax.set_title("AREA sweep: n_state3 vs D2")
    ax.legend(fontsize=8)
    ax.set_yscale("symlog")

    ax = axes[1]
    for lab, area in present_curv:
        r = curv_data[lab]
        D2s = sorted(r["D2"].keys())
        n3s = [r["D2"][d]["n3"] for d in D2s]
        ax.plot(D2s, n3s, "o-", label=lab)
    ax.set_xlabel("D2")
    ax.set_ylabel("n_state3")
    ax.set_title(f"CURVATURE sweep (area={present_curv[0][1] if present_curv else 10.5} fixed)")
    ax.legend(fontsize=8)
    ax.set_yscale("symlog")

    ax = axes[2]
    areas, thresholds = [], []
    for lab, area in present_area:
        thr = extract_threshold(area_data[lab])
        print(f"  {lab}: area={area:.2f}  D2*={thr}")
        if thr is not None:
            areas.append(area)
            thresholds.append(thr)
    if areas:
        ax.plot(areas, thresholds, "o", ms=10, color="C0", label="measured D2*")
        a_ref = np.linspace(min(areas), max(areas), 50)
        c = thresholds[0] * areas[0]
        ax.plot(a_ref, c / a_ref, "k--", alpha=0.5, label="1/Area reference")
    ax.set_xlabel("enclosed area")
    ax.set_ylabel("D2* threshold")
    ax.set_title("Threshold vs area")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("D2_threshold_vs_area_v2_figure.png", dpi=130)
    print("Saved D2_threshold_vs_area_v2_figure.png")

    for lab, area in present_curv:
        thr = extract_threshold(curv_data[lab])
        print(f"  {lab}: area={area:.2f}  D2*={thr}")


if __name__ == "__main__":
    main()

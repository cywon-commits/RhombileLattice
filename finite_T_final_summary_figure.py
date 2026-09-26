"""Final summary figure for the whole finite-T re-verification thread:
A-vs-B confinement crossover (v2, more rigorous), O's own hub/rim
order-disorder transition (fine-resolved, plus the D3=3 comparison), and
the exact T->infinity limit (dA=dB=1/3) derived from bond counting.
"""
import csv

import numpy as np
import matplotlib.pyplot as plt


def load(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


ab = load("finite_T_scan_v2_results.csv")
o_coarse = load("hub_rim_order_transition_results.csv")
o_fine = load("hub_rim_order_transition_fine_results.csv")
o_d3_3 = load("hub_rim_order_transition_D3_3.0_results.csv")

o_all = {row["T"]: row["abs_m_mean"] for row in o_coarse}
o_all.update({row["T"]: row["abs_m_mean"] for row in o_fine})
o_T = sorted(o_all.keys())
o_m = [o_all[t] for t in o_T]

d3_3_T = [row["T"] for row in o_d3_3]
d3_3_m = [row["abs_m_mean"] for row in o_d3_3]

ab_T = [row["T"] for row in ab]
dA = [row["dA_mean"] for row in ab]
dA_se = [row["dA_se"] for row in ab]
dB = [row["dB_mean"] for row in ab]
dB_se = [row["dB_se"] for row in ab]
ratio = [row["ratio"] for row in ab]
ratio_se = [row["ratio_se"] for row in ab]

Tc_theory = 2 / np.arccosh((1 + np.sqrt(3)) / 2) * 0.5
Tc_D3_1_empirical = 0.79

fig, axes = plt.subplots(3, 1, figsize=(10, 13), sharex=True)
ax1, ax2, ax3 = axes

ax1.plot(o_T, o_m, "o-", color="tab:blue", label="D3=1 (fine)", markersize=4)
ax1.plot(d3_3_T, d3_3_m, "s-", color="tab:green", label="D3=3", markersize=4)
ax1.axvline(Tc_theory, color="tab:red", ls="--", label=f"theory T_c={Tc_theory:.3f} (D3->inf limit)")
ax1.axvline(Tc_D3_1_empirical, color="tab:blue", ls=":", alpha=0.7, label=f"empirical T_c(D3=1)~{Tc_D3_1_empirical}")
ax1.set_ylabel("<|m|> (hub order parameter)")
ax1.set_title("O (pristine): hub/rim order collapses at higher T as D3 rises toward the pure-Ising limit")
ax1.legend(fontsize=8)

ax2.errorbar(ab_T, dA, yerr=dA_se, marker="o", color="tab:orange", label="dA (short, L~5)")
ax2.errorbar(ab_T, dB, yerr=dB_se, marker="s", color="tab:purple", label="dB (long, L~20, 4x)")
ax2.axhline(1/3, color="gray", ls="--", label="exact T->inf limit = 1/3 (bond-count argument)")
ax2.axvline(Tc_D3_1_empirical, color="tab:blue", ls=":", alpha=0.7)
ax2.set_yscale("log")
ax2.set_ylabel("excess energy vs pristine reference")
ax2.set_title("A vs B excess energy: both converge to the SAME exact value 1/3 as T->infinity")
ax2.legend(fontsize=8)

ax3.errorbar(ab_T, ratio, yerr=ratio_se, marker="o", color="tab:red")
ax3.axhline(4.0, color="gray", ls=":", label="T=0 exact ratio = 4")
ax3.axhline(1.0, color="gray", ls="--", label="T->inf exact ratio = 1")
ax3.axvline(Tc_D3_1_empirical, color="tab:blue", ls=":", alpha=0.7, label="O's own T_c(D3=1)")
ax3.set_xlabel("T")
ax3.set_ylabel("confinement ratio dB/dA")
ax3.set_title("Confinement ratio drops sharply right at O's own order-disorder T_c")
ax3.legend(fontsize=8)

plt.tight_layout()
fig.savefig("finite_T_final_summary.png", dpi=140, bbox_inches="tight")
print("saved finite_T_final_summary.png")

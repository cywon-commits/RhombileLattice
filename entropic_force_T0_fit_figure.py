"""Figure for the refined T=0 entropic-force functional-form check: with
5 seeds per separation (up from 2), fits log(L) and a free power law
c*L^p to d(L), and shows the power law (exponent ~0.99, i.e. close to
LINEAR in L) fits notably better than Moore et al.'s originally-claimed
log(r) form (reduced chi^2 0.46 vs 1.16).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

L = np.array([12, 16, 20, 28, 36], dtype=float)
d = np.array([5.495, 7.131, 7.133, 11.282, 15.613])
sem = np.array([1.321, 0.825, 1.383, 1.871, 1.468])


def log_model(L, a, b):
    return a + b * np.log(L)


def power_model(L, c, p):
    return c * np.power(L, p)


popt_log, _ = curve_fit(log_model, L, d, p0=[0, 3], sigma=sem, absolute_sigma=True)
popt_pow, _ = curve_fit(power_model, L, d, p0=[1, 0.5], sigma=sem, absolute_sigma=True)

chi2_log = np.sum(((d - log_model(L, *popt_log)) / sem) ** 2) / (len(L) - 2)
chi2_pow = np.sum(((d - power_model(L, *popt_pow)) / sem) ** 2) / (len(L) - 2)

L_fine = np.linspace(10, 40, 300)

fig, ax = plt.subplots(figsize=(6.5, 5))
ax.errorbar(L, d, yerr=sem, fmt="o", ms=7, color="#201d18", ecolor="#6b6558",
            capsize=3, zorder=5, label="data (5-seed avg ± SEM)")
ax.plot(L_fine, log_model(L_fine, *popt_log), color="#c5791f", lw=2, ls="--",
        label=f"log fit: a+b·ln(L)  (reduced χ²={chi2_log:.2f})")
ax.plot(L_fine, power_model(L_fine, *popt_pow), color="#3d7a4f", lw=2,
        label=f"power fit: c·L^{popt_pow[1]:.2f}  (reduced χ²={chi2_pow:.2f})")
ax.set_xlabel("Defect separation L")
ax.set_ylabel("Excess local escape-site count d(L)")
ax.set_title("T=0 entropic force: power law (~linear) fits better than log(r)")
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig("entropic_force_T0_fit_figure.png", dpi=150)
print("saved entropic_force_T0_fit_figure.png")
print(f"log fit params: a={popt_log[0]:.3f}, b={popt_log[1]:.3f}, reduced chi2={chi2_log:.3f}")
print(f"power fit params: c={popt_pow[0]:.3f}, p={popt_pow[1]:.3f}, reduced chi2={chi2_pow:.3f}")

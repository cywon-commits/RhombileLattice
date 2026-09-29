"""Clock anisotropy of the spin order parameter psi_s = |psi| e^{i theta} (minority-spin fraction per
sublattice; sre/sim of annealed_long.py). theta = 0, 120, 240 deg are the three sqrt3 x sqrt3 crystals
(minority on one sublattice); theta = 60, 180, 300 deg their spin-inverted partners.
Per (L, T): <cos 3theta>, <cos 6theta>, and |psi|-weighted versions, with time-block errors (NB blocks).
Long-range order: <cos 6theta> (and <cos 3theta>) -> const with L. Critical phase with irrelevant
anisotropy: both decay with L as a power.
Usage: [MU=10] [DISCARD=100000] [NB=8] python3 results/long/clock_analysis.py
"""
import glob
import os
import re
from collections import defaultdict

import numpy as np

MU = os.environ.get("MU", "10")
DISCARD = int(os.environ.get("DISCARD", 100000))
NB = int(os.environ.get("NB", 8))

data = defaultdict(dict)
for f in glob.glob(f"results/long/long_L*_d21_mu{MU}_T*_s*.npz"):
    m = re.search(r"L(\d+)_d21_mu[\d.]+_T([\d.]+)_s(\d+)\.npz", f)
    if not m:
        continue
    L, T = int(m[1]), float(m[2])
    d = np.load(f)
    k = DISCARD // int(d["thin"])
    z = d["sre"][k:, 0].astype(float) + 1j * d["sim"][k:, 0].astype(float)
    r, th = np.abs(z), np.angle(z)
    obs = {"c3": np.cos(3 * th), "c6": np.cos(6 * th), "rc6": r * np.cos(6 * th), "r": r}
    n = len(z) // NB * NB
    blk = {key: v[:n].reshape(NB, -1).mean(1) for key, v in obs.items()}
    out = {}
    for key in ("c3", "c6"):
        out[key] = (blk[key].mean(), blk[key].std(ddof=1) / np.sqrt(NB))
    w = blk["rc6"] / blk["r"]                       # <|psi| cos 6theta>/<|psi|>
    out["w6"] = (w.mean(), w.std(ddof=1) / np.sqrt(NB))
    data[T][L] = out

print(f"# mu={MU}: <cos 3theta>, <cos 6theta>, <|psi|cos6theta>/<|psi|>  (time-block errors, NB={NB})")
for T in sorted(data):
    print(f"T={T:.3f}")
    for L in sorted(data[T]):
        o = data[T][L]
        print("  L=%3d  c3=%+.3f(%3.0f)  c6=%+.3f(%3.0f)  w6=%+.3f(%3.0f)" % (
            L, o["c3"][0], 1e3 * o["c3"][1], o["c6"][0], 1e3 * o["c6"][1], o["w6"][0], 1e3 * o["w6"][1]))

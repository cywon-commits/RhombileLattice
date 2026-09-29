"""Lower envelope of E(D3) from every SA state found by worm_Dc_vs_density.py.

Each state found at any D3 is a line E = n3*D3 + v (v = violated bonds),
valid at every D3, so the pointwise minimum over all of them is a tighter
upper bound on the true E(D3) than any single SA run (which also fixes the
SA runs' non-monotone E(D3)). Kinks of the envelope are the D_c values.

Usage: python3 worm_Dc_envelope.py results/worm_Dc_p0.04.log [...]
"""
import re
import sys


def analyze(path):
    lines = []
    for ln in open(path):
        m = re.match(r"D3=([\d.]+)\s+E=([\d.]+)\s+n_state3_sites=(\d+)", ln)
        if m:
            d3, e, n3 = float(m[1]), float(m[2]), int(m[3])
            lines.append((d3, e, n3, e - d3 * n3))
    states = {}
    for _, _, n3, v in lines:
        states[n3] = min(v, states.get(n3, float("inf")))
    grid = [x / 100 for x in range(50, 351)]
    arg = [min(states.items(), key=lambda s: s[0] * d + s[1])[0] for d in grid]
    kinks = [(grid[i], arg[i - 1], arg[i]) for i in range(1, len(grid)) if arg[i] != arg[i - 1]]
    print(f"== {path} ==  E_inf <= {states.get(0)}")
    for d3, e, _, _ in lines:
        env = min(k * d3 + w for k, w in states.items())
        print(f"  D3={d3:.2f}  SA E={e:6.1f}  envelope E={env:6.1f}")
    print("  kinks (D3, n3 before -> after):", kinks)


if __name__ == "__main__":
    for path in sys.argv[1:]:
        analyze(path)

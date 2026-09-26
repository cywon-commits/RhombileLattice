"""Self-tests for gadget.py: known examples + random cross-check against brute force."""
import itertools, random
from gadget import Gadget, value, describe, min_gap, breakpoints

def brute(n, E, T, b, D3):
    best = None
    free = [v for v in range(n) if v not in T]
    for cs in itertools.product((1, 2, 3), repeat=len(free)):
        s = [0] * n
        for v, c in zip(T, b): s[v] = c
        for v, c in zip(free, cs): s[v] = c
        c = sum(s[u] == s[w] for u, w in E) + D3 * sum(s[v] == 3 for v in free)
        best = c if best is None else min(best, c)
    return best

# 1. triangle, no terminals: min(1, D3)
tri = Gadget(3, [(0, 1), (1, 2), (0, 2)], [])
f = tri.frontier(())
assert f == {0: 1, 1: 0}, f
# 2. wheel W4: min(2, D3), hub is the unique colour-3 choice for D3 < 2
w4 = Gadget(5, [(0, i) for i in range(1, 5)] + [(1, 2), (2, 3), (3, 4), (4, 1)], [], names={0: "hub", 1: "r1", 2: "r2", 3: "r3", 4: "r4"})
assert w4.frontier(()) == {0: 2, 1: 0}
rep = w4.choice_report((), 1.5)
assert rep["always3"] == ["hub"] and rep["sometimes3"] == [], rep
# 3. random cross-check
random.seed(0)
for trial in range(150):
    n = random.randint(3, 8)
    E = [e for e in itertools.combinations(range(n), 2) if random.random() < 0.45]
    T = random.sample(range(n), random.randint(0, min(3, n)))
    g = Gadget(n, E, T)
    for b in itertools.product((1, 2, 3), repeat=len(T)):
        fr = g.frontier(b)
        for D3 in (0.0, 0.3, 1.0, 1.7, 2.5, 4.0):
            assert abs(value(fr, D3) - brute(n, E, T, b, D3)) < 1e-9
print("gadget.py self-tests passed")

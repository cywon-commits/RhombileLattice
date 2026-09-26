"""Cost-model table of the Garey-Johnson-Stockmeyer OR (clause) gadget used in theorem_delta3.md.

Terminals: inputs c1,c2,c3 and the reference vertex ref (colour F); inputs coloured T or F.
Result (see output): every satisfied input pattern costs min(2, 1+D3, 2*D3), FFF costs
min(3, 2+D3, 1+2*D3), so on D3 in (0,1) the penalty for FFF is exactly 1.
"""
import itertools
from gadget import Gadget, describe, min_gap

names = ['c1', 'c2', 'c3', 'ref', 'a', 'b', 'd', 'e', 'f', 'c']
ix = {s: i for i, s in enumerate(names)}
E = [('c1', 'a'), ('c2', 'b'), ('c3', 'f'), ('a', 'b'), ('a', 'd'), ('b', 'd'), ('d', 'e'),
     ('e', 'f'), ('f', 'c'), ('e', 'c'), ('c', 'ref')]
g = Gadget(len(names), [(ix[u], ix[v]) for u, v in E], [ix[x] for x in ('c1', 'c2', 'c3', 'ref')],
           names=dict(enumerate(names)))
T, F = 1, 2
fronts = {}
for bits in itertools.product((T, F), repeat=3):
    lab = ''.join('T' if x == T else 'F' for x in bits)
    fronts[lab] = g.frontier(bits + (F,))
    print(f"{lab}: f = min({describe(fronts[lab])})")
gap = min(min_gap(fronts['FFF'], fronts[k], 0, 1) for k in fronts if k != 'FFF')
print(f"min over D3 in (0,1) of f(FFF) - f(satisfied) = {gap[0]:.6f}")

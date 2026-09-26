"""Do endpoint swaps + conflict-pair rotations always reach a conflict-free state?

State = (max cut t, transversal S).  Moves:
  swap  : replace an S-vertex by its M-partner (same cut)
  rotate: for a conflict edge uw (u, w in S), flip both u and w; if the cut stays maximum,
          the new uncut edges replace uu', ww'; new transversal = old one with u,w replaced by
          any choice of endpoints of the new uncut edges at u and w.
We BFS over states reachable from each non-healable max cut (with every transversal) and ask
whether a kappa = 0 state is reachable, and the minimum number of rotations needed.
"""
import itertools
from collections import deque

from heal_core import gen
from heal_moves import all_max_cuts, sat_for_cut


def run():
    stats = {"start_states": 0, "reached": 0, "max_rot": 0}
    for G in gen(21, 300):
        n, E = G.number_of_nodes(), list(G.edges())
        adj = [set(G[v]) for v in range(n)]
        fr, cuts = all_max_cuts(n, E)
        for t0 in cuts:
            if sat_for_cut(n, E, t0):
                continue
            # BFS over cuts via rotations; a cut is 'done' if healable (swaps = 2-SAT)
            key = lambda t: tuple(t) if t[0] == 0 else tuple(1 - x for x in t)
            seen = {key(t0): 0}; q = deque([t0]); found = None
            while q and found is None:
                t = q.popleft(); d = seen[key(t)]
                M = [(u, v) for u, v in E if t[u] == t[v]]
                Mv = {x for e in M for x in e}
                for u, w in E:
                    if t[u] == t[w] or u not in Mv or w not in Mv:
                        continue          # rotations act on cut edges joining two M-vertices
                    t2 = list(t); t2[u] ^= 1; t2[w] ^= 1
                    if sum(t2[a] == t2[b] for a, b in E) != fr:
                        continue
                    k = key(t2)
                    if k in seen:
                        continue
                    seen[k] = d + 1
                    if sat_for_cut(n, E, t2):
                        found = d + 1; break
                    q.append(t2)
            stats["start_states"] += 1
            if found is not None:
                stats["reached"] += 1; stats["max_rot"] = max(stats["max_rot"], found)
            else:
                print("NOT REACHED from a non-healable max cut; n =", n)
    print(stats)


if __name__ == "__main__":
    run()

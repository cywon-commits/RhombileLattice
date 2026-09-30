"""Randomised greedy on planar cubic foams (T-join fr-oracle), to reach genuinely stuck tight sets and
exercise Healer.exchange_search at sizes beyond brute force.

Random greedy: add a uniformly random augmenting vertex (any v with S+v independent and tight; found
with the oracle) until none exists.  If fr(G-S) >= 1 at that point, S is a stuck tight set; run the
exchange chain on it and record how many chain advances it needed.

    python3 reviews/constructive_random_greedy.py SEED RUNS NMIN NMAX
"""
import sys, os, random, json, collections
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from constructive_heal import Healer, fr_planar
from lattices3 import random_cubic


def main(seed, runs, nmin, nmax):
    rng = random.Random(seed)
    tally = collections.Counter(); advances = collections.Counter()
    for r in range(runs):
        n, E, pos, info = random_cubic(rng.randint(nmin, nmax), rng.randrange(10**6))
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        h = Healer(G, fr_planar)
        S = set()
        while True:
            k, _ = h.fr(S)
            if k == 0:
                tally["reached_fr0"] += 1; break
            cand = [v for v in G if h.free(S, v) and h.fr(S | {v})[0] == k - 1]
            if not cand:
                tally["stuck_k>=1"] += 1
                adv = [0]; orig = h.component_of
                def counted(T, v, orig=orig, adv=adv):
                    adv[0] += 1; return orig(T, v)
                h.component_of = counted
                try:
                    T = h.exchange_search(S)
                    assert len(T) == len(S) and h.tight(T) and h.augmentable(T)
                    advances[adv[0]] += 1
                except RuntimeError as e:
                    tally["RAISED"] += 1
                    print("RAISED", n, sorted(G.edges()), sorted(S), e, flush=True)
                break
            S.add(rng.choice(cand))
    print(json.dumps(dict(seed=seed, runs=runs, **tally, chain_advances=dict(advances))))


if __name__ == "__main__":
    main(*map(int, sys.argv[1:5]))

"""Search for genuinely stuck tight sets (non-augmentable, fr(G-S) >= 1) and run the exchange chain on them.

Characterisation used (proved in reviews/constructive_review.md, Sec. 2): a tight independent S is
non-augmentable  <=>  every non-bipartite component of G-S is an odd cycle each of whose vertices has
a neighbour in S.  Then fr(G-S) = number of such cycles, so tightness is |S| + #cycles = fr(G) and only
fr(G) itself needs an oracle call.  Every hit is then fed to constructive_heal.Healer.exchange_search.

    python3 reviews/constructive_stuck_search.py geng "-c -d3 -D3" 16      # all cubic graphs on 16 vertices
    python3 reviews/constructive_stuck_search.py k4                        # the excluded case: must raise
"""
import subprocess, sys, os
import numpy as np
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from constructive_heal import Healer, fr_brute


def fr_of(G):
    return fr_brute(G)[0]


def stuck_tight_sets(G, fr):
    nodes = list(G); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    nb = [0] * n
    for a, b in G.edges():
        nb[idx[a]] |= 1 << idx[b]; nb[idx[b]] |= 1 << idx[a]
    full = (1 << n) - 1
    out = []

    def analyse(S):
        rest = full & ~S; k = 0
        while rest:
            start = rest & -rest; comp = start; frontier = start
            while frontier:
                v = (frontier & -frontier).bit_length() - 1; frontier &= frontier - 1
                new = nb[v] & ~S & ~comp; comp |= new; frontier |= new
            rest &= ~comp
            verts = [i for i in range(n) if comp >> i & 1]
            degs = [bin(nb[i] & comp).count("1") for i in verts]
            if all(d == 2 for d in degs):                     # a cycle
                if len(verts) % 2 == 1:
                    if not all(nb[i] & S for i in verts):
                        return None
                    k += 1
                continue
            # not a cycle: must be bipartite
            col = {}; ok = True
            for s in verts:
                if s in col: continue
                col[s] = 0; st = [s]
                while st and ok:
                    v = st.pop()
                    m = nb[v] & comp
                    while m:
                        u = (m & -m).bit_length() - 1; m &= m - 1
                        if u not in col: col[u] = 1 - col[v]; st.append(u)
                        elif col[u] == col[v]: ok = False; break
                if not ok: return None
        return k

    def rec(i, S, size):
        if i == n:
            if size >= fr: return
            k = analyse(S)
            if k is not None and k >= 1 and size + k == fr:
                out.append({nodes[j] for j in range(n) if S >> j & 1})
            return
        rec(i + 1, S, size)
        if size < fr and not (nb[i] & S):
            rec(i + 1, S | 1 << i, size + 1)
    rec(0, 0, 0)
    return out


def exercise(G, S):
    h = Healer(G, fr_brute)
    adv = [0]; orig = h.component_of
    def counted(T, v):                       # called once per chain advance
        adv[0] += 1; return orig(T, v)
    h.component_of = counted
    try:
        T = h.exchange_search(S)
        return "ok", len(T), dict(h.stats, chain_advances=adv[0])
    except RuntimeError as e:
        return "RAISED: " + str(e), None, h.stats


if __name__ == "__main__":
    if sys.argv[1] == "k4":
        G = nx.complete_graph(4)
        print("K4 stuck sets:", stuck_tight_sets(G, fr_of(G)))
        print(exercise(G, {0}))
        G = nx.disjoint_union(nx.complete_graph(4), nx.petersen_graph())
        print("K4+Petersen, S={0}:", exercise(G, {0}))
        sys.exit()
    flags, n = sys.argv[2], sys.argv[3]
    lines = subprocess.run(["nauty-geng", "-q"] + flags.split() + [n], capture_output=True, text=True).stdout.split()
    graphs = hits = 0
    for g6 in lines:
        G = nx.from_graph6_bytes(g6.encode()); graphs += 1
        if G.number_of_nodes() == 4 and G.number_of_edges() == 6:
            continue
        fr = fr_of(G)
        for S in stuck_tight_sets(G, fr):
            hits += 1
            print("HIT", g6, sorted(S), exercise(G, S), flush=True)
    print(dict(n=n, flags=flags, graphs=graphs, stuck_tight_with_k_ge_1=hits))

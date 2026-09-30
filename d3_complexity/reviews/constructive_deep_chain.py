"""Exercise Healer.exchange_search on genuinely stuck tight sets whose chain must ADVANCE (stage >= 1).

Construction (small 'tree of closers', n <= NMAX): an odd cycle Z; every vertex that still has free degree gets
(with prob. q) a fresh closer s in S joined to it and to both ends of a fresh even path Q; S = set of closers.
Then Z and the |S| cycles s+Q are pairwise edge-disjoint odd cycles, so fr(G) >= |S| + 1, while
phi(S) = |S| + fr(G-S).  Hence if G-S has exactly one frustrated component and it is an odd cycle all of whose
vertices have an S-neighbour, S is tight AND stuck (characterisation in constructive_review.md Sec. 2).
Tightness is re-verified with the fr-oracle.  Leaf paths left without closers make exchanges deep in
the tree augmentable, so the chain has to walk down the tree.

    python3 reviews/constructive_deep_chain.py SEED NGRAPHS NMAX   # NMAX <= 20: brute-force oracle; else planar T-join
"""
import sys, os, random, collections, json
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from constructive_heal import Healer, fr_brute, fr_planar


def build(rng, nmax, planar):
    G = nx.Graph(); S = []; cnt = [0]
    def new():
        cnt[0] += 1; G.add_node(cnt[0] - 1); return cnt[0] - 1
    L = rng.choice([3, 3, 5]); Z = [new() for _ in range(L)]
    for i in range(L): G.add_edge(Z[i], Z[(i + 1) % L])
    queue = [(z, 0) for z in Z]; q = rng.uniform(0.3, 0.9); D = rng.randint(0, 3)
    while queue:
        x, d = queue.pop(0)                                  # breadth first: shallow closers first
        if G.degree(x) >= 3: continue
        Lq = rng.choice([2, 2, 4])
        if cnt[0] + 1 + Lq <= nmax and (d <= D or rng.random() < q):
            s = new(); S.append(s); G.add_edge(s, x)
            Q = [new() for _ in range(Lq)]
            nx.add_path(G, Q); G.add_edge(s, Q[0]); G.add_edge(s, Q[-1])
            queue += [(y, d + 1) for y in Q]
    nonS = [v for v in G if v not in S]
    for _ in range(rng.randint(0, 3)):                      # a few random extra edges
        a, b = rng.choice(nonS), rng.choice(nonS + S)
        if a != b and G.degree(a) < 3 and G.degree(b) < 3 and not G.has_edge(a, b):
            G.add_edge(a, b)
            if planar and not nx.check_planarity(G)[0]: G.remove_edge(a, b)
    return G, set(S)


def stuck_structure(G, S):
    H = G.copy(); H.remove_nodes_from(S); k = 0
    for c in nx.connected_components(H):
        C = H.subgraph(c)
        if nx.is_bipartite(C): continue
        if not (all(d == 2 for _, d in C.degree()) and all(any(u in S for u in G[v]) for v in c)):
            return None
        k += 1
    return k


def main(seed, ngraphs, nmax):
    planar = nmax > 20; oracle = fr_planar if planar else fr_brute
    rng = random.Random(seed); tally = collections.Counter(); adv_hist = collections.Counter()
    for _ in range(ngraphs):
        G, S = build(rng, nmax, planar)
        if any(G.has_edge(a, b) for a in S for b in S): continue
        if stuck_structure(G, S) != 1: tally["not_stuck_or_k!=1"] += 1; continue
        h = Healer(G, oracle)
        if not h.tight(S): tally["NOT_TIGHT(unexpected)"] += 1; continue
        tally["stuck_tight_instances"] += 1
        adv = [0]; orig = h.component_of
        def counted(T, v, orig=orig, adv=adv):
            adv[0] += 1; return orig(T, v)
        h.component_of = counted
        try:
            T = h.exchange_search(S)
            assert len(T) == len(S) and h.tight(T) and h.augmentable(T)
            adv_hist[adv[0]] += 1
            c = nx.number_connected_components(G.subgraph(set(G) - S))
            if adv[0] > c: tally["advances>c"] += 1
        except RuntimeError as e:
            tally["RAISED"] += 1; print("RAISED", sorted(G.edges()), sorted(S), e, flush=True)
    print(json.dumps(dict(seed=seed, **tally, chain_advances=dict(sorted(adv_hist.items())))))


if __name__ == "__main__":
    main(*map(int, sys.argv[1:4]))

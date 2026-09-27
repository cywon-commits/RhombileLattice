"""Greedy augmentation for the healing theorem (C.5b): start S = {} and repeatedly add an
endpoint of a frustrated edge of SOME maximum cut of G - S that is not adjacent to S.
Lemma 1 of C.5b says tightness is preserved. Question: does greedy ever get stuck with
fr(G - S) > 0 on a non-K4 graph?  (The proof only covers maximum-size tight sets.)
Max cuts are enumerated exactly with numpy (n <= ~22)."""
import sys, random
import numpy as np
import networkx as nx

def all_maxcut_frustrated_endpoints(H):
    nodes = list(H); n = len(nodes)
    if H.number_of_edges() == 0:
        return 0, set()
    idx = {v: i for i, v in enumerate(nodes)}
    E = np.array([(idx[a], idx[b]) for a, b in H.edges()])
    m = np.arange(1 << (n - 1), dtype=np.int64)
    bits = ((m[:, None] >> np.arange(n)) & 1).astype(np.int8)
    mono = (bits[:, E[:, 0]] == bits[:, E[:, 1]])
    f = mono.sum(1); k = int(f.min())
    if k == 0:
        return 0, set()
    rows = np.where(f == k)[0]
    ends = set()
    used = mono[rows].any(0)
    for e in np.where(used)[0]:
        ends.add(nodes[E[e, 0]]); ends.add(nodes[E[e, 1]])
    return k, ends

def greedy(G, rng):
    S = set()
    while True:
        H = G.copy(); H.remove_nodes_from(S)
        k, ends = all_maxcut_frustrated_endpoints(H)
        if k == 0:
            return S, 0
        cand = [v for v in ends if not any(w in S for w in G[v])]
        if not cand:
            return S, k
        S.add(rng.choice(sorted(cand)))

def rand_graph(n, rng, cubic):
    if cubic and n % 2 == 0:
        return nx.random_regular_graph(3, n, seed=rng.randrange(10**9))
    G = nx.Graph(); G.add_nodes_from(range(n))
    for _ in range(rng.randint(0, n // 3)):
        a, b, c = rng.sample(range(n), 3)
        if all(G.degree(x) <= 1 for x in (a, b, c)):
            G.add_edges_from([(a, b), (b, c), (a, c)])
    for _ in range(4 * n):
        a, b = rng.sample(range(n), 2)
        if G.degree(a) < 3 and G.degree(b) < 3:
            G.add_edge(a, b)
    return G

if __name__ == "__main__":
    seed, count, nmin, nmax = map(int, sys.argv[1:5])
    rng = random.Random(seed); stuck = 0; runs = 0
    for t in range(count):
        n = rng.randint(nmin, nmax)
        G = rand_graph(n, rng, cubic=(t % 2 == 0))
        if not nx.is_connected(G) or nx.is_isomorphic(G, nx.complete_graph(4)):
            continue
        for rep in range(3):
            S, k = greedy(G, random.Random(rng.random()))
            runs += 1
            if k:
                stuck += 1
                print("STUCK", n, sorted(G.edges()), sorted(S), k, flush=True)
        if t % 20 == 0:
            print("progress", t, "runs", runs, "stuck", stuck, flush=True)
    print("done runs", runs, "stuck", stuck)

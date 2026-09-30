"""Referee check of a strengthening of theorem_C_constructive.md: ONE optimal 2-colouring of G suffices.

Given an optimal colouring of G - S, the proof of Lemma stuck is itself an algorithm: from a monochromatic
edge uv, walk v -> v1 -> ... recolouring as in the proof. Either a vertex with no S-neighbour is met
(it is an endpoint of a monochromatic edge of an optimal colouring => augmentation by Lemma aug, and deleting
it leaves an optimal colouring of G - S - v, by (F1)), or the walk closes and certifies that the component is an
odd cycle whose vertices all have an S-neighbour. After an exchange T' = T - t + y an optimal colouring of
G - T' is written down explicitly (Y - y recoloured as a path, t placed opposite the majority of its <= 2
neighbours). So after one ground state of G no further fr-oracle call is needed.

Choices (which monochromatic edge, which p_{j+1}) are randomised so that stuck sets are actually reached.

    python3 reviews/constructive_single_groundstate.py geng "-c -d3 -D3" 12 SEEDS
    python3 reviews/constructive_single_groundstate.py lattices SEEDS
    python3 reviews/constructive_single_groundstate.py deep SEED NGRAPHS NMAX   # start from stuck tree-of-closers sets
"""
import sys, os, random, subprocess, json, collections
import numpy as np
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def ground_states_brute(G):
    nodes = list(G); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    E = np.array([(idx[a], idx[b]) for a, b in G.edges()])
    m = np.arange(1 << n, dtype=np.int64)
    bits = ((m[:, None] >> np.arange(n)) & 1).astype(np.int8)
    f = (bits[:, E[:, 0]] == bits[:, E[:, 1]]).sum(1)
    rows = np.where(f == f.min())[0]
    return int(f.min()), [{nodes[i]: int(bits[r, i]) for i in range(n)} for r in rows]


def mono(G, S, col):
    return [(a, b) for a, b in G.edges() if a not in S and b not in S and col[a] == col[b]]


def walk(G, S, col, u, v):
    """Lemma stuck as an algorithm. Returns ('aug', x) or ('cycle', vertex set). Mutates col (stays optimal)."""
    free = lambda x: not any(w in S for w in G[x])
    if free(u): return "aug", u
    prev, a, seen = u, v, [u]
    while True:
        if free(a): return "aug", a
        others = [b for b in G[a] if b not in S and b != prev]
        assert len(others) == 1, "optimality violated (deg_K 1) or deg_K > 2"
        b = others[0]
        assert col[b] != col[a], "optimality violated (recolouring gains 2)"
        col[a] ^= 1; seen.append(a)                         # now a-b is the monochromatic edge
        if b == u: return "cycle", set(seen)
        prev, a = a, b


def find_aug(G, S, col):
    """('aug', v) or ('stuck', list of odd-cycle components). Mutates col."""
    cert = set(); cycles = []
    while True:
        m = [e for e in mono(G, S, col) if e[0] not in cert]
        if not m: return "stuck", cycles
        kind, x = walk(G, S, col, *m[0])
        if kind == "aug": return "aug", x
        cert |= x; cycles.append(x)


def exchange(G, T, col, Y, y):
    t = next(w for w in G[y] if w in T)
    T2 = (T - {t}) | {y}; c2 = dict(col); del c2[y]
    path = [x for x in nx.dfs_preorder_nodes(G.subgraph(Y - {y}), next(x for x in Y if x != y and
            sum(1 for z in G[x] if z in Y and z != y) == 1))]
    for i, x in enumerate(path): c2[x] = i % 2
    nb = [c2[w] for w in G[t] if w not in T2]
    c2[t] = 0 if sum(nb) * 2 > len(nb) else 1
    return T2, c2, t


def colour_stuck(G, S):
    """Optimal colouring of G - S when every frustrated component is an odd cycle (e.g. S stuck)."""
    H = G.copy(); H.remove_nodes_from(S); col = {}
    for c in nx.connected_components(H):
        C = H.subgraph(c)
        if nx.is_bipartite(C): col.update(nx.bipartite.color(C)); continue
        assert all(d == 2 for _, d in C.degree())
        for i, x in enumerate(nx.dfs_preorder_nodes(C)): col[x] = i % 2
    return col


def heal(G, col, rng, stats, S0=()):
    S = set(S0); col = dict(col)
    fr0 = len(S) + len(mono(G, S, col))
    while True:
        m = mono(G, S, col)
        if not m: break
        rng.shuffle(m)                                       # random choice of starting edge
        kind, x = walk(G, S, col, *m[0])
        if kind == "aug":
            S.add(x); del col[x]; continue
        kind, x = find_aug(G, S, col)
        if kind == "aug":
            S.add(x); del col[x]; continue
        # stuck: exchange chain, all exchanges checked structurally, no oracle
        stats["stuck"] += 1
        H = G.copy(); H.remove_nodes_from(S); c = nx.number_connected_components(H)
        T, tc, Y, w, found = S, col, rng.choice(x), None, None
        for stage in range(c + 1):
            cands = []
            for y in sorted(Y):
                T2, c2, t = exchange(G, T, tc, Y, y)
                assert len(mono(G, T2, c2)) == len(mono(G, S, col)), "exchange colouring not optimal"
                kind, v = find_aug(G, T2, c2)
                if kind == "aug":
                    found = (T2, c2, v); break
                if y != w: cands.append((y, T2, c2, t))
            if found: break
            y, T, tc, w = rng.choice(cands)
            Hc = G.copy(); Hc.remove_nodes_from(T); Y = set(nx.node_connected_component(Hc, w))
            stats["advances"] += 1
        assert found, "chain exceeded c stages"
        stats["max_stage"] = max(stats["max_stage"], stage)
        S, col, v = found
        S = set(S); S.add(v); del col[v]
    assert all(not G.has_edge(a, b) for a in S for b in S)
    H = G.copy(); H.remove_nodes_from(S); assert nx.is_bipartite(H)
    assert len(S) == fr0
    return S


if __name__ == "__main__":
    stats = collections.Counter()
    if sys.argv[1] == "geng":
        flags, n, seeds = sys.argv[2], sys.argv[3], int(sys.argv[4])
        for g6 in subprocess.run(["nauty-geng", "-q"] + flags.split() + [n], capture_output=True, text=True).stdout.split():
            G = nx.from_graph6_bytes(g6.encode())
            if G.number_of_nodes() == 4 and G.number_of_edges() == 6: continue
            fr, gss = ground_states_brute(G); stats["graphs"] += 1
            for s in range(seeds):
                rng = random.Random(s); heal(G, rng.choice(gss), rng, stats); stats["runs"] += 1
    elif sys.argv[1] == "deep":
        from constructive_deep_chain import build, stuck_structure
        seed, ng, nmax = map(int, sys.argv[2:5]); rng0 = random.Random(seed)
        for _ in range(ng):
            G, S = build(rng0, nmax, True)
            if any(G.has_edge(a, b) for a in S for b in S) or stuck_structure(G, S) != 1: continue
            stats["stuck_start"] += 1
            heal(G, colour_stuck(G, S), random.Random(rng0.random()), stats, S0=S)
    else:
        from lattices3 import truncated_penrose, random_cubic, defect_honeycomb
        from ising_tjoin import ising_ground_state_tjoin
        seeds = int(sys.argv[2])
        for name, fn in (("voronoi_foam", lambda: random_cubic(45, 5)),
                         ("defect_honeycomb", lambda: defect_honeycomb(9, 0.12, seed=2)),
                         ("truncated_penrose", lambda: truncated_penrose(3.2, seed=0))):
            nv, E, pos, info = fn(); G = nx.Graph(); G.add_nodes_from(range(nv)); G.add_edges_from(E)
            f, c = ising_ground_state_tjoin(nv, E)
            for s in range(seeds):
                heal(G, {i: int(c[i]) for i in range(nv)}, random.Random(s), stats); stats["runs"] += 1
            print(name, nv, f, flush=True)
    print(json.dumps(dict(sys.argv[1:4] and {"args": " ".join(sys.argv[1:])}, **stats)))

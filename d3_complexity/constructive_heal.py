"""Constructive healing (theorem_C_constructive.md): build an independent odd cycle transversal
of size fr(G) for a graph of maximum degree <= 3 without K4 components, using only an fr-oracle.

    greedy augmentation            S -> S + v   whenever S + v is independent and tight
    exchange search when stuck     T -> T - t + y  (y has exactly one T-neighbour t, result tight)

Pure greedy augmentation can get stuck (smallest cases: 12 vertices, e.g. graph6 K?`CPagPagEG,
planar), so the exchange step is necessary. The proof (theorem_C_constructive.md) shows that from a stuck tight set with fr(G-S) >= 1 a tight
set of the same size that admits an augmentation is reached by a chain of at most
(#components of G - S) exchanges, plus single auxiliary exchanges. So the loop terminates after
fr(G) augmentations. A failure of the exchange search raises an error (it would contradict the proof).

fr-oracles: brute force over 2-colourings (n <= ~22) or the planar T-join (any size, planar only).

Run from d3_complexity/:
    python3 constructive_heal.py random SEED COUNT NMIN NMAX     # brute-force oracle
    python3 constructive_heal.py lattices                        # planar lattices, T-join oracle
"""
import itertools, json, random, sys, time
import numpy as np
import networkx as nx


# ---------------------------------------------------------------- fr oracles
def fr_brute(H):
    """(fr, set of endpoints of monochromatic edges over ALL optimal colourings)."""
    comps = [H.subgraph(c) for c in nx.connected_components(H)]
    tot, ends = 0, set()
    for C in comps:
        if C.number_of_edges() < C.number_of_nodes():   # tree
            continue
        if nx.is_bipartite(C):
            continue
        nodes = list(C); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
        E = np.array([(idx[a], idx[b]) for a, b in C.edges()])
        m = np.arange(1 << (n - 1), dtype=np.int64)
        bits = ((m[:, None] >> np.arange(n)) & 1).astype(np.int8)
        mono = bits[:, E[:, 0]] == bits[:, E[:, 1]]
        f = mono.sum(1); k = int(f.min()); tot += k
        used = mono[f == k].any(0)
        for e in np.where(used)[0]:
            ends.add(nodes[E[e, 0]]); ends.add(nodes[E[e, 1]])
    return tot, ends


def fr_planar(H):
    """(fr, endpoints of monochromatic edges of ONE optimal colouring) via planar T-join."""
    from ising_tjoin import ising_ground_state_tjoin
    tot, ends = 0, set()
    for c in nx.connected_components(H):
        C = H.subgraph(c)
        if C.number_of_edges() < C.number_of_nodes() or nx.is_bipartite(C):
            continue
        nodes = list(C); idx = {v: i for i, v in enumerate(nodes)}
        f, col = ising_ground_state_tjoin(len(nodes), [(idx[a], idx[b]) for a, b in C.edges()])
        tot += f
        for a, b in C.edges():
            if col[idx[a]] == col[idx[b]]:
                ends.update((a, b))
    return tot, ends


# ---------------------------------------------------------------- the algorithm
class Healer:
    def __init__(self, G, oracle):
        self.G, self.oracle = G, oracle
        self.calls = 0
        self.fr0 = self.fr(set())[0]
        self.stats = dict(augment_free=0, augment_oracle=0, exchange_searches=0, exchanges=0)

    def fr(self, S):
        self.calls += 1
        H = self.G.copy(); H.remove_nodes_from(S)
        return self.oracle(H)

    def tight(self, S):
        return len(S) + self.fr(S)[0] == self.fr0

    def free(self, S, v):
        return v not in S and not any(w in S for w in self.G[v])

    def augmentation(self, S, k, ends, count=True):
        """A vertex v with S+v independent and tight, or None."""
        for v in sorted(ends, key=str):                  # Lemma aug: no oracle call needed
            if self.free(S, v):
                self.stats["augment_free"] += count
                return v
        for v in sorted(self.G, key=str):                # other optimal colourings: ask the oracle
            if v not in ends and self.free(S, v) and self.fr(S | {v})[0] == k - 1:
                self.stats["augment_oracle"] += count
                return v
        return None

    def component_of(self, T, v):
        H = self.G.copy(); H.remove_nodes_from(T)
        return set(nx.node_connected_component(H, v))

    def augmentable(self, T):
        k, ends = self.fr(T)
        return self.augmentation(T, k, ends, count=False) is not None

    def exchange_search(self, S):
        """Chain search of theorem_C_constructive.md. S is tight, non-augmentable, fr(G-S) >= 1.
        At every stage j, try the exchange at every vertex of the odd cycle Y_j (this covers the
        auxiliary exchanges used in the proof); if none gives an augmentable set, advance the chain.
        Polynomially many exchanges; returns a tight set of size |S| that admits an augmentation."""
        self.stats["exchange_searches"] += 1
        T = set(S)
        H = self.G.copy(); H.remove_nodes_from(T)
        c = nx.number_connected_components(H)
        Y = next(set(C) for C in nx.connected_components(H)
                 if not nx.is_bipartite(H.subgraph(C)))          # an odd-cycle component (Lemma stuck)
        w = None
        for _ in range(c + 2):
            p_next = None
            for y in sorted(Y, key=str):
                nb = [x for x in self.G[y] if x in T]
                if len(nb) != 1:
                    raise RuntimeError("vertex of Y without a unique T-neighbour: contradicts Lemma stuck")
                T2 = (T - {nb[0]}) | {y}
                self.stats["exchanges"] += 1
                if not self.tight(T2):
                    raise RuntimeError("exchange not tight: contradicts Lemma exchange")
                if self.augmentable(T2):
                    return T2
                if y != w and p_next is None:
                    p_next = (y, nb[0], T2)
            y, t, T = p_next                                      # advance the chain
            Y, w = self.component_of(T, t), t                     # Y_{j+1} = R_{j+1} + w_{j+1}
        raise RuntimeError("chain longer than #components: contradicts Lemma chain (I2)")

    def run(self):
        S = set()
        while True:
            k, ends = self.fr(S)
            if k == 0:
                return S
            v = self.augmentation(S, k, ends)
            if v is None:
                S = self.exchange_search(S)
                k, ends = self.fr(S)
                v = self.augmentation(S, k, ends)
            S.add(v)


def check(G, S, fr0):
    assert all(not G.has_edge(a, b) for a, b in itertools.combinations(S, 2)), "not independent"
    H = G.copy(); H.remove_nodes_from(S)
    assert nx.is_bipartite(H), "not a transversal"
    assert len(S) == fr0, "wrong size"


def rand_subcubic(n, rng):
    G = nx.Graph(); G.add_nodes_from(range(n))
    for _ in range(rng.randint(0, n // 3)):          # seed triangles
        a, b, c = rng.sample(range(n), 3)
        if all(G.degree(x) == 0 for x in (a, b, c)):
            G.add_edges_from([(a, b), (b, c), (a, c)])
    for _ in range(4 * n):
        a, b = rng.sample(range(n), 2)
        if G.degree(a) < 3 and G.degree(b) < 3:
            G.add_edge(a, b)
    return G


def has_K4_component(G):
    return any(len(c) == 4 and G.subgraph(c).number_of_edges() == 6 for c in nx.connected_components(G))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "random":
        seed, count, nmin, nmax = map(int, sys.argv[2:6])
        rng = random.Random(seed); tot = dict(graphs=0, augment_free=0, augment_oracle=0,
                                              exchange_searches=0, exchanges=0)
        for _ in range(count):
            G = rand_subcubic(rng.randint(nmin, nmax), rng)
            if has_K4_component(G):
                continue
            h = Healer(G, fr_brute); S = h.run(); check(G, S, h.fr0)
            tot["graphs"] += 1
            for key, val in h.stats.items():
                tot[key] += val
        print(json.dumps(tot))
    elif mode == "lattices":
        from lattices3 import truncated_penrose, random_cubic, defect_honeycomb
        for name, fn in (("truncated_penrose", lambda: truncated_penrose(3.2, seed=0)),
                         ("voronoi_foam", lambda: random_cubic(45, 5)),
                         ("defect_honeycomb", lambda: defect_honeycomb(9, 0.12, seed=2)),
                         ("truncated_penrose_4", lambda: truncated_penrose(4, seed=0))):
            n, E, pos, info = fn()
            G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
            t0 = time.time(); h = Healer(G, fr_planar); S = h.run(); check(G, S, h.fr0)
            print(json.dumps(dict(lattice=name, n=n, fr=h.fr0, oracle_calls=h.calls,
                                  seconds=round(time.time() - t0, 2), **h.stats)))

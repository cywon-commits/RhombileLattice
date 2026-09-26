"""Find non-healable maximum cuts in dumbbell-rich subcubic graphs and describe the 2-SAT core."""
import itertools, random, collections
import networkx as nx
from heal_moves import all_max_cuts, sat_for_cut


def implication_core(n, E, t):
    adj = [set() for _ in range(n)]
    for u, v in E:
        adj[u].add(v); adj[v].add(u)
    M = [(u, v) for u, v in E if t[u] == t[v]]
    partner = {}
    for u, v in M:
        partner[u] = v; partner[v] = u
    # arc x -> partner(y) for every conflict edge x~y between M-vertices of different M-edges
    D = nx.DiGraph()
    for x in partner:
        D.add_node(x)
        for y in adj[x]:
            if y in partner and partner[x] != y:
                D.add_edge(x, partner[y])
    bad = []
    for u, v in M:
        if nx.has_path(D, u, v) and nx.has_path(D, v, u):
            bad.append((u, v, nx.shortest_path(D, u, v), nx.shortest_path(D, v, u)))
    return M, bad


def gen(seed, count):
    random.seed(seed)
    def dumbbell(off):
        a, y, z, b, c, d = [off + i for i in range(6)]
        return [(a, y), (a, z), (y, z), (b, c), (b, d), (c, d), (a, b)], [y, z, c, d]
    for k in range(count):
        ndb = random.choice([1, 2])
        E = []; stubs = []; off = 0
        for j in range(ndb):
            e, p = dumbbell(off); E += e; stubs += p; off += 6
        n = random.choice([6, 8, 10, 12]); H = nx.random_regular_graph(3, n, seed=random.randrange(10**9))
        M = list(nx.maximal_matching(H))[:len(stubs) // 2]
        if len(M) < len(stubs) // 2:
            continue
        H.remove_edges_from(M); ends = [v + off for e in M for v in e]; random.shuffle(ends)
        E = E + [(u + off, v + off) for u, v in H.edges()] + list(zip(stubs, ends))
        G = nx.Graph(E)
        if not nx.is_connected(G) or max(dict(G.degree()).values()) > 3 or G.number_of_nodes() > 20:
            continue
        yield nx.convert_node_labels_to_integers(G)


if __name__ == "__main__":
    shapes = collections.Counter(); total = 0; nohealable = 0
    for G in gen(21, 300):
        n, E = G.number_of_nodes(), list(G.edges())
        fr, cuts = all_max_cuts(n, E)
        ok = [t for t in cuts if sat_for_cut(n, E, t)]
        if not ok:
            nohealable += 1
        for t in cuts:
            if t in ok:
                continue
            total += 1
            M, bad = implication_core(n, E, t)
            for u, v, p1, p2 in bad:
                shapes[(len(p1) - 1, len(p2) - 1)] += 1
    print("non-healable max cuts:", total, " graphs with no healable max cut:", nohealable)
    print("implication path lengths (a=>b, b=>a) at failing M-edges:", dict(shapes))

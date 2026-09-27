"""Second referee pass on Theorem C.5: randomized test of the purely combinatorial core of the chain lemma
on graphs built so that the local hypotheses H_j hold for many steps.

H(T,Y): Y is an odd-cycle component of G-T; every y in Y has degree 3 and exactly one T-neighbour t_y;
        the t_y are pairwise distinct; each t_y has degree 3, y is its only neighbour on Y, and the other two
        neighbours of t_y are the two endpoints of one component Q_y of G-T which is a path on an even
        number of vertices.   (= Lemma stuck + Lemma exchange (a),(b) at stage (T,Y).)
The proof of the chain lemma derives (I1)-(I3) at stage j+1 from H(T_0,Y_0),...,H(T_j,Y_j) alone
(plus independence of S).  We build random graphs rich in such configurations (odd cycles, even paths, each path
with a private 'closer' in S whose third edge goes somewhere random, extra random edges, deliberate attempts to make
closers hit rotated paths / Z-p1), and check every chain step for which H holds at all previous stages.
Usage: python3 c5_review2_structural.py SEED NGRAPHS
"""
import sys, random, collections

seed, NG = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
viol = collections.Counter(); stats = collections.Counter(); ex = {}

def bits(x):
    i = 0
    while x:
        if x & 1: yield i
        x >>= 1; i += 1

def build():
    """tree-of-closers construction: every vertex that needs an S-neighbour gets (mostly) a fresh closer s
    whose other two edges close a fresh even path; with some probability it instead reuses an existing
    S-vertex with free degree (this is how coincidences with Z-p1 / rotated paths could arise)."""
    E = set(); nxt = [0]; deg = collections.Counter(); Sset = []
    def new():
        nxt[0] += 1; return nxt[0] - 1
    def edge(u, v):
        E.add(tuple(sorted((u, v)))); deg[u] += 1; deg[v] += 1
    L = rng.choice([3, 3, 5, 5, 7]); Z = [new() for _ in range(L)]
    for i in range(L): edge(Z[i], Z[(i + 1) % L])
    queue = list(Z); cap = rng.randint(25, 90); q = rng.uniform(0.6, 0.97)
    while queue:
        x = queue.pop(rng.randrange(len(queue))) if rng.random() < 0.5 else queue.pop(0)
        if deg[x] >= 3: continue
        reuse = [s for s in Sset if deg[s] < 3 and not any(tuple(sorted((s, y))) in E for y in [x])]
        if nxt[0] < cap and (rng.random() < q or not reuse):
            s = new(); Sset.append(s); edge(s, x)
            Lq = rng.choice([2, 2, 2, 4, 4, 6]); Q = [new() for _ in range(Lq)]
            for i in range(Lq - 1): edge(Q[i], Q[i + 1])
            edge(s, Q[0]); edge(s, Q[-1])
            queue += Q
        elif reuse:
            edge(rng.choice(reuse), x)
    # leftover: random extra edges non-S -- S (keeps S independent) and occasionally non-S -- non-S
    nonS = [v for v in range(nxt[0]) if v not in set(Sset)]
    for _ in range(rng.randint(0, 6)):
        a = rng.choice(nonS); b = rng.choice(nonS if rng.random() < 0.3 else Sset)
        if a != b and deg[a] < 3 and deg[b] < 3 and tuple(sorted((a, b))) not in E: edge(a, b)
    n = nxt[0]; adj = [0] * n
    for u, v in E: adj[u] |= 1 << v; adj[v] |= 1 << u
    return n, adj, sum(1 << s for s in Sset), sorted(E)

def run(n, adj, S, E):
    full = (1 << n) - 1; deg = [bin(a).count("1") for a in adj]
    assert max(deg) <= 3
    assert not any(adj[v] & S for v in bits(S))
    def comps(mask):
        res = []; rem = mask
        while rem:
            s = rem & -rem; comp = s; front = s
            while front:
                nb = 0
                for v in bits(front): nb |= adj[v]
                nb &= mask & ~comp; comp |= nb; front = nb
            res.append(comp); rem &= ~comp
        return res
    def is_odd_cycle(c):
        vs = list(bits(c))
        return len(vs) % 2 == 1 and len(vs) >= 3 and all(bin(adj[v] & c).count("1") == 2 for v in vs)
    def H(Tj, Yj):
        csj = comps(full & ~Tj)
        if Yj not in csj or not is_odd_cycle(Yj): return None
        tl = []
        for y in bits(Yj):
            tn = adj[y] & Tj
            if deg[y] != 3 or bin(tn).count("1") != 1: return None
            tl.append(tn)
        if len(set(tl)) != len(tl): return None
        Qs = {}
        for y in bits(Yj):
            t = next(bits(adj[y] & Tj))
            if deg[t] != 3 or (adj[t] & Yj) != (1 << y): return None
            others = adj[t] & ~(1 << y)
            Q = [c for c in csj if c & others]
            if len(Q) != 1: return None
            Q = Q[0]; qv = list(bits(Q)); qd = [bin(adj[v] & Q).count("1") for v in qv]
            if len(qv) % 2 or max(qd) > 2 or sum(qd) != 2 * (len(qv) - 1): return None
            if set(v for v, d in zip(qv, qd) if d <= 1) != set(bits(others)): return None
            Qs[y] = (t, Q)
        return csj, Qs
    def rep(kind, info):
        viol[kind] += 1
        if kind not in ex: ex[kind] = (E, list(bits(S)), info); print("VIOLATION", kind, ex[kind], flush=True)
    csS = comps(full & ~S)
    for Z in csS:
        if not is_odd_cycle(Z): continue
        stack = [(S, [Z], [], [])]; budget = 20000
        while stack and budget > 0:
            budget -= 1
            Tj, R, W, P = stack.pop(); j = len(W)
            Yj = R[-1] | ((1 << W[-1]) if W else 0)
            h = H(Tj, Yj)
            if h is None: stats["end_depth_%d" % j] += 1; continue
            stats["H_states_depth_%d" % min(j, 9)] += 1
            if j > len(csS) + 1: rep("chain_longer_than_#components", (W, P)); continue
            csj, Qs = h
            if j >= 1:
                expect = set(c for c in csS if c not in R)
                expect.add(R[0] & ~(1 << P[0]))
                for i in range(1, j): expect.add((R[i] | (1 << W[i - 1])) & ~(1 << P[i]))
                expect.add(Yj)
                if expect != set(csj): rep("component_list", (W, P))
            for p in bits(R[-1]):
                w, Q = Qs[p]; T2 = (Tj & ~(1 << w)) | (1 << p); stats["steps"] += 1
                if not ((S >> w) & 1) or w in W: rep("w_in_S_minus_W", (W, P, p)); continue
                if Q not in csS: rep("I2_component_of_G-S", (W, P, p)); continue
                if Q in R: rep("I2_novel", (W, P, p)); continue
                Yn = Q | (1 << w)
                if Yn not in comps(full & ~T2) or not is_odd_cycle(Yn): rep("I1_cycle", (W, P, p)); continue
                Rn = R + [Q]
                for c in csS:
                    if c in Rn: continue
                    for v in bits(c):
                        if (adj[v] & S) != (adj[v] & T2): rep("I3_untouched", (W, P, p))
                for v in bits(Q):
                    if (adj[v] & S) != ((adj[v] & T2) | ((1 << w) if (adj[v] >> w) & 1 else 0)): rep("I3_R", (W, P, p))
                stack.append((T2, Rn, W + [w], P + [p]))

for g in range(NG):
    n, adj, S, E = build(); stats["graphs"] += 1; stats["vertices"] += n
    run(n, adj, S, E)
print("STATS", dict(sorted(stats.items())))
print("VIOLATIONS", dict(viol))

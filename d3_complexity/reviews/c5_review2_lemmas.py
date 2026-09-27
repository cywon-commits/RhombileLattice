"""Second referee pass on Theorem C.5 (healing theorem): NON-VACUOUS tests of the intermediate lemmas.

The hypotheses of Lemmas stuck/exchange/chain ("maximum tight set with k>=1") are vacuous if the theorem
is true, so we test RELAXED versions whose hypotheses are satisfiable, and whose proofs in
paper/healing_proof.tex go through verbatim:

 (A) exchange-tightness: T independent tight, Y an odd-cycle component of G-T, y in Y with exactly one
     T-neighbour t  ==>  T' = T - t + y is independent and tight.   (proof of Lemma exchange, 1st para)
 (B) relaxed exchange lemma: call an independent T "blocked" if every frustrated component of G-T is an
     odd cycle whose vertices all have degree 3 in G and exactly one T-neighbour (= conclusion of
     Lemma stuck).  If T and T' are both blocked with fr(G-T') = fr(G-T), and Y u {t} is not a K4, then
     (a),(b),(c) of Lemma exchange hold.   (The proof uses maximality ONLY through these two facts.)
 (C) relaxed chain lemma: run chains from a blocked S; as long as every T_0..T_{j+1} is blocked with
     constant fr, check (I1)-(I3), w_{j+1} in S\\W_j, and the component list (i)-(iii) of G-T_j.
     Record how each chain terminates (hypothesis breaks: fr drops / T' not blocked / K4).
 (D) Lemma stuck under the weaker hypothesis "non-augmentable tight" and lemma aug, on all tight sets.

Graphs: all connected graphs with max degree <= 3 on n <= N vertices (nauty-geng -c -D3).
Usage: python3 c5_review2_lemmas.py N
"""
import sys, subprocess, itertools, collections
import numpy as np
import networkx as nx

N = int(sys.argv[1])
viol = collections.Counter(); stats = collections.Counter(); examples = {}

def report(kind, info):
    viol[kind] += 1
    if kind not in examples:
        examples[kind] = info
        print("VIOLATION", kind, info, flush=True)

_col_cache = {}
def colourings(m):
    if m not in _col_cache:
        _col_cache[m] = ((np.arange(1 << max(m - 1, 0))[:, None] >> np.arange(m)[None, :]) & 1).astype(np.int8)
    return _col_cache[m]

def graphs(n):
    out = subprocess.run(["nauty-geng", "-c", "-D3", "-q", str(n)], capture_output=True, text=True).stdout.split()
    for g6 in out:
        yield nx.from_graph6_bytes(g6.encode())

def bits(x):
    i = 0
    while x:
        if x & 1: yield i
        x >>= 1; i += 1

def run_graph(G):
    n = G.number_of_nodes(); V = list(range(n))
    adj = [0] * n
    for u, v in G.edges(): adj[u] |= 1 << v; adj[v] |= 1 << u
    deg = [bin(a).count("1") for a in adj]
    full = (1 << n) - 1
    frmemo = {}

    def comps(mask):
        res = []; rem = mask
        while rem:
            s = rem & -rem; comp = s; front = s
            while front:
                nb = 0
                for v in bits(front): nb |= adj[v]
                nb &= mask & ~comp
                comp |= nb; front = nb
            res.append(comp); rem &= ~comp
        return res

    def fr_comp(c):
        if c in frmemo: return frmemo[c]
        vs = list(bits(c)); idx = {v: i for i, v in enumerate(vs)}
        es = [(idx[u], idx[v]) for u in vs for v in bits(adj[u] & c) if u < v]
        if not es: frmemo[c] = 0; return 0
        C = colourings(len(vs))
        a = np.array([e[0] for e in es]); b = np.array([e[1] for e in es])
        f = int((C[:, a] == C[:, b]).sum(1).min()); frmemo[c] = f; return f

    def fr(mask):
        return sum(fr_comp(c) for c in comps(mask))

    def is_odd_cycle(c):
        vs = list(bits(c))
        if len(vs) % 2 == 0 or len(vs) < 3: return False
        if any(bin(adj[v] & c).count("1") != 2 for v in vs): return False
        return True  # connected + all degrees 2 => cycle

    def blocked(T):
        """returns (is_blocked, list of odd-cycle comps, list of all comps)"""
        cs = comps(full & ~T); odd = []
        for c in cs:
            if fr_comp(c) == 0: continue
            if not is_odd_cycle(c): return False, None, cs
            for v in bits(c):
                if deg[v] != 3 or bin(adj[v] & T).count("1") != 1: return False, None, cs
            odd.append(c)
        return True, odd, cs

    frG = fr(full)
    # enumerate independent sets
    indep = []
    def rec(i, T, forb):
        if i == n: indep.append(T); return
        rec(i + 1, T, forb)
        if not (forb >> i) & 1: rec(i + 1, T | (1 << i), forb | adj[i] | (1 << i))
    rec(0, 0, 0)
    isK4 = (n == 4 and G.number_of_edges() == 6)

    for T in indep:
        k = fr(full & ~T); tight = (bin(T).count("1") + k == frG)
        if bin(T).count("1") + k < frG: report("F1", (list(G.edges()), list(bits(T))))
        cs = comps(full & ~T)
        # ---------- (A) exchange tightness, (D) aug + stuck
        if tight:
            stats["tight"] += 1
            for c in cs:
                if not is_odd_cycle(c): continue
                for y in bits(c):
                    tn = adj[y] & T
                    if bin(tn).count("1") != 1: continue
                    t = next(bits(tn)); T2 = (T & ~(1 << t)) | (1 << y)
                    stats["A_tested"] += 1
                    if any(adj[v] & T2 for v in bits(T2)): report("A_indep", (list(G.edges()), list(bits(T)), y))
                    elif bin(T2).count("1") + fr(full & ~T2) != frG: report("A_tight", (list(G.edges()), list(bits(T)), y))
            if k >= 1:
                # non-augmentable? (for every optimal colouring, every endpoint of a mono edge has a T-nbr)
                aug = False
                for c in cs:
                    if fr_comp(c) == 0: continue
                    vs = list(bits(c)); idx = {v: i for i, v in enumerate(vs)}
                    es = [(u, v) for u in vs for v in bits(adj[u] & c) if u < v]
                    C = colourings(len(vs)); a = np.array([idx[e[0]] for e in es]); b = np.array([idx[e[1]] for e in es])
                    mono = (C[:, a] == C[:, b]); opt = mono.sum(1) == fr_comp(c)
                    endpoints = set()
                    for row in np.nonzero(opt)[0]:
                        for j in np.nonzero(mono[row])[0]: endpoints |= {es[j][0], es[j][1]}
                    if any(not (adj[v] & T) for v in endpoints): aug = True
                    for v in endpoints:
                        if not (adj[v] & T):
                            T3 = T | (1 << v)
                            if bin(T3).count("1") + fr(full & ~T3) != frG: report("aug_lemma", (list(G.edges()), list(bits(T)), v))
                if not aug:
                    stats["tight_nonaug_k>=1"] += 1
                    ok, _, _ = blocked(T)
                    if not ok: report("stuck_lemma", (list(G.edges()), list(bits(T))))
                    if not isK4: report("tight_stuck_nonK4", (list(G.edges()), list(bits(T))))
        # ---------- (E) purely STRUCTURAL chain core: hypotheses H_0..H_j (conclusions of Lemma stuck +
        # Lemma exchange at each stage, stated locally) ==> conclusions of the chain lemma at j+1.
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
        csS = comps(full & ~T)
        for Z in csS:
            if not is_odd_cycle(Z): continue
            stack = [(T, [Z], [], [])]
            while stack:
                Tj, R, W, P = stack.pop(); j = len(W)
                Yj = R[-1] | ((1 << W[-1]) if W else 0)
                h = H(Tj, Yj)
                if h is None: stats["E_end_depth_%d" % j] += 1; continue
                if j > n: report("E_too_long", (list(G.edges()), list(bits(T)), W, P)); continue
                csj, Qs = h
                stats["E_states_depth_%d" % j] += 1
                if j >= 1:  # component list (i)-(iii) of G - T_j
                    expect = set(c for c in csS if c not in R)
                    expect.add(R[0] & ~(1 << P[0]))
                    for i in range(1, j): expect.add((R[i] | (1 << W[i - 1])) & ~(1 << P[i]))
                    expect.add(Yj)
                    if expect != set(csj): report("E_component_list", (list(G.edges()), list(bits(T)), W, P))
                for p in bits(R[-1]):
                    w, Q = Qs[p]; T2 = (Tj & ~(1 << w)) | (1 << p)
                    stats["E_steps"] += 1
                    if not ((T >> w) & 1) or w in W: report("E_w_in_S_minus_W", (list(G.edges()), list(bits(T)), W, P, p)); continue
                    if Q not in csS: report("E_I2_comp_of_G-S", (list(G.edges()), list(bits(T)), W, P, p)); continue
                    if Q in R: report("E_I2_novel", (list(G.edges()), list(bits(T)), W, P, p)); continue
                    if any(adj[v] & T2 for v in bits(T2)): report("E_indep", 0); continue
                    Yn = Q | (1 << w)
                    if Yn not in comps(full & ~T2) or not is_odd_cycle(Yn): report("E_I1_cycle", (list(G.edges()), list(bits(T)), W, P, p)); continue
                    Rn = R + [Q]
                    for c in csS:
                        if c in Rn: continue
                        for v in bits(c):
                            if (adj[v] & T) != (adj[v] & T2): report("E_I3_untouched", (list(G.edges()), list(bits(T)), W, P, p))
                    for v in bits(Q):
                        exp = (adj[v] & T2) | ((1 << w) if (adj[v] >> w) & 1 else 0)
                        if (adj[v] & T) != exp: report("E_I3_R", (list(G.edges()), list(bits(T)), W, P, p))
                    stack.append((T2, Rn, W + [w], P + [p]))
        # ---------- (B)+(C) relaxed exchange / chain from blocked sets
        ok, odd, _ = blocked(T)
        if not ok or not odd: continue
        stats["blocked_S_k>=1"] += 1
        if tight: stats["blocked_and_tight"] += 1
        S = T; compsS = comps(full & ~S); compid = {}
        for i, c in enumerate(compsS):
            for v in bits(c): compid[v] = i
        def tnb(v, TT): return adj[v] & TT
        for Z in odd:
            # DFS over chains. state: j, T_j, R list (masks, R[0]=Z), W list, P list
            stack = [(S, [Z], [], [])]
            steps = 0
            while stack:
                Tj, R, W, P = stack.pop(); j = len(W)
                Rj = R[-1]; Yj = Rj | ((1 << W[-1]) if W else 0)
                if j > n: report("chain_too_long", (list(G.edges()), list(bits(S)))); continue
                # local hypothesis good(T_j, Y_j): every single exchange at a vertex of Y_j is again
                # blocked with the same fr (exactly the sets on which the proof invokes Lemma stuck)
                good = True
                for y in bits(Yj):
                    tn = tnb(y, Tj)
                    if bin(tn).count("1") != 1: good = False; break
                    t_ = next(bits(tn)); okx, oddx, _ = blocked((Tj & ~(1 << t_)) | (1 << y))
                    if not okx or len(oddx) != len(odd): good = False; break
                if not good:
                    stats["C_end_not_good_at_depth_%d" % j] += 1; continue
                if isK4: stats["C_K4_good"] += 1; continue
                stats["C_good_states"] += 1
                tset = [tnb(y, Tj) for y in bits(Yj)]
                if len(set(tset)) != len(tset): report("B_distinct_Tnbrs", (list(G.edges()), list(bits(S)), W, P))
                for p in bits(Rj):
                    steps += 1
                    tn = tnb(p, Tj)
                    if bin(tn).count("1") != 1: report("C_unique_Tnbr", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    w = next(bits(tn)); T2 = (Tj & ~(1 << w)) | (1 << p)
                    ok2, odd2, cs2 = blocked(T2)
                    k2 = len(odd2) if ok2 else None
                    if not ok2: report("C_good_but_Tprime_not_blocked", 0); continue
                    if k2 != len(odd): report("C_good_but_fr_changed", 0); continue
                    if isK4: stats["C_K4"] += 1; continue
                    stats["C_steps_hyp_ok"] += 1
                    # (B) exchange conclusions
                    others = adj[w] & ~(1 << p)
                    if deg[w] != 3 or (adj[w] & Yj) != (1 << p):
                        report("B_a", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    csj = comps(full & ~Tj)
                    Q = [c for c in csj if c & others]
                    if len(Q) != 1: report("B_b_onecomp", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    Q = Q[0]
                    qv = list(bits(Q)); qd = {v: bin(adj[v] & Q).count("1") for v in qv}
                    ends = [v for v in qv if qd[v] <= 1]
                    ispath = (sum(qd.values()) == 2 * (len(qv) - 1)) and max(qd.values()) <= 2 and len(qv) >= 2
                    if not ispath or len(qv) % 2 or set(ends) != set(bits(others)):
                        report("B_b_path", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    if (Q | (1 << w)) not in odd2: report("B_c", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    # (C) chain invariants
                    if not ((S >> w) & 1) or w in W: report("C_w_in_S_minus_W", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    if Q not in compsS: report("C_I2_component_of_G-S", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    if Q in R: report("C_I2_novel", (list(G.edges()), list(bits(S)), W, P, p)); continue
                    # component list of G - T_j : (i) untouched comps, (ii) Z-p1 and rho_i, (iii) Y_j
                    if j >= 1:
                        expect = set(c for c in compsS if c not in R)
                        expect.add(R[0] & ~(1 << P[0]))
                        for i in range(1, j):
                            expect.add((R[i] | (1 << W[i - 1])) & ~(1 << P[i]))
                        expect.add(Yj)
                        if expect != set(csj): report("C_component_list", (list(G.edges()), list(bits(S)), W, P))
                    # I3 at j+1
                    Rn = R + [Q]; Wn = W + [w]; Pn = P + [p]
                    for c in compsS:
                        if c in Rn: continue
                        for v in bits(c):
                            if (adj[v] & S) != (adj[v] & T2): report("C_I3_untouched", (list(G.edges()), list(bits(S)), Wn, Pn))
                    for v in bits(Q):
                        exp = (adj[v] & T2) | ((1 << w) if (adj[v] >> w) & 1 else 0)
                        if (adj[v] & S) != exp or bin(adj[v] & T2).count("1") != 1: report("C_I3_R", (list(G.edges()), list(bits(S)), Wn, Pn))
                    stats["C_depth_%d" % (j + 1)] += 1
                    stack.append((T2, Rn, Wn, Pn))
            stats["chains_from_Z"] += 1

for n in range(3, N + 1):
    cnt = 0
    for G in graphs(n):
        cnt += 1; run_graph(nx.convert_node_labels_to_integers(G))
    stats["graphs_n%d" % n] = cnt
    print(f"n={n}: {cnt} graphs; violations so far {dict(viol)}", flush=True)
print("STATS", dict(sorted(stats.items())))
print("VIOLATIONS", dict(viol))

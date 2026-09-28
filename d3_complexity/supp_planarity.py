"""Evidence for Supplementary Note S1: the reduction graph G(phi) of Johnson et al.
(Algorithmica 87 (2025), Thm 11, IOCT part) is non-planar on explicit planar 3-SAT instances.

Run from d3_complexity/:   python3 supp_planarity.py [n_random] [seed]
Writes supp_runs/planarity_evidence.json and prints a summary.

Checks, for each formula phi:
  (a) the formula conditions: clause sizes 2-3, distinct variables in a clause, every literal in
      at most two clauses, every variable in at most three clauses (strict variant);
  (b) the literal-split incidence graph I(phi) (variable x -> edge x -- not-x) is planar;
  (c) the minor H(phi) = I(phi) + one apex adjacent to every clause vertex (Lemma S1) is planar?
  (d) G(phi) itself, built exactly as in verify_delta3_reduction.build (canonical clause and
      literal order), is planar?  If not, a Kuratowski subgraph is extracted and classified.
  (e) the reduction equivalence OPT_D3 <= 2m D3 <=> phi satisfiable (exact MILP), for small phi.
Lemma S1 implies (c) non-planar => (d) non-planar for EVERY clause and literal order.
"""
import itertools, json, os, random, sys
import networkx as nx
from verify_delta3_reduction import build, opt, sat

PHI_STAR = [[1, 2], [-1, 2], [1, -2]]   # smallest: minor H(phi*) is exactly K_{3,3}
# from the referee check reviews/supp_s1_review.md: (formula, variable order for Lichtenstein's cycle)
PHI_B2 = ([[-5, -4, 2], [-6, -3, 4], [3, 5, 1], [-6, 3, 1], [4, 6, -2], [-1, -5, 2], [5, -4, -3], [-2, -1, 6]],
          [1, 2, 4, 5, 3, 6])     # exactly 3 literals per clause, every literal exactly twice
PHI_U = ([[-4, -2], [4, 5], [3, 1], [-2, 4, -5], [-5, -3], [2, -1], [1, -3]], [1, 2, 3, 4, 5])  # unsatisfiable
PHI_2U = ([[1, -3], [4, -2], [-2, 1], [-3, 4], [-4, -1], [3, 2]], None)                         # unsatisfiable 2-CNF
PHI0 = [[5, 2], [4, -5], [1, -6], [6, 3], [3, -4], [6, 4, -2], [-2, 1], [-1, 5], [-6, -1]]


def conditions(phi):
    occ_lit, occ_var = {}, {}
    for cl in phi:
        assert 2 <= len(cl) <= 3 and len({abs(l) for l in cl}) == len(cl)
        for l in cl:
            occ_lit[l] = occ_lit.get(l, 0) + 1
            occ_var[abs(l)] = occ_var.get(abs(l), 0) + 1
    return dict(max_lit=max(occ_lit.values()), max_var=max(occ_var.values()),
                sizes=sorted({len(c) for c in phi}))


def incidence(phi, split=True):
    G = nx.Graph()
    for cl in phi:
        for l in cl:
            G.add_edge(("L", l) if split else ("V", abs(l)), ("c", tuple(cl)))
            if split:
                G.add_edge(("L", abs(l)), ("L", -abs(l)))
    return G


def apex_minor(phi):
    H = incidence(phi, split=True)
    for cl in phi:
        H.add_edge(("c", tuple(cl)), "apex")
    return H


def merged_apex_minor(phi):
    """H(phi) as defined in the note: plain incidence graph (one vertex per variable) + apex."""
    H = incidence(phi, split=False)
    for cl in phi:
        H.add_edge(("c", tuple(cl)), "apex")
    return H


def variable_cycle_planar(phi, order):
    """Lichtenstein's condition for a given variable order: incidence + cycle through the variables."""
    if order is None:
        return None
    G = incidence(phi, split=False)
    if len(order) >= 3:
        for a, b in zip(order, order[1:] + order[:1]):
            G.add_edge(("V", a), ("V", b))
    return nx.check_planarity(G)[0]


def kuratowski_type(K):
    """Classify a Kuratowski subgraph: branch vertices = degree >= 3."""
    deg = [d for _, d in K.degree()]
    b = sum(1 for d in deg if d >= 3)
    return {"branch_vertices": b, "type": "K5 subdivision" if b == 5 else "K3,3 subdivision",
            "vertices": K.number_of_nodes(), "edges": K.number_of_edges()}


def analyse(phi, milp=False):
    nv = max(abs(l) for cl in phi for l in cl)
    n, E, idx = build(nv, phi)
    G = nx.Graph(); G.add_edges_from(E)
    assert max(d for _, d in G.degree()) <= 3
    planarG, cert = nx.check_planarity(G, counterexample=True)
    rec = dict(phi=phi, m=len(phi), nv=nv, cond=conditions(phi),
               incidence_planar=nx.check_planarity(incidence(phi))[0],
               apex_minor_planar=nx.check_planarity(apex_minor(phi))[0],
               H_planar=nx.check_planarity(merged_apex_minor(phi))[0],
               plain_incidence_planar=nx.check_planarity(incidence(phi, split=False))[0],
               G_n=n, G_edges=len(E), G_planar=planarG, satisfiable=sat(nv, phi))
    if not planarG:
        inv = {v: k for k, v in idx.items()}
        K = cert
        rec["kuratowski"] = kuratowski_type(K)
        rec["kuratowski_branch"] = [str(inv[v]) for v, d in K.degree() if d >= 3]
    if milp:
        m = len(phi)
        rec["milp"] = {str(D3): opt(n, E, D3) for D3 in (0.3, 0.6, 0.9)}
        rec["milp_equiv_ok"] = all((o <= 2 * m * float(D3) + 1e-7) == rec["satisfiable"]
                                   for D3, o in rec["milp"].items())
    return rec


def random_strict_phi(rng, nv, m):
    """Clause sizes 2-3, each literal <= 2 times, each variable <= 3 times."""
    for _ in range(5000):
        occ_l, occ_v, cls = {}, {}, []
        for _ in range(m):
            k = rng.choice([2, 3])
            free = [v for v in range(1, nv + 1) if occ_v.get(v, 0) < 3]
            if len(free) < k:
                break
            vs = rng.sample(free, k)
            cl = []
            for v in vs:
                opts = [s * v for s in (1, -1) if occ_l.get(s * v, 0) < 2]
                l = rng.choice(opts)
                cl.append(l); occ_l[l] = occ_l.get(l, 0) + 1; occ_v[v] = occ_v.get(v, 0) + 1
            cls.append(cl)
        else:
            if nx.check_planarity(incidence(cls))[0]:
                return cls
    return None


def random_exact3_phi(rng, nv, clause3_only):
    """Every variable in exactly 3 clauses with both polarities (one literal twice, the other once);
    clause sizes 3 only (m = nv) or 2-3. Returns a formula with planar literal-split incidence."""
    for _ in range(20000):
        lits = []
        for v in range(1, nv + 1):
            s = rng.choice([1, -1]); lits += [s * v, s * v, -s * v]
        rng.shuffle(lits)
        if clause3_only:
            if len(lits) % 3:
                continue
            cls = [lits[i:i + 3] for i in range(0, len(lits), 3)]
        else:
            cls, i = [], 0
            while i < len(lits):
                k = rng.choice([2, 3]) if len(lits) - i not in (2, 3, 4) else (2 if len(lits) - i in (2, 4) else 3)
                cls.append(lits[i:i + k]); i += k
        if any(len({abs(l) for l in c}) != len(c) for c in cls):
            continue
        if nx.check_planarity(incidence(cls))[0]:
            return cls
    return None


def smallest_variant(rng, clause3_only, nvs, tries=60):
    for nv in nvs:
        for _ in range(tries):
            phi = random_exact3_phi(rng, nv, clause3_only)
            if phi is None:
                continue
            if not nx.check_planarity(apex_minor(phi))[0]:
                return analyse(phi, milp=True)
    return None


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    out = {"phi_star": analyse(PHI_STAR, milp=True), "phi0": analyse(PHI0, milp=True)}
    for key, (phi, order) in (("phi_B2", PHI_B2), ("phi_u", PHI_U), ("phi_2u", PHI_2U)):
        out[key] = analyse(phi, milp=True)
        out[key]["variable_cycle_order"] = order
        out[key]["variable_cycle_planar"] = variable_cycle_planar(phi, order)
    # smallest strict counterexample found by random search (var <= 3 times)
    strict = []
    stats = dict(formulas=0, apex_nonplanar=0, G_nonplanar=0, lemma_violations=0)
    for t in range(N):
        nv = rng.randint(5, 14); m = rng.randint(nv, 2 * nv)
        phi = random_strict_phi(rng, nv, m)
        if phi is None:
            continue
        r = analyse(phi, milp=False)
        stats["formulas"] += 1
        stats["apex_nonplanar"] += (not r["apex_minor_planar"])
        stats["G_nonplanar"] += (not r["G_planar"])
        if (not r["apex_minor_planar"]) and r["G_planar"]:
            stats["lemma_violations"] += 1          # would contradict Lemma S1
        if not r["apex_minor_planar"]:
            strict.append(r)
    strict.sort(key=lambda r: (r["m"], r["nv"]))
    out["random_strict"] = stats
    if strict:
        best = strict[0]
        best.update(analyse(best["phi"], milp=True))
        out["smallest_strict"] = best
    out["exact3_sizes23"] = smallest_variant(rng, False, range(3, 12))
    out["exact3_sizes3"] = smallest_variant(rng, True, range(3, 16))
    os.makedirs("supp_runs", exist_ok=True)
    json.dump(out, open("supp_runs/planarity_evidence.json", "w"), indent=1)
    for k in ("phi_star", "phi_B2", "phi_u", "phi_2u", "phi0", "smallest_strict", "exact3_sizes23", "exact3_sizes3"):
        r = out.get(k)
        if r:
            print(k, r["phi"], r["cond"], "incidence planar", r["incidence_planar"], "H planar", r["H_planar"], "split-apex planar",
                  r["apex_minor_planar"], "var-cycle planar", r.get("variable_cycle_planar"), "G planar", r["G_planar"], r.get("kuratowski"), "sat", r["satisfiable"],
                  "milp ok", r.get("milp_equiv_ok"))
    print("random strict:", stats)

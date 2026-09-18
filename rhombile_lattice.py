"""Rhombile lattice (PBC) construction, string-defect insertion, visualization.

Lattice: 3-site basis (r1, r2, r3) per unit cell, bond offsets and geometry
follow rhombile_topo_report.md / rhombile_classical_proposal.md.
"""
import numpy as np
import matplotlib.pyplot as plt

A1 = np.array([1.0, 0.0])
A2 = np.array([-0.5, np.sqrt(3) / 2])
A_MAT = np.column_stack([A1, A2])
A_MAT_INV = np.linalg.inv(A_MAT)

BASIS = {
    "r1": np.array([0.0, 0.0]),
    "r2": np.array([0.5, np.sqrt(3) / 6]),
    "r3": np.array([0.0, np.sqrt(3) / 3]),
}
SUB_INDEX = {"r1": 0, "r2": 1, "r3": 2}

# (site_a, site_b) -> offsets (dn, dm) from a(n, m) to b(n+dn, m+dm)
BOND_OFFSETS = {
    ("r1", "r2"): [(0, 0), (-1, 0), (-1, -1)],
    ("r1", "r3"): [(0, 0), (0, -1), (-1, -1)],
    ("r2", "r3"): [(0, 0), (0, -1), (1, 0)],
}

# range of periodic box translations tried when resolving a wrapping string;
# +-2 comfortably covers the minimum-image displacement for any nx, ny >= 3
WRAP_SHIFT_RANGE = range(-2, 3)


class RhombileLattice:
    def __init__(self, nx, ny, j12=-1.0, j13=-1.0, j23=0.0):
        self.nx, self.ny = nx, ny
        self.j0 = {"r1r2": j12, "r1r3": j13, "r2r3": j23}
        self.n_sites = nx * ny * 3
        self.bonds = []  # dicts: i, j, type, J, p1, p2 (unwrapped positions)
        self._neighbors = None  # lazy cache, see neighbor_table()
        self._build()

    def _site_index(self, n, m, sub):
        return ((m % self.ny) * self.nx + (n % self.nx)) * 3 + SUB_INDEX[sub]

    def _position(self, n, m, sub):
        return n * A1 + m * A2 + BASIS[sub]

    def _build(self):
        for m in range(self.ny):
            for n in range(self.nx):
                for (s1, s2), offsets in BOND_OFFSETS.items():
                    i = self._site_index(n, m, s1)
                    p1 = self._position(n, m, s1)
                    for dn, dm in offsets:
                        j = self._site_index(n + dn, m + dm, s2)
                        p2 = self._position(n + dn, m + dm, s2)
                        self.bonds.append({
                            "i": i, "j": j, "type": s1 + s2,
                            "J": self.j0[s1 + s2], "p1": p1, "p2": p2,
                        })

    def site_positions(self):
        pos = np.zeros((self.n_sites, 2))
        sub_of = np.empty(self.n_sites, dtype=object)
        for m in range(self.ny):
            for n in range(self.nx):
                for sub in BASIS:
                    idx = self._site_index(n, m, sub)
                    pos[idx] = self._position(n, m, sub)
                    sub_of[idx] = sub
        return pos, sub_of

    def neighbor_table(self):
        """list[site] -> [(other_site, bond_dict), ...], each bond referenced
        (not copied) so a later string-defect J flip is seen automatically."""
        if self._neighbors is None:
            nbrs = [[] for _ in range(self.n_sites)]
            for b in self.bonds:
                nbrs[b["i"]].append((b["j"], b))
                nbrs[b["j"]].append((b["i"], b))
            self._neighbors = nbrs
        return self._neighbors


def random_states(lattice, rng):
    """Random initial Potts state (0, 1, or 2) for every site."""
    return rng.integers(0, 3, size=lattice.n_sites)


def total_energy(lattice, states, D):
    """Anisotropic 3-state Potts Hamiltonian, generalizing H = -J S_i.S_j to
    Potts variables via delta(s_i, s_j) in place of the dot product:

        H = sum_<i,j> (-J_ij) * delta(s_i, s_j)  +  sum_i D[s_i]

    J_ij < 0 (the lattice default) therefore means antiferromagnetic:
    matching neighbors cost +|J_ij|, mismatched neighbors cost 0.
    D = (D_1, D_2, D_3) is the single-site anisotropy energy of each state.
    """
    e = float(np.asarray(D, dtype=float)[states].sum())
    for b in lattice.bonds:
        if states[b["i"]] == states[b["j"]]:
            e -= b["J"]
    return e


def state_energies(lattice, states, site, D):
    """Energy of `site` for each candidate state (0, 1, 2), holding every
    other site fixed -- the effective field used to drive Metropolis/SA
    updates. The cost of flipping `site` to state k is

        delta_E = state_energies(...)[k] - state_energies(...)[states[site]]
    """
    e = np.asarray(D, dtype=float).copy()
    for j, b in lattice.neighbor_table()[site]:
        e[states[j]] -= b["J"]
    return e


def frustrated_triangles(lattice):
    """Enumerate every elementary (r1, r2, r3) triangle and return the ones
    that are topologically frustrated: all 3 edges active (nonzero J). Such
    a triangle is an odd all-AF cycle -- no spin assignment can satisfy
    every bond on it, regardless of what the rest of the lattice does.

    Returns a list of dicts {"r1", "a", "b", "centroid"} (site indices of
    the triangle's corners and its plotting centroid).
    """
    bond_by_pair = {frozenset((b["i"], b["j"])): b for b in lattice.bonds}
    result = []
    for m in range(lattice.ny):
        for n in range(lattice.nx):
            r1_idx = lattice._site_index(n, m, "r1")
            r1_pos = lattice._position(n, m, "r1")
            nbrs = []
            for other_idx, b in lattice.neighbor_table()[r1_idx]:
                other_pos = b["p2"] if b["i"] == r1_idx else b["p1"]
                ang = np.arctan2(*(other_pos - r1_pos)[::-1])
                nbrs.append((ang, other_idx, other_pos, b["J"]))
            nbrs.sort(key=lambda x: x[0])
            k = len(nbrs)
            for i in range(k):
                _, idx1, pos1, J1 = nbrs[i]
                _, idx2, pos2, J2 = nbrs[(i + 1) % k]
                closing = bond_by_pair.get(frozenset((idx1, idx2)))
                if closing is None:
                    continue
                if J1 != 0 and J2 != 0 and closing["J"] != 0:
                    centroid = (r1_pos + pos1 + pos2) / 3
                    result.append({"r1": r1_idx, "a": idx1, "b": idx2, "centroid": centroid})
    return result


def matched_bonds(lattice, states):
    """Bonds whose two endpoints share the same state: the ones actually
    contributing +|J| to total_energy (the real physical cost), as opposed
    to bonds that were merely toggled when a string defect was drawn."""
    return [b for b in lattice.bonds if b["J"] != 0.0 and states[b["i"]] == states[b["j"]]]


def heat_bath_sweep(lattice, states, D, T, rng):
    """One sweep = resample every site's state from the Gibbs distribution
    p(k) ~ exp(-state_energies(...)[k] / T), in random order, in place."""
    for site in rng.permutation(lattice.n_sites):
        e = state_energies(lattice, states, site, D)
        w = np.exp(-(e - e.min()) / T)
        states[site] = rng.choice(3, p=w / w.sum())
    return states


def simulated_annealing(lattice, D, rng, n_sweeps=400, T_start=3.0, T_end=0.01,
                         states=None, record_energy=False):
    """Anneal from T_start down to T_end over n_sweeps heat-bath sweeps
    (geometric schedule). Returns the final states, and the energy trace
    (per sweep) if record_energy=True."""
    if states is None:
        states = random_states(lattice, rng)
    else:
        states = states.copy()
    energies = [] if record_energy else None
    for T in np.geomspace(T_start, T_end, n_sweeps):
        heat_bath_sweep(lattice, states, D, T, rng)
        if record_energy:
            energies.append(total_energy(lattice, states, D))
    return states, energies


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def segments_intersect(p1, p2, p3, p4):
    d1, d2 = _cross(p3, p4, p1), _cross(p3, p4, p2)
    d3, d4 = _cross(p1, p2, p3), _cross(p1, p2, p4)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def minimum_image_endpoint(lattice, p_start, p_end):
    """Periodic image of p_end closest to p_start, i.e. the endpoint of the
    shortest string between the two points on the torus."""
    frac_start = A_MAT_INV @ p_start
    frac_end = A_MAT_INV @ p_end
    du, dv = frac_end - frac_start
    du -= lattice.nx * np.round(du / lattice.nx)
    dv -= lattice.ny * np.round(dv / lattice.ny)
    return p_start + du * A1 + dv * A2


def apply_string_defect(lattice, p_start, p_end, wrap=True):
    """Toggle J (-1 <-> 0) on every bond crossed by the string p_start -> p_end.

    With wrap=True, p_end is first replaced by its minimum-image copy (the
    shortest string on the torus), and the crossing test is repeated over a
    small set of periodic box translations, so a string that runs off one
    edge of the simulation box and reappears on the opposite edge is handled
    correctly. With wrap=False the string is treated as a plain segment in
    the plane (no PBC wrap-around).
    """
    if wrap:
        p_end = minimum_image_endpoint(lattice, p_start, p_end)
        shifts = [kx * lattice.nx * A1 + ky * lattice.ny * A2
                  for kx in WRAP_SHIFT_RANGE for ky in WRAP_SHIFT_RANGE]
    else:
        shifts = [np.zeros(2)]

    flipped = []
    for shift in shifts:
        s_start, s_end = p_start - shift, p_end - shift
        for b in lattice.bonds:
            if segments_intersect(s_start, s_end, b["p1"], b["p2"]):
                b["J"] = 0.0 if b["J"] != 0.0 else -1.0
                flipped.append(b)
    return flipped, p_end


def clip_string_for_display(lattice, p_start, p_end_image, min_length=0.05):
    """Split the (possibly out-of-box) string p_start -> p_end_image at box
    boundaries and translate each piece back into the fundamental box, so a
    wrapping string can be drawn as separate segments re-entering the plot.

    When the string grazes a corner of the box, u and v can cross their
    periodic boundaries at two slightly different parameter values, briefly
    assigning a sliver of the path to a diagonally-adjacent cell copy; that
    sliver is real (the straight line genuinely passes through it) but only
    spans a fraction of a bond length, so pieces shorter than min_length are
    dropped here to avoid a stray dash far from the rest of the string."""
    frac0 = A_MAT_INV @ p_start
    frac1 = A_MAT_INV @ p_end_image
    du, dv = frac1 - frac0
    ts = {0.0, 1.0}
    if du != 0:
        lo, hi = sorted([frac0[0], frac1[0]])
        for k in range(int(np.floor(lo / lattice.nx)), int(np.ceil(hi / lattice.nx)) + 1):
            t = (k * lattice.nx - frac0[0]) / du
            if 0 < t < 1:
                ts.add(t)
    if dv != 0:
        lo, hi = sorted([frac0[1], frac1[1]])
        for k in range(int(np.floor(lo / lattice.ny)), int(np.ceil(hi / lattice.ny)) + 1):
            t = (k * lattice.ny - frac0[1]) / dv
            if 0 < t < 1:
                ts.add(t)
    ts = sorted(ts)

    segments = []
    for t0, t1 in zip(ts[:-1], ts[1:]):
        tm = 0.5 * (t0 + t1)
        fm = frac0 + tm * np.array([du, dv])
        kx, ky = int(np.floor(fm[0] / lattice.nx)), int(np.floor(fm[1] / lattice.ny))
        shift = kx * lattice.nx * A1 + ky * lattice.ny * A2
        q0 = p_start + t0 * (p_end_image - p_start) - shift
        q1 = p_start + t1 * (p_end_image - p_start) - shift
        if np.linalg.norm(q1 - q0) >= min_length:
            segments.append((q0, q1))
    return segments


def random_site_pair(lattice, rng):
    """Two random points = two random lattice sites."""
    pos, _ = lattice.site_positions()
    i, j = rng.choice(lattice.n_sites, size=2, replace=False)
    return pos[i], pos[j]


def random_point_pair(lattice, rng):
    """Two random points = two arbitrary continuous points inside the box."""
    frac = rng.uniform(0, 1, size=(2, 2)) * [lattice.nx, lattice.ny]
    return frac[0, 0] * A1 + frac[0, 1] * A2, frac[1, 0] * A1 + frac[1, 1] * A2


def plot_lattice(lattice, string_endpoints=None, wrap=True, ax=None, title=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    for b in lattice.bonds:
        if b["J"] != 0.0:
            ax.plot([b["p1"][0], b["p2"][0]], [b["p1"][1], b["p2"][1]],
                     color="black", linewidth=1.2, zorder=1)
    pos, sub_of = lattice.site_positions()
    colors = {"r1": "tab:red", "r2": "tab:blue", "r3": "tab:green"}
    for sub, c in colors.items():
        mask = sub_of == sub
        ax.scatter(pos[mask, 0], pos[mask, 1], s=25, color=c, label=sub, zorder=2)
    if string_endpoints is not None:
        p_start, p_end = string_endpoints
        p_end_image = minimum_image_endpoint(lattice, p_start, p_end) if wrap else p_end
        for k, (q0, q1) in enumerate(clip_string_for_display(lattice, p_start, p_end_image)):
            ax.plot([q0[0], q1[0]], [q0[1], q1[1]], "--", color="orange",
                     linewidth=2, zorder=3, label="string" if k == 0 else None)
        ax.scatter(*p_start, marker="x", color="orange", s=100, zorder=4)
        ax.scatter(*p_end, marker="x", color="orange", s=100, zorder=4)
    ax.set_aspect("equal")
    ax.legend(loc="upper right")
    if title:
        ax.set_title(title)
    return ax


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    nx, ny = 8, 8

    # 1. site-based string endpoints, PBC wrap allowed
    lat = RhombileLattice(nx, ny)
    p_start, p_end = random_site_pair(lat, rng)
    flipped, _ = apply_string_defect(lat, p_start, p_end, wrap=True)
    print(f"[site pair] {p_start} -> {p_end}: flipped {len(flipped)} bonds")
    ax = plot_lattice(lat, string_endpoints=(p_start, p_end),
                       title=f"Rhombile {nx}x{ny} (PBC): site-pair string")
    plt.savefig("rhombile_site_string.png", dpi=150)
    plt.close()

    # 2. continuous random-point string endpoints
    lat2 = RhombileLattice(nx, ny)
    p_start2, p_end2 = random_point_pair(lat2, rng)
    flipped2, _ = apply_string_defect(lat2, p_start2, p_end2, wrap=True)
    print(f"[continuous pair] {p_start2} -> {p_end2}: flipped {len(flipped2)} bonds")
    ax = plot_lattice(lat2, string_endpoints=(p_start2, p_end2),
                       title=f"Rhombile {nx}x{ny} (PBC): continuous-point string")
    plt.savefig("rhombile_continuous_string.png", dpi=150)
    plt.close()

    # 3. deliberately wrapping string (endpoints near opposite edges)
    lat3 = RhombileLattice(nx, ny)
    p_a = lat3._position(1, ny // 2, "r1")
    p_b = lat3._position(nx - 2, ny // 2, "r1")
    flipped3, p_b_image = apply_string_defect(lat3, p_a, p_b, wrap=True)
    print(f"[wrap demo] {p_a} -> {p_b} (image {p_b_image}): flipped {len(flipped3)} bonds")
    ax = plot_lattice(lat3, string_endpoints=(p_a, p_b),
                       title=f"Rhombile {nx}x{ny} (PBC): wrap-around string")
    plt.savefig("rhombile_wrap_string.png", dpi=150)
    plt.close()

    print("saved rhombile_site_string.png, rhombile_continuous_string.png, rhombile_wrap_string.png")

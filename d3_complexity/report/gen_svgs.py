"""Generate the inline SVG figures for d3_explainer.html (python3 gen_svgs.py -> figs.json).

Colours are CSS classes (c1, c2, c3) styled by the page, so the figures follow the theme.
Same-colour ("frustrated") bonds are drawn with class 'frus' (thick dashed ink).
"""
import json
import math

FIG = {}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;")


def graph_svg(pos, edges, col, w, h, r=11, labels=None, extra="", aria="", digits=True, dim=None):
    """pos: {v:(x,y)}, edges: [(u,v)], col: {v: 1|2|3|None}."""
    out = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(aria)}">']
    for u, v in edges:
        (x1, y1), (x2, y2) = pos[u], pos[v]
        f = col.get(u) is not None and col.get(u) == col.get(v)
        cls = "frus" if f else "bond"
        out.append(f'<line class="{cls}" x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"/>')
    out.append(extra)
    for v, (x, y) in pos.items():
        c = col.get(v)
        cls = f"c{c}" if c else "c0"
        if dim and v in dim:
            cls += " dim"
        out.append(f'<circle class="node {cls}" cx="{x:.1f}" cy="{y:.1f}" r="{r}"/>')
        if digits and c:
            out.append(f'<text class="nd n{c}" x="{x:.1f}" y="{y + 4:.1f}">{c}</text>')
        if labels and v in labels:
            lx, ly, t = labels[v]
            out.append(f'<text class="lab" x="{x + lx:.1f}" y="{y + ly:.1f}">{esc(t)}</text>')
    out.append("</svg>")
    return "".join(out)


# ---------------------------------------------------------------- 1. the triangle
def triangle(cols, title_y=None):
    pos = {0: (60, 22), 1: (22, 88), 2: (98, 88)}
    return pos, [(0, 1), (1, 2), (0, 2)], dict(enumerate(cols))


def tri_panel(cols, ox):
    pos, E, col = triangle(cols)
    pos = {k: (x + ox, y + 10) for k, (x, y) in pos.items()}
    return pos, E, col


def fig_triangle():
    parts = []
    W, H = 560, 150
    allpos, allE, allc = {}, [], {}
    panels = [((1, 2, 1), 20, "두 색만: 같은 색 결합 1개", "에너지 1"),
              ((1, 2, 2), 200, "두 색만: 어떻게 칠해도 1개", "에너지 1"),
              ((1, 2, 3), 380, "세 번째 색을 쓰면: 0개", "에너지 D3")]
    svg = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="삼각형은 두 색으로 칠하면 반드시 같은 색 결합이 하나 생기고, 세 번째 색을 쓰면 결합 비용은 0이 되는 대신 D3를 낸다">']
    for i, (cols, ox, cap, en) in enumerate(panels):
        pos, E, col = tri_panel(cols, ox + 20)
        g = graph_svg({(i, k): p for k, p in pos.items()}, [((i, a), (i, b)) for a, b in E],
                      {(i, k): c for k, c in col.items()}, W, H)
        svg.append(g[g.index(">") + 1:-6])
        svg.append(f'<text class="cap" x="{ox + 80}" y="128">{cap}</text>')
        svg.append(f'<text class="cap strong" x="{ox + 80}" y="145">{en}</text>')
    svg.append("</svg>")
    return "".join(svg)


FIG["triangle"] = fig_triangle()


# ---------------------------------------------------------------- 2. triangular lattice: two limits
def tri_lattice(colfun, nx=9, ny=6, s=34, ox=18, oy=20):
    pos, col, E = {}, {}, []
    for j in range(ny):
        for i in range(nx):
            x = ox + s * (i + 0.5 * (j % 2))
            y = oy + s * 0.866 * j
            pos[(i, j)] = (x, y)
            col[(i, j)] = colfun(i, j)
    for j in range(ny):
        for i in range(nx):
            nbrs = [(i + 1, j)]
            if j % 2 == 0:
                nbrs += [(i, j + 1), (i - 1, j + 1)]
            else:
                nbrs += [(i, j + 1), (i + 1, j + 1)]
            for q in nbrs:
                if q in pos:
                    E.append(((i, j), q))
    return pos, E, col


def fig_trilattice():
    W, H = 330, 200
    # proper 3-colouring on the offset grid: axial coordinate a = i - floor(j/2)
    def three(i, j):
        a = i - (j - (j % 2)) // 2 + (0 if j % 2 == 0 else 0)
        return 1 + ((a * 1 + j * 2) % 3)
    # fix: use cube-ish coordinates for offset rows (odd rows shifted right)
    def three2(i, j):
        q = i - (j - (j & 1)) // 2
        return 1 + ((q - j) % 3)
    def stripes(i, j):
        return 1 if j % 2 == 0 else 2
    out = {}
    for name, f, aria in (("tri3", three2, "삼각 격자를 세 색으로 칠하면 같은 색 결합이 하나도 없고, 색 3이 사이트의 1/3을 차지한다"),
                          ("tri2", stripes, "두 색 줄무늬 배치: 모든 삼각형에 같은 색 결합이 정확히 하나씩 있다")):
        pos, E, col = tri_lattice(f)
        out[name] = graph_svg(pos, E, col, W, H, r=9, aria=aria)
    return out


FIG.update(fig_trilattice())


# ---------------------------------------------------------------- 3. the key principle: recolouring a colour-3 site
def star(center_col, nb_cols, cx=110, cy=100, R=62, start=-90, w=220, h=200, labels=None, aria=""):
    d = len(nb_cols)
    pos = {"c": (cx, cy)}
    for k in range(d):
        a = math.radians(start + 360 * k / d)
        pos[k] = (cx + R * math.cos(a), cy + R * math.sin(a))
    E = [("c", k) for k in range(d)]
    col = {"c": center_col}
    col.update({k: c for k, c in enumerate(nb_cols)})
    return graph_svg(pos, E, col, w, h, r=13, aria=aria, labels=labels)


FIG["star_before"] = star(3, [1, 2, 1, 2, 1, 2], aria="가운데가 색 3이고 이웃이 1,2,1,2,1,2: 가운데 결합은 모두 정상")
FIG["star_after1"] = star(1, [1, 2, 1, 2, 1, 2], aria="가운데를 색 1로 바꾸면 같은 색 결합이 3개 생긴다")
FIG["star_uneven"] = star(2, [1, 1, 1, 1, 2, 1], aria="이웃이 한쪽 색에 몰려 있으면 가운데를 적은 쪽 색으로 바꿔 거의 비용 없이 색 3을 없앨 수 있다")


# ---------------------------------------------------------------- 4. wheels
def wheel(D, hub, rim, w=210, h=200, aria=""):
    cx, cy, R = w / 2, 96, 70
    pos = {"h": (cx, cy)}
    for k in range(D):
        a = math.radians(-90 + 360 * k / D)
        pos[k] = (cx + R * math.cos(a), cy + R * math.sin(a))
    E = [("h", k) for k in range(D)] + [(k, (k + 1) % D) for k in range(D)]
    col = {"h": hub}
    col.update({k: c for k, c in enumerate(rim)})
    return graph_svg(pos, E, col, w, h, r=12, aria=aria)


FIG["w4_two"] = wheel(4, 1, [2, 1, 2, 1], aria="바퀴 W4를 두 색으로 칠한 최선: 같은 색 결합 2개")
FIG["w4_hub"] = wheel(4, 3, [2, 1, 2, 1], aria="W4의 중심을 색 3으로: 같은 색 결합 0개, 비용 D3")
FIG["w5_two"] = wheel(5, 1, [2, 1, 2, 1, 2], aria="바퀴 W5를 두 색으로 칠한 최선: 같은 색 결합 3개")
FIG["w5_hub"] = wheel(5, 3, [2, 1, 2, 1, 2], aria="W5의 중심을 색 3으로: 테두리의 같은 색 결합 1개와 비용 D3")


# ---------------------------------------------------------------- 5. local structure (Delta = 4, 1 < D3 < 2)
FIG["loc_ok"] = star(3, [1, 2, 1, 2], w=170, h=170, cx=85, cy=82, R=52, start=-45, aria="허용: 색 3 사이트의 네 이웃이 1,2,1,2")
FIG["loc_ok2"] = star(3, [1, 1, 2, 2], w=170, h=170, cx=85, cy=82, R=52, start=-45, aria="허용: 네 이웃이 1,1,2,2 (순서는 상관없음)")
FIG["loc_bad1"] = star(3, [1, 2, 2, 2], w=170, h=170, cx=85, cy=82, R=52, start=-45, aria="금지: 이웃이 1,2,2,2이면 가운데를 1로 바꾸는 편이 싸다")
FIG["loc_bad2"] = star(3, [1, 2, 3, 2], w=170, h=170, cx=85, cy=82, R=52, start=-45, aria="금지: 색 3끼리 이웃")
FIG["loc_bad3"] = star(3, [1, 2, 1], w=170, h=170, cx=85, cy=82, R=52, start=-90, aria="금지: 차수 3인 사이트는 색 3이 될 수 없다")


# ---------------------------------------------------------------- 6. charts
def line_chart(lines, xmax, ymax, w=520, h=260, xt=None, yt=None, xlabel="D3", ylabel="", mark=None, aria="", envelope=True):
    """lines: [(label, slope, intercept, cls)] ; draws each line and the lower envelope."""
    L, Rm, T, B = 52, 18, 16, 40
    X = lambda x: L + (w - L - Rm) * x / xmax
    Y = lambda y: h - B - (h - T - B) * y / ymax
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(aria)}">']
    for t in (yt or []):
        o.append(f'<line class="grid" x1="{L}" y1="{Y(t):.1f}" x2="{w - Rm}" y2="{Y(t):.1f}"/>')
        o.append(f'<text class="tick" x="{L - 8}" y="{Y(t) + 4:.1f}" text-anchor="end">{t:g}</text>')
    for t in (xt or []):
        o.append(f'<text class="tick" x="{X(t):.1f}" y="{h - B + 18}" text-anchor="middle">{t:g}</text>')
    o.append(f'<line class="axis" x1="{L}" y1="{h - B}" x2="{w - Rm}" y2="{h - B}"/>')
    o.append(f'<line class="axis" x1="{L}" y1="{T}" x2="{L}" y2="{h - B}"/>')
    o.append(f'<text class="tick" x="{(L + w - Rm) / 2}" y="{h - 6}" text-anchor="middle">{esc(xlabel)}</text>')
    o.append(f'<text class="tick" x="14" y="{(T + h - B) / 2}" text-anchor="middle" transform="rotate(-90 14 {(T + h - B) / 2})">{esc(ylabel)}</text>')
    for lab, a, b, cls, (lx, ly) in lines:
        x2 = min(xmax, (ymax - b) / a) if a > 0 else xmax
        o.append(f'<line class="ln {cls}" x1="{X(0):.1f}" y1="{Y(b):.1f}" x2="{X(x2):.1f}" y2="{Y(a * x2 + b):.1f}"/>')
        o.append(f'<text class="lnlab" x="{X(lx):.1f}" y="{Y(ly):.1f}">{esc(lab)}</text>')
    if envelope:
        pts = []
        N = 400
        for k in range(N + 1):
            x = xmax * k / N
            y = min(a * x + b for _, a, b, _, _ in lines)
            pts.append(f"{X(x):.1f},{Y(y):.1f}")
        o.append(f'<polyline class="env" points="{" ".join(pts)}"/>')
    if mark:
        mx, my, t = mark
        o.append(f'<circle class="kink" cx="{X(mx):.1f}" cy="{Y(my):.1f}" r="6"><title>{esc(t)}</title></circle>')
        o.append(f'<text class="lnlab strong" x="{X(mx) + 10:.1f}" y="{Y(my) + 22:.1f}">{esc(t)}</text>')
    o.append("</svg>")
    return "".join(o)


FIG["chart_w5"] = line_chart(
    [("두 색만: 3", 0, 3, "l2", (0.15, 3.12)), ("중심을 색 3으로: D3 + 1", 1, 1, "l3", (0.45, 0.85))],
    xmax=3.5, ymax=4, xt=[0, 1, 2, 3], yt=[0, 1, 2, 3, 4], ylabel="바닥상태 에너지",
    mark=(2, 3, "교차점 D3 = 2 = ⌊5/2⌋"), aria="바퀴 W5: 두 색 해의 에너지 3과 중심 색 3 해의 에너지 D3+1이 D3=2에서 교차한다")

FIG["chart_tri"] = line_chart(
    [("Wannier 두 색 배치: 1", 0, 1, "l2", (0.15, 1.05)), ("세 색 배치: D3/3", 1 / 3, 0, "l3", (1.3, 0.34))],
    xmax=4.5, ymax=1.4, xt=[0, 1, 2, 3, 4], yt=[0, 0.5, 1], ylabel="사이트당 에너지",
    mark=(3, 1, "꺾임 D3* = 3 = ⌊6/2⌋"), aria="삼각 격자: 사이트당 에너지가 min(D3/3, 1)이며 D3=3에서 꺾인다")

# Penrose dual (radius-6 patch, 140 rhombi): exact segments from potts_exact.curve
pen = [(0, 37), (5, 32), (10, 28), (66, 0)]
FIG["chart_pen"] = line_chart(
    [({(66, 0): "상태 3 없음 (단색 66)", (10, 28): "상태 3 28개 (단색 10)"}.get((m, k), ""), k / 140, m / 140, "l3" if k else "l2",
      {(0, 37): (0, 0), (5, 32): (0, 0), (10, 28): (2.12, 0.585), (66, 0): (0.05, 0.505)}[(m, k)])
     for m, k in pen],
    xmax=3, ymax=0.62, xt=[0, 0.5, 1, 1.25, 1.5, 2, 2.5, 3], yt=[0, 0.2, 0.4, 0.6], ylabel="마름모당 에너지",
    mark=(2, 66 / 140, "D3* = 2 = ⌊4/2⌋"), aria="펜로즈 쌍대 조각의 정확한 바닥상태 에너지: D3=1, 1.25에서 작게 꺾이고 D3=2에서 색 3이 모두 사라진다")


# ---------------------------------------------------------------- 7. Penrose dual ground states
def penrose_svg(state_key, w=360):
    d = json.load(open("penrose_states.json"))
    pts = [p for poly in d["polys"] for p in poly]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    sc = (w - 20) / (x1 - x0)
    h = (y1 - y0) * sc + 20
    T = lambda p: (10 + (p[0] - x0) * sc, 10 + (y1 - p[1]) * sc)
    st = d["states"][state_key]
    col = st["col"]
    o = [f'<svg viewBox="0 0 {w} {h:.0f}" role="img" aria-label="펜로즈 마름모 타일 조각의 바닥상태 (D3 = {state_key}): 색 3 마름모 {st["n3"]}개, 같은 색이 맞닿은 변 {st["mono"]}개">']
    for i, poly in enumerate(d["polys"]):
        P = " ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in poly)
        cls = f"t{col[i]}" + (" edge" if d["deg"][i] < 4 else "")
        o.append(f'<polygon class="tile {cls}" points="{P}"/>')
    for a, b, seg in d["edges"]:
        if col[a] == col[b]:
            (xa, ya), (xb, yb) = T(seg[0]), T(seg[1])
            o.append(f'<line class="frusedge" x1="{xa:.1f}" y1="{ya:.1f}" x2="{xb:.1f}" y2="{yb:.1f}"/>')
    o.append("</svg>")
    return "".join(o)


FIG["pen05"] = penrose_svg("0.5")
FIG["pen15"] = penrose_svg("1.5")
FIG["pen25"] = penrose_svg("2.5")

json.dump(FIG, open("figs.json", "w"), ensure_ascii=False)
print("figures:", len(FIG), {k: len(v) for k, v in FIG.items()})


# ---------------------------------------------------------------- 8. degree-3 lattices: healing (Theorem C)
import sys as _sys
_sys.path.insert(0, "..")
from ising_tjoin import ising_ground_state_tjoin
from heal_check import two_sat_heal


def fig_heal_prism():
    # triangular prism: outer triangle 0,1,2 and inner triangle 3,4,5 with spokes i -- i+3
    R1, R2, cx, cy = 78, 30, 110, 104
    pos = {}
    for i in range(3):
        a = math.radians(-90 + 120 * i)
        pos[i] = (cx + R1 * math.cos(a), cy + R1 * math.sin(a))
        pos[i + 3] = (cx + R2 * math.cos(a), cy + R2 * math.sin(a))
    E = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)]
    fr, col = ising_ground_state_tjoin(6, E)
    chosen, M = two_sat_heal(6, E, col)
    healed = list(col)
    for v in chosen:
        healed[v] = 3
    a = graph_svg(pos, E, dict(enumerate(col)), 220, 200, r=12,
                  aria=f"삼각 기둥의 Ising 바닥상태: 좌절 결합 {fr}개")
    b = graph_svg(pos, E, dict(enumerate(healed)), 220, 200, r=12,
                  aria="좌절 결합마다 한 끝점을 서로 이웃하지 않게 상태 3으로 바꾼 치유 배치: 같은 상태 결합 0개")
    return a, b, fr


FIG["heal_ising"], FIG["heal_done"], _fr = fig_heal_prism()


def fig_truncation():
    cx, cy, R = 100, 90, 62
    pos = {"c": (cx, cy)}
    for k in range(5):
        a = math.radians(-90 + 72 * k)
        pos[k] = (cx + R * math.cos(a), cy + R * math.sin(a))
    before = graph_svg(pos, [("c", k) for k in range(5)], {}, 200, 180, r=8,
                       aria="이웃이 5개인 꼭짓점", digits=False)
    pos2 = {}
    for k in range(5):
        a = math.radians(-90 + 72 * k)
        pos2[("in", k)] = (cx + 24 * math.cos(a), cy + 24 * math.sin(a))
        pos2[("out", k)] = (cx + R * math.cos(a), cy + R * math.sin(a))
    E2 = [(("in", k), ("in", (k + 1) % 5)) for k in range(5)] + [(("in", k), ("out", k)) for k in range(5)]
    after = graph_svg(pos2, E2, {}, 200, 180, r=8, aria="깎은 뒤: 꼭짓점이 오각형이 되고 각 사이트의 이웃은 3개", digits=False)
    return before, after


FIG["trunc_before"], FIG["trunc_after"] = fig_truncation()

FIG["chart_tpen"] = line_chart(
    [("두 색만 (Ising): fr/N = 0.2", 0, 0.2, "l2", (0.08, 0.212)), ("치유 배치: 0.2·D3", 0.2, 0, "l3", (0.45, 0.06))],
    xmax=1.6, ymax=0.26, xt=[0, 0.5, 1, 1.5], yt=[0, 0.1, 0.2], ylabel="사이트당 에너지",
    mark=(1, 0.2, "꺾임은 D3 = 1에서 한 번"), aria="깎은 펜로즈(사이트 270개, 좌절 54): 정확한 에너지는 min(D3,1)·54이며 D3=1에서만 꺾인다")

# gallery panels (class names prefixed so they do not collide with the page's graph styles)
from deg3_gallery import LATTICES, geometric_faces, svg_for
from heal_check import healed_state
_gal = {}
for key, name, desc, fn in LATTICES:
    if key not in ("tpen", "foam", "star"):
        continue
    n, E, pos, info = fn()
    fr, s = healed_state(n, E)
    faces = geometric_faces(n, E, [tuple(p) for p in pos])
    svg = svg_for(n, E, pos, faces, s, w=360)
    svg = svg.replace('class="bond"', 'class="gbond"').replace('class="odd"', 'class="godd"').replace('class="site c', 'class="gsite c')
    FIG["gal_" + key] = svg
    _gal[key] = (n, fr, sum(1 for x in s if x == 3))
FIG["_galstats"] = json.dumps(_gal)

json.dump(FIG, open("figs.json", "w"), ensure_ascii=False)
print("added degree-3 figures", _gal)

"""Gallery of complex max-degree-3 lattices with their certified ground states (0 < D3 < 1).

Lattices: star (3.12.12), truncated rhombille, truncated Penrose, defect honeycomb, 2D foam.
Each is drawn with odd faces shaded (the frustration charges), sites coloured by the
certified ground state from heal_check.healed_state (T-join matching + 2-SAT).
Writes deg3_gallery.html (self-contained, inline SVG).
"""
import math
import sys
import time

import numpy as np

sys.path.insert(0, "..")
from lattices3 import truncated, truncated_penrose, defect_honeycomb, random_cubic
from heal_check import healed_state


def geometric_faces(n, edges, pos):
    """Faces of a straight-line planar drawing (outer face dropped)."""
    nbrs = {v: [] for v in range(n)}
    for u, v in edges:
        nbrs[u].append(v); nbrs[v].append(u)
    for v in nbrs:
        x0, y0 = pos[v]
        nbrs[v].sort(key=lambda w: math.atan2(pos[w][1] - y0, pos[w][0] - x0))
    seen = set(); faces = []
    for u in range(n):
        for v in nbrs[u]:
            if (u, v) in seen:
                continue
            f = []; a, b = u, v
            while (a, b) not in seen:
                seen.add((a, b)); f.append(a)
                L = nbrs[b]; i = L.index(a)
                a, b = b, L[(i - 1) % len(L)]
            faces.append(f)
    def area(f):
        return 0.5 * sum(pos[f[i]][0] * pos[f[(i + 1) % len(f)]][1] - pos[f[(i + 1) % len(f)]][0] * pos[f[i]][1] for i in range(len(f)))
    faces = [f for f in faces if area(f) > 0]          # keep one orientation (inner faces)
    big = max(faces, key=lambda f: abs(area(f)))
    return [f for f in faces if f is not big]


def star_lattice(L=7):
    """Honeycomb (Delaunay dual of triangular-lattice points), then truncated."""
    from lattices3 import delaunay_dual
    pts = np.array([(i + 0.5 * (j % 2), j * math.sqrt(3) / 2) for j in range(L) for i in range(L)])
    n, E, pos, info = delaunay_dual(pts)
    return truncated({v: tuple(pos[v]) for v in range(n)}, E)


def truncated_rhombille(L=5):
    """Rhombille tiling: triangular-lattice points (degree 6) + triangle centres (degree 3)."""
    from scipy.spatial import Delaunay
    pts = np.array([(i + 0.5 * j, j * math.sqrt(3) / 2) for j in range(L) for i in range(L)])
    tri = Delaunay(pts)
    Vpos = {("p", k): tuple(p) for k, p in enumerate(pts)}
    E = []
    for t_idx, t in enumerate(tri.simplices):
        c = pts[t].mean(axis=0)
        Vpos[("c", t_idx)] = tuple(c)
        for k in t:
            E.append((("c", t_idx), ("p", int(k))))
    return truncated(Vpos, E)


LATTICES = [
    ("star", "별 격자 (3.12.12)", "벌집의 각 꼭짓점을 삼각형으로 바꾼 주기 격자. 모든 삼각형이 좌절 전하다.", star_lattice),
    ("trhomb", "깎은 마름모 격자", "프로젝트의 마름모(rhombille) 격자를 깎은 것. 이웃 3인 꼭짓점은 삼각형, 이웃 6인 꼭짓점은 육각형이 된다.", truncated_rhombille),
    ("tpen", "깎은 펜로즈 준결정", "펜로즈 P3 타일을 깎은 것. 이웃 3·5·7인 꼭짓점 자리에 삼각형·오각형·칠각형이 준주기적으로 놓인다.", lambda: truncated_penrose(4.5)),
    ("defect", "결함 있는 벌집", "삼각 격자에서 점 8%를 지운 뒤 쌍대를 잡은 것. 결함 주위에 오각형·칠각형·팔각형이 생긴다.", lambda: defect_honeycomb(12, 0.08, 3)),
    ("foam", "2차원 거품 (무작위 Voronoi)", "무작위 점의 Voronoi 분할. 비누 거품처럼 모든 벽이 세 갈래로 만난다.", lambda: random_cubic(90, 5)),
]


def svg_for(n, edges, pos, faces, s, w=420):
    P = np.array(pos, float)
    x0, y0 = P.min(axis=0); x1, y1 = P.max(axis=0)
    sc = (w - 24) / max(x1 - x0, y1 - y0)
    h = (y1 - y0) * sc + 24
    T = lambda v: (12 + (P[v][0] - x0) * sc, 12 + (y1 - P[v][1]) * sc)
    o = [f'<svg viewBox="0 0 {w} {h:.0f}" role="img">']
    for f in faces:
        if len(f) % 2 == 1:
            pts = " ".join(f"{T(v)[0]:.1f},{T(v)[1]:.1f}" for v in f)
            o.append(f'<polygon class="odd" points="{pts}"><title>{len(f)}각형 (홀수 면)</title></polygon>')
    for u, v in edges:
        (a, b), (c, d) = T(u), T(v)
        o.append(f'<line class="bond" x1="{a:.1f}" y1="{b:.1f}" x2="{c:.1f}" y2="{d:.1f}"/>')
    r = max(2.2, min(5.0, sc * 0.09))
    for v in range(n):
        x, y = T(v)
        rr = r * 1.45 if s[v] == 3 else r
        o.append(f'<circle class="site c{s[v]}" cx="{x:.1f}" cy="{y:.1f}" r="{rr:.1f}"/>')
    o.append("</svg>")
    return "".join(o)


def main():
    cards = []
    for key, name, desc, fn in LATTICES:
        n, E, pos, info = fn()
        t = time.time()
        fr, s = healed_state(n, E)
        dt = time.time() - t
        faces = geometric_faces(n, E, [tuple(p) for p in pos])
        odd = sum(1 for f in faces if len(f) % 2 == 1)
        sizes = sorted({len(f) for f in faces})
        ok = s is not None
        if not ok:
            s = [1] * n
        n3 = sum(1 for x in s if x == 3)
        cards.append(dict(key=key, name=name, desc=desc, svg=svg_for(n, E, pos, faces, s), n=n, fr=fr, n3=n3,
                          odd=odd, sizes=sizes, ok=ok, dt=dt))
        print(key, n, "fr", fr, "n3", n3, "odd faces", odd, "sizes", sizes, "certified", ok, f"{dt:.1f}s", flush=True)
    html = open("deg3_gallery_template.html").read()
    body = []
    for c in cards:
        sz = ", ".join(str(k) for k in c["sizes"])
        body.append(f'''<section class="card"><h2>{c["name"]}</h2><p class="desc">{c["desc"]}</p>
<figure>{c["svg"]}</figure>
<table><tr><th>사이트</th><td>{c["n"]}</td><th>면의 크기</th><td>{sz}</td></tr>
<tr><th>홀수 면(좌절 전하)</th><td>{c["odd"]}</td><th>Ising 좌절 결합 fr</th><td>{c["fr"]}</td></tr>
<tr><th>상태 3 사이트</th><td>{c["n3"]}</td><th>인증</th><td>{"성공" if c["ok"] else "실패"} ({c["dt"]:.1f}초)</td></tr></table></section>''')
    open("deg3_gallery.html", "w").write(html.replace("{{CARDS}}", "\n".join(body)))


if __name__ == "__main__":
    main()

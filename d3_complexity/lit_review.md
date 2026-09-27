# 논문화를 위한 선행 연구 점검 (1차, 2026-09-27)

방법: 웹 검색만 했다. 원문 PDF 접근은 막혀 있다(arXiv, 출판사). "[원문 필요]" 표시가 붙은 항목은 사용자 도움이 필요하다.

## 0. 원문 확인 결과 (2026-09-27 갱신)

- **Choi–Nakajima–Rim, SRC TR 87-203(= SIAM J. Discrete Math. 2, 1989)**: 원문을 확인했다.
    - 정리 1은 최대 차수 3에서 정점 이분화와 간선 이분화의 최소 크기가 같다는 것이다. 간선 집합의 각 간선에서 끝점을 하나씩 고르는 구성이고, **고른 정점들의 독립성은 다루지 않는다.** 반대 방향에서는 오히려 지운 정점 사이 간선을 길이 3 경로로 쪼개 독립성을 인위적으로 만든다.
    - **따라서 추측 C.5는 이 논문에 없다.** 정리 C.1의 되돌리기 논증은 그들의 논증과 같은 계열이므로 인용한다.
    - 정리 2는 평면 Δ = 4에서 정점 이분화가 NP-완전이라는 것이다. **이 구성으로 우리 문제의 Δ = 4, 0 < D3 < 2 NP-완전을 얻었다**(`theorem_D.md`).
- **Žukovič–Bobák, PRE 87, 032121 (2013), arXiv:1212.5447**: 원문을 확인했다.
    - 해밀토니안은 H = −J Σ SᵢSⱼ − D Σ Sᵢ²(J < 0)이고, d = D/|J|다.
    - 삼각 격자 바닥상태는 d > 0이면 Ising(Wannier)형, −3/2 < d < 0이면 (1, −1, 0) 부격자형, d < −3/2이면 모두 0이다.
- **Chen–Li–Wang, "Frustration indices of signed subcubic graphs", arXiv:2511.15226v1 (2025)**: 원문을 확인했다.
    - 주 결과는 좌절 지수의 **상계**다. 2-간선연결 단순 이웃-3 부호 그래프는 예외 5개(Γ̂₁–Γ̂₅)를 빼면 F ≤ n/3이다. 연결 그래프는 (K₄, −)를 빼면 F ≤ (3n+2)/8이고, 등호 그래프도 특징지었다.
    - **독립 횡단(추측 C.5)은 다루지 않는다.** 크기의 상계만 있고, 지울 정점 사이의 이웃 관계는 없다.
    - 우리 fr는 모든 간선이 음인 경우 F(G, −)다. 예외 중 Γ̂₁ = (K₄, −)만 모두-음과 동치임이 확실하다. Γ̂₂는 사이클 L–T–R–B의 부호가 길이 짝수와 어긋나므로 모두-음이 아니다. 어느 경우든 예외들은 8정점 이하라서 우리 전수 검사(≤ 12정점)에 포함된다.
    - **C.5에 쓸 수 있는 따름.** K₄가 아닌 2-간선연결 이웃-3 그래프는 fr ≤ n/3이다(모두-음 예외가 있더라도 전수 검사에서 치유 가능). Brooks 정리로 3색칠 가능하므로 α ≥ n/3이다. 따라서 **필요조건 fr ≤ α가 성립한다.** 추측 C.5의 수 세기 장애는 없다.
    - Lemma 2.2(최적 부호에서 홀수 유도 사이클 C = 2k+1의 경계 절단에 음 간선이 k+1개 이상 있을 수 없음)는 우리 L1(좌절 결합은 매칭)의 일반화다. C.5 증명의 충돌 분석에 쓸 수 있는 국소 도구다.
- (사용자가 보낸 arXiv:1309.3590은 Žukovič–Mižišin–Bobák의 **쌓인 삼각 격자 Ising** 논문으로, BC 논문이 아니다.)

### 정확한 BC 대응 (확인 후 유도)
|J| = 1, S = s⁻¹(3)(BC의 0), 차수 d(v)로 두면

```
H_BC = 2·E_ours − 3·e(S) + Σ_{v∈S} (d(v) + D − 2·D3) + (상수: −|E| − D·N)
```

- **z-정규 격자이고 상태 3이 독립이면 D3 = (z + D)/2에서 두 모형이 같다.** 삼각 격자(z = 6): BC의 d = 0 ↔ D3 = 3, d = −3/2 ↔ D3 = 2.25.
    - Žukovič–Bobák의 Wannier/(1, −1, 0) 경계 d = 0은 우리 D3* = 3과 정확히 일치한다.
    - 그들의 모두-0 상은 우리 모형에 없다. 상태 3끼리의 결합이 1을 내기 때문이다.
- **중요한 함의(정직하게 기록).** BC 언어로는 "d > 0이면 0 상태가 사라진다"가 **거의 자명**하다. BC에서는 ± 반대 결합이 음의 에너지라서, 0을 ±로 바꿀 때 결합 에너지 변화가 min − max ≤ 0이기 때문이다.
    - 따라서 **정리 A는 정규 격자에서는 BC의 자명한 사실을 Potts 정규화로 옮긴 것**이다.
    - 정리 A의 내용은 (i) Potts 모형의 자연스러운 변수 D3로 쓴 격자 의존 문턱 ⌊Δ/2⌋, (ii) 불규칙 격자(차수 의존 항), (iii) 빡빡함, (iv) 국소 구조에 있다. 논문에서는 정리 A를 주 결과로 내세우지 말고 **예비 정리**로 두어야 한다.
    - 새로운 무게는 **정리 C(이웃 3 치유)와 정리 D(Δ = 4 NP-완전)**에 있다.

## 1. 가장 가까운 모형: Blume–Capel (spin-1) 반강자성체

우리 모형(상태 1, 2, 3; 같은 상태 결합 1; 상태 3 비용 D3)은 **Blume–Capel(BC) 반강자성체**와 매우 가깝다. 대응은 상태 1, 2 ↔ s = ±1, 상태 3 ↔ s = 0이다.

| 결합 | 우리 모형 | BC 반강자성 (J Σ s_i s_j) |
|---|---|---|
| ±와 ± (같음) | 1 | +J |
| ±와 ∓ (다름) | 0 | −J |
| ±와 0 | 0 | 0 |
| **0과 0** | **1** | **0** |

J = 1/2로 맞추면 BC 에너지 = (우리 에너지) − |E|/2 + Σ_{v: s=0} deg(v)/2 − (0–0 결합 수)·(3/2) + (단일 이온 항)이다.

- **정규 격자(모든 사이트 이웃 z개)에서는** 차수 항이 상수 이동(D3 ↔ D_BC + z/4, 단위에 따라 다름)이 된다. 남는 차이는 **0–0 결합 에너지**뿐이다.
- 우리 모형의 바닥상태에서 상태 3이 서로 이웃하지 않는 경우에는(Δ ≤ 4의 모든 D3 > 0, 정리 A.5), 두 모형의 바닥상태가 사실상 일치한다. 정확한 조건은 [원문 확인 후 정리 필요]이다.
- **결론:** 논문은 BC 반강자성 문헌과의 관계를 명시적으로 다뤄야 한다. 기여는 "일반 평면 격자 전체에 대한 엄밀한 정리와 정확 알고리즘"으로 규정한다.

### 핵심 선행 연구
- **Žukovič & Bobák, "Phase transitions in a triangular Blume-Capel antiferromagnet", Phys. Rev. E 87, 032121 (2013), arXiv:1212.5447** [원문 필요]
    - 삼각 격자 BC 반강자성체의 바닥상태를 다룬다. 축소 이방성 −1.5 < D/J < 0에서 두 부격자는 반강자성 정렬, 세 번째 부격자는 비자성(s = 0)이다.
    - 경계 1.5는 우리 단위로 D3* = 3에 해당하는 것으로 보인다(결합 에너지 차 2J = 1 환산).
    - **따라서 삼각 격자의 D3* = 3 자체는 (BC 형태로) 알려져 있다.** 우리 쪽 새로운 점은 모든 격자에서의 ⌊Δ/2⌋ 법칙, 빡빡함, 국소 구조다.
- Žukovič, Borovský, Bobák, 선택적 희석 삼각·벌집 Ising 반강자성체, arXiv:1212.5437 (담금질 희석이라 우리와 다르다. 비교용)
- Hartmann & Rieger, *Frustrated systems: ground state properties via combinatorial optimization* (Springer 강의록). 희석 반강자성체, 매칭·흐름 알고리즘을 다룬다 [원문 필요: 빈자리 있는 평면 Ising의 정확 알고리즘이 있는지]

## 2. 그래프 이론 쪽

- **Choi, Nakajima, Rim, "Graph bipartization and via minimization", SIAM J. Discrete Math. 2(1):38–47 (1989)** [원문 필요]
    - 최대 차수 3에서 "정점 이분화 최소 = 간선 이분화 최소(fr)"를 보였다.
    - 확인할 것: 그들의 구성이 이미 **독립** 정점 집합을 주는지. 그렇다면 추측 C.5의 상당 부분이 알려져 있다.
- **Faria, Klein, Stehlík, "Odd cycle transversals and independent sets in fullerene graphs", SIAM J. Discrete Math. 26(3) (2012)** [원문 필요]. 풀러렌(평면 3-정규)에서 홀수 사이클 횡단과 독립집합의 관계를 다룬다.
- **Chen–Li–Wang, arXiv:2511.15226**: 확인 완료(§0). 상계만 다루며 C.5는 없다.
- Brooks(1941), Lovász의 분할 정리 등 고전 결과(인용용)

## 3. 복잡도 쪽 (이미 확인)
- Johnson 외, *Complexity Framework for Forbidden Subgraphs I*
    - arXiv:2211.12887v5는 확인했다(평면성 빈틈 발견).
    - **저널판: Algorithmica 87(3):429–464 (2025)** [원문 필요]. 정리 11의 평면성 논증이 고쳐졌는지 확인해야 한다. 고쳐지지 않았으면 저자들에게 알린다.
- Pilz 2019, Kazda–Kolmogorov–Rolínek, Fulla–Živný, Barahona 1982: 확인 완료

## 4. 추가 검색 키워드 (사용자 검색용)

- Blume-Capel antiferromagnet ground state kagome / honeycomb / diced / Shastry-Sutherland / arbitrary lattice
- spin-1 Ising antiferromagnet crystal field frustration relief nonmagnetic sites "ground-state phase diagram"
- annealed vacancies frustrated Ising exact ground state matching / "site-diluted" planar spin glass chemical potential
- "independent odd cycle transversal" subcubic / "stable bipartization" maximum degree 3
- "3-colouring" "colour class" minimum cubic graph / "bipartite subgraph" 3-colouring subcubic
- Potts antiferromagnet crystal field / single-ion anisotropy ground state

## 5. 지금까지의 판단

| 결과 | 선행 연구 상태 |
|---|---|
| 삼각 격자 D3* = 3 | **BC 형태로 알려짐**(Žukovič–Bobák 2013). 인용하고 재유도로 제시한다 |
| 일반 ⌊Δ/2⌋ 법칙 + 모든 Δ에서의 빡빡함 + 사이트별 형태 | 검색에서 발견 못함 (새로울 가능성 높음) |
| 국소 구조 정리 (독립성, 2:2) | 발견 못함 |
| 정리 C (이웃 3 치유, E = min(D3,1)·fr) | 발견 못함. Choi–Nakajima–Rim과의 관계 확인 필요 |
| 매칭 + 2-SAT 인증 알고리즘 | 발견 못함 |
| 추측 C.5 | 발견 못함. Choi–Nakajima–Rim, Chen–Li–Wang 모두 독립성은 다루지 않음(확인 완료) |
| 펜로즈 쌍대·깎은 펜로즈 결과 | 발견 못함 |

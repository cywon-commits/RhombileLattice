# 인수인계 노트: D3 복잡도 프로젝트 (2026-09-27 기준)

새 대화(일반 채팅 또는 새 Claude Code 세션)에서 이 작업을 이어 가기 위한 요약이다. **이 파일 전체를 새 대화의 첫 메시지로 붙여넣으면 된다.** 소통은 한국어로 하고, 제출용 원고는 영문이다.

---

## 0. 어디에 무엇이 있나

- **저장소:** GitHub `cywon-commits/RhombileLattice`
    - 작업 브랜치: **`claude/d3-complexity`** (여기에만 push한다)
    - 폴더: `d3_complexity/`
    - https://github.com/cywon-commits/RhombileLattice/tree/claude/d3-complexity/d3_complexity
- **공유 작업 문서(Claude Docs):** "D3 논문 원고", https://claude.ai/code/artifact/921480cd-99f1-4005-a8ba-63659fdec2a1
    - 한국어 작업본이다. 초록, 결과 표, 1–8절, 남은 일 체크리스트, 참고문헌이 들어 있다.
    - 규칙: 내용은 문서에서 먼저 고치고, 영문 `paper/main.tex`에 옮긴다.
- **해설 페이지(초기 산출물):** `report/d3_explainer.html`. 일반 독자용 그림 해설이며 정리 D 이전 판이다.
- **작업 제약(처음부터 유지):**
    - 부모 브랜치 `claude/modest-feynman-dtkpqn`에는 push하지 않는다. 다른 세션들이 그 브랜치를 쓴다.
    - 저장소의 Part I–IX 물리 내용은 범위 밖이다(COMPLEXITY.md에 인용된 빡빡함 예시만 예외).
    - 변경은 모두 커밋하고 push한다.

## 1. 문제

평면 격자(그래프) G에서 3상태 배치 s : V → {1, 2, 3}의 에너지를 다음과 같이 둔다.

```
E(s) = #{같은 상태 결합} + D3 · #{상태 3 사이트}
```

이는 이방성 (0, 0, D3)가 있는 반강자성 3상태 Potts 모형이다. 바닥상태 문제의 계산 복잡도가 **최대 이웃 수 Δ와 D3**에 따라 어떻게 바뀌는지 묻는다.

기호:
- fr(G): 2상태(Ising) 바닥상태의 좌절 결합 수. |E| − maxcut과 같다.
- ioct(G): 지우면 이분 그래프가 되는 **독립** 정점 집합의 최소 크기.
- 재정식화: OPT(D3) = min_S [D3|S| + e(S) + fr(G − S)]. OPT는 D3에 대해 오목한 조각별 선형 함수다.

## 2. 결과와 상태

| 결과 | 내용 | 상태 | 근거 |
|---|---|---|---|
| 보조정리(국소 재색칠) | 이웃 d개인 상태 3 사이트는 좌절을 최대 ⌊d/2⌋개 푼다 | 증명 | `theorem_A.md` |
| 정리 A(문턱) | D3 ≥ ⌊Δ/2⌋이면 OPT = fr, 평면이면 매칭으로 다항 | 증명(**예비 결과로 격하**) | `theorem_A.md` |
| 명제(빡빡함) | 바퀴 W_Δ에서 모든 Δ ≥ 3에 대해 문턱이 달성된다 | 증명 | `theorem_A.md` |
| 정리 A.5(국소 구조) | Δ ≤ 4면 상태 3이 고립된다. Δ = 4, 1 < D3 < 2면 이웃이 2:2다 | 증명(**(iv) 정정됨**, 아래 5절) | `theorem_A.md` §4, `main.tex` |
| BC 정확 대응 | H_BC = 2E − 3e(S) + Σ_{v∈S}(d(v) + D − 2D3) − \|E\| − D\|V\|. z-정규이고 상태 3이 고립이면 D3 = (z + D)/2 | 유도, Žukovič–Bobák과 일치 | `lit_review.md` §0 |
| 정리 C.1(샌드위치) | Δ ≤ 3, 0 ≤ D3 ≤ 1: D3·fr ≤ OPT ≤ min(fr, D3·ioct) | 증명 | `theorem_C.md` |
| 치유 판정(C.2) | OPT = D3·fr ⇔ ioct = fr ⇔ 어떤 Ising 바닥상태의 좌절 매칭에 독립 횡단이 있다 | 증명 | `theorem_C.md` |
| **정리 C.5(치유 정리)** | **K₄ 성분이 없는 Δ ≤ 3 그래프는 ioct = fr.** 따라서 모든 D3에서 OPT = min(D3, 1)·fr. 평면성 불필요 | **증명 초안, 적대적 검토 2회 통과** | `theorem_C.md` §C.5b, `paper/healing_proof.tex`, `reviews/c5b_review.md`, `reviews/c5b_review2.md` |
| **정리 D(Δ = 4 어려움)** | 평면 Δ = 4, 고정 0 < D3 < 2에서 NP-완전(D3 = 0은 3색칠) | 증명, 평면성은 계산으로 확인 | `theorem_D.md`, `reviews/cnr_planarity.md` |
| 인증 알고리즘 | 평면 T-join(매칭) + 2-SAT으로 바닥상태를 구성하고 최적성을 인증 | 구현, 15개 격자 전부 성공 | `heal_check.py`, `ising_tjoin.py`, `timing_cert.py` |

**복잡도 지형(평면):**

| Δ | D3 = 0 | 0 < D3 < 1 | 1 ≤ D3 < 2 | 2 ≤ D3 < ⌊Δ/2⌋ | D3 ≥ ⌊Δ/2⌋ |
|---|---|---|---|---|---|
| 3 | 다항(Brooks) | 다항(C.5) | 다항 | 해당 없음 | 다항 |
| 4, 5 | NP-완전 | NP-완전(D) | NP-완전(D) | 해당 없음 | 다항 |
| ≥ 6 | NP-완전 | NP-완전 | NP-완전 | **열림** | 다항 |

**논문 구도:** "이웃 3은 모든 D3에서 쉽고, 이웃 4는 D3 = 2에서 정확히 갈린다. 평면 Δ ≤ 5는 완전히 분류된다."
- 가제: *Coordination three heals, coordination four is hard: the antiferromagnetic Potts model with a costly third state*

## 3. 핵심 증명의 뼈대

### 정리 C.5 (치유 정리)
1. φ(X) = |X| + fr(G − X) ≥ fr이다(되살리기: 각 정점을 넣을 때 좌절이 1 이하로 는다). 독립 S가 φ(S) = fr이면 "빡빡하다". 크기가 최대인 빡빡한 S가 k = fr(G − S) = 0임을 보이면 된다. k ≥ 1이라 가정한다.
2. **덧붙이기:** G − S 최대 컷의 좌절 결합 끝점이 S와 이웃하지 않으면 S에 더해도 빡빡하다. 이는 최대성에 어긋난다.
3. **막힌 구조:** 좌절이 끝점 사이를 미끄러지므로, 좌절 있는 성분은 홀수 사이클이다. 각 정점은 차수 3이고 S-이웃이 정확히 하나다.
4. **교환:** 사이클 정점 y ↔ S-이웃 t로 바꿔도 최대 빡빡 집합이다. 그러면 t는 사이클에 이웃이 하나뿐이고, 나머지 두 이웃은 짝수 경로 성분 Q의 끝점이며, Q + t가 새 홀수 사이클이 된다. 예외는 K₄뿐이다.
5. **사슬:** 교환을 반복하면 새 경로 R_j는 항상 원래 G − S의 새로운 성분이다(회전 경로 배제는 경우 나누기로 보인다).
6. **유한성:** 사슬은 끝없이 늘릴 수 있으므로, 유한 그래프에서 모순이다. (처음 초안은 사적 결합 세기로 |Σ| ≥ 3 + 2|Σ| 모순을 냈다. 이것도 옳지만 더 길다.)

- 검토 결과(`reviews/c5b_review.md`): 치명적 오류가 없다. 결론 단순화와 표현 보완을 제안받아 반영했다.
- 2차 검토(`reviews/c5b_review2.md`, 2026-09-27): 판정 "건전". 회전 경로 배제의 경우 나누기는 빠짐이 없고, theorem_C.md 227행의 의문은 N(p_i) ⊆ R_{i−1} ∪ {w_{i−1}, w_i}로 해소된다. 사소한 보완 4개(K₄ 단계의 최대 빡빡 집합 명시, D − t가 한 성분인 이유, G − T_j 성분 목록 (i)–(iii)의 근거, j = 0 기저)와 표현 문제 몇 개가 **아직 `healing_proof.tex`에 반영되지 않았다**. 수치: n ≤ 11 전수 8,093개에서 중간 보조정리 위반 0, 무작위 사슬 구조 시험 80,000그래프·456k 단계에서 위반 0, MILP 1,386개(16–100정점, 평면 729)에서 불일치 0.
- 이 정리는 Johnson 외(Algorithmica 2025) 정리 11의 "평면 이웃-3 IOCT는 NP-완전" 주장과 양립할 수 없다(P ≠ NP 가정). 그 논문의 평면성 논증에 빈틈이 있음을 이미 찾았고, 저널판에도 그대로 있다.
- **투고 전 사람의 외부 검토를 권한다.**

### 정리 D (Δ = 4, 0 < D3 < 2 NP-완전)
- Choi–Nakajima–Rim(1989) 정리 2의 그래프 G_φ를 쓴다.
    - 변수마다 길이 4n 사이클이 있고, 각 변 위에 삼각형 꼭짓점을 얹는다.
    - 절마다 w₁…w₄ 경로가 있고, 절 사이클 길이는 13이다.
    - 크기는 28m 정점, 43m 간선이다.
- 삼각형별 비용 나눠 매기기로 보인다. 상태 3 밑변 b는 D3/2씩, 꼭짓점 a는 D3를 매긴다. D3 < 2이므로 삼각형마다 ≥ D3/2이고, 12m개이므로 OPT ≥ 6m·D3다. 등호 ⇔ φ가 충족 가능이다.
- **평면성 조건:** 변수 사이클의 출현 순서는 결합 그래프 평면 임베딩의 회전 순서를 따르고, 절의 리터럴 순서는 절 노드 회전의 **반대 방향**이어야 한다. 평면 3-SAT 307개가 이 조건에서 전부 평면이었다(같은 방향이면 67개).

## 4. 수치 근거 요약

- **C.5 전수 검사:** 연결 Δ ≤ 3 그래프 12정점까지 27,523개 전부에서 성립(K₄만 예외)(`exhaustive_heal.py`, `hub_runs/exhaustive_12.log`).
- **C.5 추가 검증(검토 에이전트):**
    - 9정점 이하 836개에서 보조정리 1과 (F1) 위반 0.
    - MILP로 큰 그래프 251개(14–60정점)에서 ioct = fr 불일치 0(`reviews/tight_exhaustive.py`, `reviews/milp_stress.py`).
- **탐욕적 덧붙이기:** 12–20정점 972회에서 한 번도 막히지 않았다(`greedy_heal.py`, `hub_runs/greedy_12_20.log`).
- **강한 판("모든 최대 컷이 치유 가능")은 거짓:** 18정점 반례(`hub_runs/strong_version_counterexample.json`).
- **정리 D:**
    - 무작위 충족 가능 48/48에서 등호.
    - 불충족 인스턴스(224사이트)에서 OPT − 6mD3 = min(D3, 2 − D3).
    - 평면 인스턴스 11개도 재확인(`cnr_reduction.py`, `cnr_unsat_test.py`, `cnr_planarity_check.py`).
- **인증 알고리즘 시간:**
    - 약 4,000사이트에서 3–35초, 최대 7,974사이트에서 340초. 시간은 대부분 매칭(networkx)이 차지한다.
    - 2-SAT은 0.01초 이하.
    - 결과: `hub_runs/timing_cert.jsonl`, 원고 표 `tab:timing`.
- **펜로즈 면-인접 그래프(Δ = 4):**
    - D3 = 1, 2의 꺾임은 원리적이다.
    - 중간 꺾임(4/3, 5/4, 7/5, 3/2, 5/3)은 가장자리에 갇혀 있지 않다.
    - 그러나 사이트당 이득이 0.004–0.024로 흔들린다. 열역학 극한은 미정이고 주기 근사체가 필요하다(`penrose_runs/SUMMARY.md`).

## 5. 철회·정정된 것 (같은 실수를 반복하지 말 것)

- **정리 B1(Δ = 3, 0 < D3 < 1 NP-완전) 철회.**
    - 인용한 Johnson 외 구성이 어려운 인스턴스에서 평면이 아니다(`theorem_delta3.md`, `planarity_gap.py`).
    - 정리 C.5가 맞다면 이 경우는 오히려 다항이다.
- **정리 A.5(iv) 정정.** "Δ = 3, 0 < D3 < 1이면 상태 3 사이트는 차수 3"은 틀렸다.
    - 차수 2 사이트도 가능하다(삼각형 + 가지 하나 반례).
    - 올바른 진술: a = 0, n₁, n₂ ≥ 1, d ∈ {2, 3}.
- **펜로즈 "벌크 효과" 결론 일부 철회.** 사이트당 이득이 크기와 무관하다는 주장은 r8(248 마름모)에서 0.004가 나와 철회했다.
- **정리 A의 무게.** BC 언어로는 거의 자명하므로 예비 결과로 낮췄다.
- **작업 교훈:**
    - arXiv 등 외부 사이트는 프록시로 막혀 있다. 논문은 사용자가 PDF를 올려 준다(pymupdf로 읽음).
    - `pkill -f`에 자기 명령줄에 들어가는 패턴을 쓰지 말 것.

## 6. 문헌 점검 결과 (`lit_review.md`)

- **Žukovič–Bobák, PRE 87, 032121 (2013):** 삼각 격자 BC 경계 D = 0이 우리 D3* = 3과 일치한다.
- **Choi–Nakajima–Rim, SIAM J. Discrete Math. 2 (1989):** 정리 1(Δ ≤ 3에서 OCT = fr)에는 독립성이 없다. 정리 2가 정리 D의 구성이다.
- **Chen–Li–Wang, arXiv:2511.15226:** 부호 이웃-3 그래프 좌절 지수 상계를 다룬다. C.5는 없다. fr ≤ n/3 ≤ α와 일치한다.
- **Johnson 외, Algorithmica 87:429–464 (2025), doi:10.1007/s00453-024-01289-2:** 정리 11의 평면성 논증은 v5와 동일하고, 빈틈이 남아 있다. 사용자는 저자에게 알리지 않기로 했다.
- **확인 완료:** Pilz 2019, Kazda–Kolmogorov–Rolínek, Fulla–Živný, Barahona 1982.

## 7. 원고 상태 (`paper/`)

- `main.tex`(revtex4-2, PRE):
    - 새 제목과 초록.
    - 서론의 결과 목록은 (i) 이웃 3 쉬움, (ii) 이웃 4 어려움, (iii) 국소 경계, (iv) 인증 알고리즘이다.
    - 절 구성: 모형·BC 대응, 문턱, 국소 구조(증명 포함), 정확 알고리즘, 이웃 3 치유(`\input{healing_proof}`), 이웃 4 어려움, 격자, 논의, 부록(수치 방법·시간표).
- 그림 7개(`paper/figs/*.pdf`, `paper/make_figs.py`로 재생성): 재색칠, 바퀴, 삼각 격자, 펜로즈, 치유, 깎기 구성, 시간.
- **2026-09-27 컴파일 성공**(TeX Live + `texlive-publishers`로 revtex4-2 설치, `latexmk -pdf main.tex`): 9쪽, 오류·미정의 참조 0, overfull 박스 2개(식 (1), 표 `tab:timing`)를 고쳐 0개.
- 남은 표시: `[TODO authors]`, `[TODO affiliation]`.

## 7a. 보충자료 (2026-09-28)

- `paper/supplementary.tex`(Note A): Johnson 외 정리 11(평면 이웃-3 IOCT NP-완전)과의 관계를 자세히 기술. 본문·초록은 "refute" 등 강한 표현을 빼고 보충자료를 가리키도록 완화했다.
- 내용: 합의점(환원 동치성, 일반 이웃-3 NP-완전성) / 쓰는 구조 (F0)–(F3) / 보조정리: H(φ) = 연결 그래프 + 모든 절에 붙은 꼭짓점이 G_φ의 minor / 가장 작은 예 φ* = (x∨y)(x̄∨y)(x∨ȳ)에서 H = K₃,₃ / 변형별 사례 표(절 크기 3, 변수 4회, 불충족 사례, Lichtenstein 변수 사이클 조건 포함) / 무작위 131개 중 130개 비평면 / 구조적 이유(co-nested 식은 다항).
- 근거 스크립트: `supp_planarity.py`(출력 `supp_runs/planarity_evidence.json`), 그림 `paper/make_supp_figs.py`.
- 적대적 검토: `reviews/supp_s1_review.md`(치명 0, 주요 3건 반영: ioct ≥ fr의 차수 조건, co-nested 논증 정밀화, Tovey 부류 사례 교체).
- **남은 확인(빨간 [VERIFY] 표시):** 정리 11 원문 인용, 구성 (F0)–(F3)과 그림 3 대조, 쓰인 평면 3-SAT 변형, Kratochvíl–Křivánek·Pilz 진술, Yannakakis 인용, 저장소 URL. **저널판 PDF가 필요하다.**

## 8. 남은 일 (우선순위 순)

1. ~~컴파일~~ 완료(2026-09-27). 남은 것은 내용 교정.
2. 저자·소속·투고 저널을 정한다(현재 PRE 형식).
3. 정리 C.5에 대한 사람의 외부 검토(출판된 주장과 배치되는 결과이므로). 에이전트 2차 검토는 완료. 그 사소한 보완을 `healing_proof.tex`에 반영하는 일이 남았다.
4. 공유 문서와 `main.tex` 동기화 점검(문서가 한국어 기준본).
5. 펜로즈 주기 근사체로 중간 꺾임의 열역학 극한을 확인한다.
6. (선택) 치유 정리의 구성적 판: 탐욕적 덧붙이기가 항상 성공하는지.
7. (선택) Δ ≥ 6, 2 ≤ D3 < ⌊Δ/2⌋ 구간. 남은 유일한 열린 문제이고 후속 연구감이다.
8. (선택) 해설 페이지 `report/d3_explainer.html`에 정리 C.5와 정리 D를 반영한다.

## 9. 코드 빠른 안내 (`d3_complexity/`)

| 파일 | 용도 |
|---|---|
| `potts_exact.py` | `solve(n, edges, D3)` MILP 정확해(HiGHS), `curve()` 에너지 포락선 |
| `ising_tjoin.py`, `heal_check.py` | 평면 T-join Ising 바닥상태, 2-SAT 치유, `healed_state()` |
| `lattices3.py`, `penrose.py` | 이웃-3 격자 생성기(깎은 펜로즈, Voronoi 거품, 결함 벌집), 펜로즈 타일링 |
| `exhaustive_heal.py`, `greedy_heal.py`, `stuck_test.py` | C.5 전수 검사, 탐욕법, 막힌 집합 시험 |
| `cnr_reduction.py`, `cnr_unsat_test.py`, `cnr_planarity_check.py` | 정리 D 구성과 검증 |
| `penrose_boundary.py` | 펜로즈 조각 꺾임 분석 |
| `timing_cert.py` | 인증 알고리즘 시간 측정 |
| `paper/make_figs.py` | 원고 그림과 시간표 생성(`paper/`에서 실행) |

필요 패키지: `requirements.txt`(networkx, numpy, scipy, matplotlib; PDF 읽기용 pymupdf는 선택). 논문 항목별 재현 명령은 `README.md`의 대응표에 있다. `reviews/`의 스크립트는 `PYTHONPATH=.`로 실행한다.

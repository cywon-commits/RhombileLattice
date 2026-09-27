# 원고 뼈대 (영문, revtex4-2 / Phys. Rev. E 형식)

> 함께 보고 고치는 한국어 작업본은 Claude Docs 문서 "D3 논문 원고"에 있다: https://claude.ai/code/artifact/921480cd-99f1-4005-a8ba-63659fdec2a1 . 내용 수정은 그 문서에서 먼저 하고, 제출용 영문은 이 `main.tex`에 반영한다.

- `main.tex`: 원고 본문. 정리 문장과 짧은 증명은 들어 있다. 서론·논의·격자 절 등은 `[TODO]`로 남겨 두었다.
- `refs.bib`: 참고문헌. `note = {check}` 항목은 원문 대조가 필요하다.
- 컴파일: `latexmk -pdf main.tex` (TeX Live + REVTeX, 예: `texlive-publishers`). 2026-09-27 컨테이너에서 경고 없이 컴파일됨.

## 절 구성과 근거 문서

| 절 | 내용 | 근거 (`d3_complexity/`) |
|---|---|---|
| 1 서론 | 좌절, 세 번째 상태, Barahona, Blume–Capel, 결과 요약 | `D3_RESEARCH_PLAN.md`, `lit_review.md` |
| 2 모형 | 해밀토니안, 정확한 BC 대응, 재정식화 | `theorem_A.md` §4 |
| 3 문턱 | 보조정리(국소 재색칠), 정리(문턱), 명제(바퀴), 삼각 격자 | `theorem_A.md` §1–3, §5 |
| 4 국소 구조 | 정리(독립성, 2:2, 구간 구조) | `theorem_A.md` §4 |
| 5 정확 알고리즘 | 평면 max-cut / T-join | `theorem_A.md` §3 |
| 6 이웃 3 치유 | 샌드위치, 치유 판정, 인증 알고리즘, 치유 정리와 증명(`healing_proof.tex`) | `theorem_C.md` |
| 7 Δ = 4 어려움 | 정리 D (CNR 구성 + 삼각형 나눠 매기기) | `theorem_D.md`, `cnr_reduction.py` |
| 8 격자 | 삼각, 펜로즈 쌍대, 깎기, 깎은 펜로즈·거품·결함 벌집 | `report/`, `penrose_scan.py`, `lattices3.py` |
| 9 논의 | 복잡도 지형(짧게), 대칭 대 장, 열린 문제 | `perfect_colouring_obstruction.md`, `theorem_delta3.md` |
| 부록 A | 수치 방법, 전수 생성, 강한 판 반례 | `potts_exact.py`, `exhaustive_heal.py`, `hub_runs/` |

## 그림 목록

`python3 make_figs.py` (이 폴더에서)가 `figs/*.pdf`(matplotlib 벡터)를 만들고, `main.tex` 부록의 시간표(`% BEGIN generated timing table` 표지 사이)를 `../hub_runs/timing_cert.jsonl`로 다시 쓴다. `python3 make_figs.py penrose timing`처럼 일부만 만들 수도 있다. 색: 상태 1 파랑 원, 상태 2 주황 사각형, 상태 3 노랑 삼각형(타일은 빗금), 좌절 결합 굵은 빨강. 작은 그래프에는 상태 숫자를 넣어 흑백에서도 읽힌다.

| 파일 | label | 위치 | 내용 | 원천 |
|---|---|---|---|---|
| `figs/fig_recolour.pdf` | `fig:recolour` | §3 보조정리 뒤 | (a) 삼각형 좌절과 세 번째 상태, (b) d = 4 국소 재색칠: 2:2 대 1:3 | `report/gen_svgs.py` `triangle`, `star_*` |
| `figs/fig_wheel.pdf` | `fig:wheel` | §3 명제(바퀴) 뒤 | W₄, W₅ 바닥상태와 D3 = 2 교차 | `w*_*`, `chart_w5` |
| `figs/fig_triangular.pdf` | `fig:triangular` | §8 삼각 격자 | 세 부격자 배치, Wannier 줄무늬, 사이트당 에너지 | `tri3`, `tri2`, `chart_tri` |
| `figs/fig_penrose.pdf` (figure*) | `fig:penrose` | §8 펜로즈 쌍대 | 140 마름모 조각의 바닥상태(D3 = 0.5, 1.5, 2.5)와 정확한 곡선 + 다른 조각 9개 | `report/penrose_states.json`, `penrose_runs/r*_s*.json` |
| `figs/fig_healing.pdf` | `fig:healing` | §6 인증 알고리즘 | (a) 삼각 기둥 치유 3단계, (b) 깎은 펜로즈(N = 270) 곡선과 MILP 점 | `ising_tjoin.py`, `heal_check.py`, `potts_exact.py` |
| `figs/fig_truncation.pdf` (figure*) | `fig:truncation` | §8 깎은 격자 | 깎기 구성, 깎은 펜로즈·거품·별 격자의 인증 바닥상태 | `lattices3.py`, `report/deg3_gallery.py` |
| `figs/fig_timing.pdf` | `fig:timing` | 부록 A 시간 | 인증 알고리즘 시간(매칭, 2-SAT) 대 사이트 수 | `hub_runs/timing_cert.jsonl` |
| 표 `tab:timing` | `tab:timing` | 부록 A 시간 | 같은 자료의 표 (make_figs.py가 생성) | `hub_runs/timing_cert.jsonl` |

계획에 있던 "국소 구조: 허용·금지 모양"(`loc_*`) 그림은 넣지 않았다(정리 문장으로 충분하다고 보고 뺌; 필요하면 추가).

## 남은 일
1. ~~[TODO-BC]~~ 완료: 정확한 BC 대응(식 BC)과 Žukovič–Bobák 확인을 반영했다.
2. ~~Choi–Nakajima–Rim 확인~~ 완료: 추측 C.5는 없고, 정리 2로 정리 D(새 절)를 얻었다. 부호 그래프 문헌과 Johnson 외 저널판은 남았다.
3. ~~[TODO-NUM] 펜로즈 가장자리 효과~~ 완료: 중간 꺾임은 벌크 재배열이다(`penrose_runs/SUMMARY.md`). 주기 근사체는 남았다.
4. 서론·격자·논의 절 본문과 부록을 쓴다.
5. (병행) 추측 C.5의 평면판 증명을 시도한다.

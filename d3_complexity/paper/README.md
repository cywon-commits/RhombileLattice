# 원고 뼈대 (영문, revtex4-2 / Phys. Rev. E 형식)

- `main.tex`: 원고 본문. 정리 문장과 짧은 증명은 들어 있다. 서론·논의·격자 절 등은 `[TODO]`로 남겨 두었다.
- `refs.bib`: 참고문헌. `note = {check}` 항목은 원문 대조가 필요하다.
- 컴파일: 이 컨테이너에는 TeX이 없다. 로컬에서 `latexmk -pdf main.tex`로 컴파일한다.

## 절 구성과 근거 문서

| 절 | 내용 | 근거 (`d3_complexity/`) |
|---|---|---|
| 1 서론 | 좌절, 세 번째 상태, Barahona, Blume–Capel, 결과 요약 | `D3_RESEARCH_PLAN.md`, `lit_review.md` |
| 2 모형 | 해밀토니안, 정확한 BC 대응, 재정식화 | `theorem_A.md` §4 |
| 3 문턱 | 보조정리(국소 재색칠), 정리(문턱), 명제(바퀴), 삼각 격자 | `theorem_A.md` §1–3, §5 |
| 4 국소 구조 | 정리(독립성, 2:2, 구간 구조) | `theorem_A.md` §4 |
| 5 정확 알고리즘 | 평면 max-cut / T-join | `theorem_A.md` §3 |
| 6 이웃 3 치유 | 샌드위치, 치유 정리, 허브 완화, 인증 알고리즘, 추측 | `theorem_C.md` |
| 7 Δ = 4 어려움 | 정리 D (CNR 구성 + 삼각형 나눠 매기기) | `theorem_D.md`, `cnr_reduction.py` |
| 8 격자 | 삼각, 펜로즈 쌍대, 깎기, 깎은 펜로즈·거품·결함 벌집 | `report/`, `penrose_scan.py`, `lattices3.py` |
| 9 논의 | 복잡도 지형(짧게), 대칭 대 장, 열린 문제 | `perfect_colouring_obstruction.md`, `theorem_delta3.md` |
| 부록 A | 수치 방법, 전수 생성, 강한 판 반례 | `potts_exact.py`, `exhaustive_heal.py`, `hub_runs/` |

## 그림 목록 (계획)

| 번호 | 내용 | 원천 |
|---|---|---|
| 1 | 삼각형 좌절과 세 번째 상태 (개념도) | `report/gen_svgs.py` `triangle` |
| 2 | 국소 재색칠: 이웃 반반 대 치우침 | `star_*` |
| 3 | 바퀴 W₄, W₅와 에너지 교차 | `w*_*`, `chart_w5` |
| 4 | 삼각 격자: 세 색 배치, Wannier 배치, 에너지 | `tri3`, `tri2`, `chart_tri` |
| 5 | 국소 구조: 허용·금지 모양 | `loc_*` |
| 6 | 펜로즈 쌍대: 바닥상태 세 장면과 정확한 에너지 곡선 [TODO-NUM 더 큰 조각] | `pen*`, `chart_pen` |
| 7 | 이웃 3 치유: 삼각 기둥, 깎은 펜로즈 에너지 곡선 | `heal_*`, `chart_tpen` |
| 8 | 깎기 구성과 갤러리(깎은 펜로즈, 거품, 별 격자) | `trunc_*`, `gal_*` |
| 9 | 인증 알고리즘 시간표 [TODO-FIG] | `heal_check.py` |

그림은 현재 SVG다. 원고용으로는 벡터 PDF로 변환하거나 matplotlib으로 다시 그린다.

## 남은 일
1. ~~[TODO-BC]~~ 완료: 정확한 BC 대응(식 BC)과 Žukovič–Bobák 확인을 반영했다.
2. ~~Choi–Nakajima–Rim 확인~~ 완료: 추측 C.5는 없고, 정리 2로 정리 D(새 절)를 얻었다. 부호 그래프 문헌과 Johnson 외 저널판은 남았다.
3. ~~[TODO-NUM] 펜로즈 가장자리 효과~~ 완료: 중간 꺾임은 벌크 재배열이다(`penrose_runs/SUMMARY.md`). 주기 근사체는 남았다.
4. 서론·격자·논의 절 본문과 부록을 쓴다.
5. (병행) 추측 C.5의 평면판 증명을 시도한다.

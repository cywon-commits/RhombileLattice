# Part IX 인계 문서: annealed 결합 삼각격자 반강자성체 (PRE 논문 준비)

작성 2026-09-29. 이 폴더 하나로 논문 작업을 이어갈 수 있게 정리했다. 계산 세션(Claude Code, 저장소 `cywon-commits/RhombileLattice`, 브랜치 `claude/modest-feynman-dtkpqn`, 마지막 커밋 `af18c72`)에서 넘어온 자료다.

**읽는 순서:** 이 문서 → `paper_draft.md`(현재 초안) → `data/`(수치 근거) → `notes/`(모형·문헌 메모).

**범위:** Part IX만 다룬다. 저장소의 Part I–VIII(quenched 결함, matching principle, D3=∞ 정확해, flux sector)은 이 논문에 넣지 않는다.

---

## 1. 한 줄 요지

> 삼각격자 이징 반강자성체(TAFM)의 결합을 annealed로 풀어 주면, 스핀 결함 하나의 에너지 비용이 정확히 두 배가 된다. 그래서 같은 조건에서 상관길이가 지수적으로 길어지고 질서가 더 높은 온도까지 유지되지만, 보편성 부류(3-state Potts)는 바뀌지 않는다. 결함 fugacity는 μ로 조절된다.

## 2. 모형과 단위

- 배위는 두 가지로 이뤄진다. 삼각형 인접 그래프(벌집 격자)의 matching과, 삼각격자 각 점의 3-state Potts 스핀이다.
- 짝지어진 두 삼각형은 공유 결합 하나를 끈다(비활성, 마름모 하나). 짝이 없는 삼각형은 세 결합이 모두 살아 있는 "좌절 삼각형(monomer)"이다. monomer 수는 제약하지 않는다.
- 해밀토니안: `H = Σ_active δ(s_i,s_j) + D2·n2 + D3·n3 + μ·n_mon`.
- 표준 설정은 D1=0, D3=10이다. 세 번째 상태가 얼어서 사실상 이징(0/1)이 된다.
- 계는 L×L 토러스이고, N_site=3L², N_tri=6L², N_bond=9L²이다.
- TAFM과의 대응: 같은 스핀끼리의 활성 결합 하나의 비용이 1이므로 이징 J=1/2, 그리고 **h/J = D2**다.

## 3. 확정된 결과 (논문 본문에 쓸 수 있는 것)

| # | 결과 | 수치 | 근거 파일 |
|---|---|---|---|
| 1 | D1=D2일 때 T=0 앙상블 = TAFM 바닥상태 앙상블 (해석적) | S(0)/N_tri = 0.16153 (Wannier 값) | notes/PART9.md, notes/PART9_literature.md |
| 2 | 0<D2<6일 때 바닥상태 = 자기장 속 TAFM의 √3×√3 상 (해석적). 바닥상태는 유일(허브 부격자 선택 3가지), E/N_tri = D2/6 | 존재 범위 0<D2<6 = TAFM의 0<h<6J | notes/PART9.md |
| 3 | T=∞ 정확값과 열역학 적분이 서로 맞음 | 좌절 비율 0.399, s(∞) = ½ln3 + 0.582 = 1.131. D2=0에서 S(0.06)=0.150 (정확값 0.1615) | notes/PART9.md |
| 4 | **결함 비용 2배 (해석적, 수치로 확인).** 결함 삼각형은 결합을 하나만 가릴 수 있어서 같은 스핀 결합이 2개 이상 드러나고, 결함 쌍의 비용이 2에서 4 이상이 된다 | 활성화 에너지 annealed 3.0–3.8, TAFM 1.5–1.6 | data/tafm_defects_L24_s1.txt, data/mechanism_summary.txt |
| 5 | 계산 절차 검증: 고정 결합 TAFM의 문헌 값 재현 | h/J=1: T_c/J = 0.9426(20), 문헌 0.9458. h/J=3: 1.3436(20), CTMRG 1.3440 | data/tafm_fss_summary.txt |
| 6 | D1=D2에는 유한 온도 전이가 없음 (L=8–32, T≥0.34, PT 시드 4개) | C_max/N_tri 0.647 → 0.636 (L=8→32), ξ/L와 U가 교차하지 않음, ξ ∝ exp(2/T) (기울기 1.98–1.99). TAFM은 exp(1/T) | data/d1d2_summary.txt, data/tafm_xi_L32_s1.txt |
| 7 | 좌절 삼각형은 많지만(∝ exp(−1/T)) 스핀 상관을 끊지 못함: 비만족 결합 string에 묶임 | T=0.46에서 0.126, T=0.34에서 0.055 | data/d1d2_summary.txt |
| 8 | μ=0, D2=1의 결정화는 **연속 3-state Potts 전이** | T_c = 0.528(2). P(E)는 봉우리 하나. 지수는 Potts 값과 맞음. 같은 절차로 본 TAFM(h/J=1)은 0.4713(10) → annealed가 약 12% 높음 | data/fss_d21_summary.txt, paper_draft.md 표 |
| 9 | T_c가 μ와 함께 단조 증가 | μ = 0 / 0.5 / 1 / 2 / 3 / 5 / 10 → 0.528 / 0.620 / 0.696 / ≈0.80 / ≈0.95 / 1.03–1.07 / 1.141(3). μ≥2 값은 ±0.02–0.03 | data/fss_mu_summary.txt, data/fss_mu_spin_summary.txt |
| 10 | μ≤2에서는 T_c 아래가 보통의 장거리 질서 | T/T_c≈0.85–0.95에서 η≈0, U → 0.5 | data/fss_mu_spin_summary.txt |
| 11 | μ≥3에서는 T_c 바로 아래에 멱법칙 감쇠 구간이 보임 (**관측 사실만**) | ψ_s² ∝ L^−η (L=32–96), η 0.1–0.4, U≈0.46 크기 무관. 구간 폭은 μ와 함께 넓어짐. μ=10도 T/T_c≈0.6에서는 장거리 질서 | data/seeds_summary.txt |
| 12 | Z6 대칭 임계상은 배제됨 (μ=10) | ⟨cos 6θ⟩ = 0.1–0.2, ⟨cos 3θ⟩ = 0.3–0.45, L=24→96에서 줄지 않음. Z6 임계상이면 L^−5 정도로 사라져야 함 | data/clock_mu10.txt |
| 13 | 비열 분해 C_n = C_μ − β²cov(E,n)²/var(n) | C_μ의 약 절반이 좌절 수 요동. corr(E,n)=0.71 (μ=0) | paper_draft.md, notes/PART9.md |
| 14 | 방법론: 시드끼리의 차이가 자기상관 오차보다 최대 약 2배 큼 (χ²/dof 최대 4) | 시드 1개짜리 결론은 신뢰도가 낮음 | data/seeds_summary.txt |

## 4. 철회한 주장 (초안이나 예전 메모에 남아 있으면 지울 것)

- "1<μ_c<2": T_c에 너무 가까운 점을 잘못 해석한 것.
- "μ=3은 낮은 T에서 장거리 질서": L=96을 더하자 사라짐.
- "μ_c≈3 다중임계점", "Z6 clock형 중간상": 결과 12로 근거가 약해짐.
- "h/J=1 전이가 1차처럼 보인다": 평형화 부족으로 생긴 착시.
- **paper_draft.md 초록의 "모든 경우의 지수는 3-state Potts와 일치"는 과장이다.** 확인된 것은 μ=0(과 고정 결합 TAFM)뿐이다. μ>0 지수는 L≤48 수준이라 오차가 크다. 초록을 고칠 것.

## 5. 미확정 (주장하지 말고 열린 문제로 둘 것)

1. 큰 μ의 멱법칙 구간이 진짜 상인지 크로스오버인지.
2. 큰 μ 전이가 Potts형인지 KT형인지.
   - **작업 가설:** Z3 비등방성은 η<4/9에서 relevant다. 그러면 이 구간은 결국 질서상으로 가는 긴 크로스오버다.
   - μ=10의 유효 η는 T=1.04–1.10에서 0.39–0.40으로 머문다. 이는 Z3 잠김 값 4/9에 가깝고 Potts 값 4/15보다 크다.
   - 결함(소용돌이)은 η=4/9에서도 relevant하다(스케일링 차원 1/(2η) = 9/8 < 2). 따라서 아주 큰 L에서는 Potts로 넘어갈 가능성이 있다.
   - 이 가설이 맞으면 "μ_c"는 날카로운 점이 아니라 크로스오버 척도다.
3. μ_c 값. 지금 말할 수 있는 것은 2 ≲ μ_c ≲ 3뿐이다(μ=2는 장거리 질서, μ=3은 멱법칙).
4. α/ν(비열 배경 때문에 미결정), μ>0 정밀 지수, μ→∞에서 T_c가 dice 격자 이징값 1.2027로 가는지.
5. 결함 활성화 기울기의 prefactor 보정. TAFM ξ는 시드 1개뿐이다.

## 6. 논문 구성 제안

- **본문의 뼈대 (확정 결과만):** 모형 → 정확한 대응(결과 1–3) → 계산 절차 검증(5) → D1=D2에서 전이 없음과 ξ∝exp(2/T)(6–7) → 결함 비용 2배 메커니즘(4) → D2>0의 Potts 전이와 TAFM 대비 T_c 상승(8) → μ 조절(9–10) → 비열 분해(13).
- **큰 μ는 짧은 절 하나로:** 결과 11–12를 관측으로 보고하고, 5절의 가설은 전망으로만 쓴다.
- **선행연구와의 차별점:**
  - 가장 가까운 연구는 Shokef–Souslov–Lubensky, PNAS 108, 11804 (2011)다. 탄성 격자의 연속 변위로 퇴화를 푸는 order-by-disorder다.
  - 우리 모형은 이산적인 결합 선택으로 결함 에너지 자체를 두 배로 만든다.
  - "좌절 결합 = dimer" 대응 자체는 알려진 것(TAFM ↔ 벌집 dimer)이다. 새로운 것은 그 dimer를 독립적인 annealed 변수로 만든 결과다.
- **남은 문헌 확인:**
  - Johnston 1994 (Phys. Lett. B, hep-th/9406138)를 아직 확인하지 않았다.
  - Yin–Gross–Chakraborty 저자 목록을 확인해야 한다.
  - arXiv:2306.09046의 저널 정보를 채워야 한다.
  - 같은 모형(matching 제약이 있는 annealed-bond TAFM)의 선행연구가 있는지 한 번 더 검색해야 한다.

## 7. 그림 목록

| 그림 | 내용 | 파일 / 데이터 |
|---|---|---|
| 1 | MC 스냅샷 4개: 바닥상태, D1=D2 무작위 타일링, √3 결정, 무질서 | figures/fig_snapshots.png |
| 2 | Binder U4(T), L=32/48/64, μ=0 | data/fss_d21_summary.txt (원 시계열은 저장소 results/long) |
| 3 | 결함 비용 2배 도식 | figures/fig_mechanism.png |
| 3b | ξ(T): annealed vs TAFM (ln ξ vs 1/T) | data/d1d2_summary.txt, data/tafm_xi_L32_s1.txt |
| 3c | 결함·좌절 밀도 vs 1/T | data/tafm_defects_L24_s1.txt |
| 4 | T_c 근처 비열, annealed vs TAFM, C_μ/C_n | paper_draft.md 표, data/fss_d21_summary.txt, data/tafm_fss_summary.txt |
| 5 | T_c(μ) | data/key_numbers.csv |
| 6 | 유효 η vs T/T_c, μ별 (큰 μ 절) | data/key_numbers.csv, data/seeds_summary.txt |

- 초안 문서의 인터랙티브 차트 6개는 Claude Docs 안에만 있다. paper_draft.md에는 `[차트/위젯: …]` 자리표시로만 남아 있다.
- 출판용 그림은 위 데이터로 새로 그려야 한다.

## 8. 폴더 구성

```
handoff/
  HANDOFF.md                 이 문서
  paper_draft.md             논문 초안 (Claude Docs rev 18을 마크다운으로 변환)
  figures/                   fig_snapshots.png, fig_mechanism.png
  data/
    key_numbers.csv          T_c(μ), η 맞춤 등 핵심 수치 한 표
    fss_d21_summary.txt      μ=0 장시간 런 FSS (T_c, 지수, 진단)
    fss_mu_summary.txt       μ>0 FSS (허브 질서변수)
    fss_mu_spin_summary.txt  μ>0 FSS (스핀 질서변수)
    seeds_summary.txt        μ=3,5 시드 평균 ψ_s² vs L, η 맞춤
    clock_mu10.txt           μ=10 ⟨cos3θ⟩, ⟨cos6θ⟩
    tafm_fss_summary.txt     고정 결합 TAFM (h/J=1, 3) FSS
    tafm_xi_L32_s1.txt       TAFM 영자기장 ξ(T)
    tafm_defects_L24_s1.txt  결함·좌절 밀도, annealed vs TAFM
    d1d2_summary.txt         D1=D2 PT 요약
    mechanism_summary.txt    메커니즘·μ 분석 작업 기록 (시간순, 수정 이력 포함)
  notes/
    PART9.md                 모형 정의와 초기 결과 메모
    PART9_literature.md      문헌 위치 정리와 인용 목록
  code/                      핵심 스크립트 (재현용, 저장소 없이는 실행 불가)
```

## 9. 재현과 추가 계산 (계산 세션에서만)

- 원 시계열(.npz, 약 340 MB)은 저장소 `results/long`, `results/tafm`, `results/d1d2`에 있다. 채팅이나 프로젝트로는 넘기지 않았다.
- 주요 명령:
  - `python3 annealed_long.py <L> <T> <n_total> <thin> <chunk> <seed> [D2] [mu]`
  - `ORDER=spin MU=<μ> python3 annealed_long_analysis.py 1 32 48 64`
  - `PYTHONPATH=. python3 results/long/seeds_analysis.py`
  - `MU=10 python3 results/long/clock_analysis.py`
- 실무 주의:
  - numba 캐시 경쟁이 생기므로 작은 작업 하나로 먼저 컴파일한다.
  - `pkill -f`로 자기 셸을 죽이지 않게 PID로 종료한다.
  - 컨테이너가 자주 재시작되니 결과는 바로 커밋·푸시한다.
- 비용 감각: L=96 한 점(5×10⁵ 스윕)에 CPU 약 4–5시간이 든다. 자기상관시간은 최대 2.5만 스윕이다.

## 10. 남은 할 일 (우선순위)

1. **원고 정리 (계산 불필요):**
   - 초록 과장 수정(4절).
   - μ_c·중간상 서술을 5절 수준으로 낮추기.
   - 출판용 그림 5–6개 만들기.
   - 부록 작성(알고리즘, 평형 진단).
2. **문헌 확인:** 6절 목록.
3. **공유 문서 §11 갱신 (선택):** https://claude.ai/artifact/9jxHSSeFMssZzDUT7Ph4eS 의 §11 Part IX는 이번 결과로 아직 갱신하지 않았다. 갱신 전에 먼저 읽을 것.
4. **추가 계산 (계산 세션으로 돌아갈 경우, 선택):**
   - (a) μ=0의 T_c에서 스핀 ψ_s의 η가 4/15인지 확인(L=32–96). 5절 가설을 가장 직접 판별한다.
   - (b) μ=10에 KT형 맞춤(기존 데이터로 가능).
   - (c) 결함 활성화의 prefactor 맞춤, TAFM ξ 시드 추가.

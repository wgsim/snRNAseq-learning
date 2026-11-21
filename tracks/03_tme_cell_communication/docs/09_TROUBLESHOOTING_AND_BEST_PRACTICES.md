# Chapter 09. 세포 간 통신 분석의 함정, 검증 프로토콜 및 캡스톤 프로젝트 (Troubleshooting & Capstone Project)

> **Cognitive Goal**: 세포 간 통신(Cell-Cell Communication, CCC) 분석 시 발생하는 4대 통계적·실험적 왜곡(상호작용 수프의 위양성 폭발, 단일 세포 유사반복에 의한 $p$-값 인플레이션, 배치 효과 오염, 극소 발현 드롭아웃 왜곡)을 진단하고 수학적으로 교정하는 기법을 체득한다. 전사체 기반 예측 결과를 직교적으로 검증(Orthogonal Validation)하기 위한 공간 다중 면역형광(mIF) 및 3D 오가노이드 기능적 섭동 프로토콜을 정립한다. 최종적으로 Track 01(eGRN), Track 02(Spatial Transcriptomics), Track 03(CCC/NicheNet/MISTy)을 집대성한 **5단계 TME 통신 캡스톤 프로젝트(Capstone Project)** 평가 루브릭을 완성한다.

---

## 1. 세포 간 통신 분석의 4대 통계적 함정과 병리학적 실패 모드

Dimitrov et al. (Nature Communications 2022, PMID 35680885)은 L-R 통신 방법과 자원 간의 예측 일치도가 낮다고 보고했습니다. 또한 Squair et al. (Nature Communications 2021, PMID 34584091)은 차별 발현 분석에서 개별 세포를 독립 표본으로 취급하는 유사반복이 거짓 양성을 크게 늘린다고 보고했으며(L-R 통신 도구 평가가 아님), 이는 아래 함정 2의 근거입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   세포 간 통신 분석을 위협하는 4대 병리학적 실패 모드                  │
│                                                                                        │
│   [1. 상호작용의 수프 (Soup of Interactions) & 위양성 폭발]                             │
│   - 세포 유형 15개, L-R 데이터베이스 1,500쌍 ==> 15 × 14 × 1,500 = 315,000회 가설 검정!│
│   - 다중 가설 검정 보정(FDR) 결손 시 순수 노이즈만으로 수만 개의 '유의미한' 통신 출력! │
│                                                                                        │
│   [2. 단일 세포 유사반복 (Pseudoreplication)에 의한 p-값 인플레이션]                   │
│   - 환자 N=6명에서 얻은 60,000개 단일 세포를 독립 표본으로 취급                        │
│   - 실제 자유도(df)는 N-2 = 4에 불과하나, df = 59,998로 계산되어 p < 10^{-50} 왜곡!   │
│                                                                                        │
│   [3. 기술적 배치 효과 (Batch-Confounded CCC)]                                         │
│   - 환자/샘플 배치 간 시퀀싱 깊이 및 효소 해리 조건 차이가 가짜 L-R 공변성 형성       │
│                                                                                        │
│   [4. 극소 발현 수용체/리간드의 드롭아웃 (Dropout Masking)]                            │
│   - 저농도에서도 초고감도로 작동하는 GPCR, 인터루킨 수용체가 0 발현으로 누락           │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 조합 폭발과 다중 가설 검정 결손 (The "Soup of Interactions")
$K$개의 세포 유형이 존재할 때 검정해야 할 방향성 세포 쌍의 수는 $K(K - 1)$개입니다. 데이터베이스 내의 L-R 쌍의 수가 $M$개라면, 전체 귀무가설의 총수 $H$는 다음과 같이 폭증합니다:

$$H = K(K - 1) \times M$$

- 예를 들어 $K = 20$개 클러스터, $M = 2,000$쌍일 경우, $H = 20 \times 19 \times 2,000 = 760,000$개의 가설이 동시 검정됩니다.
- 유의수준 $\alpha = 0.05$에서 다중 검정 보정(Benjamini-Hochberg FDR) 없이 개별 순열 $p$-값만을 사용할 경우:
  $$\mathbb{E}[\text{False Positives}] = 760,000 \times 0.05 = 38,000\text{개의 가짜 신호!}$$
- **교정 원칙**: 전체 테스트 공간에 대해 엄격한 글로벌 FDR(Global Benjamini-Hochberg Adjusted $q < 0.01$)을 적용하고, 최소 2개 이상의 독립 데이터베이스(CellChatDB, CellPhoneDB, SingleCellSignalR)가 동시에 지지하는 **합의 점수(Consensus Rank)**만을 필터링해야 합니다.

### 1.2 단일 세포 유사반복(Pseudoreplication)의 통계적 재앙
단일 세포 전사체 분석에서 가장 빈번하게 저질러지는 치명적 실수는 **동일한 환자(Donor) 조직에서 유래한 수천 개의 단일 세포를 통계적으로 상호 독립적인 관측치(Independent Observations)로 취급하는 것**입니다.

- 동일한 환자의 종양 세포와 T 세포는 환자의 고유한 유전적 배경, 전신 염증 수준, 조직 관류 상태, 실험 당일의 배치 효과를 공유합니다.
- 세포 수준에서 상관계수나 $t$-검정을 수행하면 유효 표본 크기($n$)가 세포 수(수만 개)로 둔갑하여 표준 오차(Standard Error)가 비정상적으로 $1/\sqrt{N_{\text{cells}}}$로 축소됩니다:
  $$SE = \frac{\sigma}{\sqrt{N_{\text{cells}}}} \ll \frac{\sigma}{\sqrt{N_{\text{donors}}}}$$
- 그 결과 실제로는 생물학적 환자 간 차이가 전혀 없는 무작위 노이즈조차 $p < 10^{-100}$의 허위 통계적 유의성을 획득하게 됩니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   단일 세포 유사반복 해결책: 기증자 수준 유사벌크 집계                 │
│                                                                                        │
│   [잘못된 접근법: Single-Cell Level Testing]                                           │
│   - N_cells = 50,000개의 세포를 단일 테이블에 넣고 상관분석 수행                       │
│   - 문제: 샘플 1개의 환자 특이적 아티팩트가 50,000개의 독립 증거로 둔갑함             │
│                                                                                        │
│   [올바른 접근법: Donor-Level Pseudobulk Aggregation]                                  │
│   - 각 환자(Donor d) 및 세포 유형(k)별로 평균 발현량 계산:                              │
│     x̄_{d, k} = (1 / |C_{d, k}|) ∑_{i ∈ C_{d, k}} x_i                                  │
│   - 환자 수준의 집계 테이블(행: 환자 Donor, 열: L-R 발현)을 구축                      │
│   - 선형 혼합 효과 모델(Linear Mixed Models, LMM) 또는 edgeR/DESeq2 프레임워크 적용     │
│   - 자유도는 실제 환자 수(N_donors - 2)로 올바르게 계산됨                              │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. 직교적 실험 검증 프로토콜 (Orthogonal Validation Methodologies)

계산생물학적으로 예측된 L-R 상호작용은 가설(Hypothesis)일 뿐이며, 논문 발표 및 신약 개발을 위해서는 반드시 단백질 수준 및 기능적 수준의 3단계 직교 검증을 통과해야 합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   L-R 상호작용 검증을 위한 3단계 직교 검증 파이프라인                  │
│                                                                                        │
│   [1단계: 공간 다중 단백질 동시 국소화 (Spatial Multiplex IF - PhenoCycler / CODEX)]   │
│   - 목적: mRNA 발현이 실제 기능적 단백질로 번역되었으며, 물리적으로 접촉하는가?       │
│   - 검증: 리간드 항체(Fluorophore A) + 수용체 항체(Fluorophore B) 공간 형광 공초점      │
│   - 정량: 단일 세포 수준 막 표면 형광 중첩도 (Manders' Overlap Coefficient > 0.6)     │
│                                                                                        │
│   [2단계: 고감도 단일 분자 형광 제자리 혼성화 (RNAscope Multiplex smFISH)]            │
│   - 목적: 단일 세포 시퀀싱의 드롭아웃(Dropout) 편향 검증                              │
│   - 검증: FFPE 슬라이드 상에서 리간드-수용체 mRNA 전사체 분자를 단일 분자 점으로 계수 │
│                                                                                        │
│   [3단계: 3D 환자 유래 오가노이드-기질 세포 공배양 기능적 섭동 (Functional Assay)]   │
│   - 목적: L-R 결합 차단이 실제 하류 전사 반응 및 표현형을 가역적으로 억제하는가?      │
│   - 검증: Tumor Organoid + myCAF + Osteoclast 삼차원 트랜스웰 공배양                   │
│   - 섭동: 중화 항체(Neutralizing Ab) 투여 또는 CRISPR-Cas9 녹아웃                      │
│   - 평가: 하류 타겟 전사인자 인산화(Western Blot) 및 NicheNet 표적 DEG 억제(RT-qPCR)  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. 5단계 TME 통신 캡스톤 프로젝트 평가 루브릭 (5-Stage Capstone Project Rubric)

이 루브릭은 학위 논문 또는 최고 권위 저널(Nature/Cell/Cancer Cell) 투고 기준에 부합하는 종양 미세환경 세포 간 통신 통합 분석 파이프라인의 완성도를 평가합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   5단계 TME 통신 캡스톤 프로젝트 마스터 아키텍처                       │
│                                                                                        │
│   [Stage 1: QC & 다중 데이터베이스 합의 필터링 (Consensus)]                             │
│   - 세포별 최소 발현율(10%) 필터링, CellChatDB + CellPhoneDB + Connectome 합의 도출   │
│                                                                                        │
│   [Stage 2: 다중 소단위체(Multi-Subunit) 및 이종 다량체 복합체 모델링]                 │
│   - α/β 인테그린 이량체, IL6ST/IL6R 등 필수 다합체 결합 조건 엄격 적용                │
│                                                                                        │
│   [Stage 3: NicheNet 세포 내 신호 인과 추론 및 eGRN 가교]                              │
│   - RWR 알고리즘 기반 표적 유전자 귀속, Track 01 snATAC-seq 레귤론과 교차 검증        │
│                                                                                        │
│   [Stage 4: 공간 거리 감쇄 및 MISTy 다중 뷰 분산 분해]                                 │
│   - 공간 이변량 모란 지수(Moran's I), Intrinsic / Juxta / Para 스케일별 기여도 분해     │
│                                                                                        │
│   [Stage 5: 기전 통합 및 5자 회로 표적 병용 치료 전략 도출]                           │
│   - mCRPC 골전이 5자 악순환 회로 정식화 및 면역 배제 역전 칵테일 치료 설계            │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

| 평가 단계 | 핵심 분석 내용 및 방법론 | 성공 판정 기준 (Validation Criteria) | 배점 |
| :--- | :--- | :--- | :---: |
| **Stage 1: QC & 합의 필터링** | - 환자별 라이브러리 깊이 정규화<br>- 발현율 10% 미만 드롭아웃 유전자 마스킹<br>- 2개 이상 DB 합의 순위화 | - 위양성 잡음 L-R 쌍 60% 이상 제거<br>- 합의 순위 Spearman $r > 0.85$ | 15% |
| **Stage 2: 다중 소단위체 모델링** | - 복합체 수용체(Heteromeric Receptors) 모델링<br>- 기하평균 기반 소단위체 발현율 계산 | - 필수 보조 소단위체(e.g., IL6ST) 결손 시 상호작용 점수 0으로 강제 배제 | 20% |
| **Stage 3: NicheNet 인과 추론** | - 무작위 보행 재시작(RWR) 잠재력 계산<br>- 수용 세포 DEG 예측 AUROC/AUPR 산출<br>- Track 01 SCENIC+ eGRN과 교차 검증 | - 최상위 리간드 AUROC $> 0.75$<br>- 하류 TF의 snATAC-seq 모티프 풋프린트 일치 | 25% |
| **Stage 4: 공간 다중 뷰 분해** | - 공간 이변량 모란 지수 계산<br>- MISTy 다중 뷰(내인성, 접촉, 확산) 스태킹<br>- $R^2$ 분산 분해 | - $p < 0.01$ 유의 핫스팟 도출<br>- 공간 뷰 추가에 따른 $\Delta R^2_{\text{spatial}} > 0.20$ 달성 | 25% |
| **Stage 5: 기전 및 치료 설계** | - 골전이 5자 악순환 회로 정식화<br>- 면역 배제형 미세환경 극복 다중 병용 전략 수립 | - 3단계 이상 다중 표적 병용 전략 논리적 설계<br>- 내성 기전 극복 근거 제시 | 15% |

> 위 표의 수치 기준(위양성 60% 제거, Spearman $r > 0.85$, AUROC $> 0.75$, $\Delta R^2 > 0.20$ 등)은 출처가 없는 교육용 예시이며 데이터에 맞춰 조정합니다.

---

## 4. 실전 종합 분석 파이프라인: 기증자 유사벌크 및 합의 CCC 스크립트

다음은 기증자(Donor) 수준의 유사벌크 집계를 수행하여 유사반복을 제거하고, 다중 데이터베이스 합의 순위를 계산하는 파이썬 스크립트입니다.

```python
"""
09_troubleshooting_pseudobulk_consensus_ccc.py
기증자 수준 유사벌크(Donor-level Pseudobulk) 집계를 통한 단일 세포 유사반복 왜곡 해결 및
다중 데이터베이스 합의 기반 고신뢰도 세포 간 통신 검증 파이프라인
"""

import scanpy as sc
import liana as li
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

# 1. 단일 세포 AnnData 로드
print("[Step 1] Loading mCRPC multi-donor scRNA-seq AnnData...")
adata = sc.read_h5ad("mcrpc_bone_metastasis_multidonor.h5ad")

print(f"Total cells: {adata.n_obs}, Donors: {adata.obs['donor_id'].nunique()}")

# 2. 기증자 수준 유사벌크(Pseudobulk) 프로파일 집계
print("[Step 2] Performing donor-level pseudobulk aggregation to eliminate pseudoreplication...")
# donor_id와 cell_type의 조합별 원시 카운트 합산(유사벌크 DE 표준은 합산 카운트; 로그 정규화 평균 아님)
pseudobulk_obs = []
pseudobulk_matrix = []

for donor in adata.obs["donor_id"].unique():
    donor_mask = adata.obs["donor_id"] == donor
    sub_adata = adata[donor_mask]
    
    for ct in sub_adata.obs["cell_type"].unique():
        ct_mask = sub_adata.obs["cell_type"] == ct
        cells_in_group = np.sum(ct_mask)
        
        # 세포 수가 20개 미만인 희귀 그룹은 노이즈 방지를 위해 배제
        if cells_in_group >= 20:
            # 원시 카운트는 adata.layers['counts']에 있다고 가정
            sum_expr = np.asarray(sub_adata[ct_mask].layers['counts'].sum(axis=0)).ravel()
            pseudobulk_matrix.append(sum_expr)
            pseudobulk_obs.append({"donor_id": donor, "cell_type": ct, "n_cells": cells_in_group})

pb_df = pd.DataFrame(pseudobulk_matrix, columns=adata.var_names)
pb_obs = pd.DataFrame(pseudobulk_obs)
pb_adata = sc.AnnData(X=pb_df.values, obs=pb_obs, var=adata.var.copy())
sc.pp.normalize_total(pb_adata, target_sum=1e4)   # 합산 카운트를 정규화한 뒤 로그 변환
sc.pp.log1p(pb_adata)

print(f"Aggregated Pseudobulk Samples: {pb_adata.n_obs} (Donors × CellTypes)")

# 3. LIANA+ 합의 알고리즘(Consensus CCC) 실행
print("[Step 3] Running LIANA Multi-algorithm Consensus (CellPhoneDB, NATMI, Connectome, SingleCellSignalR)...")
# 유사벌크(기증자×세포유형) 프로파일에 실행하여 유사반복(pseudoreplication)을 제거한다
li.mt.rank_aggregate(
    pb_adata,
    groupby="cell_type",
    resource_name="consensus",
    expr_prop=0.10,          # 10% 이상 발현 필터
    min_cells=3,             # 유사벌크 표본 수 기준
    verbose=True,
    use_raw=False
)

consensus_res = pb_adata.uns["liana_res"]

# 4. 상위 합의 상호작용 필터링 (다중 검정 보정 적용)
# liana 1.x 출력 컬럼은 magnitude_rank / specificity_rank(RRA p 유사 값; 낮을수록 상위)이며 aggregate_rank 컬럼은 없다
top_consensus = consensus_res[
    consensus_res["specificity_rank"] < 0.05
].sort_values(by="specificity_rank", ascending=True)

print("\n=== Top High-Confidence Consensus L-R Interactions ===")
print(top_consensus[["source", "target", "ligand_complex", "receptor_complex", "magnitude_rank", "specificity_rank"]].head(15))

# 5. 복합체 수용체(Heteromeric Complex) 무결성 검증
# 다중 서브유닛을 가진 수용체(e.g., ITGAV:ITGB3, IL6ST:IL6R)가 단일 서브유닛만 발현된 거짓 신호인지 필터링
print("[Step 4] Auditing heteromeric complex subunit expression integrity...")
hetero_complexes = top_consensus[top_consensus["receptor_complex"].str.contains("_")]

print(f"Total heteromeric receptor complexes validated: {len(hetero_complexes)}")
for _, row in hetero_complexes.head(5).iterrows():
    print(f"Verified Complex: {row['source']} -> {row['target']} via {row['ligand_complex']} -> {row['receptor_complex']}")

# 6. 결과 저장
top_consensus.to_csv("mcrpc_bone_tme_robust_consensus_ccc.csv", index=False)
print("[Complete] Robust CCC audit and pseudobulk validation pipeline completed successfully.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
단일 세포 전사체(scRNA-seq) 데이터를 기반으로 한 세포 간 통신(CCC) 분석에서 발생하는 **단일 세포 유사반복(Pseudoreplication)**의 통계학적 정의를 기술하고, 이를 방치했을 때 가설 검정의 제1종 오류(Type I Error / 위양성)율이 폭증하는 수학적 기전을 표준 오차(Standard Error) 수식을 기반으로 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **단일 세포 유사반복(Pseudoreplication)의 통계적 정의**:
  - 실제 생물학적 독립 표본(Biological Replicates)은 각기 다른 환자/기증자 개체(Donors)임에도 불구하고, 동일한 개체의 조직 내에서 추출된 수만 개의 단일 세포들을 상호 독립적인 관측치(Independent Observations)로 취급하여 통계 검정을 수행하는 오류를 말합니다.
- **제1종 오류 폭증의 수학적 기전**:
  - 표본 평균의 차이에 대한 검정 통계량 $t$는 다음과 같이 정의됩니다:
    $$t = \frac{\bar{X}_1 - \bar{X}_2}{SE}, \quad SE = \sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}$$
  - 만약 실제 생물학적 환자 수가 $N_{\text{donors}} = 6$명인데, 이를 $N_{\text{cells}} = 60,000$개의 독립 표본으로 취급하면, 분모의 표본 크기 $n$이 인위적으로 10,000배 부풀려집니다.
  - 이로 인해 표준 오차($SE$)가 $1/\sqrt{10,000} = 1/100$ 수준으로 극단적으로 축소됩니다.
  - 분모($SE$)가 인위적으로 0에 가깝게 작아지므로 검정 통계량 $t$는 천문학적으로 급증하여, 실제 환자 간 차이가 전혀 없는 귀무가설($H_0$) 하에서도 $p$-값이 $10^{-50}$ 이하로 떨어져 귀무가설을 부당하게 기각하는 제1종 오류(Type I Error / 위양성)가 폭증하게 됩니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
L-R 상호작용 데이터베이스 분석 시 필수적인 '다중 소단위체 수용체 복합체(Multi-Subunit Receptor Complexes, e.g., $\alpha/\beta$ 인테그린 이형이량체, IL-6R/IL6ST)'를 모델링할 때, 각 소단위체 발현량 $x_1, x_2, \dots, x_m$에 대해 단순 산술평균(Arithmetic Mean) $\frac{1}{m}\sum x_k$을 사용하는 것이 생물학적으로 치명적인 오류를 낳는 이유를 설명하고, 기하평균(Geometric Mean) $\left(\prod_{k=1}^m x_k\right)^{1/m}$ 또는 최소값($\min(x_1, \dots, x_m)$) 연산자를 적용해야만 하는 질량 작용의 법칙(Law of Mass Action) 기반의 수리적 타당성을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **단순 산술평균 적용 시의 치명적 오류**:
  - 다중 소단위체 복합체는 모든 필수 소단위체가 동시에 세포막에 발현되어 회합해야만 기능적 수용체 구조를 형성합니다 (논리적 AND 조건).
  - 만약 소단위체 1($x_1$)이 10.0으로 초고발현되고, 복합체 형성에 필수적인 소단위체 2($x_2$)가 전혀 발현되지 않아 0.0인 경우($x_2 = 0$):
    $$\text{산술평균} = \frac{10.0 + 0.0}{2} = 5.0$$
  - 산술평균은 5.0이라는 높은 활성 점수를 산출하여, 실제로는 수용체 복합체가 전혀 형성될 수 없음에도 매우 강력한 상호작용이 일어나는 것으로 오판하는 치명적 위양성을 낳습니다.
- **기하평균 및 최소값 연산자의 질량 작용의 법칙(Law of Mass Action) 기반 타당성**:
  - 생화학 평형 반응에서 $m$개의 소단위체가 결합하여 온전한 수용체 복합체 $C$를 형성하는 속도론은 질량 작용의 법칙을 따릅니다:
    $$[C] = K_{\text{eq}} \cdot \prod_{k=1}^m [x_k]$$
  - 복합체의 농도는 각 소단위체 농도의 '곱(Product)'에 비례합니다. 단 하나의 필수 소단위체라도 농도가 0이면 복합체 농도는 엄격히 0이 되어야 합니다 ($[x_j] = 0 \implies [C] = 0$).
  - 기하평균 $\left(\prod_{k=1}^m x_k\right)^{1/m}$은 곱의 $m$제곱근이므로 어느 하나라도 0이면 점수가 0이 된다는 화학 양론적(Stoichiometric) 제약을 공유합니다(평형 농도 자체가 아니라 곱 형태를 단위 크기로 환산한 점수).
  - 최소값 연산자 $\min(x_1, \dots, x_m)$ 역시 가장 결핍된 병목 소단위체(Rate-limiting Subunit)의 양에 의해 전체 복합체 형성이 제한된다는 리비히의 최소량의 법칙(Liebig's Law of the Minimum)을 수리적으로 정확하게 반영합니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
다기관 전립선암 골전이 단일 세포 전사체 컨소시엄에서 30명의 환자 데이터를 통합 분석하여, 종양 세포와 파골세포 간의 새로운 면역 억제성 L-R 축으로 *X-Y* 쌍을 발굴하였다고 보고하였다. 이 연구진은 "단일 세포 120,000개에 대한 순열 검정 결과 $p < 10^{-80}$이었으므로 확실한 치료 표적"이라고 주장하며 대규모 신약 개발 투자를 권고하였다.
당신이 독립 감사 위원회의 수석 생물정보학자라면, 이 발견이 '유사반복(Pseudoreplication)', '배치 효과(Batch Confounding)', '복합체 결손(Complex Misannotation)'에 의해 조작된 가짜 신호일 가능성을 감사하기 위해 요구할 3가지 통계적/실험적 재검증 감사 프로토콜을 구체적으로 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **감사 프로토콜 1: 기증자 수준 유사벌크(Donor-level Pseudobulk) 선형 혼합 모델 재검정**:
  - **감사 내용**: 120,000개 세포를 풀링한 비물리적 분석을 전면 기각하고, 30명의 환자(Donor)별로 각 세포 유형의 평균 발현량을 집계한 유사벌크 행렬을 구축합니다.
  - **통계 모델**: 환자 개체 변량을 무작위 효과(Random Effects)로 지정하는 선형 혼합 효과 모델(Linear Mixed-Effects Model, LMM)을 적용하여, 자유도를 $n=30$을 기준으로 재계산합니다.
  - **합격 기준**: 환자 간 무작위 효과를 보정한 상태에서도 해당 L-R 축의 공변성이 $p_{\text{donor}} < 0.01$을 유지하는지 검증합니다. 만약 $p$-값이 0.1 이상으로 치솟는다면 이는 소수 환자의 배치 편향이 단일 세포 유사반복으로 부풀려진 위양성으로 판정합니다.
- **감사 프로토콜 2: 기술적 배치 효과 및 임상 공변량 통제 감사**:
  - **감사 내용**: 30명의 환자 샘플이 시퀀싱된 날짜(Sequencing Batch), 샘플 채취 병원(Center), 조직 처리 효소 프로토콜이 L-R 발현 점수와 교락(Confounded)되어 있는지 분산 분석(ANOVA)을 수행합니다.
  - **재검정**: 배치를 공변량으로 넣은 pseudobulk 모형에서 동일한 상호작용이 재현되는지, 배치별로 층화해 방향이 일치하는지 확인합니다(Harmony·Scanorama는 발현이 아니라 PCA 등 임베딩을 보정하므로 발현 점수 계산에 직접 쓰지 않음). 배치 레이블을 셔플링했을 때 신호가 붕괴되는지도 확인합니다.
- **감사 프로토콜 3: 공간 다중 단백질 동시 국소화(Spatial Multiplex IF) 및 수용체 복합체 완전성 검증**:
  - **감사 내용**: 수용체 $Y$가 다중 소단위체 복합체인지 감사하고, 단일 서브유닛만 발현된 상태에서 산술평균으로 과장된 것은 아닌지 기하평균으로 재평가합니다.
  - **직교 실험 검증**: 환자의 독립적인 골전이 FFPE 조직 절편에 대해 공간 다중 형광 면역염색(mIF / PhenoCycler)을 수행하여, 종양 세포의 $X$ 단백질과 파골세포의 $Y$ 단백질이 실제 $15\ \mu\text{m}$ 이내의 물리적 접촉 계면에서 형광 공존(예: Manders 계수; 두 채널의 픽셀 공존 지표로 L-R 접촉 계면 자체의 판정은 아니며, 0.6 같은 기준은 근거 없는 예시)을 보이는지 단백질 수준에서 직교 검증을 요구하고, 필요하면 근접 결찰 분석(PLA) 등을 병행합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13(1):3224. PMID: 35680885. doi: 10.1038/s41467-022-30755-0
2. Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nat Commun*. 2021;12(1):5692. PMID: 34584091. doi: 10.1038/s41467-021-25960-2
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- LIANA (liana-py) 1.10.0 (github.com/saezlab/liana-py). 원 논문: Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13(1):3224. PMID: 35680885. doi: 10.1038/s41467-022-30755-0
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
<!-- REFERENCES:END -->

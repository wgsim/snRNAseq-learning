# 09. 트러블슈팅, 실험 설계 및 통계 검증 가이드 (Troubleshooting, Experimental Design & Statistical Benchmarking)

---

## 학습 목표 (Cognitive Goals)

1. **단일세포 실험 설계의 황금률(Golden Rule): 공여자 수($N$) 대 세포 심도($n$)의 통계적 검정력(Power) 비교 정식화**: 생물학적 분산($\sigma_D^2$)과 기술적/세포 내 분산($\sigma_c^2$)의 분산 성분 분해를 통해, 3명 $\times$ 10,000개 세포 설계 대비 30명 $\times$ 1,000개 세포 설계가 차등 발현(DEG) 검정력을 획기적으로 향상시키는 수리적 기전을 증명한다.
2. **배치 효과 교락(Batch Confounding)의 진단 및 잠재 변수 잔차 제거 기전 수립**: 시퀀싱 배치와 실험군이 뒤섞인 불완전 블록 설계에서 단순 배치 제거의 생물학적 신호 파괴 위험을 진단하고, 대리 변수 분석(SVA) 및 잔차 변동 제거(RUVseq)를 GLM 모델 행렬에 투입하는 회귀적 해결책을 도출한다.
3. **독립 가설 가중치(Independent Hypothesis Weighting, IHW)의 적응적 다중 검정 보정 유도**: 가설별 통계적 검정력의 이질성을 반영하여, 유전자 평균 발현량 공변량 기반의 볼록 최적화(Convex Optimization)를 통해 가중치 $w_g$를 배분함으로써 진정한 생물학적 표적의 발굴률을 극대화한다.
4. **5단계 단일세포 카운트 통계 및 풍부도 분석 종합 캡스톤 프로젝트 루브릭(Capstone Rubric) 완성**: 실험 설계부터 세포 수준 GLMM(NEBULA), 궤적 동역학(GAMM), 국소 풍부도(Milo), 전역 조성(scCODA)까지 전 과정을 아우르는 임상 중개 전립선암 데이터 분석의 재현성 평가 기준을 확립한다.

---

## 1. 단일세포 통계학의 황금률: 세포 수($n$) vs 공여자 수($N$) 파워 분석

단일세포 시퀀싱 프로젝트를 기획할 때 가장 자주 마주치는 자원 배분 딜레마는 다음과 같다: **"동일한 시퀀싱 예산(예: 총 30,000개 세포)이 주어졌을 때, 군당 3명(총 6명)의 환자로부터 각 5,000개 세포를 뽑을 것인가, 아니면 군당 30명(총 60명)의 환자로부터 각 500개 세포를 뽑을 것인가?"**

통계학적 답은 대체로 명확하다: 총 세포 수가 고정되고 공여자 확보 비용을 무시하는 가정에서는 **공여자 수($N$)를 늘리는 것이 세포 수를 늘리는 것보다 훨씬 효율적이다.**

```
[The Allocation Dilemma: Sample Size vs Cell Depth Power Analysis]

Design A: "Deep Cell Profiling"
  군당 3 Donors (총 6) x 5,000 Cells = 30,000 Cells
  - Degrees of Freedom (Biological): df = 6 - 2 = 4!
  - Standard Error: SE(beta) = sqrt[ 2 sigma_D^2 / 3 + 2 sigma_c^2 / (3*5,000) ] ≈ 0.5775  (sigma_D^2 = 0.5, sigma_c^2 = 1.0)
  - Power (beta* = 0.3, alpha = 0.05, 양측): ≈ 7~8%  (Underpowered, High Type II Error)

Design B: "Broad Donor Profiling"
  군당 30 Donors (총 60) x 500 Cells = 30,000 Cells
  - Degrees of Freedom (Biological): df = 60 - 2 = 58!
  - Standard Error: SE(beta) = sqrt[ 2 sigma_D^2 / 30 + 2 sigma_c^2 / (30*500) ] ≈ 0.1829
  - Power (beta* = 0.3, alpha = 0.05, 양측): ≈ 36~38%  (훨씬 높지만 여전히 충분하지는 않음)
```

### 1.1 분산 성분 분해(Variance Components Decomposition)의 수학적 유도

처리군(Treated)과 대조군(Control) 간의 평균 차이 $\hat{\beta} = \bar{Y}_{\text{Treat}} - \bar{Y}_{\text{Ctrl}}$의 추정 분산은 두 가지 독립적인 노이즈 원천으로 분해된다:
1. **공여자 간 생물학적 분산($\sigma_D^2$)**: 인간 개체 간 유전적, 면역학적, 대사적 차이.
2. **세포 간 기술적/생물학적 분산($\sigma_c^2$)**: UMI 포획 효율, 전사 버스팅(Transcriptional Bursting), 세포 내 미세 이질성.

군당 공여자 수가 $N$명이고 공여자당 포획된 세포 수가 $n$개일 때, 추정치의 표준 오차 $\text{SE}(\hat{\beta})$는 다음과 같다:

$$\text{SE}(\hat{\beta}) = \sqrt{ \frac{2 \sigma_D^2}{N} + \frac{2 \sigma_c^2}{N \cdot n} }$$

실제 포유류 조직 및 인간 종양 전사체에서 $\sigma_D^2$와 $\sigma_c^2$의 크기를 비교하면, 생물학적 환자 간 변이가 세포 간 변이보다 압도적으로 우세하다($\sigma_D^2 \gg \frac{\sigma_c^2}{n}$ for $n \ge 500$).

- **세포 수 $n \to \infty$ 극한**: 공여자당 세포 수를 아무리 무한대로 늘려도, 두 번째 항 $\frac{2\sigma_c^2}{N \cdot n}$만 0으로 수렴할 뿐 첫 번째 항 $\frac{2\sigma_D^2}{N}$는 전혀 줄어들지 않는다:
  $$\lim_{n \to \infty} \text{SE}(\hat{\beta}) = \sqrt{\frac{2\sigma_D^2}{N}}$$
- 따라서 군당 환자 수 $N$이 작으면(e.g., $N=3$), 세포를 10만 개 시퀀싱하더라도 환자 간 고유 변이의 장벽에 가로막혀 통계적 검정력(Power)이 결코 증가하지 않는다.

---

## 2. 배치 효과 교락(Batch Confounding)의 진단 및 잔차 모델링

단일세포 데이터에서 배치 효과(Batch Effect)는 서로 다른 날짜, 다른 효소 분해 시약, 서로 다른 시퀀싱 레인(Sequencing Lane) 간에 발생하는 체계적 기술적 변동이다.

```
[Confounding Scenarios in Clinical Cohorts]

Scenario 1: Complete Confounding (Fatal Flaw)
  Batch 1:  All 10 Pre-treatment Samples
  Batch 2:  All 10 Post-treatment Samples
  * Treatment effect and Batch effect are mathematically colinear.
  * Impossible to separate biological signal from technical artifact!

Scenario 2: Randomized Block Design (Golden Standard)
  Batch 1:  5 Pre-treatment  +  5 Post-treatment
  Batch 2:  5 Pre-treatment  +  5 Post-treatment
  * Full rank design matrix: log(mu) = beta_0 + beta_Treat * X_Treat + beta_Batch * X_Batch
```

### 2.1 단순 배치 보정 행렬 사용의 함정 (The Danger of Corrected Counts)

Seurat의 `IntegrateData`, Harmony, 또는 Scanpy의 `bbknn`과 같은 배치 통합 도구는 **시각화(UMAP/t-SNE) 및 이웃 그래프 구축을 위한 보정 좌표**를 제공할 뿐이다. 
- 보정된 발현 행렬(Corrected matrix)은 카운트의 음이항 정수 특성을 파괴하고, 잔차 분산을 인위적으로 왜곡시킨다.
- **원칙**: 차등 발현(DEG) 및 차등 풍부도(DA) 검정 시에는 반드시 **원시 카운트(Raw Counts)를 사용하고, 배치 변수를 GLM 설계 행렬(Design Matrix) 내부에 공변량으로 직접 포함**시켜야 한다.

### 2.2 대리 변수 분석(SVA) 및 RUVseq를 통한 미관측 배치 통제

배치 정보가 기록되지 않았거나 다중 기술 변수가 복잡하게 얽힌 경우, `RUVseq` (Remove Unwanted Variation) 또는 `sva` (Surrogate Variable Analysis)를 적용한다:
1. **음성 대조 유전자군(Negative Control Genes, $\mathcal{G}_{\text{neg}}$) 선정**: 실험 조건에 의해 변하지 않는 하우스키핑 유전자 집합 선정.
2. **원치 않는 변동 축($W$) 추출**:
   $$Y_{\text{neg}} = X \boldsymbol{\beta}_{\text{neg}} + W \boldsymbol{\alpha}_{\text{neg}} + \boldsymbol{\epsilon}$$
   잔차 행렬에 특이값 분해(SVD)를 적용하여 주요 잠재 교란 벡터 $\mathbf{w}_1, \dots, \mathbf{w}_K$를 추출.
3. **최종 모델 적합**:
   $$\log(\mu_{ij}) = \mathbf{x}_{ij}^T \boldsymbol{\beta} + \sum_{k=1}^K w_{ik} \alpha_k + \log(s_{ij})$$

---

## 3. 독립 가설 가중치 (Independent Hypothesis Weighting, IHW)의 수학적 원리

단일세포 RNA-seq에서는 통상 15,000 ~ 25,000개의 유전자를 동시에 검정한다. 표준 Benjamini-Hochberg (BH) 절차는 모든 유전자에 균일한 사전 가중치를 부여하므로, 기저 발현량이 극도로 낮아 애초에 통계적 검정력이 거의 없는 수천 개의 유전자가 유효 가설들의 $p$-값을 희석시키는 심각한 다중 검정 손실을 초래한다.

```
[Independent Hypothesis Weighting (IHW) vs Standard BH]

Standard BH Procedure:
  All 20,000 genes get equal weight (w_g = 1).
  Low-expressed noise genes consume error budget and inflate adjusted p-values!

IHW Procedure (Ignatiadis et al., Nat Methods 2016):
  Stratify genes by independent informative covariate X_g (e.g., Mean Base Expression).
  Allocate weights w_g via Convex Optimization:
    Maximize Expected Discoveries subject to sum(w_g) = G and w_g >= 0.
  High-powered strata get higher weights -> Rescues genuine biological signals!
```

### 3.1 IHW 목적 함수 및 볼록 최적화 (Convex Optimization)

각 유전자 $g$에 대해, 영가설($H_0$) 하에서 $p$-값과 통계적으로 독립이면서 대립가설($H_1$) 하의 검정력과 강하게 상관된 공변량 $X_g$(예: 공여자 전체에서의 유전자 평균 카운트 $\bar{Y}_g$)를 정의한다.

유전자들을 공변량 $X_g$에 따라 $B$개의 층(Strata)으로 분할한다 ($b = 1, \dots, B$). 각 층의 가중치를 $w_b$라 할 때, 최적화 문제는 다음과 같이 정식화된다:

$$\max_{w_1, \dots, w_B} \sum_{b=1}^B \sum_{g \in \text{Stratum } b} \mathbb{E}\left[ \mathbb{I}\left( P_g \le \frac{t \cdot w_b}{G} \right) \right]$$

$$\text{subject to } \quad \frac{1}{G} \sum_{b=1}^B G_b \, w_b = 1, \quad w_b \ge 0 \; (\forall b)$$

(여기서 $G_b$는 층 $b$의 유전자 수, $G = \sum_b G_b$이며, 유전자 단위로 환산하면 $\sum_g w_g = G$이다.)

- 검정력이 0에 가까운 초저발현 층은 $w_b \to 0$으로 배분되어 오류 예산(Error budget) 낭비를 방지한다.
- 검정력이 풍부한 중·고발현 층은 $w_b > 1$을 부여받아 유의성 기각 역치가 완화된다.
- **영가설 하 독립성 조건(Independence under the Null)**: 유전자 $g$가 실제 차등 발현되지 않는 경우($H_{0, g}$가 참), 공변량 $X_g$의 크기와 무관하게 $P_g \sim \text{Uniform}(0, 1)$이어야 한다. 이 조건이 성립하면 IHW는 전체 위발견율(FDR)을 공칭 수준 $\alpha$로 통제한다(Ignatiadis et al. 2016).

---

## 4. 단일세포 통계 분석 10대 실패 요인 및 해결 지침

| 번호 | 흔한 실패 요인 (Pitfall) | 통계적 결과 | 올바른 해결 지침 (Best Practice) |
| :---: | :--- | :--- | :--- |
| **1** | 개별 세포를 독립 표본으로 취급하여 단순 t-test/Wilcoxon 적용 | 의사복제(Pseudoreplication)로 인한 $p < 10^{-50}$ 위양성 폭발 | 반드시 생물학적 반복(공여자)을 반영하는 Pseudobulk 또는 GLMM(NEBULA) 적용 |
| **2** | UMI 데이터에 제로-인플레이션 음이항(ZINB) 모델 오용 | 모수 과다 추정 및 불필요한 모형 복잡도 증가 | UMI 기반 10x 데이터의 0은 대부분 샘플링 영(zeros)이므로 표준 NB 모델 사용 |
| **3** | 보정된 발현 행렬(Integrated/Batch-corrected)을 DEG 모델에 입력 | 잔차 분산 왜곡 및 비선형 스케일링으로 인한 모수 왜곡 | 원시 카운트(Raw counts)를 입력하고 배치는 디자인 매트릭스 공변량으로 투입 |
| **4** | 분화 궤적을 임의의 클러스터로 잘라 두 그룹 간 비율 검정 수행 | 과도기 동역학 손실 및 해상도 편향 발생 | Milo의 $k$-NN 그래프 기반 이웃 DA 또는 GAMM 궤적 연속 검정 적용 |
| **5** | 심플렉스 제약을 무시하고 세포 비율에 단순 선형 회귀 적용 | 허위 음의 상관으로 인한 가짜 면역 고갈 선언 | scCODA 베이지안 디리클레-다항 회귀 또는 CLR/ilr 변환 적용 |
| **6** | $N=3$ 소규모 코호트에서 수천 개 DEG의 절대적 임상 유효성 주장 | 환자 특이적 아티팩트의 과적합 | 발견 코호트(Discovery)와 독립 검증 코호트(Validation)의 교차 검증 필수 |
| **7** | 다중 검정 보정(FDR) 누락 또는 유의미하지 않은 유전자 자의적 선별 | 재현 불가능한 p-hacking 발생 | Benjamini-Hochberg, Spatial FDR, 또는 IHW를 통한 엄격한 다중 검정 통제 |
| **8** | 샘플 간 총 포획 세포 수 편차($N_s$)에 대한 오프셋(Offset) 누락 | 시퀀싱 수율 차이가 생물학적 풍부도 변화로 둔갑 | GLM 적합 시 반드시 $\log(\text{Total Cells}_s)$를 오프셋으로 지정 |
| **9** | 계통 분기점 분석 시 전구체 구간 동등성 제약 미부여 | 분기 이전 구간에서의 거짓 차이 및 분기점 탐색 실패 | 분기 트리 위계 모델(Shared progenitor constraint) 적용 |
| **10** | 세포 주기 및 미토콘드리아 전사체 등 세포 수준 교란 요인 방치 | 대사 스트레스 유전자가 약제 저항 마커로 오인 | 세포 수준 GLMM(NEBULA)을 통해 연속형 공변량을 동시 회귀 보정 |

---

## 5. 5단계 단일세포 통계 캡스톤 프로젝트 평가 루브릭

거세 저항성 전립선암(mCRPC) 환자 코호트의 단일세포/단일핵 전사체 데이터를 종합 분석하는 연구자를 위한 5단계 표준 평가 루브릭(Capstone Rubric)이다.

```
[5-Stage Capstone Statistical Workflow]

Stage 1: Experimental Design & Power Analysis
         Verify biological replicates (N from power analysis; >=15 is illustrative) & design matrix full rankness.
                         │
                         ▼
Stage 2: Hierarchical DEG Architecture (NEBULA / Pseudobulk)
         Select GLMM vs Pseudobulk based on continuous cell-level confounders.
                         │
                         ▼
Stage 3: Continuous Manifold Trajectory Modeling (GAMM)
         Fit penalized B-splines along DPT/Latent time & test bifurcation contrasts.
                         │
                         ▼
Stage 4: High-Resolution Differential Abundance (Milo)
         Tile kNN graph with overlapping neighborhoods & apply Spatial FDR.
                         │
                         ▼
Stage 5: Global Simplex Compositional Modeling (scCODA)
         Anchor with homeostatic reference & run Bayesian Dirichlet-Multinomial MCMC.
```

### 루브릭 상세 기준표

아래 임계값($N \ge 15$ 등)은 교육용 예시이며, 실제 기준은 예상 효과 크기와 공여자 간 분산에 기반한 검정력 분석으로 정한다.

| 단계 | 핵심 평가 항목 | 최우수 (Distinction) | 보통 (Proficient) | 불합격 (Unsatisfactory) |
| :---: | :--- | :--- | :--- | :--- |
| **1단계** | 실험 설계 및 검정력 분석 | 공여자 분산 성분 분해를 수행하고 $N \ge 15$ 환자 확보, 배치와 조건이 직교 배치됨 | 환자 수가 8~14명 수준이나 공변량 모델링으로 교락을 적절히 보정함 | $N < 5$에서 세포 수만 5만 개 이상 확보하고 의사복제 모델을 계획함 |
| **2단계** | 위계적 차등 발현 분석 | 세포 수준 공변량(G2M, mito)과 공여자 임의 효과를 결합한 NEBULA 모델 적합 성공 | 세포 수준 공변량이 중요하지 않은 경우 적절한 Pseudobulk DESeq2 워크플로우 수행 | 개별 세포에 단순 Wilcoxon/Seurat FindMarkers를 적용하여 $p$-값 인플레이션 발생 |
| **3단계** | 궤적 기반 동역학 모델링 | GAMM을 통해 환자별 Random Smooths를 분리하고 동시 신뢰 띠로 분기점 $t^*$ 정확 추정 | 단순 GAM으로 궤적 곡선을 적합하고 Wald test로 유의 유전자 도출 | 분화 궤적을 3개 클러스터로 강제 분할하여 군집 간 DEG 수행 |
| **4단계** | 국소 차등 풍부도 검출 | Milo 이웃 그래프($k=30$) 구축, 총 세포 수 오프셋 투입, Spatial FDR 적용 | Milo 분석을 수행하였으나 Spatial FDR 대신 단순 BH를 적용함 | 클러스터 분할표를 작성하여 피셔 정확 검정으로 세포 비율 비교 |
| **5단계** | 전역 조성 데이터 모델링 | scCODA 자동 앵커링으로 최적 기준 세포형 선정, Spike-and-slab MCMC 수렴 완료 | speckle의 arcsin 변환을 활용하여 환자 반복을 limma로 적절히 검정함 | 심플렉스 제약을 무시하고 종양 세포 비율 폭증을 면역 세포 고갈로 단순 결론 |

---

## 6. 실전 코드 파이프라인

### 6.1 Python 환경: 공여자 수 대 세포 수 검정력 시뮬레이션

```python
"""
Python Simulation: Power Analysis comparing Donor Replication vs Cell Depth
"""
import numpy as np
from scipy import stats

def simulate_power(
    n_sims: int = 200,
    n_donors: int = 10,
    cells_per_donor: int = 1000,
    sigma_donor: float = 0.6,
    sigma_cell: float = 1.0,
    effect_size: float = 0.5,
    alpha: float = 0.05
) -> float:
    """
    주어진 공여자 수와 세포 수 하에서 공여자 수준 슈도벌크(합산) 로그 변환 t-검정의 검정력을 계산.
    (공여자 합에 Poisson GLM을 쓰면 공여자 간 과분산을 무시해 귀무 하 기각률이 크게 팽창하므로 사용하지 않음)
    """
    significant_count = 0
    
    for _ in range(n_sims):
        # 1. 공여자 수준 임의 효과 및 조건 할당 (1:1 비율)
        half_d = n_donors // 2
        donor_conditions = np.array([0] * half_d + [1] * half_d)
        donor_effects = np.random.normal(0, sigma_donor, size=n_donors)
        
        # 2. 세포 수준 카운트 시뮬레이션
        donor_idx = np.repeat(np.arange(n_donors), cells_per_donor)
        treat_vec = donor_conditions[donor_idx]
        
        # 선형 예측자 (Poisson-Log-Normal 근사)
        eta = 2.0 + effect_size * treat_vec + donor_effects[donor_idx] + \
              np.random.normal(0, sigma_cell, size=len(donor_idx))
        mu = np.exp(eta)
        counts = np.random.poisson(mu)
        
        # 3. 슈도벌크 집계 (Donor-level Sum) 후 로그 변환
        donor_sums = np.array([np.sum(counts[donor_idx == j]) for j in range(n_donors)])
        log_sums = np.log(donor_sums + 1.0)
        
        # 4. 공여자 수준 t-검정 (공여자가 실험 단위)
        _, p_val = stats.ttest_ind(log_sums[donor_conditions == 1], log_sums[donor_conditions == 0])
        if p_val < alpha:
            significant_count += 1
            
    power = significant_count / n_sims
    return power

# 시뮬레이션 실행 (effect_size=0.0이면 귀무 하 Type I error를 확인; 값은 실행마다 달라지는 예시임):
# print("Design A (군당 3 Donors x 5,000 cells):", simulate_power(n_donors=6, cells_per_donor=5000))
# print("Design B (군당 30 Donors x 500 cells):", simulate_power(n_donors=60, cells_per_donor=500))
```

### 6.2 R 환경: IHW를 통합한 차등 발현 적응적 다중 검정

```r
# ==============================================================================
# R Pipeline: Integrating IHW with DESeq2 for Single-Cell Pseudobulk
# ==============================================================================
library(DESeq2)
library(IHW)
library(ggplot2)

# 1. DESeq2 객체 적합 (dds: Pseudobulk count matrix)
dds <- DESeq(dds)

# 2. IHW 가설 가중치 적용
# filterFun = baseMean (유전자별 평균 발현량을 공변량으로 투입)
ihw_res <- results(
  dds, 
  filterFun = ihw, 
  alpha     = 0.05
)

# 3. 표준 Benjamini-Hochberg와 IHW의 발견 유전자 수 비교
# (참고: DESeq2 results()는 기본으로 baseMean 기반 independent filtering을 적용하므로 IHW의 추가 이득은 데이터에 따라 작을 수 있다)
standard_res <- results(dds, pAdjustMethod = "BH")
n_sig_bh  <- sum(standard_res$padj < 0.05, na.rm = TRUE)
n_sig_ihw <- sum(ihw_res$padj < 0.05, na.rm = TRUE)

message(sprintf("Significant Genes (Standard BH): %d", n_sig_bh))
message(sprintf("Significant Genes (IHW Adaptive): %d (Gain: +%d genes)", 
                n_sig_ihw, n_sig_ihw - n_sig_bh))

# 4. 공변량 층(Strata)별 가중치 배분 시각화
plot(metadata(ihw_res)$ihwResult)
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 군당 3명 $\times$ 5,000개 세포 vs 군당 30명 $\times$ 500개 세포 실험군에서 분산 성분 분해 및 표준 오차 비교
연구자가 전립선암 골전이 코호트 단일핵 전사체 분석을 설계하고 있다.
- **가정**:
  - 공여자 간 생물학적 분산: $\sigma_D^2 = 0.5$
  - 세포 간 기술적/생물학적 분산: $\sigma_c^2 = 1.0$
  - 총 세포 수는 30,000개로 고정 ($= 2 N n$; 위 두 설계 모두 군당 $N \cdot n = 15,000$).
1. **설계 1** (군당 $N=3$ 환자, $n=5,000$ 세포)과 **설계 2** (군당 $N=30$ 환자, $n=500$ 세포)에 대해, 회귀 계수의 추정 분산 $\text{Var}(\hat{\beta}) = \frac{2\sigma_D^2}{N} + \frac{2\sigma_c^2}{N \cdot n}$ 및 표준 오차 $\text{SE}(\hat{\beta})$를 각각 계산하시오.
2. 실제 참 효과 크기(True Effect Size)가 $\beta^* = 0.3$일 때, 두 설계에서 Wald 통계량 $W = \frac{\beta^*}{\text{SE}(\hat{\beta})}$의 크기를 비교하고, 왜 설계 1에서 검정력이 완전히 붕괴하는지 수리적으로 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 두 설계의 분산 및 표준 오차 정량 계산**:
- 군당 세포 수는 $N \cdot n = 15,000$으로 같다.
- 세포 수준의 오차 성분은 두 설계 모두 동일하다:
  $$\frac{2\sigma_c^2}{N \cdot n} = \frac{2 \cdot 1.0}{15,000} \approx 0.000133$$

- **설계 1 (군당 $N=3, n=5,000$)**:
  $$\text{Var}(\hat{\beta})_1 = \frac{2 \cdot 0.5}{3} + 0.000133 = 0.333333 + 0.000133 = 0.333467$$
  $$\text{SE}(\hat{\beta})_1 = \sqrt{0.333467} \approx 0.5775$$

- **설계 2 (군당 $N=30, n=500$)**:
  $$\text{Var}(\hat{\beta})_2 = \frac{2 \cdot 0.5}{30} + 0.000133 = 0.033333 + 0.000133 = 0.033467$$
  $$\text{SE}(\hat{\beta})_2 = \sqrt{0.033467} \approx 0.1829$$

- **결과 비교**: 설계 2의 표준 오차는 설계 1 대비 $\frac{0.1829}{0.5775} \approx 0.3167$, 즉 **약 3.16배 작다(정밀도가 3배 이상 높다)**.

**2. Wald 통계량 및 검정력 붕괴 원인**:
- **설계 1의 Wald 통계량**:
  $$W_1 = \frac{0.3}{0.5775} \approx 0.52$$
  표준정규분포에서 $Z = 0.52$에 대한 양측 $p$-값은 약 $0.60$이다. 통계적 유의수준 $\alpha = 0.05$ ($Z_{\text{crit}} = 1.96$)에 턱없이 미치지 못하므로 귀무가설을 기각할 확률(검정력)이 바닥 수준(정규 근사 약 8%, $t$-분포 $df=4$에서는 약 7%)에 머문다.
- **설계 2의 Wald 통계량**:
  $$W_2 = \frac{0.3}{0.1829} \approx 1.64$$
  양측 검정 기준 정규 근사 검정력은 약 38%($t$-분포 $df=58$에서는 약 36%)로 설계 1보다 훨씬 높지만 아직 충분하지 않으며, 80% 검정력에는 효과 크기가 더 크거나 공여자가 더 필요하다.
- **원인 분석**:
  세포를 세포당 10,000개로 깊게 시퀀싱하더라도, 공여자 간 분산 항 $\frac{2\sigma_D^2}{N}$이 분산의 $99.96\%$(설계 1)를 차지한다. 군당 환자가 3명에 불과하면 환자 간 차이를 평균화(Averaging out)할 수 있는 표본이 부족하고 검정의 자유도도 $df=4$에 그쳐 검정력이 크게 떨어진다.

</details>

---

### 문제 2: 독립 가설 가중치(IHW)에서 공변량의 영가설 하 독립성(Null Independence) 조건의 필연성
IHW 프레임워크에서 유전자별 가중치 $w_g$를 결정하기 위한 공변량(Informativeness Covariate $X_g$)을 선정할 때, "영가설($H_0$) 하에서 공변량 $X_g$는 $p$-값 $P_g$와 통계적으로 상호 독립이어야 한다"는 조건이 필수적이다.
1. 만약 연구자가 공변량 $X_g$로 '단순 치료군과 대조군 사이의 발현량 차이의 절대값($|\bar{Y}_{\text{Treat}} - \bar{Y}_{\text{Ctrl}}|$)을 잘못 선택했다고 가정하자. 영가설($H_0: \beta_g = 0$) 하에서 이 공변량이 $P_g$와 독립성을 상실하는 수리적 이유를 설명하시오.
2. 이러한 잘못된 공변량을 IHW 알고리즘에 투입했을 때, 볼록 최적화가 가중치 $w_g$를 어떻게 왜곡 배분하며, 왜 전체 위발견율(FDR) 통제가 붕괴하는지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 영가설 하 독립성 붕괴 수리적 증명**:
- 영가설 $H_0: \mu_{\text{Treat}} = \mu_{\text{Ctrl}}$이 참일 때, 두 표본 평균의 차이 $D_g = \bar{Y}_{g, \text{Treat}} - \bar{Y}_{g, \text{Ctrl}}$는 기댓값이 0인 정규분포를 따른다:
  $$D_g \sim \mathcal{N}\left(0, \sigma_{\text{diff}}^2\right)$$
- Wald 검정 통계량의 크기는 $|Z_g| = \frac{|D_g|}{\sigma_{\text{diff}}}$이다.
- 양측 $p$-값은 $P_g = 2 \left( 1 - \Phi(|Z_g|) \right) = 2 \left( 1 - \Phi\left( \frac{|D_g|}{\sigma_{\text{diff}}} \right) \right)$로 완벽히 $|D_g|$의 단조 감소 함수이다.
- 따라서 제안된 공변량 $X_g = |D_g|$는 영가설 하에서도 $P_g$와 완벽한 음의 결정론적 상관관계를 갖는다 ($X_g$가 크면 $P_g$는 필연적으로 0에 가까움).
- 즉, $P(P_g \le \alpha \mid X_g) \ne \alpha$가 되어 영가설 하 독립성 가정(Independence under the null)이 완전히 파탄 난다.

**2. 가중치 왜곡 및 FDR 붕괴 기전**:
- IHW의 목적 함수는 기각 가능한 가설 수를 최대화하는 것이다.
- 알고리즘은 $X_g = |D_g|$가 큰 유전자 집단이 가장 작은 $p$-값을 갖고 있음을 감지하고, 이 층(Stratum)에 가중치 예산 $w_b$를 사실상 몰아주게 된다($w_{\text{high}} \gg 1$, $w_{\text{low}} \to 0$).
- 영가설이 참인 무작위 노이즈 유전자들 중에서도, 표본 추출의 우연으로 인해 $|D_g|$가 크게 튄 유전자들이 높은 가중치를 배정받게 된다.
- 보정 $p$-값 $P_g^* = P_g / w_g$가 극단적으로 축소되어, 실제로는 참인 귀무가설들이 대거 기각 임계치를 통과하게 된다.
- 이는 통계적 가설 검정의 근본인 '데이터를 이용한 p-hacking'을 알고리즘 수준에서 체계적으로 수행하는 꼴이 되며, 결과적으로 실제 위발견율(True FDR)이 공칭 수준(e.g., 5%)을 크게 초과하게 된다(초과 정도는 설명용이며 계산된 수치가 아니다).
- 반면 올바른 공변량인 '전체 표본 평균 발현량($\bar{Y}_g$)'은 영가설 하에서 평균 발현량이 높든 낮든 $p$-값은 균일분포 $\mathcal{U}(0, 1)$를 따르므로 독립성이 (NB 검정에서 근사적으로) 성립한다.

</details>

---

### 문제 3: 골전이 전립선암 코호트의 부분 배치 교락(Partial Confounding) 시 RUVseq 기반 잔차 통제 설계
전립선암 골전이 환자 24명(치료 반응군 12명, 저항군 12명)의 생검 조직을 시퀀싱하였다.
- 배치 편성:
  - 배치 1: 반응군 9명 + 저항군 3명
  - 배치 2: 반응군 3명 + 저항군 9명
  - (배치가 기록되어 있으므로 1순위는 배치를 GLM 공변량으로 포함하는 것이며, RUV는 기록되지 않은 변동을 위한 보조·민감도 분석이다.)
치료 반응과 시퀀싱 배치 사이에 통계적으로 유의한 연관성(Partial Confounding, $\chi^2$ $p < 0.05$)이 존재한다.
1. 배치 효과를 무시하고 치료 반응만으로 모델을 적합했을 때 발생하는 위험과, 반대로 배치 통합 알고리즘(예: ComBat)으로 카운트 행렬을 사전에 전역 보정한 뒤 DEG를 수행했을 때 발생하는 과다 보정(Over-correction) 위험을 비교하시오.
2. `RUVseq`의 `RUVr` (Residual RUV) 기법을 적용하여 1차 모델의 잔차(Residuals)로부터 기술적 교란 벡터 $W_1$을 추정하고, 이를 최종 음이항 GLM 모델의 설계 행렬에 통합하는 3단계 통계 절차를 수식과 함께 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 두 접근법의 통계적 위험 비교**:
- **배치 무시 모델 (Unadjusted Model)**:
  배치 1과 배치 2 사이에 존재하는 체계적 시퀀싱 깊이 차이나 효소 소화 효율 차이가 저항성군과 반응군 사이의 발현량 차이로 잘못 전이된다. 배치 1에서 높게 측정되는 모든 유전자가 반응군 특이적 유전자로 오인되어 거대한 위양성이 발생한다.
- **사전 전역 보정 (Pre-transformation via ComBat / Integration)**:
  배치 변수와 생물학적 조건이 부분적으로 상관(Partial Confounding)되어 있으므로, 배치 간 차이를 제거하는 수리적 과정에서 배치와 공변하는 **진정한 생물학적 치료 저항성 신호의 상당 부분이 배치 변이로 간주되어 함께 제거(Over-correction)**된다. 이로 인해 통계적 검정력이 급감하고 실제 핵심 유전자가 소실되는 심각한 위음성이 발생한다.

**2. RUVr 기반 잔차 잠재 변수 추출 및 통합 3단계 워크플로우**:
- **1단계: 1차 생물학적 모델 적합 및 잔차 추출**:
  원시 카운트 행렬에 대해 관심 생물학적 변수만을 포함한 1차 음이항 GLM을 적합한다:
  $$\log(\mu_{ij}^{(0)}) = \beta_0 + \beta_{\text{Treat}} X_{j}^{\text{Treat}} + \log(s_{ij})$$
  적합된 기댓값 $\hat{\mu}_{ij}^{(0)}$으로부터 로그 척도의 데비앙스 잔차(Deviance Residuals) 또는 피어슨 잔차 행렬 $R \in \mathbb{R}^{G \times N}$을 계산한다:
  $$r_{gj} = \frac{y_{gj} - \hat{\mu}_{gj}^{(0)}}{\sqrt{\hat{\mu}_{gj}^{(0)} + \hat{\phi}_g (\hat{\mu}_{gj}^{(0)})^2}}$$
  이 잔차 행렬 $R$에는 치료 효과에 의해 설명되지 않은 순수 기술적 배치 변동과 비설명 잡음만이 잔류한다.

- **2단계: 특이값 분해(SVD)를 통한 잠재 교란 벡터 $\mathbf{w}$ 도출**:
  잔차 행렬 $R$에 대해 특이값 분해(SVD)를 수행한다:
  $$R = U \Sigma V^T$$
  우특이벡터 행렬 $V$의 첫 번째 열(또는 상위 $K$개 열)을 추출하여, 각 환자 샘플의 미관측 배치 효과를 대변하는 잠재 대리 변수 벡터 $\mathbf{w}_1 \in \mathbb{R}^N$을 도출한다.

- **3단계: 최종 통합 GLM 모델 적합 및 가설 검정**:
  추출된 잠재 벡터 $\mathbf{w}_1$을 원래의 원시 카운트 음이항 GLM의 설계 행렬에 공변량으로 투입하여 최종 적합한다:
  $$\log(\mu_{gj}) = \beta_{0, g} + \beta_{\text{Treat}, g} X_{j}^{\text{Treat}} + \alpha_{1, g} w_{1, j} + \log(s_{gj})$$
  이 모델에서 $\beta_{\text{Treat}, g}$에 대한 Wald 또는 Likelihood Ratio 검정을 수행한다. 단, RUVr의 $W$는 치료 효과를 회귀로 제거한 잔차에서 추출하므로 치료와 상관된 배치 성분은 포착하지 못해 치료-배치 교락을 제거하지 못한다. 기록된 배치는 설계 행렬에 직접 포함하는 것이 우선이며(부분 교락이면 설계행렬이 풀랭크여서 추정 가능), RUV는 미관측 변동에 대한 보조 분석이다. ComBat도 condition을 공변량으로 지정하면 생물학적 신호를 보존할 수 있으므로 과다 보정은 사용 방식에 의존한다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Ignatiadis N, Klaus B, Zaugg JB, et al. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods*. 2016;13(7):577-80. PMID: 27240256. doi: 10.1038/nmeth.3885
### 검수 근거 문헌
1. Choi K, Chen Y, Skelly DA, et al. Bayesian model selection reveals biological origins of zero inflation in single-cell transcriptomics. *Genome Biol*. 2020;21(1):183. PMID: 32718323. doi: 10.1186/s13059-020-02103-2
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- DESeq2 (bioconductor.org/packages/DESeq2). 원 논문: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15(12):550. PMID: 25516281. doi: 10.1186/s13059-014-0550-8
- IHW (bioconductor.org/packages/IHW). 원 논문: Ignatiadis N, Klaus B, Zaugg JB, et al. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods*. 2016;13(7):577-80. PMID: 27240256. doi: 10.1038/nmeth.3885
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
- ggplot2 (ggplot2.tidyverse.org).
<!-- REFERENCES:END -->

# 06. 분화 궤적 기반 연속 차등 발현 분석 (Trajectory Differential Expression & Lineage Bifurcation Testing)

---

## 학습 목표 (Cognitive Goals)

1. **연속적 세포 가소성(Continuous Plasticity)과 이산 클러스터(Discrete Clusters) 가정의 한계 규명**: 약제 저항성 획득 및 계통 전이 과정에서 인위적 군집화가 초래하는 위음성 및 경계 아티팩트를 비판적으로 분석하고, 연속 유사시간(Pseudotime $t_i$) 축을 따른 비선형 동역학 모델링의 수학적 필연성을 이해한다.
2. **벌점화 3차 B-스플라인(Penalized Cubic B-splines) 기반 일반화 가법 모델(GAM) 정식화**: 조도 페널티(Roughness Penalty $\lambda \int [f''(t)]^2 dt$)와 일반화 교차 검증(GCV) / 제한 최대 우도(REML)를 통한 최적 평활화 매개변수 $\lambda$의 유도 과정을 수학적으로 정식화한다.
3. **환자 간 이질성을 흡수하는 일반화 가법 혼합 모델(GAMM)과 Random Smooths 유도**: 다중 공여자(Multi-donor) 코호트에서 환자별 기저선 편차(Random Intercept) 및 개체별 궤적 곡률 편차(Random Spline Slopes)를 분리하여 가짜 궤적 차등 발현(Spurious Trajectory DE)을 억제한다.
4. **전립선암 계통 분기점 대조 검정(Lineage Bifurcation Contrast Testing) 구축**: 안드로겐 수용체 양성(AR+) 선암종(Adenocarcinoma)에서 치료 유도성 신경내분비암(t-NEPC)과 이중 음성 전립선암(DNPC)으로 갈라지는 임계 분기점(Bifurcation Point $t^*$)을 수치적으로 탐색하고, 통계적 동시 신뢰 띠(Simultaneous Confidence Bands) 기반 가설 검정을 구현한다.

---

## 1. 이산 클러스터 기반 DEG 분석의 이론적 붕괴와 연속 궤적 가설

단일세포 생물학에서 세포 분화, 상피-중배엽 전이(EMT), 암세포의 약제 적응성 재프로그래밍은 비가역적 스텝 함수가 아니라 **연속적 매니폴드 상의 점진적 전사체 리모델링 과정**이다. 

```
[Discrete Clustering Failure vs Continuous Trajectory Modeling]

1. Discrete Clustering Paradox
   Progenitor  ───>  Transient State 1  ───>  Transient State 2  ───>  Resistant Clone
   [ Cluster 1 ]    │ [    Cluster 2    ]   │  [    Cluster 3    ]   │ [  Cluster 4  ]
                    ▲                       ▲                        ▲
             Artificial Cut 1        Artificial Cut 2         Artificial Cut 3
   * Discards intermediate transition dynamics.
   * Creates false boundary artifacts and inflates within-cluster variance.

2. Continuous Trajectory Differential Expression (GAM / GAMM)
   y_i ~ NB(mu_i, theta)
   log(mu_i) = beta_0 + f(t_i) + log(s_i)
   Smooth Spline Curve f(t) captures continuous kinetic activation, peak timing, and bifurcation!
```

### 1.1 이산 클러스터링의 세 가지 통계적 오류

1. **경계 절단 아티팩트(Boundary Cutting Artifacts)**:
   연속적인 발현 구배를 갖는 세포 집단을 $k$-평균(k-means)이나 라이덴(Leiden) 알고리즘으로 강제 분할하면, 분할 경계선 바로 양옆에 위치한 사실상 동일한 세포들이 서로 다른 군집으로 분리되어 인위적인 평균 차이가 발생한다.
2. **과도기 중간 상태(Transient Intermediate States)의 통계적 소멸**:
   약제 노출 직후 일시적으로 활성화되었다가 최종 저항 상태에서 다시 감소하는 펄스형 유전자(Pulse-like early-response genes; e.g., 초기 스트레스 반응 인자, 일시적 후성유전학적 리프로그래밍 인자)는 시작 클러스터와 종결 클러스터의 단순 비교(Two-group comparison)에서 완전히 상쇄되어 발견되지 않는다.
3. **분기점(Branching Point) 검출 불가능**:
   단일 전구세포에서 두 개 이상의 운명(Fate)으로 갈라지는 결정적 분기점에서 어떤 전사 인자가 선제적으로 활성화되는지 시간적 인과관계(Temporal Ordering)를 규명할 수 없다.

### 1.2 Track 06과의 연계: 확산 유사시간($t_i$)과 범용 잠재 시간($\tau_i$)의 입력

[Track 06 (Manifold & Trajectory)](../../06_manifold_trajectory/Manifold_Trajectory_Master_Course.html)에서 유도된 두 가지 시간 변수를 본 통계 모델의 독립 변수 $t_i$로 직접 사용한다:
- **확산 유사시간(Diffusion Pseudotime, DPT $t_i \in [0, 1]$)**: 전이 확률 행렬의 다중 스케일 무작위 보행(확산 거리) 적분으로부터 유도된, 루트 세포 기준 확산 거리 기반의 상대 순서(측지선 거리가 아님).
- **범용 잠재 시간(Universal Latent Time, scVelo/CellRank $\tau_i \in [0, 1]$)**: RNA 스플라이싱 동역학(Splicing Kinetics) 및 전사-성숙-분해 미분방정식 모형으로부터 추정된, 모형 의존적인 상대적 시간(실제 시계 시간이 아님).

---

## 2. 일반화 가법 모델(GAM)과 벌점화 3차 B-스플라인

세포 $i$의 유전자 발현 카운트 $y_i$는 음이항 분포를 따르며, 조건부 기댓값 $\mu_i$와 유사시간 $t_i$ 사이의 관계는 비모수적(Non-parametric) 평활 함수 $f(t)$를 통해 모델링된다:

$$y_i \sim \text{NB}(\mu_i, \theta)$$

$$\log(\mu_i) = \beta_0 + f(t_i) + \log(s_i)$$

여기서 $s_i$는 세포 $i$의 라이브러리 크기(Sequencing depth offset)이다.

```
[Penalized Spline Fitting Mechanics]

Data Points (t_i, y_i)
      │
      ▼
Basis Functions B_k(t)  ──>  f(t) = sum_k beta_k B_k(t)  (Linear in parameters beta)
      │
      ▼
Roughness Penalty Integral:  lambda * ∫ [f''(t)]^2 dt = lambda * beta^T S beta
      │
      ▼
Objective: Maximize Penalized Likelihood l(beta) - 0.5 * lambda * beta^T S beta
  - High lambda: Underfitting (Straight Line, High Bias)
  - Low lambda:  Overfitting (Wiggly Noise Fitting, High Variance)
  - Optimal lambda via GCV / REML: Balances Biological Dynamics vs Smoothness
```

### 2.1 3차 B-스플라인(Cubic B-splines) 기저 전개

유사시간 구간 $[0, 1]$ 사이에 $K-4$개의 내부 매듭(Interior Knots) $\xi_1 < \xi_2 < \dots < \xi_{K-4}$을 배치하면, 임의의 매끄러운 함수 $f(t)$는 $K$개의 3차 B-스플라인 기저 함수 $B_k(t)$의 선형 결합으로 정밀하게 근사된다:

$$f(t) = \sum_{k=1}^K \beta_k B_k(t) = \mathbf{B}(t) \boldsymbol{\beta}$$

여기서 $\mathbf{B}(t) = [B_1(t), B_2(t), \dots, B_K(t)]$는 기저 벡터이며, $\boldsymbol{\beta} \in \mathbb{R}^K$는 추정해야 할 스플라인 계수 벡터이다. B-스플라인 기저는 국소 지지(Compact Local Support) 특성을 가지므로 수치적으로 극히 안정적이다.

### 2.2 조도 페널티(Roughness Penalty)와 평활화 파라미터 $\lambda$

과적합(Overfitting)을 억제하고 생물학적으로 타당한 매끄러운 궤적 곡선을 얻기 위해, 함수의 2차 도함수의 제곱 적분(곡률, Curvature)을 페널티로 부과한다:

$$\mathcal{P}(f) = \int_0^1 [f''(t)]^2 dt = \int_0^1 \left( \sum_{k=1}^K \beta_k B_k''(t) \right)^2 dt = \boldsymbol{\beta}^T \mathbf{S} \boldsymbol{\beta}$$

여기서 페널티 행렬 $\mathbf{S} \in \mathbb{R}^{K \times K}$의 원소는 기저 함수의 2차 도함수 간의 내적으로 계산된다:

$$S_{kl} = \int_0^1 B_k''(t) B_l''(t) dt$$

벌점화 로그 유사우도(Penalized Log-Likelihood) 목적 함수는 다음과 같다:

$$\ell_{\text{pen}}(\boldsymbol{\beta}; \lambda) = \ell(\boldsymbol{\beta}; \mathbf{y}) - \frac{1}{2} \lambda \boldsymbol{\beta}^T \mathbf{S} \boldsymbol{\beta}$$

- $\lambda \to 0$: 데이터의 모든 노이즈에 과적합(Interpolation, 무한 자유도).
- $\lambda \to \infty$: 2차 도함수가 0으로 강제되어 단순 1차 선형 회귀($f(t) = \beta_0 + \beta_1 t$)로 수렴.
- **최적 $\lambda$ 추정**: 제한 최대 우도(Restricted Maximum Likelihood, REML) 또는 일반화 교차 검증(Generalized Cross Validation, GCV) 점수를 최소화하여 데이터 기반으로 최적 결정된다.

---

## 3. 일반화 가법 혼합 모델(GAMM): 다중 환자 코호트의 환자별 임의 곡선

단일세포 궤적 분석에서 가장 치명적인 함정은 **서로 다른 환자(Donor)에서 유래한 세포들을 하나의 풀(Pool)로 병합하여 단일 GAM 곡선을 적합하는 것**이다.

### 3.1 환자별 궤적 편차(Patient Trajectory Confounding)

실제 임상 데이터에서는 다음과 같은 현상이 발생한다:
1. **환자별 기저 발현 편차(Random Intercept)**: 특정 환자 $A$는 유전자 $X$의 기저 발현이 전반적으로 높고, 환자 $B$는 낮다.
2. **환자별 분화 곡률 편차(Random Smooths)**: 환자의 면역 미세환경이나 유전적 배경에 따라 유사시간에 따른 유전자 활성화의 기울기나 정점(Peak) 위치가 미세하게 다를 수 있다.

만약 이를 단일 GAM으로 풀면, 특정 환자의 세포가 궤적의 후반부에 우연히 몰려 있을 때(Sampling bias), 환자 간 기저 차이가 궤적에 따른 유전자 발현 변화로 둔갑하는 **가짜 궤적 DEG(Spurious Trajectory DEG)**가 대량 발생한다.

```
[Spurious DE from Uncontrolled Donor Effects along Pseudotime]

Pseudotime (t):   0.0 ──────────────> 0.5 ──────────────> 1.0
Donor A (High X): [ooo ooo ooo]       [           ]       [           ] (All in t < 0.3)
Donor B (Low X):  [           ]       [           ]       [xxx xxx xxx] (All in t > 0.7)

Naive GAM fit:
  Draws a sharp downward curve for Gene X!
  Conclusion: "Gene X is strongly downregulated during neuroendocrine differentiation."
  REALITY: Gene X is constant within both donors! Difference is 100% Patient Batch Effect!

GAMM Resolution:
  log(mu_ij) = beta_0 + f_pop(t_ij) + b_j + f_j(t_ij) + log(s_ij)
  f_pop(t) correctly estimates a flat horizontal line!
```

### 3.2 GAMM 수리적 구조

공여자 $j$ ($j=1, \dots, J$)의 세포 $i$에 대해 다음과 같이 위계적 GAMM을 수립한다:

$$\log(\mu_{ij}) = \beta_0 + f_{\text{pop}}(t_{ij}) + b_j + f_j(t_{ij}) + \log(s_{ij})$$

- $f_{\text{pop}}(t) = \sum_{k=1}^K \beta_k B_k(t)$: 전체 환자 집단(Population)을 관통하는 평균 분화 동역학 곡선.
- $b_j \sim \mathcal{N}(0, \sigma_b^2)$: 공여자 $j$의 기저선 임의 절편(Random Intercept).
- $f_j(t) = \sum_{m=1}^M u_{jm} \widetilde{B}_m(t)$: 공여자 $j$ 특이적 잔차 평활 곡선(Random Smooth). 계수 벡터 $\mathbf{u}_j \sim \mathcal{N}(\mathbf{0}, \sigma_u^2 \mathbf{S}^{-})$는 스플라인 페널티 행렬 $\mathbf{S}$를 정밀도(Precision)로 갖는 가우스 마르코프 랜덤 필드(GMRF)로 규정된다($\mathbf{S}$는 상수·선형 성분의 영공간 때문에 랭크 결손이므로 역행렬이 없고 $\mathbf{S}^{-}$는 유사역행렬/improper prior를 뜻하며, mgcv는 영공간을 별도 고정효과·패널티로 분리해 처리한다).

이 구조는 Wood (2017)의 mixed model representation of splines에 의해 선형 혼합 모델의 분산 성분(Variance Components) 추정 프레임워크로 전환될 수 있다.

---

## 4. 계통 분기점 대조 검정 (Lineage Bifurcation Contrast Testing)

전립선암에서 엔잘루타마이드 치료 저항성은 종종 계통 가소성(Lineage Plasticity)에 의해 유발된다. AR 신호전달이 차단된 선암종 세포는 공통 전구 상태를 거쳐 두 가지 치명적인 저항 경로로 분기한다고 가정한다(t-NEPC·DNPC 분류는 Labrecque 2019, Bluemn 2017 등 mCRPC 표현형 문헌에 근거하며, 이 분기 구조 자체는 분석 설계를 위한 모식적 가정이다):
1. **t-NEPC (Treatment-induced Neuroendocrine Prostate Cancer)**: 신경내분비 마커(`SYP`, `CHGA`, `ASCL1`, `NEUROD1`) 발현 획득, 급격한 증식.
2. **DNPC (Double-Negative Prostate Cancer)**: AR 음성 및 NE 음성, FGF-MAPK 수용체 과발현, EMT 전사 인자(`ZEB1`, `SNAI2`) 주도 중간엽 침윤.

```
[Lineage Bifurcation Architecture in Prostate Cancer]

                       ┌── Lineage 1: t-NEPC (ASCL1+, SYP+, Chromogranin A+)
                       │   Curve: f_1(t)
Common Progenitor ─────○ Bifurcation Point t*
(AR+ Adenocarcinoma)   │
Curve: f_prog(t)       └── Lineage 2: DNPC (FGF+, ZEB1+, VIM+, Stem-like)
                           Curve: f_2(t)

Hypothesis Testing along Pseudotime:
  Before t*: delta(t) = f_1(t) - f_2(t) = 0   (Lineages are statistically identical)
  After t*:  delta(t) = f_1(t) - f_2(t) != 0  (Significant divergence begins)
```

### 4.1 계통 간 차이 함수 $\delta(t)$ 및 동시 신뢰 띠(Simultaneous Confidence Bands)

두 분기 계통 $L_1$과 $L_2$에 대해 각각의 평활 곡선을 $f_1(t)$와 $f_2(t)$라 할 때, 차이 함수는 다음과 같다:

$$\delta(t) = f_1(t) - f_2(t) = \mathbf{B}_1(t) \hat{\boldsymbol{\beta}}_1 - \mathbf{B}_2(t) \hat{\boldsymbol{\beta}}_2 = \mathbf{C}(t) \hat{\boldsymbol{\beta}}$$

차이 함수의 점별 표준 오차(Pointwise Standard Error)는 모델의 공분산 행렬 $\mathbf{V}_{\boldsymbol{\beta}} = \text{Cov}(\hat{\boldsymbol{\beta}})$로부터 도출된다:

$$\text{SE}(\delta(t)) = \sqrt{\mathbf{C}(t) \mathbf{V}_{\boldsymbol{\beta}} \mathbf{C}(t)^T}$$

그러나 전체 유사시간 구간 $t \in [0, 1]$에 걸쳐 다중 검정(Multiple testing across continuous time)을 올바르게 수행하기 위해서는 점별(Pointwise) 신뢰구간이 아닌 **동시 신뢰 띠(Simultaneous Confidence Bands)**를 구축해야 한다.

임계값 $c_{\alpha}$는 다변량 정규 시뮬레이션(Multivariate Normal Simulation)을 통해 supremum 통계량을 통제하여 계산된다:

$$\sup_{t \in [0, 1]} \frac{|\mathbf{C}(t) (\hat{\boldsymbol{\beta}} - \boldsymbol{\beta})|}{\text{SE}(\delta(t))} \le c_{\alpha}$$

따라서 $100(1-\alpha)\%$ 동시 신뢰 띠는 다음과 같다:

$$\text{SCB}_{1-\alpha}(t) = \left[ \delta(t) - c_{\alpha} \text{SE}(\delta(t)), \; \delta(t) + c_{\alpha} \text{SE}(\delta(t)) \right]$$

### 4.2 임계 분기점 $t^*$의 수치적 추정

유전자 발현이 두 계통 간에 통계적으로 유의하게 분기하는 최초의 유사시간 $t^*$은 동시 신뢰 띠가 0을 완전히 배제하기 시작하는 지점으로 정의된다:

$$t^* = \inf \left\{ t \in [0, 1] : 0 \notin \text{SCB}_{1-\alpha}(t) \text{ and } 0 \notin \text{SCB}_{1-\alpha}(t') \; \forall t' \in [t, t+\epsilon] \right\}$$

- $t^*$의 조기 출현(Early divergence): 해당 유전자가 계통 결정을 유도하는 상위 구동 전사 인자(Lineage Driver TF; e.g., `ASCL1`)임을 시사하는 가설이다(의사시간은 발현 유사도로 추정한 순서이므로 조기 분기가 인과적 구동을 입증하지는 않으며, 후속 기능 실험이 필요하다).
- $t^*$의 후기 출현(Late divergence): 계통 결정 이후에 유도되는 하위 분화 실행 마커(Terminal Effector; e.g., `CHGA`)임을 시사.

---

## 5. 실전 R/Python 구현 파이프라인

> **주의**: 이 절의 `tradeSeq`(NB-GAM, 영과잉 가중치 제공)와 아래 Python NB-GAM은 환자 임의 효과가 없는 **풀링 모형**이므로, 3절이 경고한 공여자 효과 교란을 통제하지 못합니다. 환자 효과를 통제하려면 `mgcv`의 GAMM(예: `bam`에서 환자 임의 절편 `s(donor, bs="re")`와 환자별 평활 `s(t, donor, bs="fs")`)이나 유사시간 구간별 pseudobulk 접근을 사용해야 합니다. 또한 유사시간·계통 가중치가 같은 발현 데이터에서 추론되므로 $p$-value는 탐색적으로 해석합니다.

### 5.1 R 환경: `tradeSeq` 및 `mgcv`를 활용한 궤적 차등 발현 분석

```r
# ==============================================================================
# R Pipeline: Trajectory DE & Bifurcation Testing with tradeSeq / mgcv
# ==============================================================================
library(tradeSeq)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)

# 1. Slingshot 또는 CellRank로부터 도출된 유사시간 및 계통 가중치 준비
# sce: SingleCellExperiment 객체
# pseudotime_mat: N x 2 행렬 (Lineage 1: t-NEPC, Lineage 2: DNPC)
# cell_weights: N x 2 행렬 (각 계통에 속할 세포별 할당 확률)

# 2. 최적 매듭 수(Knots) 평가
# evaluateK 함수를 통해 BIC/AIC 기반 최적 B-스플라인 기저 수 선정 (일반적으로 K=5~7)
set.seed(42)
ic_mat <- evaluateK(
  counts      = counts(sce)[1:500, ], # 상위 500개 유전자로 고속 사전 평가
  pseudotime  = pseudotime_mat,
  cellWeights = cell_weights,
  k           = 3:8,
  nGenes      = 200
)

# 3. 음이항 GAM 모델 전체 적합 (fitGAM)
# UMI depth를 내재적 offset으로 자동 처리
gam_models <- fitGAM(
  counts      = counts(sce),
  pseudotime  = pseudotime_mat,
  cellWeights = cell_weights,
  nknots      = 6,
  parallel    = TRUE,
  BPPARAM     = BiocParallel::MulticoreParam(workers = 8)
)

# 4. 가설 검정 A: Association Test (유사시간을 따라 유의하게 발현이 변하는가?)
assoc_res <- associationTest(gam_models)
assoc_res$FDR <- p.adjust(assoc_res$pvalue, method = "BH")

# 5. 가설 검정 B: EarlyDE Test / DiffEnd Test (분기점에서 두 계통 간 발현 차이)
# knots = c(2, 4) 구간에서 t-NEPC와 DNPC 간 발현 궤적이 분기하는 유전자 검정
early_de_res <- earlyDETest(
  models     = gam_models,
  knots      = c(2, 4),
  l2fc       = 0.5      # 최소 log2 fold change 역치
)
early_de_res$FDR <- p.adjust(early_de_res$pvalue, method = "BH")

# 6. 상위 분기 유전자 시각화 (예: ASCL1 vs ZEB1)
top_branch_genes <- rownames(early_de_res[order(early_de_res$pvalue), ])[1:5]
plotSmoothers(gam_models, counts(sce), gene = top_branch_genes[1]) +
  ggtitle(paste0("Bifurcation Dynamics: ", top_branch_genes[1]))
```

### 5.2 Python 환경: `statsmodels` 및 스플라인 기저 생성을 통한 커스텀 GAM 적합

```python
"""
Python Implementation: Fitting Negative Binomial Spline GAM along Pseudotime
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrix
from scipy.interpolate import BSpline

def fit_trajectory_nb_gam(
    pseudotime: np.ndarray,
    gene_counts: np.ndarray,
    library_size: np.ndarray,
    n_knots: int = 6
) -> dict:
    """
    주어진 유전자에 대해 3차 B-스플라인 음이항 GAM을 적합하고 통계량을 반환.
    """
    # 1. 결측치 및 비유효 유사시간 필터링
    valid_mask = ~np.isnan(pseudotime) & (pseudotime >= 0) & (pseudotime <= 1)
    t_val = pseudotime[valid_mask]
    y_val = gene_counts[valid_mask]
    lib_val = library_size[valid_mask]
    
    # 2. 균일 내부 매듭(Knots) 생성 및 B-스플라인 기저 행렬(Design Matrix) 생성
    knots = np.linspace(0, 1, n_knots)[1:-1]
    # patsy의 bs를 사용하여 cubic B-spline 행렬 구축
    B_matrix = dmatrix(
        "bs(t, knots=knots, degree=3, include_intercept=False)",
        {"t": t_val, "knots": knots},
        return_type='dataframe'
    )
    
    # 3. 모델 프레임 구축 및 Offset(log library size) 추가
    df_model = B_matrix.copy()
    df_model['y'] = y_val
    df_model['offset'] = np.log(lib_val)
    
    # 4. Negative Binomial GLM 적합
    exog_cols = [c for c in df_model.columns if c not in ['y', 'offset']]
    X = sm.add_constant(df_model[exog_cols])
    
    # NB2 최대우도 적합: 과분산 alpha를 고정하지 않고 데이터로 추정 (마지막 파라미터가 alpha)
    model = sm.NegativeBinomial(
        endog=df_model['y'],
        exog=X,
        offset=df_model['offset']
    )
    fit_res = model.fit(disp=0, maxiter=200)
    
    # 5. 전역 스플라인 유의성 검정 (Wald Test on all spline coefficients)
    # H0: beta_1 = beta_2 = ... = beta_K = 0 (수평선인가?)
    n_params = len(fit_res.params) - 1  # 마지막은 alpha
    r_matrix = np.eye(len(fit_res.params))[1:n_params]  # 절편(const)과 alpha 제외
    wald_test = fit_res.wald_test(r_matrix, scalar=True)
    
    return {
        "params": fit_res.params,
        "pvalue_trajectory": wald_test.pvalue,
        "wald_statistic": wald_test.statistic,
        "aic": fit_res.aic
    }

# 실행 예제:
# res = fit_trajectory_nb_gam(adata.obs['dpt_pseudotime'].values, 
#                             adata[:, 'ASCL1'].X.toarray().flatten(), 
#                             adata.obs['total_counts'].values)
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: B-스플라인 기저 수 $K$와 평활화 페널티 $\lambda$의 과대/과소 적합(Over/Under-fitting)이 생물학적 전이 동역학 추정에 미치는 영향
단일세포 분화 궤적을 분석할 때, 유사시간 축에 배치하는 매듭의 수(Knots $K$)와 평활화 페널티 파라미터 $\lambda$의 튜닝 실패가 생물학적 해석에 미치는 통계적 영향을 논하시오.
1. 매듭의 수 $K$가 너무 작거나($K=3$) $\lambda$가 비정상적으로 과대하게 설정되었을 때, 전이 중간 상태에서 일시적으로 솟구치는 펄스형 전사 인자(Transient Pulse Dynamics)의 궤적 추정치에서 발생하는 편향(Bias)을 설명하시오.
2. 반대로 $K$가 너무 크고($K=30$) $\lambda \to 0$에 수렴할 때, 샘플링 영(zeros)에 의한 노이즈나 국소적 세포 밀도 변동에 의해 초래되는 분기점 검정의 위양성(Type I Error) 기전을 수식과 함께 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 과도한 평활화($K$ 과소, $\lambda$ 과대)에 의한 생물학적 편향**:
- 조도 페널티 $\lambda \int [f''(t)]^2 dt$가 극단적으로 커지면, 모델은 곡률을 0으로 강제하여 1차 선형 회귀 곡선($f(t) \to \beta_0 + \beta_1 t$) 또는 완만한 단조 곡선으로 축퇴된다.
- 전립선암 계통 전이 과정에서 특정 전사 인자(e.g., 일시적 스트레스 반응 인자, 조기 EMT 개시 인자)는 분화 유사시간 $t \approx 0.3 \sim 0.5$ 구간에서 급격히 상승했다가 $t \approx 0.7$ 이후 기저 상태로 복귀하는 비단조 펄스 형태($\cap$-shape)를 나타낸다.
- 과도하게 평활화된 모델은 이 고주파 국소 피크를 잔차 오차(Residual Noise)로 간주하고 평평하게 깎아버린다.
- 결과적으로 시작점($t=0$)과 종결점($t=1$)의 발현량이 유사할 경우, Wald 검정 통계량이 0에 수렴하여 해당 핵심 조절 인자가 궤적을 따라 유의하게 변동하지 않는다는 치명적인 위음성(Type II Error)을 초래한다.

**2. 과소 평활화($K$ 과대, $\lambda \to 0$)의 영향 (SE가 올바르게 계산되면 위양성이 아니라 검정력 저하)**:
- 페널티가 없는 고차원 B-스플라인은 유효 자유도(Effective Degrees of Freedom)가 기저의 수 $K$에 근접한다:
  $$\text{EDF} = \text{tr}\left( (\mathbf{B}^T W \mathbf{B} + \lambda \mathbf{S})^{-1} \mathbf{B}^T W \mathbf{B} \right) \xrightarrow{\lambda \to 0} K$$
- 단일세포 전사체 데이터는 많은 샘플링 영(zeros; UMI 데이터에서 0이 많은 것은 대부분 낮은 발현의 포아송 샘플링 결과)과 시퀀싱 깊이 변동성을 지닌다.
- 제약이 없는 스플라인 곡선은 실제 생물학적 추세가 아닌, 특정 유사시간 윈도우에 우연히 몰린 2~3개 세포의 높은 카운트 노이즈를 따라가기 위해 심하게 요동(Wiggle)친다.
- 두 계통의 발현 곡선 차이 $\delta(t) = f_1(t) - f_2(t)$를 계산할 때, 이 국소적 진동으로 인해 $f_1(t)$의 노이즈 정점과 $f_2(t)$의 노이즈 골짜기가 교차하는 지점에서 $|\delta(t)|$가 인위적으로 급증한다.
- SE를 올바르게 계산하면 과적합은 분산을 키워 동시 신뢰 띠가 넓어지므로 위양성보다는 검정력 저하로 이어진다. 허위 분기점(Spurious Bifurcation Point)이 무더기로 검출되는 위양성은 과적합에 더해 SE를 과소추정할 때(공여자 군집·과분산 미반영 등) 발생한다.

</details>

---

### 문제 2: 계통 분기 검정 시 전구체 구간($t < t^*$)에서의 동등성 제약(Constraint) 부여 여부에 따른 검정력 차이
AR 양성 선암종에서 t-NEPC와 DNPC로 갈라지는 분기 계통을 분석할 때, 두 계통의 궤적 모델을 수립하는 전략은 두 가지가 있다.
- **전략 A (비구속 독립 적합, Unconstrained Independent Fit)**: 두 계통의 데이터를 완전히 분리하여 각 계통마다 독자적으로 $f_1(t)$와 $f_2(t)$를 적합.
- **전략 B (분기 트리 위계 모델, Branching Tree Hierarchy with Shared Progenitor)**: 공통 전구 구간($t < t^*$)에서는 $f_1(t) \equiv f_2(t)$로 강제 구속하고, 분기점 $t^*$ 이후부터 독자적인 스플라인 기저를 전개.
1. 전략 A를 사용할 경우 공통 전구체 구간($t < t^*$)에서 생물학적으로 동일한 상태의 세포군에 대해 표본 분할(Data Splitting)로 인해 추정 분산이 어떻게 팽창하는지 설명하시오.
2. 전략 B가 왜 분기점($t^*$) 직후에서 계통 간 차이 $\delta(t)$를 검출하는 통계적 검정력(Statistical Power)을 월등히 향상시키는지 자유도(Degrees of Freedom) 및 모수 추정 정밀도 측면에서 논하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 전략 A의 표본 분할에 따른 추정 분산 팽창**:
- 전략 A는 $t < t^*$ 구간에 위치한 공통 전구 세포들을 계통 1과 계통 2의 할당 확률 가중치로 강제 분할하여 각각 독립된 회귀 모델에 투입한다.
- 스플라인 계수 추정치의 분산은 디자인 행렬의 정보량에 반비례한다:
  $$\text{Var}(\hat{\boldsymbol{\beta}}_1) \propto (\mathbf{B}_1^T W_1 \mathbf{B}_1)^{-1}$$
- 동일한 전구 세포 집단을 반으로 쪼개어 학습하므로 각 모델의 유효 세포 수가 절반으로 감소하며, 이는 B-스플라인 계수의 표준 오차를 $\sqrt{2}$배(약 41%) 증가시킨다.
- 더욱이 동일한 생물학적 시작점을 공유함에도 불구하고, 무작위 샘플링 오차로 인해 $t=0$에서의 추정 절편 $\hat{f}_1(0)$과 $\hat{f}_2(0)$이 불일치하는 모순(Starting Point Inconsistency)이 발생하며, 분기 이전 구간에서도 거짓 차이 신호가 감지될 위험이 커진다.

**2. 전략 B의 구속 모델을 통한 검정력 극대화 원리**:
- 전략 B는 $t \le t^*$ 구간에서 $f_1(t) = f_2(t) = f_{\text{prog}}(t)$라는 동등성 제약을 부여한다.
- 이로 인해 공통 전구체 구간의 모든 세포($N_{\text{prog}} = N_1 + N_2$)가 단일 기저 계수 $\boldsymbol{\beta}_{\text{prog}}$를 추정하는 데 온전히 통합 활용되므로, 추정량의 분산이 최소화되고 기저선의 신뢰 띠 폭이 극도로 좁아진다.
- 분기점 직후($t > t^*$)에서 계통 간 차이 $\delta(t) = f_1(t) - f_2(t)$를 검정할 때, 시작 기저선의 불확실성(Baseline Uncertainty)이 사전에 완벽히 제거되어 있으므로 분기 계수 $\boldsymbol{\beta}_{\text{diff}} = \boldsymbol{\beta}_1 - \boldsymbol{\beta}_2$의 표준 오차가 최소화된다:
  $$\text{SE}(\delta(t^* + \epsilon))_{\text{constrained}} \ll \text{SE}(\delta(t^* + \epsilon))_{\text{unconstrained}}$$
- 따라서 Wald 통계량 $W = \frac{\delta(t)}{\text{SE}(\delta(t))}$의 분모가 대폭 감소하여, 분기 초기 단계에서 발생하는 미세한 발현량 차이(Small Effect Size)를 조기에 통계적으로 유의하게 포획할 수 있는 검정력이 획기적으로 상승한다.

</details>

---

### 문제 3: 전립선암 환자 15명의 전치 궤적 데이터에서 Random Smooths 미적용 시 발생하는 위양성 분기 유전자 진단
전립선암 환자 15명으로부터 추출된 50,000개 세포의 단일세포 궤적 데이터를 분석하였다.
1. 환자 15명 중 단 1명의 특정 환자(환자 ID: PT-09)가 고유한 생식선 변이로 인해 유전자 $Y$를 전사체 전반에서 비정상적으로 과발현하고 있으며, 이 환자의 암세포들이 우연히 t-NEPC 분기 궤적의 말단에 집중 분포하고 있다. 일반화 가법 모델(단순 GAM)을 적용했을 때 유전자 $Y$가 t-NEPC 특이적 핵심 구동 유전자로 잘못 판정되는 수리적 과정을 설명하시오.
2. 일반화 가법 혼합 모델(GAMM)에서 환자별 임의 절편($b_j$) 및 임의 평활 곡선($f_j(t)$)을 포함시켰을 때, 이 위양성 오류가 어떻게 수학적으로 무력화되는지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 단순 GAM에서 샘플 불균형에 의한 위양성 생성 과정**:
- 단순 GAM은 모든 관측 세포를 동일한 확률 과정에서 독립적으로 추출된 것으로 가정한다:
  $$\log(\mu_i) = \beta_0 + f(t_i)$$
- 환자 PT-09의 세포들이 t-NEPC 궤적의 후반부($t > 0.8$)에 대거 밀집되어 있다면, 데이터셋에서 $t > 0.8$ 구간의 평균 카운트 $\bar{y}_{t > 0.8}$는 PT-09의 이상치 발현값에 의해 지배된다.
- 최적화 알고리즘은 전체 로그 유사우도를 최대화하기 위해 스플라인 함수 $f(t)$의 후반부를 강제로 위로 끌어올린다(Upward inflection).
- 모델은 세포 수($N \approx 50,000$)에 기반하여 스플라인 계수의 표준 오차를 계산하므로, $\text{SE}(\hat{\beta}_k) \propto 1/\sqrt{N}$이 되어 극도로 작은 표준 오차를 갖게 된다.
- 그 결과, 유전자 $Y$의 후반부 상승 곡선에 대해 $p < 10^{-50}$의 극단적 통계 유의성이 부여되며, 실제로는 단 1명의 환자 특이적 아티팩트(Private Patient Artifact)에 불과한 유전자가 질환 전체의 t-NEPC 핵심 드라이버로 둔갑한다.

**2. GAMM의 임의 효과 모델링을 통한 위양성 소거 기전**:
- GAMM은 세포 수준의 발현 기댓값에 환자 수준의 임의 효과를 분리 결합한다:
  $$\log(\mu_{ij}) = \beta_0 + f_{\text{pop}}(t_{ij}) + b_j + f_j(t_{ij}) + \log(s_{ij})$$
- PT-09의 전반적인 초과 발현량은 환자별 임의 절편 $b_{\text{PT-09}} \sim \mathcal{N}(0, \sigma_b^2)$에 의해 거대한 양의 값으로 우선 흡수된다.
- 나아가 궤적 상의 특이적 진동은 환자 특이적 임의 평활 함수 $f_{\text{PT-09}}(t)$로 분리된다.
- 집단 전체를 대표하는 고정 효과 궤적 곡선 $f_{\text{pop}}(t)$의 계수 $\boldsymbol{\beta}_{\text{pop}}$의 분산-공분산 행렬은 세포 수가 아니라 환자 수($J=15$)를 기준으로 하는 한계 공분산 행렬 $V$를 통해 추정된다:
  $$\text{Var}(\hat{\boldsymbol{\beta}}_{\text{pop}}) = \left( X^T V^{-1} X \right)^{-1}$$
- 단, 이 가정 시나리오처럼 PT-09의 세포가 궤적 말단($t > 0.8$)에 집중되어 다른 환자의 세포가 그 구간에 거의 없다면 말단 구간의 $f_{\text{pop}}(t)$는 PT-09 한 명으로만 식별되어 환자 효과와 궤적 효과를 분리할 수 없다(비식별/외삽). 말단 구간에 여러 환자의 세포가 존재할 때에만 나머지 환자에서 $Y$의 체계적 증가가 없음을 근거로 $f_{\text{pop}}(t)$가 평평하게 추정된다.
- 동시에 PT-09 1명에 기인한 환자 간 분산 $\sigma_b^2$이 크게 추정되어 계수의 표준 오차가 현실적인 크기로 복원되므로, Wald 검정 통계량이 임계값을 넘지 못해 위양성이 크게 완화된다(완벽한 방어는 아니며 말단 구간의 환자 다양성 확인과 민감도 분석이 필요하다).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Wood SN. Generalized Additive Models. 2017. doi: 10.1201/9781315370279
2. Bluemn EG, Coleman IM, Lucas JM, et al. Androgen Receptor Pathway-Independent Prostate Cancer Is Sustained through FGF Signaling. *Cancer Cell*. 2017;32(4):474-489.e6. PMID: 29017058. doi: 10.1016/j.ccell.2017.09.003
3. Labrecque MP, Coleman IM, Brown LG, et al. Molecular profiling stratifies diverse phenotypes of treatment-refractory metastatic castration-resistant prostate cancer. *J Clin Invest*. 2019;129(10):4492-4505. PMID: 31361600. doi: 10.1172/JCI128212
### 검수 근거 문헌
1. Van den Berge K, Roux de Bézieux H, Street K, et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat Commun*. 2020;11(1):1201. PMID: 32139671. doi: 10.1038/s41467-020-14766-3
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- patsy 1.0.3 (github.com/pydata/patsy).
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
- slingshot (bioconductor.org/packages/slingshot). 원 논문: Street K, Risso D, Fletcher RB, et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics*. 2018;19(1):477. PMID: 29914354. doi: 10.1186/s12864-018-4772-0
- statsmodels 0.15.0 (github.com/statsmodels/statsmodels).
- tradeSeq (bioconductor.org/packages/tradeSeq). 원 논문: Van den Berge K, Roux de Bézieux H, Street K, et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat Commun*. 2020;11(1):1201. PMID: 32139671. doi: 10.1038/s41467-020-14766-3
- BiocParallel (bioconductor.org/packages/BiocParallel).
- ggplot2 (ggplot2.tidyverse.org).
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- SingleCellExperiment (bioconductor.org/packages/SingleCellExperiment). 원 논문: Amezquita RA, Lun ATL, Becht E, et al. Orchestrating single-cell analysis with Bioconductor. *Nat Methods*. 2020;17(2):137-145. PMID: 31792435. doi: 10.1038/s41592-019-0654-x
<!-- REFERENCES:END -->

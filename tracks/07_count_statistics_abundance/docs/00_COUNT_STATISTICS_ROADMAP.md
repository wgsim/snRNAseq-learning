# 00. 단일세포 전사체 카운트 통계 및 차등 분석 로드맵 (Count Statistics Roadmap)

---

## 학습 목표 (Cognitive Goals)

1. **단일세포 카운트 통계학의 근본적 패러다임 전환 인식**: 단일세포 RNA 시퀀싱(scRNA-seq/snRNA-seq) 데이터에서 개별 세포의 수($n_{\text{cells}}$)와 독립적인 생물학적 개체 수($N_{\text{donors}}$)의 차이를 통계학적으로 엄밀하게 구별하고, 세포 수준 $i.i.d.$ 가정의 붕괴가 초래하는 위계적 데이터 구조의 수학적 본질을 규명한다.
2. **문헌 속 단일세포 재현성 위기(Reproducibility Crisis) 진단**: 림프구, 상피세포, 종양 미세환경 연구에서 널리 오용되어 온 세포 수준 단순 비모수 검정(Wilcoxon rank-sum test)과 포아송 기반 검정이 어떻게 수천 개의 가짜 차등 발현 유전자(False Positive DEGs)를 양산하고 통계적 검정력을 왜곡했는지 비판적으로 해체한다.
3. **대표적 임상 암 코호트 벤치마크 확립**: Kfoury et al. (Cancer Cell 2021)의 전립선암 골전이 환자 9명(각 환자 3개 분획) 및 양성 골수 대조 7명(GSE143791) 코호트와 Chen et al. (Nat Cell Biol 2021)의 원발 전립선 종양 13개 아틀라스(GSE141445)를 벤치마크 데이터셋으로 정의하고, 복잡한 임상 공변량 제어 전략을 수립한다.
4. **선행 트랙(Tracks 00~06)과의 상호 운용 체계 정립**: Track 04의 유전체 서브클론(Genomic Subclones) 특이적 발현 비교 및 Track 06의 연속적 분화 궤적(Continuous Trajectory/Pseudotime) 상에서의 동적 유전자 검정에서 카운트 기반 통계 모델이 결합되는 메커니즘을 설계한다.
5. **프로덕션 등급 카운트 분석 전산 환경 구축**: `scanpy`, `decoupler`, `statsmodels`, `edgeR`, `DESeq2`, `glmmTMB`, `nebula`로 이어지는 파이썬 및 R 하이브리드 통계 생태계를 완벽하게 구성한다.

---

## 1. 단일세포 카운트 통계학의 패러다임 전환: $n_{\text{cells}}$ vs $N_{\text{donors}}$

### 1.1 대용량 시퀀싱의 착시와 표본 크기(Sample Size)의 착각

현대 단일세포 전사체 기술(10x Chromium, Parse Biosciences, Smart-seq3 등)의 비약적 발전으로 단일 실험에서 $10^4 \sim 10^6$개의 세포를 손쉽게 프로파일링할 수 있게 되었다. 그러나 이 대규모 데이터 행렬 $X \in \mathbb{N}_0^{G \times n_{\text{cells}}}$은 통계학적으로 심각한 착시를 유발한다.

많은 연구자들이 $n = 50,000$개의 세포를 관측했다는 사실을 표본 크기 $N = 50,000$인 독립 표본으로 오인하고, 귀무가설 검정 시 표준 오차(Standard Error)가 $\sqrt{n}$에 반비례하여 0으로 수렴할 것이라 가정한다:

$$\text{SE}_{\text{naive}}(\bar{Y}) = \frac{\sigma}{\sqrt{n_{\text{cells}}}} \approx 0$$

그러나 동일한 생물학적 개체(환자, 마우스, 생검 조직)에서 유래한 세포들은 유전적 배경, 미세환경적 노출, 면역 상태, 조직 처리 배치(Batch) 효과를 공유한다. 즉, **세포들은 상호 독립적이지 않다(Non-independent)**.

```
[The Hierarchical Sampling Structure in Single-Cell Cohorts]

Cohort Level (True Biological Replicates):
       Donor 1 (Patient A)              Donor 2 (Patient B)              Donor 3 (Patient C)
      [ N_donors = 3 biological units of inference for treatment/disease effect ]
              │                                │                                │
              ├── Single-cell dissociation     ├── Single-cell dissociation     ├── Single-cell dissociation
              ▼                                ▼                                ▼
Cell Level (Observational Units / Technical Sub-samples):
   Cell 1, Cell 2, ..., Cell n_1    Cell 1, Cell 2, ..., Cell n_2    Cell 1, Cell 2, ..., Cell n_3
   (n_1 = 5,000 cells)              (n_2 = 4,200 cells)              (n_3 = 6,100 cells)
   [ Total cells = 15,300 cells: Correlated cluster measurements! ]
```

### 1.2 통계적 추론 단위의 괴리

- **실험 단위 (Experimental Unit)**: 치료제 투여, 질병 상태(Case vs Control), 유전자 녹아웃 등 연구자가 관심을 두는 처치(Treatment)가 독립적으로 할당되는 단위. 임상 및 전임상 연구에서 이는 **환자(Donor) 또는 동물 개체($N_{\text{donors}}$)**이다.
- **측정 단위 (Observational Unit)**: 전사체 발현량이 실제로 측정되는 최소 단위. 즉, **개별 단일세포($n_{\text{cells}}$)**이다.

실험 단위와 측정 단위가 일치하지 않는 위계적 구조에서 측정 단위(세포)를 독립 표본으로 취급하는 순간, 가설 검정의 유효 자유도(Effective Degrees of Freedom)는 과대평가되고, 검정 통계량은 비정상적으로 팽창하여 재앙적인 **가짜 반복(Pseudoreplication)**이 발생한다.

---

## 2. 단일세포 문헌의 재현성 위기 (The Reproducibility Crisis)

### 2.1 Squair et al. (Nature Communications, 2021)의 충격적 벤치마크

단일세포 문헌의 상당수가 Seurat의 `FindMarkers`나 Scanpy의 `sc.tl.rank_genes_groups`에서 기본 제공하는 **단일세포 수준 비모수 검정(Wilcoxon rank-sum test)** 또는 단순 포아송/음이항 회귀를 무비판적으로 사용해왔다.

2021년 Squair et al.은 matched bulk RNA-seq를 정답으로 둔 18개의 gold-standard 데이터셋에서 14개 통계 방법의 일치도(AUCC)를 비교하고, 생물학적 차이가 없는 무작위 그룹 비교에서 발견되는 DE 유전자 수를 확인하여 다음 정성적 결론을 보고하였다 (논문은 아래와 같은 Type I error 수치표를 제시하지 않으므로 수치는 인용하지 않는다):

```
[Squair et al. (2021) 정성적 결론 요약]
Test Method                      논문의 결론
──────────────────────────────────────────────────────────────────────────────────────────
세포 단위 방법(Wilcoxon 등)      생물학적 반복 간 변동을 무시 → 편향되고 위양성 발견에 취약
                                 (무작위 비교에서도 수백 개의 DE 유전자를 발견할 수 있음)
Pseudobulk (edgeR, DESeq2 등)    반복 간 변동을 반영 → 정답(bulk)과의 일치도가 가장 높음
GLMM 등 혼합모형                 약 500세포 이상에서 pseudobulk 성능에 수렴하나 계산 비용이 큼
──────────────────────────────────────────────────────────────────────────────────────────
```

단순 세포 수준 Wilcoxon 검정은 환자 간의 생물학적 변동(Donor-to-donor variability)을 무시하므로, 두 군(Case vs Control) 간에 실제 유전자 발현 차이가 없는 귀무가설($H_0$) 하에서도 **다수의 위양성 차등 발현 유전자(False DEGs)**를 보고할 수 있다(위양성의 정도는 반복 수·세포 수·발현량 등 조건에 의존한다).

### 2.2 $p$-value의 기만과 초극단적 유의성($p < 10^{-300}$)의 진실

많은 논문에서 $p = 10^{-50}$ 또는 $p = 0.0$(R 등에서 매우 작은 값이 0으로 표시되는 언더플로)과 같은 극단적인 $p$-value를 제시하며 발굴된 마커 유전자의 유의성을 과시한다. 

그러나 이 수치는 유전자의 실제 생물학적 효과 크기(Effect Size)를 증명하는 것이 아니라, 단지 **상관된 세포 수 $n$을 분모에 대입하여 표준 오차를 인위적으로 축소시킨 수치적 인공물(Artifact)**에 불과하다. 환자 3명 대 3명의 비교에서 순수한 생물학적 자유도는 단 $df = 3 + 3 - 2 = 4$에 불과하며, $df=4$에서 $p<10^{-50}$이 되려면 $|t|\approx5\times10^{12}$이 필요하므로(수학적으로는 가능하나) 현실적인 데이터로는 도달할 수 없는 수치이다.

---

## 3. 대표 임상 벤치마크 코호트

본 트랙에서는 단일세포 전사체 카운트 통계와 차등 분석의 엄밀한 검증을 위해, 실제 임상 난제와 복잡한 공변량을 내포한 두 가지 고난도 전립선암 단일세포/단일핵 코호트를 기준 데이터셋으로 채택한다.

```
[Benchmark Cohorts for Count Statistics & Differential Expression]

1. Kfoury et al. (Cancer Cell 2021) - GSE143791
   ┌─────────────────────────────────────────────────────────────┐
   │ Cohort: 9 PCa Bone-Metastasis Patients (Tumor/Involved/Distal│
   │         marrow fractions) + 7 Benign BM controls            │
   │ Platform: 10x Chromium 3' v2 scRNA-seq                      │
   │ Challenges: Severe donor-specific infiltration bias,        │
   │             extreme cell-type sparsity in osteoblasts,      │
   │             prior androgen deprivation therapy (ADT) history│
   └─────────────────────────────────────────────────────────────┘

2. Chen et al. (Nat Cell Biol 2021) - GSE141445
   ┌─────────────────────────────────────────────────────────────┐
   │ Cohort: 13 Primary Prostate Tumors (radical prostatectomy)  │
   │ Cells: 36,424 Single Cells (Epithelial & TME Atlas)         │
   │ Challenges: Patient age distribution,                       │
   │             Gleason score confounding, multi-batch designs, │
   │             luminal lineage plasticity along castration     │
   └─────────────────────────────────────────────────────────────┘
```

### 3.1 Kfoury et al. (Cancer Cell 2021) 골전이 전립선암 코호트 (GSE143791)
- **임상적 배경**: 전립선암 사망의 주원인인 골전이(Bone Metastasis) 미세환경 프로파일링.
- **통계적 도전 과제**:
  - 골 조직 특성상 환자 간 세포 수율의 불균형 가능성(예: 어떤 환자는 200개, 다른 환자는 5,000개라는 가정적 예시).
  - 환자별 이전 요법(화학요법, 호르몬 치료 등) 노출 이력의 이질성(예: Radium-223 노출 환자는 상세 분석에서 제외되었음).
  - 전이 미세환경 내 골아세포(Osteoblast), 파골세포(Osteoclast), 면역세포 간 상호작용 검정 시 환자 무작위 효과(Random Effects) 제어 필수.

### 3.2 Chen et al. (Nat Cell Biol 2021) 원발암 vs 진행성 전립선암 아틀라스 (GSE141445)
- **임상적 배경**: 근치적 전립선 절제술로 얻은 원발 종양 13개(36,424 세포)의 상피·미세환경 아틀라스(논문은 별도 11개 샘플로 CRPC 비교를 추가).
- **통계적 도전 과제**:
  - 환자 연령(Age), 병기(Stage), 글리슨 점수(Gleason Score) 등 연속형/순서형 임상 공변량 제어.
  - 관강 세포(Luminal Epithelia) 내부에서 치료 저항성에 따른 클론 확장 및 조성비 변화(Differential Abundance)와 유전자 발현 변화(Differential Expression)의 상호 분리.

---

## 4. 선행 트랙(Tracks 00~06)과의 유기적 통합 체계

본 트랙 07의 카운트 통계 모델은 독립적으로 존재하는 것이 아니라, 이전 트랙들에서 도출된 계산생물학적 객체들을 가설 검정의 단위로 수용한다.

```
[Integration Architecture with Tracks 00 ~ 06]

Track 04: Cancer Genomics & Clones
  │ (InferCNV / Copy Number Calling)
  ▼
Subclone Architecture: [Clone A (Diploid)] vs [Clone B (8q Gain, PTEN Del)]
  │
  ├─► Statistical Challenge: Cells within Clone B belong to specific patients!
  └─► Solution: Hierarchical GLMM testing Clone Effect while blocking Patient Effect!
       log(mu_ij) = beta_0 + beta_clone * Clone_ij + u_patient(i) + ...

Track 06: Manifold Trajectory & Kinetics
  │ (CellRank 2 / Waddington-OT / Latent Time)
  ▼
Continuous Differentiation Trajectory: Latent Time tau in [0, 1]
  │
  ├─► Statistical Challenge: Testing dynamic DEG along tau without pseudoreplication!
  └─► Solution: Generalized Additive Mixed Models (GAMM) or Pseudobulk windowing
       log(mu_ij) = f(tau_ij) + u_patient(i) + log(depth_ij)
```

1. **Track 04 (유전체 서브클론 가설 검정)**:
   - 종양 세포 내에서 $PTEN$ 결손 또는 $AR$ 증폭을 보유한 클론 B와 잔존 이배체 클론 A 간의 차등 발현을 비교할 때, 클론 A와 클론 B가 여러 환자에 걸쳐 분포하는 경우 환자 효과($u_i$)를 무작위 효과로 블록화(Block)하지 않으면 클론 특이적 신호가 환자 고유 발현 변이와 뒤섞인다.
2. **Track 06 (연속적 분화 궤적 상의 동역학 검정)**:
   - 전사체 유사시간(Pseudotime, $\tau$) 또는 잠재 시간(Latent Time)을 따라 점진적으로 유도/억제되는 드라이버 유전자를 규명할 때, 연속형 변수 $\tau$를 공변량으로 취급하면서 개체 간 기저 발현 변이를 제어하는 일반화 가법 혼합 모델(GAMM) 및 위계적 스플라인 모델이 필수적으로 요구된다.

---

## 5. 트랙 07 전체 아키텍처 및 학습 로드맵

본 트랙은 총 10개의 장(00~09)으로 구성되며, 단일세포 카운트 데이터의 확률적 기저부터 위계적 모델링, 궤적 DE, 차등 풍부도·조성 분석, 트러블슈팅까지 단계별로 심화된다:

```
[Track 07 Curriculum Architecture]

Part 1: Probability Foundations & Pseudobulk Framework (Docs 00 ~ 04)
 ├── 00. Count Statistics Roadmap (Current)
 ├── 01. Probability Models: Poisson to Negative Binomial (Gamma-Poisson Derivation)
 ├── 02. The Pseudoreplication Pitfall & Type I Error Explosion (ICC & N_eff)
 ├── 03. Pseudobulk Framework: edgeR & DESeq2 (Empirical Bayes Shrinkage & QLF)
 └── 04. Hierarchical GLMM: Random Effects & Cell-level Covariates (glmmTMB & Laplace)

Part 2: Cell-level Models, Trajectory DE & Differential Abundance (Docs 05 ~ 08)
 ├── 05. Fast Cell-level Mixed Models: NEBULA
 ├── 06. Trajectory Differential Expression (tradeSeq/GAM)
 ├── 07. Neighborhood-based Differential Abundance: Milo (k-NN Graph Overlapping Neighborhoods)
 └── 08. Compositional Data Analysis (CoDa): propeller/speckle, scCODA

Part 3: Best Practices (Doc 09)
 └── 09. Troubleshooting & Best Practices (IHW, QC, 루브릭)
```

---

## 6. 프로덕션 분석 환경 구축 (Environment Setup)

단일세포 카운트 통계와 혼합 모델 분석을 수행하기 위해 Python 생태계와 R Bioconductor 생태계를 긴밀하게 연동하는 전용 가상환경을 구성한다.

### 6.1 Conda/Mamba 환경 생성 및 필수 패키지 설치

```bash
# 1. 전용 환경 생성 (Python 3.11 및 R 4.3 기반; decoupler 2.x가 Python>=3.11 요구)
conda create -n sn_stats python=3.11 r-base=4.3.2 -c conda-forge -c bioconda -y
conda activate sn_stats

# 2. Python 핵심 단일세포 및 통계 패키지 설치
pip install \
    scanpy==1.10.1 \
    anndata==0.10.7 \
    "decoupler>=2.0" \
    statsmodels==0.14.2 \
    scipy==1.13.0 \
    scikit-misc==0.3.1 \
    patsy==0.5.6 \
    rpy2==3.5.15 \
    sccoda   # scCODA(Python); mgcv는 R 기본 권장 패키지로 포함

# 3. R 핵심 통계 및 Bioconductor 패키지 설치 (R 콘솔 실행)
R -e '
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")

BiocManager::install(c(
    "edgeR",
    "DESeq2",
    "limma",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "miloR",
    "scater",
    "sccomp",        # sccomp는 CRAN이 아니라 Bioconductor 패키지
    "speckle", "tradeSeq", "slingshot", "IHW", "RUVSeq", "sva"
), ask = FALSE)

install.packages(c(
    "glmmTMB",
    "lme4",
    "nebula",
    "tidyverse",
    "statmod",
    "car"
), repos="https://cloud.r-project.org")
'
```

### 6.2 Python-R 상호운용성 검증 스크립트

```python
"""
00_test_stats_environment.py
환경 구축 및 Python-R(edgeR, glmmTMB) 연동 무결성 검증
"""
import sys
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm

print(f"[Python System Info] Python Version: {sys.version.split()[0]}")
print(f"[Python Packages] Scanpy: {sc.__version__}, AnnData: {ad.__version__}")

# rpy2를 통한 R 핵심 패키지 로드 테스트
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    
    base = importr('base')
    edger = importr('edgeR')
    deseq2 = importr('DESeq2')
    glmmtmb = importr('glmmTMB')
    nebula = importr('nebula')
    
    r_ver = ro.r('R.version.string')[0]
    print(f"[R Integration] Successfully linked to: {r_ver}")
    print("[R Packages] edgeR, DESeq2, glmmTMB, nebula imported successfully!")
    print("\n>>> All Count Statistics Environment Components are READY!")
except Exception as e:
    print(f"[Warning] R integration test encountered an issue: {e}")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 표본 크기 착시와 자유도(Degrees of Freedom) 붕괴
어느 암 전사체 연구팀이 대장암 환자 4명(정상 점막 4개체 vs 종양 조직 4개체, 총 $N=8$)으로부터 총 80,000개의 단일 상피세포를 시퀀싱하였다. 연구팀은 80,000개 세포를 데이터 테이블의 개별 행으로 설정하고 두 군(정상 40,000개 vs 종양 40,000개) 간에 독립 표본 $t$-검정을 수행하여 3,500개의 유전자가 $p < 10^{-100}$의 유의수준으로 차등 발현된다고 보고하였다.
1. 이 분석에서 귀무가설 검정에 적용된 표본 크기($n$)와 실제 통계적 추론의 기초가 되는 생물학적 자유도($df$)의 차이를 명시하시오.
2. 세포 간 상관성이 존재할 때, 단순 독립 표본 검정이 산출하는 표준 오차($\text{SE}$)와 $t$-통계량($t = \frac{\bar{X}_1 - \bar{X}_2}{\text{SE}}$)이 겪게 되는 수학적 왜곡 기전을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 표본 크기와 실제 자유도의 괴리**:
- **연구팀이 착각한 표본 크기 및 자유도**: 연구팀은 개별 세포를 독립된 관측치로 간주하여 $n_1 = 40,000, n_2 = 40,000$으로 설정하였다. 이때 적용된 겉보기 자유도는 $df_{\text{naive}} = 40,000 + 40,000 - 2 = 79,998$이다.
- **실제 생물학적 표본 크기 및 유효 자유도**: 치료 및 질병 상태가 부여된 독립적인 생물학적 실험 단위는 환자(Donor) 개체 수이다 ($N_1 = 4, N_2 = 4$). 따라서 유전자가 질병 군과 정상 군 간에 진정으로 차이가 있는지를 검정하기 위한 엄밀한 통계적 자유도는 $df_{\text{true}} = N_1 + N_2 - 2 = 4 + 4 - 2 = 6$에 불과하다.

**2. 표준 오차 과소평가 및 $t$-통계량 폭발 기전**:
- 동일 환자 유래 세포들 간의 공분산(Covariance)은 양수($\text{Cov}(X_{ij}, X_{ik}) = \sigma_{\text{donor}}^2 > 0$)이다.
- 클러스터링된 데이터에서 평균의 참 분산은 다음과 같다:
  $$\text{Var}(\bar{X}) = \frac{\sigma_{\text{cell}}^2}{N \cdot n} + \frac{\sigma_{\text{donor}}^2}{N}$$
- 그러나 단순 독립 표본 공식을 사용할 경우 $\sigma_{\text{donor}}^2$ 항이 완전히 누락되고 분모가 $N \cdot n = 80,000$으로 처리되어 표준 오차($\text{SE}$)가 이론적 참값보다 최대 수십 배에서 수백 배까지 비정상적으로 작게 계산된다:
  $$\text{SE}_{\text{naive}} \ll \text{SE}_{\text{true}}$$
- 그 결과, $t$-통계량 분모의 극단적 축소로 인해 $t = \frac{\bar{X}_1 - \bar{X}_2}{\text{SE}_{\text{naive}}}$가 비정상적으로 팽창하여 수치적으로 $t > 30 \implies p < 10^{-100}$이라는 허구의 초유의성(Spurious Hyper-significance)이 도출된다.

</details>

---

### 문제 2: Kfoury et al. 골전이 데이터셋의 비대칭적 세포 분포와 가짜 상관
Kfoury et al. 골전이 코호트(환자 9명)와 유사한 구조를 가정한 가상 시나리오이다. 골아세포(Osteoblast) 군집을 분석하는 과정에서, 특정 환자 P07 한 명에게서만 전체 골아세포 4,000개 중 3,200개(80%)가 포획되었고 나머지 8명의 환자에게서는 평균 100개씩의 골아세포만 포획되었다.
1. 이 상태에서 환자 P07이 포함된 치료 저항군과 비치료군 간에 세포 수준 Wilcoxon rank-sum 검정을 수행할 경우, P07 환자 개인의 고유 유전적 변이(Private Genetic Variations)나 배치 효과가 전체 저항군 마커 발굴에 미치는 치명적 교란을 논하시오.
2. 슈도벌크(Pseudobulk) 집계 프레임워크가 이러한 단일 환자 독점 왜곡을 어떻게 방어하는지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 단일 환자 독점에 따른 개인 변이의 마커 둔갑(Private Variation Hijacking)**:
- 세포 수준 Wilcoxon 검정은 모든 세포의 발현 순위(Rank)를 동등하게 취급한다.
- 저항군 골아세포의 80%가 환자 P07 단 한 명의 세포이므로, 저항군 세포 집단의 발현량 분포는 사실상 '저항군의 보편적 생물학'이 아니라 '환자 P07의 개인적 발현 특성'을 그대로 반영하게 된다.
- 예를 들어, 환자 P07이 우연히 보유한 생식세포계 변이(Germline Variant), P07 개인의 체질, 또는 P07 샘플 채취 당일의 조직 분리 효소 농도 편차로 인해 발현이 높은 무관한 유전자 $GeneX$가 존재한다면, 세포 수준 검정에서는 3,200개의 P07 세포가 모두 높은 발현 순위를 차지하여 $GeneX$가 $p < 10^{-200}$을 기록하는 최상위 "골전이 치료 저항성 핵심 바이오마커"로 잘못 식별된다.

**2. 슈도벌크 집계의 방어 기전**:
- 슈도벌크 프레임워크는 각 환자별로 세포의 카운트를 하나의 프로파일로 합산(Summation)하므로, 환자 P07은 아무리 많은 세포(3,200개)를 가졌더라도 저항군 내에서 **단 하나의 생물학적 샘플($N=1$)**로 환원된다.
- edgeR이나 DESeq2와 같은 슈도벌크 통계 엔진은 환자 간 분산(Inter-donor variance)을 모델링하므로, $GeneX$가 P07 환자 한 명에서만 극단적으로 높고 나머지 저항군 환자들($N_2, N_3, \dots$)에서 재현되지 않는다면 높은 개체 간 분산($\phi$)으로 인해 통계적 유의성이 완전히 기각된다.
- 즉, 슈도벌크는 특정 환자의 세포 과다 추출로 인한 위양성 편향을 크게 완화한다(단, 세포 수가 적은 환자의 pseudobulk는 잡음이 커 최소 셀/카운트 필터와 가중이 필요하며, 환자 수가 적으면 환자 특이 변이의 검정력·재현성 문제는 남는다).

</details>

---

### 문제 3: Track 04 서브클론과 Track 06 연속 궤적 상의 슈도벌크 한계
환자 10명의 암 조직 단일세포 데이터에서, 각 환자의 종양 내부에는 정상 이배체 클론, 서브클론 A($TP53$ 변이), 서브클론 B($AR$ 증폭)가 서로 다른 비율로 공존하고 있다. 또한 세포들은 줄기세포 유사 상태에서 분화된 관강 세포 상태로 연속적인 분화 궤적 잠재 시간($\tau \in [0, 1]$)을 따라 이동 중이다.
1. 이 데이터셋에서 서브클론 간 유전자 발현 차이를 규명하거나 잠재 시간 $\tau$에 따른 동적 발현을 분석할 때, 환자 단위로 세포를 완전히 합산해 버리는 단순 슈도벌크(Vanilla Pseudobulk) 방식이 직면하는 두 가지 근본적 정보 손실을 서술하시오.
2. 이러한 한계를 극복하기 위해 제4장에서 다룰 위계적 일반화 선형 혼합 모델(GLMM)이 도입되어야 하는 당위성을 수학적 선형 예측자($\eta_{ij}$) 형태로 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 단순 슈도벌크의 두 가지 근본적 정보 손실**:
- **세포 내 이질성 및 서브클론 해상도 붕괴**: 환자별로 세포를 통째로 합산하면, 동일 환자 내부에서 공존하는 서브클론 A, 서브클론 B, 정상 클론의 카운트가 믹서기에 갈리듯 하나로 혼합(Bulk Averaging)되어 서브클론 고유의 발현 특이성이 희석된다.
- **연속적 공변량의 소실(Loss of Continuous Covariates)**: 각 세포가 가진 연속적 분화 잠재 시간 $\tau \in [0, 1]$이나 세포 주기 점수($G_2/M$ Score), 미토콘드리아 분율과 같은 단일세포 수준 연속 변수는 단일 스칼라 요약값(평균 등)으로 압축되면서 세포별 동역학적 궤적 정보가 파괴된다.

**2. 위계적 GLMM 도입의 당위성 및 선형 예측자 정식화**:
- GLMM은 세포 수준의 연속형/범주형 공변량을 그대로 유지하면서도, 환자 간의 상관성을 무작위 효과(Random Effects)로 완벽히 흡수한다.
- 환자 $i$의 세포 $j$에 대한 로그 연결 선형 예측자 $\eta_{ij} = \log(\mu_{ij})$:
  $$\log(\mu_{ij}) = \beta_0 + \beta_{\text{clone}} \cdot \text{Clone}_{ij} + f(\tau_{ij}) + \mathbf{X}_{ij}^{\text{cell}} \boldsymbol{\gamma} + u_i + \log(\text{depth}_{ij})$$
  - $\text{Clone}_{ij}$: 세포 $j$가 속한 유전체 서브클론 지시 변수 (Track 04 연계).
  - $f(\tau_{ij})$: 잠재 시간 $\tau$에 따른 비선형 동적 발현 스플라인 함수 (Track 06 연계).
  - $\mathbf{X}_{ij}^{\text{cell}}$: 세포 주기 및 미토콘드리아 분율 등 세포 수준 교란 변수.
  - $u_i \sim \mathcal{N}(0, \sigma_{\text{donor}}^2)$: 환자 $i$의 고유 기저 발현 변이를 흡수하는 무작위 절편(Random Intercept).
  - $\log(\text{depth}_{ij})$: 세포별 시퀀싱 깊이에 대한 오프셋(Offset).
- 이를 통해 세포 수준의 미세한 공변량을 온전히 살리면서도 환자 수준의 가짜 반복을 수학적으로 엄밀하게 차단한다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nat Commun*. 2021;12(1):5692. PMID: 34584091. doi: 10.1038/s41467-021-25960-2
2. Kfoury Y, Baryawno N, Severe N, et al. Human prostate cancer bone metastases have an actionable immunosuppressive microenvironment. *Cancer Cell*. 2021;39(11):1464-1478.e8. PMID: 34719426. doi: 10.1016/j.ccell.2021.09.005
3. Chen S, Zhu G, Yang Y, et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. *Nat Cell Biol*. 2021;23(1):87-98. PMID: 33420488. doi: 10.1038/s41556-020-00613-6
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- anndata 0.13.4 (github.com/scverse/anndata). 원 논문: Virshup I, Rybakov S, Theis FJ, et al. anndata: Access and store annotated data matrices. *JOSS*. 2024;9(101):4371. doi: 10.21105/joss.04371
- DESeq2 (bioconductor.org/packages/DESeq2). 원 논문: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15(12):550. PMID: 25516281. doi: 10.1186/s13059-014-0550-8
- edgeR (bioconductor.org/packages/edgeR). 원 논문: Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26(1):139-40. PMID: 19910308. doi: 10.1093/bioinformatics/btp616
- glmmTMB (github.com/glmmTMB/glmmTMB). 원 논문: Brooks M, Kristensen K, Benthem K, et al. glmmTMB Balances Speed and Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling. *The R Journal*. 2017;9(2):378. doi: 10.32614/rj-2017-066
- NEBULA (github.com/lhe17/nebula). 원 논문: He L, Davila-Velderrain J, Sumida TS, et al. NEBULA is a fast negative binomial mixed model for differential or co-expression analysis of large-scale multi-subject single-cell data. *Commun Biol*. 2021;4(1):629. PMID: 34040149. doi: 10.1038/s42003-021-02146-6
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- rpy2 3.6.9 (github.com/rpy2/rpy2).
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- statsmodels 0.15.0 (github.com/statsmodels/statsmodels).
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
<!-- REFERENCES:END -->

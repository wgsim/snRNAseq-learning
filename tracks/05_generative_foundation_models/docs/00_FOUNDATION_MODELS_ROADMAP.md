# 00. 단일세포 전사체 생성형 파운데이션 모델 및 딥러닝 로드맵 (Generative Foundation Models Roadmap)

---

## 학습 목표 (Cognitive Goals)

1. **고전적 선형 차원 축소 및 휴리스틱 배치 보정의 수학적 한계 규명**: 주성분 분석(PCA), Harmony, Seurat CCA 등 전통적 선형·휴리스틱 접근법이 단일세포 전사체 카운트의 비가우시안성(Non-Gaussianity)과 비선형 생물학적 다양성을 왜곡하는 통계적 결함을 이해하고, 심층 생성 모델(Deep Generative Modeling)로의 패러다임 전환 필연성을 도출할 수 있다.
2. **다기관·초대형 단일세포 아틀라스 통합의 핵심 수리적 난제 분석**: 시퀀싱 플랫폼 차이(10x Chromium v2/v3, 5', snRNA vs scRNA), 환자 간 이질성(Donor Heterogeneity), 시퀀싱 깊이(Library Size) 차이가 유발하는 비선형 배치 효과(Batch Effects)와 연속적 세포 분화 궤적(Continuous Plasticity) 보존 간의 최적화 트레이드오프를 수리생물학적으로 정식화할 수 있다.
3. **전립선암(mCRPC/t-NEPC) 및 전신 아틀라스 벤치마크 코호트 카탈로그 체계화**: 원발성 전립선암(TCGA-PRAD), 전이성 거세저항성 전립선암(SU2C/PCF), 치료유도 신경내분비 분화(t-NEPC) 단일세포 코호트와 Human Cell Atlas(HCA) Pan-tissue 아틀라스, Replogle et al. 전장 유전체 Perturb-seq 섭동 벤치마크의 유전체학적 특성을 파악한다.
4. **엔터프라이즈급 심층학습 단일세포 연산 환경 구축**: `python 3.10+`(최신 scvi-tools는 더 높은 Python 요구), `torch` (CUDA/MPS 하드웨어 가속), `scvi-tools`, `transformers`(Geneformer용), `scanpy` 생태계를 연계한 재현 가능한 딥러닝 파이프라인 아키텍처를 설계하고 환경을 검증할 수 있다.

---

## 1. 단일세포 전사체 분석의 패러다임 전환: 휴리스틱 선형 파이프라인에서 생성형 파운데이션 모델로

단일세포 RNA 시퀀싱(scRNA-seq) 및 단일핵 RNA 시퀀싱(snRNA-seq)은 종양 미세환경(Tumor Microenvironment, TME)의 고도 이질성을 해명하는 표준 기술로 확립되었다. 그러나 지난 10년간 단일세포 생물정보학을 지배해 온 고전적 파이프라인은 데이터의 물리적 측정 원리와 본질적으로 충돌하는 휴리스틱 가정들에 의존해 왔다.

```
[Evolution of Single-Cell Computational Paradigms]

Phase 1: Linear & Heuristic Era (2014 - 2018)
  Raw Counts ──> Library-size Norm ──> log(1+x) ──> Highly Variable Genes ──> PCA ──> Harmony / Seurat CCA ──> UMAP / Leiden
  * 한계: 카운트 통계학 파괴, 인위적 분산 왜곡, 배치 효과의 선형적 독립 가정

Phase 2: Deep Generative Modeling Era (2018 - 2022)
  Raw Counts ──> Probabilistic Count Likelihood (NB/ZINB) ──> Variational Autoencoder (scVI) ──> Latent Space z ~ N(0, I)
  * 혁신: 원시 카운트 직접 모델링, 비선형 배치 분리, 불확실성 정량화, 준지도 학습(scANVI)

Phase 3: Large Foundation Model & Continual Learning Era (2022 - Present)
  Massive Multi-Atlas Pretraining ──> Gene-tokenized Transformers (scGPT, Geneformer) ──> Architectural Surgery (scArches)
  * 도약: Zero-shot 세포형 분류, 유전자 발현 섭동(in silico Perturbation) 예측, 지속 학습(Lifelong Learning)
```

### 1.1 고전적 휴리스틱 선형 파이프라인의 통계적 결함

전통적인 단일세포 분석 파이프라인(`Seurat v3/v4`, 표준 `Scanpy`)은 다음과 같은 일련의 휴리스틱 변환을 거친다:

1.  **라이브러리 크기 정규화(Library Size Normalization)**:
    세포 $c$의 유전자 $g$에 대한 측정 카운트 $x_{cg}$에 대해, 총 시퀀싱 리드 수 $N_c = \sum_{g} x_{cg}$를 목표 상수 $S = 10,000$으로 나눈다:
    $$\tilde{x}_{cg} = \frac{x_{cg}}{N_c} \times 10,000$$
2.  **휴리스틱 로그 변환(Log-Transformation)**:
    $$y_{cg} = \log(\tilde{x}_{cg} + 1)$$
3.  **선형 주성분 분석(PCA)**:
    $$Z = Y W, \quad W \in \mathbb{R}^{G \times d}$$

이러한 고전적 흐름은 세 가지 치명적 수리통계적 결함을 내포한다:

*   **가산 가우시안 노이즈(Additive Gaussian Noise) 가정의 붕괴**:
    PCA는 데이터의 오차가 정규분포(Gaussian distribution)를 따르며 등분산성(Homoskedasticity)을 만족한다는 수학적 전제하에 최대 분산 직교 기저(Orthogonal Basis)를 찾는다. 그러나 scRNA-seq 카운트는 0 이상의 불연속 정수(Discrete Non-negative Integers)이며, 평균 발현량이 낮을수록 기술적 영점(Technical Zeros)의 비율이 기하급수적으로 증가하는 이분산성(Heteroskedasticity)을 띤다.
*   **$\log(1+x)$ 변환에 의한 인위적 분산 왜곡(Artifactual Variance Distortion)**:
    Townes et al. (*Genome Biology*, 2019)은 $\log(1+x)$ 정규화가 카운트 자료의 분산 구조를 왜곡해 차원 축소에서 인위적 변동을 만들 수 있음을 지적하였다(GLM-PCA 제안). 반면 Ahlmann-Eltze & Huber (*Nature Methods*, 2023)는 델타법·피어슨 잔차·잠재 상태·인자 분석 기반 변환을 비교한 결과 의사카운트를 둔 로그 변환 후 PCA가 정교한 대안과 동등하거나 더 나은 성능을 보임을 보고하였다. 따라서 로그 변환의 한계(의사카운트·라이브러리 크기 의존성)는 실제로 있으나 항상 큰 왜곡을 낳는다고 단정할 수 없으며, 확률적 생성 모델의 필요성은 이 맥락에서 평가해야 한다.
*   **Harmony 및 선형 배치 보정의 과잉 보정(Over-correction) 위험**:
    Harmony(Korsunsky et al., *Nature Methods*, 2019)는 PCA 임베딩 공간에서 최대 다양성 군집화(Maximum Diversity Clustering)를 통해 각 클러스터 내 배치 비율이 전체 표본 비율과 동일해지도록 선형 보정 벡터를 계산한다. 그러나 환자 코호트 간에 질환 진행 단계가 상이하여 특정 코호트에만 존재하는 고유한 세포 상태(예: 전이성 거세저항성 전립선암 환자에만 존재하는 신경내분비 종양세포)가 존재할 경우, Harmony는 이를 환자 간 배치 효과로 간주하여 강제로 정상 관강 상피세포와 혼합해 버리는 치명적 '생물학적 신호 소거(Biological Signal Erasure)'를 초래한다.

### 1.2 심층 생성 모델(scVI, scANVI, scArches)의 수학적 패러다임

심층 생성 모델은 원시 정수 카운트(Raw Count)를 어떠한 로그 변환이나 휴리스틱 정규화 없이 확률론적 우도 함수(Probabilistic Likelihood Function)로 직접 모델링한다:

$$x_{cg} \sim \text{NegativeBinomial}\left( \text{mean} = l_c \cdot \rho_{cg}(z_c, s_c), \ \text{dispersion} = \theta_{sg} \right)$$

*   $l_c$: 세포 $c$의 시퀀싱 깊이(라이브러리 크기 잠재 변수)
*   $z_c \in \mathbb{R}^d$: 배치 효과가 줄어들도록 조건화된 저차원 세포 상태 벡터 ($d \approx 10 \sim 30$; 완전한 배치 불변성은 모델이 보장하지 않으며 데이터에 따라 잔존한다)
*   $s_c$: 기술적 배치 식별자 (원-핫 벡터 또는 범주형 인덱스)
*   $\rho_{cg}(z_c, s_c)$: 심층 신경망(Decoder Network)을 통해 복원된 유전자 $g$의 상대적 생물학적 발현 비율 ($\sum_g \rho_{cg} = 1$)
*   $\theta_{sg}$: 유전자별/배치별 과산포 파라미터

신경망의 인코더 $q_\phi(z_c \mid x_c, s_c)$는 고차원 비선형 매핑을 통해 복잡한 기술적 노이즈를 흡수하고, 세포의 본질적 유전자 발현 조절 프로그램을 $z_c$ 공간에 보존한다.

### 1.3 대규모 트랜스포머 기반 파운데이션 모델(scGPT, Geneformer)의 지평

자연어 처리(NLP) 분야의 트랜스포머(Transformer) 아키텍처가 전사체학으로 확장되면서, 유전자를 단어(Token), 세포를 문서(Document)로 간주하는 '유전자 언어 모델(Gene Language Models)'이 대두되었다:

*   **scGPT (Cui et al., *Nature Methods*, 2024)**: 유전자 토큰 임베딩, 발현량 연속값 임베딩, 배치 토큰 임베딩을 결합한 Generative Pretrained Transformer 아키텍처를 통해 3,300만 개 이상의 단일세포 전사체를 사전학습. 유전자 간 동시 발현 주의집중(Self-Attention) 맵을 통해 유전자 조절 네트워크(GRN)를 자가 학습함.
*   **Geneformer (Theodoris et al., *Nature*, 2023)**: 3,000만 개 전사체를 기반으로 각 세포 내 유전자를 코퍼스 전체의 유전자별 발현으로 스케일한 값 기준으로 순위화(Rank-value encoding; 하우스키핑 유전자는 강등)해 입력으로 사용하며, 희귀 심근병증 및 유전자 결실 섭동의 영향을 zero-shot/few-shot으로 탐색(검증 범위는 논문 참조).

---

## 2. 대규모 단일세포 아틀라스 통합의 핵심 난제 (Core Integration Challenges)

단일 연구실 단위의 소규모 실험(1~2만 개 세포)을 넘어 전신 규모(Pan-tissue, 수백만~수천만 개 세포) 아틀라스를 구축할 때 계산생물학자는 세 가지 근본적 장벽에 직면한다.

```
┌────────────────────────────────────────────────────────────────────────┐
│            대규모 단일세포 아틀라스 통합의 3대 핵심 난제              │
└────────────────────────────────────────────────────────────────────────┘
                                    │
         ┌──────────────────────────┼──────────────────────────┐
         ▼                          ▼                          ▼
 1. 복합 배치 효과          2. 희귀 상태 및 분화        3. 계산 복잡도 및
    (Complex Batch Effects)    연속 궤적 보존             확장성 (Scalability)
  - 10x v2 vs v3 vs 5'       - t-NEPC (1% 미만)         - $O(N^2)$ 메모리 폭발
  - snRNA vs scRNA           - 중간 전이 상태 왜곡       - 전역 재학습 불가
  - 다기관 검체 동결 프로토콜 - Mode Collapse 방지       - 분산 미니배치 최적화
```

### 2.1 다기관·다중 플랫폼 배치 효과의 비선형성

단일세포 데이터의 기술적 변동은 단순한 평균값 이동(Location Shift)이 아니다:
1.  **화학 반응 및 프로토콜 차이**: 10x Chromium Single Cell 3' v2 키트는 v3 키트에 비해 역전사 효율 및 캡처 효율이 낮아 동일 유전자에 대해서도 평균 UMI 카운트와 드롭아웃 확률이 상이하다.
2.  **scRNA-seq vs snRNA-seq**: 고형암 조직의 장기 동결 검체는 단일핵(single-nucleus) 분리를 적용해야 한다. snRNA-seq은 세포질 RNA(Cytoplasmic mRNA, e.g., 미토콘드리아 유전자, 리보솜 단백질)가 대폭 결손되고 미성숙 전구전사체(Intronic pre-mRNA) 비율이 30~50%에 달해 전사체 조성 자체가 근본적으로 왜곡된다.
3.  **환자 간 유전적 배경(Genetic Background)**: 환자별 고유한 HLA 유전자형, 체세포 돌연변이 프로파일, 면역 반응 강도에 따라 동일 세포형 내에서도 기저 발현량이 비선형적으로 변동한다.

### 2.2 생물학적 희귀 아형 및 연속적 분화 궤적 보존

대규모 아틀라스 통합 알고리즘의 주요 실패 양상은 **과도한 혼합(Over-mixing)**이다. 전립선암 거세저항성(mCRPC) 환자에서 안드로겐 수용체 차단 치료(Enzalutamide) 후 발생하는 **치료유도 신경내분비 전립선암(t-NEPC)** 세포는 전체 표본의 1% 미만으로 존재할 수 있다:
*   선형 정렬 알고리즘은 다수 집단(수십만 개의 관강 상피세포 및 기질세포)의 분산을 최소화하는 방향으로 손실 함수를 최적화하므로, 희귀한 t-NEPC 집단을 정상 세포 클러스터 내부로 강제 정렬하여 소멸시키는 치명적 오류를 범한다.
*   따라서 생성 모델은 세포형의 사전 확률(Prior Probability)과 잠재 공간의 연속적 밀도 분포를 온전히 보존해야 한다.

### 2.3 계산 복잡도와 지속 학습(Lifelong Learning)의 필요성

전통적인 그래프 기반 배치 보정(MNN, BBKNN)은 모든 세포 쌍 간의 $k$-최근접 이웃(k-Nearest Neighbors) 그래프를 구축해야 하므로 시간 및 공간 복잡도가 $O(N^2)$ 또는 $O(N \log N)$에 달한다:
*   세포 수 $N = 10^7$(1천만 개)에 이르면 그래프 구축에 수 테라바이트(TB)의 RAM과 수일의 연산 시간이 소요된다.
*   더욱이 새로운 환자 코호트(예: 신규 임상시험 환자 10명, 5만 개 세포)가 추가될 때마다 기존 1천만 개 아틀라스 전체를 처음부터 다시 재학습(Full Retraining)하는 것은 막대한 탄소 배출과 계산 비용으로 인해 불가능하다.
*   이 문제를 해결하기 위해 기존 아틀라스를 완전히 동결(Freeze)하고 신규 쿼리 코호트의 기술적 파라미터만 국소적으로 수술(Architectural Surgery)하는 **scArches** 패러다임이 필수적으로 요구된다.

---

## 3. 대표 벤치마크 아틀라스 및 코호트 카탈로그

생성형 단일세포 모델의 검증 및 임상적 적용을 위해 확립된 글로벌 표준 벤치마크 데이터셋 카탈로그를 정리한다.

### 3.1 전립선암(Prostate Cancer) 단일세포 및 분자 유전체 코호트

| 코호트 / 연구 | 출판 연도 및 저널 | 표본 규모 및 시퀀싱 기술 | 분자생물학적 및 계산생물학적 핵심 의의 |
| :--- | :--- | :--- | :--- |
| **TCGA-PRAD** | *Cell* (2015) 163:1011–1025 | 원발성 전립선 선암종 333례 (Bulk WES, RNA-seq, SNP 6.0) | 원발암 분자 아형 기준선 확립; *ERG*, *ETV1* 융합 및 8p 결손, 8q 획득 등 기저 유전체 랜드마크 제공 |
| **SU2C/PCF mCRPC** | *Cell* (2015) 161:1215–1228<br>*PNAS* (2019) 116:11428–11436 | 전이성 거세저항성 전립선암 생검 코호트: Cell 2015 150례, PNAS 2019 429례 (WES + Bulk RNA-seq) | 치료 저항성 핵심 유전체 드라이버 확립 (Cell 2015 초록 기준 *AR*·ETS·*TP53*·*PTEN* 이상이 각각 40~60% 범위로 빈번) |
| **Chen et al.** | *Nature Cell Biology* (2021) 23:87–98 (doi:10.1038/s41556-020-00613-6) | 원발성 전립선 종양 13례, 36,424개 단일세포 (scRNA-seq) | 질병 공격성과 연관된 상피 세포·종양미세환경 전사체 프로그램 규명 (`Dong et al., Science 2020` 인용은 존재하지 않아 교체) |
| **He et al.** | *Nature Medicine* (2021) 27:426–433 | 전이성 거세저항성 전립선암(mCRPC) 환자 14명, 다부위 전이 검체 (snRNA-seq) | 뼈, 간, 림프절 전이 병변의 단일핵 아틀라스; 종양 내 클론 간 발현 다양성 및 화학치료 저항 서브클론 규명 |
| **Labrecque et al.** | *J Clin Invest* (2019) 129:4492–4505 | mCRPC 전이 검체 및 PDX 라인 (전장 RNA-seq, GSEA, IHC) | *AR*/신경내분비 유전자 발현에 기반한 5개 mCRPC 표현형(ARPC, ARLPC, AMPC, DNPC, SCNPC) 분류 및 REST 활성 역할 |

### 3.2 전신 규모 인간 세포 아틀라스 (Human Cell Atlas, HCA)

1.  **Tabula Sapiens (The Tabula Sapiens Consortium, *Science*, 2022)**:
    *   24개 조직에서 추출한 약 50만 개의 단일세포 전사체 아틀라스(다수 조직이 같은 기증자에서 유래).
    *   같은 기증자에서 여러 조직을 확보해 유전적 배경 영향을 줄인 채 조직 특이적 대식세포, 혈관내피세포, 섬유아세포의 장기 간 전사체 차이를 체계적으로 분석할 수 있는 절대적 레퍼런스.
2.  **Human Lung Cell Atlas (HLCA, Sikkema et al., *Nature Medicine*, 2023)**:
    *   486명의 개체에서 수집된 240만 개 세포를 scVI 및 scArches를 통해 통합한 최대 규모 호흡기 질환 표준 아틀라스.
    *   다기관 코호트의 체계적 통합 파이프라인의 글로벌 표준으로 인정받음.

### 3.3 대규모 기능 유전체학 섭동 벤치마크: Perturb-seq

*   **Replogle et al. (*Cell*, 2022)**:
    *   **발현 유전자 전체를 표적**(Genome-scale)으로 한 CRISPR 간섭(CRISPRi) 단일세포 전사체 데이터셋(총 250만 개 이상 단일세포; 전장 규모는 K562, RPE1은 필수 유전자 위주 스크린으로 알려져 있으며 유전자 수치는 직접 확인하지 않음).
    *   각 유전자 녹다운(Knockdown)에 따른 비선형 전사 반응을 단일세포 해상도로 측정하여, scGPT, Geneformer, scVI 등 생성 모델이 유전자 간 인과적 조절 네트워크(Causal Gene Networks)를 얼마나 정확히 예측하는지 평가하는 골드 스탠다드 벤치마크.

---

## 4. 계산 환경 구축 및 딥러닝 소프트웨어 스택

단일세포 생성 모델 및 파운데이션 모델 실습을 위한 딥러닝 환경은 고성능 GPU 텐서 연산 및 분산 미니배치 데이터로더를 지원해야 한다.

```
[Deep Generative Single-Cell Software Stack]

  User Application / Pipeline
    └── scvi-tools (v1.1+) / Hugging Face Transformers (v4.35+)
          │
  Core Single-Cell Data Structure
    └── AnnData (v0.10+) / Scanpy (v1.9+)
          │
  Deep Learning Compute Engine
    └── PyTorch (v2.1+) (torch.distributions, Autograd, CUDA / MPS)
          │
  Low-Level Acceleration Library
    └── NVIDIA CUDA 12.x / cuDNN / Apple Metal Performance Shaders (MPS)
```

### 4.1 핵심 의존성 패키지 명세

*   `python >= 3.10`: 단 최신 `scvi-tools`(1.5.1 기준 PyPI)는 `python >= 3.12`를 요구하므로, 3.10/3.11 환경에서는 구버전이 설치된다.
*   `torch >= 2.1.0`: 신경망 구축 및 `torch.distributions.NegativeBinomial` 기반 확률 밀도 계산.
*   `scvi-tools >= 1.1.0`: `scVI`, `scANVI`, `scArches` 모델링을 위한 표준 프로덕션 라이브러리.
*   `scanpy >= 1.9.6`: 고품질 전처리, HVG 선정, UMAP 시각화.
*   `anndata >= 0.10.3`: 메모리 매핑(`backed="r"`) 및 희소 행렬(CSR Matrix) 처리.
*   `transformers >= 4.35.0`: `Geneformer`의 HF `BertForMaskedLM` 백본 연동용(`scGPT` 패키지는 `scanpy`·`torch`·`torchtext`에 의존하며 `transformers`를 쓰지 않는다).

### 4.2 가상환경 구축 및 하드웨어 가속 검증 스크립트

다음 스크립트는 딥러닝 단일세포 연산 가상환경을 구성하고, GPU(NVIDIA CUDA 또는 Apple Silicon MPS) 가속 및 핵심 라이브러리 가용성을 자동으로 검증한다.

```python
"""
00_verify_environment.py
단일세포 생성형 파운데이션 모델 전용 하드웨어 가속 및 의존성 검증 스크립트.
"""

import sys
import platform
import importlib

def check_single_cell_dl_stack():
    print(f"=== Single-Cell Deep Learning Environment Diagnostic ===")
    print(f"Python Version : {sys.version.split()[0]}")
    print(f"Platform       : {platform.platform()}")
    print("-" * 55)

    # 1. 핵심 패키지 임포트 검증
    required_packages = [
        ("torch", "PyTorch Core Engine"),
        ("scvi", "scvi-tools (Generative Models)"),
        ("scanpy", "Scanpy (Single-Cell Toolkit)"),
        ("anndata", "AnnData Data Structure"),
        ("transformers", "Hugging Face Transformers (Geneformer용; scGPT만 쓰면 선택)"),
        ("scipy", "SciPy Scientific Computing"),
        ("numpy", "NumPy Matrix Operations")
    ]

    all_passed = True
    for pkg_name, description in required_packages:
        try:
            mod = importlib.import_module(pkg_name)
            ver = getattr(mod, "__version__", "unknown")
            print(f" [PASS] {pkg_name:<15} v{ver:<10} - {description}")
        except ImportError as e:
            print(f" [FAIL] {pkg_name:<15} NOT FOUND   - {description} ({e})")
            all_passed = False

    print("-" * 55)

    # 2. PyTorch 하드웨어 텐서 가속 검증
    import torch
    print(f"PyTorch Version: {torch.__version__}")
    
    device_found = False
    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_name(0)
        gpu_count = torch.cuda.device_count()
        vram_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        print(f" [ACCELERATOR] NVIDIA CUDA Active:")
        print(f"   - GPU Device : {gpu_name} (Total: {gpu_count} GPU(s))")
        print(f"   - Total VRAM : {vram_gb:.2f} GB")
        print(f"   - CUDA Arch  : {torch.version.cuda}")
        device_found = True
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        print(f" [ACCELERATOR] Apple Silicon MPS (Metal Performance Shaders) Active")
        device_found = True
        device = torch.device("mps")
    else:
        print(f" [WARNING] No Hardware Acceleration Detected. Fallback to CPU.")
        device = torch.device("cpu")

    # 3. 텐서 연산 벤치마크 테스트 (Negative Binomial Log-Likelihood)
    x = torch.randint(low=0, high=100, size=(1000, 2000), device=device).float()
    mu = torch.rand(1000, 2000, device=device) * 50.0 + 1e-4
    theta = torch.ones(2000, device=device) * 2.0
    
    # NB log prob 수치 테스트
    from torch.distributions import NegativeBinomial
    # torch NegativeBinomial의 평균은 total_count * p / (1 - p) 이므로 평균을 mu로 만들려면 p = mu / (theta + mu)
    probs = mu / (theta + mu)
    nb_dist = NegativeBinomial(total_count=theta, probs=probs)
    log_probs = nb_dist.log_prob(x)
    
    print(f" [TEST] Negative Binomial Likelihood Tensor Computation:")
    print(f"   - Input Shape : {tuple(x.shape)}")
    print(f"   - Device Used : {x.device}")
    print(f"   - Mean Log-Likelihood: {log_probs.mean().item():.4f}")
    print("-" * 55)
    
    if all_passed and device_found:
        print(">>> 시스템 준비 완료: 모든 생성형 파운데이션 모델 학습 가능.")
    else:
        print(">>> 환경 구성 확인 필요: 미흡한 의존성 또는 CPU 전용 모드.")

if __name__ == "__main__":
    check_single_cell_dl_stack()
```

---

## 5. 단원별 아키텍처 및 커리큘럼 로드맵 (Curriculum Map)

본 트랙은 단일세포 전사체의 수리통계적 확률 분포 모델링부터 최신 아틀라스 지속 학습 아키텍처까지 10개 챕터(00~09)로 구성된다. 아래 도식은 00~04(확률 모델·scVI 계열)만 그리고, 05~09(트랜스포머 계열·섭동·트러블슈팅)는 표에 정리한다.

```
[Track 05 Foundation Models Curriculum Architecture]

┌────────────────────────────────────────────────────────────────────────┐
│ 00_FOUNDATION_MODELS_ROADMAP.md (본 단원)                               │
│ - 선형 정렬의 붕괴, 심층 생성 모델 패러다임, HCA/Perturb-seq 벤치마크      │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│ 01_PROBABILISTIC_COUNT_MODELING.md                                     │
│ - 포아송 등분산성 실패, Gamma-Poisson (음이항, NB) 과산포 모델           │
│ - ZINB 제로 팽창 게이팅 vs UMI 프로토콜 논쟁, 라이브러리 크기 승법 분해   │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│ 02_SCVI_VARIATIONAL_AUTOENCODER.md                                     │
│ - 변분 오토인코더(VAE) 구조, 인코더/디코더 네트워크 매핑               │
│ - 재파라미터화 트릭(Reparameterization Trick) 및 ELBO 손실 함수 완전 유도│
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│ 03_SCANVI_SEMISUPERVISED_ANNOTATION.md                                 │
│ - 대형 아틀라스 주석/비주석 비대칭, 결합 확률 모델 p(x, y, z | s)        │
│ - 레이블 교차 엔트로피 + 비레이블 주변화(Marginalization) 손실, t-NEPC 전이│
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│ 04_SCARCHES_ARCHITECTURAL_SURGERY.md                                   │
│ - 파괴적 망각 극복, 동결 레퍼런스(Frozen Reference) + 쿼리 전용 어댑터    │
│ - 연합 학습 기반 프라이버시 보존 아틀라스 확장 및 신규 임상 환자 매핑     │
└────────────────────────────────────────────────────────────────────────┘
```

| 챕터 번호 및 파일명 | 핵심 다루는 주제 | 주요 수식 및 신경망 아키텍처 | 주요 실습 도구 |
| :--- | :--- | :--- | :--- |
| **00. 로드맵** | 패러다임 전환, 대형 코호트 벤치마크, 딥러닝 환경 구축 | $Z = Y W$ 선형 축소의 한계, Perturb-seq | `scvi-tools`, `PyTorch` |
| **01. 확률적 카운트 모델링** | 포아송 등분산 붕괴, 음이항(NB) 과산포, ZINB, 승법 분해 | $\text{Var}(X) = \mu + \frac{\mu^2}{\theta}$, $\mu_{cg} = l_c \rho_{cg}$ | `torch.distributions` |
| **02. scVI VAE** | 인코더-디코더, 재파라미터화 트릭, ELBO 수리 유도, 배치 보정 | $\mathcal{L}_{\text{ELBO}} = \mathbb{E}[\log p(x \mid z)] - D_{\text{KL}}$, $z = \mu + \sigma \odot \epsilon$ | `scvi.model.SCVI` |
| **03. scANVI 준지도 주석** | 주석 비대칭 극복, 결합 생성 모델, 비레이블 주변화 손실 | $\mathcal{L}_{\text{unlabeled}} = \sum_y q(y \mid x) \mathcal{L}(x, y) + \mathcal{H}(q)$, t-NEPC 전이 | `scvi.model.SCANVI` |
| **04. scArches 아키텍처 수술** | 동결 레퍼런스 아틀라스, 어댑터 레이어, 파괴적 망각 방지, 연합 학습 | $h' = \sigma(W_{\text{ref}}x + W_{\text{query}}s + b)$, $z_{\text{query}} \to z_{\text{ref}}$ | `scvi.model.SCVI.load_query_data` |
| **05. 트랜스포머 토큰화** | 단일세포-언어 모델 상동성, 유전자 토큰화, 셀프 어텐션 | $\text{softmax}(QK^\top/\sqrt{d})V$ | `torch.nn` (직접 구현) |
| **06. scGPT 임베딩** | 값 구간화(binning), 다중 태스크 사전학습, 세포 임베딩 추출 | 값 임베딩 + 유전자 임베딩 | `torch.nn` (직접 구현 예시) |
| **07. Geneformer** | 순위 기반 입력, 어텐션 기반 유전자 네트워크 분석 | rank 정렬 시퀀스 + 절대 위치 임베딩 | `torch.nn` (직접 구현 예시) |
| **08. 인 실리코 섭동** | 유전자 섭동 예측, 가상 약물 스크리닝 | 섭동 전후 임베딩 코사인 변위 | `torch` (직접 구현 예시) |
| **09. 트러블슈팅·벤치마크** | GPU 메모리, 수치 안정성, scIB 벤치마킹, 캡스톤 루브릭 | scIB 종합 점수(생물 보존/배치 제거 가중) | `torch.cuda.amp`, `sklearn` |

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: Harmony 선형 보정과 scVI 딥러닝 VAE의 수리적 메커니즘 차이 및 과잉 보정 위험
다기관에서 수집된 20만 개의 거세저항성 전립선암(mCRPC) 단일세포 전사체 데이터를 통합하고자 한다. 환자군 A(10명)는 표준 호르몬 요법에 반응하는 관강 선암종(Luminal adenocarcinoma) 환자들로만 구성되어 있고, 환자군 B(3명)는 공격적인 치료유도 신경내분비 전립선암(t-NEPC, *SYP+*, *CHGA+*)으로 분화한 말기 환자들이다. 
1. Harmony 알고리즘을 적용했을 때 t-NEPC 세포 집단에 발생할 수 있는 '과잉 보정(Over-correction)'의 수학적 원인을 최대 다양성 군집화(Maximum Diversity Clustering)의 가정에 기반하여 설명하시오.
2. scVI가 조건부 변분 오토인코더(Conditional VAE) 구조를 통해 이 문제를 어떻게 해결하며, 배치 효과와 실제 생물학적 표현형 차이를 분리해내는 메커니즘을 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. Harmony의 과잉 보정 수리적 메커니즘**:
*   Harmony는 PCA 공간상에서 세포들을 클러스터링할 때, 각 클러스터 내의 배치(환자) 분포가 전체 데이터셋의 배치 사전 비율(Prior Batch Distribution)과 동일해지도록 정보 엔트로피를 최대화하는 **최대 다양성 군집화(Maximum Diversity Clustering)**를 수행한다:
    $$R_{ki} \propto \exp\left(-\frac{\|z_i - \mu_k\|^2}{\sigma}\right) \times \left( \frac{P(\text{batch}_i)}{P(\text{batch}_i \mid k)} \right)^\theta$$
*   여기서 환자군 B에만 독점적으로 존재하는 t-NEPC 세포들은 환자군 A의 세포를 전혀 포함하지 않는다. Harmony의 페널티 함수는 이를 "특정 배치(환자군 B)에 과도하게 편향된 클러스터"로 간주하여 심각한 페널티를 부여한다.
*   그 결과, Harmony는 순수 생물학적 악성 변이(t-NEPC 특이 발현)를 환자군 간 기술적 배치 편차로 오인하여, 환자군 B의 t-NEPC 세포들을 환자군 A의 정상/선암종 관강세포 클러스터 내부로 강제로 끌어당겨 병합(Over-correction)시킨다. 이는 임상적으로 치명적인 희귀 아형 신호 소멸을 초래한다.

**2. scVI의 분리 메커니즘**:
*   scVI는 배치 식별자 $s$를 디코더 네트워크의 조건부 입력(Conditional input)으로 명시적으로 주입한다:
    $$x_{cg} \sim \text{NegativeBinomial}\left( l_c \cdot f_\theta(z_c, s_c)_g, \ \theta_{sg} \right)$$
*   디코더 $f_\theta$는 $s_c$를 통해 환자군 간의 전역적 측정 효율, 시퀀싱 깊이, 배치 특이적 기술적 절편(Intercept)을 별도로 흡수한다.
*   반면 인코더를 통해 생성되는 잠재 변수 $z_c$는 $p(z) = \mathcal{N}(0, I)$로 정규화(KL)되며, 배치 효과의 일부를 디코더의 $s_c$ 항이 흡수할 수 있어 $z_c$에 배치 정보가 줄어드는 경향이 있다. 다만 ELBO에는 $z \perp s$를 강제하는 항이 없으므로 독립성이 **보장되지는 않는다**(scvi-tools 기본값은 `encode_covariates=False`로 인코더가 $s$를 입력받지 않으며 디코더에만 깊이 주입한다).
*   위 시나리오처럼 t-NEPC가 환자군 B에만 존재하는 **완전 교란**에서는 `batch_key=환자`로 조건화할 때 모델이 t-NEPC 고유 신호의 일부를 배치 항으로 흡수할 수 있으므로 scVI도 과잉 보정에서 자유롭지 않다. 따라서 보정 후 SYP/CHGA/ASCL1 같은 마커 발현으로 군집을 검증하고, 가능하면 배치(환자)와 생물학이 교차하는 샘플을 확보해야 한다.

</details>

---

### 문제 2: mCRPC 코호트에서 t-NEPC 희귀 전이 아형 추적 시 휴리스틱 로그 정규화의 왜곡
단일세포 RNA 시퀀싱에서 일반적으로 사용되는 라이브러리 정규화 및 로그 변환 $\log(1 + \frac{x_{cg}}{N_c} \times 10^4)$은 전립선암 상피세포의 전이 과정 추적 시 심각한 통계적 왜곡을 유발한다.
1. 전립선 선암종 세포(Luminal adenocarcinoma)와 t-NEPC 세포 간에 전사 활성도(Total RNA Content, $N_c$)가 현저히 다를 때, 단순 라이브러리 크기 비례 정규화가 저발현 전사인자(Transcription Factors, e.g., *ASCL1*, *FOXA2*)와 고발현 구조 유전자(e.g., *ACTB*, *MALAT1*)의 상대적 비율에 미치는 왜곡 효과를 논하시오.
2. 확률적 생성 모델(scVI)이 라이브러리 크기 $l_c$를 상수가 아닌 잠재 변수(Latent Variable)로 모델링함으로써 이를 극복하는 방식을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 휴리스틱 로그 정규화의 왜곡 효과**:
*   **세포 크기 및 전사 출력(Transcriptional Output)의 비대칭**: 빠르게 증식하는 미분화 t-NEPC 세포는 크기가 작고 핵-세포질 비율(N:C ratio)이 높아 세포당 총 mRNA 분자 수($N_c$)가 분화된 대형 관강 상피세포의 30~50% 수준에 불과할 수 있다.
*   단순 비율 정규화 $\frac{x_{cg}}{N_c} \times 10^4$는 총 RNA 양이 적은 세포에서 모든 유전자의 측정값을 인위적으로 2~3배 팽창시킨다.
*   **$\log(1+x)$의 비선형 압축 왜곡**: 고발현 구조 유전자(e.g., *ACTB*)는 로그의 평평한 구간(Plateau)에 도달하여 실제 생물학적 차이가 축소되는 반면, *ASCL1*이나 *FOXA2*와 같은 마스터 전사인자는 1~2 카운트 수준의 미세한 발현이 0에서 급격히 증가하는 구간에 위치하여 노이즈가 과대 증폭된다. 이는 거짓 양성(False Positive) 분화 경로를 생성하거나 세포 주기 효과를 분화 궤적으로 오인하게 만든다.

**2. scVI의 잠재 변수 모델링 접근법**:
*   scVI(원 논문 모델)는 관측 카운트의 총합을 단순 고정 상수로 나누지 않고, 세포별 라이브러리 크기 $l_c$를 관측 데이터로부터 추론되는 로그 정규 잠재 변수로 설정한다 (참고: 이는 Lopez et al. 2018 원 논문의 모델 서술이다. scvi-tools의 기본값 `use_observed_lib_size=True`에서는 $l_c$를 관측 총합 $N_c$로 고정하며 잠재 변수·`l_encoder`·LogNormal 사전을 쓰지 않는다. `False`로 지정하면 이 서술의 잠재 변수 모델이 된다.):
    $$l_c \sim \text{LogNormal}(\mu_l, \sigma_l^2)$$
*   디코더는 순수 생물학적 상대 발현 비율 $\rho_{cg} \in (0, 1)$ ($\sum_g \rho_{cg} = 1$)를 $z_c$로부터 먼저 복원한 후, 이를 라이브러리 스케일 $l_c$와 곱하여 음이항 평균 $\mu_{cg} = l_c \cdot \rho_{cg}$를 산출한다.
*   이로 인해 총 RNA 함량 차이에 의한 기술적 스케일 효과가 생물학적 유전자 발현 패턴($\rho_{cg}$)과 엄격히 분리(Deconvolution)되므로, t-NEPC 특이 전사인자의 진정한 농도 변화를 왜곡 없이 추정할 수 있다.

</details>

---

### 문제 3: Perturb-seq 대규모 유전체 섭동 데이터 분석에서 파운데이션 모델의 GRN 임베딩 원리
Replogle et al. (*Cell*, 2022)의 전장 유전체 CRISPR 섭동 데이터셋(Perturb-seq)은 250만 개의 세포에 걸쳐 수천 개 유전자의 녹다운 전사체 결과를 담고 있다.
1. 고전적 차원 축소 기법(PCA)이 유전자 섭동에 의한 비선형적 하류 표적(Downstream Targets) 연쇄 반응을 포착하지 못하는 이유를 선형 중첩(Linear Superposition)의 한계 관점에서 설명하시오.
2. 트랜스포머 기반 파운데이션 모델(scGPT, Geneformer)의 자기 주의집중(Self-Attention) 메커니즘이 세포 내 유전자 조절 네트워크(Gene Regulatory Network, GRN)의 유향(Directed) 및 상위성(Epistasis) 상호작용을 학습할 수 있는 메커니즘을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 고전적 PCA의 선형 중첩 한계**:
*   PCA는 임의의 복합 섭동 상태 $Y$를 주성분 고유벡터(Eigenvectors)들의 선형 결합 $Y = \sum_i c_i v_i$으로 표현한다.
*   그러나 세포 내 신호전달 및 유전자 조절은 스위치형 온-오프(Sigmoidal/Hill kinetics), 피드백 루프, 그리고 두 유전자가 동시에 억제되었을 때 예상치 못한 제3의 경로가 활성화되는 상위성(Epistasis/Synthetic Lethality)과 같은 극단적 비선형성을 보인다.
*   선형 모델은 유전자 A 억제 효과와 유전자 B 억제 효과의 기하학적 합($v_A + v_B$)만을 예측할 수 있으므로, 임계점을 넘었을 때 발생하는 세포 사멸이나 신경내분비성 계통 전이(Lineage Plasticity)와 같은 복잡한 상태 전이를 포착하기 어렵다. 다만 이는 이론적 한계이며, 실제 섭동 예측 벤치마크에서는 단순 선형·가산 기준선이 딥러닝 파운데이션 모델과 동등하거나 더 낫다는 보고(Ahlmann-Eltze 외, 2025)가 있으므로 기준선 비교로 확인해야 한다.

**2. 트랜스포머 Self-Attention의 GRN 학습 메커니즘**:
*   트랜스포머 모델은 각 유전자 토큰 $g_i$에 대해 Query ($Q$), Key ($K$), Value ($V$) 벡터를 생성하고, 유전자 쌍 간의 주의집중 가중치를 동적으로 계산한다:
    $$\text{Attention}(Q, K, V) = \text{softmax}\left( \frac{Q K^T}{\sqrt{d_k}} \right) V$$
*   여기서 주의집중 행렬 $A_{ij} = \text{softmax}(q_i k_j^T / \sqrt{d_k})$는 유전자 $j$의 발현 상태가 유전자 $i$의 임베딩 표현을 업데이트하는 데 미치는 영향력을 나타내며, 이를 유전자 조절 네트워크의 방향성 에지(Directed Edge) 후보로 해석하려는 시도가 있으나, 어텐션 가중치는 인과·조절 관계 자체가 아니다(Geneformer 논문도 어텐션이 네트워크 위계를 인코딩한다고만 주장). 따라서 어텐션 기반 후보는 실험적·외부 데이터 검증이 필요한 가설이다.
*   다층 트랜스포머 레이어(Multi-layer Transformers)를 통과하면서, 상위 레이어는 단순 1:1 조절을 넘어 전사인자 결합 $\to$ 크로마틴 리모델러 활성화 $\to$ 하류 표적 유전자 발현으로 이어지는 고차 연쇄 조절(Higher-order Epistatic Interactions)을 표현할 수 있는 용량을 갖는다(실제로 학습·재현하는지는 별도 검증 대상이다).
*   이론적으로는 관측되지 않은 신규 유전자 조합 섭동(Unseen in silico Perturbations)에 대한 반응을 예측하는 데 쓰일 수 있으나, Ahlmann-Eltze 외(Nature Methods, 2025)는 다섯 파운데이션 모델과 두 DL 모델이 단일·이중 섭동 예측에서 단순 기준선을 이기지 못했다고 보고하였다. 따라서 섭동 예측은 기준선(평균·가산 모델) 대비 평가가 필수이다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Replogle JM, Saunders RA, Pogson AN, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*. 2022;185(14):2559-2575.e28. PMID: 35688146. doi: 10.1016/j.cell.2022.05.013
2. Townes FW, Hicks SC, Aryee MJ, et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model. *Genome Biol*. 2019;20(1):295. PMID: 31870412. doi: 10.1186/s13059-019-1861-6
3. Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. *Nat Methods*. 2023;20(5):665-672. PMID: 37037999. doi: 10.1038/s41592-023-01814-1
4. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods*. 2019;16(12):1289-1296. PMID: 31740819. doi: 10.1038/s41592-019-0619-0
5. Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nat Methods*. 2024;21(8):1470-1480. PMID: 38409223. doi: 10.1038/s41592-024-02201-0
6. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618(7965):616-624. PMID: 37258680. doi: 10.1038/s41586-023-06139-9
7. Tabula Sapiens Consortium*, Jones RC, Karkanias J, et al. The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. *Science*. 2022;376(6594):eabl4896. PMID: 35549404. doi: 10.1126/science.abl4896
8. Sikkema L, Ramírez-Suástegui C, Strobl DC, et al. An integrated cell atlas of the lung in health and disease. *Nat Med*. 2023;29(6):1563-1577. PMID: 37291214. doi: 10.1038/s41591-023-02327-2
9. Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods*. 2025;22(8):1657-1661. PMID: 40759747. doi: 10.1038/s41592-025-02772-6
10. Cancer Genome Atlas Research Network. The Molecular Taxonomy of Primary Prostate Cancer. *Cell*. 2015;163(4):1011-25. PMID: 26544944. doi: 10.1016/j.cell.2015.10.025
11. Robinson D, Van Allen EM, Wu YM, et al. Integrative clinical genomics of advanced prostate cancer. *Cell*. 2015;161(5):1215-1228. PMID: 26000489. doi: 10.1016/j.cell.2015.05.001
12. Abida W, Cyrta J, Heller G, et al. Genomic correlates of clinical outcome in advanced prostate cancer. *Proc Natl Acad Sci U S A*. 2019;116(23):11428-11436. PMID: 31061129. doi: 10.1073/pnas.1902651116
13. Chen S, Zhu G, Yang Y, et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. *Nat Cell Biol*. 2021;23(1):87-98. PMID: 33420488. doi: 10.1038/s41556-020-00613-6
14. He MX, Cuoco MS, Crowdis J, et al. Transcriptional mediators of treatment resistance in lethal prostate cancer. *Nat Med*. 2021;27(3):426-433. PMID: 33664492. doi: 10.1038/s41591-021-01244-6
15. Labrecque MP, Coleman IM, Brown LG, et al. Molecular profiling stratifies diverse phenotypes of treatment-refractory metastatic castration-resistant prostate cancer. *J Clin Invest*. 2019;129(10):4492-4505. PMID: 31361600. doi: 10.1172/JCI128212
16. Lopez R, Regier J, Cole MB, et al. Deep generative modeling for single-cell transcriptomics. *Nat Methods*. 2018;15(12):1053-1058. PMID: 30504886. doi: 10.1038/s41592-018-0229-2
### 검수 근거 문헌
1. Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods*. 2022;19(1):41-50. PMID: 34949812. doi: 10.1038/s41592-021-01336-8
2. Kedzierska KZ, Crawford L, Amini AP, et al. Zero-shot evaluation reveals limitations of single-cell foundation models. *Genome Biol*. 2025;26(1):101. PMID: 40251685. doi: 10.1186/s13059-025-03574-x
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- PyTorch 2.14.1 (github.com/pytorch/pytorch).
<!-- REFERENCES:END -->

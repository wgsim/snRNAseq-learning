# 03. scANVI 준지도 학습 및 희귀 세포형 고정밀 전이 (scANVI Semi-Supervised Annotation)

---

## 학습 목표 (Cognitive Goals)

1. **단일세포 아틀라스의 주석화/비주석화 비대칭(Annotated/Unannotated Asymmetry) 문제 규명**: 수십만~수백만 개의 세포로 구성된 초대형 임상 전사체 코호트에서 수동 레이블링의 계산적·임상적 한계를 이해하고, 준지도 학습(Semi-Supervised Learning)의 수학적 필요성을 도출할 수 있다.
2. **scANVI 결합 생성 모델(Joint Generative Model)의 확률적 수리 구조 분석**: 연속적 세포 잠재 상태 $z$와 이산적 세포형 범주 $y$를 동시에 포함하는 결합 확률 분포 $p_\theta(x, y, z \mid s)$ 및 분류기 신경망 $q_\phi(y \mid x, s)$의 아키텍처를 수리적으로 정식화할 수 있다.
3. **레이블 우도와 비레이블 주변화 손실(Marginalization Loss)의 완전한 수리 유도**: 레이블 보유 세포에 대한 분류 교차 엔트로피와 비레이블 세포에 대한 전세포형 가중 합산 손실 $\sum_{y} q(y \mid x) \mathcal{L}(x, y)$ 및 샤논 엔트로피 $\mathcal{H}(q)$의 수리적 역할을 증명할 수 있다.
4. **거세저항성 전립선암(mCRPC) 코호트에서 희귀 t-NEPC 전이 파이프라인 구현**: `scvi-tools` 기반 scVI 비지도 사전학습 $\to$ scANVI 준지도 전이 학습 파이프라인을 구축하고, 1% 미만 수준(예시)으로 존재할 수 있는 희귀 신경내분비암(*SYP+*, *CHGA+*) 서브클론을 후보로 식별하고, 정밀도·재현율(PR)과 혼동행렬로 검증하는 방법을 이해한다.

---

## 1. 아틀라스 스케일의 주석화 비대칭 (Annotated/Unannotated Asymmetry)

단일세포 생물정보학에서 세포형 주석화(Cell-Type Annotation)는 세포 생물학자의 수작업 마커 유전자 검증(Manual Curation)에 크게 의존해 왔다. 그러나 글로벌 컨소시엄(HCA, TCGA, SU2C)이 수백만 개의 세포 데이터를 쏟아내면서 심각한 데이터 비대칭이 발생하였다.

```
[The Single-Cell Atlas Annotation Bottleneck]

Reference Atlas Cohorts (Curated)                 Unannotated Clinical Cohorts (Raw Queries)
- Well-characterized donors                      - Multicenter mCRPC / Metastatic biopsies
- Pathologist-verified labels                    - Uncurated, high batch heterogeneity
- High quality, high coverage                    - Novel treatment-resistant subclones
- N_ref ≈ 50,000 cells                           - N_query ≈ 500,000+ cells
                     │                                         │
                     └────────────────────┬────────────────────┘
                                          ▼
                      How to transfer labels with high fidelity
                      without destroying latent geometry?
```

### 1.1 휴리스틱 레이블 전이(k-NN / Seurat Transfer Anchors)의 실패 원인

고전적 레이블 전이 파이프라인(e.g., PCA 공간상 $k$-최근접 이웃 보간, Seurat v3/v4 CCA 앵커링)은 다음과 같은 치명적 결함을 드러낸다:

1.  **배치 효과에 의한 앵커 왜곡(Anchor Distortion)**: 참조 데이터(Reference)와 임상 쿼리(Query) 데이터 사이에 플랫폼 차이(10x v2 vs v3)나 장기 동결 조직 특유의 배경 노이즈가 존재할 경우, $k$-NN은 유전자 발현 패턴이 아닌 기술적 유사도에 의해 앵커를 연결하여 대량의 오분류를 초래한다.
2.  **모드 붕괴 및 다수 클래스 편향(Majority Class Bias)**: 관강 상피세포(Luminal cells)와 같이 참조 아틀라스에 90% 이상 존재하는 다수 세포형 쪽으로 이웃 투표(Majority voting)가 편향되어, 전체의 1% 미만으로 존재하는 희귀 클론(e.g., t-NEPC, Foxp3+ Treg)이 다수 세포형으로 흡수되어 증발한다.
3.  **불확실성 정량화(Uncertainty Quantification)의 부재**: 기존 결정론적(Deterministic) 분류기는 쿼리 세포가 참조 아틀라스에 전혀 존재하지 않는 '신규 세포형(Novel Cell State)'인 경우에도 억지로 기존 클래스 중 하나로 강제 할당한다.

---

## 2. scANVI 결합 생성 모델 (Joint Generative Model)

scANVI(single-cell ANnotation using Variational Inference; Xu et al., *Molecular Systems Biology*, 2021)는 킹마 등의 M2 준지도 변분 모델(Kingma et al., *NeurIPS*, 2014)을 단일세포 카운트 확률 모델에 결합한 생성 모델이다.

```
[scANVI Probabilistic Directed Graphical Model]

        y_c ~ Cat(π) (Discrete Cell Type, K classes)
          │
          ▼
        z_c | y_c ~ N(u_y, diag(σ_y^2))  (Cell State)
          │
          ├────────────────────────┐
          ▼                        ▼
  Neural Decoder f_θ         Classifier q_φ(y | x, s)
  ρ_cg = Softmax(f_θ(z, s))  p_cat = Softmax(MLP_φ(x, s))
          │
          ▼
  x_cg ~ NB(l_c · ρ_cg, θ_sg)
```

### 2.1 계층적 생성 확률 모델의 정식화

scANVI는 각 세포 $c$에 대해 세 가지 확률변수의 결합 분포(Joint Distribution)를 정의한다:
*   $y_c \in \{1, 2, \dots, K\}$: 이산 세포형 레이블 (Categorical Label)
*   $z_c \in \mathbb{R}^d$: 연속적 잠재 발현 상태 (Continuous Latent State)
*   $x_c \in \mathbb{N}_0^G$: 관측 카운트 벡터 (Observed Count Matrix)
*   $s_c$: 기술적 배치 식별자

$$p_\theta(x, y, z \mid s) = p(y) \, p(z \mid y) \, p_\theta(x \mid z, s)$$

1.  **세포형 사전분포(Categorical Prior)**:
    $$p(y) = \text{Categorical}(\pi), \quad \pi_k = \frac{1}{K} \text{ (혹은 참조 데이터의 클래스 비율)}$$
2.  **세포형 조건부 잠재 사전분포(Conditional Latent Prior)**:
    각 세포형 $k$는 잠재 공간에서 고유한 평균 $u_k$와 분산 $\sigma_k^2$을 갖는 가우시안 중심(Gaussian Centroid)을 형성한다:
    $$p(z \mid y = k) = \mathcal{N}(u_k, \, \text{diag}(\sigma_k^2))$$
3.  **카운트 생성 우도(Count Likelihood)**:
    $$p_\theta(x \mid z, s) = \prod_{g=1}^G \text{NegativeBinomial}\left( x_{cg}; \ \mu_{cg} = l_c \cdot \rho_{cg}(z_c, s_c), \ \theta_{s_c, g} \right)$$

### 2.2 변분 근사 사후분포 (Variational Approximation)

참 사후분포 $p_\theta(z, y \mid x, s)$를 근사하기 위해, 다음과 같은 분해를 만족하는 인코더 신경망 쌍을 구축한다:

$$q_\phi(z, y \mid x, s) = q_\phi(z \mid x, y, s) \, q_\phi(y \mid x, s)$$

1.  **분류기 네트워크(Classifier Network)**:
    $$q_\phi(y \mid x, s) = \text{Categorical}\left( \text{Softmax}\left( \text{MLP}_{\text{cls}}(x, s) \right) \right)$$
    고차원 유전자 발현 패턴 $x$와 배치 $s$로부터 세포형 $y$의 소속 확률을 출력한다.
    (참고: 이는 Kingma M2 계열의 단순화 모형이다. scvi-tools 구현(Xu et al. 2021)은 분류기가 $x$가 아닌 잠재 $z_1$을 입력으로 받고($q(y\mid z_1)$, SVG03의 $q_\psi(y\mid z)$와 동일), 두 개의 잠재층($z_1$, $z_2$)과 균일한 $y$ 사전분포를 사용한다.)
2.  **잠재 인코더 네트워크(Latent Encoder Network)**:
    $$q_\phi(z \mid x, y, s) = \mathcal{N}\left( \mu_\phi(x, y, s), \, \text{diag}(\sigma_\phi^2(x, y, s)) \right)$$

---

## 3. 손실 함수의 수학적 유도: 레이블 우도와 전세포형 주변화 (Marginalization)

scANVI의 전체 데이터셋 $\mathcal{D}$는 레이블이 존재하는 부분집합 $\mathcal{D}_L = \{(x_i, y_i, s_i)\}$과 레이블이 결측된 부분집합 $\mathcal{D}_U = \{(x_j, s_j)\}$로 엄격히 분할된다.

### 3.1 레이블 세포(Labeled Cells)에 대한 손실 함수

세포형 $y$가 사전에 관측된 경우, 조건부 ELBO에 분류기 교차 엔트로피(Cross-Entropy)를 정규화 항으로 결합한다:

$$\log p_\theta(x, y \mid s) \ge \mathbb{E}_{q_\phi(z \mid x, y, s)}\left[ \log \frac{p_\theta(x, z, y \mid s)}{q_\phi(z \mid x, y, s)} \right]$$

$$= \mathbb{E}_{q_\phi(z \mid x, y, s)}\left[ \log p_\theta(x \mid z, s) \right] - D_{\text{KL}}\left( q_\phi(z \mid x, y, s) \parallel p(z \mid y) \right) + \log p(y) \equiv \mathcal{L}_{\text{ELBO}}(x, y \mid s)$$

여기에 실제 관측된 정답 레이블 $y$를 분류기가 정확하게 예측하도록 유도하는 지도 학습 교차 엔트로피 항을 추가하여 최종 레이블 손실을 정의한다:

$$\mathcal{L}_{\text{labeled}}(x, y, s) = -\mathcal{L}_{\text{ELBO}}(x, y \mid s) - \alpha \cdot \log q_\phi(y \mid x, s)$$

*   $\alpha > 0$: 비지도 생성 손실과 지도 분류 손실 간의 균형을 조절하는 하이퍼파라미터 (scvi-tools `classification_ratio`의 기본값은 50).

### 3.2 비레이블 세포(Unlabeled Cells)에 대한 주변화 손실 (Marginalization Loss)

세포형 $y$가 관측되지 않은 임상 쿼리 세포의 경우, $y$는 잠재 변수이므로 가능한 모든 $K$개의 세포형 범주에 대해 합산(Marginalization)을 수행해야 한다:

$$\log p_\theta(x \mid s) = \log \sum_{k=1}^K p_\theta(x, y=k \mid s)$$

$$\ge \mathbb{E}_{q_\phi(y \mid x, s)}\left[ \mathbb{E}_{q_\phi(z \mid x, y, s)}\left[ \log \frac{p_\theta(x, z, y \mid s)}{q_\phi(z, y \mid x, s)} \right] \right]$$

$$= \sum_{k=1}^K q_\phi(y=k \mid x, s) \left[ \mathcal{L}_{\text{ELBO}}(x, y=k \mid s) \right] - \sum_{k=1}^K q_\phi(y=k \mid x, s) \log q_\phi(y=k \mid x, s)$$

$$= \sum_{k=1}^K q_\phi(y=k \mid x, s) \, \mathcal{L}_{\text{ELBO}}(x, y=k \mid s) + \mathcal{H}\left( q_\phi(y \mid x, s) \right)$$

*   $\mathcal{H}(q) = - \sum_{k=1}^K q_k \log q_k$: 분류기 사후확률의 **샤논 엔트로피(Shannon Entropy)**.

따라서 비레이블 세포에 대한 음의 로그 우도 최소화 목적 함수는 다음과 같다:

$$\mathcal{L}_{\text{unlabeled}}(x, s) = -\sum_{k=1}^K q_\phi(y=k \mid x, s) \, \mathcal{L}_{\text{ELBO}}(x, y=k \mid s) - \mathcal{H}\left( q_\phi(y \mid x, s) \right)$$

### 3.3 scANVI 전역 목적 함수 (Total Objective)

전체 학습 데이터셋에 대한 목적 함수는 두 손실의 가중 합산으로 귀결된다:

$$\mathcal{L}_{\text{scANVI}}(\theta, \phi) = \sum_{(x, y, s) \in \mathcal{D}_L} \mathcal{L}_{\text{labeled}}(x, y, s) + \sum_{(x, s) \in \mathcal{D}_U} \mathcal{L}_{\text{unlabeled}}(x, s)$$

이 목적 함수는 비레이블 세포의 전사체 패턴이 잠재 공간의 기하학적 매니폴드를 보강하도록 강제하며, 레이블 정보는 잠재 공간의 세포형 클러스터 경계를 선명하게 정렬(Clustering Alignment)하는 역할을 수행한다.

---

## 4. 학습 전략: 콜드 스타트 방지와 scVI 기반 전이 학습 (Warm-Start Transfer)

scANVI 모델을 무작위 가중치로 초기화하여 처음부터 학습(Cold-start)시키면 심각한 **분류기 국소 최적화 편향(Classifier Bias)**에 빠진다:
*   학습 초기에는 잠재 공간 $z$의 매니폴드가 정립되지 않았기 때문에, 분류기 신경망이 몇몇 고발현 노이즈 유전자에 과적합되어 비레이블 세포들의 예측 확률을 특정 클래스로 쏠리게 만든다.
*   이는 비레이블 손실의 주변화 가중치 $q_\phi(y \mid x)$를 오염시켜 전체 인코더와 디코더를 파괴적인 최적화 경로로 끌고 간다.

```
[Two-Stage Warm-Start Pipeline: scVI to scANVI]

Stage 1: Pure Unsupervised scVI Pretraining (All Cells: Labeled + Unlabeled)
  Raw Counts X ──> scVI VAE ──> Stable Latent Manifold z ~ N(0, I)
  * 잠재 기하학 정립, 배치 효과 선제거, 라이브러리 스케일 분리

Stage 2: scANVI Initialization & Semi-Supervised Fine-Tuning
  scVI Weights (Encoder, Decoder, Dispersions) ──(Warm-Start)──> scANVI Model
  * Labeled Data: Cross-Entropy Supervision
  * Unlabeled Data: Marginalized ELBO Reinforcement
```

따라서 프로덕션 파이프라인에서는 항상 **2단계 전이 학습(Two-stage Warm-start)**을 수행한다:
1.  모든 데이터(레이블 + 비레이블)에 대해 비지도 scVI를 먼저 학습시켜 안정적인 기하학적 매니폴드와 음이항 파라미터를 수렴시킨다.
2.  `scvi.model.SCANVI.from_scvi_model`을 호출하여 사전학습된 가중치를 scANVI로 복사한 후, 분류기 모듈과 세포형별 사전분포를 추가하여 준지도 미세조정(Fine-tuning)을 수행한다.

---

## 5. 임상 케이스 스터디: mCRPC 코호트에서 희귀 t-NEPC 고정밀 식별

### 5.1 임상적 배경 및 문제의 특수성

전립선 선암종(PRAD) 환자에게 차세대 안드로겐 수용체 표적 치료제(Enzalutamide, Abiraterone)를 투여하면, 종양 세포는 강력한 선택압을 받아 안드로겐 신호전달계를 차단하고 신경내분비 표현형으로 탈분화(Lineage Plasticity / Transdifferentiation)하는 **치료유도 신경내분비 전립선암(t-NEPC)**으로 진화한다.

*   **주요 마커 발현 특성**:
    *   *AR*, *KLK3*(PSA), *PRAC1* 발현의 완전한 침묵(Loss).
    *   신경내분비 특이 분비 과립 마커 *SYP*(Synaptophysin), *CHGA*(Chromogranin A), 세포 접착 분자 *NCAM1*(CD56)의 급격한 상향 조절.
    *   마스터 신경 전사인자 *ASCL1*, *NEUROD1* 과발현 및 *TP53*/*RB1* 유전체 이중 소실.
*   **단일세포 분석의 난제**: 전이 생검 검체에서 t-NEPC 세포는 전체 생검 조직에서 극히 드물게 존재할 수 있으며(여기서 쓰는 0.1~0.5% 같은 비율은 출처 미확인 예시), 환자 간 전사체 편차가 극심하여 고전적 군집화로는 검출이 불가능하다.

### 5.2 프로덕션 `scANVI` 파이프라인 구현 코드

```python
"""
03_scanvi_tnepc_annotation.py
scVI 사전학습 및 scANVI 준지도 전이 학습을 통한 mCRPC 코호트 내 희귀 t-NEPC 고정밀 검출.
"""

import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def run_scanvi_tnepc_detection(adata_path: str):
    # 1. 통합 AnnData 로드
    adata = sc.read_h5ad(adata_path)
    print(f"Loaded Cohort: {adata.n_obs} cells across {adata.obs['patient_id'].nunique()} patients.")

    # 2. 레이블 비대칭 상태 점검
    # 'cell_type' 컬럼에 참조 아틀라스 세포는 실제 레이블, 임상 쿼리 세포는 'Unknown'으로 마킹
    print("Label distribution before scANVI:")
    print(adata.obs["cell_type"].value_counts(dropna=False))

    # 3. 데이터 등록 및 1단계: 비지도 scVI 사전학습 (Warm-start Stage)
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="patient_id",
        labels_key="cell_type" # 참조 라벨 등록
    )

    print(">>> Stage 1: Pretraining unsupervised scVI...")
    scvi_model = scvi.model.SCVI(
        adata,
        n_hidden=128,
        n_latent=20,
        n_layers=2,
        gene_likelihood="nb"
    )
    scvi_model.train(
        max_epochs=150,
        early_stopping=True,
        early_stopping_patience=15,
        plan_kwargs={"lr": 1e-3, "n_epochs_kl_warmup": 30}
    )

    # 4. 2단계: scVI 가중치를 상속받아 scANVI 인스턴스 구축
    print(">>> Stage 2: Initializing scANVI from pretrained scVI weights...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category="Unknown",
        labels_key="cell_type"
    )

    # 5. scANVI 준지도 미세조정 (Fine-tuning)
    scanvi_model.train(
        max_epochs=60,
        early_stopping=True,
        early_stopping_patience=10,
        plan_kwargs={"lr": 5e-4}
    )

    # 6. 임상 쿼리 비레이블 세포에 대한 세포형 예측 및 불확실성(엔트로피) 계산
    adata.obs["predicted_cell_type"] = scanvi_model.predict(adata)
    
    # 각 세포형에 대한 사후확률 소프트맥스 행렬 추출 (Cells x K)
    posterior_probs = scanvi_model.predict(adata, soft=True)
    
    # 샤논 엔트로피 계산: H(p) = - sum(p * log(p))
    probs_matrix = posterior_probs.values
    probs_clamped = np.clip(probs_matrix, 1e-8, 1.0)
    adata.obs["prediction_entropy"] = -np.sum(probs_clamped * np.log(probs_clamped), axis=1)

    # 7. 희귀 t-NEPC 예측 서브클론 검증
    tnepc_mask = (adata.obs["predicted_cell_type"] == "t-NEPC") & (adata.obs["cell_type"] == "Unknown")
    n_detected_tnepc = tnepc_mask.sum()
    # 주의: 정답 라벨이 없으므로 이 숫자는 '예측된 후보 수'일 뿐 정밀도·재현율이 아니다.
    #       검증은 독립 truth 라벨로 PR·혼동행렬을 계산해야 한다.
    print(f"Predicted {n_detected_tnepc} candidate t-NEPC cells among unlabeled queries")

    # 8. 핵심 바이오마커 발현 정량 검증
    markers = ["AR", "KLK3", "SYP", "CHGA", "NCAM1", "ASCL1"]
    print("=== Expression of Key Lineage Markers in Predicted t-NEPC ===")
    # scvi_normalized 레이어는 이 함수에서 만들지 않았으므로 모델에서 직접 계산한다.
    # get_normalized_expression은 로그 값이 아니라 (비율 x library_size) 스케일의 값이므로 expm1을 적용하지 않는다.
    present = [m for m in markers if m in adata.var_names]
    if n_detected_tnepc > 0 and present:
        norm = scanvi_model.get_normalized_expression(
            adata,
            indices=np.where(tnepc_mask.values)[0],
            gene_list=present,
            library_size=1e4,
            return_mean=True,
        )
        for marker in present:
            print(f" Marker {marker:<8} : Mean Denoised Expression = {norm[marker].mean():.3f}")

    # 9. 잠재 표현 추출 및 저장
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)
    
    return adata, scanvi_model
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: scVI 사전학습(Warm-start)과 콜드 스타트(Cold-start)의 최적화 동역학 차이
연구팀이 20만 개의 단일세포 전사체 코호트(참조 아틀라스 3만 개, 미분류 임상 코호트 17만 개)에 scANVI를 적용하고자 한다. 주니어 연구원이 "학습 시간을 단축하기 위해 scVI 사전학습 단계를 건너뛰고, 무작위 초기화 상태에서 scANVI를 곧바로 학습(Cold-start)시키겠다"고 제안하였다.
1. 콜드 스타트 학습 시 분류기 네트워크 $q_\phi(y \mid x)$의 초기 무작위 예측이 비레이블 주변화 손실 $\sum_{y} q(y \mid x) \mathcal{L}(x, y)$을 통해 전체 잠재 공간 매니폴드를 어떻게 왜곡시키는지 수학적으로 설명하시오.
2. 1단계 비지도 scVI 사전학습이 제공하는 '기하학적 앵커(Geometric Anchor)'의 역할을 정보 엔트로피와 잠재 사후분포 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 콜드 스타트 시의 왜곡 메커니즘**:
*   학습 초기 인코더와 분류기의 가중치가 무작위 상태일 때, 분류기 $q_\phi(y \mid x)$는 고유한 생물학적 마커 패턴을 인식하지 못하고 유전체 전반의 무작위 노이즈나 고발현 유전자(리보솜, 미토콘드리아 등)에 의해 임의의 클래스에 높은 가중치를 배정한다.
*   비레이블 세포에 대한 손실 함수는 각 세포형 $k$에 대한 가중 합산이다:
    $$\mathcal{L}_{\text{unlabeled}} = \sum_{k=1}^K q_\phi(y=k \mid x) \mathcal{L}_{\text{ELBO}}(x \mid y=k) + \mathcal{H}(q)$$
*   만약 $q_\phi(y \mid x)$가 특정 잘못된 클래스 $k'$로 조기에 쏠리게 되면(Premature Confidence), 디코더와 인코더는 해당 세포를 강제로 세포형 $k'$의 잠재 사전분포 $p(z \mid y=k')$에 맞추도록 그래디언트를 역전파한다.
*   이로 인해 잘못된 분류 정보가 잠재 공간의 기하학적 구조를 강제로 뒤틀고, 뒤틀린 잠재 공간이 다시 분류기의 성능을 저하시키는 **파괴적 양의 피드백 루프(Vicious Confirmation Bias Loop)**가 발생하여 심각한 모드 붕괴와 국소 최적해에 영구히 갇히게 된다.

**2. scVI 사전학습의 기하학적 앵커 역할**:
*   scVI는 세포형 레이블에 대한 어떠한 선입견 없이 오직 전사체 카운트의 복원 우도와 가우시안 정규화($z \sim \mathcal{N}(0, I)$)만을 최적화한다.
*   이를 통해 유전자 간의 비선형 공변동 네트워크와 배치 효과가 분리된 **자연스러운 연속 잠재 다양체(Continuous Latent Manifold)**가 사전에 확립된다.
*   이 안정된 매니폴드 위에서 scANVI를 초기화하면, 분류기 신경망 $q_\phi(y \mid x)$는 이미 생물학적으로 의미 있게 분리된 클러스터 경계를 학습하게 되므로 무작위 예측에 의한 왜곡이 발생하지 않는다. 따라서 비레이블 세포의 주변화 합산이 매끄러운 정규화 신호로 작동하여 빠르고 안정적인 전역 최적해로 수렴한다.

</details>

---

### 문제 2: 비레이블 주변화 손실의 계산 복잡도와 대규모 아틀라스 확장성
단일세포 아틀라스의 해상도가 극도로 정밀해져, 미세 아형(Subtypes)을 포함한 총 세포형 수 $K$가 200개에 도달하였다.
1. 비레이블 손실 함수 $\mathcal{L}_{\text{unlabeled}} = \sum_{k=1}^K q_\phi(y=k \mid x) \mathcal{L}_{\text{ELBO}}(x, y=k)$에서 $K$가 10에서 200으로 증가할 때 미니배치(Batch Size = 1,024) 순전파/역전파 연산량 및 GPU VRAM 요구량이 팽창하는 수학적 원인을 분석하시오.
2. 초대형 $K$ 환경에서 연산 복잡도를 $O(K)$에서 $O(M)$ ($M \ll K$)로 줄이기 위해 적용할 수 있는 상위 $M$개 확률 기반 가지치기(Top-$M$ Truncated Marginalization) 근사 기법을 수식으로 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 연산량 및 VRAM 팽창 원인**:
*   비레이블 세포 $x$ 하나에 대해, 모델은 $K$개의 모든 가능한 세포형 가설 $y \in \{1, \dots, K\}$에 대해 조건부 ELBO를 독립적으로 계산해야 한다:
    $$\mathcal{L}_{\text{ELBO}}(x, y=k) = \mathbb{E}_{q_\phi(z \mid x, k)}[\log p(x \mid z)] - D_{\text{KL}}(q_\phi(z \mid x, k) \parallel p(z \mid k))$$
*   즉, 미니배치 내의 세포 수 $B=1,024$에 대해 실제로 디코더와 인코더를 통과해야 하는 유효 텐서 크기는 $B \times K$가 된다.
*   $K=10$일 때는 배치당 10,240개의 가상 세포를 계산하지만, $K=200$이 되면 $1,024 \times 200 = 204,800$개의 세포에 대한 잠재 벡터 및 유전자 재구성 행렬($204,800 \times 3,000$ HVG $\approx 6.1 \times 10^8$ float 요소)을 생성해야 한다.
*   이는 순전파뿐만 아니라 역전파를 위한 계산 그래프 텐서 캐싱으로 인해 수십 기가바이트(GB)의 VRAM을 즉각 고갈시켜 `CUDA Out of Memory` 오류를 유발한다.

**2. Top-$M$ 절단 주변화 (Truncated Marginalization) 기법**:
*   분류기 $q_\phi(y \mid x)$가 출력하는 사후확률 분포에서, 대부분의 무관한 세포형(예: 상피세포 입력에 대한 T세포, B세포 확률)은 $q_k < 10^{-5}$로 거의 0에 수렴한다.
*   따라서 확률 상위 $M$개(e.g., $M = 3 \sim 5$)의 가장 유력한 세포형 집합 $\mathcal{S}_M(x) = \text{Top-M}_k \{q_\phi(y=k \mid x)\}$만을 동적으로 선별하고 확률을 재정규화한다:
    $$\tilde{q}_k = \frac{q_\phi(y=k \mid x)}{\sum_{j \in \mathcal{S}_M(x)} q_\phi(y=j \mid x)}, \quad k \in \mathcal{S}_M(x)$$
*   절단된 손실 함수는 다음과 같이 근사된다:
    $$\tilde{\mathcal{L}}_{\text{unlabeled}}(x) = \sum_{k \in \mathcal{S}_M(x)} \tilde{q}_k \, \mathcal{L}_{\text{ELBO}}(x, y=k) + \mathcal{H}(\tilde{q})$$
*   이 근사를 통해 계산 복잡도와 VRAM 사용량은 $O(K)$에서 $O(M)$으로 줄어든다. 다만 재정규화한 $\tilde{q}$와 $\mathcal{H}(\tilde{q})$를 쓰는 절단 근사는 원래 ELBO의 **편향 근사**이므로 정밀도를 보존한다고 보장할 수 없다(제안 수준의 기법이며 scvi-tools에는 구현되어 있지 않다).

</details>

---

### 문제 3: 사후 엔트로피 $\mathcal{H}(q)$의 임상생물학적 진단 가치와 Lineage Plasticity
mCRPC 환자 코호트에서 scANVI를 통해 세포형을 예측한 후, 각 세포의 분류 사후확률 분포 $q_\phi(y \mid x)$에 대한 샤논 엔트로피 $\mathcal{H}(q) = -\sum_k q_k \log q_k$를 산출하였다.
1. 엔트로피가 0에 가까운 세포군, 엔트로피가 매우 높은 세포군이 각각 임상생물학적으로 의미하는 바를 전립선암의 분화 가소성(Lineage Plasticity) 및 탈분화(Dedifferentiation) 관점에서 해석하시오.
2. 참조 아틀라스에 전혀 존재하지 않았던 신규 희귀 표현형(예: *AR* 음성이면서 신경내분비 마커도 발현하지 않는 이중 음성 전립선암, Double-Negative Prostate Cancer, DNPC)이 쿼리 데이터에 포함되어 있을 때, 엔트로피 $\mathcal{H}(q)$를 활용한 이상치 탐지(Out-of-Distribution, OOD Detection) 기준을 설계하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 엔트로피 수준별 임상생물학적 의미**:
*   **저엔트로피 ($\mathcal{H} \approx 0$)**: 특정 클래스의 사후확률이 $1.0$에 수렴하는 상태이다. 이는 분화 정체성이 확고하게 고정된 세포를 의미한다. 예를 들어, 전형적인 안드로겐 의존성 관강 선암종(*KLK3+*, *AR+*)이거나, 이미 신경내분비 분화가 완전히 종결된 성숙한 t-NEPC(*SYP+*, *CHGA+*) 세포군이다.
*   **고엔트로피 ($\mathcal{H} \gg 0$)**: 분류기가 관강 세포(Luminal)와 신경내분비 세포(t-NEPC) 사이에서 확신을 갖지 못하고 확률이 분산된 상태($q_{\text{Luminal}} \approx 0.5, q_{\text{t-NEPC}} \approx 0.5$)이다.
*   생물학적으로 이는 단순한 모델의 분류 실패가 아니라, 안드로겐 수용체 차단 치료의 압력하에 상피 정체성을 잃어가면서 신경내분비 분자 프로그램을 막 가동하기 시작한 **전이 중간 상태(Intermediate / Hybrid Transition State)**의 세포 가소성을 반영하는 후보 지표일 수 있다. 다만 고엔트로피는 배치 효과·더블릿·ambient RNA 등으로도 생기고 분류기 과신으로 신규 상태가 저엔트로피로 나올 수도 있으므로, 재구성 NLL·잠재 거리 등과 함께 해석해야 한다(결정적 바이오마커로 단정할 수 없다).

**2. OOD 탐지 및 신규 DNPC 클론 발굴 전략**:
*   **이중 음성 전립선암(DNPC)의 특성**: *AR* 경로가 차단되었으나 *SYP*나 *CHGA*와 같은 전형적 신경내분비 마커도 획득하지 않고 *FGF*, *MAPK* 신호전달계를 통해 생존하는 신규 치료 저항성 아형이다. 참조 아틀라스에 이 클래스가 등록되어 있지 않다면, 모델은 이를 강제로 기존 세포형 중 하나로 매핑하려 시도한다.
*   **OOD 탐지 기준 수립**:
    1.  **엔트로피 임계치**: 신규 세포형은 기지의 여러 세포형에 확률이 분산되어 사후 엔트로피가 높아질 수 있다: $\mathcal{H}(x) > \tau_{\text{entropy}}$ (e.g., 상위 5% 백분위수). 다만 소프트맥스 분류기는 분포 밖(OOD) 입력에도 과신하는 경우가 흔해 신규 세포형이 항상 고엔트로피가 되지는 않으며, 고엔트로피는 배치·더블릿·ambient RNA로도 생길 수 있다.
    2.  **재구성 우도 적합도 검정(Reconstruction Goodness-of-Fit)**: DNPC는 고유한 유전자 발현 시그니처를 지니므로, 기존 참조 아틀라스의 디코더 우도 $p_\theta(x \mid z)$로 복원했을 때의 음의 로그 우도(NLL, 재구성 오차)가 정상 세포군에 비해 현저히 높다: $\text{NLL}(x) > \tau_{\text{recon}}$.
*   두 조건을 동시에 만족하는 세포군($\mathcal{H} > \tau_{\text{entropy}} \land \text{NLL} > \tau_{\text{recon}}$)을 '미확인 신규 아형 후보군(Putative Novel Subtypes)'으로 격리하여 후속 무감독 클러스터링을 수행함으로써, 기존 아틀라스에 얽매이지 않고 임상적으로 치명적인 신규 저항성 클론(DNPC)을 선제적으로 발굴할 수 있다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Xu C, Lopez R, Mehlman E, et al. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Mol Syst Biol*. 2021;17(1):e9620. PMID: 33491336. doi: 10.15252/msb.20209620
2. Kingma DP, Rezende DJ, Mohamed S, et al. Semi-Supervised Learning with Deep Generative Models. *arXiv* preprint. 2014. arXiv: 1406.5298.
### 검수 근거 문헌
1. Lotfollahi M, Naghipourfar M, Luecken MD, et al. Mapping single-cell data to reference atlases by transfer learning. *Nat Biotechnol*. 2022;40(1):121-130. PMID: 34462589. doi: 10.1038/s41587-021-01001-7
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- scvi-tools 1.5.1 (github.com/scverse/scvi-tools). 원 논문: Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol*. 2022;40(2):163-166. PMID: 35132262. doi: 10.1038/s41587-021-01206-w
<!-- REFERENCES:END -->

# 06. 점진적 주성분 분석(Incremental PCA)과 무작위 SVD(Randomized SVD)

---

## 학습 목표 (Cognitive Goals)

1. **고전적 특이값 분해(Classical Full SVD)의 계산 복잡도 및 메모리 병목 수학적 증명**: 100만 개(1M) 단일세포 전립선암 아틀라스($N=10^6, G=35,000$) 규모에서 밀집 특이값 분해(Dense SVD)의 $O(N \cdot G^2)$ 시간 복잡도와 테라바이트급 메모리 요구량을 계산하고, 대규모 단일세포 분석에서 근사 선형 대수학(Approximate Linear Algebra)의 필연성을 설명할 수 있다.
2. **Ross et al. (2008)의 점진적 PCA(Incremental PCA / IPCA) 수리적 유도**: 미니배치 스트리밍 환경에서 전체 데이터의 평균 벡터를 축차 갱신($\boldsymbol{\mu}_{new}$)하고, 랭크-$k$ 부분공간에 신규 관측치 잔차를 증강한 저차원 QR 분해 및 SVD 갱신 알고리즘을 증명할 수 있다.
3. **Halko, Martinsson, & Tropp (2011)의 무작위 SVD(Randomized SVD) 정리 및 파워 반복 유도**: 가우시안 무작위 테스트 행렬 $\boldsymbol{\Omega} \in \mathbb{R}^{G \times (k + p)}$을 통한 레인지 탐색(Range Finder)과 완만한 스펙트럼 감쇠(Slow Spectral Decay)를 극복하는 파워 반복(Power Iteration, $(\mathbf{X}\mathbf{X}^T)^q \mathbf{X}\boldsymbol{\Omega}$)의 수학적 기전을 증명할 수 있다.
4. **특이값 재구성 오차 한계(Singular Value Reconstruction Error Bound) 정식화**: Eckart-Young-Mirsky 최적 저계수 근사 정리 대비 무작위 SVD의 기대 오차 상한 $\mathbb{E}\|\mathbf{X} - \mathbf{Q}\mathbf{Q}^T \mathbf{X}\|_2 \le (1 + \epsilon) \sigma_{k+1}$을 분석하고, 1M 전립선암 세포 아틀라스 벤치마크 환경에서 IPCA와 rSVD의 정확도-속도 트레이드오프를 평가할 수 있다.

---

## 1. 100만 세포 주성분 분석의 연산 장벽과 고전적 SVD의 붕괴

단일세포 RNA 시퀀싱(scRNA-seq/snRNA-seq) 파이프라인에서 주성분 분석(Principal Component Analysis, PCA)은 차원 축소, 이웃 그래프(kNN) 구축, 클러스터링 및 시각화(UMAP/t-SNE)의 필수적인 선행 단계이다. 중심화된 발현 행렬 $\tilde{\mathbf{X}} = \mathbf{X} - \mathbf{1}_N \boldsymbol{\mu}^T \in \mathbb{R}^{N \times G}$에 대해 상위 $k$개의 주성분을 추출하는 것은 $\tilde{\mathbf{X}}$의 절단된 특이값 분해(Truncated SVD)를 구하는 것과 수학적으로 동일하다:

$$\tilde{\mathbf{X}} = \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^T$$

여기서 $\mathbf{V}_k \in \mathbb{R}^{G \times k}$는 유전자의 로딩 벡터(Loadings), $\mathbf{U}_k \boldsymbol{\Sigma}_k \in \mathbb{R}^{N \times k}$는 각 세포의 주성분 임베딩 좌표(PC Embeddings)이다.

```
[The Computational Wall of Classical Full SVD on 1M Single Cells]

Matrix Dimensions:
  N = 1,000,000 cells (Prostate Cancer Atlas)
  G =    35,000 genes (or Top 3,000 HVGs)
  Target PCs (k) = 50

1. Classical Golub-Reinsch / Divide-and-Conquer SVD:
   Time Complexity : O(N * G^2) FLOPs
   For G = 35,000: 10^6 * (3.5 * 10^4)^2 = 1.225 x 10^15 FLOPs
     (지속 100 GFLOPS 가정 시 ~3.4시간, 상수 계수 무시; 메모리가 더 큰 문제)
   Even with G = 3,000 HVGs: 10^6 * (3 * 10^3)^2 = 9.0 x 10^12 FLOPs (같은 가정에서 ~90초)

2. Memory Allocation Bottleneck:
   Centering operation X - 1 * μ^T:
   Sparse Matrix (CSR, 98% sparse) ──> DENSE Centered Matrix!
   Dense Float32: 10^6 cells x 35,000 genes x 4 Bytes = 140 GB RAM!
   Left Singular Matrix U (full_matrices=True일 때만): 10^6 x 10^6 x 4 Bytes = 4 Terabytes! (Impossible)
   (경제형 SVD의 U는 N x G 크기: G=35,000이면 140 GB)
```

따라서 100만 세포 규모의 전사체 분석에서는 행렬 전체를 한 번에 메모리에 적재하여 분해하는 고전적 SVD 알고리즘이 완전히 붕괴한다. 이를 해결하기 위해 두 가지 상호 보완적인 선형대수학 알고리즘이 개발되었다:
1. **점진적 PCA (Incremental PCA, IPCA)**: 메모리 상한이 제한된 환경에서 배치 단위 스트리밍으로 부분공간을 점진적 갱신.
2. **무작위 SVD (Randomized SVD, rSVD)**: 확률론적 선형대수학(Randomized Linear Algebra)을 통해 $O(N \cdot G \cdot k)$의 최소 연산량으로 이론적 최적치에 근접하는 부분공간 고속 추출.

---

## 2. 점진적 PCA (Incremental PCA, IPCA)의 수학적 유도

Ross, Lim, Lin, and Yang (2008)이 제안한 Incremental PCA는 대용량 스트리밍 관측치를 처리하기 위해 개발된 알고리즘으로, 전체 데이터를 크기 $m$의 미니배치 단위로 순차 주입받으면서 고유공간(Eigenspace)과 평균 벡터를 점진적으로 갱신한다.

```
[Incremental PCA Online Subspace Update Pipeline]

Batch A (Prior n cells):
  Subspace: U_A (G x k), Singular Values: Σ_A (k x k), Mean: μ_A (G)
               │
               ▼
Incoming Batch B (New m cells):
  Data: X_B (m x G), Mean: μ_B (G)
               │
               ▼
[Step 1: Update Global Mean: μ_new = (n*μ_A + m*μ_B) / (n + m)]
               │
               ▼
[Step 2: Construct Augmented Residual Matrix Z ∈ R^{G x (k + m + 1)}]
  Z = [ U_A * Σ_A  |  (X_B - 1*μ_B^T)^T  |  sqrt(nm / (n+m)) * (μ_A - μ_B) ]
               │
               ▼
[Step 3: Thin QR Decomposition: Z = Q * R, where Q ∈ R^{G x (k+m+1)}]
               │
               ▼
[Step 4: Small SVD on Upper Triangular R: R = U_tilde * Σ_new * V_tilde^T]
               │
               ▼
[Step 5: Updated Subspace: U_new = Q * U_tilde[:, :k], Σ_new = Σ_new[:k, :k]]
(Memory consumption bounded by O(G * (k + m)), independent of total N!)
```

### 2.1 미니배치 평균 갱신 공식

이전에 처리된 $n$개의 세포 집합 $A$의 표본 평균을 $\boldsymbol{\mu}_A \in \mathbb{R}^G$, 새로 유입된 $m$개 세포 집합 $B$의 표본 평균을 $\boldsymbol{\mu}_B \in \mathbb{R}^G$라 하자.
통합된 $n + m$개 표본의 새로운 전역 평균 $\boldsymbol{\mu}_{\text{new}}$는 다음과 같이 정확히 갱신된다:

$$\boldsymbol{\mu}_{\text{new}} = \frac{n \boldsymbol{\mu}_A + m \boldsymbol{\mu}_B}{n + m}$$

두 집합 간의 평균 편차 벡터(Mean Shift Vector)를 다음과 같이 정의한다:

$$\Delta \boldsymbol{\mu} = \boldsymbol{\mu}_A - \boldsymbol{\mu}_B$$

### 2.2 중심화 보정 잔차 행렬과 저차원 QR 분해

집합 $A$와 $B$가 통합되었을 때의 전체 공분산 행렬 $\mathbf{C}_{\text{new}}$는 각 배치의 공분산과 평균 편차 외적의 합으로 분해된다:

$$(n + m) \mathbf{C}_{\text{new}} = n \mathbf{C}_A + m \mathbf{C}_B + \frac{nm}{n + m} \Delta \boldsymbol{\mu} \Delta \boldsymbol{\mu}^T$$

집합 $A$의 주성분 기저를 $\mathbf{U}_A \in \mathbb{R}^{G \times k}$, 특이값 대각 행렬을 $\boldsymbol{\Sigma}_A \in \mathbb{R}^{k \times k}$라 하면, $n \mathbf{C}_A \approx \mathbf{U}_A \boldsymbol{\Sigma}_A^2 \mathbf{U}_A^T$이다.
신규 배치 $B$의 중심화된 데이터 행렬을 $\tilde{\mathbf{X}}_B = \mathbf{X}_B - \mathbf{1}_m \boldsymbol{\mu}_B^T \in \mathbb{R}^{m \times G}$라 하면, $m \mathbf{C}_B = \tilde{\mathbf{X}}_B^T \tilde{\mathbf{X}}_B$이다.

따라서 증강 잔차 행렬(Augmented Matrix) $\mathbf{Z} \in \mathbb{R}^{G \times (k + m + 1)}$를 다음과 같이 수평 결합(Horizontal Concatenation)하여 구성한다:

$$\mathbf{Z} = \begin{bmatrix} \mathbf{U}_A \boldsymbol{\Sigma}_A & \tilde{\mathbf{X}}_B^T & \sqrt{\frac{nm}{n + m}} \Delta \boldsymbol{\mu} \end{bmatrix}$$

이 행렬 $\mathbf{Z}$에 대해 얇은 QR 분해(Thin QR Decomposition)를 수행한다:

$$\mathbf{Z} = \mathbf{Q} \mathbf{R}$$

- $\mathbf{Q} \in \mathbb{R}^{G \times (k + m + 1)}$: 정규직교 열벡터들을 갖는 기저 행렬 ($\mathbf{Q}^T \mathbf{Q} = \mathbf{I}$).
- $\mathbf{R} \in \mathbb{R}^{(k + m + 1) \times (k + m + 1)}$: 작은 크기의 상삼각 행렬.

### 2.3 소형 행렬 SVD 및 부분공간 절단

이제 전체 $G \times G$ 행렬 대신, 매우 작은 $(k + m + 1) \times (k + m + 1)$ 크기의 상삼각 행렬 $\mathbf{R}$에 대해서만 표준 특이값 분해를 수행한다:

$$\mathbf{R} = \tilde{\mathbf{U}} \boldsymbol{\Sigma}_{\text{new}} \tilde{\mathbf{V}}^T$$

최종 갱신된 상위 $k$개의 좌측 특이 벡터(유전자 로딩) $\mathbf{U}_{\text{new}}$와 특이값 $\boldsymbol{\Sigma}_{\text{new}}$는 다음과 같이 유도된다:

$$\mathbf{U}_{\text{new}} = \mathbf{Q} \tilde{\mathbf{U}}_{[:, 1:k]} \in \mathbb{R}^{G \times k}$$

$$\boldsymbol{\Sigma}_{\text{new}} = \boldsymbol{\Sigma}_{\text{new}[1:k, 1:k]} \in \mathbb{R}^{k \times k}$$

- **메모리 복잡도**: $O(G(k + m))$ — 전체 세포 수 $N$에 완전히 독립적이다.
- **시간 복잡도**: $N/m$개 배치에 대해 각 배치당 $O(G(k + m)^2)$이므로 총 $O(N \cdot G \cdot (k+m)^2 / m)$이며, $m \gg k$이면 약 $O(N \cdot G \cdot m)$이다($N$에 대해 선형; $O(N \cdot G \cdot k)$는 $m \approx k$일 때만 성립).

---

## 3. 무작위 SVD (Randomized SVD)의 수학적 원리와 Halko 정리

Halko, Martinsson, and Tropp (2011, *SIAM Review*)가 정립한 무작위 SVD(Randomized SVD)는 현대 대규모 머신러닝 및 단일세포 전사체학에서 가장 지배적으로 사용되는 고속 분해 기법이다. 무작위 SVD는 임의의 거대 행렬의 핵심 에너지가 소수의 지배적 고유공간에 집중되어 있다는 성질을 활용하여, 무작위 투영(Random Projection)을 통해 대상 부분공간의 작고 정밀한 정규직교 기저를 먼저 포착한다.

```
[Randomized SVD Algorithmic Architecture]

Input Matrix: X ∈ R^{N x G} (Centered or Sparse Log-normalized)
Target Rank: k (e.g., 50), Oversampling: p (e.g., 10), Dimension: l = k + p = 60

Stage A: Low-Rank Range Finder
  1. Generate Gaussian Random Test Matrix: Ω ∈ R^{G x l}, Ω_ij ~ N(0, 1)
  2. Subspace Projection: Y_0 = X * Ω ∈ R^{N x l}
  3. Power Iterations (q = 1 or 2, for slow decay):
     For j = 1 to q:
       Y_temp = X^T * Y_{j-1}  (Sparse matrix-matrix multiplication)
       Y_j    = X   * Y_temp
  4. Economic QR Decomposition: Y = Q * R ──> Q ∈ R^{N x l} (Q^T * Q = I_l)

Stage B: Direct Deterministic SVD
  1. Project X to Small Dense Subspace: B = Q^T * X ∈ R^{l x G} (l << min(N, G))
  2. Compute SVD of Tiny Matrix B: B = U_hat * Σ * V^T (U_hat ∈ R^{l x l})
  3. Reconstruct Global Left Singular Vectors: U = Q * U_hat ∈ R^{N x l}
  4. Truncate to top k components: U_k = U[:, :k], Σ_k = Σ[:k, :k], V_k = V[:, :k]
```

### 3.1 1단계: 무작위 레인지 탐색기 (Randomized Range Finder)

목표 랭크 $k$에 오버샘플링 파라미터(Oversampling Parameter) $p$ (통상 $p = 5 \sim 10$)를 더하여 $l = k + p$차원의 부분공간을 구성한다.
독립 항등 표준정규분포 원소로 구성된 무작위 가우시안 테스트 행렬 $\boldsymbol{\Omega} \in \mathbb{R}^{G \times l}$을 생성한다:

$$\Omega_{ij} \sim \mathcal{N}(0, 1)$$

표본 추출 행렬 $\mathbf{Y} \in \mathbb{R}^{N \times l}$는 행렬 곱으로 산출된다:

$$\mathbf{Y} = \mathbf{X} \boldsymbol{\Omega}$$

무작위 행렬 이론(Random Matrix Theory)에 의해, 가우시안 벡터는 등방성(Isotropic)을 가지므로 행렬 $\mathbf{X}$의 열공간(Column Space)에서 에너지가 가장 큰 상위 $k$개의 특이 벡터 방향과 직교할 확률이 측도 0(Measure Zero)이다. 따라서 $\mathbf{Y}$의 열들은 높은 확률로 상위 주성분 기저를 강력하게 포괄한다.

### 3.2 파워 반복 (Power Iteration)의 수리적 필요성

단일세포 전사체 카운트 행렬은 생물학적 전사 버스팅(Transcriptional Bursting) 및 기계적 드롭아웃(Dropout) 노이즈로 인해 특이값 스펙트럼이 지수적으로 급격히 감소하지 않고 매우 완만하게 감쇠(Slow Spectral Decay)하는 특성을 갖는다:

$$\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_k \approx \sigma_{k+1} \approx \dots \gg 0$$

이 경우 단순 투영 $\mathbf{Y} = \mathbf{X} \boldsymbol{\Omega}$는 불필요한 하위 노이즈 특이값 성분들에 오염된다. 이를 수치적으로 극복하기 위해 대칭 전치 거듭제곱 연산인 **파워 반복(Power Iteration)**을 적용한다:

$$\mathbf{Y} = (\mathbf{X} \mathbf{X}^T)^q \mathbf{X} \boldsymbol{\Omega}$$

특이값 분해 $\mathbf{X} = \sum_{j} \sigma_j \mathbf{u}_j \mathbf{v}_j^T$를 대입하면:

$$(\mathbf{X} \mathbf{X}^T)^q \mathbf{X} = \sum_{j=1}^{\min(N, G)} \sigma_j^{2q+1} \mathbf{u}_j \mathbf{v}_j^T$$

- 상위 특이값 $\sigma_1, \dots, \sigma_k$는 거듭제곱 승수 $2q+1$에 의해 기하급수적으로 증폭되는 반면,
- 하위 노이즈 특이값 $\sigma_j$ ($j > k$)는 상대 비율 $(\sigma_j / \sigma_k)^{2q+1}$에 수렴하여 급격히 소멸한다.
- 단일세포 전사체 데이터에서는 통상 $q=1$ 또는 $q=2$의 작은 반복수만으로도 기계적 노이즈를 완벽하게 분리할 수 있다.

### 3.3 정규직교 기저 $\mathbf{Q}$ 추출 및 소형 투영 행렬 $\mathbf{B}$ 구축

파워 반복을 거친 $\mathbf{Y}$에 대해 수치적으로 안정적인 얇은 QR 분해를 수행하여 정규직교 기저 $\mathbf{Q} \in \mathbb{R}^{N \times l}$를 추출한다:

$$\mathbf{Q} = \text{qr}(\mathbf{Y}) \quad (\mathbf{Q}^T \mathbf{Q} = \mathbf{I}_l)$$

이제 원본 행렬 $\mathbf{X}$를 기저 $\mathbf{Q}$가 형성하는 $l$차원 부분공간으로 직교 투영하여 축소 행렬 $\mathbf{B} \in \mathbb{R}^{l \times G}$를 형성한다:

$$\mathbf{B} = \mathbf{Q}^T \mathbf{X}$$

### 3.4 2단계: 소형 행렬 결정론적 SVD 및 복원

축소 행렬 $\mathbf{B}$의 행 크기 $l$은 통상 50~60에 불과하므로, 표준 LAPACK 결정론적 SVD를 즉각 수행할 수 있다:

$$\mathbf{B} = \hat{\mathbf{U}} \boldsymbol{\Sigma} \mathbf{V}^T$$

여기서 $\hat{\mathbf{U}} \in \mathbb{R}^{l \times l}$, $\boldsymbol{\Sigma} \in \mathbb{R}^{l \times l}$, $\mathbf{V} \in \mathbb{R}^{G \times l}$이다.
마지막으로 좌측 특이 벡터(세포 임베딩 방향)를 원본 $N$차원 공간으로 회전 복원한다:

$$\mathbf{U} = \mathbf{Q} \hat{\mathbf{U}} \in \mathbb{R}^{N \times l}$$

상위 $k$개의 열만을 슬라이싱하여 최종 Truncated SVD 결과를 도출한다.

---

## 4. 특이값 재구성 오차 한계 (Reconstruction Error Bound)의 증명

Eckart-Young-Mirsky 정리에 따르면, 임의의 행렬 $\mathbf{X}$에 대한 최적의 랭크-$k$ 근사 오차는 $k+1$번째 특이값 $\sigma_{k+1}$에 의해 하한(Lower Bound)이 결정된다:

$$\min_{\text{rank}(\mathbf{A}) \le k} \|\mathbf{X} - \mathbf{A}\|_2 = \sigma_{k+1}$$

Halko, Martinsson, and Tropp (2011)은 무작위 SVD가 산출한 투영 행렬 $\mathbf{Q}\mathbf{Q}^T \mathbf{X}$의 기대 스펙트럼 노름 오차에 대해 다음의 엄밀한 상한 정리를 증명하였다.

### 정리 (Halko et al., Theorem 1.1 & 10.6)

가우시안 무작위 행렬 $\boldsymbol{\Omega} \in \mathbb{R}^{G \times (k+p)}$ ($p \ge 2$)를 사용하여 구축한 정규직교 기저 $\mathbf{Q}$에 대해, 스펙트럼 노름(Spectral Norm) 하에서의 기대 근사 오차는 다음을 만족한다:

$$\mathbb{E} \|\mathbf{X} - \mathbf{Q}\mathbf{Q}^T \mathbf{X}\|_2 \le \left( 1 + \sqrt{\frac{k}{p - 1}} \right) \sigma_{k+1} + \frac{e\sqrt{k+p}}{p} \left( \sum_{j > k} \sigma_j^2 \right)^{1/2}$$

나아가 파워 반복 $q$회를 적용한 경우, 이 오차 상한은 다음과 같이 기하급수적으로 축소된다:

$$\mathbb{E} \|\mathbf{X} - \mathbf{Q}\mathbf{Q}^T \mathbf{X}\|_2 \le \left[ \left( 1 + \sqrt{\frac{k}{p - 1}} \right) \sigma_{k+1}^{2q+1} + \frac{e\sqrt{k+p}}{p} \left( \sum_{j > k} \sigma_j^{2(2q+1)} \right)^{1/2} \right]^{\frac{1}{2q+1}}$$

- **오버샘플링 $p$의 효과**: $p=10$만 설정해도 선행 계수 $\sqrt{\frac{k}{p-1}}$가 급격히 작아져 꼬리 특이값의 누적 영향이 억제된다.
- **파워 반복 $q$의 효과**: $q \ge 1$일 때, 괄호 내부의 상한에 $\frac{1}{2q+1}$ 승이 취해지므로 오차 상한이 $(1 + \epsilon) \sigma_{k+1}$ 쪽으로 빠르게 줄어든다(이는 기대 오차의 **상한**이며 실제 오차가 즉각 최적에 도달한다는 뜻은 아니다).

---

## 5. 1M 전립선암 세포 아틀라스 벤치마크 구현 및 비교

아래 파이썬 코드는 100만 개 전립선암 단일세포 전사체 가상 데이터셋에 대해, Scikit-learn의 `IncrementalPCA`와 Halko의 `Randomized SVD`를 직접 구현하여 시간/메모리 복잡도 및 부분공간 재구성 오차를 정밀하게 벤치마크한다.

```python
"""
High-Performance Incremental PCA & Randomized SVD Benchmark
for 1,000,000 Prostate Cancer Single-Cell Atlas
"""

import time
import gc
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import IncrementalPCA
from typing import Tuple, Dict, Any


def generate_synthetic_prostate_1m(n_cells: int = 1_000_000, n_genes: int = 3000, 
                                   rank_latent: int = 50, sparsity: float = 0.98) -> sp.csr_matrix:
    """
    1M 전립선암 단일세포 발현의 저계수 구조(Low-rank biological latent structure)를
    반영한 현실적 희소 CSR 카운트 행렬을 생성한다.
    """
    print(f"[Data Synthesis] Generating {n_cells:,} cells x {n_genes:,} genes (Sparsity: {sparsity*100:.1f}%)...")
    t0 = time.time()
    
    np.random.seed(42)
    # 기저 저차원 임베딩 생성 (세포 표현형 공간)
    U_true = np.random.randn(n_cells, rank_latent).astype(np.float32)
    V_true = np.random.randn(n_genes, rank_latent).astype(np.float32)
    
    # 청크 단위로 희소 행렬 빌드 (메모리 스필 방지)
    chunk_size = 100_000
    n_chunks = n_cells // chunk_size
    csr_list = []
    
    for i in range(n_chunks):
        u_chunk = U_true[i*chunk_size : (i+1)*chunk_size]
        # 저계수 신호 합성 + 노이즈
        signal = np.dot(u_chunk, V_true.T)
        # 포아송 버스팅 모사: 임계치 이하 마스킹
        threshold = np.percentile(signal, sparsity * 100)
        signal[signal < threshold] = 0.0
        # 희소 CSR 변환
        csr_chunk = sp.csr_matrix(signal, dtype=np.float32)
        csr_list.append(csr_chunk)
        
    X_csr = sp.vstack(csr_list, format="csr")
    print(f"[Data Synthesis] Completed in {time.time() - t0:.2f}s. Total NNZ: {X_csr.nnz:,}")
    return X_csr


def randomized_svd_sparse(X: sp.csr_matrix, k: int = 50, n_iter: int = 2, 
                          oversampling: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Halko et al. (2011) 알고리즘 기반 대규모 희소 행렬 전용 무작위 SVD.
    중심화(Centering)를 명시적 밀집 변환 없이 암묵적 연산자(Implicit Operator)로 처리.
    """
    t0 = time.time()
    N, G = X.shape
    l = k + oversampling

    # 1. 유전자별 전역 평균 계산 (희소성 유지)
    mean_g = np.asarray(X.mean(axis=0)).ravel()  # shape (G,)

    # 암묵적 중심화 행렬-벡터 곱 연산자 정의:
    # A * v = X * v - 1 * (mean_g^T * v)
    # A^T * u = X^T * u - mean_g * (1^T * u)
    def matmul_centered(v: np.ndarray) -> np.ndarray:
        # v: shape (G, l)
        return X.dot(v) - np.outer(np.ones(N, dtype=np.float32), np.dot(mean_g, v))

    def rmatmul_centered(u: np.ndarray) -> np.ndarray:
        # u: shape (N, l)
        u_sum = np.sum(u, axis=0)  # shape (l,)
        return X.T.dot(u) - np.outer(mean_g, u_sum)

    # 2. Stage A: 무작위 가우시안 테스트 행렬 생성
    Omega = np.random.randn(G, l).astype(np.float32)

    # 3. 파워 반복을 통한 레인지 탐색 (Range Finder with Power Iterations)
    print(f"[rSVD] Starting Power Iteration (q = {n_iter}, subspace dim = {l})...")
    Y = matmul_centered(Omega)  # shape (N, l)

    for it in range(n_iter):
        Z = rmatmul_centered(Y)  # shape (G, l)
        # 직교화로 부동소수점 오버플로우 방어
        Z, _ = np.linalg.qr(Z, mode="reduced")
        Y = matmul_centered(Z)
        Y, _ = np.linalg.qr(Y, mode="reduced")

    # Q 기저 행렬 확정
    Q, _ = np.linalg.qr(Y, mode="reduced")  # shape (N, l)

    # 4. Stage B: 소형 축소 행렬 B 산출
    # B = Q^T * X_centered = (X_centered^T * Q)^T
    print("[rSVD] Projecting to core matrix B...")
    B = rmatmul_centered(Q).T  # shape (l, G)

    # 5. 소형 행렬 결정론적 SVD
    U_hat, s, Vt = np.linalg.svd(B, full_matrices=False)

    # 6. 좌측 특이 벡터 회전 복원 및 k차원 절단
    U = np.dot(Q, U_hat[:, :k])
    singular_values = s[:k]
    V = Vt[:k, :].T

    elapsed = time.time() - t0
    print(f"[rSVD] Finished in {elapsed:.2f}s.")
    return U, singular_values, V


def run_incremental_pca_benchmark(X: sp.csr_matrix, k: int = 50, batch_size: int = 50_000):
    """
    Scikit-learn의 IncrementalPCA를 활용한 미니배치 스트리밍 벤치마크
    """
    print(f"\n[IPCA] Starting IncrementalPCA (Batch Size: {batch_size:,}, Target PCs: {k})...")
    t0 = time.time()
    N, G = X.shape
    ipca = IncrementalPCA(n_components=k, batch_size=batch_size)

    n_batches = int(np.ceil(N / batch_size))
    for b in range(n_batches):
        start = b * batch_size
        end = min((b + 1) * batch_size, N)
        # 배치 슬라이스를 밀집 행렬로 변환 (IPCA는 밀집 배치 요구)
        batch_dense = X[start:end].toarray()
        ipca.partial_fit(batch_dense)
        del batch_dense

    elapsed = time.time() - t0
    print(f"[IPCA] Finished fitting in {elapsed:.2f}s.")
    return ipca


def benchmark_subspace_fidelity(U_rsvd: np.ndarray, V_ipca: np.ndarray, V_rsvd: np.ndarray):
    """
    두 알고리즘이 추출한 주성분 부분공간의 일치도(Grassmannian Subspace Overlap) 산출
    """
    # 주성분 유전자 로딩 간의 부분공간 사잇각 (Principal Angles) 분석
    # Basis: V_rsvd (G x k), V_ipca (k x G) -> V_ipca.T (G x k)
    V1, _ = np.linalg.qr(V_rsvd)
    V2, _ = np.linalg.qr(V_ipca.T)
    
    # 부분공간 투영 오버랩 메트릭: ||V1^T * V2||_F^2 / k (1.0이면 완벽 일치)
    overlap = np.linalg.norm(np.dot(V1.T, V2), ord="fro")**2 / V1.shape[1]
    print(f"[Subspace Overlap] Grassmannian Fidelity Metric: {overlap:.4f} / 1.0000")
    return overlap
```

---

## 6. 전립선암 1M 세포 아틀라스 알고리즘 비교 결과

100만 개 단일세포 전립선암 아틀라스($10^6 \times 3,000$ HVG, 50 PCs)에 대한 알고리즘 비교 표이다. **아래 시간·메모리·정밀도 수치는 이 레포에서 측정한 값이 아니라 개념 비교를 위한 예시(가상) 수치이며**, 메모리 열은 $10^6 \times 3,000$ float32 밀집 행렬(= 12 GB)을 기준으로 다시 계산했다(140 GB는 $G=35,000$일 때의 값).

| 비교 분석 항목 | 고전적 밀집 Full SVD | 점진적 PCA (Incremental PCA) | 무작위 SVD (Randomized SVD) |
| :--- | :--- | :--- | :--- |
| **시간 복잡도 (Time)** | $O(N \cdot G^2)$ (수 시간 이상) | $O(N \cdot G \cdot k)$ (약 14분) | **$O(N \cdot G \cdot (k+p) \cdot q)$ (약 1분 48초)** |
| **최대 메모리 사용량 (Peak RAM)** | **12 GB 이상** (float32 밀집 중심화 행렬만, SVD 작업 공간 별도) | **1.8 GB** (배치 5만 개 기준) | **4.2 GB** (희소 행렬 상주) |
| **행렬 희소성 보존 여부** | 완전 손실 (밀집 중심화 강제) | 배치 단위 밀집화 필요 | **완벽 보존 (암묵적 연산자 적용)** |
| **수치적 정밀도 (Fidelity)** | 이론적 100% (기준점) | 99.82% (배치 평균 근사) | **99.96% (Halko 오차 상한 도달)** |
| **분산 환경 확장성** | 불가능 | 순차 배치 의존적 (직렬) | **완전 병렬화 가능 (Map-Reduce)** |

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 무작위 SVD의 파워 반복(Power Iteration) 횟수 $q$에 따른 스펙트럼 수렴 해석
단일세포 전사체 데이터에서 유전자 발현 행렬 $\mathbf{X}$의 특이값 스펙트럼이 완만한 지수 감소 패턴 $\sigma_j = j^{-0.5}$를 따를 때:
1. $q=0$ (파워 반복 미수행)인 단순 무작위 SVD를 적용할 경우, 오버샘플링 $p=10$, 목표 차원 $k=50$에서 재구성 오차 한계 수식의 꼬리 특이값 누적 항 $\left( \sum_{j > k} \sigma_j^2 \right)^{1/2}$이 발산하거나 거대해지는 수학적 원인을 조화급수(Harmonic Series) 관점에서 설명하시오.
2. $q=2$의 파워 반복을 적용했을 때, 유효 특이값 스펙트럼 $\sigma_j^{2q+1} = \sigma_j^5$이 꼬리 노이즈($j > 50$)의 기여도를 억제하고 오차 한계를 이론적 최적치 $\sigma_{k+1}$에 근접시키는 메커니즘을 증명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. $q=0$일 때 꼬리 특이값 누적 항의 팽창과 조화급수 발산성**:
특이값이 $\sigma_j = j^{-0.5}$로 천천히 감쇠할 때, 꼬리 특이값 제곱합은 다음과 같이 적분 판정법으로 근사된다:

$$\sum_{j=k+1}^G \sigma_j^2 = \sum_{j=51}^{35,000} \left( j^{-0.5} \right)^2 = \sum_{j=51}^{35,000} \frac{1}{j} \approx \ln(35,000) - \ln(50) = \ln(700) \approx 6.55$$

Halko의 오차 한계에서 $q=0$일 때의 기대 오차는 다음과 같다:

$$\mathbb{E} \|\mathbf{X} - \mathbf{Q}\mathbf{Q}^T \mathbf{X}\|_2 \le \left( 1 + \sqrt{\frac{k}{p - 1}} \right) \sigma_{k+1} + \frac{e\sqrt{k+p}}{p} \sqrt{\sum_{j > k} \sigma_j^2}$$

여기서 $k=50, p=10$을 대입하면:
- 제1항 계수: $1 + \sqrt{50/9} \approx 1 + 2.36 = 3.36$, $\sigma_{51} = 51^{-0.5} \approx 0.140 \implies 3.36 \times 0.140 \approx 0.470$
- 제2항: $\frac{2.718 \times \sqrt{60}}{10} \times \sqrt{6.55} \approx \frac{21.05}{10} \times 2.56 \approx 5.39$

제2항(꼬리 노이즈 누적)의 크기($5.39$)가 제1항(목표 특이값 크기 $0.47$) 대비 11배 이상 거대해진다. 즉, 조화급수적으로 천천히 감소하는 단일세포 전사체 노이즈가 누적되어 오차 **상한**이 크게 커지므로(상한이 크다고 실제 실패가 보장되는 것은 아니다), $q=0$인 단순 투영은 상위 50개 고유공간 포착에 부분공간 왜곡 위험이 큰 설정인다.

**2. $q=2$ 파워 반복 적용 시 꼬리 노이즈 억제 증명**:
$q=2$를 적용하면 유효 행렬 $(\mathbf{X}\mathbf{X}^T)^2 \mathbf{X}$의 특이값은 $\tilde{\sigma}_j = \sigma_j^{2q+1} = \sigma_j^5 = j^{-2.5}$로 변환된다.
변환된 꼬리 제곱합을 계산하면:

$$\sum_{j=k+1}^G \tilde{\sigma}_j^2 = \sum_{j=51}^{35,000} \left( j^{-2.5} \right)^2 = \sum_{j=51}^{35,000} j^{-5} \approx \int_{50}^\infty x^{-5} dx = \left[ -\frac{1}{4} x^{-4} \right]_{50}^\infty = \frac{1}{4 \times 50^4} = \frac{1}{25,000,000} = 4 \times 10^{-8}$$

이 값의 제곱근은 $\sqrt{4 \times 10^{-8}} = 2 \times 10^{-4}$에 불과하다!
반면 51번째 유효 특이값은 $\tilde{\sigma}_{51} = 51^{-2.5} \approx 5.3 \times 10^{-5}$이다.
최종 오차 한계는 전체에 $\frac{1}{2q+1} = \frac{1}{5}$ 승을 취하므로:

$$\mathbb{E} \|\mathbf{X} - \mathbf{Q}\mathbf{Q}^T \mathbf{X}\|_2 \le \left[ C_1 \cdot \sigma_{k+1}^5 + C_2 \cdot 10^{-4} \right]^{1/5} \approx 1.6\,\sigma_{k+1} \quad (\text{이 수치 예의 상한 값})$$

즉, 단 2회의 파워 반복만으로도 느린 스펙트럼 감쇠를 급격한 멱급수 감쇠($j^{-5}$)로 강제 변환시켜, 오차 상한의 꼬리 항을 크게 줄여 최적 근사치에 가까운 상한을 준다(계산: $q=0$ 상한 $0.470+5.39 \approx 5.86$, $q=2$ 상한 $\approx (6\times10^{-4})^{1/5} \approx 0.227 \approx 1.6\,\sigma_{51}$로 최적 $\sigma_{51}=0.140$에 완전히 도달하지는 않는다).

</details>

---

### 문제 2: 대규모 희소 행렬 무작위 SVD에서 암묵적 중심화(Implicit Centering)의 중요성
100만 세포 발현 행렬 $\mathbf{X}$ (희소도 98%, 크기 $10^6 \times 3,000$)에 대해 PCA를 수행하기 위해 평균을 0으로 맞추는 중심화(Centering) 과정이 필수적이다.
1. 명시적으로 $\tilde{\mathbf{X}} = \mathbf{X} - \mathbf{1}_N \boldsymbol{\mu}^T$를 계산하여 변수에 저장하는 순간 발생하는 시스템 메모리 붕괴 과정을 설명하시오.
2. 무작위 SVD의 핵심 연산인 $\mathbf{Y} = \tilde{\mathbf{X}} \boldsymbol{\Omega}$와 $\mathbf{Z} = \tilde{\mathbf{X}}^T \mathbf{Y}$를 수행할 때, 원본 희소 행렬 $\mathbf{X}$의 희소성을 100% 보존하면서 암묵적으로 중심화 행렬 곱을 수행하는 결합 법칙(Associative Property) 수식을 전개하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 명시적 중심화 시 시스템 메모리 붕괴 메커니즘**:
- 원본 행렬 $\mathbf{X}$는 $10^6 \times 3,000$ 크기의 희소 CSR 행렬로, 98%가 0이다. 비영 원소는 $6 \times 10^7$개이며 메모리 점유율은 약 $480 \text{ MB}$에 불과하다.
- 각 유전자의 표본 평균 벡터 $\boldsymbol{\mu} \in \mathbb{R}^{3000}$은 대부분 0이 아닌 양수 값을 갖는다.
- 따라서 $\tilde{\mathbf{X}}_{ig} = X_{ig} - \mu_g$를 계산하면, 기존에 0이었던 98%의 모든 원소가 $-\mu_g \ne 0$으로 채워지면서 행렬의 희소도가 0%인 **완전 밀집 행렬(Dense Matrix)**로 즉각 전환된다.
- 밀집 행렬의 메모리 요구량은:
  $$1,000,000 \times 3,000 \times 4 \text{ Bytes (float32)} = 12,000,000,000 \text{ Bytes} = 12.0 \text{ GB}$$
- 파이썬 연산 과정에서 중간 복사본이 생성되거나 배속 정밀도(float64)로 승격될 경우 $24 \text{ GB} \sim 36 \text{ GB}$ 이상의 RAM이 순간적으로 점유되어 워커 노드의 OOM Killer가 작동한다.

**2. 암묵적 중심화 행렬 곱의 결합 법칙 수식 전개**:
$\boldsymbol{\Omega} \in \mathbb{R}^{G \times l}$ ($l \ll G$, e.g., $l=60$)에 대해 행렬 곱의 분배 법칙 및 결합 법칙을 적용한다:

$$\mathbf{Y} = \tilde{\mathbf{X}} \boldsymbol{\Omega} = (\mathbf{X} - \mathbf{1}_N \boldsymbol{\mu}^T) \boldsymbol{\Omega} = \mathbf{X} \boldsymbol{\Omega} - \mathbf{1}_N (\boldsymbol{\mu}^T \boldsymbol{\Omega})$$

- $\mathbf{X} \boldsymbol{\Omega}$: 희소 행렬 $\mathbf{X}$와 밀집 소형 행렬 $\boldsymbol{\Omega}$의 곱으로, 희소 CSR 블록 연산(`sp.csr_matrix.dot`)을 통해 비영 원소에 대해서만 연산이 수행되므로 메모리 팽창 없이 초고속으로 $N \times l$ 밀집 행렬을 산출한다.
- $\boldsymbol{\mu}^T \boldsymbol{\Omega}$: $1 \times G$ 벡터와 $G \times l$ 행렬의 곱으로, $O(G \cdot l)$의 무시할 수 있는 연산으로 $1 \times l$ 행 벡터 $\mathbf{w}^T$를 생성한다.
- $\mathbf{1}_N \mathbf{w}^T$: 크기 $N \times l$의 브로드캐스트 외적이므로, 최종 $\mathbf{Y}$는 $\mathbf{X}\boldsymbol{\Omega}$의 각 열에서 스칼라 $w_j$를 빼주는 단순 In-place 연산으로 완료된다.

마찬가지로 전치 행렬 곱 $\mathbf{Z} = \tilde{\mathbf{X}}^T \mathbf{Y}$ 연산 역시:

$$\mathbf{Z} = (\mathbf{X} - \mathbf{1}_N \boldsymbol{\mu}^T)^T \mathbf{Y} = (\mathbf{X}^T - \boldsymbol{\mu} \mathbf{1}_N^T) \mathbf{Y} = \mathbf{X}^T \mathbf{Y} - \boldsymbol{\mu} (\mathbf{1}_N^T \mathbf{Y})$$

- $\mathbf{1}_N^T \mathbf{Y}$: $N \times l$ 행렬의 열별 합계(Column Sum)를 구하여 $1 \times l$ 벡터 $\mathbf{s}^T$를 산출 ($O(N \cdot l)$).
- $\boldsymbol{\mu} \mathbf{s}^T$: $G \times 1$ 벡터와 $1 \times l$ 벡터의 외적이므로 크기 $G \times l$ 행렬이 생성된다.
- $\mathbf{X}^T \mathbf{Y}$: 희소 전치 행렬과 $N \times l$ 밀집 행렬의 곱으로 $G \times l$ 행렬 생성.

이와 같이 결합 순서를 변경함으로써, 단 한 번도 $N \times G$ 밀집 행렬을 메모리에 생성하지 않고 오직 희소 행렬과 $N \times l$ 크기의 저차원 버퍼만을 사용하여 완벽하게 중심화된 SVD를 수행할 수 있다.

</details>

---

### 문제 3: 점진적 PCA(IPCA)와 무작위 SVD(rSVD)의 알고리즘 선택 트레이드오프
전립선암 다기관 임상 시험 코호트 데이터를 분석하는 생물정보학 엔지니어링 환경에서:
1. 전국 10개 병원에서 매주 10,000개 세포씩 실시간으로 추가 시퀀싱 데이터가 지속적으로 유입되는 연속 학습(Continual / Online Learning) 파이프라인에서 IPCA가 rSVD 대비 구조적으로 우월한 이유를 설명하시오.
2. 반면 이미 구축이 완료된 100만 세포 정적 아틀라스(Static Atlas) 데이터셋에 대해 클라우드 분산 노드(Dask 클러스터)를 활용하여 상위 50개 주성분을 최단 시간 내에 추출해야 할 때 rSVD가 압도적으로 권장되는 분산 컴퓨팅(Map-Reduce) 관점의 이유를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 실시간 온라인 유입 환경에서 IPCA의 우월성**:
- IPCA는 과거 $n$개 세포에 대해 구축된 부분공간 기저 $\mathbf{U}_A$, 특이값 $\boldsymbol{\Sigma}_A$, 평균 $\boldsymbol{\mu}_A$만을 보존하고 있다면, 신규 유입된 $m$개 세포 $\mathbf{X}_B$가 도착하는 즉시 이전 원본 데이터를 다시 읽어올 필요 없이 $O(G(k+m)^2)$ 시간 내에 고유공간을 축차 갱신(Online State Update)할 수 있다.
- 반면 rSVD는 원본 데이터 행렬 전체 $\mathbf{X}$에 접근하여 가우시안 투영 및 파워 반복($\mathbf{X} \boldsymbol{\Omega}$, $\mathbf{X}^T \mathbf{Y}$)을 여러 번 수행해야 하는 배치 알고리즘(Batch Algorithm)이다.
- 따라서 데이터가 스트리밍으로 유입될 때마다 rSVD를 재실행하려면 기존 누적 데이터와 신규 데이터를 합친 전체 코호트를 매번 처음부터 다시 분해해야 하므로 I/O 및 연산 비용이 누적 표본 수에 비례하여 $O(N_{\text{total}} \cdot G \cdot k)$로 지속 증가하는 비효율이 발생한다. 따라서 온라인 스트리밍 환경에는 IPCA가 구조적으로 최적이다.

**2. 정적 대규모 아틀라스 분산 분석에서 rSVD의 압도적 우월성 (Map-Reduce 친화성)**:
- rSVD의 가장 큰 수학적 장점은 연산의 95% 이상이 **행렬-행렬 곱(GEMM / SpGEMM)**으로 구성되어 있다는 점이다:
  $$\mathbf{Y} = \mathbf{X} \boldsymbol{\Omega} = \begin{bmatrix} \mathbf{X}_1 \\ \mathbf{X}_2 \\ \vdots \\ \mathbf{X}_P \end{bmatrix} \boldsymbol{\Omega} = \begin{bmatrix} \mathbf{X}_1 \boldsymbol{\Omega} \\ \mathbf{X}_2 \boldsymbol{\Omega} \\ \vdots \\ \mathbf{X}_P \boldsymbol{\Omega} \end{bmatrix}$$
- 행렬 $\mathbf{X}$가 $P$개의 분산 노드에 파티셔닝되어 있을 때, 각 워커 노드는 공유된 소형 행렬 $\boldsymbol{\Omega}$ ($3000 \times 60$, 수백 KB)만 브로드캐스트 받으면 통신(Communication) 없이 로컬 청크 $\mathbf{X}_p$에 대해 독립적으로 부분 곱을 즉시 수행(Embarrassingly Parallel Map)할 수 있다.
- 또한 전치 곱 $\mathbf{X}^T \mathbf{Y} = \sum_{p=1}^P \mathbf{X}_p^T \mathbf{Y}_p$ 연산 역시 각 워커가 $G \times l$ 크기의 작은 로컬 행렬을 계산한 후 전역 덧셈(All-Reduce Sum) 한 번으로 즉시 병합된다.
- 반면 IPCA는 수학적으로 직렬적(Sequential)이다. 배치 $B$의 QR 갱신은 반드시 직전 배치 $A$의 고유공간 갱신 결과가 완료되어야만 진행될 수 있으므로 분산 클러스터의 수십 개 워커 노드가 병렬로 동시에 참여할 수 없는 심각한 암달의 법칙(Amdahl's Law) 병목에 갇히게 된다. 따라서 대규모 정적 아틀라스의 초고속 병렬 추출에는 rSVD가 압도적으로 우수하다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Ross DA, Lim J, Lin RS, et al. Incremental Learning for Robust Visual Tracking. *Int J Comput Vis*. 2008;77(1-3):125-141. doi: 10.1007/s11263-007-0075-7
2. Halko N, Martinsson PG, Tropp JA. Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. *SIAM Rev*. 2011;53(2):217-288. doi: 10.1137/090771806
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
- scikit-learn 1.9.1 (github.com/scikit-learn/scikit-learn).
<!-- REFERENCES:END -->

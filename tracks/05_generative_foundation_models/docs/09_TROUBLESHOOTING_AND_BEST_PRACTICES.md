# Chapter 09. 트러블슈팅, 벤치마킹 및 캡스톤 프로젝트 루브릭 (Troubleshooting, scIB Benchmarking & Capstone Project Rubric)

> **Cognitive Goal**: 대규모 단일 세포 트랜스포머 및 생성형 파운데이션 모델(scGPT, Geneformer) 훈련 및 추론 시 필연적으로 직면하는 GPU 메모리 한계(Out-of-Memory, OOM)와 연산 병목을 극복하는 하드웨어 최적화 엔지니어링 기법(Gradient Accumulation, Mixed Precision FP16/BF16, FlashAttention-2, Dynamic Chunked Tokenization)을 체득한다. 파운데이션 모델의 잠재 공간 왜곡 현상(Transformer Hallucinations, Batch Effect Leakage, Representation Oversmoothing)의 근본 원인을 진단하고 체계적인 디버깅 솔루션을 확립한다. 단일 세포 통합 벤치마킹 표준인 `scIB`의 4대 핵심 정량 지표($\text{cLISI}, \text{iLISI}, \text{NMI}, \text{ASW}$)의 수학적 유도 과정을 이해하고, 배치 제거와 생물학적 보존성 간의 파레토 프론티어(Pareto Frontier)를 최적화한다. 최종적으로 Track 05 전체를 완성하는 5단계 딥 생성 및 파운데이션 모델 캡스톤 프로젝트의 수행 지침과 다면 평가 루브릭(Rubric)을 완성한다.

---

## 1. GPU 메모리 병목 및 딥러닝 하드웨어 가속 엔지니어링

단일 세포 트랜스포머는 입력 시퀀스 길이 $L$(보통 $1,200 \sim 2,048$개 유전자 토큰)과 다층 인코더 구조로 인해 방대한 GPU 메모리를 소비합니다. 특히 고전적 멀티헤드 셀프 어텐션의 연산 복잡도와 활성화 텐서(Activation Tensor) 메모리는 $\mathcal{O}(L^2)$로 급격히 증가합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               단일 세포 파운데이션 모델을 위한 4대 GPU 메모리 최적화 기법              │
│                                                                                        │
│   1. Mixed Precision (AMP - FP16 / BF16)                                               │
│      - 활성화를 16-bit로 계산 ──► 메모리 절감, Tensor Core 가속 (가중치 FP32 유지)     │
│      - 언더플로우 방지를 위한 동적 손실 스케일링(Dynamic Loss Scaling) 적용            │
│                                                                                        │
│   2. 그래디언트 누적 (Gradient Accumulation)                                           │
│      - 물리적 소형 미니배치(예: batch_size=8)로 순전파 후 그래디언트 누적              │
│      - k 스텝마다 optimizer.step() 실행 ──► 거대 유효 배치 크기(Effective Batch Size)│
│                                                                                        │
│   3. FlashAttention-2 (IO-Aware Exact Attention)                                       │
│      - GPU SRAM과 HBM(High Bandwidth Memory) 간의 메모리 읽기/쓰기 최소화 (Tiling)     │
│      - O(L^2) 중간 어텐션 행렬을 HBM에 쓰지 않고 온라인 소프트맥스로 실시간 계산       │
│                                                                                        │
│   4. 동적 버킷팅 및 청크 토큰화 (Dynamic Bucketing & Chunking)                         │
│      - 세포별 실제 발현 유전자 수에 따라 비슷한 길이끼리 배치 구성 (Padding 최소화)    │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 혼합 정밀도 훈련 (Automatic Mixed Precision, AMP)
표준 32비트 단정밀도(FP32) 대신 16비트 반정밀도(FP16 또는 BF16)로 순전파 연산을 수행합니다(autocast는 마스터 가중치를 FP32로 유지하므로 절감되는 것은 주로 활성화 메모리이며, 절감 폭과 가속은 모델·하드웨어에 따라 다릅니다):
- **BF16 (Bfloat16)의 생물학적 우수성**:
  - FP16은 지수부(Exponent)가 5비트에 불과하여 발현량의 스케일 차이가 클 때 언더플로우(Underflow)나 오버플로우가 발생하기 쉽습니다.
  - 반면 BF16은 FP32와 동일한 8비트 지수부를 유지하므로, 별도의 손실 스케일링(Loss Scaling) 없이도 수렴 안정성을 보장하며 GPU 텐서 코어(Tensor Cores)에서 2~3배 빠른 연산 속도를 제공합니다.

### 1.2 그래디언트 누적 (Gradient Accumulation)의 수리적 등가성
목표 유효 배치 크기를 $B_{\text{target}} = K \times B_{\text{micro}}$라 할 때, 파라미터 업데이트는 다음과 같이 수행됩니다:

$$\mathbf{g}_{\text{acc}} = \frac{1}{K} \sum_{k=1}^K \nabla_{\theta} \mathcal{L}(\mathbf{x}^{(k)}; \theta)$$

$$\theta \leftarrow \theta - \eta \cdot \text{Adam}(\mathbf{g}_{\text{acc}})$$

- 물리적인 작은 미니배치 $B_{\text{micro}}$만을 GPU 메모리에 올려 활성화 텐서 크기를 최소화하면서도, 수학적으로는 거대 배치 $B_{\text{target}}$를 한 번에 연산한 것과 완전히 동일한 그래디언트 기댓값을 달성합니다.

### 1.3 FlashAttention-2: 메모리 계층 구조 기반 최적화
Dao (2023)의 FlashAttention-2는 어텐션 연산 시 거대한 $L \times L$ 어텐션 가중치 행렬을 느린 GPU HBM에 기록하지 않고, 빠른 온칩 SRAM 내부에서 블록 단위로 타일링(Tiling) 및 점진적 정규화(Online Softmax)를 수행합니다.
- 메모리 사용량: $\mathcal{O}(L^2) \to \mathcal{O}(L)$ (시퀀스 길이에 대해 선형 복잡도 달성)
- 연산 속도: FlashAttention(v1) 논문은 최적화된 베이스라인 대비 $2 \sim 4$배를 보고했고, FlashAttention-2는 v1 대비 약 2배 추가 가속을 보고했다(하드웨어·시퀀스 길이에 따라 다름)

---

## 2. 파운데이션 모델의 잠재 공간 왜곡 및 환각 디버깅

대규모 모델을 전이학습하거나 제로샷 추론을 수행할 때 자주 마주치는 3대 대표적 실패 양상과 이에 대한 체계적인 디버깅 가이드라인입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                      단일 세포 파운데이션 모델 3대 실패 양상 및 처방                   │
│                                                                                        │
│   [증상 1: 잠재 공간 과도한 스무딩 (Representation Oversmoothing)]                     │
│   - 현상: UMAP 상에서 선암종, 기저세포, 신경내분비세포가 하나의 둥근 솜사탕처럼 뭉침   │
│   - 원인: 트랜스포머 레이어가 너무 깊어 마르코프 정상 분포로 수렴, 가중치 감쇠 과다   │
│   - 처방: Residual connection 가중치 강화, Layer 수 축소, 정규화 위치 검증            │
│                                                                                        │
│   [증상 2: 배치 효과 누수 (Batch Effect Leakage)]                                      │
│   - 현상: 환자 번호나 시퀀싱 기술(10x v2 vs v3)에 따라 세포들이 군집 분리됨           │
│   - 원인: 절대 카운트 수치가 임베딩에 노출되었거나, 배치 토큰 정규화 실패              │
│   - 처방: Normalized Rank-value Encoding으로 전환, DANN 그래디언트 반전 계층 적용     │
│                                                                                        │
│   [증상 3: 비생물학적 전사체 환각 (Biological Hallucination)]                          │
│   - 현상: GEP 헤드가 미발현 유전자(0 count)에 대해 엉뚱하게 높은 발현량을 임의 생성   │
│   - 원인: Zero-inflation 정규화 결핍, GEP 손실 가중치(lambda_GEP)의 과도한 설정        │
│   - 처방: 비영 위치 마스크 손실 또는 Huber Loss 도입 (ZINB는 raw count용)             │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. 정량적 평가 메트릭: scIB (Single-Cell Integration Benchmarking) 프레임워크

모델이 배치 효과를 성공적으로 제거하면서도 실제 세포의 생물학적 정체성을 얼마나 잘 보존하였는지를 정량화하기 위해, Luecken et al. (Nature Methods 2022)의 **scIB**의 지표(공식 scIB는 생물 보존 9개 + 배치 제거 5개, 총 14개 지표)를 산출합니다. 아래에서는 그중 대표 지표 몇 개만 골라 교육용 간이 점수로 계산합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                         scIB 벤치마킹 지표 체계 및 최적화 목표                         │
│                                                                                        │
│   1. 생물학적 보존성 (Biological Conservation): cLISI, NMI, ASW_celltype               │
│      - 목표: 동일 세포형은 잠재 공간에서 밀집되고, 상이한 계통 간 경계는 선명히 유지 │
│                                                                                        │
│   2. 배치 혼합도 (Batch Removal / Correction): iLISI, ASW_batch                       │
│      - 목표: 서로 다른 환자 및 플랫폼에서 온 동일 세포들이 완벽하게 뒤섞여 통합       │
│                                                                                        │
│   [종합 평가 점수 (Overall Integration Score)]                                         │
│   Total Score = 0.6 * Bio Conservation Score + 0.4 * Batch Removal Score               │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 3.1 국소 역 심슨 지수 (Local Inverse Simpson's Index, LISI)의 수리적 정의
세포 $i$의 잠재 공간 좌표 $\mathbf{z}_i$ 주변의 $k$-최근접 이웃(k-NN, 보통 $k=90$)에 대해 가우시안 커널 기반 전이 확률을 정의합니다:

$$p_{ij} = \frac{\exp\left( - \frac{\|\mathbf{z}_i - \mathbf{z}_j\|_2^2}{2\sigma_i^2} \right)}{\sum_{m \neq i} \exp\left( - \frac{\|\mathbf{z}_i - \mathbf{z}_m\|_2^2}{2\sigma_i^2} \right)}$$

주변 이웃 내에서 특정 범주형 라벨 $c$ (배치 또는 세포형)의 확률은 $P_{i}(c) = \sum_{j \in \text{label } c} p_{ij}$입니다.
세포 $i$의 국소 역 심슨 지수 $D_i$는 다음과 같습니다:

$$D_i = \frac{1}{\sum_{c} P_i(c)^2}$$

1. **$\text{iLISI}$ (Integration LISI - 배치 혼합도)**:
   - 라벨 $c$를 '배치(Batch ID)'로 설정하여 계산합니다.
   - 배치 수가 $B$개일 때, $1 \le \text{iLISI} \le B$의 범위를 갖습니다.
   - 모든 배치가 국소 이웃 내에 균일하게 뒤섞여 있을 때 최댓값 $B$에 도달하며, 높을수록 우수한 배치 제거를 의미합니다.
2. **$\text{cLISI}$ (Cell-type LISI - 세포형 순도)**:
   - 라벨 $c$를 '세포형(Cell Type)'으로 설정하여 계산합니다.
   - 이상적인 경우 각 세포의 국소 이웃은 오직 자신과 동일한 세포형만으로 둘러싸여 있어야 하므로, $\text{cLISI} \to 1.0$에 수렴해야 합니다.

### 3.2 정규화 상호정보량 (Normalized Mutual Information, NMI)
잠재 공간의 비지도 클러스터링 라벨 $U$와 실제 생물학적 세포형 라벨 $V$ 사이의 정보 공유량을 측정합니다:

$$\text{NMI}(U, V) = \frac{2 I(U; V)}{H(U) + H(V)}$$

여기서 $I(U; V)$는 상호정보량(Mutual Information), $H(\cdot)$는 섀넌 엔트로피(Shannon Entropy)입니다. $0 \le \text{NMI} \le 1$ 범위를 가지며, 1에 가까울수록 클러스터링이 세포형을 완벽히 분리함을 나타냅니다.

### 3.3 평균 실루엣 폭 (Average Silhouette Width, ASW)
세포 $i$에 대해 같은 세포형 내 다른 세포들과의 평균 거리 $a(i)$와 가장 가까운 타 세포형 군집과의 평균 거리 $b(i)$를 측정합니다:

$$s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}, \quad s(i) \in [-1, 1]$$

- **$\text{ASW}_{\text{bio}}$**: 세포형 라벨에 대해 계산한 실루엣 점수의 평균을 $[0, 1]$로 재스케일링: $\frac{\bar{s} + 1}{2}$.
- **$\text{ASW}_{\text{batch}}$**: 배치 라벨에 대해 계산하며, 배치가 뒤섞여 실루엣이 0에 가까울수록 높은 점수를 부여: $1 - |\bar{s}|$.

---

## 4. 실전 Python 코드: GPU 하드웨어 가속 훈련 루프 및 scIB 메트릭 자동 산출

PyTorch의 AMP(Automatic Mixed Precision) 컨텍스트를 활용한 고속 훈련 루프와, 추출된 임베딩에 대한 간이 scIB 지표 산출 파이프라인입니다.

```python
import torch
import torch.nn as nn
from torch.cuda.amp import autocast, GradScaler
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score, silhouette_score
from typing import Dict, Tuple

def train_single_cell_transformer_amp(
    model: nn.Module,
    dataloader,
    optimizer: torch.optim.Optimizer,
    num_epochs: int = 5,
    gradient_accumulation_steps: int = 4,
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
) -> nn.Module:
    """
    Automatic Mixed Precision (AMP) 및 Gradient Accumulation이 적용된 트랜스포머 훈련 루프.
    """
    model.to(device)
    model.train()
    
    # FP16/BF16 그래디언트 스케일러 초기화
    # 주의: GradScaler는 FP16 용이다. BF16 autocast를 쓸 때는 필요 없다.
    scaler = GradScaler(enabled=(device == "cuda"))
    criterion = nn.CrossEntropyLoss()

    for epoch in range(num_epochs):
        running_loss = 0.0
        optimizer.zero_grad()

        for step, (tokens, values, labels) in enumerate(dataloader):
            tokens = tokens.to(device)
            values = values.to(device)
            labels = labels.to(device)

            # 1. AMP Context Manager를 통한 반정밀도 순전파
            with autocast(enabled=(device == "cuda")):
                _, logits = model(tokens, values)
                loss = criterion(logits, labels)
                # 그래디언트 누적 스텝으로 나눔
                loss = loss / gradient_accumulation_steps

            # 2. Scaled Backward Pass
            scaler.scale(loss).backward()
            running_loss += loss.item() * gradient_accumulation_steps

            # 3. 누적 스텝 도달 시 파라미터 업데이트
            if (step + 1) % gradient_accumulation_steps == 0 or (step + 1) == len(dataloader):
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()

        avg_loss = running_loss / len(dataloader)
        print(f"Epoch [{epoch+1}/{num_epochs}] - Mean Loss: {avg_loss:.4f}")

    return model


def compute_simplified_scib_metrics(
    embeddings: np.ndarray,
    cell_types: np.ndarray,
    batches: np.ndarray,
    k_neighbors: int = 30
) -> Dict[str, float]:
    """
    추출된 세포 임베딩에 대해 NMI, ASW, 간이 iLISI, cLISI 지표를 계산합니다.
    교육용 간이 점수이며 scIB 공식 정규화와 다르다. 난수 임베딩도 중간 정도의 종합 점수(아래 __main__ 난수 입력에서 약 0.6)를 받으므로 실제 평가에서는 scib-metrics를 쓰고 난수/PCA 베이스라인과 비교한다.
    """
    n_samples = embeddings.shape[0]
    
    # 1. ASW (Silhouette Score) 계산
    asw_celltype = (silhouette_score(embeddings, cell_types) + 1.0) / 2.0
    asw_batch = 1.0 - abs(silhouette_score(embeddings, batches))

    # 2. k-NN 기반 LISI 근사 계산
    nn_model = NearestNeighbors(n_neighbors=k_neighbors).fit(embeddings)
    _, indices = nn_model.kneighbors(embeddings)

    ilisi_scores = []
    clisi_scores = []

    for i in range(n_samples):
        nbrs = indices[i]
        
        # Batch Simpson's Index
        batch_counts = pd.Series(batches[nbrs]).value_counts(normalize=True).values
        simpson_batch = (batch_counts ** 2).sum()
        ilisi_scores.append(1.0 / (simpson_batch + 1e-8))

        # Cell-type Simpson's Index
        celltype_counts = pd.Series(cell_types[nbrs]).value_counts(normalize=True).values
        simpson_celltype = (celltype_counts ** 2).sum()
        clisi_scores.append(1.0 / (simpson_celltype + 1e-8))

    mean_ilisi = float(np.mean(ilisi_scores))
    mean_clisi = float(np.mean(clisi_scores))

    # 2.5 NMI: 임베딩을 KMeans(k=세포형 수)로 군집화한 결과와 세포형 라벨 간 NMI
    km = KMeans(n_clusters=len(np.unique(cell_types)), n_init=10, random_state=0).fit_predict(embeddings)
    nmi = normalized_mutual_info_score(cell_types, km)

    # 3. 종합 점수 가중 합산
    bio_score = (asw_celltype + nmi + (1.0 / mean_clisi)) / 3.0 # clisi는 1에 가까울수록 완벽
    batch_score = 0.5 * asw_batch + 0.5 * (mean_ilisi / len(np.unique(batches)))
    total_score = 0.6 * bio_score + 0.4 * batch_score

    return {
        "ASW_CellType (Bio)": asw_celltype,
        "ASW_Batch (Batch)": asw_batch,
        "mean_iLISI (Batch)": mean_ilisi,
        "mean_cLISI (Bio)": mean_clisi,
        "NMI_CellType (Bio)": nmi,
        "Overall_Score": total_score
    }


# 테스트 및 검증 실행
if __name__ == "__main__":
    np.random.seed(42)
    N = 200
    D = 64
    
    # 모의 임베딩, 세포형(2종: Adeno vs NEPC), 배치(2개 병원)
    mock_emb = np.random.randn(N, D)
    mock_celltypes = np.array(["AR-PC"] * 100 + ["t-NEPC"] * 100)
    mock_batches = np.array(["Hospital_A", "Hospital_B"] * 100)

    metrics = compute_simplified_scib_metrics(mock_emb, mock_celltypes, mock_batches)
    print("\n[단일 세포 파운데이션 모델 scIB 벤치마킹 산출 결과]")
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")
```

---

## 5. 5단계 딥 생성 및 파운데이션 모델 캡스톤 프로젝트 루브릭 (Capstone Project Rubric)

Track 05 전체(Chapter 00~09)에서 습득한 이론과 실습을 종합하여, 다기관 전립선암 단일 세포 전사체 아틀라스를 분석하고 신약 표적을 도출하는 **5단계 종합 캡스톤 프로젝트** 수행 지침 및 채점 기준입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   5단계 딥 생성 & 파운데이션 모델 캡스톤 파이프라인                     │
│                                                                                        │
│   [Stage 1: 데이터 조화 및 토큰화] ──► 10x v3 & Smart-seq2 전처리 및 Rank 인코딩      │
│   [Stage 2: 생성 대조군 모델링]   ──► scVI / scANVI 구축 및 잠재 공간 Baseline 확보    │
│   [Stage 3: 파운데이션 전이학습]   ──► scGPT / Geneformer Fine-tuning & 임베딩 추출      │
│   [Stage 4: In Silico 섭동 시뮬레이션] ──► EZH2/ASCL1 Virtual KO & 계통 역전 스크리닝 │
│   [Stage 5: scIB 정량 평가 및 중개] ──► cLISI/iLISI 평가 및 최종 연구 보고서 완성     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 캡스톤 다면 평가 기준표 (Capstone Project Rubric)

| 평가 영역 | 미흡 (Unsatisfactory, <60점) | 보통 (Proficient, 60~85점) | 탁월 (Exemplary, 85~100점) |
| :--- | :--- | :--- | :--- |
| **Stage 1: 토큰화 및 데이터 전처리** | 결측치 및 기술적 드롭아웃 미처리, 유전자 심볼 어휘집 불일치 | Log-CP10k 정규화 수행 및 기본 유전자 순위 정렬 완료 | 코퍼스 기준 스케일링을 포함한 Rank-Value 인코딩 구현 및 배치 편향 사전 차단 |
| **Stage 2: 딥 생성 모델 (scVI)** | 잠재 변수 차원 선정 근거 부재, ELBO 손실 발산 | scVI 기본 파라미터 훈련 및 재구성 오차 모니터링 성공 | 우도 선택(`gene_likelihood` zinb/nb) 근거 제시 및 분산 파라미터 최적화 및 계통 궤적(Latent Trajectory) 시각화 |
| **Stage 3: 트랜스포머 미세조정** | GPU OOM 에러 미해결, 단순 출력 레이어만 임의 학습 | AMP 및 Gradient Accumulation 적용, 미세조정 완료 | Pre-LN 및 동분산 불확실성 손실 가중치 최적화, 제로샷 임베딩 추출 |
| **Stage 4: In Silico 유전자 섭동** | 단순 무작위 섭동, 변위 벡터 수학적 정의 오류 | Virtual KO 수행 및 목표 표현형과의 유클리드 거리 측정 | Value Zeroing/Token Deletion 비교, 치환 검정 기반 $p$-value 검증 |
| **Stage 5: scIB 벤치마킹 및 보고서** | 정량 지표 미산출, 정성적 UMAP 시각화에만 의존 | ASW, NMI 산출 및 배치 제거 효과 확인 | cLISI/iLISI 파레토 최적화 증명 및 생물학적 기전 임상 보고서 완성 |

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [개념 및 메커니즘 (Conceptual Rigor)]
FlashAttention-2가 고전적 표준 셀프 어텐션(Standard Scaled Dot-Product Attention)과 수리적으로 동일한 출력(근사 없음)을 보장하면서도 연산 속도를 크게 가속하고 메모리 복잡도를 $\mathcal{O}(L^2)$에서 $\mathcal{O}(L)$로 절감할 수 있는 핵심 알고리즘 원리(Tiling 및 Online Softmax)를 GPU 메모리 계층 구조(HBM vs SRAM) 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **GPU 메모리 병목의 근본 원인**:
  - 표준 어텐션은 $S = \frac{QK^T}{\sqrt{d_k}} \in \mathbb{R}^{L \times L}$, $P = \text{softmax}(S) \in \mathbb{R}^{L \times L}$, $O = PV \in \mathbb{R}^{L \times d}$를 순차 계산합니다.
  - 이때 크기가 $L \times L$인 거대한 중간 행렬 $S$와 $P$를 느린 GPU 글로벌 메모리(High Bandwidth Memory, HBM)에 쓰고 다시 읽어오는 반복적인 I/O(Memory Bandwidth Bottleneck)가 전체 실행 시간의 상당 부분을 점유합니다(점유 비율은 시퀀스 길이·하드웨어에 따라 다르다).
- **FlashAttention-2의 2대 핵심 원리**:
  1. **타일링 (Tiling / Block Decomposition)**:
     - 쿼리 $Q$, 키 $K$, 밸류 $V$ 행렬을 프로세서 온칩의 초고속 SRAM(Static RAM) 크기에 맞추어 작은 블록($B_r \times d$, $B_c \times d$)으로 분할합니다.
     - 입력 블록들을 SRAM에 올린 뒤 연산을 수행하여 HBM과의 왕복 통신 횟수를 획기적으로 줄입니다.
  2. **온라인 소프트맥스 (Online Softmax Scaling)**:
     - 표준 소프트맥스는 분모 $\sum_{j=1}^L e^{S_{ij}}$를 계산하기 위해 행 전체($L$개 원소)를 한 번에 조회해야 합니다.
     - FlashAttention은 점화식을 활용하여 새로운 블록이 들어올 때마다 이전 블록의 최댓값 $m_{\text{new}} = \max(m_{\text{old}}, \max(S_{\text{block}}))$과 누적 합 $l_{\text{new}} = l_{\text{old}} e^{m_{\text{old}} - m_{\text{new}}} + \sum e^{S_{\text{block}} - m_{\text{new}}}$을 점진적으로 스케일링하며 갱신합니다.
     - 따라서 거대한 $L \times L$ 어텐션 가중치 행렬을 HBM에 단 한 번도 저장하지 않고, SRAM 내부에서 블록별로 계산된 중간 가중치를 축적하여 최종 출력 $O$를 직접 산출함으로써 메모리 복잡도를 $\mathcal{O}(L)$로 축소시킵니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
Luecken et al. (scIB)의 국소 역 심슨 지수(Local Inverse Simpson's Index)에서:
1) $B$개의 독립적인 배치로 구성된 데이터셋에서, 완벽한 배치 혼합이 이루어졌을 때 $\text{iLISI}$의 이론적 최댓값이 $B$가 됨을 역 심슨 공식 $D = \frac{1}{\sum_{b=1}^B P(b)^2}$을 통해 수리적으로 유도하시오.
2) 한 분석가가 배치 통합 손실(Batch Adversarial Loss)의 가중치를 지나치게 높여 훈련한 결과, $\text{iLISI}$는 이론적 최댓값인 $B$에 도달하였으나 $\text{cLISI}$가 $1.0$에서 $4.5$로 폭증하였다. 이 현상의 기하학적 의미를 '과보정(Over-correction)' 및 세포형 경계 붕괴 관점에서 수학적으로 해석하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. $\text{iLISI} = B$의 이론적 최댓값 유도**:
  - $B$개의 배치에 대해 세포 $i$의 국소 이웃 내 배치 점유 확률을 $P(b)$라 하자 ($\sum_{b=1}^B P(b) = 1$, $P(b) \ge 0$).
  - 코시-슈바르츠 부등식(Cauchy-Schwarz Inequality) 또는 옌센 부등식(Jensen's Inequality)에 의해, 제곱합 $\sum_{b=1}^B P(b)^2$은 모든 배치의 점유 확률이 완벽하게 균등할 때($P(1) = P(2) = \dots = P(B) = \frac{1}{B}$) 최솟값을 갖습니다.
  - 최솟값 계산:
    $$\min \sum_{b=1}^B P(b)^2 = \sum_{b=1}^B \left( \frac{1}{B} \right)^2 = B \times \frac{1}{B^2} = \frac{1}{B}$$
  - 따라서 역 심슨 지수 $D$의 최댓값은 다음과 같습니다:
    $$D_{\max} = \frac{1}{\min \sum_{b=1}^B P(b)^2} = \frac{1}{1/B} = B$$
  - 즉, 국소 이웃 내에 $B$개의 배치가 동일한 비율로 완벽히 융합되어 있을 때 $\text{iLISI} = B$가 성립합니다.

- **2. $\text{cLISI} = 4.5$ 폭증 현상의 수학적/기하학적 해석 (과보정 및 경계 붕괴)**:
  - $\text{cLISI}$의 이상적인 값은 $1.0$입니다. 즉, 국소 이웃 내에 단 하나의 세포형만 존재해야 생물학적 계통 순도가 유지됩니다.
  - $\text{cLISI} = 4.5$라는 것은 임의의 한 세포 주변 국소 이웃에 평균 $4.5$개의 서로 다른 세포형(예: 관상 선암종, 기저세포, 신경내분비세포, 면역세포)이 무차별적으로 뒤섞여 있음을 의미합니다.
  - **기하학적 원인 (Over-correction)**:
    - 모델이 기술적 배치를 지우는 데만 급급하여 세포형 간의 고유한 전사체 분별 특성(Biological Variation)까지 기술적 노이즈로 오인하여 강제로 압축·소거해 버렸습니다.
    - 그 결과 고차원 매니폴드 상에서 서로 격리되어 있어야 할 세포형 클러스터들이 단일 지점으로 붕괴·융합되어, 생물학적 마커 해석 및 차등 유전자 분석이 완전히 불가능해진 치명적인 '생물학적 정보 소실 상태'를 수학적으로 증명합니다.

</details>

---

### Q3. [임상 및 생물학적 케이스 스터디 (Clinical Problem Solving)]
3개 대학병원(Site 1: 10x v3, Site 2: 10x v2, Site 3: Smart-seq2)에서 수집된 mCRPC 환자 120명의 단일 세포 전사체 아틀라스를 통합 훈련하는 과정에서, 검증 단계 손실(Validation Loss)이 3 에포크 이후 급격히 발산하며 모델의 GEP 헤드가 모든 유전자에 대해 $0.001$ 이하의 값만을 출력하는 '자명한 붕괴(Trivial Zero Collapse)' 현상이 발생하였다.
1) 손실 함수 발산 및 자명한 붕괴를 초래한 3대 핵심 원인(학습률 스케줄러, 단일 세포 스파시티, 손실 함수 스케일 불균형)을 진단하시오.
2) 이 파이프라인을 복구하여 성공적인 전이학습 모델로 전환하기 위해 취해야 할 3단계 엔지니어링 조치(Warmup Scheduler, GEP 손실 교체, Layer-wise Differential LR)를 구체적으로 처방하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 발산 및 자명한 붕괴의 3대 핵심 원인 진단**:
  - **학습률 스케줄러 미비**: 트랜스포머 초기화 직후 큰 학습률이 적용되면서 어텐션 가중치가 극단값으로 튀어 소프트맥스 그래디언트가 폭발(Gradient Explosion).
  - **극단적 스파시티(Zero-inflation)에 의한 MSE 왜곡**: 단일 세포 데이터의 $90\%$ 이상이 0 카운트이므로, 단순 MSE 손실 하에서 모델이 모든 출력을 0으로 예측할 때 전체 오차의 $90\%$가 손쉽게 줄어드는 거짓 국소 최적점(Trivial Zero Trap)에 포획됨.
  - **다중 손실 함수 간 척도 불균형**: 분류 손실(교차 엔트로피, 값 $2 \sim 5$ 가정) 대비 비정규화된 GEP(MSE, 데이터 척도에 따라 큰 값이 될 수 있다고 가정; 예시)의 그래디언트 크기가 압도하여 파라미터를 파괴적으로 업데이트함.

- **2. 파이프라인 복구를 위한 3단계 엔지니어링 처방**:
  - **조치 1: 선형 웜업(Linear Warmup) 및 코사인 어닐링 스케줄러 적용**:
    - 초기 5 에포크 동안 학습률을 $0$에서 $10^{-4}$까지 선형적으로 서서히 증가시켜 트랜스포머 어텐션 층이 안정적인 기저 기하 구조를 형성하도록 유도한 후, 코사인 곡선을 따라 감쇄시킵니다.
  - **조치 2: GEP 손실 함수를 Masked Normalized Huber Loss 등으로 교체**(scGPT의 GEP는 구간화·정규화된 값을 회귀하므로 raw count용 ZINB 교체는 호환되지 않는다):
    - 비발현된 0 값에 대해서는 직접적인 MSE를 부여하지 않고, 실제로 발현된 Non-zero 유전자 마스크 위치에서만 스케일링된 오차를 역전파하도록 제한합니다.
  - **조치 3: 계층별 차등 학습률(Layer-wise Differential LR) 및 사전학습 가중치 동결**:
    - 대규모 인간 코퍼스에서 기학습된 하위 트랜스포머 인코더 계층(예: 12층 중 1~8층)은 고정(Freeze, $\text{lr}=0$)하거나 극소 학습률($10^{-6}$)을 부여하고, 상위 계층 및 출력 헤드에만 상대적으로 높은 학습률($10^{-4}$)을 적용하여 기저의 전사체 문법 지식 붕괴(Catastrophic Forgetting)를 차단합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods*. 2022;19(1):41-50. PMID: 34949812. doi: 10.1038/s41592-021-01336-8
2. Dao T. FlashAttention-2: Faster Attention with Better Parallelism and Work Partitioning. *arXiv* preprint. 2023. arXiv: 2307.08691.
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- scikit-learn 1.9.1 (github.com/scikit-learn/scikit-learn).
- PyTorch 2.14.1 (github.com/pytorch/pytorch).
<!-- REFERENCES:END -->

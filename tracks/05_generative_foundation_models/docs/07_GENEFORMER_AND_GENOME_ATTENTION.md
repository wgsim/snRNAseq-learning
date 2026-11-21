# Chapter 07. Geneformer 아키텍처와 어텐션 기반 유전자 네트워크 분석 (Geneformer: Genome-Scale Pre-training & Attention-based Network Biology)

> **Cognitive Goal**: 약 3,000만 개(Genecorpus-30M, 악성 세포·불멸화 세포주 제외)의 인간 단일 세포 전사체 데이터로 사전학습된 유전체 파운데이션 모델 `Geneformer`(Theodoris et al., Nature 2023)의 6계층(V1, 10M 파라미터 기준) 트랜스포머 아키텍처와 핵심 생물정보학적 철학을 체득한다. 극단적인 시퀀싱 깊이 차이와 기술적 증폭 아티팩트를 비모수적으로 제거하는 정규화 순위 인코딩(Normalized Rank-Value Encoding)의 수학적 불변성(Invariance)을 증명한다. 정적 유전자 임베딩(Static Embedding)의 한계를 뛰어넘어 세포 상태와 문맥에 따라 동적으로 재구성되는 문맥 인식 유전자 표상(Context-Aware Gene Representations)을 규명한다. 다층 트랜스포머의 어텐션 행렬들을 집계하는 어텐션 롤아웃(Attention Rollout) 알고리즘을 유도하고, 어텐션 기반 중심성(Attention-based Network Centrality) 지표를 수학적으로 정의하여 전립선암 계통 가소성(t-NEPC)을 추동하는 핵심 마스터 조절자(Master Regulator: *ASCL1*, *EZH2*, *SOX2*)를 발굴하는 분석 파이프라인을 구축한다.

---

## 1. Geneformer 모델의 패러다임과 아키텍처 설계

생명체의 게놈은 모든 세포에서 거의 동일하지만, 각 세포가 어떤 표현형을 획득하는가는 유전자 조절 네트워크(Gene Regulatory Network, GRN)의 위상학적 배선(Wiring)에 의해 결정됩니다. 전통적인 네트워크 생물학(WGCNA, ARACNe)은 세포 집단 전체의 선형적 상관계수에 의존하므로 희귀 세포군이나 세포 상태별 동적 조절 관계를 포착하는 데 명확한 한계가 있었습니다.

Theodoris et al. (Nature 2023)이 발표한 **Geneformer**는 대규모 자기지도 사전학습을 통해 인간 전사체에 내재된 전사 조절 문법(Regulatory Grammar)을 스스로 학습한 유전체 파운데이션 모델입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        Geneformer 6-Layer 트랜스포머 아키텍처                           │
│                                                                                        │
│   [단일 세포 입력]                                                                     │
│   발현된 유전자들의 정규화 순위 정렬: [g_(1), g_(2), ..., g_(L)] (L <= 2048)            │
│               │                                                                        │
│               ▼ 유전자 토큰 임베딩 (25,424개 유전자 어휘집, d_model = 256)             │
│   입력 텐서 H_0: [L, 256]                                                              │
│               │                                                                        │
│   ┌───────────┴────────────────────────────────────────────────────────────────────┐   │
│   │   Transformer Layer 1 (Multi-Head Self-Attention, H=4, d_k=64, FFN=512)        │   │
│   │   Transformer Layer 2 (Residual Connection + LayerNorm)                        │   │
│   │   Transformer Layer 3                                                          │   │
│   │   Transformer Layer 4                                                          │   │
│   │   Transformer Layer 5                                                          │   │
│   │   Transformer Layer 6                                                          │   │
│   └───────────┬────────────────────────────────────────────────────────────────────┘   │
│               │                                                                        │
│       ┌───────┴────────────────────────────────────────┐                               │
│       ▼                                                ▼                               │
│   [문맥 인식 유전자 임베딩 H_6]                     [다층 어텐션 가중치 A^(1~6)]        │
│   - 특정 세포 문맥 하의 유전자 상태 표상             - Attention Rollout 집계           │
│   - In Silico Perturbation 시뮬레이션 기반          - 유전자 조절 네트워크 중심성 계산  │
│                                                     - 질환 마스터 조절자 발굴          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 하이퍼파라미터 및 모델 구조
- **사전학습 코퍼스**: 인간의 다양한 비암성(정상 및 비악성 질환) 조직에서 수집된 약 3,000만 개(30M) 단일 세포 전사체. 높은 돌연변이 부담의 악성 세포와 불멸화 세포주는 제외되었다.
- **아키텍처 구성**:
  - 인코더 전용 트랜스포머(BERT-style Encoder).
  - 레이어 수 ($L_{\text{layers}}$): 6개 계층.
  - 은닉 차원 ($d_{\text{model}}$): 256.
  - 어텐션 헤드 수 ($H$): 4개 헤드 ($d_k = d_{\text{model}} / H = 64$).
  - 피드포워드 차원 ($d_{\text{ff}}$): 512.
  - 최대 입력 시퀀스 길이: 2,048개 유전자 토큰.
  - 유전자 어휘집 크기 ($|\mathcal{V}|$): 25,424개의 Ensembl 유전자 ID.
  - 주의: 위 구성은 **V1(10M 파라미터)** 기준이다. 현재 HF 기본 모델은 V2(104M·316M; 입력 4,096, 어휘 약 20K, 코퍼스 약 1.04억 세포)이며, 암 도메인 튜닝 변형(약 1,400만 암 세포)도 별도로 있다. 사전학습 코퍼스가 악성 세포를 제외했으므로, t-NEPC 같은 종양 세포에 적용할 때는 분포 밖(OOD) 한계를 감안하고 미세조정·검증을 병행해야 한다.

---

## 2. 정규화 순위 인코딩(Normalized Rank-Value Encoding)의 수리적 엄밀성

Geneformer의 가장 독창적인 핵심 메커니즘은 절대적인 카운트 수치를 입력 벡터에 더하지 않고, 유전자별 코퍼스 기준으로 정규화한 발현값의 **순위(Rank-Value)**를 통해 시퀀스를 구성한다는 점입니다(Geneformer 모델카드: 세포 내 발현을 코퍼스 전체의 해당 유전자 발현으로 스케일한 뒤 순위화하여, 어디서나 높게 발현되는 하우스키핑 유전자는 강등하고 세포 상태를 구분하는 유전자는 승격).

### 2.1 전통적 카운트 정규화의 한계
단일 세포 RNA-seq 데이터에서 측정된 원시 카운트(Raw Count) $C_{cg}$ (세포 $c$, 유전자 $g$)는 세포의 크기, 역전사 효율, PCR 증폭 편향, 시퀀싱 심도(Sequencing Depth) $N_c = \sum_{g} C_{cg}$에 의해 심각하게 왜곡됩니다.
일반적으로 사용하는 $\log(\text{CP10k} + 1)$ 정규화:

$$\tilde{x}_{cg} = \log\left( \frac{C_{cg}}{N_c} \times 10^4 + 1 \right)$$

는 세포 간 시퀀싱 심도 차이를 선형 보정하지만, 고발현 유전자의 비선형적 포화나 드롭아웃(Zero-inflation) 확률의 차이를 근본적으로 제거하지 못합니다.

### 2.2 정규화 순위 인코딩 알고리즘
세포 $c$에서 발현된(비영) 유전자 집합을 $\mathcal{G}_c = \{g \in \mathcal{V} \mid C_{cg} > 0\}$이라 하고, 검출된 총 유전자 수를 $L_c = |\mathcal{G}_c|$라 합시다.

1. **코퍼스 기준 정규화**: 세포 $c$의 발현량을 라이브러리 크기로 정규화한 뒤, 유전자 $g$의 코퍼스 비-영(非零) 발현 중앙값 $m_g$로 나눕니다: $\hat{x}_{cg} = \dfrac{C_{cg}/N_c}{m_g}$. (이 단계가 없으면 하우스키핑 유전자가 항상 상위에 놓입니다.)
2. **순위 매김**: $\hat{x}_{cg}$를 기준으로 내림차순 정렬하여 각 유전자에 순위 $\text{Rank}_c(g) \in \{1, 2, \dots, L_c\}$를 부여합니다 (가장 높은 유전자가 1위). 순위 자체(정수)를 토큰 위치로 사용하며, 순위를 $L_c$로 나눈 백분위수 $r_c = \text{Rank}/L_c$는 설명용 지표일 뿐 실제 입력에 쓰이지 않습니다.
3. **토큰 시퀀스 구성**: 가장 지배적으로 발현하는 상위 유전자부터 차례대로 나열하여 시퀀스 $\mathbf{s}_c = [g_{(1)}, g_{(2)}, \dots, g_{(L_c)}]$를 구성합니다 ($L_c > 2048$인 경우 상위 2,048개에서 절단).

### 2.3 전역 스케일링에 대한 순위 불변성 (제한된 성질)
**정리 (Scale Invariance of Rank Encoding)**:
세포 $c$의 **모든 유전자에 동일하게** 적용되는 양의 단조 증가 변환 $f: \mathbb{R}^+ \to \mathbb{R}^+$ (예: 전역 상수 배 $f(x) = \alpha x$, $\alpha > 0$) 하에서, 유전자 시퀀스 $\mathbf{s}_c$와 순위는 불변이다. (아래 증명은 이 자명한 성질에 한정되며, 실제 깊이 차이에는 성립하지 않는다.)

**증명**:
임의의 두 발현 유전자 $g_a, g_b \in \mathcal{G}_c$에 대해, $C_{c, g_a} > C_{c, g_b}$라 하자.
단조 증가 함수의 정의에 의해:
$$f(C_{c, g_a}) > f(C_{c, g_b}) \quad (\forall \alpha > 0)$$
따라서 모든 유전자 쌍의 대소 관계가 보존되므로:
$$\text{Rank}_{f(\mathbf{C}_c)}(g) = \text{Rank}_{\mathbf{C}_c}(g) \quad (\forall g \in \mathcal{G}_c)$$
정규화 순위 분모 $L_c$ 역시 $f(C) > 0 \iff C > 0$이므로 불변이다.
$$\therefore r_c^{(f)}(g) = \frac{\text{Rank}_{f(\mathbf{C}_c)}(g)}{L_c} = r_c(g)$$
즉, 전역 상수 배 스케일링에는 입력 시퀀스 $\mathbf{s}_c$가 같습니다. $\blacksquare$

**한계**: 실제 시퀀싱 깊이 차이는 단순한 전역 상수 배가 아닙니다. 깊이가 낮으면 검출되는 유전자 집합 $\mathcal{G}_c$ 자체($L_c$)가 줄고 저발현 유전자의 순서가 샘플링 잡음으로 바뀌므로, 순위 인코딩이 깊이·배치 효과를 "완벽히" 제거한다고 할 수 없고 일부 완화할 뿐입니다.

---

## 3. 문맥 인식 유전자 표상 (Context-Aware Gene Representations)

기존 생물정보학에서 유전자를 임베딩하는 방식(예: Gene2Vec, STRING PPI 임베딩)은 정적 임베딩(Static Embedding)이었습니다. 즉, 유전자 *AR*이나 *TP53*은 항상 고정된 하나의 벡터 $\mathbf{v}_{\text{AR}} \in \mathbb{R}^d$로 표현되었습니다.

그러나 생물학적으로 유전자의 기능과 상호작용은 세포의 상태와 문맥에 따라 극단적으로 달라집니다:
- **루미날 전립선 상피세포 문맥**: *AR*은 *FOXA1*, *HOXB13*과 복합체를 이루어 *KLK3*(PSA), *TMPRSS2*를 전사 활성화하는 마스터 생존 인자로 작동합니다.
- **치료 유발 신경내분비암(t-NEPC) 문맥**: *TP53/RB1* 결손 하에서 *AR*은 *EZH2* 매개 후성유전학적 침묵에 의해 억제되며, 잔여 전사체 문맥에서도 정상적인 관상 전사 회로와의 상호작용이 단절됩니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   정적 임베딩 vs Geneformer 문맥 인식 유전자 임베딩                     │
│                                                                                        │
│   [정적 임베딩 (Static Embedding: Gene2Vec)]                                           │
│   AR Token ──► [고정 Lookup Table] ──► v_AR (모든 세포에서 동일한 벡터)               │
│                                                                                        │
│   [Geneformer 동적 문맥 임베딩 (Dynamic Contextualized Representation)]                │
│   Cell 1 (Adeno): [KLK3, FOXA1, AR, ...] ──► [Transformer 6-Layers] ──► h_AR^(Adeno)   │
│                                                                           ▲            │
│                                                                           │ d_cosine   │
│                                                                           ▼            │
│   Cell 2 (t-NEPC): [ASCL1, EZH2, AR, ...] ──► [Transformer 6-Layers] ──► h_AR^(t-NEPC) │
│                                                                                        │
│   * 결과: 동일한 AR 토큰이라도 세포 계통 문맥에 따라 잠재 공간 내 좌표가 동적으로 이동!  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

Geneformer의 6개 트랜스포머 인코더 계층을 통과한 유전자 $g_i$의 최종 은닉 벡터 $\mathbf{h}_i^{(6)} \in \mathbb{R}^{d_{\text{model}}}$는 해당 세포에 함께 공존하는 다른 모든 유전자들과의 셀프 어텐션 연산을 거쳤으므로, **세포 특이적 기능 상태가 반영된 문맥 인식 벡터(Contextualized Vector)**가 됩니다.

---

## 4. 어텐션 행렬 기반 유전자 조절 네트워크 및 중심성(Network Centrality) 유도

단일 세포 트랜스포머의 어텐션 행렬은 유전자 토큰 간 정보 전달의 한 표현이며, 이를 유전자 조절 네트워크 후보로 해석하려는 시도가 있습니다. 다만 어텐션 가중치는 인과·조절 관계 자체가 아니며(Geneformer 논문은 어텐션이 네트워크 위계를 인코딩한다고만 주장), 이하의 중심성 분석은 후보 탐색용이고 실험 검증이 필요합니다.

### 4.1 어텐션 롤아웃 (Attention Rollout) 알고리즘의 수리적 전개
트랜스포머의 각 계층 $l \in \{1, \dots, L\}$에서 헤드별 어텐션 행렬 $A^{(l, h)} \in \mathbb{R}^{K \times K}$ ($h \in \{1, \dots, H\}$)의 평균 어텐션 행렬을 $\bar{A}^{(l)}$이라 정의합시다:

$$\bar{A}^{(l)} = \frac{1}{H} \sum_{h=1}^H A^{(l, h)}$$

트랜스포머의 잔차 연결(Residual Connection: $\mathbf{x}^{(l)} = \mathbf{x}^{(l-1)} + \text{MHA}(\mathbf{x}^{(l-1)})$)을 고려하면, 입력 계층에서 다음 계층으로 전달되는 정보 흐름 행렬은 항등 행렬 $I$와 어텐션 행렬의 합으로 나타낼 수 있습니다:

$$\hat{A}^{(l)} = \frac{1}{2} I + \frac{1}{2} \bar{A}^{(l)} \quad \text{또는 가중치 정규화 형태 } \tilde{A}^{(l)} = I + \bar{A}^{(l)}$$

Abnar & Zuidema (ACL 2020)의 **어텐션 롤아웃(Attention Rollout)** 재귀 식에 따라, 계층 1부터 최종 계층 $L$까지 정보가 전파되는 누적 유전자 상호작용 행렬 $\mathcal{W}_{\text{rollout}} \in \mathbb{R}^{K \times K}$은 다음과 같은 행렬 곱으로 계산됩니다:

$$\mathcal{W}_{\text{rollout}} = \prod_{l=1}^L \left( \frac{I + \bar{A}^{(l)}}{2} \right) = \left( \frac{I + \bar{A}^{(L)}}{2} \right) \left( \frac{I + \bar{A}^{(L-1)}}{2} \right) \dots \left( \frac{I + \bar{A}^{(1)}}{2} \right)$$

- 행렬의 각 원소 $\mathcal{W}_{ij}$는 유전자 $i$의 최종 표현형에 유전자 $j$가 기여한 전체 다층 경로(Multi-layer Path)의 누적 주의력(Cumulative Attention Flow)을 나타냅니다.

### 4.2 유전자 네트워크 중심성(Centrality) 지표의 정의
누적 상호작용 행렬 $\mathcal{W}$를 가중 방향성 그래프 $G = (V, E)$의 인접 행렬(Adjacency Matrix)로 간주할 때, 각 유전자의 시스템적 중요도를 나타내는 중심성 지표를 수학적으로 정의합니다:

1. **아웃디그리 중심성 (Out-degree Centrality / Regulatory Influence)**:
   $$C_{\text{out}}(g_j) = \sum_{i=1}^K \mathcal{W}_{ij}$$
   - 유전자 $j$가 네트워크 내의 다른 모든 유전자들에게 미치는 영향력의 총합입니다.
   - 높은 $C_{\text{out}}$을 갖는 유전자는 전사 인자(Transcription Factor)나 크로마틴 리모델러와 같은 **마스터 조절자(Master Regulator)**에 해당합니다.

2. **행 합과 통합 엔트로피 (In-degree는 정보가 없음)**:
   $$C_{\text{in}}(g_i) = \sum_{j=1}^K \mathcal{W}_{ij} = 1 \quad (\forall i)$$
   - 롤아웃 행렬은 각 층 $(I+\bar{A}^{(l)})/2$(행 합이 1인 행-확률 행렬)의 곱이므로 **모든 행 합이 정확히 1**입니다. 따라서 행 합으로 정의한 인디그리 중심성은 모든 유전자에서 1.0이어서 하류 효과기를 선별하는 지표로 쓸 수 없습니다.
   - 대안으로 행의 분포가 얼마나 여러 유전자에 퍼져 있는지를 보는 **통합 엔트로피** $H_i = -\sum_j \mathcal{W}_{ij}\log \mathcal{W}_{ij}$를 쓸 수 있습니다(값이 클수록 많은 유전자의 정보를 통합).

3. **어텐션 고유벡터 중심성 (Attention Eigenvector Centrality)**:
   인접 행렬의 전치 $\mathcal{W}^T$의 최대 고유값 $\lambda_{\max}$에 대응하는 고유벡터 $\mathbf{v}$:
   $$\mathcal{W}^T \mathbf{v} = \lambda_{\max} \mathbf{v}$$
   - 단순히 많은 연결을 갖는 것에 그치지 않고, "영향력이 큰 다른 마스터 조절자들과 강력하게 결합되어 있는 핵심 허브 유전자"를 정밀하게 식별합니다.

---

## 5. 전립선암 계통 전이(t-NEPC) 마스터 조절자 발굴 사례

전립선암의 안드로겐 의존성 선암종(AR-PC)과 거세 저항성 신경내분비암(t-NEPC) 환자 시료에서 어텐션 중심성의 차이($\Delta C_{\text{out}}$)를 계산하면, 질환 진행에 관여하는 전사 인자 후보를 가려내는 데 쓰일 수 있습니다. 아래 표는 **실제 데이터로 산출한 결과가 아니라 개념 설명을 위한 가상 예시**입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│        전립선암 진행 단계별 어텐션 아웃디그리 중심성 변화 (가상 예시 수치)             │
│                                                                                        │
│   유전자 (Gene)       AR-PC 상태 (Adeno)          t-NEPC 상태 (Neuroendocrine)  ΔC_out │
│   ──────────────────────────────────────────────────────────────────────────────────   │
│   AR                  ████████████████ (0.88)     █░░░░░░░░░░░░░░░ (0.05)       -0.83  │
│   FOXA1               ██████████████░░ (0.75)     ███░░░░░░░░░░░░░ (0.18)       -0.57  │
│   ASCL1               ░░░░░░░░░░░░░░░░ (0.01)     ████████████████ (0.92)       +0.91  │
│   EZH2                ███░░░░░░░░░░░░░ (0.15)     ██████████████░░ (0.84)       +0.69  │
│   SOX2                ░░░░░░░░░░░░░░░░ (0.02)     ████████████░░░░ (0.72)       +0.70  │
│   INSM1               ░░░░░░░░░░░░░░░░ (0.00)     ███████████░░░░░ (0.68)       +0.68  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

- **해석**:
  - 치료 저항성 t-NEPC로 전환되면서, 세포의 전사 조절 중심축이 *AR* / *FOXA1*에서 *ASCL1* / *EZH2* / *SOX2*로 전면 교체(Paradigm Shift)됩니다.
  - (가상 예시이므로) 이 표 자체가 단순 발현량 차이(Differential Expression) 분석보다 민감하다는 근거가 될 수 없습니다. 실제로 어텐션 중심성이 DE보다 더 많은 정보를 주는지는 실제 모델·데이터에서 기준선과 비교해 검증해야 합니다.

---

## 6. 실전 Python 코드: Geneformer 스타일 어텐션 롤아웃 및 네트워크 중심성 산출

PyTorch로 구현된 다층 트랜스포머의 어텐션 텐서로부터 Attention Rollout을 계산하고 `networkx`를 통해 마스터 조절자 중심성을 시각화하는 파이프라인입니다.

```python
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Tuple

def compute_attention_rollout(all_layer_attentions: List[torch.Tensor]) -> np.ndarray:
    """
    모든 계층의 멀티헤드 어텐션 가중치로부터 Attention Rollout 행렬을 계산합니다.
    - all_layer_attentions: 리스트 길이 L, 각 원소의 형상 [Batch=1, Num_Heads, K, K]
    """
    num_layers = len(all_layer_attentions)
    K = all_layer_attentions[0].shape[-1]
    
    # 롤아웃 누적 행렬 초기화 (항등 행렬)
    rollout = torch.eye(K, dtype=torch.float32)

    for l in range(num_layers):
        # 헤드 평균 계산: [1, H, K, K] -> [K, K]
        attn_layer = all_layer_attentions[l].squeeze(0).mean(dim=0)
        
        # 잔차 연결(Residual Connection) 합성: 0.5 * I + 0.5 * A
        identity = torch.eye(K, dtype=torch.float32)
        attn_with_residual = 0.5 * identity + 0.5 * attn_layer
        
        # 행 정규화 (Row-normalization)
        attn_normalized = attn_with_residual / attn_with_residual.sum(dim=-1, keepdim=True)
        
        # 롤아웃 행렬 곱 갱신: Rollout = A_norm * Rollout
        rollout = torch.matmul(attn_normalized, rollout)

    return rollout.detach().cpu().numpy()


def calculate_network_centrality(rollout_matrix: np.ndarray, gene_names: List[str]) -> pd.DataFrame:
    """
    어텐션 롤아웃 행렬로부터 Out-degree, 통합 엔트로피, Eigenvector 중심성을 계산합니다.
    """
    K = len(gene_names)
    
    # 1. Out-degree (각 열의 합: 타 유전자에 미치는 영향력)
    out_degrees = rollout_matrix.sum(axis=0)
    
    # 2. 행 합은 롤아웃 행렬(행-확률 행렬의 곱)에서 항상 1이라 정보가 없으므로,
    #    통합 엔트로피(행 분포가 얼마나 퍼져 있는지)를 사용한다.
    P = np.clip(rollout_matrix, 1e-12, None)
    in_degrees = -(P * np.log(P)).sum(axis=1)

    # 3. NetworkX 방향성 가중치 그래프 생성
    G = nx.DiGraph()
    for i in range(K):
        G.add_node(gene_names[i])

    for i in range(K):
        for j in range(K):
            if i != j and rollout_matrix[i, j] > 1e-4:
                # 방향: j (조절자) -> i (표적)
                G.add_edge(gene_names[j], gene_names[i], weight=float(rollout_matrix[i, j]))

    # 고유벡터 중심성 계산
    try:
        eigenvector_dict = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
    except Exception:
        eigenvector_dict = {g: 0.0 for g in gene_names}

    df = pd.DataFrame({
        "Gene": gene_names,
        "Out_Degree": out_degrees,
        "Integration_Entropy": in_degrees,
        "Eigenvector_Centrality": [eigenvector_dict[g] for g in gene_names]
    }).sort_values(by="Out_Degree", ascending=False).reset_index(drop=True)

    return df


# 실행 및 파이프라인 무결성 검증
if __name__ == "__main__":
    # 전립선암 핵심 10대 유전자 정의
    target_genes = ["ASCL1", "EZH2", "SOX2", "INSM1", "AR", "KLK3", "FOXA1", "TP53", "RB1", "SYP"]
    K = len(target_genes)
    num_layers = 6
    num_heads = 4

    # 모의 멀티헤드 어텐션 행렬 생성. 주의: ASCL1·EZH2의 어텐션을 인위적으로 키워 넣은 장난감 입력이므로
    # 아래 결과(ASCL1·EZH2가 상위)는 입력 가정의 재현일 뿐 생물학적 발견이 아니다 (파이프라인 작동 확인용).
    torch.manual_seed(42)
    mock_attentions = []
    
    for l in range(num_layers):
        attn = torch.randn(1, num_heads, K, K).abs()
        # ASCL1(0번), EZH2(1번)에 높은 키 가중치 부여
        attn[:, :, :, 0] += 2.5
        attn[:, :, :, 1] += 2.0
        # 소프트맥스 정규화
        attn = torch.softmax(attn, dim=-1)
        mock_attentions.append(attn)

    # 1. Attention Rollout 계산
    rollout_mat = compute_attention_rollout(mock_attentions)
    print(f"계산된 Attention Rollout 행렬 크기: {rollout_mat.shape}")

    # 2. 중심성 지표 산출
    centrality_df = calculate_network_centrality(rollout_mat, target_genes)
    print("\n[Geneformer 기반 전립선암 조절 네트워크 중심성 분석 결과]")
    print(centrality_df.to_string(index=False))
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [개념 및 메커니즘 (Conceptual Rigor)]
Geneformer의 정규화 순위 인코딩(Normalized Rank-Value Encoding)이 시퀀싱 심도(Sequencing Depth)와 같은 기술적 배치 효과(Batch Artifacts)를 어느 범위까지 완화하는지 수학적 메커니즘을 설명하고, 반대로 이 인코딩 방식이 세포의 전사체 전체가 균일하게 억제되거나 활성화되는 전역적 전사체 증폭(Global Transcriptional Amplification, 예: *MYC* 유전체 증폭에 의한 전사 폭발) 현상을 포착할 때 가질 수 있는 이론적 한계점을 논하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **기술적 배치 효과 제거의 수학적 메커니즘 (순위 불변성)**:
  - 정규화 순위 인코딩은 유전자의 절대 카운트 수치를 폐기하고, 코퍼스 기준으로 정규화한 값 $\hat{x}_g$의 세포 내 내림차순 순위만을 입력으로 사용합니다.
  - 시퀀싱 심도의 차이는 모든 유전자 발현량에 대략 양의 상수 스케일러 $\alpha > 0$가 곱해지는 단조 증가 변환(Monotonic Transformation)으로 근사할 수 있습니다.
  - 순위 연산은 전역 상수 배에 대해 불변($\text{Rank}(\alpha \mathbf{x}) = \text{Rank}(\mathbf{x})$)이므로 이 부분의 깊이 차이는 완화됩니다. 그러나 실제 깊이 차이는 검출 유전자 집합과 저발현 유전자의 순서를 바꾸므로 동일 세포의 10,000 UMI와 50,000 UMI 프로파일이 정확히 같은 시퀀스를 만든다고 보장할 수 없고, 기술적 노이즈를 "근원적으로" 차단하지는 못합니다.
- **이론적 한계점: 전역적 전사 증폭(Global Amplification) 포착 불가**:
  - 암 생물학에서 *MYC* 종양유전자의 과발현은 특정 유전자 몇 개만 올리는 것이 아니라, 세포 내 RNA 중합효소 II(Pol II)를 전역적으로 자극하여 세포당 총 mRNA 분자 수를 2~3배 폭증시키는 '전사 증폭(Transcriptional Amplification)'을 일으킵니다.
  - 정규화 순위 인코딩은 전체 분자 수의 절대적 증가를 전혀 반영하지 못하며 오직 유전자들 간의 상대적 점유 비율만을 보존합니다.
  - 따라서 모든 유전자가 균일하게 3배씩 증가한 *MYC* 과발현 세포와 정상 세포가 Geneformer 입력단에서는 완전히 동일한 순위 시퀀스로 변환되어, 질환의 치명적인 전역적 전사체 과잉 상태를 감지하지 못하고 놓치는 맹점(Blind Spot)이 존재합니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
다층 트랜스포머에서 계층별 어텐션 행렬 $A^{(l)}$을 단순히 곱하지 않고 잔차 연결(Residual Connection)을 포함하여 $\tilde{A}^{(l)} = \frac{1}{2}(I + A^{(l)})$로 정의한 후 어텐션 롤아웃 $\mathcal{W} = \prod_{l=1}^L \tilde{A}^{(l)}$을 수행하는 수리적 근거를 설명하시오. 만약 항등 행렬 $I$를 더하지 않고 순수 어텐션 행렬만을 연속 곱($\prod_{l=1}^L A^{(l)}$)할 경우, 마르코프 연쇄(Markov Chain)의 에르고딕성(Ergodicity)과 정류 분포(Stationary Distribution) 관점에서 어텐션 맵에 어떠한 수학적 붕괴가 일어나는지 증명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 잔차 연결 항등 행렬 $I$의 수리적 필요성**:
  - 트랜스포머의 각 계층은 입력 벡터 $\mathbf{x}^{(l-1)}$에 셀프 어텐션 출력 $\text{MHA}(\mathbf{x}^{(l-1)})$을 더하는 잔차 구조 $\mathbf{x}^{(l)} = \mathbf{x}^{(l-1)} + \text{MHA}(\mathbf{x}^{(l-1)})$를 채택합니다.
  - 즉, 정보의 일부는 어텐션을 거치지 않고 그대로 다음 계층으로 통과(Identity Mapping)합니다.
  - 따라서 이전 계층의 유전자 표상이 자기 자신에게 보존되는 비율을 반영하기 위해 단위 행렬 $I$를 어텐션 행렬 $A^{(l)}$에 동일 가중치(0.5)로 더해주어야만 실제 신경망 내부의 잔차 정보 흐름과 수리적으로 일치하게 됩니다.

- **2. 순수 어텐션 곱($\prod_{l=1}^L A^{(l)}$) 시의 수학적 붕괴 증명 (마르코프 수렴)**:
  - 소프트맥스를 통과한 어텐션 행렬 $A$의 각 행은 합이 1이고 원소가 모두 양수($A_{ij} > 0$)인 확률 행렬(Stochastic Matrix / Transition Matrix)입니다.
  - 계층 간 어텐션 행렬의 단순 곱 $P^{(L)} = A^{(L)} A^{(L-1)} \dots A^{(1)}$은 비균일 시간 마르코프 연쇄(Non-homogeneous Markov Chain)의 다단계 전이 확률 행렬과 수학적으로 동일합니다.
  - 페론-프로베니우스 정리(Perron-Frobenius Theorem) 및 마르코프 연쇄의 수렴 정리에 의해, 원소가 모두 양수인 기약적이고 비주기적인(Irreducible & Aperiodic) 확률 행렬들을 무한히 곱하면, 행렬은 모든 행이 동일한 유일한 정류 분포(Stationary Distribution) 벡터 $\boldsymbol{\pi}^T$로 구성되는 랭크 1(Rank-1) 행렬로 수렴합니다:
    $$\lim_{L \to \infty} \prod_{l=1}^L A^{(l)} = \mathbf{1} \boldsymbol{\pi}^T = \begin{pmatrix} \pi_1 & \pi_2 & \dots & \pi_K \\ \pi_1 & \pi_2 & \dots & \pi_K \\ \vdots & \vdots & \ddots & \vdots \\ \pi_1 & \pi_2 & \dots & \pi_K \end{pmatrix}$$
  - **수학적 붕괴 결론**:
    - 계층 $L$이 깊어질수록 모든 유전자 $i$의 어텐션 분포(행 벡터)가 완전히 동일한 $\boldsymbol{\pi}$로 수렴하여, 개별 유전자 간의 차별적 상호작용 정보가 완전히 소실되는 '어텐션 스무딩 붕괴(Over-smoothing Collapse)'가 일어납니다.
    - 항등 행렬 $I$를 더해주는 것은 마르코프 전이에서 자기 자신으로의 전이 확률(Self-loop)을 강하게 유지시켜, 계층이 깊어지더라도 유전자 고유의 정체성이 소실되지 않고 보존되도록 보장하는 수학적 안전장치입니다.

</details>

---

### Q3. [임상 및 생물학적 케이스 스터디 (Clinical Problem Solving)]
Enzalutamide 치료에 실패한 전이성 거세 저항성 전립선암(mCRPC) 환자의 종양 조직 단일 세포 데이터를 Geneformer 모델로 분석하였다.
1) 차등 발현(DEG) 분석 결과에서는 *AR*의 mRNA 발현량이 치료 전과 비교하여 통계적으로 유의미한 차이가 없었으나($p > 0.05$), Geneformer 어텐션 네트워크 분석에서는 *AR*의 아웃디그리 중심성($C_{\text{out}}(\text{AR})$)이 $0.85$에서 $0.12$로 급격히 붕괴되었고, 대신 *ASCL1*과 *SOX2*의 중심성이 상위 0.1%로 폭증하였다. 이 모순적 현상이 시사하는 분자생물학적 병태생리를 진단하시오.
2) 환자의 어텐션 네트워크에서 *FOXA1*과 *HOXB13* 유전자 토큰 간 어텐션이 *AR*에서 약해지고 *ASCL1* 쪽으로 이동한 것으로 관찰되었다(가상 사례; Geneformer 입력에는 프로모터·인핸서 정보가 없으므로 "프로모터 조절부로의 재배선"은 해석할 수 없고 유전자 수준 어텐션 변화로만 읽는다). 이 환자에게 단순 AR 길항제 변경 투여의 효과가 제한적일 수 있는 이유와, 탐색적으로 고려할 수 있는 병용 전략을 논하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 모순적 현상의 분자생물학적 진단 (비기능적 AR과 계통 가소성)**:
  - *AR* mRNA의 양적 발현이 유지됨에도 아웃디그리 중심성($C_{\text{out}}$)이 붕괴되었다는 것은, AR 단백질이 여전히 전사되고 있으나 전립선 상피세포의 전사 조절 프로그램(*KLK3*, *TMPRSS2* 등 표적 유전자 활성화)을 지휘하는 마스터 조절자로서의 기능적 연결성(Functional Wiring)을 완전히 상실했음을 의미합니다.
  - 동시에 *ASCL1*과 *SOX2*의 중심성이 폭증한 것은, 종양 세포가 안드로겐 수용체 신호 전달망으로부터 완전히 이탈하여 신경내분비 및 다능성(Pluripotency) 전사인자가 유전체 전반을 장악하는 **치료 유발 신경내분비암(t-NEPC) 계통 가소성 전이(Lineage Plasticity Switch)**가 이미 진행되었음을 시사합니다(가상 사례 해석).
- **2. 치료 전략 및 최적 병용 요법**:
  - **Darolutamide 단순 변경의 무용성**:
    - AR 네트워크 영향력이 낮아진 것으로 나타났다면(가상 사례), 다른 AR 길항제로 바꿔도 *ASCL1*/*SOX2* 구동 회로에 의존하는 세포에는 효과가 제한적일 수 있습니다(반응이 "전혀" 없다고 단정할 수는 없음).
  - **네트워크 재배선 기반 최적 병용 요법**:
    - 아래는 전임상·초기 임상 수준의 **탐색적 가설**이며 어느 것도 이 분석만으로 선택해서는 안 되고, 치료 결정은 임상 가이드라인과 전문가 판단에 따릅니다(Tazemetostat은 전립선암에 승인되지 않음).
    - 1) **EZH2 억제제**: 후성유전학적 리프로그래밍 억제를 노리는 전임상 가설(AR 길항제와의 병용 포함).
    - 2) **백금계 세포독성 항암요법(Carboplatin + Etoposide)**: 소세포/신경내분비 표현형에 쓰이는 기존 요법으로 임상 가이드라인에서 고려되는 선택지.
    - 3) **DLL3 표적 BiTE(Tarlatamab)**: DLL3-CD3 이중특이성 T세포 인게이저로, ASCL1 억제제가 아니라 **DLL3를 발현하는 세포를 직접 표적**합니다. 승인 적응증은 소세포폐암(2024-05 가속승인)이며, 신경내분비 전립선암에서는 초기 임상시험(1b상)에서 DLL3 양성 선별 하에 예비 활성이 보고된 단계입니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618(7965):616-624. PMID: 37258680. doi: 10.1038/s41586-023-06139-9
2. Abnar S, Zuidema W. Quantifying Attention Flow in Transformers. *Proceedings of the 58th Annual Meeting of the Association for Computational Linguistics*. 2020:4190-4197. doi: 10.18653/v1/2020.acl-main.385
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NetworkX 3.7 (github.com/networkx/networkx).
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- PyTorch 2.14.1 (github.com/pytorch/pytorch).
<!-- REFERENCES:END -->

# Chapter 07. 세포 니치 및 공간 미세환경 클러스터링 (Cellular Niches and Microenvironments)

> **Cognitive Goal**: 개별 단일 세포의 유전체 발현을 넘어, 세포들이 물리적 조직 공간에서 형성하는 기능적 상호작용 단위인 '세포 니치(Cellular Niche / Neighborhood)'의 개념과 수학적 그래프 토폴로지를 정립한다. Delaunay 삼각분할, 반경 그래프, 공간적 $k\text{NN}$ 그래프를 구성하고, 레이블 무작위 순열 기반의 이웃 풍부도 검정(Neighborhood Enrichment Test)을 유도한다. 이어 국소 윈도우 세포 구성비 집계를 통한 비지도 니치 클러스터링 알고리즘을 체계화하고, 전립선암(mCRPC)의 4대 핵심 미세환경인 선암 코어, 섬유성 CAF 피막, 3차 림프구조(TLS), 골전이 용골 전선을 규명한다.

---

## 1. 세포 니치(Cellular Niche)의 개념: 세포 유형을 넘어선 공간 기능 단위

단일 세포 전사체 분석(scRNA-seq)은 조직을 효소 처리하여 개별 세포로 해리(Dissociation)시키므로, "어떤 세포 유형이 존재하는가(Cell Type Catalog)"는 알려주지만 **"그 세포가 누구의 곁에서 어떤 물리적 자극을 받으며 기능하는가(Sociology of Cells)"**에 대한 공간적 맥락(Spatial Context)을 완벽히 상실합니다.

생체 내 종양 조직에서 세포의 운명과 약물 반응성은 세포 자체의 유전적 변이뿐만 아니라, **인접한 1차/2차 이웃 세포들과 형성하는 국소 미세환경 니치(Cellular Microenvironment Niche)**에 의해 결정적으로 지배됩니다.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│                   단일 세포 프로파일 vs 공간적 세포 니치(Niche)                   │
│                                                                                  │
│   [단일 세포 관점]                                                               │
│   - CD8+ 세포독성 T 세포: 단독으로는 세포독성 마커(GZMB, PRF1)를 발현함         │
│                                                                                  │
│   [공간 니치 관점 (Spatial Niche)]                                               │
│   (Case A) TLS 니치 내부:                                                        │
│     - CD8+ T 세포 + CD20+ B 세포 + CXCL13+ CD4+ Tfh 세포 이웃                    │
│     ──► 활성화 항원 제시 및 증식 니치 (면역요법 반응성 우수!)                    │
│                                                                                  │
│   (Case B) 데스모플라시아 CAF 피막 니치:                                         │
│     - CD8+ T 세포 + FAP+/ACTA2+ myCAF + Dense Collagen I 섬유 이웃               │
│     ──► 물리적 포획 및 TGF-β 매개 기능 저하 (면역 배제 / 치료 저항성!)           │
└──────────────────────────────────────────────────────────────────────────────────┘
```

따라서 공간전사체학의 핵심 과제는 개별 세포의 정체성을 규명하는 것을 넘어, **인접 세포들의 조성과 상호작용 패턴으로 정의되는 공간적 니치(Cellular Neighborhoods, CNs)**를 수학적으로 모델링하고 분할하는 것입니다.

---

## 2. 공간 이웃 그래프(Spatial Neighborhood Graphs)의 기하학

세포 니치를 정의하기 위한 첫 번째 수학적 단계는 2차원 또는 3차원 조직 공간에서 세포들의 물리적 연결성을 표현하는 **공간 그래프 $G = (V, E)$**를 구축하는 것입니다:
- 정점 집합 $V = \{v_1, v_2, \dots, v_N\}$: 세포 또는 스팟의 집합.
- 간선 집합 $E \subseteq V \times V$: 물리적으로 인접하여 상호작용 가능한 세포 쌍 $(v_i, v_j)$의 연결.

```
      (A) Delaunay Triangulation        (B) Radius Graph (r)        (C) Spatial kNN (k=4)
               ●───────●                         ●───────●                    ●───────●
              / \     / \                       / \     /                    / \     /
             /   \   /   \                     /   \   /                    /   \   /
            ●─────●───────●                   ●─────●                      ●─────●
             \   / \     /                     \   /                        \   / \
              \ /   \   /                       \ /                          \ /   \
               ●─────●───●                       ●                            ●─────●
```

### 2.1 델로네 삼각분할 (Delaunay Triangulation)
- **정의**: 평면 상의 점 집합 $V$에 대해, 임의의 삼각형의 외접원(Circumcircle) 내부에 $V$의 다른 어떤 점도 포함하지 않도록 정점들을 연결하는 삼각분할 기법입니다.
- **수학적 특성**: 보로노이 다이어그램(Voronoi Diagram)의 기하학적 쌍대 그래프(Dual Graph)입니다. 두 세포 $i$와 $j$의 보로노이 셀(Voronoi Cell)이 모서리를 공유할 때에만 간선 $e_{ij} \in E$가 생성됩니다.
- **장점**: 조직 전체의 세포 밀도 변화(조밀한 종양 코어 vs 성긴 기질 구역)에 유연하게 적응하여, 밀도에 구애받지 않고 자연스러운 1차 접촉 이웃을 완벽하게 포착합니다.
- **단점**: 조직 절편의 가장자리나 세포가 비어 있는 구멍(Tear / Lumen)을 가로지르는 비정상적으로 긴 가상 간선(Spurious long edges)이 형성될 수 있으므로, 최대 허용 거리 $d_{\max}$ 임계값을 통한 가지치기(Pruning)가 필수적입니다.

### 2.2 고정 반경 그래프 (Radius Graph)
- **정의**: 두 세포 사이의 유클리드 거리가 특정 물리적 반경 $r$ 이내일 때에만 간선을 형성합니다:
  $$E = \left\{ (v_i, v_j) \;\mid\; \| s_i - s_j \|_2 \le r \right\}$$
- **임계 반경의 선택 기준**:
  - **단일 세포 해상도 플랫폼 (Xenium, CosMx, MERSCOPE)**: 단일 세포 직경($10 \sim 15\ \mu\text{m}$)과 방출된 사이토카인의 유효 확산 거리를 고려하여 **$r = 30 \sim 50\ \mu\text{m}$**를 설정합니다 (직접 접촉 및 인접 파라크라인 니치; 코드·그림·파이프라인은 이 범위의 예시값 $r = 30\ \mu\text{m}$, $k = 20$을 사용).
  - **스팟 기반 플랫폼 (10x Visium)**: 스팟 중심 간 거리($100\ \mu\text{m}$)를 감안하여 인접한 6개 스팟을 포섭하는 **$r = 120 \sim 150\ \mu\text{m}$**를 설정합니다.

### 2.3 공간적 $k$-최근접 이웃 그래프 (Spatial $k\text{NN}$)
- **정의**: 각 세포 $i$에 대해 유클리드 거리가 가장 가까운 상위 $k$개의 세포를 연결합니다 ($E = \{(v_i, v_j) \mid v_j \in \mathcal{N}_k(v_i)\}$).
- 대칭성을 확보하기 위해 통상 무방향 그래프(Undirected Graph, 즉 $i$가 $j$의 이웃이거나 $j$가 $i$의 이웃인 경우 간선 유지)로 변환합니다.

---

## 3. 공간 공기(Co-occurrence) 점수 및 이웃 풍부도 검정 (Neighborhood Enrichment Test)

공간 그래프 $G$ 상에서 특정 두 세포 유형 $A$와 $B$가 무작위 기대치보다 더 자주 물리적으로 인접해 있는가(Enrichment / Attraction), 아니면 서로를 배제하여 격리되어 있는가(Depletion / Avoidance)를 통계적으로 엄밀히 검정하는 방법이 **이웃 풍부도 검정(Neighborhood Enrichment Test)**입니다.

### 3.1 가설 검정 통계량의 유도
그래프 $G = (V, E)$에서 세포 유형 $A$에 속하는 세포 정점 집합을 $V_A$, 세포 유형 $B$에 속하는 집합을 $V_B$라 할 때, 두 세포 유형 간의 연결 간선 수 관측값 $e_{A, B}$는 다음과 같습니다:

$$e_{A, B} = \sum_{(v_i, v_j) \in E} \left[ \mathbb{I}(v_i \in V_A, v_j \in V_B) + \mathbb{I}(v_i \in V_B, v_j \in V_A) \right]$$

### 3.2 레이블 순열 귀무가설 (Label Permutation Null Model)
1. **위상 구조 보존(Topology Preservation)**: 공간 그래프 $G$의 물리적 기하 구조(정점 위치 및 간선 연결 $E$)는 그대로 고정합니다.
2. **세포 유형 무작위 셔플링**: 세포 유형 레이블 벡터 $L = [l_1, l_2, \dots, l_N]$을 정점들에 무작위로 $B$회(통상 $B = 1000$) 순열(Permutation) 배치합니다.
3. 각 순열 반복 $b \in \{1, \dots, B\}$마다 생성된 세포 간선 수 $e_{A, B}^{(b)}$를 계산하여 귀무가설 하의 표본 평균 $\mu_{A, B}^{\text{perm}}$과 표준편차 $\sigma_{A, B}^{\text{perm}}$을 산출합니다:
   $$\mu_{A, B}^{\text{perm}} = \frac{1}{B} \sum_{b=1}^B e_{A, B}^{(b)}, \quad \sigma_{A, B}^{\text{perm}} = \sqrt{\frac{1}{B-1} \sum_{b=1}^B \left( e_{A, B}^{(b)} - \mu_{A, B}^{\text{perm}} \right)^2}$$
4. **표준화 Z-점수(Z-score) 도출**:
   $$Z_{A, B} = \frac{e_{A, B} - \mu_{A, B}^{\text{perm}}}{\sigma_{A, B}^{\text{perm}}}$$

- **$Z_{A, B} > +2$ ($p < 0.05$)**: 세포 유형 $A$와 $B$가 통계적으로 유의미하게 물리적으로 공존 및 응집(Spatial Colocalization).
- **$Z_{A, B} < -2$**: 세포 유형 $A$와 $B$가 물리적 장벽에 의해 공간적으로 격리 및 배제(Spatial Exclusion).

> **해석상 한계 (전역 순열)**: 조직 전체에서 레이블을 전역 셔플하는 귀무는 차수(degree) 이질성은 통제하지만 조직의 **구획(compartment) 구조**는 통제하지 못합니다. 같은 구획에 공존하는 두 유형은 서로 끌어당기는 상호작용이 없어도 큰 $Z$를 얻을 수 있습니다(검증 시뮬레이션: 40×40 격자, 좌반부 A/B 무작위 혼합·우반부 C일 때 전역 셔플 $Z=+23.9$, 구획 내 셔플 $Z=-0.39$). 따라서 $Z>2$는 '국소 친화성'이 아니라 '구획 공존'일 수 있으며, 구획(ROI/니치) 내부 셔플 또는 계층 순열을 병행하고 해석을 '공간 공존'으로 제한해야 합니다.

---

## 4. 세포 니치 클러스터링(Cellular Neighborhood Clustering) 알고리즘

Schürch et al. (Cell, 2020)이 CODEX 데이터에 적용한 **세포 니치 클러스터링(Cellular Neighborhood Windowing & Clustering)**(각 세포의 최근접 10개 이웃을 윈도우로 삼아 구성비를 k-means로 군집화, 9~10개 CN)은 (공간 세포 관계의 쌍별 분석은 Goltsev et al., Cell 2018의 Delaunay 접촉 그래프·log-odds 분석에 기반) 세포의 개별 발현량이 아니라, **세포 주변의 국소적 세포 구성비(Local Cell Type Proportions)**를 특징 벡터로 변환하여 비지도 클러스터링을 수행하는 알고리즘입니다.

```
                      국소 윈도우 세포 구성비 산출 및 니치 클러스터링
                      
  [1. 세포 좌표 및 유형 레이블]             [2. 국소 윈도우 (반경 R 또는 k=20 이웃; 예시)]
       (A) (T) (T)                             중심 세포 i 주변의 이웃 세포 조성 계수:
     (A) [i] (F)                               - AR+ 암세포 (A) : 60%
       (F) (F)                                 - myCAF 섬유 (F) : 30%
                                               - CD8+ T세포 (T) : 10%
          │                                              │
          ▼                                              ▼
  [3. N × K 국소 조성 행렬 구축]           [4. 비지도 군집화 (Leiden / K-means)]
     세포   AR+  myCAF  CD8_T  TAM ...            │
     1   [ 0.60  0.30   0.10   0.00 ]             ▼
     2   [ 0.05  0.80   0.00   0.15 ] ──► [최종 세포 니치 (Cellular Niches)]
     :                                      - Niche 1: 선암 상피 코어 (Adeno Core)
     N   [ 0.00  0.10   0.70   0.20 ]       - Niche 2: 데스모플라시아 CAF 장벽
                                            - Niche 3: 면역 허브 (TLS Immune Hub)
```

### 4.1 알고리즘 세부 단계

1. **국소 이웃 윈도우 정의**:
   - 각 세포 $i \in \{1, \dots, N\}$에 대해 공간 그래프 상의 $k$-최근접 이웃 집합 $\mathcal{N}_k(i)$(자기 자신 포함)를 추출합니다 (Schürch 2020은 $k = 10$을 사용; 데이터 밀도에 따라 $k = 10 \sim 50$ 범위를 시도).
2. **국소 조성 벡터(Local Composition Vector) 산출**:
   - 총 $K$개의 세포 유형에 대해, 이웃 집합 내에서 세포 유형 $c \in \{1, \dots, K\}$가 차지하는 빈도 비율을 계산하여 $K$차원 벡터 $v_i \in \mathbb{R}^K$를 생성합니다:
     $$v_{i, c} = \frac{1}{|\mathcal{N}_k(i)|} \sum_{j \in \mathcal{N}_k(i)} \mathbb{I}(\text{CellType}(j) = c), \quad \sum_{c=1}^K v_{i, c} = 1$$
3. **전체 조직 조성 행렬 $M \in \mathbb{R}^{N \times K}$ 구성**:
   - 모든 세포에 대한 국소 벡터를 적층하여 행렬을 구축하고, 희소 세포 유형의 과소평가를 방지하기 위해 각 열(세포 유형)에 대해 z-점수 표준화 또는 MinMax 스케일링을 적용합니다.
4. **니치 군집화 및 공간 투영**:
   - 행렬 $M$에 대해 유클리드 거리를 기반으로 $k$-means 군집화 또는 Leiden 그래프 클러스터링을 실행하여 최종 $M$개의 공간적 니치 레이블을 할당합니다.

---

## 5. 전립선암(mCRPC)의 4대 핵심 미세환경 니치

mCRPC 공간전사체 분석에서 나타날 수 있다고 가설화한 4대 대표적 공간 니치(문헌에서 공통적으로 확립된 분류가 아닌 설명용 틀)의 병리학적 아키텍처는 다음과 같습니다:

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│              전립선암 전이 및 치료 저항성을 지배하는 4대 공간 니치               │
│                                                                                  │
│   [1. 선암 코어 니치]               [2. 데스모플라시아 CAF 피막 니치]            │
│   - AR+ 관강세포 밀집 군집          - ACTA2+/FAP+ myCAF 촘촘한 띠 형성          │
│   - 선강(Lumen) 형성 및 AMACR+      - 고밀도 Type I Collagen 물리적 차단벽       │
│   - 혈관 배제에 의한 저산소 중심부  - CD8+ T세포 침투 차단 (면역 배제 표현형)   │
│                                                                                  │
│   [3. 3차 림프구조 (TLS) 니치]      [4. 골전이 경계 니치 (혼합형, 골형성 우세)]│
│   - CD20+ B 세포 여포 (Follicle)    - 파골세포 (TRAP+/RANK+) 골 흡수 전선       │
│   - CD3E+/CD8A+ T 세포 경계 구역    - RANKL(주로 골아세포계 유래) 매개 골 흡수 요소│
│   - CXCL13/CCL19 케모카인 분비 활성 - 기질 내 TGF-β 대량 방출 ──► 종양 증식     │
└──────────────────────────────────────────────────────────────────────────────────┘
```

### 5.1 4대 니치의 분자병리학적 특성 비교

| 니치 명칭 | 주요 구성 세포 (Dominant Cell Types) | 핵심 분자 마커 및 경로 | 임상적 의의 및 치료 표적 |
| :--- | :--- | :--- | :--- |
| **1. Cribriform Adenocarcinoma Core** | AR+ Luminal 종양 상피세포, 분화 선암세포 | *AR, KLK3, AMACR, FOXA1, HOXB13* | 호르몬 반응성 중심 구역. AR 표적 치료(엔잘루타마이드)의 주 공격 표적 |
| **2. Desmoplastic CAF Rim** | 근섬유아세포성 myCAF, 혈관주위세포, ECM | *ACTA2, FAP, COL1A1, TGFB1, POSTN* | 물리적 데스모플라시아 장벽. T 세포의 종양 코어 침투를 기계적으로 차단 |
| **3. TLS Immune Hub** | CD20+ B 세포, CD4+/CD8+ T 세포, 성숙 DC | *MS4A1(CD20), CD3E, CD8A, CXCL13, LAMP3* | 국소 항종양 체액성/세포성 면역 유도. 면역관문억제제(ICI) 반응성의 양성 예측 인자 |
| **4. Bone Metastasis Interface (혼합형, 골형성 우세)** | 종양 세포(AR+ 또는 t-NEPC/DNPC), 골아세포, 파골세포(Osteoclasts) | *TNFSF11(RANKL), TNFRSF11A(RANK), CTSK, PTHLH* | 전립선암 골전이는 주로 골형성성(osteoblastic)·혼합형이며(Yin 2003) 순수 용골성은 소수입니다. 골 흡수 요소가 있는 경우 파골세포 활성화와 TGF-$\beta$ 유리(악순환 고리) 가설이 적용되며 데노수맙은 골 관련 사건 예방에 쓰입니다 |

---

## 6. 실전 파이썬 워크플로우: Squidpy 공간 이웃 및 니치 클러스터링

다음은 단일 세포/스팟 해상도 공간전사체 데이터에서 공간 그래프를 구축하고, 이웃 풍부도를 검정한 뒤 국소 윈도우 기반 니치 클러스터링을 수행하는 실전 파이프라인입니다:

```python
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------
# Step 1: 세포 유형 주석이 완료된 AnnData 로드
# ---------------------------------------------------------
adata = sc.read_h5ad("data/prostate_cancer_annotated_cells.h5ad")

# ---------------------------------------------------------
# Step 2: 공간 이웃 그래프(Delaunay Triangulation / Radius) 구축
# ---------------------------------------------------------
# 반경 30 µm 고정 반경 그래프 생성 (단일 세포 해상도 기준; 예시값, 파이프라인 02와 동일)
sq.gr.spatial_neighbors(
    adata,
    coord_type="generic",
    spatial_key="spatial",
    radius=30.0                   # 마이크로미터 단위 유효 반경 (예시값)
)

# ---------------------------------------------------------
# Step 3: 이웃 풍부도 검정 (Neighborhood Enrichment Test)
# ---------------------------------------------------------
# 세포 레이블 1,000회 순열을 통한 Z-score 행렬 계산
sq.gr.nhood_enrichment(
    adata,
    cluster_key="cell_type",
    n_perms=1000,
    seed=42
)

# Z-score 히트맵 시각화 (CAF와 T 세포의 배제 여부 확인)
sq.pl.nhood_enrichment(
    adata,
    cluster_key="cell_type",
    figsize=(6, 5),
    cmap="vlag",
    vmin=-50,
    vmax=50
)

# ---------------------------------------------------------
# Step 4: 국소 윈도우 집계를 통한 세포 니치(Cellular Neighborhood) 산출
# ---------------------------------------------------------
def compute_cellular_neighborhoods(adata, cluster_key, k_neighbors=20, n_niches=4):
    """각 세포의 국소 이웃 세포 구성비를 계산하고 니치 클러스터링을 수행합니다."""
    # 1. 공간 좌표 기반 kNN 그래프 탐색
    from sklearn.neighbors import NearestNeighbors
    coords = adata.obsm["spatial"]
    nn = NearestNeighbors(n_neighbors=k_neighbors).fit(coords)
    _, indices = nn.kneighbors(coords)
    
    # 2. 세포 유형 원핫 인코딩
    cell_types = pd.get_dummies(adata.obs[cluster_key])
    K = cell_types.shape[1]
    
    # 3. 이웃 윈도우 조성 행렬 M 계산 (N x K)
    N = adata.n_obs
    composition_matrix = np.zeros((N, K))
    for i in range(N):
        neighbor_types = cell_types.iloc[indices[i]]
        composition_matrix[i, :] = neighbor_types.mean(axis=0).values
        
    # 4. 스케일링 및 비지도 K-means 군집화
    scaled_comp = StandardScaler().fit_transform(composition_matrix)
    kmeans = MiniBatchKMeans(n_clusters=n_niches, random_state=42)
    niche_labels = kmeans.fit_predict(scaled_comp)
    
    # 5. AnnData 메타데이터에 등록
    adata.obs["cellular_niche"] = pd.Categorical([f"Niche_{l}" for l in niche_labels])
    adata.obsm["niche_composition"] = composition_matrix
    return adata

adata = compute_cellular_neighborhoods(adata, cluster_key="cell_type", k_neighbors=30, n_niches=4)

# ---------------------------------------------------------
# Step 5: 복원된 4대 세포 니치 공간 맵 시각화
# ---------------------------------------------------------
with plt.rc_context({"figure.figsize": (7, 7)}):
    sc.pl.spatial(
        adata,
        color="cellular_niche",
        palette=["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"],
        size=1.2,
        title="Prostate Cancer Cellular Niches"
    )
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
이웃 풍부도 검정(Neighborhood Enrichment Test)에서 그래프 토폴로지(정점 위치 및 간선 집합 $E$)를 고정하고 오직 세포 유형 레이블만을 무작위 셔플링(Permutation)하는 수학적 이유를 '조직 밀도 편향(Tissue Density Confounding)'의 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **조직 기하 구조 및 밀도 불균일성의 통제**:
  - 생체 조직의 세포 밀도는 균일하지 않습니다. 종양 코어 구역은 세포가 극도로 조밀하여 정점당 간선 차수(Degree, 연결된 이웃 수)가 15~20개에 달하는 반면, 섬유화 기질이나 괴사 구역은 정점당 간선 수가 2~3개에 불과합니다.
  - 만약 그래프 구조 자체를 무작위화(예: Erdős–Rényi 랜덤 그래프 생성)하면, 원래 조직이 가진 고유한 공간적 밀도 불균형과 표면적 대 부피비(Surface-to-volume ratio) 정보가 완전히 파괴됩니다.
- **순수한 세포 배치 선호도(Attraction/Repulsion) 분리**:
  - 그래프 위상($E$)을 완벽히 보존한 상태에서 레이블만 셔플링하면, "동일한 밀도와 동일한 이웃 연결 수를 가진 정점 환경에서, 두 세포 유형이 우연에 의해 인접할 확률 기댓값($\mu^{\text{perm}}$)"을 정확히 측정할 수 있습니다.
  - 이를 통해 고밀도 구역 밀집(차수 불균일)에 의한 효과는 통제할 수 있습니다. 다만 조직의 구획 구조는 통제되지 않으므로, 전역 순열의 높은 $Z$를 곧바로 케모카인 유인이나 물리적 부착에 의한 '진정한 미세환경 친화성'으로 해석할 수 없고 '공간 공존'으로 제한해 읽어야 하며, 구획 내 셔플로 확인해야 합니다.

</details>

---

### Q2. [알고리즘 비교 분석 (Algorithmic Rigor)]
단일 세포의 유전자 발현량 벡터($x_i \in \mathbb{R}^G$)를 직접 클러스터링하는 것과, 본 단원에서 다룬 국소 이웃 세포 조성비 벡터($v_i \in \mathbb{R}^K$)를 클러스터링하는 '세포 니치 클러스터링'이 제공하는 정보의 차이를 비교하고, 동일한 세포 유형(예: CD8+ T 세포)이 서로 다른 니치에 속할 때 나타나는 생물학적 상태 전이를 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **표현 공간의 근본적 차이**:
  - **발현 벡터 클러스터링 ($x_i \in \mathbb{R}^G$)**: 세포의 고유한 세포 계통 및 분화 상태(Cell Identity)를 규명합니다 (예: "이 세포는 CD8+ 세포독성 T 세포이다").
  - **니치 조성 클러스터링 ($v_i \in \mathbb{R}^K$)**: 세포가 처한 물리적 사회학적 환경(Spatial Microenvironment / Context)을 규명합니다 (예: "이 세포는 myCAF 피막 니치에 포획되어 있는가, 아니면 TLS 면역 니치에 존재하는가").
- **동일 세포 유형의 니치 의존적 상태 전이 (CD8+ T 세포 사례)**:
  - **TLS 니치 소속 CD8+ T 세포**:
    - CD20+ B 세포 및 성숙 수지상세포의 지속적인 항원 제시와 IL-2/IL-15 사이토카인 자극을 받아 높은 줄기세포성/기억 상태(*TCF7+, IL7R+*)를 유지하며 활발히 증식합니다.
  - **데스모플라시아 myCAF 피막 니치 소속 CD8+ T 세포**:
    - myCAF가 분비하는 고농도 활성 TGF-$\beta$와 콜라겐 섬유의 물리적 기계적 압박(Extracellular Matrix Stiffness)을 받아 *PDCD1, HAVCR2, LAG3*를 발현하는 말초 소진(Terminal Exhaustion) 상태로 전이되거나 암세포와의 직접 접촉이 원천 차단됩니다.
  - 즉, 니치 클러스터링을 통해서만 세포의 기능적 활성도가 국소 이웃 상호작용에 의해 어떻게 재프로그래밍되는지를 규명할 수 있습니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
mCRPC 환자의 척추 뼈 전이(Spinal Bone Metastasis) 조직 절편에 대해 고해상도 공간전사체 분석을 수행하였다. 연구팀은 비지도 니치 클러스터링을 통해 골 기질 경계부에서 독특한 'Niche 4'를 발견하였다. 
Niche 4의 세포 조성은 t-NEPC 신경내분비 암세포(40%), TRAP/CTSK 양성 다핵 파골세포(30%), 활성화 골아세포(20%), 골수성 면역억제세포(10%)로 구성되어 있었다. 이 니치에서 유의하게 과발현되는 리간드-수용체 쌍은 *TNFSF11(RANKL)-TNFRSF11A(RANK)* 및 *TGFB1-TGFBR2*였다. 
이 뼈 전이 니치(가상 시나리오; 전립선암 골전이는 주로 골형성성·혼합형이므로 여기서는 골 흡수 요소가 두드러진 혼합형 병변을 가정)에서 일어날 수 있는 **'골 흡수 악순환(Osteolytic Vicious Cycle)'**의 병태생리학적 가설을 수립하고, 이 환자에게 적용 가능한 병용 표적 치료 전략을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **골 흡수 악순환(Osteolytic Vicious Cycle)의 분자병태 기전 (가설)**:
  1. **파골세포 과활성화 촉발**: RANKL(*TNFSF11*)은 주로 골아세포계·기질 유래이며 종양 세포의 RANKL 발현은 일부 보고가 있으나 t-NEPC 특이 근거는 확인되지 않았습니다. 이 가설적 시나리오에서는 *PTHLH*(전립선암 세포에서 발현 보고)와 RANKL 신호가 미세환경 내 파골세포(Osteoclasts) 표면의 *TNFRSF11A*(RANK) 수용체를 지속적으로 자극합니다.
  2. **골 기질 파괴 및 성장인자 유리**: 과활성화된 파골세포가 카텝신 K(*CTSK*)와 염산을 분비하여 석회화된 뼈 기질을 급속히 용해(Osteolysis)시킵니다. 이 과정에서 정상 뼈 기질 내부에 불활성 상태로 고착되어 있던 막대한 양의 형질전환성장인자-베타(**TGF-$\beta$**) 및 IGF-1이 가용성 형태로 대량 방출됩니다.
  3. **종양 증식 및 면역 억제 증폭**: 방출된 TGF-$\beta$는 종양 세포의 TGFBR2에 결합하여 암세포의 침습을 촉진하고 면역억제 환경에 기여할 수 있습니다(정량적 기전은 이 사례로 입증되지 않은 가설). 이로써 자가 증폭적인 악순환 루프가 *가능*합니다.
- **병용 표적 치료 전략**:
  - **1. 골 파괴 루프 차단**: 인간화 단일클론 항체 **데노수맙(Denosumab)**은 RANKL을 표적으로 하며, CRPC에서는 골관련 사건(SRE) 발생을 지연시키는 효과로 승인되었습니다(Fizazi 2011 Lancet; 종양 억제나 생존 이득이 입증된 것은 아님). 따라서 "악순환 차단"은 기전 가설입니다.
  - **2. 신경내분비 종양 사멸 및 기질 리모델링**: t-NEPC 세포를 사멸시키기 위한 백금 기반 항암화학요법(Cisplatin/Etoposide)과 함께, TGF-$\beta$ 신호를 억제하는 **TGF-$\beta$ 수용체 키나아제 억제제**(예: Galunisertib; 임상 개발 현황은 확인이 필요하며 현재 표준 요법이 아님)의 병용을 탐색적 가설로 제시할 수 있습니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Schürch CM, Bhate SS, Barlow GL, et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. *Cell*. 2020;182(5):1341-1359.e19. PMID: 32763154. doi: 10.1016/j.cell.2020.07.005
2. Goltsev Y, Samusik N, Kennedy-Darling J, et al. Deep Profiling of Mouse Splenic Architecture with CODEX Multiplexed Imaging. *Cell*. 2018;174(4):968-981.e15. PMID: 30078711. doi: 10.1016/j.cell.2018.07.010
3. Yin JJ, Mohammad KS, Käkönen SM, et al. A causal role for endothelin-1 in the pathogenesis of osteoblastic bone metastases. *Proc Natl Acad Sci U S A*. 2003;100(19):10954-9. PMID: 12941866. doi: 10.1073/pnas.1830978100
4. Fizazi K, Carducci M, Smith M, et al. Denosumab versus zoledronic acid for treatment of bone metastases in men with castration-resistant prostate cancer: a randomised, double-blind study. *Lancet*. 2011;377(9768):813-22. PMID: 21353695. doi: 10.1016/S0140-6736(10)62344-6
### 검수 근거 문헌
1. Smith MR, Coleman RE, Klotz L, et al. Denosumab for the prevention of skeletal complications in metastatic castration-resistant prostate cancer: comparison of skeletal-related events and symptomatic skeletal events. *Ann Oncol*. 2015;26(2):368-74. PMID: 25425475. doi: 10.1093/annonc/mdu519
2. Hashimoto K, Ochi H, Sunamura S, et al. Cancer-secreted hsa-miR-940 induces an osteoblastic phenotype in the bone metastatic microenvironment via targeting ARHGAP1 and FAM134A. *Proc Natl Acad Sci U S A*. 2018;115(9):2204-2209. PMID: 29440427. doi: 10.1073/pnas.1717363115
3. Kfoury Y, Baryawno N, Severe N, et al. Human prostate cancer bone metastases have an actionable immunosuppressive microenvironment. *Cancer Cell*. 2021;39(11):1464-1478.e8. PMID: 34719426. doi: 10.1016/j.ccell.2021.09.005
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- Squidpy 1.8.3 (github.com/scverse/squidpy). 원 논문: Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nat Methods*. 2022;19(2):171-178. PMID: 35102346. doi: 10.1038/s41592-021-01358-2
- scikit-learn 1.9.1 (github.com/scikit-learn/scikit-learn).
<!-- REFERENCES:END -->

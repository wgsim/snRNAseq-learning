# Chapter 06. 공간 전사체 기반 세포 간 통신 통합 분석 (Spatial CCC & LIANA+)

> **Cognitive Goal**: 세포의 물리적 위치 정보가 소실된 단일 세포 전사체(scRNA-seq) 분석에서 발생하는 수 밀리미터($\text{mm}$) 규모의 비물리적 허위 양성(False Positives)을 극복하기 위해, Track 02의 공간 전사체(Spatial Transcriptomics) 물리 좌표계를 리간드-수용체 상호작용 추론에 결합하는 원리를 이해한다. 공간 이변량 모란 지수(Spatial Bivariate Moran's $I$) 및 국소 공간 교차 상관(Local Cross-Correlation)의 수학적 원리를 도출하고, 신호 전달 양식별 생물리학적 거리 감쇄 함수(Distance Decay Function)를 적용하는 방법론을 체득한다. 최종적으로 `LIANA+` 프레임워크를 활용하여 전립선암 골전이 미세환경의 공간 접촉형/분비형 신호 허브를 정밀 규명한다.

---

## 1. 공간 좌표 결합의 필연성: "1 mm 거리의 L-R 통신은 물리적 허구이다"

단일 세포 전사체(scRNA-seq) 분석에서는 조직을 효소로 해리하여 세포 현탁액(Single-cell Suspension)을 만든 후 시퀀싱을 진행합니다. 이로 인해 모든 세포가 동일한 시험관 내에서 자유롭게 만날 수 있다는 비현실적인 가정이 수반됩니다. 

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               해리 기반 scRNA-seq vs. 공간 전사체 기반 세포 간 통신의 차이             │
│                                                                                        │
│   [해리 기반 scRNA-seq CCC (CellPhoneDB, NATMI 등)]                                    │
│   - 물리적 좌표 상실 (x, y 좌표 없음)                                                  │
│   - 세포군 A (리간드 고발현) ──?──► 세포군 B (수용체 고발현)                           │
│   - 판정: 두 유전자의 통계적 발현량이 높으므로 유의미한 상호작용으로 판단!             │
│   - 실제 조직 내 실체: 세포 A는 생검 조직 좌측 끝, 세포 B는 3 mm 우측 끝에 위치!      │
│     ==> 확산 불가능한 원거리! 100% 생체 내 불가능한 거짓 양성(False Discovery)!        │
│                                                                                        │
│   [공간 전사체 기반 CCC (LIANA+)]                                                      │
│   - 마이크로미터(µm) 단위의 절대 공간 좌표 (s_i, s_j) 보유                             │
│   - 세포 간 유클리드 거리 d_{ij} = ||s_i - s_j||_2를 명시적 가중치로 부여              │
│   - 신호 양식별 생물리학적 확산 한계(d_0) 이내에 존재하는 유효 세포 쌍만 평가         │
│     ==> 위양성을 원천 배제하고 실제 조직 마이크로니치(Microniche) 내 통신만 포착!     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 조직 내 분자 확산의 생물리학적 한계
생체 조직의 세포외 기질(ECM)은 콜라겐, 라미닌, 피브로넥틴 및 프로테오글리칸이 얽힌 조밀한 하이드로겔 그물망입니다. 조직 내에서 분비 단백질의 확산은 픽의 확산 법칙(Fick's Law)과 세포외 분해 효소(MMPs, 기타 프로테아제) 및 세포 표면 엔도사이토시스에 의한 소모에 의해 급격히 제한됩니다.

- **유효 확산 한계**: 대부분의 가용성 사이토카인과 케모카인은 분비 지점으로부터 대략 $100 \sim 250\ \mu\text{m}$ 수준을 넘어서면 유효 농도($K_d$ 이상)를 유지하기 어렵다는 것이 통용되는 경험칙입니다(분자·조직에 따라 크게 다르며 이 수치의 1차 출처는 제시하지 못했습니다).
- **막 결합형 접촉 신호**: Notch-DLL, Ephrin-Eph, Cadherin, 면역관문(PD-L1/PD-1, CD47-SIRPA) 등은 두 세포막 사이의 물리적 간극이 $10 \sim 50\ \text{nm}$, 유효 중심 거리가 대략 $15\ \mu\text{m}$ 이내인 인접 세포 간에서만 작동하는 것으로 간주합니다(경험칙).

따라서 공간 좌표를 통합하지 않은 세포 간 통신 분석 결과는 고형암 미세환경의 구조적 구획화(Compartmentalization)를 반영하지 못합니다.

---

## 2. 생물리학적 거리 감쇄 함수 (Distance Decay Function)

두 세포 $i$와 $j$ 사이의 2차원 공간 유클리드 거리를 $d_{ij} = \| s_i - s_j \|_2$라 할 때, 리간드 분자가 물리적으로 도달하여 수용체에 결합할 확률적 가중치 $w_{ij}$는 거리 감쇄 커널(Kernel) 함수로 정의됩니다.

### 2.1 지수 감쇄 커널 (Exponential Decay Kernel)
정상 상태의 1차 반응–확산 방정식(Reaction-Diffusion Equation)에서 1차원의 경우 해가 $\exp(-x/\lambda)$ 형태의 지수 감쇄가 됩니다(2차원은 변형 베셀 함수 $K_0$, 3차원은 $\exp(-r/\lambda)/r$ 형태이므로 지수 감쇄는 근사적 선택입니다). 이 형태가 가용성 분비 리간드에 널리 적용됩니다:

$$w_{ij} = \begin{cases} \exp\left( - \frac{d_{ij}}{d_0} \right) & \text{if } d_{ij} \le d_{\max} \\ 0 & \text{if } d_{ij} > d_{\max} \end{cases}$$

- $d_0$: **특성 신호 범위(Characteristic Signaling Range)**. 신호 강도가 $1/e \approx 36.8\%$로 감소하는 물리적 거리.
- $d_{\max}$: 계산 복잡도를 제어하고 수치적 잡음을 차단하기 위한 최대 유효 컷오프 거리 (통상 $3 d_0$).

### 2.2 가우스 감쇄 커널 (Gaussian Decay Kernel)
리간드가 등방성(Isotropic)으로 확산되며 조직 내 농도 구배가 완만한 정규분포를 따를 때 사용됩니다:

$$w_{ij} = \exp\left( - \frac{d_{ij}^2}{2 d_0^2} \right)$$

### 2.3 생물학적 신호 양식별 특성 거리 $d_0$ 가이드라인
```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                 신호 전달 양식별 특성 거리 (d_0) 경험적 가이드(출처 없음)               │
│                                                                                        │
│   1. 세포 접촉 의존형 (Juxtacrine / Contact-dependent):                               │
│      - 대표 리간드-수용체: Notch-DLL1/4, CD80/CD86-CD28, CD274(PD-L1)-PDCD1(PD-1)      │
│      - 권장 d_0: 10 ~ 15 µm (직접 맞닿은 물리적 이웃 세포)                             │
│                                                                                        │
│   2. 국소 단거리 파라크라인 (Short-range Paracrine):                                   │
│      - 대표 리간드-수용체: WNT family, Hedgehog, Ephrin, 단거리 케모카인               │
│      - 권장 d_0: 50 ~ 100 µm (세포 2~5개 직경 범위)                                    │
│                                                                                        │
│   3. 원거리 파라크라인 (Long-range Paracrine):                                         │
│      - 대표 리간드-수용체: TGF-β1, VEGFA, IL-6, TNFSF11(RANKL), CXCL12                  │
│      - 권장 d_0: 150 ~ 250 µm (미세환경 내 최대 기능적 구배)                           │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. 공간 이변량 모란 지수 (Spatial Bivariate Moran's $I$)의 수학적 유도

공간 전사체 격자 상에서 리간드 $L$을 발현하는 세포들과 수용체 $R$을 발현하는 세포들이 무작위 분포를 벗어나 물리적으로 유의미하게 공간적 공존(Spatial Co-occurrence)을 이루는지를 정량화하기 위해 **공간 이변량 모란 지수(Spatial Bivariate Moran's $I$)**를 사용합니다.

### 3.1 전역 이변량 모란 지수 (Global Bivariate Moran's $I_{L, R}$)
(아래는 Anselin형 이변량 모란 $I$의 표준 정의입니다. LIANA+ `li.mt.bivariate`가 계산하는 지표는 SpatialDM 방식의 **Moran's R**(가중치 정규화도 SpatialDM 방식)이며, 이 외에 Lee's L, cosine, pearson, spearman, jaccard, norm_product, masked_spearman 등을 제공하므로 값이 정확히 일치하지는 않습니다. `li.ut.spatial_neighbors`의 커널은 gaussian·exponential·linear·misty_rbf를 지원합니다.)

총 $N$개의 공간 스팟(또는 단일 세포)에 대해, 세포 $i$의 리간드 발현량을 $L_i$, 세포 $j$의 수용체 발현량을 $R_j$라 하고, 각각의 평균을 $\bar{L}$, $\bar{R}$, 표준편차를 $\sigma_L$, $\sigma_R$이라 합니다.

표준화된 $z$-점수 벡터를 다음과 같이 정의합니다:

$$z_{L, i} = \frac{L_i - \bar{L}}{\sigma_L}, \quad z_{R, j} = \frac{R_j - \bar{R}}{\sigma_R}$$

전역 이변량 모란 지수 $I_{L, R}$은 공간 가중치 행렬 $W = [w_{ij}]$를 매개로 한 두 확률 변수의 표준화 공분산으로 정의됩니다:

$$I_{L, R} = \frac{N}{S_0} \frac{\sum_{i=1}^N \sum_{j=1}^N w_{ij} (L_i - \bar{L})(R_j - \bar{R})}{\sqrt{\sum_{i=1}^N (L_i - \bar{L})^2} \sqrt{\sum_{j=1}^N (R_j - \bar{R})^2}}$$

- $S_0 = \sum_{i=1}^N \sum_{j=1}^N w_{ij}$: 전체 공간 가중치의 총합.
- 벡터-행렬 표기법으로 나타내면 다음과 같이 우아하게 표현됩니다:
  $$I_{L, R} = \frac{N}{S_0} \frac{z_L^T W z_R}{N} = \frac{z_L^T W z_R}{S_0}$$
- **행-표준화 가중치 행렬(Row-standardized Weight Matrix, $\sum_j w_{ij} = 1$)**을 사용할 경우 $S_0 = N$이 되므로:
  $$I_{L, R} = \frac{1}{N} z_L^T W z_R$$

### 3.2 수리적 해석
- $I_{L, R} > 0$: **공간적 양의 교차 상관(Positive Cross-Correlation)**. 리간드가 풍부한 위치 주변에 수용체가 통계적으로 유의하게 밀집해 있음을 의미 (신호 전달 활성 핫스팟).
- $I_{L, R} \approx 0$: 공간적 무작위 분포(Spatial Randomness). 두 분자 간에 공간적 연관성이 전혀 없음.
- $I_{L, R} < 0$: **공간적 배제(Negative Cross-Correlation / Spatial Exclusion)**. 리간드가 높은 영역에는 수용체가 체계적으로 배제되어 있음 (예: 면역 회피 장벽).

### 3.3 국소 이변량 모란 지수 (Local Bivariate Moran's $I_i$) 및 LISA 분석
단일 세포 또는 스팟 $i$ 수준에서 신호 핫스팟의 위치를 핀포인트하기 위해 안셀린(Luc Anselin)의 **국소 공간 연관 지표(LISA)**를 도출합니다:

$$I_i(L, R) = z_{L, i} \sum_{j=1}^N w_{ij} z_{R, j}$$

각 위치 $i$는 4개의 사분면(LISA Quadrants)으로 명확히 분류됩니다:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        국소 이변량 모란 지수 사분면 (LISA Quadrants)                  │
│                                                                                        │
│   1. High-High (HH, I_i > 0):                                                          │
│      - 세포 i의 리간드 높음 (z_{L, i} > 0) AND 이웃 j들의 수용체 높음 (W z_R > 0)       │
│      - 생물학적 의미: 강력한 능동적 세포 간 통신 핫스팟(Active Signaling Hotspot)       │
│                                                                                        │
│   2. Low-Low (LL, I_i > 0):                                                            │
│      - 세포 i의 리간드 낮음 AND 이웃 j들의 수용체 낮음                                  │
│      - 생물학적 의미: 신호 침묵 영역(Signaling Desert / Cold Spot)                      │
│                                                                                        │
│   3. High-Low (HL, I_i < 0):                                                           │
│      - 세포 i의 리간드 높음 BUT 이웃 j들의 수용체 결핍                                  │
│      - 생물학적 의미: 리간드가 분비되나 수신자가 없는 '신호 낭비/고립' 영역              │
│                                                                                        │
│   4. Low-High (LH, I_i < 0):                                                           │
│      - 세포 i의 리간드 결핍 BUT 이웃 j들의 수용체 밀집                                  │
│      - 생물학적 의미: 수용체는 준비되었으나 국소 자극이 부재한 영역                    │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 4. 서브셀룰러(Subcellular) 초고해상도 공간 전사체에서의 접촉 통신

10x Xenium, Nanostring CosMx, BGI Stereo-seq과 같은 진정한 단일 세포/세브셀룰러 해상도 기술에서는 캡처 스팟 혼합(Deconvolution) 문제 없이 세포 경계(Cell Segmentation Polygon)를 직접 다룰 수 있습니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                서브셀룰러 공간 플랫폼에서의 이웃 그래프 구축 전략                     │
│                                                                                        │
│   [전략 A: 세그멘테이션 다각형 경계 접촉 (Polygon Boundary Contact)]                   │
│   - 인접 세포 간의 세포막 윤곽선이 물리적으로 맞닿아 있는 경우(공유 경계선 길이 > 0)    │
│   - 가중치: w_{ij} = Shared_Perimeter(i, j)                                            │
│   - 적용: Notch, E-cadherin, Claudin, PD-L1/PD-1 접촉 검증에 최적                      │
│                                                                                        │
│   [전략 B: 들로네 삼각분할 (Delaunay Triangulation) 및 k-NN 그래프]                    │
│   - 세포 중심(Centroid) 좌표 간의 들로네 다포체(Voronoi Diagram의 쌍대 그래프) 형성    │
│   - 물리적 직접 인접 이웃(First-order Neighbors)을 수학적으로 유일하게 결정            │
│   - 가중치: 거리 감쇄 w_{ij} = exp(-d_{ij} / d_0)                                      │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 5. LIANA+ 프레임워크: 다중 알고리즘 합의 및 공간 결합

Dimitrov et al. (Nature Communications 2022, PMID 35680885; LIANA+는 Nature Cell Biology 2024;26:1613–1622, PMID 39223377)에 의해 개발된 **LIANA+ (Ligand-Receptor Interaction Analysis Framework Plus)**는 전 세계에 흩어진 10여 개 이상의 세포 간 통신 알고리즘(CellPhoneDB, SingleCellSignalR, NATMI, Connectome, SCA 등)을 단일 표준 API로 통합하고, 여기에 공간 통계학(SpatialDM 방식의 Moran's R, Lee's L, cosine 유사도 등의 이변량 지표와 MISTy 등)을 융합한 최첨단 프레임워크입니다.

### 5.1 합의 순위화 (Consensus Ranking)의 원리
각 개별 알고리즘 $k \in \{1, \dots, K\}$가 출력한 특정 L-R 쌍의 순위를 $r_k(L, R)$이라 할 때, LIANA/LIANA+의 `rank_aggregate`는 Robust Rank Aggregation(RRA, Kolde 2012; 옵션 `mean`)으로 magnitude 순위와 specificity 순위를 각각 합의 집계합니다:

$$\text{Rank}_{\text{consensus}}(L, R) = \text{aggregate}\big( r_1(L, R), r_2(L, R), \dots, r_K(L, R) \big)$$

공간 정보는 합의 점수에 사후에 곱하는 것이 아닙니다. `rank_aggregate`에 `spatial_key`를 주면 **각 개별 메서드의 점수에 근접(proximity) 가중치**가 적용되고, 공간 이변량 지표(예: Moran's $I$ 기반)는 별도 메서드(`li.mt.bivariate`)로 계산됩니다. 두 결과를 함께 보아 "통계적 발현 강도"와 "공간 밀집도"가 모두 높은 쌍을 후보로 삼을 수 있습니다.

---

## 6. 실전 분석 파이프라인: Python LIANA+ 공간 CCC 워크플로우

다음은 10x Visium 또는 Xenium 기반 전립선암 골전이 공간 전사체 데이터에서 `liana`를 활용하여 공간 이변량 모란 지수를 계산하고 국소 핫스팟을 시각화하는 완전한 파이프라인입니다.

```python
"""
06_spatial_ccc_liana_plus.py
전립선암 골전이 공간전사체(Spatial Transcriptomics) 데이터에서
LIANA+를 활용한 공간 이변량 모란 지수(Bivariate Moran's R) 계산 및
국소 신호 핫스팟 규명 파이프라인
"""

import scanpy as sc
import liana as li
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. 공간 전사체 AnnData 로드 (log 정규화 발현이 adata.X, 좌표는 adata.obsm["spatial"])
print("[Step 1] Loading Spatial Transcriptomics AnnData (mCRPC Bone Metastasis)...")
adata = sc.read_h5ad("mcrpc_bone_metastasis_spatial.h5ad")

# 2. 공간 이웃 그래프 구축 (거리 기반 가우시안 가중치 행렬 W 생성)
# bandwidth는 obsm["spatial"]의 좌표 단위(µm 가정)이며, 골전이 미세환경의 확산 거리를 고려해 150으로 설정
print("[Step 2] Building spatial neighborhood graph with Gaussian distance decay...")
li.ut.spatial_neighbors(
    adata,
    bandwidth=150,
    cutoff=0.1,         # 가중치가 0.1 미만인 이웃은 제거
    kernel="gaussian",
    set_diag=False
)  # adata.obsp["spatial_connectivities"]에 저장

# 3-4. LIANA+ 이변량 공간 지표 계산 (Consensus 리소스, 전역 Moran's R + 국소 Moran's R)
print("[Step 3-4] Computing spatially bivariate metrics across L-R pairs...")
# li.mt.bivariate는 DataFrame이 아니라 AnnData(lrdata: spot x L-R 쌍)를 반환한다.
#   lrdata.X                : 국소 점수 (local_name)
#   lrdata.var              : 쌍별 전역 점수 "morans"와 순열 p-value "morans_pvals"
#   lrdata.layers["pvals"]  : 국소 p-value (n_perms 지정 시)
#   lrdata.layers["cats"]   : add_categories=True일 때 High-High=1, Low-Low=-1, 기타=0
lrdata = li.mt.bivariate(
    adata,
    resource_name="consensus",
    local_name="morans",
    global_name="morans",
    n_perms=1000,        # 몬테카를로 순열 검정 횟수
    add_categories=True,
    seed=42,
    use_raw=False,
    verbose=True
)

# 5. 공간적으로 유의미한 상위 L-R 상호작용 필터링
spatial_res = lrdata.var
top_spatial_lr = spatial_res[
    (spatial_res["morans_pvals"] < 0.01) & (spatial_res["morans"] > 0.15)
].sort_values(by="morans", ascending=False)

print("\n=== Top 10 Spatially Correlated L-R Pairs in Bone Metastasis ===")
print(top_spatial_lr[["ligand", "receptor", "morans", "morans_pvals"]].head(10))

# 6. 핵심 신호 축 국소 분석 (High-High 핫스팟)
# 타겟 L-R 쌍: TNFSF11 (RANKL) -> TNFRSF11A (RANK) 및 CXCL12 -> CXCR4
target_pairs = [("TNFSF11", "TNFRSF11A"), ("CXCL12", "CXCR4")]

for lig, rec in target_pairs:
    print(f"\n[Step 5] Extracting local bivariate Moran's R for {lig} -> {rec}...")
    j = lrdata.var_names.get_loc(f"{lig}^{rec}")   # interaction 이름은 "리간드^수용체"

    lisa_scores = lrdata.X[:, j].toarray().ravel() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X[:, j]).ravel()
    cats = lrdata.layers["cats"][:, j].toarray().ravel()

    adata.obs[f"{lig}_{rec}_lisa_score"] = lisa_scores
    adata.obs[f"{lig}_{rec}_quadrant"] = pd.Categorical(
        np.where(cats == 1, "High-High", "Other")
    )

# 7. 공간적 핫스팟 및 미세환경 지도 시각화
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 7.1 RANKL - RANK 핫스팟 (골용해성 파골세포 활성화 니치)
sc.pl.spatial(
    adata,
    color="TNFSF11_TNFRSF11A_quadrant",
    groups=["High-High"],
    palette={"High-High": "#E41A1C", "Other": "#D3D3D3"},
    ax=axes[0],
    show=False,
    title="RANKL - RANK Active Hotspot (Osteoclast Niche)"
)

# 7.2 CXCL12 - CXCR4 핫스팟 (CAF 매개 T세포 결속 및 면역 배제 띠)
sc.pl.spatial(
    adata,
    color="CXCL12_CXCR4_quadrant",
    groups=["High-High"],
    palette={"High-High": "#377EB8", "Other": "#D3D3D3"},
    ax=axes[1],
    show=False,
    title="CXCL12 - CXCR4 Hotspot (CAF Immune Exclusion Rim)"
)

plt.tight_layout()
plt.savefig("spatial_ccc_lisa_hotspots.png", dpi=300)
plt.close()

print("[Complete] Spatial CCC LIANA+ workflow finished successfully.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
단일 세포 전사체(scRNA-seq) 데이터를 이용한 전통적인 L-R 상호작용 분석 결과와 비교할 때, 공간 전사체(Spatial Transcriptomics) 데이터를 활용하여 거리 감쇄 함수($w_{ij} = \exp(-d_{ij}/d_0)$)를 적용했을 때 제거할 수 있는 결정적 위양성(False Positive)의 유형을 생물리학적 관점에서 2가지 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 유효 확산 거리를 초과하는 원거리 세포 간의 분비 신호 허위 판정 배제**:
  - 생체 조직 내에서 단백질 리간드는 세포외 기질(ECM)의 물리적 저항과 단백질 분해 효소로 인해 유효 확산 거리가 통상 $100 \sim 250\ \mu\text{m}$ 이내로 제한됩니다.
  - scRNA-seq은 $1 \sim 5\ \text{mm}$ 이상 격리되어 실제로는 리간드 분자가 물리적으로 도달할 수 없는 원거리 세포 쌍에 대해서도 단순 발현량이 높다는 이유만으로 유의한 상호작용으로 오판하지만, 거리 감쇄 함수를 적용하면 $d_{ij} > d_{\max}$인 쌍의 가중치가 0으로 수렴하여 이러한 장거리 위양성이 원천 배제됩니다.
- **2. 세포막 직접 접촉이 필수적인 접촉 의존형(Juxtacrine) 신호의 비물리적 오류 교정**:
  - Notch-Delta, Ephrin-Eph, 면역관문(PD-L1/PD-1, CD47-SIRPA) 등은 두 세포막이 직접 맞닿는 근접 간극($< 15\ \mu\text{m}$)에서만 수용체 결합 및 하류 신호가 촉발됩니다.
  - 세포 해리액 가정을 기반으로 한 기존 분석은 완전히 분리된 구획의 세포들 사이에서도 접촉 신호 점수를 과도하게 높게 부여하는 모순을 발생시키며, 공간 분석에서는 인접 경계를 공유하지 않는 세포 쌍의 상호작용을 정확히 기각합니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
공간 이변량 모란 지수(Spatial Bivariate Moran's $I$)의 국소 형태인 $I_i(L, R) = z_{L, i} \sum_j w_{ij} z_{R, j}$에서, LISA 4분면 중 'High-High'와 'High-Low'의 수학적 조건 및 생물학적 의미의 차이를 상세히 비교 설명하고, 단일 세포 시퀀싱의 무작위 노이즈(Dropout)가 $I_{L, R}$ 검정력에 미치는 영향을 완화하기 위해 공간 가중치 행렬 $W$를 행-표준화(Row-standardization)해야 하는 수학적 이유를 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **LISA High-High vs. High-Low의 수리생물학적 차이**:
  - **High-High ($z_{L, i} > 0$ 이고 $\sum_j w_{ij} z_{R, j} > 0$ & $I_i > 0$)**:
    - 세포 $i$의 리간드 발현량이 전체 평균보다 유의하게 높고($z_{L, i} > 0$), 그 물리적 주변 이웃 세포 $j$들의 가중 수용체 발현량 역시 평균보다 높음을 의미합니다.
    - 생물학적으로 실제 리간드가 분비되어 이웃 세포의 수용체와 성공적으로 결합하는 '활성 신호 전달 핫스팟(Active Signaling Hotspot)'을 나타냅니다.
  - **High-Low ($z_{L, i} > 0$ 이고 $\sum_j w_{ij} z_{R, j} < 0$ & $I_i < 0$)**:
    - 세포 $i$는 리간드를 과발현하고 있으나, 그 물리적 이웃 세포들의 수용체 발현량은 평균 이하(결핍 상태)임을 의미합니다.
    - 생물학적으로 리간드가 분비되더라도 신호를 받아줄 수용 세포가 주변에 존재하지 않는 '신호 낭비/고립 영역(Signaling Orphan / Waste)'을 나타냅니다.
- **공간 가중치 행렬 $W$의 행-표준화($\sum_j w_{ij} = 1$) 필요성**:
  - 조직 내 세포 밀도가 불균일할 경우(예: 종양 코어의 초고밀도 영역 vs. 기질 외곽의 저밀도 영역), 특정 세포 $i$의 이웃 세포 수 $k_i$가 크게 달라집니다.
  - 표준화하지 않으면 이웃 세포 수가 많은 위치에서 단순 합산에 의해 $\sum_j w_{ij} z_{R, j}$의 분산이 왜곡되어 국소 모란 지수가 인위적으로 부풀려집니다.
  - 각 행의 합이 1이 되도록 행-표준화($w_{ij}^* = w_{ij} / \sum_k w_{ik}$)를 적용하면, $\sum_j w_{ij}^* z_{R, j}$는 세포 $i$의 이웃들에서 관측된 수용체 $z$-점수의 '공간적 국소 가중 평균(Spatially Lagged Local Average)'이 됩니다. 이는 단일 세포 드롭아웃(Dropout) 노이즈를 주변 세포들의 평균을 통해 스무딩(Smoothing)하여 통계적 분산 안정성을 보장합니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
전립선암 골전이 환자의 10x Xenium 초고해상도 공간 전사체 데이터에 LIANA+ 공간 이변량 모란 지수 분석을 적용하였다. 분석 결과 두 가지 흥미로운 공간 패턴이 관찰되었다:
1) 골-종양 접촉 경계면(Bone-Tumor Interface)에서 *TNFSF11*(RANKL)과 *TNFRSF11A*(RANK) 쌍이 강력한 High-High 공간 핫스팟($I_i = +0.82, p < 0.001$)을 형성하였다.
2) 종양 중심부 내부에서는 종양 세포의 *CD274*(PD-L1) 발현이 관찰되었으나, CD8+ T 세포의 *PDCD1*(PD-1)과의 공간 상호작용은 High-Low 음의 상관($I_i = -0.45, p < 0.01$)을 보였고, 반면 종양 피막을 둘러싼 myCAF 띠 부위에서 *CXCL12*와 T 세포의 *CXCR4*가 초고밀도 High-High 핫스팟($I_i = +0.76, p < 0.001$)을 형성하고 있었다.

이 환자의 공간 전사체 프로파일이 시사하는 면역 미세환경의 병리학적 상태를 진단하고, 항-PD-1 면역치료 단독 요법이 실패할 수밖에 없는 기전과 이를 극복하기 위한 다학제적 치료 설계를 수립하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **병리학적 상태 진단**:
  - **1. 골-종양 접촉면의 활성 골용해성 파골세포 분화 (Active Osteolytic Niche)**:
    - RANKL-RANK의 High-High 핫스팟은 종양 침윤 전선에서 파골세포가 종양 유래 RANKL에 의해 지속적으로 활성화되어 골 기질을 파괴하고 있음을 시사합니다.
  - **2. 기질 결속에 의한 전형적인 면역 배제형(Immune-Excluded / Stromal Trapping) 종양**:
    - 종양 내부의 PD-L1/PD-1이 High-Low라는 것은 종양 코어 내부에 T 세포가 침투하지 못해 상호작용할 파트너가 없음을 뜻합니다.
    - 반면 외곽 myCAF 띠에서 CXCL12-CXCR4가 High-High 핫스팟을 형성하는 것은, 종양 코어로 향하던 CD8+ T 세포들이 myCAF가 분비하는 강력한 CXCL12 케모카인 장벽에 포획(Trapping)되어 종양 외곽 기질에 물리적으로 억류되어 있을 가능성을 시사합니다(CXCR4는 T 세포 외에 골수계 세포·종양 세포도 발현하므로 공간 동시 국소화만으로 T 세포 포획이 입증되지는 않으며, 위 수치는 가상 시나리오입니다).
- **항-PD-1 단독 요법 실패 기전**:
  - 항-PD-1 단클론 항체가 작동하려면 세포독성 CD8+ T 세포가 종양 세포와 물리적으로 접촉(시냅스 형성)해야 합니다.
  - 그러나 T 세포가 종양 외곽 myCAF 피막에 갇혀 종양 코어 내부의 신생항원에 물리적으로 접근조차 하지 못하므로, PD-1 관문을 차단하더라도 종양 살상 반응이 제한되어 단독 요법의 반응이 낮을 것으로 예상됩니다.
- **다학제적 복합 치료 설계 (탐색적 가설; 확립된 치료가 아니며 Denosumab의 근거는 SRE 지연(Fizazi 2011)·nmCRPC 골전이-무진행생존 연장(Smith 2012)이고 Plerixafor+ICI는 마우스 전임상 근거(Feig 2013)입니다)**:
  - **1단계: 골 파괴 악순환 차단**: 항-RANKL 항체 **데노수맙(Denosumab)**을 투여하여 RANKL-RANK High-High 니치의 골 흡수를 억제하고 골 기질 파괴 및 성장 인자 유리를 줄일 것으로 기대합니다.
  - **2단계: 기질 결속 해제 및 T 세포 진입 유도**: CXCR4 저해제인 **플레릭사포르(Plerixafor)** 또는 **모티틱사이드(Motixafortide)**를 투여하여 myCAF의 CXCL12 장벽에 붙잡힌 T 세포의 결속을 해제, 종양 코어 내부로의 이동(Infiltration)을 촉진합니다.
  - **3단계: 종양 코어 면역 살상 복원**: T 세포가 종양 코어 내부로 침투한 시점에 **항-PD-1/PD-L1 면역관문억제제(Nivolumab 또는 Pembrolizumab)**를 병용 투여하는 접근이 제안될 수 있으나 효과는 입증되지 않았습니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13(1):3224. PMID: 35680885. doi: 10.1038/s41467-022-30755-0
2. Dimitrov D, Schäfer PSL, Farr E, et al. LIANA+ provides an all-in-one framework for cell-cell communication inference. *Nat Cell Biol*. 2024;26(9):1613-1622. PMID: 39223377. doi: 10.1038/s41556-024-01469-w
3. Kolde R, Laur S, Adler P, et al. Robust rank aggregation for gene list integration and meta-analysis. *Bioinformatics*. 2012;28(4):573-80. PMID: 22247279. doi: 10.1093/bioinformatics/btr709
4. Fizazi K, Carducci M, Smith M, et al. Denosumab versus zoledronic acid for treatment of bone metastases in men with castration-resistant prostate cancer: a randomised, double-blind study. *Lancet*. 2011;377(9768):813-22. PMID: 21353695. doi: 10.1016/S0140-6736(10)62344-6
5. Smith MR, Saad F, Coleman R, et al. Denosumab and bone-metastasis-free survival in men with castration-resistant prostate cancer: results of a phase 3, randomised, placebo-controlled trial. *Lancet*. 2012;379(9810):39-46. PMID: 22093187. doi: 10.1016/S0140-6736(11)61226-9
6. Feig C, Jones JO, Kraman M, et al. Targeting CXCL12 from FAP-expressing carcinoma-associated fibroblasts synergizes with anti-PD-L1 immunotherapy in pancreatic cancer. *Proc Natl Acad Sci U S A*. 2013;110(50):20212-7. PMID: 24277834. doi: 10.1073/pnas.1320318110
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- LIANA (liana-py) 1.10.0 (github.com/saezlab/liana-py). 원 논문: Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13(1):3224. PMID: 35680885. doi: 10.1038/s41467-022-30755-0
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
<!-- REFERENCES:END -->

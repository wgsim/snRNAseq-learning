# Chapter 07. 치료 선택압 및 계통 가소성 진화 (Treatment Selection Bottleneck & Lineage Plasticity)

> **Cognitive Goal**: 안드로겐 박탈 요법(Androgen Deprivation Therapy, ADT) 및 차세대 안드로겐 수용체 신호 억제제(ARSI: Enzalutamide, Abiraterone)가 가하는 강력한 진화적 병목(Evolutionary Bottleneck) 하에서 암세포 집단이 겪는 클론 선택(Clonal Selection)의 수리적 동역학을 이해한다. 치료 전 기저에 잠재되어 있던 희귀 내성 서브클론의 선별적 팽창(Selection of Pre-existing Subclones)과 약제 노출에 의한 적응적 후성유전 재프로그래밍(Adaptive De Novo Reprogramming)의 집단유전학적 모델을 규명한다. 특히 *TP53*과 *RB1*의 동시 결손(Co-loss)이 신경 전사 억제자 *REST* 복합체를 해체하고, 신경내분비 파이오니어 전사인자(*ASCL1*, *NEUROD1*, *INSM1*, *SOX2*) 및 *EZH2* 매개 염색질 재구성을 촉발하여 치료 유발 신경내분비 전립선암(t-NEPC)으로의 계통 가소성(Lineage Plasticity)을 완성하는 분자 기전과 그 일부 가역성(EZH2/SOX2 의존)을 체득한다. 최종적으로 단일 세포 유전체 서브클론을 AR/NE 축 기반 4대 표현형 상태(Labrecque 2019, Bluemn 2017 참고)와 매핑하는 전이 역학 프레임워크를 정립한다.

---

## 1. 호르몬 치료 병목(Evolutionary Bottleneck)과 클론 선택 역학

전립선암은 전통적으로 안드로겐 수용체(Androgen Receptor, AR) 신호 전달 경로에 절대적으로 의존하는 호르몬 의존성 종양입니다. 테스토스테론을 고갈시키는 LHRH 작용제/길항제 기반의 1차 ADT 및 AR 리간드 결합을 강력히 저해하는 2차 ARSI(Enzalutamide, Apalutamide, Darolutamide) 치료는 종양 집단에 극단적인 선택압(Selective Pressure)을 가합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               치료 병목(Evolutionary Bottleneck) 하의 클론 선택 다이어그램              │
│                                                                                        │
│   [치료 전: 이질적 다클론성 원발암]                                                    │
│   ● Clone 1 (AR+ 민감 클론): 99% ──► [ADT / Enzalutamide 투여] ──► 대량 세포사멸(PSA 급감)│
│   ▲ Clone 2 (TP53/RB1 결손): 1%  ──► [치료 저항 및 생존]      ──► Clonal Sweeping     │
│                                                                                        │
│   [치료 후: 거세 저항성 진행 (CRPC / t-NEPC)]                                          │
│   ▲ Clone 2의 적응적 팽창 ──► 95% 점유 (AR 발현 소실, 신경내분비 마커 폭발적 발현)    │
│                                                                                        │
│   * 핵심 질문: Clone 2는 치료 전 이미 존재하였는가(Pre-existing),                      │
│                아니면 약제 스트레스가 유도한 신규 획득(De novo)인가?                   │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 상미분방정식(ODE) 기반 클론 경쟁 및 적합도(Fitness) 모델
종양 내에 $K$개의 서브클론이 존재한다고 가정합니다. 각 클론 $i \in \{1, \dots, K\}$의 세포 수를 $N_i(t)$, 총 종양 세포 수를 $N(t) = \sum_{i=1}^K N_i(t)$라 합니다. 약물 농도 $u(t)$ 하에서 각 클론의 증식률을 $r_i(u)$, 사멸률을 $d_i(u)$라 할 때, 클론별 순성장률(Net Growth Rate / Fitness) $s_i(u) = r_i(u) - d_i(u)$는 다음과 같은 미분방정식계를 따릅니다:

$$\frac{dN_i(t)}{dt} = r_i(u(t)) N_i(t) \left( 1 - \frac{N(t)}{K_{\text{carry}}} \right) - d_i(u(t)) N_i(t)$$

(수용 용량 항은 증식에만 곱해집니다. 순성장률 $s_i<0$에 로지스틱 항을 그대로 곱하면 $N\approx K$에서 사멸이 멈추는 비물리적 거동이 생기기 때문입니다. $N\ll K$인 구간에서는 $dN_i/dt\approx s_i N_i$의 지수 근사가 성립하며 아래 Q2는 이 근사를 사용합니다.)

여기서 $K_{\text{carry}}$는 종양 미세환경의 자원(산소, 영양분, 물리적 공간) 제약에 따른 수용 용량(Carrying Capacity)입니다.

- **약제 투여 전 ($u = 0$)**:
  - AR 신호가 온전한 민감 클론 $S$가 가장 높은 적합도를 갖습니다: $s_S(0) > s_R(0)$. 따라서 클론 $S$가 종양의 대부분($>99\%$)을 차지합니다.
- **치료 병목 적용 ($u(t) = u_{\text{enzalutamide}}$)**:
  - 민감 클론 $S$의 사멸률이 증식률을 압도하여 순성장률이 음수로 전락합니다: $s_S(u) \ll 0$.
  - 반면, AR 비의존적 내성 클론 $R$은 약물에 영향을 받지 않고 양의 순성장률을 유지합니다: $s_R(u) > 0$.
  - 클론 $S$가 대량 사멸하면서 가용 자원 공간($K_{\text{carry}} - N(t)$)이 개방되고, 클론 $R$이 급격히 빈 공간을 메우며 종양 전체를 독점하는 **클론 스위핑(Clonal Sweeping)**이 일어납니다.

### 1.2 사전 존재(Pre-existing) vs 신규 획득(De Novo) 판별 모델
집단유전학의 고전적인 **룰리아-델브뤽(Luria-Delbrück) 요동 검정(Fluctuation Test)** 원리는 단일 세포 계통수 해석에 직결됩니다:
1. **Pre-existing Selection (다윈적 선택 모델)**:
   - 치료 전 환자의 단일핵 전사체/유전체 데이터를 초고심도로 탐색했을 때, $0.01\% \sim 0.1\%$ 빈도로 *TP53* 결손 및 *RB1* 결손을 이미 보유한 미세 서브클론이 발견되는 경우입니다. 이 경우 치료는 단지 기존에 존재하던 클론을 체(Sieve)처럼 걸러내어 농축시키는 역할을 합니다.
2. **De Novo Adaptive Plasticity (라마르크적 적응 모델)**:
   - 치료 전에는 유전적으로 균일한 상피세포였으나, 극심한 안드로겐 차단 스트레스가 세포 내 후성유전학적 탈메틸화 및 히스톤 수식을 유도하여 중간 가소성 상태(Plastic intermediate state)를 거쳐 신경내분비 상태로 재프로그래밍되는 경우입니다.
   - 단일 세포 전사체 분석은 이 두 가설이 상호 배타적인 것이 아니라, **"Pre-existing 유전적 취약성(*TP53/RB1* co-loss)을 지닌 서브클론이 약제 선택압 하에서 후성유전학적 De Novo 가소성을 가속화한다"**는 복합 모델을 규명하고 있습니다.

---

## 2. *TP53* 및 *RB1* 동시 결손(Co-loss)과 계통 가소성 유도 기전

정상 전립선 상피 조직에서 세포는 관상세포(Luminal), 기저세포(Basal), 소수의 신경내분비세포(Neuroendocrine)로 분화 운명이 엄격히 고정되어 있습니다. 이러한 계통 정체성(Lineage Fidelity)을 유지하는 핵심 분자적 수호자가 바로 p53과 pRb입니다.

Ku et al. (Science 2017, PMID 28059767)은 *PTEN* 결손 배경의 마우스 모델에서 *RB1*과 *TP53*이 협력하여 계통 가소성·전이·항안드로겐 저항성을 억제함을(두 유전자의 동시 결손에서 가소성과 저항성이 크게 증가함을), Mu et al. (Science 2017, PMID 28059768)은 *TP53*/*RB1* 결손 전립선암에서 *SOX2*가 계통 가소성과 항안드로겐 저항성을 촉진함을 보고하였습니다. 아래 단계별 분자 연쇄(특히 *REST* 복합체 분해)는 이들 논문의 직접 결과가 아니라 후속 문헌을 종합한 *모형적 설명*으로 읽어야 합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│            TP53 / RB1 동시 결손 매개 계통 가소성(Lineage Plasticity) 신호 전달망      │
│                                                                                        │
│   [TP53 손실 (Chr 17p)]  +  [RB1 양측 결손 (Chr 13q)]                                  │
│               │                         │                                              │
│               ▼                         ▼                                              │
│      p53 사멸 회피 / 유전체 불안정   E2F1 전사인자 과활성화                             │
│               └─────────────┬───────────┘                                              │
│                             ▼                                                          │
│   [1단계: 분화 억제자 해체]  REST (RE1-Silencing TF) 복합체 단백 분해 및 전사 억제 해제 │
│                             │                                                          │
│                             ▼                                                          │
│   [2단계: 파이오니어 TF 발현] ASCL1, NEUROD1, INSM1, SOX2 폭발적 전사 개방            │
│                             │                                                          │
│                             ▼                                                          │
│   [3단계: 후성유전학적 침묵] EZH2 (PRC2 복합체) 상향 조절 ──► H3K27me3 침착            │
│                             ├──► AR, KLK3 (관상세포 마커) 프로모터 강한 후성유전 침묵    │
│                             └──► CHGA, SYP, NCAM1 (신경내분비 마커) 고도 전사 활성화   │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 분자 연쇄 반응 단계
1. **E2F1 탈억제 및 REST 축 약화 (모형적 가설)**:
   - *RB1*이 결손되면 E2F 전사인자가 항구적으로 유리되어 세포 주기를 가속할 뿐만 아니라, 신경 유전자 억제자인 **REST (RE1-Silencing Transcription Factor)** 축의 약화와 연결되는 것으로 제안되어 왔습니다(Ku/Mu 2017의 직접 결과는 아님).
2. **신경내분비 마스터 전사인자(Master TFs)의 개방**:
   - REST 억제가 풀리면서, 평소에는 강력하게 응축되어 닫혀 있던 **ASCL1 (Achaete-scute homolog 1)**, **NEUROD1**, **INSM1 (Insulinoma-associated protein 1)**, **SOX2**의 유전자 인핸서/프로모터가 개방(Chromatin accessibility 획득)되고 대량 전사됩니다.
3. **EZH2 활성화 및 루미날 전사체의 후성유전학적 억제**:
   - *TP53/RB1* 결손은 폴리콤 억제 복합체 2(PRC2)의 촉매 효소인 **EZH2 (Enhancer of zeste homolog 2)**의 과발현을 직접 유도합니다.
   - EZH2는 *AR*, *KLK3*(PSA), *NKX3-1* 등 전립선 루미날 분화 마커의 전사 조절 부위에 억제성 히스톤 마크인 **H3K27me3 (Histone H3 Lysine 27 Trimethylation)**을 대량 침착시켜 유전자를 후성유전학적으로 강하게 침묵(Silencing)시킵니다.
   - 결과적으로 종양 세포는 안드로겐 수용체 없이도 생존과 증식이 가능한 소세포 신경내분비암(t-NEPC)으로 전환될 수 있습니다. 이 전환은 EZH2/SOX2 의존적이며 일부 되돌릴 수 있다는 보고가 있습니다(Ku et al. 2017: EZH2 억제제가 AR 발현과 항안드로겐 민감성을 회복; Mu et al. 2017: TP53/RB1 기능 회복 또는 SOX2 억제로 가소성이 역전).

---

## 3. 유전체 서브클론과 4단계 표현형 상태의 시공간적 매핑

전립선암의 4대 표현형 상태(AR-PC/Amphicrine/DNPC/t-NEPC; Labrecque 2019, Bluemn 2017 등의 AR/NE 축 분류를 참고하며 Tang 2022의 CRPC 4군 분류와는 다름)는 아래와 같이 유전체 서브클론과 *경향적으로* 대응합니다. 이 대응은 단순화한 가설이며 아형별 유전형이 필수·보편적이라는 뜻이 아닙니다(예: t-NEPC도 TP53/RB1 양측 결손이 필수·보편은 아님).

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   전립선암 4대 표현형 상태와 유전체 서브클론 매핑 표                   │
│                                                                                        │
│   상태 (State)        대표 전사체 마커               주요 유전체 이상 (Genomic Driver)  │
│   ──────────────────────────────────────────────────────────────────────────────────   │
│   1. AR-PC (Adeno)    AR+, KLK3+, CHGA-, ASCL1-      TMPRSS2-ERG, PTEN loss (Truncal)   │
│   2. Amphicrine       AR+, KLK3+, CHGA+, SYP+        TP53 mut 또는 RB1 mono-loss (전이) │
│   3. DNPC             AR-, KLK3-, NE-                FGF/MAPK 활성(Bluemn), Wnt는 별개 CRPC-WNT 아형(Tang)│
│   4. t-NEPC           AR-, KLK3-, ASCL1++, INSM1++   TP53 / RB1 양측 완전 결손 (Biallelic)│
└────────────────────────────────────────────────────────────────────────────────────────┘
```

1. **AR+ 선암종 (AR-PC, Adenocarcinoma)**:
   - 전형적인 안드로겐 의존적 상태로, 트렁컬 변이(*TMPRSS2-ERG*, *PTEN* 결손)를 지니며 ARSI 약제에 초기 반응을 보입니다.
2. **양성분비 상태 (Amphicrine / Double-Positive)**:
   - 단일 세포 내에서 안드로겐 수용체와 신경내분비 마커(*CHGA*, *SYP*)가 동시에 검출되는 과도기적 전이 상태(Plastic Intermediate)입니다. 주로 *TP53* 변이를 먼저 획득한 후 *RB1* 결손이 진행되는 도중에 포착됩니다.
3. **이중 음성 전립선암 (DNPC, Double-Negative Prostate Cancer)**:
   - AR 발현과 신경내분비 마커가 둘 다 결핍된 상태($\text{AR}^- / \text{NE}^-$)입니다. 안드로겐 신호 대신 FGF/MAPK 경로 활성이 보고되었고(Bluemn 2017), Wnt/$\beta$-catenin 경로는 Tang 2022의 CRPC-WNT 아형(AR 저/음성)에서 논의된 별개의 개념입니다.
4. **치료 유발 신경내분비암 (t-NEPC)**:
   - *TP53* 및 *RB1*의 양측 완전 결손(Biallelic Co-loss)이 완결된 서브클론으로, 극단적인 신경내분비 분화(*ASCL1*, *INSM1*, *DLL3* 과발현)와 함께 안드로겐 신호로부터 100% 해방되어 임상적으로 가장 치명적인 경과를 보입니다.

---

## 4. 가소성 상태 전이 마르코프 모델 (Markov State Transition Model)

시간에 따른 단일 세포의 표현형 상태 전이는 연속 시간 마르코프 연쇄(Continuous-Time Markov Chain, CTMC)로 모델링할 수 있습니다. 4개의 이산적 상태 공간 $\mathcal{S} = \{1: \text{AR-PC}, 2: \text{Amphicrine}, 3: \text{DNPC}, 4: \text{t-NEPC}\}$에 대해, 상태 전이 속도 행렬 $Q(\mathbf{u}, \mathbf{G})$는 약제 농도 벡터 $\mathbf{u}$와 세포의 유전체 유전자형 $\mathbf{G} = (\text{TP53}, \text{RB1})$에 의해 결정됩니다:

$$Q = \begin{pmatrix} 
-q_{12} - q_{13} & q_{12}(\mathbf{u}) & q_{13}(\mathbf{u}) & 0 \\
q_{21} & -q_{21} - q_{24} & 0 & q_{24}(\mathbf{u}, \mathbf{G}) \\
q_{31} & 0 & -q_{31} - q_{34} & q_{34}(\mathbf{u}, \mathbf{G}) \\
0 & 0 & 0 & 0 
\end{pmatrix}$$

- **흡수 상태 가정과 그 한계**:
  - 모형을 단순화하기 위해 $\mathbf{G} = (\text{TP53}^{\text{loss}}, \text{RB1}^{\text{loss}})$인 유전체 배경에서 t-NEPC 상태(State 4)를 흡수 상태(Absorbing State)로 가정할 수 있습니다: 즉, $q_{4j} = 0$ ($j \neq 4$).
  - 그러나 이는 모형상의 가정일 뿐 생물학적 비가역성을 뜻하지 않습니다. Mu et al. (Science 2017)은 *TP53*/*RB1* 기능 회복 또는 *SOX2* 억제로 이 가소성이 역전될 수 있음을, Ku et al. (Science 2017)은 EZH2 억제제가 *AR* 발현과 항안드로겐 민감성을 회복시킴을 보고했습니다. 따라서 현실적인 모형에서는 EZH2/SOX2 의존적인 $q_{4j} > 0$을 허용하는 것이 타당하며, 테스토스테론 보충(Bipolar Androgen Therapy) 후 AR-PC로의 역분화 여부는 별도의 임상·실험 근거가 필요합니다.

---

## 5. Python 실습: 단일 세포 서브클론별 4-State 표현형 시그니처 점수화 파이프라인

`scanpy`를 활용하여 단일핵 전사체 데이터에서 AR-PC, Amphicrine, DNPC, t-NEPC의 유전자 시그니처 점수를 산출하고, CNA 서브클론과의 교차 매핑을 시각화하는 파이프라인입니다.

```python
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def score_prostate_lineage_states(adata: sc.AnnData) -> sc.AnnData:
    """
    단일 세포 전사체 데이터에 대해 전립선암 4대 표현형 상태 시그니처 점수를 계산합니다.
    """
    # 유전자 세트 정의
    gene_sets = {
        "AR_Luminal": ["AR", "KLK3", "KLK2", "TMPRSS2", "NKX3-1", "FOLH1"],
        "Neuroendocrine": ["CHGA", "CHGB", "SYP", "ASCL1", "NEUROD1", "INSM1", "PEG10"],
        "DNPC_Alternative": ["FGFR1", "FGFR2", "FGF2", "WNT5A", "FZD5", "SOX9"],
        "Proliferation_E2F": ["MKI67", "TOP2A", "E2F1", "CCNE1", "EZH2"]
    }
    
    # Scanpy score_genes를 이용한 모듈 점수 산출
    for set_name, genes in gene_sets.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if len(valid_genes) >= 2:
            sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"score_{set_name}")
        else:
            print(f"경고: {set_name} 세트의 검출 유전자 수가 부족합니다.")
            
    # 4대 상태 이산적 분류
    states = []
    for i in range(adata.n_obs):
        s_ar = adata.obs["score_AR_Luminal"].iloc[i]
        s_ne = adata.obs["score_Neuroendocrine"].iloc[i]
        
        if s_ar > 0.1 and s_ne <= 0.0:
            states.append("1_AR-PC")
        elif s_ar > 0.0 and s_ne > 0.1:
            states.append("2_Amphicrine")
        elif s_ar <= 0.0 and s_ne > 0.2:
            states.append("4_t-NEPC")
        elif s_ar <= 0.0 and s_ne <= 0.1:
            states.append("3_DNPC")   # AR-/NE- 명시 조건
        else:
            states.append("Unclassified")   # 경계 영역(예: AR+ 약한 NE)은 억지 분류하지 않음 (임계값 0.1/0.2는 임의 예시)
            
    adata.obs["Lineage_State"] = pd.Categorical(states)
    return adata

# --- 가상 데이터 시뮬레이션 ---
if __name__ == "__main__":
    np.random.seed(42)
    n_cells = 400
    genes = ["AR", "KLK3", "TMPRSS2", "CHGA", "SYP", "ASCL1", "INSM1", "EZH2", "FGFR1", "MKI67"]
    
    # 4개 클론 서브그룹 생성
    # Clone 1: AR-PC (200개)
    # Clone 2: t-NEPC (100개)
    # Clone 3: Amphicrine (50개)
    # Clone 4: DNPC (50개)
    X = np.zeros((n_cells, len(genes)))
    
    # Clone 1 (AR-PC)
    X[:200, :3] = np.random.normal(3.0, 0.5, (200, 3))
    # Clone 2 (t-NEPC)
    X[200:300, 3:8] = np.random.normal(3.5, 0.6, (100, 5))
    # Clone 3 (Amphicrine)
    X[300:350, :5] = np.random.normal(2.0, 0.5, (50, 5))
    # Clone 4 (DNPC)
    X[350:, 8] = np.random.normal(3.0, 0.5, (50,))
    
    # 양수화 및 희소화
    X = np.clip(X, 0, None)
    
    adata = sc.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    
    # CNA 클론 메타데이터 주입
    cna_clones = ["Clone_WT_PTEN"] * 200 + ["Clone_TP53_RB1_loss"] * 100 + ["Clone_TP53_mut"] * 50 + ["Clone_FGF_amp"] * 50
    adata.obs["CNA_Subclone"] = cna_clones
    
    # 계통 점수 계산
    adata = score_prostate_lineage_states(adata)
    
    # 교차 분할표 출력
    cross_tab = pd.crosstab(adata.obs["CNA_Subclone"], adata.obs["Lineage_State"])
    print("\n[CNA 서브클론 vs 계통 표현형 상태 교차표]")
    print(cross_tab)
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
전립선암 세포에서 *TP53* 단독 결손이나 *RB1* 단독 결손만으로는 신경내분비 전립선암(t-NEPC)으로의 계통 가소성(Lineage Plasticity) 전이가 촉발되지 않고, 반드시 두 유전자의 동시 결손(Double Co-loss)이 수반되어야 하는 분자세포생물학적 기전을 신경 억제 복합체 *REST* 및 후성유전학 조절자 *EZH2*와 연계하여 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **상보적 분화 억제 장벽의 유지**:
  - 정상 전립선 상피에서 pRb는 E2F 계열 전사인자를 억제하여 세포 분열을 통제함과 동시에, 신경 유전자 억제자인 **REST (RE1-Silencing Transcription Factor)** 복합체를 안정화하여 신경내분비 프로그램(*ASCL1*, *CHGA* 등)을 물리적으로 봉쇄합니다.
  - 한편 p53은 유전체 손상 감시뿐만 아니라 비정상적인 계통 이탈이나 줄기세포성 역분화 신호가 감지될 때 p21 매개 세포 주기 정지 또는 아포토시스(Apoptosis)를 유도하는 자살 스위치로 기능합니다.
- **단독 결손 시의 안전장치 가동**:
  - *RB1*만 단독으로 결손될 경우: E2F 활성화로 인해 가소성 유도 신호가 발생하지만, 온전한 p53 단백질이 비정상적 증식 스트레스를 감지하여 암세포를 즉각 노화(Senescence)시키거나 사멸(Apoptosis)로 유도하므로 t-NEPC 클론으로 팽창하지 못합니다.
  - *TP53*만 단독으로 결손될 경우: 사멸 장벽은 무너지지만, 온전한 pRb-REST 축이 신경 프로모터를 여전히 단단히 억제하고 있으므로 계통 전환이 일어나지 않습니다.
- **동시 결손(Double Co-loss) 시의 협동적 폭발 기전**:
  - 두 유전자가 동시에 소실되면, p53 사멸 감시망이 마비된 상태에서 pRb 소실에 따른 E2F1 과활성화가 급격히 진행됩니다.
  - 이로 인해 신경 억제자인 **REST 복합체의 기능이 소실**되어 프로모터 억제가 풀리고(이 단계는 모형적 가설), 파이오니어 전사인자인 **ASCL1, NEUROD1, INSM1, SOX2**가 폭발적으로 전사됩니다.
  - 동시에 E2F1은 PRC2 복합체의 촉매 소단위인 **EZH2**의 전사를 대량 유도하여, 전립선 루미날 분화 유전자(*AR*, *KLK3*)의 프로모터에 억제성 **H3K27me3** 마크를 침착시켜 강하게 침묵시킵니다.
  - 결과적으로 세포 사멸 회피와 후성유전학적 계통 전환이 동시에 결합되어 t-NEPC 가소성이 확립됩니다(다만 EZH2/SOX2 의존적으로 일부 가역적일 수 있음).

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
약제 농도 $u$ 하에서 증식하는 민감 암세포 클론 $S$와 내성 클론 $R$의 세포 수를 각각 $N_S(t)$, $N_R(t)$라 하자.
치료 전($t=0$) 종양의 총 세포 수는 $N(0) = N_S(0) + N_R(0) = 10^{10}$개이며, 내성 클론은 $0.01\%$ ($N_R(0) = 10^6$개)의 극미량으로 존재하였다.
약제 투여 시 민감 클론의 순성장률은 $s_S = -0.05 \text{ day}^{-1}$ (지수적 사멸)이고, 내성 클론의 순성장률은 $s_R = +0.02 \text{ day}^{-1}$ (지수적 증식)이다.
1) 종양 내에서 내성 클론의 분율 $f_R(t) = \frac{N_R(t)}{N_S(t) + N_R(t)}$이 $50\%$에 도달하여 임상적 치료 실패(PSA 최저점 통과 후 재발)가 일어나는 시점 $t^*$를 계산하시오.
2) 내성 클론이 사전 존재하는 모델(Pre-existing)과 치료에 의해 매일 민감 세포의 일부가 속도 $\mu = 10^{-6} \text{ day}^{-1}$로 내성 세포로 획득 전환되는 모델(De Novo Plasticity) 간의 분율 도함수 $\frac{df_R}{dt}$ 수식을 비교 유도하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. $50\%$ 도달 시점 $t^*$의 수리적 계산**:
  - $N_S(t) = N_S(0) e^{s_S t} \approx 10^{10} e^{-0.05 t}$
  - $N_R(t) = N_R(0) e^{s_R t} = 10^6 e^{+0.02 t}$
  - 내성 분율 $f_R(t^*) = 0.5$라는 것은 두 클론의 세포 수가 일치함을 의미합니다:
    $$N_S(t^*) = N_R(t^*) \iff 10^{10} e^{-0.05 t^*} = 10^6 e^{+0.02 t^*}$$
    $$\frac{10^{10}}{10^6} = e^{(0.02 - (-0.05)) t^*} \iff 10^4 = e^{0.07 t^*}$$
    $$0.07 t^* = \ln(10^4) = 4 \ln(10) \approx 4 \times 2.3026 = 9.2104$$
    $$t^* = \frac{9.2104}{0.07} \approx 131.58 \text{ 일}$$
  - 즉, 약제 투여 후 약 132일(약 4.4개월) 시점에 내성 클론이 종양의 과반을 점유하게 됩니다.
- **2. 분율 도함수 $\frac{df_R}{dt}$의 유도 및 비교**:
  - **사전 존재 모델 (Pre-existing Selection)**:
    $$f_R = \frac{N_R}{N_S + N_R}$$
    $$\frac{df_R}{dt} = \frac{\dot{N}_R (N_S + N_R) - N_R (\dot{N}_S + \dot{N}_R)}{(N_S + N_R)^2} = \frac{s_R N_R N_{\text{tot}} - N_R (s_S N_S + s_R N_R)}{N_{\text{tot}}^2}$$
    $$\frac{df_R}{dt} = (s_R - s_S) f_R (1 - f_R)$$
    이는 로지스틱 형태의 적합도 차이 $(s_R - s_S)$에 정비례하는 순수 다윈적 선택 방정식입니다.
  - **신규 적응 전환 모델 (De Novo Plasticity)**:
    민감 세포가 속도 $\mu$로 내성 세포로 상태 전이되므로:
    $$\dot{N}_S = s_S N_S - \mu N_S = (s_S - \mu) N_S$$
    $$\dot{N}_R = s_R N_R + \mu N_S$$
    도함수를 전개하면:
    $$\frac{df_R}{dt} = \frac{(s_R N_R + \mu N_S) N_{\text{tot}} - N_R (s_S N_S + s_R N_R)}{N_{\text{tot}}^2}$$
    $$\frac{df_R}{dt} = (s_R - s_S) f_R (1 - f_R) + \mu (1 - f_R)$$
  - **수학적 비교 결론**:
    - 사전 존재 모델에서는 $f_R(0) = 0$이면 $\frac{df_R}{dt} = 0$이 되어 내성 클론이 영원히 발생할 수 없습니다.
    - 반면 De Novo 가소성 항 $\mu(1 - f_R)$이 추가되면, 초기에 내성 클론이 전혀 없더라도($f_R(0)=0$) 초기 발생 속도 $\left.\frac{df_R}{dt}\right|_{t=0} = \mu > 0$가 보장되어 치료 스트레스 자체가 내성의 기폭제로 작동함을 증명합니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
전이성 거세 저항성 전립선암(mCRPC) 환자가 Enzalutamide 치료 9개월 후 간 전이 병변이 급격히 커지며 임상적으로 악화되었다. 혈액 검사상 혈중 PSA는 $0.5 \text{ ng/mL}$로 매우 낮았으나 혈청 뉴런 특이 에놀라아제(NSE)와 크로모그라닌 A(CgA) 수치는 정상 상한치의 15배 이상 폭증하였다.
간 생검 조직의 단일핵 전사체 및 타깃 DNA 시퀀싱 결과:
- 종양 세포의 98%가 Chr 17p (*TP53*) 및 Chr 13q (*RB1*)의 양측 완전 결손(Biallelic loss)을 보유함.
- 전사체 상에서 *AR*, *KLK3* 발현은 0에 수렴하였으며, *ASCL1*, *INSM1*, *EZH2*가 고도 발현됨.

1) 이 환자에서 발생한 질환의 정확한 병리학적/분자생물학적 진명(Diagnosis)을 기술하고, PSA가 정상 수준임에도 간 병변이 급속 진행하는 모순적 임상 양상의 분자적 이유를 설명하시오.
2) 이 환자에게 안드로겐 차단 약제(예: Abiraterone)를 지속하거나 증량 투여하는 것이 무의미한 이유를 진화적 관점에서 기술하고, 현재 임상 지침(NCCN Guideline) 및 최신 표적 치료 기전에 기반하여 즉각 권고해야 할 1차 구제 요법(Salvage Therapy)을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 최종 진단명 및 모순적 임상 양상의 원인**:
  - **최종 진단명**: 치료 유발 신경내분비 전립선암 (Treatment-induced Neuroendocrine Prostate Cancer, **t-NEPC**).
  - **모순적 임상 양상의 원인**:
    - *PSA (Prostate-Specific Antigen, KLK3)*는 안드로겐 수용체(AR)가 전사인자로서 프로모터에 직접 결합해야만 발현되는 단백질입니다.
    - 종양 세포가 *TP53* 및 *RB1* 동시 결손 하에서 EZH2 매개 H3K27me3 침착으로 *AR*과 *KLK3* 유전자를 강하게 침묵시켰기 때문에, 암세포가 체내에 수조 개 증식하더라도 혈중 PSA는 전혀 분비되지 않아 정상치($0.5\text{ ng/mL}$)를 보입니다.
    - 반면 REST 억제 해제로 인해 신경내분비 과립 단백질인 NSE와 Chromogranin A(CgA)의 합성이 폭발적으로 증가하여 혈중으로 분비되며, AR 비의존적 고증식성 소세포 암종의 특성상 간, 내장기 등으로 치명적인 공격적 전이를 일으킵니다.
- **2. 호르몬 치료 무용성 및 구제 치료 전략**:
  - **호르몬 약제 지속/증량의 무용성**:
    - t-NEPC 클론은 AR 신호 경로를 이미 완전히 이탈(Lineage Switch)하였으므로, Abiraterone을 투여하여 안드로겐 합성을 억제하거나 AR을 차단하더라도 종양 세포의 적합도(Fitness)에는 아무런 타격을 주지 못하며 불필요한 약물 독성만을 초래합니다.
  - **권고 구제 치료 전략**:
    - 1) **백금 기반 복합 항암화학요법 (Platinum-based Chemotherapy)**:
      - 소세포폐암(SCLC)과 분자 프로파일(*TP53/RB1* 결손, *ASCL1* 활성)이 동일하므로, **카보플라틴(Carboplatin) 또는 시스플라틴 + 에토포사이드(Etoposide)** 요법(소세포 전립선암/t-NEPC에서 흔히 고려되는 레지멘; 가이드라인 확인 필요)을 고려하여 DNA 손상을 유도하고 급속 분열 세포를 제거해야 합니다.
    - 2) **차세대 표적 치료 임상시험 고려 (DLL3 또는 EZH2 표적)**:
      - ASCL1의 직접 하류 표적인 델타-유사 리간드 3(**DLL3**)를 표적으로 하는 T세포 인게이저(BiTE, 예: Tarlatamab) 또는 후성유전학적 PRC2 억제제(**Tazemetostat**)의 병용은 임상시험·전임상 단계의 *시험적 선택지*로 고려할 수 있습니다(전립선 NEPC에서 확립된 근거는 아직 없음).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Bluemn EG, Coleman IM, Lucas JM, et al. Androgen Receptor Pathway-Independent Prostate Cancer Is Sustained through FGF Signaling. *Cancer Cell*. 2017;32(4):474-489.e6. PMID: 29017058. doi: 10.1016/j.ccell.2017.09.003
2. Labrecque MP, Coleman IM, Brown LG, et al. Molecular profiling stratifies diverse phenotypes of treatment-refractory metastatic castration-resistant prostate cancer. *J Clin Invest*. 2019;129(10):4492-4505. PMID: 31361600. doi: 10.1172/JCI128212
3. Tang F, Xu D, Wang S, et al. Chromatin profiles classify castration-resistant prostate cancers suggesting therapeutic targets. *Science*. 2022;376(6596):eabe1505. PMID: 35617398. doi: 10.1126/science.abe1505
4. Ku SY, Rosario S, Wang Y, et al. Rb1 and Trp53 cooperate to suppress prostate cancer lineage plasticity, metastasis, and antiandrogen resistance. *Science*. 2017;355(6320):78-83. PMID: 28059767. doi: 10.1126/science.aah4199
5. Mu P, Zhang Z, Benelli M, et al. SOX2 promotes lineage plasticity and antiandrogen resistance in TP53- and RB1-deficient prostate cancer. *Science*. 2017;355(6320):84-88. PMID: 28059768. doi: 10.1126/science.aah4307
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
<!-- REFERENCES:END -->

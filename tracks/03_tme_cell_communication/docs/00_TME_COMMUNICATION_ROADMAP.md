# Chapter 00. Track 03 로드맵: 종양 미세환경(TME) 세포 간 신호 전달 네트워크 (Cell-Cell Communication)

> **Cognitive Goal**: 종양 생물학의 패러다임이 세포 자율적(Cell-autonomous) 유전체 변이 분석에서 이질적 세포 간 동적 상호작용 네트워크(Intercellular Communication Network)로 전환된 배경을 이해하고, 전립선암 원발암 및 치명적 거세저항성 뼈 전이(mCRPC Bone Metastasis) 니치에서 작동하는 세포 간 신호 전달 체계의 생물학적 기전과 핵심 공개 벤치마크 데이터셋, 재현 가능한 분석 컴퓨팅 환경을 구축할 수 있다.

---

## 1. 패러다임 전환: 정적 전사체 프로파일링에서 동적 세포 간 신호 네트워크로

지난 10년간 단일세포 RNA 시퀀싱(scRNA-seq / snRNA-seq)은 종양 조직 내 세포의 다양성을 규명하는 강력한 도구로 자리잡았습니다. 그러나 초기 단일세포 연구는 세포 집단을 $t$-SNE나 UMAP 공간에 군집화(Clustering)하고, 세포 유형별 차별 발현 유전자(Marker DEGs)를 나열하는 **"정적 세포 목록 작성(Static Cellular Inventory)"** 수준에 머물렀습니다.

```
[패러다임의 진화: 정적 도감에서 동적 세포 사회학으로]
1세대: 벌크 시퀀싱 (Bulk RNA-seq) ──► 종양 조직 전체의 평균 전사체 (이질성 상실)
2세대: 단일세포 발현 프로파일링   ──► 개별 세포 분류 및 정적 마커 유전자 도감 (세포 간 관계 단절)
3세대: 세포 간 통신망(CCC) 추론   ──► 리간드-수용체 매개 종양-기질-면역 동적 상호작용 생태계 복원
```

### 1.1 세포 자율적(Cell-Autonomous) 분석의 본질적 한계
종양은 악성 상피세포(Malignant Epithelial Cells) 단독으로 증식하는 고립계가 아닙니다. 세포 자율적 분석에만 의존할 경우 다음과 같은 중대한 생물학적 맹점이 발생합니다:
1. **미세환경적 수용성(Microenvironmental Receptivity) 간과**: 암세포가 특정 수용체(예: *MET*, *CXCR4*)를 과발현하더라도, 기질 세포나 면역세포가 상응하는 동족 리간드(HGF, CXCL12)를 충분히 분비하지 않는다면 해당 발암 신호 경로는 실제로 활성화되지 않습니다.
2. **치료 저항성의 외인성(Extrinsic) 기전 추적 불가**: 안드로겐 차단 요법(ADT)이나 방사선 치료에 대한 저항성은 암세포 자체의 유전자 돌연변이뿐만 아니라, 암 연관 섬유아세포(CAF)가 분비하는 *WNT16B*, *SPP1*, *IL6* 등의 외인성 생존 신호에 의해 강력히 유도됩니다.
3. **면역 회피(Immune Evasion) 다이내믹스 왜곡**: 면역관문 분자($PDCD1$, $CTLA4$, $LAG3$)의 mRNA 검출 빈도는 드물고 발현량이 낮아(Dropout 현상), 단일 전사체 마커 비교만으로는 면역 억제적 접촉(Juxtacrine) 및 억제성 사이토카인 구배(Paracrine gradient)의 실제 영향력을 통계적으로 검증하기 어렵습니다.

따라서 단일세포 분석의 핵심 과제는 **"누가 그곳에 있는가(Cell Composition)"**를 넘어, **"누가 누구에게 어떤 생화학적 언어(Ligand-Receptor Pairs)로 지시하고 있는가(Cell-Cell Cross-talk)"**를 정량화하는 계산생물학적 신호망 재구성으로 전환되었습니다.

---

## 2. 전립선암 TME의 생물학적 복잡성: 원발암 vs mCRPC 뼈 전이 니치

전립선암(Prostate Cancer, PCa)은 전형적인 **"면역학적 냉성 종양(Immunologically Cold Tumor)"**으로 분류되며, 진행 단계에 따라 종양 미세환경의 세포 구성과 신호 전달 네트워크가 극적으로 재편됩니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        전립선암 미세환경(TME)의 단계별 재편                              │
└────────────────────────────────────────────────────────────────────────────────────────┘
 [원발성 선암 (Primary Adenocarcinoma)]          [mCRPC 골 전이 니치 (Bone Metastasis)]
 ┌────────────────────────────────────┐          ┌─────────────────────────────────────┐
 │ • 고밀도 반응성 간질 (myCAF 우세)    │          │ • 골아세포-파골세포 악순환 루프      │
 │ • TGF-β 중심 세포외기질(ECM) 경화   │   전이   │ • CXCL12/CXCR4 매개 내막 니치 호밍   │
 │ • CD8+ T세포 기질 침윤 장벽 형성     │ ──────► │ • SPP1+/CD206+ TAM에 의한 강력 면역억제│
 │ • 치밀한 혈관 기저막 구조           │          │ • 골수 지방세포 매개 지질 대사 리프로그래밍 │
 └────────────────────────────────────┘          └─────────────────────────────────────┘
```

### 2.1 원발성 국소 전립선암(Primary PCa) 미세환경
- **반응성 간질(Reactive Stroma)**: 정상 전립선의 평활근(Smooth Muscle) 세포가 소실되고, $\alpha$-SMA($ACTA2$)와 콜라겐을 대량 합성하는 근섬유아세포성 CAF(myCAF)가 종양 선관(Gland) 주위를 포위합니다.
- **물리적 면역 배제(Immune Exclusion)**: 기질 세포가 분비하는 고농도의 $TGF\beta1$ 및 치밀한 피브로넥틴($FN1$)/콜라겐 섬유망은 세포독성 T세포(CD8+ T)가 종양 실질로 침투하지 못하도록 물리적·생화학적 차단벽을 형성합니다.
- **안드로겐 의존적 신호 전달**: 관강 상피세포(Luminal Epithelia)는 안드로겐 수용체(AR) 신호에 종속적이며, 상피-기질 간에는 $FGF2$, $FGF10$, $IGF1$을 통한 항상성 및 증식 신호가 활발합니다.

### 2.2 전이성 거세저항성 전립선암(mCRPC) 뼈 전이 니치: 치명적 악순환(Vicious Cycle)
mCRPC 환자의 다수(보고에 따라 약 80–90%; 전체 전립선암 환자의 비율이 아니며 정확한 수치는 문헌 확인이 필요함)에서 발생하는 뼈 전이는 단순한 물리적 안착이 아니라, 골수(Bone Marrow) 특유의 세포 생태계와의 파괴적인 공모(Collusion)를 수반합니다.
1. **호밍 및 정착 (Homing & Endosteal Niche Seeding)**:
   - 골수 내피세포 및 조골세포가 분비하는 $CXCL12$ 구배를 따라 암세포의 $CXCR4$ 수용체가 결합하여 골내막(Endosteum)에 암세포가 착상합니다.
   - $\alpha_V \beta_3$, $\alpha_5 \beta_1$ 인테그린은 골기질 단백질인 오스테오폰틴($SPP1$), 본 시알로프로틴($IBSP$)과 결합하여 세포 사멸 저항성을 획득합니다.
2. **골 형성/골 흡수 악순환 (The Bone Metastasis Vicious Cycle)**:
   - 전립선암은 전통적으로 조골성(Osteoblastic) 병변을 주로 형성하지만, 실제로는 파골세포(Osteoclast) 매개 골 흡수가 선행하거나 공존합니다.
   - 암세포는 부갑상선호르몬 관련 단백질($PTHrP$), $IL6$ 등을 분비하여 조골세포 계통의 $RANKL$ 발현을 자극할 수 있다고 알려져 있습니다. (*WNT5A*는 이와 반대로 골아세포 니치에서 유래해 ROR2/SIAH2를 통해 전립선암 세포의 휴면을 유도·유지한다는 보고가 있어(Ren et al., J Exp Med 2019, PMID 30593464) 방향과 역할이 상충하므로 여기서 종양 분비 인자로 단정하지 않습니다.)
   - $RANKL$은 파골 전구세포의 $RANK$와 결합하여 파골세포 분화와 골 흡수를 촉발합니다.
   - 용해된 뼈 기질에서 방출된 활성형 $TGF\beta$, $IGF1$, $BMP$는 다시 암세포를 자극하여 암세포 증식, EMT, 추가적인 $PTHrP$ 분비를 유도하는 자가 증폭 악순환 루프를 완성합니다.

---

## 3. 실전 벤치마크 공개 단일세포 전사체 데이터셋 카탈로그

Track 03 전반에서 검증 및 분석에 활용할 수 있는 핵심 공개 데이터셋 3종의 메타데이터와 생물학적 의의를 정리합니다:

| 연구 문헌 (Citation) | 등록 번호 (Accession) | 샘플 및 조직 특성 | 세포 수 (Cells) | 핵심 분석 가치 및 L-R 검증 포인트 |
| :--- | :--- | :--- | :--- | :--- |
| **Kfoury et al.** (*Cancer Cell*, 2021) | **GSE143791** | 전립선암 골전이 환자 9명의 종양·병변 인접 골수·원위 골수 분획 및 양성 골수 대조 7명 | 총 세포 수 미확정 (논문 보고: 샘플당 평균 약 2,490 cells) | 골전이 면역억제 미세환경(T 세포 소진, 골전이 특이 대식세포), $CCL20$-$CCR6$ 축 |
| **Chen et al.** (*Nat Cell Biol*, 2021) | **GSE141445** | 근치적 전립선 절제술 원발 종양 13개 (논문은 별도 11개 샘플 추가 시퀀싱, CRPC 포함) | 36,424 cells (GEO 요약) | 미세환경의 $KLK3$ 이소성 발현, 종양과 활발히 통신하는 내피세포 하위집단(aEC)과 종양 간 세포 간 통신 |
| **Bill et al.** (*Science*, 2023;381:515, PMID 37535729) | 확인된 전립선 accession 없음 (이전 표기 `EGAS00001006653`는 확인 불가) | 다암종(pan-cancer) 인간 종양 대상 연구 — 전립선 3자 아틀라스가 아님 | 확인 불가 (이전 `~65,000`은 근거 없음) | CXCL9:SPP1 대식세포 극성이 암 전반의 세포 프로그램 네트워크를 규정한다는 보고. 전립선 특이 L-R 검증 자료로 쓰려면 별도 전립선 코호트 확인 필요 |

---

## 4. Track 03 커리큘럼 및 4단계 인지 구조 (Cognitive Progression)

Track 03은 단순한 분석 툴 사용법 나열을 지양하고, 생화학적 수용체-리간드 결합 물리화학에서 출발하여 그래프 신경망 및 다중 조건 비교 분석까지 도달하는 나선형 심화 구조로 설계되었습니다.

```
[Level 1: Biophysical Foundation & Cellular Sociology]
  ├── Ch 01. TME 생물학: 이질적 세포 생태계와 리간드-수용체 신호 양식 (Juxtacrine/Paracrine/Autocrine)
  └── Ch 02. 다중 데이터베이스 큐레이션 및 LIANA 합의 순위 집계 (Consensus Rank Aggregation)
          │
[Level 2: Statistical Modeling & Algorithmic Rigor]
  ├── Ch 03. 단일세포 드롭아웃 극복 발현 임계치 설정 및 질량작용 법칙 기반 순열 검정 (Permutation Test)
  └── Ch 04. CellChat 다중 소단위체(Multi-subunit Heteromer) 기하평균 모델링 및 4대 네트워크 중심성
          │
[Level 3: Intracellular Downstream Signaling & Target Validation]
  ├── Ch 05. NicheNet 표적인과 추론: 상류 리간드에서 하류 표적 유전자 조절망(eGRN) 연계
  └── Ch 06. 공간 세포 간 통신: LIANA+ 공간 이변량 지표(Moran's R 등)와 국소 핫스팟
          │
[Level 4: Advanced Comparative Analytics & Translational Discovery]
  ├── Ch 07. 다중 뷰 공간 모델링: MISTy (intra/juxta/para 뷰)
  ├── Ch 08. 전립선암 골전이 미세환경의 5자 신호 전달 회로
  └── Ch 09. 트러블슈팅과 모범 사례
```

---

## 5. 재현 가능한 세포 간 상호작용 분석 환경 구축 (Environment Setup)

세포 간 상호작용 분석 생태계는 Python(`liana`, `scanpy`, `decoupler`; R 패키지명은 decoupleR)과 R(`CellChat`, `NicheNet`)의 최신 패키지들이 융합되어 있습니다. 본 과정에서는 효율성과 대용량 데이터 처리 속도를 극대화하기 위해 Python 기반 통합 플랫폼인 **`LIANA+`** 및 `scanpy`를 중심으로 표준 환경을 구성합니다.

```bash
# 1. TME 세포 간 신호 전달 전용 가상환경 생성 (Python 3.11 이상 권장: decoupler 2.x가 Python>=3.11 요구, liana 1.10은 3.10~3.13 지원)
conda create -n tme_comm_env python=3.11 -y
conda activate tme_comm_env

# 2. 고성능 과학 계산 및 컴파일러 도구 설치
conda install -c conda-forge openblas libgfortran compilers -y

# 3. 단일세포 데이터 분석 코어 라이브러리 설치
pip install scanpy anndata mudata scrublet

# 4. 세포 간 통신망(CCC) 핵심 분석 프레임워크 설치 (LIANA+, decoupler)
pip install liana decoupler

# 5. 네트워크 분석 및 그래프 시각화 유틸리티
pip install networkx igraph leidenalg pyvis statsmodels seaborn matplotlib

# 6. Jupyter 커널 등록
python -m ipykernel install --user --name tme_comm_env --display-name "Python (TME Cell Communication)"
```

### 설치 무결성 검증 스크립트
환경 구성이 정상적으로 완료되었는지 확인하기 위해 아래 검증 스크립트를 실행합니다:

```python
import scanpy as sc
import anndata as ad
import liana as li
import decoupler as dc
import networkx as nx

print(f"[OK] Scanpy Version: {sc.__version__}")
print(f"[OK] AnnData Version: {ad.__version__}")
print(f"[OK] LIANA+ Version: {li.__version__}")
print(f"[OK] DecoupleR Version: {dc.__version__}")
print(f"[OK] NetworkX Version: {nx.__version__}")

# LIANA 기본 내장 합의 데이터베이스(OmniPath Consensus) 로드 테스트
consensus_db = li.rs.select_resource('consensus')
print(f"[OK] OmniPath Consensus L-R Resource Loaded: {len(consensus_db)} consensus pairs available.")
print("=== All TME Communication Core Dependencies Successfully Verified ===")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
단일세포 RNA 시퀀싱 데이터를 기반으로 세포 간 통신(Cell-Cell Communication, CCC)을 추론할 때, 단순히 개별 세포 클러스터의 마커 유전자 차별 발현(DEG)만을 분석하는 접근법이 지닌 생물학적·통계적 한계점 3가지를 제시하고, 리간드-수용체 쌍(L-R Pair) 기반 상호작용 점수 모델링이 이를 어떻게 극복하는지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **단순 DEG 분석의 3대 한계점**:
   - **생물학적 수용성(Receptivity) 결여**: 발신 세포(Sender)가 리간드를 아무리 고발현하더라도, 미세환경 내에 동족 수용체를 발현하는 표적 세포(Receiver)가 없으면 실제 신호 전달이 성립하지 않습니다. DEG 분석은 세포 쌍 간의 '양방향 매칭'을 평가하지 못합니다.
   - **다중 소단위체(Multi-subunit Complexes) 누락**: 인테그린($\alpha_V\beta_3$)이나 인터루킨 수용체($IL6R + IL6ST$)와 같이 2개 이상의 단백질이 물리적으로 결합해야만 작동하는 복합체의 경우, 단일 서브유닛 유전자 하나만 DEG로 검출된다고 해서 기능적 수용체가 형성되었다고 단정할 수 없습니다.
   - **드롭아웃(Dropout)에 의한 위음성 왜곡**: 사이토카인, 성장인자, GPCR 계열 수용체는 기저 발현량이 극도로 낮아 단일세포 데이터에서 0으로 관측되는 비율이 높습니다. 단순 발현량 필터링이나 DEG 순위 산출에서는 이러한 핵심 신호 분자들이 하위 순위로 밀려나 누락됩니다.
2. **L-R 쌍 기반 모델링의 극복 기전**:
   - 발신 세포의 리간드 평균 발현량과 수용 세포의 수용체 평균 발현량을 결합하여 세포 쌍 간의 공동 발현 강도를 정량화합니다. 결합 방식은 도구마다 다릅니다: CellPhoneDB는 두 평균의 **산술평균**(한쪽이 0이면 0), NATMI/Connectome 계열은 곱($\bar{E}_{A,L}\times\bar{E}_{B,R}$), CellChat은 질량작용 법칙 기반 Hill 포화 식입니다.
   - 다중 서브유닛 복합체의 경우 CellPhoneDB는 소단위체의 **최소값**, CellChat은 기하평균(Geometric Mean)을 적용하여 모든 필수 소단위체가 동시에 존재하는지 검증하며, 세포 라벨 순열 검정(Permutation Test)을 통해 우연에 의한 공동 발현을 통계적으로 배제합니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
전이성 거세저항성 전립선암(mCRPC)의 뼈 전이 니치(Bone Metastasis Niche)에서 발생하는 종양세포(Tumor), 골아세포(Osteoblast), 파골세포(Osteoclast) 간의 **'골 전이 악순환(Vicious Cycle)'**을 구성하는 핵심 분자 상호작용 축 4가지를 명시하고, 이 루프가 암세포의 안드로겐 비의존적 생존을 지원하는 기전을 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **악순환을 구성하는 4대 분자 축**:
   - **$PTHrP$ / $IL6$ 축 (Tumor ──► Osteoblast)**: 암세포가 분비한 부갑상선호르몬 관련 단백질($PTHrP$) 및 $IL-6$가 골아세포 표면의 수용체에 결합하여 골아세포의 $RANKL$ 발현을 상향 조절하고, 파골세포 억제 인자인 오스테오프로테게린($OPG$)을 억제합니다.
   - **$RANKL$ / $RANK$ 축 (Osteoblast ──► Osteoclast)**: 골아세포 표면의 $RANKL$이 파골 전구세포의 $RANK$와 결합하여 성숙한 파골세포로의 분화 및 골 기질 흡수(Bone Resorption) 활성을 강력히 촉진합니다.
   - **골 기질 저장 성장인자 방출 ($TGF\beta$, $IGF1$) (Bone Matrix ──► Tumor)**: 파골세포의 산(Acid) 및 단백질 분해 효소 분비로 무기질화된 뼈 기질이 용해되면서, 기질 내에 갇혀 있던 고농도의 잠재형 $TGF\beta$, $IGF1$, $BMP$가 유리되어 활성화됩니다.
   - **$TGF\beta$ / $TGFBR$ 축 (Bone Matrix ──► Tumor Autocrine/Paracrine)**: 활성화된 $TGF\beta$가 전립선암 세포의 $TGFBR1/2$ 복합체에 결합하여 $SMAD2/3$ 전사인자를 활성화하고, 이는 다시 암세포의 침습성, 상피-간엽 전이(EMT), 추가적인 $PTHrP$ 분비를 자극하는 폐쇄 회로(Closed-loop)를 형성합니다.
2. **안드로겐 비의존적 생존 지원 기전**:
   - 지속적으로 방출되는 $IGF1$과 $IL6$는 암세포의 $PI3K/AKT/mTOR$ 경로 및 $JAK/STAT3$ 경로를 교차 활성화(Cross-activation)합니다.
   - 이는 안드로겐 결핍 환경에서도 안드로겐 수용체(AR)의 리간드 결합 영역(LBD)과 무관하게 AR 전사 복합체를 인산화하여 강제 활성화하거나, AR 신호 경로를 우회하는 강력한 세포 사멸 저항성 신호를 지속적으로 공급합니다.

</details>

---

### [문제 3. 임상 사례 분석 (Clinical Interpretation Case Study)]
**[임상 시나리오]**
원발성 전립선 선암 환자군(Cohort A, n=15)과 표준 치료(Enzalutamide + Docetaxel) 후 진행된 mCRPC 뼈 전이 환자군(Cohort B, n=12)의 생검 조직에서 단일세포 전사체 데이터를 수집했습니다. 연구팀은 면역관문 억제제(Anti-PD-1 / Anti-CTLA-4) 병용 요법이 원발암에서는 일부 면역 활성을 유도했으나 mCRPC 뼈 전이 병변에서는 100% 불응성을 보인 원인을 규명하고자 세포 간 신호 전달 네트워크 분석을 기획하고 있습니다.

**[질문]**
1. Cohort A(원발암) 대비 Cohort B(뼈 전이)에서 세포독성 CD8+ T세포의 침윤 및 탈진(Exhaustion)에 관여할 것으로 예상되는 핵심 기질세포(CAF 하위 아형) 및 골수 골수계 세포(TAM/MDSC)의 특이적 리간드-수용체 상호작용 축 2가지를 가설로 제시하시오.
2. mCRPC 뼈 전이 환자에서 면역관문 단독 치료의 실패를 극복하기 위해, 세포 간 신호 전달망 분석 결과를 토대로 제안할 수 있는 혁신적인 '미세환경 재프로그래밍 병용 치료 전략'을 2가지 수립하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **예상되는 핵심 리간드-수용체 상호작용 축 2가지**:
   - **$CXCL12$ - $CXCR4$ 축 (CAF/골수간질세포 ──► CD8+ T세포)**: 뼈 전이 니치의 염증성 섬유아세포(iCAF) 및 골수 간질세포(MSCs)가 고농도의 $CXCL12$를 분비하여 종양 침윤 림프구(TIL)의 종양 세포 접촉을 물리적으로 배제(Exclusion)하거나, T세포 표면의 $CXCR4$ 신호 전달을 통해 림프구의 운동성을 마비시키고 기능적 탈진을 유도합니다.
   - **$SPP1$ - $CD44$ / $\alpha_V\beta_3$ 축 (M2-like TAM ──► CD8+ T세포 및 종양 세포)**: 골수 전이 병변에 풍부한 $SPP1^+$ 골수계 대식세포(TAM)가 오스테오폰틴($SPP1$)을 분비하여 CD8+ T세포의 증식을 직접 억제하고 대식세포 자신의 면역 억제적 M2 분극화를 자가 유지하며, 암세포의 침습 및 생존을 촉진합니다.
2. **미세환경 재프로그래밍 병용 치료 전략 2가지**:
   - **$CXCR4$ 길항제(Plerixafor / Motixafortide) + 면역관문 억제제 병용 (탐색적 가설)**: $CXCL12$ 장벽을 약화시켜 T세포 침투를 허용할 가능성이 제안됩니다. 근거(Feig 2013, 마우스 췌장암 모델)는 전임상이며 전립선 골전이에서 입증된 전략이 아닙니다.
   - **$TGF\beta$ 경로 억제 + $RANKL$ 차단제(Denosumab) 병용 (탐색적 가설)**: Denosumab의 근거는 CRPC 골전이 골격계 사건(SRE) 지연(Fizazi 2011, PMID 21353695)과 nmCRPC 골전이-무진행생존 연장(Smith 2012, PMID 22093187)이며, TGF-β '원천 봉쇄'나 면역관문 억제제 감작 근거는 없습니다. TGF-β 표적 약제 중 bintrafusp alfa는 어느 적응증에서도 승인되지 않았고 INTR@PID BTC 055 시험이 무용성으로 중단되었습니다(Merck KGaA 2021; NCT04066491). Galunisertib의 현재 개발 상태는 확인하지 못했습니다. 따라서 이 병용은 확립된 치료가 아니라 연구 가설입니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Ren D, Dai Y, Yang Q, et al. Wnt5a induces and maintains prostate cancer cells dormancy in bone. *J Exp Med*. 2019;216(2):428-449. PMID: 30593464. doi: 10.1084/jem.20180661
2. Fizazi K, Carducci M, Smith M, et al. Denosumab versus zoledronic acid for treatment of bone metastases in men with castration-resistant prostate cancer: a randomised, double-blind study. *Lancet*. 2011;377(9768):813-22. PMID: 21353695. doi: 10.1016/S0140-6736(10)62344-6
3. Smith MR, Saad F, Coleman R, et al. Denosumab and bone-metastasis-free survival in men with castration-resistant prostate cancer: results of a phase 3, randomised, placebo-controlled trial. *Lancet*. 2012;379(9810):39-46. PMID: 22093187. doi: 10.1016/S0140-6736(11)61226-9
4. Feig C, Jones JO, Kraman M, et al. Targeting CXCL12 from FAP-expressing carcinoma-associated fibroblasts synergizes with anti-PD-L1 immunotherapy in pancreatic cancer. *Proc Natl Acad Sci U S A*. 2013;110(50):20212-7. PMID: 24277834. doi: 10.1073/pnas.1320318110
5. Kfoury Y, Baryawno N, Severe N, et al. Human prostate cancer bone metastases have an actionable immunosuppressive microenvironment. *Cancer Cell*. 2021;39(11):1464-1478.e8. PMID: 34719426. doi: 10.1016/j.ccell.2021.09.005
6. Chen S, Zhu G, Yang Y, et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. *Nat Cell Biol*. 2021;23(1):87-98. PMID: 33420488. doi: 10.1038/s41556-020-00613-6
7. Bill R, Wirapati P, Messemaker M, et al. CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers. *Science*. 2023;381(6657):515-524. PMID: 37535729. doi: 10.1126/science.ade2292
8. EMD Serono Research & Development Institute, Inc. Gemcitabine Plus Cisplatin With or Without Bintrafusp Alfa (M7824) in Participants With 1L BTC. *ClinicalTrials.gov* registry (Phase 2/3; Terminated; 2019-2022). ClinicalTrials.gov: NCT04066491.
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- anndata 0.13.4 (github.com/scverse/anndata). 원 논문: Virshup I, Rybakov S, Theis FJ, et al. anndata: Access and store annotated data matrices. *JOSS*. 2024;9(101):4371. doi: 10.21105/joss.04371
- decoupler 2.2.0 (github.com/scverse/decoupler). 원 논문: Badia-I-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinform Adv*. 2022;2(1):vbac016. PMID: 36699385. doi: 10.1093/bioadv/vbac016
- LIANA (liana-py) 1.10.0 (github.com/saezlab/liana-py). 원 논문: Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun*. 2022;13(1):3224. PMID: 35680885. doi: 10.1038/s41467-022-30755-0
- NetworkX 3.7 (github.com/networkx/networkx).
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- conda (docs.conda.io).
<!-- REFERENCES:END -->

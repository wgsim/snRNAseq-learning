# Chapter 09. 암 유전체 및 클론 진화 분석 트러블슈팅 및 캡스톤 프로젝트 (Troubleshooting & Capstone Project)

> **Cognitive Goal**: 단일핵 전사체(snRNA-seq) 데이터로부터 카피수 변이(CNA) 및 체세포 변이를 추론할 때 발생하는 치명적인 기술적 아티팩트와 편향을 진단하고 해결하는 고급 트러블슈팅 역량을 확립한다. 특히 유전자 고밀도 발현 영역(Chr 6p MHC, Chr 12/17 Keratin, Chr 19q Kallikrein 클러스터)에서 기인하는 가짜 발현 파동(Expression Wave Artifacts)의 수리생물학적 원인을 규명하고 이를 필터링하는 정규화 프로토콜을 체득한다. 비종양 정상 참조 세포(Reference Cells) 선택 편향에 따른 이배체 기저선 드리프트(Baseline Drift)를 교정하고, snRNA-seq의 저조한 시퀀싱 깊이(Coverage) 하에서 체세포 돌연변이(SNV) 호출의 위양성/위음성을 통제하는 베스트 프랙티스를 정립한다. 마지막으로 학습한 전 과정을 집약하여 실제 환자 다중 오믹스 데이터를 분석하는 5단계 캡스톤 프로젝트(Capstone Project) 가이드라인과 100점 만점 평가 루브릭을 제시한다.

---

## 1. 발현 파동 아티팩트(Expression Wave Artifacts) 진단 및 제거

단일 세포 전사체 기반 CNA 추론 도구(inferCNV, CopyKAT 등)는 각 염색체를 따라 유전자 발현량의 이동 평균(Moving Average) 또는 가우스 커널 스무딩(Gaussian Kernel Smoothing)을 적용하여 연속된 유전체 영역의 상대적 전사 활성도를 측정합니다. 그러나 이 과정에서 특정 염색체 영역에 국소적으로 유전자들이 밀집되어 있거나, 생물학적으로 폭발적인 과발현이 일어나는 유전자 군집(Gene Clusters)이 존재할 경우 유전체 복제수가 실제로 변하지 않았음에도 거대한 증폭(Amplification)으로 오인되는 **발현 파동 아티팩트(Expression Wave Artifacts)**가 필연적으로 발생합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   주요 발현 파동 아티팩트 발생 유전자좌 및 생물학적 기전                │
│                                                                                        │
│   게놈 유전자좌     밀집 유전자 군집 (Gene Clusters)       주요 생물학적 활성화 원인    │
│   ──────────────────────────────────────────────────────────────────────────────────   │
│   1. Chr 6p21       HLA 복합체 (HLA-A, HLA-B, HLA-C 등)   인터페론 신호 및 면역 반응   │
│   2. Chr 12q / 17q  케라틴 패밀리 (KRT8, KRT18, KRT19 등)  상피세포 계통 분화 기본 발현 │
│   3. Chr 19q13.3-4  칼리크레인 패밀리 (KLK2, KLK3, KLK4)  안드로겐 수용체(AR) 과활성화 │
│   4. Chr 14q32      면역글로불린 중쇄 (IGH locus)          형질세포/B세포 오염 유입     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 전립선암에서의 결정적 함정: Chr 19q Kallikrein (*KLK*) 클러스터
- **현상**: 전립선암 상피세포에서 염색체 19q13.3-4 영역에 약 300kb에 걸쳐 밀집된 **KLK2, KLK3 (PSA), KLK4** 유전자는 강력한 안드로겐 수용체(AR) 전사 활성화에 의해 상당히 높은 수준으로 발현됩니다(세포 내 RNA 점유율은 검체별로 다르며 일반 수치는 제시하지 않습니다).
- **아티팩트 발생**: 이동 윈도우 스무딩($100 \sim 150$ 유전자 윈도우)을 적용하면, *KLK3* 등의 초고발현 값이 윈도우 안의 인접한 정상 유전자들로 번져나가면서 해당 국소 구간(윈도우 폭 안팎)에 가짜 피크(Pseudo-peak)를 형성할 수 있습니다. 번지는 범위는 윈도우 크기로 제한되며 19q 전체가 아닙니다. 또한 infercnvpy는 평활화 이전에 신호를 $\pm3$(`lfc_clip`)으로 클리핑하므로 영향이 더 줄어듭니다. 예시 산술: 로그 척도에서 +3으로 클리핑된 3개 유전자가 100유전자 윈도우에 있으면 윈도우 평균은 약 $3\times3/100=0.09$ 상승에 그쳐, 복제수 3(상대 $\ln 1.5\approx0.41$)보다 작은 효과입니다. 실제 크기는 발현 수준·참조군·정규화에 따라 달라지므로 KLK 구간의 피크는 '고도 증폭'이 아닌 '국소 위양성 후보'로 보고 마스킹 전후를 비교해 확인하십시오.

### 1.2 아티팩트 필터링 수학적 알고리즘: 블랙리스트 마스킹 및 유전자 밀도 가중 Z-점수
이를 방지하기 위해 윈도우 스무딩 전 단계에서 극단적 이상치 유전자를 제거하는 3단계 정규화 프로토콜을 적용합니다:

1. **사전 지식 기반 유전자 블랙리스트(Blacklist Masking)**:
   - 생물학적 과발현 클러스터(HLA, KRT, KLK, 리보솜 단백질 RPL/RPS, 미토콘드리아 유전자 MT-)를 CNA 계산 입력 대상에서 사전에 마스킹하여 배제합니다.
2. **윈도우 가중 Z-점수(Weighted Z-Score) 정규화**:
   - 각 유전자 $i$에 대해 정상 이배체 참조 세포 집단(Ref)에서의 평균 $\mu_{\text{ref}, i}$과 표준편차 $\sigma_{\text{ref}, i}$를 계산하여 국소 정규화를 수행합니다:

$$\tilde{R}_{c, i} = \frac{R_{c, i} - \mu_{\text{ref}, i}}{\sigma_{\text{ref}, i} + \epsilon}$$

   - 여기서 $\epsilon$은 분모가 0이 되는 것을 방지하는 정규화 상수입니다.
3. **윈도우 절삭 평균(Trimmed Window Mean)**:
   - 윈도우 내 $W$개의 유전자 중 상위/하위 5%의 극단값을 절삭(Trimming)한 후 평균을 산출함으로써, 단일 유전자의 폭발적 발현이 주변 윈도우 전체를 왜곡하는 엣지 파동 효과를 원천 차단합니다.

---

## 2. 정상 참조 세포(Reference Cells) 선택 편향 및 기저선 드리프트 방지

CNA 추론 알고리즘의 성패는 **"진정한 정상 이배체(2N Diad) 기준선"**을 얼마나 신뢰성 있게 설정하느냐에 전적으로 의존합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   정상 참조 세포 선택 편향에 따른 기저선 왜곡 메커니즘                 │
│                                                                                        │
│   [시나리오 A: 편향된 단일 참조군 (오직 T세포만 사용)]                                 │
│   - T세포 특이적 발현 유전자와 상피세포 발현 유전자 간의 거대한 계통적 발현 격차 발생  │
│   ──► 종양 세포 전체에서 수많은 염색체 팔이 가짜 결실/증폭으로 산출됨 (기저선 왜곡)     │
│                                                                                        │
│   [시나리오 B: 권고되는 다중 계통 캘리브레이션 (Multi-Lineage Reference)]             │
│   - 기질세포(섬유아세포/내피세포) + 면역세포(T세포/단핵구) + 정상 전립선 기저/상피세포 │
│   ──► 계통 특이적 전사 노이즈가 상쇄되어 순수한 염색체 복제수 2N 기저선 확립         │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 기저선 드리프트(Baseline Drift) 교정 프로토콜
1. **다중 세포 계통 풀링(Multi-Lineage Pooling)**:
   - 단일 세포 유형(예: T세포만 선택)을 참조군으로 사용할 경우 세포 계통 특이적 유전자 발현 패턴이 복제수 신호로 오인됩니다.
   - 따라서 조직에 실제 존재하는 비악성 계통(기질세포: Fibroblasts·Endothelial cells, 면역세포: T cells·Myeloid cells 등)을 포괄하도록 참조군을 구성하는 것이 권장됩니다(혼합 비율의 정답은 없으며 $1:1$은 예시일 뿐, 비율 변화에 대한 민감도 확인이 필요합니다).
2. **GMM 기반 기저선 센터링(Gaussian Mixture Model Centering)**:
   - 암세포 전체의 스무딩된 복제수 분포에 3-성분 가우스 혼합 모델(Loss, Neutral $2N$, Gain)을 적합(Fitting)하여, 중립 성분의 최빈값(Mode)이 $0.0$ ($\log_2 \text{Ratio}$)에 오도록 전역 기저선을 정렬하는 휴리스틱입니다(중립 성분이 최대 성분이라는 가정이 필요하며, 전장 배가(WGD)·고배수성 종양에서는 틀릴 수 있어 BAF 등 직교 정보로 확인해야 합니다).

---

## 3. snRNA-seq 체세포 SNV 호출 최적화 및 커버리지 극복 프로토콜

단일핵 전사체(snRNA-seq)는 전장 게놈 시퀀싱(WGS)에 비해 유전자 코딩 영역에 리드가 집중되지만, 포획 효율의 한계로 인해 타깃 유전자좌의 시퀀싱 깊이(Depth)가 일반적으로 낮고 유전자별 편차가 매우 큽니다(구체적 수치는 라이브러리·발현량에 따라 달라 일반값을 제시하지 않습니다). 이로 인해 다음과 같은 품질 관리 기준(QC Filters)을 엄격히 준수해야 합니다.

### 3.1 5대 체세포 변이 필터링 기준
1. **염기 품질 점수 (Base Quality $\ge 30$)**:
   - 일루미나 시퀀싱의 광학적 오류 확률이 $0.1\%$ 미만인 염기만 판독에 반영합니다.
2. **매핑 품질 점수 (Mapping Quality, 예시 임계값 $\ge 50$; 정렬기마다 MAPQ 체계가 다름)**:
   - 가위 유전자(Pseudogenes)나 반복 서열(Repeat elements)로 인한 오정렬을 배제하기 위해 고유 정렬 리드만 허용합니다.
3. **독립 UMI 수 (Unique UMI $\ge 2$)**:
   - PCR 증폭 과정 중 최초 사이클에서 발생한 Taq 중합효소 오류(Taq error)가 다수의 리드로 복제되어 단일 세포 돌연변이로 둔갑하는 현상을 방지하기 위해, 반드시 서로 다른 UMI를 가진 최소 2개 이상의 독립 분자가 동일 변이를 지지해야 합니다.
4. **리드 위치 편향 필터 (Read Position Bias Filtering)**:
   - 변이 염기가 리드의 양 끝단(5bp 이내)에만 편중되어 위치하는 경우(정렬기의 엣지 미스매치 아티팩트) 제거합니다.
5. **가닥 편향 필터 (Strand Bias Filtering)**:
   - 전사체의 정방향(Forward) 및 역방향(Reverse) 가닥 중 단 한쪽 가닥에서만 100% 검출되는 산화적 손상 손상 변이(예: 8-oxoG에 의한 $G \to T$ 치환)를 배제합니다.

---

## 4. 암 유전체학 및 클론 진화 실전 캡스톤 프로젝트 루브릭 (5-Stage Capstone Project Rubric)

본 트랙의 전 과정을 통합하여 실제 임상 암 환자의 단일 세포 다중 오믹스 데이터셋을 분석하는 실전 연구 프로젝트를 수행합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   암 유전체학 및 클론 진화 5단계 캡스톤 프로젝트 파이프라인             │
│                                                                                        │
│   [Stage 1] 전처리 및 이배체 정상 참조군 다중 계통 캘리브레이션 (20점)                 │
│        │                                                                               │
│   [Stage 2] 염색체 단위 CNA 추론 및 은닉 마르코프 모델(HMM) 분절화 (20점)              │
│        │                                                                               │
│   [Stage 3] 체세포 SNV / TMPRSS2-ERG 융합 검출 및 듀얼 클론 디컨볼루션 (20점)          │
│        │                                                                               │
│   [Stage 4] NJ / MEDICC2 계통수 복원 및 뮬러 플롯(Muller Plot) 클론 동역학 시각화 (20점)│
│        │                                                                               │
│   [Stage 5] ADT 치료 병목 선택압 및 전이 파종 마이그레이션 모델링 (20점)               │
│        │                                                                               │
│        ▼                                                                               │
│   [최종 판정] 총 100점 만점 / 80점 이상 연구 수행 합격 (Pass)                          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 4.1 정량적 평가 루브릭 (Grading Rubric)

> 아래 배점·80점 합격선은 학습 평가를 위한 *임의 설계값*입니다. 점수와 함께 재현성(코드·환경 고정), 직교 검증(DNA/FISH 등) 수행 여부, 한계 서술 같은 정성 기준을 병행해 평가하는 것을 권장합니다.

| 평가 단계 (Stage) | 핵심 평가 항목 (Criteria) | 배점 | 감점 요인 (Penalties) |
| :--- | :--- | :---: | :--- |
| **Stage 1: 전처리 및 참조군 캘리브레이션** | - QC 통과 단일핵 선별<br>- 다중 계통(면역+기질) 이배체 참조군 구성<br>- 발현 파동(HLA, KLK, KRT) 마스킹 | **20점** | - 단일 면역세포만 참조군으로 사용하여 기저선 왜곡 발생 (-10점)<br>- KLK 클러스터 미마스킹으로 Chr 19q 가짜 증폭 보고 (-10점) |
| **Stage 2: CNA 추론 및 HMM 분절화** | - 염색체별 스무딩 및 기저선 센터링<br>- HMM 기반 복제수 상태 산출(inferCNV i6 모델: 상대 복제수 {0,0.5,1,1.5,2,3}, 즉 CN 0·1·2·3·4·6의 6상태)<br>- *PTEN*(10q), *NKX3-1*(8p) 결손 정확한 감지 | **20점** | - 과도한 스무딩으로 국소 드라이버 결손 소실 (-8점)<br>- HMM 전이 확률 튜닝 실패로 인한 과다 분절화 (-10점) |
| **Stage 3: 체세포 변이 및 융합 검출** | - `CB`/`UB` 기반 *TMPRSS2-ERG* 접합부 리드 식별<br>- *TP53* R273H, *SPOP* F133V 단일 세포 VAF 계산<br>- 베타-이항 과포산 유전자형 호출 | **20점** | - PCR 중복 미제거로 단일 UMI를 다세포 변이로 과대평가 (-10점)<br>- 융합 유전자 위양성 필터 미적용 (-8점) |
| **Stage 4: 계통수 복원 및 뮬러 플롯** | - 절단점 거리 또는 `MEDICC2` 최소 사건 거리 산출<br>- Saitou-Nei NJ 또는 최소 진화 계통수 완성<br>- 시간/치료 경과별 뮬러 플롯 클론 팽창 가시화 | **20점** | - 단순 $L_2$ 유클리드 거리 사용으로 절단점·사건 기반 구조를 반영하지 못한 계통수(비균일 진화 속도로 인한 위상 오류 포함) (-10점)<br>- 비간명적(Non-parsimonious) 위상 오류 (-10점) |
| **Stage 5: 치료 선택압 및 전이 모델링** | - ADT/ARSI 병목 하 클론 스위핑 역학 ODE 적합<br>- *TP53/RB1* co-loss 매개 t-NEPC 가소성 매핑<br>- 전이 간(Met-to-Met) 마이그레이션 그래프 도출 | **20점** | - t-NEPC 표현형을 새로운 원발암으로 오분류 (-10점)<br>- 다클론 파종 증거 무시 및 단순 선형 전이 단정 (-10점) |
| **합계 (Total)** | **전체 5단계 종합 연구 결과 보고서** | **100점** | **80점 이상: 이수 (Pass) / 80점 미만: 재수행 (Fail)** |

---

## 5. Python 실습: 발현 파동 아티팩트 마스킹 및 기저선 자동 캘리브레이션 파이프라인

단일 세포 발현 행렬에서 Chr 19q *KLK* 및 Chr 6p *HLA* 클러스터를 탐색하여 마스킹하고, 가우스 혼합 모델(GMM)을 통해 복제수 기저선(Baseline)을 $0.0$으로 자동 보정하는 파이프라인입니다.

```python
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

def filter_expression_wave_artifacts(
    expr_df: pd.DataFrame, 
    gene_annotation: pd.DataFrame,
    blacklist_clusters: list = None
) -> pd.DataFrame:
    """
    발현 파동 아티팩트를 유발하는 유전자 군집(HLA, KLK, KRT 등)을 마스킹합니다.
    
    Parameters:
    -----------
    expr_df : pd.DataFrame
        세포(행) x 유전자(열) 정규화 발현 행렬
    gene_annotation : pd.DataFrame
        유전자별 염색체(chromosome) 및 물리적 위치(start, end) 정보
    blacklist_clusters : list
        배제할 유전자 심볼 리스트 (예: KLK2, KLK3, HLA-A, KRT8 등)
    """
    if blacklist_clusters is None:
        blacklist_clusters = [
            "KLK2", "KLK3", "KLK4",
            "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
            "KRT8", "KRT18", "KRT19"
        ]
        
    genes_to_mask = [g for g in blacklist_clusters if g in expr_df.columns]
    print(f"[아티팩트 필터링] 마스킹 대상 유전자 {len(genes_to_mask)}개 배제 완료: {genes_to_mask[:5]}...")
    
    # 해당 유전자 제외
    filtered_df = expr_df.drop(columns=genes_to_mask)
    return filtered_df

def calibrate_diploid_baseline_gmm(
    cna_signal_df: pd.DataFrame, 
    n_components: int = 3
) -> pd.DataFrame:
    """
    Gaussian Mixture Model (GMM)을 적용하여 복제수 중립(2N) 상태의 중심점을 0.0으로 캘리브레이션합니다.
    """
    # 전체 세포 및 유전자 윈도우의 신호를 1차원 배열로 평탄화
    flat_signals = cna_signal_df.values.flatten().reshape(-1, 1)
    
    # 3-성분 GMM 적합 (Loss, Neutral 2N, Gain)
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(flat_signals)
    
    # 가장 가중치(비율)가 높은 성분을 정상 이배체(2N) 기저선으로 추정
    diploid_comp_idx = np.argmax(gmm.weights_)
    diploid_baseline_shift = gmm.means_[diploid_comp_idx][0]
    
    print(f"[기저선 캘리브레이션] 추정된 2N 기저선 오프셋: {diploid_baseline_shift:.4f}")
    
    # 오프셋 차감으로 센터링
    calibrated_df = cna_signal_df - diploid_baseline_shift
    return calibrated_df

# --- 실행 시뮬레이션 ---
if __name__ == "__main__":
    np.random.seed(42)
    genes = [f"Gene_{i}" for i in range(50)] + ["KLK2", "KLK3", "HLA-A", "HLA-B"]
    cells = [f"Cell_{j}" for j in range(100)]
    
    # 배경 발현 데이터 생성
    data = np.random.normal(1.0, 0.3, (100, len(genes)))
    # KLK 및 HLA 유전자에 극단적인 20배 폭발 발현 주입 (아티팩트 시뮬레이션)
    data[:, -4:] = np.random.normal(20.0, 2.0, (100, 4))
    
    df_expr = pd.DataFrame(data, index=cells, columns=genes)
    
    # 1. 아티팩트 유전자 마스킹
    df_clean = filter_expression_wave_artifacts(df_expr, gene_annotation=None)
    
    # 2. 기저선 드리프트 보정
    # 기저선이 1.2로 치우쳐진 신호 시뮬레이션
    mock_cna_signal = df_clean + 0.25
    df_calibrated = calibrate_diploid_baseline_gmm(mock_cna_signal)
    
    print("\n[캘리브레이션 완료 후 요약 통계]")
    print(df_calibrated.iloc[:5, :5])
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
단일 세포 전사체 데이터에서 카피수 변이(CNA) 분석을 수행할 때, 전립선암 조직에서 염색체 19q에 위치한 칼리크레인(*KLK*) 유전자 군집과 염색체 6p에 위치한 주조직적합성복합체(*HLA*) 유전자 군집이 일으키는 '발현 파동 아티팩트(Expression Wave Artifacts)'의 수리생물학적 원인과, 이를 단일 세포 알고리즘 차원에서 사전 차단하기 위한 최적의 계산적 해결책을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **수리생물학적 발생 원인**:
  - inferCNV, CopyKAT 등의 도구는 염색체 물리적 좌표 순서대로 정렬된 수백 개의 유전자 발현량을 이동 윈도우(Moving average / Gaussian window)로 스무딩하여 카피수 신호를 추론합니다.
  - 전립선 상피암세포에서 *KLK2, KLK3(PSA), KLK4*는 안드로겐 수용체(AR)에 의해, *HLA-A, HLA-B, HLA-C*는 종양 미세환경의 인터페론 신호에 의해 인접한 일반 유전자들보다 크게 과발현될 수 있습니다.
  - 윈도우 스무딩 연산 시 이 극단적인 국소 이상치(Outlier) 값이 윈도우 너비(예: 100 유전자)만큼 좌우로 번져나가면서, 실제로는 정상 복제수($2N$)를 유지하고 있는 19q 또는 6p의 국소 구간(윈도우 폭 안팎)의 평균 발현 비율을 인위적으로 끌어올려 가짜 파동(False Positive Peak)을 형성할 수 있습니다(번지는 범위는 윈도우 크기로 제한되며, 염색체 팔 전체 수준의 효과는 아닙니다).
- **최적의 계산적 해결책**:
  - **1. 사전 지식 기반 블랙리스트 마스킹(Blacklist Masking)**: 윈도우 스무딩을 수행하기 전 단계에서 *KLK* 패밀리, *HLA* 복합체, *KRT*(케라틴), 리보솜 단백질 유전자들을 계산 대상 목록에서 영구 배제(Masking out)합니다.
  - **2. 윈도우 절삭 평균(Trimmed Mean) 및 유전자별 Z-점수 정규화**: 윈도우 내 상하위 극단값을 잘라내고, 각 유전자 발현량을 이배체 정상 참조 세포의 분산($\sigma_{\text{ref}}$)으로 나눈 표준화 점수를 활용하여 단일 유전자의 과발현이 전체 윈도우에 미치는 전파 왜곡을 원천 차단합니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
정상 참조 세포 집단(Reference)의 선택이 편향되어(특히 *영역 특이적* 참조 편향이 실제 위험이며, 02장의 세포별 median 중심화는 아래의 *전역 오프셋*을 일부 제거하므로 이 시나리오는 중심화가 불완전한 경우를 가정), 참조군의 실제 배수성이 이배체($2N$)가 아닌 림프구 집단의 전사 억제나 측정 오류로 인해 $\theta_{\text{ref}} = 2N - \delta$ ($\delta > 0$)로 과소평가되었다고 가정하자.
종양 세포 $c$의 실제 게놈 복제수가 중립($C_{c, k} = 2N$)인 영역 $k$에서, 추론된 로그 복제수 비율 $R_{c, k} = \log_2 \left( \frac{C_{c, k}}{\theta_{\text{ref}}} \right)$의 드리프트 편향을 수학적으로 유도하고, 이 왜곡이 HMM(은닉 마르코프 모델)의 유전자 상태 디코딩에서 전체 게놈 단위의 광범위한 위양성 증폭(False Positive Genome-wide Amplification)을 유발하는 메커니즘을 전이 확률 관점에서 증명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 로그 복제수 비율 $R_{c, k}$의 기저선 드리프트 편향 유도**:
  - 참된 정상 이배체 기저선을 $2N$이라 할 때, 편향된 참조군의 기준치는 $\theta_{\text{ref}} = 2N - \delta = 2N(1 - \frac{\delta}{2N})$입니다.
  - 종양 세포의 실제 복제수가 완전히 정상($C_{c, k} = 2N$)인 영역에서 관측되는 추정 비율 $R_{c, k}$는:
    $$R_{c, k} = \log_2 \left( \frac{2N}{2N(1 - \frac{\delta}{2N})} \right) = -\log_2 \left( 1 - \frac{\delta}{2N} \right)$$
  - 테일러 전개($\ln(1 - x) \approx -x$)를 적용하면:
    $$R_{c, k} \approx \frac{1}{\ln 2} \cdot \frac{\delta}{2N} > 0$$
  - 즉, 실제로는 아무런 카피수 변이가 없는 모든 정상 염색체 영역의 관측 신호가 $+ \Delta R$만큼 전역적으로 상향 편향(Positive Baseline Drift)됩니다.
- **2. HMM 디코딩에서의 위양성 증폭 유발 메커니즘 증명**:
  - HMM의 은닉 상태를 $S_k \in \{0: \text{Loss}, 1: \text{Neutral } 2N, 2: \text{Gain}\}$이라 하자.
  - 방출 확률(Emission Probability)은 가우스 분포 $P(R_k \mid S_k = s) = \frac{1}{\sqrt{2\pi}\sigma} \exp\left( -\frac{(R_k - \mu_s)^2}{2\sigma^2} \right)$를 따르며, 이론적 중심값은 $\mu_1 = 0.0$, $\mu_2 = \log_2(3/2) \approx 0.585$입니다.
  - 전역 기저선 오프셋 $\Delta R = \frac{\delta}{2N \ln 2}$가 발생하면, 중립 상태($S_k=1$)의 방출 확률 로그 우도는:
    $$\ln P(R_k \mid S_k = 1) = -\frac{(\Delta R - 0)^2}{2\sigma^2} + C$$
    반면 증폭 상태($S_k=2$)의 방출 확률 로그 우도는:
    $$\ln P(R_k \mid S_k = 2) = -\frac{(\Delta R - \mu_2)^2}{2\sigma^2} + C$$
  - 두 우도의 차이를 계산하면:
    $$\ln \frac{P(R_k \mid S_k = 2)}{P(R_k \mid S_k = 1)} = \frac{\mu_2}{\sigma^2} \left( \Delta R - \frac{\mu_2}{2} \right)$$
  - 편향 $\Delta R$이 임계값 $\frac{\mu_2}{2} \approx 0.29$를 초과하거나, HMM의 전이 지속성(Self-transition probability $T_{22} \gg T_{12}$)과 결합되면, 비터비(Viterbi) 알고리즘은 실제 정상인 염색체 전체를 증폭 상태($S_k = 2$)로 디코딩하게 됩니다.
  - 결과적으로 환자의 암세포 전체 게놈이 광범위하게 증폭된 것으로 잘못 판정되는 치명적인 전역 위양성이 발생합니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
다발성 골 전이를 동반한 전립선암 환자의 전립선 절제 검체에서 8,000개의 단일핵(snRNA-seq)을 분석한 연구원이 다음과 같은 문제에 봉착하여 프로젝트 심사(Rubric 평가)에서 불합격 판정(65점)을 받았다:
1) Chr 19q 전체가 모든 암세포에서 8-copy 증폭으로 보고되었음.
2) 환자의 골 미세환경 T세포 500개만을 유일한 정상 참조군으로 지정하여 분석한 결과, 종양 세포뿐만 아니라 검체 내 정상 혈관 내피세포마저 염색체 8번과 10번 전체가 결실된 것으로 보고됨.
3) *TP53* 점돌연변이(R273H) 분석에서 시퀀싱 리드가 단 1개만 매핑된 세포들을 모두 돌연변이 양성 세포로 판정하여 환자 종양의 95%가 *TP53* 변이 클론이라고 결론지었으나, 벌크 DNA 전장 엑솜 시퀀싱(WES) 결과 해당 변이의 VAF는 12%에 불과하였음.

본 장의 트러블슈팅 가이드라인 및 캡스톤 프로젝트 루브릭에 근거하여, 위 3가지 결함의 근본 원인을 진단하고 점수를 만점(100점)으로 개선하기 위해 즉각 적용해야 할 수정 파이프라인 단계를 상세히 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. Chr 19q 증폭 아티팩트의 진단 및 교정**:
  - **원인 진단**: 안드로겐 수용체(AR)에 의해 극단적으로 과발현된 *KLK2, KLK3(PSA), KLK4* 유전자가 윈도우 스무딩 과정에서 주변 19q 구간으로 발현 파동(Expression Wave)을 전파시킨 전형적인 위양성 피크 후보입니다.
  - **교정 단계**: 윈도우 스무딩 전처리 단계에서 *KLK* 유전자 클러스터를 블랙리스트로 지정하여 마스킹 배제하고, 절삭 평균(Trimmed mean) 윈도우 필터를 적용하여 Chr 19q의 참된 카피수 기저선을 재산출합니다.
- **2. 정상 참조군 선택 편향 및 기저선 드리프트의 진단 및 교정**:
  - **원인 진단**: 오직 활성화된 T세포만을 단일 참조군으로 사용함으로써 면역세포 특이적 전사체 패턴이 상피/내피세포의 기준선을 심각하게 왜곡시켰습니다. 정상 내피세포마저 Chr 8, 10 결실로 나타난 것은 기저선 드리프트에 의해 정상 이배체($2N$)가 결실 상태로 오분류되었음을 강하게 시사하는 단서입니다.
  - **교정 단계**: T세포 외에 검체 내의 정상 섬유아세포(Fibroblasts)와 정상 혈관 내피세포(Endothelial cells)를 동등한 비율로 병합한 **다중 계통 참조군(Multi-lineage reference pool)**을 재구성하고, GMM 센터링 알고리즘을 실행하여 비종양 정상 세포들의 복제수 상대 비율이 정확히 $0.0$ ($2N$)에 정렬되도록 기저선을 자동 보정합니다.
- **3. 단일 세포 SNV 오분류 및 과대평가의 진단 및 교정**:
  - **원인 진단**: snRNA-seq의 PCR 증폭 오류(Taq polymerase error) 및 시퀀싱 에러를 고려하지 않고, 단 1개의 UMI 리드($\text{Depth}=1, k^{\text{alt}}=1$)만을 지지 근거로 돌연변이를 확정했습니다. 이로 인해 수많은 야생형 세포가 위양성 변이 세포로 둔갑하여 WES 결과(12%)와 8배 이상의 불일치를 초래했습니다.
  - **교정 단계**:
    - 최소 2개 이상의 독립적인 고유 분자 식별자(Unique UMI $\ge 2$)가 변이를 지지하는 엄격한 분자 중복 필터를 적용합니다.
    - 염기 품질(BQ $\ge 30$) 및 매핑 품질(MQ $\ge 50$) 필터를 적용한 후, 본문에서 유도한 **베타-이항 혼합 모델(Beta-Binomial Mixture Model)**을 구동하여 사후 확률 $P(G_c = \text{Mut} \mid \text{Data}) \ge 0.85$인 세포만을 변이 클론으로 판정함으로써 벌크 WES VAF(12%)와 대략적으로 일치하는지로 교차 검증합니다(snRNA 기반 변이 양성 세포 분율은 대립유전자 특이 발현·NMD·발현 수준에 의존해 DNA VAF와 정확히 일치할 수 없으며, 순도·카피수 보정 후 세포분율 ≈ 2×VAF/순도 등으로 비교해야 합니다).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Saitou N, Nei M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol Biol Evol*. 1987;4(4):406-25. PMID: 3447015. doi: 10.1093/oxfordjournals.molbev.a040454
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- scikit-learn 1.9.1 (github.com/scikit-learn/scikit-learn).
<!-- REFERENCES:END -->

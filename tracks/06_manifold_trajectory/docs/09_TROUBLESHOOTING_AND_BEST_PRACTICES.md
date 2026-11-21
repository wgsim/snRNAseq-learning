# Chapter 09. 궤적 분석 문제 해결, 품질 진단 및 캡스톤 루브릭 (Trajectory Troubleshooting, Diagnostics & Capstone Rubric)

> **Cognitive Goal**:
> - 단일 세포 RNA velocity 및 매니폴드 동역학 궤적 추론에서 빈번히 마주치는 3대 핵심 함정(다중 동역학 체계 Multiple kinetics, 속도-단백질 불일치 Velocity-protein discordance, 루트/싱크 세포 오지정 Root/sink misassignment)의 수학적 및 생물학적 근본 원인을 진단하고, 엄밀한 디버깅 해결책을 수립한다.
> - 속도 신뢰도(Velocity confidence), 그래프 일관성(Graph coherence), 경계 교차 정확도(Cross-boundary correctness) 등 정량적 품질 진단 지표(Diagnostic metrics)의 수리적 정의와 임계값 판정 기준을 정립한다.
> - 배치 효과(Batch effect) 보정 시 세포 매니폴드의 기하학적 연속성과 스플라이싱 물리 항등식($\dot{s} = \beta u - \gamma s$)이 파괴되지 않도록 보호하는 실무 의사결정 트리(Decision Trees)를 체계화한다.
> - 전이성 거세 저항성 전립선암(mCRPC)의 계통 가소성(Lineage Plasticity) 전이를 완벽히 규명하기 위한 5단계 엔드-투-엔드 매니폴드 캡스톤 프로젝트 루브릭(5-Stage Manifold & Dynamic Trajectory Capstone Project Rubric)을 정립하고, Track 01(Multiome), Track 02(Spatial), Track 04(Genomic subclones)와의 교차 검증 파이프라인을 완성한다.

---

## 1. 궤적 분석의 3대 핵심 함정(Pitfalls)과 근본 원인 진단

단일 세포 동역학 분석은 정적인 스냅샷(Static snapshot) 데이터로부터 시간에 따른 변화율을 역추정하는 고도의 역문제(Inverse problem)입니다. 따라서 전사체 생물학과 수학적 모델의 가정이 불일치할 때 심각한 아티팩트가 발생합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   궤적 분석의 3대 주요 함정과 수학적/생물학적 해결책                    │
│                                                                                        │
│   [함정 1: 다중 동역학 체계 (Multiple Kinetics)]                                       │
│   ● 현상: 서로 다른 세포 계통에서 전사율/분해율 상이 ──► 왜곡된 위상 곡선 적합        │
│   ● 처방: 유전자별 다중 동역학 적합 또는 CytoTRACE/Pseudotime 커널과의 융합           │
│                                                                                        │
│   [함정 2: 속도-단백질 불일치 (Velocity-Protein Discordance)]                          │
│   ● 현상: mRNA 속도는 급증하나 번역 지연으로 표현형 전환 미발생                       │
│   ● 처방: CITE-seq 표면 단백질 및 다중 오믹스 통합 검증                               │
│                                                                                        │
│   [함정 3: 루트-싱크 역전 (Root-Sink Inversion)]                                       │
│   ● 현상: 약제 스트레스(Enzalutamide) 시 전사 급감 ──► (경우에 따라) 거짓 역류 (NEPC ──► Adeno) │
│   ● 처방: CellRank 2 G-PCCA 슈어 분해 및 거시상태 전이 행렬 P̃를 통한 초기/말단 재지정 │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 함정 1: 다중 동역학 체계 (Multiple Kinetic Regimes)
- **근본 원인**:
  - 표준 scVelo 모델은 전체 세포 집단에 걸쳐 유전자당 하나의 스플라이싱 속도 상수 $\beta$와 분해 속도 상수 $\gamma$가 고정되어 있다고 가정합니다 ($\frac{ds}{dt} = \beta u - \gamma s$).
  - 그러나 mCRPC에서 루미날 선암종 세포와 소세포 신경내분비암(t-NEPC) 세포는 서로 다른 크로마틴 환경과 전사 인자 풀을 가집니다. 예를 들어 세포 상태에 따라 성숙 RNA 분해율 $\gamma$가 다를 수 있다는 것은 가설이며(분열 주기 속도와 $\gamma$의 직접 연관을 보인 근거는 제시하지 못했습니다), 이 경우 단일 모델 적합이 어긋납니다.
  - 단일 위상 타원 곡선(Phase portrait)을 강제로 적합하면, 분해율이 다른 두 클러스터가 마치 단일 분화 경로의 상승기와 하강기인 것처럼 오인되어 궤적이 왜곡됩니다.
- **해결 프로토콜**:
  - `scv.tl.recover_dynamics` 실행 시 계통별 사전 확률(Priors)을 반영하거나, 유전자 공유 잠재 시간(Gene-shared latent time)을 활용합니다.
  - 계통 간 동역학 상수가 심하게 분열된 유전자는 속도 계산 피처에서 제외하고, `CellRank 2`의 `CytoTRACEKernel` 또는 `PseudotimeKernel`과의 볼록 결합을 통해 위상 왜곡을 상쇄합니다.

### 1.2 함정 2: 속도-단백질 불일치 (Velocity-Protein Discordance)
- **근본 원인**:
  - RNA velocity 벡터 $\mathbf{v}$는 단지 수 분~수십 분 단위의 즉각적인 mRNA 변화율일 뿐이며, 단백질 번역(Translation), 번역 후 변형(PTM), 단백질 분해(Degradation) 과정의 지연(Lag time, 수 시간~수일)을 반영하지 못합니다.
  - 특히 전립선암의 계통 가소성 스위치는 REST 기능 상실(REST4 대체 스플라이싱 등)과 EZH2에 의한 H3K27me3 침착 등 후성유전/RNA 처리 수준의 이벤트와 연관되는 것으로 보고됩니다.
- **해결 프로토콜**:
  - CITE-seq 또는 다중 오믹스(snATAC-seq) 데이터를 연계하여 전사체 속도와 염색질 접근성/단백질 수준의 일치도를 교차 검증합니다.

### 1.3 함정 3: 루트/싱크 역전 (Root-Sink Inversion)
- **근본 원인**:
  - Enzalutamide 처리로 AR 표적 유전자의 전사가 차단되면 해당 루미날 유전자의 속도는 음수($v<0$, 발현 감소)가 됩니다. 이는 루미날 상태에서 멀어지는 정방향 신호이지만, 하향 조절 유전자가 투영 방향을 지배하거나 $\gamma$·$u/s$ 추정이 틀리면 가소성 세포가 루미날 선암종으로 역류하는 것처럼 보이는 인위적 벡터장이 형성될 수 있습니다.
- **해결 프로토콜**:
  - 단순 국소 속도 벡터의 흐름에만 의존하지 않고, G-PCCA 실 슈어 분해(Real Schur Decomposition)를 통해 거시상태의 자체 전이 확률 $\tilde{P}_{ss}$을 검토합니다.
  - 발원지(Source)와 흡수 종착지(Sink)를 생물학적 마커(*KRT8/18* vs *CHGA/SYP*) 및 CytoTRACE 발달 잠재력으로 앵커링(Anchoring; AR+/Amphicrine/DNPC/t-NEPC 4상태 분류의 출처는 07장 5절 참조)하여 방향성을 강제 정렬합니다.

---

## 2. 정량적 품질 진단 지표 (Diagnostic Metrics)

궤적 추론 결과의 신뢰성을 평가하기 위한 정량적 지표입니다. scVelo는 `velocity_confidence`를, CellRank 2는 `Kernel.cbc`를 제공하며, 아래 2.2·2.3의 식은 이를 참고한 교육용 자체 정의 지표(두 패키지에 같은 이름의 함수 없음)입니다.

### 2.1 속도 신뢰도 (Velocity Confidence: $C_{\text{vel}}$)
세포 $i$의 속도 벡터 $\mathbf{v}_i$가 국소 이웃 세포들 $\mathcal{N}(i)$의 속도 벡터들과 얼마나 일관된 방향성을 갖는지 측정합니다. (scVelo의 `velocity_confidence`는 각 속도 벡터에서 유전자 평균을 제거한 뒤 이웃과의 Pearson 상관의 평균을 구하고 0 미만은 0으로 절삭합니다. 아래 코사인 식은 이를 단순화한 표기입니다.)

$$C_{\text{vel}}(i) = \frac{1}{|\mathcal{N}(i)|} \sum_{j \in \mathcal{N}(i)} \frac{\langle \mathbf{v}_i, \mathbf{v}_j \rangle}{\|\mathbf{v}_i\|_2 \|\mathbf{v}_j\|_2} = \frac{1}{|\mathcal{N}(i)|} \sum_{j \in \mathcal{N}(i)} \cos(\mathbf{v}_i, \mathbf{v}_j)$$

- **품질 판정 기준**:
  - (아래 0.75/0.50 기준은 scVelo가 제시한 값이 아니라 경험적 예시이며 데이터에 맞춰 조정합니다.)
  - $C_{\text{vel}}(i) \ge 0.75$: 높은 신뢰도. 주변 세포들과 통일된 동역학 흐름을 형성함.
  - $C_{\text{vel}}(i) < 0.50$: 잡음 세포. 전사 버스팅이나 기술적 드롭아웃으로 인해 무작위 방향으로 요동치고 있으므로 궤적 해석에서 제외 권고.

### 2.2 그래프 일관성 (Graph Coherence)
세포 $i$의 속도 벡터 $\mathbf{v}_i$와 전이 행렬 $\mathbf{P}$로부터 예측된 국소 기대 변위 벡터 $\mathbf{d}_i = \sum_{j \in \mathcal{N}(i)} P_{ij} (\mathbf{x}_j - \mathbf{x}_i)$ 간의 코사인 유사도를 평가합니다.

$$\text{Coh}(i) = \frac{\langle \mathbf{v}_i, \mathbf{d}_i \rangle}{\|\mathbf{v}_i\|_2 \|\mathbf{d}_i\|_2}$$

- 속도 벡터와 kNN 그래프 상의 전이 확률류가 조화를 이루는지 검증하는 척도입니다.

### 2.3 경계 교차 정확도 (Cross-Boundary Correctness: CBC)
알려진 생물학적 계통 전이 경계(Boundary, 예: 루미날 상피 $\to$ 가소성 중간체)를 가로지르는 전이 확률류의 순방향성(Net Forward Flux)을 검증합니다.

$$\text{CBC}(A \to B) = \frac{\sum_{i \in A} \sum_{j \in B} P_{ij} - \sum_{i \in B} \sum_{j \in A} P_{ji}}{\sum_{i \in A} \sum_{j \in B} P_{ij} + \sum_{i \in B} \sum_{j \in A} P_{ji}}$$

- $\text{CBC} \in [-1, 1]$: $+1$에 가까울수록 순방향 전이가 완벽히 복원된 것이며, 음수가 나타나면 역류 아티팩트가 발생한 것입니다.

---

## 3. 실무 의사결정 트리 및 배치 효과 보정 프로토콜

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   CellRank 2 실무 분석 의사결정 트리 (Decision Tree)                   │
│                                                                                        │
│   [입력 데이터 진단]                                                                   │
│   ├── 시계열 샘플(Longitudinal Timepoints)이 존재하는가?                               │
│   │   ├── YES ──► RealTimeKernel (WOT, w = 0.6) + VelocityKernel (w = 0.4)             │
│   │   └── NO  ──► VelocityKernel (w = 0.7) + CytoTRACEKernel (w = 0.3)                 │
│   │                                                                                    │
│   ├── 공간 전사체(Spatial Transcriptomics, Track 02)가 결합되어 있는가?                │
│   │   ├── YES ──► PrecomputedKernel로 공간 전이행렬 결합 (SpatialKernel 없음) │
│   │   └── NO  ──► 단일 세포 발현 공간 매니폴드 유지                                    │
│   │                                                                                    │
│   └── 샘플 간 배치 효과(Batch Effect)가 존재하는가?                                    │
│       ├── 주의: Spliced / Unspliced 원본 카운트 행렬에는 배치 보정 절대 금지!           │
│       └── 원칙: scVI / Harmony 잠재 공간(Latent space)에서만 kNN 그래프 구축           │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 배치 효과(Batch Effect) 처리 시의 절대 원칙
> **경고 (Do Not Normalize Spliced/Unspliced Across Batches)**:
> Harmony, ComBat, Scanorama 등의 배치 보정 도구를 spliced($s$) 및 unspliced($u$) 카운트 행렬 자체에 직접 적용해서는 안 됩니다.
> - 이유: 배치 보정 알고리즘은 비선형 스케일링을 가하여 실제 분자 카운트의 물리적 의미를 왜곡합니다. 이로 인해 스플라이싱 동역학의 근본 물리 항등식인 $\frac{ds}{dt} = \beta u - \gamma s$의 파라미터 적합이 완전히 붕괴됩니다.
> - **올바른 접근법**: 배치 보정은 세포 간 이웃 관계를 정의하는 저차원 임베딩 공간(scVI 잠재 공간 또는 Harmony PCA 공간)에만 적용하여 kNN 그래프 $\mathcal{G}$를 생성하고, 속도 계산은 배치 보정하지 않은 spliced/unspliced 레이어(라이브러리 정규화와 kNN 평활 모멘트 Ms, Mu를 거친 값)를 기반으로 수행해야 합니다.

---

## 4. 5단계 매니폴드 및 궤적 분석 캡스톤 프로젝트 루브릭 (Capstone Rubric)

전립선암 거세 저항성 계통 가소성 전이 분석을 완벽히 수행하기 위한 5단계 캡스톤 평가 루브릭입니다.

| 단계 (Stage) | 핵심 작업 (Task) | 정량적 합격 기준 (Success Metrics) | Track 연계 |
| :--- | :--- | :--- | :--- |
| **Stage 1: QC & Splicing Audit** | Velocyto/Alevin-fry 전사체 양적 감사 | - Unspliced 비율: 샘플 유형별 기준 설정(예시: scRNA-seq $15\% \sim 25\%$, snRNA-seq는 더 높아 $40\% \sim 60\%$)<br>- 유전자별 최소 전사체 수 $\ge 20$ | Track 00 |
| **Stage 2: Manifold & Multi-Kernel** | scVI/PCA 매니폴드 및 다중 커널 융합 | - VelocityConfidence 평균 $\ge 0.75$<br>- 볼록 결합 행-확률성 $\sum_j P_{ij} = 1.0$ | Track 02 |
| **Stage 3: G-PCCA Macrostates** | 실 슈어 분해 및 거시상태 분해 | - 스펙트럼 갭 $\Delta_k$ 극대화 지점 ($k=4 \sim 6$)<br>- 말단 거시상태 자체 전이 $\tilde{P}_{ss} \ge 0.90$ | Track 04 |
| **Stage 4: Fate Probabilities** | GMRES 기반 흡수 마르코프 체인 | - 흡수 행렬 행합 $\sum_k B_{ik} = 1.0$<br>- 잔차 수렴 허용오차 $\le 10^{-7}$ | Track 06 |
| **Stage 5: Driver Discovery** | 음이항 GAM 스플라인 및 eGRN 교차 | - FDR $q < 10^{-5}$, 상관계수 절댓값 $> 0.5$<br>- SCENIC+ regulon 활성 일치도 | Track 01 |

> 위 합격 기준(0.75, 0.90, $10^{-5}$ 등)은 출처가 없는 예시 임계이며 데이터에 맞춰 조정해야 합니다. 예컨대 이 저장소의 합성 데이터 파이프라인 01의 평균 속도 일관성은 0.161로 Stage 2 기준에 못 미칩니다.

---

## 5. 실전 진단 및 디버깅 파이썬 코드 스니펫

```python
import scvelo as scv
import cellrank as cr
import numpy as np

# 1. 속도 신뢰도 (Velocity Confidence) 계산 및 저품질 세포 필터링
scv.tl.velocity_graph(adata)        # velocity_confidence 내부(transition 단계)가 velocity_graph를 요구함
scv.tl.velocity_confidence(adata)
low_conf_cells = adata.obs["velocity_confidence"] < 0.50
print(f"Flagged Low Confidence Cells: {np.sum(low_conf_cells)} / {adata.n_obs}")

# 2. 결합 커널 구축 (배치 보정된 scVI 이웃 그래프 사용)
# scVI latent space 기반 kNN 그래프가 adata.obsp["connectivities"]에 저장되어 있어야 함
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

ck = cr.kernels.CytoTRACEKernel(adata)
ck.compute_cytotrace()
ck.compute_transition_matrix()

# 가중 볼록 결합
combined_kernel = 0.7 * vk + 0.3 * ck

# 3. 행-확률성 및 결손 행 검증
row_sums = np.array(combined_kernel.transition_matrix.sum(axis=1)).flatten()
assert np.allclose(row_sums, 1.0), "Error: Transition matrix is not row-stochastic!"

# 4. G-PCCA 거시상태 안정성 평가
g = cr.estimators.GPCCA(combined_kernel)
g.compute_schur(n_components=20)
g.compute_macrostates(n_states=6, cluster_key="cell_type")
g.predict_terminal_states()   # 흡수 확률 계산 전 말단 상태 지정(없으면 RuntimeError)

# 조대입자화 전이 행렬의 대각 성분 확인
diag_probs = np.diag(g.coarse_T)
print("Coarse-grained Self-transition Probabilities:", diag_probs)

# 5. 흡수 확률 계산 (GMRES 수치 해법)
g.compute_fate_probabilities(solver="gmres", tol=1e-7)
assert np.allclose(g.fate_probabilities.X.sum(axis=1), 1.0), "Error: Fate probabilities do not sum to 1!"
print("Trajectory Pipeline Validation Passed Successfully!")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [궤적 역전 진단 및 디버깅 (Root-Sink Inversion Debugging)]
연구원이 Enzalutamide 치료 저항성 mCRPC 환자 코호트의 snRNA-seq 데이터를 scVelo 동역학 모델로 분석하였다.
결과 시각화 상에서 세포 분화 벡터가 t-NEPC 클러스터에서 출발하여 초기 미치료 루미날 선암종 클러스터로 거꾸로 흘러 들어가는 치명적인 '루트-싱크 역전(Root-Sink Inversion)' 현상이 발생하였다.
1) 이러한 현상을 유발하는 3대 근본 원인을 [1] 약제 약리 기전(전사 감속), [2] 유전자별 상이한 스플라이싱 동역학, [3] 미성숙 RNA 검출 품질 관점에서 진단하시오.
2) `CellRank 2` 프레임워크를 적용하여 이 문제를 해결하고 올바른 순방향 궤적을 복원하기 위한 구체적인 3단계 디버깅 프로토콜을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 3대 근본 원인 진단**:
  - **[1] 약제 약리 기전에 의한 전사 급감 (Transcriptional Deceleration Artifact)**:
    - Enzalutamide는 AR 리간드 결합을 저해하여 전립선 루미날 유전자의 신규 전사를 즉각 중단시킵니다.
    - 미성숙 RNA($u$)가 성숙 RNA($s$)보다 훨씬 빠르게 고갈되어 $u < \frac{\gamma}{\beta}s$가 되고 해당 유전자의 $v < 0$이 됩니다. 이 자체는 루미날 이탈(정방향)이지만, 속도 투영이 하향 조절 유전자에 지배되거나 파라미터 추정이 틀리면 세포가 과거 루미날 상태로 되돌아가는 것처럼 보이는 거짓 궤적이 형성될 수 있습니다.
  - **[2] 유전자별 상이한 스플라이싱 동역학 (Multiple Kinetics)**:
    - 선암종과 신경내분비암종 간에 스플라이싱 반응속도($\beta$) 및 분해 속도($\gamma$)가 상이함에도 불구하고, scVelo가 단일 타원 곡선 적합을 강제하여 기울기(Steady-state ratio $\gamma/\beta$)를 잘못 추정함으로써 속도 부호가 반전되었습니다.
  - **[3] 단일핵 시퀀싱(snRNA-seq)의 미성숙 RNA 검출 편향**:
    - 핵 내에는 세포질에 비해 성숙 RNA가 적고 인트론을 함유한 미성숙 RNA 비율이 높습니다(예시 수준으로 $40\sim 60\%$; 조직·정량화 방식에 따라 다름). 만약 정량화 단계(Alevin-fry/Velocyto)에서 인트론 매핑이 왜곡되거나 드롭아웃이 심할 경우, 정상적인 스플라이싱 비율 관계가 붕괴됩니다.
- **2. CellRank 2 기반 3단계 디버깅 해결 프로토콜**:
  - **1단계 (역류 유전자 필터링 및 동역학 정제)**:
    - 루미날 표지자(*KLK3, NKX3-1*) 등 전사 급감이 발생한 유전자들의 상평형 곡선을 검사하여 위상 왜곡을 유발하는 상위 유전자들을 속도 추정 피처(`adata.var['velocity_genes']`)에서 제외합니다.
  - **2단계 (CytoTRACE / Pseudotime 커널과의 볼록 결합)**:
    - 줄기세포능과 발달 잠재력을 대변하는 `CytoTRACEKernel`을 계산하여, 미분화 가소성 세포가 높은 잠재력을 갖고 분화된 선암종이 낮은 잠재력을 갖도록 시간의 화살(Arrow of Time)을 정립합니다.
    - 결합 가중치를 $w_{\text{vel}} = 0.5, w_{\text{cyto}} = 0.5$로 설정하여 속도 벡터의 일시적 음수 왜곡을 상쇄합니다.
  - **3단계 (G-PCCA 거시상태 분석을 통한 시작점/말단점 앵커링)**:
    - 슈어 분해 후 조대입자화 전이 행렬 $\tilde{\mathbf{P}}$를 검토하여, 루미날 선암종을 초기 상태(`initial_states`)로, t-NEPC 및 DNPC를 말단 상태(`terminal_states`)로 생물학적 마커 기반 수동 앵커링(Biological Anchoring)을 수행함으로써 흡수 마르코프 체인의 흐름을 강건하게 고정합니다.

</details>

---

### Q2. [정량적 품질 지표 수리 계산 (Diagnostic Metrics Calculation)]
단일 세포 매니폴드 상의 특정 세포 $i$의 속도 벡터가 $\mathbf{v}_i = (2.0, 1.0)^T$이다.
세포 $i$의 3개 최근접 이웃 세포 $\{j_1, j_2, j_3\} \in \mathcal{N}(i)$의 속도 벡터가 다음과 같다:
- $\mathbf{v}_{j_1} = (1.5, 0.75)^T$
- $\mathbf{v}_{j_2} = (0.0, 2.0)^T$
- $\mathbf{v}_{j_3} = (-2.0, -1.0)^T$

1) 세포 $i$와 각 이웃 세포 간의 코사인 유사도 $\cos(\mathbf{v}_i, \mathbf{v}_j)$를 각각 계산하고, 속도 신뢰도 $C_{\text{vel}}(i)$의 최종 수치를 소수점 둘째 자리까지 도출하시오.
2) 계산된 $C_{\text{vel}}(i)$ 수치를 품질 합격 기준(임계값 $\ge 0.75$)과 비교하여, 이 세포를 신뢰할 수 있는 전이 세포로 채택할 수 있는지 판정하고 그 이유를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 코사인 유사도 및 속도 신뢰도 $C_{\text{vel}}(i)$ 계산**:
  - 세포 $i$의 크기(L2 norm):
    $$\|\mathbf{v}_i\|_2 = \sqrt{2.0^2 + 1.0^2} = \sqrt{4 + 1} = \sqrt{5} \approx 2.236$$
  - **이웃 $j_1$과의 코사인 유사도**:
    $$\mathbf{v}_{j_1} = (1.5, 0.75)^T = 0.75 \times (2, 1)^T$$
    방향이 완전히 일치하므로:
    $$\cos(\mathbf{v}_i, \mathbf{v}_{j_1}) = \frac{2(1.5) + 1(0.75)}{\sqrt{5} \times \sqrt{1.5^2 + 0.75^2}} = \frac{3.75}{\sqrt{5} \times \sqrt{2.8125}} = \frac{3.75}{\sqrt{5} \times 1.677} = \frac{3.75}{3.75} = +1.00$$
  - **이웃 $j_2$와의 코사인 유사도**:
    $$\|\mathbf{v}_{j_2}\|_2 = \sqrt{0^2 + 2.0^2} = 2.0$$
    내적: $\langle \mathbf{v}_i, \mathbf{v}_{j_2} \rangle = 2(0) + 1(2.0) = 2.0$
    $$\cos(\mathbf{v}_i, \mathbf{v}_{j_2}) = \frac{2.0}{\sqrt{5} \times 2.0} = \frac{1}{\sqrt{5}} \approx \frac{1}{2.236} \approx +0.45$$
  - **이웃 $j_3$와의 코사인 유사도**:
    $$\mathbf{v}_{j_3} = (-2.0, -1.0)^T = -1 \times \mathbf{v}_i$$
    방향이 정반대이므로:
    $$\cos(\mathbf{v}_i, \mathbf{v}_{j_3}) = -1.00$$
  - **속도 신뢰도 $C_{\text{vel}}(i)$ 평균**:
    $$C_{\text{vel}}(i) = \frac{1}{3} \left( \cos(\mathbf{v}_i, \mathbf{v}_{j_1}) + \cos(\mathbf{v}_i, \mathbf{v}_{j_2}) + \cos(\mathbf{v}_i, \mathbf{v}_{j_3}) \right)$$
    $$C_{\text{vel}}(i) = \frac{1.00 + 0.45 - 1.00}{3} = \frac{0.45}{3} = 0.15$$
- **2. 품질 판정 및 결론**:
  - **판정**: **불합격 (Reject / Flag as Noise)**
  - **이유**:
    - 계산된 속도 신뢰도는 $C_{\text{vel}}(i) = 0.15$로, 품질 임계값 $0.75$에 훨씬 미치지 못할 뿐만 아니라 $0.50$ 이하의 극심한 잡음 상태입니다.
    - 특히 이웃 $j_3$와의 방향이 완전히 반대($-1.0$)인 것은 국소 이웃 내에서 속도 벡터의 방향이 심각하게 상충하고 있음을 나타냅니다.
    - 이는 국소 전사 요동이나 시퀀싱 드롭아웃에 기인한 노이즈이므로, 이 세포는 신뢰할 수 있는 계통 궤적의 구성원으로 사용할 수 없으며 후속 분석에서 마스킹(Filtering)해야 합니다.

</details>

---

### Q3. [캡스톤 실전 임상 파이프라인 종합 설계 (Capstone Clinical Problem Solving)]
mCRPC 환자 30명으로부터 [1] Enzalutamide 치료 개시 전(Baseline) 및 [2] 질병 진행(Progression) 시점의 쌍 생검(Paired biopsies) 조직을 채취하여 단일핵 멀티옴(snRNA-seq + snATAC-seq, Track 01), 단일 세포 복제수 변이(scCNA, Track 04), 10x Visium 공간 전사체(Spatial Transcriptomics, Track 02) 데이터를 확보하였다.
1) 이 환자군에서 AR 저항성 계통 가소성(Lineage Plasticity) 전이를 완벽히 복원하기 위해, 본 캡스톤 루브릭의 5단계(Stage 1 ~ Stage 5)에 입각한 엔드-투-엔드 분석 파이프라인의 핵심 수행 절차를 기술하시오.
2) 질병 진행 생검에서 t-NEPC로 전환된 환자군(Cohort A)과 DNPC로 전환된 환자군(Cohort B)의 치료 전(Baseline) 상태를 역추적했을 때, 사전 존재하는 유전체 클론(Pre-existing Subclone)과 치료 유도 적응성(De Novo Plasticity)을 판별하기 위한 멀티오믹스 판정 기준을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 5단계 캡스톤 엔드-투-엔드 분석 파이프라인 수행 절차**:
  - **Stage 1 (QC 및 정량화 감사)**:
    - Alevin-fry를 통해 인트론 및 엑손 리드를 정량화하고, 각 생검 샘플별 미성숙 RNA 비율이 샘플 유형별 예시 범위(snRNA-seq는 $40\sim 60\%$ 수준, scRNA-seq는 $15\sim 25\%$ 수준)에 드는지 전수 감사. 저품질 세포 제거.
  - **Stage 2 (매니폴드 구축 및 다중 커널 융합)**:
    - 배치 효과 제거를 위해 MultiVI(`scvi.model.MULTIVI`, RNA+ATAC 결합 잠재 공간; scVI는 RNA 전용) 또는 WNN/MOFA+를 활용하여 단일핵 전사체와 snATAC-seq의 결합 잠재 공간을 구축.
    - 원시 카운트로부터 `VelocityKernel`을 계산하고, WOT 기반 `RealTimeKernel`(치료 전 vs 진행 시점)과 `CytoTRACEKernel`, 그리고 Track 02 Visium 공간 그래프에서 만든 전이행렬을 `PrecomputedKernel`로 구성한 커널(cellrank.kernels에 `SpatialKernel`은 없음)을 볼록 결합($\mathbf{P}_{\text{comb}} = 0.4 \mathbf{P}_{\text{ot}} + 0.3 \mathbf{P}_{\text{vel}} + 0.2 \mathbf{P}_{\text{cyto}} + 0.1 \mathbf{P}_{\text{sp}}$)하여 단일 행-확률 전이 행렬 도출.
  - **Stage 3 (G-PCCA 마르코프 체인 및 거시상태 분해)**:
    - 실 슈어 분해(Real Schur Decomposition)를 수행하고 스펙트럼 갭($\Delta_k$)을 평가하여 6대 거시상태(초기 선암종, 가소성 중간체, AR+ 증폭 선암종, Amphicrine, DNPC, t-NEPC)를 도출하고, 조대입자화 전이 행렬 $\tilde{\mathbf{P}}$ 대각 성분 검증($\tilde{P}_{ss} \ge 0.90$).
  - **Stage 4 (흡수 운명 확률 및 분기 엔트로피 계산)**:
    - GMRES 크릴로프 반복법을 적용하여 $(\mathbf{I} - \mathbf{Q})\mathbf{B} = \mathbf{R}$ 선형계를 풀어 t-NEPC 및 DNPC 흡수 확률 행렬 $\mathbf{B}$ 계산.
    - 가소성 체크포인트 세포의 계통 프라이밍 엔트로피($H_i$) 피크 영역 동정.
  - **Stage 5 (드라이버 발굴 및 멀티옴 교차 검증)**:
    - 음이항 GAM 스플라인 모델을 통해 $B_{\text{t-NEPC}}$ 및 $B_{\text{DNPC}}$를 구동하는 유전자 캐스케이드(EZH2/SOX2 $\to$ ASCL1/INSM1 vs FGFR1 $\to$ ETV4)를 발굴하고, Track 01 SCENIC+ regulon 활성도 및 Track 04 CNA 서브클론과 일치도 검증.
- **2. Pre-existing vs De Novo 가소성 판별을 위한 멀티오믹스 기준**:
  - **Cohort A (t-NEPC 전환군)의 판별 기준**:
    - **Pre-existing 모델 증거**:
      - 치료 전(Baseline) 생검의 Track 04 scCNA 데이터에서 극미량($0.1\sim 1\%$)이라도 Chr 17p (*TP53*) 및 Chr 13q (*RB1*) 영역의 복제수 결손을 지닌 서브클론이 이미 검출되는 경우. (scCNA는 복제수만 측정하므로 양대립 결손의 확정에는 점돌연변이 호출·LOH 분석이 별도로 필요하고, 0.1% 클론은 총 세포 수와 시퀀싱 깊이에 따라 검출되지 않을 수 있음)
      - 이 서브클론 세포들의 t-NEPC 흡수 확률이 치료 전부터 $B_{\text{t-NEPC}} > 0.85$로 기배정되어 있는 경우.
    - **De Novo 적응 모델 증거**:
      - 치료 전에는 종양 세포 전체가 유전체적으로 *TP53/RB1* 정상 야생형(Wild-type)이며, 치료 후 생검에서 비로소 새로운 복제수 결손이나 점돌연변이가 출현한 경우.
      - 치료 전에는 $B_{\text{t-NEPC}} \approx 0$이었으나, Enzalutamide 노출 후 크로마틴 리모델러(EZH2) 활성화와 함께 가소성 중간체($H_i > 1.8$)를 거쳐 단계적으로 흡수 확률이 증가한 경우.
  - **Cohort B (DNPC 전환군)의 판별 기준**:
    - *RB1* 결손 없이 *FGFR1* 증폭 등(가설적 예시)을 지닌 서브클론의 유무와, Track 01 eGRN 상에서 *ETV4* regulon의 전사적 활성화 시점을 추적하여 진화적 기원을 판정함.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- CellRank 2.3.3 (github.com/scverse/cellrank). 원 논문: Weiler P, Lange M, Klein M, et al. CellRank 2: unified fate mapping in multiview single-cell data. *Nat Methods*. 2024;21(7):1196-1205. PMID: 38871986. doi: 10.1038/s41592-024-02303-9
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- scVelo 0.3.4 (github.com/theislab/scvelo). 원 논문: Bergen V, Lange M, Peidli S, et al. Generalizing RNA velocity to transient cell states through dynamical modeling. *Nat Biotechnol*. 2020;38(12):1408-1414. PMID: 32747759. doi: 10.1038/s41587-020-0591-3
<!-- REFERENCES:END -->

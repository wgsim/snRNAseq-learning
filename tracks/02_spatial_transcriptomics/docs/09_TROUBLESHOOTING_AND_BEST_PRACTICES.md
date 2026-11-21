# Chapter 09. 트러블슈팅, 품질 제어 및 캡스톤 프로젝트 루브릭 (Troubleshooting, QC & Capstone Rubric)

> **Cognitive Goal**: 공간전사체 실험 및 데이터 분석 파이프라인 전주기에서 마주치는 4대 핵심 기술적 난제(프로브 비특이적 결합 노이즈, 조직 박리 및 주름 아티팩트, 효소 투과화에 의한 mRNA 측면 확산 블러, 디콘볼루션 다중공선성)의 근본 원인을 진단하고 통계적·실험적 해결책을 정립한다. 이어 대학원 수준의 독립 연구 역량을 검증하기 위한 '5단계 공간전사체 캡스톤 프로젝트 루브릭(5-Stage Capstone Project Rubric)'을 제시하여, 전처리부터 임상적 신호 플럭스 규명에 이르는 엄밀한 평가 기준을 확립한다.

---

## 1. 공간전사체 현장 트러블슈팅 가이드 (Field Troubleshooting)

공간전사체 기술은 기존의 용액 기반 시퀀싱과 달리, **조직 절편의 물리적 특성, 광학 현미경 이미징, 효소 반응의 동역학, 초다변량 공간 통계**가 복합적으로 얽혀 있어 독특한 기술적 노이즈와 아티팩트가 발생합니다.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│                   공간전사체 4대 주요 아티팩트 및 발생 기전                      │
│                                                                                  │
│   [1. 프로브 배경 잡음 (Background Noise)]                                       │
│   - 원인: FFPE 자가형광, 비특이적 프로브 결합, 미결합 프로브의 불완전 세척      │
│   - 증상: 조직이 없는 빈 슬라이드 영역에서 인공 카운트 검출                      │
│                                                                                  │
│   [2. 조직 박리 및 주름 (Detachment & Folding)]                                  │
│   - 원인: 동결 절편 부착 불량, 탈수 수축, 기포 형성                             │
│   - 증상: 기준 마커(Fiducial Frame) 왜곡, 2겹 중첩으로 인한 카운트 2배 폭증     │
│                                                                                  │
│   [3. mRNA 측면 확산 블러 (Lateral Diffusion Blur)]                             │
│   - 원인: 단백질분해효소(Pepsin) 과다 투과화 ──► 세포외 기질로 mRNA 유출         │
│   - 증상: 종양 특이 마커(KLK3)가 인접 기질/지방 조직 스팟으로 스며드는 흐림 현상 │
│                                                                                  │
│   [4. 디콘볼루션 다중공선성 (Deconvolution Collinearity)]                       │
│   - 원인: 유사 암세포 클론(AR+ vs Amphicrine) 간 시그니처 상관계수 r > 0.95      │
│   - 증상: 인접 스팟 간 가중치가 0과 최대치를 교대로 진동하는 불안정성            │
└──────────────────────────────────────────────────────────────────────────────────┘
```

---

### 1.1 이슈 1: 프로브 비특이적 배경 잡음 (Probe Background Noise & Non-specific Binding)

#### [원인 및 병태 기전]
- 10x Visium FFPE와 Xenium은 RNA에 하이브리드화하는 한 쌍의 올리고 프로브(LHS/RHS) 결합 및 결찰(Ligation)에 의존합니다(CosMx는 결찰 방식이 아니라 ISH 프로브와 광절단 리포터를 사용하는 별개 방식).
- 고정(Fixation) 불량 검체나 잔여 파라핀, 조직 내 리포푸신(Lipofuscin) 자가형광, 프로브 세척(Washing) 단계의 온도/염 농도 불일치로 인해 프로브가 비특이적으로 유리 슬라이드 표면에 흡착되는 현상이 발생합니다.

#### [진단 지표 (Diagnostic Metrics)]
1. **음성 대조 프로브(Negative Control Probes) 비율 검사**:
   - 인간 게놈에 존재하지 않는 염기서열(예: 세균성 프로브)로 설계된 음성 대조 프로브의 UMI 카운트가 전체 스팟에서 검출되는 비율을 계산합니다:
     $$\text{Noise Ratio} = \frac{\sum_{g \in \text{NegCtrl}} Y_{s, g}}{\sum_{g \in \text{All}} Y_{s, g}} \times 100$$
   - 예시 기준(출처 없는 경험적 임계값, 플랫폼·패널별로 설정): 정상 $< 1.0\%$, $3\%$ 초과 시 비특이적 결합이 심각.
2. **배경 스팟 발현량(Off-tissue Counts)**:
   - H&E 또는 DAPI 이미지 상 조직이 전혀 존재하지 않는 슬라이드 공백 영역의 스팟 UMI 중앙값이 50 카운트를 초과하는 경우(예시 기준; 데이터별 설정).

#### [해결 및 보정 프로토콜]
- **배경 차감(Background Subtraction)**:
  음성 대조 프로브의 스팟별 평균 카운트 $\mu_{\text{neg}}(s)$를 추정하여 원시 카운트에서 제로 절단 차감(Zero-truncated Subtraction)을 수행합니다:
  $$Y_{s, g}^{\text{corrected}} = \max\left(0, \; Y_{s, g} - \mu_{\text{neg}}(s)\right)$$
- **프로브 필터링**: 전체 조직에서 음성 대조 프로브의 최대 카운트보다 낮은 발현을 보이는 하위 $15\%$ 유전자를 SVG 및 디콘볼루션 분석 풀에서 사전 배제합니다.

---

### 1.2 이슈 2: 조직 박리 및 주름 아티팩트 (Tissue Detachment & Folding)

#### [원인 및 병태 기전]
- 극저온 박절(Cryosectioning) 또는 FFPE 수조 부유(Water bath floating) 과정에서 조직이 슬라이드 포획 영역에 불완전하게 밀착되면, 염색 및 투과화 시약 처리 과정에서 가장자리가 말려 올라가거나(Detachment) 조직 절편이 겹쳐지는(Folding/Wrinkling) 현상이 발생합니다.

#### [진단 지표 (Diagnostic Metrics)]
1. **국소 UMI 폭증 이상치(Local Library Size Outlier)**:
   - 주름이 생겨 조직이 2겹으로 겹친 구역은 단위 면적당 세포 수가 2배가 되므로, 인접한 단일 층 스팟들에 비해 총 UMI 수($N_s$)와 검출 유전자 수($G_s$)가 갑작스럽게 $2.0 \sim 2.5$배 급증하는 공간적 이상치(High-High Outlier)를 형성합니다.
2. **피듀셜 프레임 오정렬(Fiducial Misalignment)**:
   - 슬라이드 모서리의 육각 정렬 마커(Fiducial Frame)가 손상되거나 조직이 프레임을 침범하여 공간 이미지와 스팟 바코드 격자 간의 회전 각도 오차가 발생하는 경우.

#### [해결 및 보정 프로토콜]
- **주름 영역 마스킹(Spatial Masking)**:
  현미경 광학 이미지 상에서 주름 및 박리가 확인된 스팟들은 단순 정규화로 보정할 수 없습니다. 수동 어노테이션 도구(`Loupe Browser`) 또는 컴퓨터 비전 기반 임계 분할을 통해 해당 스팟 바코드를 `adata.obs['is_artifact'] = True`로 라벨링한 후 다운스트림 분석에서 완전히 필터링(Removal)해야 합니다.
- **아핀 변환(Affine Transformation) 재정렬**:
  피듀셜 정렬 실패 시, OpenCV 또는 SimpleITK를 사용하여 스팟 중심 좌표와 현미경 이미지 사이의 3점 호모그래피(Homography) 행렬을 수동 재계산합니다.

---

### 1.3 이슈 3: mRNA 측면 확산 블러 (Lateral Diffusion Blur)

#### [원인 및 병태 기전]
- 10x Visium의 가장 중대한 실험적 취약점입니다. 조직 절편 표면에 결합된 mRNA를 포획 올리고로 유도하기 위해 단백질분해효소(Pepsin/Permeabilization enzyme)를 처리합니다.
- 투과화 시간이 최적 시간(Tissue-specific optimal time)보다 과다하게 길어지면, 세포막이 완전히 붕괴된 상태에서 방출된 mRNA 분자들이 수직으로 가라앉지 않고 조직액 내부에서 측면 방향으로 브라운 운동을 하며 확산(Lateral Diffusion)됩니다.

```
[최적 투과화 (Optimal)]                   [과투과화에 의한 측면 확산 (Over-permeabilization)]
    세포 A         세포 B                      세포 A          세포 B
      │              │                           ＼          ／
      ▼ 수직 이동    ▼ 수직 이동                   ＼ 측면 확산／
 ─────────────  ─────────────               ─────────────  ─────────────
   스팟 A 포획     스팟 B 포획                  스팟 A 포획     스팟 B 포획
 (순수 시그니처 보존)                       (A와 B의 전사체가 뒤섞여 경계가 흐려짐!)
```

#### [진단 지표 (Diagnostic Metrics)]
1. **조직 경계부 유출 지수(Boundary Leakage Index)**:
   - 종양 세포가 전혀 존재하지 않는 종양 외곽 정상 기질 또는 지방 조직 스팟에서 종양 특이 고발현 유전자(예: *KLK3, AMACR*)의 분자 카운트가 검출되는 거리를 측정합니다.
   - 확산 거리가 인접 스팟 중심 간격(Visium $100\ \mu\text{m}$)에 가까워지면(예: $\gtrsim 50\ \mu\text{m}$, 예시 기준) 과투과화를 의심합니다.
2. **미토콘드리아 전사체 누출 비율(Mitochondrial Diffusion Ratio)**:
   - 과투과화·세포질 mRNA 유실 시 `%mt` 비율이 변할 수 있다는 가설이 있으나 기전은 확립되지 않았습니다(임계값 $30\%$는 예시). FFPE Visium v1 프로브 세트에는 MT 프로브가 없어 `%mt` 지표 자체가 적용되지 않을 수 있습니다(3장 참조).

#### [해결 및 보정 프로토콜]
- **공간 디블러링(Richardson-Lucy 역합성곱; 라플라시안 샤프닝과는 다른 방법)**:
  측면 확산을 2차원 가우시안 점 확산 함수(Point Spread Function, PSF: $h(x, y) = \frac{1}{2\pi \sigma_{\text{diff}}^2} e^{-\frac{x^2+y^2}{2\sigma_{\text{diff}}^2}}$)로 모델링하고, 리처드슨-루시(Richardson-Lucy) 역합성곱 알고리즘을 적용하여 인접 스팟으로 번진 카운트를 원래의 발현 중심으로 복원하려는 시도입니다. 단, Richardson-Lucy는 포아송 노이즈 최대우도 역합성곱이며 Visium 스팟 피치(100 µm)가 확산 PSF 폭($\sigma_{\text{diff}}$ 약 20~50 µm 가정)보다 커 격자 샘플링이 부족하면 복원이 불량조건이 되므로 적용 한계를 확인해야 하고, 대안으로 SpotClean 등 공간 오염 보정 도구를 고려합니다.
- **실험 조건 수정**: 동일 조직 타입의 예비 슬라이드에서 형광 cDNA 시간 적정(Tissue Optimization; 시간 예시 6·12·18·24분, 실제 범위는 제조사 프로토콜 확인)을 재수행하여 확산 없이 신호 대 잡음비가 최대화되는 지점을 재설정합니다.

---

### 1.4 이슈 4: 디콘볼루션 다중공선성 및 식별 불가능성 (Deconvolution Collinearity)

#### [원인 및 병태 기전]
- snRNA-seq 레퍼런스 데이터에 지나치게 미세하고 중복된 세포 아형(예: 'AR+ Luminal Subtype 1'과 'AR+ Luminal Subtype 2')을 독립된 기준 시그니처로 주입할 때 발생합니다.
- 두 세포 상태 간의 전사체 상관계수가 $r > 0.98$이면 시그니처 행렬이 거의 선형 종속(조건수 증가; 랭크 결손은 $r=1$에서만 발생)에 가까워져, 베이지안 MCMC 샘플링 시 사후분포가 다봉성(Multimodal ridge)을 형성하고 두 세포 유형의 존재량이 스팟마다 무작위로 뒤바뀌는 현상이 발생합니다.

#### [진단 지표 (Diagnostic Metrics)]
1. **분산 팽창 계수(Variance Inflation Factor, VIF)**:
   - 시그니처 행렬 $G$의 각 세포 유형 $k$에 대해 다른 모든 세포 유형들로 설명되는 결정계수 $R_k^2$를 구하여 VIF를 산출합니다:
     $$\text{VIF}_k = \frac{1}{1 - R_k^2}$$
   - 판정 기준(일반적 경험칙): $\text{VIF}_k > 10$이면 심각한 다중공선성, $> 50$이면 병합을 검토합니다. 두 변수만 있을 때 $\text{VIF}=1/(1-r^2)$이므로 $r=0.90 \to 5.3$, $r=0.95 \to 10.3$, $r=0.98 \to 25$, $r \approx 0.99 \to 50$입니다.
2. **사후 추정치의 음의 상관관계(Posterior Anti-correlation)**:
   - cell2location 학습 후, 스팟 전체에 걸쳐 두 세포 유형의 사후 추정치 간 상관계수가 강한 음수($r < -0.7$)를 나타내는 경우 (하나가 올라가면 다른 하나가 기계적으로 0으로 떨어지는 식별 불가 증상).

#### [해결 및 보정 프로토콜]
- **계통 트리 기반 시그니처 병합(Hierarchical Signature Pruning)**:
  계층적 군집화를 통해 상관계수가 매우 높은(예시 기준 $r \gtrsim 0.9$; 데이터별로 설정하며 VIF 10 이상에 대응하는 $r \approx 0.95$ 이상은 거의 확실한 병합 대상) 미세 아형들을 하나의 상위 계통 단위(예: 'Luminal Progenitor'로 통합)로 클러스터를 조기 병합한 후 디콘볼루션을 재수행합니다.
- **사전분포 하이퍼파라미터 조정**:
  cell2location은 감마 사전분포 기반 변분추론 모델이며 $L_2$ 능형 정규화를 쓰지 않습니다. 공선성 완화의 정석은 참조 세포유형 병합·마커 선택이고, 보조적으로 `N_cells_per_location`, `detection_alpha` 같은 사전분포 하이퍼파라미터를 조정할 수 있습니다.

---

## 2. 종합 요약: 4대 아티팩트 진단 및 조치 매트릭스

| 아티팩트 유형 | 핵심 탐지 지표 | 허용 임계값 (Threshold) | 즉각적 조치 (Action Plan) |
| :--- | :--- | :--- | :--- |
| **1. 비특이적 프로브 노이즈** | 음성 대조 프로브 비율, 빈 공간 UMI | 전체 UMI의 $< 1.0\%$ | 음성 프로브 기반 배경 차감 및 하위 노이즈 유전자 배제 |
| **2. 조직 박리 및 주름** | 스팟 면적당 UMI 2배 폭증, H&E 주름 | 단일 스팟 z-score $> 3.5$ | 문제 스팟 마스킹 및 물리적 다운스트림 분석 배제 |
| **3. mRNA 측면 확산 블러** | 종양 마커 기질 유출 거리, %mt | 유출 거리 $< 50\ \mu\text{m}$ (예시 기준; Visium 스팟 간격은 100 µm) | 투과화 시간 최적화(Tissue Optimization)와 SpotClean 등 공간 오염 보정(Richardson-Lucy 역합성곱은 스팟 피치가 PSF보다 큰 경우 불량조건이므로 한계 확인) |
| **4. 디콘볼루션 공선성** | 시그니처 VIF, 사후 추정치 음의 상관 | VIF $< 10$, 사후 상관 $r > -0.5$ | 상관 계통 클러스터 병합 및 상위 계통 단위 재분석 |

---

## 3. 공간전사체 5단계 캡스톤 프로젝트 평가 루브릭 (Capstone Project Rubric)

본 캡스톤 루브릭은 임상 전립선암(mCRPC) 또는 동등한 복합 고형암 조직의 공간전사체 데이터를 독립적으로 분석하고, 병리학적 기전을 규명할 수 있는가를 평가하는 대학원 수준의 종합 역량 평가 기준입니다.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│                   5-Stage Spatial Transcriptomics Capstone Workflow              │
│                                                                                  │
│   [Stage 1: 전처리 & 공간 QC]       ──► UMI/유전자 분포, 음성 프로브, 피듀셜 정렬 │
│   [Stage 2: 세포 분할 & 할당]       ──► Baysor / Mesmer 단일 세포 경계 및 전사체   │
│   [Stage 3: 스팟 디콘볼루션]        ──► cell2location 절대 세포수 복원 & 공선성 통제│
│   [Stage 4: SVG & 공간 니치]        ──► Moran's I, 국소 윈도우 4대 니치 클러스터링 │
│   [Stage 5: 거리 제약 신호 & 번역]  ──► COMMOT 플럭스 벡터장, 면역 배제 기전 규명   │
└──────────────────────────────────────────────────────────────────────────────────┘
```

---

### Stage 1: 데이터 전처리 및 공간 품질 제어 (Preprocessing & Spatial QC)

- **평가 초점**: 조직 절편의 공간적 결함, 오염, 기술적 저품질 스팟을 통계학적으로 정확히 식별하고 필터링할 수 있는가?
- **평가 척도**:
  - **Distinction (90~100점)**:
    - 총 UMI 수, 검출 유전자 수, 미토콘드리아 전사체 비율의 분포를 조직 해부학적 구역(H&E 매핑)과 대조하여 다차원 이상치를 검출함.
    - 음성 대조 프로브를 기반으로 비특이적 신호를 정량화하여 유효 임계값을 초과하는 노이즈를 완벽히 필터링함.
    - 피듀셜 정렬 이상 및 조직 주름 구역을 컴퓨터 비전 기법 또는 엄밀한 이상치 마스킹으로 제거하고, 필터링 전후의 공간 분포 변화를 투명하게 입증함.
  - **Merit (80~89점)**:
    - 표준 임계값(예: UMI > 500, %mt < 20%)을 적용하여 기본적인 이상치 스팟을 필터링함.
    - 조직 이미지 상의 결함 부위를 수동으로 확인하고 배제하였으나, 음성 프로브 배경 차감 수식이 다소 형식적임.
  - **Pass (60~79점)**:
    - 기본적인 scanpy QC 함수만을 기계적으로 실행하였으며, 공간적 가장자리 결함이나 조직 주름의 영향을 깊이 있게 검토하지 못함.
  - **Fail (0~59점)**:
    - 조직이 없는 빈 슬라이드 스팟이나 심각한 주름 아티팩트를 제거하지 않고 그대로 다운스트림 분석에 포함하여 왜곡을 초래함.

---

### Stage 2: 단일 세포 분할 및 전사체 할당 (Cell Segmentation & Transcript Assignment)

- **평가 초점**: 서브마이크론 해상도(Xenium/CosMx) 또는 Visium HD 데이터에서 DAPI 핵 이미지를 기반으로 정확한 세포 경계를 분할하고 전사체 포획 효율을 극대화하였는가?
- **평가 척도**:
  - **Distinction (90~100점)**:
    - 단순 유클리드 핵 팽창(Nuclear Expansion)의 한계를 인식하고, `Mesmer`(딥러닝) 또는 `Baysor`(베이지안 전사체 사전확률 결합 모델)을 적용하여 정밀 세포 분할을 수행함.
    - 다핵 세포(파골세포), 무핵 적혈구, 조밀한 암세포 덩어리에서의 분할 분리 정확도(Segmentation F1-score)를 검증함.
    - 세포질 외부 유리 전사체(Unassigned Transcripts) 비율을 $15\%$ 미만으로 제어하면서 과분할(Over-segmentation) 및 병합 오류(Under-segmentation)를 완벽히 통제함.
  - **Merit (80~89점)**:
    - DAPI 기반 셀포즈(Cellpose) 또는 스타디스트(StarDist)를 적절히 활용하여 개별 핵을 검출하고 합리적인 반경으로 세포 경계를 확장함.
    - 분할 결과에 대한 정량적 지표가 다소 부족하나 대표 시야(ROI)에서 시각적 타당성을 확보함.
  - **Pass (60~79점)**:
    - 고정 반경(예: 핵 중심에서 $5\ \mu\text{m}$) 단순 원형 팽창을 적용하여 밀집 구역에서 세포 간 전사체 교차 오염이 관찰됨.
  - **Fail (0~59점)**:
    - 세포 분할에 실패하여 다수의 세포가 하나의 거대 덩어리로 뭉개지거나 단일 세포 분석 단계로 전환하지 못함.

---

### Stage 3: 스팟 디콘볼루션 및 세포 구성 복원 (Spot Deconvolution & Abundance)

- **평가 초점**: 단일 세포 기준 데이터셋과 공간전사체 카운트 모델 간의 수학적 정합성을 확보하고 절대 세포 수를 정확히 추정하였는가?
- **평가 척도**:
  - **Distinction (90~100점)**:
    - `cell2location` 또는 `RCTD`의 생성적 확률 모델(포아송/음이항)을 정확히 이해하고, H&E 병리 기반의 스팟당 평균 세포 수 사전정보($N_{\text{cells\_per\_location}}$)를 적절히 주입함.
    - snRNA-seq 레퍼런스의 미세 클론 간 다중공선성(VIF)을 사전 점검하고 고상관 시그니처를 계통학적으로 병합함.
    - 단순 비율(Proportion)의 착시를 배제하고 사후 5% 분위수(q05) 절대 세포 수를 도출하여 실제 조직학적 밀도와 완벽히 부합시킴.
  - **Merit (80~89점)**:
    - cell2location 파이프라인을 성공적으로 구동하고 ELBO 수렴 곡선을 확인함.
    - 스팟당 세포 유형별 기여도를 합리적으로 추정하였으나, 절대 세포 수와 상대 비율 간의 비교 분석이 다소 정성적임.
  - **Pass (60~79점)**:
    - 단순 OLS 또는 기본 NNLS 패키지를 사용하여 다중공선성에 의한 가중치 진동을 방어하지 못하거나, 사전분포 하이퍼파라미터 튜닝 없이 기본값만 실행함.
  - **Fail (0~59점)**:
    - 레퍼런스 시그니처 정렬에 실패하여 대다수 유전자가 누락되었거나, 모델 학습이 발산하여 모든 스팟의 세포 추정치가 0으로 수렴함.

---

### Stage 4: 공간 가변 유전자 및 미세환경 니치 군집화 (SVGs & Cellular Niches)

- **평가 초점**: 공간 좌표를 보존하는 통계 검정을 통해 형태학적 구배를 규명하고, 국소 세포 조성 기반의 기능적 니치를 분할하였는가?
- **평가 척도**:
  - **Distinction (90~100점)**:
    - `SPARK-X` 또는 `Squidpy Moran's I`를 적용하여 다중 가설 검정(FDR $< 0.01$)을 통과한 핵심 SVG를 체계적으로 도출함.
    - Moran's 산점도 4사분면 분석을 통해 핫스팟과 경계부 침습 전선 특이 마커를 기능적으로 분류함.
    - 국소 윈도우 세포 구성비 행렬($N \times K$)을 구축하여 전립선암 4대 핵심 니치(선암 코어, CAF 피막, TLS, 골전이 침습 전선)를 명확히 분리하고 이웃 풍부도 검정($Z$-점수)으로 검증함.
  - **Merit (80~89점)**:
    - Moran's I 통계량을 성공적으로 계산하고 상위 SVG의 공간 발현 맵을 시각화함.
    - 세포 니치 클러스터링을 수행하여 주요 종양 영역과 기질 영역을 구분하였으나, TLS 등 미세 니치의 분별력이 다소 미흡함.
  - **Pass (60~79점)**:
    - 공간 좌표를 고려하지 않고 전사체 발현만으로 단순 Leiden 클러스터링을 수행한 후 이를 공간 니치라고 오판함.
  - **Fail (0~59점)**:
    - 공간 자기상관 검정을 수행하지 못했거나, 군집화 결과가 조직학적 병리 구조와 완전히 괴리됨.

---

### Stage 5: 거리 제약 리간드-수용체 신호 및 임상 중개 (Spatial Signaling & Clinical Translation)

- **평가 초점**: 생물리학적 확산 한계를 반영하여 세포 간 통신을 모델링하고, 치료 저항성을 극복할 수 있는 실행 가능한 임상 가설을 도출하였는가?
- **평가 척도**:
  - **Distinction (90~100점)**:
    - `COMMOT` 또는 거리 감쇄 함수($d_0 = 100 \sim 200\ \mu\text{m}$) 기반의 최적 수송 모델을 적용하여 허위 양성 신호를 완벽히 차단함.
    - 신호 플럭스 벡터장($\vec{v}$)을 시각화하여 종양-기질 경계부의 면역 배제 기전(예: CAF 유래 CXCL12/TGF-$\beta$ 와류) 또는 골전이 용골성 악순환(RANKL-RANK)을 규명함.
    - 관측된 공간 신호 축을 무력화하기 위한 구체적인 병용 표적 치료제(예: CXCR4 저해제 + 데노수맙 + 면역항암제)와 임상 바이오마커를 논리적으로 제시함.
  - **Merit (80~89점)**:
    - 거리 제약 조건 하에서 유의미한 리간드-수용체 상호작용 쌍을 도출하고 이웃 풍부도와의 연관성을 논의함.
    - 신호 플럭스의 방향성 해석이 포함되었으나 임상 치료 전략과의 연계가 일반론적인 수준에 머무름.
  - **Pass (60~79점)**:
    - 공간 거리를 고려하지 않은 단순 CellPhoneDB 점수를 나열하고, 물리적으로 격리된 세포 간 통신을 배제하지 못함.
  - **Fail (0~59점)**:
    - 리간드-수용체 상호작용 분석을 수행하지 못했거나 도출된 결과가 생물학적으로 불가능한 오류를 포함함.

---

## 4. 실전 파이썬 워크플로우: 품질 제어 및 아티팩트 검출 진단 스크립트

다음은 공간전사체 원시 객체에서 4대 아티팩트 지표를 자동으로 스캔하고 불량 스팟을 플래깅하는 진단 파이프라인입니다:

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def run_spatial_qc_diagnostics(adata, neg_probe_prefix="NegControl", max_mt=20.0):
    """
    공간전사체 슬라이드의 4대 아티팩트(노이즈, 주름, 확산, 미토콘드리아)를 자동 진단합니다.
    """
    # 1. 라이브러리 크기 및 미토콘드리아 비율 산출
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # 2. 음성 대조 프로브 노이즈 비율 진단
    neg_genes = [g for g in adata.var_names if g.startswith(neg_probe_prefix)]
    if len(neg_genes) > 0:
        neg_counts = np.sum(adata[:, neg_genes].X, axis=1)
        if hasattr(neg_counts, "A1"):
            neg_counts = neg_counts.A1
        adata.obs["neg_probe_ratio"] = (neg_counts / adata.obs["total_counts"]) * 100
        adata.obs["flag_high_noise"] = adata.obs["neg_probe_ratio"] > 1.0
    else:
        adata.obs["flag_high_noise"] = False

    # 3. 조직 주름 및 중첩(Folding) 이상치 진단 (Local UMI Outlier, z > 3.5)
    mean_counts = adata.obs["total_counts"].mean()
    std_counts = adata.obs["total_counts"].std()
    adata.obs["flag_folded_tissue"] = adata.obs["total_counts"] > (mean_counts + 3.5 * std_counts)
    
    # 4. 미토콘드리아 과다(투과화 불량/괴사) 진단
    adata.obs["flag_high_mt"] = adata.obs["pct_counts_mt"] > max_mt
    
    # 5. 종합 불량 스팟 플래그 (QC Fail Mask)
    adata.obs["qc_fail"] = (
        adata.obs["flag_high_noise"] | 
        adata.obs["flag_folded_tissue"] | 
        adata.obs["flag_high_mt"]
    )
    
    print("=== Spatial Transcriptomics QC Diagnostic Report ===")
    print(f"Total Spots Analyzed: {adata.n_obs}")
    print(f"High Background Noise Spots: {adata.obs['flag_high_noise'].sum()}")
    print(f"Folded/Double-layered Spots: {adata.obs['flag_folded_tissue'].sum()}")
    print(f"High Mitochondrial Spots: {adata.obs['flag_high_mt'].sum()}")
    print(f"Total QC Filtered Spots: {adata.obs['qc_fail'].sum()} ({adata.obs['qc_fail'].mean()*100:.2f}%)")
    
    return adata

# 실행 및 시각화
adata = sc.read_h5ad("data/visium_prostate_raw.h5ad")
adata = run_spatial_qc_diagnostics(adata)

# 이상치 스팟 공간 매핑 시각화
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.spatial(
        adata,
        color="qc_fail",
        palette={False: "lightgray", True: "crimson"},
        size=1.3,
        title="Artifact & QC Failure Spots"
    )
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
10x Visium 실험에서 조직 절편 효소 투과화(Permeabilization) 시간이 과도하게 길어졌을 때 발생하는 'mRNA 측면 확산(Lateral Diffusion)' 현상의 생물리학적 발생 기전과, 이로 인해 스팟 디콘볼루션 결과에서 나타나는 치명적 아티팩트(인접 조직으로의 시그니처 번짐)를 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **생물리학적 발생 기전**:
  - 단백질분해효소(Pepsin)가 조직의 세포막과 세포외 기질 단백질을 장시간 분해하면, mRNA를 세포 내에 고착시켜 주던 세포골격망(Cytoskeleton)이 완전히 액화(Liquefaction)됩니다.
  - 포획 슬라이드 표면의 올리고-dT 프로브와 결합하기 전에 유리된 mRNA 분자들이 수직으로 침강하지 못하고, 조직액 내부에서 픽의 확산 법칙에 따라 측면(Lateral) 방향으로 무작위 열운동(Brownian Motion)을 하며 인접 영역으로 이동합니다.
- **디콘볼루션 결과에서의 치명적 아티팩트**:
  - 종양 결절 중심부(예: AR+ 선암)에서만 고발현되어야 할 고유 유전자(*KLK3, TMPRSS2*)가 인접 스팟(Visium 중심 간격 $100\ \mu\text{m}$) 및 그 외곽의 정상 결합조직이나 혈관 주변부 스팟으로 흘러들어갑니다(Bleeding / Halo Effect). (이 문서군의 확산 거리는 시나리오별 **예시 가정**이다: 01장 $20 \sim 30\ \mu\text{m}$, 03장 $\lambda_d = 15\ \mu\text{m}$, 09장 $50 \sim 200\ \mu\text{m}$이며 실제 값은 Tissue Optimization 등 실측으로 정한다.)
  - 결과적으로 디콘볼루션 알고리즘(`cell2location` 등)은 이 번진 mRNA를 실제 종양 세포가 침윤한 것으로 오판하여, 순수한 기질이나 정상 지방 조직 구역에 'AR+ 암세포'가 광범위하게 존재하는 것으로 거짓 양성 세포 수를 할당함으로써 종양 침습 전선의 경계가 모호해지는 치명적 분석 오류가 발생합니다.

</details>

---

### Q2. [알고리즘 및 통계학적 진단 (Statistical Diagnostics)]
snRNA-seq 기준 데이터에 매우 밀접한 두 세포 계통(예: 성숙 관강 상피세포 Luminal AR+ 와 치료 저항성 중간 가소성 세포 Amphicrine)이 포함되어 있어 분산 팽창 계수 $\text{VIF} > 80$에 달하는 심각한 다중공선성이 존재할 때, 디콘볼루션 가중치 추정치가 불안정해지는 이유를 설명하고 이를 통계학적으로 완화하기 위한 구체적인 전처리 해결책 2가지를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **추정치 불안정의 수학적 이유**:
  - 두 세포 유형의 유전자 시그니처 상관계수가 극도로 높으면, 시그니처 행렬 $G$의 열들이 거의 평행해져 최소 제곱 공분산 행렬 $(G G^T)$의 최소 고유값이 0에 수렴합니다.
  - 가중치 추정량의 분산이 폭증하므로, 실제로는 두 세포 유형이 스팟 내에 안정적으로 5:5 비율로 존재하더라도 알고리즘은 하나의 세포를 10으로 극대화하고 다른 하나를 0으로 강제 절단하는 등의 격렬한 수치적 진동(Numerical Instability)을 일으킵니다.
- **구체적인 전처리 해결책 2가지**:
  - **1. 계통학적 시그니처 조기 병합(Coarse-graining / Hierarchical Merging)**:
    - 상관계수가 높은(예시 기준 $r \gtrsim 0.9$) 두 미세 아형을 별도로 분리하여 디콘볼루션하지 않고, Stage 1 단계에서 상위 개념인 'AR-active Plastic Epithelium' 단일 클러스터로 병합하여 레퍼런스 시그니처 행렬의 조건수를 낮춘 후 디콘볼루션을 수행합니다.
  - **2. 공통 유전자 배제 및 계통 특이적 독점 마커 유전자 기반 서브셋 추출**:
    - 두 세포 유형 간에 발현량이 동일한 수천 개의 보존성 하우스키핑 유전자를 시그니처 매트릭스에서 제거하고, 두 클러스터를 통계적으로 유의미하게 구별해 주는 핵심 차등 발현 유전자(예: Amphicrine 고유의 *CHGA, SYP* 및 Luminal 고유의 *KLK3, SPINK1*)만을 선택하여 시그니처 행렬의 직교성(Orthogonality)을 회복시킵니다.

</details>

---

### Q3. [캡스톤 루브릭 기반 종합 임상 평가 (Comprehensive Clinical Evaluation)]
본 단원의 5단계 캡스톤 루브릭을 기준으로, 한 학생 연구팀이 제출한 다음 mCRPC 골전이 공간전사체 연구 보고서의 각 Stage별 결함을 진단하고, 해당 보고서가 'Fail'을 면하고 'Distinction'으로 도약하기 위해 수정해야 할 핵심 사항을 단계별로 자문하시오:

```
[학생 연구팀 보고서 요약]
- Stage 1: UMI 수가 100 미만인 스팟만 제거하고, H&E 상 심하게 접힌 주름 구역을 그대로 분석에 포함함.
- Stage 2: 고해상도 Xenium 슬라이드에서 DAPI 핵 이미지를 보지 않고 스팟 그리드 단위(100 µm)로 임의 분할함.
- Stage 3: 기본 NNLS 알고리즘을 사용하여 세포 '비율(Proportions)'만을 산출하고 골전이 코어의 면역 밀도를 논함.
- Stage 4: 공간 좌표를 배제한 채 전사체 발현만으로 UMAP을 그리고 10개 클러스터로 나눈 뒤 이를 니치로 명명함.
- Stage 5: CellPhoneDB로 전체 세포 간 통신을 돌려 5 mm 떨어진 암세포와 파골세포 간 RANKL-RANK가 최고점이라고 보고함.
```

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **Stage 1 (QC 결함) $\to$ 수정 자문**:
  - **결함**: UMI 100 미만 기준은 지나치게 낮아 빈 공간 노이즈를 포함하며, 주름 구역 방치로 단위 면적당 카운트가 2배 폭증한 아티팩트가 생김.
  - **개선**: 음성 대조 프로브 비율($< 1.0\%$)과 미토콘드리아($< 20\%$)를 기준으로 삼고, 컴퓨터 비전 또는 수동 마스킹으로 주름 부위 스팟을 완전히 배제해야 함.
- **Stage 2 (세포 분할 결함) $\to$ 수정 자문**:
  - **결함**: 서브마이크론 해상도 Xenium 데이터를 $100\ \mu\text{m}$ 그리드로 다운샘플링하여 단일 세포 해상도의 이점을 완전히 상실함.
  - **개선**: `Baysor` 또는 `Mesmer`를 DAPI 및 세포막 염색 이미지에 적용하여 정밀 단일 세포 경계를 복원하고 세포당 전사체 할당을 수행해야 함.
- **Stage 3 (디콘볼루션 결함) $\to$ 수정 자문**:
  - **결함**: NNLS는 다중공선성에 취약하며, 단순 '비율'을 사용하여 밀도가 매우 높은 골전이 코어와 희소 기질 간의 세포 수 왜곡(Compositional Illusion)을 유발함.
  - **개선**: H&E 평균 핵 개수를 사전분포로 주입하는 `cell2location`을 적용하여 스팟당 절대 세포 수(q05 Abundance)를 복원해야 함.
- **Stage 4 (SVG 및 니치 결함) $\to$ 수정 자문**:
  - **결함**: 공간 좌표가 없는 UMAP 클러스터는 공간적 니치가 아님.
  - **개선**: `SPARK-X`/`Moran's I`로 공간 가변 유전자를 선별하고, 국소 윈도우($k=30$) 세포 구성비 행렬을 구축하여 전립선암 4대 기능적 니치를 분할해야 함.
- **Stage 5 (신호 분석 결함) $\to$ 수정 자문**:
  - **결함**: $5\ \text{mm}$ 거리는 분비 단백질 확산 한계($200\ \mu\text{m}$)를 25배 이상 초과한 완벽한 허위 양성(False Positive)임.
  - **개선**: 거리 감쇄 제약($d \le 200\ \mu\text{m}$)을 부과하는 `COMMOT`을 적용하여 암세포-파골세포 접촉 침습 전선에서의 실제 국소 유선 벡터장만을 유효 신호로 입증해야 함.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 검수 근거 문헌
1. He S, Bhatt R, Brown C, et al. High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging. *Nat Biotechnol*. 2022;40(12):1794-1806. PMID: 36203011. doi: 10.1038/s41587-022-01483-z
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
<!-- REFERENCES:END -->

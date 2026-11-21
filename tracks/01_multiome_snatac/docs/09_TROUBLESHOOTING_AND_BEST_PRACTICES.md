# Chapter 09. 멀티오믹스 실전 트러블슈팅 및 5단계 캡스톤 프로젝트 루브릭

> **Cognitive Goal**: 멀티오믹스(snRNA-seq + snATAC-seq) 연구 현장에서 빈번히 발생하는 4대 치명적 기술 결함(미토콘드리아 오염, TSS 농축도 붕괴, LSI 배치 아티팩트, 골전이 조직 탈회 손상)의 분자생물학적 메커니즘을 규명하고 즉각적인 트러블슈팅 솔루션을 도출한다. 아울러 원시 데이터 품질 관리부터 삼중 eRegulon 전사조절망 복원 및 재현 가능한 파이프라인 패키징에 이르는 '5단계 멀티오믹스 캡스톤 프로젝트 완주 루브릭'을 체계적으로 정립한다.

---

## 1. 서론: 왜 멀티오믹스 실험은 실패 확률이 두 배 높은가?

단일 세포 전사체(sc/snRNA-seq)와 후성유전체(snATAC-seq)를 단일 핵 내부에서 동시에 포획하는 멀티오믹스(Multiome) 플랫폼은 세포의 유전형-표현형 연결고리를 규명하는 가장 강력한 무기입니다. 그러나 두 가지 서로 다른 생화학 반응—**mRNA 역전사(Reverse Transcription)**와 **Tn5 트랜스포사제 전치 반응(Tagmentation)**—을 단일 반응 챔버(GEM) 안에서 완벽히 양립시켜야 하므로, 실험적 실패율과 기술적 노이즈가 단일 오믹스에 비해 기하급수적으로 증가합니다.

핵막의 미세한 손상은 snATAC의 크로마틴 구조를 붕괴시키고, 세척(Wash)의 불완전성은 미토콘드리아 오염을 유발하며, 조직 파쇄 시의 기계적 전단력은 뉴클레오솜 래더(Nucleosome ladder; 약 200 bp 주기)를 완전히 지워버립니다. 따라서 고급 생물정보학자는 원시 FASTQ부터 최종 조절망까지 각 단계의 잠재적 위험 요소를 사전에 감지하고 구제(Rescue)할 수 있는 심층적 진단 역량을 갖추어야 합니다.

---

## 2. 현장 중심 멀티오믹스 4대 치명적 결함 및 트러블슈팅 가이드

```
+---------------------------------------------------------------------------------------------------+
|                           멀티오믹스 4대 실패 모드 및 대응 아키텍처                               |
+---------------------------------------------------------------------------------------------------+
|  [Issue 1: 미토콘드리아 오염]  chrM 리드 > 20%  ──► Lysis 세척 실패, OptiPrep 밀도구배 원심분리   |
|  [Issue 2: TSS 농축도 붕괴]   TSS Score < 4.0   ──► Dounce 과파쇄, 저온 제어 및 핵 형태학 검경   |
|  [Issue 3: LSI 배치 왜곡]     LSI 1 심도 상관성 ──► LSI 1 배제, LSI 2-30 Harmony & iLISI 검증    |
|  [Issue 4: 골전이 탈회 손상]  DNA 가수분해/가교 ──► 강산 배제, 4℃ 0.5M EDTA 온화 탈회 프로토콜 |
+---------------------------------------------------------------------------------------------------+
```

### 2.1 Issue 1: snATAC에서 비정상적인 미토콘드리아 DNA 분율 ($> 20\%$)

#### 1) 분자적 근본 원인 (Root Cause)
- 본래 세포 전체를 사용하는 scATAC-seq과 달리, **단일 핵 ATAC(snATAC-seq)**은 세포막을 파쇄(Cell lysis)하여 세포질을 제거하므로 이론적으로 미토콘드리아 DNA($\text{mtDNA}$, `chrM`)가 거의 검출되지 않아야 합니다 ($\le 2 \sim 5\%$).
- 그러나 용해 완충액(Lysis buffer)의 계면활성제 강도가 부적절하여 세포질 막이 불완전하게 제거되었거나, 원심분리 후 세척(Wash) 단계가 불충분할 경우 파괴된 미토콘드리아 막 잔해(Mitochondrial debris)가 핵막 외벽에 비특이적으로 달라붙게 됩니다.
- 노출된 순수 원형 이중나선인 mtDNA는 뉴클레오솜으로 보호되지 않으므로, Tn5 효소의 무제한 표적이 되어 엄청난 양의 tagmentation 단편을 생성합니다.

#### 2) 데이터에 미치는 악영향
- 시퀀싱 리드의 $20 \sim 50\%$가 불필요한 미토콘드리아 게놈(16.5 kb)에 집중되어 핵 게놈(Nuclear genome, 3.2 Gb)의 시퀀싱 깊이(Sequencing Depth)가 급감합니다.
- 세포당 유효 핵 게놈 피크 단편(FRiP)이 바닥나 전사 인자 모티프 탐색 및 원거리 인핸서 공접근성 검출력이 붕괴됩니다.

#### 3) 진단 기준 및 해결책
$$f_{\text{mt}} = \frac{N_{\text{fragments, chrM}}}{N_{\text{total fragments}}} \times 100\%$$
- **경고 기준**: $f_{\text{mt}} > 5\%$, **심각한 결함**: $f_{\text{mt}} > 20\%$
- **Wet-lab 해결책**:
  1. **계면활성제 미세 적정(Lysis Buffer Titration)**: IGEPAL CA-630(0.1%), Tween-20(0.1%), Digitonin(0.01%)의 농도와 용해 시간을 조직 종류별로 단계별 최적화합니다.
  2. **OptiPrep 밀도 구배 원심분리(Iodixanol Density Gradient)**: $29\% / 35\%$ OptiPrep 불연속 구배를 통과시켜 밀도가 낮은 미토콘드리아 막 잔해와 세포 파편을 상층액으로 분리하고 순수한 Intact Nuclei 펠릿만을 회수합니다.
  3. **형광 검경**: DAPI 염색 후 공초점 현미경으로 핵막 외곽에 부착된 세포질/미토콘드리아 형광 잔해가 없는지 사전 확인합니다.

---

### 2.2 Issue 2: TSS Enrichment Score 붕괴 ($< 4$) 및 뉴클레오솜 래더(Nucleosome Ladder) 소실

#### 1) 분자적 근본 원인 (Root Cause)
- **TSS Enrichment Score**는 활성 유전자의 전사시작점(TSS) 주변 열린 크로마틴 영역에 Tn5 삽입이 얼마나 특이적으로 집중되었는지를 측정하는 신호 대 잡음비(Signal-to-Noise Ratio) 지표입니다.
- 조직 분쇄(Dounce homogenization) 과정에서 과도한 물리적 전단력(Over-shearing), 마찰열 발생, 혹은 검체의 반복적인 동결-해동(Freeze-thaw cycles)은 뉴클레오솜에 감겨 있던 닫힌 염색질(Heterochromatin)을 강제로 풀어헤치고 DNA 가닥을 무작위로 절단합니다.
- 이 상태에서 Tn5를 투입하면 프로모터/인핸서 특이성 없이 게놈 전역에 무작위로 삽입되는 **과전치(Over-tagmentation) 및 무작위 백그라운드 노이즈**가 발생합니다.

#### 2) 데이터에 미치는 악영향
- 유전자의 프로모터 피크가 검출되지 않고 시퀀싱 커버리지가 게놈 전체에 평평하게 퍼집니다.
- 단편 크기 분포(Fragment Size Distribution)에서 단일 뉴클레오솜(약 $180\sim220\text{ bp}$), 디뉴클레오솜(약 $380\sim420\text{ bp}$)의 뚜렷한 주기성 피크가 사라지고, 단편 크기가 $100\text{ bp}$ 미만으로 완전히 부서진 단일 붕괴 분포를 보입니다.

#### 3) 진단 기준 및 해결책
$$\text{TSS Enrichment} = \frac{\text{Mean normalized coverage at TSS center } (\pm 50\text{ bp})}{\text{Mean normalized coverage at flanking flanks } (\pm 1900 \sim 2000\text{ bp})}$$

| TSS Score 등급 | 기준 수치 | 품질 평가 및 조치 |
| :--- | :--- | :--- |
| **Fail (실패)** | $\text{TSS} < 4.0$ | 분석 불가 (백그라운드 노이즈 지배, 재실험 필수) |
| **Acceptable (통과)** | $4.0 \le \text{TSS} < 6.0$ | 기본 분석 가능 (엄격한 피크 필터링 필요) |
| **Good (우수)** | $6.0 \le \text{TSS} < 10.0$ | 표준 멀티오믹스 분석 최적 |
| **Ideal (최우수)** | $\text{TSS} \ge 10.0$ | 고해상도 TF 풋프린팅 및 SCENIC+ 분석 최적 |

- **Wet-lab 해결책**:
  1. **Dounce 스트로크 수 감소 및 완충**: 조직 분쇄 시 느슨한 Pestle A(10~15회) 후 팽팽한 Pestle B는 5~10회 이내로 최소화합니다.
  2. **철저한 $4^\circ\text{C}$ 저온 유지**: 분쇄용 튜브, 완충액, 원심분리기를 사전 냉각하고 얼음 위에서 신속히 조작하여 내인성 뉴클레아제(DNase) 활성을 원천 억제합니다.
  3. **핵 형태학(Nuclear Morphology) 품질 검증**: 현미경 하에서 핵이 둥글고 매끄러운 원형을 유지하는지 확인하며, 핵막이 쭈그러들거나 블레빙(Blebbing)이 관찰되면 즉시 분쇄 강도를 낮춥니다.

---

### 2.3 Issue 3: snATAC 잠재 시맨틱 인덱싱(LSI) 공간에서의 심각한 배치 효과 (Batch Effects)

#### 1) 통계적/생물학적 근본 원인 (Root Cause)
- snATAC-seq의 원시 카운트 매트릭스는 $98\%$ 이상이 0으로 채워진 극도의 희소(Sparse) 이진 행렬입니다. 이를 다루기 위해 자연어 처리 기법인 **TF-IDF 정규화 후 특이값 분해(SVD)**를 수행하는 **LSI(Latent Semantic Indexing)**를 표준 차원 축소로 사용합니다.
- 그러나 LSI의 첫 번째 성분(**LSI Component 1**)은 거의 전적으로 각 단일 핵의 **총 시퀀싱 깊이(Total Sequencing Depth / Total Fragments per Nucleus)**와 극단적인 양의 상관관계($r > 0.90$)를 갖습니다:
$$r(\text{LSI}_1, \log_{10}(\text{Total Fragments})) \ge 0.90$$
- 서로 다른 실험 날짜, 라이브러리 제작 배치, Tn5 로트(Lot) 번호에 따라 평균 시퀀싱 깊이가 다르면, LSI 공간에서 세포들이 생물학적 계통이 아니라 실험 배치별로 완전히 갈라지게 됩니다.

#### 2) 분석 시의 치명적 함정
- LSI 1을 제거하지 않으면 깊이 편차가 배치 구조와 섞여 보정 알고리즘(Harmony, Scanorama, BBKNN)이 깊이 신호를 무리하게 맞추려 할 수 있습니다. 한편 **과보정(Over-correction)**은 배치와 생물학적 변이가 교락된 경우(예: 환자=배치) 강하게 보정할 때 일어나며, 이는 LSI 1 포함 여부와 별개의 문제입니다(LSI 1 제외는 깊이 교란을 줄이는 독립 단계, Signac `dims = 2:30`).

#### 3) 해결책 및 벤치마크 워크플로우
1. **LSI 1 성분 배제**: 다운스트림 UMAP 차원 축소 및 이웃 그래프 구축 시 LSI 1을 강제로 제외하고 **LSI 2~30 성분만 활용**합니다.
2. **LSI 2~30 공간에 Harmony 적용**:
   - `sc.external.pp.harmony_integrate(adata, key='batch', basis='X_lsi', adjusted_basis='X_lsi_harmony')`
3. **통합 품질 정량화 (LISI 지표 산출)**:
   - **iLISI (integration LISI)**: 배치 간 혼합도로서, 동일 이웃 내에 서로 다른 배치가 얼마나 골고루 섞였는지 측정. 값의 범위는 $1 \sim B$($B$=배치 수)이며 $B$에 가까울수록 잘 섞인 것입니다(배치가 2개일 때만 목표치가 2.0에 해당). 배치 수가 다른 데이터 간 비교에는 0~1로 정규화한 iLISI를 사용하며, 정규화 값은 1에 가까울수록 좋습니다.
   - **cLISI (cell-type LISI)**: 세포 유형 분리도로서, 동일 이웃 내에 서로 다른 세포 유형이 섞이지 않고 잘 보존되었는지 측정 ($\text{목표치} \to 1.0$).

---

### 2.4 Issue 4: 진행성 전립선암(mCRPC) 골전이(Bone Metastasis) 조직의 탈회(Decalcification) 손상

#### 1) 분자적 근본 원인 (Root Cause)
- 전립선암 환자의 원격 전이 중 $80 \sim 90\%$가 골격계(Bone)로 침윤하며, 이 골전이 병변은 치명적인 치료 저항성이 발현되는 핵심 부위입니다.
- 그러나 뼈 조직 검체는 단단한 수산화인회석(Hydroxyapatite) 무기질 칼슘 결정으로 둘러싸여 있어, 전통적인 병리학적 절편 제작을 위해 **강산(Strong acids: 포름산, 염산, 질산 등)**을 처리하는 급속 탈회 과정을 거칩니다.
- 강산에 노출된 조직은 pH가 1~2 수준으로 급격히 떨어지며, DNA의 퓨린 염기 탈락(Depurination) 및 인산디에스테르 골격의 비가역적 가수분해(Hydrolysis)가 일어나 **게놈 DNA가 산산조각 파편화**됩니다. 또한 산성 환경은 포르말린 고정 시의 크로마틴 가교를 왜곡하여 Tn5 효소의 접근을 방해할 수 있습니다.

#### 2) 해결책: 온화한 저온 EDTA 탈회 프로토콜
골전이 임상 검체에서 강산 탈회는 핵산 손상이 크므로 피하고, 아래와 같은 온화한 EDTA 탈회를 대안으로 고려할 수 있습니다. 다만 이는 이 교재가 제시하는 가정적 프로토콜이며 Multiome(RNA+ATAC)에 대한 검증된 출처는 제시하지 못했습니다. 일반적으로 Multiome은 고정하지 않은 동결 조직에서 핵을 분리하며, 수일간의 EDTA 처리는 RNA 무결성을 해칠 수 있고 잔류 EDTA는 Mg²⁺ 의존적인 Tn5 반응을 억제하므로 철저한 세척과 RNA 품질(RIN/DV200) 확인이 필요합니다:
1. **중성 완충 킬레이트제 사용**: 강산 대신 칼슘 이온만을 선택적으로 포획하는 **$0.5\text{ M}$ EDTA ($\text{pH } 7.4 \sim 8.0$)** 용액을 사용합니다.
2. **저온 장시간 탈회**: DNA 단편화와 단백질 분해 효소의 활성을 차단하기 위해 반드시 **$4^\circ\text{C}$ 저온 회전 교반기**에서 $3 \sim 7\text{일}$ 동안 천천히 탈회합니다.
3. **골용해성 병변 미세 절제(Osteolytic Microdissection)**: 치밀골(Cortical bone) 전체를 처리하지 않고, CT/방사선 가이드 하에 종양 세포가 집중된 골용해성 연조직 부위만을 큐렛(Curette) 또는 레이저 캡처로 긁어내어 탈회 시간 자체를 최소화합니다.

---

## 3. 멀티오믹스 5단계 캡스톤 프로젝트 완주 루브릭 (Capstone Rubric)

단일 세포 멀티오믹스(snRNA + snATAC) 연구의 전 과정을 완성도 높은 논문급 포트폴리오로 승화시키기 위해, 전산 및 분자생물학적 검증 요건을 5단계로 체계화한 평가 루브릭입니다.

```
[ Stage 1: Dual-modal QC ]  ──►  [ Stage 2: Peak & LSI ]  ──►  [ Stage 3: WNN Integration ]
       │                                 │                                    │
       ▼                                 ▼                                    ▼
- snRNA Mito < 10%              - MACS3 Pseudobulk Peak              - MuData WNN 통합 그래프
- snATAC chrM < 5%              - TF-IDF + SVD                       - RNA/ATAC 가중치 산출
- TSS > 6.0, nucleosome ladder      - LSI 1 심도 상관성 배제             - 통합 클러스터링 최적화
- AMULET 더블렛 제거            - Harmony 배치 보정                  - 모달리티 기여도 검증
                                                                              │
                                                                              ▼
[ Stage 5: eGRN & Packaging ] ◄── [ Stage 4: Peak-to-Gene & Enhancer Dissection ]
       │                                 │
       ▼                                 ▼
- pycisTopic 토픽 모델링        - Cicero 정규화 그래피컬 라쏘
- pycistarget 모티프 탐색       - 순열 검정 기반 P2G 링크 확립
- Triplet eRegulon 복원         - AR 650 kb 원거리 인핸서 분석
- Nextflow/Snakemake 패키징     - 3D 루핑 연계 생물학적 해석
```

### [표] 멀티오믹스 캡스톤 프로젝트 5단계 상세 평가 루브릭

| 단계 (Stage) | 핵심 평가 항목 | 최우수 달성 기준 (Exemplary, 100점) | 보통 (Proficient, 80점) | 미흡 및 결함 (Needs Work, <60점) |
| :--- | :--- | :--- | :--- | :--- |
| **Stage 1: 이중 모달리티 품질 관리 (Dual-modal QC)** | • snRNA 품질<br>• snATAC 품질<br>• 이중체(Doublet) 제거 | • RNA mito $< 10\%$, ATAC chrM $< 5\%$ 달성<br>• TSS Score $\ge 6.0$(03장 등급표 우수 이상) 및 선명한 ~200 bp 뉴클레오솜 래더 확인<br>• `AMULET` 및 `Scrublet`으로 동시 이중체 완벽 제거 | • RNA mito $< 15\%$, ATAC chrM $< 10\%$<br>• TSS Score $4.0 \sim 5.9$<br>• 단일 모달리티 기반 이중체만 필터링 | • ATAC chrM $> 20\%$ 방치<br>• TSS Score $< 4.0$ 방치<br>• 이중체 필터링 생략 |
| **Stage 2: 피크 콜링 및 LSI 차원 축소** | • MACS3 피크 동정<br>• LSI 차원 축소<br>• LSI 1 심도 보정 | • 세포 유형별 Pseudobulk 분할 후 MACS3 피크 콜링<br>• 재현 가능한 Consensus Peakset ($150\text{k} \sim 250\text{k}$) 수립<br>• LSI 1 깊이 상관성($|r| > 0.75$, ArchR `corCutOff` 기본값) 확인 후 LSI 1 제외, Harmony 보정 완료 | • 전체 세포 통합 1회 MACS3 피크 콜링<br>• LSI 1 상관성은 확인했으나 LSI 1을 포함하여 보정 진행 | • 10x 기본 피크셋을 무검증 사용<br>• LSI 1 상관성 검증 생략<br>• 배치 효과 미보정 |
| **Stage 3: WNN 멀티모달 통합 및 클러스터링** | • MuData 데이터 구조화<br>• WNN 친화도 가중치 산출<br>• 모달리티 기여도 평가 | • `MuData` 기반 표준 컨테이너 구축<br>• RNA/ATAC 각각의 최근접 이웃 예측 오차 기반 WNN 가중치 산출<br>• 희귀 아형(t-NEPC 등)에서 특정 모달리티 기여도 상승 해석 완료 | • WNN 통합을 수행했으나 가중치 분포 분석 누락<br>• 단순히 두 임베딩을 물리적으로 결합(Concatenation)하여 클러스터링 | • WNN을 적용하지 않고 단일 RNA UMAP에 의존<br>• 두 모달리티 간 클러스터 불일치 방치 |
| **Stage 4: Peak-to-Gene 링크 및 원거리 인핸서 규명** | • Cicero 공접근성<br>• P2G 유의성 검정<br>• *AR* 원거리 인핸서 추적 | • Regularized Graphical Lasso 기반 거리 감쇄 보정 공접근성 산출<br>• 배경 피크 순열 검정 기반 FDR로 유의 P2G 링크 확립(06장 예시 기준 FDR $< 0.05$)<br>• *AR* 유전자 $650\text{ kb}$ 원거리 인핸서 중복 영역의 활성화 링크 증명 | • 단순 피어슨 상관계수만으로 P2G 링크 산출<br>• 순열 검정 생략<br>• 거리 감쇄 보정 미비 | • 프로모터 $\pm 2\text{ kb}$ 영역만 분석하고 원거리 인핸서 배제<br>• 통계적 유의성 검정 부재 |
| **Stage 5: eGRN 복원 및 재현 가능한 파이프라인** | • SCENIC+ 파이프라인<br>• Triplet eRegulon 검증<br>• 워크플로우 패키징 | • `pycisTopic` LDA 토픽 모델링 $\to$ `pycistarget` 모티프 탐색 $\to$ GBM 회귀 완주<br>• 4대 핵심 eRegulon(`AR`, `FOXA1`, `ASCL1`, `ETV4`) 단일 핵 활성도 점수화<br>• Nextflow 또는 Snakemake 기반 전 과정 자동화 및 MultiQC 리포트 패키징 | • SCENIC+를 구동했으나 기본 eRegulon 일부 누락<br>• 주피터 노트북 형태로만 스크립트 보존 | • 1세대 SCENIC 결과로 대체<br>• 파이프라인 자동화 및 재현성 검증 실패 |

---

## 4. 프로덕션 레벨 실전 코드: LSI 배치 보정 및 iLISI/cLISI 벤치마크

실제 연구 현장에서 snATAC 시퀀싱 깊이 편차를 극복하고 LSI 1을 올바르게 배제하여 Harmony 배치를 보정하는 Python 파이프라인 스크립트입니다:

```python
import muon as mu
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import harmonypy as hm

def process_snatac_lsi_and_harmony(mdata, batch_key='patient_id', n_comps=30):
    """
    snATAC 데이터에 대해 LSI 차원 축소를 수행하고,
    LSI 1의 시퀀싱 깊이 상관성을 검증한 뒤 LSI 2~30에 대해 Harmony를 실행합니다.
    """
    atac = mdata['atac']
    
    print("[1] Calculating Total Fragments and QC Metrics...")
    if 'total_fragments' not in atac.obs.columns:
        atac.obs['total_fragments'] = np.array(atac.X.sum(axis=1)).flatten()
    
    print("[2] Running TF-IDF and SVD (LSI)...")
    # muon TF-IDF: 기본값은 log(TF) x log(IDF)이며 Signac 기본(log1p(TF x IDF x 1e4))과는 다른 변형이다
    # sc.tl.pca(zero_center=False)는 TruncatedSVD와 유사하나 muon lsi와 달리 성분 표준화가 없다
    mu.atac.pp.tfidf(atac, scale_factor=10000)
    sc.tl.pca(atac, n_comps=n_comps, zero_center=False, svd_solver='arpack')
    atac.obsm['X_lsi'] = atac.obsm['X_pca'].copy()
    
    # LSI 1과 총 단편 수(Total fragments) 사이의 피어슨 상관계수 산출
    lsi1 = atac.obsm['X_lsi'][:, 0]
    depth = np.log10(atac.obs['total_fragments'] + 1)
    corr, pval = pearsonr(lsi1, depth)
    print(f"[*] Pearson correlation between LSI_1 and log10(Total Fragments): {corr:.4f} (p-value: {pval:.2e})")
    
    # LSI 1이 시퀀싱 심도와 강하게 상관된 경우 (|r| > 0.75, ArchR corCutOff 기본값) LSI 1 제외
    if abs(corr) > 0.75:
        print("[!] Strong depth correlation detected. Excluding LSI_1 from batch correction!")
        lsi_for_harmony = atac.obsm['X_lsi'][:, 1:n_comps]  # LSI 2 to n_comps
    else:
        print("[+] Weak depth correlation. Using all LSI components.")
        lsi_for_harmony = atac.obsm['X_lsi'][:, 0:n_comps]
        
    print(f"[3] Running Harmony batch correction across '{batch_key}'...")
    harmony_res = hm.run_harmony(
        lsi_for_harmony,
        atac.obs,
        batch_key,
        max_iter_harmony=20,
        random_state=42
    )
    
    # 보정된 임베딩 저장
    atac.obsm['X_lsi_harmony'] = harmony_res.Z_corr.T
    
    print("[4] Computing Neighbor Graph and UMAP on Harmony-corrected LSI space...")
    sc.pp.neighbors(atac, use_rep='X_lsi_harmony', n_neighbors=30, metric='cosine')
    sc.tl.umap(atac, min_dist=0.3)
    
    print("[✓] snATAC LSI Integration and Harmony Complete.")
    return mdata
```

---

## 5. 트랙 01 완주 총괄 결론 및 거대한 지평의 연결

본 멀티오믹스 트랙([Track 01: Multiome snRNA + snATAC])을 통해 여러분은 현대 단일 세포 생물학의 가장 심오한 영역을 완주하셨습니다:

1. **Chapter 01~04 (전반부)**: 단일 세포에서 단일 핵으로의 전환 패러다임, 미토콘드리아 및 이중체 QC, 10x Multiome GEM 생화학, 그리고 MACS3 기반 정밀 Pseudobulk 피크셋 구축을 달성했습니다.
2. **Chapter 05 (WNN 통합)**: 전사체와 염색질 개방성이 서로 다른 노이즈 수준을 가질 때, 각 모달리티가 이웃 세포 상태를 얼마나 잘 예측하는지(커널 기반 예측 비교)를 통해 두 모달리티를 세포별로 유연하게 가중 결합하는 Weighted Nearest Neighbors 원리를 정립했습니다.
3. **Chapter 06 (원거리 인핸서)**: 선형 게놈의 한계를 깨고 3차원 크로마틴 루핑 모델, 거리 감쇄 함수, Cicero 그래피컬 라쏘를 통해 *AR* 유전자 650 kb 상류의 슈퍼 인핸서 중복을 단일 세포 수준에서 규명했습니다.
4. **Chapter 07 (TF 풋프린팅)**: 단순 개방을 넘어 1 bp 단위의 Tn5 입체 장해(Steric hindrance) 풋프린트를 복원하고, ChromVAR 편차 모델을 통해 전립선암 파이오니어 전사 인자(*FOXA1*, *HOXB13*)의 모티프 접근성 편차(활성도)를 정량화했습니다.
5. **Chapter 08 (SCENIC+ eGRN)**: [TF] ──► [Enhancer] ──► [Target Gene]으로 이어지는 3차원 삼중 조절자 복합체를 복원하여, 거세저항성 전립선암의 4대 악성화 계통 회로(`AR`, `FOXA1`, `ASCL1`, `ETV4`)를 완성했습니다.
6. **Chapter 09 (트러블슈팅 & 루브릭)**: 현장의 4대 기술 결함을 극복하는 생물물리학적 통찰과 5단계 캡스톤 프로젝트 완주 기준을 확보했습니다.

이제 염색질의 문(Gate)을 열어 전사 조절의 인과율을 손에 쥔 여러분은, 세포들이 조직이라는 3차원 물리적 생태계에서 이웃 세포들과 어떻게 대화하고 종양 미세환경(TME)을 구축하는지를 탐구하는 **[Track 02: 공간 전사체학 (Spatial Transcriptomics: Visium & Xenium)]**으로 당당히 나아갈 준비가 되었습니다.

---

## 6. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)] snATAC-seq의 LSI(Latent Semantic Indexing) 차원 축소에서 'LSI Component 1'이 지니는 고유한 수학적/기술적 특성을 설명하고, 배치 보정(Harmony 등)을 적용할 때 LSI 1을 반드시 검증하고 제외해야 하는 이유를 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **LSI Component 1의 기술적 특성**:
  - snATAC-seq 데이터에 TF-IDF(용어 빈도-역문서 빈도) 정규화 후 특이값 분해(SVD)를 적용하는 LSI 알고리즘에서, 제1 특이벡터(Singular Vector 1, LSI 1)는 데이터의 변량(Variance) 중 가장 큰 축을 설명합니다.
  - 그런데 단일 핵 ATAC의 희소 이진 특성상, 가장 큰 변량을 차지하는 요인은 생물학적 세포 유형 차이가 아니라 각 단일 핵 라이브러리가 획득한 **총 시퀀싱 깊이(Total Sequencing Depth / Total Fragment Count per Nucleus)**입니다. 따라서 LSI 1은 단편 수의 로그값($\log_{10}(\text{Depth})$)과 거의 완벽한 양의 상관관계($r > 0.85 \sim 0.95$)를 보입니다.
- **배치 보정 시 LSI 1을 제외해야 하는 이유**:
  - 만약 서로 다른 실험 배치 간에 평균 시퀀싱 깊이 편차가 존재할 때 LSI 1을 포함하여 Harmony나 Scanorama 등의 배치 보정 도구에 투입하면, 알고리즘은 '기술적 깊이 차이'를 '배치 편차'로 오인하여 강제로 세포들을 뒤섞어 버립니다.
  - 이 과정에서 고증식성 암세포(풍부한 개방 염색질과 높은 리드 수)와 정지기/분화 세포(낮은 리드 수) 사이의 **진정한 생물학적 염색질 활성 변이가 완전히 지워지는 치명적인 과보정(Over-correction)**이 발생합니다. 따라서 배치 보정 및 이웃 그래프 구축 시에는 LSI 1을 배제하고 생물학적 구조를 담고 있는 **LSI 2~30 성분만을 선별하여 사용**해야 합니다.

</details>

---

### Q2. [능동 회상 (Active Recall)] mCRPC 환자의 골전이(Bone metastasis) 생검 조직에서 멀티오믹스(snRNA + snATAC) 라이브러리를 제작하고자 할 때, 통상적인 임상 병리 검체 탈회 기법인 '강산(포름산, 질산 등) 급속 탈회'가 초래하는 치명적인 분자생물학적 손상 기전을 서술하고, 이를 대체할 온화한 탈회 프로토콜의 3대 핵심 조건(화합물, pH, 온도)을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **강산 급속 탈회의 치명적 손상 기전**:
  - 뼈 조직의 무기질 칼슘을 신속히 제거하기 위해 사용하는 질산, 염산, 포름산 등은 반응 용액의 pH를 $1.0 \sim 2.0$의 강산성으로 급강하시킵니다.
  - 이 극단적인 산성 환경은 DNA 이중나선의 퓨린 염기를 탈락(Depurination)시키고 **인산디에스테르 결합(Phosphodiester bond)을 무작위로 가수분해(Hydrolysis)**하여 게놈 DNA를 수십 bp 단위로 산산조각 냅니다.
  - 그 결과 온전한 고분자량 크로마틴 구조가 소실되어 snATAC을 위한 뉴클레오솜 래더가 붕괴되고, mRNA 역시 분해되어 역전사 반응이 어려워져 멀티오믹스 라이브러리 제작이 사실상 어려워집니다.
- **대체 온화 탈회 프로토콜의 3대 핵심 조건**:
  1. **화합물**: 강산 대신 칼슘 이온을 온화하게 킬레이트 결합으로 포획하는 **EDTA (Ethylenediaminetetraacetic acid)**를 사용합니다.
  2. **농도 및 pH**: **$0.5\text{ M}$ EDTA**, 생리적 중성 상태인 **$\text{pH } 7.4 \sim 8.0$**으로 완충하여 산에 의한 핵산 가수분해를 방지합니다.
  3. **온도 및 환경**: 내인성 뉴클레아제와 단백질 분해 효소의 활성을 차단하기 위해 반드시 **$4^\circ\text{C}$ 저온**에서 회전 교반(Agitation)하며 수일간 천천히 진행합니다.

</details>

---

### Q3. [임상/트러블슈팅 케이스 스터디 (Troubleshooting Case Study)]

**[임상 및 실험 시나리오]**
귀하의 연구팀은 전이성 거세저항성 전립선암(mCRPC) 환자 12명의 다기관 생검 코호트 검체로 10x Genomics Multiome(snRNA + snATAC) 시퀀싱을 완료하였다. 그러나 1차 파이프라인 분석 보고서에서 다음과 같은 3가지 중대한 기술적 이상 징후가 보고되었다:

1. **환자 A, B 검체**: snATAC 시퀀싱 리드 중 미토콘드리아 게놈(`chrM`) 할당 비율이 **$34\%$**에 달하며, 핵 게놈 피크 영역 단편 비율(FRiP)이 $15\%$ 미만으로 급감함.
2. **환자 C, D 검체**: snATAC의 **TSS Enrichment Score가 $2.6$**으로 추락하였고, 단편 크기 히스토그램에서 단일(~200 bp)·디뉴클레오솜(~400 bp) 피크가 전혀 관찰되지 않고 $<100\text{ bp}$ 영역에 뭉개져 있음.
3. **전체 데이터 통합 단계**: 분석 담당자가 LSI 1부터 30까지 모든 성분을 그대로 투입하여 환자 ID를 기준으로 Harmony 배치 보정을 강하게 수행했더니, UMAP 상에서 12명 환자의 종양 세포가 완벽히 하나로 섞여버렸고 환자 고유의 클론 변이 및 *AR* 인핸서 복제 신호가 통째로 사라짐.

**[질의]**
- (1) 위 3가지 이상 징후 각각의 근본적인 원인(Root Cause)을 실험적/전산학적 관점에서 정확히 진단하시오.
- (2) 이미 시퀀싱이 완료된 데이터 단계에서 수행할 수 있는 '생물정보학적 구제(In silico Rescue) 조치'와, 향후 재실험 시 습식 실험(Wet-lab) 현장에 적용해야 할 '긴급 프로토콜 수정안'을 각각 구체적으로 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

#### (1) 이상 징후별 근본 원인(Root Cause) 진단
1. **환자 A, B (chrM 34%)**:
   - **원인**: 단일 핵 분리 과정에서 세포 용해 완충액(Lysis buffer) 세척 단계가 미흡하여, 파괴된 세포질 내 미토콘드리아 막 잔해가 핵막 외벽에 대량으로 부착된 상태로 Tn5 반응에 투입되었습니다. 비보호 상태의 원형 mtDNA에 Tn5 효소가 집중 결합하여 시퀀싱 리드를 잠식했습니다.
2. **환자 C, D (TSS 2.6, 래더 소실)**:
   - **원인**: 조직 파쇄(Dounce homogenization) 시 과도한 기계적 전단력(Over-shearing, Pestle 스트로크 과다)이 가해졌거나, 검체 채취 및 운송 과정에서 해동-동결이 반복되어 크로마틴 입체 구조가 완전히 파괴되었습니다. 염색질이 풀려 게놈 전역에 Tn5가 무작위 삽입되는 과전치(Over-tagmentation) 및 비특이적 백그라운드 노이즈가 발생했습니다.
3. **전체 데이터 통합 (환자 클론 신호 소실)**:
   - **원인**: 환자 ID를 배치 변수로 Harmony를 강하게 적용하면 환자 고유의 카피수 변이(CNA), *AR* 인핸서 증폭, 종양 특이적 클론 이질성처럼 환자 특이적인 **진정한 생물학적 변이가 배치 효과로 간주되어 제거되는 과보정(Over-correction)**이 일어날 수 있습니다(환자=배치의 교락). LSI 1을 포함했다는 사실은 깊이 교란을 추가하는 별도 문제이며 이 신호 소실의 주원인으로 단정할 수 없습니다. 해결은 theta 하향, 종양 세포는 환자 단위 보정 생략, 보정 전후 비교입니다.

#### (2) 생물정보학적 구제(In silico) 및 습식 실험(Wet-lab) 해결책
- **1. 환자 A, B 구제 및 개선책**:
  - **In silico Rescue**: 분석 파이프라인에서 `chrM` 리드를 완전 필터링한 후, 세포별 최소 핵 게놈 피크 카운트(FRiP $> 0.2$, Non-mito fragments $> 3,000$) 기준을 엄격히 상향 조정하여 온전한 핵만을 다운스트림에 통과시킵니다.
  - **Wet-lab 개선안**: 핵 분리 프로토콜에 **OptiPrep (Iodixanol 29%/35%) 밀도 구배 원심분리**를 즉각 도입하여 미토콘드리아 및 세포질 잔해를 물리적으로 완벽히 분리하고, DAPI 형광 검경으로 핵막 순도를 확인한 후 반응을 개시합니다.
- **2. 환자 C, D 구제 및 개선책**:
  - **In silico Rescue**: TSS Score 2.6 검체는 크로마틴 물리 구조가 손상되었으므로 단독 피크 콜링은 불가능합니다. 타 환자들의 고품질 검체로 생성한 'Consensus Peakset'을 참조 피크로 삼아 카운트 매트릭스만 추출하거나, 신뢰도가 유지된 snRNA-seq 모달리티를 주축으로 삼고 ATAC은 클러스터 수준 pseudobulk 평균치로 제한 활용합니다. (단, TF 풋프린팅이나 P2G 분석에서는 해당 검체를 반드시 제외해야 함).
  - **Wet-lab 개선안**: 신선 동결 검체 프로토콜을 점검하여 해동을 엄격히 방지하고, Dounce Pestle B 스트로크 수를 기존 대비 절반 이하(5~8회)로 줄이며 모든 과정을 $4^\circ\text{C}$ 얼음 위에서 초저온으로 진행합니다.
- **3. 전체 데이터 통합 구제 및 개선책**:
  - **In silico Rescue**: 즉시 기존 Harmony 결과를 파기합니다. SVD 후 LSI 1과 $\log_{10}(\text{Depth})$의 피어슨 상관계수를 확인하여 LSI 1을 강제로 제외하고, **LSI 2~30 성분만을 추출**합니다. 환자 고유의 클론 변이를 보존하기 위해 Harmony의 군집화 다양성 페널티 파라미터(`theta`)를 대폭 낮추거나, 정상 기질/면역 세포에 대해서만 배치를 보정하고 종양 세포는 환자별 고유 클러스터로 분리되도록 파라미터를 재조정합니다. 보정 후 반드시 **iLISI 및 cLISI**를 산출하여(환자 12명을 배치로 두면 iLISI의 범위는 1~12이므로 1.8 같은 값은 거의 보정되지 않은 상태에 해당합니다. 정규화 iLISI(0~1)로 보고하는 것이 해석에 유리합니다. cLISI는 1에 가까울수록 세포 유형이 잘 보존된 것입니다) 생물학적 보존도를 검증합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 검수 근거 문헌
1. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods*. 2019;16(12):1289-1296. PMID: 31740819. doi: 10.1038/s41592-019-0619-0
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- muon 0.2.0 (github.com/scverse/muon). 원 논문: Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. *Genome Biol*. 2022;23(1):42. PMID: 35105358. doi: 10.1186/s13059-021-02577-8
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
- harmonypy 2.1.0 (github.com/slowkow/harmonypy). 원 논문: Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods*. 2019;16(12):1289-1296. PMID: 31740819. doi: 10.1038/s41592-019-0619-0
<!-- REFERENCES:END -->

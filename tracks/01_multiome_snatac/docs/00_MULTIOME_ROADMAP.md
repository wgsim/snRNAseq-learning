# Chapter 00. Track 01 로드맵: 전립선암 10x Multiome (snRNA + snATAC) 학습 과정

> **Cognitive Goal**: 전립선암 후성유전체 연구에서 왜 단일 전사체(snRNA-seq)를 넘어 단일 핵 다중오믹스(10x Multiome: snRNA + snATAC)가 필수적인지 이해하고, 전체 학습 단계별 인지 구조, 공개 벤치마크 데이터셋 활용법 및 재현 가능한 멀티오믹스 컴퓨팅 환경을 구축할 수 있다.

---

## 1. 전립선암 연구에서 왜 10x Multiome (snRNA + snATAC)인가?

단일 핵 RNA 시퀀싱(snRNA-seq)은 종양 미세환경 내 이질적인 세포 계통과 악성 클론을 식별하는 강력한 도구입니다. 그러나 **전사체(Transcriptome)는 세포가 "현재 실행 중인 결과(Consequence)"만을 반영**할 뿐, **"왜 그 유전자가 전사되었는지" 또는 "앞으로 어떤 세포 상태로 전이할 잠재력을 지녔는지"를 설명하는 후성유전학적 메커니즘(Mechanism & Potential)**을 직접 보여주지 못합니다.

```
[후성유전학적 인과 사슬 (Epigenetic Causality)]
크로마틴 개방 (Chromatin Accessibility) ──► 전사인자 결합 (TF Binding) ──► 전사 개시 (Transcription) ──► 단백질 번역
         ▲                                                ▲
    snATAC-seq (조절 기전: 원인)                     snRNA-seq (발현 상태: 결과)
```

### 1.1 RNA 단독 분석의 생물학적 한계
1. **잠재(Primed) 상태 vs 침묵(Repressed) 상태의 식별 불가**:
   - 어떤 유전자의 mRNA 카운트가 $0$일 때, snRNA-seq 데이터로는 해당 프로모터/인핸서가 이질염색질(Heterochromatin)로 완전히 닫혀 있는 **불가역적 억제 상태(Repressed)**인지, 크로마틴은 이미 활짝 열려 있으나 적절한 전사인자 신호만을 기다리는 **기능적 준비 상태(Primed)**인지 구분할 수 없습니다.
2. **원격 시스-조절 인핸서(Distal Enhancers) 활성의 추적 불가**:
   - 전립선암의 핵심 종양 유전자들은 수십~수백 kb 떨어진 원격 인핸서 복합체에 의해 구동됩니다. mRNA 발현량만으로는 어떤 비코딩(Non-coding) 조절 좌위가 실제 해당 유전자의 전사를 구동하는지 추론할 수 없습니다.
3. **선구자 전사인자(Pioneer Transcription Factors)의 후성유전적 각인**:
   - *FOXA1*, *HOXB13*과 같은 선구자 인자는 응축된 염색질에 먼저 결합하여 크로마틴을 열어젖히며, 이후 안드로겐 수용체(*AR*)가 접근할 수 있는 통로를 개방합니다. 이러한 전사 네트워크의 선후 관계는 염색질 개방성 데이터 없이는 해석이 불가능합니다.

### 1.2 전산적 가상 매칭(In-silico Pairing)의 치명적 위양성
초기 단일세포 연구에서는 scRNA-seq과 scATAC-seq을 별도의 세포 샘플에서 독립적으로 수행한 뒤, 유전자 활성 점수(Gene Activity Score)나 공통 잠재 공간 투영(Co-embedding)을 통해 사후 전산 매칭(In-silico pairing)을 시도했습니다.

그러나 이는 극심한 종양 이질성을 지닌 전립선암에서 다음과 같은 치명적 오류를 낳습니다:
- **잘못된 세포 쌍 매칭(False-positive Pairing)**: 이종 세포 집단(예: *AR* 발현이 유사하지만 크로마틴 상태가 완전히 다른 두 클론)이 알고리즘에 의해 강제로 동일 세포로 묶입니다.
- **다대일/일대다 왜곡**: 크로마틴 변화는 유전자 발현 변화에 선행하는 경향이 보고되어 있으므로(시간 척도는 맥락 의존적), 발현 유사도만으로 세포 상태를 매핑하면 분화 궤적의 인과 관계가 역전됩니다.

**10x Multiome 기술**은 **동일한 물리적 단일 핵(The very same physical nucleus)**에서 snRNA와 snATAC 라이브러리를 물리적으로 1:1 바코딩하여 동시 포획하므로, 동일 GEM 내 동시 포획이므로 사후 전산 매칭에 따른 오류가 없는 결합 오믹스(Paired Multiomics)를 제공합니다. 다만 다중 핵 GEM(doublet), 앰비언트 RNA/DNA 오염, 바코드 변환 불일치 등은 여전히 QC가 필요하므로 '완전한 ground truth'는 아닙니다.

---

## 2. Track 01 학습 로드맵 및 인지 단계 (Cognitive Levels)

```
[Level 1: Biophysical Foundation]
  ├── Ch 01. 후성유전학 원리와 Tn5 트랜스포세이스 절단-접합 화학 (Tagmentation)
  └── Ch 02. 업스트림 정렬(Cell Ranger ARC, Chromap) 및 Fragment 파일 구조
          │
[Level 2: Data Cleaning & Low-Dimensional Geometry]
  ├── Ch 03. snATAC 특화 QC 지표 (TSS Enrichment Score, 147bp 뉴클레오솜 래더)
  └── Ch 04. 슈도벌크 MACS3 피크 콜링 및 LSI (TF-IDF + SVD) 차원 축소 기하학
          │
[Level 3: Dual-Modal Integration & Biological Discovery]
  ├── Ch 05. WNN (Weighted Nearest Neighbor) 멀티모달 통합 및 클러스터링
  └── Ch 06. Peak-to-Gene 연관 분석 (Chromatin-Gene Linkage & Distance Filtering)
          │
[Level 4: Expert Mastery: Gene Regulatory Networks]
  ├── Ch 07. 전사조절망 추론 (pycisTopic 모티프 디스커버리 & SCENIC+ eGRN 구축)
  └── Ch 08. 전립선암 치료 저항성 및 계통 가소성 후성유전체 마스터클래스
```

---

## 3. 선수 지식 점검 (Prerequisites)

- **후성유전학 기초**:
  - 진정염색질(Euchromatin) vs 이질염색질(Heterochromatin).
  - 뉴클레오솜(Nucleosome) 구조: 코어 히스톤 8량체 ($H2A, H2B, H3, H4 \times 2$)를 147 bp DNA가 1.65바퀴 감고 있는 형태.
  - 주요 히스톤 번역 후 수식화(PTM): 활성 마커($H3K27ac$, $H3K4me3$) vs 억제 마커($H3K27me3$, $H3K9me3$).
- **전립선암 전사 조절 기구**:
  - 마스터 선구자 인자: *FOXA1*, *HOXB13*.
  - 안드로겐 수용체 축: *AR*, *NKX3-1*, 상류 거대 인핸서 증폭 좌위.
  - 신경내분비 계통 전사인자: *ASCL1*, *NEUROD1*, *INSM1*, *REST*(억제자).
- **컴퓨팅 및 수학적 기초**:
  - 고차원 희소 행렬(Sparse Matrix) 연산 (`scipy.sparse.csr_matrix`).
  - 자연어 처리 및 정보 검색 모델: TF-IDF (Term Frequency - Inverse Document Frequency).
  - 선형대수: Truncated Singular Value Decomposition (SVD), 직교 투영.
  - 멀티모달 데이터 컨테이너: Python `muon` 라이브러리의 `MuData` 데이터 구조.

---

## 4. 재현 가능한 멀티오믹스 분석 환경 구축 (Environment Setup)

10x Multiome 분석은 대용량 지놈 정렬, Fragment 파싱, MACS3 피크 콜링 및 딥러닝/그래프 기반 멀티모달 통합 도구를 광범위하게 요구합니다. 격리된 Conda 가상환경을 구축합니다:

```bash
# 1. Multiome 전용 가상환경 생성 (Python 3.10 권장)
conda create -n multiome_env python=3.10 -y
conda activate multiome_env

# 2. 필수 바이오인포매틱스 C/C++ 기반 업스트림 유틸리티 설치
conda install -c bioconda samtools tabix bedtools ucsc-bedgraphtobigwig chromap -y

# 3. MACS3 피크 콜러 설치 (Python 3.10 호환)
pip install macs3

# 4. 단일세포 및 멀티오믹스 핵심 프레임워크 설치
pip install scanpy anndata muon episcanpy snapatac2

# 5. 후성유전학 모티프 분석 및 확장 전사조절망 도구 (SCENIC+)
pip install pycistopic pybiomart pygenometracks
# 08장 코드는 SCENIC+ v1.0.0 API 기준이므로 태그를 고정 (태그 존재는 GitHub에서 확인했으나, 이 설치 명령과
# 의존성 해결은 직접 실행해 검증하지 않았음; git main은 snakemake 중심으로 함수 시그니처가 다름)
pip install git+https://github.com/aertslab/scenicplus.git@v1.0.0

# 6. Jupyter 커널 등록
python -m ipykernel install --user --name multiome_env --display-name "Python (Multiome snATAC)"
```

### 설치 무결성 검증 스크립트
```python
import muon as mu
import scanpy as sc
import snapatac2 as snap
import pycisTopic

print(f"Muon Version: {mu.__version__}")
print(f"Scanpy Version: {sc.__version__}")
print(f"SnapATAC2 Version: {snap.__version__}")
print("All core multiome dependencies successfully verified.")
```

---

## 5. 실전 벤치마크 공개 멀티오믹스 데이터셋 카탈로그 (Public Multiome Datasets)

전립선암 10x Multiome 연구 및 방법론 검증을 위한 핵심 공개 데이터셋 카탈로그입니다:

| 연구 (논문) | 접근 번호 (Accession) | 플랫폼 / 모달리티 | 샘플 구성 및 환자 코호트 | 핵심 후성유전학적 분석 포인트 |
|:---|:---|:---|:---|:---|
| **10x Genomics Official Dataset** | `10x Genomics Dataset Hub` | Chromium Multiome (Cell Ranger ARC 2.0.0 배포본) | PBMC granulocyte-sorted 10k (전립선 종양 아님) | 10x 공개 Multiome 예제; Cell Ranger ARC 파이프라인 및 WNN 검증용 (공개 전립선 종양 Multiome 10x 데이터는 확인하지 못함) |

> **데이터셋 출처 주의**: 이전에 이 표에 있던 `GSE175634`(심근 분화 scRNA-seq)와 `GSE210358`(마우스 GEMM·오가노이드와 인간 mCRPC 시료의 RNA-seq, Chan et al. 2022)는 전립선 10x Multiome이 아님이 확인되어 삭제했습니다. 공개 전립선 Multiome accession은 아직 확정하지 못했습니다(accession 미확정: 원 논문의 Data Availability 확인 필요).

### 10x 공개 Multiome 예제(PBMC granulocyte-sorted 10k) 다운로드 가이드

> 아래 URL은 HTTP HEAD로 존재를 확인한 PBMC 예제이다(fragments 파일 약 2.4 GB). 이전에 적혀 있던 `human_prostate_tumor_multiome` 경로는 실재하지 않아(HTTP 403) 삭제했다. 전립선 종양 데이터를 쓰려면 원 논문의 Data Availability에서 확인된 accession/URL만 사용한다.
```bash
# 작업 디렉토리 생성
mkdir -p data/10x_pbmc_multiome && cd data/10x_pbmc_multiome

# 1. 필터링된 멀티오믹스 특성 매트릭스 (H5 형식: GEX + ATAC 동시 수록)
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

# 2. snATAC 골드 스탠다드 Fragment 파일 및 인덱스 (필수)
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

# 3. 종합 QC 및 정렬 요약 보고서 (Web Summary)
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_web_summary.html
```

---

## 6. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
단일 핵 전사체(snRNA-seq) 단독 분석 대비 10x Multiome (snRNA + snATAC) 기술이 종양 생물학 연구자에게 제공하는 3가지 결정적 분자생물학적 정보를 기술하고, 서로 다른 환자 검체에서 개별적으로 생성된 scRNA와 scATAC 데이터를 알고리즘으로 사후 매칭(In-silico pairing)할 때 발생하는 통계적 위험성을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **10x Multiome이 제공하는 3대 결정적 정보**:
   - **유전자 발현 잠재력(Primed vs Repressed)의 규명**: 현재 mRNA 발현이 0인 휴면 유전자에 대해 크로마틴이 열려 있어 즉각 반응 가능한 상태인지, 헤테로크로마틴으로 단단히 잠겨 영구 억제된 상태인지 명확히 구별할 수 있습니다.
   - **원격 비코딩 시스-조절 인핸서(Distal Cis-Regulatory Elements)와 표적 유전자의 1:1 결합**: 수십~수백 kb 떨어진 인핸서 피크의 개방성과 동일 물리적 핵 내 mRNA 발현 간의 상관계수(Peak-to-Gene Linkage)를 직접 연산하여 표적 암 유전자를 정확히 도출합니다.
   - **선구자 전사인자(Pioneer TFs)의 염색질 개방과 하류 유전자 전사의 인과관계 규명**: *FOXA1*, *HOXB13* 등의 결합 모티프 개방성과 하류 종양 신호 경로 활성화를 동일 세포 단위에서 직접 연계합니다.
2. **사후 전산 매칭(In-silico Pairing)의 통계적 위험성**:
   - 전립선암은 환자 간, 클론 간 유전적/후성유전적 이질성이 극심합니다. 사후 전산 매칭은 "발현 패턴이 비슷하면 염색질 구조도 비슷할 것"이라는 가정을 전제하지만, 실제 생물학에서는 염색질 개방이 전사 개시보다 수 시간~수일 선행하므로 중간 과도기 세포(Plastic/intermediate cells)가 엉뚱한 성숙 세포로 강제 매핑되는 **위양성 매칭 오류(False-positive Pairing Artifact)**를 초래합니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
Python 기반 멀티오믹스 분석 프레임워크인 `MuData` (Muon 라이브러리)의 데이터 저장 구조를 설명하고, 전통적인 단일 모달리티 `AnnData` 객체와 비교하여 다중 모달리티 간 바코드 공유 및 슬라이싱이 어떻게 작동하는지 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **`MuData` 아키텍처**:
   - `MuData` (Multimodal Data)는 서로 다른 생물학적 모달리티를 담당하는 개별 `AnnData` 객체들을 하위 모달리티 슬롯(`mdata.mod['rna']`, `mdata.mod['atac']`)으로 계층적 포함하는 상위 컨테이너 구조입니다.
2. **`AnnData`와의 비교 및 바코드 슬라이싱**:
   - 기존 `AnnData`는 단일 특성 행렬(예: 세포 × 유전자)만을 관리하지만, `MuData`는 공유된 세포 바코드(`obs`)를 축으로 하여 `rna`의 유전자 카운트 행렬(세포 × 36,000 유전자)과 `atac`의 게놈 피크 행렬(세포 × 150,000 피크)을 동시에 인덱싱합니다.
   - `mdata.obs`는 두 모달리티의 공통 메타데이터(예: 환자 ID, 세포 유형, 이중핵 판정)를 관리하며, `mdata[selected_cells, :]`와 같이 슬라이싱을 수행하면 두 모달리티(`rna`, `atac`)의 내부 행렬이 완벽한 정렬 상태를 유지하며 동시에 필터링됩니다.

</details>

---

### [문제 3. 임상 사례 분석 (Clinical Interpretation Case Study)]
**[임상 시나리오]**
원발암 수술 후 안드로겐 차단 요법(ADT: Leuprolide + Enzalutamide) 치료를 받던 중 생화학적 재발(PSA 상승)과 함께 간 전이가 발생한 mCRPC 환자의 간 생검 조직을 채취하여 10x Multiome 시퀀싱을 수행했습니다.
분석 결과, 간 전이 악성 암세포 클러스터에서 *AR*(안드로겐 수용체)과 *KLK3*(PSA)의 mRNA 발현량은 거의 $0$에 가까웠으나, 염색질 분석 결과 *AR* 유전자 상류 650 kb에 위치한 거대 증폭 인핸서(Centromeric Super-Enhancer) 영역과 *AR* 프로모터의 크로마틴 개방성(ATAC Peak)은 정상 전립선 상피세포 대비 8배 이상 강력하게 열려 있었습니다.

**[질문]**
1. 이 암세포 클러스터의 상태를 mCRPC 4대 분자 표현형 모델 관점에서 진단하고, mRNA 발현과 염색질 개방성 간의 불일치(Discrepancy)가 발생하는 분자생물학적 메커니즘을 설명하시오.
2. 만약 연구팀이 snRNA-seq 단독 분석만을 수행했을 경우 어떤 잘못된 결론을 내렸을지 비판하고, 10x Multiome 데이터를 토대로 제안할 수 있는 차세대 병용 표적 치료 전략을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **분자 진단 및 메커니즘 분석**:
   - **진단 (가설)**: 이 세포군을 '후성유전학적으로 AR 신호가 재활성화될 준비가 된 잠재 상태(Epigenetically Primed AR State)'로 볼 수 있다는 것은 하나의 *가설*이며, t-NEPC/DNPC와 구별하려면 신경내분비 마커 발현·병리·*RB1*/*TP53* 상태 같은 직교 근거가 필요합니다.
   - **불일치 메커니즘에 대한 주의**: 항안드로겐(엔잘루타마이드)은 AR의 활성(핵 내 전위, DNA 결합)을 억제할 뿐 *AR* mRNA를 0으로 만들지 않으며, 항안드로겐 저항 모델에서는 오히려 *AR* mRNA가 증가한다고 보고되었습니다(Chen et al., Nat Med 2004, PMID 14702632). 따라서 '약물이 mRNA 검출을 막았다'는 설명은 성립하지 않습니다. 또한 ATAC 신호는 DNA 카피수에 비례하므로, mCRPC에서 흔한 *AR* 상류 인핸서 중복(70–87%, Viswanathan et al., Cell 2018, PMID 29909985)이 있다면 '8배 개방'은 카피수 효과일 수 있어 카피수 보정(예: CNV 추정) 없이는 'primed' 상태로 해석할 수 없습니다. 이 시나리오는 가상이며, PSA 상승의 공급원은 이 클러스터가 아닌 다른 AR 활성 세포일 가능성도 있습니다.
2. **단독 분석의 오류 및 치료 전략**:
   - **단독 분석 시의 위험**: snRNA-seq만 보았다면 *AR* 음성이므로 "AR 경로가 결손된 AR-비의존성 종양"으로 해석할 수 있으나, 위 이유로 ATAC 단독 해석도 카피수 교란에 취약하므로 두 모달리티와 카피수·NE 마커를 함께 확인해야 합니다.
   - **치료 전략 (가설적 탐색)**: 열린 인핸서가 실제로 AR 재활성화 경로임이 확인된다면, 전사 조절 복합체를 표적하는 **BET 브로모도메인 억제제** 또는 CDK7/CDK9 억제제를 AR 길항제와 병용하는 접근이 *연구 가설*로 제안될 수 있으나, 이는 확립된 표준 치료가 아니며 이 시나리오로부터 임상 권고를 도출할 수 없습니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Chan JM, Zaidi S, Love JR, et al. Lineage plasticity in prostate cancer depends on JAK/STAT inflammatory signaling. *Science*. 2022;377(6611):1180-1191. PMID: 35981096. doi: 10.1126/science.abn0478
2. Chen CD, Welsbie DS, Tran C, et al. Molecular determinants of resistance to antiandrogen therapy. *Nat Med*. 2004;10(1):33-9. PMID: 14702632. doi: 10.1038/nm972
3. Viswanathan SR, Ha G, Hoff AM, et al. Structural Alterations Driving Castration-Resistant Prostate Cancer Revealed by Linked-Read Genome Sequencing. *Cell*. 2018;174(2):433-447.e19. PMID: 29909985. doi: 10.1016/j.cell.2018.05.036
### 검수 근거 문헌
1. Ma S, Zhang B, LaFave LM, et al. Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin. *Cell*. 2020;183(4):1103-1116.e20. PMID: 33098772. doi: 10.1016/j.cell.2020.09.056
2. Elorbany R, Popp JM, Rhodes K, et al. Single-cell sequencing reveals lineage-specific dynamic genetic regulation of gene expression during human cardiomyocyte differentiation. *PLoS Genet*. 2022;18(1):e1009666. PMID: 35061661. doi: 10.1371/journal.pgen.1009666
3. Tang F, Xu D, Wang S, et al. Chromatin profiles classify castration-resistant prostate cancers suggesting therapeutic targets. *Science*. 2022;376(6596):eabe1505. PMID: 35617398. doi: 10.1126/science.abe1505
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- muon 0.2.0 (github.com/scverse/muon). 원 논문: Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. *Genome Biol*. 2022;23(1):42. PMID: 35105358. doi: 10.1186/s13059-021-02577-8
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- SnapATAC2 2.10.0 (github.com/scverse/SnapATAC2). 원 논문: Zhang K, Zemke NR, Armand EJ, et al. A fast, scalable and versatile tool for analysis of single-cell omics data. *Nat Methods*. 2024;21(2):217-227. PMID: 38191932. doi: 10.1038/s41592-023-02139-9
- pycisTopic (github.com/aertslab/pycisTopic).
<!-- REFERENCES:END -->

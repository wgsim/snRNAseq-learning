# 00. 업스트림 정렬 및 BAM 분석 로드맵: 시퀀서 형광 신호부터 단일세포 전사체 카운트 행렬까지 (Upstream Alignment & BAM Roadmap)

---

## 학습 목표 (Cognitive Goals)

1. **차세대 시퀀싱 업스트림 블랙박스(Upstream Blackbox)의 완전한 탈신비화**: 일루미나(Illumina) 시퀀서의 4색/2색 레이저 형광 영상 신호(BCL)로부터 베이스 콜링(Base Calling), FASTQ 생성, 서열 정렬(Alignment/Pseudoalignment), UMI 중복 제거(Deduplication), 최종 유전자 발현 카운트 행렬($N \times G$) 생성에 이르는 계산생물학적 데이터 변환 파이프라인의 전 과정을 전산학적 계층 구조로 체계화한다.
2. **전립선암(Prostate Cancer) 정밀의학에서의 BAM 레벨 분석 필연성 규명**: 단순 유전자 수준(Gene-level) 카운트 행렬이 놓치는 대체 스플라이싱(Alternative Splicing, 예: 치료 저항성 안드로겐 수용체 변이체 *AR-V7*) 및 염색체 구조 변이(Structural Rearrangement, 예: *TMPRSS2-ERG* 융합)를 검출하기 위해 왜 원시 BAM 파일의 접합부 리드(Junction Reads) 및 CIGAR 스트링 분석이 필수적인지 임상 분자생물학적 기전을 이해한다.
3. **표준 벤치마크 데이터셋 구조 분석**: 10x Genomics Chromium 3' v3 PBMC 표준 데이터와 거세저항성 전립선암(CRPC) 조직 단일핵 RNA 시퀀싱(snRNA-seq) 원시 FASTQ의 리드 구조, paired-end 규격, 품질 관리(QC) 메트릭을 정량적으로 비교 분석한다.
4. **프로덕션급 생물정보학 업스트림 연산 환경 구축**: `STAR` (v2.7.11+), `salmon`/`alevin-fry`, `samtools`, `pysam`, `umi_tools` 등 최신 C/C++ 및 Rust 기반 정렬·정량 도구 스택을 설정하고, 재현 가능한 가상환경 명세와 시스템 성능 검증 파이프라인을 구축한다.

---

## 1. 업스트림 블랙박스의 탈신비화: BCL에서 카운트 행렬까지

단일세포 전사체 분석의 대다수 실무자는 Cell Ranger나 STARsolo가 출력하는 `filtered_feature_bc_matrix` (또는 `matrix.mtx`)에서 분석을 시작한다. 그러나 이 카운트 행렬은 시퀀서 하드웨어와 전산 알고리즘이 협력하여 처리한 수억 개의 연산 결과물일 뿐이다. 원시 물리 신호가 생물학적 수치로 변환되는 과정을 단계별로 분해한다.

```
[Illumina Sequencer Flowcell]
          │
          ▼ 1. Laser Excitation & Fluorescence Imaging (λ1, λ2, λ3, λ4)
[Binary Base Call (BCL) Files]
          │
          ▼ 2. bcl2fastq / bcl-convert (Demultiplexing by i7/i5 Sample Indices)
[FASTQ Files (R1: CB+UMI, R2: cDNA Insert)]
          │
          ├──────────────────────────────────────────────┐
          ▼ 3A. Spliced Genome Alignment                 ▼ 3B. Quasi-mapping / Pseudoalignment
      [STARsolo]                                     [alevin-fry]
          │                                              │
          ▼ Coordinate-Sorted BAM                        ▼ RAD (Reduced Alignment Data)
          │ (CIGAR, MD, NH, CB, UB tags)                 │ (Chunk-based binary stream)
          │                                              │
          └──────────────────────┬───────────────────────┘
                                 │
                                 ▼ 4. Error Correction & Graph Deduplication
                                    - 1-Hamming Distance CB correction against Whitelist
                                    - Directional Adjacency Graph UMI deduplication
                                 │
                                 ▼ 5. Sparse Matrix Assembly
[Count Matrix: raw_feature_bc_matrix (MEX / H5AD)]
          │
          ▼ 6. Cell Calling (EmptyDrops Dirichlet-Multinomial vs. Knee-Point)
[Filtered Count Matrix (N_cells × G_genes)]
```

### 1.1 시퀀싱 하드웨어의 물리적 신호 변환 (Physics to Bits)

현대 Illumina 시퀀서(NovaSeq 6000, NovaSeq X 등)는 플로우셀(Flowcell) 상의 나노웰(Nanowell)에 고정된 DNA 클러스터(Cluster)를 대상으로 합성 기반 시퀀싱(Sequencing-by-Synthesis, SBS)을 수행한다:

1. **형광 표지 뉴클레오타이드 혼입**: 각 사이클마다 형광단(Fluorophore)과 가역적 종결자(Reversible Terminator)가 결합된 dNTP(A, C, G, T)가 도입된다.
2. **레이저 여기 및 광학 센서 검출**: 레이저 광선에 의해 방출된 형광 파장을 고해상도 CCD/CMOS 센서가 픽셀 강도(Intensity) 데이터로 기록한다.
3. **BCL(Binary Base Call) 포맷 저장**: 각 웰 위치의 형광 강도 비율을 수학적으로 디컨볼루션(Deconvolution)하여 식별된 염기(A, C, G, T, N)와 이에 대응하는 품질 신뢰도(Phred Score)를 8비트 바이너리 구조로 압축하여 사이클별 디렉토리에 저장한다.

### 1.2 BCL에서 FASTQ로의 변환 (Demultiplexing)

`bcl2fastq` 또는 `bcl-convert` 소프트웨어는 실험 설계서(SampleSheet.csv)를 참조하여:
* 플로우셀의 레인(Lane)별 형광 데이터를 디멀티플렉싱(Demultiplexing)한다.
* 시퀀싱 라이브러리 양 말단에 부여된 인덱스(i7, i5 Index) 서열의 1-2bp 미스매치를 허용하면서 개별 샘플별 FASTQ 파일 쌍(Read 1, Read 2)을 생성한다.

### 1.3 데이터 계층별 파일 포맷과 압축 아키텍처

단일세포 분석 파이프라인에서 다루는 주요 파일 포맷의 구조적 특징과 메모리 오버헤드는 다음과 같다:

| 단계 | 파일 형식 | 주요 보관 정보 | 내부 압축 알고리즘 | 대략적 용량 (1억 리드 기준) |
| :--- | :--- | :--- | :--- | :--- |
| **원시 데이터** | `BCL / CBCL` | 사이클별 클러스터 형광 강도 및 염기 콜 | Snappy / GZIP 바이너리 | ~15 - 25 GB |
| **시퀀싱 리드** | `FASTQ.gz` | 리드 식별자, 16bp CB + 12bp UMI, 90bp cDNA 서열, Phred Q-score | Block GZIP (BGZF) | ~5 - 8 GB |
| **유전체 정렬** | `BAM` | 염색체 좌표, CIGAR 스트링, MAPQ, 플래그, SAM 태그 (CB, UB, GN) | BGZF 압축 이진 SAM | ~8 - 15 GB |
| **의사정렬 스트림** | `RAD` | 세포 바코드, UMI, 전사체 동등성 클래스 (Equivalence Class) | 청크 단위 이진 직렬화 | ~0.8 - 1.5 GB |
| **최종 발현 행렬** | `MEX / MTX` | 세포 식별자(`barcodes.tsv`), 유전자 정보(`features.tsv`), 희소 카운트(`matrix.mtx`) | 텍스트 압축 (GZIP) | ~50 - 200 MB |
| **단일세포 객체** | `H5AD` | 발현 카운트, 관측치(obs), 변수(var), 저차원 임베딩 | HDF5 계층형 압축 (Chunked/Gzip) | ~30 - 150 MB |

---

## 2. 전립선암(Prostate Cancer) 정밀의학에서의 임상적 필연성

단일세포 RNA 시퀀싱에서 유전자 발현 수준(Gene-level quantification)만 관찰하는 것은 심각한 정보 손실을 초래한다. 특히 거세저항성 전립선암(Castration-Resistant Prostate Cancer, CRPC)과 같은 난치성 암종에서는 동일 유전자좌(Locus)에서 생성되는 다양한 스플라이싱 변이체와 유전체 구조 변이가 환자의 생존과 치료 반응을 직접 결정한다.

```
[Normal Full-length Androgen Receptor (AR-FL)]
Exon:  [ 1 ]───[ 2 ]───[ 3 ]───[ 4 ]───[ 5 ]───[ 6 ]───[ 7 ]───[ 8 ]
Domain: NTD     DBD     Hinge       Ligand-Binding Domain (LBD)
                                    ▲ Enzalutamide / Abiraterone 결합 부위

[Constitutively Active Splicing Variant (AR-V7)]
Exon:  [ 1 ]───[ 2 ]───[ 3 ]───[ CE3 ] (Cryptic Exon 3, CE3: 인트론 3 안의 대체 3′ 말단 엑손)
Domain: NTD     DBD     (LBD 결손!)
                        ▲ 리간드 없이도 전사인자 상시 활성화 -> 내성 유발!
```

### 2.1 안드로겐 수용체 스플라이스 변이체 (*AR-V7*)와 치료제 내성 기전

*AR* (Androgen Receptor) 유전자는 8개의 표준 엑손으로 구성되며, 전립선암 표준 표적 치료제 중 엔잘루타마이드(Enzalutamide)는 C-말단의 리간드 결합 도메인(Ligand-Binding Domain, LBD, Exon 4-8 코딩)에 결합하는 AR 길항제이고, 아비라테론(Abiraterone)은 CYP17A1 억제제로 안드로겐 합성 자체를 줄여 LBD를 통한 신호 전달을 차단한다(둘 다 LBD 의존적 신호를 표적으로 하므로 LBD가 없는 AR-V7에는 효과가 제한적일 수 있다).

그러나 진행성 CRPC 종양 세포는 인트론 3 내부에 존재하는 잠재적 엑손(cryptic exon 3, CE3; 문헌 표기는 '3b'가 아님)을 선택적으로 포함시키는 비정상적 대체 스플라이싱을 수행하여 **AR-V7 (Androgen Receptor Splice Variant 7)** 단백질을 생성한다:
* **구조적 특성**: Exon 1-2-3 뒤에 인트론 3 안의 대체 3′ 말단 엑손 CE3(RefSeq 변이체 3의 말단 엑손 약 1.4 kb, 3′UTR 포함)가 접합되고 16개의 변이체 특이 아미노산 뒤에서 번역이 종결되어 LBD(Exon 4-8) 코딩 서열이 전사체에 포함되지 않는다(엑손 4-8 건너뛰기나 일반 PTC가 아니라 말단 엑손 교체; 이전의 '101bp'는 근거 미확인으로 삭제).
* **임상적 의미**: 안드로겐 리간드(DHT)나 약제가 없어도 핵 내부로 자가 전위(Nuclear Translocation)하여 표적 유전자(*KLK3/PSA*, *TMPRSS2*, *UBE2C*) 전사를 촉진함으로써 치료제에 즉각적인 불응성(Intrinsic Resistance)을 나타낸다.
* **정량 실패 원인**: 일반적인 `CellRanger count`나 STARsolo의 표준 Gene 카운트는 *AR* 유전자 엔트리(`ENSG00000169083`) 아래로 모든 엑손의 리드를 단일 정수로 병합한다. 따라서 종양 세포가 정상 *AR-FL*을 발현하는지 내성 변이체인 *AR-V7*을 발현하는지 전혀 식별할 수 없다. 이를 규명하기 위해서는 Exon 3과 CE3를 가로지르는 **접합부 리드(Junction Reads)**를 BAM 파일에서 직접 정량해야 한다.

### 2.2 *TMPRSS2-ERG* 융합 유전자와 유전체 재배열 검출

서양 전립선암 환자의 약 50%에서 21번 염색체의 약 3Mb 미세 결실(Interstitial Deletion) 또는 염색체 전좌(Translocation)로 인해 안드로겐 조절 유전자인 *TMPRSS2*의 프로모터/엑손 1과 발암성 전사인자 *ERG* (ETS-related gene) 코딩 영역이 융합된다.

$$\text{Chromosomal Locus 21q22.2-3: } [\text{TMPRSS2 (Promoter/Exon 1)}] \xrightarrow{\Delta \approx 3\text{Mb}} [\text{ERG (Exon 4-16)}]$$

* **분자 병리학적 결과**: 안드로겐 자극에 의해 *ERG* 종양 유전자가 강력하게 이소성 발현(Aberrant Overexpression)되어 침습 및 상피-간엽 전이(EMT)가 유도된다.
* **BAM 수준 검출**: *TMPRSS2* Exon 1과 *ERG* Exon 4 사이의 비정상적 키메라 접합부(Chimeric Junction)에 걸쳐 정렬된 Split-read 및 불일치 정합(Discordant Paired-end) 리드를 검출하기 위해서는 정렬 매개변수와 정렬 태그(`CB`, `UB`, `SA`, `chimeric.out.sam`)를 정밀하게 추적해야 한다.

---

## 3. 벤치마크 데이터셋 및 연구 설계

본 트랙 전반에 걸쳐 사용할 핵심 벤치마크 데이터셋은 두 가지 상호보완적 축으로 구성된다:

```
[Benchmark Datasets for Track 09]
   │
   ├─► 1. 10x Genomics 3' v3 PBMC (Peripheral Blood Mononuclear Cells)
   │     - 세포 수 / 시퀀싱 깊이: 데이터셋마다 다름 (예시 규모: 약 10,000 세포, 세포당 약 50,000 reads)
   │     - 목적: 표준 면역세포 프로파일링, 바코드 보정, 고속 의사정렬(Alevin-fry) 벤치마크
   │     - 출처: 10x Genomics 공개 데이터셋 페이지(사용 시 정확한 데이터셋 이름·버전·URL 기재; 규모는 데이터셋별로 상이)
   │
   └─► 2. Human Primary & Metastatic Prostate Cancer snRNA-seq (Tumor Core)
         - 임상 표본: Primary Gleason Score 4+5 및 CRPC 골 전이 조직
         - 라이브러리: 10x Genomics Single Cell Multiome / Chromium 3' v3 Nuclei
         - 목적: 인트론 보존 전구체 전사체(pre-mRNA), AR-V7 접합부, TMPRSS2-ERG 융합 식별
         - 특성: 핵 분리(Nuclear isolation)에 따른 높은 비성숙 RNA(Intron) 비율 (>40%)
```

### 3.1 10x Genomics Chromium 3' v2 vs. v3 라이브러리 사양 비교

| 사양 지표 | Chromium Single Cell 3' v2 | Chromium Single Cell 3' v3 / v3.1 | snRNA-seq 분석 시 고려사항 |
| :--- | :--- | :--- | :--- |
| **Read 1 구조** | 16bp CB + 10bp UMI (총 26bp) | 16bp CB + 12bp UMI (총 28bp) | v3의 UMI 공간이 16배 확장 ($4^{12} \approx 1.67 \times 10^7$) |
| **Read 2 구조** | 98bp (cDNA insert) | 90bp - 150bp (cDNA insert) | 스플라이스 접합부 스팬(Span) 확률 증가 |
| **바코드 화이트리스트** | 737,280 개 (`737K-august-2016.txt`) | 6,794,880 개 (`3M-february-2018.txt`) | 바코드 충돌 확률 대폭 감소 |
| **캡처 효율(Capture Rate)** | 키트·샘플에 따라 상이 (수치는 출처 미확인으로 삭제) | 키트·샘플에 따라 상이 | v3에서 세포당 검출 UMI/유전자가 증가한다는 일반 경향(정량은 데이터셋 의존) |
| **snRNA-seq 호환성** | 전장 cDNA 전사체 중심 | (v2→v3 화학 차이 효과가 아님) 핵 내 pre-mRNA 비중과 oligo-dT의 인트론 A-rich 구간 프라이밍으로 snRNA에서 인트론 리드가 많음 | 전구 전사체(pre-mRNA) 포괄적 매핑 필수 |

---

## 4. 고성능 생물정보학 업스트림 분석 환경 구축 가이드

단일세포 업스트림 연산은 대규모 지놈 인덱스 로딩(Human GRCh38 기준 32GB RAM 필요)과 초당 수십만 개의 FASTQ 레코드 처리를 요구하므로 최적화된 C/Rust 바이너리 구축이 필수적이다.

### 4.1 Conda/Mamba 환경 명세서 (`environment_upstream.yml`)

아래 환경 명세서는 `STAR`, `salmon`, `alevin-fry`, `samtools`, `pysam`, `umi_tools`의 최신 호환 버전을 격리된 환경으로 정의한다.

```yaml
name: snrna_upstream_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - pip
  # Alignment & Quantification Core
  - star=2.7.11b
  - salmon=1.10.3
  - alevin-fry=0.9.0
  - pyroe  # 03장 스크립트가 사용 (bioconda)
  # SAM/BAM Processing
  - samtools=1.19.2
  - htslib=1.19.1  # 참고: pysam 0.22.0은 자체 번들 htslib(1.18)을 사용하므로 이 핀과 무관
  - bedtools=2.31.1
  - igv-reports=1.10.1
  # Python High-Performance Libraries
  - pysam=0.22.0
  - umi_tools=1.1.5
  - numpy>=1.24.0
  - scipy>=1.10.0
  - pandas>=2.0.0
  - scanpy>=1.9.8
  - anndata>=0.10.0
  - numba>=0.58.0
  # CLI Utilities
  - pigz=2.8
  - parallel=20240122
```

### 4.2 환경 설치 및 검증 스크립트

```bash
#!/usr/bin/env bash
# Step 1: Mamba를 통한 초고속 가상환경 생성
mamba env create -f environment_upstream.yml

# Step 2: 환경 활성화
conda activate snrna_upstream_env

# Step 3: 핵심 바이너리 설치 및 링크 상태 검증
echo "=== 업스트림 핵심 바이너리 무결성 점검 ==="
STAR --version || { echo "STAR 설치 실패"; exit 1; }
salmon --version || { echo "Salmon 설치 실패"; exit 1; }
alevin-fry --version || { echo "alevin-fry 설치 실패"; exit 1; }
samtools --version | head -n 2 || { echo "samtools 설치 실패"; exit 1; }
umi_tools --version || { echo "umi_tools 설치 실패"; exit 1; }

echo "=== Python HTSlib(pysam) 및 C-API 연동 검증 ==="
python -c "
import pysam
import umi_tools
# alevin-fry는 Rust CLI이므로 Python 모듈이 없다 (위 alevin-fry --version으로 점검)
print(f'pysam version: {pysam.__version__}')
print(f'bundled samtools: {pysam.__samtools_version__}')
print(f'bundled htslib: {pysam.version.__htslib_version__}')
"
echo "환경 설정 및 바이너리 연동 완료."
```

---

## 5. 본 트랙(Track 09)의 커리큘럼 아키텍처

```
[Track 09: Upstream Alignment & BAM Processing Curriculum]
│
├── 00_ALIGNMENT_AND_BAM_ROADMAP.md (본 문서)
│
├── 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md
│     └── 일루미나 화학, Phred 점수, 1-해밍 거리 볼 탐색, EmptyDrops 수리 모델
│
├── 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md
│     └── 접미사 배열(SA), BWT, Maximal Mappable Prefix, 스플라이스 접합부 동적계획법
│
├── 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md
│     └── k-mer 준정렬, Splici 인덱스, 이진 RAD 포맷, 메모리/연산 10배 절감
│
├── 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md
│     └── PCR 증폭 바이어스, UMI 방향성 인접 그래프(Directional Adjacency Graph) 수학적 증명
│
├── 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md
│     └── SAM 비트 플래그, CIGAR 스트링, 10x 태그(CB/CR, UB/UR, GN, GX), HTSlib 바이너리
│
├── 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md
│     └── pysam 프로그래밍 접근, pileup 기반 AR 변이 UMI 카운팅
│
├── 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md
│     └── 전립선암 AR-V7 정량, Cryptic Exon 3 (CE3) 리드 추출, pysam 스플라이스 접합부 파싱
│
├── 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md
│     └── TMPRSS2-ERG 융합 분할 리드/불일치 쌍 검출과 템플릿 스위칭 아티팩트 필터
│
└── 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md
      └── 업스트림 장애 진단, 베스트 프랙티스, 캡스톤 루브릭
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 유전자 카운트 행렬 vs. 원시 BAM 정렬 파일 간 정보 손실과 *AR-V7* 변이체 발현 검출 불능 문제
한 분자유전학 연구팀이 엔잘루타마이드 치료를 받는 전립선암 환자의 종양 검체로부터 단일세포 전사체(10x Genomics 3' v3) 데이터를 수집하여 표준 Cell Ranger 파이프라인(`cellranger count`)을 구동하였다. 연구원은 출력된 `filtered_feature_bc_matrix.h5` 객체에서 *AR* 유전자(`ENSG00000169083`)의 발현량이 높은 세포 군집을 발견하고 환자가 치료 저항성 변이체인 *AR-V7*을 다량 발현하고 있다고 결론지었다.
1. 이 결론이 생물정보학적으로 왜 심각한 오류(False Premise)를 내포하고 있는지 Cell Ranger의 유전자 레벨 카운트 생성 알고리즘(`GN`/`GX` 태그 할당 방식) 관점에서 비판하시오.
2. *AR-V7* 변이체의 실제 발현 여부를 참 양성(True Positive)으로 입증하기 위해, 원시 BAM 파일(`possorted_genome_bam.bam`)에서 반드시 확인해야 하는 정렬 특징(Feature)과 스플라이스 접합부의 좌표 조건을 명시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 결론의 오류 분석 (Gene-level Count의 한계)**:
*   Cell Ranger의 표준 정량 파이프라인은 리드가 참조 전사체(Ensembl GTF) 주석 상의 어떤 엑손에 오버랩되는지를 검사하여 단일 `gene_id`(`GX`) 및 `gene_name`(`GN`) 태그를 할당한다.
*   *AR* 유전자의 정상 전장 전사체(*AR-FL*)와 치료 저항성 변이체(*AR-V7*)는 엑손 1, 2, 3을 완벽하게 공유한다.
*   따라서 엑손 1, 2, 3에 정렬된 수많은 리드는 *AR-FL*과 *AR-V7* 양쪽 모두에서 기인할 수 있음에도 불구하고, 파이프라인은 이를 단순히 *AR*(`ENSG00000169083`) 카운트로 합산한다. 환자가 약제 반응성인 *AR-FL*만 높게 발현하고 있더라도 *AR* 유전자 카운트는 극도로 높게 측정되므로, 카운트 행렬 수치만으로 *AR-V7* 매개 약제 내성을 판정하는 것은 명백한 위양성(False Positive) 추론이다.

**2. 원시 BAM 파일 기반 입증 방법**:
*   *AR-V7*에 특이적인 분자적 표지는 **Exon 3의 3' 말단과 잠재적 엑손인 Cryptic Exon 3 (CE3)의 5' 말단을 연결하는 스플라이스 접합부(Splice Junction)**이다.
*   따라서 원시 BAM 파일에서 다음 조건을 만족하는 리드를 필터링해야 한다:
    1.  **CIGAR 스트링에 'N' 연산자(스플라이스 인트론) 포함**: 리드가 인트론 3 영역을 건너뛰어 Exon 3과 CE3에 걸쳐 정렬된 Split-read여야 한다 (예: AR-V7 접합부라면 `45M8546N45M`; GRCh38에서 exon 3 끝과 CE3 시작 사이 인트론 길이).
    2.  **접합부 게놈 좌표 일치**: Exon 3 공여자 부위(Donor site, GRCh38 기준 exon 3 끝 `chrX:67,686,126`)와 Cryptic Exon 3(CE3) 수용자 부위(Acceptor site, CE3 시작 `chrX:67,694,673`; AR은 + 가닥이므로 donor < acceptor)를 정확히 연결하는 접합부 리드 수($N_{\text{junction}}$)를 확인하고, 이 리드들이 고유한 UMI(`UB` 태그)를 지니고 있는지 검증해야 한다.

</details>

---

### 문제 2: 10x Genomics 3' v3 라이브러리에서 Read 1 손상 시 복구 한계
단일세포 시퀀싱 중 NovaSeq 장비의 광학 오류로 인해 Read 1 파일(`R1.fastq.gz`)의 1번부터 16번째 사이클(Cell Barcode 영역) 전반에서 Phred Quality score가 Q5 이하로 급락하여 대규모 염기 미결정(Base N-calling)이 발생하였다. 반면 Read 2(cDNA insert)의 서열은 평균 Q38 이상의 무결성을 유지하였다.
1. 이 상황에서 고품질의 Read 2 cDNA 서열만을 이용하여 유전체 정렬(STAR)을 수행할 경우, 벌크(Bulk) RNA-seq 수준의 유전자 발현 정량은 가능한지 여부를 판정하시오.
2. 단일세포 전사체 분석 관점에서 세포별 분해능(Single-cell resolution) 복구가 기술적으로 불가능한 이유를 정보이론적(Information Theory) 엔트로피와 조합론 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 벌크 수준 발현 정량 가능 여부**:
*   **가능하다**. Read 2(cDNA insert)는 mRNA 분자의 전사체 서열 자체를 담고 있으므로, 유전체(GRCh38)에 정렬하여 어떤 유전자에서 전사된 리드인지 식별하고 조직 전체의 총 전사량을 집계(Bulk RNA-seq 모드)하는 것은 완벽히 수행할 수 있다.

**2. 단일세포 분해능 복구 불가능 이유 (정보이론 및 조합론)**:
*   단일세포 시퀀싱에서 리드의 물리적 기원(어느 단일 액적/세포에서 유래했는가)을 지정하는 유일한 메타데이터는 Read 1의 16bp Cell Barcode이다.
*   16bp 바코드가 가질 수 있는 이론적 상태 공간은 $4^{16} \approx 4.29 \times 10^9$개이며, 10x Genomics 3' v3 화이트리스트는 이 중 약 $6.79 \times 10^6$개의 사전 정의된 직교(Orthogonal) 좌표를 사용한다.
*   만약 16bp 중 N으로 손상된 염기 수 $k$가 커지면 후보가 $4^k$개로 늘어나며 각 후보가 화이트리스트에 있을 확률은 약 $6{,}794{,}880/4^{16}\approx1.6\times10^{-3}$이므로 기대 중복 후보 수는 $4^k\times1.6\times10^{-3}$입니다(예: $k=4$면 약 0.4개로 복구 가능, 수천 개 후보는 $k\approx11$ 이상에서 성립). 16 사이클 전체가 손상되면 복구가 불가능하다는 결론은 유효합니다.
*   Read 2의 서열은 세포 기원 정보를 전혀 담고 있지 않으므로(서로 다른 1,000개의 세포에서 유래한 동일 유전자 mRNA는 완벽히 동일한 Read 2 서열을 가짐), 손실된 바코드 정보를 Read 2로부터 역산하거나 복원할 수 있는 수학적 정보량(Mutual Information $I(R2; CB) = 0$)이 전무하다. 따라서 세포별 디멀티플렉싱은 회복 불가능하게 파괴된다.

</details>

---

### 문제 3: 고처리량 단일세포 업스트림 파이프라인의 디스크 I/O 병목과 공유 메모리(`/dev/shm`) 최적화
100개의 대규모 단일세포 FASTQ 샘플(총 용량 약 2TB)을 128코어 CPU와 512GB RAM, 고성능 NVMe SSD 스토리지를 갖춘 단일 계산 노드에서 병렬로 STARsolo를 구동하여 처리하고자 한다. 연구원은 16개의 샘플을 동시에 실행(`--runThreadN 8`로 16개 프로세스 동시 구동)하였으나, CPU 이용률이 30% 미만에 머물고 심각한 I/O 대기(I/O Wait) 현상이 발생함을 목격하였다.
1. 인간 유전체(GRCh38)의 접미사 배열(Suffix Array) 인덱스 크기를 고려할 때, 16개의 독립 프로세스가 동시에 인덱스를 개별 로드할 경우 발생하는 RAM 부족(OOM) 및 스왑(Swap) 디스크 스래싱(Thrashing) 메커니즘을 분석하시오.
2. STAR의 `--genomeLoad LoadAndExit`/`LoadAndKeep`/`Remove` 옵션과 Linux 공유 메모리(`/dev/shm`, POSIX Shared Memory) 아키텍처를 결합하여 이 디스크 I/O 병목을 해결하고 연산 처리량을 3배 이상 극대화하는 파이프라인 구성 전략을 수립하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. OOM 및 디스크 스래싱 메커니즘 분석**:
*   인간 게놈 GRCh38(약 3.1 Gb)의 STAR 접미사 배열(Suffix Array) 인덱스는 메모리 상에서 약 30~32 GB의 RAM을 점유한다.
*   16개의 독립된 STAR 프로세스가 기본 설정(`--genomeLoad NoSharedMemory`)으로 실행되면 각 프로세스가 독자적인 32GB 인덱스를 메모리에 할당하려 시도하므로, 요구 메모리는 $16 \times 32\text{ GB} = 512\text{ GB}$에 도달한다.
*   여기에 리드 정렬 버퍼, 동적 계획법 매트릭스, FASTQ 버퍼가 프로세스당 수 GB씩 추가로 필요하므로 시스템 물리 메모리(512GB)를 즉시 초과하여 Linux 커널의 OOM(Out-Of-Memory) Killer가 작동하거나 스왑 파티션(Swap Space)으로 메모리 페이지가 강제 이동된다.
*   NVMe SSD라 할지라도 수십 GB의 페이지 폴트(Major Page Fault)와 스왑 I/O를 감당하지 못하여 심각한 I/O 스래싱이 발생하고 CPU는 대기(I/O Wait) 상태에 빠지게 된다.

**2. 공유 메모리(`--genomeLoad`) 최적화 전략**:
*   **메모리 상주 인덱스 생성**: 파이프라인 시작 전 마스터 프로세스가 인덱스를 공유 메모리 세그먼트에 단 한 번만 적재하고 종료한다(`LoadAndExit`: 읽기 입력 없이 적재만 수행):
    ```bash
    STAR --genomeDir /path/to/GRCh38_index --genomeLoad LoadAndExit
    ```
*   **프로세스 간 인덱스 공유 (Zero-copy Shared Memory)**:
    16개의 STARsolo 프로세스가 구동될 때 모두 `--genomeLoad LoadAndKeep`을 부여한다(`LoadAndRemove`를 병렬 작업에 쓰면 먼저 끝난 작업이 인덱스를 제거해 다른 작업과 충돌할 수 있다; `SharedMemory`는 유효한 값이 아니다). 공유 메모리를 쓰는 경우 `--outSAMtype BAM SortedByCoordinate`에는 `--limitBAMsortRAM` 지정이 필수이다(기본값 0 불가).
    모든 워커 프로세스는 물리 메모리에 단 1벌(32GB)만 존재하는 인덱스를 `shmget`/`mmap` 시스템 콜을 통해 가상 메모리 공간에 읽기 전용(Read-only)으로 매핑하여 공유한다.
*   **효과**:
    *   메모리 사용량이 512GB에서 32GB + (워커당 작업 버퍼 ~3GB $\times 16) \approx 80\text{ GB}$ 수준으로 극적으로 감소한다.
    *   32GB 인덱스를 디스크에서 반복 읽기(Disk Read)할 필요가 전혀 없으므로 NVMe I/O 대기가 0에 수렴하며, 128코어 CPU가 100% 포화 상태로 연산에 집중할 수 있어 전체 분석 처리량이 수직 상승한다.
    *   마지막 16번째 작업 완료 시 인덱스를 공유 메모리에서 안전하게 언로드(`STAR --genomeLoad Remove`)하도록 파이프라인을 종료한다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 검수 근거 문헌
1. Hu R, Dunn TA, Wei S, et al. Ligand-independent androgen receptor variants derived from splicing of cryptic exons signify hormone-refractory prostate cancer. *Cancer Res*. 2009;69(1):16-22. PMID: 19117982. doi: 10.1158/0008-5472.CAN-08-2764
2. de Bono JS, Logothetis CJ, Molina A, et al. Abiraterone and increased survival in metastatic prostate cancer. *N Engl J Med*. 2011;364(21):1995-2005. PMID: 21612468. doi: 10.1056/NEJMoa1014618
3. Hong R, Koga Y, Bandyadka S, et al. Comprehensive generation, visualization, and reporting of quality control metrics for single-cell RNA sequencing data. *Nat Commun*. 2022;13(1):1688. PMID: 35354805. doi: 10.1038/s41467-022-29212-9
4. Curion F, Wu X, Heumos L, et al. hadge: a comprehensive pipeline for donor deconvolution in single-cell studies. *Genome Biol*. 2024;25(1):109. PMID: 38671451. doi: 10.1186/s13059-024-03249-z
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- alevin-fry (github.com/COMBINE-lab/alevin-fry). 원 논문: He D, Zakeri M, Sarkar H, et al. Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data. *Nat Methods*. 2022;19(3):316-322. PMID: 35277707. doi: 10.1038/s41592-022-01408-3
- SAMtools (github.com/samtools/samtools). 원 논문: Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *Gigascience*. 2021;10(2). PMID: 33590861. doi: 10.1093/gigascience/giab008
- STAR 2.7.11 (github.com/alexdobin/STAR). 원 논문: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15-21. PMID: 23104886. doi: 10.1093/bioinformatics/bts635
- Salmon (github.com/COMBINE-lab/salmon). 원 논문: Patro R, Duggal G, Love MI, et al. Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods*. 2017;14(4):417-419. PMID: 28263959. doi: 10.1038/nmeth.4197
- UMI-tools (github.com/CGATOxford/UMI-tools). 원 논문: Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Res*. 2017;27(3):491-499. PMID: 28100584. doi: 10.1101/gr.209601.116
<!-- REFERENCES:END -->

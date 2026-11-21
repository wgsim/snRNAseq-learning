# snRNA-seq 마스터 과정: 전립선암(원발암 & 전이암) 중심 분석 및 이론 체계

> **From Fundamentals to Expert**: 이론적 생물학/통계학 원리 및 분석 프로세스 인지(Cognitive Understanding) 중심의 단일 핵 전사체(Single-nucleus RNA-seq) 완성 가이드

> [!IMPORTANT]
> **면책 안내 (Disclaimer)**
> - 이 자료는 **교육용**입니다. 임상 진단, 치료 결정, 환자 관리에 사용하지 마십시오. 의학적 판단은 자격 있는 의료 전문가와 상의하십시오.
> - 파이프라인 결과(`results_pipeline_*`)는 모두 **합성(시뮬레이션) 데이터**로 만든 것이며 실제 환자나 실험의 측정값이 아닙니다(각 결과 폴더의 `SIMULATION_NOTICE.md` 참고).
> - 수치, 임계값, 벤치마크, 도구 사용법은 작성 시점의 문헌과 소프트웨어 버전을 기준으로 하며 틀리거나 오래되었을 수 있습니다. 실제 분석에서는 원 논문과 각 도구의 공식 문서를 확인하십시오. 검수 기록은 `review/REVIEW_REPORT.md`에 있습니다.
> - 자료 작성과 검수에 AI 도구를 활용했습니다. 검수는 문헌·공식 문서·코드 실행으로 확인했지만 모든 문장을 원문과 대조하지는 않았습니다.
>
> 라이선스: [MIT](LICENSE) © 2026 wgsim

---

## 🌟 단일 세포 & 공간 오믹스 마스터 커리큘럼 (Master Ecosystem)

> **전체 11개 트랙 종합 로드맵 (100% 완주)**: [MASTER_CURRICULUM_ROADMAP.md](MASTER_CURRICULUM_ROADMAP.md)에서 단일 세포 및 공간 오믹스 11대 마스터 트랙(Track 00 ~ Track 10)의 전체 구조, 학습자별 권장 경로, 완성된 교재 및 산출물 링크를 확인하실 수 있습니다.

### 📚 개설 및 완주된 마스터 트랙
- **[Track 00: snRNA-seq Foundations & Core Downstream Pipeline](README.md#2-디렉토리-구조-directory-structure)** (루트 디렉토리)
  - 전립선암 동결 조직 핵 분리, Pre-mRNA 인트론, QC & CellBender, UMAP/Leiden, inferCNV, scVI, pySCENIC
  - [마스터 교재 (Markdown)](snRNAseq_Prostate_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](snRNAseq_Prostate_Master_Course.html) | [Word 문서 (DOCX)](snRNAseq_Prostate_Master_Course.docx)
- **[Track 01: 10x Multiome (snRNA + snATAC) & Epigenomics](tracks/01_multiome_snatac/)**
  - 단일 핵 동시 전사체+후성유전체, Tn5 화학, Fragment 파일, LSI, WNN 가중 통합, Peak-to-Gene 링크, SCENIC+ 삼중 eRegulon
  - [멀티오믹스 교재 (Markdown)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.html) | [Word 문서 (DOCX)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.docx)
- **[Track 02: 공간 전사체학 (Spatial Transcriptomics)](tracks/02_spatial_transcriptomics/)**
  - 전립선암 조직 병리(Gleason/Cribriform) 아키텍처, 10x Visium/Visium HD vs Xenium/CosMx/MERSCOPE, SpatialData OME-Zarr, Cellpose vs Baysor, cell2location 베이지안 디콘볼루션, Moran's I/SpatialDE, Squidpy Niche 및 COMMOT 최적수송 L-R 신호망
  - [공간전사체 교재 (Markdown)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.html) | [Word 문서 (DOCX)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.docx)
- **[Track 03: 종양 미세환경(TME) 세포 간 신호 전달 (Cell-Cell Communication)](tracks/03_tme_cell_communication/)**
  - TME 세포 사회학, 3대 신호 양식(Juxta/Para/Auto), OmniPath/LIANA Borda 합의 집계, 발현 임계치 & 순열 검정, CellChat 다중 소단위체 기하평균 & 중심성, NicheNet RWR 조절 잠재력 표적인과 추론, LIANA+ 가우시안 거리 감쇄 커널, MISTy 다중 뷰 머신러닝, mCRPC 뼈 전이 악순환 5자 회로 및 병용 표적 치료 전략
  - [세포대화 교재 (Markdown)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.html) | [Word 문서 (DOCX)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.docx)
- **[Track 04: 단일 세포 암 유전체 & 클론 진화 (Single-Cell Cancer Genomics & Clonal Evolution)](tracks/04_cancer_genomics_clones/)**
  - 이수성(Aneuploidy) 감별, 전립선암 6대 핵심 CNA (8p, 8q, 10q, 13q, 17p, Xq), inferCNV 6상태 HMM 비터비 세그멘테이션, CopyKAT GMM 정상 이배체 기저선 탐색, Numbat 위상차 BAF 결합 베이지안 HMM 및 복제중립 LOH (CN-LOH), 단일세포 NJ/MEDICC2 계통수 복원, 체세포 SNV 및 TMPRSS2-ERG 융합, ADT 선택압 매개 t-NEPC 진화, 원발암-골전이 파종(Seeding) 모델링
  - [암유전체 교재 (Markdown)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.html) | [Word 문서 (DOCX)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.docx)
- **[Track 05: 단일 세포 생성형 딥러닝 & 파운데이션 모델 (Generative Foundation Models)](tracks/05_generative_foundation_models/)**
  - 원시 카운트 확률 모델링(Poisson vs NB vs ZINB), scVI 심층 변분 오토인코더(ELBO & 재매개변수화 트릭), scANVI 준지도 세포형 전이(희귀 t-NEPC 고정밀 주석), scArches 전역 재학습 없는 지속적 아키텍처 수술, 트랜스포머 유전자 토큰화 & 임베딩, scGPT 다중 작업 사전학습(MGM, GEP) 및 창발적 조절망(GRN), Geneformer 랭크-값 인코딩, 인 실리코 유전자 섭동(Virtual Knockout)을 통한 전립선암 t-NEPC 계통 가소성 역전 가상 약물 스크리닝 및 scIB 벤치마킹 체계
  - [생성모델 교재 (Markdown)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html) | [Word 문서 (DOCX)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.docx)
- **[Track 06: 매니폴드 기하학 & 동역학 궤적 추론 (Manifold Geometry & Dynamic Trajectory Inference)](tracks/06_manifold_trajectory/)**
  - 확산 맵(Diffusion Maps)과 Markov 확산 유사시간(DPT), 전사체 스플라이싱 반응속도론(RNA Velocity) 및 결정론적/동역학(Dynamical scVelo) 모델링과 보편 잠재 시간(Latent Time), 엔트로피 정규화 최적 수송(Waddington-OT) 싱크혼 커플링, CellRank 2 다중 커널 융합(속도·유사시간·최적수송·공간), 비가역 마르코프 연쇄의 실 슈어 분해 및 일반화 페론 클러스터 분석(G-PCCA) 거시상태 분해, 흡수 마르코프 연쇄 기본 행렬 기반 전립선암 4대 계통 가소성(AR 선암종, Amphicrine, DNPC, t-NEPC) 분화 운명 흡수 확률 및 마스터 드라이버 전사인자 발굴
  - [매니폴드 궤적 교재 (Markdown)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html) | [Word 문서 (DOCX)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.docx)
- **[Track 07: 단일 세포 카운트 통계학 & 차등 풍부도 (Single-Cell Count Statistics & Differential Abundance)](tracks/07_count_statistics_abundance/)**
  - 과포화 음이항 분포(Negative Binomial vs Poisson), 가짜 반복(Pseudoreplication)의 제1종 오류 대폭발(Type I Error Inflation) 실증, 환자 수준 슈도벌크(Pseudobulk) 집계 및 경험적 베이즈 분산 수축(DESeq2/edgeR), 위계적 일반화 선형 혼합 모델(GLMM) 2단계 아키텍처 및 고속 NEBULA, 분화 궤적 기반 연속 차등 발현(Trajectory DE), Milo kNN 이웃 그래프 기반 차등 풍부도(DA) 검정 및 Spatial FDR, 조성 데이터(CoDa) Simplex 기하학 및 scCODA 베이지안 회귀 분석
  - [카운트 통계학 교재 (Markdown)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.html) | [Word 문서 (DOCX)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.docx)
- **[Track 08: 대규모 희소 연산 & GPU 가속 (Large-Scale Sparse Computation & GPU Acceleration)](tracks/08_large_scale_sparse/)**
  - 희소 행렬 메모리 레이아웃(CSR vs CSC)과 CPU 캐시 라인 국소성, SpMV 연산 복잡도, Out-of-Core 연산과 HDF5 AnnData Backed 모드, Zarr 클라우드 네이티브 청킹 및 Blosc 다중 스레드 병렬 압축, TileDB-SOMA 단일세포 다차원 희소 텐서(타일 기반 저장) 및 Arrow 제로 카피 쿼리, Dask 지연 연산(DAG) 분산 품질 관리, 증분형 주성분 분석(Incremental PCA) 및 Halko 무작위 SVD(Randomized SVD), NVIDIA RAPIDS(cuML, cuGraph, rapids-singlecell) PCIe 스트리밍, CAGRA 근사 kNN 그래프 및 cuGraph Leiden 군집화 50~100배 GPU 가속화
  - [대규모 희소 연산 교재 (Markdown)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.html) | [Word 문서 (DOCX)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.docx)
- **[Track 09: 로우레벨 서열 정렬 & BAM 메커니즘 (Low-Level Sequence Alignment & BAM Mechanics)](tracks/09_upstream_alignment_bam/)**
  - 차세대 염기서열 분석기(NGS) 형광 신호 및 FASTQ 품질 점수(Phred Q), 737K 세포 바코드 화이트리스트 1-Hamming 거리 오류 보정 해시 알고리즘, UMI 중복 제거 그래프 수리 모델(방향성 인접 검정 $N_{\text{parent}} \ge 2N_{\text{child}} - 1$), STARsolo 접미사 배열(Suffix Array) 최대 사상 시드(MMP) 인트론 정렬, Alevin-fry 의사정렬(Pseudoalignment) 및 RAD 바이너리 스트리밍, BGZF 가상 오프셋, CIGAR 문자열, SAM 플래그 비트와 옵션 태그(CB, UB, NH)를 포함한 BAM 바이너리 구조 분해, pysam 프로그래밍 정렬 데이터 마이닝, 전립선암 *AR-V7* 잠재 엑손 3b 스플라이싱 정션(PSI) 및 *TMPRSS2-ERG* 21q22 결손 키메라 리드 단일세포 해상도 검출
  - [업스트림 정렬 교재 (Markdown)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.html) | [Word 문서 (DOCX)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.docx)
- **[Track 10: 확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps (Scalable Production Pipelines, HPC/Cloud Orchestration & MLOps)](tracks/10_production_pipelines_mlops/)**
  - 재현성 위기 극복, Docker vs Singularity/Apptainer HPC 비루트 가상화, Nextflow DSL2 리액티브 채널 및 동적 OOM 자가 치유(Self-Healing), Snakemake 타깃 중심 DAG 역방향 해결 및 벤치마킹, SLURM cgroups 자원 격리 & AWS Batch/GCP 클라우드 분산 오케스트레이션, nf-core/scrnaseq 모듈 엔지니어링, W3C PROV & RO-Crate 기반 FDA 21 CFR Part 11 암호학적 출처 추적성(Data Provenance), 단일세포 전사체 딥러닝 ONNX/Triton 실시간 추론 마이크로서비스 서빙(SLA < 50ms) 및 클라우드 스팟 인스턴스 비용 경제학
  - [프로덕션 파이프라인 교재 (Markdown)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Book.md) | [인터랙티브 웹 교재 (HTML)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.html) | [Word 문서 (DOCX)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.docx)

---

## 1. 과정 소개 (Overview)

단일 핵 RNA 시퀀싱(**snRNA-seq, Single-nucleus RNA sequencing**)은 온전한 단일 세포를 분리하기 어려운 **동결 임상 조직(frozen clinical specimens), 코어 바늘 생검(core needle biopsy), 섬유화된 종양 조직, 뼈 및 림프절 전이암(metastatic prostate cancer)**의 세포 이질성을 극복할 수 있는 전사체 분석 기술입니다.

본 저장소는 **암 생물학(특히 전립선암의 원발암 및 거세저항성 전이암, mCRPC)**을 핵심 모델로 하여:
1. **생물학적/통계학적 이론적 기저(Theoretical Foundations)**
2. **분석 파이프라인 각 단계의 인지적 이해(Analytical Cognition: Why, What, How, Parameter Interpretation)**
3. **재현 가능한 실전 파이프라인 코드(Executable Python/Scanpy Pipelines)**
4. **강의 및 자습용 인터랙티브 슬라이드 덱(Interactive Presentation Decks)**
을 체계적으로 제공합니다.

---

## 2. 디렉토리 구조 (Directory Structure)

```
.
├── README.md                                    # 본 과정 소개 및 종합 안내
├── MASTER_CURRICULUM_ROADMAP.md                 # 11개 트랙(Track 00~10) 전체 커리큘럼 종합 로드맵
├── snRNAseq_Prostate_Master_Book.md             # [Track 00 Master 1] 전체 9개 단원 통합 마크다운 교재
├── snRNAseq_Prostate_Master_Course.html         # [Track 00 Master 2] 인터랙티브 사이드바/검색/다크모드/PDF 지원 단일 HTML
├── snRNAseq_Prostate_Master_Course.docx         # [Track 00 Master 3] 수식 유니코드 & 다이어그램 콜아웃 OpenXML Word 문서
├── docs/                                        # [Track 00] 상세 챕터별 교재 및 가이드 (Markdown)
│   ├── 00_CURRICULUM_AND_ROADMAP.md             # 학습 로드맵, 전립선암 공개 벤치마크 데이터셋(GSE143791 등) 카탈로그
│   ├── 01_FOUNDATIONS_AND_WET_LAB.md            # 전립선암 조직 특성과 snRNA-seq, 핵 분리 원리, sc vs sn 비교
│   ├── 02_UPSTREAM_PROCESSING.md                # 전립선암 FASTQ, Cell Ranger (--include-introns), Sparse Matrix 인지
│   ├── 03_QC_AND_CLEANING.md                    # 종양 조직 특화 QC, Ambient RNA (SoupX/CellBender), Doublets
│   ├── 04_CORE_DOWNSTREAM_PIPELINE.md           # 정규화, HVG, PCA/UMAP의 기하학, Leiden 클러스터링, DEG 원리
│   ├── 05_CELL_TYPE_ANNOTATION.md               # 전립선 상피/기저/TME 마커, 악성 판별 (inferCNV, Numbat BAF 검증)
│   ├── 06_ADVANCED_INTEGRATION_AND_NETWORKS.md  # 원발/전이 Batch 통합 (scVI, Harmony), 4대 가소성 궤적, pySCENIC
│   ├── 07_MULTIOMICS_AND_FUTURE.md              # 전립선암 snRNA + snATAC Multiome, 공간전사체(Visium/Xenium) 플랫폼 비교
│   ├── 08_TROUBLESHOOTING_AND_BEST_PRACTICES.md # 암 데이터 특유의 오류 해결, Pseudobulk 원리, 캡스톤 프로젝트 루브릭
│   └── images/                                  # 8종 고해상도 인포그래픽 벡터 다이어그램 (SVG)
│       ├── 01_wetlab_nuclear_isolation_workflow.svg      # 전립선암 동결 조직 단일 핵 분리 및 OptiPrep 원심분리
│       ├── 02_premRNA_include_introns_mechanics.svg      # Pre-mRNA 인트론 캡처 및 Cell Ranger 카운팅 메커니즘
│       ├── 03_snRNA_qc_mitochondrial_and_cellbender.svg  # 종양 snRNA-seq 특화 QC (Mito < 2%) 및 CellBender 배경 제거
│       ├── 04_infercnv_chromosomal_aneuploidy.svg        # inferCNV 염색체 이수성(Aneuploidy) 판별 원리
│       ├── 05_mcrpc_lineage_plasticity_trajectory.svg    # mCRPC 치료 저항성 및 4대 계통 가소성(AR+, Amphicrine, t-NEPC, DNPC)
│       ├── 06_multiome_wnn_and_spatial_deconvolution.svg # 10x Multiome(snRNA+snATAC) WNN 및 공간전사체 디콘볼루션
│       ├── 07_downstream_manifold_and_leiden.svg         # 차원축소(PCA/kNN/UMAP) 매니폴드 및 Leiden 클러스터링 원리
│       └── 08_pseudobulk_vs_pseudoreplication.svg        # 가짜 반복(Pseudoreplication) vs 환자 수준 Pseudobulk 통계 비교
├── slides/                                      # 브라우저 즉시 실행 가능한 인터랙티브 슬라이드 (HTML)
│   ├── 01_Foundations_and_QC_Deck.html          # Part 1: 전립선암 배경, 실험 원리, 업스트림, 종양 특화 QC
│   ├── 02_Downstream_Analysis_Deck.html         # Part 2: 핵심 파이프라인(차원축소, 클러스터링, 어노테이션, inferCNV)
│   └── 03_Advanced_and_Multiomics_Deck.html     # Part 3: 통합, AR/NEPC 분화 궤적, 전사조절망, 멀티오믹스, Best Practices
├── pipelines/                                   # 교육용 실습 파이프라인 템플릿 (Python, 합성 데이터 기반; 환경: requirements-pipelines.txt)
│   ├── 01_prostate_standard_scanpy_pipeline.py  # 전립선암 표준 Scanpy 파이프라인
│   ├── 02_prostate_qc_and_cellbender.py         # Ambient RNA 및 Doublet 제거 정제 파이프라인
│   ├── 03_prostate_infercnv_and_pseudobulk.py   # 악성 세포 판별 및 Pseudobulk DEG 분석 파이프라인
│   └── 04_prostate_integration_trajectory_and_grn.py # 다중 배치 통합, 4대 가소성 궤적, pySCENIC 조절망 & 골전이 L-R 상호작용
├── scripts/
│   ├── build_master_documents.py                # Track 00 마스터 교재 자동 컴파일 빌더 스크립트
│   ├── generate_diagrams.py                     # Track 00 SVG 다이어그램 8종 일괄 생성 스크립트
│   ├── inject_global_track_switcher.py          # 11개 트랙 HTML 간 글로벌 스위처/배너 주입기
│   └── verify_switcher_links.py                 # 스위처 링크 무결성 자동 검증 스크립트
└── tracks/                                      # 심화 및 차세대 오믹스 트랙 모음
    ├── 01_multiome_snatac/                      # [Track 01] 10x Multiome (snRNA + snATAC) & 후성유전체학
        ├── Multiome_snATAC_Master_Book.md       # [Track 01 Master 1] 멀티오믹스 10개 챕터 통합 교재
        ├── Multiome_snATAC_Master_Course.html   # [Track 01 Master 2] 인터랙티브 웹 통합 교재
        ├── Multiome_snATAC_Master_Course.docx   # [Track 01 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 01] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 01] 멀티오믹스 QC/LSI 및 P2G/eGRN 파이프라인
        ├── slides/                              # [Track 01] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 01] 다이어그램 생성 및 마스터 빌더 스크립트
    ├── 02_spatial_transcriptomics/              # [Track 02] 공간 전사체학 (Visium & Xenium)
        ├── Spatial_Transcriptomics_Master_Book.md   # [Track 02 Master 1] 공간전사체 10개 챕터 통합 교재
        ├── Spatial_Transcriptomics_Master_Course.html # [Track 02 Master 2] 인터랙티브 웹 통합 교재
        ├── Spatial_Transcriptomics_Master_Course.docx # [Track 02 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 02] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 02] 공간 QC/디콘볼루션 및 Xenium 니치/신호망 파이프라인
        ├── slides/                              # [Track 02] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 02] 마스터 빌더 및 검증 스크립트
    ├── 03_tme_cell_communication/               # [Track 03] 종양 미세환경 & 세포 간 상호작용
        ├── TME_Cell_Communication_Master_Book.md # [Track 03 Master 1] 세포대화 10개 챕터 통합 교재
        ├── TME_Cell_Communication_Master_Course.html # [Track 03 Master 2] 인터랙티브 웹 통합 교재
        ├── TME_Cell_Communication_Master_Course.docx # [Track 03 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 03] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 03] LIANA/CellChat 및 NicheNet/공간 플럭스 파이프라인
        ├── slides/                              # [Track 03] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 03] 마스터 빌더 및 검증 스크립트
    ├── 04_cancer_genomics_clones/               # [Track 04] 단일 세포 암 유전체 & 클론 진화
        ├── Cancer_Genomics_Master_Book.md       # [Track 04 Master 1] 암유전체 10개 챕터 통합 교재
        ├── Cancer_Genomics_Master_Course.html   # [Track 04 Master 2] 인터랙티브 웹 통합 교재
        ├── Cancer_Genomics_Master_Course.docx   # [Track 04 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 04] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 04] inferCNV/CopyKAT 및 Numbat/계통수 파이프라인
        ├── slides/                              # [Track 04] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 04] 마스터 빌더 및 검증 스크립트
    ├── 05_generative_foundation_models/         # [Track 05] 단일 세포 생성형 딥러닝 & 파운데이션 모델
        ├── Generative_Foundation_Master_Book.md # [Track 05 Master 1] 생성모델 10개 챕터 통합 교재
        ├── Generative_Foundation_Master_Course.html # [Track 05 Master 2] 인터랙티브 웹 통합 교재
        ├── Generative_Foundation_Master_Course.docx # [Track 05 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 05] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 05] scVI/scANVI 및 트랜스포머/가상섭동 파이프라인
        ├── slides/                              # [Track 05] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 05] 마스터 빌더 및 검증 스크립트
    ├── 06_manifold_trajectory/                  # [Track 06] 매니폴드 기하학 & 궤적 추론
        ├── Manifold_Trajectory_Master_Book.md   # [Track 06 Master 1] 매니폴드 궤적 10개 챕터 통합 교재
        ├── Manifold_Trajectory_Master_Course.html # [Track 06 Master 2] 인터랙티브 웹 통합 교재
        ├── Manifold_Trajectory_Master_Course.docx # [Track 06 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 06] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 06] 확산/벨로시티 및 CellRank 2 운명 매핑 파이프라인
        ├── slides/                              # [Track 06] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 06] 마스터 빌더 및 검증 스크립트
    ├── 07_count_statistics_abundance/           # [Track 07] 카운트 통계학 & 차등 풍부도
        ├── Count_Statistics_Master_Book.md      # [Track 07 Master 1] 통계학/차등분석 10개 챕터 통합 교재
        ├── Count_Statistics_Master_Course.html  # [Track 07 Master 2] 인터랙티브 웹 통합 교재
        ├── Count_Statistics_Master_Course.docx  # [Track 07 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 07] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 07] 슈도벌크/GLMM 및 Milo/CoDa 파이프라인
        ├── slides/                              # [Track 07] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 07] 마스터 빌더 및 검증 스크립트
    ├── 08_large_scale_sparse/                   # [Track 08] 대규모 희소 연산 & GPU 가속
        ├── Large_Scale_Sparse_Master_Book.md    # [Track 08 Master 1] 대규모 희소 연산 10개 챕터 통합 교재
        ├── Large_Scale_Sparse_Master_Course.html # [Track 08 Master 2] 인터랙티브 웹 통합 교재
        ├── Large_Scale_Sparse_Master_Course.docx # [Track 08 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 08] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 08] 희소/Out-of-Core 및 TileDB/RAPIDS 파이프라인
        ├── slides/                              # [Track 08] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 08] 마스터 빌더 및 검증 스크립트
    ├── 09_upstream_alignment_bam/               # [Track 09] 로우레벨 서열 정렬 & BAM 메커니즘
        ├── Upstream_Alignment_Master_Book.md    # [Track 09 Master 1] 업스트림 정렬 10개 챕터 통합 교재
        ├── Upstream_Alignment_Master_Course.html # [Track 09 Master 2] 인터랙티브 웹 통합 교재
        ├── Upstream_Alignment_Master_Course.docx # [Track 09 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 09] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 09] 바코드/UMI 디덥 및 BAM 스플라이스/융합 파이프라인
        ├── slides/                              # [Track 09] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 09] 마스터 빌더 및 검증 스크립트
    └── 10_production_pipelines_mlops/           # [Track 10] 프로덕션 파이프라인, HPC/클라우드 & MLOps
        ├── Production_Pipelines_Master_Book.md  # [Track 10 Master 1] 프로덕션 파이프라인 10개 챕터 통합 교재
        ├── Production_Pipelines_Master_Course.html # [Track 10 Master 2] 인터랙티브 웹 통합 교재
        ├── Production_Pipelines_Master_Course.docx # [Track 10 Master 3] OpenXML Word 교재
        ├── docs/                                # [Track 10] 챕터별 상세 교재 (00~09) 및 SVG 다이어그램
        ├── pipelines/                           # [Track 10] DSL2 엔진 및 MLOps 서빙/감사 파이프라인
        ├── slides/                              # [Track 10] 인터랙티브 프레젠테이션 슬라이드 (01~02)
        └── scripts/                             # [Track 10] 마스터 빌더 및 검증 스크립트
```

---

## 3. 학습 단계별 로드맵 (Recommended Learning Paths)

| 단계 | 목표 | 주요 내용 | 추천 교재 |
|:---|:---|:---|:---|
| **Step 1: 기초 및 실험 원리** | 전립선암 조직에서 왜 snRNA-seq이 필수적인지 인지 | scRNA-seq vs snRNA-seq 비교, 핵 분리(Nuclear isolation) 메커니즘, pre-mRNA 인트론 캡처, 공개 벤치마크 데이터셋(GSE143791 등) 탐색 | `docs/00`, `docs/01` |
| **Step 2: 업스트림 및 데이터 정제** | 원시 서열 데이터에서 신뢰할 수 있는 카운트 매트릭스 도출 | Cell Ranger `--include-introns`, 미토콘드리아 QC(핵 특이성 1~2%), CellBender(Ambient RNA 제거), Scrublet(Doublet) | `docs/02`, `docs/03` |
| **Step 3: 핵심 다운스트림 분석** | 세포 아틀라스 구축 및 악성 상피세포 분리 | Log-norm vs SCTransform, UMAP/Leiden 기하학, 전립선 TME 마커, inferCNV 및 Numbat BAF 대립유전자 분석을 통한 Aneuploid 암세포 선별 | `docs/04`, `docs/05` |
| **Step 4: 고급 전이암 분석** | 치료 저항성 기전 및 다중 샘플 통합 | Harmony vs scVI 배치 통합 및 scIB 과보정 평가, 4대 계통 가소성(AR+, Amphicrine, t-NEPC, DNPC) 궤적, pySCENIC 마스터 레귤론 | `docs/06` |
| **Step 5: 멀티오믹스 및 실전 검증** | 차세대 오믹스 확장 및 통계적 검증 | 10x Multiome(snRNA+snATAC), 공간전사체 플랫폼(Visium, Visium HD, Xenium) 비교, Pseudobulk 기반 환자 공변량 통제 모델링 및 캡스톤 루브릭 | `docs/07`, `docs/08` |

> 💡 **단원 핵심 점검 (Cognitive Checkpoint & Case Study)**:
> `docs/01`부터 `docs/08`까지 모든 단원 말미에는 단순 암기를 넘어 실제 연구 현장에서 마주치는 문제를 해결할 수 있도록 **3문항의 능동 회상(Active Recall) 및 임상 케이스 해석 퀴즈**가 수록되어 있으며, 접이식 카드(`<details>`)를 통해 즉시 상세 모범 답안을 확인하고 자가 진단할 수 있습니다.

---

## 4. 마스터 통합 교재 및 슬라이드 활용법

### 1) 인터랙티브 마스터 HTML 교재 (`snRNAseq_Prostate_Master_Course.html`)
- **실시간 키워드 검색**: 좌측 상단 검색창에 `KLK3`, `CellBender`, `DNPC`, `inferCNV`, `Pseudobulk` 등을 입력하면 즉시 관련 섹션과 문맥으로 이동합니다.
- **다크 / 라이트 모드 전환**: 상단 토글 버튼으로 야간 및 장시간 학습 시 시각적 피로도를 최소화합니다.
- **퀴즈 정답 일괄 제어**: `💡 모든 퀴즈 정답 펼치기/접기` 버튼으로 학습 및 테스트 모드를 손쉽게 전환할 수 있습니다.
- **완벽한 PDF 인쇄 지원**: `Cmd + P` 또는 `🖨️ PDF로 저장 / 인쇄` 클릭 시 사이드바가 자동 제거되고 모든 퀴즈 해설이 펼쳐져 고품질 PDF 책자로 출력됩니다.

### 2) 프레젠테이션 슬라이드 (`slides/`)
`slides/` 디렉토리에 위치한 슬라이드는 **외부 인터넷 연결이나 추가 소프트웨어 설치 없이** 크롬, 사파리, 엣지 등의 웹 브라우저에서 바로 열 수 있는 반응형 프레젠테이션입니다.
- **슬라이드 이동**: 키보드 `방향키(← / →)` 또는 `Space`, `PageUp / PageDown`
- **전체 화면 전환**: `F` 키

---

## 5. 실습 코드 실행 환경 (Execution Environment)

```bash
# 권장 Conda / Mamba 가상환경 구축
conda create -n snrnaseq python=3.10 -y
conda activate snrnaseq

# 핵심 패키지 설치
pip install scanpy anndata scrublet cellbender scvi-tools pyscenic decoupler infercnvpy scikit-misc
pip install matplotlib seaborn leidenalg
```
상세한 설치 및 패키지 의존성 안내는 [docs/00_CURRICULUM_AND_ROADMAP.md](docs/00_CURRICULUM_AND_ROADMAP.md)에서 확인하실 수 있습니다.

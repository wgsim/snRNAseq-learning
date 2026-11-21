# Grand Curriculum Roadmap: Single-Cell & Spatial Omics Master Ecosystem

> **From Single Nucleus to Spatial Tissue Ecosystems & Production MLOps**: 이론적 분자생물학·수리통계학적 기저부터 프로덕션 MLOps 파이프라인까지 포괄하는 단일 세포(Single-Cell) 및 공간 오믹스(Spatial Omics) 11대 트랙(Track 00 ~ 10) 100% 완주 종합 마스터 로드맵

---

## 1. 커리큘럼 아키텍처 개요 (Ecosystem Architecture)

본 단일 세포 & 공간 오믹스 학습 생태계는 연구 현장에서 흔히 발생하는 "도구 사용법(API) 암기 위주의 단편적 학습"을 탈피하고, **생물학적 원인(Why) ──► 수리통계학적 모델링(What) ──► 최적의 전산 구현(How) ──► 파라미터 임상 해석(Interpretation)**으로 이어지는 인지적 도제 모델(Cognitive Apprenticeship)을 지향합니다.

원발암(Primary Adenocarcinoma)부터 거세저항성 전이암(mCRPC, 골/림프절 전이)에 이르는 전립선암 임상 모델을 척추로 삼아, 11개 트랙(Track 00 ~ 10)이 유기적으로 연결된 다차원 오믹스 전주기 학습 생태계를 100% 완성하여 제공합니다:

```
[ Track 00: snRNA-seq 기초 & 핵심 다운스트림 파이프라인 ]  <-- COMPLETED
                         │
                         ▼
[ Track 01: 10x Multiome (snRNA + snATAC) & 후성유전체학 ]  <-- COMPLETED
                         │
         ┌───────────────┴───────────────┐
         ▼                               ▼
[ Track 02: 공간 전사체학 ]  <-- COMPLETED
(Visium, Xenium, CosMx)         [ Track 04: 암 유전체 & 클론 진화 ]  <-- COMPLETED
                                (inferCNV, CopyKAT, Numbat)
         │                               │
         ▼                               ▼
[ Track 03: 종양미세환경 & 세포대화 ]  <-- COMPLETED
(CellChat, LIANA, NicheNet)     [ Track 05: 생성 모델 & 파운데이션 모델 ]  <-- COMPLETED
         │                               │
         └───────────────┬───────────────┘
                         ▼
[ Track 06: 매니폴드 기하학 & 궤적 추론 (CellRank, RNA Velocity, OT) ]  <-- COMPLETED
                         │
                         ▼
[ Track 07: 단일세포 통계학 & 차등 풍부도 (GLMM, Negative Binomial, Milo) ]  <-- COMPLETED
                         │
                         ▼
[ Track 08: 대규모 희소 연산 & GPU 가속 (TileDB, RAPIDS cuGraph) ]  <-- COMPLETED
                         │
                         ▼
[ Track 09: 로우레벨 서열 정렬 & BAM 메커니즘 (STARsolo, Alevin-fry) ]  <-- COMPLETED
                         │
                         ▼
[ Track 10: 프로덕션 파이프라인 & MLOps (Nextflow, Snakemake, Docker) ]  <-- COMPLETED
```

---

## 2. 11대 마스터 트랙 로드맵 총괄표 (Track Catalog: 100% COMPLETED)

| 트랙 ID | 트랙 명칭 및 핵심 주제 | 주요 도구 및 알고리즘 | 상태 | 핵심 산출물 및 문서 링크 |
| :---: | :--- | :--- | :---: | :--- |
| **Track 00** | **snRNA-seq 기초 & 핵심 다운스트림**<br>동결 조직 핵 분리, Pre-mRNA, QC, 차원축소, 어노테이션, inferCNV | `Cell Ranger`, `CellBender`, `Scanpy`, `inferCNV`, `scVI`, `pySCENIC` | **완료** | • [통합 교재 (Markdown)](snRNAseq_Prostate_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](snRNAseq_Prostate_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](snRNAseq_Prostate_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](docs/)<br>• [프레젠테이션 슬라이드 (slides/)](slides/)<br>• [실행 파이프라인 (pipelines/)](pipelines/) |
| **Track 01** | **10x Multiome & 후성유전체학**<br>snATAC-seq, Tn5 화학, Fragment 파일, LSI, WNN 통합, P2G, SCENIC+ | `Cell Ranger ARC`, `Chromap`, `MACS3`, `SnapATAC2`, `MuData`, `Cicero`, `SCENIC+` | **완료** | • [통합 교재 (Markdown)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/01_multiome_snatac/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/01_multiome_snatac/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/01_multiome_snatac/pipelines/) |
| **Track 02** | **공간 전사체학 (Spatial Transcriptomics)**<br>시퀀싱 기반(Visium, Visium HD) vs 영상 기반(Xenium, CosMx, MERSCOPE) | `Seurat Spatial`, `Squidpy`, `SpatialData`, `Cellpose`, `Baysor`, `RCTD`, `cell2location`, `COMMOT` | **완료** | • [통합 교재 (Markdown)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/02_spatial_transcriptomics/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/02_spatial_transcriptomics/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/02_spatial_transcriptomics/pipelines/) |
| **Track 03** | **종양미세환경 & 세포 간 상호작용 (TME & Cell-Cell Communication)**<br>리간드-수용체 복합체, 다중 소단위체, 공간 근접성, 인과 추론, 골전이 회로 | `CellChat v2`, `LIANA+`, `NicheNet v2`, `MISTy`, `CytoSig` | **완료** | • [통합 교재 (Markdown)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/03_tme_cell_communication/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/03_tme_cell_communication/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/03_tme_cell_communication/pipelines/) |
| **Track 04** | **단일 세포 암 유전체 & 클론 진화**<br>복제수 변이(CNA), 대립유전자 특이적 발현(BAF), 계통수 | `inferCNV`, `CopyKAT`, `Numbat`, `MEDICC2`, `SHAPEIT4` | **완료** | • [통합 교재 (Markdown)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/04_cancer_genomics_clones/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/04_cancer_genomics_clones/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/04_cancer_genomics_clones/pipelines/) |
| **Track 05** | **딥 생성 모델 & 파운데이션 모델 (Generative & Foundation Models)**<br>변분 오토인코더(scVI), 준지도 세포형 전이(scANVI), 지속적 아틀라스 수술(scArches), 유전자 언어 트랜스포머(scGPT, Geneformer), 인 실리코 유전자 섭동 | `scVI`, `scANVI`, `scArches`, `scGPT`, `Geneformer`, `BioGPT`, `AnnData` | **완료** | • [통합 교재 (Markdown)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/05_generative_foundation_models/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/05_generative_foundation_models/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/05_generative_foundation_models/pipelines/) |
| **Track 06** | **매니폴드 기하학 & 궤적 분화 모델링 (Manifold Geometry & Dynamic Trajectory Inference)**<br>확산 맵, RNA 벨로시티, 결정론적/동역학 모델링, 최적 수송(WOT), CellRank 2 다중 커널 융합, G-PCCA 거시상태, 흡수 마르코프 체인 운명 확률, 계통 드라이버 발굴 | `Scanpy`, `scVelo`, `CellRank 2`, `Waddington-OT`, `GPCCA` | **완료** | • [통합 교재 (Markdown)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/06_manifold_trajectory/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/06_manifold_trajectory/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/06_manifold_trajectory/pipelines/) |
| **Track 07** | **단일 세포 카운트 통계학 & 차등 분석**<br>과포화 음이항 분포, GLMM 위계 모델, 세포 풍부도 | `EdgeR-LRT`, `DESeq2-LRT`, `glmmTMB`, `NEBULA`, `Milo`, `scCODA` | **완료** | • [통합 교재 (Markdown)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/07_count_statistics_abundance/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/07_count_statistics_abundance/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/07_count_statistics_abundance/pipelines/) |
| **Track 08** | **대규모 희소 연산 & GPU 가속**<br>Out-of-Core 데이터 구조, 병렬화, GPU 가속 커널 | `AnnData (backed)`, `TileDB-SOMA`, `Zarr`, `RAPIDS cuGraph/cuml`, `GPU Scanpy` | **완료** | • [통합 교재 (Markdown)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/08_large_scale_sparse/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/08_large_scale_sparse/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/08_large_scale_sparse/pipelines/) |
| **Track 09** | **로우레벨 서열 정렬 & BAM 메커니즘**<br>바코드 복구, UMI 에러 교정, 유사 정렬, C++ 엔진 | `STARsolo`, `Alevin-fry`, `samtools`, `pysam`, `BAM parsing`, `Isoform Mining` | **완료** | • [통합 교재 (Markdown)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/09_upstream_alignment_bam/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/09_upstream_alignment_bam/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/09_upstream_alignment_bam/pipelines/) |
| **Track 10** | **프로덕션 파이프라인 & 생물정보 MLOps**<br>재현 가능한 컨테이너 파이프라인, Nextflow DSL2, Snakemake, HPC SLURM, AWS Batch, nf-core/scrnaseq, 출처 추적성(W3C PROV/RO-Crate), MLOps 모델 서빙(ONNX/Triton) | `Nextflow (DSL2)`, `Snakemake`, `Docker`, `Singularity/Apptainer`, `SLURM`, `AWS Batch`, `nf-core`, `RO-Crate`, `ONNX Runtime`, `Triton` | **완료** | • [통합 교재 (Markdown)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Book.md)<br>• [인터랙티브 웹 교재 (HTML)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.html)<br>• [인쇄/출력용 문서 (DOCX)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.docx)<br>• [챕터별 상세 문서 (docs/)](tracks/10_production_pipelines_mlops/docs/)<br>• [프레젠테이션 슬라이드 (slides/)](tracks/10_production_pipelines_mlops/slides/)<br>• [실행 파이프라인 (pipelines/)](tracks/10_production_pipelines_mlops/pipelines/) |

---

## 3. 학습자 유형별 권장 트랙 가이드 (Persona-based Learning Paths)

모든 학습자가 11개 트랙을 순차적으로 학습할 필요는 없습니다. 개인의 배경 지식과 연구 목적에 맞추어 최적화된 학습 경로를 권장합니다:

### 3.1 임상 연구자 / 종양내과·병리과 전문의 (Clinical Researchers & Pathologists)
- **추천 경로**: **Track 00 ──► Track 01 ──► Track 02 ──► Track 03 ──► Track 04**
- **학습 중점**:
  - 생검 조직의 핵 분리 품질 관리(TSS, Mito)와 환자 검체의 기술적 한계 인지.
  - 치료 저항성 계통 가소성(t-NEPC, DNPC) 및 복제수 변이(CNA)에 기반한 악성 클론 감별.
  - 공간 전사체(Visium/Xenium)를 통한 종양 미세환경 내 면역 회피 및 리간드-수용체 상호작용 파악.
  - 임상 데이터 분석 시의 환자 간 교란 변수 통제(Pseudobulk 모델링).

### 3.2 전산생물학자 / 생물정보학 데이터 과학자 (Computational Biologists & Bioinformaticians)
- **추천 경로**: **Track 00 ──► Track 01 ──► Track 05 ──► Track 06 ──► Track 07**
- **학습 중점**:
  - 차원 축소(PCA, LSI, UMAP)의 기하학적 원리와 WNN 멀티모달 가중 결합의 수학적 유도.
  - 딥 변분 오토인코더(scVI)와 파운데이션 모델(scGPT)을 활용한 다기관 배치 통합.
  - 최적 수송(Optimal Transport) 및 마르코프 체인(CellRank 2)을 이용한 세포 상태 전이 궤적 복원.
  - 위계적 일반화 선형 혼합 모델(GLMM)을 통한 위양성 없는 엄밀한 통계 검정.

### 3.3 소프트웨어 엔지니어 / 생명정보 MLOps 엔지니어 (Software & MLOps Engineers)
- **추천 경로**: **Track 00 ──► Track 01 ──► Track 08 ──► Track 09 ──► Track 10**
- **학습 중점**:
  - FASTQ부터 단편 파일까지 로우레벨 C++ 정렬 엔진(STARsolo, Chromap) 및 바코드 에러 보정.
  - TileDB-SOMA 및 Zarr 기반 Out-of-Core 대규모 세포 매트릭스 아키텍처.
  - RAPIDS GPU 가속을 통한 분산 환경 클러스터링 최적화.
  - Nextflow/Snakemake 및 컨테이너(Docker/Singularity)를 통한 페타바이트급 파이프라인 자동화.

---

## 4. 완주 트랙 심층 요약 및 주요 산출물 가이드

### 4.1 Track 00: snRNA-seq Foundations & Core Downstream Pipeline
- **교육 목표**: 온전한 단일 세포 분리가 어려운 동결 전립선암 검체에서 핵(Nucleus)을 추출하고, Pre-mRNA 인트론을 보존하여 고품질 전사체 아틀라스를 구축합니다.
- **주요 구성**:
  - `docs/00_CURRICULUM_AND_ROADMAP.md` ~ `docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (9개 챕터 완비)
  - 8종의 고해상도 인포그래픽 SVG 다이어그램 (`docs/images/`)
  - 4종의 엔드-투-엔드 Scanpy 실습 파이프라인(기본 합성 데이터, 권장 환경은 `requirements-pipelines.txt`) (`pipelines/`)
  - 3종의 인터랙티브 프레젠테이션 슬라이드 덱 (`slides/`)
- **종합 마스터 문서**:
  - [마스터 마크다운 교재 (snRNAseq_Prostate_Master_Book.md)](snRNAseq_Prostate_Master_Book.md)
  - [인터랙티브 웹 교재 (snRNAseq_Prostate_Master_Course.html)](snRNAseq_Prostate_Master_Course.html) (검색, 다크모드, 퀴즈 토글, PDF 인쇄 최적화)
  - [OpenXML 워드 문서 (snRNAseq_Prostate_Master_Course.docx)](snRNAseq_Prostate_Master_Course.docx) (유니코드 수식 & 다이어그램 콜아웃 완비)

### 4.2 Track 01: 10x Multiome (snRNA + snATAC) & Epigenomics
- **교육 목표**: 동일한 단일 핵에서 전사체와 염색질 개방성을 동시에 포획하여, 원거리 인핸서-프로모터 루핑과 전사 인자 풋프린팅을 결합한 인과적 전사조절망(SCENIC+)을 복원합니다.
- **주요 구성**:
  - `tracks/01_multiome_snatac/docs/00_MULTIOME_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 후성유전학 벡터 SVG 다이어그램 (`tracks/01_multiome_snatac/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 멀티오믹스 실습 파이프라인 (`tracks/01_multiome_snatac/pipelines/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [멀티오믹스 마스터 마크다운 교재 (Multiome_snATAC_Master_Book.md)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Book.md)
  - [인터랙티브 웹 교재 (Multiome_snATAC_Master_Course.html)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.html)
  - [OpenXML 워드 문서 (Multiome_snATAC_Master_Course.docx)](tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.docx)

### 4.3 Track 02: Spatial Transcriptomics (Visium, Visium HD, Xenium, CosMx)
- **교육 목표**: 해리 시퀀싱의 한계를 넘어 전립선암의 2차원/3차원 병리 아키텍처(Gleason 패턴, Cribriform 선암종, CAF 물리적 장벽 및 면역 배제)를 보존하고, 시퀀싱 기반 및 이미징 기반 공간전사체 데이터를 체계적으로 분석합니다.
- **주요 구성**:
  - `tracks/02_spatial_transcriptomics/docs/00_SPATIAL_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 공간생물학 벡터 SVG 다이어그램 (`tracks/02_spatial_transcriptomics/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 공간 오믹스 실습 파이프라인 (`tracks/02_spatial_transcriptomics/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/02_spatial_transcriptomics/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [공간전사체 마스터 마크다운 교재 (Spatial_Transcriptomics_Master_Book.md)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Book.md)
  - [인터랙티브 웹 교재 (Spatial_Transcriptomics_Master_Course.html)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.html)
  - [OpenXML 워드 문서 (Spatial_Transcriptomics_Master_Course.docx)](tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.docx)

### 4.4 Track 03: Tumor Microenvironment & Cell-Cell Communication (TME & CCC)
- **교육 목표**: 종양 미세환경 내 이질적 세포 간 상호작용을 단일 유전자 마커 수준을 넘어 리간드-수용체 복합체 물리화학, 다중 소단위체 기하평균 모델, 세포 내 하류 표적인과 추론(NicheNet) 및 공간 거리 감쇄 모델(LIANA+/MISTY)로 정량화하고, 전립선암 골 전이 치명적 악순환 회로의 표적 치료 전략을 규명합니다.
- **주요 구성**:
  - `tracks/03_tme_cell_communication/docs/00_TME_COMMUNICATION_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 세포 간 통신 및 악순환 회로 벡터 SVG 다이어그램 (`tracks/03_tme_cell_communication/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 세포 간 통신 실습 파이프라인 (`tracks/03_tme_cell_communication/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/03_tme_cell_communication/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [세포대화 마스터 마크다운 교재 (TME_Cell_Communication_Master_Book.md)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Book.md)
  - [인터랙티브 웹 교재 (TME_Cell_Communication_Master_Course.html)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.html)
  - [OpenXML 워드 문서 (TME_Cell_Communication_Master_Course.docx)](tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.docx)

### 4.5 Track 04: Single-Cell Cancer Genomics & Clonal Evolution
- **교육 목표**: 단순 전사체 마커 발현의 한계를 넘어 염색체 복제수 변이(CNA) 및 이수성(Aneuploidy)을 정량화하여 진정한 악성 종양 세포를 감별하고, 대립유전자 특이적 발현(BAF)과 하플로타입 위상차를 융합하여 복제중립 LOH를 검출하며, 단일 세포 계통수(NJ/MEDICC2) 복원 및 치료 저항성 계통 가소성(t-NEPC)과 다기관 전이 파종(Seeding) 모델을 규명합니다.
- **주요 구성**:
  - `tracks/04_cancer_genomics_clones/docs/00_CANCER_GENOMICS_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 암 유전체학 및 클론 진화 벡터 SVG 다이어그램 (`tracks/04_cancer_genomics_clones/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 유전체/계통수 실습 파이프라인 (`tracks/04_cancer_genomics_clones/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/04_cancer_genomics_clones/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [암유전체 마스터 마크다운 교재 (Cancer_Genomics_Master_Book.md)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Book.md)
  - [인터랙티브 웹 교재 (Cancer_Genomics_Master_Course.html)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.html)
  - [OpenXML 워드 문서 (Cancer_Genomics_Master_Course.docx)](tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.docx)

### 4.6 Track 05: Deep Generative Modeling & Single-Cell Foundation Models
- **교육 목표**: 고전적 log(1+x) 정규화와 선형 배치 보정(PCA, Harmony)의 수리통계적 결함을 극복하고, 원시 카운트 직접 모델링(NB/ZINB) 기반의 심층 변분 오토인코더(scVI)와 준지도 세포형 전이(scANVI)를 구축하며, 전신 아틀라스 재학습 없는 지속적 아키텍처 수술(scArches)을 적용합니다. 나아가 유전자 언어 모델(scGPT, Geneformer)의 다중헤드 자기주의집중(Self-Attention)으로 창발적 유전자 조절망(Emergent GRN)을 복원하고, 인 실리코 유전자 섭동(In Silico Virtual Knockout)을 통해 거세저항성 전립선암의 치료 저항성 계통 가소성(t-NEPC) 역전 가상 약물 스크리닝을 수행합니다.
- **주요 구성**:
  - `tracks/05_generative_foundation_models/docs/00_FOUNDATION_MODELS_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 생성 딥러닝 및 파운데이션 트랜스포머 벡터 SVG 다이어그램 (`tracks/05_generative_foundation_models/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 생성 모델 실습 파이프라인 (`tracks/05_generative_foundation_models/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/05_generative_foundation_models/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [생성 파운데이션 마스터 마크다운 교재 (Generative_Foundation_Master_Book.md)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Book.md)
  - [인터랙티브 웹 교재 (Generative_Foundation_Master_Course.html)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html)
  - [OpenXML 워드 문서 (Generative_Foundation_Master_Course.docx)](tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.docx)

### 4.7 Track 06: Manifold Geometry & Dynamic Trajectory Inference (CellRank 2, RNA Velocity, Waddington-OT)
- **교육 목표**: 고전적 선형 유사시간의 위상 왜곡을 극복하고, 확산 맵(Diffusion Maps)과 Markov 확산 유사시간(DPT)을 정립하며, 전사체 스플라이싱 반응속도론(RNA Velocity)의 결정론적/동역학(Dynamical scVelo) 모델링과 보편 잠재 시간(Latent Time)을 유도합니다. 나아가 시간 경과 시계열 데이터의 엔트로피 최적 수송(Waddington-OT) 및 CellRank 2 다중 커널(속도·유사시간·최적수송·공간) 융합 프레임워크를 구축하고, 비가역 마르코프 연쇄의 실 슈어 분해(G-PCCA)와 흡수 마르코프 연쇄 기본 행렬을 통해 거세저항성 전립선암(mCRPC)의 4대 계통 가소성(AR 선암종, Amphicrine, DNPC, t-NEPC) 분화 운명 흡수 확률 및 마스터 드라이버 전사인자를 발굴합니다.
- **주요 구성**:
  - `tracks/06_manifold_trajectory/docs/00_TRAJECTORY_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 매니폴드 기하학 및 동역학 궤적 벡터 SVG 다이어그램 (`tracks/06_manifold_trajectory/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 확산/벨로시티 및 CellRank 2 운명 매핑 실습 파이프라인 (`tracks/06_manifold_trajectory/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/06_manifold_trajectory/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [매니폴드 궤적 마스터 마크다운 교재 (Manifold_Trajectory_Master_Book.md)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Book.md)
  - [인터랙티브 웹 교재 (Manifold_Trajectory_Master_Course.html)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html)
  - [OpenXML 워드 문서 (Manifold_Trajectory_Master_Course.docx)](tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.docx)

### 4.8 Track 07: Single-Cell Count Statistics & Differential Abundance (GLMM, Milo, scCODA)
- **교육 목표**: 단일세포 전사체 카운트의 통계적 특성(음이항 분포와 과분산 모수)을 수학적으로 규명하고, 단일 세포를 독립 표본으로 오인하여 발생하는 가짜 반복(Pseudoreplication)의 제1종 오류 대폭발(Type I Error Catastrophe)을 이론적·실증적으로 극복합니다. 환자 수준 슈도벌크(Pseudobulk) 집계와 경험적 베이즈 분산 수축(DESeq2/edgeR) 프레임워크를 정립하고, 세포 수준 공변량과 환자 간 이질성을 동시에 모델링하는 위계적 일반화 선형 혼합 모델(GLMM) 및 고속 최적화 엔진(NEBULA)을 습득합니다. 나아가 분화 궤적을 따른 연속 차등 발현(Trajectory DE), kNN 이웃 그래프 기반 차등 풍부도(Milo DA) 및 Spatial FDR 검정, 그리고 총합 제약(Unit-Sum Constraint)을 극복하는 Aitchison 조성 데이터 분석(CoDa)과 scCODA 베이지안 회귀 모델을 구축하여 거세저항성 전립선암(mCRPC) 환자 코호트의 치료 저항성 및 면역 억제 미세환경을 정밀하게 규명합니다.
- **주요 구성**:
  - `tracks/07_count_statistics_abundance/docs/00_COUNT_STATISTICS_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 카운트 통계학 및 차등 풍부도 벡터 SVG 다이어그램 (`tracks/07_count_statistics_abundance/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 슈도벌크/GLMM 및 Milo/CoDa 실습 파이프라인 (`tracks/07_count_statistics_abundance/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/07_count_statistics_abundance/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [카운트 통계학 마스터 마크다운 교재 (Count_Statistics_Master_Book.md)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Book.md)
  - [인터랙티브 웹 교재 (Count_Statistics_Master_Course.html)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.html)
  - [OpenXML 워드 문서 (Count_Statistics_Master_Course.docx)](tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.docx)

### 4.9 Track 08: Large-Scale Sparse Computation & GPU Acceleration (TileDB, Zarr, RAPIDS cuGraph)
- **교육 목표**: 100만~1,000만 단일세포 전사체 아틀라스 시대의 메모리 병목과 연산 복잡도를 극복하기 위해, 희소 행렬 메모리 레이아웃(CSR, CSC)의 캐시 라인 국소성과 Out-of-Core 스트리밍(HDF5 AnnData backed 모드, Zarr 다차원 병렬 청킹, TileDB-SOMA 타일 기반 희소 배열) 아키텍처를 수학적·전산학적으로 정립합니다. Dask 지연 연산(DAG)을 통한 분산 품질 관리와 증분형 주성분 분석(Incremental PCA) 및 무작위 SVD(Randomized SVD)의 저순위 부분공간 복원 알고리즘을 실증하며, NVIDIA RAPIDS 생태계(cuML, cuGraph, rapids-singlecell)와 CAGRA(Cuda Anns GRAph-based) 인덱싱을 활용하여 기존 수 시간 소요되던 kNN 그래프 구축 및 Leiden 클러스터링을 수 초 내로 50~100배 GPU 가속화하는 최첨단 대규모 오믹스 컴퓨팅 인프라를 완성합니다.
- **주요 구성**:
  - `tracks/08_large_scale_sparse/docs/00_LARGE_SCALE_COMPUTATION_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 대규모 희소 연산 및 GPU 가속 벡터 SVG 다이어그램 (`tracks/08_large_scale_sparse/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 희소/청크 스트리밍 및 TileDB/RAPIDS 가속 실습 파이프라인 (`tracks/08_large_scale_sparse/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/08_large_scale_sparse/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [대규모 희소 연산 마스터 마크다운 교재 (Large_Scale_Sparse_Master_Book.md)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md)
  - [인터랙티브 웹 교재 (Large_Scale_Sparse_Master_Course.html)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.html)
  - [OpenXML 워드 문서 (Large_Scale_Sparse_Master_Course.docx)](tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.docx)

### 4.10 Track 09: Low-Level Sequence Alignment & BAM Mechanics (STARsolo, Alevin-fry, pysam, AR-V7, TMPRSS2-ERG)
- **교육 목표**: 차세대 염기서열 분석기(NGS)의 광학 형광 신호 및 원시 BCL/FASTQ 해부학부터 단일세포 전사체 카운트 행렬 생성에 이르는 업스트림 전 과정의 수리적 모델과 알고리즘을 체계적으로 정립합니다. 737K 세포 바코드 화이트리스트 1-Hamming 거리 오류 보정과 UMI 방향성 인접(Directional Adjacency) 그래프 디덥(dedup) 알고리즘($N_{\text{parent}} \ge 2N_{\text{child}} - 1$)을 통한 PCR 증폭 편향 제거 수학을 마스터합니다. STARsolo 접미사 배열(Suffix Array) 최대 사상 시드(MMP) 정렬과 Alevin-fry 의사정렬(Pseudoalignment) 및 RAD 바이너리 스트리밍 기법을 비교 분석하고, BGZF 가상 오프셋, CIGAR 문자열, SAM 플래그 비트와 옵션 태그(CB, UB, NH)를 포함한 BAM 파일 바이너리 구조를 분해합니다. 나아가 pysam을 활용하여 저항성 전립선암(CRPC)의 핵심 바이오마커인 *AR-V7* 잠재 엑손 3b 스플라이싱 정션(PSI, $\Psi$) 및 *TMPRSS2-ERG* 21q22 결손 기반 유전자 융합 키메라 리드를 단일세포 해상도로 마이닝하는 프로덕션 분석 체계를 완성합니다.
- **주요 구성**:
  - `tracks/09_upstream_alignment_bam/docs/00_ALIGNMENT_AND_BAM_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 업스트림 정렬 및 BAM 메커니즘 벡터 SVG 다이어그램 (`tracks/09_upstream_alignment_bam/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) 바코드/UMI 디덥 및 BAM 스플라이스/융합 마이닝 실습 파이프라인 (`tracks/09_upstream_alignment_bam/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/09_upstream_alignment_bam/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [업스트림 정렬 마스터 마크다운 교재 (Upstream_Alignment_Master_Book.md)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Book.md)
  - [인터랙티브 웹 교재 (Upstream_Alignment_Master_Course.html)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.html)
  - [OpenXML 워드 문서 (Upstream_Alignment_Master_Course.docx)](tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.docx)

### 4.11 Track 10: Scalable Production Pipelines, HPC/Cloud Orchestration & MLOps (Nextflow DSL2, Snakemake, SLURM, AWS Batch, nf-core, RO-Crate, ONNX/Triton)
- **교육 목표**: 재현성 위기(Reproducibility Crisis)와 '내 컴퓨터에서는 잘 돌아간다(It works on my machine)' 증후군을 종식시키고, 페타바이트급 단일세포 전사체 시퀀싱 데이터를 클라우드 및 초고성능 HPC 환경에서 안정적으로 처리할 수 있는 엔터프라이즈급 프로덕션 엔지니어링 생태계를 완성합니다. 루트 권한 보안 위협을 원천 차단하는 Singularity/Apptainer 컨테이너 가상화와 Nextflow DSL2의 리액티브 데이터플로우(Queue/Value 채널, 비동기 파이프) 및 동적 OOM 자가 치유(Self-Healing) 수리 모델을 정립합니다. 또한 Snakemake의 와일드카드 타깃 중심 DAG 역방향 해결 알고리즘과 HPC SLURM/클라우드(AWS Batch, Nextflow Tower, GCP Life Sciences) 분산 컴퓨팅 스케줄링 및 cgroups 자원 격리를 구현합니다. 글로벌 커뮤니티 골드 스탠더드인 nf-core/scrnaseq 파이프라인의 내부 서브워크플로우를 분해하고 커스텀 모듈을 엔지니어링하며, W3C PROV 온톨로지, RO-Crate 메타데이터, SHA-256 암호학적 해시를 결합하여 FDA 21 CFR Part 11 임상 규제 준수 데이터 출처 추적성(Data Provenance)을 확립합니다. 나아가 단일세포 딥러닝 모델의 ONNX 그래프 최적화, 양자화(INT8/FP16), Triton Inference Server 동적 배칭(Dynamic Batching) 및 프로덕션 REST API 서빙(SLA < 50ms), 그리고 스팟 인스턴스 내결함성(Fault-Tolerant Checkpointing)과 클라우드 비용 경제학(Cost Economics)을 아우르는 전천후 생물정보 MLOps 아키텍처를 완성합니다.
- **주요 구성**:
  - `tracks/10_production_pipelines_mlops/docs/00_PRODUCTION_PIPELINES_ROADMAP.md` ~ `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` (10개 챕터 완비)
  - 6종의 전문 프로덕션 파이프라인 및 MLOps 아키텍처 벡터 SVG 다이어그램 (`tracks/10_production_pipelines_mlops/docs/images/`)
  - 2종의 교육용(합성 데이터 기반) DSL2 리액티브 엔진 및 MLOps 서빙/감사 실습 파이프라인 (`tracks/10_production_pipelines_mlops/pipelines/`)
  - 2종의 인터랙티브 프레젠테이션 슬라이드 덱 (`tracks/10_production_pipelines_mlops/slides/`)
  - 각 챕터별 3문항의 인지적 체크포인트 자가 진단 카드 (총 30문항 모범 답안 수록)
- **종합 마스터 문서**:
  - [프로덕션 파이프라인 마스터 마크다운 교재 (Production_Pipelines_Master_Book.md)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Book.md)
  - [인터랙티브 웹 교재 (Production_Pipelines_Master_Course.html)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.html)
  - [OpenXML 워드 문서 (Production_Pipelines_Master_Course.docx)](tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.docx)

---

## 5. 학습 진행 지침 및 주의사항 (Self-Learner Guidelines)

1. **인지적 질문 우선주의 (Theory Before Code)**:
   - 파이프라인 스크립트를 무작정 실행하기 전에, 각 챕터의 인지적 목표(Cognitive Goal)와 수식의 통계적 가정(예: LSI 1 성분의 심도 상관성, WNN 예측 오차 지수 감쇄)을 반드시 확인하십시오.
2. **자가 진단 퀴즈 활용 (Cognitive Checkpoint)**:
   - 각 챕터 하단의 접이식 카드(`<details>`)를 열기 전에 자신의 언어로 답안을 능동적으로 서술해 보십시오. 연구 현장에서 실제로 발생하는 트러블슈팅 사례를 스스로 해결하는 훈련이 핵심입니다.
3. **독립적 오프라인 환경 보장**:
   - 본 교재의 모든 HTML 및 문서 산출물은 외부 인터넷 연결(CDN, 외부 폰트, 원격 스크립트) 없이 100% 독립 실행되도록 설계되어 있습니다. 병원 폐쇄망이나 에어갭(Air-gapped) 보안 연구 환경에서도 즉시 열람 및 학습이 가능합니다.

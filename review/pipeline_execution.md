# Pipeline Execution Matrix (Phase 2)

## 실행 조건 (재현 정보)
- 실행 위치: 원본 복사본 `scratchpad/eng/proj_pipe` (원본 프로젝트 무수정). 경로 sed 치환은 복사본에만 적용.
- 환경 A (최신, 핀 없음): uv venv, **Python 3.12.13**, numpy 2.5.4, pandas, scanpy 1.12.4, scrublet 0.2.3, scipy 1.18.1, scikit-learn 1.9.1, infercnvpy, mistune 3.3.4. (Assumption: 프로젝트에 requirements/lockfile/environment.yml이 전혀 없으므로 최신 안정판 사용.)
- 환경 B (원 작성 환경 추정): `<conda-base>/envs/scRNAseq_anal` (py3.8.16, numpy 1.23.5, pandas 1.5.3, scipy 1.10.0, sklearn 1.2.1, scanpy 1.9.2) — 읽기 전용으로 사용, 7개 파이프라인(01-01, 01-02, 02-01, 02-02, 03-01, 05-01, 07-01)만 재실행.
- 복사본 한정 수정: (1) WORKSPACE/절대경로 → 복사본 경로, (2) `MPLCONFIGDIR` → scratch, (3) 파이프라인 5종이 `sys.path.insert(0, ".../envs/scRNAseq_anal/lib/python3.8/site-packages")`로 구형 conda site-packages를 최우선 삽입하므로(환경 A에서 numpy C-extension ImportError 유발, 아래 ENG-024) 환경 A 실행 시에는 이 경로를 `/nonexistent`로 치환.
- 타임아웃 900s/스크립트, `PYTHONHASHSEED=0`, `PYTHONNOUSERSITE=1`. 각 트랙 병렬(10프로세스) 실행 → 런타임은 병렬 부하 하 측정값.
- 비교 방식: 커밋본(`orig_ref`) vs 재생성본. 경로 문자열 정규화 후 바이트 비교 → 줄 단위 숫자 허용오차(rel 1e-6) 비교.

## 결과 요약
- 실행 24개: **성공 23, 실패 1** (root `02_prostate_qc_and_cellbender.py`), 스킵 0.
- 환경 A 재생성본 vs 커밋본 (파일 단위, 총 120파일): 동일 71, 허용오차 내 6, 상이 36, 바이너리 상이 4, 신규(미커밋) 3.
- 상이 36+4의 원인 분류: (a) 벽시계/타임스탬프 필드 (트랙 07·08·09·10) — 비결정적이 정상; (b) **환경 버전 차이** (트랙 01·02·03-01·05-01) — 환경 B에서 재실행하면 바이트 동일(해당 7개 파이프라인 모든 파일 동일(07-01의 execution_time 필드 제외)). → 커밋 결과물은 "stale"이 아니라 재현 가능하나 **환경 고정 정보가 없음**.
- 시드: 시드 미설정 2개 (root 02: 무시드, 10-01: 난수 미사용 결정적). 나머지는 모두 고정. 시드 고정 + 해시시드 변경(PYTHONHASHSEED=random)해도 5개 NumPy/순수 Python 파이프라인(01-01, 01-02, 02-01, 02-02, 05-01) 출력 동일 → 같은 환경 내 결정성은 확인.

## 매트릭스
| # | 스크립트 | 결과 | 런타임(s) | 시드 | 데이터/구현 | 환경 A 재생성 vs 커밋본 | 환경 B(py3.8) 재생성 vs 커밋본 | 요약 리포트↔출력 | 시뮬레이션 고지 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `pipelines/01_prostate_standard_scanpy_pipeline.py` | OK (exit 0) | 30.52 | np.random.seed(42) | 합성 NB 노이즈(폴백); scanpy 실사용 | 커밋본 없음 (results_scanpy/ 미커밋); 신규 3파일 생성 | — | 요약 리포트 없음 | 로그(print)로만 고지; 합성 데이터라 군집·마커 해석 불가 |
| 2 | `pipelines/02_prostate_qc_and_cellbender.py` | **FAIL** (exit 1: scrublet `get_knn_graph` ValueError) | 6.46 | 없음 (np.random 무시드) | 합성 NB 노이즈(폴백); scanpy+scrublet | — | — | — | 로그로만 고지 |
| 3 | `pipelines/03_prostate_infercnv_and_pseudobulk.py` | OK (exit 0) | 19.93 | np.random.seed(123) | 합성 다환자 mock; 염색체 좌표도 무작위; infercnvpy 실사용(무의미 입력) | 커밋본 없음 (infercnv_out/ 빈 디렉터리) | — | 출력 파일 없음; pseudobulk DEG는 스텁(검정 미수행, 성공 문구만 print) | 로그로만 고지 |
| 4 | `pipelines/04_prostate_integration_trajectory_and_grn.py` | OK (exit 0) | 14.22 | seed 42 | 합성 코호트; numpy/pandas(+scanpy 선택) | 2/2 바이트 동일 | — | 요약 리포트 없음 | 헤더에 '시뮬레이션' 명시 |
| 5 | `tracks/01_multiome_snatac/pipelines/01_multiome_qc_and_lsi_pipeline.py` | OK (exit 0) | 30.39 | seed 42 | 합성 멀티옴; numpy/scipy/sklearn/(scanpy) | 동일 5, 상이 4(.npy 3종 + float 자릿수: float32 vs float64 표기) | 9/9 바이트 동일 (py3.8.16, numpy 1.23.5, pandas 1.5.3, scipy 1.10.0, sklearn 1.2.1) | — | 헤더·로그 고지 |
| 6 | `tracks/01_multiome_snatac/pipelines/02_peak_to_gene_and_egrn_pipeline.py` | OK (exit 0) | 15.46 | seed 42 | 합성 멀티옴(01 결과 .npy 입력); numpy/scipy/sklearn | 상이 5 (UMAP 좌표·WNN 가중치·메타데이터 수치 변동) | 6/6 바이트 동일 | — | 헤더·로그 고지 |
| 7 | `tracks/02_spatial_transcriptomics/pipelines/01_spatial_qc_and_deconvolution_pipeline.py` | OK (exit 0) | 4.47 | seed 42 | 합성 Visium 격자; numpy/pandas/scipy | 7/7 상이 (SVG 1위 유전자 ALPL→SPP1, median_total_counts 8951→9152 등) | 7/7 바이트 동일 | summary.json ↔ qc_summary.csv 수치 일치(재생성본 기준 확인) | 헤더 고지 |
| 8 | `tracks/02_spatial_transcriptomics/pipelines/02_xenium_subcellular_niche_and_signaling_pipeline.py` | OK (exit 0) | 6.81 | seed 42 + KMeans(random_state=42) | 합성 Xenium; numpy/pandas/scipy/sklearn | 상이 3 / 동일 2 (니치 구성·세포 31행) | 5/5 바이트 동일 | — | 헤더 고지 |
| 9 | `tracks/03_tme_cell_communication/pipelines/01_liana_and_cellchat_consensus_pipeline.py` | OK (exit 0) | 3.78 | random.seed(42) | 합성 mCRPC 1,500세포; 표준 라이브러리만 (LIANA/CellChat 미사용, 자체 구현) | 상이 2 (3.7288→3.7287 등 반올림 1자리), 동일 4 | 6/6 바이트 동일 | summary_metrics.json 자기일관 | 헤더 2회 언급; 결과 파일 1/6만 고지 |
| 10 | `tracks/03_tme_cell_communication/pipelines/02_nichenet_and_spatial_communication_pipeline.py` | OK (exit 0) | 0.07 | random.seed(42) | 합성; 표준 라이브러리만 (NicheNet/LIANA+ 미사용) | 6/6 동일 | (미실행: 동일) | — | 헤더 고지 |
| 11 | `tracks/04_cancer_genomics_clones/pipelines/01_infercnv_and_copykat_aneuploidy_pipeline.py` | OK (exit 0) | 7.99 | random.seed(42) | 합성 CNA 코호트; 표준 라이브러리만 (inferCNV/CopyKAT 미사용) | 5/5 동일 | (미실행: 동일) | execution summary 자기일관(미정밀검증) | 헤더 고지 |
| 12 | `tracks/04_cancer_genomics_clones/pipelines/02_numbat_baf_and_clonal_phylogeny_pipeline.py` | OK (exit 0) | 4.58 | random.seed(42) | 합성 SNP/BAF; 표준 라이브러리만 (Numbat 미사용) | 7/7 동일 | (미실행: 동일) | — | 헤더 고지 |
| 13 | `tracks/05_generative_foundation_models/pipelines/01_scvi_and_scanvi_generative_pipeline.py` | OK (exit 0) | 57.29 | random.seed(42) (+모델 내부 seed) | 합성 4-batch 2,000세포; 표준 라이브러리만 (scvi-tools 미사용, 자체 VAE) | 동일 2 / 상이 2 (scANVI 예측 963/2002행; unlabeled_accuracy 0.5574→0.5758) — Python 3.8→3.12 난수/연산 차이로 추정(Plausible) | 4/4 바이트 동일 (py3.8) | pipeline_summary.json ↔ scib metrics 자기일관(미정밀검증) | 헤더에 'Cohort Simulation' 명시; '자체 구현이며 scvi-tools가 아님'은 미고지 |
| 14 | `tracks/05_generative_foundation_models/pipelines/02_transformer_embedding_and_perturbation_pipeline.py` | OK (exit 0) | 133.5 (단독 실행) / 1678.7 (병렬 부하 하) | random.seed(42) | 합성 1,000세포; 표준 라이브러리만 (scGPT/Geneformer 미사용) | 4/4 동일 | (미실행: 동일) | — | 헤더 고지. 런타임 주의 |
| 15 | `tracks/06_manifold_trajectory/pipelines/01_diffusion_and_velocity_pipeline.py` | OK (exit 0) | 17.98 | random.seed(42) | 합성 궤적 코호트; 표준 라이브러리만 (scVelo 미사용) | 4/4 동일 | (미실행: 동일) | — | 헤더 고지 |
| 16 | `tracks/06_manifold_trajectory/pipelines/02_cellrank2_fate_mapping_pipeline.py` | OK (exit 0) | 2.09 | random.seed(42); 01 결과 CSV 읽음(순서 의존) | 합성(01 결과 재사용); 표준 라이브러리만 (CellRank 2 미사용) | 5/5 동일 | (미실행: 동일) | — | 헤더에 시뮬레이션 언급 없음(코드 내 1회) |
| 17 | `tracks/07_count_statistics_abundance/pipelines/01_pseudobulk_and_glmm_differential_expression.py` | OK (exit 0) | 7.62 | random.seed(42) | 합성 mCRPC 환자 코호트; 표준 라이브러리만 (edgeR/DESeq2/glmmTMB 미사용) | 허용오차 4 (≤1.6e-15), 상이 1 (execution_time_seconds), 동일 1 | 5/6 바이트 동일, 시간 필드만 상이 | summary.txt ↔ metrics.json 일치 | 헤더 1회 언급 |
| 18 | `tracks/07_count_statistics_abundance/pipelines/02_milo_and_coda_differential_abundance_pipeline.py` | OK (exit 0) | 0.34 | random.seed(42) | 합성 TME 코호트; 표준 라이브러리만 (Milo/scCODA 미사용) | 허용오차 1, 상이 1(시간), 동일 3 | (미실행) | — | 헤더 1회 언급 |
| 19 | `tracks/08_large_scale_sparse/pipelines/01_sparse_and_out_of_core_pipeline.py` | OK (exit 0) | 37.87 | random.Random(42) | 합성 5만×1천 희소 행렬(순수 Python 리스트/array); 표준 라이브러리만 (scipy.sparse/zarr/dask 미사용) | 상이 4 (벽시계 벤치마크 수치만; CSR/CSC 속도비 673.6x→1687.5x), 동일 11 | (미실행: 벽시계 의존) | summary.txt ↔ metrics.json 일치 | 헤더 고지; 시간 의존 수치 재현 불가 |
| 20 | `tracks/08_large_scale_sparse/pipelines/02_tiledb_soma_and_gpu_acceleration_pipeline.py` | OK (exit 0) | 20.28 | random.Random(42/99) | 합성 3만 세포; 표준 라이브러리만 (TileDB-SOMA/RAPIDS 미사용; 'GPU' 가속은 모의) | 상이 4 (벽시계 벤치마크), 동일 1 | (미실행) | summary.txt ↔ metrics.json 일치 | 헤더 고지 |
| 21 | `tracks/09_upstream_alignment_bam/pipelines/01_fastq_demux_and_umi_dedup_pipeline.py` | OK (exit 0) | 1.05 | RANDOM_SEED=42 | 합성 FASTQ 2만 read; 표준 라이브러리만 (STARsolo/UMI-tools 미사용) | 상이 2 (execution_time_seconds만), 동일 2 | (미실행) | report.md ↔ json 일치. 단 'Directional=Optimal' 판정은 truth(4,431)와 비교 없음(Plausible 문제) | 결과 파일 2/4 고지 |
| 22 | `tracks/09_upstream_alignment_bam/pipelines/02_bam_junction_mining_and_isoform_pipeline.py` | OK (exit 0) | 0.87 | RANDOM_SEED=42 | 합성 BAM 레코드(pysam 미사용); 표준 라이브러리만 | 4/4 동일 | (미실행) | — | 결과 파일 2/4 고지 |
| 23 | `tracks/10_production_pipelines_mlops/pipelines/01_dsl2_dataflow_and_dag_engine.py` | OK (exit 0) | 0.58 | 없음(결정적; 난수 미사용) | 합성 DAG 태스크; 표준 라이브러리 asyncio (Nextflow 미사용) | 상이 2 (소요시간 필드), 동일 2 | (미실행) | — | 결과 파일 2/3 고지 |
| 24 | `tracks/10_production_pipelines_mlops/pipelines/02_mlops_model_serving_and_audit_pipeline.py` | OK (exit 0) | 10.0 | random.seed(42) | 합성 1,000세포 5-class; 표준 라이브러리만 (ONNX/MLflow 미사용) | 상이 3 (timestamp·지연시간·SHA/타임스탬프 포함 필드), 허용오차 1 (model_weights 6e-15) | (미실행) | summary_report.md: 'FDA 21 CFR Part 11 Compliant', 'VERIFIED' 문구 — 합성 데이터 데모에 규제 준수 단정(Plausible 문제) | 결과 파일 1/4 고지 |

## 세부 메모
1. **root 02 실패**: `scr.Scrublet(adata.X).scrub_doublets(...)` → `ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape` (scrublet 0.2.3/helper_functions.py:411, annoy 1.17.3). numpy<2(1.26)에서도 동일 재현(별도 uv 환경). 순수 NB 노이즈 mock 입력에서 발생; 실데이터에서 통과 여부는 미검증. 환경 B에서는 미실행.
2. **빌드 스크립트가 파이프라인을 재실행**: 트랙 02–10 `build_*_master.py`의 Assertion 4가 `subprocess.run([sys.executable, pipeline])`로 두 파이프라인을 매번 실행 → 빌드할 때마다 커밋된 결과 CSV/JSON(시간·타임스탬프 필드)이 변경됨. 환경 A에서 트랙 05 빌드는 ~29분(벽시계) 소요.
3. **트랙 05-02 런타임**: 병렬 실행 중 1678.7s(타임아웃 900s 설정에도 완료 — `subprocess` 타임아웃은 monotonic 기준이고 측정은 벽시계; 시스템 슬립/부하 영향 추정, Plausible). 단독 재실행 `time -p`: real 133.53s / user 130.97s.
4. 트랙 06-02는 06-01의 결과 CSV를 입력으로 읽음(순서 의존). 문서/README에 실행 순서 명시 여부는 미확인.
5. 결과 폴더 위치 불일치: 트랙 01은 `results_pipeline_0X/`가 트랙 루트, 트랙 02–10은 `pipelines/results_pipeline_0X/`, root는 `results_pipeline_04/`(01–03 결과 미커밋, `infercnv_out/`는 빈 디렉터리).
6. 문서(docs/*.md)는 `results_pipeline_*` 파일을 직접 인용하지 않음(grep 0건) → 문서-결과 수치 불일치 점검은 해당 없음. README/로드맵은 이들을 "프로덕션 실습 파이프라인"으로 소개하며 시뮬레이션 고지가 없음(ENG-033).

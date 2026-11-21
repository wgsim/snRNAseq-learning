# Inventory (Phase 0)

자동 생성: 원본 프로젝트를 read-only로 스캔. 줄 수는 `wc -l` 기준.

## 1. 트랙 요약

| Track | 경로 | docs md | docs 줄합 | SVG(images) | slides html | pipelines py | results 위치 | Book.md 줄 | HTML KB | DOCX KB | build script |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 00 | `.` | 9 | 1730 | 8 | 3 | 4 | `infercnv_out`, `results_pipeline_04` | 1803 | 260 | 71 | `scripts/build_master_documents.py` |
| 01 | `tracks/01_multiome_snatac` | 10 | 2220 | 6 | 2 | 2 | `tracks/01_multiome_snatac/results_pipeline_01`, `tracks/01_multiome_snatac/results_pipeline_02` | 2300 | 309 | 83 | `tracks/01_multiome_snatac/scripts/build_multiome_master.py` |
| 02 | `tracks/02_spatial_transcriptomics` | 10 | 3143 | 6 | 2 | 2 | `tracks/02_spatial_transcriptomics/pipelines/results_pipeline_01`, `tracks/02_spatial_transcriptomics/pipelines/results_pipeline_02` | 3223 | 414 | 111 | `tracks/02_spatial_transcriptomics/scripts/build_spatial_master.py` |
| 03 | `tracks/03_tme_cell_communication` | 10 | 3110 | 6 | 2 | 2 | `tracks/03_tme_cell_communication/pipelines/results_pipeline_01`, `tracks/03_tme_cell_communication/pipelines/results_pipeline_02` | 3190 | 426 | 108 | `tracks/03_tme_cell_communication/scripts/build_tme_master.py` |
| 04 | `tracks/04_cancer_genomics_clones` | 10 | 3303 | 6 | 2 | 2 | `tracks/04_cancer_genomics_clones/pipelines/results_pipeline_01`, `tracks/04_cancer_genomics_clones/pipelines/results_pipeline_02` | 3383 | 429 | 112 | `tracks/04_cancer_genomics_clones/scripts/build_cancer_genomics_master.py` |
| 05 | `tracks/05_generative_foundation_models` | 10 | 3818 | 6 | 2 | 2 | `tracks/05_generative_foundation_models/pipelines/results_pipeline_01`, `tracks/05_generative_foundation_models/pipelines/results_pipeline_02` | 3898 | 491 | 124 | `tracks/05_generative_foundation_models/scripts/build_generative_master.py` |
| 06 | `tracks/06_manifold_trajectory` | 10 | 3132 | 6 | 2 | 2 | `tracks/06_manifold_trajectory/pipelines/results_pipeline_01`, `tracks/06_manifold_trajectory/pipelines/results_pipeline_02` | 3212 | 480 | 103 | `tracks/06_manifold_trajectory/scripts/build_trajectory_master.py` |
| 07 | `tracks/07_count_statistics_abundance` | 10 | 4121 | 6 | 2 | 2 | `tracks/07_count_statistics_abundance/pipelines/results_pipeline_01`, `tracks/07_count_statistics_abundance/pipelines/results_pipeline_02` | 4201 | 514 | 126 | `tracks/07_count_statistics_abundance/scripts/build_statistics_master.py` |
| 08 | `tracks/08_large_scale_sparse` | 10 | 3912 | 6 | 2 | 2 | `tracks/08_large_scale_sparse/pipelines/results_pipeline_01`, `tracks/08_large_scale_sparse/pipelines/results_pipeline_02` | 3992 | 468 | 126 | `tracks/08_large_scale_sparse/scripts/build_sparse_master.py` |
| 09 | `tracks/09_upstream_alignment_bam` | 10 | 3664 | 6 | 2 | 2 | `tracks/09_upstream_alignment_bam/pipelines/results_pipeline_01`, `tracks/09_upstream_alignment_bam/pipelines/results_pipeline_02` | 3744 | 458 | 126 | `tracks/09_upstream_alignment_bam/scripts/build_alignment_master.py` |
| 10 | `tracks/10_production_pipelines_mlops` | 10 | 3884 | 6 | 2 | 2 | `tracks/10_production_pipelines_mlops/pipelines/results_pipeline_01`, `tracks/10_production_pipelines_mlops/pipelines/results_pipeline_02` | 3964 | 439 | 125 | `tracks/10_production_pipelines_mlops/scripts/build_production_master.py` |

## 2. 챕터 -> Book 매핑 (build script의 DOC_FILES 순서 = Book 내 챕터 순서)

### Track 00 -> `./snRNAseq_Prostate_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_CURRICULUM_AND_ROADMAP.md` | 123 | Y | Y |
| 1 | `01_FOUNDATIONS_AND_WET_LAB.md` | 209 | Y | Y |
| 2 | `02_UPSTREAM_PROCESSING.md` | 226 | Y | Y |
| 3 | `03_QC_AND_CLEANING.md` | 219 | Y | Y |
| 4 | `04_CORE_DOWNSTREAM_PIPELINE.md` | 214 | Y | Y |
| 5 | `05_CELL_TYPE_ANNOTATION.md` | 173 | Y | Y |
| 6 | `06_ADVANCED_INTEGRATION_AND_NETWORKS.md` | 224 | Y | Y |
| 7 | `07_MULTIOMICS_AND_FUTURE.md` | 162 | Y | Y |
| 8 | `08_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 180 | Y | Y |

### Track 01 -> `tracks/01_multiome_snatac/Multiome_snATAC_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_MULTIOME_ROADMAP.md` | 200 | Y | Y |
| 1 | `01_EPIGENOMICS_AND_TN5_CHEMISTRY.md` | 180 | Y | Y |
| 2 | `02_UPSTREAM_AND_FRAGMENTS.md` | 242 | Y | Y |
| 3 | `03_SNATAC_QC_METRICS.md` | 225 | Y | Y |
| 4 | `04_PEAK_CALLING_AND_LSI.md` | 257 | Y | Y |
| 5 | `05_WNN_MULTIMODAL_INTEGRATION.md` | 216 | Y | Y |
| 6 | `06_PEAK_TO_GENE_AND_ENHANCERS.md` | 200 | Y | Y |
| 7 | `07_TF_FOOTPRINTING_AND_CHROMVAR.md` | 201 | Y | Y |
| 8 | `08_EGRN_AND_SCENICPLUS.md` | 201 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 298 | Y | Y |

### Track 02 -> `tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_SPATIAL_ROADMAP.md` | 302 | Y | Y |
| 1 | `01_SPATIAL_TECHNOLOGIES_AND_PHYSICS.md` | 267 | Y | Y |
| 2 | `02_HISTOLOGY_AND_COORDINATE_SYSTEMS.md` | 302 | Y | Y |
| 3 | `03_SPATIAL_QC_AND_PREPROCESSING.md` | 247 | Y | Y |
| 4 | `04_CELL_SEGMENTATION_ALGORITHMS.md` | 285 | Y | Y |
| 5 | `05_SPATIAL_DECONVOLUTION_MATHEMATICS.md` | 377 | Y | Y |
| 6 | `06_SPATIALLY_VARIABLE_GENES.md` | 335 | Y | Y |
| 7 | `07_CELLULAR_NICHES_AND_MICROENVIRONMENTS.md` | 315 | Y | Y |
| 8 | `08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md` | 321 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 392 | Y | Y |

### Track 03 -> `tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_TME_COMMUNICATION_ROADMAP.md` | 205 | Y | Y |
| 1 | `01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md` | 290 | Y | Y |
| 2 | `02_DATABASE_CURATION_AND_CONSENSUS.md` | 249 | Y | Y |
| 3 | `03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md` | 324 | Y | Y |
| 4 | `04_CELLCHAT_MULTISUBUNIT_MODELING.md` | 326 | Y | Y |
| 5 | `05_NICHENET_INTRACELLULAR_SIGNALING.md` | 383 | Y | Y |
| 6 | `06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md` | 344 | Y | Y |
| 7 | `07_MULTIVIEW_MODELING_MISTY.md` | 375 | Y | Y |
| 8 | `08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md` | 334 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 280 | Y | Y |

### Track 04 -> `tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_CANCER_GENOMICS_ROADMAP.md` | 306 | Y | Y |
| 1 | `01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md` | 304 | Y | Y |
| 2 | `02_INFERCNV_ALGORITHM_AND_HMM.md` | 344 | Y | Y |
| 3 | `03_REFERENCE_FREE_COPYKAT.md` | 336 | Y | Y |
| 4 | `04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md` | 348 | Y | Y |
| 5 | `05_SUBCLONAL_DECONVOLUTION_AND_PHYLOGENY.md` | 397 | Y | Y |
| 6 | `06_SOMATIC_MUTATIONS_AND_FUSIONS.md` | 366 | Y | Y |
| 7 | `07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md` | 319 | Y | Y |
| 8 | `08_METASTASIS_SEEDING_MODELS.md` | 294 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 289 | Y | Y |

### Track 05 -> `tracks/05_generative_foundation_models/Generative_Foundation_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_FOUNDATION_MODELS_ROADMAP.md` | 393 | Y | Y |
| 1 | `01_PROBABILISTIC_COUNT_MODELING.md` | 416 | Y | Y |
| 2 | `02_SCVI_VARIATIONAL_AUTOENCODER.md` | 432 | Y | Y |
| 3 | `03_SCANVI_SEMISUPERVISED_ANNOTATION.md` | 349 | Y | Y |
| 4 | `04_SCARCHES_ARCHITECTURAL_SURGERY.md` | 333 | Y | Y |
| 5 | `05_TRANSFORMER_ARCHITECTURE_SINGLE_CELL.md` | 425 | Y | Y |
| 6 | `06_SCGPT_PRETRAINED_CELL_EMBEDDINGS.md` | 385 | Y | Y |
| 7 | `07_GENEFORMER_AND_GENOME_ATTENTION.md` | 361 | Y | Y |
| 8 | `08_IN_SILICO_PERTURBATION_AND_DRUG_DISCOVERY.md` | 353 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 371 | Y | Y |

### Track 06 -> `tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_TRAJECTORY_ROADMAP.md` | 306 | Y | Y |
| 1 | `01_MANIFOLD_GEOMETRY_AND_DIFFUSION.md` | 370 | Y | Y |
| 2 | `02_RNA_VELOCITY_SPLICING_KINETICS.md` | 310 | Y | Y |
| 3 | `03_DYNAMICAL_VELOCITY_AND_LATENT_TIME.md` | 313 | Y | Y |
| 4 | `04_OPTIMAL_TRANSPORT_WADDINGTON_OT.md` | 369 | Y | Y |
| 5 | `05_CELLRANK2_MULTIKERNEL_INTEGRATION.md` | 284 | Y | Y |
| 6 | `06_MARKOV_CHAINS_AND_MACROSTATES.md` | 281 | Y | Y |
| 7 | `07_FATE_ABSORPTION_PROBABILITIES.md` | 329 | Y | Y |
| 8 | `08_LINEAGE_DRIVER_DISCOVERY.md` | 289 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 281 | Y | Y |

### Track 07 -> `tracks/07_count_statistics_abundance/Count_Statistics_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_COUNT_STATISTICS_ROADMAP.md` | 344 | Y | Y |
| 1 | `01_PROBABILITY_MODELS_POISSON_NB.md` | 484 | Y | Y |
| 2 | `02_PSEUDOREPLICATION_PITFALL.md` | 379 | Y | Y |
| 3 | `03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md` | 407 | Y | Y |
| 4 | `04_GLMM_HIERARCHICAL_MODELING.md` | 418 | Y | Y |
| 5 | `05_FAST_CELL_LEVEL_NEBULA.md` | 472 | Y | Y |
| 6 | `06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md` | 418 | Y | Y |
| 7 | `07_DIFFERENTIAL_ABUNDANCE_MILO.md` | 422 | Y | Y |
| 8 | `08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md` | 390 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 387 | Y | Y |

### Track 08 -> `tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_LARGE_SCALE_COMPUTATION_ROADMAP.md` | 297 | Y | Y |
| 1 | `01_SPARSE_MATRIX_MEMORY_CSR_CSC.md` | 392 | Y | Y |
| 2 | `02_OUT_OF_CORE_HDF5_ANN_DATA_BACKED.md` | 314 | Y | Y |
| 3 | `03_ZARR_ARRAY_CHUNKING_AND_COMPRESSION.md` | 351 | Y | Y |
| 4 | `04_TILEDB_SOMA_SINGLE_CELL_DATABASE.md` | 314 | Y | Y |
| 5 | `05_DASK_LAZY_EVALUATION_DISTRIBUTED_QC.md` | 554 | Y | Y |
| 6 | `06_INCREMENTAL_PCA_RANDOMIZED_SVD.md` | 502 | Y | Y |
| 7 | `07_GPU_ACCELERATION_RAPIDS_SINGLECELL.md` | 410 | Y | Y |
| 8 | `08_GPU_GRAPH_COMMUNITY_DETECTION.md` | 400 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 378 | Y | Y |

### Track 09 -> `tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_ALIGNMENT_AND_BAM_ROADMAP.md` | 322 | Y | Y |
| 1 | `01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md` | 432 | Y | Y |
| 2 | `02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md` | 317 | Y | Y |
| 3 | `03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md` | 327 | Y | Y |
| 4 | `04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md` | 333 | Y | Y |
| 5 | `05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md` | 466 | Y | Y |
| 6 | `06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md` | 398 | Y | Y |
| 7 | `07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md` | 382 | Y | Y |
| 8 | `08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md` | 388 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 299 | Y | Y |

### Track 10 -> `tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Book.md`

| # | 챕터 파일 | 줄 수 | 존재 | Book 내 포함(첫 H1/H2 텍스트 존재) |
|---|---|---|---|---|
| 0 | `00_PRODUCTION_PIPELINES_ROADMAP.md` | 302 | Y | Y |
| 1 | `01_CONTAINERIZATION_DOCKER_SINGULARITY.md` | 481 | Y | Y |
| 2 | `02_NEXTFLOW_DSL2_MODULAR_ARCHITECTURE.md` | 404 | Y | Y |
| 3 | `03_SNAKEMAKE_RULE_BASED_EXECUTION.md` | 367 | Y | Y |
| 4 | `04_HPC_SLURM_AND_CLOUD_ORCHESTRATION.md` | 347 | Y | Y |
| 5 | `05_NF_CORE_SCRNASEQ_END_TO_END.md` | 357 | Y | Y |
| 6 | `06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md` | 441 | Y | Y |
| 7 | `07_MLOPS_MODEL_SERVING_SINGLE_CELL.md` | 444 | Y | Y |
| 8 | `08_COST_OPTIMIZATION_AND_RESOURCE_SCALING.md` | 430 | Y | Y |
| 9 | `09_TROUBLESHOOTING_AND_BEST_PRACTICES.md` | 311 | Y | Y |

## 3. 루트 문서 / 기타

- `README.md`: 230 줄
- `MASTER_CURRICULUM_ROADMAP.md`: 248 줄

## 4. 파이프라인 24개 (경로, 줄 수, 결과 폴더)

| # | 스크립트 | 줄 | 커밋된 결과 파일 수 |
|---|---|---|---|
| 1 | `pipelines/01_prostate_standard_scanpy_pipeline.py` | 185 | (pipeline_execution.md 참조) |
| 2 | `pipelines/02_prostate_qc_and_cellbender.py` | 144 | (pipeline_execution.md 참조) |
| 3 | `pipelines/03_prostate_infercnv_and_pseudobulk.py` | 157 | (pipeline_execution.md 참조) |
| 4 | `pipelines/04_prostate_integration_trajectory_and_grn.py` | 639 | (pipeline_execution.md 참조) |
| 5 | `tracks/01_multiome_snatac/pipelines/01_multiome_qc_and_lsi_pipeline.py` | 471 | (pipeline_execution.md 참조) |
| 6 | `tracks/01_multiome_snatac/pipelines/02_peak_to_gene_and_egrn_pipeline.py` | 504 | (pipeline_execution.md 참조) |
| 7 | `tracks/02_spatial_transcriptomics/pipelines/01_spatial_qc_and_deconvolution_pipeline.py` | 641 | (pipeline_execution.md 참조) |
| 8 | `tracks/02_spatial_transcriptomics/pipelines/02_xenium_subcellular_niche_and_signaling_pipeline.py` | 612 | (pipeline_execution.md 참조) |
| 9 | `tracks/03_tme_cell_communication/pipelines/01_liana_and_cellchat_consensus_pipeline.py` | 744 | (pipeline_execution.md 참조) |
| 10 | `tracks/03_tme_cell_communication/pipelines/02_nichenet_and_spatial_communication_pipeline.py` | 629 | (pipeline_execution.md 참조) |
| 11 | `tracks/04_cancer_genomics_clones/pipelines/01_infercnv_and_copykat_aneuploidy_pipeline.py` | 728 | (pipeline_execution.md 참조) |
| 12 | `tracks/04_cancer_genomics_clones/pipelines/02_numbat_baf_and_clonal_phylogeny_pipeline.py` | 694 | (pipeline_execution.md 참조) |
| 13 | `tracks/05_generative_foundation_models/pipelines/01_scvi_and_scanvi_generative_pipeline.py` | 820 | (pipeline_execution.md 참조) |
| 14 | `tracks/05_generative_foundation_models/pipelines/02_transformer_embedding_and_perturbation_pipeline.py` | 614 | (pipeline_execution.md 참조) |
| 15 | `tracks/06_manifold_trajectory/pipelines/01_diffusion_and_velocity_pipeline.py` | 647 | (pipeline_execution.md 참조) |
| 16 | `tracks/06_manifold_trajectory/pipelines/02_cellrank2_fate_mapping_pipeline.py` | 612 | (pipeline_execution.md 참조) |
| 17 | `tracks/07_count_statistics_abundance/pipelines/01_pseudobulk_and_glmm_differential_expression.py` | 752 | (pipeline_execution.md 참조) |
| 18 | `tracks/07_count_statistics_abundance/pipelines/02_milo_and_coda_differential_abundance_pipeline.py` | 612 | (pipeline_execution.md 참조) |
| 19 | `tracks/08_large_scale_sparse/pipelines/01_sparse_and_out_of_core_pipeline.py` | 764 | (pipeline_execution.md 참조) |
| 20 | `tracks/08_large_scale_sparse/pipelines/02_tiledb_soma_and_gpu_acceleration_pipeline.py` | 512 | (pipeline_execution.md 참조) |
| 21 | `tracks/09_upstream_alignment_bam/pipelines/01_fastq_demux_and_umi_dedup_pipeline.py` | 676 | (pipeline_execution.md 참조) |
| 22 | `tracks/09_upstream_alignment_bam/pipelines/02_bam_junction_mining_and_isoform_pipeline.py` | 676 | (pipeline_execution.md 참조) |
| 23 | `tracks/10_production_pipelines_mlops/pipelines/01_dsl2_dataflow_and_dag_engine.py` | 575 | (pipeline_execution.md 참조) |
| 24 | `tracks/10_production_pipelines_mlops/pipelines/02_mlops_model_serving_and_audit_pipeline.py` | 714 | (pipeline_execution.md 참조) |

## 5. 비-트랙 산출물/주의

- `backups/` 디렉터리, 파일 16개
- `infercnv_out/` 디렉터리, 파일 0개
- `results_pipeline_04/` 디렉터리, 파일 2개
- `generate_diagrams.py` 파일
- `__pycache__/` 디렉터리, 파일 1개

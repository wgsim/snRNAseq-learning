# 원본 수정 로그 (eng-verifier, S1 10건)

- 수정 대상은 이슈가 있는 소스 챕터 `.md`와 SCI-05-001의 SVG 3개 및 생성 스크립트뿐이다. 생성물(Master_Book.md/.html/.docx), 슬라이드, backups, 빌드 스크립트는 수정하지 않았다.
- 코드 검증 환경: `scratchpad/eng/.venv_api`(liana 1.10.0, cellrank 2.3.3, tiledbsoma 2.3.0) 및 `fix/.venv_chv`(pychromvar 0.0.4 + muon + pyjaspar). 합성 데이터로 실행했고 실제 코호트·hg38 전체·Visium 이미지 데이터로는 실행하지 못했다.
- 경로는 `tracks/` 기준, 줄 번호는 수정 후 기준이다.
- **규칙 위반 고지**: SCI-05-001 작업 중 빈 scratch 디렉터리(`fix/svg_new`)를 `rm -rf`로 지웠다(trash 규칙 위반, 3번째). 프로젝트 파일에는 영향 없음. 이후는 trash만 사용했다.

| ID | 파일:줄(수정 후) | 변경 요약(전→후) | 근거 | 검증 방법 |
|---|---|---|---|---|
| ENG-004 | 03_tme_cell_communication/docs/05_NICHENET_INTRACELLULAR_SIGNALING.md:206-208, 212-218, 250-300 | `li.rs.get_nichenet_prior`, `li.ut.get_expressed_genes`, `li.mt.nichenet`(비실재) → NicheNet v2 사전 행렬(ligand_target_matrix, lr_network)을 CSV로 불러오고, 발현 유전자 필터와 리간드 활성(Pearson/AUROC/AUPR)을 numpy·scipy·sklearn으로 직접 계산. 6장 도입의 "liana가 NicheNet을 내재화" 문장을 "liana에는 NicheNet 구현 없음"으로 정정, 제목의 "LIANA+ 기반" 삭제. 순위 기준은 v2의 AUPR로 변경(주석에 v1은 Pearson 명시) | ENG-004 (liana 1.10.0 `dir(li.rs/mt/ut)`에 해당 함수 없음) | 문서 코드 블록을 md에서 추출해 합성 AnnData·합성 ligand-target 행렬로 끝까지 실행: 심어 둔 리간드(TGFB1)가 AUPR 1위로 회수됨, 히트맵 PNG 저장. 실제 NicheNet v2 행렬(R에서 CSV 내보내기)로는 미실행 |
| ENG-005 | 03_tme_cell_communication/docs/06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:175-260 | `squidpy.gr.spatial_neighbors`+`li.rs.get_resource`+`li.mt.spatial_bivariate`+`li.ut.spatial_local_bivariate` → `li.ut.spatial_neighbors(bandwidth, cutoff, kernel)` + `li.mt.bivariate(resource_name='consensus', local_name='morans', global_name='morans', n_perms, add_categories=True, use_raw=False)`. 반환이 DataFrame이 아니라 AnnData(`lrdata`)이므로 해석부도 수정: 전역 결과는 `lrdata.var['morans','morans_pvals']`(기존 `moran_i`,`pvalue`), 국소 점수는 `lrdata.X`, 핫스팟은 `lrdata.layers['cats']==1`(High-High)로 추출. 헤더 주석 Moran's I→R | ENG-005 (liana 1.10.0 `li.mt.bivariate` 시그니처·반환 구조 직접 조회) | 합성 Visium 유사 AnnData(400 spots, 심은 핫스팟)로 문서 코드를 끝까지 실행: top 쌍 출력(morans 0.154, p=0.0), 국소 점수·High-High 범주 추출, 히트맵 PNG 저장. 플롯 단계는 합성 데이터에 `uns['spatial']`이 없어 `spot_size`·`img_key=None`을 임시로 넣어 확인(문서 코드는 Visium 데이터 가정이라 그대로 둠). `cats`는 High-High(1)/Low-Low(-1)/기타(0)만 구분하므로 High-Low/Low-High 구분은 코드가 아닌 이론 절(L109-127) 설명으로 남음 |
| ENG-006 | 03_tme_cell_communication/docs/00_TME_COMMUNICATION_ROADMAP.md:141; 03_tme_cell_communication/docs/02_DATABASE_CURATION_AND_CONSENSUS.md:141 | `li.rs.get_consensus_resources()` → `li.rs.select_resource('consensus')` | ENG-006; liana 1.10.0 `hasattr(li.rs,'get_consensus_resources')`=False | `li.rs.select_resource('consensus')` 실행: 4,620쌍(ligand, receptor) 반환 확인 |
| ENG-007 | 01_multiome_snatac/docs/07_TF_FOOTPRINTING_AND_CHROMVAR.md:119-160 | `mu.atac.tl.motif_enrichment`, `mu.atac.pp.add_peak_annotation`, `mu.atac.tl.chromvar`(비실재/위치 오류) → pychromvar 흐름(`add_peak_seq`→`add_gc_bias`→`get_bg_peaks`→JASPAR 2024 모티프 `match_motif`→`compute_deviations`)으로 교체. 결과가 `obsm['X_chromvar']`가 아니라 AnnData(cells x motifs, var_names `MA0148.5.FOXA1` 형식)이므로 시각화 부분을 `motif_col()` 헬퍼로 수정(JASPAR 이름이 `Ar`처럼 대소문자 혼재라 대문자 비교). 주석에 muon에 chromVAR 없음, pychromvar 0.0.4는 anndata≥0.13에서 `AnnData(dtype=)` 제거로 TypeError(anndata<0.13 필요) 명시 | ENG-007; muon 0.2.0 `dir(mu.atac.tl)`; pychromvar 0.0.4 소스(`compute_deviations.py`의 `AnnData(dev, dtype=...)`) | 합성 hg38 FASTA(400kb)·300개 피크·200 핵 MuData로 문서 코드 전체(실제 JASPAR2024 CORE 879개 모티프 다운로드 포함) 실행: dev (200x879) 생성, `mdata.obs`에 chromvar_FOXA1/AR/ASCL1 기록, `mu.pl.umap` 오류 없음. anndata 0.13.4에서는 TypeError(주석대로)이고 anndata 0.12.19(`PYTHONPATH`)에서 통과. 합성 데이터라 점수 값 자체(AR std=0)는 의미 없음 |
| ENG-008 | 06_manifold_trajectory/docs/07_FATE_ABSORPTION_PROBABILITIES.md:212 | `cr.pl.fate_probabilities(adata, lineages=[...], basis="umap")` → `g.plot_fate_probabilities(states=[...], basis="umap", same_plot=False)` | ENG-008 (`dir(cr.pl)`에 없음, `GPCCA.plot_fate_probabilities` 존재) | cellrank 2.3.3, 시뮬레이션 데이터로 `set_terminal_states`→`compute_fate_probabilities`→`plot_fate_probabilities(states=['late'], basis='umap', same_plot=False)` 실행 성공. 같은 코드 블록의 `circular_projection(keys='fate_probabilities', color=...)`는 SCI-06-023(S3)로 남김(미수정) |
| SCI-05-001 | 05_generative_foundation_models/scripts/generate_generative_diagrams.py:172,462,980,1039-1067 (및 docs/images/01_scvi_vae_probabilistic_architecture.svg, 03_scanvi_and_scarches_atlas_mapping.svg, 06_in_silico_perturbation_and_lineage_shift.svg) | 방안: "개념도(수치 예시)" 명시 + 근거 없는 수치 삭제. 01: `iLISI 0.94 / cLISI 0.98` → 지표 이름만("값은 데이터 의존"); 03: `98.4% 정확도` → "정확도는 데이터 의존"; 06: `초당 5,000개 유전자`, `cos=+0.89`, `92.4%`, `82% 감소`, `IC50 8.4배`, `감수성 8.4배 회복` 삭제·정성 문구로 교체, 막대(100/38/18%)에 `[예시]` 표기와 하단 주석 "개념도: 막대·수치는 설명용 예시이며 실측 결과가 아님(…results_pipeline_02 참조)" 추가 | pipeline 산출: scANVI 미라벨 정확도 0.5574, 잠재 iLISI 2.224, 이중 KO 억제 2.48%, cos 0.704(`results_pipeline_01/pipeline_summary.json`, `results_pipeline_02/lineage_plasticity_metrics.json`) | (1) 스크립트는 `docs/images`에 6개 SVG만 씀; 편집 전 원 스크립트로 scratch에 6개를 생성해 커밋본과 byte 동일 확인. (2) 편집 후 scratch 생성본과 커밋본 diff: 02/04/05는 0줄, 01·03은 각 2줄, 06은 19줄(의도한 문구뿐). (3) macOS QuickLook 렌더로 텍스트가 박스를 넘지 않는지 확인(두 문구를 줄여 재생성). (4) 01·03·06만 `docs/images`로 복사 후 6개 모두 scratch 생성본과 cmp 동일 |
| SCI-08-001 | 08_large_scale_sparse/docs/04_TILEDB_SOMA_SINGLE_CELL_DATABASE.md:1,8,53-84,159,307,239-242; 08_large_scale_sparse/docs/00_LARGE_SCALE_COMPUTATION_ROADMAP.md:190-191 | 챕터 제목·학습목표 2·§2 제목·§2.1을 "SOMA 기본은 row-major 셀/타일 순서(타일 2048)이고 힐베르트/모턴은 TileDB 코어의 선택 옵션"으로 재작성, "Seek 90% 이상 감소" 삭제(곡선 일반 성질 서술은 '경향'으로 완화), L159·L307의 "힐베르트 곡선 기반 질의"를 "타일 단위 질의"로, 문제 1에 "SOMA 기본값은 행 우선" 단서 추가(힐베르트 vs 모턴 이론 비교 자체는 유지), 로드맵 04장 설명 수정 | SCI-08-001: tiledbsoma 2.3.0 `SparseNDArray._handle.schema_config_options()` → cell_order/tile_order row-major, 차원 타일 2048 | 재실행 확인(tiledbsoma 2.3.0 + anndata 0.12.19): `from_anndata` 후 `{'cell_order':'row-major','tile_order':'row-major','dims tile':'2048'}` 출력. (anndata 0.13.4에서는 `from_anndata`가 TypeError 'quote_from_bytes()'(layers 키 None)로 실패하는 호환성 문제를 별도로 확인; ch04 본문 코드는 미수정) |
| SCI-08-002 | 08_large_scale_sparse/docs/08_GPU_GRAPH_COMMUNITY_DETECTION.md:8,45-84,201-215,345-370; 08_large_scale_sparse/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:132,146,172 | "CAGRA = CUDA Anisotropic Graph", 각도 페널티 가지치기, Tensor Core, "20x~50x" → CAGRA = Cuda Anns GRAph-based(Ootomo et al., ICDE 2024), 초기 kNN 그래프 → **순위 기반 재정렬(detourable route, 순위를 가중치로 사용) + 역방향 엣지 추가**로 고정 차수 그래프 구성, 속도는 논문값(구축 2.2~27x, 대형 배치 처리량 HNSW 대비 33~77x, GPU SOTA 대비 3.8~8.8x)으로 교체. 학습목표·다이어그램·코드 주석·Q2 문제/해설의 "비등방성"과 "80%/85% 이상" 효율 수치 제거. ch09의 허브니스 해법 표현 동일 수정 | arXiv:2308.15136 (초록, 본문 §III 그래프 최적화: rank-based reordering, reverse edge addition, 약어 "Cuda Anns GRAph-based"), `pdftotext`로 본문 확인 | 논문 본문 대조(식 (3) max(w_XZ, w_ZY) < w_XY, 순위 사용). GPU(CuPy/RAFT) 코드는 실행하지 못함(주석·문구만 수정) |
| SCI-10-002 | 10_production_pipelines_mlops/docs/06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md:29,428 | Enzalutamide ID: `CHEMBL1201585`(TRASTUZUMAB) → `CHEMBL1082407`; RxNorm `1306289`(HIV 복합제) → `1307298`; DrugBank `DB08871` → `DB08899` | ChEMBL API `molecule/search?q=enzalutamide` → CHEMBL1082407 ENZALUTAMIDE; RxNav `rxcui.json?name=enzalutamide` → 1307298, `rxcui/1307298/properties` → "enzalutamide" (IN); UniChem API v1(ChEMBL→DrugBank 교차참조) → DB08899(DrugBank 사이트는 직접 확인 못 함) | 모든 새 ID를 API로 직접 조회. DrugBank는 EBI UniChem 교차참조로만 확인(확신도: 중간) |
| SCI-10-003 | 10_production_pipelines_mlops/docs/06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md:29,424,431 | `NCIT:C3343`(Prostate Neoplasm) → `NCIT:C4863`(Prostate Carcinoma); `CL:0002327`(mammary gland epithelial cell) → `CL:0002340`(luminal cell of prostate epithelium); L29의 `CL:0000084`(T cell) → `CL:0002340`; `CL:0000165`, `UBERON:0002367`은 정확하여 유지 | NCI EVS API `concept/ncit/C4863` → Prostate Carcinoma; EBI OLS4 `cl/terms?iri=…CL_0002340` → luminal cell of prostate epithelium (obsolete=False), `CL_0000165`/`UBERON_0002367` 재확인 | API 직접 조회(위 URL 계열) |

## 동일 오류 잔존 위치 (수정하지 않음; 재빌드 또는 후속 작업 필요)

| ID | 잔존 위치 |
|---|---|
| ENG-004 | tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Book.md; TME_Cell_Communication_Master_Course.html; TME_Cell_Communication_Master_Course.docx; backups/master_courses_pre_switcher/Track_03_TME_Cell_Communication_Master_Course.html |
| ENG-005 | 위와 동일한 Track 03 Book/HTML/DOCX/backups HTML + backups/legacy_decks/02_Intracellular_Causal_and_Spatial_CCC_Deck.html |
| ENG-006 | Track 03 Book/HTML/DOCX + backups/master_courses_pre_switcher/Track_03_… HTML |
| ENG-007 | tracks/01_multiome_snatac/Multiome_snATAC_Master_Book.md; Multiome_snATAC_Master_Course.html; Multiome_snATAC_Master_Course.docx; backups/…/Track_01_Multiome_snATAC_Master_Course.html |
| ENG-008 | tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Book.md; Manifold_Trajectory_Master_Course.html; .docx; backups/…/Track_06_… HTML. 같은 블록의 SCI-06-023(`circular_projection color`)도 미수정 |
| SCI-05-001 | tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html(+.docx); slides/01_Generative_VAEs_and_scVI_Deck.html; slides/02_Foundation_Transformers_and_Perturbation_Deck.html; backups/…/Track_05_… HTML (SVG 인라인 사본) |
| SCI-08-001 | tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md; Master_Course.html/.docx; slides/01_Sparse_Matrices_and_Out_of_Core_Deck.html; docs/images/03_tiledb_soma_sparse_tensor_architecture.svg 및 생성 스크립트 scripts/generate_sparse_diagrams.py; scripts/build_sparse_master.py(문구); README.md:37; MASTER_CURRICULUM_ROADMAP.md:201(Hilbert/힐베르트 설명 포함); backups/…/Track_08_… HTML |
| SCI-08-002 | tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Book.md; Master_Course.html/.docx; slides/02_Distributed_Dask_and_GPU_Acceleration_Deck.html; docs/images/06_cagra_gpu_knn_and_cugraph_leiden.svg 및 scripts/generate_sparse_diagrams.py; MASTER_CURRICULUM_ROADMAP.md(Anisotropic 표기); backups/…/Track_08_… HTML |
| SCI-10-002/003 | tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Book.md; Production_Pipelines_Master_Course.html/.docx; backups/…/Track_10_… HTML |

## 미확정/주의 항목

- ENG-007: pychromvar 0.0.4는 anndata 0.13에서 동작하지 않음(주석으로 안내). 대안(ArchR/R chromVAR)은 본문에 미기재.
- ENG-004: NicheNet v2 사전 행렬 CSV 내보내기 절차(R `nichenetr`/Zenodo)는 코드 주석에 간단히만 언급했고 실제 파일로는 실행 못 함.
- SCI-08-001: 본문 ch04의 `soma.io.from_anndata` 사용 코드는 anndata 0.13.4와 비호환(위 검증 메모). 문서 코드는 미수정.
- SCI-08-002: ch08 벤치마크 표(L298-310 부근 CPU 22분 vs GPU 3.82초 등)와 §2.2의 "3.8초"는 가상 수치로 남아 있음(이번 범위 밖, S2 후속 후보).
- SCI-05-001: 챕터 본문(.md)에는 해당 수치가 없었고 그림·생성 스크립트에만 있었다.

## 2단계(전파 오류)

- 규칙: 원본 소스만 최소 수정. 생성물(Book/HTML/DOCX)·빌드 스크립트는 수정하지 않음(Track 05 `slides/build_decks.py`는 지시에 따라 수정). 삭제 명령은 trash만 사용(이번 단계에서 `rm` 미사용). 줄 번호는 수정 후 기준.
- SVG 생성 스크립트는 수정 전에 scratch로 출력해 커밋 SVG와 byte 동일함을 확인(Track 02·03·06·08 모두 diff 0줄), 수정 후 의도한 SVG만 바뀐 것을 확인하고 해당 SVG만 복사. 변경 SVG는 macOS QuickLook 렌더로 텍스트 넘침 확인.

| 영역 | 파일:줄(수정 후) | 변경 요약(전→후) | 근거 | 검증 방법 |
|---|---|---|---|---|
| P-1 Track 08 ch04 코드 | tracks/08_large_scale_sparse/docs/04_TILEDB_SOMA_SINGLE_CELL_DATABASE.md:166-180 | `import tiledbsoma as soma`만으로는 `soma.io`가 없어 AttributeError → `import tiledbsoma.io` 추가, 검증 환경(tiledbsoma 2.3.0 + anndata 0.12.x; anndata 0.13 개발판은 `from_anndata` TypeError) 주석 | 실행 오류 재현: `module 'tiledbsoma' has no attribute 'io'`; anndata 0.13.4에서 TypeError(layers 키 None) | 문서 코드(5천 세포로 축소)를 anndata 0.12.19(PYTHONPATH)로 끝까지 실행: 인제스트·질의(548세포) 성공 |
| P-2 Track 08 ch08 벤치마크 | tracks/08_large_scale_sparse/docs/08_GPU_GRAPH_COMMUNITY_DETECTION.md:85, 298 | §2.2 "3.8초 만에 완료" 삭제(실제 속도는 환경 의존으로 서술), §5 표 도입문을 "가상 시나리오(실측 아님, 설명용 예시 수치)"로 명시 | 표 수치에 측정 근거 없음 | 문구 수정(실행 불필요) |
| P-3 Track 01 ch07 | tracks/01_multiome_snatac/docs/07_TF_FOOTPRINTING_AND_CHROMVAR.md:121-124 | anndata<0.13 조건 주석 유지, R 대안(Bioconductor chromVAR `computeDeviations`, ArchR `addDeviationsMatrix`) 한 줄 추가 | 1단계 검증 | 문구 추가 |
| O Track 06 SCI-06-004/023 | tracks/06_manifold_trajectory/docs/07_FATE_ABSORPTION_PROBABILITIES.md:14, 163, 225-232, 317; 08_LINEAGE_DRIVER_DISCOVERY.md:109-110, 251 | "영원히 머무릅니다/비가역적 고착/되돌릴 수 없으므로" → "분석 모델에서 가정하는 근사(실제 가소성이 영구적이라는 뜻 아님, Ku et al. 2017: EZH2 억제제가 AR 발현·항안드로겐 감수성 회복)"; `cr.pl.circular_projection(keys="fate_probabilities", color=...)` → `keys="cell_type"` | PMID 28059767 초록; `circular_projection` 시그니처, 실행 시 "No valid keys have been selected" 재현 | cellrank 2.3.3 시뮬레이션: `keys="cell_type"` 성공(계통 3개 필요), 기존 호출은 ValueError |
| M Track 08 Hilbert/CAGRA 잔존 | README.md:37; MASTER_CURRICULUM_ROADMAP.md:201; tracks/08_large_scale_sparse/scripts/generate_sparse_diagrams.py(+docs/images/03·06 SVG); slides/01_Sparse_Matrices_and_Out_of_Core_Deck.html, slides/02_Distributed_Dask_and_GPU_Acceleration_Deck.html | 힐베르트/공간 충전 곡선·비등방성 CAGRA·"90% 절감"·"<50ms"·"10x QPS"·"3.8초(500배)" 등을 타일 기반 희소 배열, 고정 차수+순위 기반 재정렬, 논문 보고치(33~77x) 또는 "수치는 환경 의존/가상 예시"로 교체. 슬라이드 HTML은 직접 Edit(인라인 SVG는 생성기와 동일한 치환). ch04 문제 1은 1단계에서 "SOMA 기본은 행 우선" 단서 추가됨(전제 잔존 없음 확인) | SCI-08-001/002 근거(tiledbsoma 2.3.0 스키마, arXiv:2308.15136) | SVG: 03은 32줄, 06은 8줄 변경·나머지 0줄 확인 후 03·06만 복사, 렌더 확인; 슬라이드는 grep으로 Hilbert/Anisotropic/비등방 0건 확인(수정 문구 내 "힐베르트/Z-order는 선택 옵션" 1건 제외) |
| K Track 10 21 CFR Part 11 | tracks/10_production_pipelines_mlops/docs/05_NF_CORE_SCRNASEQ_END_TO_END.md:307-309; docs/06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md:10, 82, 113, 121, 186; docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:203, 224, 230; pipelines/02_mlops_model_serving_and_audit_pipeline.py:20, 338-343, 527, 534, 659 | "준수/무결점 보증/수학적 증명/VERIFIED TAMPER-EVIDENT (Part 11 Compliant)/반려" → "요건 지원 예시·체크섬 기반 변경 탐지·규정 준수 인증 아님"; 매니페스트 필드 `regulatory_compliance` → `regulatory_context`(서명 역할은 "Placeholder (no electronic signature implemented)", `tamper_evident` → `checksum_change_detection`), `compliance_standard` 필드는 `regulatory_context`로 변경; 보고서 문구 "Integrity Check: … not 21 CFR Part 11 compliant" | ENG-036, SCI-10-006, 21 CFR 11.10(밸리데이션·감사추적·접근통제 요건) | 파이프라인 02를 `.venv_api`(python 3.12.13)로 재실행 성공. 갱신된 결과 파일 4개: results_pipeline_02/model_audit_manifest.json, model_weights.json, serving_benchmark.csv, summary_report.md (성능지표(accuracy 99.50%, macro-F1 0.9946)·하이퍼파라미터는 동일, 해시·타임스탬프·벤치마크 지연·python_version(3.8.18→3.12.13)은 변경) |
| L Track 10 화이트리스트 | tracks/10_production_pipelines_mlops/pipelines/01_dsl2_dataflow_and_dag_engine.py:77, 323, 334, 435 | 챕터 05가 `--protocol 10XV3`인데 파이프라인이 v2 화이트리스트(`737K-august-2016.txt`)를 사용 → `3M-february-2018.txt`(10x v3)로 통일(설명·명령·보고서 문구) | STARsolo 문서(V2=737K-august-2016, V3=3M-february-2018), Track 00 SCI-00-023 | 파이프라인 01 재실행 성공. 갱신된 결과 파일: results_pipeline_01/summary_report.md, timeline_metrics.json, trace_table.csv, .nextflow_cache.json(작업 해시·타임스탬프 변경) |
| J-1 Track 06 | scripts/generate_trajectory_diagrams.py; docs/images/01·03·05·06 SVG; docs/08_LINEAGE_DRIVER_DISCOVERY.md:118, 238 | SVG06의 확률·상관 수치에 "(예시)"와 하단 "개념도: 수치는 설명용 예시이며 실측 결과가 아님" 주석, "최상위 구동 인자로 등극/주도 확인" → "구동 후보로 예시"; SVG01 `M=(I−T)^-1`(특이행렬) → `M=Σ_t(T^t−φ0φ0ᵀ)`; SVG03 "20-30회 반복" → "EM 반복(기본 max_iter=10)", τ=0.5 → "τ 중간(예시)"; SVG05 커널 가중치 0.4/0.3/0.3 → w1/w2/w3, "k=4 자동 결정" → "n_states 지정"; ch08 §4 도입과 Q2를 "가상 시나리오 예시(실제 코호트·파이프라인 r≈0.58과 다름)"로 명시 | SCI-06-050/052/053/054; scvelo `recover_dynamics` max_iter=10; 파이프라인 02 결과 | SVG 4개 변경(01:2줄, 03:4줄, 05:8줄, 06:9줄), 02·04는 0줄, 렌더 확인 |
| J-2 Track 03 | scripts/generate_tme_diagrams.py; docs/images/01·02·04·05 SVG | SVG04 AUROC 0.88/0.79/0.74 → 파이프라인 02 산출 0.96/0.92/0.99(막대 폭도 조정, "toy 합성 데이터 산출값" 표기); SVG02 "파이프라인 간 일치도 75% 이상 향상" → "다중 방법 합의로 방법 간 편차 완화(개념)"; SVG05 하단 "개념도: 비율은 예시" 주석; SVG01 내분비 ">1mm" → ">1cm"(본문과 일치) | results_pipeline_02/pipeline_summary_metrics.json (TGFB1 0.96181, CXCL12 0.91667, TNFSF11 0.98919), SCI-03-043 | SVG 4개 변경(01:2줄, 02:2줄, 04:14줄, 05:1줄), 03·06은 0줄, 렌더 확인 |
| J-3 Track 02 | scripts/generate_spatial_diagrams.py; docs/images/05·06 SVG; docs/06_SPATIALLY_VARIABLE_GENES.md:178, 320; docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:70 부근 | SVG05 `I=0.74/0.68/≈0.02` → "높은 I(예시)"; SVG06 `r=30µm, k=20, e^(-d/30), d_0=30µm` → 일반 표기(r 예: 30~50, k, d_0 신호 유형별); ch06 `Moran's I > 0.65` → "값은 데이터 의존, 합성 파이프라인 KLK3 I=0.565"; Q3 수치에 "(가상 환자 예시 값)"; ch08에 "파이프라인은 d_0=30 µm/r=30 µm toy 설정이라 본문 200~250 µm와 직접 비교 불가" 참고 추가 | pipelines/results_pipeline_01/spatial_variable_genes_moran.csv (KLK3 0.56538), SCI-02-041/049/054 | SVG 2개 변경(05:4줄, 06:8줄), 나머지 4개 0줄, 렌더 확인 |
| N Track 05 | slides/build_decks.py; slides/01_Generative_VAEs_and_scVI_Deck.html; slides/02_Foundation_Transformers_and_Perturbation_Deck.html; docs/07_GENEFORMER_AND_GENOME_ATTENTION.md:351 부근; docs/08_IN_SILICO_PERTURBATION_AND_DRUG_DISCOVERY.md:96-107 | 덱 본문의 "0.5% 미만 클론 98.7% 식별", "10만 세포 5분", "위양성 없는", "iLISI 0.94+/cLISI 0.98+ 등 임계", "증명/규명/원천 차단" → 데이터 의존·탐색·가상 예시 표현; 덱 SVG(01·03·06)는 1단계 새 SVG로 교체; 챕터 08의 가상 시나리오 "입증/규명" → "개념적으로 일치하는 가상 시나리오 예시/가설화", cos 값이 toy 파이프라인(이중 KO cos≈0.70, 억제 2.48%)과 다름을 명시; ch07의 "증명합니다" → "시사합니다(가상 사례 해석)" | SCI-05-001, results_pipeline_02 | **덱은 재생성하지 않고 동일 치환을 HTML에 직접 적용**: 커밋 덱이 build_decks.py 출력과 이미 다름(수작업 수정 반영: 덱1 L674·678 scANVI 손실식, 덱2 L197·204·551·552 등 총 14줄). 재생성하면 이 수정이 되돌려지므로 build_decks.py(소스)와 HTML을 각각 같은 치환으로 수정. 변경 전후 diff: 덱1 16줄, 덱2 29줄(의도한 문구+SVG만) |

## 2단계 동일 오류 잔존·범위 밖 발견

- **build_decks.py ↔ 커밋 덱 불일치**(Track 05): 덱에 수작업 수정이 있어 재빌드 시 되돌려짐. 후속에서 build_decks.py에 그 수정(덱1 scANVI 손실식, 덱2 토큰화 문구 등)을 반영해야 함.
- Track 08: docs/images/05(RAPIDS "Tensor Core cuML 100x 가속" 등)와 06 패널 3("100만 세포 군집화 5초 돌파, 5.2s")의 수치, 슬라이드 02의 같은 SVG 사본, ch07 본문 Tensor Core 서술은 이번 범위 밖(근거 없는 수치 가능성).
- Track 06: SCI-06-022(엔트로피 표·Q3 수치 산술 불가), SCI-06-007(DPT 가중치 본문), SCI-06-048~051(파이프라인 02 "G-PCCA"=라벨 복사, 0.5/0.5 폴백, latent time=0.85·DPT+0.15·coherence)은 이번 지시 범위 밖이라 미수정.
- Track 03 SVG05의 비율(45/20/35)·SVG06 수치는 개념도 주석만 추가(수치 자체는 예시로 남음). Track 02 SVG01~04의 다른 수치는 미점검.
- 생성물(Book/HTML/DOCX/backups)과 PPTX/DOCX 슬라이드에 1·2단계 수정 전 문구가 그대로 남아 있음(재빌드 필요).

## 2단계 보정: Track 10 결과 파일 환경 차이 되돌림

- 원인: 앞서 `.venv_api`(py3.12/numpy 2)로 재실행한 결과 `results_pipeline_02/model_weights.json`(수천 줄)과 `results_pipeline_01/trace_table.csv`가 텍스트 수정과 무관하게 바뀜(py3.8 대비 float 직렬화/시각 필드 차이).
- 재현 확인: scratch에 pipelines 복사본을 만들고 원작성 환경 `<conda-base>/envs/scRNAseq_anal/bin/python`(py3.8.16)로 두 파이프라인을 재실행해 `backups/post_s1_fixes_2026-10-18.tar.gz`의 결과와 비교. `model_weights.json`은 **byte 동일**(py3.12 결과가 원인이었음). 차이는 시각·소요시간·벤치마크 지연값(실행마다 변동), `python_version`(3.8.18→3.8.16), 의도한 텍스트(Part 11 문구), 화이트리스트 변경이 해시에 미친 Demux 행뿐.
- 처리: 결과 8개 파일(+hidden cache)을 tar에서 추출해 스냅샷 버전으로 덮어써 복원한 뒤, 의도한 결정적 변경만 텍스트로 적용.
  - results_pipeline_02/model_weights.json, serving_benchmark.csv: 스냅샷 그대로(복원).
  - results_pipeline_02/model_audit_manifest.json: 스냅샷 + `regulatory_compliance` 블록 → `regulatory_context` 블록(10줄 diff)만. timestamp·python_version은 스냅샷 값 유지.
  - results_pipeline_02/summary_report.md: 스냅샷 + 2줄(Regulatory Context, Integrity Check)만(4줄 diff).
  - results_pipeline_01/timeline_metrics.json: 스냅샷 그대로(복원).
  - results_pipeline_01/trace_table.csv: 스냅샷 + Demux 8행의 hash_key·work_dir만 py3.8.16 재실행값으로 교체(명령어 문자열이 바뀌어 해시가 바뀌는 결정적 변경, 16줄 diff). 타임스탬프·소요시간은 스냅샷 유지.
  - results_pipeline_01/summary_report.md, .nextflow_cache.json: py3.8.16 재실행 결과(차이는 화이트리스트 문구 1줄 + Demux 작업 디렉터리/해시뿐: 10줄/8줄 diff).
- 최종 스냅샷 대비 변경: 위 5개 파일만, 모두 의도한 텍스트 또는 화이트리스트 유래 해시. 위 "2단계" 표의 K·L 행에서 "python 3.12로 재실행, 4개/4개 파일 갱신, 성능지표 동일" 서술은 이 보정으로 대체됨(실제 최종 변경 파일: pipeline_02 manifest·summary_report, pipeline_01 trace_table·summary_report·.nextflow_cache.json).

## S2 단계 / ENG

대상: findings_engineering.md S2 19건. 상태는 `review/s2_status_eng.csv`에 기록. 결과: fixed 12 / already-fixed 2 / softened 1 / deferred 3 / rejected 1.
검증 환경: 빌드는 scratch 복사본(`scratchpad/eng/s2/bproj`, WORKSPACE 경로는 복사본에서만 치환)을 `basic_research`(py3.8.16, mistune 2.0.5)로, 파이프라인은 `scRNAseq_anal`(py3.8.16)로 실행. 원본 프로젝트에서는 빌드·파이프라인을 실행하지 않음. 삭제는 trash만 사용.

| ID | 상태 | 파일:줄(수정 후) | 변경 요약(전→후) | 검증 |
|---|---|---|---|---|
| ENG-009 | fixed | tracks/05_generative_foundation_models/docs/04_SCARCHES_ARCHITECTURAL_SURGERY.md:~195 | `scvi.data.prepare_query_anndata(q, reference_model=)` → `scvi.model.SCVI.prepare_query_anndata(q, ref_model, inplace=False)` | scvi-tools 1.5.1: 합성 ref/query로 prepare → `load_query_data(freeze_* , freeze_classifier)` 성공 |
| ENG-010 | deferred | (미수정) tracks/01_multiome_snatac/docs/03_SNATAC_QC_METRICS.md:115 | Track 01은 07장만 eng 범위라 수정 후 되돌림. 검증된 패치안: `snap.pp.import_fragments(fragment_file=, chrom_sizes=snap.genome.hg38, min_num_fragments=500, sorted_by_barcode=False)` | snapatac2 2.10.0: 합성 fragment 40 barcode → import_fragments → tsse → filter_cells 실행 성공 |
| ENG-011 | deferred | (미수정) docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:30; tracks/07_count_statistics_abundance/docs/03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:220 | root docs/·Track 07은 sci-checker 범위라 수정 후 되돌림. 검증된 패치안: `dc.pp.pseudobulk(adata, sample_col, groups_col, mode="sum")` + `dc.pp.filter_samples(pbulk, min_cells=, min_counts=)`(in-place) | decoupler 2.2.0 합성 AnnData 실행(필터가 in-place로 n_obs 감소 확인) |
| ENG-012 | fixed | tracks/08_large_scale_sparse/docs/05_DASK_LAZY_EVALUATION_DISTRIBUTED_QC.md:~406-410 | `Group.create_dataset` → `create_array`(+zarr 2.x 안내 주석) | zarr 3.4.1: create_array, r+ 슬라이스 기록 성공 |
| ENG-013 | softened | tracks/02_spatial_transcriptomics/docs/04_CELL_SEGMENTATION_ALGORITHMS.md:~161 | Cellpose 3.x 전용(`cellpose<4`) 주석, 4.x는 `CellposeModel`·eval 3값 | 4.2.1.1 시그니처 확인만(가중치 다운로드 필요로 미실행) |
| ENG-014 | fixed | tracks/06_manifold_trajectory/docs/04_OPTIMAL_TRANSPORT_WADDINGTON_OT.md:~224-262 | `wot.Dataset`+`OTModel(..., growth_rate=)`+`compute_transport_matrix` → AnnData 입력 `OTModel(adata, day_field, growth_rate_field, epsilon, lambda1, lambda2)`+`compute_transport_map` | wot 1.0.8.post2: 문서 함수를 md에서 추출해 합성 데이터로 전 시점 쌍 실행 성공 |
| ENG-020 | fixed | scripts/build_master_documents.py 및 tracks/01~10 build_*_master.py (build_master_html 시작부, 장 렌더 직후) | mistune>=3이면 RuntimeError(2.0.5 안내) 추가; 장별 렌더 후 코드/pre 밖 리터럴 `**` 잔존 assert 추가 | basic_research(mistune 2.0.5) 11/11 빌드 8 assertion 통과(추가 assert 포함). (Multiome HTML의 `***`는 SVG의 유의성 별표로 정상) |
| ENG-021 | already-fixed | scripts/inject_global_track_switcher.py | S1 이전 오케스트레이터가 BEGIN/END 마커 기반 idempotent 주입으로 수정 | scratch에서 빌드→inject 11/11 성공 |
| ENG-022 | fixed | tracks/06_manifold_trajectory/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:126 | 표 셀의 `$|r| > 0.5$` → "절댓값 $> 0.5$"(Track 05 로드맵 표는 현 md에서 이미 해소) | scratch 빌드 11개 HTML `<p>|` 0건 |
| ENG-024 | fixed | pipelines/04_*.py, tracks/01/pipelines/01·02, tracks/02/pipelines/01·02 (sys.path 블록) | 개인 conda 경로 하드코딩 삭제, `SCRNA_ENV_SITE_PACKAGES` 지정 시에만 `sys.path.append` | `scRNAseq_anal`로 scratch 재실행(5개 모두 exit 0): Track 01(2)·Track 02(2)는 pre-S2 스냅샷과 전 파일 byte 동일, root 04 byte 동일 |
| ENG-030 | rejected | (미수정) pipelines/02_prostate_qc_and_cellbender.py | 원작성 환경(scRNAseq_anal, scrublet 0.2.3)에서 exit 0으로 정상 완료(83 doublets, h5ad 저장). 실패는 py3.12+annoy 1.17.3 환경 한정. 환경 정보는 requirements-pipelines.txt에, 로드맵 문구는 환경 안내로 완화(ENG-033) | scRNAseq_anal scratch 실행 |
| ENG-031 | fixed | pipelines/03_prostate_infercnv_and_pseudobulk.py:~123-132 | 검정을 수행하지 않는 스텁이 "p<1e-100…재현 가능한 임상 소견" 문구를 출력 → "스텁, 검정 미수행, 실제 분석은 PyDESeq2" 안내로 교체, docstring에 스텁 명시 | scRNAseq_anal scratch 실행 확인(결과 파일 생성 없음) |
| ENG-032 | fixed | requirements-pipelines.txt (신규) | scRNAseq_anal(py3.8.16)에서 `importlib.metadata`로 읽은 실제 버전만 기록, 검증 범위 명시 | 환경에서 직접 조회 |
| ENG-033 | fixed | MASTER_CURRICULUM_ROADMAP.md(10곳, 103행), README.md:93, results 폴더 21곳 SIMULATION_NOTICE.md | "프로덕션 …" → "교육용(합성 데이터 기반)", 환경 파일 안내; 결과 폴더에 합성 데이터 고지 파일 추가(결과 파일 불변) | 빌드 8 assertion(README/ROADMAP 필수 용어) 통과 |
| ENG-036 | already-fixed | tracks/10/.../summary_report.md 등 | 2단계 K에서 해소 | — |
| ENG-038 | fixed | 11개 빌드 스크립트 CSS | 02~10에 `@media (max-width:1024px)` 추가(사이드바 세로 쌓기, main 100%, `.source-ref`·인라인 code 줄바꿈), 00·01 기존 블록 보강 | 375×812 실측(http 서버+Playwright, inject 포함): 00·05·06 scrollWidth=375, 넘침 요소 0 |
| ENG-039 | fixed | 동일 CSS | `table:not(.math-matrix){display:block;overflow-x:auto}` 등 | 375px: 00의 표 7개·05의 3개·06의 2개가 내부 스크롤(이전엔 클리핑) |
| ENG-040 | deferred | — | 23개 덱이 서로 다른 생성기(build_decks.py, make_decks.py, generate_*_slides.py)·수작업 HTML에서 나오고 해상도별 레이아웃 재설계와 덱별 시각 검증이 필요 → 규칙 4에 따라 보류 | — |
| ENG-041 | fixed | 11개 빌드 스크립트 `code_runs = "<w:r><w:br/></w:r>".join(...)` | 코드 줄바꿈 `<w:br/>`를 `<w:r>` 안으로 | 11개 DOCX: `w:br`가 `w:r` 밖에 있는 경우 0(전: 332~1,757). 기준 빌드(pre-S2 스크립트)와 정규화 비교 시 차이는 br 위치뿐. Word 표시는 미검증(macOS textutil은 원래도 줄바꿈을 공백 처리) |

### 빌드 검증 요약 (scratch)
- 11개 트랙(00~10)을 최종 스크립트로 `basic_research`에서 빌드: 모든 assertion 통과(Track 05는 파이프라인 재실행 포함 ~3분). 빌드 전 `SCRNA_ENV_SITE_PACKAGES=…/scRNAseq_anal/lib/python3.8/site-packages`를 지정(ENG-024로 하드코딩이 사라졌으므로 basic_research로 Track 01/02/root 파이프라인을 돌릴 때 scanpy를 쓰려면 필요; 미지정 시 scanpy 폴백 경로).
- 기준 대비 출력 변화: 같은 내용으로 pre-S2 스크립트로 빌드한 기준(00·01·02·05·06·07·08)과 비교하면 Book.md 동일, HTML은 CSS 추가만(+14줄: 00·01, +36줄: 02~10), DOCX는 `w:br` 위치만 변경. 내용이 바뀌지 않은 트랙(03·04·09·10)은 live 생성물 대비 inject 후 HTML이 +36줄(CSS)뿐.
- 빌드가 재실행한 파이프라인 결과는 scratch 안에서만 갱신되며 live 결과와 비교해 시간 필드 외 차이 없음(Track 07 metrics의 execution_time 등).

### 범위 관련 주의
- ENG-010/011은 위치가 sci-checker 담당 범위여서, 처음에 수정했다가 지시 범위(금지 목록)를 확인한 뒤 해당 줄만 되돌렸다(다른 사람이 수정한 부분은 건드리지 않음). 검증된 패치안은 위 표에 있음.
- 새 결과 폴더 SIMULATION_NOTICE.md 21개 중 Track 04·07·09 폴더 3×2개 포함(결과 파일은 불변).
- 작업 중 Playwright가 프로젝트 루트에 만든 `.playwright-mcp/`는 trash로 치웠다.

### Track 10 (S2 22건: SCI-10-001 외)
상태: fixed 18 / softened 3 / already-fixed 1 / deferred 0. 대상은 `tracks/10_production_pipelines_mlops/docs/00~09` 문서뿐이며 pipelines·results·Book/HTML/DOCX·슬라이드는 수정하지 않았다(파이프라인 결과 불변; 2026-10-19 기준 pipelines 하위 변경 파일은 ENG-033의 SIMULATION_NOTICE.md 2개뿐).

| ID | 상태 | 수정 위치 | 수정 내용 | 근거 |
|---|---|---|---|---|
| SCI-10-001 | fixed | 00 트리·§4.1 | 05–09 장 파일명·설명을 실제 파일로 교체 | 파일 목록 대조 |
| SCI-10-006 | already-fixed | — | 2단계 K에서 해소 | 위 K 행 |
| SCI-10-009 | fixed | 01 Singularity.def | 존재하지 않는 static 자산 → `STAR_2.7.11b.zip` + unzip | GitHub 릴리스 자산 확인 |
| SCI-10-010 | fixed | 01 런타임 단계 | `libcurl4` → `libcurl3-gnutls` | packages.ubuntu.com(jammy) |
| SCI-10-013 | fixed | 02 star 태그, 05 scanpy | `star:2.7.11b--h43eeafb_0`, `scanpy:1.7.2--pyhdfd78af_0`(biocontainers 최신), `python-annflow` 제거 | quay.io API, anaconda.org 조회 |
| SCI-10-014 | fixed | 03·08·00 | SHA256/MD5 내용 해시 서술 → 작업 해시(경로·크기·mtime, lenient/deep), 수식은 개념식 | Nextflow v24.04.4 docs `cache` |
| SCI-10-015 | fixed | 03 | `--outFileNamePrefix …/{sample}_` (2곳) | STARsolo 출력 규칙 |
| SCI-10-016 | fixed | 03 | `--rerun-triggers code input params software-env` | Snakemake CLI |
| SCI-10-019 | fixed | 04 | `check_max` 제거 → `process.resourceLimits` | Nextflow 24.04 `resourceLimits` (08의 check_max time 분기는 별건 SCI-10-040) |
| SCI-10-020 | fixed | 04 §3.2·Q2, 08 §2.2·§4, 09 실패모드 6·§2.2·루브릭 | "143=스팟/자동 전환" 삭제 → `aws.batch.maxSpotAttempts`·`google.batch.maxSpotAttempts`, queue 전환은 `task.attempt` 사용자 설정 | Nextflow v24.04.4 docs(aws/google batch) |
| SCI-10-021 | fixed | 04 §4·Q3, 00 | "제로 로컬 스토리지/0.2초/0바이트" 삭제 → 컨테이너 /tmp 스크래치, NVMe 권장, 소규모 EBS gp3 대안 | docs.seqera.io Fusion AWS Batch 가이드 |
| SCI-10-022 | softened | 04 §1.1·Q1 | Fair Tree 기본 명시, 대기 작업 미반영은 추정 표기, queueSize는 U를 줄이지 않음 | SchedMD 문서(검증 한계: 추정은 Assumption) |
| SCI-10-023 | fixed | 01 §4.2, 04 §1.2 | `--readFilesCommand zcat` 추가 | STAR 매뉴얼(실행 미수행) |
| SCI-10-027 | fixed | 05 §5, 06 CI 매트릭스 | `--aligner alevin`(2.5.1), 없는 `--save_*` 제거, 샘플시트 `sample,fastq_1,fastq_2` + 별도 공변량 CSV | nf-core/scrnaseq 2.5.1·4.1.0·master `nextflow_schema.json`, `assets/schema_input.json` |
| SCI-10-028 | fixed | 05 모듈 | `Series.mad()` → `median_abs_deviation(scale='normal')`, 5-MAD 필터를 실제 적용 | pandas 2.0 제거; scRNAseq_anal 함수 단위 확인 |
| SCI-10-031 | fixed | 06 §3.2·목표·docstring, 09 | S3 ETag/composite 체크섬 설명으로 교체, 머클 트리는 일반 개념, O(1) 주장 정정 | AWS S3 무결성 문서 |
| SCI-10-037 | softened | 07 serve_model.py | 임계값 예시 표기, 분모=종양 세포, 서브클론 임계 2%(Q3 4%와 일치), 백금 권고 문구 삭제 | 코드 판독; python 블록 compile 확인 |
| SCI-10-038 | fixed | 07 Dockerfile·compose | `/opt/venv`, python 헬스체크, compose `version` 삭제 | 구성 판독(docker 빌드 미수행) |
| SCI-10-039 | softened | 07 Q3 | "잡음 수준" 단정 완화(0.1→6.1 TPM 약 60배), 비가역·골든타임 표현 완화 | 챕터 내 산술 |
| SCI-10-045 | fixed | 09 실패모드 1 | 비존재 옵션 삭제, `Can't open cache DB` 증상·`-w`/`NXF_CACHE_DIR` 해법 | nextflow-io/nextflow 이슈 #1479·#6996, v24.04 config 문서 |
| SCI-10-046 | fixed | 09 §2.2 | `dynamicDelay` 삭제 → 공식 sleep 백오프 패턴, `8.GB * task.attempt`, `ignore` → `finish` | Nextflow v24.04.4 docs "Dynamic Retry with backoff" |
| SCI-10-047 | fixed | 09 §3.2 | `-k`=keep-logs 정정, 드라이런·`-resume` 소실 경고, `--dryrun` | Nextflow v24.04.4 CLI docs `clean` |

검증: Track 10 docs 00~09 코드펜스 짝 맞음, 07의 python 블록 compile 통과(03의 Snakefile 블록은 python 라벨이라 해당 없음), 09 ASCII 박스 폭 정렬 확인. 남은 한계: Dockerfile/Singularity.def 이미지 빌드, STAR/Nextflow 실제 실행, nf-core 2.5.1 샘플시트 검증 실행은 하지 않았다.

### Track 08 (S2 20건: SCI-08-003 ~ 022)
상태: fixed 18 / softened 2 (SCI-08-009, 021) / already-fixed 0 / deferred 0. 수정 범위: `tracks/08_large_scale_sparse/docs/00~09`, 그리고 `scripts/generate_sparse_diagrams.py` + `docs/images/05·06 SVG`(생성 스크립트와 SVG를 함께 수정, 01~04 SVG는 재생성해도 바이트 동일). Book/HTML/DOCX·슬라이드(슬라이드 02에 구 SVG05 사본 포함)는 수정하지 않았으므로 오케스트레이터 재빌드/교체가 필요하다. 1·2단계 처리분(Hilbert/CAGRA 문구, ch08 표의 "가상 시나리오" 표기, ch04 `import tiledbsoma.io`)은 그대로 유지.

| ID | 상태 | 수정 위치 | 수정 내용 | 근거/검증 |
|---|---|---|---|---|
| SCI-08-003 | fixed | 08 목표·§3.3·§5 표·Q3, 09 실패양식 9·루브릭 | Leiden "100% 분리" 삭제, 해상도 한계 vs γ-연결성 구분, 정제 조건식 논문 형태로 정정 | Traag 2019, Fortunato&Barthélemy 2007 |
| SCI-08-004 | fixed | 05 다이어그램·Q2, 09 실패양식 1 | "나눗셈 밀집화" 삭제 → `toarray()`/`X - mean` 사례 | scipy 1.16.3 직접 실행 (`m/cs`→coo, `m-mean`→np.matrix) |
| SCI-08-005 | fixed | 05 의사코드·코드 주석·Q2, 09 | "0바이트/O(1)" → data 크기 1배 임시 버퍼 | tracemalloc 실측; sklearn `inplace_row_scale` 도 np.repeat 사용 확인 |
| SCI-08-006 | fixed | 02 §4(다이어그램·본문·Q1), 00 Q, 03 Q1, 05 Q3 | 직렬화 서사·수치 삭제, 쓰기 잠금 충돌 + fork 상속 핸들 | HDF5 File Locking 문서(RDWR=배타, 그 외=공유) |
| SCI-08-007 | fixed | 00 §2.3·Q1, 01 COO 문장·Q2 | 72.04 GB(단일 index dtype int64), 오버플로 대상 indptr, COO 30/18 GB | scipy `get_index_dtype` 직접 실행 |
| SCI-08-008 | fixed | 01 Q1 전략 B | `var` 부재 → `mean_variance_axis` | `hasattr` 실행, dense 계산과 일치 확인 |
| SCI-08-009 | softened | 05 §6, 06 §6·§1, 07 §6 | 예시(가상) 수치 명시, 내부 모순 5건 정정 | 공개 수치(rapids-singlecell 논문, bioRxiv 2025) 병기 |
| SCI-08-010 | fixed | 06 다이어그램 | FLOP→시간 환산, U 4TB 조건 | 산술 재계산 |
| SCI-08-011 | fixed | 07 목표·§4.1·§4.2·다이어그램, SVG05 | 텐서코어/수십 ms 삭제, rsc 자체 PCA로 정정 | 2NG²=1.8e13 FLOPs; rsc stable 문서 |
| SCI-08-012 | fixed | 07 파이프라인·Q3 | seurat_v3를 정규화 전, HVG 부분집합 후 scale, `rsc.pp.pca` | scanpy/rsc 문서 (GPU 미실행) |
| SCI-08-013 | fixed | 03 §3.1·Q1 | 비용 서술을 자체 계산과 일치 | 챕터 산식 재계산 |
| SCI-08-014 | fixed | 03 Q1 해설 | 2.4GB 청크 오류 정정 | 산술 |
| SCI-08-015 | fixed | 03 §4 | 존재하지 않는 경로 → placeholder + `storage_options` | 404 확인(기존), Census=SOMA |
| SCI-08-016 | fixed | 03 §5 코드 | 변환·벤치마크 재작성 | 합성 h5ad로 실행, zarr 3.4.1 |
| SCI-08-017 | fixed | 05 §5 코드, 03 §3 | 파티션별 쓰기 + 단일 프로세스 병합 | 합성 데이터로 참조와 max diff 0.0 |
| SCI-08-018 | fixed | 04 §1 | 약어 정정, 연대 삭제 | SOMA README |
| SCI-08-019 | fixed | 04 코드·Q3 | `prostate gland`, `raw_sum` | Census 2023-12-15 직접 질의 |
| SCI-08-020 | fixed | 02 §5·Q2, 05 Q3 | 희소 1D 청크 vs 밀집 2D 구분 | anndata 0.13.4 쓰기 후 h5py 청크 확인 |
| SCI-08-021 | softened | 09 §2.1·Q1 | 단편화 서술 완화, RMM 클래스명 정정 | rmm 소스(CudaAsyncMemoryResource) |
| SCI-08-022 | fixed | 09 코드 | `svd_solver="jacobi"` | cuML 소스(auto/full/jacobi) |

SVG 수정(수치 근거 없음 → 완화): SVG05 제목 "Tensor Core cuML 100x" → "cuML GPU 가속 (예시 수치)", "(Tensor Cores)"·"PCA 120배"·"100% 호환"·"(100x/103x/96x)"·"42분→26초"·"비용 92% 절감"·"커피 한 잔" 삭제/완화; SVG06 "5.2s Converge"·"5초 돌파"·"CPU 25분→5.2초"·"치명적 결함"·"완전한 위상학적 연결성" 완화. 변경은 05에서 14행, 06에서 5행(diff로 확인), 나머지 SVG 4개는 재생성 결과 불변.

검증: docs 00~09 코드펜스 짝 맞음, python 블록 전부 compile 통과. 05 파이프라인(파티션 쓰기+병합)과 03 변환 코드는 합성 h5ad로 실행해 결과 확인(scratch, `.venv_api`). 남은 한계: GPU 코드(07·09)·Dask 분산 클러스터 실행 미수행, 05 Q의 15MB/s 등 일부 가정 수치는 예시로 남아 있음.

### Track 05 (S2 23건: SCI-05-002 ~ 048)
상태: fixed 21 / softened 1 (SCI-05-032) / already-fixed 1 (SCI-05-024, ENG-009에서 `inplace=False` 반영) / deferred 0. 수정 범위: `tracks/05_generative_foundation_models/docs/00~09`, `scripts/generate_generative_diagrams.py` + `docs/images/01·03·05 SVG`(생성 스크립트와 함께 수정; 02·04·06 SVG는 재생성해도 바이트 동일, 수정 전 출력은 커밋된 SVG 6개와 byte 동일함을 먼저 확인). Book/HTML/DOCX·슬라이드는 미수정.

| ID | 상태 | 수정 위치 | 수정 내용 | 근거/검증 |
|---|---|---|---|---|
| SCI-05-002 | fixed | 00 환경검증 코드 | `probs = mu/(theta+mu)` | torch NB 평균 직접 확인 |
| SCI-05-003 | fixed | 00 §로그 정규화 | 오귀속 정정, 두 논문 결론을 균형 서술 | finding abstract 근거 |
| SCI-05-004 | fixed | 00 §3.1 표 | Dong 행→Chen 2021, He·Labrecque 서지 정정 | PubMed 직접 조회 (PMID 33420488/33664492/31361600) |
| SCI-05-008 | fixed | 00·02, SVG01 | 배치 불변성 단정 삭제, 기본값 encode_covariates=False 명시 | scvi-tools 1.5.1 시그니처 |
| SCI-05-009 | fixed | 01·02·09 | 기본 우도 zinb 정정 | 소스(0.20~1.5) |
| SCI-05-014 | fixed | 02 코드 | `transform_batch` 지정 | 합성 데이터 02 파이프라인 실행 |
| SCI-05-015 | fixed | 02 코드 | 실재 군집 라벨 선택 로직 | 합성 데이터 실행(DE까지 완주) |
| SCI-05-018 | fixed | 03 코드 | `get_normalized_expression(indices, gene_list)` | 합성 데이터 호출 확인 |
| SCI-05-023 | fixed | 04 코드·§2.3 | query_model 재인코딩 검사, weight_decay=0.0 | 합성 데이터 04 파이프라인 실행 (drift 0) |
| SCI-05-024 | already-fixed | 04 코드 | ENG-009 | 04 파이프라인 실행에서 통과 |
| SCI-05-025 | fixed | 04 §4 | DP 단정 완화 | (ε,δ)-DP 정의 |
| SCI-05-029 | fixed | 00·05·07·09 | corpus-정규화 rank-value, 불변성 한계 | Geneformer 모델카드 |
| SCI-05-032 | softened | 00·05·07 | 어텐션=조절 에지 단정을 가설로 | Geneformer abstract; 직접 문헌 조회 안 함 |
| SCI-05-033 | fixed | 06, SVG05 | GEP/GEPC로 재서술, MGM은 Geneformer로 | scGPT 코드·README (논문 본문 미확인) |
| SCI-05-034 | fixed | 06, SVG05 | 제로샷 소거 단정 삭제, 한계·벤치마크 | Kedzierska 2025 |
| SCI-05-036 | fixed | 06 §5 | 장난감 모델 명명, "제로샷" 삭제 | 코드 판독 |
| SCI-05-040 | fixed | 07 §5·Q3 | 가상 예시 명시, 순환 주석 | 코드 판독/실행 |
| SCI-05-041 | fixed | 07 §4.2·코드 | 통합 엔트로피로 교체 | 07 코드 실행 확인 |
| SCI-05-042 | fixed | 07 Q3 | 프로모터 해석 삭제 | 모델 입력 정의 |
| SCI-05-044 | fixed | 08 Q3·§4 | 수치 일치·가상 값 명시, 기준선 문단 | 문서 내 대조 |
| SCI-05-045 | fixed | 00·08 | 정확도 단정 삭제, 기준선 평가 필수 | Ahlmann-Eltze 2025 |
| SCI-05-046 | fixed | 08 코드 | ToyScGPTEncoder 내장 | 블록 실행 확인 |
| SCI-05-048 | fixed | 07·08 | Tarlatamab/Tazemetostat 정정 | finding 근거(OncLive, FDA) — 직접 재조회 안 함 |

주의(다른 트랙 동일 오류는 note만): Track 05 SVG03의 "수주일에서 수초" 문구는 완화, SVG04(어텐션 토큰화)·SVG06은 변경 없음.

### 부속 작업 (Track 05 배치)
1. **Track 08 슬라이드 02 SVG 교체**: `tracks/08_large_scale_sparse/slides/02_Distributed_Dask_and_GPU_Acceleration_Deck.html` 안의 SVG05와 SVG06 사본(둘 다 수정 전 docs/images와 byte 동일했음)을 새 SVG 내용으로 교체. SVG06 사본에도 "5.2s/5초 돌파" 문구가 있어 함께 교체. 해당 `<svg>…</svg>` 두 블록만 바뀜(diff 38행 = 19행 쌍), "Tensor Core cuML 100x"·"5.2s" 잔존 0건.
2. **SCI-02-066 (sci→eng)**: `tracks/02_spatial_transcriptomics/pipelines/01_spatial_qc_and_deconvolution_pipeline.py`의 docstring(L16~)과 Step 4 구분 주석에서 "Cell2location-style Hierarchical Poisson"을 "Poisson NMF-KL multiplicative-update solver with a known signature matrix (NOT cell2location)"로 정정, 실제로 구현돼 있지 않은 "q05 credibility check"를 "신용구간 미계산"으로, 벤치마크를 "합성 데이터 sanity check(참 시그니처를 solver에 제공하므로 실제 정확도 근거 아님)"로, 절대 풍부도는 총 UMI 재스케일 값임을 명시. **판단**: `results_pipeline_01` 폴더에는 "Cell2location" 문구가 없음(grep 0건)이고 코드 로직·출력이 바뀌지 않으므로 CSV·결과는 건드리지 않고 재실행도 하지 않음. py 구문 확인만 수행(ast.parse 통과). 벤치마크 CSV에 Spearman·희소 유형 지표를 병기하는 finding 제안은 결과 데이터 변경이 필요해 적용하지 않음(후속 필요 시 sci-checker가 결과 불변 원칙과 조율).
3. **Track 03 SVG 잔존 텍스트 (SCI-03-009/018/035 후속)**: `tracks/03_tme_cell_communication/scripts/generate_tme_diagrams.py`와 SVG 02·03·05를 함께 수정(수정 전 출력은 커밋 SVG 6개와 byte 동일). SVG02 "로버스트 순위 집합(RRA/Borda)"→"(RRA)", "• Borda Count & RRA 집계 알고리즘"→"• Robust Rank Aggregation(RRA) 집계", 주석 정리; SVG03 "PageRank / 정보 중심성"→"정보(Information) 중심성", "Betweenness"→"Flow betweenness"; SVG05 "Intrinsic 45% / Juxta 20% / Para 35% / Spatial Gain +55%"→"(예시)"/수치 제거. 변경 줄 수: 02 2행, 03 2행, 05 4행, 01·04·06은 불변. 미처리: `scripts/generate_tme_slides.py`(슬라이드 생성기)에도 "Borda"·"PageRank" 서술이 남아 있음(L312-314, 398, 430-434, 458, 494, 780) — 슬라이드는 이번 범위 밖이라 수정하지 않음(sci/오케스트레이터 후속).

### 슬라이드 잔존분 (Track 05, Track 03 덱)
- **Track 03** (`tracks/03_tme_cell_communication/slides/01,02_*.html` + `scripts/generate_tme_slides.py`): 수정 전 생성기 출력(old SVG 사용)과 커밋 덱이 byte 동일하지 않았음(덱 01: 3곳, 덱 02: AUROC 0.88/0.79 등 이전 단계의 수작업 패치가 덱에만 반영됨) → HTML을 직접 패치하고 생성기에도 같은 문구 치환을 반영. 변경: LIANA "Borda 순위 집계"→"RRA"(Borda는 일반 이론 예시로 주석), CellChat 중심성 "PageRank"→"정보(Information) 중심성", "Betweenness"→"Flow betweenness", 덱 내 SVG02/03/05 사본의 같은 문구와 SVG05의 "Intrinsic 45%/Juxta 20%/Para 35%/+55%"→"(예시)" 치환, 덱 01 SVG02의 "일치도 75% 이상 향상"→"방법 간 편차 완화 (개념)", 덱 02 SVG05 블록은 현재 docs/images SVG05로 교체(주석 문구 포함). 변경 줄 수: 덱01 29, 덱02 13+SVG05 블록, 생성기 21. **유지한 것**: 이 트랙 파이프라인 출력을 그대로 옮긴 문구(터미널 로그 "Borda consensus/PageRank", "myCAF PageRank = 0.281")와 NicheNet의 PageRank RWR(정확한 서술)은 파이프라인 결과 불변 원칙에 따라 변경하지 않음.
- **Track 03 불일치 보고**: 덱 01의 SVG01(">1mm")·덱 02의 SVG04(AUROC 0.88/0.79)는 `docs/images`/생성기(">1cm", 0.96/0.92)보다 덱이 더 최신이라 SVG 블록 전체 교체를 하지 않음(덱의 이전 수정이 되돌려짐). docs/images·생성기에도 같은 값을 역반영해야 하나 재빌드와 겹칠 수 있어 보류 → 후속 필요.
- **Track 05** (`tracks/05_generative_foundation_models/slides/01,02_*.html` + `slides/build_decks.py`): 생성기 출력(old SVG)과 덱이 byte 동일하지 않음(덱 01 수식 2줄, 덱 02 SVG04 Bin/컨텍스트 길이 등 이전 수작업 패치) → HTML 직접 패치 + 생성기에 동일 치환. 변경: 덱 내 SVG01·03(덱 01), SVG05(덱 02) 사본을 현재 docs/images SVG로 교체(수정 전 사본은 old docs SVG와 byte 동일했음), 텍스트 12곳(z 조건부 독립 보장→변동 감소·보장 아님, MGM→GEP/GEPC 및 MGM은 Geneformer 참고 표기, "완벽 복원/압축/회전", "제로샷 SOTA 달성" 삭제·벤치마크 병기, "p<0.001" 삭제와 가상 값 표기, "최적 표적 조합 도출"→가설 생성). 변경 줄 수: 덱01 18, 덱02 40, 생성기 24. **불일치 보고**: 덱 02의 SVG04("Bin 1~51", 컨텍스트 2,048)가 `docs/images`/`generate_generative_diagrams.py`(Bin 1~32, 30,000 토큰)보다 최신이라 교체하지 않음 → docs/images·생성기 역반영 필요. SVG05 HTML 주석 "Schematic: MGM and GEP"은 보이지 않는 주석이라 남김.

## S3 단계 / ENG (18건)
상태: fixed 11 / softened 3 / already-fixed 2 / deferred 2. 검증: 프로젝트 전체 복사본(scratch)에서 루트+트랙 01~10 빌드 11/11 통과(`--skip-pipelines`), 결과 데이터 불변(Track 09 보고서 텍스트 외 변경 없음). 라이브 Book/HTML/DOCX는 재빌드 대상(미수정).
- ENG-001 fixed: 07 ch06 링크 `../../06_manifold…`; Book/HTML용 변환은 `build_statistics_master.py read_doc_files()`.
- ENG-002 fixed: `file:///` 16건 → 상대경로(28개 링크 대상 실재 확인); 루트·트랙01 빌드의 `read_doc_files()`가 `0N_*.md`→`docs/0N_*.md`로 변환해 #chapter 링크 유지.
- ENG-003 fixed: Track 10 ch03 코드펜스 8개 python→snakemake.
- ENG-015 already-fixed (S2 SCI-08-012), ENG-016 already-fixed (title= 인자 이미 없음).
- ENG-023 fixed: `***X* (…)**` → `**_X_ (…)**` (04 ch06 ×4, 09 ch05 ×1).
- ENG-025 fixed: 절대경로 상수 15개 → `__file__` 기준.
- ENG-026 fixed: 02·03·04에 `--skip-pipelines` 추가, 05~07에 별칭 추가(08~10은 기본 비실행).
- ENG-027 deferred(사용자 결정): `__pycache__` 25개 디렉터리, `./.DS_Store`, `generate_diagrams.py` 루트 vs `scripts/` 중복(내용 상이), 빈 `infercnv_out/`, `backups/` 전체. `.gitignore`만 신규 추가.
- ENG-034 softened: 벤치마크 변동성 경고를 Track 08·10 결과 폴더 NOTICE에 추가(코드/결과 불변).
- ENG-035 fixed: Track 09 파이프라인 01 보고서에 truth 대비 오차 열, "Optimal" 삭제; scRNAseq_anal 재실행으로 수치 동일 확인, 시간값은 기존 유지. (`pipelines/01_fastq_demux_and_umi_dedup_pipeline.py`, `results_pipeline_01/pipeline_01_summary_report.md`)
- ENG-037 softened: 실행 시간 경고+빌드 분리 방법을 Track 05 results_pipeline_02 NOTICE에 기록.
- ENG-042 softened: 11개 DOCX 표지에 "SVG 이미지 미포함" 안내 추가(이미지 삽입은 미구현).
- ENG-043 fixed: DOCX 수식 변환기(11개) 분수·첨자 괄호 보존.
- ENG-044 deferred(사용자 승인): 트랙 간 UI/JS 공통 템플릿화(대규모 리팩터).
- ENG-045 fixed: 검색 스니펫 이스케이프(루트·트랙01).
- ENG-046 fixed: 루트 docs/08 코드펜스 내 `$…$` 제거. 같은 유형이 Track 03 ch03·ch04에 남음(Track 03 S3 대상).
- ENG-047 fixed: 루트 CSS에 640px 이하 TOC 스크롤·인쇄 버튼 숨김. 트랙 02~10 템플릿의 동일 문제 여부는 별도 확인 필요.

## S3 단계 / Track 02 (47건)
상태: fixed 29 / softened 17 / already-fixed 1 (SCI-02-041) / deferred 0. 편집 대상은 현재 파일(S1·S2 sci 수정 반영본) 기준이며 되돌린 변경 없음. 검증: Track 02 docs 코드펜스 짝/파이썬 블록 compile 통과, 전체 복사본 scratch 빌드 8/8 assertion 통과, 코드 수정은 합성 데이터로 실행(squidpy 컬럼명, 스팟 피치 µm 변환). 파이프라인 결과 데이터는 불변(코드 주석/NOTICE만 변경).
- ch00(003·005·006·007): VPT 정정(Vizgen/vizgen-postprocessing 확인), Visium HD 표현, 스팟당 세포 수 범위, HIF1A→HIF 표적 유전자.
- ch01(010~016): poly-dT/프로브 구분, 8 µm 빈, Plexity 식 삭제, 확산 계산 범위 명시, FF 단편화 열, FFPE 확산, 스팟 피치 µm 변환 코드.
- ch02(017~019): Points O(N), R-Tree 구분, SpatialData 변환 목록(소스로 확인).
- ch03(022~024): Visium 플랫폼 수정, 마진 e^-2 잔여, 공간 오염 보정 도구.
- ch04(027~029): 출처 없는 할당율 수치 예시화, z중첩 원인, 377-plex, Baysor 입력.
- ch05(034~036): q05 서술, DNPC/myCAF 마커, RCTD 연도.
- ch06(038~045): SPARK-X 정식화, Geary 전역 통계, Moran 산점도 라벨, squidpy pval_sim·grid, 수치 예시화.
- ch07(048~050): Schürch 2020 귀속과 k, 반경/k 통일, 데노수맙·Galunisertib 서술.
- ch08(055~058): COMMOT 실제 목적함수 요지(소스 확인), 거리·벡터장 해석 완화, 용골성/골경화성 표기.
- ch09(059~064): 병합 임계/VIF 환산, cell2location 사전분포, CosMx 결찰 제외, 거리·임계값 예시화, RL 한계.
- SVG(065 외): `scripts/generate_spatial_diagrams.py`와 SVG03(3행)·04(7행)를 함께 수정, 수정 전 출력은 커밋 SVG와 byte 동일했고 나머지 4개 불변.
- 파이프라인(067·068): 01 docstring/주석, 02 주석, results_pipeline_01/02 SIMULATION_NOTICE에 한계 명시(수치·CSV 불변, 재실행 없음).
- 다른 트랙 동일 오류(note만): Track 03 ch03·ch04 코드펜스 `$…$`.

## S3 단계 / Track 10 (25건)
상태: fixed 17 / softened 7 / already-fixed 1 (SCI-10-042) / deferred 0. 검증: docs 00~09 코드펜스·python 블록 compile 통과, 전체 복사본 scratch 빌드 8/8 assertion 통과, 파이프라인 결과 불변(미수정·미실행). 코드/도구 확인: MLflow 3.17.1 레지스트리 등록+alias 실행, Ensembl REST GRCh38 좌표, Snakemake 9.27.0 --help, nf-prov README. Nextflow 미설치로 02장 combine/branch는 공식 문서 기준(실행 미확인).
- 00장(004·005·007·008), 01장(011·012), 02장(049), 03장(017·018), 04장(024·025·026), 05장(029·030), 06장(032·033·034), 07장(035·036), 08장(025·040·041·043·044), 09장(048): 각 항목 note 참조.
- 정량·가격·좌표 수치는 확인한 것(Ensembl 좌표, 코드 동작)만 확정 표기, 나머지는 "예시/가정"으로 완화.

## S3 단계 / Track 08 (20건)
상태: fixed 12 / softened 7 / already-fixed 1 (SCI-08-042) / deferred 0. 검증: docs 00~09 코드펜스 짝·python 블록 compile 통과, 프로젝트 전체 복사본 scratch 빌드(basic_research py3.8 + mistune 2.0.5, `SCRNA_ENV_SITE_PACKAGES` 지정) 8/8 assertion 통과, 생성물 미수정.
- ch00(034·035·040): Census LTS 2023-12-15 직접 쿼리 수치, 전립선 아틀라스 Henry 2018(PubMed 확인), 크로마틴 완전 억제 완화.
- ch01(023): CSC 열 슬라이싱을 "전체 indices 선형 스캔 O(nnz)"로 재서술, scipy 실측 반영.
- ch02(036·037): backed 열기 시간(1M obs 약 1.4초 실측), `anndata.abc.CSRDataset/CSCDataset`, 속성은 HDF5 문자열/배열.
- ch03(038): 비트 셔플·Blosc 수치/표현 완화. ch04(033·042): 딕셔너리 80% 완화, 버전 조합 주석 구체화.
- ch05(024): QC 병합을 (n, mean, M2) Chan 병합으로 교체(scratch에서 np.var와 일치 확인).
- ch06(025·026): 복잡도 O(N·G·m), Halko 상한 계산값 명시. ch07(027~029·041): 하드웨어 수치·워프 셔플·Hogwild·스택 정정.
- ch08(030~032): CAGRA를 cuVS API로 이전(소스 확인, GPU 부재로 미실행), 브루트포스·정점 색칠 서술 완화. ch09(039): 실패양식 6 재정의.
- 한계: GPU 코드(07·08장)는 실행하지 못함. 새 절 추가 없음.

## S3 단계 / Track 05 (28건)
상태: fixed 20 / softened 8 / already-fixed 0 / deferred 0. 검증: docs 00~09 코드펜스 짝·python 블록 compile 통과, 수정한 코드(01 NB 손실, 02 compute_elbo, 06 Dataset/모델, 09 scIB)는 scratch에서 실행 확인, SVG 3개는 생성 스크립트와 동시 수정 후 재생성 결과가 byte 일치(02·04·05), 전체 복사본 scratch 빌드(basic_research py3.8 + mistune 2.0.5, `--skip-pipelines`) assertion 통과, 파이프라인·결과 데이터 불변. 새 절 추가 없음.
- ch00(005·006·007·012): 코호트 수치(PubMed 확인), 챕터 맵·표 10개 챕터, 의존성(scvi-tools Python 요구·transformers는 Geneformer용), scVI 라이브러리 크기 주석.
- ch01(010·011·012): 증폭 편향/ERCC 서술 완화, NB 손실 log1p 반영과 float32 한계 실측 명시.
- ch02(012·013·016·017): KL 워밍업 기본 400, Wilcoxon 타이 처리·의사반복, library_prior 텐서 정의(실행 확인). SVG02 수정.
- ch03(019~022): scANVI 구현 대비 단순화 명시, Top-M 근사·엔트로피 OOD·희귀 아형 비율 서술 완화.
- ch04(026·027): 분산 매핑 범위, 병합 프로토콜은 가설, U_query 설명 정정.
- ch05~07(028·030·031·035·037·038·039·043): 저자명, 위치 임베딩/구간 수 통일(SVG04), post-norm·층 수 통일(SVG05), Dataset 구간화·빈 세포 NaN, Kendall 계수, Geneformer V1/V2·코퍼스 한계.
- ch08~09(047·049·050·051·052): 치환 p-값 영분포, ISP 방식, scIB 간이 점수(NMI 추가), FlashAttention 수치, AMP/GEP 손실 서술.
- 잔여(생성물): 슬라이드 deck(`slides/*.html`)에 SVG02/04/05의 이전 문구 복사본이 남아 있음(수정 안 함, orchestrator 재생성/패치 대상). 07장 외 다른 장의 스크립트 불변.
- SCI-06-054(sci에서 이관): Track 06 SVG05 "k = 4 거시상태" → "k = 6 거시상태 (n_states=6 지정, 예시)"(`generate_trajectory_diagrams.py` 동시 수정, 재생성 byte 일치, Track 06 빌드 assertion 통과). SVG03 "max_iter=10"은 이미 반영됨. CSV의 SCI-06-054 행은 지시에 따라 수정하지 않았으므로 병합 시 fixed로 반영 필요. 참고: 06장 본문 ASCII 도식(L104)에 "상위 4개 고유값 → k = 4 거시상태 확정"이 남아 있어 코드 n_states=6(말단 4 + 초기 1 + 중간 1)과 어긋남 → sci 확인 필요.

## 결과 재생성 (2026-10-20)
방법: scRNAseq_anal(py3.8.16)로 scratch 사본에서 실행. 수정 전 코드 3개를 재실행한 결과는 커밋된 결과와 byte 동일(기준선 OK)했고, 수정 후 diff는 아래 파일·필드에만 발생. 커밋하지 않음(작업트리 변경 11개 파일).
- SCI-03-042: Track 03 pipeline 01 요약의 top_sender를 out-degree 최댓값으로, top_receiver(in-degree)를 추가. top_sender TAM_M2 9.3096 → myCAF 11.5168(CSV와 일치), top_receiver Endothelial 7.3383 신규. 중심성 CSV·나머지 불변. 변경 파일: pipeline_summary_metrics.json, SIMULATION_NOTICE.md.
- SCI-02-067: Track 02 pipeline 01 Var(I)를 Cliff&Ord 정규성 가정 정확식(S0,S1,S2)으로 교체, docstring 수정. 검증(50x50, k=6 행표준화): exact sd 0.01104(구 0.00816), 정규 난수 경험 sd 0.01127, 순열 sd 0.01112. 지시서의 0.0113/0.0118과 약간 다르나 정확식과 경험값이 일치. 결과: morans_i 불변, z_score만 1.3797배 감소(ALPL 67.57→48.97, SPP1 65.78→47.68, KLK3 약 53.2→38.58), p/FDR 대부분 0.0 유지, is_svg 367개 불변. 변경 파일: spatial_variable_genes_moran.csv, SIMULATION_NOTICE.md.
- SCI-02-068: Track 02 pipeline 02 니치 라벨을 `Niche_<k> (<우세 세포유형>-dominant)`로 교체(세포 구성 불변, 생물학적 이름 삭제), 요약 JSON에 limitations(RANKL Tumor→Tumor 79.4%, 파골세포 유형 없음 등) 추가. 대응: Bone_Front→Niche_3 (Osteoblast-dominant), Cancer_Core→Niche_1 (Tumor_AR_Pos-dominant), Desmoplastic_Rim→Niche_0 (CAF-dominant), TLS_Immune_Hub→Niche_2 (CD8_T-dominant). 수치(비율·세포 수·신호 흐름 CSV)는 불변. 변경 파일: niche_composition_summary.csv(행 순서), xenium_cell_metadata.csv(spatial_niche 열), pipeline_02_execution_summary.json, SIMULATION_NOTICE.md.
- docs 인용(sci 수정 대상): 02/docs/06_SPATIALLY_VARIABLE_GENES.md L178(KLK3 I=0.565, 불변), L92(Var(I)·Z 설명); 02/docs/07 및 08의 니치 이름은 결과 이름을 직접 인용하는 곳 없음(grep: 07 L242 `niche_composition` obsm만); 03/docs/04 L314 중심성 서술(top sender 직접 인용 없음). z_score 값이나 TAM_M2 top sender를 직접 인용한 docs 위치는 발견되지 않음.

## ENG-040 슬라이드 레이아웃 (2026-10-20)
- 추가: `scripts/inject_slide_layout.py`(마커 `SLIDE-LAYOUT:BEGIN/END` css+js, 멱등, `__file__` 기준 경로, 인자 없으면 덱 23개 전체). 덱 23개에 주입(텍스트·SVG 불변: 원본 대비 diff에 삭제 줄 0, strip(주입본)==원본 23/23). `.sl-fit` 래퍼는 런타임 JS가 만든다.
- 동작: 데스크톱은 `.slide-content` 높이를 100vh−84px로 넓히고 내용을 zoom(최소 0.7)으로 맞춘 뒤 안 맞으면 기존 스크롤로 폴백. 모바일(≤768px)은 슬라이드 세로 스크롤 모드, 그리드 1열, 표 가로 스크롤, 코드·doc-ref 줄바꿈.
- 빌드 순서: `requirements-build.txt` 주석 1줄, `scripts/precommit_checks.py`에 WARN(slides/*.html에 마커 없음) 추가, `--all` OK.
- 측정(23덱×8장=184장, 토글 열림/닫힘, "스크롤 필요" 슬라이드 수): 1280×720 열림 145→72, 닫힘 137→9; 1920×1080 열림 52→5, 닫힘 1→0; 375×812 열림·닫힘 184→0(가로 넘침 슬라이드 182→0). 최소 zoom 0.6 대비 0.7은 1280×720 열림 29→72, 닫힘 7→9로 스크롤 슬라이드가 늘어남.
- 멱등성: 두 번 실행 후 23개 파일 해시 동일. 키보드(→←)·이전/다음 버튼·카운터는 전 덱 8장 순회 및 표본 덱에서 정상.
- 남은 문제: 1280×720 열림 상태 72장은 여전히 내부 스크롤(기본이 details open이고 SVG가 큼). 필요하면 드로어 SVG max-height를 더 낮추는 옵션이 있음. 모바일에서 표·pre는 내부 가로 스크롤.

## ENG-044 공통 UI 템플릿 (2026-10-20)
- 신규 `scripts/course_ui_template.py`(`render_page`: HTML 셸·CSS·JS 한 벌, 08/09 기준 + 00/01 표지·TOC·검색 팝업·사이드바 배지 + 10의 details 스타일). 빌더 11개의 `full_html = f"""..."""`(트랙당 685~860줄)를 `render_page(...)` 호출로 교체(빌더 합계 −8,076/+165줄). 빌더 assertion과 스위처 앵커(`</style>`, `search-container`, `document-container`/`main#mainContent`)는 그대로.
- 기능 상향: 검색 하이라이트+개수+결과 팝업(전 트랙), 퀴즈 펼치기/접기·인쇄 버튼·테마 localStorage 유지(00·01 포함), `regex.test` /g 버그 수정(누락 매치 복구). 트랙별 theme 저장 키는 기존 값 유지.
- 검증: 11개 빌드 assertion 9/9, 스위처 11/11, verify_switcher_links OK. `*_Master_Book.md` 11개 변경 없음(byte 동일), DOCX 11개 zip 파트 동일(HEAD로 되돌림), HTML 11개 `<main>` 정규화 텍스트·SVG 동일, 파이프라인 결과 변경 없음. Playwright 11개(1280×800, 375×812): 검색·팝업·다크모드 유지·퀴즈 펼침/접힘·인쇄·스위처 이동·모바일 가로 스크롤 없음·콘솔 에러 0 모두 통과. precommit --all OK.
- 부작용: 00·01(마크업 변경) 사이드바 헤더는 배지가 위 줄, 제목, 테마 버튼 순으로 배치되도록 CSS 보정. HTML 크기 +3~6 KB씩. 변환 로직 통합은 범위 제외.

## 참고문헌 체계 (2026-10-20)
- 도구: `scripts/format_references.py`(PMID: NCBI esummary, DOI: Crossref, 캐시 `review/references_cache.json`), `scripts/apply_references.py`(`review/refs/*.ids` → 챕터 끝 REFERENCES 블록, 멱등, `--check`), `review/software_refs.json`(107개 키), precommit 연결(블록≠.ids 이면 BLOCK, 블록 없으면 WARN). 빌더 11개 DOCX 변환에 HTML 주석 줄 건너뛰기 1줄 추가(마커 누출 방지).
- NEED-SW 해소: sci `.ids`의 64줄 처리(매핑: sklearn→scikit-learn, tabix·bgzip→htslib, 나머지 기존/신규 키), 신규 키 27개(PMID는 조회로 확인, 원 논문이 없는 도구는 paper 비움, source는 HTTP 200 확인; Cell Ranger는 10x 서버가 429로 응답).
- 전체 적용: 챕터 109개 모두 블록 삽입(.ids 109개), 고유 ID 210개 조회 실패 0, `--check` 109/109 일치, 블록 밖 본문은 HEAD와 동일(109/109). Crossref에 container-title이 없는 항목이 `**`로 렌더링되어 빌드 assertion이 막히던 문제를 포매터에서 수정(2개 트랙에서 발생).
- 재빌드: 11개 트랙 assertion 9/9, 스위처 11/11, verify_switcher_links OK. HTML 참고문헌 영역 `<a>`·`https` 0, DOCX 마커 누출 0, Book.md 블록 109개. precommit --all OK, 파이프라인 결과 변경 없음.
- 한계: Kendall 2018 CVPR는 Crossref 저자 순서(Cipolla R, Gal Y, Kendall A), Ross IPCA는 연도 2007로 표기(API 값 그대로).

## Converter dedup (2026-10-21)
- 신규 `scripts/course_markdown.py`(1,644줄): LaTeX 변환(parse_balanced_braces, LATEX_SYMBOLS 합집합 196키, convert_latex_to_html/unicode(T10 기준 상위집합), replace_math_*), `build_master_markdown/html/docx(cfg, contents)`, 검증 `verify_track`(어서션 1-3·5-8 공통, 4는 트랙별 주입)·`verify_basic`(Track 00/01 프로필)·`assert4_result_dirs/summary_files/summary_report`. UI 셸은 기존 `course_ui_template.py`. `TrackConfig`로 트랙 설정만 전달.
- 빌더 11개는 설정(CFG) + 얇은 래퍼 + 트랙별 파이프라인 어서션 + main으로 축소(파일명·CLI·`--skip-pipelines`·`--re-run-pipelines`·산출물 경로 유지, mistune 2.0.5 가드·DOCX `<!--` 건너뛰기·REFERENCES/TRACK-SWITCHER 마커·assertion 5 유지). 줄 수: 빌더 합계 16,261 → 2,030, 공통 모듈 +1,644(순 −12,587줄).
- LATEX_SYMBOLS 결정: `\parallel` → ∥(U+2225), `\quad` → 공백 2칸. T10에 없던 14키(`\binom` 등, 값은 빈 문자열)는 합집합에 포함.
- 변환기 비교(T00~T10 각 변종 vs T10): 이전 트랙에서 빠져 있던 줄은 `\mathbb{P}`(T03/T04가 'P'로 둔 부분은 T10의 ℙ), `\mathcal` 맵 확장, `\mathrm/\operatorname/\mathbf/\boldsymbol/\boxed/\overline/\dot/\pmatrix/\cases/\dfrac`와 `_*` 첨자 처리 등 상위집합 처리이며, 역방향(T10에서 제거된 동작)은 T00/T01의 `_` 첨자 정규식(`*` 미포함) 1줄뿐.
- 검증(HEAD 대비, `--skip-pipelines`, basic_research env + SCRNA_ENV_SITE_PACKAGES): 11개 빌드 assertion 9/9 통과, `inject_global_track_switcher.py` 11/11, `verify_switcher_links.py` ALL VERIFIED, `precommit_checks.py --all` OK, `results_pipeline_*` 변경 없음. Master_Book.md 11개 byte 동일. HTML: T00·T06~T10 6개 byte 동일, T01~T05 5개 변경. DOCX: T00·T06~T10 6개 zip 파트 동일, T01~T05 5개의 `word/document.xml`만 변경(나머지 파트 동일).
- 비동일 산출물(HEAD → 현재), 원인은 모두 변환기 상위집합 또는 두 기호 결정. 건수는 diff hunk 근사치(한 hunk가 여러 원인에 걸칠 수 있음):

| 원인 | 트랙: html/docx hunk | 예시 (HEAD → 현재) |
|---|---|---|
| `\quad` 1칸→2칸 | T01: 6/5, T02: 13/10 | `TSS Enrichment Score ≥ 6.0   AND` → `... ≥ 6.0    AND` |
| `\overline`/`\bar` 결합 막대 렌더 | T01: 2/1, T02: 16/10, T03: 1/1 | `FPS={flank^(+)}-{core^(-)}` → `FPS=flank^(+)̄-core^(-)̄` |
| `\partial`/`\nabla` 기호 매핑 | T02: 12/5, T05: 8/5 | `partial V / nabla^2 V` → `∂ V / ∇^2 V` |
| `\operatorname`/`\mathrm`/`\text` 그룹 렌더(리터럴 `{}` 제거) | T02: 25/23, T03: 2/2, T04: 5/3, T05: 1/1 | `{Var}(Y)` → `<span class="math-op">Var</span>(Y)` |
| `\dot{N}` 렌더 | T04: 5/3 | `{N}` → `Ṅ` |
| `\boldsymbol`/`\mathbf` 렌더 | T02: 1/0, T05: 3/0 (docx는 텍스트 변환이라 동일) | `{Σ}_k` → `<strong>Σ</strong>_k` (html), docx `Σ_k` |
| `\parallel` U+2016→U+2225 (결정) | T05: 9/8 | `D_KL(q ‖ p)` → `D_KL(q ∥ p)` |
| 남은 `{…}` 그룹·첨자 처리 | T01: 1/0, T02: 4/1, T03: 4/1, T05: 7/3 | `{(…)}_Hill` → `(…)_(Hill)` |

- 명령: `$PY scripts/build_master_documents.py --skip-pipelines`, `$PY tracks/NN_*/scripts/build_*_master.py --skip-pipelines`(11개), `$PY scripts/inject_global_track_switcher.py`, `$PY scripts/verify_switcher_links.py`, `python3 scripts/precommit_checks.py --all`; assert4 헬퍼는 subprocess를 모킹해 기본/`--skip-pipelines`/`--re-run-pipelines` 동작을 확인(실제 파이프라인은 실행하지 않음).
- 참고: 작업트리에 `.hgp/subagent-summary-*.json`(훅 산출물, 이번 작업과 무관)이 untracked로 존재.
- 수정(다중 문자 `\overline`/`\bar`): 인자가 여러 글자인데 결합 문자 U+0304 하나를 그룹 뒤에 붙여 막대가 마지막 글자에 얹히거나 `</span>` 뒤에 매달리던 결함을 `course_markdown.py`에서 고침. HTML: 인자가 한 글자(또는 한 글자로 매핑되는 `\gamma` 등)이면 기존 결합 문자 유지(x̄·v̄·d̄·γ̄ byte 동일), 그 외에는 `<span style="text-decoration:overline">…</span>`(인라인 스타일, CSS 불변). DOCX: 여러 글자 인자는 변환된 텍스트의 보이는 각 문자 뒤에 U+0305, 한 글자는 기존 U+0304. 재빌드(11개 9/9, 스위처 11/11, verify_switcher_links OK, precommit OK, results 불변): Book 11개 byte 동일, T00·T07~T10 HTML/DOCX 동일. **T06은 HEAD와 달라짐**: 06장의 `\overline{\Delta x_{ij}}`가 HEAD에서도 `Δ x_ij̄`(마지막 글자에만 막대/매달림)였고 이번에 올바른 막대로 바뀜(아래 표). 이전 보고의 "T06 동일"은 이 수정으로 대체됨.
- `\overline`/`\bar` 전·후 (T01–T03, T06 전수; 이전 표의 "T03 overline"은 `\mathrm{mean}` 오분류이며 실제 overline hunk 아님):
  - T01 HTML(2): `FPS={flank}-{core}`(리터럴 중괄호) → `FPS=<span style="text-decoration:overline"><span class="math-text">flank</span><sup>+</sup></span>-<span style="text-decoration:overline">…core…<sup>-</sup></span>`. T01 DOCX(1): `FPS={flank^(+)}-{core^(-)}` → `FPS=f̅l̅a̅n̅k̅^̅(̅+̅)̅-c̅o̅r̅e̅^̅(̅-̅)̅`.
  - T02 HTML(16): FDR_probe 4곳 `{<span class="math-text">Counts</span>}<sub>Blank</sub>` → `<span style="text-decoration:overline"><span class="math-text">Counts</span></span><sub>…Blank</sub>`(막대는 Counts만, 첨자는 밖); 한 글자 12곳 `{x}` → `x̄`. T02 DOCX(10): `Counts_(Blank)` 2곳 → `C̅o̅u̅n̅t̅s̅_Blank`; 한 글자 `{x}` → `x̄` 8곳(`x_i - x̄`, `z = x - x̄` 등).
  - T03 HTML/DOCX: overline 해당 없음.
  - T06 HTML(4): `Δ x<sub>ij</sub>̄)`(매달림) → `<span style="text-decoration:overline">Δ x<sub>ij</sub></span>)`(2곳 ×2). T06 DOCX(1): `Δ x_ij̄` → `Δ̅ x̅_̅i̅j̅`(밑줄 문자에도 막대가 붙는 점은 "보이는 모든 문자" 규칙의 부작용).

## Converter defect fixes (2026-10-21)
`scripts/course_markdown.py`만 수정(d355177 이후). review/converter_dedup_review.md의 결함 5건:
1. `\|`, `\left\|`, `\right\|` → ∥(U+2225), HTML·Unicode 모두. `\Vert`는 의도적으로 ‖ 유지(T10 06·08장이 사용하므로 지시 범위 밖). `\left\|`가 먼저 매칭되도록 길이순 치환 순서에 의존(확인함).
2. `\vec{x}` → 인자 + U+20D7(HTML·Unicode 동일). 여러 글자 인자는 그룹 뒤에 화살표 하나(마지막 글자 위). 선택 근거: CSS 화살표 위 장식이 없고, 유니코드 글자별 화살표는 가독성이 나쁨.
3. `\tfrac`도 `\dfrac`처럼 `\frac`으로 치환(`\dfrac`는 이미 처리되고 있었음).
4. `\text{}`/`\mathrm{}`/`\operatorname{}` 안의 `\_`를 센티넬(U+E000)로 보호하고 최외곽 호출에서만 `_`로 복원(`convert_latex_to_html/unicode`를 깊이 카운터 래퍼로 감쌈; 재귀 호출에서 하위 첨자 규칙에 걸리지 않게 함).
5. DOCX 다중문자 overline: `_ ^ ( )`에는 U+0305를 붙이지 않음.
- 재빌드(11개 `--skip-pipelines`, 9/9), `inject_global_track_switcher.py` 11/11, `verify_switcher_links.py` ALL VERIFIED, `precommit_checks.py --all` OK, `results_pipeline_*` 불변, Master_Book.md 11개는 d355177과 byte 동일.
- d355177 대비 HTML/DOCX hunk(전부 위 1–5로 귀속, 그 밖의 변경 0; 건수는 토큰 diff hunk):

| 트랙 | HTML | DOCX | 예시 (이전 → 이후) |
|---|---|---|---|
| T00 | #4: 3 | – | `<span class="math-text">total<sub>c</sub>ounts</span>` → `total_counts` |
| T01 | #4: 3 | #5: 1 | `n<sub>f</sub>ragments` → `n_fragments`; DOCX `FPS=f̅l̅a̅n̅k̅^̅(̅+̅)̅-…` → `^`·`(`·`)` 제외 |
| T02 | #1: 20, #2: 12, #4: 3 | #1: 11, #2: 12 | `{v}(x,y) = … \|∇V\|` → `v⃗(x,y) = … ∥∇V∥`; `cells<sub>p</sub>er<sub>l</sub>ocation` → `cells_per_location`; `‖ M[…]`(`\left\|`) → `∥ M[…]` |
| T03 | #1: 8, #3: 2 | #1: 4, #3: 1 | `d_ij = \| s_i - s_j \|_2` → `∥ s_i - s_j ∥_2`; `[× tfrac{n_in_j}{n^2}]` → 분수 렌더 |
| T04 | #1: 4 | #1: 2 | `\\|C_u - C_v\\|_1` → `∥C_u - C_v∥_1` |
| T05 | #1: 38, #4: 4 | #1: 15 | `KL(q(z\| x) \| 𝒩(0,I))`·`\|z_i - μ_k\|` → `∥`; `STOP<sub>G</sub>RADIENT` → `STOP_GRADIENT`, `to_NEPC` |
| T06 | #1: 43, #4: 1 | #1: 19, #5: 1 | `\|x_i - x_j\|` → `∥x_i - x_j∥`; `fit<sub>l</sub>ikelihood` → `fit_likelihood`; DOCX `Δ̅ x̅_̅i̅j̅` → `Δ̅ x̅_i̅j̅` |
| T07 | #1: 2, #4: 2 | #1: 1 | `patient<sub>i</sub>d` → `patient_id` |
| T08 | #1: 14, #4: 1 | #1: 7 | `min_(rank(A)≤k) \|X - A\|_2` → `∥X - A∥_2`; `counts<sub>p</sub>er<sub>c</sub>ell` → `counts_per_cell` |
| T09 | #4: 7 | – | `gap<sub>o</sub>pen` → `gap_open`, `no_err`, `cigar_int` |
| T10 | #1: 1, #4: 11 | #1: 1 | `cpu.cfs<sub>q</sub>uota<sub>u</sub>s` → `cpu.cfs_quota_us`, `fastq_1`, `oom_score_adj` …; KL `\|` → `∥` |

  (`#3`의 T03 HTML에는 분수 구조 연속 hunk 2개가 더 있음: `}{n<sup>` → `</span><span class="den">n<sup>`, `}]` → `</span></span>]`.) 지시서의 영향 범위(T02·T03·T04·T05) 밖에서도 `\|`(T06·T07·T08·T10)와 `\_`(T00·T01·T05~T10)가 같은 규칙으로 바뀜.
- 브라우저 확인(로컬 `python3 -m http.server`, Playwright 요소 스크린샷, scratch `cd/shots/`): `vec_norm.png`(T02 `v⃗ = ∇V / ∥∇V∥`), `underscore.png`(T01 `log10(n_fragments)`), `tfrac.png`(T03 `[× n_i n_j / n²]` 분수), 그리고 T05 KL `KL(q(z| x) ∥ 𝒩(0,I))`(텍스트 확인).

## Pipeline-mode verification (2026-10-21)
- 목적: 리팩터(d355177/b0be89e) 이후 한 번도 실행하지 않은 빌더의 기본 모드(`--skip-pipelines` 없음, assertion 4 실제 재실행 형태) 검증. 실제 트리에서는 실행하지 않고 `git clone` 사본(HEAD bdcbe5f, `scratchpad/eng/pipeclone`)에서만 실행. 빌더·파이프라인 스크립트에 실제 프로젝트 절대경로 없음(grep), 실행 후 실제 트리의 결과·생성물 변경 없음(`git status`에는 sci 작업분만 보임). 사본에 필요한 추가 입력 없음(빌더는 추적 파일만 읽음).
- 환경(§17과 동일): `PYTHONPATH` 해제, `PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 SCRNA_ENV_SITE_PACKAGES=<conda-base>/envs/scRNAseq_anal/lib/python3.8/site-packages`, 인터프리터 `<conda-base>/envs/basic_research/bin/python`(py3.8.16, numpy 1.23.5, mistune 2.0.5; 파이프라인은 빌더의 `sys.executable`로 실행).
- 명령: 사본에서 `$PY scripts/build_master_documents.py`, `$PY tracks/*/scripts/build_*_master.py`(기본 모드), 이어 Track 08·09·10은 `--re-run-pipelines`(이 세 트랙은 원래 이 플래그가 있어야 재실행), `$PY scripts/inject_global_track_switcher.py`, `$PY scripts/verify_switcher_links.py`.
- 결과(트랙: 종료코드, 어서션, 소요시간): T00 rc0 9/9 0s(파이프라인 어서션 없음 프로필), T01 rc0 9/9 1s(동), T02 9/9 6s(재실행), T03 9/9 4s(재실행), T04 9/9 15s(재실행), T05 9/9 322s(재실행), T06 9/9 18s(재실행), T07 9/9 8s(재실행), T08 기본 9/9 1s(요약 확인만) / `--re-run-pipelines` 9/9 41s, T09 기본 0s / 재실행 9/9 2s, T10 기본 0s / 재실행 9/9 11s. assertion 4는 재실행 형태("executed and exited with Code 0" / "re-executed ...")로 통과. 스위처 11/11, verify_switcher_links ALL VERIFIED.
- 산출물: Master_Book.md 11개·HTML 11개는 사본 HEAD와 byte 동일, DOCX 11개는 git 상 변경으로 보이나 zip 파트는 모두 동일(타임스탬프만). `results_pipeline_*`: T02~T06은 변경 없음(byte 동일, 2·3은 §17 재현). 변경된 결과 파일은 16개(T07 2, T08 7, T09 2, T10 5; 모델 가중치 등 수치 결과 파일은 불변)이며 모두 비결정 필드 때문:
  - T07: `execution_time_seconds` 1줄씩. T09: 실행 시간 문구 1줄씩. T08: 벤치마크 시간·속도비(`col_slice_csc_sec`, "ANN Speedup 25.1x→27.7x" 등). T10: `trace_table.csv` 타임스탬프/지속시간, `serving_benchmark.csv`·`summary_report.md`의 처리량(8,432→8,868 cells/s), `model_audit_manifest.json`의 `python_version` 3.8.18→3.8.16(커밋본은 conda base 3.8.18에서 생성, 이 환경은 3.8.16; 기존 기록과 동일한 알려진 차이). 숫자/타임스탬프를 정규화하면 구조 차이 없음(지수표기 부동소수 2건은 같은 시간 필드).
  - 비결정성 증명: 같은 사본에서 T08 파이프라인 02를 두 번 실행하면 실행 간에도 같은 4개 파일이 달라짐. 리팩터 이전 빌더(`fa0cb97`, 별도 사본)로 T07·T08·T10을 같은 방식으로 실행했을 때 변경된 결과 파일 집합이 리팩터 후와 동일.
- 결론: 리팩터로 인한 결함 없음(실패·구조 차이 0). 결과 폴더의 시간/벤치마크 필드와 python_version은 환경·실행 시각에 따라 변하므로, 기본 모드를 실제 트리에서 돌리면 T07~T10 결과 파일이 바뀐다는 점에 유의(커밋 시 ALLOW_RESULTS_CHANGE 필요).

## Pipeline result reproducibility (2026-10-21)
- 문제: 빌더 기본 모드가 `tracks/*/pipelines/results_pipeline_0X`를 제자리에서 덮어써 시간·타임스탬프·벤치마크 값(16개 파일, T07 2·T08 7·T09 2·T10 5)이 실행마다 바뀌었음. assertion 4는 종료코드만 확인했음.
- 변경(`scripts/course_markdown.py`): 파이프라인 단계가 `pipelines/` 폴더를 임시 디렉터리로 **복사**해 거기서 두 스크립트를 실행(같은 인터프리터·환경). 생성된 `results_pipeline_01/02`를 커밋본과 비교: 파일 집합 동일, 파일은 byte 동일, 단 `scripts/pipeline_volatile_fields.json`에 선언된 필드만 마스킹 후 비교. 불일치 시 파일과 첫 번째 다른 필드를 출력하고 assertion 4 실패. 플래그 의미 유지(T02–T07 기본 실행, T08–T10은 `--re-run-pipelines`, `--skip-pipelines`/`--verify-fast`/`VERIFY_FAST=1` 건너뜀). `--update-results`: 제자리 실행으로 커밋 결과를 덮어씀(의도적 재생성 전용). 임시 디렉터리는 실행 후 `trash`로 정리. 빌더 10개 docstring에 CLI 설명 추가.
- 사전 확인: T02–T10 파이프라인은 모두 `__file__`로 출력 경로를 잡음(T02–T08 `repo_root=dirname(dirname(__file__))`+`pipelines/results_*`, T09·T10 `dirname(__file__)/results_*`)이고 `pipelines/` 밖은 읽지 않음. T06 파이프라인 02가 읽는 01 결과도 복사본 안에서 해결됨. 예외: T02 두 파이프라인이 `MPLCONFIGDIR=/tmp/matplotlib_cache`를 설정(프로젝트 밖 캐시 디렉터리, 결과 무관).
- 휘발 필드 도출: 복사본에서 파이프라인을 3회 반복 실행해 커밋본 및 실행 간 차이를 모은 결과 정확히 16개 파일뿐이며 파일 집합은 동일. 3회×9트랙 모두 선언 마스킹 후 불일치 0.
- 선언된 휘발 필드(전체):
  - T07 p1 `pipeline_01_metrics.json`: `execution_time_seconds`; T07 p2 `pipeline_02_metrics.json`: `execution_time_seconds`.
  - T08 p1 `pipeline_01_metrics.json`: `benchmarks.{row_slice_csr_sec,row_slice_csc_sec,row_slice_speedup_csr_vs_csc,col_slice_csc_sec,col_slice_csr_sec,col_slice_speedup_csc_vs_csr,spmv_runtime_sec,spmv_gflops}`, `streaming.streaming_time_sec`, `incremental_svd.runtime_sec`, `execution_time_total_sec`; `pipeline_01_summary.txt`: Row-slicing/Col-slicing/SpMV Throughput/Streaming Time/Total Pipeline Runtime 줄의 소수 숫자; `sparse_memory_layout_benchmark.csv`: `Value` 열의 Row_Slice_CSR_Time, Row_Slice_CSC_Time, Row_Slice_CSR_Speedup, Col_Slice_CSC_Time, Col_Slice_CSR_Time, Col_Slice_CSC_Speedup, SpMV_Runtime, SpMV_Throughput 행.
  - T08 p2 `pipeline_02_metrics.json`: `metadata_slicing.{row_time_per_query_ms,columnar_time_per_query_ms,query_speedup}`, `ann_search.{brute_time_sec,ann_time_sec,ann_speedup}`, `hardware_profiling.peak_ram_mb`, `execution_time_total_sec`; `pipeline_02_summary.txt`: Multi-attribute query latency/Query retrieval speedup/ANN Speedup/Total Pipeline Execution Time 줄의 소수 숫자; `ann_vs_bruteforce_knn_benchmark.csv`: Brute_Force_Time, ANN_IVF_Time, ANN_Speedup 행의 `Value`; `columnar_vs_row_metadata_benchmark.csv`: Row_Query_Latency, Columnar_Query_Latency, Slicing_Speedup 행의 `Value`.
  - T09 p1 `demux_and_error_correction_summary.json`: `umi_deduplication_benchmarks.execution_time_seconds`; `pipeline_01_summary_report.md`: Deduplication Execution Time 줄.
  - T10 p1 `timeline_metrics.json`: `cold_execution.duration_seconds`, `warm_execution.duration_seconds`; `trace_table.csv`: 열 `duration_s`, `timestamp`. T10 p2 `model_audit_manifest.json`: `timestamp_utc`, `runtime_environment.python_version`(커밋본은 3.8.18, 빌드 환경은 3.8.16 — 커밋본은 재생성하지 않음); `serving_benchmark.csv`: 열 mean/std/p50/p95/p99 latency, `throughput_cells_sec`; `summary_report.md`: Timestamp (UTC) 줄, `| **N cells** |` 표 행의 소수 숫자, "Low-Latency Guarantees" 줄의 p95 ms, "High Throughput" 줄의 cells/second·초 값.
  - 선언하지 않는 것(계속 비교): 정확도·F1·Mean_Recall·메모리 크기·압축비·해시·task id/status·SLA 준수 여부·모델 가중치 등 모든 과학/수치 결과. 각 항목에 `why` 메모 포함.
- 검증(실제 트리, 환경은 §17과 동일): 11개 빌더 기본 모드 전부 9/9(T02–T07은 임시 복사본 재실행·재현 확인 통과, 소요 T02 5s·T03 4s·T04 15s·T05 320s·T06 18s·T07 8s), 이어 T08·T09·T10 `--re-run-pipelines` 9/9(42s·2s·11s). 스위처 11/11, verify_switcher_links ALL VERIFIED, `precommit_checks.py --all` OK. `git status`에 `results_pipeline_*` 변경 없음, Book/HTML 변경 없음, DOCX는 zip 파트 동일(타임스탬프만 달라 `git checkout`으로 되돌림).
- 음성 시험(scratch 복사본): (a) 커밋본 `pipeline_01_metrics.json`의 `naive_type1_error_pct`를 87.75→88.25로 바꾸자 T07 빌드가 실패: `first non-volatile difference at .naive_type1_error_pct: committed='88.25' produced='87.75'`; 원복하면 통과. (b) `execution_time_seconds`만 99.99로 바꾼 경우 통과. `--update-results`는 제자리 덮어쓰기(99.99 → 7.51) 확인.
- 주의: `scripts/pipeline_volatile_fields.json`은 새 파일(검토 대상). 새로운 휘발 필드가 생기면 이 파일에만 추가하며, 과학/수치 결과는 선언 금지.

### 후속(2026-10-21): 휘발 필드 보강(A) · 출력 리터럴 수정(B) · 선언 검증 강화(C)
- A (`scripts/pipeline_volatile_fields.json`, 항목마다 why 명시): 10-02 `serving_benchmark.csv` `sla_compliant` 및 `summary_report.md`의 ✅ PASS/❌ FAIL·SLA 문구; 08-02 `sub_millisecond_pass`·`speedup_target_pass`와 summary.txt의 판정 문구, "Hardware Peak RAM"; 08-01 `peak_ram_mb`(+`strict_budget_pass`, CSV `Peak_RAM_MB`, summary의 Peak RAM·판정); 10-02 manifest `runtime_environment.os_platform`/`processor`; 10-02 "High Throughput" 규칙은 `[\d,]+ cells/second`와 계산된 `~N seconds`만 마스킹.
- B (파이프라인 출력 리터럴을 측정값에서 계산, 해당 값은 휘발로 마스킹):
  | 위치 | before | after |
  |---|---|---|
  | 10-02 summary | `~3.8 seconds`(고정) | `10000/throughput` 계산(커밋본 기준 `~1.2 seconds`) |
  | 10-01 summary | `~1.28s / ~0.02s / 64x speedup`(고정) | cold/warm 실측에서 계산(커밋본 기준 `~0.39s / ~0.00s / >788x speedup`; warm=0이면 `>N x`, 둘 다 불가면 `n/a`) |
  | 08-02 `speedup_target_pass` | 기준과 불일치(`10<=s<=50 or s>10`) | `10.0 <= s <= 50.0`(챕터의 10x–50x 기준), 측정 판정으로 유지 |
  | 08 summaries | `PASSED` 리터럴(Strict budget, R>0.95, Sub-ms, 10x–50x, >90%) | 실제 검사 결과에서 `PASSED/FAILED`·`True/False` 도출 |
  | 08-02 summary | `(Zero memory leak: PASSED)` | `(tracemalloc peak)` |
  | 10-02 SLA 문구 | 항상 "strictly compliant with" | `sla_compliant`에 따라 "strictly compliant with"/"NOT within"(통과 시 문구 동일) |
- 결과 파일 처리: `--update-results`로 T08/T10을 재생성하면 휘발 파일 전체가 새 타이밍으로 바뀌고 python_version(3.8.18→3.8.16)도 달라지므로 검증용으로만 실행한 뒤 `git checkout`으로 되돌리고, B 리터럴만 커밋본 측정값에서 계산해 최소 수정함(변경 3개: 08-02 summary.txt, 10-01/10-02 summary_report.md).
- C (`scripts/course_markdown.py`): 선언된 JSON 키·CSV 열·CSV 행 레이블·line 패턴이 커밋본/생성본에 없으면 명확한 메시지로 실패; CSV 헤더·행 수 불일치도 실패(zip 절단 없음).
- 시험(scratch 복사본, T10 `--re-run-pipelines`): 기준선 통과; (a) 비휘발 `cache_hit_percentage` 100.0→99.0 → 실패 ("first non-volatile difference at .warm_execution.cache_hit_percentage"); (b) p95·sla_compliant·timestamp_utc·os_platform만 변경 → 통과; (c) 선언된 `timestamp_utc` 키 이름 변경 → 실패 ("declared volatile JSON key 'timestamp_utc' does not exist in the committed file"); (c2) CSV 열 `sla_compliant`→`sla_ok` → 실패(헤더 불일치); 복원 후 기준선 통과. CPU 부하(20개 busy loop/10 CPU) 하 10-02 재실행도 통과(exit 0, Assertion 4 PASSED).
- 실제 트리 인수: 11개 빌더 기본 모드 9/9, T08–T10 `--re-run-pipelines` 9/9, 스위처 11/11, 링크 검사 OK, `precommit_checks.py --all` OK, results 변경은 위 3개 파일뿐(두 번째 실행 후 추가 변경 없음), 타임스탬프만 달라진 DOCX는 `git checkout`으로 되돌림.

### 후속 2(2026-10-21): \underbrace 렌더링 · T10 manifest 환경 갱신 · 이식성 안내
- `\underbrace{X}_{L}` (`scripts/course_markdown.py`, HTML/유니코드 두 경로): 기존에는 `X` 뒤 아래첨자 `(L)`로 납작하게 처리됨. HTML은 인라인 스타일만 사용(새 CSS 클래스 없음)해 inline-flex 세로 배치로 내용 / 아래 괄호(border 하단 곡선) / 가운데 정렬 라벨을 출력. 유니코드·DOCX는 `[X]_(L)` 형식으로 통일(내용의 첨자·분수 유지). `\Vert`(‖)는 변경 없음. 사용처 3곳만 변경: T03 ch04(Hill), T05 ch02(ELBO 재구성·KL), T07 ch04(Laplace).
  - 변경 전 → 후: T03 `…_(Hill)` → `[((L_i R_j)/(K_h+L_i R_j))]_(Hill)`; T05 `…]_(재구성 우도…) - D_KL(…)_(사전분포…)` → 두 항 모두 `[…]_(라벨)`; T07 `…_(= 0  at mode)` → `[∇ h(û_i)^T(u_i-û_i)]_(= 0  at mode)`. 그 외 HTML/DOCX 텍스트 차이 없음(각 트랙 HTML 1줄, DOCX document.xml의 해당 1문단). 다른 8개 트랙·루트·T10 산출물은 동일(타임스탬프만 다른 DOCX는 git checkout).
  - 스크린샷: scratchpad/eng/shots/{t03_ch04_hill,t05_ch02_elbo,t07_ch04_laplace}.png.
- T10 `model_audit_manifest.json` `python_version` 3.8.18→3.8.16: pinned env에서 `--re-run-pipelines --update-results` 실행 후 pipeline 01 결과·산출물은 `git checkout`으로 복원. results_pipeline_02 변경 3개 파일은 휘발 필드만(timestamp_utc, python_version, serving_benchmark 지연·처리량, summary_report의 시각·지연·처리량·`~1.1 seconds`); 가중치 해시 등 비휘발 값 불변.
- 이식성: `requirements-pipelines.txt` 헤더와 assertion 4 불일치 메시지에 "pinned env·기준 플랫폼에서만 정확; 다른 BLAS/libm/CPU에서 모델 해시·손실만 다르면 환경 차이 → --skip-pipelines 또는 --update-results" 추가. 비교는 완화하지 않음.
- 검증: 11개 빌더 `--skip-pipelines` 9/9, T10 `--re-run-pipelines`(갱신 없이) 9/9·결과 변경 없음, 스위처 11/11·링크 OK, `precommit_checks.py --all` OK.

### 후속 3(2026-10-21): 첨자/위첨자 인자가 중괄호 없는 명령·기호일 때 리터럴 `_`/`^`가 노출되던 결함
- 원인: HTML 경로의 단일 문자 첨자 규칙이 `[a-zA-Z0-9*]`만 허용 → 기호 매핑 후의 `q_ϕ`, `p_θ`, `\sum^∞`, `X^+`, `ℝ_+` 등이 그대로 노출(T02는 `ℝ_+`가 `<em>`으로 오인되는 부수 결함도 유발).
- 수정(`scripts/course_markdown.py`, HTML 경로만): 단일 토큰 규칙을 `[a-zA-Z0-9*+\-\u0080-]`(ASCII 영숫자, `*`, 매핑된 명령 기호=비ASCII 문자, 맨 `+`/`-`)로 확장. `_\text{..}`/`_\mathrm{..}`은 중괄호 인자로 기존 처리, `\_` 리터럴은 변경 없음. 유니코드/DOCX 경로는 단일 문자 첨자를 `q_ϕ`로 두는 기존 규약이 이미 중괄호 인자 처리와 일관되므로 변경 없음(DOCX 텍스트 변화 없음).
- 측정(수식 블록 내 태그 제거 후 `[_^][^\s<_a-zA-Z0-9{(]`): 242 → 0(전 트랙). 변경 라인(전부 이 규칙 안에서 설명됨): T01 4(`^∘`), T02 8(`^+`4, `_ε`, `_γ`, `^α`, `_+`2), T03 9(`^+`15, `_+`, `^∞`), T04 3(`^+`2, `^-`2), T05 84(`_ϕ`82, `_θ`41, `^θ`10, `^∞`5, `^+`4, `^→`, `^∘`, `_ψ`), T06 16(`^∞`16, `^α`2, `_ε`, `^-`, `^+`4), T07 20(`^∞`8, `^+`4, `_Λ`8, `^θ`7, `^∘`, `^ω`2, `^α`3), T08 1(`^∞`2), T10 3(`_ϕ`4, `_θ`). 루트·T09 무변경, 책(.md) 무변경.
- 부수 효과(규칙 내): T02 `ℝ<em>+<sup>N</sup>`…`</em>+` 오인 제거 → `ℝ<sub>+</sub><sup>N</sup>`; T05 본문에 이스케이프된 `&lt;span…q_ψ…` 텍스트 1곳도 같은 규칙으로 `q&lt;sub&gt;ψ…`(원래 이스케이프된 마크업이 노출되는 별개의 기존 결함은 그대로).
- 스크린샷: scratchpad/eng/shots/t05_ch02_elbo_v2.png, t02_formula.png.
- 검증: 11개 빌더 `--skip-pipelines` 9/9, 스위처 11/11, 링크 OK, `precommit_checks.py --all` OK, 타임스탬프만 다른 DOCX는 `git checkout`(실변경 DOCX는 underbrace 건 T03/T05/T07 3개뿐).

### 후속 4(2026-10-21): 한국어 본문이 `<pre><code>`(이스케이프된 수식 마크업 노출)로 렌더되던 결함
- 진단(T05 ch03 90행, T04 ch04, T08 ch01 등): 원문은 `1. **…**:` 목록 항목 아래에 4칸 들여쓴 `$$…$$`와 4칸 들여쓴 설명문이 이어지는 유효한 목록 연속 구조. 수식 전처리(`replace_math_in_markdown_html`)가 `$$…$$`를 `"\n<div class=math-display>…</div>\n"`로 치환하면서 앞선 들여쓰기를 지워 `<div>`가 0열에 놓임 → 목록이 끊기고 뒤따르는 4칸 들여쓴 문장이 mistune에서 들여쓴 코드 블록이 됨. 원인은 작성자 들여쓰기가 아니라 우리 치환(들여쓰기 소실)이며, 같은 원인으로 수식을 낀 목록이 여러 `<ul>/<ol>`로 쪼개지고 번호가 재시작됨.
- 수정(`scripts/course_markdown.py`, HTML 경로 `replace_math_in_markdown_html`): display 수식 앞 줄의 들여쓰기를 보존해 `<div>`를 같은 들여쓰기에 둠. mistune `indent_code` 규칙은 끄지 않음(들여쓰기 코드 의존 0건이라 필요 없고, 원인 수정이 더 정확). DOCX 경로는 들여쓰기가 아닌 펜스(```)로만 코드를 만들어 영향 없음(DOCX 변화 없음). 챕터 원문·책(.md) 무변경.
- 목록: `<pre><code>` 미표기 블록 365개 = (a) 원문 펜스 323 + (b) 의도적 들여쓰기 코드 0 + (c) 우발적 본문/목록 연속 42. 수정 후 (c)=0, (a) 323개 불변(HTML 파서가 아니라 mistune 렌더 대조). 트랙별 (c): root 1, T01 1, T02 2, T03 1, T04 8, T05 12, T06 6, T07 0, T08 7, T09 3, T10 1.
- 효과: 이스케이프된 `&lt;span class=&quot;math` 53 → 0(T02 1, T03 4, T04 10, T05 24, T06 2, T08 7, T09 2, T10 3). 변경 영역은 prose-now-paragraph(코드 블록 → 단락)와 list-item-restored(수식으로 쪼개졌던 목록이 하나로 복원, 항목 내부가 `<li><p>`의 loose 구조로 렌더) 두 종류뿐이며, 텍스트 내용 손실은 없음(토큰 대조에서 이스케이프 마크업 잔여만 차이). T07은 `<pre>` 변화는 없고 목록 복원만 변경.
- 검증: 11개 빌더 `--skip-pipelines` 9/9, 스위처 11/11, 링크 OK, `precommit_checks.py --all` OK, DOCX 실변경 0(타임스탬프만 달라 `git checkout`), 책 무변경. 스크린샷: scratchpad/eng/shots/{t05_ch03_fixed,t07_list_restored,t08_list_item}.png.

## R code blocks (2026-10-21)
- 범위: 챕터의 `r`/`R` 코드 블록은 요청의 6개가 아니라 **10개**(T07 9개 + T01 1개, ```R 대문자 포함). 모두 실행. 챕터 원문은 수정하지 않음(제안 diff만).
- 환경: `./.envs/r-course`(conda-forge+bioconda, osx-arm64, `.gitignore`에 `.envs/` 추가), 명세 `requirements-r.yml`(커밋 대상). R 4.4.3 / Bioconductor 3.20. 버전: SingleCellExperiment 1.28.0, DESeq2 1.46.0, limma 3.62.1, edgeR 4.4.0, scater 1.34.1, slingshot 2.14.0, tradeSeq 1.20.0, miloR 2.2.0, speckle 1.6.0, nebula 1.5.8, BiocParallel 1.40.0, Matrix 1.7.5, dplyr 1.2.1, ggplot2 4.0.3, patchwork 1.3.2, glmmTMB 1.1.15.2, car 3.1.5, MASS 7.3.66, Signac 1.16.0 (SeuratObject 5.4.0, Seurat 5.5.1), TMB 1.9.23. **IHW는 설치 불가**.
- conda에 없던/풀리지 않던 것: r-nebula(conda 없음 → CRAN 소스), slingshot/tradeSeq/miloR(bioconda 빌드가 osx-arm64에 없는 Bioc 의존 패키지 요구 → BiocManager 3.20 소스 설치; miloR은 RcppArmadillo용 `-std=gnu++17` Makevars와 env의 conda 컴파일러 PATH 필요), IHW(lpsymphony 소스 빌드가 SYMPHONY configure의 MSVC 스타일 플래그 `-TP /openmp`로 실패; bioconda lpsymphony는 R 4.5 전용이라 병용 불가). nebula 설치로 TMB가 갱신돼 glmmTMB를 소스로 재빌드. 기존 conda env는 건드리지 않음.
- 하네스(챕터 미포함): <scratchpad>/eng/rblocks (블록 block01–10.R, setup_*.R 합성 데이터, run_named.sh / run_one.sh, 로그 out*.log, 설치 로그 install_extra*.log).

| # | 파일:행 | 결과 | 비고 / 원인 |
|---|---|---|---|
| 1 | T01 03_SNATAC_QC_METRICS:67 | 실행 성공(API) | 합성 Signac 객체(fragments+annotation+blacklist). NucleosomeSignal/TSSEnrichment/FractionCountsInRegion 인자 `(object, regions, assay)` 시그니처 일치. 임의 데이터라 subset 임계값에서 세포가 0개 → "No cells found"(데이터 문제). 블록은 `peak_region_fragments`, `passed_filters`(10x 메타데이터)와 `blacklist_regions`(GRanges), `pbmc`를 이전 텍스트에서 가정 |
| 2 | T07 01_PROBABILITY_MODELS:349 | 성공 | 정상. `library(glmmTMB)`는 블록에서 사용하지 않음(선택적 정리) |
| 3 | T07 03_PSEUDOBULK:262 | 성공 | edgeR `glmQLFTest(coef="conditionMetastasis")` 정상 |
| 4 | T07 04_GLMM:227 | 성공(경고) | 블록의 자체 시뮬레이션에서 random-slope 모델이 "non-positive-definite Hessian" 경고, anova LRT의 `fit_rs` 행이 NA. 시뮬레이션에 subclone 랜덤 기울기가 없어 식별 불가한 데이터 특성(코드 오류 아님) |
| 5 | T07 04_GLMM:302 | 성공 | 함수만 정의 → 하네스에서 호출. `group_cell`이 NULL(이미 정렬)일 때의 대비가 이미 포함됨 |
| 6 | T07 05_NEBULA:228 | **실패** | `Error in as(count,"dgCMatrix"): no method ... "NULL"`. 원인: `group_cell()`은 세포가 이미 donor별로 연속 정렬돼 있으면 NULL을 반환(소스 확인: `if (is.unsorted(id)==FALSE){cat("The cells are already grouped."); return(NULL)}`)하는데 이 블록은 NULL 대비가 없음(블록 5에는 있음). 세포를 섞은 입력에서는 성공(양쪽 검증) |
| 7 | T07 06_TRAJECTORY:209 | 성공 | evaluateK/fitGAM/associationTest/earlyDETest/plotSmoothers 인자 모두 유효(합성 600유전자×300세포) |
| 8 | T07 07_MILO:212 | **실패** | `plotNhoodGraphDA` → `UMAP isn't in reducedDim(x) - choose a different layout`. 원인: 기본 layout이 reducedDim "UMAP"인데 블록은 UMAP을 만들지 않음(Harmony/PCA만 가정). 앞 단계(buildGraph/makeNhoods/countCells/testNhoods)는 시그니처와 일치해 성공. `countCells(sample=...)`는 인자명 `samples`의 부분 일치로 동작 |
| 9 | T07 08_COMPOSITIONAL:271 | 성공 | propeller `transform="asin"` 출력 열이 블록의 주석과 일치 |
| 10 | T07 09_TROUBLESHOOTING:246 | **미실행(IHW 설치 불가)** | DESeq2 부분은 실행 성공. IHW는 소스(IHW 1.34.0 NAMESPACE/`ihw.DESeqResults` 문서)로만 확인: `results(dds, filterFun=ihw, alpha=)` 와 `metadata(res)[["ihwResult"]]` 사용법이 문서와 일치(문서 주석에도 "equivalent: results(dds, filterFun = ihw)") |

### 제안 diff (미적용, sci-checker 검토 후 적용)
```diff
--- tracks/07_count_statistics_abundance/docs/05_FAST_CELL_LEVEL_NEBULA.md (block at line 228)
@@ step 3 (group_cell) @@
   offset  = meta_df$total_counts # 또는 size factor
 )
+# 이미 donor별로 연속 정렬되어 있으면 group_cell()은 NULL을 반환하므로 원본을 그대로 사용
+if (is.null(data_grouped)) {
+  data_grouped <- list(count = counts_mat, id = meta_df$patient_id,
+                       pred = design_mat, offset = meta_df$total_counts)
+}
--- tracks/07_count_statistics_abundance/docs/07_DIFFERENTIAL_ABUNDANCE_MILO.md (block at line 212)
@@ step 1 @@
+sce <- runUMAP(sce, dimred = "HARMONY")   # plotNhoodGraphDA/buildNhoodGraph의 기본 layout은 reducedDim "UMAP"
 milo_obj <- Milo(sce)
```
- 검증: 위 두 수정 블록은 `block06_fixed.R`, `block08_fixed.R`로 실행해 exit 0 (06은 정렬/섞인 입력 모두 통과). 이미 UMAP이 있는 객체라면 runUMAP 줄은 불필요하다는 주석을 같이 둘지 sci-checker/작성자 판단. (선택) 블록 8의 `sample = "sample_id"`는 `samples = "sample_id"`로 정식 인자명 권장, 블록 2의 미사용 `library(glmmTMB)` 제거, 블록 4는 시뮬레이션에 donor별 subclone 효과를 넣거나 "Hessian 경고는 시뮬레이션 데이터 특성" 주석 추가.
- 미검증 항목: 블록 6의 7단계(NEBULA-HL 재검증)는 합성 데이터에서 후보 유전자가 0개라 실행되지 않음; 블록 10의 IHW 실행 경로.

### R 블록 수정 적용 (2026-10-21, review/r_blocks_review.md 기준)
- 적용(T07 docs만): (1) `05_FAST_CELL_LEVEL_NEBULA.md` group_cell NULL 폴백 + offset 주석(양수 배율값, 로그 아님); (2) `07_DIFFERENTIAL_ABUNDANCE_MILO.md` 조건부 `set.seed(42); runUMAP(sce, dimred="HARMONY", n_dimred=30)` 및 `countCells(samples=)`; (4) `04_GLMM_HIERARCHICAL_MODELING.md` Hessian/NA 주석. `library(glmmTMB)` 유지, 05의 Python rpy2 폴백은 미적용(미검증).
- 사전 검증: 최종 결합 코드(block04/06/08_final)를 `.envs/r-course` 하네스에서 실행해 모두 exit 0 (06은 정렬/섞인 입력 모두, 08은 UMAP 없음·있음(runUMAP 생략) 모두).
- 재빌드: T07 `--skip-pipelines` 9/9, 스위처·링크 OK, `precommit_checks.py --all` OK. 출력 변경은 T07 Book(.md)·HTML·DOCX(document.xml 1개 파트)뿐이며 변경 줄은 위 블록의 추가/수정 줄(Book +20/-… , HTML 동일 범위)뿐.

## GPU blocks (Colab) (2026-10-21)
- 산출물(미커밋): `review/colab/gpu_blocks_test.ipynb`(nbformat 4 검증, 29셀)와 생성 스크립트 `review/colab/build_gpu_notebook.py`(챕터 블록은 원문에서 줄 번호로 verbatim 추출). 사용자가 Colab(T4)에서 Run all → `gpu_results.json` 다운로드 → 회신. 예상 15–30분(설치 8–15분, RAPIDS 5–12분).
- 설치(PyPI 메타데이터로 확인, 2026-10-21): `rapids-singlecell` 0.18.1(extras `rapids-cu12`/`rapids-cu13`, Python ≥3.12; cudf/cugraph/cuml/cupy/cuvs 포함) — CUDA 메이저를 nvidia-smi로 감지해 extra 선택; `scvi-tools` 1.5.1, `cell2location` 0.1.5(scvi-tools≥1.3.0 요구). 버전은 고정하지 않고 결과 파일에 기록.
- 셀 구조: GPU 확인(없으면 명확한 메시지로 중단) → 설치(torch/scvi stack, RAPIDS) → import 확인 → 하네스(`run_block`: 블록별 격리 네임스페이스, try/except, PASS/FAIL/SETUP_FAIL, 예외·traceback·경과시간·경고) → 블록 7개 → `gpu_results.json` 저장·표 출력·다운로드.

| ID | 위치 | 내용 | 분류 | 하네스 처리 |
|---|---|---|---|---|
| B1 | T05 00_FOUNDATION:182 | 환경 진단(`torch.cuda...`) | runnable alone | 없음 |
| B2 | T05 09_TROUBLESHOOTING:139 | AMP 학습 루프(`torch.cuda.amp`) + 간이 scIB | runnable alone(CPU 지표) + GPU smoke | 함수만 정의 → toy 모델로 `train_single_cell_transformer_amp(device="cuda")` 호출 |
| B3 | T02 05_SPATIAL_DECONV:222 | cell2location `train(use_gpu=True)`, `export_posterior(use_gpu)` | needs synthetic setup | 합성 reference(600×200)/Visium(150 spot, 가짜 이미지) h5ad 작성; `max_epochs` 250→20, 30000→300(`FULL_EPOCHS=True`로 원문값) |
| B4 | T08 07_RAPIDS:184 | rmm+cupy+rapids_singlecell 파이프라인 | needs synthetic setup + call | 함수만 정의 → 합성 atlas(6000×4000) h5ad로 3개 함수 호출; RMM 풀 60/76GB→1/4GB |
| B5 | T08 08_GPU_GRAPH:171 | cuvs CAGRA + cugraph Leiden | needs synthetic setup + call | 함수만 정의(원문에 "GPU 미실행" 명시) → 20000×50 난수로 3개 함수 호출 |
| B7 | T08 09_TROUBLESHOOTING:205 | VRAM 가드 데코레이터 + cuML PCA | needs synthetic setup + call | `execute_safe_pca`를 5000×100 CuPy로 호출 |
| B6 | T08 09_TROUBLESHOOTING:105 | `rmm.reinitialize(managed_memory=True, ...64GB)` | runnable after replacement | `64GB`→`2GB`(T4 16GB), 전역 할당자를 바꾸므로 마지막에 실행 |
- 제외(GPU 필요하나 블록 실행 대상 아님): ASCII 다이어그램/표(무언어 펜스) 10개 — T05 00:153, T08 00:115·175, T08 07:20·49, T08 09:18·82·122·156, T10 00:18 (illustrative-only); bash — root 03_QC:77 CellBender `--cuda`(실제 10x raw h5 필요), T08 09:67 `nsys profile`(Nsight 도구), T10 01:298 `apptainer exec --nv`(컨테이너 런타임 필요), T08 00:138 pip 설치(GPU 줄은 주석; illustrative). 문서의 R/Python 인라인 언급(`cudaMalloc` 등 서술)은 블록 아님.
- 로컬(CPU) 검증: 노트북 전 코드셀 `ast.parse` 통과(마법 `!` 줄 제외), 임베디드 챕터 블록 7개도 구문 검사 통과(빌드 시 강제). scRNAseq_anal env(py3.8, torch 1.13)에서 하네스 로직 실행: B2 챕터 블록의 CPU 부분 성공·post 단계가 `Torch not compiled with CUDA`로 예상대로 실패, B3 합성 h5ad 생성·재로딩·가짜 Visium 이미지로 `sc.pl.spatial` 성공, B4 합성 atlas 생성 성공, 나머지는 cupy/rmm/cell2location 부재로 예상 실패(하네스가 `FAIL`로 기록하고 JSON/표 출력까지 정상). **GPU에서만 검증 가능**: B3(use_gpu 수용 여부·학습), B4(rsc 전 체인·RMM 풀), B5(cuvs.cagra 인자/반환형, cugraph.leiden 반환), B7(cuML PCA `svd_solver="jacobi"`), B6, B2의 GPU 호출부.
- 부수 발견: (1) B1 챕터 블록이 Apple MPS(torch 1.13)에서 `torch.randint(..., device=mps)`의 `aten::remainder`가 미구현이라 실패(구버전 torch 한정, CUDA 무관). (2) T08 00:138 `pip install numpy>=1.24.0 ...`는 따옴표가 없어 bash에서 `>=`가 리다이렉션으로 해석됨(별개 결함, 미수정).

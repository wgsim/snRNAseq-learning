# 오케스트레이터 재검증 로그

서브에이전트 S1 판정을 오케스트레이터(Opus)가 독립적으로 재검증한 기록. 원본 파일은 수정하지 않음.

## Track 00 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-00-001 | ✅ 동의 | GEO GSE210358 brief: 제목 "Lineage plasticity in prostate cancer depends on JAK/STAT inflammatory signaling", Chan/Zaidi/Karthaus/Pe'er/Sawyers, PMID 35981096, Homo sapiens + Mus musculus. Tang et al. 2022 아님. |
| SCI-00-002 | ✅ 동의 | Single Cell Portal API SCP864: "Regenerative potential of prostate luminal cells revealed by single-cell analysis (Human)", 데이터 파일은 human. 원문 "마우스 … 인간 코호트 아님"은 반대. |
| SCI-00-012 | ✅ 동의 | `docs/01_FOUNDATIONS_AND_WET_LAB.md:68` "편향 없음 … 모든 핵을 균일하게 추출"은 단정적 오류. snRNA-seq 세포유형별 포획 편향은 문헌상 보고됨(에이전트 제시 PMID 32487174, 35042723). |
| SCI-00-023 | ✅ 동의 | 같은 챕터 `:12,18`은 v3(R1 28bp, UMI 12bp)를 설명하는데 `:59-67` STARsolo 명령은 v2 화이트리스트 `737K-august-2016.txt`를 쓰고 `--soloUMIlen 12`가 없음(기본값 10). v3 데이터에 적용하면 UMI와 바코드가 틀림. |
| SCI-00-047 | ✅ 동의 | infercnvpy 0.6.1 시그니처: `genomic_position_from_gtf(gtf_file, adata=None, *, ...)`. `(adata, gtf_file=...)`로 호출하면 gtf_file에 값이 두 번 들어가 TypeError. |
| SCI-00-063 | ✅ 동의 | scanpy 1.12.4 시그니처: `dpt_groups_pseudotime(adata, *, color_map, palette, show, marker, return_fig, save)`. `color=` 인자가 없어 TypeError. |

S2 샘플 재검증 (4건): SCI-00-032 (MAD = 0.6745σ, 수학적으로 옳음), SCI-00-034 (Scrublet README는 샘플별·raw counts 실행을 명시), SCI-00-054 (Azimuth 공식 참조 목록에 prostate 없음), SCI-00-077 (근거 없는 "98%" 수치) → 4건 모두 판정 타당. 단 SCI-00-054의 SingleR/celldex 부분은 에이전트가 직접 조회하지 않았으므로 Plausible 유지.

## Track 07 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-07-001 | ✅ 동의 (S1 유지) | `00_COUNT_STATISTICS_ROADMAP.md:57-72`는 Squair et al.에 "수백 개 데이터셋" 결과와 수치 표(Wilcoxon 80~98%, glmmTMB 4.9~5.3% 등)를 귀속함. Squair 2021 초록(PMID 34584091, Nat Commun 12:5692)은 "생물학적 반복 간 변동을 무시하면 위발견이 생긴다"는 정성적 결론만 제시함. 표의 수치와 glmmTMB 평가는 1차 문헌 근거 없이 출처를 붙인 것이므로, 학습자가 그대로 인용할 위험이 있음. |
| (03:292 factor 기준 수준) | ✅ 동의 | `03_*.md:312`의 `factor(c("Primary",…,"Metastasis"))`에는 levels 지정이 없어 R 기본값인 알파벳순이 적용되므로 기준 수준은 Metastasis가 되고, 계수 이름은 `conditionPrimary`가 됨. 따라서 `:292`의 `coef="conditionMetastasis"`는 오류이며 주석에 적힌 방향도 반대임. |
| (05 nebula 호출) | ✅ 동의 | lhe17/nebula `R/nebula.R:42-44` 시그니처는 `model='NBGMM'`(NBGMM/PMM/NBLMM), `method="LN"`(LN/HL), `cpc=0.005`이고 `cutoff_cpc` 인자는 없음. `R/scToNeb.R:25`의 시그니처는 `scToNeb(obj, assay, id, pred, offset, verbose)`로 `count=` 인자가 없음. 본문의 `model="LN"/"HL"`과 `cutoff_cpc`는 모두 오류임. |

참고: Track 07 요약에서 ID 번호가 빠진 항목(`SCI-07-0xx`)이 있어 sci-checker에 확정 번호를 쓰도록 요청함.

## Engineering Batch 1 (eng-verifier, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| ENG-004/005/006 | ✅ 동의 | liana 1.10.0에서 `hasattr` 확인: `li.mt.nichenet`, `li.mt.spatial_bivariate`, `li.rs.get_consensus_resources` 모두 False. 실제 존재하는 것은 `li.mt.bivariate`, `li.rs.select_resource`(True). |
| ENG-007 | ✅ 동의 | muon 0.2.0: `mu.atac.tl.motif_enrichment`, `mu.atac.tl.chromvar` 모두 False. |
| ENG-008 | ✅ 동의 | cellrank 2.3.3: `cr.pl.fate_probabilities` False. |
| ENG-021 | ✅ 동의 (영향 큼) | `scripts/inject_global_track_switcher.py:367-374`는 `backups/master_courses_pre_switcher/`에 백업이 있으면 무조건 `shutil.copy2(backup, full_path)`로 HTML을 되돌림. 따라서 챕터를 수정하고 빌드한 뒤 스위처를 주입하면 HTML이 예전 내용으로 돌아감. 수정 작업(향후 단계)의 전제 조건이므로 우선 수정 대상임. |
| ENG-024 | ✅ 동의 | `pipelines/04_*.py:33`이 `<conda-base>/envs/scRNAseq_anal/lib/python3.8/site-packages`를 하드코딩함. 같은 패턴이 Track 01·02 파이프라인 4개에도 있음. |

프로세스 이슈: eng-verifier가 scratch 내 로그 디렉터리를 `rm -rf`로 삭제함(규칙 위반, 프로젝트 영향 없음). 이후 배치에 재고지함.

## Track 04 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-04-023 | ✅ 동의 | kharchenkolab/numbat `inst/bin/pileup_and_phase.R`가 받는 인자는 `--label --samples --bams --barcodes --gmap --eagle --snpvcf --paneldir --outdir --ncores --bulk --smartseq`. 본문 `:228-239`의 `--bam --barcode --gts`는 존재하지 않음. |
| SCI-04-024 | ✅ 동의 | numbat `R/main.R:68-69` `run_numbat = function(count_mat, lambdas_ref, df_allele, ...)`. `lambdas_ref`는 기본값이 없는데, 본문 `:259-267` 호출에서 빠져 있음. |

참고: sci-checker가 Track 04 근거를 바탕으로 SCI-00-050을 S2→S3로 하향함(cBioPortal GISTIC 직접 계산). Track 00 집계 변경: S1 6 / S2 33 / S3 45.

## Engineering Batch 2 + Track 08 (eng-verifier, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-08-001 | ✅ 동의 | tiledbsoma 2.3.0의 `options/_tiledb_create_write_options.py:23`은 `DEFAULT_CELL_ORDER = DEFAULT_TILE_ORDER = "row-major"`이고, `:148-156`은 사용자가 지정하지 않으면 이 기본값을 씀. 본문 `tracks/08_large_scale_sparse/docs/04_*.md:61`의 "TileDB-SOMA는 힐베르트 곡선을 채택"은 오류임(TileDB 코어가 Hilbert 순서를 옵션으로 지원하는 것과 혼동한 것으로 보임). 같은 챕터의 학습 목표(`:8`)와 §2 전체가 이 전제 위에 있음. |
| SCI-08-002 | ✅ 동의 (근거 수용) | CAGRA 논문(arXiv 2308.15136)의 그래프 최적화는 rank 기반 재정렬과 역방향 엣지 추가임. "비등방성 각도 가지치기"는 논문에 없음. 논문 본문은 오케스트레이터가 직접 열지 않았고 에이전트가 제시한 근거를 수용함. |

렌더링 결과(ENG-038~047)는 Playwright 실측 기반으로 판정이 타당함. ENG-041(DOCX `<w:br/>`이 `<w:r>` 밖에 있음)이 Word에서 어떻게 보이는지는 확인하지 못함.
프로세스: eng-verifier는 Playwright가 프로젝트 루트에 만든 `.playwright-mcp/`를 trash로 정리함. 초기 탐색 시점(세션 시작) 프로젝트 목록에는 이 폴더가 없었으므로 원본 손실은 없음.

## Track 10 (eng-verifier, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-10-002 | ✅ 동의 | ChEMBL API: `CHEMBL1201585` = TRASTUZUMAB, `CHEMBL1082407` = ENZALUTAMIDE. RxNav: `1306289` = cobicistat/elvitegravir/emtricitabine/tenofovir, `1307298` = enzalutamide. 본문 `06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md:428`의 enzalutamide 매핑은 다른 약을 가리킴. DrugBank ID는 미검증. |
| SCI-10-003 | ✅ 동의 | OLS4: `CL:0002327` = mammary gland epithelial cell, `CL:0000084` = T cell. 본문 `:431`의 "prostate luminal cell" 등 매핑이 틀림. |

## Track 09 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-09-025 | ✅ 동의 | `03_*.md:231-252`: `generate-permit-list`에 `-k`와 `-u`를 동시에 지정함. `quant`에서는 `-r cr-like`와 `--resolution cr-like`를 중복 지정하고, 존재하지 않는 `--use-recs`를 씀. `-k`/`-u` 상호 배타성과 pyroe 산출물 이름(fl84 vs fl89)은 에이전트의 소스 대조를 수용함. 직접 실행은 하지 않음. |
| SCI-09-050 | ✅ 동의 | Ensembl REST(hg38) AR 영역 엑손: `67,685,941–67,686,126`(exon 3, 끝 67,686,126), `67,694,673–67,696,075`(CE3/AR-V7 말단 엑손), `67,711,402–67,711,689`(exon 4). 본문 `07_*.md:90-93,182-186`의 좌표 67,713,200 / 67,715,600 / 67,720,500와 fetch 시작점 67,710,000은 모두 실제 위치와 다름. 따라서 이 코드는 정션을 하나도 세지 못함. |
| SCI-09-068 | ✅ 동의 | `09_*.md:97-99`: v2 화이트리스트 `737K-august-2016`에 `--soloUMIlen 12`(v3 설정)를 조합함. Track 00 SCI-00-023과 같은 v2/v3 혼동이 재발함. |

## Track 05 (eng-verifier, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-05-001 | ✅ 동의 | `docs/images/03_scanvi_and_scarches_atlas_mapping.svg`에 "98.4% 정확도 분류", `01_scvi_vae_probabilistic_architecture.svg`에 "iLISI 0.94"가 적혀 있음. 같은 트랙 산출물 `pipelines/results_pipeline_01/*.json`의 값은 `unlabeled_accuracy: 0.5574`, `normalized_ilisi: 0.556`임. 교재 그림 수치가 자체 실습 결과와 충돌함. |

작업 재배분: eng-verifier가 담당 트랙(08·10·05)을 완료함. 남은 sci 트랙 중 03·06을 eng-verifier로 옮기고, sci-checker는 01·02를 담당함.

## Track 01 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-01-001 | ✅ 동의 | GEO GSE175634 제목은 "Single-Cell Sequencing Reveals Lineage-Specific Dynamic Genetic Regulation of Gene Expression During Human Cardiomyocyte Differentiation"(PMID 35061661)으로 전립선 연구가 아님. 본문 `00_*.md:124`는 이를 "Wang 2022 Nat Genet 전립선 Multiome"으로 소개함. `:125`의 GSE210358↔Tang 오귀속은 SCI-00-001과 같은 오류가 재발한 것임. |
| SCI-01-039 | ✅ 동의 | aertslab/scenicplus 저장소 트리(main, v1.0.0 둘 다)에 `eregulons.py` 모듈이 없음. 따라서 `08_*.md:127`의 `from scenicplus.eregulons import build_eregulons`는 ImportError. 나머지 시그니처 지적은 에이전트의 소스 대조를 수용함. |

패턴: 데이터셋 카탈로그 오류(Track 00 → 01)와 WNN 수식 오류(SCI-00-068 → Track 01 05장)가 트랙 간에 복제됨. 최종 리포트의 "트랙 간 전파 오류" 섹션에 반영 예정.

## Track 03 (eng-verifier, 2026-10-18)

S1 없음 (S2 16 / S3 27). S2는 2건을 표본으로 재검증함.

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-03-009 | ✅ 동의 | liana 1.10.0 `method/_pipe_utils/_aggregate.py:10`의 기본값은 `aggregate_method='rra'`이고, 패키지 전체를 검색해도 "borda" 문자열이 없음. 문서의 "LIANA Borda 합의"는 실제 구현과 다름. |
| SCI-03-010 | ✅ 동의 | PubMed esummary PMID 35680885: Dimitrov et al., *Nat Commun* 2022. 본문의 "Nat Methods"는 오기임. |

프로세스: eng-verifier가 scratch 임시 파일 `c09_0.py`를 `rm -f`로 삭제함(두 번째 규칙 위반, 프로젝트 영향 없음).

## Track 02 (sci-checker, 2026-10-18)

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-02-001 | ✅ 동의 | GEO GSE98544: "mRNA-seq comparison of N2 and tps1;tps2 double mutants…", organism *C. elegans*. BioStudies E-MTAB-11270: "ZXDA/B dependent transcription in HCT116 human colon carcinoma cells". 둘 다 전립선 공간전사체가 아님. |
| SCI-02-025 | ✅ 동의 (근거 수용) | Baysor CLI·TOML 규격 대조는 에이전트가 수집한 README/예시 설정(scratch `t02/baysor_*`) 기반. 직접 실행은 하지 않음. |
| SCI-02-030 | ✅ 동의 | cell2location 모듈 `__init__`(scratch 사본 `c2l__cell2location_module.py:75-89`)은 `N_cells_per_location=8.0`을 받고 `**kwargs`가 없음. 따라서 본문 `05_*.md:275`의 `N_cells_per_spot=15`는 TypeError. |
| SCI-02-051 | ✅ 동의 | COMMOT `_spatial_communication.py:251-` 시그니처의 `cost_type` 기본값은 `'euc'`이고 분기는 `'euc'`/`'euc_square'`뿐임. 본문 `08_*.md:209` `cost_type="quadratic"`은 지원하지 않는 값임. 나머지 4개 지적은 에이전트의 소스 대조를 수용함. |

## Track 06 (eng-verifier, 2026-10-18)

S1 없음 (S2 29 / S3 25). S2 표본 재검증:

| ID | 판정 | 재검증 근거 |
|---|---|---|
| SCI-06-007 | ✅ 동의 | scanpy 1.12.4 `neighbors/__init__.py:943-952` `_get_dpt_row`는 Σ (λ_j/(1−λ_j))²·(ψ_j(i)−ψ_j)²로 계산함. 본문 `tracks/06…/docs/01_*.md:315`의 가중치 λ²/(1−λ²)는 이와 다름. |

비고: SCI-06-048(pipeline02 "G-PCCA"가 정답 라벨을 복사)과 SCI-06-049(운명 확률 840/1500이 0.5 폴백)는 S2로 두되, 최종 리포트 §4-C에 "합성 파이프라인 결과의 신뢰성" 사례로 강조함.

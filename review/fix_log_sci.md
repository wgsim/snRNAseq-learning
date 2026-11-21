# fix_log_sci.md — sci-checker S1 원본 수정 로그 (17건)

- 범위: 원본 소스 챕터(.md)만 수정. 생성물(Master_Book/Course.html/docx, 슬라이드, backups, 빌드 스크립트)은 미수정.
- 새로 넣은 사실의 검증: GEO text view, Single Cell Portal API, EGA metadata API, UCSC REST(hg38), 각 도구 공식 소스/문서(STARsolo, pyroe, alevin-fry, Numbat, COMMOT, cell2location, nebula, SCENIC+ v1.0.0, infercnvpy, scanpy), 로컬 실행 검증(scanpy `scatter`, R `model.matrix`/`parse`).
- 실행 환경: eng-verifier의 `scratchpad/eng/.venv_api`. cell2location(pyro 비호환)과 commot(`np.Inf`, numpy 2)은 import가 실패하여 설치된 소스를 AST/텍스트로 대조함.
- 판단이 갈린 점: 규칙 7(Track 00 오류의 타 트랙 복제는 미수정)과 담당 목록이 겹치는 SCI-01-001, SCI-09-068은 담당 S1 목록에 명시되어 있어 수정함. 필요하면 §5 일괄 수정 단계에서 되돌려도 영향 없음.

## 수정 내역

| ID | 파일:줄(수정 후) | 변경 요약 (전→후) | 근거 | 검증 방법 |
|---|---|---|---|---|
| SCI-00-001 | docs/00_CURRICULUM_AND_ROADMAP.md:115 | `Tang et al. / MSKCohort`, `GSE210358`, snATAC+RNA/Multiome, AR/Amphicrine/DNPC/t-NEPC 4대 아형, ETV4 레귤론 → `Tang et al. (Science 2022)`, `GSE199190`(SuperSeries), **bulk** ATAC+RNA-seq(단일세포 아님), 인간 CRPC 오가노이드·PDX·세포주, 4군(AR-dependent/NE/Wnt/SCL) | GSE199190 = 'Chromatin profiles classify castration-resistant prostate cancers suggesting therapeutic targets' (PMID 35617398); GSE210358 = Chan 2022 'Lineage plasticity … JAK/STAT' (마우스/오가노이드 bulk RNA-seq, PMID 35981096) | GEO text view (Series_title/pubmed_id/type), Tang 4군 분류는 원 검토(SCI-00-001) |
| SCI-00-002 | docs/00_CURRICULUM_AND_ROADMAP.md:116 | SCP864를 '마우스'로 표기, 'scRNA-seq / snRNA-seq', 'Luminal P 줄기세포능 획득' → SCP864=인간 / SCP859=마우스(GEO GSE146811), scRNA-seq(10x), '거세 후 분화된 관강세포가 강화된 재생능 획득(희귀 Lum-P 포함)' | SCP864 name '… (Human)', SCP859 '… (Mouse)'; GSE146811 title 'Differentiated prostate luminal cells acquire enhanced regenerative potential after androgen ablation' (PMID 32355025) | Single Cell Portal API, GEO text view |
| SCI-00-012 | docs/01_FOUNDATIONS_AND_WET_LAB.md:68; 파생 :72, :163 | '편향 없음 … 모든 핵을 균일하게 추출' → '효소 해리 편향은 줄지만 핵 회수 편향은 존재(유형별 회수율, 긴 유전자·lncRNA 편향)'; 파생: '세포 손실을 근본적으로 차단' → '크게 줄입니다', '균일하게 포획' → '효소 유래 손실을 줄이되 유형별 회수율 차이 주의' | 원 검토 SCI-00-012 근거(Koenitzer 2020, Denisenko PMID 32487174, Gupta PMID 35042723 등; 본 단계에서 재조회하지 않음) | 문장 대조 |
| SCI-00-023 | docs/02_UPSTREAM_PROCESSING.md:64-65, 70; 파생 :17 | STARsolo `737K-august-2016.txt` + UMI 기본(10) → `3M-february-2018.txt`(gunzip) + `--soloUMIlen 12`; 주석 박스로 v2/v3 구분 추가; :17 '737,280개 화이트리스트'를 'v2 737,280개 / v3·v3.1은 3M-february-2018'로 | STARsolo.md: 'V2 737K-august-2016.txt', 'V3 3M-february-2018.txt.gz (gunzip it for STARsolo)', 'For V3, specify --soloUMIlen 12' | STARsolo 문서 원문 대조(scratchpad/starsolo.md). v3 화이트리스트 항목 수는 파일 미확보로 기재하지 않음 |
| SCI-00-047 | docs/05_CELL_TYPE_ANNOTATION.md:70 | `genomic_position_from_gtf(adata, gtf_file=…)` → `genomic_position_from_gtf("gencode.v46.annotation.gtf", adata=adata, gtf_gene_id="gene_name")` | infercnvpy 0.6.1 `(gtf_file, adata=None, *, gtf_gene_id='gene_name', adata_gene_id=None, inplace=True)` | `inspect.signature` (venv_api) |
| SCI-00-063 | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:109 | `sc.pl.dpt_groups_pseudotime(adata, color=[…])` → 유전자별 `sc.pl.scatter(adata, x='dpt_pseudotime', y=gene)` 루프 | scanpy 1.12.4 `dpt_groups_pseudotime(adata, *, color_map, palette, show, marker, return_fig, save)`에 `color` 없음 | `inspect.signature`; 합성 AnnData로 diffmap→dpt→`scatter(x='dpt_pseudotime', y='AR')` 실행 성공 |
| SCI-01-001 | tracks/01_multiome_snatac/docs/00_MULTIOME_ROADMAP.md:126 (표 L120-124) | 표의 Wang 2022(GSE175634)·Tang/MSK(GSE210358) 행 **삭제**, 아래에 '출처 주의' 인용문 추가. **전립선 Multiome accession 미확정**(원 논문 Data Availability 확인 필요) | GSE175634 = 심근 분화 scRNA-seq(PMID 35061661); GSE210358 = 마우스·인간 오가노이드 bulk RNA-seq(Series_type 'Expression profiling'); 'Wang 2022 Nat Genet 전립선 Multiome' 논문 자체를 확인하지 못함 | GEO text view. 삭제 근거: 두 accession 모두 Multiome이 아니고 대체 accession을 검증하지 못함 |
| SCI-01-039 | tracks/01_multiome_snatac/docs/08_EGRN_AND_SCENICPLUS.md:101-153 | SCENIC+ 코드를 v1.0.0 API로 재작성: `get_search_space(species, assembly, upstream=[1000,150000], downstream=[…])` → `calculate_regions_to_genes_relationships(search_space_key=…)`; `calculate_TFs_to_genes_relationships(scplus_obj, tf_file=…)`; `merge_cistromes` → `build_grn(min_target_genes=5)` → `format_egrns`; `get_eRegulons_as_signatures` → `make_rankings` → `score_eRegulons(ranking=…)`; 비실재 `build_eregulons`·`score_eregulons`·`ranking_key`·`search_space_upstream` 제거. 주석으로 'v1.0.0 기준, git main은 시그니처 상이' 명시 | aertslab/scenicplus 태그 v1.0.0 소스(scenicplus_class.py, enhancer_to_gene.py, TF_to_gene.py, grn_builder/gsea_approach.py, cistromes.py, utils.py, eregulon_enrichment.py) | 소스 시그니처 대조 + Python `ast.parse` 통과(실행은 하지 않음). 결과 열 이름은 소스 확인 불가하여 `uns['eRegulon_AUC']['Region_based']`만 기재 |
| SCI-02-001 | tracks/02_spatial_transcriptomics/docs/00_SPATIAL_ROADMAP.md:198, 201 | Berglund `GSE98544` → **accession 미확정**(논문 Data Availability: EGA 표기 `EGAS0000100300`은 EGA 조회 불가, 카운트 행렬 spatialtranscriptomicsresearch.org); Erickson `E-MTAB-11270`, Visium Standard, 림프절 n=14 → `EGAS00001006124`(EGA), 플랫폼 '원 논문 확인 필요', '전립선 포함 여러 장기의 양성·악성 공존 조직(120,000개 이상 영역)' | GSE98544 = C. elegans(PMID 29063832); E-MTAB-11270 = HCT116 ZXDA/B siRNA; EGAS00001006124 title 'Spatial atlas of clonal copy number alterations in co-existing benign and malignant tissue' (>120,000 regions, multiple organs); **정정:** 이전 검토에서 제시한 Berglund = EGAS00001003000은 틀림(EGA에서 소세포 전립선암 별개 연구) → findings_track_02.md·claims_ledger(_sci).csv의 해당 문구도 정정함 | GEO text view, BioStudies API, EGA metadata API(metadata.ega-archive.org), Berglund PMC 본문 |
| SCI-02-025 | tracks/02_spatial_transcriptomics/docs/04_CELL_SEGMENTATION_ALGORITHMS.md:185, 191-212 | `--prior-segmentation` 플래그 → `baysor run -c … -o … transcripts.csv dapi_cellpose_mask.tif`(위치 인자); `confidence="qv"` 삭제(qv 필터는 사전 전처리 주석); `min_molecules_per_cell`을 [data]로; Xenium 열 이름(x_location/y_location/z_location/feature_name); `exclude_genes="NegControl*,BLANK_*,antisense_*"`; Parquet → CSV(`transcripts.csv`) | Baysor docs run.md/segmentation.md(`PRIOR_SEGMENTATION` second positional; csv 예시), example_config.toml, Xenium용 baysor_xenium.toml | 공식 문서·설정 파일 대조(실행 안 함). Parquet 지원 여부는 확인 불가라 CSV로 기재 |
| SCI-02-030 | tracks/02_spatial_transcriptomics/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:275-276, 295, 307-310; 파생 :141, :194 | `N_cells_per_spot=15` → `N_cells_per_location=15`(기본 8 주석); `detection_alpha` 주석 보정; `q05_cell_abundance_w_sf_<type>` → `q05cell_abundance_w_sf_<type>`(obsm 키는 `q05_cell_abundance_w_sf`); 본문 기호 `N_{cells_per_spot}`·`q05_cell_abundance` 명칭 정정 | cell2location 0.1.5 소스: `Cell2location.__init__(adata, cell_state_df, model_class, detection_mean_per_sample, detection_mean_correction, **model_kwargs)`, module `N_cells_per_location=8.0`, `export_posterior`가 `obsm[f"{k}_cell_abundance_w_sf"]`, `sample2df_obs` 열 이름 `f"{summary_name}{name_prefix}_{site_name}_{i}"` | 설치 소스 텍스트 대조(import 실패: pyro 비호환). `use_gpu=True` 인자는 ENG 영역이라 미수정 |
| SCI-02-051 | tracks/02_spatial_transcriptomics/docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:187-245 | `signaling_type="Paracrine"` → `"Secreted Signaling"`; `df["ligand"]` → `df.iloc[:, 0]`; `cost_type="quadratic"` → `"euc_square"`; `lr_pair=("TGFB1","TGFBR2")` → `("TGFB1","TGFBR1_TGFBR2")`(플롯 제목 포함); `background="cell_type"`·`title=` → `background="cluster", clustering="cell_type"` + `ax.set_title`; `k` 주석 보정 | commot 0.0.3 소스·CellChatDB.ligrec.human.csv(열 '0'–'3'; TGFB1 수용체 4종 모두 이종복합체; TNFSF11→TNFRSF11A 'Secreted Signaling'), `plot_cell_communication` 시그니처 | 설치된 commot 소스 AST로 시그니처, 패키지 내 CSV로 DB 항목 확인(import는 numpy 2 때문에 실패) |
| SCI-04-023 | tracks/04_cancer_genomics_clones/docs/04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md:208-243 (+ :107, :119, :306) | `--bam/--barcode/--gts` → `--label --samples --bams --barcodes --gmap --snpvcf --paneldir --outdir --ncores`; 함수 인자 `snp_vcf, panel_dir, gmap, sample_id`; SHAPEIT4 → Eagle2(본문 L107·L119·문제 2) | numbat `inst/bin/pileup_and_phase.R` add_option 목록('Missing one or more always required arguments: --bams, --snpvcf, --outdir, --paneldir, --gmap'; cellsnp-lite + eagle), Numbat vignette(1000G SNP VCF `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz`, panel `1000G_hg38`) | 스크립트 소스·vignette 대조, Python `ast.parse` 통과 |
| SCI-04-024 | 동일 파일 :275-283 | `run_numbat(count_mat, df_allele, …)` → `lambdas_ref = ref_hca` 추가(필수 인자); allele 파일명 `{sample_id}_allele_counts.tsv.gz` | numbat R/main.R `run_numbat(count_mat, lambdas_ref, df_allele, gtf=NULL, genome='hg38', …, t=1e-5, …)`; 출력 `{outdir}/{sample}_allele_counts.tsv.gz`; vignette `run_numbat(count_mat, ref_hca, df_allele, …)` | 소스·vignette 대조 |
| SCI-07-001 | tracks/07_count_statistics_abundance/docs/00_COUNT_STATISTICS_ROADMAP.md:57-74 | '수백 개 데이터셋 + Type I error 수치표(80~98%, 4.9~5.3% …)' → 논문의 실제 설계(18 gold-standard, 14 방법, AUCC)와 **정성적 결론표**(수치 제거) | Squair 2021 PMID 34584091(PMC8479118 초록·본문: 'contingent on their ability to account for variation between biological replicates … can discover hundreds of differentially expressed genes in the absence of biological differences'; 18 gold-standard, 14 methods, GLMM이 500세포 이상에서 pseudobulk에 수렴) | 논문 본문 텍스트 검색 |
| SCI-07-026 | tracks/07_count_statistics_abundance/docs/03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:312-314 | `factor(c("Primary",…,"Metastasis"))` → `levels = c("Primary","Metastasis")` 지정 | R 기본 수준은 알파벳순 → 계수 `conditionPrimary` | R 실행: levels 지정 시 `conditionMetastasis`, 미지정 시 `conditionPrimary` 확인; 블록 `parse()` 통과(edgeR 미설치로 전체 실행은 안 함) |
| SCI-07-035 | tracks/07_count_statistics_abundance/docs/05_FAST_CELL_LEVEL_NEBULA.md:263-264, 271-278, 281, 288-289, 314, 360-362; tracks/07…/04_GLMM_HIERARCHICAL_MODELING.md:314-337 | `model="LN"`/`"HL"`/`"NBLG"` → `model="NBGMM", method="LN"/"HL"`; `cutoff_cpc` → `cpc`; 결과 `omega`,`kappa` 열 → `overdispersion` 요소의 `Subject`,`Cell` 열을 `cbind`; Python rpy2 인자 동일 수정; 04장 `scToNeb(count=…)` → `group_cell()`(+ 이미 정렬이면 NULL 반환 처리) | nebula man/nebula.Rd(model 'NBGMM','PMM','NBLMM'; method 'LN','HL'; cpc), `R/nebula.R`(summary: logFC_/se_/p_/gene, overdispersion 열 'Subject','Cell'), scToNeb.Rd(obj, assay, id, pred, offset), group_cell.Rd(NULL 반환 조건) | 문서·소스 대조, 3개 R 블록 `parse()` 통과(nebula 미설치로 실행 안 함) |
| SCI-09-025 | tracks/09_upstream_alignment_bam/docs/03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:183, 195, 200, 211, 245-249 | `pyroe make-splici … 89 …` + `splici_fl89.*` → read length `90`, 산출물 `splici_fl85.fa`/`splici_fl85_t2g_3col.tsv`(주석 '90−5=85'와 일치); `generate-permit-list`에서 `-k` 제거(`-u`만); `quant`에서 `--use-recs` 및 중복 `--resolution cr-like` 제거(`-r cr-like`만) | pyroe `make_splici_txome.py`(flank = read_length − flank_trim, prefix `splici_fl{flank}`, `_t2g_3col.tsv`)와 `bin/pyroe`(positional genome, gtf, read-length, output-dir; `--flank-trim-length`); alevin-fry main.rs(`filter-method` ArgGroup에 `knee-distance`·`unfiltered-pl` 상호 배타, quant `-r` 필수·`--use-mtx`만 존재, `--use-recs` 없음) | 소스 대조(명령 실행 안 함) |
| SCI-09-050 | tracks/09_upstream_alignment_bam/docs/07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:47, 90-94, 182-186, 329-346 | Exon3 donor 67,713,200→67,686,126; Exon4 acceptor 67,720,500→67,711,402; CE3 acceptor 67,715,600→67,694,673; 질의 구간 67,710,000-67,725,000→67,686,000-67,712,000; 인트론 길이 7.3/2.4 kb→25.3 kb(25,275 bp)/8.5 kb(8,546 bp); Q2 예시 POS 67,713,150·`50M2400N50M`→POS 67,686,077·`50M8546N50M`(donor 67,686,126, acceptor 67,694,673), 본문 '2.4 kb 하류'→'약 8.5 kb 하류' | UCSC hg38 ncbiRefSeq: NM_000044.6 exon3 67,686,009-67,686,126, exon4 67,711,401-67,711,689; NM_001348061.1(AR-V7, variant 3) 말단 exon 67,694,672-67,696,075 (UCSC는 0-based 반개구간 → 1-based 시작 +1) | UCSC REST API 직접 조회(2026-10-18), Q2 산술 재계산(67,686,077+50−1=67,686,126; +1+8,546=67,694,673), Python 블록 `ast.parse` |
| SCI-09-068 | tracks/09_upstream_alignment_bam/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:91 (주석 L89-90) | `--soloCBwhitelist /ref/10x_737K-august-2016.txt`(UMI 12 bp와 불일치) → `/ref/3M-february-2018.txt` + v2/v3 화이트리스트 주석 | STARsolo.md(SCI-00-023과 동일 근거) | 문서 원문 대조 |

## 동일 오류 잔존 위치 (미수정, §5 일괄 수정 대상)

원본 소스(.md/.py):
1. `docs/01_FOUNDATIONS_AND_WET_LAB.md:36` — 'Tang et al., 2022'를 단일세포/단일핵 연구·4-state 모델의 근거로 인용 (SCI-00-015).
2. `tracks/02_spatial_transcriptomics/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:198` — `N_{\text{cells\_per\_spot}}` (SCI-02-030 파생; 실제 인자명은 `N_cells_per_location`).
3. `tracks/07_count_statistics_abundance/docs/02_PSEUDOREPLICATION_PITFALL.md:9, :150` — 'Wilcoxon 제1종 오류 80~95%'. 같은 장의 자체 시뮬레이션 서술이지만 SCI-07-001과 같은 근거 없는 고정 수치 (조건 의존으로 완화 필요).
4. `tracks/10_production_pipelines_mlops/pipelines/01_dsl2_dataflow_and_dag_engine.py:323` — `ValueChannel("whitelist", "737K-august-2016.txt")` 문맥 점검 필요 (Track 10, ENG 소관).
5. `tracks/01_multiome_snatac/docs/00_MULTIOME_ROADMAP.md:97` — `pip install git+https://github.com/aertslab/scenicplus.git`(main)은 08장 수정본이 전제한 v1.0.0 API와 다름 (필요 시 태그 고정; 설치 검증은 하지 않음).

생성물(재빌드로 반영, 미수정): 00 `snRNAseq_Prostate_Master_Book.md/.html/.docx`, 01·02·04·07·09 각 `*_Master_Book.md`·`*_Master_Course.html`, `tracks/09_upstream_alignment_bam/slides/02_BAM_Mechanics_and_Isoform_Mining_Deck.html`(구 AR-V7 좌표), `backups/master_courses_pre_switcher/*`.

## 미확정으로 남긴 항목

- SCI-02-001: Berglund 2018 데이터 accession — 원 논문이 쓴 `EGAS0000100300`이 EGA 조회로 확인되지 않아 '미확정'으로 표기(카운트 행렬은 spatialtranscriptomicsresearch.org).
- SCI-01-001: 공개 전립선 10x Multiome accession — 검증된 대체물이 없어 오류 두 행을 삭제하고 '미확정' 주석을 남김. 같은 절의 10x 공식 데이터 URL(미존재, SCI-01-002 S2)은 이번 범위 밖이라 미수정.
- SCI-02-001: Erickson 2022의 플랫폼 세부(Visium/ST)·시료 수는 확인 불가라 '원 논문 확인 필요'로 표기.

## 기타 메모

- SCI-02 findings/ledger 정정: 이전 보고의 'Berglund = EGAS00001003000'은 오류였음(해당 EGA 연구는 소세포 전립선암). `review/findings_track_02.md`, `claims_ledger_sci.csv`, `claims_ledger.csv`의 해당 문구를 '미확정'으로 고쳤다.
- cell2location `train(use_gpu=True)` 등 scvi-tools 버전 의존 인자와 rpy2 변환 세부는 ENG 소관이라 손대지 않음.

## 2단계(전파 오류)

- 범위: 승인된 파일만 수정(원본 소스 .md, 수작업 슬라이드 1개). 생성물·빌드 스크립트·금지 범위 파일은 미수정.
- 새로 확인한 사실: GEO text view(GSE143791, GSE141445), Europe PMC 초록(PMID 34719426, 33420488, 31361600, 29017058), Kfoury 논문 본문(PMC8578470: '9 metastatic PCa patients … seven Benign BM controls', 샘플당 평균 2,490 cells, Radium-223 환자 2명 상세 분석 제외), GitHub 태그 목록(aertslab/scenicplus v1.0.0), v3 화이트리스트 크기(6,794,880: 앞서 Track 09 검토의 PMID 35354805/38671451 및 산술 교차 확인).

| 항목 | 파일:줄(수정 후) | 변경 요약 (전→후) | 근거 | 검증 방법 |
|---|---|---|---|---|
| A. 4대 아형 출처 (SCI-00-015) | docs/01_FOUNDATIONS_AND_WET_LAB.md:36, 174 | '단일세포/단일핵 연구(Tang et al., 2022)가 정립한 Four-State Lineage Model' → Labrecque 2019(bulk RNA-seq/IHC, ARPC·ARLPC·AMPC·DNPC·SCNPC)·Bluemn 2017(DNPC) 기반의 '교재 편의상 분류'로 표기, Tang 2022는 AR-dependent/NE/Wnt/SCL 별개 분류임을 명시; L174 '(Four-State Lineage Model)' 명칭 삭제 | Labrecque PMID 31361600, Bluemn PMID 29017058 초록, Tang PMID 35617398 | Europe PMC 초록 확인 |
| A | tracks/02_spatial_transcriptomics/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:202 | 'Track 00에서 정립된 4대 분자 계통' → 'Track 00에서 정리한 … (Labrecque 2019·Bluemn 2017 기반 표현형 분류)' | 위와 동일 | 문장 대조 |
| A | tracks/07_count_statistics_abundance/docs/06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md:153 | t-NEPC/DNPC 두 경로 분기를 사실처럼 서술 → 분류 출처(Labrecque 2019, Bluemn 2017) 명시 + '분기 구조는 모식적 가정' | 위와 동일 | 문장 대조 |
| B. v3 화이트리스트 크기 | tracks/09…/docs/00_ALIGNMENT_AND_BAM_ROADMAP.md:138, 287 | 3,793,305(3.79×10^6) → 6,794,880(6.79×10^6) | SCTK-QC PMID 35354805 표(PBMC1k V2 737,280 / V3 6,794,880), hadge PMID 38671451(Cell Ranger v3.1.0 raw 6,794,880), 10x v3 목록 파일 크기(~110 MB ≈ 6.79M×17 B)와 정합 | 산술 재계산 |
| B | …/01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:9, 92, 94, 177-181, 385-405 | M=3,793,305→6,794,880; 점유율 0.088%→0.158%; N=10,000 기대 충돌쌍 13.18(13.2쌍)→7.36(7.4쌍), v2 대비 5.14배→9.2배 감소; Q2(N=20,000) λ_v3 52.72→29.43, 기대쌍 52.7→29.4(약 105→59개 세포), '5.14배·80% 급감'→'약 9.2배·약 89% 급감'(271→29쌍), P≈1.0000 유지 | 산술: 1e8/(2×6,794,880)=7.358, 4e8/(2×6,794,880)=29.434, 67.817/7.358=9.216, 1−29.43/271.27=0.8915, 6,794,880/4^16=0.15821% | python 재계산 |
| B | …/02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:159 | 3,793,305개 → 6,794,880개 | 위와 동일 | 문장 대조 |
| C. WNN 거리식 (SCI-00-068, SCI-01-024/025) | docs/07_MULTIOMICS_AND_FUTURE.md:63-67 | d_WNN=√(w d² …) → 가중 친화도 θ_WNN(c,c')=Σ_M w_{c,M} θ_M(c,c'), 가중치는 within/cross 예측 커널 비의 softmax(세포 c 기준, 비대칭); '명확히 분리' → '분리하는 데 도움이 될 수 있습니다' | Hao 2021 Cell PMID 34062119 STAR Methods; Seurat clustering.R(원 검토 SCI-01-024/025에서 대조) | 원 검토 근거 재사용 |
| C | tracks/01_multiome_snatac/docs/05_WNN_MULTIMODAL_INTEGRATION.md:72-90, 115, 200 | §2.3: θ_M·affinity softmax → κ_within=exp(−d_within/σ), κ_cross, score=κ_within/(κ_cross+10⁻⁴), w=softmax(score), θ_M 하이퍼파라미터 삭제; §2.4: 거리 가중 제곱합 → 가중 친화도 합 θ_WNN로 이웃 선택; §3.1 항목 3(가중 제곱합에 의한 분리) 및 Q2 답안의 softmax 설명을 새 식에 맞춤 | 위와 동일 | 문장 대조. Q3의 'w_ATAC 0.65가 입증' 서술(SCI-01-026)은 범위 밖이라 미수정 |
| D. 3′ 라이브러리 융합 포획 (SCI-09-060, SCI-04-036) | tracks/09…/docs/08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:101-102 | '수십~수백 개의 융합 분할 리드가 반드시 포획' → 접합부 직접 포획은 원리적으로 거의 불가(ERG exon 4~poly(A) 약 4.7 kb), ERG 발현을 대리 지표로, 5′/롱리드/타깃 증폭 권장; 'ERG 엑손 4~17' 표현 제거 | UCSC hg38 NM_004449.4 엑손 길이(원 검토 SCI-09-060) | 원 검토 근거 재사용 |
| D | tracks/04_cancer_genomics_clones/docs/06_SOMATIC_MUTATIONS_AND_FUSIONS.md:57, 347 | §2.1에 '10x 3′ snRNA-seq에는 discordant pair 부적용, 접합부 분할 리드 원리적으로 드묾' 적용 범위 주의 추가; Q3 Cluster 1의 '분할 리드 다수 검출' → 'ERG 15배 상향(접합부 리드는 소수 가정)' | 01장 Q3(접합부 미포획)와 내부 모순 해소 | 문장 대조 |
| E. 가상 결과 단정 (SCI-00-066) | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:149, 203-205 | 'DNPC 마스터 레귤론(ETV4/ELK1/FOSL1)' → '후보 레귤론(가설)', Bluemn 2017(bulk)의 FGF/MAPK 활성 보고만 근거로 표기, '확정 문헌 없음·pySCENIC 검증 필요' | Bluemn PMID 29017058 초록 | 초록 확인 |
| E (SCI-01-041) | tracks/01_multiome_snatac/docs/08_EGRN_AND_SCENICPLUS.md:3, 67, 161, 179, 185 | '4대 마스터 eRegulon이 명확히 도출'·'정립한다'·'증명하는 회로를 확립했습니다'·'생물물리학적으로 입증된 직접 인과망' → '가정해 볼 수 있는 후보(가설적 예시)'·'추론하는 프레임워크'·'후보 연결망(인과성은 실험 검증 필요)' | SCENIC+ PMID 37443338 초록('candidate upstream TFs … candidate target genes') | 초록 대조(원 검토) |
| E (SCI-07-051) | tracks/07_count_statistics_abundance/docs/07_DIFFERENTIAL_ABUNDANCE_MILO.md:194-198 | §4.1 제목에 '(가상 예시)', n=10/14·log2FC +2.8/+3.4·Spatial FDR 수치가 '실제 공개 코호트 결과가 아닌 설명용 가상 값'임을 명시, '폭발적 변화가 포획된다'→'포착될 수 있다' | 공개 근거 없음(원 검토) | 문장 대조 |
| F. Kfoury·Chen 코호트 | docs/00_CURRICULUM_AND_ROADMAP.md:112, 114 | Kfoury: 'scRNA-seq/snRNA-seq, n=26' → scRNA-seq(10x 3′ v2), 골전이 환자 9명(종양·병변 인접·원위 골수 3분획)+양성 골수 대조 7명(GEO에는 비전립선 전이·마우스 샘플도 등록), 핵심 발견을 T 세포 소진·골전이 대식세포·CCL20-CCR6로; Chen: 'CRPC 포함 n=13' → 근치적 전립선 절제술 원발 종양 13개(36,424 세포), 별도 11개 샘플 추가(CRPC 포함), 핵심 발견을 KLK3 이소성 발현·aEC로 | GEO GSE143791 샘플 제목(BMET×9, BMM×7 + 기타), PMC8578470 본문; GEO GSE141445 요약(13 tumors, 36,424 cells, 'additional 11 samples'), PMID 33420488 초록 | GEO text view, 논문 본문/초록 |
| F | tracks/03_tme_cell_communication/docs/00_TME_COMMUNICATION_ROADMAP.md:69-70 | Kfoury '~35,000 cells'·'mCRPC'·'BM-MSC·SPP1-CD44' → 9+7 샘플 구성, 총 세포 수 미확정(논문 평균 ~2,490/샘플), CCL20-CCR6 등; Chen '~42,000 cells·원발+CRPC·CAF IL6/STAT3' → 원발 13개 36,424 cells, KLK3 이소성 발현·aEC | 위와 동일 | 위와 동일. 같은 표의 Bill et al. 행(71)은 범위 밖이라 미수정 |
| F | tracks/07…/docs/00_COUNT_STATISTICS_ROADMAP.md:9, 87-103, 106-117, 296-297 | '15인 mCRPC' → 9명+대조 7명; Chen '18 individuals, >100,000 cells' → 13개 원발 종양·36,424 cells; 연령 52~78·세포 수율 200/5,000 등 출처 없는 수치는 삭제 또는 '가정적 예시'; 문제 2의 15명 → '가상 시나리오, 환자 9명', 나머지 14명 평균 50~60개 → 8명 평균 100개 | 위와 동일 | 산술(3,200+8×100=4,000) |
| F | tracks/07…/docs/03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:169-175, 380 | Kfoury 교란 변수 'Primary(N=5) vs Bone Met(N=10), 연령 54~79' → 분획(Tumor/Involved/Distal, 환자 9명)·양성 대조 7명; Chen 'Paired Pre/Post ADT, multi-center' → 단일 병원 13개 원발 종양, 환자 내 전/후 생검 아님; 문제 3에 '가상 설계' 명시 | 위와 동일 | 문장 대조 |
| G. 종말 상태 '비가역' (SCI-04-041) | tracks/04_cancer_genomics_clones/docs/07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:3, 74, 84-87, 130-132, 245-246, 308 | '영구 봉인/비가역적 탈바꿈/q_4j=0이므로 BAT에도 역분화 불가/영구 침묵' → '후성유전학적 억제/전환될 수 있음(EZH2/SOX2 의존, 일부 가역: Ku 2017, Mu 2017)', q_4j=0은 '모형 가정이며 현실 모형은 q_4j>0 허용', BAT 역분화는 별도 근거 필요 | Ku PMID 28059767('Ezh2 inhibitors restore androgen receptor expression and sensitivity'), Mu PMID 28059768('reversed by restoring TP53 and RB1 function or by inhibiting SOX2') | 초록 확인(원 검토) |
| H. 잔존 소스 | tracks/02_spatial_transcriptomics/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:198 | N_{cells_per_spot} → N_{cells_per_location} | cell2location 소스(1단계 SCI-02-030) | 문자열 대조 |
| H | tracks/07…/docs/02_PSEUDOREPLICATION_PITFALL.md:9, 150 | 'Type I 오류 80~95%' → '크게 초과, 정도는 N·ρ·n에 의존'; '80~95%의 거짓 발견' → '다수(조건 의존)' | Squair 2021 PMID 34584091(1단계 SCI-07-001) | 문장 대조 |
| H | tracks/01…/docs/00_MULTIOME_ROADMAP.md:97-99 | `pip install git+…/scenicplus.git` → `…@v1.0.0` + 주석('태그 존재는 확인, 설치·의존성 해결은 직접 실행해 검증하지 않음') | GitHub 태그 목록 refs/tags/v1.0.0, v1.0a2 | GitHub API |
| I. 슬라이드 | tracks/09…/slides/02_BAM_Mechanics_and_Isoform_Mining_Deck.html:437-439, 505, 509, 616-617 | AR-V7 좌표 67,713,200/67,715,600 → 67,686,126/67,694,673; 인트론 길이 2.4 kb→8.5 kb(V7), 7.3 kb→25.3 kb(FL) (SVG 라벨 포함) | UCSC hg38 NM_000044.6, NM_001348061.1 (1단계 SCI-09-050) | grep 확인(HTML 구조 미변경). 이 덱에는 Q2 예시 없음 |

## 2단계 미확정·보류 항목

- Kfoury 2021의 총 세포 수: 논문이 총계를 제시하지 않아 '미확정'으로 표기(평균 약 2,490 cells/샘플만 기재).
- Chen 2021의 별도 11개 샘플이 GSE141445에 포함되는지: GEO 요약만으로 불확실(제출자 'raw data not available') — 문장에서 '논문은 별도 11개 샘플 추가'로만 서술.
- SCENIC+ `@v1.0.0` 설치 가능성: 태그 존재만 확인.
- C의 within/cross 커널 σ 정의는 Seurat 구현(세포별 이웃 거리)에 근거한 일반 서술이며 정확한 k는 명시하지 않음.

## 2단계 범위 밖에서 새로 발견한 동일 오류 위치 (미수정)

1. tracks/09…/slides/01_Raw_Reads_Alignment_and_Deduplication_Deck.html:127, 301 — v3 화이트리스트 '3.79M, 0.088%', 'M=3.79M → 기대 충돌 ≈13쌍'.
2. tracks/09…/slides/02_BAM_Mechanics_and_Isoform_Mining_Deck.html:664, 795, 1111 등 — 'TMPRSS2-ERG 50M50S 분할 리드' 검출 서술(3′ 라이브러리 한계 미언급, 한 장의 개요 문구).
3. tracks/07…/docs/02_PSEUDOREPLICATION_PITFALL.md:333 와 04_GLMM_HIERARCHICAL_MODELING.md:349-372 — 'Kfoury 코호트(환자 15명)'.
4. tracks/07…/docs/06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md:392-394, 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:310 — 환자 15명 설정(가상이나 '가상 설계' 미표기).
5. tracks/00 docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:128 — AR+/Amphicrine/… 틀을 합격 기준 맥락에서 사용(SCI-00-066 순환), ETV4/ELK1/FOSL1 표기는 없음을 확인.
6. tracks/04…/docs/06_SOMATIC_MUTATIONS_AND_FUSIONS.md 의 SAM 예시(POS 41.49 Mb에 GN:ERG, SA +1bp 불일치) 및 SPOP–야생형 ERG 서술(L357) — SCI-04-036 본문 일부이나 SCI-09-057/058/061 유형이므로 별도 수정 필요.
7. tracks/01…/docs/05_WNN_MULTIMODAL_INTEGRATION.md Q3 해설의 'w_ATAC=0.65가 입증' (SCI-01-026).
8. tracks/01…/docs/08_EGRN_AND_SCENICPLUS.md 의 AR 650 kb vs ±150 kb 검색 공간 (SCI-01-042), FGFR/MEK 병용 제안 (SCI-01-043).
9. Master Book/Course HTML/DOCX 등 생성물 전반(재빌드 필요).

### 2단계 추가분 (후속 요청)

| 항목 | 파일:줄 | 변경 요약 (전→후) | 근거 | 검증 방법 |
|---|---|---|---|---|
| 슬라이드 01 (v3 화이트리스트) | tracks/09…/slides/01_Raw_Reads_Alignment_and_Deduplication_Deck.html:127, 301 | '3.79M 화이트리스트(0.088%)' → '6.79M(6,794,880개) … (0.158%)'; 'N=10k, M=3.79M → P≈99.99%, 기대 충돌≈13쌍' → 'M=6.79M → P≈99.94%, 기대 충돌≈7.4쌍' | 챕터 재계산과 동일(6,794,880/4^16=0.158%, 1e8/(2M)=7.36, 1−e^−7.36=0.99936) | python 산술; HTML 구조 미변경 |
| 슬라이드 02 '50M50S 분할 리드' | …/02_BAM_Mechanics…Deck.html:593, 664, 795 | **미수정**: 세 곳 모두 CIGAR/분할 리드·FLAG 0x800 일반 개념 소개이며 3′ snRNA-seq에서의 포획을 단정하지 않아 오류가 확실하지 않음 | 문장 대조 | — |
| Track 07 Kfoury 15명 | tracks/07…/docs/02_PSEUDOREPLICATION_PITFALL.md:333; 04_GLMM_HIERARCHICAL_MODELING.md:349 | 'Kfoury 골전이 코호트(환자 15명)' → '실제 GSE143791은 골전이 9명(·대조 7명)이며, 이와 유사한 구조를 가정한 가상 설계(환자 15명)' | GEO GSE143791 / PMC8578470 (1단계 F항과 동일) | 문장 대조 |
| Track 07 06·08장 | 06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md:392-394, 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:310 | **미수정**: 두 문제 모두 Kfoury를 인용하지 않는 일반 가상 시나리오(환자 15명·30명) | 본문 대조 | — |

## S2 단계 / Track 00

- 대상: Track 00 sci S2 33건(상태는 `review/s2_status_sci.csv`). 결과: fixed 15 / softened 13 / already-fixed 5 / deferred 0 / rejected 0.
- 원칙: Plausible(근거 없음) 항목은 새 단정 없이 수치 삭제·"예시" 표기·완화만 수행. 새로 쓴 사실은 이미 검증된 것(MAD 상수 0.6745/1.4826 수학, Azimuth 공식 참조 목록, 6,794,880, Ch04의 layers["raw_counts"])에 한정.

| ID | 파일:줄 | 전→후 | 근거 | 검증 |
|---|---|---|---|---|
| SCI-00-006 | docs/00_CURRICULUM_AND_ROADMAP.md:16 | "해리 편향을 완전히 배제…표준 분석법", "핵막만 분리" → "크게 줄임…핵 분리, 고유 편향 존재" | 원 검토 SCI-00-006 | 문장 대조 |
| SCI-00-007 | 동:14 | "단일 세포 분리가 사실상 불가능" → "어렵고 편향 위험 큼"+Kfoury 반례 | Kfoury 2021(1~2단계에서 확인) | 문장 대조 |
| SCI-00-008 | 동:98; docs/01:141; docs/03:15-16,24-25,148,177-178; docs/08:112 | mt% "1~2% 미만" 등 → "대개 낮음/분포 기반(예시)" | 원 검토(고정 임계 근거 없음) | grep로 장별 수치 정리 |
| SCI-00-013 | docs/01:67 | "전무함(전사 억제제)" → "크게 감소" | 3.1 레시피에 전사억제제 없음 | 문장 대조 |
| SCI-00-014 | docs/01:68,72,163 | "용혈", "면역·섬유아세포만 생존", "차단" → "용해·사멸", "취약 세포 손실 가능", "감소" | 원 검토; 카탈로그의 상피 회수 scRNA 데이터 | 문장 대조 |
| SCI-00-016 | docs/01:47,183 | "20~25%" 삭제 → "코호트·정의별 상이" | 근거 미확인(Plausible) | — |
| SCI-00-017 | docs/01:121,202-203 | RIN≤2.0 단정·가수분해 기전 → "크게 저하 가능(예시)", 탈퓨린은 주로 DNA, FFPE는 DV200 | 원 검토 | 문장 대조 |
| SCI-00-024 | docs/02:109 (17은 S1에서) | "전체 70만 개 바코드" → v2 약 73.7만 / v3 약 680만 | 6,794,880(2단계 B) | 산술 |
| SCI-00-025 | docs/02:21 | 60~80%/<40% 권장 삭제 → depth-saturation 곡선 | 10x 가이드에 목표치 없음 | — |
| SCI-00-028 | docs/02:178-179 | "30~50% 하락", "0으로 소실·붕괴" → "크게 감소(on/off 비교)", "불안정" | 원 검토 | — |
| SCI-00-030 | docs/03:14-19,28-30 | UMI 수치를 예시값으로, "T/B/NK 전부 소실" → "선택적 탈락 가능" | 원 검토 | — |
| SCI-00-031 | docs/03:24-25,162,177-178 | "검출되지 않아야", 2~3% 단정 → "소량 검출/분포 기반(예시)" | 원 검토 | — |
| SCI-00-032 | docs/03:45 | k·MAD=3σ 주장 → raw MAD=0.6745σ, scaled MAD 필요, mt%는 상한+floor/cap | 정규분포 MAD 상수(0.6745, 1.4826) | 산술(k=3→2.02σ, k=5→3.37σ) |
| SCI-00-034 | docs/03:105-135 | 병합 adata에 1회 → 샘플별 루프(raw 카운트), 이중핵 비율 0.06은 예시 | 원 검토(Scrublet 샘플별 실행) | python ast 파싱 통과 |
| SCI-00-038 | docs/04:17,165-166 | "평균 낮을수록 분산 커짐/고발현 과소평가" → 이분산성(Var 식 명시) | Poisson/NB 분산 식(본문 자체) | 수식 대조 |
| SCI-00-039 | docs/04:24,169-170 | "완벽 제거·탁월" → "크게 줄임·도움" | 원 검토 | — |
| SCI-00-044 | docs/04:130-134 | 필수 규칙 → 탐색적 가이드, 저발현 TF·double dipping 주의 | 원 검토 | — |
| SCI-00-051 | docs/05:43 | "필수 조건" → "허용하는 흔한 동반 사건" | 원 검토 | — |
| SCI-00-052 | docs/05:95-96,152 | KLK→19q 팔 증폭 오인 확립 → 일반 위양성(HLA 6p21)+KLK는 가설·미검증 | 원 검토 | — |
| SCI-00-053 | docs/05:104 | 접합부 직접 탐색 → 3′ 한계 단서(ERG 대리, 5′/롱리드) | 2단계 D와 동일 | — |
| SCI-00-054 | docs/05:116 | Azimuth/SingleR 전립선 참조 → 공식 참조에 전립선 없음, 직접 구성 | 원 검토(Azimuth 15종 참조 목록) | SingleR/celldex 부재는 미확인으로 단정 안 함 |
| SCI-00-056 | docs/05:121 | "반드시…확정 분리" → "구분, CNA 음성이면 보류" | 원 검토(§3.3 한계와 일치) | — |
| SCI-00-058 | docs/06:31 | "완전 소거" → "줄이도록 유도" | scVI 구조 | — |
| SCI-00-059 | docs/06:37-45 | `categorical_covariate_keys=["tissue_source"]` 삭제 | scVI 문서(nuisance 요인) | 코드 대조 |
| SCI-00-060 | 동 | `layer="raw_counts"` 추가, HVG·batch_key 주의 주석 | Ch04 L39 layers["raw_counts"] | grep |
| SCI-00-065 | docs/06:122-123,220-221 | MAESTER 적용성 제한·enrichment 필요 명시 | 원 검토(MAESTER 방법) | — |
| SCI-00-076 | docs/08:136-146 | "유일한 표준·완벽 통제" → "권장(혼합모형도 대안)·소표본 한계" | 원 검토 | — |
| SCI-00-077 | docs/08:168-178 | "98%" 삭제, "진정한 드라이버" → "후보", 설계·검정력·독립 검증 | 원 검토 | — |
| already-fixed | SCI-00-003/004 (2단계 F), 015 (2단계 A), 066 (2단계 E), 068 (2단계 C) | 상태만 갱신 | 각 단계 로그 | — |

다른 트랙의 같은 오류(미수정, 참고): 해당 Track 00 항목의 복제본은 Track 01~09 문서에 일부 존재(예: DNPC 비율·mt% 임계·scVI 설정은 Track 02·07 일부 장에 유사 서술) — 각 트랙 S2 단계에서 처리.

## S2 단계 / Track 07

- 대상 35건(상태·사유는 `review/s2_status_sci.csv`의 note). 결과: fixed 23 / softened 9 / already-fixed 3 / deferred 0 / rejected 0.
- 검증: 재실행(scipy: t/정규 꼬리·검정력, 02장 시뮬레이션 seed 42, 09장 simulate_power 귀무·검정력), 소스(miloR graphSpatialFDR, scCODA, nebula 문서/man, statsmodels), R/Python 코드 블록 `parse`/`ast` 통과.
- 주의: 06장 `sm.NegativeBinomial` 수정은 합성 데이터로 실행해 확인. 04·05장 R(`nebula`, `glmmTMB`)은 패키지 미설치로 문법 검사만.

| ID | 파일:줄 | 전→후 | 근거 | 검증 |
|---|---|---|---|---|
| SCI-07-002 | (note 참조) | already-fixed: 2단계 F(Kfoury 9명+양성 대조 7명, 00:9/87-115, 03:168-176 정정; 02:332 및 04:349는 후속 요청에서 "가상 설계"로 표기) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-003 | (note 참조) | already-fixed: 2단계 F(Chen 원발 13개·36,424 cells, 00/03 정정; 03:380은 가상 설계 표기) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-008 | (note 참조) | softened: 01:7,60 "예외 없이" → "일반적으로, 고발현·이질 집단에서 두드러짐(저발현은 Var≈μ)" | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-009 | (note 참조) | fixed: 01:9,224-258 Svensson 단독 저자·세부 미검증 설계 삭제, Choi 2020 방향 정정(영과잉 존재 가능, NB GLM 기준 권고), "물리적 드롭아웃은 존재하지 않음" 완화 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-014 | (note 참조) | softened: 02:88 ρ=0.02~0.15를 설명용 가정으로 표기, Kfoury/Chen 귀속 삭제(해당 논문에 ICC 보고 없음) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-015 | (note 참조) | fixed: 02:140 N_eff=120은 평균 추정 정밀도 등가이며 추론 자유도는 N−2임을 명시 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-016 | (note 참조) | fixed: 02:245-259 seed 42 재실행 값(78.00/77.20/5.70%)으로 교체, 시드 의존 범위(74~78%, 4~6%) 명시(코드 재실행) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-017 | (note 참조) | fixed: 02:343-346 t(df=13)에서 p=1e-300은 t≈4e23(t=50은 p≈3e-16), Wilcoxon |Z|≈37, 분모가 커서가 아니라 순위 편향+독립 가정이 원인으로 정정(scipy 계산) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-018 | (note 참조) | softened: 02:275-278 세포주도 처리는 웰 단위 → 독립 복제 웰 집계/혼합모형, 한 웰 내 세포 검정은 가짜 반복 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-022 | (note 참조) | fixed: 03:76 "완벽 필터링·벌크와 동일" → 샷 노이즈 감소, 공여자 간 변동 별도(유도는 세포 독립 가정) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-023 | (note 참조) | softened: 03:107-118,144 edgeR WL-EB(개념적 가중평균은 단순화)와 QL squeezeVar 분리, σ0²=1 수축 삭제 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-024 | (note 참조) | fixed: 03:130 분산 이상치(dispersion outlier) 규칙과 Cook's distance(샘플 이상치) 분리 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-025 | (note 참조) | softened: 03:147 "완벽 통제" → "비교적 잘 통제하는 편(보장 아님), DESeq2 Wald도 양호" | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-031 | (note 참조) | softened: 04:348-372 Q1을 환자 내 처치 설계로 명시, SE→0은 무작위 효과 없는 세포 수준 GLM의 성질임을 분리, 무작위 기울기 식별 조건 명시 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-032 | (note 참조) | fixed: 04:245-253 random slope를 공여자 내 변수 subclone으로 변경((1+subclone|patient_id)), LRT 경계 문제 주석(condition은 공여자 수준이라 식별 불가) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-033 | (note 참조) | softened: 04:398 라플라스 오차는 영 비율이 아니라 군집당 정보량 의존, glmmTMB는 AGHQ 미지원, nebula는 AGHQ가 아님(L10도 정정) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-034 | (note 참조) | softened: 04:424 Mundlak 항·cis(CNA 구간) 유전자 제외 서술 추가, "환자 변이 완벽 흡수" 삭제 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-036 | (note 참조) | fixed: 05:41,55-74,120-135 속도 이득의 기전을 해석적 적분(대표본 근사)으로 정정, O(N³) 밀집 역행렬은 교육용 단순화로 표기, 100~500배 삭제(논문: orders of magnitude); 04:6.2도 정정 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-037 | (note 참조) | fixed: 05:145-152,156 LN/HL 표를 nebula 문서 기준으로 수정(LN=대표본 근사 우도, HL=h-likelihood, CPS<30이면 HL 자동), 근거 없는 N>50,000/J>15·5~10배·밀리초 삭제 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-038 | (note 참조) | fixed: 05:41,218 "Perfect Donor FDR Control/완벽 분리 제어" → 함께 보정, 공여자 수준 과산포 추정 민감도(κ/PMM 권고) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-041 | (note 참조) | fixed: 06:5절 도입부에 풀링 모형(tradeSeq/Python NB-GAM) 한계 주의 + mgcv GAMM(bam, s(donor,bs="re"), s(t,donor,bs="fs")) 안내, p값은 탐색적 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-042 | (note 참조) | fixed: 06:308-319 sm.GLM(NB alpha=1 고정) → sm.NegativeBinomial(alpha MLE), Wald 행렬에서 alpha 제외, scalar=True(합성 데이터로 실행 확인) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-043 | (note 참조) | fixed: 06:145,147 S^{-1} → S^{-}(랭크 결손·improper prior, mgcv가 영공간 분리), "완벽히 전환" 완화 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-044 | (note 참조) | fixed: 06:415-416 말단 구간 단일 환자 → 비식별/외삽 명시, "완벽 방어" → "크게 완화, 환자 다양성·민감도 분석 필요" | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-045 | (note 참조) | softened: 06:352-358 과적합은 SE가 맞으면 검정력 저하, 위양성은 SE 과소추정 시 발생으로 한정 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-048 | (note 참조) | fixed: 07:126-171,334-390 Spatial FDR을 miloR graphSpatialFDR 소스 기준으로 교체(w=1/connectivity, weighting 4종, adjp=Σw·p/cumsum(w)+cummin), 커널 밀도 ρ·p/w·평균 1·"가장 독창적" 삭제, BH PRDS 서술 정정, Q1·Q2 답안 재서술 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-049 | (note 참조) | fixed: 07:76,96 Fibonacci sampling → random + refinement | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-050 | (note 참조) | fixed: 07:243 이웃 크기 기준 <5/>300 → 평균 이웃 크기 > 5×N_samples(miloR 지침) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-051 | (note 참조) | already-fixed: 2단계 E(07:194-206 가상 예시 표기, 수치는 설명용 가상 값) | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-055 | (note 참조) | fixed: 08:193,337,353-361 scCODA 자동 기준 = 영 비율 낮은 세포형 중 상대 풍부도 분산/평균 최저(소스), 대수비 분산 합 삭제, 민감도 분석 권고 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-056 | (note 참조) | fixed: 08:173-181,370,382-388 scCODA 사전분포 = 연속 완화(β=τ·β̃, τ=sigmoid(50z)), horseshoe 미사용, 유의성은 FDR 기반 임계(Newton 2004), 일반 spike-and-slab 설명은 교육용으로 구분, "완벽 억제" 완화 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-061 | (note 참조) | fixed: 09:16-50,282-317 N=군당 공여자로 정의 통일(군당 3명×5,000 vs 군당 30명×500 = 30,000 cells), SE 0.5775/0.1829, 검정력 12%/94% → 정규근사 7~8%/36~38%(β=0.3, scipy 계산), "언제나 우월" 완화 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-062 | (note 참조) | fixed: 09:185-241 Poisson GLM(귀무 기각률 ≈94%) → log-pseudobulk t-검정(귀무 4.3%, 코드 실행), docstring의 "GLMM" 삭제, 실행마다 달라지는 예시값 명시 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-063 | (note 참조) | fixed: 09:351-385 기록된 배치는 공변량 포함이 1순위, RUVr은 치료 상관 배치 성분 포착 불가, ComBat 과다 보정은 사용 방식 의존, "완벽 차단" 삭제 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |
| SCI-07-068 | (note 참조) | fixed: 04:18, 05:44 "슈도벌크/혼합모형 모두 유효, GLMM은 세포 수준 공변량이 핵심일 때 필요(Squair: 500세포↑ 수렴)" 통일 문구; Track 00 08장은 S2 Track 00에서 이미 "권장(혼합모형도 대안)"으로 정정 | 원 검토 findings_track_07.md 및 위 검증 | 위 검증 항목 |

다른 트랙의 같은 오류(미수정, 참고): Track 07 외 'Kfoury 15명'(Track 03 등), nebula API는 Track 03 등 다른 트랙에 유사 서술 가능 — 해당 트랙 단계에서 처리.

## S2 단계 / Track 04

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-04-001 | 00장 L86-87 (TCGA-PRAD, SU2C/PCF 행) | TCGA 499례→원발 333례(프로젝트 전체 ~500 병기); SU2C "333례, AR+650kb 인핸서 증폭 60~70%, TP53 53%, PTEN 40%, RB1 21%"→Cell 2015 150례/PNAS 2019 444례, AR 증폭(GISTIC 기준 ~52%), 나머지는 "코호트·정의에 따라 상이", 650kb 인핸서 tandem duplication은 Viswanathan 2018 별도 | PMID 26544944, 26000489, Abida 2019, Viswanathan 2018 | OK
SCI-04-002 | 00장 L85 (Gundem 행) | "292개 검체 심층 WGS, 상호 파종을 최초 증명"→"10명 51개 종양 WGS(평균 ~55X; 원발은 5명만 포함), 전이소 간 전파가 흔함(단클론 de novo 또는 5례 다클론 전달), AR 신호 변이의 수렴적 반복" | Gundem 2015 (PMID 25830880) | OK
SCI-04-003 | 00장 L75 | "ADT 저항성은 거의 전적으로 CNA가 주도"→"CNA가 핵심적 역할(단, TCGA 분자 아형의 약 74% 등 비-CNA 요인 병기)" | TCGA 2015 (PMID 26544944) | OK
SCI-04-006 | 00장 L268-269 | "반드시 보유해야 한다"→"보유하는 경우가 많다(원발 TCGA에서 MYC gain/amp ~31%, PTEN 손실 ~32%)"; 이배체 프로파일→"반응성 정상 상피로 단정 말고 판정 보류(CNA 약/없는 종양·CN-LOH는 inferCNV 미검출), 변이·BAF 등 직교 근거로 확인" | TCGA 2015 | OK
SCI-04-007 | 00장 L286-287 | "SU2C/PCF 코호트의 증거"→"다부위 연구(예: Gundem 2015)의 증거"로 교체, Robinson 2015는 환자당 단일 생검 설계라 환자 내 분지 진화를 직접 보일 수 없음 병기 | Robinson 2015 (PMID 26000489) 설계, Gundem 2015 | OK
SCI-04-031 | 05장 L361, L366 (Q2 해설) | "장가지 유인"→"비균일 진화 속도(UPGMA류 위상 오류)"로 용어 정정(LBA는 파시모니/ML 계통수 아티팩트); NJ "유일한 해석적 최적해, 항상 올바른 자매 계통 보장"→"균형 최소진화 기준 탐욕적 휴리스틱, 가산 거리에서 일관성, 노이즈 CNA 거리에서는 근사" | NJ 이론(Saitou-Nei, Gascuel) | OK
SCI-04-032 | 05장 L55, L225, L346 | 절단점 거리 "매우 안정적, 우수한 견고성"→"덜 영향받을 수 있음, 발현 기반 절단점은 윈도우 해상도(CopyKAT ~5Mb, inferCNV 100유전자)로만 결정되어 불확실, 1윈도우 어긋남 시 대칭차집합 2개 사건"; 임계값 0.2 "임의의 예시값" 표기 | 알고리즘 정의로부터 논리 검토 | OK
SCI-04-033 | 05장 L150, L171-172 | truncal "환자 내 모든 세포 100% 공유·21q 결실"→"해당 환자 모든 세포 공유; 아래는 ERG+ 가정 예시, 융합은 서구 원발의 ~45~50%, PTEN 손실은 TCGA 원발 ~32%뿐이므로 모든 환자에서 truncal 아님"; 융합 기전에 전좌 병기; 사건 순서는 설명용 가설 | TCGA 2015, Perner 2006 | OK
SCI-04-009 | 01 L9, L131; 05 해당부 | TMPRSS2-ERG 기전에 ~60% 결실+전좌 병기 | Perner 2006 | OK
SCI-04-013/014/015 | 02 L9, L26-38, L53-66, L81-82, L95, L107-116, §4.1, §5.2, Q2/Q3 | HMM 6상태 i6(상대 0/0.5/1/1.5/2/3), 방출=spike-in 시뮬레이션·Baum-Welch는 일반이론 주석, Viterbi.dthmm.adj; 로그 차감→클리핑(±3)→이동평균; infercnvpy HMM 없음·chrX/Y 기본 제외 | inferCNV 공식 위키/소스, infercnvpy 0.6.1 시그니처 직접 확인 | OK
SCI-04-017/018/020 | 03 L7-9, L38-39, §2, §3 | Poisson-gamma=분절(MCpoissongamma+KS), 정규화=로그 Freeman-Tukey+dlm; 분산 GMM은 교육용, 실제=계층 6군집+F-test(p<0.05, 5개 이상)·GMM 대체 정의; 2% 수치 삭제 | Gao 2021 (PMID 33462507) | OK
SCI-04-025/026 | 04 §4 서두 주석, Q3 | Numbat 실제 6상태·Poisson-lognormal·beta-binomial(gamma=20)·min_LLR=5; "유일한 방법"→"가장 직접적 단서" | Numbat 공식 코드/문서 | OK (교육용 (c,m) 표는 유지; 전면 재작성은 deferred 대신 softened)
SCI-04-029 | 05 MEDICC2 절, truncal, Q1/Q2 | MEDICC2 입력(정수 CN 0-8, 하플로타입 특이)·MED 정의 정정, 단정 완화 | Kaufmann 2022 (PMID 36376909) | OK
SCI-04-036/041 | 06, 07 | already-fixed (2단계 D, G) | 2단계 로그 | 잔존 문구 없음 확인
SCI-04-042 | 07 §2 서두, 2.1, Q1 | Ku/Mu 2017 주장을 PTEN 결손 배경·SOX2로 한정, REST 분해는 모형적 가설 표기 | PMID 28059767/28059768 초록 | OK
SCI-04-047/048/049 | 08 §2.2, §3, Q3 | 간 전이=골 유래 일반화 삭제 → Gundem 2015(10명/51종양/벌크 WGS, 전이간 전파) 인용; DTC 휴면 수년~수십 년 삭제; self-seeding 가설 표기; Q3 "명백히 입증"→parsimony 가설, 협력·선택압 기전은 가설적 예시 | Gundem 2015 (PMID 25830880) | OK
SCI-04-053 | 09 §1.1, Q 해설 2곳 | 19q 전체 4-6N 삭제, 윈도우 폭 한정 국소 피크; 산술 예시(+3 클립×3유전자/100 = 0.09 < ln1.5≈0.41) | infercnvpy lfc_clip=3 | OK
SCI-04-054 | 09 §3 | 평균 1-5x 삭제, MAPQ>=50은 예시 임계값 | 근거 없는 수치 제거 | OK
SCI-04-060 | 09 §2.1, 02 문제1 해설 | 1:1 혼합 고정 비율 삭제, GMM 센터링=휴리스틱(WGD 주의), 참조군 편향 방향성 주석 | 논리 검토 | OK

다른 트랙의 같은 오류(미수정, 참고): Track 00/05 등의 KLK 19q 증폭 서술, inferCNV HMM 상태 표기(≥5x)는 해당 트랙 단계에서 처리.

## S2 단계 / Track 09

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-09-001 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:138, 287 | [already-fixed] 2단계 B(v3 화이트리스트 6,794,880) 적용 확인: 00장 L138, 01장, 02장 | PMID 35354805, 38671451 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-002 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:310-321 | [fixed] 00장 Q3: SharedMemory 삭제, LoadAndExit→LoadAndKeep→Remove, limitBAMsortRAM 필요 병기 (parametersDefault 직접 대조) | STAR parametersDefault (직접 조회) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-003 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:267-268 (Q1 해설); 02:62-67 | [fixed] 00장 Q1 해설·02장 도식: AR hg38 exon3 끝 67,686,126, CE3 시작 67,694,673(인트론 8,546), exon4 67,711,402 (UCSC ncbiRefSeq 직접 조회) | UCSC ncbiRefSeq (직접 조회) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-004 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:94 | [fixed] 00장 L94: 엔잘루타마이드=LBD 길항, 아비라테론=CYP17A1 억제로 분리 | de Bono NEJM 2011 PMID 21612468 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-005 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:122 | [fixed] 00장: SRP173516(Burkholderia) 삭제, 10x 공개 데이터셋 페이지로 대체, 규모는 예시 표기 | ENA: SRP173516=PRJNA510167 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-006 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:203-209 | [fixed] 00장 환경 스크립트: alevin_fry import 삭제, pysam __samtools_version__/__htslib_version__, yml에 pyroe 추가, htslib 핀 주석 | pysam/version.py, v0.22.0 NEWS | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-011 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:94, 171-181, 385-405 | [already-fixed] 2단계 B에서 01장 파생 수치 재계산됨(0.158%, 7.36쌍, 9.2배, 29.4쌍, ~89%) 재검산 일치; N 정의(GEM/세포 수) 주석만 추가 | 산술 재계산 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-012 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:130-135, 377 | [fixed] 01장 §3.3·Q1: cbMinP=0.975, P∝count×10^(-Q/10), 유일 이웃은 품질 무관 보정; Q1 판정을 폐기→보정(0.99>=0.975)으로 정정 | STAR ParametersSolo.cpp, SoloReadFeature_inputRecords.cpp | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-013 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:128-130 vs 304-309 | [fixed] 01장 코드: Q>=10 조건 제거(유일 이웃은 무조건 보정), 모호 처리는 단순화 주석 | STAR SoloReadBarcode_getCBandUMI.cpp | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-014 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:213-231, 429-430 | [fixed] 01장 §5.2: Good-Turing 주변 프로파일, Dirichlet-multinomial, p>=1/(R+1)(기본 R=10,000이면 ~1e-4), knee 보존; p<=1e-6 삭제 | Lun 2019 PMID 30902100 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-019 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:90-95, 299 | [fixed] 02장 페널티 표: 0/-4/-8/-8 (scoreGap/GCAG/ATAC/Noncan), Q3 P_splice=15 삭제, scoreGapNoncanUnannotated→scoreGapNoncan, alignSJstitchMismatchNmax는 기본값임 (parametersDefault 대조) | STAR parametersDefault (직접 조회) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-020 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:10, 180-186 | [softened] 02장: 벤치마크 표를 "예시(가상)"로 표기, 디스크 0 GB 삭제, 속도 배율 단정 삭제, genomeLoad는 인덱스 적재 절감 병기 | STARsolo preprint doi:10.1101/2021.05.05.442755 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-026 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:147-152 | [fixed] 03장: --use-quant 삭제(3열 t2g로 USA 자동), A=경계 flank 내 엑손 리드, 엑손 내부 리드는 S | alevin-fry quant.rs; PMID 35277707 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-027 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:10, 160-168, 322 | [fixed] 03장: splici dense~10GB/sparse~6.5GB(He 2022), 4GB 노드 시나리오 폐기, 3.8GB·kallisto 4.5GB 등 무출처 수치 삭제 또는 예시 | He 2022 PMID 35277707 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-028 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:7, 22-27, 41, 45-50 | [fixed] 03장: selective-alignment=ksw2 점수화(CIGAR 미산출), DP 없는 것은 --sketch; kallisto|bustools는 모듈형 선구 사례로 정정, 90% 낭비 삭제 | Srivastava 2020 PMID 32894187; Melsted 2021 PMID 33795888 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-029 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:80-112 | [fixed] 03장 RAD 도식: magic 삭제, is_paired/ref_count/names/num_chunks/tag desc, 레코드 num_aln→barcode→UMI→aln, 3바이트 UMI 삭제, raw 바코드 저장 (libradicl header.rs 대조) | libradicl header.rs/record.rs (docs.rs 직접 조회) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-032 | 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md:78, 82-85, 9 | [fixed] 04장: 5종(Adjacency 추가), Percentile=평균 1% 미만 제거, "완벽 동시 해결" 삭제 | Smith 2017 PMID 28100584 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-033 | 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md:10, 102-139, 284-306 | [fixed] 04장: 2n-1을 경험적 규칙으로 재서술(제목·목표·Q2), 분기 모형은 정성적 동기+반보존 복제 무시 한계 주석 | Smith 2017 PMID 28100584 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-037 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:382-387 (Q1 해설) | [fixed] 05장 Q1: coffset=0x00002A3C0000=708,575,232 (재계산), uoffset 4,176 확인 | SAM spec 4.1.1; 파이썬 재계산 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-038 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:274-283 | [fixed] 05장 reg2bin: 비교 시프트 >>14,17,20,23,26 (SAM 스펙, pysam bin=73으로 직접 검증), reg2bins O(1) 삭제 | SAMv1.tex; pysam 0.24.1 bin=73 직접 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-044 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:163-164, 186, 375 | [fixed] 06장: T878A chrX:67,723,710 (UCSC 서열 API에서 A, 코돈 ACT; 67,723,878=C 확인), 코드에 ref_fasta 참조 염기 검증 추가 | UCSC sequence API 직접 조회 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-051 | 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:107-111, 204, 225-235 | [softened] 07장: 3′ 라이브러리 정션 Ψ 한계 주석(AR-FL 정션 poly(A)에서 >=7.6kb), UB 필터 위험은 코드 주석(미확인 표기); Ψ 코드 자체는 유지 | UCSC ncbiRefSeq; Zheng 2017 PMID 28091601 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-052 | 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:380 | [softened] 07장: EPI-7386=NTD 억제제, ARV-110은 LBD 리간드 PROTAC로 AR-V7 근거 없음 명시(Plausible) | EPI-7386 Mol Cancer Ther 2019; ARV-110 PMID 39670468 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-056 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:305-306, 326-330 | [fixed] 08장 Q1: hg38 좌표(ERG exon4 38,445,404-621, TMPRSS2 exon1 41,508,081-158), − 가닥 breakpoint 공식(공여=POS, 수용=POS+L-1), 저장 서열 순서 ERG 먼저; 3.1 SA 예시·도식 방향 주석, 1.0절 도식 ERG/TMPRSS2 위치 정정(UCSC 조회) | UCSC ncbiRefSeq NM_005656.4, NM_004449.4 (직접 조회) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-057 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:306, 309, 332-334 | [fixed] 08장 Q1(3): 40+106≠150이므로 중복 아닌 4bp 삽입(N-region); 중복은 106S44M일 때(미세 상동성)로 분리 | SAM spec CIGAR 소비 산술 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-059 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:3, 37, 156-165, 365-371 (09:187-188, 198) | [fixed] 08장 §1.1·Q3, 09장 루브릭: inferCNV/CopyKAT로 21q22 3Mb Edel 판정은 해상도 밖(약 15유전자), DNA 근거 필요; "확증됨"→예시 | UCSC ncbiRefSeq 유전자 수; Track 04 SCI-04-032 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-060 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:100-102 | [already-fixed] 2단계 D(3′ 융합 포획) 적용 확인: 08장 §3.2 ERG exon4→poly(A) 4.7kb, 대리 지표 서술 | 2단계 D | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-061 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:377-380 | [fixed] 08장 Q3: SPOP-ERG 기전을 An 2015(절단 ERG의 degron 상실)로 정정, 야생형 ERG 안정화 서술 삭제, 상호배타성은 TCGA 아형 수준 | An 2015 PMID 26344096; TCGA PMID 26544944 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-063 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:167-169, 385-386 | [fixed] 08장 §1.2·Q3: 0카피 예후 단정 삭제, 문헌은 융합 중복+결실(2+Edel, Attard 2008 HR 6.10); HMGN1/DSCR3 근거 없음(DSCR3 구간 밖) | Attard 2008 PMID 17637754; UCSC | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-069 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:236-238, 253-258 | [fixed] 09장 Q2: R1/R2 역전 원인 삭제(82% 매핑과 불일치), 1순위=화이트리스트/chemistry, Barcodes.stats 확인 | STARsolo docs (Barcodes.stats) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-09-070 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:155 | [fixed] 09장 실패모드 07: UMI 충돌=과소계수, 처방은 UMI 길이/포화 보정; 해밍 병합은 오류용으로 구분 | Smith 2017; Fu 2011 PMID 21562209 | 제어문자 0, `$$` 짝 OK, python ast OK

다른 트랙의 같은 오류(미수정, 참고): AR 좌표·T878A 좌표 사본(Track 04, slides), TMPRSS2/ERG 방향·좌표 서술은 해당 트랙 단계에서 처리.

## S2 단계 / Track 01

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-01-002 | 00_MULTIOME_ROADMAP.md:131-141 | [fixed] 00장: 실재하지 않는 prostate_tumor URL(403) 삭제 → HEAD 200 확인한 PBMC granulocyte-sorted 10k 4개 URL로 교체, 표/디렉터리명 정정 | HTTP HEAD 직접 확인(200; fragments 2.4GB) | 제어문자 0, `$$` 짝 OK
SCI-01-003 | 00_MULTIOME_ROADMAP.md:183-198 (Q3 해설) | [fixed] 00장 Q3 해설: 엔잘루타마이드가 AR mRNA를 막는다는 설명 삭제(Chen 2004), 카피수 교란(Viswanathan 2018) 병기, primed 진단을 가설로, 병용요법은 가설 표기 | PMID 14702632, 29909985 | 제어문자 0, `$$` 짝 OK
SCI-01-005 | 01_EPIGENOMICS_AND_TN5_CHEMISTRY.md:82-88 | [fixed] 01장 §3: 10x 독립 scATAC도 핵 일괄 태그멘테이션 → 차이를 동일 GEM 내 RNA+ATAC 바코딩으로 재서술, C1류와 구분 | 10x CG000362 Getting Started scATAC | 제어문자 0, `$$` 짝 OK
SCI-01-009 | 02_UPSTREAM_AND_FRAGMENTS.md:130-133 | [fixed] 02장: KLK3 hg38 chr19:50,854,915-50,860,764(+), 질의 구간 50,853,000-50,857,000 | UCSC ncbiRefSeqSelect NM_001648.2 직접 조회 | 제어문자 0, `$$` 짝 OK
SCI-01-010 | 02_UPSTREAM_AND_FRAGMENTS.md:62-74 | [fixed] 02장 Chromap: --barcode-translate 추가(2열: 변환 후/원본), sort|bgzip, tabix -p bed 추가, 화이트리스트·번역 파일명은 자리표시자 | Chromap README 직접 확인 | 제어문자 0, `$$` 짝 OK
SCI-01-014 | 03_SNATAC_QC_METRICS.md:3, 57-72, 78-81, 193 (09:11, 58, 237, 270) | [fixed] 03·09장: 147 bp 간격 사다리 → 약 200 bp 주기(모노 180-220, 다이 380-420), 147은 코어 DNA 길이/경계 컷오프, 294 bp 삭제 | Buenrostro 2013 PMID 24097267; Signac vignette | 제어문자 0, `$$` 짝 OK
SCI-01-024 | 05_WNN_MULTIMODAL_INTEGRATION.md:84-88, 113 | [already-fixed] 2단계 C(WNN 가중 친화도 식) 적용 확인: 05장 §2.4 weighted affinity, 가중 제곱합 아님 |  | 제어문자 0, `$$` 짝 OK
SCI-01-025 | 05_WNN_MULTIMODAL_INTEGRATION.md:60-81 | [already-fixed] 2단계 C 적용 확인: 05장 §2.3 within/cross 커널 비+softmax, θ_M 없음(Seurat 소스) |  | 제어문자 0, `$$` 짝 OK
SCI-01-026 | 05_WNN_MULTIMODAL_INTEGRATION.md:96-108, 204-214 (Q3 해설) | [softened] 05장: 가중치 예시를 가상으로 표기, 입증→가설 생성용 지표, 임상 예측(AR 소실/PSA 위장/백금)은 직교 검증 요구 문구로 대체 | Hao 2021 PMID 34062119 | 제어문자 0, `$$` 짝 OK
SCI-01-029 | 06_PEAK_TO_GENE_AND_ENHANCERS.md:57-70, 162-167 | [fixed] 06장 Cicero: ρ=λ(1-(xmin/d)^s), xmin=1000, s=0.75, 가까울수록 페널티 작음; 기본 k=50, 윈도우 500kb, distance_constraint 250kb; 완전 제거→완화 | Cicero runCicero.R get_rho_mat; PMID 30078726 | 제어문자 0, `$$` 짝 OK
SCI-01-030 | 06_PEAK_TO_GENE_AND_ENHANCERS.md:121-141, 179-182 (00:184) | [fixed] 06장: AR 인핸서 중복 빈도 40-60%→mCRPC 70-87%/원발 <2%, AR 본체 2배체 일반화 삭제, 50배/10배/r>0.65/FDR 수치는 출처 없는 예시 표기 | Viswanathan 2018 PMID 29909985 | 제어문자 0, `$$` 짝 OK
SCI-01-031 | 06_PEAK_TO_GENE_AND_ENHANCERS.md:141 (00:184, 195) | [softened] 06장: snATAC 단편 수는 카피수에 비례하므로 접근성 증가·peak-to-gene 상관은 카피수 교란 가능, CN 보정 요구 병기 |  | 제어문자 0, `$$` 짝 OK
SCI-01-034 | 07_TF_FOOTPRINTING_AND_CHROMVAR.md:55-65, 168 | [fixed] 07장 TOBIAS: DWM 기반 편향, 보정=관측−기대(차감, 로그비 아님), footprint score=플랭크 양의 평균−코어 음의 평균 최대(TOBIAS 소스), CHG consensus 서술 삭제 | Bentsen 2020 PMID 32848148; TOBIAS signals.pyx 직접 확인 | 제어문자 0, `$$` 짝 OK
SCI-01-037 | 07_TF_FOOTPRINTING_AND_CHROMVAR.md:189-199 (Q3 해설) | [softened] 07장 Q3 해설: chromVAR는 모티프 피크 접근성 증감(점유율 아님), E-box는 TF 비특이, t-NEPC 완료 단정·치료 중단/백금/EZH2/DLL3 권고 삭제 | Schep 2017 PMID 28825706 | 제어문자 0, `$$` 짝 OK
SCI-01-040 | 08_EGRN_AND_SCENICPLUS.md:32-34, 126-133 | [fixed] 08장: SCENIC+ 기본 rho 임계 0.03 (r>0.3 삭제), 코드 주석 병기 | SCENIC+ TF_to_gene.py RHO_THRESHOLD=0.03 직접 확인 | 제어문자 0, `$$` 짝 OK
SCI-01-041 | 08_EGRN_AND_SCENICPLUS.md:65-92, 3 | [fixed] 08장: 2단계 L67 완화 위에 목표문·Q1 "3차원 물리적 인과 증명" → 후보 조절망 추론(인과성은 실험 검증 필요) | SCENIC+ PMID 37443338 | 제어문자 0, `$$` 짝 OK
SCI-01-044 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:85-87, 267-272, 286-290 | [softened] 09장: 과보정 원인을 환자=배치 교락으로 정정(LSI 1 포함 여부와 별개), 해결책은 theta 하향/종양 세포 비보정 | Harmony PMID 31740819; Signac vignette | 제어문자 0, `$$` 짝 OK
ENG-010 / ENG-011 (applied by sci) | tracks/01_multiome_snatac/docs/03_SNATAC_QC_METRICS.md:~115; docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:~30; tracks/07_count_statistics_abundance/docs/03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:~218-230 | ENG-010: snap.pp.import_data(genome=) → snap.pp.import_fragments(chrom_sizes=snap.genome.hg38, ...). ENG-011: dc.get_pseudobulk(..., min_cells, min_counts) → dc.pp.pseudobulk(adata, sample_col, groups_col, mode='sum') + dc.pp.filter_samples(pbulk, min_cells, min_counts)(in-place) | .venv_api에서 시그니처 재확인: snapatac2 2.10.0 (import_fragments(fragment_file, chrom_sizes, *, min_num_fragments, sorted_by_barcode, ...), import_data 없음, metrics.tsse(adata, gene_anno), pp.filter_cells(data, min_counts, min_tsse, max_counts)), decoupler 2.2.0 (dc.get_pseudobulk 없음; pp.pseudobulk, pp.filter_samples(adata, min_cells, min_counts, inplace=True)) | Track 07 문서 함수를 md에서 추출해 합성 AnnData 5000x2000으로 실행: 8 샘플x4 세포유형=32 pair 생성, 세포유형별 8 샘플 유지 / ENG 상태: ENG-010 fixed, ENG-011 fixed (s2_status_eng.csv는 수정하지 않음)

다른 트랙의 같은 오류(미수정, 참고): AR 650kb 인핸서 40~60% 빈도(Track 00·04 등), 147 bp 래더 서술(Track 02 snATAC 관련 문서가 있다면), WNN 식은 Track 00에서 이미 처리.

## S2 단계 / Track 02

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-02-002 | docs/00_SPATIAL_ROADMAP.md:98-100, 265 | [fixed] 00장 L98-100 표: "해상도"를 분해능/국소화 정밀도로 구분(광학 한계 ≈0.2 µm; CosMx/MERSCOPE <0.1 µm 삭제), 민감도 열은 정량 순위 아님으로 완화 | Abbe 한계 산술(01장 L107, 186 nm) | 제어문자 0, `$$` 짝 OK
SCI-02-004 | docs/00_SPATIAL_ROADMAP.md:80 | [fixed] 00장 L80, 01장 L39: Visium HD 빈 수 3×10^6/10.5×10^6 → 약 1.06×10^7 (6.5×6.5 mm ÷ 2 µm 격자 산술) | 산술 (6500/2)^2=1.056e7 | 제어문자 0, `$$` 짝 OK
SCI-02-008 | docs/01_SPATIAL_TECHNOLOGIES_AND_PHYSICS.md:239 | [fixed] 01장 L239: 수천 광자→10~20 nm 삭제, σ_PSF≈80 nm에서 N≈16~64 광자에서 10~20 nm, N=1000이면 이상적 2.5 nm | 산술; Thompson 2002 Biophys J | 제어문자 0, `$$` 짝 OK
SCI-02-009 | docs/01_SPATIAL_TECHNOLOGIES_AND_PHYSICS.md:241 | [fixed] 01장 L241: 라운드당 5~10% → MHD4(16비트·가중치4) 25%, Xenium은 별개 방식(RCA 색-순서 코드워드) | Chen 2015 Science PMID 25858977; 산술 4/16 | 제어문자 0, `$$` 짝 OK
SCI-02-020 | docs/03_SPATIAL_QC_AND_PREPROCESSING.md:56-57, 215 | [fixed] 03장: "단백질 코딩만 표적이라 MT 제외" 논리 삭제, Visium 프로브 세트 v1.0 MT 0개 / v2.0 12개(MT-ATP6 등), 사용한 CSV로 확인 지침 | 10x 프로브 CSV 직접 조회(v1.0 |MT- 0건, v2.0 12종) | 제어문자 0, `$$` 짝 OK
SCI-02-021 | docs/03_SPATIAL_QC_AND_PREPROCESSING.md:152, 167, 243-245 | [softened] 03장: 전역 총 UMI 분위수/3,500 상한은 접힘 증거가 아님(순환·종양 코어 제거 위험) 주석, H&E 마스크·국소 이상치 요구, 3,500은 예시 | 논리 점검 (Plausible) | 제어문자 0, `$$` 짝 OK
SCI-02-026 | docs/04_CELL_SEGMENTATION_ALGORITHMS.md:48-54, 226 | [fixed] 04장: Cellpose flow를 정상상태 라플라스 해가 아닌 중심 시작 확산 시뮬레이션(훈련 시 반복 2×extent)으로 정정, 추론 스텝 20~50 → niter 기본 200 | Stringer 2021 PMID 33318659; cellpose dynamics.py 직접 조회 | 제어문자 0, `$$` 짝 OK
SCI-02-031 | docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:118-128 | [fixed] 05장 cell2location 식: μ=(m_g Σ w g + s_{e,g})·y_s, m_g 유전자 수준 곱셈 보정, s_{e,g} 배치별 가산 배경, y_s 스팟별 검출 민감도, ε 항 삭제, 사전분포 기호 정정 | cell2location 소스 module; PMID 35027729 | 제어문자 0, `$$` 짝 OK
SCI-02-032 | docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:87-92, 94-104 | [fixed] 05장 RCTD: M-추정 문장 삭제, γ_j(플랫폼 효과)와 ε_{s,g} 랜덤 효과를 가진 Poisson-lognormal, "Robust"=플랫폼 효과 강건성 | Cable 2022 PMID 33603203 | 제어문자 0, `$$` 짝 OK
SCI-02-033 | docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:47, 329-330 | [fixed] 05장·09장: r>0.95의 κ≫10^5 삭제 → κ=(1+r)/(1-r) (r=0.95에서 39, 0.99에서 199), 랭크 결손은 r=1에서만 | 산술(2×2 상관 행렬 고유값 1±r) | 제어문자 0, `$$` 짝 OK
SCI-02-037 | docs/06_SPATIALLY_VARIABLE_GENES.md:144, 162, 305-307 | [fixed] 06장: SpatialDE는 길이 척도당 고유분해 1회 후 유전자별 O(N) (FaST-LMM식), "수 시간~수일" 삭제(N=5000 eigh 약 19 s 측정), 메모리 수백 TB → 약 8 TB(float64) | Svensson 2018 PMID 29553579; SpatialDE 소스; numpy eigh 측정 | 제어문자 0, `$$` 짝 OK
SCI-02-046 | docs/07_CELLULAR_NICHES_AND_MICROENVIRONMENTS.md:263-273 (Q1), 80-90 | [fixed] 07장: 전역 레이블 순열은 구획 구조를 통제하지 못함 주석(시뮬레이션 Z=+23.9 vs 구획 내 -0.39), 해석을 공간 공존으로 제한, Q1 해설 "엄밀 분리" 삭제 | 재현 시뮬레이션 nhood_sim.py | 제어문자 0, `$$` 짝 OK
SCI-02-047 | docs/07_CELLULAR_NICHES_AND_MICROENVIRONMENTS.md:144-146, 157, 300-313 | [fixed] 07장 니치 4·Q3, 08장 §5.1: 용골성 → 혼합형·골형성 우세, t-NEPC RANKL 분비 단정 삭제(RANKL 주로 골아세포계, PTHLH는 전립선암 세포 발현 보고), WNT5A 직조골 단정 → 가설, 07·08장 통일 | Yin 2003 PMID 12941866; PMID 29440427 | 제어문자 0, `$$` 짝 OK
SCI-02-052 | docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:78-92 | [fixed] 08장 §3: Squidpy 귀속 삭제(ligrec는 CellPhoneDB식 클러스터 레이블 순열, 공간 그래프 미사용), 일반 공간 순열 아이디어로 재표기, 실제 공간 도구(stLearn, CellChat v2, COMMOT) 명시, CellChat v2 병기 | squidpy gr/_ligrec.py; PMID 35102346 | 제어문자 0, `$$` 짝 OK
SCI-02-053 | docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:40-52 | [fixed] 08장 §2: 정상상태 해를 차원별(1D exp, 3D Yukawa, 2D K0)로, λ=√(D/k), 가우시안은 과도해, d0 vs 200~250 µm(Francis & Palsson 최대 영역) 구분 | 표준 PDE 해; Francis & Palsson 1997 PMID 9356436 | 제어문자 0, `$$` 짝 OK
SCI-02-054 | docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:46-70, 201, 274; SVG06; tracks/02_spatial | [already-fixed] 08장 §2 참고문(합성 파이프라인 d0=30 µm 비교 불가)과 SVG06/파이프라인은 1·2단계 eng 처리분; 본문 d0 서술은 위 053과 통일. SVG·파이프라인 파일은 수정하지 않음 |  | 제어문자 0, `$$` 짝 OK
SCI-02-066 | pipelines/results_pipeline_01/deconvolution_benchmark_metrics.csv:- | [deferred] 생성 산출물(results_pipeline_01 CSV)과 파이프라인 스크립트 docstring은 sci 범위 밖(생성 파일 수정 금지); eng ENG-033이 합성 데이터 고지(SIMULATION_NOTICE.md)를 추가함. 합성 sanity check·Spearman 병기·"Cell2location-style" 문구는 eng 파이프라인 정리 시 처리 권고 |  | 제어문자 0, `$$` 짝 OK

다른 트랙의 같은 오류(미수정, 참고): SVG01/06 수치·d0 표기와 파이프라인 02 docstring은 eng 소관. 본 단계로 sci S2 전 트랙(00,07,04,09,01,02) 완료.

## S2 단계 / Track 06 (재배정)

(원래 eng 담당 29건을 sci가 처리. 추적 파일: `review/s2_status_sci_t06.csv`; 수정 범위는 `tracks/06_manifold_trajectory/docs/*.md`만, SVG·파이프라인 코드·`s2_status_eng.csv` 미수정.)

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-06-001 | docs/00:171-176 | [fixed] 00장 §4.2: Replogle 2022는 인간 세포주 genome-scale CRISPRi Perturb-seq(전립선·NE 분석 아님), 방법론 사례로만 인용, AR/FOXA1/TP53/RB1→NE 도약은 가상 시나리오 표기 | PMID 35688146 초록 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-002 | docs/00:165-170 | [fixed] 00장 §4.1: Karthaus 2020 초록 수준(관강세포 재생 기여)으로 한정, 벨로시티/WOT "탈분화 벡터 증명"은 가상 적용 시나리오, Day 0-28 일정은 교육용 예시 | PMID 32355025 초록 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-003 | docs/04:319-321 | [fixed] 04장 문제 2: Karthaus/임상 시계열 인용 삭제 → 가상 시나리오 명시, Day 28과 모의 데이터(Day 0-14) 불일치 주석 | PMID 32355025; docs/04 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-004 | docs/07:14,163; docs/08:109-110; docs/00:102 | [already-fixed] 2단계 G 적용 확인: 07장 L14 흡수=모델 근사/비가역 아님, 08장 "(모델상) 고착", 00장 L102 모순 해소 |  | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-007 | docs/01:10,194-200,300,313-315,349,363 | [fixed] 01장: DPT 거리를 누적 확률 차이의 노름 D²=Σ(λ/(1-λ))²(ψ(x)-ψ(y))²로 정정(§5.2 유도, 목표문, Q 해설, 특이성 논의 동일 수정) | scanpy neighbors/__init__.py _get_dpt_row 직접 확인; Haghverdi 2016 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-010 | docs/02:201-204; docs/03:212; docs/05:164 | [fixed] 02·03·05장: scv.pp.filter_and_normalize(n_top_genes=) TypeError → min_shared_counts만 사용 + sc.pp.highly_variable_genes(subset=True) | scvelo 0.3.4 시그니처·실행 TypeError 및 대체 순서 실행 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-011 | docs/03:249,259-264 | [fixed] 03장: E-step "Brent 1차원 볼록 탐색"·"그래프 사전분포(λΣW(t-t_j)²)" 삭제 → assign_tau 격자 탐색(argmin), 세포 독립 할당; Q1(2) 질문·해설 재작성 | scvelo _em_model_utils.assign_tau 소스 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-012 | docs/03:178-184,188 | [fixed] 03장 §4.2: 가중식 w_g=logL·C를 개념식으로 격하, 실제 latent_time 절차(min_likelihood 0.1, root_key/속도 그래프 추정, compute_shared_time 백분위 선형성 2개 합산, 최대 4 root 평균) 병기, "완전 자율" 삭제 | scvelo _em_model_core.latent_time, compute_shared_time 소스 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-013 | docs/03:278-281 | [fixed] 03장 문제 2: 부호 교정 — fast-γ 클론은 u>γ̄ s, v=(γ_fast-γ̄)s>0 거짓 양의 벨로시티, slow 클론은 거짓 억제 | 대수 유도 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-015 | docs/04:236 (ENG-014 참조), 252-258, 269 | [already-fixed] ENG-014로 OTModel(adata, day_field, growth_rate_field, ...)·compute_transport_map 적용 확인(wot 1.0.x 시그니처 재확인); 성장률 난수 폴백은 코드 주석에 교육용 모의값으로 명시됨 | wot OTModel.__init__ 시그니처 및 compute_transport_map 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-017 | docs/05:177-179 | [fixed] 05장: VelocityKernel.compute_transition_matrix temperature → softmax_scale | cellrank 2.3.3 시그니처 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-018 | docs/05:185-189; docs/09:148-149 | [fixed] 05장: compute_cytotrace(layer="Ms"), compute_transition_matrix(threshold_scheme="soft", nu=0.5) (threshold_low/high 없음) | cellrank 2.3.3 시그니처 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-019 | docs/06:191; docs/09:164 | [fixed] 06·09장: coarse_grained_transition_matrix → coarse_T; Reuter 서지를 JCTC 2018(PMID 29812922)·JCP 2019(PMID 31067901)로 정정, CellRank 2 서지를 Nat Methods 2024(PMID 38871986)로 정정 | dir(GPCCA)에 coarse_T 확인; PubMed esummary | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-020 | docs/06:46-51,208-217 | [fixed] 06장: "세부균형 위반⇒비정규" 논리 정정 — 세부균형은 π-가중 자기수반, 가역 P도 일반적으로 PPᵀ≠PᵀP, 순환 치환은 정규이면서 위반; 에르미트 서술도 π-가중으로 한정 | 수치 반례(가역 3상태, 순환 치환) 및 선형대수 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-021 | docs/05:141-145,240-250; docs/09:28,49-54,186-193 | [fixed] 05·09장: v<0은 루미날 유전자 발현 감소=루미날 이탈 정방향임을 정정, 거짓 역류 원인을 하향 유전자 지배/γ·u/s 추정 오류로 재서술 | 대수: v=βu-γs<0 ⇒ ds/dt<0 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-022 | docs/07:170-174,306-316 | [fixed] 07장 엔트로피 표·Q3: 4번째 계통(Amphicrine) 잔여로 재계산 0.30/1.68/0.51/0.22 bits, Q3 Group 1 (0.44,0.42,0.14)=1.44 bits(3결과 최대 1.58), "2.0 bits 근접" 삭제 | numpy 재계산 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-025 | docs/08:152-175 | [fixed] 08장: compute_lineage_drivers(lineage=, cluster_key=) → lineages= (cluster_key는 clusters와 함께일 때만) | cellrank 2.3.3 시그니처 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-026 | docs/08:178-200 | [fixed] 08장: gene_trends/heatmap에 time_key 필수 추가, heatmap show_absorption_probabilities → show_fate_probabilities, GAMR은 R(mgcv)+rpy2·기본 gaussian(NB는 distribution="nb"+raw count) 주석 | cellrank 2.3.3 pl.gene_trends/pl.heatmap/GAMR 시그니처 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-033 | docs/08:122-137,238; SVG06 | [already-fixed] 08장 §4가 가상 시나리오·파이프라인(r≈0.58)과 다름을 명시(1·2단계); SVG06은 eng가 Schematic 표기 |  | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-034 | docs/08:282-287 | [softened] 08장 Q3: DLL3 판정 기준을 교육용 가상 기준으로 명시(임상 근거 없음, 임상 선별 금지), "반응 전무" 단정 완화 | 임상 근거 없음 확인(Plausible) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-038 | docs/09:99,262 | [fixed] 09장: SpatialKernel 삭제 → PrecomputedKernel(공간 그래프 기반 전이행렬) | dir(cr.kernels)에 SpatialKernel 없음 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-039 | docs/09:261 | [fixed] 09장: scVI(RNA 전용) → MultiVI(scvi.model.MULTIVI) 또는 WNN/MOFA+ | scvi-tools MULTIVI 존재 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-040 | docs/09:158-168 | [fixed] 09장 코드: compute_macrostates 뒤 predict_terminal_states() 추가(없으면 RuntimeError) | cellrank 2.3.3 시그니처·동작 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-048 | pipeline02:252-308; results_pipeline_02/macrostate_assignments.csv | [softened] 08장 §4 주석: 합성 파이프라인 02의 "G-PCCA"는 입력 계통 라벨 복사임을 명시(파이프라인 코드는 수정하지 않음) | pipeline02 L301 라벨 분기(findings) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-049 | pipeline02:395-405; results_pipeline_02/cell_fate_probabilities.csv, pipeline_02_summary.txt | [softened] 08장 §4 주석: P(t-NEPC)=P(DNPC)=0.5는 흡수 경계 유속 0일 때의 폴백 인공물(파이프라인 코드는 수정하지 않음) | pipeline02 L405 폴백(findings) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-050 | pipeline02 (expression L~140-175); results_pipeline_02/lineage_driver_correlations.csv | [already-fixed] 1·2단계(ENG-033 SIMULATION_NOTICE 및 08장 §4 주석)에서 합성 결과 고지; 이번에 순환 회수·비-TF 포함 주석 추가 | results_pipeline_02/SIMULATION_NOTICE.md 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-051 | pipeline01:335,518; results_pipeline_01/cell_trajectories_and_velocity.csv, pipeline_01_summary.txt | [softened] 08장 §4 주석: 파이프라인 01의 DPT/latent time은 Haghverdi DPT·scVelo 동역학 잠재 시간이 아닌 단순화 근사임을 명시(코드는 수정하지 않음) | pipeline01 L335, L518(findings) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-052 | SVG06 (텍스트); results_pipeline_02/pipeline_02_metrics.json | [already-fixed] SVG06 하드코딩 수치: eng가 2단계에서 Schematic/예시 표기(텍스트에 "Schematic", "개념도", "예시" 확인); 파일은 수정하지 않음 |  | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-06-053 | SVG01 (텍스트) | [already-fixed] SVG01 M 식이 eng 수정으로 M=Σ_t(T^t-φ0φ0ᵀ)로 교정됨 확인; 파일은 수정하지 않음 |  | 제어문자 0, `$$` 짝 OK, python ast OK

다른 트랙의 같은 오류(미수정, 참고): Track 05 등의 DPT 가중치 표기, scVelo n_top_genes 인자 사용은 해당 트랙에서 별도 확인 필요.

## S2 단계 / Track 03 (재배정)

(원래 eng 담당 16건을 sci가 처리. 추적: `review/s2_status_sci_t03.csv`; 수정 범위 `tracks/03_tme_cell_communication/docs/*.md`만. 편집 전 pre-S2 스냅샷과 비교한 결과 Track 03 docs에 eng의 S2 수정은 없었음(diff 0). SVG·생성 스크립트·`s2_status_eng.csv` 미수정, ENG-004/005/006 코드 유지.)

ID | 파일:줄 | 전→후 | 근거 | 검증
--- | --- | --- | --- | ---
SCI-03-001 | docs/00_...:199-203; 04_...:322-324; 05_...:376-381; 07_...:370-372; 08_...:141-153,324-332; images/06_...svg | [fixed] 00·04·05·06·07·08장: 미승인·중단 약제에 "탐색적 가설" 주의문 추가(bintrafusp alfa 미승인·INTR@PID BTC 055 중단, magrolimab ENHANCE 무용성·임상보류 후 중단, cabozantinib은 ROR2 표적 아니고 COMET-1 OS HR 0.90 P=0.213, galunisertib 상태 미확인, Denosumab 근거는 SRE 지연/nmCRPC 골전이-무진행생존), "원천 봉쇄/완전 사멸/완벽 살상" 완화. images/06 SVG는 수정하지 않음 | Fizazi 2011 PMID 21353695; Smith 2012 PMID 22093187; EMD Serono 2021; COMET-1 2015 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-002 | docs/00_...:69-71 | [fixed] 00장 표: Bill et al.을 Science 2023;381:515 (PMID 37535729, CXCL9:SPP1 다암종)로 정정, 이전 accession·65,000 cells 근거 없음 표기 (Kfoury/Chen 행은 이전 단계에서 정정됨 확인) | PubMed esummary | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-009 | docs/02_...:3,72-118,208-221; 06_...:161-166; images/02_...svg; pipelines/01:491-503 | [fixed] 02장 목표·파이프라인 박스·§3.1: LIANA 합의는 RRA(옵션 mean), magnitude/specificity 2축; Borda는 비LIANA 일반 이론으로 분리. 06장 §5.1 동일 정정. SVG02에 Borda 표기 잔존(미수정) | liana rank_aggregate 시그니처(aggregate_method rra|mean) 직접 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-010 | docs/02_...:58-60; 09_...:9; 06_...:159 | [fixed] 02·09·06장: Dimitrov 2022는 Nat Commun 13:3224 (PMID 35680885), LIANA+는 Nat Cell Biol 2024;26:1613 (PMID 39223377); 25~35%/60% 삭제 → 논문 수치(평균 고유 상호작용 10.4%, Jaccard 중앙값 0.080~0.119) | PubMed esummary; PMC9184522 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-011 | docs/02_...:230-247 | [fixed] 02장 문제 3: PTPRC–CD22는 cellphonedb·cellchatdb에 각 1행, lrdb·consensus·icellnet·guide2pharma에는 없음으로 전제 정정, LRdb=텍스트마이닝 서술 삭제, 구획 정합성 중심으로 재서술 | liana 1.10.0 select_resource 직접 조회 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-013 | docs/03_...:3,63-88,151; 00_...:161 | [fixed] 00·03장: 점수 정의를 도구별로 구분(CellPhoneDB=산술평균·복합체 최소값, NATMI/Connectome=곱, CellChat=Hill·기하평균), 곱 점수 코드는 교육용 변형으로 명시, CellPhoneDB p는 +1 보정 없음 병기 | liana _cellphonedb.py; CellChat 논문 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-017 | docs/04_...:42-70,152-182; images/03_...svg | [fixed] 04장: CellChat 보조인자 식을 구조식(Hill × 공수용체 선형 × agonist/antagonist Hill [× 세포 비율])으로 교체, (1+A)/(1+I)는 원식 아님 명시, Kh=0.5 고정 상수·trimean, P는 점수이며 1 초과 가능 주석. 예시 값 범위 주의. SVG03 미수정 | CellChat Nat Commun 2021 Methods 본문; liana _cellchat.py (kh=0.5, trimean) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-018 | docs/04_...:74-110,201-236,292-301; images/03_...svg | [fixed] 04장: 4대 중심성은 out-degree, in-degree, flow betweenness, information centrality(R sna)임을 명시, 코드의 betweenness/eigenvector는 근사로 표기, PageRank 폴백 → NaN. SVG03에 PageRank/정보 중심성 표기 잔존(미수정) | CellChat 논문 Methods | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-021 | docs/05_...:58-60,98-117,125 | [fixed] 05장: 가중치 최적화=mlrMBO(경사하강법 아님), ligand→TF는 PPR(감쇠 0.5, 상위 1%) 후 TF→표적은 GRN 행렬곱(두 번째 RWR 아님), "수천 건" 삭제 | nichenetr construct_ligand_to_target.R, parameter_optimization.R (findings 소스 확인분) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-023 | docs/00_...:57; 05_...:184-189,361-373; 07_...:168-171; 08_...:16,25,59 | [fixed] 00·05·07·08장: WNT5A 방향 통일 — 골아세포 니치 Wnt5a는 ROR2/SIAH2로 전립선암 휴면 유도(Ren 2019)이므로 종양 분비 단정 삭제, NFATc1→ASCL1/PEG10 캐스케이드 가설 표기, RANKL 공급원 맥락 의존 | Ren 2019 JEM PMID 30593464 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-024 | docs/05_...:361-373; 07_...:166-171; 08_...:50-56,270-275 | [softened] 05·08장: TGF-β→SMAD3→SOX2/ASCL1→t-NEPC 사슬을 근거 없는 가설로 표기하고 문헌상 드라이버(RB1/TP53 소실 Ku 2017, SRRM4-REST Li 2017, EZH2) 병기 | PMID 28059767, 27180064 (findings) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-026 | docs/06_...:161-166 | [fixed] 06장 §5.1: 합의는 RRA 기반, 공간 정보는 합의 점수 곱이 아니라 rank_aggregate(spatial_key)의 메서드별 근접 가중치와 별도 li.mt.bivariate로 정정 | liana rank_aggregate 시그니처(spatial_key) 직접 확인 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-032 | docs/07_...:32 | [fixed] 07장: MISTy 저자·연도·약어 — Tanevski et al., Genome Biol 2022;23:97 (PMID 35422018), Multiview Intercellular SpaTial modeling framework | PubMed esummary | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-035 | docs/07_...:157-171,350-373; images/05_spatial_misty_multiview_modeling.svg | [softened] 07장: MISTy 수치를 "가상 예시(계산 결과 아님)"로 표기, 임상 결론·삼중 병용은 가설로 약화, 코드는 단순화 재구현(in-sample 메타모델, gain.R2와 정의 상이) 주석. SVG05의 Intrinsic 45%/Spatial gain 55% 잔존(미수정) | SVG 텍스트; 본문에 수치 산출 코드 없음 | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-036 | docs/08_...:92-107,282-299 | [fixed] 08장: "고유치 우반평면 수학적 증명" 삭제 → 양의 폐루프 되먹임 존재일 뿐 불안정은 파라미터 의존(분해항 −δ; 특성방정식 g>abc), 수치 반례(안정 정상상태, T 0.90~0.99 포화) 병기, "K_T 도달 파멸 루프" 입증 삭제 | findings ODE 수치 실험(ode08.py) | 제어문자 0, `$$` 짝 OK, python ast OK
SCI-03-039 | docs/09_...:154-178,180-201,206-210 | [fixed] 09장 코드: 컬럼 ligand.complex/receptor.complex/aggregate_rank → ligand_complex/receptor_complex/magnitude_rank/specificity_rank, 유사벌크를 실제 입력으로 사용(합산 카운트→정규화→rank_aggregate) | liana 1.10.0 출력 컬럼 및 합성 유사벌크 실행(793행 specificity_rank<0.05) 확인 | 제어문자 0, `$$` 짝 OK, python ast OK

미반영 잔존(범위 밖, 참고): SVG02 'Borda', SVG03 'PageRank / 정보 중심성', SVG05 'Intrinsic 45% / Spatial Gain +55%' 텍스트는 eng 소관으로 수정하지 않음.

## S3 단계 / Track 00

(`s3_status_sci.csv` 기준 45건. 수정은 root `docs/*.md`만; 편집 후 제어문자 0, `$$` 짝 OK, python ast OK.)

ID | 파일:줄 | 조치 | 상태
--- | --- | --- | ---
SCI-00-005 | docs/00_CURRICULUM_AND_ROADMAP.md:113 | 00장 Song 행: Seq-Well S³(10x 아님), 생검·절제 혼합 병기 | fixed
SCI-00-009 | docs/00_CURRICULUM_AND_ROADMAP.md:101 | 00장: 배수체→이배체, inferCNV "필수"→활용+직교 근거 | fixed
SCI-00-010 | docs/00_CURRICULUM_AND_ROADMAP.md:13 | 00장: 린치→전이성/다기관, 냉동=생세포 불가에 DMSO 보존 예외 단서 | fixed
SCI-00-011 | docs/00_CURRICULUM_AND_ROADMAP.md:59-73 | 00장: 환경 코드에 버전 고정·Python 요구(scanpy>=3.12 등)·pySCENIC 분리 권장 주석 (environment.yml 신규 작성은 범위 밖) | softened
SCI-00-018 | docs/01_FOUNDATIONS_AND_WET_LAB.md:42 | 01장: AR-V7을 결실이 아닌 대체 스플라이싱(CE3), 항상성→구성적 | fixed
SCI-00-019 | docs/01_FOUNDATIONS_AND_WET_LAB.md:117 | 01장: NP-40 콜레스테롤 선택성 삭제, 조건 의존 서술 | fixed
SCI-00-020 | docs/01_FOUNDATIONS_AND_WET_LAB.md:137-138 | 00·01·02장: 인트론 40~60%/절반 이상 → 대략 수십 %(조직 의존), 리드 기준 단위 통일 | fixed
SCI-00-021 | docs/01_FOUNDATIONS_AND_WET_LAB.md:68, 163 | 용혈→세포 용해: 01장은 S2에서, 02장(L99,186)은 이번에 수정 | fixed
SCI-00-022 | docs/01_FOUNDATIONS_AND_WET_LAB.md:33, 36, 44 | 01장: 극심한 전이능→높은 전이 위험, 18~24개월→대략 1~2년 전후, amphicrine 과도기는 가설 표기 | softened
SCI-00-026 | docs/02_UPSTREAM_PROCESSING.md:74-89 | 02장 EmptyDrops: Cell Ranger는 EmptyDrops_CR, STARsolo 기본은 CellRanger2.2, 파라미터(45000-90000, umiMin 500, FDR 0.01 등) 정정 (STAR parametersDefault 확인) | fixed
SCI-00-027 | docs/02_UPSTREAM_PROCESSING.md:151 | 02장: 36601 출력을 레퍼런스 버전 의존 예시 출력으로 표기 | softened
SCI-00-029 | docs/02_UPSTREAM_PROCESSING.md:223 | 02장: n_genes<600 등을 예시값·03장 MAD와 다른 임의 기준으로 명시 | fixed
SCI-00-033 | docs/03_QC_AND_CLEANING.md:55, 185 | 03장: FOLH1(PSMA) 분비 단백질 서술 삭제(막 결합), 수천 배·수만 카피 수치 삭제 | fixed
SCI-00-035 | docs/03_QC_AND_CLEANING.md:139-142 | 03장 체크리스트: 고정 컷오프 → MAD 기반+안전 하한/상한, 수치는 예시 | fixed
SCI-00-036 | docs/03_QC_AND_CLEANING.md:79-85 | 03장 CellBender: --cuda 추가, expected/total droplets 설명, 기본값(fpr 0.01, epochs 150) 명시 | fixed
SCI-00-037 | docs/03_QC_AND_CLEANING.md:28, 209 | 03장: 다배수체→이수체, 이수체→전사활성 인과 단정 삭제(예시 수치) | fixed
SCI-00-040 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:12 | 04장: 1e4는 Seurat 기본, Scanpy normalize_total 기본은 median(target_sum=None) | fixed
SCI-00-041 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:33 | 04장: HVG 메커니즘을 seurat/cell_ranger(binned)와 seurat_v3(VST)로 구분 | fixed
SCI-00-042 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:74-75 | 04장: elbow 30~50 규칙 삭제 → 데이터 의존(jackstraw 등), 30~40은 예시 | fixed
SCI-00-043 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:113-118 | 04장 코드: scale·pca 명시, leiden flavor="igraph", n_iterations=2, directed=False | fixed
SCI-00-045 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:127 | 04장: Wilcoxon "zero-inflation에 가장 안정적" → 마커 탐색에서 경험적으로 양호, 조건 간 비교는 pseudobulk | fixed
SCI-00-046 | docs/04_CORE_DOWNSTREAM_PIPELINE.md:107, 183 | 04장: Louvain "치명적/빈번" → Traag 2019 수치(최대 25% badly connected, 16% disconnected), Leiden 희귀 세포 인공물 차단 과장 삭제 | fixed
SCI-00-048 | docs/05_CELL_TYPE_ANNOTATION.md:85-86 | 05장: infercnvpy 절차를 cnv.tl.pca→pp.neighbors→tl.leiden으로 정정 | fixed
SCI-00-049 | docs/05_CELL_TYPE_ANNOTATION.md:45, 70-79 | 05장: infercnvpy 기본 chrX/chrY 제외로 AR 증폭 평가 불가, FISH/WES 확인, ~60% 단정 삭제 | fixed
SCI-00-050 | docs/05_CELL_TYPE_ANNOTATION.md:40 | 05장: NKX3-1 8p 결손 ~70% → arm-level 호출 기준 60~80%(데이터셋 의존), FISH/array 보고는 낮음 | fixed
SCI-00-055 | docs/05_CELL_TYPE_ANNOTATION.md:100-102 | 05장: Numbat(cellsnp-lite+Eagle2)과 BAFExtract(CaSpER 계열) 분리 서술 | fixed
SCI-00-057 | docs/05_CELL_TYPE_ANNOTATION.md:22 | 05장 TAM 표: M2-like 표현을 연속적 상태·M1/M2 단순화로 완화 | softened
SCI-00-061 | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:26 | 06장: scVI 기본 gene_likelihood는 zinb, nb는 옵션 병기 | fixed
SCI-00-062 | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:68-72, 182 | 06장: iLISI "2.0 이상" 삭제, ARI/Isolated label F1을 scIB 정의로 정정(Q1 해설 포함) | fixed
SCI-00-064 | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:112 | 06장: scVelo 선행 조건(spliced/unspliced 레이어; STARsolo Velocyto) 추가 | fixed
SCI-00-067 | docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:217 | 06장: 선천적→기존(pre-existing) | fixed
SCI-00-069 | docs/07_MULTIOMICS_AND_FUTURE.md:40, 124 | 07장: 147 bp 주기성 → ~200 bp 간격 래더(147은 코어 DNA 길이) 2곳 | fixed
SCI-00-070 | docs/07_MULTIOMICS_AND_FUTURE.md:31-32 | 07장: FOXA1 인과 직접 입증 → 공접근성-발현 상관 기반 가설 생성 | fixed
SCI-00-071 | docs/07_MULTIOMICS_AND_FUTURE.md:98, 137 | 07장: RCTD를 비베이지안 MLE로 정정, 디콘볼루션 필수→권장 | fixed
SCI-00-072 | docs/07_MULTIOMICS_AND_FUTURE.md:78, 98, 136 | 07장: Spot당 5~10/5~15 → 대략 5~15(조직 밀도 의존)로 통일 | fixed
SCI-00-073 | docs/07_MULTIOMICS_AND_FUTURE.md:99 | 07장: Visium HD는 올리고 바코딩 2 µm 픽셀, 단일 세포는 분할·집계 필요 | fixed
SCI-00-074 | docs/07_MULTIOMICS_AND_FUTURE.md:100 | 07장: <0.1 µm는 국소화 정밀도(회절 한계와 구분), 패널 수치는 버전 의존 | fixed
SCI-00-075 | docs/07_MULTIOMICS_AND_FUTURE.md:89, 157-160 | 07장: 병용 필수→연구 중 가설, 10 µm 시냅스 기준 삭제, 면역 사막형도 병기 | softened
SCI-00-078 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:20-21 | 08장: N=1/N=2 다이어그램 → 군당 ≥3명 명시 | fixed
SCI-00-079 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:30-36 | 08장: dc.pp.pseudobulk는 S2(ENG-011)에서 적용 확인, layer="counts" 지정 추가 | fixed
SCI-00-080 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:57-62 | 08장: design_factors → design="~ ..." 문자열, pbulk_counts/sample_metadata 정의 추가 | fixed
SCI-00-081 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:80-81 | 08장: 고정효과 "강력 억제" → 완전 교란 시 분리 불가(06장 §1.3), 부분 겹침도 비보장 | fixed
SCI-00-082 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:95-99 | 08장: 참조 선택의 장단점으로 서술(정상 상피 참조는 benign 확신 시 장점) | softened
SCI-00-083 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:105 | 08장: OptiPrep 40/30/20% → 35/29%(01장 프로토콜)로 통일 | fixed
SCI-00-084 | docs/08_TROUBLESHOOTING_AND_BEST_PRACTICES.md:126-130 | 08장 루브릭: KLK3<1%/iLISI>1.8/cLISI<1.2/F1>0.85 임계값을 예시·보고 항목으로 격하, DNPC 레귤론 합격기준 삭제, 재현성 100%→허용오차 내 | softened

트랙 전체 일괄 치환 없이 챕터 내 일관 수정만 수행. 같은 오류 잔존 참고: 인트론 40~60% 서술(Track 02 등), 147 bp 주기성(Track 01 일부는 S2에서 처리).

## S3 단계 / Track 09

(43건. 수정은 `tracks/09_upstream_alignment_bam/docs/*.md`만; 제어문자 0, `$$` 짝 OK, python ast OK. 04장 방향성 그래프 코드는 실행해 3분자 복원 확인.)

ID | 파일:줄 | 조치 | 상태
--- | --- | --- | ---
SCI-09-007 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:87, 96-97 | 00장: CE3를 인트론 3 내 대체 3′ 말단 엑손으로 표기, 3b·101bp 삭제, 말단 엑손 교체 서술 | fixed
SCI-09-008 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:220-244 | 00장 트리: 05-09장 실제 파일명·내용으로 교정 | fixed
SCI-09-009 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:288 | 00장: N 개수별 후보 4^k × 1.6e-3 정량 서술 | fixed
SCI-09-010 | 00_ALIGNMENT_AND_BAM_ROADMAP.md:139-140 | 00장: 캡처율 수치 삭제, 인트론 리드 원인을 분리 서술 | fixed
SCI-09-015 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:204-206 | 01장: Cell Ranger 2.x 규칙(99백분위/10)·barcodeRanks knee·3.0+ EmptyDrops_CR 구분 | fixed
SCI-09-016 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:24-41, 35 | 01장: 라이브러리 도식을 개념도로 표기, 브리지 증폭→ExAmp 병기 (도식 자체는 재작성하지 않음) | softened
SCI-09-017 | 01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md:264-285 | 01장 코드: 메모리 한계·대안(정수 인코딩/C 확장/STARsolo) 주석 추가, 코드는 교육용 유지 | softened
SCI-09-018 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:7, 29-30, 39, 267-270 | 02장: SA O(m log n)(LCP O(m+log n)), FM O(m)로 표·목표·Q 통일 | fixed
SCI-09-021 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:107-127 | 02장: 락-프리/비트마스크/온라인 dedup을 개념도로 표기, 공식 구조 서술 | softened
SCI-09-022 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:311-315 | S2(SCI-09-019)에서 alignSJ* 기본값 명시 완료 확인 | already-fixed
SCI-09-023 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:292 | 02장: SYP/CHGA 크기 정정(약 12kb), NCAM1 317kb, GeneFull 필수→민감도 향상 | fixed
SCI-09-024 | 02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md:40-41 | 02장: FM-index 자체가 스플라이스 부적합이 아님(HISAT2), Bowtie2/BWA는 비스플라이스 정렬기로 서술 | fixed
SCI-09-030 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:286-303 | 03장: 해상도 전략명 정정(Tiered 없음), 스크립트 -r cr-like-em, 3′ 5′UTR 고유 리드 비현실 주석, STARsolo 다유전자 UMI 서술 | fixed
SCI-09-031 | 03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:140-145, 194 | 03장: flank 기본 L−5(pyroe)와 L−1 논증 구분, 스크립트 주석 정정 | fixed
SCI-09-034 | 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md:162-187 | 04장 코드: 양방향 간선 검사로 수정, 벤치마크 실행으로 3분자 복원 확인 | fixed
SCI-09-035 | 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md:63, 72, 82, 255-278 | 04장: 누적 오류 0.17~1.67%로 정정, 20~40%·150% 팽창을 조건부 예시로 표기(Smith 2017 1.3배) | fixed
SCI-09-036 | 04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md:7, 16 | 04장: 총 RNA 10~30 pg(mRNA는 일부)로 통일, PCR 사이클은 프로토콜 의존 예시 | fixed
SCI-09-039 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:351-363, 241 | 05장: xf 비트필드(Cell Ranger 전용) 처리, 비트 1 검사, head 파이프의 pipefail 대응 | fixed
SCI-09-040 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:9, 22, 146-163, 405-411 | 05장: FLAG 1107 연습을 일반 PE BAM 가정 예시로 표기(10x 3′ BAM은 R2 단일 말단), "수억 개 150bp PE" 문구 완화 | softened
SCI-09-041 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:456-459 | 05장 Q: SA 태그는 주 정렬에 남음, 100% 위음성 → 민감도 저하, STAR chimSegmentMin·WithinBAM 필요 명시 | fixed
SCI-09-042 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:462-464 | 05장: AR-V7을 엑손 4-8 건너뛰기가 아닌 CE3 말단 엑손 교체로 서술 | fixed
SCI-09-043 | 05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md:305-306, 321 | 05·06장: AR 구간 chrX:67,543,000-67,731,000 (Ensembl ENSG00000169083 67,543,353-67,730,619 및 RefSeq 67,544,020-67,730,619를 모두 덮도록 재보정; 1차 67,544,000은 5′ 약 650 bp 누락, 스크립트 AR_REGION 변수 누락분도 정정) | fixed
SCI-09-045 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:3, 161-164 | 06장: 엔잘루타마이드 작동제 전환은 F876L(F877L), T878A는 이중변이에서 증폭, L702H는 엑손 4 | fixed
SCI-09-046 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:277-285, 377-396 | 06장 코드: min_total_umi 임계 구현, Hemizygous/Mixed 용어, ADAR 편집 주석, Patient C 판정 불가로 수정 | fixed
SCI-09-047 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:39-41, 362 | 06장: fetch는 bam_dup1 독립 복사본을 반환(버퍼 재사용 아님), list 메모리 소모는 객체 보관 때문 | fixed
SCI-09-048 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:139-148 | 06장: pileup min_base_quality 기본 13, max_depth 8000(조용히 잘림) 안내, max_depth 상향; S2 변환 중 깨진 start/stop(67723877/67723710)도 정정 | fixed
SCI-09-049 | 06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md:317 (Q1) | 06장 Q1: FOXA1 좌표를 산술 연습용 임의 좌표로 명시(실제 codon 253은 chr14:37,592,025-27) | softened
SCI-09-053 | 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:146, 379 | 07장: 6개월 전 예측 문구 삭제, 탁산 다이닌 기전은 가설·재현 연구 병기 | softened
SCI-09-054 | 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:21, 41-42 vs 90-94, 249-259 | 07장: LBD 엑손 표기 통일 주석, 경험적 베이즈 적합 없는 고정 사전분포 명시, 임상 라벨→정션 라벨 | fixed
SCI-09-055 | 07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:25-26, 44-50 | 07장: Cryptic Exon 3b·엑손 4-8 건너뛰기 표기 전면 정정(CE3, 말단 엑손), 8.5 kb 유지 | fixed
SCI-09-058 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:80-81, 93-94, 317-320 | 08장: FLAG 99 = 1+2+32+64 분해 정정, SA 예시 일치, 도식 개념용 주석 | fixed
SCI-09-062 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:367-369, 381-382 | S2(SCI-09-061/063 포함 08장 Q3)에서 8p 공유 단클론/다초점 판정 불가로 정정 완료; 3′ 검출률 수치는 가상 표기 추가 | already-fixed
SCI-09-064 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:11, 16-24, 385 | S2(SCI-09-056)에서 ERG/TMPRSS2 방향·위치 도식 정정, DSCR3 삭제 완료 | already-fixed
SCI-09-065 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:49-59 | 08장: Type 번호·빈도(~65/25/5/5%)를 출처 미확인 예시로 표기, T1-E4는 N-말단 약 32 aa 결손 near-full-length ERG | softened
SCI-09-066 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:117-137, 211-286 | 08장: 미세 상동성 기준 ≥5로 통일, 100% 신뢰도·수치 삭제, 코드 구현 한계·STAR chim 옵션 명시 | softened
SCI-09-067 | 08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md:3, 9, 37 | 08장: 서구 원발 40~50%, 결실 60%/전좌 40%의 분모를 융합 양성으로 명시 | fixed
SCI-09-071 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:128-139 | 09장: Unique=다유전자 UMI 제외, Rescue 정의 교정(PropUnique 구분) | fixed
SCI-09-072 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:110 | 09장: "이중 무릎" 삭제, EmptyDrops_CR을 CellRanger flavor EmptyDrops(파라미터)로 교체 | fixed
SCI-09-073 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:43-45 | 09장: genomeSAindexNbases는 SAindex k-mer 길이, 희소도는 genomeSAsparseD로 정정 | fixed
SCI-09-074 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:156, 158 | 09장 표: 빌드/접두사 구분, ambient 처방을 SoupX/CellBender/DecontX로 교체 | fixed
SCI-09-075 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:214-223 | 09장: 44.4%를 GeneFull 추가 포획 분자 비율로 재해석, 핵 내 전스플라이싱 비율 단정 삭제 | fixed
SCI-09-076 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:187-198 | S2에서 Stage 4/5 "100%" 기준을 확인·대리지표로 교체 완료; 이번에 "완벽" 표현 완화 | already-fixed
SCI-09-077 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:97 vs Track 00 SCI-00-023/024; 06:161 vs  | 교차 트랙 반복 오류 4건은 S1/S2(SCI-09-068, 044/045, 060, 이번 045)로 해소됨을 확인 | already-fixed

비고: S2 편집 중 06장 pileup 예시의 start/stop이 `67723877/67723710`로 깨져 있던 것을 발견해 정정(SCI-09-048 행). 인트론 40~60% 등 트랙 간 용어는 일괄 치환하지 않음.

## S3 단계 / Track 04

(32건. 수정은 `tracks/04_cancer_genomics_clones/docs/*.md`만; 제어문자 0, `$$` 짝 OK, python ast OK. 06장 베타-이항 코드는 합성 데이터로 실행 확인.)

ID | 파일:줄 | 조치 | 상태
--- | --- | --- | ---
SCI-04-004 | 00_CANCER_GENOMICS_ROADMAP.md:63-71 | 00장: 3번 중립/빅뱅 모델은 Nowell이 아닌 후속 모형(Sottoriva 2015)으로 귀속 정정 | fixed
SCI-04-005 | 00_CANCER_GENOMICS_ROADMAP.md:160-194, 204-250 | 00장 스크립트: conda.sh source 추가, R 도구/cellsnp-lite/Eagle2 설치 안내, read_csv names 4열, SettingWithCopy 방지 | fixed
SCI-04-008 | 01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md:79, 87, 98 | 01장: NKX3-1/PTEN/AR 수치를 GISTIC 기준·데이터셋 병기(TCGA/SU2C), AR 인핸서 중복은 Viswanathan 70~87%로 분리 | fixed
SCI-04-010 | 01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md:142-148, 199-206, 205 | 01장: 4.1 공식의 이상 가정·유전자 용량 보상 한계, 남성 X 1카피·infercnvpy chrX 기본 제외 주석 + 코드 주석 | fixed
SCI-04-011 | 01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md:79 | 01장: Nkx3.1 haploinsufficiency는 마우스 근거로 한정, i(8q) "강력한 지표" 완화 | fixed
SCI-04-012 | 01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md:98, 280 | 01장: 인핸서 증폭을 Viswanathan/Takeda 2018로 귀속(SU2C 2015 아님), "수십 배" 삭제, 엔잘루타마이드 부분 작동제는 F877L과 구분 | fixed
SCI-04-016 | 02_INFERCNV_ALGORITHM_AND_HMM.md:296-297 | 02장 Q1: "6번 염색체 획득" 방향 서술을 참조 대비 방향 혼재로 정정(해설 방향은 S2에서 정정) | fixed
SCI-04-019 | 03_REFERENCE_FREE_COPYKAT.md:184-190, 320-334 | 03장: CopyKAT K-S 실제 구현(MCMC 사후평균 분포 적용, KS.cut 0.1, win.size 25) 주석, Q3 n≤2 ECDF 논거를 일반 이론으로 한정 | fixed
SCI-04-021 | 03_REFERENCE_FREE_COPYKAT.md:43-50 | 03장 §1.2: 영역별 기저선 산술 정정(8p는 정상 세포 gain 역전), 8q만 기존 예 | fixed
SCI-04-022 | 03_REFERENCE_FREE_COPYKAT.md:106, 313 | 03장: 분산 기하급수→평균의 이차, Anscombe/FT는 Poisson 가정 분산 안정화(NB에서 일정 분산 비보장), 90% 수치 삭제 | fixed
SCI-04-027 | 04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md:57 | 04장: 이형접합 SNP 간격 100~300 bp → 대략 1~1.5 kb당 1개(개인 차 큼) | fixed
SCI-04-028 | 04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md:306-321 | 04장: SHAPEIT4→Eagle2, 위상 스위치는 HMM 하플로타입 전이(θ)로 모형화, 베타-이항은 과산포 흡수 | fixed
SCI-04-030 | 05_SUBCLONAL_DECONVOLUTION_AND_PHYLOGENY.md:98-100 | 05장: NJ 복잡도 표준 O(n³), O(n²)는 RapidNJ/FastME류, 반복 n−3회 | fixed
SCI-04-034 | 05_SUBCLONAL_DECONVOLUTION_AND_PHYLOGENY.md:185, 94 | 05장: Enzalutamide는 ADT가 아닌 ARSI로 용어 정정 | fixed
SCI-04-035 | 05_SUBCLONAL_DECONVOLUTION_AND_PHYLOGENY.md:394-395 | 05장 Q: Alisertib(N-Myc)·Tazemetostat(PRC2) 기전 분리, NEPC 근거는 초기/전임상으로 가설 표기 | softened
SCI-04-037 | 06_SOMATIC_MUTATIONS_AND_FUSIONS.md:77 | 06장: 변이–약제 대응 정정(F877L 작동제 전환, T878A 프로게스테론 등, L702H 엑손 4 글루코코르티코이드) | fixed
SCI-04-038 | 06_SOMATIC_MUTATIONS_AND_FUSIONS.md:256-262, 277-290 | 06장 코드: 주석 beta(1,1/ε) 정정, min_total_umi(기본 2) 미만은 UNCERTAIN, 사전확률 인자화; 합성 데이터로 실행 확인 | fixed
SCI-04-039 | 06_SOMATIC_MUTATIONS_AND_FUSIONS.md:334-336 | 06장: P(k=0) 0.164→0.246, 이항 대비 약 8배 | fixed
SCI-04-040 | 06_SOMATIC_MUTATIONS_AND_FUSIONS.md:15-16, 96-98 | 06장: MAPQ 60/255, BQ·MQ·4bp·UMI≥2 임계값을 예시로 표기 | softened
SCI-04-043 | 07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:98-106, 111, 113-115 | 07장: 4-state 표를 경향적 대응으로 표기, DNPC의 FGF/MAPK(Bluemn)와 Wnt(Tang CRPC-WNT) 분리, t-NEPC 양측결손 필수 아님 | softened
SCI-04-044 | 07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:27-39 vs 262-270 | 07장: 클론 ODE를 증식항에만 로지스틱 제한(−d_i N_i 분리), Q2 지수 근사 관계 명시 | fixed
SCI-04-045 | 07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:166-179 | 07장 코드: DNPC 명시 조건(AR≤0, NE≤0.1) 추가, 경계 영역은 Unclassified | fixed
SCI-04-046 | 07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:315-317 | 07장 Q: cisplatin+docetaxel 삭제→carboplatin/cisplatin+etoposide(가이드라인 확인 필요), Tarlatamab/Tazemetostat은 시험적 선택지 | softened
SCI-04-050 | 08_METASTASIS_SEEDING_MODELS.md:25-35, 21 | 08장: 군집 효율 20~50배를 Aceto 2014(마우스 유방암)로 귀속, 가상 모형의 배율(490~8,600배)·P→1 한계 명시 | fixed
SCI-04-051 | 08_METASTASIS_SEEDING_MODELS.md:98-119 | 08장: 골전이 공간 분업 서술을 가설로 표기, 조골성/혼합 병변과의 불일치 명시 | softened
SCI-04-052 | 08_METASTASIS_SEEDING_MODELS.md:240-257 | 08장: t_split→0 극한을 이상화로 한정(founder·생식세포 변이 무시), 선형·평행 모델은 같은 식의 t_split 차이 | fixed
SCI-04-055 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:117 vs 02:9 | 09장: CN state 0~6 → inferCNV i6 6상태 {0,0.5,1,1.5,2,3} | fixed
SCI-04-056 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:237-260 vs 02:64-66 | 09장 Q2: 영역 특이적 참조 편향이 실제 위험이고 세포별 median 중심화가 전역 오프셋을 일부 제거함을 전제에 명시 | softened
SCI-04-057 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:66-67, 176-185 | S2(SCI-04-060)에서 GMM 센터링을 휴리스틱(WGD·고배수성 주의)으로 정정 완료 확인 | already-fixed
SCI-04-058 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:287 | 09장: WES VAF와 "정확히 일치" → 대략적 교차 검증(2×VAF/순도, ASE·NMD 영향) | fixed
SCI-04-059 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:114-121 | 09장: 루브릭 배점·80점 합격선은 임의 설계값, 재현성·직교 검증 정성 기준 병행 권고, 장가지 유인 용어 정정 | softened
SCI-04-061 | 07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:3, 93-115 (vs Track 00 SCI-00-015 | 07장: Track 00/01 정의 귀속 삭제 → Labrecque 2019/Bluemn 2017 참고 AR/NE 축 분류, Tang 2022 4군과 구분 | fixed

비고: Track 09 보정 — 05장 AR 구간을 chrX:67,543,000-67,731,000으로 재보정(SCI-09-043 행 및 s3 CSV 비고 갱신, 스크립트 AR_REGION 변수 포함).

## S3 단계 / Track 07

(31건. 수정은 `tracks/07_count_statistics_abundance/docs/*.md`만; 제어문자 0, `$$` 짝 OK, python ast OK. 사용 수치: F(1,4)=7.71·F(1,24)=4.26·F(1,54)=4.02, t(df=4) p=1e-50 -> t≈4.9e12, power(d=0.05,n=10,000/군)=0.942를 scipy로 재계산, propeller 열 이름은 speckle 소스로 확인.)

ID | 파일:줄 | 조치 | 상태
SCI-07-004 | 00_COUNT_STATISTICS_ROADMAP.md:160-181 | 00장 로드맵 장 구성을 실제 파일(00~09) 기준으로 재작성 | fixed
SCI-07-005 | 00_COUNT_STATISTICS_ROADMAP.md:187-233 | sccomp를 Bioconductor로 이동, tradeSeq/slingshot/speckle/IHW/RUVSeq/sva 추가; Python 3.11+decoupler>=2.0로 03장 dc.pp.pseudobulk와 정렬 | fixed
SCI-07-006 | 00_COUNT_STATISTICS_ROADMAP.md:80 | df=4 p<1e-50 서술을 "t≈4.9e12 필요, 현실 데이터로 도달 불가"로 정정 | fixed
SCI-07-007 | 00_COUNT_STATISTICS_ROADMAP.md:313-315 | pseudobulk "완벽 무력화"를 "크게 완화"로 수정 | fixed
SCI-07-010 | 01_PROBABILITY_MODELS_POISSON_NB.md:165-170; 03:76; 06:355 | 드롭아웃 표현을 샘플링 영(zeros)으로 통일(01 도식, 03, 06) | fixed
SCI-07-011 | 01_PROBABILITY_MODELS_POISSON_NB.md:371-378 | R LRT를 2*(logLik 차)와 0.5*pchisq 경계 보정으로 교체 | fixed
SCI-07-012 | 01_PROBABILITY_MODELS_POISSON_NB.md:430-431 | 경계에서 우도 평탄/식별 어려움, LRT는 혼합 χ²로 서술 정정 | fixed
SCI-07-013 | 01_PROBABILITY_MODELS_POISSON_NB.md:471-482 | "UMI 평균이 10배가 된다고 가정"으로 가정 명시(산술은 재계산 일치) | fixed
SCI-07-019 | 02_PSEUDOREPLICATION_PITFALL.md:18, 21 | Hurlbert 인용 블록을 의역 표기, "약 절반"으로 완화 | softened
SCI-07-020 | 02_PSEUDOREPLICATION_PITFALL.md:324 | Cohen d=0.05 검정력 99% 이상 -> 약 94%(scipy 0.942) | fixed
SCI-07-021 | 02_PSEUDOREPLICATION_PITFALL.md:349-353 | 850개·20~30개를 설명용 가상 수치로, 실제 GSE143791에 원발 대조 없음 명시 | softened
SCI-07-027 | 03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:354, 371-372 | F(1,4)=7.71은 d0=0일 때만; d0=20 -> F(1,24)=4.26, d0=50 -> F(1,54)=4.02(scipy)로 EB 설명 정정, 1,200 vs 45는 가상 표기, "완벽 통제" 삭제 | fixed
SCI-07-028 | 03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:399-405 | 연령 랭크 결손을 환자 내 불변 공변량(환자당 한 값 A_i)으로 한정, 시점 연령은 "거의" 종속 명시 | softened
SCI-07-029 | 03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:343-347 | NB IRLS 가중치 μ/(1+αμ)의 1/α 포화를 명시, 세포 수에 비례한 가중이 아님을 추가 | softened
SCI-07-030 | 03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:10, 192; 00:11, 189-196 | Python 패키지명 decoupler(R은 decoupleR) 표기 정정, 00장 설치 블록 decoupler>=2.0·Python 3.11로 정렬(03장 코드는 ENG-011로 이미 dc.pp.pseudobulk) | fixed
SCI-07-039 | 05_FAST_CELL_LEVEL_NEBULA.md:211 | FOXA2를 "신경내분비 분화와 연관된 전사인자"로 정정(AR 스플라이스 변이 조절 인자 삭제) | fixed
SCI-07-040 | 05_FAST_CELL_LEVEL_NEBULA.md:410-412, 469, 473-474 | ω 차감식 삭제(재적합 분산 성분 추정으로 서술), 과소추정 배율을 √VIF로 정정(√(N/J)=70.7은 ρ→1 상한, ρ=0.05 예시 약 16배) | fixed
SCI-07-046 | 06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md:198-199 | 조기 분기 t*를 상위 구동 TF로 보는 해석을 가설로 한정, 의사시간의 비인과성 명시 | softened
SCI-07-047 | 06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md:44-48 | DPT를 확산 거리 기반 상대 순서(측지선 아님), 잠재 시간을 모형 의존적 상대 시간(절대 시간 아님)으로 정정 | fixed
SCI-07-052 | 07_DIFFERENTIAL_ABUNDANCE_MILO.md:117-122, 409 | log(N_s) 오프셋 "완벽 정규화" 삭제, miloR norm.method 기본 TMM(RLE/logMS 선택)·조성 효과 잔존 명시 | fixed
SCI-07-053 | 07_DIFFERENTIAL_ABUNDANCE_MILO.md:296-325 | docstring을 유의 이웃 logFC 단순 평균(np.mean), 유의 이웃 없으면 0으로 코드와 일치시킴 | fixed
SCI-07-054 | 07_DIFFERENTIAL_ABUNDANCE_MILO.md:408-410 | "p<1e-15" 삭제, 오프셋 누락 시 유의한 전역 이동으로 서술하고 샘플 ~10개에서는 극단적 p 불가 명시 | softened
SCI-07-057 | 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:108 | ilr(p)=clr(p)Ψ^T -> clr(p)Ψ (1×D · D×(D-1)) | fixed
SCI-07-058 | 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:289-296 | propeller 결과 열(BaselineProp.*, PropMean.*, PropRatio, Tstatistic, P.Value, FDR; Log2FD 없음)·기본 transform logit을 speckle 소스로 확인해 정정 | fixed
SCI-07-059 | 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:309-327 | 조성 방법은 기준 대비 상대 변화만 추정, 절대 개체수는 IHC·비드 보정 유세포 등 외부 측정 필요 문장 추가 | softened
SCI-07-060 | 08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md:336-379 | Q2 전제를 "절대 수 불변(상대 비율은 상승)"으로 정정, β_c는 기준 대비 상대 효과임을 명시, 정규 사전 "거의 모든 CI" -> "일부" | fixed
SCI-07-064 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:107-117 | IHW 제약식을 층 크기 가중 (1/G)ΣG_b w_b=1 (Σ_g w_g=G)로 수정, "엄격 보장" 완화 | fixed
SCI-07-065 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:259-272, 340-341 | DESeq2 기본 independent filtering 때문에 IHW 이득이 데이터 의존적임을 코드 주석에 추가, FDR 50% 이상 삭제(설명용), 독립성 "완벽 보장" -> 근사적 성립 | softened
SCI-07-066 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:130 | "샘플링 드롭아웃" -> 샘플링 영(zeros), ZINB 오용의 "과도한 위음성" 삭제 | fixed
SCI-07-067 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:146-177 | N>=15 등 루브릭 임계값을 교육용 예시로 명시(검정력 분석으로 대체), "Spatial FDR 완벽 적용" 삭제 | softened
SCI-07-069 | 00_COUNT_STATISTICS_ROADMAP.md:9, 91-121 | S2에서 00장이 Kfoury 9+7, Chen 13개(36,424 세포)로 정정됨; GSE141445 GEO 요약(13 tumors, 36,424 cells)으로 재확인 | already-fixed

비고: F877L 표기 통일 — 정식 표기를 F877L(NM_000044.6 번호, T878A와 동일 체계)로 하고 "구 표기 F876L"은 장별 첫 언급에만 병기. F876L 단독 표기는 남아 있지 않으며 병기는 4곳(Track 04 01장:287, Track 04 06장:79, Track 09 06장:163, docs/01_FOUNDATIONS_AND_WET_LAB.md:41)뿐. T877A 표기는 어디에도 없음.
비고: SCI-07-030 — 00장 설치 블록이 Python 3.10 + decoupleR==1.8.0이었으나 03장 코드는 ENG-011로 2.x API(dc.pp.pseudobulk)를 쓰므로 Python 3.11 + decoupler>=2.0으로 정렬(decoupler 2.2.0 Requires-Python>=3.11 확인; 다른 핀(scanpy 1.10.1 등)과의 동시 설치는 미검증).

## S3 단계 / Track 01

(30건. 수정은 `tracks/01_multiome_snatac/docs/*.md`만; 제어문자 0, `$$` 짝 OK, python ast OK(04장 해설 코드는 목록 내 들여쓰기라 dedent 후 확인). 확인 근거: macs3 `--help`(FRAG, keep-dup, hs 2,913,022,398), Signac pbmc_vignette(임계값), muon `l2norm`/`neighbors` docstring, snapatac2 2.10.0 시그니처(frip/frag_size_distr/filter_cells), PubMed 9516433 초록.)

ID | 파일:줄 | 조치 | 상태
SCI-01-004 | 00_MULTIOME_ROADMAP.md:30-33 | "어떠한 가상 매칭 오류 없이 ground truth" -> 동일 GEM 동시 포획으로 사후 매칭 오류 없음, doublet·앰비언트 오염·바코드 변환 QC 필요 명시 | softened
SCI-01-006 | 01_EPIGENOMICS_AND_TN5_CHEMISTRY.md:74-76, 129-139 | +4/-5 시프트는 BAM 직접 사용 시에만, fragments.tsv.gz/Chromap 출력은 보정 완료(이중 시프트 방지)를 01장 본문·Q1·해설에 명시(ArchR offsetPlus/offsetMinus는 BAM 전용) | fixed
SCI-01-007 | 01_EPIGENOMICS_AND_TN5_CHEMISTRY.md:44-46 | E54K·M56A·L372P를 Goryshin & Reznikoff 1998 귀속(원문 목록 미확인 명시), "1,000배"·자가억제 도메인 기전 삭제 | softened
SCI-01-008 | 01_EPIGENOMICS_AND_TN5_CHEMISTRY.md:109-114 | 핵막 "구멍 형성" 기전·"배경 0 차단" 삭제, 투과화된 핵 내부로 Tn5 접근, 세척으로 배경 감소로 완화 | softened
SCI-01-011 | 02_UPSTREAM_AND_FRAGMENTS.md:59-60 | Chromap: 접두사 트리 -> 미니마이저 기반 근사 매핑, "완벽 보존" -> BWA 대비 10~20배 빠르고 유사한 정확도(README) | fixed
SCI-01-012 | 02_UPSTREAM_AND_FRAGMENTS.md:111-115, 188-189 | Tabix O(1)·B-tree 삭제 -> 계층 bin + 16 kb 선형 인덱스, 해당 BGZF 블록만 읽음(비용은 겹치는 bin·레코드 수에 비례); 코드 주석 포함 | fixed
SCI-01-013 | 02_UPSTREAM_AND_FRAGMENTS.md:102-107 | chromStart/chromEnd를 모두 0-based half-open(BED)으로 통일 | fixed
SCI-01-015 | 03_SNATAC_QC_METRICS.md:30-33, 170 | TSS 정의를 ArchR 기준(중심 101 bp / 플랭크 +-1,900~2,000 bp)으로 03장 식·Q1·본문 정정, 도구별 차이 명시(09장과 일치) | fixed
SCI-01-016 | 03_SNATAC_QC_METRICS.md:45-51, 126, 221; 09:63-68, 138-140 | TSS 등급을 09장 표(<4/4~6/6~10/>=10) 하나로 통일하고 ENCODE/10x 귀속 삭제(ArchR minTSS=4·Signac >4 참고 경험적 가이드), 코드 max_counts 100,000->80,000(Q3와 일치), 09장 루브릭 TSS 기준 정렬 | fixed
SCI-01-017 | 03_SNATAC_QC_METRICS.md:92-101 | mtDNA "히스톤 없음·나체 DNA" -> TFAM 핵양체, 핵 뉴클레오솜 미포장; 30~60% 삭제, chrM >5% 경고/>20% 심각(09장과 통일, 프로토콜 의존 경험적 가이드) | fixed
SCI-01-018 | 03_SNATAC_QC_METRICS.md:전체(2.3 추가) | 03장에 2.3 신설: nucleosome signal(147-294/<147 bp), FRiP, blacklist ratio 정의·Signac 예시 코드·임계(<4, >40, <0.01) 및 snapatac2 frag_size_distr/frip 안내(시그니처 확인) | fixed
SCI-01-019 | 04_PEAK_CALLING_AND_LSI.md:19-21, 79-84, 189 | "반드시 이진화"·"신호는 0,1,2뿐/4는 PCR 편향" 삭제 -> 이진화는 선택(ArchR 기본 binarize=TRUE, Signac RunTFIDF는 카운트 입력), 2배체 상한은 휴리스틱 | softened
SCI-01-020 | 04_PEAK_CALLING_AND_LSI.md:41-57 | MACS3 명령을 -f FRAG --barcodes로 교체(keep-dup 항상 all; BEDPE는 --keep-dup all 명시 필요), hs 유효 게놈 2,913,022,398(MACS2 2.7e9와 구분), ArchR cut-site 방식 병기(macs3 --help로 확인) | fixed
SCI-01-021 | 04_PEAK_CALLING_AND_LSI.md:72-76 | iterative overlap에 ArchR의 점수 정규화·pseudo-replicate 재현성(기본 2)·501 bp·q 0.1 병기 | fixed
SCI-01-022 | 04_PEAK_CALLING_AND_LSI.md:112-114 | IDF->0 수렴 삭제, 하한 ln2(Q2 해설과 일치), 이 식=Cusanovich 계열(Signac method 2)·Signac 기본 method 1·muon 기본 log(TF)xlog(IDF) 명시 | fixed
SCI-01-023 | 04_PEAK_CALLING_AND_LSI.md:160, 245; 09:184, 297 | LSI1-깊이 상관 임계를 ArchR corCutOff 0.75로 통일(04 Q3 본문·코드, 09 코드 0.85, 09 루브릭 r>0.9); 시나리오의 관측값 r>0.9는 유지 | fixed
SCI-01-027 | 05_WNN_MULTIMODAL_INTEGRATION.md:135-139, 143-148 | LSI 1 제외 성분 2~31을 별도 obsm(X_lsi_2plus)에 두고 mu.pp.l2norm 적용 후 neighbors; RNA도 l2norm 추가; muon lsi scale_embeddings=True 주석(muon 소스로 l2norm·neighbors 문서 확인) | fixed
SCI-01-028 | 05_WNN_MULTIMODAL_INTEGRATION.md:31, 110-112; 00:31 | "수 시간~수일" 삭제, 염색질 선행 경향은 Ma 2020(chromatin potential) 인용 + 시간 척도 맥락 의존 명시(00:31 포함) | softened
SCI-01-032 | 06_PEAK_TO_GENE_AND_ENHANCERS.md:79-85, 94-114, 188-198 | 윈도우를 도구 기본값(Signac 500 kb, ArchR 250 kb, Cicero 250 kb)으로 명시하고 AR 650 kb 포착엔 확대 필요 안내, 배경 100개->Signac 기본 200개, FDR<0.05·r>0.30은 예시 표기, 코드 window 1 Mb·multipletests import·헬퍼 함수 자리표시자 명시 | fixed
SCI-01-033 | 06_PEAK_TO_GENE_AND_ENHANCERS.md:38-40 | LLPS 전사 응축체를 논쟁 중인 가설로 표기, 루프 압출·CTCF 서술은 유지 | softened
SCI-01-035 | 07_TF_FOOTPRINTING_AND_CHROMVAR.md:37-39 | 풋프린트 깊이 "점유율에 정비례" -> 정성적 지표(TF 체류시간·서열 편향·깊이 의존, Bentsen 2020), 플랭킹 "구조 안정화" 설명 삭제 | softened
SCI-01-036 | 07_TF_FOOTPRINTING_AND_CHROMVAR.md:105-117 | chromVAR 편차를 Y=(O-E)/E, 편향 보정 편차=Y-mean(Y_bg), z=(Y-mean)/sd로 정정(Schep 2017 Methods) | fixed
SCI-01-038 | 07_TF_FOOTPRINTING_AND_CHROMVAR.md:93-96 | AR "FOXA1이 열어둔 인핸서에만" -> 맥락 의존적, HOXB13 파이오니어 분류는 논쟁적, ASCL1 마스터 스위치/선암 인핸서 폐쇄 -> ARPI 후 염색질 재편 보고(Nouruzi 2022) | softened
SCI-01-042 | 08_EGRN_AND_SCENICPLUS.md:51, 89, 113-120 | SCENIC+ get_search_space 기본(1~150 kb) 명시, AR 650 kb는 검색 공간 확장(upstream=[1000,700000]) 필요를 2.3절·표·코드 주석에 추가 | fixed
SCI-01-043 | 08_EGRN_AND_SCENICPLUS.md:189-199, 표 | Erdafitinib/Trametinib 병용을 전임상 수준 가설로 표기(임상 입증 아님); ASCL1/ETV4 표 서술도 연관 가설로 완화 | softened
SCI-01-045 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:93-94, 296 | iLISI 범위 1~B(배치 수), 목표 2.0은 B=2일 때만, 정규화 iLISI(0~1) 안내, 환자 12명 케이스의 1.8은 거의 미보정 상태임을 명시 | fixed
SCI-01-046 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:63-68, 41, 121-143, 217-218 | 09장 TSS 등급·chrM(>5% 경고/>20% 심각)·MACS2->MACS3 통일, Stage 4 순열 서술을 06장 기준과 정렬, WNN "정보-이론적 예측 오차" -> 커널 기반 예측 비교 | fixed
SCI-01-047 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:172-175 | muon tfidf 기본 log(TF)xlog(IDF)(Signac 기본과 다름)·PCA(zero_center=False)는 성분 표준화 없음을 코드 주석에 정정 | fixed
SCI-01-048 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:98-109, 246-259 | EDTA 탈회를 가정적 프로토콜로 표기, Multiome 검증 출처 미제시·RNA 무결성 위험·잔류 EDTA의 Mg2+ 의존 Tn5 억제 주의, "원천 봉쇄/불가능" 완화 | softened
SCI-01-049 | 00_MULTIOME_ROADMAP.md:124-128, 183-198; 05:84-88; 06:143, 184; 08 | S2에서 00장 GSE210358 등 데이터셋 정정, 05장 WNN 친화도 식 정정, 06장 AR 인핸서 복제 70-87%(Viswanathan 2018) 반영 확인; S3에서 08장 DNPC/ASCL1 마스터 레귤론 표현 완화 | already-fixed

ENG-010 / ENG-011 (applied by sci): S2 단계에서 적용 완료(import_data->import_fragments, dc.get_pseudobulk->dc.pp.pseudobulk); S3에서 변경 없음.
비고: 03장 코드 임계값 변경(max_counts 80,000)은 API 인자명 변경 없이 값만 수정. 03·09장 QC 임계값은 도구 기본값 기반 경험적 가이드로 표기(ENCODE/10x 공식 기준으로 귀속하지 않음). 이로써 Track 07까지 포함한 sci S3 전 트랙 완료.

## S3 단계 / Track 06 (재배정)

(25건. 수정은 `tracks/06_manifold_trajectory/docs/*.md`만, 새 절 추가 없음(기존 문장 수정·짧은 주석·표 아래 주석 1줄). 제어문자 0, `$$` 짝 OK, python ast OK. 확인 근거: scvelo 0.3.4 velocity_confidence 소스, PubMed 34312180 초록, PubMed esummary(31361600·29017058·26855148), pipeline_01_summary.txt. SVG·생성 스크립트·results는 수정하지 않음.)

ID | 파일:줄 | 조치 | 상태
SCI-06-005 | 06_MARKOV_CHAINS_AND_MACROSTATES.md:58 | S2에서 "Reuter et al., J. Chem. Theory Comput. 2018 (PMID 29812922); J. Chem. Phys. 2019 (PMID 31067901)"로 정정됨, 06장 L58 재확인 | already-fixed
SCI-06-006 | 05_CELLRANK2_MULTIKERNEL_INTEGRATION.md:49 | S2에서 CellRank 2를 Nature Methods 2024 (PMID 38871986)로 정정됨, 05장 L49 재확인 | already-fixed
SCI-06-008 | 01_MANIFOLD_GEOMETRY_AND_DIFFUSION.md:259 | 코드 주석으로 method=umap(기본)은 본문 가우시안 커널이 아닌 UMAP 연결성 가중임을 명시, 가우시안에 가깝게 하려면 method=gauss | fixed
SCI-06-009 | 01_MANIFOLD_GEOMETRY_AND_DIFFUSION.md:274 | 주석 정정: n_branchings 기본 0, 1은 분기 1회 검출(실험적), 다중 분기는 PAGA 권고 | fixed
SCI-06-014 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:35-36 | 분열 주기 빠른 세포의 γ 증가를 가설로 한정(직접 연관 근거 미제시) | softened
SCI-06-016 | 04_OPTIMAL_TRANSPORT_WADDINGTON_OT.md:139-143 | Schiebinger 성장률을 "증식·자멸 점수를 로지스틱 출생·사멸률에 대응", 본문 식은 선형 점수 지수식의 단순화 교육용 예시(원 논문 식 아님)로 표기 | fixed
SCI-06-023 | 07_FATE_ABSORPTION_PROBABILITIES.md:226-230 | S2에서 circular_projection이 keys="cell_type"(obs 키)로 정정되고 운명 확률은 사전 계산(lineages_fwd)으로 서술됨, 07장 L226-232 재확인 | already-fixed
SCI-06-024 | 07_FATE_ABSORPTION_PROBABILITIES.md:178-181 | Clone A/B/C와 0.90/0.85/0.92를 설명용 가상 예시로 표기(Track 04 정의 아님), AR+ 선암종이 초기/말단 양쪽 의미를 가질 수 있음 주석 | softened
SCI-06-027 | 08_LINEAGE_DRIVER_DISCOVERY.md:62,229 | "완벽히 차단" -> 줄인다, 영과잉 시 ZINB 필요 가능성, 흡수 확률이 같은 발현 kNN에서 파생되는 순환성(반보수적 p)·다중검정 한계와 순열·pseudobulk 보완 언급, link 설명도 완화 | softened
SCI-06-028 | 08_LINEAGE_DRIVER_DISCOVERY.md:46-47 | Spearman 식은 동순위 없을 때만 유효, 실제는 평균 순위 Pearson 상관임을 식 직후에 명시 | fixed
SCI-06-029 | 08_LINEAGE_DRIVER_DISCOVERY.md:99-110,245-251 | θ_slope는 임의 조정 모수, 활성화 순서는 인과 증명 아님 -> 후보 우선순위 도구로 한정하고 CRISPR·염색질 선행 검증 명시; 치료 표적 "반드시" -> 우선 검토 후 실험 검증 | softened
SCI-06-030 | 08_LINEAGE_DRIVER_DISCOVERY.md:125; 09:45 | REST 복합체 분해 -> REST 기능 상실(SRRM3/4에 의한 REST4 스플라이싱, 발현 감소; Labrecque 2021 Cancer Res PMID 34312180 초록으로 확인), 09장 유비퀴틴화 분해 서술도 정정 | fixed
SCI-06-031 | 08_LINEAGE_DRIVER_DISCOVERY.md:106-107,126 | 파이오니어: ASCL1·NEUROD1은 보고됨, INSM1은 NE 프로그램 TF로 근거 불분명 명시, "핵심 인과적 구동자" -> 구동자 후보 | softened
SCI-06-032 | 08_LINEAGE_DRIVER_DISCOVERY.md:133-139 | SOX9/TWIST1을 EMT 후보(가상 예시), ETV4/5 "완전 대체·독점" -> 대신하는 후보 시나리오로 완화하고 Bluemn 2017 Cancer Cell (PMID 29017058) 인용 | softened
SCI-06-035 | 08_LINEAGE_DRIVER_DISCOVERY.md:268 | ASCL1 모티프를 "CANNTG(모든 bHLH 공통 E-box) 중 CAGSTG 계열 선호"로 정정 | fixed
SCI-06-036 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:66-68,138-140 | velocity_graph 선행 실행 추가, 정의를 scVelo 구현(유전자 평균 제거 후 이웃과의 Pearson 평균, 0 미만 절삭; scvelo 0.3.4 소스 확인)으로 주석하고 0.75/0.50은 경험적 예시로 표기 | fixed
SCI-06-037 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:57-83 | 2.2 그래프 일관성·2.3 CBC를 두 패키지에 동명 함수 없는 교육용 자체 정의 지표로 명시, scVelo velocity_confidence·CellRank Kernel.cbc 언급 | softened
SCI-06-041 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:122,193,259 | Unspliced 비율 기준을 샘플 유형별 예시(scRNA 15-25%, snRNA 40-60%)로 루브릭·Q1 해설·원인 진단에서 통일 | fixed
SCI-06-042 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:110-112 | "원시 카운트" -> 배치 보정하지 않은 spliced/unspliced 레이어(정규화·kNN 평활 모멘트 Ms, Mu 이후) | fixed
SCI-06-043 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:122-126 | 루브릭 합격 기준(0.75, 0.90, 1e-5 등)을 출처 없는 예시 임계로 표 아래 주석 추가, 파이프라인 01 평균 일관성 0.161(결과 파일 확인)과 불일치 명시. L126 표 파이프 파손은 ENG-022 소관 | softened
SCI-06-044 | 07_FATE_ABSORPTION_PROBABILITIES.md:157; 08:116; 09:55 | 4상태(AR+/Amphicrine/DNPC/t-NEPC) 출처를 Labrecque 2019 JCI (PMID 31361600), Bluemn 2017 (PMID 29017058), Beltran 2016 (PMID 26855148)로 07장에 1회 인용(PubMed esummary로 서지 확인), 08장 §4·09장 §1.3에서 07장 5절 참조 | fixed
SCI-06-045 | 06_MARKOV_CHAINS_AND_MACROSTATES.md:196-197 | 코드 주석으로 거시상태 이름이 도식의 6개 거시상태와 일치하고 AR+_Adeno(치료 후 재활성 말단)는 초기 Luminal_Naive와 별개, Plastic_Stem은 말단 제외임을 명시(07장 4대 말단 서술과 일치) | fixed
SCI-06-046 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:273-280 | scCNA는 복제수만 측정하므로 양대립 결손 확정엔 변이 호출·LOH 필요, 0.1% 클론은 세포 수·깊이 의존 명시, KRAS 활성화 서술은 가설적 예시로 완화 | softened
SCI-06-047 | 05_CELLRANK2_MULTIKERNEL_INTEGRATION.md:146 | CytoTRACE 0.2/0.85를 예시(가상) 값으로 표기, 가소성 중간체가 높다는 가정은 검증 전임을 명시 | softened
SCI-06-054 | SVG03, SVG05 (텍스트) | eng/SVG: SVG03 max_iter=10, SVG05 가중치 기호화(w1~w3)는 이미 정정됨 확인. 잔여: SVG05의 "k = 4 거시상태 (n_states 지정)"가 06장 n_states=6(말단 4+초기+중간)과 달라 보임, "자동 검증" 문구는 eng 판단 | deferred

## S3 단계 / Track 03 (재배정)

(27건. 수정은 `tracks/03_tme_cell_communication/docs/*.md`만, 새 절 추가 없음. 제어문자 0, `$$` 짝 OK, python ast OK. 확인: PyPI squidpy 1.8.3(Python>=3.12), 기존 audit 계산·liana 소스 근거 재사용. SVG·pipelines·scripts·results는 수정하지 않음.)

ID | 파일:줄 | 조치 | 상태
SCI-03-003 | 00_TME_COMMUNICATION_ROADMAP.md:75-95 | 로드맵 Ch06~09를 실제 파일(06 공간 CCC·LIANA+, 07 MISTy, 08 골전이 5자 회로, 09 트러블슈팅)로 수정, NMF 서술 삭제 | fixed
SCI-03-004 | 00_TME_COMMUNICATION_ROADMAP.md:101-118; 06:184 | Python 3.11로 상향(decoupler 2.x는 >=3.11, liana 1.10은 3.10~3.13), decoupleR->decoupler 표기. squidpy는 ch06 코드가 li.ut.spatial_neighbors로 이미 대체되어 추가하지 않음(squidpy 1.8은 Python>=3.12 요구) | fixed
SCI-03-005 | 00_TME_COMMUNICATION_ROADMAP.md:51 | "전립선암 환자의 80% 이상" -> mCRPC 환자의 다수(보고에 따라 약 80-90%), 전체 환자 비율이 아니며 문헌 확인 필요 명시 | softened
SCI-03-006 | 01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md:169, 186-210 | VEGFR2 -> KDR(HGNC), 복합체 함수에 음수 입력 시 ValueError, 주석을 실제 계산(비음수 스케일 기하평균)과 일치 | fixed
SCI-03-007 | 01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md:29-39, 222-240 | CAF 3분류·국소화가 췌장암 모델(Ohlund 2017, Elyada 2019)에서 정립됨을 단서로 명시, apCAF "TIL 밀집 영역" 국소화 삭제 | softened
SCI-03-008 | 01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md:288; 04 L275-280 상당 | "위양성 완벽 차단" 삭제 -> 위양성 감소/드롭아웃 시 위음성 증가 트레이드오프, 도구별 방식(CellPhoneDB 최소값, CellChat trimean, LIANA expr_prop) 명시 | fixed
SCI-03-012 | 02_DATABASE_CURATION_AND_CONSENSUS.md:28, 35-41 | ICELLNET 약 380개·6계열(antigen binding/checkpoint/chemokine/cytokine/growth factor/Notch), LRdb 3,251쌍(3,085+166)과 구축 방식(FANTOM5·HPRD 등 주석 GO/Reactome 필터) 정정, "수천 편 문헌" 삭제 | fixed
SCI-03-014 | 03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md:113-130, 277 | (b+1)/(B+1)를 유효한 순열 p(Phipson-Smyth)로 표기하고 "비편향·보수적"·"완벽 방지" 삭제, CellPhoneDB/LIANA는 +1 보정 없음(기존 주석 유지); SVG02 해당 문구는 이미 없음 확인 | fixed
SCI-03-015 | 03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md:51-59 | CellChat은 비율 임계가 아니라 trimean(20% 발현이면 0, 30%면 양수)임을 분리 기술, "원천 차단/완전 소거" 삭제, 앰비언트 RNA는 비율 필터로 제거되지 않음 명시 | fixed
SCI-03-016 | 03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md:315-325 | "거의 100% 위양성" -> pi0~1 가정 시 추정 FDR~1이며 단정 불가, "모든 q=1.0" -> 동률 p 가설이 충분히 많지 않으면 q가 1에 가까움(M~130만에서 동률 약 2.6만 개 이상이면 BH 기각 가능; 기존 scipy 계산 근거) | fixed
SCI-03-019 | 04_CELLCHAT_MULTISUBUNIT_MODELING.md:17 | 최소값 연산자 표기에서 "v5~" 삭제(CellPhoneDB v2부터) | fixed
SCI-03-020 | 04_CELLCHAT_MULTISUBUNIT_MODELING.md:123-133, 308-324 | 안드로겐 항상성 루프를 CellChatDB 경로명 아닌 가설 예시로, SPP1의 T세포 IFNg/세포독성 직접 차단·CD44v6-b-catenin을 가설로 표기(RGD-avb3 파골 지원만 확립), 15배 급증은 가상 시나리오 | softened
SCI-03-022 | 05_NICHENET_INTRACELLULAR_SIGNALING.md:136-158, 242-246, 283-288 | v1(Pearson) 대 v2(보정 AUPR) 순위화 구분을 본문에 명시, "필수적 지표" 완화, 관심 유전자는 리간드 반응 유전자가 원 설계이며 세포 정체성 DEG는 계통 차이를 설명한다는 주의 추가(코드는 S2에서 이미 AUPR 정렬) | fixed
SCI-03-025 | 05_NICHENET_INTRACELLULAR_SIGNALING.md:54-56, 200-202; 09:109-110, 124 | "완벽히 교차 검증" -> 보완적 교차 비교, NicheNet prior는 범용 사전 지식이며 세포 유형 특이적 eGRN 아님 명시, "입증" -> 시사, 09장 루브릭 수치(60% 제거, r>0.85, AUROC>0.75, dR2>0.20)는 출처 없는 예시 주석 | softened
SCI-03-027 | 06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:78-103, 159 | LIANA+ 지표를 SpatialDM 방식 Moran R(+Lee L, cosine 등)로 명시하고 Anselin형 표준 정의와 값이 다름을 3.1에 주석, 커널 4종(gaussian/exponential/linear/misty_rbf) 병기, "완벽히 융합" 삭제(li.mt.bivariate.show_functions·docstring 기존 확인) | fixed
SCI-03-028 | 06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:193-201 | S2/ENG에서 ch06 코드가 li.ut.spatial_neighbors(bandwidth, kernel=gaussian)로 대체됨, 현재 코드 재확인 | already-fixed
SCI-03-029 | 06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:30-34, 44-50, 60 | 지수 감쇄는 1D 정상상태 해(2D K0, 3D exp(-r/lambda)/r)임을 명시하고 100-250 um, 10-15 um 수치를 출처 없는 경험칙으로 표기, d0 가이드 박스도 경험적 가이드로 표기 | softened
SCI-03-030 | 06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:321-342 | CXCL12-CXCR4 핫스팟의 "시각적 입증" -> 시사(CXCR4는 골수계·종양 세포도 발현, 가상 수치), "완전 불응"·"증명" 완화 | softened
SCI-03-031 | 06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md:248, 333; 08:15, 41, 183-184, 248, 261 | "용골성" -> "골용해성" 통일(ch06 2곳, ch08 7곳 포함 "파골세포 용골"), 변수명 Osteoclytic -> Osteolytic(ch08 코드 내 일관 변경) | fixed
SCI-03-033 | 07_MULTIVIEW_MODELING_MISTY.md:94-109, 244-286 | LIANA+ MISTy는 RF OOB 예측·계수 clip 후 합 1 정규화·gain.R2 사용임을 본문에 명시하고 문서 식은 교육용 정식화로 표기, 무작위 KFold의 공간 누수·메타모델 in-sample 주석 추가(SVG05는 수정하지 않음) | fixed
SCI-03-034 | 07_MULTIVIEW_MODELING_MISTY.md:214-242 | juxta/para 코드에 플랫폼 해상도 전제(Visium은 15 um 이웃 없음), one-hot 구성만 담는 한계, k=50 절단이 150 um 대역폭을 잘라냄(합성 데이터 50번째 이웃 평균 109 um; 감사 시 실측) 주석 추가 | fixed
SCI-03-037 | 08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md:7-9, 41, 58-64, 261-275 | S2에서 이미 완화된 부분(PTHrP 서술 삭제, RANKL 공급원 맥락 의존)을 확인하고, 골용해성 악순환 패러다임의 1차 근거가 유방암 모형(Guise 1996)이며 전립선암은 ET-1 조골 쪽(Nelson 1995; Yin 2003은 유방암 세포주)임을 구분 명시, TGF-b->ASCL1 단정 완화 | softened
SCI-03-038 | 08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md:181-252 | 회로 점수는 세포 유형별 해석, 회로 1은 파골세포 마커 위주라 종양 세포 상관이 제한적임, TGFB1·TGFBR2 중복을 코드 주석으로 명시하고 종양 상관 분석을 참고용으로 격하 | fixed
SCI-03-040 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:9 | Squair 2021은 차별 발현 유사반복 근거, Dimitrov 2022는 방법 간 낮은 일치로 인용 범위 정정(PMID 34584091, 35680885) | fixed
SCI-03-041 | 09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:241-256, 275-278, 122-125 | 기하평균은 곱의 m제곱근으로 0 규칙만 공유, Harmony는 임베딩 보정이라 발현에 직접 쓰지 않음, Manders 기준 0.6은 근거 없는 예시(PLA 병행), 루브릭 수치는 예시 주석 | softened
SCI-03-042 | pipelines/results_pipeline_01/pipeline_summary_metrics.json; pipelines/01_*.py | eng/pipeline: 요약 JSON의 top_sender가 out-degree가 아닌 out+in 정렬 결과(TAM_M2 vs CSV의 myCAF 11.5168)이며 pipelines/01 스크립트와 결과 파일 수정이 필요해 docs 범위 밖. ch04 Q3의 TAM Top Sender 서술은 가상 시나리오로 표기됨 | deferred
SCI-03-043 | scripts/generate_tme_diagrams.py; images/02,04,05_*.svg | eng/SVG: scripts/generate_tme_diagrams.py 하드코딩 수치(75% 향상, AUROC 0.88/0.79/0.74, 45/20/35%)와 SVG02/04/05는 eng 소관. docs 본문에는 해당 수치 없음 확인, 본문 내분비 ">1 cm"는 그대로 | deferred

비고: 코드펜스 안 `$...$` 수식(eng 통보) 정리 — ch02 L74, ch03 L37·L98, ch04 L11·L129의 ASCII 박스/코드블록에서 LaTeX를 평문으로 치환(예: $S_{\text{obs}}$ -> S_obs). 전 문서 코드펜스 재검색 결과 잔여 0건. ch04 L11 박스의 (C) 기하평균 행에는 CellChat이 trimean 후 적용함을 병기.

## Reference fixes (2026-10-21)

### A. Crossref metadata overrides (`review/reference_overrides.json`)

| Item | Identifier | Source queried | Verdict |
|---|---|---|---|
| Kendall, Gal, Cipolla, CVPR 2018 | DOI:10.1109/CVPR.2018.00781 | arXiv API (id_list=1705.07115): Kendall, Gal, Cipolla; CVF Open Access page of the paper (citation_author order Kendall, Gal, Cipolla). Crossref and OpenAlex (derived from the same deposit) both list Cipolla, Gal, Kendall. | Crossref author order is reversed. Override `author` -> Kendall A, Gal Y, Cipolla R. Chapter 06 (scGPT) block updated. |
| Ross, Lim, Lin, Yang, IJCV 77(1-3):125-141 | DOI:10.1007/s11263-007-0075-7 | Crossref fields of the same record: issued = published-online = 2007-08-17, published-print = 2008-05, journal-issue published-print 2008-05; Semantic Scholar Graph API: year 2008, publicationDate 2008-05-01, volume 77, pages 125-141. Springer landing page and DBLP are bot-blocked from this environment. | 2007 is the online date; the printed issue (vol 77) is 2008, matching the chapter text. Override `year` -> 2008. Chapter 08/06 (IPCA) block updated. |

Mechanism: `format_references.apply_override()` runs inside `format_id()` after the cache read (so apply_references, precommit and the CLI all see it). Cache entries are untouched (the cache diff is 15 added lines, no edits). Supported fields: `author` and `year` for DOI records, plus `source`/`reason` notes; any other field raises an error.

### B. UNRESOLVED items

New ID kind `ARXIV:<id>` in `scripts/format_references.py`: arXiv API `export.arxiv.org/api/query?id_list=...` (3 s gap), raw Atom XML cached under `ARXIV:<id>`, plain-text preprint citation without URL.

| Item | Identifier | Source queried | Verdict |
|---|---|---|---|
| Kingma & Welling 2013 (02_SCVI L93) | ARXIV:1312.6114 | arXiv API id_list: "Auto-Encoding Variational Bayes", Kingma, Welling, 2013-12 | Resolved. Added to [cited]. |
| Dao 2023 FlashAttention-2 (09_TROUBLESHOOTING L48) | ARXIV:2307.08691 | arXiv API id_list: "FlashAttention-2: Faster Attention with Better Parallelism and Work Partitioning", 2023-07 | Resolved. Note: the paper has a single author (Tri Dao), the chapter says "Dao et al. (2023)"; chapter text left unchanged (outside the generated block). |
| Kingma et al. 2014 M2 (03_SCANVI L45) | ARXIV:1406.5298 | arXiv API title search: "Semi-Supervised Learning with Deep Generative Models", Kingma, Rezende, Mohamed, Welling, 2014-06 | Resolved (arXiv version; the arXiv record gives no NeurIPS journal_ref/DOI, so the NeurIPS 2014 venue stays only in the chapter text). |
| Cuturi 2013 Sinkhorn (04_OPTIMAL_TRANSPORT L59) | ARXIV:1306.0895 | arXiv API title search: "Sinkhorn Distances: Lightspeed Computation of Optimal Transportation Distances", Cuturi, 2013-06; journal_ref "Advances in Neural Information Processing Systems 26, 2013" | Resolved; year and venue match the chapter. |
| Niu et al. 2011 Hogwild! (07_GPU_ACCELERATION L173) | ARXIV:1106.5730 | arXiv API title search ("lock-free approach to parallelizing stochastic gradient descent") and DataCite 10.48550/arxiv.1106.5730: "HOGWILD!: A Lock-Free Approach to Parallelizing Stochastic Gradient Descent", Niu, Recht, Re, Wright, 2011 | Resolved. |
| Merck KGaA 2021 (INTR@PID BTC 055 discontinuation; Track 03 ch00, 05, 07, 08) | none | Searched: Crossref/PubMed (no record); web search found company pages (Merck KGaA press-release PDF Aug 2021 "Bintrafusp 055 Study Media Statement", EMD Group news "Bintrafusp Alfa Update" 23-08-2021, EMD Serono topline-data release 16-03-2021 for the 2L BTC study), but the company sites are bot-blocked here, so title/date/authorship could not be verified from the page itself. ClinicalTrials.gov API v2 does confirm the study: NCT04066491, OrgStudyId MS200647_0055, "Gemcitabine Plus Cisplatin With or Without Bintrafusp Alfa (M7824) ... Biliary Tract Cancer", TERMINATED, whyStopped: IDMC review, unlikely to meet the OS primary objective. | Left UNRESOLVED (company announcement, no PMID/DOI/arXiv). The registry record NCT04066491 is a verifiable alternative if a later decision adds an `NCT:` kind; not added now. |

### Checks
- `scripts/apply_references.py` run (7 chapters' blocks updated: 05/02, 05/03, 05/06, 05/09, 06/04, 08/06, 08/07), `--check`: 109/109 OK.
- `python3 scripts/precommit_checks.py --all`: OK (590 files). Master files not rebuilt.

### Changed files
- scripts/format_references.py (ARXIV kind, `apply_override`, docs)
- review/reference_overrides.json (new)
- review/references_cache.json (+5 ARXIV entries, raw Atom; Crossref entries unchanged)
- review/refs/05_generative_foundation_models_02_SCVI_VARIATIONAL_AUTOENCODER.ids, ..._03_SCANVI_SEMISUPERVISED_ANNOTATION.ids, ..._09_TROUBLESHOOTING_AND_BEST_PRACTICES.ids, review/refs/06_manifold_trajectory_04_OPTIMAL_TRANSPORT_WADDINGTON_OT.ids, review/refs/08_large_scale_sparse_07_GPU_ACCELERATION_RAPIDS_SINGLECELL.ids (ARXIV lines added, UNRESOLVED comments removed)
- Reference blocks only in: tracks/05_generative_foundation_models/docs/{02_SCVI,03_SCANVI,06_SCGPT,09_TROUBLESHOOTING}*.md, tracks/06_manifold_trajectory/docs/04_OPTIMAL_TRANSPORT_WADDINGTON_OT.md, tracks/08_large_scale_sparse/docs/{06_INCREMENTAL_PCA,07_GPU_ACCELERATION}*.md

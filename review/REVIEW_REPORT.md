# snRNAseq_learning 학습 자료 검수 리포트

- 검수일: 2026-10-18
- 범위: Track 00(루트) + `tracks/01~10`. 원본 챕터 `docs/*.md` 약 36k줄, 생성물(Master_Book.md / HTML / DOCX) 11세트, 슬라이드 23개, SVG 68개, 파이프라인 24개
- 체계: 오케스트레이터(Opus) 1명과 Sonnet 서브에이전트 2명(`eng-verifier`, `sci-checker`)
- 원칙: 원본은 읽기 전용으로 두고 수정하지 않음(사용자 결정: 리포트만 작성). 빌드와 실행은 scratchpad 복사본에서만 함
- 상태: 검수 최종본 + 수정 이력 §11(S1)·§12(전파 오류)·§13(S2)·§14(S3) + §15 사용자 결정 항목

---

## 1. 결론

**현 상태로는 학습자에게 배포하지 않는 것을 권고합니다.**

- 교재 전체에 인용·참고문헌이 하나도 없습니다(DOI 0건). 검수 결과, 출처가 붙은 수치·데이터셋·논문 서지의 상당수가 틀렸거나 근거가 없었습니다.
- **S1(분석 결정을 잘못 이끄는 오류) 30건**은 오케스트레이터가 전부 독립 재검증했고, 모두 사실로 확인했습니다. 유형은 다음 세 가지입니다.
  - 존재하지 않거나 틀린 데이터셋 accession·DB ID: 학습자가 엉뚱한 데이터를 받게 됨
  - 실행하면 실패하는 코드: 존재하지 않는 함수·인자, 잘못된 좌표
  - 출처를 단 가공 수치
- 수식 유도, 확률모형, 산술 계산은 대체로 정확합니다(§7). 문제는 수식이 아니라 **사실 서술, 도구·API 서술, 인용**에 집중되어 있습니다. 따라서 전면 재작성까지는 필요 없고, 근거를 보강하면서 표적 수정하는 방식으로 해결할 수 있습니다.

## 2. 집계

### 2.1 트랙별 (심각도)

| 트랙 | 주제 | 담당 | S1 | S2 | S3 | 합계 |
|---|---|---|---|---|---|---|
| 00 | snRNA-seq 기초·핵심 파이프라인 | sci | 6 | 33 | 45 | 84 |
| 01 | Multiome / snATAC | sci | 2 | 17 | 30 | 49 |
| 02 | 공간 전사체 | sci | 4 | 17 | 47 | 68 |
| 03 | TME 세포 간 신호 | eng | 0 | 16 | 27 | 43 |
| 04 | 암 유전체·클론 | sci | 2 | 27 | 32 | 61 |
| 05 | 생성 모델·파운데이션 모델 | eng | 1 | 23 | 28 | 52 |
| 06 | 매니폴드·궤적 | eng | 0 | 29 | 25 | 54 |
| 07 | 카운트 통계·차등 풍부도 | sci | 3 | 35 | 31 | 69 |
| 08 | 대규모 희소 연산·GPU | eng | 2 | 20 | 20 | 42 |
| 09 | 업스트림 정렬·BAM | sci | 3 | 31 | 43 | 77 |
| 10 | 프로덕션 파이프라인·MLOps | eng | 2 | 22 | 25 | 49 |
| ENG | 정적 검사·빌드·파이프라인·렌더링 | eng | 5 | 19 | 18 | 42 |
| **합계** | | | **30** | **289** | **371** | **690** |

- ENG-017~019와 ENG-028~029는 에이전트가 번호를 건너뛴 것입니다. 누락된 이슈는 없습니다.
- 확신도: 외부 근거가 있으면 Confirmed, 없으면 Plausible로 표기했습니다. 근거 없는 S1/S2는 Plausible로 낮췄습니다. 각 findings 파일의 마지막 열을 보세요.

### 2.2 차원별

| 차원 | S1 | S2 | S3 | 합계 |
|---|---|---|---|---|
| A 사실 정확성 | 10 | 175 | 244 | 429 |
| B API·코드 | 19 | 93 | 86 | 198 |
| C 파이프라인 | 0 | 8 | 3 | 11 |
| D 생성물 정합성 | 1 | 4 | 5 | 10 |
| E 교육 설계 | 0 | 5 | 24 | 29 |
| F 렌더링 | 0 | 4 | 9 | 13 |

### 2.3 Claim 원장

`claims_ledger.csv`에 **910건**을 기록했습니다.
- 판정 분포: supported 293 / context-dependent 150 / unsupported 123 / wrong 344
- 주의: 이 원장은 의심 가는 주장을 골라 검증한 것이지 무작위 표본이 아닙니다. 따라서 wrong 비율을 교재 전체의 오류율로 해석하면 안 됩니다.

## 3. S1 전체 목록 (오케스트레이터 재검증 완료)

재검증 근거는 `orchestrator_verification.md`에 있습니다.

| ID | 위치 | 내용 | 재검증 방법 |
|---|---|---|---|
| SCI-00-001 | `docs/00_CURRICULUM_AND_ROADMAP.md:115` | GSE210358을 Tang 2022로 표기. 실제는 Chan 2022(JAK/STAT, 사람+마우스) | GEO |
| SCI-00-002 | `docs/00_CURRICULUM_AND_ROADMAP.md:116` | SCP864를 "마우스"로 표기. 실제는 "(Human)" 스터디 | Single Cell Portal API |
| SCI-00-012 | `docs/01_FOUNDATIONS_AND_WET_LAB.md:68` | "snRNA-seq 편향 없음·균일 추출" | 문헌 |
| SCI-00-023 | `docs/02_UPSTREAM_PROCESSING.md:59-67` | 본문은 v3를 설명하는데 STARsolo 명령은 v2 화이트리스트를 쓰고 UMI 길이도 기본값 10 | 동일 챕터 대조 |
| SCI-00-047 | `docs/05_CELL_TYPE_ANNOTATION.md:70` | `genomic_position_from_gtf(adata, gtf_file=…)` 호출 시 TypeError | infercnvpy 0.6.1 시그니처 |
| SCI-00-063 | `docs/06_ADVANCED_INTEGRATION_AND_NETWORKS.md:108` | `sc.pl.dpt_groups_pseudotime(color=…)` 호출 시 TypeError | scanpy 1.12.4 시그니처 |
| SCI-01-001 | `tracks/01…/docs/00_MULTIOME_ROADMAP.md:124-125` | GSE175634는 심근 분화 연구, GSE210358↔Tang 오귀속도 재발 | GEO |
| SCI-01-039 | `tracks/01…/docs/08_EGRN_AND_SCENICPLUS.md:100-143` | `scenicplus.eregulons` 모듈이 존재하지 않음 | GitHub 트리(main, v1.0.0) |
| SCI-02-001 | `tracks/02…/docs/00_SPATIAL_ROADMAP.md:198,201` | GSE98544는 *C. elegans* 연구, E-MTAB-11270은 HCT116 연구 | GEO, BioStudies |
| SCI-02-025 | `tracks/02…/docs/04_CELL_SEGMENTATION_ALGORITHMS.md:195-212` | Baysor CLI·TOML 사용법이 규격과 다름 | 에이전트 근거 수용 |
| SCI-02-030 | `tracks/02…/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:275,306-309` | `N_cells_per_spot` 사용 시 TypeError, 결과 열 이름 사용 시 KeyError | cell2location 소스 |
| SCI-02-051 | `tracks/02…/docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md` | COMMOT 호출 5곳이 무효(`cost_type="quadratic"` 등) | COMMOT 소스 |
| SCI-04-023 | `tracks/04…/docs/04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md:228-239` | `pileup_and_phase.R`에 존재하지 않는 인자 사용 | numbat 소스 |
| SCI-04-024 | 같은 파일 `:259-267` | `run_numbat`의 필수 인자 `lambdas_ref` 누락 | numbat 소스 |
| SCI-05-001 | `tracks/05…/docs/images/01,03,06_*.svg` | 그림 수치(정확도 98.4%, iLISI 0.94)가 트랙 자체 산출물(0.557, 0.556)과 충돌 | 산출물 JSON 대조 |
| SCI-07-001 | `tracks/07…/docs/00_COUNT_STATISTICS_ROADMAP.md:57-72` | Squair 2021에 귀속한 Type I error 표가 논문에 없음 | PubMed 초록 |
| SCI-07-026 | `tracks/07…/docs/03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md:292` | factor 기준 수준 때문에 `coef="conditionMetastasis"`가 존재하지 않음 | R 규칙 + 코드 |
| SCI-07-035 | `tracks/07…/docs/05_FAST_CELL_LEVEL_NEBULA.md:258-290` 외 | nebula `model`/`cutoff_cpc`/`scToNeb(count=)`가 모두 존재하지 않음 | nebula 소스 |
| SCI-08-001 | `tracks/08…/docs/04_TILEDB_SOMA_SINGLE_CELL_DATABASE.md:8,61` | "TileDB-SOMA가 Hilbert 곡선 채택"이라 서술. 기본값은 row-major | tiledbsoma 2.3.0 소스 |
| SCI-08-002 | `tracks/08…/docs/08_GPU_GRAPH_COMMUNITY_DETECTION.md` | CAGRA "비등방성 각도 가지치기"는 논문에 없음 | 에이전트 근거 수용 |
| SCI-09-025 | `tracks/09…/docs/03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md:194-252` | Alevin-fry 스크립트에 오류 3가지가 겹침 | 코드 대조 |
| SCI-09-050 | `tracks/09…/docs/07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md:88-94,182-186` | AR-V7 정션 좌표가 모두 틀려 카운트가 항상 0 | Ensembl REST hg38 |
| SCI-09-068 | `tracks/09…/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md:97-99` | v2 화이트리스트와 `--soloUMIlen 12` 조합(SCI-00-023 재발) | 코드 |
| SCI-10-002 | `tracks/10…/docs/06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md:428` | enzalutamide ChEMBL·RxNorm ID가 각각 trastuzumab, HIV 복합제를 가리킴 | ChEMBL, RxNav API |
| SCI-10-003 | 같은 파일 `:431` | `CL:0002327`은 실제로 mammary gland epithelial cell | OLS4 |
| ENG-004/005/006 | Track 03 ch00·02·05·06 | liana에 없는 API(`li.mt.nichenet`, `spatial_bivariate`, `get_consensus_resources`) | liana 1.10.0 |
| ENG-007 | Track 01 ch07 | `mu.atac.tl.motif_enrichment`/`chromvar`가 존재하지 않음 | muon 0.2.0 |
| ENG-008 | Track 06 ch07 | `cr.pl.fate_probabilities`가 존재하지 않음 | cellrank 2.3.3 |

## 4. 차원별 주요 결과

### A. 과학적 사실 정확성
- **데이터셋 카탈로그는 신뢰할 수 없습니다.** Track 00·01·02의 accession 가운데 전립선 연구가 아닌 것이 5건입니다(선충, 대장암 세포주, 심근 분화, 마우스 GEMM 등). 코호트 규모(n)와 서지(저널·연도·저자)도 여러 곳에서 GEO 기록과 다릅니다.
- **분류 체계의 출처를 잘못 붙였습니다.** 전립선암 4대 아형(AR/amphicrine/DNPC/t-NEPC)은 Labrecque 2019와 Bluemn 2017에서 나온 것인데, 교재 전반이 Tang 2022를 출처로 표기합니다. Tang 2022의 실제 아형은 AR/NE/WNT/SCL입니다.
- **가공 수치를 출처와 함께 제시합니다.** 예: Squair 표, "98%", 벤치마크 표, 골전이 logFC. 학습자가 이 수치를 그대로 인용할 위험이 큽니다.
- **과장·단정 표현이 반복됩니다.** 예: snRNA-seq "편향 없음·전무", pseudobulk "유일한 표준", t-NEPC "비가역", foundation model 섭동 예측의 우월성. 모두 반대 근거가 있는 문헌을 누락했습니다.
- **임상·약리 서술에 오류가 있습니다.** 예: 아비라테론을 LBD 결합제로 서술(실제는 CYP17A1 억제제), 미승인·임상 중단 약제를 "구원 요법"으로 제시, AR·KLK3 좌표 오류.

### B. API·코드 정확성
- S1 19건이 이 차원입니다. 존재하지 않는 함수·인자를 쓰거나 기본값을 틀리게 써서, 본문 코드를 복사해 실행하면 실패합니다.
- 공통 원인은 두 가지입니다. 하나는 그럴듯한 이름을 지어낸 경우(hallucinated API), 다른 하나는 메이저 버전 변경을 반영하지 않은 경우(ENG-009~014)입니다.
- 라이브러리 심볼 188개를 점검한 결과: 정상 136 / 존재하지 않거나 제거됨 16 / import 불가 35(`code_symbol_inventory.csv`).

### C. 파이프라인 실행·재현성 (`pipeline_execution.md`)
- 24개 중 23개 성공, 1개 실패(루트 `pipelines/02_prostate_qc_and_cellbender.py`, scrublet ValueError).
- **버전 고정 파일(requirements·lockfile)이 없습니다(ENG-032).** 원작성 환경(py3.8.16, numpy 1.23.5)에서는 커밋된 결과와 바이트 단위로 같지만, 최신 환경에서는 일부 결과가 달라집니다.
- **개인 conda 경로가 하드코딩되어 있습니다(ENG-024, 파이프라인 5개).** 다른 환경에서는 실패합니다.
- **24개 모두 합성 데이터를 씁니다.** 특히 Track 03~10의 16개는 LIANA, scVI, CellRank 같은 실제 라이브러리 대신 자체 모사 구현입니다. 그런데 README는 이를 "프로덕션 실습 파이프라인"으로 소개하고, 결과 폴더에는 시뮬레이션이라는 고지가 대부분 없습니다(ENG-033).
- 합성 데모의 결과 보고서가 "FDA 21 CFR Part 11 Compliant / VERIFIED"라고 단정합니다(ENG-036, SCI-10-006).
- 합성 결과를 분석 성과처럼 제시하는 사례가 있습니다. Track 06 pipeline02의 "G-PCCA" 거시상태는 정답 라벨을 그대로 복사한 것이고(SCI-06-048), 운명 확률은 1,500개 중 840개가 0.5/0.5 폴백 값입니다(SCI-06-049). Track 05의 그림 수치도 실제 산출물과 충돌합니다(SCI-05-001).

### D. 생성물 정합성
- `*_Master_Book.md` 11개와 `*.docx` 11개는 지금 빌드 스크립트로 다시 만들어도 커밋본과 바이트 단위로 같습니다.
- **HTML은 동기화되어 있지 않습니다.** `scripts/inject_global_track_switcher.py:367-374`가 백업이 있으면 무조건 백업으로 HTML을 되돌리기 때문에, 챕터를 수정해도 HTML에는 반영되지 않습니다(ENG-021). **원본을 수정하기 전에 반드시 먼저 고쳐야 합니다.**
- 최신 mistune으로 다시 빌드하면 HTML에 리터럴 `**`가 989곳 남는데, 빌드 스크립트의 내장 검사 8개는 모두 PASS로 나옵니다(ENG-020). 검사 항목에 빈틈이 있다는 뜻입니다.

### E. 교육 설계·일관성
- **같은 오류가 트랙 사이로 복제되어 있습니다(§5).** 한 곳만 고치면 다른 트랙에 그대로 남습니다.
- 챕터 간·트랙 간 모순: mito% 임계값, 골전이를 용골성으로 볼지 골형성으로 볼지, 신호 도달 거리 d0, WNT5A의 방향, pseudobulk의 지위, Nextflow 캐시 방식.
- 퀴즈 정답·해설 중 일부가 틀렸습니다. 예: Track 00 08 Q3, Track 07 04 Q1, Track 09 01 Q1·05 Q1.
- Track 10 로드맵의 챕터 목록이 실제 파일과 맞지 않습니다(SCI-10-001).

### F. 렌더링 (`findings_engineering.md` §3)
- 콘솔 오류 없음. 검색, 다크모드, 퀴즈 토글, 트랙 스위처, SVG, 수식은 모두 정상입니다.
- 모바일: 트랙 02~10은 반응형 CSS가 없어 375px 폭에서 가로 스크롤이 1290px까지 생깁니다. 넓은 표는 잘려서 읽을 수 없습니다(ENG-038, ENG-039).
- 슬라이드: 1280×720에서 내용이 잘립니다. Track 04/05/10 덱은 8장 중 7장이 잘립니다(ENG-040).
- DOCX: 코드 블록의 줄바꿈 태그 `<w:br/>`가 `<w:r>` 밖에 있습니다. 이미지는 0개이고, 수식에 `_`·`^`가 그대로 남아 있습니다(ENG-041~043). Word에서 실제로 어떻게 보이는지는 확인하지 않았습니다.

## 5. 트랙 간 전파 오류 (일괄 수정 대상)

| 오류 | 발생 위치 |
|---|---|
| GSE210358↔Tang 2022 오귀속, 4대 아형 출처 오귀속 | Track 00 (SCI-00-001, -015), 01 (SCI-01-001), 02 05장, 07 06장 (Track 06은 오귀속은 없으나 인용 자체가 없음: SCI-06-044) |
| v2 화이트리스트(737K)와 v3 설정(UMI 12) 혼용 | Track 00 (SCI-00-023), 09 (SCI-09-068) |
| v3 화이트리스트 크기 3,793,305로 표기(실제 6,794,880) | Track 09 00·01·02장 |
| WNN 거리식 `sqrt(Σw·d²)`(실제는 가중 친화도) | Track 00 (SCI-00-068), 01 05장 |
| 3′ 라이브러리로 TMPRSS2-ERG 접합부를 "반드시 포획"한다는 서술 | Track 04 (SCI-04-036), 09 08장 |
| 가상 결과를 "명확히 도출·입증"했다고 서술 | Track 00 (SCI-00-066), 01 08장, 05, 07 07장 |
| 그림·본문 수치가 자체 파이프라인 결과와 충돌 | Track 05 (S1), 02, 03, 06 (SCI-06-050, -052) |
| 말단 상태를 "비가역"·"영원한 흡수"로 단정(Ku 2017·Mu 2017과 상충) | Track 04 (SCI-04-041 계열), 06 (SCI-06-004) |
| Kfoury·Chen 코호트 규모가 GEO와 다름 | Track 00, 03, 07 |
| 21 CFR Part 11 "준수" 단정 | ENG-036, Track 10 (SCI-10-006) |

## 6. 수정 우선순위 제안 (리포트만 작성했으며 원본은 수정하지 않음)

1. **빌드 경로 수리**: ENG-021 스위처 백업 복원 로직, ENG-020 `**` 잔존과 검사 항목 보강. 이걸 먼저 고치지 않으면 이후 수정이 HTML에 반영되지 않습니다.
2. **S1 30건**: 데이터셋 카탈로그와 DB ID부터 고치고, 다음으로 실행이 실패하는 코드, 마지막으로 가공 수치를 고칩니다.
3. **§5 전파 오류 일괄 수정**: grep으로 모든 발생 위치를 찾아 한 번에 바꿉니다.
4. **참고문헌 체계 도입**: 챕터마다 References 섹션을 두고, 원장(`claims_ledger.csv`)의 evidence 열을 출발점으로 씁니다. 근거를 찾지 못한 수치는 삭제하거나 "가상 예시"로 명시합니다.
5. **재현성**: `requirements.txt`나 lockfile을 만들고(원작성 환경 py3.8 기준 + 최신 환경 기준), 하드코딩 경로를 제거하고, 합성 데이터라는 고지를 README·결과 폴더·요약 보고서에 넣습니다.
6. **S2 → S3** 순서로 처리하고, 렌더링(반응형 CSS, 슬라이드 크기, DOCX 줄바꿈)을 고친 뒤 다시 빌드합니다.
7. **pre-commit 훅 권고**: 링크 검사, 코드 블록 AST 파싱, 빌드 후 `**` 잔존 검사, 라이브러리 심볼 존재 검사(`code_symbol_inventory.csv` 방식). 현재 저장소는 git이 아니므로 `git init`이 먼저 필요합니다.

## 7. 잘된 점

- **수식 유도와 산술**은 대체로 정확합니다.
  - Poisson·Gamma-Poisson→NB 유도, Var=μ+μ²/θ, 샌드위치 팽창, VIF/N_eff
  - KL 폐형식, ELBO 분해, Halko randomized SVD 오차 상한, Chan–Golub–LeVeque 병합식, 모듈러리티 ΔQ
  - CoDa(ΣCov=−Var, Aitchison 거리, CLR), RRA 베타 유도, 스팟 비용·대기열 근사
- Track 05 ch01의 NB/ZINB 손실 구현은 `scvi.distributions`와 1e-5 이내로 일치합니다.
- cBioPortal 재계산 결과 Track 04의 전립선암 CNA 빈도(8p, PTEN, AR 증폭)는 방법을 명시하면 대체로 타당합니다.
- Track 07은 Track 00 4장의 평균-분산 관계 오류를 올바르게 서술합니다.
- 인터랙티브 HTML의 기능(검색, 테마, 퀴즈, 트랙 스위처)은 오류 없이 동작합니다.
- Book.md와 DOCX는 소스와 완전히 동기화되어 있습니다.
- Track 10 ch03의 Snakefile은 Snakemake 9.27.0에서 DAG까지 정상 생성됩니다.
- Track 06의 확산맵·Sinkhorn·UOT 수식과 흡수 Markov 대수(N=(I−Q)⁻¹, B=NR)는 정확하고, CellRank 커널 결합 코드는 실행해 보면 행확률성을 만족합니다.

## 8. 한계

- 실행하지 못한 것: GPU 스택(RAPIDS, cuML, cuGraph), R 패키지(NEBULA, CellChat, NicheNet, Numbat, CopyKAT), cell2location, cellxgene_census, Squidpy 일부. 이 부분은 소스·문서 대조로만 판정했습니다.
- 10x 공식 문서는 429로 차단되어 일부를 간접 확인했습니다. 유료 원문 일부(PMC 미공개)는 초록·2차 문헌으로 대신했고, 해당 항목은 Plausible입니다.
- S3와 Plausible S2는 오케스트레이터가 재검증하지 않았습니다. 재검증 범위는 S1 전체와 S2 표본(트랙당 1~4건)입니다.
- bash·R·groovy 코드 블록의 CLI 플래그는 일부만 검사했습니다.

## 9. 프로세스 기록

- 사용량 한도 때문에 서브에이전트가 세 번 중단됐습니다. 매번 저장된 진행분부터 재개했고, 생성 스크립트는 한 번만 실행하게 하고 ledger 중복이 0건임을 확인했습니다.
- 규칙 위반 2건: `eng-verifier`가 scratchpad 안의 임시 로그 디렉터리와 임시 파일을 `rm -rf`, `rm -f`로 삭제했습니다. 프로젝트 파일에는 영향이 없습니다.
- Playwright가 프로젝트 루트에 만든 `.playwright-mcp/`는 trash로 정리했습니다. 세션 시작 시점에는 없던 폴더입니다.
- Track 03과 06은 부하 분산을 위해 `sci-checker`에서 `eng-verifier`로 재배정했습니다.

## 10. 산출물 색인 (`review/`)

| 파일 | 내용 |
|---|---|
| `REVIEW_REPORT.md` | 본 리포트 |
| `orchestrator_verification.md` | 오케스트레이터의 S1·S2 표본 재검증 근거 |
| `findings_track_00.md` ~ `findings_track_10.md` | 트랙별 상세 이슈 표, 잘된 점, 검증하지 못한 항목 |
| `findings_engineering.md` | 정적 검사, 빌드 비교, 파이프라인, 렌더링 이슈(ENG-*) |
| `pipeline_execution.md` | 파이프라인 24개 실행 매트릭스 |
| `claims_ledger.csv` | 병합 claim 원장(sci + eng). 원본은 `claims_ledger_sci.csv`, `claims_ledger_eng.csv` |
| `inventory.md`, `code_symbol_inventory.csv`, `quant_claims_candidates.csv` | 인벤토리, 코드 심볼 목록, 정량 주장 후보(미검증 후보 목록) |

## 11. 수정 이력 (2026-10-18, 사용자 승인 후)

수정 순서는 리포트 §6의 1~2단계(ENG-021 → S1)입니다. 원본 스냅샷은 `backups/pre_review_fixes_2026-10-18.tar.gz`에 있습니다.

| 단계 | 결과 | 기록 |
|---|---|---|
| ENG-021 | 스위처 주입 스크립트가 백업으로 HTML을 되돌리지 않게 고쳤습니다. 이제 마커 기반으로 다시 주입하며, 여러 번 실행해도 결과가 같습니다. | `fix_log_orchestrator.md` |
| 빌드 환경 | `requirements-build.txt`를 새로 만들어 Python 3.8.16 + mistune 2.0.5로 고정했습니다. 이 환경으로 빌드하면 ENG-020(리터럴 `**` 989곳 → 코드 블록 밖 1곳)과 ENG-022(raw 표 → 0건)가 함께 해소됩니다. | `fix_log_orchestrator.md` |
| S1 30건 | **30/30 원본에 반영**했습니다. sci-checker 20건, eng-verifier 10건이며, 오케스트레이터가 스냅샷 대비 diff를 전수 검토했고 새로 들어간 accession·ID·좌표는 API로 다시 확인했습니다. | `fix_log_sci.md`, `fix_log_eng.md` |
| 재빌드 | 11개 트랙 모두 빌드 assertion 9/9를 통과했고, 스위처 주입과 `verify_switcher_links.py`도 11/11 통과했습니다. 빌드 중 파이프라인 결과 파일에서 바뀐 것은 Track 07의 `execution_time_seconds` 1개 필드뿐입니다. | — |

### 미확정으로 표기한 항목 (값을 지어내지 않음)
- Berglund 2018 데이터 accession (Track 02)
- 공개 전립선 10x Multiome accession (Track 01). 잘못된 두 행은 삭제했습니다.
- Erickson 2022의 플랫폼 세부 사항 (Track 02)

### 남은 작업 (이번 범위 밖)
1. **복제된 같은 오류가 남은 원본 5곳** (§5 일괄 수정 대상):
   - `docs/01_FOUNDATIONS_AND_WET_LAB.md:36` (Tang 귀속)
   - `tracks/02…/docs/09_*.md:198` (`N_cells_per_spot`)
   - `tracks/07…/docs/02_*.md:9,150` (Wilcoxon 80~95%)
   - `tracks/10…/pipelines/01_dsl2_dataflow_and_dag_engine.py:323` (v2 화이트리스트 문맥)
   - `tracks/01…/docs/00_MULTIOME_ROADMAP.md:97` (SCENIC+ 설치 버전 고정)
   - README.md:37과 MASTER_CURRICULUM_ROADMAP.md:201의 Hilbert 서술도 같은 유형입니다.
2. **슬라이드 덱 5개**는 빌드 스크립트 범위 밖이라 옛 오류가 그대로 있습니다: Track 05 2개, Track 08 2개, Track 09 1개.
3. **Track 08 SVG 03·06과 그 생성 스크립트**에 Hilbert·CAGRA 옛 서술이 남아 있습니다.
4. **`backups/`의 옛 HTML·덱**은 이력 보존을 위해 수정하지 않았습니다.
5. **S1 수정 중 새로 확인한 후속 항목:**
   - pychromvar 0.0.4가 anndata ≥0.13과 호환되지 않습니다(주석으로 안내함).
   - Track 08 ch04의 `soma.io.from_anndata` 코드가 anndata 0.13.4와 호환되지 않습니다.
   - Track 08 ch08의 가상 벤치마크 표가 남아 있습니다.
   - SCI-06-023이 남아 있습니다.
6. **S2 289건, S3 371건**은 아직 손대지 않았습니다.

## 12. 수정 이력 2단계: 전파 오류와 잔존분 (2026-10-18)

2단계 수정 전 스냅샷은 `backups/post_s1_fixes_2026-10-18.tar.gz`입니다. 상세 내역은 `fix_log_sci.md`, `fix_log_eng.md`의 "2단계" 섹션과 `fix_log_orchestrator.md`에 있습니다.

| §5 / §11 항목 | 결과 |
|---|---|
| 4대 아형 출처 오귀속 | Labrecque 2019·Bluemn 2017로 정정했습니다(Track 00·02·07). |
| v3 화이트리스트 크기 | 6,794,880으로 바꾸고 파생 계산을 모두 다시 했습니다(Track 09 ch00~02, 슬라이드 01). |
| WNN 거리식 | 가중 친화도 설명으로 바꿨습니다(Track 00 ch07, Track 01 ch05). |
| 3′ 라이브러리 융합 "반드시 포획" | 정정했습니다(Track 04 ch06, Track 09 ch08). |
| 가상 결과 단정 | "가설·가상 예시"로 명시했습니다(Track 00·01·05·06·07·08). |
| 그림·본문 수치 vs 파이프라인 | 산출값으로 맞추거나 "개념도(예시)"로 명시했습니다(Track 02·03·06). |
| Kfoury·Chen 코호트 | GEO 기준으로 정정했습니다. Kfoury는 골전이 9명 + 양성 골수 대조 7명, Chen은 종양 13개·36,424 cells입니다(Track 00·03·07). |
| 21 CFR Part 11 "준수" | "요건을 지원하는 예시, 인증 아님"으로 바꿨습니다(Track 10 docs, pipeline 02 결과). |
| 종말 상태 "비가역" | 모형 가정임을 밝히고 Ku 2017·Mu 2017의 가역성 보고를 반영했습니다(Track 04 ch07, Track 06 ch07·08). |
| §11 잔존 소스 | 5곳 모두 고쳤습니다. 화이트리스트, Track 10 pipeline 01, SCENIC+ `@v1.0.0` 고정, `N_cells_per_location`, "80~95%"의 조건부 서술 포함입니다. |
| 슬라이드 | Track 05(2), 08(2), 09(2) 덱을 HTML로 직접 패치했습니다. |
| Track 08 SVG 03·06 | 생성 스크립트를 고친 뒤 재생성했습니다. |
| 후속 항목 | Track 08 `from_anndata`(import 누락 수정 + anndata 0.12.x 조건 명시, 5천 세포로 실행 확인), ch08 가상 벤치마크 표 명시, pychromvar 대안 안내, SCI-06-023을 처리했습니다. |

검증: 재빌드 11/11 통과, 링크 11/11 통과. 위 오류 패턴 10종은 생성물·슬라이드·소스·SVG 전체에서 0건입니다.

### 남은 작업
- **미확정:**
  - Kfoury 총 세포 수
  - Chen의 별도 11개 샘플이 GSE141445에 포함되는지
  - SCENIC+ `@v1.0.0` 설치 가능 여부(실행하지 않음)
- **후속 후보:**
  - Track 05 `build_decks.py`와 커밋된 덱 사이의 불일치
  - Track 08 SVG05·ch07의 "Tensor Core 100x"와 SVG06 "5.2s"
  - Track 04 CopyKAT SVG·슬라이드의 "98.4%"
  - SCI-06-007·022·048~051
  - SCI-01-026·042·043, SCI-04-036 본문 일부
- **유지:** `backups/`의 옛 HTML·덱은 이력 보존을 위해 수정하지 않았습니다.
- **미착수:** S2(289건)와 S3(371건) 일반 처리.

## 13. 수정 이력 3단계: S2 (2026-10-19)

- 수정 전 스냅샷: `backups/pre_s2_fixes_2026-10-18.tar.gz`
- 이슈별 상태: `s2_status.csv` (289행. `s2_status_{sci,eng,sci_t03,sci_t06}.csv`를 병합한 것)
- 상세 기록: `fix_log_sci.md`, `fix_log_eng.md`의 "S2 단계" 절, 그리고 `fix_log_orchestrator.md`

| 상태 | 건수 | 의미 |
|---|---|---|
| fixed | 208 | 검증된 근거(소스·공식 문서·원 논문·직접 실행)로 정정함 |
| softened | 53 | 근거를 확인할 수 없는(Plausible) 단정을 완화했거나 "예시/가상"으로 표기함. 새로운 주장은 넣지 않음 |
| already-fixed | 26 | S1 단계나 2단계(전파 오류)에서 이미 해소됨 |
| deferred | 1 | ENG-040: 슬라이드 23개를 해상도별로 다시 설계해야 하는 규모라 보류 |
| rejected | 1 | ENG-030: 원작성 환경(py3.8)에서는 정상 실행됨. py3.12에서만 실패 |

**주요 변경 (트랙 공통)**
- 근거 없는 수치·단정을 모두 "예시"로 표기하거나 삭제했습니다.
- 재현되지 않던 시뮬레이션 수치는 다시 실행해 얻은 값으로 바꿨습니다. 예: Track 07 02장 78.0/77.2/5.7%
- API 호출은 설치한 패키지의 시그니처에 맞췄습니다. 대상: scvi-tools 1.5.1, cellrank 2.3.3, scvelo 0.3.4, liana 1.10.0, snapatac2 2.10.0, decoupler 2.2.0 등
- 좌표·서지·DB ID는 Ensembl, UCSC, PubMed, GEO, nf-core 스키마, quay로 다시 확인했습니다.

**엔지니어링**
- 빌드 스크립트 11개
  - mistune 3 이상이면 빌드를 막고, 리터럴 `**`가 남으면 assert로 잡습니다.
  - 반응형 CSS를 넣었습니다(375px에서 페이지 가로 스크롤 없음, 표는 내부 스크롤).
  - DOCX 줄바꿈을 유효한 OOXML로 고쳤습니다.
- 파이프라인 5개에서 개인 conda 경로 하드코딩을 없앴습니다.
- `requirements-pipelines.txt`를 새로 만들었습니다(scRNAseq_anal 환경의 실측 버전).
- 결과 폴더 21곳에 `SIMULATION_NOTICE.md`를 넣었습니다.

**재검증**
- 오케스트레이터가 트랙마다 diff를 검토했고, 새로 넣은 사실은 표본으로 다시 확인했습니다. 근거는 `fix_log_orchestrator.md` 표에 있습니다.
- 재빌드 결과 11/11 트랙이 통과했고, 링크 검증도 11/11 통과했습니다.
- 결과 데이터는 Track 07의 실행 시간 필드 1개만 바뀌었습니다.

**재빌드 방법** (`requirements-build.txt` 참조)

```bash
SCRNA_ENV_SITE_PACKAGES=<conda-base>/envs/scRNAseq_anal/lib/python3.8/site-packages \
  <conda-base>/envs/basic_research/bin/python scripts/build_master_documents.py
```
트랙별 `tracks/*/scripts/build_*_master.py`도 같은 방식으로 실행하고, 마지막에 `scripts/inject_global_track_switcher.py`를 실행합니다.

### 남은 작업
- **S3 371건**: 아직 손대지 않았습니다. 그중 S2 작업 중에 다시 확인한 것은 다음과 같습니다.
  - ENG-001: 깨진 링크 1건
  - ENG-002: `file:///` 절대경로 링크 16건
  - ENG-003: Snakemake 문법 코드 블록에 python 태그가 붙어 있음
- **슬라이드 잔존분**: 덱 안의 SVG 사본이 docs보다 오래된 경우와 새로운 경우가 섞여 있어, 통째로 교체하지 않았습니다.
  - Track 03 덱: SVG01(">1mm"), SVG04(AUROC 0.88/0.79)
  - Track 05 덱: SVG04
- **생성기·덱 불일치**: Track 03·05 슬라이드 생성기와 커밋된 덱이 원래부터 달랐습니다. 이번 치환은 양쪽에 모두 반영했지만, 다른 차이는 그대로 남아 있습니다.
- **ENG-040**: 슬라이드 레이아웃 재설계(보류)
- **§12 미확정 항목**: 그대로 유지합니다.

## 14. 수정 이력 4단계: S3 (2026-10-19 ~ 20)

- 수정 전 스냅샷: `backups/pre_s3_fixes_2026-10-19.tar.gz`
- 이슈별 상태: `s3_status.csv` (371행)

| 상태 | 건수 |
|---|---|
| fixed | 249 |
| softened | 101 |
| already-fixed | 18 |
| deferred | 3 (ENG-027, ENG-044, SCI-03-042: 모두 사용자 결정 필요) |

**주요 변경**
- **링크·경로:** 깨진 링크와 `file:///` 절대경로 링크를 모두 상대경로로 바꿨습니다. 빌드·주입·검증 스크립트는 `__file__` 기준 경로를 쓰도록 바꿔, 저장소 위치가 달라져도 동작합니다.
- **빌드 옵션:** 문서 빌드와 파이프라인 재실행을 분리하는 `--skip-pipelines` 옵션을 추가했습니다.
- **렌더링:** Markdown 파싱 오류, 코드 블록 안의 `$…$`, 검색 스니펫 이스케이프, 모바일 TOC·인쇄 버튼, DOCX 수식 변환을 고쳤습니다.
- **용어·표기:** 챕터 안에서 표기를 통일했습니다. 트랙 사이에서는 AR 변이 표기만 신 번호 기준 F877L로 통일했습니다(구 표기 F876L 병기).
- **슬라이드:** 덱 안의 SVG 사본을 docs와 동기화했습니다(Track 05·06).
- **검증:** 재빌드 11/11 통과, 링크·제어문자·`$$` 오류 0건입니다. 결과 데이터는 의도한 Track 09 보고서와 시간 필드 외에는 바뀌지 않았습니다.

## 15. 전체 수정 요약과 사용자 결정 필요 항목

| 단계 | 대상 | 결과 |
|---|---|---|
| 1 | ENG-021 + S1 | 30/30 반영 |
| 2 | §5 전파 오류 + §11 잔존분 | 반영 완료 |
| 3 | S2 | 289건 처리 (fixed 208, softened 53, already-fixed 26, deferred 1, rejected 1) |
| 4 | S3 | 371건 처리 (fixed 249, softened 101, already-fixed 18, deferred 3) |

### 사용자 결정 필요 (삭제·대규모 변경·결과 데이터 변경이 따르는 항목)

1. **파일 정리 (ENG-027).** 지금은 `.gitignore`만 추가해 두었습니다.
   - `__pycache__` 25개 디렉터리와 `.DS_Store`
   - 루트와 `scripts/`에 각각 있는 `generate_diagrams.py`: 내용이 서로 다르므로 어느 쪽을 남길지 정해야 합니다.
   - 빈 `infercnv_out/`
   - `backups/`
     - 원래 있던 legacy_decks 208K, master_courses_pre_switcher 4.5M, ghost_scripts 32K
     - 이번 수정 단계 스냅샷 4개(각 18M, 합계 약 72M)
2. **결과 데이터 재생성.** 완료했습니다(§17).
   - SCI-02-067: Track 02 Moran's I 분산이 과소 추정돼 Z가 약 1.4배 부풀려집니다.
   - SCI-02-068: Track 02 니치 이름이 순환적으로 붙습니다.
   - SCI-03-042: Track 03 pipeline 01 요약 JSON의 `top_sender`가 CSV와 맞지 않습니다.
3. **ENG-040.** 완료했습니다(§18).
4. **ENG-044.** 완료했습니다(§19).
5. **git 도입과 pre-commit 훅(권고).** 저장소가 git이 아니라서 변경 이력이 스냅샷 tar로만 남아 있습니다. `git init` 후 아래를 pre-commit 훅으로 걸어 두기를 권합니다.
   - 링크 검사
   - 코드 블록 AST 검사
   - 제어문자·`$$` 검사
   - 빌드 assertion

## 16. 저장소 정리와 git 도입 (2026-10-20)

### git 이력
수정 단계별 스냅샷 tar 4개를 차례로 커밋해 이력을 보존했습니다. 정리 전 tar는 합계 72MB였고, `git gc` 후 `.git`은 23MB입니다.

| 커밋 | 내용 |
|---|---|
| `c7425e6` | 원본(검수 수정 전). 당시의 `backups/` 레거시 디렉터리 포함 |
| `c6f00f6` | ENG-021 + S1 30건 |
| `bdc15a0` | 전파 오류·잔존분(§12) |
| `587d27a` | S2 289건(§13) |
| `dbf9c70` | S3 371건 + `review/` 기록 + 빌드 도구(§14–15) |
| `d8a6309` | 구 백업 디렉터리를 트리에서 제거, `.gitignore` 정리 |
| `3a4ebe4` | pre-commit 훅 추가 |

수정 전 상태와의 비교는 `git diff c7425e6 HEAD -- <경로>`로 할 수 있습니다.

### pre-commit 훅
- 구성: `.githooks/pre-commit`이 `scripts/precommit_checks.py`를 호출합니다. 표준 라이브러리만 사용합니다.
- 활성화: clone할 때마다 `git config core.hooksPath .githooks`를 한 번 실행합니다. 이 저장소에는 이미 설정되어 있습니다.

| 수준 | 검사 | 근거가 된 검수 이슈 |
|---|---|---|
| 차단 | 소스 챕터의 제어문자 | S2 단계에서 LaTeX가 깨진 사고 |
| 차단 | `$$` 짝 불일치 | 수식 렌더링 오류 |
| 차단 | python 코드 블록 구문 오류(dedent 후 `ast.parse`) | ENG-003, S1 코드 오류 |
| 차단 | `file://` 링크, 존재하지 않는 상대 링크 | ENG-001, ENG-002 |
| 차단 | `.py` 구문 오류 | — |
| 차단 | 5MB 초과 파일(결과 폴더 제외) | 스냅샷 tar가 실수로 커밋되는 것 방지 |
| 차단 | `results_pipeline_*` 변경. `ALLOW_RESULTS_CHANGE=1`일 때만 허용 | Track 10에서 py3.12로 재실행해 결과가 바뀐 사고 |
| 경고 | 소스 챕터는 바뀌었는데 Master Book/HTML/DOCX를 함께 커밋하지 않음 | ENG-021 유형의 생성물 불일치 |
| 경고 | 생성물만 바뀌고 소스 챕터는 그대로 | 생성물을 손으로 고친 경우 |

훅에서 제외한 것: 전체 빌드와 파이프라인 실행. 수 분이 걸리기 때문입니다. 대신 경고로 재빌드를 안내합니다. 전체 검사는 `python3 scripts/precommit_checks.py --all`로 할 수 있습니다(현재 472개 파일 OK).

### 파일 정리
- **휴지통으로 옮긴 것**
  - 빈 `infercnv_out/`: `pipelines/03_*.py`가 실행 시 `os.makedirs`로 다시 만들며, ignore 처리했습니다.
  - `backups/` 전체: 구 백업 디렉터리 3개와 스냅샷 tar 4개입니다. 모두 git 이력에 들어 있습니다(커밋 `c7425e6`부터 `587d27a`까지). 앞으로 생기는 `backups/`는 ignore됩니다.
- **유지한 것**
  - `__pycache__`, `.DS_Store`: ignore 처리만 했습니다(빌드할 때마다 다시 생깁니다).
- **`generate_diagrams.py` 두 파일 (사용자 승인 후 처리 완료)**
  - 처리: `scripts/generate_diagrams.py`를 남기고, 루트 `generate_diagrams.py`는 휴지통으로 옮겼습니다. 이전 내용은 git 이력에 남아 있습니다. README 디렉터리 트리도 함께 고쳤습니다.
  - 근거:
    - `scripts/` 버전(2026-10-15, 44KB)은 실행하면 현재 `docs/images`의 Track 00 SVG 8개를 byte 단위로 똑같이 재현합니다.
    - 루트 버전(2026-10-13, 27KB)은 SVG 6개만 만들고, 그중 SVG01은 현재 파일과 다릅니다(화살표 마커 등). 07·08 생성 함수도 없습니다.
    - 다른 파일에서 이 두 스크립트를 참조하는 곳은 없습니다.

## 17. 파이프라인 결과 재생성 (2026-10-20, 사용자 승인)

| ID | 변경 | 결과 |
|---|---|---|
| SCI-03-042 | Track 03 pipeline 01: top sender를 out-degree 기준으로 선택하고 top receiver(in-degree)를 추가 | top sender가 TAM_M2 9.3096에서 **myCAF 11.5168**로 바뀌어 CSV와 일치합니다. top receiver는 Endothelial 7.3383입니다. |
| SCI-02-067 | Track 02 pipeline 01: Moran's I 분산을 Cliff & Ord 정규성 가정 정확식으로 교체하고 docstring 정정 | I 값은 그대로이고 Z만 1.3797배 줄었습니다. 예: ALPL 67.57→48.97, KLK3 약 53.2→38.58. `is_svg` 367개는 변하지 않았습니다. sd 검증 결과 정확식 0.01104, 순열 0.01112로 일치합니다. |
| SCI-02-068 | Track 02 pipeline 02: 니치 이름을 데이터 기반 라벨로 변경하고 요약에 한계 명시 | 니치 이름이 바뀌었습니다: Bone_Front→Niche_3 (Osteoblast-dominant), Cancer_Core→Niche_1 (Tumor_AR_Pos-dominant), Desmoplastic_Rim→Niche_0 (CAF-dominant), TLS_Immune_Hub→Niche_2 (CD8_T-dominant). 세포 구성과 수치는 그대로입니다. |

**검증**
- 원작성 환경(scRNAseq_anal, py3.8.16)에서 수정 전 코드로 다시 실행한 결과가 커밋된 결과와 byte 단위로 같았습니다(기준선 확인).
- 수정 후 diff는 해당 지표, 이름, 요약 문구에만 생겼습니다.
- Track 02·03 재빌드가 같은 결과를 byte 단위로 재현했습니다. 빌드 assertion과 링크 검증도 통과했습니다.
- 본문 인용 여부를 확인했습니다. 바뀐 값을 직접 인용하는 docs는 없습니다. Track 03 04장 Q3의 "TAM Top Sender"는 SPP1 경로를 다룬 가상 시나리오로 명시돼 있어 파이프라인 전체 요약과 충돌하지 않습니다.

## 18. 슬라이드 레이아웃 개선 ENG-040 (2026-10-20, 사용자 승인)

덱 23개(184장)는 모두 같은 뼈대를 쓰고 `@media` 규칙이 없었습니다. 덱의 텍스트와 SVG는 그대로 두고, `scripts/inject_slide_layout.py`로 `SLIDE-LAYOUT:BEGIN/END` 마커가 붙은 CSS·JS 레이어를 주입했습니다.
- 데스크톱: 화면에 맞게 자동 축소합니다. 최소 배율은 0.7로, 가독성을 위해 0.6에서 올렸습니다. 0.7로도 다 들어가지 않으면 기존 스크롤 방식으로 돌아갑니다.
- 모바일(768px 이하): 축소하지 않고 세로 스크롤 모드로 바뀝니다. 그리드는 1열이 되고 표·코드는 가로 스크롤합니다.
- 덱 diff에서 삭제된 줄은 0입니다. 다시 실행해도 결과가 같습니다(오케스트레이터가 md5로 확인).

"스크롤이 필요한 슬라이드" 수(전 → 후, 총 184장):

| 해상도 | 토글 열림 | 토글 닫힘 |
|---|---|---|
| 1280×720 | 145 → 72 | 137 → 9 |
| 1920×1080 | 52 → 5 | 1 → 0 |
| 375×812 | 184 → 0 (가로 넘침 182 → 0) | 184 → 0 |

**남은 제약**
- 1280×720에서 다이어그램을 펼친 상태로 보면 72장은 여전히 스크롤해야 합니다. 원인은 기본으로 열려 있는 SVG입니다. 드로어 SVG의 최대 높이를 더 낮추면 줄일 수 있습니다.
- Firefox 126 미만은 `zoom`을 지원하지 않아 기존 스크롤 방식으로 보입니다.
- 생성기 5개로 덱을 다시 만든 뒤에는 `inject_slide_layout.py`를 반드시 다시 실행해야 합니다. 이 순서는 `requirements-build.txt`에 적어 두었고, 빠뜨리면 pre-commit이 경고합니다.

## 19. 공통 UI 템플릿 ENG-044 (2026-10-20, 사용자 승인)

빌더 11개가 각자 갖고 있던 HTML 셸, CSS, JS를 `scripts/course_ui_template.py`의 `render_page(...)` 하나로 합쳤습니다. 빌더 쪽은 −8,076줄, +165줄입니다. 빌드 assertion과 스위처 앵커는 손대지 않았습니다.

- **기능 통일:** 가장 완성도 높은 구현에 맞췄으며, 이 과정에서 빠진 기능은 없습니다.
  - 검색: 하이라이트, 결과 개수, 결과 팝업
  - 다크모드: localStorage에 저장해 새로고침 후에도 유지
  - 퀴즈 펼치기·접기, 인쇄 버튼
  - 검색 정규식의 lastIndex 버그를 고쳤습니다. Track 10에서 누락되던 검색 매치가 178개에서 185개로 복구됐습니다.
- **검증 (HEAD 대비)**
  - Master Book 11개는 byte 동일하고, DOCX 11개는 zip 내부 파트가 동일합니다.
  - HTML 11개의 `<main>` 본문 텍스트와 SVG가 모두 동일합니다.
  - assertion 9/9, 스위처·링크 검증 11/11 통과했습니다.
  - Playwright로 11개 트랙을 1280×800과 375×812에서 확인했습니다. 검색, 다크모드 유지, 퀴즈, 인쇄, 스위처가 모두 동작했고, 모바일 가로 스크롤과 콘솔 에러는 0건입니다.
  - 오케스트레이터가 Track 00·01 화면을 직접 확인했습니다.
- **주의:** 공통 모듈을 고칠 때는 11개 트랙을 모두 다시 빌드한 뒤 `<main>` 본문이 그대로인지 확인해야 합니다.
- **범위 밖 (후속 후보):** 마크다운→HTML 변환 로직(수식 변환, 파서)도 11개 빌더에 중복돼 있습니다. 공통 변환 모듈로 합치는 작업은 이번에 하지 않았습니다.

## 20. 참고문헌 체계 도입 (2026-10-20, 사용자 승인)

검수를 시작할 때 가장 큰 문제였던 "인용 0건"을 해소했습니다. 소스 챕터 109개 전부의 끝에 `## 참고문헌 (References)` 블록을 넣었습니다.

**구성**
- 블록은 `본문 인용 문헌`, `검수 근거 문헌`, `소프트웨어` 순서로 구성됩니다. 내용이 없는 구역은 생략합니다.
- 소프트웨어 항목에 적힌 버전은 검수 때 API를 확인한 버전입니다. 파이프라인 실행 환경은 `requirements-pipelines.txt`에 따로 적혀 있습니다.
- 오프라인 교재 설계와 빌드 assertion 5를 지키기 위해 하이퍼링크 없이 평문으로 표기합니다(`PMID: …`, `doi: …`, 도메인).

**생성 방식 (서지를 손으로 쓰지 않음)**
- `review/refs/*.ids`(챕터별 ID 목록)를 `scripts/apply_references.py`가 읽어, 마커(`REFERENCES:BEGIN/END`)로 감싼 블록을 넣습니다. 다시 실행해도 결과가 같습니다.
- 서지 정보는 `scripts/format_references.py`가 NCBI esummary(PMID)와 Crossref(DOI)에서만 가져옵니다. 받은 원본은 `review/references_cache.json`에 캐시해 같은 결과를 재현합니다.
- 소프트웨어 목록은 `review/software_refs.json`(107개)입니다. 원 논문 PMID는 조회로 확인했고, 출처 도메인은 HTTP 200을 확인했습니다. 10x Cell Ranger 페이지만 접속 제한(429)으로 확인하지 못했습니다.
- 수집하는 동안 에이전트가 기억으로 떠올린 PMID 중 틀린 것들이 있었습니다(LIANA, Squidpy, SciPy, scCODA 등). 모두 조회 단계에서 걸러졌습니다.

**검증**
- 고유 ID 210개를 모두 조회했고 실패는 0건입니다.
- `apply_references.py --check`로 109/109 블록이 캐시에서 다시 만든 결과와 byte 단위로 같았습니다.
- 블록 밖 본문은 HEAD와 동일합니다(109/109).
- 재빌드 결과 11/11 트랙 assertion 통과, 링크 검증 통과, HTML 외부 링크 0개, DOCX 마커 노출 0건입니다. DOCX에서 주석 단락을 건너뛰도록 빌더 11개를 고친 결과입니다.
- Crossref에 저널명이 없는 항목(책·프리프린트)이 `**`로 렌더링되던 문제를 고쳤습니다.
- pre-commit은 블록을 손으로 고치면 차단합니다(`.ids`·캐시와 불일치하면 BLOCK). 소스 챕터에 블록이 없으면 경고합니다.

**남은 한계**
- **UNRESOLVED:** 없음. 6건 중 5건은 arXiv ID, Merck KGaA 2021은 임상시험 등록번호로 보완했습니다(§22).
- **Crossref 메타데이터 오류:** Kendall 2018 저자 순서와 Ross IPCA 연도는 override로 바로잡았습니다(§22).
- **문헌 판단에 가정을 둔 3건:** `.ids` 주석에 적어 두었습니다(Wood 2017, Cusanovich 2018, Prekovic 2016).
- **본문 내 위치 표시 없음:** 본문 안에 번호식 인용 마커는 넣지 않았습니다. 본문은 기존 저자-연도 표기를 그대로 쓰고, 블록과는 저자명으로 대응됩니다.

## 21. 마크다운→HTML 변환 로직 통합 (2026-10-21, 사용자 요청)

§19에서 후속 후보로 남긴 변환 로직 중복을 정리했습니다.

**구조**
- `scripts/course_markdown.py`(신규)에 LaTeX 변환(`LATEX_SYMBOLS` 196키, HTML·Unicode 변환기), Book/HTML/DOCX 생성, 공통 검증을 모았습니다.
- 빌더 11개는 `TrackConfig`, 트랙별 파이프라인 assertion, `main()`만 남겼습니다(16,261줄 → 2,030줄).
- 파일 이름, CLI(`--skip-pipelines`), 출력 경로, mistune 가드, DOCX 주석 건너뛰기, 마커, assertion 5는 그대로입니다.

**출력 변화** (Master Book 11개는 전 과정에서 byte 동일)
- `d355177` 통합: T01–T05가 T06–T10의 LaTeX 처리를 그대로 물려받았습니다(`\partial`→∂, `\nabla`→∇, `\operatorname`·`\text` 중괄호 제거, `\dot`, `\mathbf` 등).
  - 결정 2건: `\parallel`→∥(T05 KL 9곳), `\quad`→공백 2칸.
  - 여러 글자에 걸친 `\overline`이 마지막 글자 뒤에 떨어져 붙던 문제를 고쳤습니다(T01·T02, 그리고 HEAD에도 있던 T06).
  - sci-checker 검토(`review/converter_dedup_review.md`): 바뀐 수식 54개 중 후퇴 0건.
- 후속 커밋: HEAD부터 있던 렌더링 결함 5건을 공통 모듈에서 고쳤고, 11개 트랙 전부에 반영됩니다.
  1. `\|`가 백슬래시+파이프나 `|`로 남던 것 → ∥(노름·KL 구분자). 표 안 이스케이프로 쓰인 경우가 없음을 확인했습니다.
  2. `\vec{v}`가 `{v}`로 나오던 것 → v⃗
  3. `\tfrac`가 원문 그대로 남던 것 → 분수
  4. `\text{}` 안의 `\_`가 아래첨자가 되던 것(`n_(f)ragments`) → 밑줄
  5. DOCX 다글자 윗줄이 `_ ^ ( )`에도 붙던 것 정리
- 검증: 11/11 빌드 assertion 통과, 스위처·링크 검증 통과, pre-commit `--all` OK, `results_pipeline_*` 변경 없음. 바뀐 HTML/DOCX 덩어리는 모두 위 항목으로 분류됐고 그 밖의 변화는 0건입니다. 수정별로 브라우저 화면을 확인했습니다.
- 남은 것: `\underbrace`는 아래첨자로 평탄화됩니다(의미는 유지). `\Vert`는 ‖로 둡니다.
- **주의:** 공통 모듈을 고치면 11개 트랙을 모두 다시 빌드하고 Book byte 동일성과 HTML/DOCX diff를 확인해야 합니다.

## 22. 파이프라인 모드 검증과 참고문헌 보완 (2026-10-21, 사용자 요청)

### 파이프라인 모드 빌더 검증
§21 리팩터 이후 한 번도 실행되지 않았던 파이프라인 재실행 경로를 검증했습니다. 프로젝트를 건드리지 않도록 HEAD(`bdcbe5f`)의 로컬 clone에서 실행했습니다.
- 환경은 §17과 같습니다: `basic_research` python + `SCRNA_ENV_SITE_PACKAGES`(scRNAseq_anal), `PYTHONNOUSERSITE=1`, `PYTHONHASHSEED=0`.
- 11개 빌더가 모두 9/9 assertion을 통과했습니다. assertion 4는 실제 재실행 형태로 통과했습니다(T02–T07 기본 모드, T08–T10은 `--re-run-pipelines`). Book·HTML은 byte 동일하고 DOCX는 zip 파트가 동일합니다.
- `results_pipeline_*`: T02–T06은 byte 동일합니다. T07–T10에서는 16개 파일이 바뀌었는데, 모두 실행 시각·벤치마크 시간·처리량과 `python_version`(3.8.18→3.8.16) 필드입니다.
  - 리팩터가 원인이 아님을 확인했습니다. 리팩터 전 빌더(`fa0cb97`)로 돌려도 같은 파일 집합이 바뀌고, 같은 파이프라인을 두 번 돌려도 같은 파일이 바뀝니다. 숫자 결과와 모델 가중치는 그대로입니다.
  - 커밋된 결과는 바꾸지 않았습니다.
- **주의:** 실제 트리에서 기본 모드로 빌드하면 T07–T10 결과 파일의 시간 필드가 바뀝니다. 커밋에 `ALLOW_RESULTS_CHANGE=1`이 필요하므로, 결과를 다시 만들 의도가 없을 때는 `--skip-pipelines`로 빌드합니다.

### 참고문헌 보완
- **Override 체계:** `review/reference_overrides.json`(ID별로 고친 필드와 출처·사유)을 `format_references.py`가 캐시를 읽은 뒤 적용합니다. 캐시는 API 원본 그대로 유지합니다. 허용 필드는 `author`, `year`입니다.
  - Kendall 2018 CVPR: Kendall A, Gal Y, Cipolla R(arXiv 1705.07115·CVF 페이지 기준. Crossref·OpenAlex는 역순).
  - Ross IPCA: 2008(Crossref published-print 2008-05, vol 77, issue 1-3. 2007은 온라인 게재일).
- **arXiv ID 지원:** `ARXIV:` 종류를 추가했습니다. arXiv API에서 받아 캐시하고 URL 없이 평문으로 표기합니다(assertion 5 준수). 해결 5건, 모두 arXiv API 실조회로 확인:
  - Kingma & Welling 2013 `1312.6114`, Kingma 2014 M2 `1406.5298`, Cuturi 2013 `1306.0895`, Niu 2011 Hogwild! `1106.5730`, Dao 2023 FlashAttention-2 `2307.08691`
- **본문 수정 1건:** FlashAttention-2는 단독 저자이므로 Track 05 ch09의 "Dao et al. (2023)"을 "Dao (2023)"으로 고쳤습니다.
- **Merck KGaA 2021 발표(Track 03, 4곳):** 발표문 자체는 PMID·DOI가 없고 회사 페이지는 자동 접근이 막혀 확인하지 못했습니다. 대신 해당 임상시험의 등록 기록을 넣었습니다. `NCT:` 종류(ClinicalTrials.gov API v2, URL 없는 평문)를 추가했고, `NCT04066491`(INTR@PID BTC 055, Phase 2/3, IDMC 검토 후 중단)이 본문의 무용성 중단 서술을 뒷받침합니다. 본문 4곳의 표기는 "(Merck KGaA 2021; NCT04066491)"로 바꿔 블록과 연결했습니다.
- **검증:** `apply_references.py --check` 109/109, 본문(블록 밖)은 Dao 1곳 외 HEAD와 동일, 재빌드 11/11 assertion 통과, 링크 검증 통과. 출력이 바뀐 트랙은 05·06·08뿐입니다.

## 23. 파이프라인 결과 재현성 검사 (2026-10-21, 사용자 요청)

§22에서 확인한 문제(빌더가 파이프라인을 재실행할 때마다 T07–T10 결과 파일의 시간 필드가 바뀌어 git 트리가 더러워짐)를 정리했습니다.

**동작 변경**
- 빌더는 트랙의 `pipelines/`를 임시 폴더에 복사해 그곳에서 파이프라인을 실행하고, 결과를 커밋된 `results_pipeline_*`와 비교합니다. 빌드가 git 트리를 바꾸지 않습니다. 임시 사본은 빌드 후 휴지통으로 보냅니다(최대 약 20MB/트랙, 사용자 결정).
- assertion 4는 "실행 성공"에서 "실행 성공 + 커밋된 결과 재현"으로 강화됐습니다. 불일치하면 파일과 첫 번째로 다른 필드를 출력합니다.
- 플래그 의미는 그대로입니다(T02–T07 기본 재실행, T08–T10은 `--re-run-pipelines`, `--skip-pipelines`는 건너뜀). 결과를 의도적으로 다시 만들 때만 `--update-results`를 씁니다.
- 학습자가 파이프라인을 직접 실행하는 방식은 바뀌지 않았습니다. 시간·처리량 값은 각자의 기기 측정값이 나오는 것이 정상입니다.

**휘발 필드 선언** (`scripts/pipeline_volatile_fields.json`, 항목마다 사유 기재)
- 비교에서 제외하는 값은 실행 시간, 타임스탬프, 처리량, 시간 비율로 계산한 속도비, peak RAM, 측정 시간에 따라 정해지는 통과 판정(SLA 등), 실행 환경 문자열(python 버전, OS, CPU)뿐입니다.
- sci-checker가 선언 항목마다 생성 코드를 확인했습니다(`review/volatile_fields_review.md`). 과학적·통계적 결과는 하나도 제외되지 않았습니다. 리뷰에서 나온 누락 선언(시간 기반 판정, 환경 문자열)과 과잉 마스킹 1건을 반영했습니다.
- 선언한 키·열·행이 실제 파일에 없거나 CSV 모양이 다르면 비교가 실패합니다(조용히 넘어가지 않음).

**리뷰에서 함께 찾은 출력 오류 수정**
- 파이프라인 요약문에 측정값과 맞지 않는 고정 문구가 있었습니다. 이제 코드가 측정값에서 계산합니다.
  - T10 02: 10,000세포 처리 시간 `~3.8 seconds`(고정값) → 측정 처리량에서 계산(커밋본 기준 `~1.2 seconds`)
  - T10 01: `~1.28s / ~0.02s / 64x` → cold/warm 실측에서 계산(커밋본 기준 `~0.39s / ~0.00s / >788x`)
  - T08 02: `speedup_target_pass` 조건식이 챕터 기준(10–50배)과 달랐던 것을 수정했습니다. 요약의 `PASSED` 고정 문구는 실제 검사 결과에서 만듭니다.
- 커밋된 결과는 재생성하지 않았습니다. 재생성하면 모든 시간 필드가 바뀌기 때문입니다. 바뀐 문구 3곳만 커밋된 측정값으로 계산해 고쳤으므로, 각 파일은 한 번의 실행과 계속 일치합니다.

**검증**
- 11개 빌더 기본 모드와 T08–T10 `--re-run-pipelines`가 실제 트리에서 통과했습니다. 두 번째 실행에서도 결과 변경은 0건입니다.
- 음성 시험: 비휘발 값 변경은 실패, 휘발 필드만 변경은 통과, 선언 키 이름 변경과 CSV 열 이름 변경은 실패합니다.
- CPU 부하 상태에서 T10 02를 재실행해도 통과합니다.

## 24. 남은 렌더링·환경 항목 정리 (2026-10-21, 사용자 요청)

- **`\underbrace` (3곳: T03 ch04, T05 ch02, T07 ch04):** 아래첨자로 평탄화되던 것을, HTML에서는 식 아래 곡선 괄호와 가운데 정렬 라벨로 표시합니다(인라인 스타일, 다른 출력 영향 없음). DOCX에서는 `[식]_(라벨)`로 표기합니다.
- **`\Vert`:** 변경하지 않았습니다. T10 ch06·ch08에서 해시 체인의 연결(concatenation) 연산자로 쓰이며, ‖(U+2016)가 표준 표기입니다. 노름을 뜻하던 `\|`(§21)와 용도가 다릅니다.
- **괄호 없는 단일 토큰 첨자:** `q_\phi`, `\sigma_\epsilon^2`, `\mathbb{R}_+^N`처럼 명령어·기호 하나를 첨자로 쓴 경우 HTML에 `_`·`^`가 그대로 보이던 결함(HEAD부터 존재, T05에 집중)을 고쳤습니다. LaTeX와 같이 첨자는 토큰 하나를 받습니다. 수식 안 리터럴 `_`·`^` 242곳 → 0곳. T01–T08, T10의 HTML이 바뀌었고, DOCX는 기존 표기 규약(`q_ϕ`)이 이미 일관되어 바뀌지 않았습니다.
- **T10 audit manifest:** pinned 환경에서 T10 pipeline 02 결과를 다시 만들어 `python_version`이 3.8.16이 되었습니다. 바뀐 값은 선언된 휘발 필드뿐이며 모델 가중치 해시는 그대로입니다.
- **다른 환경에서의 재현:** 비교는 pinned 환경·기준 플랫폼에서 정확하며, 다른 BLAS·libm·CPU에서는 모델 해시나 손실이 다를 수 있다는 안내를 `requirements-pipelines.txt`와 assertion 4 불일치 메시지에 넣었습니다. 비교 기준은 완화하지 않았습니다.
- 검증: 11/11 빌드 통과(Book byte 동일), T10 `--re-run-pipelines` 결과 변경 없음, 스위처·링크·pre-commit 통과. 수정별 화면을 Playwright로 확인했습니다.

## 25. 본문이 코드 블록으로 보이던 결함 (2026-10-21)

§24 검증 중 발견했습니다. 목록 항목 아래 4칸 들여쓴 `$$…$$` 수식과 설명문(정상적인 목록 연속 구조)에서, 수식 전처리가 수식 줄의 들여쓰기를 지워 목록이 끊기고, 이어지는 들여쓴 설명문이 코드 블록으로 렌더됐습니다. 학습자 화면에는 회색 코드 상자 안에 `<span class=…>` 같은 원시 마크업이 보였습니다. HEAD부터 있던 결함이며 원문 작성 방식이 아니라 변환기 문제였습니다.

- 수정: display 수식이 줄 앞에서 시작하면 그 줄의 들여쓰기를 보존합니다(`scripts/course_markdown.py`). 챕터 원문은 고치지 않았고 Master Book은 byte 동일합니다. DOCX는 영향이 없었습니다(펜스로만 코드를 만듦).
- 결과(11개 HTML 합계): 본문이 코드 블록이던 곳 42 → 0, 화면에 노출된 이스케이프 수식 마크업 53 → 0. 원문의 펜스 코드 블록(언어 미표기 323, 언어 표기 217)은 그대로입니다. 수식 때문에 쪼개지던 목록도 하나로 복원되어 번호가 이어집니다(항목 간 여백이 약간 넓어질 수 있음).
- 검증: 11/11 빌드 통과, 스위처·링크·pre-commit 통과, 수정 위치 Playwright 화면 확인(T05 ch03, T07, T08).

## 26. GPU 코드 블록 실행 검증 (Colab T4, 2026-10-27, 사용자 요청)

GPU가 필요한 챕터 코드 블록 7개를 Google Colab T4(Python 3.13.16, torch 2.11.0+cu130, scvi-tools 1.5.1, cell2location 0.1.5, rapids-singlecell 0.18.1, RAPIDS 26.06)에서 실행했습니다. Colab CLI로 세션을 열고 노트북(`review/colab/gpu_blocks_test.ipynb`)을 실행한 뒤 세션을 닫았습니다. 3회 실행에 약 0.6 compute unit이 들었습니다. 결과 파일은 `review/colab/gpu_results_2026-10-27_run{1,2,3}.json`입니다.

| 실행 | 통과 | 실패와 원인 |
|---|---|---|
| 1 (수정 전) | 4/7 | B3: `import cell2location` 실패 — pyro-ppl 1.9.2에 `deep_getattr`가 없음. B4·B6: `rmm.mr.parse_size`가 없음 |
| 2 (pyro 고정 + RMM 문자열) | 5/7 | B3: `use_gpu`를 scvi Trainer가 거부. B4: rapids-singlecell `calculate_qc_metrics`에 `percent_top` 인자가 없음 |
| 3 (모두 수정) | **7/7** | 없음 |

**수정한 챕터 원문**
- T08 ch07·ch09: `rmm.mr.parse_size("64GB")`를 `rmm.reinitialize`에 문자열로 직접 전달합니다. `rmm.reinitialize`는 `int | str`을 받고 문자열을 자체 해석합니다(GB=10^9, GiB=2^30 bytes). `rmm.mr.parse_size`는 RMM 소스에 없습니다.
- T08 ch07: `rsc.pp.calculate_qc_metrics(..., percent_top=None, ...)`에서 `percent_top`을 제거했습니다. 설치된 rapids-singlecell 0.18.1 시그니처에 이 인자가 없습니다. 다른 챕터의 `percent_top=None`은 scanpy 호출이라 유효하므로 그대로 둡니다.
- T02 ch05: `use_gpu=True`를 `accelerator="gpu"`로 바꿨습니다(`train`과 `sample_kwargs` 모두). 설치된 scvi-tools 1.5.1의 `train`·`sample_posterior`는 `accelerator`/`device`를 받고 `use_gpu`는 받지 않습니다.
- T02 ch00: 설치 줄에 `"pyro-ppl<1.9.2"`를 추가하고 이유를 주석으로 적었습니다. cell2location 0.1.5는 `pyro-ppl>=1.8.0`만 요구하고, pyro-ppl 1.9.1 소스에는 `deep_getattr`가 있으며 dev 브랜치에는 없습니다. 이 고정은 임시이며, cell2location이 새 pyro를 지원하면 제거해야 합니다.

**이 검증이 보장하는 범위와 한계**
- 합성 데이터와 축소 설정으로 돌렸습니다. B3는 `max_epochs` 250→20, 30000→300, B4는 RMM 풀 60/76GB→1/4GB, B6은 64GB→2GB입니다. 챕터 기본값(80GB GPU 대상)은 T4(15GB)에서 돌릴 수 없어 실행하지 않았습니다.
- 통과는 "코드가 오류 없이 끝났다"는 뜻이며 출력 수치의 생물학적 타당성은 검증 대상이 아닙니다.
- 패키지 버전은 pyro-ppl 외에는 고정하지 않고 2026-10-27 기준 최신으로 설치했습니다. 이후 버전이 바뀌면 다시 깨질 수 있습니다.
- B2의 `torch.cuda.amp.autocast/GradScaler` 사용은 torch 2.11에서 FutureWarning(deprecated)만 내고 동작합니다. 바꾸지 않았습니다.
- GPU 블록 이외 코드 블록은 이 검증에 포함되지 않습니다.

## 27. 80GB급 GPU 기본값 검증과 CAGRA `max_queries` 결함 (Colab, 2026-10-27, 사용자 요청)

§26의 T4 검증은 챕터의 80GB 기준 기본값을 줄여서 돌렸습니다. 이번에는 80GB급 GPU에서 챕터 기본값 그대로(RMM 풀 60GB/76GB, 64GB, 메모리 가드 70/65GB) 합성 데이터 25k·100k세포 × 35,000유전자로 돌렸습니다. 100만 규모는 호스트 RAM 위험 때문에 100k에서 멈췄습니다. 노트북은 `review/colab/gpu_blocks_80gb.ipynb`(생성기 `build_gpu_80gb_notebook.py`)이고 블록마다 별도 프로세스로 실행하며 최대 GPU·호스트 메모리를 기록합니다.

**GPU 확보 결과**

| 요청 | 결과 | 사용률 |
|---|---|---|
| A100 | 열림, 그러나 `A100-SXM4-40GB`(40GB)라 76GB 풀이 들어가지 않아 닫음 | 5.30 CU/hr |
| H100 | 이 계정에 권한 없음(백엔드가 거부), 세션 생성 안 됨 | - |
| G4 | `RTX PRO 6000 Blackwell Server Edition` 97,887 MiB, 호스트 RAM 약 176GB(CPU 48코어) | 8.90 CU/hr |

**G4(96GB급) 전체 실행: 5/6 통과** (`gpu_results_2026-10-27_80gb_g4.json`)
- 통과: B4(25k, 100k: 풀 60GB/76GB 초기화와 전체 파이프라인), B5 25k, B7 100k, B6(64GB 관리형 풀 + 4GiB 할당). 호스트 최대 RSS는 약 4.4GB였습니다.
- 실패: B5 100k. `CuvsException: cudaLaunchKernelExC … invalid argument`(cuVS CAGRA 검색 커널 실행).

**B5 100k 실패의 원인 추적**
- 같은 실패가 A100(40GB)에서도 재현되어(`…_b5_a100.json`) Blackwell 고유 문제가 아닙니다. 25k점에서는 두 GPU 모두 통과했습니다. §26 T4 실행에서는 20,000점으로 통과했기 때문에 가려져 있었습니다.
- 챕터는 `cagra.SearchParams(max_queries=100_000, …)`를 명시합니다. 같은 노트북에서 이 값만 바꿔 A100으로 재실행했습니다.
  - `max_queries=50000`: 통과(100k점이 2배치로 나뉨, `…_max_queries_50000.json`)
  - `max_queries=0`: 통과(`…_max_queries_0.json`)
  - `max_queries=100_000`(챕터 원문): 실패
- cuVS 소스의 docstring은 `max_queries`를 "Auto select when 0", 기본값 0으로 적고 상한은 적지 않습니다. 메커니즘은 확인하지 못했습니다(CUDA 그리드 한 축의 한계 65,535를 넘는 배치가 원인일 수 있다는 것은 가설입니다). 확인된 것은 100,000에서 실패하고 0과 50,000에서 통과한다는 사실뿐입니다.

**수정한 챕터 원문** (T08 ch08)
- `max_queries=100_000` → `max_queries=0`(자동 선택)과 실측 근거 주석.
- 블록 머리 주석의 "GPU 미실행"을 "합성 100k점으로 A100에서 실행 확인, 100만 규모는 미실행"으로 고쳤습니다.
- 재검증은 치환 실행(`max_queries=0`)으로 했고, 챕터 원문은 주석만 다릅니다. 수정된 원문 자체를 다시 GPU에서 돌리지는 않았습니다.

**이 검증의 한계**
- 80GB가 아니라 96GB급 Blackwell에서 돌렸습니다. A100/H100 80GB는 이 계정으로 확보하지 못했습니다. 메모리 기본값의 동작은 확인했지만 챕터의 A100 벤치마크 수치(가속 배율)와는 비교할 수 없습니다.
- B4의 "최대 GPU 사용량 57,981 MiB"는 RMM 풀 초기 예약(60GB)이며 실제 사용량이 아닙니다.
- B6의 90GB 관리형 할당 시험은 97GB GPU에서 성공했으므로 물리 메모리 초과(오버서브스크립션) 주장은 확인하지 못했습니다.
- 첫 실행 시간은 JIT 컴파일을 포함하므로(같은 설정이 17.9초와 4.4초) 단계별 시간을 비교하면 안 됩니다.
- 1M 세포·1M점 규모는 실행하지 않았습니다. B7의 메모리 가드 방어 분기(70/65GB)는 발동하지 않았습니다.
- 4개 세션 합계 약 0.4 compute unit을 썼습니다(잔량 342.04 → 341.64, 표시는 지연될 수 있음). 모든 세션은 종료했습니다.

## 28. 공개 준비 (GitHub `snRNAseq-learning`, Pages, 2026-11-20, 사용자 결정)

GitHub Pages는 개인 계정(Free, Pro)에서 공개 저장소만 지원하므로 저장소를 공개합니다. 공개 전 점검(추적 파일 606개, 커밋 28개)에서 비밀 정보, 제3자 이미지·폰트, 외부 리소스, 환자 데이터는 발견되지 않았습니다.

- **로컬 경로 제거**: `review/` 6개 파일과 `requirements-build.txt`의 개인 경로를 `<repo>`, `<scratchpad>`, `<conda-base>`, `/Users/<user>/` 자리표시자로 바꿨습니다. 검수 기록의 의미(어떤 종류의 경로였는지)는 유지됩니다.
- **연락처 이메일 제거**: `scripts/format_references.py`는 연락처를 `REF_MAILTO` 환경변수에서만 읽습니다. 설정하지 않으면 주소를 보내지 않습니다. 캐시가 있어 재실행에 네트워크가 필요 없습니다.
- **재발 방지**: `scripts/precommit_checks.py`가 로컬 경로(`/Users/<name>/`, 세션 scratchpad, conda Caskroom)와 개인 메일 주소가 든 파일의 커밋을 막습니다(`--all`에도 적용).
- **`.hgp/`**: 로컬 도구 상태 파일이라 추적을 끊고 `.gitignore`에 넣었습니다.
- **DOCX 작성자**: 메타데이터 작성자를 `wgsim`으로 바꾸고 11개 트랙을 다시 빌드했습니다. 바뀐 것은 DOCX 11개의 `docProps/core.xml`뿐이고 Book·HTML은 byte 동일합니다.
- **라이선스와 면책**: `LICENSE`(MIT, © 2026 wgsim, 파운데이션 모델 자료와 같은 문구), README 상단 면책 안내(교육용, 임상 사용 금지, 합성 데이터, AI 도구 활용).
- **Pages 진입점**: 루트 `index.html`(트랙 목록, 면책 안내, 라이선스)과 `.nojekyll`. Pages 설정은 `main` 브랜치 루트를 소스로 씁니다.
- **기록**: 공개 저장소에는 현재 트리만 담은 기록 1개짜리 브랜치를 올립니다. 이전 커밋(이메일·로컬 경로 포함)은 로컬에만 남습니다.

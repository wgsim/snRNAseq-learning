# Chapter 09. 업스트림 정렬 트러블슈팅 및 베스트 프랙티스 (Troubleshooting & Alignment Best Practices)

> **Cognitive Goal**: 대규모 단일핵 전사체(snRNA-seq) 업스트림 정렬 파이프라인에서 발생하는 시스템 수준의 장애 요소를 체계적으로 진단하고 해결하는 엔지니어링 역량을 배양한다. 인공 트랜스진(GFP, Cas9) 및 키메라 융합 유전자(*TMPRSS2-ERG*)를 포함하는 맞춤형 참조 유전체(Custom Reference Genome) 구축 및 STAR 색인화 수식을 정립한다. 단일핵 시퀀싱 고유의 핵 내 미성숙 전구체 RNA(Pre-mRNA) 인트론 리드를 온전히 포획하기 위한 STARsolo 파라미터(`--soloFeatures GeneFull`) 최적화 원리를 체득한다. 다중 매핑 리드(`NH > 1`)의 기대값 최대화(EM) 확률 분배 알고리즘과 가유전자(Pseudogenes) 노이즈 제어 기법을 도출한다. 현업에서 빈번히 발생하는 10대 실패 모드(Top 10 Failure Modes)를 분석하고, 5단계 엔드투엔드 업스트림 캡스톤 프로젝트 평가 루브릭을 완성한다.

---

## 1. 맞춤형 참조 유전체(Custom Reference Genome) 구축 가이드

단일 세포 암 연구에서는 표준 참조 유전체(GRCh38)에 존재하지 않는 외인성 서열(Exogenous sequences)을 정량화해야 하는 경우가 빈번합니다:
- **인공 리포터 및 트랜스진**: eGFP, mCherry, Firefly Luciferase, SpCas9, Cre-recombinase.
- **종양 바이러스 게놈**: 고위험군 인유두종바이러스(HPV16/18), EBV, B형/C형 간염 바이러스.
- **합성 키메라 전사체**: *TMPRSS2-ERG* 인공 엑손 융합 주석.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        맞춤형 참조 게놈 구축 및 STAR 인덱싱 워크플로우                │
│                                                                                        │
│   [1. FASTA 서열 병합]                                                                 │
│   GRCh38.primary_assembly.genome.fa + egfp_reporter.fa ──► custom_genome.fa            │
│                                                                                        │
│   [2. GTF 주석 파일 제작]                                                              │
│   gencode.v44.primary_assembly.annotation.gtf                                          │
│   + egfp_annotation.gtf (gene, transcript, exon 계층 구조 엄격 준수)                   │
│   ──► custom_annotation.gtf                                                            │
│                                                                                        │
│   [3. STAR 게놈 인덱스 생성 (STAR --runMode genomeGenerate)]                           │
│   SA Index 계산식: --genomeSAindexNbases = min(14, floor(log2(GenomeLength)/2 - 1))   │
│   인덱스 디렉토리: /ref/STAR_GRCh38_eGFP/ (사전 계산된 Suffix Array 및 Splice Junction)│
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 GTF(Gene Transfer Format) 레코드 작성 규격
단일 세포 정렬기(STARsolo, Cell Ranger)는 GTF 주석의 계층 구조와 속성(Attributes)을 엄격하게 검증합니다. `gene_id`, `transcript_id`, `gene_name`이 누락되거나 철자가 틀리면 카운팅에서 영구 배제됩니다.

```gtf
# eGFP 커스텀 GTF 주석 예시 (Tab-delimited)
eGFP_chr	transgene	gene	1	720	.	+	.	gene_id "eGFP_gene"; gene_name "eGFP"; gene_biotype "protein_coding";
eGFP_chr	transgene	transcript	1	720	.	+	.	gene_id "eGFP_gene"; transcript_id "eGFP_transcript"; gene_name "eGFP"; gene_biotype "protein_coding";
eGFP_chr	transgene	exon	1	720	.	+	.	gene_id "eGFP_gene"; transcript_id "eGFP_transcript"; exon_number "1"; gene_name "eGFP";
```

### 1.2 STAR 인덱싱 수학적 파라미터 도출
STAR의 게놈 인덱스 생성 시 가장 중요한 파라미터 중 하나는 SAindex(접두 룩업 테이블)에 쓰는 k-mer 길이인 `--genomeSAindexNbases`입니다(접미사 배열 자체의 희소도는 `--genomeSAsparseD`):

$$\text{genomeSAindexNbases} = \min\left(14, \, \left\lfloor \frac{\log_2(L_{\text{genome}})}{2} - 1 \right\rfloor \right)$$

- 인간 전장 유전체 ($L_{\text{genome}} \approx 3.1 \times 10^9\text{ bp}$):
  $$\left\lfloor \frac{\log_2(3.1 \times 10^9)}{2} - 1 \right\rfloor = \left\lfloor \frac{31.53}{2} - 1 \right\rfloor = \lfloor 14.76 \rfloor = 14$$
  따라서 기본값 14를 사용합니다.
- 만약 바이러스 게놈이나 미니 유전체(예: $1\text{ Mb}$)를 단독 인덱싱할 경우:
  $$\left\lfloor \frac{\log_2(10^6)}{2} - 1 \right\rfloor = \left\lfloor \frac{19.93}{2} - 1 \right\rfloor = 8$$
  값을 8 내지 9로 하향하지 않으면 STAR가 세그멘테이션 오류(Segmentation fault)로 비정상 종료됩니다.

---

## 2. 단일핵 전사체(snRNA-seq)를 위한 STARsolo 파라미터 최적화

단일세포(scRNA-seq)와 단일핵(snRNA-seq)은 포획되는 전사체 분자의 성상에 근본적인 생물학적 차이가 있습니다.
- 세포질(Cytoplasm)에는 스플라이싱이 완료된 성숙 전사체(Mature mRNA)가 존재합니다.
- 반면 세포핵(Nucleus) 내부에는 스플라이싱이 진행 중인 미성숙 전구체 전사체(Nascent Pre-mRNA)가 전체 RNA의 **40% ~ 50%**를 차지하며, 막대한 양의 미정제 인트론(Intronic sequences)을 포함하고 있습니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               snRNA-seq에서 STARsolo --soloFeatures 옵션별 카운팅 비교                  │
│                                                                                        │
│   [전사체 구조]    Exon 1            Intron 1 (Pre-mRNA)             Exon 2            │
│                   ┌────────┐─────────────────────────────────────┌────────┐           │
│                   │ Mature │            Un-spliced               │ Mature │           │
│                   └────────┘─────────────────────────────────────└────────┘           │
│   포획된 리드  :    Read A                      Read B                     Read C      │
│                     (Exonic)                   (Intronic)                 (Spliced)    │
│                                                                                        │
│   1. --soloFeatures Gene (Standard Exonic):                                            │
│      ├── Read A: 인정 (Count +1)                                                       │
│      ├── Read B: 폐기 (Intronic - Discarded!) ──► snRNA-seq 분자의 40~50% 증발!        │
│      └── Read C: 인정 (Count +1)                                                       │
│                                                                                        │
│   2. --soloFeatures GeneFull (Recommended for snRNA-seq):                              │
│      ├── Read A: 인정 (Count +1)                                                       │
│      ├── Read B: 인정 (Count +1, Pre-mRNA Intron 포획 성공!)                           │
│      └── Read C: 인정 (Count +1) ──► 총 유전자 검출 수 및 UMI 1.8~2.2배 대폭 상승      │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 프로덕션급 STARsolo 실행 스크립트

```bash
#!/usr/bin/env bash
set -euo pipefail

# 화이트리스트는 케미스트리와 일치해야 함: v3/v3.1(UMI 12 bp)은 3M-february-2018.txt.gz를 gunzip한 파일,
# v2(UMI 10 bp)는 737K-august-2016.txt (STARsolo 문서)
STAR \
    --runThreadN 16 \
    --genomeDir /ref/STAR_GRCh38_GENCODE44 \
    --readFilesIn data/R2_cDNA.fastq.gz data/R1_Barcodes.fastq.gz \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /ref/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures GeneFull \
    --soloCellFilter EmptyDrops_CR \
    --soloStrand Forward \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --limitBAMsortRAM 32000000000 \
    --outFileNamePrefix results/STARsolo_snRNA_
```

- `--soloFeatures GeneFull`: 엑손과 인트론을 아우르는 전체 유전자 영역에 매핑된 리드를 단일 분자로 카운팅하여 단일핵 데이터의 UMI 포획률을 극대화합니다.
- `--soloCellFilter EmptyDrops_CR`: CellRanger flavor의 EmptyDrops(기대 세포 수, UMI 순위 45,000~90,000 구간 ambient, umiMin 500, FDR 0.01, simN 10000)를 적용하여 진짜 세포와 빈 액적을 구별합니다('이중 무릎' 알고리즘은 공식 문서에서 확인되지 않음).

---

## 3. 다중 매핑 리드(`NH > 1`) 및 가유전자(Pseudogenes) 해결

참조 게놈에는 높은 서열 상동성을 지닌 가유전자(Pseudogenes)나 중복 유전자 패밀리(Paralogs)가 수천 개 존재합니다.
대표적으로 전립선암의 안드로겐 수용체(*AR*)는 전사 억제성 안티센스 유전자 및 가유전자(*AR-AS1*)와 서열을 일부 공유하며, *GAPDH*, *ACTB* 등도 다수의 가유전자를 지닙니다.

정렬 레코드에서 `NH:i` 태그는 해당 리드가 참조 게놈 내에서 정렬 가능한 위치의 총 개수를 나타냅니다:
- $\text{NH} = 1$: 고유 매핑 (Uniquely mapped read)
- $\text{NH} > 1$: 다중 매핑 (Multi-mapping read)

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        다중 매핑 리드 해결의 세 가지 전략                              │
│                                                                                        │
│   1. Unique Only (기본 필터링)                                                         │
│      (STARsolo 기본) 여러 유전자에 걸치는 다유전자 UMI를 제외. 한 유전자 내 다중       │
│      위치/겹침은 구분되나, 중복 유전자 발현량이 결손(Under-counting)되는 부작용.        │
│                                                                                        │
│   2. EM (Expectation-Maximization) 확률 분배 알고리즘                                  │
│      고유 매핑된 리드들의 유전자별 상대적 발현 밀도를 사전 확률로 사용하여,           │
│      다중 매핑 분자를 사후 확률 가중치에 따라 비례 배분:                              │
│                                                                                        │
│      P(Read r ∈ Gene g) = [ θ_g * P(r | g) ] / ∑_g' [ θ_g' * P(r | g') ]              │
│                                                                                        │
│   3. Rescue 전략 (STARsolo --soloMultiMappers Rescue)                                  │
│      다유전자 UMI를 '고유 UMI + 균등 분배 UMI'의 합에 비례해 분배(EM 1단계).          │
│      (UMI를 공유하는 고유 리드로 구제하는 방식에 가까운 것은 PropUnique)               │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 4. 단일 세포 정렬 및 BAM 제어 10대 실패 모드 (Top 10 Failure Modes)

| 순번 | 실패 모드 (Failure Mode) | 근본 원인 (Root Cause) | 임상 및 분석적 파급 효과 | 감지 신호 및 해결 솔루션 |
|---|---|---|---|---|
| **01** | **R1 / R2 FASTQ 파일 역전** | 시퀀서 demux 시 Read 1(바코드 28bp)과 Read 2(cDNA 90bp) 파일명이 뒤바뀜 | 정렬 매핑률 0.5% 미만으로 추락, 전체 실험 실패 | 파일 사이즈 확인 (cDNA가 훨씬 큼), STARsolo 입력 순서 교정 |
| **02** | **BAM 정렬 순서 불일치** | `samtools index` 실행 시 `BAI index not found / unsorted` 에러 | 게놈 슬라이싱 불가, igv 및 pysam fetch 즉시 crash | `samtools sort -@ 8 -o sorted.bam input.bam` 사전 수행 |
| **03** | **BAM 좌표 정렬 시 OOM** | `samtools sort -m 4G` 지정 시 스레드당 4GB 할당 ($4\times 16 = 64\text{GB}$) | OS 커널 OOM Killer에 의해 프로세스 SIGKILL 강제 종료 | `samtools sort -m`은 스레드당 메모리임을 명심! `-m 1G`로 축소 |
| **04** | **좌표계 0 vs 1 오프-바이-원** | IGV(1-based) 좌표를 `pysam.fetch()`(0-based half-open)에 그대로 대입 | 정션 및 단일염기변이 핫스팟의 정확한 염기 1bp 누락 | `pysam` 인자 전달 시 $start_{\text{pysam}} = pos_{\text{1-based}} - 1$ 변환 공식 적용 |
| **05** | **CIGAR 'N' 미고려 파싱 버그** | RNA-seq 리드의 게놈 도달 범위를 단순 서열 길이(`len(query)`)로 계산 | 스플라이스 인트론 건너뛰기 무시로 종료 좌표 수천 bp 오차 | 반드시 CIGAR $M, D, N, =, X$의 합으로 게놈 범위 계산 |
| **06** | **pysam 파일 핸들러 누수** | 루프 내에서 `AlignmentFile`을 열고 `close()` 누락 | `OSError: Too many open files` 시스템 충돌 발생 | 반드시 `with pysam.AlignmentFile(...) as bam:` 컨텍스트 매니저 사용 |
| **07** | **UMI 충돌 미보정** | 세포당 발현량이 매우 높은 유전자에서 서로 다른 분자가 우연히 동일 UMI 획득 | 서로 다른 분자가 하나로 합쳐져 과소계수, 상한선 포화 및 카운트 왜곡 | UMI 길이 확대(키트 선택) 또는 포화 보정(stochastic labeling 모형, Fu 2011)으로 보정. (해밍 거리 1 병합은 PCR/시퀀싱 *오류*로 인한 과대계수용이며 충돌은 오히려 악화시킨다) |
| **08** | **참조 유전체 빌드·접두사 불일치** | 빌드(GRCh37/38)와 `chr` 접두사 규약은 별개 문제(Ensembl GRCh38도 접두사 없음). 예: BAM은 `chr` 접두사, GTF는 접두사 없음 혼용 | 모든 리드가 Unassigned_NoFeatures로 분류되어 카운트 0 | BAM 헤더 `@SQ`와 GTF 염색체 명명 규칙(`chr1` vs `1`) 일치성 검증 |
| **09** | **저품질 매핑(`MAPQ=0`) 미필터** | 가유전자 및 반복 서열에 정렬된 저품질 리드가 SNV/Fusion 분석에 유입 | 위양성 돌연변이 및 융합 유전자 대량 호출 | `pysam` 분석 시 반드시 `read.mapping_quality >= 30` 필터 적용 |
| **10** | **Ambient RNA에 의한 가짜 클론** | 용해된 세포의 유리 mRNA가 빈 액적에 포획되어 가짜 종양 서브클론 형성 | 정상 세포가 *AR-V7*이나 *TMPRSS2-ERG* 양성으로 오인 | SoupX/CellBender/DecontX 등 ambient RNA 제거 및 정상 참조 기질세포 발현 대조 (EmptyDrops_CR은 빈 액적 구별일 뿐 세포 내 오염 제거가 아님) |

---

## 5. 5단계 단일 세포 시퀀스 정렬 캡스톤 프로젝트 평가 루브릭

단일핵 전사체 업스트림 정렬 및 변이 마이닝 파이프라인의 완성도를 평가하기 위한 대학원 수준의 5단계 종합 루브릭입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               5-Stage Single-Cell Alignment Capstone Evaluation Rubric                 │
│                                                                                        │
│   Stage 1: 원천 QC 및 맞춤형 레퍼런스 유전체 구축 (Raw QC & Custom Reference)           │
│   ├── FastQC/Fastp를 통한 R1(바코드 Q30 > 85%) 및 R2(cDNA 어댑터 트리밍) 무결성 검증     │
│   └── 외인성 유전자(eGFP) 및 전이성 융합 주석을 포함하는 GTF 무결성 검증 (FASTA/GTF Sync) │
│                                                                                        │
│   Stage 2: 고성능 스플라이싱 핵 정렬 (Spliced Nuclear Alignment via STARsolo)          │
│   ├── --soloFeatures GeneFull 적용으로 Pre-mRNA 인트론 리드 포획 효율 극대화           │
│   └── EmptyDrops_CR 적용으로 Ambient RNA 배경 잡음 제거 및 유효 세포 콜링             │
│                                                                                        │
│   Stage 3: 게놈 좌표 슬라이싱 및 SAM 태그 무결성 검증 (BAM Slicing & Tags)             │
│   ├── samtools view -F 3844 기반의 고속 유전자좌 슬라이싱 (AR, TMPRSS2, ERG)           │
│   └── CB(Cell Barcode), UB(UMI) 보조 태그의 정합성 검증 (xf는 Cell Ranger BAM 전용)   │
│                                                                                        │
│   Stage 4: 단일 세포 이소폼 및 융합 마이닝 (Single-Cell Isoform & Fusion Mining)       │
│   ├── CIGAR 'N' 연산자 기반 AR-V7 스플라이싱 접합부 디코딩 및 단일 세포 PSI(Ψ) 산출     │
│   └── SA 보충 정렬 태그 및 CIGAR 'S' 기반 TMPRSS2-ERG 키메라 분할 리드 확인 (3′ 라이브러리는 한계)│
│                                                                                        │
│   Stage 5: 멀티오믹스 클론 계통 결합 및 임상 리포팅 (Multi-Omic Integration)           │
│   ├── 단일 세포 대규모 CNA(Chr 8p 등)와 TMPRSS2-ERG 융합의 일치성 검증 (Edel은 DNA 근거)  │
│   └── Enzalutamide 저항성 AR-V7 고발현 서브클론의 단일 세포 진화 계통수 도출           │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

| 평가 단계 | 배점 | Pass 기준 (상위 5% 프로덕션 기준) | Fail 기준 (재작업 필요) |
|---|---|---|---|
| **Stage 1** | 20점 | R1/R2 구조 정의 정확, 커스텀 GTF 주석 문법 오류 전무, STAR SA 인덱스 정상 생성 | R1/R2 파일 뒤바뀜, GTF 내 `gene_id` 누락으로 정렬기 충돌 |
| **Stage 2** | 20점 | `GeneFull` 적용으로 총 UMI 포획률 40% 이상 향상, 고유 매핑률 > 75% 달성 | 엑손만 카운팅(`Gene`)하여 인트론 리드 50% 유실, OOM 충돌 |
| **Stage 3** | 20점 | `samtools` 비트 플래그 필터링(-F 3844) 적용, 64비트 가상 오프셋 슬라이싱 밀리초 완료 | 인덱스 미생성으로 선형 탐색, 보충 정렬 오제거로 융합 유실 |
| **Stage 4** | 20점 | `pysam` CIGAR 파싱으로 *AR-V7* PSI 정확 계산, 융합 분할 리드 UMI 중복 제거 수행·검증(3′ 라이브러리 한계 명시) | Off-by-one 버그로 정션 좌표 불일치, PCR 중복 미제거 위양성 |
| **Stage 5** | 20점 | 융합 양성 세포와 대규모 CNA(8p 등) 클론의 일치성 검토(21q22 국소 결실은 DNA 근거 필요), 임상 약제 내성 기전 규명 | CNA 클론과 돌연변이 상태 간 불일치 원인 미진단, 단순 나열 |

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [수학적 원리 규명 (Mathematical Rigor)]
단일핵 시퀀싱(snRNA-seq) 데이터에서 동일한 FASTQ 파일을 대상으로 STARsolo 정렬을 수행하였다:
- Run A (`--soloFeatures Gene`): 총 유효 검출 UMI 수 $N_A = 12,000,000$, 세포당 중앙 UMI 수 $= 2,400$.
- Run B (`--soloFeatures GeneFull`): 총 유효 검출 UMI 수 $N_B = 21,600,000$, 세포당 중앙 UMI 수 $= 4,320$.
1) Run B에서 검출된 전체 전사체 분자 중 핵 내 미성숙 전구체 전사체(Pre-mRNA) 유래 인트론 리드가 차지하는 분자 비율(Fraction of Intronic Molecules)을 수학적으로 계산하시오.
2) 만약 분석관이 snRNA-seq 데이터에 Run A 옵션을 적용하여 하류 클러스터링 및 차별발현유전자(DEG) 분석을 진행할 경우, 검출 한계(Limit of Detection) 및 세포 필터링(Cell Calling) 단계에서 발생할 수 있는 치명적 문제점을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 인트론 유래 분자 비율 계산**:
  - Run A는 순수한 성숙 엑손(Exonic) 영역에만 매핑된 분자를 카운팅합니다:
    $$N_{\text{exon}} = N_A = 12,000,000$$
  - Run B는 성숙 엑손과 미성숙 인트론(Intronic) 영역에 매핑된 모든 유전자 전장 분자를 카운팅합니다:
    $$N_{\text{total}} = N_B = 21,600,000$$
  - 따라서 순수한 인트론 영역에 매핑되어 Run B에서 새롭게 구제(Rescued)된 분자 수는:
    $$N_{\text{intron}} = N_B - N_A = 21,600,000 - 12,000,000 = 9,600,000$$
  - Run B 전체 분자 중 인트론 분자의 비율:
    $$\text{Fraction}_{\text{intron}} = \frac{N_{\text{intron}}}{N_{\text{total}}} = \frac{9,600,000}{21,600,000} = \frac{4}{9} \approx 0.444 \ (44.4\%)$$
  - 즉, `GeneFull`로 추가 포획된 분자가 전체의 $44.4\%$입니다. 이는 '핵 내 전사체의 44.4%가 스플라이싱 전'이라는 뜻이 아니며(oligo-dT의 인트론 A-rich 내부 프라이밍, 중복 유전자 겹침으로 인한 다유전자 UMI 제외 차이 등도 영향), 4/9 산술만 정확합니다.

- **2. Run A 적용 시의 치명적 문제점**:
  - **세포 손실(False Cell Dropout)**: `EmptyDrops`와 같은 세포 콜링 알고리즘은 액적 내의 총 UMI 분자 수를 기반으로 빈 액적(Ambient noise)과 진짜 세포를 구분합니다. 전체 RNA 분자의 $44.4\%$가 버려지면, 발현량이 상대적으로 낮은 희귀 세포군(예: 줄기세포, 비활성 림프구)은 총 UMI 수가 필터링 임계값 이하로 추락하여 세포가 아닌 빈 액적으로 오판되어 수천 개의 세포가 영구 소실됩니다.
  - **유전자 검출 민감도 급감**: 전사인자나 신호전달 수용체와 같이 기저 발현량이 낮은 핵심 마커 유전자들은 엑손 리드가 단 1개도 잡히지 않아 발현량이 0으로 기록되는 유전자 드롭아웃(Gene Dropout)이 극대화되어, 클러스터링 해상도가 심각하게 왜곡되고 주요 바이오마커가 발굴되지 않습니다.

</details>

---

### Q2. [알고리즘 구현 및 시스템 장애 진단 (System Failure Diagnosis)]
연구팀이 고성능 HPC 서버(64코어 CPU, 256GB RAM)에서 10x 3' snRNA-seq FASTQ 데이터(총 4억 리드)를 STARsolo로 정렬하던 중 다음과 같은 두 가지 심각한 시스템 장애가 발생하였다:
- Issue 1: 정렬은 정상 종료되었으나, 생성된 미정렬 BAM 파일을 좌표 정렬하기 위해 `samtools sort -@ 32 -m 8G -o sorted.bam Aligned.out.bam`을 실행한 직후, 시스템이 심각하게 느려지며 `Killed` 메시지와 함께 프로세스가 강제 종료됨.
- Issue 2: 메모리 문제를 해결하고 생성한 최종 카운트 매트릭스에서 모든 세포의 유전자 발현량이 정확히 0으로 출력됨. STAR의 `Log.final.out`을 확인한 결과 `Uniquely mapped reads %: 82%`, `Number of reads mapped to multiple loci: 10%`로 정렬은 매우 훌륭하였음.
1) Issue 1에서 256GB라는 거대한 물리 RAM을 보유한 서버에서 OOM Killer가 발동한 시스템적 원인을 비트/메모리 할당 공식 관점에서 규명하시오.
2) Issue 2에서 정렬 성공률이 82%에 달함에도 불구하고 세포별 유전자 발현 카운트가 100% 0으로 소실된 원인을 GTF 속성 및 바코드 화이트리스트 관점에서 진단하고 즉각적인 조치 방안을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. Issue 1: `samtools sort` 메모리 폭발 원인 분석**:
  - `samtools sort` 명령에서 `-m` 옵션은 프로세스 전체의 총 메모리 상한선이 아니라, **스레드(Thread) 1개당 할당되는 최대 메모리 버퍼**를 지정하는 옵션입니다.
  - 사용자가 지정한 스레드 수는 `-@ 32` (32개)이고, 스레드당 메모리는 `-m 8G` (8GB)입니다.
  - 따라서 `samtools`가 동적 할당을 요청한 총 가상 메모리 크기는:
    $$\text{Total RAM} = 32 \times 8\text{ GB} = 256\text{ GB}$$
  - 여기에 OS 커널 메모리, 기본 버퍼 캐시, Python/Bash 프로세스 오버헤드가 더해지면서 물리 RAM 용량인 256GB를 즉시 초과하였고, 시스템 스왑(Swap) 공간이 고갈되면서 리눅스 커널의 Out-Of-Memory (OOM) Killer가 `samtools` 프로세스에 `SIGKILL` (Signal 9)을 전송하여 강제 종료시킨 것입니다.
  - **해결책**: 스레드당 메모리를 `-m 2G` ($32 \times 2 = 64\text{ GB}$) 또는 `-m 1.5G`로 하향 조정하여 실행해야 합니다.

- **2. Issue 2: 카운트 100% 소실의 원인 진단 및 조치 방안**:
  - **원인 1 (1순위): 화이트리스트/화학(chemistry) 불일치**:
    - 고유 매핑률 82%는 cDNA(Read 2)가 게놈에 잘 정렬되었다는 뜻일 뿐, 세포 바코드(Read 1) 판독과는 무관하다. v3 데이터에 v2 화이트리스트(`737K-august-2016.txt`)를 쓰거나 `--soloUMIlen`/`--soloCBlen`이 맞지 않으면 바코드가 거의 일치하지 않아 모든 리드가 세포에 배정되지 않고 카운트가 0이 된다. `Solo.out/Barcodes.stats`의 일치율로 확인한다.
  - **원인 2: feature 설정/주석 불일치**:
    - `--soloFeatures`(Gene/GeneFull 등)와 GTF 속성(`gene_id`), `--soloStrand` 설정이 데이터와 맞지 않으면 리드가 `Unassigned_NoFeature`로 분류된다. 염색체 접두사(`chr` 유무) 불일치는 보통 `genomeGenerate` 단계에서 실패하거나 겹침이 0이 되어 나타난다.
    - (참고: FASTQ R1/R2 인자 순서 역전은 매핑률이 극단적으로 낮아지는 증상(01장 실패 모드)이며, 고유 매핑 82%인 이 사례의 원인이 될 수 없다.)
  - **즉각 조치 방안**: `Solo.out/Barcodes.stats` 파일의 바코드 일치율을 먼저 확인하고(화이트리스트·chemistry 설정 교정), 이어서 GTF 파일의 속성 헤더를 `grep -v '^#' gencode.gtf | cut -f9 | head -n 5`로 조회하여 STARsolo의 주석 매핑 설정을 일치시킵니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
전립선암 전이 환자의 골전이(Bone metastasis) 조직 검체에서 단일핵 전사체(snRNA-seq) 시퀀싱을 진행하였다.
정렬 및 변이 검출 파이프라인 결과, 분석관은 "이 환자는 전이성 종양 세포 전체에서 *AR* 유전자가 10배 이상 증폭되어 있고, 모든 종양 세포에서 *TMPRSS2-ERG* 유전자 융합과 *AR-V7* 변이체가 100% 양성이므로 초고위험 내성 암으로 분류된다"고 보고하였다.
그러나 데이터 무결성을 검증하던 수석 엔지니어는 다음 세 가지 데이터 이상 징후를 포착하였다:
1. 종양 세포뿐만 아니라 검체 내에 존재하는 정상 내피세포(Endothelial cells, *PECAM1*+)와 정상 대식세포(Macrophages, *CD68*+)에서도 *AR-V7* 정션 리드와 *TMPRSS2-ERG* 융합 리드가 각각 1~2개씩 동일하게 검출됨.
2. 세포별 총 UMI 분포 곡선(Barcode-Rank Plot)에서 배경 액적(Background droplets)의 RNA 농도가 일반 검체 대비 5배 이상 비정상적으로 높음.
3. 골 조직 탈칼슘화(Decalcification) 및 단핵 분리 과정에서 세포 용해(Cell lysis) 비율이 45%로 심각하게 보고됨.

수석 생물정보학 엔지니어의 관점에서:
1) 분석관의 "모든 세포에서 융합 및 변이체 양성"이라는 임상 보고가 어떤 기술적 오류(Artifact)에 기반한 오진인지 그 분자 메커니즘을 규명하시오.
2) 이 데이터셋에서 진짜 종양 세포 고유의 유전자형을 복원하고 위양성 오염을 제거하기 위해 적용해야 할 3단계 생물정보학적 정화 파이프라인(Decontamination Pipeline)을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 임상 보고의 기술적 오류 및 분자 메커니즘 규명**:
  - **판정: 심각한 주변 유리 전사체(Ambient RNA / Soup RNA) 오염에 의한 위양성 가짜 양성(False-positive Pseudoclonality)**입니다.
  - **분자 메커니즘**:
    - 골전이 조직 처리 중 심각한 세포 용해(45% lysis)가 발생하면서, 파괴된 종양 세포핵 내부의 고농도 전사체(초과발현된 *AR*, *AR-V7*, *TMPRSS2-ERG* mRNA)가 세포외 용액(Cell suspension buffer)으로 대량 유출되었습니다.
    - 미세유체 액적(Microfluidic droplet) 형성 시, 각 액적에는 표적 단일 세포(정상 대식세포나 내피세포) 외에도 이 유리 RNA(Ambient RNA)가 함께 무작위로 캡슐화되었습니다.
    - 결과적으로 정상 세포의 비드 바코드(`CB`)가 붙은 상태에서 주변의 유리 종양 전사체가 함께 역전사되어, 실제로는 종양 변이를 전혀 가지고 있지 않은 정상 기질/면역 세포마저 종양 유전자를 발현하는 것처럼 오인된 것입니다.
    - 분석관이 UMI 임계값과 세포별 특이성을 검증하지 않고 1개라도 리드가 잡히면 양성으로 처리하여 전체 세포가 암세포로 둔갑하는 치명적 오진을 낳았습니다.

- **2. 3단계 생물정보학적 정화 파이프라인 (Decontamination Workflow)**:
  - **단계 1: 액적 주변부 RNA 프로파일 추정 (Ambient RNA Profiling via SoupX / CellBender)**:
    - 세포가 들어있지 않은 명백한 빈 액적군(Empty droplets, UMI < 100)의 BAM 리드를 모아 'Ambient RNA Soup'의 전사체 조성을 정량적으로 모델링합니다.
    - *TMPRSS2-ERG* 및 *AR-V7* 분자가 배경 잡음 용액에서 차지하는 기여도 비율(Contamination fraction $\rho$)을 베이지안 확률로 추정합니다.
  - **단계 2: 세포별 통계적 유의성 검정 (Statistical Significance Testing)**:
    - 특정 단일 세포 $c$에서 검출된 융합/변이 UMI 수가 단순히 배경 용액의 푸아송 무작위 유입(Poisson sampling)에 의해 우연히 설명될 수 있는 확률을 계산합니다:
      $$P(X \ge k_c \mid \lambda = N_c \cdot \rho_{\text{fusion}}) < 0.01$$
    - 배경 오염 기대값보다 통계적으로 유의하게 높은 독립 UMI 분자($N_{\text{UMI}} \ge 3$)를 보유한 세포만을 진성 양성 세포로 인정합니다.
  - **단계 3: 정상 세포 음성 대조군 기반 임계값 보정 및 CNA 교차 검증**:
    - 정상 내피세포 및 대식세포에서 관찰되는 최대 UMI 수(1~2개)를 노이즈 컷오프(Baseline Noise Floor)로 설정하여 그 이하의 신호를 일괄 소거합니다.
    - 최종적으로 융합 양성 세포가 InferCNV 상에서 실제로 종양 고유의 대규모 카피수 변이(예: Chr 8p 결실)를 보유하고 있는지 매핑하여(21q22 국소 결실은 inferCNV 해상도 밖이므로 DNA 수준 근거로 확인), 종양 게놈을 가진 세포에서만 진성 변이체로 확정합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Fu GK, Hu J, Wang PH, et al. Counting individual DNA molecules by the stochastic attachment of diverse labels. *Proc Natl Acad Sci U S A*. 2011;108(22):9026-31. PMID: 21562209. doi: 10.1073/pnas.1017621108
### 검수 근거 문헌
1. Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Res*. 2017;27(3):491-499. PMID: 28100584. doi: 10.1101/gr.209601.116
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- STAR 2.7.11 (github.com/alexdobin/STAR). 원 논문: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15-21. PMID: 23104886. doi: 10.1093/bioinformatics/bts635
<!-- REFERENCES:END -->

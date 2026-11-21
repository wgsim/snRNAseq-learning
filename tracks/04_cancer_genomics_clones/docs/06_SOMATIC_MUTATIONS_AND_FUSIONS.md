# Chapter 06. 체세포 단일염기변이(SNV) 및 융합 전사체 검출 (Somatic Mutations & Gene Fusions in snRNA-seq)

> **Cognitive Goal**: 단일핵 전사체(snRNA-seq)의 정렬 데이터(BAM)로부터 세포 바코드(`CB` tag)와 고유 분자 식별자(`UB`/UMI tag)를 추출하여 단일 세포 해상도의 체세포 단일염기변이(Somatic Single Nucleotide Variants, SNVs) 및 키메라 융합 전사체(Gene Fusions)를 검출하는 생물정보학적 원리를 체득한다. 전립선암의 대표적 드라이버 구조 변이인 *TMPRSS2-ERG* 접합부 리드(Split reads 및 Discordant read pairs)를 식별하고 위양성(Artifacts)을 제거하는 판독 알고리즘을 이해한다. 나아가 주요 발암 드라이버 돌연변이(*SPOP* p.F133V, *FOXA1*, *TP53* R273H, *AR* LBD 변이)의 단일 세포 대립유전자 변이 빈도(VAF)를 산출하고, 전사 버스팅 및 모노알렐릭 발현을 보정하는 베이지안 유전자형 모델을 도출한다. 최종적으로 체세포 돌연변이와 대규모 카피수 변이(CNA) 서브클론을 결합하는 공클러스터링(Co-clustering) 프레임워크를 정립한다.

---

## 1. 단일핵 전사체 BAM 파일의 구조와 세포 바코드 기반 판독 (BAM Structure & CB/UB Extraction)

단일 세포 전사체 데이터(10x Genomics Chromium, DroNc-seq 등)의 표준 정렬 파일은 게놈 좌표뿐만 아니라 각 리드가 유래한 개별 세포와 원본 RNA 분자를 특정할 수 있는 보조 SAM 태그(Optional SAM tags)를 포함합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│              단일핵 전사체 정렬 리드(SAM/BAM)의 핵심 태그 구조                         │
│                                                                                        │
│   QNAME: A00123:456:H7K3FDSX2:1:1101:1000:2048                                         │
│   FLAG : 99 (Read paired, Read mapped in proper pair, Mate reverse strand)             │
│   RNAME: chr21                POS  : 41498119             MAPQ: 60 (예시; STAR 고유 정렬은 255)│
│   CIGAR: 54M76S (54bp Match, 76bp Soft-clipped)                                        │
│   TAGS :                                                                               │
│     CB:Z:AAACCTGAGTCGTTTG-1   <── [Cell Barcode: 단일 세포 식별자]                    │
│     UB:Z:TAGGCTAACG            <── [UMI: 단일 RNA 전사 분자 식별자]                     │
│     GN:Z:ERG                   <── [Gene Name: 정렬된 유전자]                           │
│     SA:Z:chr21,41498120,+,54M,60,0;  <── [Supplementary Alignment: 키메라 분할 접합]  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 핵심 메타데이터 태그
1. **`CB` (Cell Barcode)**: 각 비드(Bead) 또는 액적(Droplet)에 고유하게 부여된 세포 식별 바코드입니다. 셀 벤더 파이프라인(Cell Ranger 등)에서 에러 보정(Error-corrected)된 바코드 뒤에 하이픈과 젬 웰(Gem group) 번호(예: `-1`)가 붙습니다.
2. **`UB` (Unique Molecular Identifier, UMI)**: PCR 증폭 전 단일 RNA 분자에 결찰(Ligation)된 무작위 염기서열(10~12bp)입니다. 동일한 `CB`와 `UB`를 공유하는 복수의 리드는 단 1개의 원본 mRNA 분자로부터 증폭된 중복 리드(PCR duplicates)이므로 1회로 병합(Deduplication)되어야 합니다.
3. **`NH` (Number of reported alignments)**: 해당 리드가 게놈 상에 매핑된 위치의 수입니다. $\text{NH} = 1$인 고유 정렬(Uniquely mapped reads)만을 사용하여 위양성 변이 검출을 억제합니다.

---

## 2. *TMPRSS2-ERG* 키메라 융합 전사체 검출 알고리즘 (Chimeric Junction & Split-Read Analysis)

전립선암 환자의 50% 이상에서 관찰되는 *TMPRSS2-ERG* 유전자 융합은 염색체 21q22 영역의 약 3Mb 간질성 결실(Interstitial deletion) 또는 상호 전좌(Translocation)에 의해 발생합니다. 안드로겐 반응성 프로모터를 지닌 *TMPRSS2* (엑손 1 또는 2)의 비번역 영역이 ETS 전사인자 계열인 *ERG* (엑손 4 또는 2)의 단백질 코딩 서열과 연결되어 안드로겐에 의해 ERG 종양 단백질이 대량으로 전사 유도됩니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   TMPRSS2-ERG 키메라 접합부 리드 검출의 두 가지 기하학                 │
│                                                                                        │
│   [1. Split-Read (분할 리드)]                                                          │
│   TMPRSS2 Exon 1              ERG Exon 4                                               │
│   ───────────────┬───────────────────────────────────                                  │
│   ... G G C A G C│A T T C A G G A C T C A G T C C ...                                  │
│   ▲              ▲                  ▲                                                  │
│   └──────────────┼──────────────────┘                                                  │
│      Read: [ 35bp Match | 65bp Match ]  ==> CIGAR: 35M65S 및 SA 태그 생성               │
│                                                                                        │
│   [2. Discordant Read Pair (불일치 리드 쌍)]                                           │
│   R1 ──► [TMPRSS2 정렬]                  R2 ◄── [ERG 정렬]                             │
│   (정상 거리인 300bp를 초과하여 약 3Mb 떨어진 위치 또는 이종 염색체 간 불일치 정렬)       │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 키메라 접합부(Chimeric Junction) 판독의 두 가지 축
> **적용 범위 주의**: 아래 두 축은 paired-end 전장 전사체/DNA 시퀀싱에 일반적인 융합 검출 기하학입니다. 10x 3' snRNA-seq에서는 R1이 바코드/UMI이고 R2가 poly(A) 부근 cDNA의 단일 말단 읽기이므로 불일치 리드 쌍(discordant pair)은 적용할 수 없고, 5' 쪽에 있는 *TMPRSS2-ERG* 접합부를 가로지르는 분할 리드도 원리적으로 거의 포획되지 않습니다(*ERG* 발현 증가를 대리 지표로 사용; 접합부 검출은 5' 라이브러리·롱리드·타깃 증폭 한정).

1. **분할 리드(Split-Reads)**:
   - 단일 시퀀싱 리드($100 \sim 150\text{bp}$) 자체가 정확히 융합 접합부(Fusion breakpoint)를 가로질러 형성된 경우입니다.
   - 리드의 전반부는 *TMPRSS2*의 엑손 말단에 정렬되고, 후반부는 *ERG*의 엑손 시작점에 정렬됩니다. 정렬기(STAR, Minimap2)는 이를 주 정렬(Primary alignment)과 보충 정렬(Supplementary alignment, `SA` tag)로 분할 보고하며, CIGAR 문자열에 소프트 클리핑(`S`)이 발생합니다.
2. **불일치 리드 쌍(Discordant Read Pairs)**:
   - Paired-end 라이브러리에서 Read 1은 *TMPRSS2*에 매핑되고, 메이트인 Read 2는 수 메가베이스(Mb) 떨어진 *ERG*에 매핑되는 경우입니다. 정상적인 cDNA 삽입 크기(Insert size, $200 \sim 500\text{bp}$)를 극단적으로 벗어나는 이상 정렬을 통해 융합을 포착합니다.

### 2.2 단일 세포 융합 유전자 판독의 위양성(Artifact) 필터링
단일 세포 전사체 역전사(RT) 효소 반응 중 발생하는 가짜 키메라(Artifactual chimeras)를 제거하기 위해 다음과 같은 엄격한 3단계 필터링 기준을 수립합니다:
- **RT 템플릿 스위칭(Template Switching) 방지 (4bp·UMI≥2 등 아래 기준은 임의 예시)**: 접합부 양측에 4bp 이상의 미세 상동성(Microhomology) 서열이 존재하여 역전사 효소가 비특이적으로 이동하면서 합성된 인공 산물인지 검증.
- **최소 독립 UMI 수**: 동일한 `CB` 내에서 최소 2개 이상의 독립된 `UB`가 동일 접합부를 지지해야 단일 세포 융합 양성($\text{Fusion}^+$)으로 판정.
- **정상 상피 및 기질 세포 대조**: 정상 전립선 상피세포, 내피세포, 섬유아세포 등 비종양 참조 세포군에서 해당 접합부 리드가 전혀 검출되지 않아야 함.

---

## 3. 발현된 체세포 단일염기변이(SNV) 호출 및 VAF 통계 모델링

전립선암의 주요 서브클론 드라이버 변이는 특정 코돈의 점돌연변이 형태로 발생합니다.
- **_SPOP_ (p.F133V, p.W131G)**: E3 유비퀴틴 리가아제의 MATH 기질 인식 도메인 변이로, AR 복합체 및 ERG 단백질의 분해를 억제하여 종양 형성을 유도함 (*ERG* 융합과 상호 배타적).
- **_FOXA1_ (p.M253K, Wing2 도메인 변이)**: 안드로겐 수용체의 개척자(Pioneer) 인자로서 염색질을 강제로 개방.
- **_TP53_ (p.R273H, p.R248W)**: DNA 결합 도메인의 핫스팟 변이로 p53 종양 억제 기능 상실.
- **_AR_ LBD 변이**: 엔잘루타마이드를 길항제에서 작동제로 바꾸는 대표 변이는 *F877L*(구 표기 F876L)이며, *T878A*는 프로게스테론·일부 항안드로겐에 의해, *L702H*(엑손 4)는 글루코코르티코이드(프레드니손)에 의해 활성화되는 변이로 구분된다.

### 3.1 단일 세포 대립유전자 변이 빈도 (Single-Cell VAF)
특정 단일 세포 $c$에서 유전자좌 $j$에 매핑된 참조 대립유전자(Reference) UMI 카운트를 $k_{c, j}^{\text{ref}}$, 대립유전자 변이(Alternative/Mutant) UMI 카운트를 $k_{c, j}^{\text{alt}}$라 정의합니다. 총 검출 분자 수 $N_{c, j} = k_{c, j}^{\text{ref}} + k_{c, j}^{\text{alt}}$일 때, 관측된 VAF는 다음과 같습니다:

$$VAF_{c, j} = \frac{k_{c, j}^{\text{alt}}}{k_{c, j}^{\text{ref}} + k_{c, j}^{\text{alt}}}$$

### 3.2 전사체 기반 변이 호출의 한계: 전사 버스팅과 알렐릭 드롭아웃
단일 세포 전사체 데이터는 DNA 시퀀싱과 달리 다음과 같은 심각한 생물학적/기술적 교란 요인을 갖습니다:
1. **단일 대립유전자 전사 버스팅 (Monoallelic Transcriptional Bursting)**:
   - 암세포가 이형접합체(Heterozygous, Genotype $0/1$) 돌연변이를 지니고 있더라도, 전사가 발생하는 짧은 순간 동안 돌연변이 대립유전자 또는 정상 대립유전자 중 단 한 가닥만 전사 버스팅을 일으킬 수 있습니다.
   - 결과적으로 실제 이형접합체 세포임에도 불구하고 $VAF = 1$ 또는 $VAF = 0$으로 극단적인 양극화 현상이 발생합니다.
2. **시퀀싱 깊이 부족 및 드롭아웃 (Low Depth & Dropout)**:
   - 많은 단일 세포에서 $N_{c, j} = 0$이거나 $N_{c, j} \le 2$로 매우 낮아 단순 빈도주의적 접근으로는 유전자형을 확정할 수 없습니다.

### 3.3 베이지안 베타-이항 혼합 모델 (Beta-Binomial Mixture Model)
단일 세포의 실제 기저 유전자형 $G_c \in \{0, 1, 2\}$ (0: Ref 동형접합 $0/0$, 1: 이형접합 $0/1$, 2: Alt 동형접합 $1/1$)에 대한 사후 확률을 추정하기 위해 베타-이항(Beta-Binomial) 모델을 구축합니다.

관측된 변이 분자 수 $k_c^{\text{alt}}$는 주어진 총 분자 수 $N_c$와 유전자형 $g$ 하에서 대립유전자 발현 비율 $\theta \sim \text{Beta}(\alpha_g, \beta_g)$를 따르는 베타-이항 분포를 따릅니다:

$$P(k_c^{\text{alt}} \mid N_c, G_c = g) = \binom{N_c}{k_c^{\text{alt}}} \frac{\text{B}(k_c^{\text{alt}} + \alpha_g, \, N_c - k_c^{\text{alt}} + \beta_g)}{\text{B}(\alpha_g, \beta_g)}$$

여기서 $\text{B}(\cdot, \cdot)$는 베타 함수(Beta function)입니다. 유전자형별 초매개변수(Hyperparameters)는 시퀀싱 에러율 $\epsilon \approx 0.005$와 오버디스퍼전(Overdispersion) $\gamma$를 반영하여 설정됩니다:
- $g = 0$ (Homozygous Ref): 평균 대립유전자 비율 $\mu_0 = \epsilon$
- $g = 1$ (Heterozygous Mut): 평균 대립유전자 비율 $\mu_1 = 0.5$ (CNA가 2인 이배체 기준)
- $g = 2$ (Homozygous Alt / LOH): 평균 대립유전자 비율 $\mu_2 = 1 - \epsilon$

베이즈 정리에 의해 단일 세포 $c$의 유전자형 사후 확률(Posterior Probability)은 다음과 같이 계산됩니다:

$$P(G_c = g \mid k_c^{\text{alt}}, N_c) = \frac{P(k_c^{\text{alt}} \mid N_c, G_c = g) \cdot P(G_c = g)}{\sum_{g' \in \{0, 1, 2\}} P(k_c^{\text{alt}} \mid N_c, G_c = g') \cdot P(G_c = g')}$$

사후 확률 $P(G_c = 1 \text{ or } 2 \mid \text{Data}) \ge 0.85$인 세포만을 해당 돌연변이 양성 세포($\text{Mutant}^+$)로 확정 분류합니다.

---

## 4. 체세포 변이(SNV/Fusion)와 카피수 변이(CNA)의 결합 공클러스터링 (Joint Co-clustering)

CNA 프로파일은 게놈 전반의 광역적 복제수 변화를 통해 세포들을 거시적 클론(Macro-clones)으로 분리하는 데 우수하지만, 점돌연변이나 융합 유전자와 같은 국소 분자 사건을 감지하지 못합니다. 반대로 SNV/Fusion 분석은 정밀한 드라이버 사건을 포착하지만 드롭아웃으로 인해 많은 세포에서 유전자형이 결측됩니다.

따라서 두 정보를 상호 보완하는 **이중 앵커(Dual-Anchor) 통합 계통 분석**을 수행합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│             CNA 거시 클론과 체세포 변이(SNV/Fusion) 결합 공클러스터링 매트릭스        │
│                                                                                        │
│   [CNA Macro-Clones]          [Somatic Drivers]             [통합 고해상도 서브클론]    │
│   Macro-Clone A ─────────────► TMPRSS2-ERG (+), PTEN loss ─► Truncal Ancestral Subclone │
│     │                                                                                  │
│     ├─► Sub-divergence 1 ────► FOXA1 (p.M253K) ────────────► Clonal Lineage 1          │
│     │                                                                                  │
│     └─► Sub-divergence 2 ────► TP53 (p.R273H) ─────────────► Pre-plastic Subclone      │
│           │                                                                            │
│           └─► + RB1 Co-loss ─► AR (T878A LBD mutation) ────► Resistant t-NEPC Subclone │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

1. **상호 배타성 검증 (Mutual Exclusivity)**:
   - 전립선암 단일 세포 코호트에서 *TMPRSS2-ERG* 융합 양성 세포 집단과 *SPOP* 돌연변이 양성 세포 집단이 상호 배타적인 클러스터를 형성하는지 교차 검증하여 분석 파이프라인의 신뢰도를 보증합니다.
2. **트렁컬-브랜치 유전체 정렬**:
   - 모든 CNA 서브클론이 *TMPRSS2-ERG* 융합 리드를 공유하면 해당 융합은 초기 트렁컬 사건임이 확정됩니다.
   - 반면 *TP53* R273H 돌연변이 분자가 Chr 17p 결손을 지닌 특정 CNA 서브클론에만 농축되어 검출된다면, 이는 질병 진행 과정에서 분기된 후기 브랜치 사건으로 정의됩니다.

---

## 5. Python 실습: `pysam`을 활용한 단일 세포 BAM 파싱 및 SNV/융합 검출 파이프라인

`pysam` 라이브러리를 이용하여 단일핵 전사체 BAM 파일로부터 지정된 유전자좌(*TP53* R273H)의 세포 바코드별 대립유전자 카운트를 집계하고, VAF 매트릭스를 구축하는 실전 파이프라인입니다.

```python
import pysam
import numpy as np
import pandas as pd
from scipy.stats import betabinom

def extract_single_cell_snv(
    bam_path: str,
    chrom: str,
    pos_1based: int,
    ref_allele: str,
    alt_allele: str,
    min_bq: int = 30,
    min_mq: int = 50   # 예시 임계값. STAR 고유 정렬 MAPQ는 255이고 BWA/minimap2류는 60이므로 정렬기별로 조정
) -> pd.DataFrame:
    """
    단일핵 전사체 BAM 파일에서 특정 게놈 유전자좌의 세포 바코드(CB)별 Ref/Alt UMI 카운트를 추출합니다.
    
    Parameters:
    -----------
    bam_path : str
        인덱싱된 BAM 파일 경로 (.bam 및 .bai 필요)
    chrom : str
        염색체 이름 (예: "chr17")
    pos_1based : int
        1-based 게놈 위치 (예: TP53 R273H 좌표 7673802)
    ref_allele : str
        참조 염기 (예: "C")
    alt_allele : str
        변이 염기 (예: "T")
    min_bq : int
        최소 염기 품질 점수 (Base Quality)
    min_mq : int
        최소 매핑 품질 점수 (Mapping Quality)
        
    Returns:
    --------
    df_snv : pd.DataFrame
        세포 바코드(CB)를 인덱스로 갖고 Ref_UMI, Alt_UMI, VAF를 포함하는 데이터프레임
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    pos_0based = pos_1based - 1
    
    # {cell_barcode: {'ref_umis': set(), 'alt_umis': set()}}
    cell_data = {}
    
    # 해당 위치를 포함하는 리드 반복 탐색
    for pileupcolumn in samfile.pileup(chrom, pos_0based, pos_0based + 1, truncate=True, min_base_quality=min_bq):
        if pileupcolumn.pos != pos_0based:
            continue
            
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
                
            alignment = pileupread.alignment
            if alignment.mapping_quality < min_mq:
                continue
                
            # SAM 태그에서 CB 및 UB 추출
            if not alignment.has_tag("CB") or not alignment.has_tag("UB"):
                continue
                
            cb = alignment.get_tag("CB")
            ub = alignment.get_tag("UB")
            
            read_base = alignment.query_sequence[pileupread.query_position].upper()
            
            if cb not in cell_data:
                cell_data[cb] = {"ref_umis": set(), "alt_umis": set()}
                
            if read_base == ref_allele.upper():
                cell_data[cb]["ref_umis"].add(ub)
            elif read_base == alt_allele.upper():
                cell_data[cb]["alt_umis"].add(ub)
                
    samfile.close()
    
    # 통계 요약 및 VAF 계산
    records = []
    for cb, data in cell_data.items():
        n_ref = len(data["ref_umis"])
        n_alt = len(data["alt_umis"])
        total = n_ref + n_alt
        vaf = n_alt / total if total > 0 else np.nan
        records.append({
            "cell_barcode": cb,
            "ref_umi_count": n_ref,
            "alt_umi_count": n_alt,
            "total_umi_count": total,
            "vaf": vaf
        })
        
    df_snv = pd.DataFrame(records).set_index("cell_barcode")
    return df_snv

def call_genotype_beta_binomial(df_snv: pd.DataFrame, error_rate: float = 0.005,
                                min_total_umi: int = 2,
                                priors=(1/3, 1/3, 1/3)) -> pd.DataFrame:
    """
    베타-이항 혼합 모델을 적용하여 각 단일 세포의 유전자형 사후 확률을 계산합니다.
    min_total_umi 미만(기본 2)의 세포는 판정 불가(NaN→UNCERTAIN)로 둡니다(1개 UMI만으로 MUT를 호출하지 않음).
    priors=(Hom Ref, Het, Hom Alt) 사전확률: 기본 균등은 교육용 가정이며, 실제 체세포 변이는 희귀하므로
    ambient RNA·시퀀싱 오류·변이 사전확률을 반영해 조정해야 합니다.
    """
    results = df_snv.copy()
    post_mutant = []
    
    for _, row in results.iterrows():
        k = int(row["alt_umi_count"])
        N = int(row["total_umi_count"])
        
        if N < min_total_umi:
            post_mutant.append(np.nan)
            continue
            
        # 가설 0: Hom Ref (0/0) -> beta(1, 1/error_rate) (error_rate=0.005면 beta(1, 200))
        p_data_g0 = betabinom.pmf(k, N, 1, int(1 / error_rate))
        # 가설 1: Het (0/1) -> beta(10, 10)
        p_data_g1 = betabinom.pmf(k, N, 10, 10)
        # 가설 2: Hom Alt (1/1) -> beta(100, 1)
        p_data_g2 = betabinom.pmf(k, N, int(1 / error_rate), 1)
        
        # 사전 확률 (기본 균등 1/3: 교육용 가정)
        evidence = priors[0] * p_data_g0 + priors[1] * p_data_g1 + priors[2] * p_data_g2
        if evidence == 0:
            post_mutant.append(0.0)
        else:
            p_mut = (priors[1] * p_data_g1 + priors[2] * p_data_g2) / evidence
            post_mutant.append(p_mut)
            
    results["posterior_mutant_prob"] = post_mutant
    results["genotype_call"] = results["posterior_mutant_prob"].apply(
        lambda p: "MUT" if p >= 0.85 else ("WT" if p <= 0.15 else "UNCERTAIN")
    )
    return results

# --- 가상 데이터 실행 시연 ---
if __name__ == "__main__":
    mock_data = pd.DataFrame({
        "cell_barcode": ["Cell_A-1", "Cell_B-1", "Cell_C-1", "Cell_D-1", "Cell_E-1"],
        "ref_umi_count": [12, 0, 4, 1, 0],
        "alt_umi_count": [0, 8, 5, 0, 1],
        "total_umi_count": [12, 8, 9, 1, 1],
        "vaf": [0.0, 1.0, 0.555, 0.0, 1.0]
    }).set_index("cell_barcode")
    
    df_calls = call_genotype_beta_binomial(mock_data)
    print("[Single-cell SNV Genotyping Results]")
    print(df_calls[["total_umi_count", "vaf", "posterior_mutant_prob", "genotype_call"]])
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)]
단일 세포 전사체(snRNA-seq) 데이터에서 체세포 돌연변이(SNV)를 호출할 때, 단순 벌크 DNA 시퀀싱에서 유효한 기준인 "VAF $\approx 0.5 \implies$ 이형접합체(Heterozygous) 변이"라는 등식이 단일 세포 수준에서는 파탄에 이르는 두 가지 분자생물학적/기술적 기전을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 단일 대립유전자 전사 버스팅 (Monoallelic Transcriptional Bursting)**:
  - 진핵세포 유전자의 전사는 연속적이지 않고 무작위적인 온/오프(On/Off) 펄스 형태로 일어나는 '전사 버스팅'을 따릅니다.
  - 이형접합체($0/1$) 변이를 가진 세포일지라도, 세포를 용해하여 RNA를 포획하는 순간에 우연히 정상 대립유전자(Ref) 또는 변이 대립유전자(Alt) 중 오직 한쪽 대립유전자에서만 전사가 활성화되어 있을 수 있습니다.
  - 결과적으로 세포 내에는 단일 종류의 전사체만 존재하게 되어, 실제 게놈은 이형접합체임에도 관측되는 전사체 기반 VAF는 0.0 또는 1.0으로 극단적으로 왜곡(Allelic Dropout)됩니다.
- **2. 시퀀싱 깊이의 극단적 저하 및 UMI 샘플링 오차 (Low Depth & Sampling Bias)**:
  - 단일핵 시퀀싱은 포획 효율의 한계로 인해 대다수의 유전자좌에서 단 1개 내지 3개의 UMI 분자만 포획됩니다.
  - 분자 수가 1개일 경우 관측 가능한 VAF는 물리적으로 0.0 또는 1.0 외에는 존재할 수 없습니다 ($k/1 \in \{0, 1\}$).
  - 따라서 $0.5$라는 중간값은 수천 개의 세포가 섞인 벌크 평균에서만 통계적으로 관측되는 수치일 뿐, 단일 세포 전사체에서는 VAF의 이산적 표본 추출 오차가 발생하므로 베이지안 확률 모델링 없이는 이형접합성을 직접 판정할 수 없습니다.

</details>

---

### Q2. [수학적 원리 규명 (Mathematical Rigor)]
단일 세포에서 $N$개의 UMI 분자가 검출되었을 때, 이형접합체($G=1$, 기대 비율 $p=0.5$) 변이를 가진 세포가 전사 버스팅 또는 샘플링 오차로 인해 변이 분자가 단 1개도 검출되지 않아($k^{\text{alt}} = 0$) 정상 야생형 세포(WT)로 오분류될 가짜 음성(False Negative) 확률 $P(k^{\text{alt}} = 0 \mid G=1, N)$을 단순 이항 모델과 베타-이항(Beta-Binomial) 과포산 모델에서 각각 유도하고, 과포산 계수 $\rho > 0$가 가짜 음성률에 미치는 영향을 비교하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 단순 이항 분포 모델에서의 위음성 확률**:
  - 각 전사 분자가 변이 대립유전자일 확률이 독립적으로 $p = 0.5$라고 가정하면, $N$개의 분자 중 변이 분자가 0개일 확률은 다음과 같습니다:
    $$P_{\text{binom}}(k^{\text{alt}} = 0 \mid N) = \binom{N}{0} (0.5)^0 (1 - 0.5)^N = \left(\frac{1}{2}\right)^N = 2^{-N}$$
  - $N=5$일 때 위음성 확률은 $2^{-5} = 0.03125$ ($3.1\%$)에 불과합니다.
- **2. 베타-이항 과포산 모델에서의 위음성 확률**:
  - 전사 버스팅으로 인해 대립유전자 발현 비율 $\theta$ 자체가 세포 간에 요동치는 현상을 반영하여 $\theta \sim \text{Beta}(\alpha, \beta)$로 모델링합니다. 대칭적 이형접합성이므로 $\alpha = \beta$이며, 상관 계수(과포산 계수)를 $\rho = \frac{1}{2\alpha + 1}$ ($\alpha = \frac{1 - \rho}{2\rho}$)라 둘 수 있습니다.
  - 베타-이항 분포의 주변 확률(Marginal PMF)에서 $k=0$을 대입하면:
    $$P_{\text{betabinom}}(k^{\text{alt}} = 0 \mid N) = \frac{\text{B}(\alpha, N + \beta)}{\text{B}(\alpha, \beta)} = \prod_{j=0}^{N-1} \frac{\beta + j}{\alpha + \beta + j} = \prod_{j=0}^{N-1} \frac{\alpha + j}{2\alpha + j}$$
  - 이를 $\rho$로 표현하면:
    $$P_{\text{betabinom}}(k^{\text{alt}} = 0 \mid N) = \prod_{j=0}^{N-1} \frac{\frac{1-\rho}{2\rho} + j}{\frac{1-\rho}{\rho} + j} = \prod_{j=0}^{N-1} \frac{1 - \rho + 2j\rho}{2(1 - \rho + j\rho)}$$
- **과포산 계수 $\rho$의 영향 비교**:
  - $\rho \to 0$ ($\alpha \to \infty$)이면 단순 이항 모델인 $2^{-N}$으로 수렴합니다.
  - 그러나 전사 버스팅이 심하여 과포산이 커질수록 ($\rho > 0$), 모든 $j \ge 1$에 대해 $\frac{1 - \rho + 2j\rho}{2(1 - \rho + j\rho)} > \frac{1}{2}$가 성립합니다.
  - 예를 들어 강한 버스팅($\rho = 0.5 \implies \alpha = \beta = 0.5$) 환경에서 $N=5$일 때:
    $$P(k=0) = \frac{0.5 \times 1.5 \times 2.5 \times 3.5 \times 4.5}{1 \times 2 \times 3 \times 4 \times 5} \approx 0.246 \quad (24.6\%)$$
  - 즉, 전사 버스팅 과포산이 존재할 경우 단순 이항 모델($3.1\%$)에 비해 위음성 확률이 약 8배 급증하므로, 단일 세포 유전자형 추정 시 반드시 베타-이항 과포산 보정을 적용해야 함이 수학적으로 증명됩니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
원발성 전립선암 환자의 근치적 전립선 절제술(Radical Prostatectomy) 검체로부터 단일핵 전사체(snRNA-seq) 시퀀싱을 수행하여 5,000개의 종양 세포를 분석하였다.
분석 파이프라인에서 다음 두 가지 상호 배타적 분자 아형의 드라이버 변이가 보고되었다:
- Cluster 1 (2,800개 세포): *ERG* 전사체 발현 15배 상향(3' 라이브러리에서는 접합부 분할 리드 포획이 원리적으로 드물어 *ERG* 과발현이 융합의 주된 대리 지표이며, 접합부 리드는 소수만 관찰된다고 가정).
- Cluster 2 (1,400개 세포): *TMPRSS2-ERG* 리드 완전 결손, 대신 *SPOP* 유전자의 exon 6에서 c.397T>G (p.F133V) 돌연변이 UMI 다수 검출.
- 나머지 800개 세포: 두 변이 모두 미검출이나, Chr 8p 결실 및 Chr 10q (*PTEN*) 결손을 공유함.

1) 전립선암의 발암 기전(Oncogenesis)에 기반하여 Cluster 1과 Cluster 2가 단일 종양 내에서 공존하는 현상을 클론 진화(Clonal Evolution) 관점에서 설명하시오.
2) 두 변이가 모두 검출되지 않은 800개 세포의 진화적 정체를 규명하기 위해 생물정보학 연구자로서 추가로 수행해야 할 유전체 분석 단계 2가지를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 단일 종양 내 공존 현상의 클론 진화적 해석**:
  - 전립선암은 대표적인 **다중 중심성(Multifocal) 종양**입니다.
  - 대규모 유전체 연구(TCGA, Barbieri et al.)에서 입증된 바와 같이, *TMPRSS2-ERG* 융합과 *SPOP* 점돌연변이는 분자생물학적으로 엄격히 **상호 배타적(Mutually Exclusive)**인 발암 드라이버 사건입니다. SPOP 돌연변이 단백질은 야생형 ERG를 유비퀴틴화시켜 분해하지 못하는 등 두 경로는 세포 내에서 서로 양립할 수 없는 상이한 발암 경로를 구동합니다.
  - 따라서 Cluster 1과 Cluster 2가 동일 검체 내에서 공존한다는 것은, 하나의 암세포가 진화하는 과정에서 두 변이가 순차적으로 갈라진 것이 아니라, 전립선 상피 조직의 서로 다른 위치에서 독립적인 전암 병변으로부터 각각 독자적으로 발생한 **이차 독립 기원(Independent Multiclonal Origins / Collision Tumor)**의 결과임을 시사합니다.
- **2. 미검출 800개 세포의 정체 규명을 위한 추가 분석 단계**:
  - **단계 1: 전사체 드롭아웃 보정을 위한 Pseudo-bulk 및 DNA 일치성 검증**:
    - 800개 세포가 Chr 8p 및 Chr 10q 결손을 공유하고 있다면, 이들은 이미 악성화된 종양 세포입니다. 이들 세포에서 *TMPRSS2-ERG*나 *SPOP* 변이가 검출되지 않은 것이 단순 전사체 드롭아웃(Dropout) 때문인지 확인하기 위해, 800개 세포의 리드를 병합(Pseudo-bulk aggregation)하여 극미량의 접합부 리드나 변이 UMI가 존재하는지 재스캔하거나, 단일 세포 DNA 시퀀싱 또는 환자 맞춤형 정밀 타깃 패널(Targeted Capture)로 유전형을 확정해야 합니다.
  - **단계 2: 제3의 상호 배타적 드라이버 변이 탐색 (FOXA1, IDH1, ETV1/4 융합)**:
    - 전립선암 발암 분류 체계에 따라, *ERG* 음성 및 *SPOP* 음성 종양은 종종 다른 ETS 계열 융합 유전자(*ETV1*, *ETV4*, *FLI1*)나 개척자 인자인 *FOXA1* 변이, 혹은 *IDH1* 돌연변이에 의해 발생합니다.
    - 따라서 BAM 파일 전체를 대상으로 *ETV1/4* 키메라 리드 재스캔 및 *FOXA1* 코딩 영역 변이 호출을 수행하여, 이 800개 세포가 독립된 제3의 발암 클론인지 혹은 공통 조상 줄기세포로부터 미분화된 상태인지를 판별해야 합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Cancer Genome Atlas Research Network. The Molecular Taxonomy of Primary Prostate Cancer. *Cell*. 2015;163(4):1011-25. PMID: 26544944. doi: 10.1016/j.cell.2015.10.025
2. Barbieri CE, Baca SC, Lawrence MS, et al. Exome sequencing identifies recurrent SPOP, FOXA1 and MED12 mutations in prostate cancer. *Nat Genet*. 2012;44(6):685-9. PMID: 22610119. doi: 10.1038/ng.2279
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
- STAR 2.7.11 (github.com/alexdobin/STAR). 원 논문: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15-21. PMID: 23104886. doi: 10.1093/bioinformatics/bts635
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- pysam 0.24.1 (github.com/pysam-developers/pysam).
<!-- REFERENCES:END -->

# Chapter 08. 유전자 융합 검출 및 *TMPRSS2-ERG* 키메라 리드 분석 (Gene Fusion Detection & TMPRSS2-ERG Chimeric Mining)

> **Cognitive Goal**: 서구 원발 전립선암 환자의 대략 40~50%에서 발암 드라이버로 작용하는 *TMPRSS2-ERG* 유전자 융합(Gene Fusion)의 게놈 간질성 결실(Interstitial Deletion on 21q22.2) 기전과 전사체 수준의 키메라 접합부 해부학을 정립한다. 정렬 BAM 파일에서 유전자 융합을 입증하는 두 가지 핵심 신호인 **분할 리드(Split-Reads)**와 **불일치 리드 쌍(Discordant Read Pairs)**의 바이너리 CIGAR 및 플래그(FLAG 0x800 Supplementary Alignment) 구조를 디코딩한다. 역전사 템플릿 스위칭(RT Template Switching)에 의한 위양성 키메라 인공산물을 단일 세포 UMI 및 미세 상동성(Microhomology) 필터로 제거하는 통계적 판정 알고리즘을 구축한다. 최종적으로 융합 양성 단일 세포를 염색체 21q 국소 결실(CNA) 서브클론 계통(Track 04)과 정밀 결합하는 멀티오믹스 통합 파이프라인을 체득한다.

---

## 1. 전립선암의 대표적 드라이버 구조 변이: *TMPRSS2-ERG* 융합의 분자생물학

*TMPRSS2-ERG* 유전자 융합은 인간 고형암에서 발견되는 가장 흔한 종양 특이적 재배열(Genomic Rearrangement) 사건 중 하나로, 전립선 선암종(Prostate Adenocarcinoma) 환자의 대략 40~50%(코호트·인종 의존)에서 발생합니다.

인간 염색체 21번의 장완(`21q22.2~21q22.3`)에 위치한 두 유전자는 약 $3\text{ Mb}$ 떨어져 있으며 두 유전자 모두 − 가닥에 있어 같은 전사 방향(Telomere $\to$ Centromere)을 향합니다(*ERG*가 동원체 쪽, *TMPRSS2*가 텔로미어 쪽에 위치):
- **5' 파트너 유전자 (*TMPRSS2*)**: 안드로겐 수용체(AR) 신호에 의해 전립선 상피에서 강하게 발현되는 막관통 세린 프로테아제 (Androgen-regulated gene).
- **3' 파트너 유전자 (*ERG*)**: 정상 전립선 상피에서는 발현되지 않으나 혈관 내피 및 조혈계에서 세포 증식 및 분화를 조절하는 ETS 전사인자 계열 원종양유전자(Proto-oncogene).

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               염색체 21q22.2 구간의 ~3Mb 간질성 결실(Edel)에 의한 융합 기전            │
│                                                                                        │
│   Centromere ──────────────────────────────────────────────────────────► Telomere      │
│                [ERG]               (~3 Mb Interstitial Region)        [TMPRSS2]        │
│                Exon 4-17              ETS2, BRWD1, HMGN1 등         Exon 1-2           │
│                ┌─────────┐      ──────────────────────────────     ┌─────────┐         │
│                │ Coding  │ ◄─── [   3 Mb Genomic Deletion    ] ◄── │  5' UTR │         │
│                └─────────┘      ──────────────────────────────     └─────────┘         │
│                     │                                                   │              │
│                     └───────────────────────┬───────────────────────────┘              │
│                                             ▼                                          │
│   [TMPRSS2-ERG Fusion Transcript]                                                      │
│   ┌─────────────────────┬──────────────────────────────────────────────────────────┐   │
│   │ TMPRSS2 Exon 1 / 2  │                       ERG Exon 4 - 17                    │   │
│   │ 안드로겐 유도 프로모터│ 완전한 ETS 전사인자 ORF (강력한 종양 성장 및 침습 유도)     │   │
│   └─────────────────────┴──────────────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 게놈 재배열의 두 가지 메커니즘
1. **간질성 결실 (Deletion, Edel)**: 융합 양성 종양의 약 60%에서 관찰되며, *TMPRSS2*와 *ERG* 사이의 약 3Mb 게놈 구간이 완전히 잘려나가면서 두 유전자가 물리적으로 인접하게 됩니다. 이 결실은 DNA 수준(WGS, SNP 어레이, FISH)에서 21q22 국소 결실로 관측됩니다. 다만 scRNA/snRNA 기반 CNA(inferCNV, CopyKAT)로는 해상도 밖입니다(두 유전자 사이 RefSeq Select 유전자가 약 15개로 inferCNV 기본 창 약 100유전자·CopyKAT 약 5 Mb보다 작음).
2. **염색체 내 전좌 (Translocation, Etrans)**: 융합 양성 종양의 나머지(약 40%)에서 관찰되며, 결실 없이 염색체 역위(Inversion)나 전좌에 의해 융합이 형성되므로 카피수 변화가 수반되지 않습니다.

### 1.2 병리학적 결과
안드로겐(테스토스테론)이 분비되면, 정상 전립선 상피에서는 *TMPRSS2*만 발현되어야 하지만, 융합 전사체에 의해 *ERG* 종양 단백질이 수십~수백 배 폭발적으로 과발현됩니다. 과발현된 ERG는 세포 외 기질 분해 효소(MMP)와 상피-간엽 이행(EMT) 경로를 활성화하여 전립선암의 침습(Invasion), 혈관 신생, 전이 능력을 촉진합니다.

---

## 2. 융합 전사체 해부학: 주요 스플라이싱 이소폼

게놈 DNA 수준의 절단점(Breakpoint)은 대개 *TMPRSS2*의 인트론 1 또는 2와 *ERG*의 거대 인트론 3 내부에서 무작위로 발생합니다. 그러나 전사 후 RNA 스플라이싱 기전을 거치면서 성숙 전사체(Mature mRNA)는 정밀한 엑손 경계 접합부(Exon-Exon Junction)로 재편됩니다.

> **주의**: 아래 표의 'Type' 번호와 빈도(~65/25/5/5%)는 출처를 확인하지 못한 *예시*이며 문헌의 Type I/II 분류와 다릅니다.

| 전사체 유형 | 5' 스플라이스 공여부 | 3' 스플라이스 수용부 | 전사 판독 프레임 (ORF) | 임상 빈도 (예시, 출처 미확인) |
|---|---|---|---|---|
| **Type III (T1-E4)** | *TMPRSS2* Exon 1 | *ERG* Exon 4 | 정규 프레임 (ERG Met codon @ Exon 4) | **가장 흔함 (~65%)** |
| **Type VI (T2-E4)** | *TMPRSS2* Exon 2 | *ERG* Exon 4 | 정규 프레임 유지 | 흔함 (~25%) |
| **Type I (T1-E2)** | *TMPRSS2* Exon 1 | *ERG* Exon 2 | 프레임시프트 또는 대체 개시 | 드묾 (~5%) |
| **Type IV (T2-E5)** | *TMPRSS2* Exon 2 | *ERG* Exon 5 | N-말단 절단형 ERG 발현 | 드묾 (~5%) |

가장 지배적인 형태는 **Type III (*TMPRSS2* Exon 1 $\to$ *ERG* Exon 4)**입니다:
- *TMPRSS2* 엑손 1은 순수한 비번역 영역(5' UTR)입니다.
- *ERG* 엑손 4는 ERG 단백질의 첫 번째 개시 코돈(ATG, Methionine 1)을 포함하고 있습니다.
- 따라서 융합 전사체는 키메라 펩타이드 없이 ERG를 생산하되, 가장 흔한 T1-E4(Tomlins 표기, exon 4) 융합 단백은 **N-말단 약 32개 아미노산이 결손된 near-full-length ERG**입니다(완전한 천연 ERG가 아님; PMC4162137).

---

## 3. 정렬 BAM 파일 내 융합 신호: 분할 리드(Split-Reads) vs 불일치 리드 쌍(Discordant Pairs)

단일 세포 전사체 데이터에서 융합 유전자를 식별하는 결정적 증거는 정렬 레코드의 기하학적 이상 정렬입니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   BAM 파일 내 유전자 융합 검출의 두 가지 물리적 증거                   │
│                                                                                        │
│   [1. Split-Read (분할 리드 / Chimeric Alignment)]                                     │
│   단일 리드 서열(150bp)이 융합 접합점을 직접 가로지름:                                │
│                                                                                        │
│   TMPRSS2 Exon 1                   ERG Exon 4                                          │
│   ───────────────┬───────────────────────────────────                                  │
│   ... G G C A G C│A T T C A G G A C T C A G T C C ...                                  │
│   ▲              ▲                  ▲                                                  │
│   [ 45bp Match ] │  [ 105bp Soft-clip (Match on ERG) ]                                 │
│                                                                                        │
│   Primary Alignment:       POS=TMPRSS2, CIGAR=45M105S, FLAG=99                         │
│   Supplementary Alignment: POS=ERG,     CIGAR=45S105M, FLAG=2147 (0x800 set), SA 태그   │
│   (개념도: 실제 hg38에서는 − 가닥 유전자라 SAM 저장 순서상 ERG가 먼저 옴. Q1 참조)  │
│                                                                                        │
│   [2. Discordant Read Pair (불일치 리드 쌍)]                                           │
│   Read 1 (50bp) ──────► [TMPRSS2]                Read 2 (50bp) ◄────── [ERG]           │
│   (정상 삽입 단편 크기 300bp를 초과하여 ~3Mb 떨어지거나 이종 염색체 간 정렬)          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

> **방향 주의**: 위 도식은 개념 설명용으로 유전자 순서와 저장 서열 순서를 단순화했다. hg38에서 *TMPRSS2*와 *ERG*는 모두 − 가닥이므로 + 가닥 기준으로 저장되는 SAM 서열에서는 *ERG* 쪽 서열이 먼저 오고 *TMPRSS2* 쪽이 뒤에 온다(Q1 참조). 또한 *TMPRSS2*(5′ 파트너)는 *ERG*보다 텔로미어 쪽에 있다.

### 3.1 분할 리드 (Split-Reads)의 바이너리 명세
1. **주 정렬 (Primary Alignment)**:
   - 참조 게놈 상의 *TMPRSS2* 유전자좌에 매핑됩니다.
   - CIGAR 문자열에 소프트 클리핑(`S`)이 발생합니다 (예: `45M105S`).
   - 보조 태그에 **`SA:Z:` (Supplementary Alignment)** 필드가 자동 추가되어 보충 정렬이 위치한 파트너 염색체 좌표, 방향, CIGAR, MAPQ를 기록합니다:
     `SA:Z:chr21,41508081,+,110S40M,60,0;` (Q1의 hg38 예시 좌표와 일치하는 예; 위 도식은 개념용 단순화)
2. **보충 정렬 (Supplementary Alignment)**:
   - 참조 게놈 상의 *ERG* 유전자좌에 매핑됩니다.
   - SAM 플래그에 **`0x0800` (Decimal 2048)** 비트가 활성화됩니다.
   - 플래그 값 예: $99 + 2048 = 2147$.

### 3.2 10x Genomics snRNA-seq에서의 검출 특수성
표준 10x 3' snRNA-seq은 전사체의 3' 말단(PolyA 테일) 부근 약 $300\sim 500\text{bp}$ 단편만을 주로 포획합니다. *TMPRSS2-ERG* 융합 접합부는 전사체의 5' 쪽(*TMPRSS2* 엑손 1 → *ERG* 엑손 4 부근)에 위치하며 *ERG* 엑손 4에서 poly(A)까지는 약 $4.7\text{ kb}$(엑손 4~11 합계 4,726 nt, hg38 NM_004449.4 기준)나 떨어져 있으므로, 3' 편향 라이브러리에서는 접합부를 직접 통과하는 분할 리드가 원리적으로 거의 포획되지 않습니다.
따라서 3' snRNA-seq에서는 *ERG* 엑손 4 이후 전사체의 발현 증가(및 정상 전립선 상피에서의 *ERG* 부재)를 융합의 대리 지표로 사용하고, 접합부 자체를 확인하려면 5' 라이브러리, 롱리드 또는 타깃 증폭, DNA 수준 검증이 필요합니다(비정상적 내부 폴리아데닐화 프라이밍 등으로 소수의 융합 리드가 우연히 포획될 수는 있으나 기대할 수 없습니다).

---

## 4. 단일 세포 위양성(Artifact) 제거 알고리즘

단일 세포 시퀀싱의 역전사(RT) 반응에는 Moloney Murine Leukemia Virus (M-MLV) 역전사 효소가 사용됩니다. 이 효소는 RNA 주형을 복사하는 도중 서열 상동성이 있는 다른 부위로 점프하는 **주형 교환(Template Switching)** 현상을 일으키며, 이는 실제 세포 내에 존재하지 않는 인공적 키메라 cDNA(Artifactual Chimeras)를 대량 생성합니다.

이를 필터링하기 위한 4단계 정밀 검증 수식입니다:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        위양성 융합 리드 필터링 4단계 결정 트리                         │
│                                                                                        │
│   후보 분할 리드 검출 (CIGAR with 'S' & SA tag)                                        │
│          │                                                                             │
│          ▼ [검증 1: 미세 상동성 (Microhomology) 길이 검사]                             │
│   접합부 경계 상동 서열 길이 L_homo >= 5 bp ? ──► [YES] ──► 의심(RT Artifact 가능)     │
│          │ [NO]                                                                        │
│          ▼ [검증 2: 최소 오버행 (Anchor Length) 검사]                                  │
│   min(L_donor, L_acceptor) < 15 bp ? ──────────► [YES] ──► 즉시 폐기 (False Alignment)│
│          │ [NO]                                                                        │
│          ▼ [검증 3: 분자 중복성 (Independent UMI Support)]                             │
│   동일 CB 내 독립 UB 분자 수 N_umi < 2 ? ──────► [YES] ──► 보류 (Low-confidence)       │
│          │ [NO, N_umi >= 2]                                                            │
│          ▼ [검증 4: 정상 참조 세포군 대조 (Normal Negative Control)]                   │
│   정상 기질/상피 세포(Normal cells)에서 검출 ? ─► [YES] ──► 즉시 폐기 (Germline/Noise)  │
│          │ [NO]                                                                        │
│          ▼                                                                             │
│   [후보 TMPRSS2-ERG 단일 세포 융합 양성 (Fusion+); 직교 검증 필요]                    │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

1. **미세 상동성 검사**: 파트너 유전자 간의 접합부 경계에 존재하는 상동 서열 길이 $L_{\text{homo}} < 5\text{bp}$인 리드만 인정.
2. **최소 오버행 길이**: 접합부 양쪽 유전자에 걸친 리드의 매핑 길이 $\min(L_1, L_2) \ge 15\text{bp}$.
3. **독립 UMI 수**: $N_{\text{UMI}}(\text{cell}) \ge 2$.
4. **정상 기질 음성 대조군**: 비종양 내피세포 및 섬유아세포에서 검출되지 않을 것(관찰 기준; 재발성 템플릿 스위칭은 독립 UMI로도 같은 부위에 반복될 수 있어 위 기준이 위양성을 완전히 제거하지는 못함).

> **구현 한계**: 아래 코드는 TMPRSS2→ERG 키메라 SA 항목을 계수할 뿐 exon1–exon4 접합부 확인, 미세 상동성 검사, 정상 세포 대조를 구현하지 않으며, STAR에서 `--chimSegmentMin > 0`과 `--chimOutType WithinBAM`을 켠 BAM에서만 SA가 존재합니다(02장 기본 명령 BAM에는 SA 태그가 없어 빈 결과).

---

## 5. 융합 양성 단일 세포와 카피수 변이(CNA) 서브클론의 유전체 통합 (Track 04 연계)

*TMPRSS2-ERG* 융합 분석의 진정한 가치는 단일 전사체 수준의 융합 상태를 세포별 게놈 대규모 카피수 변이(CNA) 계통수와 결합할 때 폭발합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│             CNA 계통 트리와 TMPRSS2-ERG 단일 세포 융합 상태의 공클러스터링             │
│                                                                                        │
│   [정상 2배체 세포군]                                                                  │
│   Normal Prostate Epithelial / Stromal ──► Chr 21 Neutral (2 Copies) ──► Fusion (-)   │
│                                                                                        │
│   [종양 클론 진화 트리]                                                                │
│   Founder Truncal Cancer Clone                                                         │
│   ├── Chr 8p Loss (-1 Copy)                                                            │
│   ├── Chr 21q22.2 Interstitial Focal Deletion (-1 Copy, Edel) ──► TMPRSS2-ERG Fusion(+)│
│   │                                                                                    │
│   ├──► Branch Subclone 1 : + PTEN Loss (10q23 Del) ────────────► Fusion (+)           │
│   │                                                                                    │
│   └──► Branch Subclone 2 : + TP53 Mutation (17p Loss) ─────────► Fusion (+)           │
│                                                                                        │
│   결론: 모든 하위 악성 서브클론이 21q22 결실과 융합을 공유하므로,                      │
│         TMPRSS2-ERG는 이 가상 예시에서 '트렁컬(Truncal)' 조상 변이로 시사됨 (예시)   │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

- **이형 접합 결실 vs 동형 접합 결실**:
  - 대부분의 전립선암은 한쪽 대립유전자에서 3Mb 결실이 발생하여 융합이 형성됩니다 ($CNA = 1$).
  - 일부 종양에서는 융합이 중복되고 간질성 결실이 겹쳐 나타나며('2+Edel'), Attard et al. (Oncogene 2008, PMID 17637754)의 FISH 연구에서 이 군이 특히 불량한 예후와 연관되었습니다(0카피 동형 결실이 아닌 융합 중복+결실). 0카피 서술은 근거가 없어 삭제합니다.

---

## 6. Python 실습: `pysam` 기반 *TMPRSS2-ERG* 분할 리드 마이닝 파이프라인

단일핵 BAM 파일로부터 `SA` 태그를 보유한 보충 정렬 리드를 스캔하여, *TMPRSS2* 엑손 1과 *ERG* 엑손 4 접합부를 통과하는 참 양성 분할 리드를 추출하고 세포별 UMI 매트릭스를 생성하는 코드입니다.

```python
#!/usr/bin/env python3
"""
Single-Cell TMPRSS2-ERG Gene Fusion Mining Pipeline using pysam
Author: Senior Bioinformatics Engineer
Target: TMPRSS2 (chr21:41,460,000-41,530,000) - ERG (chr21:38,300,000-38,500,000)
"""

from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import pysam
import pandas as pd


def parse_sa_tag(sa_str: str) -> List[Tuple[str, int, str, str, int]]:
    """
    SA 태그 문자열을 파싱하여 (rname, pos, strand, cigar, mapq) 튜플 리스트를 반환합니다.
    Format: (chr, pos, strand, CIGAR, MAPQ, NM;)
    """
    entries = []
    for entry in sa_str.strip(";").split(";"):
        if not entry:
            continue
        tokens = entry.split(",")
        if len(tokens) >= 5:
            rname = tokens[0]
            pos = int(tokens[1])
            strand = tokens[2]
            cigar = tokens[3]
            mapq = int(tokens[4])
            entries.append((rname, pos, strand, cigar, mapq))
    return entries


def mine_tmprss2_erg_fusions(
    bam_path: str,
    tmprss2_region: Tuple[str, int, int] = ("chr21", 41460000, 41530000),
    erg_region: Tuple[str, int, int] = ("chr21", 38300000, 38500000),
    min_anchor: int = 15,
    min_mapq: int = 20
) -> pd.DataFrame:
    """
    단일핵 BAM에서 TMPRSS2-ERG 융합 지지 분할 리드를 마이닝하고 세포별 UMI를 도출합니다.
    """
    print(f"[*] Mining TMPRSS2-ERG chimeric alignments in {bam_path}...")
    
    t_chr, t_start, t_end = tmprss2_region
    e_chr, e_start, e_end = erg_region
    
    # 융합 분자 보관함: (Cell_Barcode, UMI) -> List of split-read details
    fusion_molecules: Dict[Tuple[str, str], list] = defaultdict(list)
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(t_chr, t_start, t_end):
            # 1. 매핑 품질 및 단일 세포 태그 검증
            if read.mapping_quality < min_mapq or read.is_duplicate:
                continue
            if not (read.has_tag("CB") and read.has_tag("UB")):
                continue
                
            cb = read.get_tag("CB")
            ub = read.get_tag("UB")
            
            # 2. 소프트 클리핑('S')이 최소 15bp 이상 존재하는지 검사
            cigartuples = read.cigartuples
            if not cigartuples:
                continue
            has_sufficient_clip = any(op == 4 and length >= min_anchor for op, length in cigartuples)
            if not has_sufficient_clip:
                continue
                
            # 3. Supplementary Alignment (SA 태그) 검사
            if not read.has_tag("SA"):
                continue
                
            sa_tag = read.get_tag("SA")
            sa_entries = parse_sa_tag(sa_tag)
            
            # 4. SA 대상이 ERG 유전자좌에 정확히 매핑되는지 검증
            for s_rname, s_pos, s_strand, s_cigar, s_mapq in sa_entries:
                if s_rname == e_chr and (e_start <= s_pos <= e_end) and s_mapq >= min_mapq:
                    # 진성 키메라 분할 리드 등록
                    fusion_molecules[(cb, ub)].append({
                        "qname": read.query_name,
                        "tmprss2_pos": read.reference_start + 1,
                        "tmprss2_cigar": read.cigarstring,
                        "erg_pos": s_pos,
                        "erg_cigar": s_cigar
                    })

    print(f"[*] Identified {len(fusion_molecules)} unique fusion molecules (CB-UB).")
    
    # 세포별 집계 및 신뢰도 점수 부여
    cell_fusion_stats = defaultdict(lambda: {"supported_umi_count": 0, "read_names": []})
    for (cb, ub), reads in fusion_molecules.items():
        cell_fusion_stats[cb]["supported_umi_count"] += 1
        cell_fusion_stats[cb]["read_names"].extend([r["qname"] for r in reads])
        
    records = []
    for cb, stats in cell_fusion_stats.items():
        n_umi = stats["supported_umi_count"]
        # 엄격 기준: 2개 이상의 독립 UMI 지지 시 High_Confidence_Fusion_Positive
        call = "High_Confidence_Fusion_Positive" if n_umi >= 2 else "Candidate_Low_Confidence"
        records.append({
            "cell_barcode": cb,
            "fusion_umi_count": n_umi,
            "fusion_status": call
        })
        
    return pd.DataFrame(records)


if __name__ == "__main__":
    import sys
    test_bam = sys.argv[1] if len(sys.argv) > 1 else "data/sample.bam"
    try:
        df_fusions = mine_tmprss2_erg_fusions(test_bam)
        print(df_fusions.head(10))
    except FileNotFoundError:
        print(f"[!] Target BAM {test_bam} not found.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [수학적 원리 규명 (Mathematical Rigor)]
150bp 페어드 엔드 리드 시퀀싱에서 단일핵 전사체 리드가 *TMPRSS2-ERG* 융합 접합부에 매핑되었다. (*TMPRSS2*와 *ERG*는 모두 hg38에서 − 가닥이며, 융합 전사체의 SAM 저장 서열은 + 가닥 기준이므로 *ERG* 쪽 서열이 먼저 오고 *TMPRSS2* 쪽이 뒤에 온다. 좌표는 hg38 RefSeq 기준이며 리드는 예시이다.)
주 정렬(Primary Alignment)은 `chr21:38,445,512` (*ERG* Exon 4, exon 4 = 38,445,404–38,445,621)에 위치하며, CIGAR 문자열은 `110M40S`이다.
보충 정렬(Supplementary Alignment)은 `chr21:41,508,081` (*TMPRSS2* Exon 1, exon 1 = 41,508,081–41,508,158)에 위치하며, CIGAR 문자열은 `110S40M`이다.
1) 주 정렬과 보충 정렬 레코드의 SAM `FLAG` 값에서 각각 활성화되어야 하는 고유 플래그(16진수 및 10진수)를 명시하시오.
2) 주 정렬의 110bp 매치 구간과 보충 정렬의 40bp 매치 구간이 리드 전체 150bp 서열을 빈틈없이 설명함을 증명하고, 두 유전자가 − 가닥임을 고려하여 접합부의 1-기반 게놈 좌표(Breakpoint coordinates)를 계산하시오.
3) 만약 보충 정렬의 CIGAR가 `110S40M`이 아니라 (a) `114S36M`이라면, (b) `106S44M`이라면 각각 두 정렬 사이에 무슨 일이 일어난 것인지 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 주 정렬 및 보충 정렬의 고유 SAM FLAG 명세**:
  - **주 정렬 (Primary)**:
    - 대표 정렬이므로 `0x0100` (Secondary) 및 `0x0800` (Supplementary) 플래그가 모두 0이어야 합니다.
    - 페어드 엔드 R1 순가닥 매핑 기준 기본 플래그 예: `0x0001` (PAIRED) + `0x0002` (PROPER_PAIR) + `0x0020` (MREVERSE, 짝 리드는 역가닥) + `0x0040` (READ1) = **`99` (Decimal)** (0x0001+0x0002+0x0040만이면 67).
  - **보충 정렬 (Supplementary)**:
    - 키메라 분할 리드의 잔여 정렬이므로 **`0x0800` (Decimal 2048, SUPPLEMENTARY)** 플래그가 반드시 1로 설정되어야 합니다.
    - 따라서 기본 플래그에 2048을 더한 **`2147` (Decimal, 0x0863)** 값을 갖습니다.

- **2. 리드 서열 소모 및 절단점 게놈 좌표 증명**:
  - 주 정렬: 리드의 $1 \sim 110\text{번째}$ 염기가 *ERG*에 매치(`110M`), 이후 $111 \sim 150\text{번째}$ 염기는 클리핑(`40S`).
  - 보충 정렬: 리드의 $1 \sim 110\text{번째}$ 염기가 클리핑(`110S`), 이후 $111 \sim 150\text{번째}$ 염기(총 40bp)가 *TMPRSS2*에 매치(`40M`).
  - $110 + 40 = 150\text{ bp}$로 리드의 전장 150bp가 누락이나 갭 없이 100% 설명됩니다.
  - **절단점 게놈 좌표 (두 유전자 모두 − 가닥)**: − 가닥 유전자에서는 전사 방향이 좌표 감소 방향이므로, 5′ 파트너 *TMPRSS2*의 공여점(exon 1 말단)은 해당 세그먼트의 **왼쪽 끝**이고, 3′ 파트너 *ERG*의 수용점(exon 4 시작)은 해당 세그먼트의 **오른쪽 끝**이다(+ 가닥 가정의 $\text{POS}+L-1$ 공식을 그대로 쓰면 안 된다):
    $$\text{TMPRSS2 공여점} = \text{POS}_{\text{supp}} = 41{,}508{,}081 \quad (\text{chr21:41,508,081})$$
    $$\text{ERG 수용점} = \text{POS}_{\text{prim}} + L - 1 = 38{,}445{,}512 + 110 - 1 = 38{,}445{,}621 \quad (\text{chr21:38,445,621})$$

- **3. 두 정렬 사이 서열 관계가 시사하는 현상**:
  - (a) `110M40S` + `114S36M`: 리드 111~114번째 4염기는 어느 정렬에도 포함되지 않는다($110+36=146\neq150$). 이는 중복이 아니라 접합부의 **4bp 비주형 삽입(N-region)**이며, NHEJ 말단 가공이나 말단 전이효소 유사 활성으로 설명될 수 있다.
  - (b) `110M40S` + `106S44M`: 리드 107~110번째 4염기가 두 정렬에 모두 포함되어 4bp가 **겹친다**. 이는 두 유전자의 경계에 공통 4bp 서열이 있는 **미세 상동성(Microhomology)**을 뜻하며, MMEJ/NHEJ에 의한 재결합이거나 역전사효소(RT)의 주형 교환(Template Switching) 인공산물일 가능성이 있다(둘은 이 데이터만으로 구분할 수 없고 DNA 수준 검증이 필요하다).

</details>

---

### Q2. [알고리즘 구현 및 위양성 판정 (Algorithmic Logic)]
단일핵 시퀀싱 데이터 분석 도중 전립선암 환자 검체에서 다음과 같은 두 개의 분할 리드 클러스터가 검출되었다:
- Candidate A: 종양 세포 바코드(`CB_1`)에서 동일한 UMI(`UB_A`)를 가진 8개의 PCR 중복 리드가 검출됨. 접합부 양측 매핑 길이는 각각 75bp, 75bp이며 미세 상동성이 0bp임. 그러나 해당 `CB_1`에서 다른 UMI를 가진 융합 리드는 전무함.
- Candidate B: 종양 세포 바코드(`CB_2`)에서 서로 다른 3개의 독립된 UMI(`UB_1`, `UB_2`, `UB_3`)를 가진 분할 리드가 각각 2개씩 총 6개 검출됨.
1) Candidate A와 Candidate B 중 실제 종양 세포 내에서 전사된 진성 융합 유전자(True Positive)일 확률이 압도적으로 높은 후보를 선택하고, PCR 증폭 편향(PCR Duplication Bias) 관점에서 근거를 논리적으로 입증하시오.
2) Candidate A와 같은 단일 UMI 분할 리드가 단일 세포 시퀀싱에서 흔하게 관찰되는 기술적 원인을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 참 양성 후보 선택 및 논리적 입증**:
  - **선택: Candidate B**.
  - **입증**:
    - Candidate A는 총 리드 수가 8개로 많아 보이지만, 모두 **동일한 단 1개의 UMI(`UB_A`)**를 공유하는 PCR 중복 리드(Duplicates)에 불과합니다. 즉, 라이브러리 제작 과정에서 최초 단 1개의 cDNA 분자가 증폭된 것입니다. 단 1분자만 관측된 사건은 역전사 반응 중 1회성으로 발생한 RT 주형 스위칭(Artifactual chimera)이나 시퀀싱 오류일 확률을 통계적으로 배제할 수 없습니다.
    - 반면 Candidate B는 **서로 다른 3개의 독립적인 UMI 분자(`UB_1`, `UB_2`, `UB_3`)**가 동일한 융합 접합점을 지지하고 있습니다. 세 개의 서로 다른 원본 mRNA 분자가 세포 내에서 독립적으로 역전사되어 동일한 접합부를 형성했다는 것은, 물리적 확률상 효소의 우연한 에러일 가능성이 $P < 10^{-6}$ 이하로 수렴하므로 진성 융합 전사체임이 통계학적으로 확증됩니다.

- **2. Candidate A와 같은 단일 UMI 아티팩트의 생성 원인**:
  - 역전사 효소(M-MLV RT)는 높은 가공성(Processivity)을 갖지만, 2차 구조가 강한 RNA 영역이나 유사 서열을 만났을 때 원래 주형에서 일시적으로 해리된 후 인접한 다른 cDNA/RNA 가닥에 결합하여 합성을 이어가는 '주형 전이(Template Switching)' 에러를 빈번하게 일으킵니다.
  - 이로 인해 생성된 1개의 가짜 키메라 cDNA 분자가 이후 12~14 사이클의 PCR 지수 증폭을 거치면서 수십 개의 동일 UMI 리드로 복제되어 나타나는 것입니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
원발성 전립선암 환자의 전립선 전절제술 조직 검체로부터 단일핵 전사체(snRNA-seq) 6,000개 세포를 분석하여 대규모 카피수 변이(CNA) 클론 계통수(InferCNV)와 유전자 융합 상태를 결합하였다.
데이터 분석 결과 다음과 같은 서브클론 구조가 도출되었다:
- Subclone 1 (2,500개 세포): Chr 8p 결실, Chr 21q22.2 결실 (1 copy loss; 가상 예시이며 실제 inferCNV로는 이 해상도의 결실을 판정할 수 없음), *TMPRSS2-ERG* 융합 리드 검출률 95% (가상 수치; 3′ 라이브러리에서는 비현실적).
- Subclone 2 (2,000개 세포): Chr 8p 결실, Chr 21q22.2 정상 2배체 유지, *TMPRSS2-ERG* 리드 미검출, 대신 *SPOP* p.F133V 점돌연변이 검출률 90% (가상 수치; 3′ 라이브러리에서는 비현실적).
- Subclone 3 (1,500개 세포): 정상 2배체 (CNA 변이 없음, 두 유전자 변이 모두 음성).
1) Subclone 1과 Subclone 2의 분자생물학적 상호 배타성(Mutual Exclusivity)이 시사하는 발암 기전을 설명하고, 이 환자의 전립선암이 단일 클론 기원(Monoclonal)인지 다초점성 다중 기원(Multifocal / Multiclonal)인지 판정하시오.
2) 만약 Subclone 1 내의 특정 500개 세포에서 Chr 21q22.2 영역이 0카피(Biallelic complete loss)로 관측되었다면, 이 관측을 어떻게 해석·검증해야 하는지, 그리고 문헌상 불량 예후와 연관된 21q22 변이 형태는 무엇인지 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 상호 배타성 기전 및 종양 기원 판정**:
  - **분자 기전**:
    - *SPOP*은 Cullin3 기반 E3 리가아제의 기질 인식 서브유닛으로, 변이(F133V 등)는 기질 결합을 해친다. *ERG*와의 관계에서는 TMPRSS2-ERG에서 만들어진 N-말단 절단 ERG가 SPOP 분해 degron을 잃어 SPOP 매개 분해에 저항한다는 보고가 있다(An et al., Mol Cell 2015, PMID 26344096). 즉 정상 전립선 상피에서 발현되지 않는 야생형 ERG 단백질을 SPOP 변이가 안정화시키는 것이 아니다.
    - *SPOP* 변이 종양과 *ERG* 융합 종양이 서로 다른 분자 아형으로 구분되고 상호 배타적으로 관찰되는 것은 TCGA 수준의 관찰이며(PMID 26544944), 서로 다른 드라이버 경로의 선택으로 해석한다.
  - **종양 기원 판정**:
    - 두 클론은 Chr 8p 결실을 공유하므로 공통 조상에서 분기했을 가능성과, 8p 결실이 독립적으로 두 번 일어난 다초점 기원 가능성이 모두 남는다. 드라이버(*ERG* 융합 vs *SPOP* 변이)가 상호 배타적이라는 사실만으로 단일 클론 여부를 판정할 수 없으며, 다부위 WGS/계통수 분석으로 확인해야 한다.

- **2. Chr 21q22.2 0카피(Biallelic Deletion)의 임상적 함의**:
  - 3 Mb 구간 안의 RefSeq Select 유전자는 약 15개(ETS2, PSMG1, BRWD1, HMGN1, GET1/WRB 등)이며 inferCNV 기본 창(약 100유전자)·CopyKAT(약 5 Mb)보다 작으므로, snRNA-seq CNA로 21q22.2 0카피를 판정하는 것은 해상도 밖이다. 이 판정은 WGS/SNP 어레이/FISH 등 DNA 근거가 필요하다.
  - 불량 예후와 연관된 형태는 0카피 동형 결실이 아니라 '융합 중복 + 간질성 결실(2+Edel)'이다(Attard et al., Oncogene 2008, PMID 17637754; FISH 코호트, 위험비 6.10). *HMGN1*, *DSCR3* 같은 '종양 억제 유전자 소실'은 근거가 없으며 *DSCR3*(VPS26C)는 결실 구간 밖(ERG보다 동원체 쪽, 약 37.2 Mb)에 있다. 따라서 이 500개 세포의 '초고위험' 예후를 단정할 수 없다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Attard G, Clark J, Ambroisine L, et al. Duplication of the fusion of TMPRSS2 to ERG sequences identifies fatal human prostate cancer. *Oncogene*. 2008;27(3):253-63. PMID: 17637754. doi: 10.1038/sj.onc.1210640
2. An J, Ren S, Murphy SJ, et al. Truncated ERG Oncoproteins from TMPRSS2-ERG Fusions Are Resistant to SPOP-Mediated Proteasome Degradation. *Mol Cell*. 2015;59(6):904-16. PMID: 26344096. doi: 10.1016/j.molcel.2015.07.025
3. Tomlins SA, Rhodes DR, Perner S, et al. Recurrent fusion of TMPRSS2 and ETS transcription factor genes in prostate cancer. *Science*. 2005;310(5748):644-8. PMID: 16254181. doi: 10.1126/science.1117679
### 검수 근거 문헌
1. Cancer Genome Atlas Research Network. The Molecular Taxonomy of Primary Prostate Cancer. *Cell*. 2015;163(4):1011-25. PMID: 26544944. doi: 10.1016/j.cell.2015.10.025
2. Zheng GX, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. *Nat Commun*. 2017;8:14049. PMID: 28091601. doi: 10.1038/ncomms14049
3. Rastogi A, Tan SH, Mohamed AA, et al. Functional antagonism of TMPRSS2-ERG splice variants in prostate cancer. *Genes Cancer*. 2014;5(7-8):273-84. PMID: 25221645. doi: 10.18632/genesandcancer.25
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- pysam 0.24.1 (github.com/pysam-developers/pysam).
<!-- REFERENCES:END -->

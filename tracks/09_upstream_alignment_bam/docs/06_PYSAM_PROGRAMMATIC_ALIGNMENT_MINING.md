# Chapter 06. pysam 기반 프로그래밍 방식 정렬 데이터 마이닝 (pysam Programmatic Alignment Mining)

> **Cognitive Goal**: C 언어로 작성된 고성능 HTSlib 라이브러리의 Python 바인딩인 `pysam`의 내부 아키텍처와 Cython 브리지(Bridge) 메모리 수명주기(Memory Lifecycle)를 이해한다. 수백만 개의 정렬 레코드를 스트리밍할 때 발생하는 메모리 누수(Memory Leak)를 방지하는 제너레이터 패턴을 체득하고, 0-기반 반개방(0-based half-open) 좌표계 기반의 게놈 슬라이싱 질의(`bam.fetch()`)를 마스터한다. 단일 세포 핵심 태그(`CB`, `UB`, `GN`)를 $O(1)$ 시간 복잡도로 안전하게 파싱하는 예외 처리 아키텍처를 수립한다. `bam.pileup()` 엔진을 통해 전립선암 호르몬 치료 내성 핵심 핫스팟인 *AR* LBD 점돌연변이(*AR* p.T878A 등)의 염기별 정렬 깊이와 대립유전자 카운트를 단일 세포 해상도로 역공학(Reverse Engineering)하는 프로덕션급 파이썬 파이프라인을 구축한다.

---

## 1. `pysam`과 C HTSlib 아키텍처 및 Cython 브리지

생물정보학 분석에서 대용량 BAM 파일(수십 GB)을 순수 Python 문자열로 처리하면 극심한 속도 저하와 GIL(Global Interpreter Lock) 병목이 발생합니다. `pysam`은 전 세계 표준 C 라이브러리인 **HTSlib**의 C 구조체 및 함수를 Cython C-API 바인딩으로 직접 감싸(Wrap) 순수 C 수준의 I/O 처리 속도와 Python의 유연성을 동시에 제공합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        pysam의 내부 C-API 브리지 아키텍처                             │
│                                                                                        │
│   Python Layer :   read = next(bam_iter) ──────► read.get_tag("CB")                   │
│                          │                               │                             │
│                    [Cython Wrapper]                [C Tag Parser]                      │
│                          ▼                               ▼                             │
│   C HTSlib Layer:  sam_read1(fp, hdr, b) ──────► bam_aux_get(b, "CB")                 │
│                          │                               │                             │
│                    [C Structs]                           ▼                             │
│                    bam1_t (core alignment)      uint8_t* raw_tag                       │
│                    bam_hdr_t (header)                                                  │
│                    samFile (BGZF file stream)                                          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 핵심 C 구조체와 Python 객체 매핑
1. **`samFile` $\leftrightarrow$ `pysam.AlignmentFile`**: 압축된 BGZF 스트림 디스크 핸들러 및 파일 포인터.
2. **`bam_hdr_t` $\leftrightarrow$ `pysam.AlignmentHeader`**: 참조 서열 이름, 길이, 정렬 프로그램 정보 등을 보관하는 메타데이터 헤더.
3. **`bam1_t` $\leftrightarrow$ `pysam.AlignedSegment`**: 단일 리드의 플래그, 좌표, CIGAR, 시퀀스, 퀄리티, 보조 태그가 C 연속 메모리 블록(Flat memory block)으로 할당된 핵심 구조체.

---

## 2. 대용량 스트리밍과 이터레이터 메모리 수명주기 (Memory Lifecycle)

수억 개의 리드를 처리할 때 가장 빈번하게 발생하는 주니어 개발자의 실수는 `list(bam.fetch())`와 같이 수백만 개의 `pysam.AlignedSegment` 객체를 메모리 리스트에 한 번에 담는 것입니다.

### 2.1 C 포인터 재참조와 메모리 폭발 방지
`pysam.AlignmentFile.fetch()`는 레코드를 하나씩 읽어 반환하는 **지연 평가 이터레이터(Lazy Evaluation)**를 반환합니다. 각 반환 객체는 `bam_dup1`로 독립 복사된 `AlignedSegment`이며(버퍼를 재사용하지 않음), `list(bam.fetch())`가 메모리를 소모하는 이유는 이 객체들을 모두 보관하기 때문입니다.

```python
# [위험한 안티 패턴: Out-Of-Memory (OOM) 유발]
# 1천만 개 리드 * 리드당 객체 오버헤드 (~1KB) = 10GB RAM 즉시 고갈!
reads_cache = [read for read in bam.fetch("chrX", 67543000, 67731000)]

# [권장 모범 패턴: O(1) 메모리 스트리밍 파이프라인]
# 메모리 상에는 항상 현재 처리 중인 단 1개의 리드 객체만 상주
for read in bam.fetch("chrX", 67543000, 67731000):
    process_read(read)
```

### 2.2 컨텍스트 매니저(`with`)의 엄격한 자원 반환
BGZF 파일 포인터와 메모리 버퍼의 누수(File descriptor leak)를 차단하기 위해 반드시 Python 컨텍스트 매니저를 사용해야 합니다.

```python
import pysam

with pysam.AlignmentFile("data/sample.bam", "rb") as bam:
    header = bam.header.to_dict()
    # 안전한 스트리밍 처리 블록
# 블록 탈출 시 C 함수 sam_close()가 자동 호출되어 C 버퍼 완전 해제
```

---

## 3. 좌표계 변환 수학: SAM/IGV 1-기반 폐구간 vs pysam 0-기반 반개방 구간

생물정보학 파이프라인에서 가장 흔하게 발생하는 오프-바이-원 에러(Off-by-one bug)는 좌표계 불일치에서 비롯됩니다.

| 체계 (Format) | 시작 좌표 | 종료 좌표 | 표기법 | 예시 구간 ($100\text{bp}$) | 길이 계산 수식 |
|---|---|---|---|---|---|
| **SAM / VCF / IGV** | 1-based | 1-based (Inclusive) | $[start, end]$ | `chrX:101-200` | $L = end - start + 1$ |
| **BED / Python / pysam** | 0-based | 0-based (Exclusive) | $[start, end)$ | `chrX:100-200` | $L = end - start$ |

```
IGV / SAM (1-based closed) :     1    2    3    4    5    6    7    8
                                ├───┼───┼───┼───┼───┼───┼───┼───┤
Base Sequence              :     A    C    G    T    A    C    G    T
                                ├───┼───┼───┼───┼───┼───┼───┼───┤
pysam / BED (0-based half) :    0    1    2    3    4    5    6    7    8
```

- **좌표 변환 공식**:
  - SAM/IGV 좌표 $[S_{\text{sam}}, E_{\text{sam}}]$을 `pysam.fetch()` 인자로 전달할 때:
    $$start_{\text{pysam}} = S_{\text{sam}} - 1$$
    $$stop_{\text{pysam}} = E_{\text{sam}}$$

예를 들어, IGV에서 확인한 *AR* 엑손 8의 핫스팟 좌표가 `chrX:67,723,800-67,723,900`이라면, `pysam` 질의는 반드시 `bam.fetch("chrX", 67723799, 67723900)`으로 호출되어야 단 1개의 염기도 누락되지 않습니다.

---

## 4. 단일 세포 보조 태그의 고속 파싱 및 예외 처리 아키텍처

단일핵 전사체 BAM의 모든 리드가 `CB`나 `UB` 태그를 가지고 있는 것은 아닙니다. 비드 바코드 복원에 실패했거나 품질 점수가 낮은 리드는 태그가 생략되어 있습니다.

### 4.1 `read.has_tag()` vs `try-except KeyError` 성능 비교
존재하지 않는 태그를 `read.get_tag()`로 호출하면 Python `KeyError` 예외가 발생합니다. 수천만 리드를 순회할 때 예외 처리는 스택 트레이스 생성 오버헤드로 인해 약 10~20배의 성능 저하를 초래합니다.

```python
# [최적화 기법: has_tag 사전 검사 또는 get_tag 방어적 패턴]
def extract_cell_metadata(read: pysam.AlignedSegment):
    # C 수준에서 태그 존재 여부를 신속 확인
    if not (read.has_tag("CB") and read.has_tag("UB")):
        return None
    
    cell_barcode = read.get_tag("CB")  # str: e.g. "AAACCTGAGTCGTTTG-1"
    umi = read.get_tag("UB")           # str: e.g. "TAGGCTAACG"
    gene_name = read.get_tag("GN") if read.has_tag("GN") else "Unknown"
    
    return cell_barcode, umi, gene_name
```

### 4.2 태그 타입(Tag Type) 바이너리 규격
`pysam`은 SAM 명세의 타입 문자에 따라 C 데이터를 적절한 Python 자료형으로 자동 캐스팅합니다:
- `Z`: 널 종단 ASCII 문자열 $\to$ Python `str`
- `i`, `c`, `s`, `S`, `I`: 정수형 $\to$ Python `int`
- `f`: 단정밀도 부동소수점 $\to$ Python `float`
- `A`: 단일 문자 $\to$ Python `str` (길이 1)
- `B`: 바이트 또는 정수 배열 $\to$ Python `array` 또는 `list`

---

## 5. `bam.pileup()` 엔진을 활용한 단일 염기 깊이 및 돌연변이 카운팅

게놈 상의 특정 1개 염기 위치(유전자좌)에 몇 개의 리드가 쌓여 있는지, 각 염기(A, C, G, T)의 빈도가 얼마인지 계산하기 위해 `pileup()` 엔진을 사용합니다.

```
Reference Genome:  ... T  T  C  [ A ]  C  T  T ...  (Position: chrX:67,723,710, hg38)
Read 1:                T  T  C    A    C  T
Read 2:                   T  C  [ G ]  C  T  T         <── Mutant Base (Alt)
Read 3:             G  T  T  C    A    C  T
Read 4:                   T  C  [ G ]  C  T            <── Mutant Base (Alt)
                                  ▲
                       Pileup Column (Depth = 4: Ref A=2, Alt G=2)
```

### 5.1 `pileup()` 핵심 파라미터 제어
```python
pileup_iter = bam.pileup(
    contig="chrX",
    start=67723709,      # 0-based start (1-based 67,723,710)
    stop=67723710,       # 0-based stop (단 1개 염기 위치)
    truncate=True,       # 질의 구간 외의 파일업 컬럼을 강제 절단하여 정확히 해당 위치만 반환
    stepper="samtools",  # samtools 기본 필터링 규칙 적용
    min_base_quality=20, # Phred Score < 20 (Q20, 에러율 > 1%) 염기 자동 마스킹 (pysam 기본값은 13)
    min_mapping_quality=30, # MAPQ < 30 리드 배제
    max_depth=100000     # pysam 기본 max_depth=8000: 초과 깊이는 조용히 잘리므로 깊은 좌위(고발현 AR)는 상향
)
# samtools stepper는 orphans·DUP·SECONDARY·QCFAIL 리드를 기본으로 제외합니다.
```

### 5.2 `PileupRead` 객체의 세부 속성 추출
각 파일업 컬럼(`PileupColumn`)은 그 위치를 가로지르는 리드들의 목록(`PileupRead`)을 갖습니다:
- `pileupread.is_del`: 해당 위치에서 리드가 결실(`D`) 상태인지 여부.
- `pileupread.is_refskip`: 해당 위치가 RNA 스플라이싱 인트론(`N`) 구간인지 여부 (RNA-seq에서 필수 체크!).
- `pileupread.query_position`: 리드 서열(`read.query_sequence`) 내에서 현재 염기가 위치한 0-기반 인덱스.

---

## 6. 전립선암 임상 타깃: *AR* LBD 점돌연변이 단일 세포 마이닝 파이프라인

안드로겐 수용체(AR)의 리간드 결합 도메인(LBD)에는 항안드로겐 내성을 주는 점돌연변이가 알려져 있습니다. 엔잘루타마이드를 길항제에서 작용제로 전환시키는 대표 변이는 *F877L*(구 표기 F876L; Korpal 2013은 F876L로 표기)이며, *T878A*는 단독으로는 엔잘루타마이드 작동 효과가 약하고 F877L/T878A 이중 변이에서 작동 활성을 크게 증폭시킵니다(Korpal 2013; Mol Cancer Ther 2016). *L702H*는 엑손 8이 아니라 엑손 4에 있습니다. 아래는 T878A(엑손 8)를 예시로 한 UMI 카운팅 파이프라인입니다.
- **표적 돌연변이**: *AR* c.2632A>G (p.T878A)
- **게놈 좌표 (GRCh38)**: `chrX:67,723,710` (NM_000044.6 c.2632, 코돈 `ACT`의 첫 염기; hg38 chrX:67,723,878의 참조 염기는 `C`이므로 이전 좌표를 쓰면 참조가 불일치해 모든 UMI가 other로 분류된다. 코드에 참조 염기 검증을 둘 것)
- **참조 염기 (Ref)**: `A` / **변이 염기 (Alt)**: `G`

다음은 단일핵 전사체 BAM 파일로부터 단일 세포(`CB`) 및 분자(`UB`) 해상도로 *AR* T878A 변이를 호출하고, 세포별 Ref/Alt UMI 매트릭스를 구축하는 프로덕션급 Python 코드입니다.

```python
#!/usr/bin/env python3
"""
Single-Cell AR T878A Mutation Mining Pipeline using pysam
Author: Senior Bioinformatics Engineer
Target: AR p.T878A (chrX:67,723,710 A>G, GRCh38)
"""

from typing import Dict, Tuple, Set, Optional
from collections import defaultdict
import pysam
import pandas as pd
import numpy as np


def mine_single_cell_mutation(
    bam_path: str,
    contig: str = "chrX",
    pos_1based: int = 67723710,
    ref_base: str = "A",
    alt_base: str = "G",
    min_mapq: int = 30,
    min_baseq: int = 20,
    ref_fasta: Optional[str] = None,
    min_total_umi: int = 2,
    max_depth: int = 100000   # pysam pileup 기본 max_depth=8000: 깊은 좌위(고발현 AR)에서 조용히 잘림
) -> pd.DataFrame:
    """
    단일핵 BAM 파일로부터 특정 유전자좌의 세포별 Ref/Alt UMI 분자 수를 정확히 집계합니다.
    """
    pos_0based = pos_1based - 1

    # 좌표 오류 방지: 참조 FASTA가 주어지면 해당 위치의 참조 염기가 ref_base와 같은지 확인
    if ref_fasta is not None:
        with pysam.FastaFile(ref_fasta) as fa:
            actual = fa.fetch(contig, pos_0based, pos_0based + 1).upper()
        if actual != ref_base:
            raise ValueError(f"참조 염기 불일치: {contig}:{pos_1based} = {actual}, 기대값 {ref_base}")
    
    # 분자 중복 제거를 위한 딕셔너리: (Cell_Barcode, UMI) -> List of observed bases
    # PCR 증폭 산물이 여러 개일 경우 가장 품질이 높은 염기 또는 다수결을 취함
    molecule_observations: Dict[Tuple[str, str], list] = defaultdict(list)
    
    print(f"[*] Analyzing locus {contig}:{pos_1based} ({ref_base}>{alt_base}) in {bam_path}...")
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # 파일업 엔진 구동 (truncate=True로 정확한 1개 염기 타깃)
        pileup_columns = bam.pileup(
            contig=contig,
            start=pos_0based,
            stop=pos_0based + 1,
            truncate=True,
            min_base_quality=min_baseq,
            min_mapping_quality=min_mapq,
            stepper="samtools"
        )
        
        for col in pileup_columns:
            if col.reference_pos != pos_0based:
                continue
                
            for p_read in col.pileups:
                # 1. 인트론 스킵(N) 및 결실(D) 리드 필터링
                if p_read.is_del or p_read.is_refskip:
                    continue
                
                read: pysam.AlignedSegment = p_read.alignment
                
                # 2. PCR 중복 및 2차/보충 정렬 필터링 (Primary 고유 정렬만 수용)
                if read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue
                
                # 3. 단일 세포 바코드 및 UMI 존재 여부 검증
                if not (read.has_tag("CB") and read.has_tag("UB")):
                    continue
                
                cb = read.get_tag("CB")
                ub = read.get_tag("UB")
                
                # 4. 리드 서열 내 해당 위치의 염기 및 퀄리티 추출
                qpos = p_read.query_position
                if qpos is None:
                    continue
                    
                base = read.query_sequence[qpos].upper()
                base_qual = read.query_qualities[qpos]
                
                if base_qual >= min_baseq:
                    molecule_observations[(cb, ub)].append((base, base_qual))

    print(f"[*] Processed {len(molecule_observations)} unique single-cell molecules (CB-UB pairs).")
    
    # 분자(UMI) 수준의 합의(Consensus) 염기 판정 및 세포별 집계
    # 동일 UMI 내에서 염기 일치도가 80% 이상이어야 신뢰할 수 있는 단일 분자로 인정
    cell_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"ref_umi": 0, "alt_umi": 0, "other_umi": 0})
    
    for (cb, ub), obs in molecule_observations.items():
        bases = [b for b, q in obs]
        total_obs = len(bases)
        
        ref_matches = bases.count(ref_base)
        alt_matches = bases.count(alt_base)
        
        if ref_matches / total_obs >= 0.8:
            cell_counts[cb]["ref_umi"] += 1
        elif alt_matches / total_obs >= 0.8:
            cell_counts[cb]["alt_umi"] += 1
        else:
            cell_counts[cb]["other_umi"] += 1

    # DataFrame 변환 및 통계량 산출
    rows = []
    for cb, counts in cell_counts.items():
        ref_u = counts["ref_umi"]
        alt_u = counts["alt_umi"]
        tot_u = ref_u + alt_u
        
        vaf = (alt_u / tot_u) if tot_u > 0 else np.nan
        
        # 유전자형 분류 (최소 min_total_umi UMI 이상일 때만 콜링)
        # AR은 X 연관이라 남성 종양에서는 반접합(Hemizygous): Alt만 있으면 Hemizygous_Mut,
        # Ref/Alt 공존은 AR 증폭·혼합 클론 등을 시사. (A>G는 ADAR RNA 편집으로도 나타날 수 있음)
        if tot_u == 0:
            call = "DropOut"
        elif tot_u < min_total_umi:
            call = "LowDepth"
        elif alt_u >= 1 and ref_u == 0:
            call = "Hemizygous_Mut"
        elif alt_u >= 1 and ref_u >= 1:
            call = "Mixed_Ref_Alt"
        else:
            call = "WildType_Ref"
            
        rows.append({
            "cell_barcode": cb,
            "ref_umi": ref_u,
            "alt_umi": alt_u,
            "total_informative_umi": tot_u,
            "single_cell_vaf": vaf,
            "mutation_call": call
        })
        
    df = pd.DataFrame(rows)
    return df


if __name__ == "__main__":
    import sys
    # 실행 예시: python 06_mine_ar_mut.py sample.bam
    bam_file = sys.argv[1] if len(sys.argv) > 1 else "data/sample.bam"
    try:
        res_df = mine_single_cell_mutation(bam_file)
        print(res_df.head(10))
        print(f"[*] Mutation calls summary:\n{res_df['mutation_call'].value_counts()}")
    except FileNotFoundError:
        print(f"[!] Test file {bam_file} not found. Please provide valid BAM path.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [수학적 원리 규명 (Mathematical Rigor)]
UCSC 게놈 브라우저 및 IGV 상에서 전립선암 핵심 핫스팟 변이인 *FOXA1* (p.M253K)의 유전체 좌표가 `chr14:37,590,500` (1-based, GRCh38; **좌표 변환 산술 연습용 임의 좌표**이며 실제 FOXA1 M253 위치가 아님 — 실제 코돈 253은 음성 가닥으로 chr14:37,592,025-37,592,027, 37,590,500은 3′UTR. M253이 PCa 핵심 핫스팟인지도 확인하지 못함)으로 가정한다.
1) `pysam.AlignmentFile.fetch()` 및 `pysam.AlignmentFile.pileup()`을 호출할 때 전달해야 하는 0-기반 반개방(0-based half-open) 시작(`start`) 및 종료(`stop`) 정수 인자를 수학적으로 유도하시오.
2) 만약 주니어 분석가가 실수로 `bam.fetch("chr14", 37590500, 37590500)`으로 시작과 끝에 동일한 값을 전달했을 때, `pysam`이 반환하는 리드 집합의 결과를 설명하고, 구간의 길이 공식 $L = stop - start$ 관점에서 그 이유를 증명하시오.
3) `pileup()` 호출 시 `truncate=True` 옵션을 누락했을 때 발생하는 생물정보학적 계산 오류와 파이프라인의 오버헤드를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 0-기반 반개방 좌표 변환 유도**:
  - 1-기반 폐구간(Inclusive)에서 단일 염기 좌표는 $[37590500, 37590500]$입니다.
  - 0-기반 반개방(Half-open) 변환 공식:
    $$start_{\text{0-based}} = \text{pos}_{\text{1-based}} - 1 = 37,590,500 - 1 = 37,590,499$$
    $$stop_{\text{0-based}} = \text{pos}_{\text{1-based}} = 37,590,500$$
  - 따라서 전달할 좌표는 `start=37590499`, `stop=37590500`입니다.

- **2. 동일 좌표 전달 시의 결과 증명**:
  - 반개방 구간 $[a, b)$에서 구간 내에 포함되는 게놈 염기의 총 개수는 $L = b - a$입니다.
  - $start = 37590500$, $stop = 37590500$인 경우:
    $$L = 37590500 - 37590500 = 0$$
  - 길이가 0인 공집합(Empty set) 구간이 되므로, `pysam`은 단 1개의 리드도 반환하지 않고 즉시 종료되는 빈 제너레이터(Empty Generator)를 반환합니다.

- **3. `truncate=True` 누락 시의 치명적 결함**:
  - `pileup()` 엔진은 기본적으로 질의 구간과 단 1bp라도 겹치는 모든 리드들을 검색한 후, 그 리드들이 걸쳐 있는 가장 왼쪽 염기부터 가장 오른쪽 염기까지 **리드 전체의 모든 파일업 컬럼을 순회**합니다.
  - 예를 들어 리드 길이가 150bp이고 거대 인트론(100kb)을 가로지르는 스플라이스 리드가 존재한다면, `truncate=True`가 없을 경우 단 1개 염기를 질의했음에도 불구하고 앞뒤 수만 bp에 달하는 모든 유전체 위치의 `PileupColumn`을 순차적으로 생성하여 반환합니다.
  - 이는 파이프라인에 수천 배의 불필요한 루프 오버헤드를 야기하며, `col.reference_pos`를 검사하지 않고 첫 번째 컬럼을 대상 염기로 오인하면 엉뚱한 유전자좌의 깊이와 염기를 읽게 되는 심각한 계산 오류를 발생시킵니다.

</details>

---

### Q2. [알고리즘 구현 및 성능 최적화 (Algorithmic Logic)]
단일핵 시퀀싱 데이터에서 1,000만 개의 리드가 매핑된 BAM 파일을 대상으로 유전자별 단일 세포 발현량을 계산하는 스크립트를 작성하려 한다.
1) 다음 두 파이썬 코드 패턴 (A)와 (B)의 메모리 복잡도($O(N)$ vs $O(1)$)를 비교하고, 프로덕션 환경에서 발생할 수 있는 시스템 충돌(Out-Of-Memory) 원인을 분석하시오:
   - 패턴 A: `all_reads = [r for r in bam.fetch("chrX")]`
   - 패턴 B: `for r in bam.fetch("chrX"): process(r)`
2) `pysam`에서 `read.get_tag("CB")`를 호출할 때, `KeyError` 예외 처리(`try-except`) 블록을 수천만 번 실행하는 것보다 `read.has_tag("CB")`로 사전 검증하는 것이 전산학적으로 훨씬 빠른 이유를 Python 가상머신(CPython)의 바이트코드 및 예외 프레임워크 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 메모리 복잡도 비교 및 OOM 분석**:
  - **패턴 A**: 메모리 복잡도 $O(N)$.
    - 리스트 컴프리헨션은 염색체 X에 정렬된 수백만 개의 `AlignedSegment` Python 객체를 생성하여 하나의 리스트 배열에 모두 보관합니다.
    - C 구조체 `bam1_t` 외에도 Python 객체 헤더(PyObject header, 레퍼런스 카운트, 타입 포인터 등)로 인해 리드당 약 500B ~ 1KB의 RAM이 소모됩니다. 1천만 리드의 경우 5~10GB 이상의 메모리가 순간적으로 점유되어 OS의 OOM Killer에 의해 프로세스가 강제 종료(SIGKILL)됩니다.
  - **패턴 B**: 메모리 복잡도 $O(1)$.
    - 제너레이터(Generator) 방식으로, 디스크 BGZF 블록으로부터 한 번에 1개 레코드의 바이트만 메모리에 로드하고 루프 본문에서 소비한 후 다음 단계에서 동일 버퍼를 재사용하거나 가비지 컬렉터(GC)에 즉각 반환합니다. 리드 수가 수억 개로 증가해도 상주 메모리는 수십 MB 수준으로 일정하게 유지됩니다.

- **2. `has_tag()` 사전 검증이 예외 처리보다 빠른 이유**:
  - CPython에서 `try-except` 블록 자체는 비용이 거의 들지 않지만(Zero-cost exception handling), 실제로 예외(`KeyError`)가 발생하는 순간 극심한 오버헤드가 발생합니다.
  - 예외가 발생하면 Python 가상머신은 현재 스택 프레임을 중단하고, C 백트레이스(Traceback) 객체를 동적 힙에 할당하며, 호출 스택을 역추적(Unwinding)하는 고비용의 전산 작업을 수행합니다.
  - 단일핵 전사체에서 비드 바코드가 없는 리드는 수백만 개에 달하므로, 수백만 번의 `KeyError` 발생은 초당 처리량을 수십 배 떨어뜨립니다.
  - 반면 `read.has_tag("CB")`는 Cython을 통해 HTSlib의 C 함수 `bam_aux_get(b, "CB")`를 호출하여 C 포인터가 `NULL`인지 여부만 $O(1)$로 단순 비트 비교하므로 Python 예외 객체 생성 없이 극한의 속도로 필터링할 수 있습니다.

</details>

---

### Q3. [임상 케이스 분석 (Clinical Problem Solving)]
전립선암 환자 코호트에서 전이성 병변 조직의 단일핵 전사체 BAM 파일을 분석하여 Enzalutamide 내성을 유도하는 *AR* p.T878A (chrX:67,723,710 A>G) 돌연변이를 스크리닝하였다.
분석관의 스크립트 실행 결과, 3명의 환자 검체에서 다음과 같은 단일 세포 UMI 집계 패턴이 관찰되었다:
- Patient A (세포 500개): 각 세포당 $N_{\text{tot}} \ge 4$ UMI가 검출되었고, 모든 세포에서 $VAF = 1.0$ (Alt `G` UMI만 검출).
- Patient B (세포 450개): 각 세포당 $N_{\text{tot}} \ge 5$ UMI가 검출되었고, 모든 세포에서 $VAF \approx 0.5$ (Ref `A` 2~3개, Alt `G` 2~3개 검출).
- Patient C (세포 600개): 각 세포당 $N_{\text{tot}} = 1$ UMI만 검출되었으며, 300개 세포는 Ref `A` ($VAF=0$), 300개 세포는 Alt `G` ($VAF=1.0$)로 완벽히 양극화됨.
1) 환자 A와 환자 B의 *AR* 유전체 복제수 변이(CNA) 및 염색체 상태를 비교 분석하시오 (안드로겐 수용체는 염색체 X에 위치함).
2) 환자 C에서 관찰된 $VAF = 0$과 $VAF = 1$의 양극화 현상이 독립된 2개의 서브클론(정상 클론 vs 돌연변이 클론)을 의미하는지, 아니면 단일 세포 기술적 한계에 의한 착시인지 판정하고 그 분자생물학적/기술적 메커니즘을 증명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **1. 환자 A와 환자 B의 유전체 복제수(CNA) 상태 분석**:
  - 남성 게놈에서 염색체 X는 정상적으로 반접합체(Hemizygous, 복제수 1)입니다.
  - **Patient A**: 모든 세포에서 $VAF = 1.0$으로 단일 대립유전자(Alt `G`)만 발현됩니다. 이는 남성 반접합체 *AR* 대립유전자에 T878A 돌연변이가 발생한 전형적인 동형접합(Hemizygous mutant) 상태이거나, 돌연변이 대립유전자가 국소 유전자 증폭(Focal Amplification)된 클론입니다.
  - **Patient B**: 모든 세포에서 $VAF \approx 0.5$로 Ref와 Alt가 정확히 $1:1$ 비율로 공존합니다. 정상 남성에서는 결코 $VAF=0.5$가 나올 수 없습니다. 따라서 환자 B의 종양 세포는 **전립선암 진행 과정에서 *AR* 유전자좌가 복제(Duplication) 또는 증폭(Amplification)되어 최소 2개 이상의 대립유전자를 획득한 상태에서, 그 중 하나의 대립유전자에만 T878A 돌연변이가 발생한 이형증폭(Heterozygous/Altered amplification)** 상태임을 입증합니다.

- **2. 환자 C의 VAF 양극화 원인 판정 및 메커니즘 증명**:
  - **판정: 데이터만으로는 독립된 두 클론인지 단일 대립유전자 전사 버스팅(Monoallelic Transcriptional Bursting)과 극단적 UMI 저샘플링(Low Depth)에 의한 기술적 착시인지 구별할 수 없습니다(후자 가능성이 높지만 단정할 수 없음).**
  - **증명**:
    1. 환자 C의 모든 세포는 $N_{\text{tot}} = 1$ UMI만을 갖습니다. 단 1개의 분자만 포획되었을 때, 수학적으로 관측 가능한 VAF 값의 치역은 $k/1 \in \{0.0, 1.0\}$ 외에는 존재할 수 없습니다. $0.5$라는 중간값은 물리적으로 표현 불가능합니다.
    2. 진핵생물의 전사는 연속적이지 않고 On/Off 펄스로 일어나는 '전사 버스팅'을 따릅니다. 환자 C가 실제로는 환자 B처럼 두 개의 대립유전자(Ref 1개 + Alt 1개)를 가진 이형 세포군이라 할지라도, 용해 시점에 단 1개의 mRNA 분자만 포획되면 50%의 확률로 Ref가 잡히고($VAF=0$), 50%의 확률로 Alt가 잡혀($VAF=1$) 세포 집단이 인위적으로 양분됩니다.
    3. 따라서 $N_{\text{tot}} = 1$인 상태에서는 '두 서브클론'과 '이형 상태 + 저샘플링'이 모두 이 데이터와 양립하므로 **판정 불가(추가 근거 필요)**입니다. 세포들을 병합(Pseudo-bulk)하거나 베이지안 베타-이항 모델로 분석하고 DNA 수준 근거로 확인해야 합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Korpal M, Korn JM, Gao X, et al. An F876L mutation in androgen receptor confers genetic and phenotypic resistance to MDV3100 (enzalutamide). *Cancer Discov*. 2013;3(9):1030-43. PMID: 23842682. doi: 10.1158/2159-8290.CD-13-0142
2. Prekovic S, van Royen ME, Voet AR, et al. The Effect of F877L and T878A Mutations on Androgen Receptor Response to Enzalutamide. *Mol Cancer Ther*. 2016;15(7):1702-12. PMID: 27196756. doi: 10.1158/1535-7163.MCT-15-0892
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- pysam 0.24.1 (github.com/pysam-developers/pysam).
<!-- REFERENCES:END -->

# Chapter 03. snATAC 특화 품질 관리(QC) 지표와 염색질 무결성 판별

> **Cognitive Goal**: snATAC-seq의 품질을 결정하는 타협할 수 없는 2대 핵심 지표인 전사 시작점 농축 점수(TSS Enrichment Score)와 약 200 bp 주기의 뉴클레오솜 사다리(Nucleosome Banding Periodicity; 147 bp는 코어 DNA 길이)의 생물물리학적 원리를 체득하고, 고품질 단일 핵과 저품질 파편(Debris) 및 나체 DNA 오염을 수학적·통계적으로 완벽히 분리한다.

---

## 1. snATAC 품질 판정의 2대 절대 기준 (The Two Non-Negotiable Pillars)

snRNA-seq에서는 UMI 수, 유전자 검출 수, 그리고 세포질 잔여 미토콘드리아 비율(`pct_counts_mt`)을 통해 세포의 품질을 평가했습니다.

그러나 snATAC-seq은 DNA 수준의 분석이므로 RNA 품질 지표를 그대로 적용할 수 없습니다. 단순히 "단편 수(Total fragments)가 많다"고 해서 좋은 세포가 아닙니다. 세포막과 핵막이 터져 세포 밖으로 쏟아져 나온 **나체 게놈 DNA(Naked DNA debris)** 역시 Tn5 효소에 의해 수만 개의 단편으로 잘려 나갈 수 있기 때문입니다.

진정한 고품질 단일 핵을 감별하기 위해 후성유전학 분야에서는 다음 **2대 절대 품질 기준**을 필수로 적용합니다:

```
[snATAC 2대 핵심 품질 지표]
1. TSS Enrichment Score (신호 대 잡음비: S/N Ratio)
   ──► 열려 있어야 할 핵심 유전자 프로모터에 단편이 집중적으로 농축되었는가?
2. ~200 bp Nucleosome Banding Periodicity (염색질 구조 무결성)
   ──► 히스톤 8량체가 DNA를 규칙적으로 보호하고 있는 구조적 주기성이 유지되는가?
```

---

## 2. 전사 시작점 농축 점수 (TSS Enrichment Score)

**TSS (Transcription Start Site) Enrichment Score**는 게놈 전체의 무작위 비특이적 절단 신호 대비, 활성 프로모터의 뉴클레오솜 결핍 영역(NDR)에 떨어진 진정한 생물학적 신호의 강도를 나타내는 **신호 대 잡음비(Signal-to-Noise Ratio)**입니다.

### 2.1 수학적 정의 및 계산 공식
인간 유전체 내 검증된 하우스키핑 및 표준 유전자 수만 개의 전사시작점(TSS, $0 \text{ bp}$)을 기준으로, 양쪽 $\pm 2,000 \text{ bp}$ 영역의 단편 삽입 중심점(Tn5 insertion sites)의 빈도를 정규화하여 집계합니다:

$$\text{TSS Enrichment Score} = \frac{\text{Mean Insertion Signal at TSS Center } [-50 \text{ bp}, +50 \text{ bp}]}{\text{Mean Background Signal in Flanking Regions } [-2000 \text{ bp}, -1900 \text{ bp}] \cup [+1900 \text{ bp}, +2000 \text{ bp}]}$$

위 정의는 ArchR의 구현(TSS 중심 101 bp 창을 양쪽 끝 1,900~2,000 bp의 100 bp 구간으로 정규화)에 따른 것입니다. Signac은 ENCODE의 TSS enrichment 정의를 따르며 창 크기·정규화 세부가 도구마다 조금씩 다르므로, 도구가 달라지면 같은 임계값을 그대로 옮겨 쓰지 않습니다.

```
Relative Signal
     ▲
 12x │                ┌──┐ [TSS Center Peak: 고농축 열린 프로모터]
  8x │               ┌┘  └┐
  4x │              ┌┘    └┐
  1x │══════════════┘      └══════════════ [Flanking Background Baseline = 1.0]
     └─────────────┬────────┬─────────────►
                -2000 bp   TSS   +2000 bp
```

### 2.2 TSS 농축 점수의 품질 판정 기준치 (경험적 가이드; 09장 등급표와 동일)
등급 경계는 ENCODE나 10x가 공식 제시한 표준이 아니라, 도구 기본값(ArchR `minTSS=4`·`minFragments=1000`, Signac PBMC 예시 `TSS.enrichment > 4`)을 참고한 경험적 가이드이며 조직·키트·참조 TSS 세트에 따라 달라집니다.
- **Ideal (최우수)**: $\text{TSS Score} \ge 10.0$
  - 핵막이 잘 보존되어 프로모터에 Tn5 삽입이 고도로 집중된 검체.
- **Good (우수)**: $6.0 \le \text{TSS Score} < 10.0$
  - 표준 분석에 적합한 범위.
- **Acceptable (통과)**: $4.0 \le \text{TSS Score} < 6.0$
  - 다기관 동결 생검, 골전이 검체 등 섬유화/석회화가 심한 조직에서 흔히 얻는 하한 범위. 피크 필터링을 엄격히 적용.
- **Poor / Failed (실패/폐기)**: $\text{TSS Score} < 4.0$
  - 세포질 및 핵막 파열로 인해 프로모터 특이적 신호가 사라지고 게놈 전반에 균일한 배경 잡음이 덮인 검체. **다운스트림 분석 투입 불가**.

### 2.3 TSS 외 보완 QC 지표: nucleosome signal, FRiP, blacklist ratio
TSS 농축 점수와 함께 Signac·ArchR 워크플로우에서 표준적으로 계산하는 세 지표입니다. 임계값은 데이터셋별 예시이며 분포를 보고 정합니다.

| 지표 | 정의 | Signac PBMC 예시 임계값 |
| :--- | :--- | :--- |
| **nucleosome signal** | 세포별 모노뉴클레오솜 단편(147~294 bp) 수 / 뉴클레오솜 결핍(NFR) 단편(< 147 bp) 수의 비. 과도하면 NFR 신호 부족을 의미 | `nucleosome_signal < 4` |
| **FRiP** (`pct_reads_in_peaks`) | 피크 영역에 떨어진 단편 수 / 세포의 전체 유효 단편 수 × 100 | `pct_reads_in_peaks > 40` |
| **blacklist ratio** | ENCODE blacklist 영역에 떨어진 카운트 / 피크 영역 카운트(비율). 높으면 인공 신호 | `blacklist_ratio < 0.01` |

```r
# Signac (PBMC 튜토리얼 방식; 임계값은 해당 데이터 기준 예시)
pbmc <- NucleosomeSignal(pbmc)    # nucleosome_signal
pbmc <- TSSEnrichment(pbmc)       # TSS.enrichment
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- FractionCountsInRegion(pbmc, assay = "peaks", regions = blacklist_regions)
pbmc <- subset(pbmc, subset = pct_reads_in_peaks > 40 & blacklist_ratio < 0.01 &
                              nucleosome_signal < 4 & TSS.enrichment > 4)
```

SnapATAC2에서는 `snap.metrics.frag_size_distr`(단편 길이 분포)와 `snap.metrics.frip`(영역 dict를 받아 세포별 FRiP 계산)을 사용할 수 있습니다.

---

## 3. 뉴클레오솜 사다리 주기성 (Nucleosome Banding Periodicity)

단일 핵 ATAC 라이브러리의 절단 단편 길이(Fragment Length)를 $0 \text{ bp}$부터 $800 \text{ bp}$까지 플롯하면, 염색질의 물리적 구조를 대변하는 **약 200 bp 간격(모노뉴클레오솜 ≈180–220, 다이뉴클레오솜 ≈380–420)의 주기적 사다리 패턴(Banding Ladder)**이 나타납니다. (147 bp는 히스톤 코어에 감긴 DNA 길이이며 사다리의 주기가 아니라 NFR/모노 경계 컷오프로 쓰인다.)

```
Fragment Density
     ▲
     │  [1. NFR: < 100 bp]
     │       ▲
     │      ╱ ╲
     │     ╱   ╲      [2. Mononucleosome: ~200 bp]
     │    ╱     ╲           ▲
     │   ╱       ╲         ╱ ╲      [3. Dinucleosome: ~400 bp]
     │  ╱         ╲       ╱   ╲           ▲
     │ ╱           ╲─────╱     ╲─────────╱ ╲─────────
     └────────────────────────────────────────────────► Fragment Length (bp)
       0         100     147   200       300 400   (147 = 코어 DNA 길이, 주기 아님)
```

### 3.1 단편 길이 히스토그램의 3대 핵심 영역
1. **뉴클레오솜 결핍 단편 (NFR: Nucleosome-Free Region, $< 100 \text{ bp}$)**:
   - 두 개의 인접한 히스톤 단백질 사이의 벌어진 링커 DNA나 활성 인핸서의 전사인자 결합 부위 내부에서 두 번의 Tn5 절단이 일어난 경우.
   - **전사인자 모티프 발굴 및 피크 콜링 신호의 핵심 원천**입니다.
2. **단일 뉴클레오솜 단편 (Mononucleosome Peak, $\sim 180 - 220 \text{ bp}$)**:
   - 히스톤 8량체 코어(147 bp)를 둘러싼 양쪽 링커 DNA 부위를 한 번씩 절단하여, 1개의 뉴클레오솜을 온전히 포함하는 길이의 단편.
3. **이중 뉴클레오솜 단편 (Dinucleosome Peak, $\sim 380 - 420 \text{ bp}$)**:
   - 2개의 연속된 뉴클레오솜을 포함하는 단편 ($147 \times 2 + \text{Linkers} \approx 400 \text{ bp}$).

### 3.2 주기성 붕괴에 따른 라이브러리 병리 진단
- **정상 (Intact Chromatin)**: 명확한 NFR 피크와 함께 200 bp, 400 bp의 감쇠하는 파동형 뉴클레오솜 피크가 선명히 관찰됨.
- **과절단 (Over-tagmentation)**: Tn5 효소 과다 투입으로 인해 뉴클레오솜 피크가 완전히 소멸하고 $80 \text{ bp}$ 이하의 초미세 단편만 비정상적으로 치솟음.
- **나체 DNA 오염 (Naked DNA Degradation)**: 핵막 붕괴 및 히스톤 박리로 인해 뉴클레오솜 사다리 주기가 완전히 평탄화(Flat)되어, 단순 무작위 분해 곡선인 **지수적 단조 감소 곡선(Monotonic exponential decay)**을 나타냄.

---

## 4. 미토콘드리아 DNA (mtDNA) 단편 오염 및 필터링

미토콘드리아는 자체적인 원형 이중가닥 DNA(약 16.5 kb)를 세포당 수백~수천 카피 보유하고 있습니다.

### 4.1 미토콘드리아 DNA의 특성: 보호막 없는 나체 DNA
- 미토콘드리아 DNA는 **진핵세포의 핵 뉴클레오솜(코어 히스톤)으로 포장되어 있지 않습니다**. (히스톤이 없을 뿐 TFAM 등으로 핵양체(nucleoid)에 포장되어 있어 완전한 '나체 DNA'는 아니지만, 핵 염색질에 비해 Tn5 접근성이 높습니다.)
- 따라서 mtDNA는 핵 게놈보다 Tn5에 훨씬 잘 절단되는 표적이 됩니다.
- 단일 핵 분리 과정에서 세포질이 완벽히 제거되지 않은 채 Tn5 효소에 노출되면, 고활성 Tn5 효소가 미토콘드리아 DNA를 무차별적으로 과절단(Hyper-tagmentation)합니다.
- 이로 인해 시퀀싱 리드의 상당 부분이 미토콘드리아 DNA에 낭비 매핑되어 시퀀싱 비용 손실과 게놈 피크 검출력 저하를 초래합니다. 핵 분리가 잘 된 snATAC에서는 chrM 비율이 경험적으로 수 % 이하이고 ($f_{\text{mt}} > 5\%$ 경고, $> 20\%$ 심각; 09장 기준과 동일), 세포 전체를 쓰는 scATAC이나 세포질 오염이 큰 프로토콜에서는 이보다 훨씬 높게 나타날 수 있습니다(조직·프로토콜 의존적인 경험적 가이드).

### 4.2 필터링 전략
업스트림 및 QC 단계에서 염색체 `chrM`에 정렬된 모든 단편을 집계한 뒤, 각 바코드별 `pct_counts_mitochondrial`을 계산하여 **단일 핵 ATAC 데이터셋에서 chrM 유래 단편을 전량 제거**한 후 게놈 상염색체 분석을 진행해야 합니다.

---

## 5. Python 기반 snATAC 2차원 품질 필터링 파이프라인

단일 차원(예: 단편 수만으로 필터링)의 한계를 극복하기 위해, **$\log_{10}(\text{Unique Fragments})$ vs $\text{TSS Enrichment Score}$ 2차원 산점도**를 기반으로 진정한 세포를 선별합니다:

```python
import snapatac2 as snap
import muon as mu
import matplotlib.pyplot as plt

# 1. Fragment 파일로부터 AnnData/SnapATAC2 객체 로드
data = snap.pp.import_fragments(
    fragment_file="outs/atac_fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38,
    min_num_fragments=500,
    sorted_by_barcode=False
)

# 2. TSS Enrichment Score 계산
snap.metrics.tsse(data, snap.genome.hg38)

# 3. 2차원 필터링 기준 수립
# - 예시 임계값(조직·키트 의존): TSS Score >= 6.0 AND 1,500 <= Unique Fragments <= 80,000
#   (참고: snapatac2 filter_cells 기본값 min_counts=1000, min_tsse=5.0 / ArchR 기본 minTSS=4, minFrags=1000)
snap.pp.filter_cells(
    data, 
    min_counts=1500, 
    max_counts=80000, 
    min_tsse=6.0
)

print(f"Filtered High-Quality Nuclei: {data.n_obs}")
```

```
TSS Score
     ▲
  15 │        ┌─────────────────────────┐
     │        │ ★ High-Quality Nuclei   │
  10 │        │ (TSS > 6, Frags > 1500) │
     │        └─────────────────────────┘
   6 ├───────────────────────────────────── [TSS Threshold]
     │    Low TSS / High Noise
   2 │    (Naked DNA / Broken Nuclei)
     └────────┬─────────────────────────►
             10^3   10^4   10^5   Unique Fragments
```

---

## 6. 결론 및 다음 단계 연결

TSS Enrichment Score($\ge 6$)와 뉴클레오솜 주기성을 모두 만족하는 고품질 단일 핵만을 선별함으로써, 전립선암 조직 파쇄 과정에서 발생한 나체 DNA 파편과 세포질 잔해물을 효과적으로 걸러냈습니다.

그러나 단일 핵 ATAC 데이터는 세포 1개당 두 벌의 대립유전자(Diploid)만 존재하므로 극도로 희소(Sparse)합니다. 다음 장 [Chapter 04. 피크 콜링과 LSI 차원 축소 기하학](04_PEAK_CALLING_AND_LSI.md)에서는 세포 집단을 슈도벌크(Pseudobulk)로 묶어 MACS3로 고감도 피크를 콜링하고, TF-IDF와 SVD를 결합한 Latent Semantic Indexing(LSI)으로 크로마틴 지형을 차원 축소하는 수학적 원리를 학습합니다.

---

## 7. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
snATAC-seq의 핵심 품질 관리 지표인 전사 시작점 농축 점수(TSS Enrichment Score)의 수학적 계산식을 서술하고, 온전하게 분리된 단일 핵과 파열된 세포질/핵막 파편(Debris) 사이에서 TSS 점수의 극단적인 차이가 발생하는 후성유전학적 메커니즘을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **TSS Enrichment Score 수학 공식**:
   $$\text{TSS Score} = \frac{\text{Mean Normalized Insertions at TSS Center } [-50 \text{ bp}, +50 \text{ bp}]}{\text{Mean Background Insertions at Flanking Regions } [-2000 \text{ bp}, -1900 \text{ bp}] \cup [+1900 \text{ bp}, +2000 \text{ bp}]}$$
   즉, 전사시작점 정중앙($\pm 50 \text{ bp}$)의 Tn5 절단 빈도를 프로모터 양 외곽($\pm 1,900 \sim 2,000 \text{ bp}$)의 배경 기저 절단 빈도로 나눈 비율입니다.
2. **품질 차이의 후성유전학적 메커니즘**:
   - **온전한 핵**: 살아있거나 정상 동결된 핵 내부에서는 필수 유전자의 프로모터 영역이 뉴클레오솜이 비어 있는 열린 염색질(NDR) 상태를 유지합니다. 따라서 크기가 큰 Tn5 효소는 응축된 게놈 영역을 피해 프로모터 중심부에 집중적으로 결합하여 높은 TSS 피크를 형성합니다 (위 등급표 기준 우수 이상, $\text{TSS} \ge 6$).
   - **파열된 파편(Debris) 및 나체 DNA**: 핵막이 터지거나 히스톤 단백질이 박리되어 게놈 DNA가 노출되면, 히스톤에 의한 물리적 차단막이 소실되어 Tn5가 프로모터뿐 아니라 비코딩 영역, 유전자 간 영역 등 게놈 전체를 무차별적이고 균일하게 절단합니다. 이로 인해 분모인 배경 신호가 급증하여 TSS Score가 $1 \sim 3$ 수준으로 폭락하게 됩니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
snATAC-seq 라이브러리의 절단 단편 길이(Fragment Length) 분포 곡선에서 약 200 bp 간격의 주기적인 뉴클레오솜 래더(Ladder) 패턴이 관찰되는 생물물리학적 이유를 설명하고, 라이브러리 제작 과정에서 효소 과절단(Over-tagmentation)이 일어났을 때 히스토그램 형태가 어떻게 변하는지 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **뉴클레오솜 래더(~200 bp 주기) 형성의 생물물리학적 원리**:
   - 진핵세포 DNA는 코어 히스톤 8량체를 약 147 bp 길이로 1.65바퀴 감싸고 있습니다.
   - Tn5 트랜스포세이스는 히스톤 단백질과 결합된 147 bp DNA 내부로는 침투하여 자를 수 없으며, 뉴클레오솜과 뉴클레오솜 사이의 노출된 링커 DNA(Linker DNA)만을 선택적으로 절단합니다.
   - 따라서 두 번의 절단 이벤트가 일어날 때:
     - 링커 DNA 내부에서 두 번 잘리면 뉴클레오솜 결핍 단편(NFR, $< 100 \text{ bp}$)이 생성됩니다.
     - 1개의 뉴클레오솜 양쪽 링커에서 잘리면 $147 \text{ bp} + \text{링커} \approx 200 \text{ bp}$ (모노뉴클레오솜)가 생성됩니다.
     - 2개의 뉴클레오솜을 건너뛰어 잘리면 $\sim 400 \text{ bp}$ (다이뉴클레오솜)가 생성됩니다.
   - 결과적으로 약 200 bp(모노 ≈180–220, 다이 ≈380–420)의 정수배에 가까운 주기적 파동 피크가 관찰됩니다. 147 bp는 코어 DNA 길이이고 링커가 더해진 값이 주기이므로 '147의 정수배'가 아닙니다.
2. **과절단(Over-tagmentation) 시의 히스토그램 변화**:
   - 반응 시간 과다 또는 Tn5 효소 농도가 과도하면, 링커 DNA뿐 아니라 뉴클레오솜 경계면까지 효소가 과도하게 침투하여 잘게 부숩니다.
   - 200 bp 이상의 모노/다이뉴클레오솜 피크가 완전히 소멸하고, $50 \sim 80 \text{ bp}$ 이하의 초미세 단편만이 비정상적으로 높게 치솟는 단일 첨두형 분포로 왜곡됩니다.

</details>

---

### [문제 3. 데이터 품질 트러블슈팅 사례 연구 (Case Study)]
**[임상 데이터 분석 시나리오]**
연구원이 mCRPC 환자 8명의 생검 조직 10x Multiome 데이터를 처리하면서, scRNA-seq의 통념에 익숙하여 snATAC 파이프라인에서도 단순히 "단편 수가 많은 세포가 좋은 세포"라고 판단하고 `n_fragments > 2,000` 조건만을 단일 필터로 적용하여 세포를 선별했습니다.
그 결과, 세포 수는 12,000개로 넉넉하게 확보되었으나, 이후 MACS3 피크 콜링 결과 유의한 피크가 거의 검출되지 않았고, UMAP 시각화 시 세포들이 유형별 클러스터를 형성하지 못하고 UMAP 중앙에 거대한 무정형의 노이즈 덩어리(Continuous Noise Smear)로 뭉개져 버렸습니다.

**[질문]**
1. 이 분석 실패의 근본 원인을 snATAC 품질 지표(TSS Score 및 Fragment Length) 관점에서 진단하시오.
2. 이 데이터셋에서 진정한 악성 종양 세포 및 면역세포를 구제하기 위해 즉각 적용해야 할 **2차원 결합 QC 필터링 파라미터와 진단 플롯 생성 절차**를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **실패 원인 진단**:
   - **나체 DNA 및 파열된 세포질 잔해물(Ghost Debris)의 대량 유입**: 단편 수(`n_fragments`)가 많다는 것은 세포가 온전하다는 증거가 아닙니다. 세포막과 핵막이 터져 세포질이나 용액 상에 떠돌던 나체 게놈 DNA 덩어리(Naked DNA clump)가 단일 액적에 들어가면 수만 개의 무작위 단편을 생성합니다.
   - 이 바코드들은 단편 수는 높지만 TSS Enrichment Score가 $1 \sim 2$ 수준으로 바닥이며, 뉴클레오솜 주기성이 없는 순수 노이즈입니다.
   - 연구원이 TSS Score를 검사하지 않고 단편 수만으로 필터링했기 때문에, 유효 세포보다 수배 이상 많은 가짜 파편 바코드가 데이터셋에 대량 잔류하여 게놈 전반의 신호 대 잡음비를 파괴하고 UMAP 상에서 무정형 노이즈 스미어를 형성한 것입니다.
2. **구제 및 2차원 결합 QC 해결책**:
   - **진단 플롯 생성**: $X$축을 $\log_{10}(\text{n\_fragments})$, $Y$축을 $\text{TSS Enrichment Score}$로 하는 2차원 육각형 빈 산점도(Hexbin scatter plot)를 생성하여 세포 집단과 노이즈 집단의 바이모달(Bimodal) 분리 지점을 시각적으로 확인합니다.
   - **엄격한 2차원 교차 필터링 적용**:
     $$\text{TSS Enrichment Score} \ge 6.0 \quad \text{AND} \quad 1,500 \le \text{n\_fragments} \le 80,000$$
     TSS 점수가 6.0 미만인 모든 저순도 바코드를 전량 폐기합니다.
   - **뉴클레오솜 사다리 검증**: 필터링 통과 세포들만을 대상으로 Fragment length distribution을 재플롯하여 모노/다이뉴클레오솜 주기성이 회복되었는지 최종 검증한 후 피크 콜링을 재수행합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 검수 근거 문헌
1. Buenrostro JD, Giresi PG, Zaba LC, et al. Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nat Methods*. 2013;10(12):1213-8. PMID: 24097267. doi: 10.1038/nmeth.2688
2. Granja JM, Corces MR, Pierce SE, et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nat Genet*. 2021;53(3):403-411. PMID: 33633365. doi: 10.1038/s41588-021-00790-6
3. Stuart T, Srivastava A, Madad S, et al. Single-cell chromatin state analysis with Signac. *Nat Methods*. 2021;18(11):1333-1341. PMID: 34725479. doi: 10.1038/s41592-021-01282-5
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- muon 0.2.0 (github.com/scverse/muon). 원 논문: Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. *Genome Biol*. 2022;23(1):42. PMID: 35105358. doi: 10.1186/s13059-021-02577-8
- SnapATAC2 2.10.0 (github.com/scverse/SnapATAC2). 원 논문: Zhang K, Zemke NR, Armand EJ, et al. A fast, scalable and versatile tool for analysis of single-cell omics data. *Nat Methods*. 2024;21(2):217-227. PMID: 38191932. doi: 10.1038/s41592-023-02139-9
<!-- REFERENCES:END -->

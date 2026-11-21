# Chapter 07. 전사 인자 풋프린팅 및 단일 핵 모티프 편차 분석 (ChromVAR & TOBIAS)

> **Cognitive Goal**: 전사 인자(TF)가 DNA 모티프에 물리적으로 결합할 때 발생하는 Tn5 트랜스포사제의 입체 장애(Steric Hindrance)와 풋프린트(Footprint)의 생물물리학적 원리를 이해하고, Tn5 고유의 염기서열 삽입 편향(Sequence Bias)을 수학적으로 보정하는 알고리즘을 체계화한다. 아울러 단일 핵 단위에서 전사 인자 패밀리별 활성도를 정량화하는 `ChromVAR` 편차(Deviation) 점수화 모델과 전립선암 파이오니어 전사 인자(*FOXA1*, *HOXB13*)의 염색질 개방 기전을 마스터한다.

---

## 1. 전사 인자 풋프린팅(Footprinting)의 생물물리학적 원리

ATAC-seq의 기본 원리는 Tn5 트랜스포사제가 닫힌 염색질(히스톤에 감긴 DNA)을 피해 열린 염색질 부위에 무작위로 침투하여 시퀀싱 어댑터를 삽입(Tagmentation)하는 것입니다. 그러나 열린 염색질 영역 내부에서도 **전사 인자(Transcription Factor, TF) 단백질이 자신의 동족 DNA 모티프(Cognate Motif, 10~20 bp)에 강력하게 결합해 있다면**, 거대한 전사 인자 단백질의 물리적 질량이 Tn5의 접근을 물리적으로 차단하는 **입체 장애(Steric Hindrance)**가 발생합니다.

```
                  [전사 인자 단백질 (예: FOXA1 / AR)]
                  ┌─────────────────────────────────┐
                  │        Transcription Factor     │  (물리적 차폐 장벽: ~15 bp)
                  └────────────────┬────────────────┘
                                   ▼
[열린 염색질 (Flank)]     [동족 결합 모티프 (Core)]     [열린 염색질 (Flank)]
  Tn5 침투 활발 ──►       Tn5 절단 완전 차단!        ◄── Tn5 침투 활발
 (High Cleavage)            (Zero Cleavage)             (High Cleavage)
══════════════════════╤═════════════════════╤══════════════════════
                      │  T A A A C A A A G  │
                      └─────────────────────┘
                                   │
                                   ▼
 [Tn5 절단 빈도 프로파일 (Cleavage Count Profile at 1 bp Resolution)]
 절단수
   ▲          ▲▲▲                                   ▲▲▲
   │         ▲▲▲▲▲                                 ▲▲▲▲▲
   │        ▲▲▲▲▲▲▲                               ▲▲▲▲▲▲▲
   │       ▲▲▲▲▲▲▲▲▲                             ▲▲▲▲▲▲▲▲▲
   │      ▲▲▲▲▲▲▲▲▲▲▲                           ▲▲▲▲▲▲▲▲▲▲▲
 0 ┼─────┴───────────┴──────▼▼▼▼▼▼▼▼▼▼▼▼▼──────┴───────────┴────► 게놈 좌표 (bp)
                         [깊은 계곡 (Footprint Valley)]
```

### 1.1 1 bp 해상도의 풋프린트 계곡 (Footprint Valley)
- **플랭킹 피크(Flanking Shoulders)**: 결합 부위 양옆의 노출된(열린) DNA는 Tn5 절단 빈도가 상대적으로 높게 나타납니다.
- **코어 계곡(Core Footprint Valley)**: 전사 인자가 직접 덮고 있는 모티프 중심부(Core 10~15 bp)는 절단 카운트가 0에 가깝게 푹 꺼지는 깊은 계곡(Valley/Dip)을 형성합니다.
- **생물학적 판정**: 이 계곡의 깊이(Footprint Depth)는 해당 전사 인자의 결합 점유율을 가늠하는 **정성적 지표**입니다. 정비례 관계가 아니며 TF별 DNA 체류 시간, Tn5의 서열 편향, 시퀀싱 깊이에 의존합니다(집계 신호에서 풋프린트가 검출되는 TF는 50% 수준이고 개별 TFBS 수준의 지지는 더 낮다고 보고됨; Bentsen et al. 2020).

---

## 2. Tn5 시퀀싱 서열 편향(Sequence Bias)과 수학적 디바이아싱 (TOBIAS)

풋프린팅 분석에서 가장 치명적인 함정은, **Tn5 트랜스포사제 자체가 특정 DNA 염기서열(특정 서열 문맥, 이중뉴클레오타이드 수준의 편향이 보고됨)을 선호하여 비대칭적으로 삽입되는 내인성 서열 편향(Intrinsic Sequence Preference)**을 지닌다는 점입니다.

```
[원시 절단 신호 (Raw Cuts)]    = [진정한 TF 풋프린트]  +  [Tn5 효소 서열 선호도 노이즈]
                                          │
                                          ▼
           [위양성 위험: 모티프 서열 자체가 Tn5 비선호 염기일 경우,
            단백질이 붙지 않았음에도 계곡(Valley)이 형성되어 TF 결합으로 오인!]
```

### 2.1 TOBIAS (ATACorrect) 알고리즘의 서열 보정 수식
`TOBIAS`나 `HINT-ATAC` 같은 차세대 풋프린팅 도구는 게놈 전체의 무단백질 시험관 내(In vitro) naked DNA Tn5 절단 데이터 또는 전사 인자가 결합하지 않는 비활성 영역의 리드를 기반으로 **Tn5 서열 편향 모델**을 구축합니다(TOBIAS ATACorrect는 이중뉴클레오타이드 가중치 행렬(dinucleotide weight matrix, DWM) 기반; Bentsen et al., Nat Commun 2020, PMID 32848148):

1. **예측 배경 절단율 ($E_i$) 모델링**:
   - 게놈 위치 $i$ 주변 염기서열 $S_i$에 대해, DWM으로부터 서열 편향에 따른 기대 절단 확률 $p(S_i)$를 도출하고 국소 신호량 $N_{\text{local}}$에 곱하여 기대 절단 수(개념식)를 얻습니다:
   $$E_i = N_{\text{local}} \cdot p(S_i)$$
2. **보정 절단 점수 (Corrected Cleavage Score, $C_i$) 도출**:
   - ATACorrect는 관측 신호에서 기대(편향) 신호를 **차감**하여 보정 신호를 만듭니다(Bentsen 2020: expected cuts를 uncorrected signal에서 subtract). 로그 비율이 아닙니다:
   $$C_i = O_i - E_i$$
3. **풋프린트 점수 (Footprint Score, FPS)**:
   - TOBIAS ScoreBigwig(`--score footprint`)는 보정 신호에서 좌우 플랭크의 **양의 신호 평균**과 중앙(footprint) 구간의 **음의 신호 평균**의 차, $\text{FPS}=\overline{\text{flank}^{+}}-\overline{\text{core}^{-}}$를 계산하고, 여러 플랭크/코어 너비 조합 중 최댓값을 취합니다(TOBIAS 소스 `tobias_footprint_array`). 단순한 '플랭크 평균 − 코어 최솟값'이 아닙니다.

---

## 3. 전립선암의 파이오니어 전사 인자: *FOXA1*과 *HOXB13*의 크로마틴 개방

전립선암의 발암 및 계통 가소성 전이 과정에서 가장 중요한 후성유전학적 개척자는 **파이오니어 전사 인자(Pioneer Transcription Factor)**인 **FOXA1**과 **HOXB13**입니다.

```
[닫힌 응축 염색질 (Heterochromatin)]
  히스톤 8량체에 단단히 감겨 일반 전사 인자(AR 등)가 결합 불가능한 상태
                      │
                      ▼
[Step 1. FOXA1 단백질의 파이오니어 침투]
  - Winged-helix 도메인이 히스톤 코어를 깨지 않고 누클레오솜 외곽에 직접 안착
  - SWI/SNF (BAF) 크로마틴 리모델러 및 p300/CBP 히스톤 아세틸화 효소 동원
                      │
                      ▼
[Step 2. 국소 염색질 개방 및 H3K27ac 마킹]
  - 히스톤이 물리적으로 밀려나며 광범위한 열린 염색질(Open Chromatin Peak) 형성
                      │
                      ▼
[Step 3. 안드로겐 수용체(AR) 2차 결합 및 전사 폭발]
  - 열린 공간으로 활성화된 AR 이합체가 결합하여 KLK3, TMPRSS2 전사 가동!
```

| 전사 인자 (TF) | 구조적 계열 (Family) | 파이오니어 기능 및 전립선암 특이적 역할 | Multiome 풋프린트 특성 |
|:---|:---|:---|:---|
| **FOXA1** | Forkhead / Winged-helix | 누클레오솜 침투 마스터 스위치; *AR* 및 *HOXB13* 결합 부위 사전 개방; 점돌연변이 시 결합 선호도 변경 | 뚜렷하고 대칭적인 ~12 bp 코어 풋프린트 계곡 형성 |
| **HOXB13** | Homeodomain | 전립선 상피 계통 결정 인자(파이오니어로 분류하는 것은 논쟁적이며 FOXA1과 협동하는 인자로 보는 문헌이 많음); G84E 생식세포 변이 시 전립선암 발병률 급증; FOXA1과 협동하여 인핸서 규정 | ~10 bp의 정밀한 홈 저해 신호 |
| **AR** (Androgen Receptor) | Nuclear Receptor | 호르몬 의존성 핵 수용체; 전립선 세포에서 많은 AR 결합 인핸서가 FOXA1·HOXB13이 사전에 열어둔 염색질 맥락에서 결합하지만 이는 맥락 의존적이며 FOXA1만이 유일한 결정 요인은 아님 | 리간드 결합 시에만 15 bp 인버티드 리피트 풋프린트 출현 |
| **ASCL1** | Basic Helix-Loop-Helix (bHLH) | 신경내분비(t-NEPC) 전환과 연관된 bHLH 인자; ARPI 이후 염색질 재편(AR 프로그램 약화, 신경계 프로그램 개방)을 유도한다는 보고가 있음(Nouruzi et al. 2022) | E-box 모티프(`CANNTG`) 상에서 고강도 풋프린트 검출 |

---

## 4. 단일 핵 수준의 모티프 활성도 정량화: ChromVAR 알고리즘

단일 핵 ATAC-seq 데이터는 개별 세포 단위에서 풋프린트 계곡을 그리기에는 리드 뎁스(Depth per cell)가 너무 얕습니다(보통 세포당 5,000~30,000 리드).
Greenleaf 연구실이 개발한 **`ChromVAR`**는 단일 핵 수준에서 특정 전사 인자 패밀리의 총체적 조절 활성도를 측정하기 위해, **게놈 전체에 퍼져 있는 수천 개의 동일 TF 모티프 피크를 단일 세포 내에서 집계(Aggregation)하는 통계적 편차(Deviation) 모델**을 사용합니다.

### 4.1 ChromVAR 수학적 정규화 단계
1. **관측 모티프 카운트 ($O_{i, j}$)**:
   - 세포 $j$에서 전사 인자 모티프 $i$를 포함하는 모든 피크 $k \in P_i$에 떨어진 리드 수의 합산:
   $$O_{i, j} = \sum_{k \in P_i} X_{k, j}$$
2. **기대 모티프 카운트 ($E_{i, j}$)**:
   - 피크 $k$의 전체 세포 평균 접근성과 세포 $j$의 총 시퀀싱 뎁스를 바탕으로 산출된 기대값:
   $$E_{i, j} = \sum_{k \in P_i} \left( \frac{\sum_{j'} X_{k, j'}}{\sum_{k', j'} X_{k', j'}} \right) \cdot \left( \sum_{k'} X_{k', j} \right)$$
3. **GC 함량 및 뎁스 매칭 배경 피크 집합 생성 (Background Matching)**:
   - 모티프 $i$를 포함하는 피크들과 **(1) GC 함량, (2) 평균 발현 빈도**가 동일하지만 해당 모티프는 없는 배경 피크 집합을 50세트($b = 1, \dots, 50$) 무작위 추출하여 배경 기대치 분포 $\mu_{\text{bg}}, \sigma_{\text{bg}}$를 계산합니다.
4. **Z-편차 점수 (Z-score of Deviation)**:
   - chromVAR의 원시 편차(raw deviation)는 기대값으로 나눈 상대 편차 $Y_{i,j} = (O_{i,j} - E_{i,j}) / E_{i,j}$이며, 배경 세트에서 같은 방식으로 계산한 $Y^{(b)}_{i,j}$의 평균·표준편차로 편향 보정 편차와 Z-점수를 구합니다(Schep et al. 2017). 이를 통해 시퀀싱 뎁스·GC 편향이 상당 부분 보정된 세포별 TF 모티프 활성도를 얻습니다:

$$\text{bias-corrected deviation}_{i,j} = Y_{i,j} - \mu_{\text{bg}}(i,j), \qquad Z_{i, j} = \frac{Y_{i,j} - \mu_{\text{bg}}(i, j)}{\sigma_{\text{bg}}(i, j)}$$

   - 여기서 $\mu_{\text{bg}}, \sigma_{\text{bg}}$는 배경 세트의 $Y^{(b)}_{i,j}$ 평균과 표준편차입니다.

```python
# ChromVAR 편차 계산 워크플로우 (Python / pychromVAR + Muon 인터페이스)
# 주의: muon(mu.atac)에는 chromVAR 구현이 없다(mu.atac.tl.chromvar, mu.atac.tl.motif_enrichment 없음).
# Python에서는 별도 패키지 pychromvar(pip install pychromvar)가 chromVAR 알고리즘을 구현한다.
# (pychromvar 0.0.4는 AnnData(dtype=)를 사용하므로 anndata 0.13 이상에서는 compute_deviations가 TypeError; anndata<0.13 사용)
# R 대안: Bioconductor chromVAR(computeDeviations) 또는 ArchR(addDeviationsMatrix)
import pandas as pd
import muon as mu
import scanpy as sc
import pychromvar as pc
from pyjaspar import jaspardb

atac = mdata.mod['atac']   # 피크 이름은 "chr1-1000-1500" 형식이어야 함

# 1. 피크 서열 추출 및 피크별 GC 편향 계산 (hg38 FASTA 필요)
pc.add_peak_seq(atac, genome_file="hg38.fa", delimiter="-")   # atac.uns['peak_seq']
pc.add_gc_bias(atac)                                           # atac.var['gc_bias']

# 2. GC 함량·평균 접근성이 비슷한 배경 피크 50세트 선정 (atac.varm['bg_peaks'])
pc.get_bg_peaks(atac, niterations=50)

# 3. JASPAR 2024 CORE 전사 인자 모티프 매칭 (MOODS 기반; atac.varm['motif_match'])
jdb = jaspardb(release="JASPAR2024")
motifs = jdb.fetch_motifs(collection="CORE", tax_group=["vertebrates"])
pc.match_motif(atac, motifs=motifs)

# 4. ChromVAR Z-score 편차 계산
# 반환값은 AnnData(Cells x Motifs)이며 var_names는 "MA0148.4.FOXA1"처럼 "매트릭스ID.TF이름" 형식
dev = pc.compute_deviations(atac)

def motif_col(tf):
    # JASPAR 이름은 대소문자가 섞여 있어(예: AR 모티프 MA0007.4는 "Ar") 대문자로 비교; 동일 TF의 모티프가 여럿이면 첫 번째 사용
    return dev.var_names[dev.var_names.str.split(".").str[-1].str.upper() == tf][0]

# 5. 전립선암 핵심 전사 인자 모티프 활성도 UMAP 시각화
for tf in ['FOXA1', 'AR', 'ASCL1']:
    mdata.obs[f'chromvar_{tf}'] = pd.Series(dev[:, motif_col(tf)].X.ravel(), index=dev.obs_names)

mu.pl.umap(mdata, color=['leiden_wnn', 'chromvar_FOXA1', 'chromvar_AR', 'chromvar_ASCL1'], cmap='coolwarm')
```

---

## 5. 결론 및 다음 단계 연결

풋프린팅의 입체 장애 메커니즘과 서열 편향 보정, 그리고 `ChromVAR`를 통한 단일 핵 단위의 Z-편차 점수화를 통해, 전립선암 세포 내부에서 어떤 전사 인자가 후성유전학적으로 염색질을 개방하고 활성화되었는지 정량화할 수 있게 되었습니다.

이제 개별적으로 규명된 **(1) 전사 인자 발현량, (2) 인핸서 개방성 및 풋프린트, (3) 표적 유전자 발현**의 세 요소를 하나로 묶어, 세포 운명을 결정짓는 궁극의 조절 회로를 재구성할 때입니다. 다음 장 [Chapter 08. 인핸서 주도 전사조절망 복원: SCENIC+ 및 삼중 eRegulon 아키텍처](08_EGRN_AND_SCENICPLUS.md)에서는 원거리 인핸서를 매개로 하는 차세대 유전자 조절망 복원 기술을 배웁니다.

---

## 6. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### Q1. [능동 회상 (Active Recall)] ATAC-seq 풋프린팅(Footprinting) 분석 시 Tn5 트랜스포사제의 고유한 염기서열 선호도(Sequence Insertion Bias)를 수학적으로 보정하지 않을 경우 발생하는 치명적인 위양성(False-positive) 발생 메커니즘을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **위양성 발생 메커니즘**:
  - Tn5 효소는 완전히 무작위로 절단하는 것이 아니라, 특정 염기서열 문맥에 우선적으로 삽입되는 내인성 서열 선호도를 지닙니다(이중뉴클레오타이드 수준의 편향이 보고되었으며, 단일 consensus로 환원되지 않습니다).
  - 만약 특정 전사 인자의 동족 모티프(Cognate motif) 염기서열이 우연히 Tn5가 물리적으로 삽입되기 극도로 꺼려하는 비선호 염기 조합으로 이루어져 있다면, 단백질이 전혀 결합하지 않은 벌거벗은 상태(Unbound/Naked DNA)에서도 해당 위치의 Tn5 절단 빈도가 0에 가깝게 낮아집니다.
  - 서열 편향을 보정하지 않으면 알고리즘은 이 자연적인 절단 결손을 전사 인자의 물리적 입체 장애로 오인하여, 실제로는 비활성 상태인 전사 인자가 강력하게 결합해 있다는 정반대의 가짜 풋프린트(False-positive footprint)를 보고하게 됩니다. 따라서 `TOBIAS`처럼 기대 편향 절단($E_i$)을 관측 신호에서 차감($C_i = O_i - E_i$)하는 디바이아싱을 거쳐야 합니다.

</details>

---

### Q2. [능동 회상 (Active Recall)] `ChromVAR` 알고리즘이 단일 세포 내에서 단순 모티프 카운트가 아닌 '배경 매칭 순열(Background Peak Matching Z-score)'을 사용하여 편차(Deviation)를 산출하는 이유와, 배경 피크를 선정할 때 반드시 통제해야 하는 2가지 물리적/기술적 변수를 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **배경 매칭 순열을 사용하는 이유**:
  - 단일 핵 ATAC 데이터는 세포마다 시퀀싱 뎁스(포획된 총 단편 수) 편차가 극심하며, 게놈 상에 수천 개가 존재하는 모티프와 수십 개만 존재하는 모티프 간의 원시 카운트 불균형이 큽니다. 단순 카운트를 사용하면 단순히 뎁스가 깊은 세포가 모든 TF의 활성이 높은 것처럼 왜곡됩니다.
- **반드시 통제해야 하는 2대 변수**:
  1. **GC 함량 (GC Content)**: PCR 증폭 편향 및 Tn5 침투율은 국소 GC 비율에 민감하게 좌우되므로, 모티프 피크와 동일한 GC% 분포를 갖는 배경 피크를 추출해야 염기 조성에 의한 기술적 왜곡을 상쇄할 수 있습니다.
  2. **전체 세포 평균 접근성 (Average Accessibility / Peak Depth)**: 유전체 상에서 상시 열려 있는 하우스키핑 프로모터 피크와 세포 특이적 원거리 인핸서 피크는 기저 접근성 빈도가 완전히 다릅니다. 기저 접근성 빈도가 일치하는 배경 피크 집합(50세트)을 기준으로 관측치의 표준편차 거리($Z = (O - E - \mu_{\text{bg}}) / \sigma_{\text{bg}}$)를 환산해야만 순수한 세포 유형 특이적 모티프 활성도를 추출할 수 있습니다.

</details>

---

### Q3. [임상 케이스 해석 (Clinical Interpretation)] 호르몬 치료에 저항성을 보이는 mCRPC 환자의 Multiome 검체를 분석하였다. 종양 세포의 전사체 데이터에서는 *AR*(안드로겐 수용체) mRNA 발현량이 대조군 대비 5배 높게 유지되고 있었다. 그러나 `ChromVAR` 분석 결과 AR 모티프의 Z-편차 점수는 음수($Z = -2.4$)로 급락해 있었고, 반대로 파이오니어 인자인 *ASCL1* 및 *NEUROD1* 모티프의 Z-편차 점수는 극단적인 양수($Z > +6.5$)로 치솟아 있었다. 풋프린트 프로파일에서도 *AR* 결합 부위의 계곡은 완전히 소실되어 있었고 *ASCL1* 모티프의 깊은 계곡이 관찰되었다. 이 종양 세포 내부에서 일어나고 있는 전사 조절 기전을 해석하고, 임상의가 환자에게 차세대 AR 저해제(다롤루타마이드)를 투여하려 할 때 이를 중단시켜야 하는 병리학적 근거를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

- **전사 조절 기전 해석 (가설)**:
  - chromVAR 편차는 해당 모티프를 가진 피크 집합의 접근성 증감(Schep 2017: 같은 모티프를 공유하는 피크의 접근성 gain/loss)일 뿐 TF의 단백질 점유율이 아닙니다. 따라서 AR 모티프 $Z<0$는 AR 표적 피크의 접근성 감소를, ASCL1/NEUROD1 모티프 $Z>+6.5$는 bHLH E-box(`CANNTG`) 모티프를 가진 피크의 접근성 증가를 시사할 뿐입니다. E-box는 ASCL1·NEUROD1을 포함한 여러 bHLH TF가 공유하므로 TF 특이성이 낮습니다.
  - *AR* mRNA가 5배 높은 점은 AR 증폭이나 음성 되먹임 같은 다른 설명과도 양립하므로, '신경내분비 계통 전환이 완료되었다'고 결론 내리려면 신경내분비 마커(*SYP*, *CHGA* 등) 발현, 병리 소견, *RB1*/*TP53* 상태 같은 직교 근거가 필요합니다. 이 시나리오는 가상이며, 해석은 '신경내분비 프로그램 개방 가설'로 제한됩니다.
- **치료 방침에 대하여**:
  - 이 데이터만으로 AR 표적 치료 중단이나 백금 기반 화학요법, EZH2/DLL3 표적 치료로의 전환을 권고할 수 없습니다(이들은 확립된 표준 권고가 아니며 위 직교 검증과 임상 가이드라인에 따라 결정해야 함).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Schep AN, Wu B, Buenrostro JD, et al. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nat Methods*. 2017;14(10):975-978. PMID: 28825706. doi: 10.1038/nmeth.4401
2. Bentsen M, Goymann P, Schultheis H, et al. ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation. *Nat Commun*. 2020;11(1):4267. PMID: 32848148. doi: 10.1038/s41467-020-18035-1
3. Nouruzi S, Ganguli D, Tabrizian N, et al. ASCL1 activates neuronal stem cell-like lineage programming through remodeling of the chromatin landscape in prostate cancer. *Nat Commun*. 2022;13(1):2282. PMID: 35477723. doi: 10.1038/s41467-022-29963-5
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- muon 0.2.0 (github.com/scverse/muon). 원 논문: Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. *Genome Biol*. 2022;23(1):42. PMID: 35105358. doi: 10.1186/s13059-021-02577-8
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- pychromVAR (github.com/pinellolab/pychromvar).
- pyJASPAR (github.com/asntech/pyjaspar).
<!-- REFERENCES:END -->

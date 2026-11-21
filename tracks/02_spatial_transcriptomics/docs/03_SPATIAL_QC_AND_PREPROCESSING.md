# Chapter 03. 공간전사체 품질 관리(QC)와 전처리 기법 (Spatial QC & Preprocessing)

> **Cognitive Goal**: 조직 광학 영상과 분자 카운트 간의 불일치를 판별하는 조직 마스킹(Tissue Masking) 원리를 이해하고, FFPE 특화 프로브 기반 데이터에서의 품질 메트릭스(Mitochondrial/Ribosomal 비율)를 단일세포 전사체와 비교 해석하며, 음성 대조군 프로브(Negative Control Probes)를 활용한 위발견율(FDR) 수학적 정량화 및 조직 찢김·접힘·측면 확산 블러링 아티팩트를 식별·보정할 수 있다.

---

## 1. 조직 영역 마스킹 (Tissue Masking & Background Filtering)

공간전사체 실험 슬라이드에는 실제 생체 조직이 놓인 영역과 bare glass(빈 유리 슬라이드) 영역이 공존합니다. 1차 분석 단계는 순수한 조직 내부 스팟/빈(In-tissue)만을 정확히 격리하는 **조직 마스킹(Tissue Masking)**입니다.

```
[조직 마스킹 및 아티팩트 판별 흐름]
Raw Slide Image (H&E / DAPI) ──► Otsu 이진화 / 영상 마스킹 ──► 조직 영역 추출 (In-tissue)
                                                                 │
    ┌────────────────────────────────────────────────────────────┴────────────────────────────────┐
    ▼                                                            ▼                                ▼
[정상 조직 스팟]                                              [조직 접힘 (Tissue Folding)]      [유리 영역 오염 (RNA Leakage)]
- H&E 조직 영상 정상                                           - H&E 상 짙은 보라색 겹침          - H&E 상 조직 없음 (Glass)
- UMI 밀도 및 유전자 수 정상                                   - UMI 및 유전자 수 비정상 급증     - 주변 조직으로부터 확산된 RNA 검출
```

### 1.1 영상 기반 마스킹 vs 발현량 기반 필터링
1. **광학 영상 기반 마스킹 (Image-based Thresholding)**:
   - 명시야(Brightfield) H&E 영상 또는 DAPI 형광 영상의 그레이스케일 강도 히스토그램에 **오츠(Otsu) 이진화 알고리즘** 또는 적응형 임계치(Adaptive Thresholding)를 적용하여 조직 외곽 마스크를 생성합니다.
   - 모폴로지 연산(Morphological Closing: 팽창 후 침식)을 통해 조직 내부의 미세 관강이나 미소 구멍이 잘못 배경으로 분류되는 것을 방지합니다.
2. **발현량 기반 필터링 (UMI-based Filtering)**:
   - 스팟/빈의 총 UMI 카운트($N_{\text{UMI}}$)가 특정 역치(예: $N_{\text{UMI}} \ge 500$) 이상인 지점을 조직으로 간주합니다.
3. **Space Ranger 기본 `in_tissue` 플래그의 오류 가능성**:
   - 상용 소프트웨어(10x Space Ranger)가 자동으로 부여하는 `in_tissue == 1` 플래그는 완벽하지 않습니다.
   - **조직 박리 및 찢김(Tears)**: 영상 분석 상으로는 조직 경계 내부로 판정되었으나 실제로는 미세절편 제작 중 조직이 찢겨 나가 RNA가 전혀 검출되지 않는 위양성(False-positive in-tissue) 스팟이 발생합니다.
   - **조직 외 부유 RNA 오염(Floating Ambient RNA)**: 투과화 효소 세척이 불완전할 경우 조직 외부의 순수 유리 표면에 떠돌던 RNA가 결합하여 `in_tissue == 0` 영역에서도 수백 개의 UMI가 검출되는 현상이 빈번합니다.

---

## 2. 공간적 스팟/빈 품질 지표 (Spot Metrics & FFPE Chemistry)

### 2.1 핵심 품질 관리 메트릭
각 공간 단위(스팟, 빈, 분할된 단일세포)에 대해 계산되는 기본 QC 통계량:
1. **총 UMI 수 ($N_{\text{UMI}}$ / `total_counts`)**: 캡처된 전사체 분자의 총합.
2. **검출된 유전자 수 ($N_{\text{genes}}$ / `n_genes_by_counts`)**: 1 카운트 이상 발현된 고유 유전자 수.
3. **미토콘드리아 유전자 비율 (`pct_counts_mt`)**: 전체 UMI 중 미토콘드리아 유래 mRNA($MT\text{-}CO1, MT\text{-}ND1$ 등)의 백분율.
4. **리보솜 유전자 비율 (`pct_counts_ribo`)**: 세포질 리보솜 단백질($RPL, RPS$ 계열) 유전자 발현 백분율.
5. **헤모글로빈 유전자 비율 (`pct_counts_hb`)**: 적혈구 오염 및 혈관 내강 영역 식별 ($HBA1, HBA2, HBB$).

### 2.2 신선 동결(FF) vs FFPE 조직에서의 미토콘드리아 QC 해석의 근본적 차이
단일세포 RNA 시퀀싱(scRNA-seq)에서 `pct_counts_mt > 15%`는 세포막이 천공되어 세포질 RNA는 유출되고 무거운 미토콘드리아만 남은 **사멸 세포(Dying/Apoptotic cell)**를 배제하는 절대적 QC 기준입니다. 그러나 **FFPE 기반 공간전사체에서는 이 기준을 무비판적으로 적용해서는 안 됩니다**:

```
[FF vs FFPE 미토콘드리아 QC 메커니즘 비교]

1. Fresh Frozen (Poly-A 전장 시퀀싱):
   세포막 파열 ──► 세포질 mRNA 유출 (Loss) ──► 미토콘드리아 mRNA 잔존 ──► MT 비율 인위적 급증 (사멸 신호)

2. FFPE 공간전사체 (Probe Hybridization 패널):
   - 포르말린 가교로 인해 세포막 천공에 의한 세포질 유출이 물리적으로 차단됨.
   - 프로브 패널에 미토콘드리아 유전자 프로브가 있는지는 **패널/프로브 세트 버전에 따라 다름**:
     Visium Human Transcriptome Probe Set v1.0에는 MT 프로브가 0개이나 v2.0(CytAssist FFPE)에는 12개 MT 유전자
     (MT-ATP6, -CO1/2/3, -CYB, -ND1/2/3/4/4L/5/6) 프로브가 포함됨(10x 프로브 세트 CSV 직접 확인). Xenium은 패널별로 확인 필요.
   - 따라서 사용한 프로브 세트 CSV에서 MT 프로브 유무를 먼저 확인하고, v1.0처럼 MT 프로브가 없으면 %mt는 QC 지표로 쓸 수 없으며,
     v2.0처럼 있으면 %mt를 해석하되 FFPE의 세포막 보존 특성상 임계값은 신선 동결과 달리 데이터 분포로 정함.
```

---

## 3. 광학 노이즈와 음성 대조군 기반 위발견율 (FDR) 제어

이미징 기반 공간전사체(Xenium, CosMx, MERSCOPE)는 초고감도 형광 이미징을 수행하므로 다양한 비생물학적 광학 노이즈에 노출됩니다.

### 3.1 광학 노이즈의 3대 근원
1. **내인성 자가형광 (Tissue Autofluorescence)**:
   - 전립선암 기질의 콜라겐, 엘라스틴, 세포 노화 지표인 리포푸신(Lipofuscin) 과립, 적혈구의 헴(Heme) 구조가 가시광선 영역 전반에서 자가형광을 방출하여 형광 피크 검출을 방해.
2. **광학 채널 간 블리딩 (Spectral Bleed-through / Crosstalk)**:
   - 형광체의 방출 스펙트럼 꼬리가 인접 필터 채널로 새어 들어가 가짜 바코드 비트($0 \to 1$)를 유발.
3. **비특이적 혼성화 (Non-specific Hybridization)**:
   - 프로브가 표적 외 유사 서열이나 조직 기질 단백질에 물리적으로 엉겨 붙어 신호를 발생.

### 3.2 음성 대조군 프로브(Negative Control Probes)와 코드워드(Codewords)
시스템 수준의 위양성을 측정하기 위해 두 가지 유형의 음성 대조군이 설계에 포함됩니다:
- **Negative Control Probes (Target Probes)**: 인간 게놈/전사체에 존재하지 않는 인공 서열(예: 해파리 GFP, 세균 서열, 또는 Scrambled sequence)을 표적으로 삼는 실제 올리고뉴클레오타이드 프로브.
- **Negative Control Codewords (Unassigned Codewords)**: 순환 이미징의 유효 코드북(Codebook) 상에 어떤 유전자도 할당하지 않고 남겨둔 빈 바코드 조합. (이 코드워드가 검출된다면 이는 $100\%$ 광학적 디코딩 에러를 의미함).

### 3.3 위발견율 (False Discovery Rate, FDR) 정량 수학 공식
검출된 전체 표적 전사체 신호 중 위양성 분자의 비율인 **프로브 단위 위발견율($\text{FDR}_{\text{probe}}$)**은 다음과 같이 수학적으로 정량화됩니다:

$$\text{FDR}_{\text{probe}} = \frac{\frac{\sum_{i=1}^{M_{\text{Blank}}} \text{Counts}_{\text{Blank}, i}}{M_{\text{Blank}}}}{\frac{\sum_{j=1}^{N_{\text{Target}}} \text{Counts}_{\text{Target}, j}}{N_{\text{Target}}}} = \frac{\overline{\text{Counts}}_{\text{Blank}}}{\overline{\text{Counts}}_{\text{Target}}}$$

- $M_{\text{Blank}}$: 패널에 포함된 음성 대조군(Blank) 프로브의 총 개수 (통상 $30 \sim 50$개).
- $N_{\text{Target}}$: 실제 생물학적 표적 유전자 프로브의 개수 (통상 $300 \sim 5,000$개).
- $\text{Counts}_{\text{Blank}, i}$: $i$번째 음성 대조군 프로브에서 검출된 총 신호 분자 수.
- $\text{Counts}_{\text{Target}, j}$: $j$번째 표적 유전자 프로브에서 검출된 총 신호 분자 수.

통상적으로 신뢰할 수 있는 단일 분자 공간전사체 데이터셋은 **$\text{FDR}_{\text{probe}} < 1\%$ ($0.01$)** 미만으로 유지되어야 하며, 개별 분자에 대해서는 Phred 품질 점수 $Q = -10 \log_{10}(P_{\text{error}}) \ge 20$ (신뢰도 $99\%$ 이상)을 만족하는 전사체만을 하류 분석에 투입합니다.

---

## 4. 물리적 아티팩트: 조직 찢김, 접힘 및 RNA 측면 확산 블러링

### 4.1 조직 찢김(Tearing) 및 박리(Detachment)의 진단
- **병리 영상과 UMI 히트맵의 국소 상관계수 저하**:
  - 정상적인 조직에서는 H&E 영상의 광학 밀도(Optical Density)와 UMI 카운트 밀도 간에 높은 양의 상관관계가 나타납니다.
  - 박리나 찢김이 발생한 영역에서는 H&E 영상 상에 조직 그림자가 남아 있더라도 UMI 카운트가 급격히 0으로 떨어지는 **비상관성 이상치(Discordant Outlier)**가 나타납니다.

### 4.2 조직 접힘(Tissue Folding)의 식별
- 박절(Sectioning) 과정에서 조직 절편이 겹쳐지면 H&E 영상에서 헤마톡실린 염색이 극도로 짙어지며, 스팟 당 UMI 수와 유전자 검출 수가 주변 정상 영역 대비 2~3배 이상 비정상적으로 폭증합니다.
- 이러한 접힘 영역은 단일 세포 레이어가 아닌 다층 세포가 중첩된 상태이므로, 정상 조직 분석을 왜곡하지 않도록 사전에 마스킹하여 제거해야 합니다.

### 4.3 RNA 측면 확산 블러링 (Lateral Diffusion Blur) 정량화
조직 경계면(Tissue boundary)에서 빈 슬라이드(Glass) 영역으로 방출된 RNA가 확산되어 나가는 현상은 1차원 지수 감쇄 모델(Exponential Decay Model)로 정량화할 수 있습니다:

$$C(r) = C_0 \exp\left( - \frac{r}{\lambda_d} \right) + C_{\text{bg}}$$

- $r$: 조직 최외곽 경계선으로부터 수직 방향으로 떨어진 물리적 거리 ($\mu\text{m}$).
- $C(r)$: 거리 $r$에서의 UMI 농도 밀도 ($\text{count}/\mu\text{m}^2$).
- $C_0$: 조직 가장자리 표면에서의 전사체 농도.
- $\lambda_d$: **특성 확산 길이(Characteristic Diffusion Length)**.
- $C_{\text{bg}}$: 슬라이드 전반의 기저 배경 노이즈 농도.

$\lambda_d$가 스팟 직경이나 비닝 크기(예: Visium HD $2\ \mu\text{m}$)보다 현저히 큰 경우(예: $\lambda_d > 20\ \mu\text{m}$), 심각한 투과화 과다(Over-permeabilization)로 인해 공간적 경계면이 블러링되었음을 의미합니다.

---

## 5. Python 실습: `Scanpy` & `Squidpy` 공간 QC 파이프라인

`Scanpy`와 `Squidpy`를 활용하여 공간 데이터의 이상치 스팟을 필터링하고 조직 접힘 및 배경 영역을 격리하는 표준 코드입니다.

```python
import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

def run_spatial_qc_pipeline(adata: sc.AnnData, 
                            min_counts: int = 500, 
                            max_counts_quantile: float = 0.995,
                            min_genes: int = 300) -> sc.AnnData:
    """
    공간전사체 AnnData 객체에 대한 공간 품질 관리(Spatial QC) 실행
    """
    # 1. 기본 QC 지표 계산 (Scanpy)
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPL", "RPS", "rpl", "rps"))
    adata.var["hb"] = adata.var_names.str.startswith(("HBA", "HBB", "HBD"))
    
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ribo", "hb"], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    # 2. 총 UMI 상위 분위수 이상치 후보 계산 (주의: 정의상 항상 상위 (1-q)가 표시되므로 '접힘'의 증거가 아니며,
    #    접힘 판정은 H&E 영상 마스크 또는 이웃 스팟 대비 국소 이상치로 확인해야 함. 밀도 높은 종양 코어도 UMI가 높을 수 있음)
    upper_count_limit = np.quantile(adata.obs["total_counts"], max_counts_quantile)
    
    print(f"[QC] 필터링 전 총 스팟 수: {adata.n_obs:,}")
    print(f"[QC] 99.5% 카운트 상한치 (조직 접힘 감지용): {upper_count_limit:.1f} counts")
    
    # 3. 이상치 및 비조직 스팟 마스킹
    valid_spots = (
        (adata.obs["in_tissue"] == 1) &
        (adata.obs["total_counts"] >= min_counts) &
        (adata.obs["total_counts"] <= upper_count_limit) &
        (adata.obs["n_genes_by_counts"] >= min_genes)
    )
    
    # 이상치 메타데이터 기록
    adata.obs["qc_pass"] = valid_spots
    adata.obs["is_folded"] = adata.obs["total_counts"] > upper_count_limit
    
    # 4. 공간 UMI 밀도 및 QC 지표 시각화
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    sq.pl.spatial_scatter(adata, color="total_counts", ax=axs[0], title="Total UMI Counts", cmap="viridis")
    sq.pl.spatial_scatter(adata, color="n_genes_by_counts", ax=axs[1], title="Detected Genes", cmap="plasma")
    sq.pl.spatial_scatter(adata, color="qc_pass", ax=axs[2], title="QC Pass Spots", cmap="coolwarm")
    plt.tight_layout()
    plt.show()
    
    # 필터링 적용된 정제 객체 반환
    adata_filtered = adata[adata.obs["qc_pass"]].copy()
    print(f"[QC] 필터링 후 잔존 유효 스팟 수: {adata_filtered.n_obs:,}")
    return adata_filtered
```

---

## 6. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
이미징 기반 공간전사체(Xenium, CosMx 등)에서 비특이적 광학 노이즈와 형광 채널 간 블리딩을 통제하기 위해 도입된 **Negative Control Probes**와 **Negative Control Codewords**의 설계 차이점을 규명하고, 검출된 신호로부터 **프로브 단위 위발견율($\text{FDR}_{\text{probe}}$)**을 계산하는 수학 공식을 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **Negative Control의 설계 차이점**:
   - **Negative Control Probes (실제 물리적 프로브)**: 인간 전사체에 존재하지 않는 외래 서열(예: 세균 유전자, 해파리 형광단백질, 인공 난수 서열)을 표적으로 하도록 합성된 실제 올리고 프로브입니다. 이 프로브에서 신호가 검출된다면 이는 화학적 비특이적 결합(Non-specific binding)이나 조직 단백질에 대한 물리적 점착에 기인한 에러를 반영합니다.
   - **Negative Control Codewords (순수 정보적 빈 바코드)**: 순환 이미징의 이진 바코드 체계(Codebook) 상에서 유효한 바코드 패턴 중 어떤 유전자에도 할당하지 않고 공백으로 남겨둔 코드워드입니다. 실제 프로브가 투입되지 않으므로, 이 바코드가 검출된다면 이는 형광 소광 실패, 광학적 중첩, 또는 채널 간 형광 블리딩으로 인한 **순수 광학적 디코딩 오류(Optical Decoding Error)**를 의미합니다.
2. **프로브 단위 위발견율($\text{FDR}_{\text{probe}}$) 수식**:
   $$\text{FDR}_{\text{probe}} = \frac{\frac{\sum_{i=1}^{M_{\text{Blank}}} \text{Counts}_{\text{Blank}, i}}{M_{\text{Blank}}}}{\frac{\sum_{j=1}^{N_{\text{Target}}} \text{Counts}_{\text{Target}, j}}{N_{\text{Target}}}} = \frac{\overline{\text{Counts}}_{\text{Blank}}}{\overline{\text{Counts}}_{\text{Target}}}$$
   - 여기서 분자는 음성 대조군 프로브 1개당 평균 검출 카운트이며, 분모는 생물학적 표적 프로브 1개당 평균 검출 카운트입니다. 이 비율이 통상 $1\%$ ($0.01$) 미만일 때 전체 데이터셋의 분자 신뢰성이 검증된 것으로 판정합니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
신선 동결(Fresh Frozen) 조직을 사용한 단일세포 전사체(scRNA-seq) 분석에서는 높은 미토콘드리아 RNA 비율(`pct_counts_mt > 15~20%`)을 **세포막 손상 및 사멸 세포(Apoptotic Cell)**의 핵심 배제 기준으로 사용합니다. 그러나 포르말린 고정 파라핀 포매(FFPE) 전립선암 조직에 패들록 프로브(Padlock Probe) 혼성화 방식의 공간전사체를 적용할 때, 이 미토콘드리아 QC 기준을 동일하게 적용할 수 없는 생물물리학적 및 분자생물학적 이유 2가지를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **포르말린 단백질-핵산 가교(Formaldehyde Cross-linking)에 의한 누출 차단**:
   - 신선 동결 scRNA-seq에서는 세포가 손상되면 세포막 천공 부위를 통해 세포질 mRNA가 외부로 빠져나가고, 크고 무거운 미토콘드리아만 핵 주변에 남아 상대적으로 미토콘드리아 전사체 비율이 급증합니다.
   - 반면 FFPE 조직은 고정 과정에서 포름알데히드에 의해 세포질과 세포막 단백질이 촘촘한 공유결합 메틸렌 가교망으로 굳어져 있습니다. 세포질의 용해나 선택적 유출이 원천적으로 일어날 수 없으므로 세포막 손상에 따른 미토콘드리아 비율 상승 메커니즘 자체가 성립하지 않습니다.
2. **표적 프로브 패널 설계(Panel Design)의 차이**:
   - FFPE 공간전사체(Visium FFPE, Xenium 등)는 Poly(A) 꼬리를 임의로 잡는 방식이 아니라, 패널에 포함된 특정 프로브 쌍만을 결합시키는 표적 프로브 혼성화(Targeted Probe Hybridization) 방식을 사용합니다.
   - 미토콘드리아 게놈에서 전사되는 13개 단백질 코딩 유전자($MT\text{-}CO1, MT\text{-}ND1$ 등)의 프로브 포함 여부는 프로브 세트 버전에 따라 다릅니다(Visium v1.0: 0개, v2.0: 12개 MT 유전자). MT 프로브가 없는 패널에서는 미토콘드리아 카운트 비율이 세포 생존율 지표가 될 수 없으며, 있는 패널에서는 사용한 CSV로 확인한 뒤 해석해야 합니다.

</details>

---

### [문제 3. 임상 사례 분석 (Clinical Interpretation Case Study)]
**[임상 시나리오]**
연구팀은 전립선암 환자의 수술 검체 절편을 10x Visium(55 µm 스팟) 플랫폼으로 분석하였습니다. (Visium HD의 2 µm 빈은 면적이 55 µm 스팟의 약 1/594이라 스팟당 UMI 규모가 달라 이 시나리오의 수치와 맞지 않습니다.)
QC 분석 과정에서 두 가지 기이한 공간 패턴이 관찰되었습니다:
1. **패턴 A**: H&E 병리 영상에서 진한 보라색 선관 구조를 보이는 Gleason Pattern 4 종양 영역 내부에서 스팟 당 UMI 카운트가 주변 평균($1,200$ UMI) 대비 3배 이상 높은 $4,500$ UMI로 치솟았으며, 인접 스팟과의 경계선이 뚜렷한 선형(Linear) 띠 형태로 나타났습니다.
2. **패턴 B**: H&E 영상 상에서 명백히 조직이 존재하지 않는 순수 유리 슬라이드(Bare Glass) 영역인데도 불구하고, 조직 외곽 경계선으로부터 약 $30\ \mu\text{m}$ 떨어진 지점까지 $KLK3$ 및 $ACPP$와 같은 전립선 특이적 암 유전자 전사체가 지수 함수 형태로 감쇄하며 다수 검출되었습니다.

**[질문]**
1. 패턴 A와 패턴 B가 나타난 원인을 조직 준비 및 공간 생물물리학적 아티팩트 관점에서 각각 진단하시오.
2. 패턴 B의 확산 감쇄 곡선($C(r) = C_0 \exp(-r/\lambda_d)$)에서 특성 확산 길이 $\lambda_d$가 $15\ \mu\text{m}$로 계산되었다면, 이 데이터셋을 사용하여 세포 간 직접 접촉 신호(Juxtacrine interaction, 거리 $< 5\ \mu\text{m}$)를 분석할 때 발생할 수 있는 치명적 오류를 비판하고, QC 파이프라인에서 취해야 할 정제 전략을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **패턴 A와 패턴 B의 아티팩트 진단**:
   - **패턴 A 진단 (조직 접힘, Tissue Folding)**:
     마이크로톰 박절 및 슬라이드 부착 과정에서 얇은 절편의 일부가 접히면서 두 개 이상의 조직 층이 물리적으로 겹쳐진(Overlapped) 현상입니다. H&E 상에서 헤마톡실린 염색이 두 배로 짙어지며, 스팟 당 물리적 세포 수가 2~3배로 증가하여 UMI 카운트가 비정상적으로 폭증한 것입니다.
   - **패턴 B 진단 (RNA 측면 확산 블러링, Lateral Diffusion Blurring)**:
     조직 투과화(Permeabilization) 반응 시간이 과도하게 길어져, 세포막이 붕괴된 후 용출된 고농도의 전립선 암세포 mRNA($KLK3, ACPP$)가 표면 캡처 프로브에 즉시 포획되지 못하고 수용액 상에서 브라운 운동(Brownian diffusion)을 통해 조직 밖 유리 표면으로 번져나간(Diffusion crosstalk) 아티팩트입니다.
2. **오류 비판 및 QC 정제 전략**:
   - **오류 비판**: 특성 확산 길이 $\lambda_d = 15\ \mu\text{m}$는 전사체가 원래 발생한 세포 위치로부터 평균 $15\ \mu\text{m}$ 이상 측면으로 이동하여 캡처되었음을 의미합니다. 이는 세포 간 직접 접촉 신호(Juxtacrine, 유효 거리 $< 2 \sim 5\ \mu\text{m}$)보다 훨씬 넓은 범위입니다. 따라서 종양 세포의 고발현 전사체가 인접한 면역 세포나 기질 세포 스팟으로 누출되어, 실제로는 상피세포에 국한된 리간드가 면역 세포에서도 발현되는 것처럼 보이는 **가짜 세포 간 신호(Spurious Juxtacrine Crosstalk 위양성)**를 대규모로 유발합니다.
   - **QC 정제 전략**:
     - **유리 영역 및 확산 마진 마스킹**: 형태학적 조직 경계선으로부터 확산 거리 $\lambda_d$의 2배($30\ \mu\text{m}$)에 해당하는 외곽 마진 영역의 스팟들을 분석에서 배제(Exclusion buffer zone)합니다. 지수 감쇄 $C(r)=C_0 e^{-r/\lambda_d}$에서 $r=2\lambda_d$에서도 신호가 $e^{-2}\approx13.5\%$ 남고 $3\lambda_d$에서 약 $5\%$이므로 "완전 배제"가 아니며, 마진은 허용 잔여 비율(예: $3\lambda_d$ 이상)로 정당화해야 합니다.
     - **배경 오염 감산 알고리즘 적용(SoupX·DecontX는 비공간 scRNA-seq용 도구이며 공간 적응형 모델은 아님; 공간 데이터에는 SpotClean 등 공간 특화 도구를 사용)**: 유리 영역에서 측정된 전사체 조성비를 바탕으로 조직 내부 스팟에 침투한 확산 전사체(Ambient RNA) 비율을 추정하여 카운트 행렬에서 통계적으로 감산(Subtract) 보정합니다.
     - **패턴 A 제거**: H&E 영상으로 접힘이 확인된 스팟(또는 이웃 스팟 대비 국소 이상치)을 마스킹합니다. 전역 총 UMI 상한(예: $\text{total\_counts} > 3,500$)은 이 가상 데이터에서의 예시일 뿐이며, 세포 밀도가 높은 종양 코어를 선택적으로 제거할 위험이 있어 임계는 데이터 분포를 보고하며 정합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 검수 근거 문헌
1. Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *Gigascience*. 2020;9(12). PMID: 33367645. doi: 10.1093/gigascience/giaa151
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Matplotlib 3.11.2 (github.com/matplotlib/matplotlib). 원 논문: Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput. Sci. Eng*. 2007;9(3):90-95. doi: 10.1109/mcse.2007.55
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- Squidpy 1.8.3 (github.com/scverse/squidpy). 원 논문: Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nat Methods*. 2022;19(2):171-178. PMID: 35102346. doi: 10.1038/s41592-021-01358-2
<!-- REFERENCES:END -->

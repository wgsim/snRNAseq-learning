# Chapter 00. Track 02 로드맵: 전립선암 공간전사체학 (Spatial Transcriptomics) 학습 과정

> **Cognitive Goal**: 단일세포 해리 시퀀싱(sc/snRNA-seq)의 공간적 한계를 극복하고, 전립선암 조직의 고유한 3차원 병리 아키텍처(Gleason 패턴, Cribriform 선암종, CAF 물리적 장벽)를 보존하는 공간전사체학(Spatial Transcriptomics)의 핵심 원리와 기술 분류체계를 이해하며, 재현 가능한 공간 컴퓨팅 환경을 구축하고 전체 학습 로드맵을 확립할 수 있다.

---

## 1. 왜 전립선암 연구에서 공간전사체학(Spatial Transcriptomics)인가?

단일세포/단일핵 RNA 시퀀싱(scRNA-seq / snRNA-seq)은 종양 미세환경(Tumor Microenvironment, TME)을 구성하는 세포들의 분자적 다양성을 규명하는 데 혁명적인 기여를 하였습니다. 그러나 전립선암(Prostate Adenocarcinoma, PRAD) 연구에서 전통적인 단일세포 해리 시퀀싱은 근본적인 생물학적 한계에 부딪힙니다.

```
[단일세포 해리 시퀀싱 (Dissociation-based)의 한계]
생검/수술 조직 ──► 효소 처리 및 기계적 파쇄 ──► 단일세포/핵 현탁액 ──► 시퀀싱
                     (조직 아키텍처 소실)            (해리 편향 발생)
                     - 세포 간 물리적 접촉 소실      - 단단한 기질/상피 세포 사멸
                     - 조직 침윤 경계면 붕괴          - 면역 세포 배제 패턴 관찰 불가

[공간전사체학 (Spatial Transcriptomics)의 패러다임]
온전한 조직 절편 (FFPE / Fresh Frozen) ──► 형태학적 영상 (H&E / IF) + 2차원 공간 발현 동시 획득
  - Gleason 패턴(Grade Group)과 전사체 1:1 대응
  - 암세포-기질세포-면역세포의 공간적 병치(Juxtacrine) 신호 포획
```

### 1.1 조직 해리(Dissociation)에 따른 공간적 정보 소실과 해리 편향
1. **공간적 위상 정보(Spatial Topology)의 완전한 소실**:
   - 해리된 세포 현탁액에서는 암세포가 관강(Lumen) 내부에 국한된 상태인지, 기저막(Basement Membrane)을 뚫고 기질로 침윤 중인 침윤 선단(Invasive Front)에 위치하는지 구분할 수 없습니다.
2. **효소 및 물리적 해리 편향(Dissociation Bias)**:
   - 전립선암 조직은 치밀한 결합조직성 기질(Desmoplastic Stroma)과 세포외 기질(ECM)로 둘러싸여 있습니다. 효소(Collagenase, Dispase) 처리 시 단단하게 얽힌 암세포 및 암 관련 섬유아세포(Cancer-Associated Fibroblasts, CAFs)는 저항성을 보이며 용해되거나 사멸하여, 최종 라이브러리에서 실제 비율보다 과소 대표(Under-represented)됩니다.
3. **근접 세포 간 상호작용(Juxtacrine Signaling) 단절**:
   - 리간드-수용체 상호작용은 세포 간 거리에 결정적으로 의존합니다. 해리된 상태에서는 공간적으로 수백 $\mu\text{m}$ 떨어져 있어 물리적으로 결합할 수 없는 세포 쌍 간에도 통계적 유의성이 나타나는 위양성(False-positive) 세포 간 통신 예측이 빈번하게 발생합니다.

### 1.2 전립선암의 독특한 병리학적 아키텍처와 Gleason Grading
전립선암의 예후 판정 및 치료 방침 결정은 여전히 H&E 염색 기반의 **Gleason Grading System**에 절대적으로 의존합니다:

- **Gleason Pattern 3**: 개별적이고 명확히 분리된 작은 원형 관강(Discrete, well-formed glands)을 형성하며 기질 침윤이 제한적임.
- **Gleason Pattern 4**: 관강 구조가 서로 융합(Fused glands)되거나 불완전하게 형성되며, 특히 중요한 아형으로 **체상 선암종(Cribriform Adenocarcinoma)**이 관찰됨.
  - *Cribriform 구조의 임상적 중요성*: 다공성 격자 형태를 띠며 중심부에 저산소증(Hypoxia)이 있을 것이라는 기전 가설(HIF 표적 유전자 $CA9$, $VEGFA$, $SLC2A1$ 등의 상승으로 평가; $HIF1A$ mRNA 자체는 단백질 안정화로 조절되므로 저산소 지표로 부적합)이 제시되며, 생화학적 재발(Biochemical Recurrence, BCR) 및 원격 전이 위험도가 급격히 상승함.
- **Gleason Pattern 5**: 관강 구조가 완전히 소실되어 시트상 고형 침윤(Sheet-like solid growth), 단일 세포 침윤(Single-cell cords), 또는 중심부 면포 괴사(Comedonecrosis)를 동반함.

공간전사체학은 이러한 미세 병리 구조를 파괴하지 않고, H&E 병리 영상 위에서 Gleason Pattern 3, 4, 5 영역 및 Cribriform 암소포의 전사체 프로파일을 직접 분리하여 정량화할 수 있습니다.

### 1.3 암 관련 섬유아세포(CAF)의 물리적 장벽과 면역 배제(Immune Exclusion)
전립선암은 대표적인 면역학적 "Cold Tumor"(면역 결핍/배제 종양)로 알려져 있습니다. 이는 면역관문억제제(Immune Checkpoint Inhibitors, ICI)에 대한 반응률이 극히 낮은 주요 원인입니다:
- **myCAF(Myofibroblastic CAF) 장벽**: $\alpha\text{-SMA}$($ACTA2$), $FAP$, 콜라겐 섬유($COL1A1, COL3A1$)를 고분비하는 활성화된 섬유아세포가 종양 소포 둘레에 빽빽한 피막 형태의 **결합조직 기질 장벽(Desmoplastic Barrier)**을 형성합니다.
- **CD8+ T 세포의 공간적 고립**: 세포독성 T 세포와 NK 세포는 이 기질 장벽을 뚫지 못하고 종양 외곽 기질(Peritumoral Stroma)에 물리적으로 갇히는 **면역 배제(Immune-excluded)** 패턴을 보입니다.
- 공간전사체학을 통해서만 면역 세포의 기질 내 트래핑 위치와 리간드-수용체 장벽(예: $CXCL12\text{-}CXCR4$, $TGFB1\text{-}TGFBR2$)의 공간 구배를 명확히 입증할 수 있습니다.

---

## 2. 공간전사체학 기술 분류체계 (Taxonomy of Spatial Transcriptomics)

현재 상용화된 공간전사체 기술은 분자 포획 및 판독 메커니즘에 따라 크게 **시퀀싱 기반(Sequencing-based Arrays)**과 **이미징 기반 인시츄(Imaging-based In Situ)** 기술로 양분됩니다.

```
                    ┌───────────────────────────────────────────────┐
                    │     공간전사체 기술 (Spatial Transcriptomics)   │
                    └───────────────────────┬───────────────────────┘
                                            │
            ┌───────────────────────────────┴───────────────────────────────┐
            ▼                                                               ▼
┌───────────────────────────────────────┐       ┌───────────────────────────────────────┐
│     시퀀싱 기반 (Sequencing-based)     │       │    이미징 기반 인시츄 (Imaging-based)  │
│  - 비표적 전장 전사체 (Whole-genome)   │       │  - 표적 유전자 패널 (Targeted Panel)  │
│  - NGS 시퀀싱 연계 (Poly-A / Probe)   │       │  - 형광 현미경 단일 분자 광학 검출   │
└───────────────────┬───────────────────┘       └───────────────────┬───────────────────┘
                    │                                               │
        ┌───────────┴───────────┐                       ┌───────────┴───────────┐
        ▼                       ▼                       ▼                       ▼
 10x Visium Standard      10x Visium HD            10x Xenium In Situ      NanoString CosMx SMI /
 (55 µm spot array)       (2 µm square bin)        (Padlock + RCA)         Vizgen MERSCOPE (MERFISH)
```

### 2.1 시퀀싱 기반 공간 어레이 (Sequencing-based Arrays)
- **10x Visium Standard**:
  - $6.5 \times 6.5\text{ mm}$ 캡처 영역 내 약 5,000개의 스팟 배열.
  - 각 스팟 직경 $55\ \mu\text{m}$, 스팟 중심 간 거리(Pitch) $100\ \mu\text{m}$.
  - 스팟 사이에 $45\ \mu\text{m}$의 비포획 간극(Dead space)이 존재하며, 스팟 하나당 조직 밀도에 따라 1~수십 개의 세포가 혼합(Mixture)되어 단일세포 디콘볼루션(Deconvolution)이 필수적임.
- **10x Visium HD**:
  - 슬라이드 표면에 연속적으로 패턴된 $2\ \mu\text{m} \times 2\ \mu\text{m}$ 정사각형 바코드 영역(프로브 결찰 기반 포획)으로 간극($0\ \mu\text{m}$)이 전혀 없음.
  - 슬라이드당 약 $1.06 \times 10^7$개($6.5\text{ mm}\times6.5\text{ mm}$ ÷ $2\ \mu\text{m}$ 격자 산술)의 $2\ \mu\text{m}$ 빈(Bin)을 제공하여 단일 세포 및 세포하(Subcellular) 분해능을 구현.
  - 통계적 신뢰도 향상을 위해 $2\ \mu\text{m}$, $8\ \mu\text{m}$, $16\ \mu\text{m}$ 스케일로 동적 비닝(Dynamic binning) 수행.

### 2.2 이미징 기반 인시츄 공간전사체 (Imaging-based In Situ)
- **10x Xenium In Situ**:
  - 패들록 프로브(Padlock Probe)와 롤링 서클 증폭(Rolling Circle Amplification, RCA)을 통해 생성된 형광 DNA 볼(Rolony)을 순환 광학 이미징으로 디코딩.
  - 높은 검출 효율과 신호 대 잡음비(SNR), 전립선 FFPE 조직에 최적화된 패널(300~5,000개 유전자).
- **NanoString CosMx SMI**:
  - 광절단성(Photocleavable) 형광 리포터 바코드를 사용하여 1,000~6,000개 전사체 및 단백질을 동시 검출.
- **Vizgen MERSCOPE (MERFISH)**:
  - 다중화 오류 정정 형광 인시츄 혼성화(Multiplexed Error-Robust FISH). 16-bit 이진 바코드 및 해밍 거리 4(Hamming Distance 4, HD4)를 적용하여 광학 노이즈와 형광 손실을 실시간 교정.

### 2.3 기술별 핵심 사양 비교 매트릭스

| 기술 플랫폼 | 분해능/국소화 정밀도 (Resolution) | 포획 유전자 수 (Plexity) | 검출 민감도 (Sensitivity; 정량 순위 아님) | 시료 적합성 (Sample Type) | 주요 전처리 & 분석 도구 |
|:---|:---|:---|:---|:---|:---|
| **10x Visium Standard** | $55\ \mu\text{m}$ (다세포 스팟) | 전장 전사체 (~20,000) | 중간 (Spot당 수천 UMI) | Fresh Frozen & FFPE | `Space Ranger`, `Scanpy`, `cell2location` |
| **10x Visium HD** | $2\ \mu\text{m}$ (단일세포/핵 수준) | 전장 전사체 (~20,000) | Bin당 희소 (누적 합산 필요) | FFPE 특화 | `Space Ranger 3.0+`, `SpatialData` |
| **10x Xenium In Situ** | 광학 분해능 ≈0.2 µm 한계, 분자 국소화 정밀도는 수십~수백 nm (광자수 의존) | 표적 패널 (300 ~ 5,000) | 높음 (패널·조직에 따라 상이) | FFPE & Fresh Frozen | `Xenium Ranger`, `Baysor`, `Cellpose` |
| **NanoString CosMx** | 분자 국소화 정밀도 수십~수백 nm (광학 회절 한계 ≈0.2 µm보다 좋을 수 없음) | 표적 패널 (1,000 ~ 6,000) | 높음 (패널·조직에 따라 상이) | FFPE 특화 | `CosMx AtoMx`, `Giotto`, `InSituType` |
| **Vizgen MERSCOPE** | 분자 국소화 정밀도 수십~수백 nm (회절 한계 ≈0.2 µm보다 좋을 수 없음) | 표적 패널 (500 ~ 1,000) | 높음 (보정 방식·패널에 따라 변동; 플랫폼 간 순위는 벤치마크 문헌 필요) | Fresh Frozen & FFPE | VPT(Vizgen Post-processing Tool), `Baysor`, `Squidpy` |

---

## 3. Track 02 전체 학습 단계별 인지 구조 (Cognitive Levels)

```
[Level 1: Biophysical Foundation & Data Standards]
  ├── Ch 01. 공간전사체 기술 원리와 광학·생물물리학적 기초 (RCA, MERFISH, 아베 한계, RNA 확산)
  └── Ch 02. 조직병리학 영상 정합과 공간 좌표계 체계 (H&E, Gleason, Affine 변환, SpatialData OME-Zarr)
          │
[Level 2: Signal Purification & Single-Cell Segmentation]
  ├── Ch 03. 공간전사체 품질 관리(QC)와 전처리 기법 (FFPE 프로브 특성, 광학 노이즈, FDR 계산)
  └── Ch 04. 고해상도 세포 분할 알고리즘과 전사체 할당 (Cellpose 그래디언트 플로우 vs Baysor MRF 포인트 프로세스)
          │
[Level 3: Cell Identity Mapping & Deconvolution]
  ├── Ch 05. 스팟 수준 세포 디콘볼루션 (cell2location 계층적 베이지안 모델링 & RCTD)
  └── Ch 06. 공간 도메인 클러스터링 및 조직 미세구조 분할 (BayesSpace 공간 사전확률 & SpaGCN 그래프 합성곱)
          │
[Level 4: Spatial Topology & Intercellular Crosstalk]
  ├── Ch 07. 공간적 자기상관 및 가변 유전자 탐색 (Moran's I, Geary's C, SpatialDE 가우스 과정)
  └── Ch 08. 세포-세포 상호작용 및 공간적 틈새(Niche) 분석 (Squidpy 이웃 그래프 & COMMOT 유적 최적 수송)
          │
[Level 5: Clinical Translation & 3D Histopathology]
  └── Ch 09. 다중 슬라이드 3D 정합, 클론 진화 매핑 및 전립선암 치료 저항성 공간적 임상 중개
```

---

## 4. 선수 지식 점검 (Prerequisites)

- **조직병리학 및 전립선암 임상 지식**:
  - H&E 염색 및 면역형광(IF) 기본 원리.
  - 전립선 해부학: 이행대(Transition Zone), 중심대(Central Zone), 말초대(Peripheral Zone, 암 발생의 70% 이상).
  - Gleason Score ($3+3=6$에서 $5+5=10$) 및 세계비뇨기병리학회(ISUP) Grade Group 1~5 분류.
- **공간 수학 및 기하학적 기초**:
  - 2차원 아핀 변환(Affine Transformation: 회전, 전단, 이동, 스케일링).
  - 들로네 삼각분할(Delaunay Triangulation) 및 보로노이 다이어그램(Voronoi Tessellation).
  - 공간 통계학: 공간 가중치 행렬(Spatial Weight Matrix $W$), 공간 자기상관 지수.
- **Python 공간 컴퓨팅 프레임워크**:
  - `SpatialData`: OME-Zarr 기반 멀티모달 공간 데이터 표준 컨테이너.
  - `Squidpy`: 공간 그래프 연산 및 이웃 풍부도 분석 도구.
  - `Scanpy` & `AnnData`: 고차원 전사체 행렬 조작.

---

## 5. 재현 가능한 공간전사체 분석 컴퓨팅 환경 구축 (Environment Setup)

공간전사체 분석은 기가바이트(GB) 단위의 대용량 고해상도 WSI(Whole Slide Image) 병리 영상, 수천만 개의 단일 분자 포인트 데이터, 그리고 딥러닝 기반 세포 분할 모델을 다룹니다. 따라서 최적화된 독립 가상환경이 필수적입니다.

### 5.1 Conda 가상환경 생성 및 의존성 패키지 설치
```bash
# 1. 공간전사체 전용 Python 3.10 가상환경 생성
conda create -n spatial_env python=3.10 -y
conda activate spatial_env

# 2. 필수 지리정보 및 영상 처리 C-라이브러리 설치
conda install -c conda-forge gdal geos proj libspatialindex -y

# 3. 단일세포 및 공간전사체 표준 라이브러리 설치
pip install scanpy anndata squidpy spatialdata spatialdata-io spatialdata-plot

# 4. 이미지 처리 및 딥러닝 세포 분할 프레임워크 설치
pip install tifffile zarr dask[complete] geopandas shapely
pip install cellpose torch torchvision --index-url https://download.pytorch.org/whl/cpu

# 5. 공간 통계 및 디콘볼루션 필수 도구
# pyro-ppl 1.9.2는 cell2location 0.1.5가 import하는 deep_getattr를 제거했으므로 <1.9.2로 고정한다
pip install "pyro-ppl<1.9.2" cell2location scvi-tools

# 6. Jupyter 커널 연결
python -m ipykernel install --user --name spatial_env --display-name "Python (Spatial Transcriptomics)"
```

### 5.2 환경 무결성 검증 스크립트
```python
import scanpy as sc
import squidpy as sq
import spatialdata as sd
import cellpose
import shapely
import geopandas as gpd

print(f"Scanpy Version: {sc.__version__}")
print(f"Squidpy Version: {sq.__version__}")
print(f"SpatialData Version: {sd.__version__}")
print(f"Cellpose Version: {cellpose.__version__}")
print(f"GeoPandas Version: {gpd.__version__}")
print("모든 공간전사체 핵심 라이브러리가 성공적으로 로드되었습니다.")
```

---

## 6. 실전 벤치마크 공개 공간전사체 데이터셋 카탈로그 (Public Spatial Datasets)

전립선암 공간전사체 분석 파이프라인 개발 및 검증을 위해 본 과정에서 활용하는 주요 공개 벤치마크 데이터셋입니다:

| 연구 (논문) | 접근 식별자 / 출처 | 기술 플랫폼 | 시료 특성 (Sample Cohort) | 핵심 생물학적 / 분석 포인트 |
|:---|:---|:---|:---|:---|
| **Berglund et al. (*Nat Commun*, 2018)** | accession 미확정 (원 논문 Data Availability 확인 필요: EGA 표기 `EGAS0000100300`은 EGA 조회로 확인 불가; 카운트 행렬은 spatialtranscriptomicsresearch.org) | Legacy ST ($100\ \mu\text{m}$) | 전립선 전절제술 조직 (n=12 섹션) | 역사상 최초의 전립선암 공간전사체; Gleason 3 vs 4 vs 기질 간 전사 구배 |
| **10x Genomics Official (2024)** | `10x Dataset Portal` | 10x Visium HD ($2\ \mu\text{m}$) | 인간 전립선암 FFPE 조직 (Adenocarcinoma) | $2\ \mu\text{m}$ 해상도에서 Cribriform 침윤 선단 및 미세 관강 구조 규명 |
| **10x Genomics Official (2023)** | `10x Dataset Portal` | 10x Xenium In Situ | 인간 전립선암 FFPE 조직 (377-plex) | 수백만 개의 단일 분자 전사체 기반 세포 분할 및 침윤 면역 세포 니치 분석 |
| **Erickson et al. (*Nature*, 2022)** | `EGAS00001006124` (EGA) | 공간 전사체 (플랫폼 세부는 원 논문 확인 필요) | 전립선을 포함한 여러 장기의 양성·악성 공존 조직 (120,000개 이상 영역) | 공간적 유전체 복제수 변이(CNA) 추정 및 양성-악성 조직 간 클론 관계 규명 |

### 6.1 10x Genomics Visium HD 및 Xenium 벤치마크 데이터 디렉토리 구조
실습의 재현성을 위해 프로젝트 루트에 다음과 같은 정형화된 데이터 저장 디렉토리를 구축합니다:

```bash
mkdir -p data/spatial/{visium_hd_prostate,xenium_prostate,erickson_st}
```

```
data/spatial/
├── visium_hd_prostate/
│   ├── square_002um/               # 2 µm 고해상도 특성 매트릭스
│   ├── square_008um/               # 8 µm 비닝 매트릭스 (세포 크기 수준)
│   ├── square_016um/               # 16 µm 비닝 매트릭스 (빠른 탐색용)
│   ├── spatial/
│   │   ├── tissue_positions.parquet# 비드별 공간 물리 좌표
│   │   └── tissue_lowres_image.png # 저해상도 H&E 이미지
│   └── binned_outputs/
└── xenium_prostate/
    ├── transcripts.parquet         # 단일 전사체 분자 (x, y, z, gene, q_score)
    ├── cells.parquet               # 분할된 세포 메타데이터 (면적, 중심점)
    ├── cell_boundaries.parquet     # 세포 다각형 경계선 벡터 지오메트리
    └── morphology_focus.ome.tif    # DAPI 고해상도 OME-TIFF 형광 영상
```

---

## 7. 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
단일핵 전사체 시퀀싱(snRNA-seq) 단독 연구와 비교하여, 전립선암 조직에서 공간전사체학(Spatial Transcriptomics) 기술이 제공하는 3가지 결정적 병리학적/분자적 가치를 기술하고, Gleason Pattern 4의 특이 형태인 **체상 선암종(Cribriform Adenocarcinoma)** 분석에서 조직 해리 시퀀싱이 왜 심각한 위음성 및 정보 왜곡을 초래하는지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **공간전사체학의 3대 결정적 가치**:
   - **조직학적 위상 및 병리 구조의 보존**: 병리학적 진단의 골드 스탠다드인 Gleason Pattern (3, 4, 5)과 분자 발현 프로파일을 1:1로 직접 오버레이하여, 형태학적 이상과 종양 분자 표현형 간의 인과관계를 입증할 수 있습니다.
   - **해리 편향(Dissociation Bias) 없는 전 체계적 포획**: 물리적/효소적 해리에 취약하여 사멸하기 쉬운 치밀한 결합조직 기질(Desmoplastic CAFs)과 취약한 악성 상피세포를 현장(In situ)에서 온전히 보존하여 분석합니다.
   - **진정한 국소 세포 간 상호작용(Juxtacrine/Paracrine) 규명**: 세포 간 물리적 거리(수 $\mu\text{m}$ 이내)에 기반한 리간드-수용체 상호작용 및 공간적 면역 배제 장벽을 정확히 식별합니다.
2. **Cribriform 선암종 분석에서 조직 해리의 치명적 왜곡**:
   - Cribriform 구조는 세포들이 치밀하게 융합된 다공성 격자 구조를 이룹니다. 이 구조의 중심부는 혈관 접근성이 극히 떨어져 심각한 **만성 저산소증(Hypoxia) 및 영양 결핍**에 노출되어 있습니다.
   - 조직 해리를 수행하면 중심부의 저산소 취약 세포들이 효소 처리 및 원심분리 과정에서 우선적으로 파괴되어 사멸(Lysis)합니다.
   - 설령 생존하더라도, 주변의 정상 산소 상태를 유지하는 외곽 기질 세포들과 무작위로 섞여버리므로, Cribriform 중심부 특이적인 저산소증 반응(HIF 표적) 유전자($VEGFA, CA9$ 등)와 전이 유도 신호가 평균 발현값에 희석되어 위음성(False-negative)으로 소실됩니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
시퀀싱 기반 공간전사체(Sequencing-based Arrays, 예: Visium HD)와 이미징 기반 인시츄 공간전사체(Imaging-based In Situ, 예: Xenium In Situ)의 작동 원리를 비교하고, ① 유전자 포획 범위(Plexity), ② 물리적 분해능(Resolution), ③ 세포 분할(Cell Segmentation)의 의존성 측면에서 두 기술의 근본적인 차이점과 분석 파이프라인 상의 고려사항을 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **원리 비교**:
   - **Visium HD (Sequencing-based)**: 고체 표면에 $2\ \mu\text{m}$ 단위로 고정된 공간 바코드 올리고뉴클레오타이드 어레이 위에서 조직의 RNA를 캡처한 후 차세대 시퀀싱(NGS)을 통해 염기서열을 해독합니다.
   - **Xenium (Imaging-based In Situ)**: 조직 내부에서 표적 mRNA에 패들록 프로브를 결합시키고 롤링 서클 증폭(RCA)을 통해 형광 신호를 증폭한 뒤, 현미경 하에서 형광 리포터의 순환 혼성화 이미징을 통해 단일 분자의 좌표($x, y, z$)와 유전자 정체성을 직접 해독합니다.
2. **3대 핵심 측면 비교**:
   - **유전자 포획 범위(Plexity)**:
     - Visium HD는 전장 전사체(Whole Transcriptome, 약 20,000개 유전자)를 편향 없이 탐색(Unbiased discovery)할 수 있어 신규 바이오마커 발굴에 유리합니다.
     - Xenium은 수백~수천 개(300~5,000개)로 사전에 정의된 표적 유전자 패널에 국한되지만, 검출 효율이 극도로 높습니다.
   - **물리적 분해능(Resolution)**:
     - Visium HD는 $2\ \mu\text{m}$ 정사각형 격자(Grid) 단위로 데이터를 수집하므로, 세포 경계와 비드 격자가 일치하지 않는 격자화 오차(Discretization artifact)가 존재합니다.
     - Xenium은 광학 회절 한계(≈0.2 µm) 안팎의 분해능에서 개별 RNA 분자마다 독립적인 실수형 좌표를 부여하며, 분자 국소화 정밀도는 광자수에 따라 수십~수백 nm입니다.
   - **세포 분할(Cell Segmentation)의 의존성**:
     - Visium HD는 고정된 크기의 빈($2\ \mu\text{m}, 8\ \mu\text{m}$)을 합산하여 분석할 수 있으므로 분할 없이도 무감독 분석이 가능합니다.
     - Xenium은 단일 분자 포인트 구름(Point cloud) 형태로 데이터가 생성되므로, 개별 분자를 특정 세포의 발현으로 귀속시키기 위한 DAPI/세포막 기반의 정밀한 세포 분할(Cell Segmentation, 예: Cellpose, Baysor) 알고리즘 적용이 분석의 성패를 가르는 절대적 선행 조건입니다.

</details>

---

### [문제 3. 임상 사례 분석 (Clinical Interpretation Case Study)]
**[임상 시나리오]**
Gleason Score $3+4=7$ (Grade Group 2)로 진단받은 64세 전립선암 환자의 근치적 전립선 절제술 조직을 대상으로 공간전사체 분석을 수행하였습니다.
병리 판독 결과, 종양 중심부에는 Gleason 4 Cribriform 구조가 광범위하게 관찰되었고, 그 주변은 빽빽한 섬유성 결합조직 기질이 둘러싸고 있었습니다.
공간전사체 발현 분석 결과:
1. 종양 외곽 기질 영역에서는 $\alpha\text{-SMA}$($ACTA2$), $FAP$, $COL1A1$, $TGFB1$의 발현이 매우 높게 나타났습니다.
2. 세포독성 T 세포 마커인 $CD3D, CD8A, GZMB, PRF1$을 발현하는 세포들은 기질 영역 외곽에만 고밀도로 군집해 있었으며, 종양 내부 Cribriform 암세포 소포 내부로는 전혀 진입하지 못하는 양상을 보였습니다.
3. Cribriform 암세포 군집 내부에서는 $AR, KLK3$와 함께 저산소증 반응 유전자(HIF 표적인 $CA9$, $LDHA$ 등) 및 $VEGFA$가 강하게 검출된 가상 사례입니다.

**[질문]**
1. 이 환자의 종양 미세환경을 면역학적 표현형(Immunological Phenotype) 관점에서 진단하고, T 세포가 종양 소포 내부로 침투하지 못하는 공간적/생물학적 기전을 설명하시오.
2. 이 환자에게 항-PD-1/PD-L1 면역관문억제제 단독 요법을 투여할 경우 예상되는 임상 반응성을 평가하고, 공간전사체 분석 결과를 바탕으로 T 세포의 종양 침투를 유도하기 위한 논리적인 복합 치료 병용 전략(Combination Therapy)을 2가지 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **면역학적 진단 및 기전 분석**:
   - **진단**: 전형적인 **면역 배제(Immune-Excluded / Immune Cold) 종양 표현형**입니다.
   - **기전**:
     - $\alpha\text{-SMA}^+ / FAP^+$ 활성화 myCAF가 생성한 고밀도 콜라겐($COL1A1$) 복합체가 종양 주위에 물리적 피막(Physical Desmoplastic Barrier)을 형성하여 세포독성 T 세포의 물리적 혈관외 유출 및 이동을 차단합니다.
     - myCAF에서 분비된 고농도의 $TGF\text{-}\beta 1$은 T 세포의 이동성을 억제하고 주화성을 교란합니다.
     - Cribriform 암소포 내부의 저산소증(HIF 표적 $CA9, LDHA$ 상승으로 추정되는) 환경은 젖산 축적에 따른 세포외 산성화를 유발하여 침투를 시도하는 소수의 T 세포마저 기능적으로 무력화(Anergy)시킵니다.
2. **임상 반응 평가 및 복합 치료 전략**:
   - **단독 요법 반응성**: 불응(Non-responsive). 활성화된 CD8+ T 세포가 표적 암세포와 직접 물리적으로 접촉하지 못하므로, PD-1 억제제로 T 세포의 탈진을 회복시켜도 종양 사멸 효과를 발휘할 수 없습니다.
   - **복합 치료 전략**:
     - **TGF-β 억제제(또는 FAP 표적 기질 리모델링제) + 면역관문억제제**: $TGF\text{-}\beta$ 수용체 저해제(예: Galunisertib)를 병용하여 myCAF의 콜라겐 합성 및 세포외기질 형성을 억제함으로써 기질 물리적 장벽을 해체하고 T 세포의 종양 내부 침투를 유도합니다.
     - **항-VEGF 신생혈관 억제제(예: Bevacizumab) + 면역관문억제제**: Cribriform 내부에서 분비되는 $VEGFA$ 신호를 차단하여 비정상적인 종양 미세혈관을 정상화(Vascular Normalization)하고, 내피세포의 ICAM-1/VCAM-1 발현을 증가시켜 T 세포의 종양 중심부 혈관 투과 및 호밍을 촉진합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Berglund E, Maaskola J, Schultz N, et al. Spatial maps of prostate cancer transcriptomes reveal an unexplored landscape of heterogeneity. *Nat Commun*. 2018;9(1):2419. PMID: 29925878. doi: 10.1038/s41467-018-04724-5
2. Erickson A, He M, Berglund E, et al. Spatially resolved clonal copy number alterations in benign and malignant tissue. *Nature*. 2022;608(7922):360-367. PMID: 35948708. doi: 10.1038/s41586-022-05023-2
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Cellpose 4.2.1.1 (github.com/MouseLand/cellpose). 원 논문: Stringer C, Wang T, Michaelos M, et al. Cellpose: a generalist algorithm for cellular segmentation. *Nat Methods*. 2021;18(1):100-106. PMID: 33318659. doi: 10.1038/s41592-020-01018-x
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- Squidpy 1.8.3 (github.com/scverse/squidpy). 원 논문: Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nat Methods*. 2022;19(2):171-178. PMID: 35102346. doi: 10.1038/s41592-021-01358-2
- GeoPandas 1.2.0 (github.com/geopandas/geopandas).
- Shapely 2.2.0 (github.com/shapely/shapely).
- SpatialData 0.8.0 (github.com/scverse/spatialdata). 원 논문: Marconato L, Palla G, Yamauchi KA, et al. SpatialData: an open and universal data framework for spatial omics. *Nat Methods*. 2025;22(1):58-62. PMID: 38509327. doi: 10.1038/s41592-024-02212-x
<!-- REFERENCES:END -->

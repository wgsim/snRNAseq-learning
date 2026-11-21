# 00. 암 단일세포 전사체 기반 유전체 및 클론 진화 분석 로드맵 (Cancer Single-Cell Genomics & Clonal Evolution Roadmap)

---

## 학습 목표 (Cognitive Goals)

1. **단일세포 전사체 군집화의 내재적 한계 규명**: 순수 유전자 발현량(Gene Expression Count) 기반 클러스터링(t-SNE, UMAP, Leiden/Louvain) 시 발생하는 '반응성 기질/정상 상피세포(Reactive Normal Epithelia)의 악성 세포 오분류' 메커니즘을 규명하고, 왜 유전체적 복제수 변이(Copy Number Alteration, CNA) 프로파일링이 악성 세포 판별의 절대적 표준(Ground Truth)인지 설명할 수 있다.
2. **이수성(Aneuploidy)과 Nowell 클론 진화 모델의 수리생물학적 이해**: 고형암의 결정적 특징인 대규모 이수성(Aneuploidy)의 형성 원리를 이해하고, 피터 노웰(Peter Nowell)의 클론 진화 모델(선형 진화, 분지형 진화, 중립 진화)을 단일세포 유전체 분기 트리(Phylogenetic Tree) 및 클론 파편화(Clonal Architecture) 관점에서 수학적으로 모델링할 수 있다.
3. **글로벌 암 전사체 및 유전체 벤치마크 코호트 분석**: 전이성 전립선암 클론 파종을 전장 유전체로 입증한 Gundem et al. (*Nature* 2015), 원발암 다중오믹스 표준인 TCGA-PRAD, 그리고 거세저항성 전립선암(mCRPC)의 핵심 코호트인 SU2C/PCF 데이터셋의 분자생물학적 특징을 파악한다.
4. **대규모 단일세포 전사체 복제수 분석 환경 구축**: Python 생태계의 `scanpy`, `infercnvpy`, R 연동 `numbat`, 계통수 최적화 도구 `fastme`를 연계한 고성능 연산 파이프라인을 설계하고 재현 가능한 분석 환경을 구성할 수 있다.

---

## 1. 순수 전사체 군집화(Naive Transcriptomic Clustering)의 치명적 한계

단일세포 RNA 시퀀싱(scRNA-seq 및 snRNA-seq)은 종양 미세환경(Tumor Microenvironment, TME)의 이질성을 단일세포 해상도로 조명하는 혁신적 도구로 자리잡았다. 그러나 **전사체 프로파일(Transcriptomic Profile)만으로 악성 종양 세포(Malignant Tumor Cells)와 비악성 세포(Non-malignant Cells)를 분리하는 관행**은 심각한 방법론적 오류를 수반한다.

### 1.1 전사체 가소성과 염증성 유도체에 의한 표현형 왜곡

세포의 유전자 발현 프로파일은 유전체(Genotype)뿐만 아니라 국소 미세환경 자극(Microenvironmental Stimuli)에 의해 유연하게 변화하는 표현형(Phenotype)이다. 특히 암 조직 내 침윤된 대식세포(TAMs), T세포, 섬유아세포(CAFs)가 분비하는 염증성 사이토카인($\text{TNF-}\alpha$, $\text{TGF-}\beta$, $\text{IL-6}$ 등) 및 저산소 미세환경(Hypoxia)은 정상 기저세포(Basal cells)나 관강세포(Luminal cells)에 극심한 반응성 전사 변화(Reactive Atypia)를 유도한다:

*   **상피-간엽 이행(EMT-like State)**: 기질과 인접한 정상 관강 상피세포가 $\text{TGF-}\beta$ 신호전달계 활성화에 의해 관강 특이적 분화 마커(예: *KLK3*, *PRAC1*)의 발현을 낮추고 *VIM*, *FN1*, *MMP* 계열 유전자를 상향 조절한다.
*   **스트레스 반응 및 증식 신호 모방**: 세포 밀도가 높거나 괴사(Necrosis) 경계면에 위치한 정상 세포들이 열충격 단백질(Heat Shock Proteins) 및 전사인자 *JUN*, *FOS*, *MYC* 표적 유전자를 일시적으로 과발현하여 공격적인 암세포 군집과 유사한 유전자 세트를 공유하게 된다.

### 1.2 반응성 정상 상피(Reactive Normal Epithelia) 오분류의 수리적 메커니즘

단일세포 분석의 표준 차원 축소 알고리즘(PCA $\to$ UMAP) 및 그래프 기반 군집화(Leiden/Louvain)는 고차원 유전자 발현 벡터 간의 유클리드 거리(Euclidean Distance) 또는 코사인 유사도(Cosine Similarity)를 기반으로 작동한다:

$$D(\mathbf{x}_i, \mathbf{x}_j) = \sqrt{\sum_{g=1}^{G} (x_{ig} - x_{jg})^2}$$

여기서 상위 고변이 유전자(Highly Variable Genes, HVGs, 통상 2,000~3,000개) 중 상당수가 세포 주기(Cell Cycle, e.g., *MKI67*, *TOP2A*), 대사 스트레스(Metabolic Stress), 리보솜 생합성 단백질 등 종양 특이적 돌연변이와 무관한 생리적 반응 유전자에 의해 지배된다. 그 결과:

1.  **가짜 종양 아형 생성 (False Subclonal Splitting)**: 악성 돌연변이가 전혀 없는 완전한 정상 이배체(Diploid) 세포 집단이 단지 국소 염증에 반응했다는 이유로 UMAP 상에서 독립된 클러스터로 분리되어 '고위험 신규 암 서브클론'으로 잘못 명명된다.
2.  **종양 클러스터 오염 (Tumor Cluster Contamination)**: 반대로, 증식이 느리고 분화도가 높은 저위험 암세포가 정상 관강세포와 동일한 주성분 공간(PC Space)에 투영되어 정상 세포 클러스터 내부로 흡수된다.

따라서 암 생물학에서 악성 세포를 판정하기 위한 절대적 기준은 일시적 전사체 반응이 아니라, 비가역적인 **유전체 구조적 결함(Genomic Structural Defects)**, 즉 **체세포 복제수 변이(Somatic Copy Number Alterations, SCNAs)**와 **이수성(Aneuploidy)**이어야 한다.

---

## 2. 악성의 결정적 표지: 이수성(Aneuploidy)과 Nowell 클론 진화 이론

### 2.1 테오도르 보베리 가설부터 현대 유전체학까지

1914년 독일의 세포생물학자 테오도르 보베리(Theodor Boveri)는 다극성 유사분열(Multipolar Mitosis)에 의해 발생하는 염색체 수의 비정상적 불균형이 악성 종양의 본질적 기원이라는 가설을 제시하였다. 현대 암 유전체학(Hanahan & Weinberg, *Hallmarks of Cancer*)에 이르러, 인간 고형암(Solid Tumors)의 **90% 이상**에서 광범위한 염색체 이수성(Aneuploidy)과 염색체 불안정성(Chromosomal Instability, CIN)이 관찰됨이 확립되었다.

정상 체세포는 엄격한 세포분열 검문점(Spindle Assembly Checkpoint, SAC)을 통해 46개의 염색체($2n=46$)를 온전하게 보존한다. 반면 종양 세포는 유사분열 조절 유전자 결함, 중심체 과다 복제(Centrosome Amplification), 텔로미어 위기(Telomere Crisis)를 겪으며 염색체 팔(Arm) 단위의 획득(Gain)과 소실(Loss), 혹은 전장 유전체 배가(Whole-Genome Doubling, WGD)를 겪는다. 이러한 대규모 복제수 변이는 수백~수천 개의 유전자 발현 용량 효과(Gene Dosage Effect)를 유발하여 발암 경로를 영구적으로 재편한다.

### 2.2 피터 노웰(Peter Nowell, 1976)의 클론 진화 모델

1976년 피터 노웰은 단일 전구세포에서 출발한 종양이 선택압 하에서 순차적으로 클론 진화한다는 모델을 제안하였다. 아래 3번의 중립/빅뱅 모델은 노웰의 모델이 아니라 후속 연구(예: Sottoriva 등 2015, 대장암)에서 제안된 모형이며 염색체 분쇄증·WGD와 직접 결부된 모델로 단정할 수 없다.

```
[클론 진화 토폴로지 — 1·2는 Nowell(1976)의 선택 기반 클론 진화, 3은 후속 모형(예: Sottoriva 2015 Big Bang)]

1. 선형 진화 (Linear Evolution)
   [Clone 0] --(Driver A)--> [Clone 1] --(Driver B)--> [Clone 2] (Clonal Sweep)

2. 분지형 진화 (Branched Evolution)
                              +--(Driver B)--> [Clone 1A] (Subclone A)
   [Clone 0] --(Driver A)-----+
                              +--(Driver C)--> [Clone 1B] (Subclone B)

3. 중립 진화 / 빅뱅 모델 (Neutral / Big Bang Evolution)
                                  +--> [Clone 1] (Passenger mutations)
   [Truncal WGD / Catastrophe] ---+--> [Clone 2] (Equal fitness)
                                  +--> [Clone 3] (Pervasive heterogeneity)
```

1.  **선형 진화 (Linear Evolution)**: 새로운 드라이버 돌연변이를 획득한 서브클론이 기존 클론보다 압도적인 적응도(Fitness advantage)를 가져 집단 전체를 대체하는 클론 스윕(Clonal Sweep)이 연속적으로 발생한다.
2.  **분지형 진화 (Branched Evolution)**: 공간적 격리(Spatial Isolation)나 미세환경의 차이로 인해 여러 서브클론이 동시에 독자적인 돌연변이를 축적하며 공존한다. 현재 대부분의 전이성 고형암에서 확인되는 지배적 모델이다.
3.  **중립적 진화 (Neutral Evolution / Big Bang Model; Nowell 모델이 아닌 후속 모형)**: 초기 발생 단계에서 대규모 유전체 격변(Chromothripsis, WGD)이 발생한 후, 후속 돌연변이들은 세포의 적응도에 유의미한 차이를 주지 않고 중립적으로 누적되어 극도의 종양 내 이질성(Intra-tumor Heterogeneity, ITH)을 형성한다.

### 2.3 전립선암(Prostate Cancer)의 특이적 진화 동역학: CNA-Driven Malignancy

대장암(CRC)이나 흑색종(Melanoma)과 달리, 전립선암은 체세포 점돌연변이 빈도(Tumor Mutational Burden, TMB)가 메가베이스(Mb)당 1개 미만으로 매우 낮다. 대신, **전립선암의 악성화, 침윤, 안드로겐 수용체 차단 치료(ADT) 저항성에서 대규모 염색체 구조 변이 및 복제수 변이(CNA)가 핵심적 역할을 한다(단, TCGA 분자 아형의 약 74%가 *ERG* 등 유전자 융합과 *SPOP*, *FOXA1*, *IDH1* 점돌연변이로 정의되므로 융합·SNV도 필수 드라이버이다)**. 따라서 전립선암 단일세포 전사체 연구에서 CNA를 복원하는 것은 클론의 계통수를 재구성하고 환자의 불량 예후를 예측하는 데 핵심적이다.

---

## 3. 표준 벤치마크 데이터셋 카탈로그

단일세포 복제수 추정 알고리즘을 검증하고 암 진화 계통수를 검증하기 위해서는 대규모 전장 유전체(WGS) 및 다중오믹스로 검증된 레퍼런스 데이터셋의 이해가 필수적이다.

| 코호트 / 연구 | 출판 연도 및 저널 | 표본 특성 및 분석 기법 | 유전체적 핵심 의의 |
| :--- | :--- | :--- | :--- |
| **Gundem et al.** | *Nature* (2015) 520:353–357 | 전이성 전립선암 환자 10명의 다부위(Multi-region) 51개 종양 WGS (평균 약 55X; 원발 종양은 5명에서만 포함) | **전이소 간 전파(Metastasis-to-metastasis spread)**가 흔하며(단클론 de novo 파종 또는 5례에서 다클론 전달), AR 신호 유전자 변이가 서로 다른 전이소에서 수렴적으로 반복됨을 보고 |
| **TCGA-PRAD** | *Cell* (2015) 163:1011–1025 | 원발성 전립선 선암종(Primary PRAD) 333례의 WES, RNA-seq, SNP 6.0 Array 다중오믹스 (TCGA-PRAD 프로젝트 전체는 약 500례) | 전립선암의 7대 분자 아형 정의(*ERG*, *ETV1*, *ETV4*, *FLI1* 융합 및 *SPOP*, *FOXA1*, *IDH1* 돌연변이); 8p 결손, 8q 획득 등 표준 CNA 빈도 확립 |
| **SU2C/PCF mCRPC** | *Cell* (2015) 161:1215–1228<br>*PNAS* (2019) 116:11428–11436 | Cell 2015: mCRPC 150례, PNAS 2019: 444례 생검 코호트 (WES + RNA-seq) | 치료 저항성 핵심 축 규명: *AR* 유전자 증폭(GISTIC 기준 SU2C 2015에서 약 52%), *TP53*·*PTEN*·*RB1* 불활성화(빈도는 코호트·정의에 따라 상이). *AR* 650kb 상류 인핸서 tandem duplication은 별도의 Viswanathan 2018(linked-read WGS)에서 mCRPC의 70~87%, 원발암 <2%로 보고 |

```
[Genomic Architecture Comparison: Primary vs. mCRPC]

Primary PRAD (TCGA-PRAD)
  Genome Status: Predominantly Diploid (Low CIN)
  Early Events:  TMPRSS2-ERG fusion (21q deletion), NKX3-1 loss (8p21)
  Driver Alterations: SPOP mutations (10%), PTEN hemizygous loss (20%)

                │
                │ Progression under Androgen Deprivation Therapy (ADT)
                ▼

Metastatic CRPC (SU2C/PCF, Gundem et al.)
  Genome Status: High CIN, Frequent Whole-Genome Doubling (WGD, ~50%)
  Late Events:   AR locus & upstream enhancer amplification (Xq12)
                 Biallelic loss of TP53 (17p), RB1 (13q), PTEN (10q)
                 Cross-organ clonal seeding (Bone, Liver, Lymph nodes)
```

---

## 4. 계산 생물학 분석 환경 및 파이프라인 아키텍처

단일세포 전사체 복제수 분석 파이프라인은 크게 전사체 정량화, 정상 레퍼런스 세포 식별, 염색체 위치 정렬 평활화, 세그멘테이션(HMM), 그리고 계통수 구축 단계로 구성된다.

```
[Single-Cell Cancer Genomics Pipeline Flowchart]

Raw FASTQ (snRNA-seq / scRNA-seq)
         │
         ▼ (Cell Ranger / STARsolo)
Count Matrix (Cells x Genes) + Cell Metadata
         │
         ▼ (Scanpy Quality Control & Normalization)
Filtered AnnData (Mitochondrial filter, Doublet removal)
         │
         ├───────────────────────────────────┐
         ▼                                   ▼
Transcriptomic Cell Typing             Genomic Coordinate Mapping
(Immune, Stroma, Epithelial)           (Ensembl GTF v44/GRCh38)
         │                                   │
         ▼ (Select Reference Normal)         ▼
Reference Baseline Cells               Sorted Genomic Gene Matrix
(Fibroblasts, T cells, Endothelia)                   │
         │                                           │
         └─────────────────┬─────────────────────────┘
                           ▼
             Baseline Normalization & Moving Average
             (infercnvpy / CopyKAT / Numbat)
                           │
                           ▼
             HMM Copy Number State Segmentation
                           │
                           ▼
             Subclone Assignment & Phylogeny (FastME)
```

### 4.1 핵심 의존성 패키지 및 시스템 환경

| 패키지/도구명 | 권장 버전 | 역할 및 수리적 기반 |
| :--- | :--- | :--- |
| **`python`** | `>= 3.10` | 기본 분석 런타임 환경 |
| **`scanpy`** | `>= 1.9.6` | 대규모 단일세포 전사체 전처리, 차원축소, 발현 행렬 구조체(`AnnData`) 관리 |
| **`anndata`** | `>= 0.10.0` | 메모리 효율적 희소 행렬(HDF5/zarr) 및 다중 레이어 유전체 메타데이터 저장 |
| **`infercnvpy`** | `>= 0.4.0` | Python 기반 고속 병렬 inferCNV 구현체 (이동평균 및 이수성 점수 계산) |
| **`scikit-learn`**| `>= 1.3.0` | GMM 이배체 기저선 분리 및 주성분 공간 클러스터링 |
| **`fastme`** | `>= 2.1.6` | 최소 진화 원리(Minimum Evolution Principle) 기반 대규모 계통수 고속 계산 바이너리 |
| **`numbat`** | `>= 1.4.0` | R 기반 대립유전자 특이적(BAF) 단일세포 복제수/하플로타입 HMM 파이프라인 |

### 4.2 재현 가능한 분석 환경 구축 쉘 스크립트

```bash
#!/usr/bin/env bash
# cancer_genomics_env_setup.sh: 암 단일세포 유전체 분석 환경 구성
set -euo pipefail

# 1. Conda 가상환경 생성
conda create -n sc_cancer_genomics python=3.10 -y
# 비대화형 스크립트에서는 쉘 초기화 없이 conda activate가 실패하므로 먼저 초기화 스크립트를 source
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sc_cancer_genomics

# 2. 필수 고성능 C/C++ 의존성 및 계통수 엔진 FastME 설치
conda install -c bioconda -c conda-forge fastme samtools bcftools htslib -y
# (이후 장에서 필요한 R 도구: numbat(R), copykat(R), cellsnp-lite, Eagle2는 별도 설치 필요 — 예: conda install -c bioconda cellsnp-lite eagle2, R에서 devtools::install_github('kharchenkolab/numbat'), install_github('navinlabcode/copykat'))

# 3. Python 핵심 라이브러리 및 단일세포 유전체 도구 설치
pip install --upgrade pip
pip install \
    scanpy==1.9.8 \
    anndata==0.10.5 \
    infercnvpy==0.4.2 \
    scikit-learn==1.4.0 \
    scipy==1.12.0 \
    pandas==2.2.0 \
    matplotlib==3.8.2 \
    seaborn==0.13.2

# 4. 환경 검증
python -c "
import scanpy as sc
import infercnvpy as cnp
import sklearn
print(f'Scanpy version: {sc.__version__}')
print(f'Infercnvpy version: {cnp.__version__}')
print(f'Scikit-learn version: {sklearn.__version__}')
"
echo "단일세포 암 유전체 분석 환경 구성 완료."
```

### 4.3 Scanpy 및 Infercnvpy 초기화 파이프라인 예제

```python
import scanpy as sc
import infercnvpy as cnp
import pandas as pd
import numpy as np

def initialize_cancer_genomics_pipeline(h5ad_path: str, gtf_gene_pos_path: str):
    """
    단일세포 전사체 데이터를 로드하고 유전체 좌표 기반 정렬 및 CNV 계산을 준비한다.
    
    Parameters:
        h5ad_path: 전처리된 Scanpy AnnData 경로
        gtf_gene_pos_path: 유전자별 염색체, 시작, 종료 위치가 기록된 탭 구분 파일
    
    Returns:
        adata: 유전체 위치 정보가 주석화된 AnnData 객체
    """
    # 1. 단일세포 데이터 로드
    adata = sc.read_h5ad(h5ad_path)
    adata.var_names_make_unique()
    
    # 2. 염색체 위치 메타데이터 로드
    # Format: gene_name \t chromosome \t start \t end
    gene_pos = pd.read_csv(
        gtf_gene_pos_path, 
        sep='\t', 
        names=['gene_name', 'chromosome', 'start', 'end'],   # 4열 파일이므로 이름도 4개
        index_col=0
    )
    
    # 3. 데이터셋 내 존재하는 유전자와 유전체 좌표 교집합 추출
    common_genes = adata.var_names.intersection(gene_pos.index)
    adata = adata[:, common_genes].copy()
    gene_pos_filtered = gene_pos.loc[common_genes].copy()   # SettingWithCopyWarning 방지
    
    # 4. 염색체 및 위치 정렬 (chr1 -> chr22 -> chrX -> chrY)
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    gene_pos_filtered['chromosome'] = pd.Categorical(
        gene_pos_filtered['chromosome'], 
        categories=chr_order, 
        ordered=True
    )
    gene_pos_sorted = gene_pos_filtered.sort_values(by=['chromosome', 'start'])
    
    # 정렬 순서대로 AnnData 열 재배치
    adata = adata[:, gene_pos_sorted.index].copy()
    adata.var['chromosome'] = gene_pos_sorted['chromosome']
    adata.var['start'] = gene_pos_sorted['start']
    adata.var['end'] = gene_pos_sorted['end']
    
    print(f"총 {adata.n_obs}개 세포 및 {adata.n_vars}개 유전자에 대한 유전체 좌표 매핑 완료.")
    return adata
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 반응성 비정형 정상 세포의 오분류 판별
임상 연구팀이 전립선암 환자의 전립선 전절제술 검체에서 채취한 8,000개의 단일세포 RNA 시퀀싱 데이터를 분석하였다. UMAP 상에서 고위험 관강 상피세포(Luminal epithelial cells)와 밀접하게 묶이며 증식 마커(*MKI67*, *TOP2A*) 및 EMT 마커(*VIM*)를 높게 발현하는 군집 'Cluster 4'가 발견되었다. 연구팀은 이 군집을 '신규 공격성 암 줄기세포 서브클론'으로 잠정 결론지었다. 이 결론의 타당성을 평가하고, 이를 검증하거나 기각하기 위해 수행해야 할 유전체적 분석 파이프라인을 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 결론의 타당성 평가 (기각)**:
잠정 결론은 타당하지 않다. 종양 미세환경 내에 존재하는 정상 관강 상피세포는 종양 주변 염증(IL-6, TNF-$\alpha$, TGF-$\beta$)이나 국소 허혈성 저산소증에 반응하여 일시적으로 세포 주기를 재가동하고 간엽성 마커(*VIM*)를 상향 조절하는 **반응성 비정형화(Reactive Atypia)** 현상을 겪는다. 순수 전사체 발현량만으로는 이 세포가 진정한 악성 클론인지, 단순한 반응성 정상 세포인지 구별할 수 없다.

**2. 검증 및 기각을 위한 유전체 분석 파이프라인**:
*   **염색체 복제수 변이(CNA) 재구성**: `infercnvpy` 또는 `CopyKAT`을 적용하여 정상 기저 기질세포(섬유아세포, 혈관내피세포, 침윤 면역세포)를 기준(Reference)으로 설정하고, Cluster 4 세포들의 염색체 단위 발현 배율(Expression fold-change)을 전장 염색체에 대해 매핑한다.
*   **전립선암 필수 유전체 랜드마크 확인**:
    *   Cluster 4가 악성 종양 클론이라면 전립선암의 대표적 드라이버 CNA인 **8p 결손(*NKX3-1*), 8q 획득(*MYC*), 10q 결손(*PTEN*)** 등의 대규모 염색체 구조 이상을 보유하는 경우가 많다(단, 원발 TCGA에서 *MYC* gain/amp는 약 31%, *PTEN* 손실은 약 32%에 그치므로 모든 악성 세포가 이를 갖는 것은 아니다).
    *   만약 Cluster 4가 1번부터 22번 염색체까지 이배체(Diploid, CNV magnitude $\approx 0$) 프로파일을 유지한다면, CNA 신호가 없는 것만으로 반응성 정상 상피로 단정하지 말고 판정을 보류하며(CNA가 약하거나 없는 종양과 CN-LOH는 inferCNV로 검출되지 않음), 변이·BAF 등 직교 근거로 확인한 뒤 제외 여부를 결정한다.

</details>

---

### 문제 2: Nowell 클론 진화 모델과 전립선암 치료 저항성 기전
피터 노웰의 클론 진화 모델에서 '선형 진화(Linear Evolution)'와 '분지형 진화(Branched Evolution)'의 분자계통학적 차이를 수리적/트리 구조 관점에서 비교 설명하시오. 또한, 안드로겐 수용체 차단 요법(Androgen Deprivation Therapy, ADT)을 받은 전립선암 환자에서 거세저항성(mCRPC)이 발생하는 과정은 두 모델 중 어느 쪽에 더 부합하는지 SU2C/PCF 코호트의 유전체적 증거를 들어 논하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 수리적/트리 구조 관점의 비교**:
*   **선형 진화**: 체세포 변이 세트 $S_t \subset S_{t+1}$의 포함 관계가 성립하며 계통수 상에 분기점(Branching point) 없이 단일 궤적(Trunk)만 존재한다. 새로운 변이를 획득한 클론의 적응도 $W(C_{new}) \gg W(C_{old})$이므로 이전 클론을 완전히 도태시키는 **선택적 스윕(Selective Sweep)**이 연속적으로 발생한다.
*   **분지형 진화**: 공통 조상 줄기(Trunk)에서 다수의 서브클론이 갈라져 나오며, 각 서브클론 $A$와 $B$는 공유 변이 $S_{\text{trunk}}$ 외에 서로 배타적인 독자 변이 세트($S_A \setminus S_B \neq \emptyset$, $S_B \setminus S_A \neq \emptyset$)를 축적한다. 공간적 격리나 상이한 국소 선택압에 의해 다수의 서브클론이 공존한다.

**2. 전립선암 거세저항성(mCRPC) 발생 양상과 유전체적 증거**:
*   전립선암의 진행은 전형적인 **분지형 진화(Branched Evolution)** 모델을 따른다.
*   **다부위 연구(예: Gundem 2015)의 증거**: (SU2C/PCF의 Robinson 2015는 환자당 단일 생검 설계라 환자 내 분지 진화를 직접 보일 수 없다.) ADT 치료 압력이 가해지면 안드로겐 비의존적 생존 능력을 갖춘 다수의 잠복 서브클론이 서로 다른 분자 경로를 통해 독립적으로 진화한다. 예를 들어, 동일 환자의 전이 병변 내에서 어떤 서브클론은 *AR* 유전자 좌위 및 650kb 상류 인핸서의 고도 증폭(*AR* amplification)을 획득하는 반면, 다른 전이 부위의 서브클론은 *TP53* 및 *RB1*의 이중 결손(Biallelic loss)을 통해 신경내분비 전립선암(t-NEPC, Lineage Plasticity)으로 분기한다. 이는 선형적 대체가 아닌 치료 압력에 의한 다방향성 분지 진화의 결정적 증거이다.

</details>

---

### 문제 3: 대규모 다부위 전이 전장유전체 벤치마크 설계
Gundem et al. (*Nature* 2015) 연구는 단일 환자 내 다부위 전이 전립선암 검체의 WGS를 통해 전이 파종(Metastatic Seeding)의 새로운 패러다임을 규명하였다. 이 연구가 밝혀낸 '상호 전이 파종(Cross-metastatic Seeding)'의 개념을 설명하고, 단일세포 전사체 복제수 분석 파이프라인에서 다부위 생검 데이터를 활용하여 서브클론 간 파종 방향성을 추론할 때 발생할 수 있는 주요 계산생물학적 한계점 2가지를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 상호 전이 파종 (Cross-metastatic Seeding)**:
기존의 전통적 전이 모델은 원발암(Primary tumor)에서 단방향으로 전이 세포가 림프절이나 뼈로 이동한다는 가설이었다. 그러나 Gundem et al.은 유전체 계통수 분석을 통해, 이미 형성된 전이 병변 A(예: 골 전이)의 서브클론이 또 다른 전이 병변 B(예: 간 전이)로 2차 전이를 일으키거나, 심지어 전이 병변에서 다시 원발 부위로 재파종(Self-seeding)되는 복잡한 전이 간 네트워크가 존재함을 증명하였다.

**2. 단일세포 전사체 CNA 파이프라인의 계산생물학적 한계점**:
1.  **전사체 기반 CNA의 불연속적 해상도 및 역수렴(Homoplasy/Convergence) 문제**: 단일세포 RNA 시퀀싱은 발현되는 유전자에만 의존하므로 유전자 밀도가 낮은 유전체 사막(Gene deserts)의 복제수 변화를 추적할 수 없다. 또한 서로 다른 서브클론이 독립적으로 동일한 염색체 팔 결손(예: 전립선암에서 극히 흔한 8p 결손)을 획득할 경우, 분기 계통 분석에서 이를 단일 조상 사건으로 오인하여 파종의 선후 관계를 잘못 연결할 위험이 높다.
2.  **중립 진화 및 전장 유전체 배가(WGD) 후 대립유전자 비특이성(Allelic Ambiguity)**: 전이성 암에서 흔한 WGD($2n \to 4n$)가 발생하면 단순 발현량 비율(InferCNV식 접근)만으로는 실제 유전체 카피 수가 2카피인지 4카피인지 구별할 수 없으며, 모계/부계 대립유전자의 소실 여부(Copy-Neutral LOH)를 파악하지 못해 계통수 상에서 조상 클론과 후손 클론의 분기 시간을 정확히 교정하기 어렵다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Nowell PC. The clonal evolution of tumor cell populations. *Science*. 1976;194(4260):23-8. PMID: 959840. doi: 10.1126/science.959840
2. Sottoriva A, Kang H, Ma Z, et al. A Big Bang model of human colorectal tumor growth. *Nat Genet*. 2015;47(3):209-16. PMID: 25665006. doi: 10.1038/ng.3214
3. Gundem G, Van Loo P, Kremeyer B, et al. The evolutionary history of lethal metastatic prostate cancer. *Nature*. 2015;520(7547):353-357. PMID: 25830880. doi: 10.1038/nature14347
4. Robinson D, Van Allen EM, Wu YM, et al. Integrative clinical genomics of advanced prostate cancer. *Cell*. 2015;161(5):1215-1228. PMID: 26000489. doi: 10.1016/j.cell.2015.05.001
5. Abida W, Cyrta J, Heller G, et al. Genomic correlates of clinical outcome in advanced prostate cancer. *Proc Natl Acad Sci U S A*. 2019;116(23):11428-11436. PMID: 31061129. doi: 10.1073/pnas.1902651116
6. Cancer Genome Atlas Research Network. The Molecular Taxonomy of Primary Prostate Cancer. *Cell*. 2015;163(4):1011-25. PMID: 26544944. doi: 10.1016/j.cell.2015.10.025
7. Viswanathan SR, Ha G, Hoff AM, et al. Structural Alterations Driving Castration-Resistant Prostate Cancer Revealed by Linked-Read Genome Sequencing. *Cell*. 2018;174(2):433-447.e19. PMID: 29909985. doi: 10.1016/j.cell.2018.05.036
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- infercnvpy 0.6.1 (github.com/icbi-lab/infercnvpy).
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
<!-- REFERENCES:END -->

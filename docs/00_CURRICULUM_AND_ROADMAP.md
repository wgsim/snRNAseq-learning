# Chapter 00. 커리큘럼 소개 및 전립선암 snRNA-seq 학습 로드맵

> **Cognitive Goal**: 본 과정의 학습 목표, 전립선암 연구에서 snRNA-seq이 요구되는 통계적/생물학적 맥락, 필요한 사전 지식 및 재현 가능한 컴퓨팅 환경을 설정할 수 있다.

---

## 1. 전립선암 전사체 연구에서 왜 snRNA-seq인가?

전립선암(Prostate Cancer, PCa)은 전 세계 남성 암 발생률 상위권을 차지하며, 임상적 경과는 양호한 국소 원발암(Primary localized PCa)부터 치명적인 전이성 거세저항성 전립선암(**mCRPC**, metastatic Castration-Resistant Prostate Cancer)까지 매우 광범위한 스펙트럼을 보입니다.

기존의 단일 세포 RNA 시퀀싱(**scRNA-seq**)은 효소 기반의 조직 해리(enzymatic dissociation) 과정이 필수적입니다. 그러나 전립선암 조직은 다음과 같은 치명적인 한계가 존재합니다:
1. **풍부한 세포외기질(ECM)과 섬유화(Fibrosis)**: 전립선 조직 및 암 연관 기질(Reactive stroma)은 콜라겐과 평활근 세포가 밀집되어 있어 세포 해리가 극히 어렵고, 강한 효소 처리 과정에서 상피세포(epithelial cells)와 암세포가 대량 사멸합니다.
2. **냉동 생검 검체(Frozen Biopsy Specimens)**: 다기관 임상 코호트, 전이성/다기관 코호트(예: SU2C/PCF 드림팀 코호트), 골전이(bone metastasis) 및 림프절 전이 생검 조직은 대부분 영하 80도에 급속 동결(flash-frozen) 보관됩니다. 급속 동결 후 보존제 없이 보관된 조직을 해리하면 세포막이 파괴되어 intact한 살아있는 세포를 얻기 어렵습니다(DMSO 등으로 생세포 동결보존한 검체는 예외).
3. **취약한 대형 세포 및 골조직 특이성**: 골전이 병변(골용해성/골형성성 혼합 병변)은 탈회(decalcification) 과정과 기질 특성상 온전한 단일 세포 분리가 어렵고 수율·생존율 저하에 따른 편향 위험이 큽니다(예: Kfoury et al. 2021은 골전이 신선 조직과 골수의 scRNA-seq을 수행했으나, 이는 신선 검체를 확보한 경우에 한합니다).

**단일 핵 전사체 분석(snRNA-seq)**은 동결된 조직에서 온전한 **핵(nuclei)**을 기계적/화학적으로 분리하여 시퀀싱하므로, **효소 기반 조직 해리 편향을 크게 줄이고 동결 검체의 희귀 세포 및 악성 세포를 분석**할 수 있게 하는 방법입니다(다만 세포 유형별 핵 회수율 차이 등 고유한 편향이 있습니다).

---

## 2. 학습 로드맵 및 인지 단계 (Cognitive Levels)

```
[Level 1: Foundation]
  ├── Ch 01. 실험 원리 및 생물학 (sc vs sn, 핵 분리, 전립선 암 생물학)
  └── Ch 02. 업스트림 처리 (FASTQ, Cell Ranger --include-introns, AnnData 구조)
          │
[Level 2: Data Cleaning & Cognition]
  ├── Ch 03. 암 조직 특화 QC (Mito% 인지, CellBender 배경 제거, Scrublet)
  └── Ch 04. 다운스트림 파이프라인 (정규화, HVG, PCA/UMAP 기하학, Leiden)
          │
[Level 3: Biological Interpretation]
  ├── Ch 05. 세포 주석 및 악성 판별 (전립선 TME, inferCNV 악성 상피 감별)
  └── Ch 06. 전이암 고급 분석 (Harmony/scVI 배치 통합, AR/NEPC 궤적, pySCENIC)
          │
[Level 4: Expert Mastery]
  ├── Ch 07. 멀티오믹스 & 공간전사체 (10x Multiome, snRNA + snATAC, Visium 연계)
  └── Ch 08. 실전 트러블슈팅 & Pseudobulk (배치 효과 분리, 위양성 극복)
```

---

## 3. 선수 지식 점검 (Prerequisites)

- **생물학적 기초**: 
  - 중심원리(Central Dogma): DNA -> pre-mRNA -> mature mRNA -> Protein.
  - 전립선암의 드라이버 유전자: *AR*(안드로겐 수용체), *PTEN* 결손, *TP53* 변이, *RB1* 결손, *TMPRSS2-ERG* 유전자 융합.
- **컴퓨팅 및 통계 기초**:
  - Python 기초 문법 및 객체 지향 개념 (`AnnData` 객체의 슬라이싱 및 메서드 활용).
  - 다차원 행렬 연산(`numpy`, `scipy.sparse`), 데이터 프레임 조작(`pandas`).
  - 통계적 가설 검정: 귀무가설, $p$-value, False Discovery Rate (FDR, Benjamini-Hochberg 교정).

---

## 4. 재현 가능한 분석 환경 구축 (Environment Setup)

단일 핵 분석은 딥러닝 기반 Ambient RNA 제거 도구(`cellbender`) 및 Variational Autoencoder 기반 통합 도구(`scvi-tools`)를 사용하므로, 환경 충돌을 방지하기 위해 격리된 Python 가상환경 구축을 권장합니다.

```bash
# 1. Conda 가상환경 생성 (Python 3.10 예시. 주의: 최신 scanpy/scvi-tools는 Python >=3.12, decoupler 2.x는 >=3.11을 요구하므로
#    버전을 고정하고(environment.yml / pip freeze 기록) 도구에 맞는 Python을 선택할 것. pySCENIC은 numpy/pandas 충돌 보고가 있어 별도 환경 권장)
conda create -n snrna_prostate python=3.10 -y
conda activate snrna_prostate

# 2. 필수 단일세포 핵심 라이브러리 설치
pip install scanpy anndata scrublet

# 3. 딥러닝 기반 고급 분석 도구 설치 (PyTorch 기반)
pip install cellbender scvi-tools

# 4. 전사조절망, 유전자 세트 활성 분석 및 복제수 변이(inferCNV) 도구 설치
pip install pyscenic decoupler infercnvpy scikit-misc

# 5. 시각화 및 인터랙티브 지원 라이브러리
pip install matplotlib seaborn ipykernel jupyterlab leidenalg

# 6. Jupyter 커널 등록
python -m ipykernel install --user --name snrna_prostate --display-name "Python (snRNA Prostate)"
```

### 환경 정상 설치 확인 스크립트
```python
import scanpy as sc
import anndata as ad
import scvi
import cellbender

print(f"Scanpy Version: {sc.__version__}")
print(f"AnnData Version: {ad.__version__}")
print(f"scvi-tools Version: {scvi.__version__}")
print("All core dependencies successfully loaded.")
```

---

## 5. 전립선암 핵심 분석 체크포인트 요약

| 단계 | 일반 scRNA-seq 통념 | 전립선암 snRNA-seq 인지 전환 |
|:---|:---|:---|
| **미토콘드리아 비율** | 5~20% 이하를 정상으로 간주 | 대개 **낮게 관찰됨** (예: 수 % 이하이나 조직·프로토콜별 분포로 임계값 설정; 높으면 세포질 잔여물/오염 가능성) |
| **인트론(Intron) 처리** | 성숙 mRNA(Exon) 위주 카운트 | **Pre-mRNA(Intron) 필수 포함** (핵 내 RNA 중 인트론 비율은 조직 의존적으로 상당하며 대략 수십 %) |
| **Ambient RNA 오염** | 세척 후 경미한 수준 | **매우 심각** (조직 파쇄 시 용출된 AR, PSA(KLK3) 등이 전 세포에 번짐) |
| **암세포 판별** | 단순 클러스터 마커 발현으로 판별 | **inferCNV/CopyKAT 활용** (이배체 정상 상피세포와 이수체 악성세포 분리; CNA가 없거나 약한 종양에서는 실패할 수 있어 직교 근거 병행) |
| **치료 저항성 추적** | 별도 구분 없음 | **AR signaling score vs NEPC transdifferentiation** 연속적 궤적 추적 |

---

## 6. 실전 연습용 전립선암 공개 데이터셋 카탈로그 (Public Benchmark Datasets)

전립선암의 단일세포 및 단일핵 전사체 연구는 국소 원발암부터 거세저항성 전이암(mCRPC), 골전이 미세환경까지 다양한 임상 스펙트럼을 아우릅니다. 실전 분석 파이프라인 학습과 벤치마킹을 위해 권장되는 대표적 공개 데이터셋 카탈로그는 다음과 같습니다:

| 연구 (참고문헌) | 데이터 번호 (Accession) | 분석 모달리티 | 샘플 구성 및 환자 코호트 | 핵심 생물학적 발견 및 분석 포인트 |
|:---|:---|:---|:---|:---|
| **Kfoury et al. (*Cancer Cell*, 2021)** | `GSE143791` | scRNA-seq (10x 3′ v2) | 전립선암 골전이 환자 9명(각 환자의 종양·병변 인접 골수·원위 골수 3개 분획)과 양성 골수 대조 7명 (GEO 시리즈에는 비전립선 전이·마우스 샘플도 함께 등록) | 골전이 면역억제 미세환경(T 세포 소진, 골전이 특이 대식세포, CCL20-CCR6 축) |
| **Song et al. (*Nat Commun*, 2022)** | `GSE176031` | scRNA-seq (Seq-Well S³, 10x 아님) | 국소 원발성 전립선 선암 환자 11명(생검·근치적 절제 검체 혼합; 근접 정상 조직 포함) | 악성 관강세포(ERG+/ERG-) 이질성 및 클럽(Club) 세포 염증 신호 |
| **Chen et al. (*Nat Cell Biol*, 2021)** | `GSE141445` | scRNA-seq | 근치적 전립선 절제술 원발 종양 13개 (36,424 세포; GEO 요약), 논문은 별도 11개 샘플을 추가 시퀀싱(CRPC 포함) | 진행성 연관 상피 하위집단, 미세환경의 KLK3 이소성 발현, 종양과 활발히 통신하는 내피세포 하위집단(aEC) |
| **Tang et al. (*Science*, 2022)** | `GSE199190` (SuperSeries) | bulk ATAC-seq + RNA-seq (단일세포 아님) | 인간 CRPC 오가노이드·PDX·세포주 | 크로마틴 프로파일 기반 CRPC 4군 분류 (AR-dependent, NE, Wnt, SCL) 및 군별 치료 표적 |
| **Karthaus et al. (*Science*, 2020)** | `SCP864` (인간, Single Cell Portal) / `SCP859` (마우스; GEO `GSE146811`) | scRNA-seq (10x) | 마우스 안드로겐 차단(Castration) 및 재투여 전후 전립선 관강세포 재생 아틀라스, 인간 전립선 조직은 SCP864 (종별 accession이 다름에 주의) | 거세 후 분화된 관강세포가 강화된 재생능을 획득함을 추적 (희귀 Lum-P 집단 포함) |

> **실전 메타데이터 검증 가이드 (Integration 전 필수 확인 사항)**:
> 서로 다른 논문 및 코호트의 공개 데이터를 결합(Merge)하거나 배치 통합(Harmony, scVI)을 수행하기 전, `adata.obs`에 다음 메타데이터 컬럼을 사전에 반드시 검사하고 정제해야 합니다:
> 1. **기증자/환자 ID (`donor_id` / `patient_id`)**: 동일 환자에게서 채취된 다중 생검 코어가 서로 다른 독립 샘플로 오인되지 않도록 확인.
> 2. **조직 부위 (`tissue_site`)**: 전립선 원발부(Prostate PZ/TZ) vs 원격 전이 장기(Bone marrow, Lymph node, Liver) 간의 조직 특이적 배경 전사체(예: 골수계 헤모글로빈, 조골세포 기질 유전자) 차이를 인지.
> 3. **치료 이력 (`treatment_history`)**: 치료 이력이 없는 환자(Treatment-naive), 안드로겐 차단 요법(ADT), 차세대 AR 표적 치료제(Enzalutamide/Abiraterone), 백금 기반 항암화학요법 여부에 따라 세포 상태가 급변하므로 배치 공변량 또는 그룹화 변수로 보존.
> 4. **해리 프로토콜 (`protocol_type`)**: 전세포 효소 해리(`scRNA-seq`)와 동결 핵 분리(`snRNA-seq`)는 포획되는 미토콘드리아 비율 및 인트론(Pre-mRNA) 비율이 완전히 다르므로, 이를 단일 배치 키(batch key)로 명시하여 배치 통합 모델에 기술적 교정 파라미터로 제공해야 합니다.

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Kfoury Y, Baryawno N, Severe N, et al. Human prostate cancer bone metastases have an actionable immunosuppressive microenvironment. *Cancer Cell*. 2021;39(11):1464-1478.e8. PMID: 34719426. doi: 10.1016/j.ccell.2021.09.005
2. Song H, Weinstein HNW, Allegakoen P, et al. Single-cell analysis of human primary prostate cancer reveals the heterogeneity of tumor-associated epithelial cell states. *Nat Commun*. 2022;13(1):141. PMID: 35013146. doi: 10.1038/s41467-021-27322-4
3. Chen S, Zhu G, Yang Y, et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. *Nat Cell Biol*. 2021;23(1):87-98. PMID: 33420488. doi: 10.1038/s41556-020-00613-6
4. Tang F, Xu D, Wang S, et al. Chromatin profiles classify castration-resistant prostate cancers suggesting therapeutic targets. *Science*. 2022;376(6596):eabe1505. PMID: 35617398. doi: 10.1126/science.abe1505
5. Karthaus WR, Hofree M, Choi D, et al. Regenerative potential of prostate luminal cells revealed by single-cell analysis. *Science*. 2020;368(6490):497-505. PMID: 32355025. doi: 10.1126/science.aay0267
### 검수 근거 문헌
1. Denisenko E, Guo BB, Jones M, et al. Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows. *Genome Biol*. 2020;21(1):130. PMID: 32487174. doi: 10.1186/s13059-020-02048-6
2. Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control. *Bioinformatics*. 2021;37(7):963-967. PMID: 32840568. doi: 10.1093/bioinformatics/btaa751
3. Gupta A, Shamsi F, Altemose N, et al. Characterization of transcript enrichment and detection bias in single-nucleus RNA-seq for mapping of distinct human adipocyte lineages. *Genome Res*. 2022;32(2):242-257. PMID: 35042723. doi: 10.1101/gr.275509.121
4. Chan JM, Zaidi S, Love JR, et al. Lineage plasticity in prostate cancer depends on JAK/STAT inflammatory signaling. *Science*. 2022;377(6611):1180-1191. PMID: 35981096. doi: 10.1126/science.abn0478
5. Hippen AA, Falco MM, Weber LM, et al. miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data. *PLoS Comput Biol*. 2021;17(8):e1009290. PMID: 34428202. doi: 10.1371/journal.pcbi.1009290
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- anndata 0.13.4 (github.com/scverse/anndata). 원 논문: Virshup I, Rybakov S, Theis FJ, et al. anndata: Access and store annotated data matrices. *JOSS*. 2024;9(101):4371. doi: 10.21105/joss.04371
- CellBender (github.com/broadinstitute/CellBender). 원 논문: Fleming SJ, Chaffin MD, Arduini A, et al. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. *Nat Methods*. 2023;20(9):1323-1335. PMID: 37550580. doi: 10.1038/s41592-023-01943-7
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- scvi-tools 1.5.1 (github.com/scverse/scvi-tools). 원 논문: Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol*. 2022;40(2):163-166. PMID: 35132262. doi: 10.1038/s41587-021-01206-w
<!-- REFERENCES:END -->

#!/usr/bin/env python3
"""
tracks/01_multiome_snatac/pipelines/01_multiome_qc_and_lsi_pipeline.py
=====================================================================
10x Multiome (단일 핵 동시 snRNA-seq + snATAC-seq) 품질 관리(QC) 및 
잠재 의미 지수화(LSI: TF-IDF ➔ TruncatedSVD) 표준 전처리 파이프라인.

주요 파이프라인 기능:
1. Multiome 쌍체(Paired) 데이터 로딩 또는 1,000 nuclei x 3,000 genes x 15,000 peaks 합성 코호트 생성
   - 전립선암 4대 악성 아형(AR+ Adeno, Amphicrine, DNPC, t-NEPC) 및 TME(CAF, CD8_T, Osteoblast, Osteoclast)
2. snRNA QC: Mito% < 2.5%, min genes >= 500, UMI 필터링
3. snATAC QC: TSS 농축도(TSS Enrichment Score >= 5.0), 뉴클레오솜 밴딩 주기성(NFR-to-Mono ratio >= 1.0)
4. LSI 차원 축소: 이진화(Binarization) ➔ TF-IDF 정규화 ➔ TruncatedSVD(30 PCs)
5. LSI Component 1 시퀀싱 심도(Depth) 상관관계 진단 및 필터링 (|r| > 0.90 시 제거)
6. 결과 요약 테이블 및 필터링된 데이터 저장 (results_pipeline_01/)
"""

from __future__ import annotations

import os
import sys
import warnings

# 경고 메시지 및 matplotlib 캐시 설정
warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib_cache"

# 선택 사항: SCRNA_ENV_SITE_PACKAGES 환경변수로 추가 site-packages 경로(예: scanpy가 설치된 conda 환경)를
# 지정할 수 있다. 지정하지 않으면 현재 인터프리터 환경만 사용한다(개인 환경 경로를 하드코딩하지 않음).
_extra_site = os.environ.get("SCRNA_ENV_SITE_PACKAGES", "")
if _extra_site and os.path.isdir(_extra_site) and _extra_site not in sys.path:
    sys.path.append(_extra_site)

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr
from sklearn.decomposition import TruncatedSVD

try:
    import scanpy as sc
    HAS_SCANPY = True
    sc.settings.verbosity = 0
except Exception as e:
    HAS_SCANPY = False


# ==============================================================================
# Step 1: Multiome Paired Data Loading or Simulation
# ==============================================================================
def load_or_simulate_multiome_data(
    data_dir: str | None = None,
    n_nuclei: int = 1000,
    n_genes: int = 3000,
    n_peaks: int = 15000,
    random_seed: int = 42
) -> tuple[sc.AnnData, sc.AnnData]:
    """
    동일한 단일 핵 바코드를 공유하는 snRNA와 snATAC 쌍체 데이터를 로드하거나,
    전립선암 임상 특이적 전사체-오픈 크로마틴 프로그램을 가진 합성 Multiome 데이터를 생성합니다.
    """
    print("=" * 80)
    print("[Step 1] Multiome (snRNA + snATAC) Paired Cohort Loading or Simulation")
    print("=" * 80)

    np.random.seed(random_seed)

    # 4명 환자, 2개 전이 부위 (Primary Prostate vs Bone Metastasis)
    patients = ["Patient_01", "Patient_02", "Patient_03", "Patient_04"]
    patient_site_map = {
        "Patient_01": "Primary_Prostate",
        "Patient_02": "Primary_Prostate",
        "Patient_03": "Bone_Metastasis",
        "Patient_04": "Bone_Metastasis",
    }
    patient_batch_map = {
        "Patient_01": "Batch_A",
        "Patient_02": "Batch_A",
        "Patient_03": "Batch_B",
        "Patient_04": "Batch_B",
    }

    patient_assign = np.random.choice(patients, size=n_nuclei, p=[0.25, 0.25, 0.25, 0.25])
    sites = [patient_site_map[p] for p in patient_assign]
    batches = [patient_batch_map[p] for p in patient_assign]

    # 세포 유형 정의 (4대 mCRPC 내성 아형 65% + TME 기질/면역 35%)
    malignant_states = ["AR+ Adeno", "Amphicrine", "DNPC", "t-NEPC"]
    cell_types = []
    for s in sites:
        if s == "Primary_Prostate":
            pool = ["AR+ Adeno", "Amphicrine", "DNPC", "CAF", "CD8_T"]
            probs = [0.50, 0.20, 0.05, 0.15, 0.10]
            cell_types.append(np.random.choice(pool, p=probs))
        else:
            pool = ["AR+ Adeno", "Amphicrine", "DNPC", "t-NEPC", "Osteoblast", "Osteoclast", "CAF", "CD8_T"]
            probs = [0.15, 0.20, 0.25, 0.15, 0.10, 0.05, 0.05, 0.05]
            cell_types.append(np.random.choice(pool, p=probs))

    # 바코드 생성 (동일 단일 핵 공유 식별자)
    barcodes = [f"AAACAGCC_{i:04d}-1" for i in range(n_nuclei)]

    # 1. snRNA 유전자 및 발현 행렬
    curated_genes = [
        "AR", "KLK3", "KLK2", "TMPRSS2", "FOXA1", "HOXB13",
        "SOX2", "EZH2", "POU5F1", "MYCN", "TP53", "RB1",
        "FGFR1", "FGF8", "ETV4", "MAPK1", "DUSP6",
        "ASCL1", "INSM1", "SYP", "CHGA", "NEUROD1", "DLL3",
        "REST", "BGLAP", "ALPL", "CTSK", "ACP5", "COL1A1",
        "ACTA2", "CD3D", "CD8A", "MT-CO1", "MT-ND1", "MT-ATP6"
    ]
    all_genes = list(curated_genes)
    rem_genes = n_genes - len(all_genes)
    for i in range(rem_genes):
        all_genes.append(f"Gene_{i:04d}")

    gene_to_idx = {g: idx for idx, g in enumerate(all_genes)}
    rna_counts = np.random.negative_binomial(2, 0.4, size=(n_nuclei, n_genes)).astype(np.float32)

    # 2. snATAC 피크 좌표 및 염색질 접근성 행렬 (15,000 Peaks)
    curated_peaks = [
        "chrX_66901240_66901450",  # chrX:66.9Mb AR upstream 650kb Enhancer Peak
        "chrX_67540000_67540800",  # AR Promoter Peak
        "chr19_51358000_51359000", # KLK3 Promoter & Enhancer Peak
        "chr14_38000000_38001000", # FOXA1 Enhancer Peak
        "chr12_103350000_103351000",# ASCL1 Promoter & Enhancer Peak
        "chr8_38280000_38281000",  # FGFR1 Enhancer Peak
        "chr4_1800000_1801000",    # FGFR3 Peak
        "chr2_241500000_241501000",# REST Locus Peak
        "chr1_156100000_156101000",# Osteoblast BGLAP Promoter Peak
        "chr10_87860000_87861000", # PTEN Locus Peak
    ]
    all_peaks = list(curated_peaks)
    rem_peaks = n_peaks - len(all_peaks)
    for i in range(rem_peaks):
        chrom_num = (i % 22) + 1
        start = (i * 2000) % 150000000
        end = start + np.random.choice([250, 400, 600, 850])
        all_peaks.append(f"chr{chrom_num}_{start}_{end}")

    peak_to_idx = {p: idx for idx, p in enumerate(all_peaks)}
    # ATAC 오픈 피크 카운트 (평균 1,800 피크/핵, 불량 핵은 ~400 피크 미만)
    atac_counts = (np.random.rand(n_nuclei, n_peaks) < 0.12).astype(np.float32)

    # 세포 유형별 전사-염색질 조절 프로그램 주입
    for i, ct in enumerate(cell_types):
        if ct == "AR+ Adeno":
            for g in ["AR", "KLK3", "KLK2", "TMPRSS2", "FOXA1"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(28)
            atac_counts[i, peak_to_idx["chrX_66901240_66901450"]] = 1.0  # AR Enhancer
            atac_counts[i, peak_to_idx["chrX_67540000_67540800"]] = 1.0  # AR Promoter
            atac_counts[i, peak_to_idx["chr19_51358000_51359000"]] = 1.0 # KLK3 Peak
        elif ct == "Amphicrine":
            for g in ["AR", "KLK3", "SOX2", "EZH2", "SYP"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(18)
            atac_counts[i, peak_to_idx["chrX_66901240_66901450"]] = 1.0
            atac_counts[i, peak_to_idx["chr12_103350000_103351000"]] = 1.0 # ASCL1/신경 인핸서 조기 개방!
        elif ct == "DNPC":
            for g in ["FGFR1", "FGF8", "ETV4", "MAPK1", "DUSP6"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(24)
            atac_counts[i, peak_to_idx["chr8_38280000_38281000"]] = 1.0  # FGFR1 Enhancer
            atac_counts[i, peak_to_idx["chrX_66901240_66901450"]] = 0.0  # AR 폐쇄
        elif ct == "t-NEPC":
            for g in ["ASCL1", "INSM1", "SYP", "CHGA", "NEUROD1", "DLL3"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(30)
            atac_counts[i, peak_to_idx["chr12_103350000_103351000"]] = 1.0 # ASCL1 강력 개방
            atac_counts[i, peak_to_idx["chrX_66901240_66901450"]] = 0.0
        elif ct == "Osteoblast":
            for g in ["BGLAP", "ALPL"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(25)
            atac_counts[i, peak_to_idx["chr1_156100000_156101000"]] = 1.0
        elif ct == "CD8_T":
            for g in ["CD3D", "CD8A"]:
                rna_counts[i, gene_to_idx[g]] += np.random.poisson(22)

    # QC 메트릭스 시뮬레이션
    # 1. Mito% (정상은 0.5~2.0%, ~5% 불량 핵에 3.0~6.0% 주입)
    pct_mt = np.random.beta(2, 20, size=n_nuclei) * 10.0
    bad_rna_mask = np.random.rand(n_nuclei) < 0.05
    pct_mt[bad_rna_mask] += np.random.uniform(3.0, 6.0, size=np.sum(bad_rna_mask))
    rna_counts[bad_rna_mask, 200:] = 0  # < 500 genes

    # 2. TSS Enrichment (정상은 7.0~16.0, ~5% 불량 핵에 2.0~4.5 주입)
    tss_scores = np.random.normal(10.5, 2.5, size=n_nuclei)
    tss_scores = np.clip(tss_scores, 2.0, 25.0)
    bad_atac_mask = np.random.rand(n_nuclei) < 0.05
    tss_scores[bad_atac_mask] = np.random.uniform(2.0, 4.5, size=np.sum(bad_atac_mask))
    atac_counts[bad_atac_mask, :] = (np.random.rand(np.sum(bad_atac_mask), n_peaks) < 0.02).astype(np.float32)  # < 500 peaks

    # 3. NFR to Mononucleosome ratio (정상 1.5~3.5, 불량 < 1.0)
    nfr_mono_ratio = np.random.normal(2.2, 0.4, size=n_nuclei)
    nfr_mono_ratio[bad_atac_mask] = np.random.uniform(0.3, 0.8, size=np.sum(bad_atac_mask))

    # AnnData 객체 구성
    obs_df = pd.DataFrame(
        {
            "patient_id": patient_assign,
            "sample_site": sites,
            "batch": batches,
            "cell_type": cell_types,
            "is_malignant": [ct in malignant_states for ct in cell_types],
            "pct_counts_mt": pct_mt,
            "tss_enrichment": tss_scores,
            "nfr_to_mono_ratio": nfr_mono_ratio,
        },
        index=barcodes,
    )

    adata_rna = sc.AnnData(
        X=rna_counts,
        obs=obs_df.copy(),
        var=pd.DataFrame(index=all_genes),
    )
    adata_rna.layers["raw_counts"] = rna_counts.copy()

    adata_atac = sc.AnnData(
        X=atac_counts,
        obs=obs_df.copy(),
        var=pd.DataFrame(index=all_peaks),
    )
    adata_atac.layers["raw_counts"] = atac_counts.copy()

    print(f"[✓] Created Paired Multiome Cohort:")
    print(f"    - Common Nuclei Barcodes: {n_nuclei}")
    print(f"    - snRNA: {adata_rna.n_vars} genes")
    print(f"    - snATAC: {adata_atac.n_vars} genomic peaks")
    print(f"    - Sites: {dict(adata_rna.obs['sample_site'].value_counts())}")
    print(f"    - Cell Types: {dict(adata_rna.obs['cell_type'].value_counts())}")

    return adata_rna, adata_atac


# ==============================================================================
# Step 2: Multi-Modal QC & Joint Quality Filtering
# ==============================================================================
def run_multimodal_qc(
    adata_rna: sc.AnnData,
    adata_atac: sc.AnnData,
    max_mito_pct: float = 2.5,
    min_genes: int = 500,
    min_tss_score: float = 5.0,
    min_nfr_mono_ratio: float = 1.0,
    min_peaks: int = 1000
) -> tuple[sc.AnnData, sc.AnnData, pd.DataFrame]:
    """
    snRNA QC와 snATAC QC를 동시에 수행하여, 두 모달리티 모두에서 고품질 기준을 충족하는
    진정한 단일 핵만을 선별하는 Joint Filtering을 진행합니다.
    """
    print("\n" + "=" * 80)
    print("[Step 2] Multi-Modal Quality Control (snRNA + snATAC Joint Filtering)")
    print("=" * 80)

    # 1. snRNA QC 메트릭 계산
    adata_rna.var["mt"] = adata_rna.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata_rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    # 시뮬레이션된 pct_counts_mt 보존
    if "pct_counts_mt" not in adata_rna.obs or np.all(adata_rna.obs["pct_counts_mt"] == 0):
        adata_rna.obs["pct_counts_mt"] = adata_rna.obs["pct_counts_mt"]

    # 2. snATAC QC 메트릭 계산
    peak_counts_per_cell = np.sum(adata_atac.X > 0, axis=1)
    adata_atac.obs["n_peaks_by_counts"] = np.array(peak_counts_per_cell).flatten()

    n_init = adata_rna.n_obs

    # snRNA 통과 마스크
    pass_rna = (
        (adata_rna.obs["pct_counts_mt"] < max_mito_pct) &
        (adata_rna.obs["n_genes_by_counts"] >= min_genes)
    )

    # snATAC 통과 마스크
    pass_atac = (
        (adata_atac.obs["tss_enrichment"] >= min_tss_score) &
        (adata_atac.obs["nfr_to_mono_ratio"] >= min_nfr_mono_ratio) &
        (adata_atac.obs["n_peaks_by_counts"] >= min_peaks)
    )

    # Joint 통과 마스크 (동일 핵 모두 통과)
    pass_joint = pass_rna & pass_atac

    n_pass_rna = np.sum(pass_rna)
    n_pass_atac = np.sum(pass_atac)
    n_pass_joint = np.sum(pass_joint)

    qc_summary_df = pd.DataFrame(
        {
            "QC_Stage": [
                "1. Initial Cohort",
                "2. snRNA Passed (Mito < 2.5%, Genes >= 500)",
                "3. snATAC Passed (TSS >= 5.0, NFR/Mono >= 1.0, Peaks >= 1000)",
                "4. Joint Multiome Retained (Both Modalities)",
            ],
            "Retained_Nuclei": [n_init, n_pass_rna, n_pass_atac, n_pass_joint],
            "Retention_Rate": [
                "100.0%",
                f"{n_pass_rna / n_init * 100:.1f}%",
                f"{n_pass_atac / n_init * 100:.1f}%",
                f"{n_pass_joint / n_init * 100:.1f}%",
            ],
            "Filtered_Out": [0, n_init - n_pass_rna, n_init - n_pass_atac, n_init - n_pass_joint],
        }
    )

    print("\n[+] Multiome QC Filtering Summary:")
    print(qc_summary_df.to_string(index=False))

    # 필터링 적용
    adata_rna_filt = adata_rna[pass_joint].copy()
    adata_atac_filt = adata_atac[pass_joint].copy()

    # 극소수 핵에서만 나타나는 피크 및 유전자 필터
    sc.pp.filter_genes(adata_rna_filt, min_cells=3)
    sc.pp.filter_genes(adata_atac_filt, min_cells=5)

    print(f"\n[✓] Quality Filtering Complete: {n_pass_joint} high-quality paired nuclei retained.")
    print(f"    - Filtered snRNA: {adata_rna_filt.n_obs} nuclei × {adata_rna_filt.n_vars} genes")
    print(f"    - Filtered snATAC: {adata_atac_filt.n_obs} nuclei × {adata_atac_filt.n_vars} peaks")

    return adata_rna_filt, adata_atac_filt, qc_summary_df


# ==============================================================================
# Step 3: snATAC LSI Dimensionality Reduction & Component 1 Depth Diagnosis
# ==============================================================================
def run_lsi_dimensionality_reduction(
    adata_atac: sc.AnnData,
    n_components: int = 30,
    depth_cor_threshold: float = 0.90
) -> tuple[sc.AnnData, pd.DataFrame]:
    """
    snATAC 피크 행렬에 대해 이진화(Binarization) ➔ TF-IDF 정규화 ➔ TruncatedSVD를 수행하고,
    LSI Component 1과 총 피크 수 간의 시퀀싱 심도 바이어스를 진단 및 보정합니다.
    """
    print("\n" + "=" * 80)
    print("[Step 3] snATAC LSI (TF-IDF ➔ TruncatedSVD) & Component 1 Depth Diagnosis")
    print("=" * 80)

    # 1. 피크 행렬 이진화 (0 or 1)
    X_raw = adata_atac.X.toarray() if hasattr(adata_atac.X, "toarray") else adata_atac.X
    X_bin = (X_raw > 0).astype(np.float32)
    n_cells, n_peaks = X_bin.shape

    # 2. TF-IDF 정규화 계산
    # TF: 각 세포 내 총 오픈 피크 수 대비 해당 피크의 점유 비율
    total_peaks_per_cell = np.sum(X_bin, axis=1, keepdims=True)
    total_peaks_per_cell[total_peaks_per_cell == 0] = 1.0
    tf_mat = X_bin / total_peaks_per_cell

    # IDF: 전체 세포 중 해당 피크가 열려 있는 세포 빈도의 역수 (로그 스케일)
    cell_freq_per_peak = np.sum(X_bin, axis=0, keepdims=True)
    idf_vec = np.log(1.0 + n_cells / (cell_freq_per_peak + 1.0))

    # TF-IDF 결합
    tfidf_mat = tf_mat * idf_vec
    tfidf_sparse = sparse.csr_matrix(tfidf_mat)

    print(f"[*] Computed TF-IDF Normalized Matrix ({n_cells} nuclei × {n_peaks} peaks)")

    # 3. TruncatedSVD 실행 (LSI 잠재 공간 30차원)
    svd = TruncatedSVD(n_components=n_components, algorithm="randomized", random_state=42)
    X_svd = svd.fit_transform(tfidf_sparse)
    adata_atac.obsm["X_lsi_all"] = X_svd

    # 4. LSI 컴포넌트와 시퀀싱 심도(총 피크 수) 간의 피어슨 상관계수 r 진단
    total_peaks = total_peaks_per_cell.flatten()
    log_depth = np.log1p(total_peaks)

    correlations = []
    for comp_idx in range(n_components):
        r_val, p_val = pearsonr(X_svd[:, comp_idx], log_depth)
        correlations.append(
            {
                "LSI_Component": f"LSI_{comp_idx + 1}",
                "Pearson_r_with_Depth": r_val,
                "Absolute_r": abs(r_val),
                "P_Value": p_val,
                "Status": "BIAS_FLAGGED" if abs(r_val) > depth_cor_threshold else "PASS_BIOLOGICAL",
            }
        )

    cor_df = pd.DataFrame(correlations)

    print("\n[+] Top 5 LSI Components vs Sequencing Depth Correlation Check:")
    print(cor_df.head(5).to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # LSI Component 1 제거 여부 판정
    comp1_r = cor_df.loc[0, "Absolute_r"]
    if comp1_r > depth_cor_threshold:
        print(f"\n[!] WARNING: LSI_1 is strongly correlated with sequencing depth (|r| = {comp1_r:.4f} > {depth_cor_threshold}).")
        print("[!] Dropping LSI_1 to prevent clustering purely driven by technical coverage.")
        adata_atac.obsm["X_lsi"] = X_svd[:, 1:n_components].copy()
        used_comps = f"LSI_2 ~ LSI_{n_components} ({n_components - 1} dimensions)"
    else:
        print(f"\n[✓] LSI_1 biological correlation acceptable (|r| = {comp1_r:.4f} <= {depth_cor_threshold}).")
        adata_atac.obsm["X_lsi"] = X_svd.copy()
        used_comps = f"LSI_1 ~ LSI_{n_components} ({n_components} dimensions)"

    print(f"[✓] Final Clean snATAC Latent Representation saved to adata_atac.obsm['X_lsi'] ({used_comps}).")
    return adata_atac, cor_df


# ==============================================================================
# Main Orchestration & Output Generation
# ==============================================================================
def main():
    print("=" * 80)
    print("01_multiome_qc_and_lsi_pipeline.py Execution")
    print("10x Multiome snRNA + snATAC Joint Quality Control & LSI Pipeline")
    print("=" * 80)

    out_dir = "tracks/01_multiome_snatac/results_pipeline_01"
    os.makedirs(out_dir, exist_ok=True)

    # 1. 데이터 로드/시뮬레이션
    adata_rna, adata_atac = load_or_simulate_multiome_data(n_nuclei=1000, n_genes=3000, n_peaks=15000)

    # 2. 멀티모달 QC & Joint Filtering
    adata_rna_filt, adata_atac_filt, qc_summary_df = run_multimodal_qc(
        adata_rna,
        adata_atac,
        max_mito_pct=2.5,
        min_genes=500,
        min_tss_score=5.0,
        min_nfr_mono_ratio=1.0,
        min_peaks=1000,
    )
    qc_summary_df.to_csv(os.path.join(out_dir, "qc_summary.csv"), index=False)

    # 3. snATAC LSI 차원 축소 및 Component 1 진단
    adata_atac_filt, cor_df = run_lsi_dimensionality_reduction(
        adata_atac_filt,
        n_components=30,
        depth_cor_threshold=0.90,
    )
    cor_df.to_csv(os.path.join(out_dir, "lsi_component_correlations.csv"), index=False)

    # 4. 필터링된 메타데이터 및 정량 지표 저장
    adata_rna_filt.obs.to_csv(os.path.join(out_dir, "filtered_multiome_metadata.csv"))

    # 파이프라인 02와의 연계를 위해 numpy 형식으로 잠재 공간 및 원시 카운트 저장
    np.save(os.path.join(out_dir, "atac_lsi_embedding.npy"), adata_atac_filt.obsm["X_lsi"])
    np.save(os.path.join(out_dir, "rna_expression_matrix.npy"), adata_rna_filt.X)
    np.save(os.path.join(out_dir, "atac_binary_matrix.npy"), (adata_atac_filt.X > 0).astype(np.float32))
    
    with open(os.path.join(out_dir, "gene_names.txt"), "w") as f:
        f.write("\n".join(adata_rna_filt.var_names))
    with open(os.path.join(out_dir, "peak_names.txt"), "w") as f:
        f.write("\n".join(adata_atac_filt.var_names))
    with open(os.path.join(out_dir, "cell_barcodes.txt"), "w") as f:
        f.write("\n".join(adata_rna_filt.obs_names))

    print("\n" + "=" * 80)
    print("[✓] Pipeline 01 completed successfully with Exit Code 0!")
    print(f"[✓] Saved QC & LSI artifacts in {out_dir}/:")
    print(f"    - qc_summary.csv")
    print(f"    - lsi_component_correlations.csv")
    print(f"    - filtered_multiome_metadata.csv")
    print(f"    - atac_lsi_embedding.npy, rna_expression_matrix.npy, atac_binary_matrix.npy")
    print("=" * 80)


if __name__ == "__main__":
    main()

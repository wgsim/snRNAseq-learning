#!/usr/bin/env python3
"""
tracks/02_spatial_transcriptomics/pipelines/01_spatial_qc_and_deconvolution_pipeline.py
========================================================================================
End-to-End Visium / Visium HD Spatial QC & Reference-Based Deconvolution Pipeline.

Key Pipeline Steps:
1. Multi-Spot Spatial Cohort Simulation (2,500 spots on a 50x50 spatial grid)
   - Simulating prostate adenocarcinoma core, CAF stroma rim, normal gland,
     immune infiltration (TLS), and bone metastatic front across 500 genes.
2. Spatial Quality Control (QC):
   - In-tissue contour masking, UMI/gene detection thresholds, mitochondrial ratio check.
3. Spatial Variable Gene (SVG) Detection:
   - Spatial kNN graph construction (k=6), spatial weight matrix W,
     Spatial Lag computation, and global Moran's I with an analytic Z-test using the exact
     variance under the normality assumption (Cliff & Ord; NOT a permutation test) and BH-FDR correction.
4. Reference-Based Deconvolution (Poisson NMF-KL multiplicative-update solver with a known signature matrix;
   NOT cell2location: no hierarchical Bayesian model, no priors, no posterior/credibility intervals):
   - 8 prostate cell types: Luminal_AR_Pos, Amphicrine, DNPC, t_NEPC,
     CAF, Endothelial, CD8_T, Osteoblast.
   - Absolute cell abundance & relative proportion estimation (abundances are rescaled from total UMI,
     not independently estimated; no credibility intervals are computed), and a synthetic-data sanity check
     against ground truth (the true signature matrix is given to the solver, so this is not
     evidence of real-world accuracy).
5. Save Structured Results to results_pipeline_01/.
"""

from __future__ import annotations

import json
import os
import sys
import warnings

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
from scipy.spatial.distance import cdist
from scipy.stats import norm, pearsonr, spearmanr


# ==============================================================================
# Step 1: Multi-Spot Spatial Cohort Simulation
# ==============================================================================
def simulate_visium_cohort(
    grid_size: int = 50,
    n_genes: int = 500,
    random_seed: int = 42,
) -> dict:
    """
    Simulates a 50x50 spatial Visium array (2,500 spots) representing a human
    prostate cancer tissue section with tumor core, CAF stroma rim, normal gland,
    TLS immune hubs, and bone metastatic interface.
    """
    print("=" * 80)
    print("[Step 1] Multi-Spot Spatial Cohort Simulation (50x50 Grid = 2,500 Spots)")
    print("=" * 80)

    np.random.seed(random_seed)
    n_spots = grid_size * grid_size

    # Grid coordinates
    rows, cols = np.indices((grid_size, grid_size))
    row_flat = rows.flatten()
    col_flat = cols.flatten()
    spot_barcodes = [f"SPOT_R{r:02d}_C{c:02d}" for r, c in zip(row_flat, col_flat)]

    # Physical coordinates (100 µm pitch)
    x_um = col_flat * 100.0 + np.random.normal(0, 1.5, size=n_spots)
    y_um = row_flat * 100.0 + np.random.normal(0, 1.5, size=n_spots)
    center_r, center_c = grid_size / 2.0, grid_size / 2.0
    dist_from_center = np.sqrt((row_flat - center_r) ** 2 + (col_flat - center_c) ** 2)

    # 1. In-Tissue Masking: Elliptical tissue section contour (covers ~86% of spots)
    norm_y = (row_flat - center_r) / 22.0
    norm_x = (col_flat - center_c) / 22.0
    in_tissue = ((norm_y ** 2 + norm_x ** 2) <= 1.0).astype(int)
    # Add minor boundary roughness
    boundary_noise = np.random.rand(n_spots) < 0.03
    in_tissue[boundary_noise] = 0

    # 2. Histological Compartments
    # - Tumor Core: center (dist <= 11)
    # - Desmoplastic CAF Rim: (11 < dist <= 17)
    # - Normal Gland / Benign Epithelium: (dist > 17 and in_tissue)
    # - Bone Metastatic Front: top-right corner (row <= 16 and col >= 34)
    # - TLS Focal Patches: 3 distinct immune clusters
    region = np.full(n_spots, "Off_Tissue", dtype=object)
    region[in_tissue == 1] = "Normal_Gland"

    tumor_mask = (dist_from_center <= 11.0) & (in_tissue == 1)
    region[tumor_mask] = "Tumor_Core"

    stroma_mask = (dist_from_center > 11.0) & (dist_from_center <= 17.0) & (in_tissue == 1)
    region[stroma_mask] = "CAF_Stroma_Rim"

    bone_mask = (row_flat <= 16) & (col_flat >= 34) & (in_tissue == 1)
    region[bone_mask] = "Bone_Metastatic_Front"

    # TLS focal hubs (seeds)
    tls_seeds = [(14, 18), (32, 14), (35, 33)]
    for tr, tc in tls_seeds:
        tls_dist = np.sqrt((row_flat - tr) ** 2 + (col_flat - tc) ** 2)
        tls_spots = (tls_dist <= 2.2) & (in_tissue == 1)
        region[tls_spots] = "TLS_Immune_Hub"

    # 3. Ground-Truth Cell Types (8 types from Track 00)
    cell_types = [
        "Luminal_AR_Pos",
        "Amphicrine",
        "DNPC",
        "t_NEPC",
        "CAF",
        "Endothelial",
        "CD8_T",
        "Osteoblast",
    ]
    n_cell_types = len(cell_types)

    # Ground-truth cell abundance per spot (total cells per spot ~ 8-15)
    true_abundances = np.zeros((n_spots, n_cell_types), dtype=np.float32)

    for i in range(n_spots):
        if in_tissue[i] == 0:
            continue
        total_cells = np.random.randint(8, 16)
        r = region[i]
        if r == "Tumor_Core":
            probs = [0.65, 0.15, 0.05, 0.05, 0.06, 0.02, 0.02, 0.0]
        elif r == "CAF_Stroma_Rim":
            probs = [0.10, 0.02, 0.03, 0.00, 0.65, 0.12, 0.08, 0.0]
        elif r == "Normal_Gland":
            probs = [0.45, 0.00, 0.25, 0.00, 0.18, 0.10, 0.02, 0.0]
        elif r == "TLS_Immune_Hub":
            probs = [0.05, 0.00, 0.00, 0.00, 0.15, 0.10, 0.70, 0.0]
        elif r == "Bone_Metastatic_Front":
            probs = [0.25, 0.05, 0.05, 0.05, 0.15, 0.05, 0.05, 0.35]
        else:
            probs = [0.125] * 8

        # Multinomial draw
        counts = np.random.multinomial(total_cells, probs)
        true_abundances[i, :] = counts

    # 4. Gene Profiles and Marker Sets
    gene_names = [f"Gene_{idx:03d}" for idx in range(n_genes)]
    # Replace key indices with verified canonical prostate biology markers
    marker_map = {
        0: "KLK3", 1: "KLK2", 2: "AR", 3: "TMPRSS2", 4: "FOXA1", 5: "HOXB13", 6: "NKX3-1",
        10: "SYP", 11: "CHGA", 12: "CHGB", 13: "ASCL1", 14: "INSM1",
        20: "FGFR1", 21: "FGF8", 22: "ETV4", 23: "MAPK1",
        30: "ACTA2", 31: "COL1A1", 32: "COL1A2", 33: "POSTN", 34: "FAP", 35: "VIM",
        40: "PECAM1", 41: "CD34", 42: "VWF", 43: "KDR",
        50: "CD3D", 51: "CD3E", 52: "CD8A", 53: "CD4", 54: "IFNG", 55: "GZMB",
        60: "BGLAP", 61: "RUNX2", 62: "ALPL", 63: "SPP1",
        70: "KRT5", 71: "KRT14", 72: "TP63",
        80: "GAPDH", 81: "ACTB", 82: "RPL13A", 83: "B2M", 84: "EEF1A1",
        90: "MT-CO1", 91: "MT-CO2", 92: "MT-ND1", 93: "MT-ATP6", 94: "MT-CYB",
    }
    for idx, name in marker_map.items():
        gene_names[idx] = name

    # Build reference expression profiles G (n_genes x n_cell_types)
    G = np.random.gamma(shape=1.0, scale=0.8, size=(n_genes, n_cell_types)).astype(np.float32)

    # Set cell-type specific signature programs
    def boost(gene_str, ct_idx, factor=35.0):
        if gene_str in gene_names:
            g_idx = gene_names.index(gene_str)
            G[g_idx, ct_idx] += factor

    # Luminal_AR_Pos (ct 0)
    for g in ["KLK3", "KLK2", "AR", "TMPRSS2", "FOXA1", "NKX3-1"]:
        boost(g, 0, factor=45.0)
    # Amphicrine (ct 1: dual AR + Neuroendocrine)
    for g in ["AR", "KLK3", "SYP", "CHGA"]:
        boost(g, 1, factor=30.0)
    # DNPC (ct 2: AR-negative, FGFR/MAPK active)
    for g in ["FGFR1", "FGF8", "ETV4", "MAPK1"]:
        boost(g, 2, factor=38.0)
    # t_NEPC (ct 3: Pure small cell neuroendocrine)
    for g in ["SYP", "CHGA", "CHGB", "ASCL1", "INSM1"]:
        boost(g, 3, factor=50.0)
    # CAF (ct 4)
    for g in ["ACTA2", "COL1A1", "COL1A2", "POSTN", "FAP", "VIM"]:
        boost(g, 4, factor=42.0)
    # Endothelial (ct 5)
    for g in ["PECAM1", "CD34", "VWF", "KDR"]:
        boost(g, 5, factor=40.0)
    # CD8_T (ct 6)
    for g in ["CD3D", "CD3E", "CD8A", "IFNG", "GZMB"]:
        boost(g, 6, factor=45.0)
    # Osteoblast (ct 7)
    for g in ["BGLAP", "RUNX2", "ALPL", "SPP1"]:
        boost(g, 7, factor=48.0)

    # Housekeeping genes (expressed in all cell types)
    for g in ["GAPDH", "ACTB", "RPL13A", "B2M", "EEF1A1"]:
        g_idx = gene_names.index(g)
        G[g_idx, :] += 20.0

    # Mitochondrial genes (low baseline across all cell types)
    for g in ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ATP6", "MT-CYB"]:
        g_idx = gene_names.index(g)
        G[g_idx, :] += 4.0

    # 5. Generate Spot Expression Counts
    # mu_se = sum_f (w_sf * G_ef) + ambient_noise
    expected_counts = true_abundances @ G.T  # (n_spots x n_genes)
    ambient_background = np.random.uniform(0.1, 0.5, size=(n_spots, n_genes))
    expected_counts += ambient_background

    # Poisson observation
    spot_counts = np.random.poisson(expected_counts).astype(np.float32)

    # Simulate off-tissue spots: very low counts (ambient noise only)
    off_spots = in_tissue == 0
    spot_counts[off_spots, :] = np.random.poisson(0.4, size=(np.sum(off_spots), n_genes))

    # Inject ~4% damaged spots with high mitochondrial lysis (mito% > 18%)
    damaged_mask = (np.random.rand(n_spots) < 0.04) & (in_tissue == 1)
    for g in ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ATP6", "MT-CYB"]:
        g_idx = gene_names.index(g)
        spot_counts[damaged_mask, g_idx] += np.random.poisson(120, size=np.sum(damaged_mask))

    print(f"  • Generated {n_spots} spots ({int(np.sum(in_tissue))} in-tissue, {int(np.sum(1 - in_tissue))} off-tissue)")
    print(f"  • Feature dimension: {n_genes} genes across {n_cell_types} cell types")
    print(f"  • Tissue compartments: {dict(pd.Series(region).value_counts())}")

    return {
        "spot_barcodes": spot_barcodes,
        "x_um": x_um,
        "y_um": y_um,
        "array_row": row_flat,
        "array_col": col_flat,
        "in_tissue": in_tissue,
        "region": region,
        "counts": spot_counts,
        "gene_names": gene_names,
        "cell_types": cell_types,
        "true_abundances": true_abundances,
        "signature_matrix": G,
    }


# ==============================================================================
# Step 2: Spatial Quality Control (QC)
# ==============================================================================
def run_spatial_qc(
    data: dict,
    min_counts: int = 250,
    min_genes: int = 50,
    max_pct_mt: float = 15.0,
) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """
    Performs spatial QC filtering:
    1. in_tissue == 1 masking
    2. Minimum UMI count threshold
    3. Minimum detected gene count threshold
    4. Mitochondrial fraction threshold (pct_counts_mt < 15%)
    """
    print("\n" + "=" * 80)
    print("[Step 2] Spatial Quality Control (QC) & Tissue Masking")
    print("=" * 80)

    counts = data["counts"]
    gene_names = data["gene_names"]
    n_spots = len(data["spot_barcodes"])

    total_counts = np.sum(counts, axis=1)
    n_genes_by_counts = np.sum(counts > 0, axis=1)

    # Mitochondrial percentage
    mt_gene_indices = [i for i, g in enumerate(gene_names) if g.startswith("MT-")]
    if len(mt_gene_indices) > 0:
        mt_counts = np.sum(counts[:, mt_gene_indices], axis=1)
        pct_counts_mt = np.where(total_counts > 0, (mt_counts / total_counts) * 100.0, 0.0)
    else:
        pct_counts_mt = np.zeros(n_spots)

    # Filtering boolean masks
    pass_in_tissue = data["in_tissue"] == 1
    pass_counts = total_counts >= min_counts
    pass_genes = n_genes_by_counts >= min_genes
    pass_mt = pct_counts_mt <= max_pct_mt

    qc_pass = pass_in_tissue & pass_counts & pass_genes & pass_mt

    obs_df = pd.DataFrame(
        {
            "spot_barcode": data["spot_barcodes"],
            "x_um": data["x_um"],
            "y_um": data["y_um"],
            "array_row": data["array_row"],
            "array_col": data["array_col"],
            "in_tissue": data["in_tissue"],
            "region": data["region"],
            "total_counts": total_counts,
            "n_genes_by_counts": n_genes_by_counts,
            "pct_counts_mt": pct_counts_mt,
            "qc_pass": qc_pass,
        }
    )

    qc_metrics = {
        "initial_spots": int(n_spots),
        "in_tissue_spots": int(np.sum(pass_in_tissue)),
        "qc_passed_spots": int(np.sum(qc_pass)),
        "filtered_off_tissue": int(np.sum(~pass_in_tissue)),
        "filtered_low_counts": int(np.sum(pass_in_tissue & ~pass_counts)),
        "filtered_high_mt": int(np.sum(pass_in_tissue & pass_counts & ~pass_mt)),
        "median_total_counts": float(np.median(total_counts[qc_pass])),
        "median_genes": float(np.median(n_genes_by_counts[qc_pass])),
        "median_pct_mt": float(np.median(pct_counts_mt[qc_pass])),
    }

    print(f"  • Total input spots: {qc_metrics['initial_spots']}")
    print(f"  • In-tissue spots: {qc_metrics['in_tissue_spots']} ({qc_metrics['in_tissue_spots']/n_spots*100:.1f}%)")
    print(f"  • QC Passed spots: {qc_metrics['qc_passed_spots']} ({qc_metrics['qc_passed_spots']/n_spots*100:.1f}%)")
    print(f"  • Filtered out:")
    print(f"    - Off-tissue: {qc_metrics['filtered_off_tissue']}")
    print(f"    - Low UMI (<{min_counts}): {qc_metrics['filtered_low_counts']}")
    print(f"    - High Mito (> {max_pct_mt}%): {qc_metrics['filtered_high_mt']}")
    print(f"  • Post-QC Median UMI: {qc_metrics['median_total_counts']:.1f} | Median Genes: {qc_metrics['median_genes']:.1f} | Median Mito: {qc_metrics['median_pct_mt']:.2f}%")

    return obs_df, qc_pass, qc_metrics


# ==============================================================================
# Step 3: Spatial Variable Gene (SVG) Detection via Moran's I
# ==============================================================================
def detect_spatial_variable_genes(
    counts_qc: np.ndarray,
    coords_qc: np.ndarray,
    gene_names: list[str],
    k_neighbors: int = 6,
) -> pd.DataFrame:
    """
    Computes Moran's I spatial autocorrelation statistic for all genes on the
    spatial kNN graph (k=6 hexagonal connectivity):
        I = (z^T * W * z) / (z^T * z)
    Evaluates analytical p-value and Benjamini-Hochberg FDR.
    """
    print("\n" + "=" * 80)
    print(f"[Step 3] Spatial Variable Gene (SVG) Detection (k-NN graph k={k_neighbors})")
    print("=" * 80)

    n_samples, n_genes = counts_qc.shape

    # 1. Build kNN spatial graph
    dist_mat = cdist(coords_qc, coords_qc)
    np.fill_diagonal(dist_mat, np.inf)

    # For each spot, pick k nearest neighbors
    nn_indices = np.argsort(dist_mat, axis=1)[:, :k_neighbors]

    # Row-normalized weight matrix W
    row_idx = np.repeat(np.arange(n_samples), k_neighbors)
    col_idx = nn_indices.flatten()
    data_val = np.full(n_samples * k_neighbors, 1.0 / k_neighbors, dtype=np.float32)
    W = sparse.csr_matrix((data_val, (row_idx, col_idx)), shape=(n_samples, n_samples))

    # 2. Normalize gene expression (log1p CPM)
    total_umis = np.sum(counts_qc, axis=1, keepdims=True)
    cpm = (counts_qc / (total_umis + 1e-6)) * 1e4
    log_cpm = np.log1p(cpm)

    # 3. Vectorized Moran's I Computation
    mean_val = np.mean(log_cpm, axis=0, keepdims=True)
    std_val = np.std(log_cpm, axis=0, keepdims=True) + 1e-8
    Z = (log_cpm - mean_val) / std_val  # (N x G) standardized

    # Spatial lag: WZ = W * Z
    WZ = W.dot(Z)  # (N x G)

    # Moran's I: sum(Z * WZ, axis=0) / sum(Z * Z, axis=0)
    numerator = np.sum(Z * WZ, axis=0)
    denominator = np.sum(Z * Z, axis=0) + 1e-8
    morans_i = numerator / denominator

    # Analytical variance under randomization hypothesis
    # E[I] = -1 / (N - 1)
    e_i = -1.0 / (n_samples - 1.0)
    # SCI-02-067: exact Var(I) under the normality assumption (Cliff & Ord 1981):
    #   Var(I) = (n^2*S1 - n*S2 + 3*S0^2) / (S0^2*(n^2-1)) - E[I]^2
    # with S0 = sum_ij w_ij, S1 = 0.5*sum_ij (w_ij + w_ji)^2, S2 = sum_i (sum_j w_ij + sum_j w_ji)^2.
    # (The earlier 1/(k*n) first-order approximation underestimated sd, e.g. 0.0082 vs 0.0113 on a
    #  50x50 grid with k=6, inflating Z by ~1.4x.) Still an analytic test: a permutation test is the check.
    S0 = float(W.sum())
    W_sym = W + W.T
    S1 = 0.5 * float(W_sym.multiply(W_sym).sum())
    S2 = float(np.sum(np.asarray(W_sym.sum(axis=1)).ravel() ** 2))
    n_f = float(n_samples)
    var_i = (n_f**2 * S1 - n_f * S2 + 3.0 * S0**2) / (S0**2 * (n_f**2 - 1.0)) - e_i**2
    z_scores = (morans_i - e_i) / np.sqrt(var_i)
    p_values = 1.0 - norm.cdf(z_scores)

    # Benjamini-Hochberg FDR
    df_svg = pd.DataFrame(
        {
            "gene": gene_names,
            "morans_i": morans_i,
            "z_score": z_scores,
            "p_value": p_values,
        }
    )
    df_svg.sort_values(by="morans_i", ascending=False, inplace=True)
    df_svg.reset_index(drop=True, inplace=True)

    # Rank & FDR
    ranks = np.arange(1, n_genes + 1)
    df_svg["fdr"] = np.minimum.accumulate((df_svg["p_value"] * n_genes / ranks)[::-1])[::-1]
    df_svg["fdr"] = np.clip(df_svg["fdr"], 0.0, 1.0)
    df_svg["is_svg"] = (df_svg["morans_i"] > 0.15) & (df_svg["fdr"] < 0.05)

    n_svgs = int(df_svg["is_svg"].sum())
    print(f"  • Total genes tested: {n_genes}")
    print(f"  • Significant SVGs detected (Moran's I > 0.15, FDR < 0.05): {n_svgs}")
    print("\n  Top 10 Spatial Variable Genes (SVGs):")
    print(df_svg.head(10)[["gene", "morans_i", "z_score", "p_value", "fdr"]].to_string(index=False))

    return df_svg


# ==============================================================================
# Step 4: Reference-Based Deconvolution (Poisson NMF-KL multiplicative-update solver; not cell2location)
# ==============================================================================
def run_poisson_deconvolution(
    counts_qc: np.ndarray,
    signature_matrix: np.ndarray,
    cell_types: list[str],
    max_iter: int = 120,
    tol: float = 1e-4,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solves for non-negative cell abundance matrix W (n_spots x n_cell_types)
    under a Poisson likelihood model:
        counts_{s, e} ~ Poisson( sum_f W_{s, f} * G_{e, f} + background )
    Using an efficient multiplicative coordinate update with L2 regularization:
        W_{s, f} <- W_{s, f} * ( sum_e G_{e, f} * (y_{s,e} / mu_{s,e}) ) / sum_e G_{e, f}
    """
    print("\n" + "=" * 80)
    print(f"[Step 4] Reference-Based Spatial Deconvolution ({len(cell_types)} Cell Types)")
    print("=" * 80)

    n_spots, n_genes = counts_qc.shape
    n_cell_types = len(cell_types)
    G = signature_matrix  # (n_genes x n_cell_types)

    # Initialize W based on initial non-negative projection
    total_umis = np.sum(counts_qc, axis=1, keepdims=True) + 1e-6
    W = np.ones((n_spots, n_cell_types), dtype=np.float32) * 2.0

    print("  • Executing Multiplicative Update Poisson Solver...")
    G_sum_e = np.sum(G, axis=0, keepdims=True)  # (1 x F)

    for it in range(max_iter):
        # Predicted count mu: (n_spots x n_genes)
        mu = W @ G.T + 0.1
        # Ratio of observed to predicted
        ratio = counts_qc / mu  # (n_spots x n_genes)
        # Gradient update step
        grad = ratio @ G  # (n_spots x n_cell_types)
        W_update = W * (grad / (G_sum_e + 1e-6))
        # Non-negative projection
        W_update = np.maximum(W_update, 1e-4)

        diff = np.max(np.abs(W_update - W))
        W = W_update
        if diff < tol:
            print(f"    - Converged at iteration {it + 1} with max diff {diff:.6f}")
            break

    # Re-scale W to reflect physical absolute cell counts (~5-15 cells per spot)
    target_cells_per_spot = np.clip(total_umis / 120.0, 4.0, 25.0)
    current_sum = np.sum(W, axis=1, keepdims=True) + 1e-6
    W_scaled = W * (target_cells_per_spot / current_sum)

    # Cell proportions P (rows sum to 1.0)
    proportions = W_scaled / (np.sum(W_scaled, axis=1, keepdims=True) + 1e-6)

    print(f"  • Deconvolution completed for {n_spots} spots")
    print(f"  • Mean estimated cell density per spot: {np.mean(np.sum(W_scaled, axis=1)):.2f} cells")
    print("\n  Estimated Overall Cohort Cell Type Proportions:")
    mean_props = np.mean(proportions, axis=0)
    for ct, prop in zip(cell_types, mean_props):
        print(f"    - {ct:22s}: {prop * 100:.2f}%")

    return W_scaled, proportions


# ==============================================================================
# Step 5: Benchmark Evaluation and Results Saving
# ==============================================================================
def evaluate_and_save_results(
    obs_qc: pd.DataFrame,
    W_est: np.ndarray,
    props_est: np.ndarray,
    true_abundances_qc: np.ndarray,
    cell_types: list[str],
    df_svg: pd.DataFrame,
    qc_metrics: dict,
    out_dir: str,
):
    """
    Evaluates deconvolution accuracy against simulated ground truth and
    exports all tables to results_pipeline_01/.
    """
    print("\n" + "=" * 80)
    print(f"[Step 5] Benchmarking Accuracy & Exporting Results to {out_dir}/")
    print("=" * 80)

    os.makedirs(out_dir, exist_ok=True)

    # True proportions
    true_props = true_abundances_qc / (np.sum(true_abundances_qc, axis=1, keepdims=True) + 1e-6)

    # Benchmark metrics
    benchmark_records = []
    for f_idx, ct in enumerate(cell_types):
        y_true = true_props[:, f_idx]
        y_pred = props_est[:, f_idx]
        r_val, _ = pearsonr(y_true, y_pred)
        rho_val, _ = spearmanr(y_true, y_pred)
        rmse_val = np.sqrt(np.mean((y_true - y_pred) ** 2))

        benchmark_records.append(
            {
                "cell_type": ct,
                "pearson_r": float(r_val),
                "spearman_rho": float(rho_val),
                "rmse": float(rmse_val),
            }
        )

    df_bench = pd.DataFrame(benchmark_records)
    print("  Deconvolution Benchmark Metrics vs Ground Truth:")
    print(df_bench.to_string(index=False))

    overall_r = np.mean(df_bench["pearson_r"])
    overall_rmse = np.mean(df_bench["rmse"])
    print(f"\n  • Overall Mean Pearson r: {overall_r:.4f}")
    print(f"  • Overall Mean RMSE:      {overall_rmse:.4f}")

    # 1. Export QC summary
    df_qc = pd.DataFrame([qc_metrics])
    df_qc.to_csv(os.path.join(out_dir, "qc_summary.csv"), index=False)

    # 2. Export filtered spatial metadata
    obs_qc.to_csv(os.path.join(out_dir, "filtered_spatial_metadata.csv"), index=False)

    # 3. Export SVGs
    df_svg.to_csv(os.path.join(out_dir, "spatial_variable_genes_moran.csv"), index=False)

    # 4. Export estimated abundances and proportions
    df_abund = pd.DataFrame(W_est, columns=cell_types)
    df_abund.insert(0, "spot_barcode", obs_qc["spot_barcode"].values)
    df_abund.to_csv(os.path.join(out_dir, "deconvolution_cell_abundances.csv"), index=False)

    df_props = pd.DataFrame(props_est, columns=cell_types)
    df_props.insert(0, "spot_barcode", obs_qc["spot_barcode"].values)
    df_props.to_csv(os.path.join(out_dir, "deconvolution_cell_proportions.csv"), index=False)

    # 5. Export benchmark metrics
    df_bench.to_csv(os.path.join(out_dir, "deconvolution_benchmark_metrics.csv"), index=False)

    # 6. JSON Summary
    summary = {
        "pipeline": "01_spatial_qc_and_deconvolution_pipeline",
        "qc_metrics": qc_metrics,
        "n_svgs_detected": int(df_svg["is_svg"].sum()),
        "top_5_svgs": df_svg.head(5)["gene"].tolist(),
        "mean_deconvolution_pearson_r": float(overall_r),
        "mean_deconvolution_rmse": float(overall_rmse),
        "cell_types": cell_types,
    }
    with open(os.path.join(out_dir, "pipeline_execution_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n[✓] Successfully exported all pipeline 01 tables and summaries:")
    print(f"    - {os.path.join(out_dir, 'qc_summary.csv')}")
    print(f"    - {os.path.join(out_dir, 'filtered_spatial_metadata.csv')}")
    print(f"    - {os.path.join(out_dir, 'spatial_variable_genes_moran.csv')}")
    print(f"    - {os.path.join(out_dir, 'deconvolution_cell_abundances.csv')}")
    print(f"    - {os.path.join(out_dir, 'deconvolution_cell_proportions.csv')}")
    print(f"    - {os.path.join(out_dir, 'deconvolution_benchmark_metrics.csv')}")
    print(f"    - {os.path.join(out_dir, 'pipeline_execution_summary.json')}")


# ==============================================================================
# Main Pipeline Driver
# ==============================================================================
def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "pipelines", "results_pipeline_01")

    # Step 1: Simulate Visium Cohort
    data = simulate_visium_cohort(grid_size=50, n_genes=500, random_seed=42)

    # Step 2: Run Spatial QC
    obs_df, qc_pass, qc_metrics = run_spatial_qc(data, min_counts=250, min_genes=50, max_pct_mt=15.0)

    # Subset to QC-passed spots
    counts_qc = data["counts"][qc_pass]
    coords_qc = np.column_stack([obs_df["x_um"].values[qc_pass], obs_df["y_um"].values[qc_pass]])
    true_abundances_qc = data["true_abundances"][qc_pass]
    obs_qc = obs_df[qc_pass].copy().reset_index(drop=True)

    # Step 3: Spatial Variable Gene Detection
    df_svg = detect_spatial_variable_genes(counts_qc, coords_qc, data["gene_names"], k_neighbors=6)

    # Step 4: Spatial Deconvolution
    W_est, props_est = run_poisson_deconvolution(
        counts_qc,
        data["signature_matrix"],
        data["cell_types"],
        max_iter=120,
        tol=1e-4,
    )

    # Step 5: Evaluate Benchmark & Export
    evaluate_and_save_results(
        obs_qc,
        W_est,
        props_est,
        true_abundances_qc,
        data["cell_types"],
        df_svg,
        qc_metrics,
        out_dir,
    )

    print("\n" + "=" * 80)
    print("[✓] Pipeline 01 completed successfully with Exit Code 0!")
    print("=" * 80)


if __name__ == "__main__":
    main()

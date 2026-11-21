#!/usr/bin/env python3
"""
tracks/02_spatial_transcriptomics/pipelines/02_xenium_subcellular_niche_and_signaling_pipeline.py
================================================================================================
End-to-End In Situ Subcellular Resolution (Xenium-style) Niche & Spatial Signaling Pipeline.

Key Pipeline Steps:
1. Subcellular Transcript Simulation:
   - 50,000 RNA molecules across 200 gene targets with (x, y, z) coordinates,
     and 1,500 simulated DAPI nuclear centroids in an 800x800 µm field of view.
2. Transcript Assignment to Cells:
   - KDTree-based Voronoi radius assignment (r <= 7.0 µm), cell segmentation filtering,
     and single-cell count matrix generation.
3. Cell-Type Annotation & 2D Spatial Graph Construction:
   - Marker-based cell type classification and radius spatial graph construction (r = 30 µm).
4. Cellular Neighborhood (CN) Niche Clustering:
   - Window aggregation of local cell-type compositions and K-means clustering
     into 4 distinct microenvironments: Cancer Core, Desmoplastic Rim, TLS Immune Hub, Bone Front.
5. Distance-Constrained Ligand-Receptor Signaling (COMMOT-Style):
   - Exponential distance decay S_ij = L_i * R_j * exp(-d_ij / 30) for:
     * TGFB1 -> TGFBR2 (CAF-to-Tumor EMT axis)
     * TNFSF11 (RANKL) -> TNFRSF11A (RANK) (Tumor-to-Bone osteolytic axis).
6. Export Structured Tables & Summaries to results_pipeline_02/.
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
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans


# ==============================================================================
# Step 1: Subcellular Transcript Simulation
# ==============================================================================
def simulate_xenium_subcellular_data(
    n_transcripts: int = 50000,
    n_cells: int = 1500,
    n_genes: int = 200,
    fov_size: float = 800.0,
    random_seed: int = 42,
) -> dict:
    """
    Simulates subcellular In Situ transcriptomics with exact (x, y) coordinates
    and nuclear centroids across 4 histological zones in an 800 µm x 800 µm tissue section.
    """
    print("=" * 80)
    print(f"[Step 1] Subcellular Transcript Simulation ({n_transcripts} Transcripts, {n_cells} Cells)")
    print("=" * 80)

    np.random.seed(random_seed)

    # 1. Simulate Cell Centroids & Ground Truth Cell Types
    cell_types_list = ["Tumor_AR_Pos", "CAF", "CD8_T", "Osteoblast", "Endothelial"]
    centroids_x = []
    centroids_y = []
    cell_types = []

    # Assign cell locations according to tissue histology
    # Zone 1: Cancer Core (center 300-500 um)
    n_tumor = int(n_cells * 0.45)
    t_x = np.random.normal(400.0, 75.0, size=n_tumor)
    t_y = np.random.normal(400.0, 75.0, size=n_tumor)
    centroids_x.extend(t_x)
    centroids_y.extend(t_y)
    cell_types.extend(["Tumor_AR_Pos"] * n_tumor)

    # Zone 2: Desmoplastic Rim (CAF ring around tumor core)
    n_caf = int(n_cells * 0.30)
    theta = np.random.uniform(0, 2 * np.pi, size=n_caf)
    r_caf = np.random.normal(210.0, 30.0, size=n_caf)
    c_x = 400.0 + r_caf * np.cos(theta)
    c_y = 400.0 + r_caf * np.sin(theta)
    centroids_x.extend(c_x)
    centroids_y.extend(c_y)
    cell_types.extend(["CAF"] * n_caf)

    # Zone 3: TLS Immune Hub (focal cluster at x=180, y=600)
    n_immune = int(n_cells * 0.12)
    i_x = np.random.normal(180.0, 35.0, size=n_immune)
    i_y = np.random.normal(600.0, 35.0, size=n_immune)
    centroids_x.extend(i_x)
    centroids_y.extend(i_y)
    cell_types.extend(["CD8_T"] * n_immune)

    # Zone 4: Bone Metastatic Front (top-right corner x=650, y=650)
    n_bone = int(n_cells * 0.08)
    b_x = np.random.normal(650.0, 45.0, size=n_bone)
    b_y = np.random.normal(650.0, 45.0, size=n_bone)
    centroids_x.extend(b_x)
    centroids_y.extend(b_y)
    cell_types.extend(["Osteoblast"] * n_bone)

    # Zone 5: Dispersed Endothelial cells
    n_endo = n_cells - (n_tumor + n_caf + n_immune + n_bone)
    e_x = np.random.uniform(50.0, fov_size - 50.0, size=n_endo)
    e_y = np.random.uniform(50.0, fov_size - 50.0, size=n_endo)
    centroids_x.extend(e_x)
    centroids_y.extend(e_y)
    cell_types.extend(["Endothelial"] * n_endo)

    centroids_x = np.clip(np.array(centroids_x), 10.0, fov_size - 10.0)
    centroids_y = np.clip(np.array(centroids_y), 10.0, fov_size - 10.0)
    cell_ids = [f"CELL_{idx:04d}" for idx in range(n_cells)]
    cell_types = np.array(cell_types)

    # 2. Gene Targets Definition (200 genes)
    gene_list = [f"Target_{i:03d}" for i in range(n_genes)]
    canonical_markers = {
        0: "KLK3", 1: "AR", 2: "FOLH1", 3: "NKX3-1", 4: "FOXA1",
        10: "TGFB1", 11: "TGFBR2", 12: "TNFSF11", 13: "TNFRSF11A",
        20: "ACTA2", 21: "COL1A1", 22: "COL1A2", 23: "POSTN", 24: "FAP",
        30: "CD3D", 31: "CD3E", 32: "CD8A", 33: "IFNG", 34: "GZMB",
        40: "BGLAP", 41: "RUNX2", 42: "ALPL", 43: "SPP1",
        50: "PECAM1", 51: "CD34", 52: "VWF",
    }
    for idx, name in canonical_markers.items():
        gene_list[idx] = name

    # 3. Simulate Subcellular Transcripts
    # ~85% inside cells (Gaussian spread around cell centroid sigma ~ 3.2 um)
    # ~15% extracellular ambient background
    n_in_cell = int(n_transcripts * 0.85)
    n_background = n_transcripts - n_in_cell

    # Assign transcripts to cells
    cell_assignments_true = np.random.choice(n_cells, size=n_in_cell)
    tx_x = centroids_x[cell_assignments_true] + np.random.normal(0, 3.2, size=n_in_cell)
    tx_y = centroids_y[cell_assignments_true] + np.random.normal(0, 3.2, size=n_in_cell)

    # Assign genes based on parent cell type
    tx_genes = []
    for c_idx in cell_assignments_true:
        ct = cell_types[c_idx]
        if ct == "Tumor_AR_Pos":
            p_genes = ["KLK3", "AR", "FOLH1", "NKX3-1", "TNFSF11", "TGFBR2"]
            weights = [0.25, 0.20, 0.15, 0.15, 0.12, 0.13]
        elif ct == "CAF":
            p_genes = ["ACTA2", "COL1A1", "COL1A2", "POSTN", "TGFB1"]
            weights = [0.30, 0.25, 0.20, 0.15, 0.10]
        elif ct == "CD8_T":
            p_genes = ["CD3D", "CD3E", "CD8A", "IFNG", "GZMB"]
            weights = [0.25, 0.25, 0.25, 0.15, 0.10]
        elif ct == "Osteoblast":
            p_genes = ["BGLAP", "RUNX2", "ALPL", "TNFRSF11A"]
            weights = [0.30, 0.30, 0.25, 0.15]
        else:  # Endothelial
            p_genes = ["PECAM1", "CD34", "VWF"]
            weights = [0.40, 0.35, 0.25]

        # 75% marker gene, 25% background gene
        if np.random.rand() < 0.75:
            w_norm = np.array(weights, dtype=np.float64)
            w_norm /= w_norm.sum()
            gene = np.random.choice(p_genes, p=w_norm)
        else:
            gene = np.random.choice(gene_list)
        tx_genes.append(gene)

    # Extracellular ambient transcripts
    bg_x = np.random.uniform(5.0, fov_size - 5.0, size=n_background)
    bg_y = np.random.uniform(5.0, fov_size - 5.0, size=n_background)
    bg_genes = np.random.choice(gene_list, size=n_background)

    all_tx_x = np.concatenate([tx_x, bg_x])
    all_tx_y = np.concatenate([tx_y, bg_y])
    all_tx_genes = tx_genes + list(bg_genes)
    qv_scores = np.random.uniform(22.0, 40.0, size=n_transcripts)

    print(f"  • Simulated FOV: {fov_size} x {fov_size} µm")
    print(f"  • Total transcripts: {n_transcripts} (In-cell: {n_in_cell}, Background: {n_background})")
    print(f"  • Cell centroids: {n_cells} across 5 cell types: {dict(pd.Series(cell_types).value_counts())}")

    return {
        "cell_ids": cell_ids,
        "centroids_x": centroids_x,
        "centroids_y": centroids_y,
        "true_cell_types": cell_types,
        "tx_x": all_tx_x,
        "tx_y": all_tx_y,
        "tx_genes": all_tx_genes,
        "qv_scores": qv_scores,
        "gene_list": gene_list,
    }


# ==============================================================================
# Step 2: Transcript Assignment to Cells (Voronoi / KDTree Radius Assignment)
# ==============================================================================
def assign_transcripts_to_cells(
    subcellular_data: dict,
    max_radius: float = 7.0,
    min_transcripts_per_cell: int = 15,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assigns each transcript molecule to the nearest cell centroid using a spatial
    KDTree within a cellular boundary threshold (r <= 7.0 µm).
    Filters out cells with fewer than 15 transcripts.
    """
    print("\n" + "=" * 80)
    print(f"[Step 2] Transcript Assignment to Cells (KDTree radius <= {max_radius} µm)")
    print("=" * 80)

    cell_coords = np.column_stack([subcellular_data["centroids_x"], subcellular_data["centroids_y"]])
    tx_coords = np.column_stack([subcellular_data["tx_x"], subcellular_data["tx_y"]])
    n_cells = len(cell_coords)
    n_transcripts = len(tx_coords)

    # Build KDTree on cell centroids
    tree = cKDTree(cell_coords)
    dists, indices = tree.query(tx_coords, k=1, distance_upper_bound=max_radius)

    # Assigned mask
    assigned_mask = dists <= max_radius
    assigned_indices = indices[assigned_mask]
    assigned_genes = np.array(subcellular_data["tx_genes"])[assigned_mask]

    print(f"  • Total transcripts evaluated: {n_transcripts}")
    print(f"  • Assigned to cells: {np.sum(assigned_mask)} ({np.sum(assigned_mask)/n_transcripts*100:.1f}%)")
    print(f"  • Extracellular / unassigned: {np.sum(~assigned_mask)} ({np.sum(~assigned_mask)/n_transcripts*100:.1f}%)")

    # Build count matrix (n_cells x n_genes)
    gene_to_idx = {g: i for i, g in enumerate(subcellular_data["gene_list"])}
    count_matrix = np.zeros((n_cells, len(subcellular_data["gene_list"])), dtype=np.int32)

    for c_idx, g_name in zip(assigned_indices, assigned_genes):
        count_matrix[c_idx, gene_to_idx[g_name]] += 1

    cell_counts = np.sum(count_matrix, axis=1)

    # Cell-level QC: min transcripts
    pass_cell_mask = cell_counts >= min_transcripts_per_cell
    print(f"  • Cells passing transcript filter (>= {min_transcripts_per_cell}): {np.sum(pass_cell_mask)} / {n_cells}")

    # Build cell metadata DataFrame
    df_cells = pd.DataFrame(
        {
            "cell_id": np.array(subcellular_data["cell_ids"])[pass_cell_mask],
            "x_um": subcellular_data["centroids_x"][pass_cell_mask],
            "y_um": subcellular_data["centroids_y"][pass_cell_mask],
            "true_cell_type": subcellular_data["true_cell_types"][pass_cell_mask],
            "total_transcripts": cell_counts[pass_cell_mask],
            "cell_area_um2": np.pi * (max_radius ** 2),  # circular assumption
        }
    )

    filtered_counts = count_matrix[pass_cell_mask]
    df_counts = pd.DataFrame(filtered_counts, columns=subcellular_data["gene_list"], index=df_cells["cell_id"])

    return df_cells, df_counts


# ==============================================================================
# Step 3: Cell-Type Annotation & 2D Spatial Graph Construction
# ==============================================================================
def annotate_cells_and_build_spatial_graph(
    df_cells: pd.DataFrame,
    df_counts: pd.DataFrame,
    radius_graph: float = 30.0,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Annotates cells based on marker profiles and constructs a 2D radius spatial graph
    (r = 30 µm) representing direct cell-cell contacts and paracrine signaling limits.
    """
    print("\n" + "=" * 80)
    print(f"[Step 3] Cell-Type Annotation & Spatial Graph Construction (r = {radius_graph} µm)")
    print("=" * 80)

    # Scoring markers
    tumor_score = df_counts["KLK3"] + df_counts["AR"] + df_counts["FOLH1"]
    caf_score = df_counts["ACTA2"] + df_counts["COL1A1"] + df_counts["POSTN"]
    cd8_score = df_counts["CD3D"] + df_counts["CD8A"] + df_counts["IFNG"]
    osteo_score = df_counts["BGLAP"] + df_counts["RUNX2"] + df_counts["ALPL"]
    endo_score = df_counts["PECAM1"] + df_counts["CD34"]

    score_matrix = np.column_stack([tumor_score, caf_score, cd8_score, osteo_score, endo_score])
    labels = ["Tumor_AR_Pos", "CAF", "CD8_T", "Osteoblast", "Endothelial"]
    assigned_types = [labels[idx] for idx in np.argmax(score_matrix, axis=1)]
    df_cells["annotated_cell_type"] = assigned_types

    # Agreement with ground truth
    accuracy = np.mean(df_cells["annotated_cell_type"] == df_cells["true_cell_type"]) * 100.0
    print(f"  • Cell type annotation accuracy vs ground truth: {accuracy:.2f}%")
    print(f"  • Annotated cell distribution: {dict(df_cells['annotated_cell_type'].value_counts())}")

    # Build 2D spatial graph (radius r = 30 µm)
    coords = df_cells[["x_um", "y_um"]].values
    tree = cKDTree(coords)
    adj_pairs = tree.query_pairs(r=radius_graph)

    n_cells = len(df_cells)
    adj_matrix = np.zeros((n_cells, n_cells), dtype=np.float32)
    for i, j in adj_pairs:
        adj_matrix[i, j] = 1.0
        adj_matrix[j, i] = 1.0

    degrees = np.sum(adj_matrix, axis=1)
    print(f"  • Spatial graph edges: {len(adj_pairs)} (Mean neighbors per cell: {np.mean(degrees):.1f})")

    return df_cells, adj_matrix


# ==============================================================================
# Step 4: Cellular Neighborhood (CN) Niche Clustering
# ==============================================================================
def cluster_cellular_neighborhoods(
    df_cells: pd.DataFrame,
    adj_matrix: np.ndarray,
    n_niches: int = 4,
    random_seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregates local cell-type composition vectors within the neighborhood window
    (r <= 30 µm) and identifies 4 spatial niches using K-means. Niches get DATA-DRIVEN labels
    ("Niche_<k> (<dominant cell type>-dominant)"); biological interpretations such as
    Cancer Core / Desmoplastic Rim / TLS Immune Hub / Bone Front are hypotheses only and are
    NOT written to the results (naming them from the dominant type would be circular).
    """
    print("\n" + "=" * 80)
    print(f"[Step 4] Cellular Neighborhood (CN) Niche Clustering (K={n_niches} Niches)")
    print("=" * 80)

    cell_types = ["Tumor_AR_Pos", "CAF", "CD8_T", "Osteoblast", "Endothelial"]
    n_cells = len(df_cells)
    n_types = len(cell_types)

    # One-hot encoded cell type matrix
    type_dummies = pd.get_dummies(df_cells["annotated_cell_type"])[cell_types].values.astype(np.float32)

    # Local neighborhood composition: (A + I) * OneHot / degree
    # Self included in neighborhood
    A_self = adj_matrix + np.eye(n_cells, dtype=np.float32)
    degrees_self = np.sum(A_self, axis=1, keepdims=True)
    local_compositions = (A_self @ type_dummies) / degrees_self

    # K-means clustering of composition vectors
    kmeans = KMeans(n_clusters=n_niches, random_state=random_seed, n_init=10)
    niche_labels = kmeans.fit_predict(local_compositions)

    # Interpret niches by their dominant cell type
    centers = kmeans.cluster_centers_
    # SCI-02-068: labels are data-driven (cluster id + dominant cell type of the K-means center). The simulated
    # cohort composition is unchanged; biological names (hypotheses, not findings) are intentionally not assigned.
    # The synthetic cohort has only 5 cell types (no B cells/DCs/osteoclasts, no tumor in the bone front), so these
    # niches do not reproduce the 4 niches described in the docs (ch07/ch08).
    niche_name_map = {}
    for cluster_id in range(n_niches):
        dom_ct = cell_types[int(np.argmax(centers[cluster_id]))]
        niche_name_map[cluster_id] = f"Niche_{cluster_id} ({dom_ct}-dominant)"

    # Resolve any duplicate name conflicts
    assigned_names = []
    for c_id in niche_labels:
        assigned_names.append(niche_name_map[c_id])
    df_cells["spatial_niche"] = assigned_names

    # Niche composition summary
    df_comp = pd.DataFrame(local_compositions, columns=[f"prop_{ct}" for ct in cell_types])
    df_comp["spatial_niche"] = assigned_names

    summary_df = df_comp.groupby("spatial_niche").mean()
    print("  Identified 4 Spatial Niches and their Average Cell Type Proportions:")
    for niche, row in summary_df.iterrows():
        print(f"    - [{niche}]: Tumor: {row['prop_Tumor_AR_Pos']*100:.1f}%, CAF: {row['prop_CAF']*100:.1f}%, CD8 T: {row['prop_CD8_T']*100:.1f}%, Osteo: {row['prop_Osteoblast']*100:.1f}%")

    return df_cells, summary_df


# ==============================================================================
# Step 5: Distance-Constrained Ligand-Receptor Signaling (COMMOT-Style)
# ==============================================================================
def compute_distance_decay_signaling(
    df_cells: pd.DataFrame,
    df_counts: pd.DataFrame,
    d_0: float = 30.0,
    max_dist: float = 60.0,
) -> pd.DataFrame:
    """
    Computes spatial distance-decayed ligand-receptor signaling flux:
        S_{ij} = L_i * R_j * exp(-d_{ij} / d_0)
    for two key axes in prostate cancer:
    1. TGFB1 -> TGFBR2 (CAF-to-Tumor EMT & desmoplasia)
    2. TNFSF11 (RANKL) -> TNFRSF11A (RANK) (Tumor-to-Bone osteoclast activation).
    """
    print("\n" + "=" * 80)
    print(f"[Step 5] Distance-Constrained L-R Signaling Flux (d_0 = {d_0} µm)")
    print("=" * 80)

    coords = df_cells[["x_um", "y_um"]].values
    n_cells = len(df_cells)
    tree = cKDTree(coords)

    # Pairs within max_dist
    pairs = tree.query_pairs(r=max_dist)
    print(f"  • Evaluated {len(pairs)} spatial cell pairs within {max_dist} µm")

    # Ligand-Receptor pairs to test
    lr_pairs = [
        ("TGFB1", "TGFBR2", "CAF_to_Tumor_EMT"),
        ("TNFSF11", "TNFRSF11A", "Tumor_to_Bone_Osteolysis"),
    ]

    signaling_records = []

    for lig_gene, rec_gene, biological_axis in lr_pairs:
        lig_expr = df_counts[lig_gene].values.astype(np.float32)
        rec_expr = df_counts[rec_gene].values.astype(np.float32)

        # Pairwise flux accumulator across cell type pairs
        ct_pair_fluxes = {}
        niche_pair_fluxes = {}

        for i, j in pairs:
            dist = np.sqrt(np.sum((coords[i] - coords[j]) ** 2))
            decay = np.exp(-dist / d_0)

            # Direction 1: i -> j
            flux_ij = lig_expr[i] * rec_expr[j] * decay
            if flux_ij > 0:
                ct_pair = (df_cells["annotated_cell_type"].iloc[i], df_cells["annotated_cell_type"].iloc[j])
                ct_pair_fluxes[ct_pair] = ct_pair_fluxes.get(ct_pair, 0.0) + flux_ij

                niche_pair = (df_cells["spatial_niche"].iloc[i], df_cells["spatial_niche"].iloc[j])
                niche_pair_fluxes[niche_pair] = niche_pair_fluxes.get(niche_pair, 0.0) + flux_ij

            # Direction 2: j -> i
            flux_ji = lig_expr[j] * rec_expr[i] * decay
            if flux_ji > 0:
                ct_pair_rev = (df_cells["annotated_cell_type"].iloc[j], df_cells["annotated_cell_type"].iloc[i])
                ct_pair_fluxes[ct_pair_rev] = ct_pair_fluxes.get(ct_pair_rev, 0.0) + flux_ji

                niche_pair_rev = (df_cells["spatial_niche"].iloc[j], df_cells["spatial_niche"].iloc[i])
                niche_pair_fluxes[niche_pair_rev] = niche_pair_fluxes.get(niche_pair_rev, 0.0) + flux_ji

        # Convert to records
        total_lr_flux = sum(ct_pair_fluxes.values()) + 1e-8
        for (sender_ct, rec_ct), fl in sorted(ct_pair_fluxes.items(), key=lambda x: x[1], reverse=True)[:5]:
            signaling_records.append(
                {
                    "biological_axis": biological_axis,
                    "ligand": lig_gene,
                    "receptor": rec_gene,
                    "sender_type": sender_ct,
                    "receiver_type": rec_ct,
                    "signaling_flux": float(fl),
                    "flux_percentage": float(fl / total_lr_flux * 100.0),
                }
            )

    df_signaling = pd.DataFrame(signaling_records)
    print("\n  Top Spatial Ligand-Receptor Signaling Fluxes:")
    print(df_signaling[["biological_axis", "sender_type", "receiver_type", "signaling_flux", "flux_percentage"]].to_string(index=False))

    return df_signaling


# ==============================================================================
# Step 6: Export Results to results_pipeline_02/
# ==============================================================================
def export_pipeline_02_results(
    df_cells: pd.DataFrame,
    df_counts: pd.DataFrame,
    df_niche_summary: pd.DataFrame,
    df_signaling: pd.DataFrame,
    out_dir: str,
):
    """
    Exports all subcellular, niche, and signaling tables to results_pipeline_02/.
    """
    print("\n" + "=" * 80)
    print(f"[Step 6] Exporting Results to {out_dir}/")
    print("=" * 80)

    os.makedirs(out_dir, exist_ok=True)

    # 1. Cell metadata
    df_cells.to_csv(os.path.join(out_dir, "xenium_cell_metadata.csv"), index=False)

    # 2. Key marker expression matrix
    key_markers = [
        "KLK3", "AR", "FOLH1", "NKX3-1", "TGFB1", "TGFBR2", "TNFSF11", "TNFRSF11A",
        "ACTA2", "COL1A1", "POSTN", "CD3D", "CD8A", "IFNG", "BGLAP", "RUNX2", "PECAM1",
    ]
    df_counts[key_markers].to_csv(os.path.join(out_dir, "xenium_cell_expression_matrix.csv"))

    # 3. Niche summary
    df_niche_summary.to_csv(os.path.join(out_dir, "niche_composition_summary.csv"))

    # 4. Signaling flux
    df_signaling.to_csv(os.path.join(out_dir, "distance_decay_signaling_summary.csv"), index=False)

    # 5. JSON Summary
    summary = {
        "pipeline": "02_xenium_subcellular_niche_and_signaling_pipeline",
        "total_cells_analyzed": int(len(df_cells)),
        "cell_types": {k: int(v) for k, v in df_cells["annotated_cell_type"].value_counts().items()},
        "niches": {k: int(v) for k, v in df_cells["spatial_niche"].value_counts().items()},
        "top_signaling_axes": [
            {
                "biological_axis": str(r["biological_axis"]),
                "sender_type": str(r["sender_type"]),
                "receiver_type": str(r["receiver_type"]),
                "flux_percentage": float(r["flux_percentage"]),
            }
            for _, r in df_signaling.head(4).iterrows()
        ],
        "limitations": [
            "Niche labels are data-driven (cluster id + dominant cell type); no biological niche names are assigned.",
            "Synthetic cohort has 5 cell types only (no osteoclasts, B cells or DCs); it does not reproduce the 4 niches in the docs.",
            "RANKL (TNFSF11->TNFRSF11A) axis top flux is Tumor->Tumor (79.4%) with no osteoclast type; it is not evidence of osteoclast activation.",
        ],
    }
    with open(os.path.join(out_dir, "pipeline_02_execution_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("[✓] Successfully exported all pipeline 02 tables and summaries:")
    print(f"    - {os.path.join(out_dir, 'xenium_cell_metadata.csv')}")
    print(f"    - {os.path.join(out_dir, 'xenium_cell_expression_matrix.csv')}")
    print(f"    - {os.path.join(out_dir, 'niche_composition_summary.csv')}")
    print(f"    - {os.path.join(out_dir, 'distance_decay_signaling_summary.csv')}")
    print(f"    - {os.path.join(out_dir, 'pipeline_02_execution_summary.json')}")


# ==============================================================================
# Main Pipeline Driver
# ==============================================================================
def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "pipelines", "results_pipeline_02")

    # Step 1: Simulate Subcellular Transcripts
    subcellular_data = simulate_xenium_subcellular_data(
        n_transcripts=50000,
        n_cells=1500,
        n_genes=200,
        fov_size=800.0,
        random_seed=42,
    )

    # Step 2: Assign Transcripts to Cells
    df_cells, df_counts = assign_transcripts_to_cells(
        subcellular_data,
        max_radius=7.0,
        min_transcripts_per_cell=15,
    )

    # Step 3: Annotate Cells & Build Spatial Graph
    df_cells, adj_matrix = annotate_cells_and_build_spatial_graph(
        df_cells,
        df_counts,
        radius_graph=30.0,
    )

    # Step 4: Cellular Neighborhood (CN) Niche Clustering
    df_cells, df_niche_summary = cluster_cellular_neighborhoods(
        df_cells,
        adj_matrix,
        n_niches=4,
        random_seed=42,
    )

    # Step 5: Distance-Constrained Signaling
    df_signaling = compute_distance_decay_signaling(
        df_cells,
        df_counts,
        d_0=30.0,
        max_dist=60.0,
    )

    # Step 6: Export Results
    export_pipeline_02_results(
        df_cells,
        df_counts,
        df_niche_summary,
        df_signaling,
        out_dir,
    )

    print("\n" + "=" * 80)
    print("[✓] Pipeline 02 completed successfully with Exit Code 0!")
    print("=" * 80)


if __name__ == "__main__":
    main()

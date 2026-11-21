#!/usr/bin/env python3
"""
tracks/03_tme_cell_communication/pipelines/02_nichenet_and_spatial_communication_pipeline.py
=============================================================================================
End-to-End Pipeline for NicheNet Intracellular Ligand-to-Target Causal Attribution,
Random Walk with Restart (RWR) Regulatory Potential Matrix Calculation, AUROC / Pearson
Scoring, and Track 02 Spatial Distance-Decay Signaling Flux in Prostate Bone Metastasis.
"""

import os
import sys
import math
import random
import csv
import json
import time
from collections import defaultdict


# ==============================================================================
# Step 1: Simulate Receiver Cell DEGs & Biological Programs
# ==============================================================================
def simulate_receiver_deg_programs():
    """
    Defines realistic transcriptional programs triggered by microenvironmental signaling
    in mCRPC bone metastasis:
      1. CD8+ T-cell Exhaustion Program
      2. Tumor EMT & Invasion Program
      3. Osteoclast Maturation & Resorption Program
      4. Control background non-DEGs (50 genes)
    """
    deg_programs = {
        "CD8_T_Exhaustion": [
            "PDCD1", "HAVCR2", "LAG3", "TOX", "TIGIT", "CTLA4", "ENTPD1", "PRDM1",
            "CXCL13", "LAYN", "BATF", "IRF4"
        ],
        "Tumor_EMT_Metastasis": [
            "VIM", "CDH2", "SNAI1", "SNAI2", "ZEB1", "TWIST1", "FN1", "MMP2", "MMP9",
            "ACTA2", "COL1A1", "TGFBI"
        ],
        "Osteoclast_Maturation": [
            "NFATC1", "ACP5", "CTSK", "CALCR", "ITGB3", "OSCAR", "DCSTAMP", "OCSTAMP",
            "ATP6V0D2", "TRAF6"
        ],
    }

    # 50 Control non-DEG background genes
    control_genes = [f"Ctrl_Gene_{i:02d}" for i in range(1, 51)]

    # Flatten unique target genes
    all_target_genes = []
    for prog, genes in deg_programs.items():
        for g in genes:
            if g not in all_target_genes:
                all_target_genes.append(g)
    for cg in control_genes:
        all_target_genes.append(cg)

    return deg_programs, control_genes, all_target_genes


# ==============================================================================
# Step 2: Prior Signaling-GRN Graph & RWR Regulatory Potential Matrix W
# ==============================================================================
def construct_prior_signaling_grn_graph(all_target_genes):
    """
    Builds a 3-layer prior knowledge network linking:
      Layer 1: Upstream Ligands
      Layer 2: Receptors
      Layer 3: Signaling Intracellular Kinases & Transducers
      Layer 4: Transcription Factors (TFs, Track 01 GRN link)
      Layer 5: Downstream Target Genes (DEGs + Controls)
    """
    ligands = [
        "TGFB1", "CXCL12", "TNFSF11", "IL6", "SPP1", "VEGFA", "PTHLH", "WNT5A", "BMP2"
    ]
    receptors = [
        "TGFBR1", "TGFBR2", "CXCR4", "TNFRSF11A", "IL6R", "IL6ST",
        "CD44", "KDR", "PTH1R", "FZD5", "BMPR1A"
    ]
    transducers = [
        "SMAD2", "SMAD3", "SMAD4", "JAK1", "STAT3", "TRAF6",
        "MAPK1", "MAPK3", "PIK3CA", "AKT1", "SRC", "CTNNB1", "SMAD1", "SMAD5"
    ]
    tfs = [
        "AR", "FOXA1", "NFATC1", "NFKB1", "STAT3_TF", "SMAD3_TF",
        "SNAI1_TF", "ZEB1_TF", "TOX_TF", "RUNX2_TF"
    ]

    all_nodes = list(dict.fromkeys(ligands + receptors + transducers + tfs + all_target_genes))
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}

    # Graph adjacency: source -> list of (target, weight)
    adj = defaultdict(list)

    # 1. Ligand -> Receptor edges
    lr_edges = [
        ("TGFB1", "TGFBR1", 1.0), ("TGFB1", "TGFBR2", 1.0),
        ("CXCL12", "CXCR4", 1.0),
        ("TNFSF11", "TNFRSF11A", 1.0),
        ("IL6", "IL6R", 1.0), ("IL6", "IL6ST", 0.9),
        ("SPP1", "CD44", 0.95),
        ("VEGFA", "KDR", 0.95),
        ("PTHLH", "PTH1R", 1.0),
        ("WNT5A", "FZD5", 0.9),
        ("BMP2", "BMPR1A", 0.95),
    ]
    for u, v, w in lr_edges:
        adj[u].append((v, w))

    # 2. Receptor -> Transducers
    rt_edges = [
        ("TGFBR1", "SMAD2", 0.95), ("TGFBR1", "SMAD3", 0.95), ("TGFBR2", "TGFBR1", 0.9),
        ("CXCR4", "PIK3CA", 0.85), ("CXCR4", "MAPK1", 0.8), ("CXCR4", "SRC", 0.75),
        ("TNFRSF11A", "TRAF6", 0.95), ("TRAF6", "MAPK3", 0.8),
        ("IL6R", "JAK1", 0.95), ("IL6ST", "JAK1", 0.95), ("JAK1", "STAT3", 0.95),
        ("CD44", "SRC", 0.8), ("CD44", "AKT1", 0.85),
        ("KDR", "PIK3CA", 0.85), ("KDR", "MAPK1", 0.8),
        ("PTH1R", "SRC", 0.8), ("PTH1R", "AKT1", 0.8),
        ("FZD5", "CTNNB1", 0.85),
        ("BMPR1A", "SMAD1", 0.9), ("BMPR1A", "SMAD5", 0.9),
    ]
    for u, v, w in rt_edges:
        adj[u].append((v, w))

    # 3. Transducers -> TFs
    tf_edges = [
        ("SMAD2", "SMAD4", 0.9), ("SMAD3", "SMAD4", 0.9),
        ("SMAD4", "SMAD3_TF", 0.95), ("SMAD3_TF", "SNAI1_TF", 0.85), ("SMAD3_TF", "ZEB1_TF", 0.85),
        ("TRAF6", "NFKB1", 0.9), ("TRAF6", "NFATC1", 0.95),
        ("STAT3", "STAT3_TF", 0.95), ("STAT3_TF", "AR", 0.75),
        ("MAPK1", "FOXA1", 0.7), ("MAPK3", "NFATC1", 0.8),
        ("AKT1", "AR", 0.75), ("AKT1", "ZEB1_TF", 0.8),
        ("CTNNB1", "AR", 0.8), ("CTNNB1", "RUNX2_TF", 0.85),
        ("SMAD1", "RUNX2_TF", 0.85), ("SMAD5", "RUNX2_TF", 0.85),
        ("MAPK1", "TOX_TF", 0.85), ("AKT1", "TOX_TF", 0.75),
    ]
    for u, v, w in tf_edges:
        adj[u].append((v, w))

    # 4. TFs -> Target Genes (DEGs)
    # CD8 T Exhaustion targets
    t_exh_edges = [
        ("TOX_TF", "TOX", 0.95), ("TOX_TF", "PDCD1", 0.9), ("TOX_TF", "HAVCR2", 0.9),
        ("TOX_TF", "LAG3", 0.85), ("TOX_TF", "TIGIT", 0.85), ("TOX_TF", "CTLA4", 0.8),
        ("TOX_TF", "ENTPD1", 0.8), ("TOX_TF", "PRDM1", 0.8), ("TOX_TF", "CXCL13", 0.75),
        ("TOX_TF", "LAYN", 0.75), ("TOX_TF", "BATF", 0.7), ("TOX_TF", "IRF4", 0.7),
        ("NFKB1", "PDCD1", 0.7), ("STAT3_TF", "PRDM1", 0.75),
    ]
    for u, v, w in t_exh_edges:
        adj[u].append((v, w))

    # Tumor EMT targets
    emt_edges = [
        ("SNAI1_TF", "VIM", 0.95), ("SNAI1_TF", "CDH2", 0.9), ("SNAI1_TF", "SNAI1", 0.95),
        ("SNAI1_TF", "SNAI2", 0.85), ("SNAI1_TF", "FN1", 0.9), ("SNAI1_TF", "MMP2", 0.85),
        ("ZEB1_TF", "ZEB1", 0.95), ("ZEB1_TF", "TWIST1", 0.85), ("ZEB1_TF", "MMP9", 0.85),
        ("ZEB1_TF", "ACTA2", 0.8), ("ZEB1_TF", "COL1A1", 0.85), ("ZEB1_TF", "TGFBI", 0.9),
        ("SMAD3_TF", "VIM", 0.85), ("SMAD3_TF", "TGFBI", 0.9),
    ]
    for u, v, w in emt_edges:
        adj[u].append((v, w))

    # Osteoclast Maturation targets
    osteo_edges = [
        ("NFATC1", "NFATC1", 0.95), ("NFATC1", "ACP5", 0.95), ("NFATC1", "CTSK", 0.95),
        ("NFATC1", "CALCR", 0.9), ("NFATC1", "ITGB3", 0.85), ("NFATC1", "OSCAR", 0.9),
        ("NFATC1", "DCSTAMP", 0.9), ("NFATC1", "OCSTAMP", 0.85), ("NFATC1", "ATP6V0D2", 0.85),
        ("NFKB1", "TRAF6", 0.85), ("NFKB1", "NFATC1", 0.9),
    ]
    for u, v, w in osteo_edges:
        adj[u].append((v, w))

    # Weak baseline links to control genes (to simulate biological noise)
    random.seed(42)
    for cg in all_target_genes:
        if cg.startswith("Ctrl_Gene"):
            # Random weak edge from a random TF
            rand_tf = random.choice(tfs)
            adj[rand_tf].append((cg, 0.05 + 0.05 * random.random()))

    return ligands, receptors, transducers, tfs, all_nodes, node_to_idx, adj


def compute_rwr_regulatory_potential_matrix(ligands, all_target_genes, all_nodes, node_to_idx, adj, restart_r=0.5):
    """
    Computes Random Walk with Restart (RWR) stationary probability vector for each ligand:
      p^(t+1) = (1 - r) * M * p^(t) + r * e_L
    Extracts Regulatory Potential Matrix W [L x G].
    """
    n_nodes = len(all_nodes)

    # Build transition probability matrix: P(u -> v)
    transition_weights = defaultdict(dict)
    for u in all_nodes:
        out_edges = adj.get(u, [])
        total_w = sum(w for _, w in out_edges)
        if total_w > 0:
            for v, w in out_edges:
                transition_weights[u][v] = w / total_w
        else:
            # Self-loop for absorbing sink nodes
            transition_weights[u][u] = 1.0

    W_matrix = {}  # ligand -> dict(target_gene -> score)

    for lig in ligands:
        lig_idx = node_to_idx[lig]

        # Initial probability vector p
        p = [0.0] * n_nodes
        p[lig_idx] = 1.0

        # Iterative RWR Power Iteration
        for _ in range(60):
            p_next = [0.0] * n_nodes
            # Flow propagation: p_next[v] = sum_u (p[u] * P(u -> v))
            for u_name, u_idx in node_to_idx.items():
                prob_u = p[u_idx]
                if prob_u > 0:
                    for v_name, prob_uv in transition_weights[u_name].items():
                        v_idx = node_to_idx[v_name]
                        p_next[v_idx] += prob_u * prob_uv

            # Apply restart probability
            for i in range(n_nodes):
                p_next[i] = (1.0 - restart_r) * p_next[i]
            p_next[lig_idx] += restart_r

            # Check convergence
            l1_diff = sum(abs(p_next[i] - p[i]) for i in range(n_nodes))
            p = p_next
            if l1_diff < 1e-8:
                break

        # Extract target gene scores for ligand
        scores = {}
        for tg in all_target_genes:
            tg_idx = node_to_idx[tg]
            scores[tg] = p[tg_idx]

        W_matrix[lig] = scores

    return W_matrix


# ==============================================================================
# Step 3: Ligand Activity Scoring (Pearson Correlation & AUROC)
# ==============================================================================
def compute_pearson_correlation(x_vals, y_vals):
    n = len(x_vals)
    if n < 2:
        return 0.0
    mean_x = sum(x_vals) / n
    mean_y = sum(y_vals) / n

    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_vals, y_vals))
    var_x = sum((x - mean_x) ** 2 for x in x_vals)
    var_y = sum((y - mean_y) ** 2 for y in y_vals)

    denom = math.sqrt(var_x * var_y)
    if denom == 0:
        return 0.0
    return cov / denom


def compute_auroc(score_label_pairs):
    """
    Computes AUROC via Mann-Whitney U rank statistic.
    score_label_pairs: list of (score, binary_label)
    """
    pos_count = sum(1 for _, l in score_label_pairs if l == 1)
    neg_count = sum(1 for _, l in score_label_pairs if l == 0)
    if pos_count == 0 or neg_count == 0:
        return 0.5

    # Sort ascending by score with tie-breaking rank
    sorted_pairs = sorted(score_label_pairs, key=lambda x: x[0])

    # Assign average ranks for ties
    ranks = [0.0] * len(sorted_pairs)
    i = 0
    while i < len(sorted_pairs):
        j = i
        while j < len(sorted_pairs) and sorted_pairs[j][0] == sorted_pairs[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    # Sum of ranks for positive class
    rank_sum_pos = sum(ranks[k] for k, (_, label) in enumerate(sorted_pairs) if label == 1)
    u_stat = rank_sum_pos - (pos_count * (pos_count + 1)) / 2.0
    auroc = u_stat / (pos_count * neg_count)
    return round(auroc, 5)


def score_ligand_activity(W_matrix, ligands, deg_programs, all_target_genes):
    """
    Evaluates how accurately each ligand's regulatory potential predicts
    the observed DEGs in each biological program using Pearson and AUROC.
    """
    activity_results = []

    for prog_name, deg_genes in deg_programs.items():
        deg_set = set(deg_genes)

        for lig in ligands:
            scores = W_matrix[lig]
            x_vals = []
            y_vals = []
            pairs = []

            for g in all_target_genes:
                s = scores.get(g, 0.0)
                label = 1 if g in deg_set else 0
                x_vals.append(s)
                y_vals.append(label)
                pairs.append((s, label))

            pearson_r = compute_pearson_correlation(x_vals, y_vals)
            auroc_val = compute_auroc(pairs)

            activity_results.append({
                "program": prog_name,
                "ligand": lig,
                "n_degs": len(deg_genes),
                "pearson_r": round(pearson_r, 4),
                "auroc": auroc_val,
                "regulatory_potential_sum": round(sum(scores[g] for g in deg_genes), 6),
            })

    # Sort by program and AUROC descending
    activity_results.sort(key=lambda x: (x["program"], -x["auroc"]))
    return activity_results


# ==============================================================================
# Step 4: Spatial Distance Decay Signaling Flux (Track 02 Integration)
# ==============================================================================
def simulate_spatial_tme_coordinates(n_cells=1500, random_seed=42):
    """
    Simulates a 1,000 µm x 1,000 µm histological bone metastasis tissue core.
    Cell niches are spatially clustered:
      - Bone-tumor interface: Tumor + Osteoblasts + Osteoclasts
      - Stromal desmoplastic rim: myCAF
      - Vascular core: Endothelial
      - Peritumoral invasive margin: TAM M2 + CD8 T cells
    """
    random.seed(random_seed)
    coordinates = []

    cell_types = [
        ("Luminal_AR_Pos", 350, (500, 500), 180),
        ("t_NEPC", 150, (650, 450), 120),
        ("myCAF", 250, (400, 400), 220),
        ("Endothelial", 150, (520, 550), 160),
        ("TAM_M2", 250, (380, 580), 200),
        ("CD8_T", 150, (320, 620), 150),
        ("Osteoclast", 100, (750, 300), 100),
        ("Osteoblast", 100, (700, 350), 110),
    ]

    cell_id = 0
    for ctype, count, center, spread in cell_types:
        cx, cy = center
        for _ in range(count):
            x = max(10.0, min(990.0, random.gauss(cx, spread)))
            y = max(10.0, min(990.0, random.gauss(cy, spread)))
            coordinates.append({
                "cell_id": f"Cell_{cell_id:04d}",
                "cell_type": ctype,
                "x_um": round(x, 2),
                "y_um": round(y, 2),
            })
            cell_id += 1

    return coordinates


def compute_spatial_decay_flux(coordinates, bandwidth_sigma=30.0):
    """
    Applies spatial distance decay kernel w(d) = exp(-d / sigma)
    to calculate directional signaling flux between cell types.
    """
    cells_by_type = defaultdict(list)
    for c in coordinates:
        cells_by_type[c["cell_type"]].append((c["x_um"], c["y_um"]))

    cell_types = sorted(cells_by_type.keys())

    # Pre-defined expression strengths for canonical key pairs
    pair_expression = {
        ("Luminal_AR_Pos", "Osteoblast", "PTHLH_PTH1R"): (2.5, 3.4),
        ("Osteoblast", "Osteoclast", "TNFSF11_TNFRSF11A"): (3.7, 3.6),
        ("Osteoclast", "Luminal_AR_Pos", "TGFB1_TGFBR"): (3.4, 2.2),
        ("myCAF", "CD8_T", "CXCL12_CXCR4"): (3.8, 3.1),
        ("myCAF", "CD8_T", "TGFB1_TGFBR"): (3.6, 2.1),
        ("TAM_M2", "CD8_T", "CD274_PDCD1"): (3.1, 3.2),
        ("Luminal_AR_Pos", "TAM_M2", "CD47_SIRPA"): (2.8, 3.2),
        ("myCAF", "Endothelial", "VEGFA_KDR"): (2.5, 3.0),
    }

    spatial_flux_results = []

    for (sender, receiver, pathway), (e_lig, e_rec) in pair_expression.items():
        s_coords = cells_by_type[sender]
        r_coords = cells_by_type[receiver]

        # Sample cell pairs for spatial flux calculation
        sample_size = min(60, len(s_coords))
        s_sample = s_coords[:sample_size]
        r_sample = r_coords[:sample_size]

        total_weight = 0.0
        n_pairs = 0
        total_dist = 0.0

        for sx, sy in s_sample:
            for rx, ry in r_sample:
                dist = math.sqrt((sx - rx) ** 2 + (sy - ry) ** 2)
                total_dist += dist
                # Exponential decay kernel
                weight = math.exp(-dist / bandwidth_sigma)
                total_weight += weight
                n_pairs += 1

        avg_dist = total_dist / n_pairs if n_pairs > 0 else 0.0
        avg_weight = total_weight / n_pairs if n_pairs > 0 else 0.0

        raw_flux = e_lig * e_rec
        spatial_flux = raw_flux * avg_weight

        spatial_flux_results.append({
            "sender": sender,
            "receiver": receiver,
            "pathway": pathway,
            "mean_spatial_distance_um": round(avg_dist, 2),
            "distance_decay_weight": round(avg_weight, 5),
            "raw_expression_flux": round(raw_flux, 4),
            "effective_spatial_flux": round(spatial_flux, 5),
        })

    spatial_flux_results.sort(key=lambda x: x["effective_spatial_flux"], reverse=True)
    return spatial_flux_results


# ==============================================================================
# Step 5: Prostate Bone Metastasis Vicious Cycle Quantification & Drug Screening
# ==============================================================================
def quantify_vicious_cycle_and_inhibition(spatial_flux_results):
    """
    Quantifies the 4-leg vicious cycle flux in mCRPC bone metastasis:
      Leg 1: Tumor (PTHrP) -> Osteoblast
      Leg 2: Osteoblast (RANKL) -> Osteoclast
      Leg 3: Osteoclast (TGFb) -> Tumor
      Leg 4: myCAF (CXCL12 / TGFb) -> CD8 T (Immunosuppression)
    Simulates therapeutic blockade:
      - Baseline Untreated
      - Denosumab (anti-RANKL, -95% Leg 2)
      - Galunisertib (anti-TGFb, -90% Leg 3 & Leg 4)
      - Pembrolizumab (anti-PD1, -85% immune exhaustion)
      - Tri-combination therapy
    """
    # Lookup effective spatial flux
    flux_dict = {
        (r["sender"], r["receiver"], r["pathway"]): r["effective_spatial_flux"]
        for r in spatial_flux_results
    }

    leg1 = flux_dict.get(("Luminal_AR_Pos", "Osteoblast", "PTHLH_PTH1R"), 0.02)
    leg2 = flux_dict.get(("Osteoblast", "Osteoclast", "TNFSF11_TNFRSF11A"), 0.08)
    leg3 = flux_dict.get(("Osteoclast", "Luminal_AR_Pos", "TGFB1_TGFBR"), 0.03)
    leg4 = flux_dict.get(("myCAF", "CD8_T", "CXCL12_CXCR4"), 0.04)

    regimens = [
        ("Untreated_Baseline", 1.0, 1.0, 1.0, 1.0),
        ("Denosumab_anti_RANKL", 1.0, 0.05, 0.30, 1.0),
        ("Galunisertib_anti_TGFb", 0.50, 0.80, 0.10, 0.20),
        ("Pembrolizumab_anti_PD1", 0.90, 1.0, 1.0, 0.50),
        ("Combination_Denosumab_Galunisertib", 0.40, 0.05, 0.05, 0.15),
        ("Tri_Therapy_Denosumab_Galunisertib_antiPD1", 0.20, 0.02, 0.02, 0.05),
    ]

    vicious_cycle_metrics = []

    for reg_name, f1, f2, f3, f4 in regimens:
        eff_l1 = leg1 * f1
        eff_l2 = leg2 * f2
        eff_l3 = leg3 * f3
        eff_l4 = leg4 * f4

        # Bone osteolysis vicious cycle score (geometric mean of 3 legs)
        cycle_osteolysis_score = math.pow(eff_l1 * eff_l2 * eff_l3, 1.0 / 3.0)
        # Immunosuppression index
        immunosuppression_index = eff_l4 * 100.0
        # Overall vicious pathology score
        pathology_score = cycle_osteolysis_score * (1.0 + eff_l4 * 10.0)

        # Baseline comparison
        vicious_cycle_metrics.append({
            "treatment_regimen": reg_name,
            "flux_leg1_tumor_to_osteoblast": round(eff_l1, 6),
            "flux_leg2_osteoblast_to_osteoclast": round(eff_l2, 6),
            "flux_leg3_osteoclast_to_tumor": round(eff_l3, 6),
            "flux_leg4_mycaf_to_cd8t": round(eff_l4, 6),
            "vicious_cycle_osteolysis_score": round(cycle_osteolysis_score, 6),
            "immunosuppression_index": round(immunosuppression_index, 4),
            "composite_pathology_score": round(pathology_score, 6),
        })

    # Add percent reduction relative to baseline
    baseline_score = vicious_cycle_metrics[0]["composite_pathology_score"]
    for m in vicious_cycle_metrics:
        pct_inhibition = (1.0 - m["composite_pathology_score"] / baseline_score) * 100.0 if baseline_score > 0 else 0.0
        m["vicious_cycle_percent_inhibition"] = round(pct_inhibition, 2)

    return vicious_cycle_metrics


# ==============================================================================
# Main Pipeline Runner
# ==============================================================================
def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "pipelines", "results_pipeline_02")
    os.makedirs(outdir, exist_ok=True)

    print("================================================================================")
    print("Track 03 - Pipeline 02: NicheNet Causal Attribution & Spatial Signaling Flux")
    print("================================================================================")

    # 1. Receiver DEGs
    print("[1/5] Simulating receiver cell transcriptional DEG programs (Exhaustion, EMT, Osteoclast)...")
    deg_programs, control_genes, all_target_genes = simulate_receiver_deg_programs()
    deg_path = os.path.join(outdir, "receiver_deg_signatures.csv")
    with open(deg_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["program_name", "target_gene", "is_deg"])
        for prog, genes in deg_programs.items():
            for g in genes:
                writer.writerow([prog, g, 1])
        for cg in control_genes:
            writer.writerow(["Control_NonDEG", cg, 0])
    total_degs = sum(len(v) for v in deg_programs.values())
    print(f"      [✓] Exported DEG signatures: {deg_path} ({total_degs} DEGs across 3 programs)")

    # 2. Prior Graph & RWR
    print("[2/5] Building prior signaling-GRN graph and computing RWR regulatory potential matrix W...")
    ligands, receptors, transducers, tfs, all_nodes, node_to_idx, adj = construct_prior_signaling_grn_graph(all_target_genes)
    W_matrix = compute_rwr_regulatory_potential_matrix(ligands, all_target_genes, all_nodes, node_to_idx, adj, restart_r=0.5)

    matrix_path = os.path.join(outdir, "rwr_regulatory_potential_matrix_W.csv")
    with open(matrix_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand"] + all_target_genes)
        for lig in ligands:
            row = [lig] + [round(W_matrix[lig][tg], 7) for tg in all_target_genes]
            writer.writerow(row)
    print(f"      [✓] Exported Regulatory Potential Matrix W: {matrix_path} ({len(ligands)} ligands x {len(all_target_genes)} genes)")

    # 3. Ligand Activity Scoring
    print("[3/5] Scoring ligand activity predicting receiver DEGs via Pearson correlation and AUROC...")
    activity_results = score_ligand_activity(W_matrix, ligands, deg_programs, all_target_genes)
    act_path = os.path.join(outdir, "ligand_activity_auroc_pearson_ranking.csv")
    with open(act_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "program", "ligand", "n_degs", "pearson_r", "auroc", "regulatory_potential_sum"
        ])
        writer.writeheader()
        writer.writerows(activity_results)
    print(f"      [✓] Exported ligand activity rankings: {act_path} ({len(activity_results)} evaluations)")

    # 4. Spatial Coordinates & Distance Decay
    print("[4/5] Simulating histological bone metastasis coordinates and computing distance-decay flux (sigma=30µm)...")
    coordinates = simulate_spatial_tme_coordinates(n_cells=1500, random_seed=42)
    spatial_flux_results = compute_spatial_decay_flux(coordinates, bandwidth_sigma=30.0)
    spatial_path = os.path.join(outdir, "spatial_distance_decay_flux_matrix.csv")
    with open(spatial_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sender", "receiver", "pathway", "mean_spatial_distance_um",
            "distance_decay_weight", "raw_expression_flux", "effective_spatial_flux"
        ])
        writer.writeheader()
        writer.writerows(spatial_flux_results)
    print(f"      [✓] Exported spatial distance-decay flux: {spatial_path}")

    # 5. Vicious Cycle Quantification & Intervention Screening
    print("[5/5] Quantifying mCRPC bone metastasis vicious cycle flux & simulating targeted therapies...")
    vicious_cycle_metrics = quantify_vicious_cycle_and_inhibition(spatial_flux_results)
    cycle_path = os.path.join(outdir, "vicious_cycle_pathway_flux_and_inhibition_scores.csv")
    with open(cycle_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "treatment_regimen", "flux_leg1_tumor_to_osteoblast", "flux_leg2_osteoblast_to_osteoclast",
            "flux_leg3_osteoclast_to_tumor", "flux_leg4_mycaf_to_cd8t", "vicious_cycle_osteolysis_score",
            "immunosuppression_index", "composite_pathology_score", "vicious_cycle_percent_inhibition"
        ])
        writer.writeheader()
        writer.writerows(vicious_cycle_metrics)
    print(f"      [✓] Exported vicious cycle therapeutic screening: {cycle_path}")

    # Summary Metrics JSON
    summary_metrics = {
        "n_ligands_evaluated": len(ligands),
        "n_target_genes": len(all_target_genes),
        "n_deg_programs": len(deg_programs),
        "rwr_restart_probability": 0.5,
        "spatial_bandwidth_sigma_um": 30.0,
        "top_drivers_by_program": {
            prog: [
                f"{r['ligand']} (AUROC={r['auroc']}, r={r['pearson_r']})"
                for r in activity_results if r["program"] == prog
            ][:3]
            for prog in deg_programs
        },
        "baseline_vicious_cycle_pathology_score": vicious_cycle_metrics[0]["composite_pathology_score"],
        "tri_therapy_percent_inhibition": vicious_cycle_metrics[-1]["vicious_cycle_percent_inhibition"],
    }
    json_path = os.path.join(outdir, "pipeline_summary_metrics.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary_metrics, f, indent=2)
    print(f"      [✓] Exported pipeline summary: {json_path}")

    print("\n[✓] Pipeline 02 successfully completed all steps with exit code 0.")


if __name__ == "__main__":
    main()

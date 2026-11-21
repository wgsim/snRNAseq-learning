#!/usr/bin/env python3
"""
tracks/03_tme_cell_communication/pipelines/01_liana_and_cellchat_consensus_pipeline.py
======================================================================================
End-to-End Pipeline for Multi-Subunit Ligand-Receptor Scoring, LIANA Consensus
Aggregation, Cell-Label Permutation Testing, and 2D Network Centrality Analysis.

Biological Model:
  Castration-Resistant Prostate Cancer (mCRPC) Bone Metastasis Ecosystem.
  8 Cell Types: Luminal_AR_Pos, t_NEPC, myCAF, Endothelial, TAM_M2, CD8_T, Osteoclast, Osteoblast.
"""

import os
import sys
import math
import random
import csv
import json
import time
from collections import defaultdict, deque


# ==============================================================================
# Step 1: Cohort Simulation (1,500 cells, 8 cell types, 300 genes)
# ==============================================================================
def simulate_mcrpc_cohort(n_cells=1500, random_seed=42):
    random.seed(random_seed)

    cell_type_distribution = [
        ("Luminal_AR_Pos", 350),
        ("t_NEPC", 150),
        ("myCAF", 250),
        ("Endothelial", 150),
        ("TAM_M2", 250),
        ("CD8_T", 150),
        ("Osteoclast", 100),
        ("Osteoblast", 100),
    ]

    # Validate cell count
    total_cells = sum(cnt for _, cnt in cell_type_distribution)
    assert total_cells == n_cells, f"Total cells {total_cells} != {n_cells}"

    # Build gene list: 300 genes
    # 1) Key L-R genes & markers
    lr_genes = [
        "IL6", "IL6R", "IL6ST",
        "TGFB1", "TGFBR1", "TGFBR2",
        "TNFSF11", "TNFRSF11A",
        "CD274", "PDCD1",
        "CD47", "SIRPA",
        "CXCL12", "CXCR4",
        "SPP1", "CD44",
        "VEGFA", "FLT1", "KDR",
        "PTHLH", "PTH1R",
        "COL1A1", "ITGA1", "ITGB1",
        "DLL1", "NOTCH1",
        "WNT5A", "FZD5", "LRP5",
        "BMP2", "BMPR1A",
    ]
    marker_genes = [
        "KLK3", "AR", "FOXA1", "TMPRSS2",
        "SYP", "CHGA", "NCAM1", "ASCL1", "NEUROD1",
        "ACTA2", "FAP", "POSTN", "COL1A2",
        "PECAM1", "VWF", "CDH5",
        "CD68", "MRC1", "CD163", "ARG1",
        "CD8A", "CD3D", "GZMB", "PRF1", "HAVCR2", "LAG3", "TOX",
        "ACP5", "CTSK", "NFATC1", "CALCR",
        "BGLAP", "RUNX2", "ALPL", "COL1A1_ob",
    ]
    n_background = 300 - len(lr_genes) - len(marker_genes)
    background_genes = [f"Gene_{i:03d}" for i in range(1, n_background + 1)]
    all_genes = lr_genes + marker_genes + background_genes
    assert len(all_genes) == 300, f"Total genes: {len(all_genes)} != 300"

    # Define expression profiles per cell type (mean log-expression)
    # Baseline for all genes is low: ~0.1 - 0.4
    cell_metadata = []
    expression_matrix = []  # list of cell vectors

    # Signature boosts
    cell_signatures = {
        "Luminal_AR_Pos": {
            "KLK3": 3.8, "AR": 3.5, "FOXA1": 3.2, "TMPRSS2": 3.0,
            "PTHLH": 2.5, "IL6": 1.2, "VEGFA": 2.2, "CD47": 2.8, "CD274": 1.5,
            "BMPR1A": 1.8, "FZD5": 1.6, "CXCR4": 1.4, "ITGB1": 1.9
        },
        "t_NEPC": {
            "SYP": 3.6, "CHGA": 3.4, "NCAM1": 3.1, "ASCL1": 2.9, "NEUROD1": 2.8,
            "DLL1": 3.2, "PTHLH": 2.8, "VEGFA": 2.6, "CD47": 3.0, "CD44": 2.2
        },
        "myCAF": {
            "ACTA2": 3.7, "FAP": 3.5, "POSTN": 3.4, "COL1A1": 4.1, "COL1A2": 3.9,
            "TGFB1": 3.6, "CXCL12": 3.8, "IL6": 2.6, "WNT5A": 2.7, "BMP2": 2.0,
            "TGFBR1": 2.1, "TGFBR2": 2.3, "ITGA1": 2.4, "ITGB1": 2.8
        },
        "Endothelial": {
            "PECAM1": 3.8, "VWF": 3.6, "CDH5": 3.5,
            "FLT1": 3.2, "KDR": 3.0, "NOTCH1": 2.5, "CXCR4": 2.1,
            "TGFBR1": 1.9, "TGFBR2": 1.8, "SIRPA": 1.6
        },
        "TAM_M2": {
            "CD68": 3.7, "MRC1": 3.5, "CD163": 3.4, "ARG1": 3.0,
            "SPP1": 3.9, "TGFB1": 2.9, "IL6": 2.4, "CD274": 3.1, "SIRPA": 3.2,
            "IL6R": 2.3, "IL6ST": 2.6, "CD44": 2.5, "CXCR4": 2.0
        },
        "CD8_T": {
            "CD8A": 3.8, "CD3D": 3.7, "GZMB": 2.8, "PRF1": 2.6,
            "PDCD1": 3.2, "HAVCR2": 2.9, "LAG3": 2.7, "TOX": 2.8,
            "CXCR4": 3.1, "TGFBR1": 2.0, "TGFBR2": 2.2, "CD44": 2.4
        },
        "Osteoclast": {
            "ACP5": 3.9, "CTSK": 3.8, "NFATC1": 3.4, "CALCR": 3.2,
            "TNFRSF11A": 3.6, "SIRPA": 2.2, "ITGB1": 2.1, "CD44": 2.0
        },
        "Osteoblast": {
            "BGLAP": 3.8, "RUNX2": 3.5, "ALPL": 3.6,
            "TNFSF11": 3.7, "PTH1R": 3.4, "BMP2": 2.8, "COL1A1": 3.5,
            "TGFBR1": 2.2, "TGFBR2": 2.4, "LRP5": 2.8, "FZD5": 2.2
        },
    }

    cell_idx = 0
    for ctype, count in cell_type_distribution:
        sig = cell_signatures.get(ctype, {})
        for i in range(count):
            cell_id = f"Cell_{cell_idx:04d}_{ctype}"
            cell_idx += 1
            meta = {
                "cell_id": cell_id,
                "cell_type": ctype,
                "n_counts": int(random.gauss(4500, 500)),
                "n_genes": int(random.gauss(1800, 150)),
            }
            cell_metadata.append(meta)

            # Generate expression vector for 300 genes
            expr_row = []
            for g in all_genes:
                if g in sig:
                    mean_val = sig[g]
                    # Log-normal variation with dropout probability
                    val = max(0.0, random.gauss(mean_val, 0.45))
                    # Biological dropout for 10% of cells
                    if random.random() < 0.10:
                        val = 0.0
                else:
                    # Basal noise with high dropout
                    val = max(0.0, random.gauss(0.15, 0.20)) if random.random() < 0.25 else 0.0
                expr_row.append(round(val, 4))
            expression_matrix.append(expr_row)

    return cell_metadata, all_genes, expression_matrix


# ==============================================================================
# Step 2: Multi-Subunit L-R Database Curation
# ==============================================================================
def curate_lr_complex_database():
    """
    Curates realistic multi-subunit Ligand-Receptor interactions relevant to
    prostate cancer bone metastasis and immunosuppressive TME.
    """
    lr_database = [
        {
            "pair_name": "IL6_(IL6R+IL6ST)",
            "pathway": "IL6",
            "ligand_subunits": ["IL6"],
            "receptor_subunits": ["IL6R", "IL6ST"],
            "evidence": "Heterodimer gp130 signal transducer",
            "mode": "Paracrine/Autocrine",
        },
        {
            "pair_name": "TGFB1_(TGFBR1+TGFBR2)",
            "pathway": "TGFb",
            "ligand_subunits": ["TGFB1"],
            "receptor_subunits": ["TGFBR1", "TGFBR2"],
            "evidence": "Ser/Thr kinase heterotetrameric receptor",
            "mode": "Paracrine/Autocrine",
        },
        {
            "pair_name": "TNFSF11_TNFRSF11A",
            "pathway": "RANKL",
            "ligand_subunits": ["TNFSF11"],
            "receptor_subunits": ["TNFRSF11A"],
            "evidence": "Osteoblast-osteoclast bone remodeling axis",
            "mode": "Paracrine/Juxtacrine",
        },
        {
            "pair_name": "CD274_PDCD1",
            "pathway": "Checkpoint",
            "ligand_subunits": ["CD274"],
            "receptor_subunits": ["PDCD1"],
            "evidence": "Immune checkpoint inhibitory axis",
            "mode": "Juxtacrine",
        },
        {
            "pair_name": "CD47_SIRPA",
            "pathway": "Innate_Checkpoint",
            "ligand_subunits": ["CD47"],
            "receptor_subunits": ["SIRPA"],
            "evidence": "Don't eat me macrophage phagocytosis inhibition",
            "mode": "Juxtacrine",
        },
        {
            "pair_name": "CXCL12_CXCR4",
            "pathway": "CXCL",
            "ligand_subunits": ["CXCL12"],
            "receptor_subunits": ["CXCR4"],
            "evidence": "CAF stromal exclusion and tumor homing",
            "mode": "Paracrine",
        },
        {
            "pair_name": "SPP1_CD44",
            "pathway": "SPP1",
            "ligand_subunits": ["SPP1"],
            "receptor_subunits": ["CD44"],
            "evidence": "TAM-driven stemness and metastasis",
            "mode": "Paracrine",
        },
        {
            "pair_name": "VEGFA_(FLT1+KDR)",
            "pathway": "VEGF",
            "ligand_subunits": ["VEGFA"],
            "receptor_subunits": ["FLT1", "KDR"],
            "evidence": "Endothelial angiogenesis receptor complex",
            "mode": "Paracrine",
        },
        {
            "pair_name": "PTHLH_PTH1R",
            "pathway": "PTHrP",
            "ligand_subunits": ["PTHLH"],
            "receptor_subunits": ["PTH1R"],
            "evidence": "Tumor-osteoblast feedforward vicious cycle",
            "mode": "Paracrine",
        },
        {
            "pair_name": "COL1A1_(ITGA1+ITGB1)",
            "pathway": "Collagen",
            "ligand_subunits": ["COL1A1"],
            "receptor_subunits": ["ITGA1", "ITGB1"],
            "evidence": "Integrin α1β1 collagen ECM adhesive complex",
            "mode": "Juxtacrine/ECM",
        },
        {
            "pair_name": "DLL1_NOTCH1",
            "pathway": "NOTCH",
            "ligand_subunits": ["DLL1"],
            "receptor_subunits": ["NOTCH1"],
            "evidence": "Direct cell contact juxtacrine Notch signaling",
            "mode": "Juxtacrine",
        },
        {
            "pair_name": "WNT5A_(FZD5+LRP5)",
            "pathway": "WNT",
            "ligand_subunits": ["WNT5A"],
            "receptor_subunits": ["FZD5", "LRP5"],
            "evidence": "Frizzled-LRP co-receptor non-canonical Wnt",
            "mode": "Paracrine",
        },
        {
            "pair_name": "BMP2_BMPR1A",
            "pathway": "BMP",
            "ligand_subunits": ["BMP2"],
            "receptor_subunits": ["BMPR1A"],
            "evidence": "Osteosclerotic osteoblast differentiation",
            "mode": "Paracrine",
        },
    ]
    return lr_database


# ==============================================================================
# Step 3: Geometric Mean Expression & Hill Communication Probability
# ==============================================================================
def compute_trimmed_mean(values, trim_fraction=0.10):
    if not values:
        return 0.0
    s_vals = sorted(values)
    k = int(len(s_vals) * trim_fraction)
    if k > 0 and len(s_vals) > 2 * k:
        s_vals = s_vals[k:-k]
    return sum(s_vals) / len(s_vals)


def calculate_complex_expression(subunit_means):
    """
    Computes geometric mean across all subunits in a complex:
      E_complex = (∏_{i=1}^n E_i)^(1/n)
    """
    if not subunit_means:
        return 0.0
    prod = 1.0
    for val in subunit_means:
        if val <= 0:
            return 0.0
        prod *= val
    return math.pow(prod, 1.0 / len(subunit_means))


def compute_communication_probability(ligand_expr, receptor_expr, Kh=0.5):
    """
    CellChat Hill Cooperativity Function:
      P = (L * R) / (Kh + L * R)
    """
    s = ligand_expr * receptor_expr
    if s <= 0:
        return 0.0
    return s / (Kh + s)


def evaluate_cell_type_communication(cell_metadata, all_genes, expression_matrix, lr_database, Kh=0.5):
    """
    Computes communication probabilities across all 8 cell types for each L-R pair.
    """
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    cells_by_type = defaultdict(list)
    for c_idx, meta in enumerate(cell_metadata):
        cells_by_type[meta["cell_type"]].append(c_idx)

    cell_types = sorted(cells_by_type.keys())

    # Precalculate mean expression of required genes per cell type
    gene_mean_by_type = {}
    for ctype in cell_types:
        gene_mean_by_type[ctype] = {}
        c_indices = cells_by_type[ctype]
        for lr in lr_database:
            for g in lr["ligand_subunits"] + lr["receptor_subunits"]:
                if g in gene_to_idx and g not in gene_mean_by_type[ctype]:
                    g_idx = gene_to_idx[g]
                    vals = [expression_matrix[ci][g_idx] for ci in c_indices]
                    gene_mean_by_type[ctype][g] = compute_trimmed_mean(vals, trim_fraction=0.10)

    # Compute communication probabilities
    comm_results = []
    tensor_matrix = defaultdict(dict)  # (sender, receiver, lr_name) -> prob

    for lr in lr_database:
        pair_name = lr["pair_name"]
        pathway = lr["pathway"]
        l_subunits = lr["ligand_subunits"]
        r_subunits = lr["receptor_subunits"]

        for sender in cell_types:
            l_means = [gene_mean_by_type[sender].get(g, 0.0) for g in l_subunits]
            e_ligand = calculate_complex_expression(l_means)

            for receiver in cell_types:
                r_means = [gene_mean_by_type[receiver].get(g, 0.0) for g in r_subunits]
                e_receptor = calculate_complex_expression(r_means)

                prob = compute_communication_probability(e_ligand, e_receptor, Kh=Kh)
                tensor_matrix[(sender, receiver)][pair_name] = prob

                comm_results.append({
                    "sender": sender,
                    "receiver": receiver,
                    "pair_name": pair_name,
                    "pathway": pathway,
                    "ligand_subunits": "+".join(l_subunits),
                    "receptor_subunits": "+".join(r_subunits),
                    "ligand_expr": round(e_ligand, 4),
                    "receptor_expr": round(e_receptor, 4),
                    "comm_prob": round(prob, 5),
                })

    return comm_results, gene_mean_by_type, cell_types, gene_to_idx


# ==============================================================================
# Step 4: Cell-Label Permutation Testing (1,000 iterations) & LIANA Consensus
# ==============================================================================
def run_permutation_testing(cell_metadata, all_genes, expression_matrix, lr_database,
                            comm_results, n_permutations=1000, Kh=0.5, random_seed=42):
    """
    Shuffles cell-type labels 1,000 times to construct an empirical null distribution
    and compute p-values and Benjamini-Hochberg FDR.
    """
    print(f"[*] Running Cell-Label Permutation Test with {n_permutations} iterations...")
    random.seed(random_seed)

    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    orig_labels = [m["cell_type"] for m in cell_metadata]
    n_cells = len(orig_labels)

    # Find unique genes needed for L-R pairs
    active_genes = sorted(list({g for lr in lr_database for g in (lr["ligand_subunits"] + lr["receptor_subunits"])}))
    active_gene_indices = [gene_to_idx[g] for g in active_genes]
    g_to_local_idx = {g: i for i, g in enumerate(active_genes)}

    # Pre-extract active gene expression matrix
    active_expr = [[row[gi] for gi in active_gene_indices] for row in expression_matrix]

    cell_types = sorted(list(set(orig_labels)))
    ctype_to_id = {ct: i for i, ct in enumerate(cell_types)}
    n_ctypes = len(cell_types)
    orig_type_ids = [ctype_to_id[l] for l in orig_labels]

    # Pre-calculate counts per type
    ctype_counts = [0] * n_ctypes
    for tid in orig_type_ids:
        ctype_counts[tid] += 1

    # Map observable probabilities into lookup key
    obs_map = {}
    for r in comm_results:
        key = (r["sender"], r["receiver"], r["pair_name"])
        obs_map[key] = r["comm_prob"]

    # Null exceedance counters
    exceed_counts = defaultdict(int)

    # Pre-index L-R sub-indices
    lr_specs = []
    for lr in lr_database:
        l_locs = [g_to_local_idx[g] for g in lr["ligand_subunits"]]
        r_locs = [g_to_local_idx[g] for g in lr["receptor_subunits"]]
        lr_specs.append((lr["pair_name"], l_locs, r_locs))

    t_start = time.time()
    # Fast permutation loop
    for b in range(n_permutations):
        shuffled_ids = random.sample(orig_type_ids, n_cells)

        # Accumulate sums
        sums = [[0.0] * len(active_genes) for _ in range(n_ctypes)]
        for ci, tid in enumerate(shuffled_ids):
            c_row = active_expr[ci]
            t_row = sums[tid]
            for gi, val in enumerate(c_row):
                t_row[gi] += val

        # Means per type
        means = [[sums[tid][gi] / ctype_counts[tid] for gi in range(len(active_genes))] for tid in range(n_ctypes)]

        # Check exceedance
        for pair_name, l_locs, r_locs in lr_specs:
            for s_id, sender in enumerate(cell_types):
                # Ligand geometric mean
                l_prod = 1.0
                for loc in l_locs:
                    l_prod *= means[s_id][loc]
                e_lig = math.pow(l_prod, 1.0 / len(l_locs)) if l_prod > 0 else 0.0

                for r_id, receiver in enumerate(cell_types):
                    # Receptor geometric mean
                    r_prod = 1.0
                    for loc in r_locs:
                        r_prod *= means[r_id][loc]
                    e_rec = math.pow(r_prod, 1.0 / len(r_locs)) if r_prod > 0 else 0.0

                    s_val = e_lig * e_rec
                    null_p = s_val / (Kh + s_val) if s_val > 0 else 0.0

                    key = (sender, receiver, pair_name)
                    if null_p >= obs_map[key]:
                        exceed_counts[key] += 1

        if (b + 1) % 250 == 0:
            print(f"  [>] Permutation {b + 1}/{n_permutations} completed ({round(time.time() - t_start, 2)}s)")

    # Calculate empirical p-values
    p_results = []
    for r in comm_results:
        key = (r["sender"], r["receiver"], r["pair_name"])
        exceed = exceed_counts[key]
        p_val = (exceed + 1.0) / (n_permutations + 1.0)
        p_results.append({
            **r,
            "null_exceed_count": exceed,
            "p_value": p_val,
        })

    # Benjamini-Hochberg FDR correction
    p_results.sort(key=lambda x: x["p_value"])
    m_tests = len(p_results)
    for rank_idx, item in enumerate(p_results):
        rank = rank_idx + 1
        fdr = min(1.0, item["p_value"] * m_tests / rank)
        item["rank_p"] = rank
        item["fdr"] = fdr

    # Ensure monotonic FDR
    min_fdr = 1.0
    for item in reversed(p_results):
        min_fdr = min(min_fdr, item["fdr"])
        item["fdr"] = round(min_fdr, 5)
        item["p_value"] = round(item["p_value"], 5)

    # LIANA Consensus Borda Aggregation
    # Rank by comm_prob (descending)
    p_results.sort(key=lambda x: x["comm_prob"], reverse=True)
    for rank_idx, item in enumerate(p_results):
        item["rank_prob"] = rank_idx + 1

    # Consensus score: average rank (lower is better)
    for item in p_results:
        item["borda_score"] = round((item["rank_prob"] + item["rank_p"]) / 2.0, 2)
        item["is_significant"] = (item["fdr"] < 0.05 and item["comm_prob"] >= 0.10)

    # Final sort by Borda score
    p_results.sort(key=lambda x: (x["borda_score"], -x["comm_prob"]))
    return p_results


# ==============================================================================
# Step 5: 2D Network Centrality Analysis (Sender, Receiver, Mediator, Influencer)
# ==============================================================================
def compute_network_centrality(consensus_results, cell_types):
    """
    Computes Out-degree, In-degree, Betweenness Centrality, and PageRank
    for the cell-cell communication network aggregated across significant signals.
    """
    n_nodes = len(cell_types)
    ctype_to_idx = {ct: i for i, ct in enumerate(cell_types)}

    # Weighted adjacency matrix W[u][v]
    W = [[0.0] * n_nodes for _ in range(n_nodes)]
    count_matrix = [[0] * n_nodes for _ in range(n_nodes)]

    for r in consensus_results:
        if r["is_significant"]:
            u = ctype_to_idx[r["sender"]]
            v = ctype_to_idx[r["receiver"]]
            W[u][v] += r["comm_prob"]
            count_matrix[u][v] += 1

    # 1. Out-degree (Sender Strength) & In-degree (Receiver Strength)
    out_degree = [sum(W[u][v] for v in range(n_nodes)) for u in range(n_nodes)]
    in_degree = [sum(W[u][v] for u in range(n_nodes)) for v in range(n_nodes)]

    # 2. Betweenness Centrality (Brandes' Algorithm with edge distance = 1 / (W_uv + 1e-4))
    betweenness = [0.0] * n_nodes
    dist_matrix = [[float("inf")] * n_nodes for _ in range(n_nodes)]
    for u in range(n_nodes):
        dist_matrix[u][u] = 0.0
        for v in range(n_nodes):
            if W[u][v] > 0:
                dist_matrix[u][v] = 1.0 / (W[u][v] + 0.001)

    for s in range(n_nodes):
        # Dijkstra
        S = []
        P = [[] for _ in range(n_nodes)]
        sigma = [0.0] * n_nodes
        sigma[s] = 1.0
        d = [float("inf")] * n_nodes
        d[s] = 0.0
        visited = [False] * n_nodes

        for _ in range(n_nodes):
            # pick min unvisited
            u = -1
            min_d = float("inf")
            for i in range(n_nodes):
                if not visited[i] and d[i] < min_d:
                    min_d = d[i]
                    u = i
            if u == -1 or min_d == float("inf"):
                break
            visited[u] = True
            S.append(u)

            for v in range(n_nodes):
                if W[u][v] > 0:
                    alt = d[u] + dist_matrix[u][v]
                    if alt < d[v]:
                        d[v] = alt
                        sigma[v] = sigma[u]
                        P[v] = [u]
                    elif math.isclose(alt, d[v], rel_tol=1e-5):
                        sigma[v] += sigma[u]
                        P[v].append(u)

        # Accumulation
        delta = [0.0] * n_nodes
        while S:
            w = S.pop()
            for u in P[w]:
                delta[u] += (sigma[u] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                betweenness[w] += delta[w]

    # Normalize betweenness: / ((N-1)(N-2))
    norm_factor = (n_nodes - 1) * (n_nodes - 2) if n_nodes > 2 else 1.0
    norm_betweenness = [round(b / norm_factor, 5) for b in betweenness]

    # 3. PageRank (Influencer Score) via Power Iteration
    # Stochastic matrix M
    damp = 0.85
    pr = [1.0 / n_nodes] * n_nodes
    M = [[0.0] * n_nodes for _ in range(n_nodes)]
    for u in range(n_nodes):
        row_sum = sum(W[u])
        for v in range(n_nodes):
            M[u][v] = W[u][v] / row_sum if row_sum > 0 else 1.0 / n_nodes

    for _ in range(100):
        new_pr = [(1.0 - damp) / n_nodes] * n_nodes
        for v in range(n_nodes):
            for u in range(n_nodes):
                new_pr[v] += damp * pr[u] * M[u][v]
        diff = sum(abs(new_pr[i] - pr[i]) for i in range(n_nodes))
        pr = new_pr
        if diff < 1e-6:
            break

    # Compile Centrality Table
    centrality_table = []
    for i, ctype in enumerate(cell_types):
        centrality_table.append({
            "cell_type": ctype,
            "out_degree_sender": round(out_degree[i], 4),
            "in_degree_receiver": round(in_degree[i], 4),
            "betweenness_mediator": norm_betweenness[i],
            "pagerank_influencer": round(pr[i], 5),
            "role_archetype": assign_archetype(out_degree[i], in_degree[i], norm_betweenness[i]),
        })

    centrality_table.sort(key=lambda x: x["out_degree_sender"] + x["in_degree_receiver"], reverse=True)
    return centrality_table


def assign_archetype(out_d, in_d, betw):
    if out_d > 1.5 * in_d:
        return "Dominant Sender (Signal Source)"
    elif in_d > 1.5 * out_d:
        return "Dominant Receiver (Signal Target)"
    elif betw > 0.15:
        return "Central Mediator (Hub & Relay)"
    else:
        return "Balanced Communicator"


# ==============================================================================
# Main Pipeline Runner & Result Exporter
# ==============================================================================
def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "pipelines", "results_pipeline_01")
    os.makedirs(outdir, exist_ok=True)

    print("================================================================================")
    print("Track 03 - Pipeline 01: LIANA Consensus & CellChat Communication Analysis")
    print("================================================================================")

    # 1. Simulate Cohort
    print("[1/5] Simulating 1,500 mCRPC single cells across 8 cell types...")
    cell_metadata, all_genes, expression_matrix = simulate_mcrpc_cohort(n_cells=1500, random_seed=42)
    meta_path = os.path.join(outdir, "simulated_mcrpc_cell_metadata.csv")
    with open(meta_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["cell_id", "cell_type", "n_counts", "n_genes"])
        writer.writeheader()
        writer.writerows(cell_metadata)
    print(f"      [✓] Exported metadata: {meta_path} ({len(cell_metadata)} cells)")

    # 2. Curate L-R Database
    print("[2/5] Curating multi-subunit L-R database for prostate bone metastasis...")
    lr_database = curate_lr_complex_database()
    lr_path = os.path.join(outdir, "curated_lr_complex_database.csv")
    with open(lr_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["pair_name", "pathway", "ligand_subunits", "receptor_subunits", "evidence", "mode"])
        writer.writeheader()
        for r in lr_database:
            writer.writerow({
                **r,
                "ligand_subunits": "+".join(r["ligand_subunits"]),
                "receptor_subunits": "+".join(r["receptor_subunits"]),
            })
    print(f"      [✓] Exported L-R DB: {lr_path} ({len(lr_database)} pairs)")

    # 3. Compute Communication Probabilities
    print("[3/5] Computing trimmed subunit geometric means and Hill communication probabilities (Kh=0.5)...")
    comm_results, gene_means, cell_types, gene_to_idx = evaluate_cell_type_communication(
        cell_metadata, all_genes, expression_matrix, lr_database, Kh=0.5
    )
    prob_path = os.path.join(outdir, "cellchat_lr_communication_probabilities.csv")
    with open(prob_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sender", "receiver", "pair_name", "pathway", "ligand_subunits",
            "receptor_subunits", "ligand_expr", "receptor_expr", "comm_prob"
        ])
        writer.writeheader()
        writer.writerows(comm_results)
    print(f"      [✓] Exported raw probabilities: {prob_path} ({len(comm_results)} interactions)")

    # 4. Run Permutations & Consensus
    print("[4/5] Executing 1,000 cell-label permutations for empirical FDR & LIANA Borda aggregation...")
    consensus_results = run_permutation_testing(
        cell_metadata, all_genes, expression_matrix, lr_database, comm_results,
        n_permutations=1000, Kh=0.5, random_seed=42
    )
    cons_path = os.path.join(outdir, "liana_consensus_permutation_results.csv")
    with open(cons_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sender", "receiver", "pair_name", "pathway", "comm_prob",
            "p_value", "fdr", "rank_prob", "rank_p", "borda_score", "is_significant"
        ], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(consensus_results)
    sig_count = sum(1 for r in consensus_results if r["is_significant"])
    print(f"      [✓] Exported LIANA consensus results: {cons_path} ({sig_count}/{len(consensus_results)} significant)")

    # 5. Network Centrality
    print("[5/5] Calculating 2D network centrality (Out/In-degree, Betweenness, PageRank)...")
    centrality_table = compute_network_centrality(consensus_results, cell_types)
    cent_path = os.path.join(outdir, "cell_type_network_centrality_metrics.csv")
    with open(cent_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "cell_type", "out_degree_sender", "in_degree_receiver",
            "betweenness_mediator", "pagerank_influencer", "role_archetype"
        ])
        writer.writeheader()
        writer.writerows(centrality_table)
    print(f"      [✓] Exported centrality table: {cent_path}")

    # Summary Metrics JSON
    summary_metrics = {
        "n_cells_simulated": len(cell_metadata),
        "n_cell_types": len(cell_types),
        "n_genes": len(all_genes),
        "n_lr_pairs": len(lr_database),
        "n_total_interactions_tested": len(comm_results),
        "n_significant_interactions_fdr05": sig_count,
        # top sender는 out-degree 기준(표 정렬 키인 out+in 합과 다를 수 있음), top receiver는 in-degree 기준
        "top_sender_cell_type": max(centrality_table, key=lambda x: x["out_degree_sender"])["cell_type"],
        "top_sender_strength": max(x["out_degree_sender"] for x in centrality_table),
        "top_receiver_cell_type": max(centrality_table, key=lambda x: x["in_degree_receiver"])["cell_type"],
        "top_receiver_strength": max(x["in_degree_receiver"] for x in centrality_table),
        "top_mediator_cell_type": max(centrality_table, key=lambda x: x["betweenness_mediator"])["cell_type"],
        "top_influencer_cell_type": max(centrality_table, key=lambda x: x["pagerank_influencer"])["cell_type"],
        "top_consensus_lr_pairs": [
            f"{r['sender']} -> {r['receiver']} ({r['pair_name']}): prob={r['comm_prob']}, fdr={r['fdr']}"
            for r in consensus_results[:5]
        ],
    }
    json_path = os.path.join(outdir, "pipeline_summary_metrics.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary_metrics, f, indent=2)
    print(f"      [✓] Exported pipeline summary: {json_path}")

    print("\n[✓] Pipeline 01 successfully completed all steps with exit code 0.")


if __name__ == "__main__":
    main()

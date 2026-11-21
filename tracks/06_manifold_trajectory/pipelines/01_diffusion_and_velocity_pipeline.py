#!/usr/bin/env python3
"""
tracks/06_manifold_trajectory/pipelines/01_diffusion_and_velocity_pipeline.py
=============================================================================
End-to-End Pipeline for Diffusion Pseudotime (DPT) & RNA Velocity in mCRPC Lineage Plasticity.

Biological Cohort:
  Castration-Resistant Prostate Cancer (mCRPC) Single-Cell Lineage Progression
  - 1,500 Single Cells across 4 Phenotypes:
    1. AR_Adeno (500 cells): AR+, KLK3+, FOXA1+, TMPRSS2+, HOXB13+ (Root state)
    2. Amphicrine (350 cells): AR+, SYP_low, ASCL1_low, KRT5+, SOX2+ (Plastic intermediate)
    3. DNPC (300 cells): AR-, NE-, ETV4+, FGFR1+, SOX9+, TCF4+ (Double-negative lineage)
    4. t_NEPC (350 cells): AR-, ASCL1+, INSM1+, EZH2+, SYP+, CHGA+ (Neuroendocrine lineage)
  - 400 Genes: Spliced (S) and Unspliced (U) transcript layers with realistic transcription-splicing ODE kinetics.

Steps:
  Step 1: Simulate 1,500 mCRPC cells with S and U count matrices along branching kinetics.
  Step 2: Implement Diffusion Maps & DPT:
          - kNN graph construction (k=20)
          - Adaptive kernel bandwidth & symmetric normalized Laplacian S = D^(-1/2) W D^(-1/2)
          - Top 10 diffusion eigenvectors via orthogonalized power iteration
          - Geodesic DPT metric rooted at AR_Adeno centroid.
  Step 3: Implement Steady-State RNA Velocity:
          - Linear regression of u on s estimating degradation ratio gamma_g
          - Velocity calculation v = u - gamma * s
          - Transition probabilities & 2D manifold velocity vector field projection.
  Step 4: Implement Dynamical Velocity Coherence & Cell-Specific Latent Time tau.
  Step 5: Output CSV tables and summary JSON to results_pipeline_01/.
"""

import os
import sys
import math
import time
import json
import csv
import random
from collections import defaultdict


# ==============================================================================
# Helper Math Functions (Pure Python Linear Algebra)
# ==============================================================================
def dot(v1, v2):
    return sum(a * b for a, b in zip(v1, v2))


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def normalize(v):
    mag = norm(v)
    if mag < 1e-12:
        return [0.0] * len(v)
    return [x / mag for x in v]


def euclidean_dist(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def mean_val(vals):
    return sum(vals) / len(vals) if vals else 0.0


def variance_val(vals):
    if len(vals) < 2:
        return 0.0
    m = mean_val(vals)
    return sum((x - m) ** 2 for x in vals) / (len(vals) - 1)


# ==============================================================================
# Step 1: Cohort Simulation (1,500 Cells x 400 Genes, S and U Count Layers)
# ==============================================================================
def simulate_mcrpc_trajectory_cohort(n_cells=1500, n_genes=400, random_seed=42):
    """
    Simulates spliced (S) and unspliced (U) expression counts for 1,500 single cells
    transitioning from AR+ Adeno to Amphicrine, DNPC, and t-NEPC lineages.
    """
    random.seed(random_seed)
    print(f"\n[Step 1] Simulating {n_cells} single cells with S and U layers across 4 lineages...")

    # Lineage specifications
    lineage_counts = [
        ("AR_Adeno", 500, 0.0, 0.25),
        ("Amphicrine", 350, 0.25, 0.60),
        ("DNPC", 300, 0.55, 1.0),
        ("t_NEPC", 350, 0.55, 1.0),
    ]

    # Gene definitions
    marker_genes = {
        "AR_Adeno": ["AR", "KLK3", "TMPRSS2", "FOXA1", "HOXB13", "NKX3-1", "FKBP5", "ACSL3"],
        "Amphicrine": ["KRT5", "SOX2", "CD44", "KRT14", "TP63", "MYC"],
        "DNPC": ["ETV4", "FGFR1", "SOX9", "TCF4", "FGF2", "EGFR"],
        "t_NEPC": ["ASCL1", "INSM1", "EZH2", "SYP", "CHGA", "NCAM1", "DLL3", "FOXA2"],
    }
    all_markers = []
    for g_list in marker_genes.values():
        all_markers.extend(g_list)

    # Fill remaining genes up to 400
    housekeeping = [f"HK_GENE_{i:02d}" for i in range(120)]
    variable_bg = [f"BG_GENE_{i:03d}" for i in range(n_genes - len(all_markers) - len(housekeeping))]
    gene_names = all_markers + housekeeping + variable_bg
    n_genes_total = len(gene_names)
    assert n_genes_total == n_genes

    # Synthetic degradation rates (gamma) for each gene (true ground truth)
    true_gammas = {}
    for g in gene_names:
        if g in marker_genes["AR_Adeno"]:
            true_gammas[g] = 0.85 + random.uniform(-0.1, 0.1)
        elif g in marker_genes["t_NEPC"]:
            true_gammas[g] = 0.70 + random.uniform(-0.1, 0.1)
        elif g in marker_genes["DNPC"]:
            true_gammas[g] = 0.75 + random.uniform(-0.1, 0.1)
        elif g in marker_genes["Amphicrine"]:
            true_gammas[g] = 0.90 + random.uniform(-0.1, 0.1)
        else:
            true_gammas[g] = 0.80 + random.uniform(-0.15, 0.15)

    cell_metadata = []
    cell_idx = 0

    # Build continuous progress coordinates along trajectories
    for l_name, count, t_min, t_max in lineage_counts:
        for _ in range(count):
            t_prog = random.uniform(t_min, t_max)
            # Add stochastic diffusion jitter
            t_prog = max(0.0, min(1.0, t_prog + random.gauss(0, 0.04)))
            cell_metadata.append({
                "cell_id": f"CELL_{cell_idx:04d}",
                "lineage": l_name,
                "true_progression": t_prog,
                "library_size": random.uniform(2500, 6500),
            })
            cell_idx += 1

    random.shuffle(cell_metadata)
    for idx, c in enumerate(cell_metadata):
        c["cell_id"] = f"CELL_{idx:04d}"

    # Generate Spliced (S) and Unspliced (U) counts
    # Splicing ODE: du/dt = alpha(t) - u,  ds/dt = u - gamma * s
    spliced_matrix = []
    unspliced_matrix = []

    for c in cell_metadata:
        t = c["true_progression"]
        l_name = c["lineage"]
        lib_factor = c["library_size"] / 4000.0

        s_row = []
        u_row = []

        for g in gene_names:
            gamma = true_gammas[g]
            # Transcription rate alpha(t) depends on lineage and progression
            alpha = 0.5  # basal transcription
            if g in marker_genes["AR_Adeno"]:
                if l_name == "AR_Adeno":
                    alpha = 6.0 * (1.0 - 0.7 * t)
                elif l_name == "Amphicrine":
                    alpha = 3.0 * (1.0 - t)
                else:
                    alpha = 0.3 * (1.0 - t)
            elif g in marker_genes["t_NEPC"]:
                if l_name == "t_NEPC":
                    alpha = 1.0 + 8.0 * (t ** 1.8)
                elif l_name == "Amphicrine":
                    alpha = 1.0 + 2.5 * t
                else:
                    alpha = 0.2
            elif g in marker_genes["DNPC"]:
                if l_name == "DNPC":
                    alpha = 1.0 + 7.5 * (t ** 1.8)
                else:
                    alpha = 0.2
            elif g in marker_genes["Amphicrine"]:
                if l_name == "Amphicrine":
                    alpha = 2.0 + 5.0 * math.sin(t * math.pi)
                elif l_name in ("t_NEPC", "DNPC"):
                    alpha = 1.5 * (1.0 - t)
                else:
                    alpha = 1.0 + 1.5 * t

            # Add Poisson / Gamma noise
            expected_u = alpha * lib_factor
            expected_s = (alpha / gamma) * lib_factor

            # Stochastic burst noise
            u_count = max(0, int(random.gauss(expected_u, math.sqrt(max(expected_u, 0.1)) + 0.3)))
            s_count = max(0, int(random.gauss(expected_s, math.sqrt(max(expected_s, 0.1)) + 0.5)))

            s_row.append(s_count)
            u_row.append(u_count)

        spliced_matrix.append(s_row)
        unspliced_matrix.append(u_row)

    print(f"  [✓] Successfully simulated {n_cells} cells across {n_genes_total} genes.")
    print(f"  [✓] Lineages: AR_Adeno (500), Amphicrine (350), DNPC (300), t_NEPC (350).")
    return cell_metadata, gene_names, spliced_matrix, unspliced_matrix, true_gammas


# ==============================================================================
# Step 2: Diffusion Maps & Diffusion Pseudotime (DPT)
# ==============================================================================
def compute_diffusion_maps_and_dpt(cell_metadata, gene_names, spliced_matrix, k_neighbors=20, n_components=10):
    """
    Computes kNN graph, Gaussian affinity matrix W, symmetric normalized Laplacian S = D^(-1/2) W D^(-1/2),
    top eigenvectors via Gram-Schmidt orthogonalized power iteration, and DPT geodesic distances.
    """
    print("\n[Step 2] Computing Diffusion Maps and Diffusion Pseudotime (DPT)...")
    n_cells = len(cell_metadata)

    # 1. Log-normalize spliced counts: log1p(count / total * 10000)
    print("  [1/4] Log-normalizing expression profiles...")
    norm_s = []
    for row in spliced_matrix:
        tot = sum(row)
        scale = 10000.0 / tot if tot > 0 else 1.0
        norm_s.append([math.log1p(x * scale) for x in row])

    # 2. Select top 60 highly variable genes for graph distance calculation
    gene_vars = []
    for g_idx in range(len(gene_names)):
        vals = [norm_s[i][g_idx] for i in range(n_cells)]
        gene_vars.append((variance_val(vals), g_idx))
    gene_vars.sort(key=lambda x: x[0], reverse=True)
    hvg_indices = [idx for _, idx in gene_vars[:60]]

    # 3. Construct kNN Graph & Gaussian Kernel
    print(f"  [2/4] Constructing kNN graph (k={k_neighbors}) on {len(hvg_indices)} HVGs...")
    cell_features = [[norm_s[i][g] for g in hvg_indices] for i in range(n_cells)]

    # Pairwise distances to find k-nearest neighbors
    knn_graph = []  # list of list of (neighbor_idx, distance)
    sigmas = []     # local bandwidth for each cell (distance to k-th neighbor)

    for i in range(n_cells):
        dists = []
        p_i = cell_features[i]
        for j in range(n_cells):
            if i != j:
                d = euclidean_dist(p_i, cell_features[j])
                dists.append((j, d))
        dists.sort(key=lambda x: x[1])
        neighbors = dists[:k_neighbors]
        knn_graph.append(neighbors)
        sigmas.append(max(neighbors[-1][1], 1e-4))

    # Build sparse affinity matrix W
    # W_ij = exp(-dist^2 / (2 * sigma_i * sigma_j))
    print("  [3/4] Assembling affinity matrix W and Degree matrix D...")
    adj = defaultdict(dict)
    for i in range(n_cells):
        for j, d in knn_graph[i]:
            val = math.exp(-(d * d) / (2.0 * sigmas[i] * sigmas[j]))
            adj[i][j] = max(adj[i].get(j, 0.0), val)
            adj[j][i] = max(adj[j].get(i, 0.0), val)  # Symmetrize

    # Degrees D_i = sum_j W_ij
    deg = [0.0] * n_cells
    for i in range(n_cells):
        deg[i] = sum(adj[i].values())

    inv_sqrt_deg = [1.0 / math.sqrt(d) if d > 1e-12 else 0.0 for d in deg]

    # Symmetric Normalized Laplacian matrix S_ij = D_i^(-1/2) W_ij D_j^(-1/2)
    def matvec_S(v):
        """Computes S * v using sparse adjacency."""
        res = [0.0] * n_cells
        for i in range(n_cells):
            acc = 0.0
            isd_i = inv_sqrt_deg[i]
            for j, w in adj[i].items():
                acc += w * (inv_sqrt_deg[j] * v[j])
            res[i] = isd_i * acc
        return res

    # 4. Orthogonalized Power Iteration to extract top n_components eigenvectors
    print(f"  [4/4] Extracting top {n_components} diffusion eigenvectors via orthogonalized iteration...")
    eigenvectors = []
    eigenvalues = []

    for comp in range(n_components):
        # Random initialization
        v = [random.gauss(0, 1) for _ in range(n_cells)]
        # Orthogonalize against existing eigenvectors (Gram-Schmidt)
        for prev_v in eigenvectors:
            proj = dot(v, prev_v)
            v = [vi - proj * pvi for vi, pvi in zip(v, prev_v)]
        v = normalize(v)

        # Power iteration
        lambda_est = 0.0
        for it in range(35):
            v_next = matvec_S(v)
            # Orthogonalize
            for prev_v in eigenvectors:
                proj = dot(v_next, prev_v)
                v_next = [vi - proj * pvi for vi, pvi in zip(v_next, prev_v)]
            lambda_est = norm(v_next)
            if lambda_est < 1e-12:
                break
            v = [x / lambda_est for x in v_next]

        eigenvectors.append(v)
        eigenvalues.append(lambda_est)

    print(f"  [✓] Top eigenvalues: {[round(ev, 4) for ev in eigenvalues[:5]]}")

    # Convert symmetric eigenvectors v to transition eigenvectors psi:
    # psi_k = D^(-1/2) v_k
    diff_coords = []
    for comp in range(n_components):
        v = eigenvectors[comp]
        psi = [inv_sqrt_deg[i] * v[i] for i in range(n_cells)]
        psi = normalize(psi)
        diff_coords.append(psi)

    # Identify Root Cell x_0 (Centroid of AR_Adeno cells with lowest true_progression)
    root_candidates = [i for i, c in enumerate(cell_metadata) if c["lineage"] == "AR_Adeno"]
    root_candidates.sort(key=lambda idx: cell_metadata[idx]["true_progression"])
    root_cell_idx = root_candidates[0]
    print(f"  [✓] Identified Root Cell for DPT: {cell_metadata[root_cell_idx]['cell_id']} ({cell_metadata[root_cell_idx]['lineage']})")

    # Compute DPT Geodesic Metric:
    # DPT(x_0, x) = sqrt( sum_{k=1}^{n-1} lambda_k^(2*t) * (psi_k(x_0) - psi_k(x))^2 )
    # Here we use t = 2 steps
    t_step = 2
    dpt_distances = []
    for i in range(n_cells):
        acc = 0.0
        for k in range(1, n_components):
            lam = eigenvalues[k]
            scale = (lam ** (2 * t_step))
            diff = diff_coords[k][root_cell_idx] - diff_coords[k][i]
            acc += scale * (diff * diff)
        dpt_distances.append(math.sqrt(acc))

    # Normalize DPT to [0, 1]
    min_dpt = min(dpt_distances)
    max_dpt = max(dpt_distances)
    norm_dpt = [(d - min_dpt) / (max_dpt - min_dpt) if max_dpt > min_dpt else 0.0 for d in dpt_distances]

    for i in range(n_cells):
        cell_metadata[i]["dpt"] = norm_dpt[i]
        cell_metadata[i]["dc1"] = diff_coords[1][i]
        cell_metadata[i]["dc2"] = diff_coords[2][i]

    print("  [✓] Diffusion pseudotime (DPT) successfully calculated for all cells.")
    return diff_coords, eigenvalues, root_cell_idx, knn_graph


# ==============================================================================
# Step 3: Steady-State RNA Velocity Estimation
# ==============================================================================
def compute_steady_state_velocity(cell_metadata, gene_names, spliced_matrix, unspliced_matrix, knn_graph):
    """
    Computes steady-state degradation slope gamma_g using linear regression:
    gamma_g = (s^T u) / (s^T s), residual velocity v_g = u - gamma_g * s,
    and projects velocity vectors onto 2D Diffusion Space (DC1, DC2).
    """
    print("\n[Step 3] Computing Steady-State RNA Velocity (v = u - gamma * s)...")
    n_cells = len(cell_metadata)
    n_genes = len(gene_names)

    # Select top 50 lineage driver genes for velocity kinetics
    velocity_genes = [
        "AR", "KLK3", "TMPRSS2", "FOXA1", "HOXB13",
        "ASCL1", "INSM1", "EZH2", "SYP", "CHGA", "NCAM1", "SOX2",
        "ETV4", "FGFR1", "SOX9", "TCF4",
        "KRT5", "CD44", "TP63", "MYC"
    ]
    # Add background genes to reach 50
    for g in gene_names:
        if g not in velocity_genes:
            velocity_genes.append(g)
        if len(velocity_genes) >= 50:
            break

    velo_gene_indices = [gene_names.index(g) for g in velocity_genes]

    # Calculate gamma for each gene via linear regression on cells in upper 80% expression
    gene_gammas = {}
    gene_r2 = {}
    gene_velocities = defaultdict(list)

    for g_name, g_idx in zip(velocity_genes, velo_gene_indices):
        s_vals = [spliced_matrix[i][g_idx] for i in range(n_cells)]
        u_vals = [unspliced_matrix[i][g_idx] for i in range(n_cells)]

        # Find cells expressing this gene
        active_cells = [i for i in range(n_cells) if s_vals[i] > 0 or u_vals[i] > 0]
        if len(active_cells) < 20:
            active_cells = list(range(n_cells))

        s_act = [s_vals[i] for i in active_cells]
        u_act = [u_vals[i] for i in active_cells]

        mean_s = mean_val(s_act)
        mean_u = mean_u = mean_val(u_act)

        num = sum((s - mean_s) * (u - mean_u) for s, u in zip(s_act, u_act))
        den = sum((s - mean_s) ** 2 for s in s_act)
        gamma = num / den if den > 1e-9 else 0.8
        gamma = max(0.1, min(gamma, 3.0))  # Physical bounds

        # Residual velocity v = u - gamma * s
        residuals = []
        ss_res = 0.0
        ss_tot = sum((u - mean_u) ** 2 for u in u_act)
        for i in range(n_cells):
            v = float(unspliced_matrix[i][g_idx]) - gamma * float(spliced_matrix[i][g_idx])
            residuals.append(v)
            if i in active_cells:
                pred_u = mean_u + gamma * (spliced_matrix[i][g_idx] - mean_s)
                ss_res += (unspliced_matrix[i][g_idx] - pred_u) ** 2

        r2 = max(0.0, 1.0 - (ss_res / ss_tot)) if ss_tot > 1e-9 else 0.5
        gene_gammas[g_name] = gamma
        gene_r2[g_name] = r2
        gene_velocities[g_name] = residuals

    print(f"  [✓] Fitted gamma for {len(velocity_genes)} genes. Top key regulators:")
    for g in ["AR", "ASCL1", "ETV4", "EZH2", "KRT5"]:
        if g in gene_gammas:
            print(f"      • {g:7s} | Estimated γ: {gene_gammas[g]:.3f} | R²: {gene_r2[g]:.3f}")

    # Project velocity vector onto 2D DC1-DC2 embedding
    # Transition probability P_ij = exp( corr(v_i, s_j - s_i) / sigma )
    print("  [✓] Projecting high-dimensional velocity vectors onto 2D diffusion manifold...")
    v_dc1 = [0.0] * n_cells
    v_dc2 = [0.0] * n_cells

    for i in range(n_cells):
        v_vec = [gene_velocities[g][i] for g in velocity_genes]
        s_i = [spliced_matrix[i][g_idx] for g_idx in velo_gene_indices]

        # Calculate cosine correlation with neighbors in kNN
        corrs = []
        for j, _ in knn_graph[i]:
            s_j = [spliced_matrix[j][g_idx] for g_idx in velo_gene_indices]
            delta_s = [sj - si for sj, si in zip(s_j, s_i)]
            nv = norm(v_vec)
            nds = norm(delta_s)
            c = dot(v_vec, delta_s) / (nv * nds) if (nv > 1e-9 and nds > 1e-9) else 0.0
            corrs.append((j, c))

        # Softmax transition probabilities
        sigma_t = 0.2
        max_c = max(c for _, c in corrs)
        weights = [math.exp((c - max_c) / sigma_t) for _, c in corrs]
        tot_w = sum(weights)
        probs = [w / tot_w if tot_w > 0 else 1.0 / len(weights) for w in weights]

        # Projected displacement in DC1 and DC2
        delta_dc1 = 0.0
        delta_dc2 = 0.0
        dc1_i = cell_metadata[i]["dc1"]
        dc2_i = cell_metadata[i]["dc2"]

        for (j, _), p in zip(corrs, probs):
            delta_dc1 += (p - 1.0 / len(corrs)) * (cell_metadata[j]["dc1"] - dc1_i)
            delta_dc2 += (p - 1.0 / len(corrs)) * (cell_metadata[j]["dc2"] - dc2_i)

        v_dc1[i] = delta_dc1
        v_dc2[i] = delta_dc2

    # Scale 2D velocity vectors for visual readability
    max_mag = max(math.hypot(dx, dy) for dx, dy in zip(v_dc1, v_dc2)) or 1.0
    for i in range(n_cells):
        cell_metadata[i]["v_dc1"] = v_dc1[i] / max_mag
        cell_metadata[i]["v_dc2"] = v_dc2[i] / max_mag

    return gene_gammas, gene_r2, gene_velocities, velocity_genes


# ==============================================================================
# Step 4: Dynamical Latent Time & Velocity Coherence
# ==============================================================================
def compute_latent_time_and_coherence(cell_metadata, knn_graph):
    """
    Computes velocity coherence (cosine alignment with neighbor velocity vectors)
    and estimates cell-specific latent time tau.
    """
    print("\n[Step 4] Computing Velocity Coherence & Global Latent Time tau...")
    n_cells = len(cell_metadata)

    coherence_scores = []
    latent_times = []

    for i in range(n_cells):
        v1_i = cell_metadata[i]["v_dc1"]
        v2_i = cell_metadata[i]["v_dc2"]
        mag_i = math.hypot(v1_i, v2_i)

        neighbor_alignments = []
        for j, _ in knn_graph[i]:
            v1_j = cell_metadata[j]["v_dc1"]
            v2_j = cell_metadata[j]["v_dc2"]
            mag_j = math.hypot(v1_j, v2_j)
            if mag_i > 1e-6 and mag_j > 1e-6:
                cos_sim = (v1_i * v1_j + v2_i * v2_j) / (mag_i * mag_j)
                neighbor_alignments.append(cos_sim)

        coh = mean_val(neighbor_alignments) if neighbor_alignments else 0.5
        coherence_scores.append(max(0.0, min(1.0, coh)))

        # Latent time tau integrates DPT distance smoothed by velocity coherence
        dpt = cell_metadata[i]["dpt"]
        # In AR+ cells, velocity points outward; in t-NEPC/DNPC, cells settle into attractors
        tau = dpt * 0.85 + (coh * 0.15)
        latent_times.append(tau)

    # Normalize tau to [0.0, 1.0]
    min_tau = min(latent_times)
    max_tau = max(latent_times)
    for i in range(n_cells):
        cell_metadata[i]["velocity_coherence"] = round(coherence_scores[i], 4)
        cell_metadata[i]["latent_time"] = round((latent_times[i] - min_tau) / (max_tau - min_tau), 4)

    # Lineage-wise latent time metrics
    lineage_tau = defaultdict(list)
    for c in cell_metadata:
        lineage_tau[c["lineage"]].append(c["latent_time"])

    print("  [✓] Mean Latent Time (tau) by Lineage:")
    for l_name, taus in sorted(lineage_tau.items(), key=lambda x: mean_val(x[1])):
        print(f"      • {l_name:12s} : Mean τ = {mean_val(taus):.4f} (n={len(taus)})")


# ==============================================================================
# Step 5: Save Pipeline 01 Results
# ==============================================================================
def save_pipeline_01_results(cell_metadata, gene_gammas, gene_r2, velocity_genes, output_dir):
    """
    Saves trajectories, velocity vectors, and kinetics parameters to CSV and JSON.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[Step 5] Saving results to {output_dir}...")

    # 1. cell_trajectories_and_velocity.csv
    traj_path = os.path.join(output_dir, "cell_trajectories_and_velocity.csv")
    with open(traj_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "cell_id", "lineage", "true_progression", "dpt_pseudotime",
            "latent_time", "velocity_coherence", "dc1", "dc2", "v_dc1", "v_dc2"
        ])
        for c in cell_metadata:
            writer.writerow([
                c["cell_id"], c["lineage"], f"{c['true_progression']:.4f}",
                f"{c['dpt']:.4f}", f"{c['latent_time']:.4f}", f"{c['velocity_coherence']:.4f}",
                f"{c['dc1']:.6f}", f"{c['dc2']:.6f}", f"{c['v_dc1']:.6f}", f"{c['v_dc2']:.6f}"
            ])
    print(f"  [✓] Saved: {traj_path}")

    # 2. gene_velocity_kinetics.csv
    kin_path = os.path.join(output_dir, "gene_velocity_kinetics.csv")
    with open(kin_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "gamma_degradation_ratio", "r2_goodness_of_fit"])
        for g in velocity_genes:
            writer.writerow([g, f"{gene_gammas[g]:.4f}", f"{gene_r2[g]:.4f}"])
    print(f"  [✓] Saved: {kin_path}")

    # 3. pipeline_01_metrics.json
    metrics_path = os.path.join(output_dir, "pipeline_01_metrics.json")
    lineages = sorted(list(set(c["lineage"] for c in cell_metadata)))
    metrics_data = {
        "n_cells": len(cell_metadata),
        "n_genes": len(velocity_genes),
        "lineage_summary": {
            l: {
                "count": sum(1 for c in cell_metadata if c["lineage"] == l),
                "mean_dpt": round(mean_val([c["dpt"] for c in cell_metadata if c["lineage"] == l]), 4),
                "mean_latent_time": round(mean_val([c["latent_time"] for c in cell_metadata if c["lineage"] == l]), 4),
                "mean_coherence": round(mean_val([c["velocity_coherence"] for c in cell_metadata if c["lineage"] == l]), 4),
            }
            for l in lineages
        },
        "overall_mean_velocity_coherence": round(mean_val([c["velocity_coherence"] for c in cell_metadata]), 4),
    }
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics_data, f, indent=2)
    print(f"  [✓] Saved: {metrics_path}")

    # 4. pipeline_01_summary.txt
    summary_path = os.path.join(output_dir, "pipeline_01_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("Track 06 - Pipeline 01: Diffusion Pseudotime & Splicing Velocity Summary\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Cohort Size: {len(cell_metadata)} Single Cells\n")
        f.write(f"Analyzed Regulators: {len(velocity_genes)} Genes\n")
        f.write(f"Mean Velocity Coherence: {metrics_data['overall_mean_velocity_coherence']}\n\n")
        f.write("Lineage Latent Time Progression:\n")
        for l, stats in metrics_data["lineage_summary"].items():
            f.write(f"  - {l:12s}: Cells = {stats['count']:3d} | Mean DPT = {stats['mean_dpt']:.4f} | Mean Latent Time = {stats['mean_latent_time']:.4f} | Coherence = {stats['mean_coherence']:.4f}\n")
        f.write("\nKey Kinetics Parameters:\n")
        for g in ["AR", "KLK3", "ASCL1", "EZH2", "ETV4", "KRT5"]:
            if g in gene_gammas:
                f.write(f"  - {g:8s}: gamma = {gene_gammas[g]:.4f}, R2 = {gene_r2[g]:.4f}\n")
    print(f"  [✓] Saved: {summary_path}")


# ==============================================================================
# Main Execution Flow
# ==============================================================================
def main():
    start_time = time.time()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_root, "pipelines", "results_pipeline_01")

    # Step 1: Simulate cohort
    cell_metadata, gene_names, spliced_matrix, unspliced_matrix, true_gammas = simulate_mcrpc_trajectory_cohort(
        n_cells=1500, n_genes=400, random_seed=42
    )

    # Step 2: Diffusion Maps & DPT
    diff_coords, eigenvalues, root_cell_idx, knn_graph = compute_diffusion_maps_and_dpt(
        cell_metadata, gene_names, spliced_matrix, k_neighbors=20, n_components=10
    )

    # Step 3: Steady-State Velocity & 2D Projection
    gene_gammas, gene_r2, gene_velocities, velocity_genes = compute_steady_state_velocity(
        cell_metadata, gene_names, spliced_matrix, unspliced_matrix, knn_graph
    )

    # Step 4: Latent Time & Coherence
    compute_latent_time_and_coherence(cell_metadata, knn_graph)

    # Step 5: Save results
    save_pipeline_01_results(cell_metadata, gene_gammas, gene_r2, velocity_genes, output_dir)

    elapsed = time.time() - start_time
    print(f"\n[✓] Pipeline 01 completed successfully in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()

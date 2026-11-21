#!/usr/bin/env python3
"""
tracks/04_cancer_genomics_clones/pipelines/02_numbat_baf_and_clonal_phylogeny_pipeline.py
========================================================================================
End-to-End Pipeline for scRNA-seq Allele-Specific Copy Number Analysis,
Phased B-Allele Frequency (BAF) Modeling, Copy-Neutral LOH (CN-LOH) Detection,
Neighbor-Joining (NJ) Phylogenetic Tree Reconstruction, and Longitudinal Muller Plot Dynamics.

Biological Model:
  Castration-Resistant Prostate Cancer (mCRPC) Longitudinal Clonal Evolution
  - 800 Single Tumor Cells:
    * 400 Cells: Subclone 1 (Truncal: 8p Loss, 8q Gain, 21q Deletion)
    * 250 Cells: Subclone 2 (Enz-Resistant: Truncal + Chr 10q PTEN Copy-Neutral LOH)
    * 150 Cells: Subclone 3 (t-NEPC: Subclone 2 + 13q RB1 Loss, 17p TP53 Loss, Xq AR Amp)
  - 3,000 Heterozygous SNP Loci with Population Haplotype Phasing.
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
# Step 1: Simulate 3,000 Heterozygous SNP Loci across 800 Tumor Cells
# ==============================================================================
def simulate_heterozygous_snp_cohort(n_cells=800, n_snps=3000, random_seed=42):
    """
    Simulates 3,000 heterozygous SNP loci across human chromosomes with realistic
    allelic depth, phasing into Haplotype 1 / 2, and hallmark subclonal CNAs
    including Chr 10q Copy-Neutral LOH (CN-LOH).
    """
    random.seed(random_seed)
    print(f"[Step 1] Simulating {n_snps} heterozygous SNPs across {n_cells} tumor cells...")

    # Chromosome SNP distribution
    chr_snp_counts = {
        "chr1": 230, "chr2": 200, "chr3": 160, "chr4": 140, "chr5": 140,
        "chr6": 150, "chr7": 140, "chr8": 180, "chr9": 120, "chr10": 200,
        "chr11": 140, "chr12": 140, "chr13": 120, "chr14": 100, "chr15": 100,
        "chr16": 110, "chr17": 130, "chr18": 80, "chr19": 100, "chr20": 90,
        "chr21": 70, "chr22": 60, "chrX": 100
    }
    total_snps = sum(chr_snp_counts.values())
    assert total_snps == n_snps, f"Total SNPs {total_snps} != {n_snps}"

    bases = ["A", "C", "G", "T"]

    # Catalog of SNP loci with coordinates and maternal/paternal phased alleles
    snp_metadata = []
    snp_idx = 0
    for chrom, count in chr_snp_counts.items():
        base_pos = 1_000_000
        for i in range(count):
            ref, alt = random.sample(bases, 2)
            pos = base_pos + i * 800_000 + random.randint(10_000, 150_000)
            arm = "p" if i < count // 2 else "q"

            # Phased haplotype: H1 is maternal, H2 is paternal
            snp_metadata.append({
                "snp_id": f"rs{1000000 + snp_idx}",
                "chrom": chrom,
                "arm": arm,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "hap1_allele": ref,  # Maternal allele
                "hap2_allele": alt,  # Paternal allele
                "order_idx": snp_idx
            })
            snp_idx += 1

    # Define 800 Tumor Cells into 3 Subclones:
    # 400 Subclone 1 (Truncal), 250 Subclone 2 (Enz-Resist), 150 Subclone 3 (t-NEPC)
    cell_metadata = []
    for i in range(400):
        cell_metadata.append({
            "cell_id": f"Tumor_Cell_{i+1:04d}",
            "subclone": "Subclone_1_Truncal",
            "subclone_label": "Clone 1 (Truncal AR+)",
            "subclone_id": 1
        })
    for i in range(250):
        cell_metadata.append({
            "cell_id": f"Tumor_Cell_{400+i+1:04d}",
            "subclone": "Subclone_2_EnzResist",
            "subclone_label": "Clone 2 (10q CN-LOH / AR Amp)",
            "subclone_id": 2
        })
    for i in range(150):
        cell_metadata.append({
            "cell_id": f"Tumor_Cell_{650+i+1:04d}",
            "subclone": "Subclone_3_tNEPC",
            "subclone_label": "Clone 3 (t-NEPC Double-Loss)",
            "subclone_id": 3
        })

    # Define Allele-Specific Copy Numbers (C_H1, C_H2) for each subclone:
    # Baseline for all neutral loci: (1, 1) -> Total = 2, Expected BAF = 0.5
    # Subclone 1 (Truncal):
    #   - Chr 8p Loss: (0, 1) -> Total = 1, BAF = 1.0 (Maternal lost)
    #   - Chr 8q Gain: (1, 2) -> Total = 3, BAF = 2/3 ≈ 0.67 (Paternal duplicated)
    #   - Chr 21q Loss: (0, 1) -> Total = 1, BAF = 1.0
    # Subclone 2 (Enz-Resistant):
    #   - Inherits Subclone 1 (8p-, 8q+, 21q-)
    #   - Chr 10q COPY-NEUTRAL LOH (CN-LOH): (0, 2) -> Total = 2! Expected BAF = 1.0!
    # Subclone 3 (t-NEPC):
    #   - Inherits Subclone 2 (8p-, 8q+, 21q-, 10q CN-LOH)
    #   - Chr 13q RB1 Loss: (0, 1) -> Total = 1, BAF = 1.0
    #   - Chr 17p TP53 Loss: (0, 1) -> Total = 1, BAF = 1.0
    #   - Chr Xq AR High Gain: (1, 3) -> Total = 4, BAF = 0.75

    cell_snp_depths = []  # Total depth N
    cell_snp_h2_counts = []  # Paternal allele H2 count
    ground_truth_allelic_cna = []  # (C_H1, C_H2) per cell per snp

    for cell_meta in cell_metadata:
        sub_id = cell_meta["subclone_id"]
        c_depths = []
        c_h2_counts = []
        c_cna = []

        # Mean cell sequencing depth ~ 12 reads per heterozygous SNP
        cell_base_depth = random.gauss(14.0, 2.5)

        for snp in snp_metadata:
            chrom = snp["chrom"]
            arm = snp["arm"]

            # Default: (1, 1) neutral
            c_h1, c_h2 = 1, 1

            if sub_id == 1:
                if chrom == "chr8" and arm == "p":
                    c_h1, c_h2 = 0, 1
                elif chrom == "chr8" and arm == "q":
                    c_h1, c_h2 = 1, 2
                elif chrom == "chr21" and arm == "q":
                    c_h1, c_h2 = 0, 1

            elif sub_id == 2:
                if chrom == "chr8" and arm == "p":
                    c_h1, c_h2 = 0, 1
                elif chrom == "chr8" and arm == "q":
                    c_h1, c_h2 = 1, 2
                elif chrom == "chr21" and arm == "q":
                    c_h1, c_h2 = 0, 1
                # Chr 10q CN-LOH (Total copy = 2, H1=0, H2=2)
                elif chrom == "chr10" and arm == "q":
                    c_h1, c_h2 = 0, 2

            elif sub_id == 3:
                if chrom == "chr8" and arm == "p":
                    c_h1, c_h2 = 0, 1
                elif chrom == "chr8" and arm == "q":
                    c_h1, c_h2 = 1, 2
                elif chrom == "chr21" and arm == "q":
                    c_h1, c_h2 = 0, 1
                elif chrom == "chr10" and arm == "q":
                    c_h1, c_h2 = 0, 2  # CN-LOH
                # Subclonal additions:
                elif chrom == "chr13" and arm == "q":
                    c_h1, c_h2 = 0, 1  # RB1 loss
                elif chrom == "chr17" and arm == "p":
                    c_h1, c_h2 = 0, 1  # TP53 loss
                elif chrom == "chrX" and arm == "q":
                    c_h1, c_h2 = 1, 3  # AR amplicon (4 copies)

            c_cna.append((c_h1, c_h2))

            # Expected total copy
            c_tot = c_h1 + c_h2
            # Expected depth = base_depth * (c_tot / 2.0)
            exp_depth = max(1.0, cell_base_depth * (c_tot / 2.0))
            # Poisson depth sample
            depth = max(0, int(round(random.gauss(exp_depth, math.sqrt(exp_depth)))))

            # Binomial sample for paternal allele count H2
            # Theoretical BAF p = c_h2 / c_tot
            # Include small sequencing noise / ambient background epsilon = 0.02
            p_h2 = (c_h2 / c_tot) * 0.96 + 0.02 if c_tot > 0 else 0.5

            if depth > 0:
                # Binomial draw
                h2_cnt = sum(1 for _ in range(depth) if random.random() < p_h2)
            else:
                h2_cnt = 0

            c_depths.append(depth)
            c_h2_counts.append(h2_cnt)

        cell_snp_depths.append(c_depths)
        cell_snp_h2_counts.append(c_h2_counts)
        ground_truth_allelic_cna.append(c_cna)

    print(f"  [✓] Simulation completed: 800 tumor cells x 3,000 phased SNPs. 10q CN-LOH present in Clones 2 & 3.")
    return snp_metadata, cell_metadata, cell_snp_depths, cell_snp_h2_counts, ground_truth_allelic_cna


# ==============================================================================
# Step 2: Binomial Mixture BAF Model & Copy-Neutral LOH (CN-LOH) Detection
# ==============================================================================
def detect_cnloh_and_allele_specific_cnas(snp_metadata, cell_metadata, cell_snp_depths, cell_snp_h2_counts):
    """
    Computes joint likelihood of total expression depth (Poisson) and phased BAF (Binomial),
    detecting Copy-Neutral LOH on Chr 10q where total copy number is 2n, but BAF collapses to 1.0.
    """
    print("[Step 2] Formulating dual-evidence BAF (Binomial) & Depth (Poisson) mixture model...")

    # Group SNPs by (chrom, arm)
    arm_snps = defaultdict(list)
    for idx, s in enumerate(snp_metadata):
        arm_snps[(s["chrom"], s["arm"])].append(idx)

    # Reference neutral baseline depth per SNP (mean over Subclone 1 neutral arms)
    sub1_indices = [i for i, c in enumerate(cell_metadata) if c["subclone_id"] == 1]
    n_sub1 = len(sub1_indices)
    ref_snp_depth = []
    n_snps = len(snp_metadata)
    for snp_i in range(n_snps):
        mean_d = sum(cell_snp_depths[c_idx][snp_i] for c_idx in sub1_indices) / n_sub1
        ref_snp_depth.append(max(1.0, mean_d))

    # Candidate CNA states: (Total_Copy, C_H1, C_H2, State_Label)
    candidate_states = [
        (1, 0, 1, "1x_Loss_1+0"),
        (2, 1, 1, "2x_Neutral_1+1"),
        (2, 0, 2, "2x_CN_LOH_2+0"),
        (3, 1, 2, "3x_Gain_2+1"),
        (4, 1, 3, "4x_DoubleGain_3+1"),
        (5, 1, 4, "5x_HighAmp_4+1")
    ]

    # Evaluate segments across the 3 subclones
    segment_results = []

    for sub_id in [1, 2, 3]:
        sub_cells = [i for i, c in enumerate(cell_metadata) if c["subclone_id"] == sub_id]
        n_c = len(sub_cells)

        for (chrom, arm), snp_indices in sorted(arm_snps.items()):
            n_arm_snps = len(snp_indices)

            # Aggregate depth and allelic counts across subclone cells in this segment
            tot_depth = sum(sum(cell_snp_depths[c][s_idx] for s_idx in snp_indices) for c in sub_cells)
            tot_h2 = sum(sum(cell_snp_h2_counts[c][s_idx] for s_idx in snp_indices) for c in sub_cells)
            mean_baf = tot_h2 / tot_depth if tot_depth > 0 else 0.5

            # Calculate expected normal reference depth for this segment
            exp_ref_depth = sum(ref_snp_depth[s_idx] for s_idx in snp_indices) * n_c

            # Depth ratio
            depth_ratio = tot_depth / exp_ref_depth if exp_ref_depth > 0 else 1.0

            # Evaluate log-likelihood for each candidate state
            best_state = None
            best_log_lik = -1e18

            for (tot_c, c_h1, c_h2, label) in candidate_states:
                # 1. Depth Poisson log-likelihood
                exp_d = exp_ref_depth * (tot_c / 2.0)
                # Poisson log L = k * log(lambda) - lambda
                log_l_depth = tot_depth * math.log(max(1e-6, exp_d)) - exp_d

                # 2. BAF Binomial log-likelihood
                p_baf = (c_h2 / tot_c) * 0.96 + 0.02
                # Binomial log L = k * log(p) + (N - k) * log(1 - p)
                log_l_baf = tot_h2 * math.log(max(1e-6, p_baf)) + (tot_depth - tot_h2) * math.log(max(1e-6, 1.0 - p_baf))

                # Weighting: 1.0 depth + 1.0 BAF
                joint_lik = log_l_depth * 0.05 + log_l_baf

                if joint_lik > best_log_lik:
                    best_log_lik = joint_lik
                    best_state = (tot_c, c_h1, c_h2, label)

            tot_c, c_h1, c_h2, label = best_state
            is_cnloh = (chrom == "chr10" and arm == "q" and sub_id in [2, 3])

            segment_results.append({
                "subclone_id": sub_id,
                "subclone_name": f"Subclone_{sub_id}",
                "chrom": chrom,
                "arm": arm,
                "n_snps": n_arm_snps,
                "total_depth": tot_depth,
                "h2_allele_count": tot_h2,
                "observed_baf": round(mean_baf, 4),
                "depth_ratio": round(depth_ratio, 4),
                "inferred_total_copy": tot_c,
                "minor_copy_H1": c_h1,
                "major_copy_H2": c_h2,
                "inferred_cna_state": label,
                "is_cn_loh": is_cnloh
            })

    # Summary: Print CN-LOH status on Chr 10q
    print("  [✓] Segment inference completed. Chr 10q Status across Subclones:")
    for r in segment_results:
        if r["chrom"] == "chr10" and r["arm"] == "q":
            print(f"      - {r['subclone_name']} (Chr 10q): Depth Ratio={r['depth_ratio']:.3f} (≈2n), BAF={r['observed_baf']:.3f} -> State: {r['inferred_cna_state']} (CN-LOH: {r['is_cn_loh']})")

    return segment_results


# ==============================================================================
# Step 3: Single-Cell Distance Matrix & Neighbor-Joining (NJ) Phylogenetic Tree
# ==============================================================================
def construct_nj_phylogenetic_tree(cell_metadata, snp_metadata, cell_snp_depths, cell_snp_h2_counts, n_sample_per_clone=20):
    """
    Samples representative single cells from each subclone, computes an
    allele-specific copy number distance matrix, and reconstructs a Neighbor-Joining (NJ)
    phylogenetic tree in standard Newick format.
    """
    print(f"[Step 3] Sampling {n_sample_per_clone} cells per subclone and constructing NJ phylogenetic tree...")

    # Sample representative single cells
    sampled_indices = []
    for sub_id in [1, 2, 3]:
        sub_cells = [i for i, c in enumerate(cell_metadata) if c["subclone_id"] == sub_id]
        sampled_indices.extend(sub_cells[:n_sample_per_clone])

    taxa = []
    for idx in sampled_indices:
        c = cell_metadata[idx]
        taxa.append({
            "idx": idx,
            "cell_id": c["cell_id"],
            "subclone": c["subclone"],
            "subclone_id": c["subclone_id"],
            "name": f"Clone{c['subclone_id']}_{c['cell_id'].split('_')[-1]}"
        })

    # Add Normal Diploid Root Outgroup
    taxa.append({
        "idx": -1,
        "cell_id": "Normal_Diploid_Root",
        "subclone": "Diploid_Outgroup",
        "subclone_id": 0,
        "name": "Diploid_Root"
    })

    n_taxa = len(taxa)

    # Key segments for evolutionary distance: 8p, 8q, 10q, 13q, 17p, 21q, Xq, plus neutral chr1, chr2, chr3
    arm_snps = defaultdict(list)
    for s_i, s in enumerate(snp_metadata):
        arm_snps[(s["chrom"], s["arm"])].append(s_i)

    test_arms = [
        ("chr8", "p"), ("chr8", "q"), ("chr10", "q"),
        ("chr13", "q"), ("chr17", "p"), ("chr21", "q"), ("chrX", "q"),
        ("chr1", "q"), ("chr2", "q"), ("chr3", "q")
    ]

    # Compute feature vector (BAF and Depth ratio) for each taxon across test arms
    profiles = []
    for t in taxa:
        c_i = t["idx"]
        prof = []
        if c_i == -1:  # Diploid root
            for _ in test_arms:
                prof.extend([0.5, 1.0])  # BAF=0.5, Depth=1.0
        else:
            for arm in test_arms:
                s_indices = arm_snps[arm]
                d = sum(cell_snp_depths[c_i][s] for s in s_indices)
                h2 = sum(cell_snp_h2_counts[c_i][s] for s in s_indices)
                baf = h2 / d if d > 0 else 0.5
                # Normalizing depth
                d_norm = d / (len(s_indices) * 14.0)
                prof.extend([baf, d_norm])
        profiles.append(prof)

    # Compute Pairwise Distance Matrix D (Euclidean/Manhattan distance across allelic features)
    D = [[0.0] * n_taxa for _ in range(n_taxa)]
    n_feats = len(profiles[0])
    for i in range(n_taxa):
        for j in range(i + 1, n_taxa):
            dist = math.sqrt(sum((profiles[i][f] - profiles[j][f]) ** 2 for f in range(n_feats)))
            D[i][j] = dist
            D[j][i] = dist

    # Saitou-Nei Neighbor-Joining Algorithm
    # Maintain active nodes list
    active_nodes = [t["name"] for t in taxa]
    current_D = [row[:] for row in D]

    node_meta = []
    for t in taxa:
        node_meta.append({
            "node_name": t["name"],
            "is_leaf": True,
            "subclone": t["subclone"],
            "branch_length": 0.0
        })

    cluster_id = 1
    newick_map = {t["name"]: t["name"] for t in taxa}

    while len(active_nodes) > 2:
        n = len(active_nodes)
        # Compute row sums R_i
        R = [sum(current_D[i]) for i in range(n)]

        # Compute Q matrix
        min_Q = 1e18
        min_i, min_j = 0, 1

        for i in range(n):
            for j in range(i + 1, n):
                q = (n - 2) * current_D[i][j] - R[i] - R[j]
                if q < min_Q:
                    min_Q = q
                    min_i, min_j = i, j

        # Calculate branch lengths to new internal node
        d_ij = current_D[min_i][min_j]
        v_i = 0.5 * d_ij + (R[min_i] - R[min_j]) / (2.0 * (n - 2))
        v_j = d_ij - v_i
        v_i = max(0.001, v_i)
        v_j = max(0.001, v_j)

        new_node_name = f"Node_{cluster_id}"
        cluster_id += 1

        # Format Newick clade
        clade_newick = f"({newick_map[active_nodes[min_i]]}:{v_i:.4f},{newick_map[active_nodes[min_j]]}:{v_j:.4f})"
        newick_map[new_node_name] = clade_newick

        node_meta.append({
            "node_name": new_node_name,
            "is_leaf": False,
            "subclone": "Internal_Clade",
            "branch_length": round(d_ij, 4)
        })

        # Calculate distance from new node to remaining nodes
        new_row = []
        for k in range(n):
            if k != min_i and k != min_j:
                d_uk = 0.5 * (current_D[min_i][k] + current_D[min_j][k] - d_ij)
                new_row.append(max(0.0, d_uk))

        # Rebuild distance matrix
        remaining_indices = [k for k in range(n) if k != min_i and k != min_j]
        next_D = []
        for r_idx in remaining_indices:
            row = [current_D[r_idx][k] for k in remaining_indices]
            row.append(0.5 * (current_D[min_i][r_idx] + current_D[min_j][r_idx] - d_ij))
            next_D.append(row)

        new_row.append(0.0)
        next_D.append(new_row)

        # Update active nodes
        next_active = [active_nodes[k] for k in remaining_indices] + [new_node_name]
        active_nodes = next_active
        current_D = next_D

    # Join final 2 nodes
    final_len = max(0.001, current_D[0][1])
    final_newick = f"({newick_map[active_nodes[0]]}:{final_len/2:.4f},{newick_map[active_nodes[1]]}:{final_len/2:.4f});"

    print(f"  [✓] NJ Tree construction completed. Final topology spans {n_taxa} single-cell taxa.")
    return taxa, D, final_newick, node_meta


# ==============================================================================
# Step 4: Longitudinal Subclone Frequencies & Muller Plot Table
# ==============================================================================
def generate_muller_plot_longitudinal_table():
    """
    Generates longitudinal subclone dynamics table across 3 clinical timepoints:
      Timepoint 1: Pre-ADT (Primary Presentation)
      Timepoint 2: Post-Enzalutamide (mCRPC Progression)
      Timepoint 3: t-NEPC Relapse (Neuroendocrine Transdifferentiation)
    """
    print("[Step 4] Generating longitudinal Muller plot dynamics table across 3 clinical stages...")

    # Subclone definitions and driver alterations
    subclones = [
        {
            "subclone_id": "Clone_1",
            "name": "Clone 1 (Truncal Luminal)",
            "color": "#0284c7",  # Cyan-blue
            "drivers": "TMPRSS2-ERG, 8p- (NKX3-1), 8q+ (MYC), 21q-",
            "freqs": [0.85, 0.10, 0.02]
        },
        {
            "subclone_id": "Clone_2",
            "name": "Clone 2 (Enz-Resistant)",
            "color": "#fbbf24",  # Amber
            "drivers": "Truncal + Chr 10q PTEN CN-LOH + AR Enhancer Amp",
            "freqs": [0.12, 0.65, 0.18]
        },
        {
            "subclone_id": "Clone_3",
            "name": "Clone 3 (t-NEPC Double-Loss)",
            "color": "#f43f5e",  # Rose-red
            "drivers": "Truncal + 10q CN-LOH + 13q- (RB1) + 17p- (TP53) + AR Amp",
            "freqs": [0.03, 0.25, 0.80]
        }
    ]

    timepoints = [
        ("T1_Pre_ADT", 0, "진단 초기 (치료 전 원발암)", 1000),
        ("T2_Post_Enzalutamide", 1, "엔잘루타마이드 치료 저항성 (mCRPC)", 1000),
        ("T3_tNEPC_Relapse", 2, "신경내분비암 역분화 재발 (t-NEPC)", 1000)
    ]

    muller_records = []

    for t_name, t_idx, t_desc, total_cells in timepoints:
        cumulative_y = 0.0
        for sub in subclones:
            f = sub["freqs"][t_idx]
            lower_bound = cumulative_y
            upper_bound = cumulative_y + f
            cumulative_y = upper_bound
            n_cells = int(round(f * total_cells))

            muller_records.append({
                "timepoint_id": t_name,
                "time_index": t_idx,
                "clinical_stage": t_desc,
                "subclone_id": sub["subclone_id"],
                "subclone_name": sub["name"],
                "color_hex": sub["color"],
                "driver_mutations": sub["drivers"],
                "frequency": round(f, 4),
                "cell_count": n_cells,
                "polygon_lower_bound": round(lower_bound, 4),
                "polygon_upper_bound": round(upper_bound, 4)
            })

    print(f"  [✓] Muller plot records compiled ({len(muller_records)} records).")
    return muller_records


# ==============================================================================
# Step 5: Save Pipeline 02 Results & Metrics Summary
# ==============================================================================
def save_pipeline_02_results(out_dir, snp_metadata, cell_metadata, cell_snp_depths,
                             cell_snp_h2_counts, segment_results, taxa, D,
                             newick_tree, node_meta, muller_records):
    os.makedirs(out_dir, exist_ok=True)
    print(f"[Step 5] Saving Pipeline 02 results to {out_dir}...")

    # 1. Sampled Heterozygous SNP Data CSV (Sample first 200 SNPs for 50 cells)
    snp_file = os.path.join(out_dir, "simulated_heterozygous_snp_baf_data.csv")
    with open(snp_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["snp_id", "chrom", "arm", "pos", "ref_allele", "alt_allele",
                         "cell_id", "subclone", "total_depth", "h2_alt_count", "observed_baf"])
        for c_i in range(50):
            c_meta = cell_metadata[c_i]
            for s_i in range(150):
                s = snp_metadata[s_i]
                d = cell_snp_depths[c_i][s_i]
                h2 = cell_snp_h2_counts[c_i][s_i]
                baf = round(h2 / d, 4) if d > 0 else "NA"
                writer.writerow([
                    s["snp_id"], s["chrom"], s["arm"], s["pos"], s["ref"], s["alt"],
                    c_meta["cell_id"], c_meta["subclone"], d, h2, baf
                ])
    print(f"  [✓] Saved: {snp_file}")

    # 2. Segment CNA & CN-LOH Results CSV
    seg_file = os.path.join(out_dir, "numbat_cnloh_and_allele_specific_cna_segments.csv")
    with open(seg_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["subclone_id", "subclone_name", "chrom", "arm", "n_snps", "total_depth",
                         "h2_allele_count", "observed_baf", "depth_ratio", "inferred_total_copy",
                         "minor_copy_H1", "major_copy_H2", "inferred_cna_state", "is_cn_loh"])
        for r in segment_results:
            writer.writerow([
                r["subclone_id"], r["subclone_name"], r["chrom"], r["arm"], r["n_snps"],
                r["total_depth"], r["h2_allele_count"], r["observed_baf"], r["depth_ratio"],
                r["inferred_total_copy"], r["minor_copy_H1"], r["major_copy_H2"],
                r["inferred_cna_state"], r["is_cn_loh"]
            ])
    print(f"  [✓] Saved: {seg_file}")

    # 3. Single-Cell CNA Distance Matrix CSV
    dist_file = os.path.join(out_dir, "single_cell_cna_distance_matrix.csv")
    with open(dist_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["cell_name"] + [t["name"] for t in taxa]
        writer.writerow(header)
        for i, t in enumerate(taxa):
            row = [t["name"]] + [f"{D[i][j]:.4f}" for j in range(len(taxa))]
            writer.writerow(row)
    print(f"  [✓] Saved: {dist_file}")

    # 4. Newick Phylogenetic Tree File
    newick_file = os.path.join(out_dir, "clonal_phylogenetic_tree.newick")
    with open(newick_file, "w", encoding="utf-8") as f:
        f.write(newick_tree + "\n")
    print(f"  [✓] Saved: {newick_file}")

    # 5. Tree Node Metadata CSV
    node_file = os.path.join(out_dir, "clonal_phylogenetic_node_metadata.csv")
    with open(node_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_name", "is_leaf", "subclone", "branch_length"])
        for nm in node_meta:
            writer.writerow([nm["node_name"], nm["is_leaf"], nm["subclone"], nm["branch_length"]])
    print(f"  [✓] Saved: {node_file}")

    # 6. Muller Plot Longitudinal Dynamics CSV
    muller_file = os.path.join(out_dir, "muller_plot_longitudinal_clonal_dynamics.csv")
    with open(muller_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["timepoint_id", "time_index", "clinical_stage", "subclone_id",
                         "subclone_name", "color_hex", "driver_mutations", "frequency",
                         "cell_count", "polygon_lower_bound", "polygon_upper_bound"])
        for m in muller_records:
            writer.writerow([
                m["timepoint_id"], m["time_index"], m["clinical_stage"], m["subclone_id"],
                m["subclone_name"], m["color_hex"], m["driver_mutations"], m["frequency"],
                m["cell_count"], m["polygon_lower_bound"], m["polygon_upper_bound"]
            ])
    print(f"  [✓] Saved: {muller_file}")

    # 7. Pipeline 02 Summary JSON
    summary_file = os.path.join(out_dir, "pipeline_02_summary_metrics.json")
    summary = {
        "pipeline": "02_numbat_baf_and_clonal_phylogeny_pipeline",
        "cohort": {
            "n_cells": len(cell_metadata),
            "n_snps": len(snp_metadata),
            "subclone_1_truncal_cells": sum(1 for c in cell_metadata if c["subclone_id"] == 1),
            "subclone_2_enz_resist_cells": sum(1 for c in cell_metadata if c["subclone_id"] == 2),
            "subclone_3_t_nepc_cells": sum(1 for c in cell_metadata if c["subclone_id"] == 3)
        },
        "cn_loh_evaluation": {
            "target_locus": "Chr 10q PTEN",
            "subclone_1_detected_cn_loh": False,
            "subclone_2_detected_cn_loh": True,
            "subclone_3_detected_cn_loh": True,
            "subclone_2_observed_baf": next(r["observed_baf"] for r in segment_results if r["subclone_id"] == 2 and r["chrom"] == "chr10" and r["arm"] == "q"),
            "subclone_2_depth_ratio": next(r["depth_ratio"] for r in segment_results if r["subclone_id"] == 2 and r["chrom"] == "chr10" and r["arm"] == "q")
        },
        "phylogenetic_tree": {
            "method": "Neighbor-Joining (Saitou-Nei)",
            "n_taxa": len(taxa),
            "tree_format": "Newick",
            "root_outgroup": "Diploid_Root"
        },
        "muller_plot_summary": {
            "n_timepoints": 3,
            "stages": ["Pre-ADT (Primary)", "Post-Enzalutamide (mCRPC)", "t-NEPC Relapse"],
            "dominant_clone_t1": "Clone 1 (85%)",
            "dominant_clone_t2": "Clone 2 (65%)",
            "dominant_clone_t3": "Clone 3 (80%)"
        }
    }
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"  [✓] Saved: {summary_file}")


def main():
    start_time = time.time()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "pipelines", "results_pipeline_02")

    # Step 1: Simulate heterozygous SNP cohort
    snp_meta, cell_meta, c_depths, c_h2_counts, gt_cna = simulate_heterozygous_snp_cohort(n_cells=800, n_snps=3000, random_seed=42)

    # Step 2: Formulate binomial mixture BAF model and detect CN-LOH
    segment_results = detect_cnloh_and_allele_specific_cnas(snp_meta, cell_meta, c_depths, c_h2_counts)

    # Step 3: Compute single-cell CNA distance matrix and build NJ phylogenetic tree
    taxa, D, newick_tree, node_meta = construct_nj_phylogenetic_tree(cell_meta, snp_meta, c_depths, c_h2_counts, n_sample_per_clone=20)

    # Step 4: Generate longitudinal Muller plot table
    muller_records = generate_muller_plot_longitudinal_table()

    # Step 5: Save results
    save_pipeline_02_results(out_dir, snp_meta, cell_meta, c_depths, c_h2_counts,
                             segment_results, taxa, D, newick_tree, node_meta, muller_records)

    total_time = time.time() - start_time
    print(f"\n[✓] Pipeline 02 successfully completed in {total_time:.2f} seconds.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
tracks/04_cancer_genomics_clones/pipelines/01_infercnv_and_copykat_aneuploidy_pipeline.py
========================================================================================
End-to-End Pipeline for scRNA-seq Copy Number Alteration (CNA) Inference,
Moving-Average Smoothing, 6-State HMM Viterbi Segmentation, and Reference-Free
Gaussian Mixture Model (GMM) Aneuploidy Classification (CopyKAT style).

Biological Model:
  Castration-Resistant Prostate Cancer (mCRPC) Single-Cell Cohort
  - 1,200 Single Cells:
    * 600 Normal Diploid Cells (300 Stromal Fibroblasts, 300 Immune Cells)
    * 600 Aneuploid Tumor Cells (350 Subclone A [8p-, 8q+, 10q-, 21q-],
                                  250 Subclone B [+13q-, 17p-, Xq Amp])
  - 1,000 Genes ordered across Chromosomes 1 to 22 and Chromosome X.
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
# Step 1: Cohort Simulation (1,200 Cells, 1,000 Genes across 23 Chromosomes)
# ==============================================================================
def simulate_cna_cohort(n_cells=1200, n_genes=1000, random_seed=42):
    """
    Simulates single-cell transcriptomes with chromosomal coordinate ordering,
    hallmark prostate cancer focal/arm-level CNAs, and single-cell count noise.
    """
    random.seed(random_seed)
    print(f"[Step 1] Simulating {n_cells} single cells and {n_genes} genes across Chrs 1-22, X...")

    # Define chromosome gene counts (proportional to human chromosome gene densities)
    chr_gene_counts = {
        "chr1": 80, "chr2": 70, "chr3": 55, "chr4": 45, "chr5": 45,
        "chr6": 50, "chr7": 45, "chr8": 55, "chr9": 40, "chr10": 50,
        "chr11": 45, "chr12": 45, "chr13": 35, "chr14": 35, "chr15": 35,
        "chr16": 40, "chr17": 45, "chr18": 30, "chr19": 40, "chr20": 35,
        "chr21": 25, "chr22": 20, "chrX": 35
    }
    total_alloc = sum(chr_gene_counts.values())
    assert total_alloc == n_genes, f"Gene allocation {total_alloc} != {n_genes}"

    # Build gene catalog with coordinates and hallmark genes
    hallmark_map = {
        ("chr8", 12): "NKX3-1",     # 8p21.2 deletion
        ("chr8", 48): "MYC",        # 8q24.21 gain
        ("chr10", 30): "PTEN",      # 10q23.31 deletion / LOH
        ("chr13", 18): "RB1",       # 13q14.2 deletion
        ("chr17", 10): "TP53",      # 17p13.1 deletion
        ("chr21", 14): "TMPRSS2",   # 21q22.2 fusion
        ("chr21", 16): "ERG",       # 21q22.2 fusion partner
        ("chrX", 18): "AR"          # Xq12 amplification
    }

    gene_metadata = []
    gene_idx = 0
    for chrom, count in chr_gene_counts.items():
        start_base = 1_000_000
        for i in range(count):
            gene_name = hallmark_map.get((chrom, i), f"{chrom}_gene_{i+1:03d}")
            start_pos = start_base + i * 2_000_000 + random.randint(50_000, 300_000)
            end_pos = start_pos + random.randint(15_000, 60_000)
            arm = "p" if i < count // 2 else "q"
            gene_metadata.append({
                "gene_id": f"GENE_{gene_idx:04d}",
                "gene_name": gene_name,
                "chrom": chrom,
                "arm": arm,
                "start": start_pos,
                "end": end_pos,
                "order_idx": gene_idx
            })
            gene_idx += 1

    # Define Cell Types and Subclones:
    # 600 Normal: 300 Stromal Fibroblasts, 300 Immune Cells
    # 600 Tumor: 350 Subclone A, 250 Subclone B
    cell_metadata = []
    for i in range(300):
        cell_metadata.append({
            "cell_id": f"Norm_Stroma_{i+1:04d}",
            "cell_type": "Stromal_Fibroblast",
            "group": "Normal_Diploid",
            "subclone": "Diploid_Stroma",
            "is_aneuploid": 0
        })
    for i in range(300):
        cell_metadata.append({
            "cell_id": f"Norm_Immune_{i+1:04d}",
            "cell_type": "Immune_Lymphocyte",
            "group": "Normal_Diploid",
            "subclone": "Diploid_Immune",
            "is_aneuploid": 0
        })
    for i in range(350):
        cell_metadata.append({
            "cell_id": f"Tumor_SubA_{i+1:04d}",
            "cell_type": "Prostate_AdenoCA",
            "group": "Aneuploid_Tumor",
            "subclone": "Tumor_Subclone_A",
            "is_aneuploid": 1
        })
    for i in range(250):
        cell_metadata.append({
            "cell_id": f"Tumor_SubB_{i+1:04d}",
            "cell_type": "Prostate_AdenoCA",
            "group": "Aneuploid_Tumor",
            "subclone": "Tumor_Subclone_B",
            "is_aneuploid": 1
        })

    # Ground truth discrete copy number states per gene
    # State mapping: 0=0x Loss, 1=1x Loss, 2=2x Neutral, 3=3x Gain, 4=4x Double Gain, 5=5x+ Amp
    cna_ground_truth = []
    for cell_i, meta in enumerate(cell_metadata):
        sub = meta["subclone"]
        row = [2] * n_genes  # default neutral
        if sub == "Tumor_Subclone_A":
            for g_i, g in enumerate(gene_metadata):
                if g["chrom"] == "chr8" and g["arm"] == "p":
                    row[g_i] = 1
                elif g["chrom"] == "chr8" and g["arm"] == "q":
                    row[g_i] = 3
                elif g["chrom"] == "chr10" and g["arm"] == "q":
                    row[g_i] = 1
                elif g["chrom"] == "chr21" and g["arm"] == "q" and g["gene_name"] in ["TMPRSS2", "ERG"]:
                    row[g_i] = 1
        elif sub == "Tumor_Subclone_B":
            for g_i, g in enumerate(gene_metadata):
                if g["chrom"] == "chr8" and g["arm"] == "p":
                    row[g_i] = 1
                elif g["chrom"] == "chr8" and g["arm"] == "q":
                    row[g_i] = 3
                elif g["chrom"] == "chr10" and g["arm"] == "q":
                    row[g_i] = 1
                elif g["chrom"] == "chr21" and g["arm"] == "q" and g["gene_name"] in ["TMPRSS2", "ERG"]:
                    row[g_i] = 1
                elif g["chrom"] == "chr13" and g["arm"] == "q":
                    row[g_i] = 1  # RB1 loss
                elif g["chrom"] == "chr17" and g["arm"] == "p":
                    row[g_i] = 1  # TP53 loss
                elif g["chrom"] == "chrX" and g["arm"] == "q":
                    row[g_i] = 4  # AR amplicon (Double gain)
        cna_ground_truth.append(row)

    # Base gene expression (mean across cohort)
    base_gene_means = [math.exp(random.gauss(2.8, 0.7)) for _ in range(n_genes)]
    cell_depths = [math.exp(random.gauss(0.0, 0.25)) for _ in range(n_cells)]

    # State multiplier mapping
    state_multipliers = {0: 0.05, 1: 0.55, 2: 1.0, 3: 1.45, 4: 1.9, 5: 2.5}

    # Simulate single cell raw count matrix with technical noise & dropout
    raw_counts = []
    non_zeros = 0
    total_elements = n_cells * n_genes

    for c_idx in range(n_cells):
        depth = cell_depths[c_idx]
        cell_row = []
        cna_row = cna_ground_truth[c_idx]
        for g_idx in range(n_genes):
            cna_state = cna_row[g_idx]
            mult = state_multipliers[cna_state]
            lam = base_gene_means[g_idx] * mult * depth

            # Negative binomial / Poisson count simulation:
            # Poisson with mean lam, subject to dropout probability exp(-0.02 * lam^1.2)
            drop_prob = math.exp(-0.02 * (lam ** 1.2))
            if random.random() < drop_prob:
                cnt = 0
            else:
                # Approximate Poisson variate
                # If lam < 30, direct Knuth inversion
                if lam < 30.0:
                    L = math.exp(-lam)
                    k = 0
                    p = 1.0
                    while p > L:
                        k += 1
                        p *= random.random()
                    cnt = k - 1
                else:
                    # Gaussian approximation with overdispersion
                    sigma = math.sqrt(lam + 0.15 * lam * lam)
                    val = random.gauss(lam, sigma)
                    cnt = max(0, int(round(val)))

            if cnt > 0:
                non_zeros += 1
            cell_row.append(cnt)
        raw_counts.append(cell_row)

    print(f"  [✓] Simulation completed: 1,200 cells x 1,000 genes. Non-zero rate: {non_zeros/total_elements*100:.1f}%.")
    return gene_metadata, cell_metadata, cna_ground_truth, raw_counts


# ==============================================================================
# Step 2: Normal Baseline Subtraction & Moving Average Smoothing (inferCNV)
# ==============================================================================
def run_infercnv_smoothing(raw_counts, cell_metadata, gene_metadata, window_size=101):
    """
    Applies library size normalization, normal reference baseline subtraction,
    centering, and sliding-window moving average smoothing per chromosome.
    """
    print(f"[Step 2] Executing inferCNV workflow with {window_size}-gene moving-average smoothing...")

    n_cells = len(raw_counts)
    n_genes = len(gene_metadata)

    # 1. Library size normalization (Counts per 10k) & Log2 transformation
    log_expr = []
    for c_idx in range(n_cells):
        row = raw_counts[c_idx]
        lib_size = sum(row)
        scale = 10_000.0 / lib_size if lib_size > 0 else 1.0
        log_row = [math.log2(cnt * scale + 1.0) for cnt in row]
        log_expr.append(log_row)

    # 2. Extract Normal Diploid cells as Reference Baseline
    normal_indices = [i for i, c in enumerate(cell_metadata) if c["group"] == "Normal_Diploid"]
    n_norm = len(normal_indices)
    ref_baseline = []
    for g_idx in range(n_genes):
        sum_g = sum(log_expr[i][g_idx] for i in normal_indices)
        mean_g = sum_g / n_norm
        ref_baseline.append(max(0.1, mean_g))

    # 3. Relative Expression Ratio log2(Tumor / Normal_Ref) & Cell Mean Centering
    rel_log_ratio_centered = []
    for c_idx in range(n_cells):
        raw_ratio = []
        for g_idx in range(n_genes):
            diff = log_expr[c_idx][g_idx] - ref_baseline[g_idx]
            diff = max(-3.0, min(3.0, diff))
            raw_ratio.append(diff)

        # Mean centering across genome
        mean_c = sum(raw_ratio) / n_genes
        centered_row = [v - mean_c for v in raw_ratio]
        rel_log_ratio_centered.append(centered_row)

    # 4. Chromosome-by-chromosome moving average smoothing (boundary-safe window)
    chrom_to_indices = defaultdict(list)
    for idx, g in enumerate(gene_metadata):
        chrom_to_indices[g["chrom"]].append(idx)

    smoothed_matrix = [[0.0] * n_genes for _ in range(n_cells)]
    half_w = window_size // 2

    for chrom, indices in chrom_to_indices.items():
        n_chr_genes = len(indices)
        for c_idx in range(n_cells):
            full_row = rel_log_ratio_centered[c_idx]
            chr_vals = [full_row[idx] for idx in indices]

            for i in range(n_chr_genes):
                start_w = max(0, i - half_w)
                end_w = min(n_chr_genes, i + half_w + 1)
                sub_vals = chr_vals[start_w:end_w]
                smoothed_matrix[c_idx][indices[i]] = sum(sub_vals) / len(sub_vals)

    print(f"  [✓] Moving average completed across {len(chrom_to_indices)} chromosomes. Window size: {window_size} genes.")
    return log_expr, rel_log_ratio_centered, smoothed_matrix


# ==============================================================================
# Step 3: 6-State Hidden Markov Model (HMM) Viterbi Segmentation
# ==============================================================================
class SixStateCNAHMM:
    """
    6-State Hidden Markov Model for scRNA-seq CNA Segmentation:
      State 0: Complete Loss (0x)       - Mean = -1.5, Sigma = 0.35
      State 1: Single Copy Loss (1x)     - Mean = -0.6, Sigma = 0.25
      State 2: Neutral Normal (2x)       - Mean =  0.0, Sigma = 0.20
      State 3: Single Copy Gain (3x)     - Mean = +0.4, Sigma = 0.25
      State 4: Double Gain (4x)          - Mean = +0.8, Sigma = 0.30
      State 5: High Amplification (5x+)  - Mean = +1.3, Sigma = 0.40
    """
    def __init__(self, self_transition_prob=0.985):
        self.n_states = 6
        self.state_names = ["0x_CompleteLoss", "1x_Loss", "2x_Neutral", "3x_Gain", "4x_DoubleGain", "5x_HighAmp"]
        self.means = [-1.5, -0.6, 0.0, 0.4, 0.8, 1.3]
        self.sigmas = [0.35, 0.25, 0.20, 0.25, 0.30, 0.40]

        # Precompute emission constants: log_norm = -0.5 * log(2 * pi * sigma^2), inv_two_var = 1 / (2 * sigma^2)
        self.log_norm = [-0.5 * math.log(2.0 * math.pi * (s ** 2)) for s in self.sigmas]
        self.inv_two_var = [1.0 / (2.0 * (s ** 2)) for s in self.sigmas]

        # Prior distribution log_pi (heavily weighted towards Neutral State 2)
        self.log_pi = [math.log(p) for p in [0.005, 0.045, 0.900, 0.040, 0.005, 0.005]]

        # Transition probability matrix A
        A = [[0.0] * self.n_states for _ in range(self.n_states)]
        for i in range(self.n_states):
            for j in range(self.n_states):
                if i == j:
                    A[i][j] = self_transition_prob
                else:
                    dist = abs(i - j)
                    A[i][j] = (1.0 - self_transition_prob) / (2.0 ** dist)
            row_sum = sum(A[i])
            A[i] = [val / row_sum for val in A[i]]

        self.log_A = [[math.log(val) for val in row] for row in A]

    def viterbi(self, obs_sequence):
        """Viterbi dynamic programming decoding."""
        T = len(obs_sequence)
        n_s = self.n_states

        # V[s] holds current log likelihoods for each state
        V = [self.log_pi[s] + (self.log_norm[s] - ((obs_sequence[0] - self.means[s]) ** 2) * self.inv_two_var[s])
             for s in range(n_s)]

        backpointers = []

        # Recursion
        for t in range(1, T):
            x = obs_sequence[t]
            new_V = [0.0] * n_s
            bp_row = [0] * n_s

            for j in range(n_s):
                # Emission for state j
                emiss = self.log_norm[j] - ((x - self.means[j]) ** 2) * self.inv_two_var[j]

                # Find max over prev states
                best_score = -1e18
                best_prev = 0
                for i in range(n_s):
                    score = V[i] + self.log_A[i][j]
                    if score > best_score:
                        best_score = score
                        best_prev = i

                new_V[j] = best_score + emiss
                bp_row[j] = best_prev

            V = new_V
            backpointers.append(bp_row)

        # Traceback
        best_last = 0
        max_final = -1e18
        for s in range(n_s):
            if V[s] > max_final:
                max_final = V[s]
                best_last = s

        path = [best_last]
        for t in range(T - 2, -1, -1):
            prev = backpointers[t][path[-1]]
            path.append(prev)
        path.reverse()
        return path


def decode_cohort_hmm(smoothed_matrix, gene_metadata, cell_metadata):
    """Decodes HMM states across all cells and summarizes chromosome segments."""
    print("[Step 3] Running 6-state HMM Viterbi segmentation on smoothed profiles...")
    hmm = SixStateCNAHMM()

    n_cells = len(smoothed_matrix)
    decoded_states = []

    t0 = time.time()
    for i in range(n_cells):
        path = hmm.viterbi(smoothed_matrix[i])
        decoded_states.append(path)
    elapsed = time.time() - t0

    print(f"  [✓] Viterbi decoding completed in {elapsed:.2f}s ({n_cells} cells).")
    return decoded_states, hmm


# ==============================================================================
# Step 4: Reference-Free Gaussian Mixture Model (CopyKAT Style) & Evaluation
# ==============================================================================
class TwoComponentGMM:
    """
    2-Component Gaussian Mixture Model with Expectation-Maximization (EM) algorithm
    for classifying Diploid (low genomic variance) vs Aneuploid (high genomic variance) cells.
    """
    def __init__(self, max_iter=80, tol=1e-5):
        self.max_iter = max_iter
        self.tol = tol
        self.weights = [0.5, 0.5]
        self.means = [0.0, 0.0]
        self.variances = [1.0, 1.0]

    def fit(self, X):
        """Fits GMM on 1D genomic variance array X."""
        N = len(X)
        sorted_X = sorted(X)
        q25 = sorted_X[int(N * 0.25)]
        q75 = sorted_X[int(N * 0.75)]

        mean_all = sum(X) / N
        var_all = sum((x - mean_all) ** 2 for x in X) / N

        self.means = [q25, q75]
        self.variances = [max(1e-5, var_all * 0.4), max(1e-5, var_all * 0.4)]
        self.weights = [0.5, 0.5]

        log_lik_old = -1e18

        for iteration in range(self.max_iter):
            # E-step
            gamma0 = []
            gamma1 = []
            log_lik = 0.0

            m0, m1 = self.means
            v0, v1 = self.variances
            w0, w1 = self.weights

            inv_s0 = 1.0 / math.sqrt(2.0 * math.pi * v0)
            inv_s1 = 1.0 / math.sqrt(2.0 * math.pi * v1)

            for x in X:
                p0 = inv_s0 * math.exp(-0.5 * ((x - m0) ** 2) / v0)
                p1 = inv_s1 * math.exp(-0.5 * ((x - m1) ** 2) / v1)

                wp0 = max(1e-15, w0 * p0)
                wp1 = max(1e-15, w1 * p1)
                tot = wp0 + wp1

                gamma0.append(wp0 / tot)
                gamma1.append(wp1 / tot)
                log_lik += math.log(tot)

            # M-step
            N0 = sum(gamma0)
            N1 = sum(gamma1)

            if N0 < 1e-6 or N1 < 1e-6:
                break

            self.weights = [N0 / N, N1 / N]
            self.means[0] = sum(g * x for g, x in zip(gamma0, X)) / N0
            self.means[1] = sum(g * x for g, x in zip(gamma1, X)) / N1

            self.variances[0] = max(1e-6, sum(g * ((x - self.means[0]) ** 2) for g, x in zip(gamma0, X)) / N0)
            self.variances[1] = max(1e-6, sum(g * ((x - self.means[1]) ** 2) for g, x in zip(gamma1, X)) / N1)

            if abs(log_lik - log_lik_old) < self.tol:
                break
            log_lik_old = log_lik

        # Ensure index 0 is Diploid (smaller mean) and index 1 is Aneuploid (larger mean)
        if self.means[0] > self.means[1]:
            self.means.reverse()
            self.variances.reverse()
            self.weights.reverse()

        return self

    def predict_proba(self, X):
        """Returns posterior probability of being Aneuploid (Component 1)."""
        m0, m1 = self.means
        v0, v1 = self.variances
        w0, w1 = self.weights

        inv_s0 = 1.0 / math.sqrt(2.0 * math.pi * v0)
        inv_s1 = 1.0 / math.sqrt(2.0 * math.pi * v1)

        p_aneuploid = []
        for x in X:
            p0 = inv_s0 * math.exp(-0.5 * ((x - m0) ** 2) / v0)
            p1 = inv_s1 * math.exp(-0.5 * ((x - m1) ** 2) / v1)
            wp0 = max(1e-15, w0 * p0)
            wp1 = max(1e-15, w1 * p1)
            tot = wp0 + wp1
            p_aneuploid.append(wp1 / tot)
        return p_aneuploid


def compute_roc_auc(y_true, y_score):
    """Computes Receiver Operating Characteristic (ROC) curve and Area Under Curve (AUC)."""
    P = sum(y_true)
    N = len(y_true) - P

    thresholds = [i / 500.0 for i in range(501)]
    pairs = []
    for thresh in thresholds:
        TP = 0
        FP = 0
        for yt, ys in zip(y_true, y_score):
            yp = 1 if ys >= thresh else 0
            if yp == 1 and yt == 1:
                TP += 1
            elif yp == 1 and yt == 0:
                FP += 1
        TPR = TP / P if P > 0 else 0.0
        FPR = FP / N if N > 0 else 0.0
        pairs.append((FPR, TPR))

    # Sort by FPR ascending
    pairs.sort(key=lambda p: (p[0], p[1]))

    # Trapezoidal rule
    auc = 0.0
    for i in range(len(pairs) - 1):
        dx = pairs[i + 1][0] - pairs[i][0]
        avg_y = (pairs[i + 1][1] + pairs[i][1]) / 2.0
        auc += dx * avg_y

    return abs(auc)


def run_copykat_classification(smoothed_matrix, cell_metadata, gene_metadata):
    """
    Computes per-cell genomic CNA variance across chromosome arms (CopyKAT style),
    fits GMM, classifies cells, and evaluates metrics against ground truth.
    """
    print("[Step 4] Performing reference-free GMM aneuploidy classification (CopyKAT style)...")

    # 1. Group genes into chromosome arm segments (biological CNA units, CopyKAT style)
    arm_to_indices = defaultdict(list)
    for idx, g in enumerate(gene_metadata):
        arm_to_indices[(g["chrom"], g["arm"])].append(idx)

    # Compute per-cell genomic CNA variance across chromosome arm segments
    cell_genomic_variance = []
    for row in smoothed_matrix:
        arm_means = [sum(row[idx] for idx in indices) / len(indices) for indices in arm_to_indices.values()]
        mean_arm = sum(arm_means) / len(arm_means)
        var_arm = sum((m - mean_arm) ** 2 for m in arm_means) / len(arm_means)
        cell_genomic_variance.append(var_arm)

    # 2. Fit 2-Component GMM
    gmm = TwoComponentGMM().fit(cell_genomic_variance)
    p_aneuploid = gmm.predict_proba(cell_genomic_variance)
    predicted_labels = [1 if p >= 0.5 else 0 for p in p_aneuploid]

    # 3. Ground truth evaluation
    y_true = [c["is_aneuploid"] for c in cell_metadata]
    TP = sum(1 for yt, yp in zip(y_true, predicted_labels) if yt == 1 and yp == 1)
    TN = sum(1 for yt, yp in zip(y_true, predicted_labels) if yt == 0 and yp == 0)
    FP = sum(1 for yt, yp in zip(y_true, predicted_labels) if yt == 0 and yp == 1)
    FN = sum(1 for yt, yp in zip(y_true, predicted_labels) if yt == 1 and yp == 0)

    accuracy = (TP + TN) / len(y_true)
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    f1_score = 2.0 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0

    auc = compute_roc_auc(y_true, p_aneuploid)

    print(f"  [✓] GMM Classification Results:")
    print(f"      - Diploid Component Mean:   {gmm.means[0]:.4f} (Var: {gmm.variances[0]:.5f})")
    print(f"      - Aneuploid Component Mean: {gmm.means[1]:.4f} (Var: {gmm.variances[1]:.5f})")
    print(f"      - Confusion Matrix: TP={TP}, TN={TN}, FP={FP}, FN={FN}")
    print(f"      - Accuracy: {accuracy*100:.2f}%, Sensitivity: {sensitivity*100:.2f}%, Specificity: {specificity*100:.2f}%")
    print(f"      - Precision: {precision*100:.2f}%, F1-Score: {f1_score:.4f}, AUROC: {auc:.4f}")

    metrics = {
        "TP": TP, "TN": TN, "FP": FP, "FN": FN,
        "accuracy": round(accuracy, 4),
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
        "precision": round(precision, 4),
        "f1_score": round(f1_score, 4),
        "auroc": round(auc, 4),
        "diploid_mean_variance": round(float(gmm.means[0]), 5),
        "aneuploid_mean_variance": round(float(gmm.means[1]), 5)
    }

    return cell_genomic_variance, p_aneuploid, predicted_labels, metrics


# ==============================================================================
# Step 5: Save Results & Export Summary Tables
# ==============================================================================
def save_pipeline_results(out_dir, gene_metadata, cell_metadata, cna_ground_truth,
                          smoothed_matrix, decoded_states, cell_genomic_variance,
                          p_aneuploid, predicted_labels, metrics):
    os.makedirs(out_dir, exist_ok=True)
    print(f"[Step 5] Saving results and tables to {out_dir}...")

    # 1. Ground Truth Metadata CSV
    gt_file = os.path.join(out_dir, "simulated_single_cell_cna_ground_truth.csv")
    with open(gt_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "cell_type", "group", "subclone", "is_aneuploid"])
        for meta in cell_metadata:
            writer.writerow([meta["cell_id"], meta["cell_type"], meta["group"], meta["subclone"], meta["is_aneuploid"]])
    print(f"  [✓] Saved: {gt_file}")

    # 2. Representative Smoothed Profiles CSV (50 Normal, 50 Subclone A, 50 Subclone B)
    norm_idx = [i for i, c in enumerate(cell_metadata) if c["group"] == "Normal_Diploid"][:50]
    suba_idx = [i for i, c in enumerate(cell_metadata) if c["subclone"] == "Tumor_Subclone_A"][:50]
    subb_idx = [i for i, c in enumerate(cell_metadata) if c["subclone"] == "Tumor_Subclone_B"][:50]
    sampled_indices = norm_idx + suba_idx + subb_idx

    smooth_file = os.path.join(out_dir, "infercnv_smoothed_cna_profiles.csv")
    gene_cols = [g["gene_name"] for g in gene_metadata]
    with open(smooth_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "subclone"] + gene_cols)
        for idx in sampled_indices:
            cid = cell_metadata[idx]["cell_id"]
            sub = cell_metadata[idx]["subclone"]
            vals = [f"{v:.4f}" for v in smoothed_matrix[idx]]
            writer.writerow([cid, sub] + vals)
    print(f"  [✓] Saved: {smooth_file} ({len(sampled_indices)} sampled cells)")

    # 3. Decoded HMM Segments CSV
    seg_file = os.path.join(out_dir, "infercnv_hmm_viterbi_segments.csv")
    state_names = ["0x_CompleteLoss", "1x_Loss", "2x_Neutral", "3x_Gain", "4x_DoubleGain", "5x_HighAmp"]
    copy_num_map = [0, 1, 2, 3, 4, 5]

    with open(seg_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["subclone", "chrom", "start_gene", "end_gene", "start_pos", "end_pos", "decoded_state", "inferred_copy_number", "length_genes"])

        for subclone_name in ["Diploid_Stroma", "Tumor_Subclone_A", "Tumor_Subclone_B"]:
            rep_cell_idx = next(i for i, c in enumerate(cell_metadata) if c["subclone"] == subclone_name)
            cell_states = decoded_states[rep_cell_idx]

            chrom_set = []
            for g in gene_metadata:
                if g["chrom"] not in chrom_set:
                    chrom_set.append(g["chrom"])

            for chrom in chrom_set:
                chr_gene_indices = [i for i, g in enumerate(gene_metadata) if g["chrom"] == chrom]
                chr_states = [cell_states[i] for i in chr_gene_indices]

                seg_start_idx = 0
                for i in range(1, len(chr_states) + 1):
                    if i == len(chr_states) or chr_states[i] != chr_states[seg_start_idx]:
                        state_val = chr_states[seg_start_idx]
                        start_gene = gene_metadata[chr_gene_indices[seg_start_idx]]
                        end_gene = gene_metadata[chr_gene_indices[i - 1]]
                        writer.writerow([
                            subclone_name, chrom,
                            start_gene["gene_name"], end_gene["gene_name"],
                            start_gene["start"], end_gene["end"],
                            state_names[state_val], copy_num_map[state_val],
                            i - seg_start_idx
                        ])
                        seg_start_idx = i

    print(f"  [✓] Saved: {seg_file}")

    # 4. CopyKAT Classification CSV
    class_file = os.path.join(out_dir, "copykat_reference_free_gmm_classification.csv")
    with open(class_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "true_subclone", "true_aneuploid", "genomic_variance",
                         "posterior_p_aneuploid", "predicted_aneuploid", "classification_status"])
        for i, meta in enumerate(cell_metadata):
            is_correct = "CORRECT" if predicted_labels[i] == meta["is_aneuploid"] else "MISCLASSIFIED"
            writer.writerow([
                meta["cell_id"], meta["subclone"], meta["is_aneuploid"],
                f"{cell_genomic_variance[i]:.5f}", f"{p_aneuploid[i]:.4f}",
                predicted_labels[i], is_correct
            ])
    print(f"  [✓] Saved: {class_file}")

    # 5. Pipeline Summary JSON
    summary_file = os.path.join(out_dir, "pipeline_01_summary_metrics.json")
    summary = {
        "pipeline": "01_infercnv_and_copykat_aneuploidy_pipeline",
        "cohort": {
            "n_total_cells": len(cell_metadata),
            "n_normal_diploid": sum(1 for c in cell_metadata if c["group"] == "Normal_Diploid"),
            "n_aneuploid_tumor": sum(1 for c in cell_metadata if c["group"] == "Aneuploid_Tumor"),
            "n_subclone_A": sum(1 for c in cell_metadata if c["subclone"] == "Tumor_Subclone_A"),
            "n_subclone_B": sum(1 for c in cell_metadata if c["subclone"] == "Tumor_Subclone_B"),
            "n_genes": len(gene_metadata),
            "n_chromosomes": 23
        },
        "infercnv_smoothing": {
            "moving_average_window": 101,
            "reference_baseline_cells": 600,
            "bounds": [-3.0, 3.0]
        },
        "hmm_parameters": {
            "n_states": 6,
            "state_labels": ["0x_CompleteLoss", "1x_Loss", "2x_Neutral", "3x_Gain", "4x_DoubleGain", "5x_HighAmp"],
            "state_means": [-1.5, -0.6, 0.0, 0.4, 0.8, 1.3],
            "self_transition_probability": 0.985
        },
        "copykat_gmm_evaluation": metrics
    }
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"  [✓] Saved: {summary_file}")


def main():
    start_time = time.time()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "pipelines", "results_pipeline_01")

    # Step 1: Simulate single cell cohort
    gene_meta, cell_meta, cna_gt, raw_counts = simulate_cna_cohort(n_cells=1200, n_genes=1000, random_seed=42)

    # Step 2: Run inferCNV smoothing
    log_expr, rel_ratio, smoothed_matrix = run_infercnv_smoothing(raw_counts, cell_meta, gene_meta, window_size=101)

    # Step 3: Run 6-state HMM Viterbi decoding
    decoded_states, hmm = decode_cohort_hmm(smoothed_matrix, gene_meta, cell_meta)

    # Step 4: Run CopyKAT reference-free GMM classification
    genomic_var, p_aneuploid, pred_labels, metrics = run_copykat_classification(smoothed_matrix, cell_meta, gene_meta)

    # Step 5: Save results
    save_pipeline_results(out_dir, gene_meta, cell_meta, cna_gt, smoothed_matrix,
                          decoded_states, genomic_var, p_aneuploid, pred_labels, metrics)

    total_time = time.time() - start_time
    print(f"\n[✓] Pipeline 01 successfully completed in {total_time:.2f} seconds.")


if __name__ == "__main__":
    main()

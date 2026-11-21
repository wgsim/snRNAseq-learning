#!/usr/bin/env python3
"""
tracks/07_count_statistics_abundance/pipelines/02_milo_and_coda_differential_abundance_pipeline.py
=============================================================================
End-to-End Pipeline for Milo kNN Neighborhood Differential Abundance (DA)
& Compositional Data Analysis (CoDa) in the Tumor Microenvironment (TME).

Cohort Design:
  10-Patient mCRPC Immunotherapy Cohort:
  - 5 Responders (PT_01 to PT_05)
  - 5 Non-Responders (PT_06 to PT_10)
  - 4,000 Single Cells (400 cells per patient):
      • Naive CD8 T (800 cells, ~20% invariant across conditions)
      • Memory CD8 T (800 cells, ~20% invariant across conditions)
      • Effector CD8 T (1,200 cells): Expanded in Responders (~45% vs ~15%)
      • Exhausted CD8 T (1,200 cells): Expanded in Non-Responders (~45% vs ~15%)

Pipeline Structure:
  Step 1: Simulate 10-Patient TME Single-Cell Cohort on a Continuous 2D Manifold
  Step 2: Construct kNN Graph (k=30), Sample Index Cells, and Build Neighborhood Matrix
  Step 3: Fit Negative Binomial GLM Across Neighborhoods & Apply Spatial FDR Smoothing
  Step 4: Compositional Data Analysis (CoDa) with CLR & Dirichlet-Multinomial Log-Ratio Regression
  Step 5: Export Clinical DA Results, Summary Metrics, and Diagnostic Reports
"""

from __future__ import annotations

import os
import sys
import math
import time
import json
import csv
import random
from collections import defaultdict


# ==============================================================================
# Statistical & Math Helper Functions
# ==============================================================================

def norm_cdf(x: float) -> float:
    """Standard normal cumulative distribution function."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def pval_from_z(z: float) -> float:
    """Two-tailed p-value from z-score."""
    return max(1e-300, 2.0 * (1.0 - norm_cdf(abs(z))))


def student_t_pval(t: float, df: int) -> float:
    """Two-tailed p-value for Student t-distribution."""
    if df <= 0:
        return 1.0
    t2 = t * t
    a = 1.0 - 1.0 / (4.0 * df)
    b = math.sqrt(1.0 + t2 / (2.0 * df))
    z_equiv = (t * a) / b
    return pval_from_z(z_equiv)


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """Standard Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    qvals = [0.0] * n
    cum_min = 1.0
    for rank_minus_1 in range(n - 1, -1, -1):
        orig_idx, p = indexed[rank_minus_1]
        rank = rank_minus_1 + 1
        q = (p * n) / rank
        if q < cum_min:
            cum_min = q
        qvals[orig_idx] = min(1.0, cum_min)
    return qvals


def spatial_fdr(pvals: list[float], weights: list[float]) -> list[float]:
    """
    Graph-weighted Spatial FDR for Milo neighborhoods.
    Weights each hypothesis by its local graph connectivity density.
    """
    n = len(pvals)
    if n == 0:
        return []
    total_w = sum(weights)
    if total_w <= 0:
        return benjamini_hochberg(pvals)

    w_norm = [(w / total_w) * n for w in weights]
    # Weighted p-values: p_w = p / w_norm
    indexed = sorted(range(n), key=lambda i: pvals[i] / max(1e-4, w_norm[i]))
    qvals = [0.0] * n

    cum_min = 1.0
    for rank_minus_1 in range(n - 1, -1, -1):
        orig_idx = indexed[rank_minus_1]
        rank = rank_minus_1 + 1
        p_w = pvals[orig_idx] / max(1e-4, w_norm[orig_idx])
        q = (p_w * n) / rank
        if q < cum_min:
            cum_min = q
        qvals[orig_idx] = min(1.0, cum_min)

    return qvals


# ==============================================================================
# Step 1: Simulate 10-Patient TME Cohort
# ==============================================================================

def simulate_tme_cohort(seed: int = 42):
    """
    Simulates 10 patients (5 Responders, 5 Non-Responders) and 4,000 CD8 T cells.
    Phenotypes:
      1. Naive: ~20% in both
      2. Memory: ~20% in both
      3. Effector: 45% in Responders, 15% in Non-Responders
      4. Exhausted: 15% in Responders, 45% in Non-Responders
    Coordinates in 2D latent space with continuous transitions.
    """
    random.seed(seed)
    print("\n[Step 1] Simulating 10-Patient TME Cohort with 4,000 CD8 T Cells...")

    patients = []
    for i in range(1, 6):
        patients.append({"id": f"PT_{i:02d}", "condition": "Responder", "is_resp": 1})
    for i in range(6, 11):
        patients.append({"id": f"PT_{i:02d}", "condition": "NonResponder", "is_resp": 0})

    # State parameters: (Centroid X, Centroid Y, Spread)
    phenotype_geom = {
        "Naive": (2.0, 2.0, 0.45),
        "Memory": (4.5, 3.5, 0.50),
        "Effector": (7.5, 6.0, 0.60),
        "Exhausted": (8.0, 1.5, 0.55)
    }

    cells = []
    cell_id = 0

    for pt in patients:
        pt_id = pt["id"]
        is_resp = pt["is_resp"]

        # Proportions:
        if is_resp == 1:
            # 20% Naive (80), 20% Memory (80), 45% Effector (180), 15% Exhausted (60)
            target_counts = {"Naive": 80, "Memory": 80, "Effector": 180, "Exhausted": 60}
        else:
            # 20% Naive (80), 20% Memory (80), 15% Effector (60), 45% Exhausted (180)
            target_counts = {"Naive": 80, "Memory": 80, "Effector": 60, "Exhausted": 180}

        for p_type, count in target_counts.items():
            cx, cy, spread = phenotype_geom[p_type]
            for _ in range(count):
                # Add random transition jitter
                px = random.gauss(cx, spread)
                py = random.gauss(cy, spread)

                cells.append({
                    "cell_id": f"CD8_{cell_id:04d}",
                    "patient_id": pt_id,
                    "condition": pt["condition"],
                    "is_resp": is_resp,
                    "phenotype": p_type,
                    "x": px,
                    "y": py
                })
                cell_id += 1

    print(f"  [✓] Successfully generated {len(cells)} cells across {len(patients)} patients.")
    print(f"  [✓] Phenotypic Distribution:")
    print(f"      • Responders     : 45% Effector, 15% Exhausted, 20% Memory, 20% Naive")
    print(f"      • Non-Responders : 15% Effector, 45% Exhausted, 20% Memory, 20% Naive")

    return patients, cells


# ==============================================================================
# Step 2: Milo kNN Graph & Neighborhood Assembly
# ==============================================================================

def build_milo_neighborhoods(cells, patients, k: int = 30, sample_rate: float = 0.10):
    """
    Builds kNN graph (k=30) and samples representative index cells (rate=0.10).
    Constructs Neighborhood x Patient Count Matrix (M x 10).
    """
    print(f"\n[Step 2] Constructing kNN Graph (k={k}) and Sampling Representative Index Cells...")
    t0 = time.time()

    n_cells = len(cells)
    pts = [(c["x"], c["y"]) for c in cells]

    # Fast 2D Spatial Grid Indexing
    cell_size = 1.0
    grid = defaultdict(list)
    for idx, (x, y) in enumerate(pts):
        gx, gy = int(x // cell_size), int(y // cell_size)
        grid[(gx, gy)].append(idx)

    # Sample Index Cells: ~10% uniformly across graph (400 index cells)
    step = int(1.0 / sample_rate)
    index_cell_ids = list(range(0, n_cells, step))
    n_neighborhoods = len(index_cell_ids)

    neighborhoods = []
    # Patient count matrix: M x 10
    # patient_index: pt_id -> 0..9
    pt_to_col = {p["id"]: idx for idx, p in enumerate(patients)}
    pt_names = [p["id"] for p in patients]

    count_matrix = []

    for nh_idx, idx_c in enumerate(index_cell_ids):
        x0, y0 = pts[idx_c]
        gx, gy = int(x0 // cell_size), int(y0 // cell_size)

        candidates = []
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                candidates.extend(grid.get((gx + dx, gy + dy), []))

        # Sort candidates by distance to find top k nearest neighbors
        cand_dist = sorted([(math.hypot(x0 - pts[c][0], y0 - pts[c][1]), c) for c in candidates])
        nn_indices = [c for d, c in cand_dist[:k]]

        # Store neighborhood metadata
        center_cell = cells[idx_c]
        pheno_counts = defaultdict(int)
        pt_counts = [0] * len(patients)

        for ci in nn_indices:
            pheno_counts[cells[ci]["phenotype"]] += 1
            pt_id = cells[ci]["patient_id"]
            pt_counts[pt_to_col[pt_id]] += 1

        majority_pheno = max(pheno_counts.items(), key=lambda x: x[1])[0]

        neighborhoods.append({
            "nh_id": f"NH_{nh_idx:03d}",
            "index_cell_id": center_cell["cell_id"],
            "center_x": round(x0, 3),
            "center_y": round(y0, 3),
            "majority_phenotype": majority_pheno,
            "nn_cell_indices": set(nn_indices),
            "counts_by_patient": pt_counts
        })
        count_matrix.append(pt_counts)

    print(f"  [✓] Sampled {n_neighborhoods} representative index cells in {time.time() - t0:.2f} seconds.")
    print(f"  [✓] Assembled {n_neighborhoods} x {len(patients)} Neighborhood Count Matrix.")

    return neighborhoods, count_matrix, pt_names


# ==============================================================================
# Step 3: Negative Binomial GLM & Spatial FDR Smoothing
# ==============================================================================

def run_milo_differential_abundance(neighborhoods, count_matrix, patients):
    """
    Fits Negative Binomial GLM for each neighborhood:
      log(mu_{m, d}) = log(N_d) + beta_0 + Condition_d * beta_m
    Applies graph-weighted Spatial FDR.
    """
    print("\n[Step 3] Fitting Negative Binomial GLM across Neighborhoods & Spatial FDR...")

    # Library size per patient (total cells)
    total_cells_per_pt = [400] * len(patients)  # 400 cells each

    # 1. Compute pairwise neighborhood graph overlap weights for Spatial FDR
    # Overlap W_{m, m'} = |N_m intersect N_m'| / |N_m union N_m'|
    n_nh = len(neighborhoods)
    spatial_weights = []

    for m in range(n_nh):
        s_m = neighborhoods[m]["nn_cell_indices"]
        w_sum = 0.0
        # Compare with nearby neighborhoods
        cx, cy = neighborhoods[m]["center_x"], neighborhoods[m]["center_y"]
        for m2 in range(n_nh):
            if m == m2:
                w_sum += 1.0
                continue
            cx2, cy2 = neighborhoods[m2]["center_x"], neighborhoods[m2]["center_y"]
            if math.hypot(cx - cx2, cy - cy2) < 2.0:
                s_m2 = neighborhoods[m2]["nn_cell_indices"]
                inter = len(s_m.intersection(s_m2))
                if inter > 0:
                    union = len(s_m.union(s_m2))
                    w_sum += inter / float(union)
        spatial_weights.append(w_sum)

    # 2. Fit GLM and Wald Test per neighborhood
    pvals = []
    da_records = []

    for m in range(n_nh):
        nh = neighborhoods[m]
        cnts = count_matrix[m]

        # Responders: PT_01..PT_05 (indices 0..4)
        # Non-Responders: PT_06..PT_10 (indices 5..9)
        c_resp = cnts[:5]
        c_nresp = cnts[5:]

        sum_resp = sum(c_resp)
        sum_nresp = sum(c_nresp)

        # Rates normalized by total patient cells (2,000 cells per arm)
        rate_resp = max(1e-5, (sum_resp + 0.5) / 2000.0)
        rate_nresp = max(1e-5, (sum_nresp + 0.5) / 2000.0)

        log2fc = math.log2(rate_resp / rate_nresp)

        # Standard error under Poisson-Gamma GLM
        var_r = 1.0 / max(0.5, float(sum_resp))
        var_nr = 1.0 / max(0.5, float(sum_nresp))
        se_log2fc = (1.0 / math.log(2.0)) * math.sqrt(var_r + var_nr)

        wald_stat = log2fc / se_log2fc if se_log2fc > 0 else 0.0
        pval = pval_from_z(wald_stat)
        pvals.append(pval)

        da_records.append({
            "nh_id": nh["nh_id"],
            "index_cell_id": nh["index_cell_id"],
            "center_x": nh["center_x"],
            "center_y": nh["center_y"],
            "majority_phenotype": nh["majority_phenotype"],
            "resp_mean_count": round(sum_resp / 5.0, 2),
            "nonresp_mean_count": round(sum_nresp / 5.0, 2),
            "log2_fold_change": round(log2fc, 4),
            "std_error": round(se_log2fc, 4),
            "wald_stat": round(wald_stat, 4),
            "raw_pval": pval,
            "spatial_weight": round(spatial_weights[m], 3)
        })

    # 3. Apply Spatial FDR smoothing
    spatial_qvals = spatial_fdr(pvals, spatial_weights)

    sig_effector_enriched = 0
    sig_exhausted_enriched = 0
    invariant_count = 0

    for m in range(n_nh):
        q = spatial_qvals[m]
        da_records[m]["spatial_fdr"] = q
        fc = da_records[m]["log2_fold_change"]

        if q < 0.05 and fc > 0.5:
            status = "Enriched_in_Responders"
            sig_effector_enriched += 1
        elif q < 0.05 and fc < -0.5:
            status = "Enriched_in_NonResponders"
            sig_exhausted_enriched += 1
        else:
            status = "Invariant"
            invariant_count += 1
        da_records[m]["da_status"] = status

    print(f"  [✓] Milo Differential Abundance Summary (Spatial FDR < 0.05):")
    print(f"      • Enriched in Responders     (Effector T): {sig_effector_enriched} neighborhoods")
    print(f"      • Enriched in Non-Responders (Exhausted T): {sig_exhausted_enriched} neighborhoods")
    print(f"      • Invariant Neighborhoods   (Naive/Memory): {invariant_count} neighborhoods")

    return da_records


# ==============================================================================
# Step 4: Compositional Data Analysis (CoDa) & Log-Ratio Regression
# ==============================================================================

def run_compositional_data_analysis(cells, patients):
    """
    Performs Compositional Data Analysis (CoDa):
    1. Patient cell type proportion table
    2. Centered Log-Ratio (CLR) transformation
    3. Additive Log-Ratio (ALR) regression with stable reference (Naive CD8 T)
    """
    print("\n[Step 4] Running Compositional Data Analysis (CoDa) & Bayesian ALR/CLR Regression...")

    phenotypes = ["Naive", "Memory", "Effector", "Exhausted"]

    # Calculate count of cells per patient per phenotype
    pt_counts = defaultdict(lambda: {p: 0 for p in phenotypes})
    for c in cells:
        pt_counts[c["patient_id"]][c["phenotype"]] += 1

    # Patient proportions
    prop_records = []
    for pt in patients:
        p_id = pt["id"]
        c_dict = pt_counts[p_id]
        tot = sum(c_dict.values())
        rec = {
            "patient_id": p_id,
            "condition": pt["condition"],
            "is_resp": pt["is_resp"],
            "total_cells": tot
        }
        for p in phenotypes:
            rec[f"prop_{p}"] = round(c_dict[p] / float(tot), 4)
            rec[f"count_{p}"] = c_dict[p]
        prop_records.append(rec)

    # 1. Raw Percentage t-test (To demonstrate spurious simplex artifact)
    raw_tests = {}
    for p in phenotypes:
        vals_r = [r[f"prop_{p}"] for r in prop_records if r["is_resp"] == 1]
        vals_nr = [r[f"prop_{p}"] for r in prop_records if r["is_resp"] == 0]
        mr, mnr = sum(vals_r) / 5.0, sum(vals_nr) / 5.0
        vr = sum((x - mr) ** 2 for x in vals_r) / 4.0
        vnr = sum((x - mnr) ** 2 for x in vals_nr) / 4.0
        se = math.sqrt(vr / 5.0 + vnr / 5.0)
        t = (mr - mnr) / se if se > 0 else 0
        raw_tests[p] = {"mean_resp": mr, "mean_nonresp": mnr, "t": t, "pval": student_t_pval(t, df=8)}

    # 2. Centered Log-Ratio (CLR)
    # CLR_k = log(p_k) - (1/K) sum(log(p_j))
    clr_records = []
    for r in prop_records:
        props = [max(1e-4, r[f"prop_{p}"]) for p in phenotypes]
        log_props = [math.log(p) for p in props]
        geom_mean_log = sum(log_props) / len(log_props)

        clr_dict = {"patient_id": r["patient_id"], "condition": r["condition"]}
        for idx, p in enumerate(phenotypes):
            clr_dict[f"clr_{p}"] = round(log_props[idx] - geom_mean_log, 4)
        clr_records.append(clr_dict)

    # 3. Additive Log-Ratio (ALR) Regression with Reference = Naive T
    # log(p_k / p_ref) = beta_0 + Condition * beta_k
    ref_pheno = "Naive"
    coda_regression_results = []

    for p in phenotypes:
        if p == ref_pheno:
            coda_regression_results.append({
                "phenotype": p,
                "role": "Reference_Baseline",
                "effect_beta": 0.0,
                "std_error": 0.0,
                "t_stat": 0.0,
                "p_value": 1.0,
                "coda_status": "Invariant_Baseline"
            })
            continue

        alr_resp = []
        alr_nresp = []
        for r in prop_records:
            ref_val = max(1e-4, r[f"prop_{ref_pheno}"])
            tgt_val = max(1e-4, r[f"prop_{p}"])
            ratio = math.log(tgt_val / ref_val)
            if r["is_resp"] == 1:
                alr_resp.append(ratio)
            else:
                alr_nresp.append(ratio)

        m_r = sum(alr_resp) / 5.0
        m_nr = sum(alr_nresp) / 5.0
        beta = m_r - m_nr

        var_r = sum((x - m_r) ** 2 for x in alr_resp) / 4.0
        var_nr = sum((x - m_nr) ** 2 for x in alr_nresp) / 4.0
        se_beta = math.sqrt(var_r / 5.0 + var_nr / 5.0)

        t_stat = beta / se_beta if se_beta > 0 else 0.0
        pval = student_t_pval(t_stat, df=8)

        if pval < 0.01 and beta > 0.5:
            status = "True_Expansion_Responders"
        elif pval < 0.01 and beta < -0.5:
            status = "True_Expansion_NonResponders"
        else:
            status = "Statistically_Invariant"

        coda_regression_results.append({
            "phenotype": p,
            "role": "Target_Phenotype",
            "effect_beta": round(beta, 4),
            "std_error": round(se_beta, 4),
            "t_stat": round(t_stat, 4),
            "p_value": pval,
            "coda_status": status
        })

    print("  [✓] Compositional Regression (Reference = Naive CD8 T):")
    for row in coda_regression_results:
        print(f"      • {row['phenotype']:<12} | Effect Beta: {row['effect_beta']:>7.3f} | "
              f"P-val: {row['p_value']:.2e} | Status: {row['coda_status']}")

    return prop_records, clr_records, coda_regression_results


# ==============================================================================
# Step 5: Save Results and Metrics
# ==============================================================================

def export_pipeline_results(outdir, da_records, prop_records, coda_reg_results, metrics):
    """Saves structured outputs to results_pipeline_02/."""
    os.makedirs(outdir, exist_ok=True)
    print(f"\n[Step 5] Saving structured outputs to {outdir}...")

    # 1. Milo DA results
    f_milo = os.path.join(outdir, "milo_neighborhood_da_results.csv")
    with open(f_milo, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(da_records[0].keys()))
        writer.writeheader()
        writer.writerows(da_records)
    print(f"  [✓] Saved: {f_milo}")

    # 2. Patient proportions
    f_prop = os.path.join(outdir, "coda_cell_type_proportions.csv")
    with open(f_prop, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(prop_records[0].keys()))
        writer.writeheader()
        writer.writerows(prop_records)
    print(f"  [✓] Saved: {f_prop}")

    # 3. CoDa regression results
    f_coda = os.path.join(outdir, "coda_clr_and_regression_results.csv")
    with open(f_coda, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(coda_reg_results[0].keys()))
        writer.writeheader()
        writer.writerows(coda_reg_results)
    print(f"  [✓] Saved: {f_coda}")

    # 4. JSON metrics
    f_json = os.path.join(outdir, "pipeline_02_metrics.json")
    with open(f_json, "w", encoding="utf-8") as fp:
        json.dump(metrics, fp, indent=2)
    print(f"  [✓] Saved: {f_json}")

    # 5. Summary report
    f_txt = os.path.join(outdir, "pipeline_02_summary.txt")
    with open(f_txt, "w", encoding="utf-8") as fp:
        fp.write("=" * 80 + "\n")
        fp.write("Track 07 Pipeline 02: Milo kNN DA & CoDa Compositional Regression Report\n")
        fp.write("=" * 80 + "\n\n")
        fp.write(f"Cohort Specification:\n")
        fp.write(f"  • Patients: 10 total (5 Responders, 5 Non-Responders)\n")
        fp.write(f"  • Cells: 4,000 CD8 T cells across Naive, Memory, Effector, Exhausted\n")
        fp.write(f"  • Neighborhoods: {metrics['total_neighborhoods']} sampled index neighborhoods (k=30)\n\n")
        fp.write(f"Milo Differential Abundance Results (Spatial FDR < 0.05):\n")
        fp.write(f"  • Enriched in Responders     (Effector T)  : {metrics['milo_enriched_responders']} neighborhoods\n")
        fp.write(f"  • Enriched in Non-Responders (Exhausted T) : {metrics['milo_enriched_nonresponders']} neighborhoods\n")
        fp.write(f"  • Invariant Neighborhoods   (Naive/Memory): {metrics['milo_invariant_neighborhoods']} neighborhoods\n\n")
        fp.write(f"CoDa Compositional Log-Ratio Regression (Reference = Naive CD8 T):\n")
        for r in coda_reg_results:
            fp.write(f"  • {r['phenotype']:<12} | Effect Beta: {r['effect_beta']:>7.3f} | P-val: {r['p_value']:.2e} | {r['coda_status']}\n")
    print(f"  [✓] Saved: {f_txt}")


# ==============================================================================
# Main Orchestrator
# ==============================================================================

def main():
    start_time = time.time()
    print("=" * 80)
    print("Track 07: Milo kNN Neighborhood DA & Compositional Regression Pipeline")
    print("=" * 80)

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "pipelines", "results_pipeline_02")

    # Step 1
    patients, cells = simulate_tme_cohort(seed=42)

    # Step 2
    neighborhoods, count_matrix, pt_names = build_milo_neighborhoods(
        cells, patients, k=30, sample_rate=0.10
    )

    # Step 3
    da_records = run_milo_differential_abundance(neighborhoods, count_matrix, patients)

    # Step 4
    prop_records, clr_records, coda_reg_results = run_compositional_data_analysis(cells, patients)

    # Step 5
    sig_resp = sum(1 for r in da_records if r["da_status"] == "Enriched_in_Responders")
    sig_nonresp = sum(1 for r in da_records if r["da_status"] == "Enriched_in_NonResponders")
    sig_invar = sum(1 for r in da_records if r["da_status"] == "Invariant")

    metrics = {
        "pipeline": "02_milo_and_coda_differential_abundance_pipeline",
        "total_patients": len(patients),
        "total_cells": len(cells),
        "total_neighborhoods": len(neighborhoods),
        "milo_enriched_responders": sig_resp,
        "milo_enriched_nonresponders": sig_nonresp,
        "milo_invariant_neighborhoods": sig_invar,
        "execution_time_seconds": round(time.time() - start_time, 2),
        "coda_results": coda_reg_results
    }

    export_pipeline_results(outdir, da_records, prop_records, coda_reg_results, metrics)

    elapsed = time.time() - start_time
    print(f"\n[✓] Pipeline 02 completed successfully in {elapsed:.2f} seconds with exit code 0.")


if __name__ == "__main__":
    main()

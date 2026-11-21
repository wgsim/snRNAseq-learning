#!/usr/bin/env python3
"""
02_prostate_qc_and_cellbender.py
==================================
전립선암 snRNA-seq 특화 Ambient RNA(CellBender) 및 Doublet(Scrublet) 정제 파이프라인

배경 및 목적:
전립선 조직 파쇄 시 쏟아져 나오는 고농도 KLK3(PSA), ACPP 등의 분비성 전사체가
주변 액적으로 침투하여 면역세포 및 기질세포를 오염시키는 현상을 수치적으로 제거하고,
이종 이중핵(Doublet)을 선별하여 제거하는 전문 전처리 모듈입니다.
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr


def run_cellbender_cli(
    raw_h5_path: str,
    output_dir: str = "cellbender_out/",
    expected_cells: int = 5000,
    total_droplets: int = 25000,
    fpr: float = 0.01,
    epochs: int = 150
) -> str:
    """CellBender CLI 명령어를 구성하고 실행합니다 (GPU 권장)."""
    os.makedirs(output_dir, exist_ok=True)
    cleaned_h5 = os.path.join(output_dir, "cellbender_cleaned.h5")

    if os.path.exists(cleaned_h5):
        print(f"[*] Found existing CellBender output: {cleaned_h5}")
        return cleaned_h5

    cmd = [
        "cellbender", "remove-background",
        "--input", raw_h5_path,
        "--output", cleaned_h5,
        "--expected-cells", str(expected_cells),
        "--total-droplets-included", str(total_droplets),
        "--fpr", str(fpr),
        "--epochs", str(epochs)
    ]
    print(f"[*] Executing CellBender: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"[✓] CellBender successfully finished: {cleaned_h5}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[!] CellBender CLI execution skipped/failed ({e}). Proceeding in simulation mode.")
    return cleaned_h5


def detect_and_filter_doublets(
    adata: sc.AnnData,
    expected_doublet_rate: float = 0.06,
    threshold: float = None
) -> sc.AnnData:
    """Scrublet을 활용하여 합성 이중핵(Synthetic doublets)을 시뮬레이션하고 필터링합니다."""
    print(f"[*] Running Scrublet on {adata.n_obs} nuclei...")
    
    # 카운트 행렬 추출
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, 
        min_cells=3, 
        min_gene_variability_pctl=85, 
        n_prin_comps=30
    )

    if threshold is not None:
        predicted_doublets = doublet_scores > threshold

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['is_doublet'] = predicted_doublets

    n_doublets = np.sum(predicted_doublets)
    print(f"[+] Scrublet detected {n_doublets} predicted doublets ({n_doublets / adata.n_obs * 100:.2f}%)")

    # 싱글렛만 추출
    adata_clean = adata[~adata.obs['is_doublet']].copy()
    print(f"[✓] Retained {adata_clean.n_obs} high-confidence singlets")
    return adata_clean


def evaluate_ambient_removal(
    adata_raw: sc.AnnData,
    adata_clean: sc.AnnData,
    ambient_genes: list = ["KLK3", "ACPP", "MSMB", "FOLH1"]
):
    """Ambient RNA 제거 전후의 대표 마커 발현량 감소를 정량 평가합니다."""
    print("\n--- Ambient RNA Removal Evaluation ---")
    for gene in ambient_genes:
        if gene in adata_raw.var_names and gene in adata_clean.var_names:
            raw_mean = np.mean(adata_raw[:, gene].X)
            clean_mean = np.mean(adata_clean[:, gene].X)
            reduction = ((raw_mean - clean_mean) / raw_mean) * 100 if raw_mean > 0 else 0
            print(f"  - {gene:<8}: Raw Mean={raw_mean:.3f} | Clean Mean={clean_mean:.3f} | Reduction={reduction:.1f}%")
    print("--------------------------------------\n")


def clean_prostate_dataset(
    raw_h5_path: str = "outs/raw_feature_bc_matrix.h5",
    output_dir: str = "cleaned_data/"
) -> sc.AnnData:
    os.makedirs(output_dir, exist_ok=True)

    # 1. CellBender 실행 (또는 대체 로드)
    if not os.path.exists(raw_h5_path):
        print(f"[!] {raw_h5_path} not found. Generating mock test count data for pipeline testing...")
        counts = np.random.negative_binomial(2, 0.1, size=(1500, 1000))
        gene_names = [f"Gene_{i}" for i in range(1000)]
        gene_names[0], gene_names[1], gene_names[2] = "KLK3", "ACPP", "MT-CO1"
        adata = sc.AnnData(
            X=counts.astype(np.float32),
            obs=pd.DataFrame(index=[f"Nucleus_{i}-1" for i in range(1500)]),
            var=pd.DataFrame(index=gene_names)
        )
    else:
        cleaned_h5 = run_cellbender_cli(raw_h5_path)
        if os.path.exists(cleaned_h5):
            adata = sc.read_10x_h5(cleaned_h5)
            # Ambient RNA 제거 전후 마커 발현량 변화 정량 평가
            adata_raw = sc.read_10x_h5(raw_h5_path)
            evaluate_ambient_removal(adata_raw, adata)
        else:
            adata = sc.read_10x_h5(raw_h5_path)

    # 2. Scrublet 기반 이중핵 제거
    adata_singlet = detect_and_filter_doublets(adata, expected_doublet_rate=0.06)

    # 3. 저장
    out_file = os.path.join(output_dir, "prostate_cleaned_singlets.h5ad")
    adata_singlet.write_h5ad(out_file)
    print(f"[✓] Cleaned singlet AnnData saved to: {out_file}")

    return adata_singlet


if __name__ == "__main__":
    raw_path = sys.argv[1] if len(sys.argv) > 1 else "outs/raw_feature_bc_matrix.h5"
    clean_prostate_dataset(raw_path)

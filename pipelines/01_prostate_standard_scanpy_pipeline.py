#!/usr/bin/env python3
"""
01_prostate_standard_scanpy_pipeline.py
=========================================
전립선암(원발암 & mCRPC) snRNA-seq 표준 다운스트림 분석 파이프라인

핵심 기능:
1. AnnData 객체 로딩 및 핵 특화 QC 지표(미토콘드리아 1~2% 상한, 유전자/UMI 분포) 계산
2. Shifted Log-정규화 및 고변이 유전자(HVG: 2,500개) 선별
3. 주성분 분석(PCA, 35 PCs) 및 kNN 그래프 구축
4. Leiden 군집화(연결 보장 그래프 클러스터링) 및 2D UMAP 시각화
5. 전립선암 표준 마커(Epithelial, CAF, Immune, Malignant) 발현 확인 및 Wilcoxon DEG 산출
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# 시각화 설정
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)


def run_standard_pipeline(
    matrix_dir: str = "outs/filtered_feature_bc_matrix/",
    output_dir: str = "results_scanpy/",
    sample_id: str = "Prostate_Tumor_01"
) -> sc.AnnData:
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir
    print(f"[*] Starting snRNA-seq pipeline for sample: {sample_id}")

    # ---------------------------------------------------------
    # Step 1: 데이터 로드
    # ---------------------------------------------------------
    if not os.path.exists(matrix_dir):
        print(f"[!] Warning: {matrix_dir} not found. Creating synthetic demo AnnData for validation.")
        # 파이프라인 검증용 합성 데이터셋 생성
        np.random.seed(42)
        n_cells, n_genes = 1200, 3000
        counts = np.random.negative_binomial(2, 0.1, size=(n_cells, n_genes))
        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        # 전립선 핵심 마커 유전자 주입
        prostate_markers = [
            "KLK3", "AR", "KRT8", "KRT18", "TP63", "KRT5", 
            "ACTA2", "COL1A1", "PECAM1", "CD3D", "CD68", "AMACR", "MT-CO1", "MT-ND1"
        ]
        for idx, m in enumerate(prostate_markers):
            if idx < n_genes:
                gene_names[idx] = m
                
        adata = sc.AnnData(
            X=counts.astype(np.float32),
            obs=pd.DataFrame(index=[f"Nucleus_{i}-1" for i in range(n_cells)]),
            var=pd.DataFrame(index=gene_names)
        )
        adata.var_names_make_unique()
    else:
        adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()

    print(f"[+] Initial dataset: {adata.n_obs} nuclei x {adata.n_vars} genes")

    # ---------------------------------------------------------
    # Step 2: snRNA-seq 특화 QC 계산 및 필터링
    # ---------------------------------------------------------
    # 미토콘드리아 및 리보솜 유전자 플래그
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'rps', 'rpl'))

    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt', 'ribo'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )

    print(f"[*] Pre-filter Median Mito%: {np.median(adata.obs['pct_counts_mt']):.2f}%")
    print(f"[*] Pre-filter Median Genes: {np.median(adata.obs['n_genes_by_counts']):.0f}")

    # snRNA-seq 특화 필터링 룰:
    # 1. 미토콘드리아는 2.5% 이하 (이상적 목표 2.0% 미만; 핵막에 부착된 세포질 오염 배제)
    # 2. UMI는 면역세포 보존을 위해 최소 700 이상
    n_before = adata.n_obs
    adata = adata[
        (adata.obs['pct_counts_mt'] < 2.5) &
        (adata.obs['n_genes_by_counts'] >= 500) &
        (adata.obs['n_genes_by_counts'] <= 6000) &
        (adata.obs['total_counts'] >= 700)
    ].copy()
    print(f"[+] Quality Filter: {n_before} -> {adata.n_obs} nuclei retained (Filtered out {n_before - adata.n_obs})")

    # 극소수 세포에서만 발현되는 노이즈 유전자 제거
    sc.pp.filter_genes(adata, min_cells=3)

    # ---------------------------------------------------------
    # Step 3: 정규화, 로그 변환 및 고변이 유전자(HVG) 선별
    # ---------------------------------------------------------
    # 원시 카운트 백업 (layers에 저장)
    adata.layers["raw_counts"] = adata.X.copy()

    # 라이브러리 크기 정규화 (10,000 기준) 및 log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 고변이 유전자 (HVG: 2,500개)
    try:
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=2500, 
            flavor='seurat_v3',
            layer="raw_counts"
        )
    except (ImportError, Exception):
        print("[!] 'scikit-misc' not available for flavor='seurat_v3'. Falling back to flavor='seurat'.")
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=2500, 
            flavor='seurat'
        )
    print(f"[+] Identified {np.sum(adata.var['highly_variable'])} highly variable genes")

    # ---------------------------------------------------------
    # Step 4: PCA 차원 축소 및 이웃 그래프 (kNN) 구축
    # ---------------------------------------------------------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=35, use_highly_variable=True, svd_solver='arpack')

    # kNN 그래프 구축: 전립선 조직의 복합 계통 보존을 위해 n_neighbors=25 사용
    sc.pp.neighbors(adata, n_neighbors=25, n_pcs=35)

    # ---------------------------------------------------------
    # Step 5: Leiden 클러스터링 및 UMAP 투영
    # ---------------------------------------------------------
    # Louvain의 분리 노드 오류를 방지하는 Leiden 알고리즘 사용
    sc.tl.leiden(adata, resolution=0.6, key_added="leiden_res0.6")
    sc.tl.umap(adata, min_dist=0.3)

    print(f"[+] Leiden clustering resolved {len(adata.obs['leiden_res0.6'].unique())} clusters")

    # ---------------------------------------------------------
    # Step 6: 전립선 표준 마커 발현 확인 및 DEG 계산
    # ---------------------------------------------------------
    canonical_markers = [
        "KLK3", "AR", "AMACR",     # Malignant / Luminal
        "TP63", "KRT5",            # Basal
        "ACTA2", "COL1A1",         # CAF / Stroma
        "PECAM1",                  # Endothelial
        "CD3D", "CD68"             # T cell, Macrophage
    ]
    # 실제 데이터셋에 존재하는 마커만 필터
    available_markers = [m for m in canonical_markers if m in adata.var_names]
    if available_markers:
        sc.pl.dotplot(
            adata, 
            var_names=available_markers, 
            groupby='leiden_res0.6', 
            show=False, 
            save=f"_{sample_id}_markers_dotplot.png"
        )

    # Wilcoxon Rank-Sum 검정을 통한 각 클러스터별 마커 유전자 산출
    sc.tl.rank_genes_groups(adata, groupby='leiden_res0.6', method='wilcoxon', pts=True)
    
    # 상위 DEG 테이블 추출
    deg_df = sc.get.rank_genes_groups_df(adata, group=None)
    deg_file = os.path.join(output_dir, f"{sample_id}_deg_markers.csv")
    deg_df.to_csv(deg_file, index=False)
    print(f"[+] Differential expression table saved to: {deg_file}")

    # 최종 AnnData 저장
    h5ad_path = os.path.join(output_dir, f"{sample_id}_processed.h5ad")
    adata.write_h5ad(h5ad_path)
    print(f"[✓] Pipeline complete! Output written to: {h5ad_path}")

    return adata


if __name__ == "__main__":
    matrix_path = sys.argv[1] if len(sys.argv) > 1 else "outs/filtered_feature_bc_matrix/"
    run_standard_pipeline(matrix_dir=matrix_path)

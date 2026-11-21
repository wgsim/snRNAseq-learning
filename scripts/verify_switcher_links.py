#!/usr/bin/env python3
import os
import re

WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # ENG-025: 스크립트 위치 기준(이식 가능)

TRACK_FILES = [
    ("00", "snRNAseq_Prostate_Master_Course.html"),
    ("01", "tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.html"),
    ("02", "tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.html"),
    ("03", "tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.html"),
    ("04", "tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.html"),
    ("05", "tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html"),
    ("06", "tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html"),
    ("07", "tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.html"),
    ("08", "tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.html"),
    ("09", "tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.html"),
    ("10", "tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.html")
]

print("=== Verifying Global Track Switcher Links ===")
for tid, rel_path in TRACK_FILES:
    abs_path = os.path.join(WORKSPACE_ROOT, rel_path)
    assert os.path.exists(abs_path), f"File missing: {abs_path}"
    with open(abs_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    assert "global-track-switcher" in content, f"Missing sidebar switcher in {rel_path}"
    assert "global-curriculum-banner" in content, f"Missing top banner in {rel_path}"
    
    # Extract all relative links in switcher & banner
    curr_dir = os.path.dirname(abs_path)
    # Match value="..." or href="..." ending in .html
    hrefs = re.findall(r'(?:href|value)=["\']([^"\']+\.html)["\']', content)
    
    print(f"Track {tid} ({rel_path}): found {len(hrefs)} .html links in switcher/banner")
    assert len(hrefs) >= 20, f"Expected at least 20 links, found {len(hrefs)}"
    
    for h in hrefs:
        target = os.path.normpath(os.path.join(curr_dir, h))
        assert os.path.exists(target), f"Broken link in {rel_path}: {h} -> {target}"
    
    # Check 0 external CDNs
    cdns = re.findall(r'<(?:script|link)[^>]*(?:src|href)=["\'](https?://[^"\']+)["\']', content)
    assert len(cdns) == 0, f"Found external CDNs in {rel_path}: {cdns}"
    print(f"[OK] Track {tid}: {rel_path} (All {len(hrefs)} links verified, 0 CDNs, Switcher & Banner intact)")

print("\n=== ALL 11 MASTER COURSES VERIFIED SUCCESSFULLY ===")

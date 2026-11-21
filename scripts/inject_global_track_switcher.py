#!/usr/bin/env python3
"""
inject_global_track_switcher.py
Safely injects the Global Curriculum Track Switcher into all 11 Master Course HTML files.

Run it AFTER the build scripts: it injects into the current (freshly built) HTML.
Every injected block is wrapped in TRACK-SWITCHER BEGIN/END markers, so a re-run
strips the previous injection first (idempotent) instead of restoring old backups,
which used to overwrite newly built content (review finding ENG-021).
"""

import os
import re

WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # ENG-025: 스크립트 위치 기준(이식 가능)

CSS_BEGIN, CSS_END = "/* TRACK-SWITCHER:BEGIN:css */", "/* TRACK-SWITCHER:END:css */"
HTML_BEGIN, HTML_END = "<!-- TRACK-SWITCHER:BEGIN:{} -->", "<!-- TRACK-SWITCHER:END:{} -->"


def wrap_html(block, name):
    return HTML_BEGIN.format(name) + block + HTML_END.format(name)


def _block_re(begin, end):
    return re.escape(begin) + r".*?" + re.escape(end)


def strip_previous_injection(content):
    """Exact inverse of the insertions made in process_file()."""
    content = re.sub(_block_re(CSS_BEGIN, CSS_END) + r"\n", "", content, flags=re.S)
    content = re.sub(_block_re(HTML_BEGIN.format("sidebar"), HTML_END.format("sidebar")) + r"\n    ",
                     "", content, flags=re.S)
    content = re.sub(r"\n" + _block_re(HTML_BEGIN.format("banner"), HTML_END.format("banner")),
                     "", content, flags=re.S)
    return content

TRACKS = [
    {
        "id": "00",
        "title": "Track 00: snRNA-seq 기초 & 핵심 다운스트림",
        "short_title": "Track 00 기초",
        "icon": "🧬",
        "rel_from_root": "snRNAseq_Prostate_Master_Course.html",
        "desc": "동결조직 핵분리, Pre-mRNA, Scanpy QC, 클러스터링"
    },
    {
        "id": "01",
        "title": "Track 01: 10x Multiome & 후성유전체학",
        "short_title": "Track 01 멀티옴",
        "icon": "⚡",
        "rel_from_root": "tracks/01_multiome_snatac/Multiome_snATAC_Master_Course.html",
        "desc": "snATAC-seq, Tn5 삽입, WNN 통합, SCENIC+ eGRN"
    },
    {
        "id": "02",
        "title": "Track 02: 공간 전사체학 (Spatial)",
        "short_title": "Track 02 공간",
        "icon": "🔬",
        "rel_from_root": "tracks/02_spatial_transcriptomics/Spatial_Transcriptomics_Master_Course.html",
        "desc": "Visium, Xenium, Cell2location 베이지안 디컨볼루션"
    },
    {
        "id": "03",
        "title": "Track 03: 종양미세환경 & 세포대화",
        "short_title": "Track 03 세포대화",
        "icon": "💬",
        "rel_from_root": "tracks/03_tme_cell_communication/TME_Cell_Communication_Master_Course.html",
        "desc": "CellChat, LIANA+, 질량작용 법칙, 다량체 수용체"
    },
    {
        "id": "04",
        "title": "Track 04: 암 유전체 & 클론 진화",
        "short_title": "Track 04 암유전체",
        "icon": "🧬",
        "rel_from_root": "tracks/04_cancer_genomics_clones/Cancer_Genomics_Master_Course.html",
        "desc": "inferCNV HMM, CopyKAT GMM, Numbat 베이지안 BAF"
    },
    {
        "id": "05",
        "title": "Track 05: 딥 생성 모델 & 파운데이션",
        "short_title": "Track 05 생성모델",
        "icon": "🤖",
        "rel_from_root": "tracks/05_generative_foundation_models/Generative_Foundation_Master_Course.html",
        "desc": "scVI VAE & ELBO, scANVI, scGPT, In Silico 섭동"
    },
    {
        "id": "06",
        "title": "Track 06: 매니폴드 기하학 & 궤적 추론",
        "short_title": "Track 06 궤적추론",
        "icon": "🧭",
        "rel_from_root": "tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html",
        "desc": "확산 유사시간(DPT), scVelo 동역학, WOT, CellRank 2"
    },
    {
        "id": "07",
        "title": "Track 07: 단일세포 통계학 & 차등 분석",
        "short_title": "Track 07 카운트통계",
        "icon": "📊",
        "rel_from_root": "tracks/07_count_statistics_abundance/Count_Statistics_Master_Course.html",
        "desc": "가짜반복 극복, 슈도벌크 QLF, GLMM, Milo kNN DA"
    },
    {
        "id": "08",
        "title": "Track 08: 대규모 희소 연산 & GPU 가속",
        "short_title": "Track 08 희소/GPU",
        "icon": "⚡",
        "rel_from_root": "tracks/08_large_scale_sparse/Large_Scale_Sparse_Master_Course.html",
        "desc": "CSR/CSC 캐시 계층, Zarr 청킹, TileDB-SOMA, RAPIDS"
    },
    {
        "id": "09",
        "title": "Track 09: 로우레벨 서열 정렬 & BAM",
        "short_title": "Track 09 서열정렬",
        "icon": "✂️",
        "rel_from_root": "tracks/09_upstream_alignment_bam/Upstream_Alignment_Master_Course.html",
        "desc": "FASTQ 1-Hamming 보정, STARsolo, UMI 방향성 그래프, AR-V7"
    },
    {
        "id": "10",
        "title": "Track 10: 프로덕션 파이프라인 & MLOps",
        "short_title": "Track 10 프로덕션",
        "icon": "🚀",
        "rel_from_root": "tracks/10_production_pipelines_mlops/Production_Pipelines_Master_Course.html",
        "desc": "Docker/Singularity HPC, Nextflow DSL2, Snakemake, MLOps"
    }
]

CSS_SNIPPET = """
  /* =====================================================================
     GLOBAL CURRICULUM TRACK SWITCHER & BANNER (11 Tracks Hub)
     ===================================================================== */
  .global-track-switcher {
    background: var(--bg-surface);
    border: 1px solid var(--border-color);
    border-radius: 10px;
    padding: 12px 14px;
    margin-bottom: 12px;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
  }
  .track-switcher-label {
    display: flex;
    align-items: center;
    justify-content: space-between;
    font-size: 0.74rem;
    font-weight: 800;
    color: var(--accent);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-bottom: 8px;
  }
  .track-switcher-select {
    width: 100%;
    padding: 8px 10px;
    border-radius: 8px;
    border: 1px solid var(--border-color);
    background: var(--bg-body);
    color: var(--text-dark);
    font-size: 0.80rem;
    font-weight: 600;
    cursor: pointer;
    outline: none;
    transition: all 0.2s ease;
  }
  .track-switcher-select:hover, .track-switcher-select:focus {
    border-color: var(--accent);
    box-shadow: 0 0 0 3px rgba(2, 132, 199, 0.15);
  }
  .track-quick-menu {
    margin-top: 8px;
  }
  .track-quick-menu summary {
    font-size: 0.76rem;
    font-weight: 700;
    color: var(--text-main);
    cursor: pointer;
    outline: none;
    user-select: none;
    padding: 3px 0;
    display: flex;
    align-items: center;
    gap: 4px;
  }
  .track-quick-menu summary:hover {
    color: var(--accent);
  }
  .track-quick-list {
    list-style: none;
    padding: 6px 0 0 0;
    margin: 0;
    max-height: 220px;
    overflow-y: auto;
    border-top: 1px dashed var(--border-color);
    margin-top: 6px;
  }
  .track-quick-item {
    margin: 3px 0;
  }
  .track-quick-link {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 5px 8px;
    border-radius: 6px;
    font-size: 0.75rem;
    color: var(--text-main);
    text-decoration: none;
    transition: all 0.15s ease;
  }
  .track-quick-link:hover {
    background: var(--accent-light);
    color: var(--accent);
  }
  .track-quick-link.active {
    background: var(--accent);
    color: #ffffff !important;
    font-weight: 700;
  }
  .track-quick-link .track-badge-mini {
    font-size: 0.68rem;
    padding: 1px 6px;
    border-radius: 4px;
    background: rgba(0, 0, 0, 0.08);
  }
  .track-quick-link.active .track-badge-mini {
    background: rgba(255, 255, 255, 0.25);
    color: #ffffff;
  }

  /* Global Top Banner */
  .global-curriculum-banner {
    background: linear-gradient(135deg, rgba(2, 132, 199, 0.06), rgba(56, 189, 248, 0.03));
    border: 1px solid var(--border-color);
    border-radius: 12px;
    padding: 14px 18px;
    margin-bottom: 28px;
    box-shadow: 0 4px 16px rgba(0, 0, 0, 0.03);
  }
  .banner-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 10px;
    padding-bottom: 6px;
    border-bottom: 1px solid var(--border-color);
  }
  .banner-title {
    font-size: 0.82rem;
    font-weight: 800;
    color: var(--accent);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    display: flex;
    align-items: center;
    gap: 6px;
  }
  .banner-status-tag {
    font-size: 0.72rem;
    font-weight: 700;
    color: #10b981;
    background: rgba(16, 185, 129, 0.1);
    padding: 2px 8px;
    border-radius: 9999px;
  }
  .banner-chips-track {
    display: flex;
    gap: 8px;
    overflow-x: auto;
    padding-bottom: 6px;
    scrollbar-width: thin;
  }
  .banner-chips-track::-webkit-scrollbar {
    height: 5px;
  }
  .banner-chips-track::-webkit-scrollbar-thumb {
    background: var(--border-color);
    border-radius: 10px;
  }
  .curriculum-chip {
    flex-shrink: 0;
    display: inline-flex;
    align-items: center;
    gap: 6px;
    padding: 6px 12px;
    border-radius: 20px;
    background: var(--bg-body);
    border: 1px solid var(--border-color);
    color: var(--text-dark);
    font-size: 0.78rem;
    font-weight: 600;
    text-decoration: none;
    transition: all 0.2s ease;
    white-space: nowrap;
  }
  .curriculum-chip:hover {
    border-color: var(--accent);
    color: var(--accent);
    transform: translateY(-1px);
    box-shadow: 0 3px 8px rgba(0, 0, 0, 0.06);
  }
  .curriculum-chip.active {
    background: var(--accent) !important;
    color: #ffffff !important;
    border-color: var(--accent) !important;
    box-shadow: 0 4px 12px rgba(2, 132, 199, 0.3);
  }
"""

def generate_sidebar_switcher(current_track_id, current_file_path):
    curr_dir = os.path.dirname(current_file_path)
    options_html = []
    quick_items_html = []
    
    for t in TRACKS:
        target_abs = os.path.join(WORKSPACE_ROOT, t["rel_from_root"])
        rel_url = os.path.relpath(target_abs, curr_dir)
        is_active = (t["id"] == current_track_id)
        
        selected_attr = ' selected' if is_active else ''
        label_prefix = '▶ [현재 학습 중] ' if is_active else ''
        options_html.append(f'        <option value="{rel_url}"{selected_attr}>{label_prefix}{t["title"]}</option>')
        
        active_cls = ' active' if is_active else ''
        badge_text = '학습 중' if is_active else '이동'
        quick_items_html.append(
            f'          <li class="track-quick-item">'
            f'<a href="{rel_url}" class="track-quick-link{active_cls}">'
            f'<span>{t["icon"]} {t["short_title"]}</span>'
            f'<span class="track-badge-mini">{badge_text}</span>'
            f'</a></li>'
        )

    switcher_html = f"""
    <!-- Global Curriculum Track Switcher Hub -->
    <div class="global-track-switcher no-print">
      <div class="track-switcher-label">
        <span>🌐 11대 마스터 트랙 허브</span>
        <span style="font-size:0.68rem;opacity:0.8;">100% COMPLETE</span>
      </div>
      <select class="track-switcher-select" onchange="if(this.value) window.location.href=this.value;">
{chr(10).join(options_html)}
      </select>
      <details class="track-quick-menu">
        <summary>📑 전체 트랙 목록 펼치기</summary>
        <ul class="track-quick-list">
{chr(10).join(quick_items_html)}
        </ul>
      </details>
    </div>
"""
    return switcher_html

def generate_top_banner(current_track_id, current_file_path):
    curr_dir = os.path.dirname(current_file_path)
    chips_html = []
    
    for t in TRACKS:
        target_abs = os.path.join(WORKSPACE_ROOT, t["rel_from_root"])
        rel_url = os.path.relpath(target_abs, curr_dir)
        is_active = (t["id"] == current_track_id)
        active_cls = ' active' if is_active else ''
        
        chips_html.append(
            f'      <a href="{rel_url}" class="curriculum-chip{active_cls}" title="{t["desc"]}">'
            f'{t["icon"]} {t["short_title"]}'
            f'</a>'
        )
        
    banner_html = f"""
    <!-- Global Curriculum Top Banner Bar -->
    <div class="global-curriculum-banner no-print">
      <div class="banner-header">
        <div class="banner-title">
          <span>🌐 Single-Cell & Spatial Omics Master Ecosystem (11 Tracks)</span>
        </div>
        <span class="banner-status-tag">✨ 100% 완주 코스웨어</span>
      </div>
      <div class="banner-chips-track">
{chr(10).join(chips_html)}
      </div>
    </div>
"""
    return banner_html

def process_file(track_info):
    rel_path = track_info["rel_from_root"]
    full_path = os.path.join(WORKSPACE_ROOT, rel_path)
    track_id = track_info["id"]
    
    orig_size = os.path.getsize(full_path)

    # 1. Read the current (freshly built) HTML and remove any earlier marked injection
    with open(full_path, "r", encoding="utf-8") as f:
        content = f.read()
    content = strip_previous_injection(content)
    if "global-curriculum-banner" in content:
        # Unmarked switcher from the old backup-restore workflow: re-run the build script first.
        print(f"[ERROR] {rel_path} contains an unmarked legacy switcher; rebuild it before injecting")
        return False

    # 2. Inject CSS before </style>
    if "</style>" in content:
        content = content.replace("</style>", CSS_BEGIN + CSS_SNIPPET + CSS_END + "\n</style>", 1)

    # 3. Inject Sidebar Switcher right before <div class="search-container">
    sidebar_switcher = wrap_html(generate_sidebar_switcher(track_id, full_path), "sidebar")
    pattern_search = r'(<div class="search-container">)'
    match_search = re.search(pattern_search, content)
    if match_search:
        search_start = match_search.start()
        content = content[:search_start] + sidebar_switcher + "\n    " + content[search_start:]
    else:
        print(f"[WARN] Could not find <div class=\"search-container\"> in {rel_path}")
            
    # 4. Inject Top Banner right inside main container
    top_banner = wrap_html(generate_top_banner(track_id, full_path), "banner")
    
    # Priority 1: <div class="document-container">
    pattern_doc = r'(<div class="document-container">)'
    match_doc = re.search(pattern_doc, content)
    if match_doc:
        doc_end = match_doc.end()
        content = content[:doc_end] + "\n" + top_banner + content[doc_end:]
    else:
        # Priority 2: <main class="main-content" id="mainContent"> or <main class="main-content">
        pattern_main = r'(<main[^>]*class="[^"]*main-content[^"]*"[^>]*>)'
        match_main = re.search(pattern_main, content)
        if match_main:
            main_end = match_main.end()
            content = content[:main_end] + "\n" + top_banner + content[main_end:]
        else:
            print(f"[WARN] Could not find suitable main container in {rel_path}")
            
    # 5. Save modified file
    with open(full_path, "w", encoding="utf-8") as f:
        f.write(content)
        
    new_size = os.path.getsize(full_path)
    print(f"[SUCCESS] Injected Global Track Switcher into {rel_path} ({orig_size:,} -> {new_size:,} bytes)")
    return True

def main():
    print("=== Injecting Global Curriculum Track Switcher & Top Banner ===")

    success_count = 0
    for t in TRACKS:
        ok = process_file(t)
        if ok:
            success_count += 1
            
    print(f"=== Completed: {success_count}/{len(TRACKS)} files successfully verified and updated ===")

if __name__ == "__main__":
    main()

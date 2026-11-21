#!/usr/bin/env python3
"""
build_trajectory_master.py
Consolidates the 10 Track 06 Single-Cell Manifold Geometry & Dynamic Trajectory Inference
course modules into:
1. Manifold_Trajectory_Master_Book.md (Unified Markdown Textbook)
2. Manifold_Trajectory_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Manifold_Trajectory_Master_Course.docx (Native OpenXML Microsoft Word Document)

Features:
- Zero-dependency offline LaTeX math conversion (HTML stacked fractions + Unicode text for DOCX).
- Fully inlined 6 vector SVG diagrams with zero external HTTP(S) asset dependencies.
- Interactive Study Guide UX: Sticky sidebar, real-time keyword search, dark/light theme toggle, quiz solution toggle.
- Print/PDF optimization with automatic details expansion on beforeprint.
- Automated integrity assertions verifying all 8 quality criteria:
  1. All 10 markdown chapters exist and have even code fence counts.
  2. All 6 SVG diagrams exist and parse as valid XML.
  3. Both slide decks exist and have embedded SVGs.
  4. Both Python pipelines executed and exited with Code 0.
  5. Exactly 30 Cognitive Checkpoint quiz items exist with <details> solutions.
  6. Master HTML contains zero external http/https CDN references.
  7. Master DOCX is valid OpenXML with size > 50 KB.
  8. Master roadmap and README have valid links to Track 06.

CLI (pipeline assertion 4):
  default            re-run both pipelines on a TEMP COPY of pipelines/ (same interpreter/env) and require exit 0 AND that the
                     produced results_pipeline_0X equal the committed ones, except for the fields declared in
                     scripts/pipeline_volatile_fields.json (timing/timestamp/benchmark values, python_version)
  --skip-pipelines   do not run the pipelines (also --verify-fast / VERIFY_FAST=1 / --quick)
  --update-results   run in place and overwrite the committed results (intentional regeneration only)
"""

from __future__ import annotations

import os
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('course_markdown', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'scripts', 'course_markdown.py'))
course_markdown = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_markdown)

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # ENG-025: 스크립트 위치 기준(이식 가능)
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "06_manifold_trajectory")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_TRAJECTORY_ROADMAP.md",
    "01_MANIFOLD_GEOMETRY_AND_DIFFUSION.md",
    "02_RNA_VELOCITY_SPLICING_KINETICS.md",
    "03_DYNAMICAL_VELOCITY_AND_LATENT_TIME.md",
    "04_OPTIMAL_TRANSPORT_WADDINGTON_OT.md",
    "05_CELLRANK2_MULTIKERNEL_INTEGRATION.md",
    "06_MARKOV_CHAINS_AND_MACROSTATES.md",
    "07_FATE_ABSORPTION_PROBABILITIES.md",
    "08_LINEAGE_DRIVER_DISCOVERY.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 단일세포 매니폴드 궤적 분석 및 동역학 로드맵 (Trajectory Roadmap)",
    "01. 매니폴드 기하학과 확산 유사시간 (Manifold Geometry & Diffusion Pseudotime)",
    "02. RNA 벨로시티와 스플라이싱 반응속도론 (RNA Velocity & Splicing Kinetics)",
    "03. 동역학 모델과 보편 잠재 시간 (Dynamical Velocity & Latent Time)",
    "04. 최적 수송 이론과 Waddington-OT (Optimal Transport & Waddington-OT)",
    "05. CellRank 2 다중 커널 융합 프레임워크 (CellRank 2 Multi-Kernel Integration Framework)",
    "06. 마르코프 연쇄와 거시상태의 분해 (Markov Chains & G-PCCA Macrostates)",
    "07. 흡수 마르코프 연쇄와 계통 운명 흡수 확률 (Absorbing Markov Chains & Fate Absorption Probabilities)",
    "08. 계통 운명 구동 전사인자 및 동역학 드라이버 발굴 (Lineage Fate Drivers & Dynamic Gene Discovery)",
    "09. 궤적 분석 문제 해결, 품질 진단 및 캡스톤 루브릭 (Trajectory Troubleshooting, Diagnostics & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 06",
    source_prefix="tracks/06_manifold_trajectory/docs/", track_prefix="tracks/06_manifold_trajectory/",
    md_name="Manifold_Trajectory_Master_Book.md", html_name="Manifold_Trajectory_Master_Course.html", docx_name="Manifold_Trajectory_Master_Course.docx",
    md_header="""# 단일세포 매니폴드 기하학 & 동역학 궤적 추론 (Single-Cell Manifold Geometry & Dynamic Trajectory Inference) 마스터 교재
**From Diffusion Maps & Splicing Kinetics (scVelo) to Entropic Optimal Transport (Waddington-OT) & CellRank 2 Multi-Kernel Lineage Fate Mapping**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 단일세포 매니폴드 궤적 분석 및 동역학 로드맵](#part-00-00-단일세포-매니폴드-궤적-분석-및-동역학-로드맵-trajectory-roadmap)
- [Chapter 01. 매니폴드 기하학과 확산 유사시간](#part-01-01-매니폴드-기하학과-확산-유사시간-manifold-geometry--diffusion-pseudotime)
- [Chapter 02. RNA 벨로시티와 스플라이싱 반응속도론](#part-02-02-rna-벨로시티와-스플라이싱-반응속도론-rna-velocity--splicing-kinetics)
- [Chapter 03. 동역학 모델과 보편 잠재 시간](#part-03-03-동역학-모델과-보편-잠재-시간-dynamical-velocity--latent-time)
- [Chapter 04. 최적 수송 이론과 Waddington-OT](#part-04-04-최적-수송-이론과-waddington-ot-optimal-transport--waddington-ot)
- [Chapter 05. CellRank 2 다중 커널 융합 프레임워크](#part-05-chapter-05-cellrank-2-다중-커널-융합-프레임워크-cellrank-2-multi-kernel-integration-framework)
- [Chapter 06. 마르코프 연쇄와 거시상태의 분해](#part-06-chapter-06-마르코프-연쇄와-거시상태의-분해-markov-chains--g-pcca-macrostates)
- [Chapter 07. 흡수 마르코프 연쇄와 계통 운명 흡수 확률](#part-07-chapter-07-흡수-마르코프-연쇄와-계통-운명-흡수-확률-absorbing-markov-chains--fate-absorption-probabilities)
- [Chapter 08. 계통 운명 구동 전사인자 및 동역학 드라이버 발굴](#part-08-chapter-08-계통-운명-구동-전사인자-및-동역학-드라이버-발굴-lineage-fate-drivers--dynamic-gene-discovery)
- [Chapter 09. 궤적 분석 문제 해결, 품질 진단 및 캡스톤 루브릭](#part-09-chapter-09-궤적-분석-문제-해결-품질-진단-및-캡스톤-루브릭-trajectory-troubleshooting-diagnostics--capstone-rubric)

---

""",
    figcaption_desc="단일세포 매니폴드 기하학 및 동역학 궤적 핵심 아키텍처",
    html_title=f"""단일세포 매니폴드 기하학 &amp; 동역학 궤적 추론 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""Track 06 매니폴드 궤적""",
    html_search_placeholder=f"""전체 교재 검색 (예: scVelo, CellRank, DPT)...""",
    html_theme_key='manifold_trajectory_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>단일세포 매니폴드 기하학 &amp; 동역학 궤적 추론 마스터 교재</h1>
      <p class="subtitle">From Diffusion Maps &amp; Splicing Kinetics (scVelo) to Entropic Optimal Transport (Waddington-OT) &amp; CellRank 2 Multi-Kernel Lineage Fate Mapping</p>
      <div class="hero-badges">
        <span class="hero-badge">10대 모듈 완비</span>
        <span class="hero-badge">6종 전문 SVG 다이어그램 인라인 수록</span>
        <span class="hero-badge">30문항 인지 점검 퀴즈 &amp; 해설</span>
        <span class="hero-badge">100% 무의존성 오프라인 실행</span>
      </div>
    </div>""",
    html_sidebar_badge='',
    html_main_prefix='',
    html_main_suffix='',
    docx_title="단일세포 매니폴드 기하학 & 동역학 궤적 추론 종합 교재",
    docx_subtitle="From Diffusion Maps & Splicing Kinetics (scVelo) to Entropic Optimal Transport (Waddington-OT) & CellRank 2 Multi-Kernel Lineage Fate Mapping",
    docx_meta1="임상 모델: 거세저항성 전립선암(mCRPC) 치료 저항성 및 4대 계통 가소성(AR 선암종, 양방향성 Amphicrine, 이중음성 DNPC, 치료유도 신경내분비암 t-NEPC) 분화 궤적, 확률론적 운명 전이 및 마스터 드라이버 전사인자 발굴",
    docx_meta2="주요 내용: 고전적 선형 유사시간의 위상 왜곡 극복, 확산 맵(Diffusion Maps)과 무방향성 Markov 확산 유사시간(DPT), 전사체 스플라이싱 반응속도론(RNA Velocity) 및 결정론적/동역학(Dynamical scVelo) 모델링과 보편 잠재 시간(Latent Time), 시간 경과 단일 세포 데이터의 엔트로피 정규화 최적 수송(Waddington-OT)과 싱크혼(Sinkhorn) 알고리즘, CellRank 2 다중 커널 융합(속도·유사시간·최적수송·공간), 비가역 마르코프 연쇄의 실 슈어 분해 및 일반화 페론 클러스터 분석(G-PCCA) 거시상태(Macrostates) 분해, 흡수 마르코프 연쇄(Absorbing Markov Chains)와 기본 행렬(Fundamental Matrix)을 통한 분화 운명 흡수 확률(Fate Probabilities), 유전자 발현-운명 상관성 기반 계통 동역학 드라이버(Lineage Drivers) 발굴, 궤적 분석 3대 함정 트러블슈팅 및 캡스톤 프로젝트 루브릭",
    diagram_ref_label="동역학 궤적 시각 다이어그램 참조",
    docx_core_title="단일세포 매니폴드 기하학 &amp; 동역학 궤적 추론 마스터 교재",
    docx_core_subject="Manifold Geometry &amp; Dynamic Trajectory Inference: From Diffusion Maps &amp; scVelo to Waddington-OT &amp; CellRank 2",
    svg_files=[
        "01_diffusion_maps_and_pseudotime_geometry.svg",
        "02_rna_velocity_phase_portrait_kinetics.svg",
        "03_dynamical_scvelo_latent_time_optimization.svg",
        "04_waddington_ot_sinkhorn_transport_coupling.svg",
        "05_cellrank2_multikernel_gpcca_macrostates.svg",
        "06_fate_probabilities_and_lineage_drivers.svg"
    ],
    svg_exact=False,
    roadmap_terms=[
        "Track 06",
        "Manifold_Trajectory_Master_Book.md",
        "Manifold_Trajectory_Master_Course.html",
        "Manifold_Trajectory_Master_Course.docx",
        "tracks/06_manifold_trajectory/docs/",
        "tracks/06_manifold_trajectory/slides/",
        "tracks/06_manifold_trajectory/pipelines/"
    ],
    readme_terms=[
        "Track 06",
        "Manifold_Trajectory_Master_Book.md",
        "Manifold_Trajectory_Master_Course.html",
        "Manifold_Trajectory_Master_Course.docx",
        "tracks/06_manifold_trajectory/"
    ],
)


def read_doc_files() -> list[str]:
    return course_markdown.read_doc_files(CFG)


def build_master_markdown(contents: list[str]) -> str:
    return course_markdown.build_master_markdown(CFG, contents)


def build_master_html(contents: list[str]) -> str:
    return course_markdown.build_master_html(CFG, contents)


def build_master_docx(contents: list[str]) -> str:
    return course_markdown.build_master_docx(CFG, contents)


def _assert_pipelines():
    course_markdown.assert4_summary_files(CFG, "01_diffusion_and_velocity_pipeline.py", "02_cellrank2_fate_mapping_pipeline.py", "pipeline_01_summary.txt", "pipeline_02_summary.txt",
                                          ("  Running Pipeline 01 (Diffusion & Velocity)...", "  Running Pipeline 02 (CellRank 2 Fate Mapping)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 06 Manifold Geometry & Dynamic Trajectory Inference Master Documents...")
    contents = read_doc_files()
    
    # 1. Master Markdown
    build_master_markdown(contents)
    
    # 2. Master HTML (Interactive UX, Search, Dark Mode, Print/PDF)
    html_content = build_master_html(contents)
    
    # 3. Master Word (.docx)
    docx_xml = build_master_docx(contents)
    
    # 4. Automated Verification (All 8 assertions)
    verify_build_artifacts(html_content, docx_xml)

if __name__ == "__main__":
    main()

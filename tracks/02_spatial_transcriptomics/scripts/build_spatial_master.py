#!/usr/bin/env python3
"""
build_spatial_master.py
Consolidates the 10 Track 02 Spatial Transcriptomics course modules into:
1. Spatial_Transcriptomics_Master_Book.md (Unified Markdown Textbook)
2. Spatial_Transcriptomics_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Spatial_Transcriptomics_Master_Course.docx (Native OpenXML Microsoft Word Document)

Features:
- Zero-dependency offline LaTeX math conversion (HTML stacked fractions + Unicode text for DOCX).
- Fully inlined 6 vector SVG diagrams with zero external HTTP(S) asset dependencies.
- Interactive Study Guide UX: Sticky sidebar, real-time keyword search, dark/light theme toggle, quiz solution toggle.
- Print/PDF optimization with automatic details expansion on beforeprint.
- Automated integrity assertions verifying all 8 quality criteria:
  1. All 10 markdown chapters exist and have even code fence counts.
  2. All 6 SVG diagrams exist and parse as valid XML.
  3. Both slide decks exist and have embedded SVGs.
  4. Both Python pipelines exit with code 0.
  5. Exactly 30 Cognitive Checkpoint quiz items exist with <details> solutions.
  6. Master HTML contains zero external http/https CDN references.
  7. Master DOCX is valid OpenXML with size > 50 KB.
  8. Master roadmap and README have valid links to Track 02.

CLI (pipeline assertion 4):
  default            re-run both pipelines on a TEMP COPY of pipelines/ (same interpreter/env) and require exit 0 AND that the
                     produced results_pipeline_0X equal the committed ones, except for the fields declared in
                     scripts/pipeline_volatile_fields.json (timing/timestamp/benchmark values, python_version)
  --skip-pipelines   do not run the pipelines (also --verify-fast / VERIFY_FAST=1)
  --update-results   run in place and overwrite the committed results (intentional regeneration only)
"""

from __future__ import annotations

import os
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('course_markdown', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'scripts', 'course_markdown.py'))
course_markdown = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_markdown)

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # ENG-025: 스크립트 위치 기준(이식 가능)
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "02_spatial_transcriptomics")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_SPATIAL_ROADMAP.md",
    "01_SPATIAL_TECHNOLOGIES_AND_PHYSICS.md",
    "02_HISTOLOGY_AND_COORDINATE_SYSTEMS.md",
    "03_SPATIAL_QC_AND_PREPROCESSING.md",
    "04_CELL_SEGMENTATION_ALGORITHMS.md",
    "05_SPATIAL_DECONVOLUTION_MATHEMATICS.md",
    "06_SPATIALLY_VARIABLE_GENES.md",
    "07_CELLULAR_NICHES_AND_MICROENVIRONMENTS.md",
    "08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. Track 02 로드맵: 전립선암 공간전사체학 (Spatial Transcriptomics) 학습 과정",
    "01. 공간전사체 기술 원리와 광학·생물물리학적 기초 (Spatial Technologies & Biophysics)",
    "02. 조직병리학 영상 정합과 공간 좌표계 체계 (Histology & Coordinate Systems)",
    "03. 공간전사체 품질 관리(QC)와 전처리 기법 (Spatial QC & Preprocessing)",
    "04. 고해상도 세포 분할 알고리즘과 전사체 할당 (Cell Segmentation: Cellpose & Baysor)",
    "05. 스팟 디콘볼루션 수학과 세포 구성 복원 (Spatial Deconvolution Mathematics)",
    "06. 공간 가변 유전자 규명과 통계적 검정 (Spatially Variable Genes, SVGs)",
    "07. 세포 니치 및 공간 미세환경 클러스터링 (Cellular Niches and Microenvironments)",
    "08. 공간 리간드-수용체 상호작용 및 신호 플럭스 모델링 (Spatial Ligand-Receptor Signaling & COMMOT)",
    "09. 트러블슈팅, 품질 제어 및 캡스톤 프로젝트 루브릭 (Troubleshooting, QC & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 02",
    source_prefix="tracks/02_spatial_transcriptomics/docs/", track_prefix="tracks/02_spatial_transcriptomics/",
    md_name="Spatial_Transcriptomics_Master_Book.md", html_name="Spatial_Transcriptomics_Master_Course.html", docx_name="Spatial_Transcriptomics_Master_Course.docx",
    md_header="""# 전립선암 공간전사체학 (Spatial Transcriptomics) 마스터 교재
**From Sequencing-Based Arrays (Visium & Visium HD) to Subcellular Imaging (Xenium & CosMx): Theory, Mathematics, and Clinical Implementation**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. Track 02 로드맵: 전립선암 공간전사체학 학습 과정](#part-00-chapter-00-track-02-로드맵-전립선암-공간전사체학-spatial-transcriptomics-학습-과정)
- [Chapter 01. 공간전사체 기술 원리와 광학·생물물리학적 기초](#part-01-chapter-01-공간전사체-기술-원리와-광학생물물리학적-기초-spatial-technologies--biophysics)
- [Chapter 02. 조직병리학 영상 정합과 공간 좌표계 체계](#part-02-chapter-02-조직병리학-영상-정합과-공간-좌표계-체계-histology--coordinate-systems)
- [Chapter 03. 공간전사체 품질 관리(QC)와 전처리 기법](#part-03-chapter-03-공간전사체-품질-관리qc와-전처리-기법-spatial-qc--preprocessing)
- [Chapter 04. 고해상도 세포 분할 알고리즘과 전사체 할당](#part-04-chapter-04-고해상도-세포-분할-알고리즘과-전사체-할당-cell-segmentation-cellpose--baysor)
- [Chapter 05. 스팟 디콘볼루션 수학과 세포 구성 복원](#part-05-chapter-05-스팟-디콘볼루션-수학과-세포-구성-복원-spatial-deconvolution-mathematics)
- [Chapter 06. 공간 가변 유전자 규명과 통계적 검정](#part-06-chapter-06-공간-가변-유전자-규명과-통계적-검정-spatially-variable-genes-svgs)
- [Chapter 07. 세포 니치 및 공간 미세환경 클러스터링](#part-07-chapter-07-세포-니치-및-공간-미세환경-클러스터링-cellular-niches-and-microenvironments)
- [Chapter 08. 공간 리간드-수용체 상호작용 및 신호 플럭스 모델링](#part-08-chapter-08-공간-리간드-수용체-상호작용-및-신호-플럭스-모델링-spatial-ligand-receptor-signaling--commot)
- [Chapter 09. 트러블슈팅, 품질 제어 및 캡스톤 프로젝트 루브릭](#part-09-chapter-09-트러블슈팅-품질-제어-및-캡스톤-프로젝트-루브릭-troubleshooting-qc--capstone-rubric)

---

""",
    figcaption_desc="전립선암 공간전사체학 핵심 메커니즘 및 아키텍처",
    html_title=f"""전립선암 공간전사체학 (Spatial Transcriptomics) 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""🔬 Spatial Master""",
    html_search_placeholder=f"""교재 전체 실시간 검색...""",
    html_theme_key='spatial_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>전립선암 공간전사체학 (Spatial Transcriptomics) 마스터 교재</h1>
      <p class="subtitle">Single-Cell Dissociation Limits to 2D/3D Tumor Architecture: From Visium HD &amp; Xenium Subcellular Imaging to cell2location Deconvolution and COMMOT Signaling Flux</p>
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
    docx_title="전립선암 공간전사체학 (Spatial Transcriptomics) 종합 교재",
    docx_subtitle="From Sequencing-Based Arrays (Visium & Visium HD) to Subcellular Imaging (Xenium & CosMx): Theory, Mathematics, and Clinical Implementation",
    docx_meta1="임상 모델: 전립선 선암 (Gleason 3/4/5, 체상 Cribriform), CAF 기질 피막 & 면역 배제, mCRPC 골전이",
    docx_meta2="주요 내용: 해리 편향 극복, 아베 회절 한계, SpatialData OME-Zarr 표준, Cellpose vs Baysor 세포 분할, cell2location 베이지안 디콘볼루션, Moran's I 공간 가변 유전자, Squidpy 공간 니치, COMMOT 최적 수송 신호망",
    diagram_ref_label="공간전사체 시각 다이어그램 참조",
    docx_core_title="전립선암 공간전사체학 마스터 교재",
    docx_core_subject="Spatial Transcriptomics: From Visium to Subcellular Imaging &amp; COMMOT",
    svg_files=[
        "01_spatial_technologies_resolution_comparison.svg",
        "02_coordinate_registration_and_spatialdata.svg",
        "03_cellpose_vs_baysor_segmentation.svg",
        "04_cell2location_deconvolution_model.svg",
        "05_morans_i_and_spatial_autocorrelation.svg",
        "06_spatial_niche_and_commot_signaling.svg"
    ],
    svg_exact=True,
    roadmap_terms=[
        "Track 02",
        "Spatial_Transcriptomics_Master_Book.md",
        "Spatial_Transcriptomics_Master_Course.html",
        "Spatial_Transcriptomics_Master_Course.docx",
        "tracks/02_spatial_transcriptomics/docs/",
        "tracks/02_spatial_transcriptomics/slides/",
        "tracks/02_spatial_transcriptomics/pipelines/"
    ],
    readme_terms=[
        "Track 02: 공간 전사체학",
        "Spatial_Transcriptomics_Master_Book.md",
        "Spatial_Transcriptomics_Master_Course.html",
        "Spatial_Transcriptomics_Master_Course.docx",
        "tracks/02_spatial_transcriptomics/"
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
    course_markdown.assert4_result_dirs(CFG, "01_spatial_qc_and_deconvolution_pipeline.py", "02_xenium_subcellular_niche_and_signaling_pipeline.py")


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 02 Spatial Transcriptomics Master Documents...")
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

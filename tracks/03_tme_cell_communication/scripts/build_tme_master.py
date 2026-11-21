#!/usr/bin/env python3
"""
build_tme_master.py
Consolidates the 10 Track 03 Tumor Microenvironment & Cell-Cell Communication course modules into:
1. TME_Cell_Communication_Master_Book.md (Unified Markdown Textbook)
2. TME_Cell_Communication_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. TME_Cell_Communication_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 03.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "03_tme_cell_communication")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_TME_COMMUNICATION_ROADMAP.md",
    "01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md",
    "02_DATABASE_CURATION_AND_CONSENSUS.md",
    "03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md",
    "04_CELLCHAT_MULTISUBUNIT_MODELING.md",
    "05_NICHENET_INTRACELLULAR_SIGNALING.md",
    "06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md",
    "07_MULTIVIEW_MODELING_MISTY.md",
    "08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. Track 03 로드맵: 종양 미세환경(TME) 세포 간 신호 전달 네트워크 (Cell-Cell Communication)",
    "01. 종양 미세환경(TME) 생물학: 세포 사회학과 리간드-수용체 신호 양식 (TME Biology & Signaling Modalities)",
    "02. 리간드-수용체 데이터베이스 큐레이션 및 LIANA 합의 순위 집계 (Database Curation & Consensus Rank)",
    "03. 발현 임계치 설정과 질량작용 법칙 기반 순열 검정 (Expression Thresholding & Permutation Testing)",
    "04. CellChat 다중 소단위체 모델링 및 네트워크 중심성 (CellChat Multi-subunit & Network Centrality)",
    "05. 세포 내 신호 전달 및 표적 유전자 조절망 추론 (NicheNet Causal Downstream Signaling)",
    "06. 공간 전사체 기반 세포 간 통신 통합 분석 (Spatial CCC & LIANA+)",
    "07. 다중 뷰 공간 분해 앙상블 모델링 (Multiview Modeling with MISTy)",
    "08. 전립선암 골전이 미세환경의 5자 신호 전달 회로 (mCRPC Bone Metastasis Signaling Circuits)",
    "09. 트러블슈팅, 품질 제어 및 캡스톤 프로젝트 루브릭 (Troubleshooting, QC & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 03",
    source_prefix="tracks/03_tme_cell_communication/docs/", track_prefix="tracks/03_tme_cell_communication/",
    md_name="TME_Cell_Communication_Master_Book.md", html_name="TME_Cell_Communication_Master_Course.html", docx_name="TME_Cell_Communication_Master_Course.docx",
    md_header="""# 전립선암 종양 미세환경(TME) 세포 간 신호 전달 네트워크 마스터 교재
**From Ligand-Receptor Biophysics and Multi-Subunit Modeling to NicheNet Causal Inference and Prostate Bone Metastasis Circuits**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. Track 03 로드맵: 종양 미세환경(TME) 세포 간 신호 전달 네트워크](#part-00-chapter-00-track-03-로드맵-종양-미세환경tme-세포-간-신호-전달-네트워크-cell-cell-communication)
- [Chapter 01. 종양 미세환경(TME) 생물학: 세포 사회학과 리간드-수용체 신호 양식](#part-01-chapter-01-종양-미세환경tme-생물학-세포-사회학과-리간드-수용체-신호-양식)
- [Chapter 02. 리간드-수용체 데이터베이스 큐레이션 및 LIANA 합의 순위 집계](#part-02-chapter-02-리간드-수용체-데이터베이스-큐레이션-및-liana-합의-순위-집계)
- [Chapter 03. 발현 임계치 설정과 질량작용 법칙 기반 순열 검정](#part-03-chapter-03-발현-임계치-설정과-질량작용-법칙-기반-순열-검정)
- [Chapter 04. CellChat 다중 소단위체 모델링 및 네트워크 중심성](#part-04-chapter-04-cellchat-다중-소단위체-모델링-및-네트워크-중심성)
- [Chapter 05. 세포 내 신호 전달 및 표적 유전자 조절망 추론 (NicheNet)](#part-05-chapter-05-세포-내-신호-전달-및-표적-유전자-조절망-추론-nichenet--causal-signaling-networks)
- [Chapter 06. 공간 전사체 기반 세포 간 통신 통합 분석 (LIANA+)](#part-06-chapter-06-공간-전사체-기반-세포-간-통신-통합-분석-spatial-ccc--liana)
- [Chapter 07. 다중 뷰 공간 분해 앙상블 모델링 (MISTy)](#part-07-chapter-07-다중-뷰-공간-분해-앙상블-모델링-multiview-modeling-with-misty)
- [Chapter 08. 전립선암 골전이 미세환경의 5자 신호 전달 회로](#part-08-chapter-08-전립선암-골전이-미세환경의-5자-신호-전달-회로-mcrpc-bone-metastasis-signaling-circuits)
- [Chapter 09. 트러블슈팅, 품질 제어 및 캡스톤 프로젝트 루브릭](#part-09-chapter-09-세포-간-통신-분석의-함정-검증-프로토콜-및-캡스톤-프로젝트-troubleshooting--capstone-project)

---

""",
    figcaption_desc="전립선암 TME 세포 간 신호 전달 핵심 아키텍처",
    html_title=f"""전립선암 종양 미세환경(TME) 세포 간 신호 전달 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""🧬 TME Master""",
    html_search_placeholder=f"""교재 전체 실시간 검색...""",
    html_theme_key='tme_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>전립선암 종양 미세환경(TME) 세포 간 신호 전달 마스터 교재</h1>
      <p class="subtitle">From Static Single-Cell Profiles to Dynamic Intercellular Communication Networks: CellChat Multi-Subunit Modeling, NicheNet Causal Inference, Spatial CCC &amp; Bone Metastasis Circuits</p>
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
    docx_title="전립선암 종양 미세환경(TME) 세포 간 신호 전달 종합 교재",
    docx_subtitle="From Ligand-Receptor Biophysics and Multi-Subunit Modeling to NicheNet Causal Inference and Prostate Bone Metastasis Circuits",
    docx_meta1="임상 모델: 전립선 선암 및 거세저항성 뼈 전이(mCRPC Bone Metastasis), 골아세포-파골세포 악순환 루프, myCAF ECM 피막 장벽 & 면역 억제",
    docx_meta2="주요 내용: 3대 신호 양식(Juxta/Para/Auto), OmniPath/LIANA 앙상블 Borda 합의 집계, 발현 임계치 & 순열 검정, CellChat 다중 소단위체 기하평균, NicheNet RWR 조절 잠재력 표적인과 추론, LIANA+ 가우시안 거리 감쇄 커널, MISTy 다중 뷰 머신러닝, 골전이 5자 신호망 및 표적 병용 치료 전략",
    diagram_ref_label="세포 간 통신 시각 다이어그램 참조",
    docx_core_title="전립선암 종양 미세환경 세포 간 신호 전달 마스터 교재",
    docx_core_subject="TME Cell-Cell Communication: From CellChat &amp; LIANA to NicheNet &amp; Bone Metastasis Circuits",
    svg_files=[
        "01_tme_cellular_ecosystem_and_signaling_types.svg",
        "02_liana_consensus_and_permutation_testing.svg",
        "03_cellchat_communication_probability_and_centrality.svg",
        "04_nichenet_ligand_to_target_causal_network.svg",
        "05_spatial_misty_multiview_modeling.svg",
        "06_prostate_bone_metastasis_vicious_cycle.svg"
    ],
    svg_exact=True,
    roadmap_terms=[
        "Track 03",
        "TME_Cell_Communication_Master_Book.md",
        "TME_Cell_Communication_Master_Course.html",
        "TME_Cell_Communication_Master_Course.docx",
        "tracks/03_tme_cell_communication/docs/",
        "tracks/03_tme_cell_communication/slides/",
        "tracks/03_tme_cell_communication/pipelines/"
    ],
    readme_terms=[
        "Track 03: 종양 미세환경",
        "TME_Cell_Communication_Master_Book.md",
        "TME_Cell_Communication_Master_Course.html",
        "TME_Cell_Communication_Master_Course.docx",
        "tracks/03_tme_cell_communication/"
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
    course_markdown.assert4_result_dirs(CFG, "01_liana_and_cellchat_consensus_pipeline.py", "02_nichenet_and_spatial_communication_pipeline.py")


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 03 TME & Cell-Cell Communication Master Documents...")
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

#!/usr/bin/env python3
"""
build_statistics_master.py
Consolidates the 10 Track 07 Single-Cell Count Statistics & Differential Abundance
course modules into:
1. Count_Statistics_Master_Book.md (Unified Markdown Textbook)
2. Count_Statistics_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Count_Statistics_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 07.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "07_count_statistics_abundance")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_COUNT_STATISTICS_ROADMAP.md",
    "01_PROBABILITY_MODELS_POISSON_NB.md",
    "02_PSEUDOREPLICATION_PITFALL.md",
    "03_PSEUDOBULK_FRAMEWORK_EDGER_DESEQ2.md",
    "04_GLMM_HIERARCHICAL_MODELING.md",
    "05_FAST_CELL_LEVEL_NEBULA.md",
    "06_TRAJECTORY_DIFFERENTIAL_EXPRESSION.md",
    "07_DIFFERENTIAL_ABUNDANCE_MILO.md",
    "08_COMPOSITIONAL_DATA_ANALYSIS_CODA.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 단일세포 전사체 카운트 통계 및 차등 분석 로드맵 (Count Statistics Roadmap)",
    "01. 단일세포 카운트 확률 모델: 포아송(Poisson)에서 음이항(Negative Binomial)으로 (Probability Models)",
    "02. 가짜 반복(Pseudoreplication)의 통계학적 함정과 재현성 위기 (Pseudoreplication Pitfall)",
    "03. 슈도벌크(Pseudobulk) 분석 프레임워크: edgeR과 DESeq2의 통계학 (Pseudobulk Framework)",
    "04. 위계적 일반화 선형 혼합 모델(GLMM): 세포 수준 공변량과 무작위 효과 (Hierarchical GLMM)",
    "05. 고속 세포 수준 음이항 혼합 모델 (Fast Cell-Level GLMM with NEBULA)",
    "06. 분화 궤적 기반 연속 차등 발현 분석 (Trajectory Differential Expression)",
    "07. 이웃 그래프 기반 차등 풍부도 분석과 Spatial FDR (Differential Abundance with Milo & Spatial FDR)",
    "08. 조성 데이터 분석과 scCODA 베이지안 회귀 (Compositional Data Analysis & scCODA)",
    "09. 트러블슈팅, 실험 설계 및 통계 검증 가이드 (Troubleshooting, Experimental Design & Statistical Benchmarking)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 07",
    source_prefix="tracks/07_count_statistics_abundance/docs/", track_prefix="tracks/07_count_statistics_abundance/",
    md_name="Count_Statistics_Master_Book.md", html_name="Count_Statistics_Master_Course.html", docx_name="Count_Statistics_Master_Course.docx",
    md_header="""# 단일세포 카운트 통계학 및 차등 풍부도 (Single-Cell Count Statistics & Differential Abundance) 마스터 교재
**From Negative Binomial Overdispersion & Pseudoreplication to Pseudobulk (DESeq2/edgeR), Hierarchical GLMM (NEBULA) & Milo Differential Abundance**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 단일세포 전사체 카운트 통계 및 차등 분석 로드맵](#part-00-00-단일세포-전사체-카운트-통계-및-차등-분석-로드맵-count-statistics-roadmap)
- [Chapter 01. 단일세포 카운트 확률 모델: 포아송에서 음이항으로](#part-01-01-단일세포-카운트-확률-모델-포아송poisson에서-음이항negative-binomial으로)
- [Chapter 02. 가짜 반복(Pseudoreplication)의 통계학적 함정과 재현성 위기](#part-02-02-가짜-반복pseudoreplication의-통계학적-함정과-재현성-위기)
- [Chapter 03. 슈도벌크(Pseudobulk) 분석 프레임워크: edgeR과 DESeq2의 통계학](#part-03-03-슈도벌크pseudobulk-분석-프레임워크-edger과-deseq2의-통계학)
- [Chapter 04. 위계적 일반화 선형 혼합 모델(GLMM): 세포 수준 공변량과 무작위 효과](#part-04-04-위계적-일반화-선형-혼합-모델glmm-세포-수준-공변량과-무작위-효과)
- [Chapter 05. 고속 세포 수준 음이항 혼합 모델 (Fast Cell-Level GLMM with NEBULA)](#part-05-05-고속-세포-수준-음이항-혼합-모델-fast-cell-level-glmm-with-nebula)
- [Chapter 06. 분화 궤적 기반 연속 차등 발현 분석 (Trajectory Differential Expression)](#part-06-06-분화-궤적-기반-연속-차등-발현-분석-trajectory-differential-expression--lineage-bifurcation-testing)
- [Chapter 07. 이웃 그래프 기반 차등 풍부도 분석과 Spatial FDR (Differential Abundance with Milo)](#part-07-07-이웃-그래프-기반-차등-풍부도-분석과-spatial-fdr-differential-abundance-with-milo--spatial-fdr)
- [Chapter 08. 조성 데이터 분석과 scCODA 베이지안 회귀 (Compositional Data Analysis & scCODA)](#part-08-08-조성-데이터-분석과-sccoda-베이지안-회귀-compositional-data-analysis--sccoda)
- [Chapter 09. 트러블슈팅, 실험 설계 및 통계 검증 가이드 (Troubleshooting & Best Practices)](#part-09-09-트러블슈팅-실험-설계-및-통계-검증-가이드-troubleshooting-experimental-design--statistical-benchmarking)

---

""",
    figcaption_desc="단일세포 카운트 통계학 및 차등 풍부도 핵심 아키텍처",
    html_title=f"""단일세포 카운트 통계학 &amp; 차등 풍부도 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""Track 07 카운트 통계학""",
    html_search_placeholder=f"""전체 교재 검색 (예: GLMM, Milo, DESeq2)...""",
    html_theme_key='count_statistics_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>단일세포 카운트 통계학 &amp; 차등 풍부도 마스터 교재</h1>
      <p class="subtitle">From Negative Binomial Overdispersion &amp; Pseudoreplication to Pseudobulk (DESeq2/edgeR), Hierarchical GLMM (NEBULA) &amp; Milo Differential Abundance</p>
      <div class="hero-badges">
        <span class="hero-badge">10대 핵심 모듈 완비</span>
        <span class="hero-badge">6종 전문 SVG 다이어그램 인라인 수록</span>
        <span class="hero-badge">30문항 인지 점검 퀴즈 &amp; 해설</span>
        <span class="hero-badge">100% 무의존성 오프라인 실행</span>
      </div>
    </div>""",
    html_sidebar_badge='',
    html_main_prefix='',
    html_main_suffix='',
    docx_title="단일세포 카운트 통계학 및 차등 풍부도 종합 교재",
    docx_subtitle="From Negative Binomial Overdispersion & Pseudoreplication to Pseudobulk (DESeq2/edgeR), Hierarchical GLMM (NEBULA) & Milo Differential Abundance",
    docx_meta1="임상 모델: 거세저항성 전립선암(mCRPC) 환자 코호트의 약물 반응성(치료 반응군 vs 비반응군), 면역 억제 미세환경(탈진 T세포 vs 이펙터 T세포) 및 치료 저항성 신경내분비암(t-NEPC) 분화 궤적 차등 발현 분석",
    docx_meta2="주요 내용: 포아송 분포의 등분산 제약과 단일세포 전사체 과분산(Overdispersion), 음이항(NB) 분포 유도 및 분산 모수(Dispersion θ), 단일세포 단위 통계 검정(Wilcoxon, t-test)에 의한 가짜 반복(Pseudoreplication)의 제1종 오류 대폭발(Type I Error Catastrophe) 실증, 환자 수준 슈도벌크(Pseudobulk) 집계 및 경험적 베이즈(Empirical Bayes) 분산 수축(DESeq2/edgeR LRT 및 QL-F 검정), 세포 수준 공변량과 환자 간 무작위 효과를 동시에 모델링하는 위계적 일반화 선형 혼합 모델(GLMM) 및 고속 NEBULA 엔진, 분화 궤적 잠재시간에 따른 GAM 연속 차등 발현(Trajectory DE), kNN 이웃 그래프 기반 세포 집단 차등 풍부도(Milo DA) 및 Spatial FDR, 조성 데이터(CoDa) Simplex 왜곡 극복과 scCODA 베이지안 회귀 분석, 트러블슈팅 및 캡스톤 프로젝트 루브릭",
    diagram_ref_label="카운트 통계학 시각 다이어그램 참조",
    docx_core_title="단일세포 카운트 통계학 및 차등 풍부도 마스터 교재",
    docx_core_subject="Single-Cell Count Statistics &amp; Differential Abundance: From Negative Binomial Overdispersion &amp; Pseudoreplication to Pseudobulk, Hierarchical GLMM &amp; Milo",
    svg_files=[
        "01_poisson_vs_negative_binomial_overdispersion.svg",
        "02_pseudoreplication_type1_error_catastrophe.svg",
        "03_pseudobulk_aggregation_and_dispersion_shrinkage.svg",
        "04_glmm_hierarchical_random_effects_architecture.svg",
        "05_milo_knn_neighborhood_differential_abundance.svg",
        "06_compositional_simplex_and_coda_regression.svg"
    ],
    svg_exact=False,
    roadmap_terms=[
        "Track 07",
        "Count_Statistics_Master_Book.md",
        "Count_Statistics_Master_Course.html",
        "Count_Statistics_Master_Course.docx",
        "tracks/07_count_statistics_abundance/docs/",
        "tracks/07_count_statistics_abundance/slides/",
        "tracks/07_count_statistics_abundance/pipelines/"
    ],
    readme_terms=[
        "Track 07",
        "Count_Statistics_Master_Book.md",
        "Count_Statistics_Master_Course.html",
        "Count_Statistics_Master_Course.docx",
        "tracks/07_count_statistics_abundance/"
    ],
    read_filter=lambda text: text.replace("](../../0", "](../0"),
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
    course_markdown.assert4_summary_files(CFG, "01_pseudobulk_and_glmm_differential_expression.py", "02_milo_and_coda_differential_abundance_pipeline.py", "pipeline_01_summary.txt", "pipeline_02_summary.txt",
                                          ("  Running Pipeline 01 (Pseudobulk & GLMM)...", "  Running Pipeline 02 (Milo DA & CoDa)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 07 Count Statistics & Differential Abundance Master Documents...")
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

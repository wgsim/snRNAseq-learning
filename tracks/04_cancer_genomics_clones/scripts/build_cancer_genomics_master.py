#!/usr/bin/env python3
"""
build_cancer_genomics_master.py
Consolidates the 10 Track 04 Single-Cell Cancer Genomics & Clonal Evolution course modules into:
1. Cancer_Genomics_Master_Book.md (Unified Markdown Textbook)
2. Cancer_Genomics_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Cancer_Genomics_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 04.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "04_cancer_genomics_clones")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_CANCER_GENOMICS_ROADMAP.md",
    "01_CNA_BIOLOGY_AND_CHROMOSOMAL_ALTERATIONS.md",
    "02_INFERCNV_ALGORITHM_AND_HMM.md",
    "03_REFERENCE_FREE_COPYKAT.md",
    "04_ALLELIC_IMBALANCE_AND_BAF_NUMBAT.md",
    "05_SUBCLONAL_DECONVOLUTION_AND_PHYLOGENY.md",
    "06_SOMATIC_MUTATIONS_AND_FUSIONS.md",
    "07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md",
    "08_METASTASIS_SEEDING_MODELS.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 암 단일세포 전사체 기반 유전체 및 클론 진화 분석 로드맵 (Cancer Single-Cell Genomics & Clonal Evolution Roadmap)",
    "01. 염색체 복제수 변이의 분자생물학 및 전립선암 랜드마크 (CNA Biology & Prostate Cancer Hallmarks)",
    "02. inferCNV 알고리즘과 6상태 은닉 마르코프 모델 (inferCNV Algorithm & 6-State HMM)",
    "03. 참조 비의존성 CopyKAT 알고리즘과 가우시안 혼합 모델 (Reference-Free CopyKAT & GMM)",
    "04. 대립유전자 불균형, BAF, 그리고 Numbat 알고리즘 (Allelic Imbalance, BAF & Numbat)",
    "05. 서브클론 디컨볼루션 및 단일 세포 계통수 복원 (Subclonal Deconvolution & Single-Cell Phylogeny)",
    "06. 체세포 단일염기변이(SNV) 및 융합 전사체 검출 (Somatic Mutations & Gene Fusions in snRNA-seq)",
    "07. 치료 선택압 및 계통 가소성 진화 (Treatment Selection Bottleneck & Lineage Plasticity)",
    "08. 전이 파종 모델 및 다클론성 전파 (Metastasis Seeding Models & Polyclonal Dissemination)",
    "09. 암 유전체 및 클론 진화 분석 트러블슈팅 및 캡스톤 프로젝트 (Troubleshooting & Capstone Project)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 04",
    source_prefix="tracks/04_cancer_genomics_clones/docs/", track_prefix="tracks/04_cancer_genomics_clones/",
    md_name="Cancer_Genomics_Master_Book.md", html_name="Cancer_Genomics_Master_Course.html", docx_name="Cancer_Genomics_Master_Course.docx",
    md_header="""# 단일 세포 암 유전체 & 클론 진화(Single-Cell Cancer Genomics & Clonal Evolution) 마스터 교재
**From Chromosomal Aneuploidy and Copy-Neutral LOH to Clonal Phylogeny and Treatment-Resistant Lineage Plasticity**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 암 단일세포 전사체 기반 유전체 및 클론 진화 분석 로드맵](#part-00-00-암-단일세포-전사체-기반-유전체-및-클론-진화-분석-로드맵-cancer-single-cell-genomics--clonal-evolution-roadmap)
- [Chapter 01. 염색체 복제수 변이의 분자생물학 및 전립선암 랜드마크](#part-01-01-염색체-복제수-변이의-분자생물학-및-전립선암-랜드마크-cna-biology--prostate-cancer-hallmarks)
- [Chapter 02. inferCNV 알고리즘과 6상태 은닉 마르코프 모델](#part-02-02-infercnv-알고리즘과-6상태-은닉-마르코프-모델-infercnv-algorithm--6-state-hmm)
- [Chapter 03. 참조 비의존성 CopyKAT 알고리즘과 가우시안 혼합 모델](#part-03-03-참조-비의존성-copykat-알고리즘과-가우시안-혼합-모델-reference-free-copykat--gmm)
- [Chapter 04. 대립유전자 불균형, BAF, 그리고 Numbat 알고리즘](#part-04-04-대립유전자-불균형-baf-그리고-numbat-알고리즘-allelic-imbalance-baf--numbat)
- [Chapter 05. 서브클론 디컨볼루션 및 단일 세포 계통수 복원](#part-05-chapter-05-서브클론-디컨볼루션-및-단일-세포-계통수-복원-subclonal-deconvolution--single-cell-phylogeny)
- [Chapter 06. 체세포 단일염기변이(SNV) 및 융합 전사체 검출](#part-06-chapter-06-체세포-단일염기변이snv-및-융합-전사체-검출-somatic-mutations--gene-fusions-in-snrna-seq)
- [Chapter 07. 치료 선택압 및 계통 가소성 진화](#part-07-chapter-07-치료-선택압-및-계통-가소성-진화-treatment-selection-bottleneck--lineage-plasticity)
- [Chapter 08. 전이 파종 모델 및 다클론성 전파](#part-08-chapter-08-전이-파종-모델-및-다클론성-전파-metastasis-seeding-models--polyclonal-dissemination)
- [Chapter 09. 암 유전체 및 클론 진화 분석 트러블슈팅 및 캡스톤 프로젝트](#part-09-chapter-09-암-유전체-및-클론-진화-분석-트러블슈팅-및-캡스톤-프로젝트-troubleshooting--capstone-project)

---

""",
    figcaption_desc="단일 세포 암 유전체 &amp; 클론 진화 핵심 아키텍처",
    html_title=f"""단일 세포 암 유전체 &amp; 클론 진화 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""🧬 Cancer Genomics""",
    html_search_placeholder=f"""교재 전체 실시간 검색...""",
    html_theme_key='cancer_genomics_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>단일 세포 암 유전체 &amp; 클론 진화 마스터 교재</h1>
      <p class="subtitle">From Chromosomal Aneuploidy and Copy-Neutral LOH to Single-Cell Clonal Phylogeny, Treatment Selection Bottlenecks &amp; Metastatic Seeding Models</p>
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
    docx_title="단일 세포 암 유전체 & 클론 진화 종합 교재",
    docx_subtitle="From Chromosomal Aneuploidy and Copy-Neutral LOH to Single-Cell Clonal Phylogeny and Treatment-Resistant Lineage Plasticity",
    docx_meta1="임상 모델: 전립선 선암 및 거세저항성 전이암(mCRPC), 8p/8q/10q/13q/17p/Xq 재발 변이, TMPRSS2-ERG 융합, TP53/RB1 결손 매개 t-NEPC 가소성 & 골전이 다클론 파종",
    docx_meta2="주요 내용: 단순 발현 클러스터링의 한계 극복, inferCNV 100유전자 윈도우 이동평균 & 6상태 HMM 비터비 세그멘테이션, CopyKAT 분산안정화 포아송-감마 & 2성분 GMM 참조 비의존 이배체 탐지, Numbat 1000G 위상차 BAF 결합 베이지안 HMM 및 복제중립 LOH (CN-LOH), 단일세포 NJ 및 MEDICC2 클론 계통수 복원, snRNA-seq 리드 기반 체세포 SNV/융합유전자 검출, ADT 진화 병목 및 t-NEPC 가소성 전이, 단일/다클론 골전이 파종 및 캡스톤 프로젝트 루브릭",
    diagram_ref_label="암 유전체 시각 다이어그램 참조",
    docx_core_title="단일 세포 암 유전체 및 클론 진화 마스터 교재",
    docx_core_subject="Cancer Genomics &amp; Clonal Evolution: From inferCNV &amp; CopyKAT to Numbat &amp; Clonal Phylogeny",
    svg_files=[
        "01_prostate_cancer_hallmark_cna_landscape.svg",
        "02_infercnv_hmm_smoothing_workflow.svg",
        "03_copykat_gmm_diploid_aneuploid_separation.svg",
        "04_numbat_haplotype_phased_baf_and_cnloh.svg",
        "05_clonal_phylogeny_tree_and_subclone_evolution.svg",
        "06_therapy_selection_and_metastatic_seeding.svg"
    ],
    svg_exact=True,
    roadmap_terms=[
        "Track 04",
        "Cancer_Genomics_Master_Book.md",
        "Cancer_Genomics_Master_Course.html",
        "Cancer_Genomics_Master_Course.docx",
        "tracks/04_cancer_genomics_clones/docs/",
        "tracks/04_cancer_genomics_clones/slides/",
        "tracks/04_cancer_genomics_clones/pipelines/"
    ],
    readme_terms=[
        "Track 04",
        "Cancer_Genomics_Master_Book.md",
        "Cancer_Genomics_Master_Course.html",
        "Cancer_Genomics_Master_Course.docx",
        "tracks/04_cancer_genomics_clones/"
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
    course_markdown.assert4_result_dirs(CFG, "01_infercnv_and_copykat_aneuploidy_pipeline.py", "02_numbat_baf_and_clonal_phylogeny_pipeline.py")


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 04 Cancer Genomics & Clones Master Documents...")
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

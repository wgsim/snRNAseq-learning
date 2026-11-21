#!/usr/bin/env python3
"""
build_generative_master.py
Consolidates the 10 Track 05 Single-Cell Deep Generative & Foundation Models course modules into:
1. Generative_Foundation_Master_Book.md (Unified Markdown Textbook)
2. Generative_Foundation_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Generative_Foundation_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 05.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "05_generative_foundation_models")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_FOUNDATION_MODELS_ROADMAP.md",
    "01_PROBABILISTIC_COUNT_MODELING.md",
    "02_SCVI_VARIATIONAL_AUTOENCODER.md",
    "03_SCANVI_SEMISUPERVISED_ANNOTATION.md",
    "04_SCARCHES_ARCHITECTURAL_SURGERY.md",
    "05_TRANSFORMER_ARCHITECTURE_SINGLE_CELL.md",
    "06_SCGPT_PRETRAINED_CELL_EMBEDDINGS.md",
    "07_GENEFORMER_AND_GENOME_ATTENTION.md",
    "08_IN_SILICO_PERTURBATION_AND_DRUG_DISCOVERY.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 단일세포 전사체 생성형 파운데이션 모델 및 딥러닝 로드맵 (Generative Foundation Models Roadmap)",
    "01. 단일세포 전사체 확률론적 카운트 모델링 (Probabilistic Count Modeling)",
    "02. scVI 변분 오토인코더 및 확률론적 생성 모델 (scVI Variational Autoencoder)",
    "03. scANVI 준지도 학습 및 희귀 세포형 고정밀 전이 (scANVI Semi-Supervised Annotation)",
    "04. scArches 아키텍처 수술과 지속적 아틀라스 매핑 (scArches Architectural Surgery)",
    "05. 단일 세포 트랜스포머 아키텍처 및 토큰화 메커니즘 (Transformer Architecture & Tokenization)",
    "06. scGPT 사전학습 모델과 세포 임베딩 추출 (scGPT: Pre-trained Foundation Models)",
    "07. Geneformer 아키텍처와 어텐션 기반 유전자 네트워크 분석 (Geneformer: Attention-based Network Biology)",
    "08. 인 실리코 유전자 섭동과 가상 약물 스크리닝 (In Silico Genetic Perturbation & Virtual Drug Discovery)",
    "09. 트러블슈팅, scIB 벤치마킹 및 캡스톤 프로젝트 루브릭 (Troubleshooting, scIB Benchmarking & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 05",
    source_prefix="tracks/05_generative_foundation_models/docs/", track_prefix="tracks/05_generative_foundation_models/",
    md_name="Generative_Foundation_Master_Book.md", html_name="Generative_Foundation_Master_Course.html", docx_name="Generative_Foundation_Master_Course.docx",
    md_header="""# 단일세포 전사체 생성형 파운데이션 모델 & 딥러닝(Single-Cell Generative Deep Learning & Foundation Models) 마스터 교재
**From Probabilistic Count Modeling (NB/ZINB) & scVI Variational Autoencoders to Transformer Foundation Models (scGPT, Geneformer) & In Silico Perturbation Screening**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 단일세포 전사체 생성형 파운데이션 모델 및 딥러닝 로드맵](#part-00-00-단일세포-전사체-생성형-파운데이션-모델-및-딥러닝-로드맵-generative-foundation-models-roadmap)
- [Chapter 01. 단일세포 전사체 확률론적 카운트 모델링](#part-01-01-단일세포-전사체-확률론적-카운트-모델링-probabilistic-count-modeling)
- [Chapter 02. scVI 변분 오토인코더 및 확률론적 생성 모델](#part-02-02-scvi-변분-오토인코더-및-확률론적-생성-모델-scvi-variational-autoencoder)
- [Chapter 03. scANVI 준지도 학습 및 희귀 세포형 고정밀 전이](#part-03-03-scanvi-준지도-학습-및-희귀-세포형-고정밀-전이-scanvi-semi-supervised-annotation)
- [Chapter 04. scArches 아키텍처 수술과 지속적 아틀라스 매핑](#part-04-04-scarches-아키텍처-수술과-지속적-아틀라스-매핑-scarches-architectural-surgery)
- [Chapter 05. 단일 세포 트랜스포머 아키텍처 및 토큰화 메커니즘](#part-05-chapter-05-단일-세포-트랜스포머-아키텍처-및-토큰화-메커니즘-transformer-architecture--tokenization-in-single-cell-biology)
- [Chapter 06. scGPT 사전학습 모델과 세포 임베딩 추출](#part-06-chapter-06-scgpt-사전학습-모델과-세포-임베딩-추출-scgpt-generative-pre-trained-foundation-models-for-single-cell-biology)
- [Chapter 07. Geneformer 아키텍처와 어텐션 기반 유전자 네트워크 분석](#part-07-chapter-07-geneformer-아키텍처와-어텐션-기반-유전자-네트워크-분석-geneformer-genome-scale-pre-training--attention-based-network-biology)
- [Chapter 08. 인 실리코 유전자 섭동과 가상 약물 스크리닝](#part-08-chapter-08-인-실리코-유전자-섭동과-가상-약물-스크리닝-in-silico-genetic-perturbation--virtual-drug-discovery)
- [Chapter 09. 트러블슈팅, scIB 벤치마킹 및 캡스톤 프로젝트 루브릭](#part-09-chapter-09-트러블슈팅-벤치마킹-및-캡스톤-프로젝트-루브릭-troubleshooting-scib-benchmarking--capstone-project-rubric)

---

""",
    figcaption_desc="단일세포 전사체 생성형 파운데이션 모델 핵심 아키텍처",
    html_title=f"""단일세포 전사체 생성형 파운데이션 모델 &amp; 딥러닝 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""Track 05 파운데이션""",
    html_search_placeholder=f"""전체 교재 검색 (예: scVI, ELBO, scGPT)...""",
    html_theme_key='generative_foundation_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>단일세포 전사체 생성형 파운데이션 모델 &amp; 딥러닝 마스터 교재</h1>
      <p class="subtitle">From Probabilistic Count Modeling (NB/ZINB) &amp; scVI Variational Autoencoders to Transformer Foundation Models (scGPT, Geneformer) &amp; In Silico Perturbation Screening</p>
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
    docx_title="단일세포 전사체 생성형 파운데이션 모델 & 딥러닝 종합 교재",
    docx_subtitle="From Probabilistic Count Modeling (NB/ZINB) & scVI Variational Autoencoders to Transformer Foundation Models (scGPT, Geneformer) & In Silico Perturbation Screening",
    docx_meta1="임상 모델: 거세저항성 전립선암(mCRPC) 및 치료유도 신경내분비암(t-NEPC) 계통 가소성, 다기관 코호트 배치 효과, 희귀 아형 고정밀 전이 및 가상 약물 스크리닝",
    docx_meta2="주요 내용: 고전적 log(1+x) 정규화와 선형 배치 보정의 한계 극복, 원시 카운트 음이항(NB/ZINB) 확률 모델링, scVI 변분 오토인코더(VAE)와 재매개변수화 트릭, ELBO 목적 함수 및 KL 어닐링 붕괴 방지, scANVI 준지도 세포형 전이 및 희귀 아형 구출, scArches 조건부 편향 수술을 통한 전신 아틀라스 지속 학습, 트랜스포머 유전자 토큰화 및 어텐션 기반 창발적 유전자 조절망(Emergent GRN), scGPT 다중 작업 사전학습(MGM, GEP), Geneformer 랭크-값 인코딩, 인 실리코 유전자 섭동(Virtual KO)을 통한 t-NEPC 가소성 역전 및 scIB 벤치마킹 체계",
    diagram_ref_label="생성 모델 시각 다이어그램 참조",
    docx_core_title="단일세포 전사체 생성형 파운데이션 모델 마스터 교재",
    docx_core_subject="Generative Foundation Models: From Probabilistic Count Models &amp; scVI to Transformers &amp; In Silico Perturbation",
    svg_files=[
        "01_scvi_vae_probabilistic_architecture.svg",
        "02_elbo_mathematical_optimization.svg",
        "03_scanvi_and_scarches_atlas_mapping.svg",
        "04_single_cell_transformer_tokenization.svg",
        "05_scgpt_multitask_pretraining_and_finetuning.svg",
        "06_in_silico_perturbation_and_lineage_shift.svg"
    ],
    svg_exact=True,
    roadmap_terms=[
        "Track 05",
        "Generative_Foundation_Master_Book.md",
        "Generative_Foundation_Master_Course.html",
        "Generative_Foundation_Master_Course.docx",
        "tracks/05_generative_foundation_models/docs/",
        "tracks/05_generative_foundation_models/slides/",
        "tracks/05_generative_foundation_models/pipelines/"
    ],
    readme_terms=[
        "Track 05",
        "Generative_Foundation_Master_Book.md",
        "Generative_Foundation_Master_Course.html",
        "Generative_Foundation_Master_Course.docx",
        "tracks/05_generative_foundation_models/"
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
    course_markdown.assert4_summary_files(CFG, "01_scvi_and_scanvi_generative_pipeline.py", "02_transformer_embedding_and_perturbation_pipeline.py", "pipeline_summary.json", "pipeline_summary.json",
                                          ("  Running Pipeline 01 (scVI & scANVI)...", "  Running Pipeline 02 (Transformer & Perturbation)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 05 Generative & Foundation Models Master Documents...")
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

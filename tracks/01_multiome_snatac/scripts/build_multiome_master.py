#!/usr/bin/env python3
"""
build_multiome_master.py
Consolidates the 10 Track 01 Multiome (snRNA + snATAC) course modules into:
1. Multiome_snATAC_Master_Book.md (Unified Markdown Textbook)
2. Multiome_snATAC_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Multiome_snATAC_Master_Course.docx (Native OpenXML Microsoft Word Document)

Features:
- Zero-dependency offline LaTeX math conversion (HTML stacked fractions + Unicode text for DOCX).
- Fully inlined 6 vector SVG diagrams with zero external HTTP(S) asset dependencies.
- Interactive Study Guide UX: Sticky sidebar, real-time keyword search, dark/light theme toggle, quiz solution toggle.
- Print/PDF optimization with automatic details expansion on beforeprint.
- Automated integrity assertions verifying 6 inlined SVGs, 0 raw math delimiters, 0 external links, and valid OpenXML.
"""

from __future__ import annotations

import os
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('course_markdown', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'scripts', 'course_markdown.py'))
course_markdown = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_markdown)
import re

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # ENG-025: 스크립트 위치 기준(이식 가능)
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "01_multiome_snatac")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")

DOC_FILES = [
    "00_MULTIOME_ROADMAP.md",
    "01_EPIGENOMICS_AND_TN5_CHEMISTRY.md",
    "02_UPSTREAM_AND_FRAGMENTS.md",
    "03_SNATAC_QC_METRICS.md",
    "04_PEAK_CALLING_AND_LSI.md",
    "05_WNN_MULTIMODAL_INTEGRATION.md",
    "06_PEAK_TO_GENE_AND_ENHANCERS.md",
    "07_TF_FOOTPRINTING_AND_CHROMVAR.md",
    "08_EGRN_AND_SCENICPLUS.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. Track 01 로드맵: 전립선암 10x Multiome 학습 과정 (Roadmap & Environment)",
    "01. 후성유전학 원리와 Tn5 트랜스포세이스 절단-접합 화학 (Tagmentation Chemistry)",
    "02. 업스트림 정렬 및 Fragment 파일 구조: FASTQ에서 게놈 좌표까지 (Cell Ranger ARC & Chromap)",
    "03. snATAC 특화 품질 관리(QC) 지표와 염색질 무결성 판별 (TSS Enrichment & Nucleosome Ladder)",
    "04. 피크 콜링(MACS3)과 LSI 차원 축소의 기하학 (Pseudobulk MACS3 & TF-IDF/SVD)",
    "05. WNN 멀티모달 통합 및 계통 미세 클론 규명 (Weighted Nearest Neighbor)",
    "06. 원거리 인핸서-프로모터 루핑 및 Peak-to-Gene 연계 분석 (Cicero & AR Enhancer)",
    "07. 전사 인자 풋프린팅 및 단일 핵 모티프 편차 분석 (ChromVAR & TOBIAS)",
    "08. 인핸서 주도 전사조절망 복원: SCENIC+ 및 삼중 eRegulon 아키텍처 (Triplet eRegulons)",
    "09. 멀티오믹스 실전 트러블슈팅 및 5단계 캡스톤 프로젝트 루브릭 (Troubleshooting & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    source_prefix="tracks/01_multiome_snatac/docs/", track_prefix="tracks/01_multiome_snatac/",
    md_name="Multiome_snATAC_Master_Book.md", html_name="Multiome_snATAC_Master_Course.html", docx_name="Multiome_snATAC_Master_Course.docx",
    md_header="""# 전립선암 10x Multiome (snRNA + snATAC) & 후성유전학 마스터 교재
**Single-Nucleus Multiome Sequencing for Cancer Epigenomics: From Tn5 Tagmentation to SCENIC+ eGRN**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Part 00. Chapter 00. Track 01 로드맵: 전립선암 10x Multiome (snRNA + snATAC) 학습 과정](#part-00-chapter-00-track-01-로드맵-전립선암-10x-multiome-snrna--snatac-학습-과정)
- [Part 01. Chapter 01. 후성유전학 원리와 Tn5 트랜스포세이스 절단-접합 화학 (Tagmentation Chemistry)](#part-01-chapter-01-후성유전학-원리와-tn5-트랜스포세이스-절단-접합-화학-tagmentation-chemistry)
- [Part 02. Chapter 02. 업스트림 정렬 및 Fragment 파일 구조: FASTQ에서 게놈 좌표까지](#part-02-chapter-02-업스트림-정렬-및-fragment-파일-구조-fastq에서-게놈-좌표까지)
- [Part 03. Chapter 03. snATAC 특화 품질 관리(QC) 지표와 염색질 무결성 판별](#part-03-chapter-03-snatac-특화-품질-관리qc-지표와-염색질-무결성-판별)
- [Part 04. Chapter 04. 피크 콜링(MACS3)과 LSI 차원 축소의 기하학](#part-04-chapter-04-피크-콜링macs3과-lsi-차원-축소의-기하학)
- [Part 05. Chapter 05. WNN 멀티모달 통합 및 계통 미세 클론 규명 (Weighted Nearest Neighbor)](#part-05-chapter-05-wnn-멀티모달-통합-및-계통-미세-클론-규명-weighted-nearest-neighbor)
- [Part 06. Chapter 06. 원거리 인핸서-프로모터 루핑 및 Peak-to-Gene 연계 분석 (Cicero & AR Enhancer)](#part-06-chapter-06-원거리-인핸서-프로모터-루핑-및-peak-to-gene-연계-분석-cicero--ar-enhancer)
- [Part 07. Chapter 07. 전사 인자 풋프린팅 및 단일 핵 모티프 편차 분석 (ChromVAR & TOBIAS)](#part-07-chapter-07-전사-인자-풋프린팅-및-단일-핵-모티프-편차-분석-chromvar--tobias)
- [Part 08. Chapter 08. 인핸서 주도 전사조절망 복원: SCENIC+ 및 삼중 eRegulon 아키텍처](#part-08-chapter-08-인핸서-주도-전사조절망-복원-scenic-및-삼중-eregulon-아키텍처)
- [Part 09. Chapter 09. 멀티오믹스 실전 트러블슈팅 및 5단계 캡스톤 프로젝트 루브릭](#part-09-chapter-09-멀티오믹스-실전-트러블슈팅-및-5단계-캡스톤-프로젝트-루브릭)

---

""",
    figcaption_desc="10x Multiome 생물학적 메커니즘 및 워크플로우",
    html_title=f"""전립선암 10x Multiome & 후성유전학 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""전립선암 10x Multiome 학습 가이드""",
    html_search_placeholder=f"""키워드 검색 (Tn5, LSI, WNN, SCENIC+...)""",
    html_theme_key='prostate_snrna_doc_theme',
    html_hero_html=f"""

      <!-- Cover Page -->
      <div class="cover-page">
        <div class="cover-badge">Track 01 Multiome Masterclass</div>
        <div class="cover-title">전립선암 10x Multiome & 후성유전학<br>(snRNA + snATAC) 종합 교재</div>
        <div class="cover-subtitle">Single-Nucleus Multiome Sequencing for Cancer Epigenomics:<br>From Tn5 Tagmentation to SCENIC+ eGRN</div>
        <div class="cover-meta">
          <p>• <strong>대상 질환 모델</strong>: 거세저항성 전립선암 (mCRPC, 골/간/림프절 전이 및 치료 유도 신경내분비암 t-NEPC)</p>
          <p>• <strong>핵심 학습 영역</strong>: Tn5 트랜스포세이스 절단-접합 화학, Fragment 파일 파싱, TSS Enrichment/뉴클레오솜 래더 QC, MACS3 Pseudobulk 피크셋, LSI TF-IDF/SVD 기하학 및 Component 1 심도 배제, WNN 가중 통합, Cicero 3차원 루핑 및 AR 원거리 인핸서 분석, TF 풋프린팅 & ChromVAR, SCENIC+ 삼중 eRegulon, 골전이 저온 EDTA 탈회 프로토콜</p>
          <p>• <strong>버전 및 형식</strong>: Comprehensive Master Edition (Interactive HTML / Print-Ready PDF / Word DOCX)</p>
        </div>
      </div>

      <!-- Table of Contents -->
      <div class="toc-section">
        <div class="toc-title">📖 교재 목차 (Table of Contents)</div>
        <ul class="toc-list">
          {''.join([f'<li class="toc-item"><a href="#chapter-{i:02d}"><span>{title}</span><span>Chap {i:02d}</span></a></li>' for i, title in enumerate(CHAPTER_TITLES)])}
        </ul>
      </div>

      <!-- Chapter Articles -->""",
    html_sidebar_badge=f"""Multiome Masterclass""",
    html_main_prefix='<div class="document-container">',
    html_main_suffix='</div>',
    docx_title="전립선암 10x Multiome (snRNA + snATAC) 종합 교재",
    docx_subtitle="Single-Nucleus Multiome Sequencing for Cancer Epigenomics: From Tn5 Tagmentation to SCENIC+ eGRN",
    docx_meta1="임상 모델: 거세저항성 전이암 (mCRPC, Bone/Visceral Mets, t-NEPC)",
    docx_meta2="주요 내용: Tn5 화학, Fragment 파일, TSS/Ladder QC, MACS3 피크셋, LSI 1 심도 배제, WNN 통합, Cicero 루핑, AR 원거리 인핸서, TF 풋프린팅 & ChromVAR, SCENIC+ eGRN, 골전이 EDTA 탈회",
    diagram_ref_label="시각 다이어그램 참조",
    docx_core_title="전립선암 10x Multiome &amp; 후성유전학 마스터 교재",
    docx_core_subject="Single-Nucleus Multiome Sequencing for Cancer Epigenomics",
    read_filter=lambda text: re.sub(r'\]\((0[0-9]_[A-Z0-9_]+\.md)\)', r'](docs/\1)', text),
)


def read_doc_files() -> list[str]:
    return course_markdown.read_doc_files(CFG)


def build_master_markdown(contents: list[str]) -> str:
    return course_markdown.build_master_markdown(CFG, contents)


def build_master_html(contents: list[str]) -> str:
    return course_markdown.build_master_html(CFG, contents)


def build_master_docx(contents: list[str]) -> str:
    return course_markdown.build_master_docx(CFG, contents)


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_basic(html_content, docx_xml, 6, "Track 01")


def main():
    print("🚀 Starting compilation of Track 01 Multiome (snRNA + snATAC) Master Documents...")
    contents = read_doc_files()
    
    # 1. Master Markdown
    build_master_markdown(contents)
    
    # 2. Master HTML (Interactive UX, Search, Dark Mode, Print/PDF)
    html_content = build_master_html(contents)
    
    # 3. Master Word (.docx)
    docx_xml = build_master_docx(contents)
    
    # 4. Automated Verification
    verify_build_artifacts(html_content, docx_xml)

if __name__ == "__main__":
    main()

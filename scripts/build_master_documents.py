#!/usr/bin/env python3
"""
build_master_documents.py
Consolidates the 9 snRNA-seq Prostate Cancer course modules into:
1. snRNAseq_Prostate_Master_Book.md (Unified Markdown Textbook)
2. snRNAseq_Prostate_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. snRNAseq_Prostate_Master_Course.docx (Native OpenXML Microsoft Word Document)

Features:
- Zero-dependency offline LaTeX math conversion (HTML stacked fractions + Unicode text for DOCX).
- Fully inlined 8 vector SVG diagrams with zero external HTTP(S) asset dependencies.
- Interactive UX: Sticky sidebar, real-time keyword search, dark/light theme toggle, quiz solution toggle.
- Print/PDF optimization with automatic details expansion on beforeprint.
- Automated integrity assertions verifying 0 raw math delimiters, 0 external links, and valid OpenXML.
"""

from __future__ import annotations

import os
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('course_markdown', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'course_markdown.py'))
course_markdown = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_markdown)
import re

WORKSPACE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # ENG-025: 스크립트 위치 기준(이식 가능)
DOCS_DIR = os.path.join(WORKSPACE, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")

DOC_FILES = [
    "00_CURRICULUM_AND_ROADMAP.md",
    "01_FOUNDATIONS_AND_WET_LAB.md",
    "02_UPSTREAM_PROCESSING.md",
    "03_QC_AND_CLEANING.md",
    "04_CORE_DOWNSTREAM_PIPELINE.md",
    "05_CELL_TYPE_ANNOTATION.md",
    "06_ADVANCED_INTEGRATION_AND_NETWORKS.md",
    "07_MULTIOMICS_AND_FUTURE.md",
    "08_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 커리큘럼 및 학습 로드맵 (Curriculum & Roadmap)",
    "01. snRNA-seq 기초 이론 및 전립선 조직 핵 분리 (Wet-lab)",
    "02. 업스트림 프로세싱 (Pre-mRNA & Cell Ranger Mechanics)",
    "03. 품질 관리 및 데이터 정제 (snRNA QC & Ambient RNA Removal)",
    "04. 핵심 다운스트림 분석 파이프라인 (Dimensionality & Clustering)",
    "05. 전립선암 세포 유형 동정 및 악성세포 감별 (Annotation & inferCNV)",
    "06. 전이암 고급 통합, 계통 가소성 및 전사조절망 (scVI, Trajectory, pySCENIC)",
    "07. 멀티오믹스 및 최신 확장 기술 (10x Multiome & Spatial)",
    "08. 트러블슈팅, Pseudobulk 및 연구 재현성 가이드라인"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=WORKSPACE, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    source_prefix="docs/", track_prefix="",
    md_name="snRNAseq_Prostate_Master_Book.md", html_name="snRNAseq_Prostate_Master_Course.html", docx_name="snRNAseq_Prostate_Master_Course.docx",
    md_header="""# 전립선암 단일핵 전사체 분석 (snRNA-seq) 마스터 교재
**Single-Nucleus RNA Sequencing for Prostate Cancer: From Foundations to Clinical Multiomics**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Part 00. Chapter 00. 커리큘럼 소개 및 전립선암 snRNA-seq 학습 로드맵](#part-00-chapter-00-커리큘럼-소개-및-전립선암-snrna-seq-학습-로드맵)
- [Part 01. Chapter 01. 실험 원리 및 생물학: 전립선암 조직과 snRNA-seq](#part-01-chapter-01-실험-원리-및-생물학-전립선암-조직과-snrna-seq)
- [Part 02. Chapter 02. 업스트림 데이터 처리: 전립선암 FASTQ에서 카운트 매트릭스까지](#part-02-chapter-02-업스트림-데이터-처리-전립선암-fastq에서-카운트-매트릭스까지)
- [Part 03. Chapter 03. 종양 조직 특화 QC 및 데이터 클리닝 (Ambient RNA & Doublets)](#part-03-chapter-03-종양-조직-특화-qc-및-데이터-클리닝-ambient-rna--doublets)
- [Part 04. Chapter 04. 핵심 다운스트림 분석 파이프라인: 수학적 원리와 알고리즘 인지](#part-04-chapter-04-핵심-다운스트림-분석-파이프라인-수학적-원리와-알고리즘-인지)
- [Part 05. Chapter 05. 세포 유형 어노테이션 및 악성 세포 판별 (inferCNV & TME)](#part-05-chapter-05-세포-유형-어노테이션-및-악성-세포-판별-infercnv--tme)
- [Part 06. Chapter 06. 전이암 고급 분석: 배치 통합, 계통 가소성 궤적, 전사조절망 (pySCENIC)](#part-06-chapter-06-전이암-고급-분석-배치-통합-계통-가소성-궤적-전사조절망-pyscenic)
- [Part 07. Chapter 07. 멀티오믹스 및 차세대 공간전사체 연계 (10x Multiome & Spatial)](#part-07-chapter-07-멀티오믹스-및-차세대-공간전사체-연계-10x-multiome--spatial)
- [Part 08. Chapter 08. 실전 트러블슈팅, Pseudobulk 검증 및 베스트 프랙티스](#part-08-chapter-08-실전-트러블슈팅-pseudobulk-검증-및-베스트-프랙티스)

---

""",
    figcaption_desc="상세 생물학적 메커니즘 및 파이프라인 워크플로우",
    html_title=f"""전립선암 snRNA-seq 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""전립선암 단일핵 전사체 학습 가이드""",
    html_search_placeholder=f"""키워드 검색 (KLK3, CellBender, DNPC...)""",
    html_theme_key='prostate_snrna_doc_theme',
    html_hero_html=f"""

      <!-- Cover Page -->
      <div class="cover-page">
        <div class="cover-badge">Specialized Cancer Genomics Masterclass</div>
        <div class="cover-title">전립선암 단일핵 전사체 분석<br>(snRNA-seq) 종합 교재</div>
        <div class="cover-subtitle">Single-Nucleus RNA Sequencing for Prostate Cancer:<br>From Foundations to Clinical Multiomics</div>
        <div class="cover-meta">
          <p>• <strong>대상 질환 모델</strong>: 전립선 선암 (Primary Adenocarcinoma) 및 거세저항성 전이암 (mCRPC, 골/림프절 전이)</p>
          <p>• <strong>핵심 학습 영역</strong>: 동결 조직 핵 분리 프로토콜, Pre-mRNA/Intron 생물학, 엄밀한 품질관리(QC & CellBender), 차원축소 & 클러스터링 인지적 이해, inferCNV 악성 감별, scVI 배치 통합, NEPC 계통 가소성 추적, pySCENIC 전사조절망, 10x Multiome 및 Pseudobulk 임상 검증</p>
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
    html_sidebar_badge=f"""snRNA-seq Masterclass""",
    html_main_prefix='<div class="document-container">',
    html_main_suffix='</div>',
    docx_title="전립선암 단일핵 전사체 분석 (snRNA-seq) 종합 교재",
    docx_subtitle="Single-Nucleus RNA Sequencing for Prostate Cancer: From Foundations to Clinical Multiomics",
    docx_meta1="임상 모델: 전립선 선암 (Primary Adenocarcinoma) & 거세저항성 전이암 (mCRPC, Bone/LN Mets)",
    docx_meta2="주요 내용: 핵 분리 Wet-lab, Pre-mRNA/Intron, 엄밀한 QC, 차원축소·클러스터링 원리, inferCNV, scVI 통합, NEPC 분화 궤적, pySCENIC 조절망, 10x Multiome, Pseudobulk 검증",
    diagram_ref_label="시각 다이어그램 참조",
    docx_core_title="전립선암 snRNA-seq 마스터 교재",
    docx_core_subject="Single-Nucleus RNA Sequencing for Prostate Cancer",
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
    course_markdown.verify_basic(html_content, docx_xml, 8)


def main():
    print("🚀 Starting compilation of snRNA-seq Prostate Master Documents...")
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

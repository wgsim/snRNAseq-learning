#!/usr/bin/env python3
"""
build_alignment_master.py
Consolidates the 10 Track 09 Low-Level Sequence Alignment & BAM Mechanics
course modules into:
1. Upstream_Alignment_Master_Book.md (Unified Markdown Textbook)
2. Upstream_Alignment_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Upstream_Alignment_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 09.

CLI (pipeline assertion 4):
  default            only verify that the pipeline summary artifacts exist and contain the expected text
  --re-run-pipelines re-run both pipelines on a TEMP COPY of pipelines/ (same interpreter/env) and require exit 0 AND that the
                     produced results_pipeline_0X equal the committed ones, except for the fields declared in
                     scripts/pipeline_volatile_fields.json (timing/timestamp/benchmark values, python_version)
  --update-results   (with --re-run-pipelines) run in place and overwrite the committed results (intentional regeneration only)
"""

from __future__ import annotations

import os
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('course_markdown', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'scripts', 'course_markdown.py'))
course_markdown = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_markdown)

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # ENG-025: 스크립트 위치 기준(이식 가능)
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "09_upstream_alignment_bam")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_ALIGNMENT_AND_BAM_ROADMAP.md",
    "01_FASTQ_ANATOMY_AND_BARCODE_WHITELISTS.md",
    "02_STARSOLO_SUFFIX_ARRAY_ALIGNMENT.md",
    "03_ALEVIN_FRY_PSEUDOALIGNMENT_AND_RAD.md",
    "04_UMI_DEDUPLICATION_GRAPH_ALGORITHMS.md",
    "05_BAM_FILE_MECHANICS_AND_SAMTOOLS.md",
    "06_PYSAM_PROGRAMMATIC_ALIGNMENT_MINING.md",
    "07_ALTERNATIVE_SPLICING_AND_ISOFORM_AR_V7.md",
    "08_GENE_FUSION_DETECTION_TMPRSS2_ERG.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 업스트림 정렬 및 BAM 분석 로드맵: 시퀀서 형광 신호부터 단일세포 전사체 카운트 행렬까지 (Upstream Alignment & BAM Roadmap)",
    "01. FASTQ 파일 해부학 및 세포 바코드 화이트리스트 오류 보정 수리 모델 (FASTQ Anatomy & Barcode Whitelist Correction)",
    "02. STARsolo 접미사 배열(Suffix Array) 정렬 알고리즘 및 초고속 스플라이스 정량 (STARsolo Suffix Array Alignment)",
    "03. Alevin-fry 의사정렬(Pseudoalignment) 및 RAD 바이너리 스트리밍 파이프라인 (Alevin-fry & RAD Streaming)",
    "04. UMI 중복 제거 그래프 알고리즘 및 방향성 인접 검정 수리 모델 (UMI Deduplication Graph Algorithms)",
    "05. BAM 파일의 바이너리 구조와 samtools 기반 정렬 데이터 제어 (BAM Binary Mechanics & samtools Alignment Slicing)",
    "06. pysam 기반 프로그래밍 방식 정렬 데이터 마이닝 (pysam Programmatic Alignment Mining)",
    "07. 대안적 스플라이싱 및 전립선암 AR-V7 전사체 검출 (Alternative Splicing & AR-V7 Isoform Mining)",
    "08. 유전자 융합 검출 및 TMPRSS2-ERG 키메라 리드 분석 (Gene Fusion Detection & TMPRSS2-ERG Chimeric Mining)",
    "09. 업스트림 정렬 트러블슈팅 및 베스트 프랙티스 (Troubleshooting & Alignment Best Practices)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 09",
    source_prefix="tracks/09_upstream_alignment_bam/docs/", track_prefix="tracks/09_upstream_alignment_bam/",
    md_name="Upstream_Alignment_Master_Book.md", html_name="Upstream_Alignment_Master_Course.html", docx_name="Upstream_Alignment_Master_Course.docx",
    md_header="""# 저수준 시퀀싱 정렬, 스플라이스 정션 및 BAM 조작론 (Low-Level Sequence Alignment & BAM Mechanics) 마스터 교재
**Track 09: From Raw FASTQ Anatomy & Whitelist Error Correction to Suffix Array Alignment (STARsolo), RAD Streaming (Alevin-fry), UMI Graph Deduplication, BAM Binary Mechanics, and Single-Cell Isoform/Fusion Mining (AR-V7 & TMPRSS2-ERG)**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 업스트림 정렬 및 BAM 분석 로드맵: 시퀀서 형광 신호부터 단일세포 전사체 카운트 행렬까지](#part-00-00-업스트림-정렬-및-bam-분석-로드맵-시퀀서-형광-신호부터-단일세포-전사체-카운트-행렬까지-upstream-alignment--bam-roadmap)
- [Chapter 01. FASTQ 파일 해부학 및 세포 바코드 화이트리스트 오류 보정 수리 모델](#part-01-01-fastq-파일-해부학-및-세포-바코드-화이트리스트-오류-보정-수리-모델-fastq-anatomy--barcode-whitelist-correction)
- [Chapter 02. STARsolo 접미사 배열(Suffix Array) 정렬 알고리즘 및 초고속 스플라이스 정량](#part-02-02-starsolo-접미사-배열suffix-array-정렬-알고리즘-및-초고속-스플라이스-정량-starsolo-suffix-array-alignment)
- [Chapter 03. Alevin-fry 의사정렬(Pseudoalignment) 및 RAD 바이너리 스트리밍 파이프라인](#part-03-03-alevin-fry-의사정렬pseudoalignment-및-rad-바이너리-스트리밍-파이프라인-alevin-fry--rad-streaming)
- [Chapter 04. UMI 중복 제거 그래프 알고리즘 및 방향성 인접 검정 수리 모델](#part-04-04-umi-중복-제거-그래프-알고리즘-및-방향성-인접-검정-수리-모델-umi-deduplication-graph-algorithms)
- [Chapter 05. BAM 파일의 바이너리 구조와 samtools 기반 정렬 데이터 제어](#part-05-chapter-05-bam-파일의-바이너리-구조와-samtools-기반-정렬-데이터-제어-bam-binary-mechanics--samtools-alignment-slicing)
- [Chapter 06. pysam 기반 프로그래밍 방식 정렬 데이터 마이닝](#part-06-chapter-06-pysam-기반-프로그래밍-방식-정렬-데이터-마이닝-pysam-programmatic-alignment-mining)
- [Chapter 07. 대안적 스플라이싱 및 전립선암 *AR-V7* 전사체 검출](#part-07-chapter-07-대안적-스플라이싱-및-전립선암-ar-v7-전사체-검출-alternative-splicing--ar-v7-isoform-mining)
- [Chapter 08. 유전자 융합 검출 및 *TMPRSS2-ERG* 키메라 리드 분석](#part-08-chapter-08-유전자-융합-검출-및-tmprss2-erg-키메라-리드-분석-gene-fusion-detection--tmprss2-erg-chimeric-mining)
- [Chapter 09. 업스트림 정렬 트러블슈팅 및 베스트 프랙티스](#part-09-chapter-09-업스트림-정렬-트러블슈팅-및-베스트-프랙티스-troubleshooting--alignment-best-practices)

---

""",
    figcaption_desc="저수준 서열 정렬, 스플라이스 정션 및 BAM 조작론 핵심 아키텍처",
    html_title=f"""저수준 시퀀싱 정렬, 스플라이스 정션 및 BAM 조작론 마스터 교재 (Track 09: Low-Level Sequence Alignment &amp; BAM Mechanics)""",
    html_sidebar_title=f"""Track 09 업스트림 정렬 &amp; BAM""",
    html_search_placeholder=f"""전체 교재 검색 (예: FASTQ, STARsolo, Alevin-fry, UMI, CIGAR, AR-V7, TMPRSS2)...""",
    html_theme_key='alignment_bam_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>저수준 시퀀싱 정렬, 스플라이스 정션 및 BAM 조작론 마스터 코스북</h1>
      <p class="subtitle">From Raw FASTQ &amp; Whitelist Error Correction to Suffix Array Alignment (STARsolo), RAD Streaming (Alevin-fry), UMI Graph Deduplication, BAM Binary Mechanics, and Single-Cell Isoform/Fusion Mining (AR-V7 &amp; TMPRSS2-ERG)</p>
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
    docx_title="저수준 시퀀싱 정렬, 스플라이스 정션 및 BAM 조작론 마스터 교재",
    docx_subtitle="Track 09: Low-Level Sequence Alignment & BAM Mechanics Master Course (STARsolo, Alevin-fry, pysam, AR-V7, TMPRSS2-ERG)",
    docx_meta1="임상 & 유전체 모델: 단일세포 전사체 시퀀싱(snRNA-seq/scRNA-seq) 업스트림 프로세싱 및 전립선암(CRPC) Enzalutamide 치료 저항성 세포의 대안적 스플라이싱(AR-V7)과 유전자 융합(TMPRSS2-ERG) 고해상도 마이닝",
    docx_meta2="주요 내용: 차세대 염기서열 분석기(NGS) 광학 형광 신호 및 원시 BCL/FASTQ 품질 점수(Phred Q), 737K 세포 바코드 화이트리스트 1-Hamming 거리 오류 보정 해시 알고리즘, UMI 중복 제거 그래프 수리 모델(방향성 인접 검정 N_parent >= 2*N_child - 1), STARsolo 접미사 배열(Suffix Array) 최대 사상 시드(MMP) 인트론 정렬, Alevin-fry 의사정렬(Pseudoalignment) 및 RAD 바이너리 청크 스트리밍 기법, BGZF 블록 가상 오프셋, CIGAR 문자열, SAM 플래그 비트와 옵션 태그(CB, UB, NH)를 포함한 BAM 파일 바이너리 구조 분해, pysam 프로그래밍 정렬 데이터 마이닝, 전립선암 AR-V7 잠재 엑손 3b 스플라이싱 정션(PSI) 및 TMPRSS2-ERG 21q22 결손 기반 키메라 분할 리드 단일세포 해상도 검출, 업스트림 정렬 트러블슈팅 및 프로덕션 베스트 프랙티스",
    diagram_ref_label="업스트림 정렬 시각 다이어그램 참조",
    docx_core_title="저수준 시퀀싱 정렬, 스플라이스 정션 및 BAM 조작론 마스터 교재",
    docx_core_subject="Track 09: Low-Level Sequence Alignment &amp; BAM Mechanics Master Course",
    svg_files=[
        "01_fastq_read_structure_and_barcode_correction.svg",
        "02_starsolo_suffix_array_spliced_alignment.svg",
        "03_alevin_fry_rad_file_and_splici_index.svg",
        "04_umi_directional_adjacency_graph_deduplication.svg",
        "05_bam_binary_record_cigar_and_tags.svg",
        "06_ar_v7_splice_junction_and_tmprss2_erg_fusion.svg"
    ],
    svg_exact=False,
    roadmap_terms=[
        "Track 09",
        "Upstream_Alignment_Master_Book.md",
        "Upstream_Alignment_Master_Course.html",
        "Upstream_Alignment_Master_Course.docx",
        "tracks/09_upstream_alignment_bam/docs/",
        "tracks/09_upstream_alignment_bam/slides/",
        "tracks/09_upstream_alignment_bam/pipelines/"
    ],
    readme_terms=[
        "Track 09",
        "Upstream_Alignment_Master_Book.md",
        "Upstream_Alignment_Master_Course.html",
        "Upstream_Alignment_Master_Course.docx",
        "tracks/09_upstream_alignment_bam/"
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


def _check_pipeline_reports(p1_text, p2_text):
    assert "Directional Adjacency" in p1_text or "Deduplicated UMIs" in p1_text, f"Pipeline 1 report did not record valid metrics!"
    assert "AR-V7" in p2_text and "TMPRSS2-ERG" in p2_text, f"Pipeline 2 report did not record valid metrics!"


def _assert_pipelines():
    course_markdown.assert4_summary_report(CFG, "01_fastq_demux_and_umi_dedup_pipeline.py", "02_bam_junction_mining_and_isoform_pipeline.py", "pipeline_01_summary_report.md", "pipeline_02_summary_report.md",
                                           _check_pipeline_reports, ("  Re-running Pipeline 01 (FASTQ Demux & UMI Dedup)...", "  Re-running Pipeline 02 (BAM Junction & Isoform Mining)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 09 Low-Level Sequence Alignment & BAM Mechanics Master Documents...")
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

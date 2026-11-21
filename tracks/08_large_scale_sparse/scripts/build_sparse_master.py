#!/usr/bin/env python3
"""
build_sparse_master.py
Consolidates the 10 Track 08 Large-Scale Sparse Computation & GPU Acceleration
course modules into:
1. Large_Scale_Sparse_Master_Book.md (Unified Markdown Textbook)
2. Large_Scale_Sparse_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Large_Scale_Sparse_Master_Course.docx (Native OpenXML Microsoft Word Document)

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
  8. Master roadmap and README have valid links to Track 08.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "08_large_scale_sparse")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_LARGE_SCALE_COMPUTATION_ROADMAP.md",
    "01_SPARSE_MATRIX_MEMORY_CSR_CSC.md",
    "02_OUT_OF_CORE_HDF5_ANN_DATA_BACKED.md",
    "03_ZARR_ARRAY_CHUNKING_AND_COMPRESSION.md",
    "04_TILEDB_SOMA_SINGLE_CELL_DATABASE.md",
    "05_DASK_LAZY_EVALUATION_DISTRIBUTED_QC.md",
    "06_INCREMENTAL_PCA_RANDOMIZED_SVD.md",
    "07_GPU_ACCELERATION_RAPIDS_SINGLECELL.md",
    "08_GPU_GRAPH_COMMUNITY_DETECTION.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 대규모 단일세포 전사체 희소 행렬 연산 및 저장 아키텍처 로드맵 (Large-Scale Single-Cell Sparse Matrix Computation & Storage Architecture Roadmap)",
    "01. 희소 행렬 메모리 레이아웃: CSR, CSC와 CPU 캐시 계층 아키텍처 (Sparse Matrix Memory Layout: CSR vs. CSC and CPU Cache Hierarchy)",
    "02. Out-of-Core 연산과 HDF5 AnnData Backed 모드의 내부 구조 (Out-of-Core Computation and Internal Architecture of HDF5 AnnData Backed Mode)",
    "03. Zarr 클라우드 네이티브 청킹과 Blosc 병렬 압축 아키텍처 (Zarr Cloud-Native Chunking and Blosc Parallel Compression Architecture)",
    "04. TileDB-SOMA 단일세포 데이터베이스와 공간 충전 곡선 아키텍처 (TileDB-SOMA Single-Cell Database and Space-Filling Curve Architecture)",
    "05. Dask 지연 연산(Lazy Evaluation)과 분산 품질 관리(Distributed Out-of-Core QC)",
    "06. 점진적 주성분 분석(Incremental PCA)과 무작위 SVD(Randomized SVD)",
    "07. GPU 가속 단일세포 전처리: RAPIDS 및 rapids-singlecell 아키텍처",
    "08. GPU 가속 그래프 구축 및 커뮤니티 탐지: CAGRA와 cuGraph Leiden",
    "09. 대규모 단일세포 컴퓨팅 문제 해결과 모범 사례 (Troubleshooting & Best Practices)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 08",
    source_prefix="tracks/08_large_scale_sparse/docs/", track_prefix="tracks/08_large_scale_sparse/",
    md_name="Large_Scale_Sparse_Master_Book.md", html_name="Large_Scale_Sparse_Master_Course.html", docx_name="Large_Scale_Sparse_Master_Course.docx",
    md_header="""# 대규모 단일세포 아틀라스 희소 행렬 연산 및 GPU 가속 (Large-Scale Sparse Computation & GPU Acceleration) 마스터 교재
**From Sparse Matrix Memory Layout (CSR/CSC) & Out-of-Core Streaming (HDF5/Zarr/TileDB-SOMA) to Dask Distributed QC, Incremental PCA & RAPIDS GPU Acceleration (cuGraph Leiden)**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 대규모 단일세포 전사체 희소 행렬 연산 및 저장 아키텍처 로드맵](#part-00-00-대규모-단일세포-전사체-희소-행렬-연산-및-저장-아키텍처-로드맵-large-scale-single-cell-sparse-matrix-computation--storage-architecture-roadmap)
- [Chapter 01. 희소 행렬 메모리 레이아웃: CSR, CSC와 CPU 캐시 계층 아키텍처](#part-01-01-희소-행렬-메모리-레이아웃-csr-csc와-cpu-캐시-계층-아키텍처-sparse-matrix-memory-layout-csr-vs-csc-and-cpu-cache-hierarchy)
- [Chapter 02. Out-of-Core 연산과 HDF5 AnnData Backed 모드의 내부 구조](#part-02-02-out-of-core-연산과-hdf5-anndata-backed-모드의-내부-구조-out-of-core-computation-and-internal-architecture-of-hdf5-anndata-backed-mode)
- [Chapter 03. Zarr 클라우드 네이티브 청킹과 Blosc 병렬 압축 아키텍처](#part-03-03-zarr-클라우드-네이티브-청킹과-blosc-병렬-압축-아키텍처-zarr-cloud-native-chunking-and-blosc-parallel-compression-architecture)
- [Chapter 04. TileDB-SOMA 단일세포 데이터베이스와 공간 충전 곡선 아키텍처](#part-04-04-tiledb-soma-단일세포-데이터베이스와-공간-충전-곡선-아키텍처-tiledb-soma-single-cell-database-and-space-filling-curve-architecture)
- [Chapter 05. Dask 지연 연산(Lazy Evaluation)과 분산 품질 관리(Distributed Out-of-Core QC)](#part-05-05-dask-지연-연산lazy-evaluation과-분산-품질-관리distributed-out-of-core-qc)
- [Chapter 06. 점진적 주성분 분석(Incremental PCA)과 무작위 SVD(Randomized SVD)](#part-06-06-점진적-주성분-분석incremental-pca과-무작위-svdrandomized-svd)
- [Chapter 07. GPU 가속 단일세포 전처리: RAPIDS 및 rapids-singlecell 아키텍처](#part-07-07-gpu-가속-단일세포-전처리-rapids-및-rapids-singlecell-아키텍처)
- [Chapter 08. GPU 가속 그래프 구축 및 커뮤니티 탐지: CAGRA와 cuGraph Leiden](#part-08-08-gpu-가속-그래프-구축-및-커뮤니티-탐지-cagra와-cugraph-leiden)
- [Chapter 09. 대규모 단일세포 컴퓨팅 문제 해결과 모범 사례](#part-09-09-대규모-단일세포-컴퓨팅-문제-해결과-모범-사례-troubleshooting--best-practices)

---

""",
    figcaption_desc="대규모 단일세포 희소 행렬 연산 및 GPU 가속 핵심 아키텍처",
    html_title=f"""대규모 단일세포 아틀라스 희소 행렬 연산 &amp; GPU 가속 마스터 교재 (Master Course)""",
    html_sidebar_title=f"""Track 08 대규모 희소 연산 &amp; GPU""",
    html_search_placeholder=f"""전체 교재 검색 (예: CSR, Zarr, TileDB, Dask, cuGraph)...""",
    html_theme_key='sparse_computation_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>대규모 단일세포 아틀라스 희소 행렬 연산 &amp; GPU 가속 마스터 교재</h1>
      <p class="subtitle">From Sparse Matrix Memory Layout (CSR/CSC) &amp; Out-of-Core Streaming (HDF5/Zarr/TileDB-SOMA) to Dask Distributed QC, Incremental PCA &amp; RAPIDS GPU Acceleration (cuGraph Leiden)</p>
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
    docx_title="대규모 단일세포 아틀라스 희소 행렬 연산 및 GPU 가속 종합 교재",
    docx_subtitle="From Sparse Matrix Memory Layout (CSR/CSC) & Out-of-Core Streaming (HDF5/Zarr/TileDB-SOMA) to Dask Distributed QC, Incremental PCA & RAPIDS GPU Acceleration (cuGraph Leiden)",
    docx_meta1="임상 & 빅데이터 모델: 100만~1,000만 단일세포 전사체 아틀라스(Human Cell Atlas, 전립선암 원발암/전이암 대규모 코호트)의 메모리 절약형 고성능 분산·가속 컴퓨팅",
    docx_meta2="주요 내용: 희소 행렬 포인터 압축 레이아웃(CSR, CSC)의 CPU 캐시 미스 및 SpMV(희소 행렬-벡터 곱) 계산 복잡도 O(nnz), HDF5 단일 프로세스 파일 잠금 한계 및 AnnData Backed 모드 스트리밍, Zarr 클라우드 네이티브 다차원 청킹(Chunking) 및 Blosc 다중 스레드 병렬 압축, TileDB-SOMA 단일세포 희소 텐서 데이터베이스 및 힐베르트 공간 충전 곡선(Hilbert SFC) 제로 카피 쿼리, Dask 지연 연산 방향성 비순환 그래프(DAG) 기반 Out-of-Core 분산 QC 및 청크 단위 통계량 집계, 스트리밍 증분형 주성분 분석(Incremental PCA) 및 Halko 무작위 SVD(Randomized SVD)의 저순위 부분공간 복원, NVIDIA RAPIDS 생태계(cuML, cuGraph, rapids-singlecell) 호스트-디바이스 PCIe 스트리밍, CAGRA 근사 kNN 그래프 및 GPU 가속 Leiden 커뮤니티 탐지(50~100배 속도 향상), 트러블슈팅 및 프로덕션 MLOps 벤치마크 루브릭",
    diagram_ref_label="희소 연산 시각 다이어그램 참조",
    docx_core_title="대규모 단일세포 아틀라스 희소 행렬 연산 및 GPU 가속 마스터 교재",
    docx_core_subject="Track 08: Large-Scale Sparse Computation &amp; GPU Acceleration Master Course",
    svg_files=[
        "01_sparse_matrix_csr_csc_memory_layout.svg",
        "02_hdf5_backed_mode_vs_zarr_chunking.svg",
        "03_tiledb_soma_sparse_tensor_architecture.svg",
        "04_dask_lazy_dag_out_of_core_normalization.svg",
        "05_rapids_gpu_acceleration_cuda_pipeline.svg",
        "06_cagra_gpu_knn_and_cugraph_leiden.svg"
    ],
    svg_exact=False,
    roadmap_terms=[
        "Track 08",
        "Large_Scale_Sparse_Master_Book.md",
        "Large_Scale_Sparse_Master_Course.html",
        "Large_Scale_Sparse_Master_Course.docx",
        "tracks/08_large_scale_sparse/docs/",
        "tracks/08_large_scale_sparse/slides/",
        "tracks/08_large_scale_sparse/pipelines/"
    ],
    readme_terms=[
        "Track 08",
        "Large_Scale_Sparse_Master_Book.md",
        "Large_Scale_Sparse_Master_Course.html",
        "Large_Scale_Sparse_Master_Course.docx",
        "tracks/08_large_scale_sparse/"
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
    assert "PASSED" in p1_text, f"Pipeline 1 summary did not record PASSED status!"
    assert "PASSED" in p2_text, f"Pipeline 2 summary did not record PASSED status!"


def _assert_pipelines():
    course_markdown.assert4_summary_report(CFG, "01_sparse_and_out_of_core_pipeline.py", "02_tiledb_soma_and_gpu_acceleration_pipeline.py", "pipeline_01_summary.txt", "pipeline_02_summary.txt",
                                           _check_pipeline_reports, ("  Re-running Pipeline 01 (Sparse & Out-of-Core)...", "  Re-running Pipeline 02 (TileDB & GPU Acceleration)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 08 Large-Scale Sparse & GPU Acceleration Master Documents...")
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

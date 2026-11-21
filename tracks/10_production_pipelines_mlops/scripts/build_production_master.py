#!/usr/bin/env python3
"""
build_production_master.py
Consolidates the 10 Track 10 Scalable Production Pipelines, HPC/Cloud Orchestration & MLOps
course modules into:
1. Production_Pipelines_Master_Book.md (Unified Markdown Textbook)
2. Production_Pipelines_Master_Course.html (Interactive Master Course & Print/PDF-Ready HTML)
3. Production_Pipelines_Master_Course.docx (Native OpenXML Microsoft Word Document)

Features:
- Zero-dependency offline LaTeX math conversion (HTML stacked fractions + Unicode text for DOCX).
- Fully inlined 6 vector SVG diagrams with zero external HTTP(S) asset dependencies.
- Interactive Study Guide UX: Sticky sidebar, real-time keyword search, dark/light theme toggle, quiz solution toggle.
- Print/PDF optimization with automatic details expansion on beforeprint.
- Automated integrity assertions verifying all 8 quality criteria:
  1. All 10 markdown chapters exist and have even code fence counts.
  2. All 6 SVG diagrams exist and parse as valid XML.
  3. Both slide decks exist and have embedded SVGs.
  4. Both Python pipelines executed and exited with Code 0 (verified by checking results summary reports).
  5. Exactly 30 Cognitive Checkpoint quiz items exist with <details> solutions.
  6. Master HTML contains zero external http/https CDN references.
  7. Master DOCX is valid OpenXML with size > 50 KB.
  8. Master roadmap and README have valid links to Track 10.

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
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "10_production_pipelines_mlops")
DOCS_DIR = os.path.join(TRACK_DIR, "docs")
IMAGES_DIR = os.path.join(DOCS_DIR, "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
PIPELINES_DIR = os.path.join(TRACK_DIR, "pipelines")

DOC_FILES = [
    "00_PRODUCTION_PIPELINES_ROADMAP.md",
    "01_CONTAINERIZATION_DOCKER_SINGULARITY.md",
    "02_NEXTFLOW_DSL2_MODULAR_ARCHITECTURE.md",
    "03_SNAKEMAKE_RULE_BASED_EXECUTION.md",
    "04_HPC_SLURM_AND_CLOUD_ORCHESTRATION.md",
    "05_NF_CORE_SCRNASEQ_END_TO_END.md",
    "06_DATA_PROVENANCE_AND_REPRODUCIBILITY.md",
    "07_MLOPS_MODEL_SERVING_SINGLE_CELL.md",
    "08_COST_OPTIMIZATION_AND_RESOURCE_SCALING.md",
    "09_TROUBLESHOOTING_AND_BEST_PRACTICES.md"
]

CHAPTER_TITLES = [
    "00. 프로덕션 파이프라인 및 MLOps 로드맵: 재현성 위기 극복부터 엔터프라이즈 임상 분석까지 (Production Pipelines & MLOps Roadmap)",
    "01. 컨테이너 가상화: Docker vs. Singularity/Apptainer HPC 아키텍처 (Containerization for Production Bioinformatics)",
    "02. Nextflow DSL2 모듈형 아키텍처와 리액티브 데이터플로우 (Nextflow DSL2 Reactive Architecture)",
    "03. Snakemake 규칙 기반 실행 및 타깃 중심 DAG 아키텍처 (Snakemake Rule-Based Workflow Architecture)",
    "04. HPC SLURM 클러스터 및 퍼블릭 클라우드 분산 오케스트레이션 (HPC SLURM & Cloud Batch Orchestration)",
    "05. nf-core/scrnaseq 엔드투엔드 파이프라인 분석 및 커스텀 모듈 엔지니어링 (End-to-End nf-core/scrnaseq Pipeline & Custom Module Engineering)",
    "06. 데이터 출처 추적성(Provenance), FAIR 원칙 및 암호학적 재현성 보증 (Data Provenance, FAIR Principles & Cryptographic Reproducibility)",
    "07. 단일세포 전사체 MLOps 모델 서빙 및 실시간 추론 아키텍처 (Single-Cell MLOps Model Serving & Real-Time Inference Architecture)",
    "08. 클라우드 비용 경제학, 스팟 인스턴스 내결함성 및 대규모 자원 최적화 (Cloud Cost Economics, Spot Fault-Tolerance & Resource Optimization)",
    "09. 프로덕션 파이프라인 트러블슈팅, 엔터프라이즈 모범 사례 및 캡스톤 프로젝트 루브릭 (Production Troubleshooting, Best Practices & Capstone Rubric)"
]

# Track configuration; all conversion logic lives in scripts/course_markdown.py
CFG = course_markdown.TrackConfig(
    workspace=WORKSPACE, out_dir=TRACK_DIR, docs_dir=DOCS_DIR, images_dir=IMAGES_DIR, slides_dir=SLIDES_DIR, pipelines_dir=PIPELINES_DIR,
    doc_files=DOC_FILES, chapter_titles=CHAPTER_TITLES,
    track_label="Track 10",
    source_prefix="tracks/10_production_pipelines_mlops/docs/", track_prefix="tracks/10_production_pipelines_mlops/",
    md_name="Production_Pipelines_Master_Book.md", html_name="Production_Pipelines_Master_Course.html", docx_name="Production_Pipelines_Master_Course.docx",
    md_header="""# 확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps (Scalable Production Pipelines, HPC/Cloud Orchestration & MLOps) 마스터 교재
**Track 10: From Container Virtualization (Docker/Singularity), Nextflow DSL2 Reactive Channels, Snakemake Rule-Based DAGs, HPC SLURM/Cloud Batch Orchestration, nf-core/scrnaseq Custom Engineering, Cryptographic Data Provenance (W3C PROV/RO-Crate), and Single-Cell MLOps Model Serving (ONNX/Triton) to Cloud Cost Economics**

---

## 📚 교재 전체 목차 (Table of Contents)

- [Chapter 00. 프로덕션 파이프라인 및 MLOps 로드맵: 재현성 위기 극복부터 엔터프라이즈 임상 분석까지](#part-00-00-프로덕션-파이프라인-및-mlops-로드맵-재현성-위기-극복부터-엔터프라이즈-임상-분석까지-production-pipelines--mlops-roadmap)
- [Chapter 01. 컨테이너 가상화: Docker vs. Singularity/Apptainer HPC 아키텍처](#part-01-01-컨테이너-가상화-docker-vs-singularityapptainer-hpc-아키텍처-containerization-for-production-bioinformatics)
- [Chapter 02. Nextflow DSL2 모듈형 아키텍처와 리액티브 데이터플로우](#part-02-02-nextflow-dsl2-모듈형-아키텍처와-리액티브-데이터플로우-nextflow-dsl2-reactive-architecture)
- [Chapter 03. Snakemake 규칙 기반 실행 및 타깃 중심 DAG 아키텍처](#part-03-03-snakemake-규칙-기반-실행-및-타깃-중심-dag-아키텍처-snakemake-rule-based-workflow-architecture)
- [Chapter 04. HPC SLURM 클러스터 및 퍼블릭 클라우드 분산 오케스트레이션](#part-04-04-hpc-slurm-클러스터-및-퍼블릭-클라우드-분산-오케스트레이션-hpc-slurm--cloud-batch-orchestration)
- [Chapter 05. nf-core/scrnaseq 엔드투엔드 파이프라인 분석 및 커스텀 모듈 엔지니어링](#part-05-05-nf-corescrnaseq-엔드투엔드-파이프라인-분석-및-커스텀-모듈-엔지니어링-end-to-end-nf-corescrnaseq-pipeline--custom-module-engineering)
- [Chapter 06. 데이터 출처 추적성(Provenance), FAIR 원칙 및 암호학적 재현성 보증](#part-06-06-데이터-출처-추적성provenance-fair-원칙-및-암호학적-재현성-보증-data-provenance-fair-principles--cryptographic-reproducibility)
- [Chapter 07. 단일세포 전사체 MLOps 모델 서빙 및 실시간 추론 아키텍처](#part-07-07-단일세포-전사체-mlops-모델-서빙-및-실시간-추론-아키텍처-single-cell-mlops-model-serving--real-time-inference-architecture)
- [Chapter 08. 클라우드 비용 경제학, 스팟 인스턴스 내결함성 및 대규모 자원 최적화](#part-08-08-클라우드-비용-경제학-스팟-인스턴스-내결함성-및-대규모-자원-최적화-cloud-cost-economics-spot-fault-tolerance--resource-optimization)
- [Chapter 09. 프로덕션 파이프라인 트러블슈팅, 엔터프라이즈 모범 사례 및 캡스톤 프로젝트 루브릭](#part-09-09-프로덕션-파이프라인-트러블슈팅-엔터프라이즈-모범-사례-및-캡스톤-프로젝트-루브릭-production-troubleshooting-best-practices--capstone-rubric)

---

""",
    figcaption_desc="확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps 핵심 아키텍처",
    html_title=f"""확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps 마스터 교재 (Track 10: Scalable Production Pipelines, HPC/Cloud Orchestration &amp; MLOps)""",
    html_sidebar_title=f"""Track 10 프로덕션 &amp; MLOps""",
    html_search_placeholder=f"""전체 교재 검색 (예: Docker, Nextflow, DSL2, Snakemake, SLURM, nf-core, RO-Crate, Triton)...""",
    html_theme_key='production_mlops_doc_theme',
    html_hero_html=f"""<div class="course-hero">
      <h1>확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps 마스터 코스북</h1>
      <p class="subtitle">From Container Virtualization (Docker/Singularity), Nextflow DSL2 Reactive Channels, Snakemake Rule-Based DAGs, HPC SLURM/Cloud Batch Orchestration, nf-core/scrnaseq Custom Engineering, Cryptographic Data Provenance (W3C PROV/RO-Crate), and Single-Cell MLOps Model Serving (ONNX/Triton) to Cloud Cost Economics</p>
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
    docx_title="확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps 마스터 교재",
    docx_subtitle="Track 10: Scalable Production Pipelines, HPC/Cloud Orchestration & MLOps Master Course (Nextflow DSL2, Snakemake, SLURM, AWS Batch, nf-core, RO-Crate, ONNX/Triton)",
    docx_meta1="임상 & 유전체 모델: 단일세포 전사체 시퀀싱(snRNA-seq/scRNA-seq) 엔터프라이즈 프로덕션 파이프라인 구축 및 전립선암(CRPC/t-NEPC) 코호트의 대규모 분산 오케스트레이션, 암호학적 출처 추적성 및 MLOps 실시간 추론 마이크로서비스 서빙",
    docx_meta2="주요 내용: 재현성 위기 극복 및 '내 컴퓨터 증후군' 종식, Docker 컨테이너 및 Singularity/Apptainer 비루트 HPC 가상화, Nextflow DSL2 모듈형 리액티브 데이터플로우(Queue/Value 채널, 비동기 연산자, 동적 OOM 자가 치유 수리 모델), Snakemake 타깃 중심 규칙 기반 실행 및 와일드카드 DAG 역방향 위상정렬 알고리즘, HPC SLURM 배치 스케줄러(cgroups 자원 격리, 잡 어레이, sbatch 디렉티브) 및 퍼블릭 클라우드(AWS Batch, Nextflow Tower, GCP Life Sciences) 분산 오케스트레이션, nf-core/scrnaseq 골드 스탠더드 파이프라인 분해 및 커스텀 서브워크플로우 엔지니어링, W3C PROV 온톨로지·RO-Crate 메타데이터 패키징 및 SHA-256 암호학적 해시 검증 기반 FDA 21 CFR Part 11 출처 추적성(Provenance), 단일세포 딥러닝 모델의 ONNX 그래프 최적화·INT8 양자화·Triton Inference Server 동적 배칭 기반 초저지연 실시간 추론 마이크로서비스 서빙(SLA < 50ms), 클라우드 스팟 인스턴스 내결함성 체크포인팅 및 비용 경제학, 프로덕션 트러블슈팅과 캡스톤 프로젝트 루브릭",
    diagram_ref_label="프로덕션 파이프라인 시각 다이어그램 참조",
    docx_core_title="확장 가능한 프로덕션 파이프라인, HPC/클라우드 오케스트레이션 및 MLOps 마스터 교재",
    docx_core_subject="Track 10: Scalable Production Pipelines, HPC/Cloud Orchestration &amp; MLOps Master Course",
    svg_files=[
        "01_docker_vs_singularity_hpc_container_architecture.svg",
        "02_nextflow_dsl2_reactive_dataflow_channels.svg",
        "03_snakemake_target_driven_dag_resolution.svg",
        "04_hpc_slurm_and_cloud_batch_executor.svg",
        "05_data_provenance_fair_and_cryptographic_audit.svg",
        "06_single_cell_mlops_inference_lifecycle.svg"
    ],
    svg_exact=False,
    roadmap_terms=[
        "Track 10",
        "Production_Pipelines_Master_Book.md",
        "Production_Pipelines_Master_Course.html",
        "Production_Pipelines_Master_Course.docx",
        "tracks/10_production_pipelines_mlops/docs/",
        "tracks/10_production_pipelines_mlops/slides/",
        "tracks/10_production_pipelines_mlops/pipelines/"
    ],
    readme_terms=[
        "Track 10",
        "Production_Pipelines_Master_Book.md",
        "Production_Pipelines_Master_Course.html",
        "Production_Pipelines_Master_Course.docx",
        "tracks/10_production_pipelines_mlops/"
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
    assert "Nextflow DSL2" in p1_text or "Queue/Value channels" in p1_text, f"Pipeline 1 report did not record valid metrics!"
    assert "MLOps" in p2_text or "FDA 21 CFR Part 11" in p2_text, f"Pipeline 2 report did not record valid metrics!"


def _assert_pipelines():
    course_markdown.assert4_summary_report(CFG, "01_dsl2_dataflow_and_dag_engine.py", "02_mlops_model_serving_and_audit_pipeline.py", "summary_report.md", "summary_report.md",
                                           _check_pipeline_reports, ("  Re-running Pipeline 01 (DSL2 Dataflow & DAG Engine)...", "  Re-running Pipeline 02 (Single-Cell MLOps Serving & Audit)..."))


def verify_build_artifacts(html_content: str, docx_xml: str):
    course_markdown.verify_track(CFG, html_content, docx_xml, _assert_pipelines)


def main():
    print("🚀 Starting compilation of Track 10 Scalable Production Pipelines & MLOps Master Documents...")
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

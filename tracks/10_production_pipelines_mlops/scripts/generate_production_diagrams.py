#!/usr/bin/env python3
"""
tracks/10_production_pipelines_mlops/scripts/generate_production_diagrams.py
=============================================================================
Generates 6 publication-grade, high-resolution vector SVG diagrams for Track 10:
Production Pipelines, Workflow Orchestration & Single-Cell MLOps.

Diagrams:
1. 01_docker_vs_singularity_hpc_container_architecture.svg
2. 02_nextflow_dsl2_reactive_dataflow_channels.svg
3. 03_snakemake_target_driven_dag_resolution.svg
4. 04_hpc_slurm_and_cloud_batch_executor.svg
5. 05_data_provenance_fair_and_cryptographic_audit.svg
6. 06_single_cell_mlops_inference_lifecycle.svg
"""

import os
import xml.etree.ElementTree as ET


def make_01_docker_vs_singularity_hpc_container_architecture(filepath: str):
    """
    1. Docker vs. Singularity/Apptainer HPC Container Architecture & Security
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrRose01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrCyan01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrEmerald01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">Docker vs. Singularity/Apptainer HPC 컨테이너 보안 및 실행 아키텍처</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">Root 데몬 취약점(Docker) vs HPC 멀티유저 환경 비특권 SIF 사용자 공간 격리(Singularity/Apptainer)</text>

  <!-- Panel 1: Docker Root Daemon & HPC Vulnerability -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12" font-weight="700">1. Docker Root 데몬 &amp; HPC 취약점</text>

    <!-- Schematic: Docker Architecture -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    
    <!-- User CLI -->
    <rect x="25" y="52" width="65" height="26" rx="4" fill="#1e293b" stroke="#94a3b8" stroke-width="1"/>
    <text x="57" y="65" text-anchor="middle" fill="#f1f5f9" font-size="8" font-weight="700">User CLI</text>
    <text x="57" y="74" text-anchor="middle" fill="#94a3b8" font-size="7">UID: 1002</text>

    <!-- Arrow to dockerd -->
    <line x1="90" y1="65" x2="116" y2="65" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrRose01)"/>
    <text x="103" y="60" text-anchor="middle" fill="#f43f5e" font-size="6.5">socket</text>

    <!-- dockerd (Root Daemon) -->
    <rect x="120" y="48" width="68" height="34" rx="4" fill="#881337" stroke="#f43f5e" stroke-width="1.2"/>
    <text x="154" y="62" text-anchor="middle" fill="#ffe4e6" font-size="7.5" font-weight="700">dockerd</text>
    <text x="154" y="74" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">Root (UID 0)</text>

    <!-- Arrow to container -->
    <line x1="188" y1="65" x2="208" y2="65" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrRose01)"/>

    <!-- Container -->
    <rect x="212" y="52" width="63" height="26" rx="4" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
    <text x="243" y="64" text-anchor="middle" fill="#fb7185" font-size="7.5" font-weight="700">Container</text>
    <text x="243" y="73" text-anchor="middle" fill="#f43f5e" font-size="6.5">Root 권한상승</text>

    <!-- Threat warning banner -->
    <rect x="25" y="90" width="250" height="38" rx="4" fill="#4c0519" fill-opacity="0.6" stroke="#e11d48" stroke-width="1"/>
    <text x="150" y="103" text-anchor="middle" fill="#fecdd3" font-size="7.5" font-weight="700">🚨 보안 결함: 컨테이너 탈출(Breakout) 시</text>
    <text x="150" y="114" text-anchor="middle" fill="#fda4af" font-size="7">호스트 커널 및 공유 스토리지(/shared) 전체 장악 위험</text>
    <text x="150" y="123" text-anchor="middle" fill="#fda4af" font-size="6.5">다중 사용자(Multi-tenant) 클러스터에서 사용 불가</text>

    <!-- Bullet description -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 백그라운드 Root 데몬 의존성:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">클라이언트와 분리된 데몬이 호스트 슈퍼유저 권한 행사</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• POSIX 파일 권한 왜곡:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">생성된 파일이 root:root 소유로 생성되어 권한 충돌</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• HPC 보안 정책 위반:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">기관 공유 클러스터에서 일반 연구원 sudo 권한 부여 불가</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#881337" fill-opacity="0.4"/>
    <text x="150" y="259" text-anchor="middle" fill="#fb7185" font-size="9" font-weight="600">HPC 환경 배제: 슈퍼유저 탈취 및 공유 파일 파괴 위험</text>
  </g>

  <!-- Panel 2: Singularity/Apptainer Unprivileged SIF Architecture -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">2. Singularity/Apptainer 비특권 SIF 구조</text>

    <!-- Schematic: Singularity Architecture -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- User CLI -->
    <rect x="25" y="52" width="70" height="28" rx="4" fill="#0369a1" stroke="#38bdf8" stroke-width="1"/>
    <text x="60" y="65" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">User Process</text>
    <text x="60" y="75" text-anchor="middle" fill="#bae6fd" font-size="7">UID: 1002 (Host)</text>

    <!-- Direct invocation arrow -->
    <line x1="95" y1="66" x2="120" y2="66" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan01)"/>

    <!-- SIF Image -->
    <rect x="124" y="49" width="74" height="34" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="161" y="62" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">image.sif</text>
    <text x="161" y="74" text-anchor="middle" fill="#7dd3fc" font-size="6.5">SquashFS (RO)</text>

    <!-- Exec arrow -->
    <line x1="198" y1="66" x2="220" y2="66" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan01)"/>

    <!-- Executing Container -->
    <rect x="224" y="52" width="62" height="28" rx="4" fill="#0c4a6e" stroke="#38bdf8" stroke-width="1"/>
    <text x="255" y="65" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">Apptainer</text>
    <text x="255" y="75" text-anchor="middle" fill="#38bdf8" font-size="7">UID: 1002 (유지)</text>

    <!-- Mounts & Security banner -->
    <rect x="25" y="90" width="260" height="38" rx="4" fill="#082f49" fill-opacity="0.6" stroke="#0284c7" stroke-width="1"/>
    <text x="155" y="103" text-anchor="middle" fill="#bae6fd" font-size="7.5" font-weight="700">🛡️ 비특권 사용자 공간 마운트 (Unprivileged)</text>
    <text x="155" y="114" text-anchor="middle" fill="#e0f2fe" font-size="7">데몬 없음 (Daemonless) ➔ 프로세스 포크 방식으로 직접 실행</text>
    <text x="155" y="123" text-anchor="middle" fill="#7dd3fc" font-size="6.5">호스트 사용자 ID/그룹 ID 및 파일시스템 권한 100% 보존</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 단일 파일 이식성 (SIF 포맷):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">압축된 읽기 전용 SquashFS 파일 하나로 배포 및 공유</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 네이티브 HPC 스토리지 바인드:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">--bind /scratch,/shared 옵션으로 Lustre/GPFS 직접 접근</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• GPU 드라이버 네이티브 패스스루:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">--nv 플래그를 통한 NVIDIA CUDA 라이브러리 직접 바인딩</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">HPC 표준: 일반 유저 권한 실행 및 완전한 격리 보장</text>
  </g>

  <!-- Panel 3: OCI to SIF Conversion & SLURM Workflow -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. OCI 변환 및 SLURM 파이프라인 연계</text>

    <!-- Workflow Pipeline Graphic -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Step 1: Docker registry -->
    <rect x="22" y="50" width="72" height="32" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="58" y="62" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">DockerHub/Quay</text>
    <text x="58" y="74" text-anchor="middle" fill="#94a3b8" font-size="6.5">docker://starsolo</text>

    <line x1="94" y1="66" x2="110" y2="66" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald01)"/>

    <!-- Step 2: SIF cache -->
    <rect x="114" y="50" width="76" height="32" rx="3" fill="#064e3b" stroke="#34d399" stroke-width="1"/>
    <text x="152" y="62" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">NXF_SINGULARITY</text>
    <text x="152" y="74" text-anchor="middle" fill="#6ee7b7" font-size="6.5">CACHEDIR (*.sif)</text>

    <line x1="190" y1="66" x2="206" y2="66" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald01)"/>

    <!-- Step 3: SLURM Task -->
    <rect x="210" y="50" width="75" height="32" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="247" y="62" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">SLURM Worker</text>
    <text x="247" y="74" text-anchor="middle" fill="#fde68a" font-size="6.5">cgroups 64GB</text>

    <!-- Execution code snippet -->
    <rect x="22" y="90" width="263" height="38" rx="4" fill="#022c22" stroke="#059669" stroke-width="1"/>
    <text x="30" y="103" fill="#a7f3d0" font-family="monospace" font-size="7">$ singularity exec --bind /scratch \\</text>
    <text x="38" y="114" fill="#a7f3d0" font-family="monospace" font-size="7">  starsolo_2.7.sif STAR --soloType Droplet ...</text>
    <text x="30" y="123" fill="#6ee7b7" font-size="6.5">➔ SLURM 배치 스크립트 내 비대화형(Headless) 자동 실행</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• OCI 레이어 자동 병합 및 변환:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">Nextflow가 docker:// 태그를 감지하여 SIF로 자동 빌드</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 중앙 캐시를 통한 네트워크 I/O 최소화:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">수백 개 워커 노드가 공유 스토리지 SIF를 재참조</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 재현성 및 버전 고정(Digest Pinning):</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">sha256 다이제스트를 통한 소프트웨어 빌드 불변성 보장</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">SLURM + Singularity + Nextflow 무결점 통합</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_nextflow_dsl2_reactive_dataflow_channels(filepath: str):
    """
    2. Nextflow DSL2 Reactive Dataflow Channels & Process Concurrency
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad02" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrAmber02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrPurple02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">Nextflow DSL2 반응형 데이터플로우 채널 및 비동기 동시성 엔진</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">Dataflow 패러다임: Queue/Value 채널 분리, 선언적 고차 오퍼레이터 및 이벤트 기반 프로세스 스케줄링</text>

  <!-- Panel 1: Queue vs. Value Channels -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 채널 유형: Queue vs. Value 채널</text>

    <!-- Schematic: Queue & Value channels -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Queue Channel Box -->
    <rect x="22" y="48" width="256" height="38" rx="4" fill="#0369a1" fill-opacity="0.25" stroke="#38bdf8" stroke-width="1"/>
    <text x="30" y="60" fill="#38bdf8" font-size="8" font-weight="700">Queue Channel (FIFO, 일회성 소비):</text>
    
    <!-- Packets in Queue -->
    <rect x="30" y="65" width="68" height="16" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="64" y="76" text-anchor="middle" fill="#bae6fd" font-size="7">Sample_01</text>
    <rect x="103" y="65" width="68" height="16" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="137" y="76" text-anchor="middle" fill="#bae6fd" font-size="7">Sample_02</text>
    <rect x="176" y="65" width="68" height="16" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="210" y="76" text-anchor="middle" fill="#bae6fd" font-size="7">Sample_03</text>

    <!-- Value Channel Box -->
    <rect x="22" y="92" width="256" height="38" rx="4" fill="#581c87" fill-opacity="0.25" stroke="#c084fc" stroke-width="1"/>
    <text x="30" y="104" fill="#c084fc" font-size="8" font-weight="700">Value Channel (영구 브로드캐스트):</text>
    <rect x="30" y="109" width="120" height="16" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
    <text x="90" y="120" text-anchor="middle" fill="#e9d5ff" font-size="7">Reference Genome (GRCh38)</text>
    <text x="160" y="120" fill="#c084fc" font-size="6.5">➔ 모든 태스크에 무한 재사용</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• Queue Channel (Channel.fromPath):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">선입선출(FIFO) 비동기 버퍼, 요소 도착 즉시 다운스트림 트리거</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• Value Channel (Channel.value):</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">소진되지 않는 단일 상태 유지, 공통 참조 인덱스 공유</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 반응형 데이터 바인딩:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">폴링 없는 이벤트 기반 통신으로 스레드 락 및 블로킹 제거</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">생산자-소비자(Producer-Consumer) 비동기 큐잉</text>
  </g>

  <!-- Panel 2: DSL2 Reactive Operators -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">2. DSL2 선언적 채널 오퍼레이터</text>

    <!-- Operators Schema -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Operator block 1: map & filter -->
    <rect x="22" y="49" width="125" height="36" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="28" y="60" fill="#fbbf24" font-family="monospace" font-size="7.5" font-weight="700">.map { id, r1, r2 -&gt; }</text>
    <text x="28" y="72" fill="#94a3b8" font-family="monospace" font-size="7">  [meta(id), [r1, r2]]</text>
    <text x="28" y="80" fill="#fde68a" font-size="6.5">메타데이터 튜플 주입</text>

    <!-- Operator block 2: combine & groupTuple -->
    <rect x="155" y="49" width="132" height="36" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
    <text x="160" y="60" fill="#34d399" font-family="monospace" font-size="7.5" font-weight="700">.combine(genome_idx)</text>
    <text x="160" y="72" fill="#94a3b8" font-family="monospace" font-size="7">.groupTuple(by: 0)</text>
    <text x="160" y="80" fill="#6ee7b7" font-size="6.5">데카르트 곱 &amp; 샘플 병합</text>

    <!-- Operator block 3: branch & multi-ch -->
    <rect x="22" y="90" width="265" height="42" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="28" y="102" fill="#38bdf8" font-family="monospace" font-size="7.5" font-weight="700">reads.branch {</text>
    <text x="35" y="113" fill="#94a3b8" font-family="monospace" font-size="7">  pass: it.reads &gt; 50000 ➔ FASTQC | fail: ➔ DISCARD }</text>
    <text x="28" y="125" fill="#bae6fd" font-size="6.5">조건부 데이터플로우 분기(Branching) 제어</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 스트림 함수형 조작 (Functional Stream):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">map, filter, flatMap을 통한 비파괴적 데이터 스트림 가공</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 백프레셔(Backpressure) 자동 조절:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">소비 속도에 맞추어 업스트림 생성 속도를 제어하여 OOM 방지</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 모듈식 재사용 (DSL2 Modules):</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">프로세스와 워크플로우를 독립 모듈로 캡슐화 및 재사용</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#b45309" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">선언적 함수형 합성으로 복잡한 파이프라인 단순화</text>
  </g>

  <!-- Panel 3: Asynchronous Parallel Process Execution -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12" font-weight="700">3. 비동기 병렬 프로세스 실행 엔진</text>

    <!-- Execution graphic -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Process A -->
    <rect x="22" y="48" width="76" height="36" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="60" y="60" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">FASTQC</text>
    <text x="60" y="70" text-anchor="middle" fill="#94a3b8" font-size="6">work/a1/3b4c</text>
    <text x="60" y="78" text-anchor="middle" fill="#34d399" font-size="6.5">COMPLETED</text>

    <!-- Process B -->
    <rect x="108" y="48" width="82" height="36" rx="4" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="149" y="60" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">STARSOLO</text>
    <text x="149" y="70" text-anchor="middle" fill="#94a3b8" font-size="6">work/f2/89de</text>
    <text x="149" y="78" text-anchor="middle" fill="#fbbf24" font-size="6.5">RUNNING (48 cores)</text>

    <!-- Process C -->
    <rect x="200" y="48" width="85" height="36" rx="4" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
    <text x="242" y="60" text-anchor="middle" fill="#c084fc" font-size="7" font-weight="700">EMPTYDROPS</text>
    <text x="242" y="70" text-anchor="middle" fill="#94a3b8" font-size="6">work/7c/11aa</text>
    <text x="242" y="78" text-anchor="middle" fill="#94a3b8" font-size="6.5">WAITING (Ch input)</text>

    <!-- Dynamic Retry & workdir isolation -->
    <rect x="22" y="90" width="263" height="42" rx="3" fill="#3b0764" fill-opacity="0.4" stroke="#a855f7" stroke-width="1"/>
    <text x="30" y="103" fill="#e9d5ff" font-size="7.5" font-weight="700">🔄 동적 실패 복구 (Dynamic Error Strategy):</text>
    <text x="30" y="114" fill="#c084fc" font-family="monospace" font-size="7">errorStrategy 'retry', maxRetries 3, memory { 16.GB * task.attempt }</text>
    <text x="30" y="125" fill="#d8b4fe" font-size="6.5">OOM 137 에러 감지 시 메모리 2배 자동 증설 후 재시도</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 작업 디렉토리 완전 격리 (Workdir Sandbox):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">프로세스마다 독립적 고유 해시 디렉토리에서 I/O 실행</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 스마트 캐싱 (-resume):</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">입력 파일 해시 + 스크립트 코드 검증으로 기연산 작업 100% 스킵</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 무제한 스케일아웃 (Scale-out):</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">로컬 스레드 풀, SLURM 노드, AWS Batch 워커로 투명 확장</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#581c87" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#c084fc" font-size="9" font-weight="600">자가 치유(Self-Healing) 및 대규모 분산 스케일아웃</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_snakemake_target_driven_dag_resolution(filepath: str):
    """
    3. Snakemake Target-Driven DAG Resolution & Wildcards Pattern Expansion
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad03" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrAmber03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrEmerald03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">Snakemake 타깃 주도형 DAG 해결 및 와일드카드 전파 메커니즘</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">역방향 타깃 해결(Rule all ➔ input/output 매칭), 정규표현식 와일드카드 확장 및 위상 정렬 스케줄링</text>

  <!-- Panel 1: Target-Driven Bottom-Up Resolution -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 역방향 타깃 해결 (Bottom-Up)</text>

    <!-- Schematic: Rule hierarchy -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Rule All Target -->
    <rect x="22" y="48" width="256" height="20" rx="3" fill="#0369a1" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
    <text x="30" y="62" fill="#38bdf8" font-family="monospace" font-size="7.5" font-weight="700">rule all: input: "counts/{sample}.h5ad"</text>

    <!-- Backward arrow 1 -->
    <line x1="150" y1="68" x2="150" y2="76" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan03)"/>

    <!-- Rule Count -->
    <rect x="22" y="78" width="256" height="26" rx="3" fill="#1e293b" stroke="#334155" stroke-width="1"/>
    <text x="30" y="89" fill="#94a3b8" font-family="monospace" font-size="7">rule count: output: "counts/{sample}.h5ad"</text>
    <text x="30" y="99" fill="#fbbf24" font-family="monospace" font-size="6.5">  input: "aligned/{sample}.bam"</text>

    <!-- Backward arrow 2 -->
    <line x1="150" y1="104" x2="150" y2="112" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan03)"/>

    <!-- Rule Align -->
    <rect x="22" y="113" width="256" height="21" rx="3" fill="#1e293b" stroke="#334155" stroke-width="1"/>
    <text x="30" y="123" fill="#34d399" font-family="monospace" font-size="6.5">rule align: output: "aligned/{sample}.bam" &lt;- raw FASTQ</text>
    <text x="30" y="131" fill="#64748b" font-size="6">역방향 체이닝으로 필요한 규칙만 선택적 평가</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• GNU Make의 데이터 기반 진화:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">최종 타깃에서 시작하여 output ➔ input 매칭으로 의존성 추론</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 불필요한 연산의 원천 차단:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">타깃 산출물에 도달하지 않는 중간 산출물은 DAG에서 배제</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 파일 타임스탬프 기반 캐싱:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">입력 파일보다 최신인 출력 파일은 연산 없이 자동 재사용</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">타깃 주도형(Target-Driven) 역방향 의존성 해석</text>
  </g>

  <!-- Panel 2: Wildcards Pattern Matching & Expansion -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">2. 와일드카드 패턴 매칭 및 확장</text>

    <!-- Schematic: Wildcards expansion -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Pattern Box -->
    <rect x="22" y="48" width="265" height="34" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="30" y="60" fill="#fbbf24" font-family="monospace" font-size="7.5" font-weight="700">pattern: "qc/{sample}_{batch}_{run}.html"</text>
    <text x="30" y="73" fill="#94a3b8" font-size="7">정규식 토큰화: {sample: [A-Za-z0-9]+}, {batch: B[0-9]+}</text>

    <!-- Expand function -->
    <rect x="22" y="86" width="265" height="46" rx="3" fill="#78350f" fill-opacity="0.3" stroke="#d97706" stroke-width="1"/>
    <text x="30" y="99" fill="#fde68a" font-family="monospace" font-size="7.5" font-weight="700">expand("results/{pt}_{cond}.tsv",</text>
    <text x="40" y="110" fill="#fde68a" font-family="monospace" font-size="7">       pt=["P1", "P2"], cond=["Ctrl", "Drug"])</text>
    <text x="30" y="124" fill="#fbbf24" font-size="6.5">➔ 2 × 2 데카르트 곱(Cartesian Product) = 4개 잡 자동 생성</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 문자열 템플릿 기반 메타데이터 추출:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">파일 경로 자체로부터 샘플/배치/조건 변수를 자동으로 파싱</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 와일드카드 제약(Constraints):</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">wildcard_constraints를 설정하여 탐욕적 매칭 오류 완벽 방지</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 다차원 파라메트릭 스윕:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">수백 명 환자 코호트의 조합론적 작업망을 한 줄 코드로 확장</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#b45309" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">선언적 패턴 전파로 수천 개 병렬 태스크 자동 생성</text>
  </g>

  <!-- Panel 3: Topological Sort DAG & Job Scheduling -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. 위상 정렬 DAG 및 스케줄링</text>

    <!-- DAG Schema -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- DAG Nodes -->
    <circle cx="45" cy="65" r="14" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="45" y="68" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">QC (A)</text>

    <circle cx="45" cy="112" r="14" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="45" y="115" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">QC (B)</text>

    <!-- Edges -->
    <line x1="59" y1="65" x2="105" y2="85" stroke="#34d399" stroke-width="1.2" marker-end="url(#arrEmerald03)"/>
    <line x1="59" y1="112" x2="105" y2="92" stroke="#34d399" stroke-width="1.2" marker-end="url(#arrEmerald03)"/>

    <rect x="110" y="72" width="70" height="32" rx="4" fill="#064e3b" stroke="#34d399" stroke-width="1"/>
    <text x="145" y="85" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">STARsolo</text>
    <text x="145" y="96" text-anchor="middle" fill="#6ee7b7" font-size="6.5">In-degree: 0 도달</text>

    <line x1="180" y1="88" x2="215" y2="88" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald03)"/>

    <circle cx="235" cy="88" r="18" fill="#1e293b" stroke="#fbbf24" stroke-width="1.5"/>
    <text x="235" y="87" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">Report</text>
    <text x="235" y="97" text-anchor="middle" fill="#94a3b8" font-size="6">(Target)</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• Kahn 알고리즘 기반 위상 정렬:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">진입차수(In-degree)가 0인 작업부터 순차/병렬 큐에 투입</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 순환 참조 검출 (Cycle Detection):</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">A ➔ B ➔ A 순환 발생 시 실행 전 즉각 PeriodicRuleException 방출</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 드라이런(-n) 사전 시뮬레이션:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">snakemake -n -p 로 파일 입출력 및 전체 DAG 유효성 선검증</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">위상 정렬 기반의 수학적 무결성 및 최적 자원 패킹</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_hpc_slurm_and_cloud_batch_executor(filepath: str):
    """
    4. On-Premises SLURM Cluster vs. AWS/GCP Cloud Batch Auto-scaling
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad04" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrEmerald04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrPurple04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">HPC SLURM 클러스터 vs. 퍼블릭 클라우드 배치 실행 아키텍처</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">온프레미스 고정 노드 큐잉(SLURM) vs 클라우드 탄력적 스팟 인스턴스(Cloud Batch) &amp; 오브젝트 스토리지 패브릭</text>

  <!-- Panel 1: On-Premises SLURM Cluster Executor -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 온프레미스 HPC SLURM 클러스터</text>

    <!-- Schematic: SLURM Nodes -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Controller -->
    <rect x="22" y="48" width="70" height="26" rx="3" fill="#0369a1" stroke="#38bdf8" stroke-width="1"/>
    <text x="57" y="60" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">slurmctld</text>
    <text x="57" y="70" text-anchor="middle" fill="#bae6fd" font-size="6.5">마스터 스케줄러</text>

    <line x1="92" y1="61" x2="114" y2="61" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan04)"/>

    <!-- Worker Nodes -->
    <rect x="118" y="47" width="160" height="42" rx="3" fill="#1e293b" stroke="#334155" stroke-width="1"/>
    <text x="198" y="59" text-anchor="middle" fill="#f1f5f9" font-size="7.5" font-weight="700">Compute Nodes (Partition: 'scRNA')</text>
    <rect x="124" y="65" width="70" height="18" rx="2" fill="#0f172a" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="159" y="77" text-anchor="middle" fill="#7dd3fc" font-size="6.5">Node 01 (64c/512G)</text>
    <rect x="202" y="65" width="70" height="18" rx="2" fill="#0f172a" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="237" y="77" text-anchor="middle" fill="#7dd3fc" font-size="6.5">Node 02 (64c/512G)</text>

    <!-- POSIX Storage Bar -->
    <rect x="22" y="94" width="256" height="38" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="30" y="106" fill="#38bdf8" font-size="7.5" font-weight="700">Shared POSIX Storage (/scratch):</text>
    <text x="30" y="117" fill="#94a3b8" font-size="7">Lustre / IBM GPFS 병렬 파일시스템 (직접 인피니밴드 RDMA)</text>
    <text x="30" y="127" fill="#fde68a" font-size="6.5">한계: 메타데이터 락 경합 및 대규모 동시 I/O 병목</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 고정 인프라 자원 관리:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">정적 노드 풀, 공정 점유(Fair-share) 큐 정책으로 대기 시간 발생</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• cgroups 하드웨어 엄격 제어:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">CPU/RAM 초과 시 커널 OOM Killer로 작업 강제 종료 (Code 137)</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 지속적 운영비(CapEx) 효율성:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">24/7 지속 가동하는 대규모 유전체 센터에 적합</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">대규모 상시 워크로드에 최적화된 온프레미스 HPC</text>
  </g>

  <!-- Panel 2: Cloud Batch & Serverless Auto-scaling -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">2. 클라우드 배치 &amp; 스팟 자동확장</text>

    <!-- Schematic: Cloud Batch -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Orchestrator -->
    <rect x="22" y="48" width="75" height="28" rx="3" fill="#065f46" stroke="#34d399" stroke-width="1"/>
    <text x="59" y="60" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">Nextflow Tower</text>
    <text x="59" y="71" text-anchor="middle" fill="#a7f3d0" font-size="6.5">Cloud API Hub</text>

    <line x1="97" y1="62" x2="120" y2="62" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald04)"/>

    <!-- Dynamic Fleet -->
    <rect x="124" y="47" width="158" height="42" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
    <text x="203" y="58" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">AWS Batch / GCP Life Sciences</text>
    <text x="203" y="70" text-anchor="middle" fill="#6ee7b7" font-size="7">0 ➔ 100 ➔ 0 VM 동적 프로비저닝</text>
    <text x="203" y="82" text-anchor="middle" fill="#a7f3d0" font-size="6.5">Spot / Preemptible (70~80% 비용 절감)</text>

    <!-- Fault tolerant checkpoint -->
    <rect x="22" y="94" width="260" height="38" rx="3" fill="#022c22" stroke="#059669" stroke-width="1"/>
    <text x="30" y="106" fill="#6ee7b7" font-size="7.5" font-weight="700">⚡ 스팟 회수 내성 (Preemption Resilience):</text>
    <text x="30" y="117" fill="#a7f3d0" font-size="7">인스턴스 회수 2분 경고 수신 시 작업 상태 S3 체크포인트 저장</text>
    <text x="30" y="127" fill="#d1fae5" font-size="6.5">신규 스팟 VM에서 -resume으로 즉각 무손실 복구 실행</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 제로 대기시간 (Zero Queue Time):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">큐 적체 없이 분석 요청 즉시 가상 머신 수백 대를 탄력 할당</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 초단위 과금(OpEx) 효율화:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">분석 완료 즉시 VM을 종료하여 유휴 자원 비용 완전 제거</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 최신 GPU 인프라 즉시 도입:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">H100/A100 등 고성능 가속기를 별도 하드웨어 구매 없이 활용</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">탄력적 대규모 버스트 분석 및 클라우드 비용 극대화</text>
  </g>

  <!-- Panel 3: Storage Fabric: POSIX vs. Object Store -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12" font-weight="700">3. 하이브리드 스토리지 패브릭</text>

    <!-- Schematic: S3 Staging & Compute -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Object Store -->
    <rect x="22" y="48" width="80" height="34" rx="3" fill="#581c87" stroke="#c084fc" stroke-width="1"/>
    <text x="62" y="60" text-anchor="middle" fill="#f3e8ff" font-size="7.5" font-weight="700">Cloud Storage</text>
    <text x="62" y="70" text-anchor="middle" fill="#d8b4fe" font-size="6.5">AWS S3 / GCS</text>
    <text x="62" y="78" text-anchor="middle" fill="#e9d5ff" font-size="6">페타바이트 불변 저장</text>

    <!-- Staging arrows -->
    <line x1="102" y1="60" x2="132" y2="60" stroke="#c084fc" stroke-width="1.5" marker-end="url(#arrPurple04)"/>
    <text x="117" y="55" text-anchor="middle" fill="#c084fc" font-size="6">Stage-in</text>

    <line x1="132" y1="74" x2="102" y2="74" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald04)"/>
    <text x="117" y="83" text-anchor="middle" fill="#34d399" font-size="6">Stage-out</text>

    <!-- Ephemeral Worker -->
    <rect x="136" y="48" width="145" height="38" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="208" y="59" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">Ephemeral Worker Instance</text>
    <rect x="142" y="65" width="133" height="17" rx="2" fill="#0f172a" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="208" y="77" text-anchor="middle" fill="#7dd3fc" font-size="6.5">Local NVMe SSD (/tmp 고속 연산)</text>

    <!-- Storage decoupling explanation -->
    <rect x="22" y="94" width="260" height="38" rx="3" fill="#3b0764" fill-opacity="0.3" stroke="#a855f7" stroke-width="1"/>
    <text x="30" y="106" fill="#e9d5ff" font-size="7.5" font-weight="700">💾 컴퓨팅-스토리지 완전 분리 (Decoupled Fabric):</text>
    <text x="30" y="117" fill="#d8b4fe" font-size="7">영구 보관은 99.999999999% 내구성의 오브젝트 스토리지(S3)</text>
    <text x="30" y="127" fill="#c084fc" font-size="6.5">연산 시에만 로컬 고속 NVMe에 캐싱 후 최종 결과물 역전송</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• HTTP/S3 기반 스트리밍 I/O:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">POSIX 파일 핸들 부족 및 네트워크 파일 락 문제 원천 해소</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 데이터 지역성(Locality) 최적화:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">리전(Region) 내 데이터 전송 무료 정책 활용으로 Egress 비용 0원</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 멀티 클라우드 및 하이브리드 연동:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">온프레미스 원천 데이터 ➔ 클라우드 연산 ➔ 결과 저장의 유연성</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#581c87" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#c084fc" font-size="9" font-weight="600">페타바이트급 오믹스 빅데이터의 클라우드 스토리지 패브릭</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_data_provenance_fair_and_cryptographic_audit(filepath: str):
    """
    5. Data Provenance, FAIR Principles & Cryptographic Audit Trail
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad05" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrAmber05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrEmerald05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">데이터 출처성(Provenance), FAIR 원칙 및 암호학적 감사 추적체계</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">Git 커밋 SHA, 불변 컨테이너 다이제스트, SHA-256 체크섬을 통한 규제 준수(FDA 21 CFR Part 11)</text>

  <!-- Panel 1: Code & Runtime Immutability Triad -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 코드 및 런타임 불변성 삼합 (Triad)</text>

    <!-- Schematic: Triad -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Pillar 1: Git SHA -->
    <rect x="22" y="48" width="76" height="34" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="60" y="59" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">Git Commit</text>
    <text x="60" y="70" text-anchor="middle" fill="#bae6fd" font-family="monospace" font-size="6.5">sha: e84f29a</text>
    <text x="60" y="78" text-anchor="middle" fill="#94a3b8" font-size="6">파이프라인 코드</text>

    <!-- Pillar 2: Container Digest -->
    <rect x="108" y="48" width="82" height="34" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="149" y="59" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">Image Digest</text>
    <text x="149" y="70" text-anchor="middle" fill="#fde68a" font-family="monospace" font-size="6.5">sha256: 7f8a...</text>
    <text x="149" y="78" text-anchor="middle" fill="#94a3b8" font-size="6">소프트웨어 환경</text>

    <!-- Pillar 3: Lockfile -->
    <rect x="198" y="48" width="80" height="34" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
    <text x="238" y="59" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Params Lock</text>
    <text x="238" y="70" text-anchor="middle" fill="#a7f3d0" font-family="monospace" font-size="6.5">params.json</text>
    <text x="238" y="78" text-anchor="middle" fill="#94a3b8" font-size="6">하이퍼파라미터</text>

    <!-- Bottom summary banner -->
    <rect x="22" y="90" width="256" height="38" rx="3" fill="#082f49" fill-opacity="0.4" stroke="#0284c7" stroke-width="1"/>
    <text x="150" y="103" text-anchor="middle" fill="#bae6fd" font-size="7.5" font-weight="700">🔒 절대적 재현성 보장 공식</text>
    <text x="150" y="114" text-anchor="middle" fill="#e0f2fe" font-size="7">Output = F(Code_SHA, Environment_Digest, Input_Hash)</text>
    <text x="150" y="124" text-anchor="middle" fill="#7dd3fc" font-size="6.5">태그 부동(:latest)을 엄격히 금지하고 암호학적 해시로 고정</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 부동 태그(Floating Tag) 차단 정책:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">latest 태그 사용 시 시간 경과에 따른 침묵적 변형(Drift) 발생</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 파이프라인 매니페스트 고정:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">Nextflow/Snakemake 구성 파일 내 실행 버전 명시적 잠금</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 10년 후에도 동일한 바이트 재현:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">동일 아티팩트 해시를 통해 논문 및 임상 규제 승인 영구 보증</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">완벽한 과학적 재현성을 위한 불변 3요소</text>
  </g>

  <!-- Panel 2: FAIR Data Provenance Graph (W3C PROV-O) -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">2. FAIR 데이터 출처 그래프 (W3C PROV-O)</text>

    <!-- Schematic: W3C PROV -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Entity 1 -->
    <ellipse cx="55" cy="65" rx="35" ry="16" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="55" y="65" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">Entity: FASTQ</text>
    <text x="55" y="73" text-anchor="middle" fill="#94a3b8" font-size="6">SHA256: 4a2f...</text>

    <line x1="90" y1="65" x2="118" y2="65" stroke="#fbbf24" stroke-width="1.2" marker-end="url(#arrAmber05)"/>
    <text x="104" y="58" text-anchor="middle" fill="#fbbf24" font-size="6">used</text>

    <!-- Activity -->
    <rect x="122" y="48" width="70" height="34" rx="3" fill="#78350f" fill-opacity="0.4" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="157" y="60" text-anchor="middle" fill="#fde68a" font-size="7" font-weight="700">Activity</text>
    <text x="157" y="70" text-anchor="middle" fill="#fbbf24" font-size="6.5">STARsolo_Align</text>
    <text x="157" y="78" text-anchor="middle" fill="#94a3b8" font-size="5.5">Job_ID: #4092</text>

    <line x1="192" y1="65" x2="220" y2="65" stroke="#fbbf24" stroke-width="1.2" marker-end="url(#arrAmber05)"/>
    <text x="206" y="58" text-anchor="middle" fill="#fbbf24" font-size="6">generated</text>

    <!-- Entity 2 -->
    <ellipse cx="255" cy="65" rx="32" ry="16" fill="#1e293b" stroke="#34d399" stroke-width="1.2"/>
    <text x="255" y="65" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Entity: h5ad</text>
    <text x="255" y="73" text-anchor="middle" fill="#94a3b8" font-size="6">SHA256: 9b1e...</text>

    <!-- Agent at bottom -->
    <rect x="110" y="94" width="94" height="24" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
    <text x="157" y="105" text-anchor="middle" fill="#c084fc" font-size="6.5" font-weight="700">Agent: Nextflow / nf-core</text>
    <text x="157" y="114" text-anchor="middle" fill="#e9d5ff" font-size="6">wasAssociatedWith (Activity)</text>
    <line x1="157" y1="94" x2="157" y2="82" stroke="#c084fc" stroke-width="1" stroke-dasharray="2,2"/>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• W3C PROV-O 데이터 모델 구현:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">Entity(데이터), Activity(프로세스), Agent(분석자/도구) 삼각 관계</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• RO-Crate 메타데이터 패키징:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">JSON-LD 기반 파일 간 계보 및 전사적 속성을 기계 판독형으로 기록</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• FAIR 원칙의 완벽한 충족:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">검색 가능(F), 접근 가능(A), 상호운용(I), 재사용 가능(R) 보증</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#b45309" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">글로벌 표준 출처 그래프 모델로 투명한 데이터 계보 구축</text>
  </g>

  <!-- Panel 3: Cryptographic Audit Trail & FDA 21 CFR Part 11 -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. 암호학적 감사 추적 (Part 11)</text>

    <!-- Schematic: Audit ledger -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Audit record 1 -->
    <rect x="22" y="48" width="265" height="24" rx="3" fill="#064e3b" fill-opacity="0.3" stroke="#34d399" stroke-width="0.8"/>
    <text x="30" y="58" fill="#34d399" font-family="monospace" font-size="6.5">TS: 2026-10-16T14:22:01Z | User: analyst_dr</text>
    <text x="30" y="67" fill="#6ee7b7" font-family="monospace" font-size="6">OP: MODEL_FIT | Checksum: sha256:d5a491e0... [VERIFIED]</text>

    <!-- Audit record 2 -->
    <rect x="22" y="75" width="265" height="24" rx="3" fill="#064e3b" fill-opacity="0.3" stroke="#34d399" stroke-width="0.8"/>
    <text x="30" y="85" fill="#34d399" font-family="monospace" font-size="6.5">TS: 2026-10-16T14:35:12Z | Action: PROMOTE_PROD</text>
    <text x="30" y="94" fill="#6ee7b7" font-family="monospace" font-size="6">SIG: RSA-4096 Signed by Chief Pathologist [OK]</text>

    <!-- Tamper evident callout -->
    <rect x="22" y="103" width="265" height="28" rx="3" fill="#022c22" stroke="#059669" stroke-width="1"/>
    <text x="30" y="115" fill="#a7f3d0" font-size="7" font-weight="700">🛡️ 위변조 방지 원장 (Tamper-Evident Ledger):</text>
    <text x="30" y="125" fill="#6ee7b7" font-size="6.5">블록체인/해시체인 연결로 사후 수정이나 조작 시 즉시 무효화</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• FDA 21 CFR Part 11 규제 준수:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">전자 서명(Electronic Signature), 사용자 식별 및 불변 감사 기록</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 자동화된 런타임 추적 테이블:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">CPU 사용률, peak RSS 메모리, I/O 읽기/쓰기 바이트 자동 집계</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• CLIA/CAP 임상 검사실 공인 요건:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">단일세포 암 진단 파이프라인의 결과 신뢰도 및 재현성 공식 보증</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">임상 진단 및 신약 개발을 위한 완벽한 규제 감사 추적</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_single_cell_mlops_inference_lifecycle(filepath: str):
    """
    6. Single-Cell MLOps Full Lifecycle: Training to Low-Latency REST Serving
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad06" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrAmber06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrRose06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">단일세포 전주기 MLOps 및 실시간 임상 추론 서빙 수명주기</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">FASTQ 전처리 ➔ scVI 딥러닝 잠재공간 학습 ➔ MLflow 모델 레지스트리 ➔ FastAPI 저지연 임상 진단</text>

  <!-- Panel 1: Upstream Pipeline to Latent Space -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 업스트림 정제 및 잠재공간 학습</text>

    <!-- Schematic: Upstream to scVI -->
    <rect x="15" y="42" width="270" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Fastq & nf-core -->
    <rect x="22" y="48" width="68" height="34" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="56" y="59" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">Raw FASTQ</text>
    <text x="56" y="69" text-anchor="middle" fill="#bae6fd" font-size="6.5">nf-core/scrnaseq</text>
    <text x="56" y="77" text-anchor="middle" fill="#94a3b8" font-size="6">10x Genomics</text>

    <line x1="90" y1="65" x2="108" y2="65" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan06)"/>

    <!-- Counts Matrix -->
    <rect x="112" y="48" width="68" height="34" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
    <text x="146" y="59" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Counts h5ad</text>
    <text x="146" y="69" text-anchor="middle" fill="#a7f3d0" font-size="6.5">20,000 유전자</text>
    <text x="146" y="77" text-anchor="middle" fill="#94a3b8" font-size="6">QC &amp; Ambient 제거</text>

    <line x1="180" y1="65" x2="198" y2="65" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan06)"/>

    <!-- scVI Model -->
    <rect x="202" y="48" width="75" height="34" rx="3" fill="#0c4a6e" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="239" y="59" text-anchor="middle" fill="#ffffff" font-size="7" font-weight="700">scVI / VAE</text>
    <text x="239" y="69" text-anchor="middle" fill="#7dd3fc" font-size="6.5">Z-Latent (10D)</text>
    <text x="239" y="77" text-anchor="middle" fill="#bae6fd" font-size="6">배치 효과 완벽 제거</text>

    <!-- Training Callout -->
    <rect x="22" y="90" width="256" height="38" rx="3" fill="#082f49" fill-opacity="0.4" stroke="#0284c7" stroke-width="1"/>
    <text x="150" y="103" text-anchor="middle" fill="#bae6fd" font-size="7.5" font-weight="700">🧬 심층 생성 모델(Deep Generative Model) 학습</text>
    <text x="150" y="114" text-anchor="middle" fill="#e0f2fe" font-size="7">Zero-Inflated Negative Binomial (ZINB) 분포 기반</text>
    <text x="150" y="124" text-anchor="middle" fill="#7dd3fc" font-size="6.5">기술적 잡음 분리 및 생물학적 세포 내재 상태 추출</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 자동화된 파이프라인 정제:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">FASTQ에서 AnnData 단일세포 객체까지 무인 자동화 전처리</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 전립선암 바이오마커 피처링:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">AR, SYP, CHGA, TMPRSS2-ERG 융합 스플라이싱 피처 통합</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• GPU 기반 고속 잠재공간 학습:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">100만 개 단일세포 전사체를 PyTorch 분산 환경에서 초고속 임베딩</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">원천 서열에서 고차원 잠재 표현형으로의 자동 추출</text>
  </g>

  <!-- Panel 2: Model Registry & Artifact Versioning -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">2. 모델 레지스트리 &amp; 패키징 (MLflow)</text>

    <!-- Schematic: MLflow Registry -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- Training Experiments -->
    <rect x="22" y="48" width="80" height="36" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="62" y="59" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">MLflow Tracking</text>
    <text x="62" y="70" text-anchor="middle" fill="#fde68a" font-size="6.5">Loss, F1, Val Acc</text>
    <text x="62" y="79" text-anchor="middle" fill="#94a3b8" font-size="6">Run: #9832a</text>

    <line x1="102" y1="66" x2="124" y2="66" stroke="#fbbf24" stroke-width="1.5" marker-end="url(#arrAmber06)"/>

    <!-- Model Registry Staging -> Prod -->
    <rect x="128" y="48" width="158" height="36" rx="3" fill="#78350f" fill-opacity="0.3" stroke="#fbbf24" stroke-width="1"/>
    <text x="207" y="59" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">MLflow Model Registry</text>
    <rect x="134" y="65" width="70" height="15" rx="2" fill="#1e293b" stroke="#94a3b8" stroke-width="0.8"/>
    <text x="169" y="75" text-anchor="middle" fill="#94a3b8" font-size="6">Stage: Staging</text>
    <line x1="204" y1="72" x2="216" y2="72" stroke="#34d399" stroke-width="1" marker-end="url(#arrEmerald04)"/>
    <rect x="218" y="65" width="62" height="15" rx="2" fill="#065f46" stroke="#34d399" stroke-width="0.8"/>
    <text x="249" y="75" text-anchor="middle" fill="#a7f3d0" font-size="6">Production</text>

    <!-- Artifacts & ONNX -->
    <rect x="22" y="90" width="265" height="38" rx="3" fill="#1e293b" stroke="#334155" stroke-width="1"/>
    <text x="30" y="103" fill="#fbbf24" font-size="7.5" font-weight="700">📦 모델 아티팩트 서명 및 ONNX 런타임 변환:</text>
    <text x="30" y="114" fill="#fde68a" font-family="monospace" font-size="6.5">weights.onnx (SHA-256: 8f2b...) + schema.json</text>
    <text x="30" y="124" fill="#94a3b8" font-size="6.5">Python 독립적 C++ 런타임 최적화로 추론 오버헤드 80% 감축</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 중앙화된 모델 수명주기 관리:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">실험 로깅, 파라미터 계측, 모델 버전 승격을 단일 UI에서 추적</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 암호학적 모델 카드(Model Card):</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">학습 데이터 출처, 가중치 체크섬, 한계점을 메타데이터로 고정</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 온톨로지 기반 입력 스키마 검증:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">유전자 심볼 일치성 및 카운트 정규화 조건 자동 검증</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#b45309" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">실험 관리에서 검증된 모델 아티팩트의 엔터프라이즈 승격</text>
  </g>

  <!-- Panel 3: Low-Latency Serving & Clinical Decision Support -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#f43f5e" font-size="12" font-weight="700">3. 고속 서빙 &amp; 임상 진단 대시보드</text>

    <!-- Schematic: FastAPI & Inference -->
    <rect x="15" y="42" width="280" height="96" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>

    <!-- FastAPI Container -->
    <rect x="22" y="48" width="80" height="34" rx="3" fill="#881337" fill-opacity="0.5" stroke="#f43f5e" stroke-width="1"/>
    <text x="62" y="59" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">FastAPI / Triton</text>
    <text x="62" y="70" text-anchor="middle" fill="#fca5a5" font-size="6.5">비동기 마이크로서비스</text>
    <text x="62" y="78" text-anchor="middle" fill="#fecdd3" font-size="6">&lt; 50ms 초저지연</text>

    <line x1="102" y1="65" x2="124" y2="65" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrRose06)"/>

    <!-- Diagnosis Result -->
    <rect x="128" y="48" width="158" height="34" rx="3" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
    <text x="207" y="59" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">임상 진단 추론 (Single-Cell)</text>
    <text x="207" y="69" text-anchor="middle" fill="#fecdd3" font-size="6.5">t-NEPC 신경내분비 분화 (98.2%)</text>
    <text x="207" y="77" text-anchor="middle" fill="#fda4af" font-size="6">Enzalutamide 내성 위험군 경보</text>

    <!-- Monitoring & Drift -->
    <rect x="22" y="90" width="265" height="38" rx="3" fill="#4c0519" fill-opacity="0.4" stroke="#e11d48" stroke-width="1"/>
    <text x="30" y="103" fill="#fecdd3" font-size="7.5" font-weight="700">📊 실시간 드리프트 감지 (Evidently / Prometheus):</text>
    <text x="30" y="114" fill="#fda4af" font-size="7">KS-test 유전자 발현 드리프트 감지 시 재학습 트리거</text>
    <text x="30" y="124" fill="#fca5a5" font-size="6.5">의료진 대시보드(Streamlit/React)와 즉각적 소통</text>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 고성능 비동기 배치 추론 (Async Batching):</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">세포 수천 개의 발현 벡터를 50ms 미만으로 일괄 분류</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 전립선암 거세저항성(CRPC) 전이 예측:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">AR 양성 선암종 ➔ t-NEPC 전이 서브클론을 조기 감지</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 지속적 모델 관측성 (Continuous Observability):</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">추론 레이턴시, 메모리, 예측 확신도 모니터링</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#881337" fill-opacity="0.4"/>
    <text x="155" y="259" text-anchor="middle" fill="#f43f5e" font-size="9" font-weight="600">단일세포 분해능의 실시간 정밀 종양학 의사결정 지원</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_docker_vs_singularity_hpc_container_architecture.svg", make_01_docker_vs_singularity_hpc_container_architecture),
        ("02_nextflow_dsl2_reactive_dataflow_channels.svg", make_02_nextflow_dsl2_reactive_dataflow_channels),
        ("03_snakemake_target_driven_dag_resolution.svg", make_03_snakemake_target_driven_dag_resolution),
        ("04_hpc_slurm_and_cloud_batch_executor.svg", make_04_hpc_slurm_and_cloud_batch_executor),
        ("05_data_provenance_fair_and_cryptographic_audit.svg", make_05_data_provenance_fair_and_cryptographic_audit),
        ("06_single_cell_mlops_inference_lifecycle.svg", make_06_single_cell_mlops_inference_lifecycle),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 10 Production & MLOps publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

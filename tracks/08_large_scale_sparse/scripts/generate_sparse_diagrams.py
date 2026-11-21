#!/usr/bin/env python3
"""
tracks/08_large_scale_sparse/scripts/generate_sparse_diagrams.py
=============================================================================
Generates 6 publication-grade, high-resolution vector SVG diagrams for Track 08:
Large-Scale Single-Cell Atlases, Sparse Computation & GPU Acceleration.

Diagrams:
1. 01_sparse_matrix_csr_csc_memory_layout.svg
2. 02_hdf5_backed_mode_vs_zarr_chunking.svg
3. 03_tiledb_soma_sparse_tensor_architecture.svg
4. 04_dask_lazy_dag_out_of_core_normalization.svg
5. 05_rapids_gpu_acceleration_cuda_pipeline.svg
6. 06_cagra_gpu_knn_and_cugraph_leiden.svg
"""

import os
import xml.etree.ElementTree as ET


def make_01_sparse_matrix_csr_csc_memory_layout(filepath: str):
    """
    1. Sparse Matrix Memory Layout: COO vs CSR vs CSC & 64B CPU Cache Line Hit Rate
    Pointer arrays (indptr, indices, data) and row vs col slicing memory geometry.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrCyan01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrAmber01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrEmerald01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrRose01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">희소 행렬 메모리 레이아웃: COO vs CSR vs CSC 및 CPU 캐시 라인(64B) 적중률 기하학</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">포인터 배열(indptr, indices, data) ➔ 행(세포) vs 열(유전자) 슬라이싱 성능 격차 ➔ L1/L2 캐시 스트리밍</text>

  <!-- Panel 1: COO Format -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. COO 형식 (Coordinate Format)</text>

    <!-- Schematic: COO 3 Arrays -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#e11d48" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Row array -->
      <text x="0" y="14" fill="#94a3b8" font-size="7.5" font-weight="700">row :</text>
      <rect x="35" y="4" width="26" height="14" rx="2" fill="#1e293b" stroke="#e11d48" stroke-width="0.8"/>
      <text x="48" y="14" text-anchor="middle" fill="#f43f5e" font-size="7.5">0</text>
      <rect x="65" y="4" width="26" height="14" rx="2" fill="#1e293b" stroke="#e11d48" stroke-width="0.8"/>
      <text x="78" y="14" text-anchor="middle" fill="#f43f5e" font-size="7.5">0</text>
      <rect x="95" y="4" width="26" height="14" rx="2" fill="#1e293b" stroke="#e11d48" stroke-width="0.8"/>
      <text x="108" y="14" text-anchor="middle" fill="#f43f5e" font-size="7.5">1</text>
      <rect x="125" y="4" width="26" height="14" rx="2" fill="#1e293b" stroke="#e11d48" stroke-width="0.8"/>
      <text x="138" y="14" text-anchor="middle" fill="#f43f5e" font-size="7.5">2</text>
      <text x="165" y="14" fill="#64748b" font-size="7.5">... (NNZ개)</text>

      <!-- Col array -->
      <text x="0" y="34" fill="#94a3b8" font-size="7.5" font-weight="700">col :</text>
      <rect x="35" y="24" width="26" height="14" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="48" y="34" text-anchor="middle" fill="#38bdf8" font-size="7.5">1</text>
      <rect x="65" y="24" width="26" height="14" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="78" y="34" text-anchor="middle" fill="#38bdf8" font-size="7.5">3</text>
      <rect x="95" y="24" width="26" height="14" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="108" y="34" text-anchor="middle" fill="#38bdf8" font-size="7.5">2</text>
      <rect x="125" y="24" width="26" height="14" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="138" y="34" text-anchor="middle" fill="#38bdf8" font-size="7.5">0</text>
      <text x="165" y="34" fill="#64748b" font-size="7.5">... (NNZ개)</text>

      <!-- Data array -->
      <text x="0" y="54" fill="#94a3b8" font-size="7.5" font-weight="700">data:</text>
      <rect x="35" y="44" width="26" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="48" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">4.2</text>
      <rect x="65" y="44" width="26" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="78" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">1.8</text>
      <rect x="95" y="44" width="26" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="108" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">9.1</text>
      <rect x="125" y="44" width="26" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="138" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">3.5</text>
      <text x="165" y="54" fill="#64748b" font-size="7.5">... (NNZ개)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 단순 튜플 저장소 (Triplets):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">(row, col, value) 3개 벡터를 비압축 나열</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 메모리 오버헤드: 3 × NNZ × 8 B</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">행 번호가 반복 저장되어 심각한 메모리 낭비</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• CPU 캐시 라인 미스 (Cache Miss):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">연속 인덱스 부재로 64B 캐시 프리페치 적중률 &lt; 20%</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.2" font-weight="600">용도: 10x Genomics 로우 파일 적재 전용 (연산 부적합)</text>
  </g>

  <!-- Panel 2: CSR Format -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">2. CSR 형식 (행/세포 슬라이싱 최적화)</text>

    <!-- Schematic: CSR Structure -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- indptr array -->
      <text x="0" y="14" fill="#94a3b8" font-size="7.5" font-weight="700">indptr :</text>
      <rect x="38" y="4" width="24" height="14" rx="2" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="50" y="14" text-anchor="middle" fill="#38bdf8" font-size="7.5">0</text>
      <rect x="66" y="4" width="24" height="14" rx="2" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="78" y="14" text-anchor="middle" fill="#38bdf8" font-size="7.5">2</text>
      <rect x="94" y="4" width="24" height="14" rx="2" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="106" y="14" text-anchor="middle" fill="#38bdf8" font-size="7.5">3</text>
      <rect x="122" y="4" width="24" height="14" rx="2" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="134" y="14" text-anchor="middle" fill="#38bdf8" font-size="7.5">5</text>
      <text x="155" y="14" fill="#64748b" font-size="7.5">... (M+1)</text>

      <!-- indices array -->
      <text x="0" y="34" fill="#94a3b8" font-size="7.5" font-weight="700">indices:</text>
      <rect x="38" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="50" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5">1</text>
      <rect x="66" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="78" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5">3</text>
      <rect x="94" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="106" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5">2</text>
      <rect x="122" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="134" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5">0</text>
      <rect x="150" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="162" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5">3</text>

      <!-- data array -->
      <text x="0" y="54" fill="#94a3b8" font-size="7.5" font-weight="700">data   :</text>
      <rect x="38" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="50" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">4.2</text>
      <rect x="66" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="78" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">1.8</text>
      <rect x="94" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="106" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">9.1</text>
      <rect x="122" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="134" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">3.5</text>
      <rect x="150" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="162" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">7.0</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 압축 포인터: indptr [0, 2, 3, 5, ...]</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">i번째 세포의 0이 아닌 값 범위: indptr[i] ~ indptr[i+1]</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 세포별 슬라이싱: O(1) 포인터 즉시 접근</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">연속 메모리 블록 스트리밍 ➔ 64B 캐시 적중률 &gt; 96%</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자별 슬라이싱 병목 (Catastrophe):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">단일 유전자 추출 시 M개 모든 행 순회 ➔ 50x 지연!</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">세포 필터링 / 라이브러리 크기 / 세포 간 거리 계산의 제왕</text>
  </g>

  <!-- Panel 3: CSC Format -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. CSC 형식 (열/유전자 슬라이싱 최적화)</text>

    <!-- Schematic: CSC Structure -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- indptr col array -->
      <text x="0" y="14" fill="#94a3b8" font-size="7.5" font-weight="700">indptr :</text>
      <rect x="38" y="4" width="24" height="14" rx="2" fill="#059669" fill-opacity="0.4" stroke="#34d399" stroke-width="0.8"/>
      <text x="50" y="14" text-anchor="middle" fill="#34d399" font-size="7.5">0</text>
      <rect x="66" y="4" width="24" height="14" rx="2" fill="#059669" fill-opacity="0.4" stroke="#34d399" stroke-width="0.8"/>
      <text x="78" y="14" text-anchor="middle" fill="#34d399" font-size="7.5">1</text>
      <rect x="94" y="4" width="24" height="14" rx="2" fill="#059669" fill-opacity="0.4" stroke="#34d399" stroke-width="0.8"/>
      <text x="106" y="14" text-anchor="middle" fill="#34d399" font-size="7.5">2</text>
      <rect x="122" y="4" width="24" height="14" rx="2" fill="#059669" fill-opacity="0.4" stroke="#34d399" stroke-width="0.8"/>
      <text x="134" y="14" text-anchor="middle" fill="#34d399" font-size="7.5">3</text>
      <text x="155" y="14" fill="#64748b" font-size="7.5">... (N+1)</text>

      <!-- indices row array -->
      <text x="0" y="34" fill="#94a3b8" font-size="7.5" font-weight="700">indices:</text>
      <rect x="38" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
      <text x="50" y="34" text-anchor="middle" fill="#c084fc" font-size="7.5">2</text>
      <rect x="66" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
      <text x="78" y="34" text-anchor="middle" fill="#c084fc" font-size="7.5">0</text>
      <rect x="94" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
      <text x="106" y="34" text-anchor="middle" fill="#c084fc" font-size="7.5">1</text>
      <rect x="122" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
      <text x="134" y="34" text-anchor="middle" fill="#c084fc" font-size="7.5">0</text>
      <rect x="150" y="24" width="24" height="14" rx="2" fill="#1e293b" stroke="#c084fc" stroke-width="0.8"/>
      <text x="162" y="34" text-anchor="middle" fill="#c084fc" font-size="7.5">2</text>

      <!-- data array -->
      <text x="0" y="54" fill="#94a3b8" font-size="7.5" font-weight="700">data   :</text>
      <rect x="38" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="50" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">3.5</text>
      <rect x="66" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="78" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">4.2</text>
      <rect x="94" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="106" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">9.1</text>
      <rect x="122" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="134" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">1.8</text>
      <rect x="150" y="44" width="24" height="14" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="162" y="54" text-anchor="middle" fill="#34d399" font-size="7.5">7.0</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 열 압축 포인터: indptr 크기 N + 1</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">j번째 유전자의 발현 벡터를 O(1) 시간에 즉각 슬라이싱</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 고변동성 유전자(HVG) 선별 최적화:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">유전자별 분산/평균 계산 시 연속 메모리 고속 스캔</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 형식 변환 비용: CSR ↔ CSC</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">O(NNZ) 전치 소요 ➔ 파이프라인 단계별 사전 변환 필수</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">유전자별 마커 탐색 / DEG 회귀 분석 / HVG 선별의 핵심</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_hdf5_backed_mode_vs_zarr_chunking(filepath: str):
    """
    2. HDF5 Backed Mode vs Zarr Sharded Chunking & Concurrency Architecture
    Monolithic file lock vs directory-based sharded chunks with Blosc multi-threaded decompression.
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
    <marker id="arrRose02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrEmerald02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">디스크 기반 스토리지: HDF5 (.h5ad) vs Zarr Sharded Chunking &amp; 동시성 아키텍처</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">모놀리식 단일 파일 글로벌 락(HDF5 Lock) ➔ 디렉터리 기반 비동기 청크 샤딩 및 Blosc 병렬 압축 해제</text>

  <!-- Panel 1: HDF5 Limitations -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. HDF5 (.h5ad Backed)의 병목</text>

    <!-- Schematic: Single Monolithic File with Global Lock -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#e11d48" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Single file box -->
      <rect x="10" y="8" width="100" height="50" rx="4" fill="#1e293b" stroke="#f43f5e" stroke-width="1.2"/>
      <text x="60" y="26" text-anchor="middle" fill="#f43f5e" font-size="8" font-weight="700">dataset.h5ad</text>
      <text x="60" y="42" text-anchor="middle" fill="#94a3b8" font-size="6.8">(Single File Lock)</text>

      <!-- Lock icon / barrier -->
      <rect x="130" y="16" width="36" height="34" rx="4" fill="#e11d48" fill-opacity="0.3" stroke="#f43f5e" stroke-width="1"/>
      <text x="148" y="32" text-anchor="middle" fill="#f43f5e" font-size="12">🔒</text>
      <text x="148" y="44" text-anchor="middle" fill="#f43f5e" font-size="6.5">GIL Lock</text>

      <!-- Worker threads blocked -->
      <line x1="175" y1="20" x2="235" y2="20" stroke="#f43f5e" stroke-width="1.5" stroke-dasharray="2,2"/>
      <circle cx="238" cy="20" r="3" fill="#f43f5e"/>
      <text x="215" y="14" fill="#94a3b8" font-size="6.5">Worker 1 (Wait)</text>
      <line x1="175" y1="45" x2="235" y2="45" stroke="#f43f5e" stroke-width="1.5" stroke-dasharray="2,2"/>
      <circle cx="238" cy="45" r="3" fill="#f43f5e"/>
      <text x="215" y="40" fill="#94a3b8" font-size="6.5">Worker 2 (Wait)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• C 라이브러리 전역 락 (SWMR 한계):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">복수 워커 동시 접근 시 입출력 직렬화 강제 병목</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 클라우드 오브젝트 스토리지 부적합:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">S3/GCS에서 파일 일부 갱신 불가 (전체 재업로드)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• B-tree 메타데이터 탐색 지연:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">수십만 청크 접근 시 메타데이터 캐시 경합 발생</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.2" font-weight="600">위험: Dask / Ray 병렬 워커 확장 시 심각한 I/O 동결</text>
  </g>

  <!-- Panel 2: Zarr Sharded Chunking -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">2. Zarr 디렉터리 기반 비동기 청킹</text>

    <!-- Schematic: Sharded Chunks Directory -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Directory Tree -->
      <text x="5" y="14" fill="#38bdf8" font-size="7.5" font-weight="700">adata.zarr/X/</text>
      <!-- Chunks grid -->
      <rect x="95" y="6" width="35" height="22" rx="2" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="112" y="20" text-anchor="middle" fill="#38bdf8" font-size="7">c/0/0</text>
      <rect x="135" y="6" width="35" height="22" rx="2" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="152" y="20" text-anchor="middle" fill="#38bdf8" font-size="7">c/0/1</text>
      <rect x="95" y="34" width="35" height="22" rx="2" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="112" y="48" text-anchor="middle" fill="#38bdf8" font-size="7">c/1/0</text>
      <rect x="135" y="34" width="35" height="22" rx="2" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="152" y="48" text-anchor="middle" fill="#38bdf8" font-size="7">c/1/1</text>

      <!-- Parallel workers arrows -->
      <path d="M 175 17 L 235 12" stroke="#34d399" stroke-width="1.4" fill="none" marker-end="url(#arrEmerald02)"/>
      <text x="245" y="14" fill="#34d399" font-size="6.8">Worker A</text>
      <path d="M 175 45 L 235 48" stroke="#38bdf8" stroke-width="1.4" fill="none" marker-end="url(#arrCyan02)"/>
      <text x="245" y="50" fill="#38bdf8" font-size="6.8">Worker B</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• POSIX / S3 독립 파일 샤딩 (Zarr v3):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">각 청크가 별도 오브젝트로 저장되어 글로벌 락 제로</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 최적 청크 크기(Chunk Sizing): 16~64 MB</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">행(세포: 4,000) × 열(유전자: 2,000) 2차원 최적 타일링</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• Sharding Specification:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">작은 청크를 단일 파일 내 오프셋 샤드로 묶어 inode 보호</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">완전 비동기 병렬 읽기 지원 (Lock-Free I/O Architecture)</text>
  </g>

  <!-- Panel 3: Blosc Multithreaded Decompression -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. Blosc 압축 &amp; 초고속 메모리 대역폭</text>

    <!-- Schematic: Blosc pipeline -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Disk chunk -->
      <rect x="5" y="12" width="46" height="38" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="28" y="28" text-anchor="middle" fill="#fbbf24" font-size="6.5">NVMe Disk</text>
      <text x="28" y="40" text-anchor="middle" fill="#94a3b8" font-size="5.8">3.5 GB/s</text>

      <path d="M 53 31 L 78 31" stroke="#38bdf8" stroke-width="1.5" fill="none" marker-end="url(#arrCyan02)"/>

      <!-- Blosc Core -->
      <rect x="80" y="8" width="85" height="46" rx="4" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1.2"/>
      <text x="122" y="24" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">Blosc + BitShuffle</text>
      <text x="122" y="38" text-anchor="middle" fill="#cbd5e1" font-size="6.5">zstd / lz4 (8 Threads)</text>

      <path d="M 167 31 L 192 31" stroke="#34d399" stroke-width="1.5" fill="none" marker-end="url(#arrEmerald02)"/>

      <!-- L3 Cache / RAM -->
      <rect x="194" y="12" width="55" height="38" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
      <text x="221" y="28" text-anchor="middle" fill="#34d399" font-size="6.5">L3 Cache/RAM</text>
      <text x="221" y="40" text-anchor="middle" fill="#38bdf8" font-size="5.8">&gt; 12.0 GB/s</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 비트 셔플 전처리 (BitShuffle):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">비슷한 지수/비트를 연속 배치하여 압축률 3.5x 극대화</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 멀티스레드 캐시 내 압축 해제:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">L1/L2 캐시 블록 단위로 분할하여 전송과 동시에 해제</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유효 메모리 대역폭 역전 현상:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">디스크 I/O 속도(3.5GB/s)가 RAM 대역폭 수준(&gt;12GB/s) 체감</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">대규모 아틀라스 분석(CellxGene, scverse)의 표준 스토리지</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_tiledb_soma_sparse_tensor_architecture(filepath: str):
    """
    3. TileDB-SOMA Sparse Tensor Architecture & Tile-Based Storage
    Multi-array hierarchy (SOMAExperiment, SOMADataFrame, SOMASparseNDArray), tile/fragment storage & Arrow zero-copy.
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
    <marker id="arrPurple03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
    <marker id="arrEmerald03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">TileDB-SOMA 아키텍처: 희소 텐서 공간 인덱싱과 Apache Arrow 제로-카피 메모리</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">SOMA 계층 컬렉션 ➔ 타일 기반 2D 희소 배열(타일·프래그먼트) ➔ Arrow IPC Zero-Copy 상호운용성</text>

  <!-- Panel 1: SOMA Collection Hierarchy -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">1. SOMA 계층 구조 (Experiment Hierarchy)</text>

    <!-- Schematic: SOMA tree -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#a855f7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- SOMAExperiment -->
      <rect x="70" y="4" width="120" height="18" rx="3" fill="#a855f7" fill-opacity="0.3" stroke="#c084fc" stroke-width="1"/>
      <text x="130" y="16" text-anchor="middle" fill="#c084fc" font-size="7.5" font-weight="700">SOMAExperiment</text>

      <!-- Lines branching -->
      <path d="M 90 22 L 40 38" stroke="#94a3b8" stroke-width="1"/>
      <path d="M 130 22 L 130 38" stroke="#94a3b8" stroke-width="1"/>
      <path d="M 170 22 L 210 38" stroke="#94a3b8" stroke-width="1"/>

      <!-- Children nodes -->
      <rect x="10" y="38" width="60" height="22" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="40" y="49" text-anchor="middle" fill="#38bdf8" font-size="6.5" font-weight="700">obs</text>
      <text x="40" y="57" text-anchor="middle" fill="#94a3b8" font-size="5.5">(DataFrame)</text>

      <rect x="85" y="38" width="90" height="22" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="130" y="49" text-anchor="middle" fill="#fbbf24" font-size="6.5" font-weight="700">ms['RNA']</text>
      <text x="130" y="57" text-anchor="middle" fill="#94a3b8" font-size="5.5">(SOMAMeasurement)</text>

      <rect x="185" y="38" width="65" height="22" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="217" y="49" text-anchor="middle" fill="#34d399" font-size="6.5" font-weight="700">X['raw']</text>
      <text x="217" y="57" text-anchor="middle" fill="#94a3b8" font-size="5.5">(SparseNDArray)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• Single-cell Operations for Memory-friendly Apps:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">CZI 바이오허브 &amp; TileDB 연합 개방형 단일세포 표준 규격</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 완벽한 다중 오믹스 (Multi-Modal):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">단일 실험 내 RNA, ATAC, Spatial 등 복수 계측 동시 관리</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 독립적인 스키마 진화 (Evolvable):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">obs 메타데이터 추가 시 100GB 거대 카운트 행렬 재기록 불요</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#a855f7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="9.2" font-weight="600">글로벌 1억 세포 휴먼 셀 아틀라스(HCA) 핵심 인프라</text>
  </g>

  <!-- Panel 2: Tile-Based Sparse Array -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">2. 타일 기반 2D 희소 배열</text>

    <!-- Schematic: Tiles overlapping a query box -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Grid 4x4 background -->
      <line x1="20" y1="10" x2="100" y2="10" stroke="#334155" stroke-width="0.8"/>
      <line x1="20" y1="26" x2="100" y2="26" stroke="#334155" stroke-width="0.8"/>
      <line x1="20" y1="42" x2="100" y2="42" stroke="#334155" stroke-width="0.8"/>
      <line x1="20" y1="58" x2="100" y2="58" stroke="#334155" stroke-width="0.8"/>
      <line x1="20" y1="10" x2="20" y2="58" stroke="#334155" stroke-width="0.8"/>
      <line x1="46" y1="10" x2="46" y2="58" stroke="#334155" stroke-width="0.8"/>
      <line x1="73" y1="10" x2="73" y2="58" stroke="#334155" stroke-width="0.8"/>
      <line x1="100" y1="10" x2="100" y2="58" stroke="#334155" stroke-width="0.8"/>

      <!-- Tiles overlapping the query range (highlighted) -->
      <rect x="46" y="26" width="27" height="16" fill="#38bdf8" fill-opacity="0.35" stroke="#38bdf8" stroke-width="1.2"/>
      <rect x="73" y="26" width="27" height="16" fill="#38bdf8" fill-opacity="0.35" stroke="#38bdf8" stroke-width="1.2"/>
      <circle cx="33" cy="18" r="2.5" fill="#fbbf24"/>

      <!-- Explanatory arrows to 1D disk -->
      <path d="M 115 34 L 145 34" stroke="#fbbf24" stroke-width="1.4" fill="none" marker-end="url(#arrCyan03)"/>
      <rect x="150" y="16" width="105" height="36" rx="3" fill="#1e293b" stroke="#0284c7" stroke-width="0.8"/>
      <text x="202" y="30" text-anchor="middle" fill="#38bdf8" font-size="6.8" font-weight="700">겹치는 타일만 읽기</text>
      <text x="202" y="44" text-anchor="middle" fill="#94a3b8" font-size="6">비어 있지 않은 셀만 저장</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 타일 기반 희소 배열:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">2D 세포×유전자 좌표를 고정 크기 타일(기본 2048×2048)로 분할 저장</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 2D 범위 쿼리(Sub-slicing) 최적화:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">특정 환자군 × 유전자 서브셋 슬라이싱 시 겹치는 타일만 읽음</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 희소 타일 압축:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">차원·속성별 압축 필터(ZSTD, 델타 계열)로 저장 용량 절감</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">행·열 슬라이싱 성능은 타일 크기와 데이터 정렬에 의존</text>
  </g>

  <!-- Panel 3: Apache Arrow Zero-Copy -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. Apache Arrow 제로-카피(Zero-Copy)</text>

    <!-- Schematic: Zero-copy pointer bridge -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- TileDB C++ core -->
      <rect x="5" y="12" width="65" height="40" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="37" y="28" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">TileDB Core</text>
      <text x="37" y="42" text-anchor="middle" fill="#94a3b8" font-size="6">(C++ Engine)</text>

      <!-- Zero-copy arrow bridge -->
      <path d="M 72 32 L 105 32" stroke="#34d399" stroke-width="1.8" fill="none" marker-end="url(#arrEmerald03)"/>
      <text x="88" y="26" text-anchor="middle" fill="#34d399" font-size="6" font-weight="700">0-Copy</text>

      <!-- Shared memory buffer -->
      <rect x="108" y="8" width="70" height="48" rx="4" fill="#059669" fill-opacity="0.25" stroke="#34d399" stroke-width="1.2"/>
      <text x="143" y="24" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Arrow IPC</text>
      <text x="143" y="38" text-anchor="middle" fill="#cbd5e1" font-size="6">Shared Memory</text>
      <text x="143" y="49" text-anchor="middle" fill="#94a3b8" font-size="5.5">RecordBatch</text>

      <!-- Python / R bindings -->
      <path d="M 180 32 L 205 32" stroke="#34d399" stroke-width="1.8" fill="none" marker-end="url(#arrEmerald03)"/>
      <rect x="208" y="12" width="50" height="40" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="233" y="28" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">Python / R</text>
      <text x="233" y="42" text-anchor="middle" fill="#94a3b8" font-size="6">PyArrow/SciPy</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 메모리 직렬화 오버헤드 0%:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">C++ 엔진 버퍼를 파이썬 런타임이 메모리 복사 없이 포인터 참조</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 가비지 컬렉터(GC) 압박 원천 제거:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">수천만 개 파이썬 객체 생성 방지로 메모리 누수 및 일시정지 차단</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 언어 간 무결성 상호운용:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">Python(AnnData) ↔ R(SingleCellExperiment) 완벽 동기화</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">필요한 타일만 읽어 서브셋 쿼리 저지연(수치는 환경 의존)</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_dask_lazy_dag_out_of_core_normalization(filepath: str):
    """
    4. Dask Lazy DAG Out-of-Core Normalization & Incremental SVD Memory Profile
    Disk partition streaming -> Lazy task DAG -> Library size normalization & Low-rank QR SVD (Peak RAM < 16GB).
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
    <marker id="arrAmber04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrEmerald04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">Dask 지연 연산 DAG와 메모리 상한(&lt;16GB) Out-of-Core 정규화 &amp; 점진적 SVD</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">디스크 청크 스트리밍 ➔ 태스크 그래프(DAG) ➔ 분할 라이브러리 크기 정규화 &amp; 점진적 로우 랭크 QR 분해</text>

  <!-- Panel 1: Dask Lazy Task DAG -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. Dask 비순환 유향 그래프 (Task DAG)</text>

    <!-- Schematic: Task Graph Flow -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Step A: Read chunks -->
      <rect x="0" y="10" width="48" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="24" y="22" text-anchor="middle" fill="#38bdf8" font-size="6.5">Read(C1)</text>
      <rect x="0" y="34" width="48" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="24" y="46" text-anchor="middle" fill="#38bdf8" font-size="6.5">Read(C2)</text>

      <path d="M 48 19 L 68 25" stroke="#94a3b8" stroke-width="1" marker-end="url(#arrCyan04)"/>
      <path d="M 48 43 L 68 37" stroke="#94a3b8" stroke-width="1" marker-end="url(#arrCyan04)"/>

      <!-- Step B: Sum & Scale -->
      <rect x="70" y="22" width="60" height="20" rx="3" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="100" y="35" text-anchor="middle" fill="#38bdf8" font-size="6.5" font-weight="700">LibSizeSum</text>

      <path d="M 130 32 L 150 32" stroke="#94a3b8" stroke-width="1" marker-end="url(#arrCyan04)"/>

      <!-- Step C: Norm & Log1p -->
      <rect x="152" y="14" width="46" height="34" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="175" y="28" text-anchor="middle" fill="#fbbf24" font-size="6.5">Scale &amp;</text>
      <text x="175" y="40" text-anchor="middle" fill="#fbbf24" font-size="6.5">Log1p</text>

      <path d="M 198 32 L 216 32" stroke="#94a3b8" stroke-width="1" marker-end="url(#arrEmerald04)"/>

      <!-- Step D: Incremental QR -->
      <rect x="218" y="14" width="42" height="34" rx="3" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1"/>
      <text x="239" y="28" text-anchor="middle" fill="#34d399" font-size="6.5" font-weight="700">Inc.</text>
      <text x="239" y="40" text-anchor="middle" fill="#34d399" font-size="6.5" font-weight="700">SVD</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 지연 평가 (Lazy Execution):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">코드 실행 시 계산하지 않고 계산 트리(DAG)만 정의</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 토폴로지 최적 스케줄링:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">메모리 방출(Garbage Collection) 우선순위 워커 할당</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• dask.array.map_blocks 연산:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">각 청크별 희소 연산자(SpMV)를 독립 프로세스 병렬 수행</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">수백 기가바이트 단일세포 행렬을 랩톱에서도 무결 처리</text>
  </g>

  <!-- Panel 2: Memory Profile Under 16GB Limit -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 엄격한 메모리 상한 (&lt; 16GB Flat Profile)</text>

    <!-- Schematic: RAM usage curve over time -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Axes -->
      <line x1="20" y1="58" x2="255" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="58" x2="20" y2="8" stroke="#475569" stroke-width="1"/>
      <text x="5" y="16" fill="#f43f5e" font-size="6">128GB</text>
      <text x="5" y="38" fill="#34d399" font-size="6">16GB</text>
      <text x="240" y="66" fill="#64748b" font-size="6">Time</text>

      <!-- In-memory naive crash curve (dashed red) -->
      <path d="M 20 54 Q 60 50 80 30 Q 95 10 100 8" fill="none" stroke="#f43f5e" stroke-width="1.5" stroke-dasharray="3,3"/>
      <circle cx="100" cy="8" r="3" fill="#f43f5e"/>
      <text x="105" y="12" fill="#f43f5e" font-size="6.5" font-weight="700">OOM Crash!</text>

      <!-- Out-of-core Flat curve (solid green oscillating below 16GB) -->
      <path d="M 20 55 L 40 40 L 55 50 L 75 38 L 90 48 L 110 39 L 130 49 L 150 40 L 170 48 L 195 39 L 220 48 L 245 40" fill="none" stroke="#34d399" stroke-width="2"/>
      <!-- Safe limit line -->
      <line x1="20" y1="36" x2="255" y2="36" stroke="#fbbf24" stroke-width="1" stroke-dasharray="2,2"/>
      <text x="160" y="32" fill="#fbbf24" font-size="6.5">Peak RAM Budget (16GB)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 인메모리(In-Memory) 방식의 OOM 폭발:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">100만 세포 로드 시 128GB 이상 메모리 요구 ➔ 프로세스 사망</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 청크 스트리밍(Chunk Streaming):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">한 번에 2~4GB 분할만 활성 로드 ➔ 연산 즉시 메모리 반납</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 메모리 역치 기반 백프레셔(Backpressure):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">워커 메모리가 75% 도달 시 디스크 읽기 속도 자동 조절</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">피크 RAM 사용량을 12GB 이하로 일정하게 고정 유지</text>
  </g>

  <!-- Panel 3: Incremental SVD (Low-Rank QR) -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 점진적 SVD (Incremental PCA / QR)</text>

    <!-- Schematic: QR Update step -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Existing Basis Q -->
      <rect x="5" y="10" width="45" height="42" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="27" y="26" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">기존 기저</text>
      <text x="27" y="40" text-anchor="middle" fill="#94a3b8" font-size="6.5">Q (k=50)</text>

      <text x="60" y="34" fill="#fbbf24" font-size="12" font-weight="700">⊕</text>

      <!-- New Chunk X_i -->
      <rect x="75" y="10" width="55" height="42" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="102" y="26" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">신규 청크</text>
      <text x="102" y="40" text-anchor="middle" fill="#94a3b8" font-size="6.5">X_i (ΔCells)</text>

      <path d="M 138 32 L 165 32" stroke="#34d399" stroke-width="1.5" fill="none" marker-end="url(#arrEmerald04)"/>

      <!-- QR Decomposition block -->
      <rect x="170" y="8" width="85" height="46" rx="3" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1"/>
      <text x="212" y="24" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Low-Rank QR</text>
      <text x="212" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">SVD on R (k×k)</text>
      <text x="212" y="48" text-anchor="middle" fill="#38bdf8" font-size="5.8">Q_new = [Q P]·U'</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 직교 투영 잔차 분해 (Gram-Schmidt QR):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">신규 청크를 기존 기저 Q에 투영 후 잔차 벡터 직교화</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 연산 복잡도 O(k² + k·N):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">전체 세포 수 M과 무관하며, 오직 PC 수 k와 유전자 수 N에 비례</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• Full SVD 대비 완벽한 일치도:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">상위 50개 주성분 코사인 유사도 R &gt; 0.999 검증</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">점진적 업데이트를 통해 100만 세포 차원 축소를 단일 패스로 완료</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_rapids_gpu_acceleration_cuda_pipeline(filepath: str):
    """
    5. RAPIDS cuGraph & cuML GPU Acceleration CUDA Pipeline
    Host CPU RAM -> PCIe Gen4/5 -> GPU HBM streaming and cuML GPU speedup (illustrative numbers).
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
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">RAPIDS GPU 전사체 파이프라인: PCIe Gen4/5 전송과 cuML GPU 가속 (예시 수치)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">호스트 CPU RAM ➔ 고속 PCIe 버스 ➔ GPU HBM3 스트리밍 및 cuML PCA/k-NN/UMAP 병렬 연산</text>

  <!-- Panel 1: Memory Hierarchy & PCIe Bus -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 이기종 메모리 계층 &amp; PCIe 대역폭</text>

    <!-- Schematic: CPU RAM -> PCIe -> GPU HBM -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- CPU Host RAM -->
      <rect x="0" y="8" width="60" height="46" rx="3" fill="#1e293b" stroke="#94a3b8" stroke-width="1"/>
      <text x="30" y="24" text-anchor="middle" fill="#cbd5e1" font-size="7" font-weight="700">Host RAM</text>
      <text x="30" y="36" text-anchor="middle" fill="#94a3b8" font-size="6">DDR5</text>
      <text x="30" y="46" text-anchor="middle" fill="#38bdf8" font-size="5.8">~200 GB/s</text>

      <!-- PCIe Gen5 Bridge -->
      <path d="M 64 31 L 96 31" stroke="#fbbf24" stroke-width="2" fill="none" marker-end="url(#arrAmber05)"/>
      <rect x="98" y="16" width="66" height="30" rx="3" fill="#d97706" fill-opacity="0.25" stroke="#fbbf24" stroke-width="1"/>
      <text x="131" y="29" text-anchor="middle" fill="#fbbf24" font-size="6.5" font-weight="700">PCIe Gen5 x16</text>
      <text x="131" y="40" text-anchor="middle" fill="#cbd5e1" font-size="6">64 GB/s (Uni)</text>

      <path d="M 168 31 L 198 31" stroke="#34d399" stroke-width="2" fill="none" marker-end="url(#arrEmerald05)"/>

      <!-- GPU HBM3 -->
      <rect x="202" y="8" width="64" height="46" rx="3" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1.2"/>
      <text x="234" y="24" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">GPU HBM3</text>
      <text x="234" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">H100 80GB</text>
      <text x="234" y="46" text-anchor="middle" fill="#34d399" font-size="5.8">3,350 GB/s</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 메모리 대역폭의 16배 도약:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Host DDR5(200GB/s) 대비 HBM3(3,350GB/s)의 압도적 속도</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• Pinned (Page-Locked) Memory:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">OS 페이징을 차단하고 GPU Direct DMA 전송 극대화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 비동기 CUDA 스트림 오버랩:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">청크 N+1 복사와 청크 N의 연산을 완벽 동시 수행</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">데이터 전송 병목을 연산 시간 뒤로 완전히 은폐</text>
  </g>

  <!-- Panel 2: cuML Tensor Core Workflow -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. cuML &amp; GPU 파이프라인</text>

    <!-- Schematic: cuML pipeline flow -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Stage 1: cuSparse / TruncatedSVD -->
      <rect x="5" y="10" width="70" height="44" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="40" y="24" text-anchor="middle" fill="#38bdf8" font-size="6.8" font-weight="700">cuml.PCA</text>
      <text x="40" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">TruncatedSVD</text>
      <text x="40" y="47" text-anchor="middle" fill="#94a3b8" font-size="5.5">(GEMM / eigh)</text>

      <path d="M 78 32 L 95 32" stroke="#fbbf24" stroke-width="1.5" fill="none" marker-end="url(#arrAmber05)"/>

      <!-- Stage 2: NearestNeighbors -->
      <rect x="98" y="10" width="75" height="44" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="135" y="24" text-anchor="middle" fill="#fbbf24" font-size="6.8" font-weight="700">cuml.KNN</text>
      <text x="135" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">Brute/IVF-PQ</text>
      <text x="135" y="47" text-anchor="middle" fill="#94a3b8" font-size="5.5">(k=30 L2 Graph)</text>

      <path d="M 176 32 L 193 32" stroke="#34d399" stroke-width="1.5" fill="none" marker-end="url(#arrEmerald05)"/>

      <!-- Stage 3: cuML UMAP -->
      <rect x="196" y="10" width="68" height="44" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
      <text x="230" y="24" text-anchor="middle" fill="#34d399" font-size="6.8" font-weight="700">cuml.UMAP</text>
      <text x="230" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">Parallel SGD</text>
      <text x="230" y="47" text-anchor="middle" fill="#94a3b8" font-size="5.5">(2D Manifold)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 밀집 행렬 곱셈 (GEMM) 기반 PCA:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">공분산 계산·고유분해를 GPU에서 병렬 처리</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 병렬 k-NN 그래프 생성:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">수천 개 CUDA 스레드가 세포 쌍별 거리 동시 계산</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• GPU 비동기 확률적 경사하강법 (UMAP SGD):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">원자적 잠금 없는 Hogwild 방식 가속으로 2D 투영</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">Scanpy와 유사한 API를 제공하는 GPU 가속 라이브러리</text>
  </g>

  <!-- Panel 3: Speedup Benchmark (1M Cells) -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 100만 세포 런타임 비교 (설명용 예시 수치)</text>

    <!-- Schematic: Bar chart comparison -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Item 1: PCA -->
      <text x="0" y="15" fill="#cbd5e1" font-size="6.8" font-weight="700">PCA (50 PC):</text>
      <!-- CPU bar -->
      <rect x="65" y="8" width="130" height="6" rx="1.5" fill="#f43f5e"/>
      <text x="200" y="14" fill="#f43f5e" font-size="6">CPU 180s</text>
      <!-- GPU bar -->
      <rect x="65" y="16" width="3" height="6" rx="1" fill="#34d399"/>
      <text x="75" y="22" fill="#34d399" font-size="6" font-weight="700">GPU 1.8s</text>

      <!-- Item 2: k-NN (k=30) -->
      <text x="0" y="37" fill="#cbd5e1" font-size="6.8" font-weight="700">k-NN Graph :</text>
      <!-- CPU bar -->
      <rect x="65" y="30" width="145" height="6" rx="1.5" fill="#f43f5e"/>
      <text x="215" y="36" fill="#f43f5e" font-size="6">CPU 950s</text>
      <!-- GPU bar -->
      <rect x="65" y="38" width="4" height="6" rx="1" fill="#34d399"/>
      <text x="75" y="44" fill="#34d399" font-size="6" font-weight="700">GPU 9.2s</text>

      <!-- Item 3: UMAP -->
      <text x="0" y="59" fill="#cbd5e1" font-size="6.8" font-weight="700">UMAP (2D)  :</text>
      <!-- CPU bar -->
      <rect x="65" y="52" width="160" height="6" rx="1.5" fill="#f43f5e"/>
      <text x="230" y="58" fill="#f43f5e" font-size="6">CPU 1,400s</text>
      <!-- GPU bar -->
      <rect x="65" y="60" width="5" height="6" rx="1" fill="#34d399"/>
      <text x="75" y="66" fill="#34d399" font-size="6" font-weight="700">GPU 14.5s</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 전체 분석 파이프라인 시간 단축:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">예시: CPU 수십 분 ➔ GPU 수십 초 (조건 의존)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 대화형 탐색 (Interactive Latency):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">파라미터 변경(e.g., k값, 해상도) 후 즉시 실시간 재시각화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 클라우드 비용 절감 가능성:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">소요 시간 단축에 따른 인스턴스 사용 시간 감소</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">100만 세포 분석 시간을 크게 단축하는 GPU 가속</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_cagra_gpu_knn_and_cugraph_leiden(filepath: str):
    """
    6. NVIDIA CAGRA GPU k-NN Graph Traversal & cuGraph Parallel Leiden Clustering
    Approximate graph search -> Warp-level memory coalescing -> Atomic community aggregation on 1M cells.
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
    <marker id="arrPurple06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
    <marker id="arrEmerald06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">NVIDIA CAGRA 초고속 근사 k-NN 그래프 탐색 &amp; cuGraph 원자적 병렬 라이덴 군집화</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">GPU 압축 근사 그래프(CAGRA) ➔ 워프 레벨 병렬 탐색 ➔ cuGraph 원자적 커뮤니티 병합 &amp; 100만 세포 분할</text>

  <!-- Panel 1: NVIDIA CAGRA Graph Traversal -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. NVIDIA CAGRA 근사 탐색 구조</text>

    <!-- Schematic: CAGRA Small-degree graph -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Query node -->
      <circle cx="25" cy="32" r="7" fill="#fbbf24"/>
      <text x="25" y="35" text-anchor="middle" fill="#0f172a" font-size="6.5" font-weight="700">Q</text>

      <!-- Edges to first hop -->
      <line x1="32" y1="32" x2="80" y2="16" stroke="#38bdf8" stroke-width="1.2"/>
      <line x1="32" y1="32" x2="80" y2="48" stroke="#38bdf8" stroke-width="1.2"/>

      <!-- Hop 1 nodes -->
      <circle cx="85" cy="16" r="6" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
      <circle cx="85" cy="48" r="6" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>

      <!-- Edges to hop 2 -->
      <line x1="91" y1="16" x2="140" y2="10" stroke="#34d399" stroke-width="1"/>
      <line x1="91" y1="16" x2="140" y2="28" stroke="#34d399" stroke-width="1"/>
      <line x1="91" y1="48" x2="140" y2="44" stroke="#34d399" stroke-width="1"/>
      <line x1="91" y1="48" x2="140" y2="58" stroke="#34d399" stroke-width="1"/>

      <!-- Hop 2 nodes -->
      <circle cx="145" cy="10" r="4.5" fill="#34d399"/>
      <circle cx="145" cy="28" r="4.5" fill="#34d399"/>
      <circle cx="145" cy="44" r="4.5" fill="#34d399"/>
      <circle cx="145" cy="58" r="4.5" fill="#34d399"/>

      <!-- Label block -->
      <rect x="165" y="14" width="95" height="38" rx="3" fill="#1e293b" stroke="#0284c7" stroke-width="0.8"/>
      <text x="212" y="27" text-anchor="middle" fill="#38bdf8" font-size="6.5" font-weight="700">Small-Degree Graph</text>
      <text x="212" y="38" text-anchor="middle" fill="#cbd5e1" font-size="6">Fixed Out-Degree (d=32)</text>
      <text x="212" y="47" text-anchor="middle" fill="#34d399" font-size="5.5">Rank-based pruning</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• CAGRA 고정 차수 그래프 k-NN:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">GPU 하드웨어 구조에 맞춤 설계된 최첨단 벡터 색인</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 고정 차수(Degree) 제약 설계:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">모든 노드가 동일한 수의 이웃을 가져 GPU 워프 발산 감소</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• CPU HNSW 대비 속도 혁신:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">대형 배치 처리량 HNSW 대비 33~77배 (Ootomo et al. 2024)</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">초대규모 세포 데이터베이스 실시간 벡터 검색 엔진</text>
  </g>

  <!-- Panel 2: Warp-Level Memory Coalescing -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 워프 레벨 메모리 병합 (Coalesced)</text>

    <!-- Schematic: 32 threads warp evaluating distances -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#a855f7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Warp rectangle -->
      <rect x="5" y="8" width="115" height="48" rx="3" fill="#a855f7" fill-opacity="0.25" stroke="#c084fc" stroke-width="1"/>
      <text x="62" y="22" text-anchor="middle" fill="#c084fc" font-size="7" font-weight="700">CUDA Warp (32 Threads)</text>
      <!-- Small thread indicators -->
      <g transform="translate(15, 28)">
        <rect x="0" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <rect x="11" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <rect x="22" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <rect x="33" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <text x="50" y="9" fill="#94a3b8" font-size="6">...</text>
        <rect x="62" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <rect x="73" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
        <rect x="84" y="0" width="8" height="12" rx="1" fill="#38bdf8"/>
      </g>

      <path d="M 125 32 L 155 32" stroke="#34d399" stroke-width="1.8" fill="none" marker-end="url(#arrEmerald03)"/>

      <!-- Shared memory priority queue -->
      <rect x="160" y="8" width="105" height="48" rx="3" fill="#059669" fill-opacity="0.25" stroke="#34d399" stroke-width="1"/>
      <text x="212" y="24" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Shared Memory Queue</text>
      <text x="212" y="38" text-anchor="middle" fill="#cbd5e1" font-size="6">Top-k Candidate Heap</text>
      <text x="212" y="48" text-anchor="middle" fill="#38bdf8" font-size="5.5">Zero Global Memory Stall</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 단일 쿼리 ➔ 워프(32 Threads) 합동 할당:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">후보 노드의 32차원 특성을 동시 병합 로드(Coalesced Access)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 공유 메모리(Shared Memory) 우선순위 큐:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">칩 내부 초고속 SRAM에서 k개 최단 거리 힙 정렬</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 글로벌 메모리 지연 은폐 (Latency Hiding):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">DRAM 접근 대기 시간 동안 다른 워프가 연산 파이프라인 점유</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#a855f7" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#c084fc" font-size="9.2" font-weight="600">하드웨어 아키텍처 수준의 극한 메모리 대역폭 활용</text>
  </g>

  <!-- Panel 3: cuGraph Parallel Leiden -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. cuGraph 원자적 병렬 라이덴(Leiden)</text>

    <!-- Schematic: Atomic aggregation & refinement -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Unclustered cells -->
      <g transform="translate(10, 8)">
        <circle cx="10" cy="12" r="4.5" fill="#38bdf8"/>
        <circle cx="22" cy="18" r="4.5" fill="#38bdf8"/>
        <circle cx="15" cy="30" r="4.5" fill="#38bdf8"/>
        <circle cx="32" cy="36" r="4.5" fill="#fbbf24"/>
        <circle cx="42" cy="22" r="4.5" fill="#fbbf24"/>
        <text x="24" y="46" text-anchor="middle" fill="#94a3b8" font-size="5.5">1M Cells</text>
      </g>

      <path d="M 65 32 L 95 32" stroke="#fbbf24" stroke-width="1.5" fill="none" marker-end="url(#arrAmber06)"/>

      <!-- Atomic Modularity Kernel -->
      <rect x="98" y="10" width="80" height="44" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="138" y="24" text-anchor="middle" fill="#fbbf24" font-size="6.8" font-weight="700">Atomic Kernel</text>
      <text x="138" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">atomicCAS / Add</text>
      <text x="138" y="47" text-anchor="middle" fill="#38bdf8" font-size="5.5">Modularity Move</text>

      <path d="M 183 32 L 205 32" stroke="#34d399" stroke-width="1.5" fill="none" marker-end="url(#arrEmerald06)"/>

      <!-- Partitioned communities -->
      <rect x="208" y="8" width="60" height="48" rx="3" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1.2"/>
      <text x="238" y="24" text-anchor="middle" fill="#34d399" font-size="6.8" font-weight="700">Clusters</text>
      <text x="238" y="36" text-anchor="middle" fill="#cbd5e1" font-size="6">Connected</text>
      <text x="238" y="48" text-anchor="middle" fill="#fbbf24" font-size="5.8">Converge</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 원자적 커뮤니티 이동 (Atomic CAS):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">동시 다발 노드 이동 시 충돌 없는 모듈래리티 이득 계산</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 라이덴 리파인먼트 (Refinement) 보장:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">루뱅(Louvain)에서 생길 수 있는 연결성 단절 커뮤니티 방지</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 100만 세포 GPU 군집화 가속:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">CPU 대비 크게 단축 (소요 시간은 하드웨어·그래프에 의존)</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">γ-연결성을 보장하는 병렬 파티셔닝</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_sparse_matrix_csr_csc_memory_layout.svg", make_01_sparse_matrix_csr_csc_memory_layout),
        ("02_hdf5_backed_mode_vs_zarr_chunking.svg", make_02_hdf5_backed_mode_vs_zarr_chunking),
        ("03_tiledb_soma_sparse_tensor_architecture.svg", make_03_tiledb_soma_sparse_tensor_architecture),
        ("04_dask_lazy_dag_out_of_core_normalization.svg", make_04_dask_lazy_dag_out_of_core_normalization),
        ("05_rapids_gpu_acceleration_cuda_pipeline.svg", make_05_rapids_gpu_acceleration_cuda_pipeline),
        ("06_cagra_gpu_knn_and_cugraph_leiden.svg", make_06_cagra_gpu_knn_and_cugraph_leiden),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 08 Large-Scale Sparse publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

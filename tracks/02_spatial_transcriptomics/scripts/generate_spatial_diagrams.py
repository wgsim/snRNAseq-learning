#!/usr/bin/env python3
"""
tracks/02_spatial_transcriptomics/scripts/generate_spatial_diagrams.py
=====================================================================
Generates publication-grade, high-resolution SVG diagrams for Track 02:
Spatial Transcriptomics & In Situ Subcellular Biology in Prostate Cancer.

Covers:
1. Spatial Technologies Resolution Comparison (Visium vs Visium HD vs Xenium)
2. Coordinate Registration and SpatialData Architecture
3. Single-Cell Segmentation: Watershed vs Cellpose vs Baysor
4. Cell2location Bayesian Spatial Deconvolution Model
5. Moran's I Spatial Autocorrelation & LISA Scatterplot
6. Cellular Neighborhood Niche & Distance-Decay Receptor-Ligand Signaling
"""

import os
import xml.etree.ElementTree as ET


def make_01_spatial_technologies_resolution(filepath: str):
    """
    1. Spatial Technologies Resolution Comparison
    Visium (55µm spot) vs Visium HD (2µm square bins) vs Xenium (<0.1µm subcellular single molecules)
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="visiumGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="100%" stop-color="#0284c7"/>
    </linearGradient>
    <linearGradient id="hdGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#4ade80"/>
      <stop offset="100%" stop-color="#16a34a"/>
    </linearGradient>
    <linearGradient id="xeniumGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#c084fc"/>
      <stop offset="100%" stop-color="#9333ea"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">공간 전사체 기술별 해상도 &amp; 캡처 아키텍처 비교</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">NGS 시퀀싱 기반 다세포 스팟(Visium) vs 초고해상도 비닝(Visium HD) vs 단일 분자 광학 이미징(Xenium In Situ)</text>

  <!-- Panel 1: 10x Visium Standard -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 10x Visium (55 µm Spot)</text>

    <!-- Schematic: 55um circular spot with mixed cells -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Outer glass boundary -->
    <circle cx="65" cy="79" r="28" fill="#0284c7" fill-opacity="0.2" stroke="#38bdf8" stroke-width="2" stroke-dasharray="3,2"/>
    <text x="65" y="60" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">55 µm Spot</text>
    <!-- Multi-cells inside spot -->
    <circle cx="55" cy="75" r="7" fill="#f43f5e" fill-opacity="0.75"/>
    <circle cx="72" cy="72" r="6" fill="#c084fc" fill-opacity="0.75"/>
    <circle cx="60" cy="89" r="6.5" fill="#4ade80" fill-opacity="0.75"/>
    <circle cx="77" cy="87" r="5.5" fill="#f59e0b" fill-opacity="0.75"/>
    <text x="110" y="66" fill="#f43f5e" font-size="10" font-weight="700">다세포 혼합체 (Mixture)</text>
    <text x="110" y="80" fill="#cbd5e1" font-size="9.5">스팟당 5~15개 이종 세포 포함</text>
    <text x="110" y="94" fill="#94a3b8" font-size="9">Poly-dT 포획 프로브 (전장 전사체)</text>

    <!-- Technical Specs -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 스팟 직경 및 간격(Pitch):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">Ø 55 µm, 중심 간 거리 100 µm (공백 존재)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 캡처 방식 &amp; 스루풋:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">6.5 x 6.5 mm 영역 내 약 5,000개 스팟</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 계산적 필수 단계:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">Cell2location 등 참조 기반 디컨볼루션 필수</text>

    <!-- Bottom Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">장점: 무편향 전장 유전체 (Whole Transcriptome)</text>
  </g>

  <!-- Panel 2: 10x Visium HD -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#4ade80" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#4ade80" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#4ade80" font-size="13" font-weight="700">2. Visium HD (2 µm Bins)</text>

    <!-- Schematic: Continuous gapless 2um raster -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#16a34a" stroke-width="1"/>
    <!-- 2um grid raster -->
    <g transform="translate(35, 52)">
      <rect x="0" y="0" width="54" height="54" fill="#0f172a" stroke="#4ade80" stroke-width="1.2"/>
      <line x1="18" y1="0" x2="18" y2="54" stroke="#4ade80" stroke-width="0.8" stroke-dasharray="2,2"/>
      <line x1="36" y1="0" x2="36" y2="54" stroke="#4ade80" stroke-width="0.8" stroke-dasharray="2,2"/>
      <line x1="0" y1="18" x2="54" y2="18" stroke="#4ade80" stroke-width="0.8" stroke-dasharray="2,2"/>
      <line x1="0" y1="36" x2="54" y2="36" stroke="#4ade80" stroke-width="0.8" stroke-dasharray="2,2"/>
      <!-- Cell overlay on 2um grid -->
      <circle cx="27" cy="27" r="16" fill="#4ade80" fill-opacity="0.25" stroke="#4ade80" stroke-width="1.5"/>
      <text x="27" y="30" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">Single Cell</text>
    </g>
    <text x="110" y="66" fill="#4ade80" font-size="10" font-weight="700">연속 격자 (Continuous Grid)</text>
    <text x="110" y="80" fill="#cbd5e1" font-size="9.5">빈 공간 없는 2 µm x 2 µm 마이크로 타일</text>
    <text x="110" y="94" fill="#94a3b8" font-size="9">8 µm / 16 µm 동적 비닝 및 단일 세포 정렬</text>

    <!-- Technical Specs -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 해상도 및 비닝 체계:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">기본 2 µm 격자 ➔ 8 µm (단일 세포급) 집계</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 정보 손실 최소화:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">스팟 간 공백(Gap) 제거 ➔ 미세 침윤 경계 보존</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 전사체 포획 효율:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">FFPE 조직 호환 및 1,800만+ 개 2 µm 피처</text>

    <!-- Bottom Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#16a34a" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#4ade80" font-size="10" font-weight="600">장점: 전장 전사체 + 준단일세포(Sub-cellular) 경계</text>
  </g>

  <!-- Panel 3: 10x Xenium In Situ -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. Xenium In Situ (&lt;0.1 µm)</text>

    <!-- Schematic: Subcellular puncta with cell boundary -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#9333ea" stroke-width="1"/>
    <!-- Cell membrane and nucleus -->
    <path d="M 30 75 Q 40 52 65 55 Q 90 60 85 88 Q 65 105 40 95 Z" fill="#9333ea" fill-opacity="0.15" stroke="#c084fc" stroke-width="1.2"/>
    <ellipse cx="58" cy="76" rx="12" ry="9" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="1"/>
    <!-- Individual RNA puncta -->
    <circle cx="56" cy="74" r="1.8" fill="#f43f5e"/>
    <circle cx="62" cy="78" r="1.8" fill="#f43f5e"/>
    <circle cx="45" cy="68" r="1.8" fill="#4ade80"/>
    <circle cx="72" cy="65" r="1.8" fill="#4ade80"/>
    <circle cx="75" cy="85" r="1.8" fill="#f59e0b"/>
    <circle cx="42" cy="88" r="1.8" fill="#f59e0b"/>
    <text x="110" y="66" fill="#c084fc" font-size="10" font-weight="700">단일 분자 광학 디코딩</text>
    <text x="110" y="80" fill="#cbd5e1" font-size="9.5">개별 mRNA 분자 (x, y, z) 정밀 좌표</text>
    <text x="110" y="94" fill="#94a3b8" font-size="9">Padlock 프로브 + 회전환 증폭(RCA)</text>

    <!-- Technical Specs -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 광학 회절 한계 해상도:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">&lt; 0.1 µm 수준의 세포내 국소화 (Subcellular)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 타깃 유전자 패널:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">300 ~ 5,000+ 개 큐레이션 타깃 (고감도 검출)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 세포 경계 분할(Segmentation):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">DAPI/Cellpose/Baysor 기반 단일 세포 할당</text>

    <!-- Bottom Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#9333ea" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">장점: 100% 진정한 단일 세포 + 세포내 구획화</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_coordinate_registration_spatialdata(filepath: str):
    """
    2. Coordinate Registration and SpatialData Architecture
    H&E histology pixel space -> Fiducial frame affine alignment -> Micron space -> SpatialData OME-Zarr representation
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad02" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">멀티스케일 좌표계 등록 &amp; SpatialData OME-Zarr 아키텍처</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">조직학 H&amp;E 픽셀 좌표 ➔ 기준 마커(Fiducial) Affine 정렬 ➔ 물리적 마이크론(µm) 공간 ➔ 계층적 OME-Zarr 통합</text>

  <!-- Panel 1: Pixel Space & Fiducial Alignment -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 조직학 픽셀 공간 (Pixel Space)</text>

    <!-- Schematic: Slide with fiducial corners -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <rect x="35" y="52" width="60" height="54" fill="#334155" fill-opacity="0.4" stroke="#64748b" stroke-width="1"/>
    <!-- Fiducial markers on 4 corners -->
    <circle cx="40" cy="57" r="3" fill="#f43f5e"/>
    <circle cx="90" cy="57" r="3" fill="#f43f5e"/>
    <circle cx="40" cy="101" r="3" fill="#f43f5e"/>
    <circle cx="90" cy="101" r="3" fill="#f43f5e"/>
    <!-- Tissue blob in middle -->
    <ellipse cx="65" cy="79" rx="18" ry="14" fill="#f43f5e" fill-opacity="0.25" stroke="#f43f5e" stroke-width="1"/>
    <text x="115" y="66" fill="#38bdf8" font-size="10" font-weight="700">초고해상도 WSI 이미지</text>
    <text x="115" y="80" fill="#cbd5e1" font-size="9.5">Full-res H&amp;E / DAPI (px 좌표)</text>
    <text x="115" y="94" fill="#f43f5e" font-size="9">기준 마커: 4-Corner Fiducials</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 스캔 배율 및 픽셀 피치:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">20x/40x 대물렌즈 (0.25 ~ 0.5 µm/pixel)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 광학 왜곡 및 회전 발생:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">슬라이드 장착 편차로 인한 미세 각도 회전</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 피듀셜 프레임 매핑:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">알려진 슬라이드 마커 기반 기준 좌표 추출</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">입력: (x_px, y_px) 원시 카메라 이미지</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 345 200" stroke="#38bdf8" stroke-width="2"/>
  <polygon points="345,196 350,200 345,204" fill="#38bdf8"/>

  <!-- Panel 2: Affine Transformation to Micron Space -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#f59e0b" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f59e0b" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f59e0b" font-size="13" font-weight="700">2. 어파인 변환 &amp; 마이크론 정렬</text>

    <!-- Schematic: Affine Matrix Transform -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <!-- Formula block -->
    <text x="150" y="66" text-anchor="middle" fill="#f59e0b" font-size="11" font-weight="700">Affine Transformation Matrix A</text>
    <text x="150" y="84" text-anchor="middle" fill="#cbd5e1" font-size="10" font-family="monospace">[x_µm, y_µm, 1]ᵀ = A · [x_px, y_px, 1]ᵀ</text>
    <text x="150" y="102" text-anchor="middle" fill="#4ade80" font-size="9">스케일(s), 회전(θ), 전단(sh), 이동(tx, ty) 보정</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 전역 아핀 정렬 (Global Affine):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">Fiducial 4점을 목표 좌표계로 최소자승 피팅</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 물리적 마이크론 단위 변환:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">픽셀 인덱스를 실제 조직 물리 좌표(µm)로 투영</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 다층 연속 절편(Serial Slice) 정합:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">연속 절편 간 탄성 변형(Elastic B-spline) 교정</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">출력: 공통 물리 공간 (Physical Micron Space)</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 655 200 L 670 200" stroke="#f59e0b" stroke-width="2"/>
  <polygon points="670,196 675,200 670,204" fill="#f59e0b"/>

  <!-- Panel 3: SpatialData & OME-Zarr Framework -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#4ade80" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#4ade80" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#4ade80" font-size="13" font-weight="700">3. SpatialData &amp; OME-Zarr 표현</text>

    <!-- Schematic: SpatialData 4 container elements -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#16a34a" stroke-width="1"/>
    <rect x="25" y="49" width="55" height="26" rx="3" fill="#0284c7" fill-opacity="0.3"/>
    <text x="52" y="66" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">Images</text>
    <rect x="85" y="49" width="55" height="26" rx="3" fill="#d97706" fill-opacity="0.3"/>
    <text x="112" y="66" text-anchor="middle" fill="#f59e0b" font-size="8.5" font-weight="700">Shapes</text>
    <rect x="145" y="49" width="55" height="26" rx="3" fill="#9333ea" fill-opacity="0.3"/>
    <text x="172" y="66" text-anchor="middle" fill="#c084fc" font-size="8.5" font-weight="700">Points</text>
    <rect x="205" y="49" width="65" height="26" rx="3" fill="#16a34a" fill-opacity="0.3"/>
    <text x="237" y="66" text-anchor="middle" fill="#4ade80" font-size="8.5" font-weight="700">Tables</text>
    <text x="150" y="94" text-anchor="middle" fill="#cbd5e1" font-size="9">단일 OME-Zarr 청크 저장소 및 지연 평가(Dask)</text>
    <text x="150" y="107" text-anchor="middle" fill="#94a3b8" font-size="8">Transform Graph로 모든 엘리먼트 간 실시간 동기화</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• Images (멀티스케일 피라미드):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">H&amp;E/형광 다채널 거대 이미지 블록 캐싱</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• Shapes &amp; Points (지오메트리 &amp; 분자):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">GeoPandas 폴리곤 및 1천만+ RNA 전사체 점</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• Tables (AnnData 메타데이터 연결):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">발현 행렬 + 세포 주석 + 공간 좌표 상호 링크</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#16a34a" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#4ade80" font-size="10" font-weight="600">결과: 표준화된 오픈소스 멀티모달 공간 데이터</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_cellpose_vs_baysor_segmentation(filepath: str):
    """
    3. Single-Cell Segmentation: Watershed vs Cellpose vs Baysor
    Standard DAPI watershed errors vs Cellpose gradient flow vector fields vs Baysor Bayesian MRF point process
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad03" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">단일 세포 공간 분할 알고리즘: Watershed vs Cellpose vs Baysor</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">단순 형광 임계값 분할의 한계 극복: 딥러닝 경사 흐름(Gradient Flow) vs 전사체 결합 베이지안 점 과정(MRF)</text>

  <!-- Panel 1: DAPI Watershed (Traditional) -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">1. 전통적 DAPI Watershed (오류 유발)</text>

    <!-- Schematic: Clustered nuclei merged or split incorrectly -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    <!-- Merged nuclei -->
    <ellipse cx="60" cy="79" rx="22" ry="15" fill="#f43f5e" fill-opacity="0.3" stroke="#f43f5e" stroke-width="1.5"/>
    <circle cx="53" cy="78" r="4" fill="#38bdf8"/>
    <circle cx="68" cy="78" r="4" fill="#38bdf8"/>
    <text x="60" y="103" text-anchor="middle" fill="#f43f5e" font-size="8" font-weight="700">융합 오류 (Under-seg)</text>
    <!-- Over-split nucleus -->
    <path d="M 125 65 L 145 65 L 140 92 L 120 92 Z" fill="#f43f5e" fill-opacity="0.2" stroke="#f43f5e" stroke-width="1" stroke-dasharray="2,2"/>
    <text x="135" y="103" text-anchor="middle" fill="#f59e0b" font-size="8" font-weight="700">과분할 (Over-seg)</text>
    <!-- Cytoplasmic RNA loss -->
    <circle cx="190" cy="72" r="1.8" fill="#4ade80"/>
    <circle cx="200" cy="85" r="1.8" fill="#4ade80"/>
    <text x="215" y="79" fill="#94a3b8" font-size="8.5">&gt;50% 세포질</text>
    <text x="215" y="90" fill="#f43f5e" font-size="8.5">RNA 누락</text>

    <!-- Technical Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 원리: 거리 변환 및 국소 극소점 침수:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">DAPI 형광 밝기 임계값에 전적으로 의존</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 종양 코어(Dense Core) 취약성:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">핵이 밀집된 암 조직에서 세포 경계 융합</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 세포질 전사체 확장(Expansion) 한계:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">단순 반경 확장은 인접 세포 전사체 오염 초래</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">한계: 부정확한 경계 &amp; 전사체 유실율 높음</text>
  </g>

  <!-- Panel 2: Cellpose (Deep Learning Gradient Flows) -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">2. Cellpose (경사 흐름 벡터장)</text>

    <!-- Schematic: Vector field pointing inward to center -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Irregular cell shape with flow vectors -->
    <path d="M 40 78 Q 50 56 80 60 Q 105 70 95 95 Q 65 102 45 92 Z" fill="#0284c7" fill-opacity="0.2" stroke="#38bdf8" stroke-width="1.5"/>
    <!-- Centroid dot -->
    <circle cx="72" cy="78" r="3" fill="#38bdf8"/>
    <!-- Arrows pointing to centroid -->
    <line x1="52" y1="68" x2="66" y2="74" stroke="#4ade80" stroke-width="1.5"/>
    <line x1="88" y1="68" x2="78" y2="74" stroke="#4ade80" stroke-width="1.5"/>
    <line x1="58" y1="92" x2="68" y2="83" stroke="#4ade80" stroke-width="1.5"/>
    <line x1="86" y1="92" x2="76" y2="82" stroke="#4ade80" stroke-width="1.5"/>
    <text x="120" y="66" fill="#38bdf8" font-size="10" font-weight="700">Vector Field (∇x, ∇y)</text>
    <text x="120" y="80" fill="#cbd5e1" font-size="9.5">확산 시뮬레이션 역추적 (오일러 적분)</text>
    <text x="120" y="94" fill="#4ade80" font-size="9">임의의 불규칙 형태 완벽 복원</text>

    <!-- Technical Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• U-Net 기반 위상 경사 학습:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">세포 중심을 향하는 연속 벡터 흐름 예측</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 세포질 염색(DAPI + Poly-A/IF) 통합:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">핵뿐만 아니라 실제 세포 외곽 막까지 정밀 분할</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 사전 학습 일반화 성능(Zero-shot):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">CAF 방추형 및 암세포의 다형성 형태에 강력</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">장점: 형태학적 유연성 &amp; 높은 이미지 정확도</text>
  </g>

  <!-- Panel 3: Baysor (Bayesian MRF Point Process) -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. Baysor (베이지안 분자 점 과정)</text>

    <!-- Schematic: Transcripts clustered by identity and space -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#9333ea" stroke-width="1"/>
    <!-- Weak DAPI background -->
    <circle cx="55" cy="78" r="14" fill="#0284c7" fill-opacity="0.2" stroke="#38bdf8" stroke-width="1" stroke-dasharray="2,2"/>
    <!-- Cell 1: Red transcripts (Luminal markers) -->
    <circle cx="48" cy="72" r="2.2" fill="#f43f5e"/>
    <circle cx="56" cy="70" r="2.2" fill="#f43f5e"/>
    <circle cx="52" cy="82" r="2.2" fill="#f43f5e"/>
    <circle cx="62" cy="78" r="2.2" fill="#f43f5e"/>
    <!-- Cell 2: Green transcripts (CAF markers) -->
    <circle cx="85" cy="72" r="2.2" fill="#4ade80"/>
    <circle cx="92" cy="80" r="2.2" fill="#4ade80"/>
    <circle cx="88" cy="88" r="2.2" fill="#4ade80"/>
    <!-- Inferred Bayesian boundaries -->
    <ellipse cx="55" cy="76" rx="17" ry="16" fill="none" stroke="#f43f5e" stroke-width="1.5"/>
    <ellipse cx="88" cy="80" rx="14" ry="15" fill="none" stroke="#4ade80" stroke-width="1.5"/>
    <text x="120" y="66" fill="#c084fc" font-size="10" font-weight="700">MRF 결합 확률 모델</text>
    <text x="120" y="80" fill="#cbd5e1" font-size="9.5">P(세포 | 위치좌표 + 유전자 조성)</text>
    <text x="120" y="94" fill="#f59e0b" font-size="9">핵 신호가 약해도 분자 좌표로 보완</text>

    <!-- Technical Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 좌표 + 분자 식별자 동시 추론:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">이웃 분자의 유전자형 마커 일치도를 확률 반영</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• DAPI 염색 불균일/누락 보정:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">분자 좌표·유전자 조성으로 분할 (핵은 선택적 사전정보)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 세포외 배경 노이즈 분자 필터링:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">세포 외 확산 분자를 자동 검출하여 배제</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#9333ea" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">장점: 분자 기반 할당으로 할당율 개선 (수치는 데이터 의존)</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_cell2location_deconvolution_model(filepath: str):
    """
    4. Cell2location Bayesian Spatial Deconvolution Model
    snRNA-seq cell type reference profiles + Visium spot counts -> Hierarchical Bayesian model (negative binomial likelihood) -> cell type abundance map
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad04" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">Cell2location: 계층적 베이지안 포아송 공간 디컨볼루션 모델</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">snRNA-seq 참조 서명 + Visium 스팟 카운트 ➔ 포아송 관측 모델 ➔ 세포형별 절대 풍부도(Absolute Cell Abundance) 사후 추론</text>

  <!-- Panel 1: Reference Signatures & Spatial Input -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 입력: 참조 시그니처 &amp; 공간 카운트</text>

    <!-- Schematic: Matrix reference g_ef and spot counts d_se -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- snRNA Matrix -->
    <rect x="25" y="52" width="45" height="54" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
    <text x="47" y="74" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">snRNA</text>
    <text x="47" y="86" text-anchor="middle" fill="#cbd5e1" font-size="8">g_ef (GxF)</text>
    <!-- Plus symbol -->
    <text x="80" y="82" fill="#cbd5e1" font-size="14" font-weight="700">+</text>
    <!-- Spatial Spot Matrix -->
    <rect x="95" y="52" width="45" height="54" fill="#16a34a" fill-opacity="0.3" stroke="#4ade80" stroke-width="1"/>
    <text x="117" y="74" text-anchor="middle" fill="#4ade80" font-size="8.5" font-weight="700">Visium</text>
    <text x="117" y="86" text-anchor="middle" fill="#cbd5e1" font-size="8">d_se (SxE)</text>
    <text x="155" y="66" fill="#38bdf8" font-size="10" font-weight="700">Track 00 전립선 참조</text>
    <text x="155" y="80" fill="#cbd5e1" font-size="9">8대 세포형 (AR+, t-NEPC,</text>
    <text x="155" y="93" fill="#cbd5e1" font-size="9">DNPC, CAF, CD8 T, 골아세포)</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• snRNA-seq 음이항(NB) 회귀:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">세포형 f별 유전자 e의 기저 발현 프로파일 g_ef</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 다세포 스팟 혼합 모델링:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">단일 스팟에 중첩된 5~15개 세포의 선형 합성</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 배치 효과 및 민감도 차이 보정:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">단일 핵 분리 vs 조직 포획 기술 간 스케일 차</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">입력: 유전자 발현 시그니처 행렬</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 345 200" stroke="#38bdf8" stroke-width="2"/>
  <polygon points="345,196 350,200 345,204" fill="#38bdf8"/>

  <!-- Panel 2: Hierarchical Bayesian Poisson Model -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#f59e0b" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f59e0b" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f59e0b" font-size="13" font-weight="700">2. 계층적 베이지안 생성 모델 (NB)</text>

    <!-- Schematic: Generative formula box -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <text x="150" y="62" text-anchor="middle" fill="#f59e0b" font-size="11" font-weight="700">Y_sg ~ NB(μ_sg, α_g)</text>
    <text x="150" y="80" text-anchor="middle" fill="#4ade80" font-size="9.5" font-family="monospace">μ_sg = (m_g Σ_k w_sk g_kg + s_eg) · y_s</text>
    <text x="150" y="98" text-anchor="middle" fill="#94a3b8" font-size="8.5">w_sk: 세포수 | s_eg: 가산 배경 | y_s: 스팟 검출 민감도</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 음이항(NB) 카운트 우도(Likelihood):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">비정규화 원시 카운트 데이터 직접 모델링</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 비특이적 신호 분리 (가산 배경 s_eg):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">세포 외 유리 mRNA 및 비특이적 결합 오염 격리</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 변분 추론 (Variational Inference):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">ELBO 최대화를 통해 사후 분포 파라미터 최적화</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">핵심: 절대 세포 수 w_sk 추정 (사전 의존)</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 655 200 L 670 200" stroke="#f59e0b" stroke-width="2"/>
  <polygon points="670,196 675,200 670,204" fill="#f59e0b"/>

  <!-- Panel 3: Spatial Abundance Map & q05 Filtering -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#4ade80" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#4ade80" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#4ade80" font-size="13" font-weight="700">3. 절대 세포 풍부도 지도 (Abundance)</text>

    <!-- Schematic: Spatial grid heatmaps -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#16a34a" stroke-width="1"/>
    <!-- Heatmap 1: AR+ Cancer Core -->
    <rect x="30" y="52" width="45" height="45" rx="4" fill="#be123c" fill-opacity="0.7"/>
    <text x="52" y="77" text-anchor="middle" fill="#fff" font-size="8" font-weight="700">AR+ Core</text>
    <!-- Heatmap 2: CAF rim -->
    <rect x="85" y="52" width="45" height="45" rx="4" fill="#d97706" fill-opacity="0.7"/>
    <text x="107" y="77" text-anchor="middle" fill="#fff" font-size="8" font-weight="700">CAF Rim</text>
    <!-- Heatmap 3: TLS T cells -->
    <rect x="140" y="52" width="45" height="45" rx="4" fill="#0284c7" fill-opacity="0.7"/>
    <text x="162" y="77" text-anchor="middle" fill="#fff" font-size="8" font-weight="700">CD8+ TLS</text>
    <text x="215" y="68" fill="#4ade80" font-size="9.5" font-weight="700">q05 분위수</text>
    <text x="215" y="82" fill="#cbd5e1" font-size="8.5">보수적 추정</text>
    <text x="215" y="94" fill="#94a3b8" font-size="8">위양성 95% 억제</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 5% 사후 분위수(q05) 기준 채택:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">신뢰도 95% 이상으로 존재하는 세포 수만 시각화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 전립선암 아형 공간 배타성 규명:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">AR+ 선암 코어 vs t-NEPC/DNPC 침윤 변연부</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 골전이 골아세포 틈새(Niche) 정량화:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">파골세포/골아세포와의 물리적 상호작용 검증</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#16a34a" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#4ade80" font-size="10" font-weight="600">결과: 스팟별 고신뢰도 단일세포 아형 밀도 지도</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_morans_i_spatial_autocorrelation(filepath: str):
    """
    5. Moran's I Spatial Autocorrelation & LISA Scatterplot
    Spatial weight matrix W, spatial lag, Moran's I scatterplot (High-High, Low-Low, High-Low, Low-High)
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad05" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">공간적 자기상관성 &amp; 공간 가변 유전자(SVG): Moran's I &amp; LISA 산점도</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">공간 가중치 행렬 W ➔ 공간 지연(Spatial Lag) 계산 ➔ 글로벌 Moran's I 검정 ➔ 국소 4사분면 공간 핫스팟/이상치 규명</text>

  <!-- Panel 1: Spatial Weight Matrix & Spatial Lag -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 공간 가중치 행렬 W &amp; 지연(Lag)</text>

    <!-- Schematic: k-NN / distance graph and spatial lag vector -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Central node with neighbors -->
    <circle cx="55" cy="79" r="6" fill="#38bdf8"/>
    <circle cx="35" cy="65" r="4.5" fill="#64748b"/>
    <circle cx="75" cy="65" r="4.5" fill="#64748b"/>
    <circle cx="35" cy="93" r="4.5" fill="#64748b"/>
    <circle cx="75" cy="93" r="4.5" fill="#64748b"/>
    <line x1="55" y1="79" x2="35" y2="65" stroke="#38bdf8" stroke-width="1"/>
    <line x1="55" y1="79" x2="75" y2="65" stroke="#38bdf8" stroke-width="1"/>
    <line x1="55" y1="79" x2="35" y2="93" stroke="#38bdf8" stroke-width="1"/>
    <line x1="55" y1="79" x2="75" y2="93" stroke="#38bdf8" stroke-width="1"/>
    <text x="110" y="65" fill="#38bdf8" font-size="10" font-weight="700">k-NN 그래프 (k=6)</text>
    <text x="110" y="79" fill="#cbd5e1" font-size="9">행-표준화 가중치: Σ_j w_ij = 1</text>
    <text x="110" y="94" fill="#4ade80" font-size="9" font-family="monospace">Spatial Lag: [Wz]_i = Σ w_ij z_j</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 이웃 연결 그래프 구성:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">육각형(Hex) 스팟 구조에 기반한 k=6 인접망</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 거리 감쇠 가중치 옵션:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">w_ij = exp(-d_ij² / 2σ²) 가우시안 커널 적용</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 지연(Spatial Lag)의 의미:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">i번 스팟 주변 이웃들의 평균 표준화 발현량</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">토대: 조직 내 국소 이웃 공간 구조 정량화</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 345 200" stroke="#38bdf8" stroke-width="2"/>
  <polygon points="345,196 350,200 345,204" fill="#38bdf8"/>

  <!-- Panel 2: Global Moran's I Formula & SVG Discovery -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#f59e0b" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f59e0b" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f59e0b" font-size="13" font-weight="700">2. 글로벌 Moran's I 통계 수식</text>

    <!-- Schematic: Formula box -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <text x="150" y="64" text-anchor="middle" fill="#f59e0b" font-size="11" font-weight="700">Moran's I = (zᵀ W z) / (zᵀ z)</text>
    <text x="150" y="82" text-anchor="middle" fill="#cbd5e1" font-size="9.5">귀무가설 하 기댓값: E[I] = -1 / (N - 1) ≈ 0</text>
    <text x="150" y="99" text-anchor="middle" fill="#4ade80" font-size="9">I &gt; 0: 공간 응집(SVG) | I &lt; 0: 체스판 분산</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 통계적 유의성 검정 (P-value):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">순열 검정(Permutation) ➔ Benjamini-Hochberg FDR</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 가변 유전자(SVG) 선별:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">KLK3, AR (종양), ACTA2 (CAF)에서 높은 I (예시)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 노이즈 하우스키핑 유전자 배제:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">공간 균일 발현 유전자(I ≈ 0)와 분리 (예시)</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">결과: 조직 특이적 공간 패턴 유전자 SVG 도출</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 655 200 L 670 200" stroke="#f59e0b" stroke-width="2"/>
  <polygon points="670,196 675,200 670,204" fill="#f59e0b"/>

  <!-- Panel 3: LISA 4-Quadrant Scatterplot -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#4ade80" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#4ade80" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#4ade80" font-size="13" font-weight="700">3. 국소 LISA 산점도 (4사분면 해석)</text>

    <!-- Schematic: 4-quadrant plot -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#16a34a" stroke-width="1"/>
    <!-- Crosshairs -->
    <line x1="75" y1="46" x2="75" y2="110" stroke="#475569" stroke-width="1"/>
    <line x1="25" y1="78" x2="125" y2="78" stroke="#475569" stroke-width="1"/>
    <!-- Quadrant dots -->
    <circle cx="100" cy="58" r="3.5" fill="#f43f5e"/><text x="105" y="56" fill="#f43f5e" font-size="7" font-weight="700">Q1: HH</text>
    <circle cx="50" cy="58" r="3" fill="#38bdf8"/><text x="32" y="56" fill="#38bdf8" font-size="7" font-weight="700">Q2: LH</text>
    <circle cx="50" cy="98" r="3.5" fill="#64748b"/><text x="32" y="104" fill="#64748b" font-size="7" font-weight="700">Q3: LL</text>
    <circle cx="100" cy="98" r="3" fill="#f59e0b"/><text x="105" y="104" fill="#f59e0b" font-size="7" font-weight="700">Q4: HL</text>
    <text x="140" y="64" fill="#f43f5e" font-size="9" font-weight="700">HH: 암 코어 핫스팟</text>
    <text x="140" y="78" fill="#64748b" font-size="9" font-weight="700">LL: 비발현 콜드스팟</text>
    <text x="140" y="92" fill="#f59e0b" font-size="9" font-weight="700">HL/LH: 경계부 이상치</text>
    <text x="140" y="106" fill="#94a3b8" font-size="8">침윤 마이크로 전이 검출</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• Q1 (High-High) 공간 핫스팟:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">자신과 주변 이웃 모두 고발현 ➔ 단단한 암 중심부</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• Q3 (Low-Low) 공간 콜드스팟:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">발현이 결손된 기질 및 정상 선조직 영역</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• Q2/Q4 공간 이상치 (Spatial Outliers):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">침윤 경계(Invasive Edge) 및 고립된 암세포 규명</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#16a34a" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#4ade80" font-size="10" font-weight="600">활용: 조직학적 경계와 일치하는 미세 구획 자동 분할</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_spatial_niche_and_commot_signaling(filepath: str):
    """
    6. Cellular Neighborhood Niche & Distance-Decay Receptor-Ligand Signaling
    Cellular neighborhood windowing, spatial niche clusters, and distance-decay ligand-receptor gradient
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad06" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">세포 틈새(Spatial Niche) 분석 &amp; 거리 감쇠 리간드-수용체 상호작용</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">국소 세포 이웃(CN) 윈도우 집계 ➔ 4대 공간 틈새 클러스터링 ➔ COMMOT 최적 수송 거리 제약 신호 흐름</text>

  <!-- Panel 1: Cellular Neighborhood (CN) Windowing -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 세포 이웃(CN) 윈도우 집계</text>

    <!-- Schematic: Center cell with radius r=30um circle counting neighbors -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Radius circle -->
    <circle cx="65" cy="79" r="26" fill="#0284c7" fill-opacity="0.2" stroke="#38bdf8" stroke-width="1.5" stroke-dasharray="3,2"/>
    <!-- Index cell -->
    <circle cx="65" cy="79" r="5" fill="#f43f5e"/>
    <!-- Diverse neighbor cells inside radius -->
    <circle cx="55" cy="68" r="4" fill="#f43f5e"/>
    <circle cx="75" cy="66" r="4" fill="#f59e0b"/>
    <circle cx="56" cy="92" r="4" fill="#4ade80"/>
    <circle cx="78" cy="88" r="3.5" fill="#c084fc"/>
    <text x="105" y="65" fill="#38bdf8" font-size="10" font-weight="700">국소 반경 윈도우 (r, 예: 30~50 µm)</text>
    <text x="105" y="79" fill="#cbd5e1" font-size="9">주변 k개 최근접 세포 탐색</text>
    <text x="105" y="94" fill="#4ade80" font-size="8.5">세포형 조성 벡터: v_i = [p₁, p₂, .., pₖ]</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 세포 기반 이웃 반경 정의:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">세포 간 직접 접촉 및 측분비(Paracrine) 확산 한계</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 국소 세포형 빈도 벡터 변환:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">각 단일 세포를 주변 환경 구성비(Composition)로 인코딩</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 스케일 불변 정규화:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">총 세포 수 편차 보정 (Σ p_k = 1.0)</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">토대: 세포 단위 환경을 고차원 프로파일로 사상</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 345 200" stroke="#38bdf8" stroke-width="2"/>
  <polygon points="345,196 350,200 345,204" fill="#38bdf8"/>

  <!-- Panel 2: Spatial Niche Clustering -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f59e0b" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f59e0b" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f59e0b" font-size="13" font-weight="700">2. 4대 공간 틈새(Niche) 클러스터</text>

    <!-- Schematic: 4 niche badge chips -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <!-- Niche 1 -->
    <rect x="22" y="48" width="125" height="28" rx="4" fill="#be123c" fill-opacity="0.3"/>
    <text x="27" y="61" fill="#f43f5e" font-size="8.5" font-weight="700">Niche 1: Cancer Core</text>
    <text x="27" y="71" fill="#cbd5e1" font-size="7.5">AR+ 선암 상피세포 농축</text>
    <!-- Niche 2 -->
    <rect x="153" y="48" width="125" height="28" rx="4" fill="#d97706" fill-opacity="0.3"/>
    <text x="158" y="61" fill="#f59e0b" font-size="8.5" font-weight="700">Niche 2: Desmoplastic Rim</text>
    <text x="158" y="71" fill="#cbd5e1" font-size="7.5">CAF / 기질 경계부 장벽</text>
    <!-- Niche 3 -->
    <rect x="22" y="81" width="125" height="28" rx="4" fill="#0284c7" fill-opacity="0.3"/>
    <text x="27" y="94" fill="#38bdf8" font-size="8.5" font-weight="700">Niche 3: TLS Immune Hub</text>
    <text x="27" y="104" fill="#cbd5e1" font-size="7.5">CD8 T / B / DC 면역 림프</text>
    <!-- Niche 4 -->
    <rect x="153" y="81" width="125" height="28" rx="4" fill="#16a34a" fill-opacity="0.3"/>
    <text x="158" y="94" fill="#4ade80" font-size="8.5" font-weight="700">Niche 4: Bone Front</text>
    <text x="158" y="104" fill="#cbd5e1" font-size="7.5">골아/파골세포 골기질 재형성</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• K-Means / Leiden 비지도 클러스터링:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">조성 벡터 군집화로 반복되는 미세조직 구조 발견</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 종양 침윤 변연부(Invasive Border):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">섬유아세포 기질 장벽과 암세포 간의 물리적 대치</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 골전이 특이적 전선(Front):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">전립선암의 조골성(Osteoblastic) 병변 미세환경 형성</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">결과: 생물학적으로 해석 가능한 4대 공간 틈새</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 655 200 L 670 200" stroke="#f59e0b" stroke-width="2"/>
  <polygon points="670,196 675,200 670,204" fill="#f59e0b"/>

  <!-- Panel 3: Distance-Decayed Ligand-Receptor Signaling (COMMOT) -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#4ade80" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#4ade80" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#4ade80" font-size="13" font-weight="700">3. 거리 감쇠 L-R 신호망 (COMMOT)</text>

    <!-- Schematic: Sender-Receiver signaling gradient -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#16a34a" stroke-width="1"/>
    <!-- Sender cell -->
    <circle cx="45" cy="79" r="10" fill="#d97706" fill-opacity="0.6"/>
    <text x="45" y="82" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">CAF</text>
    <!-- Diffusion gradient curve -->
    <path d="M 60 79 Q 95 65 130 79" fill="none" stroke="#f59e0b" stroke-width="2"/>
    <polygon points="130,76 135,79 130,82" fill="#f59e0b"/>
    <!-- Receiver cell -->
    <circle cx="150" cy="79" r="10" fill="#be123c" fill-opacity="0.6"/>
    <text x="150" y="82" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">Cancer</text>
    <text x="180" y="65" fill="#4ade80" font-size="8.5" font-weight="700">S_ij = L_i·R_j·e^(-d/d_0)</text>
    <text x="180" y="78" fill="#f59e0b" font-size="8">TGFB1 ➔ TGFBR2 (EMT)</text>
    <text x="180" y="91" fill="#c084fc" font-size="8">RANKL ➔ RANK (골흡수)</text>
    <text x="180" y="104" fill="#94a3b8" font-size="7.5">거리 초과 시 지수적 소멸</text>

    <!-- Technical Details -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 물리적 확산 제약 (Distance Decay):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단순 동시 발현 오류 방지 (d_0: 신호 유형별 설정)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 최적 수송(Optimal Transport) 플럭스:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">공급자(리간드)와 소비자(수용체) 간 방향성 흐름 추정</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 전립선암 골전이 악성화 축 규명:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">CAF 유래 TGFβ 신호 및 암세포-파골세포 RANKL 회로</text>

    <rect x="15" y="242" width="270" height="22" rx="4" fill="#16a34a" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#4ade80" font-size="10" font-weight="600">결과: 실제 공간 거리 기반 정밀 세포 간 통신망 규명</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_spatial_technologies_resolution_comparison.svg", make_01_spatial_technologies_resolution),
        ("02_coordinate_registration_and_spatialdata.svg", make_02_coordinate_registration_spatialdata),
        ("03_cellpose_vs_baysor_segmentation.svg", make_03_cellpose_vs_baysor_segmentation),
        ("04_cell2location_deconvolution_model.svg", make_04_cell2location_deconvolution_model),
        ("05_morans_i_and_spatial_autocorrelation.svg", make_05_morans_i_spatial_autocorrelation),
        ("06_spatial_niche_and_commot_signaling.svg", make_06_spatial_niche_and_commot_signaling),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 02 Spatial Transcriptomics publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

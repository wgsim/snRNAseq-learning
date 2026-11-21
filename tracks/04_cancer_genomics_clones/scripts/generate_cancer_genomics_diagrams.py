#!/usr/bin/env python3
"""
tracks/04_cancer_genomics_clones/scripts/generate_cancer_genomics_diagrams.py
=============================================================================
Generates publication-grade, high-resolution SVG diagrams for Track 04:
Cancer Genomics & Clonal Evolution in Prostate Cancer.

Covers:
1. Prostate Cancer Hallmark CNA Landscape & Genomic Drivers
2. inferCNV: Genomic Coordinate Smoothing & 6-State HMM Workflow
3. CopyKAT: Poisson Variance Stabilization & GMM Diploid/Aneuploid Separation
4. Numbat: Haplotype Phased BAF & Copy-Neutral LOH (CN-LOH) Detection
5. Clonal Phylogeny Tree (Neighbor-Joining) & Muller Plot Dynamics
6. Therapy Selection Bottleneck & Metastatic Seeding Patterns
"""

import os
import xml.etree.ElementTree as ET


def make_01_prostate_cancer_hallmark_cna_landscape(filepath: str):
    """
    1. Prostate Cancer Hallmark CNA Landscape & Genomic Drivers
    Chromosome Ideograms (8p/q, 10q, 13q, 17p, 21q, Xq) & Key Drivers (NKX3-1, MYC, PTEN, RB1, TP53, TMPRSS2-ERG, AR).
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrowCyan01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrowRose01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">전립선암 특징적 염색체 수 이상(CNA) 파노라마 &amp; 핵심 종양 유전자</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">염색체 이데오그램(8p/q, 10q, 13q, 17p, 21q, Xq)과 핵심 드라이버 (PTEN, TP53, RB1, MYC, TMPRSS2-ERG, AR)</text>

  <!-- Panel 1: Deletion Hotspots -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">1. 주요 결손 핫스팟 (Deletion Hotspots)</text>

    <!-- Schematic: Chromosome Ideograms with Deletions -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    
    <!-- Chr 8p -->
    <g transform="translate(25, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="1" y="2" width="10" height="24" rx="3" fill="#f43f5e" fill-opacity="0.9"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 8</text>
      <text x="18" y="16" fill="#f43f5e" font-size="8.5" font-weight="700">8p- (NKX3-1)</text>
    </g>

    <!-- Chr 10q -->
    <g transform="translate(95, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="1" y="30" width="10" height="24" rx="3" fill="#f43f5e" fill-opacity="0.9"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 10</text>
      <text x="18" y="44" fill="#f43f5e" font-size="8.5" font-weight="700">10q- (PTEN)</text>
    </g>

    <!-- Chr 13q -->
    <g transform="translate(165, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="1" y="32" width="10" height="22" rx="3" fill="#fbbf24" fill-opacity="0.9"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 13</text>
      <text x="18" y="46" fill="#fbbf24" font-size="8.5" font-weight="700">13q- (RB1)</text>
    </g>

    <!-- Chr 17p -->
    <g transform="translate(235, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="1" y="2" width="10" height="24" rx="3" fill="#e11d48" fill-opacity="0.9"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 17</text>
      <text x="-4" y="16" text-anchor="end" fill="#e11d48" font-size="8.5" font-weight="700">17p- (TP53)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 8p &amp; 10q 동시 결손 (초기 클론 드라이버):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">NKX3-1(상피분화 소실) &amp; PTEN(PI3K-AKT 활성화)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 13q &amp; 17p 결손 (후기 악성화 드라이버):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">RB1(세포주기 제동 소실) &amp; TP53(게놈 불안정성)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 2-Hit 불활성화 모델 (Knudson):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">대립유전자 결손(LOH) + 점돌연변이 또는 후성 억제</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#f43f5e" fill-opacity="0.2"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">결손 빈도: PTEN 40~60%, TP53 30~50%, RB1 25~35%</text>
  </g>

  <!-- Panel 2: Gains & Fusions -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">2. 증폭 &amp; 염색체 융합 (Gains &amp; Fusions)</text>

    <!-- Schematic: Gains & Amplicons -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <!-- Chr 8q Gain -->
    <g transform="translate(30, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="1" y="30" width="10" height="24" rx="3" fill="#38bdf8" fill-opacity="0.9"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 8</text>
      <text x="18" y="45" fill="#38bdf8" font-size="8.5" font-weight="700">8q+ (MYC)</text>
    </g>

    <!-- Chr 21q Deletion / Fusion -->
    <g transform="translate(115, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <line x1="6" y1="26" x2="6" y2="38" stroke="#f43f5e" stroke-width="3" stroke-dasharray="2,2"/>
      <circle cx="6" cy="32" r="5" fill="#c084fc"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 21</text>
      <text x="18" y="35" fill="#c084fc" font-size="8" font-weight="700">TMPRSS2-ERG</text>
    </g>

    <!-- Chr Xq High Amplicon (AR) -->
    <g transform="translate(210, 50)">
      <rect x="0" y="0" width="12" height="56" rx="6" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <rect x="0" y="16" width="12" height="20" rx="3" fill="#fbbf24" stroke="#f59e0b" stroke-width="1.5"/>
      <text x="6" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">Chr X</text>
      <text x="18" y="29" fill="#fbbf24" font-size="8.5" font-weight="700">Xq12 (AR Amp)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 8q24 증폭 (MYC 과발현):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단백질 합성 및 대사 재프로그래밍으로 성장 촉진</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 21q22.2 결손-융합 (TMPRSS2-ERG):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">안드로겐 유도 프로모터에 의한 발암성 전사인자 폭발 (~50%)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• Xq12 유전자 좌위 증폭 (AR Super-Amp):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">거세 저항성 상태에서 AR 및 인핸서 좌위 10~50배 증폭</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">진행 기전: 호르몬 결핍 후 AR 증폭 및 변이 클론 팽창</text>
  </g>

  <!-- Panel 3: Clinical Trajectory -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 임상 진행 단계별 이수체 동역학</text>

    <!-- Schematic: Trajectory Arrow -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <!-- Stages Boxes -->
    <g transform="translate(25, 52)">
      <!-- Primary -->
      <rect x="0" y="0" width="68" height="42" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="34" y="16" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">Primary PCa</text>
      <text x="34" y="30" text-anchor="middle" fill="#94a3b8" font-size="7.5">국소 결손 / ERG</text>
      <text x="34" y="38" text-anchor="middle" fill="#64748b" font-size="7">CNA &lt; 5%</text>

      <line x1="72" y1="21" x2="88" y2="21" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrowCyan01)"/>

      <!-- mCRPC -->
      <rect x="92" y="0" width="68" height="42" rx="4" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="126" y="16" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">mCRPC</text>
      <text x="126" y="30" text-anchor="middle" fill="#94a3b8" font-size="7.5">AR Amp / PTEN-</text>
      <text x="126" y="38" text-anchor="middle" fill="#64748b" font-size="7">CNA 25~40%</text>

      <line x1="164" y1="21" x2="180" y2="21" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrowRose01)"/>

      <!-- t-NEPC -->
      <rect x="184" y="0" width="66" height="42" rx="4" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
      <text x="217" y="16" text-anchor="middle" fill="#f43f5e" font-size="8.5" font-weight="700">t-NEPC</text>
      <text x="217" y="30" text-anchor="middle" fill="#94a3b8" font-size="7.5">TP53/RB1 소실</text>
      <text x="217" y="38" text-anchor="middle" fill="#64748b" font-size="7">CNA &gt; 55%</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 원발암 (Primary Adenocarcinoma):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">비교적 정적인 유전체, 국소적 8p/PTEN 결손 및 ERG 융합</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 거세저항성 전이암 (mCRPC):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">광범위 염색체 완(Arm) 단위 결손/증폭, AR 고도 증폭</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 치료 유도 신경내분비암 (t-NEPC):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">TP53 &amp; RB1 동시 소실로 인한 극단적 유전체 불안정성</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">임상 지표: CNA 부하(CNA Burden)는 전이 및 사망 위험 예측</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_infercnv_hmm_smoothing_workflow(filepath: str):
    """
    2. inferCNV: Genomic Coordinate Smoothing & 6-State HMM Workflow
    Genomic Coordinate Ordering -> Normal Baseline Subtraction -> 100-Gene Moving Average -> 6-State HMM Viterbi Segmentation.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad02" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrowCyan02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">inferCNV: 유전체 좌표 기반 평활화 &amp; 6-상태 HMM 분절화 워크플로우</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">염색체 위치 정렬 ➔ 2배체 기저선 차감 ➔ 100-유전자 이동평균 ➔ 6-상태 은닉 마르코프 모델(HMM) 비터비 분절</text>

  <!-- Panel 1: Ordering & Baseline Subtraction -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 좌표 정렬 &amp; 정상 기저선 차감</text>

    <!-- Schematic: Chromosome Tracks & Subtraction -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <!-- Chromosome blocks -->
    <g transform="translate(25, 52)">
      <rect x="0" y="0" width="40" height="18" rx="2" fill="#334155"/>
      <text x="20" y="12" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 1</text>
      <rect x="44" y="0" width="36" height="18" rx="2" fill="#1e293b"/>
      <text x="62" y="12" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 2</text>
      <text x="90" y="12" fill="#64748b" font-size="10">···</text>
      <rect x="105" y="0" width="32" height="18" rx="2" fill="#334155"/>
      <text x="121" y="12" text-anchor="middle" fill="#94a3b8" font-size="8">Chr 8</text>
      <text x="146" y="12" fill="#64748b" font-size="10">···</text>
      <rect x="160" y="0" width="28" height="18" rx="2" fill="#334155"/>
      <text x="174" y="12" text-anchor="middle" fill="#94a3b8" font-size="8">Chr X</text>
    </g>

    <!-- Formula Subtraction Box -->
    <g transform="translate(25, 78)">
      <rect x="0" y="0" width="250" height="30" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="125" y="19" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-family="monospace">log2(E_tumor + 1) - log2(E_normal_ref + 1)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 유전체 물리적 좌표 정렬:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">GENCODE 주석 기준 전사 시작점(TSS) 순 정렬</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 참조 2배체 기저선 (Normal Ref):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">T/B 림프구, 내피세포, 정상 기질 세포 풀링 기준선</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 유전자별 평균 발현 표준화:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">시퀀싱 깊이(Depth) 및 라이브러리 크기 인자 보정</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">출력: 정상 대조군 대비 상대적 log2 발현 배율 행렬</text>
  </g>

  <!-- Panel 2: Moving Average Smoothing -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 100-유전자 이동평균 평활화</text>

    <!-- Schematic: Raw Noisy vs Smoothed Signal -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Baseline zero -->
      <line x1="0" y1="32" x2="250" y2="32" stroke="#475569" stroke-width="1" stroke-dasharray="2,2"/>
      
      <!-- Raw noisy signal (grey) -->
      <polyline points="5,32 15,10 25,48 35,5 45,40 55,18 65,52 75,22 85,38 95,12 105,45 115,20 125,10 135,15 145,12 155,50 165,55 175,48 185,58 195,50 205,32 215,30 225,34 235,32"
        fill="none" stroke="#64748b" stroke-width="1" stroke-opacity="0.6"/>

      <!-- 100-gene moving window highlight box -->
      <rect x="100" y="5" width="60" height="52" fill="#34d399" fill-opacity="0.15" stroke="#34d399" stroke-width="1" stroke-dasharray="2,2"/>
      <text x="130" y="62" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">100-gene Window</text>

      <!-- Smoothed signal (cyan line) -->
      <polyline points="5,32 50,30 90,20 120,15 145,15 160,52 195,52 215,32 245,32"
        fill="none" stroke="#38bdf8" stroke-width="2.5"/>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 단일세포 전사체 드롭아웃 극복:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">개별 유전자의 고빈도 영점(Zero-count) 노이즈 여과</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 유전체 연관 신호대잡음비(SNR) 극대화:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">인접 유전자 집합의 복제수 공유 물리학 반영</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 동적 염색체 경계 처리 (Centromere Clamping):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">염색체 완(Arm) 간 신호 오염 방지 경계 절단</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">윈도우 규격: 표준 101 유전자 (Default span: 100 genes)</text>
  </g>

  <!-- Panel 3: 6-State HMM & Viterbi -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. 6-상태 HMM &amp; 비터비 분절화</text>

    <!-- Schematic: 6 discrete states -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>

    <!-- State Badges -->
    <g transform="translate(25, 50)">
      <!-- State 1: 0x Loss -->
      <rect x="0" y="0" width="36" height="24" rx="3" fill="#be123c"/>
      <text x="18" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">0x Del</text>
      <text x="18" y="21" text-anchor="middle" fill="#fecdd3" font-size="6.5">-1.5</text>

      <!-- State 2: 1x Loss -->
      <rect x="42" y="0" width="36" height="24" rx="3" fill="#f43f5e"/>
      <text x="60" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">1x Loss</text>
      <text x="60" y="21" text-anchor="middle" fill="#fecdd3" font-size="6.5">-0.6</text>

      <!-- State 3: 2x Neutral -->
      <rect x="84" y="0" width="36" height="24" rx="3" fill="#475569"/>
      <text x="102" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">2x Neut</text>
      <text x="102" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6.5">0.0</text>

      <!-- State 4: 3x Gain -->
      <rect x="126" y="0" width="36" height="24" rx="3" fill="#0284c7"/>
      <text x="144" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">3x Gain</text>
      <text x="144" y="21" text-anchor="middle" fill="#bae6fd" font-size="6.5">+0.4</text>

      <!-- State 5: 4x Gain -->
      <rect x="168" y="0" width="36" height="24" rx="3" fill="#3b82f6"/>
      <text x="186" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">4x Gain</text>
      <text x="186" y="21" text-anchor="middle" fill="#bfdbfe" font-size="6.5">+0.8</text>

      <!-- State 6: 5x Amp -->
      <rect x="210" y="0" width="40" height="24" rx="3" fill="#d97706"/>
      <text x="230" y="11" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">5x+ Amp</text>
      <text x="230" y="21" text-anchor="middle" fill="#fef3c7" font-size="6.5">+1.3</text>

      <!-- Viterbi Decoding Track -->
      <g transform="translate(0, 32)">
        <text x="0" y="18" fill="#94a3b8" font-size="8">Viterbi 경로:</text>
        <rect x="52" y="8" width="50" height="12" fill="#475569" rx="2"/>
        <text x="77" y="17" text-anchor="middle" fill="#fff" font-size="7">2x Neutral</text>

        <rect x="104" y="8" width="60" height="12" fill="#f43f5e" rx="2"/>
        <text x="134" y="17" text-anchor="middle" fill="#fff" font-size="7">1x Loss (10q)</text>

        <rect x="166" y="8" width="84" height="12" fill="#0284c7" rx="2"/>
        <text x="208" y="17" text-anchor="middle" fill="#fff" font-size="7">3x Gain (8q MYC)</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 6단계 유전체 복제수 이산화 (Discretization):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">0x(호모결손)부터 5x+(초고도 증폭)까지 모델링</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 전이 확률 행렬 (Self-Transition ~0.98):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">인접 유전자 간 동일 복제 상태 유지 관성 부여</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 비터비(Viterbi) 동적 계획법 최적화:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">전체 전사체에서 최대 사후 우도 경로를 역추적 탐색</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#7e22ce" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">최종 결과: 단일세포별 정밀 CNA 블록 분절화 지도 완성</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_copykat_gmm_diploid_aneuploid_separation(filepath: str):
    """
    3. CopyKAT: Poisson Variance Stabilization & GMM Diploid/Aneuploid Separation
    Reference-Free Poisson Variance Stabilization -> Bimodal Variance GMM -> Aneuploid Tumor vs Diploid Stroma Identification.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad03" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">CopyKAT: 포아송 분산 안정화 &amp; GMM 기반 2배체/이수체 비지도 분리</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">기준세포 비의존 분산 안정화 ➔ 이봉 게놈 분산 GMM 군집 ➔ 2배체 정상 기질 vs 이수체 종양 세포 분리</text>

  <!-- Panel 1: Reference-Free Poisson Stabilization -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 기준세포 비의존 분산 안정화</text>

    <!-- Schematic: Poisson count to Freeman-Tukey / Anscombe -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 52)">
      <rect x="0" y="0" width="105" height="48" rx="4" fill="#1e293b" stroke="#475569" stroke-width="0.8"/>
      <text x="52" y="16" text-anchor="middle" fill="#f43f5e" font-size="8" font-weight="700">Raw scRNA Count</text>
      <text x="52" y="30" text-anchor="middle" fill="#94a3b8" font-size="7">Poisson Noise: y ~ Poiss(λ)</text>
      <text x="52" y="42" text-anchor="middle" fill="#64748b" font-size="6.5">Var(y) = Mean(y)</text>

      <text x="118" y="28" fill="#38bdf8" font-size="14" font-weight="700">→</text>

      <rect x="135" y="0" width="110" height="48" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="190" y="16" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">Anscombe Transform</text>
      <text x="190" y="30" text-anchor="middle" fill="#34d399" font-size="7.5" font-family="monospace">z = 2√(y + 3/8)</text>
      <text x="190" y="42" text-anchor="middle" fill="#94a3b8" font-size="6.5">Stabilized Var(z) ≈ 1.0</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 정상 대조 세포 사전 지식 불필요:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일 종양 검체 내 자체 데이터만으로 복제수 프로파일링</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 포아송 분산 안정화 (Anscombe 변환):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">발현량 크기에 종속적인 기저 분산을 상수(1.0)로 정규화</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 5MB 게놈 윈도우 집계:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">국소 염색체 단위 윈도우 스팬으로 전사율 변동 흡수</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">특징: 종양 순도(Purity)가 낮은 생검 검체 분석 최적화</text>
  </g>

  <!-- Panel 2: GMM Diploid Discovery -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. GMM 이봉 분산 2배체 발견</text>

    <!-- Schematic: Bimodal distribution & GMM -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <line x1="10" y1="52" x2="245" y2="52" stroke="#475569" stroke-width="1"/>
      <line x1="10" y1="10" x2="10" y2="52" stroke="#475569" stroke-width="1"/>

      <!-- Diploid peak curve (emerald) -->
      <path d="M 15,52 Q 45,5 75,52" fill="#34d399" fill-opacity="0.25" stroke="#34d399" stroke-width="2"/>
      <text x="45" y="25" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">2배체 (Diploid)</text>
      <text x="45" y="34" text-anchor="middle" fill="#94a3b8" font-size="6.5">낮은 분산 σ²~0.02</text>

      <!-- Aneuploid peak curve (rose) -->
      <path d="M 85,52 Q 165,18 240,52" fill="#f43f5e" fill-opacity="0.25" stroke="#f43f5e" stroke-width="2"/>
      <text x="165" y="30" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">이수체 (Aneuploid)</text>
      <text x="165" y="40" text-anchor="middle" fill="#94a3b8" font-size="6.5">높은 분산 σ²~0.25</text>

      <text x="210" y="60" text-anchor="middle" fill="#64748b" font-size="7">게놈 분산 (σ²)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 게놈 복제수 분산 지표 (Genomic Variance):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포별 염색체 전반의 분절 발현 변동폭(Variance) 산출</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 이봉 분포(Bimodal) 자동 분할:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">EM 알고리즘 기반 2-Component 가우시안 혼합 모델 피팅</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 내부 2배체 기준군 자동 추출:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">최저 분산 성분 세포군을 고신뢰 정상 기준군(Baseline) 지정</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">분류 기준: 사후 확률 P(Aneuploid | σ²) &gt; 0.5 판정</text>
  </g>

  <!-- Panel 3: Separation & Validation -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#fbbf24" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#fbbf24" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#fbbf24" font-size="13" font-weight="700">3. 종양 vs 기질 분리 &amp; 벤치마크</text>

    <!-- Schematic: Clustered Heatmap & ROC Curve -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Heatmap rows -->
      <rect x="0" y="2" width="135" height="18" fill="#334155" rx="2"/>
      <text x="67" y="14" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">Diploid Normal (Stroma/Immune)</text>

      <rect x="0" y="23" width="135" height="15" fill="#f43f5e" fill-opacity="0.85" rx="2"/>
      <text x="67" y="34" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">Tumor Clone A (8p-, 8q+, 10q-)</text>

      <rect x="0" y="41" width="135" height="15" fill="#be123c" fill-opacity="0.9" rx="2"/>
      <text x="67" y="52" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">Tumor Clone B (+13q-, 17p-)</text>

      <!-- Inset ROC Curve -->
      <g transform="translate(155, 2)">
        <rect x="0" y="0" width="90" height="54" rx="3" fill="#1e293b" stroke="#475569" stroke-width="0.8"/>
        <line x1="8" y1="46" x2="82" y2="46" stroke="#64748b" stroke-width="1"/>
        <line x1="8" y1="8" x2="8" y2="46" stroke="#64748b" stroke-width="1"/>
        <line x1="8" y1="46" x2="82" y2="8" stroke="#475569" stroke-width="1" stroke-dasharray="2,2"/>
        <path d="M 8,46 L 10,12 L 82,8" fill="none" stroke="#fbbf24" stroke-width="2"/>
        <text x="45" y="30" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">ROC AUC</text>
        <text x="45" y="40" text-anchor="middle" fill="#fff" font-size="8.5" font-weight="700">0.992</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 계층적 클러스터링(Hierarchical Clustering):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포 복제수 프로파일 거리 기반 서브클론 계통 분기</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 정상 비종양 세포 완벽 분리:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">평평한 2배체 기저선을 갖는 섬유아세포/T세포 정확한 격리</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• DNA Flow Cytometry 검증 일치도:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">WGS 기반 DNA Ploidy 지표와 98.4% 일치하는 정확도</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#fbbf24" font-size="10" font-weight="600">성능 지표: 민감도 99.1%, 특이도 98.6%, AUC 0.992</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_numbat_haplotype_phased_baf_and_cnloh(filepath: str):
    """
    4. Numbat: Haplotype Phased BAF & Copy-Neutral LOH (CN-LOH) Detection
    Total Expression (Poisson) vs Phased BAF (Binomial) -> Population Phasing & Switch Error Mitigation -> Copy-Neutral LOH Detection.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad04" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">Numbat: 대립유전자 BAF &amp; 반수체 위상 기반 복제중립 LOH (CN-LOH) 탐지</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">총발현량(Poisson)과 위상 BAF(Binomial) 결합 ➔ 집단 연관 위상화 및 전이 에러 보정 ➔ 10q PTEN 좌위 CN-LOH 규명</text>

  <!-- Panel 1: Dual Evidence Integration -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 총발현량과 BAF의 이중 증거</text>

    <!-- Schematic: Expression Depth vs BAF at Heterozygous SNPs -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 50)">
      <!-- Expression Depth Track -->
      <rect x="0" y="0" width="115" height="24" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="57" y="11" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">총 발현량 (Total Depth)</text>
      <text x="57" y="20" text-anchor="middle" fill="#94a3b8" font-size="6.5">Poisson / Neg-Binomial: Y ~ NB(μ, θ)</text>

      <!-- Plus sign -->
      <text x="125" y="28" fill="#cbd5e1" font-size="14" font-weight="700">+</text>

      <!-- Phased BAF Track -->
      <rect x="135" y="0" width="115" height="24" rx="3" fill="#1e293b" stroke="#f43f5e" stroke-width="0.8"/>
      <text x="192" y="11" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">위상 BAF (Allelic BAF)</text>
      <text x="192" y="20" text-anchor="middle" fill="#94a3b8" font-size="6.5">Binomial / Beta-Binomial: K_b ~ Bin(N, p)</text>

      <!-- Joint Likelihood Box -->
      <rect x="0" y="32" width="250" height="24" rx="4" fill="#0284c7" fill-opacity="0.2" stroke="#38bdf8" stroke-width="1"/>
      <text x="125" y="48" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">결합 우도: log L_total = log L_depth + log L_baf</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 총발현량의 한계 (Depth Blindness):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포 드롭아웃 및 증폭 잡음으로 복제수 미세 변동 간과</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 이형접합(Heterozygous) SNP 대립유전자 비대칭:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">B-대립유전자 빈도(BAF)의 0.5 이탈을 통한 강력한 증거</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 상호 보완적 결합 확률 모델:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">발현량이 평평해도 BAF 편향을 통해 LOH 완벽 탐지</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">원리: 카피 수 변화와 대립유전자 불균형의 동시 확률 적분</text>
  </g>

  <!-- Panel 2: Haplotype Phasing & Switch Error -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 집단 위상화 &amp; 스위치 에러 보정</text>

    <!-- Schematic: Haplotype A/B & Switch Error HMM -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 50)">
      <!-- Unphased vs Phased -->
      <g transform="translate(0, 0)">
        <text x="0" y="10" fill="#94a3b8" font-size="7.5">비위상(Unphased):</text>
        <circle cx="75" cy="8" r="3" fill="#f43f5e"/>
        <circle cx="95" cy="8" r="3" fill="#38bdf8"/>
        <circle cx="115" cy="8" r="3" fill="#f43f5e"/>
        <circle cx="135" cy="8" r="3" fill="#38bdf8"/>
        <text x="175" y="11" fill="#64748b" font-size="7">산발적 노이즈</text>
      </g>

      <!-- Phased block via 1000G -->
      <g transform="translate(0, 20)">
        <text x="0" y="10" fill="#34d399" font-size="7.5" font-weight="700">1000G 위상화:</text>
        <rect x="65" y="2" width="80" height="12" fill="#34d399" fill-opacity="0.2" rx="2" stroke="#34d399" stroke-width="0.8"/>
        <circle cx="75" cy="8" r="3" fill="#38bdf8"/>
        <circle cx="95" cy="8" r="3" fill="#38bdf8"/>
        <circle cx="115" cy="8" r="3" fill="#38bdf8"/>
        <circle cx="135" cy="8" r="3" fill="#38bdf8"/>
        <text x="160" y="11" fill="#38bdf8" font-size="7" font-weight="700">Haplotype 1 집계</text>
      </g>

      <!-- Switch error mitigation -->
      <g transform="translate(0, 38)">
        <text x="0" y="12" fill="#fbbf24" font-size="7.5">스위치 에러(γ):</text>
        <text x="75" y="12" fill="#94a3b8" font-size="7">HMM 전이 확률 γ ~ 10⁻⁴ 로 위상 역전 보정</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 1000 Genomes 외부 연관불평형(LD) 패널:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포의 희소한 SNP를 모계/부계 반수체로 위상화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 위상 스위칭 에러(Switch Error) 허용:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">위상 반전 확률(γ)을 내장하여 인위적 대립 왜곡 차단</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 대립유전자 신호 집계력(Power) 10배 증대:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">수백 개 SNP를 하나의 메이저 반수체 블록으로 통합 검정</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">위상화 효과: 희소한 단일세포 전사체 BAF 통계적 검정력 달성</text>
  </g>

  <!-- Panel 3: Copy-Neutral LOH Detection -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">3. 복제중립 LOH (CN-LOH) 검출</text>

    <!-- Schematic: 10q PTEN locus CN-LOH -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Expression vs BAF comparison on Chr 10q -->
      <g transform="translate(0, 0)">
        <rect x="0" y="0" width="245" height="24" rx="3" fill="#1e293b"/>
        <text x="10" y="15" fill="#94a3b8" font-size="7.5">총 발현량 (Chr 10q):</text>
        <line x1="95" y1="12" x2="235" y2="12" stroke="#34d399" stroke-width="2"/>
        <text x="165" y="9" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">2n (Neutral, 정상 발현처럼 위장)</text>
      </g>

      <g transform="translate(0, 28)">
        <rect x="0" y="0" width="245" height="24" rx="3" fill="#1e293b"/>
        <text x="10" y="15" fill="#f43f5e" font-size="7.5" font-weight="700">위상 BAF (Chr 10q):</text>
        <line x1="95" y1="8" x2="235" y2="8" stroke="#f43f5e" stroke-width="2"/>
        <line x1="95" y1="20" x2="235" y2="20" stroke="#f43f5e" stroke-width="2"/>
        <text x="165" y="14" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">BAF = 1.0 &amp; 0.0 (완전 대립유전자 소실!)</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• CN-LOH의 분자 생물학적 메커니즘:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">모계 대립유전자 소실 후 부계 대립유전자 중복 (UPD, 2n)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 발현량 전용 도구(inferCNV)의 치명적 한계:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">총 카피 수가 2이므로 정상 2배체로 오분류하는 위음성 발생</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 전립선암 10q PTEN &amp; 17p TP53 좌위 규명:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">돌연변이된 악성 대립유전자 2배 복제 및 종양 억제 완전 차단</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">임상 의의: 숨겨진 LOH 종양 억제 유전자 결손의 정밀 진단</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_clonal_phylogeny_tree_and_subclone_evolution(filepath: str):
    """
    5. Clonal Phylogeny Tree (Neighbor-Joining) & Muller Plot Dynamics
    Truncal vs Subclonal Event Hierarchies -> Single-Cell Distance Matrix & NJ Tree -> Muller Plot Longitudinal Clonal Dynamics.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad05" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">단일세포 계통수 재건 &amp; 뮬러 플롯(Muller Plot) 클론 진화 동역학</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">트렁크(Truncal) vs 서브클론(Subclonal) 이벤트 ➔ 이웃결합(NJ) 계통수 ➔ 치료 경과별 종단 클론 팽창/수축</text>

  <!-- Panel 1: Truncal vs Subclonal -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 트렁크 vs 서브클론 변이 위계</text>

    <!-- Schematic: Clonal Trunk and Branches -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 52)">
      <!-- Trunk Line -->
      <line x1="15" y1="36" x2="95" y2="36" stroke="#38bdf8" stroke-width="4"/>
      <circle cx="15" cy="36" r="4" fill="#38bdf8"/>
      <text x="55" y="28" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">트렁크 (Trunk)</text>
      <text x="55" y="48" text-anchor="middle" fill="#94a3b8" font-size="7">ERG / 8p- / PTEN-</text>

      <!-- Fork Point -->
      <circle cx="95" cy="36" r="4" fill="#fbbf24"/>

      <!-- Branch 1 (Subclone 1: Enz-Resistant) -->
      <line x1="95" y1="36" x2="165" y2="15" stroke="#fbbf24" stroke-width="2.5"/>
      <rect x="165" y="6" width="75" height="18" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="0.8"/>
      <text x="202" y="18" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">AR Amp / 10q LOH</text>

      <!-- Branch 2 (Subclone 2: t-NEPC) -->
      <line x1="95" y1="36" x2="165" y2="52" stroke="#f43f5e" stroke-width="2.5"/>
      <rect x="165" y="44" width="75" height="18" rx="2" fill="#1e293b" stroke="#f43f5e" stroke-width="0.8"/>
      <text x="202" y="56" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">TP53- / RB1- 소실</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 트렁크(Truncal) 원형 변이 (100% 보존):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">모든 종양 세포가 공유하는 TMPRSS2-ERG 융합 및 8p 결손</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 서브클론(Subclonal) 분지 변이:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">국소 미세환경 및 치료 저항성에 따라 차별 획득된 변이</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 종양 내 이종성(ITH)의 유전체 근원:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">동일 환자 내에서도 복수 클론의 공존 및 생존 경쟁</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">계통 법칙: 트렁크 변이는 표적 치료의 가장 안전한 타겟</text>
  </g>

  <!-- Panel 2: NJ Tree Reconstruction -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 이웃결합(NJ) 계통수 재건</text>

    <!-- Schematic: Neighbor-Joining Tree -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Root: Normal Diploid -->
      <rect x="0" y="22" width="45" height="16" rx="2" fill="#475569"/>
      <text x="22" y="33" text-anchor="middle" fill="#fff" font-size="7">Normal 2n</text>

      <line x1="45" y1="30" x2="75" y2="30" stroke="#94a3b8" stroke-width="1.5"/>

      <!-- Internal Node 1 (Founder) -->
      <line x1="75" y1="12" x2="75" y2="48" stroke="#38bdf8" stroke-width="1.5"/>

      <!-- Clade 1: Clone 1 -->
      <line x1="75" y1="12" x2="135" y2="12" stroke="#38bdf8" stroke-width="1.5"/>
      <rect x="135" y="4" width="105" height="15" rx="2" fill="#0284c7"/>
      <text x="187" y="15" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">Clone 1 (Primary AR+)</text>

      <!-- Fork 2 -->
      <line x1="75" y1="48" x2="110" y2="48" stroke="#94a3b8" stroke-width="1.5"/>
      <line x1="110" y1="34" x2="110" y2="60" stroke="#94a3b8" stroke-width="1.5"/>

      <!-- Clade 2: Clone 2 -->
      <line x1="110" y1="34" x2="135" y2="34" stroke="#fbbf24" stroke-width="1.5"/>
      <rect x="135" y="26" width="105" height="15" rx="2" fill="#d97706"/>
      <text x="187" y="37" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">Clone 2 (Enz-Resist AR-Amp)</text>

      <!-- Clade 3: Clone 3 -->
      <line x1="110" y1="60" x2="135" y2="60" stroke="#f43f5e" stroke-width="1.5"/>
      <rect x="135" y="52" width="105" height="15" rx="2" fill="#be123c"/>
      <text x="187" y="63" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">Clone 3 (t-NEPC Double-Loss)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 게놈 복제수 거리 행렬 (Saitou-Nei NJ):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포 간 맨해튼/해밍 거리 기반 Q-행렬 최적화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 진화 속도 이질성(Rate Heterogeneity) 보정:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">서브클론 간 유전체 돌연변이 가속도 차이 정밀 계산</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 뉴윅(Newick) 표준 포맷 트리 출력:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">계통 발생 분기점 및 가지 길이(Branch Length) 정량화</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">수식: Q(i,j) = (n-2)d(i,j) - Σd(i,k) - Σd(j,k)</text>
  </g>

  <!-- Panel 3: Muller Plot Dynamics -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#fbbf24" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#fbbf24" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#fbbf24" font-size="13" font-weight="700">3. 뮬러 플롯 종단적 클론 동역학</text>

    <!-- Schematic: Muller Stream Plot -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>

    <g transform="translate(25, 46)">
      <!-- Timeline axis -->
      <line x1="10" y1="62" x2="240" y2="62" stroke="#475569" stroke-width="1"/>
      <text x="30" y="69" text-anchor="middle" fill="#94a3b8" font-size="6.5">Pre-ADT</text>
      <text x="125" y="69" text-anchor="middle" fill="#94a3b8" font-size="6.5">Enzalutamide</text>
      <text x="215" y="69" text-anchor="middle" fill="#94a3b8" font-size="6.5">t-NEPC</text>

      <!-- Muller Polygon Streams -->
      <!-- Clone 1 (Primary - Dominant then shrank) -->
      <path d="M 10,10 C 60,10 100,35 130,42 C 160,50 200,55 240,58 L 240,62 L 10,62 Z"
        fill="#0284c7" fill-opacity="0.85"/>
      <text x="45" y="40" fill="#fff" font-size="7.5" font-weight="700">Clone 1 (85%)</text>

      <!-- Clone 2 (AR Amp - Rose in Enz) -->
      <path d="M 50,15 C 80,12 110,8 140,10 C 180,12 210,32 240,38 L 240,58 C 200,55 160,50 130,42 Z"
        fill="#fbbf24" fill-opacity="0.85"/>
      <text x="135" y="28" fill="#1e293b" font-size="7.5" font-weight="700">Clone 2 (65%)</text>

      <!-- Clone 3 (t-NEPC - Exploded at end) -->
      <path d="M 140,10 C 170,8 200,6 240,5 L 240,38 C 210,32 180,12 140,10 Z"
        fill="#f43f5e" fill-opacity="0.9"/>
      <text x="210" y="20" fill="#fff" font-size="7.5" font-weight="700">Clone 3 (80%)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 시간적/치료별 서브클론 빈도 변화:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">진단부터 치료 반응 및 최종 재발까지의 클론 적합도 추적</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 엔잘루타마이드 치료 선택압 (Selection):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">기존 Clone 1 박멸 후 AR 증폭 Clone 2의 신속한 대체</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 신경내분비 역분화 클론의 폭발적 우점:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">TP53/RB1 소실 Clone 3가 종양의 80% 이상을 장악</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#fbbf24" font-size="10" font-weight="600">진화 법칙: 치료압력은 저항성 희귀 클론의 선택적 번식을 촉발</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_therapy_selection_and_metastatic_seeding(filepath: str):
    """
    6. Therapy Selection Bottleneck & Metastatic Seeding Patterns
    ADT Bottleneck Selection -> TP53/RB1 Double-Loss t-NEPC Plasticity -> Monoclonal vs Polyclonal Metastatic Seeding.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad06" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrowRose06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">치료 선택압(Therapy Bottleneck) &amp; 전이 파종(Metastatic Seeding) 모형</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">ADT 치료 병목 ➔ TP53/RB1 이중 결손 t-NEPC 계통 가소성 ➔ 단일클론(간) vs 다클론(골) 전이 파종</text>

  <!-- Panel 1: ADT Selection Bottleneck -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 호르몬 차단 요법(ADT) 병목 선택</text>

    <!-- Schematic: Funnel / Bottleneck -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Left: Heterogeneous Primary Tumor -->
      <g transform="translate(0, 5)">
        <rect x="0" y="0" width="55" height="52" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <circle cx="15" cy="15" r="4" fill="#38bdf8"/>
        <circle cx="35" cy="15" r="4" fill="#38bdf8"/>
        <circle cx="20" cy="35" r="4" fill="#38bdf8"/>
        <circle cx="45" cy="40" r="4" fill="#f43f5e"/>
        <text x="27" y="28" text-anchor="middle" fill="#38bdf8" font-size="6.5">감수성 95%</text>
      </g>

      <!-- Bottleneck funnel -->
      <path d="M 60,8 L 105,25 L 105,37 L 60,54 Z" fill="#334155" fill-opacity="0.5" stroke="#64748b" stroke-width="1"/>
      <text x="82" y="34" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">ADT Bottleneck</text>

      <!-- Center: MRD (Minimal Residual Disease) -->
      <g transform="translate(112, 18)">
        <rect x="0" y="0" width="25" height="26" rx="2" fill="#be123c" stroke="#f43f5e" stroke-width="1"/>
        <circle cx="12" cy="13" r="4" fill="#f43f5e"/>
        <text x="12" y="35" text-anchor="middle" fill="#f43f5e" font-size="6.5">MRD</text>
      </g>

      <!-- Arrow -->
      <line x1="142" y1="31" x2="165" y2="31" stroke="#f43f5e" stroke-width="2" marker-end="url(#arrowRose06)"/>

      <!-- Right: Expanded Resistant Clones -->
      <g transform="translate(175, 5)">
        <rect x="0" y="0" width="70" height="52" rx="4" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
        <circle cx="18" cy="15" r="4" fill="#f43f5e"/>
        <circle cx="42" cy="15" r="4" fill="#f43f5e"/>
        <circle cx="22" cy="35" r="4" fill="#f43f5e"/>
        <circle cx="52" cy="35" r="4" fill="#fbbf24"/>
        <text x="35" y="48" text-anchor="middle" fill="#f43f5e" font-size="6.5" font-weight="700">저항성 재발 (100%)</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 강력한 약리학적 선택압 (Darwinian Selection):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">안드로겐 고갈 시 호르몬 의존성 종양 세포 95% 이상 사멸</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 미세 잔존 질환(MRD) 생존 창 형성:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">사전 존재하는 극소수 비전형 클론이 치료 병목 통과</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 유전적 다양성 단순화 후 폭발:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">병목을 통과한 단일 저항성 서브클론의 적응적 재성장</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">치료 역설: 고강도 단일 치료는 가장 치명적인 변종을 선별</text>
  </g>

  <!-- Panel 2: t-NEPC Lineage Plasticity -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">2. t-NEPC 계통 가소성 &amp; 역분화</text>

    <!-- Schematic: Luminal AR+ to NEPC transdifferentiation -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>

    <g transform="translate(25, 50)">
      <!-- Luminal cell -->
      <g transform="translate(0, 0)">
        <circle cx="28" cy="24" r="18" fill="#0284c7" fill-opacity="0.7" stroke="#38bdf8" stroke-width="1.5"/>
        <text x="28" y="22" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">Luminal</text>
        <text x="28" y="32" text-anchor="middle" fill="#bae6fd" font-size="6.5">AR+/PSA+</text>
      </g>

      <!-- Double Loss Bridge -->
      <g transform="translate(68, 6)">
        <rect x="0" y="0" width="112" height="38" rx="4" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
        <text x="56" y="14" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">TP53 &amp; RB1 동시 소실</text>
        <text x="56" y="24" text-anchor="middle" fill="#fbbf24" font-size="7">EZH2 / SOX2 후성 리모델링</text>
        <text x="56" y="33" text-anchor="middle" fill="#94a3b8" font-size="6.5">계통 분화 잠금 해제 (Plasticity)</text>
      </g>

      <!-- t-NEPC cell -->
      <g transform="translate(195, 0)">
        <circle cx="28" cy="24" r="18" fill="#be123c" fill-opacity="0.85" stroke="#f43f5e" stroke-width="1.5"/>
        <text x="28" y="22" text-anchor="middle" fill="#fff" font-size="7.5" font-weight="700">t-NEPC</text>
        <text x="28" y="32" text-anchor="middle" fill="#fecdd3" font-size="6.5">SYP+/CHGA+</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• TP53 &amp; RB1 동시 소실(Double-Loss)의 격변:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">세포 분화 계통 잠금장치가 해제되어 다능성 상태 획득</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 신경내분비 형질 전환 (Neuroendocrine Shift):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">상피성 안드로겐 마커 소실 및 ASCL1/NEUROD1 발현</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 안드로겐 수용체(AR) 신호 독립성:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">AR이 결핍되어 기존 호르몬 치료제에 완전히 불응</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">임상 경고: mCRPC 환자의 15~20%에서 치명적 t-NEPC 변환</text>
  </g>

  <!-- Panel 3: Monoclonal vs Polyclonal Seeding -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 단일클론 vs 다클론 전이 파종</text>

    <!-- Schematic: Monoclonal to Liver vs Polyclonal to Bone -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Primary Tumor Source -->
      <g transform="translate(0, 10)">
        <rect x="0" y="0" width="45" height="42" rx="4" fill="#1e293b" stroke="#64748b" stroke-width="1"/>
        <circle cx="15" cy="15" r="4" fill="#f43f5e"/>
        <circle cx="30" cy="15" r="4" fill="#fbbf24"/>
        <circle cx="22" cy="30" r="4" fill="#38bdf8"/>
        <text x="22" y="49" text-anchor="middle" fill="#94a3b8" font-size="6">Primary</text>
      </g>

      <!-- Route 1: Monoclonal Seeding to Liver -->
      <g transform="translate(55, 6)">
        <line x1="0" y1="12" x2="50" y2="12" stroke="#f43f5e" stroke-width="1.5" stroke-dasharray="2,2"/>
        <rect x="55" y="0" width="130" height="24" rx="3" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
        <circle cx="68" cy="12" r="4" fill="#f43f5e"/>
        <circle cx="80" cy="12" r="4" fill="#f43f5e"/>
        <text x="135" y="11" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">간 전이: 단일클론 파종</text>
        <text x="135" y="20" text-anchor="middle" fill="#94a3b8" font-size="6">단일 저항성 클론 단독 증식</text>
      </g>

      <!-- Route 2: Polyclonal Seeding to Bone -->
      <g transform="translate(55, 34)">
        <line x1="0" y1="12" x2="50" y2="12" stroke="#34d399" stroke-width="1.5"/>
        <rect x="55" y="0" width="130" height="26" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
        <circle cx="68" cy="13" r="4" fill="#f43f5e"/>
        <circle cx="78" cy="13" r="4" fill="#fbbf24"/>
        <circle cx="88" cy="13" r="4" fill="#38bdf8"/>
        <text x="138" y="11" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">골 전이: 다클론 파종</text>
        <text x="138" y="21" text-anchor="middle" fill="#94a3b8" font-size="6">다중 서브클론 집단 협력</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 단일클론 파종 (Monoclonal Seeding - 간/내장):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일 혈중 종양 세포(CTC)가 혈류를 거쳐 단독 정착</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 다클론 집단 파종 (Polyclonal Seeding - 골 전이):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">다중 서브클론 클러스터가 골수 미세환경에 공동 파종</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 전이소 간 전이 (Metastasis-to-Metastasis):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">골전이소에서 재파종된 클론이 3차 장기로 다방향 확산</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">치료 시사점: 다클론 골전이는 복합 병용 표적 치료 필수</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_prostate_cancer_hallmark_cna_landscape.svg", make_01_prostate_cancer_hallmark_cna_landscape),
        ("02_infercnv_hmm_smoothing_workflow.svg", make_02_infercnv_hmm_smoothing_workflow),
        ("03_copykat_gmm_diploid_aneuploid_separation.svg", make_03_copykat_gmm_diploid_aneuploid_separation),
        ("04_numbat_haplotype_phased_baf_and_cnloh.svg", make_04_numbat_haplotype_phased_baf_and_cnloh),
        ("05_clonal_phylogeny_tree_and_subclone_evolution.svg", make_05_clonal_phylogeny_tree_and_subclone_evolution),
        ("06_therapy_selection_and_metastatic_seeding.svg", make_06_therapy_selection_and_metastatic_seeding),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 04 Cancer Genomics & Clones publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

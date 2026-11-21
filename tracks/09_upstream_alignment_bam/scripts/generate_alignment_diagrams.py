#!/usr/bin/env python3
"""
tracks/09_upstream_alignment_bam/scripts/generate_alignment_diagrams.py
=============================================================================
Generates 6 publication-grade, high-resolution vector SVG diagrams for Track 09:
Single-Cell Upstream Alignment, BAM Processing, and Splice Variant Mining.

Diagrams:
1. 01_fastq_read_structure_and_barcode_correction.svg
2. 02_starsolo_suffix_array_spliced_alignment.svg
3. 03_alevin_fry_rad_file_and_splici_index.svg
4. 04_umi_directional_adjacency_graph_deduplication.svg
5. 05_bam_binary_record_cigar_and_tags.svg
6. 06_ar_v7_splice_junction_and_tmprss2_erg_fusion.svg
"""

import os
import xml.etree.ElementTree as ET


def make_01_fastq_read_structure_and_barcode_correction(filepath: str):
    """
    1. FASTQ Read Structure, 1-Hamming Error Correction & EmptyDrops Cell Filtering
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
    <marker id="arrEmerald01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrAmber01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
    <marker id="arrRose01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">FASTQ 리드 구조 분석, 1-Hamming 바코드 오류 교정 및 EmptyDrops 세포 식별</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">10x R1/R2 이중 가닥 구조 ➔ 737K 화이트리스트 1-Hamming 해시 복원 ➔ Dirichlet-Multinomial 배경 액적 분리</text>

  <!-- Panel 1: FASTQ Read Architecture -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. 10x Genomics v3 R1/R2 리드 아키텍처</text>

    <!-- Schematic: Read 1 -->
    <rect x="15" y="42" width="270" height="42" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <text x="22" y="55" fill="#94a3b8" font-size="8.5" font-weight="700">Read 1 (28 bp Short Read):</text>
    <!-- Barcode 16bp -->
    <rect x="22" y="60" width="130" height="18" rx="3" fill="#059669" fill-opacity="0.25" stroke="#34d399" stroke-width="1"/>
    <text x="87" y="73" text-anchor="middle" fill="#34d399" font-size="9" font-weight="700">Cell Barcode (16 bp)</text>
    <!-- UMI 12bp -->
    <rect x="156" y="60" width="120" height="18" rx="3" fill="#d97706" fill-opacity="0.25" stroke="#fbbf24" stroke-width="1"/>
    <text x="216" y="73" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="700">UMI (12 bp)</text>

    <!-- Schematic: Read 2 -->
    <rect x="15" y="92" width="270" height="42" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <text x="22" y="105" fill="#94a3b8" font-size="8.5" font-weight="700">Read 2 (91+ bp Paired cDNA):</text>
    <rect x="22" y="110" width="254" height="18" rx="3" fill="#0284c7" fill-opacity="0.25" stroke="#38bdf8" stroke-width="1"/>
    <text x="149" y="123" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">cDNA Transcript Insert (전사항체 서열)</text>

    <!-- Bullet description -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 바코드 다양성: 4^16 = 4.29 × 10^9 조합 공간</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">737,280개 검증된 젤비즈 비드 화이트리스트 매핑</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 고유 분자 식별자: 4^12 = 1.67 × 10^7 UMI 공간</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">PCR 중복 증폭 산물과 단일 전사체 1:1 절대 계수</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• Phred 품질 점수 (Q30 &gt; 99.9% 신뢰도):</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">Illumina NovaSeq 광학적 염기 판독 오류 감지</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">Dual-index 비대칭 리드 구조 ➔ 세포 및 분자 태깅</text>
  </g>

  <!-- Panel 2: 1-Hamming Error Correction -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">2. 1-Hamming 바코드 오류 교정 해시 그래프</text>

    <!-- Graph Schema -->
    <rect x="15" y="42" width="280" height="92" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Whitelist Node -->
    <circle cx="65" cy="88" r="22" fill="#065f46" stroke="#34d399" stroke-width="1.5"/>
    <text x="65" y="85" text-anchor="middle" fill="#a7f3d0" font-size="8" font-weight="700">Whitelist</text>
    <text x="65" y="97" text-anchor="middle" fill="#34d399" font-size="7.5">ACGT...A</text>

    <!-- 1bp mutant 1 -->
    <circle cx="165" cy="62" r="15" fill="#1e293b" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="165" y="65" text-anchor="middle" fill="#fbbf24" font-size="7.5">TCGT (1bp)</text>

    <!-- 1bp mutant 2 -->
    <circle cx="165" cy="114" r="15" fill="#1e293b" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="165" y="117" text-anchor="middle" fill="#fbbf24" font-size="7.5">AGGT (1bp)</text>

    <!-- Corrected Arrow -->
    <path d="M148,65 L90,82" stroke="#34d399" stroke-width="1.5" stroke-dasharray="3,2" marker-end="url(#arrEmerald01)"/>
    <path d="M148,111 L90,94" stroke="#34d399" stroke-width="1.5" stroke-dasharray="3,2" marker-end="url(#arrEmerald01)"/>

    <!-- Unresolvable collision node -->
    <circle cx="255" cy="88" r="15" fill="#1e293b" stroke="#f43f5e" stroke-width="1.2"/>
    <text x="255" y="85" text-anchor="middle" fill="#f43f5e" font-size="7">Ambiguous</text>
    <text x="255" y="96" text-anchor="middle" fill="#f43f5e" font-size="7">Dist=1 to 2 WL</text>
    <path d="M238,88 L185,88" stroke="#f43f5e" stroke-width="1.2" stroke-dasharray="2,2"/>
    <text x="215" y="82" text-anchor="middle" fill="#f43f5e" font-size="7">충돌 폐기</text>

    <!-- Details -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 1-Hamming 근접 이웃: 3 × 16 = 48개 변이 공간</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">사전 인덱싱 해시 테이블로 O(1) 복구 탐색</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• 단일 복구율 &gt; 96.5% 달성:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">PCR 1bp 점 돌연변이 및 염기 판독 오류 완벽 복구</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• 다중 화이트리스트 충돌(Collision) 방어:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">두 개 이상 정품 바코드와 동등 거리 시 안전 폐기</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.25"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">손실 리드의 8~12%를 복원하여 세포 검출 감도 극대화</text>
  </g>

  <!-- Panel 3: EmptyDrops & Cell Calling -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">3. EmptyDrops &amp; Knee-point 세포 여과</text>

    <!-- Chart: Barcode Rank vs UMI Count -->
    <rect x="15" y="42" width="280" height="92" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Axes -->
    <line x1="40" y1="120" x2="275" y2="120" stroke="#64748b" stroke-width="1"/>
    <line x1="40" y1="120" x2="40" y2="52" stroke="#64748b" stroke-width="1"/>
    <text x="270" y="130" text-anchor="end" fill="#64748b" font-size="7">Log Barcode Rank</text>
    <text x="35" y="52" text-anchor="end" fill="#64748b" font-size="7">Log UMI</text>

    <!-- Knee Curve -->
    <path d="M42,56 C70,58 100,60 130,70 C150,80 165,100 185,114 C210,117 250,118 270,119" fill="none" stroke="#38bdf8" stroke-width="2"/>
    
    <!-- Knee & Inflection Points -->
    <circle cx="130" cy="70" r="3.5" fill="#34d399"/>
    <text x="135" y="67" fill="#34d399" font-size="7.5" font-weight="700">Knee pt (진짜 세포)</text>
    
    <circle cx="185" cy="114" r="3.5" fill="#f43f5e"/>
    <text x="190" y="110" fill="#f43f5e" font-size="7.5" font-weight="700">Inflection (배경)</text>

    <!-- Ambient Threshold Line -->
    <line x1="40" y1="114" x2="275" y2="114" stroke="#f43f5e" stroke-width="1" stroke-dasharray="2,2"/>

    <!-- Bullets -->
    <text x="14" y="152" fill="#cbd5e1" font-size="10" font-weight="600">• 단순 임계값(Knee-point)의 한계 극복:</text>
    <text x="24" y="167" fill="#94a3b8" font-size="9">소형 세포(T세포 등)의 UMI 과소 계수 탈락 방지</text>
    <text x="14" y="185" fill="#cbd5e1" font-size="10" font-weight="600">• Dirichlet-Multinomial 배경 RNA 프로파일링:</text>
    <text x="24" y="200" fill="#94a3b8" font-size="9">세포 없는 Ambient RNA 액적의 전사체 분포 모델링</text>
    <text x="14" y="218" fill="#cbd5e1" font-size="10" font-weight="600">• Benjamini-Hochberg FDR &lt; 0.01 검정:</text>
    <text x="24" y="233" fill="#94a3b8" font-size="9">유의미한 전사체 편향 액적만을 진성 세포로 판별</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#78350f" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">파편화된 Ambient RNA 오염을 배제하고 희귀 세포 구제</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_starsolo_suffix_array_spliced_alignment(filepath: str):
    """
    2. STARsolo Suffix Array (BWT) Indexing & Spliced Intron Dynamic Programming
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
    <marker id="arrEmerald02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrAmber02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">STARsolo 스플라이싱 정렬: 접미사 배열(Suffix Array) BWT 색인 및 MMP 시드 확장</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">O(L log N) 접미사 배열 탐색 ➔ 최대 매핑 접두사(MMP) 분할 ➔ 정밀 동적 계획법(DP) 스플라이스 인트론 스코어링</text>

  <!-- Panel 1: Suffix Array Indexing -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. Suffix Array (SA) &amp; BWT 유전체 색인</text>

    <!-- Schematic: Suffix Array Table -->
    <rect x="15" y="42" width="270" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="12" fill="#64748b" font-size="7.5" font-weight="700">Idx</text>
      <text x="25" y="12" fill="#64748b" font-size="7.5" font-weight="700">SA[i]</text>
      <text x="65" y="12" fill="#64748b" font-size="7.5" font-weight="700">Suffix String (사전순 정렬)</text>

      <line x1="0" y1="16" x2="250" y2="16" stroke="#334155" stroke-width="0.8"/>
      <text x="0" y="28" fill="#94a3b8" font-size="7.5">0</text>
      <text x="25" y="28" fill="#38bdf8" font-size="7.5">14,208</text>
      <text x="65" y="28" fill="#cbd5e1" font-size="7.5">AAGCTTGACGT...</text>

      <text x="0" y="42" fill="#94a3b8" font-size="7.5">1</text>
      <text x="25" y="42" fill="#38bdf8" font-size="7.5">892,104</text>
      <text x="65" y="42" fill="#34d399" font-size="7.5">ACGTTAGCTA... (MMP Match)</text>

      <text x="0" y="56" fill="#94a3b8" font-size="7.5">2</text>
      <text x="25" y="42" fill="#38bdf8" font-size="7.5">3,110,4</text>
      <text x="65" y="56" fill="#cbd5e1" font-size="7.5">AGCTAGGTCC...</text>

      <text x="0" y="70" fill="#94a3b8" font-size="7.5">...</text>
      <text x="25" y="70" fill="#64748b" font-size="7.5">3.2 × 10^9</text>
      <text x="65" y="70" fill="#64748b" font-size="7.5">TTTTTT...$ (GRCh38 인간 유전체)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 메모리 상주 접미사 배열 (~32 GB RAM):</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">BWT 압축 없이 전체 SA를 RAM에 적재해 즉각 탐색</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 이진 탐색 속도 O(L log N):</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">100 bp cDNA 리드를 수 마이크로초 내에 전역 조회</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 접미사 길이 k-mer 프리픽스 캐싱:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">첫 14bp의 SA 범위를 O(1) 직접 참조 테이블로 가속</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">전역 색인 기반 무제한 갭 스플라이스 시드 발굴</text>
  </g>

  <!-- Panel 2: Maximal Mappable Prefix (MMP) -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">2. Maximal Mappable Prefix (MMP) 시드 분할</text>

    <!-- MMP Split Diagram -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Full Read -->
    <rect x="25" y="52" width="260" height="18" rx="3" fill="#1e293b" stroke="#64748b" stroke-width="1"/>
    <text x="155" y="64" text-anchor="middle" fill="#94a3b8" font-size="8">100 bp Spliced RNA Read</text>

    <!-- MMP 1 (Exon A) -->
    <rect x="25" y="78" width="95" height="18" rx="3" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="1.2"/>
    <text x="72" y="90" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">MMP 1 (38bp)</text>
    <!-- Intron Gap -->
    <line x1="120" y1="87" x2="180" y2="87" stroke="#fbbf24" stroke-width="1.5" stroke-dasharray="3,2"/>
    <text x="150" y="82" text-anchor="middle" fill="#fbbf24" font-size="7.5">Intron 12 kb</text>
    <!-- MMP 2 (Exon B) -->
    <rect x="180" y="78" width="105" height="18" rx="3" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="232" y="90" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">MMP 2 (62bp)</text>

    <!-- Split Arrow -->
    <path d="M120,68 L120,76" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrCyan02)"/>
    <text x="120" y="112" text-anchor="middle" fill="#f43f5e" font-size="7">불일치 지점에서 분할</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 순차적 접두사 연장(Prefix Extension):</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">유전체와 완벽 일치하는 최장 부분 서열을 연속 탐색</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 엑손-엑손 접합부(Junction) 자동 검출:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">스플라이스 사이트에서 미스매치 시 즉각 차기 시드 생성</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 초장거리 키메라/구조 변이 포착:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">수십 kb 인트론뿐 아니라 염색체 간 융합 시드 동시 포섭</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.25"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">주석(Annotation) 없는 신규 스플라이싱 접합부 발굴</text>
  </g>

  <!-- Panel 3: Spliced DP Scoring -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">3. 동적 계획법(DP) 스플라이스 스코어링</text>

    <!-- Motifs & Scoring Box -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="14" fill="#fbbf24" font-size="8.5" font-weight="700">정규 모티프:</text>
      <text x="70" y="14" fill="#34d399" font-size="8.5" font-weight="700">GT (Donor) ... AG (Acceptor) [+20점]</text>
      
      <text x="0" y="32" fill="#fbbf24" font-size="8.5" font-weight="700">비정규 모티프:</text>
      <text x="70" y="32" fill="#cbd5e1" font-size="8.5">GC...AG [+15점] / AT...AC [+10점]</text>

      <text x="0" y="50" fill="#fbbf24" font-size="8.5" font-weight="700">인트론 페널티:</text>
      <text x="70" y="50" fill="#f43f5e" font-size="8.5">Score -= log(L_intron) × Penalty</text>

      <text x="0" y="68" fill="#fbbf24" font-size="8.5" font-weight="700">STARsolo 출력:</text>
      <text x="70" y="68" fill="#38bdf8" font-size="8.5" font-weight="700">Gene / GeneFull / SJ.out.tab 매트릭스</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• Needleman-Wunsch 변형 정밀 정렬:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">시드 간 공백을 정규 스플라이스 시그널로 보정 연결</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• --soloType CB_UMI_Simple 네이티브 통합:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">BAM 중간 파일 없이 단일세포 카운트 매트릭스 직출력</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• CellRanger 대비 10배 정렬 속도:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">1억 리드 정렬 및 CB/UMI 카운팅을 20분 내 완료</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#78350f" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">미세 엑손 및 대형 인트론 고정밀 스플라이싱 정량화</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_alevin_fry_rad_file_and_splici_index(filepath: str):
    """
    3. Alevin-fry Quasi-mapping, Splici Index & RAD Binary Stream to USA Matrix
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
    <marker id="arrEmerald03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrPurple03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">Alevin-fry 유사 정렬 아키텍처: Splici 색인 및 RAD 바이너리 스트림 정량화</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">k-mer Quasi-mapping ➔ Splici(전사체+인트론) 참조 유전체 ➔ RAD 바이너리 청크 및 USA 카운트 행렬</text>

  <!-- Panel 1: Splici Index Structure -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. Splici (Spliced + Intron) 인덱스 구조</text>

    <!-- Schematic: Splici Components -->
    <rect x="15" y="42" width="270" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Spliced cDNA -->
    <rect x="25" y="50" width="250" height="18" rx="3" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
    <text x="150" y="62" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">Spliced Transcripts (성숙 mRNA 엑손 연합체)</text>

    <!-- Intron with Flank -->
    <rect x="25" y="74" width="35" height="18" rx="2" fill="#d97706" fill-opacity="0.3" stroke="#fbbf24" stroke-width="1"/>
    <text x="42" y="86" text-anchor="middle" fill="#fbbf24" font-size="7">Flank</text>

    <rect x="63" y="74" width="174" height="18" rx="2" fill="#7e22ce" fill-opacity="0.3" stroke="#c084fc" stroke-width="1"/>
    <text x="150" y="86" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">Intronic Sequence (미성숙 pre-mRNA)</text>

    <rect x="240" y="74" width="35" height="18" rx="2" fill="#d97706" fill-opacity="0.3" stroke="#fbbf24" stroke-width="1"/>
    <text x="257" y="86" text-anchor="middle" fill="#fbbf24" font-size="7">Flank</text>

    <text x="150" y="106" text-anchor="middle" fill="#94a3b8" font-size="7.5">Flanking Length = Read_Length - 1 (경계 리드 포착)</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 단일핵(snRNA-seq) 최적화 참조 설계:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">핵 내부 미성숙 pre-mRNA 인트론 신호를 정확히 포섭</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• RNA Velocity 동역학 데이터 기반 마련:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">스플라이싱 전(Unspliced) vs 후(Spliced) 비율 추적</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 의사 유전자(Pseudogene) 오정렬 차단:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">게놈 전체 인트론 패딩으로 미스매핑 비율 대폭 감소</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">성숙 전사체와 인트론을 동등한 참조체로 색인화</text>
  </g>

  <!-- Panel 2: Quasi-mapping & Equivalence Classes -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12" font-weight="700">2. k-mer 유사 매핑 및 동등성 클래스 (EqClass)</text>

    <!-- Quasi Diagram -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- k-mers -->
    <rect x="25" y="52" width="45" height="16" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="47" y="63" text-anchor="middle" fill="#38bdf8" font-size="7">k1 (31bp)</text>
    <rect x="75" y="52" width="45" height="16" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="97" y="63" text-anchor="middle" fill="#38bdf8" font-size="7">k2 (31bp)</text>
    <rect x="125" y="52" width="45" height="16" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
    <text x="147" y="63" text-anchor="middle" fill="#38bdf8" font-size="7">k3 (31bp)</text>

    <!-- Target transcripts EqClass -->
    <g transform="translate(25, 78)">
      <rect width="250" height="22" rx="3" fill="#581c87" fill-opacity="0.3" stroke="#c084fc" stroke-width="1"/>
      <text x="125" y="15" text-anchor="middle" fill="#e9d5ff" font-size="8">Equivalence Class #412: {TX1_Spliced, TX1_Intron}</text>
    </g>
    <text x="155" y="118" text-anchor="middle" fill="#94a3b8" font-size="7.5">CIGAR 계산 생략 ➔ 전사체 집합 매핑 (15배 고속화)</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• CIGAR 정렬 없는 초고속 유사 매핑:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">정밀 염기 정렬 대신 k-mer 일치 전사체 라벨 집합 추출</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 다중 매핑 전사체(Multimapping) 처리:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">임의 배정 대신 동등성 클래스 단위로 모호성 유지</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 극단적인 I/O 및 CPU 부하 절감:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">기가바이트 단위 대용량 BAM 생성 없이 스트리밍</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#581c87" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#c084fc" font-size="9" font-weight="600">정렬 병목을 해소하는 경량 k-mer 해시 정량화</text>
  </g>

  <!-- Panel 3: RAD File & USA Matrix -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. RAD 스트림 및 USA (U-S-A) 카운트 행렬</text>

    <!-- RAD Format Box -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="12" fill="#38bdf8" font-size="7.5" font-weight="700">RAD Record: [CB: 64b] [UMI: 32b] [EqClass: 32b]</text>
      <!-- Binary Chunk Stream -->
      <rect x="0" y="18" width="80" height="16" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="40" y="29" text-anchor="middle" fill="#38bdf8" font-size="7">Chunk 1 (CB-A)</text>
      <rect x="85" y="18" width="80" height="16" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="125" y="29" text-anchor="middle" fill="#38bdf8" font-size="7">Chunk 2 (CB-A)</text>
      <rect x="170" y="18" width="80" height="16" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="0.8"/>
      <text x="210" y="29" text-anchor="middle" fill="#34d399" font-size="7">Chunk 3 (CB-B)</text>

      <!-- USA Matrix Breakdown -->
      <line x1="0" y1="42" x2="250" y2="42" stroke="#334155" stroke-width="0.8"/>
      <text x="0" y="55" fill="#34d399" font-size="8" font-weight="700">U (Unspliced):</text>
      <text x="80" y="55" fill="#cbd5e1" font-size="8">인트론 단독 매핑 (신생 전사체)</text>
      <text x="0" y="68" fill="#38bdf8" font-size="8" font-weight="700">S (Spliced):</text>
      <text x="80" y="68" fill="#cbd5e1" font-size="8">엑손 접합부 매핑 (성숙 mRNA)</text>
      <text x="0" y="80" fill="#fbbf24" font-size="8" font-weight="700">A (Ambiguous):</text>
      <text x="80" y="80" fill="#cbd5e1" font-size="8">엑손 내부 매핑 (U/S 공유)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• RAD (Reduced Alignment Data) 바이너리 압축:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">BAM 대비 90% 작은 파일 크기로 디스크 I/O 극소화</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 세포 바코드별 청크 스트리밍 정렬:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">Alevin-fry 엔진이 병렬 코어에서 메모리 누수 없이 집계</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• USA 3차원 카운트 행렬 동시 도출:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">스플라이싱 역동성 분석과 표준 유전자 발현 정량 동시 완성</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.25"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">BAM 불필요 ➔ 초고속 메모리 절감 차세대 업스트림 표준</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_umi_directional_adjacency_graph_deduplication(filepath: str):
    """
    4. PCR Amplification Error Branching Tree, UMI Directional Adjacency Network & 2N-1 Threshold
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
    <marker id="arrRose04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">UMI 중복 제거(Deduplication): PCR 증폭 오류 트리 및 방향성 인접 그래프</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">단순 고유 카운팅의 계수 인플레이션 한계 ➔ 1-Hamming 인접 네트워크 ➔ N_parent ≥ 2N_child - 1 임계 필터링</text>

  <!-- Panel 1: PCR Error Branching Tree -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12" font-weight="700">1. PCR 증폭 오류 및 UMI 인플레이션</text>

    <!-- Branching Tree Diagram -->
    <rect x="15" y="42" width="270" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Parent Root Node -->
    <circle cx="65" cy="86" r="20" fill="#065f46" stroke="#34d399" stroke-width="1.5"/>
    <text x="65" y="83" text-anchor="middle" fill="#a7f3d0" font-size="7.5" font-weight="700">True UMI</text>
    <text x="65" y="94" text-anchor="middle" fill="#34d399" font-size="8" font-weight="700">N = 84</text>

    <!-- Child Node 1 -->
    <circle cx="175" cy="62" r="14" fill="#1e293b" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="175" y="60" text-anchor="middle" fill="#fbbf24" font-size="7">Mutant 1</text>
    <text x="175" y="70" text-anchor="middle" fill="#fbbf24" font-size="7.5">N = 6</text>

    <!-- Child Node 2 -->
    <circle cx="175" cy="110" r="14" fill="#1e293b" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="175" y="108" text-anchor="middle" fill="#fbbf24" font-size="7">Mutant 2</text>
    <text x="175" y="118" text-anchor="middle" fill="#fbbf24" font-size="7.5">N = 2</text>

    <!-- Branch lines -->
    <path d="M85,82 L158,64" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="118" y="68" fill="#fbbf24" font-size="6.5">1bp PCR error</text>
    <path d="M85,90 L158,108" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="118" y="108" fill="#fbbf24" font-size="6.5">1bp error</text>

    <!-- Independent Collision Node -->
    <circle cx="250" cy="86" r="16" fill="#1e293b" stroke="#38bdf8" stroke-width="1.2"/>
    <text x="250" y="83" text-anchor="middle" fill="#38bdf8" font-size="7">Indep UMI</text>
    <text x="250" y="93" text-anchor="middle" fill="#38bdf8" font-size="7.5">N = 35</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• PCR 중합효소 오류율 (~10^-4 per bp):</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">12bp UMI가 증폭 사이클을 거치며 돌연변이 자손 파생</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• Naive Unique 카운팅의 치명적 결함:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">1개 분자가 3개로 카운트되어 발현량이 25~40% 과대평가</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 시퀀싱 뎁스 증가에 따른 왜곡 심화:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">시퀀싱 리드가 많아질수록 위양성 UMI 수가 폭발</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.25"/>
    <text x="150" y="259" text-anchor="middle" fill="#f43f5e" font-size="9" font-weight="600">단순 고유 계수는 전사체 정량의 신뢰도를 파괴</text>
  </g>

  <!-- Panel 2: 4 Deduplication Algorithms -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">2. UMI 디듀플리케이션 4대 알고리즘</text>

    <!-- Comparison Table / Diagram -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="12" fill="#f43f5e" font-size="7.5" font-weight="700">1. Unique (Naive):</text>
      <text x="90" y="12" fill="#cbd5e1" font-size="7.5">단순 유일 서열 계수 ➔ 인플레이션 +30%</text>

      <text x="0" y="28" fill="#fbbf24" font-size="7.5" font-weight="700">2. Cluster (CC):</text>
      <text x="90" y="28" fill="#cbd5e1" font-size="7.5">1-Hamming 연결요소 ➔ 과대 병합 위험</text>

      <text x="0" y="44" fill="#a855f7" font-size="7.5" font-weight="700">3. Adjacency:</text>
      <text x="90" y="44" fill="#cbd5e1" font-size="7.5">최다 빈도 노드 주변 1-step 통합</text>

      <text x="0" y="60" fill="#34d399" font-size="7.5" font-weight="700">4. Directional:</text>
      <text x="90" y="60" fill="#34d399" font-size="7.5" font-weight="700">N_A ≥ 2N_B - 1 판정 유향 그래프 (표준)</text>

      <line x1="0" y1="68" x2="250" y2="68" stroke="#334155" stroke-width="0.8"/>
      <text x="125" y="79" text-anchor="middle" fill="#38bdf8" font-size="7.5">UMI-tools / CellRanger / STARsolo 핵심 표준 엔진</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 연결 요소(Cluster)의 브릿징 오류:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">독립적인 두 고발현 UMI가 우연한 에러로 하나로 융합</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 방향성 인접(Directional Adjacency):</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">카운트 빈도 격차를 조건부 확률로 모델링하여 분기 차단</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 정확도 및 리콜률 최상위 수렴:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">위양성 노이즈는 99% 제거하고 진성 분자는 99.8% 보존</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="155" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">독립 분자 과소 계수 없이 PCR 증폭 오류만 정밀 제거</text>
  </g>

  <!-- Panel 3: N >= 2N - 1 Directional Edge Rule -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. N_parent ≥ 2N_child - 1 방향성 엣지 규칙</text>

    <!-- Network Graph Schema -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Node A (Parent, N=100) -->
    <circle cx="65" cy="86" r="22" fill="#065f46" stroke="#34d399" stroke-width="1.5"/>
    <text x="65" y="83" text-anchor="middle" fill="#a7f3d0" font-size="8" font-weight="700">Node A</text>
    <text x="65" y="94" text-anchor="middle" fill="#34d399" font-size="8" font-weight="700">N_A = 100</text>

    <!-- Node B (Valid child, N=12) -->
    <circle cx="165" cy="64" r="16" fill="#1e293b" stroke="#34d399" stroke-width="1.2"/>
    <text x="165" y="61" text-anchor="middle" fill="#a7f3d0" font-size="7">Node B</text>
    <text x="165" y="71" text-anchor="middle" fill="#34d399" font-size="7.5">N_B = 12</text>

    <!-- Directed Edge A -> B -->
    <path d="M88,77 L147,67" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald04)"/>
    <text x="125" y="70" fill="#34d399" font-size="6.5">100 ≥ 2(12)-1=23 [성립]</text>

    <!-- Node C (Independent, N=60) -->
    <circle cx="250" cy="86" r="20" fill="#1e293b" stroke="#f43f5e" stroke-width="1.5"/>
    <text x="250" y="83" text-anchor="middle" fill="#fca5a5" font-size="7.5">Node C</text>
    <text x="250" y="94" text-anchor="middle" fill="#f43f5e" font-size="8">N_C = 60</text>

    <!-- No Edge between B and C -->
    <line x1="183" y1="70" x2="230" y2="82" stroke="#64748b" stroke-width="1" stroke-dasharray="2,2"/>
    <text x="208" y="73" fill="#f43f5e" font-size="6.5">12 &lt; 2(60)-1=119 [절단]</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 유향 엣지 판정 규칙:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">Hamming 거리 = 1 이고 N_A ≥ 2N_B - 1 일 때만 A ➔ B</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 자식 노드의 카운트 조건부 흡수:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">자식의 발생 횟수가 부모의 절반 이하일 때만 에러로 병합</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 독립 전사체 충돌 보호:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">우연히 1bp 차이나는 고발현 독립 분자는 각각 온전히 보존</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.25"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">수학적 확률 기반 완벽한 분자 계수 정밀도 보장</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_bam_binary_record_cigar_and_tags(filepath: str):
    """
    5. BAM Binary Record Structure, CIGAR Operators & Single-Cell SAM Tags
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
    <marker id="arrEmerald05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrAmber05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">BAM 바이너리 정렬 레코드 구조: CIGAR 연산자 및 단일세포 전용 태그</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">BGZF 블록 압축 ➔ CIGAR 정렬 기하학(M, N, I, D, S) ➔ 단일세포 식별 태그(CB, UB, GN, RE)</text>

  <!-- Panel 1: BAM Binary Architecture -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">1. BAM 바이너리 레코드 메모리 구조</text>

    <!-- Binary Field Layout -->
    <rect x="15" y="42" width="270" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(22, 48)">
      <rect x="0" y="0" width="55" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="27" y="12" text-anchor="middle" fill="#38bdf8" font-size="7">refID (4B)</text>

      <rect x="60" y="0" width="55" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="87" y="12" text-anchor="middle" fill="#38bdf8" font-size="7">POS (4B)</text>

      <rect x="120" y="0" width="55" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="147" y="12" text-anchor="middle" fill="#38bdf8" font-size="7">MAPQ (1B)</text>

      <rect x="180" y="0" width="75" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
      <text x="217" y="12" text-anchor="middle" fill="#38bdf8" font-size="7">FLAG (2B)</text>

      <rect x="0" y="24" width="80" height="18" rx="2" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1"/>
      <text x="40" y="36" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">CIGAR (4B × k)</text>

      <rect x="85" y="24" width="80" height="18" rx="2" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
      <text x="125" y="36" text-anchor="middle" fill="#cbd5e1" font-size="7">SEQ (4-bit pac)</text>

      <rect x="170" y="24" width="85" height="18" rx="2" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
      <text x="212" y="36" text-anchor="middle" fill="#cbd5e1" font-size="7">QUAL (Phred)</text>

      <rect x="0" y="48" width="255" height="20" rx="2" fill="#059669" fill-opacity="0.25" stroke="#34d399" stroke-width="1"/>
      <text x="127" y="62" text-anchor="middle" fill="#34d399" font-size="8" font-weight="700">Optional Tags (CB, UB, GN, RE ... 가변 길이)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• BGZF (Blocked GNU Zip Format) 압축:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">64 KB 독립 압축 블록으로 BAI 인덱스 기반 O(1) 랜덤 점프</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 바이너리 정렬 압축 효율:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">SAM 텍스트 대비 용량 75% 절감 및 역직렬화 가속</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 0-based 게놈 좌표계 변환:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">SAM(1-based)에서 BAM 바이너리 저장 시 0-based로 인덱싱</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="150" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">인덱싱된 대용량 단일세포 전사체 표준 바이너리 포맷</text>
  </g>

  <!-- Panel 2: CIGAR String Operators -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12" font-weight="700">2. CIGAR 연산자 해부 (M, N, I, D, S)</text>

    <!-- CIGAR Schema -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="12" fill="#34d399" font-size="8" font-weight="700">M (Alignment Match):</text>
      <text x="105" y="12" fill="#cbd5e1" font-size="7.5">염기 일치 또는 단일 미스매치 (소비: Read/Ref)</text>

      <text x="0" y="28" fill="#fbbf24" font-size="8" font-weight="700">N (Skipped Region):</text>
      <text x="105" y="28" fill="#fbbf24" font-size="7.5" font-weight="700">스플라이스 인트론 건너뜀 (소비: Ref only)</text>

      <text x="0" y="44" fill="#38bdf8" font-size="8" font-weight="700">I / D (InDel):</text>
      <text x="105" y="44" fill="#cbd5e1" font-size="7.5">삽입(소비: Read) / 결실(소비: Ref)</text>

      <text x="0" y="60" fill="#f43f5e" font-size="8" font-weight="700">S (Soft Clipping):</text>
      <text x="105" y="60" fill="#cbd5e1" font-size="7.5">불일치 말단 서열 (소비: Read only)</text>

      <!-- Example CIGAR -->
      <line x1="0" y1="67" x2="250" y2="67" stroke="#334155" stroke-width="0.8"/>
      <text x="0" y="78" fill="#a7f3d0" font-size="7.5" font-weight="700">예시: 45M 12400N 55M ➔ 12.4kb 인트론을 건너뛴 스플라이싱</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• RNA-seq 스플라이싱의 핵심 연산자 N:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">게놈 상 수십만 bp 인트론을 건너뛰는 전사체 정렬 표현</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• Soft Clip(S)을 통한 융합 절단점 검출:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">유전자 융합(Fusion) 경계에서 비정상 정렬 분할점 식별</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 정렬 길이 계산 산술 규칙:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">End_POS = Start_POS + sum(M, D, N, = , X)</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#78350f" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#fbbf24" font-size="9" font-weight="600">정렬 갭과 엑손-엑손 스플라이싱 기하학을 정밀 압축 기술</text>
  </g>

  <!-- Panel 3: Single-Cell Specific Tags -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12" font-weight="700">3. 단일세포 전용 옵션 태그 (CB, UB, GN, RE)</text>

    <!-- Tags Detail Box -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <text x="0" y="12" fill="#34d399" font-size="8" font-weight="700">CB:Z:</text>
      <text x="45" y="12" fill="#cbd5e1" font-size="7.5">오류 교정된 세포 바코드 (Corrected Cell Barcode)</text>

      <text x="0" y="28" fill="#fbbf24" font-size="8" font-weight="700">UB:Z:</text>
      <text x="45" y="28" fill="#cbd5e1" font-size="7.5">오류 교정된 분자 식별자 (Corrected UMI)</text>

      <text x="0" y="44" fill="#38bdf8" font-size="8" font-weight="700">GN:Z:</text>
      <text x="45" y="44" fill="#cbd5e1" font-size="7.5">HUGO 공식 유전자 심볼 (Gene Symbol e.g. AR, ERG)</text>

      <text x="0" y="60" fill="#c084fc" font-size="8" font-weight="700">RE:A:</text>
      <text x="45" y="60" fill="#cbd5e1" font-size="7.5">영역 유형 (E: Exonic, I: Intronic, N: Intergenic)</text>

      <!-- Raw vs Corrected -->
      <line x1="0" y1="67" x2="250" y2="67" stroke="#334155" stroke-width="0.8"/>
      <text x="0" y="78" fill="#94a3b8" font-size="7.5">원시(Raw) 태그: CR:Z (Raw CB), UR:Z (Raw UMI)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 단일 리드에 다차원 메타데이터 바인딩:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">정렬 위치 + 세포 소속 + 분자 원점 + 유전자 주석 결합</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 세포별/유전자별 실시간 스트리밍 그룹화:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">BAM 전체를 스캔하며 CB/UB 기반 인메모리 행렬 축적</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 스플라이스 아이소폼 및 유전변이 마이닝:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">특정 세포 서브클론의 스플라이스 접합부 1-리드 단위 추출</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#065f46" fill-opacity="0.25"/>
    <text x="155" y="259" text-anchor="middle" fill="#34d399" font-size="9" font-weight="600">세포 수준의 발현 및 유전체 변이 추적을 가능케 하는 핵심 태그</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_ar_v7_splice_junction_and_tmprss2_erg_fusion(filepath: str):
    """
    6. AR-V7 Cryptic Exon 3b Splice Junction vs Canonical AR-FL & TMPRSS2-ERG Breakpoint
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
    <marker id="arrRose06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrEmerald06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16" font-weight="700">전립선암 단일세포 스플라이스 변이 및 염색체 융합: AR-V7 및 TMPRSS2-ERG</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">AR 잠재 엑손 3b(CE3b) 접합부 ➔ 엔잘루타마이드 저항성 서브클론 PSI(Ψ) ➔ 21q22 간질 결실 융합 유전자</text>

  <!-- Panel 1: AR-FL vs AR-V7 Splicing -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12" font-weight="700">1. AR 유전자좌: AR-FL vs AR-V7 변이체</text>

    <!-- Exon Junction Diagram -->
    <rect x="15" y="42" width="270" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Exon 3 -->
    <rect x="25" y="65" width="40" height="20" rx="3" fill="#0284c7" stroke="#38bdf8" stroke-width="1"/>
    <text x="45" y="78" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">Exon 3</text>

    <!-- Cryptic Exon 3b (AR-V7) -->
    <rect x="100" y="80" width="45" height="18" rx="3" fill="#be123c" stroke="#f43f5e" stroke-width="1.2"/>
    <text x="122" y="92" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">CE3b (V7)</text>

    <!-- Exon 4 (Canonical AR-FL) -->
    <rect x="195" y="65" width="40" height="20" rx="3" fill="#059669" stroke="#34d399" stroke-width="1"/>
    <text x="215" y="78" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">Exon 4</text>

    <!-- Exons 5-8 LBD -->
    <rect x="245" y="65" width="32" height="20" rx="3" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
    <text x="261" y="78" text-anchor="middle" fill="#94a3b8" font-size="7">Ex 5-8</text>

    <!-- Splicing lines -->
    <!-- AR-V7 Inclusion Junction -->
    <path d="M65,70 Q82,88 100,88" fill="none" stroke="#f43f5e" stroke-width="1.8"/>
    <text x="75" y="96" fill="#f43f5e" font-size="7" font-weight="700">2.3 kb (V7)</text>

    <!-- AR-FL Exclusion/Canonical Junction -->
    <path d="M65,65 Q130,42 195,65" fill="none" stroke="#34d399" stroke-width="1.8"/>
    <text x="130" y="52" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">11.9 kb (AR-FL 정규 접합)</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 잠재 엑손 3b(Cryptic Exon 3b) 조기 종결:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">인트론 3 내부 CE3b가 삽입되며 조기 정지 코돈(PTC) 발생</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 리간드 결합 도메인(LBD) 완전 결손:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">안드로겐 없이도 핵으로 이행해 하위 신호를 상시 활성화</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 엔잘루타마이드(Enzalutamide) 불응성 기전:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">호르몬 치료제 표적인 LBD 부재로 치료 저항성 획득</text>

    <!-- Badge -->
    <rect x="14" y="244" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.25"/>
    <text x="150" y="259" text-anchor="middle" fill="#f43f5e" font-size="9" font-weight="600">전립선암 표적 치료 불응의 결정적 바이오마커</text>
  </g>

  <!-- Panel 2: Single-Cell PSI Quantification -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="700">2. 단일세포 PSI (Percent Spliced-In, Ψ) 정량화</text>

    <!-- PSI Formula & Subclone Distribution -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- Formula -->
    <rect x="25" y="50" width="260" height="22" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="155" y="65" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">Ψ = I_CE3b / (I_CE3b + E_Ex4) [Inclusion Ratio]</text>

    <!-- Bimodal PSI curves -->
    <g transform="translate(25, 78)">
      <!-- Subclone Sensitive -->
      <rect x="10" y="8" width="105" height="18" rx="2" fill="#065f46" fill-opacity="0.3" stroke="#34d399" stroke-width="1"/>
      <text x="62" y="20" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">민감성 클론 (Ψ ≈ 0.02)</text>

      <!-- Subclone Resistant -->
      <rect x="135" y="8" width="115" height="18" rx="2" fill="#9f1239" fill-opacity="0.3" stroke="#f43f5e" stroke-width="1"/>
      <text x="192" y="20" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">저항성 서브클론 (Ψ ≈ 0.65)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 단일세포 스플라이스 정션 리드 마이닝:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">CIGAR 'N' 연산자의 게놈 좌표를 추출하여 I와 E 판정</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• 종양 이질성(Tumor Heterogeneity) 해명:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">벌크 분석에서 희석되던 소수 저항성 세포군을 명확히 식별</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• 최소 접합부 리드 수(Coverage ≥ 3) 필터:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">드롭아웃에 의한 PSI 과대 추정 방지 및 통계적 보정</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.2"/>
    <text x="155" y="259" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="600">세포별 아이소폼 분율의 이봉형(Bimodal) 클론 분리</text>
  </g>

  <!-- Panel 3: TMPRSS2-ERG Breakpoint -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12" font-weight="700">3. TMPRSS2-ERG 21q22 간질 결실 융합</text>

    <!-- Fusion Diagram -->
    <rect x="15" y="42" width="280" height="88" rx="6" fill="#0f172a" stroke="#334155" stroke-width="1"/>
    <!-- TMPRSS2 Exon 1 -->
    <rect x="25" y="52" width="85" height="18" rx="2" fill="#0284c7" stroke="#38bdf8" stroke-width="1"/>
    <text x="67" y="64" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">TMPRSS2 (Ex 1)</text>

    <!-- Deleted Interstitial region -->
    <line x1="110" y1="61" x2="165" y2="61" stroke="#f43f5e" stroke-width="1.5" stroke-dasharray="3,2"/>
    <text x="137" y="56" text-anchor="middle" fill="#f43f5e" font-size="6.5">~3 Mb 결실</text>

    <!-- ERG Exon 4 -->
    <rect x="165" y="52" width="85" height="18" rx="2" fill="#7e22ce" stroke="#c084fc" stroke-width="1"/>
    <text x="207" y="64" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">ERG (Ex 4)</text>

    <!-- Split Read Bar -->
    <rect x="45" y="80" width="180" height="16" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <rect x="45" y="80" width="85" height="16" rx="2" fill="#0284c7" fill-opacity="0.5"/>
    <text x="87" y="91" text-anchor="middle" fill="#38bdf8" font-size="7">50M (TMPRSS2)</text>
    <rect x="130" y="80" width="95" height="16" rx="2" fill="#7e22ce" fill-opacity="0.5"/>
    <text x="177" y="91" text-anchor="middle" fill="#c084fc" font-size="7">50S (ERG Chimeric)</text>

    <text x="155" y="116" text-anchor="middle" fill="#fbbf24" font-size="7.5">Breakpoint: chr21:41,465,100 ^ chr21:39,755,200</text>

    <!-- Bullets -->
    <text x="14" y="148" fill="#cbd5e1" font-size="10" font-weight="600">• 안드로겐 반응 프로모터에 의한 발암 구동:</text>
    <text x="24" y="163" fill="#94a3b8" font-size="9">TMPRSS2 전사 조절을 받아 ERG 전사인자가 과발현</text>
    <text x="14" y="181" fill="#cbd5e1" font-size="10" font-weight="600">• Chimeric Split Read &amp; Discordant Pairs:</text>
    <text x="24" y="196" fill="#94a3b8" font-size="9">CIGAR 50M50S 키메라 연산자를 통해 융합 절단점 검출</text>
    <text x="14" y="214" fill="#cbd5e1" font-size="10" font-weight="600">• AR-V7 발현 서브클론과의 공존성 검증:</text>
    <text x="24" y="229" fill="#94a3b8" font-size="9">동일 세포 바코드(CB) 내에서 융합 유전자와 변이체 동시 검출</text>

    <!-- Badge -->
    <rect x="14" y="244" width="282" height="22" rx="4" fill="#581c87" fill-opacity="0.3"/>
    <text x="155" y="259" text-anchor="middle" fill="#c084fc" font-size="9" font-weight="600">단일세포 분해능으로 암 유전자 재배열과 스플라이싱 통합</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_fastq_read_structure_and_barcode_correction.svg", make_01_fastq_read_structure_and_barcode_correction),
        ("02_starsolo_suffix_array_spliced_alignment.svg", make_02_starsolo_suffix_array_spliced_alignment),
        ("03_alevin_fry_rad_file_and_splici_index.svg", make_03_alevin_fry_rad_file_and_splici_index),
        ("04_umi_directional_adjacency_graph_deduplication.svg", make_04_umi_directional_adjacency_graph_deduplication),
        ("05_bam_binary_record_cigar_and_tags.svg", make_05_bam_binary_record_cigar_and_tags),
        ("06_ar_v7_splice_junction_and_tmprss2_erg_fusion.svg", make_06_ar_v7_splice_junction_and_tmprss2_erg_fusion),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 09 Upstream Alignment & BAM publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

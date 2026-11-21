#!/usr/bin/env python3
"""
tracks/01_multiome_snatac/scripts/generate_multiome_diagrams.py
==============================================================
Generates publication-grade, high-resolution SVG diagrams for Track 01:
Multiome (snRNA + snATAC) & Epigenomic Regulatory Networks in Prostate Cancer.

Covers:
1. Tn5 Tagmentation & Single-Nucleus Dual Capture Mechanism
2. 5-Column Fragment File Format & TSS Enrichment QC
3. LSI (TF-IDF + TruncatedSVD) Dimensionality Reduction Geometry
4. WNN (Weighted Nearest Neighbor) Multi-Modal Integration
5. chrX:66.9 Mb AR Enhancer Amplification & Peak-to-Gene Looping
6. SCENIC+ Triplet Architecture (TF -> Enhancer Motif -> Target Gene eRegulon)
"""

import os
import xml.etree.ElementTree as ET


def make_01_tn5_tagmentation(filepath: str):
    """
    1. Tn5 Tagmentation Mechanism
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard1" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="gradTn5" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#818cf8"/>
      <stop offset="100%" stop-color="#c084fc"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">10x Multiome: Tn5 트랜스포사제 태그멘테이션 및 단일 핵 동시 포획 메커니즘</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">Tn5 Homodimer Cut-and-Paste Transposition (9-bp Staggered Cut) ➔ 온전한 핵막 내 확산 ➔ Dual-barcoded Gel Bead 캡슐화</text>

  <!-- Panel 1: Intact Nuclear Isolation -->
  <g transform="translate(30, 70)">
    <rect width="295" height="265" rx="10" fill="url(#gradCard1)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="147" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 초저온 핵 추출 &amp; 온전한 핵막 보존</text>

    <!-- Schematic: Nuclear envelope with pore -->
    <rect x="15" y="44" width="265" height="52" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <circle cx="50" cy="70" r="18" fill="#0284c7" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="50" y="74" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">Nucleus</text>
    <path d="M 75 70 L 115 70" stroke="#38bdf8" stroke-width="2" stroke-dasharray="3,2"/>
    <text x="125" y="65" fill="#34d399" font-size="10" font-weight="700">Lamin B 핵막 보존 (0.1% Lysis)</text>
    <text x="125" y="80" fill="#94a3b8" font-size="9.5">핵공 복합체(NPC) 통한 효소 침투</text>

    <!-- Technical Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• 세포질 선택적 용해 (Cold Lysis):</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">세포막만 파괴하여 미토콘드리아 RNA 배제</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 진정염색질(Euchromatin) 접근성:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">느슨하게 열린 조절 부위만 효소 접촉 허용</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 이질염색질(Heterochromatin) 차단:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">고밀도 응축 영역의 비특이적 절단 방지</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="265" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="147" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">✓ 순수 단일 핵 무결성 100% 유지</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 350 200" stroke="#818cf8" stroke-width="2.5"/>
  <polygon points="350,195 360,200 350,205" fill="#818cf8"/>

  <!-- Panel 2: Tn5 Transposition Mechanics -->
  <g transform="translate(365, 70)">
    <rect width="300" height="265" rx="10" fill="url(#gradCard1)" stroke="#818cf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#818cf8" fill-opacity="0.2"/>
    <text x="150" y="22" text-anchor="middle" fill="#818cf8" font-size="13" font-weight="700">2. Tn5 태그멘테이션 (9-bp Staggered Cut)</text>

    <!-- Schematic: Tn5 Transposase Homodimer with Adapters -->
    <rect x="15" y="44" width="270" height="52" rx="6" fill="#0f172a" stroke="#4338ca" stroke-width="1"/>
    <ellipse cx="60" cy="70" rx="20" ry="14" fill="url(#gradTn5)" opacity="0.85"/>
    <text x="60" y="74" text-anchor="middle" fill="#fff" font-size="9" font-weight="700">Tn5-Dimer</text>
    <line x1="85" y1="64" x2="135" y2="64" stroke="#f43f5e" stroke-width="2"/>
    <line x1="94" y1="76" x2="144" y2="76" stroke="#38bdf8" stroke-width="2"/>
    <text x="155" y="66" fill="#f43f5e" font-size="9.5" font-weight="700">Read 1N</text>
    <text x="155" y="80" fill="#38bdf8" font-size="9.5" font-weight="700">Read 2N (9-bp gap)</text>

    <!-- Technical Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• Cut-and-Paste Transposition:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">염색질 절단과 동시에 19-bp ME 어댑터 결찰</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 9-bp 지그재그 돌출 말단 (Staggered Cut):</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">양 가닥에서 9-bp 편차 절단 (분석 시 +4/-5bp 보정)</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 뉴클레오솜 간격 태깅 (Banding):</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">히스톤 옥타머 회피 ➔ 147-bp 주기적 단편 생성</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="270" height="22" rx="4" fill="#4338ca" fill-opacity="0.3"/>
    <text x="150" y="250" text-anchor="middle" fill="#818cf8" font-size="10" font-weight="600">✓ 전유전체 오픈 크로마틴 즉각적 어댑터 태깅</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 670 200 L 690 200" stroke="#34d399" stroke-width="2.5"/>
  <polygon points="690,195 700,200 690,205" fill="#34d399"/>

  <!-- Panel 3: Dual-Barcoded Gel Bead Encapsulation -->
  <g transform="translate(705, 70)">
    <rect width="265" height="265" rx="10" fill="url(#gradCard1)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="265" height="32" rx="10" fill="#34d399" fill-opacity="0.2"/>
    <text x="132" y="22" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 10x Next GEM 이중 바코딩</text>

    <!-- Schematic: Gel bead with dual oligos -->
    <rect x="15" y="44" width="235" height="52" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <circle cx="45" cy="70" r="16" fill="#10b981" fill-opacity="0.3" stroke="#34d399" stroke-width="1.5"/>
    <text x="45" y="74" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">GelBead</text>
    <text x="75" y="63" fill="#38bdf8" font-size="9.5" font-weight="600">Oligo A: Poly(dT) ➔ snRNA</text>
    <text x="75" y="80" fill="#c084fc" font-size="9.5" font-weight="600">Oligo B: Spacer ➔ snATAC</text>

    <!-- Technical Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• 동일 핵(Same Cell) 1:1 결합:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">단일 유세포 액적 내 하나의 핵만 격리</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 공통 16-bp 10x 바코드 공유:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">RNA UMI와 ATAC 단편이 동일 식별자 획득</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 전사체와 후성유전체의 완벽 동기화:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">컴퓨팅 예측이 아닌 '실제 물리적 동시 계측'</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="235" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="132" y="250" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">✓ 진정한 멀티오믹스 단일 핵 융합</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_fragment_and_tss(filepath: str):
    """
    2. Fragment File Structure & TSS Enrichment QC
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard2" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">snATAC-seq 데이터 구조: 5-컬럼 Fragment 파일 &amp; TSS 농축도(Enrichment) 검증</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">Tabix 인덱싱 고속 게놈 청킹 ➔ 전사시작점(TSS ±1kb) 8× 이상 피크 프로파일 vs 배경 잡음 QC</text>

  <!-- Panel 1: Fragment File Format & Tabix -->
  <g transform="translate(30, 70)">
    <rect width="455" height="265" rx="10" fill="url(#gradCard2)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="455" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="227" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">A. 5-컬럼 단일 단편 파일 (fragments.tsv.gz) &amp; Tabix</text>

    <!-- Code Block / Fragment Table -->
    <rect x="15" y="44" width="425" height="85" rx="6" fill="#020617" stroke="#334155" stroke-width="1"/>
    <text x="25" y="62" fill="#94a3b8" font-size="10" font-family="monospace">#chrom   chromStart  chromEnd    barcode             dup_count</text>
    <text x="25" y="79" fill="#38bdf8" font-size="10.5" font-family="monospace">chrX     66901240    66901450    AAACAGCCAAGCT-1     3</text>
    <text x="25" y="96" fill="#818cf8" font-size="10.5" font-family="monospace">chr8     127735100   127735310   AAACAGCCAAGCT-1     1</text>
    <text x="25" y="113" fill="#34d399" font-size="10.5" font-family="monospace">chr10    87863200    87863410    TTTACTGAGCCTA-1     5</text>

    <!-- Key Rules Bullets -->
    <text x="15" y="150" fill="#cbd5e1" font-size="11" font-weight="600">• 0-based 반개방 좌표 &amp; Tn5 9-bp 절단점 보정:</text>
    <text x="25" y="167" fill="#94a3b8" font-size="10.5">Tn5 중심 이동: +가닥은 +4bp, -가닥은 -5bp 정밀 시프트</text>

    <text x="15" y="189" fill="#cbd5e1" font-size="11" font-weight="600">• barcode: snRNA-seq 16-bp 바코드와 100% 일치 매칭:</text>
    <text x="25" y="206" fill="#94a3b8" font-size="10.5">단일 핵 식별자 일치로 전사체와 염색질 매트릭스 직접 결합</text>

    <text x="15" y="228" fill="#cbd5e1" font-size="11" font-weight="600">• Tabix BGI 청킹 인덱스 (`tabix -p bed`):</text>
    <text x="25" y="245" fill="#94a3b8" font-size="10.5">수십 GB 단편 파일에서 임의 유전체 윈도우를 O(log N) 탐색</text>

    <!-- Bottom summary badge -->
    <rect x="15" y="260" width="425" height="2" fill="#38bdf8" fill-opacity="0.3"/>
  </g>

  <!-- Panel 2: TSS Enrichment Profile & Nucleosome Banding -->
  <g transform="translate(515, 70)">
    <rect width="455" height="265" rx="10" fill="url(#gradCard2)" stroke="#10b981" stroke-width="1.5"/>
    <rect x="0" y="0" width="455" height="32" rx="10" fill="#10b981" fill-opacity="0.2"/>
    <text x="227" y="22" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">B. 전사시작점(TSS) 농축도 및 뉴클레오솜 밴딩 QC</text>

    <!-- Graph: TSS Profile Curve -->
    <g transform="translate(25, 45)">
      <rect width="210" height="90" rx="5" fill="#0f172a" stroke="#1e293b"/>
      <!-- Axes -->
      <line x1="20" y1="80" x2="195" y2="80" stroke="#475569" stroke-width="1"/>
      <line x1="107" y1="15" x2="107" y2="80" stroke="#475569" stroke-dasharray="2,2"/>
      <text x="107" y="88" text-anchor="middle" fill="#94a3b8" font-size="8.5">TSS</text>
      <text x="35" y="88" text-anchor="middle" fill="#64748b" font-size="8">-2kb</text>
      <text x="180" y="88" text-anchor="middle" fill="#64748b" font-size="8">+2kb</text>
      <!-- TSS Sharp Peak -->
      <path d="M 25 78 Q 80 77 95 60 Q 107 20 107 20 Q 107 20 119 60 Q 134 77 190 78" fill="none" stroke="#10b981" stroke-width="2.2"/>
      <circle cx="107" cy="20" r="3.5" fill="#34d399"/>
      <text x="115" y="24" fill="#34d399" font-size="9" font-weight="700">&gt; 8× Signal</text>
      <line x1="20" y1="76" x2="195" y2="76" stroke="#f43f5e" stroke-dasharray="3,2"/>
      <text x="150" y="73" fill="#f87171" font-size="7.5">1× Background</text>
    </g>

    <!-- Inset: Nucleosome Banding Info -->
    <g transform="translate(250, 45)">
      <rect width="190" height="90" rx="5" fill="#0f172a" stroke="#047857" stroke-width="1"/>
      <text x="15" y="20" fill="#34d399" font-size="10.5" font-weight="700">뉴클레오솜 래더 주기성</text>
      <text x="15" y="38" fill="#cbd5e1" font-size="9.5">• &lt; 100 bp : NFR (무핵소체 단편)</text>
      <text x="15" y="55" fill="#cbd5e1" font-size="9.5">• ~200 bp : Mononucleosome</text>
      <text x="15" y="72" fill="#cbd5e1" font-size="9.5">• ~400 bp : Dinucleosome</text>
    </g>

    <!-- Technical QC Standards -->
    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• TSS 농축도 산출 공식:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">TSS 중심(±100bp) 평균 절단 수 / 외곽 플랭크(±1.5~2kb) 배경</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 엄격한 품질 필터 임계치 (QC Threshold):</text>
    <text x="25" y="214" fill="#34d399" font-size="10.5">TSS Score &gt; 8.0 통과 (5.0 미만은 핵 용해 및 비특이 파편 탈락)</text>

    <text x="15" y="236" fill="#cbd5e1" font-size="11" font-weight="600">• 최소 고유 단편 수 (Unique Fragments):</text>
    <text x="25" y="253" fill="#94a3b8" font-size="10.5">핵당 최소 1,000 ~ 3,000개 이상의 유효 매핑 단편 요구</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_lsi_tfidf(filepath: str):
    """
    3. LSI Dimensionality Reduction Geometry
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard3" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">snATAC-seq 차원 축소: 잠재 의미 지수화(LSI = TF-IDF ➔ TruncatedSVD)</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">초고희소성(Sparsity &gt;98%) 이진화 피크 행렬 ➔ 정보 가중 TF-IDF ➔ SVD 잠재 다양체 및 Component 1 심도 편향 진단</text>

  <!-- Panel 1: Binary Peak Matrix -->
  <g transform="translate(30, 70)">
    <rect width="290" height="265" rx="10" fill="url(#gradCard3)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="290" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="145" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 이진 피크 행렬 (Binarization)</text>

    <!-- Schematic: Binary Matrix -->
    <rect x="15" y="44" width="260" height="52" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <text x="25" y="64" fill="#38bdf8" font-size="10.5" font-family="monospace">Peak Matrix: A_ij ∈ {0, 1}</text>
    <text x="25" y="82" fill="#94a3b8" font-size="9.5">1 = Open Chromatin, 0 = Inaccessible</text>

    <!-- Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• 2배체(Diploid) 생물학적 한계:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">핵당 대립유전자는 최대 2카피 (카운트 왜곡 제거)</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 초고희소성 (Sparsity &gt; 98%):</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">10만 개 피크 중 대부분의 세포에서 0값 기록</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• PCR 중복 바이어스 배제:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">이진화로 PCR 오버앰플리피케이션 잡음 차단</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="260" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="145" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">50,000 ~ 150,000 Peaks × Nuclei</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 325 200 L 345 200" stroke="#38bdf8" stroke-width="2.5"/>
  <polygon points="345,195 355,200 345,205" fill="#38bdf8"/>

  <!-- Panel 2: TF-IDF Normalization -->
  <g transform="translate(360, 70)">
    <rect width="320" height="265" rx="10" fill="url(#gradCard3)" stroke="#818cf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#818cf8" fill-opacity="0.2"/>
    <text x="160" y="22" text-anchor="middle" fill="#818cf8" font-size="13" font-weight="700">2. TF-IDF 빈도-역문서 가중치</text>

    <!-- Formulas -->
    <rect x="15" y="44" width="290" height="60" rx="6" fill="#0f172a" stroke="#4338ca" stroke-width="1"/>
    <text x="25" y="62" fill="#818cf8" font-size="10" font-family="monospace">TF_ij = A_ij / Σ_k A_ik (세포 점유율)</text>
    <text x="25" y="78" fill="#34d399" font-size="10" font-family="monospace">IDF_j = log(1 + N_cells / Σ_i A_ij)</text>
    <text x="25" y="94" fill="#fbbf24" font-size="10" font-family="monospace">TF-IDF_ij = TF_ij × IDF_j</text>

    <!-- Bullets -->
    <text x="15" y="124" fill="#cbd5e1" font-size="11" font-weight="600">• 자연어 처리(NLP) 알고리즘 응용:</text>
    <text x="25" y="142" fill="#94a3b8" font-size="10.5">문서(세포) 내 희귀 단어(피크)의 정보 가치 증폭</text>

    <text x="15" y="164" fill="#cbd5e1" font-size="11" font-weight="600">• 유비쿼터스 하우스키핑 피크 감쇠:</text>
    <text x="25" y="182" fill="#94a3b8" font-size="10.5">모든 세포에서 열린 보편 피크는 IDF 패널티</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 희귀 암 아형 인핸서 정보 극대화:</text>
    <text x="25" y="222" fill="#94a3b8" font-size="10.5">특정 전립선암 클론만의 조절 피크 가중치 상승</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="290" height="22" rx="4" fill="#4338ca" fill-opacity="0.3"/>
    <text x="160" y="250" text-anchor="middle" fill="#818cf8" font-size="10" font-weight="600">✓ 정보 엔트로피 기반 정규화 완료</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 685 200 L 705 200" stroke="#818cf8" stroke-width="2.5"/>
  <polygon points="705,195 715,200 705,205" fill="#818cf8"/>

  <!-- Panel 3: TruncatedSVD & LSI Component 1 Warning -->
  <g transform="translate(720, 70)">
    <rect width="250" height="265" rx="10" fill="url(#gradCard3)" stroke="#ef4444" stroke-width="1.5"/>
    <rect x="0" y="0" width="250" height="32" rx="10" fill="#ef4444" fill-opacity="0.2"/>
    <text x="125" y="22" text-anchor="middle" fill="#f87171" font-size="13" font-weight="700">3. TruncatedSVD &amp; LSI-1 경고</text>

    <!-- Warning Box -->
    <rect x="12" y="44" width="226" height="60" rx="6" fill="#450a0a" stroke="#ef4444" stroke-width="1"/>
    <text x="20" y="62" fill="#fca5a5" font-size="10" font-weight="700">⚠️ LSI Component 1 시퀀싱 심도 편향</text>
    <text x="20" y="78" fill="#f87171" font-size="9.5" font-family="monospace">r(LSI_1, Total_Peaks) &gt; 0.90</text>
    <text x="20" y="93" fill="#cbd5e1" font-size="9">기술적 UMI 수와 직결된 가짜 축!</text>

    <!-- Bullets -->
    <text x="12" y="124" fill="#cbd5e1" font-size="11" font-weight="600">• LSI Component 1 진단 필수:</text>
    <text x="20" y="142" fill="#fca5a5" font-size="10">피크 총수와의 피어슨 상관계수 r 측정</text>

    <text x="12" y="164" fill="#cbd5e1" font-size="11" font-weight="600">• 표준 전처리 솔루션:</text>
    <text x="20" y="182" fill="#34d399" font-size="10">r &gt; 0.90 시 Component 1 과감히 배제 (Drop)</text>
    <text x="20" y="197" fill="#34d399" font-size="10">또는 회귀 모델(Linear Regression)로 보정</text>

    <text x="12" y="219" fill="#cbd5e1" font-size="11" font-weight="600">• 2~30번 컴포넌트 활용:</text>
    <text x="20" y="237" fill="#94a3b8" font-size="10">진정한 생물학적 후성유전 다양체 확보</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_wnn_weighting(filepath: str):
    """
    4. WNN Multi-Modal Integration
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard4" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">WNN (Weighted Nearest Neighbor) 이중 오믹스 통합 기하학</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">단일 세포별 전사체(snRNA) 및 염색질(snATAC) 정보량 동적 가중치 학습 (w_RNA + w_ATAC = 1.0) ➔ 통합 kNN 매니폴드</text>

  <!-- Panel 1: Modality Graphs -->
  <g transform="translate(30, 70)">
    <rect width="280" height="265" rx="10" fill="url(#gradCard4)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="280" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="140" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 모달리티별 국소 kNN 다양체</text>

    <!-- Schematic: Two graphs -->
    <rect x="15" y="44" width="250" height="52" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <circle cx="50" cy="70" r="14" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8"/>
    <text x="50" y="74" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">RNA</text>
    <path d="M 68 70 L 110 70" stroke="#38bdf8" stroke-width="1.5"/>
    <circle cx="128" cy="70" r="14" fill="#4338ca" fill-opacity="0.4" stroke="#818cf8"/>
    <text x="128" y="74" text-anchor="middle" fill="#818cf8" font-size="9" font-weight="700">ATAC</text>
    <text x="155" y="65" fill="#38bdf8" font-size="9.5">PCA (30 PCs)</text>
    <text x="155" y="79" fill="#818cf8" font-size="9.5">LSI (30 Dim)</text>

    <!-- Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• snRNA kNN 다양체:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">2,500 HVG 발현량 기반 표현형 유사도</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• snATAC kNN 다양체:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">유전체 오픈 크로마틴 후성유전 상태</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 모달리티의 한계:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">RNA나 ATAC 단독으로는 미세 클론 중첩 발생</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="250" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="140" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">독립적 2개 그래프 동시 구성</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 315 200 L 335 200" stroke="#818cf8" stroke-width="2.5"/>
  <polygon points="335,195 345,200 335,205" fill="#818cf8"/>

  <!-- Panel 2: Dynamic Affinity Weight Learning -->
  <g transform="translate(350, 70)">
    <rect width="335" height="265" rx="10" fill="url(#gradCard4)" stroke="#818cf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="335" height="32" rx="10" fill="#818cf8" fill-opacity="0.2"/>
    <text x="167" y="22" text-anchor="middle" fill="#818cf8" font-size="13" font-weight="700">2. 세포별 신호-잡음비 기반 가중치 학습</text>

    <!-- Formula Box -->
    <rect x="15" y="44" width="305" height="60" rx="6" fill="#0f172a" stroke="#4338ca" stroke-width="1"/>
    <text x="25" y="62" fill="#818cf8" font-size="10" font-family="monospace">x̂_i,RNA = Σ_j sim(i,j) · x_j (모달 예측)</text>
    <text x="25" y="78" fill="#34d399" font-size="10" font-family="monospace">θ_i,RNA = ||x_i - x̂_i,RNA||⁻¹ (예측 친화도)</text>
    <text x="25" y="94" fill="#fbbf24" font-size="10" font-family="monospace">w_RNA + w_ATAC = 1.0 (Softmax 가중)</text>

    <!-- Bullets -->
    <text x="15" y="124" fill="#cbd5e1" font-size="11" font-weight="600">• 교차 모달리티 예측 정확도 평가:</text>
    <text x="25" y="142" fill="#94a3b8" font-size="10.5">RNA 이웃으로 ATAC을 얼마나 잘 예측하는가</text>

    <text x="15" y="164" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 세포 수준 적응적 가중치 (Adaptive):</text>
    <text x="25" y="182" fill="#94a3b8" font-size="10.5">전체 세포 일괄 가중치가 아닌, 세포별 가중치 부여</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 저품질 모달리티 자동 감쇠:</text>
    <text x="25" y="222" fill="#94a3b8" font-size="10.5">드롭아웃 심한 축은 스스로 가중치 축소</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="305" height="22" rx="4" fill="#4338ca" fill-opacity="0.3"/>
    <text x="167" y="250" text-anchor="middle" fill="#818cf8" font-size="10" font-weight="600">✓ 정보 엔트로피 최적화 가중 결합</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 690 200 L 710 200" stroke="#34d399" stroke-width="2.5"/>
  <polygon points="710,195 720,200 710,205" fill="#34d399"/>

  <!-- Panel 3: Joint Manifold & Subtle Clones -->
  <g transform="translate(725, 70)">
    <rect width="245" height="265" rx="10" fill="url(#gradCard4)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="245" height="32" rx="10" fill="#34d399" fill-opacity="0.2"/>
    <text x="122" y="22" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 통합 WNN 다양체 &amp; 아형 분리</text>

    <!-- Graphic: Clusters -->
    <rect x="12" y="44" width="220" height="52" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <circle cx="40" cy="70" r="10" fill="#38bdf8"/>
    <text x="40" y="73" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">AR+</text>
    <circle cx="85" cy="65" r="10" fill="#f59e0b"/>
    <text x="85" y="68" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">Amphi</text>
    <circle cx="135" cy="72" r="10" fill="#818cf8"/>
    <text x="135" y="75" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">DNPC</text>
    <circle cx="185" cy="68" r="10" fill="#f43f5e"/>
    <text x="185" y="71" text-anchor="middle" fill="#fff" font-size="7" font-weight="700">t-NEPC</text>

    <!-- Bullets -->
    <text x="12" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• Amphicrine 상태 완벽 분리:</text>
    <text x="20" y="134" fill="#34d399" font-size="10">RNA는 AR 양성이나 ATAC은 신경계 개방</text>

    <text x="12" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 선행적 후성 재편 포착:</text>
    <text x="20" y="174" fill="#94a3b8" font-size="10">전사 발현 이전 염색질 개방 선행 감지</text>

    <text x="12" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 2D WNN UMAP 투영:</text>
    <text x="20" y="214" fill="#94a3b8" font-size="10">높은 해상도로 치료 저항성 궤적 가시화</text>

    <!-- Bottom Badge -->
    <rect x="12" y="235" width="220" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="122" y="250" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">✓ 임상 미세 아형 분별력 극대화</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_ar_enhancer_duplication(filepath: str):
    """
    5. AR Enhancer Amplification & Peak-to-Gene Looping
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard5" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">전립선암 핵심 후성유전 기전: chrX:66.9 Mb AR 인핸서 증폭 및 Peak-to-Gene 연계</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">원발암 정상 1-copy vs mCRPC 650 kb 원거리 인핸서 탠덤 중복(Tandem Duplication) ➔ AR 프로모터 루핑 &amp; 발현 폭발</text>

  <!-- Top Genomic Track Architecture -->
  <g transform="translate(40, 68)">
    <rect width="920" height="120" rx="8" fill="#020617" stroke="#334155" stroke-width="1.2"/>
    <text x="20" y="24" fill="#94a3b8" font-size="11" font-weight="700">Locus: chrX: 66,000,000 ~ 67,200,000 (GRCh38, 1.2 Mb Genomic Window)</text>

    <!-- DNA Backbone -->
    <line x1="40" y1="65" x2="880" y2="65" stroke="#475569" stroke-width="3"/>

    <!-- AR Gene Promoter -->
    <rect x="750" y="48" width="100" height="34" rx="5" fill="#0284c7" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="800" y="69" text-anchor="middle" fill="#fff" font-size="12" font-weight="700">AR Promoter</text>

    <!-- Normal Enhancer (1-copy) -->
    <rect x="180" y="50" width="80" height="30" rx="5" fill="#4338ca" stroke="#818cf8" stroke-width="1"/>
    <text x="220" y="69" text-anchor="middle" fill="#818cf8" font-size="10" font-weight="700">Enhancer (1×)</text>

    <!-- mCRPC Duplication Region -->
    <rect x="270" y="50" width="120" height="30" rx="5" fill="#b91c1c" fill-opacity="0.3" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,2"/>
    <text x="330" y="69" text-anchor="middle" fill="#f87171" font-size="10" font-weight="700">Tandem Duplication (2~5×)</text>

    <!-- 3D Chromatin Looping Arc -->
    <path d="M 220 50 Q 510 -15 800 48" fill="none" stroke="#f59e0b" stroke-width="2.5" stroke-dasharray="6,3"/>
    <text x="510" y="20" text-anchor="middle" fill="#fbbf24" font-size="11" font-weight="700">650-kb 원거리 크로마틴 루핑 (CTCF / Cohesin Loop)</text>

    <text x="220" y="100" text-anchor="middle" fill="#94a3b8" font-size="9.5">-650 kb Distal Enhancer</text>
    <text x="800" y="100" text-anchor="middle" fill="#38bdf8" font-size="9.5">TSS (전사시작점)</text>
  </g>

  <!-- Bottom Left: Peak-to-Gene Statistics -->
  <g transform="translate(40, 200)">
    <rect width="445" height="135" rx="8" fill="url(#gradCard5)" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="20" y="26" fill="#38bdf8" font-size="12.5" font-weight="700">Peak-to-Gene 거리 감쇠 상관 연계 (Distance-Decayed Correlation)</text>

    <rect x="15" y="38" width="415" height="36" rx="5" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <text x="25" y="54" fill="#38bdf8" font-size="9.5" font-family="monospace">Reg_Score = r(Peak_ATAC, Gene_RNA) × exp(-d / 100kb)</text>
    <text x="25" y="67" fill="#34d399" font-size="9.5" font-family="monospace">유의성 기준: |r| &gt; 0.35, p &lt; 0.01 (±250 kb Window)</text>

    <text x="15" y="93" fill="#cbd5e1" font-size="10.5">• 인핸서 개방도와 AR pre-mRNA 전사량의 완벽한 비례 관계 (r = 0.68***)</text>
    <text x="15" y="110" fill="#94a3b8" font-size="10">• snATAC 피크 증가가 실시간으로 AR RNA 발현 폭발로 이어짐을 통계적 입증</text>
  </g>

  <!-- Bottom Right: Clinical Relevance & Resistance -->
  <g transform="translate(515, 200)">
    <rect width="445" height="135" rx="8" fill="url(#gradCard5)" stroke="#f59e0b" stroke-width="1.5"/>
    <text x="20" y="26" fill="#fbbf24" font-size="12.5" font-weight="700">호르몬 요법(Enzalutamide) 저항성 분자 기전</text>

    <text x="15" y="52" fill="#cbd5e1" font-size="10.5">• 비암호화 유전체 구조 변이의 전사 활성화:</text>
    <text x="25" y="68" fill="#94a3b8" font-size="10">코딩 돌연변이 없이 650kb 상류 인핸서 탠덤 증폭만으로 AR 재활성화</text>

    <text x="15" y="90" fill="#cbd5e1" font-size="10.5">• mCRPC 환자 코호트의 70~80%에서 발견되는 범용 저항성 표적:</text>
    <text x="25" y="106" fill="#f87171" font-size="10" font-weight="600">BET 브로모도메인 저해제(JQ1) 및 p300/CBP 억제제를 통한 루핑 차단 필요</text>

    <!-- Bottom summary line -->
    <rect x="15" y="120" width="415" height="2" fill="#f59e0b" fill-opacity="0.3"/>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_scenicplus_triplet(filepath: str):
    """
    6. SCENIC+ Triplet Architecture (TF -> Enhancer Motif -> Target Gene)
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard6" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">SCENIC+ eGRN 아키텍처: TF ➔ 인핸서 모티프 ➔ 표적 유전자 3중체(Triplet) 네트워크</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">단일 핵 동시 발현-접근성 기반 전사 인자-인핸서-표적 유전자 인과 조절 모듈(eRegulon) 구축</text>

  <!-- Node 1: Master TF -->
  <g transform="translate(30, 70)">
    <rect width="280" height="265" rx="10" fill="url(#gradCard6)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="280" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="140" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 마스터 전사 인자 (TF)</text>

    <!-- Schematic: TF node -->
    <rect x="15" y="44" width="250" height="52" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <circle cx="50" cy="70" r="16" fill="#0284c7" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="50" y="74" text-anchor="middle" fill="#fff" font-size="10" font-weight="700">TF</text>
    <text x="78" y="64" fill="#38bdf8" font-size="10.5" font-weight="700">FOXA1 / AR / ASCL1</text>
    <text x="78" y="79" fill="#94a3b8" font-size="9.5">전사체 발현량 (GEX UMI)</text>

    <!-- Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• FOXA1 (Pioneer Factor):</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">닫힌 염색질을 최초로 열어내는 개척자</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• AR (Nuclear Receptor):</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">호르몬 의존적 종양 분화 및 증식 지휘</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• ASCL1 (bHLH Factor):</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">t-NEPC 신경내분비 계통 전환 마스터 TF</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="250" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="140" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">Regulator Expression Score</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 315 200 L 335 200" stroke="#818cf8" stroke-width="2.5"/>
  <polygon points="335,195 345,200 335,205" fill="#818cf8"/>

  <!-- Node 2: Enhancer Motif & ChromVAR -->
  <g transform="translate(350, 70)">
    <rect width="320" height="265" rx="10" fill="url(#gradCard6)" stroke="#818cf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#818cf8" fill-opacity="0.2"/>
    <text x="160" y="22" text-anchor="middle" fill="#818cf8" font-size="13" font-weight="700">2. 오픈 인핸서 &amp; DNA 바인딩 모티프</text>

    <!-- Schematic: DNA Motif -->
    <rect x="15" y="44" width="290" height="52" rx="6" fill="#0f172a" stroke="#4338ca" stroke-width="1"/>
    <text x="25" y="63" fill="#818cf8" font-size="10.5" font-family="monospace">Motif: AGAACAnnnTGTTCT (ARE)</text>
    <text x="25" y="80" fill="#34d399" font-size="9.5" font-family="monospace">ChromVAR Z-score: Motif Accessibility</text>

    <!-- Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• cisTopic 토픽 모델링:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">공 개방 인핸서 피크 공발현 클러스터 도출</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• pycistarget 모티프 인리치먼트:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">열린 피크 중심 실제 TF 결합 서열 물리적 검증</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• Peak-to-Gene 링크 연계:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">인핸서와 프로모터 간 기능적 연결 검증</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="290" height="22" rx="4" fill="#4338ca" fill-opacity="0.3"/>
    <text x="160" y="250" text-anchor="middle" fill="#818cf8" font-size="10" font-weight="600">✓ 물리적 결합 타당성 확보</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 675 200 L 695 200" stroke="#34d399" stroke-width="2.5"/>
  <polygon points="695,195 705,200 695,205" fill="#34d399"/>

  <!-- Node 3: Target Gene & eRegulon -->
  <g transform="translate(710, 70)">
    <rect width="260" height="265" rx="10" fill="url(#gradCard6)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="260" height="32" rx="10" fill="#34d399" fill-opacity="0.2"/>
    <text x="130" y="22" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 표적 유전자 발현 &amp; eRegulon</text>

    <!-- Formula / eRegulon Box -->
    <rect x="15" y="44" width="230" height="52" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <text x="25" y="63" fill="#34d399" font-size="10" font-weight="700">eRegulon = TF ➔ Peak ➔ Target</text>
    <text x="25" y="80" fill="#94a3b8" font-size="9">KLK3, TMPRSS2, SYP, CHGA</text>

    <!-- Bullets -->
    <text x="15" y="116" fill="#cbd5e1" font-size="11" font-weight="600">• 단순 상관관계를 넘은 인과성:</text>
    <text x="25" y="134" fill="#94a3b8" font-size="10.5">TF발현 + 인핸서개방 + 표적전사 3박자 일치</text>

    <text x="15" y="156" fill="#cbd5e1" font-size="11" font-weight="600">• 활성화(+) vs 억제(-) 분기:</text>
    <text x="25" y="174" fill="#94a3b8" font-size="10.5">REST 억제 해제 및 AR 활성 감쇠 규명</text>

    <text x="15" y="196" fill="#cbd5e1" font-size="11" font-weight="600">• 표적 발굴의 높은 재현성:</text>
    <text x="25" y="214" fill="#94a3b8" font-size="10.5">단일 모달 대비 위양성 80% 이상 제거</text>

    <!-- Bottom Badge -->
    <rect x="15" y="235" width="230" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="130" y="250" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">✓ 고신뢰도 유전자 조절망 완성</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_tn5_tagmentation_mechanism.svg", make_01_tn5_tagmentation),
        ("02_fragment_file_and_tss_enrichment.svg", make_02_fragment_and_tss),
        ("03_lsi_tfidf_svd_geometry.svg", make_03_lsi_tfidf),
        ("04_wnn_modality_weighting.svg", make_04_wnn_weighting),
        ("05_ar_enhancer_duplication_peak_to_gene.svg", make_05_ar_enhancer_duplication),
        ("06_scenicplus_triplet_egrn.svg", make_06_scenicplus_triplet),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 01 Multiome publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

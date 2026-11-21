#!/usr/bin/env python3
"""
scripts/generate_diagrams.py
============================
Generates publication-grade SVG diagrams for snRNA-seq prostate cancer learning materials.
Covers wet-lab isolation, pre-mRNA alignment, QC, inferCNV, lineage plasticity,
multiome/spatial deconvolution, downstream manifold/Leiden, and pseudobulk DEG.
"""

import os
import xml.etree.ElementTree as ET


def generate_svg_01(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="gradAccent" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="100%" stop-color="#818cf8"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="35" text-anchor="middle" fill="#38bdf8" font-size="20" font-weight="700">전립선암 조직 단일 핵 분리(Nuclear Isolation) 및 정제 메커니즘</text>
  <text x="500" y="55" text-anchor="middle" fill="#94a3b8" font-size="12">Cold Dounce Lysis (4°C) ➔ 세포질 선택적 용해 ➔ OptiPrep 밀도 구배 정제 ➔ 10x Chromium 캡슐화</text>

  <!-- Step 1 -->
  <g transform="translate(30, 80)">
    <rect width="200" height="230" rx="10" fill="url(#gradCard)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="200" height="30" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="100" y="20" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 동결 조직 생검 분쇄</text>
    
    <!-- Graphic: Frozen Biopsy Core -->
    <rect x="35" y="45" width="130" height="30" rx="6" fill="#0284c7" stroke="#38bdf8" stroke-width="1"/>
    <text x="100" y="65" text-anchor="middle" fill="#ffffff" font-size="11" font-weight="600">Flash-Frozen Biopsy (-80°C)</text>
    
    <text x="15" y="105" fill="#cbd5e1" font-size="11" font-weight="600">• 조직 파쇄:</text>
    <text x="25" y="123" fill="#94a3b8" font-size="11">Dounce Tissue Grinder</text>
    <text x="15" y="145" fill="#cbd5e1" font-size="11" font-weight="600">• 저온 조건 유지:</text>
    <text x="25" y="163" fill="#94a3b8" font-size="11">4°C 유지 (얼음 배스 상)</text>
    <text x="15" y="185" fill="#cbd5e1" font-size="11" font-weight="600">• 목적:</text>
    <text x="25" y="203" fill="#94a3b8" font-size="11">전사 동결 &amp; 스트레스 차단</text>
  </g>

  <!-- Arrow 1 -->
  <path d="M 240 195 L 265 195" stroke="#38bdf8" stroke-width="2.5"/>
  <polygon points="265,190 275,195 265,200" fill="#38bdf8"/>

  <!-- Step 2 -->
  <g transform="translate(285, 80)">
    <rect width="200" height="230" rx="10" fill="url(#gradCard)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="200" height="30" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="100" y="20" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">2. 계면활성제 화학 용해</text>
    
    <rect x="35" y="45" width="130" height="30" rx="6" fill="#4338ca" stroke="#818cf8" stroke-width="1"/>
    <text x="100" y="65" text-anchor="middle" fill="#ffffff" font-size="11" font-weight="600">Lysis Buffer (NP-40 0.1%)</text>
    
    <text x="15" y="105" fill="#cbd5e1" font-size="11" font-weight="600">• 세포막(Plasma Memb.):</text>
    <text x="25" y="123" fill="#f43f5e" font-size="11">지질 이중층 용해 (파괴)</text>
    <text x="15" y="145" fill="#cbd5e1" font-size="11" font-weight="600">• 핵막(Nuclear Memb.):</text>
    <text x="25" y="163" fill="#10b981" font-size="11">라미나 지지 (구조 보존)</text>
    <text x="15" y="185" fill="#cbd5e1" font-size="11" font-weight="600">• RNase Inhibitor:</text>
    <text x="25" y="203" fill="#94a3b8" font-size="11">0.4 U/µL (내인성 RNase 억제)</text>
  </g>

  <!-- Arrow 2 -->
  <path d="M 495 195 L 520 195" stroke="#38bdf8" stroke-width="2.5"/>
  <polygon points="520,190 530,195 520,200" fill="#38bdf8"/>

  <!-- Step 3 -->
  <g transform="translate(540, 80)">
    <rect width="200" height="230" rx="10" fill="url(#gradCard)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="200" height="30" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="100" y="20" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">3. 밀도 구배 원심분리</text>
    
    <rect x="35" y="45" width="130" height="30" rx="6" fill="#065f46" stroke="#34d399" stroke-width="1"/>
    <text x="100" y="65" text-anchor="middle" fill="#ffffff" font-size="11" font-weight="600">OptiPrep™ (Iodixanol)</text>
    
    <text x="15" y="105" fill="#cbd5e1" font-size="11" font-weight="600">• 29%/35% 구배:</text>
    <text x="25" y="123" fill="#94a3b8" font-size="11">침강 계수 기반 분리</text>
    <text x="15" y="145" fill="#cbd5e1" font-size="11" font-weight="600">• 제거 대상:</text>
    <text x="25" y="163" fill="#f43f5e" font-size="11">세포질 파편, 지질, 무기질</text>
    <text x="15" y="185" fill="#cbd5e1" font-size="11" font-weight="600">• 회수 대상:</text>
    <text x="25" y="203" fill="#10b981" font-size="11">계면 상 순수 온전한 핵 펠렛</text>
  </g>

  <!-- Arrow 3 -->
  <path d="M 750 195 L 775 195" stroke="#38bdf8" stroke-width="2.5"/>
  <polygon points="775,190 785,195 775,200" fill="#38bdf8"/>

  <!-- Step 4 -->
  <g transform="translate(795, 80)">
    <rect width="175" height="230" rx="10" fill="url(#gradCard)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="175" height="30" rx="10" fill="#34d399" fill-opacity="0.2"/>
    <text x="87" y="20" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">4. 단일 핵 캡슐화</text>
    
    <circle cx="87" cy="65" r="22" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="2"/>
    <circle cx="87" cy="65" r="10" fill="#38bdf8"/>
    
    <text x="12" y="115" fill="#cbd5e1" font-size="11" font-weight="600">• DAPI 염색 검증:</text>
    <text x="20" y="133" fill="#38bdf8" font-size="11">형광 현미경 무결성 확인</text>
    <text x="12" y="155" fill="#cbd5e1" font-size="11" font-weight="600">• 10x Chromium:</text>
    <text x="20" y="173" fill="#94a3b8" font-size="11">Next GEM 유세포 액적</text>
    <text x="12" y="195" fill="#cbd5e1" font-size="11" font-weight="600">• 포획 대상:</text>
    <text x="20" y="213" fill="#10b981" font-size="11">Nucleus + Barcoded Bead</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def generate_svg_02(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 340" width="1000" height="340" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">핵 전사체 Pre-mRNA 특성과 --include-introns 정렬 원리</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">단일 세포(scRNA-seq, Exon 위주) 대비 단일 핵(snRNA-seq, Intron 40~60%) 정량 메커니즘</text>

  <!-- Genomic DNA & Transcript Model -->
  <g transform="translate(60, 75)">
    <!-- Genomic DNA -->
    <text x="0" y="20" fill="#94a3b8" font-size="12" font-weight="700">Genomic DNA:</text>
    <line x1="120" y1="16" x2="880" y2="16" stroke="#475569" stroke-width="4"/>
    
    <!-- Exon 1 -->
    <rect x="180" y="4" width="90" height="24" rx="4" fill="#0284c7" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="225" y="20" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">Exon 1</text>
    
    <!-- Intron 1 -->
    <text x="365" y="20" text-anchor="middle" fill="#94a3b8" font-size="11">Intron 1 (Pre-mRNA)</text>
    
    <!-- Exon 2 -->
    <rect x="460" y="4" width="110" height="24" rx="4" fill="#0284c7" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="515" y="20" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">Exon 2</text>
    
    <!-- Intron 2 -->
    <text x="670" y="20" text-anchor="middle" fill="#94a3b8" font-size="11">Intron 2 (Pre-mRNA)</text>
    
    <!-- Exon 3 -->
    <rect x="760" y="4" width="80" height="24" rx="4" fill="#0284c7" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="800" y="20" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">Exon 3</text>
  </g>

  <!-- Comparison Panel -->
  <g transform="translate(60, 140)">
    <!-- Standard scRNA-seq Panel -->
    <rect x="0" y="0" width="420" height="170" rx="8" fill="#1e293b" stroke="#f43f5e" stroke-width="1.5"/>
    <text x="20" y="28" fill="#f43f5e" font-size="14" font-weight="700">A. 기존 scRNA-seq 정렬 (Exon-Only Counting)</text>
    
    <text x="20" y="58" fill="#cbd5e1" font-size="12">• 엑손 리드 (Exon Reads): <tspan fill="#10b981" font-weight="700">카운트 인정 (+1 UMI)</tspan></text>
    <text x="20" y="82" fill="#cbd5e1" font-size="12">• 인트론 리드 (Intronic Reads): <tspan fill="#f43f5e" font-weight="700">폐기 (Discarded, 0 Count)</tspan></text>
    
    <rect x="20" y="105" width="380" height="48" rx="6" fill="#0f172a" stroke="#334155"/>
    <text x="35" y="125" fill="#f43f5e" font-size="11" font-weight="700">⚠️ snRNA-seq 적용 시 치명적 결과:</text>
    <text x="35" y="142" fill="#94a3b8" font-size="11">핵 전사체의 50% 이상 소실 ➔ 유전자 검출 수 극단적 급감!</text>
  </g>

  <g transform="translate(520, 140)">
    <!-- snRNA-seq Panel -->
    <rect x="0" y="0" width="420" height="170" rx="8" fill="#1e293b" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="20" y="28" fill="#38bdf8" font-size="14" font-weight="700">B. snRNA-seq 필수 정렬 (--include-introns=true)</text>
    
    <text x="20" y="58" fill="#cbd5e1" font-size="12">• 엑손 리드 (Exon Reads): <tspan fill="#10b981" font-weight="700">카운트 인정 (+1 UMI)</tspan></text>
    <text x="20" y="82" fill="#cbd5e1" font-size="12">• 인트론 리드 (Intron Reads): <tspan fill="#10b981" font-weight="700">카운트 전면 인정 (+1 UMI)</tspan></text>
    
    <rect x="20" y="105" width="380" height="48" rx="6" fill="#0f172a" stroke="#0284c7"/>
    <text x="35" y="125" fill="#38bdf8" font-size="11" font-weight="700">✓ 검증된 성과 지표:</text>
    <text x="35" y="142" fill="#cbd5e1" font-size="11">검출 유전자 수 1.8~2.5배 증가 ➔ 미세 종양 아형 분별력 확보</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def generate_svg_03(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 350" width="1000" height="350" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">종양 snRNA-seq 특화 QC: 미토콘드리아 기준 &amp; Ambient RNA 메커니즘</text>

  <!-- Panel 1: Mito QC -->
  <g transform="translate(50, 60)">
    <rect width="420" height="260" rx="10" fill="#1e293b" stroke="#334155" stroke-width="1.5"/>
    <text x="20" y="30" fill="#38bdf8" font-size="15" font-weight="700">1. 미토콘드리아 비율(Mito%) 인지 전환</text>
    
    <!-- scRNA box -->
    <rect x="20" y="55" width="380" height="85" rx="6" fill="#0f172a"/>
    <text x="35" y="78" fill="#e2e8f0" font-size="12" font-weight="700">• 단일 세포 (scRNA-seq):</text>
    <text x="45" y="98" fill="#94a3b8" font-size="11">세포막 파열 시 세포질 탈출 ➔ 미토콘드리아 잔류 (사멸 지표)</text>
    <text x="45" y="118" fill="#fbbf24" font-size="11" font-weight="600">허용 임계치: 5% ~ 20% 이하 (암 조직 허용 범위 넓음)</text>
    
    <!-- snRNA box -->
    <rect x="20" y="155" width="380" height="90" rx="6" fill="#0f172a" stroke="#10b981" stroke-width="1"/>
    <text x="35" y="178" fill="#10b981" font-size="12" font-weight="700">• 단일 핵 (snRNA-seq):</text>
    <text x="45" y="198" fill="#94a3b8" font-size="11">미토콘드리아는 세포질 소기관 ➔ 순수 핵에는 없어야 정상!</text>
    <text x="45" y="218" fill="#f43f5e" font-size="11" font-weight="600">엄격한 임계치: 1.0% ~ 2.0% 미만 (초과 시 세포질 오염 탈락)</text>
  </g>

  <!-- Panel 2: Ambient RNA & CellBender -->
  <g transform="translate(520, 60)">
    <rect width="430" height="260" rx="10" fill="#1e293b" stroke="#334155" stroke-width="1.5"/>
    <text x="20" y="30" fill="#38bdf8" font-size="15" font-weight="700">2. Ambient RNA 오염과 CellBender VAE 분해</text>
    
    <text x="20" y="60" fill="#cbd5e1" font-size="12">• 전립선 분비성 단백질 전사체: <tspan fill="#f43f5e" font-weight="700">KLK3(PSA), ACPP, MSMB</tspan></text>
    <text x="20" y="80" fill="#94a3b8" font-size="11">조직 균질화 시 액체 버퍼(Soup) 전체로 초고농도 누출</text>

    <!-- Formula Box -->
    <rect x="20" y="100" width="390" height="50" rx="6" fill="#0f172a" stroke="#0284c7"/>
    <text x="35" y="122" fill="#38bdf8" font-size="11" font-family="monospace">Observed Count = True_Nucleus(λ) + Droplet_Soup(ρ·μ)</text>
    <text x="35" y="140" fill="#94a3b8" font-size="10">CellBender VAE가 액적별 배경 오염(ρ·μ)을 정밀 감산 추정</text>

    <text x="20" y="175" fill="#cbd5e1" font-size="12">• 미보정 방치 시:</text>
    <text x="35" y="193" fill="#f43f5e" font-size="11">T세포, 대식세포, CAF가 KLK3·AR을 발현하는 가짜 잡종세포 출현</text>
    
    <text x="20" y="218" fill="#10b981" font-size="12">• CellBender 보정 후:</text>
    <text x="35" y="236" fill="#34d399" font-size="11">면역/기질세포의 배경 신호 소멸 ➔ 순수 종양 상피세포만 KLK3 보존</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def generate_svg_04(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 340" width="1000" height="340" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">inferCNV 원리: 전립선 상피세포 염색체 이수성(Aneuploidy) 판별</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">정상 2배체 기질/면역 세포(Reference) vs 악성 선암 클론(Malignant) 염색체 복제수 변이</text>

  <!-- Chromosome Coordinates Axis -->
  <g transform="translate(60, 75)">
    <rect x="0" y="0" width="880" height="25" rx="4" fill="#1e293b"/>
    <text x="40" y="17" fill="#94a3b8" font-size="11">Chr 1</text>
    <text x="130" y="17" fill="#94a3b8" font-size="11">Chr 3</text>
    <text x="250" y="17" fill="#94a3b8" font-size="11">Chr 8</text>
    <text x="360" y="17" fill="#94a3b8" font-size="11">Chr 10</text>
    <text x="470" y="17" fill="#94a3b8" font-size="11">Chr 13</text>
    <text x="590" y="17" fill="#94a3b8" font-size="11">Chr 17</text>
    <text x="720" y="17" fill="#94a3b8" font-size="11">Chr 21</text>
    <text x="820" y="17" fill="#94a3b8" font-size="11">Chr X</text>
  </g>

  <!-- Diploid Reference Track -->
  <g transform="translate(60, 115)">
    <text x="0" y="25" fill="#38bdf8" font-size="12" font-weight="700">2배체 정상 대조군</text>
    <text x="0" y="42" fill="#94a3b8" font-size="10">(T세포, 혈관내피세포)</text>
    <rect x="150" y="10" width="730" height="35" rx="4" fill="#0f172a" stroke="#334155"/>
    <line x1="150" y1="27" x2="880" y2="27" stroke="#10b981" stroke-width="2.5" stroke-dasharray="6,2"/>
    <text x="515" y="32" text-anchor="middle" fill="#10b981" font-size="11" font-weight="600">중립 복제수 상태 (Relative Ratio = 1.0, 정상 2배체)</text>
  </g>

  <!-- Malignant Prostate Clone Track -->
  <g transform="translate(60, 185)">
    <text x="0" y="25" fill="#f43f5e" font-size="12" font-weight="700">악성 선암 클론</text>
    <text x="0" y="42" fill="#94a3b8" font-size="10">(Malignant Luminal)</text>
    <rect x="150" y="10" width="730" height="75" rx="4" fill="#0f172a" stroke="#f43f5e" stroke-width="1"/>
    
    <!-- Neutral baseline -->
    <line x1="150" y1="47" x2="880" y2="47" stroke="#475569" stroke-width="1"/>
    
    <!-- 8p Deletion -->
    <rect x="330" y="47" width="40" height="25" fill="#38bdf8" fill-opacity="0.8"/>
    <text x="350" y="85" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">8p Del</text>
    
    <!-- 8q Amplification (MYC) -->
    <rect x="375" y="18" width="45" height="29" fill="#f43f5e" fill-opacity="0.8"/>
    <text x="397" y="14" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="700">8q Gain(MYC)</text>
    
    <!-- 10q Deletion (PTEN) -->
    <rect x="440" y="47" width="50" height="25" fill="#38bdf8" fill-opacity="0.8"/>
    <text x="465" y="85" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">10q Del(PTEN)</text>
    
    <!-- 13q Deletion (RB1) -->
    <rect x="520" y="47" width="45" height="25" fill="#38bdf8" fill-opacity="0.8"/>
    <text x="542" y="85" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">13q Del(RB1)</text>

    <!-- 17p Deletion (TP53) -->
    <rect x="620" y="47" width="35" height="25" fill="#38bdf8" fill-opacity="0.8"/>
    <text x="637" y="85" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">17p Del(TP53)</text>
    
    <!-- Xq Amplification (AR in mCRPC) -->
    <rect x="830" y="18" width="40" height="29" fill="#f43f5e" fill-opacity="0.8"/>
    <text x="850" y="14" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="700">AR Amp</text>
  </g>

  <!-- Legend -->
  <g transform="translate(250, 300)">
    <rect x="0" y="0" width="18" height="12" fill="#f43f5e"/>
    <text x="25" y="10" fill="#cbd5e1" font-size="11">염색체 영역 증폭 (Gain / Amplification)</text>
    <rect x="270" y="0" width="18" height="12" fill="#38bdf8"/>
    <text x="295" y="10" fill="#cbd5e1" font-size="11">염색체 영역 결손 (Loss / Deletion)</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def generate_svg_05(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 340" width="1000" height="340" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">mCRPC 치료 저항성 및 신경내분비(t-NEPC) 계통 가소성 궤적</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">안드로겐 차단 요법(ADT) 압력 하에서의 관강 선암(AR+) ➔ 중간 가소성 상태 ➔ t-NEPC(AR-) 전환</text>

  <!-- State 1 -->
  <g transform="translate(60, 80)">
    <rect width="250" height="220" rx="10" fill="#1e293b" stroke="#38bdf8" stroke-width="2"/>
    <rect x="0" y="0" width="250" height="35" rx="10" fill="#0284c7" fill-opacity="0.3"/>
    <text x="125" y="23" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 관강 선암 (Adeno)</text>
    
    <text x="20" y="60" fill="#10b981" font-size="12" font-weight="700">• AR Signaling: <tspan fill="#38bdf8">High +++</tspan></text>
    <text x="20" y="82" fill="#cbd5e1" font-size="11">KLK3(PSA), ACPP, TMPRSS2</text>
    <text x="20" y="110" fill="#cbd5e1" font-size="12" font-weight="700">• 파이오니어 전사 인자:</text>
    <text x="20" y="130" fill="#94a3b8" font-size="11">FOXA1, HOXB13 고활성</text>
    <text x="20" y="158" fill="#cbd5e1" font-size="12" font-weight="700">• 임상 상태:</text>
    <text x="20" y="178" fill="#94a3b8" font-size="11">호르몬 치료 반응성 / 초기</text>
    <text x="20" y="198" fill="#94a3b8" font-size="11">신경 마커 음성 (SYP-, CHGA-)</text>
  </g>

  <!-- Transition Arrow 1 -->
  <g transform="translate(320, 175)">
    <path d="M 0 0 L 40 0" stroke="#f59e0b" stroke-width="3"/>
    <polygon points="40,-5 50,0 40,5" fill="#f59e0b"/>
    <text x="25" y="-10" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="700">Enzalutamide</text>
    <text x="25" y="20" text-anchor="middle" fill="#f59e0b" font-size="10">선택압 투여</text>
  </g>

  <!-- State 2 -->
  <g transform="translate(380, 80)">
    <rect width="250" height="220" rx="10" fill="#1e293b" stroke="#f59e0b" stroke-width="2"/>
    <rect x="0" y="0" width="250" height="35" rx="10" fill="#f59e0b" fill-opacity="0.3"/>
    <text x="125" y="23" text-anchor="middle" fill="#f59e0b" font-size="13" font-weight="700">2. 중간 가소성 상태 (Plastic)</text>
    
    <text x="20" y="60" fill="#f59e0b" font-size="12" font-weight="700">• AR Signaling: <tspan fill="#f43f5e">저하/소실 중</tspan></text>
    <text x="20" y="82" fill="#cbd5e1" font-size="11">AR+/-, KLK3 급감</text>
    <text x="20" y="110" fill="#cbd5e1" font-size="12" font-weight="700">• 유전체 결손 &amp; 후성 재편:</text>
    <text x="20" y="130" fill="#f43f5e" font-size="11">TP53 Loss &amp; RB1 Loss</text>
    <text x="20" y="150" fill="#38bdf8" font-size="11">EZH2 과발현, SOX2 획득</text>
    <text x="20" y="178" fill="#cbd5e1" font-size="12" font-weight="700">• REST 억제자 억제 해제:</text>
    <text x="20" y="198" fill="#94a3b8" font-size="11">신경계 유전자 탈억제 시작</text>
  </g>

  <!-- Transition Arrow 2 -->
  <g transform="translate(640, 175)">
    <path d="M 0 0 L 40 0" stroke="#f43f5e" stroke-width="3"/>
    <polygon points="40,-5 50,0 40,5" fill="#f43f5e"/>
    <text x="25" y="-10" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="700">계통 전환</text>
    <text x="25" y="20" text-anchor="middle" fill="#f43f5e" font-size="10">Lineage Switch</text>
  </g>

  <!-- State 3 -->
  <g transform="translate(700, 80)">
    <rect width="250" height="220" rx="10" fill="#1e293b" stroke="#f43f5e" stroke-width="2"/>
    <rect x="0" y="0" width="250" height="35" rx="10" fill="#f43f5e" fill-opacity="0.3"/>
    <text x="125" y="23" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">3. 신경내분비 전립선암 (t-NEPC)</text>
    
    <text x="20" y="60" fill="#f43f5e" font-size="12" font-weight="700">• AR Signaling: <tspan fill="#f43f5e">완전 소실 (AR-)</tspan></text>
    <text x="20" y="82" fill="#94a3b8" font-size="11">호르몬 요법 완전 불응</text>
    <text x="20" y="110" fill="#10b981" font-size="12" font-weight="700">• 신경계 마스터 TF:</text>
    <text x="20" y="130" fill="#34d399" font-size="11">ASCL1, NEUROD1, INSM1</text>
    <text x="20" y="158" fill="#10b981" font-size="12" font-weight="700">• 신경 분비성 마커:</text>
    <text x="20" y="178" fill="#34d399" font-size="11">SYP, CHGA, ENO2, NCAM1</text>
    <text x="20" y="198" fill="#f43f5e" font-size="11" font-weight="600">예후 극히 불량 / 백금 화학요법</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def generate_svg_06(filepath):
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 340" width="1000" height="340" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">10x Multiome (snRNA + snATAC) 및 공간전사체 디콘볼루션</text>

  <!-- Left Panel: 10x Multiome -->
  <g transform="translate(50, 60)">
    <rect width="420" height="260" rx="10" fill="#1e293b" stroke="#334155" stroke-width="1.5"/>
    <text x="20" y="30" fill="#38bdf8" font-size="14" font-weight="700">A. 10x Multiome: 단일 핵 동시 이중 오믹스</text>
    
    <circle cx="65" cy="85" r="28" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="2"/>
    <text x="65" y="90" text-anchor="middle" fill="#fff" font-size="10" font-weight="700">Single Nucleus</text>
    
    <path d="M 105 75 L 145 60" stroke="#38bdf8" stroke-width="2"/>
    <path d="M 105 95 L 145 110" stroke="#818cf8" stroke-width="2"/>
    
    <rect x="150" y="45" width="240" height="35" rx="6" fill="#0f172a" stroke="#38bdf8"/>
    <text x="160" y="67" fill="#38bdf8" font-size="11" font-weight="600">snRNA: pre-mRNA 발현량 (GEX)</text>
    
    <rect x="150" y="95" width="240" height="35" rx="6" fill="#0f172a" stroke="#818cf8"/>
    <text x="160" y="117" fill="#818cf8" font-size="11" font-weight="600">snATAC: 열린 염색질 (Tn5 Tagment)</text>
    
    <!-- Integration Box -->
    <rect x="20" y="150" width="380" height="90" rx="6" fill="#0f172a" stroke="#10b981"/>
    <text x="35" y="172" fill="#10b981" font-size="12" font-weight="700">WNN (Weighted Nearest Neighbor) 융합:</text>
    <text x="35" y="192" fill="#cbd5e1" font-size="11">• 세포별 정보량 가중치 학습: w_RNA + w_ATAC = 1.0</text>
    <text x="35" y="210" fill="#94a3b8" font-size="11">• Peak-to-Gene 링크: AR 인핸서 증폭 ➔ 전사 활성 입증</text>
    <text x="35" y="228" fill="#94a3b8" font-size="11">• ATAC QC: TSS Enrichment &gt; 4-8, 147bp 밴딩</text>
  </g>

  <!-- Right Panel: Spatial Deconvolution -->
  <g transform="translate(520, 60)">
    <rect width="430" height="260" rx="10" fill="#1e293b" stroke="#334155" stroke-width="1.5"/>
    <text x="20" y="30" fill="#38bdf8" font-size="14" font-weight="700">B. 공간전사체 (10x Visium) 디콘볼루션</text>
    
    <rect x="20" y="55" width="390" height="70" rx="6" fill="#0f172a"/>
    <text x="35" y="78" fill="#cbd5e1" font-size="12">• 공간 스팟 한계: 55 µm (스팟당 5~10개 세포 혼합)</text>
    <text x="35" y="98" fill="#94a3b8" font-size="11">• 해결책: snRNA-seq 순수 단일 핵 레퍼런스 투영</text>
    <text x="35" y="115" fill="#38bdf8" font-size="11" font-weight="600">알고리즘: cell2location / RCTD 베이지안 확률 역추정</text>
    
    <rect x="20" y="140" width="390" height="100" rx="6" fill="#0f172a" stroke="#0284c7"/>
    <text x="35" y="162" fill="#38bdf8" font-size="12" font-weight="700">전립선암 미세환경 공간 규명 성과:</text>
    <text x="35" y="182" fill="#f43f5e" font-size="11">• 면역 배제 (Immune Exclusion): CD8 T세포가 CAF 장벽에 고립</text>
    <text x="35" y="200" fill="#34d399" font-size="11">• 침습성 경계면: Malignant Luminal 세포와 myCAF 접점</text>
    <text x="35" y="218" fill="#fbbf24" font-size="11">• 골전이 니치: 조골세포와 암세포 간 Wnt/RANKL 상호작용</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def make_07_downstream_manifold(filepath):
    """
    Generates docs/images/07_downstream_manifold_and_leiden.svg (viewBox="0 0 1000 360")
    Illustrates:
    - Panel 1: Variance stabilizing normalization & 2,500 HVGs removing depth bias
    - Panel 2: PCA (35 PCs) -> kNN graph -> UMAP topological manifold (warning on Euclidean distances)
    - Panel 3: Louvain disconnected bridge bug vs Leiden guaranteed connected sub-communities (res=0.8)
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard7" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">전립선 snRNA-seq 다운스트림 분석: 정규화 ➔ 매니폴드 투영 ➔ Leiden 군집화</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">Pearson Residuals / Log1p ➔ 35 PCs &amp; kNN Graph ➔ UMAP 위상 다양체 ➔ Leiden 노드 정제(연결 보장)</text>

  <!-- Panel 1: Variance Stabilization & HVG -->
  <g transform="translate(30, 70)">
    <rect width="295" height="265" rx="10" fill="url(#gradCard7)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="147" y="22" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 분산 안정화 정규화 &amp; 2,500 HVG</text>

    <!-- Formula Box -->
    <rect x="15" y="45" width="265" height="42" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <text x="25" y="62" fill="#38bdf8" font-size="10.5" font-family="monospace">Log1p: log(Count/Total * 10⁴ + 1)</text>
    <text x="25" y="78" fill="#34d399" font-size="10" font-family="monospace">Pearson Residuals: (y - μ) / √(μ + μ²/θ)</text>

    <!-- Bullets -->
    <text x="15" y="107" fill="#cbd5e1" font-size="11" font-weight="600">• 시퀀싱 심도(Sequencing Depth) 편향 제거:</text>
    <text x="25" y="125" fill="#94a3b8" font-size="10.5">핵별 총 UMI 변동에 의한 기술적 왜곡 보정</text>

    <text x="15" y="147" fill="#cbd5e1" font-size="11" font-weight="600">• 고변이 유전자(2,500 HVG) 선별:</text>
    <text x="25" y="165" fill="#94a3b8" font-size="10.5">평균-분산 적합(seurat_v3) 기반 전립선 핵심</text>
    <text x="25" y="181" fill="#38bdf8" font-size="10.5">마커(AR, KLK3, AMACR, SYP) 보존</text>

    <text x="15" y="203" fill="#cbd5e1" font-size="11" font-weight="600">• 희소 행렬(Sparse Matrix) 최적화:</text>
    <text x="25" y="221" fill="#94a3b8" font-size="10.5">90%+ Zero 인플레이션 유지, 메모리 경량화</text>

    <!-- Badge -->
    <rect x="15" y="235" width="265" height="22" rx="4" fill="#0369a1" fill-opacity="0.3"/>
    <text x="147" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">✓ 과분산(Overdispersion) 보정 완료</text>
  </g>

  <!-- Flow Arrow 1 -->
  <path d="M 330 200 L 350 200" stroke="#38bdf8" stroke-width="2.5"/>
  <polygon points="350,195 360,200 350,205" fill="#38bdf8"/>

  <!-- Panel 2: PCA to UMAP Manifold Geometry -->
  <g transform="translate(365, 70)">
    <rect width="295" height="265" rx="10" fill="url(#gradCard7)" stroke="#818cf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#818cf8" fill-opacity="0.2"/>
    <text x="147" y="22" text-anchor="middle" fill="#818cf8" font-size="13" font-weight="700">2. PCA ➔ kNN ➔ UMAP 다양체 기하학</text>

    <!-- Visual schematic: Scree / kNN / UMAP -->
    <rect x="15" y="45" width="265" height="42" rx="6" fill="#0f172a" stroke="#4338ca" stroke-width="1"/>
    <text x="25" y="62" fill="#818cf8" font-size="11" font-weight="600">35 직교 PC ➔ kNN(k=15) ➔ 2D UMAP</text>
    <text x="25" y="78" fill="#94a3b8" font-size="10">Elbow Point: 전역 분산 85%+ 포착 다양체 구축</text>

    <!-- Bullets -->
    <text x="15" y="107" fill="#cbd5e1" font-size="11" font-weight="600">• 직교 주성분 분석(PCA, 35 PCs):</text>
    <text x="25" y="125" fill="#94a3b8" font-size="10.5">고차원 유전자 공간의 기술적 노이즈 압축</text>

    <text x="15" y="147" fill="#cbd5e1" font-size="11" font-weight="600">• 국소 연결 그래프(kNN Graph, k=15):</text>
    <text x="25" y="165" fill="#94a3b8" font-size="10.5">코사인/유클리드 거리 기반 단일 핵 이웃 결속</text>

    <!-- Crucial Warning Box -->
    <rect x="15" y="180" width="265" height="75" rx="6" fill="#450a0a" stroke="#ef4444" stroke-width="1.2"/>
    <text x="25" y="198" fill="#f87171" font-size="10.5" font-weight="700">⚠️ UMAP 위상 해석 치명적 주의점:</text>
    <text x="25" y="216" fill="#fca5a5" font-size="10">• 2D UMAP 간격 ≠ 고차원 유클리드 거리</text>
    <text x="25" y="232" fill="#fca5a5" font-size="10">• 클러스터 간 절대 거리는 유사도가 아님</text>
    <text x="25" y="248" fill="#cbd5e1" font-size="9.5">국소 연속성(Local Neighborhood)만 신뢰할 것!</text>
  </g>

  <!-- Flow Arrow 2 -->
  <path d="M 665 200 L 685 200" stroke="#818cf8" stroke-width="2.5"/>
  <polygon points="685,195 695,200 685,205" fill="#818cf8"/>

  <!-- Panel 3: Louvain Flaw vs Leiden Node Refinement -->
  <g transform="translate(700, 70)">
    <rect width="270" height="265" rx="10" fill="url(#gradCard7)" stroke="#10b981" stroke-width="1.5"/>
    <rect x="0" y="0" width="270" height="32" rx="10" fill="#10b981" fill-opacity="0.2"/>
    <text x="135" y="22" text-anchor="middle" fill="#10b981" font-size="13" font-weight="700">3. Louvain 결함 vs Leiden 정제</text>

    <!-- Contrast Diagram -->
    <rect x="12" y="42" width="118" height="60" rx="6" fill="#0f172a" stroke="#ef4444" stroke-width="1"/>
    <text x="71" y="58" text-anchor="middle" fill="#f87171" font-size="10" font-weight="700">Louvain 알고리즘</text>
    <circle cx="35" cy="80" r="10" fill="#ef4444" fill-opacity="0.3" stroke="#ef4444"/>
    <circle cx="105" cy="80" r="10" fill="#ef4444" fill-opacity="0.3" stroke="#ef4444"/>
    <line x1="45" y1="80" x2="95" y2="80" stroke="#ef4444" stroke-dasharray="3,3" stroke-width="1.5"/>
    <text x="71" y="96" text-anchor="middle" fill="#ef4444" font-size="8.5">단절된 브릿지 결함</text>

    <rect x="138" y="42" width="120" height="60" rx="6" fill="#0f172a" stroke="#10b981" stroke-width="1"/>
    <text x="198" y="58" text-anchor="middle" fill="#34d399" font-size="10" font-weight="700">Leiden (Node Refine)</text>
    <circle cx="165" cy="80" r="10" fill="#10b981" fill-opacity="0.3" stroke="#10b981"/>
    <circle cx="230" cy="80" r="10" fill="#10b981" fill-opacity="0.3" stroke="#10b981"/>
    <line x1="175" y1="80" x2="220" y2="80" stroke="#10b981" stroke-width="2"/>
    <text x="198" y="96" text-anchor="middle" fill="#34d399" font-size="8.5">연결성 수학적 보장</text>

    <!-- Bullets -->
    <text x="12" y="122" fill="#cbd5e1" font-size="11" font-weight="600">• Louvain의 치명적 결함:</text>
    <text x="20" y="139" fill="#fca5a5" font-size="10.5">모듈성 최적화 중 내부가 끊어진(Disconnected)</text>
    <text x="20" y="154" fill="#fca5a5" font-size="10.5">비연결 하위 커뮤니티가 동일 클러스터로 병합</text>

    <text x="12" y="176" fill="#cbd5e1" font-size="11" font-weight="600">• Leiden 3단계 노드 정제 혁신:</text>
    <text x="20" y="193" fill="#34d399" font-size="10.5">1) Local Move ➔ 2) Node Refine ➔</text>
    <text x="20" y="208" fill="#34d399" font-size="10.5">3) Graph Aggregate 로 연결성 100% 보장</text>

    <!-- Resolution recommendation -->
    <rect x="12" y="222" width="246" height="35" rx="5" fill="#064e3b" fill-opacity="0.4" stroke="#10b981" stroke-width="1"/>
    <text x="135" y="237" text-anchor="middle" fill="#34d399" font-size="10" font-weight="700">전립선암 최적 해상도: resolution=0.8</text>
    <text x="135" y="250" text-anchor="middle" fill="#a7f3d0" font-size="9">악성 Luminal 아형 vs TME 기질세포 정밀 분리</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


def make_08_pseudobulk_vs_pseudoreplication(filepath):
    """
    Generates docs/images/08_pseudobulk_vs_pseudoreplication.svg (viewBox="0 0 1000 360")
    Illustrates:
    - Left Panel (Red Warning): Pseudoreplication Trap with N=10,000 single nuclei inflating patient noise to false p=10^-150
    - Center Panel: Pseudobulk Sum Aggregation (patient_id x cell_type_annotation) on raw integer counts
    - Right Panel (Green Success): Patient-level Negative Binomial GLM (DESeq2/PyDESeq2, N=5 vs 5) modeling true biological dispersion alpha_g
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="gradCard8" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="32" text-anchor="middle" fill="#38bdf8" font-size="19" font-weight="700">단일세포 차등발현 검정: 위양성 가짜 반복(Pseudoreplication) vs 환자 수준 Pseudobulk</text>
  <text x="500" y="52" text-anchor="middle" fill="#94a3b8" font-size="12">단일 핵 N-인플레이션 오류(N=10,000) 극복 ➔ 환자 간 생물학적 분산(α_g)을 반영한 음이항 GLM (DESeq2)</text>

  <!-- Left Panel: Pseudoreplication Trap (Red Warning) -->
  <g transform="translate(30, 70)">
    <rect width="310" height="265" rx="10" fill="url(#gradCard8)" stroke="#ef4444" stroke-width="1.8"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#ef4444" fill-opacity="0.25"/>
    <text x="155" y="22" text-anchor="middle" fill="#f87171" font-size="12.5" font-weight="700">❌ 잘못된 단일세포 수준 검정 (가짜 반복)</text>

    <!-- Schematic: 10,000 cells treated as independent -->
    <rect x="15" y="44" width="280" height="58" rx="6" fill="#450a0a" stroke="#ef4444" stroke-width="1"/>
    <text x="25" y="63" fill="#fca5a5" font-size="11" font-weight="700">5 Responders vs 5 Non-responders</text>
    <text x="25" y="80" fill="#f87171" font-size="10" font-family="monospace">Wilcoxon test: N = 10,000 single nuclei</text>
    <text x="25" y="94" fill="#fca5a5" font-size="10" font-weight="700">위양성 폭발: p-value = 10⁻¹⁵⁰ (가짜 유의성)</text>

    <!-- Technical Critiques -->
    <text x="15" y="122" fill="#cbd5e1" font-size="11" font-weight="600">• 통계적 독립성 가정(i.i.d.) 위반:</text>
    <text x="25" y="140" fill="#94a3b8" font-size="10.5">동일 환자에게서 나온 1,000개 핵은 강한</text>
    <text x="25" y="156" fill="#94a3b8" font-size="10.5">환자 내 상관관계(Intra-donor Correlation) 보유</text>

    <text x="15" y="178" fill="#cbd5e1" font-size="11" font-weight="600">• 환자 특이적 잡음(Patient Noise) 과대 포장:</text>
    <text x="25" y="196" fill="#94a3b8" font-size="10.5">특정 환자 1명의 개체 특이적 유전자가 전체</text>
    <text x="25" y="212" fill="#f87171" font-size="10.5">치료 반응 핵심 바이오마커로 둔갑</text>

    <text x="15" y="234" fill="#cbd5e1" font-size="11" font-weight="600">• 임상 재현 실패율:</text>
    <text x="25" y="252" fill="#ef4444" font-size="10.5" font-weight="700">독립 코호트 검증 시 80~90% 이상 재현 불능</text>
  </g>

  <!-- Center Arrow / Aggregation Mechanism -->
  <g transform="translate(350, 70)">
    <rect width="280" height="265" rx="10" fill="url(#gradCard8)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="280" height="32" rx="10" fill="#38bdf8" fill-opacity="0.2"/>
    <text x="140" y="22" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">⚡ Pseudobulk 집계 메커니즘</text>

    <!-- Formula & Aggregation Card -->
    <rect x="15" y="44" width="250" height="58" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <text x="25" y="63" fill="#38bdf8" font-size="10.5" font-weight="700">patient_id × cell_type_annotation</text>
    <text x="25" y="80" fill="#34d399" font-size="10" font-family="monospace">Y_gj = Σ_{c ∈ Donor_j} Raw_UMI_gc</text>
    <text x="25" y="94" fill="#94a3b8" font-size="9.5">Malignant Luminal 세포군 대상 정수 합산</text>

    <!-- Technical Bullets -->
    <text x="15" y="122" fill="#cbd5e1" font-size="11" font-weight="600">• 원시 정수(Raw Integer UMI) 단순 합산:</text>
    <text x="25" y="140" fill="#94a3b8" font-size="10.5">Log 변환이나 정규화 전 카운트 합산 필수</text>

    <text x="15" y="162" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 핵 희소성(Dropout) 완벽 극복:</text>
    <text x="25" y="180" fill="#94a3b8" font-size="10.5">Zero 드롭아웃 극복 ➔ 신호 대 잡음비(SNR) 급증</text>

    <text x="15" y="202" fill="#cbd5e1" font-size="11" font-weight="600">• 표본 품질 컷오프(QC Filtering):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10.5">환자당 세포 수 ≥ 10~20 nuclei &amp; UMI ≥ 500</text>

    <!-- Bottom summary badge -->
    <rect x="15" y="235" width="250" height="22" rx="4" fill="#0284c7" fill-opacity="0.3"/>
    <text x="140" y="250" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">10,000 Nuclei ➔ 10 Patient Samples (5 vs 5)</text>
  </g>

  <!-- Right Panel: Patient-level DESeq2 GLM (Green Success) -->
  <g transform="translate(640, 70)">
    <rect width="330" height="265" rx="10" fill="url(#gradCard8)" stroke="#10b981" stroke-width="1.8"/>
    <rect x="0" y="0" width="330" height="32" rx="10" fill="#10b981" fill-opacity="0.25"/>
    <text x="165" y="22" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">✅ 환자 수준 음이항 GLM (DESeq2 / PyDESeq2)</text>

    <!-- Schematic: DESeq2 with true biological dispersion -->
    <rect x="15" y="44" width="300" height="58" rx="6" fill="#064e3b" fill-opacity="0.5" stroke="#10b981" stroke-width="1"/>
    <text x="25" y="63" fill="#a7f3d0" font-size="11" font-weight="700">음이항 GLM 설계: ~ batch + condition</text>
    <text x="25" y="80" fill="#34d399" font-size="10" font-family="monospace">Var(Y_gj) = μ_gj + α_g · μ_gj² (N = 5 vs 5)</text>
    <text x="25" y="94" fill="#a7f3d0" font-size="10" font-weight="600">엄밀한 환자 간 생물학적 분산(α_g) 적합</text>

    <!-- Technical Strengths -->
    <text x="15" y="122" fill="#cbd5e1" font-size="11" font-weight="600">• 진정한 통계적 자유도 확보 (N = 5 vs 5):</text>
    <text x="25" y="140" fill="#94a3b8" font-size="10.5">환자 간 개체 차이를 분산 매개변수로 정확히 모델링</text>

    <text x="15" y="162" fill="#cbd5e1" font-size="11" font-weight="600">• 임상 공변량(Covariates) 통제:</text>
    <text x="25" y="180" fill="#94a3b8" font-size="10.5">배치 효과, 치료 이력, 환자 연령 교란 변수 분리</text>

    <text x="15" y="202" fill="#cbd5e1" font-size="11" font-weight="600">• 경험적 베이즈 분산 수축(Empirical Bayes Shrinkage):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10.5">소규모 코호트에서도 안정적인 Fold Change 추정</text>

    <!-- Badge -->
    <rect x="15" y="235" width="300" height="22" rx="4" fill="#065f46" fill-opacity="0.4"/>
    <text x="165" y="250" text-anchor="middle" fill="#34d399" font-size="10" font-weight="700">✓ 임상 재현 가능한 진정한 전립선암 바이오마커 발굴</text>
  </g>
</svg>'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg)


# Aliases to support both naming styles
generate_svg_07 = make_07_downstream_manifold
generate_svg_08 = make_08_pseudobulk_vs_pseudoreplication


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    files_and_generators = [
        ("01_wetlab_nuclear_isolation_workflow.svg", generate_svg_01),
        ("02_premRNA_include_introns_mechanics.svg", generate_svg_02),
        ("03_snRNA_qc_mitochondrial_and_cellbender.svg", generate_svg_03),
        ("04_infercnv_chromosomal_aneuploidy.svg", generate_svg_04),
        ("05_mcrpc_lineage_plasticity_trajectory.svg", generate_svg_05),
        ("06_multiome_wnn_and_spatial_deconvolution.svg", generate_svg_06),
        ("07_downstream_manifold_and_leiden.svg", make_07_downstream_manifold),
        ("08_pseudobulk_vs_pseudoreplication.svg", make_08_pseudobulk_vs_pseudoreplication),
    ]

    for fname, func in files_and_generators:
        target = os.path.join(outdir, fname)
        func(target)
        # XML syntax validation
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 8 publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

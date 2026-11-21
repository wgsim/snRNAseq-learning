#!/usr/bin/env python3
"""
tracks/03_tme_cell_communication/scripts/generate_tme_diagrams.py
=================================================================
Generates publication-grade, high-resolution SVG diagrams for Track 03:
Tumor Microenvironment (TME) & Intercellular Signaling Networks in Prostate Cancer.

Covers:
1. TME Cellular Diversity & 4 Intercellular Signaling Modes
2. LIANA Consensus Framework & Cell-Label Permutation Testing
3. CellChat Communication Probability & Network Centrality
4. NicheNet Intracellular Ligand-to-Target Causal Network
5. MISTy Multiview Spatial Decomposition & Contextual Modeling
6. Prostate Cancer Bone Metastasis Vicious Cycle & Immunosuppression
"""

import os
import xml.etree.ElementTree as ET


def make_01_tme_cellular_ecosystem_and_signaling_types(filepath: str):
    """
    1. TME Cellular Diversity & 4 Intercellular Signaling Modes
    Tumor, myCAF, Endothelial, TAM M2, CD8 T + Endocrine, Paracrine, Juxtacrine, Autocrine + Heterodimer stoichiometry.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="cyanPill" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="100%" stop-color="#0284c7"/>
    </linearGradient>
    <linearGradient id="rosePill" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#f43f5e"/>
      <stop offset="100%" stop-color="#be123c"/>
    </linearGradient>
    <linearGradient id="emeraldPill" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#34d399"/>
      <stop offset="100%" stop-color="#059669"/>
    </linearGradient>
    <marker id="arrowCyan" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrowEmerald" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">종양 미세환경(TME) 세포 생태계 &amp; 4대 세포 간 신호전달 모드</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">mCRPC 이종 세포 생태계, 거리별 신호전달 물리적 분류, 및 다중 소단위체(Heterodimer) 화학량론</text>

  <!-- Panel 1: Cellular Diversity -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. TME 이종 세포 생태계 (Cellular Diversity)</text>

    <!-- Schematic: Heterogeneous Cells -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Cells icons -->
    <circle cx="45" cy="65" r="14" fill="#f43f5e" fill-opacity="0.8" stroke="#f43f5e" stroke-width="1.5"/>
    <text x="45" y="69" text-anchor="middle" fill="#ffffff" font-size="8.5" font-weight="700">Tumor</text>
    <circle cx="85" cy="88" r="12" fill="#fbbf24" fill-opacity="0.8" stroke="#fbbf24" stroke-width="1.5"/>
    <text x="85" y="92" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">myCAF</text>
    <circle cx="130" cy="62" r="11" fill="#38bdf8" fill-opacity="0.8" stroke="#38bdf8" stroke-width="1.5"/>
    <text x="130" y="66" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">Endo</text>
    <circle cx="175" cy="85" r="13" fill="#a855f7" fill-opacity="0.8" stroke="#a855f7" stroke-width="1.5"/>
    <text x="175" y="89" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">TAM M2</text>
    <circle cx="225" cy="66" r="12" fill="#34d399" fill-opacity="0.8" stroke="#34d399" stroke-width="1.5"/>
    <text x="225" y="70" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">CD8 T</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 종양 세포 (Luminal AR+ vs t-NEPC):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">AR 의존성 선암에서 신경내분비계(NEPC) 역분화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 기질 및 면역 틈새 (Stromal &amp; Immune):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">myCAF(CXCL12 장벽), M2 대식세포(면역 억제)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 혈관 및 림프 침윤 구획:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">내피세포(VEGFA 반응) 및 고갈된 CD8+ T세포</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">특징: 단일세포 전사체 기반 세포 상호작용 지도화</text>
  </g>

  <!-- Panel 2: 4 Signaling Modes -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 4대 신호전달 모드 (Signaling Modes)</text>

    <!-- Grid of 4 modes -->
    <g transform="translate(15, 42)">
      <!-- Autocrine -->
      <rect x="0" y="0" width="130" height="40" rx="5" fill="#0f172a" stroke="#34d399" stroke-width="1"/>
      <text x="10" y="16" fill="#34d399" font-size="10" font-weight="700">자가분비 (Autocrine)</text>
      <text x="10" y="30" fill="#94a3b8" font-size="8.5">자체 자극 루프 (TGFB1)</text>

      <!-- Juxtacrine -->
      <rect x="140" y="0" width="130" height="40" rx="5" fill="#0f172a" stroke="#fbbf24" stroke-width="1"/>
      <text x="150" y="16" fill="#fbbf24" font-size="10" font-weight="700">접촉분비 (Juxtacrine)</text>
      <text x="150" y="30" fill="#94a3b8" font-size="8.5">막접촉 (&lt;15µm, Notch/PD1)</text>

      <!-- Paracrine -->
      <rect x="0" y="46" width="130" height="40" rx="5" fill="#0f172a" stroke="#38bdf8" stroke-width="1"/>
      <text x="10" y="62" fill="#38bdf8" font-size="10" font-weight="700">방분비 (Paracrine)</text>
      <text x="10" y="76" fill="#94a3b8" font-size="8.5">국소 확산 (10~250µm)</text>

      <!-- Endocrine -->
      <rect x="140" y="46" width="130" height="40" rx="5" fill="#0f172a" stroke="#f43f5e" stroke-width="1"/>
      <text x="150" y="62" fill="#f43f5e" font-size="10" font-weight="700">내분비 (Endocrine)</text>
      <text x="150" y="76" fill="#94a3b8" font-size="8.5">혈류 순환계 원거리 (&gt;1cm)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="142" fill="#cbd5e1" font-size="11" font-weight="600">• 접촉 분비 (Juxtacrine):</text>
    <text x="25" y="158" fill="#94a3b8" font-size="10">리간드-수용체 세포막 직접 결합 (PD-L1/PD-1)</text>

    <text x="15" y="178" fill="#cbd5e1" font-size="11" font-weight="600">• 방분비 확산 필드 (Paracrine Field):</text>
    <text x="25" y="194" fill="#94a3b8" font-size="10">조직 간극 농도 구배 및 기질 결합(CXCL12-CXCR4)</text>

    <text x="15" y="214" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 거리 감쇠 (Distance Decay):</text>
    <text x="25" y="230" fill="#94a3b8" font-size="10">유효 신호 도달 거리에 따른 통신 확률 보정</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">거리 법칙: 국소 미세환경과 원격 골전이 신호의 분리</text>
  </g>

  <!-- Panel 3: Complex Stoichiometry -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. 복합체 화학량론 (Complex Stoichiometry)</text>

    <!-- Schematic: Heterodimer Complex -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#9333ea" stroke-width="1"/>
    <!-- Heterodimer schematic -->
    <circle cx="50" cy="79" r="16" fill="#38bdf8" fill-opacity="0.7"/>
    <text x="50" y="83" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">IL6</text>
    <line x1="68" y1="79" x2="105" y2="79" stroke="#94a3b8" stroke-width="2" stroke-dasharray="3,2"/>
    <!-- Receptor subunits -->
    <rect x="110" y="58" width="55" height="42" rx="6" fill="#c084fc" fill-opacity="0.3" stroke="#c084fc" stroke-width="1"/>
    <text x="137" y="73" text-anchor="middle" fill="#c084fc" font-size="8.5" font-weight="700">IL6R</text>
    <text x="137" y="90" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">IL6ST</text>
    <text x="195" y="72" fill="#cbd5e1" font-size="9" font-weight="700">기하평균 공식</text>
    <text x="195" y="88" fill="#fbbf24" font-size="9">R = √(R₁ · R₂)</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 다중 소단위체 필수성 (Subunit Law):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">하나라도 결손 시 결합 불가: S = min(R₁, R₂) 또는 GM</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 핵심 이종복합체 (Key Heterodimers):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">TGFB1 (TGFBR1+TGFBR2), IL6 (IL6R+IL6ST)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 유전자 과발현 오류 배제:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">공동수용체 부재 시 위양성 상호작용 필터링</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#9333ea" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">원리: 생물학적 실제 결합 구조를 반영한 엄밀 계산</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_liana_consensus_and_permutation_testing(filepath: str):
    """
    2. LIANA Consensus Framework & Cell-Label Permutation Testing
    5 DB resources -> Consensus rank aggregation -> Permutation null testing.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad02" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">LIANA 합의 프레임워크 &amp; 순열 검정을 통한 신호 유의성 검증</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">5대 L-R 데이터베이스 통합, 로버스트 순위 집합(RRA), 및 세포 라벨 순열 검정(Permutation Test)</text>

  <!-- Panel 1: DB Harmonization -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 이종 L-R 리소스 통합 (DB Resources)</text>

    <!-- Schematic: 5 DBs flowing to LIANA -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- DB blocks -->
    <rect x="25" y="50" width="60" height="18" rx="3" fill="#0284c7" fill-opacity="0.4"/>
    <text x="55" y="63" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">CellPhoneDB</text>
    <rect x="95" y="50" width="55" height="18" rx="3" fill="#059669" fill-opacity="0.4"/>
    <text x="122" y="63" text-anchor="middle" fill="#34d399" font-size="8" font-weight="700">CellChat</text>
    <rect x="160" y="50" width="55" height="18" rx="3" fill="#d97706" fill-opacity="0.4"/>
    <text x="187" y="63" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">SignalR</text>
    <rect x="225" y="50" width="50" height="18" rx="3" fill="#9333ea" fill-opacity="0.4"/>
    <text x="250" y="63" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">Connectome</text>

    <!-- Harmonized Union Box -->
    <rect x="45" y="78" width="210" height="26" rx="4" fill="#38bdf8" fill-opacity="0.15" stroke="#38bdf8" stroke-width="1"/>
    <text x="150" y="95" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">LIANA OmniPath 표준화 통합 L-R DB</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 데이터베이스 불일치 해소:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">DB마다 다른 복합체 표기법 및 명명법 단일화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 유전자 심볼 및 복합체 매핑:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">HGNC 공식 심볼 기준 다중 소단위체 자동 분해</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 단일 DB 편향성(Bias) 극복:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">문헌 수집 편향 제거 및 상호작용 검증 확대</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">핵심: 신뢰성 높은 통합 L-R 카탈로그 확보</text>
  </g>

  <!-- Panel 2: Consensus Rank Aggregation -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#fbbf24" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#fbbf24" font-size="13" font-weight="700">2. 합의 순위 집합 (Consensus Rank)</text>

    <!-- Schematic: Robust Rank Aggregation (RRA) -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <!-- Rank matrix visual -->
    <g transform="translate(25, 50)">
      <rect x="0" y="0" width="85" height="16" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
      <text x="5" y="12" fill="#cbd5e1" font-size="8">NATMI: #12</text>
      <rect x="0" y="20" width="85" height="16" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
      <text x="5" y="32" fill="#cbd5e1" font-size="8">CellChat: #3</text>
      <rect x="0" y="40" width="85" height="16" fill="#1e293b" stroke="#64748b" stroke-width="0.8"/>
      <text x="5" y="52" fill="#cbd5e1" font-size="8">SingleCellR: #5</text>
    </g>
    <!-- Arrow -->
    <path d="M125,79 L155,79" stroke="#fbbf24" stroke-width="2" marker-end="url(#arrowEmerald)"/>
    <!-- Consensus box -->
    <rect x="165" y="52" width="110" height="52" rx="5" fill="#fbbf24" fill-opacity="0.15" stroke="#fbbf24" stroke-width="1.2"/>
    <text x="220" y="70" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">Robust Rank (RRA)</text>
    <text x="220" y="85" text-anchor="middle" fill="#cbd5e1" font-size="8">순서통계량 집합</text>
    <text x="220" y="98" text-anchor="middle" fill="#34d399" font-size="8" font-weight="600">p_aggregate &lt; 0.001</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• Robust Rank Aggregation(RRA) 집계:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">순위의 일관성을 확률론적 순서통계량으로 평가</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 상위 1% 공통 상호작용 선별:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">단일 알고리즘의 편향된 이상치(Outlier) 순위 제거</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 재현성(Reproducibility) 극대화:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">다중 방법 합의로 방법 간 편차 완화 (개념)</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#fbbf24" font-size="10" font-weight="600">원리: 방법론 간 합의(Consensus)를 통한 위양성 억제</text>
  </g>

  <!-- Panel 3: Permutation Testing -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">3. 세포 라벨 순열 검정 (Permutation)</text>

    <!-- Schematic: Permutation Null curve -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    <!-- Gaussian null curve -->
    <path d="M35,102 Q85,102 110,62 Q125,48 135,48 Q145,48 160,62 Q185,102 235,102" fill="none" stroke="#64748b" stroke-width="1.5"/>
    <!-- Observed threshold line -->
    <line x1="195" y1="46" x2="195" y2="105" stroke="#f43f5e" stroke-width="2" stroke-dasharray="3,2"/>
    <text x="195" y="42" text-anchor="middle" fill="#f43f5e" font-size="8" font-weight="700">S_obs</text>
    <!-- Shaded p-value tail -->
    <path d="M195,78 Q210,102 235,102 L195,102 Z" fill="#f43f5e" fill-opacity="0.5"/>
    <text x="245" y="80" fill="#f43f5e" font-size="8" font-weight="700">p &lt; 0.001</text>
    <text x="70" y="65" fill="#94a3b8" font-size="8">영분포 (B=1000)</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 무작위 셔플링 (B = 1,000회):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">세포 유형 라벨을 섞어 귀무 가설 영분포 구축</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 경험적 P-value 산출 공식:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">P_perm = (Σ 𝕀(S_null ≥ S_obs) + 1) / (B + 1)</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• Benjamini-Hochberg FDR 보정:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">다중가설 검정 보정을 통해 위발견율(FDR &lt; 0.05) 제어</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">검증: 통계적 유의성이 검증된 신호 상호작용만 잔류</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_cellchat_communication_probability_and_centrality(filepath: str):
    """
    3. CellChat Communication Probability & Network Centrality
    Mass action law + Hill cooperativity -> 3D tensor -> 2D network centrality.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad03" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">CellChat 신호 전달 확률 모델 &amp; 네트워크 중심성 분석</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">질량작용 법칙 기반 Hill 협동성 함수, 3차원 신호 텐서, 및 2D 세포 통신 허브 역할 규명</text>

  <!-- Panel 1: Mass Action & Hill Function -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 질량작용 &amp; Hill 함수 (Probability)</text>

    <!-- Schematic: Hill Curve -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Coordinate axes -->
    <line x1="35" y1="105" x2="125" y2="105" stroke="#64748b" stroke-width="1"/>
    <line x1="35" y1="105" x2="35" y2="52" stroke="#64748b" stroke-width="1"/>
    <!-- Hill sigmoid curve -->
    <path d="M35,103 Q65,102 75,85 Q85,68 120,58" fill="none" stroke="#38bdf8" stroke-width="2"/>
    <text x="75" y="102" fill="#fbbf24" font-size="7.5">Kh=0.5</text>
    <!-- Formula box -->
    <rect x="135" y="48" width="140" height="60" rx="4" fill="#0284c7" fill-opacity="0.15" stroke="#38bdf8" stroke-width="0.8"/>
    <text x="205" y="66" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">P = (L · R) / (Kh + L · R)</text>
    <text x="205" y="82" text-anchor="middle" fill="#cbd5e1" font-size="8">L = (∏ L_i)^(1/n)</text>
    <text x="205" y="97" text-anchor="middle" fill="#cbd5e1" font-size="8">R = (∏ R_j)^(1/m)</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 비선형 포화 역동학 (Hill Kinetics):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">리간드-수용체 농도 증가 시 점근적 포화 상태 반영</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 25% 절삭 평균 (Trimmed Mean):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">단일세포 드롭아웃 및 이상치에 의한 과대추정 차단</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 보조인자(Co-factors) 효과 통합:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">촉진/억제 보조 수용체 가중치 가산/감산 계산</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">장점: 생화학적 수용체 결합 동역학 완벽 모사</text>
  </g>

  <!-- Panel 2: 3D Communication Tensor -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 3D 신호 전달 텐서 (3D Tensor)</text>

    <!-- Schematic: 3D Tensor Box -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <!-- Isometric Cube / Tensor slices -->
    <g transform="translate(60, 48)">
      <!-- Back slice -->
      <polygon points="30,10 90,10 80,45 20,45" fill="#059669" fill-opacity="0.2" stroke="#34d399" stroke-width="0.8"/>
      <!-- Mid slice -->
      <polygon points="40,18 100,18 90,53 30,53" fill="#059669" fill-opacity="0.3" stroke="#34d399" stroke-width="0.8"/>
      <!-- Front slice -->
      <polygon points="50,26 110,26 100,61 40,61" fill="#059669" fill-opacity="0.4" stroke="#34d399" stroke-width="1.2"/>
      <text x="75" y="47" text-anchor="middle" fill="#ffffff" font-size="8.5" font-weight="700">P_ijk</text>
    </g>
    <text x="185" y="65" fill="#34d399" font-size="9" font-weight="700">Tensor Dimensions</text>
    <text x="185" y="80" fill="#cbd5e1" font-size="8">i: 발신 세포 (K)</text>
    <text x="185" y="93" fill="#cbd5e1" font-size="8">j: 수신 세포 (K)</text>
    <text x="185" y="106" fill="#cbd5e1" font-size="8">k: L-R 경로 (M)</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 경로(Pathway) 수준 신호 통합:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일 L-R들을 생물학적 경로(TGFb, WNT 등)로 축약</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 조건 간 비교 분석 (Differential Flux):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">호르몬 반응성 vs 거세저항 전립선암(CRPC) 텐서 비교</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 다변량 행렬 인수분해 (Tensor NMF):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">공통 신호 패턴 모듈 자동 추출 및 클러스터링</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">구조: 세포군 간 전방위 통신을 3차원 공간에 완벽 저장</text>
  </g>

  <!-- Panel 3: Centrality Roles -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. 2D 네트워크 중심성 (Centrality)</text>

    <!-- Schematic: 4 Roles 2x2 -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#9333ea" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Sender -->
      <rect x="0" y="0" width="120" height="28" rx="4" fill="#0284c7" fill-opacity="0.3"/>
      <text x="8" y="14" fill="#38bdf8" font-size="8.5" font-weight="700">발신자 (Sender)</text>
      <text x="8" y="24" fill="#94a3b8" font-size="7">Out-degree 중심성</text>
      <!-- Receiver -->
      <rect x="130" y="0" width="120" height="28" rx="4" fill="#059669" fill-opacity="0.3"/>
      <text x="138" y="14" fill="#34d399" font-size="8.5" font-weight="700">수신자 (Receiver)</text>
      <text x="138" y="24" fill="#94a3b8" font-size="7">In-degree 중심성</text>
      <!-- Mediator -->
      <rect x="0" y="32" width="120" height="28" rx="4" fill="#d97706" fill-opacity="0.3"/>
      <text x="8" y="46" fill="#fbbf24" font-size="8.5" font-weight="700">중계자 (Mediator)</text>
      <text x="8" y="56" fill="#94a3b8" font-size="7">Flow betweenness 중심성</text>
      <!-- Influencer -->
      <rect x="130" y="32" width="120" height="28" rx="4" fill="#9333ea" fill-opacity="0.3"/>
      <text x="138" y="46" fill="#c084fc" font-size="8.5" font-weight="700">영향력자 (Influencer)</text>
      <text x="138" y="56" fill="#94a3b8" font-size="7">정보(Information) 중심성</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 발신자 vs 수신자 특이성:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">CAF(지배적 발신자), 내피세포/T세포(주요 수신자)</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 매개자(Mediator) 허브 규명:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">M2 대식세포가 암세포와 혈관/T세포 간 중계 수행</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 2D 투영 시각화 (Centrality Plot):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">Out-degree vs In-degree 산점도로 신호 흐름 한눈 파악</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#9333ea" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">응용: 네트워크 급소(Vulnerability) 타깃 치료 표적 발굴</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_nichenet_ligand_to_target_causal_network(filepath: str):
    """
    4. NicheNet Intracellular Ligand-to-Target Causal Network
    Prior signaling-GRN graph -> Random Walk with Restart (RWR) -> Target gene regulatory potential W.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad04" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrowRose" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">NicheNet 세포 내 리간드-표적 유전자 인과 조절 네트워크</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">3단계 사전 지식 그래프(L-R, 신호전달 PPI, GRN-TF), 랜덤 워크(RWR), 및 표적 유전자 귀속</text>

  <!-- Panel 1: 3-Layer Graph -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 3단계 통합 사전 지식망 (Prior Graph)</text>

    <!-- Schematic: 3 Layers -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Layer 1 -->
    <rect x="25" y="48" width="70" height="18" rx="3" fill="#0284c7" fill-opacity="0.3"/>
    <text x="60" y="61" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">L-R 결합층</text>
    <!-- Arrow -->
    <line x1="95" y1="57" x2="115" y2="57" stroke="#38bdf8" stroke-width="1.5"/>
    <!-- Layer 2 -->
    <rect x="115" y="48" width="70" height="18" rx="3" fill="#d97706" fill-opacity="0.3"/>
    <text x="150" y="61" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">신호 PPI층</text>
    <!-- Arrow -->
    <line x1="185" y1="57" x2="205" y2="57" stroke="#fbbf24" stroke-width="1.5"/>
    <!-- Layer 3 -->
    <rect x="205" y="48" width="70" height="18" rx="3" fill="#9333ea" fill-opacity="0.3"/>
    <text x="240" y="61" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">GRN-TF층</text>
    <!-- Subtext inside schematic -->
    <text x="150" y="88" text-anchor="middle" fill="#cbd5e1" font-size="8.5">리간드 ➔ 인산화 캐스케이드 ➔ TF (AR, FOXA1)</text>
    <text x="150" y="102" text-anchor="middle" fill="#34d399" font-size="8" font-weight="600">Track 01 유전자 조절망(GRN) 직접 연계</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 단순 발현 상관관계를 넘는 인과성:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단순 L-R 존재 여부가 아닌 수신세포 반응 유전자 추적</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 방대한 분자 생물학 데이터 통합:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">STRING, SignaLink, ChIP-seq, RegNetwork 융합</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 전사인자(TF) 결합 부위 매핑:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">프로모터/인핸서 모티프 기반 최종 표적 유전자 연결</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">혁신: 세포 외 신호와 핵 내 유전자 발현 변화의 결합</text>
  </g>

  <!-- Panel 2: Random Walk with Restart -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#fbbf24" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#fbbf24" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#fbbf24" font-size="13" font-weight="700">2. 재시작 랜덤 워크 (RWR Potential)</text>

    <!-- Schematic: RWR Equation & Matrix W -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <!-- Matrix W -->
    <rect x="30" y="52" width="55" height="54" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
    <text x="57" y="75" text-anchor="middle" fill="#fbbf24" font-size="11" font-weight="700">W</text>
    <text x="57" y="92" text-anchor="middle" fill="#94a3b8" font-size="7">L × Target</text>
    <!-- RWR Flow -->
    <path d="M95,79 L125,79" stroke="#fbbf24" stroke-width="2" marker-end="url(#arrowRose)"/>
    <g transform="translate(135, 52)">
      <text x="0" y="16" fill="#fbbf24" font-size="8.5" font-weight="700">p_(t+1) = (1-r) M p_t + r e_L</text>
      <text x="0" y="32" fill="#cbd5e1" font-size="8">r = 재시작 확률 (0.50)</text>
      <text x="0" y="46" fill="#34d399" font-size="8">정상상태(Steady-state) 수렴</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 조절 잠재력 행렬 (Regulatory Potential W):</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">리간드 L이 표적 유전자 G에 미치는 사전 확률 정량화</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 다중 분기 경로의 전파 감쇠 반영:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">네트워크 거리 및 노드 차수에 따른 신호 확산 계산</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 허브 단백질 편향 제거 (Degree Normalization):</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">과도한 연결을 가진 단백질에 의한 인위적 점수 왜곡 방지</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#d97706" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#fbbf24" font-size="10" font-weight="600">수학: 마르코프 연쇄(Markov Chain) 기반 신호 유전</text>
  </g>

  <!-- Panel 3: Ligand Activity Scoring -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">3. 리간드 활성도 평가 (Ligand Activity)</text>

    <!-- Schematic: Pearson & AUROC Bar -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    <!-- Top Ranked Ligand Bars -->
    <g transform="translate(25, 50)">
      <text x="0" y="12" fill="#cbd5e1" font-size="8" font-weight="600">TGFB1</text>
      <rect x="45" y="3" width="128" height="11" rx="2" fill="#f43f5e"/>
      <text x="183" y="12" fill="#f43f5e" font-size="8" font-weight="700">AUROC 0.96</text>

      <text x="0" y="29" fill="#cbd5e1" font-size="8" font-weight="600">CXCL12</text>
      <rect x="45" y="20" width="122" height="11" rx="2" fill="#fbbf24"/>
      <text x="177" y="29" fill="#fbbf24" font-size="8" font-weight="700">AUROC 0.92</text>

      <text x="0" y="46" fill="#cbd5e1" font-size="8" font-weight="600">TNFSF11</text>
      <rect x="45" y="37" width="131" height="11" rx="2" fill="#34d399"/>
      <text x="186" y="46" fill="#34d399" font-size="8" font-weight="700">AUROC 0.99</text>

      <text x="0" y="60" fill="#94a3b8" font-size="6">※ toy 합성 데이터 산출값 (pipelines/results_pipeline_02)</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 수신 세포 DEG 예측력 검증:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">CD8+ T 세포 고갈 유전자(PDCD1, TOX) 발현 예측</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• AUROC &amp; Pearson 상관분석:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">조절 잠재력과 실제 발현 차이 간의 일치도 검정</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 핵심 업스트림 드라이버 도출:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">TGFB1(EMT/고갈 유도), CXCL12(기질 장벽 형성)</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">결과: 표현형 변화를 직접 촉발한 원인 리간드 규명</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_spatial_misty_multiview_modeling(filepath: str):
    """
    5. MISTy Multiview Spatial Decomposition & Contextual Modeling
    Intrinsic (0µm) + Juxtaview (<15µm) + Paraview (100-250µm) -> Random forest ensemble -> Spatial Gain.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad05" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">MISTy 멀티뷰 공간 분해 &amp; 맥락적 세포 신호 모델링</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">고유 발현(Intrinsic), 근접 접촉(Juxtaview), 원격 확산(Paraview)의 기계학습 기반 앙상블 기여도 분해</text>

  <!-- Panel 1: Multiview Spatial Decomposition -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 3대 공간 맥락 분해 (Multiview)</text>

    <!-- Schematic: Concentric Circles -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- Circles -->
    <circle cx="70" cy="79" r="30" fill="#0284c7" fill-opacity="0.1" stroke="#38bdf8" stroke-width="1" stroke-dasharray="2,2"/>
    <circle cx="70" cy="79" r="18" fill="#fbbf24" fill-opacity="0.2" stroke="#fbbf24" stroke-width="1.2"/>
    <circle cx="70" cy="79" r="8" fill="#f43f5e" fill-opacity="0.8"/>
    <!-- Legend -->
    <g transform="translate(120, 52)">
      <circle cx="5" cy="6" r="4" fill="#f43f5e"/>
      <text x="15" y="9" fill="#cbd5e1" font-size="8.5">Intrinsic (0 µm): 세포 고유 GRN</text>
      <circle cx="5" cy="22" r="4" fill="#fbbf24"/>
      <text x="15" y="25" fill="#cbd5e1" font-size="8.5">Juxtaview (&lt;15 µm): 직접 접촉</text>
      <circle cx="5" cy="38" r="4" fill="#38bdf8"/>
      <text x="15" y="41" fill="#cbd5e1" font-size="8.5">Paraview (100~250 µm): 확산장</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 거리 계층적 뷰 생성:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">단일세포의 발현량을 물리적 반경에 따라 분리</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 지수 감쇠 커널 (Exponential Decay):</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">w(d) = exp(-d / σ) 가중치를 적용한 주변부 농도 합산</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 세포 유형 밀도 vs 리간드 농도:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">이웃 세포의 구성비와 분비 물질의 공간적 중첩 포착</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">공간 분해: 세포 내부 요인과 외부 환경 자극의 분리</text>
  </g>

  <!-- Panel 2: Random Forest Ensemble -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">2. 멀티뷰 앙상블 학습 (Ensemble RF)</text>

    <!-- Schematic: Stacked Views into Model -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 52)">
      <rect x="0" y="0" width="60" height="16" fill="#f43f5e" fill-opacity="0.4"/>
      <text x="30" y="11" text-anchor="middle" fill="#ffffff" font-size="7.5">X_intra</text>
      <rect x="0" y="19" width="60" height="16" fill="#fbbf24" fill-opacity="0.4"/>
      <text x="30" y="30" text-anchor="middle" fill="#ffffff" font-size="7.5">X_juxta</text>
      <rect x="0" y="38" width="60" height="16" fill="#38bdf8" fill-opacity="0.4"/>
      <text x="30" y="49" text-anchor="middle" fill="#ffffff" font-size="7.5">X_para</text>
    </g>
    <!-- Model Box -->
    <rect x="105" y="52" width="95" height="54" rx="5" fill="#059669" fill-opacity="0.2" stroke="#34d399" stroke-width="1.2"/>
    <text x="152" y="73" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">Random Forest</text>
    <text x="152" y="88" text-anchor="middle" fill="#cbd5e1" font-size="7.5">Meta-regressor</text>
    <text x="235" y="75" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="700">y_gene</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 2단계 스태킹(Stacking) 구조:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">각 공간 뷰별 모델 학습 후 메타 회귀기로 융합</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 비선형 특징 상호작용 포착:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">단순 선형 결합이 아닌 임계값 효과 및 포화 현상 반영</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 다변량 유전자 동시 예측:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">전체 전사체 패널에 대한 공간 기여도 병렬 계산</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">머신러닝: 해석 가능한 공간 유전자 조절 역학 규명</text>
  </g>

  <!-- Panel 3: Spatial Gain Attribution -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#c084fc" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#c084fc" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#c084fc" font-size="13" font-weight="700">3. 공간 이득 &amp; 기여도 귀속 (Spatial Gain)</text>

    <!-- Schematic: R2 Gain Partition -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#9333ea" stroke-width="1"/>
    <!-- Stacked Bar for R2 -->
    <g transform="translate(30, 60)">
      <rect x="0" y="0" width="100" height="18" fill="#f43f5e" fill-opacity="0.8"/>
      <text x="50" y="13" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">Intrinsic (예시)</text>
      <rect x="100" y="0" width="50" height="18" fill="#fbbf24" fill-opacity="0.8"/>
      <text x="125" y="13" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">Juxta</text>
      <rect x="150" y="0" width="80" height="18" fill="#38bdf8" fill-opacity="0.8"/>
      <text x="190" y="13" text-anchor="middle" fill="#ffffff" font-size="8" font-weight="700">Para</text>
    </g>
    <text x="150" y="98" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">Spatial Gain = R²(Full) - R²(Intrinsic) (예시 도식)</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 이득(Spatial Gain) 지표:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">공간 정보를 포함했을 때 설명 분산(R²)의 유의미한 증가</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 기질 장벽 vs 면역 억제 규명:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">Paraview의 높은 기여로 myCAF 유래 원격 신호 검증</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 공간 치료 타깃 선정 근거:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">주변 환경 의존성이 높은 종양 악성화 유전자 우선순위화</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#9333ea" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#c084fc" font-size="10" font-weight="600">결론: 공간적 맥락이 유전자 발현에 미치는 영향 완전 정량화</text>
  </g>
  <text x="500" y="353" text-anchor="middle" fill="#64748b" font-size="9">※ 개념도: 막대의 비율(%)은 설명용 예시이며 실측 결과가 아님</text>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_prostate_bone_metastasis_vicious_cycle(filepath: str):
    """
    6. Prostate Cancer Bone Metastasis Vicious Cycle & Immunosuppression
    Tumor - Osteoclast - Osteoblast - myCAF - CD8 T 5-way vicious cycle and immunosuppressive circuit in mCRPC.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad06" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrowCycleRose" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrowCycleAmber" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#fbbf24"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="30" text-anchor="middle" fill="#38bdf8" font-size="18" font-weight="700">전립선암 골전이 악순환 회로(Vicious Cycle) &amp; 면역억제 신호망</text>
  <text x="500" y="49" text-anchor="middle" fill="#94a3b8" font-size="11.5">mCRPC 골전이 틈새 내 암세포-파골세포-골아세포-CAF-CD8 T 5자 피드포워드 증폭 및 면역 회피</text>

  <!-- Panel 1: 5-Way Cellular Architecture -->
  <g transform="translate(25, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="13" font-weight="700">1. 골미세환경 5대 핵심 세포 (Niche)</text>

    <!-- Schematic: Pentagonal Cell Arrangement -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <!-- 5 Cells -->
    <circle cx="50" cy="65" r="13" fill="#f43f5e" fill-opacity="0.8"/>
    <text x="50" y="69" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">PCa Tumor</text>
    <circle cx="110" cy="58" r="12" fill="#fbbf24" fill-opacity="0.8"/>
    <text x="110" y="62" text-anchor="middle" fill="#ffffff" font-size="7" font-weight="700">Osteoblast</text>
    <circle cx="170" cy="58" r="12" fill="#c084fc" fill-opacity="0.8"/>
    <text x="170" y="62" text-anchor="middle" fill="#ffffff" font-size="7" font-weight="700">Osteoclast</text>
    <circle cx="230" cy="65" r="12" fill="#38bdf8" fill-opacity="0.8"/>
    <text x="230" y="69" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">myCAF</text>
    <circle cx="140" cy="94" r="11" fill="#34d399" fill-opacity="0.8"/>
    <text x="140" y="98" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">CD8 T</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• 골 친화성(Osteotropism) 특이 틈새:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">mCRPC 환자의 80% 이상에서 척추 및 골반 골전이 발생</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 골경화성(Osteosclerotic) 전이 특징:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">파골세포 흡수와 비정상적 골아세포 골생성이 동시 폭발</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 면역 사막(Immune Desert) 환경:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">골수 내 세포독성 T세포의 접근 차단 및 조기 고갈</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#0284c7" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">병리: 파골세포-골아세포 상호작용의 종양 매개 교란</text>
  </g>

  <!-- Panel 2: Molecular Vicious Cycle -->
  <g transform="translate(350, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="13" font-weight="700">2. 악순환 분자 피드포워드 루프</text>

    <!-- Schematic: Vicious Cycle Flow -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    <!-- Step 1 -> Step 2 -> Step 3 Flow -->
    <g transform="translate(25, 48)">
      <rect x="0" y="0" width="60" height="26" rx="3" fill="#f43f5e" fill-opacity="0.3"/>
      <text x="30" y="11" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">암세포</text>
      <text x="30" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6.5">PTHrP/IL6</text>

      <path d="M62,13 L88,13" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrowCycleRose)"/>

      <rect x="90" y="0" width="70" height="26" rx="3" fill="#fbbf24" fill-opacity="0.3"/>
      <text x="125" y="11" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">골아세포</text>
      <text x="125" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6.5">RANKL 분비</text>

      <path d="M162,13 L188,13" stroke="#fbbf24" stroke-width="1.5" marker-end="url(#arrowCycleAmber)"/>

      <rect x="190" y="0" width="60" height="26" rx="3" fill="#c084fc" fill-opacity="0.3"/>
      <text x="220" y="11" text-anchor="middle" fill="#c084fc" font-size="7" font-weight="700">파골세포</text>
      <text x="220" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6.5">골 기질 흡수</text>
    </g>
    <!-- Loopback arrow -->
    <path d="M225,82 Q125,115 45,82" fill="none" stroke="#34d399" stroke-width="1.5" stroke-dasharray="3,2"/>
    <text x="135" y="104" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">기질 내 TGF-β &amp; IGF-1 방출 ➔ 종양 EMT 촉진</text>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• PTHrP ➔ RANKL ➔ RANK 축:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">암세포가 골아세포를 조종해 파골세포 분화 촉진</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• 골 흡수성 TGF-β 방출 피드백:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">뼈 분해 과정에서 방출된 TGF-β가 암세포 증식 증폭</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 엔도텔린-1(ET-1) 매개 골경화증:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">동시에 골아세포 증식을 촉발하여 비정상 골 형성</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#be123c" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="10" font-weight="600">치료 타깃: Denosumab (anti-RANKL) + Galunisertib</text>
  </g>

  <!-- Panel 3: Stromal & Immune Blockade -->
  <g transform="translate(675, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.5"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#34d399" fill-opacity="0.18"/>
    <text x="150" y="21" text-anchor="middle" fill="#34d399" font-size="13" font-weight="700">3. 기질 차단 및 면역 고갈 회피</text>

    <!-- Schematic: CAF Barrier & Checkpoints -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <!-- CAF Barrier -->
    <g transform="translate(25, 52)">
      <rect x="0" y="0" width="85" height="54" rx="4" fill="#38bdf8" fill-opacity="0.2" stroke="#38bdf8" stroke-width="1"/>
      <text x="42" y="18" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">myCAF 장벽</text>
      <text x="42" y="32" text-anchor="middle" fill="#cbd5e1" font-size="7">CXCL12 과발현</text>
      <text x="42" y="46" text-anchor="middle" fill="#94a3b8" font-size="7">T세포 침투 차단</text>

      <!-- Immune Checkpoints -->
      <rect x="100" y="0" width="145" height="54" rx="4" fill="#f43f5e" fill-opacity="0.2" stroke="#f43f5e" stroke-width="1"/>
      <text x="172" y="16" text-anchor="middle" fill="#f43f5e" font-size="8" font-weight="700">면역 체크포인트 억제</text>
      <text x="172" y="30" text-anchor="middle" fill="#fbbf24" font-size="7.5">PD-L1 : PD-1 | CD47 : SIRPA</text>
      <text x="172" y="45" text-anchor="middle" fill="#34d399" font-size="7">CD8 T: TOX+, HAVCR2+, LAG3+</text>
    </g>

    <!-- Bullets -->
    <text x="15" y="132" fill="#cbd5e1" font-size="11" font-weight="600">• myCAF의 면역 세포 물리적 배제:</text>
    <text x="25" y="148" fill="#94a3b8" font-size="10">치밀한 콜라겐 기질 및 CXCL12 코팅으로 T세포 격리</text>

    <text x="15" y="168" fill="#cbd5e1" font-size="11" font-weight="600">• CD47-SIRPA 대식세포 탐식 회피:</text>
    <text x="25" y="184" fill="#94a3b8" font-size="10">"Don't eat me" 신호 전달을 통한 대식세포 탐식 무력화</text>

    <text x="15" y="204" fill="#cbd5e1" font-size="11" font-weight="600">• 복합 면역 요법 저항성 극복 전략:</text>
    <text x="25" y="220" fill="#94a3b8" font-size="10">TGF-β 차단 + 항-PD-1/PD-L1 병용으로 면역 재활성화</text>

    <!-- Badge -->
    <rect x="15" y="242" width="270" height="22" rx="4" fill="#059669" fill-opacity="0.25"/>
    <text x="150" y="257" text-anchor="middle" fill="#34d399" font-size="10" font-weight="600">임상 의의: 골전이 틈새 악순환 차단을 통한 mCRPC 치료</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_tme_cellular_ecosystem_and_signaling_types.svg", make_01_tme_cellular_ecosystem_and_signaling_types),
        ("02_liana_consensus_and_permutation_testing.svg", make_02_liana_consensus_and_permutation_testing),
        ("03_cellchat_communication_probability_and_centrality.svg", make_03_cellchat_communication_probability_and_centrality),
        ("04_nichenet_ligand_to_target_causal_network.svg", make_04_nichenet_ligand_to_target_causal_network),
        ("05_spatial_misty_multiview_modeling.svg", make_05_spatial_misty_multiview_modeling),
        ("06_prostate_bone_metastasis_vicious_cycle.svg", make_06_prostate_bone_metastasis_vicious_cycle),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 03 TME & Cell Communication publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
tracks/07_count_statistics_abundance/scripts/generate_statistics_diagrams.py
=============================================================================
Generates 6 publication-grade, high-resolution vector SVG diagrams for Track 07:
Count Statistics, Pseudoreplication & Differential Abundance in Single-Cell Biology.

Diagrams:
1. 01_poisson_vs_negative_binomial_overdispersion.svg
2. 02_pseudoreplication_type1_error_catastrophe.svg
3. 03_pseudobulk_aggregation_and_dispersion_shrinkage.svg
4. 04_glmm_hierarchical_random_effects_architecture.svg
5. 05_milo_knn_neighborhood_differential_abundance.svg
6. 06_compositional_simplex_and_coda_regression.svg
"""

import os
import xml.etree.ElementTree as ET


def make_01_poisson_vs_negative_binomial_overdispersion(filepath: str):
    """
    1. Poisson vs Negative Binomial Overdispersion & Bursting Kinetics
    Var = Mean vs Var = Mean + alpha*Mean^2, transcriptional bursting, and empirical mean-variance geometry.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="curveGrad01" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="50%" stop-color="#fbbf24"/>
      <stop offset="100%" stop-color="#f43f5e"/>
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
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">포아송(Poisson) vs 음이항(Negative Binomial) 과대산포와 전사 파열 기하학</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">평균-분산 불일치(Var &gt; Mean) ➔ 전사 파열(Transcriptional Bursting) ➔ 감마-포아송 계층 모델(Var = μ + α·μ²)</text>

  <!-- Panel 1: Poisson Model & Limitations -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 포아송 모델의 한계 (Poisson Model)</text>

    <!-- Schematic: Poisson PMF Curve -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 48)">
      <!-- Axes -->
      <line x1="15" y1="58" x2="245" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="15" y1="58" x2="15" y2="8" stroke="#475569" stroke-width="1"/>
      <!-- Narrow Poisson PMF bell -->
      <path d="M 25 58 Q 65 57 85 40 Q 115 5 130 5 Q 145 5 170 42 Q 195 58 235 58" fill="none" stroke="#38bdf8" stroke-width="2"/>
      <line x1="130" y1="58" x2="130" y2="8" stroke="#38bdf8" stroke-width="1" stroke-dasharray="2,2"/>
      <circle cx="130" cy="5" r="3.5" fill="#38bdf8"/>
      <text x="135" y="16" fill="#38bdf8" font-size="7.5" font-weight="700">μ = 10 (Mean)</text>
      <text x="200" y="25" fill="#94a3b8" font-size="7">Var = μ = 10</text>
      <text x="200" y="37" fill="#fbbf24" font-size="7">Fano Factor = 1</text>
      <text x="50" y="66" fill="#64748b" font-size="6.5">Count (Y)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 기술적 샷 노이즈(Shot Noise)의 가정:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">P(Y=k) = (λ^k · e^(-λ)) / k! (무작위 분자 포획)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 등분산성 제약 (Equidispersion):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">E[Y] = Var(Y) = λ (분산이 평균과 동일해야 함)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 실제 단일세포 전사체의 실패:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">생물학적 이질성으로 인해 분산이 평균의 수십 배 상회</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">위험: 포아송 가정 시 p-값 허위 축소 및 위양성 폭발</text>
  </g>

  <!-- Panel 2: Negative Binomial & Transcriptional Bursting -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 음이항 분포와 전사 파열 동역학</text>

    <!-- Schematic: Bursting & Overdispersed Curve -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Promoter state diagram -->
      <rect x="5" y="5" width="48" height="22" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="29" y="19" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">OFF</text>
      <path d="M 53 10 Q 68 5 83 10" stroke="#fbbf24" stroke-width="1.2" fill="none" marker-end="url(#arrAmber01)"/>
      <text x="68" y="5" text-anchor="middle" fill="#fbbf24" font-size="6.5">k_on</text>
      <path d="M 83 22 Q 68 27 53 22" stroke="#94a3b8" stroke-width="1.2" fill="none" marker-end="url(#arrCyan01)"/>
      <text x="68" y="32" text-anchor="middle" fill="#94a3b8" font-size="6.5">k_off</text>
      <rect x="83" y="5" width="48" height="22" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
      <text x="107" y="19" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">ON</text>

      <!-- Broad curve comparison -->
      <line x1="150" y1="58" x2="265" y2="58" stroke="#475569" stroke-width="1"/>
      <!-- Poisson curve (dashed) -->
      <path d="M 160 58 Q 190 20 205 10 Q 220 20 240 58" fill="none" stroke="#38bdf8" stroke-width="1.2" stroke-dasharray="2,2"/>
      <!-- Negative Binomial (broad with heavy tail and zero peak) -->
      <path d="M 152 20 Q 165 42 185 45 Q 215 48 260 58" fill="none" stroke="#fbbf24" stroke-width="2"/>
      <text x="180" y="12" fill="#38bdf8" font-size="6.5">Poisson</text>
      <text x="215" y="38" fill="#fbbf24" font-size="6.8" font-weight="700">NB (Overdispersed)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 감마-포아송 혼합 (Gamma-Poisson):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Y ~ NB(μ, α) ⟺ Y|λ ~ Poisson(λ), λ ~ Gamma(θ, μ/θ)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이차항 과대산포 공식 (Variance Equation):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">Var(Y) = μ + α·μ² = μ + μ² / θ (α: 산포, θ: 형상 모수)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 전사 파열 기전 (Transcriptional Bursting):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">프로모터 ON/OFF 전이로 인한 폭발적 mRNA 방출 모형화</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">단일세포 카운트 데이터링의 표준 통계 분포 (edgeR/sctransform)</text>
  </g>

  <!-- Panel 3: Empirical Mean-Variance Geometry -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 단일세포 평균-분산 기하학</text>

    <!-- Schematic: Log(Mean) vs Log(Var) Scatter -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Axes -->
      <line x1="20" y1="58" x2="255" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="58" x2="20" y2="8" stroke="#475569" stroke-width="1"/>
      <text x="20" y="5" fill="#64748b" font-size="6.5">log(Var)</text>
      <text x="240" y="66" fill="#64748b" font-size="6.5">log(Mean)</text>

      <!-- Poisson Line y=x -->
      <line x1="25" y1="53" x2="160" y2="15" stroke="#38bdf8" stroke-width="1.2" stroke-dasharray="3,3"/>
      <text x="95" y="42" fill="#38bdf8" font-size="6.5">Var = Mean (Shot Noise)</text>

      <!-- NB Quadratic Trend Line -->
      <path d="M 25 53 Q 120 40 230 10" fill="none" stroke="#34d399" stroke-width="2"/>
      <text x="165" y="16" fill="#34d399" font-size="6.8" font-weight="700">Var = μ + α·μ²</text>

      <!-- Data Dots -->
      <circle cx="50" cy="46" r="2" fill="#94a3b8"/>
      <circle cx="80" cy="38" r="2" fill="#94a3b8"/>
      <circle cx="110" cy="32" r="2" fill="#94a3b8"/>
      <circle cx="140" cy="24" r="2" fill="#94a3b8"/>
      <circle cx="180" cy="16" r="2.5" fill="#fbbf24"/>
      <!-- HVG outlier dots -->
      <circle cx="195" cy="8" r="3" fill="#f43f5e"/>
      <circle cx="215" cy="6" r="3" fill="#f43f5e"/>
      <text x="215" y="19" fill="#f43f5e" font-size="6.5" font-weight="700">HVGs (AR, SYP)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 생물학적 변동 계수 (BCV):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">BCV = sqrt(α) = sqrt(Var - μ) / μ (순수 생물학적 편차)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 2계층 분산 지배 원리:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">저발현: 기술적 샷 노이즈 지배 / 고발현: 생물학적 변동 지배</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 고변동성 유전자(HVG) 선별:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">음이항 추세선을 유의하게 상회하는 유전자만 생물학적 가치</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">결론: 과대산포 모델링 없이는 차등발현 검정의 유의성 왜곡</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_pseudoreplication_type1_error_catastrophe(filepath: str):
    """
    2. Pseudoreplication Type I Error Catastrophe
    Naive single-cell test (N=10,000 cells, Type I error exploding to 85%) vs Donor-level pseudobulk (N=3 donors, 5%).
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad02" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrRose02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrEmerald02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#f43f5e" font-size="16.5" font-weight="700">유사복제(Pseudoreplication)의 제1종 오류(Type I Error) 폭발 참사</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">개별 세포(N=10,000)를 독립 표본 취급 시 귀무가설 붕괴 vs 공여자 수준(N=3) 의사벌크(Pseudobulk)</text>

  <!-- Panel 1: Naive Single-Cell Test Disaster -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad02)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. 나이브 세포 수준 검정 (N=10,000)</text>

    <!-- Schematic: Exploding Volcano Plot under Pure Null -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#e11d48" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Volcano Axes -->
      <line x1="20" y1="58" x2="245" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="130" y1="58" x2="130" y2="6" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="46" x2="245" y2="46" stroke="#fbbf24" stroke-width="1" stroke-dasharray="2,2"/>
      <text x="210" y="43" fill="#fbbf24" font-size="6">p = 0.05</text>
      <!-- False positive swarm (red dots everywhere above line) -->
      <circle cx="110" cy="35" r="2" fill="#f43f5e"/>
      <circle cx="150" cy="32" r="2" fill="#f43f5e"/>
      <circle cx="95" cy="20" r="2.2" fill="#f43f5e"/>
      <circle cx="165" cy="18" r="2.2" fill="#f43f5e"/>
      <circle cx="80" cy="10" r="2.5" fill="#f43f5e"/>
      <circle cx="180" cy="8" r="2.5" fill="#f43f5e"/>
      <circle cx="120" cy="15" r="2" fill="#f43f5e"/>
      <circle cx="140" cy="12" r="2" fill="#f43f5e"/>
      <text x="130" y="4" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">85% 허위 유의 (False Positives)</text>
      <text x="35" y="20" fill="#f43f5e" font-size="6.5">Pure Null Data!</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 동일 환자 세포 간 상관성 무시:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">급내상관계수(ICC &gt; 0)를 0으로 취급하여 I.I.D. 가정 위반</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 표준오차(SE) 극단적 과소추정:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">SE_naive = σ / sqrt(N_cells) ≪ SE_true = σ_donor / sqrt(N_donors)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• z-score 인위적 10~50배 폭증:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">아무런 생물학적 차이가 없는 난수 유전자도 p &lt; 10^(-20)</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.2" font-weight="600">위양성률(Type I Error) 80~90% 폭발 (Squair et al., 2021)</text>
  </g>

  <!-- Panel 2: Effective Sample Size & Intraclass Correlation -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 유효 표본 크기(N_eff)와 ICC 왜곡</text>

    <!-- Schematic: Type I Error Comparison Bar Chart -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Baseline axis -->
      <line x1="20" y1="56" x2="255" y2="56" stroke="#475569" stroke-width="1"/>
      <!-- Nominal 5% alpha line -->
      <line x1="20" y1="48" x2="255" y2="48" stroke="#38bdf8" stroke-width="1" stroke-dasharray="2,2"/>
      <text x="220" y="45" fill="#38bdf8" font-size="6">목표 α = 0.05</text>

      <!-- Bar 1: Wilcoxon Naive -->
      <rect x="35" y="10" width="40" height="46" rx="2" fill="#f43f5e"/>
      <text x="55" y="24" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">82%</text>
      <text x="55" y="65" text-anchor="middle" fill="#94a3b8" font-size="6.5">Wilcoxon</text>

      <!-- Bar 2: t-test Naive -->
      <rect x="95" y="8" width="40" height="48" rx="2" fill="#f43f5e"/>
      <text x="115" y="22" text-anchor="middle" fill="#ffffff" font-size="7.5" font-weight="700">86%</text>
      <text x="115" y="65" text-anchor="middle" fill="#94a3b8" font-size="6.5">t-test</text>

      <!-- Bar 3: Mixed Model -->
      <rect x="155" y="47" width="40" height="9" rx="2" fill="#34d399"/>
      <text x="175" y="44" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">5.2%</text>
      <text x="175" y="65" text-anchor="middle" fill="#94a3b8" font-size="6.5">GLMM</text>

      <!-- Bar 4: Pseudobulk -->
      <rect x="210" y="47.5" width="40" height="8.5" rx="2" fill="#38bdf8"/>
      <text x="230" y="44" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">4.8%</text>
      <text x="230" y="65" text-anchor="middle" fill="#94a3b8" font-size="6.5">Pseudobulk</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유효 표본 크기 축소 공식 (Kish’s Design Effect):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">N_eff = N_total / [ 1 + (m - 1)·ρ ] (m: 환자당 세포 수)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 급내상관계수(ICC, ρ)의 치명성:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">ρ = 0.03에 불과해도 1,000세포의 실질 정보량은 30개로 급감</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 세포 복제(Cell replicate) ≠ 생물학적 복제:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">세포를 10만 개 늘려도 환자 수가 3명이면 N=3 통계</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">핵심 정리: 표본 수는 세포 방울 수가 아닌 환자 수에 종속</text>
  </g>

  <!-- Panel 3: Donor-Level Pseudobulk Solution -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 공여자 수준 의사벌크 (Pseudobulk)</text>

    <!-- Schematic: Strictly Calibrated Null Volcano Plot -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Volcano Axes -->
      <line x1="20" y1="58" x2="255" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="135" y1="58" x2="135" y2="6" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="20" x2="255" y2="20" stroke="#fbbf24" stroke-width="1" stroke-dasharray="2,2"/>
      <text x="215" y="17" fill="#fbbf24" font-size="6">p = 0.05</text>
      <!-- All null dots sit safely below p=0.05 line, uniform spread -->
      <circle cx="130" cy="52" r="1.8" fill="#94a3b8"/>
      <circle cx="140" cy="48" r="1.8" fill="#94a3b8"/>
      <circle cx="110" cy="45" r="1.8" fill="#94a3b8"/>
      <circle cx="160" cy="42" r="1.8" fill="#94a3b8"/>
      <circle cx="90" cy="38" r="1.8" fill="#94a3b8"/>
      <circle cx="180" cy="36" r="1.8" fill="#94a3b8"/>
      <circle cx="125" cy="30" r="1.8" fill="#94a3b8"/>
      <circle cx="150" cy="25" r="1.8" fill="#94a3b8"/>
      <!-- Only exactly ~5% cross line purely by chance -->
      <circle cx="170" cy="14" r="2" fill="#34d399"/>
      <text x="135" y="66" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">정확히 5% 통제 (Calibrated Null)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 분석 단위 복원 (Unit of Replication):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Y_{g,d} = Σ_{c ∈ d} Y_{g,c} (환자 d 내 세포 카운트 총합)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 환자 간 생물학적 변동 분산 완벽 반영:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">자유도 df = N_donors - p로 통계적 검정력 정규화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 임상 코호트 검증 신뢰도 100%:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">독립적인 외부 코호트에서 재현 가능한 DEG 도출</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">골드 스탠다드: 다중 환자 단일세포 DEG의 타협 없는 기준</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_pseudobulk_aggregation_and_dispersion_shrinkage(filepath: str):
    """
    3. Pseudobulk Aggregation & Empirical Bayes Dispersion Shrinkage
    Single-cell sparse matrix -> Donor-by-celltype summation tensor -> Empirical Bayes dispersion shrinkage.
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
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">의사벌크(Pseudobulk) 텐서 집계 &amp; 경험적 베이즈 산포 수축(Dispersion Shrinkage)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">단일세포 희소 행렬 ➔ 공여자×세포형 요약 ➔ edgeR/DESeq2 경험적 베이즈 사전분포(Shrinkage)</text>

  <!-- Panel 1: Single-Cell 3D Tensor Aggregation -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 단일세포 텐서 집계 (Pseudobulk Tensor)</text>

    <!-- Schematic: Sparse Matrix to 3D Tensor -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Sparse 2D Matrix -->
      <rect x="5" y="8" width="55" height="48" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="32" y="34" text-anchor="middle" fill="#94a3b8" font-size="7">Cells × Genes (85% 0)</text>

      <path d="M 68 32 L 95 32" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan03)"/>
      <text x="81" y="24" text-anchor="middle" fill="#fbbf24" font-size="6.5">Σ_{c ∈ (d,k)}</text>

      <!-- 3D Pseudobulk Tensor Blocks -->
      <rect x="105" y="14" width="70" height="42" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1.2"/>
      <rect x="115" y="8" width="70" height="42" rx="3" fill="#334155" fill-opacity="0.4" stroke="#c084fc" stroke-width="1"/>
      <text x="140" y="32" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">Genes × Donors</text>
      <text x="150" y="44" text-anchor="middle" fill="#c084fc" font-size="6.5">by CellType k</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 3차원 텐서 계수 집계:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Z_{g,d,k} = Σ_{c ∈ C_{d,k}} Y_{g,c} (정수 카운트 합산 보존)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 희소 드롭아웃(Dropout)의 극적 해소:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">환자당 수백 세포 합산으로 벌크 수준의 깊은 시퀀싱 확보</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 라이브러리 크기 인자(Size Factor) 정규화:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">TMM / RLE 기반으로 세포 수 및 시퀀싱 깊이 차이 보정</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">세포형 고유의 신호는 보존하고 기술적 노이즈만 소멸</text>
  </g>

  <!-- Panel 2: Mean-Dispersion Trend Line Estimation -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 평균-산포 공통 추세선 (Trend Estimation)</text>

    <!-- Schematic: Raw MLE scatter vs smooth trend line -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Axes -->
      <line x1="20" y1="58" x2="255" y2="58" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="58" x2="20" y2="6" stroke="#475569" stroke-width="1"/>
      <text x="20" y="4" fill="#64748b" font-size="6">Dispersion (α)</text>
      <text x="235" y="65" fill="#64748b" font-size="6">Mean (μ)</text>

      <!-- Raw MLE noisy scatter -->
      <circle cx="45" cy="18" r="2" fill="#94a3b8" fill-opacity="0.7"/>
      <circle cx="55" cy="40" r="2" fill="#94a3b8" fill-opacity="0.7"/>
      <circle cx="80" cy="22" r="2" fill="#94a3b8" fill-opacity="0.7"/>
      <circle cx="110" cy="48" r="2" fill="#94a3b8" fill-opacity="0.7"/>
      <circle cx="150" cy="35" r="2" fill="#94a3b8" fill-opacity="0.7"/>
      <circle cx="200" cy="52" r="2" fill="#94a3b8" fill-opacity="0.7"/>

      <!-- Smooth Trend Curve alpha_0(mu) -->
      <path d="M 30 15 Q 80 32 150 44 T 245 50" fill="none" stroke="#c084fc" stroke-width="2.2"/>
      <text x="145" y="32" fill="#c084fc" font-size="7.5" font-weight="700">α_0(μ) = A/μ + B (Trend)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자별 산포 MLE의 한계:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">소수 환자(N=3~6) 조건에서 개별 산포 추정값의 높은 분산</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 전체 전사체 정보 공유 (Information Sharing):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">수만 개 유전자의 평균 발현 수준별 공통 산포 추세 모델링</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 국소 회귀(Loess / Spline) 적합:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">저발현의 높은 산포와 고발현의 점근적 생물학적 편차 포착</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#c084fc" font-size="9.2" font-weight="600">공통 사전분포(Empirical Bayes Prior)의 기준 앵커 수립</text>
  </g>

  <!-- Panel 3: Empirical Bayes Shrinkage -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 경험적 베이즈 사후 수축 (EB Shrinkage)</text>

    <!-- Schematic: Shrinkage Arrows toward Trend -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Trend line -->
      <path d="M 25 20 Q 110 38 245 48" fill="none" stroke="#c084fc" stroke-width="1.5" stroke-dasharray="3,3"/>
      <text x="195" y="42" fill="#c084fc" font-size="6">Prior Trend</text>

      <!-- Shrunk point 1 (pulled down) -->
      <circle cx="65" cy="10" r="2.5" fill="#f43f5e"/>
      <path d="M 65 13 L 65 24" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald03)"/>
      <circle cx="65" cy="26" r="3" fill="#34d399"/>

      <!-- Shrunk point 2 (pulled up) -->
      <circle cx="140" cy="52" r="2.5" fill="#38bdf8"/>
      <path d="M 140 49 L 140 42" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald03)"/>
      <circle cx="140" cy="40" r="3" fill="#34d399"/>

      <text x="80" y="16" fill="#f43f5e" font-size="6">Raw MLE (Over-disp)</text>
      <text x="145" y="58" fill="#38bdf8" font-size="6">Raw MLE (Under-disp)</text>
      <text x="100" y="32" fill="#34d399" font-size="7" font-weight="700">Shrunk α~</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 경험적 베이즈 사후 산포 추정:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">α~_g = [ df_prior·α_0(μ_g) + df_data·α^_g ] / (df_prior + df_data)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 산포 과소추정으로 인한 위양성 차단:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">우연히 작은 분산을 보인 유전자가 극단적 p-값을 얻는 오류 방지</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이상치(Outlier) 유전자의 강건성 유지:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">강력한 생물학적 근거가 있는 유전자는 데이터의 자유도로 방어</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">결과: 소표본 임상 전사체 분석의 통계적 검정력 및 안정성 극대화</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_glmm_hierarchical_random_effects_architecture(filepath: str):
    """
    4. GLMM Hierarchical Random Effects Architecture & NEBULA Acceleration
    Cell-level fixed covariates + Donor-level random effects + NEBULA block-diagonal acceleration.
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
    <marker id="arrPurple04" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">2단계 계층적 일반화 선형 혼합 모델(GLMM) &amp; NEBULA 가속 기하학</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">세포 수준 고정 효과(Covariates) + 공여자 수준 무작위 효과(Random Effects) + 블록-대각 정보 행렬 가속</text>

  <!-- Panel 1: Hierarchical Formulation -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 계층적 2단계 수식 아키텍처</text>

    <!-- Schematic: Formula block -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <rect x="5" y="5" width="240" height="26" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="125" y="22" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">log(μ_ic) = log(S_c) + X_c^T β + u_i</text>

      <text x="35" y="45" fill="#34d399" font-size="7">Y_ic ~ NB(μ_ic, ϕ)</text>
      <text x="145" y="45" fill="#fbbf24" font-size="7">u_i ~ N(0, σ²_donor)</text>
      <text x="125" y="58" text-anchor="middle" fill="#94a3b8" font-size="6.5">i: Donor ID, c: Single Cell ID, S_c: Sequencing Depth</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 세포 수준 고정 효과 (Fixed Effects β):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">약물 반응군(Responder), 세포주기 점수, 미토콘드리아 %</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 공여자 무작위 절편 (Random Intercept u_i):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">환자 간 기저 전사체 편차 및 유전적 배경 차이 흡수</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 오프셋(Offset)을 통한 깊이 보정:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">세포별 총 유전자 발현량 log(S_c)를 계수 1로 고정하여 정규화</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">세포 수준 연속 공변량 보정과 환자 통계 엄밀성의 결합</text>
  </g>

  <!-- Panel 2: Cell-within-Donor Nested Hierarchy Tree -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 공여자-세포 내포 계층 트리</text>

    <!-- Schematic: Tree of Donors to Cells -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Root Cohort -->
      <circle cx="140" cy="8" r="6" fill="#38bdf8"/>
      <text x="140" y="7" text-anchor="middle" fill="#0f172a" font-size="6" font-weight="700">C</text>

      <!-- Branches to Donors -->
      <line x1="140" y1="14" x2="60" y2="28" stroke="#475569" stroke-width="1"/>
      <line x1="140" y1="14" x2="140" y2="28" stroke="#475569" stroke-width="1"/>
      <line x1="140" y1="14" x2="220" y2="28" stroke="#475569" stroke-width="1"/>

      <rect x="35" y="25" width="50" height="15" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="60" y="36" text-anchor="middle" fill="#fbbf24" font-size="6.5" font-weight="700">Donor 1 (u_1)</text>

      <rect x="115" y="25" width="50" height="15" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="140" y="36" text-anchor="middle" fill="#fbbf24" font-size="6.5" font-weight="700">Donor 2 (u_2)</text>

      <rect x="195" y="25" width="50" height="15" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="220" y="36" text-anchor="middle" fill="#fbbf24" font-size="6.5" font-weight="700">Donor M (u_M)</text>

      <!-- Branches to Single Cells -->
      <line x1="60" y1="40" x2="45" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <line x1="60" y1="40" x2="75" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <circle cx="45" cy="54" r="3" fill="#34d399"/>
      <circle cx="75" cy="54" r="3" fill="#34d399"/>

      <line x1="140" y1="40" x2="125" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <line x1="140" y1="40" x2="155" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <circle cx="125" cy="54" r="3" fill="#34d399"/>
      <circle cx="155" cy="54" r="3" fill="#34d399"/>

      <line x1="220" y1="40" x2="205" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <line x1="220" y1="40" x2="235" y2="52" stroke="#64748b" stroke-width="0.8"/>
      <circle cx="205" cy="54" r="3" fill="#34d399"/>
      <circle cx="235" cy="54" r="3" fill="#34d399"/>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 분산 분해 원리 (Variance Partitioning):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">σ²_total = σ²_donor (환자간) + σ²_cell (세포간) + 샷노이즈</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 동일 환자 세포의 잠재 상관성 흡수:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">무작위 효과 u_i를 통해 환자 내 세포 간 유사도 완벽 분리</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 의사벌크가 평균화로 잃는 정보 보존:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">개별 세포의 상태(미세 이질성, 아형 비율)를 직접 활용</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">환자 간 이질성이 높은 면역 세포 및 종양 클론 분석의 최적해</text>
  </g>

  <!-- Panel 3: NEBULA Block-Diagonal Speedup -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad04)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">3. NEBULA 블록 대각 정보 행렬 가속</text>

    <!-- Schematic: Dense vs Block-Diagonal Hessian -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Dense Hessian (standard glmm) -->
      <rect x="15" y="8" width="46" height="46" rx="2" fill="#1e293b" stroke="#f43f5e" stroke-width="1.2"/>
      <path d="M 18 12 L 58 50 M 18 50 L 58 12" stroke="#f43f5e" stroke-width="0.8"/>
      <text x="38" y="62" text-anchor="middle" fill="#f43f5e" font-size="6">Dense O(N³)</text>

      <!-- Versus Arrow -->
      <line x1="75" y1="31" x2="105" y2="31" stroke="#c084fc" stroke-width="1.5" marker-end="url(#arrPurple04)"/>
      <text x="90" y="24" text-anchor="middle" fill="#fbbf24" font-size="6.5">NEBULA</text>

      <!-- Block-Diagonal Matrix -->
      <rect x="120" y="8" width="55" height="46" rx="2" fill="#1e293b" stroke="#34d399" stroke-width="1.2"/>
      <rect x="123" y="11" width="13" height="11" fill="#34d399"/>
      <rect x="140" y="25" width="13" height="11" fill="#34d399"/>
      <rect x="158" y="40" width="13" height="11" fill="#34d399"/>
      <text x="147" y="62" text-anchor="middle" fill="#34d399" font-size="6">Block-Diag O(N)</text>

      <text x="215" y="28" fill="#34d399" font-size="8" font-weight="700">100x ~ 1000x</text>
      <text x="215" y="40" fill="#cbd5e1" font-size="6.8">Speedup!</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 전통적 GLMM 역행렬 계산 병목:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">N_cells × N_cells 조밀 헤시안 역행렬 계산으로 수 시간 소요</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• NEBULA 블록 대각 구조 분해:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">환자별 그룹화를 통해 피셔 정보 행렬 I(β,u)를 대각 블록화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 셔먼-모리슨-우드베리(SMW) 공식:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">선형 시간 O(M·m²) 닫힌 형태 역행렬로 유전자당 밀리초 달성</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#c084fc" font-size="9.2" font-weight="600">수십만 세포 코호트에서 전유전체(20,000개) 실시간 적합 가능</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_milo_knn_neighborhood_differential_abundance(filepath: str):
    """
    5. Milo kNN Neighborhood Differential Abundance (DA) Framework
    kNN graph -> Representative index sampling -> Overlapping neighborhoods -> NB GLM & Spatial FDR.
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
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">Milo kNN 이웃 그래프 차등 세포 풍부도(Differential Abundance) 프레임워크</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">고차원 kNN 그래프 ➔ 대표 인덱스 세포 추출 ➔ 중첩 이웃 계수 행렬 ➔ 음이항 GLM &amp; 공간 FDR 보정</text>

  <!-- Panel 1: Graph & Index Cell Sampling -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. kNN 그래프 &amp; 인덱스 샘플링</text>

    <!-- Schematic: Graph with overlapping neighborhoods -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Edges -->
      <line x1="40" y1="20" x2="80" y2="15" stroke="#334155" stroke-width="1"/>
      <line x1="40" y1="20" x2="60" y2="45" stroke="#334155" stroke-width="1"/>
      <line x1="80" y1="15" x2="120" y2="35" stroke="#334155" stroke-width="1"/>
      <line x1="60" y1="45" x2="120" y2="35" stroke="#334155" stroke-width="1"/>
      <line x1="120" y1="35" x2="170" y2="25" stroke="#334155" stroke-width="1"/>
      <line x1="170" y1="25" x2="210" y2="48" stroke="#334155" stroke-width="1"/>

      <!-- Neighborhood 1 (dashed circle) -->
      <circle cx="60" cy="25" r="22" fill="#38bdf8" fill-opacity="0.15" stroke="#38bdf8" stroke-width="1" stroke-dasharray="2,2"/>
      <!-- Neighborhood 2 (dashed circle overlapping) -->
      <circle cx="120" cy="35" r="24" fill="#fbbf24" fill-opacity="0.15" stroke="#fbbf24" stroke-width="1" stroke-dasharray="2,2"/>

      <!-- Nodes -->
      <circle cx="40" cy="20" r="3.5" fill="#94a3b8"/>
      <circle cx="80" cy="15" r="3.5" fill="#94a3b8"/>
      <circle cx="60" cy="45" r="3.5" fill="#94a3b8"/>
      <circle cx="170" cy="25" r="3.5" fill="#94a3b8"/>
      <circle cx="210" cy="48" r="3.5" fill="#94a3b8"/>

      <!-- Index Cells (Star/Gold) -->
      <circle cx="60" cy="25" r="5" fill="#38bdf8"/>
      <circle cx="120" cy="35" r="5" fill="#fbbf24"/>
      <text x="60" y="24" text-anchor="middle" fill="#0f172a" font-size="5" font-weight="700">v1</text>
      <text x="120" y="34" text-anchor="middle" fill="#0f172a" font-size="5" font-weight="700">v2</text>
      <text x="210" y="20" fill="#fbbf24" font-size="6.5">Overlap!</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이산적 클러스터링 경계 오류 탈피:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Leiden/Louvain의 주관적 해상도 및 경계 인위성 극복</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 그래프 기반 대표 인덱스 샘플링:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">그래프 위상 전체를 균일하게 덮는 대표 세포 추출 (~10%)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 연속적 전이 표현형(Cell State) 보존:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">각 인덱스 세포의 kNN 반경 내 중첩 이웃(Neighborhood) 정의</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">세포 상태의 연속적인 궤적과 아형 변화를 보존하는 기하학</text>
  </g>

  <!-- Panel 2: Neighborhood Count Matrix Assembly -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 이웃별 환자 계수 행렬 (Count Matrix)</text>

    <!-- Schematic: Matrix M x D -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Matrix Box -->
      <rect x="10" y="6" width="160" height="48" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="90" y="18" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">Donors (D_1 ... D_P)</text>
      <!-- Grid rows -->
      <line x1="10" y1="22" x2="170" y2="22" stroke="#475569" stroke-width="0.8"/>
      <line x1="10" y1="36" x2="170" y2="36" stroke="#475569" stroke-width="0.8"/>
      <text x="15" y="31" fill="#38bdf8" font-size="6">N_1 (Index 1)</text>
      <text x="15" y="46" fill="#38bdf8" font-size="6">N_m (Index m)</text>
      <!-- Values -->
      <text x="80" y="31" fill="#cbd5e1" font-size="6.5">14</text>
      <text x="110" y="31" fill="#cbd5e1" font-size="6.5">2</text>
      <text x="145" y="31" fill="#f43f5e" font-size="6.5">38</text>
      <text x="80" y="46" fill="#cbd5e1" font-size="6.5">5</text>
      <text x="110" y="46" fill="#34d399" font-size="6.5">42</text>
      <text x="145" y="46" fill="#cbd5e1" font-size="6.5">3</text>

      <!-- Modeling Arrow -->
      <line x1="180" y1="30" x2="205" y2="30" stroke="#fbbf24" stroke-width="1.5" marker-end="url(#arrAmber05)"/>
      <rect x="210" y="10" width="55" height="40" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
      <text x="237" y="22" text-anchor="middle" fill="#34d399" font-size="6.5" font-weight="700">NB GLM</text>
      <text x="237" y="34" text-anchor="middle" fill="#cbd5e1" font-size="5.5">Log2FC</text>
      <text x="237" y="44" text-anchor="middle" fill="#fbbf24" font-size="5.5">p-value</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이웃별 환자 세포 수 집계:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">C_{m,d} = Σ_{c ∈ d} I(c ∈ N_m) (M개 이웃 × P명 환자 행렬)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 환자별 총 세포 수 오프셋 정규화:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">환자간 총 세포 포획 효율 편차를 라이브러리 크기로 보정</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 음이항 GLM을 통한 조건별 효과 검정:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">반응군 vs 비반응군 간 세포 풍부도 변화 log2FC 및 유의성 도출</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.2" font-weight="600">세포 조성 차이를 신뢰할 수 있는 카운트 통계 검정으로 전환</text>
  </g>

  <!-- Panel 3: Spatial FDR Smoothing & Beeswarm -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 공간 FDR 평활화 &amp; Beeswarm</text>

    <!-- Schematic: Beeswarm Plot of Neighborhoods -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Central zero FC line -->
      <line x1="20" y1="31" x2="255" y2="31" stroke="#475569" stroke-width="1"/>
      <text x="25" y="25" fill="#64748b" font-size="6">Log2FC &gt; 0 (Resp Enriched)</text>
      <text x="25" y="45" fill="#64748b" font-size="6">Log2FC &lt; 0 (Non-Resp Enriched)</text>

      <!-- Neighborhood dots across phenotypes -->
      <!-- Naive T (near zero) -->
      <circle cx="65" cy="30" r="3" fill="#94a3b8"/>
      <circle cx="75" cy="33" r="2.5" fill="#94a3b8"/>
      <text x="70" y="62" text-anchor="middle" fill="#94a3b8" font-size="6">Naive T</text>

      <!-- Effector T (Enriched in Responders: Blue/Cyan) -->
      <circle cx="120" cy="18" r="4" fill="#38bdf8"/>
      <circle cx="130" cy="12" r="4.5" fill="#38bdf8"/>
      <circle cx="140" cy="16" r="3.5" fill="#38bdf8"/>
      <text x="130" y="62" text-anchor="middle" fill="#38bdf8" font-size="6" font-weight="700">Effector T</text>

      <!-- Exhausted T (Enriched in Non-Responders: Red/Rose) -->
      <circle cx="190" cy="48" r="4.5" fill="#f43f5e"/>
      <circle cx="205" cy="52" r="5" fill="#f43f5e"/>
      <circle cx="220" cy="46" r="4" fill="#f43f5e"/>
      <text x="205" y="62" text-anchor="middle" fill="#f43f5e" font-size="6" font-weight="700">Exhausted T</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 중첩 이웃 가설 검정의 다중성 보정:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">공간적 상관성을 반영하는 가중 Benjamini-Hochberg 보정</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 그래프 거리 기반 공간 가중치 (Spatial Weights):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">이웃 간 연결 강도(W_ij)를 통해 국소적 p-값 평활화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 소진(Exhaustion) 대 작동(Effector) T세포 분리:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">약물 저항 환자에서 특이적으로 팽창하는 면역 아형 확정</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">임상 시험 바이오마커 발굴을 위한 최상위 정밀 차등 풍부도 도구</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_compositional_simplex_and_coda_regression(filepath: str):
    """
    6. Compositional Simplex Constraint & scCODA Bayesian Regression
    Simplex constraint (sum p_c = 1) spurious negative correlation -> CLR projection -> scCODA Dirichlet-Multinomial.
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
    <marker id="arrEmerald06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title & Subtitle -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">조성 데이터(CoDa) 심플렉스 제약 &amp; scCODA 베이지안 회귀 분석</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">합 = 1 심플렉스 제약의 허위 음의 상관관계 ➔ 중심 로그비(CLR) 유클리드 투영 ➔ 희소 디리클레-다항 회귀</text>

  <!-- Panel 1: The Simplex Constraint & Spurious Depletion -->
  <g transform="translate(20, 65)">
    <rect width="300" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="300" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="150" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. 심플렉스 제약과 허위 감소 참사</text>

    <!-- Schematic: Proportional Stacked Bars Collapse -->
    <rect x="15" y="42" width="270" height="74" rx="6" fill="#0f172a" stroke="#e11d48" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Baseline Patient A -->
      <rect x="20" y="8" width="35" height="24" fill="#38bdf8"/>
      <rect x="20" y="32" width="35" height="24" fill="#fbbf24"/>
      <text x="37" y="22" text-anchor="middle" fill="#0f172a" font-size="6" font-weight="700">50% T</text>
      <text x="37" y="46" text-anchor="middle" fill="#0f172a" font-size="6" font-weight="700">50% M</text>
      <text x="37" y="65" text-anchor="middle" fill="#94a3b8" font-size="6">Patient A (100 T, 100 M)</text>

      <!-- Arrow -->
      <path d="M 68 32 L 95 32" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrCyan06)"/>
      <text x="82" y="24" text-anchor="middle" fill="#f43f5e" font-size="6.5">+900 T!</text>

      <!-- Patient B (T expanded 10x, M unchanged at 100) -->
      <rect x="110" y="8" width="35" height="42" fill="#38bdf8"/>
      <rect x="110" y="50" width="35" height="6" fill="#fbbf24"/>
      <text x="127" y="28" text-anchor="middle" fill="#0f172a" font-size="6.5" font-weight="700">91% T</text>
      <text x="127" y="55" text-anchor="middle" fill="#0f172a" font-size="5" font-weight="700">9%</text>
      <text x="127" y="65" text-anchor="middle" fill="#f43f5e" font-size="6">M depleted?! (False!)</text>

      <!-- Warning text -->
      <text x="185" y="25" fill="#f43f5e" font-size="7" font-weight="700">Σ p_k = 1</text>
      <text x="185" y="37" fill="#cbd5e1" font-size="6">스퓨리어스 음의 상관</text>
      <text x="185" y="47" fill="#fbbf24" font-size="6">t-test 시 위양성</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 심플렉스 닫힘 제약 (Closure Constraint):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">S^D = { p ∈ R^D : p_k &gt; 0, Σ_{k=1}^D p_k = 1 } (자유도 D-1)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 허위 음의 상관관계 (Spurious Correlation):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">한 세포형의 진정한 팽창이 다른 무변화 세포형의 비율 급감 초래</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 원시 비율(Percentage) 검정의 위험성:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">t-test나 ANOVA를 비율에 직접 적용 시 수많은 가짜 감소 세포 식별</text>

    <!-- Badge -->
    <rect x="14" y="242" width="272" height="22" rx="4" fill="#e11d48" fill-opacity="0.22"/>
    <text x="150" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.2" font-weight="600">세포 비율의 독립성 가정을 전제하는 표준 통계학의 완전한 붕괴</text>
  </g>

  <!-- Panel 2: Aitchison Geometry & CLR Projection -->
  <g transform="translate(340, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">2. 아치슨 기하학과 중심 로그비(CLR)</text>

    <!-- Schematic: Simplex to Euclidean Plane Projection -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Ternary Simplex Triangle -->
      <polygon points="40,55 75,12 110,55" fill="none" stroke="#fbbf24" stroke-width="1.5"/>
      <circle cx="75" cy="38" r="3" fill="#f43f5e"/>
      <text x="75" y="8" text-anchor="middle" fill="#fbbf24" font-size="6">S³ Simplex</text>

      <!-- Projection Arrow -->
      <line x1="125" y1="34" x2="155" y2="34" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan06)"/>
      <text x="140" y="26" text-anchor="middle" fill="#38bdf8" font-size="6.5">CLR / ilr</text>

      <!-- Euclidean 2D Plane -->
      <rect x="170" y="10" width="85" height="46" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
      <line x1="175" y1="33" x2="250" y2="33" stroke="#475569" stroke-width="0.8"/>
      <line x1="212" y1="12" x2="212" y2="53" stroke="#475569" stroke-width="0.8"/>
      <circle cx="230" cy="22" r="3" fill="#34d399"/>
      <text x="212" y="62" text-anchor="middle" fill="#34d399" font-size="6">R^(D-1) Euclidean Space</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 중심 로그비 공식 (Centered Log-Ratio):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">CLR(p) = [ log(p_1 / g(p)), ..., log(p_D / g(p)) ], g(p) = (Π p_k)^(1/D)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 아치슨 기하학(Aitchison Geometry)의 해방:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">심플렉스 곡면을 유클리드 무제약 공간 R^(D-1)으로 선형 사상</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 참조 세포형(Reference Cell Type) 선택:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">안정적 기저 세포형(예: 섬유아세포)을 분모로 설정하여 왜곡 해소</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.2" font-weight="600">조성 데이터를 위한 유일한 수학적으로 엄밀한 변환 프레임워크</text>
  </g>

  <!-- Panel 3: scCODA Sparse Bayesian Regression -->
  <g transform="translate(670, 65)">
    <rect width="310" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="310" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="155" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. scCODA 희소 베이지안 회귀 분석</text>

    <!-- Schematic: Forest plot of scCODA effects -->
    <rect x="15" y="42" width="280" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    <g transform="translate(25, 46)">
      <!-- Zero effect axis -->
      <line x1="140" y1="8" x2="140" y2="58" stroke="#475569" stroke-width="1"/>
      <text x="140" y="65" text-anchor="middle" fill="#64748b" font-size="6">Effect β = 0</text>

      <!-- Cell Type 1: CD8 Exhausted (True Expansion, Credible) -->
      <text x="20" y="18" fill="#f43f5e" font-size="6.5" font-weight="700">CD8 Exhausted</text>
      <line x1="170" y1="16" x2="230" y2="16" stroke="#f43f5e" stroke-width="1.5"/>
      <circle cx="200" cy="16" r="3.5" fill="#f43f5e"/>
      <text x="235" y="19" fill="#f43f5e" font-size="6" font-weight="700">Credible (+)</text>

      <!-- Cell Type 2: CD8 Effector (Shrunk to 0 by spike prior) -->
      <text x="20" y="32" fill="#94a3b8" font-size="6.5">CD8 Effector</text>
      <circle cx="140" cy="30" r="3" fill="#94a3b8"/>
      <text x="150" y="33" fill="#94a3b8" font-size="6">β = 0 (Shrunk)</text>

      <!-- Cell Type 3: Macrophages (Unchanged, Shrunk to 0) -->
      <text x="20" y="46" fill="#34d399" font-size="6.5">Macrophage</text>
      <circle cx="140" cy="44" r="3" fill="#34d399"/>
      <text x="150" y="47" fill="#34d399" font-size="6">β = 0 (Ref/Invariant)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 디리클레-다항 분포 (Dirichlet-Multinomial):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">Y_d ~ Multinomial(N_d, π_d), π_d ~ Dirichlet(α_d)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 스파이크 앤 슬랩(Spike-and-Slab) 희소 사전분포:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">유의하지 않은 미세 변화를 사후 확률적으로 정확히 0으로 수축</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 최적 참조 세포형 자동 탐색:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">환자 간 변동성이 가장 낮은 불변 세포형을 기준선으로 확정</text>

    <!-- Badge -->
    <rect x="14" y="242" width="282" height="22" rx="4" fill="#059669" fill-opacity="0.22"/>
    <text x="155" y="257" text-anchor="middle" fill="#34d399" font-size="9.2" font-weight="600">소표본(N=5~10) 코호트에서 조성 변화를 규명하는 궁극의 베이지안 해법</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_poisson_vs_negative_binomial_overdispersion.svg", make_01_poisson_vs_negative_binomial_overdispersion),
        ("02_pseudoreplication_type1_error_catastrophe.svg", make_02_pseudoreplication_type1_error_catastrophe),
        ("03_pseudobulk_aggregation_and_dispersion_shrinkage.svg", make_03_pseudobulk_aggregation_and_dispersion_shrinkage),
        ("04_glmm_hierarchical_random_effects_architecture.svg", make_04_glmm_hierarchical_random_effects_architecture),
        ("05_milo_knn_neighborhood_differential_abundance.svg", make_05_milo_knn_neighborhood_differential_abundance),
        ("06_compositional_simplex_and_coda_regression.svg", make_06_compositional_simplex_and_coda_regression),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 07 Count Statistics publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

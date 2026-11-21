#!/usr/bin/env python3
"""
tracks/06_manifold_trajectory/scripts/generate_trajectory_diagrams.py
=============================================================================
Generates publication-grade, high-resolution SVG diagrams for Track 06:
Manifold Learning, Pseudotime Trajectories & RNA Velocity in Single-Cell Biology.

Covers:
1. Diffusion Maps & Pseudotime Geometry (Graph Laplacian, Operator T, Eigen-coords, DPT)
2. RNA Velocity Phase Portrait Kinetics (ODE, Splicing Ratio gamma, Vector Projections)
3. Dynamical scVelo Latent Time Optimization (EM Kinetics, Parameter Estimation, Latent tau)
4. Waddington-OT Sinkhorn Transport Coupling (Push-forward, Kantorovich, Entropic Sinkhorn)
5. CellRank 2 Multi-Kernel & G-PCCA Macrostates (Kernel Fusion, Transition P, Schur Macrostates)
6. Absorbing Markov Fate Probabilities & Lineage Drivers (Matrix B = N R, Driver Ranking)
"""

import os
import xml.etree.ElementTree as ET


def make_01_diffusion_maps_and_pseudotime_geometry(filepath: str):
    """
    1. Diffusion Maps & Pseudotime Geometry
    Graph Laplacian -> Markov Transition Operator T -> Diffusion Eigen-coordinates -> DPT Geodesic Metric.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="diffGrad01" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="50%" stop-color="#c084fc"/>
      <stop offset="100%" stop-color="#34d399"/>
    </linearGradient>
    <marker id="arrCyan01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
    <marker id="arrPurple01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
    <marker id="arrEmerald01" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">확률론적 확산 맵(Diffusion Maps) &amp; 확산 의시간(DPT) 기하학</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">그래프 라플라시안 ➔ 마르코프 전이 연산자 T ➔ 고유공간 스펙트럼 분해 ➔ 분기 매니폴드 측지선 거리</text>

  <!-- Panel 1: Affinity Matrix & Markov Transition Operator -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 친화도 행렬 W &amp; 전이 연산자 T</text>

    <!-- Schematic: kNN graph to Transition Matrix -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <!-- Graph Network Nodes -->
    <g transform="translate(25, 48)">
      <line x1="20" y1="20" x2="55" y2="12" stroke="#334155" stroke-width="1.2"/>
      <line x1="20" y1="20" x2="35" y2="48" stroke="#334155" stroke-width="1.2"/>
      <line x1="55" y1="12" x2="75" y2="38" stroke="#334155" stroke-width="1.2"/>
      <line x1="35" y1="48" x2="75" y2="38" stroke="#334155" stroke-width="1.2"/>
      <circle cx="20" cy="20" r="6" fill="#38bdf8"/>
      <circle cx="55" cy="12" r="5" fill="#818cf8"/>
      <circle cx="35" cy="48" r="5" fill="#38bdf8"/>
      <circle cx="75" cy="38" r="6" fill="#0284c7"/>
      <text x="45" y="60" text-anchor="middle" fill="#94a3b8" font-size="7.5">kNN Graph (W_ij)</text>
    </g>

    <line x1="115" y1="78" x2="135" y2="78" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan01)"/>

    <!-- Normalized Transition Matrix -->
    <g transform="translate(145, 50)">
      <rect width="120" height="56" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="60" y="18" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">T = D^(-1) W</text>
      <text x="60" y="32" text-anchor="middle" fill="#cbd5e1" font-size="7.5">D_ii = Σ_j W_ij</text>
      <text x="60" y="46" text-anchor="middle" fill="#fbbf24" font-size="7.5">Σ_j T_ij = 1 (확률 보존)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 국소 적응형 가우시안 친화도 행렬:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">W_ij = exp( -||x_i - x_j||^2 / (2 σ_i σ_j) ) (밀도 적응)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 마르코프 무작위 보행 (Random Walk):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">T_ij: 세포 i에서 세포 j로 1단계 전이 확률</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 밀도 정규화 (α-Normalization):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">W^(α) = D^(-α) W D^(-α)로 샘플링 밀도 바이어스 제거</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">수학적 특징: 연속 다양체 상의 열 확산 방정식 이산화</text>
  </g>

  <!-- Panel 2: Diffusion Eigen-Coordinates -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad01)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 확산 고유공간(Diffusion Space)</text>

    <!-- Schematic: Eigenvalue Decay & 2D Projection -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>

    <!-- Eigenvalue spectrum curve -->
    <g transform="translate(25, 48)">
      <line x1="5" y1="52" x2="95" y2="52" stroke="#475569" stroke-width="1"/>
      <line x1="5" y1="52" x2="5" y2="8" stroke="#475569" stroke-width="1"/>
      <path d="M 5 10 Q 18 25 35 44 T 90 51" fill="none" stroke="#fbbf24" stroke-width="1.6"/>
      <circle cx="5" cy="10" r="2.5" fill="#38bdf8"/>
      <circle cx="20" cy="28" r="2.5" fill="#c084fc"/>
      <circle cx="35" cy="44" r="2.5" fill="#34d399"/>
      <text x="50" y="16" fill="#fbbf24" font-size="7.5">1=λ_0 &gt; λ_1 ≥ λ_2</text>
      <text x="50" y="61" fill="#94a3b8" font-size="7">Spectral Gap</text>
    </g>

    <!-- Arrow -->
    <line x1="128" y1="78" x2="148" y2="78" stroke="#c084fc" stroke-width="1.5" marker-end="url(#arrPurple01)"/>

    <!-- 2D Diffusion Space -->
    <g transform="translate(155, 47)">
      <rect width="140" height="62" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <!-- Trajectory arc -->
      <path d="M 15 48 Q 55 12 125 24" fill="none" stroke="url(#diffGrad01)" stroke-width="2.5"/>
      <circle cx="15" cy="48" r="4" fill="#38bdf8"/>
      <circle cx="68" cy="20" r="3.5" fill="#c084fc"/>
      <circle cx="125" cy="24" r="4" fill="#34d399"/>
      <text x="20" y="58" fill="#38bdf8" font-size="7" font-weight="700">Root (ψ_1)</text>
      <text x="110" y="40" fill="#34d399" font-size="7" font-weight="700">Fate Branch</text>
      <text x="70" y="12" text-anchor="middle" fill="#94a3b8" font-size="6.5">ψ_1 vs ψ_2 Manifold</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 스펙트럼 분해 (Spectral Decomposition):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">T ψ_k = λ_k ψ_k (k=1,...,d 차원 매핑)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 확산 거리 (Diffusion Distance):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">D_t^2(x, y) = Σ_k λ_k^(2t) (ψ_k(x) - ψ_k(y))^2</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 노이즈 필터링 메커니즘:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">시간 t 거듭제곱으로 단기 노이즈(작은 λ_k) 지수적 감쇠</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">스펙트럼 수렴: 리만 다양체의 본질적 저차원 측지선 보존</text>
  </g>

  <!-- Panel 3: DPT Geodesic Metric & Branching -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 확산 의시간(DPT) &amp; 분기 매니폴드</text>

    <!-- Schematic: Branching Geodesics -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Branching lines -->
      <path d="M 15 32 Q 55 32 85 32" stroke="#38bdf8" stroke-width="2.2" fill="none"/>
      <path d="M 85 32 Q 130 18 165 14" stroke="#f43f5e" stroke-width="2.2" fill="none"/>
      <path d="M 85 32 Q 130 45 165 50" stroke="#34d399" stroke-width="2.2" fill="none"/>

      <!-- Root node -->
      <circle cx="15" cy="32" r="5" fill="#38bdf8"/>
      <text x="15" y="46" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">Root x_0</text>

      <!-- Branch point -->
      <circle cx="85" cy="32" r="4" fill="#fbbf24"/>
      <text x="85" y="24" text-anchor="middle" fill="#fbbf24" font-size="6.5">Bifurcation</text>

      <!-- Terminals -->
      <circle cx="165" cy="14" r="5" fill="#f43f5e"/>
      <text x="200" y="16" fill="#f43f5e" font-size="7" font-weight="700">t-NEPC</text>

      <circle cx="165" cy="50" r="5" fill="#34d399"/>
      <text x="200" y="52" fill="#34d399" font-size="7" font-weight="700">DNPC</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 확산 의시간 공식 (Diffusion Pseudotime):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">dpt(x0,x) = ||M(e_x0 - e_x)||, M = Σ_t (T^t - φ0φ0ᵀ)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 누적 무작위 보행 (Accumulated Walk):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">모든 시간 t∈[1, ∞)에 걸친 연속 전이 확률 적분</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• mCRPC 계통 가소성 경로 정량화:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">AR+ 선암종 기원 ➔ t-NEPC 및 DNPC 2개 분기 운명 추적</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">임상적 가치: 약제 유발 신경내분비 역분화의 연속적 진행도 산출</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_rna_velocity_phase_portrait_kinetics(filepath: str):
    """
    2. RNA Velocity Phase Portrait Kinetics
    Splicing ODE (u -> s) -> Phase portrait slope gamma -> Velocity vector projection.
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
    <marker id="arrEmerald02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">RNA 속도(RNA Velocity) 기본 동역학 &amp; 위상 평면(Phase Portrait)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">스플라이싱 미분방정식(ODE) ➔ 위상 평면 정상상태 기울기 γ ➔ 2D 매니폴드 속도 벡터 투영</text>

  <!-- Panel 1: Splicing Kinetics ODE -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 전사-스플라이싱 반응 속도론</text>

    <!-- Schematic: Transcription ODE flow -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <g transform="translate(25, 52)">
      <!-- Box 1: DNA/Alpha -->
      <rect width="45" height="26" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="22" y="16" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">DNA (α)</text>

      <line x1="47" y1="13" x2="68" y2="13" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan02)"/>

      <!-- Box 2: Unspliced u -->
      <g transform="translate(72, 0)">
        <rect width="48" height="26" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
        <text x="24" y="16" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">u (미성숙)</text>
      </g>

      <line x1="122" y1="13" x2="143" y2="13" stroke="#fbbf24" stroke-width="1.5" marker-end="url(#arrAmber02)"/>

      <!-- Box 3: Spliced s -->
      <g transform="translate(147, 0)">
        <rect width="48" height="26" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
        <text x="24" y="16" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">s (성숙)</text>
      </g>

      <line x1="197" y1="13" x2="218" y2="13" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald02)"/>

      <!-- Degradation -->
      <text x="226" y="17" fill="#f43f5e" font-size="10" font-weight="700">∅ (γ)</text>

      <text x="120" y="44" text-anchor="middle" fill="#cbd5e1" font-size="8">du/dt = α - βu,   ds/dt = βu - γs</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 상미분방정식 (Kinetics ODE System):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">전사 속도 α(t), 인트론 스플라이싱 β, mRNA 분해율 γ</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• RNA 속도 정의 (Velocity Definition):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">v = ds/dt = u - γ s (단위 스플라이싱율 β ≡ 1 정규화)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 미래 상태 예측 원리 (Time Horizon):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">현재 미성숙/성숙 불균형으로 수 시간 뒤 발현 예측</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">핵심 가정: 스플라이싱 반응은 분해보다 수 배 빠른 선행 과정</text>
  </g>

  <!-- Panel 2: Phase Portrait & Extreme Slope -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 위상 평면(Phase Portrait) &amp; 정상상태</text>

    <!-- Schematic: u vs s Phase plane -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>

    <g transform="translate(30, 48)">
      <!-- Axes -->
      <line x1="15" y1="56" x2="115" y2="56" stroke="#475569" stroke-width="1"/>
      <line x1="15" y1="56" x2="15" y2="6" stroke="#475569" stroke-width="1"/>
      <text x="115" y="60" fill="#94a3b8" font-size="6.5">s (Spliced)</text>
      <text x="10" y="8" fill="#94a3b8" font-size="6.5">u</text>

      <!-- Steady line u = gamma * s -->
      <line x1="15" y1="56" x2="95" y2="12" stroke="#64748b" stroke-width="1.4" stroke-dasharray="3,2"/>
      <text x="75" y="16" fill="#fbbf24" font-size="7">기울기 γ</text>

      <!-- Induction curve (upper) & Repression curve (lower) -->
      <path d="M 15 56 Q 30 14 65 18" stroke="#38bdf8" stroke-width="1.6" fill="none"/>
      <path d="M 65 18 Q 95 24 95 48 Q 65 56 15 56" stroke="#f43f5e" stroke-width="1.6" fill="none"/>

      <text x="35" y="26" fill="#38bdf8" font-size="6.5" font-weight="700">유도(v &gt; 0)</text>
      <text x="80" y="42" fill="#f43f5e" font-size="6.5" font-weight="700">억제(v &lt; 0)</text>
    </g>

    <g transform="translate(180, 50)">
      <rect width="115" height="56" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="57" y="16" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">정상상태(Steady-State)</text>
      <text x="57" y="30" text-anchor="middle" fill="#cbd5e1" font-size="7">u = γ s ➔ v = 0</text>
      <text x="57" y="46" text-anchor="middle" fill="#38bdf8" font-size="7">γ = (s^T u) / (s^T s)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 극단 분위수 선형 회귀 (Quantile Regression):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">상위/하위 극단 발현 세포들을 통과하는 직교 기울기 γ 추정</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 속도 부호의 물리적 의미:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">v &gt; 0 (선 상단): 유도 중 ➔ mRNA 합성 가속</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 정상상태 가정의 기술적 한계:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">세포 집단 내에 완전한 평형 상태 세포가 반드시 존재해야 함</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.5" font-weight="600">진단 지표: R^2 적합도 &amp; 직교 잔차 분석으로 유효 유전자 선별</text>
  </g>

  <!-- Panel 3: Manifold Velocity Vectors -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 저차원 매니폴드 속도 벡터 장</text>

    <!-- Schematic: 2D Velocity Arrows -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Cells with arrows -->
      <circle cx="20" cy="30" r="4" fill="#38bdf8"/>
      <line x1="20" y1="30" x2="42" y2="24" stroke="#38bdf8" stroke-width="2" marker-end="url(#arrCyan02)"/>

      <circle cx="60" cy="22" r="4" fill="#818cf8"/>
      <line x1="60" y1="22" x2="88" y2="18" stroke="#818cf8" stroke-width="2" marker-end="url(#arrCyan02)"/>

      <circle cx="110" cy="20" r="4" fill="#fbbf24"/>
      <line x1="110" y1="20" x2="135" y2="35" stroke="#fbbf24" stroke-width="2" marker-end="url(#arrAmber02)"/>

      <circle cx="150" cy="42" r="4" fill="#34d399"/>
      <line x1="150" y1="42" x2="182" y2="48" stroke="#34d399" stroke-width="2" marker-end="url(#arrEmerald02)"/>

      <text x="20" y="52" fill="#38bdf8" font-size="7" font-weight="700">AR+ Adeno</text>
      <text x="100" y="10" fill="#fbbf24" font-size="7" font-weight="700">가소성 분기</text>
      <text x="160" y="60" fill="#34d399" font-size="7" font-weight="700">t-NEPC</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이웃 전이 확률 (Transition Probability):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">P_ij ∝ exp( corr(v_i, s_j - s_i) / σ )</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 2D 그리드 투영 (Vector Field):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">Δx_i = Σ_j (P_ij - 1/n) (x_j - x_i) 로 2D 화살표 산출</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 방향성 인과 추론 (Causal Trajectory):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">루미날 상피에서 신경내분비 계통으로의 비가역 전이 입증</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">결과: 정적 스냅샷 데이터에서 시간 역학의 화살(Arrow of Time) 복원</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_dynamical_scvelo_latent_time_optimization(filepath: str):
    """
    3. Dynamical scVelo Latent Time Optimization
    Steady-state failure vs Dynamical EM trajectory -> Cell-specific latent time tau.
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
    <marker id="arrRose03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrPurple03" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">scVelo 동적 모델(Dynamical Model) &amp; 전역 잠재 시간(Latent Time)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">정상상태 가정 탈피 ➔ 비선형 전사항 완전 동역학 EM 최적화 ➔ 세포 고유 전역 잠재 시간 τ</text>

  <!-- Panel 1: Steady-State Failure vs Full Dynamical Loop -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad03)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. 선형 한계 vs 완전 동역학 루프</text>

    <!-- Schematic: Closed ellipse loop vs straight line -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>

    <g transform="translate(30, 48)">
      <!-- Axes -->
      <line x1="10" y1="56" x2="95" y2="56" stroke="#475569" stroke-width="1"/>
      <line x1="10" y1="56" x2="10" y2="6" stroke="#475569" stroke-width="1"/>

      <!-- False linear fit -->
      <line x1="10" y1="56" x2="85" y2="15" stroke="#64748b" stroke-width="1.4" stroke-dasharray="2,2"/>
      <text x="75" y="12" fill="#64748b" font-size="6">선형 오적합</text>

      <!-- True 4-phase loop -->
      <path d="M 10 56 C 25 10, 60 12, 85 28 C 90 45, 60 55, 10 56 Z" fill="none" stroke="#38bdf8" stroke-width="2"/>
      <text x="42" y="24" fill="#38bdf8" font-size="7" font-weight="700">동역학 곡선</text>
    </g>

    <g transform="translate(150, 48)">
      <rect width="120" height="60" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="60" y="16" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">4대 동역학 상태</text>
      <text x="60" y="29" text-anchor="middle" fill="#38bdf8" font-size="6.8">1. 유도 (Induction, α&gt;0)</text>
      <text x="60" y="41" text-anchor="middle" fill="#fbbf24" font-size="6.8">2. 정상 (Steady, u=γs)</text>
      <text x="60" y="53" text-anchor="middle" fill="#c084fc" font-size="6.8">3. 억제 (Repression, α=0)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 정상상태 가정의 치명적 결함:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">과도기(transient) 유전자는 닫힌 루프 형성 ➔ 기울기 왜곡</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 완전 해석적 미분방정식 해 (Exact ODE):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">u(t) = u_0 e^(-βt) + (α/β)(1 - e^(-βt))</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자별 4개 반응 파라미터 최적화:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">θ_g = {α^(k), β, γ, t_on, t_off} 동시 추정</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#be123c" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.5" font-weight="600">장점: 부분적 스플라이싱 및 발현 펄스 유전자 완벽 복원</text>
  </g>

  <!-- Panel 2: EM Parameter Estimation -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 기댓값 최대화(EM) 최적화</text>

    <!-- Schematic: EM Loop -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- E-step box -->
      <rect width="105" height="30" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="52" y="15" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">E-단계: 잠재 시간</text>
      <text x="52" y="25" text-anchor="middle" fill="#94a3b8" font-size="6.8">t_i = argmin dist(c_i, γ(t))</text>

      <!-- Cycle arrows -->
      <path d="M 115 56 C 135 48, 145 48, 155 56" stroke="#fbbf24" stroke-width="1.5" fill="none" marker-end="url(#arrCyan03)"/>
      <path d="M 155 70 C 145 78, 135 78, 115 70" stroke="#c084fc" stroke-width="1.5" fill="none" marker-end="url(#arrPurple03)"/>

      <!-- M-step box -->
      <g transform="translate(160, 0)">
        <rect width="105" height="30" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
        <text x="52" y="15" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">M-단계: 모수 갱신</text>
        <text x="52" y="25" text-anchor="middle" fill="#94a3b8" font-size="6.8">max L(α, β, γ, t_on, t_off)</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 잠재 시간 직교 사상 (Orthogonal Projection):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">관측 세포 (u_i, s_i)에서 위상 곡선으로 최단 거리 시점 t_i 할당</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 전사 우도 극대화 (Likelihood Maximization):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">log L(θ) = -Σ_i || (u_i, s_i) - (u(t_i), s(t_i)) ||^2 / (2 σ^2)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자 신뢰도 필터링:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">우도비 검정(Likelihood Ratio)으로 동역학 유의 유전자 선택</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">수렴 안정성: EM 반복(기본 max_iter=10)으로 파라미터 최적화</text>
  </g>

  <!-- Panel 3: Global Latent Time tau -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">3. 전역 잠재 시간(τ) &amp; 계통 축</text>

    <!-- Schematic: Latent time gradient -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 52)">
      <!-- Gradient bar -->
      <defs>
        <linearGradient id="tauBarGrad" x1="0%" y1="0%" x2="100%" y2="0%">
          <stop offset="0%" stop-color="#38bdf8"/>
          <stop offset="50%" stop-color="#fbbf24"/>
          <stop offset="100%" stop-color="#f43f5e"/>
        </linearGradient>
      </defs>
      <rect width="250" height="18" rx="4" fill="url(#tauBarGrad)"/>
      <line x1="20" y1="28" x2="20" y2="34" stroke="#38bdf8" stroke-width="1.5"/>
      <text x="20" y="44" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">τ = 0.0 (AR+)</text>

      <line x1="125" y1="28" x2="125" y2="34" stroke="#fbbf24" stroke-width="1.5"/>
      <text x="125" y="44" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">τ 중간 (Amphi, 예시)</text>

      <line x1="230" y1="28" x2="230" y2="34" stroke="#f43f5e" stroke-width="1.5"/>
      <text x="230" y="44" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">τ = 1.0 (t-NEPC)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 다중 유전자 시간 앙상블 정렬:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">수백 개 유전자의 시간축 {t_ig}을 공통 전역 시간 τ_i로 통합</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 일관성 가중 결합 (Coherence Weighting):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">유전자 간 속도 상관도 기반 가중치 부여로 비생물학적 노이즈 제거</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 절대 생물학적 시간축 복원:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">DPT와 달리 거리 왜곡 없이 실제 전사 사이클 진행도 제공</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">결과: 계통 가소성 전이의 정확한 시간적 순서 및 분기 시점 특정</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_waddington_ot_sinkhorn_transport_coupling(filepath: str):
    """
    4. Waddington-OT Sinkhorn Transport Coupling
    Timepoint t1 -> t2 push-forward -> Entropic Kantorovich -> Sinkhorn scaling.
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

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">와딩턴 최적 수송(Waddington-OT) &amp; 엔트로피 싱크혼(Sinkhorn) 결합</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">시점별 확률 측도 밀어내기(Push-forward) ➔ 칸토로비치 수송 계획 행렬 Γ* ➔ 싱크혼 반복 초고속 스케일링</text>

  <!-- Panel 1: Multi-timepoint Measures & Growth Rates -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 다중 시점 측도 &amp; 증식/사멸 보정</text>

    <!-- Schematic: t1 to t2 cloud mapping -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Cloud t1 -->
      <circle cx="20" cy="22" r="3.5" fill="#38bdf8"/>
      <circle cx="32" cy="14" r="4.5" fill="#38bdf8"/>
      <circle cx="24" cy="36" r="3" fill="#38bdf8"/>
      <text x="25" y="54" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">시점 t_1 (μ)</text>

      <!-- Transport lines with cost C_ij -->
      <path d="M 38 20 C 70 15, 110 18, 145 22" stroke="#64748b" stroke-width="1.2" stroke-dasharray="2,2"/>
      <path d="M 38 30 C 70 35, 110 32, 145 35" stroke="#38bdf8" stroke-width="1.5"/>

      <!-- Cloud t2 -->
      <g transform="translate(130, 0)">
        <circle cx="25" cy="18" r="4" fill="#34d399"/>
        <circle cx="38" cy="28" r="5" fill="#34d399"/>
        <circle cx="20" cy="38" r="3.5" fill="#34d399"/>
        <text x="30" y="54" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">시점 t_2 (ν)</text>
      </g>

      <text x="88" y="14" text-anchor="middle" fill="#fbbf24" font-size="6.5">비용 C_ij = ||x_i - y_j||^2</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 비평형 세포 집단 동역학:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">시점 간 단순 매칭이 아닌 세포 분열(g_i) 및 사멸율 모델링</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전체 표현형 수송 비용 (Cost Matrix):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">C_ij = ||x_i - y_j||_2^2 (고차원 PCA 발현 공간 거리)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 와딩턴 지형(Landscape) 역학 재구성:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">파괴적 단일세포 측정의 시간적 단절을 수학적 수송으로 연결</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">핵심 원리: 최소 에너지 투입으로 전체 세포 집단 밀도 이동 설명</text>
  </g>

  <!-- Panel 2: Entropic Regularized Kantorovich Problem -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad04)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 엔트로피 정규화 칸토로비치 문제</text>

    <!-- Schematic: Coupling Matrix Gamma & Formula -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <rect width="90" height="58" rx="4" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="45" y="16" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">수송 계획 행렬 Γ*</text>
      <!-- Matrix grid representation -->
      <circle cx="25" cy="30" r="3" fill="#38bdf8" fill-opacity="0.8"/>
      <circle cx="45" cy="30" r="1.5" fill="#64748b"/>
      <circle cx="65" cy="30" r="4.5" fill="#34d399"/>
      <circle cx="25" cy="46" r="2" fill="#64748b"/>
      <circle cx="45" cy="46" r="5" fill="#f43f5e"/>
      <circle cx="65" cy="46" r="2.5" fill="#fbbf24"/>
    </g>

    <g transform="translate(130, 48)">
      <rect width="165" height="58" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="82" y="16" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">min_Γ &lt;C, Γ&gt; - ε H(Γ)</text>
      <text x="82" y="30" text-anchor="middle" fill="#cbd5e1" font-size="6.8">s.t.  Γ 1 = μ,  Γ^T 1 = ν</text>
      <text x="82" y="46" text-anchor="middle" fill="#34d399" font-size="6.8">H(Γ) = -Σ Γ_ij (log Γ_ij - 1)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 칸토로비치 완화 (Kantorovich Relaxation):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">1:1 단사 함수 한계를 극복하고 확률적 1:다 분기 허용</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 엔트로피 정규화의 수학적 효과:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">강한 볼록성(Strong convexity) 부여 ➔ 유일한 최적해 보장</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 깁스 커널 (Gibbs Kernel Form):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">최적해 구조: Γ_ij^* = u_i K_ij v_j,  K_ij = exp( -C_ij / ε )</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.5" font-weight="600">계산 혁신: O(N^3 log N) 선형 계획법을 O(N^2) 행렬 연산으로 단축</text>
  </g>

  <!-- Panel 3: Sinkhorn-Knopp Scaling Algorithm -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 싱크혼-크놉(Sinkhorn) 반복 알고리즘</text>

    <!-- Schematic: Sinkhorn Iteration -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <rect width="115" height="26" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="57" y="17" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">u^(k+1) = μ ./ (K v^(k))</text>

      <path d="M 125 61 C 145 54, 145 54, 155 61" stroke="#fbbf24" stroke-width="1.5" fill="none" marker-end="url(#arrCyan04)"/>
      <path d="M 155 72 C 145 78, 145 78, 125 72" stroke="#34d399" stroke-width="1.5" fill="none" marker-end="url(#arrEmerald04)"/>

      <g transform="translate(145, 0)">
        <rect width="115" height="26" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
        <text x="57" y="17" text-anchor="middle" fill="#34d399" font-size="8" font-weight="700">v^(k+1) = ν ./ (K^T u^(k+1))</text>
      </g>
      <text x="130" y="52" text-anchor="middle" fill="#94a3b8" font-size="6.8">교대 대각 스케일링 (Diagonal Scaling)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 교대 행렬-벡터 곱셈 (Matrix-Vector Ops):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">GPU 병렬 연산에 극도로 최적화된 대규모 수송 최적화</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 선조 및 후손 확률 쿼리 (Ancestor/Descendant):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">특정 t-NEPC 종말 세포의 기원 선조 클론을 t_1에서 역추적</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 엔잘루타마이드 치료 유도 계통 전이 검증:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">약물 노출 전 루미날 세포 중 어떤 서브클론이 저항성으로 가는지 규명</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">임상 적용: 치료 저항성 획득 클론의 조상 집단 및 유도 인자 사전 포착</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_cellrank2_multikernel_gpcca_macrostates(filepath: str):
    """
    5. CellRank 2 Multi-Kernel & G-PCCA Macrostates
    Multi-kernel fusion -> Transition matrix P -> Schur decomposition & G-PCCA macrostates.
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
    <marker id="arrPurple05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
    <marker id="arrEmerald05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">CellRank 2 멀티커널 융합 &amp; 일반화 PCCA+(G-PCCA) 거시상태</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">RNA 속도 + 의시간 + 실시간 최적수송 커널 결합 ➔ 비가역 전이 행렬 P ➔ 슈어(Schur) 분해 &amp; 정상 거시상태</text>

  <!-- Panel 1: Multi-Kernel Fusion -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 다중 커널 융합 (Multi-Kernel)</text>

    <!-- Schematic: Kernel Additions -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(20, 50)">
      <rect width="60" height="24" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="30" y="15" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">K_velo (w1)</text>

      <text x="68" y="16" fill="#cbd5e1" font-size="10" font-weight="700">+</text>

      <g transform="translate(76, 0)">
        <rect width="60" height="24" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
        <text x="30" y="15" text-anchor="middle" fill="#c084fc" font-size="7.5" font-weight="700">K_dpt (w2)</text>
      </g>

      <text x="144" y="16" fill="#cbd5e1" font-size="10" font-weight="700">+</text>

      <g transform="translate(152, 0)">
        <rect width="60" height="24" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
        <text x="30" y="15" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">K_ot (w3)</text>
      </g>

      <line x1="216" y1="12" x2="232" y2="12" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan05)"/>
      <text x="245" y="16" fill="#fbbf24" font-size="11" font-weight="700">P</text>

      <text x="125" y="42" text-anchor="middle" fill="#94a3b8" font-size="7.2">전이 행렬: P_ij ≥ 0,  Σ_j P_ij = 1 (행 정규화)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 단일 모달리티 한계 극복:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">RNA 속도의 국소 노이즈를 전역 DPT 및 OT 방향성으로 상호 보완</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 가중치 볼록 결합 (Convex Combination):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">P = w_1 T_velo + w_2 T_dpt + w_3 T_ot (Σ w_k = 1)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 세포별 신뢰도 기반 동적 적응:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">속도 벡터 품질이 낮은 미분화 세포는 DPT 가중치 자동 증대</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">견고성 달성: 생물학적 모순 신호 억제 및 안정적 궤적 형성</text>
  </g>

  <!-- Panel 2: Real Schur Decomposition & Spectrum -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad05)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 실수 슈어(Schur) 분해 &amp; 거시상태</text>

    <!-- Schematic: Schur Form -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <rect width="115" height="58" rx="4" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
      <text x="57" y="16" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">슈어 분해 (Schur Form)</text>
      <text x="57" y="32" text-anchor="middle" fill="#cbd5e1" font-size="7.5">P = Q R Q^T</text>
      <text x="57" y="48" text-anchor="middle" fill="#38bdf8" font-size="6.8">Q: 직교, R: 준상삼각</text>
    </g>

    <g transform="translate(155, 48)">
      <rect width="135" height="58" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="67" y="16" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">스펙트럼 갭 (Gap)</text>
      <!-- Real eigenvalues plot -->
      <line x1="15" y1="46" x2="120" y2="46" stroke="#475569" stroke-width="1"/>
      <circle cx="25" cy="28" r="3" fill="#38bdf8"/>
      <circle cx="45" cy="30" r="3" fill="#38bdf8"/>
      <circle cx="65" cy="31" r="3" fill="#38bdf8"/>
      <circle cx="85" cy="33" r="3" fill="#38bdf8"/>
      <line x1="95" y1="20" x2="95" y2="46" stroke="#f43f5e" stroke-width="1.2" stroke-dasharray="2,2"/>
      <circle cx="108" cy="42" r="2.5" fill="#64748b"/>
      <text x="67" y="55" text-anchor="middle" fill="#94a3b8" font-size="6.5">k = 6 거시상태 (n_states=6 지정, 예시)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 비가역(Non-reversible) 연산자의 안정적 분해:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">전이 행렬 P는 비대칭이므로 복소수 고유값 문제 내포</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 실수 슈어 벡터 (Real Schur Vectors):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">불변 부분공간(Invariant subspace)을 수치적으로 완벽 보존</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 일반화 페론 클러스터 분석 (G-PCCA):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">슈어 기저를 최적 심플렉스(Simplex) 상의 소속 확률로 사상</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">수학적 진보: 복소 고유벡터 섭동(Perturbation) 결함 원천 차단</text>
  </g>

  <!-- Panel 3: 4 Macrostates & Lineage Branches -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. mCRPC 4대 거시상태 분류</text>

    <!-- Schematic: 4 Macrostate Cards -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Root: AR+ Adeno -->
      <rect width="58" height="24" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="29" y="15" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">AR+ Adeno</text>
      <text x="29" y="32" text-anchor="middle" fill="#94a3b8" font-size="6">초기 상태</text>

      <path d="M 60 12 L 80 12" stroke="#38bdf8" stroke-width="1.4" marker-end="url(#arrCyan05)"/>

      <!-- Intermediate: Amphicrine -->
      <g transform="translate(85, 0)">
        <rect width="58" height="24" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
        <text x="29" y="15" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">Amphicrine</text>
        <text x="29" y="32" text-anchor="middle" fill="#94a3b8" font-size="6">과도 전이체</text>
      </g>

      <!-- Fork arrows -->
      <path d="M 145 12 L 165 4" stroke="#c084fc" stroke-width="1.4" marker-end="url(#arrPurple05)"/>
      <path d="M 145 12 L 165 20" stroke="#34d399" stroke-width="1.4" marker-end="url(#arrEmerald05)"/>

      <!-- Terminal 1: DNPC -->
      <g transform="translate(170, -6)">
        <rect width="55" height="20" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
        <text x="27" y="13" text-anchor="middle" fill="#c084fc" font-size="6.8" font-weight="700">DNPC</text>
      </g>

      <!-- Terminal 2: t-NEPC -->
      <g transform="translate(170, 16)">
        <rect width="55" height="20" rx="3" fill="#1e293b" stroke="#34d399" stroke-width="1"/>
        <text x="27" y="13" text-anchor="middle" fill="#34d399" font-size="6.8" font-weight="700">t-NEPC</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 초기(Initial) 및 종말(Terminal) 상태 자동 검증:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">마르코프 연쇄의 유출/유입 플럭스 비율로 생물학적 역할 분별</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 종말 거시상태 1: t-NEPC (ASCL1/SYP/CHGA):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">완전 신경내분비 역분화 (강한 분해 속도 및 비가역 트랩)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 종말 거시상태 2: DNPC (ETV4/FGFR1):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">AR 음성 / NE 음성 이중 음성 표현형 독립 분기</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">임상적 함의: 호르몬 내성 전립선암의 2대 가소성 종착역 정량 규명</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_fate_probabilities_and_lineage_drivers(filepath: str):
    """
    6. Fate Probabilities & Lineage Drivers
    Absorbing Markov chain -> Single-cell fate probability B = N R -> Driver correlation curve.
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

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="16.5" font-weight="700">흡수 마르코프 운명 확률(Fate Probabilities) &amp; 계통 구동자(Drivers)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">흡수 전이 행렬 분할 ➔ 기본 행렬 N &amp; 운명 흡수 확률 B = N R ➔ 전이 구동 마스터 전사인자 발굴</text>

  <!-- Panel 1: Absorbing Markov Chain Formulation -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 흡수 마르코프 연쇄 분할</text>

    <!-- Schematic: Matrix Block Partition -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>

    <g transform="translate(30, 48)">
      <!-- Partitioned matrix -->
      <rect width="105" height="58" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <line x1="65" y1="0" x2="65" y2="58" stroke="#475569" stroke-width="1" stroke-dasharray="2,2"/>
      <line x1="0" y1="36" x2="105" y2="36" stroke="#475569" stroke-width="1" stroke-dasharray="2,2"/>

      <text x="32" y="22" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">Q (과도기)</text>
      <text x="85" y="22" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">R (흡수)</text>
      <text x="32" y="49" text-anchor="middle" fill="#64748b" font-size="8.5">0</text>
      <text x="85" y="49" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">I (항등)</text>
    </g>

    <g transform="translate(150, 48)">
      <rect width="120" height="58" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="60" y="16" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">기본 행렬 N</text>
      <text x="60" y="30" text-anchor="middle" fill="#cbd5e1" font-size="7">N = (I - Q)^(-1)</text>
      <text x="60" y="46" text-anchor="middle" fill="#fbbf24" font-size="6.8">N = Σ_(k=0)^∞ Q^k</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 마르코프 정준형 분할 (Canonical Form):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">상태를 과도 전이 상태(Transient)와 흡수 상태(Absorbing)로 분리</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 기본 행렬 N의 물리적 의미:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">N_ij: 세포 i가 흡수되기 전 상태 j에 머무르는 기대 방문 횟수</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 완전 폐쇄형 대수해 (Closed-Form):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">몬테카를로 경로 무작위 샘플링 없이 역행렬 대수로 엄밀 산출</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">수학적 엄밀성: 무한 급수 전이 경로의 정확한 총합 도출</text>
  </g>

  <!-- Panel 2: Absorption Probability Matrix B -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad06)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 운명 흡수 확률 행렬 B = N R</text>

    <!-- Schematic: Single cell fate bar plot -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <rect width="110" height="58" rx="4" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="55" y="16" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">B = N · R</text>
      <text x="55" y="32" text-anchor="middle" fill="#cbd5e1" font-size="7.5">B_ij ∈ [0, 1]</text>
      <text x="55" y="48" text-anchor="middle" fill="#34d399" font-size="7">Σ_j B_ij = 1.0</text>
    </g>

    <g transform="translate(150, 48)">
      <rect width="140" height="58" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="70" y="14" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">단일세포 i 운명 확률</text>
      <!-- Horizontal stacked bar -->
      <rect x="15" y="22" width="85" height="14" rx="2" fill="#34d399"/>
      <rect x="100" y="22" width="25" height="14" rx="2" fill="#c084fc"/>
      <text x="57" y="33" text-anchor="middle" fill="#0f172a" font-size="7" font-weight="700">t-NEPC: 78% (예시)</text>
      <text x="112" y="33" text-anchor="middle" fill="#f8fafc" font-size="6.5">DN: 22%</text>
      <text x="70" y="49" text-anchor="middle" fill="#94a3b8" font-size="6.5">세포별 절대 분화 편향성</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 단일세포 수준의 흡수 운명 정량화:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">각 세포가 궁극적으로 t-NEPC에 도달할 절대 확률 산출</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 분기 가소성 윈도우 (Plasticity Window):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">B_t-NEPC ≈ 0.5인 임계 세포군이 치료 개입의 최적 타깃</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 엔잘루타마이드 치료 유도 가소성 측정:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">루미날 마커 발현 세포 중 t-NEPC 운명 기약 세포 사전 선별</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.5" font-weight="600">진단적 가치: 형태학적 변화 이전의 전사체 운명 결정 상태 포착</text>
  </g>

  <!-- Panel 3: Lineage Drivers Correlation -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 계통 구동 유전자(Drivers) 발굴</text>

    <!-- Schematic: Expression vs Fate Correlation curve -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <g transform="translate(25, 48)">
      <!-- Axes -->
      <line x1="15" y1="56" x2="115" y2="56" stroke="#475569" stroke-width="1"/>
      <line x1="15" y1="56" x2="15" y2="6" stroke="#475569" stroke-width="1"/>
      <text x="65" y="60" text-anchor="middle" fill="#94a3b8" font-size="6.5">B_t-NEPC 운명 확률 (0 ➔ 1)</text>

      <!-- Up curve (ASCL1) -->
      <path d="M 15 52 Q 55 50 110 12" stroke="#34d399" stroke-width="2" fill="none"/>
      <text x="110" y="10" fill="#34d399" font-size="7" font-weight="700">ASCL1 (r=0.92)</text>

      <!-- Down curve (AR) -->
      <path d="M 15 15 Q 55 20 110 54" stroke="#f43f5e" stroke-width="1.6" fill="none"/>
      <text x="75" y="38" fill="#f43f5e" font-size="6.5">AR (r=-0.89)</text>
    </g>

    <g transform="translate(170, 48)">
      <rect width="115" height="58" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <text x="57" y="14" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">Top t-NEPC Drivers</text>
      <text x="57" y="27" text-anchor="middle" fill="#cbd5e1" font-size="6.8">1. ASCL1 (r = 0.92, 예시)</text>
      <text x="57" y="38" text-anchor="middle" fill="#cbd5e1" font-size="6.8">2. INSM1 (r = 0.88)</text>
      <text x="57" y="49" text-anchor="middle" fill="#cbd5e1" font-size="6.8">3. EZH2  (r = 0.85)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 운명 상관 분석 (Fate Correlation):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">r_g = Corr( X_g, B_(:, target) ) (피어슨/스피어만 상관계수)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 핵심 신경내분비 구동자 규명:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">ASCL1, INSM1, EZH2, SOX2를 구동 후보로 예시</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이중 음성(DNPC) 구동자 규명:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">ETV4, FGFR1을 DNPC 구동 후보로 예시</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">신약 개발 타깃: 후성유전학적 리프로그래밍(EZH2/ASCL1) 표적 억제</text>
  </g>
  <text x="500" y="353" text-anchor="middle" fill="#64748b" font-size="9">※ 개념도: 확률·상관계수 수치는 설명용 예시이며 실측 결과가 아님 (toy 파이프라인 산출은 pipelines/results_pipeline_02 참조)</text>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_diffusion_maps_and_pseudotime_geometry.svg", make_01_diffusion_maps_and_pseudotime_geometry),
        ("02_rna_velocity_phase_portrait_kinetics.svg", make_02_rna_velocity_phase_portrait_kinetics),
        ("03_dynamical_scvelo_latent_time_optimization.svg", make_03_dynamical_scvelo_latent_time_optimization),
        ("04_waddington_ot_sinkhorn_transport_coupling.svg", make_04_waddington_ot_sinkhorn_transport_coupling),
        ("05_cellrank2_multikernel_gpcca_macrostates.svg", make_05_cellrank2_multikernel_gpcca_macrostates),
        ("06_fate_probabilities_and_lineage_drivers.svg", make_06_fate_probabilities_and_lineage_drivers),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 06 Manifold Trajectory publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

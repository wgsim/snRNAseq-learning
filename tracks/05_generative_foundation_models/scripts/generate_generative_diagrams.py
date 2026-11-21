#!/usr/bin/env python3
"""
tracks/05_generative_foundation_models/scripts/generate_generative_diagrams.py
=============================================================================
Generates publication-grade, high-resolution SVG diagrams for Track 05:
Generative Deep Learning & Foundation Models in Single-Cell Biology.

Covers:
1. scVI VAE Probabilistic Architecture & Reparameterization Trick
2. ELBO Mathematical Optimization & Regularization Tradeoff
3. scANVI Semi-Supervised Classification & scArches Surgery
4. Single-Cell Transformer Tokenization & Self-Attention Matrix
5. scGPT Multi-Task Pre-training (GEP, GEPC) & Downstream Adaptation
6. In Silico Virtual Perturbation & Prostate Plasticity Reversal
"""

import os
import xml.etree.ElementTree as ET


def make_01_scvi_vae_probabilistic_architecture(filepath: str):
    """
    1. scVI VAE Probabilistic Architecture
    Encoder q_phi(z|x,s) -> Reparameterization (mu, sigma) -> Latent z ->
    Decoder p_theta(x|z,s) -> Negative Binomial parameters (rho, theta, l).
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad01" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="encGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="100%" stop-color="#818cf8"/>
    </linearGradient>
    <linearGradient id="decGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#818cf8"/>
      <stop offset="100%" stop-color="#34d399"/>
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

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">scVI 심층 확률론적 생성 모델 아키텍처 &amp; 변분 추론 파이프라인</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">인코더 q_φ(z|x,s) ➔ 재파라미터화 트릭(μ, σ) ➔ 잠재 공간 z ➔ 디코더 p_θ(x|z,s) ➔ 음이항(NB) 파라미터(ρ, θ, l)</text>

  <!-- Panel 1: Inference Network (Encoder) -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad01)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 변분 인코더 &amp; 재파라미터화</text>

    <!-- Schematic: Encoder Flow -->
    <rect x="15" y="42" width="265" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <!-- Input box -->
    <g transform="translate(25, 52)">
      <rect width="52" height="32" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="26" y="15" text-anchor="middle" fill="#f8fafc" font-size="9" font-weight="700">x_n, s_n</text>
      <text x="26" y="26" text-anchor="middle" fill="#94a3b8" font-size="7">카운트+배치</text>
    </g>

    <!-- Arrow -->
    <line x1="82" y1="68" x2="104" y2="68" stroke="#38bdf8" stroke-width="1.5" marker-end="url(#arrCyan01)"/>

    <!-- Encoder MLP box -->
    <g transform="translate(108, 48)">
      <rect width="64" height="40" rx="4" fill="#0369a1" fill-opacity="0.4" stroke="#38bdf8" stroke-width="1"/>
      <text x="32" y="17" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">q_φ(z|x,s)</text>
      <text x="32" y="30" text-anchor="middle" fill="#cbd5e1" font-size="7.5">Dense + LN</text>
    </g>

    <!-- Arrow split -->
    <path d="M 176 68 L 194 58" stroke="#38bdf8" stroke-width="1.2" marker-end="url(#arrCyan01)"/>
    <path d="M 176 68 L 194 78" stroke="#38bdf8" stroke-width="1.2" marker-end="url(#arrCyan01)"/>

    <!-- Mu and Sigma -->
    <g transform="translate(198, 46)">
      <rect width="36" height="18" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="18" y="12" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">μ_z</text>
    </g>
    <g transform="translate(198, 70)">
      <rect width="36" height="18" rx="3" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
      <text x="18" y="12" text-anchor="middle" fill="#f43f5e" font-size="8.5" font-weight="700">log σ_z</text>
    </g>

    <!-- Sampling box -->
    <g transform="translate(242, 54)">
      <circle cx="16" cy="14" r="14" fill="#0f172a" stroke="#c084fc" stroke-width="1.5"/>
      <text x="16" y="12" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">ε~N(0,I)</text>
      <text x="16" y="22" text-anchor="middle" fill="#cbd5e1" font-size="7">z=μ+σε</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 입력 특징: 비정규화 원시 카운트 x_n ∈ ℕ^G</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">시퀀싱 깊이(l_n) 및 배치 지시자(s_n)를 명시적 조건부 결합</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 변분 사후분포 근사 (Variational Posterior):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">q_φ(z|x,s) = N(μ_φ(x,s), diag(σ_φ^2(x,s)))</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 재파라미터화 트릭 (Kingma &amp; Welling):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">확률적 샘플링 노드를 통과하는 역전파 기울기 전파 보장</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">수학적 특징: 몬테카를로 경로별(Pathwise) 저분산 미분</text>
  </g>

  <!-- Panel 2: Latent Manifold & Invariance -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad01)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 잠재 공간 정규화 &amp; 배치 조건화</text>

    <!-- Schematic: Latent Space UMAP & Prior -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>

    <!-- Prior Gaussian curve -->
    <g transform="translate(25, 48)">
      <path d="M 0 52 Q 22 52 35 15 Q 45 3 55 15 Q 68 52 90 52" fill="none" stroke="#fbbf24" stroke-width="1.5"/>
      <line x1="45" y1="3" x2="45" y2="52" stroke="#64748b" stroke-dasharray="2,2"/>
      <text x="45" y="62" text-anchor="middle" fill="#fbbf24" font-size="8">Prior p(z)=N(0,I)</text>
    </g>

    <!-- Arrow to manifold -->
    <line x1="122" y1="78" x2="148" y2="78" stroke="#c084fc" stroke-width="1.5" marker-end="url(#arrCyan01)"/>

    <!-- 2D Latent clustering schematic -->
    <g transform="translate(154, 46)">
      <rect width="140" height="64" rx="4" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <!-- Cluster 1: Luminal AR+ (Batch 1,2,3 mixed) -->
      <circle cx="35" cy="24" r="3.5" fill="#38bdf8"/>
      <circle cx="42" cy="28" r="3.5" fill="#fbbf24"/>
      <circle cx="28" cy="32" r="3.5" fill="#34d399"/>
      <circle cx="38" cy="36" r="3.5" fill="#f43f5e"/>
      <text x="35" y="48" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">Luminal AR+</text>

      <!-- Cluster 2: t-NEPC (Batch 1,2,3 mixed) -->
      <circle cx="105" cy="22" r="3.5" fill="#38bdf8"/>
      <circle cx="112" cy="26" r="3.5" fill="#fbbf24"/>
      <circle cx="98" cy="30" r="3.5" fill="#34d399"/>
      <circle cx="108" cy="35" r="3.5" fill="#f43f5e"/>
      <text x="105" y="48" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">t-NEPC</text>

      <text x="70" y="14" text-anchor="middle" fill="#94a3b8" font-size="6.5">배치 혼합 &amp; 생물학적 보존</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 잠재 차원 압축: z_n ∈ ℝ^d (d = 10~30)</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">세포 고유의 전사체 매니폴드 및 연속적 분화 상태 포착</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 배치 효과 완화 (Batch Conditioning):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">배치 s 조건화로 z의 기술적 변동 감소 (보장 아님)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 사후분포 정규화 (Kullback-Leibler):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">D_KL(q_φ(z|x,s) ∥ p(z))를 통해 표준정규분포로 밀집 제어</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">scIB 지표: iLISI(배치 통합) &amp; cLISI(라벨 보존)로 평가 — 값은 데이터 의존</text>
  </g>

  <!-- Panel 3: Generative Decoder -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad01)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 생성적 디코더 &amp; 음이항(NB) 모델</text>

    <!-- Schematic: Decoder Network -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>

    <!-- Latent Input -->
    <g transform="translate(25, 52)">
      <rect width="45" height="32" rx="4" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
      <text x="22" y="15" text-anchor="middle" fill="#c084fc" font-size="9" font-weight="700">z_n, s_n</text>
      <text x="22" y="26" text-anchor="middle" fill="#94a3b8" font-size="7">잠재+배치</text>
    </g>

    <!-- Arrow -->
    <line x1="74" y1="68" x2="96" y2="68" stroke="#34d399" stroke-width="1.5" marker-end="url(#arrEmerald01)"/>

    <!-- Decoder MLP -->
    <g transform="translate(100, 48)">
      <rect width="64" height="40" rx="4" fill="#065f46" fill-opacity="0.4" stroke="#34d399" stroke-width="1"/>
      <text x="32" y="17" text-anchor="middle" fill="#34d399" font-size="9" font-weight="700">p_θ(x|z,s)</text>
      <text x="32" y="30" text-anchor="middle" fill="#cbd5e1" font-size="7.5">Softmax MLP</text>
    </g>

    <!-- Arrow split -->
    <path d="M 168 68 L 186 56" stroke="#34d399" stroke-width="1.2" marker-end="url(#arrEmerald01)"/>
    <path d="M 168 68 L 186 80" stroke="#34d399" stroke-width="1.2" marker-end="url(#arrEmerald01)"/>

    <!-- Parameters output -->
    <g transform="translate(190, 44)">
      <rect width="96" height="22" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="48" y="15" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="700">ρ_ng (평균 발현율)</text>
    </g>
    <g transform="translate(190, 70)">
      <rect width="96" height="22" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
      <text x="48" y="15" text-anchor="middle" fill="#fbbf24" font-size="8.5" font-weight="700">θ_g (분산 매개변수)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자별 발현 비율: ρ_ng = Softmax(f_θ(z, s))</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">세포별 평균 기대값: μ_ng = l_n · ρ_ng (l_n: 시퀀싱 깊이)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 과분산 관측 분포 (Negative Binomial):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">x_ng ~ NB(Mean = μ_ng, Dispersion = θ_g)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 노이즈 제거 및 순수 발현 복원 (Denoising):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">기술적 드롭아웃 및 증폭 바이어스가 제거된 대립유전자 복원</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">출력: 배치 보정 발현량 ρ_ng &amp; 차등 발현(DE) 베이즈 팩터 산출</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_02_elbo_mathematical_optimization(filepath: str):
    """
    2. ELBO Mathematical Optimization
    Reconstruction loss (Log-Likelihood) vs Regularization loss (KL divergence)
    tradeoff curve, beta-annealing warmup, and ELBO maximization geometry.
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
    <marker id="arrRose02" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">증거 하한(ELBO) 수학적 목적함수 &amp; 정규화 손실 트레이드오프</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">재구성 우도 E_q[log p(x|z)] vs KL 발산 손실 D_KL(q||p) &amp; 사후분포 붕괴 방지 β-어닐링 웜업</text>

  <!-- Panel 1: Mathematical Decomposition & Jensen -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad02)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 증거 하한(ELBO) 수학적 분해</text>

    <!-- Equation Schematic Box -->
    <rect x="14" y="42" width="267" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <!-- Big Bar comparison -->
    <g transform="translate(24, 50)">
      <!-- Total Marginal Log-Likelihood log p(x) -->
      <rect x="0" y="0" width="246" height="20" rx="3" fill="#0284c7" fill-opacity="0.25" stroke="#38bdf8" stroke-width="1"/>
      <text x="123" y="14" text-anchor="middle" fill="#38bdf8" font-size="9" font-weight="700">주변 로그 우도: log p(x) (난해한 참값)</text>

      <!-- Subdivided bar: ELBO + KL -->
      <rect x="0" y="26" width="166" height="24" rx="3" fill="#059669" fill-opacity="0.35" stroke="#34d399" stroke-width="1"/>
      <text x="83" y="42" text-anchor="middle" fill="#34d399" font-size="9" font-weight="700">ELBO(φ, θ) [최대화]</text>

      <rect x="170" y="26" width="76" height="24" rx="3" fill="#f43f5e" fill-opacity="0.25" stroke="#f43f5e" stroke-width="1"/>
      <text x="208" y="42" text-anchor="middle" fill="#f43f5e" font-size="8.5" font-weight="700">D_KL(q∥p) ≥ 0</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 옌센 부등식(Jensen's Inequality) 적용:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">log p(x) = log ∫ p(x,z) dz ≥ E_q[log (p(x,z)/q(z))]</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 두 가지 상충 손실의 합 (Objective Formulation):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">ELBO = E_q[log p_θ(x|z)] - D_KL(q_φ(z|x) ∥ p(z))</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 참 사후분포와의 괴리 최소화:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">ELBO를 최대화하는 것은 D_KL(q_φ ∥ p(z|x))를 0으로 수렴시킴</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">증명: log p(x) = ELBO + D_KL(q(z|x) ∥ p(z|x)) ≥ ELBO</text>
  </g>

  <!-- Panel 2: Loss Tradeoff Curves & Beta-Annealing -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad02)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 손실 트레이드오프 &amp; β-어닐링 곡선</text>

    <!-- Schematic: Line Chart -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    
    <!-- Chart axes -->
    <g transform="translate(25, 48)">
      <line x1="20" y1="56" x2="260" y2="56" stroke="#475569" stroke-width="1"/>
      <line x1="20" y1="8" x2="20" y2="56" stroke="#475569" stroke-width="1"/>
      <text x="260" y="66" text-anchor="end" fill="#64748b" font-size="7">Epochs (0~100)</text>

      <!-- Curve 1: Reconstruction Error (Green, dropping) -->
      <path d="M 20 12 Q 80 40 160 48 T 255 51" fill="none" stroke="#34d399" stroke-width="2"/>
      <text x="160" y="44" fill="#34d399" font-size="7" font-weight="700">재구성 손실 (NB LogLik ↓)</text>

      <!-- Curve 2: KL Divergence (Rose, starting small then leveling) -->
      <path d="M 20 54 Q 60 52 120 38 T 255 32" fill="none" stroke="#f43f5e" stroke-width="2"/>
      <text x="210" y="26" fill="#f43f5e" font-size="7" font-weight="700">KL 발산 (D_KL ↑)</text>

      <!-- Curve 3: Beta schedule (Amber dashed, ramping 0 to 1) -->
      <path d="M 20 54 L 110 14 L 255 14" fill="none" stroke="#fbbf24" stroke-width="1.5" stroke-dasharray="3,2"/>
      <text x="115" y="11" fill="#fbbf24" font-size="7" font-weight="700">β-어닐링 웜업 (0 ➔ 1.0)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 사후분포 붕괴(Posterior Collapse) 방지:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">학습 초기에 과도한 KL 정규화로 z 정보가 소멸되는 참사 억제</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 선형 β-Warmup 스케줄러:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">β(t) = min(1.0, t / T_warmup), 초기엔 재구성 학습에 집중 유도</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 정규화 강도와 매니폴드 표현력의 균형:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">β=1 완전 베이지안 수렴 시 노이즈 강인성 및 보간 능력 극대화</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.5" font-weight="600">학습 전략: KL 페널티를 T_warmup 에포크 동안 점진 인가</text>
  </g>

  <!-- Panel 3: Stochastic Gradient Descent Geometry -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad02)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 확률적 최적화 곡면 &amp; 경사 전파</text>

    <!-- Schematic: Contour Map with SGD Trajectory -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    
    <g transform="translate(25, 46)">
      <!-- Concentric ellipses (Contour lines of ELBO surface) -->
      <ellipse cx="140" cy="33" rx="110" ry="26" fill="none" stroke="#334155" stroke-width="1"/>
      <ellipse cx="140" cy="33" rx="75" ry="18" fill="none" stroke="#475569" stroke-width="1"/>
      <ellipse cx="140" cy="33" rx="40" ry="10" fill="none" stroke="#059669" stroke-width="1"/>
      <circle cx="140" cy="33" r="3" fill="#34d399"/>
      <text x="140" y="29" text-anchor="middle" fill="#34d399" font-size="6.5" font-weight="700">최적 파라미터 (θ*, φ*)</text>

      <!-- SGD Trajectory -->
      <path d="M 45 12 L 65 24 L 80 18 L 105 28 L 120 29 L 138 33" fill="none" stroke="#38bdf8" stroke-width="1.8" marker-end="url(#arrCyan02)"/>
      <circle cx="45" cy="12" r="3" fill="#f43f5e"/>
      <text x="45" y="8" fill="#f43f5e" font-size="6.5" font-weight="700">초기화</text>
      <text x="80" y="42" fill="#38bdf8" font-size="6.5">Adam + Mini-batch SGD</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 파라미터 공동 최적화 (Joint Optimization):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">인코더 φ와 디코더 θ를 Adam(lr=1e-3, weight_decay=1e-6) 동시 갱신</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 미니배치 샘플링 (128~512 세포):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">O(N) 확장성으로 수백만 세포 단일세포 아틀라스 초고속 수렴</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 음이항 로그 우도 닫힌 형태 미분:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">log Γ(x + θ) - log Γ(θ) 감마 함수 기반 역전파 고속 계산</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">수렴 특성: 50 에포크 이내 Loss 평탄화 및 잠재 공간 안정화</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_03_scanvi_and_scarches_atlas_mapping(filepath: str):
    """
    3. scANVI and scArches Atlas Mapping
    scANVI semi-supervised label transfer + scArches architectural surgery
    attaching new query batches to a frozen reference atlas.
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

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">scANVI 준지도 라벨 전이 &amp; scArches 동적 수술적 아틀라스 매핑</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">희소 라벨 활용 혼합 모델 분류 + 참조 가중치 동결 및 쿼리 어댑터 수술을 통한 파괴적 망각 없는 지속적 통합</text>

  <!-- Panel 1: scANVI Semi-Supervised Model -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad03)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. scANVI 준지도 라벨 전이 엔진</text>

    <!-- Schematic: Dual Latent (z + y) -->
    <rect x="14" y="42" width="267" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <g transform="translate(22, 50)">
      <!-- Labeled data (30%) -->
      <rect x="0" y="0" width="70" height="24" rx="3" fill="#0369a1" fill-opacity="0.5" stroke="#38bdf8" stroke-width="1"/>
      <text x="35" y="12" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="700">30% 라벨 세포</text>
      <text x="35" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6.5">y 관측 (지도 학습)</text>

      <!-- Unlabeled data (70%) -->
      <rect x="0" y="32" width="70" height="24" rx="3" fill="#334155" stroke="#64748b" stroke-width="1"/>
      <text x="35" y="44" text-anchor="middle" fill="#cbd5e1" font-size="8" font-weight="700">70% 미지 세포</text>
      <text x="35" y="53" text-anchor="middle" fill="#94a3b8" font-size="6.5">y 잠재 (한계화 연산)</text>

      <!-- Arrows -->
      <line x1="72" y1="12" x2="102" y2="24" stroke="#38bdf8" stroke-width="1.2" marker-end="url(#arrCyan03)"/>
      <line x1="72" y1="44" x2="102" y2="34" stroke="#94a3b8" stroke-width="1.2" marker-end="url(#arrCyan03)"/>

      <!-- Classifier MLP q_psi(y|z) -->
      <rect x="104" y="16" width="65" height="30" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
      <text x="136" y="29" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">q_ψ(y|z)</text>
      <text x="136" y="40" text-anchor="middle" fill="#94a3b8" font-size="6.5">분류 네트워크</text>

      <!-- Arrow to output -->
      <line x1="171" y1="31" x2="188" y2="31" stroke="#c084fc" stroke-width="1.2" marker-end="url(#arrPurple03)"/>

      <!-- Transferred labels -->
      <rect x="190" y="15" width="58" height="32" rx="3" fill="#065f46" fill-opacity="0.5" stroke="#34d399" stroke-width="1"/>
      <text x="219" y="27" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">예측 라벨 ŷ</text>
      <text x="219" y="39" text-anchor="middle" fill="#cbd5e1" font-size="6.5">확률 사후분포</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 이중 잠재 구조: 연속 상태 z + 이산 라벨 y</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">p(z|y) ~ N(μ_y, Σ_y) 가우시안 혼합 모델(GMM) 사전분포</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 준지도 통합 목적함수 (Semi-Supervised ELBO):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">ℒ = ℒ_unlabeled(x) + α · ℒ_labeled(x,y) + β · ℒ_classifier(y,z)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 세포 유형 기반 가이던스:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">배치 효과를 강력히 제거하면서 생물학적 경계를 엄격히 보존</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">성능: 일부 라벨로 나머지 라벨 전이 (정확도는 데이터 의존)</text>
  </g>

  <!-- Panel 2: scArches Architectural Surgery -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad03)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. scArches 아키텍처 수술 (Surgery)</text>

    <!-- Schematic: Frozen weights + surgical adapter -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>
    
    <g transform="translate(25, 48)">
      <!-- Frozen Reference Network -->
      <rect x="0" y="0" width="130" height="36" rx="4" fill="#1e293b" stroke="#38bdf8" stroke-width="1.5"/>
      <rect x="3" y="3" width="124" height="12" rx="2" fill="#0284c7" fill-opacity="0.3"/>
      <text x="65" y="12" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">참조 아틀라스 가중치 W_ref</text>
      <text x="65" y="28" text-anchor="middle" fill="#fbbf24" font-size="8" font-weight="700">🔒 동결 (FROZEN - 고정)</text>

      <!-- Surgical Query Adapter -->
      <rect x="0" y="42" width="130" height="22" rx="4" fill="#065f46" fill-opacity="0.4" stroke="#34d399" stroke-width="1.5"/>
      <text x="65" y="56" text-anchor="middle" fill="#34d399" font-size="7.5" font-weight="700">✨ 쿼리 어댑터 가중치 W_query (수술)</text>

      <!-- Arrow -->
      <line x1="135" y1="28" x2="162" y2="28" stroke="#c084fc" stroke-width="1.5" marker-end="url(#arrPurple03)"/>

      <!-- Output Latent Mapping -->
      <rect x="165" y="6" width="95" height="52" rx="4" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
      <text x="212" y="22" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="700">통합 잠재 공간 z</text>
      <text x="212" y="34" text-anchor="middle" fill="#cbd5e1" font-size="6.5">Zero-Init 어댑터</text>
      <text x="212" y="46" text-anchor="middle" fill="#34d399" font-size="6.5">망각 최소화</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 아키텍처 수술 원리 (Architectural Surgery):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">신규 쿼리 배치가 들어오면 인코더 첫 층에 전용 노드/가중치만 접목</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 제로 초기화(Zero-Initialization) 기법:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">초기 W_query = 0 설정 ➔ 첫 스텝에서 참조 모델 출력과 100% 동일</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 파괴적 망각(Catastrophic Forgetting) 배제:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">수십만 세포 기준 아틀라스를 재학습하지 않고 신규 환자 즉시 투영</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">장점: 프라이버시 보호(참조 원시 데이터 공유 불필요) &amp; 초경량 학습</text>
  </g>

  <!-- Panel 3: Unified Atlas Integration & Novel Clone Discovery -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad03)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 통합 아틀라스 &amp; 신규 클론 발견</text>

    <!-- Schematic: UMAP with Reference Atlas + Query Biopsies -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    
    <g transform="translate(25, 46)">
      <!-- Main Reference cluster (Grey dots) -->
      <circle cx="45" cy="35" r="2.5" fill="#64748b"/>
      <circle cx="52" cy="40" r="2.5" fill="#64748b"/>
      <circle cx="38" cy="44" r="2.5" fill="#64748b"/>
      <circle cx="46" cy="48" r="2.5" fill="#64748b"/>
      <!-- Query 1 (Cyan dots) mapped seamlessly -->
      <circle cx="48" cy="38" r="3" fill="#38bdf8"/>
      <circle cx="42" cy="42" r="3" fill="#38bdf8"/>
      <text x="45" y="60" text-anchor="middle" fill="#38bdf8" font-size="6.5" font-weight="700">Luminal AR+ 매핑</text>

      <!-- Immune cluster -->
      <circle cx="115" cy="20" r="2.5" fill="#64748b"/>
      <circle cx="122" cy="25" r="2.5" fill="#64748b"/>
      <circle cx="118" cy="22" r="3" fill="#38bdf8"/>
      <text x="118" y="36" text-anchor="middle" fill="#38bdf8" font-size="6.5" font-weight="700">TAM/면역 세포</text>

      <!-- Out-of-Distribution Novel Cluster (Rose dots, separated) -->
      <circle cx="205" cy="38" r="3.5" fill="#f43f5e"/>
      <circle cx="214" cy="42" r="3.5" fill="#f43f5e"/>
      <circle cx="200" cy="46" r="3.5" fill="#f43f5e"/>
      <circle cx="210" cy="50" r="3.5" fill="#f43f5e"/>
      <rect x="175" y="10" width="75" height="18" rx="3" fill="#f43f5e" fill-opacity="0.25" stroke="#f43f5e" stroke-width="1"/>
      <text x="212" y="22" text-anchor="middle" fill="#f43f5e" font-size="6.5" font-weight="700">신규 t-NEPC 클론 발견!</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 기준 아틀라스 투영 (mCRPC Reference):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">골전이/원발암 10만 세포 기준 좌표계에 간전이 환자 생검 정합</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• OOD(분포 외) 신규 세포군 자동 탐지:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">기준 아틀라스에 없던 희귀 저분화/이중음성 클론을 독립 군집화</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 임상 라벨 정밀 전달:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">수작업 주석 대비 세포 주석 시간을 크게 단축</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">결과: 기준 좌표계 일관성 유지 + 신규 치료 저항 서브클론 탐색</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_04_single_cell_transformer_tokenization(filepath: str):
    """
    4. Single-Cell Transformer Tokenization
    Rank-value tokenization -> Multi-head self-attention matrix ->
    Context-aware gene and cell embeddings.
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
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">단일세포 트랜스포머 토큰화 &amp; 멀티헤드 자기주의(Self-Attention) 메커니즘</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">발현량 순위-값(Rank-Value) 토큰화 ➔ 자기주의 가중치 행렬 ➔ 문맥 인식 유전자 조절망 및 [CLS] 세포 임베딩</text>

  <!-- Panel 1: Rank-Value Tokenization -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad04)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 발현 순위-값(Rank-Value) 토큰화</text>

    <!-- Schematic: Gene sorting & embedding tokens -->
    <rect x="14" y="42" width="267" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <g transform="translate(22, 50)">
      <!-- Expression bar ranks -->
      <text x="0" y="12" fill="#94a3b8" font-size="7">발현량 내림차순 정렬:</text>
      <!-- Tokens -->
      <g transform="translate(0, 18)">
        <!-- Token 1: AR -->
        <rect x="0" y="0" width="46" height="36" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <text x="23" y="13" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">AR</text>
        <rect x="3" y="18" width="40" height="12" rx="2" fill="#0284c7" fill-opacity="0.3"/>
        <text x="23" y="27" text-anchor="middle" fill="#cbd5e1" font-size="6.5">Val: 31 (Rank 1)</text>

        <!-- Token 2: KLK3 -->
        <rect x="52" y="0" width="48" height="36" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <text x="76" y="13" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">KLK3</text>
        <rect x="55" y="18" width="42" height="12" rx="2" fill="#0284c7" fill-opacity="0.3"/>
        <text x="76" y="27" text-anchor="middle" fill="#cbd5e1" font-size="6.5">Val: 28 (Rank 2)</text>

        <!-- Token 3: EZH2 -->
        <rect x="106" y="0" width="48" height="36" rx="3" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
        <text x="130" y="13" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">EZH2</text>
        <rect x="109" y="18" width="42" height="12" rx="2" fill="#d97706" fill-opacity="0.3"/>
        <text x="130" y="27" text-anchor="middle" fill="#cbd5e1" font-size="6.5">Val: 24 (Rank 3)</text>

        <!-- Truncation / Pad -->
        <text x="162" y="22" fill="#64748b" font-size="10">···</text>

        <!-- Token K: Top 128 -->
        <rect x="180" y="0" width="58" height="36" rx="3" fill="#1e293b" stroke="#64748b" stroke-width="1"/>
        <text x="209" y="13" text-anchor="middle" fill="#cbd5e1" font-size="7.5">Top-K</text>
        <text x="209" y="27" text-anchor="middle" fill="#94a3b8" font-size="6.5">Non-Zero 유전자</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 스파스성 극복: 0-카운트 유전자 제외</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">세포당 1,000~2,000개 활성 발현 유전자만 선별하여 시퀀스 압축</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 3중 임베딩 합성 (Embedding Triplet):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">x_i = E_gene(ID) + E_value(Bin, B=51) + E_batch</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 순위·값 정보 보존 (모델마다 방식 상이):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">Geneformer: 순위 정렬 시퀀스 / scGPT: 값 구간 임베딩</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">어휘: 유전자 토큰 + 51 단계 발현량 구간 (scGPT 기준)</text>
  </g>

  <!-- Panel 2: Multi-Head Self-Attention Matrix -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad04)" stroke="#fbbf24" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#fbbf24" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#fbbf24" font-size="12.5" font-weight="700">2. 자기주의 행렬 &amp; 전사 조절망</text>

    <!-- Schematic: Attention Heatmap & Formula -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#d97706" stroke-width="1"/>
    
    <g transform="translate(25, 46)">
      <!-- 5x5 Attention Grid -->
      <rect x="0" y="0" width="65" height="65" fill="#1e293b" stroke="#334155" stroke-width="1"/>
      <!-- Strong weights (AR-KLK3 block) -->
      <rect x="0" y="0" width="26" height="26" fill="#f43f5e" fill-opacity="0.7"/>
      <rect x="26" y="0" width="13" height="13" fill="#fbbf24" fill-opacity="0.8"/>
      <!-- Strong weights (EZH2-ASCL1 block) -->
      <rect x="39" y="39" width="26" height="26" fill="#38bdf8" fill-opacity="0.75"/>
      <text x="13" y="16" text-anchor="middle" fill="#fff" font-size="6" font-weight="700">AR 모듈</text>
      <text x="52" y="55" text-anchor="middle" fill="#fff" font-size="6" font-weight="700">NE 모듈</text>

      <!-- Labels -->
      <text x="80" y="16" fill="#fbbf24" font-size="8" font-weight="700">Attn = Softmax(QK^T / √d_k) V</text>
      <text x="80" y="30" fill="#cbd5e1" font-size="7">• Q, K, V 투영 행렬 (Head = 8~16)</text>
      <text x="80" y="42" fill="#cbd5e1" font-size="7">• 유전자 간 전사 공발현 결합도 포착</text>
      <text x="80" y="54" fill="#34d399" font-size="7">• 동적 유전자 조절망(GRN) 자동 구축</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자-유전자 전역 상호작용 (All-to-All Interaction):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">전사인자(AR, EZH2)와 하류 타깃(KLK3, ASCL1) 간 주의 집중</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 문맥 의존적 전사 네트워크 해석:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">세포 상태(mCRPC vs t-NEPC)에 따라 자기주의 가중치 재구성</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 고속화 최적화 (FlashAttention &amp; Flash-Decoding):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">메모리 I/O 병목 해소로 O(N^2) 연산의 GPU 타일링 가속</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#d97706" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#fbbf24" font-size="9.5" font-weight="600">핵심 기능: 사전 지식 없이 전사체 발현만으로 인과 조절망 역공학</text>
  </g>

  <!-- Panel 3: Contextual Gene & Cell Embeddings -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad04)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 문맥 인식 임베딩 &amp; [CLS] 세포 벡터</text>

    <!-- Schematic: CLS Token & Gene vectors -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    
    <g transform="translate(25, 48)">
      <!-- CLS Token output -->
      <rect x="0" y="0" width="70" height="32" rx="3" fill="#047857" fill-opacity="0.5" stroke="#34d399" stroke-width="1.5"/>
      <text x="35" y="14" text-anchor="middle" fill="#34d399" font-size="8.5" font-weight="700">h_[CLS]</text>
      <text x="35" y="25" text-anchor="middle" fill="#cbd5e1" font-size="6.5">전체 세포 대표 벡터</text>

      <!-- Arrow -->
      <line x1="74" y1="16" x2="100" y2="16" stroke="#34d399" stroke-width="1.2" marker-end="url(#arrEmerald04)"/>

      <!-- Downstream tasks -->
      <rect x="104" y="0" width="165" height="32" rx="3" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="186" y="13" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">하류 과제: 세포 주석 &amp; 가상 섭동</text>
      <text x="186" y="24" text-anchor="middle" fill="#94a3b8" font-size="6.5">MLP Classifier / Cosine Similarity</text>

      <!-- Gene embeddings below -->
      <g transform="translate(0, 38)">
        <rect x="0" y="0" width="55" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <text x="27" y="12" text-anchor="middle" fill="#38bdf8" font-size="7">h_AR (문맥)</text>
        <rect x="60" y="0" width="55" height="18" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
        <text x="87" y="12" text-anchor="middle" fill="#fbbf24" font-size="7">h_EZH2 (문맥)</text>
        <text x="125" y="12" fill="#64748b" font-size="8">➔ 유전자 기능 클러스터링</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 단일세포 전역 표현 (Cell-Level Representation):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">시퀀스 선두 [CLS] 토큰 또는 전체 유전자 벡터 평균 풀링</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 동적 유전자 임베딩 (Contextual Gene Embeddings):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">동일한 EZH2 유전자라도 t-NEPC 문맥과 루미날 문맥에서 벡터 변이</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 다목적 전이 학습 역량 (Generalizability):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">512차원 잠재 공간에서 세포 상태 분류 및 약물 스크리닝 즉시 수행</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">차원: d_model = 512, Layer = 12 층의 다층 인코더 심층 표현</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_05_scgpt_multitask_pretraining_and_finetuning(filepath: str):
    """
    5. scGPT Multi-Task Pre-training and Fine-tuning
    Masked expression-value prediction (GEP) + cell-embedding value prediction (GEPC) ->
    Pre-trained foundation model -> Downstream fine-tuning.
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
    <marker id="arrPurple05" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#c084fc"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">scGPT 다중 작업 자기지도 사전학습 &amp; 하류 과제 미세조정(Fine-Tuning)</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">마스킹된 발현값 예측(GEP) + 셀 임베딩 기반 값 예측(GEPC) ➔ 대규모 파운데이션 백본 ➔ 세포 주석·배치 보정·유전자 조절망</text>

  <!-- Panel 1: Pre-training Objectives -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad05)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">1. 자기지도 사전학습 (Pre-training)</text>

    <!-- Schematic: MGM and GEP -->
    <rect x="14" y="42" width="267" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <g transform="translate(22, 48)">
      <!-- Input with Mask -->
      <g transform="translate(0, 4)">
        <rect x="0" y="0" width="38" height="24" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <text x="19" y="15" text-anchor="middle" fill="#38bdf8" font-size="7.5">AR</text>

        <!-- Masked Token -->
        <rect x="42" y="0" width="48" height="24" rx="2" fill="#f43f5e" fill-opacity="0.3" stroke="#f43f5e" stroke-width="1.5"/>
        <text x="66" y="15" text-anchor="middle" fill="#f43f5e" font-size="7.5" font-weight="700">[MASK]</text>

        <rect x="94" y="0" width="38" height="24" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
        <text x="113" y="15" text-anchor="middle" fill="#38bdf8" font-size="7.5">FOXA1</text>
      </g>

      <!-- Target Prediction -->
      <g transform="translate(142, 0)">
        <rect x="0" y="0" width="105" height="32" rx="3" fill="#0369a1" fill-opacity="0.4" stroke="#38bdf8" stroke-width="1"/>
        <text x="52" y="13" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">GEP 손실:</text>
        <text x="52" y="25" text-anchor="middle" fill="#34d399" font-size="7">Target = KLK3 발현값</text>
      </g>

      <!-- GEP Continuous regression -->
      <g transform="translate(0, 36)">
        <rect x="0" y="0" width="248" height="20" rx="2" fill="#1e293b" stroke="#fbbf24" stroke-width="1"/>
        <text x="124" y="14" text-anchor="middle" fill="#fbbf24" font-size="7.5" font-weight="700">GEPC 손실: 셀 임베딩 기반 값 회귀 (MSE)</text>
      </g>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 마스킹된 발현값 예측 (GEP):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">유전자 토큰은 주어지고 가려진 발현값을 문맥으로 예측</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 셀 임베딩 기반 값 예측 (GEPC):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">세포 임베딩이 발현 크기 정보를 담도록 유도</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 다중 작업 결합 목적함수:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">ℒ_pretrain = ℒ_GEP + λ ℒ_GEPC (+ ECS 등)</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">데이터: 수천만 정상 인간 세포 (공개 체크포인트)</text>
  </g>

  <!-- Panel 2: Pre-trained Foundation Backbone -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad05)" stroke="#c084fc" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#c084fc" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#c084fc" font-size="12.5" font-weight="700">2. 대규모 트랜스포머 백본 구조</text>

    <!-- Schematic: Stack of Transformer Blocks -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#7e22ce" stroke-width="1"/>
    
    <g transform="translate(25, 48)">
      <!-- Layer 1 -->
      <rect x="0" y="38" width="265" height="18" rx="3" fill="#1e293b" stroke="#64748b" stroke-width="1"/>
      <text x="132" y="50" text-anchor="middle" fill="#94a3b8" font-size="7">Layer 1: 어휘 임베딩 및 기초 조절 결합</text>

      <!-- Middle dots -->
      <text x="132" y="32" text-anchor="middle" fill="#c084fc" font-size="9">· · · (12 인코더 레이어) · · ·</text>

      <!-- Layer 12 -->
      <rect x="0" y="2" width="265" height="18" rx="3" fill="#581c87" fill-opacity="0.5" stroke="#c084fc" stroke-width="1.2"/>
      <text x="132" y="14" text-anchor="middle" fill="#c084fc" font-size="7.5" font-weight="700">Layer 12: 고차원 세포 운명 및 전사 상태 공간 (d=512)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 생물학적 유전체 문법 (Biological Syntax) 체득:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">조직 간 불변 전사 규칙 및 세포 주기, 대사 경로 규칙 내재화</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 모델 파라미터 스케일링:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">BERT/GPT형 트랜스포머 구조를 단일세포에 전이</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• Zero-Shot 일반화 능력:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">제로샷 활용은 가능하나 신뢰성은 벤치마크로 검증 필요</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#7e22ce" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#c084fc" font-size="9.5" font-weight="600">아키텍처: FlashAttention 탑재 12계층 트랜스포머 인코더</text>
  </g>

  <!-- Panel 3: Downstream Adaptation & Transfer Tasks -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad05)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 하류 과제 전이 &amp; 미세조정(PEFT)</text>

    <!-- Schematic: Task branching -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    
    <g transform="translate(25, 48)">
      <!-- Task 1: Annotation -->
      <rect x="0" y="0" width="80" height="26" rx="3" fill="#047857" fill-opacity="0.4" stroke="#34d399" stroke-width="1"/>
      <text x="40" y="11" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">과제 1: 세포 주석</text>
      <text x="40" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6">Linear Head</text>

      <!-- Task 2: Integration -->
      <rect x="90" y="0" width="80" height="26" rx="3" fill="#0284c7" fill-opacity="0.4" stroke="#38bdf8" stroke-width="1"/>
      <text x="130" y="11" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">과제 2: 배치 통합</text>
      <text x="130" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6">Batch Token 교정</text>

      <!-- Task 3: In Silico Perturbation -->
      <rect x="180" y="0" width="85" height="26" rx="3" fill="#d97706" fill-opacity="0.4" stroke="#fbbf24" stroke-width="1"/>
      <text x="222" y="11" text-anchor="middle" fill="#fbbf24" font-size="7" font-weight="700">과제 3: 가상 섭동</text>
      <text x="222" y="21" text-anchor="middle" fill="#cbd5e1" font-size="6">디지털 녹아웃 예측</text>

      <!-- LoRA / Adapter badge -->
      <rect x="35" y="34" width="195" height="20" rx="3" fill="#1e293b" stroke="#c084fc" stroke-width="1"/>
      <text x="132" y="47" text-anchor="middle" fill="#c084fc" font-size="7.5" font-weight="700">⚡ LoRA 어댑터 미세조정: 어댑터 파라미터만 갱신</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 파라미터 효율적 미세조정 (PEFT / LoRA):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">기반 백본을 동결하고 저순위 어댑터만 조정하여 단일 GPU로 완수</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 소량 표본 전이 학습 (Few-Shot Transfer):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">라벨 세포가 100개 미만인 희귀 암 서브클론도 정밀 분류</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 유전자 조절망(GRN) 동적 추출:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">파인튜닝 전후 어텐션 맵 차이를 분석해 약물 작용 표적 동정</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">효율성: 풀 미세조정 대비 학습 파라미터·메모리 절감</text>
  </g>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def make_06_in_silico_perturbation_and_lineage_shift(filepath: str):
    """
    6. In Silico Perturbation and Lineage Shift
    Virtual knockout of master regulators (EZH2, ASCL1) ->
    Transformer latent vector displacement Delta z ->
    t-NEPC to therapy-sensitive reprogramming trajectory.
    """
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" width="1000" height="360" style="background:#0f172a; border-radius:12px; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
  <defs>
    <linearGradient id="cardGrad06" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#1e293b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <marker id="arrRose06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#f43f5e"/>
    </marker>
    <marker id="arrEmerald06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#34d399"/>
    </marker>
    <marker id="arrCyan06" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">
      <path d="M0,0 L6,3 L0,6 Z" fill="#38bdf8"/>
    </marker>
  </defs>

  <!-- Title -->
  <text x="500" y="28" text-anchor="middle" fill="#38bdf8" font-size="17" font-weight="700">인 실리코(In Silico) 가상 유전자 섭동 &amp; 전립선암 계통 가소성 역전환</text>
  <text x="500" y="47" text-anchor="middle" fill="#94a3b8" font-size="11.5">EZH2·ASCL1 가상 결손 ➔ 잠재 공간 변위 벡터 Δz 산출 ➔ 악성 t-NEPC에서 치료 감수성 루미날 분화 회귀</text>

  <!-- Panel 1: Digital Knockout Protocol -->
  <g transform="translate(20, 65)">
    <rect width="295" height="275" rx="10" fill="url(#cardGrad06)" stroke="#f43f5e" stroke-width="1.4"/>
    <rect x="0" y="0" width="295" height="32" rx="10" fill="#f43f5e" fill-opacity="0.16"/>
    <text x="147" y="21" text-anchor="middle" fill="#f43f5e" font-size="12.5" font-weight="700">1. 디지털 유전자 녹아웃 (Virtual KO)</text>

    <!-- Schematic: Vector clamping -->
    <rect x="14" y="42" width="267" height="74" rx="6" fill="#0f172a" stroke="#be123c" stroke-width="1"/>
    
    <g transform="translate(22, 48)">
      <!-- Baseline Profile -->
      <g transform="translate(0, 0)">
        <text x="0" y="10" fill="#94a3b8" font-size="7">기준 t-NEPC 세포 벡터:</text>
        <rect x="0" y="14" width="90" height="20" rx="2" fill="#1e293b" stroke="#f43f5e" stroke-width="1"/>
        <text x="45" y="27" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">EZH2 = High (Val 31)</text>
      </g>

      <!-- Digital Clamping Action -->
      <line x1="95" y1="24" x2="135" y2="24" stroke="#f43f5e" stroke-width="1.5" marker-end="url(#arrRose06)"/>
      <text x="115" y="18" text-anchor="middle" fill="#fbbf24" font-size="6.5">KO: 0 강제</text>

      <!-- Perturbed Profile -->
      <g transform="translate(140, 0)">
        <text x="0" y="10" fill="#94a3b8" font-size="7">가상 섭동 세포 벡터:</text>
        <rect x="0" y="14" width="95" height="20" rx="2" fill="#065f46" fill-opacity="0.5" stroke="#34d399" stroke-width="1"/>
        <text x="47" y="27" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">EZH2 = 0 (Val 0)</text>
      </g>

      <!-- Formula below -->
      <rect x="0" y="40" width="245" height="18" rx="2" fill="#1e293b" stroke="#38bdf8" stroke-width="1"/>
      <text x="122" y="52" text-anchor="middle" fill="#38bdf8" font-size="7.5" font-weight="700">변위 벡터: Δz = z_pert(EZH2=0) - z_baseline</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 디지털 유전체 간섭 (In Silico Intervention):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">CRISPR 스크리닝 없이 입력 토큰의 발현 순위/값을 0으로 강제 클램핑</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 다중 타깃 조합 시뮬레이션 (Combinatorial):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">단일 EZH2 KO vs ASCL1 KO vs Dual (AR 재활성화 + EZH2 억제)</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 피드포워드 고속 추론 (Zero-Latency):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">사전학습 백본의 순전파 1회만으로 전체 유전자 반응 동시 예측</text>

    <!-- Badge -->
    <rect x="14" y="242" width="267" height="22" rx="4" fill="#be123c" fill-opacity="0.22"/>
    <text x="147" y="257" text-anchor="middle" fill="#f43f5e" font-size="9.5" font-weight="600">스크리닝: 순전파 1회로 가상 결손 반응 예측 (속도는 모델·HW 의존)</text>
  </g>

  <!-- Panel 2: Latent Vector Displacement & Velocity -->
  <g transform="translate(330, 65)">
    <rect width="320" height="275" rx="10" fill="url(#cardGrad06)" stroke="#38bdf8" stroke-width="1.4"/>
    <rect x="0" y="0" width="320" height="32" rx="10" fill="#38bdf8" fill-opacity="0.16"/>
    <text x="160" y="21" text-anchor="middle" fill="#38bdf8" font-size="12.5" font-weight="700">2. 잠재 변위 벡터 Δz &amp; 계통 궤적</text>

    <!-- Schematic: 2D Manifold Vector Shift -->
    <rect x="15" y="42" width="290" height="74" rx="6" fill="#0f172a" stroke="#0284c7" stroke-width="1"/>
    
    <g transform="translate(25, 46)">
      <!-- Target Attractor: Luminal AR+ (Top right) -->
      <circle cx="215" cy="20" r="14" fill="#0369a1" fill-opacity="0.3" stroke="#38bdf8" stroke-width="1.5"/>
      <text x="215" y="22" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">루미날 AR+</text>
      <text x="215" y="30" text-anchor="middle" fill="#cbd5e1" font-size="5.5">(치료 감수성)</text>

      <!-- Source Attractor: t-NEPC (Bottom left) -->
      <circle cx="45" cy="50" r="14" fill="#881337" fill-opacity="0.4" stroke="#f43f5e" stroke-width="1.5"/>
      <text x="45" y="52" text-anchor="middle" fill="#f43f5e" font-size="7" font-weight="700">t-NEPC</text>
      <text x="45" y="60" text-anchor="middle" fill="#cbd5e1" font-size="5.5">(치료 저항성)</text>

      <!-- Trajectory 1: EZH2 KO (Cyan vector) -->
      <path d="M 58 44 Q 105 28 140 28" fill="none" stroke="#38bdf8" stroke-width="2" marker-end="url(#arrCyan06)"/>
      <text x="100" y="24" fill="#38bdf8" font-size="6.5" font-weight="700">EZH2 KO (중간 가소성)</text>

      <!-- Trajectory 2: Dual KO (Emerald vector) -->
      <path d="M 58 42 L 195 22" fill="none" stroke="#34d399" stroke-width="2" stroke-dasharray="3,2" marker-end="url(#arrEmerald06)"/>
      <text x="125" y="44" fill="#34d399" font-size="6.5" font-weight="700">Dual AR+EZH2 (완전 회귀)</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 코사인 유사도 벡터 방향성 분석:</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">cos(Δz, v_luminal): 목표 분화 경로와의 정렬도 (개념 예시)</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 후성유전학적 와딩턴(Waddington) 에너지 극복:</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">EZH2(PRC2 복합체) 억제로 크로마틴 응축 장벽을 해제하고 계통 전환</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 가상 세포 속도(In Silico Velocity):</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">RNA velocity 원리를 트랜스포머 잠재 공간 미분으로 일반화</text>

    <!-- Badge -->
    <rect x="14" y="242" width="292" height="22" rx="4" fill="#0284c7" fill-opacity="0.22"/>
    <text x="160" y="257" text-anchor="middle" fill="#38bdf8" font-size="9.5" font-weight="600">계통 변위 평가: 목표 경로와의 코사인 정렬로 정량화</text>
  </g>

  <!-- Panel 3: Plasticity Suppression & Resensitization -->
  <g transform="translate(665, 65)">
    <rect width="315" height="275" rx="10" fill="url(#cardGrad06)" stroke="#34d399" stroke-width="1.4"/>
    <rect x="0" y="0" width="315" height="32" rx="10" fill="#34d399" fill-opacity="0.16"/>
    <text x="157" y="21" text-anchor="middle" fill="#34d399" font-size="12.5" font-weight="700">3. 가소성 점수 억제 &amp; 치료 재감작</text>

    <!-- Schematic: Quantitative Bar Chart -->
    <rect x="15" y="42" width="285" height="74" rx="6" fill="#0f172a" stroke="#059669" stroke-width="1"/>
    
    <g transform="translate(25, 48)">
      <!-- Baseline Bar: NE score 100% -->
      <rect x="0" y="0" width="180" height="15" rx="2" fill="#f43f5e" fill-opacity="0.8"/>
      <text x="5" y="11" fill="#fff" font-size="6.5" font-weight="700">대조군 t-NEPC (SYP, CHGA, ASCL1 고발현): 100% [예시]</text>

      <!-- EZH2 KO Bar: NE score 38% -->
      <rect x="0" y="18" width="68" height="15" rx="2" fill="#fbbf24" fill-opacity="0.85"/>
      <text x="5" y="29" fill="#0f172a" font-size="6.5" font-weight="700">가상 EZH2 KO: 38% (-62%) [예시]</text>

      <!-- Dual KO Bar: NE score 18% -->
      <rect x="0" y="36" width="32" height="15" rx="2" fill="#34d399" fill-opacity="0.9"/>
      <text x="38" y="47" fill="#34d399" font-size="6.5" font-weight="700">Dual KO/AR: 18% (-82%) [예시]</text>

      <text x="210" y="28" text-anchor="middle" fill="#38bdf8" font-size="7" font-weight="700">엔잘루타마이드</text>
      <text x="210" y="39" text-anchor="middle" fill="#34d399" font-size="7" font-weight="700">감수성 회복: 실험 검증 필요</text>
    </g>

    <!-- Bullets -->
    <text x="14" y="132" fill="#cbd5e1" font-size="10.5" font-weight="600">• 신경내분비 시그니처 억제 (Suppression):</text>
    <text x="24" y="148" fill="#94a3b8" font-size="9.5">SYP, CHGA, ASCL1 점수 감소 및 루미날 마커(KLK3, FOXA1) 회복</text>

    <text x="14" y="168" fill="#cbd5e1" font-size="10.5" font-weight="600">• 2차 호르몬 요법 재감작 (Resensitization):</text>
    <text x="24" y="184" fill="#94a3b8" font-size="9.5">엔잘루타마이드(Enzalutamide) 저항성 극복 여부는 실험적 검증 필요</text>

    <text x="14" y="204" fill="#cbd5e1" font-size="10.5" font-weight="600">• 신약 개발 병용 전략 확립:</text>
    <text x="24" y="220" fill="#94a3b8" font-size="9.5">EZH2 억제제(타제메토스타트) + AR 표적 병용 치료의 전임상 이론 검증</text>

    <!-- Badge -->
    <rect x="14" y="242" width="287" height="22" rx="4" fill="#047857" fill-opacity="0.22"/>
    <text x="157" y="257" text-anchor="middle" fill="#34d399" font-size="9.5" font-weight="600">임상적 의의: 난치성 전이암 계통 가소성 차단 및 표적 치료 재개</text>
  </g>
  <text x="500" y="353" text-anchor="middle" fill="#64748b" font-size="9">※ 개념도: 막대·수치는 설명용 예시이며 실측 결과가 아님 (toy 파이프라인 산출은 pipelines/results_pipeline_02 참조)</text>
</svg>'''
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(svg)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo_root, "docs", "images")
    os.makedirs(outdir, exist_ok=True)

    generators = [
        ("01_scvi_vae_probabilistic_architecture.svg", make_01_scvi_vae_probabilistic_architecture),
        ("02_elbo_mathematical_optimization.svg", make_02_elbo_mathematical_optimization),
        ("03_scanvi_and_scarches_atlas_mapping.svg", make_03_scanvi_and_scarches_atlas_mapping),
        ("04_single_cell_transformer_tokenization.svg", make_04_single_cell_transformer_tokenization),
        ("05_scgpt_multitask_pretraining_and_finetuning.svg", make_05_scgpt_multitask_pretraining_and_finetuning),
        ("06_in_silico_perturbation_and_lineage_shift.svg", make_06_in_silico_perturbation_and_lineage_shift),
    ]

    for fname, func in generators:
        target = os.path.join(outdir, fname)
        func(target)
        # Parse XML to verify syntax validity
        ET.parse(target)
        print(f"[✓] Successfully generated and verified XML: {target}")

    print(f"\n[✓] All 6 Track 05 Generative Foundation Models publication-grade SVG diagrams successfully generated in {outdir}")


if __name__ == "__main__":
    main()

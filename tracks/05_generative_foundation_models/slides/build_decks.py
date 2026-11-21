import os

svg_dir = 'tracks/05_generative_foundation_models/docs/images'

def get_svg(filename):
    with open(os.path.join(svg_dir, filename), 'r', encoding='utf-8') as f:
        return f.read().strip()

svg01 = get_svg('01_scvi_vae_probabilistic_architecture.svg')
svg02 = get_svg('02_elbo_mathematical_optimization.svg')
svg03 = get_svg('03_scanvi_and_scarches_atlas_mapping.svg')
svg04 = get_svg('04_single_cell_transformer_tokenization.svg')
svg05 = get_svg('05_scgpt_multitask_pretraining_and_finetuning.svg')
svg06 = get_svg('06_in_silico_perturbation_and_lineage_shift.svg')

css = """
  :root {
    --bg-primary: #0f172a;
    --bg-secondary: #1e293b;
    --card-bg: rgba(30, 41, 59, 0.85);
    --accent: #38bdf8;
    --accent-glow: rgba(56, 189, 248, 0.25);
    --text-main: #f8fafc;
    --text-muted: #94a3b8;
    --highlight: #f43f5e;
    --success: #34d399;
    --warning: #fbbf24;
    --purple: #c084fc;
    --code-bg: #090d16;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }
  body { background: var(--bg-primary); color: var(--text-main); overflow: hidden; height: 100vh; display: flex; flex-direction: column; }
  
  /* Slide Container */
  #presentation { position: relative; flex: 1; width: 100%; height: 100%; overflow: hidden; }
  .slide { position: absolute; top: 0; left: 0; width: 100%; height: 100%; opacity: 0; visibility: hidden; transition: opacity 0.35s ease, transform 0.35s ease; transform: scale(0.97); display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 40px 60px; }
  .slide.active { opacity: 1; visibility: visible; transform: scale(1); }
  
  /* Slide Layout */
  .slide-content { max-width: 1200px; width: 100%; max-height: 85vh; overflow-y: auto; background: var(--card-bg); border: 1px solid rgba(255,255,255,0.08); border-radius: 16px; padding: 36px 44px; box-shadow: 0 20px 50px rgba(0,0,0,0.5); backdrop-filter: blur(12px); }
  .slide-header { border-bottom: 2px solid var(--accent-glow); padding-bottom: 14px; margin-bottom: 20px; display: flex; justify-content: space-between; align-items: baseline; }
  .slide-header h2 { font-size: 1.85rem; color: var(--accent); font-weight: 700; letter-spacing: -0.5px; }
  .slide-header .tag { font-size: 0.85rem; padding: 4px 12px; border-radius: 20px; background: rgba(56, 189, 248, 0.15); color: var(--accent); font-weight: 600; text-transform: uppercase; }
  .doc-ref { font-size: 0.85rem; color: #38bdf8; background: rgba(56, 189, 248, 0.12); border: 1px solid rgba(56, 189, 248, 0.3); border-radius: 6px; padding: 5px 12px; margin-bottom: 14px; display: inline-flex; align-items: center; gap: 8px; font-weight: 600; }
  
  /* Typography & Components */
  p, li { font-size: 1.02rem; line-height: 1.62; color: #cbd5e1; margin-bottom: 10px; }
  strong { color: #ffffff; }
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 16px; }
  .grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 18px; margin-top: 16px; }
  .grid-4 { display: grid; grid-template-columns: 1fr 1fr 1fr 1fr; gap: 16px; margin-top: 16px; }
  .card { background: rgba(15, 23, 42, 0.6); border: 1px solid rgba(255,255,255,0.06); border-radius: 12px; padding: 18px 20px; }
  .card.accent { border-left: 4px solid var(--accent); }
  .card.highlight { border-left: 4px solid var(--highlight); }
  .card.warning { border-left: 4px solid var(--warning); }
  .card.success { border-left: 4px solid var(--success); }
  .card.purple { border-left: 4px solid var(--purple); }
  .card h3 { font-size: 1.18rem; color: #fff; margin-bottom: 8px; }
  
  table { width: 100%; border-collapse: collapse; margin-top: 14px; font-size: 0.92rem; }
  th, td { padding: 10px 14px; text-align: left; border-bottom: 1px solid rgba(255,255,255,0.08); }
  th { background: rgba(56, 189, 248, 0.1); color: var(--accent); font-weight: 600; }
  
  pre { background: var(--code-bg); border: 1px solid rgba(255,255,255,0.1); border-radius: 8px; padding: 14px 16px; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size: 0.88rem; color: #38bdf8; overflow-x: auto; margin-top: 12px; }
  code { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 4px; color: #38bdf8; font-size: 0.92em; }
  
  /* Semantic Math Formatting */
  .math-block { background: rgba(15, 23, 42, 0.7); border: 1px solid rgba(56, 189, 248, 0.25); border-radius: 8px; padding: 8px 14px; margin: 8px 0; text-align: center; font-family: "Cambria Math", "Latin Modern Math", "STIX Two Math", Georgia, serif; font-size: 1.08rem; color: #38bdf8; }
  .math-frac { display: inline-flex; flex-direction: column; vertical-align: middle; text-align: center; padding: 0 4px; font-size: 0.95em; line-height: 1.2; }
  .math-frac .num { border-bottom: 1.5px solid currentColor; padding-bottom: 2px; }
  .math-frac .den { padding-top: 2px; }

  /* Footer Controls */
  #controls { height: 60px; background: var(--bg-secondary); border-top: 1px solid rgba(255,255,255,0.08); display: flex; justify-content: space-between; align-items: center; padding: 0 30px; }
  .nav-btn { background: rgba(255,255,255,0.08); border: 1px solid rgba(255,255,255,0.12); color: #fff; padding: 8px 18px; border-radius: 8px; cursor: pointer; font-size: 0.9rem; transition: all 0.2s; }
  .nav-btn:hover { background: var(--accent); color: var(--bg-primary); }
  #counter { font-size: 0.95rem; color: var(--text-muted); font-weight: 600; }
  .progress-bar { position: absolute; top: 0; left: 0; height: 3px; background: var(--accent); width: 0%; transition: width 0.3s ease; }

  /* Diagram Drawer */
  .diagram-drawer {
    margin-top: 14px;
    background: rgba(15, 23, 42, 0.75);
    border: 1px solid rgba(56, 189, 248, 0.3);
    border-radius: 10px;
    padding: 10px 16px;
  }
  .diagram-drawer summary {
    cursor: pointer;
    font-weight: 600;
    color: var(--accent);
    font-size: 0.92rem;
  }
  .diagram-drawer .svg-container svg {
    width: 100%;
    height: auto;
    max-height: 275px;
    margin-top: 10px;
    border-radius: 8px;
  }
"""

js_code = """
  const slides = document.querySelectorAll('.slide');
  const counter = document.getElementById('counter');
  const progressBar = document.getElementById('progressBar');
  let current = 0;

  function showSlide(index) {
    slides.forEach((s, i) => {
      s.classList.toggle('active', i === index);
    });
    counter.textContent = `Slide ${index + 1} of ${slides.length}`;
    progressBar.style.width = `${((index + 1) / slides.length) * 100}%`;
  }

  function nextSlide() {
    if (current < slides.length - 1) {
      current++;
      showSlide(current);
    }
  }

  function prevSlide() {
    if (current > 0) {
      current--;
      showSlide(current);
    }
  }

  document.getElementById('nextBtn').addEventListener('click', nextSlide);
  document.getElementById('prevBtn').addEventListener('click', prevSlide);

  document.addEventListener('keydown', (e) => {
    if (e.key === ' ') e.preventDefault();
    if (e.key === 'ArrowRight' || e.key === ' ' || e.key === 'PageDown') {
      nextSlide();
    } else if (e.key === 'ArrowLeft' || e.key === 'PageUp') {
      prevSlide();
    } else if (e.key === 'f' || e.key === 'F') {
      if (!document.fullscreenElement) {
        document.documentElement.requestFullscreen();
      } else {
        document.exitFullscreen();
      }
    }
  });

  showSlide(0);
"""

# =========================================================================
# DECK 1: 01_Generative_VAEs_and_scVI_Deck.html
# =========================================================================

deck1_slides = [
    # Slide 1: Title
    """
  <!-- Slide 1: Title -->
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 50px 36px;">
      <span class="tag" style="font-size: 0.95rem; padding: 6px 16px;">Track 05 Masterclass • Part 1</span>
      <div style="margin-top: 14px;"><span class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/00_FOUNDATION_MODELS_ROADMAP.md &amp; 01</span></div>
      <h1 style="font-size: 2.45rem; color: var(--accent); margin: 18px 0 14px; font-weight: 800; letter-spacing: -0.5px;">단일 세포 딥 생성 모델 &amp; scVI 기초</h1>
      <h3 style="font-size: 1.2rem; color: #94a3b8; font-weight: 400; margin-bottom: 26px;">베이지안 심층 잠재 변수 모델, 음이항(NB/ZINB) 전사 버스팅 계층 확률 과정, 증거 하한(ELBO) 변분 추론 및 아틀라스 규모 준지도 학습(scANVI/scArches)</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 28px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08); max-width: 1020px;">
        <p style="margin-bottom: 6px;">🎯 <strong>타깃 질환 &amp; 병리 표상</strong>: 거세저항성 전립선암(mCRPC) AR/t-NEPC 표현형 전이, 기술적 배치 왜곡(10x Chromium v2/v3, Smart-seq2, CEL-seq2), 시퀀싱 깊이 편차</p>
        <p style="margin-bottom: 6px;">🧬 <strong>핵심 분석 아키텍처</strong>: <code>scVI</code> (NB/ZINB VAE, 킹마-웰링 재매개변수화, 베이지안 차별 발현) • <code>scANVI</code> (GMM 사전분포, 한계화 로그우도 준지도 라벨 전이) • <code>scArches</code> (동결된 글로벌 참조 아틀라스에 수술적 어댑터 결합)</p>
        <p style="margin-bottom: 0;">🔬 <strong>학습 중점</strong>: 선형 PCA/Harmony의 과도 보정 및 곡률 손실 극복, 물리적 카운트 노이즈와 생물학적 변이 분리, 잠재공간 z의 연속적 매니폴드 보존, scIB 벤치마킹 (iLISI, cLISI 등)</p>
      </div>
    </div>
  </div>""",

    # Slide 2: Why Linear Integration Fails
    """
  <!-- Slide 2: Why Linear Integration Fails -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/00_FOUNDATION_MODELS_ROADMAP.md § 1.1 &amp; 01_PROBABILISTIC_COUNT_MODELING.md § 1</div>
      <div class="slide-header">
        <h2>1. 왜 선형 통합(PCA/Harmony)은 한계에 부딪히는가? (과도 보정, 비선형 카운트 노이즈, 희귀 세포 손실)</h2>
        <span class="tag">Statistical Limitations</span>
      </div>
      <p>전통적인 단일세포 파이프라인(Seurat, Scanpy)은 원시 카운트에 임의의 대수 변환 <code>log(x + 1)</code>을 가한 후 선형 PCA 투영을 수행합니다. 그러나 이는 카운트 데이터의 고유한 통계물리학적 특성을 무시하여 치명적인 아티팩트를 양산합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 가우시안 가정과 log(x+1)의 결함</h3>
          <ul>
            <li>이산형 음이항 카운트 ℕ<sup>G</sup>를 연속형 정규분포로 무리하게 왜곡하여 평균-분산 비선형 종속성을 가짜 차별 발현(DEG)으로 오인</li>
            <li>저발현 유전자 영역에서 log(0+1)=0과 log(1+1)=0.693 사이의 급격한 스케일 점프로 인해 인공적인 이형성(Heteroscedasticity) 노이즈 폭발</li>
            <li>라이브러리 크기 차이를 단순 상수로 나눈 뒤 대수 변환을 취할 경우, 유전자 발현량에 따라 실제 농도비가 심각하게 왜곡됨</li>
          </ul>
        </div>
        <div class="card warning">
          <h3>2. 선형 투영(PCA)의 곡률 손실 &amp; 과도 보정</h3>
          <ul>
            <li>PCA는 평면 선형 초평면으로만 차원을 축소하므로 생체 세포의 굽은 비선형 매니폴드(Curved Manifold)를 압착하여 위상 구조를 파괴</li>
            <li>Harmony 등 휴리스틱 선형 배치 보정은 배치 간 차이를 억지로 제거하려다 실제 환자 고유의 병리적 이질성까지 함께 삭제</li>
            <li><strong>과도 보정(Over-correction)</strong>: 거세저항성 전립선암(mCRPC) 환자 특이적 약제 내성 클론이 다른 배치의 정상 상피세포로 강제 병합되는 참사 발생</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>3. 통계적 카운트 노이즈와 희귀 세포 소멸</h3>
          <ul>
            <li>scRNA-seq의 90% 이상 제로(Zero)는 기술적 포획 실패(드롭아웃)와 생물학적 미발현이 복합된 확률적 현상</li>
            <li>선형 클러스터링은 고발현 하우스키핑 유전자의 높은 분산에 압도당하여, 0.1~0.5% 미만의 희귀 신경내분비(t-NEPC) 클론을 감지하지 못함</li>
            <li>기술적 배치 편향(Batch Bias)과 순수 생물학적 발현을 직교 분리하는 심층 비선형 확률론적 생성 모델이 필수적</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 3: Probabilistic Count Modeling
    """
  <!-- Slide 3: Probabilistic Count Modeling -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/01_PROBABILISTIC_COUNT_MODELING.md § 2 &amp; § 3</div>
      <div class="slide-header">
        <h2>2. 음이항(NB/ZINB) 카운트 모델링의 수리통계학 (전사 버스팅과 과분산 수식화)</h2>
        <span class="tag">Probabilistic Physics</span>
      </div>
      <p>단일 세포 전사체는 연속형 정규분포가 아닌, 프로모터의 확률적 On/Off 전이(Transcriptional Bursting)에 의해 지배받는 이산형 확률 과정(Discrete Stochastic Process)입니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 전사 버스팅과 포아송 등분산성 붕괴</h3>
          <ul>
            <li>순수 포아송 분포는 평균과 분산이 동일함을 가정: <i>E[X] = Var(X) = &lambda;</i></li>
            <li>단일 세포 mRNA 합성은 일정한 속도로 일어나지 않고 순간 폭발적으로 전사된 후 휴지기에 돌입 (텔레그래프 모델)</li>
            <li>실제 관측 카운트의 분산은 평균보다 비약적으로 큰 극심한 <strong>과분산(Overdispersion)</strong>을 가짐:</li>
          </ul>
          <div class="math-block">Var(x<sub>ng</sub>) = &mu;<sub>ng</sub> + <span class="math-frac"><span class="num">&mu;<sub>ng</sub><sup>2</sup></span><span class="den">&theta;<sub>g</sub></span></span> &gt;&gt; &mu;<sub>ng</sub></div>
          <ul>
            <li>과분산 매개변수 &theta;<sub>g</sub>가 유전자 고유의 전사 노이즈 폭을 결정함</li>
          </ul>
        </div>
        <div class="card success">
          <h3>2. 감마-포아송 혼합과 음이항(NB) 유도</h3>
          <ul>
            <li>전사 속도 &lambda;가 감마 분포를 따를 때, 포아송 관측과의 복합 분포는 음이항(Negative Binomial) 분포로 정확히 수렴:</li>
          </ul>
          <div class="math-block">p(x<sub>ng</sub> | &mu;<sub>ng</sub>, &theta;<sub>g</sub>) = <span class="math-frac"><span class="num">&Gamma;(x<sub>ng</sub> + &theta;<sub>g</sub>)</span><span class="den">&Gamma;(&theta;<sub>g</sub>) x<sub>ng</sub>!</span></span> (<span class="math-frac"><span class="num">&theta;<sub>g</sub></span><span class="den">&theta;<sub>g</sub> + &mu;<sub>ng</sub></span></span>)<sup>&theta;<sub>g</sub></sup> (<span class="math-frac"><span class="num">&mu;<sub>ng</sub></span><span class="den">&theta;<sub>g</sub> + &mu;<sub>ng</sub></span></span>)<sup>x<sub>ng</sub></sup></div>
          <ul>
            <li>시퀀싱 깊이(l<sub>n</sub>)와 순수 유전자 발현 비율(&rho;<sub>ng</sub>)의 엄격한 곱셈적 분해:</li>
          </ul>
          <div class="math-block">&mu;<sub>ng</sub> = l<sub>n</sub> &middot; &rho;<sub>ng</sub>, &emsp; &sum;<sub>g</sub> &rho;<sub>ng</sub> = 1</div>
        </div>
        <div class="card warning">
          <h3>3. 제로 팽창(ZINB) vs 순수 NB: UMI의 진실</h3>
          <ul>
            <li>초기 Smart-seq 등 비-UMI 기술에서는 기술적 드롭아웃 확률 &pi;<sub>ng</sub>를 고려한 ZINB(Zero-Inflated NB)가 널리 사용됨:</li>
          </ul>
          <div class="math-block">p(x=0) = &pi; + (1 - &pi;) ( <span class="math-frac"><span class="num">&theta;</span><span class="den">&theta; + &mu;</span></span> )<sup>&theta;</sup></div>
          <ul>
            <li>현대 10x Chromium UMI 데이터에서는 제로의 대다수가 순수 감마-포아송(NB)의 낮은 발현 확률 및 전사 휴지기에 의한 생물학적 영으로 증명됨 (Svensson et al., 2020)</li>
            <li>scVI는 UMI 데이터셋에서 순수 NB 디코더를 기본 채택하여 과도 파라미터화 및 최적화 불안정을 원천 방지</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 4: scVI Architecture & Reparameterization
    f"""
  <!-- Slide 4: scVI Architecture & Reparameterization -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/02_SCVI_VARIATIONAL_AUTOENCODER.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>3. scVI 변분 오토인코더 아키텍처 &amp; 재매개변수화 트릭 (수학적 생성 파이프라인)</h2>
        <span class="tag">Generative Architecture</span>
      </div>
      <p><code>scVI</code>는 비선형 심층 인공신경망과 계층적 베이지안 확률 모델을 결합하여, 기술적 배치 노이즈를 직교 분리하고 세포의 순수한 생물학적 전사체 매니폴드를 저차원 잠재 공간 <code>z</code>에 투영합니다.</p>
      
      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (scVI VAE Probabilistic Architecture)</summary>
        <div class="svg-container">
{svg01}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card accent">
          <h3>1. 계층적 생성 프로세스 (Decoder)</h3>
          <ul>
            <li>저차원 생물학적 잠재 변수: <i>z<sub>n</sub> ~ N(0, I) &isin; ℝ<sup>d</sup></i> (d = 10~30)</li>
            <li>관측된 기술적 배치 범주 지시자: <i>s<sub>n</sub> &isin; {{1, ..., B}}</i></li>
            <li>신경망 디코더: <i>&rho;<sub>ng</sub> = Softmax(f<sub>&theta;</sub>(z<sub>n</sub>, s<sub>n</sub>))</i> (순수 발현 비율)</li>
            <li>관측 카운트 생성: <i>x<sub>ng</sub> ~ NB(Mean = l<sub>n</sub> &middot; &rho;<sub>ng</sub>, Dispersion = &theta;<sub>g</sub>)</i></li>
          </ul>
        </div>
        <div class="card purple">
          <h3>2. 변분 인코더 (Inference Network)</h3>
          <ul>
            <li>참 사후분포 <i>p(z|x, s)</i>의 적분이 해석적으로 불가능하므로 신경망 기반 가우시안 변분 근사 도입:</li>
          </ul>
          <div class="math-block">q<sub>&phi;</sub>(z<sub>n</sub> | x<sub>n</sub>, s<sub>n</sub>) = N(&mu;<sub>&phi;</sub>(x<sub>n</sub>, s<sub>n</sub>), diag(&sigma;<sub>&phi;</sub><sup>2</sup>(x<sub>n</sub>, s<sub>n</sub>)))</div>
          <ul>
            <li>배치를 조건화(scvi-tools 기본은 디코더에 주입)하여 잠재 공간 <i>z</i>의 배치 변동을 줄이도록 유도 (독립성은 보장되지 않음)</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 킹마-웰링 재매개변수화 트릭</h3>
          <ul>
            <li>확률 분포로부터의 직접 샘플링 노드는 역전파 기울기 전파를 차단</li>
            <li>표준정규분포 잡음 <i>&epsilon; ~ N(0, I)</i>을 주입하여 샘플링을 미분 가능한 결정론적 변환으로 전환:</li>
          </ul>
          <div class="math-block">z<sub>n</sub> = &mu;<sub>&phi;</sub>(x<sub>n</sub>, s<sub>n</sub>) + &sigma;<sub>&phi;</sub>(x<sub>n</sub>, s<sub>n</sub>) &odot; &epsilon;, &emsp; &epsilon; ~ N(0, I)</div>
          <ul>
            <li>몬테카를로 경로별(Pathwise) 기울기 계산으로 분산이 극히 낮은 엔드투엔드 역전파 학습 보장</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 5: ELBO Optimization
    f"""
  <!-- Slide 5: ELBO Optimization -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/02_SCVI_VARIATIONAL_AUTOENCODER.md § 3</div>
      <div class="slide-header">
        <h2>4. 증거 하한(ELBO) 수학적 최적화의 기하학 (재구성 손실 vs 정규화 트레이드오프)</h2>
        <span class="tag">Variational Optimization</span>
      </div>
      <p><code>scVI</code>의 훈련은 불가능한 주변 로그 우도 <code>log p(x)</code>를 직접 최대화하는 대신, 옌센 부등식(Jensen's Inequality)을 통해 유도된 수학적 하한선(ELBO)을 경사하강법으로 최대화합니다.</p>

      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (ELBO Mathematical Optimization)</summary>
        <div class="svg-container">
{svg02}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card accent">
          <h3>1. ELBO의 수리적 분해와 옌센 부등식</h3>
          <ul>
            <li>주변 로그 우도와 변분 사후분포 사이의 쿨백-라이블러(KL) 발산 관계:</li>
          </ul>
          <div class="math-block">log p(x) = E<sub>q<sub>&phi;</sub>(z|x)</sub>[log <span class="math-frac"><span class="num">p(x, z)</span><span class="den">q<sub>&phi;</sub>(z|x)</span></span>] + D<sub>KL</sub>(q<sub>&phi;</sub>(z|x) || p(z|x))</div>
          <ul>
            <li>KL 발산은 항상 0 이상이므로 첫 번째 항이 수학적 하한(ELBO)을 형성:</li>
          </ul>
          <div class="math-block">L<sub>ELBO</sub>(&phi;, &theta;) = E<sub>q<sub>&phi;</sub>(z|x)</sub>[log p<sub>&theta;</sub>(x|z)] - D<sub>KL</sub>(q<sub>&phi;</sub>(z|x) || p(z))</div>
        </div>
        <div class="card highlight">
          <h3>2. 상충하는 두 힘: 재구성 우도 vs KL 정규화</h3>
          <ul>
            <li><strong>재구성 우도 (Reconstruction)</strong>: 인코더가 압축한 잠재 벡터로부터 실제 음이항 카운트를 잘 복원하도록 유도 (과적합 위험)</li>
            <li><strong>KL 정규화 (Regularization)</strong>: 잠재 사후분포를 사전분포 <i>p(z) = N(0, I)</i>로 끌어당겨 연속적이고 매끄러운 공간 형성 유도</li>
            <li>두 손실의 균형을 통해 잠재 공간의 거리가 생물학적 분화 거리를 반영하도록 유도 (완전한 비례를 보장하지는 않음)</li>
          </ul>
        </div>
        <div class="card warning">
          <h3>3. 사후분포 붕괴(Posterior Collapse) &amp; &beta;-어닐링</h3>
          <ul>
            <li>학습 초기에 디코더가 최적화되기 전 강한 KL 손실로 인해 인코더가 입력을 무시하고 단순 <i>N(0, I)</i>로 사상하는 붕괴 현상 발생</li>
            <li>KL 항에 점진적 가중치 <i>&beta; &isin; [0, 1]</i>를 부여하는 선형 웜업(KL Warm-up) 스케줄링 적용 (전체 에포크 초기 10% 구간)</li>
            <li>초기에는 재구성에 집중하고 점진적으로 정규화를 적용하여 안정적인 수렴 및 높은 군집 분리능 보장</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 6: scANVI Semi-supervised Transfer
    """
  <!-- Slide 6: scANVI Semi-supervised Transfer -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/03_SCANVI_SEMISUPERVISED_ANNOTATION.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>5. scANVI 준지도 학습: 미지의 환자 코호트로의 정밀 레이블 전이 (한계화 손실)</h2>
        <span class="tag">Semi-Supervised Transfer</span>
      </div>
      <p>고품질로 주석화된 참조 아틀라스(Reference Atlas)의 세포 유형 레이블을, 기술적 배치가 완전히 다른 대규모 미주석 신규 환자 코호트로 고정밀 전이하는 베이지안 준지도 학습 프레임워크입니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 아틀라스 스케일의 레이블 비대칭성</h3>
          <ul>
            <li>수십만 세포의 건강한 참조 아틀라스는 수작업 검증된 고정밀 세포형 레이블 <i>c<sub>n</sub></i>을 보유</li>
            <li>반면 임상 종양 코호트(예: mCRPC 환자 침습 생검)는 수만 세포가 완전히 비주석(Unannotated) 상태</li>
            <li>기존 k-NN이나 Seurat Transfer Anchors는 배치 왜곡 하에서 희귀 악성 세포를 주변의 흔한 정상 상피로 오분류하는 한계 노출</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>2. scANVI 결합 생성 모델 (GMM Prior)</h3>
          <ul>
            <li>잠재 변수 <i>z</i>가 단일 가우시안이 아닌, 세포 유형 <i>c</i>에 조건화된 가우시안 혼합 모델(GMM)을 형성:</li>
          </ul>
          <div class="math-block">p(z | c) = N(&mu;<sub>c</sub>, diag(&sigma;<sub>c</sub><sup>2</sup>)), &emsp; c ~ Categorical(&pi;)</div>
          <ul>
            <li>인코더는 잠재 벡터 <i>z</i>를 생성하는 동시에 레이블 예측 분류기 신경망 <i>q<sub>&phi;</sub>(c|x)</i>를 공유하여 공동 최적화</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 주변화(Marginalization) 목적함수</h3>
          <ul>
            <li><strong>라벨 세포 (Supervised)</strong>: ELBO에 교차 엔트로피 분류 손실을 직접 가산:</li>
          </ul>
          <div class="math-block">L<sub>labeled</sub> = L<sub>ELBO</sub>(x, c) + &alpha; &middot; log q<sub>&phi;</sub>(c | x)</div>
          <ul>
            <li><strong>비라벨 세포 (Unsupervised)</strong>: 모든 가능한 레이블 <i>c</i>에 대해 수학적 기대값(한계화) 수행:</li>
          </ul>
          <div class="math-block">L<sub>unlabeled</sub> = &sum;<sub>c</sub> q<sub>&phi;</sub>(c | x) [ L<sub>ELBO</sub>(x, c) - log q<sub>&phi;</sub>(c | x) ]</div>
          <ul>
            <li>희귀 신경내분비(t-NEPC) 클론 식별 정확도를 F1-score 등으로 평가 (성능은 데이터 의존)</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 7: scArches Architectural Surgery
    f"""
  <!-- Slide 7: scArches Architectural Surgery -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/04_SCARCHES_ARCHITECTURAL_SURGERY.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>6. scArches 아키텍처 수술: 동결된 글로벌 아틀라스에 신규 배치 어댑터 결합 (지속 학습)</h2>
        <span class="tag">Continual Learning</span>
      </div>
      <p>새로운 환자 코호트나 신규 기술 배치가 추가될 때마다 거대한 참조 아틀라스를 처음부터 다시 학습할 필요 없이, 신경망 아키텍처에 외과적 수술을 가하여 신규 배치 어댑터만을 초고속 학습합니다.</p>

      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (scANVI &amp; scArches Atlas Mapping)</summary>
        <div class="svg-container">
{svg03}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card warning">
          <h3>1. 전역 재학습의 파멸적 한계</h3>
          <ul>
            <li>수백만 세포로 구성된 휴먼 셀 아틀라스(HCA)는 단 1회 훈련에 수천 시간의 고비용 GPU 연산 소요</li>
            <li>신규 데이터 추가 시 전체 모델을 재학습하면 기존 참조 매니폴드의 좌표계가 뒤흔들려 기존 주석과의 호환성 상실</li>
            <li><strong>파괴적 망각 (Catastrophic Forgetting)</strong>: 신규 환자 데이터에 과적합되면서 기존 참조의 미세 면역 아형 표현 능력이 퇴화</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>2. 외과적 어댑터 확장 (Architectural Surgery)</h3>
          <ul>
            <li>기존 참조 모델의 핵심 신경망 가중치 <i>W<sub>ref</sub></i>를 영구 동결(Freeze)</li>
            <li>인코더와 디코더의 첫 번째 레이어에 신규 배치 지시자 <i>s<sub>new</sub></i>만을 수용하는 전용 조건부 가중치 <i>W<sub>adapter</sub></i>를 수술적으로 증설</li>
            <li>오직 <i>W<sub>adapter</sub></i>만을 역전파로 업데이트하므로 참조 매니폴드의 좌표는 100% 보존되면서 신규 배치는 참조 공간으로 즉각 투영</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 벤치마크 성능 &amp; 아틀라스 확장성</h3>
          <ul>
            <li>신규 환자 코호트를 단일 GPU에서 어댑터만 학습하여 참조 아틀라스로 정렬 (소요 시간은 규모·하드웨어 의존)</li>
            <li>환자 원시 카운트 데이터를 병원 외부로 반출할 필요 없이 사전학습된 가중치 파일만으로 독립 실행 (HIPAA/GDPR 준수)</li>
            <li>scIB 벤치마크: 배치 통합 지표 iLISI = 0.95 달성 및 참조 레이블 보존율 cLISI = 0.99 기록</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 8: Part 1 Milestone & Part 2 Preview
    """
  <!-- Slide 8: Part 1 Milestone & Part 2 Preview -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/00_FOUNDATION_MODELS_ROADMAP.md &amp; 05_TRANSFORMER_ARCHITECTURE_SINGLE_CELL.md</div>
      <div class="slide-header">
        <h2>7. Part 1 마일스톤 총괄 &amp; Part 2 트랜스포머 파운데이션 모델 도약</h2>
        <span class="tag">Part 1 Synthesis</span>
      </div>
      <p>선형 통계 모델의 붕괴에서 시작하여 베이지안 VAE 생성 모델, 준지도 레이블 전이, 지속 학습 아틀라스 매핑까지 완성했습니다. 이제 정적 잠재공간을 넘어 유전자 조절 어텐션을 해석하는 파운데이션 모델로 도약합니다.</p>
      <div class="grid-2">
        <div class="card success">
          <h3>🏆 Part 1 핵심 성과 및 딥러닝 잠재 공간이 시사하는 암 생물학(개념 예시)</h3>
          <ul>
            <li><strong>통계물리학적 카운트 분해</strong>: 인위적 log(x+1) 변환 폐기 &rarr; 전사 버스팅 물리 기반 음이항(NB) 분해 (<i>&mu; = l &middot; &rho;</i>, 과분산 <i>&theta;</i>)</li>
            <li><strong>비선형 매니폴드 복원</strong>: <code>scVI</code> 변분 추론을 통해 기술적 배치 노이즈를 직교 분리하고, 거세저항성 전립선암의 관강 상피 &rarr; 중간 전이 &rarr; t-NEPC 연속적 가소성 매니폴드 획득</li>
            <li><strong>아틀라스 정밀 이식</strong>: <code>scANVI</code> 한계화 손실로 희귀 악성 클론 식별 가능성을 탐색, <code>scArches</code> 동결 어댑터로 신규 쿼리 세포를 참조 아틀라스에 통합 (성능 수치는 데이터 의존)</li>
            <li><strong>베이지안 차별 발현(DE)</strong>: MCMC/변분 샘플링 기반 사후분포 로그 폴드 체인지(LFC) 및 베이즈 팩터(Bayes Factor &ge; 3)를 산출하여 치료 표적 후보 우선순위화</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>🚀 Part 2 프리뷰: 대규모 트랜스포머 &amp; 가상 섭동 시뮬레이션</h3>
          <ul>
            <li><strong>세포 언어 모델링(Language of Cells)</strong>: 유전자를 단어로, 세포를 문장으로 모델링하는 Rank-Value 토큰화 및 무순서 집합 어텐션</li>
            <li><strong><code>scGPT</code> &amp; <code>Geneformer</code> 파운데이션 모델</strong>: 3,000만 단일세포 전사체 코퍼스로부터 자기지도 학습된 유전자 공발현 네트워크 및 [CLS] 세포 임베딩</li>
            <li><strong>인 실리코 유전자 섭동(In Silico Perturbation)</strong>: 컴퓨터상에서 <i>AR</i>, <i>EZH2</i>, <i>ASCL1</i>을 가상 결손(Virtual Knockout)하여 잠재 변위 벡터 <i>&Delta;z</i> 및 세포 상태 전환 궤적 시뮬레이션</li>
            <li><strong>하드웨어 가속 및 엔지니어링</strong>: FlashAttention-2, AMP FP16, 그래디언트 누적 및 100점 만점 캡스톤 프로젝트 루브릭 완결</li>
          </ul>
        </div>
      </div>
    </div>
  </div>"""
]

deck1_html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Track 05: Generative VAEs &amp; scVI Deck</title>
<style>
{css}
</style>
</head>
<body>

<div class="progress-bar" id="progressBar"></div>

<div id="presentation">
{''.join(deck1_slides)}
</div>

<!-- Navigation Controls -->
<div id="controls">
  <button class="nav-btn" id="prevBtn">◀ 이전 (Prev)</button>
  <div id="counter">Slide 1 of 8</div>
  <button class="nav-btn" id="nextBtn">다음 (Next) ▶</button>
</div>

<script>
{js_code}
</script>
</body>
</html>
"""

# =========================================================================
# DECK 2: 02_Foundation_Transformers_and_Perturbation_Deck.html
# =========================================================================

deck2_slides = [
    # Slide 1: Title
    """
  <!-- Slide 1: Title -->
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 50px 36px;">
      <span class="tag" style="font-size: 0.95rem; padding: 6px 16px;">Track 05 Masterclass • Part 2</span>
      <div style="margin-top: 14px;"><span class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/05_TRANSFORMER_ARCHITECTURE_SINGLE_CELL.md &amp; 06</span></div>
      <h1 style="font-size: 2.45rem; color: var(--accent); margin: 18px 0 14px; font-weight: 800; letter-spacing: -0.5px;">단일 세포 트랜스포머 &amp; 가상 섭동 예측</h1>
      <h3 style="font-size: 1.2rem; color: #94a3b8; font-weight: 400; margin-bottom: 26px;">단일세포 언어 모델링, Rank-Value 토큰화, scGPT/Geneformer 사전학습, 가상 유전자 섭동(In Silico Knockout) 및 전립선암 t-NEPC 역전 분화 시뮬레이션</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 28px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08); max-width: 1020px;">
        <p style="margin-bottom: 6px;">🎯 <strong>타깃 질환 &amp; 섭동 병리</strong>: 거세저항성 전립선암(mCRPC) AR 억제 가소성, <i>TP53/RB1</i> 결손 하 <i>EZH2/ASCL1</i> 드라이버 활성화, 신경내분비(t-NEPC) 전환 억제 및 가상 약물 스크리닝</p>
        <p style="margin-bottom: 6px;">🧬 <strong>핵심 트랜스포머 아키텍처</strong>: <code>scGPT</code> (GEP 마스킹 발현값 예측, GEPC 셀 임베딩 기반 값 예측, 3중 임베딩) • <code>Geneformer</code> (3,000만 세포 사전학습, 정규화 순위 인코딩, 네트워크 어텐션 중심성) • <code>FlashAttention-2</code> &amp; <code>AMP</code></p>
        <p style="margin-bottom: 0;">🔬 <strong>학습 중점</strong>: 무순서 집합(Unordered Set) 토큰화 기하학, 세포 임베딩 [CLS] 풀링, 가상 섭동 잠재 변위 벡터 &Delta;z 및 코사인 유사도, scIB 정량적 벤치마킹 및 5단계 캡스톤 프로젝트 루브릭</p>
      </div>
    </div>
  </div>""",

    # Slide 2: Single-Cell Tokenization & Attention Geometry
    f"""
  <!-- Slide 2: Single-Cell Tokenization & Attention Geometry -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/05_TRANSFORMER_ARCHITECTURE_SINGLE_CELL.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>1. 세포를 문장으로, 유전자를 단어로: Rank-Value 토큰화 및 어텐션 기하학</h2>
        <span class="tag">Tokenization &amp; Attention</span>
      </div>
      <p>단일 세포 전사체는 순서가 고정된 자연어 문장과 달리, 유전자들이 공간적 순서 없이 공존하는 무순서 집합(Unordered Set)입니다. 트랜스포머는 이를 순위 기반 토큰화와 자기주의 집중(Self-Attention)으로 모델링합니다.</p>

      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (Single-Cell Transformer Tokenization)</summary>
        <div class="svg-container">
{svg04}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card accent">
          <h3>1. 자연어와 단일 세포의 결정적 차이</h3>
          <ul>
            <li>자연어 문장이 단어의 문맥을 결정하듯, 세포 상태는 수천 개 유전자 간의 공발현 문맥(Co-expression Context)에 의해 결정됨</li>
            <li>단일세포 전사체는 유전자 간 고유한 선후 순서가 없는 <strong>무순서 집합(Unordered Set)</strong></li>
            <li>기존 절대 위치 인코딩(Absolute Positional Encoding)을 그대로 적용하면 인위적인 유전자 순서 편향이 유발되므로 위치 독립적 어텐션 설계가 필수</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>2. Rank-Value 정규화 토큰화 알고리즘</h3>
          <ul>
            <li>극심한 제로 스파시티(90%+)를 해결하기 위해 비-제로(Non-zero) 발현 유전자만을 선별</li>
            <li>세포 내 발현량을 기준으로 유전자를 내림차순 정렬한 뒤, 발현 강도를 정규화 순위(Rank)로 변환:</li>
          </ul>
          <div class="math-block">r<sub>ng</sub> = <span class="math-frac"><span class="num">Rank(x<sub>ng</sub>)</span><span class="den">|V<sub>n</sub>|</span></span> &isin; (0, 1]</div>
          <ul>
            <li>시퀀싱 깊이 변동에 완전히 불변(Invariant)하면서 유전자의 상대적 중요도를 보존하는 비모수적 토큰 시퀀스 생성</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 자기주의 집중(Self-Attention) 기하학</h3>
          <ul>
            <li>쿼리(Q), 키(K), 밸류(V) 사상 행렬을 통한 유전자 간 상호작용 에너지 계산:</li>
          </ul>
          <div class="math-block">Attention(Q, K, V) = Softmax(<span class="math-frac"><span class="num">Q K<sup>T</sup></span><span class="den">&radic;d<sub>k</sub></span></span>) V</div>
          <ul>
            <li>어텐션 가중치 행렬 <i>A<sub>ij</sub></i>는 유전자 <i>i</i>와 <i>j</i> 사이의 실시간 전사 복합체 형성 및 동적 조절 네트워크(Dynamic GRN) 강도를 직접 표현</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 3: scGPT Multitask Pretraining
    f"""
  <!-- Slide 3: scGPT Multitask Pretraining -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/06_SCGPT_PRETRAINED_CELL_EMBEDDINGS.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>2. scGPT 생성적 사전학습(GEP/GEPC) 및 파운데이션 세포 임베딩</h2>
        <span class="tag">scGPT Pretraining</span>
      </div>
      <p><code>scGPT</code>는 수천만 세포의 대규모 전사체 코퍼스로부터 3중 임베딩 합성과 다중 태스크 마스크 학습을 통해 유전자 조절 문법을 스스로 터득한 범용 파운데이션 모델입니다.</p>

      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (scGPT Multitask Pretraining &amp; Finetuning)</summary>
        <div class="svg-container">
{svg05}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card accent">
          <h3>1. 3중 임베딩 (Tri-partite) 합성</h3>
          <ul>
            <li>각 유전자 입력에 대해 세 가지 독립된 임베딩 벡터를 요소별 합산:</li>
          </ul>
          <div class="math-block">h<sub>i</sub><sup>(0)</sup> = e<sub>gene</sub>(g<sub>i</sub>) + e<sub>expr</sub>(bin(x<sub>i</sub>)) + e<sub>cond</sub>(s)</div>
          <ul>
            <li>유전자 식별자 임베딩(<i>e<sub>gene</sub></i>), 양자화된 발현량 임베딩(<i>e<sub>expr</sub></i>), 배치/질환 조건 임베딩(<i>e<sub>cond</sub></i>)의 상호작용 표상</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>2. 듀얼 자기지도 학습 목적함수</h3>
          <ul>
            <li><strong>(참고) 마스크 유전자 모델링 (MGM)</strong>: 유전자 정체성 복원은 Geneformer의 사전학습 방식이며 scGPT에서는 유전자 토큰이 주어지고 발현 값이 마스킹됨:</li>
          </ul>
          <div class="math-block">L<sub>MGM</sub> = - &sum;<sub>i &isin; M</sub> log P(g<sub>i</sub> | g<sub>&setminus;M</sub>, x<sub>&setminus;M</sub>)</div>
          <ul>
            <li><strong>유전자 발현량 예측 (GEP)</strong>: 마스킹된 위치의 발현 값을 회귀 예측 (MSE 손실); 세포 임베딩 기반 값 예측(GEPC)을 보조 손실로 사용:</li>
          </ul>
          <div class="math-block">L<sub>GEP</sub> = <span class="math-frac"><span class="num">1</span><span class="den">|M|</span></span> &sum;<sub>i &isin; M</sub> (x̂<sub>i</sub> - x<sub>i</sub>)<sup>2</sup></div>
        </div>
        <div class="card success">
          <h3>3. [CLS] 세포 임베딩과 다운스트림 전이</h3>
          <ul>
            <li>시퀀스 맨 앞에 특수 토큰 <code>[CLS]</code>를 배치하여 전체 유전자 어텐션 정보를 통합 집계</li>
            <li>최종 레이어의 [CLS] 출력 벡터 <i>h<sub>cell</sub> &isin; ℝ<sup>512</sup></i>는 세포 전체의 정보를 요약하는 표현으로 사용</li>
            <li>단일 선형 분류기를 추가하여 미세조정(Fine-tuning)하는 활용이 가능하나, 제로샷 임베딩은 더 단순한 방법에 뒤질 수 있다는 독립 벤치마크(Kedzierska 외 2025)가 있어 성능은 직접 검증해야 함</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 4: Geneformer & Network Attention
    """
  <!-- Slide 4: Geneformer & Network Attention -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/07_GENEFORMER_AND_GENOME_ATTENTION.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>3. Geneformer 3,000만 세포 사전학습 모델과 네트워크 어텐션 중심성</h2>
        <span class="tag">Geneformer &amp; Network Biology</span>
      </div>
      <p><code>Geneformer</code>는 3,000만 단일세포 전사체 코퍼스로 사전학습된 6계층 트랜스포머 모델로, 정규화 순위 인코딩과 어텐션 가중치 기반 그래프 이론을 결합하여 유전체 규모의 핵심 조절 인자를 발굴합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 대규모 사전학습 코퍼스와 아키텍처</h3>
          <ul>
            <li>다양한 인체 조직에서 수집된 3,000만(30M+) 단일세포 전사체를 학습 코퍼스로 구축</li>
            <li>6계층 트랜스포머 인코더, 254-토큰 컨텍스트 길이, 256 임베딩 차원 채택</li>
            <li>발현량이 높은 상위 254개 핵심 유전자만을 집중 분석하여 연산 복잡도를 최소화하면서 세포 정체성의 핵심 보존</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>2. 비모수적 정규화 순위 인코딩</h3>
          <ul>
            <li>절대적 발현 카운트 대신 세포 내 발현 순위(Rank)만을 입력으로 사용</li>
            <li>시퀀싱 깊이, 배치 차이, UMI 증폭 효율의 비선형 왜곡에 완전히 불변(Invariant)한 강건성 확보</li>
            <li>저발현 노이즈 유전자를 자연스럽게 꼬리 영역으로 밀어내고 고발현 드라이버 유전자의 어텐션 기여도 극대화</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 어텐션 행렬 기반 유전자 중심성 분석</h3>
          <ul>
            <li>12개 멀티헤드 어텐션의 가중치 행렬 <i>A</i>를 유전자 간 방향성 조절 네트워크 그래프의 인접 행렬로 변환</li>
            <li>어텐션 차수 중심성(Degree Centrality) 수리적 산출:</li>
          </ul>
          <div class="math-block">C<sub>attn</sub>(g<sub>i</sub>) = &sum;<sub>j &ne; i</sub> A<sub>ji</sub></div>
          <ul>
            <li>질환 특이적 세포에서 어텐션 중심성이 급상승하는 유전자를 탐색하여 알려지지 않은 마스터 전사인자(Master Regulator) 후보 탐색</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 5: In Silico Perturbation
    """
  <!-- Slide 5: In Silico Perturbation -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/08_IN_SILICO_PERTURBATION_AND_DRUG_DISCOVERY.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>4. 가상 유전자 섭동(In Silico Perturbation)의 수학적 메커니즘</h2>
        <span class="tag">In Silico Perturbation</span>
      </div>
      <p>생체 실험실(Wet-lab)의 막대한 비용과 시간 제약을 극복하고, 컴퓨터상에서 특정 유전자의 결손(KO) 또는 과발현(OE)을 가상으로 주입하여 세포의 잠재 상태 전이 궤적을 실시간 시뮬레이션합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 가상 유전자 결손의 두 가지 구현</h3>
          <ul>
            <li><strong>방식 A (토큰 삭제 / Masking)</strong>: 입력 시퀀스에서 대상 유전자 토큰 <i>g<sub>target</sub></i>을 완전히 제거하여 자기주의 어텐션 계산에서 영구 배제</li>
            <li><strong>방식 B (발현량 제로 클램핑)</strong>: 유전자 토큰은 유지하되 발현량 임베딩 값을 강제로 0으로 고정 (<i>x<sub>target</sub> = 0</i>)</li>
            <li>트랜스포머는 비선형 어텐션을 통해 해당 유전자와 상호작용하던 수천 개 하위 타깃 유전자의 발현을 연쇄적으로 재조정</li>
          </ul>
        </div>
        <div class="card warning">
          <h3>2. 잠재 공간 변위 벡터 (Perturbation Shift)</h3>
          <ul>
            <li>야생형(WT) 세포의 임베딩 <i>z<sub>wt</sub></i>와 가상 섭동 후의 임베딩 <i>z<sub>pert</sub></i> 사이의 다차원 변위 벡터 정의:</li>
          </ul>
          <div class="math-block">&Delta;z = z<sub>pert</sub> - z<sub>wt</sub> &isin; ℝ<sup>d</sup></div>
          <ul>
            <li>변위 벡터의 크기 <i>||&Delta;z||</i>는 해당 유전자가 세포의 글로벌 전사체 상태에 미치는 총 기능적 충격량(Impact Magnitude)을 정량화</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>3. 방향성 코사인 유사도와 표현형 유도력</h3>
          <ul>
            <li>목표하는 치료 반응 또는 분화 역전 방향을 나타내는 기준 벡터 <i>&Delta;z<sub>target</sub> = z<sub>luminal</sub> - z<sub>t-NEPC</sub></i> 설정</li>
            <li>섭동 변위 벡터와 목표 벡터 간의 방향성 코사인 유사도 산출:</li>
          </ul>
          <div class="math-block">cos &theta; = <span class="math-frac"><span class="num">&Delta;z &middot; &Delta;z<sub>target</sub></span><span class="den">||&Delta;z|| &middot; ||&Delta;z<sub>target</sub>||</span></span> &isin; [-1, 1]</div>
          <ul>
            <li><i>cos &theta; &gt; 0</i>이면 목표 표현형으로의 성공적 전환, <i>cos &theta; &lt; 0</i>이면 질환 악화 방향으로의 이탈을 판정</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 6: Prostate Cancer Lineage Shift Simulation
    f"""
  <!-- Slide 6: Prostate Cancer Lineage Shift Simulation -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/08_IN_SILICO_PERTURBATION_AND_DRUG_DISCOVERY.md § 3</div>
      <div class="slide-header">
        <h2>5. 전립선암 AR 억제 가소성 및 EZH2/ASCL1 가상 결손을 통한 t-NEPC 억제 시뮬레이션</h2>
        <span class="tag">Cancer Lineage Shift Simulation</span>
      </div>
      <p>거세저항성 전립선암에서 항안드로겐 제제 처방 시 발생하는 치명적 신경내분비(t-NEPC) 분화 전이를, 파운데이션 모델의 가상 섭동 시뮬레이션을 통해 분자 수준에서 차단하고 관강 표현형으로 역전시킵니다.</p>

      <details class="diagram-drawer" open>
        <summary>📐 핵심 메커니즘 인포그래픽 아키텍처 벡터 다이어그램 보기 (In Silico Perturbation &amp; Lineage Shift)</summary>
        <div class="svg-container">
{svg06}
        </div>
      </details>

      <div class="grid-3" style="margin-top: 14px;">
        <div class="card highlight">
          <h3>1. Enzalutamide 치료 병목 하 AR 결손의 비극</h3>
          <ul>
            <li>안드로겐 수용체(AR) 억제제 투여 시 관강 선암세포는 초기 억제 반응을 보이나, <i>TP53/RB1</i> 결손 배경 하에서 비가역적 계통 가소성(Lineage Plasticity) 발동</li>
            <li><i>AR</i> 가상 결손 시뮬레이션: 잠재 변위 벡터 <i>&Delta;z</i>가 신경내분비 악성화 벡터 <i>&Delta;z<sub>NEPC</sub></i>와 강한 양의 일치도(<i>cos &theta; = +0.88</i>)를 보이며 <i>SYP, CHGA, ASCL1</i> 급상승</li>
            <li>단일 표적 억제제의 한계와 대체 우회 경로 활성화의 분자적 기전 확인</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>2. EZH2 &amp; ASCL1 가상 이중 녹아웃</h3>
          <ul>
            <li>후성유전학적 리프로그래밍 효소인 <i>EZH2</i>와 신경 분화 마스터 전사인자 <i>ASCL1</i>을 가상 결손 대상으로 지정</li>
            <li><i>EZH2</i> KO 시: 염색질 억제 복합체(PRC2) 침묵 해제로 관강 분화 억제 기전 파괴 (<i>cos &theta; = -0.62</i>)</li>
            <li><i>ASCL1</i> + <i>EZH2</i> 동시 가상 결손 시: 변위 벡터가 관강(Luminal) 표현형 회복 방향으로 회전하는 가상 시나리오 (<i>cos &theta; = -0.84</i>, 설명용 가상 값)</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. 합성 치사 및 병용 요법의 임상적 번역</h3>
          <ul>
            <li>파운데이션 모델의 가상 섭동 스크리닝을 통해 <strong>[Enzalutamide (AR-i) + Tazemetostat (EZH2-i)]</strong> 병용 요법이 t-NEPC 계통 전이를 억제할 가능성을 탐색(가상 예시, 실험 검증 필요)</li>
            <li>동물 전임상 실험 전 계산으로 후보 표적 조합을 좁히는 가설 생성 용도 (기준선 대비 검증 필요)</li>
            <li>환자 맞춤형 치료 저항성 극복을 위한 계산 종양학(Computational Oncology)의 신기원</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 7: Engineering Troubleshooting & Benchmarking
    """
  <!-- Slide 7: Engineering Troubleshooting & Benchmarking -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md § 1 &amp; § 2</div>
      <div class="slide-header">
        <h2>6. 단일 세포 딥러닝 현장 트러블슈팅: FlashAttention-2, AMP &amp; scIB 벤치마킹</h2>
        <span class="tag">Engineering &amp; Benchmarking</span>
      </div>
      <p>수십만 세포 및 2,000+ 유전자 시퀀스를 처리할 때 직면하는 GPU 메모리 OOM, IO 병목, 환각 임베딩을 극복하기 위한 최신 딥러닝 시스템 엔지니어링 기법과 검증 체계입니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. FlashAttention-2: 메모리 계층 최적화</h3>
          <ul>
            <li>표준 어텐션은 <i>N &times; N</i> 어텐션 가중치 행렬을 느린 GPU HBM에 기록/로드하여 메모리 OOM 발생</li>
            <li>FlashAttention-2는 타일링(Tiling) 기법을 통해 고속 온칩 SRAM에서 Softmax를 블록 단위로 온라인 계산</li>
            <li>메모리 복잡도를 <i>O(N<sup>2</sup>)</i>에서 <i>O(N)</i>으로 획기적으로 축소, 4,096 토큰 이상에서 3~5배 연산 가속 달성</li>
          </ul>
        </div>
        <div class="card warning">
          <h3>2. FP16/BF16 자동 혼합 정밀도(AMP)</h3>
          <ul>
            <li>FP32 대비 메모리 소모를 50% 절감하고 텐서 코어(Tensor Core)를 100% 활용하는 AMP 훈련 적용</li>
            <li>언더플로우(Underflow) 방지를 위한 동적 그래디언트 스케일러(GradScaler) 탑재</li>
            <li>하드웨어 메모리 한계를 극복하기 위해 <code>Gradient Accumulation (step = 4~8)</code>을 적용하여 유효 배치 크기를 수천 세포로 확장</li>
          </ul>
        </div>
        <div class="card success">
          <h3>3. scIB (Single-cell Benchmark) 표준 검증</h3>
          <ul>
            <li>모델의 편향과 환각을 방지하기 위해 scIB 14대 표준 지표로 성능을 엄격히 계량화:</li>
            <li><strong>배치 통합성 (Batch Correction)</strong>: <code>iLISI</code>, <code>Graph Connectivity</code> (높을수록 좋음; 임계값은 데이터 의존)</li>
            <li><strong>생물학적 보존도 (Bio Conservation)</strong>: <code>cLISI</code>, <code>NMI</code>, <code>ARI</code>, <code>Isolated Label F1</code> (임계값은 데이터 의존)</li>
            <li>전체 지표의 가중 조화평균(Overall Score)을 산출하여 최적 체크포인트 결정</li>
          </ul>
        </div>
      </div>
    </div>
  </div>""",

    # Slide 8: Capstone Project Rubric & Track Milestones
    """
  <!-- Slide 8: Capstone Project Rubric & Track Milestones -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/05_generative_foundation_models/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md § 3</div>
      <div class="slide-header">
        <h2>7. 5단계 딥 생성 &amp; 파운데이션 모델 캡스톤 프로젝트 루브릭 (100점 만점)</h2>
        <span class="tag">Track 05 Capstone Rubric</span>
      </div>
      <p>5단계 딥 생성 및 파운데이션 모델 전 과정을 집대성하는 엔드투엔드 캡스톤 프로젝트의 평가 루브릭과 마스터클래스 총괄 마일스톤입니다.</p>
      <div class="grid-2">
        <div class="card accent">
          <h3>📋 5단계 엔드투엔드 캡스톤 프로젝트 루브릭 (100점 만점)</h3>
          <ul>
            <li><strong>Step 1: 통계물리학적 카운트 모델링 &amp; 전처리 (20점)</strong>: 가우시안 왜곡 분석, 전사 버스팅 음이항(NB) 파라미터 추정 및 비-제로 마스킹</li>
            <li><strong>Step 2: scVI 심층 생성 모델 구축 &amp; 변분 최적화 (20점)</strong>: 인코더-디코더 구축, 킹마-웰링 재매개변수화, ELBO 손실 수렴 및 &beta;-어닐링 검증</li>
            <li><strong>Step 3: scANVI 준지도 전이 &amp; scArches 아틀라스 수술 (20점)</strong>: 희귀 t-NEPC 클론 라벨 전이 성능을 평가 및 가중치 동결 어댑터 결합</li>
            <li><strong>Step 4: scGPT 파운데이션 토큰화 &amp; 파인튜닝 (20점)</strong>: Rank-Value 인코딩, GEP/GEPC 목적함수 학습, FlashAttention-2 가속 및 [CLS] 임베딩</li>
            <li><strong>Step 5: In Silico 유전자 섭동 &amp; scIB 벤치마킹 (20점)</strong>: <i>AR/EZH2/ASCL1</i> 가상 녹아웃 시뮬레이션, 잠재 변위 &Delta;z 분석 및 scIB 종합 0.88+ 달성</li>
          </ul>
        </div>
        <div class="card success">
          <h3>🎓 Track 05 마스터클래스 5대 핵심 엔지니어링 역량 인증</h3>
          <ul>
            <li><strong>1. 베이지안 심층 생성 모델 설계</strong>: 카운트 데이터 물리 법칙(NB/ZINB)과 변분 오토인코더(VAE)를 결합한 차세대 비선형 전사체 분석 능력</li>
            <li><strong>2. 아틀라스 규모 준지도 지속 학습</strong>: 수백만 세포 참조 모델을 파괴하지 않고 신규 임상 환자를 수분 내 통합하는 외과적 어댑터 확장 전문성</li>
            <li><strong>3. 단일세포 트랜스포머 아키텍처 구현</strong>: 무순서 집합 토큰화, 자기지도 마스크 학습, 멀티헤드 어텐션 기반 유전자 조절망 역공학 역량</li>
            <li><strong>4. 가상 유전자 섭동 예측</strong>: 컴퓨터상에서 유전자 결손 효과를 정량 시뮬레이션하여 신규 합성치사 병용 타깃을 발굴하는 예측 유전체학 실무력</li>
            <li><strong>5. 초대형 모델 가속 및 신뢰성 검증</strong>: FlashAttention-2, AMP, scIB 벤치마킹을 능숙하게 구사하는 최고 수준의 AI 생물정보학 엔지니어 완성</li>
          </ul>
        </div>
      </div>
    </div>
  </div>"""
]

deck2_html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Track 05: Foundation Transformers &amp; Perturbation Deck</title>
<style>
{css}
</style>
</head>
<body>

<div class="progress-bar" id="progressBar"></div>

<div id="presentation">
{''.join(deck2_slides)}
</div>

<!-- Navigation Controls -->
<div id="controls">
  <button class="nav-btn" id="prevBtn">◀ 이전 (Prev)</button>
  <div id="counter">Slide 1 of 8</div>
  <button class="nav-btn" id="nextBtn">다음 (Next) ▶</button>
</div>

<script>
{js_code}
</script>
</body>
</html>
"""

with open('tracks/05_generative_foundation_models/slides/01_Generative_VAEs_and_scVI_Deck.html', 'w', encoding='utf-8') as f:
    f.write(deck1_html)

with open('tracks/05_generative_foundation_models/slides/02_Foundation_Transformers_and_Perturbation_Deck.html', 'w', encoding='utf-8') as f:
    f.write(deck2_html)

print("Both decks generated successfully!")

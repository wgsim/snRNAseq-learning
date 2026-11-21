#!/usr/bin/env python3
"""
generate_tme_slides.py
Generates the two standalone interactive HTML slide decks for Track 03:
1. slides/01_TME_Communication_Foundations_and_CellChat_Deck.html
2. slides/02_NicheNet_Spatial_and_Translational_Circuits_Deck.html
"""

import os
import re

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # ENG-025: 스크립트 위치 기준(이식 가능)
TRACK_DIR = os.path.join(WORKSPACE, "tracks", "03_tme_cell_communication")
IMAGES_DIR = os.path.join(TRACK_DIR, "docs", "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")
os.makedirs(SLIDES_DIR, exist_ok=True)

def load_svg(filename: str) -> str:
    path = os.path.join(IMAGES_DIR, filename)
    with open(path, "r", encoding="utf-8") as fp:
        raw = fp.read()
    clean = re.sub(r'<\?xml[^>]*\?>', '', raw).strip()
    return clean

svg1 = load_svg("01_tme_cellular_ecosystem_and_signaling_types.svg")
svg2 = load_svg("02_liana_consensus_and_permutation_testing.svg")
svg3 = load_svg("03_cellchat_communication_probability_and_centrality.svg")
svg4 = load_svg("04_nichenet_ligand_to_target_causal_network.svg")
svg5 = load_svg("05_spatial_misty_multiview_modeling.svg")
svg6 = load_svg("06_prostate_bone_metastasis_vicious_cycle.svg")

HTML_HEAD = """<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
  :root {{
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
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }}
  body {{ background: var(--bg-primary); color: var(--text-main); overflow: hidden; height: 100vh; display: flex; flex-direction: column; }}
  
  #presentation {{ position: relative; flex: 1; width: 100%; height: 100%; overflow: hidden; }}
  .slide {{ position: absolute; top: 0; left: 0; width: 100%; height: 100%; opacity: 0; visibility: hidden; transition: opacity 0.35s ease, transform 0.35s ease; transform: scale(0.97); display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 40px 60px; }}
  .slide.active {{ opacity: 1; visibility: visible; transform: scale(1); }}
  
  .slide-content {{ max-width: 1240px; width: 100%; max-height: 86vh; overflow-y: auto; background: var(--card-bg); border: 1px solid rgba(255,255,255,0.08); border-radius: 16px; padding: 36px 44px; box-shadow: 0 20px 50px rgba(0,0,0,0.5); backdrop-filter: blur(12px); }}
  .slide-header {{ border-bottom: 2px solid var(--accent-glow); padding-bottom: 14px; margin-bottom: 20px; display: flex; justify-content: space-between; align-items: baseline; }}
  .slide-header h2 {{ font-size: 1.85rem; color: var(--accent); font-weight: 700; letter-spacing: -0.5px; }}
  .slide-header .tag {{ font-size: 0.85rem; padding: 4px 12px; border-radius: 20px; background: rgba(56, 189, 248, 0.15); color: var(--accent); font-weight: 600; text-transform: uppercase; }}
  .doc-ref {{ font-size: 0.85rem; color: #38bdf8; background: rgba(56, 189, 248, 0.12); border: 1px solid rgba(56, 189, 248, 0.3); border-radius: 6px; padding: 5px 12px; margin-bottom: 14px; display: inline-flex; align-items: center; gap: 8px; font-weight: 600; }}
  
  p, li {{ font-size: 1.02rem; line-height: 1.65; color: #cbd5e1; margin-bottom: 10px; }}
  strong {{ color: #ffffff; }}
  .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 16px; }}
  .grid-3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 18px; margin-top: 16px; }}
  .card {{ background: rgba(15, 23, 42, 0.6); border: 1px solid rgba(255,255,255,0.06); border-radius: 12px; padding: 18px; }}
  .card.accent {{ border-left: 4px solid var(--accent); }}
  .card.highlight {{ border-left: 4px solid var(--highlight); }}
  .card.warning {{ border-left: 4px solid var(--warning); }}
  .card.success {{ border-left: 4px solid var(--success); }}
  .card.purple {{ border-left: 4px solid var(--purple); }}
  .card h3 {{ font-size: 1.15rem; color: #fff; margin-bottom: 8px; }}
  
  table {{ width: 100%; border-collapse: collapse; margin-top: 14px; font-size: 0.92rem; }}
  th, td {{ padding: 10px 14px; text-align: left; border-bottom: 1px solid rgba(255,255,255,0.08); }}
  th {{ background: rgba(56, 189, 248, 0.1); color: var(--accent); font-weight: 600; }}
  
  pre {{ background: var(--code-bg); border: 1px solid rgba(255,255,255,0.1); border-radius: 8px; padding: 14px; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size: 0.88rem; color: #38bdf8; overflow-x: auto; margin-top: 10px; }}
  code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 4px; color: #38bdf8; font-size: 0.95em; }}
  
  .math-block {{ background: rgba(15, 23, 42, 0.7); border: 1px solid rgba(56, 189, 248, 0.25); border-radius: 8px; padding: 10px 16px; margin: 10px 0; text-align: center; font-family: "Cambria Math", "Latin Modern Math", "STIX Two Math", Georgia, serif; font-size: 1.12rem; color: #38bdf8; }}
  .math-frac {{ display: inline-flex; flex-direction: column; vertical-align: middle; text-align: center; padding: 0 4px; font-size: 0.95em; line-height: 1.2; }}
  .math-frac .num {{ border-bottom: 1.5px solid currentColor; padding-bottom: 2px; }}
  .math-frac .den {{ padding-top: 2px; }}

  #controls {{ height: 58px; background: var(--bg-secondary); border-top: 1px solid rgba(255,255,255,0.08); display: flex; justify-content: space-between; align-items: center; padding: 0 30px; }}
  .nav-btn {{ background: rgba(255,255,255,0.08); border: 1px solid rgba(255,255,255,0.12); color: #fff; padding: 8px 18px; border-radius: 8px; cursor: pointer; font-size: 0.9rem; transition: all 0.2s; }}
  .nav-btn:hover {{ background: var(--accent); color: var(--bg-primary); }}
  #counter {{ font-size: 0.95rem; color: var(--text-muted); font-weight: 600; }}
  .progress-bar {{ position: absolute; top: 0; left: 0; height: 3px; background: var(--accent); width: 0%; transition: width 0.3s ease; }}

  .diagram-drawer {{
    margin-top: 14px;
    background: rgba(15, 23, 42, 0.75);
    border: 1px solid rgba(56, 189, 248, 0.3);
    border-radius: 10px;
    padding: 10px 16px;
  }}
  .diagram-drawer summary {{
    cursor: pointer;
    font-weight: 600;
    color: var(--accent);
    font-size: 0.95rem;
  }}
  .diagram-drawer .svg-container svg {{
    width: 100%;
    height: auto;
    max-height: 300px;
    margin-top: 10px;
    border-radius: 8px;
  }}
</style>
</head>
<body>
<div class="progress-bar" id="progressBar"></div>
<div id="presentation">
"""

HTML_FOOT = """
</div>

<div id="controls">
  <button class="nav-btn" id="prevBtn">◀ 이전 (Prev)</button>
  <div id="counter">Slide 1 / {slide_count}</div>
  <button class="nav-btn" id="nextBtn">다음 (Next) ▶</button>
</div>

<script>
  const slides = document.querySelectorAll('.slide');
  const counter = document.getElementById('counter');
  const progressBar = document.getElementById('progressBar');
  let current = 0;

  function showSlide(index) {{
    slides.forEach((s, i) => {{
      s.classList.toggle('active', i === index);
    }});
    counter.textContent = `Slide ${{index + 1}} / ${{slides.length}}`;
    progressBar.style.width = `${{((index + 1) / slides.length) * 100}}%`;
  }}

  function nextSlide() {{
    if (current < slides.length - 1) {{
      current++;
      showSlide(current);
    }}
  }}

  function prevSlide() {{
    if (current > 0) {{
      current--;
      showSlide(current);
    }}
  }}

  document.getElementById('nextBtn').addEventListener('click', nextSlide);
  document.getElementById('prevBtn').addEventListener('click', prevSlide);

  document.addEventListener('keydown', (e) => {{
    if (e.key === ' ') e.preventDefault();
    if (e.key === 'ArrowRight' || e.key === ' ' || e.key === 'PageDown') {{
      nextSlide();
    }} else if (e.key === 'ArrowLeft' || e.key === 'PageUp') {{
      prevSlide();
    }} else if (e.key === 'f' || e.key === 'F') {{
      if (!document.fullscreenElement) {{
        document.documentElement.requestFullscreen();
      }} else {{
        document.exitFullscreen();
      }}
    }}
  }});

  showSlide(0);
</script>
</body>
</html>
"""

# ==============================================================================
# Deck 1: Foundations & CellChat Modeling
# ==============================================================================
deck1_slides = [
    # Slide 1: Title
    """
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 60px 40px;">
      <span class="tag" style="font-size: 1rem; padding: 6px 16px;">Track 03 Masterclass • Part 1</span>
      <div style="margin-top: 15px;"><span class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/00 ~ 04</span></div>
      <h1 style="font-size: 2.6rem; color: var(--accent); margin: 20px 0 16px; font-weight: 800;">종양 미세환경(TME) 세포 간 신호 전달 생태계 &amp; CellChat/LIANA 모델링</h1>
      <h3 style="font-size: 1.35rem; color: #94a3b8; font-weight: 400; margin-bottom: 28px;">세포 사회학 패러다임, 리간드-수용체 DB 합의 앙상블, 다중 소단위체 기하평균 &amp; 4대 네트워크 중심성</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 30px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08);">
        <p style="margin-bottom: 6px;">🎯 <strong>임상 타깃 모델</strong>: 거세저항성 전립선암(mCRPC) 골전이 미세환경 및 안드로겐 차단 저항성 기질 장벽</p>
        <p style="margin-bottom: 6px;">🧬 <strong>핵심 계산 알고리즘</strong>: LIANA 앙상블 RRA 순위 집계 • 1,000회 라벨 순열 검정 • CellChat 다중 소단위체 Hill 함수</p>
        <p style="margin-bottom: 0;">🔬 <strong>학습 중점</strong>: 단일세포 드롭아웃 극복 발현 임계치(10%), 비대칭 헤테로머 병목 효과, 정보(Information) &amp; flow betweenness 중심성 통신 허브</p>
      </div>
    </div>
  </div>
    """,

    # Slide 2: Paradigm Shift
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/00_TME_COMMUNICATION_ROADMAP.md § 1, 2</div>
      <div class="slide-header">
        <h2>1. 패러다임 전환: 정적 전사체 도감에서 동적 세포 사회학으로</h2>
        <span class="tag">Paradigm Shift</span>
      </div>
      <p>단일세포 RNA-seq 분석의 지향점은 <strong>"누가 그곳에 있는가(Cell Inventory)"</strong>에서 <strong>"누가 누구에게 어떤 생화학적 언어로 지시하는가(Cell-Cell Cross-talk)"</strong>로 진화했습니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 세포 자율적(Cell-Autonomous) 한계</h3>
          <ul>
            <li>암세포가 <i>MET</i>, <i>CXCR4</i> 수용체를 발현해도 동족 리간드(HGF, CXCL12) 부재 시 경로 비활성</li>
            <li>정적 마커 DEG 나열만으로는 종양-기질-면역 동적 상호작용 추적 불가</li>
            <li>미세환경 수용성(Microenvironmental Receptivity) 간과</li>
          </ul>
        </div>
        <div class="card warning">
          <h3>2. 치료 저항성의 외인성(Extrinsic) 기전</h3>
          <ul>
            <li>ADT(Enzalutamide) 저항성은 유전 변이 외에 CAF 분비 <i>WNT16B</i>, <i>SPP1</i>에 의해 강력 유도</li>
            <li>골수 미세환경 간질세포(BM-MSCs)의 paracrine 생존 신호</li>
            <li>면역 배제 기질 피막(myCAF Fibronectin/Collagen) 장벽</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>3. 희소 면역관문 Dropout 극복</h3>
          <ul>
            <li><i>PDCD1</i>, <i>CTLA4</i> 등 면역 수용체의 극심한 전사체 드롭아웃</li>
            <li>단일 유전자 발현 비교 대신 복합 리간드-수용체 상호작용 확률 및 네트워크 기반 추론 필수</li>
            <li>위양성/위음성 제어를 위한 다중 DB 합의 집계</li>
          </ul>
        </div>
      </div>
    </div>
  </div>
    """,

    # Slide 3: Signaling Modalities & SVG 01
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/01_TME_BIOLOGY_AND_LIGAND_RECEPTORS.md § 1, 2</div>
      <div class="slide-header">
        <h2>2. 종양 미세환경의 3대 신호 전달 양식과 세포 사회학</h2>
        <span class="tag">Signaling Modalities</span>
      </div>
      <p>종양 세포와 TME 기질·면역세포는 물리적 거리 및 분자 결합 형태에 따라 3가지 주요 신호 양식으로 상호작용합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 병렬접촉 (Juxtacrine)</h3>
          <p><strong>유효 거리: &lt; 0.1 µm (직접 접촉)</strong></p>
          <ul>
            <li>세포막 표면 단백질 간 직접 결합: <i>Notch-Delta</i>, <i>Ephrin-Eph</i>, <i>PD-L1/PD-1</i></li>
            <li>간극결합(Gap Junction) 코넥신 채널을 통한 2차 메신저(cAMP, Ca<sup>2+</sup>) 수송</li>
          </ul>
        </div>
        <div class="card success">
          <h3>2. 근분비 (Paracrine)</h3>
          <p><strong>유효 거리: 1 ~ 200 µm (국소 확산)</strong></p>
          <ul>
            <li>수용성 사이토카인/케모카인 분비: <i>CXCL12-CXCR4</i>, <i>TGFB1-TGFBR1/2</i></li>
            <li>세포외기질(ECM) 상호작용: <i>SPP1-CD44</i>, 인테그린-콜라겐 결합</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>3. 자가분비 (Autocrine)</h3>
          <p><strong>유효 거리: 자가 세포막 표면 (루프 피드백)</strong></p>
          <ul>
            <li>세포 자신이 분비한 리간드가 자신의 수용체에 결합: <i>IL6-IL6R</i>, <i>VEGFA-VEGFR2</i></li>
            <li>자가 증폭 루프를 통한 줄기능 유지 및 화학요법 내성 획득</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 01] TME 세포 생태계와 3대 신호 전달 양식 벡터 구조 (Click to toggle)</summary>
        <div class="svg-container">
          {svg1}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 4: Database Curation & LIANA Consensus (SVG 02)
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/02_DATABASE_CURATION_AND_CONSENSUS.md § 2, 3</div>
      <div class="slide-header">
        <h2>3. 리간드-수용체 DB 이질성 극복 &amp; LIANA 합의 순위 집계</h2>
        <span class="tag">Database Consensus</span>
      </div>
      <p>단일 DB(CellChat, CellPhoneDB 등)에만 의존할 경우 데이터베이스 간 L-R 쌍 중첩률이 28% 미만에 불과해 심각한 편향이 발생합니다.</p>
      <div class="grid-2">
        <div class="card accent">
          <h3>주요 L-R 데이터베이스 한계 및 편향</h3>
          <ul>
            <li><strong>CellPhoneDB v4.0</strong>: 다중 소단위체 복합체 강점이나 구조 단백질 누락</li>
            <li><strong>CellChatDB v2.0</strong>: 분비성 리간드 및 보조인자 충실하나 세포외기질 커버리지 편차</li>
            <li><strong>SingleCellSignalR / ICELLNET</strong>: 점수 체계 상이, 문헌 검증 기준 불일치</li>
            <li><strong>해결책</strong>: OmniPath 기반 <strong>LIANA 앙상블 프레임워크</strong> 도입</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>합의 순위 집계 (Consensus Rank; LIANA는 RRA)</h3>
          <div class="math-block">
            Rank<sub>평균</sub>(i, j) = <span class="math-frac"><span class="num">1</span><span class="den">M</span></span> ∑<sub>m=1</sub><sup>M</sup> rank<sub>m</sub>(i, j)
          </div>
          <ul>
            <li>5대 CCC 알고리즘(NatMI, SingleCellSignalR, Connectome, CellPhoneDB, LogFC) 동시 평가</li>
            <li>최소 3개 이상 알고리즘 합의 통과 상호작용만 엄밀 선별</li>
            <li>위 평균 순위(Borda식)는 일반 순위 집계 예시이며, LIANA는 Robust Rank Aggregation(RRA; 옵션 mean)으로 집계 (Borda는 미구현)</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 02] LIANA 다중 알고리즘 앙상블 및 순열 검정 워크플로우 (Click to toggle)</summary>
        <div class="svg-container">
          {svg2}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 5: Expression Thresholding & Permutations
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>4. 드롭아웃 극복 발현 임계치 &amp; 질량작용 법칙 기반 순열 검정</h2>
        <span class="tag">Statistical Rigor</span>
      </div>
      <p>평균 발현량(Mean Expression)의 함정을 피하기 위해 <strong>비영(Non-zero) 발현 비율 임계치(10%)</strong>와 <strong>세포 라벨 순열 검정</strong>을 결합합니다.</p>
      <div class="grid-2">
        <div class="card warning">
          <h3>발현 비율 임계치 (Expr Proportion &ge; 10%)</h3>
          <ul>
            <li>단 1개의 초고발현 이상치(Outlier) 세포가 클러스터 전체 평균을 왜곡하는 현상 차단</li>
            <li>리간드 분비 세포와 수용체 표출 세포 모두 클러스터 내 10% 이상 발현 필수</li>
            <li>극소수 세포 분비 리간드에 의한 거짓 양성(False Positive) 통신 원천 제거</li>
          </ul>
        </div>
        <div class="card success">
          <h3>세포 라벨 1,000회 몬테카를로 순열 검정</h3>
          <div class="math-block">
            p<sub>emp</sub> = <span class="math-frac"><span class="num">∑<sub>b=1</sub><sup>B</sup> 𝕀(S<sub>perm</sub><sup>(b)</sup> &ge; S<sub>obs</sub>) + 1</span><span class="den">B + 1</span></span>
          </div>
          <ul>
            <li>세포 유형 라벨을 무작위 셔플링하여 귀무가설(Null Hypothesis) 분포 생성</li>
            <li>Benjamini-Hochberg FDR 보정 (q &lt; 0.05) 엄격 적용</li>
          </ul>
        </div>
      </div>
      <div class="card accent" style="margin-top: 14px;">
        <h3>질량작용의 법칙(Law of Mass Action) 수학적 정식화</h3>
        <p>신호 세기 S = [L] &times; [R] 은 리간드와 수용체의 농도 곱에 비례하며, 해리 상수 K<sub>d</sub>와 힐 계수 n에 따라 수용체 결합 포화 반응을 정량화합니다.</p>
      </div>
    </div>
  </div>
    """,

    # Slide 6: CellChat Multi-subunit & Centrality (SVG 03)
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/04_CELLCHAT_MULTISUBUNIT_MODELING.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>5. CellChat 다중 소단위체 기하평균 모델 &amp; 4대 네트워크 중심성</h2>
        <span class="tag">Multi-subunit Architecture</span>
      </div>
      <p>단백질 수용체 복합체는 여러 소단위체로 구성되며, 단 하나의 소단위체라도 결핍되면 신호가 차단되는 <strong>병목 원리</strong>를 기하평균으로 모델링합니다.</p>
      <div class="grid-2">
        <div class="card purple">
          <h3>다중 소단위체 기하평균 (Geometric Mean)</h3>
          <div class="math-block">
            R<sub>eff</sub> = ( ∏<sub>k=1</sub><sup>K</sup> R<sub>k</sub> )<sup>1/K</sup>
          </div>
          <ul>
            <li>산술평균의 치명적 결함(0 + 10 = 5) 방지: 한 소단위체가 0이면 R<sub>eff</sub> = 0</li>
            <li>IL-6 수용체: IL6R + IL6ST(gp130) 이종 복합체 모델링</li>
            <li>TGF-β 수용체: TGFBR1 + TGFBR2 이종 이량체 모델링</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>4대 네트워크 중심성 (Graph Centrality)</h3>
          <ul>
            <li><strong>송신 강도 (Out-degree)</strong>: 총 분비 리간드 플럭스 (주요 발신 세포)</li>
            <li><strong>수신 강도 (In-degree)</strong>: 총 수신 리간드 민감도 (표적 반응 세포)</li>
            <li><strong>매개 중심성 (Flow betweenness)</strong>: 다세포 신호 중계 허브 세포</li>
            <li><strong>정보 중심성 (Information)</strong>: 네트워크 전역 영향력 중심 세포</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 03] CellChat 통신 확률 언덕 함수 및 2D 네트워크 중심성 공간 (Click to toggle)</summary>
        <div class="svg-container">
          {svg3}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 7: Interactive Lab & Pipeline 01
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 실습: tracks/03_tme_cell_communication/pipelines/01_liana_and_cellchat_consensus_pipeline.py</div>
      <div class="slide-header">
        <h2>6. 대화형 코딩 실습: Pipeline 01 실행 및 핵심 결과</h2>
        <span class="tag">Hands-on Lab</span>
      </div>
      <p>1,500개 mCRPC 단일세포 전사체 시뮬레이션 데이터를 바탕으로 LIANA 합의 및 CellChat 중심성을 실시간 연산합니다.</p>
      <pre><code># Pipeline 01 실행 커맨드
python tracks/03_tme_cell_communication/pipelines/01_liana_and_cellchat_consensus_pipeline.py

# 주요 실행 단계 및 산출물:
[1/5] Simulating 1,500 mCRPC single cells across 8 cell types... (metadata.csv)
[2/5] Curating multi-subunit L-R database for prostate bone metastasis... (13 curated complexes)
[3/5] Computing trimmed geometric means & Hill probabilities (Kh=0.5)... (832 interactions)
[4/5] Executing 1,000 cell-label permutations for empirical FDR & Borda consensus... (39 significant, FDR &lt; 0.05)
[5/5] Calculating 2D network centrality (Out/In-degree, Betweenness, PageRank)... (centrality_metrics.csv)</code></pre>
      <div class="grid-2" style="margin-top: 14px;">
        <div class="card success">
          <h3>핵심 발견 1: myCAF의 압도적 송신 중심성</h3>
          <p>myCAF는 Out-degree = 0.842, PageRank = 0.281로 종양 미세환경 전체에서 가장 강력한 리간드(<i>COL1A1</i>, <i>FN1</i>, <i>TGFB1</i>, <i>CXCL12</i>) 송신 허브로 군림합니다.</p>
        </div>
        <div class="card highlight">
          <h3>핵심 발견 2: CD8+ T세포와 파골 전구세포의 표적 수신</h3>
          <p>Exhausted CD8+ T세포는 <i>PDCD1</i>, <i>CTLA4</i>, <i>TGFBR</i>을 집중 수신하며 면역 억제되고, Pre-Osteoclast는 <i>RANK</i> 수신을 통해 골 흡수 분화 신호를 수신합니다.</p>
        </div>
      </div>
    </div>
  </div>
    """,

    # Slide 8: Part 1 Summary & Part 2 Preview
    """
  <div class="slide">
    <div class="slide-content" style="text-align: center; padding: 50px 40px;">
      <div class="doc-ref">📖 Track 03 Masterclass • Part 1 요약 및 Part 2 로드맵</div>
      <div class="slide-header" style="justify-content: center;">
        <h2>Part 1 마스터 요약 &amp; Part 2 예고</h2>
      </div>
      <div style="text-align: left; max-width: 960px; margin: 0 auto;">
        <div class="card accent" style="margin-bottom: 18px;">
          <h3>Part 1 핵심 달성 마일스톤</h3>
          <p>✔ 종양 미세환경 세포 사회학 패러다임 정립 및 외인성 치료 저항성 기전 이해</p>
          <p>✔ 3대 신호 양식(Juxtacrine/Paracrine/Autocrine) 및 물리화학적 유효 거리 정량화</p>
          <p>✔ OmniPath/LIANA 앙상블 합의 순위 집계(RRA)로 DB 편향 완화</p>
          <p>✔ 10% 발현 임계치 &amp; 1,000회 몬테카를로 세포 라벨 순열 검정으로 FDR &lt; 0.05 엄격 검증</p>
          <p>✔ CellChat 다중 소단위체 기하평균 모델 및 4대 그래프 중심성(Out/In, Flow betweenness, Information) 학습</p>
        </div>
        <div class="card highlight">
          <h3>Part 2 예고: NicheNet 표적 인과 추론, 공간 전사체 &amp; 골 전이 악순환 회로</h3>
          <p>상류 리간드 ──► RWR 세포 내 신호망 ──► 하류 표적 유전자 조절망(eGRN) 인과 추론 ──► LIANA+ 공간 거리 감쇄 &amp; MISTY 다중 뷰 기계학습 ──► 전립선암 뼈 전이 치명적 악순환 회로 및 임상 타깃 도출</p>
        </div>
      </div>
    </div>
  </div>
    """
]

# Write Deck 1
deck1_content = HTML_HEAD.format(title="Track 03 Masterclass Part 1: TME Foundations & CellChat Modeling") + "".join(deck1_slides) + HTML_FOOT.format(slide_count=len(deck1_slides))
deck1_path = os.path.join(SLIDES_DIR, "01_TME_Communication_Foundations_and_CellChat_Deck.html")
with open(deck1_path, "w", encoding="utf-8") as fp:
    fp.write(deck1_content)
print(f"✅ Generated Slide Deck 1: {deck1_path} ({len(deck1_content):,} chars)")


# ==============================================================================
# Deck 2: NicheNet, Spatial CCC & Translational Circuits
# ==============================================================================
deck2_slides = [
    # Slide 1: Title
    """
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 60px 40px;">
      <span class="tag" style="font-size: 1rem; padding: 6px 16px;">Track 03 Masterclass • Part 2</span>
      <div style="margin-top: 15px;"><span class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/05 ~ 09</span></div>
      <h1 style="font-size: 2.6rem; color: var(--accent); margin: 20px 0 16px; font-weight: 800;">NicheNet 세포 내 신호 인과 추론, 공간 전사체 통신 &amp; 골 전이 표적 회로</h1>
      <h3 style="font-size: 1.35rem; color: #94a3b8; font-weight: 400; margin-bottom: 28px;">무작위 보행(RWR) 조절 잠재력, LIANA+ 공간 거리 감쇄, MISTY 다중 뷰 &amp; 치명적 악순환 5자 회로</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 30px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08);">
        <p style="margin-bottom: 6px;">🎯 <strong>임상 모델</strong>: mCRPC 골수 침윤, 골아세포-파골세포 골 흡수/신생 불균형, SPP1+ TAM 면역 억제 니치</p>
        <p style="margin-bottom: 6px;">🧬 <strong>핵심 계산 알고리즘</strong>: NicheNet PageRank RWR • AUROC/AUPR 리간드 활성도 • LIANA+ 가우시안 거리 커널 • MISTY 랜덤 포레스트</p>
        <p style="margin-bottom: 0;">🔬 <strong>학습 중점</strong>: 상류 신호와 하류 전사체 인과 연계, 세포 내/인접/근거리 미세환경 기여도 분해, 이중/삼중 억제제 시너지</p>
      </div>
    </div>
  </div>
    """,

    # Slide 2: NicheNet Causal Attribution (SVG 04)
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/05_NICHENET_INTRACELLULAR_SIGNALING.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>1. NicheNet 세포 내 신호 인과 추론: RWR 조절 잠재력 &amp; 리간드 활성도</h2>
        <span class="tag">Causal Inference</span>
      </div>
      <p>단순 L-R 결합 가능성을 넘어, 수용체 결합이 수신 세포 내부의 <strong>표적 유전자 발현(DEG)을 실제로 유도했는지</strong> 인과적으로 증명합니다.</p>
      <div class="grid-2">
        <div class="card accent">
          <h3>사전 지식 통합 그래프 무작위 보행 (RWR)</h3>
          <div class="math-block">
            p<sub>t+1</sub> = (1 - α) A<sub>norm</sub> p<sub>t</sub> + α p<sub>0</sub>
          </div>
          <ul>
            <li>신호 전달 네트워크(Signaling) + 유전자 조절망(GRN) 사전 지식 통합</li>
            <li>리간드 노드에서 재시작 확률 α=0.15로 무작위 보행 수행하여 정상 상태 확률 수렴</li>
            <li>조절 잠재력 행렬 W(Ligand &times; Target) 사전 구축</li>
          </ul>
        </div>
        <div class="card success">
          <h3>리간드 활성도 예측 (AUROC &amp; Pearson r)</h3>
          <ul>
            <li>수신 세포의 차별 발현 유전자(DEG) 시그니처와 조절 잠재력 W 벡터 간 상관성 분석</li>
            <li><strong>AUROC &gt; 0.85</strong>: 수신 세포 DEG 유도에 가장 기여도가 높은 핵심 리간드 랭킹</li>
            <li>mCRPC T세포 고갈: <i>TGFB1</i>, <i>LGALS9</i> 상위 랭크 검증</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 04] NicheNet 리간드-표적 인과 신호망 및 RWR 조절 잠재력 아키텍처 (Click to toggle)</summary>
        <div class="svg-container">
          {svg4}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 3: Spatial CCC & LIANA+
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/06_SPATIAL_CELL_COMMUNICATION_LIANA_PLUS.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>2. 공간 전사체 기반 세포 간 통신: LIANA+ 가우시안 거리 감쇄 커널</h2>
        <span class="tag">Spatial CCC</span>
      </div>
      <p>단일세포 해리 데이터에서 발생하는 물리적 거리 왜곡을 극복하기 위해 <strong>조직 공간 좌표(x, y)</strong>와 <strong>가우시안 확산 감쇄 커널</strong>을 통합합니다.</p>
      <div class="grid-2">
        <div class="card highlight">
          <h3>가우시안 확산 거리 감쇄 커널 (Distance-Decay Kernel)</h3>
          <div class="math-block">
            w<sub>ij</sub> = exp( - <span class="math-frac"><span class="num">d(i, j)<sup>2</sup></span><span class="den">2 σ<sup>2</sup></span></span> )
          </div>
          <ul>
            <li>d(i, j): 발신 세포 i와 수신 세포 j 사이의 유클리드 공간 거리 (µm)</li>
            <li>σ: 리간드 확산 길이 상수 (Paracrine: 30~50 µm, Juxtacrine: 5~10 µm)</li>
            <li>거리가 멀어질수록 신호 가중치가 지수함수적으로 0으로 수렴</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>공간 제약 국소 유효 플럭스 (Local Effective Flux)</h3>
          <div class="math-block">
            F<sub>local</sub>(i &rarr; j) = L<sub>i</sub> &times; R<sub>j</sub> &times; w<sub>ij</sub>
          </div>
          <ul>
            <li>동일 공간 니치(Cellular Niche) 내 물리적으로 인접한 세포 쌍의 실제 접촉 통신만 정밀 계산</li>
            <li>수백 µm 떨어진 세포 간의 허구적 L-R 통신(False-positive cross-talk) 90% 이상 제거</li>
          </ul>
        </div>
      </div>
      <div class="card warning" style="margin-top: 14px;">
        <h3>10x Visium HD (2 µm) 및 Xenium (&lt;0.1 µm) 인시츄 데이터 적용</h3>
        <p>세포 분할 경계(Cell Segmentation Polygon) 및 서브셀룰러 전사체 좌표를 바탕으로 진정한 단일 세포 간 직접 접촉 신호 정량화 실현.</p>
      </div>
    </div>
  </div>
    """,

    # Slide 4: Multiview Modeling MISTy (SVG 05)
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/07_MULTIVIEW_MODELING_MISTY.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>3. MISTY 다중 뷰 기계학습: Intra / Juxta / Para 기여도 분해</h2>
        <span class="tag">MISTy Multiview</span>
      </div>
      <p>표적 마커의 발현 변량을 세포 내부(Intra), 직접 접촉(Juxta), 주변 근분비(Para) 미세환경 뷰로 나누어 기계학습 앙상블로 모델링합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>View 1: Intraview (세포 자율적)</h3>
          <ul>
            <li>세포 자체 내부의 다른 유전자/단백질 발현</li>
            <li>고유 전사 인자 및 세포 주기 상태</li>
            <li>내부 기저 발현 변량 포획</li>
          </ul>
        </div>
        <div class="card success">
          <h3>View 2: Juxtaview (직접 접촉)</h3>
          <ul>
            <li>반경 r &le; 15 µm 이내 인접 세포 마커</li>
            <li>물리적 막 접촉(Notch, Cadherin, Integrin)</li>
            <li>즉각적 접촉 신호 영향력 정량화</li>
          </ul>
        </div>
        <div class="card purple">
          <h3>View 3: Paraview (광역 미세환경)</h3>
          <ul>
            <li>반경 r &le; 100 µm 범위 가우시안 평활화</li>
            <li>사이토카인 확산 구배 및 저산소 미세환경</li>
            <li>기질 피막 및 면역 침윤 구배 효과</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 05] MISTy 다중 뷰 앙상블 아키텍처 및 변량 기여도 분해 (Click to toggle)</summary>
        <div class="svg-container">
          {svg5}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 5: Bone Metastasis Vicious Cycle (SVG 06)
    f"""
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>4. 전립선암 거세저항성 뼈 전이(mCRPC) 치명적 악순환 5자 회로</h2>
        <span class="tag">Vicious Cycle Circuit</span>
      </div>
      <p>mCRPC 골전이 병변은 5대 핵심 세포 집단이 상호 되먹임(Feedback)을 형성하며 자가 증폭 악순환을 지속합니다.</p>
      <div class="grid-2">
        <div class="card highlight">
          <h3>악순환 루프(Vicious Cycle) 4단계 동역학</h3>
          <ol>
            <li><strong>호밍 &amp; 안착</strong>: 골수 내피세포 <i>CXCL12</i> ➔ 암세포 <i>CXCR4</i> 축으로 골내막 착상</li>
            <li><strong>골 흡수 촉진</strong>: 암세포가 <i>PTHrP</i>, <i>IL6</i> 분비 ➔ 조골세포 <i>RANKL</i> 상향 ➔ 파골세포 분화 촉발</li>
            <li><strong>기질 인자 방출</strong>: 뼈 기질 용해로 유리된 활성형 <i>TGFB1</i>, <i>IGF1</i>이 암세포 재자극</li>
            <li><strong>면역 억제 피막</strong>: <i>SPP1+</i> TAM 및 myCAF가 <i>PD-L1</i>, <i>TGF-β</i> 장벽으로 CD8+ T세포 무력화</li>
          </ol>
        </div>
        <div class="card accent">
          <h3>5자 신호 상호작용 참여 세포군</h3>
          <ul>
            <li><strong>mCRPC 암세포</strong>: <i>PTHrP</i>, <i>TGFB1</i>, <i>CXCR4</i>, <i>AR</i> 과발현</li>
            <li><strong>조골세포 (Osteoblast)</strong>: <i>RANKL</i>, <i>COL1A1</i>, <i>CXCL12</i> 발현</li>
            <li><strong>파골세포 (Osteoclast)</strong>: <i>RANK</i>, <i>CTSK</i>, <i>TRAP</i> 매개 골 용해</li>
            <li><strong>골수 myCAF</strong>: <i>WNT16B</i>, <i>IL6</i>, 치밀한 ECM 장벽 형성</li>
            <li><strong>SPP1+ TAM / Treg</strong>: <i>SPP1</i>, <i>IL10</i> 매개 T세포 고갈 유도</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📐 [다이어그램 06] 전립선암 뼈 전이 니치 5자 악순환 회로 및 표적 치료 표적 (Click to toggle)</summary>
        <div class="svg-container">
          {svg6}
        </div>
      </details>
    </div>
  </div>
    """,

    # Slide 6: Targeted Therapy & In Silico Screening
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 실습: tracks/03_tme_cell_communication/pipelines/02_nichenet_and_spatial_communication_pipeline.py</div>
      <div class="slide-header">
        <h2>5. 이중/삼중 표적 치료 전략 &amp; 컴퓨터 가상 스크리닝 (Pipeline 02)</h2>
        <span class="tag">Translational Discovery</span>
      </div>
      <p>Pipeline 02의 악순환 신호 플럭스 정량 모델을 통해 단일 표적 억제 대비 <strong>복합 병용 요법의 시너지 억제 효과</strong>를 입증했습니다.</p>
      <table>
        <thead>
          <tr>
            <th>치료 전략 (Therapeutic Strategy)</th>
            <th>표적 경로 및 분자 기전</th>
            <th>악순환 플럭스 억제율 (%)</th>
            <th>임상적 의의 및 저항성 극복</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>단일 요법 1 (Denosumab)</strong></td>
            <td>RANKL 중화 항체 (골 흡수 차단)</td>
            <td>41.2%</td>
            <td>골 흡수 억제되나 TGF-β/myCAF 매개 암 증식 지속</td>
          </tr>
          <tr>
            <td><strong>단일 요법 2 (Galunisertib)</strong></td>
            <td>TGF-β 수용체 키나제 억제제</td>
            <td>36.8%</td>
            <td>EMT 억제되나 RANKL-골용해 루프 잔존</td>
          </tr>
          <tr>
            <td><strong>이중 병용 (Denosumab + Galunisertib)</strong></td>
            <td>RANKL + TGF-β 동시 차단</td>
            <td>78.5%</td>
            <td>골흡수와 기질 인자 재방출 동시 차단 (강력한 시너지)</td>
          </tr>
          <tr style="background: rgba(56, 189, 248, 0.15); font-weight: bold;">
            <td><strong>삼중 복합 (Denosumab + Galunisertib + Enzalutamide)</strong></td>
            <td>RANKL + TGF-β + AR 축 동시 타깃</td>
            <td>94.7%</td>
            <td>골전이 악순환 루프 완전 와해 및 거세저항성 암세포 고갈</td>
          </tr>
        </tbody>
      </table>
      <div class="card success" style="margin-top: 14px;">
        <h3>인실리코(In Silico) 약물 섭동 시뮬레이션의 혁신성</h3>
        <p>단일세포 네트워크 플럭스 행렬을 통해 환자 맞춤형 미세환경 최적 병용 치료법을 사전에 예측하여 불필요한 약물 독성 방지 및 임상 시험 성공률 극대화.</p>
      </div>
    </div>
  </div>
    """,

    # Slide 7: Troubleshooting & Capstone Rubric
    """
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/03_tme_cell_communication/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md § 1 ~ 3</div>
      <div class="slide-header">
        <h2>6. 실전 트러블슈팅, 벤치마킹 품질 제어 &amp; 캡스톤 루브릭</h2>
        <span class="tag">QC &amp; Capstone</span>
      </div>
      <p>세포 간 통신 분석 시 마주치는 3대 기술적 결함과 해결책, 그리고 캡스톤 프로젝트 평가 루브릭을 준수합니다.</p>
      <div class="grid-3">
        <div class="card warning">
          <h3>1. 세포 수 편향 (Cell-size bias)</h3>
          <ul>
            <li>거대 클러스터(예: 종양세포 5,000개 vs 희소 수지상세포 50개) 간 순열 검정 왜곡</li>
            <li>클러스터별 서브샘플링(100~500 cells) 적용 필수</li>
          </ul>
        </div>
        <div class="card highlight">
          <h3>2. 이종 다중 복합체 누락</h3>
          <ul>
            <li>α 및 β 소단위체 중 하나라도 누락되면 거짓 양성 발생</li>
            <li>CellChatDB 다중 소단위체 기하평균 및 OmniPath 큐레이션 검증</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>3. 단일세포 전사체 해리 왜곡</h3>
          <ul>
            <li>원거리 세포 간 허구적 통신을 공간 좌표 기반 가우시안 감쇄 필터링으로 제거</li>
            <li>Spatial CCC / LIANA+ 상호 검증 필수</li>
          </ul>
        </div>
      </div>
      <div class="card purple" style="margin-top: 14px;">
        <h3>🎓 Track 03 캡스톤 프로젝트 루브릭 (100점 만점)</h3>
        <p>• <strong>데이터 전처리 &amp; 임계치 검증 (25점)</strong>: 10% 발현 임계치 및 드롭아웃 보정 타당성</p>
        <p>• <strong>LIANA 앙상블 &amp; CellChat 모델링 (25점)</strong>: 다중 소단위체 기하평균 및 RRA 순위 집계 정확성</p>
        <p>• <strong>NicheNet &amp; 공간 플럭스 분석 (25점)</strong>: RWR 기반 표적 인과 추론 및 가우시안 거리 감쇄 필터링</p>
        <p>• <strong>임상 기전 해석 &amp; 치료 타깃 도출 (25점)</strong>: 전립선암 뼈 전이 악순환 회로 차단 전략 및 시너지 정량화</p>
      </div>
    </div>
  </div>
    """,

    # Slide 8: Full Track Conclusion & Curriculum Integration
    """
  <div class="slide">
    <div class="slide-content" style="text-align: center; padding: 50px 40px;">
      <div class="doc-ref">📖 Track 03: TME Cell-Cell Communication 완주 요약</div>
      <div class="slide-header" style="justify-content: center;">
        <h2>Track 03 전과정 완주 &amp; 마스터 커리큘럼 로드맵 통합</h2>
      </div>
      <div style="text-align: left; max-width: 960px; margin: 0 auto;">
        <div class="card accent" style="margin-bottom: 16px;">
          <h3>전체 10개 챕터 마일스톤 완비</h3>
          <p>✔ <strong>Ch 00 ~ 01</strong>: TME 세포 사회학 패러다임 및 3대 신호 양식(Juxtacrine/Paracrine/Autocrine) 물리화학</p>
          <p>✔ <strong>Ch 02 ~ 04</strong>: DB 큐레이션, LIANA RRA 합의 집계, 발현 임계치, CellChat 다중 소단위체 기하평균 &amp; 중심성</p>
          <p>✔ <strong>Ch 05 ~ 07</strong>: NicheNet RWR 표적 인과 추론, LIANA+ 공간 거리 감쇄, MISTY 다중 뷰 앙상블 머신러닝</p>
          <p>✔ <strong>Ch 08 ~ 09</strong>: mCRPC 뼈 전이 악순환 5자 회로, 삼중 표적 치료 시뮬레이션, 트러블슈팅 및 캡스톤 루브릭</p>
        </div>
        <div class="card success">
          <h3>전체 커리큘럼과의 시너지 연계</h3>
          <p>• <strong>Track 00 &amp; 01 (단일세포/멀티옴 eGRN)</strong> ➔ 세포 유형별 하류 전사인자 조절망 제공</p>
          <p>• <strong>Track 02 (공간 전사체학)</strong> ➔ 2D/3D 조직 내 실제 세포 간 유클리드 좌표 및 공간 니치 제공</p>
          <p>• <strong>Track 03 (세포 간 통신망)</strong> ➔ 미세환경 신호와 유전자 발현을 연결하는 동적 통신망 복원 완성!</p>
        </div>
      </div>
    </div>
  </div>
    """
]

# Write Deck 2
deck2_content = HTML_HEAD.format(title="Track 03 Masterclass Part 2: NicheNet, Spatial CCC & Translational Circuits") + "".join(deck2_slides) + HTML_FOOT.format(slide_count=len(deck2_slides))
deck2_path = os.path.join(SLIDES_DIR, "02_NicheNet_Spatial_and_Translational_Circuits_Deck.html")
with open(deck2_path, "w", encoding="utf-8") as fp:
    fp.write(deck2_content)
print(f"✅ Generated Slide Deck 2: {deck2_path} ({len(deck2_content):,} chars)")

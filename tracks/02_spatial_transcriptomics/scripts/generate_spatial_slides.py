#!/usr/bin/env python3
"""
generate_spatial_slides.py
Generates the two interactive HTML presentation slide decks for Track 02:
[LEGACY NOTICE] 이 스크립트는 구세대 덱(01_Spatial_Foundations_Deck / 02_Spatial_Deconvolution_and_Niche_Deck)을
생성합니다. 현재 정본 덱은 01_Spatial_Foundations_and_Segmentation_Deck.html 및
02_Spatial_Deconvolution_and_Niches_Deck.html이며, 구세대 덱은 backups/legacy_decks/에 보관되어 있습니다.
1. 01_Spatial_Foundations_Deck.html (Chapters 00-04, SVGs 01-03)
2. 02_Spatial_Deconvolution_and_Niche_Deck.html (Chapters 05-09, SVGs 04-06)
"""

import os
import re

TRACK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # ENG-025: 스크립트 위치 기준(이식 가능)
IMAGES_DIR = os.path.join(TRACK_DIR, "docs", "images")
SLIDES_DIR = os.path.join(TRACK_DIR, "slides")

os.makedirs(SLIDES_DIR, exist_ok=True)

def load_svg(filename: str) -> str:
    path = os.path.join(IMAGES_DIR, filename)
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    clean_svg = re.sub(r'<\?xml[^>]*\?>', '', content).strip()
    return clean_svg

svg_01 = load_svg("01_spatial_technologies_resolution_comparison.svg")
svg_02 = load_svg("02_coordinate_registration_and_spatialdata.svg")
svg_03 = load_svg("03_cellpose_vs_baysor_segmentation.svg")
svg_04 = load_svg("04_cell2location_deconvolution_model.svg")
svg_05 = load_svg("05_morans_i_and_spatial_autocorrelation.svg")
svg_06 = load_svg("06_spatial_niche_and_commot_signaling.svg")

COMMON_CSS = """
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
    --code-bg: #090d16;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }
  body { background: var(--bg-primary); color: var(--text-main); overflow: hidden; height: 100vh; display: flex; flex-direction: column; }
  
  #presentation { position: relative; flex: 1; width: 100%; height: 100%; overflow: hidden; }
  .slide { position: absolute; top: 0; left: 0; width: 100%; height: 100%; opacity: 0; visibility: hidden; transition: opacity 0.35s ease, transform 0.35s ease; transform: scale(0.97); display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 40px 60px; }
  .slide.active { opacity: 1; visibility: visible; transform: scale(1); }
  
  .slide-content { max-width: 1200px; width: 100%; max-height: 85vh; overflow-y: auto; background: var(--card-bg); border: 1px solid rgba(255,255,255,0.08); border-radius: 16px; padding: 40px 48px; box-shadow: 0 20px 50px rgba(0,0,0,0.5); backdrop-filter: blur(12px); }
  .slide-header { border-bottom: 2px solid var(--accent-glow); padding-bottom: 16px; margin-bottom: 24px; display: flex; justify-content: space-between; align-items: baseline; }
  .slide-header h2 { font-size: 2rem; color: var(--accent); font-weight: 700; letter-spacing: -0.5px; }
  .slide-header .tag { font-size: 0.85rem; padding: 4px 12px; border-radius: 20px; background: rgba(56, 189, 248, 0.15); color: var(--accent); font-weight: 600; text-transform: uppercase; }
  .doc-ref { font-size: 0.85rem; color: #38bdf8; background: rgba(56, 189, 248, 0.12); border: 1px solid rgba(56, 189, 248, 0.3); border-radius: 6px; padding: 5px 12px; margin-bottom: 16px; display: inline-flex; align-items: center; gap: 8px; font-weight: 600; }
  
  p, li { font-size: 1.1rem; line-height: 1.65; color: #cbd5e1; margin-bottom: 12px; }
  strong { color: #ffffff; }
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; margin-top: 20px; }
  .grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; margin-top: 20px; }
  .card { background: rgba(15, 23, 42, 0.6); border: 1px solid rgba(255,255,255,0.06); border-radius: 12px; padding: 20px; }
  .card.accent { border-left: 4px solid var(--accent); }
  .card.highlight { border-left: 4px solid var(--highlight); }
  .card.warning { border-left: 4px solid var(--warning); }
  .card.success { border-left: 4px solid var(--success); }
  .card h3 { font-size: 1.25rem; color: #fff; margin-bottom: 10px; }
  
  table { width: 100%; border-collapse: collapse; margin-top: 16px; font-size: 0.95rem; }
  th, td { padding: 12px 16px; text-align: left; border-bottom: 1px solid rgba(255,255,255,0.08); }
  th { background: rgba(56, 189, 248, 0.1); color: var(--accent); font-weight: 600; }
  
  pre { background: var(--code-bg); border: 1px solid rgba(255,255,255,0.1); border-radius: 8px; padding: 16px; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size: 0.9rem; color: #38bdf8; overflow-x: auto; margin-top: 12px; }
  code { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 4px; color: #38bdf8; font-size: 0.95em; }
  
  #controls { height: 60px; background: var(--bg-secondary); border-top: 1px solid rgba(255,255,255,0.08); display: flex; justify-content: space-between; align-items: center; padding: 0 30px; }
  .nav-btn { background: rgba(255,255,255,0.08); border: 1px solid rgba(255,255,255,0.12); color: #fff; padding: 8px 18px; border-radius: 8px; cursor: pointer; font-size: 0.9rem; transition: all 0.2s; }
  .nav-btn:hover { background: var(--accent); color: var(--bg-primary); }
  #counter { font-size: 0.95rem; color: var(--text-muted); font-weight: 600; }
  .progress-bar { position: absolute; top: 0; left: 0; height: 3px; background: var(--accent); width: 0%; transition: width 0.3s ease; }

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
    font-size: 0.95rem;
  }
  .diagram-drawer .svg-container svg {
    width: 100%;
    height: auto;
    max-height: 380px;
    margin-top: 10px;
    border-radius: 8px;
  }
"""

COMMON_JS = """
<script>
  const slides = document.querySelectorAll('.slide');
  const counter = document.getElementById('counter');
  const progressBar = document.getElementById('progressBar');
  let current = 0;

  function showSlide(index) {
    slides.forEach((s, i) => {
      s.classList.toggle('active', i === index);
    });
    counter.textContent = `Slide ${index + 1} / ${slides.length}`;
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
</script>
"""

# Build Deck 1
deck_1_html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Track 02: Spatial Transcriptomics Foundations & Technologies</title>
<style>{COMMON_CSS}</style>
</head>
<body>
<div class="progress-bar" id="progressBar"></div>
<div id="presentation">

  <!-- Slide 1: Title -->
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 60px 40px;">
      <span class="tag" style="font-size: 1rem; padding: 6px 16px;">Track 02 Masterclass • Part 1</span>
      <div style="margin-top: 15px;"><span class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/00_SPATIAL_ROADMAP.md &amp; 01</span></div>
      <h1 style="font-size: 2.8rem; color: var(--accent); margin: 20px 0 16px; font-weight: 800;">전립선암 공간전사체학 기초 및 물리·생체 광학 원리</h1>
      <h3 style="font-size: 1.4rem; color: #94a3b8; font-weight: 400; margin-bottom: 30px;">조직 해리 편향 극복, Gleason 병리 아키텍처, 아베 회절 한계, OME-Zarr 및 세포 분할</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 30px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08);">
        <p style="margin-bottom: 6px;">🎯 <strong>타깃 질환</strong>: 전립선 선암 (Primary PCa Gleason 3/4/5) &amp; 거세저항성 전이암 (mCRPC)</p>
        <p style="margin-bottom: 6px;">🔬 <strong>핵심 모달리티</strong>: 10x Visium, Visium HD, 10x Xenium, NanoString CosMx, MERSCOPE</p>
        <p style="margin-bottom: 0;">🧬 <strong>학습 중점</strong>: 조직 위상 보존, Affine 병리 정합, SpatialData 컨테이너, Cellpose vs Baysor 분할</p>
      </div>
    </div>
  </div>

  <!-- Slide 2: Why Spatial Transcriptomics? -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/00_SPATIAL_ROADMAP.md § 1</div>
      <div class="slide-header">
        <h2>1. 전립선암 연구에서 공간전사체학의 필수성</h2>
        <span class="tag">Clinical Rationale</span>
      </div>
      <p>단일세포 해리 시퀀싱(sc/snRNA)은 세포를 파쇄하여 물리적 접촉과 조직 구조를 완전히 파괴합니다. 전립선암의 진정한 악성도는 <strong>공간적 구조(Spatial Topology)</strong>에 내재되어 있습니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. Gleason 패턴 &amp; Cribriform 보존</h3>
          <p>Gleason 3(분리 관강), Gleason 4(체상 구조 Cribriform), Gleason 5(시트상 침윤) 병리 소포의 전사체를 형태학적 위치와 1:1로 직접 매핑합니다.</p>
        </div>
        <div class="card warning">
          <h3>2. 해리 편향 및 CAF 기질 장벽 극복</h3>
          <p>효소 분해 시 파괴되기 쉬운 치밀한 myCAF 결합조직 피막과 ECM 장벽을 온전히 보존하여 실제 조직 내 세포 비율을 왜곡 없이 포획합니다.</p>
        </div>
        <div class="card highlight">
          <h3>3. 면역 배제 (Immune Exclusion) 입증</h3>
          <p>CD8+ T 세포와 NK 세포가 종양 내부로 침투하지 못하고 기질 외곽에 갇히는 면역 콜드(Cold) 종양의 물리적 차단 기전을 시각적으로 검증합니다.</p>
        </div>
      </div>
    </div>
  </div>

  <!-- Slide 3: Spatial Technologies & Physics -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/01_SPATIAL_TECHNOLOGIES_AND_PHYSICS.md § 1~3</div>
      <div class="slide-header">
        <h2>2. 공간전사체 기술 분류 및 생물물리학적 한계</h2>
        <span class="tag">Biophysics &amp; Optics</span>
      </div>
      <p>공간전사체 기술은 <strong>시퀀싱 기반 어레이(Sequencing-based)</strong>와 <strong>이미징 기반 인시츄(Imaging-based)</strong>로 양분되며, 각각 분해능과 유전자 커버리지 간의 고유한 트레이드오프를 가집니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>시퀀싱 기반 (Visium / Visium HD)</h3>
          <ul>
            <li><strong>Visium Standard</strong>: 55 µm 스팟 (스팟당 5~15개 세포 혼합, 디콘볼루션 필수)</li>
            <li><strong>Visium HD</strong>: 2 µm 무간극 정사각 비드 어레이 (연속적 전장 전사체 포획)</li>
            <li><strong>장점</strong>: 편향 없는 전장 전사체(~20,000 유전자), 신규 바이오마커 발굴 가능</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>이미징 기반 (Xenium / CosMx / MERSCOPE)</h3>
          <ul>
            <li><strong>Xenium</strong>: 패들록 프로브 + 롤링 서클 증폭(RCA) 형광 볼(Rolony) 디코딩</li>
            <li><strong>광학 한계</strong>: 아베 회절 한계 (d = λ / 2NA ≈ 250 nm), 광학 크라우딩(Crowding) 현상</li>
            <li><strong>장점</strong>: 단일 분자 광학 분해능, 극도로 높은 검출 민감도와 단일 세포 할당력</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 공간전사체 기술별 해상도 및 검출 원리 비교 다이어그램 (Click to toggle)</summary>
        <div class="svg-container">
          {svg_01}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 4: Histopathology & Coordinate Systems -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/02_HISTOLOGY_AND_COORDINATE_SYSTEMS.md § 1~4</div>
      <div class="slide-header">
        <h2>3. 병리 영상 정합과 멀티모달 공간 좌표계</h2>
        <span class="tag">SpatialData &amp; OME-Zarr</span>
      </div>
      <p>고해상도 H&E/IF 병리 영상 픽셀 좌표계(Pixel Space)와 분자 발현 물리 좌표계(Physical Micron Space) 간의 정밀 정합은 공간 분석의 기초입니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>2차원 아핀 변환 (Affine Transformation)</h3>
          <ul>
            <li>스케일링, 회전, 이동(Translation), 전단(Shearing)을 3x3 동차 행렬로 표현</li>
            <li>Space Ranger <code>scalefactors_json.json</code>: 스팟 중심 좌표와 축소 영상 1:1 보정</li>
            <li>단일 세포 중심 좌표와 조직 윤곽선 폴리곤의 실시간 좌표 변환</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>SpatialData &amp; OME-Zarr 표준</h3>
          <ul>
            <li><strong>Images</strong>: 피라미드형 다해상도 래스터 영상 (Zarr 청크 저장)</li>
            <li><strong>Points</strong>: 분자 수준 단일 전사체 2D/3D 좌표 (Parquet/Arrow 백엔드)</li>
            <li><strong>Shapes</strong>: 세포 경계 분할 폴리곤 (Shapely/GeoPandas)</li>
            <li><strong>Tables</strong>: AnnData 호환 전사체 발현 카운트 행렬</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 병리 영상 정합 및 SpatialData OME-Zarr 표준 아키텍처 다이어그램</summary>
        <div class="svg-container">
          {svg_02}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 5: Spatial QC & Preprocessing -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/03_SPATIAL_QC_AND_PREPROCESSING.md § 1~3</div>
      <div class="slide-header">
        <h2>4. 공간전사체 특화 품질 관리(QC)와 노이즈 정제</h2>
        <span class="tag">Quality Control</span>
      </div>
      <p>공간전사체 데이터는 FFPE 조직 손상, 비특이적 프로브 결합, 광학 노이즈로 인해 표준 scRNA-seq과 다른 엄밀한 QC 잣대가 요구됩니다.</p>
      <div class="grid-3">
        <div class="card warning">
          <h3>1. 조직 마스킹 &amp; 유리 결합 제거</h3>
          <p>H&E Otsu 역치 분할을 통해 유리 슬라이드 표면에 비특이적으로 결합한 비조직(Off-tissue) 스팟 및 경계면 조직 파편을 완벽히 필터링합니다.</p>
        </div>
        <div class="card accent">
          <h3>2. 음성 대조군 프로브 &amp; FDR 계산</h3>
          <p>Xenium Negative Control 프로브 및 블랭크 바코드 카운트를 기반으로 분자 검출 위양성률(FDR &lt; 0.05)을 동적으로 산출하고 저신뢰 분자를 제거합니다.</p>
        </div>
        <div class="card highlight">
          <h3>3. FFPE 미토콘드리아 한계치</h3>
          <p>프로브 기반 FFPE 시퀀싱에서는 Poly-A 효소 역전사가 생략되므로 Mito%가 15~20% 미만이어야 하며, 과도한 카운트는 세포 자멸사 부위를 반영합니다.</p>
        </div>
      </div>
    </div>
  </div>

  <!-- Slide 6: Cell Segmentation Algorithms -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/04_CELL_SEGMENTATION_ALGORITHMS.md § 1~4</div>
      <div class="slide-header">
        <h2>5. 단일 세포 분할: 영상 기반 vs 전사체 기반</h2>
        <span class="tag">Segmentation Algorithms</span>
      </div>
      <p>하위 단일 세포 해석의 정확도는 분할 경계(Cell Boundary)의 무결성에 직결됩니다. 최신 알고리즘은 영상 단독 분할에서 분자 포인트 결합 모델로 진화하고 있습니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>Cellpose (딥러닝 그래디언트 플로우)</h3>
          <ul>
            <li>DAPI 핵 염색 영상에서 세포 중심을 향한 확산 벡터장(Vector Field) 예측</li>
            <li>형태가 불규칙하거나 밀집된 전립선 종양 상피세포 분할에 강력한 일반화 성능</li>
            <li><strong>한계</strong>: 핵 외곽 세포질(Cytoplasm) 경계 추정이 부정확하여 전사체 누락 발생</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>Baysor (마르코프 랜덤 필드 &amp; 포인트 프로세스)</h3>
          <ul>
            <li>개별 전사체의 2D 공간 좌표와 유전자 신원(Gene Identity)을 직접 모델링</li>
            <li>DAPI 핵 사전확률을 시작점으로, 주변 분자 구성의 유사도를 학습하여 경계 최적화</li>
            <li><strong>장점</strong>: 광학적으로 염색되지 않은 세포질 및 수지상 돌기까지 전사체 할당 성공</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 Cellpose vs Baysor 세포 분할 알고리즘 비교 다이어그램</summary>
        <div class="svg-container">
          {svg_03}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 7: Technical Comparison Matrix -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/01 § 4 &amp; 03 § 4</div>
      <div class="slide-header">
        <h2>6. 공간전사체 플랫폼별 핵심 스펙 비교 매트릭스</h2>
        <span class="tag">Platform Benchmark</span>
      </div>
      <table>
        <thead>
          <tr>
            <th>기술 플랫폼</th>
            <th>해상도 (Resolution)</th>
            <th>포획 유전자 수</th>
            <th>검출 민감도</th>
            <th>시료 적합성</th>
            <th>주요 분석 도구</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>10x Visium Standard</strong></td>
            <td>55 µm (5~15 세포 혼합)</td>
            <td>전장 전사체 (~20,000)</td>
            <td>중간 (Spot당 수천 UMI)</td>
            <td>Fresh Frozen &amp; FFPE</td>
            <td><code>Space Ranger</code>, <code>cell2location</code></td>
          </tr>
          <tr>
            <td><strong>10x Visium HD</strong></td>
            <td>2 µm 무간극 정사각 빈</td>
            <td>전장 전사체 (~20,000)</td>
            <td>Bin당 희소 (비닝 필요)</td>
            <td>FFPE 특화</td>
            <td><code>Space Ranger 3.0</code>, <code>SpatialData</code></td>
          </tr>
          <tr>
            <td><strong>10x Xenium In Situ</strong></td>
            <td>&lt; 0.2 µm (단일 분자)</td>
            <td>표적 패널 (300 ~ 5,000)</td>
            <td>극도로 높음 (단일 Rolony)</td>
            <td>FFPE &amp; Fresh Frozen</td>
            <td><code>Baysor</code>, <code>Cellpose</code>, <code>Squidpy</code></td>
          </tr>
          <tr>
            <td><strong>NanoString CosMx</strong></td>
            <td>&lt; 0.1 µm (아세포 해상도)</td>
            <td>표적 패널 (1,000 ~ 6,000)</td>
            <td>매우 높음 (RNA + 단백질)</td>
            <td>FFPE 특화</td>
            <td><code>CosMx AtoMx</code>, <code>InSituType</code></td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>

  <!-- Slide 8: Part 1 Summary & Next Step -->
  <div class="slide">
    <div class="slide-content" style="text-align: center;">
      <div class="slide-header" style="justify-content: center;">
        <h2>Track 02 Part 1 총괄 요약 및 Part 2 예고</h2>
      </div>
      <p style="font-size: 1.25rem; max-width: 900px; margin: 0 auto 30px;">
        Part 1에서는 전립선암 조직 병리의 공간적 보존, 플랫폼별 물리적 분해능, SpatialData 멀티모달 컨테이너, 그리고 Cellpose/Baysor 세포 분할 메커니즘을 마스터하였습니다.
      </p>
      <div class="grid-2" style="text-align: left; max-width: 960px; margin: 0 auto;">
        <div class="card success">
          <h3>Part 1 완주 역량</h3>
          <p>✓ 시퀀싱(Visium) vs 이미징(Xenium) 물리적 분해능 한계 인지<br>✓ SpatialData OME-Zarr 표준 다해상도 영상 및 포인트 정합<br>✓ DAPI 그래디언트 플로우 및 Baysor 공간 포인트 프로세스 분할</p>
        </div>
        <div class="card accent">
          <h3>Part 2 예고: 디콘볼루션, 니치 및 신호망</h3>
          <p>cell2location 음이항 베이지안 역중첩 ──► Moran's I 공간 가변 유전자 ──► Squidpy 세포 니치 클러스터링 ──► COMMOT 최적 수송 리간드-수용체 신호 플럭스</p>
        </div>
      </div>
    </div>
  </div>

</div>

<div id="controls">
  <button class="nav-btn" id="prevBtn">◀ 이전 (Prev)</button>
  <div id="counter">Slide 1 / 8</div>
  <button class="nav-btn" id="nextBtn">다음 (Next) ▶</button>
</div>

{COMMON_JS}
</body>
</html>
"""

# Build Deck 2
deck_2_html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Track 02: Spatial Deconvolution, Microenvironments & Signaling</title>
<style>{COMMON_CSS}</style>
</head>
<body>
<div class="progress-bar" id="progressBar"></div>
<div id="presentation">

  <!-- Slide 1: Title -->
  <div class="slide active">
    <div class="slide-content" style="text-align: center; padding: 60px 40px;">
      <span class="tag" style="font-size: 1rem; padding: 6px 16px;">Track 02 Masterclass • Part 2</span>
      <div style="margin-top: 15px;"><span class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md &amp; 08</span></div>
      <h1 style="font-size: 2.8rem; color: var(--accent); margin: 20px 0 16px; font-weight: 800;">공간 디콘볼루션, 세포 니치 및 리간드-수용체 신호망</h1>
      <h3 style="font-size: 1.4rem; color: #94a3b8; font-weight: 400; margin-bottom: 30px;">cell2location 베이지안 추론, Moran's I 공간 가변 유전자, 이웃 그래프 니치 및 COMMOT 최적 수송</h3>
      <div style="display: inline-block; text-align: left; background: rgba(15,23,42,0.6); padding: 20px 30px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.08);">
        <p style="margin-bottom: 6px;">🎯 <strong>타깃 질환</strong>: 전립선암 기질-면역 장벽, 골전이 침윤 전단(Invasive Front), 림프절 전이</p>
        <p style="margin-bottom: 6px;">🔬 <strong>핵심 알고리즘</strong>: cell2location, SpatialDE, Squidpy Neighborhood Graph, COMMOT</p>
        <p style="margin-bottom: 0;">🧬 <strong>학습 중점</strong>: 스팟 세포 조성 복원, 공간 자기상관 지수, 니치 미세환경, 거리 감쇄 신호 플럭스</p>
      </div>
    </div>
  </div>

  <!-- Slide 2: Spatial Deconvolution Mathematics -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md § 1~3</div>
      <div class="slide-header">
        <h2>1. 스팟 디콘볼루션 수학: cell2location 계층적 모델</h2>
        <span class="tag">Bayesian Deconvolution</span>
      </div>
      <p>10x Visium의 각 스팟은 여러 세포가 섞인 혼합체(Mixture)입니다. 단일핵 참조 전사체(snRNA-seq)를 사전확률로 활용하여 스팟당 절대 세포 밀도를 역추정합니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>음이항(Negative Binomial) 관측 모델</h3>
          <ul>
            <li>y<sub>sg</sub> ~ NB(μ<sub>sg</sub>, α<sub>g</sub>): 유전자별 과분산(Overdispersion) 완벽 반영</li>
            <li>기대 발현량: μ<sub>sg</sub> = m<sub>s</sub> ∑<sub>k</sub> w<sub>sk</sub> g<sub>kg</sub> + l<sub>s</sub></li>
            <li>w<sub>sk</sub>: 스팟 s 내 세포 유형 k의 절대 세포 밀도(Cell Abundance)</li>
            <li>g<sub>kg</sub>: sc/snRNA 참조 아틀라스에서 사전 추정된 세포 유형별 참조 발현 프로파일</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>베이지안 추론 및 기술적 노이즈 분리</h3>
          <ul>
            <li><strong>m<sub>s</sub> (스팟 민감도)</strong>: 조직 두께 및 RNA 역전사 효율의 기술적 편차 보정</li>
            <li><strong>l<sub>s</sub> (배경 노이즈)</strong>: 조직 외 비특이적 결합 및 광학 노이즈 분리</li>
            <li>변분 추론(Variational Inference, Pyro)을 통해 수만 개 스팟의 사후확률 분포 고속 수렴</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 cell2location 계층적 베이지안 스팟 디콘볼루션 모델 다이어그램</summary>
        <div class="svg-container">
          {svg_04}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 3: Spatially Variable Genes -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/06_SPATIALLY_VARIABLE_GENES.md § 1~4</div>
      <div class="slide-header">
        <h2>2. 공간적 자기상관 및 가변 유전자(SVG) 탐색</h2>
        <span class="tag">Spatial Autocorrelation</span>
      </div>
      <p>단순 발현량 변동(HVG)을 넘어, 2차원 조직 평면에서 <strong>공간적 구조(Spatial Pattern)</strong>를 형성하는 핵심 기능 유전자를 통계적으로 규명합니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>공간 가중치 행렬 W 와 Moran's I</h3>
          <ul>
            <li><strong>Moran's I</strong>: 글로벌 공간 자기상관 지수 (-1 ~ +1)</li>
            <li>I &gt; 0: 인접 스팟 간 발현 유사도 높음 (공간적 군집화 클러스터링)</li>
            <li>k-NN 기반 또는 가우스 반경 가중치: W<sub>ij</sub> = exp(-d<sub>ij</sub>² / 2σ²)</li>
            <li>Z-score 및 순열 검정(Permutation Test) 기반 통계적 유의성(FDR) 산출</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>SpatialDE (가우스 과정 회귀 모델)</h3>
          <ul>
            <li>공간 공분산(K<sub>spatial</sub>)과 독립 노이즈(I)의 분산 성분 분해</li>
            <li>p-값 우도비 검정(Log-likelihood Ratio Test)으로 무작위 노이즈 배제</li>
            <li>전립선암 Cribriform 선암종 마커(HIF1A, VEGFA) 및 myCAF 콜라겐의 강한 공간 자기상관 포획</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 공간 자기상관 Moran's I 및 공간 가중치 기하학 다이어그램</summary>
        <div class="svg-container">
          {svg_05}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 4: Cellular Niches & Microenvironments -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/07_CELLULAR_NICHES_AND_MICROENVIRONMENTS.md § 1~3</div>
      <div class="slide-header">
        <h2>3. 세포 니치(Cellular Niche)와 공간 미세환경</h2>
        <span class="tag">Niche Clustering</span>
      </div>
      <p>단일 세포 유형 자체보다, 세포가 <strong>"어떤 이웃 세포들과 물리적으로 둘러싸여 있는가"</strong>가 치료 반응성과 예후를 결정합니다.</p>
      <div class="grid-3">
        <div class="card accent">
          <h3>1. 이웃 윈도우 벡터화</h3>
          <p>각 세포를 중심으로 반경 r(예: 30 µm) 내 존재하는 주변 세포 유형의 빈도 벡터를 계산하여 국소 미세환경 프로파일을 구축합니다.</p>
        </div>
        <div class="card warning">
          <h3>2. 전립선암 4대 핵심 니치</h3>
          <p>• <strong>Cancer Core</strong>: 순수 Luminal AR+ 암세포 고밀도<br>• <strong>Desmoplastic Rim</strong>: myCAF 피막 장벽<br>• <strong>TLS Hub</strong>: CD8 T/B 세포 3차 림프구 구조<br>• <strong>Bone Front</strong>: 파골/조골 세포 침윤 경계면</p>
        </div>
        <div class="card highlight">
          <h3>3. 임상 예후 직결성</h3>
          <p>동일한 AR+ 암세포라도 Desmoplastic Rim에 인접할수록 TGF-β 유도 EMT 표현형을 획득하며 안드로겐 차단 요법에 조기 내성을 나타냅니다.</p>
        </div>
      </div>
    </div>
  </div>

  <!-- Slide 5: Spatial Ligand-Receptor Signaling -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md § 1~4</div>
      <div class="slide-header">
        <h2>4. 거리 감쇄 리간드-수용체 신호 플럭스 (COMMOT)</h2>
        <span class="tag">Optimal Transport &amp; Flux</span>
      </div>
      <p>공간 정보가 없는 기존 도구(CellChat)의 거리 무시 위양성을 극복하고, <strong>유한 확산 거리(Diffusion Limit)</strong>를 제약 조건으로 하는 최적 수송 모델을 적용합니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>COMMOT 최적 수송 (Optimal Transport) 수학</h3>
          <ul>
            <li>송신 세포 발현(L<sub>i</sub>)에서 수신 세포 수용체(R<sub>j</sub>)로의 질량 수송 비용 최소화</li>
            <li>거리 비용 행렬: C<sub>ij</sub> = d<sub>ij</sub> (유효 확산 반경 d<sub>0</sub> 초과 시 무한대 비용)</li>
            <li><strong>엔트로피 정규화</strong>: Sinkhorn 알고리즘을 통한 고속 수렴</li>
            <li>공간 신호 방향성 벡터장(Directional Signaling Vector Field) 복원</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>전립선암 핵심 신호 축 검증</h3>
          <ul>
            <li><strong>CAF ──► 암세포 (TGFB1-TGFBR2)</strong>: 침윤 선단 암세포의 침습 촉진</li>
            <li><strong>CAF ──► CD8 T 세포 (CXCL12-CXCR4)</strong>: T 세포 종양 침투 차단 기질 트래핑</li>
            <li><strong>암세포 ──► 조골세포 (BMP6, WNT7B)</strong>: 골전이성 조골성 병변(Osteoblastic lesion) 유도</li>
          </ul>
        </div>
      </div>
      <details class="diagram-drawer" open>
        <summary>📊 공간 미세환경 니치 클러스터링 및 COMMOT 신호 플럭스 다이어그램</summary>
        <div class="svg-container">
          {svg_06}
        </div>
      </details>
    </div>
  </div>

  <!-- Slide 6: Prostate Cancer 3D Pathology & Clinical Translation -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md § 1~3</div>
      <div class="slide-header">
        <h2>5. 3차원 병리 재구성 및 임상 중개 연구</h2>
        <span class="tag">Clinical Translation</span>
      </div>
      <p>단일 슬라이드 2D 분석에서 연속 절편(Serial Sections)의 3D 공간 복원으로 확장하여 종양 체적 내 클론 분기와 치료 저항성 아키텍처를 규명합니다.</p>
      <div class="grid-3">
        <div class="card">
          <h3>1. 3D 연속 절편 비강체 정합</h3>
          <p>조직 절편 회전 및 국소 신장 왜곡을 B-spline 탄성 변환으로 정밀 교정하여 3차원 Gleason 기하학을 복원합니다.</p>
        </div>
        <div class="card accent">
          <h3>2. Cribriform 저산소 중심부</h3>
          <p>체상 구조 선암종 중심부의 저산소증(HIF1A) 및 대사 재프로그래밍 구배를 3차원 등고선 맵으로 시각화합니다.</p>
        </div>
        <div class="card warning">
          <h3>3. 골전이 틈새(Niche) 공략</h3>
          <p>뼈 전이 검체에서 골수 기질-조골세포-신경내분비 암세포 간의 공간 병치 및 이중 저항성 클론의 공간적 탈출구를 규명합니다.</p>
        </div>
      </div>
    </div>
  </div>

  <!-- Slide 7: Production Pipelines & Capstone Rubric -->
  <div class="slide">
    <div class="slide-content">
      <div class="doc-ref">📖 대응 교재: tracks/02_spatial_transcriptomics/pipelines/ &amp; docs/09 § 4</div>
      <div class="slide-header">
        <h2>6. 재현 가능한 프로덕션 파이프라인 및 실습</h2>
        <span class="tag">Executable Pipelines</span>
      </div>
      <p>본 트랙은 이론에 그치지 않고, 100% 독립 실행 가능한 2종의 표준 파이프라인을 제공하여 즉시 실제 데이터셋 분석에 적용할 수 있습니다.</p>
      <div class="grid-2">
        <div class="card">
          <h3>Pipeline 01: 스팟 QC 및 디콘볼루션</h3>
          <ul>
            <li><code>01_spatial_qc_and_deconvolution_pipeline.py</code></li>
            <li>50x50 그리드 2,500 스팟 코호트 시뮬레이션 및 Otsu 조직 마스킹</li>
            <li>Moran's I 통계 검정 기반 상위 367개 공간 가변 유전자 규명</li>
            <li>8대 세포 유형 Multiplicative Poisson solver 역중첩 (Pearson r &gt; 0.998)</li>
          </ul>
        </div>
        <div class="card accent">
          <h3>Pipeline 02: Xenium 세포 분할 &amp; COMMOT 신호망</h3>
          <ul>
            <li><code>02_xenium_subcellular_niche_and_signaling_pipeline.py</code></li>
            <li>50,000개 서브셀룰러 전사체 및 KDTree 공간 할당 (85.1% 유효 할당)</li>
            <li>r=30 µm 반경 기반 4대 공간 니치 클러스터링</li>
            <li>지수 거리 감쇄 d<sub>0</sub>=30 µm 제약 하의 CAF-암세포 EMT 신호 플럭스 계산</li>
          </ul>
        </div>
      </div>
    </div>
  </div>

  <!-- Slide 8: Grand Ecosystem Alignment -->
  <div class="slide">
    <div class="slide-content" style="text-align: center;">
      <div class="slide-header" style="justify-content: center;">
        <h2>Track 02 완주 및 그랜드 커리큘럼 로드맵 정렬</h2>
      </div>
      <p style="font-size: 1.25rem; max-width: 900px; margin: 0 auto 30px;">
        축하합니다! 단일 세포 해리 시퀀싱(Track 00)과 후성유전학 Multiome(Track 01)에 이어, 온전한 조직 아키텍처를 복원하는 <strong>공간전사체학(Track 02) 마스터 과정</strong>을 성공적으로 완주하였습니다.
      </p>
      <div class="grid-3" style="text-align: left; max-width: 1050px; margin: 0 auto;">
        <div class="card success">
          <h3>Track 00 (완료)</h3>
          <p>동결 조직 snRNA-seq, Pre-mRNA 인트론, UMAP 매니폴드, inferCNV 악성 선별</p>
        </div>
        <div class="card success">
          <h3>Track 01 (완료)</h3>
          <p>10x Multiome, Tn5 화학, WNN 멀티모달 통합, SCENIC+ Triplet eGRN</p>
        </div>
        <div class="card accent">
          <h3>Track 02 (완료)</h3>
          <p>공간전사체(Visium/Xenium), SpatialData, cell2location, Moran's I, COMMOT 신호 플럭스</p>
        </div>
      </div>
      <div style="margin-top: 30px;">
        <span class="tag" style="font-size: 0.95rem; padding: 8px 20px;">Next Track ──► Track 03: 종양미세환경 &amp; 세포 간 상호작용 (CellChat v2, LIANA+, NicheNet v2)</span>
      </div>
    </div>
  </div>

</div>

<div id="controls">
  <button class="nav-btn" id="prevBtn">◀ 이전 (Prev)</button>
  <div id="counter">Slide 1 / 8</div>
  <button class="nav-btn" id="nextBtn">다음 (Next) ▶</button>
</div>

{COMMON_JS}
</body>
</html>
"""

# Write out the slide decks
deck_1_path = os.path.join(SLIDES_DIR, "01_Spatial_Foundations_Deck.html")
deck_2_path = os.path.join(SLIDES_DIR, "02_Spatial_Deconvolution_and_Niche_Deck.html")

with open(deck_1_path, "w", encoding="utf-8") as f:
    f.write(deck_1_html)
print(f"✓ Generated Deck 1: {deck_1_path} ({len(deck_1_html):,} bytes)")

with open(deck_2_path, "w", encoding="utf-8") as f:
    f.write(deck_2_html)
print(f"✓ Generated Deck 2: {deck_2_path} ({len(deck_2_html):,} bytes)")

import os
import re

svg_dir = 'tracks/04_cancer_genomics_clones/docs/images'

def get_svg(filename):
    with open(os.path.join(svg_dir, filename), 'r', encoding='utf-8') as f:
        return f.read().strip()

svg01 = get_svg('01_prostate_cancer_hallmark_cna_landscape.svg')
svg02 = get_svg('02_infercnv_hmm_smoothing_workflow.svg')
svg03 = get_svg('03_copykat_gmm_diploid_aneuploid_separation.svg')
svg04 = get_svg('04_numbat_haplotype_phased_baf_and_cnloh.svg')
svg05 = get_svg('05_clonal_phylogeny_tree_and_subclone_evolution.svg')
svg06 = get_svg('06_therapy_selection_and_metastatic_seeding.svg')

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

js_template = """
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

print("Setup template strings ready.")

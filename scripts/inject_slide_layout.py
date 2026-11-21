#!/usr/bin/env python3
"""
inject_slide_layout.py (ENG-040)
Marker-based, idempotent injection of a common layout layer into slide decks:
  - desktop: fit-to-viewport (zoom the slide content so it fits, min 0.7, then scrolls)
  - mobile (<=768px): vertical scroll mode (no scaling), 1-column grids, scrollable tables/code
  - SVG: max-width 100%, height auto, viewport-based max-height
Deck text/SVG are not modified. Re-running strips the previous block first (idempotent).
Run it AFTER the deck generators (generate_spatial_slides.py, generate_tme_slides.py, tracks/04 make_decks /
build_all_decks.py, tracks/05 build_decks.py) and after hand edits of decks.
Usage: python3 scripts/inject_slide_layout.py [deck.html ...]   (default: all slides/*.html and tracks/*/slides/*.html)
"""
import re, sys, glob, os
CSS_B, CSS_E = "/* SLIDE-LAYOUT:BEGIN:css */", "/* SLIDE-LAYOUT:END:css */"
JS_B, JS_E = "/* SLIDE-LAYOUT:BEGIN:js */", "/* SLIDE-LAYOUT:END:js */"
CSS = r"""
/* ---- desktop fit-to-viewport: .sl-fit wraps the slide content and is zoomed by the JS below ---- */
.slide { overflow: hidden; }
.slide-content { max-height: calc(100vh - 60px - 24px) !important; }   /* 60px = #controls */
.slide-content .sl-fit { transform-origin: top center; }
.svg-container svg, .slide-content svg { max-width: 100%; height: auto; }
.diagram-drawer .svg-container svg { max-height: min(310px, 42vh); }
.slide-content table { max-width: 100%; }
.slide-content pre { max-width: 100%; overflow-x: auto; }
/* ---- mobile: vertical scroll mode, no scaling ---- */
@media (max-width: 768px) {
  body { height: 100dvh; }
  .slide { padding: 12px 10px !important; justify-content: flex-start !important; overflow-y: auto; }
  .slide-content { max-height: none !important; overflow-y: visible !important; padding: 18px 14px !important; }
  .slide-content .sl-fit { zoom: 1 !important; }
  .slide-header { flex-direction: column; align-items: flex-start !important; gap: 6px; }
  .slide-header h2 { font-size: 1.3rem !important; }
  .grid-2, .grid-3, .grid-4 { grid-template-columns: 1fr !important; }
  p, li { font-size: 0.95rem !important; }
  .slide-content [style*="grid-template-columns"] { grid-template-columns: 1fr !important; }
  .slide-content [style*="display: flex"], .slide-content [style*="display:flex"] { flex-wrap: wrap !important; }
  .slide-content :not(pre) > code { white-space: normal; overflow-wrap: anywhere; }
  .slide-content .card { max-width: 100%; min-width: 0; }
  .slide-content table { display: block; overflow-x: auto; font-size: 0.8rem; }
  .slide-content pre { font-size: 0.74rem; }
  .slide-content h1 { font-size: 1.6rem !important; }
  .doc-ref { display: flex !important; max-width: 100%; white-space: normal; word-break: break-all; font-size: 0.75rem !important; }
  #controls { padding: 0 12px; }
  .diagram-drawer .svg-container svg { max-height: none; }
}
"""
JS = r"""
(function () {
  var MIN = 0.7;
  function wrap(sc) {
    if (sc.querySelector(':scope > .sl-fit')) return sc.querySelector(':scope > .sl-fit');
    var w = document.createElement('div'); w.className = 'sl-fit';
    while (sc.firstChild) w.appendChild(sc.firstChild);
    sc.appendChild(w); return w;
  }
  function fit() {
    var mobile = window.matchMedia('(max-width: 768px)').matches;
    document.querySelectorAll('.slide-content').forEach(function (sc) {
      var w = wrap(sc);
      if (mobile) { w.style.zoom = ''; return; }
      var slide = sc.closest('.slide');
      if (slide && !slide.classList.contains('active')) return;
      w.style.zoom = 1;
      var cs = getComputedStyle(sc);
      var avail = sc.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
      var need = w.scrollHeight;
      w.style.zoom = need > avail ? Math.max(MIN, avail / need) : 1;
    });
  }
  window.addEventListener('resize', fit);
  document.addEventListener('keydown', function () { setTimeout(fit, 0); });
  document.addEventListener('click', function () { setTimeout(fit, 0); setTimeout(fit, 400); });
  document.addEventListener('toggle', function () { setTimeout(fit, 0); }, true);
  if (document.readyState !== 'loading') setTimeout(fit, 0); else document.addEventListener('DOMContentLoaded', function () { setTimeout(fit, 0); });
  setTimeout(fit, 450);
})();
"""
def strip(t):
    t = re.sub(r"\n?" + re.escape(CSS_B) + r".*?" + re.escape(CSS_E) + r"\n?", "\n", t, flags=re.S)
    t = re.sub(r"\n?<script>\s*" + re.escape(JS_B) + r".*?" + re.escape(JS_E) + r"\s*</script>\n?", "\n", t, flags=re.S)
    return t
def inject(t):
    t = strip(t)
    t = t.replace("</style>", CSS_B + CSS + CSS_E + "\n</style>", 1) if False else re.sub(r"</style>(?![\s\S]*</style>)", lambda m: CSS_B + CSS + CSS_E + "\n</style>", t, count=1)
    t = t.replace("</body>", "<script>\n" + JS_B + JS + JS_E + "\n</script>\n</body>", 1)
    return t
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(ROOT, "slides", "*.html")) +
                                  glob.glob(os.path.join(ROOT, "tracks", "*", "slides", "*.html")))
    for f in files:
        s = open(f, encoding="utf-8").read(); n = inject(s)
        assert inject(n) == n, "not idempotent"
        open(f, "w", encoding="utf-8").write(n); print("injected", f)

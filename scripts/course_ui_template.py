"""
course_ui_template.py (PROTOTYPE, ENG-044)
One HTML shell + CSS + JS for the 11 Master Course builders. Builders pass only track configuration.

    render_page(title, sidebar_title, search_placeholder, theme_key, hero_html, sidebar_nav, chapters_html,
                sidebar_badge="", main_prefix="", main_suffix="")

Anchors kept for scripts/inject_global_track_switcher.py: a "</style>" (CSS), <div class="search-container">,
and <main class="main-content" id="mainContent"> / <div class="document-container">.
"""

CSS = r'''
  :root {
    --bg-body: #ffffff;
    --text-main: #334155;
    --text-dark: #0f172a;
    --bg-surface: #f8fafc;
    --border-color: #e2e8f0;
    --code-bg: #0f172a;
    --code-text: #e2e8f0;
    --inline-code-bg: #f1f5f9;
    --inline-code-text: #0284c7;
    --accent: #0284c7;
    --accent-hover: #0369a1;
    --accent-light: #e0f2fe;
    --sidebar-bg: #ffffff;
    --sidebar-border: #e2e8f0;
    --card-bg: #f8fafc;
    --card-bg-open: #f0f9ff;
    --highlight-bg: #fef08a;
    --math-bg: #f8fafc;
  }

  [data-theme="dark"] {
    --bg-body: #0b1120;
    --text-main: #cbd5e1;
    --text-dark: #f8fafc;
    --bg-surface: #1e293b;
    --border-color: #334155;
    --code-bg: #030712;
    --code-text: #f3f4f6;
    --inline-code-bg: #1e293b;
    --inline-code-text: #38bdf8;
    --accent: #38bdf8;
    --accent-hover: #7dd3fc;
    --accent-light: #0c4a6e;
    --sidebar-bg: #0f172a;
    --sidebar-border: #1e293b;
    --card-bg: #1e293b;
    --card-bg-open: #0f2847;
    --highlight-bg: #854d0e;
    --math-bg: #1e293b;
  }

  @page {
    size: A4;
    margin: 20mm 15mm 20mm 15mm;
    @bottom-center {
      content: counter(page);
      font-size: 9pt;
      color: #64748b;
    }
  }

  * { box-sizing: border-box; }
  html { scroll-behavior: smooth; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Apple SD Gothic Neo", "Malgun Gothic", sans-serif;
    line-height: 1.75;
    color: var(--text-main);
    background-color: var(--bg-body);
    margin: 0;
    padding: 0;
    transition: background-color 0.25s ease, color 0.25s ease;
  }

  /* Layout */
  .layout-wrapper {
    display: flex;
    min-height: 100vh;
  }

  /* Sidebar */
  aside.sidebar {
    width: 330px;
    min-width: 330px;
    background: var(--sidebar-bg);
    border-right: 1px solid var(--sidebar-border);
    position: sticky;
    top: 0;
    height: 100vh;
    overflow-y: auto;
    padding: 24px 20px;
    display: flex;
    flex-direction: column;
    gap: 18px;
    z-index: 100;
  }

  .sidebar-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding-bottom: 14px;
    border-bottom: 2px solid var(--accent);
  }

  .sidebar-title {
    font-size: 1.05rem;
    font-weight: 800;
    color: var(--text-dark);
    margin: 0;
    letter-spacing: -0.5px;
  }

  .theme-toggle-btn {
    background: var(--bg-surface);
    border: 1px solid var(--border-color);
    color: var(--text-main);
    padding: 6px 10px;
    border-radius: 8px;
    cursor: pointer;
    font-size: 0.85rem;
    font-weight: 600;
    transition: all 0.2s;
  }
  .theme-toggle-btn:hover {
    background: var(--accent-light);
    color: var(--accent);
  }

  /* Search Box */
  .search-container {
    position: relative;
    width: 100%;
  }
  .search-input {
    width: 100%;
    padding: 10px 14px 10px 36px;
    background: var(--bg-surface);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    color: var(--text-dark);
    font-size: 0.9rem;
    outline: none;
    transition: border-color 0.2s;
  }
  .search-input:focus {
    border-color: var(--accent);
  }
  .search-icon {
    position: absolute;
    left: 12px;
    top: 50%;
    transform: translateY(-50%);
    font-size: 0.9rem;
    color: #94a3b8;
  }
  .search-results-count {
    font-size: 0.8rem;
    color: var(--accent);
    margin-top: 4px;
    font-weight: 600;
    min-height: 16px;
  }

  /* Global Actions */
  .global-actions {
    display: flex;
    gap: 8px;
    flex-wrap: wrap;
  }
  .action-btn {
    flex: 1;
    min-width: 85px;
    padding: 7px 10px;
    background: var(--bg-surface);
    border: 1px solid var(--border-color);
    border-radius: 6px;
    color: var(--text-main);
    font-size: 0.8rem;
    font-weight: 600;
    cursor: pointer;
    transition: all 0.2s;
    text-align: center;
  }
  .action-btn:hover {
    background: var(--accent-light);
    color: var(--accent);
    border-color: var(--accent);
  }

  /* Sidebar Navigation Tree */
  nav.sidebar-nav {
    display: flex;
    flex-direction: column;
    gap: 12px;
  }
  .sidebar-chap-group {
    display: flex;
    flex-direction: column;
    gap: 6px;
  }
  .sidebar-chap-title {
    font-size: 0.88rem;
    font-weight: 700;
    color: var(--text-dark);
    text-decoration: none;
    line-height: 1.4;
    padding: 6px 8px;
    border-radius: 6px;
    transition: background 0.15s, color 0.15s;
  }
  .sidebar-chap-title:hover {
    background: var(--accent-light);
    color: var(--accent);
  }
  .sidebar-sec-list {
    list-style: none;
    margin: 0;
    padding: 0 0 0 12px;
    border-left: 2px solid var(--border-color);
    display: flex;
    flex-direction: column;
    gap: 4px;
  }
  .sidebar-sec-list li a {
    font-size: 0.8rem;
    color: #64748b;
    text-decoration: none;
    display: block;
    padding: 3px 6px;
    border-radius: 4px;
    transition: color 0.15s, background 0.15s;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  .sidebar-sec-list li a:hover {
    color: var(--accent);
    background: var(--bg-surface);
  }

  /* Main Content Area */
  main.main-content {
    flex: 1;
    max-width: 1040px;
    margin: 0 auto;
    padding: 40px 48px;
    overflow-x: hidden;
  }

  /* Hero Section */
  .course-hero {
    background: linear-gradient(135deg, #0284c7 0%, #0369a1 100%);
    color: #ffffff;
    padding: 48px 40px;
    border-radius: 16px;
    margin-bottom: 48px;
    box-shadow: 0 10px 25px -5px rgba(2, 132, 199, 0.25);
  }
  .course-hero h1 {
    font-size: 2.2rem;
    font-weight: 800;
    margin: 0 0 14px;
    letter-spacing: -0.5px;
    color: #ffffff;
  }
  .course-hero p.subtitle {
    font-size: 1.05rem;
    color: #e0f2fe;
    margin: 0 0 24px;
    line-height: 1.6;
  }
  .hero-badges {
    display: flex;
    gap: 10px;
    flex-wrap: wrap;
  }
  .hero-badge {
    background: rgba(255, 255, 255, 0.18);
    backdrop-filter: blur(8px);
    border: 1px solid rgba(255, 255, 255, 0.3);
    color: #ffffff;
    font-size: 0.82rem;
    font-weight: 600;
    padding: 5px 12px;
    border-radius: 20px;
  }

  /* Chapter Article Containers */
  article.chapter-container {
    margin-bottom: 64px;
    padding-bottom: 48px;
    border-bottom: 2px solid var(--border-color);
  }
  article.chapter-container:last-child {
    border-bottom: none;
  }
  .chapter-meta {
    display: flex;
    align-items: center;
    gap: 12px;
    margin-bottom: 16px;
  }
  .chapter-badge {
    background: var(--accent);
    color: #ffffff;
    font-size: 0.8rem;
    font-weight: 800;
    padding: 4px 12px;
    border-radius: 14px;
    text-transform: uppercase;
  }
  .source-ref {
    font-size: 0.8rem;
    color: #94a3b8;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
  }

  /* Typography */
  h1 { font-size: 2.0rem; font-weight: 800; color: var(--text-dark); margin: 30px 0 18px; letter-spacing: -0.5px; }
  h2 { font-size: 1.5rem; font-weight: 700; color: var(--text-dark); margin: 36px 0 16px; padding-bottom: 8px; border-bottom: 1px solid var(--border-color); scroll-margin-top: 24px; }
  h3 { font-size: 1.2rem; font-weight: 700; color: var(--text-dark); margin: 26px 0 12px; }
  h4 { font-size: 1.05rem; font-weight: 600; color: var(--text-dark); margin: 20px 0 10px; }
  p { margin: 0 0 16px; }
  ul, ol { margin: 0 0 16px; padding-left: 24px; }
  li { margin-bottom: 6px; }
  hr { border: none; border-top: 1px solid var(--border-color); margin: 40px 0; }
  blockquote {
    margin: 0 0 20px;
    padding: 16px 20px;
    background: var(--bg-surface);
    border-left: 4px solid var(--accent);
    border-radius: 0 8px 8px 0;
    color: var(--text-main);
  }
  blockquote p:last-child { margin-bottom: 0; }

  /* Tables */
  table {
    width: 100%;
    border-collapse: collapse;
    margin: 20px 0 28px;
    font-size: 0.9rem;
  }
  th, td {
    padding: 10px 14px;
    text-align: left;
    border: 1px solid var(--border-color);
  }
  th {
    background: var(--bg-surface);
    color: var(--text-dark);
    font-weight: 700;
  }
  tr:nth-child(even) {
    background: var(--bg-surface);
  }

  /* Code Blocks & Inline Code */
  pre {
    background: var(--code-bg);
    color: var(--code-text);
    padding: 18px 20px;
    border-radius: 10px;
    overflow-x: auto;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    font-size: 0.88rem;
    line-height: 1.6;
    margin: 16px 0 24px;
  }
  code {
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    font-size: 0.88em;
  }
  p code, li code, td code {
    background: var(--inline-code-bg);
    color: var(--inline-code-text);
    padding: 2px 6px;
    border-radius: 4px;
    font-weight: 600;
  }

  /* Cognitive Checkpoint Quizzes (Details / Summary) */
  details {
    background: var(--card-bg);
    border: 1px solid var(--border-color);
    border-radius: 10px;
    padding: 14px 18px;
    margin: 20px 0 24px;
    transition: background 0.2s;
  }
  details[open] {
    background: var(--card-bg-open);
    border-color: var(--accent);
  }
  summary {
    font-weight: 700;
    color: var(--accent);
    cursor: pointer;
    outline: none;
    list-style: none;
    display: flex;
    align-items: center;
    justify-content: space-between;
  }
  summary::-webkit-details-marker { display: none; }
  summary::after {
    content: "▼";
    font-size: 0.75rem;
    transition: transform 0.2s;
  }
  details[open] summary::after {
    transform: rotate(180deg);
  }
  details > *:not(summary) {
    margin-top: 14px;
  }

  /* Math Rendering Styles */
  .math-inline {
    font-family: "Cambria Math", "Latin Modern Math", "STIX Two Math", Georgia, serif;
    font-style: italic;
    color: var(--text-dark);
    padding: 0 3px;
  }
  .math-display {
    font-family: "Cambria Math", "Latin Modern Math", "STIX Two Math", Georgia, serif;
    background: var(--math-bg);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    padding: 14px 20px;
    margin: 16px 0 22px;
    text-align: center;
    overflow-x: auto;
    color: var(--text-dark);
    font-size: 1.05rem;
  }
  .math-frac {
    display: inline-flex;
    flex-direction: column;
    vertical-align: middle;
    text-align: center;
    padding: 0 4px;
    font-size: 0.92em;
    line-height: 1.15;
  }
  .math-frac .num {
    border-bottom: 1.5px solid currentColor;
    padding-bottom: 2px;
  }
  .math-frac .den {
    padding-top: 2px;
  }
  .math-text {
    font-style: normal;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  }
  .math-cases {
    display: inline-block;
    vertical-align: middle;
    font-style: normal;
  }
  .math-matrix {
    display: inline-table;
    vertical-align: middle;
    border-collapse: collapse;
    margin: 0 4px;
    border-left: 2px solid currentColor;
    border-right: 2px solid currentColor;
    border-radius: 4px;
    padding: 0 2px;
  }
  .math-matrix td {
    padding: 2px 6px;
    text-align: center;
    border: none;
    background: none;
  }

  /* Embedded Vector Diagrams */
  figure.embedded-diagram {
    margin: 32px 0;
    text-align: center;
  }
  figure.embedded-diagram .svg-wrapper {
    background: var(--bg-surface);
    border: 1px solid var(--border-color);
    border-radius: 12px;
    padding: 24px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.03);
    overflow-x: auto;
  }
  figure.embedded-diagram svg {
    max-width: 100%;
    height: auto;
    display: block;
    margin: 0 auto;
  }
  figure.embedded-diagram figcaption {
    margin-top: 12px;
    font-size: 0.86rem;
    color: #64748b;
  }

  /* Search Highlight */
  mark.search-highlight {
    background: var(--highlight-bg);
    color: inherit;
    padding: 2px 4px;
    border-radius: 3px;
  }

  /* Responsive layout (ENG-038/039): stack the sidebar; wide tables/code scroll inside the page */
  @media (max-width: 1024px) {
    .layout-wrapper {
      flex-direction: column;
    }
    aside.sidebar {
      width: 100%;
      min-width: 100%;
      height: auto;
      position: static;
      border-right: none;
      border-bottom: 1px solid var(--sidebar-border);
    }
    main.main-content {
      max-width: 100%;
      min-width: 0;
      padding: 24px 20px;
    }
    table:not(.math-matrix) {
      display: block;
      max-width: 100%;
      overflow-x: auto;
    }
    pre {
      max-width: 100%;
      overflow-x: auto;
    }
    .source-ref {
      display: inline-block;
      max-width: 100%;
      overflow-wrap: anywhere;
    }
    :not(pre) > code {
      overflow-wrap: anywhere;
    }
  }
  /* Print/PDF Mode Media Queries */
  @media print {
    body {
      background: #ffffff !important;
      color: #000000 !important;
    }
    aside.sidebar {
      display: none !important;
    }
    main.main-content {
      max-width: 100% !important;
      padding: 0 !important;
      margin: 0 !important;
    }
    .course-hero {
      background: none !important;
      border: 2px solid #000000 !important;
      padding: 24px !important;
      page-break-after: always;
    }
    article.chapter-container {
      page-break-after: always;
      border-bottom: none !important;
      margin-bottom: 0 !important;
      padding-bottom: 0 !important;
    }
    details {
      border: 1px solid #ccc !important;
      background: #fdfdfd !important;
      page-break-inside: avoid;
    }
    details summary::after {
      display: none !important;
    }
    figure.embedded-diagram {
      page-break-inside: avoid;
    }
    figure.embedded-diagram .svg-wrapper {
      border: 1px solid #ddd !important;
      background: #ffffff !important;
    }
    pre, table {
      page-break-inside: avoid;
    }
  }

/* ---- features carried over from other generations (Track 00/01 cover page + popup results, Track 10 details) ---- */
pre code { background-color: transparent; color: inherit; padding: 0; }
details .details-content, details > *:not(summary) { padding: 0 18px 16px; }
.sidebar-badge { background: var(--accent); border-radius: 9999px; color: #ffffff; display: inline-block; font-size: 0.75rem; font-weight: 700; letter-spacing: 0.05em; margin-bottom: 8px; padding: 3px 10px; text-transform: uppercase; }
.search-results { background: var(--bg-surface); border-radius: 8px; border: 1px solid var(--border-color); box-shadow: 0 10px 25px rgba(0,0,0,0.15); display: none; left: 0; margin-top: 6px; max-height: 300px; overflow-y: auto; position: absolute; right: 0; top: 100%; z-index: 100; }
.search-result-item { border-bottom: 1px solid var(--border-color); cursor: pointer; font-size: 0.82rem; padding: 10px 12px; }
.search-result-item:hover { background: var(--accent-light); }
.search-result-title { color: var(--accent); font-weight: 700; margin-bottom: 2px; }
.search-result-snippet { color: var(--text-main); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.action-btn.primary { background: var(--accent); border-color: var(--accent); color: #ffffff; }
.action-btn.primary:hover { background: var(--accent-hover); }
.document-container { margin: 0 auto; max-width: 900px; }
.cover-page { border-bottom: 2px solid var(--border-color); padding: 60px 20px 80px; page-break-after: always; text-align: center; }
.cover-badge { background: var(--accent); border-radius: 9999px; color: #ffffff; display: inline-block; font-size: 0.85rem; font-weight: 700; letter-spacing: 0.05em; margin-bottom: 24px; padding: 6px 16px; text-transform: uppercase; }
.cover-title { color: var(--text-dark); font-size: 2.6rem; font-weight: 800; line-height: 1.25; margin-bottom: 18px; }
.cover-subtitle { color: var(--text-main); font-size: 1.25rem; font-weight: 400; margin-bottom: 40px; }
.cover-meta { border-top: 1px solid var(--border-color); color: var(--text-main); font-size: 0.95rem; line-height: 1.8; margin: 0 auto; max-width: 650px; padding-top: 24px; text-align: left; }
.toc-section { background: var(--bg-surface); border-radius: 12px; border: 1px solid var(--border-color); margin: 40px 0; padding: 30px; page-break-after: always; }
.toc-title { border-bottom: 2px solid var(--accent); color: var(--text-dark); font-size: 1.4rem; font-weight: 700; margin-bottom: 20px; margin-top: 0; padding-bottom: 8px; }
.toc-list { list-style: none; margin: 0; padding: 0; }
.toc-item { font-size: 1rem; margin-bottom: 10px; }
.toc-item a { color: var(--text-dark); display: flex; font-weight: 600; justify-content: space-between; text-decoration: none; }
.toc-item a:hover { color: var(--accent); }
.print-btn { background: var(--accent); border-radius: 50px; border: none; bottom: 24px; box-shadow: 0 6px 16px rgba(2, 132, 199, 0.4); color: #ffffff; cursor: pointer; font-size: 0.95rem; font-weight: 700; padding: 12px 22px; position: fixed; right: 24px; transition: all 0.2s ease; z-index: 1000; }
.sidebar-header { flex-wrap: wrap; }
.sidebar-header .sidebar-badge { flex: 0 0 100%; border-radius: 6px; padding: 3px 8px; margin-bottom: 6px; text-align: left; }
.print-btn:hover { background: var(--accent-hover); transform: translateY(-2px); }
@media print { .no-print, .print-btn { display: none !important; } .document-container { max-width: none; padding: 0; } }
@media (max-width: 640px) { .print-btn { top: 8px; right: 8px; } }
'''

JS = r'''
  // Theme Toggle
  const themeToggle = document.getElementById('themeToggle');
  const currentTheme = localStorage.getItem('__THEME_KEY__') || 'light';
  document.documentElement.setAttribute('data-theme', currentTheme);

  themeToggle.addEventListener('click', () => {
    const activeTheme = document.documentElement.getAttribute('data-theme');
    const newTheme = activeTheme === 'light' ? 'dark' : 'light';
    document.documentElement.setAttribute('data-theme', newTheme);
    localStorage.setItem('__THEME_KEY__', newTheme);
  });

  // Collapsible Quiz Controls
  const expandBtn = document.getElementById('expandQuizzesBtn');
  const collapseBtn = document.getElementById('collapseQuizzesBtn');

  expandBtn.addEventListener('click', () => {
    document.querySelectorAll('details').forEach(d => d.open = true);
  });

  collapseBtn.addEventListener('click', () => {
    document.querySelectorAll('details').forEach(d => d.open = false);
  });

  // Print Action
  document.getElementById('printDocBtn').addEventListener('click', () => {
    document.querySelectorAll('details').forEach(d => d.open = true);
    window.print();
  });

  // Expand all details on beforeprint event
  window.addEventListener('beforeprint', () => {
    document.querySelectorAll('details').forEach(d => d.open = true);
  });

  // Real-Time Search & Highlight
  const searchInput = document.getElementById('searchInput');
  const searchCount = document.getElementById('searchCount');
  const mainContent = document.getElementById('mainContent');
  let originalArticlesHTML = new Map();

  document.querySelectorAll('article.chapter-container').forEach(art => {
    originalArticlesHTML.set(art.id, art.innerHTML);
  });

  searchInput.addEventListener('input', (e) => {
    const query = e.target.value.trim();
    if (!query) {
      searchCount.textContent = '';
      originalArticlesHTML.forEach((html, id) => {
        document.getElementById(id).innerHTML = html;
      });
      return;
    }

    if (query.length < 2) {
      searchCount.textContent = '2자 이상 입력...';
      return;
    }

    let totalMatches = 0;
    const esc = query.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    const regex = new RegExp(`(${esc})`, 'gi');
    const probe = new RegExp(esc, 'i');   // non-global: avoids the lastIndex bug of regex.test() with /g

    document.querySelectorAll('article.chapter-container').forEach(art => {
      const orig = originalArticlesHTML.get(art.id);
      art.innerHTML = orig;
      
      const walker = document.createTreeWalker(art, NodeFilter.SHOW_TEXT, null, false);
      const nodesToReplace = [];
      let node;
      while ((node = walker.nextNode())) {
        if (node.parentElement && ['SCRIPT', 'STYLE', 'SVG'].includes(node.parentElement.tagName)) continue;
        if (probe.test(node.nodeValue)) {
          nodesToReplace.push(node);
        }
      }

      nodesToReplace.forEach(textNode => {
        const matches = textNode.nodeValue.match(regex);
        if (matches) totalMatches += matches.length;
        // Build highlight nodes with DOM APIs (no innerHTML re-parse of decoded text)
        const text = textNode.nodeValue;
        const frag = document.createDocumentFragment();
        let last = 0, m;
        regex.lastIndex = 0;
        while ((m = regex.exec(text)) !== null) {
          if (m[0] === '') { regex.lastIndex++; continue; }
          frag.appendChild(document.createTextNode(text.slice(last, m.index)));
          const mark = document.createElement('mark');
          mark.className = 'search-highlight';
          mark.textContent = m[0];
          frag.appendChild(mark);
          last = m.index + m[0].length;
        }
        frag.appendChild(document.createTextNode(text.slice(last)));
        textNode.parentNode.replaceChild(frag, textNode);
      });
    });

    searchCount.textContent = `${totalMatches}개 발견됨`;
  });

  // Result list popup (ported from the Track 00/01 search): click a result to jump to the highlighted match
  const searchResults = document.getElementById('searchResults');
  function renderResultList() {
    if (!searchResults) return;
    const marks = Array.from(document.querySelectorAll('mark.search-highlight'));
    if (!marks.length) { searchResults.style.display = 'none'; searchResults.innerHTML = ''; return; }
    const seen = new Set(), items = [];
    for (const m of marks) {
      const host = m.closest('p, li, h2, h3, td') || m.parentElement;
      if (seen.has(host)) continue; seen.add(host); items.push([m, host]);
      if (items.length >= 10) break;
    }
    searchResults.innerHTML = '';
    items.forEach(([m, host]) => {
      const chap = host.closest('article.chapter-container');
      const badge = chap && chap.querySelector('.chapter-badge') ? chap.querySelector('.chapter-badge').textContent : '본문';
      const row = document.createElement('div'); row.className = 'search-result-item';
      const t = document.createElement('div'); t.className = 'search-result-title'; t.textContent = badge;
      const s = document.createElement('div'); s.className = 'search-result-snippet'; s.textContent = (host.textContent || '').trim().slice(0, 80);
      row.appendChild(t); row.appendChild(s);
      row.addEventListener('click', () => { m.scrollIntoView({ behavior: 'smooth', block: 'center' }); });
      searchResults.appendChild(row);
    });
    searchResults.style.display = 'block';
  }
  searchInput.addEventListener('input', () => setTimeout(renderResultList, 0));

'''

def render_page(title, sidebar_title, search_placeholder, theme_key, hero_html, sidebar_nav, chapters_html,
                sidebar_badge="", main_prefix="", main_suffix=""):
    badge = f'<div class="sidebar-badge">{sidebar_badge}</div>' if sidebar_badge else ""
    js = JS.replace("__THEME_KEY__", theme_key)
    return f"""<!DOCTYPE html>
<html lang="ko" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
{CSS}
</style>
</head>
<body>

<div class="layout-wrapper">
  <aside class="sidebar">
    <div class="sidebar-header">
      {badge}
      <h3 class="sidebar-title">{sidebar_title}</h3>
      <button class="theme-toggle-btn" id="themeToggle">🌓 테마</button>
    </div>

    <div class="search-container">
      <span class="search-icon">🔍</span>
      <input type="text" class="search-input" id="searchInput" placeholder="{search_placeholder}" autocomplete="off" />
      <div class="search-results-count" id="searchCount"></div>
      <div id="searchResults" class="search-results" style="display:none"></div>
    </div>

    <div class="global-actions">
      <button class="action-btn" id="expandQuizzesBtn">💡 퀴즈 펼치기</button>
      <button class="action-btn" id="collapseQuizzesBtn">🔒 퀴즈 접기</button>
      <button class="action-btn" id="printDocBtn">🖨️ PDF 인쇄</button>
    </div>

    <nav class="sidebar-nav">
      {sidebar_nav}
    </nav>
  </aside>

  <main class="main-content" id="mainContent">
    {main_prefix}
    {hero_html}
    {chapters_html}
    {main_suffix}
  </main>
</div>

<script>
{js}
</script>
</body>
</html>
"""

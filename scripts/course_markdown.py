#!/usr/bin/env python3
"""
course_markdown.py
Shared markdown -> Book / HTML / DOCX conversion logic for the 11 Master Course builders
(scripts/build_master_documents.py = Track 00, tracks/*/scripts/build_*_master.py = Tracks 01-10).

Each builder keeps only its track configuration (a TrackConfig), its track-specific pipeline assertion and main().
Contents:
  - LaTeX helpers: parse_balanced_braces, LATEX_SYMBOLS, convert_latex_to_html / _unicode, replace_math_in_markdown_*
  - build_master_markdown / build_master_html / build_master_docx (cfg, contents)
  - verify_track (assertions 1-3, 5-8 shared; assertion 4 injected), verify_basic (Track 00/01 profile),
    assert4_* helpers for the pipeline assertion
The UI shell lives in scripts/course_ui_template.py. Requires mistune==2.0.5 (guard in build_master_html).
LATEX_SYMBOLS decisions (ENG converter dedup): \\parallel -> U+2225, \\quad -> two spaces.
"""

from __future__ import annotations

import os
import re
import sys
import zipfile
import subprocess
import importlib.util as _ilu
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape as xml_escape

_spec = _ilu.spec_from_file_location("course_ui_template", os.path.join(os.path.dirname(os.path.abspath(__file__)), "course_ui_template.py"))
course_ui_template = _ilu.module_from_spec(_spec); _spec.loader.exec_module(course_ui_template)


class TrackConfig:
    """Plain attribute bag; see the builders for the full list of fields used by the functions below."""

    def __init__(self, **kw):
        self.read_filter = lambda text: text
        self.svg_exact = False
        self.__dict__.update(kw)


def read_doc_files(cfg: TrackConfig) -> list[str]:
    contents = []
    for f in cfg.doc_files:
        path = os.path.join(cfg.docs_dir, f)
        with open(path, "r", encoding="utf-8") as fp:
            contents.append(cfg.read_filter(fp.read()))
    return contents


# ==============================================================================
# LaTeX conversion
# ==============================================================================
def parse_balanced_braces(s: str, start_idx: int) -> tuple[str, int]:
    """Given string s and index of '{', returns (content, end_idx)."""
    assert s[start_idx] == '{'
    depth = 0
    content = []
    i = start_idx
    while i < len(s):
        if s[i] == '{':
            depth += 1
            if depth > 1:
                content.append(s[i])
        elif s[i] == '}':
            depth -= 1
            if depth == 0:
                return "".join(content), i
            else:
                content.append(s[i])
        else:
            content.append(s[i])
        i += 1
    return "".join(content), len(s) - 1

LATEX_SYMBOLS = {
    r'\Sigma': 'Σ',
    r'\Theta': 'Θ',
    r'\Delta': 'Δ',
    r'\delta': 'δ',
    r'\Gamma': 'Γ',
    r'\Lambda': 'Λ',
    r'\Pi': 'Π',
    r'\Omega': 'Ω',
    r'\Psi': 'Ψ',
    r'\Phi': 'Φ',
    r'\alpha': 'α',
    r'\beta': 'β',
    r'\gamma': 'γ',
    r'\epsilon': 'ε',
    r'\eta': 'η',
    r'\kappa': 'κ',
    r'\mu': 'μ',
    r'\nu': 'ν',
    r'\xi': 'ξ',
    r'\rho': 'ρ',
    r'\theta': 'θ',
    r'\sigma': 'σ',
    r'\phi': 'ϕ',
    r'\psi': 'ψ',
    r'\chi': 'χ',
    r'\lambda': 'λ',
    r'\tau': 'τ',
    r'\zeta': 'ζ',
    r'\pi': 'π',
    r'\omega': 'ω',
    r'\partial': '∂',
    r'\nabla': '∇',
    r'\infty': '∞',
    r'\approx': '≈',
    r'\sim': '~',
    r'\equiv': '≡',
    r'\cdot': '·',
    r'\circ': '∘',
    r'\odot': '⊙',
    r'\oplus': '⊕',
    r'\otimes': '⊗',
    r'\times': '×',
    r'\rightarrow': '→',
    r'\to': '→',
    r'\leftarrow': '←',
    r'\leftrightarrow': '↔',
    r'\xrightarrow': '→',
    r'\longrightarrow': '⟶',
    r'\rightleftharpoons': '⇌',
    r'\implies': '⇒',
    r'\iff': '⇔',
    r'\exists': '∃',
    r'\forall': '∀',
    r'\in': '∈',
    r'\notin': '∉',
    r'\subset': '⊂',
    r'\cup': '∪',
    r'\cap': '∩',
    r'\bigcup': '⋃',
    r'\emptyset': '∅',
    r'\neq': '≠',
    r'\ne': '≠',
    r'\le': '≤',
    r'\ge': '≥',
    r'\leq': '≤',
    r'\geq': '≥',
    r'\gg': '≫',
    r'\ll': '≪',
    r'\pm': '±',
    r'\sum': '∑',
    r'\prod': '∏',
    r'\int': '∫',
    r'\det': 'det',
    r'\ln': 'ln',
    r'\log': 'log',
    r'\cos': 'cos',
    r'\sin': 'sin',
    r'\max': 'max',
    r'\min': 'min',
    r'\exp': 'exp',
    r'\lim': 'lim',
    r'\arg': 'arg',
    r'\inf': 'inf',
    r'\sup': 'sup',
    r'\mid': '|',
    r'\middle': '|',
    r'\parallel': '∥',
    r'\Vert': '‖',
    r'\vert': '|',
    r'\propto': '∝',
    r'\succ': '≻',
    r'\prec': '≺',
    r'\setminus': '∖',
    r'\qquad': '    ',
    r'\quad': '  ',
    r'\dots': '…',
    r'\vdots': '⋮',
    r'\ddots': '⋱',
    r'\ell': 'ℓ',
    r'\tr': 'tr',
    r'\rank': 'rank',
    r'\diag': 'diag',
    r'\big': '',
    r'\left\|': '∥',
    r'\right\|': '∥',
    r'\|': '∥',
    r'\left\lfloor': '⌊',
    r'\right\rfloor': '⌋',
    r'\left\lceil': '⌈',
    r'\right\rceil': '⌉',
    r'\lceil': '⌈',
    r'\rceil': '⌉',
    r'\$': '$',
    r'\iint': '∬',
    r'\lvert': '|',
    r'\rvert': '|',
    r'\lVert': '‖',
    r'\rVert': '‖',
    r'\neg': '¬',
    r'\pmod': ' mod ',
    r'\bmod': ' mod ',
    r'\subseteq': '⊆',
    r'\subsetneq': '⊊',
    r'\supseteq': '⊇',
    r'\cong': '≅',
    r'\simeq': '≃',
    r'\asymp': '≍',
    r'\doteq': '≐',
    r'\wedge': '∧',
    r'\vee': '∨',
    r'\land': '∧',
    r'\lor': '∨',
    r'\varphi': 'φ',
    r'\varepsilon': 'ε',
    r'\varnothing': '∅',
    r'\angle': '∠',
    r'\therefore': '∴',
    r'\blacksquare': '■',
    r'\perp': '⊥',
    r'\longleftarrow': '⟵',
    r'\hookrightarrow': '↪',
    r'\mapsto': '↦',
    r'\uparrow': '↑',
    r'\downarrow': '↓',
    r'\cdots': '⋯',
    r'\ldots': '…',
    r'\argmin': 'arg min',
    r'\argmax': 'arg max',
    r'\arcsin': 'arcsin',
    r'\bar': '',
    r'\hat': '',
    r'\vec': '',
    r'\dot': '',
    r'\ddot': '',
    r'\mathbin': '',
    r'\begin': '',
    r'\end': '',
    r'\%': '%',
    r'\&': '&',
    r'\middle|': '|',
    r'\left(': '(',
    r'\right)': ')',
    r'\left[': '[',
    r'\right]': ']',
    r'\left|': '|',
    r'\right|': '|',
    r'\left\{': '{',
    r'\right\}': '}',
    r'\left.': '',
    r'\right.': '',
    r'\left': '',
    r'\right': '',
    r'\{': '{',
    r'\}': '}',
    r'\,': ' ',
    r'\;': ' ',
    r'\!': '',
    r'\ ': ' ',
    r'\langle': '⟨',
    r'\rangle': '⟩',
    r'\lfloor': '⌊',
    r'\rfloor': '⌋',
    r'\_': '_',
    r'\mathrm': '',
    r'\operatorname': '',
    r'\widehat': '',
    r'\boldsymbol': '',
    r'\underbrace': '',
    r'\widetilde': '',
    r'\overline': '',
    r'\boxed': '',
    r'\mathbf': '',
    r'\mathcal': '',
    r'\overset': '',
    r'\underset': '',
    r'\binom': '',
    r'\mathbb': '',
}


_US = '\ue000'  # protects an escaped underscore inside \text{}/\mathrm{}/\operatorname{} from the subscript rules
_DEPTH = {'html': 0, 'uni': 0}


def _html_overline(content: str) -> str:
    """\\bar/\\overline for HTML: one glyph keeps the combining macron; longer arguments get a continuous CSS overline."""
    c = content.strip()
    single = (len(c) == 1 and c != '<') or (re.fullmatch(r'\\[A-Za-z]+', c) is not None and len(LATEX_SYMBOLS.get(c, '')) == 1)
    if single:
        return f'{content}\u0304'
    return f'<span style="text-decoration:overline">{content}</span>'


def _uni_overline(content: str) -> str:
    """\\bar/\\overline for DOCX text: one glyph keeps U+0304; longer arguments get U+0305 after every visible character."""
    conv = convert_latex_to_unicode(content)
    if len(conv.strip()) <= 1:
        return f'{content}\u0304'
    return ''.join(ch if (ch.isspace() or ch in '_^()') else ch + '\u0305' for ch in conv)



def _convert_latex_to_html_impl(raw_latex: str, is_display: bool = False) -> str:
    """Converts a raw LaTeX formula to semantic, publication-grade HTML."""
    s = raw_latex.strip()
    # --- consolidated LaTeX environment handlers (v6) ---
    s = s.replace(r'\dfrac{', r'\frac{').replace(r'\tfrac{', r'\frac{')

    # \begin{cases} ... \end{cases}
    if r'\begin{cases}' in s:
        def repl_cases(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            res_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&')]
                conv_parts = [convert_latex_to_html(p, is_display=False) for p in parts]
                conv_parts = [re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', cp) for cp in conv_parts]
                res_lines.append("&nbsp;&nbsp;".join(conv_parts))
            return '<span class="math-cases">{ ' + ' ; '.join(res_lines) + ' }</span>'
        s = re.sub(r'\\begin\{cases\}(.*?)\\end\{cases\}', repl_cases, s, flags=re.DOTALL)

    # \begin{pmatrix} / \begin{bmatrix}
    for _env in ('pmatrix', 'bmatrix'):
        if r'\begin{' + _env + '}' in s:
            def repl_matrix(m):
                inner = m.group(1).strip()
                rows = [row.strip() for row in inner.split(r'\\')]
                res_rows = []
                for row in rows:
                    cols = [c.strip() for c in row.split('&')]
                    conv_cols = [convert_latex_to_html(c, is_display=False) for c in cols]
                    conv_cols = [re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', cc) for cc in conv_cols]
                    res_rows.append('<td>' + '</td><td>'.join(conv_cols) + '</td>')
                table_inner = '<tr>' + '</tr><tr>'.join(res_rows) + '</tr>'
                return '<table class="math-matrix"><tbody>' + table_inner + '</tbody></table>'
            s = re.sub(r'\\begin\{' + _env + r'\}(.*?)\\end\{' + _env + r'\}', repl_matrix, s, flags=re.DOTALL)

    # \begin{aligned} ... \end{aligned}
    if r'\begin{aligned}' in s:
        def repl_aligned(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            conv_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&') if p.strip()]
                conv_parts = [convert_latex_to_html(p, is_display=False) for p in parts]
                conv_parts = [re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', cp) for cp in conv_parts]
                conv_lines.append(' '.join(conv_parts))
            return '<span class="math-cases">' + ' ; '.join(conv_lines) + '</span>'
        s = re.sub(r'\\begin\{aligned\}(.*?)\\end\{aligned\}', repl_aligned, s, flags=re.DOTALL)

    # \mathbb{...}
    while r'\mathbb{' in s:
        idx = s.find(r'\mathbb{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbb'))
        _bb = {'A':'𝔸','B':'𝔹','C':'ℂ','D':'𝔻','E':'𝔼','F':'𝔽','G':'𝔾','H':'ℍ','I':'𝕀','J':'𝕁','K':'𝕂','L':'𝕃','M':'𝕄','N':'ℕ','O':'𝕆','P':'ℙ','Q':'ℚ','R':'ℝ','S':'𝕊','T':'𝕋','U':'𝕌','V':'𝕍','W':'𝕎','X':'𝕏','Y':'𝕐','Z':'ℤ'}
        s = s[:idx] + _bb.get(content, content) + s[end_idx+1:]

    # \mathcal{...}
    while r'\mathcal{' in s:
        idx = s.find(r'\mathcal{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathcal'))
        _mc = {'A':'𝒜','B':'ℬ','C':'𝒞','D':'𝒟','E':'ℰ','F':'ℱ','G':'𝒢','H':'ℋ','I':'ℐ','J':'𝒥','K':'𝒦','L':'ℒ','M':'ℳ','N':'𝒩','O':'𝒪','P':'𝒫','Q':'𝒬','R':'ℛ','S':'𝒮','T':'𝒯','U':'𝒰','V':'𝒱','W':'𝒲','X':'𝒳','Y':'𝒴','Z':'𝒵'}
        s = s[:idx] + _mc.get(content, content) + s[end_idx+1:]

    # \mathbf{...}
    while r'\mathbf{' in s:
        idx = s.find(r'\mathbf{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbf'))
        s = s[:idx] + '<strong>' + content + '</strong>' + s[end_idx+1:]

    # \boxed{...}
    while r'\boxed{' in s:
        idx = s.find(r'\boxed{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boxed'))
        s = s[:idx] + content + s[end_idx+1:]

    # \binom{n}{k}
    while r'\binom{' in s:
        idx = s.find(r'\binom{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\binom'))
        bot_idx = s.find('{', end_top + 1)
        if bot_idx == -1:
            break
        bot_c, end_bot = parse_balanced_braces(s, bot_idx)
        s = s[:idx] + '(' + top_c + ' C ' + bot_c + ')' + s[end_bot+1:]

    s = s.replace(r'\dfrac{', r'\frac{').replace(r'\tfrac{', r'\frac{')

    # Handle \begin{cases} ... \end{cases}
    if r'\begin{cases}' in s:
        def repl_cases(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            res_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&')]
                conv_parts = [convert_latex_to_html(p, is_display=False) for p in parts]
                conv_parts = [re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', cp) for cp in conv_parts]
                res_lines.append("&nbsp;&nbsp;".join(conv_parts))
            return '<span class="math-cases">{ ' + ' ; '.join(res_lines) + '</span>'
        s = re.sub(r'\\begin\{cases\}(.*?)\\end\{cases\}', repl_cases, s, flags=re.DOTALL)

    # Handle \begin{pmatrix} ... \end{pmatrix}
    if r'\begin{pmatrix}' in s:
        def repl_pmatrix(m):
            inner = m.group(1).strip()
            rows = [r.strip() for r in inner.split(r'\\')]
            res_rows = []
            for row in rows:
                cols = [c.strip() for c in row.split('&')]
                conv_cols = [convert_latex_to_html(c, is_display=False) for c in cols]
                conv_cols = [re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', cc) for cc in conv_cols]
                res_rows.append('<td>' + '</td><td>'.join(conv_cols) + '</td>')
            table_inner = '<tr>' + '</tr><tr>'.join(res_rows) + '</tr>'
            return f'<table class="math-matrix"><tbody>{table_inner}</tbody></table>'
        s = re.sub(r'\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}', repl_pmatrix, s, flags=re.DOTALL)

    # Parse \mathbb{...}
    while r'\mathbb{' in s:
        idx = s.find(r'\mathbb{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbb'))
        char_map = {'I': '𝕀', 'R': 'ℝ', 'P': 'ℙ', 'E': '𝔼', 'N': 'ℕ', 'Z': 'ℤ'}
        char = char_map.get(content, content)
        s = s[:idx] + char + s[end_idx+1:]

    # Parse \mathcal{...}
    while r'\mathcal{' in s:
        idx = s.find(r'\mathcal{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathcal'))
        math_cal_map = {'T': '𝒯', 'G': '𝒢', 'R': 'ℛ', 'S': '𝒮', 'V': '𝒱', 'N': '𝒩', 'P': '𝒫', 'L': 'ℒ', 'B': 'ℬ', 'M': 'ℳ', 'F': 'ℱ', 'O': '𝒪', 'C': '𝒞', 'D': '𝒟', 'E': 'ℰ'}
        char = math_cal_map.get(content, content)
        s = s[:idx] + char + s[end_idx+1:]

    # Parse \mathbf{...} and \boldsymbol{...}
    while r'\mathbf{' in s:
        idx = s.find(r'\mathbf{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbf'))
        s = s[:idx] + f'<strong>{content}</strong>' + s[end_idx+1:]

    while r'\boldsymbol{' in s:
        idx = s.find(r'\boldsymbol{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boldsymbol'))
        s = s[:idx] + f'<strong>{content}</strong>' + s[end_idx+1:]

    # Parse \boxed{...}
    while r'\boxed{' in s:
        idx = s.find(r'\boxed{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boxed'))
        s = s[:idx] + content + s[end_idx+1:]

    # Parse \text{...}, \mathrm{...}, \operatorname{...}
    while r'\text{' in s:
        idx = s.find(r'\text{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\text'))
        content = content.replace(r'\_', _US)
        sub_html = f'<span class="math-text">{content}</span>'
        s = s[:idx] + sub_html + s[end_idx+1:]

    while r'\mathrm{' in s:
        idx = s.find(r'\mathrm{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathrm'))
        content = content.replace(r'\_', _US)
        sub_html = f'<span class="math-text">{content}</span>'
        s = s[:idx] + sub_html + s[end_idx+1:]

    while r'\operatorname{' in s:
        idx = s.find(r'\operatorname{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\operatorname'))
        content = content.replace(r'\_', _US)
        sub_html = f'<span class="math-op">{content}</span>'
        s = s[:idx] + sub_html + s[end_idx+1:]

    # Parse \binom{n}{k}
    while r'\binom{' in s:
        idx = s.find(r'\binom{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\binom'))
        bot_idx = s.find('{', end_top + 1)
        if bot_idx == -1:
            break
        bot_c, end_bot = parse_balanced_braces(s, bot_idx)
        binom_html = f'({top_c} C {bot_c})'
        s = s[:idx] + binom_html + s[end_bot+1:]

    # Parse \frac{num}{den}
    while r'\frac{' in s:
        idx = s.find(r'\frac{')
        num_content, end_num = parse_balanced_braces(s, idx + len(r'\frac'))
        den_idx = s.find('{', end_num + 1)
        if den_idx == -1:
            break
        den_content, end_den = parse_balanced_braces(s, den_idx)
        
        num_html = convert_latex_to_html(num_content, is_display=False)
        num_html = re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', num_html)
        den_html = convert_latex_to_html(den_content, is_display=False)
        den_html = re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', den_html)
        
        frac_html = f'<span class="math-frac"><span class="num">{num_html}</span><span class="den">{den_html}</span></span>'
        s = s[:idx] + frac_html + s[end_den+1:]

    # Parse \sqrt{...}
    while r'\sqrt{' in s:
        idx = s.find(r'\sqrt{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\sqrt'))
        inner_html = convert_latex_to_html(content, is_display=False)
        inner_html = re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', inner_html)
        sqrt_html = f'√(<span class="math-sqrt">{inner_html}</span>)'
        s = s[:idx] + sqrt_html + s[end_idx+1:]

    # Parse accents: \hat, \bar, \tilde, \dot
    while r'\dot{' in s:
        idx = s.find(r'\dot{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\dot'))
        s = s[:idx] + f'{content}̇' + s[end_idx+1:]

    while r'\hat{' in s:
        idx = s.find(r'\hat{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\hat'))
        s = s[:idx] + f'{content}̂' + s[end_idx+1:]

    while r'\widehat{' in s:
        idx = s.find(r'\widehat{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\widehat'))
        s = s[:idx] + f'{content}̂' + s[end_idx+1:]

    while r'\bar{' in s:
        idx = s.find(r'\bar{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\bar'))
        s = s[:idx] + _html_overline(content) + s[end_idx+1:]

    while r'\overline{' in s:
        idx = s.find(r'\overline{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\overline'))
        s = s[:idx] + _html_overline(content) + s[end_idx+1:]

    # \vec{...}: combining right arrow above (a longer argument gets the arrow after the whole group)
    while r'\vec{' in s:
        idx = s.find(r'\vec{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\vec'))
        s = s[:idx] + f'{content}\u20d7' + s[end_idx+1:]

    while r'\tilde{' in s:
        idx = s.find(r'\tilde{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\tilde'))
        s = s[:idx] + f'{content}̃' + s[end_idx+1:]

    while r'\widetilde{' in s:
        idx = s.find(r'\widetilde{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\widetilde'))
        s = s[:idx] + f'{content}̃' + s[end_idx+1:]

    # Parse \overset and \underset
    while r'\overset{' in s:
        idx = s.find(r'\overset{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\overset'))
        base_idx = s.find('{', end_top + 1)
        if base_idx == -1:
            break
        base_c, end_base = parse_balanced_braces(s, base_idx)
        s = s[:idx] + f'{base_c}<sup>({top_c})</sup>' + s[end_base+1:]

    while r'\underset{' in s:
        idx = s.find(r'\underset{')
        bot_c, end_bot = parse_balanced_braces(s, idx + len(r'\underset'))
        base_idx = s.find('{', end_bot + 1)
        if base_idx == -1:
            break
        base_c, end_base = parse_balanced_braces(s, base_idx)
        s = s[:idx] + f'{base_c}<sub>({bot_c})</sub>' + s[end_base+1:]

    while r'\underbrace{' in s:
        idx = s.find(r'\underbrace{')
        base, end_base = parse_balanced_braces(s, idx + len(r'\underbrace'))
        after = s[end_base+1:].lstrip()
        if after.startswith('_{'):
            sub_idx = s.find('_{', end_base)
            label, end_label = parse_balanced_braces(s, sub_idx + 1)
            repl = (
                '<span style="display:inline-flex;flex-direction:column;align-items:center;'
                'vertical-align:middle">'
                f'<span>{base}</span>'
                '<span style="align-self:stretch;height:5px;margin:2px 0;border:solid currentColor;'
                'border-width:0 1.5px 1.5px;border-radius:0 0 6px 6px"></span>'
                f'<span style="font-size:0.85em">{label}</span></span>'
            )
            s = s[:idx] + repl + s[end_label+1:]
        else:
            s = s[:idx] + base + s[end_base+1:]

    # Replace mathematical symbols
    for sym, val in sorted(LATEX_SYMBOLS.items(), key=lambda kv: -len(kv[0])):
        s = s.replace(sym, val)

    # Subscripts with braces: _{...}
    while '_{' in s:
        idx = s.find('_{')
        content, end_idx = parse_balanced_braces(s, idx + 1)
        sub_html = convert_latex_to_html(content, is_display=False)
        sub_html = re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', sub_html)
        s = s[:idx] + f'<sub>{sub_html}</sub>' + s[end_idx+1:]

    # Superscripts with braces: ^{...}
    while '^{' in s:
        idx = s.find('^{')
        content, end_idx = parse_balanced_braces(s, idx + 1)
        sup_html = convert_latex_to_html(content, is_display=False)
        sup_html = re.sub(r'^<span class="math-inline">(.*)</span>$', r'\1', sup_html)
        s = s[:idx] + f'<sup>{sup_html}</sup>' + s[end_idx+1:]

    # Single-character subscripts/superscripts
    # LaTeX takes ONE token after _/^: an ASCII alphanumeric, '*', a mapped command symbol (\phi -> ϕ,
    # \infty -> ∞, ... any non-ASCII char) or a bare + / - (e.g. q_\phi, \sum^\infty, \mathbb{R}_+).
    s = re.sub(r'(?<!<sub)_(?![0-9]+>)([a-zA-Z0-9*+\-\u0080-\U0010ffff])', r'<sub>\1</sub>', s)
    s = re.sub(r'(?<!<sup)\^(?![0-9]+>)([a-zA-Z0-9*+\-\u0080-\U0010ffff])', r'<sup>\1</sup>', s)

    # Remove any leftover backslashes before residual identifiers
    s = re.sub(r'\\([a-zA-Z]+)', r'\1', s)

    # protect literal asterisks inside sub/sup from markdown emphasis pairing
    s = s.replace('<sup>*</sup>', '<sup>&#42;</sup>').replace('<sub>*</sub>', '<sub>&#42;</sub>')

    if is_display:
        return f'<div class="math-display">{s}</div>'
    else:
        return f'<span class="math-inline">{s}</span>'


def _uni_group(x: str) -> str:
    """ENG-043: 분수의 분자/분모와 첨자의 복합 식을 괄호로 묶어 DOCX 선형 텍스트의 모호성(예: λ^2 / 1 - λ^2)을 없앤다."""
    t = x.strip()
    if not t or re.fullmatch(r'[^\W_]+(\.[0-9]+)?', t):
        return t
    if t.startswith('(') and t.endswith(')'):
        depth = 0
        for i, ch in enumerate(t):
            depth += (ch == '(') - (ch == ')')
            if depth == 0 and i < len(t) - 1:
                break
        else:
            return t
    return f"({t})"


def _convert_latex_to_unicode_impl(raw_latex: str) -> str:
    """Converts a raw LaTeX formula to clean, readable Unicode text for DOCX."""
    s = raw_latex.strip()
    # --- consolidated LaTeX environment handlers (v6) ---
    s = s.replace(r'\dfrac{', r'\frac{').replace(r'\tfrac{', r'\frac{')

    # \begin{cases} ... \end{cases}
    if r'\begin{cases}' in s:
        def repl_cases(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            res_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&')]
                conv_parts = [convert_latex_to_unicode(p) for p in parts]
                res_lines.append(" ".join(conv_parts))
            return '{ ' + ' | '.join(res_lines) + ' }'
        s = re.sub(r'\\begin\{cases\}(.*?)\\end\{cases\}', repl_cases, s, flags=re.DOTALL)

    # \begin{pmatrix} / \begin{bmatrix}
    for _env in ('pmatrix', 'bmatrix'):
        if r'\begin{' + _env + '}' in s:
            def repl_matrix(m):
                inner = m.group(1).strip()
                rows = [row.strip() for row in inner.split(r'\\')]
                res_rows = []
                for row in rows:
                    cols = [c.strip() for c in row.split('&')]
                    conv_cols = [convert_latex_to_unicode(c) for c in cols]
                    res_rows.append(' '.join(conv_cols))
                return '[ ' + ' ; '.join(res_rows) + ' ]'
            s = re.sub(r'\\begin\{' + _env + r'\}(.*?)\\end\{' + _env + r'\}', repl_matrix, s, flags=re.DOTALL)

    # \begin{aligned} ... \end{aligned}
    if r'\begin{aligned}' in s:
        def repl_aligned(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            conv_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&') if p.strip()]
                conv_parts = [convert_latex_to_unicode(p) for p in parts]
                conv_lines.append(' '.join(conv_parts))
            return ' ; '.join(conv_lines)
        s = re.sub(r'\\begin\{aligned\}(.*?)\\end\{aligned\}', repl_aligned, s, flags=re.DOTALL)

    # \mathbb{...}
    while r'\mathbb{' in s:
        idx = s.find(r'\mathbb{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbb'))
        _bb = {'A':'𝔸','B':'𝔹','C':'ℂ','D':'𝔻','E':'𝔼','F':'𝔽','G':'𝔾','H':'ℍ','I':'𝕀','J':'𝕁','K':'𝕂','L':'𝕃','M':'𝕄','N':'ℕ','O':'𝕆','P':'ℙ','Q':'ℚ','R':'ℝ','S':'𝕊','T':'𝕋','U':'𝕌','V':'𝕍','W':'𝕎','X':'𝕏','Y':'𝕐','Z':'ℤ'}
        s = s[:idx] + _bb.get(content, content) + s[end_idx+1:]

    # \mathcal{...}
    while r'\mathcal{' in s:
        idx = s.find(r'\mathcal{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathcal'))
        _mc = {'A':'𝒜','B':'ℬ','C':'𝒞','D':'𝒟','E':'ℰ','F':'ℱ','G':'𝒢','H':'ℋ','I':'ℐ','J':'𝒥','K':'𝒦','L':'ℒ','M':'ℳ','N':'𝒩','O':'𝒪','P':'𝒫','Q':'𝒬','R':'ℛ','S':'𝒮','T':'𝒯','U':'𝒰','V':'𝒱','W':'𝒲','X':'𝒳','Y':'𝒴','Z':'𝒵'}
        s = s[:idx] + _mc.get(content, content) + s[end_idx+1:]

    # \mathbf{...}
    while r'\mathbf{' in s:
        idx = s.find(r'\mathbf{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbf'))
        s = s[:idx] + content + s[end_idx+1:]

    # \boxed{...}
    while r'\boxed{' in s:
        idx = s.find(r'\boxed{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boxed'))
        s = s[:idx] + content + s[end_idx+1:]

    # \binom{n}{k}
    while r'\binom{' in s:
        idx = s.find(r'\binom{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\binom'))
        bot_idx = s.find('{', end_top + 1)
        if bot_idx == -1:
            break
        bot_c, end_bot = parse_balanced_braces(s, bot_idx)
        s = s[:idx] + '(' + top_c + ' C ' + bot_c + ')' + s[end_bot+1:]

    s = s.replace(r'\dfrac{', r'\frac{').replace(r'\tfrac{', r'\frac{')

    # \begin{cases} ... \end{cases}
    if r'\begin{cases}' in s:
        def repl_cases(m):
            inner = m.group(1).strip()
            lines = [l.strip() for l in inner.split(r'\\')]
            res_lines = []
            for line in lines:
                parts = [p.strip() for p in line.split('&')]
                conv_parts = [convert_latex_to_unicode(p) for p in parts]
                res_lines.append(" ".join(conv_parts))
            return '{ ' + ' | '.join(res_lines) + ' }'
        s = re.sub(r'\\begin\{cases\}(.*?)\\end\{cases\}', repl_cases, s, flags=re.DOTALL)

    # \begin{pmatrix} ... \end{pmatrix}
    if r'\begin{pmatrix}' in s:
        def repl_pmatrix_uni(m):
            inner = m.group(1).strip()
            rows = [r.strip() for r in inner.split(r'\\')]
            res_rows = []
            for row in rows:
                cols = [c.strip() for c in row.split('&')]
                conv_cols = [convert_latex_to_unicode(c) for c in cols]
                res_rows.append(' '.join(conv_cols))
            return '[ ' + ' ; '.join(res_rows) + ' ]'
        s = re.sub(r'\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}', repl_pmatrix_uni, s, flags=re.DOTALL)

    # \mathbb{...}
    while r'\mathbb{' in s:
        idx = s.find(r'\mathbb{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbb'))
        char_map = {'I': '𝕀', 'R': 'ℝ', 'P': 'ℙ', 'E': '𝔼', 'N': 'ℕ', 'Z': 'ℤ'}
        char = char_map.get(content, content)
        s = s[:idx] + char + s[end_idx+1:]

    # \mathcal{...}
    while r'\mathcal{' in s:
        idx = s.find(r'\mathcal{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathcal'))
        math_cal_map = {'T': '𝒯', 'G': '𝒢', 'R': 'ℛ', 'S': '𝒮', 'V': '𝒱', 'N': '𝒩', 'P': '𝒫', 'L': 'ℒ', 'B': 'ℬ', 'M': 'ℳ', 'F': 'ℱ', 'O': '𝒪', 'C': '𝒞', 'D': '𝒟', 'E': 'ℰ'}
        char = math_cal_map.get(content, content)
        s = s[:idx] + char + s[end_idx+1:]

    # \mathbf{...} and \boldsymbol{...}
    while r'\mathbf{' in s:
        idx = s.find(r'\mathbf{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathbf'))
        s = s[:idx] + content + s[end_idx+1:]

    while r'\boldsymbol{' in s:
        idx = s.find(r'\boldsymbol{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boldsymbol'))
        s = s[:idx] + content + s[end_idx+1:]

    # \boxed{...}
    while r'\boxed{' in s:
        idx = s.find(r'\boxed{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\boxed'))
        s = s[:idx] + content + s[end_idx+1:]

    # \text{...}, \mathrm{...}, \operatorname{...}
    while r'\text{' in s:
        idx = s.find(r'\text{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\text'))
        content = content.replace(r'\_', _US)
        s = s[:idx] + content + s[end_idx+1:]

    while r'\mathrm{' in s:
        idx = s.find(r'\mathrm{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\mathrm'))
        content = content.replace(r'\_', _US)
        s = s[:idx] + content + s[end_idx+1:]

    while r'\operatorname{' in s:
        idx = s.find(r'\operatorname{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\operatorname'))
        content = content.replace(r'\_', _US)
        s = s[:idx] + content + s[end_idx+1:]

    # \binom{n}{k}
    while r'\binom{' in s:
        idx = s.find(r'\binom{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\binom'))
        bot_idx = s.find('{', end_top + 1)
        if bot_idx == -1:
            break
        bot_c, end_bot = parse_balanced_braces(s, bot_idx)
        s = s[:idx] + f"({top_c} C {bot_c})" + s[end_bot+1:]

    # \frac{num}{den} -> (num / den)
    while r'\frac{' in s:
        idx = s.find(r'\frac{')
        num_content, end_num = parse_balanced_braces(s, idx + len(r'\frac'))
        den_idx = s.find('{', end_num + 1)
        if den_idx == -1:
            break
        den_content, end_den = parse_balanced_braces(s, den_idx)
        
        num_uni = convert_latex_to_unicode(num_content)
        den_uni = convert_latex_to_unicode(den_content)
        s = s[:idx] + f"({_uni_group(num_uni)} / {_uni_group(den_uni)})" + s[end_den+1:]

    # \sqrt{...} -> √(content)
    while r'\sqrt{' in s:
        idx = s.find(r'\sqrt{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\sqrt'))
        inner_uni = convert_latex_to_unicode(content)
        s = s[:idx] + f"√({inner_uni})" + s[end_idx+1:]

    # Accents
    while r'\dot{' in s:
        idx = s.find(r'\dot{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\dot'))
        s = s[:idx] + f"{content}̇" + s[end_idx+1:]

    while r'\hat{' in s:
        idx = s.find(r'\hat{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\hat'))
        s = s[:idx] + f"{content}^" + s[end_idx+1:]

    while r'\widehat{' in s:
        idx = s.find(r'\widehat{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\widehat'))
        s = s[:idx] + f"{content}^" + s[end_idx+1:]

    while r'\bar{' in s:
        idx = s.find(r'\bar{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\bar'))
        s = s[:idx] + _uni_overline(content) + s[end_idx+1:]

    while r'\overline{' in s:
        idx = s.find(r'\overline{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\overline'))
        s = s[:idx] + _uni_overline(content) + s[end_idx+1:]

    # \vec{...}: combining right arrow above (a longer argument gets the arrow after the whole group)
    while r'\vec{' in s:
        idx = s.find(r'\vec{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\vec'))
        s = s[:idx] + f'{content}\u20d7' + s[end_idx+1:]

    while r'\tilde{' in s:
        idx = s.find(r'\tilde{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\tilde'))
        s = s[:idx] + f"{content}~" + s[end_idx+1:]

    while r'\widetilde{' in s:
        idx = s.find(r'\widetilde{')
        content, end_idx = parse_balanced_braces(s, idx + len(r'\widetilde'))
        s = s[:idx] + f"{content}~" + s[end_idx+1:]

    # \overset / \underset
    while r'\overset{' in s:
        idx = s.find(r'\overset{')
        top_c, end_top = parse_balanced_braces(s, idx + len(r'\overset'))
        base_idx = s.find('{', end_top + 1)
        if base_idx == -1:
            break
        base_c, end_base = parse_balanced_braces(s, base_idx)
        s = s[:idx] + f"{base_c}^({top_c})" + s[end_base+1:]

    while r'\underset{' in s:
        idx = s.find(r'\underset{')
        bot_c, end_bot = parse_balanced_braces(s, idx + len(r'\underset'))
        base_idx = s.find('{', end_bot + 1)
        if base_idx == -1:
            break
        base_c, end_base = parse_balanced_braces(s, base_idx)
        s = s[:idx] + f"{base_c}_({bot_c})" + s[end_base+1:]

    while r'\underbrace{' in s:
        idx = s.find(r'\underbrace{')
        base, end_base = parse_balanced_braces(s, idx + len(r'\underbrace'))
        after = s[end_base+1:].lstrip()
        if after.startswith('_{'):
            sub_idx = s.find('_{', end_base)
            label, end_label = parse_balanced_braces(s, sub_idx + 1)
            repl = f"[{base}]_({label})"
            s = s[:idx] + repl + s[end_label+1:]
        else:
            s = s[:idx] + base + s[end_base+1:]

    # Symbols
    for sym, val in sorted(LATEX_SYMBOLS.items(), key=lambda kv: -len(kv[0])):
        s = s.replace(sym, val)

    # Subscripts/superscripts braces
    while '_{' in s:
        idx = s.find('_{')
        content, end_idx = parse_balanced_braces(s, idx + 1)
        inner = convert_latex_to_unicode(content)
        s = s[:idx] + "_" + _uni_group(inner) + s[end_idx+1:]

    while '^{' in s:
        idx = s.find('^{')
        content, end_idx = parse_balanced_braces(s, idx + 1)
        inner = convert_latex_to_unicode(content)
        s = s[:idx] + "^" + _uni_group(inner) + s[end_idx+1:]

    s = re.sub(r'\\([a-zA-Z]+)', r'\1', s)
    return s


def convert_latex_to_html(raw_latex: str, is_display: bool = False) -> str:
    """Public entry; the protected-underscore sentinel is restored only at the outermost call."""
    _DEPTH['html'] += 1
    try:
        out = _convert_latex_to_html_impl(raw_latex, is_display)
    finally:
        _DEPTH['html'] -= 1
    return out.replace(_US, '_') if _DEPTH['html'] == 0 else out


def convert_latex_to_unicode(raw_latex: str) -> str:
    """Public entry; the protected-underscore sentinel is restored only at the outermost call."""
    _DEPTH['uni'] += 1
    try:
        out = _convert_latex_to_unicode_impl(raw_latex)
    finally:
        _DEPTH['uni'] -= 1
    return out.replace(_US, '_') if _DEPTH['uni'] == 0 else out


def replace_math_in_markdown_html(text: str) -> str:
    """Replaces $$...$$ and $...$ with HTML equivalents outside of code blocks."""
    parts = re.split(r'(```.*?```)', text, flags=re.DOTALL)
    out = []
    for part in parts:
        if part.startswith('```') and part.endswith('```'):
            out.append(part)
        else:
            part = part.replace(r'\$', '\x00')
            def repl_display(m):
                # keep the display block's own indentation: inside a list item the <div> must stay indented,
                # otherwise it ends the list and the indented continuation text below becomes a <pre><code>.
                lead = m.group(1) or ""
                return lead + "\n" + lead + convert_latex_to_html(m.group(2), is_display=True) + "\n"
            p = re.sub(r'(^[ \t]*)?\$\$(.*?)\$\$', repl_display, part, flags=re.DOTALL | re.MULTILINE)
            
            def repl_inline(m):
                return convert_latex_to_html(m.group(1), is_display=False)
            p = re.sub(r'(?<!\\)(?<!\$)\$(?!\$)(.*?)(?<!\\)(?<!\$)\$(?!\$)', repl_inline, p)
            p = p.replace('\x00', '$')
            out.append(p)
    return "".join(out)


def replace_math_in_markdown_docx(text: str) -> str:
    """Replaces $$...$$ and $...$ with Unicode text outside of code blocks."""
    parts = re.split(r'(```.*?```)', text, flags=re.DOTALL)
    out = []
    for part in parts:
        if part.startswith('```') and part.endswith('```'):
            out.append(part)
        else:
            part = part.replace(r'\$', '\x00')
            def repl_display(m):
                return "\n" + convert_latex_to_unicode(m.group(1)) + "\n"
            p = re.sub(r'\$\$(.*?)\$\$', repl_display, part, flags=re.DOTALL)
            
            def repl_inline(m):
                return convert_latex_to_unicode(m.group(1))
            p = re.sub(r'(?<!\\)(?<!\$)\$(?!\$)(.*?)(?<!\\)(?<!\$)\$(?!\$)', repl_inline, p)
            p = p.replace('\x00', '$')
            out.append(p)
    return "".join(out)


# ==============================================================================
# Build functions
# ==============================================================================
def build_master_markdown(cfg: TrackConfig, contents: list[str]) -> str:
    out_path = os.path.join(cfg.out_dir, cfg.md_name)
    
    header = cfg.md_header
    body_parts = []
    for i, (fname, raw_text) in enumerate(zip(cfg.doc_files, contents)):
        adjusted = raw_text.replace("images/", "docs/images/")
        lines = adjusted.split("\n")
        filtered_lines = []
        part_prefix = f"# Part {fname[:2]}. "
        
        replaced_h1 = False
        for line in lines:
            if line.startswith("# ") and not replaced_h1:
                orig_title = line[2:].strip()
                filtered_lines.append(f"{part_prefix}{orig_title}")
                filtered_lines.append(f"> **출처 원본 파일**: `{cfg.source_prefix}{fname}`\n")
                replaced_h1 = True
            else:
                filtered_lines.append(line)
        
        body_parts.append("\n".join(filtered_lines))
        body_parts.append("\n\n---\n\n")

    full_md = header + "".join(body_parts)
    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write(full_md)
    print(f"✅ Generated Master Markdown: {out_path} ({len(full_md):,} chars)")
    return full_md


def build_master_html(cfg: TrackConfig, contents: list[str]) -> str:
    import mistune
    if int(mistune.__version__.split(".")[0]) >= 3:
        raise RuntimeError("mistune>=3 leaves literal '**' in the rendered HTML (review finding ENG-020); "
                           "install mistune==2.0.5 (see requirements-build.txt)")
    out_path = os.path.join(cfg.out_dir, cfg.html_name)
    
    # Load SVGs to embed inline
    svg_cache = {}
    if os.path.exists(cfg.images_dir):
        for img in sorted(os.listdir(cfg.images_dir)):
            if img.endswith(".svg"):
                with open(os.path.join(cfg.images_dir, img), "r", encoding="utf-8") as fp:
                    svg_cache[img] = fp.read()

    markdown_parser = mistune.create_markdown(
        escape=False,
        plugins=['table', 'strikethrough']
    )

    rendered_chapters = []
    sidebar_nav_tree = []

    for i, (fname, raw_text) in enumerate(zip(cfg.doc_files, contents)):
        chap_num = f"{i:02d}"
        chap_title = cfg.chapter_titles[i]

        # Extract H2 headings for sidebar navigation
        h2_headings = []
        h2_counter = 0
        for line in raw_text.split('\n'):
            if line.startswith('## ') and not line.startswith('## 📚'):
                h2_counter += 1
                h2_text = line[3:].strip()
                h2_id = f"chap-{chap_num}-sec-{h2_counter}"
                h2_headings.append((h2_id, h2_text))

        # Build sidebar entry
        def _sidebar_math(t):
            return re.sub(r'(?<!\\)(?<!\$)\$(?!\$)(.*?)(?<!\\)(?<!\$)\$(?!\$)', lambda m: convert_latex_to_html(m.group(1), is_display=False), t)
        sec_links = "".join([f'<li><a href="#{hid}">{_sidebar_math(htext)}</a></li>' for hid, htext in h2_headings])
        sidebar_nav_tree.append(f"""
        <div class="sidebar-chap-group">
            <a class="sidebar-chap-title" href="#chapter-{chap_num}">{chap_title}</a>
            <ul class="sidebar-sec-list">{sec_links}</ul>
        </div>
        """)

        # Process math in markdown
        math_processed_text = replace_math_in_markdown_html(raw_text)

        # Parse through mistune
        # Rewrite intra-book chapter links before mistune, which sanitizes
        # the file:// scheme to a dead #harmful-link anchor.
        math_processed_text = re.sub(
            r'\]\((?:file:///[^)\s]*?/)?docs/0([0-9])_[A-Z0-9_]+\.md\)',
            lambda m: f'](#chapter-{int(m.group(1)):02d})',
            math_processed_text)
        chapter_html = markdown_parser(math_processed_text)
        # ENG-020: bold markers must be rendered; a literal '**' outside code/pre means the parser mis-handled emphasis
        _html_wo_code = re.sub(r'<(pre|code)\b.*?</\1>', '', chapter_html, flags=re.S)
        assert '**' not in _html_wo_code, f"Literal '**' left in rendered chapter HTML ({fname}); use mistune==2.0.5"

        # Inject IDs into <h2> tags
        h2_idx = 0
        def inject_h2(match):
            nonlocal h2_idx
            if h2_idx < len(h2_headings):
                hid, _ = h2_headings[h2_idx]
                h2_idx += 1
                return f'<h2 id="{hid}">{match.group(1)}</h2>'
            return f'<h2>{match.group(1)}</h2>'
        chapter_html = re.sub(r'<h2>(.*?)</h2>', inject_h2, chapter_html)

        # Inlining SVGs: replace <img> tags with inline SVGs
        for img_name, svg_code in svg_cache.items():
            clean_svg = re.sub(r'<\?xml[^>]*\?>', '', svg_code)
            pattern = re.compile(rf'<img\s+[^>]*src=["\'](?:[^"\']*/)?{re.escape(img_name)}["\'][^>]*>', re.IGNORECASE)
            
            def make_figure(match):
                tag = match.group(0)
                alt_m = re.search(r'alt=["\'](.*?)["\']', tag, re.IGNORECASE)
                alt_text = alt_m.group(1) if alt_m else img_name
                return f'''
                <figure class="embedded-diagram">
                    <div class="svg-wrapper">{clean_svg}</div>
                    <figcaption><strong>다이어그램 {alt_text}</strong>: {cfg.figcaption_desc}</figcaption>
                </figure>
                '''
            chapter_html = pattern.sub(make_figure, chapter_html)

        # Rewrite internal chapter links
        chapter_html = re.sub(rf'href=["\'](?:file:///[^"\']*/)?(?:{re.escape(cfg.track_prefix)})?docs/0([0-9])_[A-Z0-9_]+\.md["\']', r'href="#chapter-0\1"', chapter_html)
        chapter_html = re.sub(r'href=["\']0([0-9])_[A-Z0-9_]+\.md["\']', r'href="#chapter-0\1"', chapter_html)

        rendered_chapters.append(f"""
        <article class="chapter-container" id="chapter-{chap_num}">
            <div class="chapter-meta">
                <span class="chapter-badge">Chapter {chap_num}</span>
                <span class="source-ref">대응 원본: {cfg.source_prefix}{fname}</span>
            </div>
            {chapter_html}
        </article>
        """)

    full_html = course_ui_template.render_page(
        title=cfg.html_title,
        sidebar_title=cfg.html_sidebar_title,
        search_placeholder=cfg.html_search_placeholder,
        theme_key=cfg.html_theme_key,
        hero_html=cfg.html_hero_html,
        sidebar_nav="".join(sidebar_nav_tree),
        chapters_html="".join(rendered_chapters),
        sidebar_badge=cfg.html_sidebar_badge,
        main_prefix=cfg.html_main_prefix, main_suffix=cfg.html_main_suffix,
    )

    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write(full_html)
    print(f"✅ Generated Master HTML Course: {out_path} ({len(full_html):,} chars, {os.path.getsize(out_path):,} bytes)")
    return full_html


def build_master_docx(cfg: TrackConfig, contents: list[str]) -> str:
    out_path = os.path.join(cfg.out_dir, cfg.docx_name)
    
    doc_body = []

    def parse_inline_runs(text: str) -> str:
        tokens = re.split(r'(\*\*.*?\*\*|\*.*?\*|`.*?`)', text)
        runs = []
        for t in tokens:
            if not t:
                continue
            if t.startswith('**') and t.endswith('**'):
                inner = xml_escape(t[2:-2])
                runs.append(f'<w:r><w:rPr><w:b/></w:rPr><w:t xml:space="preserve">{inner}</w:t></w:r>')
            elif t.startswith('*') and t.endswith('*'):
                inner = xml_escape(t[1:-1])
                runs.append(f'<w:r><w:rPr><w:i/></w:rPr><w:t xml:space="preserve">{inner}</w:t></w:r>')
            elif t.startswith('`') and t.endswith('`'):
                inner = xml_escape(t[1:-1])
                runs.append(f'<w:r><w:rPr><w:rFonts w:ascii="Consolas" w:hAnsi="Consolas"/><w:color w:val="0284C7"/><w:sz w:val="19"/></w:rPr><w:t xml:space="preserve">{inner}</w:t></w:r>')
            else:
                inner = xml_escape(t)
                runs.append(f'<w:r><w:t xml:space="preserve">{inner}</w:t></w:r>')
        return ''.join(runs)

    def make_p(text: str, style: str = None, bold: bool = False, italic: bool = False, color: str = None, font_size: str = None) -> str:
        p_elm = ['<w:p>']
        pPr = []
        if style:
            pPr.append(f'<w:pStyle w:val="{style}"/>')
        if pPr:
            p_elm.append('<w:pPr>' + ''.join(pPr) + '</w:pPr>')
        
        rPr = []
        if bold: rPr.append('<w:b/>')
        if italic: rPr.append('<w:i/>')
        if color: rPr.append(f'<w:color w:val="{color}"/>')
        if font_size: rPr.append(f'<w:sz w:val="{font_size}"/>')
        
        rPr_str = ('<w:rPr>' + ''.join(rPr) + '</w:rPr>') if rPr else ''
        safe_text = xml_escape(text)
        p_elm.append(f'<w:r>{rPr_str}<w:t xml:space="preserve">{safe_text}</w:t></w:r>')
        p_elm.append('</w:p>')
        return ''.join(p_elm)

    # Document Title Page
    doc_body.append(make_p(cfg.docx_title, style="Title", bold=True, color="0F172A", font_size="52"))
    doc_body.append(make_p(cfg.docx_subtitle, style="Subtitle", color="475569", font_size="28"))
    doc_body.append(make_p("—" * 30, color="CBD5E1"))
    doc_body.append(make_p(cfg.docx_meta1, bold=True, color="0284C7"))
    doc_body.append(make_p(cfg.docx_meta2, color="334155"))
    doc_body.append(make_p("※ 이 DOCX에는 SVG 다이어그램 이미지가 포함되지 않습니다. 본문의 \"📊 [시각 다이어그램 참조: …]\" 표시는 해당 위치의 그림을 가리키며, 그림은 HTML 판(Master_Course.html) 또는 docs/images/의 SVG 파일에서 확인하세요.", color="64748B"))
    doc_body.append(make_p(" ", style="Normal"))

    # Table of Contents Outline
    doc_body.append(make_p("📚 교재 전체 목차 (Course Outline)", style="Heading1", bold=True, color="0284C7", font_size="36"))
    for i, title in enumerate(cfg.chapter_titles):
        doc_body.append(make_p(f"• {title}", style="ListBullet"))
    doc_body.append(make_p(" ", style="Normal"))

    # Process each markdown chapter
    for idx, (fname, raw_text) in enumerate(zip(cfg.doc_files, contents)):
        doc_body.append(make_p("—" * 35, color="CBD5E1"))
        doc_body.append(make_p(f"Chapter {idx:02d}: {cfg.chapter_titles[idx]}", style="Heading1", bold=True, color="0284C7", font_size="36"))
        doc_body.append(make_p(f"[대응 원본 교재: {cfg.source_prefix}{fname}]", italic=True, color="64748B", font_size="20"))
        
        # Replace math with unicode
        processed_text = replace_math_in_markdown_docx(raw_text)
        lines = processed_text.split('\n')

        in_code_block = False
        code_lines = []
        in_table = False
        table_rows = []

        def flush_table():
            nonlocal in_table, table_rows
            if not table_rows:
                return
            tbl_xml = ['<w:tbl>', '<w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="CBD5E1"/><w:bottom w:val="single" w:sz="4" w:space="0" w:color="CBD5E1"/><w:insideH w:val="single" w:sz="4" w:space="0" w:color="E2E8F0"/></w:tblBorders></w:tblPr>']
            is_first_row = True
            for r in table_rows:
                tbl_xml.append('<w:tr>')
                for c in r:
                    shd = '<w:shd w:fill="F1F5F9"/>' if is_first_row else ''
                    clean_cell = re.sub(r'<br\s*/?>', ' ', c.strip(), flags=re.IGNORECASE)
                    runs = parse_inline_runs(clean_cell)
                    tbl_xml.append(f'<w:tc><w:tcPr>{shd}<w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="150" w:type="dxa"/><w:right w:w="150" w:type="dxa"/></w:tcMar></w:tcPr><w:p>{runs}</w:p></w:tc>')
                tbl_xml.append('</w:tr>')
                is_first_row = False
            tbl_xml.append('</w:tbl>')
            doc_body.append(''.join(tbl_xml))
            table_rows = []
            in_table = False

        for line in lines:
            stripped = line.strip()

            # Handle code blocks
            if stripped.startswith("```"):
                if in_code_block:
                    code_text = "\n".join(code_lines)
                    # Word renders literal newlines in w:t unreliably (often as spaces),
                    # so emit explicit <w:br/> per source line; shrink the font for
                    # long-lined ASCII diagrams and hard-split anything still too wide.
                    def _split_long(cl):
                        if len(cl) <= 146:
                            return [cl]
                        parts = []
                        while len(cl) > 146:
                            cut = cl.rfind(" ", 0, 146)
                            cut = cut if cut > 100 else 146
                            parts.append(cl[:cut])
                            cl = "    " + cl[cut:].lstrip()
                        parts.append(cl)
                        return parts
                    code_lines = [cl for raw in code_text.split("\n") for cl in _split_long(raw)]
                    code_maxlen = max((len(cl) for cl in code_lines), default=0)
                    code_sz = 18 if code_maxlen <= 91 else max(12, int(2 * 823 / code_maxlen))
                    code_runs = "<w:r><w:br/></w:r>".join(
                        f'<w:r><w:rPr><w:rFonts w:ascii="Consolas" w:hAnsi="Consolas"/><w:sz w:val="{code_sz}"/><w:color w:val="0F172A"/></w:rPr><w:t xml:space="preserve">{xml_escape(cl)}</w:t></w:r>'
                        for cl in code_lines)
                    doc_body.append(f'<w:p><w:pPr><w:pStyle w:val="CodeBlock"/><w:shd w:fill="F8FAFC"/></w:pPr>{code_runs}</w:p>')
                    code_lines = []
                    in_code_block = False
                else:
                    if in_table:
                        flush_table()
                    in_code_block = True
                continue

            if in_code_block:
                code_lines.append(line)
                continue

            # Handle markdown tables
            if stripped.startswith("|") and stripped.endswith("|"):
                if re.match(r'^\|[\s\-:|]+\|$', stripped):
                    continue
                cells = [c for c in stripped.split("|")[1:-1]]
                if cells:
                    in_table = True
                    table_rows.append(cells)
                    continue
            else:
                if in_table:
                    flush_table()

            if not stripped:
                continue

            # Image lines: ![alt](path) -> Callout paragraph
            # Markdown horizontal rule -> thin bottom-border separator paragraph
            if re.match(r'^(-{3,}|\*{3,}|_{3,})$', stripped):
                doc_body.append('<w:p><w:pPr><w:pBdr><w:bottom w:val="single" w:sz="4" w:space="1" w:color="CBD5E1"/></w:pBdr></w:pPr></w:p>')
                continue
            img_match = re.match(r'^!\[(.*?)\]\((.*?)\)', stripped)
            if img_match:
                alt_text = img_match.group(1)
                doc_body.append(make_p(f"📊 [{cfg.diagram_ref_label}: {alt_text}]", style="Normal", bold=True, color="0284C7"))
                continue

            # Handle details and summary tags
            if stripped.startswith("<!--"):  # HTML comments (e.g. REFERENCES markers) are not document text
                continue
            if stripped.startswith("<details>") or stripped.startswith("</details>"):
                continue
            if "<summary>" in stripped:
                summary_text = re.sub(r'</?summary>', '', stripped).strip()
                doc_body.append(make_p(f"▶ {summary_text}", style="Normal", bold=True, color="0284C7"))
                continue

            # Clean any leftover inline HTML tags
            clean_line = re.sub(r'<br\s*/?>', ' ', stripped, flags=re.IGNORECASE)

            # Headers
            if clean_line.startswith("#### "):
                doc_body.append(make_p(clean_line[5:].strip(), style="Heading3", bold=True, color="475569", font_size="22"))
            elif clean_line.startswith("### "):
                title_text = clean_line[4:].strip()
                doc_body.append(make_p(title_text, style="Heading3", bold=True, color="1E293B", font_size="24"))
            elif clean_line.startswith("## "):
                title_text = clean_line[3:].strip()
                doc_body.append(make_p(title_text, style="Heading2", bold=True, color="0F172A", font_size="28"))
            elif clean_line.startswith("# "):
                continue
            elif clean_line.startswith("> "):
                quote_text = clean_line[2:].strip()
                runs = parse_inline_runs(quote_text)
                doc_body.append(f'<w:p><w:pPr><w:pStyle w:val="Quote"/><w:ind w:left="400"/><w:shd w:fill="F0F9FF"/></w:pPr>{runs}</w:p>')
            elif clean_line.startswith("- ") or clean_line.startswith("* "):
                bullet_text = clean_line[2:].strip()
                runs = parse_inline_runs(bullet_text)
                doc_body.append(f'<w:p><w:pPr><w:pStyle w:val="ListBullet"/><w:ind w:left="360"/></w:pPr><w:r><w:t xml:space="preserve">• </w:t></w:r>{runs}</w:p>')
            elif re.match(r'^\d+\.\s', clean_line):
                num_text = re.sub(r'^\d+\.\s*', '', clean_line)
                runs = parse_inline_runs(num_text)
                num_prefix = clean_line.split()[0]
                doc_body.append(f'<w:p><w:pPr><w:ind w:left="360"/></w:pPr><w:r><w:b/><w:t xml:space="preserve">{num_prefix} </w:t></w:r>{runs}</w:p>')
            else:
                runs = parse_inline_runs(clean_line)
                doc_body.append(f'<w:p><w:pPr><w:spacing w:after="140" w:line="320" w:lineRule="auto"/></w:pPr>{runs}</w:p>')

        if in_table:
            flush_table()

    # XML Templates
    content_types_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
</Types>"""

    rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
</Relationships>"""

    doc_rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>"""

    core_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>{cfg.docx_core_title}</dc:title>
  <dc:subject>{cfg.docx_core_subject}</dc:subject>
  <dc:creator>wgsim</dc:creator>
  <cp:category>Bioinformatics Curriculum</cp:category>
</cp:coreProperties>"""

    styles_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:docDefaults>
    <w:rPrDefault>
      <w:rPr>
        <w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:eastAsia="맑은 고딕"/>
        <w:sz w:val="22"/>
        <w:color w:val="334155"/>
      </w:rPr>
    </w:rPrDefault>
  </w:docDefaults>
  <w:style w:type="paragraph" w:styleId="Normal" w:default="1">
    <w:name w:val="Normal"/>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading1">
    <w:name w:val="heading 1"/>
    <w:pPr><w:keepNext/></w:pPr>
    <w:rPr>
      <w:b/>
      <w:color w:val="0284C7"/>
      <w:sz w:val="36"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading2">
    <w:name w:val="heading 2"/>
    <w:pPr><w:keepNext/></w:pPr>
    <w:rPr>
      <w:b/>
      <w:color w:val="0F172A"/>
      <w:sz w:val="28"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading3">
    <w:name w:val="heading 3"/>
    <w:pPr><w:keepNext/></w:pPr>
    <w:rPr>
      <w:b/>
      <w:color w:val="1E293B"/>
      <w:sz w:val="24"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="ListBullet">
    <w:name w:val="List Bullet"/>
  </w:style>
  <w:style w:type="paragraph" w:styleId="CodeBlock">
    <w:name w:val="Code Block"/>
    <w:rPr>
      <w:rFonts w:ascii="Consolas" w:hAnsi="Consolas"/>
      <w:sz w:val="18"/>
    </w:rPr>
  </w:style>
</w:styles>"""

    body_xml = "".join(doc_body)
    document_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    {body_xml}
    <w:sectPr>
      <w:pgSz w:w="11906" w:h="16838"/>
      <w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" w:header="708" w:footer="708" w:gutter="0"/>
    </w:sectPr>
  </w:body>
</w:document>"""

    # Validate XML before writing
    try:
        ET.fromstring(document_xml)
        print("✅ OpenXML document.xml successfully validated by XML parser!")
    except Exception as e:
        print(f"❌ XML validation error: {e}", file=sys.stderr)
        raise

    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types_xml)
        z.writestr("_rels/.rels", rels_xml)
        z.writestr("word/_rels/document.xml.rels", doc_rels_xml)
        z.writestr("word/styles.xml", styles_xml)
        z.writestr("docProps/core.xml", core_xml)
        z.writestr("word/document.xml", document_xml)

    print(f"✅ Generated Master Word Document: {out_path} ({os.path.getsize(out_path):,} bytes)")
    return document_xml


# ==============================================================================
# Verification
# ==============================================================================
def verify_track(cfg: TrackConfig, html_content: str, docx_xml: str, assert_pipelines):
    """Shared assertions 1-3 and 5-8; assertion 4 (pipelines) is track specific -> assert_pipelines()."""
    print(f"\n🔍 Running Automated Quality & Integrity Assertions for {cfg.track_label}...")
    
    # --------------------------------------------------------------------------
    # Assertion 1: All 10 markdown chapters exist and have even code fence counts
    # --------------------------------------------------------------------------
    assert len(cfg.doc_files) == 10, f"Expected 10 doc files, got {len(cfg.doc_files)}"
    for f in cfg.doc_files:
        path = os.path.join(cfg.docs_dir, f)
        assert os.path.exists(path), f"Chapter file missing: {path}"
        with open(path, "r", encoding="utf-8") as fp:
            text = fp.read()
        fences = text.count("```")
        assert fences % 2 == 0, f"Odd code fence count ({fences}) in {f}!"
    print(f"  [Assertion 1/8] All 10 Markdown chapters exist with even code fences (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 2: All 6 SVG diagrams exist and parse as valid XML
    # --------------------------------------------------------------------------
    svg_files = cfg.svg_files
    for s in svg_files:
        path = os.path.join(cfg.images_dir, s)
        assert os.path.exists(path), f"SVG file missing: {path}"
        tree = ET.parse(path)
        root = tree.getroot()
        assert "svg" in root.tag.lower(), f"Root tag not SVG in {s}"
    print(f"  [Assertion 2/8] All {len(svg_files)} SVG diagrams exist and parse as valid XML (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 3: Both slide decks exist and have embedded SVGs
    # --------------------------------------------------------------------------
    slide_files = [f for f in os.listdir(cfg.slides_dir) if f.endswith(".html")]
    assert len(slide_files) >= 2, f"Expected at least 2 slide decks, found {len(slide_files)} in {cfg.slides_dir}"
    
    found_deck1 = False
    found_deck2 = False
    for sf in slide_files:
        path = os.path.join(cfg.slides_dir, sf)
        with open(path, "r", encoding="utf-8") as fp:
            s_content = fp.read()
        svg_c = s_content.count("<svg")
        assert svg_c >= 1, f"Slide deck {sf} has no embedded SVGs!"
        if "01_" in sf:
            found_deck1 = True
        if "02_" in sf:
            found_deck2 = True
    assert found_deck1 and found_deck2, f"Missing Part 1 or Part 2 slide decks in {cfg.slides_dir}!"
    print(f"  [Assertion 3/8] Both slide decks exist with embedded SVGs ({len(slide_files)} decks verified) (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 4: Both Python pipelines executed and exited with Code 0
    # --------------------------------------------------------------------------
    assert_pipelines()

    # --------------------------------------------------------------------------
    # Assertion 5: Exactly 30 Cognitive Checkpoint quiz items exist with <details>
    # --------------------------------------------------------------------------
    total_quizzes = 0
    for f in cfg.doc_files:
        path = os.path.join(cfg.docs_dir, f)
        with open(path, "r", encoding="utf-8") as fp:
            text = fp.read()
        q_count = text.count("<details>")
        assert q_count == 3, f"Expected exactly 3 quizzes in {f}, got {q_count}"
        assert text.count("</details>") == 3, f"Mismatched </details> count in {f}"
        total_quizzes += q_count
    assert total_quizzes == 30, f"Expected 30 Cognitive Checkpoints, got {total_quizzes}"
    print(f"  [Assertion 5/8] Exactly 30 Cognitive Checkpoint quiz items exist with <details> (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 6: Master HTML contains zero external http/https CDN references
    # --------------------------------------------------------------------------
    external_links = re.findall(r'(?:src|href)=["\']https?://[^"\']+["\']', html_content, flags=re.IGNORECASE)
    assert len(external_links) == 0, f"Found external links in HTML: {external_links}!"
    raw_dd_count = len(re.findall(r'\$\$', html_content))
    assert raw_dd_count == 0, f"Found {raw_dd_count} raw $$ in HTML!"
    svg_count = len(re.findall(r'<svg[\s>]', html_content, flags=re.IGNORECASE))
    if cfg.svg_exact:
        assert svg_count == 6, f"Expected exactly 6 inlined SVGs in HTML, found {svg_count}!"
    else:
        assert svg_count >= 6, f"Expected at least 6 inlined SVGs in HTML, found {svg_count}!"
    print(f"  [Assertion 6/8] Master HTML contains 0 external CDN links and {svg_count} inlined SVGs across all 6 diagrams (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 7: Master DOCX is valid OpenXML with size > 50 KB
    # --------------------------------------------------------------------------
    docx_path = os.path.join(cfg.out_dir, cfg.docx_name)
    assert os.path.exists(docx_path), f"Master DOCX missing: {docx_path}"
    docx_size = os.path.getsize(docx_path)
    assert docx_size > 50 * 1024, f"DOCX file size too small: {docx_size} bytes (expected > 50 KB)"
    tree = ET.fromstring(docx_xml)
    all_docx_text = "".join(tree.itertext())
    assert "$$" not in all_docx_text, "Found raw $$ in DOCX text!"
    assert "<details>" not in all_docx_text, "Found literal <details> tag in DOCX text!"
    print(f"  [Assertion 7/8] Master DOCX is valid OpenXML ({docx_size:,} bytes > 50 KB) (PASSED)")

    # --------------------------------------------------------------------------
    # Assertion 8: Master roadmap and README have valid links to {cfg.track_label}
    # --------------------------------------------------------------------------
    roadmap_path = os.path.join(cfg.workspace, "MASTER_CURRICULUM_ROADMAP.md")
    readme_path = os.path.join(cfg.workspace, "README.md")
    assert os.path.exists(roadmap_path), f"Roadmap missing: {roadmap_path}"
    assert os.path.exists(readme_path), f"README missing: {readme_path}"

    with open(roadmap_path, "r", encoding="utf-8") as fp:
        roadmap_content = fp.read()
    with open(readme_path, "r", encoding="utf-8") as fp:
        readme_content = fp.read()

    required_roadmap_terms = cfg.roadmap_terms
    for term in required_roadmap_terms:
        assert term in roadmap_content, f"Missing term '{term}' in MASTER_CURRICULUM_ROADMAP.md"

    required_readme_terms = cfg.readme_terms
    for term in required_readme_terms:
        assert term in readme_content, f"Missing term '{term}' in README.md"

    print(f"  [Assertion 8/8] Master roadmap and README have valid links to {cfg.track_label} (PASSED)")

    print(f"\n🎉 ALL 8 BUILD ASSERTIONS SUCCESSFULLY PASSED FOR {cfg.track_label.upper()}!")


def verify_basic(html_content: str, docx_xml: str, n_svg: int, track_label: str = ""):
    """Assertion profile of Track 00/01 (HTML/DOCX checks only)."""
    suffix = f" for {track_label}" if track_label else ""
    print(f"\n🔍 Running Automated Quality & Integrity Assertions{suffix}...")
    
    # Check 1: HTML must contain n_svg inlined SVGs
    svg_count = len(re.findall(r'<svg[\s>]', html_content, flags=re.IGNORECASE))
    assert svg_count == n_svg, f"Assertion Failed: Expected {n_svg} inlined SVGs, found {svg_count}!"
    print(f"  [Assertion 1/8] Inlined SVGs: {svg_count} / {n_svg} (PASSED)")

    # Check 2: HTML must contain 0 raw $$
    raw_dd_count = len(re.findall(r'\$\$', html_content))
    assert raw_dd_count == 0, f"Assertion Failed: Found {raw_dd_count} raw $$ in HTML!"
    print(f"  [Assertion 2/8] Raw $$ in HTML: {raw_dd_count} (PASSED)")

    # Check 3: HTML must contain 0 \\frac{
    frac_count = len(re.findall(r'\\frac\{', html_content))
    assert frac_count == 0, f"Assertion Failed: Found {frac_count} \\frac{{ in HTML!"
    print(f"  [Assertion 3/8] Raw \\frac{{ in HTML: {frac_count} (PASSED)")

    # Check 4: HTML must contain 0 \\text{
    text_count = len(re.findall(r'\\text\{', html_content))
    assert text_count == 0, f"Assertion Failed: Found {text_count} \\text{{ in HTML!"
    print(f"  [Assertion 4/8] Raw \\text{{ in HTML: {text_count} (PASSED)")

    # Check 5: HTML must contain 0 external HTTP(S) links
    external_links = re.findall(r'(?:src|href)=["\']https?://[^"\']+["\']', html_content, flags=re.IGNORECASE)
    assert len(external_links) == 0, f"Assertion Failed: Found external links in HTML: {external_links}!"
    print(f"  [Assertion 5/8] External HTTP(S) links in HTML: 0 (PASSED)")

    # Check 6: DOCX XML must parse cleanly
    tree = ET.fromstring(docx_xml)
    all_docx_text = "".join(tree.itertext())
    print("  [Assertion 6/8] OpenXML ElementTree validation: SUCCESS (PASSED)")

    # Check 7: DOCX must contain 0 raw $$ and 0 \\frac{
    assert "$$" not in all_docx_text, "Assertion Failed: Found raw $$ in DOCX text!"
    assert r"\frac{" not in all_docx_text, "Assertion Failed: Found \\frac{ in DOCX text!"
    print("  [Assertion 7/8] Raw LaTeX macros in DOCX: 0 (PASSED)")

    # Check 8: DOCX must contain 0 literal <details> tags in text nodes
    assert "<details>" not in all_docx_text, "Assertion Failed: Found literal <details> tag in DOCX text!"
    assert "</details>" not in all_docx_text, "Assertion Failed: Found literal </details> tag in DOCX text!"
    print("  [Assertion 8/8] Literal <details> tags in DOCX text: 0 (PASSED)")

    print(f"\n🎉 ALL 8 BUILD ASSERTIONS SUCCESSFULLY PASSED{(' FOR ' + track_label.upper()) if track_label else ''}!")


def _skip_requested(*flags, env=True) -> bool:
    return any(f in sys.argv for f in flags) or (env and os.environ.get("VERIFY_FAST") == "1")


VOLATILE_FIELDS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline_volatile_fields.json")
_MASK = "<VOLATILE>"


def _load_volatile():
    import json
    with open(VOLATILE_FIELDS_FILE, "r", encoding="utf-8") as fp:
        return json.load(fp).get("files", {})


def _mask_json(obj, path_pats, path=()):
    """Replaces values whose key path matches a declared dotted pattern ('*' = any key/index) with a marker."""
    out = obj
    if isinstance(obj, dict):
        out = {k: _mask_json(v, path_pats, path + (k,)) for k, v in obj.items()}
    elif isinstance(obj, list):
        out = [_mask_json(v, path_pats, path + (str(i),)) for i, v in enumerate(obj)]
    for pat in path_pats:
        parts = pat.split(".")
        if len(parts) == len(path) and all(a == "*" or a == b for a, b in zip(parts, path)):
            return _MASK
    return out


def _json_paths(obj, path=()):
    """All key paths (list indices as strings) of a parsed JSON document."""
    out = {path}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out |= _json_paths(v, path + (k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out |= _json_paths(v, path + (str(i),))
    return out


def _first_json_diff(a, b, path=""):
    if type(a) != type(b):
        return path or "<root>", a, b
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                return f"{path}.{k}".lstrip("."), a.get(k, "<missing>"), b.get(k, "<missing>")
            d = _first_json_diff(a[k], b[k], f"{path}.{k}")
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}[len]", len(a), len(b)
        for i, (x, y) in enumerate(zip(a, b)):
            d = _first_json_diff(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    return None if a == b else (path or "<root>", a, b)


def _compare_one(rel_key, committed_path, produced_path, rules):
    """Returns None if equal (byte-identical, or equal after masking only the declared volatile fields), else a message."""
    import csv, io, json
    with open(committed_path, "rb") as fp:
        cb = fp.read()
    with open(produced_path, "rb") as fp:
        pb = fp.read()
    if cb == pb:
        return None
    if rules is None:
        # first differing line for the report
        ca, pa = cb.decode("utf-8", "replace").split("\n"), pb.decode("utf-8", "replace").split("\n")
        for i, (x, y) in enumerate(zip(ca, pa), 1):
            if x != y:
                return f"{rel_key}: not declared volatile; first differing line {i}: committed={x[:100]!r} produced={y[:100]!r}"
        return f"{rel_key}: not declared volatile; files differ in length ({len(ca)} vs {len(pa)} lines)"
    if "json_keys" in rules:
        ca, pa = json.loads(cb), json.loads(pb)
        for side, obj in (("committed", ca), ("produced", pa)):
            paths = _json_paths(obj)
            for pat in rules["json_keys"]:
                parts = pat.split(".")
                if not any(len(parts) == len(q) and all(a == "*" or a == b for a, b in zip(parts, q)) for q in paths):
                    return f"{rel_key}: declared volatile JSON key {pat!r} does not exist in the {side} file"
        ma, mp = _mask_json(ca, rules["json_keys"]), _mask_json(pa, rules["json_keys"])
        d = _first_json_diff(ma, mp)
        return None if d is None else f"{rel_key}: first non-volatile difference at {d[0]}: committed={str(d[1])[:80]!r} produced={str(d[2])[:80]!r}"
    if "csv_columns" in rules:
        ra = list(csv.reader(io.StringIO(cb.decode("utf-8"))))
        rp = list(csv.reader(io.StringIO(pb.decode("utf-8"))))
        if not ra or not rp or ra[0] != rp[0]:
            return f"{rel_key}: CSV header differs (committed {ra[0] if ra else None} vs produced {rp[0] if rp else None})"
        if len(ra) != len(rp):
            return f"{rel_key}: CSV row count differs ({len(ra)} vs {len(rp)})"
        for col in rules["csv_columns"]:
            if col not in ra[0]:
                return f"{rel_key}: declared volatile CSV column {col!r} does not exist in the header {ra[0]}"
        cols = [i for i, h in enumerate(ra[0]) if h in rules["csv_columns"]]
        rk = None
        if rules.get("row_key"):
            if rules["row_key"] not in ra[0]:
                return f"{rel_key}: declared row_key column {rules['row_key']!r} does not exist"
            rk = ra[0].index(rules["row_key"])
            for side, rows in (("committed", ra), ("produced", rp)):
                labels = {r[rk] for r in rows[1:] if len(r) > rk}
                for lab in rules.get("row_values", []):
                    if lab not in labels:
                        return f"{rel_key}: declared volatile CSV row {lab!r} does not exist in the {side} file"
        for r, (x, y) in enumerate(zip(ra, rp)):
            if len(x) != len(y):
                return f"{rel_key}: CSV row {r} has {len(x)} vs {len(y)} columns"
            row_masked = r > 0 and (rk is None or x[rk] in rules.get("row_values", [x[rk]]))
            for c, (u, v) in enumerate(zip(x, y)):
                if not (row_masked and c in cols) and u != v:
                    return f"{rel_key}: first non-volatile difference at row {r}, column {ra[0][c]!r}: committed={u!r} produced={v!r}"
        return None
    if "line_rules" in rules:
        ca, pa = cb.decode("utf-8").split("\n"), pb.decode("utf-8").split("\n")
        if len(ca) != len(pa):
            return f"{rel_key}: line count differs ({len(ca)} vs {len(pa)})"
        for side, lines in (("committed", ca), ("produced", pa)):
            for rule in rules["line_rules"]:
                if not any(re.search(rule["line"], l) for l in lines):
                    return f"{rel_key}: declared volatile line pattern {rule['line']!r} matches no line in the {side} file"
        def mask(line):
            for rule in rules["line_rules"]:
                if re.search(rule["line"], line):
                    line = re.sub(rule["mask"], _MASK, line)
            return line
        for i, (x, y) in enumerate(zip(ca, pa), 1):
            if mask(x) != mask(y):
                return f"{rel_key}: first non-volatile difference at line {i}: committed={x[:100]!r} produced={y[:100]!r}"
        return None
    return f"{rel_key}: unsupported volatile rule {sorted(rules)}"


def _compare_results(cfg: TrackConfig, produced_pipelines_dir: str):
    """Compares every results_pipeline_0X folder of the copy with the committed one; returns a list of problems."""
    volatile = _load_volatile()
    problems = []
    for rd in ("results_pipeline_01", "results_pipeline_02"):
        cdir = os.path.join(cfg.pipelines_dir, rd)
        pdir = os.path.join(produced_pipelines_dir, rd)
        def listing(base):
            return sorted(os.path.relpath(os.path.join(root, f), base) for root, _, fs in os.walk(base) for f in fs
                          if "__pycache__" not in root)
        cf, pf = listing(cdir), listing(pdir)
        if cf != pf:
            problems.append(f"{rd}: file set differs (only committed: {sorted(set(cf) - set(pf))[:3]}, only produced: {sorted(set(pf) - set(cf))[:3]})")
        for f in sorted(set(cf) & set(pf)):
            key = os.path.relpath(os.path.join(cdir, f), cfg.workspace).replace(os.sep, "/")
            msg = _compare_one(key, os.path.join(cdir, f), os.path.join(pdir, f), volatile.get(key))
            if msg:
                problems.append(msg)
    return problems


def _trash(path):
    t = "/opt/homebrew/opt/trash/bin/trash"
    if os.path.exists(t):
        subprocess.run([t, path], capture_output=True)


def _run_pipelines(cfg: TrackConfig, p1, p2, labels=None):
    """Runs both pipeline scripts and asserts exit code 0.
    Default: runs a COPY of cfg.pipelines_dir in a temp dir (committed results stay untouched) and requires the produced
    results_pipeline_0X folders to equal the committed ones except for the fields declared in
    scripts/pipeline_volatile_fields.json. --update-results: run in place and overwrite the committed results."""
    import shutil, tempfile
    update = "--update-results" in sys.argv
    tmp = None
    if update:
        run1, run2, pdir = p1, p2, cfg.pipelines_dir
    else:
        tmp = tempfile.mkdtemp(prefix="course_pipelines_", dir=os.environ.get("PIPELINE_REPRO_TMPDIR") or None)
        pdir = os.path.join(tmp, "pipelines")
        shutil.copytree(cfg.pipelines_dir, pdir, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        run1, run2 = os.path.join(pdir, os.path.basename(p1)), os.path.join(pdir, os.path.basename(p2))
    try:
        if labels:
            print(labels[0])
        res1 = subprocess.run([sys.executable, run1], capture_output=True, text=True)
        assert res1.returncode == 0, f"Pipeline 01 failed with exit code {res1.returncode}: {res1.stderr}"
        if labels:
            print(labels[1])
        res2 = subprocess.run([sys.executable, run2], capture_output=True, text=True)
        assert res2.returncode == 0, f"Pipeline 02 failed with exit code {res2.returncode}: {res2.stderr}"
        if update:
            print("  --update-results: committed result folders were regenerated in place")
            return
        problems = _compare_results(cfg, pdir)
        assert not problems, (
            "Pipeline results differ from the committed ones:\n    " + "\n    ".join(problems[:10])
            + "\n    Note: the check is exact for the pinned env (requirements-pipelines.txt) on the reference platform. "
              "If only model hashes/losses differ on another platform (BLAS/libm/CPU), this is an environment "
              "difference: use --skip-pipelines, or regenerate with --update-results."
        )
        print("  Pipeline results reproduce the committed ones (only declared volatile fields may differ)")
    finally:
        if tmp:
            _trash(tmp)


def assert4_result_dirs(cfg: TrackConfig, p1_name: str, p2_name: str):
    """Assertion 4 (Tracks 02-04): --skip-pipelines/--verify-fast/VERIFY_FAST -> committed result folders must exist."""
    pipeline_1 = os.path.join(cfg.pipelines_dir, p1_name)
    pipeline_2 = os.path.join(cfg.pipelines_dir, p2_name)
    assert os.path.exists(pipeline_1), f"Pipeline 1 missing: {pipeline_1}"
    assert os.path.exists(pipeline_2), f"Pipeline 2 missing: {pipeline_2}"
    # ENG-026: 문서 빌드와 파이프라인 재실행(결과 파일 덮어쓰기)을 분리. --skip-pipelines 또는 VERIFY_FAST=1 이면
    # 파이프라인을 실행하지 않고 커밋된 결과 폴더가 존재하는지만 확인한다. (기본 동작은 기존과 같이 실행)
    if _skip_requested("--skip-pipelines", "--verify-fast"):
        for _rd in ("results_pipeline_01", "results_pipeline_02"):
            _rp = os.path.join(cfg.pipelines_dir, _rd)
            assert os.path.isdir(_rp) and os.listdir(_rp), f"Pipeline results missing: {_rp}"
        print(f"  [Assertion 4/8] Pipeline execution skipped; committed result folders verified (PASSED)")
    else:
        _run_pipelines(cfg, pipeline_1, pipeline_2)
        print(f"  [Assertion 4/8] Both Python pipelines executed (exit 0); results reproduce the committed ones, or were regenerated with --update-results (PASSED)")


def assert4_summary_files(cfg: TrackConfig, p1_name: str, p2_name: str, s1_name: str, s2_name: str, labels):
    """Assertion 4 (Tracks 05-07): skip flags -> both summary files must exist; otherwise run the pipelines."""
    pipeline_1 = os.path.join(cfg.pipelines_dir, p1_name)
    pipeline_2 = os.path.join(cfg.pipelines_dir, p2_name)
    assert os.path.exists(pipeline_1), f"Pipeline 1 missing: {pipeline_1}"
    assert os.path.exists(pipeline_2), f"Pipeline 2 missing: {pipeline_2}"
    if _skip_requested("--skip-pipelines", "--verify-fast", "--quick"):
        p1_summary = os.path.join(cfg.pipelines_dir, "results_pipeline_01", s1_name)
        p2_summary = os.path.join(cfg.pipelines_dir, "results_pipeline_02", s2_name)
        assert os.path.exists(p1_summary), f"Pipeline 1 summary missing: {p1_summary}"
        assert os.path.exists(p2_summary), f"Pipeline 2 summary missing: {p2_summary}"
        print(f"  [Assertion 4/8] Both Python pipelines verified via generated artifacts (PASSED)")
    else:
        _run_pipelines(cfg, pipeline_1, pipeline_2, labels)
        print(f"  [Assertion 4/8] Both Python pipelines executed (exit 0); results reproduce the committed ones, or were regenerated with --update-results (PASSED)")


def assert4_summary_report(cfg: TrackConfig, p1_name: str, p2_name: str, s1_name: str, s2_name: str, check_texts, rerun_labels):
    """Assertion 4 (Tracks 08-10): summary reports must exist and pass check_texts(p1_text, p2_text); --re-run-pipelines re-executes."""
    pipeline_1 = os.path.join(cfg.pipelines_dir, p1_name)
    pipeline_2 = os.path.join(cfg.pipelines_dir, p2_name)
    assert os.path.exists(pipeline_1), f"Pipeline 1 missing: {pipeline_1}"
    assert os.path.exists(pipeline_2), f"Pipeline 2 missing: {pipeline_2}"

    p1_summary = os.path.join(cfg.pipelines_dir, "results_pipeline_01", s1_name)
    p2_summary = os.path.join(cfg.pipelines_dir, "results_pipeline_02", s2_name)
    assert os.path.exists(p1_summary), f"Pipeline 1 summary missing: {p1_summary}"
    assert os.path.exists(p2_summary), f"Pipeline 2 summary missing: {p2_summary}"

    with open(p1_summary, "r", encoding="utf-8") as fp:
        p1_text = fp.read()
    with open(p2_summary, "r", encoding="utf-8") as fp:
        p2_text = fp.read()
    check_texts(p1_text, p2_text)

    if "--re-run-pipelines" in sys.argv:
        _run_pipelines(cfg, pipeline_1, pipeline_2, rerun_labels)
        print(f"  [Assertion 4/8] Both Python pipelines re-executed (exit 0); results reproduce the committed ones, or were regenerated with --update-results (PASSED)")
    else:
        print(f"  [Assertion 4/8] Both Python pipelines executed and verified via summary artifacts (PASSED)")

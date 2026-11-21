#!/usr/bin/env python3
"""
precommit_checks.py
Pre-commit checks for the course sources (Python standard library only).

Each check guards a problem class found in the 2026-10 review (see review/REVIEW_REPORT.md):
  BLOCK  source .md: control characters, unbalanced $$, python code blocks that do not parse,
         file:/// links, relative links whose target does not exist
  BLOCK  .py files that do not parse
  BLOCK  files larger than 5 MB outside results_pipeline_* (snapshots/tarballs belong outside git)
  BLOCK  results_pipeline_* changes unless ALLOW_RESULTS_CHANGE=1
         (results must be regenerated in the scRNAseq_anal env, see requirements-pipelines.txt)
  WARN   source chapter staged without its generated *_Master_Book.md (rebuild needed)
  WARN   generated Book/HTML/DOCX staged without any source chapter (possible hand edit);
         silent when that track's builder or a shared build module (course_markdown/course_ui_template) is staged too
  BLOCK  chapter reference block (<!-- REFERENCES:BEGIN -->) that differs from review/refs/*.ids + cache (hand edit)
  WARN   source chapter without a reference block (rollout period)
  WARN   slide deck HTML (slides/*.html) without the SLIDE-LAYOUT:BEGIN marker (run scripts/inject_slide_layout.py)
  BLOCK  text files that contain a local machine path (/Users/<name>/, a session scratchpad, a conda Caskroom path) or a
         personal mail address (the repository is public; write <repo>, <scratchpad>, <conda-base> instead)

Usage:
  python3 scripts/precommit_checks.py            # staged files (used by .githooks/pre-commit)
  python3 scripts/precommit_checks.py --all      # every tracked file
"""

import ast
import os
import re
import subprocess
import sys
import textwrap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAX_BYTES = 5 * 1024 * 1024
SOURCE_MD = re.compile(r"^(docs/[^/]+\.md|tracks/[^/]+/docs/[^/]+\.md|README\.md|MASTER_CURRICULUM_ROADMAP\.md)$")
GENERATED = re.compile(r"_Master_(Book\.md|Course\.html|Course\.docx)$")
SHARED_BUILD = {"scripts/course_markdown.py", "scripts/course_ui_template.py"}
TRACK_BUILDER = re.compile(r"^(scripts/build_master_documents\.py|tracks/[^/]+/scripts/build_[^/]+_master\.py)$")
CTRL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
PY_BLOCK = re.compile(r"```python\n(.*?)```", re.S)
MD_LINK = re.compile(r"\]\(([^)\s#]+)(?:#[^)]*)?\)")
# Local paths / personal addresses must not reach the public repository (patterns do not match their own source text).
PRIVATE = re.compile(r"/Users/(?!<)[A-Za-z0-9_.-]+/|-Users-[A-Za-z0-9_.-]+-|/private/tmp/claude-\d|Caskroom/(?:miniconda|anaconda)"
                     r"|[A-Za-z0-9._%+-]+@(?:gmail|naver|daum|hanmail|kakao|outlook|hotmail|icloud|yahoo)\.(?:com|net)")


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def target_files(all_files):
    if all_files:
        return [f for f in git("ls-files").splitlines() if f]
    return [f for f in git("diff", "--cached", "--name-only", "--diff-filter=ACMR").splitlines() if f]


def read_text(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return fh.read()


def check_markdown(path, errors):
    text = read_text(path)
    if CTRL.search(text):
        errors.append(f"{path}: control character found")
    if text.count("$$") % 2:
        errors.append(f"{path}: unbalanced $$ delimiters")
    for i, block in enumerate(PY_BLOCK.findall(text), 1):
        try:
            ast.parse(textwrap.dedent(block))
        except SyntaxError as exc:
            errors.append(f"{path}: python code block #{i} does not parse (line {exc.lineno}: {exc.msg})")
    for link in MD_LINK.findall(text):
        if link.startswith("file://"):
            errors.append(f"{path}: absolute file:// link {link}")
        elif not link.startswith(("http://", "https://", "mailto:")):
            if not os.path.exists(os.path.normpath(os.path.join(ROOT, os.path.dirname(path), link))):
                errors.append(f"{path}: broken relative link {link}")


def main():
    files = [f for f in target_files("--all" in sys.argv) if os.path.isfile(os.path.join(ROOT, f))]
    errors, warnings = [], []

    for f in files:
        try:
            with open(os.path.join(ROOT, f), encoding="utf-8") as fh:
                m = PRIVATE.search(fh.read())
        except (UnicodeDecodeError, OSError):
            m = None  # binary file
        if m:
            errors.append(f"{f}: local path or personal address ({m.group(0)[:40]!r}); replace with <repo>/<scratchpad>/<conda-base>")
        is_result = "/results_pipeline_" in f or f.startswith("results_pipeline_")
        if not is_result and os.path.getsize(os.path.join(ROOT, f)) > MAX_BYTES:
            errors.append(f"{f}: larger than 5 MB")
        if SOURCE_MD.match(f):
            check_markdown(f, errors)
        elif f.endswith(".py"):
            try:
                ast.parse(read_text(f), filename=f)
            except SyntaxError as exc:
                errors.append(f"{f}: does not parse (line {exc.lineno}: {exc.msg})")

    if "--all" not in sys.argv:
        results = [f for f in files if "/results_pipeline_" in f or f.startswith("results_pipeline_")]
        if results and os.environ.get("ALLOW_RESULTS_CHANGE") != "1":
            errors.append("pipeline results staged (" + ", ".join(results[:3]) + ("…" if len(results) > 3 else "") +
                          "); regenerate in the scRNAseq_anal env and commit with ALLOW_RESULTS_CHANGE=1")
        src_tracks = {f.split("/docs/")[0] if f.startswith("tracks/") else "." for f in files
                      if SOURCE_MD.match(f) and "/docs/" in "/" + f}
        gen_tracks = {os.path.dirname(f) or "." for f in files if GENERATED.search(f)}
        for t in sorted(src_tracks - gen_tracks):
            warnings.append(f"{t}: source chapters changed but the Master Book/HTML/DOCX were not rebuilt (see requirements-build.txt)")
        # Builder changes legitimately regenerate outputs: a shared module covers every track,
        # a track builder covers its own track.
        if any(f in SHARED_BUILD for f in files):
            gen_tracks = set()
        gen_tracks -= {"." if f == "scripts/build_master_documents.py" else f.split("/scripts/")[0]
                       for f in files if TRACK_BUILDER.match(f)}
        for t in sorted(gen_tracks - src_tracks):
            warnings.append(f"{t}: generated Master files changed without source chapters; edit docs/*.md and rebuild instead of hand-editing")

    import importlib.util
    _s = importlib.util.spec_from_file_location("apply_references", os.path.join(ROOT, "scripts", "apply_references.py"))
    _ar = importlib.util.module_from_spec(_s); _s.loader.exec_module(_ar)
    for f in files:
        if SOURCE_MD.match(f) and "/docs/" in "/" + f:
            try:
                st, msg = _ar.check_chapter(f)
            except Exception as exc:
                st, msg = "mismatch", f"{type(exc).__name__}: {exc}"
            if st in ("mismatch", "orphan"):
                errors.append(f"{f}: {msg}")
            elif st == "missing" and "--all" not in sys.argv:
                warnings.append(f"{f}: {msg}")

    for f in files:
        if re.search(r"(^|/)slides/[^/]+\.html$", f) and "SLIDE-LAYOUT:BEGIN" not in read_text(f):
            warnings.append(f"{f}: no SLIDE-LAYOUT marker; run scripts/inject_slide_layout.py after (re)generating decks")

    for w in warnings:
        print(f"[pre-commit] WARN  {w}")
    for e in errors:
        print(f"[pre-commit] ERROR {e}")
    if errors:
        print(f"[pre-commit] {len(errors)} error(s); commit blocked (bypass only with explicit intent: git commit --no-verify)")
        return 1
    print(f"[pre-commit] OK ({len(files)} file(s) checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

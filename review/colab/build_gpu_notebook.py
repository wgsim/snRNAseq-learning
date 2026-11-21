#!/usr/bin/env python3
"""Builds review/colab/gpu_blocks_test.ipynb (nbformat 4).

Every chapter block is extracted VERBATIM from the markdown source by (file, first-line) and embedded as a
string; synthetic setups / post-calls are harness code written here and labelled in the notebook.
Run:  python3 review/colab/build_gpu_notebook.py
"""
import ast, glob, json, os, re, sys
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(os.path.dirname(__file__), "gpu_blocks_test.ipynb")


def block(relglob, start_line):
    """Return (rel_path, start, end, text) for the fenced block whose opening fence is at 1-based `start_line`."""
    path = glob.glob(os.path.join(ROOT, relglob))
    assert len(path) == 1, relglob
    L = open(path[0], encoding="utf-8").read().split("\n")
    assert L[start_line - 1].strip().startswith("```"), (relglob, start_line, L[start_line - 1])
    j = start_line
    while not L[j].strip().startswith("```"):
        j += 1
    text = "\n".join(L[start_line:j]) + "\n"
    ast.parse(text)  # chapter block must be valid Python
    return os.path.relpath(path[0], ROOT), start_line, j + 1, text


BLOCKS = {
    "B1": block("tracks/05_*/docs/00_FOUNDATION_MODELS_ROADMAP.md", 182),
    "B2": block("tracks/05_*/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md", 139),
    "B3": block("tracks/02_*/docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md", 222),
    "B4": block("tracks/08_*/docs/07_GPU_ACCELERATION_RAPIDS_SINGLECELL.md", 184),
    "B5": block("tracks/08_*/docs/08_GPU_GRAPH_COMMUNITY_DETECTION.md", 171),
    "B7": block("tracks/08_*/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md", 205),
    "B6": block("tracks/08_*/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md", 105),
}
for k, (_, _, _, t) in BLOCKS.items():
    assert "'''" not in t, k  # safe to embed in r'''...'''

cells = []
md = lambda s: cells.append(new_markdown_cell(s.strip("\n")))
code = lambda s: cells.append(new_code_cell(s.strip("\n")))

md("""
# GPU code-block test (Colab)

**What this does:** runs the 7 Python code blocks of the course that need an NVIDIA GPU
(cell2location, rapids-singlecell, cuGraph/CAGRA, RMM/CuPy/cuML, PyTorch AMP) and records PASS / FAIL, the exception
and the elapsed time for each block in `gpu_results.json`.

**How to run**
1. Menu **Runtime -> Change runtime type -> Hardware accelerator: T4 GPU** (an L4/A100 also works), then Save.
2. **Runtime -> Run all.** Expected duration: about **15-30 min** (installs 8-15 min, RAPIDS being the slowest; the blocks themselves 5-10 min).
3. If Colab says *"Restart session"* after the install cells, click **Restart**, then **Run all** again (installs are idempotent).
4. At the end, download **`gpu_results.json`** (the last cell tries to start the download; otherwise use the Files pane on the left) and send it back.

**Notes**
- Chapter blocks are executed verbatim from the course sources. Everything marked *HARNESS* (synthetic data, calls into functions that
  the block only defines, shrunk epochs / pool sizes for a 16 GB T4) is test scaffolding, not part of the chapters.
- Each block is isolated (own namespace) and wrapped in try/except; one failure does not stop the others.
- Source of truth for the inventory: `review/fix_log_eng.md`, section "GPU blocks (Colab) (2026-10-21)".
""")

md("## 0. GPU check")
code('''
import subprocess, sys
print(subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout or "nvidia-smi: not found")
try:
    import torch
    HAS_GPU = torch.cuda.is_available()
    print("torch", torch.__version__, "| torch CUDA", torch.version.cuda, "| cuda available:", HAS_GPU)
except Exception as e:  # torch is preinstalled on Colab
    HAS_GPU = False
    print("torch import failed:", e)
if not HAS_GPU:
    raise RuntimeError("NO GPU DETECTED. Use Runtime -> Change runtime type -> T4 GPU, then Run all again. Stopping here.")
GPU_NAME = torch.cuda.get_device_name(0)
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1024**3
print(f"GPU: {GPU_NAME} ({GPU_MEM_GB:.1f} GB)")
''')

md("""
## 1. Installs
Grouped by stack. Package names / extras were checked against PyPI metadata on 2026-10-21:
`rapids-singlecell` (extras `rapids-cu12` / `rapids-cu13`; Python >= 3.12), `cell2location` 0.1.5 (requires `scvi-tools>=1.3.0`),
`scvi-tools` 1.5.1. **Run 2 (2026-10-27): `pyro-ppl<1.9.2` is pinned** because run 1 failed at `import cell2location` (pyro-ppl 1.9.2 no longer defines `deep_getattr`). Other versions are left unpinned on purpose so the blocks are tested against what a learner gets today;
the exact versions are recorded in `gpu_results.json`.
""")
md("### 1a. torch / scvi-tools stack (cell2location, scanpy, transformers)")
code('''
!pip install -q scanpy scvi-tools cell2location "pyro-ppl<1.9.2" transformers psutil 2>&1 | tail -3
''')
md("""
### 1b. RAPIDS stack (cuDF, cuML, cuGraph, cuVS, RMM, CuPy) via rapids-singlecell
Official route: the `rapids-singlecell` extra for your CUDA major version (pulls `cudf`, `cuml`, `cugraph`, `cuvs`, `cupy`, `rmm` from PyPI).
**Expected time: 5-12 minutes** (large wheels). Colab's CUDA major is detected automatically.
""")
code('''
import re, subprocess
smi = subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout
m = re.search(r"CUDA Version:\\s*(\\d+)\\.", smi)
cuda_major = int(m.group(1)) if m else 12
extra = "rapids-cu13" if cuda_major >= 13 else "rapids-cu12"
print("driver CUDA major:", cuda_major, "-> extra:", extra)
!pip install -q "rapids-singlecell[{extra}]" --extra-index-url=https://pypi.nvidia.com 2>&1 | tail -5
''')
md("### 1c. Import check (stack versions are written into the results file)")
code('''
import importlib
VERSIONS = {}
for name in ["torch", "scvi", "cell2location", "scanpy", "anndata", "transformers", "cupy", "cudf", "cuml", "cugraph", "cuvs", "rmm", "rapids_singlecell"]:
    try:
        mod = importlib.import_module(name)
        VERSIONS[name] = getattr(mod, "__version__", "unknown")
    except Exception as e:
        VERSIONS[name] = f"IMPORT FAILED: {type(e).__name__}: {str(e)[:150]}"
for k, v in VERSIONS.items():
    print(f"{k:18s} {v}")
''')

md("## 2. Harness (runs each block in isolation and records the result)")
code('''
import time, traceback, json, warnings, platform
RESULTS = {"meta": {}, "blocks": {}}

def run_block(bid, location, src, setup=None, post=None, replace=(), notes=""):
    """Run chapter block `src` verbatim (optionally with documented textual replacements) in a fresh namespace.
    setup(): HARNESS synthetic data;  post(ns): HARNESS calls into functions the block only defines."""
    rec = {"id": bid, "location": location, "status": None, "phase": None, "error": None,
           "traceback": None, "elapsed_s": None, "notes": notes, "replacements": [list(r) for r in replace],
           "warnings": []}
    ns = {"__name__": "__main__"}
    phase = "setup"
    t0 = time.time()
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        try:
            if setup:
                setup()
            phase = "block"
            s = src
            for a, b in replace:
                assert a in s, f"replacement target not found: {a!r}"
                s = s.replace(a, b)
            exec(compile(s, f"<{bid} {location}>", "exec"), ns)
            if post:
                phase = "post"
                post(ns)
            rec["status"] = "PASS"
        except (Exception, SystemExit) as e:
            rec["status"] = "SETUP_FAIL" if phase == "setup" else "FAIL"
            rec["error"] = f"{type(e).__name__}: {str(e)[:600]}"
            rec["traceback"] = traceback.format_exc()[-2500:]
        finally:
            rec["phase"] = phase
            rec["elapsed_s"] = round(time.time() - t0, 2)
            seen = set()
            for w in wlist:
                msg = f"{w.category.__name__}: {str(w.message)[:200]}"
                if msg not in seen and len(seen) < 12:
                    seen.add(msg)
            rec["warnings"] = sorted(seen)
    RESULTS["blocks"][bid] = rec
    print(f"[{rec['status']}] {bid} ({location}) phase={rec['phase']} {rec['elapsed_s']}s" + (f"\\n    {rec['error']}" if rec["error"] else ""))
    return rec
''')

def block_cell(bid, title, note, setup_src, post_src, replace=None, cls=""):
    rel, a, b, text = BLOCKS[bid]
    loc = f"{rel}:{a}"
    md(f"### {bid}  `{loc}` ({title})\n\n{note}")
    parts = [f"{bid}_SRC = r'''{text}'''"]
    if setup_src:
        parts.append(setup_src.strip("\n"))
    if post_src:
        parts.append(post_src.strip("\n"))
    args = [f'"{bid}"', f'"{loc}"', f"{bid}_SRC"]
    if setup_src:
        args.append(f"setup=_setup_{bid}")
    if post_src:
        args.append(f"post=_post_{bid}")
    if replace:
        args.append(f"replace={replace!r}")
    parts.append(f"run_block({', '.join(args)})")
    code("\n\n".join(parts))


md("## 3. Blocks")

block_cell("B1", "Track 05 ch00 — environment diagnostic; runnable alone",
           "Runnable alone (needs torch, scvi, scanpy, transformers). The chapter's own `[PASS]/[FAIL]` printout is shown, and the cell passes if the script runs to the end.",
           None, None)

block_cell("B2", "Track 05 ch09 — AMP training loop + simplified scIB; runnable alone, plus GPU smoke test",
           "The chapter block runs alone (CPU metrics in its `__main__`). The block only *defines* `train_single_cell_transformer_amp`, so the HARNESS "
           "additionally calls it on a toy model on the GPU (this exercises `torch.cuda.amp.GradScaler/autocast`, which newer PyTorch deprecates).",
           None,
           '''
def _post_B2(ns):
    import torch, torch.nn as nn
    class Toy(nn.Module):  # HARNESS toy model: forward(tokens, values) -> (embedding, logits)
        def __init__(self):
            super().__init__(); self.e = nn.Embedding(50, 16); self.l = nn.Linear(16, 3)
        def forward(self, t, v):
            h = (self.e(t) * v.unsqueeze(-1)).mean(1); return h, self.l(h)
    data = [(torch.randint(0, 50, (8, 10)), torch.rand(8, 10), torch.randint(0, 3, (8,))) for _ in range(8)]
    m = Toy(); opt = torch.optim.Adam(m.parameters(), 1e-3)
    out = ns["train_single_cell_transformer_amp"](m, data, opt, num_epochs=2, gradient_accumulation_steps=4, device="cuda")
    assert next(out.parameters()).is_cuda
''')

block_cell("B3", "Track 02 ch05 — cell2location (use_gpu=True); needs synthetic setup",
           "The block reads `data/mcrpc_reference_snrnaseq.h5ad` and `data/visium_prostate_cancer.h5ad`; the HARNESS writes small synthetic versions "
           "(600 reference cells, 150 Visium spots, 200 genes, fake H&E image). **Documented replacement:** `max_epochs` 250 -> 20 and 30000 -> 300 so the cell finishes in minutes "
           "(set `FULL_EPOCHS = True` to run the chapter values; that can take 15+ min). Run 1 showed `use_gpu=True` is rejected by the installed scvi-tools Trainer; the chapter now uses `accelerator=\"gpu\"` (also inside `sample_kwargs`).",
           '''
FULL_EPOCHS = False
def _setup_B3():
    import numpy as np, pandas as pd, anndata as ad, os, scipy.sparse as sp
    rng = np.random.default_rng(0); os.makedirs("data", exist_ok=True)
    G = 200; genes = [f"GENE{i}" for i in range(G)]
    types = ["Luminal_AR_plus", "t_NEPC", "myCAF", "CD8_T_cell", "Basal"]
    sig = rng.gamma(1.0, 1.0, (len(types), G)) + 0.2
    ct = rng.integers(0, len(types), 600)
    X = rng.poisson(sig[ct] * 3).astype("float32")
    ref = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({
        "patient_id": rng.choice(["P1", "P2", "P3"], 600),
        "lineage_4states": pd.Categorical(np.array(types)[ct], categories=types),
        "technology": rng.choice(["A", "B"], 600)}, index=[f"r{i}" for i in range(600)]),
        var=pd.DataFrame(index=genes))
    ref.write_h5ad("data/mcrpc_reference_snrnaseq.h5ad")
    S = 150; mix = rng.dirichlet(np.ones(len(types)), S)
    V = rng.poisson(mix @ sig * 15).astype("float32")
    vis = ad.AnnData(sp.csr_matrix(V), obs=pd.DataFrame({"sample_slide_id": rng.choice(["S1", "S2"], S)}, index=[f"s{i}" for i in range(S)]),
                     var=pd.DataFrame(index=genes))
    vis.obsm["spatial"] = rng.uniform(10, 90, (S, 2))
    vis.uns["spatial"] = {"lib": {"images": {"hires": rng.uniform(0.6, 1, (100, 100, 3)).astype("float32")},
                                  "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 4.0}}}
    vis.write_h5ad("data/visium_prostate_cancer.h5ad")
''',
           None,
           replace=None)
# B3 needs the FULL_EPOCHS switch: patch the generated call after creation
cells[-1].source = cells[-1].source.replace(
    'run_block("B3", "' + BLOCKS["B3"][0] + ':' + str(BLOCKS["B3"][1]) + '", B3_SRC, setup=_setup_B3)',
    'run_block("B3", "' + BLOCKS["B3"][0] + ':' + str(BLOCKS["B3"][1]) + '", B3_SRC, setup=_setup_B3,\n'
    '          replace=[] if FULL_EPOCHS else [("max_epochs=250", "max_epochs=20"), ("max_epochs=30000", "max_epochs=300")])')
assert "FULL_EPOCHS else" in cells[-1].source

block_cell("B4", "Track 08 ch07 — rapids-singlecell pipeline (functions only); needs synthetic setup + call",
           "The block only defines `init_gpu_memory_pool`, `load_1m_atlas_to_gpu`, `run_gpu_accelerated_pipeline`. The HARNESS builds a small synthetic atlas "
           "(6000 cells x 4000 genes, with `MT-` genes), writes it to h5ad, and calls the three functions in order. **Documented replacement of arguments:** RMM pool "
           "`60GB/76GB` -> `1GB/4GB` (the chapter defaults target an 80 GB GPU; a T4 has 16 GB). If this block fails with the *default* pool sizes only, that is expected on a small GPU.",
           None,
           '''
def _post_B4(ns):
    import numpy as np, scipy.sparse as sp, anndata as ad, pandas as pd
    rng = np.random.default_rng(1); n, g = 6000, 4000
    X = sp.csr_matrix(rng.poisson(rng.gamma(0.5, 1.0, g) * 1.5, (n, g)).astype("float32"))
    names = [f"MT-{i}" if i < 13 else f"G{i}" for i in range(g)]
    a = ad.AnnData(X, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=names))
    a.write_h5ad("synthetic_atlas.h5ad")
    ns["init_gpu_memory_pool"](initial_pool_size="1GB", max_pool_size="4GB")   # HARNESS: smaller pool
    a = ns["load_1m_atlas_to_gpu"]("synthetic_atlas.h5ad")
    a = ns["run_gpu_accelerated_pipeline"](a)
    assert "X_umap" in a.obsm
''')

block_cell("B5", "Track 08 ch08 — CAGRA kNN + cuGraph Leiden (functions only); needs synthetic setup + call",
           "The block itself says the CAGRA part was *not run on a GPU* (`GPU 미실행`), so this is the first real check of `cuvs.neighbors.cagra` argument names and return types. "
           "The HARNESS calls the three defined functions on 20,000 random 50-d points (float32).",
           None,
           '''
def _post_B5(ns):
    import cupy as cp
    cp.random.seed(0)
    X = cp.random.standard_normal((20000, 50), dtype=cp.float32)   # HARNESS synthetic PCs
    idx, dist = ns["build_cagra_knn_graph_1m"](X, k=15)
    assert idx.shape == (20000, 15) and dist.shape == (20000, 15)
    G = ns["convert_knn_to_cugraph"](idx, dist)
    parts = ns["run_cugraph_leiden_clustering"](G, resolution=0.8, max_iter=50)
    assert len(parts) == 20000
''')

block_cell("B7", "Track 08 ch09 — memory guard + cuML PCA (functions only); needs synthetic setup + call",
           "The block defines a decorator and `execute_safe_pca`; the HARNESS calls `execute_safe_pca` on a random 5000 x 100 CuPy matrix.",
           None,
           '''
def _post_B7(ns):
    import cupy as cp
    X = cp.random.standard_normal((5000, 100), dtype=cp.float32)   # HARNESS synthetic matrix
    out = ns["execute_safe_pca"](X, n_components=10)
    assert out.shape == (5000, 10)
''')

block_cell("B6", "Track 08 ch09 — RMM managed memory; runnable after one documented replacement",
           "Run last because it re-initialises the RMM allocator for the whole process. **Documented replacement:** `initial_pool_size=\"64GB\"` -> `\"2GB\"` "
           "(chapter value targets an 80 GB GPU). Verbatim run is not attempted on a 16 GB GPU.",
           None, None,
           replace=[('initial_pool_size="64GB"', 'initial_pool_size="2GB"')])

md("## 4. Write results")
code('''
import json, platform, subprocess
RESULTS["meta"] = {
    "date": time.strftime("%Y-%m-%d %H:%M:%S"), "python": platform.python_version(),
    "gpu": GPU_NAME, "gpu_mem_gb": round(GPU_MEM_GB, 1), "versions": VERSIONS,
}
try:
    RESULTS["meta"]["nvidia_smi_head"] = subprocess.run(["nvidia-smi", "--query-gpu=driver_version,name", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
except Exception as e:
    RESULTS["meta"]["nvidia_smi_head"] = f"unavailable: {e}"
with open("gpu_results.json", "w") as f:
    json.dump(RESULTS, f, indent=2, ensure_ascii=False)
print(f"{'block':5s} {'status':10s} {'sec':>8s}  location")
for bid, r in RESULTS["blocks"].items():
    print(f"{bid:5s} {r['status']:10s} {r['elapsed_s']:>8}  {r['location']}" + (f"  -> {r['error'][:90]}" if r["error"] else ""))
n_pass = sum(r["status"] == "PASS" for r in RESULTS["blocks"].values())
print(f"\\n{n_pass}/{len(RESULTS['blocks'])} blocks PASS. Results saved to gpu_results.json")
print(">>> Download gpu_results.json (Files pane on the left, or the download starting now) and send it back. <<<")
try:
    from google.colab import files
    files.download("gpu_results.json")
except Exception as e:
    print("(automatic download unavailable:", type(e).__name__, ") -> use the Files pane.")
''')

nb = new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
    "accelerator": "GPU",
    "colab": {"provenance": [], "gpuType": "T4"},
})
nbformat.validate(nb)
nbformat.write(nb, OUT)
# syntax check: every code cell (magics lines removed) and every embedded block string
bad = 0
for i, c in enumerate(nb.cells):
    if c.cell_type != "code":
        continue
    src = "\n".join(l for l in c.source.split("\n") if not l.lstrip().startswith("!"))
    try:
        ast.parse(src)
    except SyntaxError as e:
        bad += 1
        print("SYNTAX ERROR in cell", i, e)
print(f"wrote {OUT}: {len(nb.cells)} cells, syntax errors: {bad}")
sys.exit(1 if bad else 0)

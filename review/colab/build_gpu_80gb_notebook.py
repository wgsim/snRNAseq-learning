#!/usr/bin/env python3
"""Builds review/colab/gpu_blocks_80gb.ipynb (nbformat 4): the 80 GB-GPU profile of the GPU code-block test.

Differences from gpu_blocks_test.ipynb (T4 profile, see build_gpu_notebook.py):
  * the RMM pool sizes are the chapter defaults (60GB/76GB in T08 ch07, 64GB in T08 ch09) - no downsizing;
    the notebook refuses to run on a GPU with < 77e9 bytes;
  * data is synthetic but sized 25k and 100k cells x 35,000 genes (~5% non-zero), i.e. a ramp that stops at 100k;
  * every block runs in its own subprocess (RMM re-initialisation is process-wide) and records peak GPU memory
    (nvidia-smi polling) and peak host RSS;
  * only the blocks whose defaults depend on GPU memory (B4, B5, B7, B6) are run; B1-B3 are in the T4 profile.
Chapter blocks are embedded verbatim, extracted by (file, opening-fence line) like build_gpu_notebook.py.
Run:  uv run --with nbformat python review/colab/build_gpu_80gb_notebook.py
Options: --only B5_25k,B5_100k --min-gpu-bytes 0 --tag b5_a100 --out <path>  builds a subset for another GPU (results in gpu_results_<tag>.json);
         --replace "OLD=>NEW" (repeatable) applies a documented textual replacement to every block that contains OLD (recorded in the results).
"""
import argparse, ast, glob, os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--only", default="", help="comma-separated run ids (B4_25k,B4_100k,B5_25k,B5_100k,B7_100k,B6); default all")
ap.add_argument("--min-gpu-bytes", type=int, default=77 * 10**9, help="refuse to run below this GPU memory (0 = no gate)")
ap.add_argument("--tag", default="80gb", help="results file is gpu_results_<tag>.json")
ap.add_argument("--replace", action="append", default=[], help="OLD=>NEW textual replacement applied to chapter blocks that contain OLD")
ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "gpu_blocks_80gb.ipynb"))
ARGS = ap.parse_args()
OUT = ARGS.out
REPLACE = [tuple(r.split("=>", 1)) for r in ARGS.replace]
ONLY = {x for x in ARGS.only.split(",") if x}
want = lambda rid: not ONLY or rid in ONLY


def block(relglob, start_line):
    path = glob.glob(os.path.join(ROOT, relglob))
    assert len(path) == 1, relglob
    L = open(path[0], encoding="utf-8").read().split("\n")
    assert L[start_line - 1].strip().startswith("```"), (relglob, start_line, L[start_line - 1])
    j = start_line
    while not L[j].strip().startswith("```"):
        j += 1
    text = "\n".join(L[start_line:j]) + "\n"
    ast.parse(text)
    return os.path.relpath(path[0], ROOT), start_line, j + 1, text


BLOCKS = {
    "B4": block("tracks/08_*/docs/07_GPU_ACCELERATION_RAPIDS_SINGLECELL.md", 184),
    "B5": block("tracks/08_*/docs/08_GPU_GRAPH_COMMUNITY_DETECTION.md", 171),
    "B7": block("tracks/08_*/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md", 205),
    "B6": block("tracks/08_*/docs/09_TROUBLESHOOTING_AND_BEST_PRACTICES.md", 105),
}
for k, (_, _, _, t) in BLOCKS.items():
    assert "'''" not in t, k

cells = []
md = lambda s: cells.append(new_markdown_cell(s.strip("\n")))
code = lambda s: cells.append(new_code_cell(s.strip("\n")))

md("""
# GPU code-block test, 80 GB profile (Colab A100/H100 80 GB)

Runs the 4 Track 08 blocks whose defaults assume an 80 GB GPU **with the chapter defaults** (RMM pool 60GB/76GB and 64GB,
memory-guard thresholds 70/65 GB) on synthetic data of 25k and 100k cells x 35,000 genes. Each block runs in its own process and
the peak GPU memory (nvidia-smi polling, includes ~0.4 GB of the parent CUDA context) and the peak host RSS are recorded in `gpu_results_80gb.json`.
The ramp stops at 100k cells on purpose (host-RAM risk of a 1M x 35k synthetic atlas was not taken).
""")
md(f"## 0. GPU check (refuses to run below {ARGS.min_gpu_bytes:,} bytes)")
code('''
import subprocess, os, torch
print(subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout)
assert torch.cuda.is_available(), "NO GPU DETECTED"
props = torch.cuda.get_device_properties(0)
GPU_NAME, GPU_MEM_BYTES = props.name, props.total_memory
print(f"GPU: {GPU_NAME} | total_memory = {GPU_MEM_BYTES} bytes ({GPU_MEM_BYTES/1024**3:.1f} GiB)")
if GPU_MEM_BYTES < __MIN_BYTES__:
    raise RuntimeError(f"GPU has only {GPU_MEM_BYTES/1e9:.1f} GB; the chapter default RMM max pool is 76GB (decimal). Stopping before the installs.")
'''.replace("__MIN_BYTES__", str(ARGS.min_gpu_bytes)))
md("## 1. Install (RAPIDS only; scvi-tools / cell2location are not needed for these blocks)")
code('''
import re, subprocess
smi = subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout
m = re.search(r"CUDA Version:\\s*(\\d+)\\.", smi)
cuda_major = int(m.group(1)) if m else 12
extra = "rapids-cu13" if cuda_major >= 13 else "rapids-cu12"
print("driver CUDA major:", cuda_major, "-> extra:", extra)
!pip install -q "rapids-singlecell[{extra}]" psutil --extra-index-url=https://pypi.nvidia.com 2>&1 | tail -5
''')
code('''
import importlib, os, psutil
VERSIONS = {}
for name in ["torch", "scanpy", "anndata", "cupy", "cudf", "cuml", "cugraph", "cuvs", "rmm", "rapids_singlecell"]:
    try:
        VERSIONS[name] = getattr(importlib.import_module(name), "__version__", "unknown")
    except Exception as e:
        VERSIONS[name] = f"IMPORT FAILED: {type(e).__name__}: {str(e)[:150]}"
HOST = {"cpus": os.cpu_count(), "ram_gb": round(psutil.virtual_memory().total / 1e9, 1)}
print(VERSIONS); print(HOST)
''')
md("## 2. Harness: every block in its own subprocess, with peak GPU / host memory")
code("""
import json, re, subprocess, sys, time
RESULTS = {"meta": {}, "blocks": {}}
CHILD = r'''
import json, sys, time, threading, subprocess, traceback, resource, warnings
payload = json.load(open(sys.argv[1]))
bid, loc = payload["bid"], payload["loc"]
rec = {"id": bid, "location": loc, "status": None, "error": None, "traceback": None}
peak = {"gpu_mib": 0}
stop = threading.Event()
def poll():
    while not stop.is_set():
        try:
            out = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                                 capture_output=True, text=True, timeout=5).stdout.split()[0]
            peak["gpu_mib"] = max(peak["gpu_mib"], int(out))
        except Exception:
            pass
        stop.wait(0.5)
threading.Thread(target=poll, daemon=True).start()
ns = {"__name__": "__main__"}
t0 = time.time()
with warnings.catch_warnings(record=True) as wl:
    warnings.simplefilter("always")
    try:
        exec(compile(payload["block"], "<%s %s>" % (bid, loc), "exec"), ns)
        exec(compile(payload["post"], "<%s HARNESS>" % bid, "exec"), {"ns": ns, "__name__": "__main__"})
        rec["status"] = "PASS"
    except BaseException as e:
        rec["status"] = "FAIL"
        rec["error"] = "%s: %s" % (type(e).__name__, str(e)[:600])
        rec["traceback"] = traceback.format_exc()[-2500:]
    rec["warnings"] = sorted({"%s: %s" % (w.category.__name__, str(w.message)[:160]) for w in wl})[:8]
rec["elapsed_s"] = round(time.time() - t0, 2)
rec["peak_gpu_used_mib"] = peak["gpu_mib"]
rec["peak_host_rss_mib"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
stop.set()
print("@@RESULT@@" + json.dumps(rec))
'''

def run_block(bid, location, src, post, notes="", timeout=1800, replace=()):
    applied = []
    for a, b in replace:
        if a in src:
            src = src.replace(a, b)
            applied.append([a, b])
    pj = f"/content/_blk_{bid}.json"
    json.dump({"bid": bid, "loc": location, "block": src, "post": post}, open(pj, "w"))
    open("/content/_child.py", "w").write(CHILD)
    try:
        p = subprocess.run([sys.executable, "/content/_child.py", pj], capture_output=True, text=True, timeout=timeout)
        out, err, rc = p.stdout, p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        tx = lambda b: b.decode() if isinstance(b, bytes) else (b or "")
        out, err, rc = tx(e.stdout), tx(e.stderr), "TIMEOUT"
    m = re.search(r"@@RESULT@@(.*)", out)
    rec = json.loads(m.group(1)) if m else {"id": bid, "location": location, "status": "CRASH" if rc != "TIMEOUT" else "TIMEOUT",
                                             "error": f"no result line; returncode={rc}", "elapsed_s": None}
    rec["stdout_tail"] = [l for l in out.splitlines() if "@@RESULT@@" not in l][-30:]
    rec["stderr_tail"] = err.splitlines()[-8:]
    rec["notes"] = notes
    rec["replacements"] = applied
    RESULTS["blocks"][bid] = rec
    print(f"[{rec['status']}] {bid} {rec.get('elapsed_s')}s  peak GPU {rec.get('peak_gpu_used_mib')} MiB  peak host RSS {rec.get('peak_host_rss_mib')} MiB"
          + (f"\\n    {rec['error']}" if rec.get("error") else ""))
    print("    " + "\\n    ".join(rec["stdout_tail"][-14:]))
    return rec
""")
NEED_ATLAS = [n for n in (25_000, 100_000) if want(f"B4_{n//1000}k")]
if not NEED_ATLAS:
    code("NNZ = {}  # no B4 run requested -> no synthetic atlas")
else:
  md("## 3. Synthetic atlases (35,000 genes, ~5% non-zero, first 13 genes named `MT-*`)")
  code('''
import numpy as np, scipy.sparse as sp, anndata as ad, pandas as pd
G = 35000
def make_atlas(n, path, chunk=2500, seed=1):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    lam = rng.gamma(0.3, 0.2, G)                       # per-gene Poisson rate -> ~5% non-zero entries
    parts = [sp.csr_matrix(rng.poisson(lam, (min(chunk, n - i), G)).astype("float32")) for i in range(0, n, chunk)]
    X = sp.vstack(parts, format="csr")
    names = [f"MT-{j}" if j < 13 else f"G{j}" for j in range(G)]
    ad.AnnData(X, obs=pd.DataFrame(index=[f"c{j}" for j in range(n)]), var=pd.DataFrame(index=names)).write_h5ad(path)
    print(f"{path}: {n} x {G}, nnz={X.nnz:,} ({X.nnz/(n*G):.1%}), {time.time()-t0:.0f}s, host RAM used {psutil.virtual_memory().used/1e9:.1f} GB")
    return X.nnz
SIZES = __SIZES__
NNZ = {n: make_atlas(n, f"/content/atlas_{n//1000}k.h5ad") for n in SIZES}
'''.replace("__SIZES__", repr(NEED_ATLAS)))
md("## 4. Blocks (chapter defaults; the HARNESS post code only builds inputs and calls the functions the block defines)")

POST = {
    "B4": '''
ns["init_gpu_memory_pool"]()          # chapter defaults: initial 60GB, max 76GB (no override)
a = ns["load_1m_atlas_to_gpu"]("__PATH__")
a = ns["run_gpu_accelerated_pipeline"](a)
assert "X_umap" in a.obsm
print("HARNESS: cells after QC =", a.n_obs, "| genes after filter =", a.n_vars)
''',
    "B5": '''
import cupy as cp
cp.random.seed(0)
n = __N__
X = cp.random.standard_normal((n, 50), dtype=cp.float32)      # HARNESS synthetic PCs
idx, dist = ns["build_cagra_knn_graph_1m"](X, k=15)
assert idx.shape == (n, 15) and dist.shape == (n, 15)
G = ns["convert_knn_to_cugraph"](idx, dist)
parts = ns["run_cugraph_leiden_clustering"](G, resolution=0.8, max_iter=50)
assert len(parts) == n
print("HARNESS: leiden returned", len(parts), "rows")
''',
    "B7": '''
import cupy as cp
n = __N__
X = cp.random.standard_normal((n, 3000), dtype=cp.float32)    # HARNESS synthetic matrix (3000 HVGs)
out = ns["execute_safe_pca"](X, n_components=50)
assert out.shape == (n, 50)
print("HARNESS: PCA output", out.shape)
''',
    "B6": '''
import cupy as cp
from rmm.allocators.cupy import rmm_cupy_allocator
cp.cuda.set_allocator(rmm_cupy_allocator)
x = cp.ones(4 * 1024**3 // 4, dtype=cp.float32)               # HARNESS: 4 GiB inside the 64GB managed pool
cp.cuda.runtime.deviceSynchronize()
print("HARNESS: 4 GiB allocation inside the managed pool OK")
del x
try:                                                          # HARNESS probe of the chapter claim "exceeding physical 80GB VRAM"
    y = cp.empty(90 * 10**9 // 4, dtype=cp.float32)
    print("HARNESS PROBE: 90 GB managed allocation SUCCEEDED")
    del y
except Exception as e:
    print("HARNESS PROBE: 90 GB managed allocation FAILED:", type(e).__name__, str(e)[:200])
''',
}


def block_cell(bid, runid, title, note, post, notes=""):
    rel, start, end, text = BLOCKS[bid]
    loc = f"{rel}:{start}"
    md(f"### {runid}  `{loc}` ({title})\n{note}")
    code(f"{bid}_SRC = r'''{text}'''\n"
         f"POST = {post!r}\n"
         f"run_block({runid!r}, {loc!r}, {bid}_SRC, POST, notes={notes!r}, replace={REPLACE!r})")


for n in NEED_ATLAS:
    block_cell("B4", f"B4_{n//1000}k", "T08 ch07 rapids-singlecell pipeline",
               f"Chapter-default RMM pool 60GB/76GB, synthetic atlas {n:,} x 35,000.",
               POST["B4"].replace("__PATH__", f"/content/atlas_{n//1000}k.h5ad"))
for n in (25_000, 100_000):
    if want(f"B5_{n//1000}k"):
      block_cell("B5", f"B5_{n//1000}k", "T08 ch08 CAGRA + cuGraph Leiden",
                 f"{n:,} random 50-d points (the chapter targets 1M).", POST["B5"].replace("__N__", str(n)))
if want("B7_100k"):
    block_cell("B7", "B7_100k", "T08 ch09 memory guard + cuML PCA",
               "100,000 x 3,000 float32; guard thresholds are the chapter's (70/65 GB), so the defragmentation branch is expected NOT to trigger.",
               POST["B7"].replace("__N__", "100000"))
if want("B6"):
    block_cell("B6", "B6", "T08 ch09 RMM managed memory",
               "Verbatim `initial_pool_size=\"64GB\"` (run last: RMM re-initialisation). Plus a 4 GiB allocation and an informational 90 GB oversubscription probe.",
               POST["B6"])

md("## 5. Write results")
code('''
import json
RESULTS["meta"] = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "gpu": GPU_NAME, "gpu_mem_bytes": GPU_MEM_BYTES,
                   "host": HOST, "versions": VERSIONS, "atlas_nnz": {str(k): v for k, v in NNZ.items()}}
json.dump(RESULTS, open("gpu_results___TAG__.json", "w"), indent=2, ensure_ascii=False)
print(f"{'block':9s} {'status':8s} {'sec':>8s} {'peakGPU MiB':>12s} {'peakRSS MiB':>12s}")
for k, r in RESULTS["blocks"].items():
    print(f"{k:9s} {r['status']:8s} {str(r.get('elapsed_s')):>8s} {str(r.get('peak_gpu_used_mib')):>12s} {str(r.get('peak_host_rss_mib')):>12s}")
print(sum(r["status"] == "PASS" for r in RESULTS["blocks"].values()), "/", len(RESULTS["blocks"]), "PASS; saved gpu_results___TAG__.json")
'''.replace("__TAG__", ARGS.tag))

nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
bad = 0
for c in nb.cells:
    if c.cell_type == "code":
        src = "\n".join(l for l in c.source.split("\n") if not l.startswith("!"))
        try:
            ast.parse(src)
        except SyntaxError as e:
            bad += 1
            print("SYNTAX ERROR", e, "\n", c.source[:200])
nbformat.write(nb, OUT)
print(f"wrote {OUT}: {len(cells)} cells, syntax errors: {bad}")

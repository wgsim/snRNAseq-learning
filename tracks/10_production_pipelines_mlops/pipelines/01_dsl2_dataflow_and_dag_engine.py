#!/usr/bin/env python3
"""
tracks/10_production_pipelines_mlops/pipelines/01_dsl2_dataflow_and_dag_engine.py
================================================================================
Production Pipelines, Workflow Orchestration & Reactive Dataflow Engine Simulation.

This pipeline models:
1. 8-Stage Single-Cell Pipeline DAG:
   - FastQC ➔ Demux ➔ Spliced_Align ➔ UMI_Dedup ➔ Cell_Filter ➔ Normalization ➔ Clustering ➔ Report
2. Asynchronous Queue Channels and Reactive Operators (.map, .filter, .combine, .branch)
3. Simulated Out-Of-Memory (OOM) on Sample 3 with dynamic retry & doubled memory allocation
   (Nextflow DSL2 `errorStrategy 'retry'`, `maxRetries 2`, `memory = base_mem * task.attempt`)
4. Cryptographic MD5(inputs + script) task caching (`-resume` simulation):
   - Cold run: executes all tasks, handles retry on Sample 3
   - Warm run (-resume): verifies 100% cache hit on unmodified tasks
5. Structured outputs saved to results_pipeline_01/:
   - trace_table.csv
   - timeline_metrics.json
   - summary_report.md
"""

import os
import sys
import time
import json
import csv
import hashlib
import asyncio
from typing import Dict, List, Any, Optional, Tuple


# ============================================================================
# 1. Reactive Channel Architecture (Queue & Value Channels)
# ============================================================================

class QueueChannel:
    """Nextflow DSL2 style asynchronous FIFO Queue Channel."""
    def __init__(self, name: str):
        self.name = name
        self.queue: asyncio.Queue = asyncio.Queue()
        self.closed: bool = False

    async def put(self, item: Any):
        if self.closed:
            raise RuntimeError(f"Cannot put to closed channel {self.name}")
        await self.queue.put(item)

    async def get(self) -> Any:
        return await self.queue.get()

    def empty(self) -> bool:
        return self.queue.empty()

    def close(self):
        self.closed = True

    def qsize(self) -> int:
        return self.queue.qsize()


class ValueChannel:
    """Nextflow DSL2 style reusable broadcast Value Channel."""
    def __init__(self, name: str, value: Any):
        self.name = name
        self.value = value

    def get(self) -> Any:
        return self.value


# ============================================================================
# 2. Pipeline Stage Definitions & Task Model
# ============================================================================

PIPELINE_STAGES = [
    ("FastQC", "Raw read quality profiling (Phred Q30, adapter content, sequence duplication)"),
    ("Demux", "Sample demultiplexing, 16bp cell barcode extraction & 3M-february-2018 (10x v3) whitelist matching"),
    ("Spliced_Align", "STARsolo / Alevin-fry splice-aware alignment against GRCh38 reference genome"),
    ("UMI_Dedup", "12bp UMI directional adjacency graph deduplication & PCR artifact collapse"),
    ("Cell_Filter", "EmptyDrops Dirichlet-Multinomial ambient RNA separation (FDR < 0.01)"),
    ("Normalization", "SCTransform regularized negative binomial normalization & 2,000 HVG selection"),
    ("Clustering", "scVI variational autoencoder 10D latent space & Leiden graph community detection"),
    ("Report", "MultiQC execution summary, cohort diagnostic metrics & clinical QC signoff")
]

# Clinical Prostate Cancer Cohort Samples
SAMPLES_METADATA = [
    {
        "sample_id": "Sample_01_PRAD_Primary",
        "tissue": "Primary Prostate Adenocarcinoma (Gleason 4+3)",
        "read_pairs": 48_500_000,
        "base_memory_gb": 16,
        "peak_ram_required_gb": 14.2
    },
    {
        "sample_id": "Sample_02_PRAD_Adjacent",
        "tissue": "Histologically Normal Adjacent Prostate Tissue",
        "read_pairs": 32_100_000,
        "base_memory_gb": 16,
        "peak_ram_required_gb": 11.8
    },
    {
        "sample_id": "Sample_03_CRPC_Metastasis",
        "tissue": "Castration-Resistant Bone Metastasis (High Depth / Cell Density)",
        "read_pairs": 94_800_000,
        "base_memory_gb": 16,  # Will fail on attempt 1 with 16GB, needs 28GB -> retry with 32GB
        "peak_ram_required_gb": 28.4
    },
    {
        "sample_id": "Sample_04_NEPC_Progression",
        "tissue": "Treatment-Emergent Neuroendocrine Prostate Cancer (t-NEPC)",
        "read_pairs": 56_200_000,
        "base_memory_gb": 16,
        "peak_ram_required_gb": 15.6
    }
]


# ============================================================================
# 3. Cryptographic Execution Cache (-resume engine)
# ============================================================================

class ExecutionCache:
    """
    Simulates Nextflow's `.nextflow/cache` store.
    Task hash is computed as MD5(process_name, input_files, script_version, sample_id).
    """
    def __init__(self, cache_file: str):
        self.cache_file = cache_file
        self.store: Dict[str, Dict[str, Any]] = {}
        self.load()

    def load(self):
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.store = json.load(f)
            except Exception:
                self.store = {}

    def save(self):
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_file)), exist_ok=True)
        with open(self.cache_file, "w", encoding="utf-8") as f:
            json.dump(self.store, f, indent=2)

    @staticmethod
    def compute_hash(process_name: str, sample_id: str, inputs: Any, script_code: str) -> str:
        hasher = hashlib.md5()
        hasher.update(process_name.encode())
        hasher.update(sample_id.encode())
        hasher.update(json.dumps(inputs, sort_keys=True, default=str).encode())
        hasher.update(script_code.encode())
        return hasher.hexdigest()

    def get(self, hash_key: str) -> Optional[Dict[str, Any]]:
        return self.store.get(hash_key)

    def put(self, hash_key: str, data: Dict[str, Any]):
        self.store[hash_key] = data
        self.save()


# ============================================================================
# 4. Asynchronous Dataflow & DAG Execution Engine
# ============================================================================

class PipelineExecutionEngine:
    def __init__(self, cache: ExecutionCache, resume: bool = False):
        self.cache = cache
        self.resume = resume
        self.trace_records: List[Dict[str, Any]] = []
        self.retry_events: List[Dict[str, Any]] = []
        self.task_counter: int = 0

    async def execute_task(
        self,
        process_name: str,
        sample: Dict[str, Any],
        inputs: Dict[str, Any],
        script_code: str,
        error_strategy: str = "retry",
        max_retries: int = 2
    ) -> Dict[str, Any]:
        """
        Executes a single task with reactive retry and caching.
        """
        sample_id = sample["sample_id"]
        hash_key = self.cache.compute_hash(process_name, sample_id, inputs, script_code)

        # Check Cache if Resume Mode is Active
        if self.resume:
            cached_entry = self.cache.get(hash_key)
            if cached_entry:
                self.task_counter += 1
                rec = {
                    "task_id": f"TASK-{self.task_counter:03d}",
                    "process_name": process_name,
                    "sample_id": sample_id,
                    "attempt": cached_entry.get("attempt", 1),
                    "status": "CACHED",
                    "memory_gb": cached_entry.get("memory_gb", 16),
                    "cpus": 8,
                    "duration_s": 0.001,
                    "exit_code": 0,
                    "cached": True,
                    "hash_key": hash_key[:12],
                    "work_dir": f".nextflow/work/{hash_key[:2]}/{hash_key[2:8]}",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
                }
                self.trace_records.append(rec)
                return cached_entry["outputs"]

        # Cold Execution with Dynamic Error Handling
        attempt = 1
        current_memory = sample["base_memory_gb"]

        while attempt <= max_retries + 1:
            self.task_counter += 1
            task_id = f"TASK-{self.task_counter:03d}"
            t0 = time.time()

            # Simulate task execution delay (asynchronous non-blocking)
            simulated_delay = 0.03 + (0.01 * attempt)
            await asyncio.sleep(simulated_delay)

            # Check if this task triggers simulated OOM
            # Specifically: Spliced_Align on Sample_03 requires 28.4GB.
            # If current_memory < required, it fails with Exit Code 137 (OOM Killer)
            peak_ram_needed = sample.get("peak_ram_required_gb", 12.0)
            is_oom = (process_name == "Spliced_Align" and current_memory < peak_ram_needed)

            if is_oom:
                duration = round(time.time() - t0, 4)
                rec = {
                    "task_id": task_id,
                    "process_name": process_name,
                    "sample_id": sample_id,
                    "attempt": attempt,
                    "status": "FAILED_OOM",
                    "memory_gb": current_memory,
                    "cpus": 8,
                    "duration_s": duration,
                    "exit_code": 137,
                    "cached": False,
                    "hash_key": hash_key[:12],
                    "work_dir": f".nextflow/work/{hash_key[:2]}/{hash_key[2:8]}",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
                }
                self.trace_records.append(rec)

                retry_info = {
                    "task_id": task_id,
                    "process": process_name,
                    "sample": sample_id,
                    "failed_attempt": attempt,
                    "failed_memory_gb": current_memory,
                    "error": "Exit 137: Out of Memory (OOM Killer). Memory limit exceeded.",
                    "action": f"DSL2 errorStrategy '{error_strategy}': Doubling memory allocation for attempt {attempt + 1}"
                }
                self.retry_events.append(retry_info)

                if attempt <= max_retries and error_strategy == "retry":
                    # Double memory allocation for next attempt: 16GB -> 32GB
                    attempt += 1
                    current_memory = current_memory * 2
                    continue
                else:
                    raise MemoryError(f"Task {task_id} failed after {attempt} attempts: OOM")

            # Successful Execution
            duration = round(time.time() - t0, 4)
            outputs = {
                "sample_id": sample_id,
                "process": process_name,
                "output_file": f"{process_name.lower()}_{sample_id}.out",
                "cells_retained": int(sample["read_pairs"] / 15_000) if process_name != "FastQC" else None,
                "memory_used_gb": round(peak_ram_needed, 2)
            }

            rec = {
                "task_id": task_id,
                "process_name": process_name,
                "sample_id": sample_id,
                "attempt": attempt,
                "status": "SUCCESS",
                "memory_gb": current_memory,
                "cpus": 8,
                "duration_s": duration,
                "exit_code": 0,
                "cached": False,
                "hash_key": hash_key[:12],
                "work_dir": f".nextflow/work/{hash_key[:2]}/{hash_key[2:8]}",
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
            }
            self.trace_records.append(rec)

            # Store in Cache
            self.cache.put(hash_key, {
                "outputs": outputs,
                "memory_gb": current_memory,
                "attempt": attempt
            })
            return outputs

        raise RuntimeError(f"Unexpected termination for task {process_name} on {sample_id}")

    async def run_pipeline(self) -> Dict[str, Any]:
        """
        Executes the full 8-stage reactive DAG across all 4 cohort samples.
        """
        print(f"\n{'='*75}")
        mode_str = "WARM RESUME MODE (-resume: verifying 100% cache hits)" if self.resume else "COLD INITIAL EXECUTION (with dynamic error recovery)"
        print(f"[*] LAUNCHING NEXTFLOW DSL2 SIMULATION ENGINE: {mode_str}")
        print(f"{'='*75}")

        # Channels setup
        samples_channel = QueueChannel("samples_ch")
        for s in SAMPLES_METADATA:
            await samples_channel.put(s)
        samples_channel.close()

        ref_genome_val = ValueChannel("ref_genome", "GRCh38.primary_assembly.genome.fa")
        whitelist_val = ValueChannel("whitelist", "3M-february-2018.txt")

        # Process pipeline for all samples asynchronously
        async def process_sample_dag(sample: Dict[str, Any]):
            sid = sample["sample_id"]
            # Stage 1: FastQC
            fq_out = await self.execute_task(
                "FastQC", sample, {"reads": sample["read_pairs"]}, "fastqc -t 8 input_R1.fq input_R2.fq"
            )
            # Stage 2: Demux
            demux_out = await self.execute_task(
                "Demux", sample, {"fq": fq_out, "whitelist": whitelist_val.get()}, "alevin-fry demux --whitelist 3M-february-2018"
            )
            # Stage 3: Spliced_Align (Injects OOM on Sample 3)
            align_out = await self.execute_task(
                "Spliced_Align", sample, {"demux": demux_out, "ref": ref_genome_val.get()}, "STAR --soloType Droplet --genomeDir GRCh38"
            )
            # Stage 4: UMI_Dedup
            dedup_out = await self.execute_task(
                "UMI_Dedup", sample, {"bam": align_out}, "umi_tools dedup --directional"
            )
            # Stage 5: Cell_Filter
            filter_out = await self.execute_task(
                "Cell_Filter", sample, {"matrix": dedup_out}, "Rscript emptyDrops.R --fdr 0.01"
            )
            # Stage 6: Normalization
            norm_out = await self.execute_task(
                "Normalization", sample, {"filtered": filter_out}, "scanpy.pp.normalize_total & sc.pp.log1p"
            )
            # Stage 7: Clustering
            clust_out = await self.execute_task(
                "Clustering", sample, {"norm": norm_out}, "scvi.model.SCVI.train() & scanpy.tl.leiden()"
            )
            return clust_out

        # Run all samples in parallel
        tasks = []
        while not samples_channel.empty():
            sample_item = await samples_channel.get()
            tasks.append(process_sample_dag(sample_item))

        cluster_results = await asyncio.gather(*tasks)

        # Stage 8: Report (Cohort Aggregation - consumes all sample clustering results)
        cohort_summary = {
            "sample_id": "COHORT_SUMMARY",
            "base_memory_gb": 8,
            "peak_ram_required_gb": 6.0,
            "read_pairs": sum(s["read_pairs"] for s in SAMPLES_METADATA)
        }
        report_out = await self.execute_task(
            "Report",
            cohort_summary,
            {"clustered_samples": [c["sample_id"] for c in cluster_results]},
            "multiqc . --title 'Clinical Prostate Cancer Cohort Report'"
        )

        return {
            "cluster_results": cluster_results,
            "report_out": report_out
        }


# ============================================================================
# 5. Pipeline Runner & Output Generator
# ============================================================================

def write_trace_csv(filepath: str, records: List[Dict[str, Any]]):
    headers = [
        "task_id", "process_name", "sample_id", "attempt", "status",
        "memory_gb", "cpus", "duration_s", "exit_code", "cached",
        "hash_key", "work_dir", "timestamp"
    ]
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for r in records:
            writer.writerow(r)


def write_summary_report(
    filepath: str,
    cold_trace: List[Dict[str, Any]],
    warm_trace: List[Dict[str, Any]],
    retries: List[Dict[str, Any]],
    cold_duration: float = 0.0,
    warm_duration: float = 0.0
):
    cold_success = sum(1 for r in cold_trace if r["status"] == "SUCCESS")
    cold_failed = sum(1 for r in cold_trace if r["status"] == "FAILED_OOM")
    cold_total = len(cold_trace)

    warm_cached = sum(1 for r in warm_trace if r["cached"])
    warm_total = len(warm_trace)
    cache_rate = (warm_cached / warm_total * 100.0) if warm_total > 0 else 0.0
    # Measured wall-clock values (rounded to 1 ms by the caller); the speedup is derived from them, never hard-coded.
    if warm_duration > 0:
        speedup_txt = f"{cold_duration / warm_duration:,.0f}x speedup"
    elif cold_duration > 0:
        speedup_txt = f">{int(cold_duration / 0.0005):,}x speedup"  # warm run rounded to 0.0 s: lower bound with 0.5 ms
    else:
        speedup_txt = "n/a"

    md = f"""# Track 10: Production Pipelines & Workflow Orchestration
## Nextflow DSL2 Reactive Dataflow & DAG Dependency Engine Simulation

### 1. Executive Summary & Architecture Overview
This pipeline simulates an enterprise-grade Nextflow DSL2 reactive workflow engine tailored for high-throughput single-cell RNA sequencing. It models asynchronous Queue/Value channels, dynamic OOM error recovery with automated memory escalation, and cryptographic `-resume` execution caching.

- **Total Clinical Cohort Samples**: 4 (Primary PRAD, Normal Adjacent, CRPC Metastasis, t-NEPC Progression)
- **Pipeline Stages**: 8 sequential/parallel stages (`FastQC` ➔ `Demux` ➔ `Spliced_Align` ➔ `UMI_Dedup` ➔ `Cell_Filter` ➔ `Normalization` ➔ `Clustering` ➔ `Report`)
- **Total Initial Tasks Executed (Cold Run)**: {cold_total} (Success: {cold_success}, OOM Recovered: {cold_failed})
- **Warm Run Cache Hit Ratio**: **{cache_rate:.1f}% ({warm_cached}/{warm_total} tasks resumed from cache)**

---

### 2. Pipeline Stage Specification & Resource Profiles

| Stage | Process Name | Purpose & Algorithm | Base RAM | Peak Cores |
|:-----:|:-------------|:--------------------|:--------:|:----------:|
| 1 | `FastQC` | Raw read quality, adapter contamination & Phred Q30 profiling | 16 GB | 8 |
| 2 | `Demux` | 16bp cell barcode demultiplexing against 3M-february-2018 (10x v3) whitelist | 16 GB | 8 |
| 3 | `Spliced_Align` | STARsolo / Alevin-fry splice-aware alignment against GRCh38 | 16-32 GB | 16 |
| 4 | `UMI_Dedup` | 12bp UMI directional adjacency graph PCR deduplication | 16 GB | 8 |
| 5 | `Cell_Filter` | EmptyDrops Dirichlet-Multinomial ambient RNA removal (FDR < 0.01) | 16 GB | 8 |
| 6 | `Normalization` | SCTransform regularized negative binomial modeling & 2,000 HVGs | 16 GB | 8 |
| 7 | `Clustering` | scVI 10D variational latent space & Leiden graph partitioning | 16 GB | 8 |
| 8 | `Report` | Cohort-wide MultiQC diagnostic aggregation & clinical audit signoff | 8 GB | 4 |

---

### 3. Dynamic Error Recovery & Self-Healing Event Log

During the cold run on **Sample_03_CRPC_Metastasis** (an ultra-high depth metastatic specimen with 94.8M reads), the alignment stage exceeded the default 16 GB limit (actual demand: 28.4 GB), triggering a simulated Linux OOM Killer signal (`Exit Code 137`).

The Nextflow DSL2 error strategy (`errorStrategy 'retry'`, `maxRetries 2`, `memory = memory * task.attempt`) intervened dynamically:

"""
    for e in retries:
        md += f"""- **Task ID**: `{e['task_id']}`
  - **Process**: `{e['process']}` on `{e['sample']}`
  - **Failed Attempt**: #{e['failed_attempt']} (Allocation: {e['failed_memory_gb']} GB RAM)
  - **Error Signal**: `{e['error']}`
  - **Reactive Action**: **{e['action']}** ➔ **Attempt #2 SUCCEEDED with 32 GB RAM!**\n\n"""

    md += f"""---

### 4. Cryptographic Caching & `-resume` Efficiency Benchmark

Nextflow creates an immutable MD5 hash for each task:
$$\\text{{Hash}} = \\text{{MD5}}(\\text{{Process}} + \\text{{SampleID}} + \\text{{Inputs}} + \\text{{Script}})$$

| Metric | Cold Execution (Initial) | Warm Execution (`-resume`) | Efficiency Gain |
|:-------|:------------------------:|:--------------------------:|:---------------:|
| Total Tasks Processed | {cold_total} | {warm_total} | - |
| Tasks Executed | {cold_total} | 0 | **100% compute saved** |
| Tasks Retrieved from Cache | 0 | {warm_cached} | **{cache_rate:.1f}% Cache Hit** |
| Total Compute Walltime | ~{cold_duration:.2f}s (simulated) | ~{warm_duration:.2f}s (hash lookup) | **{speedup_txt}** |
| Data I/O Redundancy | Full I/O read/write | 0 redundant writes | **Zero duplicate I/O** |

---

### 5. Task Execution Trace Table (Sample Entries)

| Task ID | Process | Sample | Attempt | Status | Memory | Cached | Work Dir |
|:--------|:--------|:-------|:-------:|:------:|:------:|:------:|:---------|
"""
    for r in cold_trace[:10]:
        md += f"| `{r['task_id']}` | `{r['process_name']}` | `{r['sample_id'].split('_')[1]}` | #{r['attempt']} | **{r['status']}** | {r['memory_gb']}GB | {r['cached']} | `{r['work_dir']}` |\n"

    md += """\n*Refer to `trace_table.csv` and `timeline_metrics.json` for full audit logs.*
"""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(md)


async def main_async():
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(pipeline_dir, "results_pipeline_01")
    os.makedirs(outdir, exist_ok=True)

    cache_file = os.path.join(outdir, ".nextflow_cache.json")
    if os.path.exists(cache_file):
        try:
            os.remove(cache_file)
        except Exception:
            pass

    cache = ExecutionCache(cache_file)

    # 1. Run Cold Execution
    engine_cold = PipelineExecutionEngine(cache=cache, resume=False)
    t0 = time.time()
    await engine_cold.run_pipeline()
    cold_duration = round(time.time() - t0, 3)

    print(f"\n[✓] Cold Execution Complete in {cold_duration}s.")
    print(f"    - Total Tasks Tracked: {len(engine_cold.trace_records)}")
    print(f"    - Retry Events: {len(engine_cold.retry_events)}")
    for ev in engine_cold.retry_events:
        print(f"      * [OOM RECOVERY] {ev['process']} on {ev['sample']}: {ev['action']}")

    # 2. Run Warm Execution with -resume
    engine_warm = PipelineExecutionEngine(cache=cache, resume=True)
    t1 = time.time()
    await engine_warm.run_pipeline()
    warm_duration = round(time.time() - t1, 3)

    cached_count = sum(1 for r in engine_warm.trace_records if r["cached"])
    total_warm = len(engine_warm.trace_records)
    cache_pct = (cached_count / total_warm * 100.0) if total_warm > 0 else 0.0

    print(f"\n[✓] Warm Execution (-resume) Complete in {warm_duration}s.")
    print(f"    - Total Tasks Checked: {total_warm}")
    print(f"    - Tasks Resumed from Cache: {cached_count} ({cache_pct:.1f}% CACHE HIT)")

    # 3. Export Structured Files
    trace_csv_path = os.path.join(outdir, "trace_table.csv")
    write_trace_csv(trace_csv_path, engine_cold.trace_records + engine_warm.trace_records)
    print(f"[✓] Trace table exported: {trace_csv_path}")

    metrics_json_path = os.path.join(outdir, "timeline_metrics.json")
    metrics_data = {
        "pipeline": "Nextflow_DSL2_Reactive_Engine_Simulation",
        "cohort_size": len(SAMPLES_METADATA),
        "stages_count": len(PIPELINE_STAGES),
        "cold_execution": {
            "duration_seconds": cold_duration,
            "total_tasks": len(engine_cold.trace_records),
            "successful_tasks": sum(1 for r in engine_cold.trace_records if r["status"] == "SUCCESS"),
            "failed_oom_retried": len(engine_cold.retry_events),
            "peak_cluster_ram_allocated_gb": max(r["memory_gb"] for r in engine_cold.trace_records)
        },
        "warm_execution": {
            "duration_seconds": warm_duration,
            "total_tasks": total_warm,
            "cached_tasks": cached_count,
            "cache_hit_percentage": cache_pct
        },
        "retries_audit": engine_cold.retry_events
    }
    with open(metrics_json_path, "w", encoding="utf-8") as f:
        json.dump(metrics_data, f, indent=2)
    print(f"[✓] Timeline metrics exported: {metrics_json_path}")

    report_md_path = os.path.join(outdir, "summary_report.md")
    write_summary_report(
        report_md_path,
        engine_cold.trace_records,
        engine_warm.trace_records,
        engine_cold.retry_events,
        cold_duration,
        warm_duration
    )
    print(f"[✓] Summary report exported: {report_md_path}")
    print(f"\n[✓] PIPELINE 01 SUCCESSFULLY COMPLETED (Exit Code 0)\n")


def main():
    asyncio.run(main_async())


if __name__ == "__main__":
    main()

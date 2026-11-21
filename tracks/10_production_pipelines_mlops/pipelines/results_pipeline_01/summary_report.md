# Track 10: Production Pipelines & Workflow Orchestration
## Nextflow DSL2 Reactive Dataflow & DAG Dependency Engine Simulation

### 1. Executive Summary & Architecture Overview
This pipeline simulates an enterprise-grade Nextflow DSL2 reactive workflow engine tailored for high-throughput single-cell RNA sequencing. It models asynchronous Queue/Value channels, dynamic OOM error recovery with automated memory escalation, and cryptographic `-resume` execution caching.

- **Total Clinical Cohort Samples**: 4 (Primary PRAD, Normal Adjacent, CRPC Metastasis, t-NEPC Progression)
- **Pipeline Stages**: 8 sequential/parallel stages (`FastQC` ➔ `Demux` ➔ `Spliced_Align` ➔ `UMI_Dedup` ➔ `Cell_Filter` ➔ `Normalization` ➔ `Clustering` ➔ `Report`)
- **Total Initial Tasks Executed (Cold Run)**: 30 (Success: 29, OOM Recovered: 1)
- **Warm Run Cache Hit Ratio**: **100.0% (29/29 tasks resumed from cache)**

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

- **Task ID**: `TASK-011`
  - **Process**: `Spliced_Align` on `Sample_03_CRPC_Metastasis`
  - **Failed Attempt**: #1 (Allocation: 16 GB RAM)
  - **Error Signal**: `Exit 137: Out of Memory (OOM Killer). Memory limit exceeded.`
  - **Reactive Action**: **DSL2 errorStrategy 'retry': Doubling memory allocation for attempt 2** ➔ **Attempt #2 SUCCEEDED with 32 GB RAM!**

---

### 4. Cryptographic Caching & `-resume` Efficiency Benchmark

Nextflow creates an immutable MD5 hash for each task:
$$\text{Hash} = \text{MD5}(\text{Process} + \text{SampleID} + \text{Inputs} + \text{Script})$$

| Metric | Cold Execution (Initial) | Warm Execution (`-resume`) | Efficiency Gain |
|:-------|:------------------------:|:--------------------------:|:---------------:|
| Total Tasks Processed | 30 | 29 | - |
| Tasks Executed | 30 | 0 | **100% compute saved** |
| Tasks Retrieved from Cache | 0 | 29 | **100.0% Cache Hit** |
| Total Compute Walltime | ~0.39s (simulated) | ~0.00s (hash lookup) | **>788x speedup** |
| Data I/O Redundancy | Full I/O read/write | 0 redundant writes | **Zero duplicate I/O** |

---

### 5. Task Execution Trace Table (Sample Entries)

| Task ID | Process | Sample | Attempt | Status | Memory | Cached | Work Dir |
|:--------|:--------|:-------|:-------:|:------:|:------:|:------:|:---------|
| `TASK-001` | `FastQC` | `01` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/d0/7f0d23` |
| `TASK-002` | `FastQC` | `02` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/4d/d3a253` |
| `TASK-003` | `FastQC` | `03` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/10/0a964d` |
| `TASK-004` | `FastQC` | `04` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/c1/f4f5ad` |
| `TASK-005` | `Demux` | `01` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/0d/b00d7c` |
| `TASK-006` | `Demux` | `02` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/0a/c73d43` |
| `TASK-007` | `Demux` | `03` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/9a/1f6476` |
| `TASK-008` | `Demux` | `04` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/f9/2a2cbf` |
| `TASK-009` | `Spliced_Align` | `01` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/27/acb7f8` |
| `TASK-010` | `Spliced_Align` | `02` | #1 | **SUCCESS** | 16GB | False | `.nextflow/work/94/860314` |

*Refer to `trace_table.csv` and `timeline_metrics.json` for full audit logs.*

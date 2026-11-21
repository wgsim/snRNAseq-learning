# 09. 대규모 단일세포 컴퓨팅 문제 해결과 모범 사례 (Troubleshooting & Best Practices)

---

## 학습 목표 (Cognitive Goals)

1. **시스템 레벨 CPU/GPU 프로파일링 도구의 원리와 실전 진단 기법 체득**: `tracemalloc` 및 `memory_profiler`를 통한 파이썬 힙 메모리 누수 추적, 그리고 NVIDIA Nsight Systems (`nsys`), Nsight Compute (`ncu`), `nvtop`을 활용한 호스트-디바이스 동기화 지연 및 SM 점유율(Occupancy) 병목을 정밀 진단할 수 있다.
2. **CUDA Out of Memory (OOM)의 수학적 근원 분석 및 다층 방어 전략 수립**: 캐시 단편화(Memory Fragmentation), 최대 순간 할당 스파이크(Peak Allocation Spikes), 통합 가상 메모리(Unified Memory, UVM) 페이징 오버헤드를 분석하고, RMM 풀링, 마이크로배칭(Micro-batching), 명시적 GC 연동을 통해 OOM을 원천 차단할 수 있다.
3. **대규모 단일세포 엔지니어링 10대 실패 양식(Top 10 Failure Modes)과 방어적 프로그래밍 패턴 습득**: 암묵적 밀집화(Implicit Dense Coercion), H5AD 청킹 불일치, 호스트-디바이스 핑퐁, 분산 직렬화 오버헤드 등 100만 세포 파이프라인 마비 요인을 체계적으로 식별하고 리팩토링할 수 있다.
4. **엔터프라이즈급 5단계 대규모 단일세포 컴퓨팅 캡스톤 프로젝트 루브릭(Capstone Project Rubric) 적용**: 100만 개 전립선암 단일세포 아틀라스의 I/O 스트리밍부터 GPU 그래프 군집화 및 프로파일링까지 이어지는 프로덕션 파이프라인의 완성도를 평가 및 검증할 수 있다.

---

## 1. 시스템 레벨 하드웨어 및 메모리 자원 프로파일링

대규모 단일세포 전사체 파이프라인에서 발생하는 버그나 성능 저하는 단순한 코드 문법 오류가 아니라, 하드웨어 계층(CPU 캐시, PCIe 버스, GPU VRAM)과 런타임 메모리 관리자(cPython GC, CUDA Driver) 간의 복잡한 상호작용에서 비롯된다.

```
[Full-Stack Single-Cell Profiling Architecture]

+-------------------------------------------------------------------------+
| User Python Layer                                                       |
|   - tracemalloc: Object-level allocation tracing & snapshot diffs       |
|   - memory_profiler: Line-by-line RAM monitoring (@profile)             |
|   - NVTX Markers: nvtx.annotate("Stage: PCA", color="blue")             |
+-------------------------------------------------------------------------+
                                    │
                                    ▼
+-------------------------------------------------------------------------+
| Hardware & System Profiling Tools                                       |
|   - nvtop: Real-time GPU SM%, VRAM, Power, PCIe throughput CLI          |
|   - Nsight Systems (nsys): Timeline tracing (CUDA Streams, CPU Wait)    |
|   - Nsight Compute (ncu): Kernel-level Roofline, Warp Occupancy, Stalls |
+-------------------------------------------------------------------------+
```

### 1.1 CPU 메모리 프로파일링: `tracemalloc` 및 순환 참조(Cyclic References) 탐지

파이썬의 메모리 관리는 참조 횟수 계산(Reference Counting)과 세대별 가비지 컬렉터(Generational GC)에 의존한다. AnnData 객체는 `.obs`, `.var`, `.obsm['X_pca']`, `.uns['neighbors']` 간의 복잡한 양방향 포인터를 포함하므로 순환 참조가 빈번히 발생한다:

```python
import tracemalloc
import gc
import linecache

def profile_memory_leak_snapshot():
    tracemalloc.start()
    snapshot1 = tracemalloc.take_snapshot()

    # 의심 파이프라인 실행
    # ... adata_chunk = process_chunk() ...

    # 명시적 GC 강제 수집
    gc.collect()
    snapshot2 = tracemalloc.take_snapshot()

    top_stats = snapshot2.compare_to(snapshot1, 'lineno')
    print("\n[Tracemalloc] Top 5 Memory Allocation Spikes:")
    for stat in top_stats[:5]:
        print(stat)
```

### 1.2 GPU 프로파일링: NVIDIA Nsight Systems 및 Nsight Compute

1. **NVIDIA Nsight Systems (`nsys`)**:
   호스트 CPU 스레드와 GPU CUDA 스트림 간의 시간축 타임라인을 가시화한다. PCIe 데이터 전송과 커널 실행이 겹치지 않고 순차적으로 차단(Blocking)되는 병목을 식별한다.
   ```bash
   nsys profile --trace=cuda,nvtx,osrt --output=prostate_1m_profile python run_gpu_pipeline.py
   ```
2. **NVIDIA Nsight Compute (`ncu`)**:
   특정 CUDA 커널(e.g., custom reduction kernel, CAGRA search)의 루프라인 분석(Roofline Model), 워프 점유율(Warp Occupancy), 공유 메모리 뱅크 충돌(Bank Conflicts)을 나노초 단위로 계측한다.
   ```bash
   ncu --set full --target-processes all -o kernel_metrics python run_cagra.py
   ```

---

## 2. CUDA Out of Memory (OOM) 근원 분석 및 다층 방어 기전

GPU VRAM(e.g., A100 80GB) 환경에서도 100만 세포 분석 도중 `torch.cuda.OutOfMemoryError` 또는 `cupy.cuda.memory.OutOfMemoryError`가 빈번히 발생한다. 이는 실제 데이터 크기가 80GB를 초과해서가 아니라 다음과 같은 구조적 요인 때문이다.

```
[Anatomy of CUDA Out of Memory (OOM)]

Total Physical GPU VRAM: 80 GB
-------------------------------------------------------------------------
[ In-Use Data: 12 GB ] [ Cached/Fragmented Pool Holes: 25 GB ] [ Reserved: 38 GB ] [ Free: 5 GB ]  (illustrative)
-------------------------------------------------------------------------
                                    ▲
                                    │
Allocation Request: 6.0 GB (e.g., for UMAP Graph Adjacency)
System Check: Total Free VRAM reported is 5 GB + 25 GB = 30 GB!
Result: CRASH possible (no usable 6.0 GB block inside the cached pools / free VA range)
```

### 2.1 메모리 단편화(Memory Fragmentation) 완화: 캐싱 할당자 풀과 RMM

GPU는 가상 주소 매핑을 쓰므로 **물리적** 단편화 자체가 큰 `cudaMalloc`을 막지는 않는다. 실무의 "단편화"는 주로 CuPy/PyTorch 같은 **프레임워크의 캐싱 할당자**가 예약해 둔 풀(세그먼트) 안에서 크기가 다양한 텐서를 반복 생성/삭제하면서 생기며, 예약은 되어 있지만 요청 크기에 맞는 연속 구간이 없어 새 4GB 블록 할당이 실패하는 형태로 나타난다.
- **완화책(보장은 아님)**: 불필요해진 캐시를 반환(`cp.get_default_memory_pool().free_all_blocks()`, `torch.cuda.empty_cache()`)하고, RMM의 `rmm.mr.CudaAsyncMemoryResource`(`cudaMallocAsync` 기반 풀) 또는 `PoolMemoryResource`로 할당자를 통일하며, PyTorch는 `expandable_segments` 옵션을 검토한다. 이들은 실패 확률을 낮출 뿐 모든 할당 성공을 보장하지 않는다.

### 2.2 통합 메모리(Unified Memory, UVM) 페이징 초과 할당(Oversubscription)

물리적 80GB VRAM을 일시적으로 초과하는 거대 연산(e.g., 1M 세포의 전체 유전자 스케일링)을 수행할 때는 시스템 호스트 RAM(512GB)을 백업 스토리지로 활용하는 UVM 오버서브스크립션을 구성한다:

```python
import rmm

# VRAM 80GB 초과 시 시스템 DDR5 RAM으로 자동 투명 페이징(Paging)
rmm.reinitialize(
    managed_memory=True,  # cudaMallocManaged 활성화
    pool_allocator=True,
    initial_pool_size="64GB"  # 문자열은 RMM이 직접 해석 (GB=10^9 bytes)
)
```

---

## 3. 대규모 단일세포 엔지니어링 10대 실패 양식 (Top 10 Failure Modes)

현대 바이오인포매틱스 엔지니어링 현장에서 100만 세포 이상의 대규모 단일세포 전사체 파이프라인을 구축할 때 마주치는 가장 치명적인 10대 실패 양식과 방어적 프로그래밍 해법은 다음과 같다.

```
[Top 10 Failure Modes Matrix]

1. Implicit Dense Conversion   ──> Sparse in-place log1p / scaling (no densify)
2. H5AD Chunk Misalignment     ──> 1D CSR row-aligned Zarr chunking
3. Host-Device Ping-Pong       ──> Complete in-device AnnData chaining
4. Scale-Up vs Out Mismatch    ──> Single-node GPU preferred over distributed CPU
5. Dask Serialization Flood    ──> 100MB-300MB chunk fusion (P = 20~50)
6. Local Normalization Drift   ──> Chan-Golub-LeVeque global reduction
7. rSVD Power Iteration Defect ──> q=2 power iteration for slow decay
8. Graph Hubness              ──> CAGRA rank-based reordering + reverse edges (fixed degree 32)
9. Louvain Disconnection       ──> Leiden 3-stage connectivity refinement
10. Uninitialized VRAM Pool    ──> Pre-allocated RMM Pool + gc.collect()
```

| 실패 양식 (Failure Mode) | 근본 메커니즘 (Root Cause) | 증상 및 오류 코드 | 방어적 엔지니어링 해법 |
| :--- | :--- | :--- | :--- |
| **1. 암묵적 밀집 변환 (Implicit Dense Conversion)** | 희소 행렬에 `X.toarray()` 호출 또는 평균 중심화(`X - mean` → 밀집 `np.matrix`) | 즉각적인 메모리 140GB 폭증 및 `SIGKILL (Exit Code 137)` | CSR `data` 벡터 직접 조작 In-place 브로드캐스팅 (`np.repeat`, 임시 버퍼 = data 크기 1배) |
| **2. H5AD 청킹 불일치 (Chunk Misalignment)** | HDF5 청크가 유전자 축으로 저장되어 세포별 읽기 시 캐시 스래싱 유발 | NVMe 3GB/s 환경에서 15MB/s로 디스크 I/O 99% 급락 | 스트리밍 크기(5만 세포)에 정렬된 Zarr 1차원 CSR 청킹 |
| **3. 호스트-디바이스 핑퐁 (Host-Device Ping-Pong)** | 매 분석 스텝마다 CPU와 GPU 간 PCIe 왕복 데이터 복사 | GPU 가속에도 불구하고 전체 실행 시간이 CPU보다 느림 | AnnData 객체 전체를 GPU VRAM에 상주시키고 체이닝 실행 |
| **4. 스케일업 vs 아웃 미스매치 (Architecture Mismatch)** | 단일 GPU 80GB로 30초면 끝날 작업을 50노드 CPU Dask 클러스터로 구축 | 네트워크 직렬화 및 노드 간 셔플로 시간/비용 100배 낭비 | 1M~5M 세포는 단일 고성능 GPU(A100/H100) 스케일업 우선 채택 |
| **5. Dask 직렬화 오버헤드 (Serialization Flood)** | 청크 크기를 너무 작게 설정(<10MB)하여 수십만 태스크 생성 | Dask Scheduler 프로세스가 GIL 락에 걸려 CPU 100% 프리징 | 태스크 퓨전 및 청크당 100MB~300MB(5만 세포 단위) 정합 |
| **6. 청크별 통계 불일치 (Chunk-wise Statistics Drift)** | 라이브러리 크기 정규화는 세포 단위라 청크와 무관하나, QC 임계값·HVG 평균/분산 등 유전자·전역 통계를 청크별로 따로 계산해 적용 | 청크마다 다른 임계값/HVG 기준으로 세포·유전자가 선택되어 청크 경계 의존적 결과 발생 | 전역 평균 및 분산의 병렬 리덕션(Chan-Golub-LeVeque)으로 전역 통계 후 적용 |
| **7. rSVD 파워 반복 결함 (Power Iteration Defect)** | 단일세포 전사체의 완만한 스펙트럼 감쇠에서 $q=0$ 단순 투영 수행 | 하위 노이즈 특이값이 고유공간을 오염시켜 클러스터 분리 실패 | 파워 반복 횟수 $q=1$ 또는 $q=2$ 강제 지정 |
| **8. kNN 그래프 허브 왜곡 (Graph Hubness)** | 브루트 포스 또는 비정규 그래프에서 고차원 허브 노드로 엣지 쏠림 | UMAP 시각화에서 모든 세포가 하나의 거대한 덩어리로 뭉침 | 순위 기반 재정렬·역방향 엣지 추가로 고정 차수 그래프 구성(CAGRA) 적용 |
| **9. Louvain 연결성 결함 / 해상도 한계** | Louvain 집계 과정에서 연결되지 않은 커뮤니티 발생 가능(Leiden은 γ-연결성 보장); 모듈러리티 해상도 한계는 Leiden도 공유 | 희귀 암세포(NEPC 0.1%)가 일반 선암 클러스터로 병합될 수 있음 | Leiden 사용 + 해상도 $\gamma$ 조정 + 마커 기반 후처리 |
| **10. RMM 미할당 단편화 (VRAM Fragmentation)** | 표준 `cudaMalloc`의 무작위 할당으로 인해 가용 VRAM 파편화 | 30GB VRAM이 남아있음에도 4GB 연속 할당 요청 시 OOM 크래시 | RMM 사전 할당 메모리 풀(Pool Allocator) 필수 바인딩 |

---

## 4. 5단계 대규모 단일세포 컴퓨팅 캡스톤 프로젝트 루브릭 (Capstone Rubric)

100만 개 단일세포 전립선암 아틀라스 데이터 파이프라인의 완성도를 정량적으로 평가하고 프로덕션 배포 적합성을 검증하기 위한 5단계 평가 루브릭은 다음과 같다.

```
[5-Stage Single-Cell Engineering Capstone Pipeline]

[ Stage 1: Out-of-Core I/O & Streaming QC ]
  ├── Zarr Chunk Alignment (50k cells / chunk)
  ├── Single-pass Chan-Golub-LeVeque Reduction
  └── Sparsity-preserving In-place Normalization
             │
             ▼
[ Stage 2: Streaming / Randomized Dimensionality Reduction ]
  ├── Implicit Centered Randomized SVD (q=2)
  ├── Subspace Grassmannian Fidelity > 99.5%
  └── Peak Memory Bounded < 8 GB
             │
             ▼
[ Stage 3: GPU-Accelerated Manifold & Embedding ]
  ├── NVIDIA CAGRA fixed-degree kNN (k=15, < 5 sec)
  ├── Warp-collaborative Coalesced Traversal
  └── cuML Asynchronous SGD UMAP 2D
             │
             ▼
[ Stage 4: High-Throughput Community Detection & Marker Discovery ]
  ├── cuGraph Parallel Leiden with Graph Coloring
  ├── Guaranteed Community Connectivity (No Louvain artifact)
  └── Fast GPU Logistic Regression Biomarker Extraction
             │
             ▼
[ Stage 5: Profiling, Memory Guarantee & Production Deployment ]
  ├── Nsight Systems Zero-Stall PCIe Validation
  ├── RMM Pool Fragmentation Defense (< 5% overhead)
  └── Reproducible Docker/Singularity Image Delivery
```

### 상세 평가 기준표 (Grading Rubric)

| 평가 단계 (Stage) | 핵심 검증 항목 (Key Validation Criteria) | 우수 (Exemplary - 5점) | 보통 (Proficient - 3점) | 미흡 (Unacceptable - 1점) |
| :--- | :--- | :--- | :--- | :--- |
| **Stage 1: Out-of-Core I/O 및 스트리밍 전처리** | Zarr 청킹 정합성, Chan-Golub-LeVeque 전역 QC, 희소성 보존 $\log(1+x)$ | 밀집 변환 없이 data 크기 이내의 임시 버퍼만 사용, 1M 세포 QC 2분 이내 완주 | 계산은 정확하나 중간 복사본 발생으로 메모리 15GB 이상 점유 | 밀집 행렬 변환으로 워커 OOM 발생 또는 순차 단순 평균 오차 |
| **Stage 2: 점진적/무작위 차원 축소** | rSVD 암묵적 중심화 연산자 구현, 파워 반복 $q=2$, Grassmannian 충실도 | 암묵적 연산자 완벽 적용, 고유공간 일치도 > 99.8%, 시간 2분 이내 | Scikit-learn IPCA 사용(메모리는 안전하나 15분 이상 소요) | $N \times G$ 밀집 행렬 직접 중심화 시도 또는 $q=0$ 노이즈 오염 |
| **Stage 3: GPU 가속 매니폴드 및 임베딩** | NVIDIA CAGRA kNN 그래프 구축, 워프 메모리 응집, cuML UMAP | kNN 그래프 5초 이내 구축, Recall@15 > 98%, UMAP 5초 이내 수렴 | CPU pynndescent 사용으로 이웃 탐색에 20분 이상 소요 | Brute-force 거리 계산 시도로 GPU VRAM OOM 발생 |
| **Stage 4: 병렬 군집화 및 바이오마커 발굴** | cuGraph Leiden 정점 색칠, 해상도 튜닝, 0.1% 희귀 NEPC 후보 탐색 | Leiden 수렴 시간 측정 기록, 해상도별 NEPC 마커 점수 집중 클러스터 확인 | Louvain 사용으로 연결되지 않은 커뮤니티가 생기거나 단일 해상도만 사용 | 멀티스레드 레이스 컨디션으로 모듈러리티 루프 진동/발산 |
| **Stage 5: 시스템 프로파일링 및 무장애 배포** | Nsight Systems 타임라인 계측, RMM 풀 단편화 방어, 재현성 패키징 | RMM 풀 완벽 바인딩, PCIe 핑퐁 0건, 타임라인 오버랩 입증 | 가속은 성공했으나 수동 메모리 해제 누락으로 재실행 시 OOM | 하드웨어 프로파일링 부재, 메모리 누수 발생, 컨테이너 미제공 |

---

## 5. 실전 트러블슈팅 파이썬 진단 모듈 구현

아래 코드는 대규모 단일세포 분석 파이프라인에서 런타임 VRAM 누수 및 힙 메모리 스파이크를 실시간으로 감시하고, OOM 임계치 도달 시 메모리를 자동으로 방어 정리하는 엔터프라이즈 모니터링 데코레이터이다.

```python
"""
Enterprise Single-Cell Pipeline Diagnostic & Memory Guard Module
Real-time VRAM/RAM Monitoring and Automatic Fragmentation Defense
"""

import functools
import gc
import time
import os
import psutil
import cupy as cp
import rmm


def single_cell_memory_guard(stage_name: str, warn_vram_gb: float = 70.0):
    """
    단일세포 파이프라인의 각 단계 전후로 CPU RAM과 GPU VRAM의 누수 및 단편화를 정밀 추적하고
    비정상적인 메모리 잔류를 강제 회수하는 방어적 데코레이터.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # 1. 실행 전 메모리 상태 스냅샷
            process = psutil.Process(os.getpid())
            ram_before_mb = process.memory_info().rss / (1024 * 1024)
            
            # GPU VRAM 상태 쿼리
            free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()
            vram_used_before_gb = (total_bytes - free_bytes) / (1024**3)
            total_vram_gb = total_bytes / (1024**3)

            print(f"\n[MemoryGuard: START] >>> Stage: {stage_name}")
            print(f"  CPU RAM Baseline : {ram_before_mb:.2f} MB")
            print(f"  GPU VRAM Used    : {vram_used_before_gb:.2f} GB / {total_vram_gb:.2f} GB")

            t0 = time.time()
            try:
                result = func(*args, **kwargs)
            except cp.cuda.memory.OutOfMemoryError as oom:
                print(f"\n[CRITICAL OOM DETECTED] Stage '{stage_name}' crashed with CUDA Out of Memory!")
                print("[Emergency Recovery] Flushing CuPy Memory Pool and forcing Python GC...")
                # 긴급 메모리 회수
                gc.collect()
                cp.get_default_memory_pool().free_all_blocks()
                raise oom

            elapsed = time.time() - t0

            # 2. 실행 후 메모리 상태 스냅샷
            ram_after_mb = process.memory_info().rss / (1024 * 1024)
            free_bytes, _ = cp.cuda.runtime.memGetInfo()
            vram_used_after_gb = (total_bytes - free_bytes) / (1024**3)
            vram_delta_gb = vram_used_after_gb - vram_used_before_gb

            print(f"[MemoryGuard: END] <<< Stage: {stage_name} (Duration: {elapsed:.2f}s)")
            print(f"  CPU RAM Delta    : {ram_after_mb - ram_before_mb:+.2f} MB (Current: {ram_after_mb:.2f} MB)")
            print(f"  GPU VRAM Delta   : {vram_delta_gb:+.2f} GB (Current: {vram_used_after_gb:.2f} GB)")

            # 3. 위험 임계치 경고 및 방어 조치
            if vram_used_after_gb >= warn_vram_gb:
                print(f"  [WARNING] VRAM utilization exceeds {warn_vram_gb} GB! Performing defensive defragmentation...")
                gc.collect()
                cp.get_default_memory_pool().free_all_blocks()
                free_bytes_post, _ = cp.cuda.runtime.memGetInfo()
                reclaimed_gb = ((total_bytes - free_bytes_post) - vram_used_after_gb)
                print(f"  [Recovery Complete] Reclaimed {abs(reclaimed_gb):.2f} GB of fragmented VRAM.")

            return result
        return wrapper
    return decorator


# 사용 예시
@single_cell_memory_guard(stage_name="cuML PCA 50 PCs", warn_vram_gb=65.0)
def execute_safe_pca(X_gpu, n_components=50):
    import cuml
    # cuML PCA의 svd_solver는 'full' | 'jacobi' | 'auto' 만 지원한다 ('randomized'는 없음; 랜덤화 SVD는 rsc.pp.pca(..., svd_solver="randomized") 사용)
    pca = cuml.PCA(n_components=n_components, svd_solver="jacobi")
    return pca.fit_transform(X_gpu)
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: GPU VRAM 가용량 허위 보고와 캐시 메모리 단편화(Memory Fragmentation) 진단
100만 개 단일세포 전사체 분석 파이프라인에서 다음과 같은 오류가 발생하였다:
`cp.cuda.runtime.memGetInfo()`를 호출했을 때 GPU의 남은 여유 메모리는 34.2 GB로 보고되었으나, 크기 4.8 GB인 kNN 그래프 인접 행렬(`cupyx.scipy.sparse.csr_matrix`)을 생성하려는 순간 즉시 `OutOfMemoryError`가 발생하며 프로세스가 중단되었다.
1. 가용 메모리가 34.2 GB나 존재함에도 불구하고 4.8 GB 할당 요청이 실패할 수 있는 이유를, 프레임워크 캐싱 할당자의 풀·세그먼트 구조 관점에서 설명하시오.
2. 이 문제를 완화하기 위해 CuPy 및 PyTorch의 캐싱 할당자(Caching Allocator) 내부 풀 비우기(`free_all_blocks()`)와 RMM 풀/Async 리소스가 어떻게 작동하는지, 그리고 이들이 성공을 보장하지 않는 이유를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. VRAM 단편화와 연속 할당(Contiguous Allocation) 실패 메커니즘**:
- `cp.cuda.runtime.memGetInfo()`가 보고하는 여유 메모리(Free Memory)는 VRAM 전체에서 할당되지 않은 모든 빈 바이트들의 단순 산술 합계(Total Sum of Free Chunks)이다.
- GPU는 가상 주소 매핑을 사용하므로 물리 메모리의 단편화만으로 큰 할당이 불가능해지지는 않는다. 그러나 CuPy/PyTorch의 캐싱 할당자는 할당 요청을 자기 풀의 **연속 구간(블록/세그먼트)** 안에서 먼저 충족시키려 하며, 프레임워크가 예약해 둔 풀이 크기가 다른 블록들로 쪼개져 있으면 새 4.8 GB 요청에 맞는 단일 블록이 없을 수 있다.
- 이전 파이프라인 단계(HVG 선별, 마스킹, 슬라이싱)에서 크기가 서로 다른 많은 텐서들이 반복적으로 할당·해제되면 이런 풀 내부 조각화가 생긴다. 또한 `memGetInfo`의 여유량은 풀에 캐시된 메모리를 "사용 중"으로 보므로 프레임워크 입장의 가용량과 다를 수 있다. 어느 요인이 실제 원인인지는 풀 통계(`mempool.used_bytes()`, `total_bytes()`)로 확인해야 한다.

**2. 메모리 풀 비우기 및 RMM 가상 주소 재매핑 해결 기전**:
- `cp.get_default_memory_pool().free_all_blocks()` 및 `torch.cuda.empty_cache()`:
  - 프레임워크가 추후 재사용을 위해 캐싱하고 있던 비활성 메모리 블록들을 드라이버로 반환해, 이후의 큰 할당이 새로 예약될 수 있게 한다.
- **RMM 풀 / `CudaAsyncMemoryResource`**:
  - 하나의 할당자(풀 또는 `cudaMallocAsync` 기반 메모리 풀)로 할당을 통일해 라이브러리별 캐시가 서로의 메모리를 점유하는 상황을 줄이고, 필요 시 풀을 키워 큰 블록을 확보한다.
  - 다만 이는 실패 확률을 낮추는 완화책이며, 풀 내부 조각화나 총 요구량 초과 상황에서는 여전히 OOM이 날 수 있다("100% 할당 성공 보장"은 성립하지 않는다).

</details>

---

### 문제 2: Dask 분산 파이프라인에서 청크 크기(Chunk Size) 미스매치로 인한 스케줄러 GIL 프리징
100만 세포 전립선암 H5AD 파일을 분석하기 위해 Dask 클러스터를 구성하면서, 세분화된 병렬성을 극대화하기 위해 청크 크기를 세포 100개 단위($100 \text{ cells} \times 35,000 \text{ genes}$, 청크당 약 560 KB 희소 크기)로 분할하였다 ($P = 10,000$개 파티션).
1. 연산 파이프라인(QC -> 정규화 -> log1p -> 유전자 분산)을 지연 그래프로 컴파일할 때, Dask 스케줄러 프로세스 내부에서 발생하는 태스크 그래프 노드 폭증과 파이썬 GIL(Global Interpreter Lock) 락 경합 메커니즘을 설명하시오.
2. 이로 인해 수십 개의 고성능 워커 노드가 연산을 수행하지 못하고 장시간 유휴(CPU 0%) 상태로 정체되는 병목을 해소하기 위한 최적 청크 크기 산정 수학적 공식을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 태스크 노드 폭증과 스케줄러 GIL 프리징 메커니즘**:
- 파티션 수가 $P = 10,000$개일 때, 5단계의 연산 체인이 존재하면 Dask DAG의 전체 태스크 노드 수는:
  $$N_{\text{tasks}} = 10,000 \text{ chunks} \times 5 \text{ stages} = 50,000 \text{ tasks}$$
- Dask의 중앙 스케줄러(`dask-scheduler`)는 단일 파이썬 프로세스 내에서 실행된다.
- 스케줄러는 50,000개 태스크의 의존성 상태 머신(Dependency State Machine), 워커 노드로의 직렬화(Cloudpickle Serialization), 네트워크 소켓 I/O 이벤트 루프(Tornado)를 단일 스레드로 추적해야 한다.
- 개별 청크의 크기가 560 KB에 불과하여 워커 노드의 실제 연산 수행 시간(수 밀리초)이 스케줄러의 태스크 디스패치 및 네트워크 패킷 처리 시간(수십 밀리초)보다 훨씬 짧아진다.
- 결과적으로 스케줄러 프로세스는 파이썬 인터프리터의 GIL 락에 갇혀 CPU 100%를 소모하며 병목이 되고, 워커 노드들은 스케줄러로부터 다음 일감을 받지 못해 네트워크 I/O 대기 상태에서 유휴(Idle, CPU 0%)로 멈춰 서게 된다.

**2. 최적 청크 크기 산정 수학적 가이드라인**:
태스크 스케줄링 오버헤드 $T_{\text{sched}}$(태스크당 약 1ms)가 실제 연산 시간 $T_{\text{compute}}$의 1% 미만이 되도록 보장해야 한다:

$$\frac{T_{\text{sched}}}{T_{\text{compute}}} < 0.01 \implies T_{\text{compute}} > 100 \times T_{\text{sched}} \approx 100 \text{ ms} = 0.1 \text{ 초}$$

희소 행렬의 I/O 대역폭 및 CPU 벡터 연산 속도를 고려할 때, 단일 청크가 최소 0.5초~2초 동안 지속적으로 연산될 수 있는 크기여야 한다.
비압축 희소 CSR 기준으로 최적의 청크 메모리 크기 $S_{\text{opt}}$는:

$$S_{\text{opt}} \in [100 \text{ MB}, \quad 300 \text{ MB}]$$

100만 세포($10^6 \times 35,000$, 총 희소 크기 5.6 GB)의 경우:
- 최적 파티션 수: $P_{\text{opt}} = \frac{5.6 \text{ GB}}{280 \text{ MB}} \approx 20 \sim 30 \text{ 개}$
- 세포 축 청크 크기: $n_{\text{chunk}} = \frac{1,000,000}{20} = 50,000 \text{ cells}$
청크를 50,000개 단위로 병합하면 태스크 수가 50,000개에서 단 100개로 급감하여 스케줄러 오버헤드가 0%로 소멸하고 워커의 연산 가동률이 100%로 복원된다.

</details>

---

### 문제 3: Nsight Systems 타임라인 분석을 통한 호스트-디바이스 동기화 지연(Host-Device Sync Stall) 진단
NVIDIA A100 GPU 상에서 실행한 100만 세포 파이프라인의 Nsight Systems (`nsys`) 프로파일링 결과 타임라인에서 다음과 같은 병목 현상이 관측되었다:
GPU의 커널 실행 바(Kernel Execution Bar) 사이사이에 수 밀리초에 달하는 빈 공백(Idle Gaps)이 지속적으로 나타나고 있으며, CPU 타임라인에서는 `cudaStreamSynchronize()` 및 `cudaMemcpy`가 메인 스레드를 반복적으로 차단하고 있다.
1. 파이썬 코드 상에서 사용자가 무심코 호출한 어떤 관용적 표현(Idiomatic expressions, e.g., `.item()`, `print(tensor.shape)`, `len()`, CPU 조건문)이 GPU 비동기 파이프라인에 암묵적 강제 동기화(Implicit Hard Synchronization)를 유발하는지 분석하시오.
2. GPU 연산 커널과 PCIe 데이터 스트리밍을 완벽하게 오버랩(Overlap)시키기 위해 다중 CUDA 스트림(`cp.cuda.Stream`)과 고정 호스트 메모리(Pinned Host Memory, `cudaHostAlloc`)를 적용하는 파이프라인 리팩토링 원리를 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 암묵적 강제 동기화를 유발하는 파이썬 안티패턴**:
- CUDA 연산(CuPy, PyTorch, cuML)은 기본적으로 비동기(Asynchronous) 방식으로 동작한다. 즉, CPU가 GPU 드라이버 큐에 커널 실행 명령을 넣은 후 커널 완료를 기다리지 않고 다음 파이썬 코드로 즉시 진행한다.
- 그러나 파이썬 코드 내에서 GPU에 상주하는 텐서의 **스칼라 값이나 메타데이터를 확인하려는 순간**, CPU는 연산 결과를 알아야 분기를 진행할 수 있으므로 GPU가 해당 시점까지의 모든 작업을 완료할 때까지 메인 스레드를 멈추고 기다리는 **하드 동기화(Hard Synchronization Stall)**가 강제된다:
  1. `if (cp.sum(mask) > 0):` — GPU 리덕션 결과를 CPU 불리언으로 평가하기 위해 동기화 발생.
  2. `print(f"Current loss: {loss.item()}")` 또는 `float(metric)` — 스칼라 변환을 위한 D2H 복사 및 동기화.
  3. `adata.obs["qc"] = cp.asnumpy(gpu_array)` — 매 스텝마다 넘파이 변환 호출.
- 루프 내부에서 이러한 코드가 수백 번 실행되면, GPU는 연산 1ms 수행 후 동기화 대기로 5ms 동안 멈추는 극심한 지연 공백(Idle Gap)이 발생하여 GPU 점유율이 20% 미만으로 추락한다.

**2. CUDA 스트림 및 Pinned Memory를 통한 제로 스톨(Zero-Stall) 리팩토링**:
1. **페이지 고정 메모리 (Pinned / Page-locked Host Memory)**:
   - 일반 호스트 메모리는 OS가 디스크로 페이징할 수 있으므로 GPU DMA 컨트롤러가 직접 접근하지 못하고 임시 버퍼로 복사된 후 전송된다.
   - `cp.cuda.alloc_pinned_memory()`로 할당된 Pinned Memory를 사용하면 CPU 개입 없이 GPU DMA 컨트롤러가 PCIe 최대 대역폭(31.5 GB/s)으로 직접 데이터를 전송한다.
2. **다중 CUDA 스트림을 통한 전송과 연산의 오버랩 (Compute-Transfer Overlap)**:
   - 독립적인 두 개의 CUDA 스트림(`stream_io`, `stream_compute`)을 생성한다.
   - 청크 $i+1$의 Pinned Memory -> GPU VRAM 복사(`cudaMemcpyAsync`)를 `stream_io`에서 실행하는 동시에, 이전 청크 $i$의 PCA/정규화 연산 커널을 `stream_compute`에서 동시에 실행한다.
   - 이를 통해 데이터 전송 지연 시간이 연산 시간 뒤로 완전히 숨겨져(Latency Hiding), GPU 커널 실행 바 사이의 간극이 0이 되는 완전 포화 파이프라인이 완성된다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Chan TF, Golub GH, Leveque RJ. Algorithms for Computing the Sample Variance: Analysis and Recommendations. *The American Statistician*. 1983;37(3):242-247. doi: 10.1080/00031305.1983.10483115
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- cuML (github.com/rapidsai/cuml).
- CuPy (github.com/cupy/cupy).
- psutil 7.2.2 (github.com/giampaolo/psutil).
- RAPIDS RMM (github.com/rapidsai/rmm).
<!-- REFERENCES:END -->

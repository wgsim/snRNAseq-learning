# 02. Out-of-Core 연산과 HDF5 AnnData Backed 모드의 내부 구조 (Out-of-Core Computation and Internal Architecture of HDF5 AnnData Backed Mode)

---

## 학습 목표 (Cognitive Goals)

1. **Out-of-Core 메모리 매핑 아키텍처 원리 이해**: 워크스테이션의 물리 주기억장치(RAM) 크기를 초과하는 초대형 단일세포 전사체 데이터를 가상 메모리 및 디스크 블록 단위 스트리밍을 통해 처리하는 Out-of-Core 연산의 개념과 철학을 정립할 수 있다.
2. **HDF5 (`.h5ad`) 바이너리 컨테이너의 내부 계층 구조 해부**: HDF5 포맷의 Groups, Datasets, Attributes 계층 구조와 `.h5ad` 명세에서 희소 행렬(`data`, `indices`, `indptr`)이 저장되는 방식 및 B-Tree 기반 청크 인덱싱 메커니즘을 상세히 분석할 수 있다.
3. **`anndata.read_h5ad(path, backed='r')` 지연 로딩(Lazy Evaluation) 분석**: Backed 모드 호출 시 발생하는 파일 디스크립터 바인딩, 슬라이싱 시점의 부분 데이터 구체화(Materialization), 그리고 인메모리(In-Memory) 모드 대비 메모리 풋프린트 격차를 규명할 수 있다.
4. **HDF5 전역 파일 잠금(Global File Lock / POSIX Advisory Lock)과 멀티프로세스 병목 진단**: 딥러닝 DataLoader 및 멀티프로세싱 환경에서 HDF5 쓰기 잠금 충돌과 `fork()` 상속 핸들 문제를 구분해 진단하고, 이를 방지하기 위한 배치 스트리밍 최적화 기법을 구현할 수 있다.
5. **랜덤 셀 접근 vs 순차 청크 스트리밍 간 Throughput 격차(밀집 2D 청크 가정 시 최대 약 1,000배) 수리적 검증**: 저장 매체(NVMe SSD)의 IOPS 및 지연 시간 한계에 따른 랜덤 셀 페칭(Random Fetching)의 읽기 증폭(Read Amplification)을 규명하고 고속 스트리밍 제너레이터를 설계할 수 있다.

---

## 1. Out-of-Core 연산과 메모리 매핑(mmap) 철학

전통적인 데이터 분석 파이프라인은 모든 데이터를 주기억장치(DRAM)에 사전 로드(Eager Loading)하는 것을 전제로 설계되었다. 그러나 100만 개($10^6$) 이상의 세포를 포함하는 인간 세포 지도(HCA) 데이터셋의 경우, 전처리 및 중간 결과 버퍼를 고려할 때 수백 기가바이트의 메모리가 요구된다.

```
[In-Memory vs. Out-of-Core Processing Paradigm]

1. In-Memory (Eager Loading):
   Disk (.h5ad) ──[Full Read: 60 GB]──► DRAM (100% Occupied) ──► CPU Cache
   * Risk: OOM-Killer termination if dataset > physical RAM.

2. Out-of-Core (Lazy Chunk Streaming / Backed Mode):
   Disk (.h5ad) ──[On-Demand Slice: 50 MB]──► DRAM Buffer ──► CPU Cache
   * Advantage: Constant memory footprint O(Chunk Size), handling 100M+ cells.
```

**Out-of-Core 연산**은 전체 데이터셋 중 현재 CPU/GPU가 연산에 즉시 필요로 하는 국소 슬라이스(Slice)만을 영구 저장 매체(NVMe SSD)에서 DRAM으로 온디맨드(On-Demand) 페칭하고, 연산이 종료되면 메모리를 즉시 반환하거나 덮어쓰는 가상화 기법이다. 이를 통해 시스템은 물리 RAM 용량의 제약에서 벗어나 스토리지 용량($\sim \text{수 TB}$) 한계까지 선형적으로 확장할 수 있다.

---

## 2. HDF5 (`.h5ad`) 바이너리 컨테이너의 내부 구조

`AnnData`의 표준 영구 저장 포맷인 `.h5ad` 파일은 **HDF5 (Hierarchical Data Format version 5)** 사양을 기반으로 구축된 자체 완비형(Self-describing) 바이너리 컨테이너이다.

### 2.1 HDF5 트리 계층 구조 (Group, Dataset, Attribute)

HDF5는 파일 내부에 유닉스 파일 시스템과 유사한 계층적 디렉토리 트리를 구성한다:

```
example_atlas.h5ad (HDF5 Root Group "/")
├── X (Group or Dataset: 발현 행렬)
│   ├── data (1D Dataset, shape: [nnz], dtype: float32)
│   ├── indices (1D Dataset, shape: [nnz], dtype: int32 or int64)
│   ├── indptr (1D Dataset, shape: [N + 1], dtype: int64)
│   └── (Attributes: encoding-type="csr_matrix", shape=[N, G])
├── obs (Group: 세포 메타데이터)
│   ├── _index (1D Dataset: 세포 바코드 식별자)
│   ├── cell_type (1D Dataset: 범주형 인덱스)
│   ├── batch (1D Dataset: 배치 레이블)
│   └── (Attributes: encoding-type="dataframe")
├── var (Group: 유전자 메타데이터)
│   ├── _index (1D Dataset: Ensembl/Gene Symbol)
│   ├── highly_variable (1D Dataset: bool)
│   └── (Attributes: encoding-type="dataframe")
├── obsm (Group: 다차원 임베딩 - e.g., X_pca, X_umap)
└── uns (Group: 비정형 메타데이터, 색상 코드, 계통수)
```

### 2.2 희소 행렬의 인코딩 명세

`adata.X`가 희소 행렬(CSR)일 경우, HDF5 루트 아래의 `/X`는 단일 데이터셋이 아니라 **그룹(Group)**으로 생성되며, 내부에는 CSR의 3대 핵심 배열인 `data`, `indices`, `indptr`가 1차원 데이터셋으로 분할 저장된다. HDF5 객체의 메타데이터 속성(Attribute)에는 다음과 같은 정보가 HDF5 문자열/배열 속성으로 기록된다(JSON 문서가 아님):
*   `encoding-type`: `"csr_matrix"`
*   `encoding-version`: `"0.1.0"`
*   `shape`: `[1000000, 30000]`

### 2.3 HDF5 B-Tree 인덱싱과 청킹(Chunking) 메커니즘

HDF5 내부의 대용량 데이터셋은 연속된 단일 바이트 스트림(Contiguous layout)으로 저장되거나, 고정된 크기의 **청크(Chunk)** 단위로 분할(Chunked layout)되어 저장된다.
*   청크 단위 저장은 각 청크 블록마다 독립적인 압축(GZIP, LZF 등)을 허용한다.
*   HDF5 라이브러리는 파일 헤더 내부에 **B-Tree (또는 B+ Tree)** 인덱스를 유지하여, 임의의 다차원 좌표 슬라이스가 주어졌을 때 디스크 상의 어떤 물리적 청크 바이트 오프셋에 접근해야 하는지 $O(\log K)$ 시간 복잡도로 탐색한다.

---

## 3. `anndata.read_h5ad(path, backed='r')`의 내부 동작 원리

인메모리 로딩(`backed=None`)과 Out-of-Core Backed 로딩(`backed='r'`)의 시스템 수준 동작 방식은 근본적으로 다르다.

```python
import anndata as ad

# 1. Backed 모드로 100만 세포 아틀라스 열기
adata = ad.read_h5ad("massive_atlas.h5ad", backed="r")
```

### 3.1 지연 바인딩(Lazy Binding) 메커니즘

`backed='r'`(Read-Only Backed Mode)을 지정하면:
1.  **DRAM 소비 제로**: 발현 행렬 `X`의 원시 데이터(수십 GB)를 RAM으로 복사하지 않는다.
2.  **파일 핸들 오픈**: 운영체제 시스템 콜(`open()`)을 통해 `.h5ad` 파일의 파일 디스크립터(File Descriptor)를 획득하고, HDF5 파일 슈퍼블록(Superblock) 및 메타데이터 트리만 메모리에 캐싱한다. `obs`/`var` 데이터프레임은 열 때 읽히므로 소요 시간이 `obs` 크기에 비례한다(100만 행·단순 3개 열 `obs`의 `read_h5ad(backed='r')`가 약 1.4초로 실측됨; 수백 ms~수 초 규모).
3.  **가상 배열 래핑**: `adata.X`는 `scipy.sparse.csr_matrix`가 아니라 backed 희소 래퍼 객체(현재 `anndata.abc.CSRDataset`/`CSCDataset`; 구버전의 `SparseDataset`는 제거됨)로 바인딩되며, 이는 내부의 `h5py.Dataset` 포인터를 유지한다.

### 3.2 온디맨드 구체화(Materialization)

데이터의 실제 디스크 I/O는 사용자가 명시적인 인덱싱이나 연산을 요구하는 순간에만 발생한다:

```python
# 1,000개 세포의 서브셋 슬라이싱 (On-Demand Materialization)
subset = adata[0:1000, :].to_memory()
```

이 슬라이싱 구문이 실행되는 순간:
1.  HDF5 라이브러리는 `/X/indptr` 데이터셋에서 인덱스 $0$부터 $1000$까지의 포인터를 읽어 $1,000$개 세포의 시작/끝 바이트 오프셋을 계산한다.
2.  해당 범위에 속하는 `indices`와 `data` 디스크 블록만을 부분 읽기(Partial Read)한다.
3.  읽어온 1차원 버퍼를 결합하여 단 $1,000$개 세포 크기의 `scipy.sparse.csr_matrix`를 메모리에 인스턴스화한다. 메모리 점유율은 수십 메가바이트 수준에 머문다.

---

## 4. HDF5 파일 잠금과 멀티프로세스 사용 시 문제

HDF5 Backed 모드는 단일 스레드 순차 분석에서는 탁월한 메모리 절감 효과를 발휘하지만, **멀티프로세스 파이프라인(PyTorch DataLoader, Ray, Dask Multiprocessing)과 결합할 때 파일 핸들 공유와 잠금 설정에 주의**해야 한다.

```
[Two Failure Causes of Backed HDF5 in Multi-Worker DataLoaders]

   (A) Write-lock conflict                (B) Handle inherited by fork()
   Process W: opens file for writing      Parent: opens .h5ad (backed='r')
        │  (exclusive lock held)               │ fork()
        ▼                                      ▼
   Process R: opens read-only            Worker 0..N: share the parent's
        │  (shared lock requested)        open HDF5/h5py handle
        ▼                                      ▼
   "unable to lock file, errno = 11"     Undefined behavior / read errors
   (EAGAIN: lock refused)                Fix: each worker opens its own handle
```

HDF5는 파일을 열 때 접근 플래그에 따라 잠금을 건다. 쓰기(`H5F_ACC_RDWR`)로 열면 **배타적 잠금**, 그 외(읽기 전용)에는 **공유 잠금**을 사용한다(HDF5 File Locking 문서). 따라서 **읽기 전용으로 여는 여러 프로세스는 서로 잠금 때문에 직렬화되지 않는다.** 아래 두 가지가 실제 원인이다.

### 4.1 원인 A: 쓰기 잠금 충돌 (`unable to lock file, errno = 11`)

다른 프로세스가 같은 파일을 쓰기 모드로 열고 있으면(예: `r+`로 열린 `adata`, 또는 파일을 갱신하는 작업) 읽기 전용 열기가 `unable to lock file, errno = 11`(`EAGAIN`)로 거부될 수 있다. 읽기 전용 접근만 있는 상황이거나 단일 프로세스가 쓰기로 여는 상황이라면 `HDF5_USE_FILE_LOCKING=FALSE`로 잠금을 끌 수 있다고 HDF5 문서는 설명하지만, 쓰기 프로세스와 다른 프로세스가 동시에 접근하면 파일 손상 위험이 있다. HDF5 라이브러리의 전역 뮤텍스(`--enable-threadsafe` 빌드)는 **한 프로세스 내부의 스레드 간** 직렬화이며 서로 다른 프로세스는 이를 공유하지 않는다.

### 4.2 원인 B: `fork()`로 상속된 열린 핸들

PyTorch `DataLoader(num_workers=8)`는 `fork()`로 자식 프로세스를 만든다. 부모가 이미 연 `adata`(h5py/HDF5 핸들, 내부 상태 포함)를 자식들이 상속해 공유하면 예측할 수 없는 오류(잘못된 읽기, 세그멘테이션 폴트 등)가 발생할 수 있다. h5py 문서는 멀티프로세싱에서 **각 프로세스가 파일을 독립적으로 열 것**을 권고한다. 해법은 `worker_init_fn`(또는 첫 `__getitem__`)에서 워커마다 파일을 여는 것이다. 이때 워커 수만큼의 독립 핸들이 생기므로 I/O 대기나 GPU 활용률 저하가 있다면 원인은 잠금보다는 랜덤 접근 패턴(5절)·압축 청크 읽기 증폭·스토리지 처리량에서 먼저 찾아야 한다(구체 수치는 환경에 따라 다르며 본 문서에서 측정한 값이 아님).

---

## 5. 랜덤 셀 접근 vs 순차 청크 스트리밍 벤치마크 (Throughput 격차)

대규모 Backed HDF5 데이터셋에서 데이터를 추출할 때, **인덱싱 순서(Access Pattern)**는 처리량을 크게 좌우한다(격차 배수는 청크 크기·압축·스토리지에 의존하며, 아래 1,000배는 밀집 2D 청크 가정의 최악 사례 산술이다).

### 5.1 랜덤 셀 페칭(Random Fetching)의 읽기 증폭(Read Amplification)

딥러닝의 `shuffle=True` 옵션에 따라 임의의 무작위 인덱스 순서로 단일 세포를 하나씩 조회할 경우:
1.  **디스크 탐색 오버헤드**: NVMe SSD라 할지라도 무작위 4K 읽기 요청이 폭증하여 IOPS 한계에 도달한다.
2.  **청크 압축 해제의 낭비**: 밀집(dense) 2D 데이터셋이 예를 들어 $1,000 \times 1,000$ 크기의 청크로 압축 저장되어 있다면, 단 1개 세포의 단 1개 원소를 읽기 위해 $1,000 \times 1,000$ 크기의 압축 청크 전체를 디스크에서 읽어 압축 해제해야 한다. 이로 인해 유효 데이터 4바이트를 얻기 위해 수 메가바이트를 읽는 **극단적 읽기 증폭(Read Amplification > 1,000x)**이 발생한다. 단, 희소 H5AD의 `X/data`·`X/indices`·`X/indptr`는 **1차원** 배열이라 청크가 nnz 원소 축으로만 나뉜다(anndata 0.13.4로 5000×2000, 5% 희소 행렬을 쓰면 `X/data`·`X/indices` 청크가 `(3907,)`의 1차원으로 기록됨을 확인; 크기는 h5py 자동 추정에 따라 데이터 크기마다 다름). 이 경우 한 세포(행)는 `indptr`로 정해진 연속 구간(평균 수백~수천 원소)이어서 1~2개의 청크만 읽으면 되며, 증폭은 대략 (청크 원소 수 ÷ 행당 nnz)배 수준(예: 3,907 ÷ 100 ≈ 40배, 7,813 ÷ 1,500 ≈ 5배)에 그친다. 랜덤 접근의 비용은 주로 세포마다 별도의 `indptr` 조회와 청크 읽기·압축 해제 호출이 발생하는 점에서 나온다.

### 5.2 순차 청크 스트리밍(Sequential Chunk Streaming)

반면, 데이터를 디스크에 저장된 물리적 순서(Physical Layout)대로 일정 크기(예: 5,000~10,000개 세포)의 슬라이스로 연속 추출할 경우:
1.  HDF5 청크가 단 한 번 압축 해제되어 수천 개의 세포 데이터가 메모리로 직행한다.
2.  운영체제의 선행 읽기(Read-Ahead / Prefetching) 버퍼가 100% 작동하여 디스크 대역폭을 최대 $3,000 \sim 5,000\text{ MB/s}$까지 완전히 활용한다.

---

## 6. 실전 재현 코드: 고성능 Out-of-Core 배치 제너레이터 구현

다음은 HDF5 Backed 모드의 단점을 극복하고, 메모리 사용량을 엄격히 제한하면서 최대의 I/O 처리량을 달성하는 프로덕션급 청크 스트리밍 제너레이터 구현이다:

```python
import time
import numpy as np
import anndata as ad
import scipy.sparse as sp

class BackedChunkStreamer:
    """
    HDF5 Backed AnnData 객체로부터 I/O 병목 없이
    순차적 고속 배치 데이터를 스트리밍하는 Out-of-Core 제너레이터
    """
    def __init__(self, h5ad_path: str, chunk_size: int = 10_000, shuffle_buffer: int = 20_000):
        self.h5ad_path = h5ad_path
        self.chunk_size = chunk_size
        self.shuffle_buffer_size = shuffle_buffer
        # Backed 모드로 가상 오픈
        self.adata = ad.read_h5ad(h5ad_path, backed='r')
        self.n_cells = self.adata.n_obs
        self.n_vars = self.adata.n_vars
        
    def stream_batches(self, batch_size: int = 256):
        """
        대용량 청크를 디스크에서 순차 로드한 후, 인메모리 버퍼 내에서만
        미니 배치를 무작위 추출하여 I/O 증폭을 방지함.
        """
        indices = np.arange(self.n_cells)
        # 전체를 뒤섞지 않고 청크 단위 순회
        for chunk_start in range(0, self.n_cells, self.chunk_size):
            chunk_end = min(chunk_start + self.chunk_size, self.n_cells)
            
            # 1. 디스크에서 연속된 대용량 청크 슬라이스 로딩 (Sequential I/O)
            # to_memory()를 호출하여 해당 청크만 RAM의 순수 scipy.sparse.csr_matrix로 구체화
            chunk_adata = self.adata[chunk_start:chunk_end, :].to_memory()
            chunk_X = chunk_adata.X
            
            # 2. 청크 내부 로컬 셔플링 (In-Memory Buffer Shuffling)
            local_indices = np.random.permutation(chunk_end - chunk_start)
            
            # 3. 미니 배치 단위로 제너레이팅
            for batch_start in range(0, len(local_indices), batch_size):
                batch_end = min(batch_start + batch_size, len(local_indices))
                batch_sub_idx = local_indices[batch_start:batch_end]
                
                # 메모리 상의 CSR에서 미니 배치 추출 (완전 연속 행 슬라이싱)
                yield chunk_X[batch_sub_idx, :]
                
    def close(self):
        # Backed 파일 핸들 명시적 해제
        if hasattr(self.adata, 'file') and self.adata.file is not None:
            self.adata.file.close()

# 벤치마크 및 검증 시뮬레이션
if __name__ == "__main__":
    print("[HDF5 Backed Streaming Generator 초기화]")
    # 실제 환경에서는 대규모 H5AD 파일 경로 전달
    # streamer = BackedChunkStreamer("massive_atlas.h5ad", chunk_size=5000)
    # for batch in streamer.stream_batches(batch_size=256):
    #     pass
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: HDF5 Locking 환경에서 PyTorch DataLoader 데드락 해결
연구팀이 150만 세포 크기의 `.h5ad` 파일로 단일세포 전사체 오토인코더(Autoencoder)를 학습시키기 위해 PyTorch `DataLoader(num_workers=16)`를 가동하였다. 그러나 프로세스 시작 직후 다음과 같은 치명적 예외가 간헐적으로 발생하며 훈련이 중단되었다:
`BlockingIOError: [Errno 11] Unable to create file (unable to lock file, errno = 11, error message = 'Resource temporarily unavailable')`
1. 이 에러가 운영체제 시스템 콜(`fcntl` / `flock`) 및 HDF5 파일 락 수준에서 왜 발생하는지 규명하시오.
2. 환경 변수 `HDF5_USE_FILE_LOCKING=FALSE` 설정이 가져올 수 있는 잠재적 데이터 무결성(Data Integrity) 위험을 설명하고, 멀티프로세싱 환경에서 안전하게 HDF5 데이터를 병렬 로딩하기 위한 `worker_init_fn` 기반의 아키텍처적 해법을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 에러 발생 메커니즘 (쓰기 잠금 충돌)**:
*   HDF5 라이브러리는 파일을 열 때 접근 플래그에 따라 잠금을 건다: 쓰기로 열면 배타적 잠금, 읽기 전용이면 공유 잠금(HDF5 File Locking 문서). 읽기 전용으로 연 여러 프로세스는 서로를 막지 않는다.
*   `unable to lock file, errno = 11`(`EAGAIN`, Resource temporarily unavailable)은 **다른 프로세스(또는 같은 파일을 쓰기로 연 핸들)가 배타적 잠금을 쥐고 있을 때** 발생한다. 예: 부모가 `r+`/쓰기로 연 `adata`가 있거나, 파일을 갱신하는 작업이 동시에 돈다.
*   별개로, 부모가 열어 둔 핸들을 `fork()`된 16개 자식이 상속해 공유하면 잠금과 무관한 오동작(잘못된 읽기, 세그멘테이션 폴트)이 생길 수 있다. 요청 수가 "락 테이블 큐를 초과"해서 거부된다는 설명은 근거가 없어 삭제했다.

**2. 잠재적 위험 및 안전한 worker_init_fn 아키텍처 해법**:
*   **환경 변수 강제 비활성화의 위험**: `HDF5_USE_FILE_LOCKING=FALSE`를 적용하면 파일 락 검사는 건너뛰지만, 여러 프로세스가 동시에 파일을 수정(`r+` 모드)하거나 쓰기 작업이 일어날 경우 파일 내부 B-Tree 인덱스가 파괴되어 파일 전체가 복구 불가능하게 손상된다. 읽기 전용 접근만 있는 경우에는 HDF5 문서도 잠금을 끌 수 있다고 설명하지만, `fork()`로 상속한 핸들을 공유하는 문제는 잠금 설정과 별개로 남는다.
*   **안전한 아키텍처 해법 (`worker_init_fn` 지연 초기화)**:
    *   부모 프로세스에서는 HDF5 파일을 절대 열지 않고 메타데이터(세포 수 $N$)만 전달한다.
    *   PyTorch의 `worker_init_fn(worker_id)` 내부에서 각 자식 워커가 독립적인 파일 핸들을 독립적으로 `read_h5ad`로 열도록 지연 오픈(Lazy Open)한다.
    *   각 워커 프로세스가 독립적인 핸들과 고유 파일 디스크립터를 소유하게 되므로 상속 핸들 공유 문제가 피해진다(h5py 권고: 각 프로세스에서 독립적으로 열기). 쓰기 잠금 충돌은 쓰기 프로세스와 읽기 프로세스를 분리해야 해소된다.

</details>

---

### 문제 2: HDF5 Backed 모드에서의 Random Slicing vs Contiguous Slicing I/O 계산
100만 세포($N=10^6$), 30,000 유전자($G=30,000$) 데이터셋이 저장된 H5AD 파일이 있다. `X`가 **밀집(dense) 행렬**로 저장되어 있고 $1,000 \times 30,000$ 단위의 청크 크기(Chunk shape)로 GZIP 압축되어 있다고 가정한다(희소 CSR의 `data`/`indices`는 1차원 배열이므로 세포 단위 2D 청크가 존재하지 않는다).
어떤 알고리즘이 임의로 추출된 1,000개의 무작위 세포(Random Non-contiguous Cells)를 하나씩 반복문으로 호출(`adata[idx, :].X`)하여 데이터를 읽으려 한다.
1. 최악의 경우(Worst-case Scenario) 이 연산이 디스크에서 읽어 압축 해제해야 하는 이론적 HDF5 청크의 총 개수를 산출하시오.
2. 반면, 세포들이 디스크 상에 연속된 구간(`adata[0:1000, :].X`)으로 정렬되어 있을 때 읽어야 하는 청크의 개수를 산출하고, 두 방식 간의 **읽기 증폭률(Read Amplification Factor)**을 비교하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 무작위 세포 접근 시의 청크 압축 해제 수 (Worst-case)**:
*   전체 세포 수는 $1,000,000$개이며 청크는 $1,000$개 세포 단위로 나뉘어 있으므로 총 청크 수는 $1,000,000 / 1,000 = 1,000$개 청크이다.
*   무작위로 선택된 1,000개의 세포가 서로 다른 1,000개의 청크에 각각 1개씩 분산되어 있을 경우, 최악의 경우 **1,000개의 서로 다른 청크 전체**를 디스크에서 읽어 압축 해제해야 한다.
*   유효한 1,000개 세포의 데이터를 얻기 위해 $1,000 \times 1,000 = 1,000,000$개 세포 분량의 데이터를 압축 해제하는 파국이 발생한다.

**2. 연속 구간 접근 시의 청크 수 및 읽기 증폭률**:
*   연속된 구간 `adata[0:1000, :]`은 정확히 **단 1개의 청크(Chunk 0)** 내부에 완벽히 포함된다.
*   따라서 단 1개의 청크만 디스크에서 읽어 압축 해제하면 1,000개 세포 데이터 로딩이 완결된다.
*   **읽기 증폭률 (Read Amplification Factor)**:
    $$\text{Amplification} = \frac{1,000 \text{ chunks}}{1 \text{ chunk}} = 1,000 \times$$
*   밀집 2D 청크 가정에서 랜덤 접근은 연속 접근 대비 디스크 I/O 및 CPU 압축 해제 부하를 최대 **1,000배 증폭**시킨다(희소 CSR의 1차원 청크에서는 이보다 훨씬 작다). 이것이 Backed 모드에서 무작위 셔플링 대신 청크 버퍼 기반 스트리밍을 사용해야만 하는 수학적 이유이다.

</details>

---

### 문제 3: HDF5 SWMR (Single-Writer Multi-Reader) 모드의 한계와 실시간 파이프라인
대규모 유전체 분석 센터에서 단일세포 전사체 시퀀싱 데이터가 시퀀서에서 실시간으로 산출되는 동안, 다른 분석 프로세스들이 동시에 품질 관리(QC) 메트릭을 계산할 수 있는 파이프라인을 구축하고자 한다. 엔지니어가 HDF5의 SWMR (Single-Writer Multi-Reader) 기능을 적용하여 단일 `.h5ad` 파일에 실시간 쓰기와 동시 읽기를 수행하려 한다.
1. HDF5의 SWMR 모드가 정상 작동하기 위해 만족해야 하는 파일 시스템(POSIX) 및 데이터셋 구조상의 기술적 전제 조건을 기술하시오.
2. 희소 행렬(CSR)의 `indices`, `data`, `indptr` 데이터셋을 동적으로 확장(Append)할 때 SWMR 모드에서 발생할 수 있는 데이터 불일치(Data Inconsistency) 위험을 분석하고, 왜 산업계가 실시간 스트리밍 아키텍처로 Zarr나 TileDB를 채택하는지 논하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. HDF5 SWMR 모드의 기술적 전제 조건**:
*   **파일 시스템 전제 조건**: POSIX 호환 파일 시스템이어야 하며, 원자적 쓰기(Atomic write)와 운영체제 버퍼 플러시(`fsync()`)가 완벽히 보장되어야 한다. 분산 객체 저장소(AWS S3, NFS 등)에서는 SWMR이 동작하지 않는다.
*   **데이터셋 구조 전제 조건**: SWMR 모드에서는 새로운 Group이나 Dataset을 동적으로 생성할 수 없으며, 사전에 정의된 고정 차원 데이터셋의 크기만을 청크 단위로 무한 확장(Unlimited chunked dataset, `maxshape=(None,)`)하는 것만 허용된다.

**2. 희소 행렬 동적 확장 시의 데이터 불일치 위험 및 대안 포맷의 우위**:
*   **희소 행렬 동시성 파괴 위험**:
    *   CSR 형식은 `data`, `indices`, `indptr`라는 3개의 데이터셋이 서로 완벽하게 동기화되어야만 유효한 행렬로 성립한다.
    *   새로운 세포가 추가될 때 라이터(Writer)는 `data`와 `indices`에 $\text{nnz}$ 원소를 추가하고 `indptr`에 오프셋을 갱신해야 한다.
    *   HDF5는 다중 데이터셋 간의 원자적 트랜잭션(ACID Multi-dataset Transaction)을 지원하지 않는다. 따라서 라이터가 `data`만 플러시하고 `indptr`를 아직 플러시하지 않은 미세한 시간 간격 동안 리더(Reader)가 파일을 읽으면, 포인터 범위를 벗어난 메모리 참조(Out-of-bounds index) 또는 데이터 오염 크래시가 발생한다.
*   **Zarr 및 TileDB의 우위**:
    *   **Zarr**: 각 청크가 독립된 불변 파일(Immutable object)로 기록되므로 락 경합 없이 클라우드 객체 저장소에서도 안전하게 동시 읽기가 가능하다.
    *   **TileDB**: 데이터베이스 수준의 타임스탬프 기반 원자적 커밋(ACID Multi-writer/Multi-reader) 및 다차원 희소 인덱싱을 자체 지원하므로 실시간 스트리밍 및 동시 질의에 완벽히 부합한다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- anndata 0.13.4 (github.com/scverse/anndata). 원 논문: Virshup I, Rybakov S, Theis FJ, et al. anndata: Access and store annotated data matrices. *JOSS*. 2024;9(101):4371. doi: 10.21105/joss.04371
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
<!-- REFERENCES:END -->

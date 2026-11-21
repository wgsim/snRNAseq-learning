# 09. 프로덕션 파이프라인 트러블슈팅, 엔터프라이즈 모범 사례 및 캡스톤 프로젝트 루브릭 (Production Troubleshooting, Best Practices & Capstone Rubric)

---

## 학습 목표 (Cognitive Goals)

1. **프로덕션 유전체학 파이프라인 상위 10대 실패 모드(Top 10 Failure Modes)의 근본 원인 규명**: 공유 파일시스템 락 타임아웃, `publishDir` 심볼릭 링크 단절, 디스크 공간 고갈(ENOSPC), 컨테이너 UID/GID 권한 불일치, 리눅스 OOM 킬러(Exit Code 137)의 커널 및 스토리지 메커니즘을 해부한다.
2. **지수 백오프(Exponential Backoff) 및 동적 재시도 수리 모델을 적용한 엔터프라이즈 에러 복구 체계 구축**: 실패 확률과 네트워크 지연 시간에 따른 최적 대기 시간 $T_{\text{wait}}(k) = \min(T_{\text{max}}, T_{\text{base}} \cdot 2^{k-1} + \mathcal{U}(0, J))$을 정식화하고, `retry`, `ignore`, `terminate` 복구 전략을 프로덕션 Nextflow 코드로 구현한다.
3. **스크래치 작업 디렉토리 자동 정리 및 가비지 컬렉션(GC) 수명주기 통제**: `nextflow clean` 알고리즘과 클라우드 임시 블록 스토리지(EBS/Local NVMe) 수명주기 정책을 연계하여 스토리지 점유를 $O(N \cdot S_{\text{BAM}})$에서 $O(N \cdot S_{\text{h5ad}})$로(샘플 수에 비례하되 상수 계수를 크게 줄여) 억제한다.
4. **전이성 전립선암(mCRPC) 환자 코호트 대상 프로덕션 파이프라인 & MLOps 5단계 캡스톤 프로젝트 루브릭 체계화**: 100점 만점 기준의 상세 평가 지표를 정립하고, 임상 생검 검체의 무결성, 실시간 서빙 지연시간, 클라우드 비용 효율성을 총체적으로 검증한다.

---

## 1. 프로덕션 단일세포 파이프라인 상위 10대 실패 모드 (Top 10 Failure Modes)

엔터프라이즈 유전체학 파이프라인 운영 중 발생하는 치명적 중단 사고는 단순한 코드 버그보다 인프라, 커널, 스토리지, 컨테이너 계층의 상호작용 장애에서 기인한다.

```
       [Production Pipeline Infrastructure Failure Taxonomy]

           Host Operating System & Kernel Layer
  ┌───────────────────────────────────────────────────────────────┐
  │ [Failure 5] Silent Linux OOM Killer (SIGKILL / Exit 137)      │
  │ [Failure 4] Docker Container UID/GID Permission Mismatch      │
  │ [Failure 10] Zombie Java JVM Daemons & .nextflow Cache Locks  │
  └───────────────────────────────┬───────────────────────────────┘
                                  │
           Distributed Storage & File System Layer
  ┌───────────────────────────────┴───────────────────────────────┐
  │ [Failure 1] Lustre / NFS POSIX File Locking Deadlock Timeouts │
  │ [Failure 2] Broken Dangling Symlinks in publishDir Staging    │
  │ [Failure 3] Out-of-Disk Space (ENOSPC) in work/ Scratch Dir   │
  └───────────────────────────────┬───────────────────────────────┘
                                  │
           Cloud Compute & Container Runtime Layer
  ┌───────────────────────────────┴───────────────────────────────┐
  │ [Failure 6] AWS Batch Spot Reclaim (SIGTERM; 143 possible)    │
  │ [Failure 8] HTTP/TCP Remote Reference Streaming Timeouts      │
  │ [Failure 9] Unpinned Docker Floating Tags (:latest) Drift     │
  │ [Failure 7] FASTQ Barcode Corruptions & Truncated gzip Blocks │
  └───────────────────────────────────────────────────────────────┘
```

### 1.1 실패 모드 정밀 분석 및 근본 원인 (Root Cause)

#### 1. 공유 파일시스템 분산 락 타임아웃 (POSIX/NFS File Locking Deadlock)
* **증상**: NFS·Lustre·GPFS·BeeGFS 등 POSIX 파일 잠금을 지원하지 않는 공유 파일시스템에서 실행 디렉터리를 사용하면 시작 시 `Can't open cache DB: .../.nextflow/cache/<uuid>/db` 오류가 발생한다(GitHub nextflow-io/nextflow 이슈 #1479, #6996에 보고된 증상).
* **근본 원인**: Nextflow의 캐시 DB(`.nextflow/cache`, LevelDB 기반)가 파일 잠금을 필요로 하는데 해당 파일시스템이 이를 지원하지 않거나 락 데몬이 동작하지 않음.
* **해결책**: 캐시 DB를 잠금을 지원하는 로컬 디렉터리에 두고 대용량 work만 공유 스토리지에 둔다. 로컬 디렉터리에서 `nextflow run ... -w /shared/work`로 실행하거나, 공식 환경변수 `NXF_CACHE_DIR`(24.02.0-edge부터, 기본 `$launchDir/.nextflow`)로 캐시 위치를 지정한다.
  ```bash
  export NXF_CACHE_DIR=/local/scratch/$USER/nxf_cache
  nextflow run main.nf -w /shared/work
  ```
  (참고: 파일 잠금을 끄는 `nextflow.fileSys.disableFileLock` 같은 설정 옵션은 공식 설정 레퍼런스에 없다.)

#### 2. `publishDir` 복사 실패 및 끊어진 심볼릭 링크 (Dangling Symlinks)
* **증상**: 파이프라인은 성공했다고 표시되나 최종 결과 디렉토리의 BAM 또는 h5ad 파일이 0바이트이거나 링크 대상이 존재하지 않아 다운스트림 분석 스크립트가 크래시됨.
* **근본 원인**: `publishDir mode: 'symlink'`를 사용한 상태에서 스크래치 디렉토리(`work/`)를 조기 삭제하거나, 서로 다른 마운트 포인트 간에 상대 경로가 왜곡됨.
* **해결책**: 프로덕션에서는 반드시 `mode: 'copy'` 또는 안전한 하드링크 `mode: 'link'`를 강제 적용:
  ```groovy
  publishDir path: "${params.outdir}/anndata", mode: 'copy', overwrite: true
  ```

#### 3. 스크래치 작업 디렉토리 디스크 공간 고갈 (ENOSPC - Error No Space Left)
* **증상**: 얼라인먼트 중 `write error: No space left on device` 발생.
* **근본 원인**: STARsolo 정렬 과정에서 수십 기가바이트의 압축 해제 FIFO 파이프 및 정렬되지 않은 중간 임시 BAM 파일이 `/scratch` 볼륨에 누적됨.
* **해결책**: 파이프라인 완료 직후 선행 중간 산출물을 자동으로 정리하거나, 클라우드 워커 노드에 500GB 이상의 고속 임시 NVMe 인스턴스 스토리지를 마운트하고 파이프라인 실행 옵션에 `--save_align_intermeds false` 명시.

#### 4. 컨테이너 권한 불일치 (UID/GID Permission Mismatch)
* **증상**: 도커 컨테이너 내부에서 `Permission denied: /app/work/...` 에러 발생.
* **근본 원인**: 호스트 사용자의 UID(예: 1002)와 컨테이너 내부 기본 실행 사용자(`root` UID 0 또는 `user` UID 1000)가 일치하지 않아 호스트 볼륨에 결과물을 쓰지 못함.
* **해결책**: Docker 런타임에 호스트 사용자 ID를 동적 주입:
  ```groovy
  docker.runOptions = "-u \$(id -u):\$(id -g) --userns=host"
  ```

#### 5. 침묵의 리눅스 OOM 킬러 (Silent Linux OOM Killer - Exit Code 137)
* **증상**: 파이프라인 로그에 명시적인 소프트웨어 스택 트레이스 없이 작업이 즉각 사망하며 `Command error: ... Exit code: 137` 기록.
* **근본 원인**: 프로세스의 물리적 메모리(RSS)가 cgroups에 할당된 `memory.limit_in_bytes`를 초과하여 리눅스 커널 OOM 킬러가 프로세스에 `SIGKILL` (신호 9, $128 + 9 = 137$)을 강제 전송함.
* **해결책**: Nextflow의 `task.attempt`를 활용하여 OOM 감지 시 메모리를 지능적으로 2배 증설하는 동적 재시도 전략 구현.

#### 6. 클라우드 스팟 인스턴스 강제 회수 (SIGTERM / 종료 코드 143일 수 있음)
* **증상**: AWS Batch 또는 Google Cloud 환경에서 작업 실행 중 호스트 VM이 회수되어 작업이 실패함. 컨테이너 프로세스가 SIGTERM($128 + 15 = 143$)을 받으면 `Exit code: 143`이 보일 수 있으나, 스팟 회수 시 종료 코드가 항상 기록된다고 보장되지는 않는다.
* **근본 원인**: 클라우드 공급자가 온디맨드 용량 확보를 위해 스팟 VM을 회수함(AWS는 2분 전 중단 통지).
* **해결책**: 종료 코드로 스팟 회수를 판별하려 하지 말고 Batch 수준 스팟 재시도 옵션(`aws.batch.maxSpotAttempts`, `google.batch.maxSpotAttempts`)을 사용한다. 온디맨드 큐로의 폴백은 `task.attempt` 기반 `queue` 클로저로 직접 구성해야 하는 사용자 설정이다. 143은 타임아웃·수동 종료 등 다른 원인으로도 발생하므로 `errorStrategy`의 일반 재시도 대상 코드의 하나로만 취급한다.

#### 7. FASTQ 파일 헤더 포맷 불일치 및 손상된 gzip 블록
* **증상**: 시퀀싱 리드 수억 개 중 8천만 번째 리드에서 정렬기가 돌연 중단(`gzip: invalid block type`).
* **근본 원인**: 네트워크 전송(FTP/S3 rsync) 중단으로 인한 불완전 파일 전송, 또는 다중 레인(Multi-lane) 병합 과정에서 Read 1과 Read 2의 레코드 순서 동기화 파괴.
* **해결책**: 업스트림 진입 단계에서 `gzip -t` 무결성 검사 및 바코드 길이 정규성 사전 검증 스크립트 강제 실행.

#### 8. 원격 참조 유전체 스트리밍 소켓 타임아웃
* **증상**: 수십 개 태스크가 동시에 Ensembl/NCBI 서버로부터 FASTA 또는 GTF를 다운로드하려다 연결 거부(Connection Refused) 발생.
* **근본 원인**: 외부 웹 서버의 Rate Limiting 및 임상 파이프라인의 외부 인터넷 직접 의존성.
* **해결책**: 참조 유전자 세트는 파이프라인 실행 전 사전에 로컬 공유 스토리지 또는 사내 S3 미러 버킷에 완전히 동기화하고 캐시 경로(`--fasta /shared/ref/GRCh38.fa`)를 정적으로 전달.

#### 9. 고정되지 않은 플로팅 도커 태그 (`latest`) 사용으로 인한 재현성 파괴
* **증상**: 1년 전 완벽히 동작하던 동일한 파이프라인 코드가 새로운 서버에서 실행 시 모듈 문법 에러를 일으킴.
* **근본 원인**: Dockerfile이나 Nextflow 설정에서 `scanpy:latest`를 참조하여, 베이스 이미지 내부의 파이썬 라이브러리 메이저 버전이 업그레이드됨.
* **해결책**: 반드시 OCI 불변 다이제스트 해시를 사용하여 이미지를 고정:
  ```groovy
  container = 'quay.io/biocontainers/scanpy@sha256:<조회한_다이제스트>'   // 실제 값은 `skopeo inspect docker://quay.io/biocontainers/scanpy:<태그>`로 확인 (플레이스홀더)
  ```

#### 10. 좀비 Nextflow JVM 데몬 프로세스의 캐시 락 점유
* **증상**: 터미널이 비정상 종료된 후 재실행 시 `Unable to acquire lock on session with ID ...` 에러 발생(Nextflow 캐시 DB 세션 잠금).
* **근본 원인**: 같은 디렉터리에서 다른 Nextflow 실행이 아직 돌고 있거나, 비정상 종료된 실행이 `.nextflow/cache`의 세션 잠금을 남겨 둔 경우.
* **해결책**: 먼저 해당 실행이 실제로 동작 중인지 확인한다(`ps aux | grep -i nextflow`). 의도치 않은 잔존 프로세스라면 해당 PID만 종료한다(`kill <PID>`, 필요 시에만 `kill -9 <PID>`). `pgrep -f`로 모든 Nextflow를 일괄 `kill -9` 하면 다른 사용자/실행의 작업까지 종료될 수 있으므로 피한다. 프로세스가 없는데 잠금만 남아 있다면 해당 세션의 잠금 파일을 확인 후 정리한다.

---

## 2. 엔터프라이즈 에러 전략 및 지수 백오프(Exponential Backoff) 수리 모델

네트워크 순단, 임시 I/O 병목, 또는 클라우드 스팟 노드 회수와 같은 일시적 결함(Transient Failures)에 대해 즉각적인 재시도를 수행하면, 경합 상태를 더욱 악화시키는 재시도 폭풍(Retry Storm)을 유발할 수 있다.

### 2.1 지수 백오프와 지터(Jitter)의 수학적 정식화

재시도 시도 횟수를 $k \in \{1, 2, \dots, K_{\text{max}}\}$, 기본 대기 시간을 $T_{\text{base}}$, 최대 상한 대기 시간을 $T_{\text{max}}$, 균등 무작위 지터 상한을 $J$라 할 때, $k$번째 재시도 직전의 대기 시간 $T_{\text{wait}}(k)$는 다음과 같이 산출된다:

$$T_{\text{wait}}(k) = \min\left( T_{\text{max}}, \; T_{\text{base}} \times 2^{k-1} \right) + U_k, \quad U_k \sim \operatorname{Uniform}(0, J)$$

* **지터($U_k$)의 당위성**:
  동일한 순간에 네트워크 타임아웃을 겪은 100개의 태스크가 정확히 동일한 지수 간격($2\text{초}, 4\text{초}, 8\text{초}$)으로 재시도하면 서버에 주기적인 충격파(Resonance Congestion)가 가해진다. 무작위 지터를 부가함으로써 재시도 시점을 시간 축 상에 균일하게 분산(Decorrelation)시킨다.

### 2.2 Nextflow 지능형 에러 전략 블록

```groovy
// conf/error_strategies.config

process {
    // 1. 종료 코드별 에러 전략 분기 (미지의 종료 코드는 조용히 무시하지 않고 finish로 처리)
    //    137/139/143 등은 OOM·시그널·인프라 이벤트에 흔한 코드이나 원인을 단정하지 않고 재시도 대상으로만 분류
    errorStrategy = { task.exitStatus in [143, 137, 104, 134, 139, 255] ? 'retry' : 'finish' }
    maxRetries    = 3   // 최초 1회 + 재시도 3회 = 최대 4회 시도

    // 2. 시도 횟수 기반 자원 점진적 확장 (공식 문서의 dynamic directive 패턴)
    //    자원 클로저에서 자기 자신(task.memory)을 참조하면 순환 참조가 되므로 기준값 × task.attempt 사용
    memory = { 8.GB * task.attempt }
}

// 3. 재시도 전 지수 백오프+지터: 별도 'dynamicDelay' 디렉티브는 없으며,
//    공식 문서 "Dynamic Retry with backoff" 패턴처럼 errorStrategy 클로저 안에서 sleep 한다.
process {
    withLabel: 'network_io' {
        errorStrategy = {
            long backoff = Math.min(300_000L, 10_000L * (1L << (task.attempt - 1)))
            sleep(backoff + new Random().nextInt(15_000))
            return 'retry'
        }
        maxRetries = 5
    }
}
```

---

## 3. 스크래치 작업 디렉토리 가비지 컬렉션(GC) 및 정리

대규모 단일세포 전사체 파이프라인에서 500개 샘플을 처리할 경우, `work/` 디렉토리의 용량은 수십 테라바이트에 도달하여 막대한 클라우드 EBS/S3 스토리지 비용을 발생시킨다.

### 3.1 디스크 증가율의 점근적 한계와 정리 수학

샘플당 FASTQ 크기를 $S_{\text{FASTQ}}$, 생성되는 중간 정렬 BAM 크기를 $S_{\text{BAM}}$, 최종 필터링된 h5ad 크기를 $S_{\text{h5ad}}$라 하자 ($S_{\text{h5ad}} \ll S_{\text{BAM}} \ll S_{\text{FASTQ}}$).
- **정리 미수행 시 스토리지 공간 복잡도**:
  $$S_{\text{disk}}(N) = \sum_{i=1}^N \left( S_{\text{FASTQ}}^{(i)} + S_{\text{BAM}}^{(i)} + S_{\text{temp}}^{(i)} + S_{\text{h5ad}}^{(i)} \right) = O(N \cdot S_{\text{BAM}})$$
- **파이프라인 완료 후 가비지 컬렉션 적용 시**:
  $$S_{\text{disk, final}}(N) = \sum_{i=1}^N S_{\text{h5ad}}^{(i)} = O(N \cdot S_{\text{h5ad}})$$
  BAM과 중간 파일이 제거되면 최종 유지 스토리지 용량은 전체의 **$5\% \sim 10\%$ 이하로 급감**한다.

### 3.2 안전한 자동 정리 명령어

Nextflow의 `clean` 명령으로 이전 실행의 work 파일을 삭제한다(선별 삭제가 아니라 해당 실행의 work 전체 대상이므로, 최종 결과가 `publishDir`로 발행된 뒤에 수행한다):

```bash
# 0. 먼저 드라이런(-n)으로 삭제 대상 확인. 실행명을 생략하면 가장 최근 실행이 대상이다
nextflow clean -n

# 1. 가장 최근 실행의 work 파일 삭제 (-k: 임시 파일만 삭제하고 로그·메타데이터는 유지, -f: 강제)
#    주의: work를 지우므로 해당 실행의 -resume 캐시도 사라진다
nextflow clean -f -k

# 2. 특정 실행 ID(Run Name)에 속한 임시 파일 전부 삭제
nextflow clean -f nostalgic_euler

# 3. S3 work 접두사 전체 삭제 (파괴적 명령: --dryrun으로 먼저 확인. 장기 관리는 S3 수명주기 규칙 권장)
aws s3 rm s3://clinical-genomics-bucket/work/ --recursive --exclude "*.filtered.h5ad" --dryrun
```

---

## 4. 단일세포 프로덕션 파이프라인 & MLOps 캡스톤 프로젝트 평가 루브릭 (100점 만점)

본 과정의 최종 역량 평가를 위한 캡스톤 프로젝트 루브릭은 5대 핵심 기술 역량을 종합적으로 검증하도록 설계되었다.

```
       [Production Pipelines & MLOps 100-Point Capstone Rubric]

  ┌──────────────────────────────────────────────────────────────┐
  │ Stage 1: Pipeline Architecture & Scalable Alignment (20 pts) │
  │ - Nextflow DSL2 channel typing & metadata decoupling         │
  │ - STARsolo vs Alevin-fry asymptotic resource optimization    │
  ├──────────────────────────────────────────────────────────────┤
  │ Stage 2: Data Provenance, FAIR & Cryptographic Audit (20 pts)│
  │ - SHA-256 hash-chain chain-of-custody                        │
  │ - Software versions manifest automation (audit support)      │
  ├──────────────────────────────────────────────────────────────┤
  │ Stage 3: Automated Single-Cell MLOps Model Serving (20 pts)  │
  │ - scVI/scANVI latent inference & ONNX operator fusion        │
  │ - FastAPI asynchronous microservice & t-NEPC clinical alert  │
  ├──────────────────────────────────────────────────────────────┤
  │ Stage 4: Cloud Spot Cost Optimization & Fault Tolerance(20)  │
  │ - Spot interruption resilience & task-hash resume check      │
  │ - Dynamic memory escalation (137 OOM recovery)               │
  ├──────────────────────────────────────────────────────────────┤
  │ Stage 5: Enterprise Hardening, CI/CD & Clinical Test (20 pts)│
  │ - GitHub Actions automated lint & minimal test dataset CI    │
  │ - mCRPC clinical cohort end-to-end TAT < 24h SLA compliance  │
  └──────────────────────────────────────────────────────────────┘
```

| 평가 단계 | 세부 평가 항목 | 배점 | 채점 기준 (Grading Criteria) |
| :--- | :--- | :--- | :--- |
| **Stage 1: 파이프라인 아키텍처 및 대규모 정렬 (20점)** | DSL2 채널 설계 및 모듈화 | 10점 | `meta` 맵과 파일 페이로드가 완벽히 결합되어 크로스 컨태미네이션 없이 비동기 라우팅되는가? |
| | 얼라이너 벤치마크 및 리소스 튜닝 | 10점 | STARsolo와 Alevin-fry의 메모리 복잡도($O(G)$ vs $O(K)$)를 분석하고 데이터 특성에 맞게 선택했는가? |
| **Stage 2: 데이터 출처 추적성 및 암호학적 감사 (20점)** | SHA-256 무결성 검증 체인 | 10점 | FASTQ부터 BAM, h5ad까지 단방향 해시 체인 및 무결성 검증 보고서가 자동 생성되는가? |
| | FAIR 원칙 및 규제 요건 지원 리포트 | 10점 | 온톨로지(UBERON, NCIT) 정규화 및 `versions.yml` 소프트웨어 버전 명세가 완벽한가? |
| **Stage 3: 단일세포 MLOps 서빙 및 실시간 추론 (20점)** | 딥러닝 모델 직렬화 및 최적화 | 10점 | PyTorch 인코더를 ONNX Runtime으로 변환하고 연산자 융합을 통해 지연시간을 5배 단축시켰는가? |
| | 실시간 마이크로서비스 및 임상 경고 | 10점 | FastAPI를 통해 $p_{99} < 50\text{ms}$로 세포 유형을 분류하고 t-NEPC 위험 지수를 정확히 산출하는가? |
| **Stage 4: 클라우드 비용 최적화 및 내결함성 (20점)** | 스팟 인스턴스 내결함성 설계 | 10점 | 스팟 회수(`maxSpotAttempts` 재시도) 및 OOM(137) 발생 시 무손실 캐시 재개(`-resume`)가 동작하는가? |
| | 전산 자원 라이트사이징 및 TCO 절감 | 10점 | $Q_{0.99}$ 기반 동적 메모리 스케일링을 통해 온디맨드 대비 60% 이상의 비용 절감을 입증했는가? |
| **Stage 5: 엔터프라이즈 하드닝 및 CI/CD (20점)** | GitOps CI/CD 파이프라인 자동화 | 10점 | GitHub Actions에서 린팅 및 최소 테스트 데이터셋 엔드투엔드 검증이 자동 실행되는가? |
| | 임상 코호트 SLA 검증 및 운영성 | 10점 | 환자 생검 코호트에 대해 24시간 이내의 처리 완료 시간(TAT)과 누락 없는 감사 기록(로그·체크섬)을 확보했는가? |

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 리눅스 커널 OOM 킬러(Exit Code 137)의 cgroups 메모리 회수 메커니즘과 Nextflow 동적 메모리 확장 수리 모델
전이성 전립선암 환자의 대용량 골수 생검 검체(세포 수 35,000개, 총 3억 리드)를 STARsolo로 정렬하던 중, 작업 시작 45분 만에 아무런 경고 로그 없이 프로세스가 급작스럽게 종료되며 Nextflow 콘솔에 `Exit code: 137`이 기록되었다.
1. 호스트 머신에 128 GB의 물리 RAM이 남아있었음에도 불구하고 도커 컨테이너 내부의 프로세스가 왜 리눅스 커널에 의해 강제 사살(`SIGKILL`)되었는지 cgroups `memory.limit_in_bytes`와 OOM 스코어 수리식을 통해 설명하시오.
2. 이 문제를 해결하기 위해 엔지니어가 `nextflow.config`에 `errorStrategy = 'retry'`와 함께 시도 횟수($k$)에 따라 메모리를 $M(k) = M_0 \times 2^{k-1}$로 동적 증설하는 로직을 작성하였다. 초기 할당 메모리가 $M_0 = 32\text{ GB}$이고 태스크의 실제 최대 물리 피크 메모리가 $48\text{ GB}$일 때, 몇 번째 시도에서 작업이 성공 완주되는지 계산하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. cgroups OOM 킬러 작동 메커니즘 및 수리식**:
* 호스트 머신의 전체 여유 RAM과 무관하게, 도커 컨테이너는 리눅스 커널의 **cgroups(Control Groups)** 제어 서브시스템에 의해 프로세스 그룹별 메모리 상한선($M_{\text{limit}} = \text{memory.limit\_in\_bytes}$, 예: 32 GB)이 강제된다.
* 프로세스의 물리 점유 메모리(RSS)와 익명 페이지(Anonymous Pages)의 합이 $M_{\text{limit}}$에 도달하면 커널은 먼저 페이지 캐시 회수(Page Reclamation)를 시도한다.
* 회수에 실패하면 커널은 컨테이너 네임스페이스 내부의 프로세스들에 대해 OOM 점수(Out-Of-Memory Score)를 계산한다:
  $$\text{oom\_score} \approx \frac{\text{RSS} + \text{swap}}{M_{\text{limit}}} \times 1000 + \text{oom\_score\_adj}$$
* STARsolo 프로세스가 최대 점수(1000점)를 기록하게 되므로, 커널의 `mem_cgroup_out_of_memory()` 함수가 트리거되어 프로세스에 방어 불가능한 `SIGKILL` (신호 9)을 전송한다. 프로세스는 즉각 사살되며 쉘 종료 코드는 $128 + 9 = 137$을 반환하게 된다.

**2. 동적 메모리 증설에 따른 완주 시도 계산**:
* 초기 할당량 $M_0 = 32\text{ GB}$, 실제 필요 피크 메모리 $M_{\text{peak}} = 48\text{ GB}$.
  - **시도 1 ($k=1$)**: $M(1) = 32 \times 2^0 = 32\text{ GB} < 48\text{ GB}$ $\implies$ 메모리 부족으로 OOM 발생 (Exit Code 137 실패).
  - **시도 2 ($k=2$)**: Nextflow가 실패를 감지하고 메모리를 지능적으로 2배 증설:
    $$M(2) = 32 \times 2^{2-1} = 32 \times 2^1 = \mathbf{64\text{ GB}}$$
    $64\text{ GB} > 48\text{ GB}$를 만족하므로, cgroups 한계선 이내에서 안전하게 연산이 진행되어 **2번째 시도(Attempt 2)**에서 성공적으로 완주된다.

</details>

---

### 문제 2: 임상 시험 유전체 파이프라인에서 `publishDir mode: 'symlink'` 사용이 유발하는 데이터 영구 소실 참사 및 규제 위반 분석
한 임상 유전체학 연구소에서 비용 절감을 위해 S3 스토리지 사용량을 줄이고자 `publishDir` 지시자에 `mode: 'symlink'`를 적용하여 100명의 전립선암 환자 분석을 완료하였다. 이후 정기 유지보수 스크립트가 실행되어 1주일이 경과한 Nextflow 작업 임시 디렉토리(`work/`)를 `rm -rf work/`로 영구 삭제하였다. 다음 날 병리과 의사가 환자 치료 방침 결정을 위해 결과 폴더(`results/anndata/`)의 h5ad 파일을 열려고 했으나 모든 파일이 `broken symbolic link` 에러를 내며 열리지 않았다.
1. `symlink` 모드가 파일의 실제 데이터 블록(Data Inode/Blocks)을 어떻게 참조하기에 이와 같은 데이터 영구 유실 참사가 발생했는지 POSIX 파일시스템 구조 관점에서 설명하시오.
2. 임상 진단 및 동반진단(CDx) 파이프라인에서 데이터의 무결성과 독립성을 보장하기 위해 반드시 채택해야 하는 `publishDir` 설정과 스토리지 아키텍처의 모범 사례를 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 심볼릭 링크(Symlink) 참조 구조와 데이터 유실 원인**:
* POSIX 파일시스템에서 심볼릭 링크는 실제 파일 데이터의 아이노드(Inode)나 디스크 블록을 소유하지 않고, 단지 원본 파일의 경로 문자열(Path String, 예: `../../work/ab/c1d2e3/sample.h5ad`)만을 담고 있는 특수 포인터 파일이다.
* 파이프라인 실행 중 실제 연산 결과물(AnnData h5ad 파일의 이진 바이트 블록)은 `work/` 디렉토리 내부의 고유 해시 서브폴더에 저장된다.
* `rm -rf work/`가 실행되는 순간 원본 파일의 아이노드와 디스크 데이터 블록이 파일시스템 상에서 영구 해제(Deallocated)된다.
* 결과 디렉토리(`results/anndata/`)에 남아있던 심볼릭 링크들은 가리키는 대상(Target Path)이 소멸된 **단절 링크(Dangling/Broken Symlink)**로 전락하며, 원본 데이터가 존재하지 않으므로 어떠한 방법으로도 데이터를 복구할 수 없는 영구 소실이 발생한다. 환자의 치료 결정에 사용될 유전체 데이터가 소멸되어 심각한 의료 사고 및 규제 위반에 직면하게 된다.

**2. 임상 규제 준수를 위한 모범 사례**:
1. **`mode: 'copy'` 강제 적용**:
   임상 결과 디렉토리에는 포인터가 아닌 실제 바이트 데이터를 물리적으로 복사해야 한다:
   ```groovy
   publishDir path: "${params.outdir}/clinical_reports", mode: 'copy', overwrite: true
   ```
2. **단일 파일시스템 내 하드링크(`mode: 'link'`) 활용**:
   로컬 고성능 파일시스템 내에서 스토리지를 절약하면서도 안전성을 확보하려면 하드링크를 사용할 수 있다. 하드링크는 동일한 Inode를 공유하며 참조 카운트(Link Count)를 1 증가시키므로, `work/` 디렉토리가 삭제되더라도 `results/`의 링크 카운트가 1로 유지되어 실제 디스크 데이터 블록이 절대 소멸되지 않는다.
3. **클라우드 객체 스토리지 이중화 영속화**:
   임상 최종 산출물은 로컬 디스크가 아닌 AWS S3의 버전 관리(Versioning) 및 객체 잠금(Object Lock/WORM: Write Once Read Many)이 활성화된 규제 준수 버킷으로 즉시 복제하여 원천적인 삭제를 방지해야 한다.

</details>

---

### 문제 3: 지수 백오프와 지터(Jitter) 알고리즘이 클라우드 API 경합(Throttling) 및 재시도 폭풍(Retry Storm)을 방지하는 수리적 기전
AWS Batch 클라우드 환경에서 1,000개의 단일세포 전처리 작업이 동시에 완료되는 순간, 모든 작업이 중앙 S3 버킷의 같은 접두사에 동시에 접근하여 결과를 업로드하려다 S3 접두사당 요청 한도(3,500 PUT/COPY/POST/DELETE requests/sec)를 초과하여 전원 `HTTP 503 Slow Down (Throttling)` 에러를 수신하였다. (가상 시나리오: 각 태스크가 같은 접두사 아래 약 10개 객체를 동시에 PUT하여 순간 약 10,000건/초가 발생한다고 가정. 단순히 1,000건/초는 3,500 PUT/s 한도 미만이다.)
1. 만약 1,000개의 태스크가 고정된 대기 시간(Fixed Delay, 예: 정확히 5초 대기 후 전원 재시도)을 취할 경우, 클라우드 시스템 전체에 어떤 주기적 공진 붕괴 현상(Resonance Congestion/Retry Storm)이 발생하는지 설명하시오.
2. 완전 무작위 지터를 포함한 완전 지수 백오프(Full Jitter Exponential Backoff) 공식 $T_{\text{wait}}(k) = \operatorname{Uniform}(0, \min(T_{\text{max}}, T_{\text{base}} \cdot 2^{k-1}))$을 적용했을 때, 1,000개 태스크의 재시도 트래픽이 시간 축 상에서 어떻게 분산되어 API 스로틀링을 즉각 해소하는지 확률 밀도 함수(PDF) 관점에서 증명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 고정 지연 시간의 주기적 공진 붕괴(Retry Storm) 메커니즘**:
* $t=0$ 시점에 1,000개의 태스크(각 약 10개 PUT, 합 약 10,000건)가 S3 PUT 요청을 보내 Throttling(503)을 받았다.
* 모든 태스크가 고정 지연 $\Delta t = 5.0\text{초}$를 취한다면, $t=5.0\text{초}$가 되는 순간 1,000개의 태스크가 **위상 동기화(Phase Synchronization)**되어 정확히 동일한 밀리초에 다시 수천 건의 동시 요청을 S3 엔드포인트에 쏟아붓는다.
* 이 충격파는 필연적으로 2차 Throttling을 유발하고, 다시 $t=10.0\text{초}$, $t=15.0\text{초}$에 동일한 규모의 스파이크가 주기적으로 반복된다. 결과적으로 API 서버는 회복할 여유를 전혀 얻지 못하고 전체 파이프라인이 집단 타임아웃으로 동반 전멸하는 재시도 폭풍(Retry Storm) 재앙에 빠진다.

**2. Full Jitter 지수 백오프의 확률 밀도 분산 증명**:
* 1회차 실패 후 $k=1$일 때, 상한선은 $M_1 = T_{\text{base}} \times 2^0 = 10\text{초}$이다. 각 태스크 $i \in \{1, \dots, 1000\}$의 대기 시간 $T_{\text{wait}}^{(i)}$는 독립 균등분포를 따른다:
  $$T_{\text{wait}}^{(i)} \sim \operatorname{Uniform}(0, 10)$$
* 단위 시간당 도착 태스크 수(재시도 트래픽 밀도)의 기댓값은 균등 확률 밀도 함수 $f(t) = \frac{1}{10}$ ($0 \le t \le 10$)에 의해 다음과 같이 일정하게 분산된다:
  $$\mathbb{E}[\text{Requests per second}] = \frac{10{,}000\text{ requests}}{10\text{ seconds}} = 1{,}000\text{ requests/sec}$$
* 순간 약 10,000건의 임펄스 스파이크가 10초 구간에 걸쳐 초당 약 1,000건의 완만한 트래픽으로 평탄화된다. 이는 접두사당 한도(3,500 req/sec) 아래이므로 재시도 분산만으로 스로틀링이 크게 완화될 것으로 기대된다(S3는 요청률에 따라 내부적으로 확장되며 즉시 100% 해소를 보장하지는 않는다).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Nextflow 24.04.4 (nextflow.io). 원 논문: Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35(4):316-319. PMID: 28398311. doi: 10.1038/nbt.3820
- Docker (github.com/moby/moby).
- ONNX Runtime 1.31.0 (github.com/microsoft/onnxruntime).
- FastAPI 0.143.3 (github.com/fastapi/fastapi).
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- scvi-tools 1.5.1 (github.com/scverse/scvi-tools). 원 논문: Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol*. 2022;40(2):163-166. PMID: 35132262. doi: 10.1038/s41587-021-01206-w
<!-- REFERENCES:END -->

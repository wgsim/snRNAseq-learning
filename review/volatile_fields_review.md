# scripts/pipeline_volatile_fields.json 검토

범위는 `scripts/pipeline_volatile_fields.json`의 선언 전부와, Track 07–10 결과 폴더(`tracks/0[7-9]_*/pipelines/results_pipeline_0*`, `tracks/10_*/pipelines/results_pipeline_0*`)에서 선언되지 않은 시간 의존 값입니다. 검토는 코드 읽기에 더해, 실제 트리가 아닌 scratchpad 사본에서 파이프라인을 재실행해 확인했습니다.

- 재실행 환경: `scRNAseq_anal` conda 환경(Python 3.8.16). 파이프라인 7개를 실행했습니다(07-01, 07-02, 08-01, 08-02, 09-01, 10-01, 10-02). 09-02는 시간 값이 출력(print)뿐이라 실행하지 않았습니다.
- `scripts/course_markdown.py`의 `_compare_one`으로 비교한 결과, 재실행한 15개 파일이 모두 "선언된 필드만 마스킹하면 통과"였습니다. 같은 환경과 같은 머신 한 번 비교이므로, 아래 flaky-risk 항목은 재실행 비교를 통과해도 사라지지 않습니다.
- 부하 실험: 바쁜 루프 프로세스 10~24개를 띄운 상태에서 10-01을 3회, 10-02를 1회 재실행했습니다.
- 줄 번호는 모두 `tracks/<트랙>/pipelines/<스크립트>` 기준입니다.
- 판정: **OK** = 선언한 값이 시간/환경에서만 나옴, **over-masked** = 시간이 아닌 값까지 가림, **under-masked** = 시간 의존 값이 선언에서 빠짐(다른 머신에서 flaky), **flaky-risk** = 시간에서 파생된 판정 값(불리언 등)이 선언되지 않음.

## 1. 선언 항목 검토

| 파일 | 필드 | 생성 코드 | 판정 | 비고 |
|---|---|---|---|---|
| 07 `results_pipeline_01/pipeline_01_metrics.json` | `execution_time_seconds` | `07/pipelines/01_pseudobulk…py:701, 741` (`round(time.time()-start_time,2)`) | OK | 시간 측정만 담습니다. 같은 파일의 `fdr`·`log2fc`·오류율 등은 시드 고정 계산이라 비교 대상으로 남아 있고, 재실행에서도 동일했습니다. |
| 07 `results_pipeline_02/pipeline_02_metrics.json` | `execution_time_seconds` | `07/pipelines/02_milo_and_coda…py:566, 601` | OK | 위와 같습니다. |
| 08 `results_pipeline_01/pipeline_01_metrics.json` | `benchmarks.row_slice_csr_sec`, `row_slice_csc_sec`, `col_slice_csc_sec`, `col_slice_csr_sec` | `08/pipelines/01_sparse…py:153–187` (`time.time()` 차이) | OK | 측정 시간입니다. |
| 〃 | `benchmarks.row_slice_speedup_csr_vs_csc`, `col_slice_speedup_csc_vs_csr` | `01…py` 위 시간들의 비 | OK | 시간 비율입니다. |
| 〃 | `benchmarks.spmv_runtime_sec`, `spmv_gflops` | `01…py:194–201` (spmv 시간, FLOP/시간) | OK | `spmv_flops`와 `dense_flops`는 nnz 기반이라 선언하지 않았고 재실행에서 같았습니다. |
| 〃 | `streaming.streaming_time_sec` | `01…py:287, 368` | OK | |
| 〃 | `incremental_svd.runtime_sec` | `01…py:400, 495` | OK | `mean_cosine_similarity`, `concordance_pass`(`>0.95`)는 시간과 무관하며 비교 대상입니다. |
| 〃 | `execution_time_total_sec` | `01…py:577, 728` | OK | |
| 08 `…/pipeline_01_summary.txt` | 줄 규칙 `Row-slicing|Col-slicing|SpMV Throughput|Streaming Time|Total Pipeline Runtime`, 마스크 `\d[\d,]*\.\d+` | `01…py:~745–757` (f-string 서식) | OK | 해당 5줄에 있는 소수는 전부 시간이나 시간 파생값이며, 비율(`673.6x`)도 시간 비입니다. 비슷한 줄 `Peak RAM: 4.97 MB (… PASSED)`는 규칙에 걸리지 않아 비교됩니다. 이 값은 결정적이며 재실행에서 동일했습니다. |
| 08 `…/sparse_memory_layout_benchmark.csv` | `Value` 열, `Metric`이 `Row_Slice_*`·`Col_Slice_*`·`SpMV_*`인 8행 | `01…py:~660–680` | OK | 마스크 대상이 이름으로 정확히 한정됩니다. 메모리 크기, `Compression_Ratio`, `FLOP_Reduction` 행은 비교 대상으로 남습니다. |
| 08 `results_pipeline_02/pipeline_02_metrics.json` | `metadata_slicing.row_time_per_query_ms`, `columnar_time_per_query_ms`, `query_speedup` | `08/pipelines/02_tiledb…py:148–168` | OK | |
| 〃 | `ann_search.brute_time_sec`, `ann_time_sec`, `ann_speedup` | `02…py:230–260` | OK | `mean_recall`(1.0), `recall_target_pass`는 시드 고정이고 선언되지 않아 올바릅니다. |
| 〃 | `hardware_profiling.peak_ram_mb` | `02…py:353, 402` (`tracemalloc` 피크) | OK (주의) | 시간이 아닌 메모리 측정이지만 허용 범주(peak RAM)입니다. 재실행에서 `39.51887`이 `39.51923`으로 실제로 달라졌으므로 선언은 필요합니다. 같은 값이 `summary.txt`에도 있는데 거기서는 마스킹되지 않았습니다(아래 2번 표). |
| 〃 | `execution_time_total_sec` | `02…py:348, 477` | OK | |
| 08 `…/pipeline_02_summary.txt` | 줄 규칙 `Multi-attribute query latency|Query retrieval speedup|ANN Speedup|Total Pipeline Execution Time`, 마스크 소수 | `02…py:~480–505` | OK (일부 under-masked) | 마스킹되는 소수는 모두 시간이나 시간 비입니다. 그런데 `Multi-attribute query latency` 줄에는 `Sub-millisecond pass: True`가 같이 있고, 이 불리언은 비교됩니다(아래 2번 표). |
| 08 `…/ann_vs_bruteforce_knn_benchmark.csv` | `Value`, 행 `Brute_Force_Time`, `ANN_IVF_Time`, `ANN_Speedup` | `02…py:~440–460` | OK | `Mean_Recall`, `Queries_Evaluated`는 비교 대상으로 남습니다. |
| 08 `…/columnar_vs_row_metadata_benchmark.csv` | `Value`, 행 `Row_Query_Latency`, `Columnar_Query_Latency`, `Slicing_Speedup` | `02…py:~420–440` | OK | 메모리 행(`Columnar_Memory`, `Row_Memory`, `Memory_Reduction`)은 비교됩니다. 이 값은 `buffer_info`와 `n_cells*320`으로 계산되어 결정적입니다(`02…py:118–130`). |
| 09 `results_pipeline_01/demux_and_error_correction_summary.json` | `umi_deduplication_benchmarks.execution_time_seconds` | `09/pipelines/01_fastq…py:455, 496, 509` | OK | |
| 09 `…/pipeline_01_summary_report.md` | 줄 규칙 `Deduplication Execution Time`, 마스크 소수 | `01…py:591` | OK | 이 줄에는 시간 값 하나뿐입니다. |
| 10 `results_pipeline_01/timeline_metrics.json` | `cold_execution.duration_seconds`, `warm_execution.duration_seconds` | `10/pipelines/01_dsl2…py:506–548` | OK | `cached_tasks`, `cache_hit_percentage`는 해시 조회로 결정되며 비교됩니다. |
| 10 `…/trace_table.csv` | 열 `duration_s`, `timestamp` (모든 행) | `01…py:203–208, 233–247, 271–293` | OK | `row_key`가 없어 모든 행에 적용되는데, 의도에 맞습니다. `task_id`, `attempt`, `status`, `hash_key`, `work_dir`는 비교됩니다. |
| 10 `results_pipeline_02/model_audit_manifest.json` | `timestamp_utc` | `10/pipelines/02_mlops…py:337` | OK | |
| 〃 | `runtime_environment.python_version` | `02…py:359` | OK | |
| 10 `…/serving_benchmark.csv` | 열 `mean_latency_ms`, `std_latency_ms`, `p50/p95/p99_latency_ms`, `throughput_cells_sec` | `02…py:469–498` (`time.perf_counter`) | OK | 선언된 열 자체는 모두 순수한 측정값입니다. |
| 10 `…/summary_report.md` | 줄 규칙 `Timestamp \(UTC\)` / ISO 시각 | `02…py:~560` | OK | |
| 〃 | 줄 규칙 `^\| \*\*\d+ cells\*\* \|` / 소수 | `02…py:582` | OK (일부 under-masked) | 해당 4행의 소수는 모두 지연·처리량입니다. 같은 행의 `✅ PASS`/`❌ FAIL`은 비교됩니다(아래 2번 표). |
| 〃 | 줄 규칙 `Low-Latency Guarantees` / `\d+\.\d+ ms` | `02…py:586` | OK | `p95 … ms`만 가립니다. |
| 〃 | 줄 규칙 `High Throughput` / `[\d,]+ cells/second|~[\d.]+ seconds` | `02…py:587` | **over-masked** | `cells/second` 숫자(측정값)는 정당한 마스크입니다. 그러나 `~3.8 seconds`는 코드에 고정된 문자열 상수입니다. 측정값이 아니므로 마스킹할 필요가 없고, 오히려 상수의 오류(10,000셀 ÷ 8,432 셀/초 ≈ 1.2초이지 3.8초가 아님)를 가립니다. 제안: 마스크를 `[\d,]+ cells/second`로 줄이고, 3.8초 상수는 코드에서 측정값으로 계산하거나 수정합니다. |

## 2. 선언되지 않은 시간 의존 값 (역방향 검토)

| 파일 | 필드 | 생성 코드 | 판정 | 비고 |
|---|---|---|---|---|
| 10 `…/serving_benchmark.csv` | `sla_compliant` (`p95 < 50 ms`) | `02…py:488` | **flaky-risk (실측 확인)** | p95 지연에서 파생된 불리언이며 선언되지 않았습니다. 이 머신에서 batch 256의 p95는 28.7~36.1 ms로 한계의 59~72%에 있습니다. 바쁜 루프 24개를 띄운 부하 재실행에서는 batch 64(p95 85.5 ms)와 batch 256(p95 129.1 ms)이 `False`로 뒤집혔습니다. 더 느린 머신이나 CI에서도 같은 일이 일어납니다. |
| 10 `…/summary_report.md` | 표 행 `✅ PASS` / `❌ FAIL` | `02…py:581` | **flaky-risk** | 위 `sla_compliant`에서 그대로 파생되며, 줄 규칙은 소수만 마스킹합니다. 같은 행 안에서 숫자는 가려지는데 판정 문자는 비교되는 불일치입니다. |
| 10 `…/summary_report.md` | "strictly compliant with the 50 ms clinical SLA" | `02…py:586` | 참고 | 문구가 SLA 결과와 무관하게 고정되어 있어, 느린 머신에서는 본문이 거짓이 됩니다. 비교는 이를 잡지 못하며, 생성기를 수정할 사안입니다. |
| 10 `…/model_audit_manifest.json` | `runtime_environment.os_platform`, `processor` | `02…py:360–361` (`platform.platform()`, `platform.processor()`) | **under-masked (환경)** | `python_version`은 선언되었지만 같은 부류인 OS 문자열(`macOS-27.0.1-arm64-arm-64bit`)과 프로세서는 선언되지 않았습니다. 다른 OS(Linux CI 등)에서는 바로 실패합니다. 제안: `json_keys`에 두 키를 추가합니다. |
| 08 `…/pipeline_02_metrics.json` | `metadata_slicing.sub_millisecond_pass` (`columnar_time_per_query_ms < 1.0`) | `08/…/02_tiledb…py:178` | **flaky-risk** | 지연에서 파생된 불리언입니다. 현재 0.27~0.28 ms로 3.6배 여유가 있어 느린 머신에서 뒤집힐 수 있습니다. |
| 08 `…/pipeline_02_summary.txt` | `Sub-millisecond pass: True` (같은 줄) | `02…py:~491` | **flaky-risk** | 줄 규칙은 소수만 마스킹하므로 `True`가 비교됩니다. |
| 08 `…/pipeline_02_metrics.json` | `ann_search.speedup_target_pass` | `02…py:268` | **flaky-risk** | `10.0 <= speedup <= 50.0 or speedup > 10.0`은 사실상 `speedup > 10`이며, 시간 비율에서 파생된 불리언입니다. 실측 비율은 21–27배(커밋본 25.07, 재실행 27.30)로 여유가 2배 남짓입니다. 상한 50이 무의미하다는 점도 코드 오류입니다. |
| 08 `…/pipeline_02_summary.txt` | `Hardware Peak RAM: 39.52 MB` | `02…py:504` | **under-masked** | JSON의 `hardware_profiling.peak_ram_mb`는 선언했는데 같은 값을 쓰는 텍스트 줄은 선언하지 않았습니다. 재실행에서 값이 실제로 달라졌고(`39.5189` 대 `39.5192`), 지금은 반올림 결과가 우연히 같지만 경계값 근처에서는 `39.51`/`39.53`으로 달라질 수 있습니다. 제안: 줄 규칙에 `Hardware Peak RAM`을 추가합니다. |
| 08 `…/pipeline_01_metrics.json` | `streaming.peak_ram_mb`, `strict_budget_pass`, CSV `Peak_RAM_MB` | `08/…/01_sparse…py:286–379, 328` | OK (미선언이 맞음) | 같은 `tracemalloc` 부류지만 재실행에서 소수점 15자리까지 같았습니다(`4.973285675048828`). 08-02의 `peak_ram_mb`와 처리 방침이 다르므로, 정책을 한쪽으로 통일하면 좋습니다. 이 값을 선언하지 않는 쪽은 안전합니다. |
| 08 `…/pipeline_01_summary.txt`, `pipeline_02_summary.txt` | `PASSED` 문구 | `01…py:751, 755`; `02…py:497–498, 504` | 참고 | 리터럴 `PASSED`이므로 결과와 무관하게 항상 같습니다. 가독성 문제일 뿐 flaky는 아닙니다. |
| 08 `…/pipeline_02_metrics.json` `gpu_acceleration.*`, `gpu_acceleration_profiling_metrics.csv` | 단계별 CPU/GPU 시간 | `02…py:295–340` | OK (미선언이 맞음) | 상수와 `matrix_size_mb/(64·1024)`로 계산되는 시뮬레이션 값이며 측정값이 아닙니다. 재실행에서 바이트 동일했습니다. |
| 10 `results_pipeline_01/trace_table.csv` | `task_id` 순서 (`asyncio.gather` + `sleep`) | `10/…/01_dsl2…py:176–225, 364` | flaky-risk (낮음) | `task_id`는 task 시작 시점의 카운터로 부여되므로 타이머 순서에 의존합니다. 지연 상수는 0.03+0.01·attempt초이고 asyncio는 마감 시각 순으로 깨웁니다. 바쁜 루프 10개 부하에서 3회 재실행해도 마스킹 밖의 차이가 없었습니다(비교 열 전부 동일). 지금은 안전합니다. |
| 10 `results_pipeline_01/summary_report.md` | "Total Compute Walltime ~1.28s / ~0.02s / 64x speedup" | `01…py:471` | 참고 | 하드코딩된 문자열이라 flaky하지 않습니다(재실행에서 동일). 다만 같은 폴더의 `timeline_metrics.json`(cold 0.394 s, warm 0.0 s)과 서로 맞지 않는 문서 불일치입니다. |
| 07, 09-02 | — | `07/…/01:615`, `09/…/02:594–607` | OK | 경과 시간은 `print`로만 쓰이고 결과 파일에는 `execution_time_seconds`만 들어갑니다. 결과 파일에 다른 시간 값은 없습니다. |
| 10 `results_pipeline_02/model_weights.json`, `model_audit_manifest.json` | `weights_sha256`, `final loss`, `accuracy` | `02…py:~337–360, 629–646` | 참고 | 시간 의존은 아닙니다. 다만 `math.exp`/`math.log` 결과는 libm 구현(macOS와 glibc)에 따라 마지막 자리가 달라질 수 있어, 다른 OS에서는 해시와 손실이 달라질 수 있습니다. 위 `os_platform`과 같은 부류의 환경 의존입니다. |

## 3. 비교 코드 (`scripts/course_markdown.py:1635–1760`)

- 마스킹은 필드 범위입니다. 선언이 없는 파일은 바이트 일치를 요구하고(`rules is None`), 파일 전체를 건너뛰는 경로는 없습니다.
- JSON: 경로가 정확히 일치하는 값만 `<VOLATILE>`로 바꿉니다(`*` 허용). 비교 후에도 키 누락과 타입 차이는 `_first_json_diff`가 잡습니다. 컨테이너를 마스킹할 수 있으나, 선언은 모두 스칼라입니다.
- CSV: 헤더 동일성과 행 수를 확인하고, 선언 열이면서 `row_values`에 속한 행만 건너뜁니다. 그 밖의 열과 행은 모두 비교합니다. `Metric` 열 자체는 마스킹 대상이 아니라서, 행 이름이 달라지면 걸립니다.
- 텍스트 줄 규칙: 줄 수가 같아야 하고, 줄 규칙에 걸린 줄에서 마스크 정규식에 맞는 부분만 가립니다. 같은 줄의 나머지 텍스트(예: `True`, `PASS`)는 비교됩니다. 그래서 위 2번 표의 `Sub-millisecond pass`와 `✅ PASS`가 flaky-risk로 남습니다.
- 사소한 약점 1: CSV 비교는 `zip(x, y)`라서 두 행의 열 개수가 다르면 짧은 쪽에서 잘립니다(헤더와 행 수는 확인함).
- 사소한 약점 2: 선언한 키나 열이 실제 파일에 없어도 오류가 나지 않습니다(오타가 있으면 아무것도 마스킹하지 않아 더 엄격한 방향으로만 틀립니다).
- 사소한 약점 3: JSON은 파싱 후 비교하므로 공백과 키 순서 차이는 무시됩니다.

## 요약

- 선언 항목 중 약 30개(JSON 키, CSV 열/행, 텍스트 줄 규칙)를 코드 줄까지 확인했고, **1건만 over-masked**입니다(10-02 `summary_report.md`의 `~[\d.]+ seconds`가 코드에 고정된 상수 `~3.8 seconds`를 가림). 나머지 선언은 모두 시간/환경에서만 나온 값입니다.
- 과소 선언은 7건입니다.
  - **flaky-risk**: 10-02 `sla_compliant`(CSV)와 `✅ PASS`/`❌ FAIL`(md) 2건(**실측으로 확인**: 부하 상태에서 batch 64·256이 FAIL로 뒤집힘), 08-02 `sub_millisecond_pass`(JSON과 summary.txt) 2건, 08-02 `speedup_target_pass` 1건.
  - **under-masked**: 08-02 `summary.txt`의 `Hardware Peak RAM` 줄, 10-02 manifest의 `os_platform`/`processor` 2건.
- 시드 고정 과학 결과(recall, 정확도, FDR, 오류율, 해시)를 가리는 선언은 없습니다. 시간에서 파생된 판정 값(SLA, 서브밀리초, 속도 목표)이 선언되지 않은 것이 가장 큰 구멍입니다.
- 제안 선언(수정은 하지 않았습니다):
  - `serving_benchmark.csv`의 `csv_columns`에 `sla_compliant`를 추가합니다.
  - `summary_report.md` 둘째 규칙의 마스크를 `\d[\d,]*\.\d+|✅ PASS|❌ FAIL`로 확장합니다.
  - `model_audit_manifest.json`의 `json_keys`에 `runtime_environment.os_platform`, `runtime_environment.processor`를 추가합니다.
  - `08 pipeline_02_metrics.json`의 `json_keys`에 `metadata_slicing.sub_millisecond_pass`, `ann_search.speedup_target_pass`를 추가합니다.
  - `08 pipeline_02_summary.txt`에 `Hardware Peak RAM` 줄 규칙과 `Sub-millisecond pass: (True|False)` 마스크를 추가합니다.
  - `10-02 summary_report.md`의 `High Throughput` 규칙에서 `~[\d.]+ seconds`를 제거합니다.
- 코드 쪽으로 같이 정리할 것: `speedup_target_pass`의 `10 <= s <= 50 or s > 10`(실질 `s > 10`), 하드코딩된 "strictly compliant", `~3.8 seconds`, 10-01 `summary_report.md`의 고정 walltime 문구(실측 0.394 s와 불일치).
- 비교 코드는 필드 단위로 마스킹하며 파일 전체를 건너뛰지 않습니다. 사소한 약점은 CSV 열 개수 불일치를 `zip`이 가리는 것뿐입니다.
- 재실행 환경과 실제 환경이 같다면 지금 선언은 현 머신에서 통과합니다. 위 flaky-risk 항목은 부하나 다른 머신에서 깨집니다.

# 08. 클라우드 비용 경제학, 스팟 인스턴스 내결함성 및 대규모 자원 최적화 (Cloud Cost Economics, Spot Fault-Tolerance & Resource Optimization)

---

## 학습 목표 (Cognitive Goals)

1. **클라우드 컴퓨팅 비용 경제학 및 스팟 인스턴스(Spot Instance) 수학적 모델링**: 온디맨드 대비 크게 할인되는(AWS 공개 상한 최대 90%; 본 장 예시는 약 75% 가정) 스팟 마켓의 가격 결정 메커니즘과 중단 확률(Interruption Probability)을 포아송 점 과정(Poisson Point Process)으로 정식화한다.
2. **Nextflow 캐시 무결성과 작업 디렉토리 해시 알고리즘 해부**: `-resume` 지시자가 선행 작업의 성공 여부를 검증하는 수리적 조건식($H_{\text{task}} = \text{Hash}(S_{\text{script}}, \sum H(I), C_{\text{container}}, E_{\text{env}})$, 개념식)을 완벽히 이해하고, 스팟 중단 발생 시 제로 중복 재실행 메커니즘을 규명한다.
3. **자원 과다 프로비저닝(Over-provisioning) 진단 및 동적 적응형 스케일링**: Nextflow 추적 지표로부터 경험적 분위수($Q_{0.99}$) 기반 메모리/CPU 라이트사이징(Right-sizing) 알고리즘을 도출하고, 지수적 백오프 재시도 전략(`errorStrategy 'retry'`)을 구현한다.
4. **500명 대규모 전립선암(mCRPC) 환자 코호트 비용 최적화 아키텍처 구축**: 500억 개 리드 쌍의 대규모 유전체 데이터를 24시간 이내의 임상 소요 시간(TAT)으로 처리하면서 샘플당 단가를 낮추는(예시 목표 \$75→\$15는 저장·전송·재처리 등을 포함한 총비용 가정이며, Q1의 컴퓨팅 비용만 계산한 샘플당 \$1.5→\$0.48과는 전제가 다름) AWS Batch/하이브리드 HPC 인프라를 배포한다.

---

## 1. 클라우드 비용 경제학과 스팟 인스턴스 동역학

대규모 단일세포 전사체 연구에서 컴퓨팅 인프라 비용은 시퀀싱 시약 비용에 필적하는 주요 예산 항목이다. 상용 퍼블릭 클라우드(AWS, Google Cloud, Azure)는 미사용 여유 컴퓨팅 용량을 시장가로 제공하며(과거의 입찰/경매 방식은 더 이상 쓰이지 않음), 이를 **스팟 인스턴스(Spot Instances)** 또는 선점형 VM(Preemptible VMs)이라 부른다.

```
       [Cloud Cost & Spot Interruption Optimization Topology]

     On-Demand Node Pool                   Spot Instance Fleet
     (Stable / High Cost)                 (~75% 예시 가정; 상한 90%)
  ┌─────────────────────────┐           ┌─────────────────────────┐
  │ Nextflow Head Node      │           │ AWS Batch Worker Tasks  │
  │ • Workflow Orchestrator │           │ • FastQC / MultiQC      │
  │ • State Tracking (.next)│           │ • Alevin-fry / STARsolo │
  │ • Zero Interruption Req │           │ • AnnData QC Filtering  │
  └────────────┬────────────┘           └────────────┬────────────┘
               │                                     │
               │ HTTP API Orchestration              │ Spot 2-min Warning
               ▼                                     ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ Distributed Shared Object Storage (AWS S3 / FSx Lustre)       │
  │ - Task Work Dirs: s3://bucket/work/ab/c1d2e3...               │
  │ - Task hash: Hash(inputs metadata, script, container, env)    │
  └───────────────────────────────────────────────────────────────┘
```

### 1.1 스팟 인스턴스 중단 확률의 수리적 모델링

클라우드 공급자는 데이터센터의 실시간 온디맨드 수요가 급증할 경우, 스팟 인스턴스를 2분의 유예 경고(2-minute Termination Notice) 후 즉각 회수(Interruption)한다.

노드 회수 사건이 모수 $\lambda$를 갖는 비동차 포아송 과정(Inhomogeneous Poisson Process)을 따른다고 가정할 때, 실행 시간 $t$ 동안 인스턴스가 적어도 한 번 중단될 확률 $P_{\text{interrupt}}(t)$는 다음과 같다:

$$P_{\text{interrupt}}(t) = 1 - \exp\left( -\int_0^t \lambda(u) \, du \right) \approx 1 - e^{-\lambda t}$$

여기서 $\lambda$는 특정 인스턴스 패밀리(예: `c6i.4xlarge`, `r6i.4xlarge`)와 가용 영역(Availability Zone, AZ)의 공급-수요 긴장도에 의해 결정되는 위험률(Hazard Rate)이다.

### 1.2 기대 전산 비용 비교 수리식

작업 실행 시간 $t$, 온디맨드 시간당 단가 $C_o$, 스팟 시간당 단가 $C_s$ ($C_s \approx 0.25 C_o$ (75% 할인 예시 가정)), 스팟 중단 발생 시 손실되는 평균 시간 $t_{\text{lost}}$, 그리고 재시도 시의 추가 비용을 고려할 때, 기대 비용 모델은 다음과 같이 정식화된다:

$$\mathbb{E}[\text{Cost}_{\text{On-Demand}}] = C_o \cdot t$$

$$\mathbb{E}[\text{Cost}_{\text{Spot}}] = C_s \cdot t + P_{\text{interrupt}}(t) \cdot \left( C_s \cdot t_{\text{lost}} + \mathbb{E}[\text{Cost}_{\text{Retry}}] \right)$$

* **임계 조건**:
  단일 작업의 소요 시간 $t$가 40시간 이상으로 극도로 길고 체크포인트가 없다면 $P_{\text{interrupt}} \to 1.0$이 되어 재시도 비용이 절감분을 초과할 수 있다.
  그러나 Nextflow의 경우 모든 개별 프로세스가 수십 분~수 시간 단위의 독립된 태스크로 세분화되어 있고, S3 객체 스토리지에 작업 결과가 즉각 영속화(Persisted)되므로 $\lambda t \ll 1$이 성립하여 스팟 전략이 항상 압도적인 경제적 우위($60\% \sim 75\%$ 순비용 절감)를 점한다.

---

## 2. Nextflow 작업 해시 알고리즘과 `-resume` 캐시 무결성

Nextflow의 무손실 내결함성(Fault-tolerance)과 비용 절감의 핵심은 **고유 작업 디렉토리 해싱(Task Work Directory Hashing)** 알고리즘에 있다.

### 2.1 태스크 해시 계산 공식

Nextflow 엔진은 프로세스가 큐에 디스패치되기 직전, 해당 태스크의 상태 변수(공식 문서의 해시 구성: 세션 ID, 태스크 이름, 컨테이너·conda 환경, 입력, 스크립트, 스크립트가 참조하는 전역 변수 등)를 수집해 해시를 계산한다. 해시 함수의 종류는 공식 문서에 명시되어 있지 않으며, 아래 수식은 이해를 위한 개념식이다:

$$H_{\text{task}} = \operatorname{Hash}\left( S_{\text{script}} \mathbin{\Vert} \sum_{i} H(I_i) \mathbin{\Vert} C_{\text{container}} \mathbin{\Vert} E_{\text{env}} \mathbin{\Vert} V_{\text{params}} \right)$$

1. **$S_{\text{script}}$ (스크립트 바이트열)**: 프로세스의 `script:` 또는 `shell:` 블록에 기술된 쉘 명령 및 템플릿 변수가 치환된 최종 문자열.
2. **$\sum H(I_i)$ (입력 파일 메타데이터)**: 기본(`cache true`)에서는 입력 파일의 이름(경로), 크기(Size), 최종 수정 시각(Last Modified Timestamp)을 사용한다. `cache 'lenient'`는 이름·크기만, `cache 'deep'`은 파일 내용까지 포함한다.
3. **$C_{\text{container}}$ (컨테이너 식별자)**: 실행에 사용된 Docker/Singularity 이미지의 전체 태그 또는 OCI SHA-256 다이제스트.
4. **$E_{\text{env}}$ (실행 환경)**: conda 환경·환경 모듈 등 태스크의 소프트웨어 환경 정보.
5. **$V_{\text{params}}$ (입력 값 파라미터)**: `val(meta)`로 전달된 환자 ID, 샘플 유형 등의 원시 데이터 값.

### 2.2 해시 캐시 적중(Cache Hit) 및 스팟 중단 극복 메커니즘

계산된 16진수 해시 문자열이 `a3b89c...`라면, Nextflow는 작업 디렉토리 경로를 다음과 같이 생성한다:

$$\text{Work Dir} = \text{work}/a3/b89c\dots$$

- **스팟 인스턴스 강제 종료 발생 시**:
  AWS Batch 노드가 중간에 회수되어 작업이 비정상 종료되면(종료 코드는 SIGTERM 기반 143일 수도, 없을 수도 있다), Batch의 스팟 재시도(`aws.batch.maxSpotAttempts`) 또는 Nextflow의 `errorStrategy` 재시도로 해당 태스크가 다시 실행된다.
- **`-resume` 재실행 동작**:
  새로운 스팟 노드에서 파이프라인이 재개될 때, 이미 완료된 선행 태스크들(예: FastQC, STARsolo 등)은 동일한 $H_{\text{task}}$를 산출하므로 작업을 다시 실행하지 않고 작업 디렉토리의 캐시 파일 출력을 즉시 다운스트림 채널로 전달한다.
  따라서 손실되는 비용은 **중단된 그 순간에 실행 중이던 단일 태스크의 미완료 시간 분량**에 불과하며, 이미 완료된 수백 개의 환자 태스크는 단 1센트의 중복 비용도 발생시키지 않는다.

---

## 3. 자원 과다 프로비저닝 진단 및 적응형 라이트사이징(Right-sizing)

생명정보학 파이프라인의 전형적인 비용 낭비는 실패(OOM)를 두려워하여 모든 작업에 최대 사양(예: 128 GB RAM, 32 코어)을 일괄 할당하는 관행에서 비롯된다.

### 3.1 경험적 분위수 기반 자원 최적화 수리 모델

$N$개의 환자 샘플 처리 과정에서 관측된 프로세스 $P$의 물리적 피크 메모리 사용량 표본을 $M = \{m_1, m_2, \dots, m_N\}$이라 하자.
연속형 경험적 분포 함수 $\hat{F}(m)$에 대해 $1-\alpha$ 신뢰 수준(통상 $99\%$, 즉 $\alpha=0.01$)의 분위수 $Q_{1-\alpha}(M)$을 정의한다:

$$Q_{1-\alpha}(M) = \inf \{ m \mid \hat{F}(m) \ge 1-\alpha \}$$

안전 마진(Safety Margin) $\gamma = 1.25$ ($25\%$ 버퍼)를 적용할 때, 1차 시도에 할당해야 할 최적화된 기준 메모리 $M_{\text{alloc}}^{(1)}$은 다음과 같다:

$$M_{\text{alloc}}^{(1)} = \max\left( M_{\text{min}}, \; \lceil Q_{0.99}(M) \times 1.25 \rceil \right)$$

만약 $1\%$의 극단적 대용량 검체(초고심도 시퀀싱 검체)로 인해 메모리 초과가 발생하면, 고정된 자원을 낭비하는 대신 **시도 횟수($\text{task.attempt}$)에 비례하여 동적으로 메모리를 증설**하는 적응형 에러 핸들러를 가동한다:

$$M_{\text{alloc}}^{(k)} = \min\left( M_{\text{max}}, \; M_{\text{alloc}}^{(1)} \times 2^{k-1} \right) \quad \text{for attempt } k \in \{1, 2, \dots, 1+\text{maxRetries}\}$$

---

## 4. 프로덕션 동적 자원 확장 Nextflow 설정

다음은 AWS Batch 및 로컬 HPC 클러스터에서 스팟 중단 재시도(`maxSpotAttempts`)와 메모리 부족(137) 재시도를 구성하는 엔터프라이즈 `nextflow.config` 구현이다:

```groovy
// conf/aws_batch_resources.config

// 1. 최대 허용 자원 한계선 정의 (Nextflow 24.04+ resourceLimits; 구 check_max 함수 대체)
process.resourceLimits = [ cpus: 64, memory: 256.GB, time: 48.h ]

// 3. 프로세스별 스팟 내결함성 및 지능형 재시도 전략
process {
    // 기본 에러 전략: 재시도 3회 허용
    errorStrategy = { task.exitStatus in [143, 137, 104, 134, 139] ? 'retry' : 'finish' }
    maxRetries    = 3
    maxErrors     = '-1'

    // 글로벌 기본 경량 프로파일
    cpus   = { 2 * task.attempt }
    memory = { 4.GB * task.attempt }
    time   = { 2.h * task.attempt }

    // 고연산/고메모리 얼라이너 레이블 (STARsolo 등)
    withLabel: 'process_high_memory' {
        cpus   = 8
        // 1회차: 40GB -> OOM 재시도마다 2배: 80GB -> 160GB -> 4회차 320GB는 resourceLimits로 256GB에 제한 (maxRetries=3이면 최대 4회 시도)
        memory = { 40.GB * (2 ** (task.attempt - 1)) }
        time   = { 6.h * task.attempt }
    }

    // 초경량 준정렬 레이블 (Alevin-fry 등)
    withLabel: 'process_light_quant' {
        cpus   = 4
        memory = { 8.GB * task.attempt }
        time   = { 1.h * task.attempt }
    }
}

// 4. AWS Batch 실행자 프로파일 설정
profiles {
    awsbatch {
        process.executor = 'awsbatch'
        process.queue    = 'arn:aws:batch:us-east-1:123456789012:job-queue/mcrpc-spot-queue'
        aws.region       = 'us-east-1'
        aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
        workDir          = 's3://clinical-genomics-bucket/nextflow-workdir'
    }
}
```

---

## 5. Nextflow 추적 분석 및 클라우드 비용 감사 도구 (`cost_analyzer.py`)

다음은 파이프라인 완료 후 `trace.tsv`를 파싱하여 각 프로세스의 자원 낭비율을 진단하고, 온디맨드 대비 스팟 인스턴스 적용 시 절감된 실제 달러(\$) 금액을 정량 산출하는 분석 스크립트이다.

```python
#!/usr/bin/env python3
"""
cost_analyzer.py
Nextflow trace.tsv 분석을 통한 자원 과다 프로비저닝 진단 및 클라우드 비용 감사 도구.
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path

# vCPU 및 GB-hour 가상 예시 단가 ($/hour). 아래 두 값은 EC2 c6i 단가가 아니라 AWS Fargate(x86) 시간당 단가로 알려진 값(미조회)이며
# 가격은 리전·시점에 따라 변하므로 실제 감사에는 현재 요금표로 교체해야 한다.
PRICE_VCPU_ONDEMAND = 0.04048
PRICE_GB_ONDEMAND   = 0.004445
SPOT_DISCOUNT_RATIO = 0.75      # 가정: 스팟 75% 할인 (Q1과 동일한 예시 가정; AWS 공개 상한은 최대 90%, 실제는 시점·인스턴스별 상이). 재시도 비용은 반영되지 않음

def parse_memory_to_gb(mem_str: str) -> float:
    """Nextflow 메모리 문자열(예: '38.2 GB', '4120 MB')을 기가바이트 실수로 변환한다."""
    if pd.isna(mem_str) or not mem_str:
        return 0.0
    mem_str = str(mem_str).strip()
    parts = mem_str.split()
    if len(parts) != 2:
        return 0.0
    val, unit = float(parts[0]), parts[1].upper()
    if "GB" in unit:
        return val
    elif "MB" in unit:
        return val / 1024.0
    elif "KB" in unit:
        return val / (1024.0 ** 2)
    elif "B" in unit:
        return val / (1024.0 ** 3)
    return val

def parse_duration_to_hours(dur_str: str) -> float:
    """Nextflow 소요 시간 문자열(예: '1h 23m 10s', '45m 12s')을 시간(Hour) 실수로 변환한다."""
    if pd.isna(dur_str) or not dur_str:
        return 0.0
    dur_str = str(dur_str).strip()
    # 밀리초 또는 단순 숫자 처리
    if dur_str.endswith("ms"):
        return float(dur_str[:-2]) / (1000.0 * 3600.0)
    
    total_seconds = 0.0
    tokens = dur_str.split()
    for token in tokens:
        if token.endswith("h"):
            total_seconds += float(token[:-1]) * 3600.0
        elif token.endswith("m"):
            total_seconds += float(token[:-1]) * 60.0
        elif token.endswith("s"):
            total_seconds += float(token[:-1])
    return total_seconds / 3600.0

def main():
    parser = argparse.ArgumentParser(description="Analyze Nextflow Execution Costs and Over-provisioning")
    # 주의: 이 스크립트는 process, status, duration, peak_rss, memory, cpus 열을 사용한다. 기본 trace 필드에 모두 포함된다고 가정하지 말고
    # nextflow.config에 trace.fields = 'task_id,process,status,duration,peak_rss,memory,cpus' 를 명시한다.
    parser.add_argument("--trace", required=True, type=Path, help="Nextflow trace.tsv 경로")
    args = parser.parse_args()

    df = pd.read_csv(args.trace, sep="\t")

    # 성공한 태스크만 필터링
    df_success = df[df["status"] == "COMPLETED"].copy()
    if df_success.empty:
        print("[ERROR] 완료된(COMPLETED) 태스크가 없습니다.")
        sys.exit(1)

    df_success["duration_hours"] = df_success["duration"].apply(parse_duration_to_hours)
    df_success["peak_rss_gb"] = df_success["peak_rss"].apply(parse_memory_to_gb)
    df_success["req_mem_gb"] = df_success["memory"].apply(parse_memory_to_gb)
    df_success["cpus_num"] = pd.to_numeric(df_success["cpus"], errors="coerce").fillna(1)

    # 태스크별 온디맨드 비용 및 스팟 비용 계산
    df_success["cost_ondemand"] = df_success["duration_hours"] * (
        (df_success["cpus_num"] * PRICE_VCPU_ONDEMAND) +
        (df_success["req_mem_gb"] * PRICE_GB_ONDEMAND)
    )
    df_success["cost_spot"] = df_success["cost_ondemand"] * (1.0 - SPOT_DISCOUNT_RATIO)

    total_ondemand = df_success["cost_ondemand"].sum()
    total_spot = df_success["cost_spot"].sum()
    total_saved = total_ondemand - total_spot

    # 프로세스별 자원 낭비 통계
    summary = df_success.groupby("process").agg(
        task_count=("task_id", "count"),
        mean_duration_hr=("duration_hours", "mean"),
        mean_req_mem_gb=("req_mem_gb", "mean"),
        max_peak_rss_gb=("peak_rss_gb", "max"),
        mean_peak_rss_gb=("peak_rss_gb", "mean"),
        total_ondemand_usd=("cost_ondemand", "sum"),
        total_spot_usd=("cost_spot", "sum")
    ).reset_index()

    summary["mem_utilization_pct"] = (summary["mean_peak_rss_gb"] / summary["mean_req_mem_gb"]) * 100.0

    print("================================================================================")
    print("                     NEXTFLOW CLOUD COST & EFFICIENCY AUDIT                    ")
    print("================================================================================")
    print(f"총 분석 태스크 수: {len(df_success)} 개")
    print(f"온디맨드 인스턴스 환산 총 비용: ${total_ondemand:,.2f} USD")
    print(f"스팟 인스턴스(Spot) 실제 총 비용: ${total_spot:,.2f} USD")
    print(f"스팟 활용에 따른 순 절감액:       ${total_saved:,.2f} USD (절감율: {SPOT_DISCOUNT_RATIO*100:.1f}%)")
    print("--------------------------------------------------------------------------------")
    print("프로세스별 자원 효율성 및 메모리 활용률 분석:")
    for _, row in summary.iterrows():
        print(f"\n[Process: {row['process']}] (실행 횟수: {row['task_count']})")
        print(f"  - 평균 실행 시간: {row['mean_duration_hr']:.2f} 시간")
        print(f"  - 할당 메모리: {row['mean_req_mem_gb']:.1f} GB | 실제 평균 피크 RAM: {row['mean_peak_rss_gb']:.1f} GB (최대 피크: {row['max_peak_rss_gb']:.1f} GB)")
        print(f"  - 메모리 활용도: {row['mem_utilization_pct']:.1f}%")
        if row['mem_utilization_pct'] < 25.0:
            rec_mem = max(8.0, row['max_peak_rss_gb'] * 1.25)   # 3.1절의 안전 마진 γ=1.25와 통일 (표본이 작으면 Q0.99 대신 최댓값 사용)
            print(f"  -> [경고: 심각한 자원 과다 프로비저닝] 메모리 요청량을 {row['mean_req_mem_gb']:.1f} GB에서 {rec_mem:.1f} GB로 하향 조정을 권고합니다.")
    print("================================================================================")

if __name__ == "__main__":
    main()
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 500명 전립선암(mCRPC) 환자 코호트의 온디맨드 vs 스팟 인스턴스 총소유비용(TCO) 및 스팟 중단 리스크 계산
임상시험 연구소에서 500명의 mCRPC 환자 코호트(총 500개 샘플)의 snRNA-seq 데이터를 클라우드 인프라에서 분석하고자 한다.
- 샘플당 전처리 및 얼라인먼트 파이프라인 평균 실행 시간: $T = 3.0\text{ 시간}$
- 사용 인스턴스 규격: `r6i.2xlarge` (8 vCPU, 64 GB RAM)
- 가격 조건: 온디맨드 단가 $C_o = \$0.504\text{ / hour}$, 스팟 단가 $C_s = \$0.126\text{ / hour}$ ($75\%$ 할인)
- 과거 가용 영역(AZ) 통계에 따른 3시간 런타임 내 스팟 중단 확률: $P_{\text{interrupt}} = 0.08$ ($8\%$)
- 중단 발생 시 Nextflow `-resume` 캐시 메커니즘을 통해 손실되는 시간은 해당 태스크의 평균 절반($1.5\text{ 시간}$)이며, 재시도 작업은 100% 온디맨드로 성공 완주된다고 가정한다.
1. 500개 샘플 전체를 순수 온디맨드 인스턴스로만 처리할 때의 총 컴퓨팅 비용 $C_{\text{total, ondemand}}$를 계산하시오.
2. 스팟 인스턴스 기반 재시도 전략을 도입했을 때의 기대 총비용 $\mathbb{E}[C_{\text{total, spot}}]$과 최종 절감액(\$)을 산출하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 순수 온디맨드 총 비용 계산**:
* 총 컴퓨팅 시간: $500\text{ 샘플} \times 3.0\text{ 시간} = 1,500\text{ 인스턴스-시간}$
* 온디맨드 총 비용:
  $$C_{\text{total, ondemand}} = 1,500 \times \$0.504 = \mathbf{\$756.00\text{ USD}}$$
  (샘플당 평균 비용: $\$1.512$)

**2. 스팟 전략 하에서의 기대 총비용 및 절감액 계산**:
* 500개 작업 중 성공 및 중단 예상 수:
  - 첫 시도 성공 태스크 수: $500 \times (1 - 0.08) = 460\text{ 개}$
  - 스팟 중단 발생 태스크 수: $500 \times 0.08 = 40\text{ 개}$
* 세부 비용 요소 계산:
  1. 첫 시도 성공 작업 비용 (460개):
     $$C_1 = 460 \times 3.0\text{ 시간} \times \$0.126 = \$173.88$$
  2. 중단된 40개 작업의 손실 시간 비용 (평균 1.5시간 스팟 요금 소모):
     $$C_{\text{lost}} = 40 \times 1.5\text{ 시간} \times \$0.126 = \$7.56$$
  3. 중단 후 재시도 작업 비용 (40개 작업을 온디맨드로 3시간 완주):
     $$C_{\text{retry}} = 40 \times 3.0\text{ 시간} \times \$0.504 = \$60.48$$
* **기대 총비용 ($\mathbb{E}[C_{\text{total, spot}}]$)**:
  $$\mathbb{E}[C_{\text{total, spot}}] = C_1 + C_{\text{lost}} + C_{\text{retry}} = \$173.88 + \$7.56 + \$60.48 = \mathbf{\$241.92\text{ USD}}$$
  (샘플당 평균 비용: $\approx \$0.484$)
* **순 절감액**:
  $$\text{총 절감액} = \$756.00 - \$241.92 = \mathbf{\$514.08\text{ USD}} \quad (\mathbf{68.0\%}\text{ 순비용 절감})$$
* **결론**: 8%의 스팟 노드 회수가 발생하더라도(참고: AWS Spot Instance Advisor는 월간 중단 빈도를 <5%, 5–10%, 10–15%, 15–20%, >20% 구간으로 제공하므로 3시간 내 8%는 높은 편인 스트레스 시나리오 가정) Nextflow 태스크 캐싱 덕분에 온디맨드 대비 68%의 비용을 절감할 수 있다.

</details>

---

### 문제 2: Nextflow 작업 해시 계산 시 스크립트 공백 수정 및 파일 타임스탬프 왜곡으로 인한 거짓 캐시 미스(False Cache Miss)
수백 명의 환자 유전체 데이터를 정렬 완료한 후, 연구원이 보고서 출력을 개선하기 위해 `SCRNASEQ_FILTER_ANNDATA` 모듈의 bash 스크립트 내부에서 단순 주석(`# Filter cells`) 1줄을 추가하고 들여쓰기를 정돈하였다. 이후 `nextflow run main.nf -resume`을 실행했으나, 파이프라인이 이전 결과를 재사용하지 않고 100개 샘플 전체에 대해 40시간에 걸친 재연산을 시작해버리는 치명적인 클라우드 과금 참사가 발생했다.
1. Nextflow가 작업 디렉토리 고유 해시 $H_{\text{task}}$를 계산할 때 어떤 수학적 구성 요소를 포함하기에 주석 변경만으로 캐시가 무효화(Cache Invalidation)되었는지 원인을 규명하시오.
2. 클라우드 S3 버킷에서 로컬로 파일을 동기화할 때, 파일 내용의 1바이트도 변하지 않았음에도 불구하고 `rsync`나 `aws s3 sync`의 어떤 동작 특성이 Nextflow 캐시 미스를 유발할 수 있는지 설명하고 이를 방지하는 `nextflow.config` 튜닝 기법을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 주석 수정에 의한 캐시 무효화 원인 규명**:
* Nextflow의 작업 해시 개념식은 $H_{\text{task}} = \operatorname{Hash}(S_{\text{script}} \mathbin{\Vert} \dots)$이다.
* 여기서 $S_{\text{script}}$는 템플릿 변수가 치환된 후 실행될 **쉘 스크립트의 전체 텍스트 바이트 스트림**이다. Nextflow 엔진은 코드의 구문론적 의미(Abstract Syntax Tree)를 파싱하여 주석과 실행 코드를 구별하지 않고, 문자열 전체의 바이트 단위 해시를 계산한다.
* 주석 한 줄 추가 및 공백 변경만으로도 $S_{\text{script}}$의 이진 데이터가 변하고, 해시 함수의 눈사태 효과(Avalanche Effect)로 인해 완전히 상이한 해시 다이제스트(예: `work/7f/2a...` $\to$ `work/e4/1b...`)가 생성된다. 결과적으로 새 작업 디렉토리가 생성되어 선행 결과물이 완전히 무시되고 전면 재실행된 것이다.

**2. S3 동기화 시 파일 타임스탬프 왜곡 및 해결책**:
* **캐시 미스 원인**: Nextflow는 입력 파일의 해시 $\sum H(I_i)$를 계산할 때 대용량 유전체 파일 전체를 읽는 부하를 피하기 위해 파일의 `(Path, Size, LastModifiedTimestamp)` 튜플을 기본 검사 대상으로 삼는다. `aws s3 sync`나 `scp`로 파일을 복사하면 파일 내용이 동일하더라도 복사 방식(옵션)에 따라 로컬 파일의 수정 시각(mtime)이 원본과 달라질 수 있다. 타임스탬프가 달라졌으므로 Nextflow는 입력 데이터가 변경된 것으로 오판하여 캐시 미스를 선언한다.
* **방지 대책 (`nextflow.config`)**:
  `cache`는 프로세스 디렉티브이므로 `nextflow.config`의 `process` 스코프에서 지정한다:
  ```groovy
  process.cache = 'lenient'   // 또는 'deep'
  ```
  - `cache = 'deep'`: 파일 내용을 캐시 키에 포함한다(타임스탬프에는 영향받지 않으나 대용량 파일 읽기 비용이 든다).
  - `cache = 'lenient'`: 파일 크기(Size)와 이름만 일치하고 타임스탬프가 더 최신이더라도 캐시 적중을 허용한다.

</details>

---

### 문제 3: 임상 긴급도(TAT < 24시간)와 예산 제약(\$15/sample)을 동시 만족하는 하이브리드 HPC-Cloud 버스팅(Cloud Bursting) 설계
신약 임상 1b/2상 시험 중 전이성 전립선암 환자가 급성 질병 진행(Rapid Progression)을 보여 24시간 이내에 전사체 분석을 완수하여 표적 치료제 투여 여부를 결정해야 하는 응급 검체가 정기적으로 유입된다.
- 원내 온프레미스 HPC 클러스터: 128코어, 상시 무료 이용 가능, 그러나 타 연구자들의 대기열(Queue Wait Time)이 평균 18시간 발생함.
- AWS 클라우드 인프라: 즉각적인 오토스케일링 가능(대기열 0분), 그러나 온디맨드 단가가 고가임.
1. 임상 검체의 턴어라운드 시간(TAT = Queue Wait Time + Execution Time)을 24시간 이내로 보증하면서도 평시 연구 예산(\$15/샘플)을 초과하지 않도록 하는 **우선순위 기반 하이브리드 스케줄링 라우팅 규칙**을 수립하시오.
2. Nextflow의 복수 엑시큐터(`executor = 'slurm'` vs `executor = 'awsbatch'`) 프로파일 동적 전환 전략을 Groovy 조건문으로 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 우선순위 기반 하이브리드 클라우드 버스팅 라우팅 규칙**:
* **검체 분류 (Triage)**:
  - **Class A (임상 긴급 검체, Urgent TAT < 24h)**: 표적 치료제 결정용 mCRPC 바늘 생검 검체.
  - **Class B (정기 연구 검체, Batch TAT < 7일)**: 후향적 바이오뱅크 및 모델링용 검체.
* **라우팅 알고리즘**:
  1. Class B 검체는 100% 온프레미스 HPC(SLURM)로 라우팅하여 클라우드 비용을 \$0으로 유지.
  2. Class A 검체 유입 시 온프레미스 대기 시간($T_{\text{queue}}$)을 쿼리:
     - 만약 $T_{\text{queue}} + T_{\text{exec}} \le 18\text{시간}$(안전 마진 포함 임계, 24시간 TAT 대비 6시간 여유)이면 온프레미스 HPC의 최상위 우선순위 큐(`qos=urgent`)로 디스패치.
     - 그렇지 않으면(18시간 초과 — 구간 누락 없이 단일 임계로 정의) 즉각 AWS Cloud로 버스팅(Cloud Bursting) 실행.
  3. 클라우드 버스팅 시에도 기본 정렬은 고성능 스팟 플릿으로 처리하고 실패 시에만 온디맨드로 승격시켜 샘플당 클라우드 지출을 \$12 이내로 통제.

**2. Nextflow 조건부 엑시큐터 프로파일 구현**:
(주의: 이는 태스크별 동적 전환이 아니라 **실행(run) 단위** 전환이다. `awsbatch` 실행은 `workDir`가 S3여야 하고 입력 데이터가 클라우드에 있어야 하며, 임상 검체(PHI)의 클라우드 전송에는 기관의 데이터 거버넌스·BAA 등 규정 검토가 선행되어야 한다.)
```groovy
// nextflow.config 내의 동적 프로파일 라우팅
params.urgent_clinical = false // 실행 시 --urgent_clinical 플래그로 제어

process {
    executor = params.urgent_clinical ? 'awsbatch' : 'slurm'
    queue    = params.urgent_clinical ? 'arn:aws:batch:...:job-queue/urgent-spot-queue' : 'normal'

    // 긴급 검체일 경우 클라우드 다중 코어 고속 병렬 모드로 확장
    cpus   = params.urgent_clinical ? 16 : 8
    memory = params.urgent_clinical ? 64.GB : 32.GB

    errorStrategy = { task.exitStatus in [143, 137] ? 'retry' : 'finish' }
    maxRetries    = 2
}
```
* 이 설계를 통해 긴급 환자 검체는 대기 시간 0분으로 AWS Batch에 즉각 투입되어 4~6시간 만에 h5ad 분석 및 바이오마커 리포트 발행을 완주(TAT $\approx 5$시간)하면서도, 전체 연간 분석 예산을 안전하게 준수할 수 있다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Nextflow 24.04.4 (nextflow.io). 원 논문: Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35(4):316-319. PMID: 28398311. doi: 10.1038/nbt.3820
- Slurm (slurm.schedmd.com).
- FastQC (github.com/s-andrews/FastQC).
- MultiQC (github.com/MultiQC/MultiQC). 원 논문: Ewels P, Magnusson M, Lundin S, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016;32(19):3047-8. PMID: 27312411. doi: 10.1093/bioinformatics/btw354
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
<!-- REFERENCES:END -->

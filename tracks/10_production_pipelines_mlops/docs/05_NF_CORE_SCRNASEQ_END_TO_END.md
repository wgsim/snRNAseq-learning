# 05. nf-core/scrnaseq 엔드투엔드 파이프라인 분석 및 커스텀 모듈 엔지니어링 (End-to-End nf-core/scrnaseq Pipeline & Custom Module Engineering)

---

## 학습 목표 (Cognitive Goals)

1. **골드 스탠다드 `nf-core/scrnaseq` 아키텍처의 채널 토폴로지 분석**: Nextflow DSL2 기반의 멀티오믹스 파이프라인 구조를 분해하고, FastQC/MultiQC 입력 큐, 전처리 채널 라우팅, 그리고 업스트림 얼라이너 모듈로의 메타데이터(`meta` 맵) 전파 메커니즘을 규명한다.
2. **단일세포 얼라이너 모듈 트레이드오프와 벤치마크 수학적 정식화**: STARsolo, Alevin-fry, Cell Ranger 간의 메모리 복잡도 $O(N)$, 계산 시간 복잡도, 슈도얼라인먼트(Pseudoalignment) vs 스플라이스 인식 전장 정렬(Splice-aware Full Alignment)의 알고리즘적 차이를 정량화한다.
3. **nf-core 호환 커스텀 모듈 설계 및 소프트웨어 추적성 보장**: `conda`, Docker/Singularity OCI 컨테이너 지시자, 입출력 튜플 규격화, 그리고 다중 도구 버전 정보를 단일 YAML 스트림으로 방출하는 `versions.yml` 추적 패턴을 완벽히 구현한다.
4. **전이성 거세저항성 전립선암(mCRPC) 다중 환자 코호트 샘플시트 설계 및 자동화 QC**: 골 및 내장 전이 생검 검체의 이질성을 반영하는 샘플시트(`sample,fastq_1,fastq_2`; 임상 공변량은 별도 메타데이터 파일로 분리)를 구축하고, 미토콘드리아 RNA 비율 및 세포 생존도 기반의 자동화 AnnData/Seurat 필터링 파이프라인을 구축한다.

---

## 1. nf-core/scrnaseq 파이프라인 아키텍처 및 채널 라우팅

`nf-core/scrnaseq`는 전 세계 생명정보학 커뮤니티가 공동 개발하고 유지보수하는 검증된 단일세포 전사체(scRNA-seq/snRNA-seq) 업스트림 프로세싱 표준 파이프라인이다. 수백 개의 환자 검체로부터 생산된 테라바이트 규모의 FASTQ 데이터를 안정적이고 재현 가능하게 처리할 수 있도록 Nextflow DSL2 모듈형 서브워크플로우(Subworkflows) 체계로 구성되어 있다.

```
                     [Clinical mCRPC FASTQ Samplesheet]
                                     │
                                     ▼
                      ┌──────────────────────────────┐
                      │    SAMPLESHEET_CHECK Module  │
                      └──────────────┬───────────────┘
                                     │ Channel: [meta, [r1, r2]]
            ┌────────────────────────┼────────────────────────┐
            ▼                                                 ▼
┌───────────────────────┐                         ┌───────────────────────┐
│     FASTQC Module     │                         │  Aligner Subworkflow  │
│  (Read Quality Base)  │                         │ (STARsolo/Alevin/CR)  │
└───────────┬───────────┘                         └───────────┬───────────┘
            │                                                 │
            │ qc_files                                        │ BAM + Raw Matrix
            ▼                                                 ▼
┌───────────────────────┐                         ┌───────────────────────┐
│    MULTIQC Module     │                         │ MTX_TO_H5AD / SEURAT  │
│ (Cohort-level Report) │                         │  (AnnData QC Filter)  │
└───────────────────────┘                         └───────────┬───────────┘
                                                              │
                                                              ▼
                                                  [Final Filtered AnnData]
                                                  (Clinical Biomarker Ready)
```

### 1.1 입력 채널 생성 및 `meta` 맵 전파 메커니즘

nf-core 아키텍처의 핵심 원칙은 **데이터 페이로드(Files)와 메타데이터 컨텍스트(Map)의 원자적 결합(Atomic Tuple Coupling)**이다. 파이프라인은 CSV 샘플시트로부터 레코드를 읽어들인 뒤 다음과 같은 형태의 채널 요소를 생성한다:

$$\text{ch\_reads} = \left[ \mathbf{meta}, [ \text{fastq\_1}, \text{fastq\_2} ] \right]$$

여기서 $\mathbf{meta}$는 불변(Immutable) 해시 맵으로 정의되며, 다운스트림의 모든 프로세스가 샘플의 임상적/기술적 정체성을 즉각 식별할 수 있도록 보장한다:

```groovy
meta = [
    id:           "mCRPC_Patient042_BoneMet",
    single_end:   false,
    strandedness: "forward",
    patient_id:   "P042",
    biopsy_site:  "bone_marrow",
    therapy:      "enzalutamide_resistant",
    chemistry:    "10x_v3"
]
```

이 구조 덕분에 특정 환자의 FASTQ 리드가 복수의 병렬 얼라이너 모듈로 분기되거나 다중 QC 프로세스로 전달되더라도, 결과 파일들이 서로 섞이거나(Cross-contamination) 메타데이터가 유실되는 위험이 원천적으로 차단된다.

### 1.2 FastQC 및 MultiQC를 통한 코호트 품질 관리

원시 시퀀싱 리드는 얼라인먼트 단계로 진입하기 전, 병렬 비동기 채널을 통해 `FASTQC` 모듈로 전달된다:
1. **Per-base Sequence Quality**: 각 시퀀싱 사이클별 Phred 점수($Q$)의 분포를 박스플롯으로 집계한다.
2. **Per-base Sequence Content**: 10x Chromium 3' v3 Read 1의 경우 첫 16bp(Cell Barcode) 및 다음 12bp(UMI)에서 합성 특유의 염기 편향이 확인되는지, Read 2의 시작부에서 TSO(Template Switching Oligo) 아티팩트가 존재하는지 검증한다.
3. **Adapter Content**: cDNA 인서트 길이가 짧아 3' 말단에서 Illumina TruSeq 어댑터나 Poly-A 서열이 리드 내부로 읽혀 들어갔는지 확인한다.

개별 샘플의 FastQC zip 아카이브, 얼라이너 통계 요약표(Alignment summary metrics), 유전자 검출 통계는 수집 연산자(`.collect()`)를 통해 단일 채널로 모인 후 `MULTIQC` 모듈에 주입되어 단일 상호작용형 HTML 리포트로 융합된다.

---

## 2. 얼라이너 모듈 선택 및 알고리즘적 트레이드오프

단일세포 전사체 분석의 전산 비용과 메모리 요구량의 80% 이상은 FASTQ 리드를 참조 유전체에 정렬하고 유전자별 UMI 카운트 매트릭스를 생성하는 업스트림 단계에서 발생한다. `nf-core/scrnaseq`는 세 가지 주요 얼라이너 엔진(`--aligner`)을 지원하며, 각 엔진은 고유한 수학적/알고리즘적 특성을 지닌다.

| 얼라이너 엔진 | 정렬 알고리즘 | 메모리 요구량 (Peak RAM) | 연산 속도 (1억 리드 기준) | 인트론 매핑 (snRNA-seq) |
| :--- | :--- | :--- | :--- | :--- |
| **STARsolo** | Suffix Array 접두사 매칭 + 동적 계획법 | $\approx 32\text{ GB} \sim 40\text{ GB}$ (인간 유전체) | $\approx 40\text{ 분} \sim 1\text{ 시간}$ | 완벽 지원 (`--soloFeatures Gene GeneFull`) |
| **Alevin-fry** | K-mer 슈도얼라인먼트 + RAD 파일 그래프 분해 | $\approx 4\text{ GB} \sim 8\text{ GB}$ | $\approx 15\text{ 분} \sim 25\text{ 분}$ | Splici (Spliced + Intron) 인덱스 활용 |
| **Cell Ranger** | STAR 정렬 래퍼 + 독자적 UMI 그래프 보정 | $\approx 32\text{ GB} \sim 64\text{ GB}$ | $\approx 1.5\text{ 시간} \sim 2.5\text{ 시간}$ | `include-introns true` 지원 |

(위 표의 메모리·속도 수치는 출처 없는 대략적 예시값이며 데이터·파라미터·하드웨어에 따라 크게 다릅니다.)

### 2.1 STARsolo의 메모리 복잡도와 Asymptotics

STARsolo는 접두사 배열(Suffix Array, SA)을 메모리에 상주시킨 상태에서 Read 2의 전사체 서열을 탐색한다. 길이 $G$인 인간 참조 유전체($G \approx 3.1 \times 10^9\text{ bp}$)에 대해 Suffix Array 인덱스의 크기는 다음과 같이 결정된다:

$$M_{\text{SA}} = G \times \lceil \log_2 G \rceil \text{ bits} \approx 3.1 \times 10^9 \times 32\text{ bits} \approx 12.4\text{ GB}$$

위 식은 한 가닥 기준의 단순 추정이며 STAR는 양 가닥($2G$)의 접미사를 담으므로 실제 인덱스는 더 큽니다. STAR 매뉴얼의 경험칙은 **게놈 크기의 약 10배 RAM(인간 게놈 약 30 GB)**입니다(본 문서에서 재확인하지 않음). 여기에 유전체 서열 자체($G/4\text{ bytes} \approx 0.78\text{ GB}$)와 접두사 해시 테이블($4^{p} \times 4\text{ bytes}$, $p=14$일 때 약 $1\text{ GB}$), 그리고 세포 바코드 및 UMI 중복 제거 해시 맵이 추가되므로, 런타임 최소 메모리는 약 $32\text{ GB}$ 이상이 필요하다 ($O(G)$ 공간 복잡도).

### 2.2 Alevin-fry의 초경량화 메모리 기전 (RAD 중간 표현형)

반면 Alevin-fry는 전장 정렬 대신 서열 결정적 유사성 매핑(Selective Alignment/Pseudoalignment)을 수행하여 리드를 직접 BAM으로 쓰지 않고 경량화된 이진 포맷인 **RAD(Reduced Alignment Data)** 파일로 스트리밍한다:

1. **K-mer 준정렬**: 전사체 인덱스(Splici)를 해시 테이블로 구축하여 $O(K \cdot L)$ 시간 복잡도로 동등 클래스(Equivalence Class)를 추출한다.
2. **2단계 파이프라인 분리**:
   - `salmon alevin`: FASTQ로부터 세포 바코드, UMI, 유전자 동등 클래스 식별자만을 추출하여 RAD 파일로 직렬화 (RAM $< 4\text{ GB}$).
   - `alevin-fry quant`: RAD 파일을 메모리 매핑(mmap)하여 세포별 바코드 화이트리스트 보정 및 유전자 단위 UMI 그래프 중복 제거(Graph Deduplication)를 인메모리에서 수행 (RAM $< 4\text{ GB}$).

따라서 대규모 mCRPC 환자 코호트(예: 100명 이상의 생검 검체)를 AWS Batch나 온프레미스 HPC 클러스터에서 고병렬로 동시 처리할 때, Alevin-fry는 노드당 동시 실행 가능 작업 수($N_{\text{tasks}} = \frac{\text{Node RAM}}{\text{Task RAM}}$)를 늘릴 수 있다(표의 예시 RAM 값을 쓰면 문제 1의 계산처럼 약 2.7배; 4~8배는 출처 없는 상한 추정) 파이프라인 전체 완료 시간(Turnaround Time, TAT)과 클라우드 컴퓨팅 비용을 극적으로 절감시킨다.

---

## 3. 다운스트림 자동화 AnnData 생성 및 미토콘드리아 QC

정렬 및 분자 카운팅이 완료된 원시 희소 행렬(Raw Feature-Barcode Matrix, MEX 포맷: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`)은 파이프라인 내에서 즉각 표준 단일세포 전사체 구조체(AnnData `.h5ad` 또는 Seurat `.rds`)로 변환되고 1차 필터링이 수행된다.

### 3.1 미토콘드리아 및 리보솜 전사체 분율 수학적 모델

(전체 세포 scRNA-seq 시나리오) 전립선암 전이 생검 조직(특히 뼈 및 간 전이)은 바늘 생검(Core Needle Biopsy) 및 효소적 해리 과정에서 기계적 전단력(Shear Stress)으로 인해 세포막이 부분 파열된 사멸 세포(Apoptotic/Necrotic Cells)가 다량 발생한다. 세포막이 손상되면 세포질 RNA(Cytoplasmic mRNA)는 외부로 유출되고 상대적으로 견고한 미토콘드리아 막 내부의 전사체만 잔존하여 비정상적으로 높은 미토콘드리아 리드 비율이 관측된다.

세포 $i$의 총 검출 유전자 집합을 $\mathcal{G}$, 미토콘드리아 염색체에 암호화된 유전자 집합을 $\mathcal{G}_{\text{MT}} = \{g \in \mathcal{G} \mid \text{symbol}(g) \text{ starts with } \text{"MT-"}\}$라 할 때, 미토콘드리아 분율 $f_{\text{MT}}^{(i)}$는 다음과 같이 정의된다:

$$f_{\text{MT}}^{(i)} = \frac{\sum_{g \in \mathcal{G}_{\text{MT}}} X_{i, g}}{\sum_{g \in \mathcal{G}} X_{i, g}} = \frac{N_{\text{MT}}^{(i)}}{N_{\text{total}}^{(i)}}$$

통상 세포주(Cell line)나 신선 절제 조직에서는 $f_{\text{MT}} < 10\%$가 표준적인 통과 임계치이나(주의: 이 트랙의 snRNA-seq에서는 핵 분리 과정에서 세포질 mRNA가 제거되어 $f_{\text{MT}}$가 매우 낮고, 루트 3장처럼 약 2% 안팎(조직·프로토콜 의존 예시)을 상한으로 두며 높은 값은 세포질 오염/핵 분리 불량의 지표로 해석한다), 전이성 전립선암 생검 조직의 경우 종양 미세환경 내의 심각한 저산소증(Hypoxia) 및 괴사율을 감안하여 중앙값 절대편차(Median Absolute Deviation, MAD) 기반의 적응형 이상치 제거 모델을 적용한다:

$$\text{Threshold}_{\text{MT}} = \operatorname{median}(f_{\text{MT}}) + 3 \times \operatorname{MAD}(f_{\text{MT}})$$

$$\operatorname{MAD}(f_{\text{MT}}) = \operatorname{median}\left( \left| f_{\text{MT}} - \operatorname{median}(f_{\text{MT}}) \right| \right)$$

이를 통해 전립선 기저/내강 상피세포(Basal/Luminal Epithelial), 신경내분비 세포(NEPC), 침윤 T세포, 거핵세포 간의 생물학적 고유 발현 차이를 보존하면서 기술적 괴사 세포만을 정확히 배제한다.

---

## 4. nf-core 호환 커스텀 모듈 작성 엔지니어링

임상 현장에서 상용화된 `nf-core/scrnaseq` 표준 파이프라인에 연구기관 고유의 알고리즘(예: *AR-V7* 대체 스플라이싱 검출 또는 초기 품질 기반 자동 더블렛 필터링 모듈)을 추가하려면, nf-core 커뮤니티의 엄격한 모듈 가이드라인을 준수해야 한다.

### 4.1 nf-core 모듈의 3대 불변 법칙

1. **격리된 소프트웨어 환경 선언**: `conda` 환경 파일(`environment.yml`)과 Docker Hub/Quay.io의 OCI 다이제스트(`container`)를 반드시 병행 명시해야 한다.
2. **`meta` 맵 입출력 불변성**: 프로세스로 진입한 `meta` 튜플은 변경되지 않고 그대로 출력 채널로 방출되어야 한다.
3. **버전 추적성 (`versions.yml`)**: 모듈 내부에서 실행된 모든 CLI 도구의 바이너리 버전을 추출하여 단일 YAML 파일로 기록해야 한다.

### 4.2 프로덕션 커스텀 모듈 구현: `SCRNASEQ_FILTER_ANNDATA`

다음은 원시 카운트 매트릭스를 수신하여 빈 액적 제거, 미토콘드리아 QC 필터링, 그리고 h5ad 아티팩트를 방출하는 프로덕션 레벨의 Nextflow DSL2 모듈 정의이다:

```groovy
// modules/local/scrnaseq_filter_anndata.nf

process SCRNASEQ_FILTER_ANNDATA {
    tag "$meta.id"
    label 'process_medium'

    // scanpy가 anndata·numpy·scipy·pandas를 함께 설치한다. biocontainers에 존재하는 최신 scanpy 이미지는 1.7.2이므로 conda와 컨테이너 버전을 일치시킴
    conda "conda-forge::scanpy=1.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(raw_matrix_dir)

    output:
    tuple val(meta), path("*.filtered.h5ad"), emit: h5ad
    tuple val(meta), path("*.qc_metrics.csv"), emit: metrics
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'EOF'
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import sys
    from scipy.stats import median_abs_deviation

    # 1. 원시 희소 행렬 로드 (10x MEX 포맷)
    adata = sc.read_10x_mtx(
        "${raw_matrix_dir}",
        var_names='gene_symbols',
        cache=False
    )
    adata.var_names_make_unique()

    # 2. 메타데이터 주입 (meta map 동기화)
    adata.obs['sample_id'] = "${meta.id}"
    adata.obs['patient_id'] = "${meta.patient_id}"
    adata.obs['biopsy_site'] = "${meta.biopsy_site}"

    # 3. 미토콘드리아 및 리보솜 전사체 검출
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'rps', 'rpl'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # 4. 적응형 MAD 기반 이상치 판정 (pandas Series.mad()는 평균절대편차이며 pandas 2.0에서 제거됨 →
    #    scipy의 median_abs_deviation(scale='normal') 사용)
    def is_outlier(df, metric, nmads):
        M = df[metric]
        mad = median_abs_deviation(M, scale='normal')
        return (M < M.median() - nmads * mad) | (M > M.median() + nmads * mad)

    adata.obs['log1p_total_counts'] = np.log1p(adata.obs['total_counts'])
    adata.obs['log1p_n_genes'] = np.log1p(adata.obs['n_genes_by_counts'])
    mad_outlier = (
        is_outlier(adata.obs, 'log1p_total_counts', 5) |
        is_outlier(adata.obs, 'log1p_n_genes', 5)
    )

    # 임상 전립선암 생검 기준 엄격한 하한선(고정값) + 적응형 MAD 이상치 제거를 함께 적용
    min_genes = 300
    min_counts = 500
    max_mt_pct = 20.0

    valid_mask = (
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['pct_counts_mt'] <= max_mt_pct) &
        (~mad_outlier)
    )

    qc_summary = pd.DataFrame({
        'total_cells_raw': [adata.n_obs],
        'total_cells_filtered': [int(valid_mask.sum())],
        'median_genes': [float(adata.obs['n_genes_by_counts'].median())],
        'median_counts': [float(adata.obs['total_counts'].median())],
        'median_mt_pct': [float(adata.obs['pct_counts_mt'].median())]
    })
    qc_summary.to_csv("${prefix}.qc_metrics.csv", index=False)

    # 5. 필터링된 객체 저장
    adata_filtered = adata[valid_mask].copy()
    adata_filtered.write_h5ad("${prefix}.filtered.h5ad", compression="gzip")
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
```

---

## 5. 임상 mCRPC 다중 환자 코호트 샘플시트 설계 및 실행 전략

실제 정밀 종양학(Precision Oncology) 임상 연구 환경에서는 다양한 전이 부위(골수, 간, 림프절)와 치료 이력(엔잘루타마이드, 아비라테론 내성)을 가진 환자 코호트가 연속적으로 유입된다.

### 5.1 프로덕션 샘플시트 명세 (`samplesheet.csv`)

nf-core/scrnaseq의 샘플시트 스키마(`assets/schema_input.json`, 릴리스 2.5.1 기준)는 `sample`, `fastq_1`, `fastq_2`와 선택 열(`expected_cells`, `seq_center` 등)만 허용하며, `strandedness`는 nf-core/rnaseq의 열이다. 환자·조직·치료 이력 같은 임상 공변량을 샘플시트에 직접 넣으면 스키마 검증과 충돌할 수 있으므로, 샘플시트는 규격 열만 두고 공변량은 `sample`을 키로 하는 별도 파일로 관리한다:

```csv
sample,fastq_1,fastq_2
mCRPC_P001_Bone,s3://clinical-genomics-bucket/mCRPC/P001/P001_Bone_R1.fastq.gz,s3://clinical-genomics-bucket/mCRPC/P001/P001_Bone_R2.fastq.gz
mCRPC_P001_Liver,s3://clinical-genomics-bucket/mCRPC/P001/P001_Liver_R1.fastq.gz,s3://clinical-genomics-bucket/mCRPC/P001/P001_Liver_R2.fastq.gz
mCRPC_P002_Lymph,s3://clinical-genomics-bucket/mCRPC/P002/P002_Lymph_R1.fastq.gz,s3://clinical-genomics-bucket/mCRPC/P002/P002_Lymph_R2.fastq.gz
mCRPC_P003_Bone,s3://clinical-genomics-bucket/mCRPC/P003/P003_Bone_R1.fastq.gz,s3://clinical-genomics-bucket/mCRPC/P003/P003_Bone_R2.fastq.gz
```

```csv
sample,patient_id,tissue_type,therapy_status
mCRPC_P001_Bone,P001,bone_marrow,enzalutamide_resistant
mCRPC_P001_Liver,P001,liver_metastasis,enzalutamide_resistant
mCRPC_P002_Lymph,P002,lymph_node,abiraterone_naive
mCRPC_P003_Bone,P003,bone_marrow,taxane_chemo_refractory
```

### 5.2 프로덕션 실행 CLI 및 프로파일 조합

임상 코호트 분석 시, 재현성과 추적성을 확보하기 위해 GitHub 태그로 고정된 파이프라인 릴리스 버전을 사용하고 AWS Batch 또는 로컬 고성능 Docker 엔진 프로파일을 활성화하여 실행한다:

```bash
# 2.5.1의 --aligner 유효값: kallisto, star, alevin, cellranger, cellrangerarc, universc
# (4.x 이후 'alevin'은 'simpleaf'로 대체: kallisto, star, simpleaf, cellranger, ...)
nextflow run nf-core/scrnaseq \
    -r 2.5.1 \
    -profile docker \
    --input ./cohort_mcrpc_samplesheet.csv \
    --outdir s3://clinical-genomics-bucket/results/mCRPC_Cohort_202610 \
    --genome GRCh38 \
    --aligner alevin \
    --protocol 10XV3 \
    -with-report s3://clinical-genomics-bucket/reports/execution_report.html \
    -with-trace s3://clinical-genomics-bucket/reports/trace.tsv \
    -with-timeline s3://clinical-genomics-bucket/reports/timeline.html \
    -resume
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: STARsolo와 Alevin-fry의 메모리 복잡도 격차 및 클라우드 배치 병렬 처리 최적화
전립선암 전이 환자 코호트 120명의 snRNA-seq 라이브러리(샘플당 평균 1억 2천만 리드 쌍)를 AWS 클라우드 환경에서 처리해야 한다. 현재 가용한 단일 클라우드 가상 머신(EC2 인스턴스) 사양은 64개 vCPU와 128 GB RAM을 탑재하고 있다.
1. 인간 유전체 인덱스(GRCh38) 기반 정렬 시, STARsolo(작업당 피크 RAM 36 GB, 8 vCPU 할당)와 Alevin-fry(작업당 피크 RAM 6 GB, 8 vCPU 할당)를 사용할 때 해당 단일 머신에서 안전하게 동시 실행(Concurrent Tasks) 가능한 최대 프로세스 수를 각각 계산하시오.
2. 두 방식 간의 메모리 병목 원인을 접두사 배열(Suffix Array) 상주 방식과 RAD 파일 스트리밍 방식의 알고리즘적 차이 관점에서 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 동시 실행 가능 최대 프로세스 수 계산**:
* **STARsolo**:
  * CPU 기준 병렬 한계: $\lfloor 64 \text{ vCPU} / 8 \text{ vCPU} \rfloor = 8\text{ 작업}$
  * 메모리 기준 병렬 한계: $\lfloor 128 \text{ GB} / 36 \text{ GB} \rfloor = 3\text{ 작업}$
  * 병목 제약 조건: $\min(8, 3) = \mathbf{3\text{ 작업}}$. 8개의 vCPU를 갖는 3개 작업 실행 시 CPU는 24개만 사용되고 40개 vCPU가 유휴(Idle) 상태로 방치되어 하드웨어 낭비율이 극심하다.
* **Alevin-fry**:
  * CPU 기준 병렬 한계: $\lfloor 64 \text{ vCPU} / 8 \text{ vCPU} \rfloor = 8\text{ 작업}$
  * 메모리 기준 병렬 한계: $\lfloor 128 \text{ GB} / 6 \text{ GB} \rfloor = 21\text{ 작업}$
  * 병목 제약 조건: $\min(8, 21) = \mathbf{8\text{ 작업}}$. 8개 프로세스가 64 vCPU를 100% 포화시키며 총 필요 메모리는 $8 \times 6\text{ GB} = 48\text{ GB}$로 128 GB RAM 용량 이내에서 OOM 위험 없이 안전하게 병렬 완주된다.
  * **결론**: Alevin-fry는 동일 하드웨어 노드에서 동시 처리량을 $2.67\text{배}$ ($3 \to 8$) 향상시킨다.

**2. 알고리즘적 메모리 병목 원인 규명**:
* **STARsolo**: 인간 유전체(약 3.1 Gb) 전체에 대한 접두사 배열(Suffix Array, SA) 및 염색체 서열 전체를 공유 메모리 또는 프로세스 전용 힙에 상주시켜야 하므로 본질적으로 $O(G)$의 거대한 메모리 풋프린트($\ge 32\text{ GB}$)가 고정된다.
* **Alevin-fry**: 전장 염색체 인덱스를 사용하지 않고, 전사체 서열(Splici Index)에 대한 K-mer 해시 기반 준정렬을 수행하여 매핑 정보를 즉각 경량화된 중간 표현형인 RAD(Reduced Alignment Data) 바이너리로 디스크에 스트리밍 방출한다. 대용량 BAM 정렬 연산을 건너뛰고 UMI 중복 제거를 RAD 파일 메모리 매핑(mmap)을 통해 인메모리 국소 영역에서 처리하므로 RAM 사용량이 $4 \sim 6\text{ GB}$ 수준으로 극히 작다.

</details>

---

### 문제 2: nf-core 커스텀 모듈 작성 시 `versions.yml` 누락이 유발하는 규제 감사(Regulatory Audit) 위험 및 해결책
임상시험 소프트웨어 적격성 평가(FDA 21 CFR Part 11 및 CLIA/CAP 가이드라인 등 규제 요건 대응)를 준비 중인 바이오텍 연구소에서 `nf-core/scrnaseq` 파이프라인 내부에 독자적인 전립선암 *AR-V7* 스플라이싱 분석 모듈을 Nextflow DSL2로 구현하여 배포하였다. 그러나 해당 모듈의 `script:` 블록에서 소프트웨어 버전을 방출하는 `versions.yml` 생성 로직을 누락하였다.
1. 이 모듈이 도커 컨테이너 내부에서 정상 실행되어 결과 h5ad와 통계 CSV를 성공적으로 산출했음에도 불구하고, 품질 감사·적격성 평가에서 "소프트웨어 버전 추적성(provenance) 미비"로 지적받을 수 있는 구체적 메커니즘을 설명하시오(규정에 `versions.yml`이 명시된 것은 아니며, 재현성·추적성 요건을 지원하는 품질 시스템 관점의 위험이다).
2. nf-core 표준에 부합하도록 다중 종속성(Python 패키지, Samtools 바이너리) 버전을 정확히 파싱하여 단일 다중 레벨 YAML 포맷으로 방출하는 Nextflow bash heredoc 스크립트를 작성하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 규제 감사 지적 메커니즘**:
* 임상 진단 및 동반진단(CDx) 파이프라인은 결과물(h5ad/CSV)뿐만 아니라 **결과물을 생성한 정확한 알고리즘 및 도구 바이너리의 고유 버전(Software Provenance Manifest)**을 입증할 수 있어야 한다.
* 컨테이너 태그(예: `scanpy:latest` 또는 `v1.0`)는 빌드 시점이나 레지스트리 재업로드 과정에서 내부 라이브러리(NumPy, SciPy, Pysam 등)의 서브 버전이 변경될 가능성이 상존한다.
* `versions.yml`이 누락되면 상위 파이프라인의 `MULTIQC` 및 `WORKFLOW_SUMMARY`가 최종 종합 감사 보고서(`software_versions.yml`)를 집계할 때 해당 커스텀 단계의 의존성 계통을 누락하게 된다. 결과적으로 2년 후 환자 생존 분석 재검증 요청 시 동일한 계산 환경을 100% 동일한 비트 단위로 복원할 수 없게 되어 재현성 보증 결함으로 최종 승인이 기각된다.

**2. 표준 `versions.yml` 방출 스크립트 구현**:
```bash
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1 | sed -n '1s/^samtools //p'))
        python: \$(python -c "import sys; print(sys.version.split()[0])")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
```
* `$task.process`를 최상위 키로 매핑하고, 정규식 또는 파이썬 런타임 검사를 통해 각 도구의 순수 버전 문자열만을 추출함으로써 다중 모듈 병합 시 YAML 구문 파싱 에러를 원천 방지한다.

</details>

---

### 문제 3: 골전이 전립선암 생검 조직의 과도한 미토콘드리아 RNA 비율에 대한 적응형 QC 필터링 설계
골전이(Bone Metastasis)를 겪고 있는 전립선암 환자의 골수 생검 코어 조직을 단일세포 전사체 시퀀싱하였다(효소 해리된 전체 세포 scRNA-seq 가정; snRNA-seq에서는 해당하지 않음). 조직학적 특성상 단단한 골기질(Bone Matrix)을 분해하기 위한 강력한 콜라게나아제 처리로 인해 세포막 손상이 발생하였다. FastQC 및 1차 메트릭 분석 결과, 검출된 바코드의 평균 미토콘드리아 전사체 비율($f_{\text{MT}}$)이 $18.5\%$에 달했다.
1. 통상적인 단일세포 분석 기준인 고정 임계값 $f_{\text{MT}} \le 10\%$를 기계적으로 적용할 경우, 이 코호트에서 어떤 임상적 왜곡(Clinical Bias) 및 세포 군집 결손이 발생하는지 설명하시오.
2. 전립선암 상피세포(Adenocarcinoma)와 침윤 골수 유래 면역세포(CD8+ T세포, 대식세포)의 고유한 세포 크기 및 미토콘드리아 함량을 고려하여, 단일 컷오프의 한계를 극복하는 이중 모달리티(Bimodal) 적응형 MAD 기반 필터링 수식을 제안하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 고정 컷오프(10%) 적용 시 발생하는 임상적 왜곡**:
* 표본 전체의 평균 미토콘드리아 비율이 $18.5\%$인 고난도 골전이 검체에 일괄적으로 $f_{\text{MT}} \le 10\%$를 적용하면 전체 세포의 상당 부분(이 가상 시나리오에서는 예시로 70~80%)이 위양성 사멸 세포(False-positive Dead Cells)로 판정되어 폐기될 수 있다.
* 특히 암세포 중에서도 높은 대사 활성과 산화적 인산화(OXPHOS)를 영위하는 악성 공격성 종양 클론(AR-V7 고발현 클론)과 호흡 폭발(Respiratory Burst) 활성이 높은 종양 관련 대식세포(TAMs)가 우선적으로 제거된다. 결과적으로 치료 저항성을 유발하는 핵심 표적 군집이 전사체 아틀라스에서 증발하는 치명적인 선별 편향(Selection Bias)이 유발된다.

**2. 적응형 군집별 MAD 필터링 수식 설계**:
세포 크기(총 UMI 수 $N_{\text{total}}^{(i)}$)와 미토콘드리아 분율 $f_{\text{MT}}^{(i)}$의 결합 분포를 평가하여 1차 조대 필터링(Coarse Filtering) 후 비지도 군집화 기반의 적응형 임계값을 구축한다:
1. **1차 조대 이상치 제거**:
   진정한 파열 세포(Debris)만을 제거하기 위한 완화된 절대 상한선 적용:
   $$\text{Filter}_1: N_{\text{total}}^{(i)} \ge 400 \quad \land \quad f_{\text{MT}}^{(i)} \le 35\%$$
2. **세포 계통 군집별 적응형 MAD 임계값**:
   조대 필터링을 통과한 세포들을 저차원 임베딩하여 상피성 종양 클러스터($C_{\text{Tumor}}$)와 면역 클러스터($C_{\text{Immune}}$)로 1차 분할한 뒤, 각 군집 $k \in \{C_{\text{Tumor}}, C_{\text{Immune}}\}$별로 별도의 컷오프를 동적 계산한다:
   $$\text{Cutoff}_{\text{MT}}^{(k)} = \operatorname{median}_{i \in k}(f_{\text{MT}}^{(i)}) + 3 \times \operatorname{MAD}_{i \in k}(f_{\text{MT}}^{(i)})$$
* 예시로, 면역세포 군집은 대사량이 낮아 약 $12\%$, 상피 종양 세포군은 약 $24\%$의 컷오프가 산출될 수 있다(가상 수치). 군집별 MAD 컷오프는 생물학적 신호와 사멸 노이즈를 분리하는 데 도움이 되지만 완벽한 분리를 보장하지는 않으며 군집 정의 자체가 순환적일 수 있다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- Nextflow 24.04.4 (nextflow.io). 원 논문: Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35(4):316-319. PMID: 28398311. doi: 10.1038/nbt.3820
- Docker (github.com/moby/moby).
- Apptainer (Singularity) (github.com/apptainer/apptainer). 원 논문: Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. *PLoS One*. 2017;12(5):e0177459. PMID: 28494014. doi: 10.1371/journal.pone.0177459
- SAMtools (github.com/samtools/samtools). 원 논문: Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *Gigascience*. 2021;10(2). PMID: 33590861. doi: 10.1093/gigascience/giab008
- FastQC (github.com/s-andrews/FastQC).
- MultiQC (github.com/MultiQC/MultiQC). 원 논문: Ewels P, Magnusson M, Lundin S, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016;32(19):3047-8. PMID: 27312411. doi: 10.1093/bioinformatics/btw354
- nf-core/scrnaseq (nf-co.re/scrnaseq). 원 논문: Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol*. 2020;38(3):276-278. PMID: 32055031. doi: 10.1038/s41587-020-0439-x
- conda (docs.conda.io).
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- alevin-fry (github.com/COMBINE-lab/alevin-fry). 원 논문: He D, Zakeri M, Sarkar H, et al. Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data. *Nat Methods*. 2022;19(3):316-322. PMID: 35277707. doi: 10.1038/s41592-022-01408-3
<!-- REFERENCES:END -->

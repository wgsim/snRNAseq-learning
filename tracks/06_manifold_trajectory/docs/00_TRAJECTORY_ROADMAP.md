# 00. 단일세포 매니폴드 궤적 분석 및 동역학 로드맵 (Trajectory Roadmap)

---

## 학습 목표 (Cognitive Goals)

1. **와딩턴 후성유전학적 풍경(Waddington's Epigenetic Landscape)과 '시간의 화살' 패러독스 규명**: 파괴적 측정(Destructive Measurement)이라는 단일세포 시퀀싱의 물리적 한계 하에서 정적 스냅샷(Static Snapshot)들로부터 시간의 연속성과 비가역적 방향성(Arrow of Time)을 복원하는 계산생물학적 정식화 원리를 이해한다.
2. **정적 그래프 유사시간(Static Graph Pseudotime)의 근본적 결함 분석**: 유클리드/그래프 최단 거리 기반 알고리즘(Monocle, DPT)이 가역적 상태 전이, 탈분화(Dedifferentiation), 교차분화(Transdifferentiation) 및 암 세포 계통 가소성(Lineage Plasticity)을 왜곡하는 수학적 기전을 증명한다.
3. **단일세포 동역학 4대 계산 패러다임의 유기적 연계성 정립**:
   - 매니폴드 확산 기하학(Diffusion Geometry & DPT)
   - 전사체 스플라이싱 반응속도론(RNA Velocity Splicing Kinetics)
   - 완전 동역학 EM 최적화 및 보편 잠재 시간(Dynamical scVelo & Latent Time)
   - 엔트로피 정규화 최적 수송 모델(Waddington-OT)
   위 방법론들의 수학적 가정, 적용 한계 및 상호보완적 역할을 체계화한다.
4. **전립선암 계통 가소성 벤치마크 및 분석 환경 구축**: Karthaus et al. (Science 2020)의 안드로겐 결핍/재생 시계열 모델 및 Replogle et al. (Cell 2022)의 Perturb-seq 섭동 모델을 통한 궤적 분석 실습을 준비하고, 전용 파이썬 환경(`scvelo`, `cellrank`, `wot`, `scanpy`)을 구성한다.

---

## 1. 와딩턴 후성유전학적 풍경과 '시간의 화살' (Arrow of Time)

### 1.1 파괴적 측정(Destructive Measurement)과 단일세포 스냅샷의 한계

단일세포 RNA 시퀀싱(scRNA-seq/snRNA-seq)은 세포를 화학적으로 용해(Lysis)하여 전사체를 포획하므로, **동일한 세포를 시간에 따라 연속적으로 추적하는 것이 물리적으로 불가능**하다. 우리가 획득하는 데이터 행렬 $X \in \mathbb{R}^{N \times G}$는 특정 시점에 정지된 개별 세포들의 정적 스냅샷(Static Snapshot)들의 집합에 불과하다.

```
[The Destructive Measurement Paradox in Single-Cell Biology]

True Continuous Biological Time (t):
t_0 ──────────────────> t_1 ──────────────────> t_2 ──────────────────> t_final
[Cell A] ──(Differentiates)──> [Cell A'] ──(Transdifferentiates)──> [Cell A'']

Observed High-Throughput scRNA-seq:
t_0: {Cell 1, Cell 2, ...} (Cells Destroyed by Lysis!)
t_1: {Cell 101, Cell 102, ...} (Completely Different Individual Cells!)
t_2: {Cell 201, Cell 202, ...} (No Direct Ancestor-Descendant Tracking!)
```

이러한 물리적 제약은 통계역학과 동역학계 이론(Dynamical Systems Theory)에서 다루는 **시간의 화살(Arrow of Time)** 문제를 야기한다:
- 불연속적인 세포 집단 사이에서 어느 세포가 '조상(Ancestor)'이고 어느 세포가 '후손(Descendant)'인가?
- 세포의 상태 전이는 엔트로피가 증가하는 방향으로만 진행하는 비가역적 과정인가, 아니면 병리적 미세환경에서 역류(Reverse Flow)할 수 있는가?

### 1.2 콘래드 와딩턴의 고전적 풍경(Waddington's Landscape)과 현대적 재해석

1957년 발생생물학자 콘래드 와딩턴(Conrad Waddington)은 세포 분화 과정을 높은 고도에서 깊은 골짜기로 굴러떨어지는 구슬에 비유하였다.

```
[Classical Waddington Landscape]         [Modern Cancer Plasticity Landscape]

         Stem/Progenitor                             Stem/Progenitor
             (High PE)                                  (Unstable)
               ○                                          ○   ▲
              / \                                        / \  │ (Dedifferentiation:
             /   \                                      /   \ │  e.g., Castration,
            /     \                                    /     \│  TP53/RB1 Loss)
           /   │   \                                  /       \
          ▼    │    ▼                                ▼         ▼
      Branch 1 │ Branch 2                        Basal/     Luminal
       (Valley)│  (Valley)                      Neuroendo   (Reversible)
               ▼                                     ▲─────────┘
        Terminally Differentiated                    (Transdifferentiation)
           (Global Minimum)
```

1. **고전적 발생학 모델(Classical Development)**:
   - 에너지 지형(Potential Energy Surface $U(x)$)은 정상 세포에서 엄격한 단조 감소(Monotonically Decreasing) 함수로 취급된다:
     $$\frac{dx}{dt} = -\nabla U(x) + \xi(t)$$
   - 여기서 $\xi(t)$는 백색 소음(White Noise) 형태의 전사 파열 요동이며, 구슬은 국소 최소점(Local Minima, 분화 종착점)에 도달하면 외부의 극단적인 자극 없이는 계곡을 거슬러 올라갈 수 없다.
2. **현대 종양생물학 및 계통 가소성(Cancer Lineage Plasticity)**:
   - 유전적 결함(e.g., $TP53, RB1$ 결손) 또는 치료적 스트레스(e.g., 항호르몬 치료, 면역항암제) 하에서 세포는 골짜기 사이의 능선(Barrier)을 넘는 **교차분화(Transdifferentiation)**를 일으키거나, 계곡 바닥에서 다시 상류로 거슬러 올라가는 **탈분화(Dedifferentiation)**를 감행한다.
   - 따라서 정적인 포텐셜 경사만으로는 암 조직 내의 비평형(Non-equilibrium) 동역학을 포착할 수 없다.

---

## 2. 정적 그래프 유사시간(Static Graph Pseudotime)의 붕괴

초기 궤적 추정 알고리즘들(Wanderlust, Monocle 1/2, Diffusion Pseudotime)은 세포 간의 발현 프로파일 유사도를 바탕으로 $k$-최근접 이웃 그래프(k-NN Graph)를 구축한 뒤, 그래프 상의 최단 경로(Shortest Geodesic Path)를 유사시간(Pseudotime, $t_{\text{pseudo}}$)으로 환산하였다.

### 2.1 정적 알고리즘의 네 가지 구조적 한계

```
[Pitfalls of Distance-Based Static Pseudotime]

1. Boundary Condition Artifact:
   Selected Root (User) ───[Shortest Path Distance]───> Artifact Endpoints
   (Misidentification of root inverts biological time entirely!)

2. Irreversibility Assumption:
   A (Luminal) ───────────────> B (Basal)  [Assumed One-way]
   Real Biology: A <──────────> B (Plasticity & Oscillation)

3. Density & Sampling Bias:
   High Cell Density Region ===> Graph distance artificially compressed!
   Low Cell Density Region  ===> Bridge edges span spurious transitions!

4. Rate Distortion (Pseudotime ≠ Physical Time):
   Linear Pseudotime Δτ = 1 does NOT mean Δt = 1 hour.
   Cells can arrest in quiescent states for weeks without gene changes.
```

1. **방향성의 부재(Lack of Intrinsic Directionality)**:
   그래프 상의 최단 경로는 대칭적($d(u, v) = d(v, u)$)이다. 분석가가 사전에 시작 세포(Root cell)를 수동으로 지정하지 않으면, 시스템은 탈분화인지 정상 분화인지 자체적으로 판단할 수 없다.
2. **비가역성 공리의 오류(Violation of Irreversibility)**:
   정적 그래프는 순환 고리(Cycle)나 궤적 역행을 다룰 수 없으며, 강제적으로 최소 신장 트리(Minimum Spanning Tree, MST)나 단방향 유향 비순환 그래프(DAG)로 매니폴드를 투영한다.
3. **세포 밀도 불균일성에 따른 거리 왜곡**:
   세포 증식 속도와 생존율 차이로 인해 중간 상태의 세포 밀도가 극히 낮을 경우, 그래프는 연결성이 끊어지거나 잘못된 지름길(Short-circuit Bridge)을 생성한다.
4. **물리적 시간과의 괴리**:
   유사시간은 유전자 발현 공간의 '거리'일 뿐이며, 세포 주기 정지(Quiescence)나 휴지기 동안의 물리적 절대 시간($t$)의 흐름을 반영하지 못한다.

---

## 3. 단일세포 동역학 4대 계산 패러다임

단일세포 연구는 정적 그래프 거리를 탈피하여, 세포 내부의 물리적/통계적 원리를 이용해 내재적 방향성 벡터를 복원하는 동역학 패러다임으로 진화하였다.

```
[Evolution of Single-Cell Dynamical Trajectory Inference]

+───────────────────────────────────────────────────────────────────────────+
| 1. Manifold Diffusion Geometry (Haghverdi et al., Nature Methods 2016)   |
|    Graph Laplacian & Diffusion Random Walk:                               |
|    D^2(x, y) = Σ [T^t(x, ·) - T^t(y, ·)]^2                               |
|    - Random-walk transition noise filtering, multi-scale geodesics        |
+───────────────────────────────────────────────────────────────────────────+
                                    │
                                    ▼
+───────────────────────────────────────────────────────────────────────────+
| 2. RNA Velocity Splicing Kinetics (La Manno et al., Nature 2018)          |
|    Biophysical pre-mRNA ODEs:                                             |
|    du/dt = α - β*u,   ds/dt = β*u - γ*s                                   |
|    - First intrinsic "Arrow of Time" using splicing displacement          |
+───────────────────────────────────────────────────────────────────────────+
                                    │
                                    ▼
+───────────────────────────────────────────────────────────────────────────+
| 3. Dynamical scVelo & Latent Time (Bergen et al., Nature Biotech 2020)    |
|    Full Kinetic EM Optimization:                                          |
|    Optimizes (α, β, γ, t) without steady-state equilibrium assumption     |
|    - Reconstructs universal latent time (τ) resolving transient regimes   |
+───────────────────────────────────────────────────────────────────────────+
                                    │
                                    ▼
+───────────────────────────────────────────────────────────────────────────+
| 4. Entropic Optimal Transport: Waddington-OT (Schiebinger et al., 2019)  |
|    Coupling discrete longitudinal populations:                            |
|    min <C, Γ> - ε H(Γ), with proliferation/death growth rates g(x)        |
|    - Rigorous mathematical push-forward across experimental timepoints    |
+───────────────────────────────────────────────────────────────────────────+
```

### 3.1 방법론별 수리물리학적 핵심 비교

| 패러다임 | 핵심 수리 도구 | 방향성 정보원 | 주요 강점 | 적용 시 주의점 |
| :--- | :--- | :--- | :--- | :--- |
| **Diffusion Pseudotime (DPT)** | 그래프 라플라시안 스펙트럼 분해, 마르코프 무작위 보행 | 외부 지정 (Root Cell 필요) | 노이즈 강건성, 다분기(Multibranching) 기하 구조 보존 | 역행/순환 궤적 불가능, 분기 방향 지정 오류 위험 |
| **RNA Velocity (Steady-state)** | 1차 선형 상미분방정식(ODE) 극단 위상 회귀 | 스플라이스되지 않은 pre-mRNA 잔차 ($u - \tilde{\gamma}s$) | 단일 시점 스냅샷에서 세포별 순간 속도 벡터 도출 | 전사 정상상태(Equilibrium) 가정 붕괴 시 왜곡 |
| **Dynamical scVelo** | 완전 동역학 우도 기반 EM 알고리즘 | 유전자별 위상 궤도 최적화 및 보편 잠재 시간 ($\tau$) | 비평형 전사 유도/억제 구간 모델링, 다중 반응속도 복원 | 계산 비용 큼, 스플라이싱이 적은 저발현 유전자 피팅 불안정 |
| **Waddington-OT** | 엔트로피 정규화 최적 수송, 싱크혼-크놉(Sinkhorn-Knopp) 알고리즘 | 실제 실험적 시계열 ($t_1 < t_2 < \dots < t_k$) | 세포 증식/사멸($g(x)$) 통합, 대규모 섭동 계통 추적 | 시계열 데이터 필수, 극단적 배치 효과 존재 시 수송 왜곡 |

---

## 4. 종양 생물학 벤치마크: 전립선암 계통 가소성

본 매니폴드 궤적 분석 트랙의 모든 방법론은 전립선암(Prostate Cancer)의 대표적 난제인 **계통 가소성(Lineage Plasticity) 및 거세 저항성 신경내분비 분화(t-NEPC)** 벤치마크 데이터를 관통하여 검증된다.

### 4.1 Karthaus et al. (Science 2020) 거세/재생 시계열 벤치마크
- **실험 디자인(교육용 요약)**: 마우스 전립선에서 거세(Castration)와 재생(Regeneration) 과정을 단일세포 전사체로 분석한 연구입니다(Karthaus et al., Science 2020). 아래 날짜별 샘플링 일정(Day 0/1/3/7/14/28)은 원 논문에서 확인한 것이 아니라 이 트랙의 *교육용 예시 일정*입니다(본 트랙 모의 데이터는 Day 0/1/3/7/14).
- **핵심 생물학적 발견**:
  - 논문 초록은 마우스와 인간의 관강 세포(Luminal cells)가 재생에 기여하며, 희귀한 줄기세포 유사 관강 세포와 분화된 관강 세포가 모두 재생에 기여함을 보고합니다(PMID 32355025). 즉 '분화 관강세포도 재생에 기여' 수준이 원 논문에서 확인된 내용입니다.
  - 이 시계열에 RNA 벨로시티나 Waddington-OT를 적용해 '탈분화 벡터'를 증명했다는 서술은 원 논문의 주장이 아니라, 이 트랙에서 방법론을 설명하기 위한 *가상의 적용 시나리오*입니다. 정적 유사시간이 방향성을 모호하게 만든다는 점은 일반론입니다.

### 4.2 Replogle et al. (Cell 2022) 게놈 규모 Perturb-seq 벤치마크
- **실험 디자인**: 인간 세포주를 대상으로 genome-scale CRISPRi Perturb-seq(250만 세포 이상)를 수행한 데이터셋입니다(Replogle et al., Cell 2022, PMID 35688146). 전립선 조직이나 신경내분비 전이를 분석한 연구가 아닙니다.
- **본 트랙에서의 위치**: 섭동 데이터에서 상태 이동을 정량하는 *방법론 사례*로만 인용합니다. 아래 '전사인자 결손(AR, FOXA1, TP53, RB1)이 선암 매니폴드에서 신경내분비 상태로 도약하는 전이 확률' 분석은 원 논문의 결과가 아니라 *가상 시나리오*입니다.

---

## 5. 계산 환경 구축 (Python Setup)

궤적 분석 및 동역학 추론을 위해 고성능 파이썬 단일세포 라이브러리 스택을 구성한다.

```bash
#!/usr/bin/env bash
# ==============================================================================
# Single-Cell Manifold Trajectory & Dynamics Environment Setup Script
# ==============================================================================
set -euo pipefail

ENV_NAME="sc_trajectory"
PYTHON_VERSION="3.10"

echo "[1/4] Conda/Mamba 환경 생성: ${ENV_NAME} (Python ${PYTHON_VERSION})"
conda create -y -n "${ENV_NAME}" python="${PYTHON_VERSION}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

echo "[2/4] 핵심 바이오인포매틱스 및 그래프 라이브러리 설치"
conda install -y -c conda-forge -c bioconda \
    scanpy \
    python-igraph \
    leidenalg \
    louvain \
    scikit-learn \
    scipy \
    pandas \
    matplotlib \
    seaborn

echo "[3/4] RNA 벨로시티 및 동역학 전용 패키지 설치"
pip install --no-cache-dir \
    scvelo>=0.3.1 \
    cellrank>=2.0.0 \
    wot>=1.0.8 \
    anndata>=0.9.0

echo "[4/4] 환경 검증 및 버전 출력"
python -c "
import scanpy as sc
import scvelo as scv
import cellrank as cr
import wot

print(f'Scanpy Version:   {sc.__version__}')
print(f'scVelo Version:   {scv.__version__}')
print(f'CellRank Version: {cr.__version__}')
print('Trajectory analysis environment configured successfully!')
"
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: 시간의 화살과 파괴적 측정(Destructive Measurement)의 수학적 정식화
단일세포 RNA 시퀀싱 데이터셋에서 $K$개의 실험 시점 $t_1 < t_2 < \dots < t_K$에서 세포 군집 $P(t_k)$를 각각 샘플링하였다.
1. 왜 임의의 세포 $i \in P(t_k)$와 세포 $j \in P(t_{k+1})$ 사이의 유클리드 거리 $\|x_i - x_j\|_2$를 최소화하는 최근접 이웃 연결이 실제 조상-후손 계통(Lineage) 관계를 보장하지 못하는지 고차원 측도 집중(Concentration of Measure) 및 세포 사멸/증식 관점에서 수식적으로 논증하시오.
2. Waddington-OT가 이 문제를 해결하기 위해 도입한 통계역학적 결합 확률 분포(Coupling Distribution $\Gamma$)의 수학적 정의를 기술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 유클리드 최근접 이웃 연결의 계통 보장 실패 기전**:
- **고차원 측도 집중 현상 (Concentration of Measure)**:
  유전자 발현 공간의 차원 $G \sim 20,000$ (고가변 유전자 사용 시에도 $G \sim 2,000$)에서, 무작위로 분포하는 점들 간의 거리 비율은 $\lim_{G \to \infty} \frac{D_{\max} - D_{\min}}{D_{\min}} \to 0$으로 수렴한다. 따라서 고차원 유클리드 거리는 생물학적 조상 관계에 의한 유사성과 단순한 기술적 드롭아웃/시퀀싱 노이즈에 의한 유사성을 통계적으로 유의미하게 구분하지 못한다.
- **세포 증식 및 사멸에 따른 질량 비보존(Non-conservation of Mass)**:
  유클리드 최근접 이웃(1-NN)은 $t_k$의 모든 세포가 $t_{k+1}$로 반드시 일대일 생존 전이한다고 암묵적으로 가정한다. 그러나 실제 생체 조직에서는 세포 자멸사(Apoptosis)로 인해 후손 없이 소멸하는 세포가 존재하며, 특정 클론의 폭발적 증식(Clonal Expansion)으로 인해 하나의 조상 세포가 수십 개의 후손 세포를 생성한다. 따라서 거리 기반의 결정론적 매칭은 사멸하는 세포에 엉뚱한 후손을 강제 할당하는 치명적 위양성(False Positive) 궤적을 생성한다.

**2. Waddington-OT의 결합 확률 분포 $\Gamma$ 정식화**:
Waddington-OT는 개별 세포 간의 점대점 매칭 대신, 두 시점의 확률 측도 $\mu \in \mathcal{P}(X_{t_k})$와 $\nu \in \mathcal{P}(X_{t_{k+1}})$ 사이의 결합 측도(Joint Coupling Distribution) $\Gamma \in \mathbb{R}_{+}^{N_1 \times N_2}$를 최적화한다:
$$\Gamma_{ij} = P\left(\text{Cell } i \text{ at } t_k \text{ transitions to Cell } j \text{ at } t_{k+1}\right)$$
이때 질량 보존 법칙은 세포의 분열/사멸 증식률 $g(x_i)$를 반영한 주변부 제약 조건(Marginal Constraints)으로 일반화된다:
$$\sum_{j=1}^{N_2} \Gamma_{ij} = g(x_i) \mu_i, \quad \sum_{i=1}^{N_1} \Gamma_{ij} = \nu_j$$
이를 통해 증식하는 세포와 사멸하는 세포의 질량 변화를 엄밀하게 보정하며 거시적 전이 확률을 복원한다.

</details>

---

### 문제 2: 정적 그래프 유사시간의 탈분화(Dedifferentiation) 판별 실패
전립선암 환자 유래 조직에서 안드로겐 수용체 양성($AR^+$) 관강 상피세포(Luminal)가 호르몬 결핍 치료 후 줄기세포 유사 상태(Stem-like) 및 기저 상피세포(Basal) 표현형을 획득하는 탈분화 과정이 발생하였다.
1. Monocle이나 일반 DPT와 같은 정적 그래프 알고리즘을 적용할 때, 분석가가 발생학적 상식에 따라 '줄기세포 유사 세포'를 시작 세포(Root)로 지정할 경우 도출되는 치명적인 생물학적 결론 왜곡을 설명하시오.
2. 이러한 왜곡을 극복하기 위해 왜 RNA 스플라이싱 반응속도론(RNA Velocity) 또는 종단적 최적 수송(OT)이 필수적인지 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 시작 세포 수동 지정에 따른 생물학적 결론 왜곡**:
- 정적 그래프 알고리즘은 대칭적 거리 행렬을 기반으로 하므로, 방향성의 원천이 알고리즘 내부에 존재하지 않고 전적으로 분석가가 지정하는 시작점(Root Cell)에 종속된다.
- 일반적인 발생학적 편향(Developmental Bias)에 따라 줄기세포 유사 세포(Stem-like)를 Root로 지정하면, 알고리즘은 자동으로:
  $$\text{Stem-like} \longrightarrow \text{Luminal}$$
  방향으로 유사시간 화살표를 강제 부여한다.
- 그 결과, 실제로는 호르몬 치료 저항성 획득을 위해 성숙 관강 세포가 탈분화($\text{Luminal} \to \text{Stem-like}$)하는 병리적 역행 프로세스가 발생했음에도 불구하고, 분석 결과는 "줄기세포가 정상적으로 관강 세포로 분화하고 있으며 치료에 잘 반응하고 있다"는 정반대의 치명적인 임상적 오판을 초래하게 된다.

**2. 동역학적 방법론(RNA Velocity / OT)의 필수성**:
- **RNA Velocity**: 세포 내부에서 새로 합성 중인 미성숙 pre-mRNA($u$)와 성숙 mRNA($s$)의 비율 편차를 측정하므로, 외부의 Root 지정 없이도 세포 내부의 유전자 전사 가속도($v = \frac{ds}{dt}$)를 통해 세포가 현재 어느 상태를 향해 이동하고 있는지 고유한 시간의 화살을 직접 계산한다. 따라서 $\text{Luminal} \to \text{Stem-like}$로 향하는 역류 벡터를 자율적으로 검출한다.
- **Waddington-OT**: 치료 전($t_0$)과 치료 후($t_1$)의 실제 물리적 시간 정보를 이용하여, $t_0$의 관강 세포 집단이 $t_1$의 줄기세포 유사 집단으로 수송되는 확률 질량($\Gamma_{\text{Luminal} \to \text{Stem-like}} > 0$)을 객관적으로 입증한다.

</details>

---

### 문제 3: 암 계통 가소성 추론에서 세포 주기(Cell Cycle)와 정지 상태(Quiescence)의 시간 왜곡
Karthaus et al. 거세 시계열 데이터에서, 거세 직후 잔존하는 관강 상피세포(Castration-tolerant Luminal Cells)는 장기간 세포 분열을 멈추고 $G_0$ 정지 상태(Quiescence)에 머무른다.
1. 그래프 상의 유사시간 거리 $\Delta \tau$와 실제 생체 물리적 시간 $\Delta t$ 사이의 불일치가 발생하는 근본 기전을 설명하시오.
2. 세포 주기 유전자 발현 변동이 매니폴드 기하학에 미치는 순환적 교란(Periodic Confounding)을 제어하기 위한 계산적 전처리 기법을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. $\Delta \tau$와 $\Delta t$ 불일치의 근본 기전**:
- **유사시간($\tau$)의 정의**: 전사체 매니폴드 상에서 세포 간의 누적된 전사체 변화량(Transcriptional Distance)의 적분값이다:
  $$\tau(x) \propto \int_{\gamma} \|\nabla X\| ds$$
- **정지 상태에서의 시간 왜곡**:
  세포가 $G_0$ 휴지기(Quiescence)에 진입하면, 물리적 시간 $t$는 수주~수개월 동안 지속적으로 경과($\Delta t \gg 0$)함에도 불구하고 대사 활동 저하 및 전사 억제로 인해 전사체 프로파일의 변화는 거의 0에 수렴한다($\Delta X \approx 0$).
- 따라서 정적 유사시간 상에서는 이 세포들이 단 몇 분 만에 통과한 극히 짧은 궤적으로 압축 표현되어, 휴지기 세포의 실제 생체 내 생존 시간과 치료 저항성 기여도가 과소평가된다.

**2. 세포 주기 순환 교란 제어 기법**:
- **세포 주기 점수화(Cell Cycle Scoring)**:
  Tirosh et al.의 $S$기 및 $G_2/M$기 특이적 마커 유전자 세트(e.g., $MKI67, TOP2A, PCNA, MCM6$)를 이용하여 각 세포의 $S$-Score 및 $G_2/M$-Score를 산출한다.
- **선형 회귀 제거(Regression out of Cell Cycle)**:
  스케일링 전 단계에서 일반 선형 모델(GLM)을 통해 세포 주기 점수를 독립변수로 설정하고, 이에 기인한 발현량 변동 잔차를 제거한다:
  $$x_{cg} = \beta_0 + \beta_1 S_c + \beta_2 G2M_c + \epsilon_{cg} \implies \tilde{x}_{cg} = \epsilon_{cg}$$
- **분화 특이적 분산 보존 주의점**:
  단, 줄기세포-전구세포 증식 과정 자체가 계통 분화와 불가분하게 얽혀 있는 경우, 무차별적인 세포 주기 제거는 분화 궤적 신호 자체를 파괴할 수 있다. 따라서 분화와 관련된 전사인자 네트워크와 순수 유사분열(Mitotic) 사이클 유전자를 분리하여 선택적으로 교란 요인을 제거해야 한다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Karthaus WR, Hofree M, Choi D, et al. Regenerative potential of prostate luminal cells revealed by single-cell analysis. *Science*. 2020;368(6490):497-505. PMID: 32355025. doi: 10.1126/science.aay0267
2. Replogle JM, Saunders RA, Pogson AN, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*. 2022;185(14):2559-2575.e28. PMID: 35688146. doi: 10.1016/j.cell.2022.05.013
3. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*. 2016;352(6282):189-96. PMID: 27124452. doi: 10.1126/science.aad0501
4. Haghverdi L, Büttner M, Wolf FA, et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods*. 2016;13(10):845-8. PMID: 27571553. doi: 10.1038/nmeth.3971
5. La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells. *Nature*. 2018;560(7719):494-498. PMID: 30089906. doi: 10.1038/s41586-018-0414-6
6. Bergen V, Lange M, Peidli S, et al. Generalizing RNA velocity to transient cell states through dynamical modeling. *Nat Biotechnol*. 2020;38(12):1408-1414. PMID: 32747759. doi: 10.1038/s41587-020-0591-3
7. Schiebinger G, Shu J, Tabaka M, et al. Optimal-Transport Analysis of Single-Cell Gene Expression Identifies Developmental Trajectories in Reprogramming. *Cell*. 2019;176(4):928-943.e22. PMID: 30712874. doi: 10.1016/j.cell.2019.01.006
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- conda (docs.conda.io).
<!-- REFERENCES:END -->

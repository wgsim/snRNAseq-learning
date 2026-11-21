# Chapter 03. 발현 임계치 설정과 질량작용 법칙 기반 순열 검정

> **Cognitive Goal**: 단일세포 전사체 시퀀싱(scRNA-seq / snRNA-seq)의 기술적 한계인 제로 팽창(Zero-Inflation) 및 드롭아웃(Dropout) 현상이 리간드-수용체 상호작용 추론에 미치는 통계적 왜곡을 이해하고, 발현 백분율 임계치($\text{min.pct}$) 필터링의 생물학적 타당성, 질량작용 법칙(Law of Mass Action)에 기반한 상호작용 점수 모델링, 세포 라벨 무작위 순열 검정(Permutation Test)과 벤자미니-호흐바르크(Benjamini-Hochberg FDR) 다중 가설 보정의 수학적 원리를 유도하고 구현할 수 있다.

---

## 1. 단일세포 시퀀싱의 제로 팽창(Zero-Inflation)과 드롭아웃(Dropout) 도전 과제

단일세포 RNA 시퀀싱 데이터의 발현 행렬(Count Matrix)은 전체 엔트리의 $80\%\sim95\%$ 이상이 $0$으로 채워진 극도의 희소 행렬(Sparse Matrix)입니다. 이러한 제로 값은 생물학적 원인과 기술적 결함이 혼재되어 발생합니다:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        단일세포 데이터의 제로(0) 발생 이중 기전                           │
├────────────────────────────────────────────────────────────────────────────────────────┤
│ 1. 생물학적 제로 (Biological Zero) : 해당 유전자가 실제로 전사되지 않아 mRNA가 없음     │
│ 2. 기술적 드롭아웃 (Technical Dropout) : 세포 내에 mRNA 분자가 소량 존재했으나          │
│    - 역전사(Reverse Transcription) 효소 효율 한계 (10% ~ 20% 포획률)                   │
│    - 마이크로플루이딕 캡슐화 및 PCR 증폭 손실                                           │
│    - 얕은 시퀀싱 깊이(Shallow Sequencing Depth, 세포당 20,000~50,000 reads)             │
│    로 인해 cDNA 라이브러리에 포획되지 못하고 0으로 관측되는 측정 오류                   │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 1.1 신호 전달 분자의 발현 특수성과 검출 한계
세포 간 신호 전달에 관여하는 분자들은 하우스키핑 유전자(액틴, 리보솜 단백질)와 근본적으로 다른 발현 역학을 지닙니다:
- **리간드(사이토카인, 성장인자)**: 정상 상태에서는 거의 전사되지 않다가 특정 자극 시 일시적으로 폭발 발현(Bursting)되거나, 세포당 불과 수 개~수십 개의 mRNA 분자만 유지됩니다.
- **수용체(GPCR, 인산화효소 수용체)**: 증폭 신호 전달 체계의 특성상 세포 표면에 소수만 존재해도 강력한 하류 캐스케이드를 유발할 수 있으므로, mRNA 카운트가 매우 낮습니다.

만약 아무런 전처리나 필터링 없이 원시 카운트 행렬을 분석에 투입할 경우, **기술적 드롭아웃으로 인한 막대한 위음성(False Negative)**과, 반대로 **소수의 이상치 세포(Outlier)에 의한 평균 발현량 왜곡(False Positive)**이 동시에 발생하는 진퇴양난에 빠지게 됩니다.

---

## 2. 발현 백분율 임계치 설정의 생물학적·통계적 당위성

세포 클러스터의 대표 발현량을 산출할 때, 단순 산술평균(Mean Expression)만을 사용하면 심각한 통계적 착시가 발생합니다.

```
[산술평균의 맹점: 단일 이중핵(Doublet) 또는 이상치 세포에 의한 왜곡]
세포 집단 A (100개 세포): 99개 세포의 발현량 = 0, 1개 세포의 발현량 = 100
  ──► 산술평균 E_bar = 1.0 (마치 집단 전체가 유의미하게 발현하는 것처럼 오판)

세포 집단 B (100개 세포): 50개 세포의 발현량 = 2, 50개 세포의 발현량 = 0
  ──► 산술평균 E_bar = 1.0 (실제 절반의 세포가 균일하게 신호를 생산)
```

### 2.1 발현 백분율 ($\text{Pct}_{c, g}$) 정의 및 임계치 필터링
세포 유형 $c$에 속한 전체 세포 집단을 $C_c$, 세포 $i$에서의 유전자 $g$의 발현량을 $E_{i, g}$라 할 때, 세포 클러스터 내 발현 세포 비율(Percentage of expressing cells)은 지시 함수(Indicator function)를 통해 다음과 같이 정의됩니다:

$$\text{Pct}_{c, g} = \frac{1}{|C_c|} \sum_{i \in C_c} \mathbb{I}(E_{i, g} > 0)$$

CellPhoneDB(`threshold` 기본 0.1)와 LIANA(`expr_prop` 기본 0.1) 등은 신호 전달의 생물학적 항상성을 담보하기 위해 **발현 비율 임계치 $\tau$ (일반적으로 $10\%\sim20\%$)**를 적용합니다. (CellChat은 비율 임계 대신 발현 평균에 trimean $(Q_1 + 2\,\text{median} + Q_3)/4$을 쓰며, 20% 발현이면 0, 30% 발현이면 양수로 바뀌어 사실상 비율 임계처럼 작동합니다.)

$$\text{Active}(L, R \mid A, B) = \begin{cases} 
\text{True}, & \text{if } \text{Pct}_{A, L} \ge \tau \quad \text{and} \quad \text{Pct}_{B, R} \ge \tau \\ 
\text{False}, & \text{otherwise} 
\end{cases}$$

- **생물학적 근거**: 한 클러스터 내에서 최소 $10\%\sim20\%$ 이상의 세포가 지속적으로 리간드를 분비하거나 수용체를 발현해야만, 세포 간질액 내에 유효 결합 농도($[L] \ge K_d$)를 형성할 수 있는 생물학적 집단 효과(Quorum/Paracrine pool)가 완성됩니다.
- **통계적 근거**: 소수의 이중핵(Doublet)이나 이상치 세포로 인해 비정상적으로 부풀려진 평균값의 영향을 줄입니다. 다만 앰비언트 RNA는 비율 필터로 제거되지 않으므로 별도 보정이 필요합니다.

---

## 3. 질량작용 법칙 (Law of Mass Action) 기반 상호작용 점수 모델링

세포 간 신호의 세기는 어떻게 수학적으로 정량화할 수 있는가? 가장 널리 사용되는 기본 물리화학적 모델은 **질량작용 법칙(Law of Mass Action)**입니다.

### 3.1 물리화학적 결합 평형과 점수 유도
단일 리간드 $L$과 단일 수용체 $R$이 용액 내에서 결합하여 복합체 $[LR]$을 형성하는 가역 반응을 고려합니다:

$$L + R \underset{k_{\text{off}}}{\overset{k_{\text{on}}}{\rightleftharpoons}} [LR]$$

화학 평형 상태에서 정반응 속도와 역반응 속도는 동일하므로:
$$k_{\text{on}} [L][R] = k_{\text{off}} [LR]$$

해리 상수(Dissociation Constant) $K_d$를 $K_d = \frac{k_{\text{off}}}{k_{\text{on}}}$로 정의하면, 형성된 복합체의 평형 농도 $[LR]$은 다음과 같습니다:
$$[LR] = \frac{[L][R]}{K_d}$$

즉, 생성되는 신호 복합체의 농도(신호 전달 강도)는 **리간드의 유효 농도와 수용체의 유효 농도의 곱에 정비례**합니다.

### 3.2 세포 유형 간 상호작용 점수 ($S_{ij}$)
단일세포 데이터에서 세포 유형 $A$가 분비하는 리간드 $L$의 평균 발현량을 $\bar{E}_{A, L}$, 세포 유형 $B$가 발현하는 수용체 $R$의 평균 발현량을 $\bar{E}_{B, R}$이라 할 때, 세포 쌍 $(A, B)$ 간의 리간드-수용체 상호작용 점수 $S_{AB}^{L-R}$은 다음과 같이 정의됩니다:

$$S_{AB}^{L-R} = \bar{E}_{A, L} \times \bar{E}_{B, R}$$

여기서 각 클러스터의 평균 발현량은 로그 정규화된 발현값(Log-normalized expression)을 원본 스케일로 역변환하여 계산하거나, 정규화된 발현 수준의 산술평균을 취합니다:
$$\bar{E}_{A, L} = \frac{1}{|C_A|} \sum_{i \in C_A} E_{i, L}, \quad \bar{E}_{B, R} = \frac{1}{|C_B|} \sum_{j \in C_B} E_{j, R}$$

이 곱셈형(Multiplicative) 점수는 두 분자 중 어느 하나라도 결손되면($\bar{E} \to 0$) 전체 상호작용 점수가 0으로 수렴합니다.

> **도구별 점수 정의 (소스 확인)**: 위 곱 $\bar{E}_{A,L}\times\bar{E}_{B,R}$은 NATMI/Connectome 계열(`expr_prod`)의 정의이며 질량작용 법칙의 단순 근사입니다. **CellPhoneDB**의 점수는 두 평균의 **산술평균** $\mathrm{mean}(\bar{E}_L,\bar{E}_R)$(한쪽이 0이면 0)이고 다중 소단위체 복합체는 소단위체 평균 발현의 **최솟값**을 씁니다(liana `_cellphonedb.py`; CellChat 논문). **CellChat**은 질량작용 법칙에 Hill 포화(Kh)를 적용한 식을 쓰고 복합체는 기하평균을 씁니다.

---

## 4. 세포 라벨 무작위 순열 검정 (Permutation Testing)

상호작용 점수 $S_{AB}^{L-R}$이 양수라는 사실만으로는 이 신호가 생물학적으로 유의미한지, 아니면 두 유전자가 전체 세포에 걸쳐 비특이적으로 흔하게 발현되어 우연히 높게 측정된 것인지 알 수 없습니다. 이를 검증하기 위해 **비모수적 순열 검정(Non-parametric Permutation Test)**을 수행합니다.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        세포 라벨 순열 검정(Permutation Test) 알고리즘                    │
├────────────────────────────────────────────────────────────────────────────────────────┤
│ [1. 관측값 계산] : 실제 세포 라벨 (A, B) 상태에서 관측 상호작용 점수 S_obs 산출        │
│ [2. 라벨 셔플링] : 모든 세포의 유전자 발현 프로파일은 그대로 보존하고,                  │
│                   세포 유형 라벨(Cluster ID)만을 무작위 순열 pi로 셔플링              │
│ [3. 영분포 생성] : 셔플된 가상 클러스터 A*, B*에 대해 영점수 S_b* 재계산 (B회 반복)      │
│ [4. p-value 계산] : 관측 점수 S_obs 이상을 기록한 순열 횟수의 비율 정량화              │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 4.1 귀무가설 ($H_0$)과 대립가설 ($H_1$)
- **귀무가설 ($H_0$)**: 세포 유형 $A$와 $B$ 사이의 리간드-수용체 공동 발현 강도는, 전체 세포 집단에서 두 세포 집단의 크기만큼 임의로 세포들을 무작위 추출했을 때 기대되는 공동 발현 강도와 통계적으로 차이가 없다 (즉, 세포 유형 특이적 상호작용이 아니다).
- **대립가설 ($H_1$)**: 세포 유형 $A$와 $B$ 사이의 리간드-수용체 공동 발현 강도는 무작위 추출 기댓값보다 유의하게 높다 (즉, 해당 세포 쌍 간에 특이적으로 풍부화된 신호 전달이다).

### 4.2 경험적 유의확률 (Empirical $p$-value)의 수학적 정식화
$B$회의 무작위 순열 반복(Permutation iterations, 예: $B=1,000$ 또는 $B=10,000$)을 수행하여 얻은 영분포 상호작용 점수를 $\{S_1^*, S_2^*, \dots, S_B^*\}$라 할 때, 경험적 $p$-value는 다음과 같이 계산됩니다:

$$p = \frac{\sum_{b=1}^B \mathbb{I}(S_b^* \ge S_{\text{obs}}) + 1}{B + 1}$$

여기서 $\mathbb{I}(\cdot)$는 조건이 참이면 1, 거짓이면 0을 반환하는 지시 함수입니다.

> **도구별 p-value 정의**: 아래 $+1$ 보정은 일반적으로 권장되는 순열 p-value(Phipson & Smyth 2010)이며, **CellPhoneDB**의 p-value는 "순열 점수가 관측 점수 이상인 횟수 / 순열 횟수"로 **+1 보정이 없습니다**(그래서 0이 나올 수 있음).

#### 왜 분자와 분모에 $+1$ (라플라스 유사 유사 카운트)을 더하는가?
- 만약 $S_{\text{obs}}$가 $B=1,000$번의 모든 순열 점수보다 엄격히 커서 $\sum \mathbb{I} = 0$이 될 경우, 단순 나눗셈을 적용하면 $p = \frac{0}{1000} = 0$이 됩니다.
- 그러나 유한한 표본 추출에서 $p=0$은 통계학적으로 정의될 수 없는 오류(불가능성 단정)입니다. 관측된 실제 데이터 자체를 하나의 유효한 순열 표본으로 간주하여 분자와 분모에 $+1$을 부여함으로써:
  1. $p$-value가 항상 엄격하게 양수($p \ge \frac{1}{B+1} > 0$)로 유지됩니다.
  2. 영가설 하에서 $P(p \le \alpha) \le \alpha$를 만족하는 유효한(valid) 순열 $p$-value가 되어 제1종 오류를 통제합니다(Phipson & Smyth 2010; 약간 보수적인 추정이며 '비편향'이라는 표현은 부정확).

---

## 5. 다중 가설 검정 보정 (Multiple Testing Correction: Benjamini-Hochberg FDR)

단일세포 통신 분석에서는 수십 개의 세포 유형 간에 수천 개의 L-R 쌍을 동시 검정하므로, 총 가설 검정의 수 $M$이 수만~수십만 개에 달합니다:
$$M = (\text{세포 유형 수})^2 \times (\text{L-R 쌍 수})$$

만약 유의수준 $\alpha = 0.05$로 다중 검정 보정 없이 개별 가설을 채택한다면, $M=20,000$개 검정 중 우연에 의해 거짓으로 유의 판정되는 위양성(False Positives)의 수가 무려 $20,000 \times 0.05 = 1,000$개에 이르게 됩니다.

### 5.1 본페로니 보정(FWER) vs 벤자미니-호흐바르크(FDR)
- **본페로니 보정 (Bonferroni FWER)**: 적어도 하나의 위양성을 허용하지 않는 극단적 보정 기준 ($\alpha' = \alpha / M$). 검정력(Statistical Power)이 극도로 저하되어 실제 중요한 생물학적 신호를 모조리 위음성으로 누락시킵니다.
- **벤자미니-호흐바르크 (Benjamini-Hochberg False Discovery Rate, BH-FDR)**: 유의하다고 판정한 전체 발견(Discoveries) 중 위양성이 차지하는 비율의 기댓값을 $q$ 이하로 통제합니다.

### 5.2 Benjamini-Hochberg (BH) 알고리즘 단계
1. $M$개의 독립 검정에서 도출된 $p$-value들을 오름차순으로 정렬합니다:
   $$p_{(1)} \le p_{(2)} \le \dots \le p_{(M)}$$
2. 주어진 FDR 허용 임계치 $\alpha$ (예: $0.05$)에 대해 다음 부등식을 만족하는 가장 큰 순위 $k$를 탐색합니다:
   $$k = \max \left\{ i \in \{1, \dots, M\} : p_{(i)} \le \frac{i}{M} \alpha \right\}$$
3. 정렬된 순위 $i \le k$에 해당하는 모든 가설에 대해 귀무가설을 기각하고 유의한 상호작용으로 최종 판정합니다.
4. 개별 상호작용의 보정된 유의확률($q$-value)은 다음과 같이 단조 증가성을 보장하도록 변환됩니다:
   $$q_{(i)} = \min_{j \ge i} \left( \min \left( 1, \frac{M}{j} p_{(j)} \right) \right)$$

---

## 6. 실전 Python 코드: AnnData 기반 순열 검정 및 FDR 계산 엔진 구현

순열 검정의 개념을 순수 Python(NumPy/SciPy/Scanpy)으로 구현한 *교육용 변형*입니다. 점수는 곱, p-value는 $+1$ 보정을 사용하므로 CellPhoneDB(평균 점수, 복합체 최소값, $+1$ 없는 p-value)의 정의를 그대로 재현한 것이 아닙니다:

```python
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import scanpy as sc

def run_tme_permutation_test(
    adata,
    source_type,
    target_type,
    ligand_gene,
    receptor_gene,
    cluster_key='cell_type',
    n_perms=1000,
    min_pct=0.10,
    seed=42
):
    """
    단일 L-R 쌍 및 세포 쌍에 대한 질량작용 법칙 점수 산출 및 무작위 라벨 순열 검정을 수행합니다.
    """
    np.random.seed(seed)
    
    # 1. 유전자 존재 여부 확인
    for g in [ligand_gene, receptor_gene]:
        if g not in adata.var_names:
            raise ValueError(f"Gene '{g}' not found in AnnData object.")
            
    # 발현 행렬 추출 (dense array 변환)
    L_expr = adata[:, ligand_gene].X
    R_expr = adata[:, receptor_gene].X
    if hasattr(L_expr, "toarray"):
        L_expr = L_expr.toarray().flatten()
        R_expr = R_expr.toarray().flatten()
    else:
        L_expr = L_expr.flatten()
        R_expr = R_expr.flatten()
        
    labels = adata.obs[cluster_key].values
    
    # 2. 발현 백분율(min.pct) 필터링 검증
    src_mask = (labels == source_type)
    tgt_mask = (labels == target_type)
    
    pct_L = np.mean(L_expr[src_mask] > 0)
    pct_R = np.mean(R_expr[tgt_mask] > 0)
    
    if pct_L < min_pct or pct_R < min_pct:
        return {
            "interaction": f"{ligand_gene} -> {receptor_gene}",
            "source": source_type,
            "target": target_type,
            "score_obs": 0.0,
            "pct_L": pct_L,
            "pct_R": pct_R,
            "p_value": 1.0,
            "passed_pct_filter": False
        }
        
    # 3. 관측 상호작용 점수 (질량작용 법칙: Mean_L * Mean_R)
    mean_L_obs = np.mean(L_expr[src_mask])
    mean_R_obs = np.mean(R_expr[tgt_mask])
    score_obs = mean_L_obs * mean_R_obs
    
    # 4. 세포 라벨 무작위 순열 검정 (Permutation Loop)
    perm_scores = np.empty(n_perms, dtype=np.float64)
    shuffled_labels = labels.copy()
    
    for b in range(n_perms):
        np.random.shuffle(shuffled_labels)
        shuf_src_mask = (shuffled_labels == source_type)
        shuf_tgt_mask = (shuffled_labels == target_type)
        
        mean_L_b = np.mean(L_expr[shuf_src_mask])
        mean_R_b = np.mean(R_expr[shuf_tgt_mask])
        perm_scores[b] = mean_L_b * mean_R_b
        
    # 5. 엄격한 경험적 p-value 산출 (+1 Laplace pseudocount correction)
    count_greater_or_equal = np.sum(perm_scores >= score_obs)
    empirical_p = (count_greater_or_equal + 1.0) / (n_perms + 1.0)
    
    return {
        "interaction": f"{ligand_gene} -> {receptor_gene}",
        "source": source_type,
        "target": target_type,
        "score_obs": float(score_obs),
        "pct_L": float(pct_L),
        "pct_R": float(pct_R),
        "p_value": float(empirical_p),
        "passed_pct_filter": True
    }

# 6. 다중 가설 결과에 대한 Benjamini-Hochberg FDR 보정 유틸리티
def apply_fdr_correction(results_df, alpha=0.05):
    """
    모든 검정된 L-R 결과 DataFrame에 Benjamini-Hochberg FDR을 적용하여 q-value를 계산합니다.
    """
    p_vals = results_df['p_value'].values
    q_vals = false_discovery_control(p_vals, method='bh')
    results_df['q_value'] = q_vals
    results_df['significant'] = (results_df['q_value'] < alpha) & (results_df['passed_pct_filter'])
    return results_df.sort_values('q_value')

print("[Verified] TME Permutation Test Engine and BH-FDR correction successfully implemented.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### [문제 1. 핵심 개념 능동 회상 (Active Recall)]
세포 라벨 무작위 순열 검정(Permutation Testing)을 수행하여 경험적 $p$-value를 계산할 때, 단순 분수 $\frac{\sum \mathbb{I}(S^* \ge S_{\text{obs}})}{B}$ 대신 분자와 분모에 각각 1을 더한 공식인 $p = \frac{\sum_{b=1}^B \mathbb{I}(S_b^* \ge S_{\text{obs}}) + 1}{B + 1}$을 엄격히 적용해야 하는 통계학적 이유 2가지를 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **$p=0$의 수학적 모순 배제 (양수성 보장)**:
   - 만약 어떤 관측 상호작용 점수 $S_{\text{obs}}$가 $B=1,000$번의 순열에서 단 한 번도 나타나지 않았다고 해서($\sum \mathbb{I} = 0$), $p$-value를 $0$으로 판정하는 것은 표본 추출 이론상 명백한 오류입니다.
   - 영가설 하에서 관측값보다 더 극단적인 값이 나타날 확률이 '절대적으로 0'이라는 것은 불가능하며, 이는 단지 유한한 $B$번의 시도에서 관측되지 않았을 뿐입니다. 분자와 분모에 $+1$을 더함으로써 최소 $p$-value가 $\frac{1}{B+1}$ (예: $B=1000$일 때 약 $0.000999$)로 한정되어 수학적 모순을 방지합니다.
2. **관측된 실제 데이터의 귀무 표본 편입 및 제1종 오류 통제 (비편향성)**:
   - 영가설($H_0$)이 참이라는 가정 하에서, 실제 실험에서 얻은 '관측된 데이터' 역시 영분포(Null distribution)로부터 추출된 유효한 하나의 무작위 표본(Realization)에 해당합니다.
   - 따라서 총 표본의 수는 순열로 생성한 $B$개와 실제 관측 표본 1개를 합친 $B+1$개가 되며, 관측값 자신은 항상 자기 자신보다 크거나 같으므로($S_{\text{obs}} \ge S_{\text{obs}}$) 분자의 지시 함수 카운트에 $+1$이 자연스럽게 추가됩니다. 이를 통해 유한 표본에서 제1종 오류율(Type I Error Rate)이 명목 유의수준 $\alpha$를 초과하는 현상을 방지하는 보수적 비편향 추정량이 완성됩니다.

</details>

---

### [문제 2. 핵심 개념 능동 회상 (Active Recall)]
단일세포 RNA 시퀀싱 데이터에서 리간드와 수용체의 상호작용 세기를 점수화할 때, 단순 클러스터 평균 발현량 $\bar{E}$만을 계산하는 방식에 비해 **최소 발현 백분율 임계치($\text{min.pct} \ge 10\%\sim20\%$) 필터링**을 사전에 필수적으로 적용해야 하는 생물학적 이유와 통계적 이유를 각각 서술하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **생물학적 이유 (세포간질액 내 유효 농도 형성)**:
   - 주변분비(Paracrine) 신호 전달이 성립하려면 분비된 리간드가 세포간질액에 일정 농도 구배를 형성하여 표적 세포의 수용체 해리 상수($K_d$) 이상으로 도달해야 합니다.
   - 클러스터 내에서 극소수(예: $1\%$)의 세포만이 고립되어 리간드를 분비할 경우, 방출된 분자들은 광대한 세포외기질 속으로 희석되거나 신속히 분해되어 수용체를 점유할 수 있는 유효 생물학적 풀(Effective biological pool)을 형성하지 못합니다. 최소 $10\%\sim20\%$ 이상의 세포가 집단적으로 분비에 동참해야만 기능적 신호 도메인이 수립됩니다.
2. **통계적 이유 (이상치 세포 및 이중핵에 의한 평균치 왜곡 차단)**:
   - 단일세포 전사체 데이터에는 2개 이상의 서로 다른 세포가 하나의 액적에 포획된 **이중핵(Doublet)**이나, 비정상적으로 전사 활성이 폭증한 단일 이상치(Outlier) 세포가 포함될 수 있습니다.
   - 만약 1,000개 세포 중 단 1개의 이상치 세포가 유전자를 $1,000$ 카운트로 극단 발현하고 나머지 999개 세포는 $0$ 카운트라면, 단순 산술평균은 $\bar{E} = 1.0$으로 계산되어 마치 해당 클러스터 전체의 특징적인 신호 분자인 것처럼 심각한 위양성(False Positive)을 초래합니다. $\text{min.pct}$ 필터는 이러한 편향을 완전히 소거합니다.

</details>

---

### [문제 3. 임상 사례 분석 (Clinical Interpretation Case Study)]
**[임상 시나리오]**
골 전이성 전립선암(mCRPC) 환자 10명의 골수 단일세포 전사체 데이터(총 18개 세포 유형, 4,000개 공공 L-R 쌍)를 분석한 생명정보학 연구원이, 세포 쌍 간의 모든 가능한 상호작용(총 $18 \times 18 \times 4000 = 1,296,000$개 가설 검정)에 대해 순열 검정($B=1,000$)을 수행했습니다.
연구원은 개별 가설의 명목 $p$-value가 $p < 0.01$인 상호작용 12,000개를 도출한 뒤, "이 중 상위 500개는 뼈 전이 미세환경을 특이적으로 구동하는 통계적으로 결정적인 신호 전달 경로"라고 학술지에 투고했습니다.
그러나 심사위원(Reviewer)은 "다중 가설 검정(Multiple Hypothesis Testing)의 기본 원리를 무시한 심각한 방법론적 오류가 있어 신뢰할 수 없다"며 전면 기각(Reject)을 통보했습니다.

**[질문]**
1. 심사위원이 지적한 치명적 통계적 오류를 **위발견율(False Discovery Rate)과 위양성 개수 기댓값** 계산을 통해 정량적으로 증명하시오.
2. 이 연구원이 벤자미니-호흐바르크(Benjamini-Hochberg) FDR 보정을 전체 파이프라인에 올바르게 적용하기 위해 수행해야 하는 단계별 통계 처리 절차를 기술하고, 만약 $B=1,000$ 순열 검정 구조 하에서 발생할 수 있는 $p$-value 해상도 한계(Resolution Limit) 문제와 그 해결책을 제시하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

1. **통계적 오류의 정량적 증명**:
   - 총 검정 횟수 $M = 1,296,000$개입니다.
   - 만약 모든 세포 간 상호작용이 실제로는 아무런 의미가 없는 완전한 영가설($H_0$) 상태라 하더라도, 개별 검정에서 명목 유의수준 $\alpha = 0.01$을 적용하면 우연히 $p < 0.01$로 잘못 검출되는 **위양성(False Positive)의 기댓값**은 다음과 같습니다:
     $$\mathbb{E}[\text{False Positives}] = M \times \alpha = 1,296,000 \times 0.01 = 12,960\text{개}$$
   - 연구원이 유의하다고 도출한 12,000개라는 숫자는 귀무가설 하에서 우연히 발생할 수 있는 위양성 기댓값($12,960$개)보다도 적습니다. 모든 검정이 영가설($\pi_0 \approx 1$)이라면 추정 FDR $= \pi_0 \cdot 12{,}960 / 12{,}000 \approx 1$이므로 12,000개 신호의 대부분이 우연일 가능성이 높으며(실제 $\pi_0$는 미지이므로 단정은 불가), 이를 근거로 상위 500개를 결론짓는 것은 완전히 무효합니다.
2. **올바른 FDR 보정 절차 및 해상도 한계 해결책**:
   - **단계별 BH-FDR 처리 절차**:
     1. 발현 백분율 임계치($\text{min.pct} \ge 0.10$) 필터링을 사전에 적용하여 발현되지도 않는 무의미한 검정 후보(수십만 개)를 사전 제거함으로써 총 검정 수 $M'$을 대폭 축소시킵니다.
     2. 유효한 $M'$개 가설의 $p$-value를 오름차순으로 정렬 ($p_{(1)} \le p_{(2)} \dots$).
     3. $p_{(i)} \le \frac{i}{M'} \times 0.05$를 만족하는 최대 순위 $k$를 찾아 $q < 0.05$인 상호작용만을 통제된 유의 신호로 채택합니다.
   - **순열 검정 해상도 한계 및 해결책**:
     - $B=1,000$ 순열 검정에서 얻을 수 있는 최소 $p$-value는 $\frac{1}{1001} \approx 0.001$입니다.
     - 그러나 검정 수가 수만~수십만 개에 달할 경우, 엄격한 BH-FDR 임계치를 통과하려면 $p$-value가 $10^{-5}\sim10^{-7}$ 수준으로 극히 작아야 합니다. 최소값이 $0.001$에 불과하고 그 값으로 동률인 가설이 충분히 많지 않으면(전체 검정 수가 약 130만 개일 때는 동률이 약 2만 6천 개 이상이어야 BH가 기각 가능) $q$-value가 1에 가깝게 남아 실제 존재하는 신호가 탈락하는 해상도 결핍이 발생합니다.
     - **해결책**: 순열 반복 횟수를 $B=10,000\sim100,000$으로 대폭 늘리거나, 순열 영분포를 표준 정규분포 또는 감마 분포로 모수적 근사(Parametric approximation via z-score fitting)하여 연속적인 미세 $p$-value를 외삽 추정하는 하이브리드 통계 모델을 채택해야 합니다.

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Stat Appl Genet Mol Biol*. 2010;9:Article39. PMID: 21044043. doi: 10.2202/1544-6115.1585
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- NumPy 2.5.4 (numpy.org). 원 논문: Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. PMID: 32939066. doi: 10.1038/s41586-020-2649-2
- pandas 2.3.3 (github.com/pandas-dev/pandas). 원 논문: McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*. 2010:56-61. doi: 10.25080/majora-92bf1922-00a
- Scanpy 1.12.4 (github.com/scverse/scanpy). 원 논문: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15. PMID: 29409532. doi: 10.1186/s13059-017-1382-0
- SciPy 1.16.3 (github.com/scipy/scipy). 원 논문: Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. PMID: 32015543. doi: 10.1038/s41592-019-0686-2
<!-- REFERENCES:END -->

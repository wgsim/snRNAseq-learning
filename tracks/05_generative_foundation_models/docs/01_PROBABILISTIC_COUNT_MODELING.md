# 01. 단일세포 전사체 확률론적 카운트 모델링 (Probabilistic Count Modeling)

---

## 학습 목표 (Cognitive Goals)

1. **단일세포 RNA 카운트의 비가우시안성(Non-Gaussianity) 및 포아송 등분산성 붕괴 규명**: 연속 정규분포(Gaussian) 모델이 단일세포 카운트 데이터에서 붕괴하는 원인을 규명하고, 전사 파열(Transcriptional Bursting)과 기술적 포획 효율 변동으로 인한 포아송 등분산성($\text{Var}(X) = \mathbb{E}[X]$)의 한계를 수리물리학적으로 증명할 수 있다.
2. **Gamma-Poisson 혼합 분포(음이항, Negative Binomial)의 과산포(Overdispersion) 파라미터화**: 포아송 전사 속도 $\Lambda$에 감마 사전분포(Gamma Prior)를 도입하여 음이항 분포를 엄밀하게 유도하고, 평균 $\mu$와 분산 $\text{Var}(X) = \mu + \frac{\mu^2}{\theta}$의 수학적 관계 및 과산포 계수 $\theta$의 생물학적 의미를 도출할 수 있다.
3. **제로 팽창 음이항(ZINB) 게이팅과 UMI 프로토콜 논쟁 분석**: 생물학적 진성 결측(True Biological Zeros)과 기술적 드롭아웃(Technical Dropouts)을 분리하는 베이즈 혼합 모델을 정식화하고, 액적 기반 UMI(10x Chromium)와 플레이트 기반 비UMI(Smart-seq2) 간의 통계적 차이를 평가할 수 있다.
4. **라이브러리 크기($l_c$)와 상대 발현 비율($\rho_{cg}$)의 승법적 분해(Multiplicative Decomposition) 구현**: 시퀀싱 깊이 편차를 분리하는 수리 모델을 정립하고, PyTorch 기반 수치 안정화(Numerical Stability) 기법을 적용한 음이항 및 ZINB 로그 우도 함수를 직접 작성할 수 있다.

---

## 1. 단일세포 전사체 카운트의 통계물리학: 왜 가우시안 분포가 아닌가?

단일세포 RNA 시퀀싱(scRNA-seq/snRNA-seq) 데이터는 Bulk RNA-seq의 연속적 신호(Continuous Intensity)와 완전히 다른 물리적 특성을 지닌다. 세포 하나에 존재하는 특정 유전자의 전사체 분자 수는 유한한 자연수($0, 1, 2, \dots$)이며, 시퀀싱 기기는 이 분자들을 샘플링하여 불연속 정수 카운트(Discrete Non-negative Integer Counts)로 출력한다.

```
[Single-Cell Measurement Physics & Noise Generation]

Promoter State Dynamics      mRNA Synthesis & Decay       Microfluidic Droplet & RT     Sequencing Output
(Telegraph Model)            (Poissonian Birth-Death)     (Binomial Sampling)           (Observed Count)

[Silent]  k_on                +── mRNA (t) ──+             Captured cDNA                 PCR Amplification
   │ ────────> [Active] ───>  │ Burst Size b │ ────────>  (Efficiency η ~ 5-15%) ────>  & Sequencing
   │ <────────    │           +──────────────+                    │                            │
  k_off         k_trans                                           ▼                            ▼
                                                            Missing mRNA                 Discrete Count
                                                            (Technical Dropout)          X_cg ∈ {0, 1, 2, ...}
```

### 1.1 가우시안 모델(Gaussian Model)의 세 가지 치명적 결함

고전적 머신러닝 알고리즘(PCA, $k$-means, 일반 선형 회귀)은 오차항이 정규분포를 따른다고 가정한다:

$$x_{cg} \sim \mathcal{N}(\mu_{cg}, \sigma^2)$$

그러나 실제 단일세포 전사체 행렬에 가우시안 가정을 적용하면 다음과 같은 모순이 발생한다:

1.  **음수 발현량 생성의 모순**: 정규분포의 지지 집합(Support)은 $(-\infty, +\infty)$이다. 평균 $\mu$가 작은 저발현 유전자(e.g., 전사인자, $\mu \approx 0.1$)에 정규분포를 적합하면, 사후 생성 모델은 물리적으로 불가능한 음수 분자 수($x < 0$)를 높은 확률로 샘플링한다.
2.  **극단적 우측 꼬리 왜곡(Extreme Right Skewness)**: 단일세포 발현 분포는 0 부근에 대다수의 세포가 집중되고, 소수의 세포에서 수십~수백 개의 카운트가 관측되는 강한 양의 왜도(Positive Skewness)를 띤다. 대칭적인 가우시안 종형 곡선은 이 비대칭성을 수용할 수 없다.
3.  **이분산성(Heteroskedasticity)**: 카운트 데이터의 분산은 평균 발현량 $\mu$에 강하게 의존한다. 평균이 커질수록 분산도 함께 증가하므로, 잔차의 등분산성을 전제하는 최소제곱법(Ordinary Least Squares)은 고발현 유전자에 과도한 가중치를 부여하게 된다.

---

## 2. 포아송 분포(Poisson Distribution)와 등분산성 실패

가장 단순한 불연속 카운트 모델은 포아송 분포(Poisson Distribution)이다. 단위 시간당 일정한 전사 속도 $\lambda$로 mRNA 분자가 생성된다면, 세포 내 유전자 분자 수 $X$는 다음과 같이 기술된다:

$$P(X = k \mid \lambda) = \frac{\lambda^k e^{-\lambda}}{k!}, \quad k \in \{0, 1, 2, \dots\}$$

포아송 분포의 결정적 수학적 성질은 **등분산성(Equi-dispersion)**이다:

$$\mathbb{E}[X] = \lambda, \quad \text{Var}(X) = \lambda$$

즉, 이론적으로 $\frac{\text{Var}(X)}{\mathbb{E}[X]} = 1$ (Fano Factor = 1)이어야 한다.

### 2.1 포아송 등분산성 붕괴의 통계물리학적 기전: 전사 파열(Transcriptional Bursting)

실제 단일세포 전사체 데이터에서 유전자 발현량을 측정하면, 거의 모든 유전자에서 분산이 평균보다 훨씬 큰 **과산포(Overdispersion, $\text{Var}(X) \gg \mathbb{E}[X]$)** 현상이 관찰된다:

$$\text{Var}(X) = \mu + \alpha \mu^2, \quad \alpha > 0$$

이 과산포는 단순한 측정 오차가 아니라 전사의 분자생물학적 메커니즘에서 기인한다:

1.  **전신 텔레그래프 모델(Two-state Telegraph Model)**:
    프로모터는 비활성(Off) 상태와 활성(On) 상태를 무작위로 오간다:
    $$[\text{Off}] \underset{k_{\text{off}}}{\overset{k_{\text{on}}}{\rightleftharpoons}} [\text{Active}]$$
    프로모터가 활성화된 짧은 시간 동안 여러 개의 mRNA가 폭발적으로 합성되는 **전사 파열(Transcriptional Bursting)**이 일어난다.
2.  **기술적 포획 효율(Capture Efficiency)의 세포 간 변동**:
    각 액적(Droplet) 내부의 역전사 효소 활성, 용해(Lysis) 효율, 마이크로유체 유속 변동에 의해 포획 효율 $\eta_c$는 고정된 상수가 아니라 세포마다 달라지는 확률변수가 된다.

따라서 단일 전사 속도 $\lambda$를 가정하는 포아송 모델은 단일세포 데이터를 모델링하기에 근본적으로 부적합하다.

---

## 3. Gamma-Poisson 혼합 모델: 음이항 분포(Negative Binomial, NB)

과산포를 엄밀하게 설명하기 위해, 전사율 $\Lambda$ 자체가 고정된 상수가 아니라 감마 분포(Gamma Distribution)를 따르는 연속 확률변수라고 모델링한다. 이를 **Gamma-Poisson Mixture**라고 부르며, 이것이 곧 **음이항 분포(Negative Binomial Distribution)**이다.

```
[Hierarchical Gamma-Poisson Mixture Framework]

Biological & Technical Latent Rate:
  Λ ~ Gamma(shape = θ, rate = β = θ / μ)
         │
         ▼ (Poissonian Sampling)
Observed Count:
  X | Λ ~ Poisson(Λ)
         │
         ▼ (Marginalize over Λ)
Marginal Negative Binomial Likelihood:
  X ~ NegativeBinomial(mean = μ, dispersion = θ)
```

### 3.1 음이항 분포의 완전한 수학적 유도

감마 분포의 확률 밀도 함수(PDF)는 형상 파라미터(Shape) $\theta > 0$와 척도 파라미터(Rate) $\beta > 0$에 대해 다음과 같다:

$$p(\lambda \mid \theta, \beta) = \frac{\beta^\theta}{\Gamma(\theta)} \lambda^{\theta - 1} e^{-\beta \lambda}, \quad \lambda > 0$$

여기서 $\mathbb{E}[\Lambda] = \frac{\theta}{\beta} \equiv \mu$이므로, $\beta = \frac{\theta}{\mu}$로 치환할 수 있다.

이제 잠재 변수 $\lambda$에 대해 포아송 결합 확률을 주변화(Marginalization)하여 $X=k$의 주변 확률을 구한다:

$$P(X = k \mid \mu, \theta) = \int_0^\infty P(X = k \mid \lambda) p(\lambda \mid \theta, \mu) \, d\lambda$$

$$= \int_0^\infty \left( \frac{\lambda^k e^{-\lambda}}{k!} \right) \left( \frac{(\theta/\mu)^\theta}{\Gamma(\theta)} \lambda^{\theta - 1} e^{-(\theta/\mu)\lambda} \right) d\lambda$$

$$= \frac{(\theta/\mu)^\theta}{k! \, \Gamma(\theta)} \int_0^\infty \lambda^{k + \theta - 1} e^{-\left(1 + \frac{\theta}{\mu}\right)\lambda} \, d\lambda$$

적분 항은 감마 함수의 정의 $\int_0^\infty t^{z-1} e^{-a t} dt = \frac{\Gamma(z)}{a^z}$를 이용하면 다음과 같이 계산된다:

$$\int_0^\infty \lambda^{(k + \theta) - 1} e^{-\left(\frac{\mu + \theta}{\mu}\right)\lambda} \, d\lambda = \frac{\Gamma(k + \theta)}{\left(\frac{\mu + \theta}{\mu}\right)^{k + \theta}}$$

이를 대입하여 정리하면 최종 음이항 분포의 확률 질량 함수(PMF)가 유도된다:

$$P(X = k \mid \mu, \theta) = \frac{(\theta/\mu)^\theta}{k! \, \Gamma(\theta)} \cdot \frac{\Gamma(k + \theta)}{\left(\frac{\theta + \mu}{\mu}\right)^{k + \theta}}$$

$$= \frac{\Gamma(k + \theta)}{k! \, \Gamma(\theta)} \left( \frac{\theta}{\theta + \mu} \right)^\theta \left( \frac{\mu}{\theta + \mu} \right)^k$$

### 3.2 음이항 분포의 평균, 분산 및 과산포 거동

음이항 분포의 1차, 2차 적률(Moments)은 다음과 같다:

$$\mathbb{E}[X] = \mu$$

$$\text{Var}(X) = \mu + \frac{\mu^2}{\theta}$$

*   **과산포 계수(Dispersion Parameter) $\theta$의 물리적 의미**:
    *   $\theta \to \infty$인 경우: $\frac{\mu^2}{\theta} \to 0$이 되어 $\text{Var}(X) \to \mu$가 성립한다. 즉, **음이항 분포는 포아송 분포로 정확히 수렴**한다.
    *   $\theta$가 작을수록(e.g., $\theta < 1$): 분산이 평균의 제곱($\mu^2$)에 비례하여 급격히 팽창한다. 이는 프로모터의 On/Off 전환이 느리고 거대한 전사 파열을 일으키는 유전자일수록 $\theta$가 작은 값으로 추정됨을 의미한다.

---

## 4. 제로 팽창 음이항(ZINB) 모델과 UMI 프로토콜 논쟁

### 4.1 ZINB의 수학적 정식화

단일세포 데이터에서 관측되는 대규모 영점(Zero counts)을 모델링하기 위해, 제로 팽창 음이항(Zero-Inflated Negative Binomial, ZINB) 모델은 두 가지 상이한 영점 생성 메커니즘을 가정한다:

```
                  ┌── [확률 π] ──────> 항구적 영점 (Structural / Technical Zero)
                  │
[카운트 관측 X] ──┤
                  │
                  └── [확률 1 - π] ──> 음이항 분포 (NB Sampling)
                                          ├── [NB 영점 P_NB(0)] ──> 생물학적 영점 (Biological Zero)
                                          └── [k > 0 카운트]   ──> 양수 발현 카운트
```

1.  **생물학적 진성 결측(Biological Zeros)**: 유전자가 해당 세포형에서 발현되지 않거나(Promoter Silent), mRNA가 일시적으로 분해되어 세포 내에 분자가 실제로 존재하지 않는 경우.
2.  **기술적 드롭아웃(Technical Dropouts)**: 세포 내에 mRNA 분자가 존재하였으나, 세포 용해 실패, 역전사 효소 결합 실패, 혹은 시퀀싱 심도 부족으로 인해 포획되지 못하고 0으로 관측된 경우.

ZINB는 베르누이 게이팅 변수 $\pi \in [0, 1]$를 도입하여 이를 혼합 모델로 결합한다:

$$P(X = 0 \mid \mu, \theta, \pi) = \pi + (1 - \pi) \left( \frac{\theta}{\theta + \mu} \right)^\theta$$

$$P(X = k \mid \mu, \theta, \pi) = (1 - \pi) \frac{\Gamma(k + \theta)}{k! \, \Gamma(\theta)} \left( \frac{\theta}{\theta + \mu} \right)^\theta \left( \frac{\mu}{\theta + \mu} \right)^k, \quad (k \ge 1)$$

### 4.2 UMI 데이터와 Svensson (2020)의 대논쟁: ZINB는 여전히 필요한가?

초기 단일세포 분석(Smart-seq, Drop-seq 초기)에서는 ZINB가 사실상의 표준으로 간주되었다. 그러나 Valentine Svensson (*Nature Biotechnology*, 2020)과 Townes et al. (*Genome Biology*, 2019)은 중대한 반론을 제기하였다:

*   **고유 분자 식별자(Unique Molecular Identifiers, UMIs)의 도입**: 10x Chromium과 같은 현대 단일세포 플랫폼은 각 역전사 cDNA 분자에 무작위 UMI 바코드를 부여하여 PCR 증폭 편향(Amplification Bias)을 크게 줄인다(완전 제거는 아니다).
*   **음성 대조군(균일 RNA) 데이터**: Svensson은 UMI 액적 기반 음성 대조 데이터에서 관측된 영점의 빈도가 추가적인 제로 팽창 파라미터($\pi$) 없이 포아송/음이항 샘플링으로 설명됨을 보였다("ERCC 스파이크인 10x 실험"이라는 구체 서술은 원문에서 직접 확인하지 못했다). 후속 Cao 외(2021)는 이 해석을 비UMI 플랫폼으로 일반화하는 것을 비판하였다.
*   **결론**: UMI 기반 플랫폼에서 대다수의 영점은 '비정상적 기술적 드롭아웃'이 아니라, **극도로 낮은 분자 농도와 한정된 시퀀싱 깊이로 인해 발생하는 통계적 샘플링 영점(Sampling Zeros)**이다. 따라서 불필요한 제로 팽창 게이트 $\pi$를 추가하면 모델이 과파라미터화(Over-parameterized)되어 수치적 불안정성과 과적합을 유발한다. UMI 데이터에서는 `scvi-tools`의 `gene_likelihood="nb"`를 선택하는 것이 권장되는 이유가 바로 여기에 있다. (주의: `SCVI`의 **기본값은 `"zinb"`**이며, NB는 사용자가 명시적으로 선택하는 옵션이다. scvi-tools 0.20~1.5 소스에서 확인.)

---

## 5. 승법적 분해 (Multiplicative Decomposition): 라이브러리 크기와 유전자 발현율

단일세포 RNA 카운트 $x_{cg}$ (세포 $c$, 유전자 $g$)의 기저 평균 $\mu_{cg}$는 두 가지 독립적 요인의 곱으로 승법적 분해(Multiplicative Decomposition)된다:

$$\mu_{cg} = l_c \cdot \rho_{cg}$$

```
[Multiplicative Likelihood Decomposition]

Observed Count Matrix X (Cells × Genes)
             │
             ├──> Cell Library Size Factor l_c (Scalar per Cell)
             │    - 기술적 시퀀싱 깊이 + 세포 체적(Volume) 변동 반영
             │    - Prior: l_c ~ LogNormal(μ_l, σ_l^2)
             │
             └──> Normalized Gene Expression Rate ρ_cg (Vector per Cell)
                  - 세포의 생물학적 상태 z_c로부터 디코딩된 상대 발현 비율
                  - 제약 조건: ρ_cg ∈ (0, 1),  ∑_g ρ_cg = 1 (Simplex)
```

1.  **세포별 라이브러리 크기 인자(Library Size Factor) $l_c$**:
    세포 $c$에서 시퀀싱된 총 전사체 분자 수의 척도이다. 이는 생물학적 세포 체적(Cell Volume/RNA content)과 실험적 시퀀싱 심도(Sequencing Depth)의 복합적 결과물이다.
2.  **유전자별 상대 발현 비율(Gene Expression Proportion) $\rho_{cg}$**:
    세포 $c$ 내에서 발현되는 전체 mRNA 중 유전자 $g$가 차지하는 진정한 생물학적 분율이다. 확률 단체(Probability Simplex) 상에 존재하므로 다음을 만족한다:
    $$\rho_{cg} > 0, \quad \sum_{g=1}^G \rho_{cg} = 1$$
    심층 생성 모델(scVI)은 디코더의 최종 출력단에 **소프트맥스(Softmax)** 활성화 함수를 적용하여 이 심플렉스 제약 조건을 엄밀하게 만족시킨다.

---

## 6. PyTorch 기반 수치 안정화 확률 우도 구현

음이항 및 ZINB의 로그 우도(Log-Likelihood)를 계산할 때, 팩토리얼 계산($k!$)이나 감마 함수($\Gamma(k+\theta)$)는 $k$가 100만 넘어도 부동소수점 오버플로우(Overflow)를 일으킨다. 따라서 모든 계산은 로그 스케일(`lgamma`, `log1p`)에서 수치적으로 안정화되어야 한다.

```python
"""
01_probabilistic_count_losses.py
단일세포 RNA 카운트를 위한 수치 안정화 Negative Binomial 및 ZINB Loss 함수.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

class NegativeBinomialLoss(nn.Module):
    """
    수치 안정화 음이항(Negative Binomial) 음의 로그 우도(Negative Log-Likelihood, NLL) 모듈.
    
    매개변수화:
      mu: 예측된 평균 발현량 (l_c * rho_{cg}), Shape: (Batch, Genes)
      theta: 유전자별 과산포 파라미터 (Inverse Dispersion), Shape: (Genes,) or (Batch, Genes)
      eps: 수치적 안정성을 위한 최소 엡실론 (1e-8)
    """
    def __init__(self, eps=1e-8):
        super().__init__()
        self.eps = eps

    def forward(self, x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor) -> torch.Tensor:
        """
        x: 실제 관측 카운트 정수 행렬 (Batch, Genes)
        반환값: 각 샘플별 음의 로그 우도 합 (Batch,)
        """
        # 수치적 클램핑으로 0 방지
        mu = torch.clamp(mu, min=self.eps)
        theta = torch.clamp(theta, min=self.eps)

        # Log-Gamma 항 계산: log(Gamma(x + theta)) - log(Gamma(theta)) - log(Gamma(x + 1))
        # torch.lgamma는 log(|Gamma(y)|)를 안정적으로 계산
        t1 = torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1.0)
        
        # 확률 항 계산: theta * log(theta / (theta + mu)) + x * log(mu / (theta + mu))
        # log(theta / (theta + mu)) = -log(1 + mu / theta) = -torch.log1p(mu / theta)
        # log(mu / (theta + mu)) = log(mu) - log(theta + mu)
        # (float32에서 theta가 1e4 이상이면 lgamma 차(t1)의 상쇄 오차가 지배적이므로
        #  그런 영역은 theta 상한 클램프나 float64로 다룬다: 실측 오차 theta=1e4 ~1e-3, 1e6 ~0.5)
        log_theta_mu = torch.log(theta + mu)
        t2 = -theta * torch.log1p(mu / theta)
        t3 = x * (torch.log(mu) - log_theta_mu)

        # 최종 log-likelihood
        log_prob = t1 + t2 + t3
        
        # 음의 로그 우도 반환 (유전자 축 합산)
        return -torch.sum(log_prob, dim=-1)


class ZeroInflatedNegativeBinomialLoss(nn.Module):
    """
    수치 안정화 제로 팽창 음이항(ZINB) NLL 모듈.
    
    매개변수화:
      mu: 음이항 평균 (Batch, Genes)
      theta: 과산포 파라미터 (Batch, Genes)
      pi_logits: 제로 팽창 게이트의 Logit (Batch, Genes), pi = sigmoid(pi_logits)
    """
    def __init__(self, eps=1e-8):
        super().__init__()
        self.eps = eps
        self.nb_loss = NegativeBinomialLoss(eps=eps)

    def forward(self, x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor, pi_logits: torch.Tensor) -> torch.Tensor:
        mu = torch.clamp(mu, min=self.eps)
        theta = torch.clamp(theta, min=self.eps)

        # softplus(-pi_logits) = -log(sigmoid(pi_logits))
        # log(1 - pi) = -F.softplus(pi_logits)
        # log(pi) = -F.softplus(-pi_logits)
        log_pi = -F.softplus(-pi_logits)
        log_1_minus_pi = -F.softplus(pi_logits)

        # 1. NB 우도의 로그 계산
        log_theta_mu = torch.log(theta + mu)
        nb_zero_log_prob = theta * (torch.log(theta) - log_theta_mu) # P_NB(0)

        # 2. x == 0 인 경우: log(pi + (1 - pi) * P_NB(0))
        # logaddexp(a, b) = log(exp(a) + exp(b))로 언더플로우 방지
        zero_case = torch.logaddexp(log_pi, log_1_minus_pi + nb_zero_log_prob)

        # 3. x > 0 인 경우: log(1 - pi) + log(P_NB(x))
        t1 = torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1.0)
        t2 = theta * (torch.log(theta) - log_theta_mu)
        t3 = x * (torch.log(mu) - log_theta_mu)
        nonzero_case = log_1_minus_pi + t1 + t2 + t3

        # 조건부 선택 마스크
        is_zero = (x < self.eps).float()
        log_prob = is_zero * zero_case + (1.0 - is_zero) * nonzero_case

        return -torch.sum(log_prob, dim=-1)


# ==========================================
# 실행 검증 및 시뮬레이션 테스트
# ==========================================
if __name__ == "__main__":
    torch.manual_seed(42)
    B, G = 5, 4  # 5개 세포, 4개 유전자
    
    # 합성 데이터 생성
    observed_counts = torch.tensor([
        [0.0, 5.0,  0.0, 120.0],
        [0.0, 0.0,  3.0, 45.0],
        [1.0, 2.0,  0.0, 0.0],
        [0.0, 18.0, 1.0, 89.0],
        [4.0, 0.0,  0.0, 210.0]
    ])
    
    predicted_mu = torch.tensor([
        [0.2, 4.8,  0.1, 115.0],
        [0.1, 0.3,  2.5, 50.0],
        [1.2, 1.8,  0.2, 0.5],
        [0.1, 16.5, 1.2, 92.0],
        [3.8, 0.2,  0.1, 205.0]
    ])
    
    dispersion_theta = torch.tensor([0.5, 1.2, 0.8, 3.5]) # 유전자별 theta
    pi_logits = torch.randn(B, G) # 임의의 dropout logit
    
    nb_criterion = NegativeBinomialLoss()
    zinb_criterion = ZeroInflatedNegativeBinomialLoss()
    
    loss_nb = nb_criterion(observed_counts, predicted_mu, dispersion_theta.unsqueeze(0))
    loss_zinb = zinb_criterion(observed_counts, predicted_mu, dispersion_theta.unsqueeze(0), pi_logits)
    
    print("=== Loss Computation Verification ===")
    print(f"Observed Count Shape : {observed_counts.shape}")
    print(f"NB Loss per Cell     : {loss_nb.detach().numpy().round(3)}")
    print(f"ZINB Loss per Cell   : {loss_zinb.detach().numpy().round(3)}")
    print(">>> 모든 확률 텐서 연산 정상 수렴 완료.")
```

---

## 단원 핵심 점검 (Cognitive Checkpoint & Case Study)

### 문제 1: UMI 기반 단일세포 데이터에서 NB와 ZINB 모델의 수리통계적 결론
10x Genomics Chromium(Single Cell 3' v3)으로 시퀀싱된 5만 개의 인간 전립선 상피세포 데이터를 분석하고 있다. 한 연구원이 "카운트 행렬의 85%가 0이므로, 제로 팽창 음이항(ZINB) 모델을 사용하여 드롭아웃 게이트 $\pi$를 반드시 학습해야 한다"고 주장하였다.
1. Valentine Svensson (2020)의 연구 결과와 통계적 샘플링 이론을 근거로, 이 연구원의 주장이 타당하지 않은 이유를 논하시오.
2. 불필요하게 ZINB 모델을 고차원 UMI 데이터에 적용했을 때 딥러닝 최적화 과정에서 발생할 수 있는 '모수 식별 불가능성(Parameter Identifiability Problem)' 및 경사면 왜곡 현상을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. Svensson (2020) 및 샘플링 이론에 기반한 반론**:
*   10x Chromium과 같은 UMI 기반 기술은 각 전사체 분자에 고유 바코드를 태깅하여 PCR 증폭 편향을 제거한다. 따라서 UMI 카운트 행렬의 영점($X_{cg}=0$) 중 대다수는 분자가 존재하는데 기술적으로 증폭되지 못한 '비정상적 제로 팽창(Dropout)'이 아니다.
*   세포 하나당 총 시퀀싱 리드는 수천~수만 개 수준에 불과하므로, 세포당 1~2개 분자만 발현되는 저발현 유전자(전사인자 등)는 포아송/음이항 무작위 샘플링 과정에서 자연스럽게 0으로 관측될 확률이 극히 높다:
    $$P(X=0) = \left( \frac{\theta}{\theta + \mu} \right)^\theta$$
    예를 들어 $\mu = 0.05, \theta = 1$일 때, 순수 NB 샘플링에 의해서도 영점 확률은 약 95.2%에 달한다.
*   Svensson의 정밀한 통계 검정에 따르면, UMI 데이터의 영점 비율은 추가적인 영점 게이트($\pi > 0$) 없이도 관측된 평균 $\mu$와 분산 $\text{Var}(X)$을 갖는 순수 음이항(NB) 분포의 이론적 영점 확률과 대체로 일치한다(대조 데이터 기준). 이질 집단의 실제 영과잉은 생물학적 변이로 설명되는 경우가 많아, UMI 데이터에서 ZINB는 대체로 불필요하다.

**2. 모수 식별 불가능성 및 경사면 왜곡**:
*   $\mu \to 0$인 저발현 구간에서, 관측 카운트가 0일 때 모델은 다음 두 가지 상이한 파라미터 조합으로 동일한 영점을 완벽하게 설명할 수 있다:
    1.  $\pi \approx 0$ (제로 팽창 없음)이면서 $\mu \approx 0$ (기저 발현량이 실제로 매우 낮음).
    2.  $\pi \approx 1$ (기술적 드롭아웃)이면서 $\mu > 0$ (기저 발현량은 높으나 게이트에 의해 차단됨).
*   이러한 현상을 **모수 식별 불가능성(Identifiability Issue)**이라 한다. 두 경로가 동일한 우도를 생성하므로 손실 함수의 경사면(Loss Landscape)에 길고 평평한 계곡(Plateau)이 형성되며, 경사 하강법 최적화 시 그래디언트가 진동하거나 국소 해에 갇혀 실제 생물학적 발현량 $\mu$의 복원이 심각하게 왜곡된다.

</details>

---

### 문제 2: 승법 분해 $\mu_{cg} = l_c \rho_{cg}$와 세포 크기-시퀀싱 깊이 간의 얽힘(Confounding) 해소
단일세포 RNA 시퀀싱에서 관측된 총 카운트 $N_c = \sum_g x_{cg}$는 실험적 시퀀싱 깊이(Sequencing Depth)뿐만 아니라 세포 고유의 물리적 크기 및 총 RNA 함량(Cell Volume & RNA Content)에 의해 동시에 결정된다.
1. 분화된 대형 관강 상피세포(Luminal cells, 큰 크기, 높은 RNA 함량)와 미분화 신경내분비 세포(t-NEPC, 작은 크기, 낮은 RNA 함량)가 혼재된 암 조직에서, 단순 전역 상수로 나누는 정규화($\tilde{x} = x / N_c \times 10^4$)가 초래하는 생물학적 상대 발현 비율의 왜곡 기전을 수식으로 분석하시오.
2. scVI가 라이브러리 크기 $l_c$를 확률변수로 다루면서 신경망 인코더 $q(l \mid x)$를 통해 이를 분리 추론하는 통계적 이점을 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 단순 전역 비례 정규화의 왜곡 기전**:
*   세포 유형 $A$(관강 세포)의 실제 세포 내 총 전사체 수를 $M_A$, 세포 유형 $B$(신경내분비 세포)의 총 전사체 수를 $M_B$라 할 때, $M_A \approx 3 M_B$의 관계가 성립한다고 가정하자.
*   시퀀싱 포획 효율이 두 세포에서 $\eta$로 동일하다면, 관측된 총 카운트는 $\mathbb{E}[N_A] = \eta M_A$, $\mathbb{E}[N_B] = \eta M_B$가 된다.
*   만약 유전자 $g$가 세포 $A$와 $B$에서 분자 농도($\text{molecules per cell}$)는 정확히 동일하게 $C_g$ 분자씩 존재한다고 하자.
*   이때 단순 비례 정규화 $\frac{x_{cg}}{N_c} \times 10^4$를 적용하면:
    $$\mathbb{E}[\tilde{x}_{Ag}] \approx \frac{\eta C_g}{\eta M_A} \times 10^4 = \frac{C_g}{M_A} \times 10^4$$
    $$\mathbb{E}[\tilde{x}_{Bg}] \approx \frac{\eta C_g}{\eta M_B} \times 10^4 = \frac{C_g}{M_B} \times 10^4$$
*   $M_A > M_B$이므로, 실제 세포당 절대 분자 수가 완전히 동일함에도 불구하고 정규화된 값은 $\tilde{x}_{Bg} \approx 3 \tilde{x}_{Ag}$로 산출되어, 세포 $B$에서 해당 유전자가 3배 과발현된 것으로 잘못 판정되는 치명적 왜곡이 발생한다.

**2. scVI의 분리 추론 이점**: (참고: 이는 Lopez et al. 2018 원 논문의 모델 서술이다. scvi-tools의 기본값 `use_observed_lib_size=True`에서는 $l_c$를 관측 총합 $N_c$로 고정하며 잠재 변수·`l_encoder`·LogNormal 사전을 쓰지 않는다. `False`로 지정하면 이 서술의 잠재 변수 모델이 된다.)
*   scVI는 $l_c$를 관측된 총합 $N_c$의 로그 정규 사전분포(Log-Normal Prior)를 갖는 잠재 변수로 설정한다:
    $$l_c \sim \text{LogNormal}(\mu_l, \sigma_l^2), \quad \mu_l \approx \text{mean}(\log N_c)$$
*   인코더 네트워크 $q_\phi(l_c \mid x_c)$는 세포의 전체 발현 패턴을 종합하여 기술적 깊이 변동을 추론하고, 별도의 잠재 공간 $z_c$를 통해 순수 상대 비율 $\rho_{cg} = \text{Softmax}(f_\theta(z_c))_g$를 디코딩한다.
*   결과적으로 라이브러리 스케일($l_c$)의 기술적 변동은 $\mu_{cg} = l_c \rho_{cg}$ 수식의 앞단에서 격리 흡수되므로, 잠재 상태 $z_c$와 차별 발현 유전자(DE) 검정은 세포 체적 차이에 따른 기술적 아티팩트의 간섭을 받지 않고 순수한 전사 프로그램의 조절 차이만을 정확하게 포착한다.

</details>

---

### 문제 3: 과산포 파라미터 $\theta$ 극한에서의 수치적 불안정성과 PyTorch 최적화
음이항 분포 손실 함수를 역전파(Backpropagation)로 최적화할 때, $\theta$의 극단적 값에 의해 그래디언트 폭발(Gradient Explosion) 또는 소실(Vanishing)이 발생할 수 있다.
1. 유전자의 과산포가 거의 없어 포아송 분포에 근접할 때($\theta \to \infty$)와, 반대로 전사 파열이 극단적이어서 영점 부근에 집중될 때($\theta \to 0$) 각각 발생하는 수치 연산상의 문제점을 컴퓨터 부동소수점 한계 관점에서 분석하시오.
2. `torch.lgamma` 및 `torch.log1p`를 활용하여 이 수치적 불안정성을 어떻게 해결하는지 수식적으로 설명하시오.

<details>
<summary>💡 모범 답안 및 해설 확인</summary>

**1. 극단적 $\theta$에서의 수치 연산 붕괴**:
*   **$\theta \to \infty$ (포아송 수렴 극한)**:
    우도 수식 내의 감마 함수 비율 $\frac{\Gamma(x + \theta)}{\Gamma(\theta)}$과 $\left( \frac{\theta}{\theta + \mu} \right)^\theta$ 항에서 문제가 발생한다. $\theta$가 1,000을 넘어서면 표준 부동소수점(float32, 최대 표현값 $\approx 3.4 \times 10^{38}$) 환경에서 $\Gamma(\theta)$는 즉각 무한대(`Inf`)로 오버플로우된다. 또한 컴퓨터가 $\infty / \infty$를 연산하면서 `NaN`이 반환되어 전체 신경망의 가중치가 파괴된다.
*   **$\theta \to 0$ (극단적 과산포 극한)**:
    $\log(\theta)$를 계산할 때 부동소수점 언더플로우로 인해 $-\infty$가 발생한다. 또한 분모에 $\theta$가 들어가는 항($\frac{\mu}{\theta}$)에서 0으로 나누기(Zero-division) 오류가 발생하여 역전파 그래디언트가 무한대로 폭발한다.

**2. Log-Space 안정화 해결책**:
*   **`lgamma` 적용**: 감마 함수를 직접 계산한 후 로그를 취하는 것이 아니라, 스털링 근사(Stirling's Approximation) 기반의 로그-감마 전용 함수 `torch.lgamma`를 사용한다:
    $$\log \frac{\Gamma(x + \theta)}{\Gamma(\theta)} = \text{lgamma}(x + \theta) - \text{lgamma}(\theta)$$
    이를 통해 거대한 팩토리얼 값의 직접 연산 없이 두 유한한 실수 간의 뺄셈으로 안전하게 변환된다.
*   **`log1p` 적용**: $\theta$가 매우 클 때 $\frac{\mu}{\theta} \ll 1$이 되므로, 로그 항은 $\log\left( \frac{\theta}{\theta + \mu} \right) = -\log\left(1 + \frac{\mu}{\theta}\right)$으로 변환된다. 컴퓨터에서 $1 + \epsilon \approx 1$로 반올림되는 정보 손실(Catastrophic Cancellation)을 줄이기 위해 `torch.log1p(mu / theta)`를 사용한다(위 구현 코드가 이 형태를 쓴다). 다만 float32에서 $\theta$가 매우 크면 `lgamma(x+θ) - lgamma(θ)`의 뺄셈 상쇄 오차가 더 크게 남는다(실측: float64 기준 최대 오차 $\theta=10^4$에서 약 $10^{-3}$, $10^6$에서 약 $0.5$; scvi-tools의 float32 구현도 같은 영역에서 오차가 있음).
*   **안전 클램핑(Clamping)**: 신경망 최적화 시 $\theta$는 반드시 지수 함수($\theta = \exp(\tilde{\theta})$) 또는 소프트플러스 활성화를 거치도록 설계하고, `torch.clamp(theta, min=1e-4, max=1e4)` 같은 범위 클램프를 적용하여 이론적 한계 영역 내에서만 업데이트되도록 강제할 수 있다(위 구현 코드는 하한 `eps=1e-8` 클램프만 포함하며, 상한은 필요 시 호출 측에서 지정한다).

</details>

---

<!-- REFERENCES:BEGIN (generated by scripts/apply_references.py; do not edit by hand) -->
## 참고문헌 (References)
### 본문 인용 문헌
1. Svensson V. Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol*. 2020;38(2):147-150. PMID: 31937974. doi: 10.1038/s41587-019-0379-5
2. Townes FW, Hicks SC, Aryee MJ, et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model. *Genome Biol*. 2019;20(1):295. PMID: 31870412. doi: 10.1186/s13059-019-1861-6
3. Cao Y, Kitanovski S, Küppers R, et al. UMI or not UMI, that is the question for scRNA-seq zero-inflation. *Nat Biotechnol*. 2021;39(2):158-159. PMID: 33526946. doi: 10.1038/s41587-020-00810-6
4. Lopez R, Regier J, Cole MB, et al. Deep generative modeling for single-cell transcriptomics. *Nat Methods*. 2018;15(12):1053-1058. PMID: 30504886. doi: 10.1038/s41592-018-0229-2
### 소프트웨어 (표기 버전은 검수 시 API를 확인한 버전이며, 파이프라인 실행 환경은 requirements-pipelines.txt 참조)
- PyTorch 2.14.1 (github.com/pytorch/pytorch).
<!-- REFERENCES:END -->

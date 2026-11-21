# scripts/course_markdown.py 통합 후 수식 렌더링 검토 (T01–T05)

- 비교 기준: HEAD `fa0cb97` 대 작업 트리. 검토 시점의 T01–T05 HTML을 스냅샷(`scratchpad/sci/snap/`)으로 떠 놓고 비교했고, 끝에 작업 트리와 `cmp`로 일치를 재확인했습니다.
- 수식 단위 비교: `math-inline`/`math-display` 요소를 HTML 파서로 뽑아 old/new를 1:1 대응시켰습니다. 요소 수는 old와 new가 같습니다(T01 328, T02 809, T03 1045, T04 836, T05 1035).
- 내용이 바뀐 수식은 **54개**입니다. 이 밖에 공백만 바뀐 요소가 15개 있으며(`\quad`), 아래 6번 행에서 다룹니다.
- 변경 파일은 T01–T05 HTML 5개이고, 그 밖의 요소는 byte 동일입니다.
- DOCX는 T02·T05의 `word/document.xml` 텍스트를 HEAD와 비교했습니다.
- 판정: **OK** = 새 렌더링이 LaTeX와 일치하고 old보다 낫거나 같음, **defect** = 새 출력에서도 여전히 틀림(old에도 있던 것), **regression** = old보다 나빠짐.

| # | 트랙 | 파일:위치 (HTML 줄 / 소스 md) | LaTeX 소스 | old | new | 판정 | 제안 |
|---|---|---|---|---|---|---|---|
| 1 | T01 | Multiome HTML L3223 / `docs/07_TF_FOOTPRINTING_AND_CHROMVAR.md:65` | `\text{FPS}=\overline{\text{flank}^{+}}-\overline{\text{core}^{-}}` | `FPS={flank⁺}-{core⁻}` (리터럴 중괄호) | `FPS=` + `flank⁺` 전체와 `core⁻` 전체 위에 overline | OK | overline이 `flank⁺`(위첨자 포함)를 덮어 소스와 일치합니다. 다중문자 `\overline`이지만 `\text` 그룹 하나라 알려진 결함의 영향이 없습니다. |
| 2 | T02 | Spatial L2410, L2518 / `docs/03_SPATIAL_QC_AND_PREPROCESSING.md:85,199` (HTML 2곳, 수식 2개) | `\frac{\overline{\text{Counts}}_{\text{Blank}}}{\overline{\text{Counts}}_{\text{Target}}}` | `{Counts}_Blank / {Counts}_Target` (중괄호 리터럴) | `Counts`에만 overline, 아래첨자는 overline 밖 (`Counts̄_Blank / Counts̄_Target`) | OK | `\overline{X}_{Y}`의 아래첨자는 overline 밖이 맞습니다. DOCX는 `C̅o̅u̅n̅t̅s̅_Blank`(글자별 U+0305)로 같은 의미입니다. |
| 3 | T02 | Spatial L2772 / `docs/04_CELL_SEGMENTATION_ALGORITHMS.md:50` | `\frac{\partial V}{\partial t} = \nabla^2 V` | `partial V / partial t = nabla² V` (명령어 이름 노출) | `∂V / ∂t = ∇²V` | OK | 확산 방정식이 정상으로 읽힙니다. |
| 4 | T02 | Spatial L4485, L4490 (수식 2개) / `docs/08_SPATIAL_LIGAND_RECEPTOR_SIGNALING.md:42,47` | `\frac{\partial C}{\partial t} = D\nabla^2 C - k_{\deg}C` 및 `\frac{\partial C}{\partial t}=0` | `partial C / partial t`, `nabla²` | `∂C/∂t = D∇²C − k_deg C`, `∂C/∂t = 0` | OK | — |
| 5 | T02 | Spatial L2774, L2932 (수식 2개) / `docs/04_CELL_SEGMENTATION_ALGORITHMS.md:54` | `\vec{v}(x,y)=\frac{\nabla V(x,y)}{\|\nabla V(x,y)\|}` | `{v}(x,y) = nabla V / \|nabla V\|` | `{v}(x,y) = ∇V / \|∇V\|` | **defect** (old에도 있음, 미해결) | (a) `\nabla`는 고쳐졌으나 노름 `\|`가 백슬래시 그대로 남습니다. `\|`→`∥`(U+2225)로 매핑해야 합니다. (b) `\vec{v}`가 `{v}`로 남습니다. `\vec`→`v⃗`(v + U+20D7)로 처리해야 합니다. (c) HTML 줄 L2932 본문(`\frac{\nabla V}{|\nabla V|}`)도 같은 (b) 문제입니다. |
| 6 | T02·T01 외 전체 | 공백만 바뀐 15개(예: Multiome L1888, L2192, L2974; Spatial L1616, L3060, L3745, L4678…) | `\quad`, `\;` 등 | 공백 3칸 | 공백 4칸 | OK (실효 없음) | HTML은 연속 공백을 합치므로 `\quad` 간격은 시각적으로 사라집니다. `.math-display`가 `white-space: pre`가 아니라서 old와 화면상 차이는 없습니다. 간격이 필요하면 `&emsp;`(또는 `&nbsp;&nbsp;`)로 출력해야 합니다. |
| 7 | T02 | Spatial L2962 / `docs/04_CELL_SEGMENTATION_ALGORITHMS.md` (Baysor 공분산) | `\boldsymbol{\Sigma}_k` | `{Σ}_k` | **Σ**<sub>k</sub> (굵게) | OK | 소스의 굵은체와 일치합니다. |
| 8 | T02 | Spatial L3260, L3508 ×2, L3793, L3795 (수식 5개) / `docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md:129,350`, `docs/06_SPATIALLY_VARIABLE_GENES.md:92,93` | `\operatorname{Var}(Y_{s,g})`, `\operatorname{Var}(I)` | `{Var}(Y)` | `Var(Y)` (math-op 클래스) | OK | 중괄호가 사라지고 정체로 렌더됩니다. |
| 9 | T02 | Spatial L3264–3266 (수식 3개) / `docs/05_SPATIAL_DECONVOLUTION_MATHEMATICS.md` (Gamma 사전분포) | `\text{Gamma}(a_k, a_k/\mu_{w,k})` 등 | `{Gamma}(a_k, …)` | `Gamma(a_k, a_k/μ_{w,k})` | OK | — |
| 10 | T02 | Spatial L3285, L3491, L3856, L4010 (수식 4개) / `docs/05…:162,334`, `docs/06…:167,312` | `\operatorname{ELBO}`, `\operatorname{Cov}`, `\operatorname{tr}` | `{ELBO}(q)`, `{Cov}(ŵ)`, `{tr}(E_C Σ_C)` | `ELBO(q)`, `Cov(ŵ)`, `tr(E_C Σ_C)/n` | OK | 괄호 안 인자 범위와 분수 `1/n` 뒤 위치가 소스와 같습니다. |
| 11 | T02 | Spatial L3763–3770, L3799, L3820, L3964, L3969 (수식 11개, Moran's I) / `docs/06_SPATIALLY_VARIABLE_GENES.md:62–74,100,124,278,284` | `\bar{x}`, `(x_i-\bar{x})(x_j-\bar{x})`, `z = x-\bar{x}\mathbf{1}` | `{x}`·`(x_i - {x})(x_j - {x})` (**막대가 사라져 평균 x와 구별 불가**) | `x̄`, `(x_i − x̄)(x_j − x̄)` | OK (old보다 크게 개선) | 가장 중요한 수정입니다. old는 `x̄`가 일반 `x`로 보여 수식 의미가 틀렸습니다. DOCX도 `x̄`로 정상입니다. |
| 12 | T03 | TME L2160 / `docs/03_EXPRESSION_THRESHOLDING_AND_PERMUTATIONS.md:90` | `\operatorname{mean}(\bar{E}_L,\bar{E}_R)` | `{mean}(Ē_L,Ē_R)` | `mean(Ē_L,Ē_R)` | OK | — |
| 13 | T03 | TME L2613 / `docs/04_CELLCHAT_MULTISUBUNIT_MODELING.md:66` | `\underbrace{\frac{L_iR_j}{K_h+L_iR_j}}_{\text{Hill}} \times M_{\text{co-receptor}}(RA,RI)…` | `{ L_iR_j/(K_h+L_iR_j) }_Hill` (바깥 중괄호 노출) | `L_iR_j/(K_h+L_iR_j)` 뒤 `_(Hill)` 아래첨자 | OK (경미한 한계) | `\underbrace`가 라벨 아래첨자로 단순화됩니다(위치는 분수 바로 뒤). 의미는 보존됩니다. 필요하면 `\underbrace` 전용 렌더링(`display:inline-block` + 하단 라벨)을 추가할 수 있습니다. |
| 14 | T03 | TME L2613 (같은 수식 끝부분) / `docs/04…:66` | `\;\big[\times \tfrac{n_in_j}{n^2}\big]` | `[× tfrac{n_i n_j}{n^2}]` | **동일** (`tfrac{n_i n_j}{n^2}` 리터럴) | **defect** (old에도 있음) | `\tfrac`/`\dfrac`를 `\frac`와 같이 `math-frac`으로 처리해야 합니다. |
| 15 | T03 | TME L4579 / `docs/08_PROSTATE_BONE_METASTASIS_TME_CIRCUITS.md:107` | `\operatorname{Re}\lambda=-0.04` | `{Re} λ=-0.04` | `Re λ=-0.04` | OK | — |
| 16 | T04 | Cancer L2599, L2600 (수식 2개) / `docs/03_REFERENCE_FREE_COPYKAT.md:145,147` | `\bar{\mathbf{z}}_i`, `\text{median}(\bar{\mathbf{z}}_i)` | `z̄_i`(굵기 소실) | **z**(굵게)에 CSS overline + `_i` | OK | `\bar{\mathbf z}`의 굵은체와 막대가 모두 보존됩니다. 같은 수식의 `\bar{z}_{i1}`은 결합 매크론(z̄)이라 막대 방식이 둘로 섞여 보이는 정도입니다. |
| 17 | T04 | Cancer L4601, L4609, L4610 (수식 3개) / `docs/07_TREATMENT_SELECTION_AND_LINEAGE_PLASTICITY.md:278,283,284` | `\dot{N}_R`, `\dot{N}_S` | `{N}_R` (점 표기 소실, N과 구별 불가) | `Ṅ_R`, `Ṅ_S` (U+1E44) | OK (old보다 크게 개선) | 시간미분 표기가 복원되어 `df_R/dt` 유도가 읽힙니다. |
| 18 | T05 | Generative L2145, L2162, L3427, L4627, L5748 (수식 5개, 6곳) / `docs/02_SCVI…:77,101`, `docs/04_SCARCHES…`, `docs/07_GENEFORMER…` 등 | `\nabla_\phi`, `\nabla_{\theta_{\text{ref}}}`, `\nabla_\theta\mathcal{L}` | `nabla_ϕ …`, `nabla_(θ)` | `∇_ϕ`, `∇_(θ_ref)`, `∇_θ ℒ` | OK | 하첨자 범위(ϕ, θ, θ_ref)가 정확합니다. |
| 19 | T05 | Generative L2176, L2178, L2181, L2340, L2556, L2564, L2740, L2939 (수식 8개) / `docs/02_SCVI…:129–162`, 같은 문서 KL 식들 | `D_{\text{KL}}(q \parallel p)` | `D_KL(q ‖ p)` (U+2016) | `D_KL(q ∥ p)` (U+2225) | OK | 평행 기호가 KL 발산의 `∥`로 정확합니다. DOCX도 같습니다. |
| 20 | T05 | Generative L2176 / `docs/02_SCVI_VARIATIONAL_AUTOENCODER.md:129` | `\underbrace{\mathbb E_q[\log p]}_{\text{재구성 우도 (Reconstruction Accuracy)}} - \underbrace{D_{KL}(\cdot\parallel\cdot)}_{\text{사전분포 정규화 (Latent Regularization)}}` | `{𝔼…}_(재구성…) - {D_KL(…)}_(사전분포…)` (바깥 중괄호 노출) | 중괄호 제거, 라벨은 `]`·`)` 뒤 아래첨자 | OK (경미한 한계) | 13번과 같은 `\underbrace` 단순화입니다. 라벨이 아래첨자로 붙어 읽을 수는 있습니다. |
| 21 | T05 | Generative L2586 / `docs/02_SCVI…:418` | `\mathrm{KL}(q(z\mid x)\,\|\,\mathcal{N}(0,I))` | `{KL}(q(z\| x) \| 𝒩(0,I))` → 렌더에서 단일 막대 `|` | `KL(q(z\| x) \| 𝒩(0,I))` → 단일 막대 `|` | **defect** (old에도 있음) | 중괄호는 해결됐으나 `\,\|\,`가 조건부 확률 `|`처럼 보여 KL 발산이 아닌 것처럼 읽힙니다(같은 문서의 `\parallel`은 `∥`). `\|`→`∥`로 통일해야 합니다. DOCX에는 `\|`가 백슬래시 포함으로 남습니다. |
| 22 | T05 | Generative L5049, L5053, L5056 (수식 3개) / `docs/07_GENEFORMER_AND_GENOME_ATTENTION.md:336–339` | `\boldsymbol{\pi}^T`, `\boldsymbol{\pi}` | `{π}^T`, `{π}` | **π**(굵게)`^T` | OK | 벡터 굵은체 표기가 소스와 같습니다. |
| 23 | T01 | Multiome L2129 · L2461 · L2467 (수식 3개, 검토 중 확인) / `docs/03_SNATAC_QC_METRICS.md:246,248`, `docs/04…:234` | `\text{n\_fragments}` | `n_(f)ragments` | **동일** (`n_fragments`가 `n`+아래첨자 `f`+`ragments`로 깨짐) | **defect** (old에도 있음) | `\text{…}` 안의 `\_`를 리터럴 `_`로 처리해야 합니다. 같은 부류가 T02에도 있습니다(Spatial `total_counts`, `n_genes`, `N_{cells\_per\_location}`; 수식 4개). |

## 요약

| 구분 | 수식 요소 수 |
|---|---|
| 내용이 바뀐 수식 | 54 (T01 1, T02 28, T03 3, T04 5, T05 17) |
| 　OK (개선 또는 동등) | 50 |
| 　defect (old에도 있고 이번에도 미해결, 일부는 부분 개선) | 4: 5번 2개 요소(Spatial L2774, L2932: `\|`·`\vec`), 14번 1개(`\tfrac`), 21번 1개(KL 단일 막대) |
| 　**regression** (old보다 악화) | **0** |
| 공백만 변경 (6번) | 15, 시각 영향 없음 |

- 판정 합계는 50 + 4 = 54입니다.
- 23번(`\_` in `\text`)은 이번 diff 대상이 아니라 새 출력 전수 검사에서 나온 동류 결함이라 54개에 넣지 않았습니다.
- 새 출력(T01–T05) 전체에서 남은 이스케이프/리터럴 잔재는 다음 4부류뿐입니다. 이 밖에 중괄호는 `\{…\}` 집합 표기에서 의도된 것만 남아 있습니다.
  - `\|`: T02 8, T03 2, T04 4, T05 16개 요소
  - `\vec{v}`→`{v}`: T02 12개 요소
  - `\tfrac`: T03 1개 요소
  - `\text{…\_…}`: T01 3, T02 4개 요소

### 제안 우선순위 (모두 `scripts/course_markdown.py` LaTeX 단계)
1. `\|`, `\Vert` → `∥`. 영향이 가장 커서 T02–T05에서 약 30개 요소입니다. 노름과 KL 발산 표기가 달라집니다.
2. `\vec{x}` → `x⃗`(결합 화살표 U+20D7)입니다. `\bar`·`\dot`는 이미 처리되므로 같은 경로에 추가하면 됩니다.
3. `\tfrac`, `\dfrac`를 `\frac`와 동일하게 처리합니다.
4. `\text{…}` 안의 `\_`를 리터럴 `_`로 처리합니다.
5. (선택) `\underbrace{A}_{B}`를 하단 라벨 렌더링으로 처리합니다.
6. (선택) `\quad`는 HTML에서 간격이 합쳐지므로 `&emsp;`를 씁니다.

### 부수 관찰 (검토 범위 밖)
- 검토 중 `tracks/06_manifold_trajectory/Manifold_Trajectory_Master_Course.html`도 수정됨으로 보였습니다(`git status`). 이번 지시 범위는 T01–T05라 검토하지 않았습니다.
- 트랙 스위처 마커(`TRACK-SWITCHER`)는 T01–T05 모두 old/new 6개로 동일합니다(초기 diff에서 보였던 CSS 줄 이동은 줄 정렬 잡음입니다).
- 다중문자 `\overline`은 지시대로 판단에서 제외했습니다. 다만 위 1·2번은 `\text` 그룹 하나를 덮는 경우라 올바르게 렌더되고 있습니다.

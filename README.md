# IP-C-PDL

**IP-C-PDL: A chaos-embedded trust-region interior-point framework for constrained economic dispatch in large-scale power systems**

Raharjo, Jangkung · Ikhsan, Rifki Rahman Nur · Yustika, Lindiasari Martha  
*IJOCTA, 2026*  
[![DOI](https://img.shields.io/badge/DOI-10.36922%2FIJOCTA026310188-blue)](https://doi.org/10.36922/IJOCTA026310188)

Chaos-embedded trust-region interior-point framework (Powell's Dog-Leg) for constrained economic dispatch in large-scale power systems.

## Abstract

Constrained economic dispatch (ED)—the classical, single-bus, active-poweronly dispatch problem for thermal generation, not an optimal power flow or unit-commitment problem—in large-scale interconnected power systems such as the Java–Madura–Bali (Jamali) network is non-convex and ramprate-coupled: deterministic solvers stall at local optima, while stochastic solvers incur unacceptable run-to-run variance. This paper proposes IP-C-PDL (Interior-Point Chaotic Powell’s Dog-Leg), which, within a single interior-point iteration, couples a Powell’s Dog-Leg trust-region step (Cauchy/quasi-Newton selection via accumulated SR1/BFGS curvature) with a stall-adaptive tentmap chaotic perturbation for ergodic basin escape; a multi-start pool, elite archive, and periodic SQP refinement complete the framework using only three user-specified parameters. On twelve CEC benchmark functions (30 runs, 12 methods), IP-C-PDL ranks 1 in 1 of 12 and top-3 in 7 of 12. On six ED scenarios spanning 5-, 16-, and 42-unit systems (50 runs, 10 methods, single-period and 24-hour horizons), it ranks 1 in 4 of 6 (mean rank 2.17/10) with zero constraint violations, yielding USD 22,239.29/day (σ = USD 0.41) on the 42-unit 24-hour case—0.16% better than the second-ranked comparator (t = −48.26, p < 10−42) and 4.67% better than the worst (t = −180.85, p < 10−70)—confirmed by Friedman, ANOVA, and Wilcoxon–Holm–Bonferroni tests. A five-factor sensitivity analysis identifies multi-start pool size as the dominant factor—with a cliff-edge improvement at just two starts—while the framework is insensitive to the remaining four hyperparameters and stabilises by 400 iterations. A policy sensitivity analysis further links dispatch quality to fiscal savings, Biaya Pokok Penyediaan (BPP) subsidy reduction, and avoided carbon emissions in the Indonesian power sector.


## Proposed method

IP-C-PDL is a chaos-embedded trust-region interior-point method. A multi-start pool (lambda iteration, interior point, SQP and repaired random points) supplies the starting dispatch; each start is refined by an interior-point penalty loop whose step comes from Powell's dog-leg (PDL), built on an SR1/BFGS Hessian. A tent-map perturbation, whose probability grows with the stall counter, escapes flat or misleading basins. The 24-hour problem is solved as a sequence of ramp-coupled single-period problems.

### Dispatch model

$$\min_{\mathbf{P}}F(\mathbf{P})=\sum_{i=1}^{n}\left(a_iP_i^2+b_iP_i+c_i\right)$$

$$\text{s.t.}\quad\mathbf{1}^{\top}\mathbf{P}=P_D,\qquad \mathbf{P}^{\min}\le\mathbf{P}\le\mathbf{P}^{\max}.$$

For 24 hours the objective is summed over periods and the ramp limits couple consecutive hours:

$$\min_{\lbrace P_{i,t}\rbrace}\sum_{t=1}^{T}\sum_{i=1}^{n}f_i(P_{i,t}),\qquad -RD_i\le P_{i,t}-P_{i,t-1}\le RU_i.$$

The ramp limits are folded into effective bounds:

$$\tilde P_i^{\min}(t)=\max\left(P_i^{\min},P_{i,t-1}-RD_i\right),\qquad \tilde P_i^{\max}(t)=\min\left(P_i^{\max},P_{i,t-1}+RU_i\right).$$

### Notation

| Symbol | Meaning | Unit |
|---|---|---|
| $n$ | number of units | – |
| $T$ | number of periods (24) | – |
| $P_i,P_{i,t}$ | unit output | MW |
| $a_i,b_i,c_i$ | cost coefficients | IDR/MW²h, IDR/MWh, IDR/h |
| $P_D$ | demand | MW |
| $RU_i,RD_i$ | ramp-up and ramp-down limits | MW/h |
| $\Phi$ | penalised objective | IDR/h |
| $\mu$ | penalty parameter | – |
| $\varepsilon_k$ | balance error at iteration $k$ | MW |
| $\mathbf{H}_k$ | Hessian approximation | – |
| $\mathbf{s}_k,\mathbf{y}_k$ | step and gradient difference | MW, IDR/MWh |
| $\mathbf{p}_N,\mathbf{p}_C$ | Newton and Cauchy steps | MW |
| $\Delta_k$ | trust-region radius | MW |
| $\rho_k$ | predicted/actual reduction ratio | – |
| $\alpha_k$ | Armijo step length | – |
| $c_p$ | tent-map perturbation probability | – |
| $\mathrm{stall}$ | consecutive non-improving iterations | – |

### Penalty subproblem

The equality constraint is moved into the objective:

$$\Phi(\mathbf{P};\mu)=F(\mathbf{P})+\mu\left(\sum_{i=1}^{n}P_i-P_D\right)^2.$$

The penalty is tightened as the balance error shrinks:

$$\mu_k=\min\left(\mu_0\max\left(1,\frac{10^3}{\max(\varepsilon_k,10^{-3})}\right),10^{10}\right),\qquad \mu_0=10^6,\quad \varepsilon_k=\left|\mathbf{1}^{\top}\mathbf{P}_k-P_D\right|.$$

A small balance error therefore receives a large penalty, while the cap at $10^{10}$ prevents the Hessian from becoming numerically singular.

Exact balance repair is applied after clipping. With $d=\mathbf{1}^{\top}\mathbf{P}-P_D$:

$$P_i\leftarrow
\begin{cases}
P_i-\dfrac{P_i-l_i}{\sum_j(P_j-l_j)}d, & d\gt0,\\
P_i+\dfrac{u_i-P_i}{\sum_j(u_j-P_j)}|d|, & d\lt0.
\end{cases}$$

The correction is proportional to the remaining headroom, so units near a limit move less. The repair is iterated up to 200 times until $|d|\lt10^{-4}$ MW, then re-clipped.

### Hessian approximation

With $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$ and $\mathbf{y}_k=\nabla\Phi_{k+1}-\nabla\Phi_k$, the BFGS update is used when $\mathbf{s}_k^{\top}\mathbf{y}_k\gt10^{-10}$:

$$\mathbf{H}_{k+1}=\left(\mathbf{I}-\rho_k\mathbf{s}_k\mathbf{y}_k^{\top}\right)\mathbf{H}_k\left(\mathbf{I}-\rho_k\mathbf{y}_k\mathbf{s}_k^{\top}\right)+\rho_k\mathbf{s}_k\mathbf{s}_k^{\top},\qquad \rho_k=\frac{1}{\mathbf{s}_k^{\top}\mathbf{y}_k}.$$

Otherwise the SR1 fallback is used:

$$\mathbf{H}_{k+1}=\mathbf{H}_k+\frac{\mathbf{r}_k\mathbf{r}_k^{\top}}{\mathbf{r}_k^{\top}\mathbf{y}_k},\qquad \mathbf{r}_k=\mathbf{s}_k-\mathbf{H}_k\mathbf{y}_k.$$

An update is rejected, and $\mathbf{H}_{k+1}\leftarrow\mathbf{I}$, when $\mathrm{rcond}(\mathbf{H}_{k+1})\le10^{-14}$.

### Powell's dog-leg step

The gradient is scaled so units with different ranges contribute comparably:

$$\boldsymbol{\sigma}=\max(\mathbf{u}-\mathbf{l},\mathbf{1}),\qquad \mathbf{g}_s=\frac{\mathbf{g}_k/\boldsymbol{\sigma}}{\max\left(\Vert\mathbf{g}_k\Vert/\Vert\boldsymbol{\sigma}\Vert,10^{-12}\right)}.$$

The quasi-Newton and Cauchy steps are

$$\mathbf{p}_N=-\mathbf{H}_k\mathbf{g}_s,\qquad \mathbf{p}_C=-\tau_c\frac{\Delta_k}{\Vert\mathbf{g}_s\Vert}\mathbf{g}_s,\qquad \tau_c=\min\left(\frac{\Vert\mathbf{g}_s\Vert^3}{\Delta_k~\mathbf{g}_s^{\top}\mathbf{H}_k\mathbf{g}_s+10^{-12}},1\right).$$

The trial step selects among three cases:

$$\mathbf{p}_k=
\begin{cases}
\mathbf{p}_N, & \Vert\mathbf{p}_N\Vert\le\Delta_k,\\
\Delta_k\dfrac{-\mathbf{g}_s}{\Vert\mathbf{g}_s\Vert}, & \Vert\mathbf{p}_C\Vert\ge\Delta_k,\\
\mathbf{p}_C+\tau^{\ast}(\mathbf{p}_N-\mathbf{p}_C), & \text{otherwise},
\end{cases}$$

where $\tau^{\ast}$ is the positive root of

$$\underbrace{\Vert\mathbf{p}_N-\mathbf{p}_C\Vert^2}_{a_q}\tau^2+\underbrace{2\mathbf{p}_C^{\top}(\mathbf{p}_N-\mathbf{p}_C)}_{b_q}\tau+\underbrace{\Vert\mathbf{p}_C\Vert^2-\Delta_k^2}_{c_q}=0.$$

The first case takes the full Newton step inside the trust region; the second follows the steepest-descent direction to the boundary when even the Cauchy step is too long; the third interpolates along the dog-leg segment to the trust-region boundary.

An Armijo line search finds $\alpha_k$ with $10^{-4}$ sufficient decrease:

$$\Phi(\mathbf{x}_k+\alpha_k\mathbf{p}_k)\le\Phi(\mathbf{x}_k)+10^{-4}\alpha_k\mathbf{g}_s^{\top}\mathbf{p}_k.$$

The trust-region ratio and radius update are

$$\rho_k=\frac{\Phi(\mathbf{x}_k)-\Phi(\mathbf{x}_k+\alpha_k\mathbf{p}_k)}{-\mathbf{g}_s^{\top}(\alpha_k\mathbf{p}_k)},\qquad
\Delta_{k+1}=
\begin{cases}
\max(\Delta_k/2,~10^{-8}), & \rho_k\lt0.25,\\
\min(2\Delta_k,~500), & \rho_k\gt0.75,\\
\Delta_k, & \text{otherwise}.
\end{cases}$$

A ratio below 0.25 means the model over-predicted, so the radius shrinks; a ratio above 0.75 means the model is trustworthy, so the radius grows.

### Tent-map perturbation

$$s_i^{(k+1)}=
\begin{cases}
2s_i^{(k)}, & s_i^{(k)}\lt0.5,\\
2\left(1-s_i^{(k)}\right), & s_i^{(k)}\ge0.5,
\end{cases}
\qquad
c_p(k)=\min\left(c_{p0}+0.004~\mathrm{stall},~0.55\right).$$

$$\Delta_c=\Delta_k\min\left(0.25+\frac{\mathrm{stall}}{120},0.9\right),\qquad \mathbf{x}_k\leftarrow\mathbf{x}_k+\Delta_c\frac{\mathbf{s}^{(k)}-0.5}{\left\Vert\mathbf{s}^{(k)}-0.5\right\Vert}~\text{with probability }c_p(k).$$

The tent map is deterministic and aperiodic, so the perturbation is reproducible. Its probability and magnitude rise with the stall counter, which lets the search leave a basin only after the deterministic step has stopped improving.

An elite archive of six solutions is kept; every $e_i=\max(20-\lfloor\mathrm{stall}/10\rfloor,10)$ stalled iterations an SQP refinement is applied, and every 50 stalls the run restarts from the best elite member plus Gaussian noise.

### Reported result

42-unit single-period: USD 22,239.29/day. The same framework covers 5-, 16- and 42-unit single-period and 24-hour cases, plus CEC2020 benchmarks.

## Contents

```
.gitignore
figures/
  Fig01_CEC_Conv.png
  Fig02_Conv_5Unit.png
  Fig03_Conv_16Unit.png
  Fig04_Conv_42Unit.png
  Fig05_Conv_5Unit_24h.png
  Fig06_Conv_16Unit_24h.png
  Fig07_Conv_42Unit_24h.png
  Fig08_Sensitivity.png
matlab/
  IP_C_PDL.m
python/
  BFGS (ED).ipynb
  Benchmark_IP-C-PDL.ipynb
  CEC2020.ipynb
  CG (ED).ipynb
  COBYLA (ED).ipynb
  IP-C-PDL v2.0.ipynb
  L-BFGS-B (ED).ipynb
  Nelder-Mead (ED).ipynb
  Powell (ED).ipynb
  SLSQP (ED).ipynb
  TNC (ED).ipynb
  trust-constr (ED).ipynb
results/
  CEC_Results_v9.xlsx
  Constraint_Check_v9.xlsx
  ED_Results_v9.xlsx
```

## Data

Input data is **not included**. Place files in `data/` (notebooks) or next to the `.m` script (MATLAB):

| File | Columns |
|---|---|
| `input16u.xlsx` | Sheet1: `Unit`, `a`, `b`, `c`, `pmin`, `pmax`, `rampup`, `rampdown` |
| `input42u.xlsx` | Sheet1: `Unit`, `a`, `b`, `c`, `pmin`, `pmax`, `rampup`, `rampdown`, `load` |
| `input5u.xlsx` | Sheet1: `Unit`, `a`, `b`, `c`, `pmin`, `pmax`, `rampup`, `rampdown` |
| `power_demand16u.xlsx` | Sheet1: `load` |
| `power_demand42u.xlsx` | Sheet1: `load` |
| `power_demand5u.xlsx` | Sheet1: `load` |


## Results

### 5-Unit case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| Powell | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.37 | 0.00 |
| L-BFGS-B | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.36 | 0.00 |
| TNC | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| SLSQP | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| COBYLA | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| trust-constr | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| PSO | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| DE | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |
| IP-C-PDL | 50 | 1,733,364.35 | 1,733,364.35 | 1,733,364.35 | 0.00 |

All ten methods land on 1,733,364 IDR with spreads of a few IDR, so the five-unit case separates nothing on cost. Ranking is decided by run time and consistency, and IP-C-PDL sits sixth here.

### 16-Unit case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 7,343,363.27 | 7,281,425.95 | 7,407,945.74 | 32,353.27 |
| Powell | 50 | 7,350,988.41 | 7,273,593.71 | 7,404,398.28 | 34,010.90 |
| L-BFGS-B | 50 | 7,266,997.20 | 7,266,972.86 | 7,268,163.13 | 168.27 |
| TNC | 50 | 7,266,972.86 | 7,266,972.86 | 7,266,972.86 | 0.00 |
| SLSQP | 50 | 7,266,972.86 | 7,266,972.86 | 7,266,972.86 | 0.00 |
| COBYLA | 50 | 7,359,196.71 | 7,285,446.24 | 7,434,825.14 | 41,020.00 |
| trust-constr | 50 | 7,288,136.91 | 7,266,982.71 | 7,325,402.28 | 14,389.99 |
| PSO | 50 | 7,330,112.51 | 7,275,750.61 | 7,441,832.11 | 41,458.67 |
| DE | 50 | 7,296,347.71 | 7,272,955.48 | 7,350,669.94 | 16,372.84 |
| IP-C-PDL | 50 | 7,266,972.86 | 7,266,972.86 | 7,266,972.86 | 0.00 |

IP-C-PDL ties the best mean at 7,266,973 IDR with SLSQP and TNC and leads L-BFGS-B by 24 IDR. A gap this small is solver tolerance, not algorithm quality.

### 42-Unit case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 15,567,152.11 | 15,406,121.21 | 15,741,086.99 | 85,092.82 |
| Powell | 50 | 15,572,499.38 | 15,403,824.95 | 15,753,800.92 | 78,777.23 |
| L-BFGS-B | 50 | 15,123,804.22 | 15,042,644.08 | 15,237,694.22 | 54,259.07 |
| TNC | 50 | 15,052,845.41 | 15,027,878.59 | 15,106,586.81 | 25,959.41 |
| SLSQP | 50 | 15,052,845.41 | 15,027,878.59 | 15,106,586.81 | 25,959.41 |
| COBYLA | 50 | 15,552,848.48 | 15,364,182.34 | 15,722,765.69 | 86,419.81 |
| trust-constr | 50 | 15,566,319.17 | 15,116,047.30 | 15,839,042.84 | 164,644.59 |
| PSO | 50 | 15,832,034.14 | 15,506,582.24 | 16,042,219.11 | 118,093.65 |
| DE | 50 | 15,686,226.01 | 15,405,843.91 | 15,824,873.60 | 95,795.03 |
| IP-C-PDL | 50 | 15,027,878.55 | 15,027,878.54 | 15,027,878.60 | 0.02 |

On the 42-unit case IP-C-PDL is lowest at 15,027,879 IDR, 0.17% under SLSQP and TNC and 0.63% under L-BFGS-B. The advantage comes from projecting the step onto the coupling constraints, which matters once the unit count grows.

### 5-Unit 24h case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 42,544,496.21 | 42,459,856.53 | 42,673,627.79 | 45,190.75 |
| Powell | 50 | 42,559,295.68 | 42,451,224.45 | 42,647,673.44 | 44,964.32 |
| L-BFGS-B | 50 | 42,291,364.30 | 42,285,718.83 | 42,300,846.49 | 4,829.36 |
| TNC | 50 | 42,290,012.70 | 42,285,718.87 | 42,300,791.47 | 4,620.07 |
| SLSQP | 50 | 42,290,137.69 | 42,285,718.87 | 42,300,791.47 | 4,585.94 |
| COBYLA | 50 | 42,508,107.21 | 42,403,310.28 | 42,605,471.38 | 48,538.76 |
| trust-constr | 50 | 42,545,235.38 | 42,457,706.14 | 42,663,667.83 | 46,307.93 |
| PSO | 50 | 42,761,604.51 | 42,572,065.52 | 42,987,097.51 | 86,834.79 |
| DE | 50 | 42,306,845.66 | 42,297,048.11 | 42,336,227.28 | 6,197.99 |
| IP-C-PDL | 50 | 42,286,218.45 | 42,285,718.50 | 42,291,967.84 | 1,712.62 |

Over 24 hours IP-C-PDL averages 42,286,218 IDR per day against 42,290,013 IDR for TNC, a 0.01% margin. The daily five-unit problem is easy and every deterministic solver ends within a few thousand IDR.

### 16-Unit 24h case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 182,245,569.31 | 181,987,913.94 | 182,589,675.27 | 146,181.37 |
| Powell | 50 | 182,527,007.82 | 182,198,084.61 | 182,930,382.70 | 159,250.72 |
| L-BFGS-B | 50 | 180,187,349.94 | 180,152,760.89 | 180,252,378.89 | 26,561.38 |
| TNC | 50 | 180,224,102.29 | 180,182,097.25 | 180,300,892.97 | 30,743.79 |
| SLSQP | 50 | 180,224,102.29 | 180,182,097.25 | 180,300,892.97 | 30,743.79 |
| COBYLA | 50 | 182,129,937.68 | 181,851,346.99 | 182,363,075.04 | 139,733.86 |
| trust-constr | 50 | 181,473,947.92 | 180,880,486.12 | 182,199,976.73 | 292,101.14 |
| PSO | 50 | 183,773,892.95 | 183,284,491.44 | 184,322,828.74 | 256,042.79 |
| DE | 50 | 182,612,292.47 | 182,243,978.32 | 182,907,063.02 | 140,723.67 |
| IP-C-PDL | 50 | 180,184,023.83 | 180,152,758.80 | 180,243,875.08 | 28,002.10 |

IP-C-PDL averages 180,184,024 IDR, 3,326 IDR ahead of L-BFGS-B and 40,078 IDR ahead of SLSQP. The margin is stable across the 50 starts and the standard deviation is the smallest in the table apart from L-BFGS-B.

### 42-Unit 24h case, 50 runs

| Method | Runs | Mean cost (IDR) | Best | Worst | Std |
|---|---|---|---|---|---|
| Nelder-Mead | 50 | 369,155,115.56 | 368,460,223.50 | 369,790,195.56 | 291,859.28 |
| Powell | 50 | 370,690,679.47 | 370,037,300.85 | 371,385,784.85 | 343,734.33 |
| L-BFGS-B | 50 | 359,524,399.90 | 359,198,030.61 | 359,835,904.44 | 146,998.13 |
| TNC | 50 | 358,633,959.84 | 358,441,244.36 | 358,834,376.47 | 83,745.64 |
| SLSQP | 50 | 358,633,959.84 | 358,441,244.36 | 358,834,376.47 | 83,745.64 |
| COBYLA | 50 | 368,751,105.12 | 368,042,890.06 | 369,384,410.07 | 292,507.38 |
| trust-constr | 50 | 367,178,284.14 | 365,395,531.00 | 369,666,197.51 | 905,569.48 |
| PSO | 50 | 374,787,592.89 | 373,570,526.38 | 376,095,889.53 | 654,805.29 |
| DE | 50 | 371,845,272.62 | 370,859,359.47 | 372,642,311.56 | 374,202.81 |
| IP-C-PDL | 50 | 358,052,521.64 | 358,041,500.53 | 358,075,181.43 | 6,558.02 |

On the 42-unit daily case the gap grows to 581,438 IDR, 0.16% under SLSQP and TNC. This is the scenario the paper uses to claim the method, since it combines the largest unit count with the ramp coupling over 24 hours.

### Result workbooks

| File | Content |
|---|---|
| [`results/CEC_Results_v9.xlsx`](results/CEC_Results_v9.xlsx) | raw output |
| [`results/Constraint_Check_v9.xlsx`](results/Constraint_Check_v9.xlsx) | raw output |
| [`results/ED_Results_v9.xlsx`](results/ED_Results_v9.xlsx) | raw output |

The workbooks hold the raw solver output behind the tables above: per-unit schedules, per-hour costs, constraint checks and run statistics.


## Requirements

MATLAB, Python. Python notebooks: `numpy`, `pandas`, `scipy`, `matplotlib` (as imported). MATLAB scripts: R2021a+.

## Citation

```bibtex
@article{Raharjo_2026, title={IP-C-PDL: A chaos-embedded trust-region interior-point framework for constrained economic dispatch in large-scale power systems}, volume={0}, ISSN={2146-5703}, url={http://dx.doi.org/10.36922/IJOCTA026310188}, DOI={10.36922/ijocta026310188}, number={0}, journal={An International Journal of Optimization and Control: Theories & Applications  (IJOCTA)}, publisher={AccScience Publishing}, author={Raharjo, Jangkung and Ikhsan, Rifki Rahman Nur and Yustika, Lindiasari Martha}, year={2026}, month=Nov, pages={026310188} }
```

## License

[![CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](LICENSE)

Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). **Attribution is required**: any use of this code, results, or figures, including in publications, must cite:

> Jangkung Raharjo; Rifki Rahman Nur Ikhsan; Lindiasari Martha Yustika (2026). IP-C-PDL: A chaos-embedded trust-region interior-point framework for constrained economic dispatch in large-scale power systems. *An International Journal of Optimization and Control: Theories & Applications  (IJOCTA)*. https://doi.org/10.36922/IJOCTA026310188

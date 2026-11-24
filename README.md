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

$$\min_{\lbraceP_{i,t}\rbrace}\sum_{t=1}^{T}\sum_{i=1}^{n}f_i(P_{i,t}),\qquad -RD_i\le P_{i,t}-P_{i,t-1}\le RU_i.$$

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

An update is rejected, and $\mathbf{H}_{k+1}\leftarrow\mathbf{I}$, when $\operatorname{rcond}(\mathbf{H}_{k+1})\le10^{-14}$.

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

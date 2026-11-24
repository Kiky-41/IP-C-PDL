# IP-C-PDL

**IP-C-PDL: A chaos-embedded trust-region interior-point framework for constrained economic dispatch in large-scale power systems**

Raharjo, Jangkung · Ikhsan, Rifki Rahman Nur · Yustika, Lindiasari Martha  
*IJOCTA, 2026*  
[![DOI](https://img.shields.io/badge/DOI-10.36922%2FIJOCTA026310188-blue)](https://doi.org/10.36922/IJOCTA026310188)

Chaos-embedded trust-region interior-point framework (Powell's Dog-Leg) for constrained economic dispatch in large-scale power systems.

## Abstract

Constrained economic dispatch (ED)—the classical, single-bus, active-poweronly dispatch problem for thermal generation, not an optimal power flow or unit-commitment problem—in large-scale interconnected power systems such as the Java–Madura–Bali (Jamali) network is non-convex and ramprate-coupled: deterministic solvers stall at local optima, while stochastic solvers incur unacceptable run-to-run variance. This paper proposes IP-C-PDL (Interior-Point Chaotic Powell’s Dog-Leg), which, within a single interior-point iteration, couples a Powell’s Dog-Leg trust-region step (Cauchy/quasi-Newton selection via accumulated SR1/BFGS curvature) with a stall-adaptive tentmap chaotic perturbation for ergodic basin escape; a multi-start pool, elite archive, and periodic SQP refinement complete the framework using only three user-specified parameters. On twelve CEC benchmark functions (30 runs, 12 methods), IP-C-PDL ranks 1 in 1 of 12 and top-3 in 7 of 12. On six ED scenarios spanning 5-, 16-, and 42-unit systems (50 runs, 10 methods, single-period and 24-hour horizons), it ranks 1 in 4 of 6 (mean rank 2.17/10) with zero constraint violations, yielding USD 22,239.29/day (σ = USD 0.41) on the 42-unit 24-hour case—0.16% better than the second-ranked comparator (t = −48.26, p < 10−42) and 4.67% better than the worst (t = −180.85, p < 10−70)—confirmed by Friedman, ANOVA, and Wilcoxon–Holm–Bonferroni tests. A five-factor sensitivity analysis identifies multi-start pool size as the dominant factor—with a cliff-edge improvement at just two starts—while the framework is insensitive to the remaining four hyperparameters and stabilises by 400 iterations. A policy sensitivity analysis further links dispatch quality to fiscal savings, Biaya Pokok Penyediaan (BPP) subsidy reduction, and avoided carbon emissions in the Indonesian power sector.

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

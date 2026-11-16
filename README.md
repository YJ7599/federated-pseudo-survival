# Federated Survival Analysis with Site-Level Heterogeneity Adjustment

Code to reproduce the simulation studies in

> Hyojung Jang, Malcolm Risk, Yaojie Wang, Norrina Bai Allen, Xu Shi and Lili Zhao.
> *Federated Survival Analysis with Site-Level Heterogeneity Adjustment.*

## Contents

- `R/`: implementation of the proposed method (federated Kaplan–Meier splines and pseudo-observations,
  renewable GEE with sandwich variance, GSURE-based heterogeneity adjustment) and the ODAC comparator
- `simulation/`: simulation scripts for Scenarios 1–3 (Section 3)
- `visualization/`: scripts that produce Figures 1–4 from the simulation results

## Requirements

R (≥ 4.3). Install the required packages with

```bash
Rscript install_packages.R
```

## Reproducing the results

| Figure | Simulation | Plot |
|---|---|---|
| 1, 2 | `simulation/scenario1_ph/run_ph_*.R` (four settings) | `visualization/figure1_ph_bias.R`, `visualization/figure2_ph_se_calibration.R` |
| 3 | `simulation/scenario2_time_varying/run_time_varying.R` | `visualization/figure3_time_varying_effect.R` |
| 4 | `simulation/scenario3_sparse_heterogeneity/run_sparse_heterogeneity.R` | `visualization/figure4_heterogeneity_rmse.R` |

Run the scripts from the repository root, for example

```bash
FEDSURV_NCORE=50 Rscript simulation/scenario2_time_varying/run_time_varying.R
Rscript visualization/figure3_time_varying_effect.R
```

Results are written to `results/` and figures to `figures/`. The number of cores is set with
`FEDSURV_NCORE` (default 50). The full simulations are computationally intensive; with 50 cores,
each Scenario 1 setting takes about 1–1.5 hours and Scenario 2 about 25 minutes. Scenario 3 can be run
in parallel by (π, τ) cell, e.g. `FEDSURV_PI=0.1 FEDSURV_TAU=0.3 Rscript simulation/scenario3_sparse_heterogeneity/run_sparse_heterogeneity.R`.
All simulations are seeded, so re-running them reproduces the reported results.

## License

`R/odac.R` contains functions from the [pda](https://github.com/Penncil/pda) package
(Apache License 2.0) with modifications.

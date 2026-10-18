# Simulation study (Section 4)

Scripts for the Monte Carlo study of the $\mathrm{MZMCMP}_\mu$ regression model in
*Zero-Modified Multivariate Conway-Maxwell-Poisson Regression with Baseline-Mean
Parameterization* (Santana, Conceição, Diniz, Andrade).

## Design

* Trivariate model with one covariate, `log(mu_ij) = beta_i0 + beta_i1 x_j`,
  constant hurdle probability `omega_0 = P(Y != 0)` and constant dispersion
  `kappa_i = log(nu_i)`; `x_j ~ U(0, 1)`, generated once for each `n` and kept fixed.
* Scenarios (true values in `sim_01_cenarios.R`):
  * **S1** zero inflation, overdispersed margins (`p_j` about 0.70);
  * **S2** zero deflation, overdispersed margins, small baseline means (`p_j` from 1.02 to 1.32);
  * **S3** zero inflation, underdispersed margins (`p_j` about 0.80, `nu_i` about 2).
* `n` in {100, 200, 300}; `N = 100` replications; 95% parametric bootstrap
  percentile intervals with `B = 500`.
* Measures: Mean, MAPE, MSE, interval width, empirical coverage (CP95) and
  convergence rates.

## Requirements

* `MultDiscreteZM` >= 0.2.1 (exact baseline mean and analytic gradient)
* `parallel` (base R)

## How to run

```r
# copy the scripts to a writable folder
file.copy(system.file("simulation", package = "MultDiscreteZM"), "C:/my/folder", recursive = TRUE)
setwd("C:/my/folder/simulation")

source("sim_01_cenarios.R")      # true values; scenario figures Fig_cenario_S1-S3.pdf

# quick test: S1, n = 100, 4 replications, B = 20
Sys.setenv(MZM_QUICK = 1, MZM_NSIM = 4, MZM_BSIM = 20)
source("sim_02_monte_carlo.R"); source("sim_03_medidas.R")

# full study (delete results/ from the test first)
Sys.unsetenv("MZM_QUICK"); Sys.setenv(MZM_NSIM = 100, MZM_BSIM = 500)
source("sim_02_monte_carlo.R"); source("sim_03_medidas.R")   # tables of results
```

Results are written to `./results` (or to `MZM_SIM_OUT`).

| Script | Content |
|---|---|
| `sim_00_setup.R` | packages, paths, seed (2026) |
| `utils_sim.R` | model quantities (`p_j`, GDI, correlations), data generation, refits |
| `sim_01_cenarios.R` | scenarios S1–S3; behaviour of `p_j`, GDI and `rho_ij` along `x` |
| `sim_02_monte_carlo.R` | replications, MLE and bootstrap intervals; records convergence |
| `sim_03_medidas.R` | performance measures and LaTeX tables |

### Options of `sim_02_monte_carlo.R` (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `MZM_NSIM` | 100 | replications per scenario and sample size |
| `MZM_BSIM` | 500 | bootstrap samples per replication |
| `MZM_SIZES` | `"100,200,300"` | sample sizes |
| `MZM_NCORES` | cores − 1 | parallel workers |
| `MZM_JOBS` | all | subset of combinations, e.g. `"S1:300,S2:300"`, to split the work among computers |
| `MZM_QUICK` | unset | `1` runs only S1 with `n = 100` |

Each replication is a separate task saved in `results/parts/` as soon as it
ends, so an interrupted run resumes from the replications not yet saved.
Time limits: 1800 s for the main fit and 120 s for each bootstrap refit; a
refit that exceeds the limit or does not converge is recorded as missing.

## Reproducibility

The covariate uses `set.seed(2026 + n)`. Each replication has its own
L'Ecuyer-CMRG stream, derived from `set.seed(2026 + 1000 * s + n)`
(`s` = 1, 2, 3 for S1–S3), so the results do not depend on the number of
cores or on how the combinations are split among computers.

## Results reported in the paper

`paper_results/` contains the results of the full run:

* `sim_results.csv`: Mean, MAPE (%), MSE, Width and CP95 for each scenario,
  sample size and parameter (simulation tables of the paper);
* `tab_sim_S1_S3.tex`: the LaTeX tables.

All 900 replications converged, as did more than 99.9% of the bootstrap
refits. Median running time per replication (minutes, `n` = 100/200/300):
S1 6.4/22.7/34.8; S2 1.7/5.8/13.6; S3 8.3/18.5/26.6.

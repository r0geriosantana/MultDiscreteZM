# Reproducing the paper

Scripts for the applications (Section 5) of *Zero-Modified Multivariate Conway-Maxwell-Poisson Regression
with Baseline-Mean Parameterization* (Santana, Conceição, Diniz, Andrade).
The simulation study (Section 4) is in `inst/simulation` (see the README there).

## Requirements

* `MultDiscreteZM` >= 0.2.1 (exact baseline-mean parameterization)
* `readxl`; `parallel` (base R)
* the data files in `inst/extdata` (see the README there)

## How to run

```r
# copy the scripts to a writable folder
file.copy(system.file("paper", package = "MultDiscreteZM"), "C:/my/folder", recursive = TRUE)
setwd("C:/my/folder/paper")
source("01_alagoas_no_covariates.R")
```

Results (`.rds` and figures) are written to `./results`. Environment variables:
`MZM_B` (bootstrap replicates, default 5000) and `MZM_OUT` (output folder).
For a quick test: `Sys.setenv(MZM_B = 50)`.

| Script | Paper (LaTeX label) | Content |
|---|---|---|
| `00_setup.R` | — | packages, data (Alagoas without Maceió; Minas Gerais joint frequencies, `freq_triv2`) |
| `utils_paper.R` | — | moments, zero probabilities, GDI, regression quantities, parallel bootstrap |
| `01_alagoas_no_covariates.R` | `Ajuste1`, `Descritiva1` | MZMP, MZMNB and MZMCMP fits; sample vs fitted quantities |
| `02_alagoas_bootstrap.R` | `tab6` | bootstrap CIs, MZMCMP and MZMNB |
| `03_alagoas_hdi_regression.R` | `tab:AjusteRegAlagoas` | MZMCMP_mu regression on the HDI, LR test, implied quantities, p_j, bootstrap |
| `04_alagoas_figures.R` | `FitedMarg`, `Vertical` | p_j versus HDI; fitted means with bootstrap bands |
| `05_minas_gerais.R` | `AjusteMG`, `Descritiva2`, `tab:BootMG` | fits, sample vs fitted quantities, bootstrap CIs |
| `06_minas_selection.R` | Section 5.2 | inference conditional on the cohort selection: omega and p estimated from the municipalities with notifications in 2020; bootstrap under the selection design |

## Notes

* The bootstrap runs in parallel (`detectCores() - 1` workers) with
  `clusterSetRNGStream(seed = 2026)`. Results are reproducible for a given
  number of workers; with a different number of cores the percentiles can
  differ in the third decimal place.
* At the geometric boundary (nu -> 0, Minas Gerais) the Hessian of the MZMCMP
  model is singular; Wald standard errors are `NaN` and inference is based on
  the bootstrap.
* `06_minas_selection.R` uses `inst/extdata/coorte_MG_2020.csv`. Municipalities
  without notifications in 2020 (group `A0`) entered the cohort only because of
  their 2023 counts; they contribute only the zero-truncated baseline pmf, so
  omega and p are estimated from the municipalities in group `Ay`.
* Approximate running times (8 cores): bootstrap of the models without
  covariates, a few minutes per model; HDI regression, about 1h15 for 5000
  replicates (version 0.2.0, before the analytic gradient).

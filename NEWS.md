# MultDiscreteZM 0.2.4

* `TrizmpoisI()`, `TrizmnbI()`, `reg3zmpoisI()`, `reg3zmnbI()`, `fitreg3zmpoisI()`,
  `fitreg3zmnbI()`: the count part of the likelihood is computed on the log
  scale. The previous `log(density + 1e-15)` floor biased the fit when an
  observation had a very small density (e.g. the count 21 in the Minas Gerais
  data, where `TrizmpoisI()` returned a non-maximum and an incorrect
  log-likelihood).
* `inst/paper/utils_paper.R`, `inst/simulation/utils_sim.R`: the bootstrap
  generator of the MZMNB model now passes the hurdle probability
  omega = p (1 - P0) to `r3zmnbI()`, whose argument `p` is P(Y != 0).
  Previously the zero-modification parameter was passed, which is wrong for
  p < 1 and fails for p > 1.
* `boot_par()` and `boot_funs_*()`: the MLE, the simulation/refit functions
  and the data they use (`n`, `x`) are evaluated before being sent to the
  parallel workers (previously every replicate failed on the workers).
* New script `inst/paper/06_minas_selection.R` and data file
  `inst/extdata/coorte_MG_2020.csv`: inference for the Minas Gerais
  application conditional on the cohort selection (omega and p estimated
  from the municipalities with notifications in 2020; bootstrap under the
  selection design).

# MultDiscreteZM 0.2.3

* Regression functions (`reg3zmcmpuI()`, `reg3zmcmpI()`, `reg3zmnbI()`,
  `reg3zmpoisI()`): the 95% confidence limits on the original scale are now
  the Wald limits on the optimization scale mapped by the transformation
  (inverse logit for the hurdle probability `w0`/`omega0`, exponential for
  the dispersion parameters), so that they respect the parameter space; the
  interval for `w0` stays in (0, 1). Standard errors on the original scale
  are still given by the delta method.
* `reg3zmpoisI()`: the upper confidence limit was computed with the wrong sign.
* Output tables: column `MLS` renamed to `MLE`; section headers printed with
  `cat()` (no `[1]` prefix).
* The fitting functions no longer print the Q-Q plot of the residuals; the
  plot is still returned as `qq_plot` in the output list.

# MultDiscreteZM 0.2.2

* Simulation study of the paper added in `inst/simulation` (scenarios S1-S3,
  Monte Carlo replications with parametric bootstrap percentile intervals,
  performance measures), with the results reported in the paper in
  `inst/simulation/paper_results`.
* Data of the applications in `inst/extdata`.

# MultDiscreteZM 0.2.1

* Analytic gradient of the MZMCMP_mu regression log-likelihood
  (`src/cmp_grad.h`, `logzmRegIu_vg()`), using
  d log(lambda)/d mu = 1/Var0 and d log(lambda)/d nu = Cov0(Y, log Y!)/Var0.
  `fitreg3zmcmpuI()` and `reg3zmcmpuI()` now pass it to `optim()`, which
  removes the numerical differentiation (about 20 likelihood evaluations per
  gradient). Estimates are unchanged up to the optimizer tolerance.

# MultDiscreteZM 0.2.0

## Exact baseline-mean parameterization (MZMCMP_mu)

* The rate parameter is now obtained exactly as lambda = g^{-1}(mu; nu), the
  solution of E_0(Y) = mu, by Newton's method on log(lambda)
  (dg/dlog(lambda) = Var_0(Y)), with step control and a bracketing fallback.
  The previous asymptotic approximation of Shmueli et al. (2005), which
  required mu > (1 - nu)/(2 nu) and is inaccurate for small lambda, was removed.
  Without covariates the MZMCMP and MZMCMP_mu likelihoods now coincide.
* New exported functions `lambda_u(mu, kappa)` and `cmp_moments(lambda, nu)`.
* Normalizing constant, mean and variance computed on the log scale with
  adaptive truncation (`src/cmp_core.h`); used by the likelihoods, densities
  and random generators.
* Guards in the lambda inversion against invalid or extreme inputs
  (no infinite loops for nu -> 0 or very large mu).

## Other changes

* `fitreg3zmcmpuI()` now returns the attributes `converged` and `optimizer`.


* Hurdle initial values: `y0` is now the indicator I(Y != 0), consistent with
  omega = P(Y != 0) (all MZMCMP fitting functions).
* Initial value of the dispersion parameters set to kappa = 0 (nu = 1).
* `TrizmcmpI()`: estimates on the original scale are computed outside the
  Hessian check (the function failed when the Hessian was singular).

## Paper material

* `inst/paper`: scripts that reproduce the tables and figures of the applications in the paper.
* `inst/extdata`: data files (see the README there).

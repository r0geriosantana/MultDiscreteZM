#' Rapid MLE Fitting for Trivariate MZMCMP Regression (Lambda Parameterization - Simulation)
#'
#' @description A streamlined version of the fitting function designed specifically
#' for Monte Carlo and Bootstrap simulation studies. Bypasses Hessian calculation
#' and console outputs to maximize computational speed.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#' @param x Numeric vector or matrix representing the covariates.
#' @param kappa_init Numeric scalar. The initial value for the dispersion parameters (log scale). Default is -0.5.
#'
#' @return A named numeric vector of the Maximum Likelihood Estimates (MLEs) on the original scale.
#' @author Rogerio A. Santana
#' @export
fitreg3zmcmpI <- function(yx, ys, yz, x, kappa_init = -0.5) {

  # --- 1. INITIALIZATION ---
  pois1 <- stats::glm(yx ~ x, family = stats::poisson())
  pois2 <- stats::glm(ys ~ x, family = stats::poisson())
  pois3 <- stats::glm(yz ~ x, family = stats::poisson())

  xmat <- stats::model.matrix(pois1)
  ka <- ncol(xmat) # Dynamically capture number of regressors

  zero1 <- ifelse(yx == 0, 1, 0)
  zero2 <- ifelse(ys == 0, 1, 0)
  zero3 <- ifelse(yz == 0, 1, 0)
  y0 <- zero1 * zero2 * zero3

  bin0 <- stats::glm.fit(xmat[, 1, drop = FALSE], y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  # Dynamic starting values based on user input
  chute <- c(gama0, beta1, beta2, beta3,
             kappa1 = kappa_init, kappa2 = kappa_init, kappa3 = kappa_init)

  likHood <- function(par) {
    logzmRegTriIcmp(
      par = par,
      xmat = xmat,
      y1 = yx,
      y2 = ys,
      y3 = yz
    )
  }

  # --- 2. OPTIMIZATION WITH FALLBACK ---
  resultado <- try(
    stats::optim(
      par = chute,
      fn = likHood,
      method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12),
      hessian = FALSE # Hessian deactivated for simulation speed
    ),
    silent = TRUE
  )

  if (!inherits(resultado, "try-error") && resultado$convergence == 0) {
    param <- resultado$par
  } else {
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1.0e-7, "maxeval" = 5000)

    resultado_nlopt <- nloptr::nloptr(
      x0 = chute,
      eval_f = likHood,
      lb = rep(-Inf, length(chute)),
      ub = rep(Inf, length(chute)),
      opts = opts
    )
    param <- resultado_nlopt$solution
  }

  # --- 3. PARAMETER EXTRACTION & TRANSFORMATION ---
  g_transform <- function(par_vec) {
    transf_par <- par_vec
    transf_par[1] <- stats::plogis(par_vec[1]) # Transform gamma0 to w0

    # Transform kappas to nu (nu = exp(kappa))
    idx_kappas <- (length(par_vec) - 2):length(par_vec)
    transf_par[idx_kappas] <- exp(par_vec[idx_kappas])

    return(transf_par)
  }

  MLEs <- round(g_transform(param), 4)

  # Build dynamic names to avoid hardcoding
  beta_names <- colnames(xmat)
  if(beta_names[1] == "(Intercept)") beta_names[1] <- "0"
  if(length(beta_names) > 1 && beta_names[2] == "x") beta_names[2] <- "1"

  param_names <- c('w0',
                   paste0('beta1', beta_names),
                   paste0('beta2', beta_names),
                   paste0('beta3', beta_names),
                   'nu1', 'nu2', 'nu3')

  names(MLEs) <- param_names

  return(MLEs)
}

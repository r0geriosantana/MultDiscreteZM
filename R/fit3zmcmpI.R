#' Rapid MLE Fitting for Trivariate MZMCMP without covariates (Simulation Wrapper)
#'
#' @description A streamlined version of the fitting function designed specifically
#' for Monte Carlo and Bootstrap simulation studies. It bypasses Hessian calculation
#' and console outputs to maximize computational speed, returning only the point estimates.
#'
#' @param y1,y2,y3 Numeric vectors representing the count responses.
#' @param zeta_init Numeric scalar. The initial value for the dispersion parameters on the log scale. Default is 0 (implying nu = 1, Poisson start).
#'
#' @return A named numeric vector of the Maximum Likelihood Estimates (MLEs) on the original scale.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#' @export
fit3zmcmpI <- function(y1, y2, y3, zeta_init = 0) {

  # --- 1. INITIALIZATION (Poisson & Binomial Guesses) ---
  pois1 <- stats::glm(y1 ~ 1, family = stats::poisson())
  pois2 <- stats::glm(y2 ~ 1, family = stats::poisson())
  pois3 <- stats::glm(y3 ~ 1, family = stats::poisson())

  xmat <- stats::model.matrix(pois1)
  ka <- ncol(xmat) # Dynamically capture number of regressors

  zero1 <- ifelse(y1 == 0, 1, 0)
  zero2 <- ifelse(y2 == 0, 1, 0)
  zero3 <- ifelse(y3 == 0, 1, 0)
  y0 <- zero1 * zero2 * zero3

  bin0 <- stats::glm.fit(xmat[, 1, drop = FALSE], y0, family = stats::binomial())

  # Extracting and stripping dirty names for a clean numeric vector
  gama0 <- unname(stats::coefficients(bin0))
  beta1 <- unname(stats::coefficients(pois1))
  beta2 <- unname(stats::coefficients(pois2))
  beta3 <- unname(stats::coefficients(pois3))

  # Dynamic starting values based on user input (zeta = log(nu))
  chute <- c(gama0, beta1, beta2, beta3,
             zeta1 = zeta_init, zeta2 = zeta_init, zeta3 = zeta_init)

  # Likelihood wrapper (using the compiled DLL from the package)
  likHood <- function(par) {
    likHoodcmpzmI(par = par, y1 = y1, y2 = y2, y3 = y3)
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
  g_transform <- function(par) {
    gama0 <- par[1]; beta1 <- par[2]; beta2 <- par[3]; beta3 <- par[4]
    zeta1 <- par[5]; zeta2 <- par[6]; zeta3 <- par[7]

    omega0 <- stats::plogis(gama0)
    lambda1 <- exp(beta1); lambda2 <- exp(beta2); lambda3 <- exp(beta3)
    nu1 <- exp(zeta1); nu2 <- exp(zeta2); nu3 <- exp(zeta3)

    dzero <- dcmp_sp(0, lambda1, nu1) * dcmp_sp(0, lambda2, nu2) * dcmp_sp(0, lambda3, nu3)
    denominador_p0 <- 1 - dzero

    p0 <- ifelse(abs(denominador_p0) < 1e-9, NA, omega0 / denominador_p0)
    return(c(p0 = p0, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, nu1 = nu1, nu2 = nu2, nu3 = nu3))
  }

  MLEs <- round(g_transform(param), 4)

  # Build dynamic names
  param_names <- c('p', 'lambda1', 'lambda2', 'lambda3', 'nu1', 'nu2', 'nu3')
  names(MLEs) <- param_names

  return(MLEs)
}

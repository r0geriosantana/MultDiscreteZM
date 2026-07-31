#' Rapid MLE Fitting for Trivariate Poisson Model without covariates (Simulation Wrapper)
#'
#' @description A streamlined version of the fitting function designed specifically
#' for Monte Carlo and Bootstrap simulation studies. It bypasses Hessian calculation
#' and console outputs to maximize computational speed, returning only the point estimates.
#'
#' @param y1,y2,y3 Numeric vectors representing the count responses.
#'
#' @return A named numeric vector of the Maximum Likelihood Estimates (MLEs) on the original scale.
#' @author Rogerio A. Santana
#' @export
fit3zmpoisI <- function(y1, y2, y3) {

  # --- 1. INITIALIZATION (Poisson & Binomial Guesses) ---
  pois1 <- stats::glm(y1 ~ 1, family = stats::poisson())
  pois2 <- stats::glm(y2 ~ 1, family = stats::poisson())
  pois3 <- stats::glm(y3 ~ 1, family = stats::poisson())

  xmat <- stats::model.matrix(pois1)

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

  # Dynamic starting values based on user input
  chute <- c(gama0, beta1, beta2, beta3)

  # --- PRE-CALCULATIONS FOR EXTREME SPEED ---
  y00 <- (y0 == 1) # All zeros
  y0c <- (y0 == 0) # At least one non-zero

  N0 <- sum(y00) # Total number of purely zero occurrences
  Nc <- sum(y0c) # Total number of non-zero occurrences

  # Slicing the count vectors to include only observations > 0
  y1_c <- y1[y0c]
  y2_c <- y2[y0c]
  y3_c <- y3[y0c]

  # --- 2. LIKELIHOOD FUNCTION (Optimized Algebraic Implementation) ---
  likHood <- function(par) {
    omega0 <- stats::plogis(par[1])

    # Lambda parameters
    lambda1 <- exp(par[2])
    lambda2 <- exp(par[3])
    lambda3 <- exp(par[4])

    # Probability (0,0,0) in the Poisson Base
    dz1 <- stats::dpois(0, lambda = lambda1)
    dz2 <- stats::dpois(0, lambda = lambda2)
    dz3 <- stats::dpois(0, lambda = lambda3)
    dzero <- dz1 * dz2 * dz3

    # Logarithmic blocks algebraically simplified
    log_zeros <- N0 * log(1 - omega0 + 1e-15)
    log_omega <- Nc * log(omega0 + 1e-15)
    log_denom <- Nc * log(1 - dzero + 1e-15)

    # dpois(..., log = TRUE) processed internally in C
    log_pois1 <- sum(stats::dpois(y1_c, lambda = lambda1, log = TRUE))
    log_pois2 <- sum(stats::dpois(y2_c, lambda = lambda2, log = TRUE))
    log_pois3 <- sum(stats::dpois(y3_c, lambda = lambda3, log = TRUE))

    logverossimilhanca <- -(log_zeros + log_omega + log_pois1 + log_pois2 + log_pois3 - log_denom)
    return(logverossimilhanca)
  }

  # --- 3. OPTIMIZATION WITH FALLBACK ---
  resultado <- try(
    stats::optim(
      par = chute,
      fn = likHood,
      method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12),
      hessian = FALSE
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

  # --- 4. PARAMETER EXTRACTION & TRANSFORMATION ---
  g_transform <- function(par) {
    omega0 <- stats::plogis(par[1])
    mu1 <- exp(par[2])
    mu2 <- exp(par[3])
    mu3 <- exp(par[4])

    dz1 <- stats::dpois(0, lambda = mu1)
    dz2 <- stats::dpois(0, lambda = mu2)
    dz3 <- stats::dpois(0, lambda = mu3)

    denominador_p0 <- 1 - (dz1 * dz2 * dz3)
    p <- ifelse(abs(denominador_p0) < 1e-9, NA, omega0 / denominador_p0)

    return(c(p = p, mu1 = mu1, mu2 = mu2, mu3 = mu3))
  }

  MLEs <- round(g_transform(param), 4)

  param_names <- c('p', 'mu1', 'mu2', 'mu3')
  names(MLEs) <- param_names

  return(MLEs)
}

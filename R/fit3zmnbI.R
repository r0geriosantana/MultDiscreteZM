#' Rapid MLE Fitting for Trivariate Negative Binomial Model without covariates (Simulation Wrapper)
#'
#' @description A streamlined version of the fitting function designed specifically
#' for Monte Carlo and Bootstrap simulation studies. It bypasses Hessian calculation
#' and console outputs to maximize computational speed, returning only the point estimates.
#'
#' @param y1,y2,y3 Numeric vectors representing the count responses.
#'
#' @return A named numeric vector of the Maximum Likelihood Estimates (MLEs) on the original scale.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#' @export
fit3zmnbI <- function(y1, y2, y3) {

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
  chute <- c(gama0, beta1, beta2, beta3, zeta1 = 0, zeta2 = 0, zeta3 = 0)

  # --- PRÉ-CÁLCULOS PARA VELOCIDADE EXTREMA ---
  # Isolando vetores e contagens fora do loop de otimização
  y00 <- (y0 == 1) # All zeros
  y0c <- (y0 == 0) # At least one non-zero

  N0 <- sum(y00) # Número total de ocorrências puramente zero
  Nc <- sum(y0c) # Número total de ocorrências não-zero

  # Fatiando os vetores de contagem apenas para as observações > 0
  y1_c <- y1[y0c]
  y2_c <- y2[y0c]
  y3_c <- y3[y0c]

  # --- 2. LIKELIHOOD FUNCTION (Optimized Algebraic Implementation) ---
  likHood <- function(par) {
    omega0 <- stats::plogis(par[1])

    # Parâmetros Escalares
    lambda1 <- exp(par[2]); lambda2 <- exp(par[3]); lambda3 <- exp(par[4])
    theta1 <- exp(par[5]); theta2 <- exp(par[6]); theta3 <- exp(par[7])

    # Probabilidade (0,0,0) na Binomial Negativa Base
    dz1 <- stats::dnbinom(0, size = theta1, mu = lambda1)
    dz2 <- stats::dnbinom(0, size = theta2, mu = lambda2)
    dz3 <- stats::dnbinom(0, size = theta3, mu = lambda3)
    dzero <- dz1 * dz2 * dz3

    # Blocos Logarítmicos
    log_zeros <- N0 * log(1 - omega0 + 1e-15)
    log_omega <- Nc * log(omega0 + 1e-15)
    log_denom <- Nc * log(1 - dzero + 1e-15)

    # dnbinom(..., log = TRUE) é processado em C no R, garantindo estabilidade
    log_nb1 <- sum(stats::dnbinom(y1_c, size = theta1, mu = lambda1, log = TRUE))
    log_nb2 <- sum(stats::dnbinom(y2_c, size = theta2, mu = lambda2, log = TRUE))
    log_nb3 <- sum(stats::dnbinom(y3_c, size = theta3, mu = lambda3, log = TRUE))

    logverossimilhanca <- -(log_zeros + log_omega + log_nb1 + log_nb2 + log_nb3 - log_denom)
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
    mu1 <- exp(par[2]); mu2 <- exp(par[3]); mu3 <- exp(par[4])
    theta1 <- exp(par[5]); theta2 <- exp(par[6]); theta3 <- exp(par[7])

    dz1 <- stats::dnbinom(0, size = theta1, mu = mu1)
    dz2 <- stats::dnbinom(0, size = theta2, mu = mu2)
    dz3 <- stats::dnbinom(0, size = theta3, mu = mu3)

    denominador_p0 <- 1 - (dz1 * dz2 * dz3)
    p <- ifelse(abs(denominador_p0) < 1e-9, NA, omega0 / denominador_p0)

    return(c(p = p, mu1 = mu1, mu2 = mu2, mu3 = mu3, theta1 = theta1, theta2 = theta2, theta3 = theta3))
  }

  MLEs <- round(g_transform(param), 4)

  param_names <- c('p', 'mu1', 'mu2', 'mu3', 'theta1', 'theta2', 'theta3')
  names(MLEs) <- param_names

  return(MLEs)
}

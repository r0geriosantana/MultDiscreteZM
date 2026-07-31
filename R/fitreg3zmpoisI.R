#' Extract MLEs for Trivariate ZM-Poisson Regression (Simulation Wrapper)
#'
#' @description Fast point-estimate extractor for ZMP Monte Carlo simulations.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#' @param x Numeric vector or matrix representing the covariates.
#'
#' @return A named numeric vector of Maximum Likelihood Estimates (original scale).
#' @export
fitreg3zmpoisI <- function(yx, ys, yz, x) {

  # --- 1. Initialization ---
  pois1 <- stats::glm(yx ~ x, family = stats::poisson())
  pois2 <- stats::glm(ys ~ x, family = stats::poisson())
  pois3 <- stats::glm(yz ~ x, family = stats::poisson())

  y1 <- pois1$y; y2 <- pois2$y; y3 <- pois3$y
  Xmat <- stats::model.matrix(pois1)
  ka <- ncol(Xmat)

  zero1 <- ifelse(y1 == 0, 1, 0)
  zero2 <- ifelse(y2 == 0, 1, 0)
  zero3 <- ifelse(y3 == 0, 1, 0)
  y0 <- zero1 * zero2 * zero3

  X0 <- Xmat[, 1, drop = FALSE]
  bin0 <- stats::glm.fit(X0, y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  chute <- c(gama0, beta1, beta2, beta3)

  # --- 2. Likelihood Function ---
  likHood <- function(par) {
    gama_0 <- par[1]
    beta_1 <- par[2:(1 + ka)]
    beta_2 <- par[(2 + ka):(1 + 2 * ka)]
    beta_3 <- par[(2 + 2 * ka):(1 + 3 * ka)]

    lambda1 <- as.numeric(exp(Xmat %*% beta_1))
    lambda2 <- as.numeric(exp(Xmat %*% beta_2))
    lambda3 <- as.numeric(exp(Xmat %*% beta_3))

    omega0 <- as.numeric(exp(X0 %*% gama_0) / (1 + exp(X0 %*% gama_0)))

    y00 <- (y0 == 1)
    y0c <- (y0 == 0)

    dzero <- stats::dpois(0, lambda = lambda1) * stats::dpois(0, lambda = lambda2) * stats::dpois(0, lambda = lambda3)

    vetor_dzeros <- (1 - omega0)[y00]

    cmp_contagem <- (stats::dpois(y1, lambda = lambda1) * stats::dpois(y2, lambda = lambda2) * stats::dpois(y3, lambda = lambda3))[y0c]

    logvetor_zeros <- sum(log(vetor_dzeros + 1e-15))
    logcmp_contagem <- sum(log(cmp_contagem + 1e-15))
    logomega_zero <- sum(log(omega0[y0c] + 1e-15))
    logdenom <- sum(log(1 - dzero[y0c] + 1e-15))

    return(-(logvetor_zeros + logomega_zero + logcmp_contagem - logdenom))
  }

  # --- 3. High-Speed Optimization (No Hessian) ---
  resultado <- try(
    stats::optim(
      par = chute, fn = likHood, method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12), hessian = FALSE
    ), silent = TRUE
  )

  if (!inherits(resultado, "try-error") && resultado$convergence == 0) {
    param <- resultado$par
  } else {
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1e-7, "maxeval" = 5000)
    resultado_nlopt <- nloptr::nloptr(
      x0 = chute, eval_f = likHood,
      lb = rep(-Inf, length(chute)), ub = rep(Inf, length(chute)), opts = opts
    )
    param <- resultado_nlopt$solution
  }

  # --- 4. Transform to Original Scale & Format Output ---
  w0 <- as.numeric(stats::plogis(param[1]))
  betas <- param[2:length(param)]

  MLE <- c(w0 = w0, betas)

  # Dynamic naming for betas
  beta_names <- colnames(Xmat)
  if(beta_names[1] == "(Intercept)") beta_names[1] <- "0"
  if(length(beta_names) > 1 && beta_names[2] == "x") beta_names[2] <- "1"

  b_names <- c(paste0("beta1", beta_names), paste0("beta2", beta_names), paste0("beta3", beta_names))
  names(MLE) <- c("w0", b_names)

  return(MLE)
}

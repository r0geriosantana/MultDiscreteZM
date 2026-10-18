#' Maximum Likelihood Estimation for Trivariate ZMP Regression (Type I)
#'
#' @description Fits the trivariate Zero-Modified Poisson regression
#' model incorporating covariates for the mean components, and generates
#' randomized quantile residuals for diagnostics.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#' @param x Numeric vector or matrix representing the covariates.
#'
#' @return A list containing parameter estimates, standard errors, model information, and residuals.
#' @export
reg3zmpoisI <- function(yx, ys, yz, x) {

  # --- 1. INITIALIZATION ---
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

  X0 <- Xmat[, 1, drop = FALSE] # Intercept only for the zero-modification
  bin0 <- stats::glm.fit(X0, y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  chute <- c(gama0, beta1, beta2, beta3)

  # --- 2. LIKELIHOOD FUNCTION ---
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

    zero_1 <- stats::dpois(0, lambda = lambda1)
    zero_2 <- stats::dpois(0, lambda = lambda2)
    zero_3 <- stats::dpois(0, lambda = lambda3)

    dzero <- zero_1 * zero_2 * zero_3
    cmpzero <- (1 - omega0)

    vetor_dzeros <- cmpzero[y00]

    log_parte_naozero <- stats::dpois(y1, lambda = lambda1, log = TRUE) + stats::dpois(y2, lambda = lambda2, log = TRUE) + stats::dpois(y3, lambda = lambda3, log = TRUE)

    log_cmp_contagem <- log_parte_naozero[y0c]

    logvetor_zeros <- sum(log(vetor_dzeros + 1e-15))
    logcmp_contagem <- sum(log_cmp_contagem)   # log densities: no floor
    logomega_zero <- sum(log(omega0[y0c] + 1e-15))
    logdenom <- sum(log(1 - dzero[y0c] + 1e-15))

    return(-(logvetor_zeros + logomega_zero + logcmp_contagem - logdenom))
  }

  # --- 3. OPTIMIZATION ---
  resultado <- try(
    stats::optim(
      par = chute, fn = likHood, method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12), hessian = TRUE
    ), silent = TRUE
  )

  opt_method <- "optim"
  if (!inherits(resultado, "try-error") && resultado$convergence == 0) {
    cat("Optimization successful using 'optim'.\n")
    parametros <- resultado$par
    log_vero <- -resultado$value
    hessiana <- resultado$hessian
    convergencia <- resultado$convergence
  } else {
    opt_method <- "nloptr"
    cat("'optim' failed or did not converge. Falling back to 'nloptr'...\n")
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1e-7, "maxeval" = 5000)
    resultado_nlopt <- nloptr::nloptr(
      x0 = chute, eval_f = likHood,
      lb = rep(-Inf, length(chute)), ub = rep(Inf, length(chute)), opts = opts
    )
    parametros <- resultado_nlopt$solution
    log_vero <- -resultado_nlopt$objective
    convergencia <- resultado_nlopt$status
    cat("Calculating Hessian at the nloptr solution...\n")
    hessiana <- numDeriv::hessian(func = likHood, x = parametros)
  }

  # --- 4. STANDARD ERRORS ---
  erro_pad <- tryCatch(
    { sqrt(diag(solve(hessiana))) },
    error = function(e) {
      cat("\n[Warning] Singular Hessian detected. Using Moore-Penrose Pseudoinverse...\n")
      tryCatch(
        { sqrt(diag(MASS::ginv(hessiana))) },
        error = function(e2) { rep(NA, length(parametros)) }
      )
    }
  )

  beta_names <- colnames(Xmat)
  if(beta_names[1] == "(Intercept)") beta_names[1] <- "Intercept"
  Parameters_names <- c('gamma0', paste0('beta1_', beta_names), paste0('beta2_', beta_names), paste0('beta3_', beta_names) )

  Z <- parametros / erro_pad
  p_valor <- 2 * stats::pnorm(abs(Z), lower.tail = FALSE)

  tabela_resumo <- data.frame(Parameters = Parameters_names, MLE = round(parametros, 4), SE = round(erro_pad, 4), Z_Value = round(Z, 4), P_Value = round(p_valor, 4))

  posto <- length(parametros)
  tabela_info <- data.frame(AIC = round(2 * (-log_vero) + 2 * posto, 3), BIC = round(2 * (-log_vero) + posto * log(length(y1)), 3), Log_likelihood = round(log_vero, 3), N_Obs = length(y1))

  cat("--- Estimated Parameters Summary (Optimization Scale) ---\n")
  print(tabela_resumo, row.names = FALSE)

  # --- 5. DELTA METHOD ---
  tabela_ic <- NULL
  if (!any(is.na(erro_pad))) {
    V_gamma_beta <- solve(hessiana)
    g_transform <- function(par) {
      transf_par <- par
      transf_par[1] <- stats::plogis(par[1]) # omega0 transform
      # NOTE: No 'thetas/zetas' to transform in Poisson model
      return(transf_par)
    }
    estimativas_originais <- g_transform(parametros)
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = parametros)
    V_p_mu <- G_jacobiana %*% V_gamma_beta %*% t(G_jacobiana)
    erros_padrao_orig <- sqrt(diag(V_p_mu))

    z_critico <- stats::qnorm(0.975)
    # Wald limits on the optimization scale, mapped by the (componentwise,
    # increasing) transformation: the interval for w0 stays in (0, 1)
    ic_inferior <- g_transform(parametros - z_critico * erro_pad)
    ic_superior <- g_transform(parametros + z_critico * erro_pad)

    Parameters_orig <- Parameters_names
    Parameters_orig[1] <- "omega0"

    tabela_ic <- data.frame(Parameters = Parameters_orig, MLE = round(estimativas_originais, 4), SE = round(erros_padrao_orig, 4), LI_95 = round(ic_inferior, 4), LS_95 = round(ic_superior, 4))
    cat("\n--- Original Scale (SE: delta method; 95% CI: transformed Wald limits) ---\n")
    print(tabela_ic, row.names = FALSE)
  }

  cat("\n--- Model Information ---\n")
  print(tabela_info, row.names = FALSE)

  # --- 6. RANDOMIZED QUANTILE RESIDUALS ---
  omega0_est <- stats::plogis(parametros[1])
  beta1_est <- parametros[2:(1 + ka)]
  beta2_est <- parametros[(2 + ka):(1 + 2 * ka)]
  beta3_est <- parametros[(2 + 2 * ka):(1 + 3 * ka)]

  mu1_vec <- as.numeric(exp(Xmat %*% beta1_est))
  mu2_vec <- as.numeric(exp(Xmat %*% beta2_est))
  mu3_vec <- as.numeric(exp(Xmat %*% beta3_est))

  # Calculate observation-specific Hurdle probability (p_mod)
  d0_1 <- stats::dpois(0, lambda = mu1_vec)
  d0_2 <- stats::dpois(0, lambda = mu2_vec)
  d0_3 <- stats::dpois(0, lambda = mu3_vec)

  denominador_p0 <- 1 - (d0_1 * d0_2 * d0_3)
  p0_vec <- ifelse(abs(denominador_p0) < 1e-9, NA, omega0_est / denominador_p0)

  # Vectorized CDF and PMF
  Fy <- p3zmpoisI(yx - 1, ys - 1, yz - 1, mu1_vec, mu2_vec, mu3_vec, p0_vec)
  fy <- d3zmpoisI(yx, ys, yz, mu1_vec, mu2_vec, mu3_vec, p0_vec)

  residuo_uniforme <- Fy + fy * stats::runif(length(yx))
  residuo_uniforme <- pmax(1e-10, pmin(1 - 1e-10, residuo_uniforme))
  residuos_normais <- stats::qnorm(residuo_uniforme)

  plot_residual <- qq_plot_zm(residuos_normais)

  return(invisible(list(
    estimates_optim = tabela_resumo,
    estimates_original = tabela_ic,
    model_info = tabela_info,
    hessian = hessiana,
    residuals = residuos_normais,
    qq_plot = plot_residual,
    convergence = convergencia,
    optimizer = opt_method
  )))
}

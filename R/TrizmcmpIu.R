#' Maximum Likelihood Estimation for the Trivariate MZMCMPu (Type I)
#'
#' @description Fits the trivariate Zero-Modified COM-Poisson distribution (mean-parameterized)
#' without covariates. Automatically handles optimization fallbacks and diagnostic plots.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#'
#' @return A list containing parameter estimates, standard errors, model information, and the Q-Q plot.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#' @export
TrizmcmpIu <- function(yx, ys, yz) {

  # --- 1. INITIALIZATION (Poisson Guesses) ---
  pois1 <- stats::glm(yx ~ 1, family = stats::poisson())
  pois2 <- stats::glm(ys ~ 1, family = stats::poisson())
  pois3 <- stats::glm(yz ~ 1, family = stats::poisson())

  # Intersections of zeros
  zero1 <- ifelse(pois1$y == 0, 1, 0)
  zero2 <- ifelse(pois2$y == 0, 1, 0)
  zero3 <- ifelse(pois3$y == 0, 1, 0)

  y0 <- 1 - zero1 * zero2 * zero3   # indicator I(Y != 0): omega = P(Y != 0)

  bin0 <- stats::glm.fit(stats::model.matrix(pois1), y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  chute <- c(gama0, beta1, beta2, beta3, kappa1 = 0, kappa2 = 0, kappa3 = 0)

  # Wrapper for the C++ Likelihood
  likHood <- function(par) {
    likHoodcmpuzmI(par = par, y1 = yx, y2 = ys, y3 = yz)
  }

  opt_method <- "optim"

  # --- 2. OPTIMIZATION WITH FALLBACK ---
  resultado <- try(
    stats::optim(
      par = chute,
      fn = likHood,
      method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12),
      hessian = TRUE
    ),
    silent = TRUE
  )

  if (!inherits(resultado, "try-error") && resultado$convergence == 0) {
    cat("Optimization successful using 'optim'.\n")
    parametros <- resultado$par
    log_vero <- -resultado$value
    hessiana <- resultado$hessian
    convergencia <- resultado$convergence

  } else {
    opt_method <- "nloptr"
    cat("'optim' failed or did not converge. Falling back to 'nloptr'...\n")

    opts <- list(
      "algorithm" = "NLOPT_LN_NELDERMEAD",
      "xtol_rel" = 1.0e-8,
      "maxeval" = 15000  # Giving much more breathing room to complex data.
    )

    resultado_nlopt <- nloptr::nloptr(
      x0 = chute, eval_f = likHood,
      lb = rep(-Inf, length(chute)), ub = rep(Inf, length(chute)),
      opts = opts
    )

    parametros <- resultado_nlopt$solution
    log_vero <- -resultado_nlopt$objective
    convergencia <- resultado_nlopt$status

    cat("Calculating the Hessian in the nloptr solution...\n")
    # We reduced the step size to avoid accidentally hitting the C++ barrier
    hessiana <- numDeriv::hessian(func = likHood, x = parametros, method.args = list(eps = 1e-6))
  }

  # --- 3. CONVERGENCE CHECK ---
  converged <- FALSE
  if (opt_method == "optim" && convergencia == 0) {
    converged <- TRUE
  } else if (opt_method == "nloptr" && convergencia >= 1 && convergencia <= 4) {
    converged <- TRUE
  }

  cat("\n------------------------ Convergence Check ----------------------\n")
  cat("Optimizer:", opt_method, "| Status Code:", convergencia, "\n")
  if (!converged) {
    cat("ATTENTION: Optimization may not have converged correctly!\n")
  } else {
    cat("Optimization converged successfully!\n")
  }
  cat("-----------------------------------------------------------------\n\n")

  # --- 4. STANDARD ERRORS & WALD TEST ---
  erro_pad <- tryCatch(
    sqrt(diag(solve(hessiana))),
    error = function(e) {
      warning("Could not calculate standard errors: ", e$message)
      rep(NA, length(parametros))
    }
  )

  Z <- parametros / erro_pad
  p_valor <- 2 * stats::pnorm(abs(Z), lower.tail = FALSE)

  Parameters <- c('gamma', 'beta1', 'beta2', 'beta3', 'kappa1', 'kappa2', 'kappa3')
  tabela_resumo <- data.frame(
    Parameters = Parameters,
    MLE = round(parametros, 4),
    SE = round(erro_pad, 4),
    Z_Value = round(Z, 4),
    P_Value = round(p_valor, 4)
  )

  posto <- length(parametros)
  tabela_info <- data.frame(
    AIC = round(2 * (-log_vero) + 2 * posto, 4),
    BIC = round(2 * (-log_vero) + posto * log(length(yx)), 4),
    Log_likelihood = round(log_vero, 4),
    N_Obs = length(yx)
  )

  cat("--- Estimated Parameters Summary ---\n")
  print(tabela_resumo, row.names = FALSE)
  cat("\n--- Model Information ---\n")
  print(tabela_info, row.names = FALSE)

  # --- 5. DELTA METHOD (ORIGINAL SCALE) ---
  g_transform <- function(par) {
    gama0 <- par[1]; beta1 <- par[2]; beta2 <- par[3]; beta3 <- par[4]
    kappa1 <- par[5]; kappa2 <- par[6]; kappa3 <- par[7]

    omega0 <- stats::plogis(gama0)
    mu1 <- exp(beta1); mu2 <- exp(beta2); mu3 <- exp(beta3)

    dz1 <- dcmp_u(0, mu1, kappa1)
    dz2 <- dcmp_u(0, mu2, kappa2)
    dz3 <- dcmp_u(0, mu3, kappa3)
    dzero <- dz1 * dz2 * dz3

    denominador_p0 <- 1 - dzero
    p0 <- ifelse(abs(denominador_p0) < 1e-9, NA, omega0 / denominador_p0)

    return(c(p0 = p0, mu1 = mu1, mu2 = mu2, mu3 = mu3, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3))
  }

  # The extraction of estimates MUST be left out of the if() statement to avoid breaking the function.
  estimativas_originais <- g_transform(parametros)

  tabela_ic <- NULL
  if (!any(is.na(erro_pad))) {
    V_gamma_beta <- solve(hessiana)
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = parametros)
    V_p_mu <- G_jacobiana %*% V_gamma_beta %*% t(G_jacobiana)
    erros_padrao_orig <- sqrt(diag(V_p_mu))

    z_critico <- stats::qnorm(0.975)
    ic_inferior <- estimativas_originais - z_critico * erros_padrao_orig
    ic_superior <- estimativas_originais + z_critico * erros_padrao_orig

    Parameters_orig <- c('p', 'mu1', 'mu2', 'mu3', 'kappa1', 'kappa2', 'kappa3')
    tabela_ic <- data.frame(
      Parameters = Parameters_orig,
      MLE = round(estimativas_originais, 4),
      SE = round(erros_padrao_orig, 4),
      LI_95 = round(ic_inferior, 4),
      LS_95 = round(ic_superior, 4)
    )
    cat("\n--- Confidence Intervals (95%) via Delta Method ---\n")
    print(tabela_ic, row.names = FALSE)
  } else {
    cat("\n--- Confidence Intervals could not be calculated (Singular Hessian) ---\n")
  }

  # --- 6. THEORETICAL MOMENTS ---
  p0 <- as.numeric(estimativas_originais[1])
  mu1 <- as.numeric(estimativas_originais[2]); mu2 <- as.numeric(estimativas_originais[3]); mu3 <- as.numeric(estimativas_originais[4])
  kappa1 <- as.numeric(estimativas_originais[5]); kappa2 <- as.numeric(estimativas_originais[6]); kappa3 <- as.numeric(estimativas_originais[7])

  # Translation to original scale for internal helper functions
  nu1 <- exp(kappa1); nu2 <- exp(kappa2); nu3 <- exp(kappa3)
  lambda1 <- lambda_u(mu1, kappa1)   # exact: E0(Y1) = mu1
  lambda2 <- lambda_u(mu2, kappa2)
  lambda3 <- lambda_u(mu3, kappa3)

  E1 <- p0 * EXcmp(lambda1, nu1); E2 <- p0 * EXcmp(lambda2, nu2); E3 <- p0 * EXcmp(lambda3, nu3)

  var_func <- function(p, lambda, nu) { p * (dz2cmp(lambda, nu) + EXcmp(lambda, nu)) - (p * EXcmp(lambda, nu))^2 }
  var1 <- var_func(p0, lambda1, nu1); var2 <- var_func(p0, lambda2, nu2); var3 <- var_func(p0, lambda3, nu3)

  cor1x2 <- corzmI(p0, lambda1, lambda2, nu1, nu2)
  cor1x3 <- corzmI(p0, lambda1, lambda3, nu1, nu3)
  cor2x3 <- corzmI(p0, lambda2, lambda3, nu2, nu3)

  ID <- GDIzmcmpI(p0, lambda1, lambda2, lambda3, nu1, nu2, nu3)

  cat("\n------------ Expected Values ------------\n")
  cat("E(Y1):", round(E1, 4), "| E(Y2):", round(E2, 4), "| E(Y3):", round(E3, 4), "\n")
  cat("-------------- Variances ----------------\n")
  cat("Var(Y1):", round(var1, 4), "| Var(Y2):", round(var2, 4), "| Var(Y3):", round(var3, 4), "\n")
  cat("------------- Correlations --------------\n")
  cat("Cor(Y1, Y2):", round(cor1x2, 4), "| Cor(Y1, Y3):", round(cor1x3, 4), "| Cor(Y2, Y3):", round(cor2x3, 4), "\n")
  cat("------- Generalized Dispersion Index ----\n")
  cat("GDI:", round(ID, 4), "\n")

  # --- 7. RESIDUALS & QQ PLOT ---
  n <- length(yx)
  residuos_normais <- numeric(n)

  for (i in 1:n) {
    # Utilizing the optimized CDF and PMF functions directly
    Fy <- p3zmcmpuI(yx[i] - 1, ys[i] - 1, yz[i] - 1, mu1, mu2, mu3, kappa1, kappa2, kappa3, p0)
    fy <- d3zmcmpuI(yx[i], ys[i], yz[i], mu1, mu2, mu3, kappa1, kappa2, kappa3, p0)

    residuo_uniforme <- Fy + fy * stats::runif(1)
    residuo_uniforme <- pmax(1e-10, pmin(1 - 1e-10, residuo_uniforme))
    residuos_normais[i] <- stats::qnorm(residuo_uniforme)
  }

  plot_residual <- qq_plot_zm(residuos_normais)

  # --- 8. RETURN STRUCTURED OBJECT ---
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

#' Maximum Likelihood Estimation for the Trivariate MZMCMP (Type I)
#'
#' @description Fits the trivariate Zero-Modified COM-Poisson distribution without covariates
#' using Maximum Likelihood. Automatically handles optimization fallbacks and diagnostic plots.
#'
#' Maximum Likelihood Estimation for the Trivariate MZMCMP (Type I)
#'
#' @description Fits the trivariate Zero-Modified COM-Poisson distribution without covariates.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#'
#' @return A list containing parameter estimates, standard errors, model information, and the Q-Q plot.
#' @author Rogerio A. Santana
#' @export
TrizmcmpI <- function(yx, ys, yz) {

  # --- 1. INITIALIZATION (Poisson Guesses) ---
  pois1 <- stats::glm(yx ~ 1, family = stats::poisson())
  pois2 <- stats::glm(ys ~ 1, family = stats::poisson())
  pois3 <- stats::glm(yz ~ 1, family = stats::poisson())

  # Intersections of zeros
  zero1 <- ifelse(pois1$y == 0, 1, 0)
  zero2 <- ifelse(pois2$y == 0, 1, 0)
  zero3 <- ifelse(pois3$y == 0, 1, 0)

  y0 <- zero1 * zero2 * zero3

  bin0 <- stats::glm.fit(stats::model.matrix(pois1), y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  chute <- c(gama0, beta1, beta2, beta3, nu1 = 1, nu2 = 1, nu3 = 1)

  # Wrapper for the C++ Likelihood
  likHood <- function(par) {
    likHoodcmpzmI(par = par, y1 = yx, y2 = ys, y3 = yz)
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

    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1.0e-7, "maxeval" = 5000)

    resultado_nlopt <- nloptr::nloptr(
      x0 = chute, eval_f = likHood,
      lb = rep(-Inf, length(chute)), ub = rep(Inf, length(chute)),
      opts = opts
    )

    parametros <- resultado_nlopt$solution
    log_vero <- -resultado_nlopt$objective
    convergencia <- resultado_nlopt$status

    cat("Calculating Hessian at the nloptr solution...\n")
    hessiana <- numDeriv::hessian(func = likHood, x = parametros)
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
    MLS = round(parametros, 4),
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

  print("--- Estimated Parameters Summary ---")
  print(tabela_resumo, row.names = FALSE)
  print("\n--- Model Information ---")
  print(tabela_info, row.names = FALSE)

  # --- 5. DELTA METHOD (ORIGINAL SCALE) ---
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

  tabela_ic <- NULL
  if (!any(is.na(erro_pad))) {
    V_gamma_beta <- solve(hessiana)
    estimativas_originais <- g_transform(parametros)
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = parametros)
    V_p_lambda <- G_jacobiana %*% V_gamma_beta %*% t(G_jacobiana)
    erros_padrao_orig <- sqrt(diag(V_p_lambda))

    z_critico <- stats::qnorm(0.975)
    ic_inferior <- estimativas_originais - z_critico * erros_padrao_orig
    ic_superior <- estimativas_originais + z_critico * erros_padrao_orig

    Parameters_orig <- c('p', 'lambda1', 'lambda2', 'lambda3', 'nu1', 'nu2', 'nu3')
    tabela_ic <- data.frame(
      Parameters = Parameters_orig,
      MLS = round(estimativas_originais, 4),
      SE = round(erros_padrao_orig, 4),
      LI_95 = round(ic_inferior, 4),
      LS_95 = round(ic_superior, 4)
    )
    cat("\n--- Confidence Intervals (95%) via Delta Method ---\n")
    print(tabela_ic, row.names = FALSE)
  }

  # --- 6. THEORETICAL MOMENTS ---
  p0 <- as.numeric(estimativas_originais[1])
  lambda1 <- as.numeric(estimativas_originais[2]); lambda2 <- as.numeric(estimativas_originais[3]); lambda3 <- as.numeric(estimativas_originais[4])
  nu1 <- as.numeric(estimativas_originais[5]); nu2 <- as.numeric(estimativas_originais[6]); nu3 <- as.numeric(estimativas_originais[7])

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
    Fy <- p3zmcmpI(yx[i] - 1, ys[i] - 1, yz[i] - 1, lambda1, lambda2, lambda3, nu1, nu2, nu3, p0)
    fy <- d3zmcmpI(yx[i], ys[i], yz[i], lambda1, lambda2, lambda3, nu1, nu2, nu3, p0)

    residuo_uniforme <- Fy + fy * stats::runif(1)
    residuo_uniforme <- pmax(1e-10, pmin(1 - 1e-10, residuo_uniforme))
    residuos_normais[i] <- stats::qnorm(residuo_uniforme)
  }

  plot_residual <- qq_plot_zm(residuos_normais)
  print(plot_residual)

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

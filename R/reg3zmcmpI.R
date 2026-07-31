#' Maximum Likelihood Estimation for the Trivariate MZMCMP Regression (Lambda Parameterization)
#'
#' @description Fits the trivariate Zero-Modified COM-Poisson regression model
#' (lambda-parameterized) incorporating covariates for the intensity components.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#' @param x Numeric vector or matrix representing the covariates.
#'
#' @return A list containing parameter estimates, standard errors, and model information.
#' @author Rogerio A. Santana
#' @export
reg3zmcmpI <- function(yx, ys, yz, x) {

  # --- 1. INITIALIZATION (Poisson & Binomial Guesses) ---
  pois1 <- stats::glm(yx ~ x, family = stats::poisson())
  pois2 <- stats::glm(ys ~ x, family = stats::poisson())
  pois3 <- stats::glm(yz ~ x, family = stats::poisson())

  y1 <- pois1$y
  y2 <- pois2$y
  y3 <- pois3$y
  xmat <- stats::model.matrix(pois1)
  ka <- ncol(xmat) # Dynamically capture number of regressors

  zero1 <- ifelse(y1 == 0, 1, 0)
  zero2 <- ifelse(y2 == 0, 1, 0)
  zero3 <- ifelse(y3 == 0, 1, 0)
  y0 <- zero1 * zero2 * zero3

  bin0 <- stats::glm.fit(xmat[, 1, drop = FALSE], y0, family = stats::binomial())

  gama0 <- stats::coefficients(bin0)
  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)

  # Initialize kappas slightly negative to account for common overdispersion
  chute <- c(gama0, beta1, beta2, beta3, kappa1 = -0.5, kappa2 = -0.5, kappa3 = -0.5)

  likHood <- function(par) {
    logzmRegTriIcmp(par = par, xmat = xmat, y1 = y1, y2 = y2, y3 = y3)
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
    param <- resultado$par
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

    param <- resultado_nlopt$solution
    log_vero <- -resultado_nlopt$objective
    convergencia <- resultado_nlopt$status

    cat("Calculating Hessian at the nloptr solution...\n")
    hessiana <- numDeriv::hessian(func = likHood, x = param)
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

  # --- 4. STANDARD ERRORS & WALD TEST (Optimization Scale) ---
  erro_pad <- tryCatch(
    sqrt(diag(solve(hessiana))),
    error = function(e) {
      warning("Could not calculate standard errors: ", e$message)
      rep(NA, length(param))
    }
  )

  beta_names <- colnames(xmat)
  if(beta_names[1] == "(Intercept)") beta_names[1] <- "Intercept"

  Parameters <- c('gamma0',
                  paste0('beta1_', beta_names),
                  paste0('beta2_', beta_names),
                  paste0('beta3_', beta_names),
                  'nu1', 'nu2', 'nu3')

  Z <- param / erro_pad
  p_valor <- 2 * stats::pnorm(abs(Z), lower.tail = FALSE)

  tabela_resumo <- data.frame(
    Parameters = Parameters,
    MLS = round(param, 4),
    SE = round(erro_pad, 4),
    Z_Value = round(Z, 4),
    P_Value = round(p_valor, 4)
  )

  posto <- length(param)
  tabela_info <- data.frame(
    AIC = round(2 * (-log_vero) + 2 * posto, 4),
    BIC = round(2 * (-log_vero) + posto * log(length(y1)), 4),
    Log_likelihood = round(log_vero, 4),
    N_Obs = length(y1)
  )

  print("--- Estimated Parameters Summary (Optimization Scale) ---")
  print(tabela_resumo, row.names = FALSE)
  print("\n--- Model Information ---")
  print(tabela_info, row.names = FALSE)

  # --- 5. DELTA METHOD (Original Scale) ---
  g_transform <- function(par_vec) {
    transf_par <- par_vec
    transf_par[1] <- stats::plogis(par_vec[1]) # gamma0 -> w0
    # kappas to nu (nu = exp(kappa))
    idx_kappas <- (length(par_vec) - 2):length(par_vec)
    transf_par[idx_kappas] <- exp(par_vec[idx_kappas])
    return(transf_par)
  }

  estimativas_originais <- g_transform(param)
  tabela_ic <- NULL

  if (!any(is.na(erro_pad))) {
    V_gamma_beta <- solve(hessiana)
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = param)
    V_p_mu <- G_jacobiana %*% V_gamma_beta %*% t(G_jacobiana)
    erros_padrao_orig <- sqrt(diag(V_p_mu))

    z_critico <- stats::qnorm(0.975)
    ic_inferior <- estimativas_originais - z_critico * erros_padrao_orig
    ic_superior <- estimativas_originais + z_critico * erros_padrao_orig

    Parameters_orig <- Parameters
    Parameters_orig[1] <- "w0"
    idx_kappas <- (length(Parameters_orig) - 2):length(Parameters_orig)
    Parameters_orig[idx_kappas] <- c("nu1", "nu2", "nu3")

    tabela_ic <- data.frame(
      Parameters = Parameters_orig,
      MLS = round(estimativas_originais, 4),
      SE = round(erros_padrao_orig, 4),
      LI_95 = round(ic_inferior, 4),
      LS_95 = round(ic_superior, 4)
    )
    cat("\n--- Confidence Intervals (95%) via Delta Method (Original Scale) ---\n")
    print(tabela_ic, row.names = FALSE)
  }
  # --- 6. RESIDUALS & QQ PLOT ---
  gama0 <- param[1]
  beta1 <- param[2:(1+ka)]
  beta2 <- param[(2+ka):(1+2*ka)]
  beta3 <- param[(2+2*ka):(1+3*ka)]

  # CORREÇÃO: O 'param' armazena os kappas. Precisamos aplicar o exp() para obter o nu.
  kappa1 <- param[length(param) - 2]
  kappa2 <- param[length(param) - 1]
  kappa3 <- param[length(param)]

  nu1 <- exp(kappa1)
  nu2 <- exp(kappa2)
  nu3 <- exp(kappa3)

  lambda1 <- as.numeric(exp(xmat %*% beta1))
  lambda2 <- as.numeric(exp(xmat %*% beta2))
  lambda3 <- as.numeric(exp(xmat %*% beta3))
  w0 <- stats::plogis(gama0)

  n <- length(y1)
  residuos_normais <- numeric(n)

  for (i in 1:n) {
    # Calcula a probabilidade de zero-modificação específica da observação (p_mod)
    d0_1 <- dcmp_sp(0, lambda1[i], nu1)
    d0_2 <- dcmp_sp(0, lambda2[i], nu2)
    d0_3 <- dcmp_sp(0, lambda3[i], nu3)

    denominador <- 1 - (d0_1 * d0_2 * d0_3)
    p_mod_i <- ifelse(abs(denominador) < 1e-9, NA, w0 / denominador)

    # Calcula a FDA e FMP para os resíduos quantílicos randomizados
    Fy <- p3zmcmpI(yx[i] - 1, ys[i] - 1, yz[i] - 1, lambda1[i], lambda2[i], lambda3[i], nu1, nu2, nu3, p_mod_i)
    fy <- d3zmcmpI(yx[i], ys[i], yz[i], lambda1[i], lambda2[i], lambda3[i], nu1, nu2, nu3, p_mod_i)

    residuo_uniforme <- Fy + fy * stats::runif(1)
    residuo_uniforme <- pmax(1e-10, pmin(1 - 1e-10, residuo_uniforme))
    residuos_normais[i] <- stats::qnorm(residuo_uniforme)
  }

  plot_residual <- qq_plot_zm(residuos_normais)
  print(plot_residual)


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

#' Fit Trivariate Zero-Modified Poisson Model (Intercept Only)
#'
#' @description Fits a trivariate ZMP model without covariates, estimating
#' theoretical moments, generalized dispersion index, and quantile residuals.
#'
#' @param yx,ys,yz Numeric vectors representing the count responses.
#'
#' @return A list containing parameter estimates, standard errors, model metrics, and theoretical moments.
#' @export
TrizmpoisI <- function(yx, ys, yz) {

  # --- 1. INITIALIZATION (Poisson & Binomial Guesses) ---
  pois1 <- stats::glm(yx ~ 1, family = stats::poisson())
  pois2 <- stats::glm(ys ~ 1, family = stats::poisson())
  pois3 <- stats::glm(yz ~ 1, family = stats::poisson())

  y1 <- pois1$y; Xmat1 <- stats::model.matrix(pois1)
  y2 <- pois2$y; Xmat2 <- stats::model.matrix(pois2)
  y3 <- pois3$y; Xmat3 <- stats::model.matrix(pois3)

  # Structural zeros intersection
  zero1 <- ifelse(y1 == 0, 1, 0)
  zero2 <- ifelse(y2 == 0, 1, 0)
  zero3 <- ifelse(y3 == 0, 1, 0)
  y0 <- zero1 * zero2 * zero3

  bin0 <- stats::glm.fit(Xmat1, y0, family = stats::binomial())

  beta1 <- stats::coefficients(pois1)
  beta2 <- stats::coefficients(pois2)
  beta3 <- stats::coefficients(pois3)
  gama0 <- stats::coefficients(bin0)

  chute <- c(gama0, beta1, beta2, beta3)

  # --- 2. LIKELIHOOD FUNCTION (R-based implementation) ---
  likHood <- function(par) {
    gama_0 <- par[1]
    beta_1 <- par[2]
    beta_2 <- par[3]
    beta_3 <- par[4]

    # Lambda parameters
    lambda1 <- exp(Xmat1 %*% beta_1)
    lambda2 <- exp(Xmat2 %*% beta_2)
    lambda3 <- exp(Xmat3 %*% beta_3)

    # Zero-modification parameter (Probability of Y > 0)
    omega0 <- exp(Xmat1 %*% gama_0) / (1 + exp(Xmat1 %*% gama_0))

    y00 <- (y0 == 1) # All zeros
    y0c <- (y0 == 0) # At least one non-zero

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

    logverossimilhanca <- -(logvetor_zeros + logomega_zero + logcmp_contagem - logdenom)
    return(logverossimilhanca)
  }

  # --- 3. OPTIMIZATION WITH FALLBACK ---
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
  } else {
    cat("'optim' failed or did not converge. Falling back to 'nloptr'...\n")
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1.0e-7, "maxeval" = 5000)
    resultado_nlopt <- nloptr::nloptr(
      x0 = chute, eval_f = likHood,
      lb = rep(-Inf, length(chute)), ub = rep(Inf, length(chute)),
      opts = opts
    )
    parametros <- resultado_nlopt$solution
    log_vero <- -resultado_nlopt$objective
    cat("Calculating Hessian at the nloptr solution...\n")
    hessiana <- numDeriv::hessian(func = likHood, x = parametros)
  }

  # --- 4. STANDARD ERRORS & WALD TEST ---
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

  parameters_names <- c('gamma0', 'beta1', 'beta2', 'beta3')
  Z <- parametros / erro_pad
  p_valor <- 2 * stats::pnorm(abs(Z), lower.tail = FALSE)

  tabela_resumo <- data.frame(
    Parameters = parameters_names,
    MLE = round(parametros, 4),
    SE = round(erro_pad, 4),
    Z_Value = round(Z, 4),
    P_Value = round(p_valor, 4)
  )

  posto <- length(parametros)
  tabela_info <- data.frame(
    AIC = round(2 * (-log_vero) + 2 * posto, 4),
    BIC = round(2 * (-log_vero) + posto * log(length(ys)), 4),
    Log_likelihood = round(log_vero, 4),
    N_Obs = length(ys)
  )

  cat("--- Estimated Parameters Summary ---\n")
  print(tabela_resumo, row.names = FALSE)
  cat("\n--- Model Information ---\n")
  print(tabela_info, row.names = FALSE)

  # --- 5. DELTA METHOD (Original Scale) ---
  tabela_ic <- NULL
  if (!any(is.na(erro_pad))) {
    V_gamma_beta <- solve(hessiana)

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

    estimativas_originais <- g_transform(parametros)
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = parametros)
    V_p_mu <- G_jacobiana %*% V_gamma_beta %*% t(G_jacobiana)
    erros_padrao_orig <- sqrt(diag(V_p_mu))

    z_critico <- stats::qnorm(0.975)
    ic_inferior <- estimativas_originais - z_critico * erros_padrao_orig
    ic_superior <- estimativas_originais + z_critico * erros_padrao_orig

    tabela_ic <- data.frame(
      Parameters = c('p', 'mu1', 'mu2', 'mu3'),
      MLE = round(estimativas_originais, 4),
      SE = round(erros_padrao_orig, 4),
      LI_95 = round(ic_inferior, 4),
      LS_95 = round(ic_superior, 4)
    )

    cat("\n--- Confidence Intervals (95%) via Delta Method ---\n")
    print(tabela_ic, row.names = FALSE)
  }

  # --- 6. THEORETICAL MOMENTS ---
  p <- as.numeric(estimativas_originais[1])
  mu1 <- as.numeric(estimativas_originais[2])
  mu2 <- as.numeric(estimativas_originais[3])
  mu3 <- as.numeric(estimativas_originais[4])

  E1 <- p * mu1; E2 <- p * mu2; E3 <- p * mu3

  vari <- function(p_val, mu_val) {
    E_mean <- p_val * mu_val
    M2 <- mu_val + (mu_val^2)
    return((p_val * M2) - (E_mean)^2)
  }

  var1 <- vari(p, mu1)
  var2 <- vari(p, mu2)
  var3 <- vari(p, mu3)

  cor1x2 <- corzmpoisI(p, mu1, mu2)
  cor1x3 <- corzmpoisI(p, mu1, mu3)
  cor2x3 <- corzmpoisI(p, mu2, mu3)
  ID <- GDIzmpoisI(p, mu1, mu2, mu3)

  cat("\n------------ Theoretical Moments -------------------\n")
  cat(sprintf("E(Y1): %.4f | Var(Y1): %.4f\n", E1, var1))
  cat(sprintf("E(Y2): %.4f | Var(Y2): %.4f\n", E2, var2))
  cat(sprintf("E(Y3): %.4f | Var(Y3): %.4f\n", E3, var3))
  cat("\n------------ Theoretical Correlations --------------\n")
  cat(sprintf("Cor(Y1, Y2): %.4f\n", cor1x2))
  cat(sprintf("Cor(Y1, Y3): %.4f\n", cor1x3))
  cat(sprintf("Cor(Y2, Y3): %.4f\n", cor2x3))
  cat("----------------------------------------------------\n")
  cat(sprintf("Generalized Dispersion Index (GDI): %.4f\n", ID))
  cat("----------------------------------------------------\n")

  # --- 7. RESIDUALS & QQ PLOT (Vectorized) ---
  n <- length(yx)

  Fy <- p3zmpoisI(yx - 1, ys - 1, yz - 1, mu1, mu2, mu3, p)
  fy <- d3zmpoisI(yx, ys, yz, mu1, mu2, mu3, p)

  residuo_uniforme <- Fy + fy * stats::runif(n)
  residuo_uniforme <- pmax(1e-10, pmin(1 - 1e-10, residuo_uniforme))
  residuos_normais <- stats::qnorm(residuo_uniforme)

  plot_residual <- qq_plot_zm(residuos_normais)

  return(invisible(list(
    estimates_optim = tabela_resumo,
    estimates_original = tabela_ic,
    model_info = tabela_info,
    residuals = residuos_normais,
    qq_plot = plot_residual
  )))
}

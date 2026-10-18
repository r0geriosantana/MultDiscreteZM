## =============================================================================
## Helper functions for the paper
## "Zero-Modified Multivariate Conway-Maxwell-Poisson Regression with
##  Baseline-Mean Parameterization" (Santana, Conceicao, Diniz, Andrade)
##
## Sourced by sim_00_setup.R and by the parallel workers.
## =============================================================================

## ---- Baseline zero probabilities and variances -----------------------------
## COM-Poisson with exact-mean parameterization (mu, kappa), nu = exp(kappa)
cmp_f0_var <- function(mu, kappa) {
  mom <- cmp_moments(lambda_u(mu, kappa), exp(kappa))
  list(f0 = exp(-mom[, "logS"]), var = mom[, "var"])
}
## Negative binomial with mean mu and size theta: Var = mu + mu^2 / theta
nb_f0_var <- function(mu, theta) {
  list(f0 = (theta / (theta + mu))^theta, var = mu + mu^2 / theta)
}
## Poisson
pois_f0_var <- function(mu) list(f0 = exp(-mu), var = mu)

## ---- Moments and zero probabilities of the Type I model ---------------------
## p: zero-modification parameter; mu, var0, f0: baseline means, variances
## and zero probabilities of the three independent margins.
typeI_summary <- function(p, mu, var0, f0) {
  EY <- p * mu
  VY <- p * var0 + p * (1 - p) * mu^2
  S  <- p * (1 - p) * outer(mu, mu); diag(S) <- VY
  R  <- stats::cov2cor(S)
  pr <- rbind(c(1, 2), c(1, 3), c(2, 3))
  s  <- sqrt(EY)
  list(mean    = EY,
       var     = VY,
       cor     = setNames(R[pr], c("12", "13", "23")),
       P0_marg = 1 - p + p * f0,
       P0_pair = setNames(1 - p + p * f0[pr[, 1]] * f0[pr[, 2]], c("12", "13", "23")),
       P0_all  = 1 - p + p * prod(f0),
       GDI     = drop(t(s) %*% S %*% s) / sum(EY^2))
}

## Sample counterparts
sample_summary <- function(Y) {
  Y  <- as.matrix(Y); Z <- Y == 0
  pr <- rbind(c(1, 2), c(1, 3), c(2, 3))
  R  <- stats::cor(Y)
  m  <- colMeans(Y); S <- stats::cov(Y); s <- sqrt(m)
  list(mean    = m,
       var     = diag(S),
       cor     = setNames(R[pr], c("12", "13", "23")),
       P0_marg = colMeans(Z),
       P0_pair = setNames(c(mean(Z[, 1] & Z[, 2]), mean(Z[, 1] & Z[, 3]),
                            mean(Z[, 2] & Z[, 3])), c("12", "13", "23")),
       P0_all  = mean(rowSums(Z) == 3),
       GDI     = drop(t(s) %*% S %*% s) / sum(m^2))
}

## Model summaries from the fitted parameter vectors
summary_cmp  <- function(est) {
  mu <- est[c("mu1", "mu2", "mu3")]; ka <- est[c("kappa1", "kappa2", "kappa3")]
  b  <- lapply(1:3, function(k) cmp_f0_var(mu[k], ka[k]))
  typeI_summary(est["p"], mu, sapply(b, `[[`, "var"), sapply(b, `[[`, "f0"))
}
summary_nb   <- function(est) {
  mu <- est[c("mu1", "mu2", "mu3")]; th <- est[c("theta1", "theta2", "theta3")]
  b  <- nb_f0_var(mu, th)
  typeI_summary(est["p"], mu, b$var, b$f0)
}
summary_pois <- function(est) {
  mu <- est[c("mu1", "mu2", "mu3")]; b <- pois_f0_var(mu)
  typeI_summary(est["p"], mu, b$var, b$f0)
}

## ---- Regression with one covariate x in the means (MZMCMP_mu) ---------------
## est: named vector (w0, beta10, beta11, beta20, beta21, beta30, beta31,
##      kappa1, kappa2, kappa3) as returned by fitreg3zmcmpuI()
reg_margins <- function(est, x) {
  X  <- cbind(1, x)
  mu <- sapply(list(c("beta10", "beta11"), c("beta20", "beta21"), c("beta30", "beta31")),
               function(b) as.numeric(X %*% est[b]))
  mu <- exp(matrix(mu, ncol = 3))
  ka <- est[c("kappa1", "kappa2", "kappa3")]
  f0 <- v0 <- mu
  for (k in 1:3) { b <- cmp_f0_var(mu[, k], ka[k]); f0[, k] <- b$f0; v0[, k] <- b$var }
  P0 <- f0[, 1] * f0[, 2] * f0[, 3]
  list(mu = mu, f0 = f0, v0 = v0, P0 = P0, p = as.numeric(est["w0"] / (1 - P0)))
}
## Municipality-specific zero-modification parameter p_j = omega / (1 - P0j)
pj_reg <- function(est, x) reg_margins(est, x)$p
## Fitted marginal means E(Y_ij | x_j) = p_j mu_ij
EY_reg <- function(est, x) { m <- reg_margins(est, x); m$p * m$mu }

## Unconditional quantities implied by the regression model, averaging over x
implied_reg <- function(est, x) {
  m  <- reg_margins(est, x); w <- as.numeric(est["w0"])
  EY <- m$p * m$mu
  VY <- m$p * m$v0 + m$p * (1 - m$p) * m$mu^2
  Vt <- sapply(1:3, function(k) mean(VY[, k]) + stats::var(EY[, k]) * (nrow(EY) - 1) / nrow(EY))
  pr <- rbind(c(1, 2), c(1, 3), c(2, 3)); cr <- zp <- numeric(3)
  for (r in 1:3) {
    i <- pr[r, 1]; j <- pr[r, 2]
    cv <- mean(m$p * (1 - m$p) * m$mu[, i] * m$mu[, j]) +
          mean(EY[, i] * EY[, j]) - mean(EY[, i]) * mean(EY[, j])
    cr[r] <- cv / sqrt(Vt[i] * Vt[j])
    zp[r] <- mean(1 - w + w * (m$f0[, i] * m$f0[, j] - m$P0) / (1 - m$P0))
  }
  list(mean = colMeans(EY), var = Vt, cor = setNames(cr, c("12", "13", "23")),
       P0_marg = colMeans(1 - w + w * (m$f0 - m$P0) / (1 - m$P0)),
       P0_pair = setNames(zp, c("12", "13", "23")))
}

## Fast refit for the bootstrap: starts at the original MLE (warm start);
## falls back to the package function if BFGS does not converge.
fit_reg_warm <- function(Y, start, x) {
  xmat <- cbind(1, x)
  par0 <- c(stats::qlogis(start["w0"]), start[-1])
  fg <- MultDiscreteZM:::zmreg_fg(xmat, Y[, 1], Y[, 2], Y[, 3])   # analytic gradient
  o <- try(stats::optim(par0, fg$fn, fg$gr, method = "BFGS",
                        control = list(maxit = 1000, reltol = 1e-10)), silent = TRUE)
  if (inherits(o, "try-error") || o$convergence != 0)
    return(fitreg3zmcmpuI(Y[, 1], Y[, 2], Y[, 3], x))
  out <- o$par; out[1] <- stats::plogis(out[1]); names(out) <- names(start)
  attr(out, "converged") <- TRUE; attr(out, "optimizer") <- "optim (warm start)"
  out
}

## Hurdle probability omega = p (1 - P0) implied by a fitted MZMNB vector
omega_nb <- function(est) {
  P0 <- prod(stats::dnbinom(0, size = est[c("theta1", "theta2", "theta3")],
                            mu = est[c("mu1", "mu2", "mu3")]))
  as.numeric(est["p"] * (1 - P0))
}

## ---- Simulation / refit functions used by the bootstrap ---------------------
## Each returns list(sim = function(est), fit = function(Y, est)).
## Data needed by the workers (n, x) travel inside the closures.
boot_funs_cmp <- function(n) list(n = force(n),   # forced: the closures go to the workers
  sim = function(est) as.matrix(r3zmcmpIu(n, est["mu1"], est["mu2"], est["mu3"],
                                          est["kappa1"], est["kappa2"], est["kappa3"], est["p"])),
  fit = function(Y, est) fit3zmcmpIu(Y[, 1], Y[, 2], Y[, 3]))

## Note: the 'p' argument of r3zmnbI() is the hurdle probability
## omega = P(Y != 0) = p (1 - P0), not the zero-modification parameter p.
boot_funs_nb <- function(n) list(n = force(n),
  sim = function(est) as.matrix(r3zmnbI(n, est["mu1"], est["mu2"], est["mu3"],
                                        est["theta1"], est["theta2"], est["theta3"],
                                        omega_nb(est))),
  fit = function(Y, est) fit3zmnbI(Y[, 1], Y[, 2], Y[, 3]))

boot_funs_reg <- function(x) list(x = force(x),
  sim = function(est) as.matrix(r3zmcmpuI_reg0(x,
          beta1 = est[c("beta10", "beta11")], beta2 = est[c("beta20", "beta21")],
          beta3 = est[c("beta30", "beta31")],
          kappa1 = est["kappa1"], kappa2 = est["kappa2"], kappa3 = est["kappa3"],
          w0 = est["w0"])),
  fit = function(Y, est) fit_reg_warm(Y, est, x))

## ---- Parallel parametric bootstrap -----------------------------------------
## est: MLE (named vector); funs: list(sim, fit) from boot_funs_*().
## Progress is saved to 'arquivo' after each block of replicates.
boot_par <- function(est, funs, B = 5000, seed = 2026, arquivo,
                     ncores = max(1, parallel::detectCores() - 1), bloco = 100,
                     utils_file = UTILS_FILE) {
  force(est); force(funs)   # evaluate before the closures are sent to the workers
  cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, suppressMessages(library(MultDiscreteZM)))
  parallel::clusterCall(cl, function(f) { source(f); NULL }, utils_file)
  parallel::clusterSetRNGStream(cl, seed)
  one <- function(b) {
    fb <- try(funs$fit(funs$sim(est), est), silent = TRUE)
    if (inherits(fb, "try-error") || length(fb) != length(est)) rep(NA_real_, length(est))
    else as.numeric(fb)
  }
  res <- NULL; t0 <- Sys.time()
  for (ini in seq(1, B, by = bloco)) {
    idx <- ini:min(B, ini + bloco - 1)
    res <- rbind(res, do.call(rbind, parallel::parLapply(cl, idx, one)))
    colnames(res) <- names(est)
    saveRDS(list(est = est, boot = res, seed = seed), arquivo)
    s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("%d/%d replicates | %.2f s/replicate | ~%.0f min left\n",
                nrow(res), B, s / nrow(res), s / nrow(res) * (B - nrow(res)) / 60))
  }
  invisible(res)
}

## Percentile intervals; kappa columns are also reported as nu = exp(kappa)
boot_ci <- function(est, res, extra = NULL) {
  tudo <- res; estim <- est
  kc <- grep("^kappa", colnames(res), value = TRUE)
  if (length(kc)) {
    nb <- exp(res[, kc, drop = FALSE]); colnames(nb) <- sub("kappa", "nu", kc)
    tudo <- cbind(tudo, nb); estim <- c(estim, setNames(exp(est[kc]), colnames(nb)))
  }
  if (!is.null(extra)) {               # derived quantities: list(name = function(theta))
    for (nm in names(extra)) {
      tudo <- cbind(tudo, apply(res, 1, extra[[nm]])); colnames(tudo)[ncol(tudo)] <- nm
      estim <- c(estim, setNames(extra[[nm]](est), nm))
    }
  }
  ic <- t(apply(tudo, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  data.frame(MLE = round(estim, 3), CI_low = round(ic[, 1], 3), CI_up = round(ic[, 2], 3),
             failures = sum(!stats::complete.cases(res)))
}

## Small formatter for side-by-side summaries
tab_summary <- function(...) {
  L <- list(...)
  do.call(cbind, lapply(L, function(s) round(unlist(s[c("mean", "var", "cor",
                                                         "P0_marg", "P0_pair", "P0_all", "GDI")]), 3)))
}

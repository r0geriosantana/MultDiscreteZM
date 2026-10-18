## =============================================================================
## sim_01_cenarios.R - Simulation study (Section 4), step 1
## Model (Eq. sim_model), scenarios and the behaviour of p_j, GDI and rho_ij
## as functions of the covariate (scenario figures).
##
## Model: Y_j ~ MZMCMP_mu(p_j, mu_1j, mu_2j, mu_3j, kappa_1, kappa_2, kappa_3)
##   logit(omega_0) = gamma_0                       (constant hurdle P(Y != 0))
##   p_j            = omega_0 / (1 - prod_i P0(mu_ij, kappa_i))
##   log(mu_ij)     = beta_i0 + beta_i1 x_j,  x_j ~ U(0, 1)
##   kappa_ij       = kappa_i
## With the exact mean parameterization there is no restriction on mu or
## kappa, and any omega_0 in (0, 1) defines a valid distribution.
## =============================================================================
source("sim_00_setup.R")

## ---- Scenarios (true values) ------------------------------------------------
## est vector in the format used by fitreg3zmcmpuI():
## (w0, beta10, beta11, beta20, beta21, beta30, beta31, kappa1, kappa2, kappa3)
mk <- function(w0, b, k) setNames(c(w0, t(b), k),
  c("w0", "beta10", "beta11", "beta20", "beta21", "beta30", "beta31",
    "kappa1", "kappa2", "kappa3"))

scenarios <- list(
  ## S1: zero inflation, overdispersed margins (p ~ 0.70)
  S1 = list(true = mk(0.70, rbind(c(1.0, -0.5), c(0.9, -0.5), c(1.1, -0.6)),
                      c(-0.7, -0.5, -0.9)), xlim = c(0, 1)),
  ## S2: zero deflation, overdispersed margins (p from ~1.02 to ~1.32)
  S2 = list(true = mk(0.85, rbind(c(-0.4, -0.6), c(-0.5, -0.5), c(-0.3, -0.7)),
                      c(-0.7, -0.5, -0.9)), xlim = c(0, 1)),
  ## S3: zero inflation, underdispersed margins (p ~ 0.80, nu ~ 2)
  S3 = list(true = mk(0.80, rbind(c(0.9, -0.5), c(0.8, -0.4), c(1.0, -0.6)),
                      c(0.7, 0.6, 0.8)), xlim = c(0, 1))
)

## ---- Behaviour along the covariate ----------------------------------------
behaviour <- function(true, x) {
  m <- reg_margins(true, x)
  t(sapply(seq_along(x), function(j) {
    s <- typeI_summary(m$p[j], m$mu[j, ], m$v0[j, ], m$f0[j, ])
    c(x = x[j], p = m$p[j], P0_base = m$P0[j], GDI = s$GDI, s$cor,
      mu1 = m$mu[j, 1], mu2 = m$mu[j, 2], mu3 = m$mu[j, 3])
  }))
}

res <- lapply(scenarios, function(s)
  behaviour(s$true, seq(s$xlim[1], s$xlim[2], length.out = 101)))

## Checks: lambda finite (valid configuration) and range of p_j, GDI, rho
for (nm in names(res)) {
  r <- res[[nm]]
  cat(sprintf("%s: p in [%.3f, %.3f] | GDI in [%.3f, %.3f] | rho12 in [%.3f, %.3f] | finite: %s\n",
              nm, min(r[, "p"]), max(r[, "p"]), min(r[, "GDI"]), max(r[, "GDI"]),
              min(r[, "12"]), max(r[, "12"]), all(is.finite(r))))
}

## ---- Figures (one per scenario): (A) p and GDI, (B) correlations ------------
for (nm in names(res)) {
  r <- res[[nm]]
  pdf(file.path(OUT_DIR, sprintf("Fig_cenario_%s.pdf", nm)), width = 10, height = 4)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 4.5))
  plot(r[, "x"], r[, "p"], type = "l", lwd = 2, xlab = expression(x[1]),
       ylab = expression(p[j]), main = "(A)")
  abline(h = 1, lty = 3)
  par(new = TRUE)
  plot(r[, "x"], r[, "GDI"], type = "l", lwd = 2, lty = 2, col = "steelblue4",
       axes = FALSE, xlab = "", ylab = "")
  axis(4, col.axis = "steelblue4"); mtext("GDI", side = 4, line = 2.5, col = "steelblue4")
  legend("topright", c(expression(p[j]), "GDI"), lty = 1:2, lwd = 2,
         col = c("black", "steelblue4"), bty = "n")
  matplot(r[, "x"], r[, c("12", "13", "23")], type = "l", lty = 1:3, lwd = 2,
          col = 1, xlab = expression(x[1]), ylab = expression(rho[ij]), main = "(B)")
  abline(h = 0, lty = 3)
  legend("topright", c(expression(rho[12]), expression(rho[13]), expression(rho[23])),
         lty = 1:3, lwd = 2, bty = "n")
  dev.off()
}
saveRDS(list(scenarios = scenarios, behaviour = res), file.path(OUT_DIR, "sim_scenarios.rds"))

## =============================================================================
## 05 - Minas Gerais without covariates: Tables 8, 9, 10 and 11 (Section 5.2)
## =============================================================================
source("00_setup.R")
y1 <- Y_MG[, 1]; y2 <- Y_MG[, 2]; y3 <- Y_MG[, 3]; n <- nrow(Y_MG)

## Fits (Table 9). At nu -> 0 the Hessian of the MZMCMP model is singular:
## Wald standard errors are NaN and inference relies on the bootstrap.
fit_P   <- TrizmpoisI(y1, y2, y3)
fit_NB  <- TrizmnbI(y1, y2, y3)
fit_CMP <- TrizmcmpIu(y1, y2, y3)
est_P   <- fit3zmpoisI(y1, y2, y3)
est_NB  <- fit3zmnbI(y1, y2, y3)
est_CMP <- fit3zmcmpIu(y1, y2, y3)
ka  <- est_CMP[c("kappa1", "kappa2", "kappa3")]
tab9 <- list(
  MZMP   = est_P, MZMNB = est_NB,
  MZMCMP = c(est_CMP[c("p", "mu1", "mu2", "mu3")],
             setNames(exp(ka), c("nu1", "nu2", "nu3")),
             setNames(lambda_u(est_CMP[c("mu1", "mu2", "mu3")], ka),
                      c("lambda1", "lambda2", "lambda3"))),
  info   = rbind(MZMP = fit_P$model_info, MZMNB = fit_NB$model_info,
                 MZMCMP = fit_CMP$model_info))
print(lapply(tab9, function(z) if (is.numeric(z)) round(z, 3) else z))

## Table 10: sample versus fitted quantities
tab10 <- tab_summary(Sample = sample_summary(Y_MG),
                     MZMNB  = summary_nb(est_NB),
                     MZMCMP = summary_cmp(est_CMP))
print(tab10)

## Table 11: parametric bootstrap for both models
boot_par(est_NB, boot_funs_nb(n), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_minas_nb.rds"))
boot_par(est_CMP, boot_funs_cmp(n), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_minas_cmp.rds"))
b_nb  <- readRDS(file.path(OUT_DIR, "boot_minas_nb.rds"))
b_cmp <- readRDS(file.path(OUT_DIR, "boot_minas_cmp.rds"))
tab11 <- list(MZMNB  = boot_ci(b_nb$est,  b_nb$boot),
              MZMCMP = boot_ci(b_cmp$est, b_cmp$boot),
              GDI    = c(MZMNB = summary_nb(est_NB)$GDI, MZMCMP = summary_cmp(est_CMP)$GDI))
print(tab11)

saveRDS(list(tab9 = tab9, tab10 = tab10, tab11 = tab11),
        file.path(OUT_DIR, "minas_gerais.rds"))

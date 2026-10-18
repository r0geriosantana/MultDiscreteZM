## =============================================================================
## 01 - Alagoas without covariates: Tables 3, 4 and 5 (Section 5.1.1)
## =============================================================================
source("00_setup.R")
y1 <- Y_AL[, 1]; y2 <- Y_AL[, 2]; y3 <- Y_AL[, 3]

## Table 3: observed joint frequencies
print(table(apply(Y_AL, 1, paste, collapse = ",")))

## Full fits (standard errors, AIC, BIC) and point estimates
fit_P   <- TrizmpoisI(y1, y2, y3)
fit_NB  <- TrizmnbI(y1, y2, y3)
fit_CMP <- TrizmcmpIu(y1, y2, y3)            # exact-mean parameterization
est_P   <- fit3zmpoisI(y1, y2, y3)
est_NB  <- fit3zmnbI(y1, y2, y3)
est_CMP <- fit3zmcmpIu(y1, y2, y3)

## Table 4: MLEs, nu = exp(kappa), lambda = g^{-1}(mu; nu), AIC and BIC
ka  <- est_CMP[c("kappa1", "kappa2", "kappa3")]
tab4 <- list(
  MZMP   = est_P,
  MZMNB  = est_NB,
  MZMCMP = c(est_CMP[c("p", "mu1", "mu2", "mu3")],
             setNames(exp(ka), c("nu1", "nu2", "nu3")),
             setNames(lambda_u(est_CMP[c("mu1", "mu2", "mu3")], ka),
                      c("lambda1", "lambda2", "lambda3"))),
  info   = rbind(MZMP = fit_P$model_info, MZMNB = fit_NB$model_info,
                 MZMCMP = fit_CMP$model_info))
print(lapply(tab4, function(z) if (is.numeric(z)) round(z, 3) else z))

## Table 5: sample versus fitted moments, correlations and zero probabilities
tab5 <- tab_summary(Sample = sample_summary(Y_AL),
                    MZMNB  = summary_nb(est_NB),
                    MZMCMP = summary_cmp(est_CMP))
print(tab5)
## Rank correlations (robustness check mentioned in the notes of Table 5)
print(round(cor(Y_AL, method = "kendall"), 3))
print(round(cor(Y_AL, method = "spearman"), 3))

saveRDS(list(est_P = est_P, est_NB = est_NB, est_CMP = est_CMP,
             tab4 = tab4, tab5 = tab5,
             logL_CMP = fit_CMP$model_info$Log_likelihood),
        file.path(OUT_DIR, "alagoas_no_covariates.rds"))

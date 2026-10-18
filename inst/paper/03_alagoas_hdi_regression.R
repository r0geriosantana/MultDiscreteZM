## =============================================================================
## 03 - Alagoas: MZMCMP_mu regression on the HDI, Table 7 (Section 5.1.2)
## =============================================================================
source("00_setup.R")
y1 <- Y_AL[, 1]; y2 <- Y_AL[, 2]; y3 <- Y_AL[, 3]; x <- HDI_AL

## Exploratory correlations with the HDI
print(round(sapply(1:3, function(k) cor(Y_AL[, k], x)), 3))                       # Pearson
print(round(sapply(1:3, function(k) cor(Y_AL[, k], x, method = "kendall")), 3))   # Kendall

## Full fit (standard errors, AIC, BIC) and point estimates
fit_reg <- reg3zmcmpuI(y1, y2, y3, x)
est_reg <- fitreg3zmcmpuI(y1, y2, y3, x)

## Likelihood ratio test for H0: beta11 = beta21 = beta31 = 0
logL0 <- TrizmcmpIu(y1, y2, y3)$model_info$Log_likelihood
LR <- 2 * (fit_reg$model_info$Log_likelihood - logL0)
LRtest <- c(LR = LR, df = 3, p_value = pchisq(LR, 3, lower.tail = FALSE))
print(LRtest)

## Effect of a 0.1 increase in the HDI on the baseline means
print(round(exp(0.1 * est_reg[c("beta11", "beta21", "beta31")]), 2))

## Quantities implied by the fitted model (averaged over the observed HDI)
print(round(unlist(implied_reg(est_reg, x)), 3))

## Municipality-specific zero-modification parameter p_j
alagoas$p_hat <- pj_reg(est_reg, x)
print(summary(alagoas$p_hat))
print(alagoas[alagoas$p_hat > 1, c("Municipio", "IDH", "p_hat")])
hdi_p1 <- tryCatch(uniroot(function(h) pj_reg(est_reg, h) - 1, c(0.40, 0.75))$root,
                   error = function(e) NA)
cat("HDI at which p = 1:", round(hdi_p1, 3), "\n")

## Parametric bootstrap with the HDI fixed (Table 7)
boot_par(est_reg, boot_funs_reg(x), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_alagoas_reg.rds"))
b <- readRDS(file.path(OUT_DIR, "boot_alagoas_reg.rds"))
tab7 <- boot_ci(b$est, b$boot,
                extra = list(dif_b21_b11 = function(t) t["beta21"] - t["beta11"]))
print(tab7)

saveRDS(list(est_reg = est_reg, info = fit_reg$model_info, LRtest = LRtest,
             implied = implied_reg(est_reg, x), p_hat = alagoas[, c("Municipio", "IDH", "p_hat")],
             hdi_p1 = hdi_p1, tab7 = tab7),
        file.path(OUT_DIR, "alagoas_hdi_regression.rds"))

## =============================================================================
## 02 - Alagoas without covariates: parametric bootstrap, Table 6
## =============================================================================
source("00_setup.R")
n <- nrow(Y_AL)
est_CMP <- fit3zmcmpIu(Y_AL[, 1], Y_AL[, 2], Y_AL[, 3])
est_NB  <- fit3zmnbI(Y_AL[, 1], Y_AL[, 2], Y_AL[, 3])

boot_par(est_CMP, boot_funs_cmp(n), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_alagoas_cmp.rds"))
boot_par(est_NB, boot_funs_nb(n), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_alagoas_nb.rds"))

b_cmp <- readRDS(file.path(OUT_DIR, "boot_alagoas_cmp.rds"))
b_nb  <- readRDS(file.path(OUT_DIR, "boot_alagoas_nb.rds"))
tab6 <- list(MZMCMP = boot_ci(b_cmp$est, b_cmp$boot),
             MZMNB  = boot_ci(b_nb$est,  b_nb$boot),
             GDI    = c(MZMCMP = summary_cmp(est_CMP)$GDI, MZMNB = summary_nb(est_NB)$GDI))
print(tab6)
saveRDS(tab6, file.path(OUT_DIR, "table6_alagoas.rds"))

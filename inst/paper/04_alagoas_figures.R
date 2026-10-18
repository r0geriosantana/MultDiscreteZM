## =============================================================================
## 04 - Alagoas: Figures 5 (p_j versus HDI) and 6 (fitted marginal means)
## Requires the bootstrap file produced by 03_alagoas_hdi_regression.R
## =============================================================================
source("00_setup.R")
x <- HDI_AL
b <- readRDS(file.path(OUT_DIR, "boot_alagoas_reg.rds"))
est_reg <- b$est
hgrid <- seq(min(x), max(x), length.out = 150)

## Figure 5
p_grid <- pj_reg(est_reg, hgrid)
pdf(file.path(OUT_DIR, "Fig_pmod_HDI.pdf"), width = 6, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
plot(hgrid, p_grid, type = "l", lwd = 2, xlab = "HDI",
     ylab = expression(hat(p)[j]), ylim = c(min(p_grid), 1.02))
points(x, pj_reg(est_reg, x), pch = 16, cex = 0.6)
abline(h = 1, lty = 2)
dev.off()

## Figure 6: fitted means with 95% pointwise bootstrap bands (1000 replicates)
set.seed(1)
idx   <- sample(which(complete.cases(b$boot)), 1000)
bands <- array(NA, c(length(idx), length(hgrid), 3))
for (r in seq_along(idx)) {
  v <- try(EY_reg(b$boot[idx[r], ], hgrid), silent = TRUE)
  if (!inherits(v, "try-error")) bands[r, , ] <- v
}
lo  <- apply(bands, c(2, 3), quantile, 0.025, na.rm = TRUE)
hi  <- apply(bands, c(2, 3), quantile, 0.975, na.rm = TRUE)
fit <- EY_reg(est_reg, hgrid)
titles <- c("(A) HIV-positive pregnant women", "(B) Children exposed to HIV",
            "(C) Congenital syphilis")
pdf(file.path(OUT_DIR, "Vertical_index2015.pdf"), width = 11, height = 3.8)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
for (k in 1:3) {
  plot(x, Y_AL[, k], pch = 16, cex = 0.7, col = "grey40",
       ylim = c(0, max(Y_AL[, k], hi[, k])), xlab = "HDI",
       ylab = "Number of cases", main = titles[k], cex.main = 1)
  polygon(c(hgrid, rev(hgrid)), c(lo[, k], rev(hi[, k])),
          col = adjustcolor("steelblue", 0.25), border = NA)
  lines(hgrid, fit[, k], lwd = 2, col = "steelblue4")
}
dev.off()

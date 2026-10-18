## =============================================================================
## sim_03_medidas.R - Simulation study (Section 4), step 3
## Performance measures and LaTeX tables (one table per scenario).
##   Mean   = average of the MLEs over the N replications
##   MAPE   = 100/N * sum |(theta_hat - theta) / theta|
##   MSE    = 1/N * sum (theta_hat - theta)^2
##   Width  = average width of the 95% bootstrap percentile intervals
##   CP95   = proportion of intervals that contain the true value
## Measures use the replications in which the fit converged; the convergence
## rate is reported in the table caption and a sensitivity table uses all
## replications.
## Requires the files results/sim_<scenario>_n<n>.rds from sim_02_monte_carlo.R.
## =============================================================================
source("sim_00_setup.R")

labels_tex <- c(w0 = "$\\omega_0$",
                beta10 = "$\\beta_{10}$", beta11 = "$\\beta_{11}$",
                beta20 = "$\\beta_{20}$", beta21 = "$\\beta_{21}$",
                beta30 = "$\\beta_{30}$", beta31 = "$\\beta_{31}$",
                kappa1 = "$\\kappa_1$", kappa2 = "$\\kappa_2$", kappa3 = "$\\kappa_3$")

## keep = which replications enter the measures (default: converged ones)
measures <- function(res, keep = res$converged) {
  th  <- res$true
  E   <- res$est[keep, names(th), drop = FALSE]
  L   <- res$low[keep, names(th), drop = FALSE]; U <- res$up[keep, names(th), drop = FALSE]
  err <- sweep(E, 2, th)
  data.frame(n        = res$n,
             par      = names(th),
             Value    = th,
             Mean     = colMeans(E),
             MAPE     = 100 * colMeans(abs(sweep(err, 2, th, "/"))),
             MSE      = colMeans(err^2),
             Width    = colMeans(U - L, na.rm = TRUE),
             CP95     = colMeans(sweep(L, 2, th, "<=") & sweep(U, 2, th, ">="), na.rm = TRUE),
             conv     = sum(res$converged) / res$N,      # convergence rate over all N
             row.names = NULL)
}

files <- list.files(OUT_DIR, pattern = "^sim_S[0-9]+_n[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("No sim_S*_n*.rds files in ", OUT_DIR, ": run sim_02_monte_carlo.R first.")
all_res <- lapply(files, readRDS)
tab <- do.call(rbind, lapply(all_res, function(r) cbind(scenario = r$scenario, measures(r))))
tab <- tab[order(tab$scenario, tab$n, match(tab$par, names(labels_tex))), ]

## Convergence: rate of the main fits and of the bootstrap refits
conv_tab <- do.call(rbind, lapply(all_res, function(r)
  data.frame(scenario = r$scenario, n = r$n,
             conv_rate = sum(r$converged) / r$N, errors = r$n_error,
             boot_conv_mean = round(mean(r$boot_conv), 3),
             minutes = round(r$minutes, 1))))
print(conv_tab)

## Sensitivity: measures using all replications that returned estimates
## (if they differ noticeably from the main table, report it in the text)
tab_all <- do.call(rbind, lapply(all_res, function(r)
  cbind(scenario = r$scenario, measures(r, keep = rep(TRUE, nrow(r$est))))))
print(format(tab, digits = 3), row.names = FALSE)
saveRDS(list(main = tab, all = tab_all, convergence = conv_tab),
        file.path(OUT_DIR, "sim_measures.rds"))

## ---- LaTeX tables (paste into Section 4) -----------------------------------
fmt <- function(v) ifelse(v < 0, sprintf("$%.3f$", v), sprintf("%.3f", v))
for (s in unique(tab$scenario)) {
  ts <- tab[tab$scenario == s, ]
  linhas <- character(0)
  for (nn in unique(ts$n)) {
    tn <- ts[ts$n == nn, ]
    mid <- ceiling(nrow(tn) / 2)
    for (k in seq_len(nrow(tn))) {
      linhas <- c(linhas, sprintf("\t\t%s & %s & %s & %s & %.3f & %.3f & %.3f & %.3f \\\\",
        if (k == mid) nn else "", labels_tex[tn$par[k]], fmt(tn$Value[k]), fmt(tn$Mean[k]),
        tn$MAPE[k], tn$MSE[k], tn$Width[k], tn$CP95[k]))
    }
    if (nn != tail(unique(ts$n), 1)) linhas <- c(linhas, "\t\t\\midrule")
  }
  tex <- c("\\begin{table}[pos=h]", "\t\\centering",
           sprintf("\t\\caption{Simulation results for the $\\mathrm{MZMCMP}_{\\mu}$ regression model, scenario %s. Convergence rates: %s.}",
                   s, paste(sprintf("%.2f ($n=%d$)", conv_tab$conv_rate[conv_tab$scenario == s],
                                    conv_tab$n[conv_tab$scenario == s]), collapse = ", ")),
           sprintf("\t\\label{tab:sim_%s}", s),
           "\t\\begin{tabular}{@{}cccccccc@{}}", "\t\t\\toprule",
           "\t\t$n$ & $\\bm{\\theta}$ & Value & Mean & MAPE & MSE & Width & $\\mathrm{CP}_{95\\%}$ \\\\",
           "\t\t\\midrule", linhas, "\t\t\\bottomrule", "\t\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(OUT_DIR, sprintf("tab_sim_%s.tex", s)))
  cat("\nwritten", sprintf("tab_sim_%s.tex", s), "\n")
}

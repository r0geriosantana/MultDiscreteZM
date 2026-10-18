## =============================================================================
## 06 - Minas Gerais: inference conditional on the cohort selection (Section 5.2)
##
## The cohort is A(y) U B(y): municipalities with notifications in 2020 (A(y))
## or in 2023 (B(y)). Municipalities in A(0) (no notifications in 2020) enter
## the cohort only because they had cases in 2023, so a joint zero in 2023 is
## impossible for them. The likelihood is therefore specified conditionally:
##   * A(y): full Type I likelihood;
##   * A(0): zero-truncated baseline pmf f+(y) only.
## Under the hurdle representation, the baseline parameters (mu, theta / kappa)
## keep their MLE from the full fit; only omega changes, being estimated from
## the municipalities in A(y). The zero-modification parameter p then refers to
## the municipalities with reporting activity in 2020.
##
## Classification of the 76 municipalities: inst/extdata/coorte_MG_2020.csv
## (notifications in 2020 from SINAN/DATASUS, extracted in 2026).
## =============================================================================
source("00_setup.R")

## ---- Data by municipality ---------------------------------------------------
mg_file <- system.file("extdata", "dados_prontosMG.xlsx", package = "MultDiscreteZM")
if (!nzchar(mg_file)) mg_file <- file.path(PAPER_DIR, "..", "extdata", "dados_prontosMG.xlsx")
coorte_file <- system.file("extdata", "coorte_MG_2020.csv", package = "MultDiscreteZM")
if (!nzchar(coorte_file)) coorte_file <- file.path(PAPER_DIR, "..", "extdata", "coorte_MG_2020.csv")

mg     <- as.data.frame(readxl::read_excel(mg_file))
coorte <- utils::read.csv(coorte_file)
mg     <- merge(mg, coorte[, c("cod_ibge", "Grupo")], by = "cod_ibge")
stopifnot(nrow(mg) == 76, all(mg$Grupo %in% c("Ay", "A0")))

## Rows ordered with A(y) first (needed by the bootstrap generator)
mg  <- mg[order(mg$Grupo != "Ay"), ]
Y   <- as.matrix(mg[, c("Esquistossomose", "LeishmanioseVisceral", "Hanseniase")])
colnames(Y) <- c("Y1", "Y2", "Y3"); rownames(Y) <- NULL
inA <- mg$Grupo == "Ay"
nA  <- sum(inA); n0 <- sum(!inA)

## Same data as Table 8 (only the order of the rows differs)
key <- function(M) sort(unname(apply(M, 1, paste, collapse = "-")))
stopifnot(identical(key(Y), key(Y_MG)))
## Selection rule: no joint zero in 2023 among the municipalities in A(0)
stopifnot(all(rowSums(Y[!inA, , drop = FALSE]) > 0))

cat(sprintf("A(y): %d municipalities (%d joint zeros) | A(0): %d municipalities\n",
            nA, sum(rowSums(Y[inA, ]) == 0), n0))

## ---- Fits: full-sample MLE and selection-corrected p ------------------------
y1 <- Y[, 1]; y2 <- Y[, 2]; y3 <- Y[, 3]
est_P   <- fit3zmpoisI(y1, y2, y3)
est_NB  <- fit3zmnbI(y1, y2, y3)
est_CMP <- fit3zmcmpIu(y1, y2, y3)

sel_P   <- sel_correct(est_P,   Y, inA, "pois")
sel_NB  <- sel_correct(est_NB,  Y, inA, "nb")
sel_CMP <- sel_correct(est_CMP, Y, inA, "cmp")

omega_full <- mean(rowSums(Y) > 0)
omega_A    <- mean(rowSums(Y[inA, ]) > 0)
tab_p <- rbind(
  omega = c(full = omega_full, conditional = omega_A),
  p_MZMP   = c(est_P["p"],   sel_P["p"]),
  p_MZMNB  = c(est_NB["p"],  sel_NB["p"]),
  p_MZMCMP = c(est_CMP["p"], sel_CMP["p"]))
print(round(tab_p, 3))

## Log-likelihood of the conditional model: only the Bernoulli part changes,
## so AIC and BIC differ from Table 9 by the same constant for all models.
ll_bern <- function(w, n1, n) n1 * log(w) + (n - n1) * log(1 - w)
delta_ll <- ll_bern(omega_A, sum(rowSums(Y[inA, ]) > 0), nA) -
            ll_bern(omega_full, sum(rowSums(Y) > 0), nrow(Y))
cat(sprintf("Change in the maximized log-likelihood (all models): %.3f\n", delta_ll))

## ---- Parametric bootstrap under the selection design -------------------------
boot_par(sel_NB, boot_funs_sel(nA, n0, "nb"), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_minas_sel_nb.rds"))
boot_par(sel_CMP, boot_funs_sel(nA, n0, "cmp"), B = B_BOOT, seed = SEED,
         arquivo = file.path(OUT_DIR, "boot_minas_sel_cmp.rds"))
b_nb  <- readRDS(file.path(OUT_DIR, "boot_minas_sel_nb.rds"))
b_cmp <- readRDS(file.path(OUT_DIR, "boot_minas_sel_cmp.rds"))

tab_sel <- list(MZMNB  = boot_ci(b_nb$est,  b_nb$boot),
                MZMCMP = boot_ci(b_cmp$est, b_cmp$boot))
print(tab_sel)

## Cohort-level fitted moments (means, variances, correlations, zero
## probabilities, GDI) are unchanged: they are linear in p, and the average of
## the selection-corrected p over A(y) and A(0) equals the full-sample p.

saveRDS(list(tab_p = tab_p, delta_ll = delta_ll, tab_sel = tab_sel,
             nA = nA, n0 = n0),
        file.path(OUT_DIR, "minas_gerais_selection.rds"))

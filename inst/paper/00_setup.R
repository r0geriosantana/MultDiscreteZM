## =============================================================================
## 00_setup.R - packages, paths and data for all paper scripts
## Copy the folder with the scripts to a working directory and run from it:
##   file.copy(system.file("paper", package = "MultDiscreteZM"), "C:/my/folder",
##             recursive = TRUE); setwd("C:/my/folder/paper")
## Results (.rds, .pdf) are written to OUT_DIR (default: ./results).
## =============================================================================
suppressMessages({
  library(MultDiscreteZM)
  library(parallel)
})

PAPER_DIR  <- getwd()
UTILS_FILE <- normalizePath(file.path(PAPER_DIR, "utils_paper.R"))
OUT_DIR    <- Sys.getenv("MZM_OUT", file.path(PAPER_DIR, "results"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(UTILS_FILE)

B_BOOT <- as.integer(Sys.getenv("MZM_B", "5000"))   # bootstrap replicates
SEED   <- 2026

## ---- Alagoas, 2015 (Section 5.1) -------------------------------------------
## Y1 = HIV-positive pregnant women, Y2 = children exposed to HIV,
## Y3 = congenital syphilis; covariate = municipal HDI.
## The capital, Maceio, is excluded (101 of the 102 municipalities).
alagoas_file <- system.file("extdata", "Vertical_2015_Alagoas.xlsx",
                            package = "MultDiscreteZM")
if (!nzchar(alagoas_file))                      # local copy of the repository
  alagoas_file <- file.path(PAPER_DIR, "..", "extdata", "Vertical_2015_Alagoas.xlsx")
alagoas_raw <- as.data.frame(readxl::read_excel(alagoas_file))
alagoas     <- alagoas_raw[!grepl("^Macei", alagoas_raw$Municipio), ]   # exclude Maceio
Y_AL   <- as.matrix(alagoas[, c("Gestantes", "Criancas", "Congenita")])
colnames(Y_AL) <- c("Y1", "Y2", "Y3")
HDI_AL <- alagoas$IDH
stopifnot(nrow(Y_AL) == 101, sum(rowSums(Y_AL) > 0) == 60)

## ---- Minas Gerais, 2023 (Section 5.2) ---------------------------------------
## Y1 = schistosomiasis, Y2 = human visceral leishmaniasis, Y3 = leprosy.
## Joint frequencies of Table 8: each row is (Y1, Y2, Y3, frequency).
tab_mg <- matrix(c(
  0,0,0,6,  0,0,1,9,  0,0,2,7,  0,0,3,7,  0,0,4,4,  0,0,5,1,
  0,0,8,1,  0,0,10,1, 0,0,22,1, 0,1,0,3,  0,1,2,1,  0,2,0,1,
  0,2,1,1,  0,3,2,1,  1,0,0,2,  1,0,1,4,  1,0,2,1,  1,0,3,1,
  1,0,4,2,  1,0,6,1,  1,0,20,1, 1,1,7,1,  2,0,10,1, 2,1,0,1,
  2,1,1,1,  2,2,2,1,  2,0,0,3,  3,0,0,1,  2,0,1,2,  3,0,1,2,
  3,0,2,2,  5,0,0,1,  5,0,2,1,  8,1,1,1,  9,0,1,1,  21,0,1,1),
  ncol = 4, byrow = TRUE)
Y_MG <- tab_mg[rep(seq_len(nrow(tab_mg)), tab_mg[, 4]), 1:3]
colnames(Y_MG) <- c("Y1", "Y2", "Y3")
stopifnot(nrow(Y_MG) == 76)

cat("Setup OK | results in:", OUT_DIR, "| B =", B_BOOT,
    "| cores =", max(1, detectCores() - 1), "\n")

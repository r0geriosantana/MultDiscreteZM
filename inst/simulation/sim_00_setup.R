## =============================================================================
## sim_00_setup.R - setup for the simulation study (Section 4)
## Independent of the application data. Run the scripts from this folder:
##   file.copy(system.file("simulation", package = "MultDiscreteZM"), "C:/my/folder",
##             recursive = TRUE); setwd("C:/my/folder/simulation")
## Requires MultDiscreteZM >= 0.2.1 (exact mean parameterization and the
## 'converged' attribute of fitreg3zmcmpuI()).
## =============================================================================
suppressMessages({
  library(MultDiscreteZM)
  library(parallel)
})
stopifnot(packageVersion("MultDiscreteZM") >= "0.2.1")   # analytic gradient

SIM_DIR    <- getwd()
UTILS_FILE <- normalizePath(file.path(SIM_DIR, "utils_sim.R"))
OUT_DIR    <- Sys.getenv("MZM_SIM_OUT", file.path(SIM_DIR, "results"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(UTILS_FILE)
SEED <- 2026
cat("Simulation setup OK | results in:", OUT_DIR,
    "| cores =", max(1, detectCores() - 1), "\n")

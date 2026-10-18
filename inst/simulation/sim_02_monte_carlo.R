## =============================================================================
## sim_02_monte_carlo.R - Simulation study (Section 4), step 2
## Monte Carlo replications: generate, fit and bootstrap.
##
## Design: scenarios S1-S3 (sim_01_cenarios.R), n in {100, 200, 300},
##         N = 100 replications, B = 500 bootstrap samples per replication,
##         95% parametric bootstrap percentile intervals.
##
## How it runs (version 2):
## * every replication is an independent task sent to a free core
##   (load balancing, no waiting for slow replications);
## * each replication is saved as soon as it ends, in results/parts/
##   (progress = number of files in that folder; a stopped run resumes
##   from the replications not yet saved);
## * time limits: main fit MAIN_LIMIT seconds, each bootstrap refit
##   BOOT_LIMIT seconds; a fit that exceeds the limit or does not converge
##   is recorded as non-converged (it is NOT replaced by a slower method);
## * at the end, the replications of each (scenario, n) are combined into
##   results/sim_<scenario>_n<n>.rds, the input of sim_03_medidas.R.
##
## Options (environment variables):
##   MZM_NSIM (default 100), MZM_BSIM (500), MZM_SIZES ("100,200,300"),
##   MZM_QUICK = 1 (only S1 with n = 100, for testing)
##   MZM_NCORES = 12 (number of parallel processes; default: logical cores - 1)
##   MZM_JOBS = "S1:300,S1:200" (run only these scenario:n combinations;
##              used to split the work among several computers)
## =============================================================================
source("sim_00_setup.R")
RNGkind("L'Ecuyer-CMRG")

sc         <- readRDS(file.path(OUT_DIR, "sim_scenarios.rds"))$scenarios
N_SIM      <- as.integer(Sys.getenv("MZM_NSIM", "100"))
B_SIM      <- as.integer(Sys.getenv("MZM_BSIM", "500"))
N_SIZES    <- as.numeric(strsplit(Sys.getenv("MZM_SIZES", "100,200,300"), ",")[[1]])
NCORES     <- as.integer(Sys.getenv("MZM_NCORES", max(1, detectCores() - 1)))   # e.g. MZM_NCORES = 12
MAIN_LIMIT <- 1800      # seconds for the main fit of a replication
BOOT_LIMIT <- 120       # seconds for each bootstrap refit
if (nzchar(Sys.getenv("MZM_QUICK"))) { sc <- sc["S1"]; N_SIZES <- 100 }
PARTS <- file.path(OUT_DIR, "parts")
dir.create(PARTS, showWarnings = FALSE)

## ---- Helpers (also exported to the workers) ---------------------------------
with_limit <- function(expr, sec) {
  setTimeLimit(elapsed = sec, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
  try(expr, silent = TRUE)
}
## Bootstrap refit: BFGS from the replication MLE; no slow fallback
fit_reg_boot <- function(Y, start, x) {
  xmat <- cbind(1, x)
  par0 <- c(stats::qlogis(start["w0"]), start[-1])
  fg <- MultDiscreteZM:::zmreg_fg(xmat, Y[, 1], Y[, 2], Y[, 3])   # analytic gradient
  o <- stats::optim(par0, fg$fn, fg$gr, method = "BFGS",
                    control = list(maxit = 500, reltol = 1e-10))
  out <- o$par; out[1] <- stats::plogis(out[1]); names(out) <- names(start)
  attr(out, "converged") <- (o$convergence == 0)
  out
}

## ---- One replication (one task) ---------------------------------------------
sim_task <- function(task) {
  if (file.exists(task$file)) return(TRUE)
  x <- task$xcov
  assign(".Random.seed", task$seed, envir = .GlobalEnv)
  gen <- boot_funs_reg(x)$sim
  Y   <- gen(task$true)
  t0  <- proc.time()[["elapsed"]]
  est <- with_limit(fitreg3zmcmpuI(Y[, 1], Y[, 2], Y[, 3], x), MAIN_LIMIT)
  if (inherits(est, "try-error") || any(!is.finite(est))) {
    res <- list(r = task$r, ran = FALSE, converged = FALSE,
                error = if (inherits(est, "try-error")) as.character(est) else "non-finite",
                seconds = proc.time()[["elapsed"]] - t0)
  } else {
    bt <- matrix(NA_real_, task$B, length(est), dimnames = list(NULL, names(est)))
    n_timeout <- 0
    for (b in seq_len(task$B)) {
      fb <- with_limit(fit_reg_boot(gen(est), est, x), BOOT_LIMIT)
      if (inherits(fb, "try-error")) { n_timeout <- n_timeout + 1; next }
      if (all(is.finite(fb)) && isTRUE(attr(fb, "converged"))) bt[b, ] <- fb
    }
    ci <- apply(bt, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    res <- list(r = task$r, ran = TRUE, converged = isTRUE(attr(est, "converged")),
                est = as.numeric(est), low = ci[1, ], up = ci[2, ],
                boot_conv = mean(stats::complete.cases(bt)), boot_timeout = n_timeout,
                seconds = proc.time()[["elapsed"]] - t0)
  }
  saveRDS(res, task$file)
  TRUE
}

## ---- Build the list of tasks (all replications not yet saved) ---------------
tasks <- list(); combos <- list()
jobs <- expand.grid(s = names(sc), n = N_SIZES, stringsAsFactors = FALSE)
if (nzchar(Sys.getenv("MZM_JOBS"))) {
  jj   <- strsplit(strsplit(Sys.getenv("MZM_JOBS"), ",")[[1]], ":")
  jobs <- data.frame(s = trimws(sapply(jj, `[`, 1)), n = as.numeric(sapply(jj, `[`, 2)),
                     stringsAsFactors = FALSE)
}
for (k in seq_len(nrow(jobs))) {
  s <- jobs$s[k]; n <- jobs$n[k]
  final <- file.path(OUT_DIR, sprintf("sim_%s_n%d.rds", s, n))
  set.seed(SEED + n);  x <- runif(n)                   # fixed covariate
  set.seed(SEED + 1000 * match(s, c("S1", "S2", "S3")) + n)   # streams for this (s, n)
  seeds <- vector("list", N_SIM); seeds[[1]] <- .Random.seed
  if (N_SIM > 1) for (r in 2:N_SIM) seeds[[r]] <- nextRNGStream(seeds[[r - 1]])
  files <- file.path(PARTS, sprintf("%s_n%d_r%03d.rds", s, n, seq_len(N_SIM)))
  combos[[length(combos) + 1]] <- list(s = s, n = n, x = x, files = files, final = final)
  if (file.exists(final)) next
  for (r in seq_len(N_SIM)) if (!file.exists(files[r]))
    tasks[[length(tasks) + 1]] <- list(r = r, seed = seeds[[r]], true = sc[[s]]$true,
                                       xcov = x, B = B_SIM, file = files[r])
}
cat(sprintf("%d replications to run on %d cores. Progress: count the files in\n  %s\n",
            length(tasks), NCORES, PARTS))
cat("(in another R window: length(list.files(\"", PARTS, "\")) )\n", sep = "")

## ---- Run --------------------------------------------------------------------
if (length(tasks)) {
  t0 <- Sys.time()
  cl <- makeCluster(NCORES)
  invisible(clusterEvalQ(cl, suppressMessages(library(MultDiscreteZM))))
  invisible(clusterCall(cl, function(f) { source(f); NULL }, UTILS_FILE))
  clusterExport(cl, c("with_limit", "fit_reg_boot", "MAIN_LIMIT", "BOOT_LIMIT"))
  invisible(parLapplyLB(cl, tasks, sim_task, chunk.size = 1))
  stopCluster(cl)
  cat(sprintf("Finished in %.1f hours\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}

## ---- Combine the replications of each (scenario, n) -------------------------
for (cb in combos) {
  if (file.exists(cb$final) || !all(file.exists(cb$files))) next
  parts <- lapply(cb$files, readRDS)
  ran   <- vapply(parts, `[[`, logical(1), "ran")
  nm    <- names(sc[[cb$s]]$true)
  pick  <- function(what) { m <- do.call(rbind, lapply(parts[ran], `[[`, what)); colnames(m) <- nm; m }
  res <- list(scenario = cb$s, n = cb$n, true = sc[[cb$s]]$true, x = cb$x,
              N = N_SIM, B = B_SIM,
              est = pick("est"), low = pick("low"), up = pick("up"),
              converged    = vapply(parts[ran], `[[`, logical(1), "converged"),
              n_error      = sum(!ran),
              boot_conv    = vapply(parts[ran], `[[`, numeric(1), "boot_conv"),
              boot_timeout = vapply(parts[ran], `[[`, numeric(1), "boot_timeout"),
              seconds      = vapply(parts, `[[`, numeric(1), "seconds"),
              minutes      = sum(vapply(parts, `[[`, numeric(1), "seconds")) / 60)
  saveRDS(res, cb$final)
  cat(sprintf("saved %s | convergence rate %.2f | median %.1f min per replication\n",
              basename(cb$final), sum(res$converged) / N_SIM, median(res$seconds) / 60))
}

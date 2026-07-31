#' Random Generation for the Trivariate Zero-Modified COM-Poisson (Original Parametrization)
#'
#' @description Generates random vectors from the trivariate Zero-Modified
#' COM-Poisson (MZMCMP) distribution without covariates, based on the original lambda-nu parametrization.
#'
#' @param n Integer. The number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param lambda1,lambda2,lambda3 Numeric vectors. The rate/intensity parameters for each margin.
#' @param nu1,nu2,nu3 Numeric vectors. The dispersion parameters for each margin.
#' @param p Numeric vector. The zero-modification parameter.
#'
#' @return A data frame with \code{n} rows and 3 columns (\code{Y1}, \code{Y2}, \code{Y3}) containing the simulated values.
#' @author Rogerio A. Santana
#'
#' @examples
#' set.seed(123)
#' # Simulating 10 observations with scalar parameters
#' sim_data <- r3zmcmpI(n = 10, lambda1 = 3, lambda2 = 2, lambda3 = 4,
#'                      nu1 = 0.8, nu2 = 1.2, nu3 = 0.9, p = 0.7)
#' head(sim_data)
#'
#' @export
r3zmcmpI <- function(n, lambda1, lambda2, lambda3, nu1, nu2, nu3, p) {

  # --- CRAN Standard for n ---
  if (length(n) > 1) n <- length(n)

  # --- Parameter Recycling ---
  # Ensures all inputs have length 'n' to avoid mismatch errors
  lambda1 <- rep_len(lambda1, n); nu1 <- rep_len(nu1, n)
  lambda2 <- rep_len(lambda2, n); nu2 <- rep_len(nu2, n)
  lambda3 <- rep_len(lambda3, n); nu3 <- rep_len(nu3, n)
  p <- rep_len(p, n)

  # --- Pre-calculate Joint Probability of Zero ---
  # Computes the marginal probabilities of zero using the univariate function
  d0_1 <- dcmp_sp(y = rep(0, n), lambda = lambda1, nu = nu1)
  d0_2 <- dcmp_sp(y = rep(0, n), lambda = lambda2, nu = nu2)
  d0_3 <- dcmp_sp(y = rep(0, n), lambda = lambda3, nu = nu3)

  d0_joint <- d0_1 * d0_2 * d0_3

  # --- Hurdle Probability (Omega) ---
  # Protection Layer: Ensures omega stays strictly bounded [0, 1].
  omega <- pmax(0, pmin(1, p * (1 - d0_joint)))

  # Initialize the output matrix
  Y <- matrix(0, nrow = n, ncol = 3)

  # Determine which observations cross the hurdle
  hurdle <- stats::rbinom(n, size = 1, prob = omega)
  simula <- which(hurdle == 1)
  n_sim <- length(simula)

  if (n_sim > 0) {
    # Extract parameters only for the observations that crossed the hurdle
    l1_s <- lambda1[simula]; n1_s <- nu1[simula]
    l2_s <- lambda2[simula]; n2_s <- nu2[simula]
    l3_s <- lambda3[simula]; n3_s <- nu3[simula]

    cand1 <- numeric(n_sim)
    cand2 <- numeric(n_sim)
    cand3 <- numeric(n_sim)

    # --- Vectorized Rejection Sampling ---
    # Boolean mask indicating which rows still need valid values
    needs_update <- rep(TRUE, n_sim)

    while (any(needs_update)) {
      idx <- which(needs_update)
      n_update <- length(idx)

      # Generate candidate values simultaneously in chunks
      cand1[idx] <- rcmp_sp(n_update, l1_s[idx], n1_s[idx])
      cand2[idx] <- rcmp_sp(n_update, l2_s[idx], n2_s[idx])
      cand3[idx] <- rcmp_sp(n_update, l3_s[idx], n3_s[idx])

      # Check truncation condition (sum strictly greater than 0)
      valid <- (cand1[idx] + cand2[idx] + cand3[idx]) > 0

      # Update the mask: valid candidates no longer need updating
      needs_update[idx[valid]] <- FALSE
    }

    # Assign validated candidates to the output matrix
    Y[simula, 1] <- cand1
    Y[simula, 2] <- cand2
    Y[simula, 3] <- cand3
  }

  return(data.frame(Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3]))
}

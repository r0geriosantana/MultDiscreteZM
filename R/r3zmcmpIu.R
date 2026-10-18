#' Random Generation for the Trivariate Zero-Modified COM-Poisson
#'
#' @description Generates random vectors from the trivariate Zero-Modified
#' COM-Poisson (MZMCMP) distribution without covariates, parameterized by the mean.
#'
#' @param n Integer. The number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu1,mu2,mu3 Numeric vectors. The mean parameters for each margin.
#' @param kappa1,kappa2,kappa3 Numeric vectors. The dispersion parameters on the log scale for each margin.
#' @param p Numeric vector. The zero-modification parameter.
#'
#' @return A data frame with \code{n} rows and 3 columns (\code{Y1}, \code{Y2}, \code{Y3}) containing the simulated values.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' set.seed(123)
#' # Simulating 10 observations with scalar parameters
#' sim_data <- r3zmcmpIu(n = 10, mu1 = 2, mu2 = 3, mu3 = 1.5,
#'                       kappa1 = -0.5, kappa2 = 0.2, kappa3 = -0.1, p = 0.8)
#' head(sim_data)
#'
#' @export
r3zmcmpIu <- function(n, mu1, mu2, mu3, kappa1, kappa2, kappa3, p) {

  # --- CRAN Standard for n ---
  if (length(n) > 1) n <- length(n)

  # --- Parameter Recycling ---
  # Ensures all inputs have length 'n' to avoid mismatch errors
  mu1 <- rep_len(mu1, n); kappa1 <- rep_len(kappa1, n)
  mu2 <- rep_len(mu2, n); kappa2 <- rep_len(kappa2, n)
  mu3 <- rep_len(mu3, n); kappa3 <- rep_len(kappa3, n)
  p <- rep_len(p, n)

  # --- Pre-calculate Joint Probability of Zero ---
  # Computes the marginal probabilities of zero using the univariate function
  d0_1 <- dcmp_u(y = rep(0, n), mu = mu1, kappa = kappa1)
  d0_2 <- dcmp_u(y = rep(0, n), mu = mu2, kappa = kappa2)
  d0_3 <- dcmp_u(y = rep(0, n), mu = mu3, kappa = kappa3)

  d0_joint <- d0_1 * d0_2 * d0_3

  # --- Hurdle Probability (Omega) ---
  # Protection Layer: Ensures omega stays strictly bounded [0, 1].
  # This prevents stats::rbinom failures during boundary testing in regression models.
  omega <- pmax(0, pmin(1, p * (1 - d0_joint)))

  # Initialize the output matrix
  Y <- matrix(0, nrow = n, ncol = 3)

  # Determine which observations cross the hurdle (structural vs non-zero states)
  hurdle <- stats::rbinom(n, size = 1, prob = omega)
  simula <- which(hurdle == 1)
  n_sim <- length(simula)

  if (n_sim > 0) {
    # Extract parameters only for the observations that crossed the hurdle
    m1_s <- mu1[simula]; k1_s <- kappa1[simula]
    m2_s <- mu2[simula]; k2_s <- kappa2[simula]
    m3_s <- mu3[simula]; k3_s <- kappa3[simula]

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
      cand1[idx] <- rcmp_u(n_update, m1_s[idx], k1_s[idx])
      cand2[idx] <- rcmp_u(n_update, m2_s[idx], k2_s[idx])
      cand3[idx] <- rcmp_u(n_update, m3_s[idx], k3_s[idx])

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

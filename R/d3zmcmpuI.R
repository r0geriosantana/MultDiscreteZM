#' Probability mass function of the Trivariate Zero-Modified COM-Poisson
#'
#' @description Calculates the probability mass function (PMF) of the trivariate
#' Zero-Modified COM-Poisson (MZMCMP) distribution parameterized by the mean.
#'
#' @param y1,y2,y3 Numeric vectors of non-negative integers (quantiles).
#' @param mu1,mu2,mu3 Numeric vectors. The mean parameters for each margin.
#' @param kappa1,kappa2,kappa3 Numeric vectors. The dispersion parameters on the log scale.
#' @param p Numeric vector. The zero-modification parameter.
#' @param log Logical; if TRUE, probabilities are returned as log-probabilities. Defaults to FALSE.
#'
#' @return A numeric vector of probabilities or log-probabilities.
#' @author Rogerio A. Santana
#'
#' @examples
#' d3zmcmpuI(y1 = c(0, 1), y2 = c(0, 2), y3 = c(0, 1),
#'           mu1 = 2, mu2 = 3, mu3 = 1.5,
#'           kappa1 = -0.5, kappa2 = 0.2, kappa3 = -0.1, p = 0.8)
#'
#' @export
d3zmcmpuI <- function(y1, y2, y3, mu1, mu2, mu3, kappa1, kappa2, kappa3, p, log = FALSE) {

  # --- CRAN Standard: Vector Recycling ---
  # Ensures all inputs (quantiles and parameters) have the same length
  n_obs <- max(length(y1), length(y2), length(y3),
               length(mu1), length(mu2), length(mu3),
               length(kappa1), length(kappa2), length(kappa3), length(p))

  y1 <- rep_len(y1, n_obs); y2 <- rep_len(y2, n_obs); y3 <- rep_len(y3, n_obs)
  mu1 <- rep_len(mu1, n_obs); mu2 <- rep_len(mu2, n_obs); mu3 <- rep_len(mu3, n_obs)
  kappa1 <- rep_len(kappa1, n_obs); kappa2 <- rep_len(kappa2, n_obs); kappa3 <- rep_len(kappa3, n_obs)
  p <- rep_len(p, n_obs)

  # --- Step 1: Calculate probabilities for each component ---
  # Joint probability of zero for the base COM-Poisson part
  dzero1 <- dcmp_u(rep(0, n_obs), mu1, kappa1)
  dzero2 <- dcmp_u(rep(0, n_obs), mu2, kappa2)
  dzero3 <- dcmp_u(rep(0, n_obs), mu3, kappa3)
  prob_z0 <- dzero1 * dzero2 * dzero3

  # --- Step 2: Calculate the base joint PMF for all observations ---
  # This computes P(Y1=y1, Y2=y2, Y3=y3) ignoring the global zero-modification
  term1 <- dcmp_u(y1, mu1, kappa1)
  term2 <- dcmp_u(y2, mu2, kappa2)
  term3 <- dcmp_u(y3, mu3, kappa3)
  prod_term <- term1 * term2 * term3

  # --- Step 3: Apply the global zero-inflation structure ---
  is_all_zeros <- (y1 == 0) & (y2 == 0) & (y3 == 0)

  densities <- numeric(n_obs)

  # For observations where all components are zeros:
  # P(Y = 0) = 1 - omega = 1 - p * (1 - P(Y_base = 0))
  pz0 <- p * (1 - prob_z0)
  densities[is_all_zeros] <- 1 - pz0[is_all_zeros]

  # For the remaining observations (at least one component is > 0):
  # P(Y = y) = omega * [f(y) / (1 - P(Y_base = 0))]
  # This simplifies mathematically to: p * f(y)
  densities[!is_all_zeros] <- p[!is_all_zeros] * prod_term[!is_all_zeros]

  # --- Step 4: Return conditional on the 'log' argument ---
  if (log) {
    return(log(densities))
  } else {
    return(densities)
  }
}

#' Probability mass function of the Trivariate Zero-Modified COM-Poisson (Original Parametrization)
#'
#' @description Calculates the probability mass function (PMF) of the trivariate
#' Zero-Modified COM-Poisson (MZMCMP) distribution based on the original lambda-nu parametrization.
#'
#' @param y1,y2,y3 Numeric vectors of non-negative integers (quantiles).
#' @param lambda1,lambda2,lambda3 Numeric vectors. The rate/intensity parameters for each margin.
#' @param nu1,nu2,nu3 Numeric vectors. The dispersion parameters for each margin.
#' @param p Numeric vector. The zero-modification parameter.
#' @param log Logical; if TRUE, probabilities are returned as log-probabilities. Defaults to FALSE.
#'
#' @return A numeric vector of probabilities or log-probabilities.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' d3zmcmpI(y1 = c(0, 1), y2 = c(0, 2), y3 = c(0, 1),
#'          lambda1 = 3, lambda2 = 2, lambda3 = 4,
#'          nu1 = 0.8, nu2 = 1.2, nu3 = 0.9, p = 0.7)
#'
#' @export
d3zmcmpI <- function(y1, y2, y3, lambda1, lambda2, lambda3, nu1, nu2, nu3, p, log = FALSE) {

  # --- CRAN Standard: Vector Recycling ---
  # Ensures all inputs (quantiles and parameters) have the same length
  n_obs <- max(length(y1), length(y2), length(y3),
               length(lambda1), length(lambda2), length(lambda3),
               length(nu1), length(nu2), length(nu3), length(p))

  y1 <- rep_len(y1, n_obs); y2 <- rep_len(y2, n_obs); y3 <- rep_len(y3, n_obs)
  lambda1 <- rep_len(lambda1, n_obs); lambda2 <- rep_len(lambda2, n_obs); lambda3 <- rep_len(lambda3, n_obs)
  nu1 <- rep_len(nu1, n_obs); nu2 <- rep_len(nu2, n_obs); nu3 <- rep_len(nu3, n_obs)
  p <- rep_len(p, n_obs)

  # --- Step 1: Calculate probabilities for each component ---
  # Joint probability of zero for the base COM-Poisson part
  dzero1 <- dcmp_sp(rep(0, n_obs), lambda1, nu1)
  dzero2 <- dcmp_sp(rep(0, n_obs), lambda2, nu2)
  dzero3 <- dcmp_sp(rep(0, n_obs), lambda3, nu3)
  prob_z0 <- dzero1 * dzero2 * dzero3

  # --- Step 2: Calculate the base joint PMF for all observations ---
  # This computes P(Y1=y1, Y2=y2, Y3=y3) ignoring the global zero-modification
  term1 <- dcmp_sp(y1, lambda1, nu1)
  term2 <- dcmp_sp(y2, lambda2, nu2)
  term3 <- dcmp_sp(y3, lambda3, nu3)
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

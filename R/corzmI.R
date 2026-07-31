#' Theoretical Correlation of the Zero-Modified COM-Poisson (Type I)
#'
#' @description Internal helper function to calculate the theoretical correlation
#' between two margins of the MZMCMP distribution.
#'
#' @param p Numeric vector. The zero-modification parameter.
#' @param lambda1 Numeric vector. The rate parameter for the first margin.
#' @param lambda2 Numeric vector. The rate parameter for the second margin.
#' @param nu1 Numeric vector. The dispersion parameter for the first margin.
#' @param nu2 Numeric vector. The dispersion parameter for the second margin.
#'
#' @return A numeric vector representing the theoretical correlation.
#' @export
corzmI <- function(p, lambda1, lambda2, nu1, nu2) {

  # --- Vector Recycling ---
  # Ensures all parameters have the same length
  n <- max(length(p), length(lambda1), length(lambda2), length(nu1), length(nu2))
  p <- rep_len(p, n)
  lambda1 <- rep_len(lambda1, n); lambda2 <- rep_len(lambda2, n)
  nu1 <- rep_len(nu1, n); nu2 <- rep_len(nu2, n)

  # --- Base Expectations ---
  # E[X] of the base COM-Poisson
  mu1 <- EXcmp(lambda1, nu1)
  mu2 <- EXcmp(lambda2, nu2)

  # --- Base Second Factorial Moments E[X(X-1)] ---
  dz2_1 <- dz2cmp(lambda1, nu1)
  dz2_2 <- dz2cmp(lambda2, nu2)

  # --- Variances of the Zero-Modified Margins ---
  # E_zm[X^2] = p * (E[X(X-1)] + E[X])
  E2X_1 <- dz2_1 + mu1
  var1 <- (p * E2X_1) - (p * mu1)^2

  E2X_2 <- dz2_2 + mu2
  var2 <- (p * E2X_2) - (p * mu2)^2

  # --- Covariance ---
  # Due to the Type I hurdle structure, covariance arises purely from the shared zero state
  cov12 <- p * (1 - p) * mu1 * mu2

  # --- Correlation ---
  # Protect against division by zero
  denom <- sqrt(var1) * sqrt(var2)

  # Calculates correlation, returning 0 if the variance is effectively zero
  corr <- ifelse(denom > 1e-15, cov12 / denom, 0)

  return(corr)
}

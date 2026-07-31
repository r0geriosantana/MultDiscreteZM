#' Theoretical Correlation for the Trivariate ZMNB Model
#'
#' @description Calculates the theoretical population correlation between two
#' variables in a Zero-Modified Negative Binomial (ZMNB) distribution based on
#' their marginal parameters and the shared zero-modification parameter.
#'
#' @param p Numeric vector or scalar. The zero-modification parameter (Probability of Y > 0).
#' @param mu1 Numeric vector or scalar. The mean parameter of the first variable.
#' @param mu2 Numeric vector or scalar. The mean parameter of the second variable.
#' @param theta1 Numeric vector or scalar. The dispersion parameter (size) of the first variable.
#' @param theta2 Numeric vector or scalar. The dispersion parameter (size) of the second variable.
#'
#' @return A numeric vector representing the theoretical correlation.
#' @export
corzmnbI <- function(p, mu1, mu2, theta1, theta2) {

  # --- 1. Expected Values (First Moments) ---
  E1 <- p * mu1
  E2 <- p * mu2

  # Joint expectation assuming conditional independence given Y > 0
  E12 <- p * mu1 * mu2

  # --- 2. Variances (Second Central Moments) ---
  # Internal helper function for marginal variance
  calc_variance <- function(p_val, mu_val, theta_val) {
    E_mean <- p_val * mu_val

    # Second uncentered moment of the non-zero component
    M2 <- mu_val + (mu_val^2) * (1 + 1 / theta_val)

    # Law of total variance / Unconditional variance
    return((p_val * M2) - (E_mean)^2)
  }

  V1 <- calc_variance(p, mu1, theta1)
  V2 <- calc_variance(p, mu2, theta2)

  # --- 3. Covariance and Correlation ---
  cov12 <- E12 - (E1 * E2)

  # Numerical protection against division by zero in highly degenerate simulated cases
  denominator <- sqrt(V1 * V2)
  correlation <- ifelse(denominator <= 1e-15, 0, cov12 / denominator)

  return(correlation)
}

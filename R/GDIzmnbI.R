#' Generalized Dispersion Index for Trivariate ZMNB Type I
#'
#' @description Computes the Generalized Dispersion Index (GDI) for the
#' trivariate Zero-Modified Negative Binomial (Type I) distribution.
#' A GDI > 1 indicates overdispersion relative to the multivariate Poisson.
#'
#' @param p Numeric vector or scalar. The zero-modification parameter (Probability of Y > 0).
#' @param mu1,mu2,mu3 Numeric vectors or scalars. Mean parameters for each component.
#' @param theta1,theta2,theta3 Numeric vectors or scalars. Dispersion parameters (size).
#'
#' @return A numeric vector representing the Generalized Dispersion Index.
#' @export
GDIzmnbI <- function(p, mu1, mu2, mu3, theta1, theta2, theta3) {

  # --- 1. Expected Values (First Moments) ---
  E1 <- p * mu1
  E2 <- p * mu2
  E3 <- p * mu3

  # --- 2. Variances ---
  # Internal helper function for marginal variance
  calc_variance <- function(p_val, mu_val, theta_val) {
    E_mean <- p_val * mu_val
    M2 <- mu_val + (mu_val^2) * (1 + 1 / theta_val)
    return((p_val * M2) - E_mean^2)
  }

  var1 <- calc_variance(p, mu1, theta1)
  var2 <- calc_variance(p, mu2, theta2)
  var3 <- calc_variance(p, mu3, theta3)

  # --- 3. Covariances (Due to the Hurdle/Zero-Modified structure) ---
  cov12 <- p * (1 - p) * mu1 * mu2
  cov13 <- p * (1 - p) * mu1 * mu3
  cov23 <- p * (1 - p) * mu2 * mu3

  # --- 4. Algebraic Vectorization of the GDI ---
  # Expanding the quadratic form: Er %*% M %*% t(Er) analytically.
  # This avoids matrix multiplication loops and natively supports long vectors.
  numerador <- (E1 * var1) + (E2 * var2) + (E3 * var3) +
    (2 * sqrt(E1 * E2) * cov12) +
    (2 * sqrt(E1 * E3) * cov13) +
    (2 * sqrt(E2 * E3) * cov23)

  denominador <- (E1^2) + (E2^2) + (E3^2)

  # Numerical protection: If denominator is practically zero, return NA to avoid Inf
  gdi <- ifelse(denominador <= 1e-15, NA, numerador / denominador)

  return(as.numeric(gdi))
}

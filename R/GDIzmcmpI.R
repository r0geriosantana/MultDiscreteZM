#' Generalized Dispersion Index (GDI) for Trivariate MZMCMP (Type I)
#'
#' @description Internal helper function to calculate the Generalized Dispersion Index.
#'
#' @param p Numeric vector. The zero-modification parameter.
#' @param lambda1,lambda2,lambda3 Numeric vectors. Rate parameters for each margin.
#' @param nu1,nu2,nu3 Numeric vectors. Dispersion parameters for each margin.
#'
#' @return A numeric vector representing the dispersion index.
#' @export
GDIzmcmpI <- function(p, lambda1, lambda2, lambda3, nu1, nu2, nu3) {

  # --- Vector Recycling ---
  n <- max(length(p), length(lambda1), length(lambda2), length(lambda3),
           length(nu1), length(nu2), length(nu3))

  p <- rep_len(p, n)
  lambda1 <- rep_len(lambda1, n); lambda2 <- rep_len(lambda2, n); lambda3 <- rep_len(lambda3, n)
  nu1 <- rep_len(nu1, n); nu2 <- rep_len(nu2, n); nu3 <- rep_len(nu3, n)

  # --- Base Expectations ---
  mu1 <- EXcmp(lambda1, nu1)
  mu2 <- EXcmp(lambda2, nu2)
  mu3 <- EXcmp(lambda3, nu3)

  # --- Modified Expectations ---
  E1 <- p * mu1
  E2 <- p * mu2
  E3 <- p * mu3

  # --- Variances ---
  # Reusing the internal helper dz2cmp to avoid duplicate computations
  var1 <- (p * (dz2cmp(lambda1, nu1) + mu1)) - E1^2
  var2 <- (p * (dz2cmp(lambda2, nu2) + mu2)) - E2^2
  var3 <- (p * (dz2cmp(lambda3, nu3) + mu3)) - E3^2

  # --- Covariances (Type I hurdle structure) ---
  cov12 <- p * (1 - p) * mu1 * mu2
  cov13 <- p * (1 - p) * mu1 * mu3
  cov23 <- p * (1 - p) * mu2 * mu3

  # --- Vectorized GDI Calculation (Algebraic Expansion) ---
  # Numerator: expansion of the quadratic form Er %*% M %*% t(Er)
  numerador <- (E1 * var1) + (E2 * var2) + (E3 * var3) +
    (2 * sqrt(E1 * E2) * cov12) +
    (2 * sqrt(E1 * E3) * cov13) +
    (2 * sqrt(E2 * E3) * cov23)

  # Denominator: expansion of E %*% t(E)
  denominador <- E1^2 + E2^2 + E3^2

  # --- Protection Layer ---
  # Prevents division by zero if all expectations are zero
  gdi <- ifelse(denominador > 1e-15, numerador / denominador, NA_real_)

  return(gdi)
}

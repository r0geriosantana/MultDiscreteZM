#' Probability Mass Function for the Trivariate ZMNB Distribution
#'
#' @description Computes the joint probability mass function of the
#' trivariate Zero-Modified Negative Binomial distribution.
#'
#' @param y1,y2,y3 Numeric vectors of count responses.
#' @param mu1,mu2,mu3 Numeric vectors of mean parameters.
#' @param theta1,theta2,theta3 Numeric vectors of dispersion parameters (size).
#' @param p Numeric vector of the zero-modification parameter (probability of Y > 0).
#'
#' @return A numeric vector of probabilities.
#' @export
d3zmnbI <- function(y1, y2, y3, mu1, mu2, mu3, theta1, theta2, theta3, p) {

  # --- Vector size verification ---
  n_obs <- length(y1)

  # Ensure all parameters are vectors of the same length or scalars
  if (length(mu1) == 1) mu1 <- rep(mu1, n_obs)
  if (length(mu2) == 1) mu2 <- rep(mu2, n_obs)
  if (length(mu3) == 1) mu3 <- rep(mu3, n_obs)
  if (length(theta1) == 1) theta1 <- rep(theta1, n_obs)
  if (length(theta2) == 1) theta2 <- rep(theta2, n_obs)
  if (length(theta3) == 1) theta3 <- rep(theta3, n_obs)
  if (length(p) == 1) p <- rep(p, n_obs)

  # --- Step 1: Calculate the baseline zero probabilities ---
  # dnbinom is inherently vectorized and handles 0 perfectly
  prob_z1 <- stats::dnbinom(0, mu = mu1, size = theta1)
  prob_z2 <- stats::dnbinom(0, mu = mu2, size = theta2)
  prob_z3 <- stats::dnbinom(0, mu = mu3, size = theta3)

  prob_z0 <- prob_z1 * prob_z2 * prob_z3

  # --- Step 2: Calculate the baseline joint PMF for the actual observations ---
  prod_term <- stats::dnbinom(y1, mu = mu1, size = theta1) * stats::dnbinom(y2, mu = mu2, size = theta2) * stats::dnbinom(y3, mu = mu3, size = theta3)

  # --- Step 3: Apply the Zero-Modified structure ---
  is_all_zeros <- (y1 == 0) & (y2 == 0) & (y3 == 0)
  densities <- numeric(n_obs)

  # For observations where all components are simultaneously zero
  densities[is_all_zeros] <- 1 - p[is_all_zeros]

  # For all other observations (at least one component > 0)
  denominador <- 1 - prob_z0[!is_all_zeros]

  # Extreme numerical protection (prevents division by zero if prob_z0 approaches 1)
  denominador <- pmax(denominador, 1e-15)

  densities[!is_all_zeros] <- p[!is_all_zeros] * (prod_term[!is_all_zeros] / denominador)

  return(densities)
}

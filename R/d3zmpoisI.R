#' Probability Mass Function for the Trivariate ZM-Poisson Distribution
#'
#' @description Computes the joint probability mass function of the
#' trivariate Zero-Modified Poisson Type I distribution.
#'
#' @param y1,y2,y3 Numeric vectors of count responses.
#' @param mu1,mu2,mu3 Numeric vectors of mean parameters.
#' @param p Numeric vector of the zero-modification parameter (probability of Y > 0).
#'
#' @return A numeric vector of probabilities.
#' @export
d3zmpoisI <- function(y1, y2, y3, mu1, mu2, mu3, p) {

  # --- Vector size verification ---
  n_obs <- length(y1)

  # Ensure all parameters are vectors of the same length or scalars
  if (length(mu1) == 1) mu1 <- rep(mu1, n_obs)
  if (length(mu2) == 1) mu2 <- rep(mu2, n_obs)
  if (length(mu3) == 1) mu3 <- rep(mu3, n_obs)
  if (length(p) == 1) p <- rep(p, n_obs)

  # --- Step 1: Calculate the baseline zero probabilities ---
  # dpois is inherently vectorized and handles 0 perfectly
  prob_z1 <- stats::dpois(0, lambda = mu1)
  prob_z2 <- stats::dpois(0, lambda = mu2)
  prob_z3 <- stats::dpois(0, lambda = mu3)

  prob_z0 <- prob_z1 * prob_z2 * prob_z3

  # --- Step 2: Calculate the baseline joint PMF for the actual observations ---
  prod_term <- stats::dpois(y1, lambda = mu1) * stats::dpois(y2, lambda = mu2) * stats::dpois(y3, lambda = mu3)

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

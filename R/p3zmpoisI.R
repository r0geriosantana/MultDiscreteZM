#' Cumulative Distribution Function for the Trivariate ZM-Poisson Distribution
#'
#' @description Computes the joint cumulative distribution function (CDF) of the
#' trivariate Zero-Modified Poisson Type I distribution analytically.
#'
#' @param y1,y2,y3 Numeric vectors of count quantiles.
#' @param mu1,mu2,mu3 Numeric vectors of mean parameters.
#' @param p Numeric vector of the zero-modification parameter (probability of Y > 0).
#'
#' @return A numeric vector of cumulative probabilities.
#' @export
p3zmpoisI <- function(y1, y2, y3, mu1, mu2, mu3, p) {

  n_obs <- max(length(y1), length(y2), length(y3))

  # --- Vector size verification (Crucial for Regression loops) ---
  if (length(mu1) == 1) mu1 <- rep(mu1, n_obs)
  if (length(mu2) == 1) mu2 <- rep(mu2, n_obs)
  if (length(mu3) == 1) mu3 <- rep(mu3, n_obs)
  if (length(p) == 1) p <- rep(p, n_obs)

  # --- Step 1: Base case protection ---
  valid_obs <- (y1 >= 0) & (y2 >= 0) & (y3 >= 0)
  cdf_values <- numeric(n_obs)

  if (!any(valid_obs)) {
    return(cdf_values)
  }

  # --- Step 2: Calculate baseline PMF at zero ---
  prob_z1 <- stats::dpois(0, lambda = mu1)
  prob_z2 <- stats::dpois(0, lambda = mu2)
  prob_z3 <- stats::dpois(0, lambda = mu3)
  prob_z0 <- prob_z1 * prob_z2 * prob_z3

  # --- Step 3: Calculate baseline joint CDF ---
  # CORREÇÃO APLICADA: Uso de ppois (Acumulada) ao invés de dpois (Densidade)
  cdf_base1 <- stats::ppois(y1, lambda = mu1)
  cdf_base2 <- stats::ppois(y2, lambda = mu2)
  cdf_base3 <- stats::ppois(y3, lambda = mu3)
  cdf_base_joint <- cdf_base1 * cdf_base2 * cdf_base3

  # --- Step 4: Apply Hurdle transformation ---
  denominador <- pmax(1 - prob_z0, 1e-15)
  term_prob <- p * ((cdf_base_joint - prob_z0) / denominador)

  # Aplicar o filtro [valid_obs] também no parâmetro 'p'
  cdf_values[valid_obs] <- (1 - p[valid_obs]) + term_prob[valid_obs]

  cdf_values <- pmax(0, pmin(1, cdf_values))

  return(cdf_values)
}

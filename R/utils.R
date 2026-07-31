#' Expected Value of the Original COM-Poisson Distribution
#'
#' @description Internal helper function to calculate E(X).
#'
#' @param lambda Numeric vector. The rate/intensity parameter.
#' @param nu Numeric vector. The dispersion parameter.
#'
#' @keywords internal
#' @noRd
EXcmp <- function(lambda, nu) {

  # Determine max length for vector recycling support
  n <- max(length(lambda), length(nu))
  lambda <- rep_len(lambda, n)
  nu <- rep_len(nu, n)

  # --- Dynamic Truncation Heuristic ---
  # Prevents truncation errors for large lambda and small nu
  approx_mu <- max(lambda)^(1 / min(nu))
  proxy_sd <- sqrt(approx_mu / min(nu))
  max_iter <- ceiling(max(100, approx_mu + 20 * proxy_sd))

  j <- 0:max_iter

  # 1. Log-terms matrix: row = observation, col = j
  # log_term = j * log(lambda) - nu * log(j!)
  log_terms <- outer(log(lambda), j, "*") - outer(nu, lgamma(j + 1), "*")

  # 2. Numerical Stability (Log-Sum-Exp Trick)
  max_log <- apply(log_terms, 1, max)
  terms <- exp(log_terms - max_log)

  # 3. Calculation of the Mean E(X)
  numerator <- rowSums(sweep(terms, 2, j, "*"))
  denominator <- rowSums(terms)

  return(numerator / denominator)
}

#' Second Factorial Moment of the Original COM-Poisson
#'
#' @description Internal helper function to calculate E[X(X-1)].
#' Used to derive the theoretical variance.
#'
#' @param lambda Numeric vector. The rate/intensity parameter.
#' @param nu Numeric vector. The dispersion parameter.
#'
#' @keywords internal
#' @noRd
dz2cmp <- function(lambda, nu) {

  # Determine max length for vector recycling support
  n <- max(length(lambda), length(nu))
  lambda <- rep_len(lambda, n)
  nu <- rep_len(nu, n)

  # --- Dynamic Truncation Heuristic ---
  approx_mu <- max(lambda)^(1 / min(nu))
  proxy_sd <- sqrt(approx_mu / min(nu))
  max_iter <- ceiling(max(100, approx_mu + 20 * proxy_sd))

  j <- 0:max_iter

  # 1. Log-terms matrix
  log_terms <- outer(log(lambda), j, "*") - outer(nu, lgamma(j + 1), "*")

  # 2. Numerical Stability (Log-Sum-Exp Trick)
  max_log <- apply(log_terms, 1, max)
  terms <- exp(log_terms - max_log)

  # 3. Calculation of E[X(X-1)]
  denom <- rowSums(terms)

  # numerator: sum of j*(j-1) * exp(log_term)
  j_factorial_2 <- j * (j - 1)
  num <- rowSums(sweep(terms, 2, j_factorial_2, "*"))

  return(num / denom)
}

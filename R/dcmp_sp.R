#' Probability mass function of the original COM-Poisson distribution.
#'
#' @description Calculate the probability mass function (PMF) of the univariate
#' COM-Poisson distribution in its original parametrization.
#'
#' @param y A numeric vector of non-negative integers (quantiles).
#' @param lambda A numeric vector. The rate/intensity parameter of the distribution (lambda > 0).
#' @param nu A numeric vector. The dispersion parameter (nu >= 0).
#' @param log Logical; if TRUE, probabilities p are given as log(p). Defaults to FALSE.
#'
#' @details
#' The function computes the density using the normalization constant \code{\link{Zconst}}.
#' It includes vector recycling to align parameters and safeguards against numerical instability.
#'
#' @return A numeric vector of probabilities or log-probabilities.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' # Probabilidade pontual
#' dcmp_sp(y = 2, lambda = 3, nu = 0.76)
#'
#' # Log-probabilidade (ideal para MLE)
#' dcmp_sp(y = 2, lambda = 3, nu = 0.76, log = TRUE)
#'
#' @export
dcmp_sp <- function(y, lambda, nu, log = FALSE) {

  # --- Vector Recycling  ---
  # Ensures that y, lambda, and nu have the same length for C++.
  n_max <- max(length(y), length(lambda), length(nu))
  y <- rep_len(y, n_max)
  lambda <- rep_len(lambda, n_max)
  nu <- rep_len(nu, n_max)

  # 1. Numerical Protection
  lambda_safe <- pmax(lambda, 1e-10)

  # 2. Components of Log-Likelihood
  # Correction: log(lambda^y) = y * log(lambda)
  parte1 <- y * log(lambda_safe)
  parte2 <- -nu * lgamma(y + 1)

  # 3. C++ call for the constant Z
  log_Z <- log(Zconst(lambda_safe, nu))

  log_res <- parte1 + parte2 - log_Z

  # 4. Conditional return
  if (log) {
    return(log_res)
  } else {
    return(exp(log_res))
  }
}

#' Probability mass function of the univariate COM-Poisson distribution parameterized by the mean.
#'
#' @description Calculate the probability mass function (PMF) of the univariate COM-Poisson distribution
#' parameterized by the mean, as proposed by Ribeiro Jr et al. (2020).
#'
#' @param y A numeric vector of non-negative integers (quantiles).
#' @param mu A numeric vector. The mean parameter of the distribution (mu > 0).
#' @param kappa A numeric vector. The dispersion parameter on the log scale (nu = exp(kappa)).
#' @param log Logical; if TRUE, probabilities p are given as log(p). Defaults to FALSE.
#'
#' @details
#' The function computes the density using the normalization constant \code{\link{Zconst}}.
#' It includes safeguards against numerical instability in the lambda approximation.
#'
#' @return A numeric vector of probabilities or log-probabilities.
#' @author Rogerio A. Santana
#'
#' @examples
#' # Probabilidade pontual
#' dcmp_u(y = 2, mu = 2, kappa = -0.36)
#'
#' # Log-probabilidade (ideal para MLE)
#' dcmp_u(y = 2, mu = 2, kappa = -0.36, log = TRUE)
#'
#' @export
dcmp_u <- function(y, mu, kappa, log = FALSE) {

  # 1. Basic definitions
  nu <- exp(kappa)

  # 2. Calculation of Lambda (with the same protection against NaN used in Zconst)
  base_term <- mu + (nu - 1) / (2 * nu)

  # Protection: If the base is <= 0, we use an infinitesimal value to allow for logging.
  lambda_base <- pmax(base_term, 1e-10)

  # Note: log(lambda^y) = y * log(lambda) = y * log(base^nu) = y * nu * log(base)
  #We calculate log(lambda) directly for stability.
  log_lambda <- nu * log(lambda_base)

  # 3. Calculation of the components of the Log-Likelihood
  # Log(P) = y*log(lambda) - nu*log(y!) - log(Z)
  parte1 <- y * log_lambda
  parte2 <- -nu * lgamma(y + 1)

  # Call to the Zuconst function (which should already be in the package)
  #  Importante: Zuconst returns Z, so we apply log(Z).
  log_Z <- log(Zuconst(mu, kappa))

  log_res <- parte1 + parte2 - log_Z

  # 4. Conditional return
  if (log) {
    return(log_res)
  } else {
    return(exp(log_res))
  }
}

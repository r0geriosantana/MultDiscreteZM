#' Random Generation for the COM-Poisson with the mean-parametrization
#'
#' @description Generates random values from the COM-Poisson (CMP)
#' distribution based on the mean-parametrization.
#'
#' @param n Integer. The number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu A numeric vector. The mean parameter of the distribution.
#' @param kappa A numeric vector. The dispersion parameter on the log scale.
#'
#' @return A vector of random integers.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' # 1. Simple Simulation (scalar parameters)
#' set.seed(123)
#' rcmp_u(n = 10, mu = 3, kappa = -0.37)
#'
#' @export
rcmp_u <- function(n, mu, kappa) {

  # --- CRAN Standard for n ---
  if (length(n) > 1) n <- length(n)

  # --- Parameter Recycling ---
  if (length(mu) != n) mu <- rep_len(mu, n)
  if (length(kappa) != n) kappa <- rep_len(kappa, n)

  # Pre-calculate Z (using the new Rcpp function Zuconst)
  Z_vec <- Zuconst(mu, kappa)
  nu_vec <- exp(kappa)
  lambda_vec <- lambda_u(mu, kappa)   # exact-mean rate parameters

  rands <- numeric(n)

  for(i in 1:n) {
    z_i <- Z_vec[i]

    # --- PROTECTION LAYER 1: Z Verification ---
    if(is.na(z_i) || !is.finite(z_i) || z_i <= 0) {
      rands[i] <- NA
      next
    }

    # Calculate base lambda dynamically for recurrence
    nu_i <- nu_vec[i]
    lambda_i <- lambda_vec[i]

    U <- stats::runif(1)
    y <- 0

    # Fast calculation of P(Y=0) using pre-calculated Z
    p <- 1 / z_i
    if(is.na(p)) p <- 0

    Fp <- p

    # --- PROTECTION LAYER 2: The While Loop ---
    while(!is.na(Fp) && U > Fp) {
      y <- y + 1

      # Fast recurrence relation: P(Y=y) = P(Y=y-1) * (lambda / y^nu)
      p <- p * (lambda_i / (y^nu_i))

      if(is.na(p)) p <- 0

      Fp <- Fp + p

      # Protection against infinite loop (very long tail)
      if(y > 10000) break
    }

    rands[i] <- y
  }

  return(rands)
}

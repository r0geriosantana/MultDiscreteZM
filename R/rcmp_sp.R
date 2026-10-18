#' Random Generation for the original COM-Poisson distribution.
#'
#' @description Generates random values from the COM-Poisson (CMP)
#' distribution based on the original lambda-nu parametrization.
#'
#' @param n Integer. The number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param lambda A numeric vector. The rate/intensity parameter of the distribution.
#' @param nu A numeric vector. The dispersion parameter.
#'
#' @return A vector of random integers.
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' # 1. Simple Simulation (scalar parameters)
#' set.seed(123)
#' rcmp_sp(n = 10, lambda = 3, nu = 0.76)
#'
#' @export
rcmp_sp <- function(n, lambda, nu) {

  # --- CRAN standard for size 'n' ---
  if (length(n) > 1) n <- length(n)

  # --- Parameter Recycling ---
  if (length(lambda) != n) lambda <- rep_len(lambda, n)
  if (length(nu) != n) nu <- rep_len(nu, n)

  # Pre-calculate Z (using the C++ function)
  Z_vec <- Zconst(lambda, nu)

  rands <- numeric(n)

  for(i in 1:n) {
    z_i <- Z_vec[i]
    lambda_i <- lambda[i]
    nu_i <- nu[i]

    # --- LAYER OF PROTECTION 1: Z Validation ---
    if(is.na(z_i) || !is.finite(z_i) || z_i <= 0) {
      rands[i] <- NA
      next
    }

    U <- stats::runif(1)
    y <- 0

    # Ultra-fast calculation of P(Y=0)
    p <- 1 / z_i
    if(is.na(p)) p <- 0

    Fp <- p

    # --- PROTECTION LAYER 2: Optimized While Loop ---
    while(!is.na(Fp) && U > Fp) {
      y <- y + 1

      # Mathematical recurrence to avoid recalculating Z.
      p <- p * (lambda_i / (y^nu_i))

      if(is.na(p)) p <- 0

      Fp <- Fp + p

      # Protection against excessively long tails
      if(y > 10000) break
    }

    rands[i] <- y
  }

  return(rands)
}

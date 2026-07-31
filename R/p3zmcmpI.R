#' Cumulative Distribution Function of the Trivariate Zero-Modified COM-Poisson (Original Parametrization)
#'
#' @description Calculates the cumulative distribution function (CDF) of the
#' trivariate Zero-Modified COM-Poisson (MZMCMP) distribution based on the original lambda-nu parametrization.
#'
#' @param y1,y2,y3 Numeric vectors of quantiles.
#' @param lambda1,lambda2,lambda3 Numeric vectors. The rate/intensity parameters for each margin.
#' @param nu1,nu2,nu3 Numeric vectors. The dispersion parameters for each margin.
#' @param p Numeric vector. The zero-modification parameter.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p Logical; if TRUE, probabilities are given as log(p). Defaults to FALSE.
#'
#' @return A numeric vector of cumulative probabilities.
#' @author Rogerio A. Santana
#'
#' @examples
#' p3zmcmpI(y1 = 1, y2 = 2, y3 = 1,
#'          lambda1 = 3, lambda2 = 2, lambda3 = 4,
#'          nu1 = 0.8, nu2 = 1.2, nu3 = 0.9, p = 0.7)
#'
#' @export
p3zmcmpI <- function(y1, y2, y3, lambda1, lambda2, lambda3, nu1, nu2, nu3, p, lower.tail = TRUE, log.p = FALSE) {

  # --- CRAN Standard: Vector Recycling ---
  # Ensures all inputs have the same length so the user can pass vectors of quantiles
  n_max <- max(length(y1), length(y2), length(y3),
               length(lambda1), length(lambda2), length(lambda3),
               length(nu1), length(nu2), length(nu3), length(p))

  y1 <- rep_len(y1, n_max); y2 <- rep_len(y2, n_max); y3 <- rep_len(y3, n_max)
  lambda1 <- rep_len(lambda1, n_max); lambda2 <- rep_len(lambda2, n_max); lambda3 <- rep_len(lambda3, n_max)
  nu1 <- rep_len(nu1, n_max); nu2 <- rep_len(nu2, n_max); nu3 <- rep_len(nu3, n_max)
  p <- rep_len(p, n_max)

  res <- numeric(n_max)

  # Process each observation
  for (idx in seq_len(n_max)) {

    # Base case: CDF is 0 for negative values
    if (y1[idx] < 0 || y2[idx] < 0 || y3[idx] < 0) {
      res[idx] <- 0
      next
    }

    # Ensure quantiles are evaluated as integers (floor)
    Y1 <- floor(y1[idx])
    Y2 <- floor(y2[idx])
    Y3 <- floor(y3[idx])

    # Creates a grid of all integer points to sum over (Highly optimized approach)
    i <- rep(0:Y1, times = (Y2 + 1) * (Y3 + 1))
    j <- rep(0:Y2, each = (Y1 + 1), times = (Y3 + 1))
    k <- rep(0:Y3, each = (Y1 + 1) * (Y2 + 1))

    grid_size <- length(i)

    # Calls the vectorized density function once for the entire grid
    densities <- d3zmcmpI(
      y1 = i, y2 = j, y3 = k,
      lambda1 = rep(lambda1[idx], grid_size), lambda2 = rep(lambda2[idx], grid_size), lambda3 = rep(lambda3[idx], grid_size),
      nu1 = rep(nu1[idx], grid_size), nu2 = rep(nu2[idx], grid_size), nu3 = rep(nu3[idx], grid_size),
      p = rep(p[idx], grid_size)
    )

    # Sums all calculated densities in the grid
    res[idx] <- sum(densities)
  }

  # --- CRAN Standard: Tail and Log Adjustments ---
  if (!lower.tail) {
    res <- 1 - res
  }

  if (log.p) {
    res <- log(res)
  }

  return(res)
}

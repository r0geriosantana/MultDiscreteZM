#' Cumulative Distribution Function of the Trivariate Zero-Modified COM-Poisson
#'
#' @description Calculates the cumulative distribution function (CDF) of the
#' trivariate Zero-Modified COM-Poisson (MZMCMP) distribution parameterized by the mean.
#'
#' @param y1,y2,y3 Numeric vectors of quantiles.
#' @param mu1,mu2,mu3 Numeric vectors. The mean parameters for each margin.
#' @param kappa1,kappa2,kappa3 Numeric vectors. The dispersion parameters on the log scale.
#' @param p Numeric vector. The zero-modification parameter.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p Logical; if TRUE, probabilities are given as log(p). Defaults to FALSE.
#'
#' @return A numeric vector of cumulative probabilities.
#' @author Rogerio A. Santana
#'
#' @examples
#' p3zmcmpuI(y1 = 1, y2 = 2, y3 = 1,
#'           mu1 = 2, mu2 = 3, mu3 = 1.5,
#'           kappa1 = -0.5, kappa2 = 0.2, kappa3 = -0.1, p = 0.8)
#'
#' @export
p3zmcmpuI <- function(y1, y2, y3, mu1, mu2, mu3, kappa1, kappa2, kappa3, p, lower.tail = TRUE, log.p = FALSE) {

  # --- CRAN Standard: Vector Recycling ---
  # Ensures all inputs have the same length so the user can pass vectors of quantiles
  n_max <- max(length(y1), length(y2), length(y3),
               length(mu1), length(mu2), length(mu3),
               length(kappa1), length(kappa2), length(kappa3), length(p))

  y1 <- rep_len(y1, n_max); y2 <- rep_len(y2, n_max); y3 <- rep_len(y3, n_max)
  mu1 <- rep_len(mu1, n_max); mu2 <- rep_len(mu2, n_max); mu3 <- rep_len(mu3, n_max)
  kappa1 <- rep_len(kappa1, n_max); kappa2 <- rep_len(kappa2, n_max); kappa3 <- rep_len(kappa3, n_max)
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
    densities <- d3zmcmpuI(
      y1 = i, y2 = j, y3 = k,
      mu1 = rep(mu1[idx], grid_size), mu2 = rep(mu2[idx], grid_size), mu3 = rep(mu3[idx], grid_size),
      kappa1 = rep(kappa1[idx], grid_size), kappa2 = rep(kappa2[idx], grid_size), kappa3 = rep(kappa3[idx], grid_size),
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

#' Random Number Generation for the Trivariate ZM-Poisson Distribution
#'
#' @description Generates pseudo-random counts from the trivariate
#' Zero-Modified Poisson (Type I) distribution using
#' high-speed vectorized rejection sampling for the truncated component.
#'
#' @param n Integer. Number of observations to generate.
#' @param mu1,mu2,mu3 Numeric scalars. Mean parameters for each component.
#' @param p Numeric scalar. The zero-modification parameter (Probability of Y > 0).
#'
#' @return A data frame with n rows and 3 columns (Y1, Y2, Y3).
#' @export
r3zmpoisI <- function(n, mu1, mu2, mu3, p) {

  # --- 1. Parameter Validation ---
  if (p < 0 || p > 1) stop("Parameter 'p' must be a probability between 0 and 1.")

  # Initialize the output matrix with zeros
  Y <- matrix(0, nrow = n, ncol = 3)

  # --- 2. Hurdle Crossing Step ---
  # In a Zero-Modified (Hurdle) model, the probability of entering
  # the truncated count component is exactly 'p'.
  crossed_hurdle <- stats::rbinom(n, size = 1, prob = p)
  n_nonzero <- sum(crossed_hurdle == 1)

  # --- 3. Vectorized Rejection Sampling for Truncated Part ---
  if (n_nonzero > 0) {
    # Matrix to store valid non-zero trivariate counts
    mat_nonzero <- matrix(0, nrow = n_nonzero, ncol = 3)

    # Track which rows still need valid (sum > 0) candidates
    pending_indices <- 1:n_nonzero

    while (length(pending_indices) > 0) {
      n_pending <- length(pending_indices)

      # Generate Poisson candidates in blocks for all pending rows
      # NOTE: stats::rpois uses 'lambda', not 'mu'
      y1_cand <- stats::rpois(n_pending, lambda = mu1)
      y2_cand <- stats::rpois(n_pending, lambda = mu2)
      y3_cand <- stats::rpois(n_pending, lambda = mu3)

      # Condition for the truncated baseline distribution: sum must be strictly positive
      is_valid <- (y1_cand + y2_cand + y3_cand) > 0

      if (any(is_valid)) {
        # Extract rows that successfully passed the truncation hurdle
        valid_rows <- pending_indices[is_valid]

        # Store the successful candidates
        mat_nonzero[valid_rows, ] <- cbind(y1_cand[is_valid],
                                           y2_cand[is_valid],
                                           y3_cand[is_valid])

        # Remove successful rows from the pending tracking vector
        pending_indices <- pending_indices[!is_valid]
      }
    }

    # Inject the truncated counts back into the final matrix
    Y[crossed_hurdle == 1, ] <- mat_nonzero
  }

  # --- 4. Format Output ---
  return(data.frame(Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3]))
}

#' Random Number Generation for Trivariate ZMNB Regression (Type I)
#'
#' @description Generates pseudo-random counts from the trivariate
#' Zero-Modified Negative Binomial (Type I) distribution under a regression
#' framework, where the mean parameters depend on covariates.
#'
#' @param x Numeric vector or matrix of covariates.
#' @param beta1,beta2,beta3 Numeric vectors of regression coefficients for the means.
#' @param theta1,theta2,theta3 Numeric scalars. Dispersion parameters (size).
#' @param w0 Numeric scalar. The constant zero-modification parameter (Probability of Y > 0).
#'
#' @return A data frame with n rows and 3 columns (Y1, Y2, Y3).
#' @export
rmzmpoisI_reg0 <- function(x, beta1, beta2, beta3, theta1, theta2, theta3, w0) {

  # --- 1. Design Matrix & Dimensions ---
  # Safely handle both single vectors and pre-formatted model matrices
  n <- if (is.matrix(x)) nrow(x) else length(x)

  if (is.matrix(x) && colnames(x)[1] == "(Intercept)") {
    X <- x
  } else {
    X <- stats::model.matrix(~ x)
  }

  # --- 2. Calculate varying means (Forced to numeric vectors) ---
  mu1 <- as.numeric(exp(X %*% beta1))
  mu2 <- as.numeric(exp(X %*% beta2))
  mu3 <- as.numeric(exp(X %*% beta3))

  # --- 3. Hurdle Crossing Step ---
  if (w0 < 0 || w0 > 1) stop("Parameter 'w0' must be a probability between 0 and 1.")

  Y <- matrix(0, nrow = n, ncol = 3)
  hurdle <- stats::rbinom(n, size = 1, prob = w0)

  # Track indices of observations that entered the truncated component
  pending_indices <- which(hurdle == 1)

  # --- 4. Vectorized Block Rejection Sampling ---
  while (length(pending_indices) > 0) {
    n_pending <- length(pending_indices)

    # Generate candidates in bulk using the specific 'lambda=mumu' for each pending row
    cand1 <- stats::rpois(n_pending, lambda = mu1[pending_indices])
    cand2 <- stats::rpois(n_pending, lambda = mu2[pending_indices])
    cand3 <- stats::rpois(n_pending, lambda = mu3[pending_indices])

    # Validation: the joint observation cannot be (0,0,0)
    is_valid <- (cand1 + cand2 + cand3) > 0

    if (any(is_valid)) {
      # Extract rows that successfully passed the truncation hurdle
      valid_rows <- pending_indices[is_valid]

      # Insert accepted candidates into the final matrix
      Y[valid_rows, 1] <- cand1[is_valid]
      Y[valid_rows, 2] <- cand2[is_valid]
      Y[valid_rows, 3] <- cand3[is_valid]

      # Remove successful rows from the pending tracking vector
      pending_indices <- pending_indices[!is_valid]
    }
  }

  return(data.frame(Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3]))
}

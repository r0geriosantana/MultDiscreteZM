#' Random Generation for the Trivariate MZMCMPu Regression Model
#'
#' @description Generates random vectors from the trivariate Zero-Modified
#' COM-Poisson (MZMCMP) regression model (mean-parameterized), incorporating a single covariate.
#'
#' @param x Numeric vector. The covariate for the regression model.
#' @param beta1,beta2,beta3 Numeric vectors. The regression coefficients (intercept and slope) for each margin.
#' @param kappa1,kappa2,kappa3 Numeric scalars or vectors. The dispersion parameters on the log scale.
#' @param w0 Numeric scalar or vector. The hurdle probability (omega), representing P(Y > 0).
#'
#' @return A data frame with \code{n} rows and 3 columns (\code{Y1}, \code{Y2}, \code{Y3}).
#' @author Rogerio A. Santana
#'
#' @examples
#' set.seed(123)
#' x_cov <- rnorm(100)
#' sim_reg <- r3zmcmpuI_reg0(x = x_cov,
#'                           beta1 = c(0.5, 0.2), beta2 = c(0.1, -0.1), beta3 = c(0.8, 0.3),
#'                           kappa1 = -0.5, kappa2 = 0.2, kappa3 = -0.1,
#'                           w0 = 0.75)
#' head(sim_reg)
#'
#' @export
r3zmcmpuI_reg0 <- function(x, beta1, beta2, beta3, kappa1, kappa2, kappa3, w0) {

  n <- length(x)

  # Build the model matrix (Intercept + Covariate)
  X <- stats::model.matrix(~ x)

  # --- Regression Link (Log-link) ---
  # Convert to numeric vectors to ensure proper subsetting later
  mu1 <- as.numeric(exp(X %*% beta1))
  mu2 <- as.numeric(exp(X %*% beta2))
  mu3 <- as.numeric(exp(X %*% beta3))

  # Parameter recycling for dispersion and hurdle
  kappa1 <- rep_len(kappa1, n)
  kappa2 <- rep_len(kappa2, n)
  kappa3 <- rep_len(kappa3, n)
  omega <- rep_len(w0, n)

  # Initialize output matrix
  Y <- matrix(0, nrow = n, ncol = 3)

  # Determine which observations cross the hurdle
  hurdle <- stats::rbinom(n, size = 1, prob = omega)
  simula_zm <- which(hurdle == 1)
  n_sim <- length(simula_zm)

  if (n_sim > 0) {
    # Extract parameters only for the non-zero structural components
    m1_s <- mu1[simula_zm]; k1_s <- kappa1[simula_zm]
    m2_s <- mu2[simula_zm]; k2_s <- kappa2[simula_zm]
    m3_s <- mu3[simula_zm]; k3_s <- kappa3[simula_zm]

    cand1 <- numeric(n_sim)
    cand2 <- numeric(n_sim)
    cand3 <- numeric(n_sim)

    # --- Vectorized Rejection Sampling ---
    needs_update <- rep(TRUE, n_sim)

    while (any(needs_update)) {
      idx <- which(needs_update)
      n_update <- length(idx)

      # Generate candidates simultaneously
      cand1[idx] <- rcmp_u(n_update, m1_s[idx], k1_s[idx])
      cand2[idx] <- rcmp_u(n_update, m2_s[idx], k2_s[idx])
      cand3[idx] <- rcmp_u(n_update, m3_s[idx], k3_s[idx])

      # Check truncation condition (sum > 0)
      valid <- (cand1[idx] + cand2[idx] + cand3[idx]) > 0

      # Update mask
      needs_update[idx[valid]] <- FALSE
    }

    # Assign validated candidates
    Y[simula_zm, 1] <- cand1
    Y[simula_zm, 2] <- cand2
    Y[simula_zm, 3] <- cand3
  }

  return(data.frame(Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3]))
}

#' Random Generation for the Trivariate MZMCMP Regression Model (Original Parametrization)
#'
#' @description Generates random vectors from the trivariate Zero-Modified
#' COM-Poisson (MZMCMP) regression model based on the original lambda-nu parametrization,
#' incorporating a single covariate.
#'
#' @param x Numeric vector. The covariate for the regression model.
#' @param beta1,beta2,beta3 Numeric vectors. The regression coefficients (intercept and slope) for each margin.
#' @param nu1,nu2,nu3 Numeric scalars or vectors. The dispersion parameters.
#' @param w0 Numeric scalar or vector. The hurdle probability (omega), representing P(Y > 0).
#'
#' @return A data frame with \code{n} rows and 3 columns (\code{Y1}, \code{Y2}, \code{Y3}).
#' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
#'
#' @examples
#' set.seed(123)
#' x_cov <- stats::rnorm(100)
#' sim_reg <- r3zmcmpI_reg0(x = x_cov,
#'                          beta1 = c(0.5, 0.2), beta2 = c(0.1, -0.1), beta3 = c(0.8, 0.3),
#'                          nu1 = 0.8, nu2 = 1.2, nu3 = 0.9,
#'                          w0 = 0.75)
#' head(sim_reg)
#'
#' @export
r3zmcmpI_reg0 <- function(x, beta1, beta2, beta3, nu1, nu2, nu3, w0) {

  n <- length(x)

  # Build the model matrix (Intercept + Covariate)
  X <- stats::model.matrix(~ x)

  # --- Regression Link (Log-link) ---
  # Convert to numeric vectors to ensure proper subsetting later
  lambda1 <- as.numeric(exp(X %*% beta1))
  lambda2 <- as.numeric(exp(X %*% beta2))
  lambda3 <- as.numeric(exp(X %*% beta3))

  # Parameter recycling for dispersion and hurdle
  nu1 <- rep_len(nu1, n)
  nu2 <- rep_len(nu2, n)
  nu3 <- rep_len(nu3, n)
  omega <- rep_len(w0, n)

  # Initialize output matrix
  Y <- matrix(0, nrow = n, ncol = 3)

  # Determine which observations cross the hurdle
  hurdle <- stats::rbinom(n, size = 1, prob = omega)
  simula_zm <- which(hurdle == 1)
  n_sim <- length(simula_zm)

  if (n_sim > 0) {
    # Extract parameters only for the non-zero structural components
    l1_s <- lambda1[simula_zm]; n1_s <- nu1[simula_zm]
    l2_s <- lambda2[simula_zm]; n2_s <- nu2[simula_zm]
    l3_s <- lambda3[simula_zm]; n3_s <- nu3[simula_zm]

    cand1 <- numeric(n_sim)
    cand2 <- numeric(n_sim)
    cand3 <- numeric(n_sim)

    # --- Vectorized Rejection Sampling ---
    needs_update <- rep(TRUE, n_sim)

    while (any(needs_update)) {
      idx <- which(needs_update)
      n_update <- length(idx)

      # Generate candidates simultaneously using the correct base function
      cand1[idx] <- rcmp_sp(n_update, l1_s[idx], n1_s[idx])
      cand2[idx] <- rcmp_sp(n_update, l2_s[idx], n2_s[idx])
      cand3[idx] <- rcmp_sp(n_update, l3_s[idx], n3_s[idx])

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

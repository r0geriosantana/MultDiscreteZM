#ifndef UTILS_CMP_H
#define UTILS_CMP_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include <cmath>

using namespace Rcpp;

// ===================================================================
// FUNÇÕES AUXILIARES COMPARTILHADAS (Mestre)
// ===================================================================

inline double compute_log_S_scalar(double lambda, double nu, int max_iter = 100) {
  if (std::isnan(lambda)) return R_NaN;
  if (lambda == 0.0) return 0.0;

  NumericVector log_terms(max_iter + 1);
  log_terms[0] = 0.0;
  double log_lambda = std::log(lambda);
  for (int j = 1; j <= max_iter; ++j) {
    log_terms[j] = j * log_lambda - nu * std::lgamma(j + 1.0);
  }
  double max_log = Rcpp::max(log_terms);
  return max_log + std::log(Rcpp::sum(Rcpp::exp(log_terms - max_log)));
}

inline double dcmpu(double y, double mu, double kappa, int max_iter = 100) {
  double nu = std::exp(kappa);
  double base = mu + (nu - 1.0) / (2.0 * nu);
  if (base < 0.0) base = 0.0;
  double lambda = std::pow(base, nu);

  double log_S = compute_log_S_scalar(lambda, nu, max_iter);
  double log_prob = y * std::log(lambda) - nu * std::lgamma(y + 1.0) - log_S;

  if (!std::isfinite(log_prob)) return 0.0;
  return std::exp(log_prob);
}

#endif // UTILS_CMP_H

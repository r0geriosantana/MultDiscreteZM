#ifndef UTILS_CMP_H
#define UTILS_CMP_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include <cmath>
#include "cmp_core.h"   // truncamento adaptativo e inversao exata da media

using namespace Rcpp;

// ===================================================================
// FUNCOES AUXILIARES COMPARTILHADAS
// ===================================================================

// log S(lambda, nu) com truncamento adaptativo (Shmueli et al., 2005, Ap. B.1).
// O argumento max_iter foi mantido apenas por compatibilidade e e ignorado.
inline double compute_log_S_scalar(double lambda, double nu, int max_iter = 100) {
  if (std::isnan(lambda)) return R_NaN;
  if (lambda == 0.0) return 0.0;
  return cmp_series(lambda, nu).logS;
}

// pmf na parametrizacao pela media EXATA: lambda = g^{-1}(mu; nu)
inline double dcmpu(double y, double mu, double kappa, int max_iter = 100) {
  double lp = cmp_logpmf_u(y, mu, kappa);
  if (!std::isfinite(lp)) return 0.0;
  return std::exp(lp);
}

#endif // UTILS_CMP_H

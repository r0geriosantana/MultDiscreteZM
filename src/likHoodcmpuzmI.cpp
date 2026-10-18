#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include "utils_cmp.h" // Puxa as funcoes do arquivo mestre

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double likHoodcmpuzmI(NumericVector par, NumericVector y1,
                      NumericVector y2, NumericVector y3){
  int n = y1.size();

  double omega = 1.0 / (1.0 + exp(-par[0]));
  double mu1 = exp(par[1]);
  double mu2 = exp(par[2]);
  double mu3 = exp(par[3]);
  double kappa1 = par[4];
  double kappa2 = par[5];
  double kappa3 = par[6];

  double nu1 = exp(kappa1), nu2 = exp(kappa2), nu3 = exp(kappa3);
  // sem covariaveis: lambda e log S calculados uma unica vez
  double l1 = cmp_lambda_from_mu(mu1, nu1);
  double l2 = cmp_lambda_from_mu(mu2, nu2);
  double l3 = cmp_lambda_from_mu(mu3, nu3);
  double logS_tot = cmp_series(l1, nu1).logS + cmp_series(l2, nu2).logS + cmp_series(l3, nu3).logS;
  double logdenom = logS_tot + std::log1p(-std::exp(-logS_tot)); // log(prod S - 1)

  double loglik = 0.0;

  for (int i = 0; i < n; i++) {
    if (y1[i] == 0 && y2[i] == 0 && y3[i] == 0) {
      loglik += std::log1p(-omega);
    } else {
      loglik += std::log(omega);
      loglik += y1[i] * std::log(l1) - nu1 * std::lgamma(y1[i] + 1.0)
              + y2[i] * std::log(l2) - nu2 * std::lgamma(y2[i] + 1.0)
              + y3[i] * std::log(l3) - nu3 * std::lgamma(y3[i] + 1.0)
              - logdenom;
    }
  }

  if (!std::isfinite(loglik)) {
    return R_PosInf;
  }

  return -loglik;
}

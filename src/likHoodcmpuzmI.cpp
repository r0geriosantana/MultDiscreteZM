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

  double loglik = 0.0;

  for (int i = 0; i < n; i++) {

    if (y1[i] == 0 && y2[i] == 0 && y3[i] == 0) {

      double p_joint_zeros = (1.0 - omega);

      loglik += log(p_joint_zeros + 1e-15);
    } else {
      loglik += log(omega + 1e-15);

      double logcontagem = log(dcmpu(y1[i], mu1, kappa1)) + log(dcmpu(y2[i], mu2, kappa2)) + log(dcmpu(y3[i], mu3, kappa3)) + 1e-15;

      double den0 = dcmpu(0, mu1, kappa1) * dcmpu(0, mu2, kappa2) * dcmpu(0, mu3, kappa3);

      double logdenom = log((1.0 - den0) + 1e-15);

      loglik += logcontagem - logdenom;
    }
  }

  if (!std::isfinite(loglik)) {
    return R_PosInf;
  }

  return -loglik;
}

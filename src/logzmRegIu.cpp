#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include "utils_cmp.h" // Puxa as funcoes do arquivo mestre

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double logzmRegIu(NumericVector par,
                  NumericMatrix xmat,
                  NumericVector y1,
                  NumericVector y2,
                  NumericVector y3) {

  int n = y1.size();
  int ka = xmat.ncol();

  mat X(xmat.begin(), n, ka, false);
  arma::vec par_arma = Rcpp::as<arma::vec>(par);

  double gama = par_arma[0];

  int start_beta = 1;
  vec beta1 = par_arma.subvec(start_beta, start_beta + ka - 1);
  vec beta2 = par_arma.subvec(start_beta + ka, start_beta + 2 * ka - 1);
  vec beta3 = par_arma.subvec(start_beta + 2 * ka, start_beta + 3 * ka - 1);

  int start_kappa = start_beta + 3 * ka;
  double kappa1 = par_arma[start_kappa];
  double kappa2 = par_arma[start_kappa + 1];
  double kappa3 = par_arma[start_kappa + 2];

  vec mu1 = arma::exp(X * beta1);
  vec mu2 = arma::exp(X * beta2);
  vec mu3 = arma::exp(X * beta3);

  double omega = 1.0 / (1.0 + std::exp(-gama));
  double loglik = 0.0;

  for (int i = 0; i < n; i++) {
    int yy1 = y1[i]; int yy2 = y2[i]; int yy3 = y3[i];

    if (yy1 == 0 && yy2 == 0 && yy3 == 0) {
      double p_joint_zeros = (1.0 - omega);
      loglik += std::log(p_joint_zeros + 1e-15);
    } else {
      loglik += std::log(omega + 1e-15);

      double contagem = dcmpu(yy1, mu1[i], kappa1) * dcmpu(yy2, mu2[i], kappa2) * dcmpu(yy3, mu3[i], kappa3);
      double den0 = dcmpu(0, mu1[i], kappa1) * dcmpu(0, mu2[i], kappa2) * dcmpu(0, mu3[i], kappa3);

      loglik += std::log(contagem + 1e-15) - std::log((1.0 - den0) + 1e-15);
    }
  }

  if (!std::isfinite(loglik)) return R_PosInf;
  return -loglik;
}

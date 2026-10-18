#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include "utils_cmp.h" // Puxa as funcoes do arquivo mestre
#include "cmp_grad.h"  // verossimilhanca + gradiente analitico

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
  double nu1 = std::exp(kappa1), nu2 = std::exp(kappa2), nu3 = std::exp(kappa3);
  double loglik = 0.0;

  for (int i = 0; i < n; i++) {
    int yy1 = y1[i]; int yy2 = y2[i]; int yy3 = y3[i];

    if (yy1 == 0 && yy2 == 0 && yy3 == 0) {
      loglik += std::log1p(-omega);
    } else {
      loglik += std::log(omega);

      // lambda = g^{-1}(mu; nu) e log S calculados uma unica vez por componente
      double l1 = cmp_lambda_from_mu(mu1[i], nu1);
      double l2 = cmp_lambda_from_mu(mu2[i], nu2);
      double l3 = cmp_lambda_from_mu(mu3[i], nu3);
      double lS1 = cmp_series(l1, nu1).logS;
      double lS2 = cmp_series(l2, nu2).logS;
      double lS3 = cmp_series(l3, nu3).logS;

      double logcont = yy1 * std::log(l1) - nu1 * std::lgamma(yy1 + 1.0)
                     + yy2 * std::log(l2) - nu2 * std::lgamma(yy2 + 1.0)
                     + yy3 * std::log(l3) - nu3 * std::lgamma(yy3 + 1.0);
      double logS_tot = lS1 + lS2 + lS3;          // log prod S
      // log f+ = logcont - log(prod S - 1)
      loglik += logcont - (logS_tot + std::log1p(-std::exp(-logS_tot)));
    }
  }

  if (!std::isfinite(loglik)) return R_PosInf;
  return -loglik;
}

// Negative log-likelihood and its analytic gradient in one pass.
// Returns c(value, gradient). Used by optim(fn, gr) in the fitting functions.
// [[Rcpp::export]]
NumericVector logzmRegIu_vg(NumericVector par,
                            NumericMatrix xmat,
                            NumericVector y1,
                            NumericVector y2,
                            NumericVector y3) {
  int n = y1.size();
  int ka = xmat.ncol();
  int p = par.size();
  NumericVector out(p + 1);
  std::vector<double> g(p);
  out[0] = zmreg_negloglik_grad(par.begin(), xmat.begin(), n, ka,
                                y1.begin(), y2.begin(), y3.begin(), g.data());
  for (int k = 0; k < p; ++k) out[k + 1] = g[k];
  return out;
}

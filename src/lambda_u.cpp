#include <Rcpp.h>
#include "cmp_core.h"
using namespace Rcpp;

//' Rate parameter from the exact baseline mean
//'
//' @description Computes lambda = g^{-1}(mu; nu), the COM-Poisson rate whose
//' exact mean equals \code{mu}, with \code{nu = exp(kappa)}. Uses Newton's method
//' on log(lambda) with adaptive series truncation.
//'
//' @param mu Numeric vector of baseline means (mu > 0).
//' @param kappa Numeric vector of log-dispersion parameters (nu = exp(kappa)).
//' @return Numeric vector of rate parameters lambda.
//' @export
// [[Rcpp::export]]
NumericVector lambda_u(NumericVector mu, NumericVector kappa) {
  int n = std::max(mu.size(), kappa.size());
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double m = mu[i % mu.size()];
    double nu = std::exp(kappa[i % kappa.size()]);
    out[i] = cmp_lambda_from_mu(m, nu);
  }
  return out;
}

//' Exact moments of the COM-Poisson distribution
//'
//' @param lambda Numeric vector of rate parameters.
//' @param nu Numeric vector of dispersion parameters.
//' @return A matrix with columns logS, mean, var and K (truncation point).
//' @export
// [[Rcpp::export]]
NumericMatrix cmp_moments(NumericVector lambda, NumericVector nu) {
  int n = std::max(lambda.size(), nu.size());
  NumericMatrix out(n, 4);
  colnames(out) = CharacterVector::create("logS", "mean", "var", "K");
  for (int i = 0; i < n; i++) {
    CmpSeries s = cmp_series(lambda[i % lambda.size()], nu[i % nu.size()]);
    out(i, 0) = s.logS; out(i, 1) = s.mean; out(i, 2) = s.var; out(i, 3) = s.K;
  }
  return out;
}

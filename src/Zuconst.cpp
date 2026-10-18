#include <Rcpp.h>
#include "cmp_core.h"
using namespace Rcpp;

//' Normalization constant for Mean-Parameterized COM-Poisson (Rcpp)
//'
//' @name Zuconst
//' @description Normalization constant S(lambda, nu) with lambda = g^{-1}(mu; nu),
//' i.e. the rate whose exact mean equals mu. Adaptive truncation.
//'
//' @param mu A numeric vector. The exact baseline mean (mu > 0).
//' @param kappa A numeric vector. The dispersion parameter on the log scale (nu = exp(kappa)).
//' @param sum_limit Ignored; kept for backward compatibility.
//' @return A numeric vector containing the normalization constants.
//' @author Rogerio A. Santana,Katiane S. Conceicao, Carlos A. R. Diniz, Marinho G. Andrade
//' @export
// [[Rcpp::export]]
NumericVector Zuconst(NumericVector mu, NumericVector kappa, int sum_limit = -1) {
  int n = mu.size();
  NumericVector Z(n);
  for (int i = 0; i < n; i++) {
    double nu = std::exp(kappa[i % kappa.size()]);
    double lambda = cmp_lambda_from_mu(mu[i], nu);
    Z[i] = std::exp(cmp_series(lambda, nu).logS);
  }
  return Z;
}

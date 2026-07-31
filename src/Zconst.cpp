#include <Rcpp.h>
using namespace Rcpp;

//' Normalization constant for Original COM-Poisson (Rcpp)
//' 
//' @name Zconst
//' @description Calculates the normalization constant (Z) of the COM-Poisson distribution
//' using the original parameters lambda and nu, implemented in Rcpp for high performance.
//'
//' @param lambda A numeric vector. The rate/intensity parameter of the distribution (lambda > 0).
//' @param nu A numeric vector. The dispersion parameter (nu > 0).
//' @param sum_limit Integer. A safety buffer to ensure the summation series covers the distribution mass. Default is -1 (dynamic calculation).
//'
//' @return A numeric vector containing the normalization constants.
//' @author Rogerio A. Santana
//' @export
// [[Rcpp::export]]
NumericVector Zconst(NumericVector lambda, NumericVector nu, int sum_limit = -1) {
  int n = lambda.size();
  NumericVector Z(n);
  
  // 1. Heurística de Truncamento Dinâmico (Adaptada para lambda e nu)
  double max_lambda = max(lambda);
  double min_nu = min(nu);
  
  // Aproximação da média na parametrização original: lambda^(1/nu)
  double approx_mu = std::pow(std::max(max_lambda, 1e-10), 1.0 / min_nu);
  double proxy_sd = std::sqrt(approx_mu / min_nu);
  
  int iter_max = sum_limit;
  if (sum_limit <= 0) {
    iter_max = std::ceil(std::max(100.0, approx_mu + 20.0 * proxy_sd));
  }

  // 2. Loop de Cálculo
  for(int i = 0; i < n; i++) {
    // Proteção: lambda deve ser estritamente positivo para o log não retornar NaN
    double current_lambda = std::max(lambda[i], 1e-10);
    double current_nu = nu[i];

    double z_val = 0.0;
    double log_lambda = std::log(current_lambda);

    for(int j = 0; j <= iter_max; j++) {
      double log_term = j * log_lambda - current_nu * R::lgammafn(j + 1.0);
      z_val += std::exp(log_term);
    }

    Z[i] = z_val;
  }

  return Z;
}
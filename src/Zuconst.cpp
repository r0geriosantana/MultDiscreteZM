#include <Rcpp.h>
using namespace Rcpp;

//' Normalization constant for Mean-Parameterized COM-Poisson (Rcpp)
//' 
//' @name Zuconst
//' @description Calculates the normalization constant (Z) of the COM-Poisson distribution
//' using Rcpp for high performance.
 //'
 //' @param mu A numeric vector. The mean parameter of the distribution.
 //' @param kappa A numeric vector. The dispersion parameter on the log scale (nu = exp(kappa)).
 //' @param sum_limit Integer. A safety buffer to ensure the summation series covers the distribution mass. Default is -1 (dynamic calculation).
 //'
 //' @return A numeric vector containing the normalization constants.
 //' @author Rogerio A. Santana
 //' @export
 // [[Rcpp::export]]
 NumericVector Zuconst(NumericVector mu, NumericVector kappa, int sum_limit = -1) {
   int n = mu.size();
   NumericVector Z(n);
   NumericVector nu = exp(kappa);

   double max_mu = max(mu);
   double min_nu = min(nu);
   double proxy_sd = std::sqrt(max_mu / min_nu);

   int iter_max = sum_limit;
   if (sum_limit <= 0) {
     iter_max = std::ceil(std::max(100.0, max_mu + 20.0 * proxy_sd));
   }

   for(int i = 0; i < n; i++) {
     double current_mu = mu[i];
     double current_nu = nu[i];

     double base_term = current_mu + (current_nu - 1.0) / (2.0 * current_nu);
     double current_lambda = std::pow(std::max(base_term, 1e-10), current_nu);

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

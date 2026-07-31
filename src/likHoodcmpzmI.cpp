#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#include "utils_cmp.h" // Puxa as funcoes do arquivo mestre

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
 double likHoodcmpzmI(NumericVector par, NumericVector y1, NumericVector y2, NumericVector y3) {
   int n = y1.size();

   // Transform parameters back to their original constrained scales
   double omega = 1.0 / (1.0 + std::exp(-par[0]));
   double lambda1 = std::exp(par[1]);
   double lambda2 = std::exp(par[2]);
   double lambda3 = std::exp(par[3]);
   double nu1 = std::exp(par[4]);
   double nu2 = std::exp(par[5]);
   double nu3 = std::exp(par[6]);

   // --- CRITICAL OPTIMIZATION: CALCULATE CONSTANTS OUTSIDE THE LOOP ---
   // Since there are no covariates, parameters are identical for all observations.
   double log_S1 = compute_log_S_scalar(lambda1, nu1);
   double log_S2 = compute_log_S_scalar(lambda2, nu2);
   double log_S3 = compute_log_S_scalar(lambda3, nu3);

   // Log-probabilities of zero for each margin: log_p(0) = 0 - 0 - log_S = -log_S
   double p0_1 = std::exp(-log_S1);
   double p0_2 = std::exp(-log_S2);
   double p0_3 = std::exp(-log_S3);

   double den0 = p0_1 * p0_2 * p0_3;
   double log_denom = std::log(1.0 - den0 + 1e-15);

   double log_omega = std::log(omega + 1e-15);
   double log_1_minus_omega = std::log(1.0 - omega + 1e-15);

   double log_lambda1 = std::log(lambda1);
   double log_lambda2 = std::log(lambda2);
   double log_lambda3 = std::log(lambda3);

   double loglik = 0.0;

   // --- LIKELIHOOD EVALUATION LOOP ---
   for (int i = 0; i < n; i++) {

     if (y1[i] == 0.0 && y2[i] == 0.0 && y3[i] == 0.0) {
       // Hurdle baseline state
       loglik += log_1_minus_omega;
     } else {
       // Direct calculation of joint log-probability inside the loop using pre-calculated log_S
       double log_p1 = y1[i] * log_lambda1 - nu1 * R::lgammafn(y1[i] + 1.0) - log_S1;
       double log_p2 = y2[i] * log_lambda2 - nu2 * R::lgammafn(y2[i] + 1.0) - log_S2;
       double log_p3 = y3[i] * log_lambda3 - nu3 * R::lgammafn(y3[i] + 1.0) - log_S3;

       double log_contagem = log_p1 + log_p2 + log_p3;

       // log( omega * P_joint / (1 - P(0,0,0)) ) = log_omega + log_joint - log_denom
       loglik += log_omega + log_contagem - log_denom;
     }
   }

   // Penalty for optimization stability if log-likelihood blows up
   if (!std::isfinite(loglik)) {
     return 1e9;
   }

   return -loglik; // Return negative for minimization routines (e.g., optim)
 }

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>
#include <cmath>
#include "utils_cmp.h" // Puxa a compute_log_S_scalar do nosso arquivo mestre

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ===================================================================
// FUNCAO AUXILIAR: Densidade dcmp (Versao ESCALAR para Lambda)
// ===================================================================
// 'static inline' blinda a funcao para nao dar conflito com outras no pacote
static inline double dcmp(double y, double lambda, double nu, int max_iter = 100) {

  double log_S = compute_log_S_scalar(lambda, nu, max_iter);

  double log_prob = y * std::log(lambda) - nu * std::lgamma(y + 1.0) - log_S;

  if (!std::isfinite(log_prob)) {
    return 0.0; // Retorna 0 se a prob for -Inf ou NaN
  }

  return std::exp(log_prob);
}

// ===================================================================
// LIKELIHOOD PRINCIPAL (PARAMETRIZACAO EM LAMBDA)
// ===================================================================
// [[Rcpp::export]]
double logzmRegTriIcmp(NumericVector par,
                       NumericMatrix xmat,
                       NumericVector y1,
                       NumericVector y2,
                       NumericVector y3) {

  int n = y1.size();
  int ka = xmat.ncol();

  // Converte para objetos Armadillo
  mat X(xmat.begin(), n, ka, false);
  arma::vec par_arma = Rcpp::as<arma::vec>(par);

  // 1. Gamas (Interceptos Escalares)
  double gama = par_arma[0];

  // 2. Betas (Coeficientes da Regressao)
  int start_beta = 1;
  vec beta1 = par_arma.subvec(start_beta, start_beta + ka - 1);
  vec beta2 = par_arma.subvec(start_beta + ka, start_beta + 2 * ka - 1);
  vec beta3 = par_arma.subvec(start_beta + 2 * ka, start_beta + 3 * ka - 1);

  // 3. Kappas (Escalares)
  int start_kappa = start_beta + 3 * ka;
  double kappa1 = par_arma[start_kappa];
  double kappa2 = par_arma[start_kappa + 1];
  double kappa3 = par_arma[start_kappa + 2];

  // ==========================================================
  // CALCULO DOS PARAMETROS
  // ==========================================================

  // Calcula lambdas (Vetores) - Mantive os nomes mu1, mu2, mu3 do seu original,
  // mas aqui representam o preditor linear para o lambda
  vec mu1 = arma::exp(X * beta1);
  vec mu2 = arma::exp(X * beta2);
  vec mu3 = arma::exp(X * beta3);

  // Calcula o parametro de dispersao nu
  double nu1 = std::exp(kappa1);
  double nu2 = std::exp(kappa2);
  double nu3 = std::exp(kappa3);

  // Calcula probabilidades logisticas (w) (Escalares)
  double omega = 1.0 / (1.0 + std::exp(-gama));

  double loglik = 0.0;

  for (int i = 0; i < n; i++) {
    int yy1 = y1[i];
    int yy2 = y2[i];
    int yy3 = y3[i];

    double lambda1 = mu1[i];
    double lambda2 = mu2[i];
    double lambda3 = mu3[i];

    if (yy1 == 0 && yy2 == 0 && yy3 == 0) {
      double p_joint_zeros = (1.0 - omega);
      loglik += std::log(p_joint_zeros + 1e-15);
    } else {
      loglik += std::log(omega + 1e-15);

      double contagem = dcmp(yy1, lambda1, nu1) * dcmp(yy2, lambda2, nu2) * dcmp(yy3, lambda3, nu3);
      double den0 = dcmp(0, lambda1, nu1) * dcmp(0, lambda2, nu2) * dcmp(0, lambda3, nu3);

      double denom = (1.0 - den0);

      loglik += std::log(contagem + 1e-15) - std::log(denom + 1e-15);
    }
  }

  if (!std::isfinite(loglik)) {
    return R_PosInf; // Protecao para o otimizador
  }

  return -loglik;
}

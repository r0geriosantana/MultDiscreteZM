// cmp_grad.h -- negative log-likelihood of the trivariate MZMCMP_mu regression
// (hurdle form) and its analytic gradient. Pure C++ (no Rcpp), so that it can
// be tested on its own.
//
// Parameters: par = (gamma, beta_1 (ka), beta_2 (ka), beta_3 (ka), kappa_1..3)
//   omega = logit^{-1}(gamma) = P(Y != 0)
//   mu_ij = exp(x_j' beta_i), nu_i = exp(kappa_i), lambda_ij = g^{-1}(mu_ij; nu_i)
// Log-likelihood contribution:
//   Y_j = 0 : log(1 - omega)
//   Y_j != 0: log(omega) + sum_i [y_ij log(lambda_ij) - nu_i log(y_ij!)]
//             - log(prod_i S_ij - 1)
// Derivatives (mu fixed when differentiating w.r.t. nu):
//   d log(lambda)/d mu      = 1 / Var0
//   d log(lambda)/d nu      = Cov0(Y, log Y!) / Var0
//   d log S / d mu          = mu / Var0
//   d log S / d nu          = mu Cov0(Y, log Y!) / Var0 - E0(log Y!)
#ifndef CMP_GRAD_H
#define CMP_GRAD_H

#include <cmath>
#include <vector>
#include <limits>
#include "cmp_core.h"

struct CmpMoments {
  double lambda, logS, mean, var, ElogF, CovYlogF;   // logF = log Y!
};

// Series with the extra sums needed for the derivatives w.r.t. nu.
inline CmpMoments cmp_moments_ext(double lambda, double nu,
                                  double delta = 1e-12, int Kmax = 2000000) {
  CmpMoments out;
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  out.lambda = lambda;
  if (!(lambda > 0.0) || !std::isfinite(lambda) || !(nu > 0.0) || !std::isfinite(nu)) {
    out.logS = out.mean = out.var = out.ElogF = out.CovYlogF = NaN;
    return out;
  }
  const double ll = std::log(lambda), logdelta = std::log(delta);
  double M = 0.0, s0 = 1.0, s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
  double la = 0.0, lf = 0.0;           // log a_r and log r!
  int r = 0;
  for (;;) {
    double la_next = la + ll - nu * std::log(r + 1.0);
    double epsK    = std::exp(ll - nu * std::log(r + 1.0));
    if (epsK < 1.0) {
      // tail bound for r^2 a_r (dominates r log(r!) a_r while log(r!) <= r^2)
      double q = epsK * std::pow((r + 2.0) / (r + 1.0), 2);
      double logS_now = M + std::log(s0);
      bool ok0 = la_next - std::log1p(-epsK) - logS_now < logdelta;
      bool ok2 = (q < 1.0) &&
                 (la_next + 2.0 * std::log(r + 1.0) - std::log1p(-q)
                  - (M + std::log(s2 > 0 ? s2 : s0)) < logdelta);
      if ((ok0 && ok2) || r >= Kmax) break;
    } else if (r >= Kmax) break;
    ++r; la = la_next; lf += std::log((double) r);
    if (la > M) {
      double f = std::exp(M - la);
      s0 *= f; s1 *= f; s2 *= f; s3 *= f; s4 *= f; M = la;
    }
    double w = std::exp(la - M);
    s0 += w; s1 += r * w; s2 += (double) r * r * w; s3 += lf * w; s4 += r * lf * w;
  }
  out.logS     = M + std::log(s0);
  out.mean     = s1 / s0;
  out.var      = s2 / s0 - out.mean * out.mean;
  out.ElogF    = s3 / s0;
  out.CovYlogF = s4 / s0 - out.mean * out.ElogF;
  return out;
}

// Negative log-likelihood and gradient. X is n x ka in column-major order.
// Returns the value; fills grad (size 1 + 3 ka + 3). Non-finite -> +Inf.
inline double zmreg_negloglik_grad(const double* par, const double* X, int n, int ka,
                                   const double* y1, const double* y2, const double* y3,
                                   double* grad) {
  const int p = 1 + 3 * ka + 3;
  const double Inf = std::numeric_limits<double>::infinity();
  for (int k = 0; k < p; ++k) grad[k] = 0.0;
  const double gama  = par[0];
  const double omega = 1.0 / (1.0 + std::exp(-gama));
  double nu[3], kap[3];
  for (int i = 0; i < 3; ++i) { kap[i] = par[1 + 3 * ka + i]; nu[i] = std::exp(kap[i]); }
  const double* Y[3] = {y1, y2, y3};
  double ll = 0.0;
  for (int j = 0; j < n; ++j) {
    double yy[3] = {Y[0][j], Y[1][j], Y[2][j]};
    if (yy[0] == 0 && yy[1] == 0 && yy[2] == 0) {
      ll += std::log1p(-omega);
      grad[0] += -omega;                       // d/d gamma of log(1 - omega)
      continue;
    }
    ll += std::log(omega);
    grad[0] += 1.0 - omega;                    // d/d gamma of log(omega)
    CmpMoments m[3]; double mu[3];
    double logP = 0.0;
    for (int i = 0; i < 3; ++i) {
      double eta = 0.0;
      for (int c = 0; c < ka; ++c) eta += X[j + c * n] * par[1 + i * ka + c];
      mu[i] = std::exp(eta);
      double lam = cmp_lambda_from_mu(mu[i], nu[i]);
      m[i] = cmp_moments_ext(lam, nu[i]);
      if (!std::isfinite(m[i].logS) || !(m[i].var > 0.0)) return Inf;
      logP += m[i].logS;
      ll += yy[i] * std::log(lam) - nu[i] * std::lgamma(yy[i] + 1.0);
    }
    // - log(P - 1) = -(logP + log(1 - 1/P))
    ll += -(logP + std::log1p(-std::exp(-logP)));
    const double cfac = 1.0 / (-std::expm1(-logP));          // P / (P - 1)
    for (int i = 0; i < 3; ++i) {
      double dmu  = (yy[i] - cfac * mu[i]) / m[i].var;        // d ll / d mu
      double dnu  = yy[i] * m[i].CovYlogF / m[i].var - std::lgamma(yy[i] + 1.0)
                  - cfac * (mu[i] * m[i].CovYlogF / m[i].var - m[i].ElogF);
      for (int c = 0; c < ka; ++c) grad[1 + i * ka + c] += dmu * mu[i] * X[j + c * n];
      grad[1 + 3 * ka + i] += dnu * nu[i];                     // d / d kappa
    }
  }
  if (!std::isfinite(ll)) return Inf;
  for (int k = 0; k < p; ++k) {
    grad[k] = -grad[k];                         // negative log-likelihood
    if (!std::isfinite(grad[k])) grad[k] = 0.0;
  }
  return -ll;
}

#endif

// cmp_core.h -- nucleo numerico da COM-Poisson (C++ puro, sem Rcpp)
// Truncamento adaptativo (Shmueli et al., 2005, Ap. B.1) e
// inversao exata da media: lambda = g^{-1}(mu; nu).
#ifndef CMP_CORE_H
#define CMP_CORE_H

#include <cmath>
#include <limits>

struct CmpSeries {
  double logS;   // log S(lambda, nu)
  double mean;   // E0(Y)
  double var;    // Var0(Y)
  int    K;      // ultimo indice usado
};

// Soma log a_r, r*a_r, r^2*a_r com reescala online (log-sum-exp).
// Para quando eps_K = lambda/(K+1)^nu < 1 e o limite relativo da cauda
// (incluindo a cauda de r^2 a_r) e menor que delta.
inline CmpSeries cmp_series(double lambda, double nu,
                            double delta = 1e-12, int Kmax = 1000000) {
  CmpSeries out;
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  if (!(lambda > 0.0) || !std::isfinite(lambda) || !(nu > 0.0) || !std::isfinite(nu)) {
    out.logS = NaN; out.mean = NaN; out.var = NaN; out.K = 0;
    return out;
  }
  const double ll = std::log(lambda);
  double M  = 0.0;                    // maximo corrente de log a_r (a_0 = 1)
  double s0 = 1.0, s1 = 0.0, s2 = 0.0; // somas reescaladas por exp(-M)
  double la = 0.0;                    // log a_r corrente
  int r = 0;
  const double logdelta = std::log(delta);
  for (;;) {
    double la_next = la + ll - nu * std::log(r + 1.0);  // log a_{r+1}
    double epsK    = std::exp(ll - nu * std::log(r + 1.0));
    if (epsK < 1.0) {
      // cauda de a_r e de r^2 a_r (razao desta ultima <= epsK*((r+2)/(r+1))^2)
      double q   = epsK * std::pow((r + 2.0) / (r + 1.0), 2);
      double logS_now = M + std::log(s0);
      bool ok0 = la_next - std::log1p(-epsK) - logS_now < logdelta;
      bool ok2 = (q < 1.0) &&
                 (la_next + 2.0 * std::log(r + 1.0) - std::log1p(-q)
                   - (M + std::log(s2 > 0 ? s2 : s0)) < logdelta);
      if ((ok0 && ok2) || r >= Kmax) break;
    } else if (r >= Kmax) break;
    ++r; la = la_next;
    if (la > M) { double f = std::exp(M - la); s0 *= f; s1 *= f; s2 *= f; M = la; }
    double w = std::exp(la - M);
    s0 += w; s1 += r * w; s2 += (double)r * r * w;
  }
  out.logS = M + std::log(s0);
  out.mean = s1 / s0;
  out.var  = s2 / s0 - out.mean * out.mean;
  out.K = r;
  return out;
}

// lambda = g^{-1}(mu; nu). eta0: valor inicial de log(lambda) (NaN = automatico).
inline double cmp_lambda_from_mu(double mu, double nu,
                                 double eta0 = std::numeric_limits<double>::quiet_NaN(),
                                 double tol = 1e-10, int maxit = 50) {
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  // Entradas fora do dominio util: devolve NaN (a verossimilhanca vira penalidade)
  if (!std::isfinite(mu) || !(mu > 0.0) || mu > 1e5 ||
      !std::isfinite(nu) || !(nu > 0.0)) return NaN;
  // Limite superior de log(lambda) para que a media nao exceda ~1e5:
  // para nu fixo, g(lambda) >= lambda^(1/nu) - (1 - nu)/(2 nu) aproximadamente,
  // e usamos um teto generoso em eta para evitar series gigantes.
  const double eta_max = nu * std::log(1e5 + 10.0) + 5.0;
  double eta = eta0;
  if (std::isnan(eta)) {
    double b = mu + (nu - 1.0) / (2.0 * nu);
    eta = (b > 1.0) ? nu * std::log(b) : std::log(mu / (1.0 + mu));
  }
  if (eta > eta_max) eta = eta_max;
  const double crit = tol * (mu > 1.0 ? mu : 1.0);
  for (int k = 0; k < maxit; ++k) {
    CmpSeries s = cmp_series(std::exp(eta), nu, 1e-12, 2000000);
    if (!std::isfinite(s.mean) || !std::isfinite(s.var) || !(s.var > 0.0)) break;
    double d = s.mean - mu;
    if (std::fabs(d) < crit) return std::exp(eta);
    double step = d / s.var;
    if (step >  1.0) step =  1.0;
    if (step < -1.0) step = -1.0;
    eta -= step;
    if (eta > eta_max) eta = eta_max;
  }
  // fallback: bissecao em eta (g crescente), com numero maximo de passos
  double lo = std::log(mu / (1.0 + mu)) - 5.0, hi = eta_max;
  CmpSeries shi = cmp_series(std::exp(hi), nu, 1e-12, 2000000);
  if (!std::isfinite(shi.mean) || shi.mean < mu) return NaN;
  for (int k = 0; k < 40 && cmp_series(std::exp(lo), nu).mean > mu; ++k) lo -= 5.0;
  for (int k = 0; k < 200 && hi - lo > 1e-13; ++k) {
    double mid = 0.5 * (lo + hi);
    double m = cmp_series(std::exp(mid), nu).mean;
    if (!std::isfinite(m)) return NaN;
    if (m < mu) lo = mid; else hi = mid;
  }
  return std::exp(0.5 * (lo + hi));
}

// log f(y | mu, kappa) na parametrizacao pela media exata
inline double cmp_logpmf_u(double y, double mu, double kappa) {
  double nu = std::exp(kappa);
  double lambda = cmp_lambda_from_mu(mu, nu);
  CmpSeries s = cmp_series(lambda, nu);
  return y * std::log(lambda) - nu * std::lgamma(y + 1.0) - s.logS;
}

#endif

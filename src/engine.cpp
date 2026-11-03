#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export(name = ".LogLik_cpp")]]
List LogLik_cpp(NumericVector omega, NumericMatrix pi_mat, NumericMatrix R,
                NumericVector m, NumericMatrix xi_mat, NumericMatrix uniform_j) {

  int n = R.nrow();
  int J = R.ncol();
  int K = omega.size();
  const double NEG_INF = -std::numeric_limits<double>::infinity();

  NumericMatrix tau(n, K);
  NumericVector logf_i(n);          // log della verosimiglianza marginale del soggetto i
  double LL = 0.0;
  std::vector<double> la(K);

  for (int i = 0; i < n; i++) {
    double mx = NEG_INF;

    for (int k = 0; k < K; k++) {
      double lf = 0.0;
      for (int j = 0; j < J; j++) {
        double binom_val = R::dbinom(R(i, j) - 1.0, m[j] - 1.0, 1.0 - xi_mat(k, j), 0);
        double dens = pi_mat(k, j) * binom_val + (1.0 - pi_mat(k, j)) * uniform_j(i, j);
        lf += std::log(dens);       // log(0) = -Inf, gestito sotto
      }
      la[k] = (omega[k] > 0.0 ? std::log(omega[k]) : NEG_INF) + lf;
      if (la[k] > mx) mx = la[k];
    }

    if (!std::isfinite(mx)) {       // configurazione davvero impossibile
      for (int k = 0; k < K; k++) tau(i, k) = NA_REAL;
      logf_i[i] = NEG_INF;
      LL += NEG_INF;
      continue;
    }

    double s = 0.0;
    for (int k = 0; k < K; k++) s += std::exp(la[k] - mx);
    double lse = mx + std::log(s);

    for (int k = 0; k < K; k++) tau(i, k) = std::exp(la[k] - lse);  // righe somma a 1
    logf_i[i] = lse;
    LL += lse;
  }

  return List::create(Named("tau") = tau,
                      Named("logMLCCUBk") = logf_i,
                      Named("LL") = LL);
}

// Diciamo a Rcpp di esportare questa funzione in R come ".eta_ijk_cpp"
// [[Rcpp::export(name = ".eta_ijk_cpp")]]
NumericVector eta_ijk_cpp(NumericMatrix pi_mat, NumericMatrix R, NumericVector m, 
                          NumericMatrix xi_mat, NumericMatrix uniform_j, 
                          int K, int J, int n) {
  
  NumericVector eta(K * J * n);
  eta.attr("dim") = IntegerVector::create(K, J, n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < J; j++) {
      for(int k = 0; k < K; k++) {
        double binom_val = R::dbinom(R(i, j) - 1.0, m[j] - 1.0, 1.0 - xi_mat(k, j), 0);
        double num = pi_mat(k, j) * binom_val;
        double den = num + (1.0 - pi_mat(k, j)) * uniform_j(i, j);
        
        int idx = k + j * K + i * K * J; 
        
        // Avoid division by zero
        if (den > 0) {
          eta[idx] = num / den;
        } else {
          eta[idx] = 0.0;
        }
      }
    }
  }
  return eta;
}

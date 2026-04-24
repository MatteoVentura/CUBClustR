#include <Rcpp.h>
using namespace Rcpp;

// Diciamo a Rcpp di esportare questa funzione in R come ".LogLik_cpp"
// [[Rcpp::export(name = ".LogLik_cpp")]]
List LogLik_cpp(NumericVector omega, NumericMatrix pi_mat, NumericMatrix R, 
                NumericVector m, NumericMatrix xi_mat, NumericMatrix uniform_j) {
  
  int n = R.nrow(); 
  int J = R.ncol(); 
  int K = omega.size();
  
  NumericMatrix multiCUB_jk(n, K);
  NumericVector MLCCUBk(n);
  double LL = 0.0;
  
  for(int i = 0; i < n; i++) {
    double sum_k = 0.0;
    for(int k = 0; k < K; k++) {
      double prod_j = 1.0;
      for(int j = 0; j < J; j++) {
        // Compute Binomial density
        double binom_val = R::dbinom(R(i, j) - 1.0, m[j] - 1.0, 1.0 - xi_mat(k, j), 0);
        double dens = pi_mat(k, j) * binom_val + (1.0 - pi_mat(k, j)) * uniform_j(i, j);
        prod_j *= dens;
      }
      multiCUB_jk(i, k) = omega[k] * prod_j;
      sum_k += multiCUB_jk(i, k);
    }
    
    MLCCUBk[i] = sum_k;
    
    // Add to LogLikelihood, avoiding log(0) crashes
    if (sum_k > 0) {
      LL += std::log(sum_k);
    } else {
      LL += -1e10; // Heavy penalty for impossible configurations
    }
  }
  
  return List::create(Named("multiCUB_w_k") = multiCUB_jk, 
                      Named("MLCCUBk") = MLCCUBk, 
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

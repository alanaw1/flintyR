#include <RcppArmadillo.h>

using namespace Rcpp;

#include "Rcpp.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "algorithm"
#include "climits"



//' Binary Vector to Bit
//'
//' Converts the first column of an integer matrix X
//' into its base 2 or bit representation. Used to speed up
//' Hamming distance vector computation.
//'
//' @param x (IntegerMatrix)
//' @return bit representation of first column (unsigned long long)
unsigned long long binvec2int(IntegerMatrix x) {
  unsigned long long u = 0;
  int n = x.nrow();
  for (int i = 0; i < n; i++) {
    u += x(i, 0) * (unsigned long long) pow(2, n - 1 - i);
  }
  return u;
}


//' Fast Bitwise Hamming Distance Vector Computation
//'
//' Takes in a binary matrix X, whose transpose t(X)
//' has N rows, and computes a vector recording all
//' {N choose 2} pairwise Hamming distances of t(X),
//' ordered lexicographically.
//'
//' @param X binary matrix (IntegerMatrix class )
//' @return vector of Hamming distances (NumericVector class)
//' @examples # t(X) = [[1,0], [0,1], [1,1]] --> output = [2,1,1]
// [[Rcpp::export]]
NumericVector hamming_bitwise(IntegerMatrix X) {
  // initialize some variables
  int i, j, k, from, to, h;
  int stepsize = sizeof(unsigned long long) * CHAR_BIT;
  int n = ceil(((double) X.nrow()) / stepsize);
  IntegerMatrix x;

  // initialize matrix H, to be modified below
  NumericMatrix H(X.ncol(), X.ncol());

  // initialize output vector HamVec, to be modified from H and returned
  NumericVector HamVec(X.ncol() * (X.ncol() - 1) / 2);

  // convert X to unsigned long long
  unsigned long long *xx =
    (unsigned long long *) malloc(n * X.ncol() * sizeof(unsigned long long));
  i = 0;
  while ((from = i * stepsize) < X.nrow()) {
    for (j = 0; j < X.ncol(); j++) {
      to = std::min(X.nrow(), (i + 1) * stepsize);
      x = X(Range(from, to - 1), Range(j, j));
      xx[i + j * n] = binvec2int(x);
    }
    i++;
  }

  // fill in entries of H
  for (i = 0; i < X.ncol(); i++) {
    for (j = i + 1; j < X.ncol(); j++) {
      for (k = 0; k < n; k++) {
        h = __builtin_popcountll(xx[k + i * n] ^ xx[k + j * n]);
        H(i, j) += h;
        H(j, i) += h;
      }
    }
  }
  free(xx);
  // convert from N x N distance matrix, H, to {N choose 2} vector
  for (i = 0; i < X.ncol(); i++) {
    for (j = i + 1; j < X.ncol(); j++) {
      h = H(i,j);
      k = X.ncol() * i - (i + 2) * (i + 1) / 2 + j;
      HamVec(k) = h;
    }
  }
  return HamVec;
}

//' Fast \eqn{l_p^p} Distance Vector Computation
//'
//' Takes in a double matrix X, whose transpose t(X)
//' has N rows, and computes a vector recording all
//' \eqn{{N \choose 2}} pairwise \eqn{l_p^p} distances of t(X),
//' ordered lexicographically.
//'
//' @param X double matrix (arma::mat class)
//' @param p numeric Minkowski power (double class)
//' @return vector of \eqn{l_p^p} distances (arma::vec class)
//' @examples # X = [[0.5,0.5],[0,1],[0.3,0.7]] --> lPVec = [x,y,z]
//' # with x = (0.5^p + 0.5^p)
// [[Rcpp::export]]
arma::vec lp_distance(arma::mat X, double p) {
  // initialize some variables
  int i, j, k;
  double h;
  arma::vec xi, xj;

  // initialize matrix H, to be modified below
  arma::mat H(X.n_cols, X.n_cols);

  // initialize output vector lPVec, to be modified from H and returned
  arma::vec lPVec(X.n_cols * (X.n_cols - 1) / 2);

  // fill in entries of Hs
  for (i = 0; i < X.n_cols; i++) {
    xi = X.col(i);
    for (j = i + 1; j < X.n_cols; j++) {
      h = pow(arma::norm(xi - X.col(j), p), p);
      H(i,j) = h;
      H(j,i) = h;
    }
  }

  // convert from N x N distance matrix, H, to {N choose 2} vector
  for (i = 0; i < X.n_cols; i++) {
    for (j = i + 1; j < X.n_cols; j++) {
      h = H(i,j);
      k = X.n_cols * i - (i + 2) * (i + 1) / 2 + j;
      lPVec(k) = h;
    }
  }
  return(lPVec);
}



#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP Matprod(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP Matprodt(Eigen::MatrixXd A, Eigen::MatrixXd B){
      Eigen::MatrixXd C = A * B.transpose();
    return Rcpp::wrap(C);
}


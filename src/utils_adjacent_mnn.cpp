#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::sp_mat convertSparse(Rcpp::S4& mat) {

      // obtain dim, i, p. x from S4 object
      Rcpp::IntegerVector dims = mat.slot("Dim");
      arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
      arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
      arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

      int nrow = dims[0], ncol = dims[1];

      // use Armadillo sparse matrix constructor
      arma::sp_mat res(i, p, x, nrow, ncol);
      return res;
}

#include <RcppArmadillo.h>
// #include <RcppParallel.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]

//[[Rcpp::export]]
arma::sp_mat calAdjacent(arma::mat& position_,
                         arma::Col<int>& batch_vec_,
                         std::string& platform,
                         const double& cutoff=50,
                         const int& ncores = 4) {
    using arma::uword;
    using arma::uvec;
    omp_set_num_threads(ncores);

    // number of spots in all batches
    const uword n(position_.n_rows);
    // allocate space for adjacent matrix
    arma::sp_mat adj_mat(n,n);
    uvec inx, inx2;
    // parallel
#pragma omp parallel for schedule(static) private(inx, inx2)
    for(uword rw=0; rw<n-1; ++rw) {
        if(platform == "visium") {
            // find index satisfing conditions in same batch
            inx = arma::find((abs(position_(rw,0) - position_.col(0))<=2) &&
                (abs(position_(rw,1) - position_.col(1))<=1) && (batch_vec_==batch_vec_(rw)));
            inx2 = arma::find(inx > rw);
            for(uword cl=0; cl<inx2.n_elem; ++cl) {
                adj_mat(rw, inx(inx2(cl))) = 1;
                adj_mat(inx(inx2(cl)), rw) = 1;
            }
        }else if(platform == "st") {
            inx = arma::find((abs(position_(rw,0) - position_.col(0))<=1) &&
                (abs(position_(rw,1) - position_.col(1))<=1) && (batch_vec_==batch_vec_(rw)));
            inx2 = arma::find(inx > rw);
            for(uword cl=0; cl<inx2.n_elem; ++cl) {
                adj_mat(rw, inx(inx2(cl))) = 1;
                adj_mat(inx(inx2(cl)), rw) = 1;
            }
        }else{
            inx = arma::find((abs(position_(rw,0) - position_.col(0))<cutoff) &&
                (abs(position_(rw,1) - position_.col(1))<cutoff) && (batch_vec_==batch_vec_(rw)));
            // select rows > given rows, only fill upper triangular
            inx2 = arma::find(inx > rw);
            for(uword cl=0; cl<inx2.n_elem; ++cl) {
                adj_mat(rw, inx(inx2(cl))) = (norm(position_.row(inx(inx2(cl))) - position_.row(rw), 2) < cutoff);
                adj_mat(inx(inx2(cl)), rw) = adj_mat(rw, inx(inx2(cl)));
            }
        }
    }
    return adj_mat;
}



// //[[Rcpp::export]]
// arma::sp_mat calAdjacent(Rcpp::NumericMatrix& position,
//                          Rcpp::IntegerVector& batch_vec,
//                          std::string& platform,
//                          const double& cutoff=2.1,
//                          const int& ncores = 4) {
//     using arma::uword;
//     using arma::uvec;
//     omp_set_num_threads(ncores);
//     // convert position to arma object
//     arma::mat position_(position.begin(), position.nrow(), position.ncol(), false);
//     arma::Col<int> batch_vec_(batch_vec.begin(), batch_vec.size(), false);
//
//     const uword n(position_.n_rows);
//     arma::sp_mat adj_mat(n,n);
//     uvec inx, inx2;
//     // parallel
//     #pragma omp parallel for schedule(static) private(inx, inx2)
//     for(uword rw=0; rw<n-1; ++rw) {
//         if(platform == "visium") {
//           // find index satisfing conditions in same batch
//           inx = arma::find((abs(position_(rw,0) - position_.col(0))<=2) &&
//             (abs(position_(rw,1) - position_.col(1))<=1) && (batch_vec_==batch_vec_(rw)));
//           inx2 = arma::find(inx > rw);
//           for(uword cl=0; cl<inx2.n_elem; ++cl) {
//             adj_mat(rw, inx(inx2(cl))) = 1;
//             adj_mat(inx(inx2(cl)), rw) = 1;
//           }
//         }else if(platform == "st") {
//           inx = arma::find((abs(position_(rw,0) - position_.col(0))<=1) &&
//             (abs(position_(rw,1) - position_.col(1))<=1) && (batch_vec_==batch_vec_(rw)));
//           inx2 = arma::find(inx > rw);
//           for(uword cl=0; cl<inx2.n_elem; ++cl) {
//             adj_mat(rw, inx(inx2(cl))) = 1;
//             adj_mat(inx(inx2(cl)), rw) = 1;
//           }
//         }else{
//           inx = arma::find((abs(position_(rw,0) - position_.col(0))<cutoff) &&
//             (abs(position_(rw,1) - position_.col(1))<cutoff) && (batch_vec_==batch_vec_(rw)));
//           // select rows > given rows, only fill upper triangular
//           inx2 = arma::find(inx > rw);
//           for(uword cl=0; cl<inx2.n_elem; ++cl) {
//             adj_mat(rw, inx(inx2(cl))) = (norm(position_.row(inx(inx2(cl))) - position_.row(rw), 2) < cutoff);
//             adj_mat(inx(inx2(cl)), rw) = adj_mat(rw, inx(inx2(cl)));
//           }
//         }
//     }
//     return adj_mat;
// }
//

#ifndef utils_em_intern
#define utils_em_intern

#include <RcppArmadillo.h>
#include <omp.h>

//----------------------------------------------------
Rcpp::List calPhi(const arma::mat& T, const arma::mat& L,
                  const arma::mat& Lambda, const arma::vec& omega_vec,
                  const arma::Col<int>& batch_vec);
arma::mat calEnergy_latent(const arma::sp_mat& adj_mat, const arma::Col<int>& c_vec,
                           const int& K, const arma::vec& eta_vec,
                           const arma::Col<int>& batch_vec);
Rcpp::List calEnergy_observed(const arma::mat& x, const arma::mat& B, const arma::mat& L,
                              const arma::mat& Mu, const arma::mat& T, const arma::mat& Lambda,
                              const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                              const int& ncores=4);
arma::mat extractDiag(const arma::cube& Efft_li);
arma::mat extractJcol(const arma::cube& Efft_li, const int& j);
arma::mat contT2mat(const arma::mat& T, const arma::Col<int>& batch_vec);
arma::mat contB2mat(const arma::mat& B, const arma::Col<int>& batch_vec);
arma::mat solveQ(const arma::mat& A, const arma::mat& B, const arma::mat& C);
// arma::mat solveQ(const arma::mat& A, const arma::mat& B,
//                  const arma::mat& C, const int& ncores);
//----------------------------------------------------

#endif//utils_em_intern

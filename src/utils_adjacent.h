#ifndef utils_adjacent
#define utils_adjacent

#include <RcppArmadillo.h>
// #include <RcppParallel.h>
#include <omp.h>

//----------------------------------------------------
arma::sp_mat calAdjacent(arma::mat& position,
                         arma::Col<int>& batch_vec,
                         std::string& platform,
                         const double& cutoff=2.1,
                         const int& ncores = 4);
//----------------------------------------------------

#endif//utils_adjacent

#ifndef utils_em
#define utils_em

#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <omp.h>

//----------------------------------------------------
Rcpp::List update_c_sp(const arma::mat& x, const arma::sp_mat& adj_mat,
                       arma::Col<int>& c_vec, const arma::mat& B,
                       const arma::mat& L, const arma::mat& Mu,
                       const arma::mat& T, const arma::mat& Lambda,
                       const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                       arma::vec& eta_vec, const int& icm_maxiter = 10,
                       const bool& verbose = true, const int& ncores=4);
Rcpp::List update_c_sc(const arma::mat& x, const arma::mat& B,
                       const arma::mat& L, const arma::mat& Mu,
                       const arma::mat& T, const arma::mat& Lambda,
                       const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                       arma::mat& pk_mat, const int& ncores=4);
Rcpp::List expect_f(const arma::mat& x, const Rcpp::List& inv_Phis_list,
                    const arma::mat& B, const arma::mat& Lambda,
                    const arma::mat& Mu, const arma::mat& T, const arma::mat& L,
                    const arma::mat& S_lik ,const arma::vec& omega_vec,
                    const arma::Col<int>& batch_vec,
                    const int& ncores);
arma::mat expect_gamma(const arma::mat& L, const double& lambda0,
                       const double& lambda1, const arma::vec& p_vec,
                       const int& ncores);
arma::mat update_T(const arma::mat& x, const arma::mat& B,  const arma::mat& L,
                   const arma::cube& Efft_li, const arma::mat& Ef_li,
                   const arma::Col<int> batch_vec, const double& nu_tau,
                   const int& ncores);
arma::mat update_B(const arma::mat& x, const arma::mat& T,
                   const arma::mat& L, const arma::mat& Ef_li,
                   const arma::Col<int>& batch_vec, const arma::uword& n_b,
                   const int& ncores);
arma::mat update_mu(const arma::mat& S_lik, const arma::mat& Lambda,
                    const arma::vec& mu_mu, const arma::mat& sigma_mu,
                    const arma::mat& Ef_k, const arma::vec& omega_vec,
                    const int& ncores);
arma::vec update_omega(const arma::mat& S_lik, const arma::mat& Mu,
                       const arma::mat& Lambda, const arma::cube& Efft_li,
                       const arma::cube& weightEf_liks, const double& nu_omega,
                       const int& ncores);
arma::mat update_Lambda(const arma::cube& Efft_li, const arma::mat& muEf,
                        const arma::vec& omega_vec, const arma::mat& S_lik,
                        const arma::mat& Mu, const int& n_Lambda,
                        const arma::mat& Sigma_Lambda, const int& ncores);
arma::vec update_p(const arma::mat& Egamma, const double& alpha_p,
                   const double& beta_p, const int& ncores);
arma::mat update_L(const arma::mat& x, const arma::mat& B,
                   const arma::mat& Egamma, const arma::Col<int>& batch_vec,
                   const arma::mat& Ef_li, const arma::cube& Efft_li,
                   const arma::mat& T, const arma::mat& L,
                   const double& lambda0, const double& lambda1,
                   const int& ncores);
double logPost(const double& loglik, const arma::mat& L,
               const arma::mat& B, const arma::mat& mu,
               const arma::vec& p_vec, const arma::vec& omega_vec,
               const arma::mat& T, const arma::mat& Lambda,
               const double& lambda0, const double& lambda1,
               const arma::vec& mu_mu,const arma::mat& Sigma_mu,
               const double& nu_omega, const double& nu_tau,
               const int& n_Lambda, const arma::mat& Sigma_Lambda,
               const double& alpha_p, const double& beta_p);
//----------------------------------------------------

#endif//utils_em

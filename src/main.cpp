#include<RcppArmadillo.h>
#include<RcppDist.h>
#include <omp.h>
#include "utils_adjacent.h"
#include "utils_adjacent_mnn.h"
#include "utils_em_intern.h"
#include "utils_em.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::depends(RcppArmadillo,RcppDist)]]
//[[Rcpp::plugins(openmp)]]

 //' @title stadia_EM_SP
 //'
 //' @description main function to run stadia algorithm for spatial transcriptomics
 //'
 //' @param x a pxn matrix, preprocessed count data.
 //' @param position a nx2 matrix containing row & col.
 //' @param batch_vec vector of batch information.
 //' @param adj_mat_mnn a sparse matrix for mnn adjacent information.
 //' @param hyper a list of all hyperParameters.
 //' @param init a list of initialization of all parametres.
 //' @param platform string to idicate platform generating data.
 //' @param d dimension of hidden space.
 //' @param K number of subtypes/domains.
 //' @param adj_cutoff cutoff when calculating adjacent matrix.
 //' @param icm_maxiter number of iteration of icm.
 //' @param em_maxiter number of iteration of em.
 //' @param verbose print information or not.
 //' @param verbose_icm print icm information or not.
 //' @param ncores number of cores used in parallel.
 //' @return a list of all refined parameters.
 //'
 //' @export
 // [[Rcpp::export]]
 Rcpp::List stadia_EM_SP(Rcpp::NumericMatrix& x, Rcpp::NumericMatrix& position,
                        Rcpp::IntegerVector& batch_vec,
                        Rcpp::S4& adj_mat_mnn,
                        Rcpp::List hyper, Rcpp::List init,
                        std::string& platform,
                        const int& d = 15, const int& K = 7,
                        const double& adj_cutoff = 50., const int& icm_maxiter = 10,
                        const int& em_maxiter = 5, const bool& verbose = true,
                        const bool& verbose_icm = false,
                        const int& ncores=4) {
     // number of spots in all batches
     const uword n(x.ncol());

     // convert Rcpp to arma::mat object
     arma::mat xx(x.rows(), x.cols());
     xx = Rcpp::as<arma::mat>(x);
     arma::mat position_(position.begin(), position.nrow(), position.ncol(), false);

     // transform data type from Rcpp to arma
     arma::Col<int> batch_vec_(batch_vec.begin(), batch_vec.size(), false);
     arma::Col<int> aus_batch_vec = arma::unique(batch_vec_);
     arma::uword n_b = aus_batch_vec.n_elem;  // number of batches
     // adjacent matrix
     arma::sp_mat adj_mat(n,n);
     if(verbose)
         Rcpp::Rcout << "Calculating adjacent matrix..." << endl;
     adj_mat = calAdjacent(position_, batch_vec_, platform, adj_cutoff, ncores);
     adj_mat = adj_mat + convertSparse(adj_mat_mnn);

     // hyperparameters
     if(verbose)
         Rcpp::Rcout << "Extracting hyperparameters..." << endl;
     int n_Lambda = hyper["n_Lambda"];
     arma::vec eta_vec = hyper["eta"];
     double alpha_p = hyper["alpha_p"], beta_p = hyper["beta_p"];
     double lambda0 = hyper["lambda0"], lambda1 = hyper["lambda1"];
     double nu_tau = hyper["nu_tau"], nu_omega = hyper["nu_omega"];
     arma::vec mu_mu = hyper["mu_mu"];
     arma::mat Sigma_mu = hyper["Sigma_mu"];
     arma::mat Sigma_Lambda = hyper["Sigma_Lambda"];

     // initialization
     if(verbose)
         Rcpp::Rcout << "Extracting initializations..." << endl;
     arma::Col<int> c_current = init["c"];
     arma::mat L_current = init["L"];
     arma::mat gamma_current = init["gamma"];

     arma::mat B_current = init["B"];
     arma::mat F_current = init["F"];
     arma::mat T_current = init["T"];
     arma::vec omega_current = init["omega"];
     arma::mat mu_current = init["mu"];
     arma::mat Lambda_current = init["Lambda"];
     arma::vec p_current = init["p"];

     // traces
     arma::vec logpost(em_maxiter);

     // EM algorithms
     int em_iter;
     arma::mat reduction(d,n); // low-dimension representations
     arma::mat Egamma(L_current.n_rows, L_current.n_cols);
     double increaseRatio;
     Rcpp::Rcout << "Main algorithm..." << endl;
     Rcpp::Rcout << "==============================================" << endl;
     for(em_iter = 0; em_iter < em_maxiter; ++em_iter) {
         // ---------------------- E step ----------------------------

         // update C
         // Rcpp::Rcout <<  "update c" << endl;
         Rcpp::List icm_result = update_c_sp(xx, adj_mat, c_current, B_current, L_current, mu_current,
                                             T_current, Lambda_current, batch_vec_, omega_current,
                                             eta_vec, icm_maxiter, verbose_icm, ncores);
         c_current = Rcpp::as<arma::Col<int>>(icm_result["c_vec"]);
         // arma::vec eta_vec = icm_result["eta_vec"];
         arma::mat S_lik = icm_result["S_lik"];
         Rcpp::List inv_Phis = icm_result["inv_Phis"];
         double loglik = icm_result["logLikelihood"];
         // calculate expectation of latent variables *f*
         // Rcpp::Rcout <<  "calculate f" << endl;
         Rcpp::List Ef_result = expect_f(xx, inv_Phis, B_current, Lambda_current,
                                         mu_current, T_current, L_current, S_lik, omega_current,
                                         batch_vec_, ncores);
         arma::mat Ef_li = Ef_result["Ef_li"];
         reduction = Ef_li;
         arma::cube Efft_li = Ef_result["Efft_li"];
         arma::mat muEf = Ef_result["muEf"];
         arma::cube weightEf_liks = Ef_result["weightEf_liks"];
         arma::mat Ef_k = Ef_result["Ef_k"];
         // calculate expectation of latent variable *gamma*
         // Rcpp::Rcout <<  "calculate gamma" << endl;
         Egamma = expect_gamma(L_current, lambda0, lambda1, p_current, ncores);

         // ---------------------- M step ----------------------------
         // update T
         // Rcpp::Rcout <<  "update T" << endl;
         T_current = update_T(xx, B_current,  L_current, Efft_li, Ef_li,
                              batch_vec_, nu_tau, ncores);
         // update B
         // Rcpp::Rcout <<  "update B" << endl;
         B_current = update_B(xx, T_current, L_current, Ef_li,
                              batch_vec_, n_b, ncores);
         // update p
         // Rcpp::Rcout <<  "update P" << endl;
         p_current = update_p(Egamma, alpha_p, beta_p, ncores);
         // update L
         // Rcpp::Rcout <<  "update L" << endl;
         L_current = update_L(xx, B_current, Egamma, batch_vec_, Ef_li,
                              Efft_li, T_current, L_current, lambda0,
                              lambda1, ncores);
         // update Lambda
         // Rcpp::Rcout <<  "update Lambda" << endl;
         Lambda_current = update_Lambda(Efft_li, muEf, omega_current, S_lik,
                                        mu_current, n_Lambda, Sigma_Lambda, ncores);
         // update mu
         // Rcpp::Rcout <<  "update mu" << endl;
         mu_current = update_mu(S_lik, Lambda_current, mu_mu, Sigma_mu,
                                Ef_k, omega_current, ncores);
         // update omega
         // Rcpp::Rcout <<  "update omega" << endl;
         omega_current = update_omega(S_lik, mu_current, Lambda_current,
                                      Efft_li, weightEf_liks, nu_omega, ncores);

         // ---------------------- log posterior ------------------------
         // Rcpp::Rcout <<  "calculate log posterior" << endl;
         logpost(em_iter) = logPost(loglik, L_current, B_current, mu_current,
                 p_current, omega_current, T_current, Lambda_current,
                 lambda0, lambda1, mu_mu, Sigma_mu, nu_omega, nu_tau,
                 n_Lambda, Sigma_Lambda, alpha_p,  beta_p);

         if(em_iter == 0) {
             Rcpp::Rcout << "logPosterior"<< em_iter <<"th iteration: " << logpost(em_iter) << endl;
             continue;
         }
         increaseRatio = (logpost(em_iter) - logpost(em_iter-1))/std::abs(logpost(em_iter-1));
         if(verbose){
             Rcpp::Rcout << "logPosterior"<< em_iter <<"th iteration: " << logpost(em_iter);
             Rcpp::Rcout << " with increaseRatio: " <<  increaseRatio << endl;
         }
         if(abs(increaseRatio) < 1e-5) break;
     }

     return Rcpp::List::create(Rcpp::_["logpost"] = logpost(span(1, em_iter-1)),
                               Rcpp::_["L"] = L_current,
                               Rcpp::_["B"] = B_current,
                               Rcpp::_["T"] = T_current,
                               Rcpp::_["mu"] = mu_current,
                               Rcpp::_["omega_vec"] = omega_current,
                               Rcpp::_["Lambda_current"] = Lambda_current,
                               Rcpp::_["p_vec"] = p_current,
                               Rcpp::_["c_vec"] = c_current,
                               Rcpp::_["factors"] = reduction,
                               Rcpp::_["gamma"] = Egamma);

 }

// //' @title stadia_EM_SC
//  //' @description main function to run stadia algorithm for scRNAs
//  //'
//  //' @param x preprocessed count data.
//  //' @param batch_vec vector of batch information.
//  //' @param adj_mat_mnn a sparse matrix for mnn adjacent information.
//  //' @param hyper a list of all hyperParameters.
//  //' @param init a list of initialization of all parametres.
//  //' @param d dimension of hidden space.
//  //' @param K number of subtypes/domains.
//  //' @param icm_maxiter number of iteration of icm.
//  //' @param em_maxiter number of iteration of em.
//  //' @param verbose print information or not.
//  //' @param verbose_icm print icm information or not.
//  //' @param ncores number of cores used in parallel.
//  //' @return a list of all refined parameters.
//  //'
//  //' @export
//  // [[Rcpp::export]]
//  Rcpp::List stadia_EM_SC(Rcpp::NumericMatrix& x, Rcpp::IntegerVector& batch_vec,
//                         Rcpp::S4& adj_mat_mnn,
//                         Rcpp::List hyper, Rcpp::List init,
//                         const int& d = 15, const int& K = 7,
//                         const int& icm_maxiter=10, const int& em_maxiter = 5,
//                         const bool& verbose = true, const bool& verbose_icm = false,
//                         const int& ncores=4) {
//      // dimension
//      const uword n(x.ncol());
//
//      // Rcpp object to arma object
//      arma::mat xx(x.rows(), x.cols());
//      xx = Rcpp::as<arma::mat>(x);
//
//      // transform data type from Rcpp to arma
//      arma::Col<int> batch_vec_(batch_vec.begin(), batch_vec.size(), false);
//      arma::Col<int> aus_batch_vec = arma::unique(batch_vec_);
//      arma::uword n_b = aus_batch_vec.n_elem;
//      // arma::uword K = arma::conv_to<arma::uword>::from(K_);
//
//      arma::sp_mat adj_mat(n,n);
//      adj_mat = convertSparse(adj_mat_mnn);
//
//      // hyperparameters
//      if(verbose)
//          Rcpp::Rcout << "Extracting hyperparameters..." << endl;
//      arma::vec eta_vec = hyper["eta"];
//      int n_Lambda = hyper["n_Lambda"];
//      double alpha_p = hyper["alpha_p"], beta_p = hyper["beta_p"];
//      double lambda0 = hyper["lambda0"], lambda1 = hyper["lambda1"];
//      double nu_tau = hyper["nu_tau"], nu_omega = hyper["nu_omega"];
//      arma::vec mu_mu = hyper["mu_mu"];
//      arma::mat Sigma_mu = hyper["Sigma_mu"];
//      arma::mat Sigma_Lambda = hyper["Sigma_Lambda"];
//
//      // initialization
//      if(verbose)
//          Rcpp::Rcout << "Extracting initializations..." << endl;
//      arma::Col<int> c_current = init["c"];
//      arma::mat L_current = init["L"];
//      arma::mat gamma_current = init["gamma"];
//
//      arma::mat B_current = init["B"];
//      arma::mat F_current = init["F"];
//      arma::mat T_current = init["T"];
//      arma::vec omega_current = init["omega"];
//      arma::mat mu_current = init["mu"];
//      arma::mat Lambda_current = init["Lambda"];
//      arma::vec p_current = init["p"];
//
//      // traces
//      arma::vec logpost(em_maxiter);
//
//      // EM algorithms
//      int em_iter;
//      arma::mat reduction(d,n);
//      arma::mat Egamma(L_current.n_rows, L_current.n_cols);
//      //double increaseRatio_old = 100;
//      double increaseRatio;
//      Rcpp::Rcout << "Main algorithm..." << endl;
//      Rcpp::Rcout << "==============================================" << endl;
//      for(em_iter = 0; em_iter < em_maxiter; ++em_iter) {
//          // ---------------------- E step ----------------------------
//
//          // update C
//          // Rcpp::Rcout <<  "update c" << endl;
//          Rcpp::List icm_result = update_c_sp(xx, adj_mat, c_current, B_current, L_current, mu_current,
//                                              T_current, Lambda_current, batch_vec_, omega_current,
//                                              eta_vec, icm_maxiter, verbose_icm, ncores);
//          c_current = Rcpp::as<arma::Col<int>>(icm_result["c_vec"]);
//          arma::mat S_lik = icm_result["S_lik"];
//          Rcpp::List inv_Phis = icm_result["inv_Phis"];
//          double loglik = icm_result["logLikelihood"];
//          // calculate expectation of latent variables *f*
//          // Rcpp::Rcout <<  "calculate f" << endl;
//          Rcpp::List Ef_result = expect_f(xx, inv_Phis, B_current, Lambda_current,
//                                          mu_current, T_current, L_current, S_lik, omega_current,
//                                          batch_vec_, ncores);
//          arma::mat Ef_li = Ef_result["Ef_li"];
//          reduction = Ef_li;
//          arma::cube Efft_li = Ef_result["Efft_li"];
//          arma::mat muEf = Ef_result["muEf"];
//          arma::cube weightEf_liks = Ef_result["weightEf_liks"];
//          arma::mat Ef_k = Ef_result["Ef_k"];
//          // calculate expectation of latent variable *gamma*
//          // Rcpp::Rcout <<  "calculate gamma" << endl;
//          Egamma = expect_gamma(L_current, lambda0, lambda1, p_current, ncores);
//
//          // ---------------------- M step ----------------------------
//          // update T
//          // Rcpp::Rcout <<  "update T" << endl;
//          T_current = update_T(xx, B_current,  L_current, Efft_li, Ef_li,
//                               batch_vec_, nu_tau, ncores);
//          // update B
//          // Rcpp::Rcout <<  "update B" << endl;
//          B_current = update_B(xx, T_current, L_current, Ef_li,
//                               batch_vec_, n_b, ncores);
//          // update p
//          // Rcpp::Rcout <<  "update P" << endl;
//          p_current = update_p(Egamma, alpha_p, beta_p, ncores);
//          // update L
//          // Rcpp::Rcout <<  "update L" << endl;
//          L_current = update_L(xx, B_current, Egamma, batch_vec_, Ef_li,
//                               Efft_li, T_current, L_current, lambda0,
//                               lambda1, ncores);
//          // update Lambda
//          // Rcpp::Rcout <<  "update Lambda" << endl;
//          Lambda_current = update_Lambda(Efft_li, muEf, omega_current, S_lik,
//                                         mu_current, n_Lambda, Sigma_Lambda, ncores);
//          // update mu
//          // Rcpp::Rcout <<  "update mu" << endl;
//          mu_current = update_mu(S_lik, Lambda_current, mu_mu, Sigma_mu,
//                                 Ef_k, omega_current, ncores);
//          // update omega
//          // Rcpp::Rcout <<  "update omega" << endl;
//          omega_current = update_omega(S_lik, mu_current, Lambda_current,
//                                       Efft_li, weightEf_liks, nu_omega, ncores);
//
//          // ---------------------- log posterior ------------------------
//          // Rcpp::Rcout <<  "calculate log posterior" << endl;
//          logpost(em_iter) = logPost(loglik, L_current, B_current, mu_current,
//                  p_current, omega_current, T_current, Lambda_current,
//                  lambda0, lambda1, mu_mu, Sigma_mu, nu_omega, nu_tau,
//                  n_Lambda, Sigma_Lambda, alpha_p,  beta_p);
//
//          if(em_iter == 0) {
//              Rcpp::Rcout << "logPosterior"<< em_iter <<"th iteration: " << logpost(em_iter) << endl;
//              continue;
//          }
//          increaseRatio = (logpost(em_iter) - logpost(em_iter-1))/std::abs(logpost(em_iter-1));
//          if(verbose){
//              Rcpp::Rcout << "logPosterior"<< em_iter <<"th iteration: " << logpost(em_iter);
//              Rcpp::Rcout << " with increaseRatio: " <<  increaseRatio << endl;
//          }
//          if(abs(increaseRatio) < 1e-5) break;
//          // if(abs(increaseRatio) < 1e-5 || (logpost(em_iter-1) - logpost(em_iter)) > 1e-7) break;
//          //if(increaseRatio < 1e-5 || increaseRatio > increaseRatio_old) break;  // add condition
//          //increaseRatio_old = increaseRatio;
//
//
//      }
//
//      return Rcpp::List::create(Rcpp::_["logpost"] = logpost(span(1, em_iter-1)),
//                                Rcpp::_["L"] = L_current,
//                                Rcpp::_["B"] = B_current,
//                                Rcpp::_["T"] = T_current,
//                                Rcpp::_["mu"] = mu_current,
//                                Rcpp::_["omega_vec"] = omega_current,
//                                Rcpp::_["Lambda_current"] = Lambda_current,
//                                Rcpp::_["p_vec"] = p_current,
//                                Rcpp::_["c_vec"] = c_current,
//                                Rcpp::_["factors"] = reduction,
//                                Rcpp::_["gamma"] = Egamma);
//
//  }

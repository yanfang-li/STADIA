#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <omp.h>
#include "utils_em_intern.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::depends(RcppArmadillo, RcppDist)]]
//[[Rcpp::plugins(openmp)]]

// =============================================================================
//                        E-step in EM-based algorithm
// =============================================================================
// ########## icm to update c
// [[Rcpp::export]]
Rcpp::List update_c_sp(const arma::mat& x, const arma::sp_mat& adj_mat,
                       arma::Col<int>& c_vec, const arma::mat& B,
                       const arma::mat& L, const arma::mat& Mu,
                       const arma::mat& T, const arma::mat& Lambda,
                       const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                       arma::vec& eta_vec, const int& icm_maxiter = 10,
                       const bool& verbose = true, const int& ncores=4) {
    const int n(x.n_cols), K(Mu.n_cols);

    // allocate space for energys
    arma::mat energy_all(n,K);
    arma::mat energy_lat(n,K);
    arma::vec energy(icm_maxiter);

    // calculate energy corresponding to data
    Rcpp::List energy_obs_list;
    energy_obs_list = calEnergy_observed(x, B, L, Mu, T, Lambda, batch_vec, omega_vec, ncores);
    arma::mat energy_obs = energy_obs_list["energy_obs"];

    for(int i=0; i<icm_maxiter; ++i){
        // calculate energy corresponding to hat_c
        energy_lat = calEnergy_latent(adj_mat, c_vec, K, eta_vec, batch_vec);
        // overall energy
        energy_all = energy_lat + energy_obs;

        // updating c with lowest energy
        c_vec = conv_to< Col<int> >::from(index_min(energy_all, 1)) + 1;

        // full energy correspondin to new c: negative log likelihood
        energy(i) = sum(min(energy_all, 1));

        if(verbose) {
            Rcpp::Rcout << "energy of "<< i <<"th iteration: " << energy(i) << endl;
        }
        if(i==0) continue;
        if((energy(i-1) - energy(i)) < 1e-5 && verbose) {
            // Rcpp::Rcout << "ICM Converged at Iteration = " << i << endl;
            break;
        }
    }

    // calculate logLikelihood of x and c
    // to avoid obstacles to computation (log(small)=-infinity), some tricks are used like DR-SC algorithm
    arma::vec rowmin_energy = min(energy_all, 1);
    arma::mat aus_energy = energy_all - repmat(rowmin_energy, 1, K);
    arma::mat logLikMat = arma::exp(-aus_energy);
    double logLikelihood = arma::accu(arma::log(arma::sum(logLikMat, 1)) - rowmin_energy);  // sum_i \log \sum_k

    // calculate S_lik used later: probability of each spot to each subtype
    arma::mat S_lik(n,K);
    S_lik = arma::normalise(logLikMat, 1, 1);

    // update eta
    arma::Col<int> aus_batch_vec = arma::unique(batch_vec);
    arma::uword n_b = aus_batch_vec.n_elem;
    arma::vec etas = arma::regspace(0.1, 0.1,  4);
    arma::mat objEta(n_b, etas.n_elem);
    for(uword i=0; i < etas.n_elem; ++i){
        energy_lat = calEnergy_latent(adj_mat, c_vec, K, ones(n_b)*etas(i), batch_vec);
        for(uword j=0; j<n_b; ++j){
            objEta(j, i) = -arma::accu(S_lik.rows(find(batch_vec==(j+1))) % energy_lat.rows(find(batch_vec==(j+1))));
        }
    }
    eta_vec = etas(arma::index_max(objEta, 1));

    // return
    return Rcpp::List::create(Rcpp::_["c_vec"] = c_vec,
                              Rcpp::_["S_lik"] = S_lik,
                              Rcpp::_["eta_vec"] = eta_vec,
                              Rcpp::_["logLikelihood"] = logLikelihood,
                              Rcpp::_["objEta"] = objEta,
                              Rcpp::_["inv_Phis"] = energy_obs_list["inv_Phis"]);
}

// [[Rcpp::export]]
Rcpp::List update_c_sc(const arma::mat& x, const arma::mat& B,
                       const arma::mat& L, const arma::mat& Mu,
                       const arma::mat& T, const arma::mat& Lambda,
                       const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                       arma::mat& pk_mat, const int& ncores=4) { // pk_mat: bxK
    const int n(x.n_cols), K(Mu.n_cols);

    // calculate energy corresponding to data
    Rcpp::List energy_obs_list;
    energy_obs_list = calEnergy_observed(x, B, L, Mu, T, Lambda, batch_vec, omega_vec, ncores);
    arma::mat energy_obs = energy_obs_list["energy_obs"];
    arma::uvec batch_inx = arma::conv_to<arma::uvec>::from(batch_vec);
    arma::mat energy_all = energy_obs - arma::log(pk_mat.rows(batch_inx - 1));
    // arma::mat energy_all = energy_obs - repmat(log(pk_vec), 1, energy_obs.n_rows).t();

    // update c
    arma::Col<int> c_vec = conv_to< Col<int> >::from(index_min(energy_all, 1)) + 1;
    // calculate logLikelihood of x
    arma::vec rowmin_energy = min(energy_obs, 1);
    arma::mat aus_energy = energy_obs - repmat(rowmin_energy, 1, K);
    arma::mat logLikMat = arma::exp(-aus_energy);
    double logLikelihood = arma::accu(arma::log(arma::sum(logLikMat, 1)) - rowmin_energy);

    // calculate S_lik used later: probability of each spot to each subtype
    arma::mat S_lik(n,K);
    S_lik = arma::normalise(logLikMat, 1, 1);
    arma::mat sub_S;
    for(uword i = 0; i<pk_mat.n_rows; ++i) {
        sub_S = S_lik.rows(batch_vec == (i+1));
        pk_mat.row(i) = sum(sub_S, 0)/sub_S.n_rows;
    }
    // pk_vec = sum(S_lik.t(), 1)/n;

    // return
    return Rcpp::List::create(Rcpp::_["c_vec"] = c_vec,
                              Rcpp::_["S_lik"] = S_lik,
                              Rcpp::_["pk_mat"] = pk_mat,
                              Rcpp::_["logLikelihood"] = logLikelihood,
                              Rcpp::_["inv_Phis"] = energy_obs_list["inv_Phis"]);
}

// ##########expectation of latent variables f and ff^t
// [[Rcpp::export]]
Rcpp::List expect_f(const arma::mat& x, const Rcpp::List& inv_Phis_list,
                    const arma::mat& B, const arma::mat& Lambda,
                    const arma::mat& Mu, const arma::mat& T, const arma::mat& L,
                    const arma::mat& S_lik ,const arma::vec& omega_vec,
                    const arma::Col<int>& batch_vec,
                    const int& ncores) {
    using arma::uword;
    omp_set_num_threads(ncores);

    // dimensions and matrices used later
    const uword n(x.n_cols), K(Mu.n_cols), d(Lambda.n_cols), p(x.n_rows);
    arma::cube inv_Phis = inv_Phis_list["inv_Phis"];
    // arma::cube sqTLs = inv_Phis_list["sqTLs"];

    // allocates space for output
    arma::cube weightEf_liks(d,K,n); // S_lik Ef_lik: used in updating omega
    arma::mat Ef_li(d,n);      // sum_k S_lik Ef_lik
    arma::cube Efft_li(d,d,n); // sum_k S_lik Ef_lik f_lik^T
    arma::mat muEf(d, d);      // sum_i omega_i sum_k S_lik mu_k Ef_lik^T: used in update_Lambda
    arma::mat Ef_k(d, K);      // sum_i omega_i S_lik Ef_lik: used in update_mu

    // note that due to dependence of updating Lambda on mu & omega, and
    // updating mu on omega, the order of updating these three parameter is
    // Lambda -> mu -> omega

    // // aus matrices
    // arma::mat sqT(T.n_rows, T.n_cols);
    // sqT = arma::sqrt(T);  // T^1/2

    // iteration across all samples
#pragma omp parallel
{
    // private variables
    arma::mat Ef_lik(d,K);        // Ef_lik
    arma::mat weightEf_lik(d,K);  // S_lik Ef_lik
    arma::vec LTx(p);             // L^T T \tilde x
    arma::mat x_spec(d,K);        // \phi_{\ell i}
    arma::mat muEf_private(d, d); // omega_i sum_k S_lik mu_k Ef_lik^T
    arma::mat Ef_k_private(d, K); // omega_li S_lik Ef_lik

#pragma omp for nowait
    for(uword i=0; i<n; i++) {
        int batch_inx = batch_vec(i) - 1;  // batch information
        LTx = L.t() * (T.col(batch_inx) % (x.col(i) - B.col(batch_inx))); // L^T T \tilde x
        // LTx = sqTLs.slice(batch_inx).t() *arma::diagmat(sqT.col(batch_inx)) * (x.col(i) - B.col(batch_inx)); // dx1
        x_spec = repmat(LTx, 1, K) + omega_vec(i)*Lambda*Mu; // dxK: \phi_{\ell i}

        Ef_lik = inv_Phis.slice(i)*x_spec; // dxK: Ef_lik
        weightEf_lik = Ef_lik.each_row() % S_lik.row(i); // dxK: S_lik Ef_lik

        // output
        weightEf_liks.slice(i) = weightEf_lik;

        // outputs
        Ef_li.col(i) = sum(weightEf_lik, 1); // dx1: sum_k S_lik Ef_lik
        Efft_li.slice(i) = weightEf_lik * (Ef_lik.t()) + inv_Phis.slice(i); // dxd: sum_k S_lik Ef_lik f_lik^T

        // outputs
        muEf_private += omega_vec(i) * Mu * (weightEf_lik.t()); // omega_i sum_k S_lik mu_k Ef_lik^T
        Ef_k_private += omega_vec(i) * weightEf_lik;            // omega_i S_lik Ef_lik
    }
#pragma omp critical
    muEf += muEf_private; // dxd: sum_i omega_i sum_k S_lik mu_k Ef_lik^T
    Ef_k += Ef_k_private; // dx1: sum_i omega_i S_lik Ef_lik
}
return Rcpp::List::create(Rcpp::_["Ef_li"] = Ef_li,     // dxn
                          Rcpp::_["Efft_li"] = Efft_li, // dxdxn
                          Rcpp::_["muEf"] = muEf, // dxd
                          Rcpp::_["weightEf_liks"] = weightEf_liks, // dxKxn
                          Rcpp::_["Ef_k"] = Ef_k // dxK
);
}

// ########## expectation of latent variable gamma
// [[Rcpp::export]]
arma::mat expect_gamma(const arma::mat& L, const double& lambda0,
                       const double& lambda1, const arma::vec& p_vec,
                       const int& ncores) {
    omp_set_num_threads(ncores);
    const uword p(L.n_rows), d(L.n_cols);

    arma::mat hMat(p,d);
    // private variables
    double aus_cont1;
    double aus_cont2;
    double aus_cont3;
    arma::vec sqLj(p);

    aus_cont1 = exp(1.5*log(lambda1))/exp(0.5*log(lambda0));
    aus_cont2 = 0.5 * (1/lambda0 - 1/lambda1);

#pragma omp parallel for schedule(static) private(aus_cont3,sqLj)
    for(uword j=0; j<d; ++j) {
        aus_cont3 = aus_cont1 * (1-p_vec(j))/p_vec(j);
        // aus_cont3 = ((1-p_vec(j)) * exp(1.5*log(lambda1))) / (p_vec(j)*exp(0.5*log(lambda0)));
        // aus_cont2 = 0.5 * (1/lambda0 - 1/lambda1);
        sqLj = arma::square(L.col(j));
        hMat.col(j) = 1/(1+(aus_cont3/sqLj) % arma::exp(- aus_cont2* sqLj));
    }
    return hMat;
}

// =============================================================================
//                        M-step in EM-based algorithm
// =============================================================================
// update T
// [[Rcpp::export]]
arma::mat update_T(const arma::mat& x, const arma::mat& B,  const arma::mat& L,
                   const arma::cube& Efft_li, const arma::mat& Ef_li,
                   const arma::Col<int> batch_vec, const double& nu_tau,
                   const int& ncores) {
    omp_set_num_threads(ncores);
    const uword n_b(B.n_cols), p(x.n_rows), d(L.n_cols);

    // mat Ip(p, p, fill::eye); // identity matrix
    // arma::mat aus_T(p,p);
    arma::vec aus_T(p);
    arma::mat T(p, n_b);
    arma::mat aus_Efft(d,d);
    arma::mat aus_LsqEfft(p, d);

#pragma omp parallel for schedule(static) private(aus_Efft, aus_T, aus_LsqEfft)
    for(uword l=0; l<n_b; ++l) {
        arma::uvec inx_b = find(batch_vec == (l+1));
        int n_l = inx_b.n_elem;

        aus_Efft = sum(Efft_li.slices(inx_b), 2); // dxd
        arma::mat tildeX = x.cols(inx_b) - repmat(B.col(l), 1, n_l); // pxn_l
        arma::mat aus_Ef_l = Ef_li.cols(inx_b); // dxn_l

        // aus_T = tildeX  * (tildeX.t()) + L * (aus_Efft) * L.t() -
        //       tildeX * (aus_Ef_l.t()) * L.t()  - L * aus_Ef_l * (tildeX.t()) + nu_tau*Ip;
        // T1.col(l) = (n_l + nu_tau - 2) * (1/aus_T.diag());

        arma::mat U, V;
        arma::vec s;
        svd(U, s, V, aus_Efft, "dc");
        aus_LsqEfft = L * (U * arma::diagmat(arma::sqrt(s))); // pxd

        aus_T = arma::sum(tildeX % tildeX, 1)  -
            2 * arma::sum(tildeX % (L * aus_Ef_l), 1) +
            nu_tau * arma::ones<vec>(p);//px1
        aus_T += arma::sum(aus_LsqEfft % aus_LsqEfft, 1);
        // arma::mat tmp = L * (aus_Efft) * L.t();
        // aus_T += tmp.diag();
        T.col(l) = (n_l + nu_tau - 2) * (1/aus_T);
    }
    return T;
}

// update B
// [[Rcpp::export]]
arma::mat update_B(const arma::mat& x, const arma::mat& T,
                   const arma::mat& L, const arma::mat& Ef_li,
                   const arma::Col<int>& batch_vec, const arma::uword& n_b,
                   const int& ncores) {
    omp_set_num_threads(ncores);
    const uword p(x.n_rows);

    mat Ip(p, p, fill::eye);
    arma::mat B(p, n_b);
    arma::mat Tl(p,p);

#pragma omp parallel for schedule(static) private(Tl)
    for(uword l=0; l<n_b; ++l){
        arma::uvec inx_b = find(batch_vec == (l+1));
        uword n_l = inx_b.n_elem;

        // Tl = arma::diagmat(T.col(l)); // pxp
        // B.col(l) = inv_sympd(n_l * Tl + Ip) * (Tl * sum(x.cols(inx_b), 1)- Tl * L * sum(Ef_li.cols(inx_b), 1));
        B.col(l) = (1/(n_l * T.col(l) + 1)) % (T.col(l) % (sum(x.cols(inx_b), 1)- (L * sum(Ef_li.cols(inx_b), 1))));
    }

    return B;
}

// update Mu
// [[Rcpp::export]]
arma::mat update_mu(const arma::mat& S_lik, const arma::mat& Lambda,
                    const arma::vec& mu_mu, const arma::mat& sigma_mu,
                    const arma::mat& Ef_k, const arma::vec& omega_vec,
                    const int& ncores) {
    omp_set_num_threads(ncores);
    const int d(Lambda.n_cols), K(S_lik.n_cols);

    arma::mat mu(d,K);

    arma::mat aux_mu(d,d);
#pragma omp parallel for schedule(static) private(aux_mu)
    for(int k=0; k<K; ++k) {
        aux_mu = inv_sympd(accu(S_lik.col(k) % omega_vec) * Lambda + inv_sympd(sigma_mu)); // dxd
        mu.col(k) = aux_mu*(Lambda * Ef_k.col(k) + inv_sympd(sigma_mu)*mu_mu); // dx1
    }

    return mu;
}

// update omega
// [[Rcpp::export]]
arma::vec update_omega(const arma::mat& S_lik, const arma::mat& Mu,
                       const arma::mat& Lambda, const arma::cube& Efft_li,
                       const arma::cube& weightEf_liks, const double& nu_omega,
                       const int& ncores) {
    omp_set_num_threads(ncores);
    int K(S_lik.n_cols), n(S_lik.n_rows), d(Lambda.n_cols);

    arma::vec aus_mu(K);
    for(int k=0; k<K; ++k) {
        aus_mu(k) = arma::as_scalar((Mu.col(k).t()) * Lambda * (Mu.col(k)));
    }
    arma::vec Smu(n);
    Smu = S_lik * aus_mu;

    arma::vec omega(n);
    arma::vec muLamEf_li(n);
    arma::mat weightEf_lik(d, K);
#pragma omp parallel for schedule(static) private(weightEf_lik)
    for(int i=0; i<n; ++i) {
        weightEf_lik = weightEf_liks.slice(i);
        for(int j=0; j<K; ++j) {
            // 1x1: sum_k S_lik mu_k^t Lambda Ef_lik
            muLamEf_li(i) += arma::as_scalar(Mu.col(j).t() * Lambda * (weightEf_lik.col(j)));
        }
        omega(i) = (d+nu_omega - 2) / (Smu(i) + arma::trace(Lambda*Efft_li.slice(i)) -
            2*muLamEf_li(i) + nu_omega);
    }
    return omega;
}

// update Lambda
// [[Rcpp::export]]
arma::mat update_Lambda(const arma::cube& Efft_li, const arma::mat& muEf,
                        const arma::vec& omega_vec, const arma::mat& S_lik,
                        const arma::mat& Mu, const int& n_Lambda,
                        const arma::mat& Sigma_Lambda, const int& ncores) {
    omp_set_num_threads(ncores);
    const uword n(S_lik.n_rows), d(Sigma_Lambda.n_rows);

    arma::mat inv_Lambda(d,d);
#pragma omp parallel
{
    // private variables
    arma::mat inv_Lambda_private(d,d);
#pragma omp for nowait
    for(uword i = 0; i < n; ++i) {
        inv_Lambda_private += omega_vec(i) * (Efft_li.slice(i) +
            (Mu.each_row() % S_lik.row(i)) * Mu.t()); // dxd
    }
#pragma omp critical
    inv_Lambda += inv_Lambda_private;
}
inv_Lambda += -muEf - muEf.t() + inv_sympd(Sigma_Lambda);
arma::mat Lambda = inv_sympd(1.0/(n+n_Lambda - d - 1)*inv_Lambda);
return Lambda;
}

// update p
// [[Rcpp::export]]
arma::vec update_p(const arma::mat& Egamma, const double& alpha_p,
                   const double& beta_p, const int& ncores) {
    omp_set_num_threads(ncores);
    const uword d(Egamma.n_cols), p(Egamma.n_rows);

    arma::vec p_vec(d);
#pragma omp parallel for schedule(static)
    for(uword j=0; j<d; ++j) {
        // p_vec(j) = (arma::accu(Egamma.col(j)) + alpha_p/(j+1) - 1) / (p + alpha_p/(j+1) + beta_p - 2);
        p_vec(j) = (arma::accu(Egamma.col(j)) + alpha_p - 1) / (p + alpha_p + beta_p - 2);
    }

    return p_vec.clamp(0,1);
}

// update L
// [[Rcpp::export]]
arma::mat update_L(const arma::mat& x, const arma::mat& B,
                   const arma::mat& Egamma, const arma::Col<int>& batch_vec,
                   const arma::mat& Ef_li, const arma::cube& Efft_li,
                   const arma::mat& T, const arma::mat& L,
                   const double& lambda0, const double& lambda1,
                   const int& ncores) {
    omp_set_num_threads(ncores);
    arma::mat A(size(Egamma)), B_(size(Egamma)), C(size(Egamma));

    arma::mat tilX = x - contB2mat(B, batch_vec);//pxn
    arma::mat aus_T = contT2mat(T, batch_vec);   // pxn
    arma::mat aus_Efft = extractDiag(Efft_li);   // dxn
    arma::mat aus_TX = aus_T % tilX; //pxn
    const uword d(L.n_cols), p(L.n_rows);

    arma::mat aus_A(p,d);
    arma::rowvec aus_B(batch_vec.n_elem);
    B_ = aus_TX * Ef_li.t();

    // 尝试
#pragma omp parallel for schedule(static)
    for(uword j = 0; j < d; ++j) {
        B_.col(j) = B_.col(j) - sum(aus_T % (L * extractJcol(Efft_li, j)), 1);
    }

    aus_A = aus_T * aus_Efft.t();
    A = -(Egamma/lambda1 + (1-Egamma)/lambda0 + aus_A); // pxd
    C = 2*Egamma; // pxd

    arma::mat out = solveQ(A, B_, C);
    return out;
}
// // [[Rcpp::export]]
// arma::mat update_L(const arma::mat& x, const arma::mat& B,
//                    const arma::mat& Egamma, const arma::Col<int>& batch_vec,
//                    const arma::mat& Ef_li, const arma::cube& Efft_li,
//                    const arma::mat& T, const arma::mat& L,
//                    const double& lambda0, const double& lambda1,
//                    const int& ncores) {
//       omp_set_num_threads(ncores);
//       arma::mat A(size(Egamma)), B_(size(Egamma)), C(size(Egamma));
//
//       arma::mat tilX = x - contB2mat(B, batch_vec);//pxn
//       arma::mat aus_T = contT2mat(T, batch_vec);   // pxn
//       arma::mat aus_Efft = extractDiag(Efft_li);   // dxn
//       arma::mat aus_TX = aus_T % tilX; //pxn
//       const uword d(L.n_cols), p(L.n_rows);
//
//       arma::mat aus_A(p,d);
// #pragma omp parallel for schedule(static)
//       for(uword i = 0; i < p; ++i) {
//             for(uword j = 0; j < d; ++j) {
//                   aus_A(i,j) = accu(aus_Efft.row(j) % aus_T.row(i));
//             }
//       }
//
//       A = -(Egamma/lambda1 + (1-Egamma)/lambda0 + aus_A); // pxd
//       C = 2*Egamma; // pxd
//
//       arma::rowvec aus_B(batch_vec.n_elem);
// #pragma omp parallel for schedule(static) private(aus_B)
//       for(uword i = 0; i < p; ++i) {
//             for(uword j = 0; j < d; ++j) {
//                   aus_B = L.row(i) * extractJcol(Efft_li, j);//1xn
//                   B_(i,j) = arma::accu(aus_TX.row(i) % Ef_li.row(j) - aus_T.row(i) % aus_B);//sum(1xn)
//             }
//       }
//
//       arma::mat out = solveQ(A, B_, C);
//       return out;
// }

// =============================================================================
//                        log posterior
// =============================================================================
// prior of L
//[[Rcpp::export]]
double logPrior_L(const arma::mat& L, const arma::vec& p_vec,
                  const double& lambda0, const double& lambda1) {
    double out;
    arma::mat aus_p = repmat(p_vec, 1, L.n_rows).t();

    // out = arma::accu(arma::log(aus_p/lambda1 % (arma::square(L) % arma::normpdf(L, 0, std::sqrt(lambda1))) +
    // (1-aus_p) % arma::normpdf(L, 0, std::sqrt(lambda0))));
    out = arma::accu(aus_p %  (arma::log(arma::square(L)/lambda1) +
        arma::log_normpdf(L, 0, std::sqrt(lambda1))) +
        (1-aus_p) % arma::log_normpdf(L, 0, std::sqrt(lambda0)));
    return out;
}

// prior of B
//[[Rcpp::export]]
double logPrior_B(const arma::mat& B) {
    double out;
    out = arma::accu(arma::log_normpdf(B, 0, 1));
    return out;
}

// prior of mu
//[[Rcpp::export]]
double logPrior_mu(const arma::mat& mu, const arma::vec& mu_mu,
                   const arma::mat& Sigma_mu) {
    double out=0;
    arma::uword K(mu.n_cols), d(mu.n_rows);
    double aus_const = - (d/2.0*std::log(2.0*M_PI)+0.5*arma::log_det_sympd(Sigma_mu));

    out = K*aus_const;
    for(uword j = 0; j < K; ++j){
        out+= - 0.5 * arma::as_scalar((mu.col(j) - mu_mu).t() * arma::inv_sympd(Sigma_mu) * (mu.col(j) - mu_mu));
    }

    return out;
}

// omega
//[[Rcpp::export]]
double logPrior_omega(const arma::vec& omega_vec, const double& nu_omega) {
    double out;
    double aus_const = 0.5*nu_omega;
    Rcpp::NumericVector aus_omega = Rcpp::NumericVector(omega_vec.begin(), omega_vec.end());

    out = Rcpp::sum(Rcpp::dgamma(aus_omega, aus_const, aus_const, true));
    return out;
}

// T
//[[Rcpp::export]]
double logPrior_T(const arma::mat& T, const double& nu_tau) {
    double out;
    double aus_const = 0.5*nu_tau;
    Rcpp::NumericMatrix aus_T = Rcpp::wrap(T);

    out = Rcpp::sum(Rcpp::dgamma(aus_T + 1e-5, aus_const, aus_const, true));
    return out;

}

// Lambda
//[[Rcpp::export]]
double logPrior_Lambda(const arma::mat& Lambda, const int& n_Lambda,
                       const arma::mat& Sigma_Lambda) {
    double out;
    out = dwish(Lambda, n_Lambda, Sigma_Lambda, true);
    return out;
}

// prior of p
//[[Rcpp::export]]
double logPrior_p(const arma::vec& p_vec, const double& alpha_p,
                  const double& beta_p) {
    double out=0;
    arma::uword d(p_vec.n_elem);

    for(uword j = 0; j < d; ++j) {
        // out += R::dbeta(p_vec(j), alpha_p/(j+1), beta_p, true);
        out += R::dbeta(p_vec(j), alpha_p, beta_p, true);
    }

    return out;
}

// log posterior
//[[Rcpp::export]]
double logPost(const double& loglik, const arma::mat& L,
               const arma::mat& B, const arma::mat& mu,
               const arma::vec& p_vec, const arma::vec& omega_vec,
               const arma::mat& T, const arma::mat& Lambda,
               const double& lambda0, const double& lambda1,
               const arma::vec& mu_mu,const arma::mat& Sigma_mu,
               const double& nu_omega, const double& nu_tau,
               const int& n_Lambda, const arma::mat& Sigma_Lambda,
               const double& alpha_p, const double& beta_p
) {

    double logpriorL = logPrior_L(L, p_vec, lambda0, lambda1);
    double logpriorB = logPrior_B(B);
    double logpriorMu = logPrior_mu(mu, mu_mu, Sigma_mu);
    double logpriorOmega = logPrior_omega(omega_vec, nu_omega);
    double logpriorT = logPrior_T(T, nu_tau);
    double logpriorLambda = logPrior_Lambda(Lambda, n_Lambda, Sigma_Lambda);
    double logpriorP = logPrior_p(p_vec, alpha_p, beta_p);

    double logpost = loglik + logpriorL + logpriorB + logpriorMu +
        logpriorOmega + logpriorT + logpriorLambda + logpriorP;
    return logpost;
}

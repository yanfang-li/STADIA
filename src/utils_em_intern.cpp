#include <RcppArmadillo.h>
#include <omp.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// =============================================================================
//                      calculate intern matrix used in E-step
// =============================================================================
// matrices used in E.f and E.ff^T: overall function
//[[Rcpp::export]]
Rcpp::List calPhi(const arma::mat& T, const arma::mat& L,
                  const arma::mat& Lambda, const arma::vec& omega_vec,
                  const arma::Col<int>& batch_vec) {
      // dimensions
      const uword n(omega_vec.n_elem), p(L.n_rows), d(L.n_cols), b(T.n_cols);

      // firstly calculate parts shared by Phi's
      arma::cube sqTLs(p,d,b);
      arma::cube LtTLs(d,d,b);
      arma::cube inv_Phis(d,d,n);

      for(uword i=0; i<b; ++i) {
            sqTLs.slice(i) = L;
            arma::vec v = arma::sqrt(T.col(i));
            sqTLs.slice(i).each_col() %= v; // T^{1/2}*L pxd
            LtTLs.slice(i) = sqTLs.slice(i).t()* sqTLs.slice(i); // L^tTL dxd
      }

      // calculate Phi's
      for(uword i=0; i<n; ++i) {
            inv_Phis.slice(i) = inv_sympd(LtTLs.slice(batch_vec(i)-1) + omega_vec(i)*Lambda);
      }
      return Rcpp::List::create(Rcpp::_["sqTLs"] = sqTLs,        // T^{1/2}*L pxd
                                Rcpp::_["inv_Phis"] = inv_Phis); // (L^tTL + omega*Lambda)^{-1}
}

// =============================================================================
//                      energy of each spot used in icm
// =============================================================================
// // inline function used in updating c: calculating energy corresponding to hat_c
// //[[Rcpp::export]]
// arma::rowvec calEnergy_latent_each(const arma::sp_mat& adj_mat,
//                                    const arma::Col<int>& c_vec,
//                                    const int& K, const int& j,
//                                    const double& eta) {
//       using arma::uword;
//       // allocate space for output: row for spot and column for subtype
//       arma::rowvec energy_latent(K);
//
//       // read-only iterator of adj_mat: all non-zero elements: j from 1
//       arma::sp_mat::const_col_iterator it = adj_mat.begin_col(j-1);
//       for(; it!=adj_mat.end_col(j-1); ++it) {
//             energy_latent(c_vec(it.row())-1) -=eta;
//       }
//       // log likelihood of latent subtype (normalized)
//       arma::rowvec energy_latent_final = -arma::log(exp(-energy_latent)/arma::sum(exp(-energy_latent)));
//       return energy_latent_final; // energy = - log(loglik_c)
// }
//

//[[Rcpp::export]]
arma::mat calEnergy_latent(const arma::sp_mat& adj_mat, const arma::Col<int>& c_vec,
                           const int& K, const arma::vec& eta_vec,
                           const arma::Col<int>& batch_vec) {
      using arma::uword;
      // allocate space for output: row for spot and column for subtype
      arma::mat energy_latent(adj_mat.n_rows, K);  // nxK

      // read-only iterator of adj_mat: all non-zero elements
      const uword n(adj_mat.n_nonzero);                  // no. of nonzero elements
      arma::sp_mat::const_iterator it = adj_mat.begin(); // iterator
      for(uword i=0; i<n; ++i) {
            energy_latent(it.row(), c_vec(it.col())-1) -=eta_vec(batch_vec(it.row())-1);
            ++it;
      }

      // negative log likelihood of latent subtype (normalized)
      arma::mat energy_latent_final = -arma::log(arma::normalise(exp(-energy_latent), 1, 1));
      return energy_latent_final; // energy = - log(loglik_c)
}

// ========================== observed energy of each spot =====================
// inline function used in updating c: calculating energy corresponding to data
//[[Rcpp::export]]
Rcpp::List calEnergy_observed(const arma::mat& x, const arma::mat& B, const arma::mat& L,
                              const arma::mat& Mu, const arma::mat& T, const arma::mat& Lambda,
                              const arma::Col<int>& batch_vec, const arma::vec& omega_vec,
                              const int& ncores=4) {
      using arma::uword;
      omp_set_num_threads(ncores);
      // dimensions: no. of spots; orig dim; no. of subtypes; low dim; no. of batches
      const uword n(x.n_cols), p(x.n_rows), K(Mu.n_cols), d(L.n_cols), b(B.n_cols);

      // declare U s V
      arma::mat U, V;
      arma::vec s;

      // extract matrices to use
      Rcpp::List Phis = calPhi(T, L, Lambda, omega_vec, batch_vec);
      arma::cube sqTLs = Phis["sqTLs"];       // T^1/2 L
      arma::cube inv_Phis = Phis["inv_Phis"]; //(L^T T_\ell L + \omega_ \Lambda)^{-1}

      // calculate utilities for log|Psi_{\ell i}|
      // from Woodbury formula:
      // Psi_{\ell i} = T^1/2 (I + \omega_{\ell i}^{-1} T^1/2 L \Lambda^{-1} L^T T^1/2)^{-1} T^1/2
      svd(U, s, V, Lambda.i());

      arma::mat aus_sqTLUsqS(d,d);
      arma::mat sqTLUsqS(p, d);
      aus_sqTLUsqS = U * arma::diagmat(arma::sqrt(s));
      arma::mat aus_logDetPsi_D2(d, b);
      arma::vec aus_logDetPsi_logDetT(b);
      for(uword i = 0; i < b; ++i){
            sqTLUsqS = sqTLs.slice(i) * aus_sqTLUsqS; // pxd: T^1/2 L U sqrt(S)
            aus_logDetPsi_D2.col(i) = arma::square(svd(sqTLUsqS));     // dx1
            aus_logDetPsi_logDetT(i) = arma::accu(arma::log(T.col(i)));// log(det(T_ell))
      }

      // calculate (x - Bm - Lmu) ~ N(0, Psi^{-1})
      arma::vec logDetPsi(n);
      arma::mat energy_obs(n,K);
      energy_obs.fill((p/2.0)*std::log(2.0*M_PI));
#pragma omp parallel for schedule(static) private(U, V, sqTLUsqS, s)
      for(uword i=0; i<n; ++i) {
            // log(det(I + \omega_{\ell i}^{-1} T^1/2 L \Lambda^{-1} L^T T^1/2))
            double aus_logDetPsi_logDetD2 =
                  arma::accu(arma::log(1+(1.0/omega_vec(i))*aus_logDetPsi_D2.col(batch_vec(i)-1)));
            // log(det(Psi_{\ell i}))
            logDetPsi(i) = aus_logDetPsi_logDetT(batch_vec(i) - 1) - aus_logDetPsi_logDetD2;

            svd(U, s, V, inv_Phis.slice(i));
            // T^1/2 L inv_Phis^{1/2}:T^1/2 L U S^1/2
            sqTLUsqS = sqTLs.slice(batch_vec(i) - 1) * (U * arma::diagmat(arma::sqrt(s)));
            arma::vec z(p);
            arma::vec zz(d);
            // p/2 log(2*pi) - 1/2*log(det(Psi_{\ell i}))
            energy_obs.row(i) -= 0.5*logDetPsi(i);
            for(uword k = 0; k < K; k++) {
                  // T^1/2 (\tilde x - L \mu_k)  px1
                  z = (arma::sqrt(T.col(batch_vec(i) - 1))) % (x.col(i) - B.col(batch_vec(i)-1) - L*Mu.col(k));
                  // z = arma::diagmat(arma::sqrt(T.col(batch_vec(i) - 1))) *
                  //       (x.col(i) - B.col(batch_vec(i)-1) - L*Mu.col(k));
                  // S^1/2 U^T L^T T^1/2 z   dx1
                  zz = sqTLUsqS.t() * z;
                  // calculate 1/2*(x - Bm - Lmu)^t Psi (x - Bm - Lmu)
                  energy_obs(i, k) += 0.5*(arma::dot(z, z) - arma::dot(zz, zz));
                  // // calculate 1/2*(x - Bm - Lmu)^t Psi (x - Bm - Lmu)
                  // energy_obs(i, k) = 0.5*(arma::dot(z, z) - arma::dot(zz, zz));
                  // // calculate final energy = - log(liklihood(x))
                  // energy_obs(i, k) += (p/2.0)*std::log(2.0*M_PI)  - 0.5*logDetPsi(i);
            }
      }

      return Rcpp::List::create(Rcpp::_["energy_obs"] = energy_obs,
                                Rcpp::_["inv_Phis"] = Phis);
}

// =============================================================================
//                      Extract elements from cubs
// =============================================================================
arma::mat extractDiag(const arma::cube& Efft_li) {
      const uword n(Efft_li.n_slices), d(Efft_li.n_rows);

      arma::mat  out(d, n);
      for(uword i=0; i<n; ++i) {
            out.col(i) = Efft_li.slice(i).diag();
      }
      return out;
}

arma::mat extractJcol(const arma::cube& Efft_li, const int& j) {
      arma::mat out(Efft_li.n_rows, Efft_li.n_slices);
      out = Efft_li(span::all, span(j), span::all); // dxn
      // out.shed_row(j);
      out.row(j) = zeros<arma::rowvec>( out.n_cols ); // dxn with jth row == 0
      return out;
}

// =============================================================================
//                      convert T/B to pxn matrix
// =============================================================================
arma::mat contT2mat(const arma::mat& T,         // pxn_b
                    const arma::Col<int>& batch_vec) {  // n
      const uword n_b(T.n_cols);

      arma::mat out(T.n_rows, batch_vec.n_elem); //pxn
      for(uword i = 0; i < n_b; ++i) {
            arma::uvec inx = arma::find(batch_vec == (i+1));
            out.cols(inx) = repmat(T.col(i), 1, inx.n_elem);
      }

      return out;
}

arma::mat contB2mat(const arma::mat& B,         // pxn_b
                    const arma::Col<int>& batch_vec) {  // n
      const uword n_b(B.n_cols);

      arma::mat out(B.n_rows, batch_vec.n_elem); //pxn
      for(uword i = 0; i < n_b; ++i) {
            arma::uvec inx = arma::find(batch_vec == (i+1));
            out.cols(inx) = repmat(B.col(i), 1, inx.n_elem);
      }

      return out;
}

// =============================================================================
//                      solver in updating L in M-step
// =============================================================================
arma::mat solveQ(const arma::mat& A, const arma::mat& B,
                 const arma::mat& C) {

      arma::mat out(size(B));
      arma::mat aus_out(size(B));

      arma::mat Delta = arma::square(B) - 4*(A%C);
      out = -B/(2*A);
      aus_out = arma::sqrt(Delta)/(2*A);

      arma::uvec posInx = find(B>=0);
      arma::uvec negInx = find(B<0);

      out(posInx) = out(posInx) - aus_out(posInx);
      out(negInx) = out(negInx) + aus_out(negInx);

      return out;
}

// arma::mat solveQ(const arma::mat& A, const arma::mat& B,
//                  const arma::mat& C) {
//
//       arma::mat out(size(B));
//
//       arma::mat Delta = arma::square(B) - 4*(A%C);
//
//       arma::uvec posInx = find(B>=0);
//       arma::uvec negInx = find(B<0);
//
//       out(posInx) = (-B(posInx) - arma::sqrt(Delta(posInx)))/(2*A(posInx));
//       out(negInx) = (-B(negInx) + arma::sqrt(Delta(negInx)))/(2*A(negInx));
// //
// //       omp_set_num_threads(ncores);
// // #pragma omp parallel for schedule(static)
// //       for(uword i: posInx) {
// //               out(i) = (-B(i) - std::sqrt(exp(2*log(B(i))) - 4*A(i)*C(i)))/(2*A(i));
// //       }
// // #pragma omp parallel for schedule(static)
// //       for(uword i: negInx) {
// //               out(i) = (-B(i) + std::sqrt(exp(2*log(B(i))) - 4*A(i)*C(i)))/(2*A(i));
// //       }
//
//       return out;
// }

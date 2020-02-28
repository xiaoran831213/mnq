#define ARMA_NO_DEBUG 
#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
#include <stdio.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// MINQUE-0
//
// y vector of response variable
// V list of kernels
// w vector of initial weights (variance components).
// X matrix of covariate
SEXP amo_mn0(SEXP _Y, SEXP _K, SEXP _X, SEXP _opt)
{
    // // [[Rcpp::export]]
    fmat Y = as<fmat>(_Y);	// column vector of response
    int N = Y.n_rows;		// sample size
    List ret;			// list of results
    List opt = as<List>(_opt);  // options

    /* list of kernels */
    List _l = as<List>(_K);	// list of kernels
    int M = _l.length();	// number of kernels
    std::vector<fmat> K(M);
    for(int i = 0; i < M; i++)
	K[i] = as<fmat>(_l[i]);

    /* covariate */
    fmat X = as<fmat>(_X);
    // Rcpp::Rcout << "X=" << X << std::endl;

    // w_0 = 1, w_1 = 0, ... , w_L = 0
    // V = w_0 * K_1 + w_1 * K_2 ... w_L * K_L = I  (since K_0 = I)

    // Q = V^{-1} - V^{-1}X (X'V^{-1}X)^{+} X'V^{-1}
    // get Q K_1, Q K_2, ... Q K_L, and Qy

    /* Get P, and Q = I - P, then R V_i and R Y, where R = V^Q */
    std::vector<fmat> RV(M);
    fmat Ry;
    fmat B, C, P;
    if(X.n_cols == 0)
    {
        /* No X, Q = V^{-1} = I ==> Q K_l = K_l, l = 1 ... L */
	for(int i = 0; i < M; i++)
	    RV[i] = K[i];
	Ry = Y;
    }
    else
    {
        // Q = V^{-1} - V^{-1}X (X'V^{-1}X)^+ X'V^{-1}
        //   = I - X (X'X)^{+} X'           (V^{-1} = I)
        //
        // Q K = K - X (X'X)^{+} (X' K)     (K_1 ... K_L)
        //
        // Q <- diag(N) - X %*% .ginv(t(X) %*% X) %*% t(X)
	B = X * pinv(X.t() * X);
	for(int i = 0; i < M; i++)
	    RV[i] = K[i] - B * (X.t() * K[i]);
	Ry = Y - B * (X.t() * Y);
    }
    // Rcpp::Rcout << "Ry=" << Ry << std::endl;

    // get y' Q K_l Q y = (Qy)' K_l (Qy), for l=0 ... L
    fmat u(M, 1);
    for(int i = 0; i < M; i++)
	u[i] = as_scalar((Ry.t() * K[i] * Ry));
    // Rcpp::Rcout << "u=" << u << std::endl;
	    
    // get Tr(QK_l QK_m) = | (QK_l)*(QK_m)'|_f
    fmat F(M, M);
    for(int i = 0; i < M; i++)
    {
	for(int j = i + 1; j < M; j++)
	{
	    // F(i, j) = accu(RV[i] % RV[j].t());
            F(i, j) = trace(RV[i] * RV[j]);
	    F(j, i) = F(i, j);
	}
	F(i, i) = accu(RV[i] % RV[i]);
    }
    // Rcpp::Rcout << "F=" << F << std::endl;

    /* solve:
    [Tr(QK_0 QK_0) ... Tr(QK_0 QK_L)]  [s_0^2] = [y'Q K_0 Q y]
    [Tr(QK_1 QK_0) ... Tr(QK_1 QK_L)]  [s_1^2] = [y'Q K_1 Q y]
    ...            ...            ...    ...          ...
    [Tr(QK_L QK_0) ... Tr(QK_L QK_L)]  [s_L^2] = [y'Q K_L Q y]
    for s_0^2, s_1^2, ... s_L^2 */
    fmat w(M, 1);
    w = pinv(F) * u;

    int gls = as<int>(opt["gls"]);
    // Rcpp::Rcout << "GLS=" << gls << std::endl;
    if(X.n_cols > 0 && gls > 0)
    {
        if(gls > 1)
        {
            fmat V(K[0] * w[0]);
            for(int i = 1; i < M; i++)
                V += K[i] * w[i];
            B = solve(V, X, solve_opts::fast + solve_opts::likely_sympd);
            ret["fix"] = pinv(X.t() * B) * B.t() * Y;
        }
        else
        {
            fmat D(N, 1, fill::zeros);
            for(int i = 0; i < M; i++)
                D += w[i] * K[i].diag();
            fmat C = X.each_col() / D;
            // Rcpp::Rcout << "D=" << D << std::endl;
            // Rcpp::Rcout << "C=" << C << std::endl;
            ret["fix"] = pinv(X.t() * C) * C.t() * Y;
        }
    }
    ret["vcs"] = w;
    ret["err"] = 0;
    return(Rcpp::wrap(ret));
}

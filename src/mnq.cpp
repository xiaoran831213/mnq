#define ARMA_NO_DEBUG 
#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
#include <stdio.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// MINQUE
//
// y vector of response variable
// V list of kernels
// w vector of initial weights (variance components).
// X matrix of covariate
SEXP amo_mnq(SEXP _Y, SEXP _K, SEXP _X, SEXP _w, SEXP _opt)
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

    /* initial weights */
    fmat w = as<fmat>(_w);
    if(w.n_cols == 0)
	w = fmat(M, 1, fill::ones);
    // Rcpp::Rcout << "w=" << w << std::endl;

    /* covariate */
    fmat X = as<fmat>(_X);
    // Rcpp::Rcout << "X=" << X << std::endl;

    // V: weighted sum of kernels, and A=V^ (V inverse)
    fmat V(N, N, fill::zeros);
    for(int i = 0; i < M; i++)
	V += K[i] * w[i];
    // Rcpp::Rcout << "V = " << V << std::endl;
    fmat A;
    if(!inv_sympd(A, V))
    {
	fmat b;
	fmat r = Y;
	if(X.n_cols > 0)
	{
	    b = pinv(X.t() * X) * X.t() * Y;
	    r = Y - X * b;
	}
        w = w * 0;
        w[0] = as_scalar(accu(r % r) / (N - X.n_cols));
	ret["vcs"] = w;
	ret["fix"] = b;
        ret["err"] = 1;
	return(Rcpp::wrap(ret));
    }
    // Rcpp::Rcout << "A = " << A << std::endl;

    /* Get P, and Q = I - P, then R V_i and R Y, where R = V^Q */
    std::vector<fmat> RV(M);
    fmat Ry;
    fmat B, C, P;
    if(X.n_cols == 0)
    {
	for(int i = 0; i < M; i++)
	    RV[i] = A * K[i];
	Ry = A * Y;
    }
    else
    {
	// B = V^X = solve(V, X)
	// P = X (X'V^X)^ (V^X)' = X (X'B)^ B'
	B = A * X;
	C = pinv(X.t() * B) * B.t();
	P = X * C;
	for(int i = 0; i < M; i++)
	    RV[i] = A * (K[i] - P * K[i]);
	Ry = A * (Y - P * Y);
    }
    // Rcpp::Rcout << "Ry=" << Ry << std::endl;
    
    // u_i = e' V_i e = y' R V_i R y
    fmat u(M, 1);
    for(int i = 0; i < M; i++)
	u[i] = as_scalar((Ry.t() * K[i] * Ry));
    // Rcpp::Rcout << "u=" << u << std::endl;
	    
    // F_{i,j} = Tr(R V_i R V_j)
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

    /* Solve variance components */
    ret["vcs"] = pinv(F) * u;

    int gls = as<int>(opt["gls"]);
    // Rcpp::Rcout << "GLS=" << gls << std::endl;
    if(X.n_cols > 0 && gls > 0)
    {
        if(gls > 1)
        {
            V = K[0] * w[0];
            for(int i = 1; i < M; i++)
                V += K[i] * w[i];
            B = solve(V, X, solve_opts::fast + solve_opts::likely_sympd);
            ret["fix"] = pinv(X.t() * B) * B.t() * Y;
        }
        else
            ret["fix"] = C * Y;
    }
    ret["err"] = 0;
    return(Rcpp::wrap(ret));
}

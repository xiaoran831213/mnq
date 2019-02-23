#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// MINQUE
//
// y vector of response variable
// V list of kernels
// w vector of initial weights (variance components).
// X matrix of covariate
SEXP amo_mnq(SEXP _y, SEXP _V, SEXP _X, SEXP _w)
{
    // // [[Rcpp::export]]
    fmat y = as<fmat>(_y);	// column vector of response
    int N = y.n_rows;		// sample size
    List ret;			// list of results

    /* list of kernels */
    List _l = as<List>(_V);	// list of kernels
    int M = _l.length();	// number of kernels
    std::vector<fmat> v(M);
    for(int i = 0; i < M; i++)
	v[i] = as<fmat>(_l[i]);

    /* initial weights */
    fmat w = as<fmat>(_w);
    if(w.n_cols == 0)
	w = fmat(M, 1, fill::ones) / M;
    // Rcpp::Rcout << "w=" << w << std::endl;

    /* covariate */
    fmat X = as<fmat>(_X);
    // Rcpp::Rcout << "X=" << X << std::endl;

    // V: weighted sum of kernels, and A=V^ (V inverse)
    fmat V(N, N, fill::zeros);
    for(int i = 0; i < M; i++)
	V += v[i] * w[i];
    // Rcpp::Rcout << "V = " << V << std::endl;
    fmat A;
    bool psd=true;
    try
    {
	A = inv_sympd(V);
    }
    catch(std::runtime_error ex)
    {
	psd = false;
    }
    if(!psd) // fall back to pinv
    {
	// Rcpp::Rcout << "fallback to pinv" << std::endl;
	A = pinv(V);
    }
    ret["psd"] = psd;
    // Rcpp::Rcout << "A = " << A << std::endl;

    /* Get P, and Q = I - P, then R V_i and R y, where R = V^Q */
    std::vector<fmat> RV(M);
    fmat Ry;
    fmat B, C, P;
    if(X.n_cols == 0)
    {
	for(int i = 0; i < M; i++)
	    RV[i] = A * v[i];
	Ry = A * y;
    }
    else
    {
	// B = V^X = solve(V, X)
	// P = X (X'V^X)^ (V^X)' = X (X'B)^ B'
	B = A * X;
	C = pinv(X.t() * B) * B.t();
	P = X * C;
	for(int i = 0; i < M; i++)
	    RV[i] = A * (v[i] - P * v[i]);
	Ry = A * (y - P * y);
    }
    // Rcpp::Rcout << "Ry=" << Ry << std::endl;
    
    // u_i = e' V_i e = y' R V_i R y
    fmat u(M, 1);
    for(int i = 0; i < M; i++)
	u[i] = as_scalar((Ry.t() * v[i] * Ry));
    // Rcpp::Rcout << "u=" << u << std::endl;
	    
    // F_{i,j} = Tr(R V_i R V_j)
    fmat F(M, M);
    for(int i = 0; i < M; i++)
    {
	for(int j = i + 1; j < M; j++)
	{
	    F(i, j) = accu(RV[i] % RV[j]);
	    F(j, i) = F(i, j);
	}
	F(i, i) = accu(RV[i] % RV[i]);
    }
    // Rcpp::Rcout << "F=" << F << std::endl;

    /* Solve variance components */
    // [s2_1, s2_2, ..., s2_K]^T =
    // [yA1y, yA2y, ..., yAMy]^T = solve(F, u) = pinv(F) u
    ret["vcs"] = pinv(F) * u;
    // Rcpp::Rcout << "w=" << w << std::endl;
    // [A_1, A_2, ..., A_M] is not directly calculated

    /* GLS for fixed effects */
    // beta(s) = (X'V^X)^ (V^X)' y
    if(X.n_cols > 0)
	ret["fix"] = C * y;

    return(Rcpp::wrap(ret));
}


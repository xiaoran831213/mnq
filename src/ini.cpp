#include <Rcpp.h>

using namespace Rcpp;

SEXP amo_mnq(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP amo_mn0(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] =
{
    {"amo_mnq", (DL_FUNC) &amo_mnq, 5},
    {"amo_mn0", (DL_FUNC) &amo_mn0, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mnq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

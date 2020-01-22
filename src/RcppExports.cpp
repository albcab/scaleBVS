// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// wTGS
List wTGS(SEXP X_, SEXP y_, SEXP n_, SEXP p_, SEXP n_iter, SEXP burnin_, SEXP h1_, SEXP h2_, SEXP c_, SEXP k_, SEXP weighted);
RcppExport SEXP _scaleBVS_wTGS(SEXP X_SEXP, SEXP y_SEXP, SEXP n_SEXP, SEXP p_SEXP, SEXP n_iterSEXP, SEXP burnin_SEXP, SEXP h1_SEXP, SEXP h2_SEXP, SEXP c_SEXP, SEXP k_SEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type burnin_(burnin_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type h1_(h1_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type h2_(h2_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type c_(c_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(wTGS(X_, y_, n_, p_, n_iter, burnin_, h1_, h2_, c_, k_, weighted));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scaleBVS_wTGS", (DL_FUNC) &_scaleBVS_wTGS, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_scaleBVS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

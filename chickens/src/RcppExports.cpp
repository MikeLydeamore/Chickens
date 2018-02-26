// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// chickens_model
List chickens_model(List parameters_patch, NumericMatrix betas, double max_time, double dt, int solver_type);
RcppExport SEXP _chickens_chickens_model(SEXP parameters_patchSEXP, SEXP betasSEXP, SEXP max_timeSEXP, SEXP dtSEXP, SEXP solver_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type parameters_patch(parameters_patchSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type solver_type(solver_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(chickens_model(parameters_patch, betas, max_time, dt, solver_type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chickens_chickens_model", (DL_FUNC) &_chickens_chickens_model, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_chickens(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

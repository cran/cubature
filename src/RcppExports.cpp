// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/cubature.h"
#include <Rcpp.h>

using namespace Rcpp;

// doHCubature
Rcpp::List doHCubature(int fDim, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL, int maxEval, double absErr, double tol, int vectorInterface, unsigned norm);
RcppExport SEXP _cubature_doHCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type fDim(fDimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type f(fSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xLL(xLLSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xUL(xULSEXP);
    Rcpp::traits::input_parameter< int >::type maxEval(maxEvalSEXP);
    Rcpp::traits::input_parameter< double >::type absErr(absErrSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type vectorInterface(vectorInterfaceSEXP);
    Rcpp::traits::input_parameter< unsigned >::type norm(normSEXP);
    rcpp_result_gen = Rcpp::wrap(doHCubature(fDim, f, xLL, xUL, maxEval, absErr, tol, vectorInterface, norm));
    return rcpp_result_gen;
END_RCPP
}
// doPCubature
Rcpp::List doPCubature(int fDim, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL, int maxEval, double absErr, double tol, int vectorInterface, unsigned norm);
RcppExport SEXP _cubature_doPCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type fDim(fDimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type f(fSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xLL(xLLSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xUL(xULSEXP);
    Rcpp::traits::input_parameter< int >::type maxEval(maxEvalSEXP);
    Rcpp::traits::input_parameter< double >::type absErr(absErrSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type vectorInterface(vectorInterfaceSEXP);
    Rcpp::traits::input_parameter< unsigned >::type norm(normSEXP);
    rcpp_result_gen = Rcpp::wrap(doPCubature(fDim, f, xLL, xUL, maxEval, absErr, tol, vectorInterface, norm));
    return rcpp_result_gen;
END_RCPP
}
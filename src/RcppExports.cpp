// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Counts
SEXP Counts(arma::colvec anno, arma::mat geno, const int n_anno, const int min_mac);
RcppExport SEXP _AllelicSeries_Counts(SEXP annoSEXP, SEXP genoSEXP, SEXP n_annoSEXP, SEXP min_macSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type anno(annoSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const int >::type n_anno(n_annoSEXP);
    Rcpp::traits::input_parameter< const int >::type min_mac(min_macSEXP);
    rcpp_result_gen = Rcpp::wrap(Counts(anno, geno, n_anno, min_mac));
    return rcpp_result_gen;
END_RCPP
}
// OLS
SEXP OLS(const arma::colvec y, const arma::mat X);
RcppExport SEXP _AllelicSeries_OLS(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(OLS(y, X));
    return rcpp_result_gen;
END_RCPP
}
// ResidVar
SEXP ResidVar(const arma::colvec y, const arma::mat X);
RcppExport SEXP _AllelicSeries_ResidVar(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ResidVar(y, X));
    return rcpp_result_gen;
END_RCPP
}
// Score
SEXP Score(const arma::colvec y, const arma::mat G, const arma::mat X, const double v);
RcppExport SEXP _AllelicSeries_Score(SEXP ySEXP, SEXP GSEXP, SEXP XSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Score(y, G, X, v));
    return rcpp_result_gen;
END_RCPP
}
// CalcSkatWeights
SEXP CalcSkatWeights(const arma::vec& anno, const arma::vec& maf, const arma::vec& weights);
RcppExport SEXP _AllelicSeries_CalcSkatWeights(SEXP annoSEXP, SEXP mafSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type anno(annoSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcSkatWeights(anno, maf, weights));
    return rcpp_result_gen;
END_RCPP
}
// GetLambda
SEXP GetLambda(const arma::mat& kernel);
RcppExport SEXP _AllelicSeries_GetLambda(SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLambda(kernel));
    return rcpp_result_gen;
END_RCPP
}
// SkatOptimalParam
SEXP SkatOptimalParam(const arma::mat& cov_z, const arma::vec& rhos);
RcppExport SEXP _AllelicSeries_SkatOptimalParam(SEXP cov_zSEXP, SEXP rhosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cov_z(cov_zSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rhos(rhosSEXP);
    rcpp_result_gen = Rcpp::wrap(SkatOptimalParam(cov_z, rhos));
    return rcpp_result_gen;
END_RCPP
}
// CorCpp
SEXP CorCpp(const arma::mat& x);
RcppExport SEXP _AllelicSeries_CorCpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CorCpp(x));
    return rcpp_result_gen;
END_RCPP
}
// isPD
SEXP isPD(const arma::mat& x, bool force_symmetry, double tau);
RcppExport SEXP _AllelicSeries_isPD(SEXP xSEXP, SEXP force_symmetrySEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type force_symmetry(force_symmetrySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(isPD(x, force_symmetry, tau));
    return rcpp_result_gen;
END_RCPP
}
// AnnoMat
arma::mat AnnoMat(const arma::colvec& anno, const int n_anno);
RcppExport SEXP _AllelicSeries_AnnoMat(SEXP annoSEXP, SEXP n_annoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type anno(annoSEXP);
    Rcpp::traits::input_parameter< const int >::type n_anno(n_annoSEXP);
    rcpp_result_gen = Rcpp::wrap(AnnoMat(anno, n_anno));
    return rcpp_result_gen;
END_RCPP
}
// BaselineSS
SEXP BaselineSS(const arma::colvec& anno, const arma::colvec& beta, const arma::mat& ld, const arma::colvec& se, const int n_anno);
RcppExport SEXP _AllelicSeries_BaselineSS(SEXP annoSEXP, SEXP betaSEXP, SEXP ldSEXP, SEXP seSEXP, SEXP n_annoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type anno(annoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ld(ldSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type se(seSEXP);
    Rcpp::traits::input_parameter< const int >::type n_anno(n_annoSEXP);
    rcpp_result_gen = Rcpp::wrap(BaselineSS(anno, beta, ld, se, n_anno));
    return rcpp_result_gen;
END_RCPP
}
// SumCountSS
SEXP SumCountSS(const arma::colvec& anno, const arma::colvec& beta, const arma::mat& ld, const arma::colvec& se, const arma::colvec& weights);
RcppExport SEXP _AllelicSeries_SumCountSS(SEXP annoSEXP, SEXP betaSEXP, SEXP ldSEXP, SEXP seSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type anno(annoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ld(ldSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(SumCountSS(anno, beta, ld, se, weights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_AllelicSeries_Counts", (DL_FUNC) &_AllelicSeries_Counts, 4},
    {"_AllelicSeries_OLS", (DL_FUNC) &_AllelicSeries_OLS, 2},
    {"_AllelicSeries_ResidVar", (DL_FUNC) &_AllelicSeries_ResidVar, 2},
    {"_AllelicSeries_Score", (DL_FUNC) &_AllelicSeries_Score, 4},
    {"_AllelicSeries_CalcSkatWeights", (DL_FUNC) &_AllelicSeries_CalcSkatWeights, 3},
    {"_AllelicSeries_GetLambda", (DL_FUNC) &_AllelicSeries_GetLambda, 1},
    {"_AllelicSeries_SkatOptimalParam", (DL_FUNC) &_AllelicSeries_SkatOptimalParam, 2},
    {"_AllelicSeries_CorCpp", (DL_FUNC) &_AllelicSeries_CorCpp, 1},
    {"_AllelicSeries_isPD", (DL_FUNC) &_AllelicSeries_isPD, 3},
    {"_AllelicSeries_AnnoMat", (DL_FUNC) &_AllelicSeries_AnnoMat, 2},
    {"_AllelicSeries_BaselineSS", (DL_FUNC) &_AllelicSeries_BaselineSS, 5},
    {"_AllelicSeries_SumCountSS", (DL_FUNC) &_AllelicSeries_SumCountSS, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_AllelicSeries(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

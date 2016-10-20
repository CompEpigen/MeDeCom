// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cppTAfact
List cppTAfact(const RMatrixIn& mDt, const RMatrixIn& mTtinit, const RMatrixIn& mAinit, double lambda, int itersMax, double tol, double tolA, double tolT);
RcppExport SEXP MeDeCom_cppTAfact(SEXP mDtSEXP, SEXP mTtinitSEXP, SEXP mAinitSEXP, SEXP lambdaSEXP, SEXP itersMaxSEXP, SEXP tolSEXP, SEXP tolASEXP, SEXP tolTSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const RMatrixIn& >::type mDt(mDtSEXP);
    Rcpp::traits::input_parameter< const RMatrixIn& >::type mTtinit(mTtinitSEXP);
    Rcpp::traits::input_parameter< const RMatrixIn& >::type mAinit(mAinitSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type itersMax(itersMaxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type tolA(tolASEXP);
    Rcpp::traits::input_parameter< double >::type tolT(tolTSEXP);
    __result = Rcpp::wrap(cppTAfact(mDt, mTtinit, mAinit, lambda, itersMax, tol, tolA, tolT));
    return __result;
END_RCPP
}
// RHLasso
List RHLasso(NumericMatrix Ginp, NumericMatrix Winp, NumericMatrix Ainp, NumericVector l);
RcppExport SEXP MeDeCom_RHLasso(SEXP GinpSEXP, SEXP WinpSEXP, SEXP AinpSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ginp(GinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Winp(WinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ainp(AinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    __result = Rcpp::wrap(RHLasso(Ginp, Winp, Ainp, l));
    return __result;
END_RCPP
}
// RQuadHC
List RQuadHC(NumericMatrix Ginp, NumericMatrix Winp, NumericMatrix Ainp, NumericVector otol, NumericVector lconstr, NumericVector uconstr);
RcppExport SEXP MeDeCom_RQuadHC(SEXP GinpSEXP, SEXP WinpSEXP, SEXP AinpSEXP, SEXP otolSEXP, SEXP lconstrSEXP, SEXP uconstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ginp(GinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Winp(WinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ainp(AinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type otol(otolSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lconstr(lconstrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type uconstr(uconstrSEXP);
    __result = Rcpp::wrap(RQuadHC(Ginp, Winp, Ainp, otol, lconstr, uconstr));
    return __result;
END_RCPP
}
// RProjSplxBox
NumericMatrix RProjSplxBox(NumericMatrix Xinp, NumericVector linp, NumericVector uinp);
RcppExport SEXP MeDeCom_RProjSplxBox(SEXP XinpSEXP, SEXP linpSEXP, SEXP uinpSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Xinp(XinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type linp(linpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type uinp(uinpSEXP);
    __result = Rcpp::wrap(RProjSplxBox(Xinp, linp, uinp));
    return __result;
END_RCPP
}
// RQuadSimplex
List RQuadSimplex(NumericMatrix Ginp, NumericMatrix Winp, NumericMatrix Ainp, NumericVector ot);
RcppExport SEXP MeDeCom_RQuadSimplex(SEXP GinpSEXP, SEXP WinpSEXP, SEXP AinpSEXP, SEXP otSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ginp(GinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Winp(WinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ainp(AinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ot(otSEXP);
    __result = Rcpp::wrap(RQuadSimplex(Ginp, Winp, Ainp, ot));
    return __result;
END_RCPP
}
// RQuadSimplexBox
List RQuadSimplexBox(NumericMatrix Ginp, NumericMatrix Winp, NumericMatrix Ainp, NumericVector linp, NumericVector uinp, NumericVector ot);
RcppExport SEXP MeDeCom_RQuadSimplexBox(SEXP GinpSEXP, SEXP WinpSEXP, SEXP AinpSEXP, SEXP linpSEXP, SEXP uinpSEXP, SEXP otSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ginp(GinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Winp(WinpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ainp(AinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type linp(linpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type uinp(uinpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ot(otSEXP);
    __result = Rcpp::wrap(RQuadSimplexBox(Ginp, Winp, Ainp, linp, uinp, ot));
    return __result;
END_RCPP
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_stan_fit4hmm_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4hmm_deltas_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mm_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mm10d01_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mm10d9_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4hmm_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4hmm_mod, 0},
    {"_rcpp_module_boot_stan_fit4hmm_deltas_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4hmm_deltas_mod, 0},
    {"_rcpp_module_boot_stan_fit4mm_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mm_mod, 0},
    {"_rcpp_module_boot_stan_fit4mm10d01_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mm10d01_mod, 0},
    {"_rcpp_module_boot_stan_fit4mm10d9_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mm10d9_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_sampinfo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

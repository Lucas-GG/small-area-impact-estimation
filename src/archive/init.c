#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _rcpp_module_boot_combined();
extern SEXP _rcpp_module_boot_ST_evaluation_case();

extern SEXP _ST_evaluation_case_shuffle_cl_cpp(SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_cluster_boot_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_shuffle_start_cpp(SEXP, SEXP);
extern SEXP _ST_evaluation_case_gmean_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_reset_y0_cpp(SEXP);
extern SEXP _ST_evaluation_case_expandx_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_rselect_cpp(SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_fit_single_tree_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ST_evaluation_case_forest_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_combined", (DL_FUNC) &_rcpp_module_boot_combined, 0},
    {"_rcpp_module_boot_ST_evaluation_case", (DL_FUNC) &_rcpp_module_boot_ST_evaluation_case, 0},
    {"_ST_evaluation_case_shuffle_cl_cpp", (DL_FUNC) &_ST_evaluation_case_shuffle_cl_cpp, 3},
    {"_ST_evaluation_case_cluster_boot_cpp", (DL_FUNC) &_ST_evaluation_case_cluster_boot_cpp, 4},
    {"_ST_evaluation_case_shuffle_start_cpp", (DL_FUNC) &_ST_evaluation_case_shuffle_start_cpp, 2},
    {"_ST_evaluation_case_gmean_cpp", (DL_FUNC) &_ST_evaluation_case_gmean_cpp, 5},
    {"_ST_evaluation_case_reset_y0_cpp", (DL_FUNC) &_ST_evaluation_case_reset_y0_cpp, 1},
    {"_ST_evaluation_case_expandx_cpp", (DL_FUNC) &_ST_evaluation_case_expandx_cpp, 4},
    {"_ST_evaluation_case_rselect_cpp", (DL_FUNC) &_ST_evaluation_case_rselect_cpp, 3},
    {"_ST_evaluation_case_fit_single_tree_cpp", (DL_FUNC) &_ST_evaluation_case_fit_single_tree_cpp, 4},
    {"_ST_evaluation_case_forest_rcpp", (DL_FUNC) &_ST_evaluation_case_forest_rcpp, 10},
    {NULL, NULL, 0}
};

void R_init_ST_evaluation_case(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

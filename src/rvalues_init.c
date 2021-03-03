#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP Vcut(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ZeroIn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Vcut",   (DL_FUNC) &Vcut,   5},
  {"ZeroIn", (DL_FUNC) &ZeroIn, 9},
  {NULL, NULL, 0}
};

void R_init_rvalues(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

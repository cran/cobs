#include <stdlib.h> // for NULL

#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#include "splines.h"
#include "SR.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

// ./splines.c :
static R_NativePrimitiveArgType spline_basis_typ[] = {
    REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP };

static R_NativePrimitiveArgType spline_value_typ[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP, REALSXP };

static R_NativePrimitiveArgType SR_R_typ[] = {
    INTSXP, REALSXP, INTSXP, INTSXP,
    REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP };

static const R_CMethodDef CEntries[] = {
    CDEF(spline_basis),
    CDEF(spline_value),
    CDEF(SR_R),
    {NULL, NULL, 0}
};

void R_init_cobs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

/*  Routines for manipulating B-splines.  These are intended for use with
 *  the New S language described in Becker, Chambers, and Wilks (Wadsworth,
 *  1988).
 *  Copyright 1989 by Douglas M. Bates.  Permission to use, copy,
 *  modify and distribute is freely granted to all.
 *  The routines are loosely based on the pseudo-code in Schumaker (Wiley,
 *  1981) and the CMLIB library DBSPLINES.
*/

/* ============
 * REPLACE THIS by calling much nicer code in standard R "splines" package !!
 * ============
 */

#include <R.h>
/* defines `Sint' as `int' for R; (and define it as `long' for S-plus !!)
 *          ====      ~~~                            ~~~~      !!!!!!
 * and  S_alloc() */
#ifdef not_in_R
# ifdef __STDC__
  extern void *S_alloc();
# else
  extern char *S_alloc();
# endif
#endif

/* GLOBALS VARIABLES for difference tables : */
static double *ldel, *rdel;
static Sint orderm1;

static void diff_table(double *ti, double x, register int ndiff)
{
  register double *r = rdel, *l = ldel, *dpt = ti;

  while (ndiff--) {
    *r++ = *dpt++ - x;
    *l++ = x - *--ti;
  }
}


static void basis_funcs(double *ti, double x, double *b)
{
  int j, r;
  double saved, term;

  diff_table(ti, x, orderm1);
  b[0] = 1.;
  for (j = 1; j <= orderm1; j++) {
    saved = 0.;
    for (r = 0; r < j; r++) {
      term = b[r]/(rdel[r] + ldel[j - 1 - r]);
      b[r] = saved + rdel[r] * term;
      saved = ldel[j - 1 - r] * term;
    }
    b[j] = saved;
  }
}

static double evaluate(double *ti, double x, double *a, int nder)
{
  register double *lpt, *rpt, *apt;
  register int inner;
  int outer = orderm1;

  while(nder--) {
    for(inner = outer, apt = a, lpt = ti - outer; inner--; apt++, lpt++)
      *apt = outer * (*(apt + 1) - *apt)/(*(lpt + outer) - *lpt);
    outer--;
  }
  diff_table(ti, x, outer);
  while(outer--)
    for(apt = a, lpt = ldel + outer, rpt = rdel, inner = outer + 1;
	inner--; lpt--, rpt++, apt++)
      *apt = (*(apt + 1) * *lpt + *apt * *rpt)/(*rpt + *lpt);
  return(*a);
}

void spline_value(double *knots, double *coeff, Sint *ncoeff, Sint *order,
		  double *x, Sint *nx, Sint *deriv, double *y)
{
  Sint n = *nx;
  double *a, *last = knots + *ncoeff;

  a = (double *) S_alloc(*order, sizeof(double));

  /* allocate the GLOBAL difference tables : */
  orderm1 = *order - 1;
  rdel = (double *) S_alloc(orderm1, sizeof(double));
  ldel = (double *) S_alloc(orderm1, sizeof(double));

  knots += *order;		/* First *order knots must be < all x's */
  while(n--) {
    while(knots < last && *knots <= *x) {knots++; coeff++;}
    memcpy((char *) a, (char *) coeff, (int) *order * sizeof(double));
    *y++ = evaluate(knots, *x++, a, (int) *deriv);
  }
}

void spline_basis(double *knots, Sint *ncoeff, Sint *order,
		  double *xvals, Sint *derivs, Sint *nx,
		  double *basis,
		  Sint *offsets)/* offsets : pointers for bin number */
{
/* evaluate the non-zero B-spline basis
 * functions (or their derivatives) at xvals.
 */
  int n = *nx, i, j;
  double *dpt, *coeff, *last = knots + *ncoeff;

  dpt = (knots += *order);	/* first order knots must be <= all xvals */

  /* allocate the GLOBAL difference tables : */
  orderm1 = *order - 1;
  rdel = (double *) S_alloc(orderm1, sizeof(double));
  ldel = (double *) S_alloc(orderm1, sizeof(double));

  coeff = (double *) S_alloc(*order, sizeof(double));
  for( ; n--; xvals++, derivs++) {
    while(dpt < last && *dpt <= *xvals) dpt++;
    if (*derivs) {		/* slow method for derivatives */
      for(i = 0; i < *order; i++) {
	for(j = 0; j < *order; j++) coeff[j] = 0;
	coeff[i] = 1;
	*basis++ = evaluate(dpt, *xvals, coeff, (int) *derivs);
      }
    }
    else {
      basis_funcs(dpt, *xvals, basis); /* fast method for value */
      basis += *order;
    }
    *offsets++ = (Sint)(dpt - knots);
  }
}

#ifdef never_used_in_COBS
void lin_interp(double *x, double *y, double *x0, double *y0, Sint *nvals)
{
  Sint n = *nvals;

  while(n--) {
    while (*x < *x0) {x++; y++;}
    if (*x > *x0) *y0++ = *y + (*(y+1) - *y)*(*x0 - *x)/(*(x+1) - *x);
    else *y0++ = *y;
    x0++;
  }
}
#endif

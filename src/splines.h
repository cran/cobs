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

void spline_basis(double *knots, Sint *ncoeff, Sint *order,
		  double *xvals, Sint *derivs, Sint *nx,
		  double *basis, Sint *offsets);

void spline_value(double *knots, double *coeff, Sint *ncoeff, Sint *order,
		  double *x, Sint *nx, Sint *deriv, double *y);



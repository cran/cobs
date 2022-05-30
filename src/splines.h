#include <R.h>
typedef int Sint;

void spline_basis(double *knots, Sint *ncoeff, Sint *order,
		  double *xvals, Sint *derivs, Sint *nx,
		  double *basis, Sint *offsets);

void spline_value(double *knots, double *coeff, Sint *ncoeff, Sint *order,
		  double *x, Sint *nx, Sint *deriv, double *y);



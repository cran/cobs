#include "nrutil.h"
/* --> dvector():       using malloc()   and
 *     free_dvector() : using free() ... FIXME */

void SR(int n, double cc, int *m1, int ind[], double x[], double y[], double r[],
	double R[], double H[], double S[], double Y[], double D[],
	double tol_[], int it_max, int verbose, // <- new
	double *phi_BL, int *NumIt);

void SR_R(int* n, double* cc, int *m1, int ind[],
	  double x[], double y[],
	  double r[], double R[], double H[],
	  double S[], double Y[], double D[],
	  double tol[], int *it_max, int *verbose,
	  double *phi_BL, int *NumIt);

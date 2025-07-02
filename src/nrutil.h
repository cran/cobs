/* nrutil.h. From Numerical Recipes. Floats are changed to doubles. */
/* --------- MM: modified to eliminate compiler warnings */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/* MM: quite a few IMIN(), DMAX(), .. macros removed which used static variables
   --  and were nowhere called, hence the ==> -Wunused  warnings (in gcc, macOS, 2025)
*/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#ifdef _not_for_R_
void nrerror(char error_text[]);
#else
# include <R_ext/Error.h>
# define nrerror Rf_error
#endif
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
short *shvector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
char *charvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_shvector(short *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_charvector(char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
double ran1(long *idum);
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);

#endif /* _NR_UTILS_H_ */

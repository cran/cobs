// SR: Support Reduction algorithm by Piet Groeneboom -- ~= a C version of conreg()

#include <R.h>

#include "SR.h"
/* --> dvector():       using malloc()   and
 *     free_dvector() : using free()
 * FIXME: in R: preferably using R_alloc() ! */

void	argmin(int m, int ind[], double x[], double *min, int *indmin);
void	rindtor(int m, double x[], int ind[], double rind[], double r[]);
void    cumsum(int n, const double x[], /* --> */ double Sx[]);
void 	s3D(int n, double **a, double b[], double r[]);
void 	detr(int m, double x[], double y[], int ind[], double rind[]);
int CheckPositivity(double x[], double y[], double tol, int m,
		    int ind[], double delta[], double rind[]);

void 	int_sort(int m, int ind[]);
void	indexx_int(int n, int arr[], int indx[]);

// SR : Support Reduction ------
void SR(int n, double cc, int *m1, int ind[], double x[], double y[], double r[],
	double R[], double H[], double S[], double Y[], double D[],
	double tol_[], int it_max, int verbose, // <- new
	double *phi_BL, int *NumIt)
{
    short cont;
    int i,m,indmin,nIt;
    double phi,min,sum,
	tol     = tol_[0],
	tol_pos = tol_[1],
	*delta = dvector(1,n),
	*rind = dvector(1,n);


    // all these arrays are 1-indexed in code below
    // originally: created in main.c with nrutil's  dvector(1,n)
    x--; y--; r--;
    R--; H--; S--; Y--; D--;

    m = 2;

    ind[1] = 1;
    ind[2] = n;

    detr(m,x,y,ind, /* -> */ rind);
    rindtor(m,x,ind,rind, /* -> */ r);

    S[1] = cc*y[1]/n;
    for (i = 2; i <= n; i++)
	S[i]= S[i-1] + cc*y[i]/n;

    Y[1] = 0;
    for (i = 2; i <= n; i++)
	Y[i]= Y[i-1] + S[i-1]*(x[i]-x[i-1]);

    delta[1] = x[n]-x[1];

    if(verbose >= 2) Rprintf("nIt   m  m'     min     iMin\n");
			 //   123 123 123 123456789012 123
			 //    3d 3d   3d     12g       3d
    indmin = 1;
    nIt = 0;
    cont = 1;
    while (cont && nIt < it_max)
    {
	nIt++;
	if(verbose >= 2) Rprintf("%3d %3d ", nIt, m);

	if (m > 2)
	{
	    m = CheckPositivity(x, y, tol_pos, m,
				// input *and* output:
				ind, delta, rind);
	    rindtor(m,x,ind,rind, /* -> */ r);
	}
	if(verbose >= 2) Rprintf("%3d ", m);


        /* Now a feasible function is found.
	   Next step is to find a profitable direction
	*/

	R[1] = cc*r[1]/n;
	for (i = 2; i <= n; i++)
	    R[i]= R[i-1] + cc*r[i]/n;

	H[1] = 0;
	for (i = 2; i <= n; i++)
	    H[i]= H[i-1] + R[i-1]*(x[i]-x[i-1]);

	for (i = 1; i <= n; i++)
	    D[i] = H[i]-Y[i];

	argmin(m,ind,D,&min,&indmin);

	if(verbose >= 2) Rprintf("%12g %3d\n", min, indmin);

	if (min < -tol)
	{
	    m++;
	    ind[m] = indmin;
	    int_sort(m,ind);
	    detr(m,x,y,ind, /* -> */ rind);
	    for (i = 1; i < m; i++)
		delta[i] = x[ind[i+1]]-x[ind[i]];
	}
	else cont = 0;
    }

    for (i = 1,sum = 0; i <= n; i++)
	sum += SQR(r[i]-y[i]);
    phi = 0.5*sum;

    if(verbose)
	Rprintf("c(nIt=%4d, phi=%12.6g, min=%14.8g, iMin=%6d, m=%4d)\n",
		nIt,phi,min,indmin, m);

    *m1 = m;
    *NumIt = nIt;

    *phi_BL = phi;

    free_dvector(delta,1,n);
    free_dvector(rind,1,n);
}

// For R interface
void SR_R(int* n, double* cc, int *m1, int ind[],
	  double x[], double y[],
	  double r[], double R[], double H[],
	  double S[], double Y[], double D[],
	  double tol[], int* it_max, int* verbose,
	  double *phi_BL, int *NumIt)
{
    SR(n[0], cc[0], m1, ind, x,y, r, R,H, S,Y,D, tol,
       it_max[0], verbose[0],
       // result :
       phi_BL, NumIt);
}

void argmin(int m, int ind[], double x[], double *min, int *indmin)
{
    double min1 = 0.;
    int indmin1 = 1;

    for (int i = 1; i < m; i++)
    {
	for (int k = ind[i]+1; k < ind[i+1]; k++)
	{
	    if (x[k] < min1)
	    {
		min1 = x[k];
		indmin1 = k;
	    }
	}
    }

    *min = min1;
    *indmin = indmin1;
}

/* {comment from Geurt's R code:}
 Needed function: s3D
 The function s3D solves a tridiagonal system  A r = b,
 i.e.,   r :=  A^{-1} b
 Here a[,] has three rows, where the first row
 consists of a 0 and the superdiagonal elements, the
 second row contains the diagonal elements and the
 third row the subdiagonal entries (last element 0)

 a[,] is modified;  b[] is not
*/
void s3D(int m, double **a, double b[], double r[])
{
    int i;
    double *x = dvector(1,m);

    for (i = 1; i <= m; i++)
	x[i] = b[i];

    for (i = 2; i <= m; i++)
    {
	double xmult = a[3][i-1]/a[2][i-1];
	a[2][i] -= xmult*a[1][i];
	x[i] 	-= xmult*x[i-1];
    }

    x[m] /=a[2][m];

    for (i = m-1; i >= 1; i--)
	x[i] = (x[i]-a[1][i+1]*x[i+1])/a[2][i];

    for (i = 1; i <= m; i++)
	r[i] = x[i];

    free_dvector(x,1,m);
}


/* {comment from Geurt's R code:}
 function rindtor, which extends the definition
 of a function r on the x determined by ind
 to all x
*/
void rindtor(int m, double x[], int ind[], double rind[],
	     /* --> */ double r[])
{
    for (int i = 1; i <= m; i++)
	r[ind[i]] = rind[i];

    for (int i = 1; i <= m-1; i++)
    {
	for (int j = ind[i]+1; j < ind[i+1]; j++)
	{
	    double
		a= (x[ind[i+1]]-x[j]) / (x[ind[i+1]] - x[ind[i]]),
		b= (x[j] - x[ind[i]]) / (x[ind[i+1]] - x[ind[i]]);
	    r[j]= a*rind[i] + b*rind[i+1];
	}
    }
}


void cumsum(int n, const double x[], /* --> */ double Sx[])
{
    Sx[1] = x[1];
    for (int i = 2; i <= n; i++)
	Sx[i] = Sx[i-1]+x[i];
}

/* {comment from Geurt's R code:}
 Needed function: detr
 Performance can be optimized a bit further, because some
 objects are recomputed. Given the x, y vector and the vector
 of kink-indices, this function computes the unrestricted
 least squares linear spline
*/
void detr(int m,
	  double x[], double y[],
	  int ind[],
	  // result:
	  double rind[])
{
    int i,j;
    double sum1, sum2,
	*rl  = dvector(1,m),
	*rl1 = dvector(1,m),
	*rl2 = dvector(1,m),
	*D1  = dvector(1,m),
	*D2  = dvector(1,m),
	**a  = dmatrix(1,3,1,m);

    for (i = 1; i <= 3; i++)
	for (j = 1; j <= m; j++)
	    a[i][j] = 0;

    for (i = 1; i <= m; i++)
	rl[i] = rl1[i] = rl2[i] = 0;

    for (i = 1; i <= m-1; i++) {
	D1[i] = x[ind[i+1]]-x[ind[i]];
	D2[i] = SQR(D1[i]);
    }

    for (j = ind[1]; j < ind[2]; j++)
	a[2][1]+= SQR(x[ind[2]]-x[j])/D2[1];

    for (j = ind[1]+1; j < ind[2]; j++)
	a[3][1]+= (x[ind[2]]-x[j])*(x[j]-x[ind[1]])/D2[1];

    for (j = ind[1],rl[1] = 0; j < ind[2]; j++)
	rl[1] += (x[ind[2]]-x[j])*y[j]/D1[1];

    for (i = 2; i <= m-1; i++)
    {
	for (j = ind[i], sum1 = 0; j < ind[i+1]; j++)
	    sum1 += SQR(x[ind[i+1]]-x[j])/D2[i];
	for (j = ind[i-1]+1, sum2 = 0; j <= ind[i]; j++)
	    sum2 += SQR(x[j]-x[ind[i-1]])/D2[i-1];

	a[2][i] = sum1+sum2-1;

	for (j = ind[i],a[3][i] = 0; j < ind[i+1]; j++)
	    a[3][i] += (x[ind[i+1]]-x[j])*(x[j]-x[ind[i]])/D2[i];
	for (j = ind[i-1]+1,rl1[i] = 0; j <= ind[i]; j++)
	    rl1[i] += y[j]*(x[j]-x[ind[i-1]])/D1[i-1];
	for (j = ind[i],rl2[i] = 0; j < ind[i+1]; j++)
	    rl2[i] += y[j]*(x[ind[i+1]]-x[j])/D1[i];

	rl[i] = rl1[i]+rl2[i]-y[ind[i]];
    }

    a[1][1] = 0;
    for (i = 2; i <= m; i++)
	a[1][i] = a[3][i-1];

    for (j = ind[m-1]+1,a[2][m] = 0; j <= ind[m]; j++)
	a[2][m] += SQR(x[j]-x[ind[m-1]])/D2[m-1];

    for (j = ind[m-1]+1,rl[m] = 0; j <= ind[m]; j++)
	rl[m] += y[j]*(x[j]-x[ind[m-1]])/D1[m-1];

    // solve  A r_ind = rl,  i.e.,  rind := A^{-1} rl :
    s3D(m,a,rl,rind);

    free_dvector(rl,1,m);
    free_dvector(rl1,1,m);
    free_dvector(rl2,1,m);
    free_dvector(D1,1,m);
    free_dvector(D2,1,m);
    free_dmatrix(a,1,3,1,m);
}


int CheckPositivity(double x[], double y[], double tol, int m,
		    // input _and_ output :
		    int ind[], double delta[], double rind[])
{
    int imin = 0, i;
    double min = 0;

    for (i = 2; i < m; i++) {
	double a =
	    (rind[i+1] - rind[i]) / delta[i] -
	    (rind[i] - rind[i-1]) / delta[i-1];
	if (a < min) {
	    min = a;
	    imin = i;
	}
    }

    if (min < -tol)
    {
	while (min < 0)
	{
	    for (int k= imin; k < m; k++) {
		ind[k] = ind[k+1];
		delta[k] = delta[k+1];
		rind[k] = rind[k+1];
	    }

	    m--;

	    detr(m,x,y,ind, /* -> */ rind);
	    for (i = 1; i < m; i++)
		delta[i] = x[ind[i+1]]-x[ind[i]];

	    double min1 = 0.;

	    for (i = 2; i < m; i++) {
		double a =
		    (rind[i+1] - rind[i]) / delta[i] -
		    (rind[i] - rind[i-1]) / delta[i-1];
		if (a < min1) {
		    min1 = a;
		    imin = i;
		}
	    }
	    min = min1;
	}
    }
    return(m);
}


//----------------------- Sort Routines ------------------------------
// FIXME: Use R's instead

void 	int_sort(int m, int ind[])
{
    int j,*intarray;
    int *intarray2;

    intarray = ivector(1,m);
    intarray2 = ivector(1,m);

    indexx_int(m,ind,intarray);
    for (j = 1; j <= m; j++) intarray2[j] = ind[j];
    for (j = 1; j <= m; j++) ind[j] = intarray2[intarray[j]];

    free_ivector(intarray,1,m);
    free_ivector(intarray2,1,m);
}


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx_int(int n, int arr[], int indx[])
{
    int i,indxt,ir = n,itemp,j,k,l = 1;
    int jstack = 0,*istack;
    int a;

    istack = ivector(1,NSTACK);
    for (j = 1; j <= n; j++) indx[j] = j;
    for (;;) {
	if (ir-l < M) {
	    for (j = l+1; j <= ir; j++) {
		indxt = indx[j];
		a = arr[indxt];
		for (i = j-1; i >= 1; i--) {
		    if (arr[indx[i]] <= a) break;
		    indx[i+1] = indx[i];
		}
		indx[i+1] = indxt;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];
	} else {
	    k=(l+ir) >> 1;
	    SWAP(indx[k],indx[l+1]);
	    if (arr[indx[l+1]] > arr[indx[ir]]) {
		SWAP(indx[l+1],indx[ir])
		    }
	    if (arr[indx[l]] > arr[indx[ir]]) {
		SWAP(indx[l],indx[ir])
		    }
	    if (arr[indx[l+1]] > arr[indx[l]]) {
		SWAP(indx[l+1],indx[l])
		    }
	    i = l+1;
	    j = ir;
	    indxt = indx[l];
	    a = arr[indxt];
	    for (;;) {
		do i++; while (arr[indx[i]] < a);
		do j--; while (arr[indx[j]] > a);
		if (j < i) break;
		SWAP(indx[i],indx[j])
		    }
	    indx[l] = indx[j];
	    indx[j] = indxt;
	    jstack += 2;
	    if (jstack > NSTACK)
#ifdef _not_for_R_
		nrerror("NSTACK too small in indexx.");
#else
	        Rf_error("NSTACK too small in indexx.");
#endif
	    if (ir-i+1 >= j-l) {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    } else {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
    free_ivector(istack,1,NSTACK);
}

#undef M
#undef NSTACK
#undef SWAP

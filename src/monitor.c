/* Printing of informative messages if  (pswitch == TRUE) : */
#include <R.h>

/* called for  trace >= 1 : ----------------------------------------------- */

int F77_SUB(monit0)(int *nrq, int *nl1, int *neqc, int *niqc, int *nvars,
		    double *t, double *lam, int *itend, int *ilend, int *trace)
{
    int up = !(*itend && *ilend);
    Rprintf(" %s N[rq,L1]= (%d, %d); ([ei]qc= %d,%d, vars=%d)%s"
	    " %s%s%s (tau,lam)= (%.14g, %.14g):\n",
	    (*trace >= 2)? "dcrqL1lt():" : "",
  	    *nrq, *nl1, *neqc, *niqc, *nvars, (up || *trace >= 2)? "\n\t" : " ",
	    (*itend)?"":"TAU",
	    (*ilend)?"":"LAMBDA",
	    (up)? " will be updated;" : "_fixed_",
	    *t, *lam);
    return 0;
}

/* called for  trace >= 2 : ----------------------------------------------- */

int F77_SUB(monit11)(int *nt, int *nact,
		     double *amag, double *cgmag, int *ifl, int *trace)
{
    Rprintf("after setup() & newpen(): ifl = %d", *ifl);

    if(*ifl != 0)
	Rprintf(" *** != 0 ***\n");
    else {
	Rprintf(" (a|cg)mag = (%g, %g)\n", *amag, *cgmag);
    }
    return 0;
} /* m..11 */

int F77_SUB(monit12)(int *icyc, int *ifl, int *trace)
{
    Rprintf("%s icyc = %d iter.; ifl =%d%s",
	    (*trace >= 3)? "end of inner loop, needed" : "search used",
	    *icyc, *ifl, (*trace >= 3)? "\n" : "; ");
    return 0;
} /* m..12 */

int F77_SUB(monit13)(int *nt, int *nact, int *itend, int *ilend,
		    double *tnxt, double *lnxt, int *ifl, int *trace)
{
    Rprintf("%s %s%s= %g\n",
	    (*trace >= 3)? "after `Next Lambda & Tau': new" : "NLT:",
	    (*itend)?"":"TAU ",
	    (*ilend)?"":"LAMBDA ",
	    (*itend)? *lnxt : *tnxt);
    return 0;
} /* m..13 */

/* called for  trace >= 3 : ----------------------------------------------- */

/* monit2() : *INNER* loop */
int F77_SUB(monit2)(int *nact, int *icyc, int *trace,
		    double *x, double *alpha, double *erql1n,
		    double *pen, double *penpar, int *indx)
{
    Rprintf("  [ic= %2d] nact=%d : "
	    " pen(par) = %g, %g;  alpha, rqL1n = %g, %g\n",
	    *icyc, *nact, *pen, *penpar, *alpha, *erql1n);
    return 0;
}

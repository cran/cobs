/* just for testing out the  dunif01() generator in ./drqssbc.f */
#include <R.h>

/* IMPORT from ./drqssbc.f : */
extern double F77_NAME(dunif01)(int *, int *);

void COBS_dunif01(int *iseed, int *n, double *u)
{
    int i, idummy = 0;
    for(i = 0; i < *n; i++) {
	if(i > 0) *iseed = 0;
	u[i] = F77_CALL(dunif01)(iseed, &idummy);
	/* Rprintf("i=%d, iseed= %d ", i, iseed); */
    }
    return;
}

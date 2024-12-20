Calling tree :
=============


	cobs()
		|
		|-> shat()
		|   |
		|   \-> dn
		|
		|--> getdim2()
		|
		|-> qbsks()
		|   |
		|   |
		|   v
		\-> drqssbc2()
			|   /^
			|  /  (possibly calls itself once)
			|--
			|
			|-> l1.design2()
			|
			|-> loo.design2()
			|
			|-> rq.fit.sfn()   {quantreg}
			|
			\-> rq.fit.sfnc()  {quantreg}



qbsks(x,y,w,pw, knots,nknots, degree,Tlambda,constraint,
      equal,smaller,greater,gradient, coef,maxiter,
      pswitch, n.equal,n.smaller,n.greater,n.gradient,
      nrq,nl1, neqc, nj0, tau,lam,tmin,kmax,lstart,
      ks,mk.flag, knots.add, ic, print.mesg, factor,print.warn)

drqssbc(x,y,w,pw, knots, degree,Tlambda,constraint,
	equal,smaller,greater,gradient, coef,maxiter,
	pswitch, n.equal,n.smaller,n.greater,n.gradient,
	nrq,nl1, neqc,niqc, nvar,nj0, tau,lam,tmin,kmax,lstart,
	factor,print.warn)

-------------------------------------------------------------------------------

The following relates to the older COBS version, nowadays in package  `cobs99`
=============================================================================

COBS -- Constrained B-splines (Version 1.0)

COBS is a constrained B-splines algorithm that computes constrained
quantile curves using linear or quadratic splines.
The median spline may be used as an attractive alternative to constrained
smoothing.

A postscript version of the paper that describes the detail of COBS
can be downloaded from  http://www.cba.nau.edu/pin-ng/cobs.html

This is a *beta* version (no longer the case)

We do appreciate comments and suggestions.  Please send
bugs or problems to Martin Maechler <maechler@stat.math.ethz.ch>
and also  pin-ng@nau.edu or he@bahadur.stat.uiuc.edu

Martin Maechler <maechler@stat.math.ethz.ch>	http://stat.ethz.ch/~maechler/
Seminar fuer Statistik, ETH-Zentrum  LEO C16	Leonhardstr. 27
ETH (Federal Inst. Technology)	8092 Zurich	SWITZERLAND

Pin T. Ng
College of Business Administration
PO Box 15066
Flagstaff, AZ 86011-5066
vox: (928) 523-8726
fax: (928) 523-7331
e-mail: Pin.Ng@nau.edu
url: http://www.cba.nau.edu/pin-ng


old affiliations :

	Pin T. Ng			Xuming He
	Department of Economics		Department of Statistics
	University of Illinois		University of Illinois
	Champaign, IL 61820		Champaign, IL 61820

------------------------------------------------------------------------

Dec 20, 2001:
    Prof. Pin T. Ng explicitly allowed Martin Maechler
    to port and distribute this for R.

--> First public release to CRAN, 18 April 2002, version 0.9-3


See also the following files

name		       content
--------	       ------------------------------------------------
./TODO		       Martin's "to do" list (he has also ./00-FIXME)
./ChangeLog	       Martin's ``log book'' on changes done.
./R_README	       calling tree of R functions, ..
./inst/scripts/README  about Ng&He (1994)'s 3 examples
./tests/README	       explaining a bit about current "instability" problem

------------------------------------------------------------------------------
Martin Maechler <maechler@stat.math.ethz.ch>	http://stat.ethz.ch/~maechler/
Seminar fuer Statistik, ETH-Zentrum
ETH (Federal Inst. Technology)	8092 Zurich	SWITZERLAND


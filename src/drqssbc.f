C -*- mode: fortran; kept-old-versions: 12;  kept-new-versions: 20; -*-

C--- FIXME [MM] : The BLAS routine names are "s...." for single precision
C--- =====        but they are *coded* as double precision in  dblas1.f
C--- this is *not* helpful --> replaced all the following,
C--- i.e. the R-internal blas routines will be called.
C s/sscal1/dscal/  1 x
C s/sasum1/dasum/  many times
C s/saxpy1/daxpy/  many times
C s/scopy1/dcopy/  many times
C s/sdot1/ddot/    many times
C ----  ------------------------
C NOTE: Could NOT replace the srotm() and srotmg1()

C Output from Public domain Ratfor, version 1.0
C From original file drqssbc.r by Ng and He (1997)
C Pretty edited (and ratfor "bugs" fixed):
C Copyright © 2001 by Martin Maechler <maechler@stat.math.ethz.ch>
c
      subroutine drqssbc(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,
     *	   trace, e,ner,x,f,erql1n,res,indx,w,nt, nsol,sol, tl,
     *	   toler,big,eps, it, tmin,k,k0,lstart,factor)

      implicit none

C
C This is a modification of Bartels and Conn (1980) as described in
C	Koenker, R. and Ng, P. (1996)
C       "A Remark on Bartels and Conn's Linearly Constrained L1 Algorithm",
C       ACM Transaction on Mathematical Software 22, 493-495.
C
C It also contains the Parametric Linear Programming (=: PLP)
C on `tau' or `lambda' as described in
C	Ng, P. (1996)
C       "An Algorithm for Quantile Smoothing Splines",
C	Computational Statistics & Data Analysis 22, 99-118.
C
C
C     ***************
C     front end interface
C     ***************
C
C +++++ parameters +++++
C ----------------------------------------------------------
C			 input
C name	 type  subscrpt	 output	       description
C			 scratch
C ..........................................................
C nrq	 int.	 none	   in	   number of observations in the rq norm that
C				   correspond to the fidelity (may be zero)
C
C nl1	 int.	 none	   in	   number of observations in the l1 norm that
C				   correspond to roughness measure (may be zero)
C
C neqc	 int.	 none	   in	   number of equality constraints (may be zero)
C
C niqc	 int.	 none	   in	   number of inequality constraints
C				   (may be zero)
C
C niqc1	 int.	 none	   in	   part of niqc that belongs to
C				   the loo roughness measure (may be zero)
C
C nvars	 int.	 none	   in	   number of variables
C
C nact	 int.	 none	   out	   number of active equations/constraints at
C				   termination (if any, their associated column
C				   positions in  e  will be listed in	indx(1)
C				   through  indx(nact) )
C
C ifl	 int.	 none	   out	   termination code (see below)
C
C mxs	 int.	 none	   in	   maximum number of steps allowed = `maxiter'
C
C trace	 int.	 none	   in	   trace level (for info printing, see below)
C				   originally was logical `psw'
C
C e	 real	  2	   in	   equation/constraint matrix
C				   the first  nrq+nl1  columns (see note below)
C				   specify equations, the remaining
C				   columns (if any) specify constraints.
C
C ner	 int.	 none	   in	   row dimension of e
C
C x	 real	  1	   in	   starting values for the unknowns
C					(use zeros if no guess is available)
C			   out	   termination values for the unknowns
C
C f	 real	  1	   in	   equation/constraint right-hand sides
C
C erql1n real	 none	  out	   rq-l1 norm of equation
C				   residuals at termination
C
C res	 real	  1	   out	   equation/constraint
C				   residuals at termination
C
C indx	 int.	  1	   out	   index vector used to record the order in
C				   which the columns of	 e  are being processed
C
C w	 real	  1	   scr.	   working storage
C nt	 int.	  none	   out	   number of unique tau or lambda solutions
C				   while performing PLP in tau or lambda
C nsol	 int.	  none	   in	   upper limit for the number of unique
C				   tau or lambda solutions
C sol	 real	  2	   out	   matrix of solutions when performing PLP
C				   in tau or lambda
C tl	 real	  1	   in	   values of initial tau and lambda
C toler	 real	  none	   in	   tolerance used in PLP
C big	 real	  none	   in	   largest representable floating point number
C eps	 real	  none	   in	   smallest number satisfying (1. + eps) > 1.
C icyc	 int.	  none	   out	   number of cycles to achieve convergence
C tmin	 real	  none	   in	   smallest value of tau to begin PLP in tau
C k	 int.	  none	   out	   effective dimension of the model
C k0	 int.	  none	   in	   the largest effective dimension of the
C				   model allowed during PLP in lambda
C lstart real	 none	  in	   largest value of lambda to begin PLP {lambda}
C factor real	 none	  in	   factor to determine how big a step
C				   to take to the next smaller lambda
C				   during PLP {lambda}
C ----------------------------------------------------------
C
C +++++ purpose +++++
C ----------------------------------------------------------
C this subroutine solves the   nrq+nl1 by nvars
C system of equations
C
C		    (a-transpose) * x	==   b
C
C subject to the  neqc	 constraints
C
C		    (g-transpose) * x	==  h
C
C and the  niqc	 inequality constraints
C
C		    (c-transpose) * x	>=  d
C
C for the unknowns  x(1),...,x(nvars).
C
C the problem must be well-posed, nontrivial
C and overdetermined in the sense that
C
C		       nvars   >= 1
C		       nrq+nl1 >= 0
C		       neqc    >= 0
C		       niqc    >= 0
C	    nrq+nl1+neqc+niqc  >= nvars.
C
C further, no column of	 a, g  or  c  should be zero.
C if these conditions are not met, the program
C will terminate without performing any substantive
C computations.
C
C a point  x  is a solution if it minimizes the equation
C residuals from among all points which satisfy the
C constraints.	at any (nondegenerate) solution
C there will be	 nact  equations and constraints
C whose residuals
C
C      (a(i)-transpose) * x - b(i)
C
C      (g(i)-transpose) * x - h(i)
C
C and
C
C      (c(i)-transpose) * x - d(i)
C
C are zero.
C
C the columns of  (a,g,c)  corresponding to the zero residuals
C are referred to as  active columns  throughout this listing.
C the numbers of the active columns are maintained as the
C entries  1,...,nact  of the array  indx.
C
C a solution  x	 is found by minimizing a piecewise
C linear penalty function formed from the  l1
C norm of the equation residuals and the sum of the
C infeasibilities in the constraints.
C the minimization proceeds in a step-by-step
C fashion, terminating after a finite number of steps.
C
C Note that  a, g  and	c  appear transposed in the
C problem formulation.	hence it is the columns of  (a,g,c)
C which define the equations and constraints respectively.
C
C The array  e	is a composite of  a, g  and  c
C and        f  is a composite of  b, h  and  d.
C e  must contain  a  as its first    nrq+nl1  columns,
C         contain  g  as its next       neqc   columns and
C         contain  c  as its remaining	niqc   columns.  Similarly
C f  should contain  b  as its first  nrq+nl1  components,
C                    h  as its next	neqc   components, and
C                    d  as its last     niqc  components.
C ----------------------------------------------------------
C
C +++++ arrays +++++
C ----------------------------------------------------------
C     e	 is to be dimensioned	at least    n by  m,
C     x				at least    n,
C     f				at least    m,
C     res			at least    m,
C     indx			at least    m,
C     w				at least    ((3*n*n+11*n+2)/2) + (2*m).
C
C	where  n = nvars  and  m = nrq+nl1+neqc+niqc
C ----------------------------------------------------------
C
C +++++ initialization +++++
C ----------------------------------------------------------
C the user must initialize
C
C      nrq,nl1,neqc,niqc,nvars,mxs,trace,e,ner,x,f .
C
C the following are set
C and do not require initialization
C
C      nact,indx,res .
C
C the array  w	is used as scratch space.
C ----------------------------------------------------------
C
C +++++ termination codes and intermediate printing +++++
C ----------------------------------------------------------
C mxs  sets a limit on the number of minimization steps to be
C taken.
C
C upon termination  ifl	 will be set according to
C the following code ...
C
C	  ifl = 1 .... successful termination.
C
C	  ifl = 2 .... unsuccessful termination.
C		       constraints cannot be satisfied.
C		       problem is infeasible.
C
C	  ifl = 3 .... limit imposed by	 mxs  reached
C		       without finding a solution.
C
C	  ifl = 4 .... program aborted.
C		       numerical difficulties due to ill-conditioning.
C     		       (MM: might to differentiate the `kinds of ifl = 4')
C
C	  ifl = 5 .... nrq, nl1, nvars, neqc and/or
C		       niqc  have improper values
C		       or  e  contains a zero column.
C
c-- comment for next ones added by MM (from reading source):
C
C	  ifl = 6 .... nsol was too small (too many lambda's in LPL).
C                      need larger sol[] & nsol = ncol(sol)
C
C	  ifl = 7 .... both lambda < 0  & tau outside [0,1];
C                      cannot find both!
C
C
C in all cases the output parameters  x,erql1n and res
C will contain the values which they reached at termination.
C
C intermediate printing will be turned off if  trace = 0 (.false.).;
C on the other hand, increasingly more details of each minimization cycle
C will be printed if  trace  is set > 0.
C
C Args
      integer nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs, ner, trace
      integer it(2), indx(1), nt,nsol, k,k0
      double precision e(ner,1),x(1),f(1), res(1),w(1),sol(nvars+6,nsol)
      double precision tl(2),erql1n, toler,big,eps, tmin, lstart,factor

      call dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,
     *	   trace, e,ner,x,f,erql1n,res,indx,w,nt,nsol,sol, tl(1),tl(2),
     *	   toler,big,eps, it(1),it(2), tmin,k,k0,lstart,factor)
      return
      end

      subroutine dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,
     *	   trace, e,ner,x,f,erql1n,res, indx,w,nt,nsol,sol, t,lam,
     *	   toler,big,eps, icyc,totcyc, tmin,k,k0,lstart,factor)

      implicit none

C     Args
      integer nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs, ner, trace
      integer indx(1), nt,nsol, icyc, totcyc, k,k0
      double precision e(ner,1),x(1),f(1), res(1),w(1),sol(nvars+6,nsol)
      double precision erql1n, t,lam,toler,big,eps, tmin, lstart,factor

C     VAR
      logical itend,ilend,ilfix
      integer i, nrql1, ddx,grdx,grd1x,iaddc,idelc
      integer px,ptex,rrx,topx,zzx
      double precision tmax,l0,l1, tnxt,lnxt
      double precision alpha,amag,cgmag,pen,penpar,told
      double precision zero,one
C
      data zero/0.d00/
      data one/1.d00/
C

C     ////////////////	begin program  /////////////////////////
C
C     initialize ifl to 0
      ifl = 0

      nt=1
      tnxt=t
      lnxt=lam
      sol(1,nt)=t
      sol(2,nt)=lam
      itend = zero .le. t .and. t .le. one
      if(.not. itend) then
C        find tau.  Note here that tmin is passed into the subroutine
	 tmax = one - toler
	 tnxt=tmin
	 told=zero
	 sol(1,nt)=tmin
      endif
      ilend = lam .ge. zero
      ilfix = ilend
      if(.not. ilend)then
c        find lambda
	 l0 = toler
c--      FIXME?  l0 = toler means that very small (< toler) lambda won't be ?
	 l1 = (big-toler)
	 lnxt=lstart
	 told=t
	 sol(2,nt)=lstart
      endif

c                           monit0(.) --> ./monitor.c
      if(trace .gt. 0) call monit0(nrq, nl1, neqc, niqc, nvars,
     *     t, lam, itend,ilend, trace)

      if(.not.itend .and. .not.ilend)then
C        don't allow both t and lam to vary
	 ifl = 7
	 return
      endif
C
      nrql1 = nrq+nl1
C     Note: penpar is assigned outside the loop
      penpar=one/lnxt**.5

      totcyc = 0
C REPEAT
 100  continue
      call drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,px,ptex,rrx,
     *	   topx,zzx,icyc,ifl,e,ner,amag,cgmag,lnxt)
c     penpar not touched above, but ifl = 2 (or = 5 if "fail"), icyc= -1;
c      -->  penpar /= 8  below
      call dnewpen(iaddc,idelc,nact,nrql1,neqc,niqc,nvars,ifl,e,ner,x,f,
     *	   res,w(ptex),alpha,penpar,indx)

      if(trace .ge. 2) call monit11(nt, nact,
     *     amag, cgmag, ifl, trace)


c  Repeat
 20   continue
      call drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,ifl,
     *	   mxs,e,ner,x,f,res,w(grdx),erql1n,pen,penpar,indx,w(zzx),
     *	   nvars,w(ddx),w(rrx),w(topx),tnxt,eps,w(grd1x))

      if(trace .ge. 3) call monit2(nact,icyc,trace,
     *     x,alpha, erql1n,pen,penpar,indx)
c               monit2(.) --> ./monitor.c

c     Find descent direction (or discover optimality):
      call drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,f,res,
     *	   w(grdx),w(px),erql1n,amag,cgmag,penpar,indx,w(zzx),nvars,
     *	   w(ddx),w(rrx),w(topx),tnxt,big,eps)
c     Piecewise line search :
      call drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,res,
     *	   w(grdx),w(px),w(ptex),alpha,penpar,indx,w(topx),tnxt,
     *	   big,eps,w(grd1x),idelc)
      if(ifl .eq. 0) goto 20
c  Until (ifl != 0)
      continue

      totcyc = totcyc + icyc
      if(trace .ge. 2) call monit12(icyc, ifl, trace)

      if(.not.(itend.and.ilend) .and.
     *	   (ifl .ne.2 .or. cgmag+penpar*amag .eq. cgmag)) then
c        Next Lambda or Tau, i.e. number [nt + 1]
	 call drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,nact,nvars,
     *	      nrq,nl1,neqc,niqc,niqc1,w(zzx),nvars,w(rrx),w(grdx),
     *	      w(px), w(topx),w(topx+nvars),ifl,idelc,iaddc,icyc,
     *	      alpha,amag, cgmag,trace,penpar,nsol,sol,x,ilfix,l0,l1,
     *	      tnxt,lnxt,toler, erql1n,eps,big,told,k0,factor)

         if(trace .ge. 2) call monit13(nt, nact, itend, ilend,
     *        tnxt, lnxt, ifl, trace)

C        update penpar to the sqrt of next lambda
	 penpar=one/lnxt**.5
	 if(ifl .ne. 0)then
	    goto 30
	 endif
      else
	 if(ifl .ne.2 .or. cgmag+penpar*amag .eq. cgmag)then
	    goto 30
	 endif
      endif

      goto 100
C UNTIL ......

C
C     compute the effective dimensionality:
 30   continue
      k = 0
      do 40 i=1,nact
	 if(indx(i) .le. nrq .or.
     1	   (indx(i) .gt. nrq+nl1 .and. indx(i) .le. nrq+nl1+neqc) .or.
     2	    indx(i) .gt. nrq+nl1+neqc+niqc1) then
	    k = k+1
	 endif
 40   continue
      return
      end
c     --- end dcrql1lt { = only main routine of drqssbc() }

      subroutine drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,
     1	   nact,nvars,nrq, nl1,neqc,niqc,niqc1, zz,nzzr,rr,
     2	   a,aa,b,bb, ifl,idelc,iaddc,icyc, alpha,amag,cgmag,
     3	   trace,penpar,nsol,sol, x,ilfix, l0,l1,tnxt,lnxt,
     4	   toler,erql1n,eps,big,told,k0,factor)
C
C     ***************
C     perform parametric programming in "lambda" and "tau",
C     compute `Next Lambda and Tau' (= `nlt'), i.e. number [nt + 1]
C     ***************
C
      implicit none

      integer nt, ner,indx(1),nact,nvars,nrq,nl1,neqc,niqc,niqc1, nzzr
      integer ifl,idelc,iaddc,icyc, nsol,  k0, trace
      logical ilfix

      double precision tmin,tmax,res(1),e(ner,1),f(1),zz(nzzr,1),rr(1)
      double precision a(1),aa(1),b(1),bb(1), alpha,amag,cgmag
      double precision penpar, sol(nvars+6,nsol), x(1),l0,l1
      double precision toler,erql1n,eps,big,told, factor
C VAR
      logical fail
      integer nrql1, isave, k,nactp1,nallq,nalqp1,ncols,nqnp1,i,j,ix
      double precision lnxt,lamb, thet,tnxt,tmp
      double precision test,prod, fidel,penal, sgn, wgt
      double precision zero,one,two

C
C     ////////////////	begin program  /////////////////////////
C
      data one/1.d00/
      data zero/0.d00/
      data two/2.d00/
C
      nrql1=nrq+nl1
      nallq=nrql1+neqc
      nalqp1=nallq+1
      ncols=nallq+niqc
      nqnp1=nrql1+1
      nactp1=nact+1
      thet = tnxt
      lamb = lnxt
c     next tau *or* next lambda :
      if(ilfix)then
	 tnxt = one+eps
      else
	 lnxt = l0
      endif
      if(ifl.eq.1.or.ifl.eq.3)then
	 call dcopy(nvars,zero,0,a,1)
	 call dcopy(nvars,zero,0,b,1)
c ORIGINAL code had `nacpt' which is *not* declared (i.e. := 0 for some...)
c	  if(nacpt.le.ncols)then
         if(nactp1 .le. ncols)then
	    do 30 i = nactp1,ncols
	       ix = indx(i)
	       sgn = dsign(one,res(ix))
	       test = dabs(f(ix))
	       do 20 j = 1,nvars
		  prod = dabs(e(j,ix)*x(j))
		  if(prod .gt. test)then
		     test = prod
		  endif
 20	       continue

	       test = eps*dsqrt(dfloat(nvars))*test
	       if(dabs(res(ix)) .lt. test)then
		  sgn = zero
	       endif
	       if(ilfix)then
		  if(ix.le.nrq)then
c     a := a + (one+sgn)*e(1,ix) :
		     call daxpy(nvars,(one+sgn),e(1,ix),1,a,1)
c     b := b -two* e(1,ix) :
		     call daxpy(nvars,-two,e(1,ix),1,b,1)
		  else
		     if(ix.le.nallq.or.sgn.le.zero)then
c     a := a + sgn*e(1,ix) :
			call daxpy(nvars,sgn,e(1,ix),1,a,1)
		     endif
		  endif
	       else if(ix.le.nrq)then
                  call daxpy(nvars,(one-two*thet+sgn),e(1,ix),1,a,1)
               else if(ix.le.nrql1)then
                  call daxpy(nvars,sgn/lamb,e(1,ix),1,b,1)
c ORIGINAL code had `allq' which is *not* declared (i.e. := 0 for some...)
C              else if(ix.le. allq.or.sgn.le.zero)then
               else if(ix.le.nallq.or.sgn.le.zero)then
                  call daxpy(nvars,sgn,	   e(1,ix),1,a,1)
	       endif
 30	    continue
	 endif

	 call dzdrgnv(nvars,nact,zz,nzzr,rr,a,aa,fail,big)
	 if(fail)then
	    ifl = 4
	 else
	    call dzdrgnv(nvars,nact,zz,nzzr,rr,b,bb,fail,big)
	    if(fail)then
	       ifl = 4
	    else
	       do 50 i = 1,nact
		  ix = indx(i)
C     a check for small bb(i) is implemented
C     to avoid floating point overflow
		  test = dabs(f(ix))
		  do 56 j = 1,nvars
		     prod = dabs(e(j,ix)*x(j))
		     if(prod .gt. test)then
			test = prod
		     endif
 56		  continue

		  test = eps*dsqrt(dfloat(nvars))*test
		  if(ix.le.nrq)then
		     if(ilfix)then
			tmp = (two+aa(i))/(two-bb(i))
			if(tmp.lt.tnxt.and.tmp.ge.thet)then
			   tnxt = tmp
			   isave = i
			else
			   tmp = aa(i)/(two-bb(i))
			   if(tmp.lt.tnxt.and.tmp.ge.thet)then
			      tnxt = tmp
			   endif
			endif
		     else
			tmp = (two*thet - aa(i))/bb(i)
			if(dabs(bb(i)).lt.test)then
C				avoid bb near zero
			   tmp = big
			endif
			if(tmp.gt.lnxt.and.tmp.lt.lamb)then
			   lnxt = tmp
			   isave = i
			else
			   tmp = (two*thet - two - aa(i))/bb(i)
			   if(dabs(bb(i)).lt.test)then
C				avoid bb near zero
			      tmp = big
			   endif
			   if(tmp.gt.lnxt.and.tmp.lt.lamb)then
			      lnxt = tmp
			      isave = i
			   endif
			endif
		     endif

		  else if(ix.le.nallq)then
		     if(ilfix)then
			tmp = (one-aa(i))/bb(i)
			if(dabs(bb(i)).lt.test)then
C			   avoid bb near zero
			   tmp = big
			endif
			if(tmp.lt.tnxt.and.tmp.ge.thet)then
			   tnxt = tmp
			   isave = i
			else
			   tmp = -(aa(i)+one)/bb(i)
			   if(dabs(bb(i)).lt.test)then
C			      avoid bb near zero
			      tmp = big
			   endif
			   if(tmp.lt.tnxt.and.tmp.ge.thet)then
			      tnxt = tmp
			      isave = i
			   endif
			endif
		     else
			tmp = -aa(i)*lamb/(bb(i)*lamb+one)
			if(tmp.gt.lnxt.and.tmp.lt.lamb)then
			   lnxt = tmp
			   isave = i
			else
			   tmp = -aa(i)*lamb/(bb(i)*lamb-one)
			   if(tmp.gt.lnxt.and.tmp.lt.lamb)then
			      lnxt = tmp
			      isave = i
			   endif
			endif
		     endif
		  else if(ilfix)then
		     tmp = -aa(i)/bb(i)
		     if(dabs(bb(i)).lt.test)then
C			avoid bb near zero
			tmp = big
		     endif
		     if(tmp.lt.tnxt.and.tmp.ge.thet)then
			tnxt = tmp
			isave = i
		     endif
		  else
		     tmp = -aa(i)/bb(i)
		     if(dabs(bb(i)).lt.test)then
C			avoid bb near zero
			tmp = big
		     endif
		     if(tmp.gt.lnxt.and.tmp.lt.lamb)then
			lnxt = tmp
			isave = i
		     endif
		  endif
 50	       continue
	    endif
	 endif
      endif
C
C     compute the effective dimensionalty, fidelity and penalty
C
      k = 0
      do 2 i=1,nact
	 if(indx(i).le. nrq .or.
     *	   (indx(i).gt. nrq+nl1 .and. indx(i).le. nrq+nl1+neqc) .or.
     *	    indx(i).gt. nrq+nl1+neqc+niqc1) then
	    k = k+1
	 endif
 2    continue

C     set the lower stopping criterion for lambda to be either k>=k0
C     or when lnxt < 0
      fidel = zero
      penal = zero
      if(k.ge.k0 .or. lnxt-lnxt*10.0d0**(factor-4.0d0) .lt.zero)then
	 l0 = lnxt+eps
      endif
      do 8 i=nactp1,ncols
	 ix = indx(i)
	 tmp = res(ix)
	 wgt = dsign(one,tmp)
	 if(ix.le.nrq)then
	    wgt = wgt+ one-two*told
	    fidel = fidel+wgt*tmp
	 else
	    if(ix.le.nrql1)then
	       penal = penal+dabs(tmp)
	    endif
	 endif
 8    continue

      nt = nt+1
c     --------- compute next tau or lambda; check for nt >= nsol below
      if(ilfix)then
c        ----  new  tau
	 if((ifl.eq.1.or.ifl.eq.3) .and. tnxt.lt.tmax)then
	    sol(1,nt) = tnxt
	    sol(2,nt) = lamb
	    sol(3,nt-1)=dble(ifl)
	    sol(4,nt-1) = fidel
	    sol(5,nt-1) = penal/lamb
	    sol(6,nt-1) = k
	    told = tnxt
	    if(indx(isave) .le. nrql1)then
	       tnxt = tnxt+10.0d0**(factor-7.0d0)*amag
	    else
	       tnxt = tnxt+10.0d0**(factor-7.0d0)*cgmag
	    endif
	    call dcopy(nvars,x,1,sol(7,nt-1),1)
	    idelc = 0
	    iaddc = nact
	    icyc = -1
	    ifl = 0
	    alpha = zero
	 else
	    sol(1,nt) = tnxt
	    sol(2,nt) = lamb
	    sol(3,nt-1)=dble(ifl)
	    sol(4,nt-1) = fidel
	    sol(5,nt-1) = penal/lamb
	    sol(6,nt-1) = k
	    if((ifl.eq.1.or.ifl.eq.3) .and. tnxt .ge.tmax)then
	       sol(1,nt) = one
	       sol(1,1) = zero
	       sol(4,1) = zero
	       sol(5,1) = zero
	       sol(6,1) = two
	       sol(2,nt) = lamb
	       sol(3,nt) = sol(3,nt-1)
	       sol(4,nt) = zero
	       sol(5,nt) = zero
	       sol(6,nt) = two
	       call dcopy(nvars,x,1,sol(7,nt),1)
	    endif
	    call dcopy(nvars,x,1,sol(7,nt-1),1)
	 endif

      else
c     ---- ilfix is FALSE : new  lambda
         if((ifl.eq.1.or.ifl.eq.3) .and. lnxt.gt.l0)then
            sol(1,nt) = thet
            sol(2,nt) = lnxt
            sol(3,nt-1)=dble(ifl)
            sol(4,nt-1) = fidel
            sol(5,nt-1) = penal/lamb
            sol(6,nt-1) = k
            lnxt = lnxt-lnxt*10.0d0**(factor-4.0d0)
            call dcopy(nvars,x,1,sol(7,nt-1),1)
            idelc = 0
            iaddc = nact
            icyc = -1
            ifl = 0
            alpha = zero
         else
            sol(1,nt) = thet
            sol(2,nt) = lnxt
            sol(3,nt-1)=dble(ifl)
            sol(4,nt-1) = fidel
            sol(5,nt-1) = penal/lamb
            sol(6,nt-1) = k
            if((ifl.eq.1.or.ifl.eq.3) .and. lnxt .le.l0)then
               sol(1,nt) = thet
               sol(2,nt) = l0
               sol(3,nt) = sol(3,nt-1)
               sol(4,nt) = sol(4,nt-1)
               sol(5,nt) = sol(5,nt-1)
               sol(6,nt) = sol(6,nt-1)
               call dcopy(nvars,x,1,sol(7,nt),1)
            endif
            call dcopy(nvars,x,1,sol(7,nt-1),1)
         endif
      endif

      if(nt.ge.nsol)then
c        have no extra space in sol(,)
         ifl = 6
      endif

C     remove lambda from e for next iteration
      if(nrql1.ge.nrq+1)then
	 do 80 i=nrq+1,nrql1
	    do 84 j=1,ner
	       e(j,i) = e(j,i)/lamb
 84	    continue
 80	 continue
      endif
      return
      end
c     --- drql1nlt()

      subroutine drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,
     *	   px,ptex,rrx,topx,zzx,icyc,ifl,e,ner,amag,cgmag,lam)
C
      implicit none

      integer ddx,grdx,icyc,ifl,neqc,nrql1,ner,nrq,nl1,grd1x
      integer niqc,nvars,px,ptex,rrx,topx,zzx
      double precision amag,cgmag,e(ner,1),lam
C
C     ***************
C     crql1  version.
C
C     set up the program
C     parameters and indices.
C     ***************
C
C     +++++++++++++++
C     system routines  dabs
C     +++++++++++++++
C
      integer i,j,ncols,nqnp1
      double precision tmp,zero
C
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
C     ***************
C     check validity of problem dimensions
C     ***************
C
      nrql1=nrq+nl1
      ncols = nrql1+neqc+niqc
      if(nvars.lt.1 .or. neqc.lt.0 .or. niqc.lt.0 .or. nrql1.lt.0 .or.
     *     ncols.lt.nvars .or. ner.lt.nvars) then
	 ifl = 5
         return
      endif
c-    ELSE

C     ***************
C     set up indices for the temporary storage vector  w.
C     ***************
C
      nqnp1 = nrql1+1
      grdx = 1
      grd1x = grdx+nvars
      px = grd1x+nvars
      ptex = px+nvars
      ddx = ptex+ncols
      rrx = ddx+nvars
      zzx = rrx+(((nvars+1)*(nvars+2))/2)
      topx = zzx+nvars*nvars
C
C     ***************
C     update e with lambda only if ifl!=2, i.e. update only for the new lambda
C     ***************
C
      if( ifl.ne.2)then
         do 40 i=nrq+1,nrql1
            do 42 j=1,ner
               e(j,i)=e(j,i)*lam
 42         continue
 40      continue
      endif

C
C     ***************
C     amag  is a rough estimate of the norm of	a.
C     cgmag  is a rough estimate of the norm of	 (g,c).
C     together they are used to determine when the
C     penalty parameter is too small and when the
C     restricted gradient is zero.
C     ***************
C
      amag = zero
      cgmag = zero
      if(1.le.nrql1)then
         do 46 j = 1,nrql1
            tmp = zero
            do 48 i = 1,nvars
               tmp = tmp+dabs(e(i,j))
 48         continue
            if(tmp.le.zero)then
               ifl = 5
               return
            endif
            if(tmp.gt.amag)then
               amag = tmp
            endif
 46      continue
      endif

      if(nqnp1.le.ncols)then
         do 56 j = nqnp1,ncols
            tmp = zero
            do 58 i = 1,nvars
               tmp = tmp+dabs(e(i,j))
 58         continue
            if(tmp.le.zero)then
               ifl = 5
               return
            endif
            if(tmp.gt.cgmag)then
               cgmag = tmp
            endif
 56      continue
      endif

C     ***************
C     initialize  ifl,icyc
C     ***************
C
      ifl = 2
      icyc = -1

      return
      end

      subroutine dnewpen(iaddc,idelc,nact,neqns,neqc,niqc,nvars,ifl,
     *	   e,ner,x,f,res,pte,alpha,penpar,indx)
C
      implicit none

      integer iaddc,idelc,ifl,indx(1),nact
      integer neqc,neqns,ner,niqc,nvars
      double precision alpha,e(ner,1),f(1),penpar,pte(1),res(1),x(1)
C
C     ***************
C     cl1  version.
C
C     begin a round of minimization steps
C     with a new penalty parameter value.
C     ***************
C
C     +++++++++++++++
C     blas  ddot
C     +++++++++++++++
C
      integer i,ncols
      double precision oct,zero
C
      double precision ddot
C
      data zero/0.0d+00/
      data oct/8.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
C     ***************
C     set penalty parameter value.
C     erase record of active equation/constraints.
C     ***************
C
      if(ifl.eq.2)then
	 ncols = neqns+neqc+niqc
	 ifl = 0
	 nact = 0
	 iaddc = 0
	 idelc = 0
	 alpha = zero
	 penpar = penpar/oct
C
C     ***************
C     initialize  indx,res,pte,indx
C     ***************
C
	 do 10 i = 1,ncols
	    res(i) = ddot(nvars,e(1,i),1, x,1) - f(i)
	    pte(i) = zero
	    indx(i) = i
 10	 continue
      endif
      return
      end

      subroutine drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,
     *	   ifl,mxs,e,ner,x,f,res,grd,erql1n,pen,penpar,indx,zz,nzzr,
     *	   dd, rr,w, theta,eps,grd1)
C
C     ***************
C     crql1  version.
C
C     preparation for next minimization step.
C     ***************
C
      implicit none
c Args
      integer iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc
      integer ifl,mxs, ner, indx(1), nzzr
      double precision e(ner,1),x(1),f(1),res(1),grd(1)
      double precision erql1n,pen,penpar, zz(nzzr,1), dd(1),rr(1),w(1)
      double precision theta,eps,grd1
C VAR
      integer nrql1,nallq,ncols
C
C     ***************
C     determine the active equations and active
C     constraints.  compute residuals and function value.
C     update the  z*d*r	 decomposition.
C     ***************
C
      nrql1 = nrq+nl1
      nallq = nrql1+neqc
      ncols = nallq+niqc
      if(ifl.eq.0)then
	 icyc = icyc+1
	 if(icyc.gt.mxs)then
	    ifl = 3
	 else
	    call ddelcol1(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,indx)
	    call dresid	 (iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
	    call daddcol (iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,e,ner,
     *		 indx,w,eps)
	    call drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,
     *		 res,grd,erql1n,pen,penpar,indx,theta,grd1)
	 endif
      endif
      return
      end

      subroutine drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,
     *	   e,ner,x,f,res,grd,p,erql1n,amag,cgmag,penpar,indx,zz,nzzr,
     *	   dd,rr,w, theta,big,eps)
C
      integer idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl, ner,indx(1),nzzr
      double precision e(ner,1),x(1),f(1),res(1),grd(1),p(1),
     *	   erql1n,amag,cgmag,penpar
      double precision zz(nzzr,1),dd(1),rr(1),w(1), theta,big,eps
C
C     ***************
C     crql1  version.
C
C     determine descent direction  p
C     (or discover optimality)
C     ***************
C
C     +++++++++++++++
C     system routines  dabs
C
C     blas  dasum,dcopy,dscal
C
C     eps  is the smallest positive number which
C     satisfies	  (1.0 + eps) .gt. 1.0	 in the
C     precision of the arithmetic being used.
C     (alternatively, for less strict zero checking,
C     eps  can be set to a user-specified tolerance.)
C     +++++++++++++++
C
      logical fail
      integer coefx,i,ix,nallq,nalqp1,ncols,nqnp1,topx,nrql1
      double precision grdnrm,one,pnrm,prod,test,zero
C
      double precision dasum
C
C     +++++++++++++++
C     the following declarations are necessary
C     for portability when  dcopy  is used, as
C     it is below, to fill arrays with a single value
C     (one=unity  and  zero=zip	 in this case).
C     +++++++++++++++
C
      double precision unity(1),zip(1)
      equivalence(one,unity(1)),(zero,zip(1))
C
      data one/1.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      idelc = 0
      if(ifl .ne. 0) return
      nrql1 = nrq+nl1

      nallq = nrql1+neqc
      nalqp1 = nallq+1
      ncols = nallq+niqc
      nqnp1 = nrql1+1
      coefx = 1
      topx = coefx+nvars
C
C     ***************
C     project the negative of the restricted gradient
C     onto the orthogonal complement of the space
C     spanned by the active columns.
C     ***************
C
      call dzdrpoc(nvars,nact,zz,nzzr,dd,grd,p,fail)
      if(fail)then
         ifl = 4
         return
      endif
C    else
      call dscal(nvars,-one,p,1)
      pnrm = dasum(nvars,p,1)
      grdnrm = dasum(nvars,grd,1)
C
C     ***************
C     if the projection is not zero,
C     it will serve as a descent direction.
C
C     otherwise find the representation of
C     the restricted gradient as a linear
C     combination of the active columns.
C     the coefficients of the linear combination
C     are to be stored in the array  coef
C     (that is, in  w(coefx),...,w(coefx+nact-1)).
C     ***************
C
      if(pnrm.le.eps*(amag*penpar+cgmag))then
         if(nact.ne.0)then
            call dzdrgnv(nvars,nact,zz,nzzr,rr,grd,w(coefx), fail,big)
            if(fail)then
               ifl = 4
               return
            endif
c--	    else
C
C     ***************
C     convert the coefficients of the linear
C     combination into a descent direction  p ,
C     or determine optimality.
C
C     if the optimality test is not satisfied,
C     drql1gv  will indicate an equation/constraint
C     to be deleted from activity by the value
C     of  idelc.  for optimality,  idelc=0.
C     ***************
C
            call drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,
     *           grd,w(coefx),penpar,indx,theta,eps)
            pnrm = zero
            if(idelc.ne.0)then
               call dzdrgit(nvars,nact,zz,nzzr,rr,w(coefx),p,
     *              fail,w(topx),big,eps)
               if(fail)then
                  ifl = 4
                  return
               else
                  pnrm = dasum(nvars,p,1)
               endif
            endif
C
C     ***************
C     if a descent direction  p	 could have been found,
C     it has been obtained by this point in the program.
C
C     check for optimality.
C
C     pnrm  has been set exactly zero
C     after the call to subroutine  drql1gv
C     if the optimality conditions are satisfied.
C     the check below has been made somewhat
C     complicated to allow for the rare event that
C     the restricted gradient is zero and no
C     columns are active,  or that the	rq  norm of
C     (a-transpose) * x - f
C     is computationally zero.
C     (the call to the subroutine  refine
C     may be omitted, if desired.)
C     ***************
            if(pnrm.gt.eps*(amag*penpar+cgmag))then
               do 10 i = 1,nrql1
                  test = dabs(f(i))
                  do 12 ix = 1,nvars
                     prod = dabs(e(ix,i)*x(ix))
                     if(prod.gt.test)then
                        test = prod
                     endif
 12               continue
                  if(dabs(res(i)).gt.eps*test)then
                     return
                  endif
 10            continue
            endif
         endif

         ifl = 1
C     call drql1rf(nact,nrq,nl1,ncols,nvars,ifl,e,ner,x,f,erql1n,res,indx,zz,
C		nzzr, rr ,w,theta,big,eps)
         if(ifl.eq.1)then
C
C     ***************
C     if the problem has constraints,
C     check feasibility.
C     ***************
C
            if(nqnp1.le.ncols)then
               do 20 i = nqnp1,ncols
                  test = dabs(f(i))
                  do 24 ix = 1,nvars
                     prod = dabs(e(ix,i)*x(ix))
                     if(prod.gt.test)then
                        test = prod
                     endif
 24               continue
C NOTE: the criterion for checking feasibility is relaxed
C      by (eps * test)^.5 rather than eps*test
                  test = (eps*test)**.5
C     test = eps*test
                  if(i.gt.nallq)then
                     if(res(i).lt.(-test))then
                        ifl = 2
                        return
                     endif
                  else
                     if(dabs(res(i)).gt.test)then
                        ifl = 2
                        return
                     endif
                  endif
 20            continue
            endif
         endif
      endif
      return
      end

      subroutine drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,
     *	   x,res,grd,p,pte,alpha,penpar,indx,alf,theta,
     *	   big,eps,grd1,idelc)
C
      integer iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl, ner,indx(1),idelc
      double precision e(ner,1),x(1),res(1),grd(1),p(1),pte(1)
      double precision alpha,penpar, alf(1),theta,big,eps, grd1(1)

      integer nrql1
      double precision sgn1
C
C     ***************
C     cl1  version.
C
C     piecewise linear line search.
C     ***************
C
C     +++++++++++++++
C     system routines dabs,dsign
C
C     blas  dasum,daxpy,ddot
C
C     eps  is the smallest positive number which
C     satisfies	  (1.0 + eps) .gt. 1.0	 in the
C     precision of the arithmetic being used.
C     (alternatively, for less strict zero checking,
C     eps  can be set to a user-specified tolerance.)
C
C     big  is the largest positive number
C     which can be represented in the
C     precision of the arithmetic being used.
C     +++++++++++++++
C
      integer i,iin,ix,jx,nactp1,nallq,ncols,num,numnac
      double precision den,grdnrm,one,pnrm,ptg,ptg1
      double precision ratio,resid,tmp,two,zero
C
      double precision dasum,ddot,etp
C
      data one/1.0d+00/
      data two/2.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
C     ***************
C     this routine determines all of the ratios	 alf
C     of the form
C     -res(i)/((e(.,i)-transp)*p),
C     for  i = k+1,...,mpl
C     which are nonnegative and hence indicate distances
C     from the point  x	 to breakpoints which will
C     be encountered in travel along direction	p.
C     the index vector	indx  is rearranged so that
C     its  k+1	through	 num  components correspond to
C     these nonnegative ratios.
C     the results are heaped so that the  alf  values can
C     be inspected in order from smallest to largest.
C     the breakpoint  alpha  giving the minimum objective
C     function value is found, and  x  is
C     adjusted to  x + alpha*p .
C
C     the inner products  (e(.,i)-transpose)*p	are saved
C     for later use in updating the residual values.
C     ***************
C
      alpha = zero
      if(ifl.eq.0)then
	 nrql1 = nrq+nl1
	 nallq = nrql1+neqc
	 ncols = nallq+niqc
	 nactp1 = nact+1
	 num = 0
	 if(1.le.nact)then
	    do 18 i = 1,nact
	       ix = indx(i)
	       pte(ix) = ddot(nvars,e(1,ix),1,p,1)
C     update the correct gradient
 18	    continue
	 endif
	 if(nactp1.le.iaddc)then
	    do 22 i = nactp1,iaddc
	       ix = indx(i)
	       etp = ddot(nvars,e(1,ix),1,p,1)
	       sgn1 = dsign(one,etp)
	       if(ix.le.nrq)then
		  sgn1 = one-two*theta + sgn1
	       endif
	       if(ix.le.nallq .or. sgn1.le.zero)then
		  if(ix.le.nrql1)then
		     sgn1 = sgn1*penpar
		  endif
		  call daxpy(nvars,sgn1,e(1,ix),1,grd1,1)
	       endif
 22	    continue
	 endif
	 if(idelc.ne.0)then
	    ix=indx(idelc)
	    etp=ddot(nvars,e(1,ix),1,p,1)
	    sgn1=dsign(one,etp)
	    if(ix.le.nrq)then
	       sgn1 = one-two*theta + sgn1
	    endif
	    if(ix.le.nallq .or. sgn1.le.zero)then
	       if(ix.le.nrql1)then
		  sgn1=sgn1*penpar
	       endif
	       call daxpy(nvars,sgn1,e(1,ix),1,grd1,1)
	    endif
	 endif
	 if(nactp1.gt.ncols)then
	    ifl = 1
	 else
	    do 40 i = nactp1,ncols
	       ix = indx(i)
	       resid = res(ix)
	       den = ddot(nvars,e(1,ix),1,p,1)
	       pte(ix) = den
	       if(dsign(one,resid).ne.dsign(one,den) .or.
     +		    resid.eq.zero) then
		  resid = dabs(resid)
		  den = dabs(den)
		  if(den.lt.one)then
		     if(resid.ge.den*big)then
			goto 40
		     endif
		  endif
		  ratio = resid/den
		  num = num+1
		  numnac = num+nact
		  jx = indx(numnac)
		  indx(numnac) = ix
		  indx(i) = jx
		  alf(num) = ratio
	       endif
 40	    continue

	    if(num.le.0)then
	       ifl = 2
	    else
C
C     ***************
C     heap the positive ratios
C     ***************
C
	       call ddkheap(.true.,num,indx(nactp1),alf)
C
C     ***************
C     travel along  p  until no further decrease in the
C     penalty function is possible
C     ***************
C
	       iin = num
	       ptg = ddot(nvars,grd,1,p,1)
	       ptg1 = ddot(nvars,grd1,1,p,1)
	       pnrm = dasum(nvars,p,1)
	       grdnrm = dasum(nvars,grd1,1)
	       do 50 i = 1,num
		  ix = indx(nactp1)
		  if(res(ix).eq.zero)then
		     tmp = zero
		  else
		     tmp = -dsign(one,res(ix))
		  endif
		  if(ix.le.nallq)then
		     tmp = tmp*two
		  endif
		  if(ix.le.nrql1)then
		     tmp = tmp*penpar
		  endif
		  ptg1 = ptg1+tmp*pte(ix)
		  if(ptg1.ge.(-eps)*grdnrm*pnrm)then
		     go to 140
		  endif
		  call ddkheap(.false.,iin,indx(nactp1),alf)
 50	       continue
	       ifl = 2
	       return

 140	       iaddc = nactp1
C
C     ***************
C     adjust  x	 to  x + alpha*p
C     ***************
C
	       alpha = alf(1)
	       call daxpy(nvars,alpha,p,1,x,1)
	    endif
	 endif
      endif
      return
      end

c Used to be called from drql1fp() --- NOT USED ANYMORE :
      subroutine drql1rf(nact,nrq,nl1,ncols,nvars,ifl, e,ner, x,f,
     *	   erql1n,res,indx,zz,nzzr,rr,w,theta,big,eps)
C
      integer nact,nrq,nl1,ncols,nvars,ifl,ner, indx(1),nzzr
      double precision e(ner,1),x(1),f(1),erql1n,res(1),zz(nzzr,1)
      double precision rr(1),w(1),theta, big,eps
C
C     ***************
C     a routine for refining the solution
C     produced by  crql1.
C
C     (this routine may be omitted if desired.)
C     ***************
C
C     +++++++++++++++
C     system routines  dabs
C
C     blas  ddot
C     +++++++++++++++
C
      logical fail
      integer i, nrql1
      double precision zero,one,two, wgt,tmp,ddot
C
      data zero/0.0d+00/
      data one/1.0d+00/
      data two/2.0d+00/
C
C     /////////////// begin program ///////////////
C
      fail = .false.
      nrql1 = nrq+nl1
      if(nact.ne.0)then
	 if(fail)then
	    ifl = 4
	 else
	    erql1n = zero
	    do 10 i = 1,ncols
	       tmp = ddot(nvars,e(1,i),1,x,1)-f(i)
	       wgt = dsign(one,tmp)
	       if(i.le.nrq)then
		  wgt = one-two*theta+wgt
	       endif
	       res(i) = tmp
	       if(i.le.nrql1)then
		  erql1n = erql1n+wgt*tmp
	       endif
 10         continue
	 endif
      endif
      return
      end
C----------------------------------------------------------------------

C
C     ---------------
C     third level subroutines --
C	   ddelcol1,dresid,addcol,drql1obj,getv
C     ---------------
C
      subroutine ddelcol1(iaddc,idelc,nact,nrow,zz,nzzr,dd,rr,indx)
C
      integer indx(1),nzzr,nact,idelc,iaddc,nrow
      double precision dd(1),rr(1),zz(nzzr,1)
C
C     ***************
C     cl1  version of  ddelcol1.
C
C     this routine administers the deletion of the column
C     indicated by the value of idelc
C     from an  nrow by nact   z*d*r   decomposition.
C     note that the value of idelc
C     is the number of a column in the decomposition
C     rather than a number which refers to
C     a column in the matrix  e.
C     (the  e-column  numbers corresponding to
C     the columns of the factorization are to be
C     found in	 indx(1),...,indx(nact) .
C     the contents of	indx(nact+1),...,indx(iaddc)
C     indicate columns of  e  which are slated for
C     addition to the decomposition.)
C     the vector  indx	 is rearranged by
C     permuting the element which corresponds to
C     the deletion out to the	iaddc-th  position.
C     nact  and	 iaddc	are decreased accordingly.
C     ***************
C
      integer i,idlp1,ixdlc
      logical fail
C
C     /////////////////	 begin program	//////////////////
C
      if(idelc.ne.0)then
	 idlp1 = idelc+1
	 ixdlc = indx(idelc)
	 do 10 i = idlp1,iaddc
	    indx(i-1) = indx(i)
 10	 continue
	 indx(iaddc) = ixdlc
	 iaddc = iaddc-1
	 call dzdrcou(nrow,nact,zz,nzzr,dd,rr,idelc,fail)
	 idelc = ixdlc
      endif
      return
      end

      subroutine dresid(iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
C
      integer iaddc,indx(1),nact,ncols,ner,nvars
      double precision e(ner,1),f(1),res(1),x(1)
C
C
C     ***************
C     compute the residuals
C     (e(.,ix)-transp)*x - f(ix)  .
C     the residuals are stored in the array  res.
C     indx  is rearranged so that the zero residuals
C     correspond to  indx(1),...,indx(iaddc)  .
C     ***************
C
C     +++++++++++++++
C     system routines  dabs,idint,dfloat,dsqrt
C
C     blas  ddot
C
C     eps  is the smallest positive number which
C     satisfies	  (1.0 + eps) .gt. 1.0	 in the
C     precision of the arithmetic being used.
C     (alternatively, for less strict zero checking,
C     eps  can be set to a user-specified tolerance.)
C     +++++++++++++++
C
      integer i,iadp1,idummy,irand,ix,j,nactp1
      double precision eps,prod,temp,test,tol,zero
C
      double precision ddot,dunif01
C
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
C     tol = eps*dsqrt(dfloat(nvars))
      tol = eps
      nactp1 = nact+1
      if (1.le.iaddc) then
C
C     ***************
C     zero out all residuals known to be zero.
C     ***************
C
	 do 10 i = 1,iaddc
	    ix = indx(i)
	    res(ix) = zero
 10	 continue
      endif
C
C     ***************
C     compute the remaining residuals.
C     detect any more residuals which
C     are computationally zero, and
C     set them exactly zero.  their
C     associated indices are permuted
C     so that they are stored in
C     indx(nact+1),...,nact(iaddc).
C
C     (a fairly tight zero check is used.
C     it is far less expensive in running
C     time to neglect an extra zero
C     residual than to accept it and risk
C     invoking the anti-cycling
C     mechanisms in the program.
C     the accuracy of the solution as
C     finally determined is not affected.)
C     ***************
C
      iadp1 = iaddc+1
      if(iadp1.le.ncols)then
	 do 20 i = iadp1,ncols
	    ix = indx(i)
	    temp = ddot(nvars,e(1,ix),1,x,1)-f(ix)
	    test = dabs(f(ix))
	    do 22 j = 1,nvars
	       prod = dabs(e(j,ix)*x(j))
	       if(prod.gt.test)then
		  test = prod
	       endif
 22	    continue
	    test = tol*test
	    if(dabs(temp).gt.test)then
	       res(ix) = temp
	    else
c              set residual to 0
	       iaddc = iaddc+1
	       indx(i) = indx(iaddc)
	       indx(iaddc) = ix
	       res(ix) = zero
	    endif
 20	 continue
      endif

C
C     ***************
C     if any new zero residuals have been found, randomize their
C     ordering as an anti-cycling device for addcol() :
C     ***************
C
      if(iaddc.gt.nactp1)then
	 do 90 i = nactp1,iaddc
	    irand = i+ifix(float(iaddc-i+1)*sngl(dunif01(0,idummy)))
	    ix = indx(irand)
	    indx(irand) = indx(i)
	    indx(i) = ix
 90	 continue
      endif
      return
      end

      subroutine daddcol(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,
     *	   e,ner,indx,w,eps)
C
      integer iaddc,idelc,indx(1),nact,ner,nvars,nzzr
      double precision dd(1),e(ner,1),rr(1)
      double precision w(1),zz(nzzr,1)
      double precision eps
C
C     ***************
C     cl1 version of addcol.
C
C     this routine administers the adjustment of the
C     z*d*r   decomposition for any new zero residuals.
C     the data corresponding to the zero residuals is indexed
C     in  indx(nact+1),...,indx(iaddc).
C     ***************
C
C     +++++++++++++++
C     blas  dasum
C
C     eps  is the smallest positive number which
C     satisfies	  (1.0 + eps) .gt. 1.0	 in the
C     precision of the arithmetic being used.
C     (alternatively, for less strict zero checking,
C     eps  can be set to a user-specified tolerance.)
C     +++++++++++++++
C
      integer i,istrt,ix,nactp1,topx
      logical fail
      double precision colnrm,prjnrm
C
      double precision dasum
C
C
C     /////////////////	 begin program	//////////////////
C
      topx = nvars+1
      istrt = nact+1
      if(istrt.le.iaddc) then
C
C     ***************
C     candidates for addition to the  z*d*r  factorization
C     are inspected in random order to hinder cycling.
C     the randomization was carried out by  resid.
C
C     if a candidate has just been released from the factorization
C     or is dependent upon the columns in the factorization,
C     then it is omitted from addition.
C
C     upon exit, indices of such omitted columns are to be found in
C     indx(nact+1),...,indx(iaddc) .
C     ***************
C
	 do 10 i = istrt,iaddc
	    nactp1 = nact+1
	    ix = indx(i)
	    call dzdrpoc(nvars,nact,zz,nzzr,dd,e(1,ix),w,fail)
	    colnrm = dasum(nvars,e(1,ix),1)
	    prjnrm = dasum(nvars,w,1)
	    if(prjnrm.gt.eps*colnrm .and. ix.ne.idelc)then
	       indx(i) = indx(nactp1)
	       indx(nactp1) = ix
	       call dzdrcin(nvars,nact,zz,nzzr,dd,rr,e(1,ix),fail,w)
	    endif
 10	 continue
      endif
      return
      end

      subroutine drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,
     *	   res,grd,erql1n,pen,penpar,indx,theta,grd1)
C
      integer iaddc,indx(1),nact,nallq,nrq,nl1,ncols,ner,nvars,nrql1
      double precision e(ner,1),erql1n,grd(1),pen,penpar,res(1),theta,
     +	   grd1(1)
C
C     ***************
C     crql1 version of object.
C
C     this routine administers the evaluation of the
C     penalty (objective) function given the equation
C     and constraint residuals.	 it also computes the
C     restricted gradient of the function.
C
C     columns which are not in the  z*d*r factorization
C     but which are associated with zero residuals must
C     be included in the restricted gradient with random
C     signs as an anti-cycling device.
C     the indices of these columns are to be
C     found in	indx(nact+1),...,indx(iaddc)
C     ***************
C
C     +++++++++++++++
C     system routines  dabs,dsign
C
C     blas  daxpy,dcopy
C     +++++++++++++++
C
      integer i,idummy,ix,nactp1
      double precision one,two,three,tmp,zero,wgt,wgt1
C
      double precision dunif01
C
C     +++++++++++++++
C     the following declarations are necessary
C     for portability when  dcopy  is used, as
C     it is below, to fill arrays with a single
C     value  (zero=zip	in this case).
C     +++++++++++++++
C
      double precision zip(1)
      equivalence(zero,zip(1))
C
      data one/1.0d+00/
      data two/2.0d+00/
      data three/3.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      nrql1 = nrq+nl1
      nactp1 = nact+1
      erql1n = zero
      pen = zero
      call dcopy(nvars,zip,0,grd,1)
      call dcopy(nvars,zip,0,grd1,1)

      if(nactp1.le.ncols) then
	 do 10 i = nactp1,ncols
	    ix = indx(i)
	    tmp = res(ix)
	    wgt = dsign(one,tmp)
	    wgt1 = wgt
	    if(i.le.iaddc)then
	       wgt1 = dunif01(0,idummy)
	       if(wgt1 .lt. one/three)then
		  wgt = -one
	       else
		  if(wgt1 .gt. two/three)then
		     wgt = one
		  else
		     wgt = zero
		  endif
	       endif
	    endif
	    if(ix.le.nrq)then
	       wgt = one-two*theta+wgt
	       wgt1 = wgt
	    endif
	    if(i.le.iaddc)then
	       wgt1 = zero
	    endif
	    if(ix.le.nallq .or. wgt.le.zero)then
C     why <= zero?
	       if(wgt1 .ne. zero)then
		  if(ix.le.nrql1)then
		     erql1n = erql1n+tmp*wgt
		     tmp = tmp*penpar
		  endif
		  pen = pen+tmp*wgt
	       endif
	       if(ix.le.nrql1)then
		  wgt = wgt*penpar
		  wgt1 = wgt1*penpar
	       endif
	       call daxpy(nvars,wgt,e(1,ix),1,grd,1)
	       call daxpy(nvars,wgt1,e(1,ix),1,grd1,1)
	    endif
 10	 continue
      endif
      return
      end

      subroutine drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,grd,
     *	   coef,penpar,indx,theta,eps)
C
      integer idelc,indx(1),nact,nallq,nrq,nl1,ner,nvars,nrql1
      double precision coef(1),e(ner,1),grd(1),penpar,theta
C
C     ***************
C     crql1  version.
C
C     set up the right-hand-side vector
C     (and store in the array  coef)
C     for the linear problem which determines
C     a descent direction  p  in the case where
C     the projection of the restricted gradient is zero.
C     ***************
C
C     +++++++++++++++
C     system routines  dabs,dfloat,idint,dsign
C
C     blas  daxpy
C
C     eps  is the smallest positive number which
C     satisfies	  (1.0 + eps) .gt. 1.0	 in the
C     precision of the arithmetic being used.
C     (alternatively, for less strict zero checking,
C     eps  can be set to a user-specified tolerance.)
C     +++++++++++++++
C
      integer i,idummy,irand,ix
      double precision cf,eps,one,ope,s,tmp,tmpsav,zero,two
C
      double precision dunif01
C
      data one/1.0d+00/
      data two/2.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
C     ***************
C     find the most out-of-kilter
C     coefficient.  begin inspecting
C     the coefficients at a random index
C     to hinder cycling.  set  coef
C     to zero on the fly.
C     ***************
C
      nrql1 = nrq+nl1
      ope = one+eps
      idelc = 0
      tmpsav = zero
      s = zero
      if(1.le.nact)then
	 irand = ifix(float(nact)*sngl(dunif01(0,idummy)))
	 do 10 i = 1,nact
	    irand = irand+1
	    if(irand.gt.nact)then
	       irand = 1
	    endif
	    ix = indx(irand)
	    cf = coef(irand)
	    coef(irand) = zero
	    if(ix.gt.nallq)then
	       tmp = cf+eps
	    else
	       if(ix.le.nrql1)then
		  cf = cf/penpar
	       endif
	       tmp = ope-dabs(cf)
	       if(ix.le.nrq)then
		  tmp = tmp+dsign(one,cf)*(-one+two*theta)
	       endif
	    endif
	    if(tmp.lt.tmpsav)then
C     ? what about w_nu >1
	       idelc = irand
	       s = dsign(one,cf)
	       tmpsav = tmp
	    endif
 10	 continue
C
C     ***************
C     if no coefficients are out of kilter,
C     then return.  otherwise set a
C     value in an appropriate component
C     (indicated by  idelc)  of	 coef
C     and adjust the restricted gradient
C     if necessary.
C     ***************
C
	 if(idelc.ne.0)then
	    coef(idelc) = -s
	    ix = indx(idelc)
	    if(ix.le.nallq)then
	       tmp = -s
	       if(ix.le.nrql1)then
		  tmp = tmp*penpar
	       endif
	       if(ix.le.nrq)then
		  tmp = tmp+(one-two*theta)*penpar
	       endif
	       call daxpy(nvars,tmp,e(1,ix),1,grd,1)
	    endif
	 endif
      endif
      return
      end
C-----------------------------------------------------------------------

C
C     ---------------
C     fourth level subroutines --
C     ddkheap,dunif01,dzdrcin,dzdrcou,
C     dzdrgit,dzdrgnv,dzdrpoc
C     ---------------
C
      subroutine ddkheap(make,ir,indx,aray)
C
      integer indx(1),ir
      logical make
      double precision aray(1)
C
C     ***************
C     an adaptation of d. e. knuth,s heaping
C     routines (see volume 3 of
C     the art of computer programming  ).
C     if  make	is  .true.,  the full heap building
C     process is carried out on
C     aray(1),...,aray(ir) ,
C     and the value of	ir  is unchanged.
C     if  make	is  .false.,  one step of the sorting
C     process is carried out to provide the next
C     element of  aray	in order,  and the variable
C     ir  is decreased by one.	the interruption of the
C     sorting phase is built in via the flag  once.
C     indx  is an index vector associated with
C     aray  which must be rearranged in parallel
C     with it.
C     ***************
C
      integer i,il,it,j
      logical once
      double precision t
C
C     /////////////////	 begin program	//////////////////
C
      if(ir.gt.1)then
C
C     ***************
C     test whether or not the initial
C     heap is to be built
C     ***************
C
	 il = 1
	 if(make)then
	    il = (ir/2)+1
	 endif
	 once = .false.
c REPEAT
 100	 continue
	 if(il.gt.1)then
C           ***************
C           the heap-building phase uses this branch
C           ***************
	    il = il-1
	    it = indx(il)
	    t = aray(il)
	 else
C           ***************
C           the sorting phase uses this branch
C           ***************
	    if(make .or. once)then
	       return
	    endif
	    once = .true.
	    it = indx(ir)
	    t = aray(ir)
	    indx(ir) = indx(1)
	    aray(ir) = aray(1)
	    ir = ir-1
	    if(ir.le.1)then
	       goto 200
	    endif
	 endif
C
C     ***************
C     the remaining statements are common
C     to both phases and embody the
C     heap-rectifying (sifting) section
C     ***************
C
	 j = il
c    Repeat
 50	 continue
	 i = j
	 j = 2*j
	 if(j.lt.ir)then
	    if(aray(j).gt.aray(j+1))then
	       j = j+1
	    endif
	 else
	    if(j.ne.ir)then
	       goto 60
	    endif
	 endif
	 if(t.le.aray(j))then
	    goto 60
	 endif
	 indx(i) = indx(j)
	 aray(i) = aray(j)
	 goto 50
c    end repeat
 60	 continue
	 indx(i) = it
	 aray(i) = t
	 goto 100
c end REPEAT

 200	 continue

	 indx(1) = it
	 aray(1) = t
      else
	 if(.not.make)then
	    ir = 0
	 endif
      endif
      return
      end

      double precision function dunif01(iseed,ix)
C
      integer iseed,ix,ix0
C
      data ix0/2/
C
C     +++++++++++++++
C     system routines  dfloat,mod
C     +++++++++++++++
C
C     --------------------------------------------------------------
C     --------------------------------------------------------------
C
C     *****purpose-
C     this function returns a pseudo-random number distributed
C     uniformly in the interval (0,1).
C
C     *****parameter description-
C     on input-
C
C     iseed,  if it is nonzero modulo 9973, becomes the
C	   new seed, i.e. it replaces the internally stored
C	   value of ix0.  on machines where fortran variables
C	   retain their values between calls, the internally
C	   stored value if ix0 is the value assigned to	 ix  in
C	   the previous invocation of  dunif01.	 otherwise -- and
C	   in the first call to	 dunif01 --  ix0=2.
C
C     on output-
C
C     ix is the next integer in a pseudo-random sequence of
C	   integers between  1	and  9972  and is generated from its
C	   predecessor	ix0  (i.e.  from  iseed,  if  iseed  is nonzero
C	   modulo 9973).  ix  is the value which  iseed	 should have
C	   in the next invocation of  dunif01  to get the next
C	   pseudo-random number.  the caller will often pass the
C	   same variable for  iseed  as for  ix,
C	   e.g.	 x = dunif01(ix,ix).
C
C     *****application and usage restrictions-
C     dunif01  should only be used when portability is important and a
C     course random number generator suffices.	applications requiring
C     a fine, high precison generator should use one with a much
C     larger modulus.
C
C     *****algorithm notes-
C     dunif01 will run on any machine having at least 20 bits of ac-
C     curacy for fixed-point arithmitic.  it is based on a generator
C     recommended in (3), which passes the spectral test with flying
C     colors -- see (1) and (2).
C
C     references-
C     (1) hoaglin, d.c. (1976), theoretical properties of congruential
C     random-number generators-	 an empirical view,
C     memorandum ns-340, dept. of statistics, harvard univ.
C
C     (2) knuth, d.e. (1969), the art of computer programming, vol. 2
C     (seminumerical algorithms), addison-wesley, reading, mass.
C
C     (3) smith, c.s. (1971), multiplicative pseudo-random number
C     generators with prime modulus, j. assoc. comput. mach. 18,
C     pp. 586-593.
C
C     *****general-
C
C     this subroutine was written in connection with research
C     supported by the national science foundation under grants
C     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
C
C     permission for the use of	 dunif01  in  cl1  was
C     generously given by  v. klema  and  d. hoaglin.
C
C     --------------------------------------------------------------
C     --------------------------------------------------------------
C
      if(iseed.ne.0)then
	 ix = mod(iseed,99730)
	 if(ix.ne.0)then
	    ix0 = ix
	 endif
C  ***
C  in order that all fixed-point calculations require only 20 bit
C  arithmetic, we use two calls to  mod	 to compute
C  ix0 = mod(3432*ix0, 9973).
C  ***
      endif
      ix0 = mod(52*mod(66*ix0,99730),99730)
      ix = ix0
      dunif01 = dfloat(ix0)/99730.0d+00
      return
      end

      subroutine dzdrcin(n,k,zz,nzzr,dd,rr,col,fail,w)
C
      integer k,n,nzzr
      logical fail
      double precision col(1),dd(1),rr(1),w(1),zz(nzzr,1)
C
C     ***************
C     prepared by richard bartels
C     the university of waterloo
C     computer science department
C     latest update .... 30 november, 1979.
C
C     given the factorization
C
C     zz*dd*rr
C
C     of some  n by k  matrix
C
C     (0 .le. k .lt. n)
C     (n .ge. 1),
C
C     where
C
C     (zz-transp)*(zz) = (dd-inv),
C     dd  is diagonal and nonsingular,
C     and
C     rr  has zeros below the diagonal,
C
C     and given a  (k+1)th  column
C     to be addedto the original matrix,
C     this program updates  zz,dd and rr.
C
C     the value of  k  is increased by one.
C
C     w	 is a scratch array.
C
C     use is made of routines from the library
C     of basic linear algebra subroutines (blas).
C
C     parameters...
C
C                   input/
C     name   type   output/   sub-    description
C                   scratch  scripts
C     -------------------------------------------
C     n	     int.      i	      number of rows
C
C     k	     int.     i/o	      number of columns
C
C     zz     double   i/o	2     scaled orthogonal matrix
C
C     nzzr   int.      i	      row dimension of zz
C
C     dd     double   i/o	1     diagonal scaling matrix
C      				       (diagonal elements only)
C
C     rr     double   i/o	1     right-triangular
C     					matrix in compact form.
C
C     col    double    i	1     column to be added to rr
C
C     fail   log.      o	     .true.  if	 k,n are improper
C
C     w	     double   scr       1	  workspace
C     -------------------------------------------
C
C     the  i-th	 segment of the array  rr  is  n-i+2 spaces
C     long and contains	 1  work space followed by the
C     k-i+1  elements of row  i	 followed by  n-k
C     scratch spaces.
C     ***************
C
C     +++++++++++++++
C     blas  dcopy,ddot,srotm1,srotmg1
C     +++++++++++++++
C
      integer i,j,jdel,kp1,kp2
      double precision di,one,param(5),wi,zero
C
      double precision ddot
C
C     +++++++++++++++
C     the following declarations are necessary
C     for portability when  dcopy  is used, as
C     it is below, to fill arrays with a single value
C     (one=unity  and  zero=zip	 in this case).
C     +++++++++++++++
C
      double precision unity(1),zip(1)
      equivalence(one,unity(1)),(zero,zip(1))
C
      data one/1.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      if(k.lt.0 .or. k.ge.n .or. n.gt.nzzr)then
	 fail = .true.
      else
         if(k.le.0) then
C
C     ***************
C     for the special case that the factorization was vacuous,
C     reset the arrays	zz and dd  to represent the identity.
C     ***************
C
	    k = 0
	    call dcopy(n,unity,0,dd,1)
	    call dcopy(n*n,zip,0,zz,1)
	    do 10 i = 1,n
	       zz(i,i) = one
 10	    continue
	 endif
	 kp1 = k+1
	 kp2 = k+2
C
C     ***************
C     transform the incoming column,
C     and store the result in  w.
C     ***************
C
	 do 20 i = 1,n
	    w(i) = ddot(n,zz(1,i),1,col,1)
 20	 continue
C
C     ***************
C     zero out the spike which would result from
C     storing  w  in  rr.   update  zz	and  dd.
C     ***************
C
	 if(kp2.le.n)then
	    do 30 i = kp2,n
	       di = dd(i)
	       wi = w(i)
	       call srotmg1(dd(kp1),di,w(kp1),wi,param)
	       w(i) = wi
	       dd(i) = di
	       call srotm1(n,zz(1,kp1),1, zz(1,i),1,param)
 30	    continue
	 endif
C
C     ***************
C     store the new column, which is still
C     in the array  w,	into  rr.
C     ***************
C
	 j = kp2
	 jdel = n
	 do 40 i = 1,kp1
	    rr(j) = w(i)
	    j = j+jdel
	    jdel = jdel-1
 40	 continue
	 k = kp1
	 fail = .false.
      endif
      return
      end

      subroutine dzdrcou(n,k,zz,nzzr,dd,rr,ic,fail)
C
      integer ic,k,n,nzzr
      logical fail
      double precision dd(1),rr(1),zz(nzzr,1)
C
C     ***************
C     prepared by richard bartels
C     the university of waterloo
C     computer science department
C     latest update .... 30 november, 1979.
C
C     given the factorization
C
C          zz * dd * rr
C
C     of some  n by k  matrix
C
C     (1 <= k <= n, i.e., n >= 1),
C
C     where
C
C     (zz-transp)*(zz) = (dd-inv),
C     dd  is diagonal and nonsingular,
C     and
C     rr  has zeros below the diagonal,
C
C     and given the index  ic  of a column
C     to be removed  (1 .le. ic .le. k),
C     this program updates  zz,dd and rr .
C
C     the value of  k  is decreased by one, and
C     the column ordering in  rr  is changed.
C
C     use is made of routines from the library
C     of basic linear algebra subroutines (blas).
C
C     parameters...
C
C     name   type    input/   sub-    description
C		     output  scripts
C     -------------------------------------------
C     n	     int.      i	      number of rows
C
C     k	     int.     i/o	      number of columns
C
C     zz     double   i/o    2	      scaled orthogonal matrix
C
C     nzzr   int.      i	      row dimension of zz
C
C     dd     double   i/o    1	      diagonal scaling matrix
C				      (diagonal elements only)
C
C     rr     double   i/o    1	      right-triangular matrix in compact form.
C
C     ic     int.      i	      index of column to be removed
C
C     fail   log.      o	      .true.  if  k,n,ic are improper
C     -------------------------------------------
C
C     the  i-th	 segment of the array  rr  is  n-i+2 spaces
C     long and contains	 1  work space followed by the
C     k-i+1  elements of row  i	 followed by  n-k
C     scratch spaces.
C     ***************
C
C     +++++++++++++++
C     blas  srotm1,srotmg1
C     +++++++++++++++
C
      integer i,im1,j,jend,jinc,jstrt,km1,lstrt
      double precision param(5)
c-    double precision di,rj
C
C     /////////////////	 begin program	//////////////////
C
      fail = k.lt.1 .or. k.gt.n .or. n.gt.nzzr
      if(fail) return

C     ELSE `ok' :
      if(k.le.1)then
C     ***************
C     special cases are handled first.
C     1.  k=1 and the factorization becomes null.
C     2.  ic=k and the updating is trivial.
C     ***************
	 k = 0
      else if(ic.ge.k)then
	 k = k-1
      else
	 km1 = k-1
C
C     ***************
C     general updating step.
C     the column to be deleted must be permuted
C     to the right, and subdiagonal elements
C     which result in  rr  have to be
C     transformed to zero.
C     ***************
C
	 jstrt = ic+1
	 jend = k
	 jinc = n
	 do 100 i = 1,k
C
C     ***************
C     permutation of the  i-th	row of rr.
C     ***************
C
	    do 94 j = jstrt,jend
	       rr(j) = rr(j+1)
 94	    continue
	    if(i.gt.ic)then
C
C     ***************
C     transformation of the current and last
C     rows  (i and i-1)	 of rr	as well as
C     corresponding changes to	zz and dd.
C
C     the extra variables  di  and  rj
C     are used to avoid an error message
C     from the	pfort verifier, and they
C     may be removed, if desired, so that
C     the call to  srotmg1  would be
C
C     call srotmg1(dd(im1),dd(i),rr(lstrt),rr(jstrt),param)
C
C     ***************
C
	       im1 = i-1
c-	       di = dd(i)
c-	       rj = rr(jstrt)
c-	       call srotmg1(dd(im1),di,rr(lstrt),rj,param)
               call srotmg1(dd(im1),dd(i),rr(lstrt),rr(jstrt),param)
c-	       rr(jstrt) = rj
c-	       dd(i) = di
	       call srotm1(jend-jstrt+1,rr(lstrt+1),1,
     +		    rr(jstrt+1),1,param)
	       call srotm1(n,zz(1,im1),1,zz(1,i),1,param)
	       jstrt = jstrt+1
	    endif
C
C     ***************
C     index updating
C     ***************
C
	    lstrt = jstrt
	    jstrt = jstrt+jinc
	    jend = jend+jinc
	    jinc = jinc-1
 100	 continue
	 k = km1
      endif
      return
      end

      subroutine dzdrgit(n,k,zz,nzzr,rr,gv,sol,fail,w,big,eps)
C
      integer k,n,nzzr
      logical fail
      double precision gv(1),rr(1),w(1),sol(1),zz(nzzr,1)
C
C     ***************
C     prepared by richard bartels
C     the university of waterloo
C     computer science department
C     latest update .... 30 november, 1979.
C
C     given the factorization
C
C     zz*dd*rr
C
C     of some  n by k  matrix
C
C     (1 .le. k .le. n)
C     (n .ge. 1),
C
C     where
C
C     (zz-transp)*(zz) = (dd-inv),
C     dd  is diagonal and nonsingular,
C     and
C     rr  has zeros below the diagonal,
C
C     and given an arbitrary vector  gv	 of
C     appropriate dimension, this routine finds the
C     vector  sol  satisfying the underdetermined system
C
C     (zz*dd*rr-transp.)*(sol) = (gv).
C
C     that is,
C
C     (sol) = ((zz*dd*rr)-gen.inv.-transp.)*(gv).
C
C     the array	 dd  is not needed by  dzdrgit.
C
C     use is made of routines from the library
C     of basic linear algebra subroutines (blas).
C
C     w	 is a scratch array.
C
C     parameters...
C
C                   input/
C     name   type   output/   sub-    description
C                   scratch  scripts
C     -------------------------------------------
C     n	     int.      i	      number of rows
C
C     k	     int.     i/o	      number of columns
C
C     zz     double   i/o      2      scaled orthogonal matrix
C
C     nzzr   int.      i	      row dimension of zz
C
C     rr     double   i/o      1      right-triangular matrix in compact form.
C
C     gv     double    i       1      given vector
C
C     sol    double    o       1      solution
C
C     fail   log.      o	      .true. if	 n,k are improper, or if
C     				             rr is singular
C
C     w	     double   scr      1      workspace
C     -------------------------------------------
C
C     the  i-th	 segment of the array  rr  is  n-i+2 spaces
C     long and contains	 1  work space followed by the
C     k-i+1  elements of row  i	 followed by  n-k
C     scratch spaces.
C
C     if  gv  and  sol	are dimensioned to the
C     maximum of  n  and  k , then the same
C     storage array may be used for both of
C     these vectors.
C     ***************
C
C     +++++++++++++++
C     system routines  dabs
C
C     blas  daxpy,dcopy
C
C     big  is the largest positive number
C     which can be represented in the
C     precision of the arithmetic being used.
C     +++++++++++++++
C
      integer i,j,jdel
      double precision big,wi,one,rrj,zero,eps
C
C     +++++++++++++++
C     the following declarations are necessary
C     for portability when  dcopy  is used, as
C     it is below, to fill arrays with a single value
C     (zero=zip	 in this case).
C     +++++++++++++++
C
      double precision zip(1)
      equivalence(zero,zip(1))
C
      data one/1.0d+00/
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      if(k.lt.1 .or. k.gt.n .or. n.gt.nzzr)then
	 fail = .true.
      else
C
C     ***************
C     first solve  (rr-transp.)*(w) = (gv)
C     ***************
C
	 call dcopy(k,gv,1,w,1)
	 j = 2
	 jdel = n+1
	 do 400 i = 1,k
	    rrj = rr(j)
	    wi = w(i)
C* Here the check for ill-condition is NOT changed to use eps instead of 1/big
C     if (dabs(rrj)<one)
C     if (eps*dabs(wi)>=dabs(rrj))
	    if(dabs(rrj).lt.one)then
	       if(dabs(wi).ge.dabs(rrj)*big)then
                  fail = .true.
                  return
	       endif
	    endif
	    w(i) = wi/rrj
	    if(i.lt.k)then
	       call daxpy(k-i,(-w(i)),rr(j+1),1,w(i+1),1)
	    endif
	    j = j+jdel
	    jdel = jdel-1
C
C     ***************
C     now  (sol) = (zz)*(w)
C     ***************
C
 400	 continue

	 call dcopy(n,zip,0,sol,1)
	 do 408 i = 1,k
	    call daxpy(n,w(i),zz(1,i),1,sol,1)
 408	 continue
	 fail = .false.
      endif
      return
      end

      subroutine dzdrgnv(n,k,zz,nzzr,rr,gv,sol,fail,big)
C
      integer k,n,nzzr
      logical fail
      double precision gv(1),rr(1),sol(1),zz(nzzr,1)
C
C     ***************
C     prepared by richard bartels
C     the university of waterloo
C     computer science department
C     latest update .... 30 november, 1979.
C
C     given the factorization
C
C     zz*dd*rr
C
C     of some  n by k  matrix
C
C     (1 .le. k .le. n)
C     (n .ge. 1),
C
C     where
C
C     (zz-transp)*(zz) = (dd-inv),
C     dd  is diagonal and nonsingular,
C     and
C     rr  has zeros below the diagonal,
C
C     and given an arbitrary vector  gv	 of
C     appropriate dimension, this routine finds the
C     vector  sol  given by
C
C     (sol) = ((zz*dd*rr)-gen.inv.)*(gv),
C
C     which represents the least squares problem
C
C     (zz*dd*rr)*(sol) = (gv).
C
C     the array	 dd  is not needed by  dzdrgnv.
C
C     use is made of routines from the library
C     of basic linear algebra subroutines (blas).
C
C     parameters...
C
C     name   type   input/    sub-    description
C                   output/  scripts
C     -------------------------------------------
C     n	     int.      i	      number of rows
C
C     k	     int.     i/o	      number of columns
C
C     zz     double   i/o      2      scaled orthogonal matrix
C
C     nzzr   int.      i	      row dimension of zz
C
C     rr     double   i/o      1      right-triangular matrix in compact form.
C
C     gv     double    i       1      given vector
C
C     sol    double    o       1      solution
C
C     fail   log.      o	      .true. if	 n,k are improper, or if
C                                            rr  is singular
C     -------------------------------------------
C
C     the  i-th	 segment of the array  rr  is  n-i+2 spaces
C     long and contains	 1  work space followed by the
C     k-i+1  elements of row  i	 followed by  n-k
C     scratch spaces.
C     ***************
C
C     +++++++++++++++
C     system routines  dabs
C
C     blas  ddot
C
C     big  is the largest positive number
C     which can be represented in the
C     precision of the arithmetic being used.
C     +++++++++++++++
C
      integer i,ix,j,jdel
      double precision big,one,tden,tnum
C
      double precision ddot
C
      data one/1.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      if(k.lt.1 .or. k.gt.n .or. n.gt.nzzr)then
	 fail = .true.
      else
C
C     ***************
C     form   (v) = (zz(1)-transp)*(gv),	  where	 zz(1)
C     is the matrix of the first  k  columns of	 zz
C
C     v	 can be stored in the array  sol.
C     ***************
C
	 do 12 i = 1,k
	    sol(i) = ddot(n,zz(1,i),1,gv,1)
 12	 continue
C
C     ***************
C     backsolve the system
C     (rr)*(sol) = (v)
C     for the vector  sol
C
C     note that	 sol  and  v
C     are stored in the same array.
C     ***************
C
	 j = (((n+1)*(n+2)-(n-k+3)*(n-k+2))/2)+2
	 jdel = n-k+3
	 do 14 ix = 1,k
	    i = k-ix+1
	    tden = rr(j)
	    tnum = sol(i)
	    if(ix.gt.1)then
	       tnum = tnum-ddot(ix-1,rr(j+1),1,sol(i+1),1)
	    endif
	    if(dabs(tden).lt.one)then
	       if(dabs(tnum).ge.dabs(tden)*big)then
		  go to 160
	       endif
	    endif
	    sol(i) = tnum/tden
	    j = j-jdel
	    jdel = jdel+1
 14	 continue
	 fail = .false.
	 return
 160	 fail = .true.
      endif
      return
      end

      subroutine dzdrpoc(n,k,zz,nzzr,dd,gv,poc,fail)
C
      integer k,n,nzzr
      logical fail
      double precision dd(1),poc(1),gv(1),zz(nzzr,1)
C
C     ***************
C     prepared by richard bartels
C     the university of waterloo
C     computer science department
C     latest update .... 30 november, 1979.
C
C     zz is an	n by n	(n .ge. 1)  scaled
C     orthogonal matrix.  dd  contains the
C     diagonal elements of a diagonal scaling
C     matrix.  gv  is a given vector of length	n.
C
C     we have
C
C     (zz-transp.)*(zz) = (dd-inv.)
C
C     and
C
C     zz*dd*rr = mat
C
C     for some	n by k	(0 .le. k .le. n)
C     matrix  rr  with zeros below the diagonal
C     and some given matrix  mat.  (neither  rr
C     nor  mat	are needed by  dzdrpoc.)
C
C     then
C
C     (proj(oc)) = (zz(2))*(dd(2))*(zz(2)-transp.)
C
C     is the (orthogonal) projector on the
C     complement of the range space of	mat,
C     where  zz(2)  represents the last	 n-k
C     columns of  zz  and  dd(2)  represents the
C     lower-right-hand	n-k  order submatrix of	 dd.
C
C     dzdrpoc  produces the vector
C
C     poc = (proj(oc))*gv .
C
C     use is made of routines from the library
C     of basic linear algebra subroutines (blas).
C
C     parameters...
C
C                    input/
C     name   type   output/   sub-    description
C                   scratch  scripts
C     -------------------------------------------
C     n	     int.      i	      order of	zz,dd
C
C     k	     int.     i/o	      number of columns of  zz  defining
C     				      range of	mat
C
C     zz     double   i/o      2      scaled orthogonal matrix
C
C     nzzr   int.      i	      row dimension of zz
C
C     dd     double   i/o      1      diagonal scaling matrix
C     				      (diagonal elements only)
C
C     gv     double    i       1      vector to be projected
C
C     poc    double    o       1      projection
C
C     fail   log.      o	      .true.  if  n,k are improper
C
C     -------------------------------------------
C
C     ***************
C
C     +++++++++++++++
C     blas  daxpy,dcopy,ddot
C     +++++++++++++++
C
      integer i,kp1
      double precision wi,zero
C
      double precision ddot
C
C     +++++++++++++++
C     the following declarations are necessary
C     for portability when  dcopy  is used, as
C     it is below, to fill arrays with a single value
C     (zero=zip	 in this case).
C     +++++++++++++++
C
      double precision zip(1)
      equivalence(zero,zip(1))
C
      data zero/0.0d+00/
C
C     /////////////////	 begin program	//////////////////
C
      kp1 = k+1
      fail = .false.
c            ------- in all but the following :
      if(k.lt.0 .or. k.gt.n .or. n.lt.1 .or. n.gt.nzzr)then
	 fail = .true.

      else if(k.le.0)then
C
C     ***************
C     case 1 ... zz(2)=zz  (k=0)
C     ***************
C
	 call dcopy(n,gv,1,poc,1)

      else if(k.ge.n)then
C
C     ***************
C     case 2 ... zz(2) is vacuous  (k=n)
C     ***************
C
	 call dcopy(n,zip,0,poc,1)

      else
C
C     ***************
C     case 3 ... zz(2)	is intermediate
C     between the other two cases
C     (0 .lt. k .lt. n)
C     ***************
C
	 call dcopy(n,zip,0,poc,1)
	 do 20 i = kp1,n
	    wi = ddot(n,zz(1,i),1,gv,1)*dd(i)
	    call daxpy(n,wi,zz(1,i),1,poc,1)
 20	 continue
      endif
      return
      end

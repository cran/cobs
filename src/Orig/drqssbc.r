subroutine drqssbc(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,psw,e,ner,x,f,erql1n,res,indx,w,nt,nsol,sol,tl,toler,big,eps,icyc,tmin,k,k0,lstart,factor)
#
# This is a modification of Bartels and Conn (1980) as described in 
# Koenker and Ng (1997), "A Remark on Bartels and Conn's Linearly Constrained
# L1 Algorithm", ACM Transaction on Mathematical Software, forthcoming.
#
# It also contains the parametric linear programming on `tau' and `lambda' as 
# described in Ng (1996), "An Algorithm for Quantile Smoothing Splines",
# Computational Statistics & Data Analysis, 22, 99-118.
#
#
#     ***************
#     front end interface
#     ***************
#
#     +++++ parameters +++++
#     ----------------------------------------------------------
#                           input
#     name   type  subscrpt  output        description
#                           scratch
#     ..........................................................
#     nrq    int.    none      in      number of observations
#                                      in the rq norm that correspond
#                                      to the fidelity
#                                      (may be zero)
#
#     nl1    int.    none      in      number of observations
#                                      in the l1 norm that correspond
#                                      to the roughness measure
#                                      (may be zero)
#
#     neqc   int.    none      in      number of equality
#                                      constraints
#                                      (may be zero)
#
#     niqc   int.    none      in      number of inequality
#                                      constraints
#                                      (may be zero)
#
#     niqc1  int.    none      in      part of niqc that belongs to
#                                      the loo roughness measure
#                                      (may be zero)
#
#     nvars  int.    none      in      number of variables
#
#     nact   int.    none      out     number of active
#                                      equations/constraints
#                                      at termination
#                                      (if any, their associated
#                                      column positions in  e  will
#                                      be listed in  indx(1)
#                                      through  indx(nact) )
#
#     ifl    int.    none      out     termination code
#                                      (see below)
#
#     mxs    int.    none      in      maximum number of steps
#                                      allowed
#
#     psw    logic.  none      in      print switch
#                                      (see below)
#
#     e      real     2        in      equation/constraint matrix
#                                      the first  nrq+nl1  columns
#                                      (see note below) specify
#                                      equations, the remaining
#                                      columns (if any) specify
#                                      constraints.
#
#     ner    int.    none      in      row dimension of e
#
#     x      real     1        in      starting values for the
#                                      unknowns (use zeros if no
#                                      guess is available)
#                              out     termination values for
#                                      the unknowns
#
#     f      real     1        in      equation/constraint
#                                      right-hand sides
#
#     erql1n real    none     out      rq-l1 norm of equation
#                                      residuals at termination
#
#     res    real     1        out     equation/constraint
#                                      residuals at termination
#
#     indx   int.     1        out     index vector used to record
#                                      the order in which the columns
#                                      of  e  are being processed
#
#     w      real     1        scr.    working storage
#     nt     int.     none     out     number of unique tau or lambda
#                                      solutions while performing parametric
#                                      in tau or lambda
#     nsol   int.     none     in      upper limit for the number of unique
#                                      tau or lambda solutions
#     sol    real     2        out     matrix of solutions when performing
#                                      parametric programming in tau or lambda
#     tl     real     1        in      values of initial tau and lambda
#     toler  real     none     in      tolerance used in parametric programming
#     big    real     none     in      largest representable floating point
#                                      number
#     eps    real     none     in      least positive number satisfying  
#                                      (1.0 + eps) .gt. 1.0
#     icyc   int.     none     out     number of cycles to achieve convergence
#     tmin   real     none     in      smallest value of tau to begin 
#                                      parametric programming in tau
#     k      int.     none     out     effective dimension of the model
#     k0     int.     none     in      the largest effective dimension of the
#                                      model allowed during parametric
#                                      programming in lambda
#     lstart real    none     in       largest value of lambda to begin 
#                                      parametric programming in lambda
#     factor real    none     in       factor to determine the how big a step
#                                      to take to the next smaller lambda 
#                                      during parametric programming in lambda
#     ----------------------------------------------------------
#
#     +++++ purpose +++++
#     ----------------------------------------------------------
#     this subroutine solves the   nrq+nl1 by nvars
#     system of equations
#
#                       (a-transpose) * x   ==   b
#
#     subject to the  neqc   constraints
#
#                       (g-transpose) * x  .eq.  h
#
#     and the  niqc  inequality constraints
#
#                       (c-transpose) * x  .ge.  d
#
#     for the unknowns  x(1),...,x(nvars).
#
#     the problem must be well-posed, nontrivial
#     and overdetermined in the sense that
#
#                          nvars .ge. 1
#                          nrq+nl1 .ge. 0
#                          neqc  .ge. 0
#                          niqc  .ge. 0
#               nrq+nl1+neqc+niqc  .ge. nvars.
#
#     further, no column of  a, g  or  c  should be zero.
#     if these conditions are not met, the program
#     will terminate without performing any substantive
#     computations.
#
#     a point  x  is a solution if it minimizes the equation
#     residuals from among all points which satisfy the
#     constraints.  at any (nondegenerate) solution
#     there will be  nact  equations and constraints
#     whose residuals
#
#          (a(i)-transpose) * x - b(i)
#
#          (g(i)-transpose) * x - h(i)
#
#     and
#
#          (c(i)-transpose) * x - d(i)
#
#     are zero.
#
#     the columns of  (a,g,c)  corresponding to the zero residuals
#     are referred to as  active columns  throughout this listing.
#     the numbers of the active columns are maintained as the
#     entries  1,...,nact  of the array  indx.
#
#     a solution  x  is found by minimizing a piecewise
#     linear penalty function formed from the  l1
#     norm of the equation residuals and the sum of the
#     infeasibilities in the constraints.
#     the minimization proceeds in a step-by-step
#     fashion, terminating after a finite number of steps.
#
#     note that  a, g  and  c  appear transposed in the
#     problem formulation.  hence it is the columns of  (a,g,c)
#     which define the equations and constraints respectively.
#
#     the array  e  is a composite of   a, g and c
#     and  f  is a composite of  b, h  and  d.
#     e  should contain  a  as its first  nrq+nl1  columns.
#     it should contain  g  as its next  neqc  columns and
#     contain  c  as its remaining  niqc  columns.
#     similarly  f  should contain  b  as its first
#     nrq+nl1  components,  h  as its next  neqc  components
#     and  d  as its last  niqc  components.
#     ----------------------------------------------------------
#
#     +++++ arrays +++++
#     ----------------------------------------------------------
#     e  is to be dimensioned at least    n  by  m,
#     x                       at least    n,
#     f                       at least    m,
#     res                     at least    m,
#     indx                    at least    m,
#     w                       at least    ((3*n*n+11*n+2)/2) + (2*m).
#
#                                         where  n = nvars  and
#                                         m = nrq+nl1+neqc+niqc
#     ----------------------------------------------------------
#
#     +++++ initialization +++++
#     ----------------------------------------------------------
#     the user must initialize
#
#          nrq,nl1,neqc,niqc,nvars,mxs,psw,e,ner,x,f .
#
#     the following are set 
#     and do not require initialization
#
#          nact,indx,res .
#
#     the array  w  is used as scratch space.
#     ----------------------------------------------------------
#
#     +++++ termination codes and intermediate printing +++++
#     ----------------------------------------------------------
#     mxs  sets a limit on the number of minimization steps to be
#     taken.
#
#     upon termination  ifl  will be set according to
#     the following code ...
#
#             ifl = 1 .... successful termination.
#
#             ifl = 2 .... unsuccessful termination.
#                          constraints cannot be satisfied.
#                          problem is infeasible.
#
#             ifl = 3 .... limit imposed by  mxs  reached
#                          without finding a solution.
#
#             ifl = 4 .... program aborted.
#                          numerical difficulties
#                          due to ill-conditioning.
#
#             ifl = 5 .... nrq, nl1, nvars, neqc and/or
#                          niqc  have improper values
#                          or  e  contains a zero column.
#
#     in all cases the output parameters  x,erql1n and res
#     will contain the values which they reached at termination.
#
#     intermediate printing will be turned off if  psw = .false. .
#     on the other hand,  details of each minimization cycle
#     will be printed if  psw  is set to  .true.
#
integer ifl,indx(1),mxs,nact,ner,nt
integer nrq,nl1,neqc,niqc,niqc1,nvars,nsol,icyc,k,k0
logical psw
double precision e(ner,1),erql1n(1),f(1),res(1),w(1),x(1)
double precision sol(nvars+6,nsol),tl(1),toler,big,eps
double precision tmin,lstart,factor
#
#     ////////////////  begin program  /////////////////////////
#
call dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,psw,e,ner,x,f,
erql1n,res,indx,w,nt,nsol,sol,tl(1),tl(2),toler,big,eps,icyc,tmin,k,k0,lstart,factor)
return
end



subroutine dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,psw,e,ner,x,f,erql1n,res,
indx,w,nt,nsol,sol,t,lam,toler,big,eps,icyc,tmin,k,k0,lstart,factor)
#
#
#     ***************
#     main body
#     ***************
#
integer ifl,indx(1),mxs,nact,neqc,nrq,nl1,ner,niqc,niqc1,nvars,nrql1,nt
integer nsol,k,k0
logical psw,itend,ilend,ilfix
double precision e(ner,1),erql1n,f(1),res(1),w(1),x(1)
double precision eps,tmin,tmax,sol(nvars+6,nsol),t,lam,zero,one,toler,big,l0,l1
double precision tnxt,lnxt,lstart,factor
#
#
integer ddx,grdx,grd1x,icyc,iaddc,idelc
integer px,ptex,rrx,topx,zzx
double precision alpha,amag,cgmag,pen,penpar,told
#
data zero/0.d00/
data one/1.d00/
#
#     ////////////////  begin program  /////////////////////////
#
ifl = 0  #initialize ifl to 0
nrql1 = nrq+nl1
itend = .true.
ilend = .true.
ilfix = .true.
nt=1
tnxt=t
sol(1,nt)=t
lnxt=lam
sol(2,nt)=lam
if (t<zero||t>one){
# Note here that tmin is passed into the subroutine
	tmax = one - toler
	itend = .false.
	tnxt=tmin
	told=zero
	sol(1,nt)=tmin
	}
if (lam < zero) {
	l0 = toler
	l1 = (big-toler)
	ilend = .false.
	ilfix = .false.
	lnxt=lstart
	told=t
	sol(2,nt)=lstart
	}
if(!itend&!ilend){ #don't allow both t and lam to vary
	ifl = 7
	return
	}
#
#Note: penpar is assigned outside the loop as well as drql1sup
#
penpar=one/lnxt**.5
repeat {
	call drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,px,ptex,rrx,topx,zzx,icyc,ifl,e,ner,amag,cgmag,penpar,lnxt)
	call dnewpen(iaddc,idelc,nact,nrql1,neqc,niqc,nvars,ifl,e,ner,x,f,res,w(ptex),alpha,penpar,indx)
	repeat {
		call drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,ifl,mxs,e,ner,x,f,res,w(grdx),erql1n,pen,penpar,indx,w(zzx),
		  nvars,w(ddx),w(rrx),w(topx),tnxt,eps,w(grd1x))
#		call dmonit(nact,neqc,niqc,nvars,icyc,psw,x,alpha,erql1n,pen,penpar,indx)
		call drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,f,res,w(grdx),w(px),erql1n,amag,cgmag,penpar,indx,w(zzx),nvars,
		  w(ddx),w(rrx),w(topx),tnxt,big,eps)
		call drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,res,w(grdx),w(px),w(ptex),alpha,penpar,indx,w(topx),tnxt,big,eps,w(grd1x),idelc)
		}
		until(ifl!=0)
	if (!(itend&&ilend) && (ifl !=2 || cgmag+penpar*amag == cgmag)) {
		call drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,nact,nvars,
			nrq,nl1,neqc,niqc,niqc1,w(zzx),nvars,w(rrx),w(grdx),w(px),
			w(topx),w(topx+nvars),ifl,idelc,iaddc,icyc,alpha,amag,
			cgmag,psw,penpar,nsol,sol,x,ilfix,l0,l1,tnxt,lnxt,toler,
			erql1n,eps,big,told,k0,factor)
#update penpar to the sqrt of next lambda
		penpar=one/lnxt**.5 
		if (ifl != 0)
			break 1
		}
	else if (ifl !=2 || cgmag+penpar*amag == cgmag){ 
		break 1
		}
	}
#
#compute the effective dimensionality
k = 0
do i=1,nact
	if (indx(i) <= nrq || (indx(i) > nrq+nl1 && indx(i) <= nrq+nl1+neqc) || indx(i) > nrq+nl1+neqc+niqc1)
		k = k+1
return
end



subroutine drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,nact,nvars,nrq,nl1,
neqc,niqc,niqc1,zz,nzzr,rr,a,aa,b,bb,ifl,idelc,iaddc,icyc,alpha,amag,cgmag,psw,penpar,
nsol,sol,x,ilfix,l0,l1,tnxt,lnxt,toler,erql1n,eps,big,told,k0,factor)
#
#
#     ***************
#     perform parametric programming in "lambda" and "tau"
#     ***************
#
logical fail,psw,ilfix
integer nt,nact,nvars,nrq,nl1,ner,nzzr,ifl,k,nactp1
integer nrql1,neqc,niqc,niqc1,nallq,nalqp1,ncols,nqnp1,ix
integer idelc,iaddc,icyc,indx(1),isave,nsol,k0
double precision tmin,tmax,sgn,one,res(1),zero,e(ner,1),zz(nzzr,1)
double precision rr(1),a(1),aa(1),b(1),bb(1),thet,tnxt,tmp,sol(nvars+6,nsol)
double precision eps,two,penpar,x(1),lnxt,lamb,l0,l1,big,toler,erql1n
double precision test,prod,f(1),amag,cgmag,fidel,penal,wgt,told,factor
#
#
#     ////////////////  begin program  /////////////////////////
#
data one/1.d00/
data zero/0.d00/
data two/2.d00/
#
nrql1=nrq+nl1
nallq=nrql1+neqc
nalqp1=nallq+1
ncols=nallq+niqc
nqnp1=nrql1+1
nactp1=nact+1
thet = tnxt
if (ilfix)
	tnxt = one+eps
lamb = lnxt
if (!ilfix)
	lnxt = l0
if (ifl==1|ifl==3) {

	call scopy1(nvars,zero,0,a,1)
	call scopy1(nvars,zero,0,b,1)
	if (nacpt<=ncols)
		do i = nactp1,ncols {
			ix = indx(i)
			sgn = dsign(one,res(ix))
			test = dabs(f(ix))
			do j = 1,nvars {
				prod = dabs(e(j,ix)*x(j))
				if (prod > test)
					test = prod
				}
			test = eps*dsqrt(dfloat(nvars))*test
			if (dabs(res(ix)) < test)
				sgn = zero
			if (ilfix){
				if (ix<=nrq){
					call saxpy1(nvars,(one+sgn),e(1,ix),1,a,1)
					call saxpy1(nvars,-two,e(1,ix),1,b,1)
					}
				else if (ix<=nallq||sgn<=zero){
					call saxpy1(nvars,sgn,e(1,ix),1,a,1)
					}
				}
			else {
				if (ix<=nrq){
					call saxpy1(nvars,(one-two*thet+sgn),e(1,ix),1,a,1)
					}
				else if (ix<=nrql1){
					call saxpy1(nvars,sgn/lamb,e(1,ix),1,b,1)
					}
					else if (ix<=allq||sgn<=zero){
						call saxpy1(nvars,sgn,e(1,ix),1,a,1)
						}
				}
			}
	call dzdrgnv(nvars,nact,zz,nzzr,rr,a,aa,fail,big)
	if (fail)
		ifl = 4
	else {
		call dzdrgnv(nvars,nact,zz,nzzr,rr,b,bb,fail,big)
		if (fail)
			ifl = 4
		else {
			do i = 1,nact {
				ix = indx(i)
#a check for small bb(i) is implemented to avoid floating point overflow
				test = dabs(f(ix))
				do j = 1,nvars {
					prod = dabs(e(j,ix)*x(j))
					if (prod > test)
						test = prod
					}
				test = eps*dsqrt(dfloat(nvars))*test
				if (ix<=nrq){
					if (ilfix){
						tmp = (two+aa(i))/(two-bb(i))
						if (tmp<tnxt&&tmp>=thet) {
							tnxt = tmp
							isave = i
							}
						else {
							tmp = aa(i)/(two-bb(i))
							if (tmp<tnxt&&tmp>=thet) {
								tnxt = tmp
								}
							}
						}
					else {
						tmp = (two*thet - aa(i))/bb(i)
						if(dabs(bb(i))<test) #avoid bb near zero
							tmp = big
						if (tmp>lnxt&&tmp<lamb) {
							lnxt = tmp
							isave = i
							}
						else {
							tmp = (two*thet-two-aa(i))/bb(i)
							if(dabs(bb(i))<test) #avoid bb near zero
								tmp = big
							if(tmp>lnxt&&tmp<lamb){
								lnxt = tmp
								isave = i
								}
							}
						}
					}
				else if (ix<=nallq){
					if (ilfix){
						tmp = (one-aa(i))/bb(i)
						if(dabs(bb(i))<test) #avoid bb near zero
							tmp = big
						if (tmp<tnxt&&tmp>=thet){
							tnxt = tmp
							isave = i
							}
						else {
							tmp = -(aa(i)+one)/bb(i)
							if(dabs(bb(i))<test) #avoid bb near zero
								tmp = big
							if (tmp<tnxt&&tmp>=thet){
								tnxt = tmp
								isave = i
								}
							}
						}
					else {
						tmp = -aa(i)*lamb/(bb(i)*lamb+one)
						if (tmp>lnxt&&tmp<lamb){
							lnxt = tmp
							isave = i
							}
						else {
							tmp = -aa(i)*lamb/(bb(i)*lamb-one)
							if (tmp>lnxt&&tmp<lamb){
								lnxt = tmp
								isave = i
								}
							}
						}
					}
				else {
					if (ilfix){
						tmp = -aa(i)/bb(i)
						if(dabs(bb(i))<test) #avoid bb near zero
							tmp = big
						if (tmp<tnxt&&tmp>=thet){
							tnxt = tmp
							isave = i
							}
						}
					else {
						tmp = -aa(i)/bb(i)
						if(dabs(bb(i))<test) #avoid bb near zero
							tmp = big
						if (tmp>lnxt&&tmp<lamb){
							lnxt = tmp
							isave = i
							}
						}
					}
				}
			}
		}
	}
#
#compute the effective dimensionalty, fidelity and penalty
#
k = 0
fidel = zero
penal = zero
do i=1,nact
	if (indx(i) <= nrq || (indx(i) > nrq+nl1 && indx(i) <= nrq+nl1+neqc) || indx(i) > nrq+nl1+neqc+niqc1)
		k = k+1
#set the lower stopping criterion for lambda to be either k>=k0 or when
#lnxt < 0
tmp = lnxt-lnxt*10.0d0**(factor-4.0d0)
if(k>=k0|tmp<zero)  
	l0 = lnxt+eps
do i=nactp1,ncols {
	ix = indx(i)
	tmp = res(ix)
	wgt = dsign(one,tmp)
	if (ix<=nrq){
		wgt = one-two*told+wgt
		fidel = fidel+wgt*tmp
		}
	else if (ix<=nrql1)
			penal = penal+dabs(tmp)
	}
if (ilfix){
	if ((ifl==1|ifl==3) && tnxt<tmax){
		nt = nt+1
		sol(1,nt) = tnxt
		sol(2,nt) = lamb
		sol(3,nt-1)=dble(ifl)
		sol(4,nt-1) = fidel
		sol(5,nt-1) = penal/lamb
		sol(6,nt-1) = k
		told = tnxt
		if (indx(isave) <= nrql1)
			tnxt = tnxt+10.0d0**(factor-7.0d0)*amag
		else
			tnxt = tnxt+10.0d0**(factor-7.0d0)*cgmag
		call scopy1(nvars,x,1,sol(7,nt-1),1)
		idelc = 0
		iaddc = nact
		icyc = -1
		ifl = 0
		alpha = zero
		}
	else{
		nt = nt+1
		if  (nt>=nsol)
			ifl = 6
		sol(1,nt) = tnxt
		sol(2,nt) = lamb
		sol(3,nt-1)=dble(ifl)
		sol(4,nt-1) = fidel
		sol(5,nt-1) = penal/lamb
		sol(6,nt-1) = k
		if ((ifl==1|ifl==3) && tnxt >=tmax){
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
			call scopy1(nvars,x,1,sol(7,nt),1)
			}
		call scopy1(nvars,x,1,sol(7,nt-1),1)
		}
	}
else {
	if ((ifl==1|ifl==3) && lnxt>l0){
		nt = nt+1
		sol(1,nt) = thet
		sol(2,nt) = lnxt
		sol(3,nt-1)=dble(ifl)
		sol(4,nt-1) = fidel
		sol(5,nt-1) = penal/lamb
		sol(6,nt-1) = k
		lnxt = lnxt-lnxt*10.0d0**(factor-4.0d0)
		call scopy1(nvars,x,1,sol(7,nt-1),1)
		idelc = 0
		iaddc = nact
		icyc = -1
		ifl = 0
		alpha = zero
		}
	else{
		nt = nt+1
		if  (nt>=nsol)
			ifl = 6
		sol(1,nt) = thet
		sol(2,nt) = lnxt
		sol(3,nt-1)=dble(ifl)
		sol(4,nt-1) = fidel
		sol(5,nt-1) = penal/lamb
		sol(6,nt-1) = k
		if ((ifl==1|ifl==3) && lnxt <=l0){
			sol(1,nt) = thet
			sol(2,nt) = l0
			sol(3,nt) = sol(3,nt-1)
			sol(4,nt) = sol(4,nt-1)
			sol(5,nt) = sol(5,nt-1)
			sol(6,nt) = sol(6,nt-1)
			call scopy1(nvars,x,1,sol(7,nt),1)
			}
		call scopy1(nvars,x,1,sol(7,nt-1),1)
		}
	}
# remove lambda from e for next iteration
if(nrql1>=nrq+1)
	do i=nrq+1,nrql1
		do j=1,ner
			e(j,i) = e(j,i)/lamb
return
end





subroutine drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,px,ptex,rrx,topx,zzx,icyc,ifl,e,ner,amag,cgmag,penpar,lam)
#
integer ddx,grdx,icyc,ifl,neqc,nrql1,ner,nrq,nl1,grd1x
integer niqc,nvars,px,ptex,rrx,topx,zzx
double precision amag,cgmag,e(ner,1),penpar,lam
#
#     ***************
#     crql1  version.
#
#     set up the program
#     parameters and indices.
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#     +++++++++++++++
#
integer i,j,ncols,nqnp1
double precision oct,tmp,zero
#
data oct/8.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#     ***************
#     check validity of problem dimensions
#     ***************
#
nrql1=nrq+nl1
ncols = nrql1+neqc+niqc
if (nvars<1||neqc<0||niqc<0||nrql1<0||ncols<nvars||ner<nvars)
	ifl = 5
else {
#
#     ***************
#     set up indices for the temporary storage vector  w.
#     ***************
#
	nqnp1 = nrql1+1
	grdx = 1
	grd1x = grdx+nvars
	px = grd1x+nvars
	ptex = px+nvars
	ddx = ptex+ncols
	rrx = ddx+nvars
	zzx = rrx+(((nvars+1)*(nvars+2))/2)
	topx = zzx+nvars*nvars
#
#     ***************
#     update e with lambda only if ifl!=2, i.e. update only for the new lambda
#     ***************
#
	if( ifl!=2)
		do i=nrq+1,nrql1
			do j=1,ner
				e(j,i)=e(j,i)*lam
#
#     ***************
#     amag  is a rough estimate of the norm of  a.
#     cgmag  is a rough estimate of the norm of  (g,c).
#     together they are used to determine when the
#     penalty parameter is too small and when the
#     restricted gradient is zero.
#     ***************
#
	amag = zero
	cgmag = zero
	if (1<=nrql1) {
		do j = 1,nrql1 {
			tmp = zero
			do i = 1,nvars
				tmp = tmp+dabs(e(i,j))
			if (tmp<=zero)  
				go to 10
			if (tmp>amag)
				amag = tmp
			}
		go to 20
		10  ifl = 5
		return
		}
	20  if (nqnp1<=ncols) {
		do j = nqnp1,ncols {
			tmp = zero
			do i = 1,nvars
				tmp = tmp+dabs(e(i,j))
			if (tmp<=zero)
				go to 30
			if (tmp>cgmag)
				cgmag = tmp
			}
		go to 40
		30  ifl = 5
		return
		}
#
#     ***************
#     initialize  ifl,icyc,penpar
#     ***************
#
	40  ifl = 2
	icyc = -1
	}
return
end



subroutine dnewpen(iaddc,idelc,nact,neqns,neqc,niqc,nvars,ifl,e,ner,x,f,res,pte,alpha,penpar,indx)
#
integer iaddc,idelc,ifl,indx(1),nact
integer neqc,neqns,ner,niqc,nvars
double precision alpha,e(ner,1),f(1),penpar,pte(1),res(1),x(1)
#
#     ***************
#     cl1  version.
#
#     begin a round of minimization steps
#     with a new penalty parameter value.
#     ***************
#
#     +++++++++++++++
#     blas  sdot1
#     +++++++++++++++
#
integer i,ncols
double precision oct,one,zero
#
double precision sdot1
#
data zero/0.0d+00/
data one/1.0d+00/
data oct/8.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#     ***************
#     set penalty parameter value.
#     erase record of active equation/constraints.
#     ***************
#
if (ifl==2) {
	ncols = neqns+neqc+niqc
	ifl = 0
	nact = 0
	iaddc = 0
	idelc = 0
	alpha = zero
	penpar = penpar/oct
#
#     ***************
#     initialize  indx,res,pte,indx
#     ***************
#
	do i = 1,ncols {
		res(i) = sdot1(nvars,e(1,i),1,x,1)-f(i)
		pte(i) = zero
		indx(i) = i
		}
	}
return
end


subroutine drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,ifl,mxs,e,ner,x,f,res,grd,erql1n,pen,penpar,indx,zz,nzzr,dd,
  rr,w,theta,eps,grd1)
#
integer iaddc,idelc,icyc,ifl,indx(1),mxs
integer nact,neqc,nrq,nl1,ner,niqc,nvars,nzzr
double precision dd(1),e(ner,1),erql1n,f(1),grd(1),pen,penpar,theta,grd1(1)
double precision res(1),rr(1),w(1),x(1),zz(nzzr,1),eps
#
#     ***************
#     crql1  version.
#
#     preparation for next minimization step.
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#     +++++++++++++++
#
integer nallq,ncols,nrql1
double precision one,zero
#
data one/1.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#     ***************
#     determine the active equations and active
#     constraints.  compute residuals and function value.
#     update the  z*d*r  decomposition.
#     ***************
#
nrql1 = nrq+nl1
nallq = nrql1+neqc
ncols = nallq+niqc
if (ifl==0) {
	icyc = icyc+1
	if (icyc>mxs)
		ifl = 3
	else {
		call ddelcol1(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,indx)
		call dresid(iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
		call daddcol(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,e,ner,indx,w,eps)
		call drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,res,grd,erql1n,pen,penpar,indx,theta,grd1)
		}
	}
return
end


subroutine drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,f,res,grd,p,erql1n,amag,cgmag,penpar,indx,zz,nzzr,dd,rr,w,
  theta,big,eps)
#
integer idelc,ifl,indx(1),nact,neqc,nrq,nl1,ner,niqc,nvars,nzzr,nrql1
double precision amag,cgmag,dd(1),e(ner,1),erql1n,f(1),grd(1),p(1),penpar
double precision res(1),rr(1),w(1),x(1),zz(nzzr,1),theta
#
#     ***************
#     crql1  version.
#
#     determine descent direction  p
#     (or discover optimality)
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#
#     blas  sasum1,scopy1,sscal1
#
#     eps  is the smallest positive number which
#     satisfies   (1.0 + eps) .gt. 1.0   in the
#     precision of the arithmetic being used.
#     (alternatively, for less strict zero checking,
#      eps  can be set to a user-specified tolerance.)
#     +++++++++++++++
#
integer coefx,i,ix,nallq,nalqp1,ncols,nqnp1,topx,nrql1
logical fail
double precision grdnrm,one,pnrm,prod,test,zero
double precision eps,big
#
double precision sasum1
#
#     +++++++++++++++
#     the following declarations are necessary
#     for portability when  scopy1  is used, as
#     it is below, to fill arrays with a single value
#     (one=unity  and  zero=zip  in this case).
#     +++++++++++++++
#
double precision unity(1),zip(1)
equivalence(one,unity(1)),(zero,zip(1))
#
data one/1.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
nrql1 = nrq+nl1
idelc = 0
if (ifl==0) {
	nallq = nrql1+neqc
	nalqp1 = nallq+1
	ncols = nallq+niqc
	nqnp1 = nrql1+1
	coefx = 1
	topx = coefx+nvars
#
#     ***************
#     project the negative of the restricted gradient
#     onto the orthogonal complement of the space
#     spanned by the active columns.
#     ***************
#
	call dzdrpoc(nvars,nact,zz,nzzr,dd,grd,p,fail)
	if (fail){
		ifl = 4
		}
	else {
		call sscal1(nvars,-one,p,1)
		pnrm = sasum1(nvars,p,1)
		grdnrm = sasum1(nvars,grd,1)
#
#     ***************
#     if the projection is not zero,
#     it will serve as a descent direction.
#
#     otherwise find the representation of
#     the restricted gradient as a linear
#     combination of the active columns.
#     the coefficients of the linear combination
#     are to be stored in the array  coef
#     (that is, in  w(coefx),...,w(coefx+nact-1)).
#     ***************
#
		if (pnrm<=eps*(amag*penpar+cgmag)) {
			if (nact!=0) {
				call dzdrgnv(nvars,nact,zz,nzzr,rr,grd,w(coefx),fail,big)
				if (fail) {
					ifl = 4
					return
					}
				else {
#
#     ***************
#     convert the coefficients of the linear
#     combination into a descent direction  p ,
#     or determine optimality.
#
#     if the optimality test is not satisfied,
#     drql1gv  will indicate an equation/constraint
#     to be deleted from activity by the value
#     of  idelc.  for optimality,  idelc=0.
#     ***************
#
					call drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,grd,w(coefx),penpar,indx,theta,eps)
					pnrm = zero
					if (idelc!=0) {
						call dzdrgit(nvars,nact,zz,nzzr,rr,w(coefx),p,fail,w(topx),big,eps)
						if (!fail)
							pnrm = sasum1(nvars,p,1)
						if (fail) {
							ifl = 4
							return
							}
						}
#
#     ***************
#     if a descent direction  p  could have been found,
#     it has been obtained by this point in the program.
#
#     check for optimality.
#
#     pnrm  has been set exactly zero
#     after the call to subroutine  drql1gv
#     if the optimality conditions are satisfied.
#     the check below has been made somewhat
#     complicated to allow for the rare event that
#     the restricted gradient is zero and no
#     columns are active,  or that the  rq  norm of
#               (a-transpose) * x - f
#     is computationally zero.
#     (the call to the subroutine  refine
#      may be omitted, if desired.)
#     ***************
#
					if (pnrm>eps*(amag*penpar+cgmag)){
						do i = 1,nrql1 {
							test = dabs(f(i))
							do ix = 1,nvars {
								prod = dabs(e(ix,i)*x(ix))
								if (prod>test)
									test = prod
								}
							if (dabs(res(i))>eps*test){
								return
								}
							}
						}
					}
				}
			ifl = 1
#			call drql1rf(nact,nrq,nl1,ncols,nvars,ifl,e,ner,x,f,erql1n,res,indx,zz,nzzr,rr,w,theta,big,eps)
			if (ifl==1)
#
#     ***************
#     if the problem has constraints,
#     check feasibility.
#     ***************
#
				if (nqnp1<=ncols) {
					do i = nqnp1,ncols {
						test = dabs(f(i))
						do ix = 1,nvars {
							prod = dabs(e(ix,i)*x(ix))
							if (prod>test)
								test = prod
							}
# NOTE: the criterion for checking feasibility is relaxed by (eps * test)^.5
# rather than eps*test
						test = (eps*test)**.5
#						test = eps*test
						if (i>nallq) {
							if (res(i)<(-test)){
								go to 20
								}
							}
						else if (dabs(res(i))>test)
							go to 10
						}
					return
					10  ifl = 2
					return
					20  ifl = 2
					}
			}
		}
	}
return
end


subroutine drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,res,grd,p,pte,alpha,penpar,indx,alf,theta,big,eps,grd1,idelc)
#
integer iaddc,ifl,indx(1),nact,neqc,nrq,nl1,nrql1,ner,niqc,nvars,idelc
double precision alpha,alf(1),e(ner,1),grd(1),p(1),grd1(1)
double precision penpar,pte(1),res(1),x(1)
double precision theta,sgn1
#
#     ***************
#     cl1  version.
#
#     piecewise linear line search.
#     ***************
#
#     +++++++++++++++
#     system routines dabs,dsign
#
#     blas  sasum1,saxpy1,sdot1
#
#     eps  is the smallest positive number which
#     satisfies   (1.0 + eps) .gt. 1.0   in the
#     precision of the arithmetic being used.
#     (alternatively, for less strict zero checking,
#      eps  can be set to a user-specified tolerance.)
#
#     big  is the largest positive number
#     which can be represented in the
#     precision of the arithmetic being used.
#     +++++++++++++++
#
integer i,iin,ix,jx,nactp1,nallq,ncols,num,numnac
double precision big,den,eps,grdnrm,one,pnrm,ptg,ptg1
double precision ratio,resid,tmp,two,zero
#
double precision sasum1,sdot1,etp
#
data one/1.0d+00/
data two/2.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#     ***************
#     this routine determines all of the ratios  alf
#     of the form
#        -res(i)/((e(.,i)-transp)*p),
#              for  i = k+1,...,mpl
#     which are nonnegative and hence indicate distances
#     from the point  x  to breakpoints which will
#     be encountered in travel along direction  p.
#     the index vector  indx  is rearranged so that
#     its  k+1  through  num  components correspond to
#     these nonnegative ratios.
#     the results are heaped so that the  alf  values can
#     be inspected in order from smallest to largest.
#     the breakpoint  alpha  giving the minimum objective
#     function value is found, and  x  is
#     adjusted to  x + alpha*p .
#
#     the inner products  (e(.,i)-transpose)*p  are saved
#     for later use in updating the residual values.
#     ***************
#
alpha = zero
if (ifl==0) {
	nrql1 = nrq+nl1
	nallq = nrql1+neqc
	ncols = nallq+niqc
	nactp1 = nact+1
	num = 0
	if (1<=nact)
		do i = 1,nact {
			ix = indx(i)
			pte(ix) = sdot1(nvars,e(1,ix),1,p,1)
			}
#update the correct gradient
	if (nactp1<=iaddc)
		do i = nactp1,iaddc{
			ix = indx(i)
			etp = sdot1(nvars,e(1,ix),1,p,1)
			sgn1 = dsign(one,etp)
			if (ix<=nrq)
				sgn1 = one-two*theta + sgn1
			if (ix<=nallq||sgn1<=zero){
				if (ix<=nrql1){
					sgn1 = sgn1*penpar
					}
				call saxpy1(nvars,sgn1,e(1,ix),1,grd1,1)
				}
			}
	if (idelc!=0){
		ix=indx(idelc)
		etp=sdot1(nvars,e(1,ix),1,p,1)
		sgn1=dsign(one,etp)
		if (ix<=nrq)
			sgn1 = one-two*theta + sgn1
		if (ix<=nallq||sgn1<=zero){
			if(ix<=nrql1){
				sgn1=sgn1*penpar
				}
			call saxpy1(nvars,sgn1,e(1,ix),1,grd1,1)
			}
		}
	if (nactp1>ncols)
		ifl = 1
	else {
		do i = nactp1,ncols {
			ix = indx(i)
			resid = res(ix)
			den = sdot1(nvars,e(1,ix),1,p,1)
			pte(ix) = den
			if (dsign(one,resid)!=dsign(one,den)||resid==zero) {
				resid = dabs(resid)
				den = dabs(den)
				if (den<one)
					if (resid>=den*big)
						next 1
				ratio = resid/den
				num = num+1
				numnac = num+nact
				jx = indx(numnac)
				indx(numnac) = ix
				indx(i) = jx
				alf(num) = ratio
				}
			}
		if (num<=0)
			ifl = 2
		else {
#
#     ***************
#     heap the positive ratios
#     ***************
#
			call ddkheap(.true.,num,indx(nactp1),alf)
#
#     ***************
#     travel along  p  until no further decrease in the
#     penalty function is possible
#     ***************
#
			iin = num
			ptg = sdot1(nvars,grd,1,p,1)
			ptg1 = sdot1(nvars,grd1,1,p,1)
			pnrm = sasum1(nvars,p,1)
			grdnrm = sasum1(nvars,grd1,1)
			do i = 1,num {
				ix = indx(nactp1)
				if(res(ix)==zero)
					tmp = zero
				else
				tmp = -dsign(one,res(ix))
				if (ix<=nallq)
					tmp = tmp*two
				if (ix<=nrql1)
					tmp = tmp*penpar
				ptg1 = ptg1+tmp*pte(ix)
				if (ptg1>=(-eps)*grdnrm*pnrm)
					go to 140
				call ddkheap(.false.,iin,indx(nactp1),alf)
				}
			ifl = 2
			return
			140  iaddc = nactp1
#
#     ***************
#     adjust  x  to  x + alpha*p
#     ***************
#
			alpha = alf(1)
			call saxpy1(nvars,alpha,p,1,x,1)
			}
		}
	}
return
end


subroutine drql1rf(nact,nrq,nl1,ncols,nvars,ifl,e,ner,x,f,erql1n,res,indx,zz,nzzr,rr,w,theta,big,eps)
#
integer ifl,indx(1),nact,ncols,nrq,nl1,ner,nvars,nzzr,nrql1
double precision e(ner,1),erql1n,f(1),res(1),rr(1),w(1),x(1),zz(nzzr,1)
double precision theta, wgt
#
#     ***************
#     a routine for refining the solution
#     produced by  crql1.
#
#     (this routine may be omitted if desired.)
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#
#     blas  sdot1
#     +++++++++++++++
#
##integer i,ix
logical fail
double precision tmp,zero,one,two
#
double precision sdot1,big,eps
#
data zero/0.0d+00/
data one/1.0d+00/
data two/2.0d+00/
#
#     /////////////// begin program ///////////////
#
nrql1 = nrq+nl1
if (nact!=0) {
	if (fail){
		ifl = 4
		}
	else {
		erql1n = zero
		do i = 1,ncols {
			tmp = sdot1(nvars,e(1,i),1,x,1)-f(i)
			wgt = dsign(one,tmp)
			if (i<=nrq)
				wgt = one-two*theta+wgt
			res(i) = tmp
			if (i<=nrql1)
				erql1n = erql1n+wgt*tmp
			}
		}
	}
return
end


#
#     ---------------
#     third level subroutines --
#          ddelcol1,dresid,addcol,drql1obj,getv
#     ---------------
#
subroutine ddelcol1(iaddc,idelc,nact,nrow,zz,nzzr,dd,rr,indx)
#
integer indx(1),nzzr,nact,idelc,iaddc,nrow
double precision dd(1),rr(1),zz(nzzr,1)
#
#     ***************
#     cl1  version of  ddelcol1.
#
#     this routine administers the deletion of the column
#     indicated by the value of idelc
#     from an  nrow by nact   z*d*r   decomposition.
#     note that the value of idelc
#     is the number of a column in the decomposition
#     rather than a number which refers to
#     a column in the matrix  e.
#     (the  e-column  numbers corresponding to
#     the columns of the factorization are to be
#     found in   indx(1),...,indx(nact) .
#     the contents of   indx(nact+1),...,indx(iaddc)
#     indicate columns of  e  which are slated for
#     addition to the decomposition.)
#     the vector  indx   is rearranged by
#     permuting the element which corresponds to
#     the deletion out to the   iaddc-th  position.
#     nact  and  iaddc  are decreased accordingly.
#     ***************
#
integer i,idlp1,ixdlc
logical fail
#
#     /////////////////  begin program  //////////////////
#
if (idelc!=0) {
	idlp1 = idelc+1
	ixdlc = indx(idelc)
	do i = idlp1,iaddc
		indx(i-1) = indx(i)
	indx(iaddc) = ixdlc
	iaddc = iaddc-1
	call dzdrcou(nrow,nact,zz,nzzr,dd,rr,idelc,fail)
	idelc = ixdlc
	}
return
end



subroutine dresid(iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
#
integer iaddc,indx(1),nact,ncols,ner,nvars
double precision e(ner,1),f(1),res(1),x(1)
#
#
#     ***************
#     compute the residuals
#          (e(.,ix)-transp)*x - f(ix)  .
#     the residuals are stored in the array  res.
#     indx  is rearranged so that the zero residuals
#     correspond to  indx(1),...,indx(iaddc)  .
#     ***************
#
#     +++++++++++++++
#     system routines  dabs,idint,dfloat,dsqrt
#
#     blas  sdot1
#
#     eps  is the smallest positive number which
#     satisfies   (1.0 + eps) .gt. 1.0   in the
#     precision of the arithmetic being used.
#     (alternatively, for less strict zero checking,
#      eps  can be set to a user-specified tolerance.)
#     +++++++++++++++
#
integer i,iadp1,idummy,irand,ix,j,nactp1
double precision eps,prod,temp,test,tol,zero
#
double precision sdot1,dunif01
#
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#tol = eps*dsqrt(dfloat(nvars))
tol = eps
nactp1 = nact+1
if (1<=iaddc)
#
#     ***************
#     zero out all residuals known to be zero.
#     ***************
#
	do i = 1,iaddc {
		ix = indx(i)
		res(ix) = zero
		}
#
#     ***************
#     compute the remaining residuals.
#     detect any more residuals which
#     are computationally zero, and
#     set them exactly zero.  their
#     associated indices are permuted
#     so that they are stored in
#     indx(nact+1),...,nact(iaddc).
#
#     (a fairly tight zero check is used.
#     it is far less expensive in running
#     time to neglect an extra zero
#     residual than to accept it and risk
#     invoking the anti-cycling
#     mechanisms in the program.
#     the accuracy of the solution as
#     finally determined is not affected.)
#     ***************
#
iadp1 = iaddc+1
if (iadp1<=ncols)
	do i = iadp1,ncols {
		ix = indx(i)
		temp = sdot1(nvars,e(1,ix),1,x,1)-f(ix)
		test = dabs(f(ix))
		do j = 1,nvars {
			prod = dabs(e(j,ix)*x(j))
			if (prod>test)
				test = prod
			}
		test = tol*test
		if (dabs(temp)>test)
			res(ix) = temp
		else {
			iaddc = iaddc+1
			indx(i) = indx(iaddc)
			indx(iaddc) = ix
			res(ix) = zero
			}
		}
#
#     ***************
#     if any new zero residuals have
#     been found, randomize their
#     ordering as an anti-cycling
#     device for  addcol.
#     ***************
#
if (iaddc>nactp1)
	do i = nactp1,iaddc {
		irand = i+ifix(float(iaddc-i+1)*sngl(dunif01(0,idummy)))
		ix = indx(irand)
		indx(irand) = indx(i)
		indx(i) = ix
		}
return
end



subroutine daddcol(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,e,ner,indx,w,eps)
#
integer iaddc,idelc,indx(1),nact,ner,nvars,nzzr
double precision dd(1),e(ner,1),rr(1)
double precision w(1),zz(nzzr,1)
double precision eps
#
#     ***************
#     cl1 version of addcol.
#
#     this routine administers the adjustment of the
#     z*d*r   decomposition for any new zero residuals.
#     the data corresponding to the zero residuals is indexed
#     in  indx(nact+1),...,indx(iaddc).
#     ***************
#
#     +++++++++++++++
#     blas  sasum1
#
#     eps  is the smallest positive number which
#     satisfies   (1.0 + eps) .gt. 1.0   in the
#     precision of the arithmetic being used.
#     (alternatively, for less strict zero checking,
#      eps  can be set to a user-specified tolerance.)
#     +++++++++++++++
#
integer i,istrt,ix,nactp1,topx
logical fail
double precision colnrm,prjnrm
#
double precision sasum1
#
#
#     /////////////////  begin program  //////////////////
#
topx = nvars+1
istrt = nact+1
if (istrt<=iaddc)
#
#     ***************
#     candidates for addition to the  z*d*r
#     factorization are inspected in random
#     order to hinder cycling.
#     the randomization was carried out by  resid.
#
#     if a candidate has just been released
#     from the factorization or is dependent upon the
#     columns in the factorization,
#     then it is omitted from addition.
#
#     upon exit, indices of such omitted
#     columns are to be found in
#          indx(nact+1),...,indx(iaddc) .
#     ***************
#
	do i = istrt,iaddc {
		nactp1 = nact+1
		ix = indx(i)
		call dzdrpoc(nvars,nact,zz,nzzr,dd,e(1,ix),w,fail)
		colnrm = sasum1(nvars,e(1,ix),1)
		prjnrm = sasum1(nvars,w,1)
		if (prjnrm>eps*colnrm&&ix!=idelc) {
			indx(i) = indx(nactp1)
			indx(nactp1) = ix
			call dzdrcin(nvars,nact,zz,nzzr,dd,rr,e(1,ix),fail,w)
			}
		}
return
end



subroutine drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,res,grd,erql1n,pen,penpar,indx,theta,grd1)
#
integer iaddc,indx(1),nact,nallq,nrq,nl1,ncols,ner,nvars,nrql1
double precision e(ner,1),erql1n,grd(1),pen,penpar,res(1),theta,grd1(1)
#
#     ***************
#     crql1 version of object.
#
#     this routine administers the evaluation of the
#     penalty (objective) function given the equation
#     and constraint residuals.  it also computes the
#     restricted gradient of the function.
#
#     columns which are not in the  z*d*r factorization
#     but which are associated with zero residuals must
#     be included in the restricted gradient with random
#     signs as an anti-cycling device.
#     the indices of these columns are to be
#     found in  indx(nact+1),...,indx(iaddc)
#     ***************
#
#     +++++++++++++++
#     system routines  dabs,dsign
#
#     blas  saxpy1,scopy1
#     +++++++++++++++
#
integer i,idummy,ix,nactp1
double precision half,one,two,three,tmp,zero,wgt,wgt1
#
double precision dunif01
#
#     +++++++++++++++
#     the following declarations are necessary
#     for portability when  scopy1  is used, as
#     it is below, to fill arrays with a single
#     value  (zero=zip  in this case).
#     +++++++++++++++
#
double precision zip(1)
equivalence(zero,zip(1))
#
data half/0.5d+00/
data one/1.0d+00/
data two/2.0d+00/
data three/3.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
nrql1 = nrq+nl1
nactp1 = nact+1
erql1n = zero
pen = zero
call scopy1(nvars,zip,0,grd,1)
call scopy1(nvars,zip,0,grd1,1)
if (nactp1<=ncols)
	do i = nactp1,ncols {
		ix = indx(i)
		tmp = res(ix)
		wgt = dsign(one,tmp)
		wgt1 = wgt
		if (i<=iaddc){
			wgt1 = dunif01(0,idummy)
			if (wgt1 < one/three)
				wgt = -one
			else if (wgt1 > two/three)
				wgt = one
				else
					wgt = zero
			}
		if (ix<=nrq){
			wgt = one-two*theta+wgt
			wgt1 = wgt
			}
		if (i<=iaddc)	
			wgt1 = zero
		if (ix<=nallq||wgt<=zero) {  #why <= zero?
			if (wgt1 != zero){
				if (ix<=nrql1){
					erql1n = erql1n+tmp*wgt
					tmp = tmp*penpar
					}
				pen = pen+tmp*wgt
				}
			if (ix<=nrql1) {
				wgt = wgt*penpar
				wgt1 = wgt1*penpar
				}
			call saxpy1(nvars,wgt,e(1,ix),1,grd,1)
			call saxpy1(nvars,wgt1,e(1,ix),1,grd1,1)
			}
		}
return
end





subroutine drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,grd,coef,penpar,indx,theta,eps)
#
integer idelc,indx(1),nact,nallq,nrq,nl1,ner,nvars,nrql1
double precision coef(1),e(ner,1),grd(1),penpar,theta
#
#     ***************
#     crql1  version.
#
#     set up the right-hand-side vector
#     (and store in the array  coef)
#     for the linear problem which determines
#     a descent direction  p  in the case where
#     the projection of the restricted gradient is zero.
#     ***************
#
#     +++++++++++++++
#     system routines  dabs,dfloat,idint,dsign
#
#     blas  saxpy1
#
#     eps  is the smallest positive number which
#     satisfies   (1.0 + eps) .gt. 1.0   in the
#     precision of the arithmetic being used.
#     (alternatively, for less strict zero checking,
#      eps  can be set to a user-specified tolerance.)
#     +++++++++++++++
#
integer i,idummy,irand,ix
double precision cf,eps,one,ope,s,tmp,tmpsav,zero,two
#
double precision dunif01
#
data one/1.0d+00/
data two/2.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
#     ***************
#     find the most out-of-kilter
#     coefficient.  begin inspecting
#     the coefficients at a random index
#     to hinder cycling.  set  coef
#     to zero on the fly.
#     ***************
#
nrql1 = nrq+nl1
ope = one+eps
idelc = 0
tmpsav = zero
if (1<=nact) {
	irand = ifix(float(nact)*sngl(dunif01(0,idummy)))
	do i = 1,nact {
		irand = irand+1
		if (irand>nact)
			irand = 1
		ix = indx(irand)
		cf = coef(irand)
		coef(irand) = zero
		if (ix>nallq)
			tmp = cf+eps
		else {
			if (ix<=nrql1)
				cf = cf/penpar
			tmp = ope-dabs(cf)
			if (ix<=nrq)
				tmp = tmp+dsign(one,cf)*(-one+two*theta)
			}
		if (tmp<tmpsav) { #? what about w_nu >1
			idelc = irand
			s = dsign(one,cf)
			tmpsav = tmp
			}
		}
#
#     ***************
#     if no coefficients are out of kilter,
#     then return.  otherwise set a
#     value in an appropriate component
#     (indicated by  idelc)  of  coef
#     and adjust the restricted gradient
#     if necessary.
#     ***************
#
	if (idelc!=0) {
		coef(idelc) = -s
		ix = indx(idelc)
		if (ix<=nallq) {
			tmp = -s
			if (ix<=nrql1)
				tmp = tmp*penpar
			if (ix<=nrq)
				tmp = tmp+(one-two*theta)*penpar
			call saxpy1(nvars,tmp,e(1,ix),1,grd,1)
			}
		}
	}
return
end


#
#     ---------------
#     fourth level subroutines --
#               ddkheap,dunif01,dzdrcin,dzdrcou,
#               dzdrgit,dzdrgnv,dzdrpoc
#     ---------------
#
subroutine ddkheap(make,ir,indx,aray)
#
integer indx(1),ir
logical make
double precision aray(1)
#
#     ***************
#     an adaptation of d. e. knuth,s heaping
#     routines (see volume 3 of
#          the art of computer programming  ).
#     if  make  is  .true.,  the full heap building
#     process is carried out on
#          aray(1),...,aray(ir) ,
#     and the value of  ir  is unchanged.
#     if  make  is  .false.,  one step of the sorting
#     process is carried out to provide the next
#     element of  aray  in order,  and the variable
#     ir  is decreased by one.  the interruption of the
#     sorting phase is built in via the flag  once.
#     indx  is an index vector associated with
#     aray  which must be rearranged in parallel
#     with it.
#     ***************
#
integer i,il,it,j
logical once
double precision t
#
#     /////////////////  begin program  //////////////////
#
if (ir>1) {
#
#     ***************
#     test whether or not the initial
#     heap is to be built
#     ***************
#
	il = 1
	if (make)
		il = (ir/2)+1
	once = .false.
	repeat {
		if (il>1) {
#
#     ***************
#     the heap-building phase uses this branch
#     ***************
#
			il = il-1
			it = indx(il)
			t = aray(il)
			}
		else {
#
#     ***************
#     the sorting phase uses this branch
#     ***************
#
			if (make||once)
				return
			once = .true.
			it = indx(ir)
			t = aray(ir)
			indx(ir) = indx(1)
			aray(ir) = aray(1)
			ir = ir-1
			if (ir<=1)
				break 1
			}
#
#     ***************
#     the remaining statements are common
#     to both phases and embody the
#     heap-rectifying (sifting) section
#     ***************
#
		j = il
		repeat {
			i = j
			j = 2*j
			if (j<ir) {
				if (aray(j)>aray(j+1))
					j = j+1
				}
			else if (j!=ir)
				break 1
			if (t<=aray(j))
				break 1
			indx(i) = indx(j)
			aray(i) = aray(j)
			}
		indx(i) = it
		aray(i) = t
		}
	indx(1) = it
	aray(1) = t
	}
else if (!make)
	ir = 0
return
end



double precision function dunif01(iseed,ix)
#
integer iseed,ix,ix0
#
data ix0/2/
#
#     +++++++++++++++
#     system routines  dfloat,mod
#     +++++++++++++++
#
#     --------------------------------------------------------------
#     --------------------------------------------------------------
#
#     *****purpose-
#     this function returns a pseudo-random number distributed
#     uniformly in the interval (0,1).
#
#     *****parameter description-
#     on input-
#
#     iseed,  if it is nonzero modulo 9973, becomes the
#          new seed, i.e. it replaces the internally stored
#          value of ix0.  on machines where fortran variables
#          retain their values between calls, the internally
#          stored value if ix0 is the value assigned to  ix  in
#          the previous invocation of  dunif01.  otherwise -- and
#          in the first call to  dunif01 --  ix0=2.
#
#     on output-
#
#     ix is the next integer in a pseudo-random sequence of
#          integers between  1  and  9972  and is generated from its
#          predecessor  ix0  (i.e.  from  iseed,  if  iseed  is nonzero
#          modulo 9973).  ix  is the value which  iseed  should have
#          in the next invocation of  dunif01  to get the next
#          pseudo-random number.  the caller will often pass the
#          same variable for  iseed  as for  ix,
#          e.g.  x = dunif01(ix,ix).
#
#     *****application and usage restrictions-
#     dunif01  should only be used when portability is important and a
#     course random number generator suffices.  applications requiring
#     a fine, high precison generator should use one with a much
#     larger modulus.
#
#     *****algorithm notes-
#     dunif01 will run on any machine having at least 20 bits of ac-
#     curacy for fixed-point arithmitic.  it is based on a generator
#     recommended in (3), which passes the spectral test with flying
#     colors -- see (1) and (2).
#
#     references-
#     (1) hoaglin, d.c. (1976), theoretical properties of congruential
#     random-number generators-  an empirical view,
#     memorandum ns-340, dept. of statistics, harvard univ.
#
#     (2) knuth, d.e. (1969), the art of computer programming, vol. 2
#     (seminumerical algorithms), addison-wesley, reading, mass.
#
#     (3) smith, c.s. (1971), multiplicative pseudo-random number
#     generators with prime modulus, j. assoc. comput. mach. 18,
#     pp. 586-593.
#
#     *****general-
#
#     this subroutine was written in connection with research
#     supported by the national science foundation under grants
#     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
#
#     permission for the use of  dunif01  in  cl1  was
#     generously given by  v. klema  and  d. hoaglin.
#
#     --------------------------------------------------------------
#     --------------------------------------------------------------
#
if (iseed!=0) {
	ix = mod(iseed,99730)
	if (ix!=0)
		ix0 = ix
	}
#  ***
#  in order that all fixed-point calculations require only 20 bit
#  arithmetic, we use two calls to  mod  to compute
#  ix0 = mod(3432*ix0, 9973).
#  ***
ix0 = mod(52*mod(66*ix0,99730),99730)
ix = ix0
dunif01 = dfloat(ix0)/99730.0d+00
return
end



subroutine dzdrcin(n,k,zz,nzzr,dd,rr,col,fail,w)
#
integer k,n,nzzr
logical fail
double precision col(1),dd(1),rr(1),w(1),zz(nzzr,1)
#
#     ***************
#     prepared by richard bartels
#     the university of waterloo
#     computer science department
#     latest update .... 30 november, 1979.
#
#     given the factorization
#
#          zz*dd*rr
#
#     of some  n by k  matrix
#
#       (0 .le. k .lt. n)
#         (n .ge. 1),
#
#     where
#
#          (zz-transp)*(zz) = (dd-inv),
#          dd  is diagonal and nonsingular,
#     and
#          rr  has zeros below the diagonal,
#
#     and given a  (k+1)th  column
#     to be addedto the original matrix,
#     this program updates  zz,dd and rr.
#
#     the value of  k  is increased by one.
#
#     w  is a scratch array.
#
#     use is made of routines from the library
#     of basic linear algebra subroutines (blas).
#
#     parameters...
#
#                     input/
#       name   type   output/   sub-    description
#                     scratch  scripts
#       -------------------------------------------
#       n      int.      i              number of rows
#
#       k      int.     i/o             number of columns
#
#       zz     double precision     i/o       2     scaled orthogonal
#                                       matrix
#
#       nzzr   int.      i              row dimension of zz
#
#       dd     double precision     i/o       1     diagonal scaling
#                                       matrix (diagonal
#                                       elements only)
#
#       rr     double precision     i/o       1     right-triangular
#                                       matrix in compact form.
#
#       col    double precision      i        1     column to be
#                                       added to  rr
#
#       fail   log.      o             .true.  if  k,n
#                                       are improper
#
#       w      double precision     scr       1     workspace
#       -------------------------------------------
#
#     the  i-th  segment of the array  rr  is  n-i+2 spaces
#     long and contains  1  work space followed by the
#     k-i+1  elements of row  i  followed by  n-k
#     scratch spaces.
#     ***************
#
#     +++++++++++++++
#     blas  scopy1,sdot1,srotm1,srotmg1
#     +++++++++++++++
#
integer i,j,jdel,kp1,kp2
double precision di,one,param(5),wi,zero
#
double precision sdot1
#
#     +++++++++++++++
#     the following declarations are necessary
#     for portability when  scopy1  is used, as
#     it is below, to fill arrays with a single value
#     (one=unity  and  zero=zip  in this case).
#     +++++++++++++++
#
double precision unity(1),zip(1)
equivalence(one,unity(1)),(zero,zip(1))
#
data one/1.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
if (k<0||k>=n||n>nzzr)
	fail = .true.
else {
	if (k<=0) {
#
#     ***************
#     for the special case that the
#     factorization was vacuous,
#     reset the arrays  zz and dd
#     to represent the  identity.
#     ***************
#
		k = 0
		call scopy1(n,unity,0,dd,1)
		call scopy1(n*n,zip,0,zz,1)
		do i = 1,n
			zz(i,i) = one
		}
	kp1 = k+1
	kp2 = k+2
#
#     ***************
#     transform the incoming column,
#     and store the result in  w.
#     ***************
#
	do i = 1,n
		w(i) = sdot1(n,zz(1,i),1,col,1)
#
#     ***************
#     zero out the spike which would result from
#     storing  w  in  rr.   update  zz  and  dd.
#     ***************
#
	if (kp2<=n)
		do i = kp2,n {
			di = dd(i)
			wi = w(i)
			call srotmg1(dd(kp1),di,w(kp1),wi,param)
			w(i) = wi
			dd(i) = di
			call srotm1(n,zz(1,kp1),1,zz(1,i),1,param)
			}
#
#     ***************
#     store the new column, which is still
#     in the array  w,  into  rr.
#     ***************
#
	j = kp2
	jdel = n
	do i = 1,kp1 {
		rr(j) = w(i)
		j = j+jdel
		jdel = jdel-1
		}
	k = kp1
	fail = .false.
	}
return
end



subroutine dzdrcou(n,k,zz,nzzr,dd,rr,ic,fail)
#
integer ic,k,n,nzzr
logical fail
double precision dd(1),rr(1),zz(nzzr,1)
#
#     ***************
#     prepared by richard bartels
#     the university of waterloo
#     computer science department
#     latest update .... 30 november, 1979.
#
#     given the factorization
#
#          zz*dd*rr
#
#     of some  n by k  matrix
#
#       (1 .le. k .le. n)
#          (n .ge. 1),
#
#     where
#
#          (zz-transp)*(zz) = (dd-inv),
#          dd  is diagonal and nonsingular,
#     and
#          rr  has zeros below the diagonal,
#
#     and given the index  ic  of a column
#     to be removed  (1 .le. ic .le. k),
#     this program updates  zz,dd and rr .
#
#     the value of  k  is decreased by one, and
#     the column ordering in  rr  is changed.
#
#     use is made of routines from the library
#     of basic linear algebra subroutines (blas).
#
#     parameters...
#
#       name   type    input/   sub-    description
#                      output  scripts
#       -------------------------------------------
#       n      int.      i              number of rows
#
#       k      int.     i/o             number of columns
#
#       zz     double precision     i/o       2     scaled orthogonal
#                                       matrix
#
#       nzzr   int.      i              row dimension of zz
#
#       dd     double precision     i/o       1     diagonal scaling
#                                       matrix (diagonal
#                                       elements only)
#
#       rr     double precision     i/o       1     right-triangular
#                                       matrix in compact form.
#
#       ic     int.      i              index of column
#                                       to be removed
#
#       fail   log.      o              .true.  if  k,n,ic
#                                       are improper
#       -------------------------------------------
#
#     the  i-th  segment of the array  rr  is  n-i+2 spaces
#     long and contains  1  work space followed by the
#     k-i+1  elements of row  i  followed by  n-k
#     scratch spaces.
#     ***************
#
#     +++++++++++++++
#     blas  srotm1,srotmg1
#     +++++++++++++++
#
integer i,im1,j,jend,jinc,jstrt,km1,lstrt
double precision di,param(5),rj
#
#     /////////////////  begin program  //////////////////
#
if (k<1||k>n||n>nzzr)
	fail = .true.
else {
	km1 = k-1
#
#     ***************
#     special cases are handled first.
#     1.  k=1 and the factorization becomes null.
#     2.  ic=k and the updating is trivial.
#     ***************
#
	if (k<=1) {
		k = 0
		fail = .false.
		}
	else if (ic>=k) {
		k = km1
		fail = .false.
		}
	else {
#
#     ***************
#     general updating step.
#     the column to be deleted must be permuted
#     to the right, and subdiagonal elements
#     which result in  rr  have to be
#     transformed to zero.
#     ***************
#
		jstrt = ic+1
		jend = k
		jinc = n
		do i = 1,k {
#
#     ***************
#     permutation of the  i-th  row of rr.
#     ***************
#
			do j = jstrt,jend
				rr(j) = rr(j+1)
			if (i>ic) {
#
#     ***************
#     transformation of the current and last
#     rows  (i and i-1)  of rr  as well as
#     corresponding changes to  zz and dd.
#
#     the extra variables  di  and  rj
#     are used to avoid an error message
#     from the  pfort verifier, and they
#     may be removed, if desired, so that
#     the call to  srotmg1  would be
#
#     call srotmg1(dd(im1),dd(i),rr(lstrt),rr(jstrt),param)
#
#     ***************
#
				im1 = i-1
				di = dd(i)
				rj = rr(jstrt)
				call srotmg1(dd(im1),di,rr(lstrt),rj,param)
				rr(jstrt) = rj
				dd(i) = di
				call srotm1(jend-jstrt+1,rr(lstrt+1),1,rr(jstrt+1),1,param)
				call srotm1(n,zz(1,im1),1,zz(1,i),1,param)
				jstrt = jstrt+1
				}
#
#     ***************
#     index updating
#     ***************
#
			lstrt = jstrt
			jstrt = jstrt+jinc
			jend = jend+jinc
			jinc = jinc-1
			}
		k = km1
		fail = .false.
		}
	}
return
end



subroutine dzdrgit(n,k,zz,nzzr,rr,gv,sol,fail,w,big,eps)
#
integer k,n,nzzr
logical fail
double precision gv(1),rr(1),w(1),sol(1),zz(nzzr,1)
#
#     ***************
#     prepared by richard bartels
#     the university of waterloo
#     computer science department
#     latest update .... 30 november, 1979.
#
#     given the factorization
#
#          zz*dd*rr
#
#     of some  n by k  matrix
#
#       (1 .le. k .le. n)
#         (n .ge. 1),
#
#     where
#
#          (zz-transp)*(zz) = (dd-inv),
#          dd  is diagonal and nonsingular,
#     and
#          rr  has zeros below the diagonal,
#
#     and given an arbitrary vector  gv  of
#     appropriate dimension, this routine finds the
#     vector  sol  satisfying the underdetermined system
#
#          (zz*dd*rr-transp.)*(sol) = (gv).
#
#     that is,
#
#          (sol) = ((zz*dd*rr)-gen.inv.-transp.)*(gv).
#
#     the array  dd  is not needed by  dzdrgit.
#
#     use is made of routines from the library
#     of basic linear algebra subroutines (blas).
#
#     w  is a scratch array.
#
#     parameters...
#
#                     input/
#       name   type   output/   sub-    description
#                     scratch  scripts
#       -------------------------------------------
#       n      int.      i              number of rows
#
#       k      int.     i/o             number of columns
#
#       zz     double precision     i/o       2     scaled orthogonal
#                                       matrix
#
#       nzzr   int.      i              row dimension of zz
#
#       rr     double precision     i/o       1     right-triangular
#                                       matrix in compact form.
#
#       gv     double precision      i        1     given vector
#
#       sol    double precision      o        1     solution
#
#       fail   log.      o              .true. if  n,k
#                                       are improper, or if
#                                       rr  is singular
#
#       w      double precision     scr       1     workspace
#       -------------------------------------------
#
#     the  i-th  segment of the array  rr  is  n-i+2 spaces
#     long and contains  1  work space followed by the
#     k-i+1  elements of row  i  followed by  n-k
#     scratch spaces.
#
#     if  gv  and  sol  are dimensioned to the
#     maximum of  n  and  k , then the same
#     storage array may be used for both of
#     these vectors.
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#
#     blas  saxpy1,scopy1
#
#     big  is the largest positive number
#     which can be represented in the
#     precision of the arithmetic being used.
#     +++++++++++++++
#
integer i,j,jdel
double precision big,wi,one,rrj,zero,eps
#
#     +++++++++++++++
#     the following declarations are necessary
#     for portability when  scopy1  is used, as
#     it is below, to fill arrays with a single value
#     (zero=zip  in this case).
#     +++++++++++++++
#
double precision zip(1)
equivalence(zero,zip(1))
#
data one/1.0d+00/
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
if (k<1||k>n||n>nzzr)
	fail = .true.
else {
#
#     ***************
#     first solve  (rr-transp.)*(w) = (gv)
#     ***************
#
	call scopy1(k,gv,1,w,1)
	j = 2
	jdel = n+1
	do i = 1,k {
		rrj = rr(j)
		wi = w(i)
#* Here the check for ill-condition is NOT changed to use eps instead of big
#		if (dabs(rrj)<one)
#			if (eps*dabs(wi)>=dabs(rrj))
		if (dabs(rrj)<one)
			if (dabs(wi)>=dabs(rrj)*big)
				go to 150
		w(i) = wi/rrj
		if (i<k)
			call saxpy1(k-i,(-w(i)),rr(j+1),1,w(i+1),1)
		j = j+jdel
		jdel = jdel-1
		}
#
#     ***************
#     now  (sol) = (zz)*(w)
#     ***************
#
	call scopy1(n,zip,0,sol,1)
	do i = 1,k
		call saxpy1(n,w(i),zz(1,i),1,sol,1)
	fail = .false.
	return
	150  fail = .true.
	}
return
end



subroutine dzdrgnv(n,k,zz,nzzr,rr,gv,sol,fail,big)
#
integer k,n,nzzr
logical fail
double precision gv(1),rr(1),sol(1),zz(nzzr,1)
#
#     ***************
#     prepared by richard bartels
#     the university of waterloo
#     computer science department
#     latest update .... 30 november, 1979.
#
#     given the factorization
#
#          zz*dd*rr
#
#     of some  n by k  matrix
#
#       (1 .le. k .le. n)
#         (n .ge. 1),
#
#     where
#
#          (zz-transp)*(zz) = (dd-inv),
#          dd  is diagonal and nonsingular,
#     and
#          rr  has zeros below the diagonal,
#
#     and given an arbitrary vector  gv  of
#     appropriate dimension, this routine finds the
#     vector  sol  given by
#
#          (sol) = ((zz*dd*rr)-gen.inv.)*(gv),
#
#     which represents the least squares problem
#
#        (zz*dd*rr)*(sol) = (gv).
#
#     the array  dd  is not needed by  dzdrgnv.
#
#     use is made of routines from the library
#     of basic linear algebra subroutines (blas).
#
#     parameters...
#
#       name   type   input/    sub-    description
#                     output/  scripts
#       -------------------------------------------
#       n      int.      i              number of rows
#
#       k      int.     i/o             number of columns
#
#       zz     double precision     i/o       2     scaled orthogonal
#                                       matrix
#
#       nzzr   int.      i              row dimension of zz
#
#       rr     double precision     i/o       1     right-triangular
#                                       matrix in compact form.
#
#       gv     double precision      i        1     given vector
#
#       sol    double precision      o        1     solution
#
#       fail   log.      o              .true. if  n,k
#                                       are improper, or if
#                                       rr  is singular
#       -------------------------------------------
#
#     the  i-th  segment of the array  rr  is  n-i+2 spaces
#     long and contains  1  work space followed by the
#     k-i+1  elements of row  i  followed by  n-k
#     scratch spaces.
#     ***************
#
#     +++++++++++++++
#     system routines  dabs
#
#     blas  sdot1
#
#     big  is the largest positive number
#     which can be represented in the
#     precision of the arithmetic being used.
#     +++++++++++++++
#
integer i,ix,j,jdel
double precision big,one,tden,tnum
#
double precision sdot1
#
data one/1.0d+00/
#
#     /////////////////  begin program  //////////////////
#
if (k<1||k>n||n>nzzr)
	fail = .true.
else {
#
#     ***************
#     form   (v) = (zz(1)-transp)*(gv),   where  zz(1)
#     is the matrix of the first  k  columns of  zz
#
#     v  can be stored in the array  sol.
#     ***************
#
	do i = 1,k
		sol(i) = sdot1(n,zz(1,i),1,gv,1)
#
#     ***************
#     backsolve the system
#          (rr)*(sol) = (v)
#     for the vector  sol
#
#     note that  sol  and  v
#     are stored in the same array.
#     ***************
#
	j = (((n+1)*(n+2)-(n-k+3)*(n-k+2))/2)+2
	jdel = n-k+3
	do ix = 1,k {
		i = k-ix+1
		tden = rr(j)
		tnum = sol(i)
		if (ix>1)
			tnum = tnum-sdot1(ix-1,rr(j+1),1,sol(i+1),1)
		if (dabs(tden)<one)
			if (dabs(tnum)>=dabs(tden)*big)
				go to 160
		sol(i) = tnum/tden
		j = j-jdel
		jdel = jdel+1
		}
	fail = .false.
	return
	160  fail = .true.
	}
return
end



subroutine dzdrpoc(n,k,zz,nzzr,dd,gv,poc,fail)
#
integer k,n,nzzr
logical fail
double precision dd(1),poc(1),gv(1),zz(nzzr,1)
#
#     ***************
#     prepared by richard bartels
#     the university of waterloo
#     computer science department
#     latest update .... 30 november, 1979.
#
#     zz is an  n by n  (n .ge. 1)  scaled
#     orthogonal matrix.  dd  contains the
#     diagonal elements of a diagonal scaling
#     matrix.  gv  is a given vector of length  n.
#
#     we have
#
#          (zz-transp.)*(zz) = (dd-inv.)
#
#     and
#
#               zz*dd*rr = mat
#
#     for some  n by k  (0 .le. k .le. n)
#     matrix  rr  with zeros below the diagonal
#     and some given matrix  mat.  (niether  rr
#     nor  mat  are needed by  dzdrpoc.)
#
#     then
#
#    (proj(oc)) = (zz(2))*(dd(2))*(zz(2)-transp.)
#
#     is the (orthogonal) projector on the
#     complement of the range space of  mat,
#     where  zz(2)  represents the last  n-k
#     columns of  zz  and  dd(2)  represents the
#     lower-right-hand  n-k  order submatrix of  dd.
#
#     dzdrpoc  produces the vector
#
#               poc = (proj(oc))*gv .
#
#     use is made of routines from the library
#     of basic linear algebra subroutines (blas).
#
#     parameters...
#
#                     input/
#       name   type   output/   sub-    description
#                     scratch  scripts
#       -------------------------------------------
#       n      int.      i              order of  zz,dd
#
#       k      int.     i/o             number of columns
#                                       of  zz  defining
#                                       range of  mat
#
#       zz     double precision     i/o       2     scaled orthogonal
#                                       matrix
#
#       nzzr   int.      i              row dimension of zz
#
#       dd     double precision     i/o       1     diagonal scaling
#                                       matrix (diagonal
#                                       elements only)
#
#       gv     double precision      i        1     vector to be projected
#
#       poc    double precision      o        1     projection
#
#       fail   log.      o              .true.  if  n,k
#                                       are improper
#
#       -------------------------------------------
#
#     ***************
#
#     +++++++++++++++
#     blas  saxpy1,scopy1,sdot1
#     +++++++++++++++
#
integer i,kp1
double precision wi,zero
#
double precision sdot1
#
#     +++++++++++++++
#     the following declarations are necessary
#     for portability when  scopy1  is used, as
#     it is below, to fill arrays with a single value
#     (zero=zip  in this case).
#     +++++++++++++++
#
double precision zip(1)
equivalence(zero,zip(1))
#
data zero/0.0d+00/
#
#     /////////////////  begin program  //////////////////
#
kp1 = k+1
if (k<0||k>n||n<1||n>nzzr)
	fail = .true.
else if (k<=0) {
#
#     ***************
#     case 1 ... zz(2)=zz  (k=0)
#     ***************
#
	call scopy1(n,gv,1,poc,1)
	fail = .false.
	}
else if (k>=n) {
#
#     ***************
#     case 2 ... zz(2) is vacuous  (k=n)
#     ***************
#
	call scopy1(n,zip,0,poc,1)
	fail = .false.
	}
else {
#
#     ***************
#     case 3 ... zz(2)  is intermediate
#     between the other two cases
#     (0 .lt. k .lt. n)
#     ***************
#
	call scopy1(n,zip,0,poc,1)
	do i = kp1,n {
		wi = sdot1(n,zz(1,i),1,gv,1)*dd(i)
		call saxpy1(n,wi,zz(1,i),1,poc,1)
		}
	fail = .false.
	}
return
end

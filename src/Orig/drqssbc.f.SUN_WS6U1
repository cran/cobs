      subroutine drqssbc(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,psw,
&      e,ner,x,f,erql1n,res,indx,w,nt,nsol,sol,tl,toler,big,eps,icyc,
&      tmin,k,k0,lstart,factor)
c
c This is a modification of Bartels and Conn (1980) as described in 
c Koenker and Ng (1997), "A Remark on Bartels and Conn's Linearly Constrained
c L1 Algorithm", ACM Transaction on Mathematical Software, forthcoming.
c
c It also contains the parametric linear programming on `tau' and `lambda' as 
c described in Ng (1996), "An Algorithm for Quantile Smoothing Splines",
c Computational Statistics & Data Analysis, 22, 99-118.
c
c
c     ***************
c     front end interface
c     ***************
c
c     +++++ parameters +++++
c     ----------------------------------------------------------
c                           input
c     name   type  subscrpt  output        description
c                           scratch
c     ..........................................................
c     nrq    int.    none      in      number of observations
c                                      in the rq norm that correspond
c                                      to the fidelity
c                                      (may be zero)
c
c     nl1    int.    none      in      number of observations
c                                      in the l1 norm that correspond
c                                      to the roughness measure
c                                      (may be zero)
c
c     neqc   int.    none      in      number of equality
c                                      constraints
c                                      (may be zero)
c
c     niqc   int.    none      in      number of inequality
c                                      constraints
c                                      (may be zero)
c
c     niqc1  int.    none      in      part of niqc that belongs to
c                                      the loo roughness measure
c                                      (may be zero)
c
c     nvars  int.    none      in      number of variables
c
c     nact   int.    none      out     number of active
c                                      equations/constraints
c                                      at termination
c                                      (if any, their associated
c                                      column positions in  e  will
c                                      be listed in  indx(1)
c                                      through  indx(nact) )
c
c     ifl    int.    none      out     termination code
c                                      (see below)
c
c     mxs    int.    none      in      maximum number of steps
c                                      allowed
c
c     psw    logic.  none      in      print switch
c                                      (see below)
c
c     e      real     2        in      equation/constraint matrix
c                                      the first  nrq+nl1  columns
c                                      (see note below) specify
c                                      equations, the remaining
c                                      columns (if any) specify
c                                      constraints.
c
c     ner    int.    none      in      row dimension of e
c
c     x      real     1        in      starting values for the
c                                      unknowns (use zeros if no
c                                      guess is available)
c                              out     termination values for
c                                      the unknowns
c
c     f      real     1        in      equation/constraint
c                                      right-hand sides
c
c     erql1n real    none     out      rq-l1 norm of equation
c                                      residuals at termination
c
c     res    real     1        out     equation/constraint
c                                      residuals at termination
c
c     indx   int.     1        out     index vector used to record
c                                      the order in which the columns
c                                      of  e  are being processed
c
c     w      real     1        scr.    working storage
c     nt     int.     none     out     number of unique tau or lambda
c                                      solutions while performing parametric
c                                      in tau or lambda
c     nsol   int.     none     in      upper limit for the number of unique
c                                      tau or lambda solutions
c     sol    real     2        out     matrix of solutions when performing
c                                      parametric programming in tau or lambda
c     tl     real     1        in      values of initial tau and lambda
c     toler  real     none     in      tolerance used in parametric programming
c     big    real     none     in      largest representable floating point
c                                      number
c     eps    real     none     in      least positive number satisfying  
c                                      (1.0 + eps) .gt. 1.0
c     icyc   int.     none     out     number of cycles to achieve convergence
c     tmin   real     none     in      smallest value of tau to begin 
c                                      parametric programming in tau
c     k      int.     none     out     effective dimension of the model
c     k0     int.     none     in      the largest effective dimension of the
c                                      model allowed during parametric
c                                      programming in lambda
c     lstart real    none     in       largest value of lambda to begin 
c                                      parametric programming in lambda
c     factor real    none     in       factor to determine the how big a step
c                                      to take to the next smaller lambda 
c                                      during parametric programming in lambda
c     ----------------------------------------------------------
c
c     +++++ purpose +++++
c     ----------------------------------------------------------
c     this subroutine solves the   nrq+nl1 by nvars
c     system of equations
c
c                       (a-transpose) * x   ==   b
c
c     subject to the  neqc   constraints
c
c                       (g-transpose) * x  .eq.  h
c
c     and the  niqc  inequality constraints
c
c                       (c-transpose) * x  .ge.  d
c
c     for the unknowns  x(1),...,x(nvars).
c
c     the problem must be well-posed, nontrivial
c     and overdetermined in the sense that
c
c                          nvars .ge. 1
c                          nrq+nl1 .ge. 0
c                          neqc  .ge. 0
c                          niqc  .ge. 0
c               nrq+nl1+neqc+niqc  .ge. nvars.
c
c     further, no column of  a, g  or  c  should be zero.
c     if these conditions are not met, the program
c     will terminate without performing any substantive
c     computations.
c
c     a point  x  is a solution if it minimizes the equation
c     residuals from among all points which satisfy the
c     constraints.  at any (nondegenerate) solution
c     there will be  nact  equations and constraints
c     whose residuals
c
c          (a(i)-transpose) * x - b(i)
c
c          (g(i)-transpose) * x - h(i)
c
c     and
c
c          (c(i)-transpose) * x - d(i)
c
c     are zero.
c
c     the columns of  (a,g,c)  corresponding to the zero residuals
c     are referred to as  active columns  throughout this listing.
c     the numbers of the active columns are maintained as the
c     entries  1,...,nact  of the array  indx.
c
c     a solution  x  is found by minimizing a piecewise
c     linear penalty function formed from the  l1
c     norm of the equation residuals and the sum of the
c     infeasibilities in the constraints.
c     the minimization proceeds in a step-by-step
c     fashion, terminating after a finite number of steps.
c
c     note that  a, g  and  c  appear transposed in the
c     problem formulation.  hence it is the columns of  (a,g,c)
c     which define the equations and constraints respectively.
c
c     the array  e  is a composite of   a, g and c
c     and  f  is a composite of  b, h  and  d.
c     e  should contain  a  as its first  nrq+nl1  columns.
c     it should contain  g  as its next  neqc  columns and
c     contain  c  as its remaining  niqc  columns.
c     similarly  f  should contain  b  as its first
c     nrq+nl1  components,  h  as its next  neqc  components
c     and  d  as its last  niqc  components.
c     ----------------------------------------------------------
c
c     +++++ arrays +++++
c     ----------------------------------------------------------
c     e  is to be dimensioned at least    n  by  m,
c     x                       at least    n,
c     f                       at least    m,
c     res                     at least    m,
c     indx                    at least    m,
c     w                       at least    ((3*n*n+11*n+2)/2) + (2*m).
c
c                                         where  n = nvars  and
c                                         m = nrq+nl1+neqc+niqc
c     ----------------------------------------------------------
c
c     +++++ initialization +++++
c     ----------------------------------------------------------
c     the user must initialize
c
c          nrq,nl1,neqc,niqc,nvars,mxs,psw,e,ner,x,f .
c
c     the following are set 
c     and do not require initialization
c
c          nact,indx,res .
c
c     the array  w  is used as scratch space.
c     ----------------------------------------------------------
c
c     +++++ termination codes and intermediate printing +++++
c     ----------------------------------------------------------
c     mxs  sets a limit on the number of minimization steps to be
c     taken.
c
c     upon termination  ifl  will be set according to
c     the following code ...
c
c             ifl = 1 .... successful termination.
c
c             ifl = 2 .... unsuccessful termination.
c                          constraints cannot be satisfied.
c                          problem is infeasible.
c
c             ifl = 3 .... limit imposed by  mxs  reached
c                          without finding a solution.
c
c             ifl = 4 .... program aborted.
c                          numerical difficulties
c                          due to ill-conditioning.
c
c             ifl = 5 .... nrq, nl1, nvars, neqc and/or
c                          niqc  have improper values
c                          or  e  contains a zero column.
c
c     in all cases the output parameters  x,erql1n and res
c     will contain the values which they reached at termination.
c
c     intermediate printing will be turned off if  psw = .false. .
c     on the other hand,  details of each minimization cycle
c     will be printed if  psw  is set to  .true.
c
      integer ifl,indx(1),mxs,nact,ner,nt
      integer nrq,nl1,neqc,niqc,niqc1,nvars,nsol,icyc,k,k0
      logical psw
      double precision e(ner,1),erql1n(1),f(1),res(1),w(1),x(1)
      double precision sol(nvars+6,nsol),tl(1),toler,big,eps
      double precision tmin,lstart,factor
c
c     ////////////////  begin program  /////////////////////////
c
      call dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,psw,e,
&      ner,x,f,erql1n,res,indx,w,nt,nsol,sol,tl(1),tl(2),toler,big,eps,
&      icyc,tmin,k,k0,lstart,factor)
      return
      end
      subroutine dcrql1lt(nrq,nl1,neqc,niqc,niqc1,nvars,nact,ifl,mxs,
&      psw,e,ner,x,f,erql1n,res,indx,w,nt,nsol,sol,t,lam,toler,big,eps,
&      icyc,tmin,k,k0,lstart,factor)
c
c
c     ***************
c     main body
c     ***************
c
      integer ifl,indx(1),mxs,nact,neqc,nrq,nl1,ner,niqc,niqc1,nvars,
&      nrql1,nt
      integer nsol,k,k0
      logical psw,itend,ilend,ilfix
      double precision e(ner,1),erql1n,f(1),res(1),w(1),x(1)
      double precision eps,tmin,tmax,sol(nvars+6,nsol),t,lam,zero,one,
&      toler,big,l0,l1
      double precision tnxt,lnxt,lstart,factor
c
c
      integer ddx,grdx,grd1x,icyc,iaddc,idelc
      integer px,ptex,rrx,topx,zzx
      double precision alpha,amag,cgmag,pen,penpar,told
c
      data zero/0.d00/
      data one/1.d00/
c
c     ////////////////  begin program  /////////////////////////
c
cinitialize ifl to 0
      ifl = 0
      nrql1 = nrq+nl1
      itend = .true.
      ilend = .true.
      ilfix = .true.
      nt=1
      tnxt=t
      sol(1,nt)=t
      lnxt=lam
      sol(2,nt)=lam
      if(.not.(t.lt.zero.or.t.gt.one))goto 23000
c Note here that tmin is passed into the subroutine
         tmax = one - toler
         itend = .false.
         tnxt=tmin
         told=zero
         sol(1,nt)=tmin
23000 continue
      if(.not.(lam .lt. zero))goto 23002
         l0 = toler
         l1 = (big-toler)
         ilend = .false.
         ilfix = .false.
         lnxt=lstart
         told=t
         sol(2,nt)=lstart
23002 continue
      if(.not.(.not.itend.and..not.ilend))goto 23004
cdon't allow both t and lam to vary
         ifl = 7
         return
c
cNote: penpar is assigned outside the loop as well as drql1sup
c
23004 continue
      penpar=one/lnxt**.5
c     repeat
23006    continue
         call drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,px,ptex,
&         rrx,topx,zzx,icyc,ifl,e,ner,amag,cgmag,penpar,lnxt)
         call dnewpen(iaddc,idelc,nact,nrql1,neqc,niqc,nvars,ifl,e,ner,
&         x,f,res,w(ptex),alpha,penpar,indx)
c        repeat
23009       continue
            call drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,
&            ifl,mxs,e,ner,x,f,res,w(grdx),erql1n,pen,penpar,indx,w(zzx)
&            ,nvars,w(ddx),w(rrx),w(topx),tnxt,eps,w(grd1x))
c		call dmonit(nact,neqc,niqc,nvars,icyc,psw,x,alpha,erql1n,pen,penpar,indx)
            call drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,
&            f,res,w(grdx),w(px),erql1n,amag,cgmag,penpar,indx,w(zzx),
&            nvars,w(ddx),w(rrx),w(topx),tnxt,big,eps)
            call drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,
&            x,res,w(grdx),w(px),w(ptex),alpha,penpar,indx,w(topx),tnxt,
&            big,eps,w(grd1x),idelc)
23010       if(.not.(ifl.ne.0))goto 23009
         if(.not.(.not.(itend.and.ilend) .and. (ifl .ne.2 .or. cgmag+
&         penpar*amag .eq. cgmag)))goto 23012
            call drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,nact,nvars,nrq,
&            nl1,neqc,niqc,niqc1,w(zzx),nvars,w(rrx),w(grdx),w(px),w(
&            topx),w(topx+nvars),ifl,idelc,iaddc,icyc,alpha,amag,cgmag,
&            psw,penpar,nsol,sol,x,ilfix,l0,l1,tnxt,lnxt,toler,erql1n,
&            eps,big,told,k0,factor)
cupdate penpar to the sqrt of next lambda
            penpar=one/lnxt**.5
            if(.not.(ifl .ne. 0))goto 23014
               goto 23008
23014       continue
            goto 23013
c        else
23012       continue
            if(.not.(ifl .ne.2 .or. cgmag+penpar*amag .eq. cgmag))goto 2
&            3016
               goto 23008
23016       continue
23013    continue
c
ccompute the effective dimensionality
23007    goto 23006
23008 continue
      k = 0
      do 23018 i=1,nact
         if(.not.(indx(i) .le. nrq .or. (indx(i) .gt. nrq+nl1 .and. 
&         indx(i) .le. nrq+nl1+neqc) .or. indx(i) .gt. nrq+nl1+neqc+
&         niqc1))goto 23020
            k = k+1
23020    continue
23018    continue
      return
      end
      subroutine drql1nlt(nt,tmin,tmax,res,e,f,ner,indx,nact,nvars,nrq,
&      nl1,neqc,niqc,niqc1,zz,nzzr,rr,a,aa,b,bb,ifl,idelc,iaddc,icyc,
&      alpha,amag,cgmag,psw,penpar,nsol,sol,x,ilfix,l0,l1,tnxt,lnxt,
&      toler,erql1n,eps,big,told,k0,factor)
c
c
c     ***************
c     perform parametric programming in "lambda" and "tau"
c     ***************
c
      logical fail,psw,ilfix
      integer nt,nact,nvars,nrq,nl1,ner,nzzr,ifl,k,nactp1
      integer nrql1,neqc,niqc,niqc1,nallq,nalqp1,ncols,nqnp1,ix
      integer idelc,iaddc,icyc,indx(1),isave,nsol,k0
      double precision tmin,tmax,sgn,one,res(1),zero,e(ner,1),zz(nzzr,1)
      double precision rr(1),a(1),aa(1),b(1),bb(1),thet,tnxt,tmp,sol(
&      nvars+6,nsol)
      double precision eps,two,penpar,x(1),lnxt,lamb,l0,l1,big,toler,
&      erql1n
      double precision test,prod,f(1),amag,cgmag,fidel,penal,wgt,told,
&      factor
c
c
c     ////////////////  begin program  /////////////////////////
c
      data one/1.d00/
      data zero/0.d00/
      data two/2.d00/
c
      nrql1=nrq+nl1
      nallq=nrql1+neqc
      nalqp1=nallq+1
      ncols=nallq+niqc
      nqnp1=nrql1+1
      nactp1=nact+1
      thet = tnxt
      if(.not.(ilfix))goto 23022
         tnxt = one+eps
23022 continue
      lamb = lnxt
      if(.not.(.not.ilfix))goto 23024
         lnxt = l0
23024 continue
      if(.not.(ifl.eq.1.or.ifl.eq.3))goto 23026
         call scopy1(nvars,zero,0,a,1)
         call scopy1(nvars,zero,0,b,1)
         if(.not.(nacpt.le.ncols))goto 23028
            do 23030 i = nactp1,ncols 
               ix = indx(i)
               sgn = dsign(one,res(ix))
               test = dabs(f(ix))
               do 23032 j = 1,nvars 
                  prod = dabs(e(j,ix)*x(j))
                  if(.not.(prod .gt. test))goto 23034
                     test = prod
23034             continue
23032             continue
               test = eps*dsqrt(dfloat(nvars))*test
               if(.not.(dabs(res(ix)) .lt. test))goto 23036
                  sgn = zero
23036          continue
               if(.not.(ilfix))goto 23038
                  if(.not.(ix.le.nrq))goto 23040
                     call saxpy1(nvars,(one+sgn),e(1,ix),1,a,1)
                     call saxpy1(nvars,-two,e(1,ix),1,b,1)
                     goto 23041
c                 else
23040                continue
                     if(.not.(ix.le.nallq.or.sgn.le.zero))goto 23042
                        call saxpy1(nvars,sgn,e(1,ix),1,a,1)
23042                continue
23041             continue
                  goto 23039
c              else
23038             continue
                  if(.not.(ix.le.nrq))goto 23044
                     call saxpy1(nvars,(one-two*thet+sgn),e(1,ix),1,a,1)
                     goto 23045
c                 else
23044                continue
                     if(.not.(ix.le.nrql1))goto 23046
                        call saxpy1(nvars,sgn/lamb,e(1,ix),1,b,1)
                        goto 23047
c                    else
23046                   continue
                        if(.not.(ix.le.allq.or.sgn.le.zero))goto 23048
                           call saxpy1(nvars,sgn,e(1,ix),1,a,1)
23048                   continue
23047                continue
23045             continue
23039          continue
23030          continue
23028    continue
         call dzdrgnv(nvars,nact,zz,nzzr,rr,a,aa,fail,big)
         if(.not.(fail))goto 23050
            ifl = 4
            goto 23051
c        else
23050       continue
            call dzdrgnv(nvars,nact,zz,nzzr,rr,b,bb,fail,big)
            if(.not.(fail))goto 23052
               ifl = 4
               goto 23053
c           else
23052          continue
               do 23054 i = 1,nact 
                  ix = indx(i)
ca check for small bb(i) is implemented to avoid floating point overflow
                  test = dabs(f(ix))
                  do 23056 j = 1,nvars 
                     prod = dabs(e(j,ix)*x(j))
                     if(.not.(prod .gt. test))goto 23058
                        test = prod
23058                continue
23056                continue
                  test = eps*dsqrt(dfloat(nvars))*test
                  if(.not.(ix.le.nrq))goto 23060
                     if(.not.(ilfix))goto 23062
                        tmp = (two+aa(i))/(two-bb(i))
                        if(.not.(tmp.lt.tnxt.and.tmp.ge.thet))goto 23064
                           tnxt = tmp
                           isave = i
                           goto 23065
c                       else
23064                      continue
                           tmp = aa(i)/(two-bb(i))
                           if(.not.(tmp.lt.tnxt.and.tmp.ge.thet))goto 23
&                           066
                              tnxt = tmp
23066                      continue
23065                   continue
                        goto 23063
c                    else
23062                   continue
                        tmp = (two*thet - aa(i))/bb(i)
                        if(.not.(dabs(bb(i)).lt.test))goto 23068
cavoid bb near zero
                           tmp = big
23068                   continue
                        if(.not.(tmp.gt.lnxt.and.tmp.lt.lamb))goto 23070
                           lnxt = tmp
                           isave = i
                           goto 23071
c                       else
23070                      continue
                           tmp = (two*thet-two-aa(i))/bb(i)
                           if(.not.(dabs(bb(i)).lt.test))goto 23072
cavoid bb near zero
                              tmp = big
23072                      continue
                           if(.not.(tmp.gt.lnxt.and.tmp.lt.lamb))goto 23
&                           074
                              lnxt = tmp
                              isave = i
23074                      continue
23071                   continue
23063                continue
                     goto 23061
c                 else
23060                continue
                     if(.not.(ix.le.nallq))goto 23076
                        if(.not.(ilfix))goto 23078
                           tmp = (one-aa(i))/bb(i)
                           if(.not.(dabs(bb(i)).lt.test))goto 23080
cavoid bb near zero
                              tmp = big
23080                      continue
                           if(.not.(tmp.lt.tnxt.and.tmp.ge.thet))goto 23
&                           082
                              tnxt = tmp
                              isave = i
                              goto 23083
c                          else
23082                         continue
                              tmp = -(aa(i)+one)/bb(i)
                              if(.not.(dabs(bb(i)).lt.test))goto 23084
cavoid bb near zero
                                 tmp = big
23084                         continue
                              if(.not.(tmp.lt.tnxt.and.tmp.ge.thet))
&                              goto 23086
                                 tnxt = tmp
                                 isave = i
23086                         continue
23083                      continue
                           goto 23079
c                       else
23078                      continue
                           tmp = -aa(i)*lamb/(bb(i)*lamb+one)
                           if(.not.(tmp.gt.lnxt.and.tmp.lt.lamb))goto 23
&                           088
                              lnxt = tmp
                              isave = i
                              goto 23089
c                          else
23088                         continue
                              tmp = -aa(i)*lamb/(bb(i)*lamb-one)
                              if(.not.(tmp.gt.lnxt.and.tmp.lt.lamb))
&                              goto 23090
                                 lnxt = tmp
                                 isave = i
23090                         continue
23089                      continue
23079                   continue
                        goto 23077
c                    else
23076                   continue
                        if(.not.(ilfix))goto 23092
                           tmp = -aa(i)/bb(i)
                           if(.not.(dabs(bb(i)).lt.test))goto 23094
cavoid bb near zero
                              tmp = big
23094                      continue
                           if(.not.(tmp.lt.tnxt.and.tmp.ge.thet))goto 23
&                           096
                              tnxt = tmp
                              isave = i
23096                      continue
                           goto 23093
c                       else
23092                      continue
                           tmp = -aa(i)/bb(i)
                           if(.not.(dabs(bb(i)).lt.test))goto 23098
cavoid bb near zero
                              tmp = big
23098                      continue
                           if(.not.(tmp.gt.lnxt.and.tmp.lt.lamb))goto 23
&                           100
                              lnxt = tmp
                              isave = i
23100                      continue
23093                   continue
23077                continue
23061             continue
23054             continue
23053       continue
23051    continue
c
ccompute the effective dimensionalty, fidelity and penalty
c
23026 continue
      k = 0
      fidel = zero
      penal = zero
      do 23102 i=1,nact
         if(.not.(indx(i) .le. nrq .or. (indx(i) .gt. nrq+nl1 .and. 
&         indx(i) .le. nrq+nl1+neqc) .or. indx(i) .gt. nrq+nl1+neqc+
&         niqc1))goto 23104
            k = k+1
cset the lower stopping criterion for lambda to be either k>=k0 or when
clnxt < 0
23104    continue
23102    continue
      tmp = lnxt-lnxt*10.0d0**(factor-4.0d0)
      if(.not.(k.ge.k0.or.tmp.lt.zero))goto 23106
         l0 = lnxt+eps
23106 continue
      do 23108 i=nactp1,ncols 
         ix = indx(i)
         tmp = res(ix)
         wgt = dsign(one,tmp)
         if(.not.(ix.le.nrq))goto 23110
            wgt = one-two*told+wgt
            fidel = fidel+wgt*tmp
            goto 23111
c        else
23110       continue
            if(.not.(ix.le.nrql1))goto 23112
               penal = penal+dabs(tmp)
23112       continue
23111    continue
23108    continue
      if(.not.(ilfix))goto 23114
         if(.not.((ifl.eq.1.or.ifl.eq.3) .and. tnxt.lt.tmax))goto 23116
            nt = nt+1
            sol(1,nt) = tnxt
            sol(2,nt) = lamb
            sol(3,nt-1)=dble(ifl)
            sol(4,nt-1) = fidel
            sol(5,nt-1) = penal/lamb
            sol(6,nt-1) = k
            told = tnxt
            if(.not.(indx(isave) .le. nrql1))goto 23118
               tnxt = tnxt+10.0d0**(factor-7.0d0)*amag
               goto 23119
c           else
23118          continue
               tnxt = tnxt+10.0d0**(factor-7.0d0)*cgmag
23119       continue
            call scopy1(nvars,x,1,sol(7,nt-1),1)
            idelc = 0
            iaddc = nact
            icyc = -1
            ifl = 0
            alpha = zero
            goto 23117
c        else
23116       continue
            nt = nt+1
            if(.not.(nt.ge.nsol))goto 23120
               ifl = 6
23120       continue
            sol(1,nt) = tnxt
            sol(2,nt) = lamb
            sol(3,nt-1)=dble(ifl)
            sol(4,nt-1) = fidel
            sol(5,nt-1) = penal/lamb
            sol(6,nt-1) = k
            if(.not.((ifl.eq.1.or.ifl.eq.3) .and. tnxt .ge.tmax))goto 23
&            122
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
23122       continue
            call scopy1(nvars,x,1,sol(7,nt-1),1)
23117    continue
         goto 23115
c     else
23114    continue
         if(.not.((ifl.eq.1.or.ifl.eq.3) .and. lnxt.gt.l0))goto 23124
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
            goto 23125
c        else
23124       continue
            nt = nt+1
            if(.not.(nt.ge.nsol))goto 23126
               ifl = 6
23126       continue
            sol(1,nt) = thet
            sol(2,nt) = lnxt
            sol(3,nt-1)=dble(ifl)
            sol(4,nt-1) = fidel
            sol(5,nt-1) = penal/lamb
            sol(6,nt-1) = k
            if(.not.((ifl.eq.1.or.ifl.eq.3) .and. lnxt .le.l0))goto 2312
&            8
               sol(1,nt) = thet
               sol(2,nt) = l0
               sol(3,nt) = sol(3,nt-1)
               sol(4,nt) = sol(4,nt-1)
               sol(5,nt) = sol(5,nt-1)
               sol(6,nt) = sol(6,nt-1)
               call scopy1(nvars,x,1,sol(7,nt),1)
23128       continue
            call scopy1(nvars,x,1,sol(7,nt-1),1)
23125    continue
23115 continue
c remove lambda from e for next iteration
      if(.not.(nrql1.ge.nrq+1))goto 23130
         do 23132 i=nrq+1,nrql1
            do 23134 j=1,ner
               e(j,i) = e(j,i)/lamb
23134          continue
23132       continue
23130 continue
      return
      end
      subroutine drql1sup(nrq,nl1,neqc,niqc,nvars,ddx,grdx,grd1x,px,
&      ptex,rrx,topx,zzx,icyc,ifl,e,ner,amag,cgmag,penpar,lam)
c
      integer ddx,grdx,icyc,ifl,neqc,nrql1,ner,nrq,nl1,grd1x
      integer niqc,nvars,px,ptex,rrx,topx,zzx
      double precision amag,cgmag,e(ner,1),penpar,lam
c
c     ***************
c     crql1  version.
c
c     set up the program
c     parameters and indices.
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c     +++++++++++++++
c
      integer i,j,ncols,nqnp1
      double precision oct,tmp,zero
c
      data oct/8.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
c     ***************
c     check validity of problem dimensions
c     ***************
c
      nrql1=nrq+nl1
      ncols = nrql1+neqc+niqc
      if(.not.(nvars.lt.1.or.neqc.lt.0.or.niqc.lt.0.or.nrql1.lt.0.or.
&      ncols.lt.nvars.or.ner.lt.nvars))goto 23136
         ifl = 5
         goto 23137
c     else
23136    continue
c
c     ***************
c     set up indices for the temporary storage vector  w.
c     ***************
c
         nqnp1 = nrql1+1
         grdx = 1
         grd1x = grdx+nvars
         px = grd1x+nvars
         ptex = px+nvars
         ddx = ptex+ncols
         rrx = ddx+nvars
         zzx = rrx+(((nvars+1)*(nvars+2))/2)
         topx = zzx+nvars*nvars
c
c     ***************
c     update e with lambda only if ifl!=2, i.e. update only for the new lambda
c     ***************
c
         if(.not.( ifl.ne.2))goto 23138
            do 23140 i=nrq+1,nrql1
               do 23142 j=1,ner
                  e(j,i)=e(j,i)*lam
23142             continue
23140          continue
c
c     ***************
c     amag  is a rough estimate of the norm of  a.
c     cgmag  is a rough estimate of the norm of  (g,c).
c     together they are used to determine when the
c     penalty parameter is too small and when the
c     restricted gradient is zero.
c     ***************
c
23138    continue
         amag = zero
         cgmag = zero
         if(.not.(1.le.nrql1))goto 23144
            do 23146 j = 1,nrql1 
               tmp = zero
               do 23148 i = 1,nvars
                  tmp = tmp+dabs(e(i,j))
23148             continue
               if(.not.(tmp.le.zero))goto 23150
                  go to 10
23150          continue
               if(.not.(tmp.gt.amag))goto 23152
                  amag = tmp
23152          continue
23146          continue
            go to 20
10          ifl = 5
            return
23144    continue
20       if(.not.(nqnp1.le.ncols))goto 23154
            do 23156 j = nqnp1,ncols 
               tmp = zero
               do 23158 i = 1,nvars
                  tmp = tmp+dabs(e(i,j))
23158             continue
               if(.not.(tmp.le.zero))goto 23160
                  go to 30
23160          continue
               if(.not.(tmp.gt.cgmag))goto 23162
                  cgmag = tmp
23162          continue
23156          continue
            go to 40
30          ifl = 5
            return
c
c     ***************
c     initialize  ifl,icyc,penpar
c     ***************
c
23154    continue
40       ifl = 2
         icyc = -1
23137 continue
      return
      end
      subroutine dnewpen(iaddc,idelc,nact,neqns,neqc,niqc,nvars,ifl,e,
&      ner,x,f,res,pte,alpha,penpar,indx)
c
      integer iaddc,idelc,ifl,indx(1),nact
      integer neqc,neqns,ner,niqc,nvars
      double precision alpha,e(ner,1),f(1),penpar,pte(1),res(1),x(1)
c
c     ***************
c     cl1  version.
c
c     begin a round of minimization steps
c     with a new penalty parameter value.
c     ***************
c
c     +++++++++++++++
c     blas  sdot1
c     +++++++++++++++
c
      integer i,ncols
      double precision oct,one,zero
c
      double precision sdot1
c
      data zero/0.0d+00/
      data one/1.0d+00/
      data oct/8.0d+00/
c
c     /////////////////  begin program  //////////////////
c
c     ***************
c     set penalty parameter value.
c     erase record of active equation/constraints.
c     ***************
c
      if(.not.(ifl.eq.2))goto 23164
         ncols = neqns+neqc+niqc
         ifl = 0
         nact = 0
         iaddc = 0
         idelc = 0
         alpha = zero
         penpar = penpar/oct
c
c     ***************
c     initialize  indx,res,pte,indx
c     ***************
c
         do 23166 i = 1,ncols 
            res(i) = sdot1(nvars,e(1,i),1,x,1)-f(i)
            pte(i) = zero
            indx(i) = i
23166       continue
23164 continue
      return
      end
      subroutine drql1up(iaddc,idelc,nact,nrq,nl1,neqc,niqc,nvars,icyc,
&      ifl,mxs,e,ner,x,f,res,grd,erql1n,pen,penpar,indx,zz,nzzr,dd, rr,
&      w,theta,eps,grd1)
c
      integer iaddc,idelc,icyc,ifl,indx(1),mxs
      integer nact,neqc,nrq,nl1,ner,niqc,nvars,nzzr
      double precision dd(1),e(ner,1),erql1n,f(1),grd(1),pen,penpar,
&      theta,grd1(1)
      double precision res(1),rr(1),w(1),x(1),zz(nzzr,1),eps
c
c     ***************
c     crql1  version.
c
c     preparation for next minimization step.
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c     +++++++++++++++
c
      integer nallq,ncols,nrql1
      double precision one,zero
c
      data one/1.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
c     ***************
c     determine the active equations and active
c     constraints.  compute residuals and function value.
c     update the  z*d*r  decomposition.
c     ***************
c
      nrql1 = nrq+nl1
      nallq = nrql1+neqc
      ncols = nallq+niqc
      if(.not.(ifl.eq.0))goto 23168
         icyc = icyc+1
         if(.not.(icyc.gt.mxs))goto 23170
            ifl = 3
            goto 23171
c        else
23170       continue
            call ddelcol1(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,indx)
            call dresid(iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
            call daddcol(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,e,ner,
&            indx,w,eps)
            call drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,
&            res,grd,erql1n,pen,penpar,indx,theta,grd1)
23171    continue
23168 continue
      return
      end
      subroutine drql1fp(idelc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,x,
&      f,res,grd,p,erql1n,amag,cgmag,penpar,indx,zz,nzzr,dd,rr,w, theta,
&      big,eps)
c
      integer idelc,ifl,indx(1),nact,neqc,nrq,nl1,ner,niqc,nvars,nzzr,
&      nrql1
      double precision amag,cgmag,dd(1),e(ner,1),erql1n,f(1),grd(1),p(1)
&      ,penpar
      double precision res(1),rr(1),w(1),x(1),zz(nzzr,1),theta
c
c     ***************
c     crql1  version.
c
c     determine descent direction  p
c     (or discover optimality)
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c
c     blas  sasum1,scopy1,sscal1
c
c     eps  is the smallest positive number which
c     satisfies   (1.0 + eps) .gt. 1.0   in the
c     precision of the arithmetic being used.
c     (alternatively, for less strict zero checking,
c      eps  can be set to a user-specified tolerance.)
c     +++++++++++++++
c
      integer coefx,i,ix,nallq,nalqp1,ncols,nqnp1,topx,nrql1
      logical fail
      double precision grdnrm,one,pnrm,prod,test,zero
      double precision eps,big
c
      double precision sasum1
c
c     +++++++++++++++
c     the following declarations are necessary
c     for portability when  scopy1  is used, as
c     it is below, to fill arrays with a single value
c     (one=unity  and  zero=zip  in this case).
c     +++++++++++++++
c
      double precision unity(1),zip(1)
      equivalence(one,unity(1)),(zero,zip(1))
c
      data one/1.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      nrql1 = nrq+nl1
      idelc = 0
      if(.not.(ifl.eq.0))goto 23172
         nallq = nrql1+neqc
         nalqp1 = nallq+1
         ncols = nallq+niqc
         nqnp1 = nrql1+1
         coefx = 1
         topx = coefx+nvars
c
c     ***************
c     project the negative of the restricted gradient
c     onto the orthogonal complement of the space
c     spanned by the active columns.
c     ***************
c
         call dzdrpoc(nvars,nact,zz,nzzr,dd,grd,p,fail)
         if(.not.(fail))goto 23174
            ifl = 4
            goto 23175
c        else
23174       continue
            call sscal1(nvars,-one,p,1)
            pnrm = sasum1(nvars,p,1)
            grdnrm = sasum1(nvars,grd,1)
c
c     ***************
c     if the projection is not zero,
c     it will serve as a descent direction.
c
c     otherwise find the representation of
c     the restricted gradient as a linear
c     combination of the active columns.
c     the coefficients of the linear combination
c     are to be stored in the array  coef
c     (that is, in  w(coefx),...,w(coefx+nact-1)).
c     ***************
c
            if(.not.(pnrm.le.eps*(amag*penpar+cgmag)))goto 23176
               if(.not.(nact.ne.0))goto 23178
                  call dzdrgnv(nvars,nact,zz,nzzr,rr,grd,w(coefx),fail,
&                  big)
                  if(.not.(fail))goto 23180
                     ifl = 4
                     return
c                 else
23180                continue
c
c     ***************
c     convert the coefficients of the linear
c     combination into a descent direction  p ,
c     or determine optimality.
c
c     if the optimality test is not satisfied,
c     drql1gv  will indicate an equation/constraint
c     to be deleted from activity by the value
c     of  idelc.  for optimality,  idelc=0.
c     ***************
c
                     call drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,
&                     grd,w(coefx),penpar,indx,theta,eps)
                     pnrm = zero
                     if(.not.(idelc.ne.0))goto 23182
                        call dzdrgit(nvars,nact,zz,nzzr,rr,w(coefx),p,
&                        fail,w(topx),big,eps)
                        if(.not.(.not.fail))goto 23184
                           pnrm = sasum1(nvars,p,1)
23184                   continue
                        if(.not.(fail))goto 23186
                           ifl = 4
                           return
23186                   continue
c
c     ***************
c     if a descent direction  p  could have been found,
c     it has been obtained by this point in the program.
c
c     check for optimality.
c
c     pnrm  has been set exactly zero
c     after the call to subroutine  drql1gv
c     if the optimality conditions are satisfied.
c     the check below has been made somewhat
c     complicated to allow for the rare event that
c     the restricted gradient is zero and no
c     columns are active,  or that the  rq  norm of
c               (a-transpose) * x - f
c     is computationally zero.
c     (the call to the subroutine  refine
c      may be omitted, if desired.)
c     ***************
c
23182                continue
                     if(.not.(pnrm.gt.eps*(amag*penpar+cgmag)))goto 2318
&                     8
                        do 23190 i = 1,nrql1 
                           test = dabs(f(i))
                           do 23192 ix = 1,nvars 
                              prod = dabs(e(ix,i)*x(ix))
                              if(.not.(prod.gt.test))goto 23194
                                 test = prod
23194                         continue
23192                         continue
                           if(.not.(dabs(res(i)).gt.eps*test))goto 23196
                              return
23196                      continue
23190                      continue
23188                continue
23181             continue
23178          continue
               ifl = 1
c			call drql1rf(nact,nrq,nl1,ncols,nvars,ifl,e,ner,x,f,erql1n,res,indx,zz,nzzr,rr,w,theta,big,eps)
               if(.not.(ifl.eq.1))goto 23198
c
c     ***************
c     if the problem has constraints,
c     check feasibility.
c     ***************
c
                  if(.not.(nqnp1.le.ncols))goto 23200
                     do 23202 i = nqnp1,ncols 
                        test = dabs(f(i))
                        do 23204 ix = 1,nvars 
                           prod = dabs(e(ix,i)*x(ix))
                           if(.not.(prod.gt.test))goto 23206
                              test = prod
23206                      continue
23204                      continue
c NOTE: the criterion for checking feasibility is relaxed by (eps * test)^.5
c rather than eps*test
                        test = (eps*test)**.5
c						test = eps*test
                        if(.not.(i.gt.nallq))goto 23208
                           if(.not.(res(i).lt.(-test)))goto 23210
                              go to 20
23210                      continue
                           goto 23209
c                       else
23208                      continue
                           if(.not.(dabs(res(i)).gt.test))goto 23212
                              go to 10
23212                      continue
23209                   continue
23202                   continue
                     return
10                   ifl = 2
                     return
20                   ifl = 2
23200             continue
23198          continue
23176       continue
23175    continue
23172 continue
      return
      end
      subroutine drql1stp(iaddc,nact,nrq,nl1,neqc,niqc,nvars,ifl,e,ner,
&      x,res,grd,p,pte,alpha,penpar,indx,alf,theta,big,eps,grd1,idelc)
c
      integer iaddc,ifl,indx(1),nact,neqc,nrq,nl1,nrql1,ner,niqc,nvars,
&      idelc
      double precision alpha,alf(1),e(ner,1),grd(1),p(1),grd1(1)
      double precision penpar,pte(1),res(1),x(1)
      double precision theta,sgn1
c
c     ***************
c     cl1  version.
c
c     piecewise linear line search.
c     ***************
c
c     +++++++++++++++
c     system routines dabs,dsign
c
c     blas  sasum1,saxpy1,sdot1
c
c     eps  is the smallest positive number which
c     satisfies   (1.0 + eps) .gt. 1.0   in the
c     precision of the arithmetic being used.
c     (alternatively, for less strict zero checking,
c      eps  can be set to a user-specified tolerance.)
c
c     big  is the largest positive number
c     which can be represented in the
c     precision of the arithmetic being used.
c     +++++++++++++++
c
      integer i,iin,ix,jx,nactp1,nallq,ncols,num,numnac
      double precision big,den,eps,grdnrm,one,pnrm,ptg,ptg1
      double precision ratio,resid,tmp,two,zero
c
      double precision sasum1,sdot1,etp
c
      data one/1.0d+00/
      data two/2.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
c     ***************
c     this routine determines all of the ratios  alf
c     of the form
c        -res(i)/((e(.,i)-transp)*p),
c              for  i = k+1,...,mpl
c     which are nonnegative and hence indicate distances
c     from the point  x  to breakpoints which will
c     be encountered in travel along direction  p.
c     the index vector  indx  is rearranged so that
c     its  k+1  through  num  components correspond to
c     these nonnegative ratios.
c     the results are heaped so that the  alf  values can
c     be inspected in order from smallest to largest.
c     the breakpoint  alpha  giving the minimum objective
c     function value is found, and  x  is
c     adjusted to  x + alpha*p .
c
c     the inner products  (e(.,i)-transpose)*p  are saved
c     for later use in updating the residual values.
c     ***************
c
      alpha = zero
      if(.not.(ifl.eq.0))goto 23214
         nrql1 = nrq+nl1
         nallq = nrql1+neqc
         ncols = nallq+niqc
         nactp1 = nact+1
         num = 0
         if(.not.(1.le.nact))goto 23216
            do 23218 i = 1,nact 
               ix = indx(i)
               pte(ix) = sdot1(nvars,e(1,ix),1,p,1)
23218          continue
cupdate the correct gradient
23216    continue
         if(.not.(nactp1.le.iaddc))goto 23220
            do 23222 i = nactp1,iaddc
               ix = indx(i)
               etp = sdot1(nvars,e(1,ix),1,p,1)
               sgn1 = dsign(one,etp)
               if(.not.(ix.le.nrq))goto 23224
                  sgn1 = one-two*theta + sgn1
23224          continue
               if(.not.(ix.le.nallq.or.sgn1.le.zero))goto 23226
                  if(.not.(ix.le.nrql1))goto 23228
                     sgn1 = sgn1*penpar
23228             continue
                  call saxpy1(nvars,sgn1,e(1,ix),1,grd1,1)
23226          continue
23222          continue
23220    continue
         if(.not.(idelc.ne.0))goto 23230
            ix=indx(idelc)
            etp=sdot1(nvars,e(1,ix),1,p,1)
            sgn1=dsign(one,etp)
            if(.not.(ix.le.nrq))goto 23232
               sgn1 = one-two*theta + sgn1
23232       continue
            if(.not.(ix.le.nallq.or.sgn1.le.zero))goto 23234
               if(.not.(ix.le.nrql1))goto 23236
                  sgn1=sgn1*penpar
23236          continue
               call saxpy1(nvars,sgn1,e(1,ix),1,grd1,1)
23234       continue
23230    continue
         if(.not.(nactp1.gt.ncols))goto 23238
            ifl = 1
            goto 23239
c        else
23238       continue
            do 23240 i = nactp1,ncols 
               ix = indx(i)
               resid = res(ix)
               den = sdot1(nvars,e(1,ix),1,p,1)
               pte(ix) = den
               if(.not.(dsign(one,resid).ne.dsign(one,den).or.resid.eq.
&               zero))goto 23242
                  resid = dabs(resid)
                  den = dabs(den)
                  if(.not.(den.lt.one))goto 23244
                     if(.not.(resid.ge.den*big))goto 23246
                        goto 23240
23246                continue
23244             continue
                  ratio = resid/den
                  num = num+1
                  numnac = num+nact
                  jx = indx(numnac)
                  indx(numnac) = ix
                  indx(i) = jx
                  alf(num) = ratio
23242          continue
23240          continue
            if(.not.(num.le.0))goto 23248
               ifl = 2
               goto 23249
c           else
23248          continue
c
c     ***************
c     heap the positive ratios
c     ***************
c
               call ddkheap(.true.,num,indx(nactp1),alf)
c
c     ***************
c     travel along  p  until no further decrease in the
c     penalty function is possible
c     ***************
c
               iin = num
               ptg = sdot1(nvars,grd,1,p,1)
               ptg1 = sdot1(nvars,grd1,1,p,1)
               pnrm = sasum1(nvars,p,1)
               grdnrm = sasum1(nvars,grd1,1)
               do 23250 i = 1,num 
                  ix = indx(nactp1)
                  if(.not.(res(ix).eq.zero))goto 23252
                     tmp = zero
                     goto 23253
c                 else
23252                continue
                     tmp = -dsign(one,res(ix))
23253             continue
                  if(.not.(ix.le.nallq))goto 23254
                     tmp = tmp*two
23254             continue
                  if(.not.(ix.le.nrql1))goto 23256
                     tmp = tmp*penpar
23256             continue
                  ptg1 = ptg1+tmp*pte(ix)
                  if(.not.(ptg1.ge.(-eps)*grdnrm*pnrm))goto 23258
                     go to 140
23258             continue
                  call ddkheap(.false.,iin,indx(nactp1),alf)
23250             continue
               ifl = 2
               return
140            iaddc = nactp1
c
c     ***************
c     adjust  x  to  x + alpha*p
c     ***************
c
               alpha = alf(1)
               call saxpy1(nvars,alpha,p,1,x,1)
23249       continue
23239    continue
23214 continue
      return
      end
      subroutine drql1rf(nact,nrq,nl1,ncols,nvars,ifl,e,ner,x,f,erql1n,
&      res,indx,zz,nzzr,rr,w,theta,big,eps)
c
      integer ifl,indx(1),nact,ncols,nrq,nl1,ner,nvars,nzzr,nrql1
      double precision e(ner,1),erql1n,f(1),res(1),rr(1),w(1),x(1),zz(
&      nzzr,1)
      double precision theta, wgt
c
c     ***************
c     a routine for refining the solution
c     produced by  crql1.
c
c     (this routine may be omitted if desired.)
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c
c     blas  sdot1
c     +++++++++++++++
c
c#integer i,ix
      logical fail
      double precision tmp,zero,one,two
c
      double precision sdot1,big,eps
c
      data zero/0.0d+00/
      data one/1.0d+00/
      data two/2.0d+00/
c
c     /////////////// begin program ///////////////
c
      nrql1 = nrq+nl1
      if(.not.(nact.ne.0))goto 23260
         if(.not.(fail))goto 23262
            ifl = 4
            goto 23263
c        else
23262       continue
            erql1n = zero
            do 23264 i = 1,ncols 
               tmp = sdot1(nvars,e(1,i),1,x,1)-f(i)
               wgt = dsign(one,tmp)
               if(.not.(i.le.nrq))goto 23266
                  wgt = one-two*theta+wgt
23266          continue
               res(i) = tmp
               if(.not.(i.le.nrql1))goto 23268
                  erql1n = erql1n+wgt*tmp
23268          continue
23264          continue
23263    continue
23260 continue
      return
      end
c
c     ---------------
c     third level subroutines --
c          ddelcol1,dresid,addcol,drql1obj,getv
c     ---------------
c
      subroutine ddelcol1(iaddc,idelc,nact,nrow,zz,nzzr,dd,rr,indx)
c
      integer indx(1),nzzr,nact,idelc,iaddc,nrow
      double precision dd(1),rr(1),zz(nzzr,1)
c
c     ***************
c     cl1  version of  ddelcol1.
c
c     this routine administers the deletion of the column
c     indicated by the value of idelc
c     from an  nrow by nact   z*d*r   decomposition.
c     note that the value of idelc
c     is the number of a column in the decomposition
c     rather than a number which refers to
c     a column in the matrix  e.
c     (the  e-column  numbers corresponding to
c     the columns of the factorization are to be
c     found in   indx(1),...,indx(nact) .
c     the contents of   indx(nact+1),...,indx(iaddc)
c     indicate columns of  e  which are slated for
c     addition to the decomposition.)
c     the vector  indx   is rearranged by
c     permuting the element which corresponds to
c     the deletion out to the   iaddc-th  position.
c     nact  and  iaddc  are decreased accordingly.
c     ***************
c
      integer i,idlp1,ixdlc
      logical fail
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(idelc.ne.0))goto 23270
         idlp1 = idelc+1
         ixdlc = indx(idelc)
         do 23272 i = idlp1,iaddc
            indx(i-1) = indx(i)
23272       continue
         indx(iaddc) = ixdlc
         iaddc = iaddc-1
         call dzdrcou(nrow,nact,zz,nzzr,dd,rr,idelc,fail)
         idelc = ixdlc
23270 continue
      return
      end
      subroutine dresid(iaddc,nact,ncols,nvars,e,ner,x,f,res,indx,eps)
c
      integer iaddc,indx(1),nact,ncols,ner,nvars
      double precision e(ner,1),f(1),res(1),x(1)
c
c
c     ***************
c     compute the residuals
c          (e(.,ix)-transp)*x - f(ix)  .
c     the residuals are stored in the array  res.
c     indx  is rearranged so that the zero residuals
c     correspond to  indx(1),...,indx(iaddc)  .
c     ***************
c
c     +++++++++++++++
c     system routines  dabs,idint,dfloat,dsqrt
c
c     blas  sdot1
c
c     eps  is the smallest positive number which
c     satisfies   (1.0 + eps) .gt. 1.0   in the
c     precision of the arithmetic being used.
c     (alternatively, for less strict zero checking,
c      eps  can be set to a user-specified tolerance.)
c     +++++++++++++++
c
      integer i,iadp1,idummy,irand,ix,j,nactp1
      double precision eps,prod,temp,test,tol,zero
c
      double precision sdot1,dunif01
c
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
ctol = eps*dsqrt(dfloat(nvars))
      tol = eps
      nactp1 = nact+1
      if(.not.(1.le.iaddc))goto 23274
c
c     ***************
c     zero out all residuals known to be zero.
c     ***************
c
         do 23276 i = 1,iaddc 
            ix = indx(i)
            res(ix) = zero
23276       continue
c
c     ***************
c     compute the remaining residuals.
c     detect any more residuals which
c     are computationally zero, and
c     set them exactly zero.  their
c     associated indices are permuted
c     so that they are stored in
c     indx(nact+1),...,nact(iaddc).
c
c     (a fairly tight zero check is used.
c     it is far less expensive in running
c     time to neglect an extra zero
c     residual than to accept it and risk
c     invoking the anti-cycling
c     mechanisms in the program.
c     the accuracy of the solution as
c     finally determined is not affected.)
c     ***************
c
23274 continue
      iadp1 = iaddc+1
      if(.not.(iadp1.le.ncols))goto 23278
         do 23280 i = iadp1,ncols 
            ix = indx(i)
            temp = sdot1(nvars,e(1,ix),1,x,1)-f(ix)
            test = dabs(f(ix))
            do 23282 j = 1,nvars 
               prod = dabs(e(j,ix)*x(j))
               if(.not.(prod.gt.test))goto 23284
                  test = prod
23284          continue
23282          continue
            test = tol*test
            if(.not.(dabs(temp).gt.test))goto 23286
               res(ix) = temp
               goto 23287
c           else
23286          continue
               iaddc = iaddc+1
               indx(i) = indx(iaddc)
               indx(iaddc) = ix
               res(ix) = zero
23287       continue
23280       continue
c
c     ***************
c     if any new zero residuals have
c     been found, randomize their
c     ordering as an anti-cycling
c     device for  addcol.
c     ***************
c
23278 continue
      if(.not.(iaddc.gt.nactp1))goto 23288
         do 23290 i = nactp1,iaddc 
            irand = i+ifix(float(iaddc-i+1)*sngl(dunif01(0,idummy)))
            ix = indx(irand)
            indx(irand) = indx(i)
            indx(i) = ix
23290       continue
23288 continue
      return
      end
      subroutine daddcol(iaddc,idelc,nact,nvars,zz,nzzr,dd,rr,e,ner,
&      indx,w,eps)
c
      integer iaddc,idelc,indx(1),nact,ner,nvars,nzzr
      double precision dd(1),e(ner,1),rr(1)
      double precision w(1),zz(nzzr,1)
      double precision eps
c
c     ***************
c     cl1 version of addcol.
c
c     this routine administers the adjustment of the
c     z*d*r   decomposition for any new zero residuals.
c     the data corresponding to the zero residuals is indexed
c     in  indx(nact+1),...,indx(iaddc).
c     ***************
c
c     +++++++++++++++
c     blas  sasum1
c
c     eps  is the smallest positive number which
c     satisfies   (1.0 + eps) .gt. 1.0   in the
c     precision of the arithmetic being used.
c     (alternatively, for less strict zero checking,
c      eps  can be set to a user-specified tolerance.)
c     +++++++++++++++
c
      integer i,istrt,ix,nactp1,topx
      logical fail
      double precision colnrm,prjnrm
c
      double precision sasum1
c
c
c     /////////////////  begin program  //////////////////
c
      topx = nvars+1
      istrt = nact+1
      if(.not.(istrt.le.iaddc))goto 23292
c
c     ***************
c     candidates for addition to the  z*d*r
c     factorization are inspected in random
c     order to hinder cycling.
c     the randomization was carried out by  resid.
c
c     if a candidate has just been released
c     from the factorization or is dependent upon the
c     columns in the factorization,
c     then it is omitted from addition.
c
c     upon exit, indices of such omitted
c     columns are to be found in
c          indx(nact+1),...,indx(iaddc) .
c     ***************
c
         do 23294 i = istrt,iaddc 
            nactp1 = nact+1
            ix = indx(i)
            call dzdrpoc(nvars,nact,zz,nzzr,dd,e(1,ix),w,fail)
            colnrm = sasum1(nvars,e(1,ix),1)
            prjnrm = sasum1(nvars,w,1)
            if(.not.(prjnrm.gt.eps*colnrm.and.ix.ne.idelc))goto 23296
               indx(i) = indx(nactp1)
               indx(nactp1) = ix
               call dzdrcin(nvars,nact,zz,nzzr,dd,rr,e(1,ix),fail,w)
23296       continue
23294       continue
23292 continue
      return
      end
      subroutine drql1obj(iaddc,nact,nrq,nl1,nallq,ncols,nvars,e,ner,
&      res,grd,erql1n,pen,penpar,indx,theta,grd1)
c
      integer iaddc,indx(1),nact,nallq,nrq,nl1,ncols,ner,nvars,nrql1
      double precision e(ner,1),erql1n,grd(1),pen,penpar,res(1),theta,
&      grd1(1)
c
c     ***************
c     crql1 version of object.
c
c     this routine administers the evaluation of the
c     penalty (objective) function given the equation
c     and constraint residuals.  it also computes the
c     restricted gradient of the function.
c
c     columns which are not in the  z*d*r factorization
c     but which are associated with zero residuals must
c     be included in the restricted gradient with random
c     signs as an anti-cycling device.
c     the indices of these columns are to be
c     found in  indx(nact+1),...,indx(iaddc)
c     ***************
c
c     +++++++++++++++
c     system routines  dabs,dsign
c
c     blas  saxpy1,scopy1
c     +++++++++++++++
c
      integer i,idummy,ix,nactp1
      double precision half,one,two,three,tmp,zero,wgt,wgt1
c
      double precision dunif01
c
c     +++++++++++++++
c     the following declarations are necessary
c     for portability when  scopy1  is used, as
c     it is below, to fill arrays with a single
c     value  (zero=zip  in this case).
c     +++++++++++++++
c
      double precision zip(1)
      equivalence(zero,zip(1))
c
      data half/0.5d+00/
      data one/1.0d+00/
      data two/2.0d+00/
      data three/3.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      nrql1 = nrq+nl1
      nactp1 = nact+1
      erql1n = zero
      pen = zero
      call scopy1(nvars,zip,0,grd,1)
      call scopy1(nvars,zip,0,grd1,1)
      if(.not.(nactp1.le.ncols))goto 23298
         do 23300 i = nactp1,ncols 
            ix = indx(i)
            tmp = res(ix)
            wgt = dsign(one,tmp)
            wgt1 = wgt
            if(.not.(i.le.iaddc))goto 23302
               wgt1 = dunif01(0,idummy)
               if(.not.(wgt1 .lt. one/three))goto 23304
                  wgt = -one
                  goto 23305
c              else
23304             continue
                  if(.not.(wgt1 .gt. two/three))goto 23306
                     wgt = one
                     goto 23307
c                 else
23306                continue
                     wgt = zero
23307             continue
23305          continue
23302       continue
            if(.not.(ix.le.nrq))goto 23308
               wgt = one-two*theta+wgt
               wgt1 = wgt
23308       continue
            if(.not.(i.le.iaddc))goto 23310
               wgt1 = zero
23310       continue
            if(.not.(ix.le.nallq.or.wgt.le.zero))goto 23312
cwhy <= zero?
               if(.not.(wgt1 .ne. zero))goto 23314
                  if(.not.(ix.le.nrql1))goto 23316
                     erql1n = erql1n+tmp*wgt
                     tmp = tmp*penpar
23316             continue
                  pen = pen+tmp*wgt
23314          continue
               if(.not.(ix.le.nrql1))goto 23318
                  wgt = wgt*penpar
                  wgt1 = wgt1*penpar
23318          continue
               call saxpy1(nvars,wgt,e(1,ix),1,grd,1)
               call saxpy1(nvars,wgt1,e(1,ix),1,grd1,1)
23312       continue
23300       continue
23298 continue
      return
      end
      subroutine drql1gv(idelc,nact,nvars,nrq,nl1,nallq,e,ner,grd,coef,
&      penpar,indx,theta,eps)
c
      integer idelc,indx(1),nact,nallq,nrq,nl1,ner,nvars,nrql1
      double precision coef(1),e(ner,1),grd(1),penpar,theta
c
c     ***************
c     crql1  version.
c
c     set up the right-hand-side vector
c     (and store in the array  coef)
c     for the linear problem which determines
c     a descent direction  p  in the case where
c     the projection of the restricted gradient is zero.
c     ***************
c
c     +++++++++++++++
c     system routines  dabs,dfloat,idint,dsign
c
c     blas  saxpy1
c
c     eps  is the smallest positive number which
c     satisfies   (1.0 + eps) .gt. 1.0   in the
c     precision of the arithmetic being used.
c     (alternatively, for less strict zero checking,
c      eps  can be set to a user-specified tolerance.)
c     +++++++++++++++
c
      integer i,idummy,irand,ix
      double precision cf,eps,one,ope,s,tmp,tmpsav,zero,two
c
      double precision dunif01
c
      data one/1.0d+00/
      data two/2.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
c     ***************
c     find the most out-of-kilter
c     coefficient.  begin inspecting
c     the coefficients at a random index
c     to hinder cycling.  set  coef
c     to zero on the fly.
c     ***************
c
      nrql1 = nrq+nl1
      ope = one+eps
      idelc = 0
      tmpsav = zero
      if(.not.(1.le.nact))goto 23320
         irand = ifix(float(nact)*sngl(dunif01(0,idummy)))
         do 23322 i = 1,nact 
            irand = irand+1
            if(.not.(irand.gt.nact))goto 23324
               irand = 1
23324       continue
            ix = indx(irand)
            cf = coef(irand)
            coef(irand) = zero
            if(.not.(ix.gt.nallq))goto 23326
               tmp = cf+eps
               goto 23327
c           else
23326          continue
               if(.not.(ix.le.nrql1))goto 23328
                  cf = cf/penpar
23328          continue
               tmp = ope-dabs(cf)
               if(.not.(ix.le.nrq))goto 23330
                  tmp = tmp+dsign(one,cf)*(-one+two*theta)
23330          continue
23327       continue
            if(.not.(tmp.lt.tmpsav))goto 23332
c? what about w_nu >1
               idelc = irand
               s = dsign(one,cf)
               tmpsav = tmp
23332       continue
23322       continue
c
c     ***************
c     if no coefficients are out of kilter,
c     then return.  otherwise set a
c     value in an appropriate component
c     (indicated by  idelc)  of  coef
c     and adjust the restricted gradient
c     if necessary.
c     ***************
c
         if(.not.(idelc.ne.0))goto 23334
            coef(idelc) = -s
            ix = indx(idelc)
            if(.not.(ix.le.nallq))goto 23336
               tmp = -s
               if(.not.(ix.le.nrql1))goto 23338
                  tmp = tmp*penpar
23338          continue
               if(.not.(ix.le.nrq))goto 23340
                  tmp = tmp+(one-two*theta)*penpar
23340          continue
               call saxpy1(nvars,tmp,e(1,ix),1,grd,1)
23336       continue
23334    continue
23320 continue
      return
      end
c
c     ---------------
c     fourth level subroutines --
c               ddkheap,dunif01,dzdrcin,dzdrcou,
c               dzdrgit,dzdrgnv,dzdrpoc
c     ---------------
c
      subroutine ddkheap(make,ir,indx,aray)
c
      integer indx(1),ir
      logical make
      double precision aray(1)
c
c     ***************
c     an adaptation of d. e. knuth,s heaping
c     routines (see volume 3 of
c          the art of computer programming  ).
c     if  make  is  .true.,  the full heap building
c     process is carried out on
c          aray(1),...,aray(ir) ,
c     and the value of  ir  is unchanged.
c     if  make  is  .false.,  one step of the sorting
c     process is carried out to provide the next
c     element of  aray  in order,  and the variable
c     ir  is decreased by one.  the interruption of the
c     sorting phase is built in via the flag  once.
c     indx  is an index vector associated with
c     aray  which must be rearranged in parallel
c     with it.
c     ***************
c
      integer i,il,it,j
      logical once
      double precision t
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(ir.gt.1))goto 23342
c
c     ***************
c     test whether or not the initial
c     heap is to be built
c     ***************
c
         il = 1
         if(.not.(make))goto 23344
            il = (ir/2)+1
23344    continue
         once = .false.
c        repeat
23346       continue
            if(.not.(il.gt.1))goto 23349
c
c     ***************
c     the heap-building phase uses this branch
c     ***************
c
               il = il-1
               it = indx(il)
               t = aray(il)
               goto 23350
c           else
23349          continue
c
c     ***************
c     the sorting phase uses this branch
c     ***************
c
               if(.not.(make.or.once))goto 23351
                  return
23351          continue
               once = .true.
               it = indx(ir)
               t = aray(ir)
               indx(ir) = indx(1)
               aray(ir) = aray(1)
               ir = ir-1
               if(.not.(ir.le.1))goto 23353
                  goto 23348
23353          continue
23350       continue
c
c     ***************
c     the remaining statements are common
c     to both phases and embody the
c     heap-rectifying (sifting) section
c     ***************
c
            j = il
c           repeat
23355          continue
               i = j
               j = 2*j
               if(.not.(j.lt.ir))goto 23358
                  if(.not.(aray(j).gt.aray(j+1)))goto 23360
                     j = j+1
23360             continue
                  goto 23359
c              else
23358             continue
                  if(.not.(j.ne.ir))goto 23362
                     goto 23357
23362             continue
23359          continue
               if(.not.(t.le.aray(j)))goto 23364
                  goto 23357
23364          continue
               indx(i) = indx(j)
               aray(i) = aray(j)
23356          goto 23355
23357       continue
            indx(i) = it
            aray(i) = t
23347       goto 23346
23348    continue
         indx(1) = it
         aray(1) = t
         goto 23343
c     else
23342    continue
         if(.not.(.not.make))goto 23366
            ir = 0
23366    continue
23343 continue
      return
      end
      double precision function dunif01(iseed,ix)
c
      integer iseed,ix,ix0
c
      data ix0/2/
c
c     +++++++++++++++
c     system routines  dfloat,mod
c     +++++++++++++++
c
c     --------------------------------------------------------------
c     --------------------------------------------------------------
c
c     *****purpose-
c     this function returns a pseudo-random number distributed
c     uniformly in the interval (0,1).
c
c     *****parameter description-
c     on input-
c
c     iseed,  if it is nonzero modulo 9973, becomes the
c          new seed, i.e. it replaces the internally stored
c          value of ix0.  on machines where fortran variables
c          retain their values between calls, the internally
c          stored value if ix0 is the value assigned to  ix  in
c          the previous invocation of  dunif01.  otherwise -- and
c          in the first call to  dunif01 --  ix0=2.
c
c     on output-
c
c     ix is the next integer in a pseudo-random sequence of
c          integers between  1  and  9972  and is generated from its
c          predecessor  ix0  (i.e.  from  iseed,  if  iseed  is nonzero
c          modulo 9973).  ix  is the value which  iseed  should have
c          in the next invocation of  dunif01  to get the next
c          pseudo-random number.  the caller will often pass the
c          same variable for  iseed  as for  ix,
c          e.g.  x = dunif01(ix,ix).
c
c     *****application and usage restrictions-
c     dunif01  should only be used when portability is important and a
c     course random number generator suffices.  applications requiring
c     a fine, high precison generator should use one with a much
c     larger modulus.
c
c     *****algorithm notes-
c     dunif01 will run on any machine having at least 20 bits of ac-
c     curacy for fixed-point arithmitic.  it is based on a generator
c     recommended in (3), which passes the spectral test with flying
c     colors -- see (1) and (2).
c
c     references-
c     (1) hoaglin, d.c. (1976), theoretical properties of congruential
c     random-number generators-  an empirical view,
c     memorandum ns-340, dept. of statistics, harvard univ.
c
c     (2) knuth, d.e. (1969), the art of computer programming, vol. 2
c     (seminumerical algorithms), addison-wesley, reading, mass.
c
c     (3) smith, c.s. (1971), multiplicative pseudo-random number
c     generators with prime modulus, j. assoc. comput. mach. 18,
c     pp. 586-593.
c
c     *****general-
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c     permission for the use of  dunif01  in  cl1  was
c     generously given by  v. klema  and  d. hoaglin.
c
c     --------------------------------------------------------------
c     --------------------------------------------------------------
c
      if(.not.(iseed.ne.0))goto 23368
         ix = mod(iseed,99730)
         if(.not.(ix.ne.0))goto 23370
            ix0 = ix
23370    continue
c  ***
c  in order that all fixed-point calculations require only 20 bit
c  arithmetic, we use two calls to  mod  to compute
c  ix0 = mod(3432*ix0, 9973).
c  ***
23368 continue
      ix0 = mod(52*mod(66*ix0,99730),99730)
      ix = ix0
      dunif01 = dfloat(ix0)/99730.0d+00
      return
      end
      subroutine dzdrcin(n,k,zz,nzzr,dd,rr,col,fail,w)
c
      integer k,n,nzzr
      logical fail
      double precision col(1),dd(1),rr(1),w(1),zz(nzzr,1)
c
c     ***************
c     prepared by richard bartels
c     the university of waterloo
c     computer science department
c     latest update .... 30 november, 1979.
c
c     given the factorization
c
c          zz*dd*rr
c
c     of some  n by k  matrix
c
c       (0 .le. k .lt. n)
c         (n .ge. 1),
c
c     where
c
c          (zz-transp)*(zz) = (dd-inv),
c          dd  is diagonal and nonsingular,
c     and
c          rr  has zeros below the diagonal,
c
c     and given a  (k+1)th  column
c     to be addedto the original matrix,
c     this program updates  zz,dd and rr.
c
c     the value of  k  is increased by one.
c
c     w  is a scratch array.
c
c     use is made of routines from the library
c     of basic linear algebra subroutines (blas).
c
c     parameters...
c
c                     input/
c       name   type   output/   sub-    description
c                     scratch  scripts
c       -------------------------------------------
c       n      int.      i              number of rows
c
c       k      int.     i/o             number of columns
c
c       zz     double precision     i/o       2     scaled orthogonal
c                                       matrix
c
c       nzzr   int.      i              row dimension of zz
c
c       dd     double precision     i/o       1     diagonal scaling
c                                       matrix (diagonal
c                                       elements only)
c
c       rr     double precision     i/o       1     right-triangular
c                                       matrix in compact form.
c
c       col    double precision      i        1     column to be
c                                       added to  rr
c
c       fail   log.      o             .true.  if  k,n
c                                       are improper
c
c       w      double precision     scr       1     workspace
c       -------------------------------------------
c
c     the  i-th  segment of the array  rr  is  n-i+2 spaces
c     long and contains  1  work space followed by the
c     k-i+1  elements of row  i  followed by  n-k
c     scratch spaces.
c     ***************
c
c     +++++++++++++++
c     blas  scopy1,sdot1,srotm1,srotmg1
c     +++++++++++++++
c
      integer i,j,jdel,kp1,kp2
      double precision di,one,param(5),wi,zero
c
      double precision sdot1
c
c     +++++++++++++++
c     the following declarations are necessary
c     for portability when  scopy1  is used, as
c     it is below, to fill arrays with a single value
c     (one=unity  and  zero=zip  in this case).
c     +++++++++++++++
c
      double precision unity(1),zip(1)
      equivalence(one,unity(1)),(zero,zip(1))
c
      data one/1.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(k.lt.0.or.k.ge.n.or.n.gt.nzzr))goto 23372
         fail = .true.
         goto 23373
c     else
23372    continue
         if(.not.(k.le.0))goto 23374
c
c     ***************
c     for the special case that the
c     factorization was vacuous,
c     reset the arrays  zz and dd
c     to represent the  identity.
c     ***************
c
            k = 0
            call scopy1(n,unity,0,dd,1)
            call scopy1(n*n,zip,0,zz,1)
            do 23376 i = 1,n
               zz(i,i) = one
23376          continue
23374    continue
         kp1 = k+1
         kp2 = k+2
c
c     ***************
c     transform the incoming column,
c     and store the result in  w.
c     ***************
c
         do 23378 i = 1,n
            w(i) = sdot1(n,zz(1,i),1,col,1)
23378       continue
c
c     ***************
c     zero out the spike which would result from
c     storing  w  in  rr.   update  zz  and  dd.
c     ***************
c
         if(.not.(kp2.le.n))goto 23380
            do 23382 i = kp2,n 
               di = dd(i)
               wi = w(i)
               call srotmg1(dd(kp1),di,w(kp1),wi,param)
               w(i) = wi
               dd(i) = di
               call srotm1(n,zz(1,kp1),1,zz(1,i),1,param)
23382          continue
c
c     ***************
c     store the new column, which is still
c     in the array  w,  into  rr.
c     ***************
c
23380    continue
         j = kp2
         jdel = n
         do 23384 i = 1,kp1 
            rr(j) = w(i)
            j = j+jdel
            jdel = jdel-1
23384       continue
         k = kp1
         fail = .false.
23373 continue
      return
      end
      subroutine dzdrcou(n,k,zz,nzzr,dd,rr,ic,fail)
c
      integer ic,k,n,nzzr
      logical fail
      double precision dd(1),rr(1),zz(nzzr,1)
c
c     ***************
c     prepared by richard bartels
c     the university of waterloo
c     computer science department
c     latest update .... 30 november, 1979.
c
c     given the factorization
c
c          zz*dd*rr
c
c     of some  n by k  matrix
c
c       (1 .le. k .le. n)
c          (n .ge. 1),
c
c     where
c
c          (zz-transp)*(zz) = (dd-inv),
c          dd  is diagonal and nonsingular,
c     and
c          rr  has zeros below the diagonal,
c
c     and given the index  ic  of a column
c     to be removed  (1 .le. ic .le. k),
c     this program updates  zz,dd and rr .
c
c     the value of  k  is decreased by one, and
c     the column ordering in  rr  is changed.
c
c     use is made of routines from the library
c     of basic linear algebra subroutines (blas).
c
c     parameters...
c
c       name   type    input/   sub-    description
c                      output  scripts
c       -------------------------------------------
c       n      int.      i              number of rows
c
c       k      int.     i/o             number of columns
c
c       zz     double precision     i/o       2     scaled orthogonal
c                                       matrix
c
c       nzzr   int.      i              row dimension of zz
c
c       dd     double precision     i/o       1     diagonal scaling
c                                       matrix (diagonal
c                                       elements only)
c
c       rr     double precision     i/o       1     right-triangular
c                                       matrix in compact form.
c
c       ic     int.      i              index of column
c                                       to be removed
c
c       fail   log.      o              .true.  if  k,n,ic
c                                       are improper
c       -------------------------------------------
c
c     the  i-th  segment of the array  rr  is  n-i+2 spaces
c     long and contains  1  work space followed by the
c     k-i+1  elements of row  i  followed by  n-k
c     scratch spaces.
c     ***************
c
c     +++++++++++++++
c     blas  srotm1,srotmg1
c     +++++++++++++++
c
      integer i,im1,j,jend,jinc,jstrt,km1,lstrt
      double precision di,param(5),rj
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(k.lt.1.or.k.gt.n.or.n.gt.nzzr))goto 23386
         fail = .true.
         goto 23387
c     else
23386    continue
         km1 = k-1
c
c     ***************
c     special cases are handled first.
c     1.  k=1 and the factorization becomes null.
c     2.  ic=k and the updating is trivial.
c     ***************
c
         if(.not.(k.le.1))goto 23388
            k = 0
            fail = .false.
            goto 23389
c        else
23388       continue
            if(.not.(ic.ge.k))goto 23390
               k = km1
               fail = .false.
               goto 23391
c           else
23390          continue
c
c     ***************
c     general updating step.
c     the column to be deleted must be permuted
c     to the right, and subdiagonal elements
c     which result in  rr  have to be
c     transformed to zero.
c     ***************
c
               jstrt = ic+1
               jend = k
               jinc = n
               do 23392 i = 1,k 
c
c     ***************
c     permutation of the  i-th  row of rr.
c     ***************
c
                  do 23394 j = jstrt,jend
                     rr(j) = rr(j+1)
23394                continue
                  if(.not.(i.gt.ic))goto 23396
c
c     ***************
c     transformation of the current and last
c     rows  (i and i-1)  of rr  as well as
c     corresponding changes to  zz and dd.
c
c     the extra variables  di  and  rj
c     are used to avoid an error message
c     from the  pfort verifier, and they
c     may be removed, if desired, so that
c     the call to  srotmg1  would be
c
c     call srotmg1(dd(im1),dd(i),rr(lstrt),rr(jstrt),param)
c
c     ***************
c
                     im1 = i-1
                     di = dd(i)
                     rj = rr(jstrt)
                     call srotmg1(dd(im1),di,rr(lstrt),rj,param)
                     rr(jstrt) = rj
                     dd(i) = di
                     call srotm1(jend-jstrt+1,rr(lstrt+1),1,rr(jstrt+1),
&                     1,param)
                     call srotm1(n,zz(1,im1),1,zz(1,i),1,param)
                     jstrt = jstrt+1
c
c     ***************
c     index updating
c     ***************
c
23396             continue
                  lstrt = jstrt
                  jstrt = jstrt+jinc
                  jend = jend+jinc
                  jinc = jinc-1
23392             continue
               k = km1
               fail = .false.
23391       continue
23389    continue
23387 continue
      return
      end
      subroutine dzdrgit(n,k,zz,nzzr,rr,gv,sol,fail,w,big,eps)
c
      integer k,n,nzzr
      logical fail
      double precision gv(1),rr(1),w(1),sol(1),zz(nzzr,1)
c
c     ***************
c     prepared by richard bartels
c     the university of waterloo
c     computer science department
c     latest update .... 30 november, 1979.
c
c     given the factorization
c
c          zz*dd*rr
c
c     of some  n by k  matrix
c
c       (1 .le. k .le. n)
c         (n .ge. 1),
c
c     where
c
c          (zz-transp)*(zz) = (dd-inv),
c          dd  is diagonal and nonsingular,
c     and
c          rr  has zeros below the diagonal,
c
c     and given an arbitrary vector  gv  of
c     appropriate dimension, this routine finds the
c     vector  sol  satisfying the underdetermined system
c
c          (zz*dd*rr-transp.)*(sol) = (gv).
c
c     that is,
c
c          (sol) = ((zz*dd*rr)-gen.inv.-transp.)*(gv).
c
c     the array  dd  is not needed by  dzdrgit.
c
c     use is made of routines from the library
c     of basic linear algebra subroutines (blas).
c
c     w  is a scratch array.
c
c     parameters...
c
c                     input/
c       name   type   output/   sub-    description
c                     scratch  scripts
c       -------------------------------------------
c       n      int.      i              number of rows
c
c       k      int.     i/o             number of columns
c
c       zz     double precision     i/o       2     scaled orthogonal
c                                       matrix
c
c       nzzr   int.      i              row dimension of zz
c
c       rr     double precision     i/o       1     right-triangular
c                                       matrix in compact form.
c
c       gv     double precision      i        1     given vector
c
c       sol    double precision      o        1     solution
c
c       fail   log.      o              .true. if  n,k
c                                       are improper, or if
c                                       rr  is singular
c
c       w      double precision     scr       1     workspace
c       -------------------------------------------
c
c     the  i-th  segment of the array  rr  is  n-i+2 spaces
c     long and contains  1  work space followed by the
c     k-i+1  elements of row  i  followed by  n-k
c     scratch spaces.
c
c     if  gv  and  sol  are dimensioned to the
c     maximum of  n  and  k , then the same
c     storage array may be used for both of
c     these vectors.
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c
c     blas  saxpy1,scopy1
c
c     big  is the largest positive number
c     which can be represented in the
c     precision of the arithmetic being used.
c     +++++++++++++++
c
      integer i,j,jdel
      double precision big,wi,one,rrj,zero,eps
c
c     +++++++++++++++
c     the following declarations are necessary
c     for portability when  scopy1  is used, as
c     it is below, to fill arrays with a single value
c     (zero=zip  in this case).
c     +++++++++++++++
c
      double precision zip(1)
      equivalence(zero,zip(1))
c
      data one/1.0d+00/
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(k.lt.1.or.k.gt.n.or.n.gt.nzzr))goto 23398
         fail = .true.
         goto 23399
c     else
23398    continue
c
c     ***************
c     first solve  (rr-transp.)*(w) = (gv)
c     ***************
c
         call scopy1(k,gv,1,w,1)
         j = 2
         jdel = n+1
         do 23400 i = 1,k 
            rrj = rr(j)
            wi = w(i)
c* Here the check for ill-condition is NOT changed to use eps instead of big
c		if (dabs(rrj)<one)
c			if (eps*dabs(wi)>=dabs(rrj))
            if(.not.(dabs(rrj).lt.one))goto 23402
               if(.not.(dabs(wi).ge.dabs(rrj)*big))goto 23404
                  go to 150
23404          continue
23402       continue
            w(i) = wi/rrj
            if(.not.(i.lt.k))goto 23406
               call saxpy1(k-i,(-w(i)),rr(j+1),1,w(i+1),1)
23406       continue
            j = j+jdel
            jdel = jdel-1
23400       continue
c
c     ***************
c     now  (sol) = (zz)*(w)
c     ***************
c
         call scopy1(n,zip,0,sol,1)
         do 23408 i = 1,k
            call saxpy1(n,w(i),zz(1,i),1,sol,1)
23408       continue
         fail = .false.
         return
150      fail = .true.
23399 continue
      return
      end
      subroutine dzdrgnv(n,k,zz,nzzr,rr,gv,sol,fail,big)
c
      integer k,n,nzzr
      logical fail
      double precision gv(1),rr(1),sol(1),zz(nzzr,1)
c
c     ***************
c     prepared by richard bartels
c     the university of waterloo
c     computer science department
c     latest update .... 30 november, 1979.
c
c     given the factorization
c
c          zz*dd*rr
c
c     of some  n by k  matrix
c
c       (1 .le. k .le. n)
c         (n .ge. 1),
c
c     where
c
c          (zz-transp)*(zz) = (dd-inv),
c          dd  is diagonal and nonsingular,
c     and
c          rr  has zeros below the diagonal,
c
c     and given an arbitrary vector  gv  of
c     appropriate dimension, this routine finds the
c     vector  sol  given by
c
c          (sol) = ((zz*dd*rr)-gen.inv.)*(gv),
c
c     which represents the least squares problem
c
c        (zz*dd*rr)*(sol) = (gv).
c
c     the array  dd  is not needed by  dzdrgnv.
c
c     use is made of routines from the library
c     of basic linear algebra subroutines (blas).
c
c     parameters...
c
c       name   type   input/    sub-    description
c                     output/  scripts
c       -------------------------------------------
c       n      int.      i              number of rows
c
c       k      int.     i/o             number of columns
c
c       zz     double precision     i/o       2     scaled orthogonal
c                                       matrix
c
c       nzzr   int.      i              row dimension of zz
c
c       rr     double precision     i/o       1     right-triangular
c                                       matrix in compact form.
c
c       gv     double precision      i        1     given vector
c
c       sol    double precision      o        1     solution
c
c       fail   log.      o              .true. if  n,k
c                                       are improper, or if
c                                       rr  is singular
c       -------------------------------------------
c
c     the  i-th  segment of the array  rr  is  n-i+2 spaces
c     long and contains  1  work space followed by the
c     k-i+1  elements of row  i  followed by  n-k
c     scratch spaces.
c     ***************
c
c     +++++++++++++++
c     system routines  dabs
c
c     blas  sdot1
c
c     big  is the largest positive number
c     which can be represented in the
c     precision of the arithmetic being used.
c     +++++++++++++++
c
      integer i,ix,j,jdel
      double precision big,one,tden,tnum
c
      double precision sdot1
c
      data one/1.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      if(.not.(k.lt.1.or.k.gt.n.or.n.gt.nzzr))goto 23410
         fail = .true.
         goto 23411
c     else
23410    continue
c
c     ***************
c     form   (v) = (zz(1)-transp)*(gv),   where  zz(1)
c     is the matrix of the first  k  columns of  zz
c
c     v  can be stored in the array  sol.
c     ***************
c
         do 23412 i = 1,k
            sol(i) = sdot1(n,zz(1,i),1,gv,1)
23412       continue
c
c     ***************
c     backsolve the system
c          (rr)*(sol) = (v)
c     for the vector  sol
c
c     note that  sol  and  v
c     are stored in the same array.
c     ***************
c
         j = (((n+1)*(n+2)-(n-k+3)*(n-k+2))/2)+2
         jdel = n-k+3
         do 23414 ix = 1,k 
            i = k-ix+1
            tden = rr(j)
            tnum = sol(i)
            if(.not.(ix.gt.1))goto 23416
               tnum = tnum-sdot1(ix-1,rr(j+1),1,sol(i+1),1)
23416       continue
            if(.not.(dabs(tden).lt.one))goto 23418
               if(.not.(dabs(tnum).ge.dabs(tden)*big))goto 23420
                  go to 160
23420          continue
23418       continue
            sol(i) = tnum/tden
            j = j-jdel
            jdel = jdel+1
23414       continue
         fail = .false.
         return
160      fail = .true.
23411 continue
      return
      end
      subroutine dzdrpoc(n,k,zz,nzzr,dd,gv,poc,fail)
c
      integer k,n,nzzr
      logical fail
      double precision dd(1),poc(1),gv(1),zz(nzzr,1)
c
c     ***************
c     prepared by richard bartels
c     the university of waterloo
c     computer science department
c     latest update .... 30 november, 1979.
c
c     zz is an  n by n  (n .ge. 1)  scaled
c     orthogonal matrix.  dd  contains the
c     diagonal elements of a diagonal scaling
c     matrix.  gv  is a given vector of length  n.
c
c     we have
c
c          (zz-transp.)*(zz) = (dd-inv.)
c
c     and
c
c               zz*dd*rr = mat
c
c     for some  n by k  (0 .le. k .le. n)
c     matrix  rr  with zeros below the diagonal
c     and some given matrix  mat.  (niether  rr
c     nor  mat  are needed by  dzdrpoc.)
c
c     then
c
c    (proj(oc)) = (zz(2))*(dd(2))*(zz(2)-transp.)
c
c     is the (orthogonal) projector on the
c     complement of the range space of  mat,
c     where  zz(2)  represents the last  n-k
c     columns of  zz  and  dd(2)  represents the
c     lower-right-hand  n-k  order submatrix of  dd.
c
c     dzdrpoc  produces the vector
c
c               poc = (proj(oc))*gv .
c
c     use is made of routines from the library
c     of basic linear algebra subroutines (blas).
c
c     parameters...
c
c                     input/
c       name   type   output/   sub-    description
c                     scratch  scripts
c       -------------------------------------------
c       n      int.      i              order of  zz,dd
c
c       k      int.     i/o             number of columns
c                                       of  zz  defining
c                                       range of  mat
c
c       zz     double precision     i/o       2     scaled orthogonal
c                                       matrix
c
c       nzzr   int.      i              row dimension of zz
c
c       dd     double precision     i/o       1     diagonal scaling
c                                       matrix (diagonal
c                                       elements only)
c
c       gv     double precision      i        1     vector to be projected
c
c       poc    double precision      o        1     projection
c
c       fail   log.      o              .true.  if  n,k
c                                       are improper
c
c       -------------------------------------------
c
c     ***************
c
c     +++++++++++++++
c     blas  saxpy1,scopy1,sdot1
c     +++++++++++++++
c
      integer i,kp1
      double precision wi,zero
c
      double precision sdot1
c
c     +++++++++++++++
c     the following declarations are necessary
c     for portability when  scopy1  is used, as
c     it is below, to fill arrays with a single value
c     (zero=zip  in this case).
c     +++++++++++++++
c
      double precision zip(1)
      equivalence(zero,zip(1))
c
      data zero/0.0d+00/
c
c     /////////////////  begin program  //////////////////
c
      kp1 = k+1
      if(.not.(k.lt.0.or.k.gt.n.or.n.lt.1.or.n.gt.nzzr))goto 23422
         fail = .true.
         goto 23423
c     else
23422    continue
         if(.not.(k.le.0))goto 23424
c
c     ***************
c     case 1 ... zz(2)=zz  (k=0)
c     ***************
c
            call scopy1(n,gv,1,poc,1)
            fail = .false.
            goto 23425
c        else
23424       continue
            if(.not.(k.ge.n))goto 23426
c
c     ***************
c     case 2 ... zz(2) is vacuous  (k=n)
c     ***************
c
               call scopy1(n,zip,0,poc,1)
               fail = .false.
               goto 23427
c           else
23426          continue
c
c     ***************
c     case 3 ... zz(2)  is intermediate
c     between the other two cases
c     (0 .lt. k .lt. n)
c     ***************
c
               call scopy1(n,zip,0,poc,1)
               do 23428 i = kp1,n 
                  wi = sdot1(n,zz(1,i),1,gv,1)*dd(i)
                  call saxpy1(n,wi,zz(1,i),1,poc,1)
23428             continue
               fail = .false.
23427       continue
23425    continue
23423 continue
      return
      end

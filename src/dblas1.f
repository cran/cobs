C--- FIXME [MM] : The BLAS routine names are "s...." for single precision
C--- =====        but they are *coded* as double precision here and in drqssbc.f
C--- this is *not* helpful ==> fix mostly in ./drqssbc.f
C            ----------------------------
C===> Now removed all these unneeded single precision versions;
C===> Original still in ./Orig/dblas1.f.~1~

C--- Problem : BLAS (Linux :  nm -g /usr/lib/libblas.so)
C    =======
C    only exports  drot() and drotg() {and srot[g], c*, }.
C    their argument list looks quite different than the
C    srotm1() and srotmg1() here !
C
C ==> Keep only  srotm1() and srotmg1()  {which MM thinks are NOT quite BLAS!}
C                ========     =========
C-------------------------------------------------------------------------

c                  ** appendix **
c
c     basic linear algebra subroutines (blas) --
c               sasum1,saxpy1,scopy1,sdot1,
c               srotm1,srotmg1,sscal1
c
c     the  blas  are described in
c
c              basic linear algebra subprograms
c                    for fortran usage
c
c                           by
c              c. l. lawson, r. j. hanson,
c              d. r. kincaid, f. t. krogh.
c              acm trans. on math software
c              vol. 5, no. 3, 308-325 (1979)
c
c              the  blas  are taken as
c              primitives.  their code is not
c              part of the definition of  cl1.
c              they are included for reference
c              only.  please refer to the above
c              publication for current information.
c
c     if the  blas  are available independantly
c     on the target machine, this appendix should
c     be removed from the code for  cl1.
c     ---------------------------------------------
c

      subroutine srotm1 (n,sx,incx,sy,incy,sparam)
c
      integer incx,incy,n
      double precision sparam(5),sx(1),sy(1)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c
c     (sx(1)     sx(n))
c     (      ...      )
c     (sy(1)     sy(n))
c
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c
      integer i,kx,ky,nsteps
      double precision sflag,sh11,sh12,sh21,sh22,two,w,z,zero
c
      data zero,two /0.0d+00,2.0d+00/
c
      sflag=sparam(1)
      if(n .le. 0 .or.(sflag+two.eq.zero)) go to 140
          if(.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if(sflag) 50,10,30
   10          continue
               sh12=sparam(4)
               sh21=sparam(3)
                    do 20 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w+z*sh12
                    sy(i)=w*sh21+z
   20               continue
               go to 140
   30          continue
               sh11=sparam(2)
               sh22=sparam(5)
                    do 40 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z
                    sy(i)=-w+sh22*z
   40               continue
               go to 140
   50          continue
               sh11=sparam(2)
               sh12=sparam(4)
               sh21=sparam(3)
               sh22=sparam(5)
                    do 60 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z*sh12
                    sy(i)=w*sh21+z*sh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if(incx .lt. 0) kx=1+(1-n)*incx
          if(incy .lt. 0) ky=1+(1-n)*incy
c
          if(sflag)120,80,100
   80     continue
          sh12=sparam(4)
          sh21=sparam(3)
               do 90 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w+z*sh12
               sy(ky)=w*sh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          sh11=sparam(2)
          sh22=sparam(5)
               do 110 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z
               sy(ky)=-w+sh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          sh11=sparam(2)
          sh12=sparam(4)
          sh21=sparam(3)
          sh22=sparam(5)
               do 130 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z*sh12
               sy(ky)=w*sh21+z*sh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
          end

      subroutine srotmg1 (sd1,sd2,sx1,sy1,sparam)
c
      double precision sd1,sd2,sparam(5),sx1,sy1
c
c     construct the modified givens transformation matrix h
c     which zeros the second component of the 2-vector
c         (sqrt(sd1)*sx1,sqrt(sd2)*sy2)**t.
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c
c     +++++++++++++++
c     system routines  dabs
c     +++++++++++++++
c
      integer iflag,igo
      double precision gam,gamsq,one,rgam,rgamsq,sflag,sh11,sh12
      double precision sh21,sh22,sp1,sp2,sq1,sq2,stemp,su,two,zero
c
      data zero,one,two /0.0d+00,1.0d+00,2.0d+00/,iflag/1/
      data gam/4096.0d+00/
      data gamsq /1.678d+07/
      data rgam /2.441e-04/
      data rgamsq /5.960e-08/
c
      if(.not. sd1 .lt. zero) go to 10
c       go zero-h-d-and-sx1..
          go to 60
   10 continue
c     case-sd1-nonnegative
      sp2=sd2*sy1
      if(.not. sp2 .eq. zero) go to 20
          sflag=-two
          go to 260
c     regular-case..
   20 continue
      sp1=sd1*sx1
      sq2=sp2*sy1
      sq1=sp1*sx1
c
      if(.not. dabs(sq1) .gt. dabs(sq2)) go to 40
          sh21=-sy1/sx1
          sh12=sp2/sp1
c
          su=one-sh12*sh21
c
          if(.not. su .le. zero) go to 30
c         go zero-h-d-and-sx1..
               go to 60
   30     continue
               sflag=zero
               sd1=sd1/su
               sd2=sd2/su
               sx1=sx1*su
c         go scale-check..
               go to 100
   40 continue
          if(.not. sq2 .lt. zero) go to 50
c         go zero-h-d-and-sx1..
               go to 60
   50     continue
               sflag=one
               sh11=sp1/sp2
               sh22=sx1/sy1
               su=one+sh11*sh22
               stemp=sd2/su
               sd2=sd1/su
               sd1=stemp
               sx1=sy1*su
c         go scale-check
               go to 100
c     procedure..zero-h-d-and-sx1..
   60 continue
          sflag=-one
          sh11=zero
          sh12=zero
          sh21=zero
          sh22=zero
c
          sd1=zero
          sd2=zero
          sx1=zero
c         return..
          go to 220
c     procedure..fix-h..
   70 continue
      if(.not. sflag .ge. zero) go to 90
c
          if(.not. sflag .eq. zero) go to 80
          sh11=one
          sh22=one
          sflag=-one
          go to 90
   80     continue
          sh21=-one
          sh12=one
          sflag=-one
   90 continue
      go to igo,(120,150,180,210)
c     procedure..scale-check
  100 continue
  110     continue
          if(.not. sd1 .le. rgamsq) go to 130
               if(.not. iflag.eq.1) go to 105
c
c                   recompute rescaling parameters
c                   more accurately..
c
                    rgam = one/gam
                    gamsq = gam**2
                    rgamsq = rgam**2
                    iflag = 2
  105          continue
               if(sd1 .eq. zero) go to 160
               assign 120 to igo
c              fix-h..
               go to 70
  120          continue
               sd1=sd1*gamsq
               sx1=sx1*rgam
               sh11=sh11*rgam
               sh12=sh12*rgam
          go to 110
  130 continue
  140     continue
          if(.not. sd1 .ge. gamsq) go to 160
               assign 150 to igo
c              fix-h..
               go to 70
  150          continue
               sd1=sd1*rgamsq
               sx1=sx1*gam
               sh11=sh11*gam
               sh12=sh12*gam
          go to 140
  160 continue
  170     continue
          if(.not. dabs(sd2) .le. rgamsq) go to 190
               if(sd2 .eq. zero) go to 220
               assign 180 to igo
c              fix-h..
               go to 70
  180          continue
               sd2=sd2*gamsq
               sh21=sh21*rgam
               sh22=sh22*rgam
          go to 170
  190 continue
  200     continue
          if(.not. dabs(sd2) .ge. gamsq) go to 220
               assign 210 to igo
c              fix-h..
               go to 70
  210          continue
               sd2=sd2*rgamsq
               sh21=sh21*gam
               sh22=sh22*gam
          go to 200
  220 continue
          if(sflag)250,230,240
  230     continue
               sparam(3)=sh21
               sparam(4)=sh12
               go to 260
  240     continue
               sparam(2)=sh11
               sparam(5)=sh22
               go to 260
  250     continue
               sparam(2)=sh11
               sparam(3)=sh21
               sparam(4)=sh12
               sparam(5)=sh22
  260 continue
          sparam(1)=sflag
          return
      end

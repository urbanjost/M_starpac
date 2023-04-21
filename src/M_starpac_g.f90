      module M_starpac_g
      implicit none
      private
      public  ::  ssifa,  dswap,  dtrdi,  dscal,   strco
      public  ::  dasum,  dnrm2,  dsifa,  sscal,   isamax
      public  ::  ddot,   daxpy,  sdot,   idamax,  sswap
      public  ::  dsidi,  dtrco,  ssidi,  strdi,   scopy
      public  ::  saxpy,  sasum,  dcopy,  snrm2

      public  ::  r9lgic,  gamma,   dlngam,  erfc,    eprint
      public  ::  dgamr,   dlbeta,  dcsevl,  gamlim,  r9lgmc
      public  ::  alngam,  xerprt,  inits,   csevl,   xgetua
      public  ::  e9rint,  xerror,  d9gmit,  xsetf,   r9gmit
      public  ::  dlnrel,  gami,    d9lgic,  derf,    xerrwv
      public  ::  d9lgit,  seterr,  fdump,   r9lgit,  xerctl
      public  ::  s88fmt,  erf,     xerabt,  alnrel,  initds
      public  ::  dbetai,  i8save,  dgamma,  derfc,   dgamlm
      public  ::  albeta,  xersav,  dgamit,  gamit,   dlgams
      public  ::  gamr,    xerclr,  betai,   j4save,  xgetf
      public  ::  d9lgmc,  algams,  dgami

      public  ::  i1mach,  r1mach,  d1mach

      public  ::  setiv,   fixprt,  backop,  iprint,  fftlen
      public  ::  ecvf,    amehdr,  llhdrp,  eiagep,  eiveq
      public  ::  acfer,   nchose,  inperl,  eisii,   eiage
      public  ::  nlskl,   bfser,   aov1hd,  eiseq,   pltsym
      public  ::  pline,   modsum,  dckhdr,  nlerr,   eivii
      public  ::  stkrel,  stphdr,  corrhd,  eiveo,   ccfer
      public  ::  acfdtl,  eriodd,  amfer,   nlhdra,  ehdr
      public  ::  icnti,   llhdrg,  erdf,    stkst,   eisge
      public  ::  geni,    eisle,   stkclr,  lstlag,  eisrng
      public  ::  setesl,  enfft,   ufser,   factor,  msgx
      public  ::  cpyvii,  stkget,  setlag,  ldscmp,  prtcnt
      public  ::  nlhdrn,  amfhdr,  icopy,   stkset,  versp

      public  ::  imdcon,  stopx,  ufparm
      integer,save,public :: G_IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF G_IERR .EQ. 0, NO ERRORS WERE DETECTED
!        IF G_IERR .EQ. 1, ERRORS HAVE BEEN DETECTED

      contains

!SSIFA
      subroutine ssifa(a,lda,n,kpvt,info)
!
!     LATEST REVISION  -  JANUARY 24, 1990
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer info,lda,n
!
!  ARRAY ARGUMENTS
      real a(lda,*)
      integer kpvt(*)
!
!  LOCAL SCALARS
     real absakk,ak,akm1,alpha,bk,bkm1,colmax,denom,mulk,mulkm1,rowmax,&
     &   t
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep
      logical swap
!
!  EXTERNAL FUNCTIONS
!      INTEGER ISAMAX
!       EXTERNAL ISAMAX
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL SAXPY,SSWAP
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,amax1,sqrt
!
!
!     SSIFA FACTORS A REAL SYMMETRIC MATRIX BY ELIMINATION
!     WITH SYMMETRIC PIVOTING.
!
!     TO SOLVE  A*X = B , FOLLOW SSIFA BY SSISL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SSIFA BY SSISL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SSIFA BY SSIDI.
!     TO COMPUTE  INERTIA(A) , FOLLOW SSIFA BY SSIDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW SSIFA BY SSIDI.
!
!     ON ENTRY
!
!        A       REAL(LDA,N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
!                     BUT IT DOES INDICATE THAT SSISL OR SSIDI MAY
!                     DIVIDE BY ZERO IF CALLED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SSWAP,ISAMAX
!     FORTRAN ABS,AMAX1,SQRT
!
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
!
      info = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      k = n
   10 continue
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0e0) info = 1
!     ......EXIT
            go to 200
   20    continue
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         km1 = k - 1
         absakk = abs(a(k,k))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         imax = isamax(k-1,a(1,k),1)
         colmax = abs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            rowmax = 0.0e0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,abs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = isamax(imax-1,a(1,imax),1)
               rowmax = amax1(rowmax,abs(a(jmax,imax)))
   50       continue
            if (abs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
!
!           1 X 1 PIVOT BLOCK.
!
            if (.not.swap) go to 120
!
!              PERFORM AN INTERCHANGE.
!
               call sswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
!
!           PERFORM THE ELIMINATION.
!
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call saxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
!
!           SET THE PIVOT ARRAY.
!
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
!
!           2 X 2 PIVOT BLOCK.
!
            if (.not.swap) go to 160
!
!              PERFORM AN INTERCHANGE.
!
               call sswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
!
!           PERFORM THE ELIMINATION.
!
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0e0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call saxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call saxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
!
!           SET THE PIVOT ARRAY.
!
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
!DSWAP
      subroutine dswap(n,dx,incx,dy,incy)
!
!     INTERCHANGE DOUBLE PRECISION DX AND DOUBLE PRECISION DY.
!     FOR I = 0 TO N-1, INTERCHANGE  DX(LX+I*INCX) AND DY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      double precision dx(*),dy(*)
!
!  LOCAL SCALARS
      double precision dtemp1,dtemp2,dtemp3
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp1 = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp1
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp1 = dx(i)
        dtemp2 = dx(i+1)
        dtemp3 = dx(i+2)
        dx(i) = dy(i)
        dx(i+1) = dy(i+1)
        dx(i+2) = dy(i+2)
        dy(i) = dtemp1
        dy(i+1) = dtemp2
        dy(i+2) = dtemp3
   50 continue
      return
   60 continue
!
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
      ns = n*incx
        do 70 i=1,ns,incx
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   70   continue
      return
      end
!DTRDI
      subroutine dtrdi(t,ldt,n,det,job,info)
!
!     DTRDI COMPUTES THE DETERMINANT AND INVERSE OF A DOUBLE PRECISION
!     TRIANGULAR MATRIX.
!
!     ON ENTRY
!
!        T       DOUBLE PRECISION(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!
!        JOB     INTEGER
!                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.
!                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.
!                = 100       DET, NO INVERSE.
!                = 110       DET, INVERSE OF LOWER TRIANGULAR.
!                = 111       DET, INVERSE OF UPPER TRIANGULAR.
!
!     ON RETURN
!
!        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     DOUBLE PRECISION(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!        INFO    INTEGER
!                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR
!                AND THE INVERSE IS REQUESTED.
!                OTHERWISE INFO CONTAINS THE INDEX OF
!                A ZERO DIAGONAL ELEMENT OF T.
!
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL
!     FORTRAN DABS,MOD
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer info,job,ldt,n
!
!  ARRAY ARGUMENTS
      double precision det(2),t(ldt,*)
!
!  LOCAL SCALARS
      double precision temp,ten
      integer i,j,k,kb,km1,kp1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DAXPY,DSCAL
!
!  INTRINSIC FUNCTIONS
      intrinsic dabs,mod
!
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 180
!
!        COMPUTE DETERMINANT
!
         if (job/100 .eq. 0) go to 70
            det(1) = 1.0d0
            det(2) = 0.0d0
            ten = 10.0d0
            do 50 i = 1, n
               det(1) = t(i,i)*det(1)
!           ...EXIT
               if (det(1) .eq. 0.0d0) go to 60
   10          if (dabs(det(1)) .ge. 1.0d0) go to 20
                  det(1) = ten*det(1)
                  det(2) = det(2) - 1.0d0
               go to 10
   20          continue
   30          if (dabs(det(1)) .lt. ten) go to 40
                  det(1) = det(1)/ten
                  det(2) = det(2) + 1.0d0
               go to 30
   40          continue
   50       continue
   60       continue
   70    continue
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
         if (mod(job/10,10) .eq. 0) go to 170
            if (mod(job,10) .eq. 0) go to 120
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  do 100 k = 1, n
                     info = k
!              ......EXIT
                     if (t(k,k) .eq. 0.0d0) go to 110
                     t(k,k) = 1.0d0/t(k,k)
                     temp = -t(k,k)
                     call dscal(k-1,temp,t(1,k),1)
                     kp1 = k + 1
                     if (n .lt. kp1) go to 90
                     do 80 j = kp1, n
                        temp = t(k,j)
                        t(k,j) = 0.0d0
                        call daxpy(k,temp,t(1,k),1,t(1,j),1)
   80                continue
   90                continue
  100             continue
                  info = 0
  110          continue
            go to 160
  120       continue
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
               do 150 kb = 1, n
                  k = n + 1 - kb
                  info = k
!     ............EXIT
                  if (t(k,k) .eq. 0.0d0) go to 180
                  t(k,k) = 1.0d0/t(k,k)
                  temp = -t(k,k)
                  if (k .ne. n) call dscal(n-k,temp,t(k+1,k),1)
                  km1 = k - 1
                  if (km1 .lt. 1) go to 140
                  do 130 j = 1, km1
                     temp = t(k,j)
                     t(k,j) = 0.0d0
                     call daxpy(n-k+1,temp,t(k,k),1,t(k,j),1)
  130             continue
  140             continue
  150          continue
               info = 0
  160       continue
  170    continue
  180 continue
      return
      end
!DSCAL
      subroutine dscal(n,da,dx,incx)
!
!     REPLACE DOUBLE PRECISION DX BY DOUBLE PRECISION DA*DX.
!     FOR I = 0 TO N-1, REPLACE DX(1+I*INCX) WITH  DA * DX(1+I*INCX)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision da
      integer incx,n
!
!  ARRAY ARGUMENTS
      double precision dx(*)
!
!  LOCAL SCALARS
      integer i,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      ns = n*incx
          do 10 i = 1,ns,incx
          dx(i) = da*dx(i)
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
!STRCO
      subroutine strco(t,ldt,n,rcond,z,job)
!***BEGIN PROLOGUE  STRCO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A3
!***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,TRIANGULAR
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  ESTIMATES THE CONDITION OF A REAL TRIANGULAR MATRIX.
!***DESCRIPTION
!     STRCO ESTIMATES THE CONDITION OF A REAL TRIANGULAR MATRIX.
!     ON ENTRY
!        T       REAL(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX.  THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!        JOB     INTEGER
!                = 0         T  IS LOWER TRIANGULAR.
!                = NONZERO   T  IS UPPER TRIANGULAR.
!     ON RETURN
!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .
!                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS
!                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  SASUM,SAXPY,SSCAL
!***END PROLOGUE  STRCO

!...SCALAR ARGUMENTS
      real rcond
      integer job,ldt,n

!...ARRAY ARGUMENTS
      real t(ldt,*),z(*)

!...LOCAL SCALARS
      real ek,s,sm,tnorm,w,wk,wkm,ynorm
      integer i1,j,j1,j2,k,kk,l
      logical lower

!...EXTERNAL FUNCTIONS
!      REAL,external :: SASUM
!
!...EXTERNAL SUBROUTINES
!       EXTERNAL SAXPY,SSCAL

!...INTRINSIC FUNCTIONS
      intrinsic abs,amax1,sign

!***FIRST EXECUTABLE STATEMENT  STRCO

      lower = job .eq. 0

!     COMPUTE 1-NORM OF T

      tnorm = 0.0e0
      do 10 j = 1, n
         l = j
         if (lower) l = n + 1 - j
         i1 = 1
         if (lower) i1 = j
         tnorm = amax1(tnorm,sasum(l,t(i1,j),1))
   10 continue

!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.

!     SOLVE TRANS(T)*Y = E

      ek = 1.0e0
      do 20 j = 1, n
         z(j) = 0.0e0
   20 continue
      do 100 kk = 1, n
         k = kk
         if (lower) k = n + 1 - kk
         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .le. abs(t(k,k))) go to 30
            s = abs(t(k,k))/abs(ek-z(k))
            call sscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (t(k,k) .eq. 0.0e0) go to 40
            wk = wk/t(k,k)
            wkm = wkm/t(k,k)
         go to 50
   40    continue
            wk = 1.0e0
            wkm = 1.0e0
   50    continue
         if (kk .eq. n) go to 90
            j1 = k + 1
            if (lower) j1 = 1
            j2 = n
            if (lower) j2 = k - 1
            do 60 j = j1, j2
               sm = sm + abs(z(j)+wkm*t(k,j))
               z(j) = z(j) + wk*t(k,j)
               s = s + abs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               w = wkm - wk
               wk = wkm
               do 70 j = j1, j2
                  z(j) = z(j) + w*t(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)

      ynorm = 1.0e0

!     SOLVE T*Z = Y

      do 130 kk = 1, n
         k = n + 1 - kk
         if (lower) k = kk
         if (abs(z(k)) .le. abs(t(k,k))) go to 110
            s = abs(t(k,k))/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
  110    continue
         if (t(k,k) .ne. 0.0e0) z(k) = z(k)/t(k,k)
         if (t(k,k) .eq. 0.0e0) z(k) = 1.0e0
         i1 = 1
         if (lower) i1 = k + 1
         if (kk .ge. n) go to 120
            w = -z(k)
            call saxpy(n-kk,w,t(i1,k),1,z(i1),1)
  120    continue
  130 continue
!     MAKE ZNORM = 1.0
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm

      if (tnorm .ne. 0.0e0) rcond = ynorm/tnorm
      if (tnorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
!DASUM
      double precision function dasum(n,dx,incx)
!***BEGIN PROLOGUE  DASUM
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A3A
!***KEYWORDS  ADD,BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAGNITUDE,SUM,
!             VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  SUM OF MAGNITUDES OF D.P. VECTOR COMPONENTS
!***DESCRIPTION
!                B L A S  SUBPROGRAM
!    DESCRIPTION OF PARAMETERS
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX
!     --OUTPUT--
!    DASUM  DOUBLE PRECISION RESULT (ZERO IF N .LE. 0)
!     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX.
!     DASUM = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX))
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DASUM

!...SCALAR ARGUMENTS
      integer incx,n

!...ARRAY ARGUMENTS
      double precision dx(*)

!...LOCAL SCALARS
      integer i,m,mp1,ns

!...INTRINSIC FUNCTIONS
      intrinsic dabs,mod

!***FIRST EXECUTABLE STATEMENT  DASUM

      dasum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      ns = n*incx
          do 10 i=1,ns,incx
          dasum = dasum + dabs(dx(i))
   10     continue
      return

!        CODE FOR INCREMENTS EQUAL TO 1.

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dasum = dasum + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))&
     &   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
   50 continue
      return
      end
!DNRM2
      double precision function dnrm2 ( n, dx, incx)
!
!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N .LE. 0 RETURN WITH RESULT = 0.
!     IF N .GE. 1 THEN INCX MUST BE .GE. 1
!
!           C.L.LAWSON, 1978 JAN 08
!
!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,n
!
!  ARRAY ARGUMENTS
      double precision dx(*)
!
!  LOCAL SCALARS
      double precision cuthi,cutlo,hitest,one,sum,xmax,zero
      integer i,j,next,nn
!
!  INTRINSIC FUNCTIONS
      intrinsic dabs,dsqrt,float
!
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
      data   zero, one /0.0d0, 1.0d0/
!
      xmax = zero
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
!
   10 assign 30 to next
      sum = zero
      nn = n * incx
!                                                 BEGIN MAIN LOOP
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        PHASE 1.  SUM IS ZERO
!
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
!
!                                PREPARE FOR PHASE 2.
      assign 70 to next
      go to 105
!
!                                PREPARE FOR PHASE 4.
!
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
!
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 sum = (sum * xmax) * xmax
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 hitest = cuthi/float( n )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
!DSIFA
      subroutine dsifa(a,lda,n,kpvt,info)
!
!     DSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION
!     WITH SYMMETRIC PIVOTING.
!
!     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.
!     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.
!
!     ON ENTRY
!
!        A       DOUBLE PRECISION(LDA,N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
!                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY
!                     DIVIDE BY ZERO IF CALLED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSWAP,IDAMAX
!     FORTRAN DABS,DMAX1,DSQRT
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer info,lda,n
!
!  ARRAY ARGUMENTS
      double precision a(lda,*)
      integer kpvt(*)
!
!  LOCAL SCALARS
     double precision absakk,ak,akm1,alpha,bk,bkm1,colmax,denom,mulk,&
     &   mulkm1,rowmax,t
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep
      logical swap
!
!  EXTERNAL FUNCTIONS
!      INTEGER IDAMAX
!       EXTERNAL IDAMAX
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DAXPY,DSWAP
!
!  INTRINSIC FUNCTIONS
      intrinsic dabs,dmax1,dsqrt
!
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
!
      info = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      k = n
   10 continue
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
!     ......EXIT
            go to 200
   20    continue
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         km1 = k - 1
         absakk = dabs(a(k,k))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
!
!           1 X 1 PIVOT BLOCK.
!
            if (.not.swap) go to 120
!
!              PERFORM AN INTERCHANGE.
!
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
!
!           PERFORM THE ELIMINATION.
!
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
!
!           SET THE PIVOT ARRAY.
!
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
!
!           2 X 2 PIVOT BLOCK.
!
            if (.not.swap) go to 160
!
!              PERFORM AN INTERCHANGE.
!
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
!
!           PERFORM THE ELIMINATION.
!
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call daxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
!
!           SET THE PIVOT ARRAY.
!
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
!SSCAL
      subroutine sscal(n,sa,sx,incx)
!
!     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
!     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real sa
      integer incx,n
!
!  ARRAY ARGUMENTS
      real sx(*)
!
!  LOCAL SCALARS
      integer i,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      ns = n*incx
          do 10 i = 1,ns,incx
          sx(i) = sa*sx(i)
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
!ISAMAX
      integer function isamax(n,sx,incx)
!
!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
!     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX))
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,n
!
!  ARRAY ARGUMENTS
      real sx(*)
!
!  LOCAL SCALARS
      real smax,xmag
      integer i,ii,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic abs
!
      isamax = 0
      if(n.le.0) return
      isamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      smax = abs(sx(1))
      ns = n*incx
      ii = 1
          do 10 i=1,ns,incx
          xmag = abs(sx(i))
          if(xmag.le.smax) go to 5
          isamax = ii
          smax = xmag
    5     ii = ii + 1
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 smax = abs(sx(1))
      do 30 i = 2,n
         xmag = abs(sx(i))
         if(xmag.le.smax) go to 30
         isamax = i
         smax = xmag
   30 continue
      return
      end
!DDOT
      double precision function ddot(n,dx,incx,dy,incy)
!
!     LATEST REVISION  -  OCTOBER 3, 1983  (JRD)
!
!     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
!     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      double precision dx(*),dy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         ddot = ddot + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +&
     &    dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
   50 continue
      return
!
!         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          ddot = ddot + dx(i)*dy(i)
   70     continue
      return
      end
!DAXPY
      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY.
!     FOR I = 0 TO N-1, REPLACE  DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
!       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision da
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      double precision dx(*),dy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0.or.da.eq.0.d0) return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          dy(i) = da*dx(i) + dy(i)
   70     continue
      return
      end
!SDOT
      real function sdot(n,sx,incx,sy,incy)
!
!     LATEST REVISION  -  OCTOBER 3, 1983  (JRD)
!
!     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
!     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      real sx(*),sy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1)5,20,60
    5 continue
!
!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sdot = sdot + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sdot = sdot + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
       sdot = sdot + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +&
     &    sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 continue
      return
!
!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 continue
      ns=n*incx
      do 70 i=1,ns,incx
        sdot = sdot + sx(i)*sy(i)
   70   continue
      return
      end
!IDAMAX
      integer function idamax(n,dx,incx)
!
!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF DOUBLE PRECISION DX.
!     IDAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(DX(1-INCX+I*INCX))
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,n
!
!  ARRAY ARGUMENTS
      double precision dx(*)
!
!  LOCAL SCALARS
      double precision dmax,xmag
      integer i,ii,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic dabs
!
      idamax = 0
      if(n.le.0) return
      idamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      dmax = dabs(dx(1))
      ns = n*incx
      ii = 1
          do 10 i = 1,ns,incx
          xmag = dabs(dx(i))
          if(xmag.le.dmax) go to 5
          idamax = ii
          dmax = xmag
    5     ii = ii + 1
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
          xmag = dabs(dx(i))
          if(xmag.le.dmax) go to 30
          idamax = i
          dmax = xmag
   30 continue
      return
      end
!SSWAP
      subroutine sswap (n,sx,incx,sy,incy)
!
!     INTERCHANGE SINGLE PRECISION SX AND SINGLE PRECISION SY.
!     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      real sx(*),sy(*)
!
!  LOCAL SCALARS
      real stemp1,stemp2,stemp3
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp1 = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp1
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp1 = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp1
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp1 = sx(i)
        stemp2 = sx(i+1)
        stemp3 = sx(i+2)
        sx(i) = sy(i)
        sx(i+1) = sy(i+1)
        sx(i+2) = sy(i+2)
        sy(i) = stemp1
        sy(i+1) = stemp2
        sy(i+2) = stemp3
   50 continue
      return
   60 continue
!
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
      ns = n*incx
        do 70 i=1,ns,incx
        stemp1 = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp1
   70   continue
      return
      end
!DSIDI
      subroutine dsidi(a,lda,n,kpvt,det,inert,work,job)
!
!     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
!     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM
!     DSIFA.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer job,lda,n
!
!  ARRAY ARGUMENTS
      double precision a(lda,*),det(2),work(*)
      integer inert(3),kpvt(*)
!
!  LOCAL SCALARS
      double precision ak,akkp1,akp1,d,t,temp,ten
      integer j,jb,k,km1,ks,kstep
      logical nodet,noert,noinv
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION DDOT
!       EXTERNAL DDOT
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DAXPY,DCOPY,DSWAP
!
!  INTRINSIC FUNCTIONS
      intrinsic dabs,iabs,mod
!
!
!     ON ENTRY
!
!        A       DOUBLE PRECISION(LDA,N)
!                THE OUTPUT FROM DSIFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY A.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX A.
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DSIFA.
!
!        WORK    DOUBLE PRECISION(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
!
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
!
!     ON RETURN
!
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
!
!        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
!               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE
!               IS NEVER REFERENCED.
!
!        DET    DOUBLE PRECISION(2)
!               DETERMINANT OF ORIGINAL MATRIX.
!               DETERMINANT = DET(1) * 10.0**DET(2)
!               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0
!               OR DET(1) = 0.0.
!
!        INERT  INTEGER(3)
!               THE INERTIA OF THE ORIGINAL MATRIX.
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED
!        AND  DSICO  HAS SET RCOND .EQ. 0.0
!        OR  DSIFA  HAS SET  INFO .NE. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DCOPY,DDOT,DSWAP
!     FORTRAN DABS,IABS,MOD
!
!
      ten = 10.0d0
!
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
!
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0d0
            det(2) = 0.0d0
   20    continue
         t = 0.0d0
         do 130 k = 1, n
            d = a(k,k)
!
!           CHECK IF 1 BY 1
!
            if (kpvt(k) .gt. 0) go to 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               if (t .ne. 0.0d0) go to 30
                  t = dabs(a(k,k+1))
                  d = (d/t)*a(k+1,k+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0d0
   40          continue
   50       continue
!
            if (noert) go to 60
               if (d .gt. 0.0d0) inert(1) = inert(1) + 1
               if (d .lt. 0.0d0) inert(2) = inert(2) + 1
               if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60       continue
!
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0d0) go to 110
   70             if (dabs(det(1)) .ge. 1.0d0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0d0
                  go to 70
   80             continue
   90             if (dabs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0d0
                  go to 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
!
!     COMPUTE INVERSE(A)
!
      if (noinv) go to 270
         k = 1
  150    if (k .gt. n) go to 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 180
!
!              1 BY 1
!
               a(k,k) = 1.0d0/a(k,k)
               if (km1 .lt. 1) go to 170
                  call dcopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
!
!              2 BY 2
!
               t = dabs(a(k,k+1))
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0d0)
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call dcopy(km1,a(1,k+1),1,work,1)
                  do 190 j = 1, km1
                     a(j,k+1) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1) + ddot(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + ddot(km1,a(1,k),1,a(1,k+1),1)
                  call dcopy(km1,a(1,k),1,work,1)
                  do 200 j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  210          continue
               kstep = 2
  220       continue
!
!           SWAP
!
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               call dswap(ks,a(1,ks),1,a(1,k),1)
               do 230 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) go to 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
!DTRCO
      subroutine dtrco(t,ldt,n,rcond,z,job)
!***BEGIN PROLOGUE  DTRCO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A3
!***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
!             MATRIX,TRIANGULAR
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR
!            MATRIX.
!***DESCRIPTION
!     DTRCO ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR
!     MATRIX.
!     ON ENTRY
!        T       DOUBLE PRECISION(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX.  THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!        JOB     INTEGER
!                = 0         T  IS LOWER TRIANGULAR.
!                = NONZERO   T  IS UPPER TRIANGULAR.
!     ON RETURN
!        RCOND   DOUBLE PRECISION
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .
!                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS
!                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!        Z       DOUBLE PRECISION(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DASUM,DAXPY,DSCAL
!***END PROLOGUE  DTRCO

!...SCALAR ARGUMENTS
      double precision rcond
      integer job,ldt,n

!...ARRAY ARGUMENTS
      double precision t(ldt,*),z(*)

!...LOCAL SCALARS
      double precision ek,s,sm,tnorm,w,wk,wkm,ynorm
      integer i1,j,j1,j2,k,kk,l
      logical lower

!...EXTERNAL FUNCTIONS
!      DOUBLE PRECISION DASUM
!       EXTERNAL DASUM

!...EXTERNAL SUBROUTINES
!       EXTERNAL DAXPY,DSCAL

!...INTRINSIC FUNCTIONS
      intrinsic dabs,dmax1,dsign

!***FIRST EXECUTABLE STATEMENT  DTRCO

      lower = job .eq. 0

!     COMPUTE 1-NORM OF T

      tnorm = 0.0d0
      do 10 j = 1, n
         l = j
         if (lower) l = n + 1 - j
         i1 = 1
         if (lower) i1 = j
         tnorm = dmax1(tnorm,dasum(l,t(i1,j),1))
   10 continue

!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.

!     SOLVE TRANS(T)*Y = E

      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 kk = 1, n
         k = kk
         if (lower) k = n + 1 - kk
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(t(k,k))) go to 30
            s = dabs(t(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (t(k,k) .eq. 0.0d0) go to 40
            wk = wk/t(k,k)
            wkm = wkm/t(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         if (kk .eq. n) go to 90
            j1 = k + 1
            if (lower) j1 = 1
            j2 = n
            if (lower) j2 = k - 1
            do 60 j = j1, j2
               sm = sm + dabs(z(j)+wkm*t(k,j))
               z(j) = z(j) + wk*t(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               w = wkm - wk
               wk = wkm
               do 70 j = j1, j2
                  z(j) = z(j) + w*t(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)

      ynorm = 1.0d0

!     SOLVE T*Z = Y

      do 130 kk = 1, n
         k = n + 1 - kk
         if (lower) k = kk
         if (dabs(z(k)) .le. dabs(t(k,k))) go to 110
            s = dabs(t(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  110    continue
         if (t(k,k) .ne. 0.0d0) z(k) = z(k)/t(k,k)
         if (t(k,k) .eq. 0.0d0) z(k) = 1.0d0
         i1 = 1
         if (lower) i1 = k + 1
         if (kk .ge. n) go to 120
            w = -z(k)
            call daxpy(n-kk,w,t(i1,k),1,z(i1),1)
  120    continue
  130 continue
!     MAKE ZNORM = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm

      if (tnorm .ne. 0.0d0) rcond = ynorm/tnorm
      if (tnorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
!SSIDI
      subroutine ssidi(a,lda,n,kpvt,det,inert,work,job)
!
!     LATEST REVISION  -  JANUARY 24, 1990  (JRD)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer job,lda,n
!
!  ARRAY ARGUMENTS
      real a(lda,*),det(2),work(*)
      integer inert(3),kpvt(*)
!
!  LOCAL SCALARS
      real ak,akkp1,akp1,d,t,temp,ten
      integer j,jb,k,km1,ks,kstep
      logical nodet,noert,noinv
!
!  EXTERNAL FUNCTIONS
!      REAL SDOT
!       EXTERNAL SDOT
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL SAXPY,SSWAP
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,iabs,mod
!
!     SSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
!     OF A REAL SYMMETRIC MATRIX USING THE FACTORS FROM SSIFA.
!
!     ON ENTRY
!
!        A       REAL(LDA,N)
!                THE OUTPUT FROM SSIFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY A.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX A.
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SSIFA.
!
!        WORK    REAL(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
!
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
!
!     ON RETURN
!
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
!
!        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
!               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE
!               IS NEVER REFERENCED.
!
!        DET    REAL(2)
!               DETERMINANT OF ORIGINAL MATRIX.
!               DETERMINANT = DET(1) * 10.0**DET(2)
!               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0
!               OR DET(1) = 0.0.
!
!        INERT  INTEGER(3)
!               THE INERTIA OF THE ORIGINAL MATRIX.
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED
!        AND  SSICO  HAS SET RCOND .EQ. 0.0
!        OR  SSIFA  HAS SET  INFO .NE. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SDOT,SSWAP
!     FORTRAN ABS,IABS,MOD
!
!
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
!
      ten = 10.0e0
!
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0e0
            det(2) = 0.0e0
   20    continue
         t = 0.0e0
         do 130 k = 1, n
            d = a(k,k)
!
!           CHECK IF 1 BY 1
!
            if (kpvt(k) .gt. 0) go to 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               if (t .ne. 0.0e0) go to 30
                  t = abs(a(k,k+1))
                  d = (d/t)*a(k+1,k+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0e0
   40          continue
   50       continue
!
            if (noert) go to 60
               if (d .gt. 0.0e0) inert(1) = inert(1) + 1
               if (d .lt. 0.0e0) inert(2) = inert(2) + 1
               if (d .eq. 0.0e0) inert(3) = inert(3) + 1
   60       continue
!
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0e0) go to 110
   70             if (abs(det(1)) .ge. 1.0e0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0e0
                  go to 70
   80             continue
   90             if (abs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0e0
                  go to 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
!
!     COMPUTE INVERSE(A)
!
      if (noinv) go to 270
         k = 1
  150    if (k .gt. n) go to 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 180
!
!              1 BY 1
!
               a(k,k) = 1.0e0/a(k,k)
               if (km1 .lt. 1) go to 170
                  call scopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = sdot(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + sdot(km1,work,1,a(1,k),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
!
!              2 BY 2
!
               t = abs(a(k,k+1))
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0e0)
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call scopy(km1,a(1,k+1),1,work,1)
                  do 190 j = 1, km1
                     a(j,k+1) = sdot(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1) + sdot(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + sdot(km1,a(1,k),1,a(1,k+1),1)
                  call scopy(km1,a(1,k),1,work,1)
                  do 200 j = 1, km1
                     a(j,k) = sdot(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k) + sdot(km1,work,1,a(1,k),1)
  210          continue
               kstep = 2
  220       continue
!
!           SWAP
!
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               call sswap(ks,a(1,ks),1,a(1,k),1)
               do 230 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) go to 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
!STRDI
      subroutine strdi(t,ldt,n,det,job,info)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer info,job,ldt,n
!
!  ARRAY ARGUMENTS
      real det(2),t(ldt,*)
!
!  LOCAL SCALARS
      real temp,ten
      integer i,j,k,kb,km1,kp1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL SAXPY,SSCAL
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,mod
!
!
!     STRDI COMPUTES THE DETERMINANT AND INVERSE OF A REAL
!     TRIANGULAR MATRIX.
!
!     ON ENTRY
!
!        T       REAL(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!
!        JOB     INTEGER
!                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.
!                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.
!                = 100       DET, NO INVERSE.
!                = 110       DET, INVERSE OF LOWER TRIANGULAR.
!                = 111       DET, INVERSE OF UPPER TRIANGULAR.
!
!     ON RETURN
!
!        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     REAL(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!        INFO    INTEGER
!                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR
!                AND THE INVERSE IS REQUESTED.
!                OTHERWISE INFO CONTAINS THE INDEX OF
!                A ZERO DIAGONAL ELEMENT OF T.
!
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SSCAL
!     FORTRAN ABS,MOD
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 180
!
!        COMPUTE DETERMINANT
!
         if (job/100 .eq. 0) go to 70
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
            do 50 i = 1, n
               det(1) = t(i,i)*det(1)
!           ...EXIT
               if (det(1) .eq. 0.0e0) go to 60
   10          if (abs(det(1)) .ge. 1.0e0) go to 20
                  det(1) = ten*det(1)
                  det(2) = det(2) - 1.0e0
               go to 10
   20          continue
   30          if (abs(det(1)) .lt. ten) go to 40
                  det(1) = det(1)/ten
                  det(2) = det(2) + 1.0e0
               go to 30
   40          continue
   50       continue
   60       continue
   70    continue
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
         if (mod(job/10,10) .eq. 0) go to 170
            if (mod(job,10) .eq. 0) go to 120
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  do 100 k = 1, n
                     info = k
!              ......EXIT
                     if (t(k,k) .eq. 0.0e0) go to 110
                     t(k,k) = 1.0e0/t(k,k)
                     temp = -t(k,k)
                     call sscal(k-1,temp,t(1,k),1)
                     kp1 = k + 1
                     if (n .lt. kp1) go to 90
                     do 80 j = kp1, n
                        temp = t(k,j)
                        t(k,j) = 0.0e0
                        call saxpy(k,temp,t(1,k),1,t(1,j),1)
   80                continue
   90                continue
  100             continue
                  info = 0
  110          continue
            go to 160
  120       continue
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
               do 150 kb = 1, n
                  k = n + 1 - kb
                  info = k
!     ............EXIT
                  if (t(k,k) .eq. 0.0e0) go to 180
                  t(k,k) = 1.0e0/t(k,k)
                  temp = -t(k,k)
                  if (k .ne. n) call sscal(n-k,temp,t(k+1,k),1)
                  km1 = k - 1
                  if (km1 .lt. 1) go to 140
                  do 130 j = 1, km1
                     temp = t(k,j)
                     t(k,j) = 0.0e0
                     call saxpy(n-k+1,temp,t(k,k),1,t(k,j),1)
  130             continue
  140             continue
  150          continue
               info = 0
  160       continue
  170    continue
  180 continue
      return
      end
!SCOPY
      subroutine scopy(n,sx,incx,sy,incy)
!
!     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY.
!     FOR I = 0 TO N-1, COPY  SX(LX+I*INCX) TO SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      real sx(*),sy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sx(i)
   70     continue
      return
      end
!SAXPY
      subroutine saxpy(n,sa,sx,incx,sy,incy)
!
!     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
!     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
!       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real sa
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      real sx(*),sy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0.or.sa.eq.0.e0) return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sa*sx(i) + sy(i)
   70     continue
      return
      end
!SASUM
      real function sasum(n,sx,incx)
!***BEGIN PROLOGUE  SASUM
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A3A
!***KEYWORDS  ADD,BLAS,LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  SUM OF MAGNITUDES OF S.P VECTOR COMPONENTS
!***DESCRIPTION
!                B L A S  SUBPROGRAM
!    DESCRIPTION OF PARAMETERS
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!     --OUTPUT--
!    SASUM  SINGLE PRECISION RESULT (ZERO IF N .LE. 0)
!     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
!     SASUM = SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  SASUM

!...SCALAR ARGUMENTS
     integer&
     &   incx,n

!...ARRAY ARGUMENTS
      real sx(*)

!...LOCAL SCALARS
     integer&
     &   i,m,mp1,ns

!...INTRINSIC FUNCTIONS
     intrinsic&
     &   abs,mod

!***FIRST EXECUTABLE STATEMENT  SASUM

      sasum = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)goto 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      ns = n*incx
          do 10 i=1,ns,incx
          sasum = sasum + abs(sx(i))
   10     continue
      return

!        CODE FOR INCREMENTS EQUAL TO 1.

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sasum = sasum + abs(sx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
       sasum = sasum + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))&
     &  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
      return
      end
!DCOPY
      subroutine dcopy(n,dx,incx,dy,incy)
!
!     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.
!     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,incy,n
!
!  ARRAY ARGUMENTS
      double precision dx(*),dy(*)
!
!  LOCAL SCALARS
      integer i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns=n*incx
          do 70 i=1,ns,incx
          dy(i) = dx(i)
   70     continue
      return
      end
!SNRM2
      real function snrm2 ( n, sx, incx)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer incx,n
!
!  ARRAY ARGUMENTS
      real sx(*)
!
!  LOCAL SCALARS
      real cuthi,cutlo,hitest,one,sum,xmax,zero
      integer i,j,next,nn
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,float,sqrt
!
      data   zero, one /0.0e0, 1.0e0/
!
!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N .LE. 0 RETURN WITH RESULT = 0.
!     IF N .GE. 1 THEN INCX MUST BE .GE. 1
!
!           C.L.LAWSON, 1978 JAN 08
!
!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
!
      xmax = zero
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
!
   10 assign 30 to next
      sum = zero
      nn = n * incx
!                                                 BEGIN MAIN LOOP
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        PHASE 1.  SUM IS ZERO
!
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
!
!                                PREPARE FOR PHASE 2.
      assign 70 to next
      go to 105
!
!                                PREPARE FOR PHASE 4.
!
  100 i = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
!
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 sum = (sum * xmax) * xmax
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 hitest = cuthi/float( n )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      do 95 j =i,nn,incx
      if(abs(sx(j)) .ge. hitest) go to 100
   95    sum = sum + sx(j)**2
      snrm2 = sqrt( sum )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
      end

!R9LGIC
      real function r9lgic (a, x, alx)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG COMPLEMENTARY INCOMPLETE GAMMA FUNCTION FOR LARGE X
! AND FOR A .LE. X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,alx,x
!
!  LOCAL SCALARS
      real eps,fk,p,r,s,t,xma,xpa
      integer k
!
!  EXTERNAL FUNCTIONS
!      REAL R1MACH
!       EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log
!
      data eps / 0.0 /
!
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
!
      xpa = x + 1.0 - a
      xma = x - 1.0 - a
!
      r = 0.0
      p = 1.0
      s = p
      do 10 k=1,200
        fk = k
        t = fk*(a-fk)*(1.0+r)
        r = -t/((xma+2.0*fk)*(xpa+2.0*fk)+t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 20
 10   continue
     call xerror (  'R9LGIC  NO CONVERGENCE IN 200 TERMS OF CONTINUED F&
     &RACTION', 57, 1, 2)
!
 20   r9lgic = a*alx - x + log(s/xpa)
!
      return
      end
!GAMMA
      real function gamma (x)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real dxrel,pi,sinpiy,sq2pil,xmax,xmin,y
      integer i,n,ngcs
!
!  LOCAL ARRAYS
      real gcs(23)
!
!  EXTERNAL FUNCTIONS
!      REAL CSEVL,R1MACH,R9LGMC
!      INTEGER INITS
!       EXTERNAL CSEVL,R1MACH,R9LGMC,INITS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL GAMLIM,XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,aint,exp,float,log,sin,sqrt
!
!
      data gcs   ( 1) / .008571195590989331e0/
      data gcs   ( 2) / .004415381324841007e0/
      data gcs   ( 3) / .05685043681599363e0/
      data gcs   ( 4) /-.004219835396418561e0/
      data gcs   ( 5) / .001326808181212460e0/
      data gcs   ( 6) /-.0001893024529798880e0/
      data gcs   ( 7) / .0000360692532744124e0/
      data gcs   ( 8) /-.0000060567619044608e0/
      data gcs   ( 9) / .0000010558295463022e0/
      data gcs   (10) /-.0000001811967365542e0/
      data gcs   (11) / .0000000311772496471e0/
      data gcs   (12) /-.0000000053542196390e0/
      data gcs   (13) / .0000000009193275519e0/
      data gcs   (14) /-.0000000001577941280e0/
      data gcs   (15) / .0000000000270798062e0/
      data gcs   (16) /-.0000000000046468186e0/
      data gcs   (17) / .0000000000007973350e0/
      data gcs   (18) /-.0000000000001368078e0/
      data gcs   (19) / .0000000000000234731e0/
      data gcs   (20) /-.0000000000000040274e0/
      data gcs   (21) / .0000000000000006910e0/
      data gcs   (22) /-.0000000000000001185e0/
      data gcs   (23) / .0000000000000000203e0/
!
      data pi /3.14159265358979324e0/
! SQ2PIL IS LOG (SQRT (2.*PI) )
      data sq2pil /0.91893853320467274e0/
      data ngcs, xmin, xmax, dxrel /0, 3*0.0 /
!
      if (ngcs.ne.0) go to 10
!
! ---------------------------------------------------------------------
! INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF
! TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
! THAN MACHINE PRECISION.
!
      ngcs = inits (gcs, 23, 0.1*r1mach(3))
!
      call gamlim (xmin, xmax)
      dxrel = sqrt (r1mach(4))
!
! ---------------------------------------------------------------------
! FINISH INITIALIZATION.  START EVALUATING GAMMA(X).
!
 10   y = abs(x)
      if (y.gt.10.0) go to 50
!
! COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND
! FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL.
!
      n = x
      if (x.lt.0.) n = n - 1
      y = x - float(n)
      n = n - 1
      gamma = 0.9375 + csevl(2.*y-1., gcs, ngcs)
      if (n.eq.0) return
!
      if (n.gt.0) go to 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.
!
      n = -n
      if (x.eq.0.) call xerror ('GAMMA   X IS 0', 14, 4, 2)
      if (x.lt.0. .and. x+float(n-2).eq.0.) call xerror ( 'GAMMA   X IS A NEGATIVE INTEGER', 31, 4, 2)
      if (x.lt.(-0.5) .and. abs((x-aint(x-0.5))/x).lt.dxrel) &
      & call xerror (  'GAMMA   ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 68, 1, 1)
!
      do i=1,n
        gamma = gamma / (x+float(i-1))
      enddo
      return
!
! GAMMA(X) FOR X .GE. 2.
!
 30   do 40 i=1,n
        gamma = (y+float(i))*gamma
 40   continue
      return
!
! COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   if (x.gt.xmax) call xerror ('GAMMA   X SO BIG GAMMA OVERFLOWS', 32, 3, 2)
!
      gamma = 0.
      if (x.lt.xmin) call xerror ('GAMMA   X SO SMALL GAMMA UNDERFLOWS', 35, 2, 1)
      if (x.lt.xmin) return
!
      gamma = exp((y-0.5)*log(y) - y + sq2pil + r9lgmc(y) )
      if (x.gt.0.) return
!
      if (abs((x-aint(x-0.5))/x).lt.dxrel) call xerror ( 'GAMMA   ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 61, 1, 1)
!
      sinpiy = sin (pi*y)
      if (sinpiy.eq.0.) call xerror ( 'GAMMA   X IS A NEGATIVE INTEGER', 31, 4, 2)
!
      gamma = -pi / (y*sinpiy*gamma)
!
      return
      end
!DLNGAM
      double precision function dlngam (x)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision dxrel,pi,sinpiy,sq2pil,sqpi2l,xmax,y
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,D9LGMC,DGAMMA
!       EXTERNAL D1MACH,D9LGMC,DGAMMA
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsqrt,int,log,sin
!
!
      data sq2pil / 0.91893853320467274178032973640562d0 /
! SQ2PIL = LOG (SQRT(2*PI)),  SQPI2L = LOG(SQRT(PI/2))
      data sqpi2l / +.225791352644727432363097614947441d+0    /
      data pi / 3.14159265358979323846264338327950d0 /
!
      data xmax, dxrel / 2*0.d0 /
!
      if (xmax.ne.0.d0) go to 10
      xmax = d1mach(2)/log(d1mach(2))
      dxrel = dsqrt (d1mach(4))
!
 10   y = abs (x)
      if (y.gt.10.d0) go to 20
!
! LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
!
      dlngam = log (abs (dgamma(x)) )
      return
!
! LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
!
20   if (y.gt.xmax) call xerror (&
     &  'DLNGAM  ABS(X) SO BIG DLNGAM OVERFLOWS', 39, 2, 2)
!
      if (x.gt.0.d0) then
         dlngam = sq2pil + (x-0.5d0)*log(x) - x + d9lgmc(y)
         return
      end if
!
      sinpiy = abs (sin(pi*y))
     if (sinpiy.eq.0.d0) call xerror (&
     &  'DLNGAM  X IS A NEGATIVE INTEGER', 31, 3, 2)
!
     if (abs ((x-int(x-0.5d0))/x).lt.dxrel) call xerror (&
    &    'DLNGAM  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE&
     &INTEGER', 68, 1, 1)
!
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
      return
!
      end
!ERFC
      real function erfc (x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real eta,sqeps,sqrtpi,xmax,xsml,y
      integer nterc2,nterf,nterfc
!
!  LOCAL ARRAYS
      real erc2cs(23),erfccs(24),erfcs(13)
!
!  EXTERNAL FUNCTIONS
!      REAL CSEVL,R1MACH
!      INTEGER INITS
!       EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,exp,log,sqrt
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000D+00
!                                        WITH WEIGHTED ERROR   7.10E-18
!                                         LOG WEIGHTED ERROR  17.15
!                               SIGNIFICANT FIGURES REQUIRED  16.31
!                                    DECIMAL PLACES REQUIRED  17.71
!
      data erfcs( 1) /   -.049046121234691808e0 /
      data erfcs( 2) /   -.14226120510371364e0 /
      data erfcs( 3) /    .010035582187599796e0 /
      data erfcs( 4) /   -.000576876469976748e0 /
      data erfcs( 5) /    .000027419931252196e0 /
      data erfcs( 6) /   -.000001104317550734e0 /
      data erfcs( 7) /    .000000038488755420e0 /
      data erfcs( 8) /   -.000000001180858253e0 /
      data erfcs( 9) /    .000000000032334215e0 /
      data erfcs(10) /   -.000000000000799101e0 /
      data erfcs(11) /    .000000000000017990e0 /
      data erfcs(12) /   -.000000000000000371e0 /
      data erfcs(13) /    .000000000000000007e0 /
!
! SERIES FOR ERFC       ON THE INTERVAL  0.          TO  2.50000D-01
!                                        WITH WEIGHTED ERROR   4.81E-17
!                                         LOG WEIGHTED ERROR  16.32
!                        APPROX SIGNIFICANT FIGURES REQUIRED  15.0
! SERIES FOR ERC2       ON THE INTERVAL  2.50000D-01 TO  1.00000D+00
!                                        WITH WEIGHTED ERROR   5.22E-17
!                                         LOG WEIGHTED ERROR  16.28
!                        APPROX SIGNIFICANT FIGURES REQUIRED  15.0
!                                    DECIMAL PLACES REQUIRED  16.96
!
      data erc2cs( 1) /   -.069601346602309501e0 /
      data erc2cs( 2) /   -.041101339362620893e0 /
      data erc2cs( 3) /    .003914495866689626e0 /
      data erc2cs( 4) /   -.000490639565054897e0 /
      data erc2cs( 5) /    .000071574790013770e0 /
      data erc2cs( 6) /   -.000011530716341312e0 /
      data erc2cs( 7) /    .000001994670590201e0 /
      data erc2cs( 8) /   -.000000364266647159e0 /
      data erc2cs( 9) /    .000000069443726100e0 /
      data erc2cs(10) /   -.000000013712209021e0 /
      data erc2cs(11) /    .000000002788389661e0 /
      data erc2cs(12) /   -.000000000581416472e0 /
      data erc2cs(13) /    .000000000123892049e0 /
      data erc2cs(14) /   -.000000000026906391e0 /
      data erc2cs(15) /    .000000000005942614e0 /
      data erc2cs(16) /   -.000000000001332386e0 /
      data erc2cs(17) /    .000000000000302804e0 /
      data erc2cs(18) /   -.000000000000069666e0 /
      data erc2cs(19) /    .000000000000016208e0 /
      data erc2cs(20) /   -.000000000000003809e0 /
      data erc2cs(21) /    .000000000000000904e0 /
      data erc2cs(22) /   -.000000000000000216e0 /
      data erc2cs(23) /    .000000000000000052e0 /
!
!                                    DECIMAL PLACES REQUIRED  17.01
!
      data erfccs( 1) /   0.0715179310202925e0 /
      data erfccs( 2) /   -.026532434337606719e0 /
      data erfccs( 3) /    .001711153977920853e0 /
      data erfccs( 4) /   -.000163751663458512e0 /
      data erfccs( 5) /    .000019871293500549e0 /
      data erfccs( 6) /   -.000002843712412769e0 /
      data erfccs( 7) /    .000000460616130901e0 /
      data erfccs( 8) /   -.000000082277530261e0 /
      data erfccs( 9) /    .000000015921418724e0 /
      data erfccs(10) /   -.000000003295071356e0 /
      data erfccs(11) /    .000000000722343973e0 /
      data erfccs(12) /   -.000000000166485584e0 /
      data erfccs(13) /    .000000000040103931e0 /
      data erfccs(14) /   -.000000000010048164e0 /
      data erfccs(15) /    .000000000002608272e0 /
      data erfccs(16) /   -.000000000000699105e0 /
      data erfccs(17) /    .000000000000192946e0 /
      data erfccs(18) /   -.000000000000054704e0 /
      data erfccs(19) /    .000000000000015901e0 /
      data erfccs(20) /   -.000000000000004729e0 /
      data erfccs(21) /    .000000000000001432e0 /
      data erfccs(22) /   -.000000000000000439e0 /
      data erfccs(23) /    .000000000000000138e0 /
      data erfccs(24) /   -.000000000000000048e0 /
!
      data sqrtpi /1.7724538509055160e0/
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0./
!
      if (nterf.ne.0) go to 10
      eta = 0.1*r1mach(3)
      nterf = inits (erfcs, 13, eta)
      nterfc = inits (erfccs, 24, eta)
      nterc2 = inits (erc2cs, 23, eta)
!
      xsml = -sqrt (-log(sqrtpi*r1mach(3)))
      xmax = sqrt (-log(sqrtpi*r1mach(1)))
      xmax = xmax - 0.5*log(xmax)/xmax - 0.01
      sqeps = sqrt (2.0*r1mach(3))
!
 10   if (x.gt.xsml) go to 20
!
! ERFC(X) = 1.0 - ERF(X) FOR X .LT. XSML
!
      erfc = 2.0
      return
!
 20   if (x.gt.xmax) go to 40
      y = abs(x)
      if (y.gt.1.0) go to 30
!
! ERFC(X) = 1.0 - ERF(X) FOR -1. .LE. X .LE. 1.
!
      if (y.lt.sqeps) then
         erfc = 1.0 - 2.0*x/sqrtpi
      else
         erfc = 1.0 - x*(1.0 + csevl (2.*x*x-1., erfcs, nterf) )
      end if
!
      return
!
! ERFC(X) = 1.0 - ERF(X) FOR 1. .LT. ABS(X) .LE. XMAX
!
 30   y = y*y
      if (y.le.4.) then
        erfc = exp(-y)/abs(x) *&
     &         (0.5 + csevl ((8.0/y-5.0)/3.0, erc2cs, nterc2) )
      else
        erfc = exp(-y)/abs(x) *&
     &          (0.5 + csevl (8.0/y-1.0, erfccs, nterfc) )
      end if
      if (x.lt.0.0) erfc = 2.0 - erfc
      return
!
 40   call xerror ('ERFC    X SO BIG ERFC UNDERFLOWS', 32, 1, 1)
      erfc = 0.0
      return
!
      end
!EPRINT
      subroutine eprint
!
!  THIS SUBROUTINE PRINTS THE LAST ERROR MESSAGE, IF ANY.
!
!
!  VARIABLE DECLARATIONS
!
!  LOCAL ARRAYS
      character messg(1)*4
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL E9RINT
!
!
      call e9rint(messg,1,1,.false.)
      return
!
      end
!DGAMR
      double precision function dgamr (x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! THIS ROUTINE, NOT DGAMMA(X), SHOULD BE THE FUNDAMENTAL ONE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision alngx,sgngx
      integer irold
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION DGAMMA
!       EXTERNAL DGAMMA
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DLGAMS,XERCLR,XGETF,XSETF
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,exp,int
!
!
      dgamr = 0.0d0
      if (x.le.0.0d0 .and. int(x).eq.x) return
!
      call xgetf (irold)
      call xsetf (1)
      if (abs(x).gt.10.0d0) go to 10
      dgamr = 1.0d0/dgamma(x)
      call xerclr
      call xsetf (irold)
      return
!
 10   call dlgams (x, alngx, sgngx)
      call xerclr
      call xsetf (irold)
      dgamr = sgngx * exp(-alngx)
      return
!
      end
!DLBETA
      double precision function dlbeta (a, b)
! JULY 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,b
!
!  LOCAL SCALARS
      double precision corr,p,q,sq2pil
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D9LGMC,DGAMMA,DLNGAM,DLNREL
!       EXTERNAL D9LGMC,DGAMMA,DLNGAM,DLNREL
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic log,max,min
!
      data sq2pil / 0.91893853320467274178032973640562d0 /
!
      p = min (a, b)
      q = max (a, b)
!
     if (p.le.0.d0) call xerror (&
     &  'DLBETA  BOTH ARGUMENTS MUST BE GT ZERO', 38, 1, 2)
!
      if (p.ge.10.d0) go to 30
      if (q.ge.10.d0) go to 20
!
! P AND Q ARE SMALL.
!
      dlbeta = log (dgamma(p) * (dgamma(q)/dgamma(p+q)) )
      return
!
! P IS SMALL, BUT Q IS BIG.
!
 20   corr = d9lgmc(q) - d9lgmc(p+q)
     dlbeta = dlngam(p) + corr + p - p*log(p+q)&
     &  + (q-0.5d0)*dlnrel(-p/(p+q))
      return
!
! P AND Q ARE BIG.
!
 30   corr = d9lgmc(p) + d9lgmc(q) - d9lgmc(p+q)
     dlbeta = -0.5d0*log(q) + sq2pil + corr + (p-0.5d0)*log(p/(p+q))&
     &  + q*dlnrel(-p/(p+q))
      return
!
      end
!DCSEVL
      double precision function dcsevl (x, a, n)
!
!     LATEST REVISION  -  OCTOBER 3, 1983  (JRD)
!
! EVALUATE THE N-TERM CHEBYSHEV SERIES A AT X.  ADAPTED FROM
! R. BROUCKE, ALGORITHM 446, C.A.C.M., 16, 254 (1973).
!
!             INPUT ARGUMENTS --
! X      DBLE PREC VALUE AT WHICH THE SERIES IS TO BE EVALUATED.
! A      DBLE PREC ARRAY OF N TERMS OF A CHEBYSHEV SERIES.  IN EVAL-
!        UATING A, ONLY HALF THE FIRST COEF IS SUMMED.
! N      NUMBER OF TERMS IN ARRAY A.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
      integer n
!
!  ARRAY ARGUMENTS
      double precision a(n)
!
!  LOCAL SCALARS
      double precision b0,b1,b2,twox
      integer i,ni1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!
      if (n.lt.1) call xerror ('DCSEVL  NUMBER OF TERMS LE 0', 28, 2,2)
     if (n.gt.1000) call xerror ('DCSEVL  NUMBER OF TERMS GT 1000',&
     &  31, 3, 2)
     if (x.lt.(-1.d0) .or. x.gt.1.d0) call xerror (&
     &  'DCSEVL  X OUTSIDE (-1,+1)', 25, 1, 1)
!
      twox = 2.0d0*x
      b0 = 0.0d0
      b1 = 0.0d0
      b2 = 0.0d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni1 = n-i+1
        b0 = twox*b1 - b2 + a(ni1)
 10   continue
!
      dcsevl = 0.5d0 * (b0-b2)
!
      return
      end
!GAMLIM
      subroutine gamlim (xmin, xmax)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! CALCULATE THE MINIMUM AND MAXIMUM LEGAL BOUNDS FOR X IN GAMMA(X).
! XMIN AND XMAX ARE NOT THE ONLY BOUNDS, BUT THEY ARE THE ONLY NON-
! TRIVIAL ONES TO CALCULATE.
!
!             OUTPUT ARGUMENTS --
! XMIN   MINIMUM LEGAL VALUE OF X IN GAMMA(X).  ANY SMALLER VALUE OF
!        X MIGHT RESULT IN UNDERFLOW.
! XMAX   MAXIMUM LEGAL VALUE OF X IN GAMMA(X).  ANY LARGER VALUE WILL
!        CAUSE OVERFLOW.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real xmax,xmin
!
!  LOCAL SCALARS
      real alnbig,alnsml,xln,xold
      integer i
!
!  EXTERNAL FUNCTIONS
!      REAL R1MACH
!       EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log,max
!
      alnsml = log(r1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
       xmin = xmin - xmin*((xmin+0.5)*xln - xmin - 0.2258 + alnsml)&
     &    / (xmin*xln + 0.5)
        if (abs(xmin-xold).lt.0.005) go to 20
 10   continue
      call xerror ('GAMLIM  UNABLE TO FIND XMIN', 27, 1, 2)
!
 20   xmin = -xmin + 0.01
!
      alnbig = log(r1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
       xmax = xmax - xmax*((xmax-0.5)*xln - xmax + 0.9189 - alnbig)&
     &    / (xmax*xln - 0.5)
        if (abs(xmax-xold).lt.0.005) go to 40
 30   continue
      call xerror ('GAMLIM  UNABLE TO FIND XMAX', 27, 2, 2)
!
 40   xmax = xmax - 0.01
      xmin = max (xmin, -xmax+1.)
!
      return
      end
!R9LGMC
      real function r9lgmc (x)
! AUGUST 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG GAMMA CORRECTION FACTOR FOR X .GE. 10.0 SO THAT
!  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real xbig,xmax
      integer nalgm
!
!  LOCAL ARRAYS
      real algmcs(6)
!
!  EXTERNAL FUNCTIONS
!      REAL CSEVL,R1MACH
!      INTEGER INITS
!       EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic exp,log,min,sqrt
!
!
! SERIES FOR ALGM       ON THE INTERVAL  0.          TO  1.00000D-02
!                                        WITH WEIGHTED ERROR   3.40E-16
!                                         LOG WEIGHTED ERROR  15.47
!                               SIGNIFICANT FIGURES REQUIRED  14.39
!                                    DECIMAL PLACES REQUIRED  15.86
!
      data algmcs( 1) /    .166638948045186e0 /
      data algmcs( 2) /   -.0000138494817606e0 /
      data algmcs( 3) /    .0000000098108256e0 /
      data algmcs( 4) /   -.0000000000180912e0 /
      data algmcs( 5) /    .0000000000000622e0 /
      data algmcs( 6) /   -.0000000000000003e0 /
!
      data nalgm, xbig, xmax / 0, 2*0.0 /
!
      if (nalgm.ne.0) go to 10
      nalgm = inits (algmcs, 6, r1mach(3))
      xbig = 1.0/sqrt(r1mach(3))
      xmax = exp (min(log(r1mach(2)/12.0), -log(12.0*r1mach(1))) )
!
 10   if (x.lt.10.0) call xerror ('R9LGMC  X MUST BE GE 10', 23, 1, 2)
      if (x.ge.xmax) go to 20
!
      r9lgmc = 1.0/(12.0*x)
      if (x.lt.xbig) r9lgmc = csevl (2.0*(10./x)**2-1., algmcs, nalgm)/x
      return
!
 20   r9lgmc = 0.0
      call xerror ('R9LGMC  X SO BIG R9LGMC UNDERFLOWS', 34, 2, 1)
      return
!
      end
!ALNGAM
      real function alngam (x)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real dxrel,pi,sinpiy,sq2pil,sqpi2l,xmax,y
!
!  EXTERNAL FUNCTIONS
!      REAL GAMMA,R1MACH,R9LGMC
!       EXTERNAL GAMMA,R1MACH,R9LGMC
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,aint,log,sin,sqrt
!
      data sq2pil / 0.91893853320467274e0/
! SQ2PIL = LOG(SQRT(2.*PI)),  SQPI2L = LOG (SQRT(PI/2.))
      data sqpi2l / 0.22579135264472743e0/
      data pi     / 3.14159265358979324e0/
!
      data xmax, dxrel / 0., 0. /
!
      if (xmax.ne.0.) go to 10
      xmax = r1mach(2)/log(r1mach(2))
      dxrel = sqrt (r1mach(4))
!
 10   y = abs(x)
      if (y.gt.10.0) go to 20
!
! LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
!
      alngam = log (abs (gamma(x)))
      return
!
! LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
!
20   if (y.gt.xmax) call xerror (&
     &  'ALNGAM  ABS(X) SO BIG ALNGAM OVERFLOWS', 38, 2, 2)
!
      if (x.gt.0.0) then
         alngam = sq2pil + (x-0.5)*log(x) - x + r9lgmc(y)
         return
      end if
!
      sinpiy = abs (sin(pi*y))
     if (sinpiy.eq.0.) call xerror ('ALNGAM  X IS A NEGATIVE INTEGER',&
     &  31, 3, 2)
!
     if (abs((x-aint(x-0.5))/x).lt.dxrel) call xerror (&
    &    'ALNGAM  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE&
     &INTEGER', 68, 1, 1)
!
      alngam = sqpi2l + (x-0.5)*log(y) - x - log(sinpiy) - r9lgmc(y)
      return
!
      end
!XERPRT
      subroutine xerprt(messg,nmessg)
!
!     ABSTRACT
!        PRINT THE HOLLERITH MESSAGE IN MESSG, OF LENGTH MESSG,
!        ON EACH FILE INDICATED BY XGETUA.
!        THIS VERSION PRINTS EXACTLY THE RIGHT NUMBER OF CHARACTERS,
!        NOT A NUMBER OF WORDS, AND THUS SHOULD WORK ON MACHINES
!        WHICH DO NOT BLANK FILL THE LAST WORD OF THE HOLLERITH.
!
!     RON JONES, JUNE 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer nmessg
!
!  ARRAY ARGUMENTS
      character messg(nmessg)*4
!
!  LOCAL SCALARS
     integer i,iunit,kunit,nchar,ncharl,nchlst,nchrem,nfield,nlines,&
     &   nunit,nword,nword1,nword2
      character la*1,lblank*1,lcom*1
!
!  LOCAL ARRAYS
      integer lun(5)
      character f(10)*1,g(14)*1
!
!  EXTERNAL FUNCTIONS
!      INTEGER I1MACH
!       EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL S88FMT,XGETUA
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
     data f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10)&
     &   / '(' ,'1' ,'X' ,',' ,' ' ,' ' ,'A' ,' ' ,' ' ,')' /
     data g(1),g(2),g(3),g(4),g(5),g(6),g(7),g(8),g(9),g(10)&
     &   / '(' ,'1' ,'X' ,' ' ,' ' ,' ' ,' ' ,' ' ,' ' ,' ' /
     data g(11),g(12),g(13),g(14)&
     &   / ' '  ,' '  ,' '  ,')'  /
      data la/'A'/,lcom/','/,lblank/' '/
!     PREPARE FORMAT FOR WHOLE LINES
      nchar = i1mach(6)
      nfield = 72/nchar
      call s88fmt(2,nfield,f(5))
      call s88fmt(2,nchar,f(8))
!     PREPARE FORMAT FOR LAST, PARTIAL LINE, IF NEEDED
      ncharl = nfield*nchar
      nlines = nmessg/ncharl
      nword  = nlines*nfield
      nchrem = nmessg - nlines*ncharl
      if (nchrem.le.0) go to 40
         do 10 i=4,13
10          g(i) = lblank
         nfield = nchrem/nchar
         if (nfield.le.0) go to 20
!        PREPARE WHOLE WORD FIELDS
            g(4) = lcom
            call s88fmt(2,nfield,g(5))
            g(7) = la
            call s88fmt(2,nchar,g(8))
20       continue
         nchlst = mod(nchrem,nchar)
         if (nchlst.le.0) go to 30
!        PREPARE PARTIAL WORD FIELD
            g(10) = lcom
            g(11) = la
            call s88fmt(2,nchlst,g(12))
30       continue
40    continue
!     PRINT THE MESSAGE
      nword1 = nword+1
      nword2 = (nmessg+nchar-1)/nchar
      call xgetua(lun,nunit)
      do 50 kunit = 1,nunit
         iunit = lun(kunit)
         if (iunit.eq.0) iunit = i1mach(4)
         if (nword.gt.0) write (iunit,f) (messg(i),i=1,nword)
         if (nchrem.gt.0) write (iunit,g) (messg(i),i=nword1,nword2)
50    continue
      return
      end
!INITS
      integer function inits (os, nos, eta)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! INITIALIZE THE ORTHOGONAL SERIES SO THAT INITS IS THE NUMBER OF TERMS
! NEEDED TO INSURE THE ERROR IS NO LARGER THAN ETA.  ORDINARILY, ETA
! WILL BE CHOSEN TO BE ONE-TENTH MACHINE PRECISION.
!
!             INPUT ARGUMENTS --
! OS     ARRAY OF NOS COEFFICIENTS IN AN ORTHOGONAL SERIES.
! NOS    NUMBER OF COEFFICIENTS IN OS.
! ETA    REQUESTED ACCURACY OF SERIES.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real eta
      integer nos
!
!  ARRAY ARGUMENTS
      real os(nos)
!
!  LOCAL SCALARS
      real err
      integer i,ii
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs
!
!
     if (nos.lt.1) call xerror (&
     &  'INITS   NUMBER OF COEFFICIENTS LT 1', 35, 2, 2)
!
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
!
20   if (i.eq.nos) call xerror ('INITS   ETA MAY BE TOO SMALL', 28,&
     &  1, 2)
      inits = i
!
      return
      end
!CSEVL
      real function csevl (x, cs, n)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE THE N-TERM CHEBYSHEV SERIES CS AT X.  ADAPTED FROM
! R. BROUCKE, ALGORITHM 446, C.A.C.M., 16, 254 (1973).  ALSO SEE FOX
! AND PARKER, CHEBYSHEV POLYS IN NUMERICAL ANALYSIS, OXFORD PRESS, P.56.
!
!             INPUT ARGUMENTS --
! X      VALUE AT WHICH THE SERIES IS TO BE EVALUATED.
! CS     ARRAY OF N TERMS OF A CHEBYSHEV SERIES.  IN EVAL-
!        UATING CS, ONLY HALF THE FIRST COEF IS SUMMED.
! N      NUMBER OF TERMS IN ARRAY CS.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
      integer n
!
!  ARRAY ARGUMENTS
      real cs(n)
!
!  LOCAL SCALARS
      real b0,b1,b2,twox
      integer i,ni
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!
      if (n.lt.1) call xerror ('CSEVL   NUMBER OF TERMS LE 0', 28, 2,2)
     if (n.gt.1000) call xerror ('CSEVL   NUMBER OF TERMS GT 1000',&
     &  31, 3, 2)
     if (x.lt.(-1.0) .or. x.gt.1.0) call xerror (&
     &  'CSEVL   X OUTSIDE (-1,+1)', 25, 1, 1)
!
      b0 = 0.0
      b1 = 0.0
      b2 = 0.0
      twox = 2.0*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
!
      csevl = 0.5 * (b0-b2)
!
      return
      end
!XGETUA
      subroutine xgetua(iunit,n)
!
!     ABSTRACT
!        XGETUA MAY BE CALLED TO DETERMINE THE UNIT NUMBER OR NUMBERS
!        TO WHICH ERROR MESSAGES ARE BEING SENT.
!        THESE UNIT NUMBERS MAY HAVE BEEN SET BY A CALL TO XSETUN,
!        OR A CALL TO XSETUA, OR MAY BE A DEFAULT VALUE.
!
!     DESCRIPTION OF PARAMETERS
!      --OUTPUT--
!        IUNIT - AN ARRAY OF ONE TO FIVE UNIT NUMBERS, DEPENDING
!                ON THE VALUE OF N.  A VALUE OF ZERO REFERS TO THE
!                DEFAULT UNIT, AS DEFINED BY THE I1MACH MACHINE
!                CONSTANT ROUTINE.  ONLY IUNIT(1),...,IUNIT(N) ARE
!                DEFINED BY XGETUA.  THE VALUES OF IUNIT(N+1),...,
!                IUNIT(5) ARE NOT DEFINED (FOR N.LT.5) OR ALTERED
!                IN ANY WAY BY XGETUA.
!        N     - THE NUMBER OF UNITS TO WHICH COPIES OF THE
!                ERROR MESSAGES ARE BEING SENT.  N WILL BE IN THE
!                RANGE FROM 1 TO 5.
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer n
!
!  ARRAY ARGUMENTS
      integer iunit(5)
!
!  LOCAL SCALARS
      integer i,index
!
!  EXTERNAL FUNCTIONS
!      INTEGER J4SAVE
!       EXTERNAL J4SAVE
!
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunit(i) = j4save(index,0,.false.)
   30 continue
      return
      end
!E9RINT
      subroutine e9rint(messg,nw,nerr,save)
!
!  THIS ROUTINE STORES THE CURRENT ERROR MESSAGE OR PRINTS THE OLD ONE,
!  IF ANY, DEPENDING ON WHETHER OR NOT SAVE = .TRUE. .
!
!     CHARACTER*4 MESSG(NW)
!     LOGICAL SAVE
!
!  MESSGP STORES AT LEAST THE FIRST 72 CHARACTERS OF THE PREVIOUS
!  MESSAGE. ITS LENGTH IS MACHINE DEPENDENT AND MUST BE AT LEAST
!
!       1 + 71/(THE NUMBER OF CHARACTERS STORED PER INTEGER WORD).
!
!     CHARACTER*4 MESSGP(36),FMT(14),CCPLUS
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer nerr,nw
      logical save
!
!  ARRAY ARGUMENTS
      character messg(nw)*4
!
!  LOCAL SCALARS
      integer i,iwunit,nerrp,nwp
      character ccplus*4
!
!  LOCAL ARRAYS
      character fmt(14)*4,messgp(36)*4
!
!  EXTERNAL FUNCTIONS
!      INTEGER I1MACH,I8SAVE
!       EXTERNAL I1MACH,I8SAVE
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL S88FMT
!
!
!  START WITH NO PREVIOUS MESSAGE.
!
      data messgp(1)/'1'/, nwp/0/, nerrp/0/
!
!  SET UP THE FORMAT FOR PRINTING THE ERROR MESSAGE.
!  THE FORMAT IS SIMPLY (A1,14X,72AXX) WHERE XX=I1MACH(6) IS THE
!  NUMBER OF CHARACTERS STORED PER INTEGER WORD.
!
      data ccplus  / '+' /
!
      data fmt( 1) / '(' /
      data fmt( 2) / 'A' /
      data fmt( 3) / '1' /
      data fmt( 4) / ',' /
      data fmt( 5) / '1' /
      data fmt( 6) / '4' /
      data fmt( 7) / 'X' /
      data fmt( 8) / ',' /
      data fmt( 9) / '7' /
      data fmt(10) / '2' /
      data fmt(11) / 'A' /
      data fmt(12) / 'X' /
      data fmt(13) / 'X' /
      data fmt(14) / ')' /
!
      if (.not.save) go to 20
!
!  SAVE THE MESSAGE.
!
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
!
        go to 30
!
 20   if (i8save(1,0,.false.).eq.0) go to 30
!
!  PRINT THE MESSAGE.
!
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(' ERROR ',i4,' IN ')
!
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
!
 30   return
!
      end
!XERROR
      subroutine xerror(messg,nmessg,nerr,level)
!
!     ABSTRACT
!        XERROR PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER
!        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE
!        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL.
!        (SEE SUBROUTINE XSETF FOR DETAILS.)
!
!     DESCRIPTION OF PARAMETERS
!      --INPUT--
!        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED, CONTAINING
!                NO MORE THAN 72 CHARACTERS.
!        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG.
!        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE.
!                NERR MUST NOT BE ZERO.
!        LEVEL - ERROR CATEGORY.
!                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR.
!                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS
!                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.)
!                =0 MEANS THIS IS A WARNING MESSAGE ONLY.
!                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE
!                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY
!                   TIMES THIS CALL IS EXECUTED.
!
!     EXAMPLES
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
!                    43,2,1)
!        CALL XERROR(  'ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL
!    1 FULLY COLLAPSED.',65,3,0)
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  7 FEB 1979
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer level,nerr,nmessg
!
!  ARRAY ARGUMENTS
      character messg(nmessg)*4
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERRWV
!
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)
      return
      end
!D9GMIT
      double precision function d9gmit(a,x,algap1,sgngam,alx)
!***BEGIN PROLOGUE  D9GMIT
!***DATE WRITTEN   770701   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  C7E
!***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
!             DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,SPECIAL FUNCTION,
!             TRICOMI
!***AUTHOR  FULLERTON, W., (LANL)
!***PURPOSE  COMPUTES D.P. TRICOMI-S INCOMPLETE GAMMA FUNCTION FOR
!            SMALL X.
!***DESCRIPTION
!
! COMPUTE TRICOMI'S INCOMPLETE GAMMA FUNCTION FOR SMALL X.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH,DLNGAM,XERROR
!***END PROLOGUE  D9GMIT
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,algap1,alx,sgngam,x
!
!  LOCAL SCALARS
      double precision ae,aeps,alg2,algs,bot,eps,fk,s,sgng2,t,te
      integer k,m,ma
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DLNGAM
!       EXTERNAL D1MACH,DLNGAM
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dble,dsign,exp,float,log
!
      data eps, bot / 2*0.d0 /
!***FIRST EXECUTABLE STATEMENT  D9GMIT
      if (eps.ne.0.d0) go to 10
      eps = 0.5d0*d1mach(3)
      bot = log (d1mach(1))
!
 10   if (x.le.0.d0) call xerror ( 'D9GMIT  X SHOULD BE GT 0', 24, 1, 2)
!
      ma = a + 0.5d0
      if (a.lt.0.d0) ma = a - 0.5d0
      aeps = a - dble(float(ma))
!
      ae = a
      if (a.lt.(-0.5d0)) ae = aeps
!
      t = 1.d0
      te = ae
      s = t
      do 20 k=1,200
        fk = k
        te = -x*te/fk
        t = te/(ae+fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 30
 20   continue
     call xerror&
     &('D9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES',54,2,2)
!
 30   if (a.ge.(-0.5d0)) then
         algs = -algap1 + log(s)
      else
         algs = -dlngam(1.d0+aeps) + log(s)
         s = 1.0d0
         m = -ma - 1
         if (m.eq.0) go to 50
         t = 1.0d0
         do 40 k=1,m
            t = x*t/(aeps-dble(float(m+1-k)))
            s = s + t
            if (abs(t).lt.eps*abs(s)) go to 50
 40      continue
!
 50      d9gmit = 0.0d0
         algs = -dble(float(ma))*log(x) + algs
         if (s.ne.0.0d0 .and. aeps.ne.0.0d0) then
            sgng2 = sgngam * dsign (1.0d0, s)
            alg2 = -x - algap1 + log(abs(s))
!
            if (alg2.gt.bot) d9gmit = sgng2 * exp(alg2)
            if (algs.gt.bot) d9gmit = d9gmit + exp(algs)
            return
         end if
      end if
!
      d9gmit = exp (algs)
      return
!
      end
!XSETF
      subroutine xsetf(kontrl)
!
!     ABSTRACT
!        XSETF SETS THE ERROR CONTROL FLAG VALUE TO KONTRL.
!        (KONTRL IS AN INPUT PARAMETER ONLY.)
!        THE FOLLOWING TABLE SHOWS HOW EACH MESSAGE IS TREATED,
!        DEPENDING ON THE VALUES OF KONTRL AND LEVEL.  (SEE XERROR
!        FOR DESCRIPTION OF LEVEL.)
!
!        IF KONTRL IS ZERO OR NEGATIVE, NO INFORMATION OTHER THAN THE
!        MESSAGE ITSELF (INCLUDING NUMERIC VALUES, IF ANY) WILL BE
!        PRINTED.  IF KONTRL IS POSITIVE, INTRODUCTORY MESSAGES,
!        TRACE-BACKS, ETC., WILL BE PRINTED IN ADDITION TO THE MESSAGE.
!
!              IABS(KONTRL)
!        LEVEL        0              1              2
!        VALUE
!          2        FATAL          FATAL          FATAL
!
!          1     NOT PRINTED      PRINTED         FATAL
!
!          0     NOT PRINTED      PRINTED        PRINTED
!
!         -1     NOT PRINTED      PRINTED        PRINTED
!                                  ONLY           ONLY
!                                  ONCE           ONCE
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  23 MAY 1979
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer kontrl
!
!  LOCAL SCALARS
      integer junk
!
!  EXTERNAL FUNCTIONS
!      INTEGER J4SAVE
!       EXTERNAL J4SAVE
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERRWV
!
      if ((kontrl.ge.(-2)).and.(kontrl.le.2)) go to 10
        call xerrwv('XSETF  -- INVALID VALUE OF KONTRL (I1).',33,1,2,&
     &   1,kontrl,0,0,0.,0.)
         return
   10 junk = j4save(2,kontrl,.true.)
      return
      end
!R9GMIT
      real function r9gmit (a, x, algap1, sgngam, alx)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE TRICOMI-S INCOMPLETE GAMMA FUNCTION FOR SMALL X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,algap1,alx,sgngam,x
!
!  LOCAL SCALARS
      real ae,aeps,alg2,algs,bot,eps,fk,s,sgng2,t,te
      integer k,m,ma
!
!  EXTERNAL FUNCTIONS
!      REAL ALNGAM,R1MACH
!       EXTERNAL ALNGAM,R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,exp,float,log,sign
!
      data eps, bot / 2*0.0 /
!
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
      if (bot.eq.0.0) bot = log(r1mach(1))
!
      if (x.le.0.0) call xerror ('R9GMIT  X SHOULD BE GT 0', 24, 1, 2)
!
      ma = a + 0.5
      if (a.lt.0.0) ma = a - 0.5
      aeps = a - float(ma)
!
      ae = a
      if (a.lt.(-0.5)) ae = aeps
!
      t = 1.0
      te = ae
      s = t
      do 20 k=1,200
        fk = k
        te = -x*te/fk
        t = te/(ae+fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 30
 20   continue
     call xerror (  'R9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SE&
     &RIES', 54, 2, 2)
!
 30   if (a.ge.(-0.5)) then
         algs = -algap1 + log(s)
      else
!
         algs = -alngam(1.0+aeps) + log(s)
         s = 1.0
         m = -ma - 1
         if (m.eq.0) go to 50
         t = 1.0
         do 40 k=1,m
            t = x*t/(aeps-float(m+1-k))
            s = s + t
            if (abs(t).lt.eps*abs(s)) go to 50
 40      continue
!
 50      r9gmit = 0.0
         algs = -float(ma)*log(x) + algs
         if (s.ne.0.0 .and. aeps.ne.0.0) then
            sgng2 = sgngam*sign(1.0,s)
            alg2 = -x - algap1 + log(abs(s))
!
            if (alg2.gt.bot) r9gmit = sgng2*exp(alg2)
            if (algs.gt.bot) r9gmit = r9gmit + exp(algs)
            return
         end if
      end if
      r9gmit = exp(algs)
      return
!
      end
!DLNREL
      double precision function dlnrel (x)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision xmin
      integer nlnrel
!
!  LOCAL ARRAYS
      double precision alnrcs(43)
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DCSEVL
!      INTEGER INITDS
!       EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsqrt,log,sngl
!
!
! SERIES FOR ALNR       ON THE INTERVAL -3.75000E-01 TO  3.75000E-01
!                                        WITH WEIGHTED ERROR   6.35E-32
!                                         LOG WEIGHTED ERROR  31.20
!                               SIGNIFICANT FIGURES REQUIRED  30.93
!                                    DECIMAL PLACES REQUIRED  32.01
!
      data alnrcs(  1) / +.10378693562743769800686267719098d+1     /
      data alnrcs(  2) / -.13364301504908918098766041553133d+0     /
      data alnrcs(  3) / +.19408249135520563357926199374750d-1     /
      data alnrcs(  4) / -.30107551127535777690376537776592d-2     /
      data alnrcs(  5) / +.48694614797154850090456366509137d-3     /
      data alnrcs(  6) / -.81054881893175356066809943008622d-4     /
      data alnrcs(  7) / +.13778847799559524782938251496059d-4     /
      data alnrcs(  8) / -.23802210894358970251369992914935d-5     /
      data alnrcs(  9) / +.41640416213865183476391859901989d-6     /
      data alnrcs( 10) / -.73595828378075994984266837031998d-7     /
      data alnrcs( 11) / +.13117611876241674949152294345011d-7     /
      data alnrcs( 12) / -.23546709317742425136696092330175d-8     /
      data alnrcs( 13) / +.42522773276034997775638052962567d-9     /
      data alnrcs( 14) / -.77190894134840796826108107493300d-10    /
      data alnrcs( 15) / +.14075746481359069909215356472191d-10    /
      data alnrcs( 16) / -.25769072058024680627537078627584d-11    /
      data alnrcs( 17) / +.47342406666294421849154395005938d-12    /
      data alnrcs( 18) / -.87249012674742641745301263292675d-13    /
      data alnrcs( 19) / +.16124614902740551465739833119115d-13    /
      data alnrcs( 20) / -.29875652015665773006710792416815d-14    /
      data alnrcs( 21) / +.55480701209082887983041321697279d-15    /
      data alnrcs( 22) / -.10324619158271569595141333961932d-15    /
      data alnrcs( 23) / +.19250239203049851177878503244868d-16    /
      data alnrcs( 24) / -.35955073465265150011189707844266d-17    /
      data alnrcs( 25) / +.67264542537876857892194574226773d-18    /
      data alnrcs( 26) / -.12602624168735219252082425637546d-18    /
      data alnrcs( 27) / +.23644884408606210044916158955519d-19    /
      data alnrcs( 28) / -.44419377050807936898878389179733d-20    /
      data alnrcs( 29) / +.83546594464034259016241293994666d-21    /
      data alnrcs( 30) / -.15731559416479562574899253521066d-21    /
      data alnrcs( 31) / +.29653128740247422686154369706666d-22    /
      data alnrcs( 32) / -.55949583481815947292156013226666d-23    /
      data alnrcs( 33) / +.10566354268835681048187284138666d-23    /
      data alnrcs( 34) / -.19972483680670204548314999466666d-24    /
      data alnrcs( 35) / +.37782977818839361421049855999999d-25    /
      data alnrcs( 36) / -.71531586889081740345038165333333d-26    /
      data alnrcs( 37) / +.13552488463674213646502024533333d-26    /
      data alnrcs( 38) / -.25694673048487567430079829333333d-27    /
      data alnrcs( 39) / +.48747756066216949076459519999999d-28    /
      data alnrcs( 40) / -.92542112530849715321132373333333d-29    /
      data alnrcs( 41) / +.17578597841760239233269760000000d-29    /
      data alnrcs( 42) / -.33410026677731010351377066666666d-30    /
      data alnrcs( 43) / +.63533936180236187354180266666666d-31    /
!
      data nlnrel, xmin / 0, 0.d0 /
!
      if (nlnrel.ne.0) go to 10
      nlnrel = initds (alnrcs, 43, 0.1*sngl(d1mach(3)))
      xmin = -1.0d0 + dsqrt(d1mach(4))
!
 10   if (x.le.(-1.d0)) call xerror ('DLNREL  X IS LE -1', 18, 2, 2)
     if (x.lt.xmin) call xerror (&
    &  'DLNREL  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 54,&
     &  1, 1)
!
      if (abs(x).le.0.375d0) then
         dlnrel = x*(1.0d0 - x*dcsevl (x/0.375d0, alnrcs, nlnrel))
      else
         dlnrel = log (1.0d0+x)
      end if
!
      return
      end
!GAMI
      real function gami (a, x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE THE INCOMPLETE GAMMA FUNCTION DEFINED BY
!
! GAMI = INTEGRAL FROM T = 0 TO X OF EXP(-T) * T**(A-1.0) .
!
! GAMI IS EVALUATED FOR POSITIVE VALUES OF A AND NON-NEGATIVE VALUES
! OF X.  A SLIGHT DETERIORATION OF 2 OR 3 DIGITS ACCURACY WILL OCCUR
! WHEN GAMI IS VERY LARGE OR VERY SMALL, BECAUSE LOGARITHMIC VARIABLES
! ARE USED.
!
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,x
!
!  LOCAL SCALARS
      real factor
!
!  EXTERNAL FUNCTIONS
!      REAL ALNGAM,GAMIT,R1MACH
!       EXTERNAL ALNGAM,GAMIT,R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic exp,log
!
      if (a.le.0.0) call xerror ('GAMI    A MUST BE GT ZERO', 25, 1, 2)
      if (x.lt.0.0) call xerror ('GAMI    X MUST BE GE ZERO', 25, 2, 2)
!
      gami = 0.0
      if (x.eq.0.0) return
!
! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
!
      factor = alngam(a) + a*log(x)
      if (factor.gt.log(r1mach(2))) then
         gami = r1mach(2)
      else
         gami = exp(factor) * gamit(a,x)
      end if
!
      return
      end
!D9LGIC
      double precision function d9lgic(a,x,alx)
!***BEGIN PROLOGUE  D9LGIC
!***DATE WRITTEN   770701   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  C7E
!***KEYWORDS  DOUBLE PRECISION,GAMMA,INCOMPLETE GAMMA FUNCTION,
!             LOGARITHM INCOMPLETE GAMMA FUNCTION,SPECIAL FUNCTION
!***AUTHOR  FULLERTON, W., (LANL)
!***PURPOSE  COMPUTES THE D.P. LOG INCOMPLETE GAMMA FUNCTION FOR LARGE X
!            AND FOR A .LE. X.
!***DESCRIPTION
!
! COMPUTE THE LOG COMPLEMENTARY INCOMPLETE GAMMA FUNCTION FOR LARGE X
! AND FOR A .LE. X.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH,XERROR
!***END PROLOGUE  D9LGIC
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,alx,x
!
!  LOCAL SCALARS
      double precision eps,fk,p,r,s,t,xma,xpa
      integer k
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH
!       EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log
!
      data eps / 0.d0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIC
      if (eps.eq.0.d0) eps = 0.5d0*d1mach(3)
!
      xpa = x + 1.0d0 - a
      xma = x - 1.d0 - a
!
      r = 0.d0
      p = 1.d0
      s = p
      do 10 k=1,300
        fk = k
        t = fk*(a-fk)*(1.d0+r)
        r = -t/((xma+2.d0*fk)*(xpa+2.d0*fk)+t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 20
 10   continue
     call xerror&
    &('D9LGIC  NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION',&
     & 57,1,2)
!
 20   d9lgic = a*alx - x + log(s/xpa)
!
      return
      end
!DERF
      double precision function derf (x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision sqeps,sqrtpi,xbig,y
      integer nterf
!
!  LOCAL ARRAYS
      double precision erfcs(21)
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DCSEVL,DERFC
!      INTEGER INITDS
!       EXTERNAL D1MACH,DCSEVL,DERFC,INITDS
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsign,dsqrt,log,sngl
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   1.28E-32
!                                         LOG WEIGHTED ERROR  31.89
!                               SIGNIFICANT FIGURES REQUIRED  31.05
!                                    DECIMAL PLACES REQUIRED  32.55
!
      data erfcs(  1) / -.49046121234691808039984544033376d-1     /
      data erfcs(  2) / -.14226120510371364237824741899631d+0     /
      data erfcs(  3) / +.10035582187599795575754676712933d-1     /
      data erfcs(  4) / -.57687646997674847650827025509167d-3     /
      data erfcs(  5) / +.27419931252196061034422160791471d-4     /
      data erfcs(  6) / -.11043175507344507604135381295905d-5     /
      data erfcs(  7) / +.38488755420345036949961311498174d-7     /
      data erfcs(  8) / -.11808582533875466969631751801581d-8     /
      data erfcs(  9) / +.32334215826050909646402930953354d-10    /
      data erfcs( 10) / -.79910159470045487581607374708595d-12    /
      data erfcs( 11) / +.17990725113961455611967245486634d-13    /
      data erfcs( 12) / -.37186354878186926382316828209493d-15    /
      data erfcs( 13) / +.71035990037142529711689908394666d-17    /
      data erfcs( 14) / -.12612455119155225832495424853333d-18    /
      data erfcs( 15) / +.20916406941769294369170500266666d-20    /
      data erfcs( 16) / -.32539731029314072982364160000000d-22    /
      data erfcs( 17) / +.47668672097976748332373333333333d-24    /
      data erfcs( 18) / -.65980120782851343155199999999999d-26    /
      data erfcs( 19) / +.86550114699637626197333333333333d-28    /
      data erfcs( 20) / -.10788925177498064213333333333333d-29    /
      data erfcs( 21) / +.12811883993017002666666666666666d-31    /
!
      data sqrtpi / 1.77245385090551602729816748334115d0 /
      data nterf, xbig, sqeps / 0, 2*0.d0 /
!
      if (nterf.ne.0) go to 10
      nterf = initds (erfcs, 21, 0.1*sngl(d1mach(3)))
      xbig = dsqrt (-log(sqrtpi*d1mach(3)))
      sqeps = dsqrt (2.0d0*d1mach(3))
!
 10   y = abs(x)
      if (y.gt.1.d0) go to 20
!
! ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
!
      if (y.le.sqeps) then
         derf = 2.0d0*x*x/sqrtpi
      else
         derf = x*(1.0d0+dcsevl(2.d0*x*x-1.d0,erfcs,nterf))
      end if
!
      return
!
! ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
!
 20   if (y.le.xbig) then
         derf = dsign (1.0d0-derfc(y), x)
      else
         derf = dsign (1.0d0, x)
      end if
!
      return
      end
!XERRWV
      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
!
!     ABSTRACT
!        XERRWV PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER
!        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE
!        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL.
!        (SEE SUBROUTINE XSETF FOR DETAILS.)
!        IN ADDITION, UP TO TWO INTEGER VALUES AND TWO REAL
!        VALUES MAY BE PRINTED ALONG WITH THE MESSAGE.
!
!     DESCRIPTION OF PARAMETERS
!      --INPUT--
!        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED.
!        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG.
!        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE.
!                NERR MUST NOT BE ZERO.
!        LEVEL - ERROR CATEGORY.
!                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR.
!                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS
!                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.)
!                =0 MEANS THIS IS A WARNING MESSAGE ONLY.
!                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE
!                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY
!                   TIMES THIS CALL IS EXECUTED.
!        NI    - NUMBER OF INTEGER VALUES TO BE PRINTED. (O TO 2)
!        I1    - FIRST INTEGER VALUE.
!        I2    - SECOND INTEGER VALUE.
!        NR    - NUMBER OF REAL VALUES TO BE PRINTED. (0 TO 2)
!        R1    - FIRST REAL VALUE.
!        R2    - SECOND REAL VALUE.
!
!     EXAMPLES
!        CALL XERROR('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
!    &   1,NUM,0,0,0.,0.)
!        CALL XERRWV(  'QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM
!    & (R2).',54,77,1,0,0,0,2,ERRREQ,ERRMIN)
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  19 MAR 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real r1,r2
      integer i1,i2,level,nerr,ni,nmessg,nr
!
!  ARRAY ARGUMENTS
      character messg(nmessg)*4
!
!  LOCAL SCALARS
     integer ifatal,iunit,junk,kdummy,kount,kunit,lerr,lkntrl,llevel,&
     &   lmessg,maxmes,mkntrl,nunit
      character lfirst*4
!
!  LOCAL ARRAYS
      integer lun(5)
!
!  EXTERNAL FUNCTIONS
!      INTEGER I1MACH,J4SAVE
!       EXTERNAL I1MACH,J4SAVE
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL FDUMP,XERABT,XERCTL,XERPRT,XERSAV,XGETUA
!
!  INTRINSIC FUNCTIONS
      intrinsic iabs,max,min
!
!     GET FLAGS
      lkntrl = j4save(2,0,.false.)
      maxmes = j4save(4,0,.false.)
!     CHECK FOR VALID INPUT
     if ((nmessg.gt.0).and.(nerr.ne.0).and.&
     &    (level.ge.(-1)).and.(level.le.2)) go to 10
         if (lkntrl.gt.0) call xerprt('FATAL ERROR IN...',17)
         call xerprt('XERROR -- INVALID INPUT',23)
         if (lkntrl.gt.0) call fdump
        if (lkntrl.gt.0) call xerprt('JOB ABORT DUE TO FATAL ERROR.',&
     &   29)
         if (lkntrl.gt.0) call xersav('    ',0,0,0,kdummy)
         call xerabt('XERROR -- INVALID INPUT',23)
         return
   10 continue
!     RECORD MESSAGE
      junk = j4save(1,nerr,.true.)
      call xersav(messg,nmessg,nerr,level,kount)
!     LET USER OVERRIDE
      lfirst = messg(1)
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
!     RESET TO ORIGINAL VALUES
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max(-2,min(2,lkntrl))
      mkntrl = iabs(lkntrl)
!     DECIDE WHETHER TO PRINT MESSAGE
      if ((llevel.lt.2).and.(lkntrl.eq.0)) go to 100
     if (((llevel.eq.(-1)).and.(kount.gt.min(1,maxmes)))&
    &.or.((llevel.eq.0)   .and.(kount.gt.maxmes))&
    &.or.((llevel.eq.1)   .and.(kount.gt.maxmes).and.(mkntrl.eq.1))&
     &.or.((llevel.eq.2)   .and.(kount.gt.max(1,maxmes)))) go to 100
         if (lkntrl.le.0) go to 20
            call xerprt('    ',1)
!           INTRODUCTION
           if (llevel.eq.(-1)) call xerprt&
     &('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            if (llevel.eq.0) call xerprt('WARNING IN...',13)
           if (llevel.eq.1) call xerprt&
     &      ('RECOVERABLE ERROR IN...',23)
            if (llevel.eq.2) call xerprt('FATAL ERROR IN...',17)
   20    continue
!        MESSAGE
         call xerprt(messg,lmessg)
         call xgetua(lun,nunit)
         do 50 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
            if (ni.ge.1) write (iunit,22) i1
            if (ni.ge.2) write (iunit,23) i2
            if (nr.ge.1) write (iunit,24) r1
            if (nr.ge.2) write (iunit,25) r2
   22       format (11x,'IN ABOVE MESSAGE, I1=',i10)
   23       format (11x,'IN ABOVE MESSAGE, I2=',i10)
   24       format (11x,'IN ABOVE MESSAGE, R1=',e20.10)
   25       format (11x,'IN ABOVE MESSAGE, R2=',e20.10)
            if (lkntrl.le.0) go to 40
!              ERROR NUMBER
               write (iunit,30) lerr
   30          format (' ERROR NUMBER =',i10)
   40       continue
   50    continue
!        TRACE-BACK
         if (lkntrl.gt.0) call fdump
  100 continue
      ifatal = 0
     if ((llevel.eq.2).or.((llevel.eq.1).and.(mkntrl.eq.2)))&
     &ifatal = 1
!     QUIT HERE IF MESSAGE IS NOT FATAL
      if (ifatal.le.0) return
      if ((lkntrl.le.0).or.(kount.gt.max(1,maxmes))) go to 120
!        PRINT REASON FOR ABORT
        if (llevel.eq.1) call xerprt&
     &   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
        if (llevel.eq.2) call xerprt&
     &   ('JOB ABORT DUE TO FATAL ERROR.',29)
!        PRINT ERROR SUMMARY
         call xersav('    ',-1,0,0,kdummy)
  120 continue
!     ABORT
      if ((llevel.eq.2).and.(kount.gt.max(1,maxmes))) lmessg = 0
      call xerabt(messg,lmessg)
      return
      end
!D9LGIT
      double precision function d9lgit(a,x,algap1)
!***BEGIN PROLOGUE  D9LGIT
!***DATE WRITTEN   770701   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  C7E
!***KEYWORDS  DOUBLE PRECISION,GAMMA,INCOMPLETE GAMMA FUNCTION,
!             LOGARITHM,SPECIAL FUNCTION,TRICOMI
!***AUTHOR  FULLERTON, W., (LANL)
!***PURPOSE  COMPUTES  THE LOG OF TRICOMI'S INCOMPLETE GAMMA FUNCTION
!            WITH PERRON'S CONTINUED FRACTION FOR LARGE X AND A .GE. X.
!***DESCRIPTION
!
! COMPUTE THE LOG OF TRICOMI'S INCOMPLETE GAMMA FUNCTION WITH PERRON'S
! CONTINUED FRACTION FOR LARGE X AND FOR A .GE. X.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH,XERROR
!***END PROLOGUE  D9LGIT
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,algap1,x
!
!  LOCAL SCALARS
      double precision a1x,ax,eps,fk,hstar,p,r,s,sqeps,t
      integer k
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH
!       EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsqrt,log
!
      data eps, sqeps / 2*0.d0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIT
      if (eps.ne.0.d0) go to 10
      eps = 0.5d0*d1mach(3)
      sqeps = dsqrt (d1mach(4))
!
10   if (x.le.0.d0 .or. a.lt.x) call xerror ( 'D9LGIT  X SHOULD BE GT 0&
     &.0 AND LE A', 35, 2, 2)
!
      ax = a + x
      a1x = ax + 1.0d0
      r = 0.d0
      p = 1.d0
      s = p
      do 20 k=1,200
        fk = k
        t = (a+fk)*x*(1.d0+r)
        r = t/((ax+fk)*(a1x+fk)-t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 30
 20   continue
     call xerror ( 'D9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR&
     &ACTION', 57, 3, 2)
!
 30   hstar = 1.0d0 - x*s/a1x
     if (hstar.lt.sqeps) call xerror ( 'D9LGIT  RESULT LESS THAN HALF P&
     &RECISION', 39, 1, 1)
!
      d9lgit = -x - algap1 - log(hstar)
      return
!
      end
!SETERR
      subroutine seterr(messg,nmessg,nerr,iopt)
!
!  SETERR SETS LERROR = NERR, OPTIONALLY PRINTS THE MESSAGE AND DUMPS
!  ACCORDING TO THE FOLLOWING RULES...
!
!    IF IOPT = 1 AND RECOVERING      - JUST REMEMBER THE ERROR.
!    IF IOPT = 1 AND NOT RECOVERING  - PRINT AND STOP.
!    IF IOPT = 2                     - PRINT, DUMP AND STOP.
!
!  INPUT
!
!    MESSG  - THE ERROR MESSAGE.
!    NMESSG - THE LENGTH OF THE MESSAGE, IN CHARACTERS.
!    NERR   - THE ERROR NUMBER. MUST HAVE NERR NON-ZERO.
!    IOPT   - THE OPTION. MUST HAVE IOPT=1 OR 2.
!
!  ERROR STATES -
!
!    1 - MESSAGE LENGTH NOT POSITIVE.
!    2 - CANNOT HAVE NERR=0.
!    3 - AN UNRECOVERED ERROR FOLLOWED BY ANOTHER ERROR.
!    4 - BAD VALUE FOR IOPT.
!
!  ONLY THE FIRST 72 CHARACTERS OF THE MESSAGE ARE PRINTED.
!
!  THE ERROR HANDLER CALLS A SUBROUTINE NAMED FDUMP TO PRODUCE A
!  SYMBOLIC DUMP. TO COMPLETE THE PACKAGE, A DUMMY VERSION OF FDUMP
!  IS SUPPLIED, BUT IT SHOULD BE REPLACED BY A LOCALLY WRITTEN VERSION
!  WHICH AT LEAST GIVES A TRACE-BACK.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer iopt,nerr,nmessg
!
!  ARRAY ARGUMENTS
      character messg(nmessg)*4
!
!  LOCAL SCALARS
      integer itemp,iwunit,nw
!
!  EXTERNAL FUNCTIONS
!      INTEGER I1MACH,I8SAVE
!       EXTERNAL I1MACH,I8SAVE
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL E9RINT,EPRINT,FDUMP
!
!  INTRINSIC FUNCTIONS
      intrinsic min
!
!
!  THE UNIT FOR ERROR MESSAGES.
!
      iwunit=i1mach(4)
!
      if (nmessg.ge.1) go to 10
!
!  A MESSAGE OF NON-POSITIVE LENGTH IS FATAL.
!
        write(iwunit,9000)
 9000   format('1ERROR    1 IN SETERR - MESSAGE LENGTH NOT POSITIVE.')
        go to 60
!
!  NW IS THE NUMBER OF WORDS THE MESSAGE OCCUPIES.
!
 10   nw=(min(nmessg,72)-1)/i1mach(6)+1
!
      if (nerr.ne.0) go to 20
!
!  CANNOT TURN THE ERROR STATE OFF USING SETERR.
!
        write(iwunit,9001)
9001   format('1ERROR    2 IN SETERR - CANNOT HAVE NERR=0'//&
     &         ' THE CURRENT ERROR MESSAGE FOLLOWS'///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
!
!  SET LERROR AND TEST FOR A PREVIOUS UNRECOVERED ERROR.
!
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
!
        write(iwunit,9002)
9002   format('1ERROR    3 IN SETERR -',&
    &         ' AN UNRECOVERED ERROR FOLLOWED BY ANOTHER ERROR.'//&
     &         ' THE PREVIOUS AND CURRENT ERROR MESSAGES FOLLOW.'///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
!
!  SAVE THIS MESSAGE IN CASE IT IS NOT RECOVERED FROM PROPERLY.
!
 30   call e9rint(messg,nw,nerr,.true.)
!
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
!
!  MUST HAVE IOPT = 1 OR 2.
!
        write(iwunit,9003)
9003   format('1ERROR    4 IN SETERR - BAD VALUE FOR IOPT'//&
     &         ' THE CURRENT ERROR MESSAGE FOLLOWS'///)
        go to 50
!
!  TEST FOR RECOVERY.
!
 40   if (iopt.eq.2) go to 50
!
      if (i8save(2,0,.false.).eq.1) return
!
      call eprint
      stop
!
 50   call eprint
 60   call fdump
      stop
!
      end
!FDUMP
      subroutine fdump
!  THIS IS A DUMMY ROUTINE TO BE SENT OUT ON
!  THE PORT SEDIT TAPE
!
      return
      end
!R9LGIT
      real function r9lgit (a, x, algap1)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG OF TRICOMI-S INCOMPLETE GAMMA FUNCTION WITH PERRON-S
! CONTINUED FRACTION FOR LARGE X AND FOR A .GE. X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,algap1,x
!
!  LOCAL SCALARS
      real a1x,ax,eps,fk,hstar,p,r,s,sqeps,t
      integer k
!
!  EXTERNAL FUNCTIONS
!      REAL R1MACH
!       EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log,sqrt
!
      data eps, sqeps / 2*0.0 /
!
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
      if (sqeps.eq.0.0) sqeps = sqrt(r1mach(4))
!
     if (x.le.0.0 .or. a.lt.x) call xerror (&
     &  'R9LGIT  X SHOULD BE GT 0.0 AND LE A', 35, 2, 2)
!
      ax = a + x
      a1x = ax + 1.0
      r = 0.0
      p = 1.0
      s = p
      do 20 k=1,200
        fk = k
        t = (a+fk)*x*(1.0+r)
        r = t/((ax+fk)*(a1x+fk)-t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 30
 20   continue
     call xerror (  'R9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED F&
     &RACTION', 57, 3, 2)
!
 30   hstar = 1.0 - x*s/a1x
     if (hstar.lt.sqeps) call xerror (&
     &  'R9LGIT  RESULT LESS THAN HALF PRECISION', 39, 1, 1)
!
      r9lgit = -x - algap1 - log(hstar)
!
      return
      end
!XERCTL
      subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
!
!     ABSTRACT
!        ALLOWS USER CONTROL OVER HANDLING OF INDIVIDUAL ERRORS.
!        JUST AFTER EACH MESSAGE IS RECORDED, BUT BEFORE IT IS
!        PROCESSED ANY FURTHER (I.E., BEFORE IT IS PRINTED OR
!        A DECISION TO ABORT IS MADE) A CALL IS MADE TO XERCTL.
!        IF THE USER HAS PROVIDED HIS OWN VERSION OF XERCTL, HE
!        CAN THEN OVERRIDE THE VALUE OF KONTROL USED IN PROCESSING
!        THIS MESSAGE BY REDEFINING ITS VALUE.
!        KONTRL MAY BE SET TO ANY VALUE FROM -2 TO 2.
!        THE MEANINGS FOR KONTRL ARE THE SAME AS IN XSETF, EXCEPT
!        THAT THE VALUE OF KONTRL CHANGES ONLY FOR THIS MESSAGE.
!        IF KONTRL IS SET TO A VALUE OUTSIDE THE RANGE FROM -2 TO 2,
!        IT WILL BE MOVED BACK INTO THAT RANGE.
!
!     DESCRIPTION OF PARAMETERS
!
!      --INPUT--
!        MESSG1 - THE FIRST WORD (ONLY) OF THE ERROR MESSAGE.
!        NMESSG - SAME AS IN THE CALL TO XERROR OR XERRWV.
!        NERR   - SAME AS IN THE CALL TO XERROR OR XERRWV.
!        LEVEL  - SAME AS IN THE CALL TO XERROR OR XERRWV.
!        KONTRL - THE CURRENT VALUE OF THE CONTROL FLAG AS SET
!                 BY A CALL TO XSETF.
!
!      --OUTPUT--
!        KONTRL - THE NEW VALUE OF KONTRL.  IF KONTRL IS NOT
!                 DEFINED, IT WILL REMAIN AT ITS ORIGINAL VALUE.
!                 THIS CHANGED VALUE OF CONTROL AFFECTS ONLY
!                 THE CURRENT OCCURRENCE OF THE CURRENT MESSAGE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer kontrl,level,nerr,nmessg
      character(len=4),intent(in) :: messg1
!
      return
      end
!S88FMT
      subroutine s88fmt( n, w, ifmt )
!
!     LATEST REVISION  -  OCTOBER 3, 1983  (JRD)
!
!  S88FMT  REPLACES IFMT(1), ... , IFMT(N) WITH
!  THE CHARACTERS CORRESPONDING TO THE N LEAST SIGNIFICANT
!  DIGITS OF W.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer n,w
!
!  ARRAY ARGUMENTS
      character ifmt(n)*4
!
!  LOCAL SCALARS
      integer idigit,nt,wt
!
!  LOCAL ARRAYS
      character digits(10)*4
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!
      data digits( 1) / '0' /
      data digits( 2) / '1' /
      data digits( 3) / '2' /
      data digits( 4) / '3' /
      data digits( 5) / '4' /
      data digits( 6) / '5' /
      data digits( 7) / '6' /
      data digits( 8) / '7' /
      data digits( 9) / '8' /
      data digits(10) / '9' /
!
      nt = n
      wt = w
!
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
!
      end
!ERF
      real function erf (x)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real sqeps,sqrtpi,xbig,y
      integer nterf
!
!  LOCAL ARRAYS
      real erfcs(13)
!
!  EXTERNAL FUNCTIONS
!      REAL CSEVL,ERFC,R1MACH
!      INTEGER INITS
!       EXTERNAL CSEVL,ERFC,R1MACH,INITS
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log,sign,sqrt
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000D+00
!                                        WITH WEIGHTED ERROR   7.10E-18
!                                         LOG WEIGHTED ERROR  17.15
!                               SIGNIFICANT FIGURES REQUIRED  16.31
!                                    DECIMAL PLACES REQUIRED  17.71
!
      data erfcs( 1) /   -.049046121234691808e0 /
      data erfcs( 2) /   -.14226120510371364e0 /
      data erfcs( 3) /    .010035582187599796e0 /
      data erfcs( 4) /   -.000576876469976748e0 /
      data erfcs( 5) /    .000027419931252196e0 /
      data erfcs( 6) /   -.000001104317550734e0 /
      data erfcs( 7) /    .000000038488755420e0 /
      data erfcs( 8) /   -.000000001180858253e0 /
      data erfcs( 9) /    .000000000032334215e0 /
      data erfcs(10) /   -.000000000000799101e0 /
      data erfcs(11) /    .000000000000017990e0 /
      data erfcs(12) /   -.000000000000000371e0 /
      data erfcs(13) /    .000000000000000007e0 /
!
      data sqrtpi /1.7724538509055160e0/
      data nterf, xbig, sqeps / 0, 0., 0./
!
      if (nterf.ne.0) go to 10
      nterf = inits (erfcs, 13, 0.1*r1mach(3))
      xbig = sqrt(-log(sqrtpi*r1mach(3)))
      sqeps = sqrt(2.0*r1mach(3))
!
 10   y = abs(x)
      if (y.gt.1.) go to 20
!
! ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1.
!
      if (y.le.sqeps) then
         erf = 2.0*x/sqrtpi
      else
         erf = x*(1.0 + csevl(2.*x**2-1., erfcs, nterf))
      end if
!
      return
!
! ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1.
!
 20   if (y.le.xbig) then
         erf = sign (1.0-erfc(y), x)
      else
         erf = sign (1.0, x)
      end if
!
      return
      end
!XERABT
      subroutine xerabt(messg,nmessg)
!
!     LATEST REVISION  -  JANUARY 24, 1990 (JRD)
!
!     ABSTRACT
!        ***NOTE*** MACHINE DEPENDENT ROUTINE
!        XERABT ABORTS THE EXECUTION OF THE PROGRAM.
!        THE ERROR MESSAGE CAUSING THE ABORT IS GIVEN IN THE CALLING
!        SEQUENCE IN CASE ONE NEEDS IT FOR PRINTING ON A DAYFILE,
!        FOR EXAMPLE.
!
!     DESCRIPTION OF PARAMETERS
!        MESSG AND NMESSG ARE AS IN XERROR, EXCEPT THAT NMESSG MAY
!        BE ZERO, IN WHICH CASE NO MESSAGE IS BEING SUPPLIED.
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  7 JUNE 1978
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer nmessg
!
!  ARRAY ARGUMENTS
      character messg(1)*4
!
      stop
      end
!ALNREL
      real function alnrel (x)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real xmin
      integer nlnrel
!
!  LOCAL ARRAYS
      real alnrcs(23)
!
!  EXTERNAL FUNCTIONS
!      REAL CSEVL,R1MACH
!      INTEGER INITS
!       EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log,sqrt
!
! SERIES FOR ALNR       ON THE INTERVAL -3.75000D-01 TO  3.75000D-01
!                                        WITH WEIGHTED ERROR   1.93E-17
!                                         LOG WEIGHTED ERROR  16.72
!                               SIGNIFICANT FIGURES REQUIRED  16.44
!                                    DECIMAL PLACES REQUIRED  17.40
!
      data alnrcs( 1) /   1.0378693562743770e0 /
      data alnrcs( 2) /   -.13364301504908918e0 /
      data alnrcs( 3) /    .019408249135520563e0 /
      data alnrcs( 4) /   -.003010755112753577e0 /
      data alnrcs( 5) /    .000486946147971548e0 /
      data alnrcs( 6) /   -.000081054881893175e0 /
      data alnrcs( 7) /    .000013778847799559e0 /
      data alnrcs( 8) /   -.000002380221089435e0 /
      data alnrcs( 9) /    .000000416404162138e0 /
      data alnrcs(10) /   -.000000073595828378e0 /
      data alnrcs(11) /    .000000013117611876e0 /
      data alnrcs(12) /   -.000000002354670931e0 /
      data alnrcs(13) /    .000000000425227732e0 /
      data alnrcs(14) /   -.000000000077190894e0 /
      data alnrcs(15) /    .000000000014075746e0 /
      data alnrcs(16) /   -.000000000002576907e0 /
      data alnrcs(17) /    .000000000000473424e0 /
      data alnrcs(18) /   -.000000000000087249e0 /
      data alnrcs(19) /    .000000000000016124e0 /
      data alnrcs(20) /   -.000000000000002987e0 /
      data alnrcs(21) /    .000000000000000554e0 /
      data alnrcs(22) /   -.000000000000000103e0 /
      data alnrcs(23) /    .000000000000000019e0 /
!
      data nlnrel, xmin /0, 0./
!
      if (nlnrel.ne.0) go to 10
      nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))
      xmin = -1.0 + sqrt(r1mach(4))
!
10   if (x.le.(-1.0)) call xerror (&
     &  'ALNREL  X IS LE -1', 18, 2, 2)
     if (x.lt.xmin) call xerror (&
    &  'ALNREL  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 54,&
     &  1, 1)
!
      if (abs(x).le.0.375) then
         alnrel = x*(1.0-x*csevl(x/0.375,alnrcs,nlnrel))
      else
         alnrel = log (1.0+x)
      end if
!
      return
      end
!INITDS
      integer function initds (dos, nos, eta)
!
! INITIALIZE THE DOUBLE PRECISION ORTHOGONAL SERIES DOS SO THAT INITDS
! IS THE NUMBER OF TERMS NEEDED TO INSURE THE ERROR IS NO LARGER THAN
! ETA.  ORDINARILY ETA WILL BE CHOSEN TO BE ONE-TENTH MACHINE PRECISION.
!
!             INPUT ARGUMENTS --
! DOS    DBLE PREC ARRAY OF NOS COEFFICIENTS IN AN ORTHOGONAL SERIES.
! NOS    NUMBER OF COEFFICIENTS IN DOS.
! ETA    REQUESTED ACCURACY OF SERIES.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real eta
      integer nos
!
!  ARRAY ARGUMENTS
      double precision dos(nos)
!
!  LOCAL SCALARS
      real err
      integer i,ii
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,sngl
!
!
     if (nos.lt.1) call xerror (&
     &  'INITDS  NUMBER OF COEFFICIENTS LT 1', 35, 2, 2)
!
      err = 0.0
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
!
20   if (i.eq.nos) call xerror ('INITDS  ETA MAY BE TOO SMALL', 28,&
     &  1, 2)
      initds = i
!
      return
      end
!DBETAI
      double precision function dbetai (x, pin, qin)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! BASED ON BOSTEN AND BATTISTE, REMARK ON ALGORITHM 179, COMM. ACM,
! V 17, P 156, (1974).
!
!             INPUT ARGUMENTS --
! X      UPPER LIMIT OF INTEGRATION.  X MUST BE IN (0,1) INCLUSIVE.
! P      FIRST BETA DISTRIBUTION PARAMETER.  P MUST BE GT 0.0.
! Q      SECOND BETA DISTRIBUTION PARAMETER.  Q MUST BE GT 0.0.
! BETAI  THE INCOMPLETE BETA FUNCTION RATIO IS THE PROBABILITY THAT A
!        RANDOM VARIABLE FROM A BETA DISTRIBUTION HAVING PARAMETERS
!        P AND Q WILL BE LESS THAN OR EQUAL TO X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision pin,qin,x
!
!  LOCAL SCALARS
     double precision alneps,alnsml,c,eps,fac1,fac2,finsum,p,ps,q,sml,&
     &   term,xb,y
      real p1
      integer i,ib,n
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DLBETA
!       EXTERNAL D1MACH,DLBETA
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dble,exp,float,int,log,max,min,sngl
!
      data             eps, alneps, sml, alnsml / 4*0.0d0 /
!
      if (eps.ne.0.0d0) go to 10
      eps = d1mach(3)
      alneps = log (eps)
      sml = d1mach(1)
      alnsml = log (sml)
!
10   if (x.lt.0.d0 .or. x.gt.1.d0) call xerror (&
     &  'DBETAI  X IS NOT IN THE RANGE (0,1)', 35, 1, 2)
     if (pin.le.0.d0 .or. qin.le.0.d0) call xerror (&
     &  'DBETAI  P AND/OR Q IS LE ZERO', 29, 2, 2)
!
      y = x
      p = pin
      q = qin
      if (q.le.p .and. x.lt.0.8d0) go to 20
      if (x.lt.0.2d0) go to 20
      y = 1.0d0 - y
      p = qin
      q = pin
!
 20   if ((p+q)*y/(p+1.d0).lt.eps) go to 80
!
! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
!
      ps = q - int(q)
      if (ps.eq.0.d0) ps = 1.0d0
      xb = p*log(y) - dlbeta(ps,p) - log(p)
      dbetai = 0.0d0
      if (xb.ge.alnsml) then
         dbetai = exp(xb)
         fac2 = 1.0
         if (ps.ne.1.0d0) then
            fac1 = 1.0
            n = max(alneps/log(y), 4.0d0)
            do 30 i=1,n
               if ((i-ps.eq.0.0d0) .or. (fac1.eq.0.0d0)) then
                  fac1 = 0.0d0
               else
                 if (log(abs(fac1)) + log(abs(i-ps)) + log(y) -&
     &                log(dble(i)) .lt. alnsml) then
                     fac1 = 0.0d0
                  else
                     fac1 = fac1 * (i-ps)*y/i
                  end if
               end if
               fac2 = fac2 + fac1*p/(p+i)
 30         continue
         end if
         dbetai = dbetai*fac2
      end if
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
      if (q.le.1.0d0) go to 70
!
      xb = p*log(y) + q*log(1.0d0-y) - dlbeta(p,q) - log(q)
      ib = max(sngl(xb/alnsml), 0.0)
      term = exp (xb - dble(float(ib))*alnsml )
      c = 1.0d0/(1.d0-y)
      p1 = q*c/(p+q-1.d0)
!
      finsum = 0.0d0
      n = q
      if (q.eq.dble(float(n))) n = n - 1
      do 50 i=1,n
        if (p1.le.1.0d0 .and. term/eps.le.finsum) go to 60
        if (q-i+1.0d0 .eq. 0.0d0) then
          term = 0.0d0
        else
         if (log(abs(q-i+1.0d0)) + log(abs(c)) + log(abs(term)) -&
     &        log(abs(p+q-i)) .lt. alnsml) then
            term = 0.0d0
          else
            term = (q-i+1.0d0)*c*term/(p+q-i)
          end if
        end if
!
        if (term.gt.1.0d0) ib = ib - 1
        if (term.gt.1.0d0) term = term*sml
!
        if (ib.eq.0) finsum = finsum + term
 50   continue
!
 60   dbetai = dbetai + finsum
 70   if (y.ne.x .or. p.ne.pin) dbetai = 1.0d0 - dbetai
      dbetai = max (min (dbetai, 1.0d0), 0.0d0)
      return
!
 80   dbetai = 0.0d0
      xb = p*log(max(y,sml)) - log(p) - dlbeta(p,q)
      if (xb.gt.alnsml .and. y.ne.0.0d0) dbetai = exp(xb)
      if (y.ne.x .or. p.ne.pin) dbetai = 1.0d0 - dbetai
!
      return
      end
!I8SAVE
      integer function i8save(isw,ivalue,set)
!
!  IF (ISW = 1) I8SAVE RETURNS THE CURRENT ERROR NUMBER AND
!               SETS IT TO IVALUE IF SET = .TRUE. .
!
!  IF (ISW = 2) I8SAVE RETURNS THE CURRENT RECOVERY SWITCH AND
!               SETS IT TO IVALUE IF SET = .TRUE. .
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer isw,ivalue
      logical set
!
!  LOCAL SCALARS
      integer lerror,lrecov
!
!  LOCAL ARRAYS
      integer iparam(2)
!
!  EQUIVALENCES
      equivalence (iparam(1),lerror), (iparam(2),lrecov)
!
!
!  START EXECUTION ERROR FREE AND WITH RECOVERY TURNED OFF.
!
      data lerror/0/ , lrecov/2/
!
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
!
      return
!
      end
!DGAMMA
      double precision function dgamma (x)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision dxrel,pi,sinpiy,sq2pil,xmax,xmin,y
      integer i,n,ngam
!
!  LOCAL ARRAYS
      double precision gamcs(42)
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,D9LGMC,DCSEVL
!      INTEGER INITDS
!       EXTERNAL D1MACH,D9LGMC,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DGAMLM,XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dble,dsqrt,exp,float,int,log,sin,sngl
!
!
! SERIES FOR GAM        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   5.79E-32
!                                         LOG WEIGHTED ERROR  31.24
!                               SIGNIFICANT FIGURES REQUIRED  30.00
!                                    DECIMAL PLACES REQUIRED  32.05
!
      data gamcs(  1) / +.8571195590989331421920062399942d-2      /
      data gamcs(  2) / +.4415381324841006757191315771652d-2      /
      data gamcs(  3) / +.5685043681599363378632664588789d-1      /
      data gamcs(  4) / -.4219835396418560501012500186624d-2      /
      data gamcs(  5) / +.1326808181212460220584006796352d-2      /
      data gamcs(  6) / -.1893024529798880432523947023886d-3      /
      data gamcs(  7) / +.3606925327441245256578082217225d-4      /
      data gamcs(  8) / -.6056761904460864218485548290365d-5      /
      data gamcs(  9) / +.1055829546302283344731823509093d-5      /
      data gamcs( 10) / -.1811967365542384048291855891166d-6      /
      data gamcs( 11) / +.3117724964715322277790254593169d-7      /
      data gamcs( 12) / -.5354219639019687140874081024347d-8      /
      data gamcs( 13) / +.9193275519859588946887786825940d-9      /
      data gamcs( 14) / -.1577941280288339761767423273953d-9      /
      data gamcs( 15) / +.2707980622934954543266540433089d-10     /
      data gamcs( 16) / -.4646818653825730144081661058933d-11     /
      data gamcs( 17) / +.7973350192007419656460767175359d-12     /
      data gamcs( 18) / -.1368078209830916025799499172309d-12     /
      data gamcs( 19) / +.2347319486563800657233471771688d-13     /
      data gamcs( 20) / -.4027432614949066932766570534699d-14     /
      data gamcs( 21) / +.6910051747372100912138336975257d-15     /
      data gamcs( 22) / -.1185584500221992907052387126192d-15     /
      data gamcs( 23) / +.2034148542496373955201026051932d-16     /
      data gamcs( 24) / -.3490054341717405849274012949108d-17     /
      data gamcs( 25) / +.5987993856485305567135051066026d-18     /
      data gamcs( 26) / -.1027378057872228074490069778431d-18     /
      data gamcs( 27) / +.1762702816060529824942759660748d-19     /
      data gamcs( 28) / -.3024320653735306260958772112042d-20     /
      data gamcs( 29) / +.5188914660218397839717833550506d-21     /
      data gamcs( 30) / -.8902770842456576692449251601066d-22     /
      data gamcs( 31) / +.1527474068493342602274596891306d-22     /
      data gamcs( 32) / -.2620731256187362900257328332799d-23     /
      data gamcs( 33) / +.4496464047830538670331046570666d-24     /
      data gamcs( 34) / -.7714712731336877911703901525333d-25     /
      data gamcs( 35) / +.1323635453126044036486572714666d-25     /
      data gamcs( 36) / -.2270999412942928816702313813333d-26     /
      data gamcs( 37) / +.3896418998003991449320816639999d-27     /
      data gamcs( 38) / -.6685198115125953327792127999999d-28     /
      data gamcs( 39) / +.1146998663140024384347613866666d-28     /
      data gamcs( 40) / -.1967938586345134677295103999999d-29     /
      data gamcs( 41) / +.3376448816585338090334890666666d-30     /
      data gamcs( 42) / -.5793070335782135784625493333333d-31     /
!
      data pi / 3.14159265358979323846264338327950d0 /
! SQ2PIL IS 0.5*LOG(2*PI) = LOG(SQRT(2*PI))
      data sq2pil / 0.91893853320467274178032973640562d0 /
      data ngam, xmin, xmax, dxrel / 0, 3*0.d0 /
!
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
!
      call dgamlm (xmin, xmax)
      dxrel = dsqrt (d1mach(4))
!
 10   y = abs(x)
      if (y.gt.10.d0) go to 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
!
      if (n.gt.0) go to 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      n = -n
      if (x.eq.0.d0) call xerror ('DGAMMA  X IS 0', 14, 4, 2)
     if (x.lt.0.0 .and. x+dble(float(n-2)).eq.0.d0) call xerror (&
     &  'DGAMMA  X IS A NEGATIVE INTEGER', 31, 4, 2)
     if (x.lt.(-0.5d0) .and. abs((x-int(x-0.5d0))/x).lt.dxrel) call&
    &  xerror (  'DGAMMA  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR N&
     &EGATIVE INTEGER', 68, 1, 1)
!
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
50   if (x.gt.xmax) call xerror ('DGAMMA  X SO BIG GAMMA OVERFLOWS',&
     &  32, 3, 2)
!
      dgamma = 0.d0
     if (x.lt.xmin) call xerror ('DGAMMA  X SO SMALL GAMMA UNDERFLOWS',&
     &  35, 2, 1)
      if (x.lt.xmin) return
!
      dgamma = exp ((y-0.5d0)*log(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
!
     if (abs((x-int(x-0.5d0))/x).lt.dxrel) call xerror (&
    &  'DGAMMA  ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER',&
     &  61, 1, 1)
!
      sinpiy = sin (pi*y)
     if (sinpiy.eq.0.d0) call xerror (&
     &  'DGAMMA  X IS A NEGATIVE INTEGER', 31, 4, 2)
!
      dgamma = -pi/(y*sinpiy*dgamma)
!
      return
      end
!DERFC
      double precision function derfc (x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision sqeps,sqrtpi,xmax,xsml,y
      real eta
      integer nterc2,nterf,nterfc
!
!  LOCAL ARRAYS
      double precision erc2cs(49),erfccs(59),erfcs(21)
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DCSEVL
!      INTEGER INITDS
!       EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsqrt,exp,log,sngl
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   1.28E-32
!                                         LOG WEIGHTED ERROR  31.89
!                               SIGNIFICANT FIGURES REQUIRED  31.05
!                                    DECIMAL PLACES REQUIRED  32.55
!
      data erfcs(  1) / -.49046121234691808039984544033376d-1     /
      data erfcs(  2) / -.14226120510371364237824741899631d+0     /
      data erfcs(  3) / +.10035582187599795575754676712933d-1     /
      data erfcs(  4) / -.57687646997674847650827025509167d-3     /
      data erfcs(  5) / +.27419931252196061034422160791471d-4     /
      data erfcs(  6) / -.11043175507344507604135381295905d-5     /
      data erfcs(  7) / +.38488755420345036949961311498174d-7     /
      data erfcs(  8) / -.11808582533875466969631751801581d-8     /
      data erfcs(  9) / +.32334215826050909646402930953354d-10    /
      data erfcs( 10) / -.79910159470045487581607374708595d-12    /
      data erfcs( 11) / +.17990725113961455611967245486634d-13    /
      data erfcs( 12) / -.37186354878186926382316828209493d-15    /
      data erfcs( 13) / +.71035990037142529711689908394666d-17    /
      data erfcs( 14) / -.12612455119155225832495424853333d-18    /
      data erfcs( 15) / +.20916406941769294369170500266666d-20    /
      data erfcs( 16) / -.32539731029314072982364160000000d-22    /
      data erfcs( 17) / +.47668672097976748332373333333333d-24    /
      data erfcs( 18) / -.65980120782851343155199999999999d-26    /
      data erfcs( 19) / +.86550114699637626197333333333333d-28    /
      data erfcs( 20) / -.10788925177498064213333333333333d-29    /
      data erfcs( 21) / +.12811883993017002666666666666666d-31    /
!
! SERIES FOR ERC2       ON THE INTERVAL  2.50000E-01 TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   2.67E-32
!                                         LOG WEIGHTED ERROR  31.57
!                               SIGNIFICANT FIGURES REQUIRED  30.31
!                                    DECIMAL PLACES REQUIRED  32.42
!
      data erc2cs(  1) / -.6960134660230950112739150826197d-1      /
      data erc2cs(  2) / -.4110133936262089348982212084666d-1      /
      data erc2cs(  3) / +.3914495866689626881561143705244d-2      /
      data erc2cs(  4) / -.4906395650548979161280935450774d-3      /
      data erc2cs(  5) / +.7157479001377036380760894141825d-4      /
      data erc2cs(  6) / -.1153071634131232833808232847912d-4      /
      data erc2cs(  7) / +.1994670590201997635052314867709d-5      /
      data erc2cs(  8) / -.3642666471599222873936118430711d-6      /
      data erc2cs(  9) / +.6944372610005012589931277214633d-7      /
      data erc2cs( 10) / -.1371220902104366019534605141210d-7      /
      data erc2cs( 11) / +.2788389661007137131963860348087d-8      /
      data erc2cs( 12) / -.5814164724331161551864791050316d-9      /
      data erc2cs( 13) / +.1238920491752753181180168817950d-9      /
      data erc2cs( 14) / -.2690639145306743432390424937889d-10     /
      data erc2cs( 15) / +.5942614350847910982444709683840d-11     /
      data erc2cs( 16) / -.1332386735758119579287754420570d-11     /
      data erc2cs( 17) / +.3028046806177132017173697243304d-12     /
      data erc2cs( 18) / -.6966648814941032588795867588954d-13     /
      data erc2cs( 19) / +.1620854541053922969812893227628d-13     /
      data erc2cs( 20) / -.3809934465250491999876913057729d-14     /
      data erc2cs( 21) / +.9040487815978831149368971012975d-15     /
      data erc2cs( 22) / -.2164006195089607347809812047003d-15     /
      data erc2cs( 23) / +.5222102233995854984607980244172d-16     /
      data erc2cs( 24) / -.1269729602364555336372415527780d-16     /
      data erc2cs( 25) / +.3109145504276197583836227412951d-17     /
      data erc2cs( 26) / -.7663762920320385524009566714811d-18     /
      data erc2cs( 27) / +.1900819251362745202536929733290d-18     /
      data erc2cs( 28) / -.4742207279069039545225655999965d-19     /
      data erc2cs( 29) / +.1189649200076528382880683078451d-19     /
      data erc2cs( 30) / -.3000035590325780256845271313066d-20     /
      data erc2cs( 31) / +.7602993453043246173019385277098d-21     /
      data erc2cs( 32) / -.1935909447606872881569811049130d-21     /
      data erc2cs( 33) / +.4951399124773337881000042386773d-22     /
      data erc2cs( 34) / -.1271807481336371879608621989888d-22     /
      data erc2cs( 35) / +.3280049600469513043315841652053d-23     /
      data erc2cs( 36) / -.8492320176822896568924792422399d-24     /
      data erc2cs( 37) / +.2206917892807560223519879987199d-24     /
      data erc2cs( 38) / -.5755617245696528498312819507199d-25     /
      data erc2cs( 39) / +.1506191533639234250354144051199d-25     /
      data erc2cs( 40) / -.3954502959018796953104285695999d-26     /
      data erc2cs( 41) / +.1041529704151500979984645051733d-26     /
      data erc2cs( 42) / -.2751487795278765079450178901333d-27     /
      data erc2cs( 43) / +.7290058205497557408997703680000d-28     /
      data erc2cs( 44) / -.1936939645915947804077501098666d-28     /
      data erc2cs( 45) / +.5160357112051487298370054826666d-29     /
      data erc2cs( 46) / -.1378419322193094099389644800000d-29     /
      data erc2cs( 47) / +.3691326793107069042251093333333d-30     /
      data erc2cs( 48) / -.9909389590624365420653226666666d-31     /
      data erc2cs( 49) / +.2666491705195388413323946666666d-31     /
!
! SERIES FOR ERFC       ON THE INTERVAL  0.          TO  2.50000E-01
!                                        WITH WEIGHTED ERROR   1.53E-31
!                                         LOG WEIGHTED ERROR  30.82
!                               SIGNIFICANT FIGURES REQUIRED  29.47
!                                    DECIMAL PLACES REQUIRED  31.70
!
      data erfccs(  1) / +.715179310202924774503697709496d-1        /
      data erfccs(  2) / -.265324343376067157558893386681d-1        /
      data erfccs(  3) / +.171115397792085588332699194606d-2        /
      data erfccs(  4) / -.163751663458517884163746404749d-3        /
      data erfccs(  5) / +.198712935005520364995974806758d-4        /
      data erfccs(  6) / -.284371241276655508750175183152d-5        /
      data erfccs(  7) / +.460616130896313036969379968464d-6        /
      data erfccs(  8) / -.822775302587920842057766536366d-7        /
      data erfccs(  9) / +.159214187277090112989358340826d-7        /
      data erfccs( 10) / -.329507136225284321486631665072d-8        /
      data erfccs( 11) / +.722343976040055546581261153890d-9        /
      data erfccs( 12) / -.166485581339872959344695966886d-9        /
      data erfccs( 13) / +.401039258823766482077671768814d-10       /
      data erfccs( 14) / -.100481621442573113272170176283d-10       /
      data erfccs( 15) / +.260827591330033380859341009439d-11       /
      data erfccs( 16) / -.699111056040402486557697812476d-12       /
      data erfccs( 17) / +.192949233326170708624205749803d-12       /
      data erfccs( 18) / -.547013118875433106490125085271d-13       /
      data erfccs( 19) / +.158966330976269744839084032762d-13       /
      data erfccs( 20) / -.472689398019755483920369584290d-14       /
      data erfccs( 21) / +.143587337678498478672873997840d-14       /
      data erfccs( 22) / -.444951056181735839417250062829d-15       /
      data erfccs( 23) / +.140481088476823343737305537466d-15       /
      data erfccs( 24) / -.451381838776421089625963281623d-16       /
      data erfccs( 25) / +.147452154104513307787018713262d-16       /
      data erfccs( 26) / -.489262140694577615436841552532d-17       /
      data erfccs( 27) / +.164761214141064673895301522827d-17       /
      data erfccs( 28) / -.562681717632940809299928521323d-18       /
      data erfccs( 29) / +.194744338223207851429197867821d-18       /
      data erfccs( 30) / -.682630564294842072956664144723d-19       /
      data erfccs( 31) / +.242198888729864924018301125438d-19       /
      data erfccs( 32) / -.869341413350307042563800861857d-20       /
      data erfccs( 33) / +.315518034622808557122363401262d-20       /
      data erfccs( 34) / -.115737232404960874261239486742d-20       /
      data erfccs( 35) / +.428894716160565394623737097442d-21       /
      data erfccs( 36) / -.160503074205761685005737770964d-21       /
      data erfccs( 37) / +.606329875745380264495069923027d-22       /
      data erfccs( 38) / -.231140425169795849098840801367d-22       /
      data erfccs( 39) / +.888877854066188552554702955697d-23       /
      data erfccs( 40) / -.344726057665137652230718495566d-23       /
      data erfccs( 41) / +.134786546020696506827582774181d-23       /
      data erfccs( 42) / -.531179407112502173645873201807d-24       /
      data erfccs( 43) / +.210934105861978316828954734537d-24       /
      data erfccs( 44) / -.843836558792378911598133256738d-25       /
      data erfccs( 45) / +.339998252494520890627359576337d-25       /
      data erfccs( 46) / -.137945238807324209002238377110d-25       /
      data erfccs( 47) / +.563449031183325261513392634811d-26       /
      data erfccs( 48) / -.231649043447706544823427752700d-26       /
      data erfccs( 49) / +.958446284460181015263158381226d-27       /
      data erfccs( 50) / -.399072288033010972624224850193d-27       /
      data erfccs( 51) / +.167212922594447736017228709669d-27       /
      data erfccs( 52) / -.704599152276601385638803782587d-28       /
      data erfccs( 53) / +.297976840286420635412357989444d-28       /
      data erfccs( 54) / -.126252246646061929722422632994d-28       /
      data erfccs( 55) / +.539543870454248793985299653154d-29       /
      data erfccs( 56) / -.238099288253145918675346190062d-29       /
      data erfccs( 57) / +.109905283010276157359726683750d-29       /
      data erfccs( 58) / -.486771374164496572732518677435d-30       /
      data erfccs( 59) / +.152587726411035756763200828211d-30       /
!
      data sqrtpi / 1.77245385090551602729816748334115d0 /
      data nterf, nterfc, nterc2, xsml, xmax, sqeps / 3*0, 3*0.d0 /
!
      if (nterf.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nterf = initds (erfcs, 21, eta)
      nterfc = initds (erfccs, 59, eta)
      nterc2 = initds (erc2cs, 49, eta)
!
      xsml = -dsqrt (-log(sqrtpi*d1mach(3)))
      xmax = dsqrt (-log(sqrtpi*d1mach(1)) )
      xmax = xmax - 0.5d0*log(xmax)/xmax - 0.01d0
      sqeps = dsqrt (2.0d0*d1mach(3))
!
 10   if (x.gt.xsml) go to 20
!
! ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
!
      derfc = 2.0d0
      return
!
 20   if (x.gt.xmax) go to 40
      y = abs(x)
      if (y.gt.1.0d0) go to 30
!
! ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
!
      if (y.lt.sqeps) then
         derfc = 1.0d0 - 2.0d0*x/sqrtpi
      else
         derfc = 1.0d0 - x*(1.0d0+dcsevl(2.0d0*x*x-1.0d0,erfcs,nterf))
      end if
!
      return
!
! ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
!
 30   y = y*y
      if (y.le.4.0d0) then
        derfc = exp(-y)/abs(x) *&
     &           (0.5d0 + dcsevl((8.0d0/y-5.0d0)/3.0d0,erc2cs,nterc2))
      else
        derfc = exp(-y)/abs(x) *&
     &           (0.5d0 + dcsevl(8.0d0/y-1.0d0,erfccs,nterfc))
      end if
      if (x.lt.0.d0) derfc = 2.0d0 - derfc
      return
!
 40   call xerror ('DERFC   X SO BIG ERFC UNDERFLOWS', 32, 1, 1)
      derfc = 0.d0
      return
!
      end
!DGAMLM
      subroutine dgamlm (xmin, xmax)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! CALCULATE THE MINIMUM AND MAXIMUM LEGAL BOUNDS FOR X IN GAMMA(X).
! XMIN AND XMAX ARE NOT THE ONLY BOUNDS, BUT THEY ARE THE ONLY NON-
! TRIVIAL ONES TO CALCULATE.
!
!             OUTPUT ARGUMENTS --
! XMIN   DBLE PREC MINIMUM LEGAL VALUE OF X IN GAMMA(X).  ANY SMALLER
!        VALUE OF X MIGHT RESULT IN UNDERFLOW.
! XMAX   DBLE PREC MAXIMUM LEGAL VALUE OF X IN GAMMA(X).  ANY LARGER
!        VALUE OF X MIGHT CAUSE OVERFLOW.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision xmax,xmin
!
!  LOCAL SCALARS
      double precision alnbig,alnsml,xln,xold
      integer i
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH
!       EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,log,max
!
!
      alnsml = log(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
       xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)&
     &    / (xmin*xln+0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
      call xerror ('DGAMLM  UNABLE TO FIND XMIN', 27, 1, 2)
!
 20   xmin = -xmin + 0.01d0
!
      alnbig = log (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
       xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)&
     &    / (xmax*xln-0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
      call xerror ('DGAMLM  UNABLE TO FIND XMAX', 27, 2, 2)
!
 40   xmax = xmax - 0.01d0
      xmin = max (xmin, -xmax+1.d0)
!
      return
      end
!ALBETA
      real function albeta (a, b)
! JULY 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,b
!
!  LOCAL SCALARS
      real corr,p,q,sq2pil
!
!  EXTERNAL FUNCTIONS
!      REAL ALNGAM,ALNREL,GAMMA,R9LGMC
!       EXTERNAL ALNGAM,ALNREL,GAMMA,R9LGMC
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic log,max,min
!
      data sq2pil / 0.91893853320467274e0 /
!
      p = min (a, b)
      q = max (a, b)
!
     if (p.le.0.0) call xerror (&
     &  'ALBETA  BOTH ARGUMENTS MUST BE GT ZERO', 38, 1, 2)
      if (p.ge.10.0) go to 30
      if (q.ge.10.0) go to 20
!
! P AND Q ARE SMALL.
!
      albeta = log(gamma(p) * (gamma(q)/gamma(p+q)) )
      return
!
! P IS SMALL, BUT Q IS BIG.
!
 20   corr = r9lgmc(q) - r9lgmc(p+q)
     albeta = alngam(p) + corr + p - p*log(p+q) +&
     &  (q-0.5)*alnrel(-p/(p+q))
      return
!
! P AND Q ARE BIG.
!
 30   corr = r9lgmc(p) + r9lgmc(q) - r9lgmc(p+q)
     albeta = -0.5*log(q) + sq2pil + corr + (p-0.5)*log(p/(p+q))&
     &  + q*alnrel(-p/(p+q))
      return
!
      end
!XERSAV
      subroutine xersav(messg,nmessg,nerr,level,icount)
!
!     ABSTRACT
!        RECORD THAT THIS ERROR OCCURRED.
!
!     DESCRIPTION OF PARAMETERS
!     --INPUT--
!       MESSG, NMESSG, NERR, LEVEL ARE AS IN XERROR,
!       EXCEPT THAT WHEN NMESSG=0 THE TABLES WILL BE
!       DUMPED AND CLEARED, AND WHEN NMESSG IS LESS THAN ZERO THE
!       TABLES WILL BE DUMPED AND NOT CLEARED.
!     --OUTPUT--
!       ICOUNT WILL BE THE NUMBER OF TIMES THIS MESSAGE HAS
!       BEEN SEEN, OR ZERO IF THE TABLE HAS OVERFLOWED AND
!       DOES NOT CONTAIN THIS MESSAGE SPECIFICALLY.
!       WHEN NMESSG=0, ICOUNT WILL NOT BE ALTERED.
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  19 MAR 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer icount,level,nerr,nmessg
!
!  ARRAY ARGUMENTS
      character messg(nmessg)*4
!
!  LOCAL SCALARS
      integer i,ii,iunit,kountx,kunit,nchar,ncol,nunit
!
!  LOCAL ARRAYS
      integer kount(10),levtab(10),lun(5),nertab(10)
      character f(17)*1,mestab(10)*4
!
!  EXTERNAL FUNCTIONS
!      INTEGER I1MACH
!       EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL S88FMT,XGETUA
!
!     NEXT THREE DATA STATEMENTS ARE NEEDED MERELY TO SATISFY
!     CERTAIN CONVENTIONS FOR COMPILERS WHICH DYNAMICALLY
!     ALLOCATE STORAGE.
     data mestab(1),mestab(2),mestab(3),mestab(4),mestab(5),&
    &     mestab(6),mestab(7),mestab(8),mestab(9),mestab(10)&
     &     /'0','0','0','0','0','0','0','0','0','0'/
     data nertab(1),nertab(2),nertab(3),nertab(4),nertab(5),&
    &     nertab(6),nertab(7),nertab(8),nertab(9),nertab(10)&
     &     /0,0,0,0,0,0,0,0,0,0/
     data levtab(1),levtab(2),levtab(3),levtab(4),levtab(5),&
    &     levtab(6),levtab(7),levtab(8),levtab(9),levtab(10)&
     &     /0,0,0,0,0,0,0,0,0,0/
!     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
!     ERROR TABLE INITIALLY
     data kount(1),kount(2),kount(3),kount(4),kount(5),&
    &     kount(6),kount(7),kount(8),kount(9),kount(10)&
     &     /0,0,0,0,0,0,0,0,0,0/
      data kountx/0/
!     NEXT DATA STATEMENT SETS UP OUTPUT FORMAT
     data f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),&
    &     f(11),f(12),f(13),f(14),f(15),f(16),f(17)&
    &     /'(' ,'1' ,'X' ,',' ,'A' ,' ' ,' ' ,',' ,'I' ,' ' ,&
     &      ' ' ,',' ,'2' ,'I' ,'1' ,'0' ,')' /
      if (nmessg.gt.0) go to 80
!     DUMP THE TABLE
         if (kount(1).eq.0) return
!        PREPARE FORMAT
         nchar = i1mach(6)
         call s88fmt(2,nchar,f(6))
         ncol = 20 - nchar
         call s88fmt(2,ncol,f(10))
!        PRINT TO EACH UNIT
         call xgetua(lun,nunit)
         do 60 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
!           PRINT TABLE HEADER
            write (iunit,10)
  10       format ('0          ERROR MESSAGE SUMMARY'/&
     &              ' FIRST WORD      NERR     LEVEL     COUNT')
!           PRINT BODY OF TABLE
            do 20 i=1,10
               if (kount(i).eq.0) go to 30
               write (iunit,f) mestab(i),nertab(i),levtab(i),kount(i)
   20       continue
   30       continue
!           PRINT NUMBER OF OTHER ERRORS
            if (kountx.ne.0) write (iunit,40) kountx
   40       format (/' OTHER ERRORS NOT INDIVIDUALLY TABULATED=',i10)
            write (iunit,50)
   50       format (1x)
   60    continue
         if (nmessg.lt.0) return
!        CLEAR THE ERROR TABLES
         do 70 i=1,10
   70       kount(i) = 0
         kountx = 0
         return
   80 continue
!     PROCESS A MESSAGE...
!     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      do 90 i=1,10
         ii = i
         if (kount(i).eq.0) go to 110
         if (messg(1).ne.mestab(i)) go to 90
         if (nerr.ne.nertab(i)) go to 90
         if (level.ne.levtab(i)) go to 90
         go to 100
   90 continue
!     THREE POSSIBLE CASES...
!     TABLE IS FULL
         kountx = kountx+1
         icount = 1
         return
!     MESSAGE FOUND IN TABLE
  100    kount(ii) = kount(ii) + 1
         icount = kount(ii)
         return
!     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    mestab(ii) = messg(1)
         nertab(ii) = nerr
         levtab(ii) = level
         kount(ii)  = 1
         icount = 1
         return
      end
!DGAMIT
      double precision function dgamit (a, x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE TRICOMI-S INCOMPLETE GAMMA FUNCTION DEFINED BY
!
! GAMIT = X**(-A)/GAMMA(A) * INTEGRAL T = 0 TO X OF EXP(-T) * T**(A-1.)
!
! AND ANALYTIC CONTINUATION FOR A .LE. 0.0.  GAMMA(X) IS THE COMPLETE
! GAMMA FUNCTION OF X.  GAMIT IS EVALUATED FOR ARBITRARY REAL VALUES OF
! A AND FOR NON-NEGATIVE VALUES OF X (EVEN THOUGH GAMIT IS DEFINED FOR
! X .LT. 0.0), EXCEPT THAT FOR X = 0 AND A .LE. 0.0, GAMIT IS INFINITE,
! A FATAL ERROR.
!
!      A SLIGHT DETERIORATION OF 2 OR 3 DIGITS ACCURACY WILL OCCUR WHEN
! GAMIT IS VERY LARGE OR VERY SMALL IN ABSOLUTE VALUE, BECAUSE LOG-
! ARITHMIC VARIABLES ARE USED.  ALSO, IF THE PARAMETER A IS VERY CLOSE
! TO A NEGATIVE INTEGER (BUT NOT A NEGATIVE INTEGER), THERE IS A LOSS
! OF ACCURACY, WHICH IS REPORTED IF THE RESULT IS LESS THAN HALF
! MACHINE PRECISION.
!
! REF. -- W. GAUTSCHI, AN EVALUATION PROCEDURE FOR INCOMPLETE GAMMA
! FUNCTIONS, ACM TRANS. MATH. SOFTWARE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,x
!
!  LOCAL SCALARS
     double precision aeps,ainta,algap1,alneps,alng,alx,bot,h,sga,&
     &   sgngam,sqeps,t
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,D9GMIT,D9LGIC,D9LGIT,DGAMR,DLNGAM
!       EXTERNAL D1MACH,D9GMIT,D9LGIC,D9LGIT,DGAMR,DLNGAM
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL DLGAMS,XERCLR,XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,dsign,dsqrt,exp,int,log
!
!
      data alneps, sqeps, bot / 3*0.d0 /
!
      if (alneps.ne.0.d0) go to 10
      alneps = -log (d1mach(3))
      sqeps = dsqrt (d1mach(4))
      bot = log (d1mach(1))
!
 10   if (x.lt.0.d0) call xerror ('DGAMIT  X IS NEGATIVE', 21, 2, 2)
!
      if (x.ne.0.d0) alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = dsign (1.0d0, a)
      ainta = int (a + 0.5d0*sga)
      aeps = a - ainta
!
      if (x.gt.0.d0) go to 20
      dgamit = 0.0d0
      if (ainta.gt.0.d0 .or. aeps.ne.0.d0) dgamit = dgamr(a+1.0d0)
      return
!
 20   if (x.gt.1.d0) go to 30
     if (a.ge.(-0.5d0) .or. aeps.ne.0.d0) call dlgams (a+1.0d0, algap1,&
     &  sgngam)
      dgamit = d9gmit (a, x, algap1, sgngam, alx)
      return
!
 30   if (a.lt.x) go to 40
      t = d9lgit (a, x, dlngam(a+1.0d0))
      if (t.lt.bot) call xerclr
      dgamit = exp (t)
      return
!
 40   alng = d9lgic (a, x, alx)
!
! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
!
      h = 1.0d0
      if (aeps.eq.0.d0 .and. ainta.le.0.d0) go to 50
!
      call dlgams (a+1.0d0, algap1, sgngam)
      t = log (abs(a)) + alng - algap1
      if (t.gt.alneps) go to 60
!
      if (t.gt.(-alneps)) h = 1.0d0 - sga * sgngam * exp(t)
      if (abs(h).gt.sqeps) go to 50
!
      call xerclr
      call xerror ('DGAMIT  RESULT LT HALF PRECISION', 32, 1, 1)
!
 50   t = -a*alx + log(abs(h))
      if (t.lt.bot) call xerclr
      dgamit = dsign (exp(t), h)
      return
!
 60   t = t - a*alx
      if (t.lt.bot) call xerclr
      dgamit = -sga * sgngam * exp(t)
      return
!
      end
!GAMIT
      real function gamit (a, x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE TRICOMI-S INCOMPLETE GAMMA FUNCTION DEFINED BY
!
! GAMIT = X**(-A)/GAMMA(A) * INTEGRAL T = 0 TO X OF EXP(-T) * T**(A-1.)
!
! AND ANALYTIC CONTINUATION FOR A .LE. 0.0.  GAMMA(X) IS THE COMPLETE
! GAMMA FUNCTION OF X.  GAMIT IS EVALUATED FOR ARBITRARY REAL VALUES OF
! A AND FOR NON-NEGATIVE VALUES OF X (EVEN THOUGH GAMIT IS DEFINED FOR
! X .LT. 0.0), EXCEPT THAT FOR X = 0 AND A .LE. 0.0, GAMIT IS INFINITE,
! A FATAL ERROR.
!
!      A SLIGHT DETERIORATION OF 2 OR 3 DIGITS ACCURACY WILL OCCUR WHEN
! GAMIT IS VERY LARGE OR VERY SMALL IN ABSOLUTE VALUE, BECAUSE LOG-
! ARITHMIC VARIABLES ARE USED.  ALSO, IF THE PARAMETER A IS VERY CLOSE
! TO A NEGATIVE INTEGER (BUT NOT A NEGATIVE INTEGER), THERE IS A LOSS
! OF ACCURACY, WHICH IS REPORTED IF THE RESULT IS LESS THAN HALF
! MACHINE PRECISION.
!
! REF. -- W. GAUTSCHI, AN EVALUATION PROCEDURE FOR INCOMPLETE GAMMA
! FUNCTIONS, ACM TRANS. MATH. SOFTWARE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real a,x
!
!  LOCAL SCALARS
      real aeps,ainta,algap1,alneps,alng,alx,bot,h,sga,sgngam,sqeps,t
!
!  EXTERNAL FUNCTIONS
!      REAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
!       EXTERNAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL ALGAMS,XERCLR,XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,aint,exp,log,sign,sqrt
!
      data alneps, sqeps, bot / 3*0.0 /
!
      if (alneps.ne.0.0) go to 10
      alneps = -log(r1mach(3))
      sqeps = sqrt(r1mach(4))
      bot = log(r1mach(1))
!
 10   if (x.lt.0.0) call xerror ('GAMIT   X IS NEGATIVE', 21, 2, 2)
!
      if (x.ne.0.0) alx = log(x)
      sga = 1.0
      if (a.ne.0.0) sga = sign (1.0, a)
      ainta = aint (a+0.5*sga)
      aeps = a - ainta
!
      if (x.gt.0.0) go to 20
      gamit = 0.0
      if (ainta.gt.0.0 .or. aeps.ne.0.0) gamit = gamr(a+1.0)
      return
!
 20   if (x.gt.1.0) go to 40
     if (a.ge.(-0.5) .or. aeps.ne.0.0) call algams (a+1.0, algap1,&
     &  sgngam)
      gamit = r9gmit (a, x, algap1, sgngam, alx)
      return
!
 40   if (a.lt.x) go to 50
      t = r9lgit (a, x, alngam(a+1.0))
      if (t.lt.bot) call xerclr
      gamit = exp(t)
      return
!
 50   alng = r9lgic (a, x, alx)
!
! EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
!
      h = 1.0
      if (aeps.eq.0.0 .and. ainta.le.0.0) go to 60
      call algams (a+1.0, algap1, sgngam)
      t = log(abs(a)) + alng - algap1
      if (t.gt.alneps) go to 70
      if (t.gt.(-alneps)) h = 1.0 - sga*sgngam*exp(t)
      if (abs(h).gt.sqeps) go to 60
      call xerclr
      call xerror ('GAMIT   RESULT LT HALF PRECISION', 32, 1, 1)
!
 60   t = -a*alx + log(abs(h))
      if (t.lt.bot) call xerclr
      gamit = sign (exp(t), h)
      return
!
 70   t = t - a*alx
      if (t.lt.bot) call xerclr
      gamit = -sga*sgngam*exp(t)
      return
!
      end
!DLGAMS
      subroutine dlgams (x, dlgam, sgngam)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE LOG ABS (GAMMA(X)) AND RETURN THE SIGN OF GAMMA(X) IN SGNGAM.
! SGNGAM IS EITHER +1.0 OR -1.0.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision dlgam,sgngam,x
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION DLNGAM
!       EXTERNAL DLNGAM
!
!  INTRINSIC FUNCTIONS
      intrinsic int,mod
!
!
      dlgam = dlngam(x)
      sgngam = 1.0d0
      if (x.gt.0.d0) return
!
!     INT = DMOD (-INT(X), 2.0D0) + 0.1D0
      if (int(mod(-int(x),2)+0.1d0).eq.0) sgngam = -1.0d0
!
      return
      end
!GAMR
      real function gamr (x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! THIS ROUTINE, NOT GAMMA(X), SHOULD BE THE FUNDAMENTAL ONE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real x
!
!  LOCAL SCALARS
      real alngx,sgngx
      integer irold
!
!  EXTERNAL FUNCTIONS
!      REAL GAMMA
!       EXTERNAL GAMMA
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL ALGAMS,XERCLR,XGETF,XSETF
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,aint,exp
!
      gamr = 0.0
      if (x.le.0.0 .and. aint(x).eq.x) return
!
      call xgetf (irold)
      call xsetf (1)
      if (abs(x).gt.10.0) go to 10
      gamr = 1.0/gamma(x)
      call xerclr
      call xsetf (irold)
      return
!
 10   call algams (x, alngx, sgngx)
      call xerclr
      call xsetf (irold)
      gamr = sgngx * exp(-alngx)
      return
!
      end
!XERCLR
      subroutine xerclr
!
!     ABSTRACT
!        THIS ROUTINE SIMPLY RESETS THE CURRENT ERROR NUMBER TO ZERO.
!        THIS MAY BE NECESSARY TO DO IN ORDER TO DETERMINE THAT
!        A CERTAIN ERROR HAS OCCURRED AGAIN SINCE THE LAST TIME
!        NUMXER WAS REFERENCED.
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  7 JUNE 1978
!
!
!  VARIABLE DECLARATIONS
!
!  LOCAL SCALARS
      integer junk
!
!  EXTERNAL FUNCTIONS
!      INTEGER J4SAVE
!       EXTERNAL J4SAVE
!
      junk = j4save(1,0,.true.)
      return
      end
!BETAI
      real function betai (x, pin, qin)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! BASED ON BOSTEN AND BATTISTE, REMARK ON ALGORITHM 179, COMM. ACM,
! V 17, P 153, (1974).
!
! X   VALUE TO WHICH FUNCTION IS TO BE INTEGRATED. X MUST BE IN (0,1).
! P   INPUT (1ST) PARAMETER (MUST BE GREATER THAN 0)
! Q   INPUT (2ND) PARAMETER (MUST BE GREATER THAN 0)
! BETAI  INCOMPLETE BETA FUNCTION RATIO, THE PROBABILITY THAT A RANDOM
!        VARIABLE FROM A BETA DISTRIBUTION HAVING PARAMETERS P AND Q
!        WILL BE LESS THAN OR EQUAL TO X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real pin,qin,x
!
!  LOCAL SCALARS
      real alneps,alnsml,c,eps,fac1,fac2,finsum,p,p1,ps,q,sml,term,xb,y
      integer i,ib,n
!
!  EXTERNAL FUNCTIONS
!      REAL ALBETA,R1MACH
!       EXTERNAL ALBETA,R1MACH
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,aint,exp,float,log,max,min,real
!
      data             eps, alneps, sml, alnsml / 4*0.0 /
!
      if (eps.ne.0.) go to 10
      eps = r1mach(3)
      alneps = log(eps)
      sml = r1mach(1)
      alnsml = log(sml)
!
10   if (x.lt.0. .or. x.gt.1.0) call xerror (&
     &  'BETAI   X IS NOT IN THE RANGE (0,1)', 35, 1, 2)
     if (pin.le.0. .or. qin.le.0.) call xerror (&
     &  'BETAI   P AND/OR Q IS LE ZERO', 29, 2, 2)
!
      y = x
      p = pin
      q = qin
      if (q.le.p .and. x.lt.0.8) go to 20
      if (x.lt.0.2) go to 20
      y = 1.0 - y
      p = qin
      q = pin
!
 20   if ((p+q)*y/(p+1.).lt.eps) go to 80
!
! EVALUATE THE INFINITE SUM FIRST.
! TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
!
      ps = q - aint(q)
      if (ps.eq.0.) ps = 1.0
      xb = p*log(y) -  albeta(ps, p) - log(p)
      betai = 0.0
      if (xb.ge.alnsml) then
         betai = exp(xb)
         fac2 = 1.0
         if (ps.ne.1.0e0) then
            fac1 = 1.0
            n = max(alneps/log(y), 4.0e0)
            do 30 i=1,n
               if ((i-ps.eq.0.0e0) .or. (fac1.eq.0.0e0)) then
                  fac1 = 0.0e0
               else
                 if (log(abs(fac1)) + log(abs(i-ps)) + log(y) -&
     &                log(real(i)) .lt. alnsml) then
                     fac1 = 0.0e0
                  else
                     fac1 = fac1 * (i-ps)*y/i
                  end if
               end if
               fac2 = fac2 + fac1*p/(p+i)
 30         continue
         end if
         betai = betai*fac2
      end if
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
      if (q.le.1.0) go to 70
!
      xb = p*log(y) + q*log(1.0-y) - albeta(p,q) - log(q)
      ib = max (xb/alnsml, 0.0)
      term = exp (xb - float(ib)*alnsml)
      c = 1.0/(1.0-y)
      p1 = q*c/(p+q-1.)
!
      finsum = 0.0
      n = q
      if (q.eq.float(n)) n = n - 1
      do 50 i=1,n
        if (p1.le.1.0 .and. term/eps.le.finsum) go to 60
        if (q-i+1.0e0 .eq. 0.0e0) then
          term = 0.0e0
        else
         if (log(abs(q-i+1.0e0)) + log(abs(c)) + log(abs(term)) -&
     &        log(abs(p+q-i)) .lt. alnsml) then
            term = 0.0e0
          else
            term = (q-i+1.0e0)*c*term/(p+q-i)
          end if
        end if
!
        if (term.gt.1.0) ib = ib - 1
        if (term.gt.1.0) term = term*sml
!
        if (ib.eq.0) finsum = finsum + term
 50   continue
!
 60   betai = betai + finsum
 70   if (y.ne.x .or. p.ne.pin) betai = 1.0 - betai
      betai = max (min (betai, 1.0), 0.0)
      return
!
 80   betai = 0.0
      xb = p*log(max(y,sml)) - log(p) - albeta(p,q)
      if (xb.gt.alnsml .and. y.ne.0.) betai = exp (xb)
      if (y.ne.x .or. p.ne.pin) betai = 1.0 - betai
      return
!
      end
!J4SAVE
      integer function j4save(iwhich,ivalue,iset)
!
!     ABSTRACT
!        J4SAVE SAVES AND RECALLS SEVERAL GLOBAL VARIABLES NEEDED
!        BY THE LIBRARY ERROR HANDLING ROUTINES.
!
!     DESCRIPTION OF PARAMETERS
!      --INPUT--
!        IWHICH - INDEX OF ITEM DESIRED.
!                 = 1 REFERS TO CURRENT ERROR NUMBER.
!                 = 2 REFERS TO CURRENT ERROR CONTROL FLAG.
!                 = 3 REFERS TO CURRENT UNIT NUMBER TO WHICH ERROR
!                     MESSAGES ARE TO BE SENT.  (0 MEANS USE STANDARD.)
!                 = 4 REFERS TO THE MAXIMUM NUMBER OF TIMES ANY
!                     MESSAGE IS TO BE PRINTED (AS SET BY XERMAX).
!                 = 5 REFERS TO THE TOTAL NUMBER OF UNITS TO WHICH
!                     EACH ERROR MESSAGE IS TO BE WRITTEN.
!                 = 6 REFERS TO THE 2ND UNIT FOR ERROR MESSAGES
!                 = 7 REFERS TO THE 3RD UNIT FOR ERROR MESSAGES
!                 = 8 REFERS TO THE 4TH UNIT FOR ERROR MESSAGES
!                 = 9 REFERS TO THE 5TH UNIT FOR ERROR MESSAGES
!        IVALUE - THE VALUE TO BE SET FOR THE IWHICH-TH PARAMETER,
!                 IF ISET IS .TRUE. .
!        ISET   - IF ISET=.TRUE., THE IWHICH-TH PARAMETER WILL BE
!                 GIVEN THE VALUE, IVALUE.  IF ISET=.FALSE., THE
!                 IWHICH-TH PARAMETER WILL BE UNCHANGED, AND IVALUE
!                 IS A DUMMY PARAMETER.
!      --OUTPUT--
!        THE (OLD) VALUE OF THE IWHICH-TH PARAMETER WILL BE RETURNED
!        IN THE FUNCTION VALUE, J4SAVE.
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     ADAPTED FROM BELL LABORATORIES PORT LIBRARY ERROR HANDLER
!     LATEST REVISION ---  23 MAY 1979
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer ivalue,iwhich
      logical iset
!
!  LOCAL ARRAYS
      integer iparam(9)
!
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,1,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
!XGETF
      subroutine xgetf(kontrl)
!
!     ABSTRACT
!        XGETF RETURNS THE CURRENT VALUE OF THE ERROR CONTROL FLAG
!        IN KONTRL.  SEE SUBROUTINE XSETF FOR FLAG VALUE MEANINGS.
!        (KONTRL IS AN OUTPUT PARAMETER ONLY.)
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  7 JUNE 1978
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer kontrl
!
!  EXTERNAL FUNCTIONS
!      INTEGER J4SAVE
!       EXTERNAL J4SAVE
!
      kontrl = j4save(2,0,.false.)
      return
      end
!D9LGMC
      double precision function d9lgmc (x)
! AUGUST 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG GAMMA CORRECTION FACTOR FOR X .GE. 10. SO THAT
! LOG (DGAMMA(X)) = LOG(DSQRT(2*PI)) + (X-.5)*LOG(X) - X + D9LGMC(X)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision x
!
!  LOCAL SCALARS
      double precision xbig,xmax
      integer nalgm
!
!  LOCAL ARRAYS
      double precision algmcs(15)
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DCSEVL
!      INTEGER INITDS
!       EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic dsqrt,exp,log,min,sngl
!
!
! SERIES FOR ALGM       ON THE INTERVAL  0.          TO  1.00000E-02
!                                        WITH WEIGHTED ERROR   1.28E-31
!                                         LOG WEIGHTED ERROR  30.89
!                               SIGNIFICANT FIGURES REQUIRED  29.81
!                                    DECIMAL PLACES REQUIRED  31.48
!
      data algmcs(  1) / +.1666389480451863247205729650822d+0      /
      data algmcs(  2) / -.1384948176067563840732986059135d-4      /
      data algmcs(  3) / +.9810825646924729426157171547487d-8      /
      data algmcs(  4) / -.1809129475572494194263306266719d-10     /
      data algmcs(  5) / +.6221098041892605227126015543416d-13     /
      data algmcs(  6) / -.3399615005417721944303330599666d-15     /
      data algmcs(  7) / +.2683181998482698748957538846666d-17     /
      data algmcs(  8) / -.2868042435334643284144622399999d-19     /
      data algmcs(  9) / +.3962837061046434803679306666666d-21     /
      data algmcs( 10) / -.6831888753985766870111999999999d-23     /
      data algmcs( 11) / +.1429227355942498147573333333333d-24     /
      data algmcs( 12) / -.3547598158101070547199999999999d-26     /
      data algmcs( 13) / +.1025680058010470912000000000000d-27     /
      data algmcs( 14) / -.3401102254316748799999999999999d-29     /
      data algmcs( 15) / +.1276642195630062933333333333333d-30     /
!
      data nalgm, xbig, xmax / 0, 2*0.d0 /
!
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, 15, sngl(d1mach(3)) )
      xbig = 1.0d0/dsqrt(d1mach(3))
      xmax = exp (min(log(d1mach(2)/12.d0), -log(12.d0*d1mach(1))))
!
 10   if (x.lt.10.d0) call xerror ('D9LGMC  X MUST BE GE 10', 23, 1, 2)
      if (x.ge.xmax) go to 20
!
      d9lgmc = 1.d0/(12.d0*x)
     if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs,&
     &  nalgm) / x
      return
!
 20   d9lgmc = 0.d0
      call xerror ('D9LGMC  X SO BIG D9LGMC UNDERFLOWS', 34, 2, 1)
      return
!
      end
!ALGAMS
      subroutine algams (x, algam, sgngam)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE LOG ABS (GAMMA(X)) AND RETURN THE SIGN OF GAMMA(X) IN SGNGAM.
! SGNGAM IS EITHER +1.0 OR -1.0.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      real algam,sgngam,x
!
!  EXTERNAL FUNCTIONS
!      REAL ALNGAM
!       EXTERNAL ALNGAM
!
!  INTRINSIC FUNCTIONS
      intrinsic int,mod
!
      algam = alngam(x)
      sgngam = 1.0
      if (x.gt.0.0) return
!
!     INT = AMOD (-AINT(X), 2.0) + 0.1
      if (int(mod(-int(x),2)+0.1).eq.0) sgngam = -1.0
!
      return
      end
!DGAMI
      double precision function dgami (a, x)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE THE INCOMPLETE GAMMA FUNCTION DEFINED BY
!
! GAMI = INTEGRAL FROM T = 0 TO X OF EXP(-T) * T**(A-1.0) .
!
! GAMI IS EVALUATED FOR POSITIVE VALUES OF A AND NON-NEGATIVE VALUES
! OF X.  A SLIGHT DETERIORATION OF 2 OR 3 DIGITS ACCURACY WILL OCCUR
! WHEN GAMI IS VERY LARGE OR VERY SMALL, BECAUSE LOGARITHMIC VARIABLES
! ARE USED.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      double precision a,x
!
!  LOCAL SCALARS
      double precision factor
!
!  EXTERNAL FUNCTIONS
!      DOUBLE PRECISION D1MACH,DGAMIT,DLNGAM
!       EXTERNAL D1MACH,DGAMIT,DLNGAM
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      intrinsic exp,log
!
!
      if (a.le.0.d0) call xerror ('DGAMI   A MUST BE GT ZERO', 25, 1,2)
      if (x.lt.0.d0) call xerror ('DGAMI   X MUST BE GE ZERO', 25, 2,2)
!
      dgami = 0.d0
      if (x.eq.0.0d0) return
!
! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
!
      factor = dlngam(a) + a*log(x)
      if (factor.gt.log(d1mach(2))) then
         dgami = d1mach(2)
      else
         dgami = exp(factor) * dgamit(a,x)
      end if
!
      return
      end

!I1MACH
      integer function i1mach(i)
      use,intrinsic :: iso_fortran_env, only : stdin=>input_unit
      use,intrinsic :: iso_fortran_env, only : stdout=>output_unit
      use,intrinsic :: iso_fortran_env, only : stderr=>error_unit
!
!     MODIFIED JANUARY 24, 1990 TO ACCORD WITH CMLIB AND PORT VERSIONS
!
!  I/O UNIT NUMBERS.
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS.
!
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!
!  INTEGERS.
!
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!
!    I1MACH( 7) = A, THE BASE.
!
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!
!    I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.   FOR FORTRAN 77, YOU MAY WISH
!  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
!  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
!
!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
!  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
!  FOR IMACH(1) - IMACH(4).
!
! intrinsics
! exponent     -  Exponent function
! fraction     -  Fractional part of the model representation
! nearest      -  Nearest representable number
! rrspacing    -  Reciprocal of the relative spacing
! scale        -  Scale a real value
! set_exponent -  Set the exponent of the model
! spacing      -  Smallest distance between two numbers of a given type
! digits       -  Significant digits function
! epsilon      -  Epsilon function
! maxexponent  -  Maximum exponent of a real kind
! minexponent  -  Minimum exponent of a real kind
! precision    -  Decimal precision of a real kind
! radix        -  Base of a model number
! range        -  Decimal exponent range of a real kind
! tiny         -  Smallest positive number of a real kind
! huge         -  Largest number of a kind
! bit_size     -  Bit size inquiry function
! storage_size - Storage size in bits
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer i
!
!  LOCAL SCALARS
      integer output,sanity
!
!  LOCAL ARRAYS
      integer imach(16)
!
!  EXTERNAL SUBROUTINES
!      external fdump
!
!  EQUIVALENCES
      equivalence (imach(4),output)
!
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -125 /
!      DATA IMACH(13) /  128 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA IMACH( 1) /    7 /
!     DATA IMACH( 2) /    2 /
!     DATA IMACH( 3) /    2 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -256 /
!     DATA IMACH(13) /  255 /
!     DATA IMACH(14) /   60 /
!     DATA IMACH(15) / -256 /
!     DATA IMACH(16) /  255 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  48 /
!     DATA IMACH( 6) /   6 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /   8 /
!     DATA IMACH(11) /  13 /
!     DATA IMACH(12) / -50 /
!     DATA IMACH(13) /  76 /
!     DATA IMACH(14) /  26 /
!     DATA IMACH(15) / -50 /
!     DATA IMACH(16) /  76 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  48 /
!     DATA IMACH( 6) /   6 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /   8 /
!     DATA IMACH(11) /  13 /
!     DATA IMACH(12) / -50 /
!     DATA IMACH(13) /  76 /
!     DATA IMACH(14) /  26 /
!     DATA IMACH(15) / -32754 /
!     DATA IMACH(16) /  32780 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    7 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   60 /
!     DATA IMACH( 6) /   10 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   48 /
!     DATA IMACH( 9) / O"00007777777777777777" /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   48 /
!     DATA IMACH(12) / -974 /
!     DATA IMACH(13) / 1070 /
!     DATA IMACH(14) /   96 /
!     DATA IMACH(15) / -927 /
!     DATA IMACH(16) / 1070 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS/VE
!
!     DATA IMACH( 1) /     5 /
!     DATA IMACH( 2) /     6 /
!     DATA IMACH( 3) /     7 /
!     DATA IMACH( 4) /     6 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -4095 /
!     DATA IMACH(13) /  4094 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -4095 /
!     DATA IMACH(16) /  4094 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA IMACH( 1) /      5 /
!     DATA IMACH( 2) /      6 /
!     DATA IMACH( 3) /      7 /
!     DATA IMACH( 4) /      6 /
!     DATA IMACH( 5) /     64 /
!     DATA IMACH( 6) /      8 /
!     DATA IMACH( 7) /      2 /
!     DATA IMACH( 8) /     47 /
!     DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
!     DATA IMACH(10) /      2 /
!     DATA IMACH(11) /     47 /
!     DATA IMACH(12) / -28625 /
!     DATA IMACH(13) /  28718 /
!     DATA IMACH(14) /     94 /
!     DATA IMACH(15) / -28625 /
!     DATA IMACH(16) /  28718 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN5 COMPILER)
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /   60 /
!     DATA IMACH( 6) /   10 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   48 /
!     DATA IMACH( 9) / O"00007777777777777777" /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   48 /
!     DATA IMACH(12) / -974 /
!     DATA IMACH(13) / 1070 /
!     DATA IMACH(14) /   96 /
!     DATA IMACH(15) / -927 /
!     DATA IMACH(16) / 1070 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN COMPILER)
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    7 /
!     DATA IMACH( 4) /6LOUTPUT/
!     DATA IMACH( 5) /   60 /
!     DATA IMACH( 6) /   10 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   47 /
!     DATA IMACH(12) / -929 /
!     DATA IMACH(13) / 1070 /
!     DATA IMACH(14) /   94 /
!     DATA IMACH(15) / -929 /
!     DATA IMACH(16) / 1069 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR CONVEX C-1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP2 AND XMP3
!
!     DATA IMACH( 1) /     5 /
!     DATA IMACH( 2) /     6 /
!     DATA IMACH( 3) /   102 /
!     DATA IMACH( 4) /     6 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    63 /
!     DATA IMACH( 9) /  777777777777777777777B /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -8189 /
!     DATA IMACH(13) /  8190 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -8099 /
!     DATA IMACH(16) /  8190 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     DATA IMACH( 1) /   11 /
!     DATA IMACH( 2) /   12 /
!     DATA IMACH( 3) /    8 /
!     DATA IMACH( 4) /   10 /
!     DATA IMACH( 5) /   16 /
!     DATA IMACH( 6) /    2 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   15 /
!     DATA IMACH( 9) /32767 /
!     DATA IMACH(10) /   16 /
!     DATA IMACH(11) /    6 /
!     DATA IMACH(12) /  -64 /
!     DATA IMACH(13) /   63 /
!     DATA IMACH(14) /   14 /
!     DATA IMACH(15) /  -64 /
!     DATA IMACH(16) /   63 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA IMACH( 1) /       5 /
!     DATA IMACH( 2) /       6 /
!     DATA IMACH( 3) /       0 /
!     DATA IMACH( 4) /       6 /
!     DATA IMACH( 5) /      24 /
!     DATA IMACH( 6) /       3 /
!     DATA IMACH( 7) /       2 /
!     DATA IMACH( 8) /      23 /
!     DATA IMACH( 9) / 8388607 /
!     DATA IMACH(10) /       2 /
!     DATA IMACH(11) /      23 /
!     DATA IMACH(12) /    -127 /
!     DATA IMACH(13) /     127 /
!     DATA IMACH(14) /      38 /
!     DATA IMACH(15) /    -127 /
!     DATA IMACH(16) /     127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70
!                           THE HONEYWELL 600/6000 SERIES
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /   43 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   63 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH(1) /      5/
!     DATA IMACH(2) /      6 /
!     DATA IMACH(3) /      4 /
!     DATA IMACH(4) /      1 /
!     DATA IMACH(5) /     16 /
!     DATA IMACH(6) /      2 /
!     DATA IMACH(7) /      2 /
!     DATA IMACH(8) /     15 /
!     DATA IMACH(9) /  32767 /
!     DATA IMACH(10)/      2 /
!     DATA IMACH(11)/     23 /
!     DATA IMACH(12)/   -128 /
!     DATA IMACH(13)/    127 /
!     DATA IMACH(14)/     39 /
!     DATA IMACH(15)/   -128 /
!     DATA IMACH(16)/    127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH(1) /      5 /
!     DATA IMACH(2) /      6 /
!     DATA IMACH(3) /      4 /
!     DATA IMACH(4) /      1 /
!     DATA IMACH(5) /     16 /
!     DATA IMACH(6) /      2 /
!     DATA IMACH(7) /      2 /
!     DATA IMACH(8) /     15 /
!     DATA IMACH(9) /  32767 /
!     DATA IMACH(10)/      2 /
!     DATA IMACH(11)/     23 /
!     DATA IMACH(12)/   -128 /
!     DATA IMACH(13)/    127 /
!     DATA IMACH(14)/     55 /
!     DATA IMACH(15)/   -128 /
!     DATA IMACH(16)/    127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES
!                           THE XEROX SIGMA 5/7/9
!                           THE SEL SYSTEMS 85/86
!                           THE PERKIN ELMER 3230
!                           THE PERKIN ELMER (INTERDATA) 7/32
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  32 /
!     DATA IMACH( 6) /   4 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  31 /
!     DATA IMACH( 9) / Z7FFFFFFF /
!     DATA IMACH(10) /  16 /
!     DATA IMACH(11) /   6 /
!     DATA IMACH(12) / -64 /
!     DATA IMACH(13) /  63 /
!     DATA IMACH(14) /  14 /
!     DATA IMACH(15) / -64 /
!     DATA IMACH(16) /  63 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32 WITH THE UNIX SYSTEM
!     FORTRAN 77 COMPILER
!
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE THE Z'S
!     SPECIFYING HEX CONSTANTS WITH Y'S
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   6 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  32 /
!     DATA IMACH( 6) /   4 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  31 /
!     DATA IMACH( 9) / Z'7FFFFFFF' /
!     DATA IMACH(10) /  16 /
!     DATA IMACH(11) /   6 /
!     DATA IMACH(12) / -64 /
!     DATA IMACH(13) /  62 /
!     DATA IMACH(14) /  14 /
!     DATA IMACH(15) / -64 /
!     DATA IMACH(16) /  62 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    5 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   54 /
!     DATA IMACH(15) / -101 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    5 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   62 /
!     DATA IMACH(15) / -128 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     32-BIT INTEGER ARITHMETIC
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   32 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   56 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     16-BIT INTEGER ARITHMETIC
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   16 /
!     DATA IMACH( 6) /    2 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   15 /
!     DATA IMACH( 9) / 32767 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   56 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!     SUPPLIED BY IGOR BRAY.
!
!      DATA IMACH( 1) /            1 /
!      DATA IMACH( 2) /            1 /
!      DATA IMACH( 3) /            2 /
!      DATA IMACH( 4) /            1 /
!      DATA IMACH( 5) /           32 /
!      DATA IMACH( 6) /            4 /
!      DATA IMACH( 7) /            2 /
!      DATA IMACH( 8) /           31 /
!      DATA IMACH( 9) / :17777777777 /
!      DATA IMACH(10) /            2 /
!      DATA IMACH(11) /           23 /
!      DATA IMACH(12) /         -127 /
!      DATA IMACH(13) /         +127 /
!      DATA IMACH(14) /           47 /
!      DATA IMACH(15) /       -32895 /
!      DATA IMACH(16) /       +32637 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!
!      DATA IMACH( 1) /     0 /
!      DATA IMACH( 2) /     0 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     0 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE SUN-3/160
!     (SEE ALSO IEEE CONSTANTS ABOVE)
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    6 /
!     DATA IMACH( 4) /    0 /
!     DATA IMACH( 5) /   32 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -125 /
!     DATA IMACH(13) /  128 /
!     DATA IMACH(14) /   53 /
!     DATA IMACH(15) / -1021 /
!     DATA IMACH(16) /  1024 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    7 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    6 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   60 /
!     DATA IMACH(15) /-1024 /
!     DATA IMACH(16) / 1023 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE VAX 11/780 WITH FORTRAN IV-PLUS COMPILER
!                   AND FOR THE VAX/VMS VERSION 2.2 WITHOUT G_FLOATING
!
!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   32 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   56 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /, SANITY/987/
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA IMACH( 1) /     1/
!     DATA IMACH( 2) /     1/
!     DATA IMACH( 3) /     0/
!     DATA IMACH( 4) /     1/
!     DATA IMACH( 5) /    16/
!     DATA IMACH( 6) /     2/
!     DATA IMACH( 7) /     2/
!     DATA IMACH( 8) /    15/
!     DATA IMACH( 9) / 32767/
!     DATA IMACH(10) /     2/
!     DATA IMACH(11) /    24/
!     DATA IMACH(12) /  -127/
!     DATA IMACH(13) /   127/
!     DATA IMACH(14) /    56/
!     DATA IMACH(15) /  -127/
!     DATA IMACH(16) /   127/, SANITY/987/
      integer,parameter :: bpi=bit_size(0)
      integer,parameter :: cpi=storage_size(0)/storage_size('a')
      integer,parameter :: ia=radix(0)
      integer,parameter :: is=digits(0)
      integer,parameter :: ibig=huge(0)
      integer,parameter :: rb=radix(0.0e0)
      integer,parameter :: rt=digits(0.0e0)
      integer,parameter :: rmine=minexponent(0.0e0)
      integer,parameter :: rmaxe=maxexponent(0.0e0)
      integer,parameter :: dt=digits(0.0d0)
      integer,parameter :: dmine=minexponent(0.0d0)
      integer,parameter :: dmaxe=maxexponent(0.0d0)
      data imach( 1)/stdin/                               ! I/O UNIT NUMBERS : THE STANDARD INPUT UNIT.
      data imach( 2)/stdout/                              ! I/O UNIT NUMBERS : THE STANDARD OUTPUT UNIT.
      data imach( 3)/7/                                   ! I/O UNIT NUMBERS : THE STANDARD PUNCH UNIT.
      data imach( 4)/stderr/                              ! I/O UNIT NUMBERS : THE STANDARD ERROR MESSAGE UNIT.
      data imach( 5)/bpi/                                 ! WORDS            : THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
      data imach( 6)/cpi/                                 ! WORDS            : THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!  INTEGERS.
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
      data imach( 7)/ia/                ! A, THE BASE.
      data imach( 8)/is/                ! S, THE NUMBER OF BASE-A DIGITS.
      data imach( 9)/ibig/              ! A**S - 1, THE LARGEST MAGNITUDE.
!  FLOATING-POINT NUMBERS.
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT, BASE-B FORM
!         SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!         WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!         0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
      data imach(10)/rb/                 ! B, THE BASE.
!  SINGLE-PRECISION
      data imach(11)/rt/                 ! T, THE NUMBER OF BASE-B DIGITS.
      data imach(12)/rmine/              ! EMIN, THE SMALLEST EXPONENT E.
      data imach(13)/rmaxe/              ! EMAX, THE LARGEST EXPONENT E.
!  DOUBLE-PRECISION
      data imach(14)/dt/                 ! T, THE NUMBER OF BASE-B DIGITS.
      data imach(15)/dmine/              ! EMIN, THE SMALLEST EXPONENT E.
      data imach(16)/dmaxe/, sanity/987/ ! EMAX, THE LARGEST EXPONENT E.
!
!  ***  ISSUE STOP IF ALL DATA STATEMENTS ARE COMMENTED...
      if (sanity .ne. 987) then
         stop 'I1MACH, D1MACH AND R1MACH HAVE NOT BEEN INITIALIZED'
      else
!
         if (i .lt. 1  .or.  i .gt. 16) then
            write(output,9000)
!            call fdump()
            stop
         else
            i1mach=imach(i)
         end if
      end if
!
      return
!
 9000 format('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
!
      end function i1mach
!R1MACH
      real function r1mach(i)
!
!     MODIFIED JANUARY 22, 1990 TO ACCORD WITH CMLIB AND PORT VERSIONS
!
!  SINGLE-PRECISION MACHINE CONSTANTS
!
!  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!
!  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!
!  R1MACH(5) = LOG10(B)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.
!
!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
!  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
!
!  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
!  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
!  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer i
!
!  LOCAL ARRAYS
      real rmach(5)
      integer diver(2),large(2),log10(2),right(2),small(2)
!
!  EXTERNAL FUNCTIONS
!      integer i1mach
!      external i1mach
!
!  EXTERNAL SUBROUTINES
!      external seterr
!
!  EQUIVALENCES
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
!      DATA SMALL(1) /     8388608 /
!      DATA LARGE(1) /  2139095039 /
!      DATA RIGHT(1) /   864026624 /
!      DATA DIVER(1) /   872415232 /
!      DATA LOG10(1) /  1050288283 /
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!      DATA SMALL(1) /    1048576 /
!      DATA LARGE(1) / 2147483647 /
!      DATA RIGHT(1) /  990904320 /
!      DATA DIVER(1) / 1007681536 /
!      DATA LOG10(1) / 1091781651 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA RMACH(1) / Z400800000 /
!     DATA RMACH(2) / Z5FFFFFFFF /
!     DATA RMACH(3) / Z4E9800000 /
!     DATA RMACH(4) / Z4EA800000 /
!     DATA RMACH(5) / Z500E730E8 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
!
!     DATA RMACH(1) / O1771000000000000 /
!     DATA RMACH(2) / O0777777777777777 /
!     DATA RMACH(3) / O1311000000000000 /
!     DATA RMACH(4) / O1301000000000000 /
!     DATA RMACH(5) / O1157163034761675 /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS
!
!     DATA RMACH(1) / O"00014000000000000000" /
!     DATA RMACH(2) / O"37767777777777777777" /
!     DATA RMACH(3) / O"16404000000000000000" /
!     DATA RMACH(4) / O"16414000000000000000" /
!     DATA RMACH(5) / O"17164642023241175720" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS/VE
!
!     DATA RMACH(1) / Z"3001800000000000" /
!     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
!     DATA RMACH(3) / Z"3FD2800000000000" /
!     DATA RMACH(4) / Z"3FD3800000000000" /
!     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA RMACH(1) / X'9000400000000000' /
!     DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /
!     DATA RMACH(3) / X'FFA3400000000000' /
!     DATA RMACH(4) / X'FFA4400000000000' /
!     DATA RMACH(5) / X'FFD04D104D427DE8' /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN5 COMPILER)
!
!     DATA RMACH(1) / O"00014000000000000000" /
!     DATA RMACH(2) / O"37767777777777777777" /
!     DATA RMACH(3) / O"16404000000000000000" /
!     DATA RMACH(4) / O"16414000000000000000" /
!     DATA RMACH(5) / O"17164642023241175720" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN COMPILER)
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR CONVEX C-1.
!
!      DATA RMACH(1) / '00800000'X /
!      DATA RMACH(2) / '7FFFFFFF'X /
!      DATA RMACH(3) / '34800000'X /
!      DATA RMACH(4) / '35000000'X /
!      DATA RMACH(5) / '3F9A209B'X /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP2 AND XMP3
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC RMACH(5)
!
!     DATA SMALL/20K,0/,LARGE/77777K,177777K/
!     DATA RIGHT/35420K,0/,DIVER/36020K,0/
!     DATA LOG10/40423K,42023K/
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
!     DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
!     DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
!     DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES
!                           THE HONEYWELL 600/6000 SERIES
!
!     DATA RMACH(1) / O402400000000 /
!     DATA RMACH(2) / O376777777777 /
!     DATA RMACH(3) / O714400000000 /
!     DATA RMACH(4) / O716400000000 /
!     DATA RMACH(5) / O776464202324 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE91), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES
!                           THE XEROX SIGMA 5/7/9
!                           THE SEL SYSTEMS 85/86
!                           THE PERKIN ELMER 3230
!                           THE PERKIN ELMER (INTERDATA) 3230
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32 WITH THE UNIX SYSTEM
!     FORTRAN 77 COMPILER
!
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE THE Z'S
!     SPECIFYING HEX CONSTANTS WITH Y'S
!
!     DATA RMACH(1) / Z'00100000' /
!     DATA RMACH(2) / Z'7EFFFFFF' /
!     DATA RMACH(3) / Z'3B100000' /
!     DATA RMACH(4) / Z'3C100000' /
!     DATA RMACH(5) / Z'41134413' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL)
!
!     DATA SMALL(1) /    8388608 /
!     DATA LARGE(1) / 2147483647 /
!     DATA RIGHT(1) /  880803840 /
!     DATA DIVER(1) /  889192448 /
!     DATA LOG10(1) / 1067065499 /
!
!     DATA RMACH(1) / O00040000000 /
!     DATA RMACH(2) / O17777777777 /
!     DATA RMACH(3) / O06440000000 /
!     DATA RMACH(4) / O06500000000 /
!     DATA RMACH(5) / O07746420233 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL)
!
!     DATA SMALL(1),SMALL(2) /   128,     0 /
!     DATA LARGE(1),LARGE(2) / 32767,    -1 /
!     DATA RIGHT(1),RIGHT(2) / 13440,     0 /
!     DATA DIVER(1),DIVER(2) / 13568,     0 /
!     DATA LOG10(1),LOG10(2) / 16282,  8347 /
!
!     DATA SMALL(1),SMALL(2) / O000200, O000000 /
!     DATA LARGE(1),LARGE(2) / O077777, O177777 /
!     DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
!     DATA DIVER(1),DIVER(2) / O032400, O000000 /
!     DATA LOG10(1),LOG10(2) / O037632, O020233 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!
!      DATA SMALL(1) / $00800000 /
!      DATA LARGE(1) / $7F7FFFFF /
!      DATA RIGHT(1) / $33800000 /
!      DATA DIVER(1) / $34000000 /
!      DATA LOG10(1) / $3E9A209B /
!
!     MACHINE CONSTANTS FOR THE SUN-3/160
!     (SEE ALSO IEEE CONSTANTS ABOVE)
!
!     DATA SMALL(1) / X'00800000' /
!     DATA LARGE(1) / X'7F7FFFFF' /
!     DATA RIGHT(1) / X'33800000' /
!     DATA DIVER(1) / X'34000000' /
!     DATA LOG10(1) / X'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     DATA RMACH(1) / O000400000000 /
!     DATA RMACH(2) / O377777777777 /
!     DATA RMACH(3) / O146400000000 /
!     DATA RMACH(4) / O147400000000 /
!     DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE VAX 11/780
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
!     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR THE VAX/VMS VERSION 2.2
!
!     DATA RMACH(1) / '00000080'X /
!     DATA RMACH(2) / 'FFFF7FFF'X /
!     DATA RMACH(3) / '00003480'X /
!     DATA RMACH(4) / '00003500'X /
!     DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA SMALL(1),SMALL(2) /     0,    256/
!     DATA LARGE(1),LARGE(2) /    -1,   -129/
!     DATA RIGHT(1),RIGHT(2) /     0,  26880/
!     DATA DIVER(1),DIVER(2) /     0,  27136/
!     DATA LOG10(1),LOG10(2) /  8347,  32538/
!
      integer,parameter       :: b = radix(0.0e0)
      integer,parameter       :: t = digits (0.0e0 )
      real,parameter :: emin=minexponent(0.0e0)
      real,parameter :: emax=maxexponent(0.0e0)

      select case(i)
      case(1); r1mach = tiny(0.0e0)                ! B**(EMIN-1), the smallest positive magnitude.
      case(2); r1mach = huge(0.0e0)                ! B**EMAX*(1-B**(-T)), the largest magnitude.
                                                   ! calculating this by formula could cause overflow without using a larger type
      case(3); r1mach = real(b)**(-t)              ! B**(-T), the smallest relative spacing.
      case(4); r1mach = epsilon(0.0e0)             ! B**(1-T), the largest relative spacing.
      case(5);
      block
      intrinsic log10
      r1mach = log10(real(b))             ! log10(B).
      endblock
      case default
         call seterr('R1MACH - I OUT OF BOUNDS',24,1,2)
      end select
!
      end function r1mach
!D1MACH
      double precision function d1mach(i)
!
!     MODIFIED JANUARY 24, 1990 TO ACCORD WITH CMLIB AND PORT VERSIONS
!
!  DOUBLE-PRECISION MACHINE CONSTANTS
!
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!
!  D1MACH( 5) = LOG10(B)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.
!
!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
!  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
!
!  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
!  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
!  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer i
!
!  LOCAL ARRAYS
      double precision dmach(5)
      integer diver(4),large(4),log10(4),right(4),small(4)
!
!  EXTERNAL FUNCTIONS
!
!  EXTERNAL SUBROUTINES
!      external seterr
!
!  EQUIVALENCES
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
!
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
!
!      DATA SMALL(1),SMALL(2) /    1048576,          0 /
!      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
!      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
!      DATA DIVER(1),DIVER(2) / 1018167296,          0 /
!      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
!     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
!     SIGNIFICANT BYTE IS STORED FIRST.
!
!      DATA SMALL(1),SMALL(2) /          0,    1048576 /
!      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
!      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
!      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
!      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!      DATA SMALL(1),SMALL(2) /    1048576,          0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
!      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
!      DATA DIVER(1),DIVER(2) /  873463808,          0 /
!      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA SMALL(1) / ZC00800000 /
!     DATA SMALL(2) / Z000000000 /
!
!     DATA LARGE(1) / ZDFFFFFFFF /
!     DATA LARGE(2) / ZFFFFFFFFF /
!
!     DATA RIGHT(1) / ZCC5800000 /
!     DATA RIGHT(2) / Z000000000 /
!
!     DATA DIVER(1) / ZCC6800000 /
!     DATA DIVER(2) / Z000000000 /
!
!     DATA LOG10(1) / ZD00E730E7 /
!     DATA LOG10(2) / ZC77800DC0 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O0000000000000000 /
!
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O0007777777777777 /
!
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O7770000000000000 /
!
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O7777777777777777 /
!
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS
!
!     DATA SMALL(1) / O"00604000000000000000" /
!     DATA SMALL(2) / O"00000000000000000000" /
!
!     DATA LARGE(1) / O"37767777777777777777" /
!     DATA LARGE(2) / O"37167777777777777777" /
!
!     DATA RIGHT(1) / O"15604000000000000000" /
!     DATA RIGHT(2) / O"15000000000000000000" /
!
!     DATA DIVER(1) / O"15614000000000000000" /
!     DATA DIVER(2) / O"15010000000000000000" /
!
!     DATA LOG10(1) / O"17164642023241175717" /
!     DATA LOG10(2) / O"16367571421742254654" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170/180 SERIES USING NOS/VE
!
!     DATA SMALL(1) / Z"3001800000000000" /
!     DATA SMALL(2) / Z"3001000000000000" /
!
!     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
!     DATA LARGE(2) / Z"4FFE000000000000" /
!
!     DATA RIGHT(1) / Z"3FD2800000000000" /
!     DATA RIGHT(2) / Z"3FD2000000000000" /
!
!     DATA DIVER(1) / Z"3FD3800000000000" /
!     DATA DIVER(2) / Z"3FD3000000000000" /
!
!     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
!     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA SMALL(1) / X'9000400000000000' /
!     DATA SMALL(2) / X'8FD1000000000000' /
!
!     DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
!     DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
!
!     DATA RIGHT(1) / X'FF74400000000000' /
!     DATA RIGHT(2) / X'FF45000000000000' /
!
!     DATA DIVER(1) / X'FF75400000000000' /
!     DATA DIVER(2) / X'FF46000000000000' /
!
!     DATA LOG10(1) / X'FFD04D104D427DE7' /
!     DATA LOG10(2) / X'FFA17DE623E2566A' /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN5 COMPILER)
!
!     DATA SMALL(1) / O"00604000000000000000" /
!     DATA SMALL(2) / O"00000000000000000000" /
!
!     DATA LARGE(1) / O"37767777777777777777" /
!     DATA LARGE(2) / O"37167777777777777777" /
!
!     DATA RIGHT(1) / O"15604000000000000000" /
!     DATA RIGHT(2) / O"15000000000000000000" /
!
!     DATA DIVER(1) / O"15614000000000000000" /
!     DATA DIVER(2) / O"15010000000000000000" /
!
!     DATA LOG10(1) / O"17164642023241175717" /
!     DATA LOG10(2) / O"16367571421742254654" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES (FTN COMPILER)
!
!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /
!
!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /
!
!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /
!
!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /
!
!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /
!
!     MACHINE CONSTANTS FOR CONVEX C-1
!
!      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
!      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
!      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
!      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
!      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP2 AND XMP3
!
!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /
!
!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777776B /
!
!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /
!
!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /
!
!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(5)
!
!     DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
!     DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
!     DATA LOG10/40423K,42023K,50237K,74776K/
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
!     DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
!     DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
!     DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES
!                           THE HONEYWELL 600/6000 SERIES
!
!     DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!     DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!     DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!     DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!     DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!     DATA SMALL(3), SMALL(4) /       0,       1 /
!     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!     DATA DIVER(3), DIVER(4) /       0,    227B /
!     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES
!                           THE XEROX SIGMA 5/7/9
!                           THE SEL SYSTEMS 85/86
!                           THE PERKIN-ELMER 3230
!                           THE PERKIN-ELMER (INTERDATA) 7/32
!
!     DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
!     DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!     DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
!     DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
!     DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
!
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32 WITH THE UNIX SYSTEM
!     FORTRAN 77 COMPILER
!
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE THE Z'S
!     SPECIFYING HEX CONSTANTS WITH Y'S
!
!     DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
!     DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
!     DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL)
!
!     DATA SMALL(1),SMALL(2) /    8388608,           0 /
!     DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!     DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!     DATA DIVER(1),DIVER(2) /  620756992,           0 /
!     DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
!
!     DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
!     DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
!     DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
!     DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
!     DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL)
!
!     DATA SMALL(1),SMALL(2) /    128,      0 /
!     DATA SMALL(3),SMALL(4) /      0,      0 /
!
!     DATA LARGE(1),LARGE(2) /  32767,     -1 /
!     DATA LARGE(3),LARGE(4) /     -1,     -1 /
!
!     DATA RIGHT(1),RIGHT(2) /   9344,      0 /
!     DATA RIGHT(3),RIGHT(4) /      0,      0 /
!
!     DATA DIVER(1),DIVER(2) /   9472,      0 /
!     DATA DIVER(3),DIVER(4) /      0,      0 /
!
!     DATA LOG10(1),LOG10(2) /  16282,   8346 /
!     DATA LOG10(3),LOG10(4) / -31493, -12296 /
!
!     DATA SMALL(1),SMALL(2) / O000200, O000000 /
!     DATA SMALL(3),SMALL(4) / O000000, O000000 /
!
!     DATA LARGE(1),LARGE(2) / O077777, O177777 /
!     DATA LARGE(3),LARGE(4) / O177777, O177777 /
!
!     DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
!     DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
!
!     DATA DIVER(1),DIVER(2) / O022400, O000000 /
!     DATA DIVER(3),DIVER(4) / O000000, O000000 /
!
!     DATA LOG10(1),LOG10(2) / O037632, O020232 /
!     DATA LOG10(3),LOG10(4) / O102373, O147770 /
!
!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!     SUPPLIED BY IGOR BRAY.
!
!      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
!      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
!      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
!      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
!      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
!
!      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
!      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
!      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
!      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
!      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
!
!     MACHINE CONSTANTS FOR THE SUN-3/160
!     (SEE ALSO IEEE CONSTANTS ABOVE)
!
!     DATA SMALL(1),SMALL(2) / X'00100000', X'00000000' /
!     DATA LARGE(1),LARGE(2) / X'7FEFFFFF', X'FFFFFFFF' /
!     DATA RIGHT(1),RIGHT(2) / X'3CA00000', X'00000000' /
!     DATA DIVER(1),DIVER(2) / X'3CB00000', X'00000000' /
!     DATA LOG10(1),LOG10(2) / X'3FD34413', X'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES (FTN COMPILER)
!
!     DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!     DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!     DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!     DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!     DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
!
!     MACHINE CONSTANTS FOR VAX 11/780
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!
!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!
!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!
!     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2 COMPILER
!
!     DATA SMALL(1),SMALL(2) / '00000080'X, '00000000'X /
!     DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
!     DATA RIGHT(1),RIGHT(2) / '00002480'X, '00000000'X /
!     DATA DIVER(1),DIVER(2) / '00002500'X, '00000000'X /
!     DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /
!
      integer,parameter       :: b = radix(0.0d0)
      integer,parameter       :: t = digits (0.0d0)
      double precision,parameter :: emin=minexponent(0.0d0)
      double precision,parameter :: emax=maxexponent(0.0d0)

      select case(i)
      case(1); d1mach = tiny(0.0d0)                        ! B**(EMIN-1), the smallest positive magnitude.
      case(2); d1mach = huge(0.0d0)                        ! B**EMAX*(1-B**(-T)), the largest magnitude.
                                                           ! calculating this by formula could cause overflow without using a larger type
      case(3); d1mach = real(b,kind=kind(0.0d0))**(-t)     ! B**(-T), the smallest relative spacing.
      case(4); d1mach = epsilon(0.0d0)                     ! B**(1-T), the largest relative spacing.
      case(5)
      block
      intrinsic log10
      d1mach = log10(real(b,kind=kind(0.0d0)))    ! log10(B).
      endblock
      case default
         call seterr('D1MACH - I OUT OF BOUNDS',24,1,2)
      end select
!
      end function d1mach

!SETIV
      subroutine setiv(vector, n, value)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS THE FIRST N ELEMENTS OF AN INTEGER VECTOR
!
!     WRITTEN BY  -  JOHN E. KOONTZ
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!        ADAPTED FROM SETRV, WRITTEN BY LINDA L. MITCHELL
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   n,value
!
!  ARRAY ARGUMENTS
     integer&
     &   vector(n)
!
!  LOCAL SCALARS
     integer&
     &   i
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        *
!     INTEGER N
!        NUMBER OF ELEMENTS TO SET
!     INTEGER VALUE
!        VALUE TO WHICH THE ELEMENTS ARE TO BE SET
!     INTEGER VECTOR(N)
!        VECTOR WHOSE FIRST N ELEMENTS ARE TO BE SET.
!
      do 10 i=1,n
         vector(i) = value
   10 continue
!
      return
!
      end
!FIXPRT
      subroutine fixprt(ifix, fixed)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS THE CHARACTER ARRAY FIXED.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  APRIL 2, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ifix
!
!  ARRAY ARGUMENTS
     character&
     &   fixed(3)*1
!
!  LOCAL SCALARS
     integer&
     &   i
!
!  LOCAL ARRAYS
     character&
     &   no(3)*1,yes(3)*1
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     CHARACTER*1 FIXED(3)
!        THE CHARACTERS USED TO LABEL THE PARAMETERS FIXED OR NOT.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IFIX
!        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
!        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
!        IF IFIX.EQ.0, THEN FIXED WILL BE SET TO NO.
!        IF IFIX.NE.0, THEN FIXED WILL BE SET TO YES.
!     CHARACTER*1 NO(3)
!        THE CHARACTERS BLANK, N, AND O
!     CHARACTER*1 YES(3)
!        THE CHARACTERS Y, E, AND S
!
      data no(1)/' '/, no(2)/'N'/, no(3)/'O'/
      data yes(1)/'Y'/, yes(2)/'E'/, yes(3)/'S'/
!
      if (ifix.ne.0) then
!
!     SET FIXED TO YES
!
         do 10 i = 1, 3
            fixed(i) = yes(i)
   10    continue
!
      else
!
!     SET FIXED TO NO
!
         do 20 i = 1, 3
            fixed(i) = no(i)
   20    continue
      end if
!
      return
!
      end
!BACKOP
      subroutine backop (mspec, nfac, npardf, mbol, mbo, nparma, nparar)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     COMPUTE NUMBER OF BACK ORDER TERMS FOR ARIMA MODEL
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 4, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   mbo,mbol,nfac,nparar,npardf,nparma
!
!  ARRAY ARGUMENTS
     integer&
     &   mspec(4,*)
!
!  LOCAL SCALARS
     integer&
     &   j
!
!  INTRINSIC FUNCTIONS
      intrinsic max
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER J
!        AN INDEX VARIABLE.
!     INTEGER MBO
!        THE MAXIMUM BACK ORDER OPERATOR.
!     INTEGER MBOL
!        THE MAXIMUM BACK ORDER ON THE LEFT
!     INTEGER MSPEC(4,NFAC)
!        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
!     INTEGER NFAC
!        THE NUMBER OF FACTORS IN THE MODEL
!     INTEGER NPARAR
!        THE NUMBER OF AUTOREGRESSIVE PARAMETERS
!     INTEGER NPARDF
!        THE ORDER OF THE EXPANDED DIFFERENCE FILTER.
!     INTEGER NPARMA
!        THE LENGTH OF THE VECTOR PARMA
!
!     COMPUTE DEGREE OF BACK OPERATOR RESULTING FROM THE NDF
!     DIFFERENCING FACTORS (= ND DOT IOD).
!
      nparar = 0
      npardf = 0
      nparma = 0
      if (nfac .eq. 0) go to 20
      do 10 j = 1, nfac
         nparar = nparar + mspec(1,j)*mspec(4,j)
         npardf = npardf + mspec(2,j)*mspec(4,j)
         nparma = nparma + mspec(3,j)*mspec(4,j)
   10 continue
!
   20 continue
!
      mbol = npardf + nparar
      mbo = max(mbol,nparma)
!
      return
!
      end
!IPRINT
      subroutine iprint(iprt)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS THE LOGICAL UNIT FOR PRINTED OUTPUT.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iprt
!
!  EXTERNAL FUNCTIONS
!      INTEGER
!     +   I1MACH
!      EXTERNAL I1MACH
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR OUTPUT.
!
      iprt = i1mach(2)
      return
      end
!FFTLEN
      subroutine fftlen(n, ndiv, nfft)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE COMPUTES THE SMALLEST VALUE OF NFFT WHICH
!     EQUALS OR EXCEEDS N + 2, SUCH THAT NFFT - 2 IS DIVISIBLE BY
!     NDIV AND HAS NO PRIME FACTORS GREATER THAN 23, AND THE
!     PRODUCT OF THE NON SQUARE PRIME FACTORS OF NFFT - 2 DO NOT
!     EXCEED 209.  THE VALUE OF NFFT THUS MEET THE REQUIREMENTS OF
!     THE EXTENDED LENGTH OF THE SERIES REQUIRED FOR ANY ROUTINE
!     USING THE SINGLETON FFT PROVIDING THE PROPER VALUE OF NDIV
!     IS CHOSEN.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   n,ndiv,nfft
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   iprt
     logical&
     &   err01,err02,head
!
!  LOCAL ARRAYS
     character&
     &   ln(8)*1,lndiv(8)*1,nmsub(6)*1
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EISGE,IPRINT,SETESL
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERR01, ERR02
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
!        IF IERR .EQ. 1, ERRORS WERE DETECTED.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     CHARACTER*1 LN(8), LNDIV(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
!        CHECKED FOR ERRORS.
!     INTEGER N
!        THE NUMBER UPON WHICH NFFT IS BASED.
!     INTEGER NDIV
!        A REQUIRED FACTOR OF NFFT - 2.
!     INTEGER NFFT
!        THE RETURNED VALUE WHICH MEETS THE ABOVE DESCRIPTION.
!     CHARACTER*1 NMSUB(6)
!        THE NAME OF THE ROUTINE CALLING THE ERROR CHECKING
!        SUBROUTINE.
!
!     SET UP NAME ARRAYS
!
     data&
    &  nmsub(1),  nmsub(2),  nmsub(3),  nmsub(4),  nmsub(5),  nmsub(6)&
     & /     'F',       'F',       'T',       'L',       'E',       'N'/
     data&
    &     ln(1),     ln(2),     ln(3),     ln(4),     ln(5),     ln(6)&
     & /     'N',       ' ',       ' ',       ' ',       ' ',       ' '/
     data&
    &     ln(7),     ln(8)&
     & /     ' ',       ' '/
     data&
    &  lndiv(1),  lndiv(2),  lndiv(3),  lndiv(4),  lndiv(5),  lndiv(6)&
     & /     'N',       'D',       'I',       'V',       ' ',       ' '/
     data&
    &  lndiv(7),  lndiv(8)&
     & /     ' ',       ' '/
!
!     ERROR CHECKING
!
      ierr = 0
      head = .true.
!
      call eisge(nmsub, ln, n, 1, 1, head, err01, ln)
!
      call eisge(nmsub, lndiv, ndiv, 1, 1, head, err02, lndiv)
!
!
      if ((.not. err01) .and. (.not. err02)) go to 10
!
!     PRINT PROPER CALL SEQUENCE
!
      ierr = 1
      call iprint (iprt)
      write (iprt, 1000)
!
      return
!
   10 continue
!
      call setesl(n, ndiv, nfft)
!
      return
!
!     FORMAT STATEMENTS
!
1000 format (/42h the correct form of the call statement is//&
     &   34h       call fftlen (n, ndiv, nfft))
!
      end
!ECVF
      subroutine ecvf(nmsub)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS AN ERROR MESSAGE WHEN THE LAG VALUE OF
!     THE LAST COVARIANCE COMPUTED BEFORE ONE WAS NOT COMPUTED
!     DUE TO MISSING DATA DOES NOT EXCEED ZERO.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JUNE 10, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1
!
!  LOCAL SCALARS
     integer&
     &   iprt
     logical&
     &   head
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IPRT
!        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
!     CHARACTER*1 NMSUB(6)
!        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
!
      call iprint(iprt)
!
      head = .true.
!
      call ehdr(nmsub, head)
!
      write(iprt, 1010)
      return
!
!     FORMAT STATEMENTS
!
1010 format (/&
    &   46h the covariances at lags zero and/or one could,&
    &   16h not be computed/&
    &   49h because of missing data.  no further analysis is,&
     &   10h possible.)
!
      end
!AMEHDR
      subroutine amehdr(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE NONLINEAR
!     LEAST SQUARES ESTIMATION ROUTINES FOR ARIMA MODELS THAT USE
!     NUMERICAL APPROXIMATIONS TO THE DERIVATIVES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  AUGUST 1, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format ('+NONLINEAR LEAST SQUARES ESTIMATION',&
     &   ' FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED')
1010 format ('+', 77(1h*)/&
    &   1x, 37h*  nonlinear least squares estimation,&
    &   40h for the parameters of an arima model  */&
    &   2h *, 16x, 45h             using backforecasts             ,&
     &   14x, 1h*/1x, 77(1h*))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!LLHDRP
      subroutine llhdrp(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE POLYNOMIAL LINEAR
!     LEAST SQUARES LLSTING ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
      if (page) write (iprt,1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt,1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (32h+linear least squares estimation,&
     &   33h with polynomial model, continued)
1010 format ('+', 59('*')/&
    &   1x, 34h*  linear least squares estimation,&
    &   25h with polynomial model  */ 1x,&
     &   59('*'))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!EIAGEP
     subroutine eiagep (nmsub, nmvar, ymmn, nvmx, head, msgtyp, nv,&
     &   nmmin)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS THE ERROR MESSAGES FOR ERAGT AND ERAGTM.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JUNE 10, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   msgtyp,nv,nvmx,ymmn
     logical&
     &   head
!
!  ARRAY ARGUMENTS
     character&
     &   nmmin(8)*1,nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER MSGTYP
!        THE INDICATOR ARGUMENT FOR THE TYPE OF MESSAGE.
!        IF (MSGTYP.GE.3) THE MESSAGE PRINTED WILL USE NMMIN
!        OTHERWISE IT WILL USE YMMN.
!        IF (MSGTYP = 1 OR 3) NO VIOLATIONS ARE ALLOWED.
!        IF (MSGTYP = 2 OR 4) THE NUMBER OF VIOLATIONS MUST
!                             BE LESS THAN   NVMX   .
!     CHARACTER*1 NMMIN(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING ROUTINES NAME.
!     CHARACTER*1 NMVAR(8)
!        THE CHARACTERS OF THE PARAMETERS NAME.
!     INTEGER NV
!        THE NUMBER OF VIOLATIONS FOUND.
!     INTEGER NVMX
!        THE LARGEST NUMBER OF VIOLATIONS ALLOWED.
!     INTEGER YMMN
!        THE MINIMUM ACCEPTABLE VALUE.
!
      call iprint(iprt)
      call ehdr(nmsub, head)
!
     if (msgtyp.le.2)&
     &   write (iprt, 1000) (nmvar(i),i=1,6), ymmn, nv
     if (msgtyp.ge.3)&
     &   write (iprt, 1005) (nmvar(i),i=1,6), (nmmin(i),i=1,8), nv
!
      go to (10, 20, 30, 40), msgtyp
!
   10 write(iprt, 1010) (nmvar(i),i=1,6), ymmn
      return
!
   20 write(iprt, 1020) (nmvar(i),i=1,6), ymmn, nvmx
      return
!
   30 write(iprt, 1030) (nmvar(i),i=1,6), (nmmin(i),i=1,8)
      return
!
   40 write(iprt, 1040) (nmvar(i),i=1,6), (nmmin(i),i=1,8), nvmx
      return
!
!     FORMAT STATEMENTS
!
1000 format (/&
    &   31h the number of values in array , 6a1,&
     &   ' LESS THAN ', i5, 4h is , i6, '.')
1005 format (/&
    &   31h the number of values in array , 6a1,&
     &   ' LESS THAN ', 8a1, 4h is , i6, '.')
1010 format(&
    &   25h the values in the array , 6a1,&
     &   ' MUST ALL BE GREATER THAN OR EQUAL TO ', i5, '.')
1020 format(&
    &   35h the number of values in the array , 6a1,&
    &   ' LESS THAN ', 8a1/&
     &   19h must be less than , i5, '.')
1030 format(&
    &   25h the values in the array , 6a1,&
     &   ' MUST ALL BE GREATER THAN OR EQUAL TO ', i5, '.')
1040 format(&
    &   35h the number of values in the array , 6a1,&
    &   ' LESS THAN ', 8a1/&
     &   19h must be less than , i5, '.')
!
      end
!EIVEQ
     subroutine eiveq (nmsub, nmvar1, ivec, n, ival, neqmn, head, neq,&
     &                  nne, msgtyp, error, nmvar2, nmvar3)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE NUMBER OF ELEMENTS OF IVEC EQUAL
!     TO IVAL IS AT LEAST NEQMN.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  AUGUST 3, 1987
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ival,msgtyp,n,neq,neqmn,nne
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     integer&
     &   ivec(*)
     character&
     &   nmsub(6)*1,nmvar1(8)*1,nmvar2(8)*1,nmvar3(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER IVEC(N)
!        THE VECTOR BEING CHECKED.
!     INTEGER MSGTYP
!        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF ERROR IS TRUE AND
!        MSGTYP = 1, THE INPUT VALUE WAS TOO SMALL BASED ON LIMITS
!                    IMPOSED BY STARPAC.
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR1(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED.
!     CHARACTER*1 NMVAR2(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED
!        AGAINST.
!     CHARACTER*1 NMVAR3(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT THAT THE ELEMENTS
!        MUST BE EQUAL TO.
!     INTEGER NEQ
!        THE NUMBER OF ELEMENTS EQUAL TO IVAL.
!     INTEGER NEQMN
!        THE MINIMUM NUMBER OF ELEMENTS EQUAL TO IVAL WHICH IS OK.
!     INTEGER NNE
!        THE NUMBER OF ELEMENTS NOT EQUAL TO IVAL.
!
      error = .false.
!
      if (n.le.0) return
!
!     CHECK FOR VALUES EQUAL TO IVAL
!
      neq = 0
      do 10 i = 1, n
         if (ivec(i) .eq. ival) neq = neq + 1
   10 continue
!
      nne = n - neq
      if (neq .ge. neqmn) return
!
!     INSUFFICIENT NUMBER OF ELEMENTS EQUAL TO IVAL.
!
      error = .true.
!
      call iprint(iprt)
!
      call ehdr(nmsub, head)
!
     if (msgtyp.eq.1) write(iprt, 1000)&
    &   (nmvar1(i),i=1,8), (nmvar2(i),i=1,8), neq,&
     &   (nmvar2(i),i=1,8), (nmvar3(i),i=1,8)
!
      return
!
!     FORMAT STATEMENTS
!
1000 format(&
    &   ' THE NUMBER OF ELEMENTS IN ', 8a1,&
    &   ' EQUAL TO ', 8a1, ' IS ', i6, '.'/&
    &   ' THE NUMBER OF ELEMENTS EQUAL TO ', 8a1,&
     &   ' MUST BE GREATER THAN OR EQUAL TO ', 8a1, '.')
!
      end
!ACFER
     subroutine acfer(nmsub, n, lagmax, lacov, ldstak, ldsmin,&
     &  differ, nfac, nd, iod, isfft, lyfft, nfft)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS IS THE ERROR CHECKING ROUTINE FOR THE ACF FAMILY
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 2, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   lacov,lagmax,ldsmin,ldstak,lyfft,n,nfac,nfft
     logical&
     &   differ,isfft
!
!  ARRAY ARGUMENTS
     integer&
     &   iod(*),nd(*)
     character&
     &   nmsub(6)*1
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i
     logical&
     &   head
!
!  LOCAL ARRAYS
     logical&
     &   err(15)
     character&
    &   llacov(8)*1,llagmx(8)*1,llds(8)*1,llgmx1(8)*1,&
    &   llyfft(8)*1,ln(8)*1,lnfft(8)*1,lnm1(8)*1,lone(8)*1,&
     &   lthree(8)*1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EISGE,EISII,ERDF
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL DIFFER
!        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
!        ROUTINE IS ACFD (DIFFER = TRUE) OR NOT (DIFFER = FALSE)
!     LOGICAL ERR(15)
!        VALUES INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
!     INTEGER IOD(NFAC)
!        THE ORDER OF EACH OF THE DIFFERENCE VACTORS
!     LOGICAL ISFFT
!        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
!        ROUTINE HAS SUFFIX F (ISFFT = TRUE) OR NOT (ISFFT = FALSE)
!     INTEGER LACOV
!        THE LENGTH OF THE VECTOR ACOV.
!     INTEGER LAGMAX
!        THE MAXIMUM LAG VALUE REQUESTED.
!     INTEGER LDSMIN
!        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
!     INTEGER LDSTAK
!        THE LENGTH OF THE ARRAY DSTAK.
!     CHARACTER*1 LLACOV(8), LLAGMX(8), LLDS(8), LLGMX1(8), LLYFFT(8),
!    *  LN(8), LNFFT(8), LNM1(8), LONE(8), LTHREE(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
!        CHECKED FOR ERRORS.
!     INTEGER LYFFT
!        THE LENGTH OF THE VECTOR YFFT.
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
!     INTEGER ND(NFAC)
!        THE ARRAY CONTAINING THE NUMBER OF TIMES THE DIFFERENCE FACTORS
!        ARE TO BE APPLIED
!     INTEGER NFAC
!        THE NUMBER OF FACTORS.
!     INTEGER NFFT
!        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
!        SUBROUTINE.
!
!
!     SET UP NAME ARRAYS
!
     data&
    & llacov(1), llacov(2), llacov(3), llacov(4), llacov(5),&
     & llacov(6), llacov(7), llacov(8) /'L','A','C','O','V',' ',' ',' '/
     data&
    & llagmx(1), llagmx(2), llagmx(3), llagmx(4), llagmx(5),&
     & llagmx(6), llagmx(7), llagmx(8) /'L','A','G','M','A','X',' ',' '/
     data&
    & llds(1), llds(2), llds(3), llds(4), llds(5),&
     & llds(6), llds(7), llds(8) /'L','D','S','T','A','K',' ',' '/
     data&
    & llgmx1(1), llgmx1(2), llgmx1(3), llgmx1(4), llgmx1(5),&
     & llgmx1(6), llgmx1(7), llgmx1(8) /'L','A','G','M','A','X','+','1'/
     data&
    & llyfft(1), llyfft(2), llyfft(3), llyfft(4), llyfft(5),&
     & llyfft(6), llyfft(7), llyfft(8) /'L','Y','F','F','T',' ',' ',' '/
     data&
    & ln(1), ln(2), ln(3), ln(4), ln(5),&
     & ln(6), ln(7), ln(8) /'N',' ',' ',' ',' ',' ',' ',' '/
     data&
    & lnm1(1), lnm1(2), lnm1(3), lnm1(4), lnm1(5),&
     & lnm1(6), lnm1(7), lnm1(8) /'(','N','-','1',')',' ',' ',' '/
     data&
    & lnfft(1), lnfft(2), lnfft(3), lnfft(4), lnfft(5),&
     & lnfft(6), lnfft(7), lnfft(8) /'N','F','F','T',' ',' ',' ',' '/
     data&
    & lone(1), lone(2), lone(3), lone(4), lone(5),&
     & lone(6), lone(7), lone(8) /'O','N','E',' ',' ',' ',' ',' '/
     data&
    & lthree(1), lthree(2), lthree(3), lthree(4), lthree(5),&
     & lthree(6), lthree(7), lthree(8) /'T','H','R','E','E',' ',' ',' '/
!
!     SET UP FOR ERROR CHECKING
!
      ierr = 0
      head = .true.
      do 10 i = 1, 15
        err(i) = .false.
   10 continue
!
!     CALL ERROR CHECKING ROUTINES
!
      call eisge(nmsub, ln, n, 3, 2, head, err(1), lthree)
!
      if (.not.err(1)) then
!
       call eisii(nmsub, llagmx, lagmax, 1, n-1, 1, head, err(2), lone,&
     &    lnm1)
!
        if (differ) call erdf(nmsub, nfac, nd, iod, n, head, err(3))
!
        if (.not.err(2)) then
!
         call eisge(nmsub, llacov, lacov, lagmax+1, 2, head, err(4),&
     &      llgmx1)
!
          call eisge(nmsub, llds, ldstak, ldsmin, 9, head, err(5), llds)
!
         if (isfft)&
    &      call eisge(nmsub, llyfft, lyfft, nfft, 2, head, err(6),&
     &      lnfft)
        end if
      end if
!
      do 20 i = 1, 15
        if (err(i)) ierr = 1
   20 continue
!
      return
!
      end
!NCHOSE
      integer function nchose(n,k)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE IS USED TO COMBINE THE DIFFERENCE FACTORS FROM A
!     (BOX-JENKINS) TIME SERIES MODEL.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   k,n
!
!  LOCAL SCALARS
     integer&
     &   i,kk,nn
!
!  INTRINSIC FUNCTIONS
      intrinsic min
!
!
      if (n .gt. k) go to 10
      nchose = 1
      return
!
   10 kk = min(k, n - k)
      nn = 1
      do 20 i = 1, kk
         nn = (nn*(n - i + 1))/i
   20 continue
      nchose = nn
      return
      end
!INPERL
      integer function inperl (idum)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE COMPUTES THE NUMBER OF VECTOR ELEMENTS THAT CAN
!     BE PRINTED IN A LINE OF OUTPUT ON THE STANDARD OUTPUT FILE.
!
!     ASSUMPTIONS RE -
!
!        1) MAXIMUM WIDTH OF LINE TO USE (IMAXW) IS 132.
!        2) NUMBER OF CHARACTERS NOT VECTOR ELEMENTS PER LINE
!                (IOCPL) IS 15.
!        2) WIDTH OF FIELD FOR AN ELEMENT, INCLUDING SPACING
!                BETWEEN ELEMENTS (IEW) IS 15.
!        4) MAXIMUM ELEMENTS PER LINE (IMAXE) IS 7.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!                       EXTRACTED FROM EARLIER LSTVEC.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   idum
!
!  LOCAL SCALARS
     integer&
     &   iew,imaxe,imaxw,iocpl,iwidth
!
!  INTRINSIC FUNCTIONS
      intrinsic min
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IDUM
!        INPUT PARAMETER.  UNUSED ARGUMENT.
!     INTEGER IEW
!        WIDTH OF A FIELD FOR PRINTING OUT A VECTOR ELEMENT,
!        INCLUDING SPACES BETWEEN ADJACENT ELEMENTS.
!     INTEGER IMAXE
!        MAXIMUM NUMBER OF ARRAY ELEMENTS PER LINE.
!     INTEGER IMAXW
!        MAXIMUM NUMBER OF CHARACTERS TO ALLOW PER LINE.
!     INTEGER IOCPL
!        NUMBER OF CHARACTERS TO BE INTRODUCED TO LINE IN ADDITION
!        TO CHARACTERS IN THE ELEMENT FIELDS.
!     INTEGER IWIDTH
!        NUMBER OF CHARACTERS IN A LINE ON THE STANDARD OUTPUT FILE.
!
!
!     INITIALIZATIONS
!
      data iew /15/, imaxe /7/, imaxw /132/, iocpl /15/
!
!     COMMENCE BODY OF ROUTINE
!
      iwidth = 132
      inperl = (min(iwidth, imaxw) - iocpl)/iew
      inperl = min(inperl, imaxe)
      return
      end
!EISII
     subroutine eisii(nmsub, nmvar, ival, ivalmn, ivalmx, msgtyp,&
     &   head, error, nmmin, nmmax)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THE ROUTINE CHECKS WHETHER THE VALUE   IVAL   IS WITHIN THE
!     THE RANGE IVALMN (INCLUSIVE) TO IVALMX (INCLUSIVE), AND PRINTS A
!     DIAGNOSTIC IF IT IS NOT.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ival,ivalmn,ivalmx,msgtyp
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmmax(8)*1,nmmin(8)*1,nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER IVAL
!        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
!     INTEGER IVALMN, IVALMX
!        THE MINIMUM AND MAXIMUM OF THE RANGE WITHIN WHICH THE
!        ARGUMENT MUST LIE.
!     INTEGER MSGTYP
!        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF ERROR IS .TRUE. AND
!        MSGTYP = 1 THE INPUT VALUE WAS OUTSIDE THE RANGE DETERMINED
!                   FROM OTHER INPUT ARGUMENTS
!        MSGTYP = 2 THE INPUT VALUE WAS OUTSIDE THE RANGE IMPOSED BY
!                   STARPAC
!     CHARACTER*1 NMMAX(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MAXIMUM.
!     CHARACTER*1 NMMIN(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR(8)
!        THE CHARACTERS OF THE ARGUMENTS NAME.
!
      error = .false.
!
     if (((ivalmn.le.ival) .and. (ival.le.ivalmx)) .or.&
     &   (ivalmx.lt.ivalmn)) return
!
      error = .true.
      call iprint(iprt)
      call ehdr(nmsub, head)
!
      if (msgtyp.le.2) write (iprt, 1000) (nmvar(i),i=1,6), ival
!
!     PRINT MESSAGE FOR VALUE OUTSIDE OF RANGE DETERMINED FROM
!     OTHER INPUT ARGUMENTS.
!
     if (msgtyp .eq. 1)&
    &   write (iprt, 1010) (nmvar(i),i=1,6), (nmmin(i),i=1,8),&
     &      (nmmax(i),i=1,8)
!
!     PRINT MESSAGE FOR VALUE OUTSIDE OF RANGE IMPOSED BY STARPAC
!
     if (msgtyp .eq. 2)&
     &   write (iprt, 1020) (nmvar(i),i=1,6), ivalmn, ivalmx
!
!     PRINT MESSAGE FOR AOV ROUTINES
!
     if (msgtyp .eq. 3)&
     &   write (iprt, 1030)
      return
!
!     FORMAT STATEMENTS
!
 1000 format (/20h the input value of , 6a1, 4h is , i6, '.')
1010 format(&
    &   27h the value of the argument , 6a1,&
    &   16h must be between, 1x, 8a1,&
     &   5h and , 8a1, 12h, inclusive.)
1020 format(&
    &   27h the value of the argument , 6a1,&
    &   16h must be between, 1x, i6,&
     &   5h and , i6, 12h, inclusive.)
1030 format(/' THE NUMBER OF DISTINCT GROUPS (NG) MUST BE BETWEEN'/&
     &  ' TWO AND ONE LESS THAN THE NUMBER OF POSITIVE TAG VALUES.')
!
      end
!EIAGE
     subroutine eiage (nmsub, nmvar, ym, n, m, iym, ymmn, nvmx,&
     &   head, msgtyp, nv, error, nmmin)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS TO ENSURE THAT NO VALUES, OR ONLY A MAXIMUM
!     OF NVMX, ARE NOT GREATER THAN A SPECIFIED LOWER BOUND YMMN,
!     WITH NAME NMMIN.   THE CHECKING OPTION IS SPECIFIED
!     WITH MSGTYP.  IF AN ERROR IS FOUND, THE ERROR IS PRINTED AND
!     AN ERROR FLAG AND THE NUMBER OF VIOLATINS ARE RETURNED.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JUNE 10, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iym,m,msgtyp,n,nv,nvmx,ymmn
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     integer&
     &   ym(*)
     character&
     &   nmmin(8)*1,nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,j
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EIAGEP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IYM
!        THE FIRST DIMENSION OF THE ARRAY YM.
!     INTEGER J
!        AN INDEXING VARIABLE.
!     INTEGER M
!        THE NUMBER OF COLUMNS OF DATA IN YM.
!     INTEGER MSGTYP
!        THE INDICATOR ARGUMENT FOR THE TYPE OF MESSAGE.
!        IF (MSGTYP.GE.3) THE MESSAGE PRINTED WILL USE NMMIN
!        OTHERWISE IT WILL USE YMMN.
!        IF (MSGTYP = 1 OR 3) NO VIOLATIONS ARE ALLOWED.
!        IF (MSGTYP = 2 OR 4) THE NUMBER OF VIOLATIONS MUST
!                             BE LESS THAN   NVMX   .
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     CHARACTER*1 NMMIN(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING ROUTINES NAME.
!     CHARACTER*1 NMVAR(8)
!        THE CHARACTERS OF THE PARAMETERS NAME.
!     INTEGER NV
!        THE NUMBER OF VIOLATIONS FOUND.
!     INTEGER NVMX
!        THE LARGEST NUMBER OF VIOLATIONS ALLOWED.
!     INTEGER YM(IYM,M)
!        THE ARRAY BEING TESTED.
!     INTEGER YMMN
!        THE MINIMUM ACCEPTABLE VALUE.
!
      error = .false.
!
      if ((n.le.0) .or. (m.le.0)) return
!
!     CHECK FOR VIOLATIONS
!
      nv = 0
      do 5 i = 1, n
         do 1 j = 1, m
            if (ym(i+(j-1)*iym) .lt. ymmn) nv = nv + 1
    1    continue
    5 continue
!
      if (nv .le. nvmx) return
!
!     VIOLATIONS FOUND
!
      error = .true.
!
     call eiagep (nmsub, nmvar, ymmn, nvmx, head, msgtyp, nv,&
     &   nmmin)
!
      return
!
      end
!NLSKL
      subroutine nlskl(iskull, page, wide, nlhdr)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS A HEADING AND WARNING MESSAGES FOR
!     SERIOUS ERRORS DETECTED BY THE NONLINEAR LEAST SQUARES ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  APRIL 2, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     logical&
     &   page,wide
!
!  ARRAY ARGUMENTS
     integer&
     &   iskull(10)
!
!  SUBROUTINE ARGUMENTS
       external nlhdr
!
!  LOCAL SCALARS
     integer&
     &   iprt,isubhd
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     EXTERNAL NLHDR
!        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISKULL(10)
!        AN ERROR MESSAGE INDICATOR VARIABLE.
!     INTEGER ISUBHD
!        AN INTEGER VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER OR NOT THE OUTPUT
!        IS TO BEGIN ON A NEW PAGE.
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
!
      isubhd = 0
      call nlhdr(page, wide, isubhd)
!
      if (wide) then
         write (iprt,1010)
         write (iprt,1020)
!        WRITE (IPRT,1030)
!        WRITE (IPRT,1040)
!        WRITE (IPRT,1050)
         write (iprt,1000)
      end if
      write (iprt,1060)
!
!     VCV COMPUTATION NOT COMPLETED
!
      if (iskull(7).ne.0) write (iprt,1120)
!
!     MAXIMUM NUMBER OF ITERATIONS REACHED BEFORE CONVERGENCE
!
      if (iskull(6).ne.0) write (iprt,1100)
!
!     FALSE CONVERGENCE
!
      if (iskull(5).ne.0) write (iprt,1090)
!
!     MEANINGLESS VCV MATRIX
!
      if (iskull(4).ne.0) write (iprt,1080)
!
!     PROBLEM IS COMPUTATIONALLY SINGULAR
!
      if (iskull(3).ne.0) write (iprt,1070)
!
!     INITIAL RESIDUAL SUM OF SQUARES COMPUTATION OVERFLOWED
!
      if (iskull(2).ne.0) write (iprt,1110)
!
      return
!
!     FORMAT STATEMENTS
!
 1000 format (///)
1010 format (/48h  w      w     aa     rrrrrrr   n      n    iiii,&
    &   19h    n      n    ggg/31h  w      w    a  a    r     rr ,&
    &   38h nn     n     ii     nn     n   g    g/12h  w      w  ,&
    &   51h  a  a    r      r  n n    n     ii     n n    n  g/&
    &   59h  ww    ww   aa  aa   r     rr  n  n   n     ii     n  n   ,&
    &   4hn  g/47h   w    w    aaaaaa   rrrrrrr   n  nn  n     ii,&
     &   23h     n  nn  n  g  ggggg)
1020 format (49h   w ww w    a    a   r r       n   n  n     ii  ,&
    &   21h   n   n  n  g      g/29h   w ww w    a    a   r  r   ,&
    &   41h   n    n n     ii     n    n n  g      g/9h    w  w ,&
    &   59h   aa    aa  r   r     n     nn     ii     n     nn   g    ,&
    &   2hgg/49h    w  w    a      a  r    r    n      n    iiii ,&
     &   21h   n      n    gggg g/)
!1010 FORMAT (/30X, 48H  W      W     AA     RRRRRRR   N      N    IIII,
!    *   19H    N      N    GGG/30X, 31H  W      W    A  A    R     RR ,
!    *   38H NN     N     II     NN     N   G    G/30X, 12H  W      W  ,
!    *   51H  A  A    R      R  N N    N     II     N N    N  G/30X,
!    *   59H  WW    WW   AA  AA   R     RR  N  N   N     II     N  N   ,
!    *   4HN  G/30X, 47H   W    W    AAAAAA   RRRRRRR   N  NN  N     II,
!    *   23H     N  NN  N  G  GGGGG)
!1020 FORMAT (30X, 49H   W WW W    A    A   R R       N   N  N     II  ,
!    *   21H   N   N  N  G      G/30X, 29H   W WW W    A    A   R  R   ,
!    *   41H   N    N N     II     N    N N  G      G/30X, 9H    W  W ,
!    *   59H   AA    AA  R   R     N     NN     II     N     NN   G    ,
!    *   2HGG/30X, 49H    W  W    A      A  R    R    N      N    IIII ,
!    *   21H   N      N    GGGG G/)
!1030 FORMAT (1(34X, 3HXXX, 58X, 3HXXX/), 31X, 6('X'), 58X, 6('X')/31X,
!    *   7('X'), 56X, 7('X')/31X, 9('X'), 52X, 9('X')/36X, 5('X'), 17X,
!    *   '(', 14('-'), ')', 17X, 5('X')/38X, 5('X'), 14X, 2H((, 14X,
!    *   2H)), 14X, 5('X')/40X, 5('X'), 10X, 2H((, 18X, 2H)), 10X,
!    *   5('X')/41X, 5('X'), 8X, 2H((, 20X, 2H)), 8X, 5('X')/43X,
!    *   5('X'), 5X, 2H((, 22X, 2H)), 5X, 5('X')/44X, 5('X'), 3X, 2H((,
!    *   24X, 2H)), 3X, 5('X'))
!1040 FORMAT (46X, 7HXXXXX (, 26X, 7H) XXXXX/48X,
!    *   5HXXX((, 7X, 2HOO, 8X, 2HOO, 7X, 5H))XXX/49X, 3HXX(, 7X,
!    *   4HO  O, 6X, 4HO  O, 7X, 3H)XX/50X, 2HX(, 7X, 4HO  O, 6X,
!    *   4HO  O, 7X, 2H)X/51X, '(', 8X, 2HOO, 8X, 2HOO, 8X, ')'/2(51X,
!    *   '(', 28X, ')'/), 51X, '(', 11X, 6HOO  OO, 11X, ')'/51X, 2H((,
!    *   10X, 6HOO  OO, 10X, 2H))/52X, 2H((, 24X, 2H))/53X, '(', 24X,
!    *   ')'/54X, '(', 22X, ')')
!1050 FORMAT (55X, 4H(--(, 14X, 4H)--)/59X, '(', 12X, ')'/58X,
!    *   3HX((, 10X, 3H))X/56X, 5HXXXX(, 10X, 5H)XXXX/54X, 9HXXXXX (II,
!    *   15HIIIIIIII) XXXXX/53X, 5('X'), 2X, 12H(IIIIIIIIII), 2X, 5('X')
!    *   /51X, 5('X'), 4X, '(', 10X, ')', 4X, 5('X')/49X, 5('X'), 6X,
!    *   2H((, 8X, 2H)), 6X, 5('X')/48X, 5('X'), 8X, 10H(--------), 8X,
!    *   5('X')/46X, 5('X'), 30X, 5('X')/44X, 5('X'), 34X, 5('X')/43X,
!    *   5('X'), 36X, 5('X')/41X, 5('X'), 40X, 5('X')/40X, 4HXXXX, 44X,
!    *   4HXXXX/38X, 5('X'), 46X, 5('X')/36X, 5('X'), 50X, 5('X')/31X,
!    *   9('X'), 52X, 9('X')/31X, 7('X'), 56X, 7('X')/31X, 6('X'), 58X,
!    *   6('X')/1(34X, 3HXXX, 58X, 3HXXX))
 1060 format (22h **  error summary  **)
1070 format (/50h this model and data are computationally singular.,&
     &   29h check your input for errors.)
1080 format (/43h at least one of the standardized residuals, 6h could,&
    &   47h not be computed because the standard deviation, 8h of the ,&
    &   18hresidual was zero./37h the validity of the covariance matri,&
     &   18hx is questionable.)
1090 format (/46h the iterations do not appear to be converging,&
    &   13h to a minimum, 41h (false convergence), indicating that the,&
    &   12h convergence, 16h criteria stopss/22h and stopp may be too ,&
    &   35hsmall for the accuracy of the model, 17h and derivatives,,&
    &   52h that there is an error in the derivative matrix, or/&
    &   15h that the model, 39h is discontinuous near the current coef,&
     &   18hficient estimates.)
1100 format (/53h program did not converge in the number of iterations,&
     &   13h or number of, 32h model subroutine calls allowed.)
1110 format (/50h the residual sum of squares could not be computed,&
     &   19h using the starting, 26h model coefficient values.)
1120 format (/44h the variance-covariance matrix could not be,&
     &   26h computed at the solution.)
      end
!BFSER
     subroutine bfser(nmsub, n, lagmax, iccov, jccov, inlppc, jnlppc,&
    &   m, index1, index2, icspc2, iphas, nf, nw, lags,&
     &   ldstak, ldsmin, lyfft, nfft, option)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS IS THE ERROR CHECKING ROUTINE FOR THE TIME SERIES
!     FOURIER UNIVARIATE SPECTRUM ANALYSIS ROUTINES.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 2, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
    &   iccov,icspc2,index1,index2,inlppc,iphas,jccov,jnlppc,&
     &   lagmax,ldsmin,ldstak,lyfft,m,n,nf,nfft,nw
!
!  ARRAY ARGUMENTS
     integer&
     &   lags(*)
     logical&
     &   option(4)
     character&
     &   nmsub(6)*1
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i,nv
     logical&
     &   head
!
!  LOCAL ARRAYS
     logical&
     &   error(30)
     character&
    &   l1(8)*1,liccov(8)*1,licspc(8)*1,lindx1(8)*1,lindx2(8)*1,&
    &   linlpp(8)*1,liphas(8)*1,ljccov(8)*1,ljnlpp(8)*1,&
    &   llagmx(8)*1,llags(8)*1,llds(8)*1,llgmx1(8)*1,&
     &   llyfft(8)*1,lm(8)*1,ln(8)*1,lnf(8)*1,lnm1(8)*1,lnw(8)*1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EISGE,EISII,EISLE,EIVII
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR(30)
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VALUE.
!     INTEGER ICCOV
!        THE FIRST DIMENSION OF THE ARRAY CCOV.
!     INTEGER ICSPC2
!        THE FIRST DIMENSION OF THE ARRAY CSPC2.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF ERR01, NO ERRORS WERE DETECTED
!        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
!     INTEGER INDEX1, INDEX2
!        THE INDICES OF THE COVARIANCES OF THE TWO SERIES.
!     INTEGER INLPPC
!        THE FIRST DIMENSION OF THE ARRAY NLPPC.
!     INTEGER IPHAS
!        THE FIRST DIMENSION OF THE ARRAY PHAS.
!     INTEGER JCCOV
!        THE FIRST DIMENSION OF THE ARRAY CCOV.
!     INTEGER JNLPPC
!        THE FIRST DIMENSION OF THE ARRAY NLPPC.
!     INTEGER LAGMAX
!        THE MAXIMUM LAG VALUE TO BE USED.
!     INTEGER LAGS(NW)
!        THE ARRAY USED TO SPECIFY THE LAG WINDOW TRUNCATION
!        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
!     INTEGER LDSTAK
!        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
!     INTEGER LDSMIN
!        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
!     CHARACTER*1 LICCOV(8), LICSPC(8), LINDX1(8),
!    *   LINDX2(8), LINLPP(8), LIPHAS(8), LJCCOV(8), LJNLPP(8),
!    *   LLAGMX(8), LLAGS(8), LLDS(8), LLGMX1(8), LLYFFT(8), LM(8),
!    *   LN(8), LNF(8), LNM1(8), LNW(8), L1(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF THE ARGUMENT(S)
!        CHECKED FOR ERRORS.
!     INTEGER LYFFT
!        THE LENGTH OF VECTOR YFFT.
!     INTEGER M
!        THE NUMBER OF SERIES FOR WHICH THE COVARIANCES WERE
!        COMPUTED
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
!     INTEGER NF
!        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
!        TO BE COMPUTED.
!     INTEGER NFFT
!        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE ARRAY CONTAINING THE NAME OF THE USER CALLED SUBROUTINE.
!     INTEGER NV
!        THE NUMBER OF VIOLATIONS FOUND WHEN CHECKING VECTOR LAGS.
!     INTEGER NW
!        THE ARGUMENT USED TO DETERMINE THE NUMBER OF DIFFERENT
!        BANDWIDTHS TO BE USED.
!     LOGICAL OPTION(4)
!        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
!        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
!        OR NOT (FALSE).
!
!     SET UP NAME ARRAYS
!
     data liccov(1), liccov(2), liccov(3), liccov(4), liccov(5),&
    &   liccov(6), liccov(7), liccov(8) /'I','C','C','O','V',' ',' ',&
     &   ' '/
     data licspc(1), licspc(2), licspc(3), licspc(4), licspc(5),&
    &   licspc(6), licspc(7), licspc(8) /'I','C','S','P','C','2',' ',&
     &   ' '/
     data lindx1(1), lindx1(2), lindx1(3), lindx1(4), lindx1(5),&
    &   lindx1(6), lindx1(7), lindx1(8) /'I','N','D','E','X','1',' ',&
     &   ' '/
     data lindx2(1), lindx2(2), lindx2(3), lindx2(4), lindx2(5),&
    &   lindx2(6), lindx2(7), lindx2(8) /'I','N','D','E','X','2',' ',&
     &   ' '/
     data liphas(1), liphas(2), liphas(3), liphas(4), liphas(5),&
    &   liphas(6), liphas(7), liphas(8) /'I','P','H','A','S',' ',' ',&
     &   ' '/
     data linlpp(1), linlpp(2), linlpp(3), linlpp(4), linlpp(5),&
    &   linlpp(6), linlpp(7), linlpp(8) /'I','N','L','P','P','C',' ',&
     &   ' '/
     data ljccov(1), ljccov(2), ljccov(3), ljccov(4), ljccov(5),&
    &   ljccov(6), ljccov(7), ljccov(8) /'J','C','C','O','V',' ',' ',&
     &   ' '/
     data ljnlpp(1), ljnlpp(2), ljnlpp(3), ljnlpp(4), ljnlpp(5),&
    &   ljnlpp(6), ljnlpp(7), ljnlpp(8) /'J','N','L','P','P','C',' ',&
     &   ' '/
     data llagmx(1), llagmx(2), llagmx(3), llagmx(4), llagmx(5),&
    &   llagmx(6), llagmx(7), llagmx(8) /'L','A','G','M','A','X',' ',&
     &   ' '/
     data llags(1), llags(2), llags(3), llags(4), llags(5), llags(6),&
     &   llags(7), llags(8) /'L','A','G','S',' ',' ',' ',' '/
     data llgmx1(1), llgmx1(2), llgmx1(3), llgmx1(4), llgmx1(5),&
    &   llgmx1(6), llgmx1(7), llgmx1(8) /'L','A','G','M','A','X','+',&
     &   '1'/
     data llds(1), llds(2), llds(3), llds(4), llds(5), llds(6),&
     &   llds(7), llds(8) /'L','D','S','T','A','K',' ',' '/
     data ln(1), ln(2), ln(3), ln(4), ln(5), ln(6), ln(7), ln(8) /'N',&
     &   ' ',' ',' ',' ',' ',' ',' '/
     data lm(1), lm(2), lm(3), lm(4), lm(5), lm(6), lm(7), lm(8) /'M',&
     &   ' ',' ',' ',' ',' ',' ',' '/
     data lnf(1), lnf(2), lnf(3), lnf(4), lnf(5), lnf(6), lnf(7),&
     &   lnf(8) /'N','F',' ',' ',' ',' ',' ',' '/
     data lnm1(1), lnm1(2), lnm1(3), lnm1(4), lnm1(5), lnm1(6),&
     &   lnm1(7), lnm1(8) /'N','-','1',' ',' ',' ',' ',' '/
     data lnw(1), lnw(2), lnw(3), lnw(4), lnw(5), lnw(6), lnw(7),&
     &   lnw(8) /'N','W',' ',' ',' ',' ',' ',' '/
     data llyfft(1), llyfft(2), llyfft(3), llyfft(4), llyfft(5),&
    &   llyfft(6), llyfft(7), llyfft(8) /'L','Y','F','F','T',' ',' ',&
     &   ' '/
     data l1(1), l1(2), l1(3), l1(4), l1(5), l1(6), l1(7), l1(8) /'1',&
     &   ' ',' ',' ',' ',' ',' ',' '/
!
!     SET UP FOR ERROR CHECKING
!
      ierr = 0
      head = .true.
!
      do 10 i=1,30
         error(i) = .false.
   10 continue
!
!     CALL ERROR CHECKING ROUTINES
!
      call eisge(nmsub, ln, n, 17, 1, head, error(1), ln)
!
      if ((.not.option(3))) go to 20
!
     call eisii(nmsub, llagmx, lagmax, 1, n-1, 1, head, error(2), l1,&
     &   lnm1)
!
      call eisge(nmsub, lm, m, 2, 1, head, error(3), lm)
!
     call eisge(nmsub, liccov, iccov, lagmax+1, 3, head, error(4),&
     &   llgmx1)
!
      call eisge(nmsub, ljccov, jccov, m, 4, head, error(5), lm)
!
      if (option(2)) then
       call eisge(nmsub, linlpp, inlppc, lagmax+1, 3, head, error(6),&
     &     llgmx1)
!
        call eisge(nmsub, ljnlpp, jnlppc, m, 4, head, error(7), lm)
      end if
!
      call eisle(nmsub, lindx1, index1, m, 2, head, error(8), lm)
!
      call eisle(nmsub, lindx2, index2, m, 2, head, error(9), lm)
!
   20 call eisge(nmsub, llyfft, lyfft, nfft, 9, head, error(10), llyfft)
!
     if (option(1) .and. (.not.option(4))) call eisge(nmsub, llds,&
     &   ldstak, ldsmin, 9, head, error(15), llds)
!
      if (option(4)) go to 40
!
      do 30 i=1,15
         if (error(i)) go to 70
   30 continue
!
      return
!
   40 continue
!
      call eisge(nmsub, lnf, nf, 1, 1, head, error(16), lnf)
!
      call eisge(nmsub, lnw, nw, 1, 1, head, error(18), lnw)
!
      if (error(18)) go to 50
      if (option(3)) then
        call eivii(nmsub, llags, lags, nw, 1, lagmax, 0,&
     &      head, 4, nv, error(19), l1, llagmx)
      else
        call eivii(nmsub, llags, lags, nw, 1, n-1, 0,&
     &      head, 4, nv, error(19), l1, lnm1)
      end if
!
   50 continue
!
      call eisge(nmsub, licspc, icspc2, nf, 3, head, error(24), lnf)
!
      call eisge(nmsub, liphas, iphas, nf, 3, head, error(25), lnf)
!
     if (error(2) .or. error(16) .or. error(18) .or. error(19)) go to&
     &   70
!
      call eisge(nmsub, llds, ldstak, ldsmin, 9, head, error(30), llds)
!
      do 60 i=1,30
         if (error(i)) go to 70
   60 continue
!
      return
!
   70 continue
      ierr = 1
      return
!
      end
!AOV1HD
      subroutine aov1hd(iprt)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     A SUBROUTINE TO PRINT OUT THE HEADING FOR THE ONEWAY ANOVA
!     FAMILY, AND IS THE ONLY SOURCE FOR HEADINGS IN THAT FAMILY
!
!     AUTHOR -
!        JOHN E. KOONTZ
!        STATISTICAL ENGINEERING DIVISION
!        NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  MAY 17, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL VERSP
!
!  VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE OUTPUT LOGICAL UNIT NUMBER
!
      call versp(.true.)
      write (iprt,1000)
      return
 1000 format(///48x, 20hanalysis of variance//)
      end
!EISEQ
     subroutine eiseq(nmsub, nmvar1, nval, neq, msgtyp, head, error,&
     &   nmvar2)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE VALUE   NVAL   IS
!     OQUAL TO   NEQ  AND PRINTS A DIAGNOSTIC IF IT IS NOT.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 2, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   msgtyp,neq,nval
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1,nmvar1(8)*1,nmvar2(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER MSGTYP
!        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF ERROR IS TRUE AND
!        MSGTYP = 1 THE INPUT VALUE WAS NOT EQUAL TO THE NUMBER OF PARAM
!                   SPECIFIED BY MSPEC (ARIMA ESTIMATION AND FORECASTING
!     INTEGER NEQ
!        THE ACCEPTABLE VALUE FOR THE ARGUMENT BEING TESTED.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR1(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED.
!     CHARACTER*1 NMVAR2(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED
!        AGAINST.
!     INTEGER NVAL
!        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
!
      error = .false.
!
      if (nval .eq. neq) return
!
      error = .true.
!
      call iprint (iprt)
!
      call ehdr(nmsub, head)
!
      write (iprt, 1000) (nmvar1(i), i=1,6), nval
!
!     PRINT MESSAGE FOR ARIMA ROUTINES
!
      write (iprt, 1010) (nmvar1(i), i=1,6), neq
      return
!
!     FORMAT STATEMENTS
!
 1000 format (/20h the input value of , 6a1, 4h is , i5, '.')
1010 format(&
    &   27h the value of the argument , 6a1,&
    &   ' MUST BE GREATER THAN OR EQUAL TO'/&
    &   1x, i5, ' = ONE PLUS THE SUM OF MSPEC(1,J)+MSPEC(3,J) FOR',&
    &   ' J = 1, ..., NFAC,'/&
    &   6x, ' = ONE PLUS THE NUMBER OF AUTOREGRESSIVE PARAMETERS PLUS'/&
     &   9x, ' THE NUMBER OF MOVING AVERAGE PARAMETERS.')
!
      end
!PLTSYM
      subroutine pltsym(iptsym, i, j, isym, n, ipoint, line, icount)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SUPPLIES THE APPROPRIATE PLOT SYMBOL FOR
!     THE PLOT LINE.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 21, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   i,ipoint,iptsym,j,n
!
!  ARRAY ARGUMENTS
     integer&
     &   icount(103),isym(n)
     character&
     &   line(103)*1
!
!  LOCAL SCALARS
     integer&
     &   isymbl
!
!  LOCAL ARRAYS
     character&
     &   sym(30)*1,sym1(10)*1
!
!  INTRINSIC FUNCTIONS
      intrinsic max,min
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEXING VARIABLE.
!     INTEGER ICOUNT(103)
!        THE NUMBER OF PLOT SYMBOLS AT EACH LOCATION.
!     INTEGER IPOINT
!        THE LOCATION IN THE PLOT STRING OF THE VALUE BEING PLOTTED.
!     INTEGER IPTSYM
!        AN INDICATOR VARIABLE USED TO DESIGNATE THE TYPE
!        OF PLOT.  IF IPTSYM = 1, THE PLOT IS A SYMPLE PAGE
!        OR VERTICAL PLOT.  IF IPTSYM = 2, THE PLOT IS A SYMBOL
!        PLOT.  IF IPTSYM = 3, THE PLOT IS A MULTIVARIATE PLOT.
!     INTEGER ISYM(N)
!        VECTOR CONTAINING SYMBOL DESIGNATIONS FOR PLOTTING
!     INTEGER ISYMBL
!        THE INDEX OF THE PLOT SYMBOL TO BE USED.
!     INTEGER J
!        AN INDEX VARIABLE.
!     CHARACTER*1 LINE(103)
!        THE VECTOR USED FOR THE PLOT STRING.
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     CHARACTER*1 SYM(30), SYM1(10)
!        THE PLOT SYMBOLS.
!
     data sym( 1)/'+'/,sym( 2)/'.'/,sym( 3)/'*'/,sym( 4)/'-'/,&
    &     sym( 5)/'A'/,sym( 6)/'B'/,sym( 7)/'C'/,sym( 8)/'D'/,&
    &     sym( 9)/'E'/,sym(10)/'F'/,sym(11)/'G'/,sym(12)/'H'/,&
    &     sym(13)/'I'/,sym(14)/'J'/,sym(15)/'K'/,sym(16)/'L'/,&
    &     sym(17)/'M'/,sym(18)/'N'/,sym(19)/'O'/,sym(20)/'P'/,&
    &     sym(21)/'Q'/,sym(22)/'R'/,sym(23)/'S'/,sym(24)/'T'/,&
    &     sym(25)/'U'/,sym(26)/'V'/,sym(27)/'W'/,sym(28)/'Y'/,&
     &     sym(29)/'Z'/,sym(30)/'Z'/
     data sym1(1)/'1'/,sym1(2)/'2'/,sym1(3)/'3'/,sym1(4)/'4'/,&
    &     sym1(5)/'5'/,sym1(6)/'6'/,sym1(7)/'7'/,sym1(8)/'8'/,&
     &     sym1(9)/'9'/,sym1(10)/'X'/
!
      icount(ipoint) = icount(ipoint) + 1
      if (icount(ipoint) .eq. 1) go to 5
!
      isymbl = min(icount(ipoint), 10)
      line(ipoint) = sym1(isymbl)
      return
!
    5 continue
      go to (10, 20, 30), iptsym
!
   10 line(ipoint) = sym(1)
      return
!
   20 isymbl = min(29, max(1, isym(i)))
      line(ipoint) = sym(isymbl)
      return
!
   30 isymbl = min(29, max(1, j+4))
      line(ipoint) = sym(isymbl)
!
      return
      end
!PLINE
      subroutine pline(imin, imax, isymbl, line)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE DEFINES ONE LINE OF A PLOT STRING FOR THE
!     VERTICAL PLOT ROUTINES.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 21, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   imax,imin
     character&
     &   isymbl*1
!
!  ARRAY ARGUMENTS
     character&
     &   line(103)*1
!
!  LOCAL SCALARS
     integer&
     &   i
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEXING VARIABLE.
!     INTEGER IMAX
!        THE LARGEST LOCATION IN THE PLOT STRING BEING DEFINED.
!     INTEGER IMIN
!        THE SMALLEST LOCATION IN THE PLOT STRING BEING DEFINED.
!     CHARACTER*1 ISYMBL
!        THE PLOTTING SYMBOL BEING USED.
!     CHARACTER*1 LINE(103)
!        THE VECTOR USED FOR THE PLOT STRING.
!
      do 10 i = imin, imax
         line(i) = isymbl
   10 continue
      return
      end
!MODSUM
      subroutine modsum(nfac, mspect)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS THE MODEL SUMMARY FOR THE ARIMA ROUTINES
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 4, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   nfac
!
!  ARRAY ARGUMENTS
     integer&
     &   mspect(nfac,4)
!
!  LOCAL SCALARS
     integer&
     &   i,iprt,j
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER NFAC
!        THE NUMBER OF FACTORS IN THE MODEL
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER J
!        AN INDEX VARIABLE.
!     INTEGER MSPECT(NFAC,4)
!        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
!
!
      call iprint(iprt)
!
!     PRINT MODEL SPECIFICATION
!
      write(iprt, 1002) (i, (mspect(i,j),j=1,4), i=1,nfac)
!
      return
!
!     FORMAT STATEMENTS
!
1002 format(//&
    &   '    MODEL SPECIFICATION'//&
    &   '       FACTOR          (P     D     Q)    S'//&
     &   (7x, i6, 6x, 4i6))
      end
!DCKHDR
      subroutine dckhdr(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE
!     DERIVATIVE CHECKING ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!       THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (21h+derivative checking,,&
     &   10h continued)
 1010 format ('+', 23(1h*)/ 24h * derivative checking */ 1x, 23(1h*))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!NLERR
      subroutine nlerr (icnvcd, iskull)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS THE ERROR FLAG IERR BASED ON THE CONVERGENCE
!     CODE RETURNED BY NL2.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  APRIL 2, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   icnvcd
!
!  ARRAY ARGUMENTS
     integer&
     &   iskull(10)
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER ICNVCD
!        THE CONVERGENCE CODE FROM NL2.
!     INTEGER IERR
!        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
!        IF IERR .GE. 1, ERRORS WERE DETECTED.
!     INTEGER ISKULL(10)
!        AN ERROR MESSAGE INDICATOR VARIABLE.
!
!     INITIALIZE MESSAGE INDICATOR VARIABLE
!
      do 5 i = 1, 10
         iskull(i) = 0
    5 continue
!
!     SET ERROR FLAG
!
     go to (10, 10, 20, 20, 20, 20, 40, 50, 60, 60, 10, 30, 10, 10,&
     &   10), icnvcd
!
!     BAD VALUE
!
   10 ierr = 1
      return
!
!     ACCEPTABLE STOPPING CONDITION
!
   20 ierr = 0
      return
!
!     INITIAL VARIANCE COMPUTATION OVERFLOWS
!
   30 ierr = 2
      iskull(2) = 1
      return
!
!     SINGULAR CONVERGENCE
!
   40 ierr = 3
      iskull(3) = 1
      return
!
!     FALSE CONVERGENCE
!
   50 ierr = 5
      iskull(5) = 1
      return
!
!     ITERATION OR FUNCTION EVALUATION LIMIT
!
   60 ierr = 6
      iskull(6) = 1
      return
!
      end
!EIVII
     subroutine eivii (nmsub, nmvar, ivec, n, iveclb, ivecub, nvmx,&
     &   head, msgtyp, nv, error, nmmin, nmmax)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS FOR VALUES IN THE INPUT VECTOR IVEC
!     WHICH ARE OUTSIDE THE (INCLUSIVE) LIMITS IVECLB TO IVECUB, PRINTS
!     AN ERROR MESSAGE IF THE NUMBER OF VIOLATIONS EXCEEDS THE LARGEST
!     NUMBER OF VIOLATIONS ALLOWED, AND RETURNS THE NUMBER OF
!     VIOLATIONS AND AN ERROR FLAG INDICATING THE RESULTS.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JUNE 10, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iveclb,ivecub,msgtyp,n,nv,nvmx
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     integer&
     &   ivec(*)
     character&
     &   nmmax(8)*1,nmmin(8)*1,nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        THE VALUE RETURNED FROM THE ERROR CHECKING ROUTINES TO INDICATE
!        WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED (TRUE)
!        OR NOT (FALSE).
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER IVEC(N)
!        THE VECTOR BEING TESTED.
!     INTEGER IVECLB, IVECUB
!        THE (INCLUSIVE) RANGE THAT THE VECTOR IS BEING TESTED
!        AGAINST.
!     INTEGER MSGTYP
!        THE INDICATOR ARGUMENT FOR THE TYPE OF MESSAGE.
!        IF (MSGTYP.GE.4) THE MESSAGE PRINTED WILL USE NMMIN AND
!        NMMAX, OTHERWISE IT WILL USE IVECLB AND IVECUB.
!        IF (MSGTYP = 1 OR 4) NO VIOLATIONS ARE ALLOWED.
!        IF (MSGTYP = 2 OR 5) THE NUMBER OF VIOLATIONS MUST
!                             BE LESS THAN   NVMX   .
!        IF (MSGTYP = 3 OR 6) VIOLATIONS ARE COUNTED ONLY IF THE
!                             THE FIRST ELEMENT IS NOT IN VIOLATION.
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     CHARACTER*1 NMMAX(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MAXIMUM.
!     CHARACTER*1 NMMIN(8)
!        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR(8)
!        THE CHARACTERS OF THE ARGUMENTS NAME.
!     INTEGER NV
!        THE NUMBER OF VIOLATIONS FOUND.
!     INTEGER NVMX
!        THE LARGEST NUMBER OF VIOLATIONS ALLOWED.
!
      error = .false.
!
      if (n.le.0) return
      if (ivecub.lt.iveclb) return
!
!     TEST WHETHER TESTING IS NECESSRY
!
     if ((mod(msgtyp,3) .eq. 0) .and.&
     &    ((ivec(1) .lt. iveclb) .or. (ivec(1) .gt. ivecub))) return
!
!     CHECK FOR VIOLATIONS
!
      nv = 0
      do 5 i = 1, n
         if ((ivec(i).lt.iveclb) .or. (ivec(i).gt.ivecub)) nv = nv + 1
    5 continue
!
      if (nv .le. nvmx) return
!
!     VIOLATIONS FOUND
!
      error = .true.
      call iprint(iprt)
      call ehdr(nmsub, head)
!
     if (msgtyp.le.3)&
     &   write (iprt, 1000) (nmvar(i),i=1,6), iveclb, ivecub, nv
     if (msgtyp.ge.4)&
    &   write (iprt, 1005) (nmvar(i),i=1,6), (nmmin(i),i=1,8),&
     &   (nmmax(i),i=1,8), nv
!
      go to (10, 20, 30, 10, 20, 30), msgtyp
!
   10 write(iprt, 1010) (nmvar(i),i=1,6)
      return
!
   20 write(iprt, 1020) (nmvar(i),i=1,6), nvmx
      return
!
   30 write(iprt, 1030) (nmvar(i),i=1,6)
      return
!
!     FORMAT STATEMENTS
!
1000 format (/&
    &   32h the number of values in vector , 6a1,&
    &   19h outside the range , i6, 3h to/&
     &   1x, i6, 16h, inclusive, is , i6, '.')
1005 format (/&
    &   32h the number of values in vector , 6a1,&
    &   19h outside the range , 8a1, 3h to/&
     &   1x, 8a1, 16h, inclusive, is , i6, '.')
1010 format(&
    &   26h the values in the vector , 6a1,&
     &   31h must all be within this range.)
1020 format(&
    &   36h the number of values in the vector , 6a1,&
    &   19h outside this range/&
     &   19h must be less than , i5, '.')
1030 format(&
    &   34h if the first value of the vector , 6a1,&
    &   21h is within this range/&
     &   45h all of the values must be within this range.)
!
      end
!STKREL
      subroutine stkrel(number)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!  DE-ALLOCATES THE LAST (NUMBER) ALLOCATIONS MADE IN THE STACK
!  BY STKGET.
!
!  ERROR STATES -
!
!    1 - NUMBER .LT. 0
!    2 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
!    3 - ATTEMPT TO DE-ALLOCATE NON-EXISTENT ALLOCATION
!    4 - THE POINTER AT ISTAK(LNOW) OVERWRITTEN
!
!     THIS FUNCTION WAS ADAPTED FROM THE FRAMEWORK FUNCTION ISTKGT
!
!     ADAPTED BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 26, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   number
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  ARRAYS IN COMMON
      double precision dstak(12)
!
!  LOCAL SCALARS
     integer&
     &   in,iprt,lbook,lmax,lnow,lout,lused
!
!  LOCAL ARRAYS
     integer&
     &   istak(12)
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!  COMMON BLOCKS
      common /cstak/dstak
      common /errchk/ierr
!
!  EQUIVALENCES
      equivalence (dstak(1),istak(1))
      equivalence (istak(1),lout)
      equivalence (istak(2),lnow)
      equivalence (istak(3),lused)
      equivalence (istak(4),lmax)
      equivalence (istak(5),lbook)
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     DOUBLE PRECISION DSTAK(12)
!        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
!        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
!     INTEGER IN
!        ...
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISTAK(12)
!        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER LBOOK
!        THE NUMBER OF WORDS USED FOR BOOKEEPING.
!     INTEGER LMAX
!        THE MAXIMUM LENGTH OF THE STACK.
!     INTEGER LNOW
!        THE CURRENT ACTIVE LENGTH OF THE STACK.
!     INTEGER LOUT
!        THE NUMBER OF CURRENT ALLOCATIONS.
!     INTEGER LUSED
!        THE MAXIMUM VALUE OF ISTAK(2) ACHIEVED.
!     INTEGER NUMBER
!        THE NUMBER OF ALLOCATIONS TO BE FREED FROM THE STACK.
!
!
      if (lnow.lt.lbook.or.lnow.gt.lused.or.lused.gt.lmax) go to 20
!
      in = number
 10      if (in.eq.0) return
!
         if (lnow.le.lbook) go to 30
!
!     CHECK TO MAKE SURE THE BACK POINTERS ARE MONOTONE.
!
         if (istak(lnow).lt.lbook.or.istak(lnow).ge.lnow-1) go to 40
!
         lout = lout-1
         lnow = istak(lnow)
         in = in-1
         go to 10
!
!     PRINT ERROR MESSAGES
!
   20 ierr = 1
      call iprint(iprt)
      write (iprt, 1000)
      return
!
   30 ierr = 1
      call iprint(iprt)
      write (iprt, 1010)
      return
!
   40 ierr = 1
      call iprint(iprt)
      write (iprt, 1020) lout
      return
!
!     FORMAT STATEMENTS
!
1000 format (///18h ***** error *****//&
     &   50h dstak bookkeeping elements have been overwritten.)
1010 format (///18h ***** error *****//&
    &   52h attempt has been made to de-allocate a non-existant,&
     &   21h allocation in dstak.)
1020 format (///18h ***** error *****//&
    &   35h the pointer for allocation number , i3, 9h has been,&
     &   13h overwritten.)
!
      end
!STPHDR
      subroutine stphdr(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE
!     STEP SIZE SELECTION ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!       THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (32h+derivative step size selection,,&
     &   10h continued)
1010 format ('+', 34(1h*)/ 35h * derivative step size selection */&
     &   1x, 34(1h*))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!CORRHD
      subroutine corrhd(iprt, m, n)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     A SUBROUTINE TO PRINT OUT THE HEADING FOR THE CORRELATION FAMILY.
!
!     AUTHOR -
!        JOHN E. KOONTZ
!        STATISTICAL ENGINEERING DIVISION
!        NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  MAY 17, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iprt,m,n
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE OUTPUT LOGICAL UNIT NUMBER
!     INTEGER M
!        THE NUMBER OF VARIABLES
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS FOR EACH VARIABLE
!
      call versp(.true.)
      write (iprt,1000) m, n
      return
!
!     FORMAT STATEMENTS
!
1000 format (/25h correlation analysis for, i3, 15h variables with,&
     &   i5, 13h observations/)
      end
!EIVEO
      subroutine eiveo (nmsub, nmvar, ivec, n, even, head)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER EACH OF THE VALUES IN THE INPUT
!     VECTOR IVEC ARE EVEN (OR ODD) AND PRINTS A
!     DIAGNOSTIC IF THEY ARE NOT.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   n
     logical&
     &   even,head
!
!  ARRAY ARGUMENTS
     integer&
     &   ivec(*)
     character&
     &   nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL EVEN
!        AN INDICATOR VARIABLE DESIGNATING WHETHER THE VALUES OF IVEC
!        SHOULD BE EVEN (TRUE) OR NOT (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER IVEC(N)
!        THE VECTOR BEING TESTED.
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR(8)
!        THE CHARACTERS OF THE PARAMETERS NAME.
!
!     CHECK FOR VIOLATIONS
!
      do 10 i = 1, n
        if ((even .and. (mod(ivec(i), 2) .eq. 1)) .or.&
     &       ((.not.even) .and. (mod(ivec(i), 2) .eq. 1))) go to 20
   10 continue
!
      return
!
!     VIOLATIONS FOUND
!
   20 continue
!
      call iprint(iprt)
!
      call ehdr(nmsub, head)
!
      if (even) go to 40
!
      write (iprt, 1010) (nmvar(i), i = 1, 6)
      return
!
   40 continue
      write (iprt, 1020) (nmvar(i), i = 1, 6)
      return
!
!     FORMAT STATEMENTS
!
1010 format(/&
    &   26h the values in the vector , 6a1,&
    &   27h must all be odd.  the next/&
     &   53h larger integer will be used in place of even values.)
1020 format(/&
    &   26h the values in the vector , 6a1,&
    &   28h must all be even.  the next/&
     &   52h larger integer will be used in place of odd values.)
!
      end
!CCFER
     subroutine ccfer(nmsub, n, lagmax, ldstak, ldsmin, iccov, jccov,&
     &  inlppc, jnlppc, m, lyfft, nfft, iym, iymfft, isfft, islong)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS IS THE ERROR CHECKING ROUTINE FOR THE CCF FAMILY
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 2, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
    &   iccov,inlppc,iym,iymfft,jccov,jnlppc,lagmax,ldsmin,ldstak,&
     &   lyfft,m,n,nfft
     logical&
     &   isfft,islong
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i
     logical&
     &   head
!
!  LOCAL ARRAYS
     logical&
     &   err(15)
     character&
    &   liccov(8)*1,linlpp(8)*1,liym(8)*1,liymff(8)*1,&
    &   ljccov(8)*1,ljnlpp(8)*1,llagmx(8)*1,llds(8)*1,&
    &   llgmx1(8)*1,llyfft(8)*1,lm(8)*1,ln(8)*1,lnfft(8)*1,&
     &   lnm1(8)*1,lone(8)*1,lthree(8)*1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EISGE,EISII
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERR(15)
!        VALUES INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER ICCOV
!        THE FIRST DIMENSION OF THE ARRAY CCOV.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
!     INTEGER INLPPC
!        THE FIRST DIMENSION OF THE ARRAY NLPPC.
!     LOGICAL ISFFT
!        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
!        ROUTINE HAS SUFFIX F (ISFFT = TRUE) OR NOT (ISFFT = FALSE)
!     LOGICAL ISLONG
!        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
!        ROUTINE HAS SUFFIX S (ISLONG = TRUE) OR NOT (ISLONG = FALSE)
!     INTEGER IYM, IYMFFT
!        THE FIRST DIMENSION OF THE ARRAYS YM AND YMFFT, RESPECTIVELY.
!     INTEGER JCCOV, JNLPPC
!        THE SECOND DIMENSIONS OF THE ARRAYS CCOV AND NLPPC,
!        RESPECTIVELY.
!     INTEGER LAGMAX
!        THE MAXIMUM LAG VALUE REQUESTED.
!     INTEGER LDSMIN
!        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
!     INTEGER LDSTAK
!        THE LENGTH OF THE ARRAY DSTAK.
!     CHARACTER*1 LICCOV(8), LINLPP(8), LIYM(8), LIYMFF(8), LJCCOV(8),
!    *  LJNLPP(8), LLAGMX(8), LLDS(8), LLGMX1(8), LLYFFT(8),
!    *  LM(8), LN(8), LNFFT(8), LNM1(8), LONE(8), LTHREE(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
!        CHECKED FOR ERRORS.
!     INTEGER M
!        THE NUMBER OF SERIES BEING ANALYZED
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
!     INTEGER NFFT
!        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
!        SUBROUTINE.
!
!     SET UP NAME ARRAYS
!
     data&
    & liccov(1), liccov(2), liccov(3), liccov(4), liccov(5),&
     & liccov(6), liccov(7), liccov(8) /'I','C','C','O','V',' ',' ',' '/
     data&
    & linlpp(1), linlpp(2), linlpp(3), linlpp(4), linlpp(5),&
     & linlpp(6), linlpp(7), linlpp(8) /'I','N','L','P','P','C',' ',' '/
     data&
    & liym(1), liym(2), liym(3), liym(4), liym(5),&
     & liym(6), liym(7), liym(8) /'I','Y','M',' ',' ',' ',' ',' '/
     data&
    & liymff(1), liymff(2), liymff(3), liymff(4), liymff(5),&
     & liymff(6), liymff(7), liymff(8) /'I','Y','M','F','F','T',' ',' '/
     data&
    & ljccov(1), ljccov(2), ljccov(3), ljccov(4), ljccov(5),&
     & ljccov(6), ljccov(7), ljccov(8) /'J','C','C','O','V',' ',' ',' '/
     data&
    & ljnlpp(1), ljnlpp(2), ljnlpp(3), ljnlpp(4), ljnlpp(5),&
     & ljnlpp(6), ljnlpp(7), ljnlpp(8) /'J','N','L','P','P','C',' ',' '/
     data&
    & llagmx(1), llagmx(2), llagmx(3), llagmx(4), llagmx(5),&
     & llagmx(6), llagmx(7), llagmx(8) /'L','A','G','M','A','X',' ',' '/
     data&
    & llds(1), llds(2), llds(3), llds(4), llds(5),&
     & llds(6), llds(7), llds(8) /'L','D','S','T','A','K',' ',' '/
     data&
    & llgmx1(1), llgmx1(2), llgmx1(3), llgmx1(4), llgmx1(5),&
     & llgmx1(6), llgmx1(7), llgmx1(8) /'L','A','G','M','A','X','+','1'/
     data&
    & llyfft(1), llyfft(2), llyfft(3), llyfft(4), llyfft(5),&
     & llyfft(6), llyfft(7), llyfft(8) /'L','Y','F','F','T',' ',' ',' '/
     data&
    & lm(1), lm(2), lm(3), lm(4), lm(5),&
     & lm(6), lm(7), lm(8) /'M',' ',' ',' ',' ',' ',' ',' '/
     data&
    & ln(1), ln(2), ln(3), ln(4), ln(5),&
     & ln(6), ln(7), ln(8) /'N',' ',' ',' ',' ',' ',' ',' '/
     data&
    & lnm1(1), lnm1(2), lnm1(3), lnm1(4), lnm1(5),&
     & lnm1(6), lnm1(7), lnm1(8) /'(','N','-','1',')',' ',' ',' '/
     data&
    & lnfft(1), lnfft(2), lnfft(3), lnfft(4), lnfft(5),&
     & lnfft(6), lnfft(7), lnfft(8) /'N','F','F','T',' ',' ',' ',' '/
     data&
    & lone(1), lone(2), lone(3), lone(4), lone(5),&
     & lone(6), lone(7), lone(8) /'O','N','E',' ',' ',' ',' ',' '/
     data&
    & lthree(1), lthree(2), lthree(3), lthree(4), lthree(5),&
     & lthree(6), lthree(7), lthree(8) /'T','H','R','E','E',' ',' ',' '/
!
!     SET UP FOR ERROR CHECKING
!
      ierr = 0
      head = .true.
      do 10 i = 1, 15
        err(i) = .false.
   10 continue
!
!     CALL ERROR CHECKING ROUTINES
!
      call eisge(nmsub, ln, n, 3, 2, head, err(1), lthree)
!
      call eisge(nmsub, lm, m, 1, 2, head, err(2), lone)
!
      if (.not.err(1)) then
!
       call eisii(nmsub, llagmx, lagmax, 1, n-1, 1, head, err(3), lone,&
     &    lnm1)
!
        if (isfft) then
          if (islong) then
           call eisge(nmsub, liymff, iymfft, nfft, 3, head, err(4),&
     &        lnfft)
          else
           call eisge(nmsub, llyfft, lyfft, nfft, 3, head, err(4),&
     &        lnfft)
          end if
        else
          call eisge(nmsub, liym, iym, n, 3, head, err(4), ln)
        end if
!
        if (.not.err(3)) then
!
          if (islong) then
           call eisge(nmsub, liccov, iccov, lagmax+1, 3, head, err(5),&
     &        llgmx1)
           call eisge(nmsub, ljccov, jccov, m, 3, head, err(6),&
     &        llgmx1)
           call eisge(nmsub, linlpp, inlppc, lagmax+1, 3, head, err(7),&
     &        llgmx1)
           call eisge(nmsub, ljnlpp, jnlppc, m, 3, head, err(8),&
     &        llgmx1)
          end if
!
          call eisge(nmsub, llds, ldstak, ldsmin, 9, head, err(9), llds)
!
        end if
      end if
!
      do 20 i = 1, 15
        if (err(i)) ierr = 1
   20 continue
!
      return
!
      end
!ACFDTL
      subroutine acfdtl (ndf, nd, iod, ntimes)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS TITLING FOR ACORRD.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ndf,ntimes
!
!  ARRAY ARGUMENTS
     integer&
     &   iod(*),nd(*)
!
!  LOCAL SCALARS
     integer&
     &   i,iprt,istop
     character&
     &   icom*1,iper*1,ipunct*1
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEX VARIABLE.
!     CHARACTER*1 ICOM
!        THE HOLLERITH VALUE -,- (COMMA)
!     INTEGER IOD(NDF)
!        THE ORDER OF EACH OF THE DIFFERENCE FACTORS.
!     CHARACTER*1 IPER
!        THE HOLLERITH VALUE -.- (PERIOD)
!     INTEGER IPRT
!        THE UNIT NUMBER OF THE DEVICE USED FOR PRINTED
!        OUTPUT.
!     CHARACTER*1 IPUNCT
!        THE HOLLERITH VALUE OF EITHER COMMA OR PERIOD.
!     INTEGER ISTOP
!        ONE LESS THAN THE NUMBER OF DIFFERENCE FACTORS.
!     INTEGER ND(NDF)
!        THE ARRAY CONTAINING THE NUMBER OF TIMES THE DIFFERENCE
!        FACTORS ARE TO BE APPLIED.
!     INTEGER NDF
!        THE NUMBER OF DIFFERENCE FACTORS.
!     INTEGER NTIMES
!        THE NUMBER OF TIMES THE DIFFERENCING FACTOR HAS BEEN APPLIED.
!
      data icom/','/, iper/'.'/
!
      call iprint (iprt)
!
      if (ndf .le. 1) go to 10
!
      istop = ndf - 1
      ipunct = iper
      if (ntimes .ge. 1) ipunct = icom
      write(iprt, 1000)
      if (ndf .eq. 2)  write(iprt, 1001) nd(2), iod(2), iper
     if (ndf .ge. 3) write(iprt, 1001)&
     &   (nd(i), iod(i), icom, i = 1, istop), nd(ndf), iod(ndf), ipunct
      go to 20
!
   10 write(iprt, 1002)
!
   20 if (ntimes .eq. 0) return
!
      if (ndf .ge. 2) write(iprt, 1003) ntimes, iod(1)
      if (ndf .eq. 1) write(iprt, 1004) ntimes, iod(1)
      return
!
!     FORMAT STATEMENTS
!
 1000 format(//47h series analyzed is input series differenced by/)
 1001 format(3x, 3(i3, ' FACTOR(S) OF ORDER ', i3, a1, 1x)/)
 1002 format(//' SERIES ANALYZED IS ORIGINAL INPUT SERIES'/)
1003 format(4x, 34h and, in addition, differenced by , i3,&
     &   18h factors of order , i3, '.'//)
1004 format(4x, 16h differenced by , i3, 18h factors of order ,&
     &   i3, '.'//)
      end
!ERIODD
      subroutine eriodd(nmsub, nmvar, nval, msgtyp, head, error)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS ERROR TO TRUE IF THE VALUE   NVAL   IS NOT EVEN
!     OR ODD, AS SPECIFIED BY THE PARAMETER ODD.  IN ADDITION, IF THIS
!     IS THE FIRST ERROR FOUND FOR THE CALLING SUBROUTINE   NMSUB   , IE
!     IF   HEAD   IS TRUE, THEN A HEADING FOR THE CALLING SUBROUTINE
!     IS ALSO PRINTED OUT.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 26, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   msgtyp,nval
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1,nmvar(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER MSGTYP
!        A VARIABLE USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF
!        MSGTYP = 1, THE INPUT VALUE SHOULD BE ODD AND
!        MSGTYP = 2, THE INPUT VALUE SHOULD BE EVEN.
!     CHARACTER*1 NMSUB(6)
!        THE ARRAY CONTAINING THE NAME OF THE CALLING SUBROUTINE.
!     CHARACTER*1 NMVAR(8)
!        THE ARRAY CONTAINING THE NAME OF THE VARIABLE BEING CHECKED.
!     INTEGER NVAL
!        THE VALUE OF THE VARIABLE BEING CHECKED.
!
      error = .false.
!
      if (msgtyp .eq. 2) go to 10
!
!     CHECK FOR ODD
!
      if (mod(nval, 2) .eq. 1) return
!
      call iprint(iprt)
      call ehdr(nmsub, head)
      write(iprt, 1010) (nmvar(i), i = 1, 6), (nmvar(i), i = 1, 6), nval
      error = .true.
      return
!
   10 continue
!
!     CHECK FOR EVEN
!
      if (mod(nval, 2) .eq. 0) return
!
      call iprint(iprt)
      call ehdr(nmsub, head)
      write(iprt, 1020) (nmvar(i), i = 1, 6), (nmvar(i), i = 1, 6), nval
      error = .true.
      return
!
!     FORMAT STATEMENTS
!
1010 format(/&
    &   27h the value of the variable , 6a1,&
    &   34h must be odd.  the input value of , 6a1/&
     &    4h is , i5, '.')
1020 format(/&
    &   27h the value of the variable , 6a1,&
    &   35h must be even.  the input value of , 6a1/&
     &    4h is , i5, '.')
!
      end
!AMFER
     subroutine amfer(nmsub, n, npar, ldstak, ldsmin,&
     &  save, mspec, nfac, ifcst, nfcst)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS IS THE ERROR CHECKING ROUTINE FOR NONLINEAR LEAST SQUARES
!     ESTIMATION ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 2, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ifcst,ldsmin,ldstak,n,nfac,nfcst,npar
     logical&
     &   save
!
!  ARRAY ARGUMENTS
     integer&
     &   mspec(4,*)
     character&
     &   nmsub(6)*1
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i,np,nv
     logical&
     &   head
!
!  LOCAL ARRAYS
     logical&
     &   error(20)
     character&
    &   lifcst(8)*1,llds(8)*1,lmspec(8)*1,ln(8)*1,lnfac(8)*1,&
     &   lnfcst(8)*1,lnpar(8)*1,lone(8)*1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EIAGE,EISEQ,EISGE
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR(20)
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        THE VARIABLE USED TO INDICATE WHETHER A HEADING IS TO BE
!        PRINTED DURING A GIVEN CALL TO THE ITERATION REPORT (TRUE)
!        OR NOT (FALSE).
!     INTEGER IERR
!        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
!        IF IERR .GE. 1, ERRORS WERE DETECTED.
!     INTEGER LDSMIN
!        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
!     INTEGER LDSTAK
!        THE LENGTH OF THE ARRAY DSTAK.
!     CHARACTER*1 LIFCST(8), LLDS(8), LMSPEC(8), LN(8), LNFAC(8),
!    *  LNPAR(8), LNFCST(8), LONE(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
!        CHECKED FOR ERRORS.
!     INTEGER MSPEC(4,NFAC)
!        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS.
!     INTEGER NFAC
!        THE NUMBER OF FACTORS IN THE MODEL
!     INTEGER NFCST
!        THE NUMBER OF FORECASTS.
!     CHARACTER*1 NMSUB(6)
!        THE NAME OF THE ROUTINE CALLING THE ERROR CHECKING ROUTINE
!     INTEGER NPAR
!        THE NUMBER OF PARAMETERS IN THE MODEL.
!     INTEGER NV
!        *
!     LOGICAL SAVE
!        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
!        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
!        (FALSE).
!
!     SET UP NAME ARRAYS
!
     data lifcst(1), lifcst(2), lifcst(3), lifcst(4), lifcst(5),&
    &   lifcst(6), lifcst(7), lifcst(8)&
     &  /'I','F','C','S','T',' ',' ',' '/
     data llds(1), llds(2), llds(3), llds(4), llds(5), llds(6),&
     &   llds(7), llds(8) /'L','D','S','T','A','K',' ',' '/
     data lmspec(1), lmspec(2), lmspec(3), lmspec(4), lmspec(5),&
    &   lmspec(6), lmspec(7), lmspec(8)&
     &  /'M','S','P','C',' ',' ',' ',' '/
     data ln(1), ln(2), ln(3), ln(4), ln(5), ln(6), ln(7), ln(8) /'N',&
     &   ' ',' ',' ',' ',' ',' ',' '/
     data lnfac(1), lnfac(2), lnfac(3), lnfac(4), lnfac(5),&
     &   lnfac(6), lnfac(7), lnfac(8) /'N','F','A','C',' ',' ',' ',' '/
     data lnfcst(1), lnfcst(2), lnfcst(3), lnfcst(4), lnfcst(5),&
    &   lnfcst(6), lnfcst(7), lnfcst(8)&
     &  /'N','F','C','S','T',' ',' ',' '/
     data lnpar(1), lnpar(2), lnpar(3), lnpar(4), lnpar(5),&
    &   lnpar(6), lnpar(7), lnpar(8) /'N','P','A','R',' ',' ',' ',&
     &   ' '/
     data lone(1), lone(2), lone(3), lone(4), lone(5),&
     &   lone(6), lone(7), lone(8) /'1',' ',' ',' ',' ',' ',' ',' '/
!
!     ERROR CHECKING
!
      do 10 i=1,20
         error(i) = .false.
   10 continue
!
      ierr = 0
      head = .true.
!
      call eisge(nmsub, ln, n, 1, 2, head, error(1), lone)
!
      call eisge(nmsub, lnfac, nfac, 1, 2, head, error(2), lone)
!
     if (.not. error(2))&
    &  call eiage(nmsub, lmspec, mspec, 4, nfac, 4, 0, 0, head, 1, nv,&
     &  error(3), lmspec)
!
      if ((.not. error(2)) .and. (.not. error(3))) then
        np = 1
         do 15 i = 1, nfac
           np = np + mspec(1,i) + mspec(3,i)
   15   continue
        call eiseq(nmsub, lnpar, npar, np, 1, head, error(4), lnpar)
      end if
!
     if ((.not.error(1)) .and. (.not.error(2)) .and. (.not.error(3))&
    &   .and. (.not.error(4)) .and. (.not.error(5)))&
    &   call eisge(nmsub, llds, ldstak, ldsmin, 9, head, error(6),&
     &   llds)
!
     if (save)&
    &   call eisge(nmsub, lifcst, ifcst, nfcst, 3, head, error(15),&
     &   lnfcst)
!
      do 20 i=1,20
         if (error(i)) go to 30
   20 continue
      return
!
   30 continue
      ierr = 1
      return
!
      end
!NLHDRA
      subroutine nlhdra(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE NONLINEAR
!     LEAST SQUARES ESTIMATION ROUTINES THAT USE ANALYTIC
!     (USER-SUPPLIED) DERIVATIVES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (35h+nonlinear least squares estimation,&
     &   42h with user-supplied derivatives, continued)
1010 format ('+', 71(1h*)/&
    &   1x, 37h*  nonlinear least squares estimation,&
     &   34h with user-supplied derivatives  */ 1x, 71(1h*))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!EHDR
      subroutine ehdr(nmsub, head)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS THE HEADING FOR THE ERROR CHECKING ROUTINES.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     logical&
     &   head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING ROUTINES NAME.
!
      if (.not.head) return
!
      call iprint(iprt)
!
      call versp(.false.)
      write(iprt,1010)
      write (iprt, 1000) (nmsub(i), i=1,6)
      head = .false.
!
      return
!
!     FORMAT STATEMENTS
!
 1000 format (/31h error checking for subroutine , 6a1/ 1x, 37('-'))
 1010 format ('+', 18(1h*)/19h * error messages */1x, 18(1h*))
!
      end
!ICNTI
      integer function icnti (iv, niv, i)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE COUNTS THE NUMBER OF OCCURENCES OF I IN IV.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING LAB/BOULDER
!                  NATIONAL BUREAU OF STANDARDS
!
!     CREATION DATE  -  APRIL 20, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   i,niv
!
!  ARRAY ARGUMENTS
     integer&
     &   iv(niv)
!
!  LOCAL SCALARS
     integer&
     &   j
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        INPUT PARAMETER.  THE INTEGER TO COUNT OCCURENCES OF.
!     INTEGER IV(NIV)
!        INPUT PARAMETER.  THE VECTOR IN WHICH TO COUNT.
!     INTEGER J
!        LOOP PARAMETER.
!     INTEGER NIV
!        INPUT PARAMETER.  THE LENGTH OF IV.
!
!     COMMENCE BODY OF ROUTINE
!
      icnti = 0
      do 10 j = 1, niv
         if (iv(j) .eq. i) icnti = icnti + 1
   10 continue
      return
      end
!LLHDRG
      subroutine llhdrg(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE UNRESTRICTED
!     LINEAR LEAST SQUARES ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
!
      call iprint(iprt)
      if (page) write (iprt,1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt,1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (32h+linear least squares estimation,&
     &  ' WITH USER-SPECIFIED MODEL, CONTINUED')
1010 format ('+', 63('*')/&
    &   1x, 34h*  linear least squares estimation,&
     &   ' WITH USER-SPECIFIED MODEL  *'/ 1x, 63('*'))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!ERDF
      subroutine erdf(nmsub, ndf, iod, nd, n, head, error)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PERFORMS ERROR CHECKING FOR THE INPUT
!     VALUES USED TO SPECIFY DIFFERENCING ON A TIME SERIES
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   n,ndf
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     integer&
     &   iod(*),nd(*)
     character&
     &   nmsub(6)*1
!
!  LOCAL SCALARS
     integer&
     &   i,ier,iprt,mbod
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX VARIABLE.
!     INTEGER IER
!        AN ERROR INDICATOR.
!     INTEGER IOD(NDF)
!        THE VECTOR CONTAINING THE ORDERS OF EACH DIFFERENCE FACTOR.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER MBOD
!        THE MAXIMUM BACKORDER DUE TO DIFFERENCING.
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
!     INTEGER ND(NDF)
!        THE VECTOR CONTAINING THE NUMBER OF TIMES EACH DIFFERENCE
!        FACTOR IS APPLIED.
!     INTEGER NDF
!        THE NUMBER OF DIFFERENCE FACTORS TO BE APPLIED TO THE SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINE NAME.
!
      error = .false.
!
      if (ndf .ge. 0) go to 10
      call iprint(iprt)
      call ehdr(nmsub, head)
      write (iprt, 1001) ndf
      error = .true.
      return
!
   10 if (ndf .eq. 0) return
!
      ier = 0
      mbod = 0
      do 30 i = 1, ndf
         if (iod(i) .ge. 1 .and. nd(i) .ge. 1) go to 20
         ier = 1
         go to 40
   20    mbod = mbod + iod(i) * nd(i)
   30 continue
      if (mbod .le. n - 1) return
!
   40 continue
      call iprint(iprt)
      call ehdr(nmsub, head)
     if (ier .eq. 1)&
     &   write (iprt, 1002) (i, nd(i), iod(i), i = 1, ndf)
      if (ier .eq. 0 .and. mbod .ge. n) write (iprt, 1003) mbod, n
      error = .true.
      return
!
!     FORMAT STATEMENTS
!
1001 format(/44h the number of difference factors (ndf) must/&
    &   54h be greater than or equal to zero.  the input value of/&
     &   8h ndf is , i6, '.')
1002 format (/46h the order of each difference factor (iod) and/&
    &   56h number of times it is applied (nd) must be greater than/&
    &   52h equal to one.  the input values of these arrays are/&
    &   25h    dif. fact.   nd   iod/&
     &   (1x, i13, i5, i6))
1003 format (/50h the maximum backorder due to differencing (mbod),&
    &  /54h that is, the sum of nd(i)*iod(i), i = 1, 2, ..., ndf,/&
    &   59h must be less than or equal to n-1.  the computed value for/&
     &   9h mbod is , i6, 33h, while the input value for n is , i6, '.')
      end
!STKST
      integer function stkst (nfact)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE REPLACES INTEGER FUNCTION ISTKST IN THE FRAMEWORK
!     FOR USE WITH STARPAC.  RETURNS ONE OF FOUR STATISTICS ON THE
!     STATE OF THE CSTAK STACK.
!
!     IMPORTANT - THIS ROUTINE ASSUMES THAT THE STACK IS INITIALIZED.
!                 IT DOES NOT CHECK TO SEE IF IT IS.  IN FACT, THERE
!                 IS NO WAY THAT IT COULD CHECK.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 14, 1983
!        BASED ON FRAMEWORK ROUTINE ISTKST.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   nfact
!
!  ARRAYS IN COMMON
      double precision dstak(12)
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  LOCAL ARRAYS
     integer&
     &   istak(12),istats(4)
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!  COMMON BLOCKS
      common /cstak/dstak
!
!  EQUIVALENCES
      equivalence (dstak(1),istak(1))
      equivalence (istak(1),istats(1))
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     DOUBLE PRECISION DSTAK(12)
!        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER IPRT
!        THE NUMBER OF THE STANDARD OUTPUT UNIT.
!     INTEGER ISTAK(12)
!        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER ISTATS(4)
!        INTEGER ARRAY INCLUDING THE FOUR STACK STATISTICS.
!     INTEGER NFACT
!
!
!     COMMENCE BODY OF ROUTINE
!
      if (nfact .gt. 0 .and. nfact .lt. 6) go to 10
!
!     REPORT ERROR STATUS
!
      call iprint (iprt)
      write (iprt, 1000) iprt
      stkst = 0
      return
!
!     REPORT TRUE VALUE OF A STATISTIC, ASSUMING STACK IS
!     DEFINED.
!
   10 stkst = istats(nfact)
      return
!
!     FORMAT STATEMENTS
!
1000 format (///18h ***** error *****//&
     &   24h illegal stack statistic, i5, 11h requested.)
      end
!EISGE
     subroutine eisge(nmsub, nmvar1, nval, nmin, msgtyp, head, error,&
     &   nmvar2)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE VALUE   NVAL   IS GREATER THAN
!     OR EQUAL TO   NMIN   AND PRINTS A DIAGNOSTIC IF IT IS NOT.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   msgtyp,nmin,nval
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1,nmvar1(8)*1,nmvar2(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER I
!        AN INDEX ARGUMENT.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER MSGTYP
!        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF ERROR IS TRUE AND
!        MSGTYP = 1 THE INPUT VALUE WAS TOO SMALL BASED
!                   ON LIMITS IMPOSED BY STARPAC
!        MSGTYP = 2 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
!                   ARGUMENTS.
!        MSGTYP = 3 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
!                   ARGUMENTS, WHERE THE VALUE INDICATES THE FIRST
!                   DIMENSION OF A DIMENSIONED ARRAY
!                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
!                         ARRAY NAME PRECEDED BY THE LETTER I.  IF THE
!                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
!                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
!                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
!                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
!        MSGTYP = 4 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
!                   ARGUMENTS, WHERE THE VALUE INDICATES THE SECOND
!                   DIMENSION OF A DIMENSIONED ARRAY
!                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
!                         ARRAY NAME PRECEDED BY THE LETTER J.  IF THE
!                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
!                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
!                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
!                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
!        MSGTYP = 5 THE ARGUMENT BEING CHECKED IS LDSTAK.
!                   NO LONGER USED.
!        MSGTYP = 6 THE ARGUMENT INDICATES THE FIRST DIMENSION OF
!                   AN ARRAY BEING CHECKED AGAINST THE NUMBER OF
!                   UNFIXED PARAMETERS.
!        MSGTYP = 7 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
!                   ARGUMENTS, WHERE THE VALUE INDICATES THE
!                   DIMENSION OF A VECTOR.
!                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
!                         ARRAY NAME PRECEDED BY THE LETTER L.  IF THE
!                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
!                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
!                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
!                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
!        MSGTYP = 8 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
!                   ARGUMENTS, WHERE THE VALUE INDICATES THE
!                   DIMENSION OF THE VECTORS ACOV AND NLPPA.
!        MSGTYP = 9 THE INPUT VALUE WAS TOO SMALL BASED ON LIMITS
!                   IMPOSED BY STARPAC, WHERE THE VALUE INDICATES THE
!                   DIMENSION OF A VECTOR.
!                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
!                         ARRAY NAME PRECEDED BY THE LETTER L.  IF THE
!                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
!                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
!                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
!                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
!     INTEGER NMIN
!        THE MINIMUM ACCEPTABLE VALUE FOR THE ARGUMENT BEING TESTED.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR1(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED.
!     CHARACTER*1 NMVAR2(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED
!        AGAINST.
!     INTEGER NVAL
!        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
!
      error = .false.
!
      if (nval .ge. nmin) return
!
      error = .true.
!
      call iprint (iprt)
!
      call ehdr(nmsub, head)
!
      write (iprt, 1000) (nmvar1(i), i=1,6), nval
!
      go to (20, 30, 40, 50, 60, 70, 80, 90, 100), msgtyp
!
!     PRINT MESSAGE FOR VALUE TOO SMALL BASED ON LIMITS IMPOSED
!     BY STARPAC.
!
   20 write (iprt, 1010) (nmvar1(i), i=1,6), nmin
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL BASED ON OTHER INPUT
!     ARGUMENTS.
!
   30 write (iprt, 1020) (nmvar1(i), i=1,6), (nmvar2(i), i=1,8)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     FIRST DIMENSION OF A DIMENSIONED ARRAY.
!
  40 write (iprt, 1030) (nmvar1(i), i=2,7), (nmvar1(i), i=1,6),&
     &   (nmvar2(i), i=1,8)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     SECOND DIMENSION OF A DIMENSIONED ARRAY.
!
  50 write (iprt, 1040) (nmvar1(i), i=2,7), (nmvar1(i), i=1,6),&
     &   (nmvar2(i), i=1,8)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHEN ARGUMENT IS LDSTAK.
!
   60 write(iprt, 1050) nmin
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     FIRST DIMENSION OF A DIMENSIONED ARRAY CHECK AGAINST THE NUMBER OF
!     UNFIXED PARAMETERS.
!
   70 write (iprt, 1060) (nmvar1(i), i=2,7), (nmvar1(i), i=1,6)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     DIMENSION OF A VECTOR.
!
  80 write (iprt, 1070) (nmvar1(i), i=2,7), (nmvar1(i), i=1,6),&
     &   (nmvar2(i), i=1,8)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     DIMENSION OF THE VECTORS ACOV AND NLPPA.
!
   90 write (iprt, 1080) (nmvar1(i), i=1,6), (nmvar2(i), i=1,8)
      return
!
!     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
!     DIMENSION OF A VECTOR.
!
 100 write (iprt, 1090) (nmvar1(i), i=2,7), (nmvar1(i), i=1,6),&
     &   nmin
      return
!
!     FORMAT STATEMENTS
!
 1000 format (/20h the input value of , 6a1, 4h is , i5, '.')
1010 format(&
    &   27h the value of the argument , 6a1,&
     &   34h must be greater than or equal to , i5, '.')
1020 format(&
    &   27h the value of the argument , 6a1,&
     &   34h must be greater than or equal to , 8a1, '.')
1030 format(&
    &   24h the first dimension of , 6a1,&
    &   30h, as indicated by the argument/&
     &    1x, 6a1, 35h, must be greater than or equal to , 8a1, '.')
1040 format(&
    &   25h the second dimension of , 6a1,&
    &   30h, as indicated by the argument/&
     &    1x, 6a1, 35h, must be greater than or equal to , 8a1, '.')
1050 format(&
    &   55h the dimension of the double precision vector dstak, as,&
    &   13h indicated by/&
    &   54h the argument ldstak, must be greater than or equal to,&
     &   i5, '.')
1060 format(&
    &   24h the first dimension of , 6a1,&
    &   30h, as indicated by the argument/&
    &    1x, 6a1, 34h, must be greater than or equal to,&
     &   34h the number of unfixed parameters.)
1070 format(&
    &   15h the length of , 6a1,&
    &   30h, as indicated by the argument/&
     &    1x, 6a1, 35h, must be greater than or equal to , 8a1, '.')
1080 format(&
    &   29h the length of acov and nlppa,&
    &   30h, as indicated by the argument/&
     &    1x, 6a1, 35h, must be greater than or equal to , 8a1, '.')
1090 format(&
    &   15h the length of , 6a1,&
    &   30h, as indicated by the argument/&
     &    1x, 6a1, 35h, must be greater than or equal to , i6, '.')
!
      end
!GENI
      subroutine geni(ivect, n, iinit, istp)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     PUT VALUES IINIT STEP ISTP THROUGH IINIT + (N - 1)*ISTP INTO
!     A VECTOR IVECT OF LENGTH N.  NO ERROR CHECKING IS DONE.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING LAB/BOULDER
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  MAY 17, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iinit,istp,n
!
!  ARRAY ARGUMENTS
     integer&
     &   ivect(n)
!
!  LOCAL SCALARS
     integer&
     &   i,j
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        INITIALIZATION VALUE.
!     INTEGER IINIT, ISTP
!        INPUT PARAMETERS.  THE INITIAL VALUE AND THE INCREMENT USED
!        IN CREATING THE INITIALIZATION VALUES.
!     INTEGER IVECT(N)
!        OUTPUT PARAMETER.  THE VECTOR INTO WHICH TO PUT THE VALUES
!        IINIT, IINIT + ISTP, ..., IINIT + (N - 1)*ISTP.
!     INTEGER J
!        LOOP PARAMETER.
!     INTEGER N
!        INPUT PARAMETER.  THE LENGTH OF IVECT.
!
      i = iinit
      do 10 j=1,n
         ivect(j) = i
         i = i + istp
   10 continue
      return
      end
!EISLE
     subroutine eisle(nmsub, nmvar1, nval, nmax, msgtyp, head, error,&
     &   nmvar2)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE VALUE   NVAL   IS LESS THAN
!     OR EQUAL TO   NMAX   AND PRINTS A DIAGNOSTIC IF IT IS NOT.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   msgtyp,nmax,nval
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1,nmvar1(8)*1,nmvar2(8)*1
!
!  LOCAL SCALARS
     integer&
     &   i,iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER NMAX
!        THE MAXIMUM ACCEPTABLE VALUE FOR THE ARGUMENT BEING TESTED.
!     INTEGER MSGTYP
!        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
!        PRINTED, WHERE IF ERROR IS TRUE AND
!        MSGTYP = 1 THE INPUT VALUE WAS TOO LARGE BASED
!                   ON LIMITS IMPOSED BY STARPAC
!        MSGTYP = 2 THE INPUT VALUE WAS TOO LARGE BASED ON OTHER INPUT
!                   ARGUMENTS.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     CHARACTER*1 NMVAR1(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED.
!     CHARACTER*1 NMVAR2(8)
!        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED
!        AGAINST.
!     INTEGER NVAL
!        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
!
      error = .false.
!
      if (nval .le. nmax) return
!
      error = .true.
!
      call iprint (iprt)
!
      call ehdr(nmsub, head)
!
      write (iprt, 1000) (nmvar1(i), i=1,6), nval
!
      go to (10, 20), msgtyp
!
!     PRINT MESSAGE FOR VALUE TOO LARGE BASED ON LIMITS IMPOSED
!     BY STARPAC.
!
   10 write (iprt, 1010) (nmvar1(i), i=1,6), nmax
      return
!
!     PRINT MESSAGE FOR VALUE TOO LARGE BASED ON OTHER INPUT
!     ARGUMENTS.
!
   20 write (iprt, 1020) (nmvar1(i), i=1,6), (nmvar2(i), i=1,8)
      return
!
!     FORMAT STATEMENTS
!
 1000 format (/20h the input value of , 6a1, 4h is , i5, '.')
1010 format(&
    &   27h the value of the argument , 6a1,&
     &   31h must be less than or equal to , i5, '.')
1020 format(&
    &   27h the value of the argument , 6a1,&
     &   31h must be less than or equal to , 8a1, '.')
!
      end
!STKCLR
      subroutine stkclr (nall0)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE IS AN ADDITION TO THE FRAMEWORK AREA MANIPULATION
!     ROUTINES.  IT CLEARS ALL ALLOCATIONS MADE SINCE THE FIRST NALL0.
!     IT IS INTENDED FOR USE DURING ERROR OR FINAL EXITS FROM STARPAC
!     ROUTINES WHICH MAKE ALLOCATIONS, TO RELEASE ALL ALLOCATIONS
!     MADE SINCE THE NALL0 EXISTING ON ENTRY TO THE STARPAC ROUTINE,
!     WITHOUT KNOWING HOW MANY ALLOCATIONS MUST BE RELEASED.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   nall0
!
!  LOCAL SCALARS
     integer&
     &   nalln
!
!  EXTERNAL FUNCTIONS
     integer&
     &   stkst
!       EXTERNAL STKST
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL STKREL
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER NALL0
!        INPUT PARAMETER.  THE NUMBER OF ALLOCATIONS TO BE PRESERVED
!        WHEN ALL LATER ONES ARE RELEASED.
!     INTEGER NALLN
!        THE TOTAL NUMBER OF ALLOCATIONS EXISTING BEFORE ANY ARE
!        RELEASED.
!
!     COMMENCE BODY OF ROUTINE
!
      nalln = stkst(1)
      call stkrel (nalln - nall0)
      return
      end
!LSTLAG
      integer function lstlag (nlppa, lagmax, lacov)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE FINDS THE LAG VALUE OF THE LAST AUTOCOVARIANCE
!     COMPUTED BEFORE ONE COULD NOT BE COMPUTED DUE TO MISSING DATA.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   lacov,lagmax
!
!  ARRAY ARGUMENTS
     integer&
     &   nlppa(lacov)
!
!  LOCAL SCALARS
     integer&
     &   lag
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER LACOV
!        THE LENGTH OF THE VECTOR ACOV.
!     INTEGER LAG
!        THE INDEXING VARIABLE INDICATING THE LAG VALUE OF THE
!        AUTOCORRELATION.
!     INTEGER NLPPA(LACOV)
!        THE ARRAY CONTAINING THE NUMBERS OF LAGGED PRODUCT PAIRS
!        USED TO COMPUTE THE ACVF AT EACH LAG.
!
!     FIND THE LAST AUTOCORRELATION TO BE COMPUTED BEFORE
!     ONE COULD NOT BE COMPUTED DUE TO MISSING DATA
!
      lstlag = -1
      if (nlppa(1) .le. 0) return
      do 20 lag = 1, lagmax
         if (nlppa(lag + 1) .ge. 1) go to 20
         lstlag = lag - 1
         return
   20 continue
      lstlag = lagmax
      return
      end
!EISRNG
      subroutine eisrng (nmsub, iseed, iseedu, head)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE INPUT VARIABLE ISEED IS
!     WITHIN [0, 2**MDIG], AND, IF NONZERO, IS ODD.
!
!     IF ISEED IS WITHIN [0, 2**MDIG] THEN
!        ISEEDU = ISEED-MOD(ISEED,2)+1
!     ELSE
!        ISEEDU = MIN[ ABS(ISEED)-MOD(ABS(ISEED),2)+1, 2**(MDIG-1)-1]
!                 AND AN ERROR MESSAGE IS PRINTED.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    CENTER FOR COMPUTING AND APPLIED MATHEMATICS
!                    NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
!                    BOULDER, COLORADO
!
!     CREATION DATE  -  JANUARY 17, 1990
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   iseed,iseedu
     logical&
     &   head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1
!
!  LOCAL SCALARS
     integer&
     &   iprt,mdig
!
!  EXTERNAL FUNCTIONS
!      INTEGER
!     +   I1MACH
!      EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT
!
!  INTRINSIC FUNCTIONS
      intrinsic abs,min,mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISEED
!        THE VALUE OF THE SEED BEING TESTED.
!     INTEGER ISEEDU
!        THE VALUE OF THE SEED ACTUALLY USED BY NRAND AND NRANDC.
!     INTEGER MDIG
!        A LOWER BOUND ON THE NUMBER OF BINARY DIGITS AVAILABLE
!        FOR REPRESENTING INTEGERS, INCLUDING THE SIGN BIT.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!
      mdig = min(i1mach(8)+1,32)
!
!     CHECK FOR VIOLATIONS
!
     if ((iseed.eq.0) .or.&
    &    ((iseed.ge.1) .and.&
    &     (iseed.le.2**(mdig-1)-1) .and.&
     &     (mod(iseed,2).eq.1))) then
!
!     SUPPLIED SEED WILL BE USED
!
         iseedu = iseed
      else
!
!     VIOLATIONS FOUND
!
         iseedu = min( abs(iseed)+mod(abs(iseed),2)-1, 2**(mdig-1)-1)
         call iprint(iprt)
         call ehdr(nmsub, head)
         write (iprt, 1010) mdig-1,iseedu
      end if
!
      return
!
!     FORMAT STATEMENTS
!
1010 format(/&
    &   ' THE VALUE OF ISEED MUST BE BETWEEN 0 AND 2**',i2,' - 1,'/&
    &   ' INCLUSIVE, AND, IF ISEED IS NOT 0, ISEED MUST BE ODD.  THE'/&
    &   ' SEED ACTUALLY USED BY THE RANDOM NUMBER GENERATOR HAS BEEN'/&
     &   ' SET TO', i10,'.')
!
      end
!SETESL
      subroutine setesl(n, ndiv, nfft)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE COMPUTES THE SMALLEST VALUE OF NFFT WHICH
!     EQUALS OR EXCEEDS N + 2, SUCH THAT NFFT - 2 IS
!     1. DIVISIBLE BY NDIV,
!     2. HAS NO MORE THAN 11 PRIME FACTORS,
!     3. HAS NO PRIME FACTOR GREATER THAN 23, AND
!     4. THE PRODUCT OF THE SQUARE FREE PRIME FACTORS OF
!        (NFFT-2)/NDIV DO NOT EXCEED 210 IF NDIV = 2, AND
!                                    105 IF NDIV = 4.
!     THE VALUE OF NFFT THUS MEET THE REQUIREMENTS OF
!     THE EXTENDED LENGTH OF THE SERIES REQUIRED FOR ANY ROUTINE
!     USING THE SINGLETON FFT PROVIDING THE PROPER VALUE OF NDIV
!     IS CHOSEN.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   n,ndiv,nfft
!
!  LOCAL SCALARS
     integer&
     &   i,npf,nsfp
!
!  LOCAL ARRAYS
     integer&
     &   ipf(50),ipfexp(50)
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL FACTOR
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!       AN INDEX VARIABLE.
!     INTEGER IPF(50), IPFEXP(50)
!        THE VECTORS OF PRIME FACTORS OF NFFT AND THEIR EXPONENTS,
!        RESPECTIVELY, WHERE THE LENGTH OF THESE VECTORS IS
!        SUFFICIENT TO ACCOMODATE THE PRIME FACTORS OF AN INTEGER
!        UP TO 2 ** 128 (APPROXIMATELY 10 ** 40).
!     INTEGER N
!        THE NUMBER UPON WHICH NFFT IS BASED.
!     INTEGER NDIV
!        A REQUIRED FACTOR OF NFFT - 2.
!     INTEGER NFFT
!        THE RETURNED VALUE WHICH MEETS THE ABOVE DESCRIPTION.
!     INTEGER NPF
!        THE NUMBER OF PRIME FACTORS IN NFFT.
!     INTEGER NSFP
!        THE PRODUCT OF THE NON SQUARE FACTORS.
!
      nfft = n
      if (nfft.le.0) return
      if (mod(nfft, ndiv) .ne. 0) nfft = nfft + ndiv - mod(nfft, ndiv)
      nfft = nfft - ndiv
   20 nfft = nfft + ndiv
      call factor(nfft/ndiv, npf, ipf, ipfexp)
      if ((npf.ge.11) .or. (ipf(npf).gt.23)) go to 20
      nsfp = 1
      if (ndiv.eq.4) nsfp = 2
      do 30 i = 1, npf
         if (mod(ipfexp(i), 2).eq.1) nsfp = nsfp * ipf(i)
   30 continue
      if (nsfp .ge. 210) go to 20
!
      nfft = nfft + 2
!
      return
!
      end
!ENFFT
      subroutine enfft(nmsub, nfft, ndiv, n, lyfft, nfft2, head, error)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE CHECKS WHETHER THE VALUE NFFT IS SUCH THAT NFFT-2 IS
!     DIVISIBLE BY NDIV AND HAS NO PRIME FACTORS GREATER THAN 23, AND
!     THE PRODUCT OF THE SQUARE FREE PRIME FACTORS OF NFFT - 2 DO NOT
!     EXCEED 209, I.E., THE VALUE OF NFFT MEETS THE REQUIREMENTS OF
!     THE EXTENDED LENGTH OF THE SERIES REQUIRED FOR ANY ROUTINE
!     USING THE SINGLETON FFT PROVIDING THE PROPER VALUE OF NDIV
!     IS CHOSEN.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   lyfft,n,ndiv,nfft,nfft2
     logical&
     &   error,head
!
!  ARRAY ARGUMENTS
     character&
     &   nmsub(6)*1
!
!  LOCAL SCALARS
     integer&
     &   iprt,nfft1
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL EHDR,IPRINT,SETESL
!
!  INTRINSIC FUNCTIONS
      intrinsic max
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERROR
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER LYFFT
!        THE LENGTH OF THE VECTOR CONTAINING THE SERIES TO BE EXTENDED.
!     INTEGER N
!        THE ACTUAL NUMBER OF OBSERVATIONS IN THE SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
!     INTEGER NFFT
!        THE USER SUPPLIED EXTENDED SERIES LENGTH.
!     INTEGER NFFT1
!        THE MAXIMUM OF NFFT AND N+2.
!     INTEGER NFFT2
!        THE SMALLEST EXTENDED SERIES LENGTH WHICH EQUALS OR
!        EXCEEDS NFFT AND WHICH MEETS THE REQUIREMENTS OF
!        SINGLETONS FFT CODE.
!
      error = .false.
      call iprint (iprt)
!
      if (nfft .ge. n+2) go to 20
!
!     PRINT WARNING
!
      call ehdr(nmsub, head)
!
      write (iprt, 1050) n
!
   20 continue
      nfft1 = max(nfft, n+2)
      call setesl(nfft1-2, ndiv, nfft2)
!
      if (nfft .eq. nfft2) return
!
!     PRINT WARNING
!
      call ehdr(nmsub, head)
!
      write (iprt, 1020) nfft, nfft2
!
      if (nfft .gt. lyfft) go to 40
!
      write (iprt, 1030) nfft2
      return
!
   40 continue
!
      error = .true.
!
      write (iprt, 1040) nfft2, lyfft
      return
!
!     FORMAT STATEMENTS
!
1020 format (/&
    &   40h the input value of the parameter nfft (, i5,&
    &   15h) does not meet/&
    &   51h the requirements of singletons fft code.  the next,&
    &   13h larger value/&
     &   15h which does is , i5, '.')
1030 format (/&
    &   11h the value , i5, 37h will be used for the extended series,&
     &    8h length.)
1040 format (/&
    &   20h however, the value , i5, 27h exceeds the length lyfft (,&
    &   i5, 8h) of the/&
    &   58h vector yfft, and therefore cannot be used as the extended/&
     &   43h series length without redimensioning yfft.)
1050 format (/&
    &   56h the extended series length (nfft) must equal or exceed,/&
    &   45h the number of observations in the series (n=, i5,&
     &    9h) plus 2.)
!
      end
!UFSER
     subroutine ufser(nmsub, n, lagmax, lacov, nf, ispcf, nw,&
     &    lags, ldstak, ldsmin, lyfft, nfft, option)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS IS THE ERROR CHECKING ROUTINE FOR THE TIME SERIES
!     FOURIER UNIVARIATE SPECTRUM ANALYSIS ROUTINES.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  JUNE 10, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ispcf,lacov,lagmax,ldsmin,ldstak,lyfft,n,nf,nfft,nw
!
!  ARRAY ARGUMENTS
     integer&
     &   lags(*)
     logical&
     &   option(4)
     character&
     &   nmsub(6)*1
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  LOCAL SCALARS
     integer&
     &   i,nv
     logical&
     &   head
!
!  LOCAL ARRAYS
     logical&
     &   err(15)
     character&
    &   l1(8)*1,lispcf(8)*1,llacov(8)*1,llagmx(8)*1,llags(8)*1,&
    &   llds(8)*1,llgmx1(8)*1,llyfft(8)*1,ln(8)*1,lnf(8)*1,&
     &   lnm1(8)*1,lnw(8)*1
!
!  EXTERNAL SUBROUTINES
!       EXTERNAL EISGE,EISII,EIVII
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     LOGICAL ERR(15)
!        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
!        (FALSE).
!     LOGICAL HEAD
!        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
!        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
!        OF HEAD WILL BE CHANGED TO FALSE.
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF ERR01, NO ERRORS WERE DETECTED
!        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
!     INTEGER ISPCF
!         THE ACTUAL FIRST DIMENSION OF THE SPECTRUM ARRAYS.
!     INTEGER LACOV
!        THE LENGTH OF THE VECTOR ACOV.
!     INTEGER LAGMAX
!        THE MAXIMUM LAG VALUE TO BE USED.
!     INTEGER LAGS(NW)
!        THE ARRAY USED TO SPECIFY THE LAG WINDOW TRUNCATION
!        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
!     INTEGER LDSMIN
!        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
!     INTEGER LDSTAK
!        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
!     CHARACTER*1 LISPCF(8), LLACOV(8), LLAGMX(8),
!    *   LLAGS(8), LLGMX1(8), LLDS(8), LN(8), LNF(8), LNM1(8),
!    *   LNW(8), LLYFFT(8), L1(8)
!        THE ARRAY(S) CONTAINING THE NAME(S) OF THE ARGUMENT(S)
!        CHECKED FOR ERRORS.
!     INTEGER LYFFT
!        THE LENGTH OF THE VECTOR YFFT.
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
!     INTEGER NF
!        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
!        TO BE COMPUTED.
!     INTEGER NFFT
!        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
!     CHARACTER*1 NMSUB(6)
!        THE ARRAY CONTAINING THE NAME OF THE USER CALLED SUBROUTINE.
!     INTEGER NV
!        THE NUMBER OF VIOLATIONS FOUND WHEN CHECKING VECTOR LAGS.
!     INTEGER NW
!        THE ARGUMENT USED TO DETERMINE THE NUMBER OF DIFFERENT
!        BANDWIDTHS TO BE USED.
!     LOGICAL OPTION(4)
!        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
!        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
!        OR NOT (FALSE).
!
!     SET UP NAME ARRAYS
!
     data lispcf(1), lispcf(2), lispcf(3), lispcf(4), lispcf(5),&
    &   lispcf(6), lispcf(7), lispcf(8) /'I','S','P','C','F',' ',' ',&
     &   ' '/
     data llacov(1), llacov(2), llacov(3), llacov(4), llacov(5),&
    &   llacov(6), llacov(7), llacov(8) /'L','A','C','O','V',' ',' ',&
     &   ' '/
     data llagmx(1), llagmx(2), llagmx(3), llagmx(4), llagmx(5),&
    &   llagmx(6), llagmx(7), llagmx(8) /'L','A','G','M','A','X',' ',&
     &   ' '/
     data llags(1), llags(2), llags(3), llags(4), llags(5), llags(6),&
     &   llags(7), llags(8) /'L','A','G','S',' ',' ',' ',' '/
     data llgmx1(1), llgmx1(2), llgmx1(3), llgmx1(4), llgmx1(5),&
    &   llgmx1(6), llgmx1(7), llgmx1(8) /'L','A','G','M','A','X','+',&
     &   '1'/
     data llds(1), llds(2), llds(3), llds(4), llds(5), llds(6),&
     &   llds(7), llds(8) /'L','D','S','T','A','K',' ',' '/
     data ln(1), ln(2), ln(3), ln(4), ln(5), ln(6), ln(7), ln(8) /'N',&
     &   ' ',' ',' ',' ',' ',' ',' '/
     data lnf(1), lnf(2), lnf(3), lnf(4), lnf(5), lnf(6), lnf(7),&
     &   lnf(8) /'N','F',' ',' ',' ',' ',' ',' '/
     data lnm1(1), lnm1(2), lnm1(3), lnm1(4), lnm1(5), lnm1(6),&
     &   lnm1(7), lnm1(8) /'N','-','1',' ',' ',' ',' ',' '/
     data lnw(1), lnw(2), lnw(3), lnw(4), lnw(5), lnw(6), lnw(7),&
     &   lnw(8) /'N','W',' ',' ',' ',' ',' ',' '/
     data llyfft(1), llyfft(2), llyfft(3), llyfft(4), llyfft(5),&
    &   llyfft(6), llyfft(7), llyfft(8) /'L','Y','F','F','T',' ',' ',&
     &   ' '/
     data l1(1), l1(2), l1(3), l1(4), l1(5), l1(6), l1(7), l1(8) /'1',&
     &   ' ',' ',' ',' ',' ',' ',' '/
!
!     SET UP FOR ERROR CHECKING
!
      ierr = 0
      head = .true.
!
      do 10 i=1,15
         err(i) = .false.
   10 continue
!
!     CALL ERROR CHECKING ROUTINES
!
      call eisge(nmsub, ln, n, 17, 1, head, err(1), ln)
!
      if (option(4)) then
        call eisge(nmsub, lnf, nf, 1, 1, head, err(6), lnf)
       if (.not.err(6))&
     &     call eisge(nmsub, lispcf, ispcf, nf, 3, head, err(7), lnf)
        call eisge(nmsub, lnw, nw, 1, 1, head, err(8), lnw)
      end if
!
      if (.not.err(1)) then
        if (option(3)) then
         call eisii(nmsub, llagmx, lagmax, 1, n-1, 1, head, err(2),&
     &       l1, lnm1)
          if (.not.err(2)) then
            if (option(2)) then
             call eisge(nmsub, llacov, lacov, lagmax+1, 8, head,&
     &          err(3), llgmx1)
            else
             call eisge(nmsub, llacov, lacov, lagmax+1, 7, head,&
     &          err(3), llgmx1)
            end if
          end if
        end if
        if (.not.err(2)) then
         if (option(1))&
    &      call eisge(nmsub, llyfft, lyfft, nfft, 9, head, err(4),&
     &        llyfft)
!
          if (.not.err(8)) then
           if (option(4)) then
            if (option(3)) then
            call eivii(nmsub, llags, lags, nw, 1, lagmax, 0, head, 3,&
     &         nv, err(9), l1, llagmx)
            else
            call eivii(nmsub, llags, lags, nw, 1, n-1, 0, head, 3, nv,&
     &         err(9), l1, lnm1)
            end if
           end if
!
           if ((.not.err(6)) .and. (.not.err(9)))&
    &         call eisge(nmsub, llds, ldstak, ldsmin, 9, head, err(14),&
     &            llds)
          end if
        end if
      end if
!
      do 40 i=1,15
         if (err(i)) ierr = 1
   40 continue
!
      return
!
      end
!FACTOR
      subroutine factor(n, npf, ipf, ipfexp)
!
!     Latest revision  -  03/15/90  (JRD)
!
!     This routine factors an input integer  "N"  and returns
!     the number of prime factors in  "NPF"  , The value of the
!     prime factors in the vector   "PF"  , and the exponent
!     of each of the prime factors in the vector  "IPFEXP"  .
!     the elements of   "IPF"  are stored in increasing order.
!     the length of the vectors is sufficient to accomodate
!     the prime factors of an integer up to 2 ** 128 (approximately
!     10 ** 40).
!
!     This routine is adapted from the factoring routine given
!     in ACM Algorithm 467 (CACM, 1973, Vol. 16, No. 11, Page 692-694).
!
!     Adapted by  -  Janet R. Donaldson
!                    Statistical Engineering Division
!                    National Bureau of Standards, Boulder, Colorado
!
!     Creation Date  -  October 23, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer n,npf
!
!  ARRAY ARGUMENTS
      integer ipf(50),ipfexp(50)
!
!  LOCAL SCALARS
      integer idiv,ifcur,iquot,npart
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IDIV, IFCUR
!        VARIOUS VARIABLES USED TO FACTOR    N   .
!     INTEGER IPF(50), IPFEXP(50)
!        THE VECTORS OF PRIME FACTORS OF   N   , AND THEIR EXPONENTS,
!        RESPECTIVELY.
!     INTEGER IQUOT
!        A VARIABLE USED TO FACTOR   N   .
!     INTEGER N
!        THE VALUE TO BE FACTORED.
!     INTEGER NPART
!        A VARIABLE USED TO FACTOR   N   .
!     INTEGER NPF
!        THE NUMBER OF FACTORS FOUND IN    N   .
!
!  DETERMINE THE FACTORS OF N
!
      npf = 0
      ifcur = 0
      npart = n
      idiv = 2
   10 continue
      iquot = npart/idiv
      if (npart.ne.idiv*iquot) go to 40
      if (idiv.le.ifcur) go to 20
      npf = npf + 1
      ipf(npf) = idiv
      ifcur = idiv
      ipfexp(npf) = 1
      go to 30
   20 continue
      ipfexp(npf) = ipfexp(npf) + 1
   30 continue
      npart = iquot
      go to 10
   40 continue
      if (iquot.le.idiv) go to 60
      if (idiv.ge.3) go to 50
      idiv = 3
      go to 10
   50 continue
      idiv = idiv + 2
      go to 10
   60 continue
      if (npart.le.1) return
      if (npart.le.ifcur) go to 70
      npf = npf + 1
      ipf(npf) = npart
      ipfexp(npf) = 1
      return
   70 continue
      ipfexp(npf) = ipfexp(npf) + 1
!
      end subroutine factor
!MSGX
      subroutine msgx(ier, iprt)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS THE RETURNED AND EXPECTED VALUES FOR THE
!     ERROR FLAG IERR
!
!     WRITTEN BY -
!        LINDA MITCHELL
!        STATISTICAL ENGINEERING DIVISION
!        NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  MAY 17, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ier,iprt
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  COMMON BLOCKS
      common /errchk/ierr
!
!     VARIBLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IER
!        EXPECTED VALUE OF ERROR FLAG IERR
!     INTEGER IERR
!        RETURNED ERROR FLAG FOUND IN THE COMMON ERRCHK
!     INTEGER IPRT
!        LOGICAL OUTPUT DEVICE
!
!
!     PRINT MESSAGE
      write (iprt,1000) ier, ierr
!
      if (ier.ne.ierr) write (iprt,1010)
!
      return
!
!     FORMAT STATEMENT
!
1000 format(/28h expected value for ierr is , i1/15h returned value,&
     &   12h for ierr is, i2)
 1010 format(48h possible error, unexpected value for error flag)
      end
!CPYVII
      subroutine cpyvii(n,x,incx,y,incy)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     COPY INTEGER X TO INTEGER Y.
!     FOR I = 0 TO N-1, COPY  X(LX+I*INCX) TO Y(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!     MODELED AFTER BLAS COPY ROUTINES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   incx,incy,n
!
!  ARRAY ARGUMENTS
     integer&
     &   x(n),y(n)
!
!  LOCAL SCALARS
     integer&
     &   i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEX VALUE.
!     INTEGER INCX
!        THE INCREMENT FOR THE MATRIX X.
!     INTEGER INCY
!        THE INCREMENT FOR THE MATRIX Y.
!     INTEGER N
!        THE NUMBER OF ROWS OF DATA TO BE COPIED FROM MATRIX X.
!     INTEGER X(N)
!        THE MATRIX TO BE COPIED FROM.
!     INTEGER Y(N)
!        THE MATRIX TO BE COPIED TO.
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        y(i) = x(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          y(i) = x(i)
   70     continue
      return
      end
!STKGET
      integer function stkget(nitems, itype)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!  ALLOCATES SPACE OUT OF THE INTEGER ARRAY ISTAK (IN COMMON
!  BLOCK CSTAK) FOR AN ARRAY OF LENGTH NITEMS AND OF TYPE
!  DETERMINED BY ITYPE AS FOLLOWS
!
!    1 - LOGICAL
!    2 - INTEGER
!    3 - REAL
!    4 - DOUBLE PRECISION
!    5 - COMPLEX
!
!  ON RETURN, THE ARRAY WILL OCCUPY
!
!    STAK(STKGET), STAK(STKGET+1), ..., STAK(STKGET-NITEMS+1)
!
!  WHERE STAK IS AN ARRAY OF TYPE ITYPE EQUIVALENCED TO ISTAK.
!
!  (FOR THOSE WANTING TO MAKE MACHINE DEPENDENT MODIFICATIONS
!  TO SUPPORT OTHER TYPES, CODES 6, 7, 8, 9, 10, 11 AND 12 HAVE
!  BEEN RESERVED FOR 1/4 LOGICAL, 1/2 LOGICAL, 1/4 INTEGER,
!  1/2 INTEGER, QUAD PRECISION, DOUBLE COMPLEX AND QUAD
!  COMPLEX, RESPECTIVELY.)
!
!  THE USE OF THE FIRST FIVE WORDS IS DESCRIBED BELOW.
!
!    ISTAK( 1) - LOUT,  THE NUMBER OF CURRENT ALLOCATIONS.
!    ISTAK( 2) - LNOW,  THE CURRENT ACTIVE LENGTH OF THE STACK.
!    ISTAK( 3) - LUSED, THE MAXIMUM VALUE OF ISTAK(2) ACHIEVED.
!    ISTAK( 4) - LMAX,  THE MAXIMUM LENGTH THE STACK.
!    ISTAK( 5) - LBOOK, THE NUMBER OF WORDS USED FOR BOOKEEPING.
!
!  THE NEXT FIVE WORDS CONTAIN INTEGERS DESCRIBING THE AMOUNT
!  OF STORAGE ALLOCATED BY THE FORTRAN SYSTEM TO THE VARIOUS
!  DATA TYPES.  THE UNIT OF MEASUREMENT IS ARBITRARY AND MAY
!  BE WORDS, BYTES OR BITS OR WHATEVER IS CONVENIENT.  THE
!  VALUES CURRENTLY ASSUMED CORRESPOND TO AN ANS FORTRAN
!  ENVIRONMENT.  FOR SOME MINI-COMPUTER SYSTEMS THE VALUES MAY
!  HAVE TO BE CHANGED (SEE I0TK00).
!
!    ISTAK( 6) - THE NUMBER OF UNITS ALLOCATED TO LOGICAL
!    ISTAK( 7) - THE NUMBER OF UNITS ALLOCATED TO INTEGER
!    ISTAK( 8) - THE NUMBER OF UNITS ALLOCATED TO REAL
!    ISTAK( 9) - THE NUMBER OF UNITS ALLOCATED TO DOUBLE PRECISION
!    ISTAK(10) - THE NUMBER OF UNITS ALLOCATED TO COMPLEX
!
!     THIS FUNCTION WAS ADAPTED FROM THE FRAMEWORK FUNCTION ISTKGT
!
!     ADAPTED BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 26, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   itype,nitems
!
!  SCALARS IN COMMON
     integer&
     &   ierr
!
!  ARRAYS IN COMMON
      double precision dstak(12)
!
!  LOCAL SCALARS
     integer&
     &   i,iprt,lbook,lmax,lnow,lout,lused
!
!  LOCAL ARRAYS
     integer&
     &   isize(5),istak(12)
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!  INTRINSIC FUNCTIONS
      intrinsic max
!
!  COMMON BLOCKS
      common /cstak/dstak
      common /errchk/ierr
!
!  EQUIVALENCES
      equivalence (dstak(1),istak(1))
      equivalence (istak(1),lout)
      equivalence (istak(2),lnow)
      equivalence (istak(3),lused)
      equivalence (istak(4),lmax)
      equivalence (istak(5),lbook)
      equivalence (istak(6),isize(1))
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     DOUBLE PRECISION DSTAK(12)
!        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER I
!        THE LOCATION OF A POINTER TO THE END OF THE PREVIOUS ALLOCATION
!     INTEGER IERR
!        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
!        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
!        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
!        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISIZE(5)
!        THE NUMBER OF WORDS IN EACH OF THE VARIOUS DATA TYPES.
!     INTEGER ISTAK(12)
!        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER ITYPE
!        THE TYPE OF ARRAY OF LENGTH NITEMS TO BE ALLOCATED.
!     INTEGER LBOOK
!        THE NUMBER OF WORDS USED FOR BOOKEEPING.
!     INTEGER LMAX
!        THE MAXIMUM LENGTH OF THE STACK.
!     INTEGER LNOW
!        THE CURRENT ACTIVE LENGTH OF THE STACK.
!     INTEGER LOUT
!        THE NUMBER OF CURRENT ALLOCATIONS.
!     INTEGER LUSED
!        THE MAXIMUM VALUE OF ISTAK(2) ACHIEVED.
!     INTEGER NITEMS
!        THE LENGTH OF THE ARRAY OF ITYPE TO BE ALLOCATED.
!
!
      stkget = (lnow*isize(2)-1)/isize(itype) + 2
      i = ( (stkget-1+nitems)*isize(itype) - 1 )/isize(2) + 3
!
!  STACK OVERFLOW IS AN UNRECOVERABLE ERROR.
!
      if (i .le. lmax) go to 10
!
      ierr = 1
      call iprint(iprt)
      write(iprt, 1000)
      return
!
   10 continue
!
!  ISTAK(I-1) CONTAINS THE TYPE FOR THIS ALLOCATION.
!  ISTAK(I  ) CONTAINS A POINTER TO THE END OF THE PREVIOUS
!             ALLOCATION.
!
      istak(i-1) = itype
      istak(i  ) = lnow
      lout = lout+1
      lnow = i
      lused = max(lused, lnow)
!
      return
!
!     FORMAT STATEMENTS
!
 1000 format(20h dstak is too short.)
!
      end
!SETLAG
      subroutine setlag (n, lagmax)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS THE NUMBER OF AUTOCORRELATIONS TO BE
!     COMPUTED.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 21, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   lagmax,n
!
!  INTRINSIC FUNCTIONS
      intrinsic min
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER LAGMAX
!        THE NUMBER OF LAGS AT WHICH THE AUTOCOVARIANCES ARE TO BE
!        COMPUTED.
!     INTEGER N
!        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
!
      if (n .ge. 96)                 lagmax = min(n / 3, 100)
      if (33 .le. n .and. n .le. 95) lagmax = 32
      if (n .le. 32)                 lagmax = n - 1
      return
      end
!LDSCMP
     subroutine ldscmp (narr, nlog, nint, nreal, ndbl, ncmp,&
     &   flag, nfp, ldsmin)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     COMPUTES LDSMIN, THE MINIMUM NUMBER OF DOUBLE PRECISION LOCATIONS
!     NEEDED BY THE FRAMEWORK TO STORE NARR ARRAYS, COMPRISING NLOG
!     LOGICAL LOCATIONS, NINT INTEGER LOCATIONS, NREAL REAL LOCATIONS,
!     NDBL DOUBLE PRECISION LOCATIONS, AND NCMP COMPLEX LOCATIONS,
!     TOGETHER WITH THE NOVER OVERHEAD INTEGER LOCATIONS THAT THE
!     FRAMEWORK ALWAYS USES AND THE 3 OVERHEAD LOCATIONS THAT IT USES
!     PER ARRAY STORED.  (ALL THE LOCATIONS ARE ASSIGNED OUT OF THE
!     LABELED COMMON CSTAK, USING A STACK DISCIPLINE.)
!
!     IT IS ASSUMED, BASED UPON THE FORTRAN STANDARD (ANSI X3.9 1966),
!     THAT DOUBLE PRECISION AND COMPLEX DATA ELEMENTS ARE TWICE AS LONG
!     AS INTEGER AND LOGICAL ELEMENTS.
!
!     WRITTEN BY - JOHN E. KOONTZ
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 7, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ldsmin,narr,ncmp,ndbl,nfp,nint,nlog,nreal
     character&
     &   flag*1
!
!  LOCAL SCALARS
     integer&
     &   nover
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     CHARACTER*1 FLAG
!        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE NFP
!        ELEMENTS ARE REAL OR DOUBLE PRECISION, WHERE FLAG=S INDICATES
!        THE NFP ELEMENTS ARE REAL (SINGLE PRECISION), AND FLAG=D
!        INDICATES THE ELEMENTS ARE DOUBLE PRECISION.
!     INTEGER LDSMIN
!        OUTPUT PARAMETER.  THE MINIMUM NUMBER OF DOUBLE PRECISION
!        LOCATIONS IN CSTAK REQUIRED FOR THE QUANTITIES OF ARRAY
!        ELEMENTS AND ARRAYS SPECIFIED BY THE INPUT PARAMETERS.
!     INTEGER NARR
!        INPUT PARAMETER.  THE NUMBER OF ARRAYS TO BE STORED IN CSTAK.
!     INTEGER NCMP
!        INPUT PARAMETER.  THE NUMBER OF COMPLEX ELEMENTS IN THE
!        ARRAYS TO BE STORED IN CSTAK.
!     INTEGER NDBL
!        INPUT PARAMETER.  THE NUMBER OF DOUBLE PRECISION ELEMENTS IN
!        THE ARRAYS TO BE STORED, IN CSTAK.
!     INTEGER NFP
!        THE NUMBER OF ELEMENTS WHICH DEPEND ON THE PRECISION OF THE
!        VERSION OF STARPAC BEING USED.
!     INTEGER NINT
!        INPUT PARAMETER.  THE NUMBER OF INTEGER ELEMENTS IN THE
!        ARRAYS TO BE STORED IN CSTAK.
!     INTEGER NLOG
!        INPUT PARAMETER.  THE NUMBER OF LOGICAL ELEMENTS IN THE
!        ARRAYS TO BE STORED IN CSTAK.
!     INTEGER NOVER
!        THE NUMBER OF INTEGER LOCATIONS THAT THE FRAMEWORK ALWAYS
!        USES FOR OVERHEAD PURPOSES.
!     INTEGER NREAL
!        INPUT PARAMETER.  THE NUMBER OF REAL ELEMENTS IN THE ARRAYS
!        TO BE STORED IN CSTAK.
!
!     DEFINE CONSTANTS
!
      data nover /10/
!
!     COMMENCE BODY OF ROUTINE
!
     ldsmin = (nlog + nint + nreal + 3*narr + nover + 1)/2&
     &       + ndbl + ncmp
      if (flag.eq.'S') then
         ldsmin = ldsmin + (nfp+1)/2
      else
         ldsmin = ldsmin + nfp
      end if
      return
      end
!PRTCNT
      subroutine prtcnt(nprt, ndigit, iptout)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE SETS UP THE PRINT CONTROL PARAMETERS.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  DECEMBER 29, 1982
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   ndigit,nprt
!
!  ARRAY ARGUMENTS
     integer&
     &   iptout(ndigit)
!
!  LOCAL SCALARS
     integer&
     &   i,ifac1,ifac2
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I, IFAC1, IFAC2
!     INTEGER IPTOUT(NDIGIT)
!        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
!     INTEGER NDIGIT
!        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
!     INTEGER NPRT
!        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
!        TO BE PROVIDED.
!
!
      if (nprt.le.-1) go to 20
!
      ifac1 = 10 ** (ndigit)
      do 10 i = 1, ndigit
         ifac2 = ifac1/10
         iptout(i) = mod(nprt, ifac1) / ifac2
         ifac1 = ifac2
   10 continue
      return
!
   20 do 30 i = 1, ndigit
         iptout(i) = 1
   30 continue
      iptout (ndigit) = 2
!
      return
!
      end
!NLHDRN
      subroutine nlhdrn(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE NONLINEAR
!     LEAST SQUARES ESTIMATION ROUTINES THAT USE NUMERICAL
!     APPROXIMATIONS TO THE DERIVATIVES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 3, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
1000 format (35h+nonlinear least squares estimation,&
     &   53h with numerically approximated derivatives, continued)
1010 format ('+', 82(1h*)/&
    &   1x, 37h*  nonlinear least squares estimation,&
     &   45h with numerically approximated derivatives  */ 1x, 82(1h*))
 1020 format ('1')
 1030 format (//30h summary of initial conditions/ 1x, 30('-'))
      end
!AMFHDR
      subroutine amfhdr(page, wide, isubhd)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE NONLINEAR
!     LEAST SQUARES ESTIMATION ROUTINES FOR ARIMA MODELS THAT USE
!     NUMERICAL APPROXIMATIONS TO THE DERIVATIVES.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  AUGUST 1, 1985
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   isubhd
     logical&
     &   page,wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT,VERSP
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER FOR PRINTED OUTPUT.
!     INTEGER ISUBHD
!        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
!     LOGICAL PAGE
!        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
!        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
!     LOGICAL WIDE
!        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
!        BE FULL WIDTH (TRUE) OR NOT (FALSE).
!
      call iprint(iprt)
      if (page) write (iprt, 1020)
      call versp(wide)
      if (page) write (iprt,1000)
      if (.not.page) write (iprt,1010)
      page = .true.
!
      if (isubhd.eq.0) return
!
      go to (10), isubhd
!
   10 write (iprt, 1030)
!
      return
!
!     FORMAT STATEMENTS FOR PAGE HEADINGS
!
 1000 format ('+ARIMA FORECASTING, CONTINUED')
 1010 format ('+', 23(1h*)/ ' *  ARIMA FORECASTING  *', /1x, 23(1h*))
 1020 format ('1')
 1030 format (//' MODEL SUMMARY'/' -------------')
      end
!ICOPY
      subroutine icopy(n,isx,incx,isy,incy)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE IS A ADAPTATION OF THE BLAS SUBROUTINE SCOPY,
!     MODIFIED TO HANDLE INTEGER ARRAYS.
!
!     COPY INTEGER ISX TO INTEGER ISY.
!     FOR I = 0 TO N-1, COPY  ISX(LX+I*INCX) TO ISY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
!     WRITTEN BY  -  JANET R. DONALDSON
!                    STATISTICAL ENGINEERING DIVISION
!                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  APRIL 2, 1981
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   incx,incy,n
!
!  ARRAY ARGUMENTS
     integer&
     &   isx(n),isy(n)
!
!  LOCAL SCALARS
     integer&
     &   i,ix,iy,m,mp1,ns
!
!  INTRINSIC FUNCTIONS
      intrinsic mod
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER I
!        AN INDEXING VARIABLE.
!     INTEGER INCX, INCY
!        THE INCREMENT USED FOR THE COPY FROM ONE VARIABLE TO THE OTHER.
!     INTEGER ISX(N)
!        THE ARRAY TO BE COPIED FROM.
!     INTEGER ISY(N)
!        THE ARRAY TO BE COPIED TO.
!     INTEGER IX, IY
!        INDEX VARIABLES.
!     INTEGER M
!        THE VALUE OF N MODULO 7.
!     INTEGER MP1
!        THE VALUE OF M + 1.
!     INTEGER N
!        THE NUMBER OF OBSERVATIONS IN THE ARRAYS ISX AND ISY.
!     INTEGER NS
!        THE VALUE OF N * INCX.
!
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        isy(iy) = isx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        isy(i) = isx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        isy(i) = isx(i)
        isy(i + 1) = isx(i + 1)
        isy(i + 2) = isx(i + 2)
        isy(i + 3) = isx(i + 3)
        isy(i + 4) = isx(i + 4)
        isy(i + 5) = isx(i + 5)
        isy(i + 6) = isx(i + 6)
   50 continue
      return
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
      do 70 i=1,ns,incx
          isy(i) = isx(i)
   70 continue
      return
      end
!STKSET
      subroutine stkset (nitems, itype)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!  INITIALIZES THE STACK TO NITEMS OF TYPE ITYPE
!
!     THIS FUNCTION WAS ADAPTED FROM THE FRAMEWORK SUBROUTINE ISTKIN
!
!     ADAPTED BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  NOVEMBER 26, 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     integer&
     &   itype,nitems
!
!  ARRAYS IN COMMON
      double precision dstak(12)
!
!  LOCAL SCALARS
     integer&
     &   lbook,lmax,lnow,lout,lused
!
!  LOCAL ARRAYS
     integer&
     &   isize(5),istak(12)
!
!  INTRINSIC FUNCTIONS
      intrinsic max
!
!  COMMON BLOCKS
      common /cstak/dstak
!
!  EQUIVALENCES
      equivalence (dstak(1),istak(1))
      equivalence (istak(1),lout)
      equivalence (istak(2),lnow)
      equivalence (istak(3),lused)
      equivalence (istak(4),lmax)
      equivalence (istak(5),lbook)
      equivalence (istak(6),isize(1))
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     DOUBLE PRECISION DSTAK(12)
!        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER ISIZE(5)
!        THE NUMBER OF WORDS IN EACH OF THE VARIOUS DATA TYPES.
!     INTEGER ISTAK(12)
!        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
!     INTEGER ITYPE
!        THE TYPE OF ARRAY OF LENGTH NITEMS TO BE ALLOCATED.
!     INTEGER LBOOK
!        THE NUMBER OF WORDS USED FOR BOOKEEPING.
!     INTEGER LMAX
!        THE MAXIMUM LENGTH OF THE STACK.
!     INTEGER LNOW
!        THE CURRENT ACTIVE LENGTH OF THE STACK.
!     INTEGER LOUT
!        THE NUMBER OF CURRENT ALLOCATIONS.
!     INTEGER LUSED
!        THE MAXIMUM VALUE OF ISTAK(2) ACHIEVED.
!     INTEGER NITEMS
!        THE LENGTH OF THE ARRAY OF ITYPE TO BE ALLOCATED.
!
!  HERE TO INITIALIZE
!
!  SET DATA SIZES APPROPRIATE FOR A STANDARD CONFORMING
!  FORTRAN SYSTEM USING THE FORTRAN "STORAGE UNIT" AS THE
!  MEASURE OF SIZE.
!
!  LOGICAL
      isize(1) = 1
!  INTEGER
      isize(2) = 1
!  REAL
      isize(3) = 1
!  DOUBLE PRECISION
      isize(4) = 2
!  COMPLEX
      isize(5) = 2
!
      lbook = 10
      lnow  = lbook
      lused = lbook
      lmax  = max( (nitems*isize(itype))/isize(2), 12 )
      lout  = 0
!
      return
!
      end
!VERSP
      subroutine versp (wide)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!     THIS ROUTINE PRINTS THE VERSION NUMBER.
!
!     WRITTEN BY - JANET R. DONALDSON
!                  STATISTICAL ENGINEERING DIVISION
!                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
!
!     CREATION DATE  -  OCTOBER 4, 1983
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
     logical&
     &   wide
!
!  LOCAL SCALARS
     integer&
     &   iprt
!
!  EXTERNAL SUBROUTINES
!      EXTERNAL IPRINT
!
!     VARIABLE DEFINITIONS (ALPHABETICALLY)
!
!     INTEGER IPRT
!        THE UNIT NUMBER OF THE DEVICE USED FOR PRINTED OUTPUT.
!     LOGICAL WIDE
!        THE MAXIMUM NUMBER OF COLUMNS THE PRINTED OUTPUT CAN USE.
!
      call iprint(iprt)
!
      if (wide) then
         write(iprt, 1000)
      else
         write(iprt, 1010)
      end if
!
      return
!
!     FORMAT STATEMENTS
!
 1000 format (105x, 'STARPAC 3.00 (05/15/2022)')
 1010 format (54x, 'STARPAC 3.00 (05/15/2022)')
      end

!IMDCON
      integer function imdcon(k)
!
!     LATEST REVISION  -  03/15/90  (JRD)
!
!  ***  RETURN INTEGER MACHINE-DEPENDENT CONSTANTS  ***
!
!     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   ***
!     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  ***
!     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            ***
!          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer k
!
!  LOCAL ARRAYS
      integer mdcon(3)
!
!  EXTERNAL FUNCTIONS
!      integer i1mach
!      external i1mach
!
      mdcon(1) = i1mach(2)
      mdcon(2) = i1mach(3)
      mdcon(3) = i1mach(1)
!
      imdcon = mdcon(k)
      end function imdcon
!STOPX
      logical function stopx(idummy)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      integer :: idummy
!
!     *****PURPOSE...
!     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION)
!     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT
!     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A
!     DYNAMIC STOPX.
!
!     *****ALGORITHM NOTES...
!     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED
!     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A
!     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT
!     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX.
!
!
      stopx = .false.

      end function stopx
!UFPARM
      subroutine ufparm

      end subroutine ufparm

      end module M_starpac_g

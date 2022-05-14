      module M_starpac
      use M_starpac_s
      use M_starpac_d
      implicit none
      private

      public :: ssifa,  dswap,  dtrdi,  dscal,   strco
      public :: dasum,  dnrm2,  dsifa,  sscal,   isamax
      public :: ddot,   daxpy,  sdot,   idamax,  sswap
      public :: dsidi,  dtrco,  ssidi,  strdi,   scopy
      public :: saxpy,  sasum,  dcopy,  snrm2

      public :: i1mach, d1mach, r1mach

      public :: r9lgic,  gamma,   dlngam,  erfc,    eprint
      public :: dgamr,   dlbeta,  dcsevl,  gamlim,  r9lgmc
      public :: alngam,  xerprt,  inits,   csevl,   xgetua
      public :: e9rint,  xerror,  d9gmit,  xsetf,   r9gmit
      public :: dlnrel,  gami,    d9lgic,  derf,    xerrwv
      public :: d9lgit,  seterr,  fdump,   r9lgit,  xerctl
      public :: s88fmt,  erf,     xerabt,  alnrel,  initds
      public :: dbetai,  i8save,  dgamma,  derfc,   dgamlm
      public :: albeta,  xersav,  dgamit,  gamit,   dlgams
      public :: gamr,    xerclr,  betai,   j4save,  xgetf
      public :: d9lgmc,  algams,  dgami

      contains

!SSIFA
      SUBROUTINE SSIFA(A,LDA,N,KPVT,INFO)
!
!     LATEST REVISION  -  JANUARY 24, 1990
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER INFO,LDA,N
!
!  ARRAY ARGUMENTS
      REAL A(LDA,*)
      INTEGER KPVT(*)
!
!  LOCAL SCALARS
      REAL ABSAKK,AK,AKM1,ALPHA,BK,BKM1,COLMAX,DENOM,MULK,MULKM1,ROWMAX,
     +   T
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP
      LOGICAL SWAP
!
!  EXTERNAL FUNCTIONS
      INTEGER ISAMAX
      EXTERNAL ISAMAX
!
!  EXTERNAL SUBROUTINES
      EXTERNAL SAXPY,SSWAP
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AMAX1,SQRT
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
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
!
      INFO = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      K = N
   10 CONTINUE
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0E0) INFO = 1
!     ......EXIT
            GO TO 200
   20    CONTINUE
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         KM1 = K - 1
         ABSAKK = ABS(A(K,K))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         IMAX = ISAMAX(K-1,A(1,K),1)
         COLMAX = ABS(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = AMAX1(ROWMAX,ABS(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ISAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = AMAX1(ROWMAX,ABS(A(JMAX,IMAX)))
   50       CONTINUE
            IF (ABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (AMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
!
!           1 X 1 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 120
!
!              PERFORM AN INTERCHANGE.
!
               CALL SSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL SAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 160
!
!              PERFORM AN INTERCHANGE.
!
               CALL SSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0E0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL SAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL SAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
!DSWAP
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*),DY(*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION DTEMP1,DTEMP2,DTEMP3
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
!
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
!
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
!DTRDI
      SUBROUTINE DTRDI(T,LDT,N,DET,JOB,INFO)
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
      INTEGER INFO,JOB,LDT,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DET(2),T(LDT,*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION TEMP,TEN
      INTEGER I,J,K,KB,KM1,KP1
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DAXPY,DSCAL
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DABS,MOD
!
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 180
!
!        COMPUTE DETERMINANT
!
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
!           ...EXIT
               IF (DET(1) .EQ. 0.0D0) GO TO 60
   10          IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
                  DET(1) = TEN*DET(1)
                  DET(2) = DET(2) - 1.0D0
               GO TO 10
   20          CONTINUE
   30          IF (DABS(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/TEN
                  DET(2) = DET(2) + 1.0D0
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  DO 100 K = 1, N
                     INFO = K
!              ......EXIT
                     IF (T(K,K) .EQ. 0.0D0) GO TO 110
                     T(K,K) = 1.0D0/T(K,K)
                     TEMP = -T(K,K)
                     CALL DSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = 0.0D0
                        CALL DAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
!     ............EXIT
                  IF (T(K,K) .EQ. 0.0D0) GO TO 180
                  T(K,K) = 1.0D0/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL DSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = 0.0D0
                     CALL DAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
      END
!DSCAL
      SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!     REPLACE DOUBLE PRECISION DX BY DOUBLE PRECISION DA*DX.
!     FOR I = 0 TO N-1, REPLACE DX(1+I*INCX) WITH  DA * DX(1+I*INCX)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION DA
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*)
!
!  LOCAL SCALARS
      INTEGER I,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
!STRCO
      SUBROUTINE STRCO(T,LDT,N,RCOND,Z,JOB)
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
      REAL RCOND
      INTEGER
     +   JOB,LDT,N

!...ARRAY ARGUMENTS
      REAL T(LDT,*),Z(*)

!...LOCAL SCALARS
      REAL EK,S,SM,TNORM,W,WK,WKM,YNORM
      INTEGER
     +   I1,J,J1,J2,K,KK,L
      LOGICAL
     +   LOWER

!...EXTERNAL FUNCTIONS
      REAL SASUM
      EXTERNAL
     +   SASUM

!...EXTERNAL SUBROUTINES
      EXTERNAL
     +   SAXPY,SSCAL

!...INTRINSIC FUNCTIONS
      INTRINSIC
     +   ABS,AMAX1,SIGN


!***FIRST EXECUTABLE STATEMENT  STRCO


      LOWER = JOB .EQ. 0

!     COMPUTE 1-NORM OF T

      TNORM = 0.0E0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = AMAX1(TNORM,SASUM(L,T(I1,J),1))
   10 CONTINUE

!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.

!     SOLVE TRANS(T)*Y = E

      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(T(K,K))) GO TO 30
            S = ABS(T(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (T(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + ABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)

      YNORM = 1.0E0

!     SOLVE T*Z = Y

      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (ABS(Z(K)) .LE. ABS(T(K,K))) GO TO 110
            S = ABS(T(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0E0) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL SAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM

      IF (TNORM .NE. 0.0E0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
!DASUM
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
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
      INTEGER
     +   INCX,N

!...ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   DX(*)

!...LOCAL SCALARS
      INTEGER
     +   I,M,MP1,NS

!...INTRINSIC FUNCTIONS
      INTRINSIC
     +   DABS,MOD


!***FIRST EXECUTABLE STATEMENT  DASUM


      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
!DNRM2
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
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
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION CUTHI,CUTLO,HITEST,ONE,SUM,XMAX,ZERO
      INTEGER I,J,NEXT,NN
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DABS,DSQRT,FLOAT
!
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
      DATA   ZERO, ONE /0.0D0, 1.0D0/
!
      XMAX = ZERO
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
!
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
!                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
!
!                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
!
!                                PREPARE FOR PHASE 4.
!
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
!
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 SUM = (SUM * XMAX) * XMAX
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 HITEST = CUTHI/FLOAT( N )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
!
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
!DSIFA
      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)
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
      INTEGER INFO,LDA,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION A(LDA,*)
      INTEGER KPVT(*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION ABSAKK,AK,AKM1,ALPHA,BK,BKM1,COLMAX,DENOM,MULK,
     +   MULKM1,ROWMAX,T
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP
      LOGICAL SWAP
!
!  EXTERNAL FUNCTIONS
      INTEGER IDAMAX
      EXTERNAL IDAMAX
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DAXPY,DSWAP
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DABS,DMAX1,DSQRT
!
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
!
      INFO = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      K = N
   10 CONTINUE
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0D0) INFO = 1
!     ......EXIT
            GO TO 200
   20    CONTINUE
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         KM1 = K - 1
         ABSAKK = DABS(A(K,K))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         IMAX = IDAMAX(K-1,A(1,K),1)
         COLMAX = DABS(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))
   50       CONTINUE
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
!
!           1 X 1 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 120
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 160
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0D0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
!SSCAL
      SUBROUTINE SSCAL(N,SA,SX,INCX)
!
!     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
!     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL SA
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      REAL SX(*)
!
!  LOCAL SCALARS
      INTEGER I,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
!ISAMAX
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
!
!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
!     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX))
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      REAL SX(*)
!
!  LOCAL SCALARS
      REAL SMAX,XMAG
      INTEGER I,II,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS
!
      ISAMAX = 0
      IF(N.LE.0) RETURN
      ISAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          ISAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END
!DDOT
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*),DY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     *    DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
!
!         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
!DAXPY
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
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
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*),DY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
!SDOT
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      REAL SX(*),SY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
!
!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     *    SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
      RETURN
!
!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
!IDAMAX
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF DOUBLE PRECISION DX.
!     IDAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(DX(1-INCX+I*INCX))
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION DMAX,XMAG
      INTEGER I,II,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DABS
!
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
!SSWAP
      SUBROUTINE SSWAP (N,SX,INCX,SY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      REAL SX(*),SY(*)
!
!  LOCAL SCALARS
      REAL STEMP1,STEMP2,STEMP3
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP1 = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
!
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        STEMP1 = SX(I)
        STEMP2 = SX(I+1)
        STEMP3 = SX(I+2)
        SX(I) = SY(I)
        SX(I+1) = SY(I+1)
        SX(I+2) = SY(I+2)
        SY(I) = STEMP1
        SY(I+1) = STEMP2
        SY(I+2) = STEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
!
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
      NS = N*INCX
        DO 70 I=1,NS,INCX
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   70   CONTINUE
      RETURN
      END
!DSIDI
      SUBROUTINE DSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)
!
!     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
!     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM
!     DSIFA.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER JOB,LDA,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
      INTEGER INERT(3),KPVT(*)
!
!  LOCAL SCALARS
      DOUBLE PRECISION AK,AKKP1,AKP1,D,T,TEMP,TEN
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NODET,NOERT,NOINV
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DAXPY,DCOPY,DSWAP
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DABS,IABS,MOD
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
      TEN = 10.0D0
!
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
!
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0D0
            DET(2) = 0.0D0
   20    CONTINUE
         T = 0.0D0
         DO 130 K = 1, N
            D = A(K,K)
!
!           CHECK IF 1 BY 1
!
            IF (KPVT(K) .GT. 0) GO TO 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               IF (T .NE. 0.0D0) GO TO 30
                  T = DABS(A(K,K+1))
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0D0
   40          CONTINUE
   50       CONTINUE
!
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
!
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0D0) GO TO 110
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0D0
                  GO TO 70
   80             CONTINUE
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0D0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
!
!              1 BY 1
!
               A(K,K) = 1.0D0/A(K,K)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
!
!              2 BY 2
!
               T = DABS(A(K,K+1))
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0D0)
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL DSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
!DTRCO
      SUBROUTINE DTRCO(T,LDT,N,RCOND,Z,JOB)
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
      DOUBLE PRECISION
     +   RCOND
      INTEGER
     +   JOB,LDT,N

!...ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   T(LDT,*),Z(*)

!...LOCAL SCALARS
      DOUBLE PRECISION
     +   EK,S,SM,TNORM,W,WK,WKM,YNORM
      INTEGER
     +   I1,J,J1,J2,K,KK,L
      LOGICAL
     +   LOWER

!...EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DASUM
      EXTERNAL
     +   DASUM

!...EXTERNAL SUBROUTINES
      EXTERNAL
     +   DAXPY,DSCAL

!...INTRINSIC FUNCTIONS
      INTRINSIC
     +   DABS,DMAX1,DSIGN


!***FIRST EXECUTABLE STATEMENT  DTRCO


      LOWER = JOB .EQ. 0

!     COMPUTE 1-NORM OF T

      TNORM = 0.0D0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = DMAX1(TNORM,DASUM(L,T(I1,J),1))
   10 CONTINUE

!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.

!     SOLVE TRANS(T)*Y = E

      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(T(K,K))) GO TO 30
            S = DABS(T(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (T(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + DABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)

      YNORM = 1.0D0

!     SOLVE T*Z = Y

      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (DABS(Z(K)) .LE. DABS(T(K,K))) GO TO 110
            S = DABS(T(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0D0) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL DAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM

      IF (TNORM .NE. 0.0D0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
!SSIDI
      SUBROUTINE SSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)
!
!     LATEST REVISION  -  JANUARY 24, 1990  (JRD)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER JOB,LDA,N
!
!  ARRAY ARGUMENTS
      REAL A(LDA,*),DET(2),WORK(*)
      INTEGER INERT(3),KPVT(*)
!
!  LOCAL SCALARS
      REAL AK,AKKP1,AKP1,D,T,TEMP,TEN
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NODET,NOERT,NOINV
!
!  EXTERNAL FUNCTIONS
      REAL SDOT
      EXTERNAL SDOT
!
!  EXTERNAL SUBROUTINES
      EXTERNAL SAXPY,SCOPY,SSWAP
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,IABS,MOD
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
!     BLAS SAXPY,SCOPY,SDOT,SSWAP
!     FORTRAN ABS,IABS,MOD
!
!
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
!
      TEN = 10.0E0
!
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0E0
            DET(2) = 0.0E0
   20    CONTINUE
         T = 0.0E0
         DO 130 K = 1, N
            D = A(K,K)
!
!           CHECK IF 1 BY 1
!
            IF (KPVT(K) .GT. 0) GO TO 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               IF (T .NE. 0.0E0) GO TO 30
                  T = ABS(A(K,K+1))
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0E0
   40          CONTINUE
   50       CONTINUE
!
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0E0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0E0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0E0) INERT(3) = INERT(3) + 1
   60       CONTINUE
!
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0E0) GO TO 110
   70             IF (ABS(DET(1)) .GE. 1.0E0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0E0
                  GO TO 70
   80             CONTINUE
   90             IF (ABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0E0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
!
!              1 BY 1
!
               A(K,K) = 1.0E0/A(K,K)
               IF (KM1 .LT. 1) GO TO 170
                  CALL SCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = SDOT(J,A(1,J),1,WORK,1)
                     CALL SAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + SDOT(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
!
!              2 BY 2
!
               T = ABS(A(K,K+1))
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0E0)
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL SCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = SDOT(J,A(1,J),1,WORK,1)
                     CALL SAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1) + SDOT(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + SDOT(KM1,A(1,K),1,A(1,K+1),1)
                  CALL SCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = SDOT(J,A(1,J),1,WORK,1)
                     CALL SAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K) + SDOT(KM1,WORK,1,A(1,K),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL SSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
!STRDI
      SUBROUTINE STRDI(T,LDT,N,DET,JOB,INFO)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER INFO,JOB,LDT,N
!
!  ARRAY ARGUMENTS
      REAL DET(2),T(LDT,*)
!
!  LOCAL SCALARS
      REAL TEMP,TEN
      INTEGER I,J,K,KB,KM1,KP1
!
!  EXTERNAL SUBROUTINES
      EXTERNAL SAXPY,SSCAL
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MOD
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
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = 1.0E0
            DET(2) = 0.0E0
            TEN = 10.0E0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
!           ...EXIT
               IF (DET(1) .EQ. 0.0E0) GO TO 60
   10          IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
                  DET(1) = TEN*DET(1)
                  DET(2) = DET(2) - 1.0E0
               GO TO 10
   20          CONTINUE
   30          IF (ABS(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/TEN
                  DET(2) = DET(2) + 1.0E0
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  DO 100 K = 1, N
                     INFO = K
!              ......EXIT
                     IF (T(K,K) .EQ. 0.0E0) GO TO 110
                     T(K,K) = 1.0E0/T(K,K)
                     TEMP = -T(K,K)
                     CALL SSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = 0.0E0
                        CALL SAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
!     ............EXIT
                  IF (T(K,K) .EQ. 0.0E0) GO TO 180
                  T(K,K) = 1.0E0/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL SSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = 0.0E0
                     CALL SAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
      END
!SCOPY
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      REAL SX(*),SY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
!SAXPY
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
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
      REAL SA
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      REAL SX(*),SY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
!SASUM
      REAL FUNCTION SASUM(N,SX,INCX)
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
      INTEGER
     +   INCX,N

!...ARRAY ARGUMENTS
      REAL SX(*)

!...LOCAL SCALARS
      INTEGER
     +   I,M,MP1,NS

!...INTRINSIC FUNCTIONS
      INTRINSIC
     +   ABS,MOD


!***FIRST EXECUTABLE STATEMENT  SASUM


      SASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
!DCOPY
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
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
      INTEGER INCX,INCY,N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DX(*),DY(*)
!
!  LOCAL SCALARS
      INTEGER I,IX,IY,M,MP1,NS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END
!SNRM2
      REAL FUNCTION SNRM2 ( N, SX, INCX)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER INCX,N
!
!  ARRAY ARGUMENTS
      REAL SX(*)
!
!  LOCAL SCALARS
      REAL CUTHI,CUTLO,HITEST,ONE,SUM,XMAX,ZERO
      INTEGER I,J,NEXT,NN
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,FLOAT,SQRT
!
      DATA   ZERO, ONE /0.0E0, 1.0E0/
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
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
!
      XMAX = ZERO
      IF(N .GT. 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
!
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
!                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
   50 IF( SX(I) .EQ. ZERO) GO TO 200
      IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
!
!                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
!
!                                PREPARE FOR PHASE 4.
!
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 IF( ABS(SX(I)) .GT. CUTLO ) GO TO 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 IF( ABS(SX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
!
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 SUM = (SUM * XMAX) * XMAX
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 HITEST = CUTHI/FLOAT( N )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      DO 95 J =I,NN,INCX
      IF(ABS(SX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2
      SNRM2 = SQRT( SUM )
      GO TO 300
!
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END

!R9LGIC
      REAL FUNCTION R9LGIC (A, X, ALX)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG COMPLEMENTARY INCOMPLETE GAMMA FUNCTION FOR LARGE X
! AND FOR A .LE. X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL A,ALX,X
!
!  LOCAL SCALARS
      REAL EPS,FK,P,R,S,T,XMA,XPA
      INTEGER K
!
!  EXTERNAL FUNCTIONS
      REAL R1MACH
      EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG
!
      DATA EPS / 0.0 /
!
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
!
      XPA = X + 1.0 - A
      XMA = X - 1.0 - A
!
      R = 0.0
      P = 1.0
      S = P
      DO 10 K=1,200
        FK = K
        T = FK*(A-FK)*(1.0+R)
        R = -T/((XMA+2.0*FK)*(XPA+2.0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERROR (  'R9LGIC  NO CONVERGENCE IN 200 TERMS OF CONTINUED F
     1RACTION', 57, 1, 2)
!
 20   R9LGIC = A*ALX - X + LOG(S/XPA)
!
      RETURN
      END
!GAMMA
      REAL FUNCTION GAMMA (X)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL DXREL,PI,SINPIY,SQ2PIL,XMAX,XMIN,Y
      INTEGER I,N,NGCS
!
!  LOCAL ARRAYS
      REAL GCS(23)
!
!  EXTERNAL FUNCTIONS
      REAL CSEVL,R1MACH,R9LGMC
      INTEGER INITS
      EXTERNAL CSEVL,R1MACH,R9LGMC,INITS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL GAMLIM,XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP,FLOAT,LOG,SIN,SQRT
!
!
      DATA GCS   ( 1) / .0085711955 90989331E0/
      DATA GCS   ( 2) / .0044153813 24841007E0/
      DATA GCS   ( 3) / .0568504368 1599363E0/
      DATA GCS   ( 4) /-.0042198353 96418561E0/
      DATA GCS   ( 5) / .0013268081 81212460E0/
      DATA GCS   ( 6) /-.0001893024 529798880E0/
      DATA GCS   ( 7) / .0000360692 532744124E0/
      DATA GCS   ( 8) /-.0000060567 619044608E0/
      DATA GCS   ( 9) / .0000010558 295463022E0/
      DATA GCS   (10) /-.0000001811 967365542E0/
      DATA GCS   (11) / .0000000311 772496471E0/
      DATA GCS   (12) /-.0000000053 542196390E0/
      DATA GCS   (13) / .0000000009 193275519E0/
      DATA GCS   (14) /-.0000000001 577941280E0/
      DATA GCS   (15) / .0000000000 270798062E0/
      DATA GCS   (16) /-.0000000000 046468186E0/
      DATA GCS   (17) / .0000000000 007973350E0/
      DATA GCS   (18) /-.0000000000 001368078E0/
      DATA GCS   (19) / .0000000000 000234731E0/
      DATA GCS   (20) /-.0000000000 000040274E0/
      DATA GCS   (21) / .0000000000 000006910E0/
      DATA GCS   (22) /-.0000000000 000001185E0/
      DATA GCS   (23) / .0000000000 000000203E0/
!
      DATA PI /3.14159 26535 89793 24E0/
! SQ2PIL IS LOG (SQRT (2.*PI) )
      DATA SQ2PIL /0.91893 85332 04672 74E0/
      DATA NGCS, XMIN, XMAX, DXREL /0, 3*0.0 /
!
      IF (NGCS.NE.0) GO TO 10
!
! ---------------------------------------------------------------------
! INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF
! TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
! THAN MACHINE PRECISION.
!
      NGCS = INITS (GCS, 23, 0.1*R1MACH(3))
!
      CALL GAMLIM (XMIN, XMAX)
      DXREL = SQRT (R1MACH(4))
!
! ---------------------------------------------------------------------
! FINISH INITIALIZATION.  START EVALUATING GAMMA(X).
!
 10   Y = ABS(X)
      IF (Y.GT.10.0) GO TO 50
!
! COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND
! FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL.
!
      N = X
      IF (X.LT.0.) N = N - 1
      Y = X - FLOAT(N)
      N = N - 1
      GAMMA = 0.9375 + CSEVL(2.*Y-1., GCS, NGCS)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.
!
      N = -N
      IF (X.EQ.0.) CALL XERROR ('GAMMA   X IS 0', 14, 4, 2)
      IF (X.LT.0. .AND. X+FLOAT(N-2).EQ.0.) CALL XERROR (
     1  'GAMMA   X IS A NEGATIVE INTEGER', 31, 4, 2)
      IF (X.LT.(-0.5) .AND. ABS((X-AINT(X-0.5))/X).LT.DXREL) CALL
     1  XERROR (  'GAMMA   ANSWER LT HALF PRECISION BECAUSE X TOO NEAR N
     2EGATIVE INTEGER', 68, 1, 1)
!
      DO 20 I=1,N
        GAMMA = GAMMA / (X+FLOAT(I-1))
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.
!
 30   DO 40 I=1,N
        GAMMA = (Y+FLOAT(I))*GAMMA
 40   CONTINUE
      RETURN
!
! COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X.GT.XMAX) CALL XERROR ('GAMMA   X SO BIG GAMMA OVERFLOWS',
     1  32, 3, 2)
!
      GAMMA = 0.
      IF (X.LT.XMIN) CALL XERROR ('GAMMA   X SO SMALL GAMMA UNDERFLOWS',
     1  35, 2, 1)
      IF (X.LT.XMIN) RETURN
!
      GAMMA = EXP((Y-0.5)*LOG(Y) - Y + SQ2PIL + R9LGMC(Y) )
      IF (X.GT.0.) RETURN
!
      IF (ABS((X-AINT(X-0.5))/X).LT.DXREL) CALL XERROR (
     1  'GAMMA   ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER',
     2  61, 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY.EQ.0.) CALL XERROR (
     1  'GAMMA   X IS A NEGATIVE INTEGER', 31, 4, 2)
!
      GAMMA = -PI / (Y*SINPIY*GAMMA)
!
      RETURN
      END
!DLNGAM
      DOUBLE PRECISION FUNCTION DLNGAM (X)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION DXREL,PI,SINPIY,SQ2PIL,SQPI2L,XMAX,Y
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,D9LGMC,DGAMMA
      EXTERNAL D1MACH,D9LGMC,DGAMMA
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSQRT,INT,LOG,SIN
!
!
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
! SQ2PIL = LOG (SQRT(2*PI)),  SQPI2L = LOG(SQRT(PI/2))
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
!
      DATA XMAX, DXREL / 2*0.D0 /
!
      IF (XMAX.NE.0.D0) GO TO 10
      XMAX = D1MACH(2)/LOG(D1MACH(2))
      DXREL = DSQRT (D1MACH(4))
!
 10   Y = ABS (X)
      IF (Y.GT.10.D0) GO TO 20
!
! LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
!
      DLNGAM = LOG (ABS (DGAMMA(X)) )
      RETURN
!
! LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
!
 20   IF (Y.GT.XMAX) CALL XERROR (
     1  'DLNGAM  ABS(X) SO BIG DLNGAM OVERFLOWS', 39, 2, 2)
!
      IF (X.GT.0.D0) THEN
         DLNGAM = SQ2PIL + (X-0.5D0)*LOG(X) - X + D9LGMC(Y)
         RETURN
      END IF
!
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY.EQ.0.D0) CALL XERROR (
     1  'DLNGAM  X IS A NEGATIVE INTEGER', 31, 3, 2)
!
      IF (ABS ((X-INT(X-0.5D0))/X).LT.DXREL) CALL XERROR (
     1    'DLNGAM  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE
     2INTEGER', 68, 1, 1)
!
      DLNGAM = SQPI2L + (X-0.5D0)*LOG(Y) - X - LOG(SINPIY) - D9LGMC(Y)
      RETURN
!
      END
!ERFC
      REAL FUNCTION ERFC (X)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL ETA,SQEPS,SQRTPI,XMAX,XSML,Y
      INTEGER NTERC2,NTERF,NTERFC
!
!  LOCAL ARRAYS
      REAL ERC2CS(23),ERFCCS(24),ERFCS(13)
!
!  EXTERNAL FUNCTIONS
      REAL CSEVL,R1MACH
      INTEGER INITS
      EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,EXP,LOG,SQRT
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000D+00
!                                        WITH WEIGHTED ERROR   7.10E-18
!                                         LOG WEIGHTED ERROR  17.15
!                               SIGNIFICANT FIGURES REQUIRED  16.31
!                                    DECIMAL PLACES REQUIRED  17.71
!
      DATA ERF CS( 1) /   -.0490461212 34691808E0 /
      DATA ERF CS( 2) /   -.1422612051 0371364E0 /
      DATA ERF CS( 3) /    .0100355821 87599796E0 /
      DATA ERF CS( 4) /   -.0005768764 69976748E0 /
      DATA ERF CS( 5) /    .0000274199 31252196E0 /
      DATA ERF CS( 6) /   -.0000011043 17550734E0 /
      DATA ERF CS( 7) /    .0000000384 88755420E0 /
      DATA ERF CS( 8) /   -.0000000011 80858253E0 /
      DATA ERF CS( 9) /    .0000000000 32334215E0 /
      DATA ERF CS(10) /   -.0000000000 00799101E0 /
      DATA ERF CS(11) /    .0000000000 00017990E0 /
      DATA ERF CS(12) /   -.0000000000 00000371E0 /
      DATA ERF CS(13) /    .0000000000 00000007E0 /
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
      DATA ERC2CS( 1) /   -.0696013466 02309501E0 /
      DATA ERC2CS( 2) /   -.0411013393 62620893E0 /
      DATA ERC2CS( 3) /    .0039144958 66689626E0 /
      DATA ERC2CS( 4) /   -.0004906395 65054897E0 /
      DATA ERC2CS( 5) /    .0000715747 90013770E0 /
      DATA ERC2CS( 6) /   -.0000115307 16341312E0 /
      DATA ERC2CS( 7) /    .0000019946 70590201E0 /
      DATA ERC2CS( 8) /   -.0000003642 66647159E0 /
      DATA ERC2CS( 9) /    .0000000694 43726100E0 /
      DATA ERC2CS(10) /   -.0000000137 12209021E0 /
      DATA ERC2CS(11) /    .0000000027 88389661E0 /
      DATA ERC2CS(12) /   -.0000000005 81416472E0 /
      DATA ERC2CS(13) /    .0000000001 23892049E0 /
      DATA ERC2CS(14) /   -.0000000000 26906391E0 /
      DATA ERC2CS(15) /    .0000000000 05942614E0 /
      DATA ERC2CS(16) /   -.0000000000 01332386E0 /
      DATA ERC2CS(17) /    .0000000000 00302804E0 /
      DATA ERC2CS(18) /   -.0000000000 00069666E0 /
      DATA ERC2CS(19) /    .0000000000 00016208E0 /
      DATA ERC2CS(20) /   -.0000000000 00003809E0 /
      DATA ERC2CS(21) /    .0000000000 00000904E0 /
      DATA ERC2CS(22) /   -.0000000000 00000216E0 /
      DATA ERC2CS(23) /    .0000000000 00000052E0 /
!
!                                    DECIMAL PLACES REQUIRED  17.01
!
      DATA ERFCCS( 1) /   0.0715179310 202925E0 /
      DATA ERFCCS( 2) /   -.0265324343 37606719E0 /
      DATA ERFCCS( 3) /    .0017111539 77920853E0 /
      DATA ERFCCS( 4) /   -.0001637516 63458512E0 /
      DATA ERFCCS( 5) /    .0000198712 93500549E0 /
      DATA ERFCCS( 6) /   -.0000028437 12412769E0 /
      DATA ERFCCS( 7) /    .0000004606 16130901E0 /
      DATA ERFCCS( 8) /   -.0000000822 77530261E0 /
      DATA ERFCCS( 9) /    .0000000159 21418724E0 /
      DATA ERFCCS(10) /   -.0000000032 95071356E0 /
      DATA ERFCCS(11) /    .0000000007 22343973E0 /
      DATA ERFCCS(12) /   -.0000000001 66485584E0 /
      DATA ERFCCS(13) /    .0000000000 40103931E0 /
      DATA ERFCCS(14) /   -.0000000000 10048164E0 /
      DATA ERFCCS(15) /    .0000000000 02608272E0 /
      DATA ERFCCS(16) /   -.0000000000 00699105E0 /
      DATA ERFCCS(17) /    .0000000000 00192946E0 /
      DATA ERFCCS(18) /   -.0000000000 00054704E0 /
      DATA ERFCCS(19) /    .0000000000 00015901E0 /
      DATA ERFCCS(20) /   -.0000000000 00004729E0 /
      DATA ERFCCS(21) /    .0000000000 00001432E0 /
      DATA ERFCCS(22) /   -.0000000000 00000439E0 /
      DATA ERFCCS(23) /    .0000000000 00000138E0 /
      DATA ERFCCS(24) /   -.0000000000 00000048E0 /
!
      DATA SQRTPI /1.772453850 9055160E0/
      DATA NTERF, NTERFC, NTERC2, XSML, XMAX, SQEPS /3*0, 3*0./
!
      IF (NTERF.NE.0) GO TO 10
      ETA = 0.1*R1MACH(3)
      NTERF = INITS (ERFCS, 13, ETA)
      NTERFC = INITS (ERFCCS, 24, ETA)
      NTERC2 = INITS (ERC2CS, 23, ETA)
!
      XSML = -SQRT (-LOG(SQRTPI*R1MACH(3)))
      XMAX = SQRT (-LOG(SQRTPI*R1MACH(1)))
      XMAX = XMAX - 0.5*LOG(XMAX)/XMAX - 0.01
      SQEPS = SQRT (2.0*R1MACH(3))
!
 10   IF (X.GT.XSML) GO TO 20
!
! ERFC(X) = 1.0 - ERF(X) FOR X .LT. XSML
!
      ERFC = 2.0
      RETURN
!
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0) GO TO 30
!
! ERFC(X) = 1.0 - ERF(X) FOR -1. .LE. X .LE. 1.
!
      IF (Y.LT.SQEPS) THEN
         ERFC = 1.0 - 2.0*X/SQRTPI
      ELSE
         ERFC = 1.0 - X*(1.0 + CSEVL (2.*X*X-1., ERFCS, NTERF) )
      END IF
!
      RETURN
!
! ERFC(X) = 1.0 - ERF(X) FOR 1. .LT. ABS(X) .LE. XMAX
!
 30   Y = Y*Y
      IF (Y.LE.4.) THEN
         ERFC = EXP(-Y)/ABS(X) *
     +         (0.5 + CSEVL ((8.0/Y-5.0)/3.0, ERC2CS, NTERC2) )
      ELSE
         ERFC = EXP(-Y)/ABS(X) *
     +          (0.5 + CSEVL (8.0/Y-1.0, ERFCCS, NTERFC) )
      END IF
      IF (X.LT.0.0) ERFC = 2.0 - ERFC
      RETURN
!
 40   CALL XERROR ('ERFC    X SO BIG ERFC UNDERFLOWS', 32, 1, 1)
      ERFC = 0.0
      RETURN
!
      END
!EPRINT
      SUBROUTINE EPRINT
!
!  THIS SUBROUTINE PRINTS THE LAST ERROR MESSAGE, IF ANY.
!
!
!  VARIABLE DECLARATIONS
!
!  LOCAL ARRAYS
      CHARACTER MESSG(1)*4
!
!  EXTERNAL SUBROUTINES
      EXTERNAL E9RINT
!
!
      CALL E9RINT(MESSG,1,1,.FALSE.)
      RETURN
!
      END
!DGAMR
      DOUBLE PRECISION FUNCTION DGAMR (X)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! THIS ROUTINE, NOT DGAMMA(X), SHOULD BE THE FUNDAMENTAL ONE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION ALNGX,SGNGX
      INTEGER IROLD
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION DGAMMA
      EXTERNAL DGAMMA
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DLGAMS,XERCLR,XGETF,XSETF
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,EXP,INT
!
!
      DGAMR = 0.0D0
      IF (X.LE.0.0D0 .AND. INT(X).EQ.X) RETURN
!
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0D0) GO TO 10
      DGAMR = 1.0D0/DGAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
!
 10   CALL DLGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      DGAMR = SGNGX * EXP(-ALNGX)
      RETURN
!
      END
!DLBETA
      DOUBLE PRECISION FUNCTION DLBETA (A, B)
! JULY 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION A,B
!
!  LOCAL SCALARS
      DOUBLE PRECISION CORR,P,Q,SQ2PIL
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D9LGMC,DGAMMA,DLNGAM,DLNREL
      EXTERNAL D9LGMC,DGAMMA,DLNGAM,DLNREL
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC LOG,MAX,MIN
!
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
!
      P = MIN (A, B)
      Q = MAX (A, B)
!
      IF (P.LE.0.D0) CALL XERROR (
     1  'DLBETA  BOTH ARGUMENTS MUST BE GT ZERO', 38, 1, 2)
!
      IF (P.GE.10.D0) GO TO 30
      IF (Q.GE.10.D0) GO TO 20
!
! P AND Q ARE SMALL.
!
      DLBETA = LOG (DGAMMA(P) * (DGAMMA(Q)/DGAMMA(P+Q)) )
      RETURN
!
! P IS SMALL, BUT Q IS BIG.
!
 20   CORR = D9LGMC(Q) - D9LGMC(P+Q)
      DLBETA = DLNGAM(P) + CORR + P - P*LOG(P+Q)
     1  + (Q-0.5D0)*DLNREL(-P/(P+Q))
      RETURN
!
! P AND Q ARE BIG.
!
 30   CORR = D9LGMC(P) + D9LGMC(Q) - D9LGMC(P+Q)
      DLBETA = -0.5D0*LOG(Q) + SQ2PIL + CORR + (P-0.5D0)*LOG(P/(P+Q))
     1  + Q*DLNREL(-P/(P+Q))
      RETURN
!
      END
!DCSEVL
      DOUBLE PRECISION FUNCTION DCSEVL (X, A, N)
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
      DOUBLE PRECISION X
      INTEGER N
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION A(N)
!
!  LOCAL SCALARS
      DOUBLE PRECISION B0,B1,B2,TWOX
      INTEGER I,NI1
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!
      IF (N.LT.1) CALL XERROR ('DCSEVL  NUMBER OF TERMS LE 0', 28, 2,2)
      IF (N.GT.1000) CALL XERROR ('DCSEVL  NUMBER OF TERMS GT 1000',
     1  31, 3, 2)
      IF (X.LT.(-1.D0) .OR. X.GT.1.D0) CALL XERROR (
     1  'DCSEVL  X OUTSIDE (-1,+1)', 25, 1, 1)
!
      TWOX = 2.0D0*X
      B0 = 0.0D0
      B1 = 0.0D0
      B2 = 0.0D0
      DO 10 I=1,N
        B2 = B1
        B1 = B0
        NI1 = N-I+1
        B0 = TWOX*B1 - B2 + A(NI1)
 10   CONTINUE
!
      DCSEVL = 0.5D0 * (B0-B2)
!
      RETURN
      END
!GAMLIM
      SUBROUTINE GAMLIM (XMIN, XMAX)
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
      REAL XMAX,XMIN
!
!  LOCAL SCALARS
      REAL ALNBIG,ALNSML,XLN,XOLD
      INTEGER I
!
!  EXTERNAL FUNCTIONS
      REAL R1MACH
      EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,MAX
!
      ALNSML = LOG(R1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5)*XLN - XMIN - 0.2258 + ALNSML)
     1    / (XMIN*XLN + 0.5)
        IF (ABS(XMIN-XOLD).LT.0.005) GO TO 20
 10   CONTINUE
      CALL XERROR ('GAMLIM  UNABLE TO FIND XMIN', 27, 1, 2)
!
 20   XMIN = -XMIN + 0.01
!
      ALNBIG = LOG(R1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5)*XLN - XMAX + 0.9189 - ALNBIG)
     1    / (XMAX*XLN - 0.5)
        IF (ABS(XMAX-XOLD).LT.0.005) GO TO 40
 30   CONTINUE
      CALL XERROR ('GAMLIM  UNABLE TO FIND XMAX', 27, 2, 2)
!
 40   XMAX = XMAX - 0.01
      XMIN = MAX (XMIN, -XMAX+1.)
!
      RETURN
      END
!R9LGMC
      REAL FUNCTION R9LGMC (X)
! AUGUST 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG GAMMA CORRECTION FACTOR FOR X .GE. 10.0 SO THAT
!  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL XBIG,XMAX
      INTEGER NALGM
!
!  LOCAL ARRAYS
      REAL ALGMCS(6)
!
!  EXTERNAL FUNCTIONS
      REAL CSEVL,R1MACH
      INTEGER INITS
      EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC EXP,LOG,MIN,SQRT
!
!
! SERIES FOR ALGM       ON THE INTERVAL  0.          TO  1.00000D-02
!                                        WITH WEIGHTED ERROR   3.40E-16
!                                         LOG WEIGHTED ERROR  15.47
!                               SIGNIFICANT FIGURES REQUIRED  14.39
!                                    DECIMAL PLACES REQUIRED  15.86
!
      DATA ALGMCS( 1) /    .1666389480 45186E0 /
      DATA ALGMCS( 2) /   -.0000138494 817606E0 /
      DATA ALGMCS( 3) /    .0000000098 108256E0 /
      DATA ALGMCS( 4) /   -.0000000000 180912E0 /
      DATA ALGMCS( 5) /    .0000000000 000622E0 /
      DATA ALGMCS( 6) /   -.0000000000 000003E0 /
!
      DATA NALGM, XBIG, XMAX / 0, 2*0.0 /
!
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITS (ALGMCS, 6, R1MACH(3))
      XBIG = 1.0/SQRT(R1MACH(3))
      XMAX = EXP (MIN(LOG(R1MACH(2)/12.0), -LOG(12.0*R1MACH(1))) )
!
 10   IF (X.LT.10.0) CALL XERROR ('R9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      R9LGMC = 1.0/(12.0*X)
      IF (X.LT.XBIG) R9LGMC = CSEVL (2.0*(10./X)**2-1., ALGMCS, NALGM)/X
      RETURN
!
 20   R9LGMC = 0.0
      CALL XERROR ('R9LGMC  X SO BIG R9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
!
      END
!ALNGAM
      REAL FUNCTION ALNGAM (X)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL DXREL,PI,SINPIY,SQ2PIL,SQPI2L,XMAX,Y
!
!  EXTERNAL FUNCTIONS
      REAL GAMMA,R1MACH,R9LGMC
      EXTERNAL GAMMA,R1MACH,R9LGMC
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,LOG,SIN,SQRT
!
      DATA SQ2PIL / 0.9189385332 0467274E0/
! SQ2PIL = LOG(SQRT(2.*PI)),  SQPI2L = LOG (SQRT(PI/2.))
      DATA SQPI2L / 0.2257913526 4472743E0/
      DATA PI     / 3.1415926535 8979324E0/
!
      DATA XMAX, DXREL / 0., 0. /
!
      IF (XMAX.NE.0.) GO TO 10
      XMAX = R1MACH(2)/LOG(R1MACH(2))
      DXREL = SQRT (R1MACH(4))
!
 10   Y = ABS(X)
      IF (Y.GT.10.0) GO TO 20
!
! LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
!
      ALNGAM = LOG (ABS (GAMMA(X)))
      RETURN
!
! LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
!
 20   IF (Y.GT.XMAX) CALL XERROR (
     1  'ALNGAM  ABS(X) SO BIG ALNGAM OVERFLOWS', 38, 2, 2)
!
      IF (X.GT.0.0) THEN
         ALNGAM = SQ2PIL + (X-0.5)*LOG(X) - X + R9LGMC(Y)
         RETURN
      END IF
!
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY.EQ.0.) CALL XERROR ('ALNGAM  X IS A NEGATIVE INTEGER',
     1  31, 3, 2)
!
      IF (ABS((X-AINT(X-0.5))/X).LT.DXREL) CALL XERROR (
     1    'ALNGAM  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE
     2INTEGER', 68, 1, 1)
!
      ALNGAM = SQPI2L + (X-0.5)*LOG(Y) - X - LOG(SINPIY) - R9LGMC(Y)
      RETURN
!
      END
!XERPRT
      SUBROUTINE XERPRT(MESSG,NMESSG)
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
      INTEGER NMESSG
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NMESSG)*4
!
!  LOCAL SCALARS
      INTEGER I,IUNIT,KUNIT,NCHAR,NCHARL,NCHLST,NCHREM,NFIELD,NLINES,
     +   NUNIT,NWORD,NWORD1,NWORD2
      CHARACTER LA*1,LBLANK*1,LCOM*1
!
!  LOCAL ARRAYS
      INTEGER LUN(5)
      CHARACTER F(10)*1,G(14)*1
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH
      EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL S88FMT,XGETUA
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
      DATA F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),F(9),F(10)
     1   / '(' ,'1' ,'X' ,',' ,' ' ,' ' ,'A' ,' ' ,' ' ,')' /
      DATA G(1),G(2),G(3),G(4),G(5),G(6),G(7),G(8),G(9),G(10)
     1   / '(' ,'1' ,'X' ,' ' ,' ' ,' ' ,' ' ,' ' ,' ' ,' ' /
      DATA G(11),G(12),G(13),G(14)
     1   / ' '  ,' '  ,' '  ,')'  /
      DATA LA/'A'/,LCOM/','/,LBLANK/' '/
!     PREPARE FORMAT FOR WHOLE LINES
      NCHAR = I1MACH(6)
      NFIELD = 72/NCHAR
      CALL S88FMT(2,NFIELD,F(5))
      CALL S88FMT(2,NCHAR,F(8))
!     PREPARE FORMAT FOR LAST, PARTIAL LINE, IF NEEDED
      NCHARL = NFIELD*NCHAR
      NLINES = NMESSG/NCHARL
      NWORD  = NLINES*NFIELD
      NCHREM = NMESSG - NLINES*NCHARL
      IF (NCHREM.LE.0) GO TO 40
         DO 10 I=4,13
10          G(I) = LBLANK
         NFIELD = NCHREM/NCHAR
         IF (NFIELD.LE.0) GO TO 20
!        PREPARE WHOLE WORD FIELDS
            G(4) = LCOM
            CALL S88FMT(2,NFIELD,G(5))
            G(7) = LA
            CALL S88FMT(2,NCHAR,G(8))
20       CONTINUE
         NCHLST = MOD(NCHREM,NCHAR)
         IF (NCHLST.LE.0) GO TO 30
!        PREPARE PARTIAL WORD FIELD
            G(10) = LCOM
            G(11) = LA
            CALL S88FMT(2,NCHLST,G(12))
30       CONTINUE
40    CONTINUE
!     PRINT THE MESSAGE
      NWORD1 = NWORD+1
      NWORD2 = (NMESSG+NCHAR-1)/NCHAR
      CALL XGETUA(LUN,NUNIT)
      DO 50 KUNIT = 1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         IF (NWORD.GT.0) WRITE (IUNIT,F) (MESSG(I),I=1,NWORD)
         IF (NCHREM.GT.0) WRITE (IUNIT,G) (MESSG(I),I=NWORD1,NWORD2)
50    CONTINUE
      RETURN
      END
!INITS
      INTEGER FUNCTION INITS (OS, NOS, ETA)
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
      REAL ETA
      INTEGER NOS
!
!  ARRAY ARGUMENTS
      REAL OS(NOS)
!
!  LOCAL SCALARS
      REAL ERR
      INTEGER I,II
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS
!
!
      IF (NOS.LT.1) CALL XERROR (
     1  'INITS   NUMBER OF COEFFICIENTS LT 1', 35, 2, 2)
!
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(OS(I))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
!
 20   IF (I.EQ.NOS) CALL XERROR ('INITS   ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITS = I
!
      RETURN
      END
!CSEVL
      REAL FUNCTION CSEVL (X, CS, N)
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
      REAL X
      INTEGER N
!
!  ARRAY ARGUMENTS
      REAL CS(N)
!
!  LOCAL SCALARS
      REAL B0,B1,B2,TWOX
      INTEGER I,NI
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!
      IF (N.LT.1) CALL XERROR ('CSEVL   NUMBER OF TERMS LE 0', 28, 2,2)
      IF (N.GT.1000) CALL XERROR ('CSEVL   NUMBER OF TERMS GT 1000',
     1  31, 3, 2)
      IF (X.LT.(-1.0) .OR. X.GT.1.0) CALL XERROR (
     1  'CSEVL   X OUTSIDE (-1,+1)', 25, 1, 1)
!
      B0 = 0.0
      B1 = 0.0
      B2 = 0.0
      TWOX = 2.0*X
      DO 10 I=1,N
        B2 = B1
        B1 = B0
        NI = N + 1 - I
        B0 = TWOX*B1 - B2 + CS(NI)
 10   CONTINUE
!
      CSEVL = 0.5 * (B0-B2)
!
      RETURN
      END
!XGETUA
      SUBROUTINE XGETUA(IUNIT,N)
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
      INTEGER N
!
!  ARRAY ARGUMENTS
      INTEGER IUNIT(5)
!
!  LOCAL SCALARS
      INTEGER I,INDEX
!
!  EXTERNAL FUNCTIONS
      INTEGER J4SAVE
      EXTERNAL J4SAVE
!
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNIT(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
!E9RINT
      SUBROUTINE E9RINT(MESSG,NW,NERR,SAVE)
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
      INTEGER NERR,NW
      LOGICAL SAVE
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NW)*4
!
!  LOCAL SCALARS
      INTEGER I,IWUNIT,NERRP,NWP
      CHARACTER CCPLUS*4
!
!  LOCAL ARRAYS
      CHARACTER FMT(14)*4,MESSGP(36)*4
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH,I8SAVE
      EXTERNAL I1MACH,I8SAVE
!
!  EXTERNAL SUBROUTINES
      EXTERNAL S88FMT
!
!
!  START WITH NO PREVIOUS MESSAGE.
!
      DATA MESSGP(1)/'1'/, NWP/0/, NERRP/0/
!
!  SET UP THE FORMAT FOR PRINTING THE ERROR MESSAGE.
!  THE FORMAT IS SIMPLY (A1,14X,72AXX) WHERE XX=I1MACH(6) IS THE
!  NUMBER OF CHARACTERS STORED PER INTEGER WORD.
!
      DATA CCPLUS  / '+' /
!
      DATA FMT( 1) / '(' /
      DATA FMT( 2) / 'A' /
      DATA FMT( 3) / '1' /
      DATA FMT( 4) / ',' /
      DATA FMT( 5) / '1' /
      DATA FMT( 6) / '4' /
      DATA FMT( 7) / 'X' /
      DATA FMT( 8) / ',' /
      DATA FMT( 9) / '7' /
      DATA FMT(10) / '2' /
      DATA FMT(11) / 'A' /
      DATA FMT(12) / 'X' /
      DATA FMT(13) / 'X' /
      DATA FMT(14) / ')' /
!
      IF (.NOT.SAVE) GO TO 20
!
!  SAVE THE MESSAGE.
!
        NWP=NW
        NERRP=NERR
        DO 10 I=1,NW
 10     MESSGP(I)=MESSG(I)
!
        GO TO 30
!
 20   IF (I8SAVE(1,0,.FALSE.).EQ.0) GO TO 30
!
!  PRINT THE MESSAGE.
!
        IWUNIT=I1MACH(4)
        WRITE(IWUNIT,9000) NERRP
 9000   FORMAT(' ERROR ',I4,' IN ')
!
        CALL S88FMT(2,I1MACH(6),FMT(12))
        WRITE(IWUNIT,FMT) CCPLUS,(MESSGP(I),I=1,NWP)
!
 30   RETURN
!
      END
!XERROR
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
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
      INTEGER LEVEL,NERR,NMESSG
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NMESSG)*4
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERRWV
!
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
!D9GMIT
      DOUBLE PRECISION FUNCTION D9GMIT(A,X,ALGAP1,SGNGAM,ALX)
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
      DOUBLE PRECISION A,ALGAP1,ALX,SGNGAM,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION AE,AEPS,ALG2,ALGS,BOT,EPS,FK,S,SGNG2,T,TE
      INTEGER K,M,MA
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DLNGAM
      EXTERNAL D1MACH,DLNGAM
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DBLE,DSIGN,EXP,FLOAT,LOG
!
      DATA EPS, BOT / 2*0.D0 /
!***FIRST EXECUTABLE STATEMENT  D9GMIT
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      BOT = LOG (D1MACH(1))
!
 10   IF (X.LE.0.D0) CALL XERROR ( 'D9GMIT  X SHOULD BE GT 0', 24, 1, 2)
!
      MA = A + 0.5D0
      IF (A.LT.0.D0) MA = A - 0.5D0
      AEPS = A - DBLE(FLOAT(MA))
!
      AE = A
      IF (A.LT.(-0.5D0)) AE = AEPS
!
      T = 1.D0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SER
     1IES', 54, 2, 2)
!
 30   IF (A.GE.(-0.5D0)) THEN
         ALGS = -ALGAP1 + LOG(S)
      ELSE
         ALGS = -DLNGAM(1.D0+AEPS) + LOG(S)
         S = 1.0D0
         M = -MA - 1
         IF (M.EQ.0) GO TO 50
         T = 1.0D0
         DO 40 K=1,M
            T = X*T/(AEPS-DBLE(FLOAT(M+1-K)))
            S = S + T
            IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40      CONTINUE
!
 50      D9GMIT = 0.0D0
         ALGS = -DBLE(FLOAT(MA))*LOG(X) + ALGS
         IF (S.NE.0.0D0 .AND. AEPS.NE.0.0D0) THEN
            SGNG2 = SGNGAM * DSIGN (1.0D0, S)
            ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
            IF (ALG2.GT.BOT) D9GMIT = SGNG2 * EXP(ALG2)
            IF (ALGS.GT.BOT) D9GMIT = D9GMIT + EXP(ALGS)
            RETURN
         END IF
      END IF
!
      D9GMIT = EXP (ALGS)
      RETURN
!
      END
!XSETF
      SUBROUTINE XSETF(KONTRL)
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
      INTEGER KONTRL
!
!  LOCAL SCALARS
      INTEGER JUNK
!
!  EXTERNAL FUNCTIONS
      INTEGER J4SAVE
      EXTERNAL J4SAVE
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERRWV
!
      IF ((KONTRL.GE.(-2)).AND.(KONTRL.LE.2)) GO TO 10
         CALL XERRWV('XSETF  -- INVALID VALUE OF KONTRL (I1).',33,1,2,
     1   1,KONTRL,0,0,0.,0.)
         RETURN
   10 JUNK = J4SAVE(2,KONTRL,.TRUE.)
      RETURN
      END
!R9GMIT
      REAL FUNCTION R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE TRICOMI-S INCOMPLETE GAMMA FUNCTION FOR SMALL X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL A,ALGAP1,ALX,SGNGAM,X
!
!  LOCAL SCALARS
      REAL AE,AEPS,ALG2,ALGS,BOT,EPS,FK,S,SGNG2,T,TE
      INTEGER K,M,MA
!
!  EXTERNAL FUNCTIONS
      REAL ALNGAM,R1MACH
      EXTERNAL ALNGAM,R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,EXP,FLOAT,LOG,SIGN
!
      DATA EPS, BOT / 2*0.0 /
!
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (BOT.EQ.0.0) BOT = LOG(R1MACH(1))
!
      IF (X.LE.0.0) CALL XERROR ('R9GMIT  X SHOULD BE GT 0', 24, 1, 2)
!
      MA = A + 0.5
      IF (A.LT.0.0) MA = A - 0.5
      AEPS = A - FLOAT(MA)
!
      AE = A
      IF (A.LT.(-0.5)) AE = AEPS
!
      T = 1.0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERROR (  'R9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SE
     1RIES', 54, 2, 2)
!
 30   IF (A.GE.(-0.5)) THEN
         ALGS = -ALGAP1 + LOG(S)
      ELSE
!
         ALGS = -ALNGAM(1.0+AEPS) + LOG(S)
         S = 1.0
         M = -MA - 1
         IF (M.EQ.0) GO TO 50
         T = 1.0
         DO 40 K=1,M
            T = X*T/(AEPS-FLOAT(M+1-K))
            S = S + T
            IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40      CONTINUE
!
 50      R9GMIT = 0.0
         ALGS = -FLOAT(MA)*LOG(X) + ALGS
         IF (S.NE.0.0 .AND. AEPS.NE.0.0) THEN
            SGNG2 = SGNGAM*SIGN(1.0,S)
            ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
            IF (ALG2.GT.BOT) R9GMIT = SGNG2*EXP(ALG2)
            IF (ALGS.GT.BOT) R9GMIT = R9GMIT + EXP(ALGS)
            RETURN
         END IF
      END IF
      R9GMIT = EXP(ALGS)
      RETURN
!
      END
!DLNREL
      DOUBLE PRECISION FUNCTION DLNREL (X)
! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION XMIN
      INTEGER NLNREL
!
!  LOCAL ARRAYS
      DOUBLE PRECISION ALNRCS(43)
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DCSEVL
      INTEGER INITDS
      EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSQRT,LOG,SNGL
!
!
! SERIES FOR ALNR       ON THE INTERVAL -3.75000E-01 TO  3.75000E-01
!                                        WITH WEIGHTED ERROR   6.35E-32
!                                         LOG WEIGHTED ERROR  31.20
!                               SIGNIFICANT FIGURES REQUIRED  30.93
!                                    DECIMAL PLACES REQUIRED  32.01
!
      DATA ALNRCS(  1) / +.1037869356 2743769800 6862677190 98 D+1     /
      DATA ALNRCS(  2) / -.1336430150 4908918098 7660415531 33 D+0     /
      DATA ALNRCS(  3) / +.1940824913 5520563357 9261993747 50 D-1     /
      DATA ALNRCS(  4) / -.3010755112 7535777690 3765377765 92 D-2     /
      DATA ALNRCS(  5) / +.4869461479 7154850090 4563665091 37 D-3     /
      DATA ALNRCS(  6) / -.8105488189 3175356066 8099430086 22 D-4     /
      DATA ALNRCS(  7) / +.1377884779 9559524782 9382514960 59 D-4     /
      DATA ALNRCS(  8) / -.2380221089 4358970251 3699929149 35 D-5     /
      DATA ALNRCS(  9) / +.4164041621 3865183476 3918599019 89 D-6     /
      DATA ALNRCS( 10) / -.7359582837 8075994984 2668370319 98 D-7     /
      DATA ALNRCS( 11) / +.1311761187 6241674949 1522943450 11 D-7     /
      DATA ALNRCS( 12) / -.2354670931 7742425136 6960923301 75 D-8     /
      DATA ALNRCS( 13) / +.4252277327 6034997775 6380529625 67 D-9     /
      DATA ALNRCS( 14) / -.7719089413 4840796826 1081074933 00 D-10    /
      DATA ALNRCS( 15) / +.1407574648 1359069909 2153564721 91 D-10    /
      DATA ALNRCS( 16) / -.2576907205 8024680627 5370786275 84 D-11    /
      DATA ALNRCS( 17) / +.4734240666 6294421849 1543950059 38 D-12    /
      DATA ALNRCS( 18) / -.8724901267 4742641745 3012632926 75 D-13    /
      DATA ALNRCS( 19) / +.1612461490 2740551465 7398331191 15 D-13    /
      DATA ALNRCS( 20) / -.2987565201 5665773006 7107924168 15 D-14    /
      DATA ALNRCS( 21) / +.5548070120 9082887983 0413216972 79 D-15    /
      DATA ALNRCS( 22) / -.1032461915 8271569595 1413339619 32 D-15    /
      DATA ALNRCS( 23) / +.1925023920 3049851177 8785032448 68 D-16    /
      DATA ALNRCS( 24) / -.3595507346 5265150011 1897078442 66 D-17    /
      DATA ALNRCS( 25) / +.6726454253 7876857892 1945742267 73 D-18    /
      DATA ALNRCS( 26) / -.1260262416 8735219252 0824256375 46 D-18    /
      DATA ALNRCS( 27) / +.2364488440 8606210044 9161589555 19 D-19    /
      DATA ALNRCS( 28) / -.4441937705 0807936898 8783891797 33 D-20    /
      DATA ALNRCS( 29) / +.8354659446 4034259016 2412939946 66 D-21    /
      DATA ALNRCS( 30) / -.1573155941 6479562574 8992535210 66 D-21    /
      DATA ALNRCS( 31) / +.2965312874 0247422686 1543697066 66 D-22    /
      DATA ALNRCS( 32) / -.5594958348 1815947292 1560132266 66 D-23    /
      DATA ALNRCS( 33) / +.1056635426 8835681048 1872841386 66 D-23    /
      DATA ALNRCS( 34) / -.1997248368 0670204548 3149994666 66 D-24    /
      DATA ALNRCS( 35) / +.3778297781 8839361421 0498559999 99 D-25    /
      DATA ALNRCS( 36) / -.7153158688 9081740345 0381653333 33 D-26    /
      DATA ALNRCS( 37) / +.1355248846 3674213646 5020245333 33 D-26    /
      DATA ALNRCS( 38) / -.2569467304 8487567430 0798293333 33 D-27    /
      DATA ALNRCS( 39) / +.4874775606 6216949076 4595199999 99 D-28    /
      DATA ALNRCS( 40) / -.9254211253 0849715321 1323733333 33 D-29    /
      DATA ALNRCS( 41) / +.1757859784 1760239233 2697600000 00 D-29    /
      DATA ALNRCS( 42) / -.3341002667 7731010351 3770666666 66 D-30    /
      DATA ALNRCS( 43) / +.6353393618 0236187354 1802666666 66 D-31    /
!
      DATA NLNREL, XMIN / 0, 0.D0 /
!
      IF (NLNREL.NE.0) GO TO 10
      NLNREL = INITDS (ALNRCS, 43, 0.1*SNGL(D1MACH(3)))
      XMIN = -1.0D0 + DSQRT(D1MACH(4))
!
 10   IF (X.LE.(-1.D0)) CALL XERROR ('DLNREL  X IS LE -1', 18, 2, 2)
      IF (X.LT.XMIN) CALL XERROR (
     1  'DLNREL  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 54,
     2  1, 1)
!
      IF (ABS(X).LE.0.375D0) THEN
         DLNREL = X*(1.0D0 - X*DCSEVL (X/0.375D0, ALNRCS, NLNREL))
      ELSE
         DLNREL = LOG (1.0D0+X)
      END IF
!
      RETURN
      END
!GAMI
      REAL FUNCTION GAMI (A, X)
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
      REAL A,X
!
!  LOCAL SCALARS
      REAL FACTOR
!
!  EXTERNAL FUNCTIONS
      REAL ALNGAM,GAMIT,R1MACH
      EXTERNAL ALNGAM,GAMIT,R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC EXP,LOG
!
      IF (A.LE.0.0) CALL XERROR ('GAMI    A MUST BE GT ZERO', 25, 1, 2)
      IF (X.LT.0.0) CALL XERROR ('GAMI    X MUST BE GE ZERO', 25, 2, 2)
!
      GAMI = 0.0
      IF (X.EQ.0.0) RETURN
!
! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
!
      FACTOR = ALNGAM(A) + A*LOG(X)
      IF (FACTOR.GT.LOG(R1MACH(2))) THEN
         GAMI = R1MACH(2)
      ELSE
         GAMI = EXP(FACTOR) * GAMIT(A,X)
      END IF
!
      RETURN
      END
!D9LGIC
      DOUBLE PRECISION FUNCTION D9LGIC(A,X,ALX)
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
      DOUBLE PRECISION A,ALX,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION EPS,FK,P,R,S,T,XMA,XPA
      INTEGER K
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG
!
      DATA EPS / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIC
      IF (EPS.EQ.0.D0) EPS = 0.5D0*D1MACH(3)
!
      XPA = X + 1.0D0 - A
      XMA = X - 1.D0 - A
!
      R = 0.D0
      P = 1.D0
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.D0+R)
        R = -T/((XMA+2.D0*FK)*(XPA+2.D0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'D9LGIC  NO CONVERGENCE IN 300 TERMS OF CONTINUED FR
     1ACTION', 57, 1, 2)
!
 20   D9LGIC = A*ALX - X + LOG(S/XPA)
!
      RETURN
      END
!DERF
      DOUBLE PRECISION FUNCTION DERF (X)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION SQEPS,SQRTPI,XBIG,Y
      INTEGER NTERF
!
!  LOCAL ARRAYS
      DOUBLE PRECISION ERFCS(21)
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DCSEVL,DERFC
      INTEGER INITDS
      EXTERNAL D1MACH,DCSEVL,DERFC,INITDS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSIGN,DSQRT,LOG,SNGL
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   1.28E-32
!                                         LOG WEIGHTED ERROR  31.89
!                               SIGNIFICANT FIGURES REQUIRED  31.05
!                                    DECIMAL PLACES REQUIRED  32.55
!
      DATA ERF CS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERF CS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERF CS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERF CS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERF CS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERF CS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERF CS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERF CS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERF CS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERF CS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERF CS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERF CS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERF CS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERF CS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERF CS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERF CS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERF CS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERF CS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERF CS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERF CS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERF CS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
!
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA NTERF, XBIG, SQEPS / 0, 2*0.D0 /
!
      IF (NTERF.NE.0) GO TO 10
      NTERF = INITDS (ERFCS, 21, 0.1*SNGL(D1MACH(3)))
      XBIG = DSQRT (-LOG(SQRTPI*D1MACH(3)))
      SQEPS = DSQRT (2.0D0*D1MACH(3))
!
 10   Y = ABS(X)
      IF (Y.GT.1.D0) GO TO 20
!
! ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
!
      IF (Y.LE.SQEPS) THEN
         DERF = 2.0D0*X*X/SQRTPI
      ELSE
         DERF = X*(1.0D0+DCSEVL(2.D0*X*X-1.D0,ERFCS,NTERF))
      END IF
!
      RETURN
!
! ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
!
 20   IF (Y.LE.XBIG) THEN
         DERF = DSIGN (1.0D0-DERFC(Y), X)
      ELSE
         DERF = DSIGN (1.0D0, X)
      END IF
!
      RETURN
      END
!XERRWV
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
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
!    1   1,NUM,0,0,0.,0.)
!        CALL XERRWV(  'QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM
!    1 (R2).',54,77,1,0,0,0,2,ERRREQ,ERRMIN)
!
!     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
!     LATEST REVISION ---  19 MAR 1980
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL R1,R2
      INTEGER I1,I2,LEVEL,NERR,NI,NMESSG,NR
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NMESSG)*4
!
!  LOCAL SCALARS
      INTEGER IFATAL,IUNIT,JUNK,KDUMMY,KOUNT,KUNIT,LERR,LKNTRL,LLEVEL,
     +   LMESSG,MAXMES,MKNTRL,NUNIT
      CHARACTER LFIRST*4
!
!  LOCAL ARRAYS
      INTEGER LUN(5)
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH,J4SAVE
      EXTERNAL I1MACH,J4SAVE
!
!  EXTERNAL SUBROUTINES
      EXTERNAL FDUMP,XERABT,XERCTL,XERPRT,XERSAV,XGETUA
!
!  INTRINSIC FUNCTIONS
      INTRINSIC IABS,MAX,MIN
!
!     GET FLAGS
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
!     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1   29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
!     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
!     LET USER OVERRIDE
      LFIRST = MESSG(1)
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
!     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX(-2,MIN(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
!     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
!           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
!        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            IF (NI.GE.1) WRITE (IUNIT,22) I1
            IF (NI.GE.2) WRITE (IUNIT,23) I2
            IF (NR.GE.1) WRITE (IUNIT,24) R1
            IF (NR.GE.2) WRITE (IUNIT,25) R2
   22       FORMAT (11X,'IN ABOVE MESSAGE, I1=',I10)
   23       FORMAT (11X,'IN ABOVE MESSAGE, I2=',I10)
   24       FORMAT (11X,'IN ABOVE MESSAGE, R1=',E20.10)
   25       FORMAT (11X,'IN ABOVE MESSAGE, R2=',E20.10)
            IF (LKNTRL.LE.0) GO TO 40
!              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (' ERROR NUMBER =',I10)
   40       CONTINUE
   50    CONTINUE
!        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
!     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX(1,MAXMES))) GO TO 120
!        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
!        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
!     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
!D9LGIT
      DOUBLE PRECISION FUNCTION D9LGIT(A,X,ALGAP1)
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
      DOUBLE PRECISION A,ALGAP1,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION A1X,AX,EPS,FK,HSTAR,P,R,S,SQEPS,T
      INTEGER K
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSQRT,LOG
!
      DATA EPS, SQEPS / 2*0.D0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIT
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      SQEPS = DSQRT (D1MACH(4))
!
 10   IF (X.LE.0.D0 .OR. A.LT.X) CALL XERROR ( 'D9LGIT  X SHOULD BE GT 0
     1.0 AND LE A', 35, 2, 2)
!
      AX = A + X
      A1X = AX + 1.0D0
      R = 0.D0
      P = 1.D0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.D0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 3, 2)
!
 30   HSTAR = 1.0D0 - X*S/A1X
      IF (HSTAR.LT.SQEPS) CALL XERROR ( 'D9LGIT  RESULT LESS THAN HALF P
     1RECISION', 39, 1, 1)
!
      D9LGIT = -X - ALGAP1 - LOG(HSTAR)
      RETURN
!
      END
!SETERR
      SUBROUTINE SETERR(MESSG,NMESSG,NERR,IOPT)
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
      INTEGER IOPT,NERR,NMESSG
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NMESSG)*4
!
!  LOCAL SCALARS
      INTEGER ITEMP,IWUNIT,NW
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH,I8SAVE
      EXTERNAL I1MACH,I8SAVE
!
!  EXTERNAL SUBROUTINES
      EXTERNAL E9RINT,EPRINT,FDUMP
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MIN
!
!
!  THE UNIT FOR ERROR MESSAGES.
!
      IWUNIT=I1MACH(4)
!
      IF (NMESSG.GE.1) GO TO 10
!
!  A MESSAGE OF NON-POSITIVE LENGTH IS FATAL.
!
        WRITE(IWUNIT,9000)
 9000   FORMAT('1ERROR    1 IN SETERR - MESSAGE LENGTH NOT POSITIVE.')
        GO TO 60
!
!  NW IS THE NUMBER OF WORDS THE MESSAGE OCCUPIES.
!
 10   NW=(MIN(NMESSG,72)-1)/I1MACH(6)+1
!
      IF (NERR.NE.0) GO TO 20
!
!  CANNOT TURN THE ERROR STATE OFF USING SETERR.
!
        WRITE(IWUNIT,9001)
 9001   FORMAT('1ERROR    2 IN SETERR - CANNOT HAVE NERR=0'//
     1         ' THE CURRENT ERROR MESSAGE FOLLOWS'///)
        CALL E9RINT(MESSG,NW,NERR,.TRUE.)
        ITEMP=I8SAVE(1,1,.TRUE.)
        GO TO 50
!
!  SET LERROR AND TEST FOR A PREVIOUS UNRECOVERED ERROR.
!
 20   IF (I8SAVE(1,NERR,.TRUE.).EQ.0) GO TO 30
!
        WRITE(IWUNIT,9002)
 9002   FORMAT('1ERROR    3 IN SETERR -',
     1         ' AN UNRECOVERED ERROR FOLLOWED BY ANOTHER ERROR.'//
     2         ' THE PREVIOUS AND CURRENT ERROR MESSAGES FOLLOW.'///)
        CALL EPRINT
        CALL E9RINT(MESSG,NW,NERR,.TRUE.)
        GO TO 50
!
!  SAVE THIS MESSAGE IN CASE IT IS NOT RECOVERED FROM PROPERLY.
!
 30   CALL E9RINT(MESSG,NW,NERR,.TRUE.)
!
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) GO TO 40
!
!  MUST HAVE IOPT = 1 OR 2.
!
        WRITE(IWUNIT,9003)
 9003   FORMAT('1ERROR    4 IN SETERR - BAD VALUE FOR IOPT'//
     1         ' THE CURRENT ERROR MESSAGE FOLLOWS'///)
        GO TO 50
!
!  TEST FOR RECOVERY.
!
 40   IF (IOPT.EQ.2) GO TO 50
!
      IF (I8SAVE(2,0,.FALSE.).EQ.1) RETURN
!
      CALL EPRINT
      STOP
!
 50   CALL EPRINT
 60   CALL FDUMP
      STOP
!
      END
!FDUMP
      SUBROUTINE FDUMP
!  THIS IS A DUMMY ROUTINE TO BE SENT OUT ON
!  THE PORT SEDIT TAPE
!
      RETURN
      END
!R9LGIT
      REAL FUNCTION R9LGIT (A, X, ALGAP1)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG OF TRICOMI-S INCOMPLETE GAMMA FUNCTION WITH PERRON-S
! CONTINUED FRACTION FOR LARGE X AND FOR A .GE. X.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL A,ALGAP1,X
!
!  LOCAL SCALARS
      REAL A1X,AX,EPS,FK,HSTAR,P,R,S,SQEPS,T
      INTEGER K
!
!  EXTERNAL FUNCTIONS
      REAL R1MACH
      EXTERNAL R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,SQRT
!
      DATA EPS, SQEPS / 2*0.0 /
!
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (SQEPS.EQ.0.0) SQEPS = SQRT(R1MACH(4))
!
      IF (X.LE.0.0 .OR. A.LT.X) CALL XERROR (
     1  'R9LGIT  X SHOULD BE GT 0.0 AND LE A', 35, 2, 2)
!
      AX = A + X
      A1X = AX + 1.0
      R = 0.0
      P = 1.0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR (  'R9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED F
     1RACTION', 57, 3, 2)
!
 30   HSTAR = 1.0 - X*S/A1X
      IF (HSTAR.LT.SQEPS) CALL XERROR (
     1  'R9LGIT  RESULT LESS THAN HALF PRECISION', 39, 1, 1)
!
      R9LGIT = -X - ALGAP1 - LOG(HSTAR)
!
      RETURN
      END
!XERCTL
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
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
      INTEGER KONTRL,LEVEL,MESSG1,NERR,NMESSG
!
      RETURN
      END
!S88FMT
      SUBROUTINE S88FMT( N, W, IFMT )
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
      INTEGER N,W
!
!  ARRAY ARGUMENTS
      CHARACTER IFMT(N)*4
!
!  LOCAL SCALARS
      INTEGER IDIGIT,NT,WT
!
!  LOCAL ARRAYS
      CHARACTER DIGITS(10)*4
!
!  INTRINSIC FUNCTIONS
      INTRINSIC MOD
!
!
      DATA DIGITS( 1) / '0' /
      DATA DIGITS( 2) / '1' /
      DATA DIGITS( 3) / '2' /
      DATA DIGITS( 4) / '3' /
      DATA DIGITS( 5) / '4' /
      DATA DIGITS( 6) / '5' /
      DATA DIGITS( 7) / '6' /
      DATA DIGITS( 8) / '7' /
      DATA DIGITS( 9) / '8' /
      DATA DIGITS(10) / '9' /
!
      NT = N
      WT = W
!
 10   IF (NT .LE. 0) RETURN
        IDIGIT = MOD( WT, 10 )
        IFMT(NT) = DIGITS(IDIGIT+1)
        WT = WT/10
        NT = NT - 1
        GO TO 10
!
      END
!ERF
      REAL FUNCTION ERF (X)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL SQEPS,SQRTPI,XBIG,Y
      INTEGER NTERF
!
!  LOCAL ARRAYS
      REAL ERFCS(13)
!
!  EXTERNAL FUNCTIONS
      REAL CSEVL,ERFC,R1MACH
      INTEGER INITS
      EXTERNAL CSEVL,ERFC,R1MACH,INITS
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,SIGN,SQRT
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000D+00
!                                        WITH WEIGHTED ERROR   7.10E-18
!                                         LOG WEIGHTED ERROR  17.15
!                               SIGNIFICANT FIGURES REQUIRED  16.31
!                                    DECIMAL PLACES REQUIRED  17.71
!
      DATA ERF CS( 1) /   -.0490461212 34691808E0 /
      DATA ERF CS( 2) /   -.1422612051 0371364E0 /
      DATA ERF CS( 3) /    .0100355821 87599796E0 /
      DATA ERF CS( 4) /   -.0005768764 69976748E0 /
      DATA ERF CS( 5) /    .0000274199 31252196E0 /
      DATA ERF CS( 6) /   -.0000011043 17550734E0 /
      DATA ERF CS( 7) /    .0000000384 88755420E0 /
      DATA ERF CS( 8) /   -.0000000011 80858253E0 /
      DATA ERF CS( 9) /    .0000000000 32334215E0 /
      DATA ERF CS(10) /   -.0000000000 00799101E0 /
      DATA ERF CS(11) /    .0000000000 00017990E0 /
      DATA ERF CS(12) /   -.0000000000 00000371E0 /
      DATA ERF CS(13) /    .0000000000 00000007E0 /
!
      DATA SQRTPI /1.772453850 9055160E0/
      DATA NTERF, XBIG, SQEPS / 0, 0., 0./
!
      IF (NTERF.NE.0) GO TO 10
      NTERF = INITS (ERFCS, 13, 0.1*R1MACH(3))
      XBIG = SQRT(-LOG(SQRTPI*R1MACH(3)))
      SQEPS = SQRT(2.0*R1MACH(3))
!
 10   Y = ABS(X)
      IF (Y.GT.1.) GO TO 20
!
! ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1.
!
      IF (Y.LE.SQEPS) THEN
         ERF = 2.0*X/SQRTPI
      ELSE
         ERF = X*(1.0 + CSEVL(2.*X**2-1., ERFCS, NTERF))
      END IF
!
      RETURN
!
! ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1.
!
 20   IF (Y.LE.XBIG) THEN
         ERF = SIGN (1.0-ERFC(Y), X)
      ELSE
         ERF = SIGN (1.0, X)
      END IF
!
      RETURN
      END
!XERABT
      SUBROUTINE XERABT(MESSG,NMESSG)
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
      INTEGER NMESSG
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(1)*4
!
      STOP
      END
!ALNREL
      REAL FUNCTION ALNREL (X)
! APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL XMIN
      INTEGER NLNREL
!
!  LOCAL ARRAYS
      REAL ALNRCS(23)
!
!  EXTERNAL FUNCTIONS
      REAL CSEVL,R1MACH
      INTEGER INITS
      EXTERNAL CSEVL,R1MACH,INITS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,SQRT
!
! SERIES FOR ALNR       ON THE INTERVAL -3.75000D-01 TO  3.75000D-01
!                                        WITH WEIGHTED ERROR   1.93E-17
!                                         LOG WEIGHTED ERROR  16.72
!                               SIGNIFICANT FIGURES REQUIRED  16.44
!                                    DECIMAL PLACES REQUIRED  17.40
!
      DATA ALNRCS( 1) /   1.0378693562 743770E0 /
      DATA ALNRCS( 2) /   -.1336430150 4908918E0 /
      DATA ALNRCS( 3) /    .0194082491 35520563E0 /
      DATA ALNRCS( 4) /   -.0030107551 12753577E0 /
      DATA ALNRCS( 5) /    .0004869461 47971548E0 /
      DATA ALNRCS( 6) /   -.0000810548 81893175E0 /
      DATA ALNRCS( 7) /    .0000137788 47799559E0 /
      DATA ALNRCS( 8) /   -.0000023802 21089435E0 /
      DATA ALNRCS( 9) /    .0000004164 04162138E0 /
      DATA ALNRCS(10) /   -.0000000735 95828378E0 /
      DATA ALNRCS(11) /    .0000000131 17611876E0 /
      DATA ALNRCS(12) /   -.0000000023 54670931E0 /
      DATA ALNRCS(13) /    .0000000004 25227732E0 /
      DATA ALNRCS(14) /   -.0000000000 77190894E0 /
      DATA ALNRCS(15) /    .0000000000 14075746E0 /
      DATA ALNRCS(16) /   -.0000000000 02576907E0 /
      DATA ALNRCS(17) /    .0000000000 00473424E0 /
      DATA ALNRCS(18) /   -.0000000000 00087249E0 /
      DATA ALNRCS(19) /    .0000000000 00016124E0 /
      DATA ALNRCS(20) /   -.0000000000 00002987E0 /
      DATA ALNRCS(21) /    .0000000000 00000554E0 /
      DATA ALNRCS(22) /   -.0000000000 00000103E0 /
      DATA ALNRCS(23) /    .0000000000 00000019E0 /
!
      DATA NLNREL, XMIN /0, 0./
!
      IF (NLNREL.NE.0) GO TO 10
      NLNREL = INITS (ALNRCS, 23, 0.1*R1MACH(3))
      XMIN = -1.0 + SQRT(R1MACH(4))
!
 10   IF (X.LE.(-1.0)) CALL XERROR (
     1  'ALNREL  X IS LE -1', 18, 2, 2)
      IF (X.LT.XMIN) CALL XERROR (
     1  'ALNREL  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 54,
     2  1, 1)
!
      IF (ABS(X).LE.0.375) THEN
         ALNREL = X*(1.0-X*CSEVL(X/0.375,ALNRCS,NLNREL))
      ELSE
         ALNREL = LOG (1.0+X)
      END IF
!
      RETURN
      END
!INITDS
      INTEGER FUNCTION INITDS (DOS, NOS, ETA)
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
      REAL ETA
      INTEGER NOS
!
!  ARRAY ARGUMENTS
      DOUBLE PRECISION DOS(NOS)
!
!  LOCAL SCALARS
      REAL ERR
      INTEGER I,II
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,SNGL
!
!
      IF (NOS.LT.1) CALL XERROR (
     1  'INITDS  NUMBER OF COEFFICIENTS LT 1', 35, 2, 2)
!
      ERR = 0.0
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(SNGL(DOS(I)))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
!
 20   IF (I.EQ.NOS) CALL XERROR ('INITDS  ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITDS = I
!
      RETURN
      END
!DBETAI
      DOUBLE PRECISION FUNCTION DBETAI (X, PIN, QIN)
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
      DOUBLE PRECISION PIN,QIN,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION ALNEPS,ALNSML,C,EPS,FAC1,FAC2,FINSUM,P,PS,Q,SML,
     +   TERM,XB,Y
      REAL P1
      INTEGER I,IB,N
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DLBETA
      EXTERNAL D1MACH,DLBETA
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DBLE,EXP,FLOAT,INT,LOG,MAX,MIN,SNGL
!
      DATA             EPS, ALNEPS, SML, ALNSML / 4*0.0D0 /
!
      IF (EPS.NE.0.0D0) GO TO 10
      EPS = D1MACH(3)
      ALNEPS = LOG (EPS)
      SML = D1MACH(1)
      ALNSML = LOG (SML)
!
 10   IF (X.LT.0.D0 .OR. X.GT.1.D0) CALL XERROR (
     1  'DBETAI  X IS NOT IN THE RANGE (0,1)', 35, 1, 2)
      IF (PIN.LE.0.D0 .OR. QIN.LE.0.D0) CALL XERROR (
     1  'DBETAI  P AND/OR Q IS LE ZERO', 29, 2, 2)
!
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8D0) GO TO 20
      IF (X.LT.0.2D0) GO TO 20
      Y = 1.0D0 - Y
      P = QIN
      Q = PIN
!
 20   IF ((P+Q)*Y/(P+1.D0).LT.EPS) GO TO 80
!
! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
!
      PS = Q - INT(Q)
      IF (PS.EQ.0.D0) PS = 1.0D0
      XB = P*LOG(Y) - DLBETA(PS,P) - LOG(P)
      DBETAI = 0.0D0
      IF (XB.GE.ALNSML) THEN
         DBETAI = EXP(XB)
         FAC2 = 1.0
         IF (PS.NE.1.0D0) THEN
            FAC1 = 1.0
            N = MAX(ALNEPS/LOG(Y), 4.0D0)
            DO 30 I=1,N
               IF ((I-PS.EQ.0.0D0) .OR. (FAC1.EQ.0.0D0)) THEN
                  FAC1 = 0.0D0
               ELSE
                  IF (LOG(ABS(FAC1)) + LOG(ABS(I-PS)) + LOG(Y) -
     +                LOG(DBLE(I)) .LT. ALNSML) THEN
                     FAC1 = 0.0D0
                  ELSE
                     FAC1 = FAC1 * (I-PS)*Y/I
                  END IF
               END IF
               FAC2 = FAC2 + FAC1*P/(P+I)
 30         CONTINUE
         END IF
         DBETAI = DBETAI*FAC2
      END IF
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
      IF (Q.LE.1.0D0) GO TO 70
!
      XB = P*LOG(Y) + Q*LOG(1.0D0-Y) - DLBETA(P,Q) - LOG(Q)
      IB = MAX(SNGL(XB/ALNSML), 0.0)
      TERM = EXP (XB - DBLE(FLOAT(IB))*ALNSML )
      C = 1.0D0/(1.D0-Y)
      P1 = Q*C/(P+Q-1.D0)
!
      FINSUM = 0.0D0
      N = Q
      IF (Q.EQ.DBLE(FLOAT(N))) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0D0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
        IF (Q-I+1.0D0 .EQ. 0.0D0) THEN
          TERM = 0.0D0
        ELSE
          IF (LOG(ABS(Q-I+1.0D0)) + LOG(ABS(C)) + LOG(ABS(TERM)) -
     +        LOG(ABS(P+Q-I)) .LT. ALNSML) THEN
            TERM = 0.0D0
          ELSE
            TERM = (Q-I+1.0D0)*C*TERM/(P+Q-I)
          END IF
        END IF
!
        IF (TERM.GT.1.0D0) IB = IB - 1
        IF (TERM.GT.1.0D0) TERM = TERM*SML
!
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
!
 60   DBETAI = DBETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0D0 - DBETAI
      DBETAI = MAX (MIN (DBETAI, 1.0D0), 0.0D0)
      RETURN
!
 80   DBETAI = 0.0D0
      XB = P*LOG(MAX(Y,SML)) - LOG(P) - DLBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.0D0) DBETAI = EXP(XB)
      IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0D0 - DBETAI
!
      RETURN
      END
!I8SAVE
      INTEGER FUNCTION I8SAVE(ISW,IVALUE,SET)
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
      INTEGER ISW,IVALUE
      LOGICAL SET
!
!  LOCAL SCALARS
      INTEGER LERROR,LRECOV
!
!  LOCAL ARRAYS
      INTEGER IPARAM(2)
!
!  EQUIVALENCES
      EQUIVALENCE (IPARAM(1),LERROR), (IPARAM(2),LRECOV)
!
!
!  START EXECUTION ERROR FREE AND WITH RECOVERY TURNED OFF.
!
      DATA LERROR/0/ , LRECOV/2/
!
      I8SAVE=IPARAM(ISW)
      IF (SET) IPARAM(ISW)=IVALUE
!
      RETURN
!
      END
!DGAMMA
      DOUBLE PRECISION FUNCTION DGAMMA (X)
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION DXREL,PI,SINPIY,SQ2PIL,XMAX,XMIN,Y
      INTEGER I,N,NGAM
!
!  LOCAL ARRAYS
      DOUBLE PRECISION GAMCS(42)
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,D9LGMC,DCSEVL
      INTEGER INITDS
      EXTERNAL D1MACH,D9LGMC,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DGAMLM,XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DBLE,DSQRT,EXP,FLOAT,INT,LOG,SIN,SNGL
!
!
! SERIES FOR GAM        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   5.79E-32
!                                         LOG WEIGHTED ERROR  31.24
!                               SIGNIFICANT FIGURES REQUIRED  30.00
!                                    DECIMAL PLACES REQUIRED  32.05
!
      DATA GAM CS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAM CS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAM CS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAM CS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAM CS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAM CS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAM CS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAM CS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAM CS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAM CS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAM CS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAM CS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAM CS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAM CS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAM CS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAM CS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAM CS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAM CS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAM CS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAM CS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAM CS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAM CS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAM CS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAM CS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAM CS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAM CS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAM CS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAM CS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAM CS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAM CS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAM CS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAM CS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAM CS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAM CS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAM CS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAM CS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAM CS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAM CS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAM CS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAM CS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAM CS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAM CS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
!
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
! SQ2PIL IS 0.5*LOG(2*PI) = LOG(SQRT(2*PI))
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA NGAM, XMIN, XMAX, DXREL / 0, 3*0.D0 /
!
      IF (NGAM.NE.0) GO TO 10
      NGAM = INITDS (GAMCS, 42, 0.1*SNGL(D1MACH(3)) )
!
      CALL DGAMLM (XMIN, XMAX)
      DXREL = DSQRT (D1MACH(4))
!
 10   Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - DBLE(FLOAT(N))
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X.EQ.0.D0) CALL XERROR ('DGAMMA  X IS 0', 14, 4, 2)
      IF (X.LT.0.0 .AND. X+DBLE(FLOAT(N-2)).EQ.0.D0) CALL XERROR (
     1  'DGAMMA  X IS A NEGATIVE INTEGER', 31, 4, 2)
      IF (X.LT.(-0.5D0) .AND. ABS((X-INT(X-0.5D0))/X).LT.DXREL) CALL
     1  XERROR (  'DGAMMA  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR N
     2EGATIVE INTEGER', 68, 1, 1)
!
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+DBLE(FLOAT(I-1)) )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        DGAMMA = (Y+DBLE(FLOAT(I))) * DGAMMA
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X.GT.XMAX) CALL XERROR ('DGAMMA  X SO BIG GAMMA OVERFLOWS',
     1  32, 3, 2)
!
      DGAMMA = 0.D0
      IF (X.LT.XMIN) CALL XERROR ('DGAMMA  X SO SMALL GAMMA UNDERFLOWS',
     1  35, 2, 1)
      IF (X.LT.XMIN) RETURN
!
      DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
!
      IF (ABS((X-INT(X-0.5D0))/X).LT.DXREL) CALL XERROR (
     1  'DGAMMA  ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER',
     2  61, 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY.EQ.0.D0) CALL XERROR (
     1  'DGAMMA  X IS A NEGATIVE INTEGER', 31, 4, 2)
!
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
!
      RETURN
      END
!DERFC
      DOUBLE PRECISION FUNCTION DERFC (X)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION SQEPS,SQRTPI,XMAX,XSML,Y
      REAL ETA
      INTEGER NTERC2,NTERF,NTERFC
!
!  LOCAL ARRAYS
      DOUBLE PRECISION ERC2CS(49),ERFCCS(59),ERFCS(21)
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DCSEVL
      INTEGER INITDS
      EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSQRT,EXP,LOG,SNGL
!
!
! SERIES FOR ERF        ON THE INTERVAL  0.          TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   1.28E-32
!                                         LOG WEIGHTED ERROR  31.89
!                               SIGNIFICANT FIGURES REQUIRED  31.05
!                                    DECIMAL PLACES REQUIRED  32.55
!
      DATA ERF CS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERF CS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERF CS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERF CS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERF CS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERF CS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERF CS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERF CS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERF CS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERF CS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERF CS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERF CS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERF CS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERF CS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERF CS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERF CS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERF CS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERF CS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERF CS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERF CS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERF CS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
!
! SERIES FOR ERC2       ON THE INTERVAL  2.50000E-01 TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   2.67E-32
!                                         LOG WEIGHTED ERROR  31.57
!                               SIGNIFICANT FIGURES REQUIRED  30.31
!                                    DECIMAL PLACES REQUIRED  32.42
!
      DATA ERC2CS(  1) / -.6960134660 2309501127 3915082619 7 D-1      /
      DATA ERC2CS(  2) / -.4110133936 2620893489 8221208466 6 D-1      /
      DATA ERC2CS(  3) / +.3914495866 6896268815 6114370524 4 D-2      /
      DATA ERC2CS(  4) / -.4906395650 5489791612 8093545077 4 D-3      /
      DATA ERC2CS(  5) / +.7157479001 3770363807 6089414182 5 D-4      /
      DATA ERC2CS(  6) / -.1153071634 1312328338 0823284791 2 D-4      /
      DATA ERC2CS(  7) / +.1994670590 2019976350 5231486770 9 D-5      /
      DATA ERC2CS(  8) / -.3642666471 5992228739 3611843071 1 D-6      /
      DATA ERC2CS(  9) / +.6944372610 0050125899 3127721463 3 D-7      /
      DATA ERC2CS( 10) / -.1371220902 1043660195 3460514121 0 D-7      /
      DATA ERC2CS( 11) / +.2788389661 0071371319 6386034808 7 D-8      /
      DATA ERC2CS( 12) / -.5814164724 3311615518 6479105031 6 D-9      /
      DATA ERC2CS( 13) / +.1238920491 7527531811 8016881795 0 D-9      /
      DATA ERC2CS( 14) / -.2690639145 3067434323 9042493788 9 D-10     /
      DATA ERC2CS( 15) / +.5942614350 8479109824 4470968384 0 D-11     /
      DATA ERC2CS( 16) / -.1332386735 7581195792 8775442057 0 D-11     /
      DATA ERC2CS( 17) / +.3028046806 1771320171 7369724330 4 D-12     /
      DATA ERC2CS( 18) / -.6966648814 9410325887 9586758895 4 D-13     /
      DATA ERC2CS( 19) / +.1620854541 0539229698 1289322762 8 D-13     /
      DATA ERC2CS( 20) / -.3809934465 2504919998 7691305772 9 D-14     /
      DATA ERC2CS( 21) / +.9040487815 9788311493 6897101297 5 D-15     /
      DATA ERC2CS( 22) / -.2164006195 0896073478 0981204700 3 D-15     /
      DATA ERC2CS( 23) / +.5222102233 9958549846 0798024417 2 D-16     /
      DATA ERC2CS( 24) / -.1269729602 3645553363 7241552778 0 D-16     /
      DATA ERC2CS( 25) / +.3109145504 2761975838 3622741295 1 D-17     /
      DATA ERC2CS( 26) / -.7663762920 3203855240 0956671481 1 D-18     /
      DATA ERC2CS( 27) / +.1900819251 3627452025 3692973329 0 D-18     /
      DATA ERC2CS( 28) / -.4742207279 0690395452 2565599996 5 D-19     /
      DATA ERC2CS( 29) / +.1189649200 0765283828 8068307845 1 D-19     /
      DATA ERC2CS( 30) / -.3000035590 3257802568 4527131306 6 D-20     /
      DATA ERC2CS( 31) / +.7602993453 0432461730 1938527709 8 D-21     /
      DATA ERC2CS( 32) / -.1935909447 6068728815 6981104913 0 D-21     /
      DATA ERC2CS( 33) / +.4951399124 7733378810 0004238677 3 D-22     /
      DATA ERC2CS( 34) / -.1271807481 3363718796 0862198988 8 D-22     /
      DATA ERC2CS( 35) / +.3280049600 4695130433 1584165205 3 D-23     /
      DATA ERC2CS( 36) / -.8492320176 8228965689 2479242239 9 D-24     /
      DATA ERC2CS( 37) / +.2206917892 8075602235 1987998719 9 D-24     /
      DATA ERC2CS( 38) / -.5755617245 6965284983 1281950719 9 D-25     /
      DATA ERC2CS( 39) / +.1506191533 6392342503 5414405119 9 D-25     /
      DATA ERC2CS( 40) / -.3954502959 0187969531 0428569599 9 D-26     /
      DATA ERC2CS( 41) / +.1041529704 1515009799 8464505173 3 D-26     /
      DATA ERC2CS( 42) / -.2751487795 2787650794 5017890133 3 D-27     /
      DATA ERC2CS( 43) / +.7290058205 4975574089 9770368000 0 D-28     /
      DATA ERC2CS( 44) / -.1936939645 9159478040 7750109866 6 D-28     /
      DATA ERC2CS( 45) / +.5160357112 0514872983 7005482666 6 D-29     /
      DATA ERC2CS( 46) / -.1378419322 1930940993 8964480000 0 D-29     /
      DATA ERC2CS( 47) / +.3691326793 1070690422 5109333333 3 D-30     /
      DATA ERC2CS( 48) / -.9909389590 6243654206 5322666666 6 D-31     /
      DATA ERC2CS( 49) / +.2666491705 1953884133 2394666666 6 D-31     /
!
! SERIES FOR ERFC       ON THE INTERVAL  0.          TO  2.50000E-01
!                                        WITH WEIGHTED ERROR   1.53E-31
!                                         LOG WEIGHTED ERROR  30.82
!                               SIGNIFICANT FIGURES REQUIRED  29.47
!                                    DECIMAL PLACES REQUIRED  31.70
!
      DATA ERFCCS(  1) / +.7151793102 0292477450 3697709496 D-1        /
      DATA ERFCCS(  2) / -.2653243433 7606715755 8893386681 D-1        /
      DATA ERFCCS(  3) / +.1711153977 9208558833 2699194606 D-2        /
      DATA ERFCCS(  4) / -.1637516634 5851788416 3746404749 D-3        /
      DATA ERFCCS(  5) / +.1987129350 0552036499 5974806758 D-4        /
      DATA ERFCCS(  6) / -.2843712412 7665550875 0175183152 D-5        /
      DATA ERFCCS(  7) / +.4606161308 9631303696 9379968464 D-6        /
      DATA ERFCCS(  8) / -.8227753025 8792084205 7766536366 D-7        /
      DATA ERFCCS(  9) / +.1592141872 7709011298 9358340826 D-7        /
      DATA ERFCCS( 10) / -.3295071362 2528432148 6631665072 D-8        /
      DATA ERFCCS( 11) / +.7223439760 4005554658 1261153890 D-9        /
      DATA ERFCCS( 12) / -.1664855813 3987295934 4695966886 D-9        /
      DATA ERFCCS( 13) / +.4010392588 2376648207 7671768814 D-10       /
      DATA ERFCCS( 14) / -.1004816214 4257311327 2170176283 D-10       /
      DATA ERFCCS( 15) / +.2608275913 3003338085 9341009439 D-11       /
      DATA ERFCCS( 16) / -.6991110560 4040248655 7697812476 D-12       /
      DATA ERFCCS( 17) / +.1929492333 2617070862 4205749803 D-12       /
      DATA ERFCCS( 18) / -.5470131188 7543310649 0125085271 D-13       /
      DATA ERFCCS( 19) / +.1589663309 7626974483 9084032762 D-13       /
      DATA ERFCCS( 20) / -.4726893980 1975548392 0369584290 D-14       /
      DATA ERFCCS( 21) / +.1435873376 7849847867 2873997840 D-14       /
      DATA ERFCCS( 22) / -.4449510561 8173583941 7250062829 D-15       /
      DATA ERFCCS( 23) / +.1404810884 7682334373 7305537466 D-15       /
      DATA ERFCCS( 24) / -.4513818387 7642108962 5963281623 D-16       /
      DATA ERFCCS( 25) / +.1474521541 0451330778 7018713262 D-16       /
      DATA ERFCCS( 26) / -.4892621406 9457761543 6841552532 D-17       /
      DATA ERFCCS( 27) / +.1647612141 4106467389 5301522827 D-17       /
      DATA ERFCCS( 28) / -.5626817176 3294080929 9928521323 D-18       /
      DATA ERFCCS( 29) / +.1947443382 2320785142 9197867821 D-18       /
      DATA ERFCCS( 30) / -.6826305642 9484207295 6664144723 D-19       /
      DATA ERFCCS( 31) / +.2421988887 2986492401 8301125438 D-19       /
      DATA ERFCCS( 32) / -.8693414133 5030704256 3800861857 D-20       /
      DATA ERFCCS( 33) / +.3155180346 2280855712 2363401262 D-20       /
      DATA ERFCCS( 34) / -.1157372324 0496087426 1239486742 D-20       /
      DATA ERFCCS( 35) / +.4288947161 6056539462 3737097442 D-21       /
      DATA ERFCCS( 36) / -.1605030742 0576168500 5737770964 D-21       /
      DATA ERFCCS( 37) / +.6063298757 4538026449 5069923027 D-22       /
      DATA ERFCCS( 38) / -.2311404251 6979584909 8840801367 D-22       /
      DATA ERFCCS( 39) / +.8888778540 6618855255 4702955697 D-23       /
      DATA ERFCCS( 40) / -.3447260576 6513765223 0718495566 D-23       /
      DATA ERFCCS( 41) / +.1347865460 2069650682 7582774181 D-23       /
      DATA ERFCCS( 42) / -.5311794071 1250217364 5873201807 D-24       /
      DATA ERFCCS( 43) / +.2109341058 6197831682 8954734537 D-24       /
      DATA ERFCCS( 44) / -.8438365587 9237891159 8133256738 D-25       /
      DATA ERFCCS( 45) / +.3399982524 9452089062 7359576337 D-25       /
      DATA ERFCCS( 46) / -.1379452388 0732420900 2238377110 D-25       /
      DATA ERFCCS( 47) / +.5634490311 8332526151 3392634811 D-26       /
      DATA ERFCCS( 48) / -.2316490434 4770654482 3427752700 D-26       /
      DATA ERFCCS( 49) / +.9584462844 6018101526 3158381226 D-27       /
      DATA ERFCCS( 50) / -.3990722880 3301097262 4224850193 D-27       /
      DATA ERFCCS( 51) / +.1672129225 9444773601 7228709669 D-27       /
      DATA ERFCCS( 52) / -.7045991522 7660138563 8803782587 D-28       /
      DATA ERFCCS( 53) / +.2979768402 8642063541 2357989444 D-28       /
      DATA ERFCCS( 54) / -.1262522466 4606192972 2422632994 D-28       /
      DATA ERFCCS( 55) / +.5395438704 5424879398 5299653154 D-29       /
      DATA ERFCCS( 56) / -.2380992882 5314591867 5346190062 D-29       /
      DATA ERFCCS( 57) / +.1099052830 1027615735 9726683750 D-29       /
      DATA ERFCCS( 58) / -.4867713741 6449657273 2518677435 D-30       /
      DATA ERFCCS( 59) / +.1525877264 1103575676 3200828211 D-30       /
!
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA NTERF, NTERFC, NTERC2, XSML, XMAX, SQEPS / 3*0, 3*0.D0 /
!
      IF (NTERF.NE.0) GO TO 10
      ETA = 0.1*SNGL(D1MACH(3))
      NTERF = INITDS (ERFCS, 21, ETA)
      NTERFC = INITDS (ERFCCS, 59, ETA)
      NTERC2 = INITDS (ERC2CS, 49, ETA)
!
      XSML = -DSQRT (-LOG(SQRTPI*D1MACH(3)))
      XMAX = DSQRT (-LOG(SQRTPI*D1MACH(1)) )
      XMAX = XMAX - 0.5D0*LOG(XMAX)/XMAX - 0.01D0
      SQEPS = DSQRT (2.0D0*D1MACH(3))
!
 10   IF (X.GT.XSML) GO TO 20
!
! ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
!
      DERFC = 2.0D0
      RETURN
!
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0D0) GO TO 30
!
! ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
!
      IF (Y.LT.SQEPS) THEN
         DERFC = 1.0D0 - 2.0D0*X/SQRTPI
      ELSE
         DERFC = 1.0D0 - X*(1.0D0+DCSEVL(2.0D0*X*X-1.0D0,ERFCS,NTERF))
      END IF
!
      RETURN
!
! ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
!
 30   Y = Y*Y
      IF (Y.LE.4.0D0) THEN
         DERFC = EXP(-Y)/ABS(X) *
     +           (0.5D0 + DCSEVL((8.0D0/Y-5.0D0)/3.0D0,ERC2CS,NTERC2))
      ELSE
         DERFC = EXP(-Y)/ABS(X) *
     +           (0.5D0 + DCSEVL(8.0D0/Y-1.0D0,ERFCCS,NTERFC))
      END IF
      IF (X.LT.0.D0) DERFC = 2.0D0 - DERFC
      RETURN
!
 40   CALL XERROR ('DERFC   X SO BIG ERFC UNDERFLOWS', 32, 1, 1)
      DERFC = 0.D0
      RETURN
!
      END
!DGAMLM
      SUBROUTINE DGAMLM (XMIN, XMAX)
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
      DOUBLE PRECISION XMAX,XMIN
!
!  LOCAL SCALARS
      DOUBLE PRECISION ALNBIG,ALNSML,XLN,XOLD
      INTEGER I
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,MAX
!
!
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERROR ('DGAMLM  UNABLE TO FIND XMIN', 27, 1, 2)
!
 20   XMIN = -XMIN + 0.01D0
!
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERROR ('DGAMLM  UNABLE TO FIND XMAX', 27, 2, 2)
!
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
!
      RETURN
      END
!ALBETA
      REAL FUNCTION ALBETA (A, B)
! JULY 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL A,B
!
!  LOCAL SCALARS
      REAL CORR,P,Q,SQ2PIL
!
!  EXTERNAL FUNCTIONS
      REAL ALNGAM,ALNREL,GAMMA,R9LGMC
      EXTERNAL ALNGAM,ALNREL,GAMMA,R9LGMC
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC LOG,MAX,MIN
!
      DATA SQ2PIL / 0.91893853320467274E0 /
!
      P = MIN (A, B)
      Q = MAX (A, B)
!
      IF (P.LE.0.0) CALL XERROR (
     1  'ALBETA  BOTH ARGUMENTS MUST BE GT ZERO', 38, 1, 2)
      IF (P.GE.10.0) GO TO 30
      IF (Q.GE.10.0) GO TO 20
!
! P AND Q ARE SMALL.
!
      ALBETA = LOG(GAMMA(P) * (GAMMA(Q)/GAMMA(P+Q)) )
      RETURN
!
! P IS SMALL, BUT Q IS BIG.
!
 20   CORR = R9LGMC(Q) - R9LGMC(P+Q)
      ALBETA = ALNGAM(P) + CORR + P - P*LOG(P+Q) +
     1  (Q-0.5)*ALNREL(-P/(P+Q))
      RETURN
!
! P AND Q ARE BIG.
!
 30   CORR = R9LGMC(P) + R9LGMC(Q) - R9LGMC(P+Q)
      ALBETA = -0.5*LOG(Q) + SQ2PIL + CORR + (P-0.5)*LOG(P/(P+Q))
     1  + Q*ALNREL(-P/(P+Q))
      RETURN
!
      END
!XERSAV
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
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
      INTEGER ICOUNT,LEVEL,NERR,NMESSG
!
!  ARRAY ARGUMENTS
      CHARACTER MESSG(NMESSG)*4
!
!  LOCAL SCALARS
      INTEGER I,II,IUNIT,KOUNTX,KUNIT,NCHAR,NCOL,NUNIT
!
!  LOCAL ARRAYS
      INTEGER KOUNT(10),LEVTAB(10),LUN(5),NERTAB(10)
      CHARACTER F(17)*1,MESTAB(10)*4
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH
      EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL S88FMT,XGETUA
!
!     NEXT THREE DATA STATEMENTS ARE NEEDED MERELY TO SATISFY
!     CERTAIN CONVENTIONS FOR COMPILERS WHICH DYNAMICALLY
!     ALLOCATE STORAGE.
      DATA MESTAB(1),MESTAB(2),MESTAB(3),MESTAB(4),MESTAB(5),
     1     MESTAB(6),MESTAB(7),MESTAB(8),MESTAB(9),MESTAB(10)
     2     /'0','0','0','0','0','0','0','0','0','0'/
      DATA NERTAB(1),NERTAB(2),NERTAB(3),NERTAB(4),NERTAB(5),
     1     NERTAB(6),NERTAB(7),NERTAB(8),NERTAB(9),NERTAB(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA LEVTAB(1),LEVTAB(2),LEVTAB(3),LEVTAB(4),LEVTAB(5),
     1     LEVTAB(6),LEVTAB(7),LEVTAB(8),LEVTAB(9),LEVTAB(10)
     2     /0,0,0,0,0,0,0,0,0,0/
!     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
!     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
!     NEXT DATA STATEMENT SETS UP OUTPUT FORMAT
      DATA F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),F(9),F(10),
     1     F(11),F(12),F(13),F(14),F(15),F(16),F(17)
     2     /'(' ,'1' ,'X' ,',' ,'A' ,' ' ,' ' ,',' ,'I' ,' ' ,
     3      ' ' ,',' ,'2' ,'I' ,'1' ,'0' ,')' /
      IF (NMESSG.GT.0) GO TO 80
!     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
!        PREPARE FORMAT
         NCHAR = I1MACH(6)
         CALL S88FMT(2,NCHAR,F(6))
         NCOL = 20 - NCHAR
         CALL S88FMT(2,NCOL,F(10))
!        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT ('0          ERROR MESSAGE SUMMARY'/
     1              ' FIRST WORD      NERR     LEVEL     COUNT')
!           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,F) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   20       CONTINUE
   30       CONTINUE
!           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (/' OTHER ERRORS NOT INDIVIDUALLY TABULATED=',I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
!        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
!     PROCESS A MESSAGE...
!     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MESSG(1).NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
!     THREE POSSIBLE CASES...
!     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
!     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
!     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MESSG(1)
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
!DGAMIT
      DOUBLE PRECISION FUNCTION DGAMIT (A, X)
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
      DOUBLE PRECISION A,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION AEPS,AINTA,ALGAP1,ALNEPS,ALNG,ALX,BOT,H,SGA,
     +   SGNGAM,SQEPS,T
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,D9GMIT,D9LGIC,D9LGIT,DGAMR,DLNGAM
      EXTERNAL D1MACH,D9GMIT,D9LGIC,D9LGIT,DGAMR,DLNGAM
!
!  EXTERNAL SUBROUTINES
      EXTERNAL DLGAMS,XERCLR,XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DSIGN,DSQRT,EXP,INT,LOG
!
!
      DATA ALNEPS, SQEPS, BOT / 3*0.D0 /
!
      IF (ALNEPS.NE.0.D0) GO TO 10
      ALNEPS = -LOG (D1MACH(3))
      SQEPS = DSQRT (D1MACH(4))
      BOT = LOG (D1MACH(1))
!
 10   IF (X.LT.0.D0) CALL XERROR ('DGAMIT  X IS NEGATIVE', 21, 2, 2)
!
      IF (X.NE.0.D0) ALX = LOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = DSIGN (1.0D0, A)
      AINTA = INT (A + 0.5D0*SGA)
      AEPS = A - AINTA
!
      IF (X.GT.0.D0) GO TO 20
      DGAMIT = 0.0D0
      IF (AINTA.GT.0.D0 .OR. AEPS.NE.0.D0) DGAMIT = DGAMR(A+1.0D0)
      RETURN
!
 20   IF (X.GT.1.D0) GO TO 30
      IF (A.GE.(-0.5D0) .OR. AEPS.NE.0.D0) CALL DLGAMS (A+1.0D0, ALGAP1,
     1  SGNGAM)
      DGAMIT = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
!
 30   IF (A.LT.X) GO TO 40
      T = D9LGIT (A, X, DLNGAM(A+1.0D0))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = EXP (T)
      RETURN
!
 40   ALNG = D9LGIC (A, X, ALX)
!
! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
!
      H = 1.0D0
      IF (AEPS.EQ.0.D0 .AND. AINTA.LE.0.D0) GO TO 50
!
      CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      T = LOG (ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 60
!
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGA * SGNGAM * EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 50
!
      CALL XERCLR
      CALL XERROR ('DGAMIT  RESULT LT HALF PRECISION', 32, 1, 1)
!
 50   T = -A*ALX + LOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = DSIGN (EXP(T), H)
      RETURN
!
 60   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = -SGA * SGNGAM * EXP(T)
      RETURN
!
      END
!GAMIT
      REAL FUNCTION GAMIT (A, X)
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
      REAL A,X
!
!  LOCAL SCALARS
      REAL AEPS,AINTA,ALGAP1,ALNEPS,ALNG,ALX,BOT,H,SGA,SGNGAM,SQEPS,T
!
!  EXTERNAL FUNCTIONS
      REAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
      EXTERNAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
!
!  EXTERNAL SUBROUTINES
      EXTERNAL ALGAMS,XERCLR,XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP,LOG,SIGN,SQRT
!
      DATA ALNEPS, SQEPS, BOT / 3*0.0 /
!
      IF (ALNEPS.NE.0.0) GO TO 10
      ALNEPS = -LOG(R1MACH(3))
      SQEPS = SQRT(R1MACH(4))
      BOT = LOG(R1MACH(1))
!
 10   IF (X.LT.0.0) CALL XERROR ('GAMIT   X IS NEGATIVE', 21, 2, 2)
!
      IF (X.NE.0.0) ALX = LOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      AINTA = AINT (A+0.5*SGA)
      AEPS = A - AINTA
!
      IF (X.GT.0.0) GO TO 20
      GAMIT = 0.0
      IF (AINTA.GT.0.0 .OR. AEPS.NE.0.0) GAMIT = GAMR(A+1.0)
      RETURN
!
 20   IF (X.GT.1.0) GO TO 40
      IF (A.GE.(-0.5) .OR. AEPS.NE.0.0) CALL ALGAMS (A+1.0, ALGAP1,
     1  SGNGAM)
      GAMIT = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
!
 40   IF (A.LT.X) GO TO 50
      T = R9LGIT (A, X, ALNGAM(A+1.0))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = EXP(T)
      RETURN
!
 50   ALNG = R9LGIC (A, X, ALX)
!
! EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
!
      H = 1.0
      IF (AEPS.EQ.0.0 .AND. AINTA.LE.0.0) GO TO 60
      CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      T = LOG(ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGA*SGNGAM*EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 60
      CALL XERCLR
      CALL XERROR ('GAMIT   RESULT LT HALF PRECISION', 32, 1, 1)
!
 60   T = -A*ALX + LOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = SIGN (EXP(T), H)
      RETURN
!
 70   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = -SGA*SGNGAM*EXP(T)
      RETURN
!
      END
!DLGAMS
      SUBROUTINE DLGAMS (X, DLGAM, SGNGAM)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE LOG ABS (GAMMA(X)) AND RETURN THE SIGN OF GAMMA(X) IN SGNGAM.
! SGNGAM IS EITHER +1.0 OR -1.0.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION DLGAM,SGNGAM,X
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION DLNGAM
      EXTERNAL DLNGAM
!
!  INTRINSIC FUNCTIONS
      INTRINSIC INT,MOD
!
!
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0D0
      IF (X.GT.0.D0) RETURN
!
!     INT = DMOD (-INT(X), 2.0D0) + 0.1D0
      IF (INT(MOD(-INT(X),2)+0.1D0).EQ.0) SGNGAM = -1.0D0
!
      RETURN
      END
!GAMR
      REAL FUNCTION GAMR (X)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
! THIS ROUTINE, NOT GAMMA(X), SHOULD BE THE FUNDAMENTAL ONE.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL X
!
!  LOCAL SCALARS
      REAL ALNGX,SGNGX
      INTEGER IROLD
!
!  EXTERNAL FUNCTIONS
      REAL GAMMA
      EXTERNAL GAMMA
!
!  EXTERNAL SUBROUTINES
      EXTERNAL ALGAMS,XERCLR,XGETF,XSETF
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP
!
      GAMR = 0.0
      IF (X.LE.0.0 .AND. AINT(X).EQ.X) RETURN
!
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0) GO TO 10
      GAMR = 1.0/GAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
!
 10   CALL ALGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      GAMR = SGNGX * EXP(-ALNGX)
      RETURN
!
      END
!XERCLR
      SUBROUTINE XERCLR
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
      INTEGER JUNK
!
!  EXTERNAL FUNCTIONS
      INTEGER J4SAVE
      EXTERNAL J4SAVE
!
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
!BETAI
      REAL FUNCTION BETAI (X, PIN, QIN)
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
      REAL PIN,QIN,X
!
!  LOCAL SCALARS
      REAL ALNEPS,ALNSML,C,EPS,FAC1,FAC2,FINSUM,P,P1,PS,Q,SML,TERM,XB,Y
      INTEGER I,IB,N
!
!  EXTERNAL FUNCTIONS
      REAL ALBETA,R1MACH
      EXTERNAL ALBETA,R1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP,FLOAT,LOG,MAX,MIN,REAL
!
      DATA             EPS, ALNEPS, SML, ALNSML / 4*0.0 /
!
      IF (EPS.NE.0.) GO TO 10
      EPS = R1MACH(3)
      ALNEPS = LOG(EPS)
      SML = R1MACH(1)
      ALNSML = LOG(SML)
!
 10   IF (X.LT.0. .OR. X.GT.1.0) CALL XERROR (
     1  'BETAI   X IS NOT IN THE RANGE (0,1)', 35, 1, 2)
      IF (PIN.LE.0. .OR. QIN.LE.0.) CALL XERROR (
     1  'BETAI   P AND/OR Q IS LE ZERO', 29, 2, 2)
!
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8) GO TO 20
      IF (X.LT.0.2) GO TO 20
      Y = 1.0 - Y
      P = QIN
      Q = PIN
!
 20   IF ((P+Q)*Y/(P+1.).LT.EPS) GO TO 80
!
! EVALUATE THE INFINITE SUM FIRST.
! TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
!
      PS = Q - AINT(Q)
      IF (PS.EQ.0.) PS = 1.0
      XB = P*LOG(Y) -  ALBETA(PS, P) - LOG(P)
      BETAI = 0.0
      IF (XB.GE.ALNSML) THEN
         BETAI = EXP(XB)
         FAC2 = 1.0
         IF (PS.NE.1.0E0) THEN
            FAC1 = 1.0
            N = MAX(ALNEPS/LOG(Y), 4.0E0)
            DO 30 I=1,N
               IF ((I-PS.EQ.0.0E0) .OR. (FAC1.EQ.0.0E0)) THEN
                  FAC1 = 0.0E0
               ELSE
                  IF (LOG(ABS(FAC1)) + LOG(ABS(I-PS)) + LOG(Y) -
     +                LOG(REAL(I)) .LT. ALNSML) THEN
                     FAC1 = 0.0E0
                  ELSE
                     FAC1 = FAC1 * (I-PS)*Y/I
                  END IF
               END IF
               FAC2 = FAC2 + FAC1*P/(P+I)
 30         CONTINUE
         END IF
         BETAI = BETAI*FAC2
      END IF
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
      IF (Q.LE.1.0) GO TO 70
!
      XB = P*LOG(Y) + Q*LOG(1.0-Y) - ALBETA(P,Q) - LOG(Q)
      IB = MAX (XB/ALNSML, 0.0)
      TERM = EXP (XB - FLOAT(IB)*ALNSML)
      C = 1.0/(1.0-Y)
      P1 = Q*C/(P+Q-1.)
!
      FINSUM = 0.0
      N = Q
      IF (Q.EQ.FLOAT(N)) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
        IF (Q-I+1.0E0 .EQ. 0.0E0) THEN
          TERM = 0.0E0
        ELSE
          IF (LOG(ABS(Q-I+1.0E0)) + LOG(ABS(C)) + LOG(ABS(TERM)) -
     +        LOG(ABS(P+Q-I)) .LT. ALNSML) THEN
            TERM = 0.0E0
          ELSE
            TERM = (Q-I+1.0E0)*C*TERM/(P+Q-I)
          END IF
        END IF
!
        IF (TERM.GT.1.0) IB = IB - 1
        IF (TERM.GT.1.0) TERM = TERM*SML
!
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
!
 60   BETAI = BETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0 - BETAI
      BETAI = MAX (MIN (BETAI, 1.0), 0.0)
      RETURN
!
 80   BETAI = 0.0
      XB = P*LOG(MAX(Y,SML)) - LOG(P) - ALBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.) BETAI = EXP (XB)
      IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0 - BETAI
      RETURN
!
      END
!J4SAVE
      INTEGER FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
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
      INTEGER IVALUE,IWHICH
      LOGICAL ISET
!
!  LOCAL ARRAYS
      INTEGER IPARAM(9)
!
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,1,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
!XGETF
      SUBROUTINE XGETF(KONTRL)
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
      INTEGER KONTRL
!
!  EXTERNAL FUNCTIONS
      INTEGER J4SAVE
      EXTERNAL J4SAVE
!
      KONTRL = J4SAVE(2,0,.FALSE.)
      RETURN
      END
!D9LGMC
      DOUBLE PRECISION FUNCTION D9LGMC (X)
! AUGUST 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! COMPUTE THE LOG GAMMA CORRECTION FACTOR FOR X .GE. 10. SO THAT
! LOG (DGAMMA(X)) = LOG(DSQRT(2*PI)) + (X-.5)*LOG(X) - X + D9LGMC(X)
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      DOUBLE PRECISION X
!
!  LOCAL SCALARS
      DOUBLE PRECISION XBIG,XMAX
      INTEGER NALGM
!
!  LOCAL ARRAYS
      DOUBLE PRECISION ALGMCS(15)
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DCSEVL
      INTEGER INITDS
      EXTERNAL D1MACH,DCSEVL,INITDS
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC DSQRT,EXP,LOG,MIN,SNGL
!
!
! SERIES FOR ALGM       ON THE INTERVAL  0.          TO  1.00000E-02
!                                        WITH WEIGHTED ERROR   1.28E-31
!                                         LOG WEIGHTED ERROR  30.89
!                               SIGNIFICANT FIGURES REQUIRED  29.81
!                                    DECIMAL PLACES REQUIRED  31.48
!
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
!
      DATA NALGM, XBIG, XMAX / 0, 2*0.D0 /
!
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITDS (ALGMCS, 15, SNGL(D1MACH(3)) )
      XBIG = 1.0D0/DSQRT(D1MACH(3))
      XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
!
 10   IF (X.LT.10.D0) CALL XERROR ('D9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
!
 20   D9LGMC = 0.D0
      CALL XERROR ('D9LGMC  X SO BIG D9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
!
      END
!ALGAMS
      SUBROUTINE ALGAMS (X, ALGAM, SGNGAM)
! JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
!
! EVALUATE LOG ABS (GAMMA(X)) AND RETURN THE SIGN OF GAMMA(X) IN SGNGAM.
! SGNGAM IS EITHER +1.0 OR -1.0.
!
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      REAL ALGAM,SGNGAM,X
!
!  EXTERNAL FUNCTIONS
      REAL ALNGAM
      EXTERNAL ALNGAM
!
!  INTRINSIC FUNCTIONS
      INTRINSIC INT,MOD
!
      ALGAM = ALNGAM(X)
      SGNGAM = 1.0
      IF (X.GT.0.0) RETURN
!
!     INT = AMOD (-AINT(X), 2.0) + 0.1
      IF (INT(MOD(-INT(X),2)+0.1).EQ.0) SGNGAM = -1.0
!
      RETURN
      END
!DGAMI
      DOUBLE PRECISION FUNCTION DGAMI (A, X)
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
      DOUBLE PRECISION A,X
!
!  LOCAL SCALARS
      DOUBLE PRECISION FACTOR
!
!  EXTERNAL FUNCTIONS
      DOUBLE PRECISION D1MACH,DGAMIT,DLNGAM
      EXTERNAL D1MACH,DGAMIT,DLNGAM
!
!  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
!
!  INTRINSIC FUNCTIONS
      INTRINSIC EXP,LOG
!
!
      IF (A.LE.0.D0) CALL XERROR ('DGAMI   A MUST BE GT ZERO', 25, 1,2)
      IF (X.LT.0.D0) CALL XERROR ('DGAMI   X MUST BE GE ZERO', 25, 2,2)
!
      DGAMI = 0.D0
      IF (X.EQ.0.0D0) RETURN
!
! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
!
      FACTOR = DLNGAM(A) + A*LOG(X)
      IF (FACTOR.GT.LOG(D1MACH(2))) THEN
         DGAMI = D1MACH(2)
      ELSE
         DGAMI = EXP(FACTOR) * DGAMIT(A,X)
      END IF
!
      RETURN
      END

!I1MACH
      INTEGER FUNCTION I1MACH(I)
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
!
!  VARIABLE DECLARATIONS
!
!  SCALAR ARGUMENTS
      INTEGER I
!
!  LOCAL SCALARS
      INTEGER OUTPUT,SANITY
!
!  LOCAL ARRAYS
      INTEGER IMACH(16)
!
!  EXTERNAL SUBROUTINES
      EXTERNAL FDUMP
!
!  EQUIVALENCES
      EQUIVALENCE (IMACH(4),OUTPUT)
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
!
!  ***  ISSUE STOP IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) THEN
         STOP 'I1MACH, D1MACH AND R1MACH HAVE NOT BEEN INITIALIZED'
      ELSE
!
         IF (I .LT. 1  .OR.  I .GT. 16) THEN
            WRITE(OUTPUT,9000)
            CALL FDUMP
            STOP
         ELSE
            I1MACH=IMACH(I)
         END IF
      END IF
!
      RETURN
!
 9000 FORMAT('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
!
      END
!R1MACH
      REAL FUNCTION R1MACH(I)
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
      INTEGER I
!
!  LOCAL ARRAYS
      REAL RMACH(5)
      INTEGER DIVER(2),LARGE(2),LOG10(2),RIGHT(2),SMALL(2)
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH
      EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL SETERR
!
!  EQUIVALENCES
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
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
!     CHECK FOR INITIALIZATION WITH DUMMY CALL TO I1MACH
!
      R1MACH = I1MACH(I)
!
      IF (I .LT. 1  .OR.  I .GT. 5)
     +   CALL SETERR('R1MACH - I OUT OF BOUNDS',24,1,2)
!
      R1MACH = RMACH(I)
      RETURN
!
      END
!D1MACH
      DOUBLE PRECISION FUNCTION D1MACH(I)
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
      INTEGER I
!
!  LOCAL ARRAYS
      DOUBLE PRECISION DMACH(5)
      INTEGER DIVER(4),LARGE(4),LOG10(4),RIGHT(4),SMALL(4)
!
!  EXTERNAL FUNCTIONS
      INTEGER I1MACH
      EXTERNAL I1MACH
!
!  EXTERNAL SUBROUTINES
      EXTERNAL SETERR
!
!  EQUIVALENCES
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
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
!
!     CHECK FOR INITIALIZATION WITH DUMMY CALL TO I1MACH
!
      D1MACH = I1MACH(I)
!
!
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL SETERR('D1MACH - I OUT OF BOUNDS',24,1,2)
!
      D1MACH = DMACH(I)
      RETURN
!
      END

      end module M_starpac

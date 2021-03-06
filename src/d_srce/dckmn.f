*DCKMN
      SUBROUTINE DCKMN(J, D, PAR, SCALE, NPAR, ETA, TAU, MDL, XM,
     +   N, NROW, M, IXM, PV, PVTEMP, MSG, LMSG)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN SUBROUTINE FOR CHECKING USER SUPPLIED
C     ANALYTIC DERIVATIVES AGAINST NUMERICAL DERIVATIVES
C
C     WRITTEN BY  -  ROBERT B. SCHNABEL (CODED BY JANET R. DONALDSON)
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   D,ETA,PV,SCALE,TAU
      INTEGER
     +   IXM,J,LMSG,M,N,NPAR,NROW
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(NPAR),PVTEMP(N),XM(IXM,M)
      INTEGER
     +   MSG(LMSG)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL MDL
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FD,PARMX,PVPSTP,STP,TEMP
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DCKCRV,DCKZRO
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX,SIGN,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION D
C        THE SCALAR IN WHICH ROW   NROW   OF THE DERIVATIVE
C        MATRIX WITH RESPECT TO THE JTH UNKNOWN PARAMETER
C        IS STORED.
C     DOUBLE PRECISION ETA
C        THE RELATIVE NOISE IN THE MODEL
C     DOUBLE PRECISION FD
C        THE FORWARD DIFFERENCE QUOTIENT DERIVATIVE WITH RESPECT TO THE
C        JTH PARAMETER
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER J
C        THE INDEX OF THE PARTIAL DERIVATIVE BEING EXAMINED.
C     INTEGER LMSG
C        THE LENGTH OF THE VECTOR MSG.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER MSG(LMSG)
C        AN ARRAY USED TO STORE MESSAGE PARAMETERS.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NROW
C        THE NUMBER OF THE ROW OF THE INDEPENDENT VARIABLE ARRAY AT
C        WHICH THE DERIVATIVE IS TO BE CHECKED.
C     DOUBLE PRECISION PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     DOUBLE PRECISION PARMX
C        THE MAXIMUM OF THE CURRENT PARAMETER ESTIMATE AND THE
C        TYPICAL VALUE OF THAT PARAMETER
C     DOUBLE PRECISION PV
C        THE SCALAR IN WHICH THE PREDICTED VALUE FROM THE MODEL FOR
C        ROW   NROW   IS STORED.
C     DOUBLE PRECISION PVPSTP
C        THE PREDICTED VALUE FOR ROW    NROW   OF THE MODEL
C        BASED ON THE CURRENT PARAMETER ESTIMATES
C        FOR ALL BUT THE JTH PARAMETER VALUE, WHICH IS PAR(J) + STP.
C     DOUBLE PRECISION PVTEMP(N)
C        THE VECTOR OF PREDICTED VALUE FROM THE MODEL.
C     DOUBLE PRECISION SCALE
C        THE TYPICAL SIZE OF THE JTH PARAMETER.
C     DOUBLE PRECISION STP
C        THE STEP SIZE CURRENTLY BEING EXAMINED FOR THE FINITE DIFFERENC
C        DERIVATIVE
C     DOUBLE PRECISION TAU
C        THE AGREEMENT TOLERANCE.
C     DOUBLE PRECISION TEMP
C        A TEMPORARY LOCATION IN WHICH THE CURRENT ESTIMATE OF THE JTH
C        PARAMETER IS STORED.
C     DOUBLE PRECISION XM(IXM,M)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C
C     CALCULATE THE JTH PARTIAL DERIVATIVE USING FORWARD DIFFERENCE
C     QUOTIENTS AND DECIDE IF IT AGREES WITH USER SUPPLIED VALUES
C
      MSG(J+1) = 0
C
      PARMX = MAX(ABS(PAR(J)),ABS(SCALE))
      IF (PARMX .EQ. 0.0D0) PARMX = 1.0D0
C
C     COMPUTE INITIAL STEP SIZE
C
      STP = (SQRT(ETA)*PARMX*SIGN(1.0D0,PAR(J))+PAR(J)) - PAR(J)
C
C     COMPUTE PREDICTED VALUES
C
      TEMP = PAR(J)
      PAR(J) = PAR(J) + STP
      CALL MDL(PAR, NPAR, XM, N, M, IXM, PVTEMP)
      PAR(J) = TEMP
C
      PVPSTP = PVTEMP(NROW)
C
      FD = (PVPSTP-PV)/STP
C
C     CHECK FOR DISAGREEMENT
C
      IF (ABS(FD-D) .GT. TAU*ABS(D)) GO TO 10
C
C     NUMERICAL AND ANALYTIC DERIVATIVES AGREE
C
C     CHECK IF ANALYTIC DERIVATIVE IS IDENTICALLY ZERO, INDICATING
C     THE POSSIBILITY THAT THE DERIVATIVE SHOULD BE RECHECKED AT
C     ANOTHER POINT.
C
      IF (D.NE.0.0D0) RETURN
C
C     JTH ANALYTIC AND NUMERICAL DERIVATIVES BOTH ARE ZERO.
C
      IF (MSG(1).EQ.0) MSG(1) = 1
      MSG(J+1) = 3
      RETURN
C
   10 CONTINUE
C
C     NUMERICAL AND ANALYTIC DERIVATIVES DISAGREE
C
C     CHECK WHY
C
      IF (D.EQ.0.0D0) THEN
         CALL DCKZRO(J, PAR, NPAR, MDL, XM, N,
     +      NROW, M, IXM, PV, PVTEMP, MSG, LMSG, FD, PARMX, PVPSTP,
     +      STP)
      ELSE
         CALL DCKCRV(J, D, PAR, NPAR, ETA, TAU, MDL, XM,
     +      N, NROW, M, IXM, PV, PVTEMP, MSG, LMSG, FD, PARMX,
     +      PVPSTP, STP)
      END IF
C
      RETURN
      END

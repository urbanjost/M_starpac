*NLSX1
      SUBROUTINE NLSX1(MOD, PAR, NPAR, PV, SDPV, RES, SDRES, VCV, N,
     +   IVCV, NNZW, NPARE, RSD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     SET THE STARTING PARAMETER VALUES FOR NLSX
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   RSD
      INTEGER
     +   IVCV,MOD,N,NNZW,NPAR,NPARE
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),PV(N),RES(N),SDPV(N),SDRES(N),VCV(IVCV,IVCV)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,J
C
C  EXTERNAL SUBROUTINES
      EXTERNAL SETRV
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IVCV
C        THE ACTUAL FIRST DIMENSION OF VCV.
C     INTEGER MOD
C        AN INDICATOR VALUE USED TO DESIGNATE THE MODEL FOR WHICH
C        THE PARAMETERS ARE TO BE SET.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C        TO BE PROVIDED.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS ESTIMATED BY THE ROUTINE.
C     INTEGER NNZW
C        THE NUMBER OF NONZERO WEIGHTS.
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     REAL PV(N)
C        THE PREDICTED VALUES.
C     REAL RES(N)
C        THE RESIDUALS.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     REAL SDPV(N)
C        THE STANDARD DEVIATION OF THE PREDICTED VALUES.
C     REAL SDRES(N)
C        THE STANDARDIZED RESIDUALS.
C     REAL VCV(IVCV,IVCV)
C        THE VARIANCE COVARIANCE MATRIX.
C
C
C
      GO TO (10, 20, 30, 40, 50, 60), MOD
C
   10 PAR(1) = 0.725E0
      PAR(2) = 4.0E0
C
      GO TO 70
C
C
   20 PAR(1) = 1.0E0
      PAR(2) = 2.0E0
      PAR(3) = 3.0E0
C
      GO TO 70
C
C
   30 PAR(1) = 6.0E0
      PAR(2) = 5.0E0
      PAR(3) = 4.0E0
      PAR(4) = 3.0E0
      PAR(5) = 2.0E0
C
      GO TO 70
C
C
   40 CALL SETRV(PAR, NPAR, 0.0E0)
C
      GO TO 70
C
C
   50 CALL SETRV(PAR, NPAR, 0.5E0)
C
      GO TO 70
C
C
   60 PAR(1) = 100.0E0
      PAR(2) = 15.0E0
C
   70 CONTINUE
C
      DO 80 I=1,N
         RES(I) = -1.0E0
         PV(I) = -1.0E0
         SDPV(I) = -1.0E0
         SDRES(I) = -1.0E0
   80 CONTINUE
C
      DO 100 I=1,IVCV
         DO 90 J=1,IVCV
            VCV(I,J) = -1.0E0
   90    CONTINUE
  100 CONTINUE
C
      NNZW = -1
      NPARE = -1
      RSD = -1.0E0
C
      IERR = -1
C
      RETURN
C
      END

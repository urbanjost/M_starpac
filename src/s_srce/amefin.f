*AMEFIN
      SUBROUTINE AMEFIN(Y, WEIGHT, NNZW, WT, LWT, XM, N, M, IXM, IFIXD,
     +   PAR, NPAR, NPARE, RES, PAGE, WIDE, IPTOUT, NDIGIT, RSSHLF, RSD,
     +   PVT, SDPVT, SDREST, RD, VCVL, LVCVL, D, AMEHDR, IVCVPT, ISKULL,
     +   NRESTS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPLETES THE ANALYSIS FOR THE NONLINEAR
C     LEAST SQUARES ESTIMATION ROUTINES ONCE THE ESTIMATES
C     HAVE BEEN FOUND.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   RSD,RSSHLF
      INTEGER
     +   IVCVPT,IXM,LVCVL,LWT,M,N,NDIGIT,NNZW,NPAR,NPARE,NRESTS
      LOGICAL
     +   PAGE,WEIGHT,WIDE
C
C  ARRAY ARGUMENTS
      REAL
     +   D(N,*),PAR(*),PVT(*),RD(*),RES(*),SDPVT(*),SDREST(*),VCVL(*),
     +   WT(*),XM(IXM,*),Y(*)
      INTEGER
     +   IFIXD(*),IPTOUT(*),ISKULL(10)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL AMEHDR
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   COND,RSS,YSS
      INTEGER
     +   I,IDF
      LOGICAL
     +   EXACT,PRTFSM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMEOUT,NLCMP
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL COND
C        THE CONDITION NUMBER OF D.
C     REAL D(N,NPAR)
C        THE FIRST DERIVATIVE OF THE MODEL (JACOBIAN).
C     LOGICAL EXACT
C        AN INDICATOR VALUE USED TO DESIGNATE WHETHER THE FIT
C        WAS EXACT TO MACHINE PRECISION (TRUE) OR NOT (FALSE).
C     EXTERNAL AMEHDR
C        THE ROUTINE USED TO PRINT THE HEADING
C     INTEGER IDF
C        THE DEGREES OF FREEDOM IN THE FIT.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IPTOUT(NDIGIT)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER ISKULL(10)
C        AN ERROR MESSAGE INDICATOR VARIABLE.
C     INTEGER IVCVPT
C        AN INDICATOR VALUE USED TO DESIGNATE WHICH FORM OF THE
C        VARIANCE COVARIANCE MATRIX (VCV) IS BEING PRINTED, WHERE
C        IVCVPT = 1 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C        IVCVPT = 2 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(HESSIAN)
C        IVCVPT = 3 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                       *INVERSE(HESSIAN)
C     INTEGER IXM
C        THE FIRST DIMENSION OF MATRIX XM.
C     INTEGER LVCVL
C        THE DIMENSION OF VECTOR VCVL.
C     INTEGER LWT
C        THE DIMENSION OF VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE OPTIMIZED.
C     INTEGER NRESTS
C        THE MAXIMUM NUMBER OF RESIDUALS TO BE COMPUTED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     LOGICAL PRTFSM
C        THE VARIABLE USED TO INDICATE WHETHER ANY OF THE SUMMARY
C        INFORMATION IS TO BE PRINTED (TRUE) OR NOT (FALSE).
C     REAL PVT(N)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES.
C     REAL RD(N)
C        THE DIAGONAL ELEMENTS OF THE R MATRIX OF THE Q - R
C        FACTORIZATION OF D.
C     REAL RES(N)
C        THE RESIDUALS FROM THE FIT.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     REAL RSS
C        THE RESIDUAL SUM OF SQUARES.
C     REAL RSSHLF
C        HALF THE RESIDUAL SUM OF SQUARES.
C     REAL SDPVT(N)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     REAL SDREST(N)
C        THE STANDARDIZED RESIDUALS.
C     REAL VCVL(LVCVL)
C        THE LOWER HALF OF THE VARIANCE-COVARIANCE MATRIX, STORED
C        ROW WISE.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     REAL WT(LWT)
C        THE USER SUPPLIED WEIGHTS.
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE.
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C     REAL YSS
C        THE SUM OF THE SQUARES ABOUT THE MEAN Y VALUE.
C
C     MODIFY VCV TO REFLECT PROPER DEGREES OF FREEDOM
C
      DO 10 I=1,LVCVL
         VCVL(I) = (NRESTS-NPAR)*VCVL(I)/(N-NPAR)
   10 CONTINUE
C
C     COMPUTE RETURNED AND/OR PRINTED VALUES.
C
      CALL NLCMP (Y, WEIGHT, WT, LWT, N, NPAR, NPARE, RES,
     +   D, RD, COND, VCVL, LVCVL, NNZW, IDF, RSSHLF, RSS, RSD, YSS,
     +   EXACT, PVT, SDPVT, SDREST, ISKULL)
C
      PRTFSM = ((IPTOUT(3).NE.0) .OR. (IPTOUT(4).NE.0) .OR.
     +   (IPTOUT(5).NE.0) .OR. (IERR.NE.0))
C
C     PRINT SUMMARY INFORMATION IF DESIRED OR IF AN ERROR FLAG
C     HAS BEEN SET.
C
      IF (PRTFSM) CALL AMEOUT(Y, N,
     +   IFIXD, PAR, NPAR, NPARE, RES, IPTOUT, NDIGIT, PAGE, IDF, COND,
     +   RSS, RSD, YSS, EXACT, PVT, SDPVT, SDREST, VCVL, LVCVL, IVCVPT,
     +   ISKULL, AMEHDR, WIDE)
      RETURN
C
      END

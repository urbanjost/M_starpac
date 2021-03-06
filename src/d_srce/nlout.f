*NLOUT
      SUBROUTINE NLOUT(Y, WEIGHT, NNZW, WT, LWT, XM, N, M, IXM, IFIXD,
     +   PAR, NPAR, NPARE, RES, IPTOUT, NDIGIT, PAGE, IDF, COND, RSS,
     +   RSD, YSS, EXACT, PV, SDPV, SDRES, VCVL, LVCVL, IVCVPT, ISKULL,
     +   NLHDR, WIDE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE FINAL SUMMARY OUTPUT FROM THE
C     NONLINEAR LEAST SQUARES SUBOUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
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
     +   COND,RSD,RSS,YSS
      INTEGER
     +   IDF,IVCVPT,IXM,LVCVL,LWT,M,N,NDIGIT,NNZW,NPAR,NPARE
      LOGICAL
     +   EXACT,PAGE,WEIGHT,WIDE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(NPAR),PV(N),RES(N),SDPV(N),SDRES(N),VCVL(LVCVL),WT(LWT),
     +   XM(IXM,M),Y(N)
      INTEGER
     +   IFIXD(NPAR),IPTOUT(NDIGIT),ISKULL(10)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL NLHDR
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FPLM,PLL,PUL,RATIO,SDPAR,T
      INTEGER
     +   I,IPK,IPK2,IPRT,ISUBHD
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   PPFT,D1MACH
      EXTERNAL PPFT,D1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FITPT1,FITPT2,IPRINT,NLSKL,VCVOTF
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX,SQRT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION COND
C        THE CONDITION NUMBER OF D.
C     LOGICAL EXACT
C        AN INDICATOR VALUE USED TO DESIGNATE WHETHER THE FIT
C        WAS EXACT TO MACHINE PRECISION (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER IDF
C        THE DEGREES OF FREEDOM IN THE FIT.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IPK
C        AN INDEX.
C     INTEGER IPK2
C        THE INDEX OF THE "DIAGONAL" ELEMENT OF THE VCV
C        MATRIX.
C     INTEGER IPRT
C        THE LOGICAL UNIT FOR PRINTED OUTPUT.
C     INTEGER IPTOUT(NDIGIT)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER ISKULL(10)
C        AN ERROR MESSAGE INDICATOR VARIABLE.
C     INTEGER ISUBHD
C        AN INTEGER VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
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
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER LVCVL
C        THE LENGTH OF THE VECTOR CONTANING
C        THE LOWER HALF OF THE VCV MATRIX, STORED ROW WISE.
C     INTEGER LWT
C        THE ACTUAL LENGTH OF THE VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     EXTERNAL NLHDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE OPTIMIZED.
C     DOUBLE PRECISION PAR(NPAR)
C        THE PARAMETER ESTIMATES.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER OR NOT THE OUTPUT
C        IS TO BEGIN ON A NEW PAGE.
C     DOUBLE PRECISION PLL, PUL
C        THE LOWER AND UPPER CONFIDENCE LIMITS FOR A GIVEN PARAMETER.
C     DOUBLE PRECISION PV(N)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES
C     DOUBLE PRECISION RATIO
C        THE RATIO OF A GIVEN PARAMETER VALUE TO ITS STANDARD ERROR.
C     DOUBLE PRECISION RES(N)
C        THE RESIDUALS FROM THE FIT.
C     DOUBLE PRECISION RSD
C        THE VALUE OF THE RESIDUAL STANDARD DEVIATION AT THE SOLUTION.
C     DOUBLE PRECISION RSS
C        THE RESIDUAL SUM OF SQUARES.
C     DOUBLE PRECISION SDPAR
C        THE STANDARD DEVIATION OF A GIVEN PARAMETER VALUE.
C     DOUBLE PRECISION SDPV(N)
C        THE STANDARD DEVIATION OF THE PREDICTED VALUE.
C     DOUBLE PRECISION SDRES(N)
C        THE STANDARD DEVIATIONS OF THE RESIDUALS.
C     DOUBLE PRECISION T
C        THE VALUE OF THE 97.5 PERCENT POINT FUNCTION FOR THE
C        T DISTRIBUTION.
C     DOUBLE PRECISION VCVL(LVCVL)
C        THE LOWER HALF OF THE VCV MATRIX, STORED ROW WISE.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION WT(LWT)
C        THE USER SUPPLIED WEIGHTS.
C     DOUBLE PRECISION XM(IXM,M)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C     DOUBLE PRECISION YSS
C        THE SUM OF SQUARES OF THE DEPENDENT VARIABLE Y.
C
C
      FPLM = D1MACH(2)
C
      CALL IPRINT(IPRT)
C
      IF ((IERR.GE.1) .AND. (IERR.NE.4)) GO TO 60
C
C     TEST FOR EXACT FIT
C
      IF ((IDF.LE.0) .OR. EXACT) GO TO 70
C
C     PRINT ERROR HEADING IF NECESSARY
C
      IF (IERR.EQ.4) CALL NLSKL(ISKULL, PAGE, WIDE, NLHDR)
C
C     PRINT PRIMARY REPORT
C
      IF ((IERR.EQ.0) .AND. (IPTOUT(3).EQ.0)) GO TO 10
      ISUBHD = 0
      CALL NLHDR(PAGE, WIDE, ISUBHD)
      CALL FITPT1(N, M, XM, Y, PV, SDPV, RES, SDRES, WT, IXM, NNZW,
     +   WEIGHT, MAX(1, IPTOUT(3)))
C
C     PRINT STANDARDIZED RESIDUAL PLOTS
C
   10 IF (IPTOUT(4).EQ.0) GO TO 20
      ISUBHD = 0
      CALL NLHDR(PAGE, WIDE, ISUBHD)
C
      CALL FITPT2 (SDRES, PV, WT, N, NNZW, WEIGHT, RES, RSS)
C
C     PRINT THE COVARIANCE AND CORRELATION MATRIX
C
   20 IF ((IERR.EQ.0) .AND. (IPTOUT(5).EQ.0)) RETURN
      ISUBHD = 0
      CALL NLHDR(PAGE, WIDE, ISUBHD)
C
      IF ((IERR.EQ.0) .AND. (IPTOUT(5).LE.1)) GO TO 30
C
      CALL VCVOTF(NPARE, VCVL, LVCVL, .TRUE., NPAR, IFIXD, IVCVPT)
C
C     PRINT ANALYSIS SUMMARY
C
   30 WRITE (IPRT,1000)
      WRITE (IPRT,1010)
C
      T = PPFT(0.975D0,IDF)
      IPK = 0
      DO 50 I=1,NPAR
C
         IF (IFIXD(I).EQ.0) GO TO 40
C
         WRITE (IPRT,1030) I, PAR(I)
         GO TO 50
C
   40    IPK = IPK + 1
         IPK2 = IPK*(IPK-1)/2 + IPK
         RATIO = FPLM
         SDPAR = SQRT(VCVL(IPK2))
         IF (SDPAR.GT.0.0D0) RATIO = PAR(I)/SDPAR
         PLL = PAR(I) - T*SDPAR
         PUL = PAR(I) + T*SDPAR
         WRITE (IPRT,1020) I, PAR(I), SDPAR, RATIO, PLL, PUL
   50 CONTINUE
      WRITE (IPRT,1040) RSS, RSD, NNZW, NPARE, IDF
      WRITE (IPRT,1050) COND
C
      IF (RSS.GT.YSS) WRITE (IPRT,1060)
C
      RETURN
C
C     PRINT OUT ERROR HEADING
C
   60 CALL NLSKL(ISKULL, PAGE, WIDE, NLHDR)
C
      IF (IERR.LE.2) RETURN
C
C     PRINT SECONDARY REPORT
C
   70 CONTINUE
      ISUBHD = 0
      CALL NLHDR(PAGE, WIDE, ISUBHD)
      IF (IERR.NE.0) WRITE (IPRT,1080)
      WRITE (IPRT,1000)
      WRITE (IPRT,1100)
      DO 80 I=1,NPAR
         IF (IFIXD(I).EQ.0) WRITE (IPRT,1110) I, PAR(I)
         IF (IFIXD(I).NE.0) WRITE (IPRT,1120) I, PAR(I)
   80 CONTINUE
      WRITE (IPRT,1040) RSS, RSD, NNZW, NPARE, IDF
      IF (IERR.NE.3) WRITE (IPRT,1050) COND
C
      IF ((IERR.EQ.0) .AND. (.NOT.EXACT) .AND. (IDF.LE.0)) WRITE
     +   (IPRT,1070)
      IF ((IERR.EQ.0) .AND. EXACT) WRITE (IPRT,1090)
C
      IF (IERR.NE.0) GO TO 100
C
      DO 90 I=1,N
         SDRES(I) = 0.0D0
         SDPV(I) = 0.0D0
   90 CONTINUE
C
      RETURN
C
  100 CONTINUE
C
      DO 110 I=1,N
         SDRES(I) = FPLM
         SDPV(I) = FPLM
  110 CONTINUE
C
C     PRINT OUT ERROR EXIT STATISTICS
C
      CALL FITPT1(N, M, XM, Y, PV, SDPV, RES, SDRES, WT, IXM, NNZW,
     +   WEIGHT, MAX(IPTOUT(3),1))
C
C
C     WIPE OUT SDRES VECTOR
C
      DO 120 I=1,N
         SDRES(I) = FPLM
  120 CONTINUE
C
C     WIPE OUT VCV MATRIX
C
      DO 140 I=1,LVCVL
         VCVL(I) = FPLM
  140 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (///33H ESTIMATES FROM LEAST SQUARES FIT/1X, 33('-'))
 1010 FORMAT (//69X, 11HAPPROXIMATE/61X, 28H95 PERCENT CONFIDENCE LIMITS
     +   /1X, 5HINDEX, 2X, 5HFIXED, 3X, 9HPARAMETER, 8X, 9HSD OF PAR,
     +   7X, 5HRATIO, 12X, 5HLOWER, 12X, 5HUPPER/)
 1020 FORMAT (1X, I3, 5X, 3H NO, 2G17.8, 2X, G10.4, 2G17.8)
 1030 FORMAT (1X, I3, 5X, 3HYES, G17.8, 10X, 3H---, 8X, 3H---, 15X,
     +   3H---, 14X, 3H---)
 1040 FORMAT (//31H RESIDUAL SUM OF SQUARES       , 8X, G15.7//6H RESID,
     +   25HUAL STANDARD DEVIATION   , 8X, G15.7/19H BASED ON DEGREES O,
     +   9HF FREEDOM, 5X, I4, 3H - , I3, 3H = , I4)
 1050 FORMAT (/29H APPROXIMATE CONDITION NUMBER, 10X, G15.7)
 1060 FORMAT (52H THE RESIDUAL SUM OF SQUARES AFTER THE LEAST SQUARES,
     +   20H FIT IS GREATER THAN/35H THE SUM OF SQUARES ABOUT THE MEAN ,
     +   19HY OBSERVATION.  THE, 14H MODEL IS LESS/17H REPRESENTATIVE O,
     +   39HF THE DATA THAN A SIMPLE AVERAGE.  DATA, 15H AND MODEL SHOU,
     +   2HLD/48H BE CHECKED TO BE SURE THAT THEY ARE COMPATABLE.)
 1070 FORMAT (/49H THE DEGREES OF FREEDOM FOR THIS PROBLEM IS ZERO.,
     +   54H  STATISTICAL ANALYSIS OF THE RESULTS IS NOT POSSIBLE.)
 1080 FORMAT (//40H THE FOLLOWING SUMMARY SHOULD BE USED TO, 8H ANALYZE,
     +   30H THE ABOVE MENTIONED PROBLEMS.)
 1090 FORMAT (/50H THE LEAST SQUARES FIT OF THE DATA TO THE MODEL IS,
     +   35H EXACT TO WITHIN MACHINE PRECISION./20H STATISTICAL ANALYSI,
     +   33HS OF THE RESULTS IS NOT POSSIBLE.)
 1100 FORMAT (//1X, 5HINDEX, 2X, 5HFIXED, 3X, 9HPARAMETER)
 1110 FORMAT (1X, I3, 5X, 3H NO, G17.8)
 1120 FORMAT (1X, I3, 5X, 3HYES, G17.8)
      END

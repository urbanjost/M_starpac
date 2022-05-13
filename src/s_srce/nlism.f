*NLISM
      SUBROUTINE NLISM (NLHDR, PAGE, WIDE, HLFRPT, NPAR, M, N, NNZW,
     +   WEIGHT, IFIXD, PAR, SCALE, IWORK, IIWORK, RWORK, IRWORK, RES,
     +   APRXDV, STP, LSTP, NPARE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS AN INITIAL SUMMARY OF THE STARTING
C     ESTIMATES AND THE CONTROL PARAMETERS FOR THE NONLINEAR
C     LEAST SQUARES SUBROUTINES.
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
      INTEGER
     +   IIWORK,IRWORK,LSTP,M,N,NNZW,NPAR,NPARE
      LOGICAL
     +   APRXDV,HLFRPT,PAGE,WEIGHT,WIDE
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),RES(N),RWORK(IRWORK),SCALE(NPAR),STP(LSTP)
      INTEGER
     +   IFIXD(NPAR),IWORK(IIWORK)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL NLHDR
C
C  LOCAL SCALARS
      REAL
     +   RSD,RSS
      INTEGER
     +   I,IPRT,ISUBHD,LMAX0,MXFCAL,MXITER,RFCTOL,XCTOL
C
C  EXTERNAL FUNCTIONS
      REAL
     +   SNRM2
      EXTERNAL SNRM2
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC REAL,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL APRXDV
C        THE VARIABLE USED TO INDICATE WHETHER NUMERICAL
C        APPROXIMATIONS TO THE DERIVATIVE WERE USED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HLFRPT
C        THE VARIABLE WHICH INDICATES WHETHER THE DERIVATIVE
C        CHECKING ROUTINE HAS ALREADY PRINTED PART OF THE
C        INITIAL SUMMARY (TRUE) OR NOT (FALSE).
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IIWORK
C        THE DIMENSION OF THE INTEGER WORK VECTOR IWORK.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IRWORK
C        THE DIMENSION OF THE REAL WORK VECTOR RWORK.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     INTEGER IWORK(IIWORK)
C        THE INTEGER WORK SPACE VECTOR USED BY THE NL2 SUBROUTINES.
C     INTEGER LMAX0
C        THE LOCATION IN RWORK OF THE VALUE INDICATING THE
C        MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     INTEGER LSTP
C        THE ACTUAL LENGTH OF THE VECTOR STP.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MXFCAL
C        THE LOCATION IN IWORK OF THE VARIABLE DESIGNATING THE
C        MAXIMUM NUMBER OF FUNCTION CALLS ALLOWED, EXCLUDING
C        CALLS NECESSARY TO COMPUTE THE DERIVATIVES AND COVARIANCE
C        MATRIX.
C     INTEGER MXITER
C        THE LOCATION IN IWORK OF THE VARIABLE DESIGNATING THE
C        MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     EXTERNAL NLHDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        NUMBER OF PARAMETERS ESTIMATED BY ROUTINE.
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     REAL RES(N)
C        THE RESIDUALS FROM THE FIT.
C     INTEGER RFCTOL
C        THE LOCATION IN RWORK OF THE RELATIVE FUNCTION CONVERGENCE
C        TOLERANCE.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION AT THE INITIAL PARAMETER VALUES
C     REAL RSS
C        THE RESIDUAL SUM OF SQUARES AT THE INITIAL PARAMETER VALUES
C     REAL RWORK(IRWORK)
C        THE REAL WORK VECTOR USED BY THE NL2 SUBROUTINES.
C     REAL SCALE(NPAR)
C        THE TYPICAL SIZE OF THE UNKNOWN PARAMETERS.
C     REAL STP(LSTP)
C        THE SELECTED RELATIVE STEP SIZES.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C       THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        FULL WIDTH (TRUE) OR NOT (FALSE).
C     INTEGER XCTOL
C        THE LOCATION IN RSTAK/DSTAK OF THE P CONVERGENCE TOLERANCE.
C
C     IWORK SUBSCRIPT VALUES
C
      DATA MXFCAL/17/, MXITER/18/
C
C     RWORK SUBSCRIPT VALUES
C
      DATA LMAX0/35/, RFCTOL/32/, XCTOL/33/
C
      CALL IPRINT(IPRT)
C
      IF (.NOT.HLFRPT) THEN
C
         ISUBHD = 1
         CALL NLHDR(PAGE, WIDE, ISUBHD)
C
         IF (APRXDV) THEN
            WRITE (IPRT, 1260)
         ELSE
            WRITE (IPRT, 1270)
         END IF
C
         DO 40 I = 1, NPAR
            IF (IFIXD(I).EQ.0) THEN
               IF (SCALE(1).GT.0.0E0) THEN
                  IF (APRXDV) THEN
                     WRITE (IPRT, 1410) I, PAR(I), SCALE(I), STP(I)
                  ELSE
                     WRITE (IPRT, 1410) I, PAR(I), SCALE(I)
                  END IF
               ELSE
                  IF (APRXDV) THEN
                     WRITE (IPRT, 1310) I, PAR(I), STP(I)
                  ELSE
                     WRITE (IPRT, 1310) I, PAR(I)
                  END IF
               END IF
            ELSE
               IF (APRXDV) THEN
                  WRITE (IPRT, 1510) I, PAR(I)
               ELSE
                  WRITE (IPRT, 1610) I, PAR(I)
               END IF
            END IF
   40    CONTINUE
C
         WRITE (IPRT, 1160) N
C
      END IF
C
      IF (WEIGHT) WRITE (IPRT, 1170) NNZW
      WRITE (IPRT, 1180) M
      WRITE(IPRT, 1070) IWORK(MXITER)
      WRITE(IPRT, 1090) IWORK(MXFCAL)
      WRITE(IPRT, 1080)
      WRITE(IPRT, 1100) RWORK(RFCTOL)
      WRITE(IPRT, 1110) RWORK(XCTOL)
      WRITE(IPRT, 1120) RWORK(LMAX0)
C
      RSD = SNRM2(N, RES, 1)
      RSS = RSD * RSD
      IF (NNZW-NPARE.GE.1) RSD = RSD /SQRT(REAL(NNZW-NPARE))
      WRITE (IPRT, 1200) RSS
      WRITE (IPRT, 1210) RSD
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1070 FORMAT (/37H MAXIMUM NUMBER OF ITERATIONS ALLOWED, 32X, 5H(MIT),
     +   1X, I5)
 1080 FORMAT(/44H CONVERGENCE CRITERION FOR TEST BASED ON THE/)
 1090 FORMAT(/41H MAXIMUM NUMBER OF MODEL SUBROUTINE CALLS,
     +   8H ALLOWED, 26X, I5)
 1100 FORMAT (5X, 39H FORECASTED RELATIVE CHANGE IN RESIDUAL,
     +   15H SUM OF SQUARES, 7X, 8H(STOPSS), 1X, G11.4)
 1110 FORMAT(5X, 49H MAXIMUM SCALED RELATIVE CHANGE IN THE PARAMETERS,
     +   13X, 7H(STOPP), 1X, G11.4)
 1120 FORMAT(//41H MAXIMUM CHANGE ALLOWED IN THE PARAMETERS,
     +   23H AT THE FIRST ITERATION, 3X, 7H(DELTA), 1X, G11.4)
 1160 FORMAT (/23H NUMBER OF OBSERVATIONS, 48X, 3H(N), 1X, I5)
 1170 FORMAT (/41H NUMBER OF NON ZERO WEIGHTED OBSERVATIONS, 27X,
     +   6H(NNZW), 1X, I5)
 1180 FORMAT (/32H NUMBER OF INDEPENDENT VARIABLES, 39X, 3H(M), 1X, I5)
 1200 FORMAT (/44H RESIDUAL SUM OF SQUARES FOR INPUT PARAMETER,
     +   7H VALUES, 24X, G11.4)
 1210 FORMAT (/48H RESIDUAL STANDARD DEVIATION FOR INPUT PARAMETER,
     +   7H VALUES, 14X, 5H(RSD), 1X, G11.4)
 1260 FORMAT (//50X, 13HSTEP SIZE FOR/
     +   50X, 13HAPPROXIMATING/
     +   7X, 24HPARAMETER STARTING VALUE, 6X, 5HSCALE, 10X,
     +   10HDERIVATIVE/
     +   1X, 5HINDEX, 2X, 5HFIXED, 6X, 5H(PAR), 12X, 7H(SCALE), 11X,
     +   5H(STP)/)
 1270 FORMAT (//6X, 25HPARAMETER STARTING VALUES, 5X, 5HSCALE/
     +   1X, 5HINDEX, 2X, 5HFIXED, 6X, 5H(PAR), 11X, 7H(SCALE)/)
 1310 FORMAT (1X, I3, 5X, ' NO', G17.8, 7X, 7HDEFAULT, 3X, G17.8)
 1410 FORMAT (1X, I3, 5X, ' NO', 3G17.8)
 1510 FORMAT (1X, I3, 5X, 'YES', G17.8, 9X, '---', 14X, '---')
 1610 FORMAT (1X, I3, 5X, 'YES', G17.8, 9X, '---')
      END

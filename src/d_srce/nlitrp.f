*NLITRP
      SUBROUTINE NLITRP(NLHDR, HEAD, PAGE, WIDE, IPTOUT, NPAR, NNZW,
     +   IWORK, IIWORK, RWORK, IRWORK, IFIXD, PARE, NPARE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE ITERATION REPORTS FOR THE
C     NONLINEAR LEAST SQUARES REGRESSION SUBROUTINES.
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
     +   IIWORK,IPTOUT,IRWORK,NNZW,NPAR,NPARE
      LOGICAL
     +   HEAD,PAGE,WIDE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PARE(NPAR),RWORK(IRWORK)
      INTEGER
     +   IFIXD(NPAR),IWORK(IIWORK)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL NLHDR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   RSD,RSS,RSSC,RSSPC
      INTEGER
     +   DST0,F,F0,FDIF,ICASE,IPRT,ISUBHD,MXITER,NFCALL,NITER,
     +   NREDUC,PREDUC,RELDX,STPPAR
      CHARACTER
     +   LETTRN*1,LETTRY*1
C
C  LOCAL ARRAYS
      CHARACTER
     +   ISCHKD(2)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,LSTVCF
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MOD,DBLE,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER DST0
C        THE LOCATION IN RWORK OF THE VALUE OF THE 2 NORM OF D TIMES
C        THE  NEWTON STEP.
C     INTEGER F
C        THE LOCATION IN RWORK OF THE VALUE OF HALF THE RESIDUAL
C        SUM OF SQUARES AT THE CURRENT PARAMETER VALUES.
C     INTEGER FDIF
C        THE LOCATION IN RWORK OF THE DIFFERENCE BETWEEN THE
C        RESIDUAL SUM OF SQUARES AT THE BEGINNING AND END OF THE
C        CURRENT ITERATION.
C     INTEGER F0
C        THE LOCATION IN RWORK OF THE VALUE OF HALF THE RESIDUAL
C        VARIANCE AT THE BEGINNING OF THE CURRENT ITERATION.
C     LOGICAL HEAD
C        THE VARIABLE USED TO INDICATE WHETHER A HEADING IS TO BE
C        PRINTED DURING A GIVEN CALL TO THE ITERATION REPORT (TRUE)
C        OR NOT (FALSE).
C     INTEGER ICASE
C        AN INDICATER VARIABLE USED TO DESIGNATE THE MESSAGE TO BE
C        PRINTED.
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
C        THE DIMENSION OF THE DOUBLE PRECISION WORK VECTOR RWORK.
C     CHARACTER*1 ISCHKD(2)
C        THE INDICATOR USED TO DESIGNATE WHETHER THE
C        TEST VALUE WAS CHECKED FOR CONVERGENCE (Y) OR NOT (N).
C     INTEGER ISUBHD
C        AN INTEGER VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     INTEGER IWORK(IIWORK)
C        THE INTEGER WORK SPACE VECTOR USED BY THE NL2 SUBROUTINES.
C     CHARACTER*1 LETTRN, LETTRY
C        THE LETTERS N AND Y, RESPECTIVELY.
C     INTEGER MXITER
C        THE LOCATION IN IWORK OF THE VARIABLE DESIGNATING THE
C        MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER NFCALL
C        THE LOCATION IN IWORK OF THE NUMBER OF FUNCTION EVALUATIONS.
C     INTEGER NITER
C        THE LOCATION IN IWORK OF THE NUMBER OF THE CURRENT ITERATION.
C     EXTERNAL NLHDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF UNKNOWN PARAMETERS TO BE OPTIMIZED.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NREDUC
C        THE LOCATION IN RWORK OF THE VALUE USED TO CHECK IF THE
C        HESSIAN APPROXIMATION IS POSITIVE DEFINITE.  IF
C        IF RWORK(NREDUC) .EQ. 0, THE HESSIAN IS SINGULAR, OTHERWISE
C        IT IS NOT.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION PARE(NPAR)
C        THE CURRENT ESTIMATES OF THE UNKNOWN PARAMETERS, BUT ONLY
C        THOSE TO BE OPTIMIZED (NOT THOSE WHOSE VALUES ARE FIXED).
C     INTEGER PREDUC
C        THE LOCATION IN RWORK OF THE PREDICTED FUNCTION REDUCTION
C        FOR THE CURRENT STEP.
C     INTEGER RELDX
C        THE LOCATION IN RWORK OF THE SCALED RELATIVE CHANGE IN
C        THE PARAMETER VALUES CAUSED BY THE CURRENT ITERATION.
C     DOUBLE PRECISION RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     DOUBLE PRECISION RSS
C        THE RESIDUAL SUM OF SQUARES.
C     DOUBLE PRECISION RSSC
C        THE CHANGE IN THE RESIDUAL SUM OF SQUARES CAUSED BY THIS
C        ITERATION.
C     DOUBLE PRECISION RSSPC
C        THE PREDICTED CHANGE IN THE RESIDUAL SUM OF SQUARES AT THIS
C        ITERATION.
C     DOUBLE PRECISION RWORK(IRWORK)
C        THE DOUBLE PRECISION WORK VECTOR USED BY THE NL2 SUBROUTINES.
C     INTEGER STPPAR
C        THE LOCATION IN RWORK OF THE MARQUARDT LAMBDA PARAMETER.
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C
C
      DATA LETTRN /'N'/, LETTRY /'Y'/
C
C     IWORK SUBSCRIPT VALUES
C
      DATA MXITER /18/, NFCALL /6/, NITER /31/
C
C     RWORK SUBSCRIPT VALUES
C
      DATA DST0 /3/, F /10/, FDIF /11/, F0 /13/, NREDUC /6/, PREDUC
     +   /7/, RELDX /17/, STPPAR /5/
C
      CALL IPRINT(IPRT)
C
      IF (IWORK(1).EQ.10) GO TO 90
      IF ((IPTOUT.EQ.1) .AND. (IWORK(NITER).NE.1) .AND.
     +   (IWORK(NITER).NE.IWORK(MXITER)) .AND. (IWORK(1).LE.2)) RETURN
C
      ISUBHD = 0
      IF (HEAD) CALL NLHDR(PAGE, WIDE, ISUBHD)
      HEAD = .FALSE.
      IF (MOD(IWORK(NITER),4).EQ.0) HEAD = .TRUE.
C
      WRITE (IPRT,1000) IWORK(NITER)
C
C     COMPUTE STATISTICS TO BE PRINTED
C
      RSS = 2.0D0*RWORK(F)
      RSD = SQRT(RSS)
      IF (NNZW-NPARE.GE.1) RSD = RSD/SQRT(DBLE(NNZW-NPARE))
C
      RSSC = 0.0D0
      IF (RWORK(F0).GT.0.0D0) RSSC = RWORK(FDIF)/RWORK(F0)
C
      RSSPC = 0.0D0
      IF (RWORK(F0).GT.0.0D0) RSSPC = RWORK(NREDUC)/RWORK(F0)
C
C     REFERENCE NL2 SUBROUTINE ASSESS, STATEMENT LABEL 300 TO 320
C
      ISCHKD(1) = LETTRN
      ISCHKD(2) = LETTRN
      IF (RWORK(FDIF).GT.2.0D0*RWORK(PREDUC)) GO TO 10
      IF (RWORK(DST0).LT.0.0D0) GO TO 10
      IF (RWORK(NREDUC).GE.0.0D0) ISCHKD(1) = LETTRY
      IF (RWORK(STPPAR).EQ.0.0D0) ISCHKD(2) = LETTRY
   10 CONTINUE
C
      WRITE (IPRT,1010) IWORK(NFCALL), RSD, RSS, RSSC, RSSPC,
     +   ISCHKD(1), RWORK(RELDX), ISCHKD(2)
      IF (NPARE.LT.NPAR) WRITE (IPRT,1020)
      IF (NPARE.GE.NPAR) WRITE (IPRT,1150)
      CALL LSTVCF(NPARE, PARE, NPAR, IFIXD)
C
      IF (IWORK(1).LE.2) RETURN
C
C     PRINT FINAL ITERATION MESSAGE
C
      ICASE = IWORK(1) - 2
      GO TO (20, 30, 40, 50, 60, 70, 80, 90, 100, 140, 110, 120, 130),
     +   ICASE
C
C     ***** PARAMETER CONVERGENCE *****
C
   20 WRITE (IPRT,1030)
      RETURN
C
C     ***** RESIDUAL SUM OF SQUARES CONVERGENCE *****
C
   30 WRITE (IPRT,1040)
      RETURN
C
C     ***** PARAMETER AND RESIDUAL SUM OF SQUARES CONVERGENCE ****
C
   40 WRITE (IPRT,1050)
      RETURN
C
C     ***** RESIDUAL SUM OF SQUARES IS EXACTLY ZERO *****
C
   50 WRITE (IPRT,1060)
      RETURN
C
C     ***** SINGULAR CONVERGENCE *****
C
   60 WRITE (IPRT,1070)
      RETURN
C
C     ***** FALSE CONVERGENCE *****
C
   70 WRITE (IPRT,1080)
      RETURN
C
C     ***** LIMIT ON NUM. OF CALLS TO THE MODEL SUBROUTINE REACHED *****
C
   80 WRITE (IPRT,1090)
      RETURN
C
C     ***** ITERATION LIMIT REACHED *****
C
   90 WRITE (IPRT,1100)
      RETURN
C
C     ***** STOPX *****
C
  100 WRITE (IPRT,1110)
      RETURN
C
C     ***** INITIAL RESIDUAL SUM OF SQUARES OVERFLOWS *****
C
  110 WRITE (IPRT,1120)
      RETURN
C
C     ***** BAD PARAMETERS TO ASSESS *****
C
  120 WRITE (IPRT,1130)
      RETURN
C
C     ***** J COULD NOT BE COMPUTED *****
C
  130 WRITE (IPRT,1140)
      RETURN
C
  140 RETURN
C
C      FORMAT STATEMENTS
C
 1000 FORMAT (//17H ITERATION NUMBER, I5/1X, 22('-'))
 1010 FORMAT (5X, 5HMODEL, 53X, 10HFORECASTED/5X, 5HCALLS, 9X, 3HRSD,
     +   13X, 3HRSS, 8X, 12HREL CHNG RSS, 4X, 12HREL CHNG RSS, 4X,
     +   12HREL CHNG PAR/62X, 5HVALUE, 3X, 4HCHKD, 4X, 5HVALUE, 3X,
     +   4HCHKD/3X, I7, 3(2X, G14.4), 2(G12.4, 3X, A1))
 1020 FORMAT (/5X, 25H CURRENT PARAMETER VALUES, 19H (ONLY UNFIXED PARA,
     +   18HMETERS ARE LISTED))
 1030 FORMAT (/34H ***** PARAMETER CONVERGENCE *****)
 1040 FORMAT (/48H ***** RESIDUAL SUM OF SQUARES CONVERGENCE *****)
 1050 FORMAT (/44H ***** PARAMETER AND RESIDUAL SUM OF SQUARES,
     +   18H CONVERGENCE *****)
 1060 FORMAT (/50H ***** THE RESIDUAL SUM OF SQUARES IS EXACTLY ZERO,
     +   6H *****)
 1070 FORMAT (/33H ***** SINGULAR CONVERGENCE *****)
 1080 FORMAT (/30H ***** FALSE CONVERGENCE *****)
 1090 FORMAT (/44H ***** LIMIT ON NUMBER OF CALLS TO THE MODEL,
     +   25H SUBROUTINE REACHED *****)
 1100 FORMAT (/36H ***** ITERATION LIMIT REACHED *****)
 1110 FORMAT (/18H ***** STOPX *****)
 1120 FORMAT (/53H ***** INITIAL RESIDUAL SUM OF SQUARES OVERFLOWS ****,
     +   1H*)
 1130 FORMAT (/37H ***** BAD PARAMETERS TO ASSESS *****)
 1140 FORMAT (/52H ***** DERIVATIVE MATRIX COULD NOT BE COMPUTED *****)
 1150 FORMAT (/5X, 25H CURRENT PARAMETER VALUES)
      END

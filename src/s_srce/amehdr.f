*AMEHDR
      SUBROUTINE AMEHDR(PAGE, WIDE, ISUBHD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE NONLINEAR
C     LEAST SQUARES ESTIMATION ROUTINES FOR ARIMA MODELS THAT USE
C     NUMERICAL APPROXIMATIONS TO THE DERIVATIVES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  AUGUST 1, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   ISUBHD
      LOGICAL
     +   PAGE,WIDE
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,VERSP
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C
      CALL IPRINT(IPRT)
      IF (PAGE) WRITE (IPRT, 1020)
      CALL VERSP(WIDE)
      IF (PAGE) WRITE (IPRT,1000)
      IF (.NOT.PAGE) WRITE (IPRT,1010)
      PAGE = .TRUE.
C
      IF (ISUBHD.EQ.0) RETURN
C
      GO TO (10), ISUBHD
C
   10 WRITE (IPRT, 1030)
C
      RETURN
C
C     FORMAT STATEMENTS FOR PAGE HEADINGS
C
 1000 FORMAT ('+NONLINEAR LEAST SQUARES ESTIMATION',
     +   ' FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED')
 1010 FORMAT ('+', 77(1H*)/
     +   1X, 37H*  NONLINEAR LEAST SQUARES ESTIMATION,
     +   40H FOR THE PARAMETERS OF AN ARIMA MODEL  */
     +   2H *, 16X, 45H             USING BACKFORECASTS             ,
     +   14X, 1H*/1X, 77(1H*))
 1020 FORMAT ('1')
 1030 FORMAT (//30H SUMMARY OF INITIAL CONDITIONS/ 1X, 30('-'))
      END

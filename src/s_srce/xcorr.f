*XCORR
      SUBROUTINE XCORR(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE EXERCISES ALL ASPECTS OF THE CORRELATION
C     FAMILY ROUTINES
C
C     WRITTEN BY -
C        LINDA MITCHELL
C        STATISTICAL ENGINEERING DIVISION
C        NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LDSTAK
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT,IVCV,IYM,J,LDSMIN,M,N
C
C  LOCAL ARRAYS
      REAL
     +   VCV(4,4),YM(10,4),Z(10,4)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CORR,CORRS,CORRXP,GENR,IPRINT,LDSCMP,MSGX,SETRA
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C         AN INDEX VARIABLE.
C     INTEGER IERR
C        COMMON FLAG INDICATING IF ANY ERRORS WERE DETECTED
C        IF IERR = 0, THEN NO ERRORS WERE FOUND
C     INTEGER IPRT
C        THE OUTPUT LOGICAL UNIT NUMBER
C     INTEGER IVCV
C        THE ROW DIMENSION OF VCV
C     INTEGER IYM
C        THE ROW DIMENSION OF YM
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER LDSMIN
C        THE SMALLEST ACCEPTABLE SIZE OF COMMON AREA CSTAK
C     INTEGER LDSTAK
C        THE SIZE OF THE COMMON AREA CSTAK
C     INTEGER M
C        THE NUMBER OF VARIABLES
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS
C     REAL VCV(4,4)
C        THE VARIANCE COVARIANCE MATRIX
C     REAL YM(10,4)
C        GENERAL DATA SET, FROM DRAPER AND SMITH
C     REAL Z(10,4)
C        TEST OBSERVATION MATRIX
C
C     THIS DATA SET TAKEN FROM
C        APPLIED REGRESSION ANALYSIS
C        DRAPER AND SMITH
C        PAGE 216
C
      DATA     YM(1,1),   YM(1,2),   YM(1,3),   YM(1,4)
     +    /      42.2E0,  11.2E0,  31.9E0, 167.1E0/
      DATA     YM(2,1),   YM(2,2),   YM(2,3),   YM(2,4)
     +    /      48.6E0,  10.6E0,  13.2E0, 174.4E0/
      DATA     YM(3,1),   YM(3,2),   YM(3,3),   YM(3,4)
     +    /      42.6E0,  10.6E0,  28.7E0, 160.8E0/
      DATA     YM(4,1),   YM(4,2),   YM(4,3),   YM(4,4)
     +    /      39.0E0,  10.4E0,  26.1E0, 162.0E0/
      DATA     YM(5,1),   YM(5,2),   YM(5,3),   YM(5,4)
     +    /      34.7E0,   9.3E0,  30.1E0, 140.8E0/
      DATA     YM(6,1),   YM(6,2),   YM(6,3),   YM(6,4)
     +    /      44.5E0,  10.8E0,   8.5E0, 174.6E0/
      DATA     YM(7,1),   YM(7,2),   YM(7,3),   YM(7,4)
     +    /      39.1E0,  10.7E0,  24.3E0, 163.7E0/
      DATA     YM(8,1),   YM(8,2),   YM(8,3),   YM(8,4)
     +    /      40.1E0,  10.0E0,  18.6E0, 174.5E0/
      DATA     YM(9,1),   YM(9,2),   YM(9,3),   YM(9,4)
     +    /      45.9E0,  12.0E0,  20.4E0, 185.7E0/
C
C     DETERMINE OUTPUT UNIT
C
      CALL IPRINT(IPRT)
C
      IVCV = 4
      IYM = 10
      M = 4
      N = 9
      IERR = 0
C
C**** TEST ROUTINES WITH CORRECT CALL STATEMENT *****
C
      WRITE (IPRT,1000)
      WRITE (IPRT,1010)
C
C     TEST CORR
C
      WRITE (IPRT,1020)
      WRITE (IPRT,1060)
      CALL CORR(YM, N, M, IYM, LDSTAK)
      CALL MSGX(0, IPRT)
C
C     TEST CORRS
C
C     PRINTOUT SUPPRESSED
C
      WRITE (IPRT,1030)
      WRITE (IPRT,1040)
      WRITE (IPRT,1060)
      CALL CORRS(YM, N, M, IYM, LDSTAK, 0, VCV, IVCV)
      CALL MSGX(0, IPRT)
C
C     PRINT STORED OUTPUT AND ZERO ARRAYS
C
      CALL CORRXP(M, VCV, IVCV, IPRT)
C
C     WITH PRINTOUT
C
      WRITE (IPRT,1050)
      WRITE (IPRT,1060)
      CALL CORRS(YM, N, M, IYM, LDSTAK, 1, VCV, IVCV)
      CALL MSGX(0, IPRT)
C
C     PRINT STORED OUTPUT
C
      CALL CORRXP(M, VCV, IVCV, IPRT)
C
C**** SPECIAL 2 COLUMN MATRIX ****
C
      WRITE (IPRT,1070)
      WRITE (IPRT,1060)
      CALL CORR(YM, N, 2, IYM, LDSTAK)
      CALL MSGX(0, IPRT)
C
C**** TEST WORK AREA REQUIREMENTS ****
C
C     TEST CORR
C
      CALL LDSCMP(12, 0, MAX(N,M), 0, 0, 0, 'S',
     +   M*M + (MAX(N,M)+M+N*(M+3)+6*M*M), LDSMIN)
      WRITE (IPRT,1090)
      CALL CORR(YM, N, M, IYM, LDSMIN-1)
      CALL MSGX(1, IPRT)
      WRITE (IPRT,1100)
      CALL CORR(YM, N, M, IYM, LDSMIN)
      CALL MSGX(0, IPRT)
C
C     TEST CORRS WITH PRINTOUT
C
      CALL LDSCMP(12, 0, MAX(N,M), 0, 0, 0, 'S',
     +   MAX(N,M)+M+N*(M+3)+6*M*M, LDSMIN)
      WRITE (IPRT,1090)
      CALL CORRS(YM, N, M, IYM, LDSMIN-1, 1, VCV, IVCV)
      CALL MSGX(1, IPRT)
      WRITE (IPRT,1100)
      CALL CORRS(YM, N, M, IYM, LDSMIN, 1, VCV, IVCV)
      CALL CORRXP(M, VCV, IVCV, IPRT)
      CALL MSGX(0, IPRT)
C
C     TEST CORRS WITHOUT PRINTOUT
C
      CALL LDSCMP(12, 0, 0, 0, 0, 0, 'S', 0, LDSMIN)
      WRITE (IPRT,1090)
      CALL CORRS(YM, N, M, IYM, LDSMIN-1, 0, VCV, IVCV)
      CALL MSGX(1, IPRT)
      WRITE (IPRT,1100)
      CALL CORRS(YM, N, M, IYM, LDSMIN, 0, VCV, IVCV)
      CALL CORRXP(M, VCV, IVCV, IPRT)
      CALL MSGX(0, IPRT)
C
C**** NUMBER OF VARIABLES LESS THAN 2 ****
C
      WRITE (IPRT,1110)
C
C     TEST CORR
C
      CALL CORR(YM, N, 1, IYM, LDSTAK)
      CALL MSGX(1, IPRT)
C
C     TEST CORRS
C
      CALL CORRS(YM, N, 1, IYM, LDSTAK, 1, VCV, IVCV)
      CALL MSGX(1, IPRT)
C
C**** NUMBER OF OBSERVATIONS LESS THAN 3 ****
C
      WRITE (IPRT,1120)
C
C     TEST CORR
C
      CALL CORR(YM, 2, 4, IYM, LDSTAK)
      CALL MSGX(1, IPRT)
C
C     TEST CORRS
C
      CALL CORRS(YM, 2, 4, IYM, LDSTAK, 1, VCV, IVCV)
      CALL MSGX(1, IPRT)
C
C**** OBSERVATION MATRIX DIMENSIONED LESS THAN N ****
C
      WRITE (IPRT,1150)
C
C     TEST CORR
C
      CALL CORR(YM, N, M, 8, LDSTAK)
      CALL MSGX(1, IPRT)
C
C     TEST CORRS
C
      CALL CORRS(YM, N, M, 8, LDSTAK, 1, VCV, IVCV)
      CALL MSGX(1, IPRT)
C
C**** VCV MATRIX DIMENSIONED LESS THAN M ****
C
      WRITE (IPRT,1130)
      CALL CORRS(YM, N, M, IYM, LDSTAK, 1, VCV, 2)
      CALL MSGX(1, IPRT)
C
C**** ALL OBSERVATIONS ON A SINGLE VARIABLE EQUAL TO ZERO ****
C
      WRITE (IPRT,1140)
      CALL SETRA(Z, 10, 4, 10, 0.0E0)
      CALL CORR(Z, 9, 4, 10, LDSTAK)
      CALL MSGX(1, IPRT)
      CALL CORRS(Z, 9, 4, 10, LDSTAK, 1, VCV, IVCV)
      CALL CORRXP(M, VCV, IVCV, IPRT)
      CALL MSGX(1, IPRT)
C
      DO 10 I=1,10
         Z(I,1) = I
         Z(I,2) = 0.0E0
   10 CONTINUE
      CALL CORR(Z, 10, 4, 10, LDSTAK)
      CALL MSGX(1, IPRT)
C
C**** ARRAY FILLED WITH A SINGLE VALUE ****
C
      WRITE (IPRT,1160)
      CALL SETRA(Z, 10, 4, 10, 4.0E0)
      CALL CORR(Z, 4, 10, 4, LDSTAK)
      CALL MSGX(1, IPRT)
C
C**** 2 COLUMNS THE SAME ****
C
      DO 20 I=1,3
         CALL GENR(Z(1,I), 5, 5.0E0*I, 5.0E0*I)
   20 CONTINUE
      DO 30 I=1,5
         Z(I,4) = Z(I,3)
   30 CONTINUE
      WRITE (IPRT,1170)
      CALL CORR(Z, 5, 4, 10, LDSTAK)
      CALL MSGX(1, IPRT)
C
C**** 2 COLUMNS INVERSELY RELATED ****
C
      J = 5
      DO 40 I=1,5
         J = J - 1
         Z(J,4) = Z(I,3)
   40 CONTINUE
      WRITE (IPRT,1170)
      CALL CORR(Z, 5, 4, 10, LDSTAK)
      CALL MSGX(1, IPRT)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT('1')
 1010 FORMAT(' ****TEST ROUTINES WITH CORRECT CALL****')
 1020 FORMAT(' TEST OF CORR')
 1030 FORMAT('1TEST OF CORRS')
 1040 FORMAT(' PRINTOUT SUPRESSED.')
 1050 FORMAT('1PRINTOUT NOT SUPRESSED.')
 1060 FORMAT(' DRAPER AND SMITH DATA SET (PAGE 216).')
 1070 FORMAT('1****SPECIAL CASE 2 COLUMN MATRIX****')
 1090 FORMAT('1****TEST WITH INSUFFICIENT WORK AREA****')
 1100 FORMAT('1****TEST WITH EXACTLY THE RIGHT AMOUNT OF WORK AREA****')
 1110 FORMAT('1****NUMBER OF VARIABLES LESS THAN 2****')
 1120 FORMAT(' ****NUMBER OF OBSERVATIONS LESS THAN 3****')
 1130 FORMAT(' ****INADEQUATE SPACE IN STORAGE ARRAYS****')
 1140 FORMAT('1****ALL OBSERVATIONS ON A VARIABLE EQUAL TO ZERO****')
 1150 FORMAT(' ****OBSERVATION MATRIX DIMENSIONED LESS THAN NUMBER',
     +       ' OF OBSERVATIONS DESIGNATED****')
 1160 FORMAT('1****ARRAY CONTAINING A SINGLE VALUE****')
 1170 FORMAT('1****2 COLUMNS RELATED****')
      END

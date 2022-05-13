*XXCH8
      SUBROUTINE XXCH8(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     TEST SUBPROGRAM FOR SIMPLE TEST OF
C     THE LINEAR LEAST SQUARES FAMILY OF ROUTINES.
C
C     LLS PROBLEM IS FROM DANIAL AND WOOD [1971], PAGES 61-65.
C
C     LLSP PROBLEM IS FROM MILLER AND FREUND [1977], PAGE 311.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  AUGUST 3, 1987
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
     +   IPRT,IXM,N1,N2,NDEG,NPAR
C
C  LOCAL ARRAYS
      REAL
     +   RES(25),X(25),XM(25,5),Y1(25),Y2(25)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,LLS,LLSP
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IERR
C        THE INTEGER VALUE DESIGNATING WHETHER ANY ERRORS WERE
C        DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        LOGICAL OUTPUT UNIT.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE MATRIX X.
C     INTEGER LDSTAK
C        THE LENGTH OF DSTAK IN COMMON /CSTAK/.
C     INTEGER NDEG
C        THE DEGREE OF THE POLYNOMIAL MODEL TO BE FIT.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS TO BE ESTIMATED.
C     INTEGER N1, N2
C        THE NUMBER OF OBSERVATIONS IN EACH PROBLEM.
C     REAL RES(25)
C        THE RESIDUALS.
C     REAL X(25)
C        THE INDEPENDENT VARIABLE.
C     REAL XM(25,5)
C        THE INDEPENDENT VARIABLE.
C     REAL Y1(25), Y2(25)
C        THE DEPENDENT VARIABLE.
C
C
      DATA      XM(1,1),  XM(1,2),  XM(1,3),  XM(1,4)
     +    /      1.0E0, 80.0E0, 27.0E0, 89.0E0/
      DATA      XM(2,1),  XM(2,2),  XM(2,3),  XM(2,4)
     +    /      1.0E0, 80.0E0, 27.0E0, 88.0E0/
      DATA      XM(3,1),  XM(3,2),  XM(3,3),  XM(3,4)
     +    /      1.0E0, 75.0E0, 25.0E0, 90.0E0/
      DATA      XM(4,1),  XM(4,2),  XM(4,3),  XM(4,4)
     +    /      1.0E0, 62.0E0, 24.0E0, 87.0E0/
      DATA      XM(5,1),  XM(5,2),  XM(5,3),  XM(5,4)
     +    /      1.0E0, 62.0E0, 22.0E0, 87.0E0/
      DATA      XM(6,1),  XM(6,2),  XM(6,3),  XM(6,4)
     +    /      1.0E0, 62.0E0, 23.0E0, 87.0E0/
      DATA      XM(7,1),  XM(7,2),  XM(7,3),  XM(7,4)
     +    /      1.0E0, 62.0E0, 24.0E0, 93.0E0/
      DATA      XM(8,1),  XM(8,2),  XM(8,3),  XM(8,4)
     +    /      1.0E0, 62.0E0, 24.0E0, 93.0E0/
      DATA      XM(9,1),  XM(9,2),  XM(9,3),  XM(9,4)
     +    /      1.0E0, 58.0E0, 23.0E0, 87.0E0/
      DATA     XM(10,1), XM(10,2), XM(10,3), XM(10,4)
     +    /      1.0E0, 58.0E0, 18.0E0, 80.0E0/
      DATA     XM(11,1), XM(11,2), XM(11,3), XM(11,4)
     +    /      1.0E0, 58.0E0, 18.0E0, 89.0E0/
      DATA     XM(12,1), XM(12,2), XM(12,3), XM(12,4)
     +    /      1.0E0, 58.0E0, 17.0E0, 88.0E0/
      DATA     XM(13,1), XM(13,2), XM(13,3), XM(13,4)
     +    /      1.0E0, 58.0E0, 18.0E0, 82.0E0/
      DATA     XM(14,1), XM(14,2), XM(14,3), XM(14,4)
     +    /      1.0E0, 58.0E0, 19.0E0, 93.0E0/
      DATA     XM(15,1), XM(15,2), XM(15,3), XM(15,4)
     +    /      1.0E0, 50.0E0, 18.0E0, 89.0E0/
      DATA     XM(16,1), XM(16,2), XM(16,3), XM(16,4)
     +    /      1.0E0, 50.0E0, 18.0E0, 86.0E0/
      DATA     XM(17,1), XM(17,2), XM(17,3), XM(17,4)
     +    /      1.0E0, 50.0E0, 19.0E0, 72.0E0/
      DATA     XM(18,1), XM(18,2), XM(18,3), XM(18,4)
     +    /      1.0E0, 50.0E0, 19.0E0, 79.0E0/
      DATA     XM(19,1), XM(19,2), XM(19,3), XM(19,4)
     +    /      1.0E0, 50.0E0, 20.0E0, 80.0E0/
      DATA     XM(20,1), XM(20,2), XM(20,3), XM(20,4)
     +    /      1.0E0, 56.0E0, 20.0E0, 82.0E0/
      DATA     XM(21,1), XM(21,2), XM(21,3), XM(21,4)
     +    /      1.0E0, 70.0E0, 20.0E0, 91.0E0/
C
      DATA        Y1(1),    Y1(2),    Y1(3)
     +    /     42.0E0, 37.0E0, 37.0E0/
      DATA        Y1(4),    Y1(5),    Y1(6)
     +    /     28.0E0, 18.0E0, 18.0E0/
      DATA        Y1(7),    Y1(8),    Y1(9)
     +    /     19.0E0, 20.0E0, 15.0E0/
      DATA       Y1(10),   Y1(11),   Y1(12)
     +    /     14.0E0, 14.0E0, 13.0E0/
      DATA       Y1(13),   Y1(14),   Y1(15)
     +    /     11.0E0, 12.0E0,  8.0E0/
      DATA       Y1(16),   Y1(17),   Y1(18)
     +    /      7.0E0,  8.0E0,  8.0E0/
      DATA       Y1(19),   Y1(20),   Y1(21)
     +    /      9.0E0, 15.0E0, 15.0E0/
C
      DATA         X(1),     X(2),     X(3)
     +    /      0.0E0,  1.0E0,  2.0E0/
      DATA         X(4),     X(5),     X(6)
     +    /      3.0E0,  4.0E0,  5.0E0/
      DATA         X(7),     X(8),     X(9)
     +    /      6.0E0,  7.0E0,  8.0E0/
C
      DATA        Y2(1),    Y2(2),    Y2(3)
     +    /     12.0E0, 10.5E0, 10.0E0/
      DATA        Y2(4),    Y2(5),    Y2(6)
     +    /      8.0E0,  7.0E0,  8.0E0/
      DATA        Y2(7),    Y2(8),    Y2(9)
     +    /      7.5E0,  8.5E0,  9.0E0/
C
C     SET PARAMETERS NECESSARY FOR THE COMPUTATIONS
C
      CALL IPRINT(IPRT)
      IXM = 25
      N1 = 21
      N2 = 9
      NPAR = 4
      NDEG = 2
C
C     PRINT HEADER
C
      WRITE (IPRT,1000)
C
C     RUN SIMPLE EXAMPLE OF LLS
C
      WRITE (IPRT,1100)
      CALL LLS(Y1, XM, N1, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,2000) IERR
C
C     RUN SIMPLE EXAMPLE OF LLSP
C
      WRITE (IPRT,1200)
      CALL LLSP(Y2, X, N2, NDEG, RES, LDSTAK)
      WRITE (IPRT,2000) IERR
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1*CH8')
 1100 FORMAT (' SIMPLE TEST OF LLS')
 1200 FORMAT ('1SIMPLE TEST OF LLSP')
 2000 FORMAT (/' THE VALUE OF IERR IS ', I4)
C
      END

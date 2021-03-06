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
      DOUBLE PRECISION
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
C     DOUBLE PRECISION RES(25)
C        THE RESIDUALS.
C     DOUBLE PRECISION X(25)
C        THE INDEPENDENT VARIABLE.
C     DOUBLE PRECISION XM(25,5)
C        THE INDEPENDENT VARIABLE.
C     DOUBLE PRECISION Y1(25), Y2(25)
C        THE DEPENDENT VARIABLE.
C
C
      DATA      XM(1,1),  XM(1,2),  XM(1,3),  XM(1,4)
     +    /      1.0D0, 80.0D0, 27.0D0, 89.0D0/
      DATA      XM(2,1),  XM(2,2),  XM(2,3),  XM(2,4)
     +    /      1.0D0, 80.0D0, 27.0D0, 88.0D0/
      DATA      XM(3,1),  XM(3,2),  XM(3,3),  XM(3,4)
     +    /      1.0D0, 75.0D0, 25.0D0, 90.0D0/
      DATA      XM(4,1),  XM(4,2),  XM(4,3),  XM(4,4)
     +    /      1.0D0, 62.0D0, 24.0D0, 87.0D0/
      DATA      XM(5,1),  XM(5,2),  XM(5,3),  XM(5,4)
     +    /      1.0D0, 62.0D0, 22.0D0, 87.0D0/
      DATA      XM(6,1),  XM(6,2),  XM(6,3),  XM(6,4)
     +    /      1.0D0, 62.0D0, 23.0D0, 87.0D0/
      DATA      XM(7,1),  XM(7,2),  XM(7,3),  XM(7,4)
     +    /      1.0D0, 62.0D0, 24.0D0, 93.0D0/
      DATA      XM(8,1),  XM(8,2),  XM(8,3),  XM(8,4)
     +    /      1.0D0, 62.0D0, 24.0D0, 93.0D0/
      DATA      XM(9,1),  XM(9,2),  XM(9,3),  XM(9,4)
     +    /      1.0D0, 58.0D0, 23.0D0, 87.0D0/
      DATA     XM(10,1), XM(10,2), XM(10,3), XM(10,4)
     +    /      1.0D0, 58.0D0, 18.0D0, 80.0D0/
      DATA     XM(11,1), XM(11,2), XM(11,3), XM(11,4)
     +    /      1.0D0, 58.0D0, 18.0D0, 89.0D0/
      DATA     XM(12,1), XM(12,2), XM(12,3), XM(12,4)
     +    /      1.0D0, 58.0D0, 17.0D0, 88.0D0/
      DATA     XM(13,1), XM(13,2), XM(13,3), XM(13,4)
     +    /      1.0D0, 58.0D0, 18.0D0, 82.0D0/
      DATA     XM(14,1), XM(14,2), XM(14,3), XM(14,4)
     +    /      1.0D0, 58.0D0, 19.0D0, 93.0D0/
      DATA     XM(15,1), XM(15,2), XM(15,3), XM(15,4)
     +    /      1.0D0, 50.0D0, 18.0D0, 89.0D0/
      DATA     XM(16,1), XM(16,2), XM(16,3), XM(16,4)
     +    /      1.0D0, 50.0D0, 18.0D0, 86.0D0/
      DATA     XM(17,1), XM(17,2), XM(17,3), XM(17,4)
     +    /      1.0D0, 50.0D0, 19.0D0, 72.0D0/
      DATA     XM(18,1), XM(18,2), XM(18,3), XM(18,4)
     +    /      1.0D0, 50.0D0, 19.0D0, 79.0D0/
      DATA     XM(19,1), XM(19,2), XM(19,3), XM(19,4)
     +    /      1.0D0, 50.0D0, 20.0D0, 80.0D0/
      DATA     XM(20,1), XM(20,2), XM(20,3), XM(20,4)
     +    /      1.0D0, 56.0D0, 20.0D0, 82.0D0/
      DATA     XM(21,1), XM(21,2), XM(21,3), XM(21,4)
     +    /      1.0D0, 70.0D0, 20.0D0, 91.0D0/
C
      DATA        Y1(1),    Y1(2),    Y1(3)
     +    /     42.0D0, 37.0D0, 37.0D0/
      DATA        Y1(4),    Y1(5),    Y1(6)
     +    /     28.0D0, 18.0D0, 18.0D0/
      DATA        Y1(7),    Y1(8),    Y1(9)
     +    /     19.0D0, 20.0D0, 15.0D0/
      DATA       Y1(10),   Y1(11),   Y1(12)
     +    /     14.0D0, 14.0D0, 13.0D0/
      DATA       Y1(13),   Y1(14),   Y1(15)
     +    /     11.0D0, 12.0D0,  8.0D0/
      DATA       Y1(16),   Y1(17),   Y1(18)
     +    /      7.0D0,  8.0D0,  8.0D0/
      DATA       Y1(19),   Y1(20),   Y1(21)
     +    /      9.0D0, 15.0D0, 15.0D0/
C
      DATA         X(1),     X(2),     X(3)
     +    /      0.0D0,  1.0D0,  2.0D0/
      DATA         X(4),     X(5),     X(6)
     +    /      3.0D0,  4.0D0,  5.0D0/
      DATA         X(7),     X(8),     X(9)
     +    /      6.0D0,  7.0D0,  8.0D0/
C
      DATA        Y2(1),    Y2(2),    Y2(3)
     +    /     12.0D0, 10.5D0, 10.0D0/
      DATA        Y2(4),    Y2(5),    Y2(6)
     +    /      8.0D0,  7.0D0,  8.0D0/
      DATA        Y2(7),    Y2(8),    Y2(9)
     +    /      7.5D0,  8.5D0,  9.0D0/
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

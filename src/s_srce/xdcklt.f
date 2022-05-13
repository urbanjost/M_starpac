*XDCKLT
      SUBROUTINE XDCKLT(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST DERIVATIVE CHECKING ROUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 21, 1983
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
      REAL
     +   DELTA
      INTEGER
     +   I,IPRT,IXM,J,JSTOP,LDSMIN,M,N,NETA,NPAR,NPRT,NROW,NTAU,
     +   NTEST
C
C  LOCAL ARRAYS
      REAL
     +   PAR(10),SCALE(10),XM(200,2)
      INTEGER
     +   NETTST(6),NROTST(5),NTATST(6)
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DCKLS,DCKLS1,DCKLSC,DRV4A,DRV4B,IPRINT,LDSCMP,MDL4,SETRV
C
C  INTRINSIC FUNCTIONS
      INTRINSIC LOG10
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL DELTA
C        *
C     EXTERNAL DRV4A, DRV4B
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        ANALYTIC DERIVATIVES (JACOBIAN MATRIX) OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY XM.
C     INTEGER J, JSTOP
C        *
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH OF THE ARRAY DSTAK ALLOWED.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL4
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS OF DATA.
C     INTEGER NETA
C        THE NUMBER OF RELIABLE DIGITS IN THE MODEL.
C     INTEGER NETTST(6)
C        VARIOUS TEST VALUES FOR NETA.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE PROVIDED, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO PRINTED OUTPUT IS GIVEN.
C     INTEGER NROTST(5)
C        VARIOUS TEST VALUES FOR NROW.
C     INTEGER NROW
C        THE NUMBER OF THE ROW OF THE INDEPENDENT VARIABLE ARRAY AT
C        WHICH THE DERIVATIVE IS TO BE CHECKED.
C     INTEGER NTATST(6)
C         VARIOUS TEST VALUES FOR NTAU.
C     INTEGER NTAU
C        THE NUMBER OF DIGITS OF AGREEMENT REQUIRED BETWEEN THE
C        NUMERICAL DERIVATIVES AND THE USER SUPPLIED DERIVATIVES.
C     INTEGER NTEST
C        THE NUMBER OF THE CURRENT TEST.
C     REAL PAR(10)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     REAL SCALE(10)
C        A DUMMY ARRAY, INDICATING USE OF DEFAULT VALUES FOR
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL XM(200,2)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C
      CALL IPRINT(IPRT)
C
C     SET PARAMETER VALUES
C
      CALL DCKLS1(N, M, IXM, PAR, NPAR, NETA, NTAU, NROW, SCALE, NPRT)
      CALL LDSCMP(5, 0, 2*NPAR+1, 0, 0, 0, 'S',
     +            N*NPAR+NPAR+N, LDSMIN)
C
      IF (LDSMIN.LE.LDSTAK) GO TO 5
C
      WRITE (IPRT, 1020) LDSMIN
      RETURN
C
    5 CONTINUE
C
C     CREATE INDEPENDENT VARIABLE
C
      DELTA = 0.0625E0
      XM(1,1) = 0.0E0
      DO 10 I=2,N
         XM(I,1) = XM(I-1,1) + DELTA
   10 CONTINUE
C
      NTEST = 0
C
C
C     TEST VARIOUS VALUES OF NETA AND NTAU
C
      SCALE(1) = 0.0E0
C
      NETTST(1) = -1
      NETTST(2) = 0
      NETTST(3) = 1
      NETTST(4) = 2
C
      NETTST(5) = -LOG10(R1MACH(4))
      NETTST(6) = NETTST(5) + 1
C
      NTATST(1) = -1
      NTATST(2) = 0
      NTATST(3) = 1
C
      JSTOP = 3
C
      DO 30 I=1,6
C
         NTATST(4) = NETTST(I)/4
         IF (I.LE.5) THEN
            NTATST(5) = (NETTST(I)-1)/2
            NTATST(6) = NTATST(5) + 1
         END IF
C
         IF (I.EQ.5) JSTOP = 6
C
         DO 20 J=1,JSTOP
C
            NTEST = NTEST + 1
            WRITE (IPRT,1130) NTEST
            WRITE (IPRT,1100)
            WRITE (IPRT,1040)
            WRITE (IPRT,1060) NETTST(I), NTATST(J), SCALE(1), NROW, NPRT
            WRITE (IPRT,1000)
            IERR = -1
            CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4A, PAR, NPAR, LDSTAK,
     +         NETTST(I), NTATST(J), SCALE, NROW, NPRT)
            WRITE (IPRT,1050) IERR
            WRITE (IPRT,1140) NETTST(I), NTATST(J), SCALE(1), NROW, NPRT
C
   20    CONTINUE
C
   30 CONTINUE
C
C     SUPPRESS OUTPUT
C
      NPRT = 0
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1100)
      WRITE (IPRT,1040)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4A, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     LARGE CALCULATION ERROR PROBLEM
C
      CALL DCKLS1(N, M, IXM, PAR, NPAR, NETA, NTAU, NROW, SCALE, NPRT)
      PAR(3) = 10.0E0**NTATST(5)
      SCALE(1) = 0.0E0
      NROW = 51
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1100)
      WRITE (IPRT,1070)
      WRITE (IPRT,1080)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4A, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     NEARLY ZERO DERIVATIVE
C
      NROW = 50
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1100)
      WRITE (IPRT,1070)
      WRITE (IPRT,1090)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4A, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     INCORRECTLY CODED DERIVATIVE
C
C     SIMPLE EXAMPLE
C
C     SET PARAMETER VALUES
C
      CALL DCKLS1(N, M, IXM, PAR, NPAR, NETA, NTAU, NROW, SCALE, NPRT)
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1110)
      WRITE (IPRT,1040)
      WRITE (IPRT,1000)
      IERR = -1
      CALL DCKLS(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK)
      WRITE (IPRT,1050) IERR
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1110)
      WRITE (IPRT,1040)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     SUPPRESS OUTPUT
C
      NPRT = 0
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1110)
      WRITE (IPRT,1040)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     LARGE CALCULATION ERROR PROBLEM
C
      CALL DCKLS1(N, M, IXM, PAR, NPAR, NETA, NTAU, NROW, SCALE, NPRT)
C
      PAR(3) = 10.0E0**NTATST(5)
      NROW = 26
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1110)
      WRITE (IPRT,1070)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
      PAR(4) = 0.75E0
      NROW = 1
C
      NTEST = NTEST + 1
      WRITE (IPRT,1130) NTEST
      WRITE (IPRT,1110)
      WRITE (IPRT,1070)
      WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROW, NPRT
      WRITE (IPRT,1010)
      IERR = -1
      CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK, NETA,
     +   NTAU, SCALE, NROW, NPRT)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROW, NPRT
C
C     CHECK VARIOUS VALUES OF NROW
C
      CALL DCKLS1(N, M, IXM, PAR, NPAR, NETA, NTAU, NROW, SCALE, NPRT)
C
      CALL SETRV(XM(1,1), N, 0.0E0)
      NROTST(1) = -1
      NROTST(2) = 0
      NROTST(3) = 1
      NROTST(4) = N
      NROTST(5) = N + 1
C
      DO 40 I=1,5
C
         NTEST = NTEST + 1
         WRITE (IPRT,1130) NTEST
         WRITE (IPRT,1110)
         WRITE (IPRT,1120)
         WRITE (IPRT,1060) NETA, NTAU, SCALE(1), NROTST(I), NPRT
         WRITE (IPRT,1010)
         IERR = -1
         CALL DCKLSC(XM, N, M, IXM, MDL4, DRV4B, PAR, NPAR, LDSTAK,
     +      NETA, NTAU, SCALE, NROTST(I), NPRT)
         WRITE (IPRT,1050) IERR
         WRITE (IPRT,1140) NETA, NTAU, SCALE(1), NROTST(I), NPRT
C
   40 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (15H TEST OF DCKLS )
 1010 FORMAT (15H TEST OF DCKLSC)
 1020 FORMAT (45H1 *** LDSTAK MUST BE GREATER THAN OR EQUAL TO , I6)
 1040 FORMAT (15H SIMPLE EXAMPLE)
 1050 FORMAT (29H ***** RETURNED RESULTS *****, 5X, 15H (-1 INDICATES ,
     +   39HVALUE NOT CHANGED BY CALLED SUBROUTINE)//9H IERR IS , I3)
 1060 FORMAT (19H INPUT   -  NETA = , I5, 9H, NTAU = , I5,
     +   13H, SCALE(1) = , G15.8, 9H, NROW = , I5, 9H, NPRT = , I5)
 1070 FORMAT (32H LARGE CALCULATION ERROR PROBLEM)
 1080 FORMAT (16H ZERO DERIVATIVE)
 1090 FORMAT (23H NEARLY ZERO DERIVATIVE)
 1100 FORMAT (27H CORRECTLY CODED DERIVATIVE)
 1110 FORMAT (' INCORRECTLY CODED DERIVATIVE FOR PARAMETERS 1, 2 AND 4')
 1120 FORMAT (' ALL INDEPENDENT VARIABLES EQUAL TO ZERO')
 1130 FORMAT (43H1DERIVATIVE CHECKING SUBROUTINE TEST NUMBER, I5)
 1140 FORMAT (19H OUTPUT  -  NETA = , I5, 9H, NTAU = , I5,
     +   13H, SCALE(1) = , G15.8, 9H, NROW = , I5, 9H, NPRT = , I5//)
      END

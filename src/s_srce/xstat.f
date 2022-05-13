*XSTAT
      SUBROUTINE XSTAT(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS PROGRAM TESTS FEATURES OF THE STAT FAMILY TO ENSURE THAT
C     ALL ASPECTS OF THE STAT FAMILY ROUTINES WORK CORRECTLY.
C
C     WRITTEN BY  -  JOHN E. KOONTZ AND JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
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
      REAL
     +   FPLM,WTEMP,YTEMP1,YTEMPN
      INTEGER
     +   I,IPRT,N,NCONST,NPRTOF,NPRTON
C
C  LOCAL ARRAYS
      REAL
     +   STS(53),WT(84),WTALL0(10),WTALL1(84),Y(84),YCONST(10),
     +   YPATH(10)
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,STAT,STATS,STATW,STATWS
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        A LOOP INDEX.
C     INTEGER IERR
C        FLAG TO INDICATE PRESENCE OF ERROR DETECTED BY PRECEDING
C        STARPAC CALL.  (0 IS OK, 1 IS ERROR)
C     INTEGER IPRT
C        LOGICAL OUTPUT UNIT.
C     INTEGER LDSTAK
C        AMOUNT OF WORK AREA.  SIZE OF DSTAK.
C     INTEGER N
C        THE LENGTH OF THE VECTOR Y.
C     INTEGER NCONST
C        LENGTH OF THE VECTOR YCONST.
C     INTEGER NPRTOF
C        FLAG FOR NO OUTPUT (EXCEPT ERROR MESSAGES).
C     INTEGER NPRTON
C        FLAG FOR FULL PRINTOUT.
C     REAL STS(53)
C        VECTOR OF STATISTICS.
C     REAL WT(84)
C        WEIGHTS VECTOR.
C     REAL WTALL0(10)
C        N VECTOR OF 0 WEIGHTS.
C     REAL WTALL1(84)
C        N VECTOR OF 1 WEIGHTS.
C     REAL WTEMP
C        TEMPORARY STORAGE FOR ONE OF THE WEIGHTS.
C     REAL Y(84)
C        DATA VECTOR FOR TESTS.
C     REAL YCONST(10)
C        VECTOR OF CONSTANT DATA.
C     REAL YPATH(10)
C        A VECTOR OF Y VALUES DESIGNED TO FORCE DIFFERENT PATHS
C        THROUGH THE SUMMATION ROUTINES.
C     REAL YTEMPN, YTEMP1
C        TEMPORARY STORAGE FOR THE FIRST AND LAST Y VALUE.
C
C     DATA INITIALIZATIONS.
C
      DATA N /84/
      DATA NCONST /10/
      DATA NPRTON /1/
      DATA NPRTOF /0/
C
C     DAVIS-HARRISON R.H. DATA, PIKES PEAK.
C
C     THIS IS AN ARBITRARILY CHOSEN DATA SET.
C
      DATA Y( 1), Y( 2), Y( 3), Y( 4)
     +    / 0.6067E0, 0.6087E0, 0.6086E0, 0.6134E0/
      DATA Y( 5), Y( 6), Y( 7)
     +    / 0.6108E0, 0.6138E0, 0.6125E0/
      DATA Y( 8), Y( 9), Y(10), Y(11)
     +    / 0.6122E0, 0.6110E0, 0.6104E0, 0.7213E0/
      DATA Y(12), Y(13), Y(14)
     +    / 0.7078E0, 0.7021E0, 0.7004E0/
      DATA Y(15), Y(16), Y(17), Y(18)
     +    / 0.6981E0, 0.7242E0, 0.7268E0, 0.7418E0/
      DATA Y(19), Y(20), Y(21)
     +    / 0.7407E0, 0.7199E0, 0.6225E0/
      DATA Y(22), Y(23), Y(24), Y(25)
     +    / 0.6254E0, 0.6252E0, 0.6267E0, 0.6218E0/
      DATA Y(26), Y(27), Y(28)
     +    / 0.6178E0, 0.6216E0, 0.6192E0/
      DATA Y(29), Y(30), Y(31), Y(32)
     +    / 0.6191E0, 0.6250E0, 0.6188E0, 0.6233E0/
      DATA Y(33), Y(34), Y(35)
     +    / 0.6225E0, 0.6204E0, 0.6207E0/
      DATA Y(36), Y(37), Y(38), Y(39)
     +    / 0.6168E0, 0.6141E0, 0.6291E0, 0.6231E0/
      DATA Y(40), Y(41), Y(42)
     +    / 0.6222E0, 0.6252E0, 0.6308E0/
      DATA Y(43), Y(44), Y(45), Y(46)
     +    / 0.6376E0, 0.6330E0, 0.6303E0, 0.6301E0/
      DATA Y(47), Y(48), Y(49)
     +    / 0.6390E0, 0.6423E0, 0.6300E0/
      DATA Y(50), Y(51), Y(52), Y(53)
     +    / 0.6260E0, 0.6292E0, 0.6298E0, 0.6290E0/
      DATA Y(54), Y(55), Y(56)
     +    / 0.6262E0, 0.5952E0, 0.5951E0/
      DATA Y(57), Y(58), Y(59), Y(60)
     +    / 0.6314E0, 0.6440E0, 0.6439E0, 0.6326E0/
      DATA Y(61), Y(62), Y(63)
     +    / 0.6392E0, 0.6417E0, 0.6412E0/
      DATA Y(64), Y(65), Y(66), Y(67)
     +    / 0.6530E0, 0.6411E0, 0.6355E0, 0.6344E0/
      DATA Y(68), Y(69), Y(70)
     +    / 0.6623E0, 0.6276E0, 0.6307E0/
      DATA Y(71), Y(72), Y(73), Y(74)
     +    / 0.6354E0, 0.6197E0, 0.6153E0, 0.6340E0/
      DATA Y(75), Y(76), Y(77)
     +    / 0.6338E0, 0.6284E0, 0.6162E0/
      DATA Y(78), Y(79), Y(80), Y(81)
     +    / 0.6252E0, 0.6349E0, 0.6344E0, 0.6361E0/
      DATA Y(82), Y(83), Y(84)
     +    / 0.6373E0, 0.6337E0, 0.6383E0/
      DATA WT( 1), WT( 2), WT( 3), WT( 4), WT( 5), WT( 6), WT( 7)
     +   / 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0/
      DATA WT( 8), WT( 9), WT(10), WT(11), WT(12), WT(13), WT(14)
     +   / 0.5E0, 0.5E0, 0.5E0, 0.0E0, 0.0E0, 0.0E0, 0.0E0/
      DATA WT(15), WT(16), WT(17), WT(18), WT(19), WT(20), WT(21)
     +   / 0.0E0, 0.0E0, 0.0E0, 0.0E0, 0.0E0, 0.0E0, 0.5E0/
      DATA WT(22), WT(23), WT(24), WT(25), WT(26), WT(27), WT(28)
     +   / 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0/
      DATA WT(29), WT(30), WT(31), WT(32), WT(33), WT(34), WT(35)
     +   / 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0/
      DATA WT(36), WT(37), WT(38), WT(39), WT(40), WT(41), WT(42)
     +   / 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0, 0.5E0/
      DATA WT(43), WT(44), WT(45), WT(46), WT(47), WT(48), WT(49)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0/
      DATA WT(50), WT(51), WT(52), WT(53), WT(54), WT(55), WT(56)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 0.0E0, 0.0E0/
      DATA WT(57), WT(58), WT(59), WT(60), WT(61), WT(62), WT(63)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0/
      DATA WT(64), WT(65), WT(66), WT(67), WT(68), WT(69), WT(70)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0/
      DATA WT(71), WT(72), WT(73), WT(74), WT(75), WT(76), WT(77)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0/
      DATA WT(78), WT(79), WT(80), WT(81), WT(82), WT(83), WT(84)
     +   / 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0, 1.0E0/
C
C     DEFINE IPRT, THE CURRENT OUTPUT UNIT.
C
      CALL IPRINT(IPRT)
C
      FPLM = R1MACH(2)
C
C     SET UP THE WEIGHTS VECTORS.
C
      DO 10 I=1,N
         WTALL1(I) = 1.0E0
   10 CONTINUE
      DO 20 I=1,NCONST
         YCONST(I) = 1.0E0
         WTALL0(I) = 0.0E0
   20 CONTINUE
C
C     HEADING.
C
      WRITE (IPRT,1150)
C
C     TEST 1.  CHECK ALL ERROR MESSAGES.
C
C     ERROR 1, TWO OR FEWER ELEMENTS.
C
      WRITE (IPRT,1180)
      WRITE(IPRT,1230)
      WRITE(IPRT,1240)
      CALL STAT(Y, 2, LDSTAK)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1250)
      CALL STATS(Y, 2, LDSTAK, STS, NPRTON)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1400)
      CALL STATW(Y, WT, 2, LDSTAK)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, 2, LDSTAK, STS, NPRTON)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
C
C     ERROR 2, NOT ENOUGH SPACE IN CSTAK.
C
      WRITE (IPRT,1190)
      WRITE(IPRT,1230)
      WRITE(IPRT,1240)
      CALL STAT(Y, N, N/4)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1250)
      CALL STATS(Y, N, N/4, STS, NPRTON)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1400)
      CALL STATW(Y, WT, N, N/4)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, N, N/4, STS, NPRTON)
      WRITE (IPRT,1170) IERR
C
C     ERROR 4, NEGATIVE WEIGHTS.
C
      WRITE (IPRT,1210)
      WTEMP = WT(2)
      WT(2) = -1.0E0
      WRITE(IPRT,1230)
      WRITE(IPRT,1400)
      CALL STATW(Y, WT, N, LDSTAK)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, N, LDSTAK, STS, NPRTON)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
      WT(2) = WTEMP
C
C     ERROR 5, ALL WEIGHTS ZERO (PLUS CONSTANT Y).
C
      WRITE (IPRT,1220)
      WRITE(IPRT,1230)
      WRITE(IPRT,1400)
      CALL STATW(YCONST, WTALL0, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1230)
      WRITE(IPRT,1410)
      CALL STATWS(YCONST, WTALL0, NCONST, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1170) IERR
C
C     TEST 2.  CHECK FOR READING OUTSIDE OF DATA ARRAY.
C
      WRITE (IPRT,1160)
      YTEMP1 = YCONST(1)
      YCONST(1) = FPLM
      YTEMPN = YCONST(NCONST)
      YCONST(NCONST) = FPLM
      WRITE(IPRT,1440)
      WRITE(IPRT,1240)
      CALL STAT(YCONST(2), NCONST-2, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1250)
      CALL STATS(YCONST(2), NCONST-2, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1400)
      CALL STATW(YCONST(2), WT, NCONST-2, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1410)
      CALL STATWS(YCONST(2), WT, NCONST-2, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1170) IERR
      YCONST(1) = YTEMP1
      YCONST(NCONST) = YTEMPN
C
C     TEST 3.  CONSTANT Y.
C
      WRITE (IPRT,1200)
      WRITE(IPRT,1440)
      WRITE(IPRT,1240)
      CALL STAT(YCONST, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1250)
      CALL STATS(YCONST, NCONST, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1400)
      CALL STATW(YCONST, WT, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE(IPRT,1440)
      WRITE(IPRT,1410)
      CALL STATWS(YCONST, WT, NCONST, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1170) IERR
C
C     TEST 4.  SEE IF TURNING OFF THE PRINTOUT WORKS.
C
      WRITE (IPRT,1260)
      WRITE (IPRT,1270)
      WRITE(IPRT,1230)
      WRITE(IPRT,1250)
      CALL STATS(Y, N, LDSTAK, STS, NPRTOF)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1280)
      WRITE(IPRT,1230)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, N, LDSTAK, STS, NPRTOF)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
C
C     TEST 5.  MAKE A WORKING RUN OF EACH ROUTINE  FIRST WITH
C              N=2 (THE MINIMUN VALID VALUE) AND THEN FOR THE WHOLE
C              DATA SET TO CHECK THE OUTPUT.
C
      WRITE (IPRT,1300)
      WRITE (IPRT,1310)
      WRITE(IPRT,1240)
      CALL STAT(Y, 3, LDSTAK)
      WRITE (IPRT,1310)
      WRITE(IPRT,1240)
      CALL STAT(Y, N, LDSTAK)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
C
      WRITE (IPRT,1320)
      WRITE(IPRT,1400)
      CALL STATW(Y, WT, 3, LDSTAK)
      WRITE (IPRT,1320)
      WRITE(IPRT,1400)
      CALL STATW(Y, WT, N, LDSTAK)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
C
      WRITE (IPRT,1340)
      WRITE(IPRT,1250)
      CALL STATS(Y, 3, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1340)
      WRITE(IPRT,1250)
      CALL STATS(Y, N, LDSTAK, STS, NPRTON)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
C
      WRITE (IPRT,1350)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, 3, LDSTAK, STS, NPRTON)
      WRITE (IPRT,1350)
      WRITE(IPRT,1410)
      CALL STATWS(Y, WT, N, LDSTAK, STS, NPRTON)
      WRITE(IPRT,1390) (Y(I), I = 1, 10)
      WRITE (IPRT,1170) IERR
C
C     TEST 5.  CHECK RESULTS OF WEIGHTING ALL OBSERVATIONS
C              WITH 1.0E0.  COMPARE WITH STAT EXECUTION.
C
      WRITE (IPRT,1370)
      WRITE(IPRT,1400)
      CALL STATW(Y, WTALL1, N, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     TEST 6.  CHECK RESULTS OF FORCING DIFFERENCE PATHS THROUGH
C              THE SUMMATION ROUTINES, USING SMALL, SIMPLE DATA SETS.
C
      WRITE (IPRT,1000)
C
C     RUN DATA SET 6.1
C
      DO 30 I=1,10
         YPATH(I) = I
   30 CONTINUE
      WRITE (IPRT,1010)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1020)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.2
C
      DO 40 I=1,10
         YPATH(I) = -I
   40 CONTINUE
      WRITE (IPRT,1030)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1040)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.3
C
      DO 50 I=1,10
         YPATH(I) = I-1
   50 CONTINUE
      WRITE (IPRT,1050)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1060)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.4
C
      DO 60 I=1,10
         YPATH(I) = 1-I
   60 CONTINUE
      WRITE (IPRT,1070)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1080)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.5
C
      DO 70 I=1,10
         YPATH(I) = I-6
   70 CONTINUE
      WRITE (IPRT,1090)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1100)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.6
C
      DO 80 I=1,10
         YPATH(I) = I-5
   80 CONTINUE
      WRITE (IPRT,1110)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1120)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.7
C
      DO 90 I=1,10
         YPATH(I) = 0.0E0
   90 CONTINUE
      YPATH(1) = -5.0E0
      YPATH(10) = 5.0E0
      WRITE (IPRT,1130)
      WRITE(IPRT,1240)
      CALL STAT(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1140)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.8
C
      DO 100 I=1,10
         YPATH(I) = 0.0E0
  100 CONTINUE
      YPATH(1) = -5.0E0
      WTALL1(1) = 0.0E0
      YPATH(10) = 5.0E0
      WTALL1(10) = 0.0E0
      WRITE (IPRT,1380)
      WRITE(IPRT,1400)
      CALL STATW(YPATH, WTALL1, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      RETURN
C
C     FORMATS
C
 1000 FORMAT(51H1TEST 6.  TRY DIFFERENT PATHS THROUGH THE SUMMATION,
     +   6H CODE.)
 1010 FORMAT('1RUN STAT ON 1, ..., 10.')
 1020 FORMAT('1RUN STATW ON 1, ..., 10.  WEIGHTS ARE ALL 1.')
 1030 FORMAT('1RUN STAT ON -1, ..., -10.')
 1040 FORMAT('1RUN STATW ON -1, ..., -10.  WEIGHTS ARE ALL 1.')
 1050 FORMAT('1RUN STAT ON 0, ..., 9.')
 1060 FORMAT('1RUN STATW ON 0, ..., 9.  WEIGHTS ARE ALL 1.')
 1070 FORMAT('1RUN STAT ON 0, ..., -9.')
 1080 FORMAT('1RUN STATW ON 0, ..., -9.  WEIGHTS ARE ALL 1.')
 1090 FORMAT('1STAT ON -5, ..., 4.')
 1100 FORMAT('1RUN STATW ON -5, ..., 4.  WEIGHTS ARE ALL 1.')
 1110 FORMAT('1RUN STAT ON -4, ..., 5.')
 1120 FORMAT('1RUN STATW ON -4, ..., 5.  WEIGHTS ARE ALL 1.')
 1130 FORMAT('1RUN STAT ON -1, 8*0, 1.')
 1140 FORMAT('1RUN STATW ON -1, 8*0, 1.  WEIGHTS ARE ALL 1.')
 1150 FORMAT('1TEST RUNS FOR THE STATISTICAL ANALYSIS FAMILY ROUTINES.')
 1160 FORMAT('1TEST RUNS TO BE SURE CODE IS NOT READING OUTSIDE',
     +       ' DATA ARRAY.')
 1170 FORMAT(/' THE VALUE OF IERR IS ', I4)
 1180 FORMAT('1TRY TWO OR FEWER ELEMENTS.')
 1190 FORMAT('1TRY INSUFFICIENT WORK AREA.')
 1200 FORMAT('1TRY CONSTANT Y.')
 1210 FORMAT('1TRY NEGATIVE WEIGHTS.')
 1220 FORMAT('1TRY ALL WEIGHTS ZERO (AND CONSTANT Y).')
 1230 FORMAT (///)
 1240 FORMAT (' CALL TO STAT')
 1250 FORMAT (' CALL TO STATS')
 1260 FORMAT(45H1TEST3.  TRY TURNING OFF THE PRINT FOR THOSE ,
     +   24HROUTINES WHICH ALLOW IT.)
 1270 FORMAT(37H TRY TURNING THE PRINT OFF FOR STATS.)
 1280 FORMAT(38H TRY TURNING THE PRINT OFF FOR STATWS.)
 1300 FORMAT(52H1TEST 4.  MAKE WORKING RUNS OF ALL ROUTINES TO CHECK,
     +   16H THE STATISTICS.)
 1310 FORMAT('1RUN STAT ON THE DAVIS-HARRISON PIKES PEAK DATA.')
 1320 FORMAT('1RUN STATW ON THE DAVIS-HARRISON PIKES PEAK DATA.')
 1340 FORMAT('1RUN STATS ON THE DAVIS-HARRISON PIKES PEAK DATA.')
 1350 FORMAT('1RUN STATWS ON THE DAVIS-HARRISON PIKES PEAK DATA.')
 1370 FORMAT('1RUN STATW ON THE DAVIS-HARRISON PIKES PEAK DATA.',
     +  '  WEIGHTS ALL EQUAL TO ONE.  COMPARE TO STAT ABOVE, NOT TO',
     +  ' STATW.')
 1380 FORMAT(42H SERIES WITH NONZERO VALUES WEIGHTED ZERO.)
 1390 FORMAT(/8H DATA = , 10F7.4)
 1400 FORMAT (14H CALL TO STATW)
 1410 FORMAT (15H CALL TO STATWS)
 1440 FORMAT ('1')
      END

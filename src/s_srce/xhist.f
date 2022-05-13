*XHIST
      SUBROUTINE XHIST(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS PROGRAM TESTS FEATURES OF THE HIST FAMILY TO ENSURE THAT
C     ALL ASPECTS OF THE HIST FAMILY ROUTINES WORK CORRECTLY.
C
C     WRITTEN BY  -  JOHN E. KOONTZ, JANET R. DONALDSON
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
     +   YLB,YUB
      INTEGER
     +   I,IPRT,LDSMIN,N,NCELL,NCONST,NPRTOF,NPRTON
C
C  LOCAL ARRAYS
      REAL
     +   Y(84),YCONST(10),YLONG(200),YPATH(10)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL HIST,HISTC,IPRINT,LDSCMP,NRAND
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ANINT,LOG10,MIN,NINT,REAL
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
C        A LOOP INDEX.
C     INTEGER IERR
C        FLAG TO INDICATE PRESENCE OF ERROR DETECTED BY PRECEDING
C        STARPAC CALL.  (0 IS OK, 1 IS ERROR)
C     INTEGER IPRT
C        LOGICAL OUTPUT UNIT.
C     INTEGER LDSMIN
C        THE MINIMUM AMOUNT OF WORK AREA NEEDED FOR A GIVEN PROBLEM.
C     INTEGER LDSTAK
C        AMOUNT OF WORK AREA.  SIZE OF DSTAK.
C     INTEGER N
C        THE LENGTH OF THE VECTOR Y.
C     INTEGER NCELL
C        THE USER SUPPLIED VALUE FOR THE NUMBER OF CELLS IN THE
C        HISTOGRAM.  IF NCELL IS LESS THAN OR EQUAL TO ZERO, THE
C        NUMBER OF CELLS TO BE USED (NCELLS) WILL BE CALCULATED FROM
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NCONST
C        LENGTH OF THE VECTOR YCONST.
C     INTEGER NPRTOF
C        FLAG FOR NO OUTPUT (EXCEPT ERROR MESSAGES).
C     INTEGER NPRTON
C        FLAG FOR FULL PRINTOUT.
C     REAL Y(84)
C        DATA VECTOR FOR TESTS.
C     REAL YCONST(10)
C        VECTOR OF CONSTANT DATA.
C     REAL YLB
C        THE LOWER BOUND FOR SELECTING DATA FROM Y FOR THE HISTOGRAM.
C     REAL YLONG(200)
C        LONG VECTOR OF DATA
C     REAL YPATH(10)
C        A VECTOR OF Y VALUES DESIGNED TO FORCE DIFFERENT PATHS
C        THROUGH THE SUMMATION ROUTINES.
C     REAL YUB
C        THE UPPER BOUND FOR SELECTING DATA FROM Y FOR THE HISTOGRAM.
C
C     DATA INITIALIZATIONS.
C
      DATA N /84/
      DATA NCONST /10/
      DATA NPRTON /1/
      DATA NPRTOF /0/
      DATA NCELL/10/
      DATA YLB/0.60E0/, YUB/0.63E0/
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
C
C     DEFINE IPRT, THE CURRENT OUTPUT UNIT.
C
      CALL IPRINT(IPRT)
C
C     CHECK FOR SUFFICIENT WORK AREA LENGTH.
C
      IF (LDSTAK.LT.300) THEN
        WRITE (IPRT, 1000)
         RETURN
      END IF
C
      DO 20 I=1,NCONST
         YCONST(I) = 1.0E0
   20 CONTINUE
C
C     HEADING.
C
      WRITE (IPRT,1150)
C
C     TEST 1.  CHECK ALL ERROR MESSAGES.
C
      WRITE (IPRT,1160)
C
C     ERROR 1, ZERO OR FEWER ELEMENTS.
C
      WRITE (IPRT,1180)
      CALL HIST(Y, 0, LDSTAK)
      WRITE (IPRT, 1350)
      WRITE (IPRT, 1360) (Y(I), I = 1, N)
      WRITE (IPRT,1170) IERR
      CALL HISTC(Y, 0, NCELL, YLB, YUB, LDSTAK)
      WRITE (IPRT, 1350)
      WRITE (IPRT, 1360) (Y(I), I = 1, N)
      WRITE (IPRT,1170) IERR
C
C     ERROR 2, NOT ENOUGH SPACE IN CSTAK.
C
      WRITE (IPRT,1190)
      CALL LDSCMP(2, 0, N, 0, 0, 0, 'S',
     +            MIN(NINT(5.5+1.5*ANINT(LOG10(REAL(N)))),25),LDSMIN)
      CALL HIST(Y, N, LDSMIN-1)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1195)
      CALL HIST(Y, N, LDSMIN)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1190)
      CALL LDSCMP(2, 0, N, 0, 0, 0, 'S', NCELL, LDSMIN)
      CALL HISTC(Y, N, NCELL, YLB, YUB, LDSMIN-1)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1195)
      CALL HISTC(Y, N, NCELL, YLB, YUB, LDSMIN)
      WRITE (IPRT,1170) IERR
C
C     CONSTANT Y. (NOT AN ERROR)
C
      WRITE (IPRT,1200)
      CALL HIST(YCONST, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT,1200)
      CALL HISTC(YCONST, NCONST, NCELL, YLB, YUB, LDSTAK)
      WRITE (IPRT,1170) IERR
C
C     ERROR 4, NO DATA WITHIN USER SUPPLIED LIMITS
C
      WRITE (IPRT, 1110)
      CALL HISTC(Y, N, 0, 4.0E0, 10.0E0, LDSTAK)
      WRITE (IPRT, 1170) IERR
C
C     TEST 2.  MAKE A WORKING RUN OF EACH ROUTINE TO CHECK
C     THE OUTPUT.
C
      WRITE (IPRT,1300)
      WRITE (IPRT,1310)
      CALL HIST(Y, N, LDSTAK)
      WRITE (IPRT, 1350)
      WRITE (IPRT, 1360) (Y(I), I = 1, N)
      WRITE (IPRT,1170) IERR
C
      WRITE (IPRT,1340)
      CALL HISTC(Y, N, NCELL, YLB, YUB, LDSTAK)
      WRITE (IPRT, 1350)
      WRITE (IPRT, 1360) (Y(I), I = 1, N)
      WRITE (IPRT,1170) IERR
C
C     RUN DATA SET 6.7.
C
      DO 90 I=1,10
         YPATH(I) = 0.0E0
   90 CONTINUE
      YPATH(1) = -1.0E0
      YPATH(10) = 1.0E0
      WRITE (IPRT,1130)
      CALL HIST(YPATH, NCONST, LDSTAK)
      WRITE (IPRT,1170) IERR
      WRITE (IPRT, 1130)
      CALL HISTC(YPATH, NCONST, 0, 0.0E0, 0.0E0, LDSTAK)
      WRITE (IPRT, 1130)
      CALL HISTC(YPATH, NCONST, 1, 0.0E0, 0.0E0, LDSTAK)
      WRITE (IPRT, 1130)
      CALL HISTC(YPATH, NCONST, 0, -0.5E0, 0.5E0, LDSTAK)
      WRITE (IPRT, 1130)
      CALL HISTC(YPATH, NCONST, 0, 1.0E0, 4.0E0, LDSTAK)
C
C     RUN DATA SET 6.8
C
      WRITE (IPRT, 1120)
      CALL NRAND (YLONG, 200, 3254767)
      CALL HIST (YLONG, 200, LDSTAK)
      RETURN
C
C     FORMATS
C
 1000 FORMAT ('1THE DIMENSION OF DSTAK AND THE VALUE OF LDSTAK NEEDED'/
     +  ' FOR HISTX MUST EQUAL OR EXCEED 300.  CHANGE DRIVER'/
     +  ' AND RECALL HISTX.')
 1110 FORMAT (41H1TRY NO DATA WITHIN USER SUPPLIED LIMITS.)
 1120 FORMAT (38H1RUN HIST ON 200 PSEUDO-RANDON NUMBERS)
 1130 FORMAT(24H1RUN HIST ON -1, 8*0, 1.)
 1150 FORMAT (48H1TEST RUNS FOR THE HISTOGRAM FAMILY OF ROUTINES.)
 1160 FORMAT(47H TEST 1.  GENERATE ONE OF EACH OF THE POSSIBLE ,
     +   15HERROR MESSAGES.)
 1170 FORMAT(22H THE VALUE OF IERR IS , I4)
 1180 FORMAT(28H TRY ZERO OR FEWER ELEMENTS.)
 1190 FORMAT('1TEST WITH INSUFFICIENT WORK AREA')
 1195 FORMAT(' TEST WITH EXACTLY THE RIGHT AMOUNT OF WORK AREA.')
 1200 FORMAT('1TRY CONSTANT Y. (NOT AN ERROR)')
 1300 FORMAT(52H1TEST 4.  MAKE WORKING RUNS OF ALL ROUTINES TO CHECK,
     +   12H THE OUTPUT.)
 1310 FORMAT(48H1RUN HIST ON THE DAVIS-HARRISON PIKES PEAK DATA.)
 1340 FORMAT(49H1RUN HISTC ON THE DAVIS-HARRISON PIKES PEAK DATA.)
 1350 FORMAT(/48H PRINT THE DATA TO INSURE THE ORIGINAL ORDER HAS,
     +   15H BEEN RESTORED.)
 1360 FORMAT (7F10.5)
      END

*FITPT2
      SUBROUTINE FITPT2 (SDRES, PV, WT, N, NNZW, WEIGHT, RES, RSS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE, ADAPTED FROM OMNITAB II, PRINTS
C     THE FOUR STANDARDIZED RESIDUAL PLOTS.
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
      REAL
     +   RSS
      INTEGER
     +   N,NNZW
      LOGICAL
     +   WEIGHT
C
C  ARRAY ARGUMENTS
      REAL
     +   PV(N),RES(N),SDRES(N),WT(N)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   ANNZW,DOT,FAC1,FAC2,FPLM,GAMMA,PI,PVDIV,PVMAX,PVMID,PVMIN,
     +   RATIO,ROWDIV,ROWMAX,ROWMID,ROWMIN,W,XDIV,XMAX,XMIN,YLABEL,
     +   YMAX,YMIN
      INTEGER
     +   I,I1,I2,IDOT,IFIRST,IMID,IPLOT,IPRB,IPRT,IPV,IROW,IX,K,L,
     +   NCOL,NCOLP1,NCOLPL,NCOLT2,NDOT,NROW
      CHARACTER
     +   IBLANK*1,IMINUS*1,IPLUS*1,ISTAR*1
C
C  LOCAL ARRAYS
      CHARACTER
     +   LINE(102)*1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      LOGICAL
     +   MVCHK
      EXTERNAL R1MACH,MVCHK
C
C  EXTERNAL SUBROUTINES
      EXTERNAL GETPI,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC INT,MAX,MIN,MOD
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ANNZW
C        THE NUMBER OF NONZERO WEIGHTS, USED IN COMPUTING
C        THE NORMAL PROBABILITY PLOT.
C     REAL DOT
C        ...
C     REAL FAC1, FAC2
C        FACTORS USED IN COMPUTING THE NORMAL PROBABILITY PLOT.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     REAL GAMMA
C        A VALUE USED IN COMPUTING THE NORMAL PROBABILITY PLOT.
C     INTEGER I
C        AN INDEX VARIABLE.
C     CHARACTER*1 IBLANK
C        THE VALUE OF THE CHARACTER -BLANK-.
C     INTEGER IERR
C        THE INTEGER VALUE DESIGNATING WHETHER ANY ERRORS HAVE
C        BEEN DETECTED.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .NE. 0, ERRORS HAVE BEEN DETECTED.
C     INTEGER IFIRST
C        THE FIRST ROW OF THE VARIABLES TO BE PLOTTED.
C     INTEGER IMID
C        THE MIDPOINT OF THE FIRST PLOT OF THE SECOND SET
C     CHARACTER*1 IMINUS
C        THE CHARACTER MINUS.
C     INTEGER IPLOT
C        AN INDICATOR VARIABLE DESIGNATING WHETHER THE FIRST OR
C        SECOND SET OF TWO PLOTS ARE BEING PRINTED.
C     CHARACTER*1 IPLUS
C        THE CHARACTER PLUS.
C     INTEGER IPRB
C        THE LOCATION IN THE PLOT STRING OF THE SYMBOL FOR THE
C        PROBABILITY PLOT.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IPV
C        THE LOCATION IN THE PLOT STRING OF THE SYMBOL FOR THE PLOT
C        VERSUS PREDICTED VALUE.
C     INTEGER IROW
C        THE ROW OF THE VARIABLES BEING PLOTTED.
C     CHARACTER*1 ISTAR
C        THE CHARACTER STAR.
C     INTEGER IX
C        THE LOCATION IN THE PLOT STRING OF THE SYMBOL FOR THE PLOTS
C        VERSUS THE INDEPENDENT VARIABLE.
C     INTEGER I1, I2
C        ...
C     INTEGER K, L
C        INDEX VARIABLES.
C     CHARACTER*1 LINE(102)
C        THE SYMBOLS (BLANKS AND CHARACTERS) FOR A GIVEN LINE
C        OF THE PLOT.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN EACH COLUMN OF DATA.
C     INTEGER NCOL, NCOLPL, NCOLP1, NCOLT2
C        THE NUMBER OF COLUMNS IN THE PLOT, NCOL+L, NCOL+1,
C        AND NCOL * 2.
C     INTEGER NDOT
C        ...
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NROW
C        THE NUMBER OF COLUMNS IN THE PLOT.
C     REAL PI
C        THE VALUE OF PI.
C     REAL PV(N)
C        THE PREDICTED VALUE BASED ON THE CURRENT COEFFICIENT ESTIMATES
C     REAL PVDIV
C        THE VALUE OF A DIVISION ALONG THE -PREDICTED VALUE- AXIS.
C     REAL PVMAX
C        THE LARGEST VALUE IN THE VECTOR PV.
C     REAL PVMID
C        THE MIDPOINT OF THE RANGE OF VALUES IN THE VECTOR PV.
C     REAL PVMIN
C        THE SMALLEST VALUE IN THE VECTOR PV.
C     REAL RATIO
C        A VALUE USED TO PRODUCE THE NORMAL PROBABILITY PLOT.
C     REAL RES(N)
C        THE RESIDUALS FROM THE FIT.
C     REAL ROWDIV
C        THE VALUE OF A DIVISION ALONG THE -ROW- AXIS.
C     REAL ROWMAX
C        THE LARGEST ROW VALUE.
C     REAL ROWMID
C        THE MIDPOINT OF THE RANGE OF THE ROWS PLOTTED.
C     REAL ROWMIN
C        THE SMALLEST ROW VALUE PLOTTED.
C     REAL RSS
C        THE RESIDUAL SUM OF SQUARES.
C     REAL SDRES(N)
C        THE STANDARD DEVIATIONS OF THE RESIDUALS.
C     REAL W
C        THE VALUE OF THE WEIGHT FOR THE CURRENT VALUE BEING PLOTTED.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     REAL WT(N)
C        THE USER SUPPLIED WEIGHTS.
C     REAL XDIV
C        THE VALUE OF A DIVISION ALONG THE X AXIS.
C     REAL XMAX
C        THE LARGEST VALUE ALONG THE X AXIS.
C     REAL XMIN
C        THE SMALLEST VALUE ALONG THE X AXIS.
C     REAL YLABEL
C        THE LABEL TO BE PRINTED ALONG THE Y AXIS.
C     REAL YMAX
C        THE LARGEST VALUE ALONG THE Y AXIS
C     REAL YMIN
C        THE SMALLEST VALUE ALONG THE Y AXIS.
C
      DATA IPLUS/'+'/, IMINUS/'-'/, ISTAR/'*'/, IBLANK/' '/
C
      CALL IPRINT(IPRT)
C
      FPLM = R1MACH(2)
C
C     CHECK FOR INSUFFICIENT POINTS TO PLOT
C
      IF (IERR.EQ.4) THEN
         DO 1 I = 1, N
            IF (SDRES(I).NE.FPLM) GO TO 5
    1    CONTINUE
         WRITE (IPRT, 1090)
         RETURN
      END IF
    5 CONTINUE
C
C     INITIALIZE VARIABLES FOR PROBABILITY PLOT
C
      CALL GETPI(PI)
      GAMMA = PI/8.0E0
      ANNZW = NNZW
      FAC1 = 1.0E0 / (ANNZW - 2.0E0*GAMMA + 1.0E0)
      FAC2 = 10.0E0
C
C     INITIALIZE THE PLOT SIZE (IN PLOT UNITS)
C
      NROW = 26
      NCOL = 51
      NCOLP1 = NCOL + 1
      NCOLT2 = 2*NCOL
      IMID = (NCOL-1)/2
C
C     FIND THE FIRST ROW OF OBSERVATIONS WITH NONZERO WEIGHTS
C
      IFIRST = 1
      IF (.NOT. WEIGHT) GO TO 20
      DO 10 I=1,N
         IF (WT(I).LE.0.0E0) GO TO 10
         IFIRST = I
         GO TO 20
   10 CONTINUE
C
C     BEGIN COMPUTATIONS FOR FIRST SET OF PLOTS
C
   20 IPLOT = 1
C
C     SET X AXIS LIMITS FOR STANDARDIZED RESIDUAL VS ROW PLOT,
C     AND STANDARDIZED RESIDUALS VS PREDICTED VALUES PLOT.
C
      ROWMIN = IFIRST
      PVMIN = PV(IFIRST)
      PVMAX = PV(IFIRST)
      ROWMAX = IFIRST
      DO 30 I=IFIRST,N
         W = 1.0E0
         IF (WEIGHT) W = WT(I)
         IF (W.GT.0.0E0) THEN
            ROWMAX = I
            IF (PV(I).LT.PVMIN) PVMIN = PV(I)
            IF (PV(I).GT.PVMAX) PVMAX = PV(I)
         END IF
   30 CONTINUE
C
      IF (PVMIN.LT.PVMAX) GO TO 35
         IF (PVMIN.EQ.0.0E0) GO TO 33
            PVMIN = PVMIN - PVMIN/2.0E0
            PVMAX = PVMAX + PVMAX/2.0E0
         GO TO 35
   33    CONTINUE
            PVMIN = -0.5E0
            PVMAX = 0.5E0
   35 CONTINUE
C
      ROWMID = (ROWMAX+ROWMIN)/2.0E0
      ROWDIV = (ROWMAX-ROWMIN)/(NCOL-1)
      PVMID = (PVMAX+PVMIN)/2.0E0
      PVDIV = (PVMAX-PVMIN)/(NCOL-1)
C
C     PRINT TITLES FOR FIRST PLOTS
C
      WRITE (IPRT,1000)
      GO TO 90
C
C     BEGIN COMPUTATIONS FOR SECOND SET OF PLOTS
C
   40 IPLOT = 2
C
C     SET AXIS LIMITS FOR THE STANDARDIZED RESIDUALS VS
C     STANDARDIZED RESIDUALS LAGED BY ONE AND FOR PROBABILITY PLOT
C
      XMIN = -3.75E0
      XMAX = 3.75E0
      XDIV = (XMAX-XMIN)/(NCOL-1)
C
C     PRINT TITLES FOR SECOND PLOTS
C
      WRITE (IPRT,1050)
C
C     WRITE FIRST LINE OF PLOTS
C
   90 CONTINUE
C
C     PRINT PLOTS, ONE LINE AT A TIME
C
      YLABEL = 3.75E0
      YMAX = FPLM
      YMIN = 4.05E0
      DO 160 K=1,NROW
         YMIN = YMIN - 0.3E0
         IF (-3.70E0.GE.YMIN) YMIN = -FPLM
         DO 100 L=1,NCOL
            NCOLPL = L + NCOL
            LINE(L) = IBLANK
            LINE(NCOLPL) = IBLANK
            IF ((K.NE.1) .AND. (K.NE.NROW)) GO TO 100
               LINE(L) = IMINUS
               LINE(NCOLPL) = IMINUS
               IF ((MOD(L,10).NE.1) .AND. (L.NE.1+NCOL/2)) GO TO 100
                  LINE(L) = IPLUS
                  LINE(NCOLPL) = IPLUS
  100    CONTINUE
         DO 110 I=1,N
            IF (WEIGHT) THEN
               W = WT(I)
            ELSE
               W = 1.0E0
            END IF
            IF ((W.NE.0.0E0) .AND. (.NOT.MVCHK(SDRES(I),FPLM))) THEN
               IF ((SDRES(I).GT.YMIN) .AND. (SDRES(I).LE.YMAX)) THEN
                  IF (IPLOT.EQ.1) THEN
C
C     SET PLOT LINE FOR FIRST SET OF PLOTS
C
                     IROW = INT(((I-ROWMIN)/ROWDIV)+1.5E0)
                     LINE(IROW) = ISTAR
                     IPV = INT((PV(I)-PVMIN)/PVDIV+1.5E0) + NCOL
                     LINE(IPV) = ISTAR
                  ELSE
C
C     SET PLOT LINE FOR PROBABILITY PLOT
C
                     RATIO = (ANNZW-GAMMA) * FAC1
                     IPRB = INT(4.91E0*(RATIO**0.14E0-
     +                         (1.0E0-RATIO)**0.14E0)*FAC2) + 77
                     IF (IPRB.LE.NCOL) IPRB = NCOL+1
                     IF (IPRB.GE.103) IPRB = 102
                     LINE(IPRB) = ISTAR
                     ANNZW = ANNZW - 1.0E0
                     IF ((ANNZW.LT.2.0E0) .AND. (NNZW.LE.10)) THEN
                        GAMMA = 1.0E0/3.0E0
                     END IF
                  END IF
               END IF
            END IF
  110    CONTINUE
C
C     SET PLOT LINE FOR CORRELATION PLOT
C
         IF (IPLOT.EQ.2) THEN
            IF (K.LE.N-1) THEN
              DOT = 0.0E0
              IF (WEIGHT) THEN
                NDOT = 0
                DO 120 IDOT = 1, N-K
                  IF ((WT(IDOT).GT.0.0E0) .AND.
     +                (WT(IDOT+K).GT.0.0E0)) THEN
                    NDOT = NDOT + 1
                    DOT = DOT + RES(IDOT)*RES(IDOT+K)
                  END IF
  120           CONTINUE
                IF (NDOT.GE.1) THEN
                   DOT = DOT * (N-K) / NDOT
                END IF
              ELSE
                DO 130 IDOT = 1, N-K
                  DOT = DOT + RES(IDOT)*RES(IDOT+K)
  130           CONTINUE
              END IF
              IX = INT(IMID*DOT/RSS) + IMID + 1
              I1 = MIN(IX,IMID+1)
              I2 = MAX(IX,IMID+1)
              DO 140 IX=I1,I2
                LINE(IX) = ISTAR
  140         CONTINUE
            END IF
         END IF
         IF (MOD(K,5).EQ.1) THEN
            IF (IPLOT.EQ.1) THEN
               WRITE (IPRT,2020) YLABEL, (LINE(L),L=1,NCOL), YLABEL,
     +         (LINE(L),L=NCOLP1,NCOLT2)
            ELSE
               WRITE (IPRT,1020) K, (LINE(L),L=1,NCOL), YLABEL,
     +         (LINE(L),L=NCOLP1,NCOLT2)
            END IF
            YLABEL = YLABEL - 1.5
         ELSE
            WRITE (IPRT,1030) (LINE(L),L=1,102)
         END IF
         YMAX = YMIN
  160 CONTINUE
C
C     PRINT BOTTOM LINE OF GRAPHS
C
      IF (IPLOT.EQ.1) THEN
C
C     PRINT X AXIS LABELS FOR FIRST SET OF PLOTS
C
         WRITE (IPRT,1040) ROWMIN, ROWMID, ROWMAX, PVMIN, PVMID, PVMAX
         GO TO 40
      ELSE
C
C     PRINT X AXIS LABELS FOR SECOND SET OF PLOTS
C
         WRITE (IPRT,1070)
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/20X, 23H STD RES VS ROW NUMBER , 35X,
     +   29H STD RES VS PREDICTED VALUES )
C1010 FORMAT (7X, 2('+', 9A1), '+', 4A1, 'X', 4A1, 2('+', 9A1), '+',
C    *   10X, 2('+', 9A1), '+', 4A1, 'X', 4A1, 2('+', 9A1), '+')
 1020 FORMAT (1X, I5, '+', 51A1, '+', 3X, F5.2, '+', 51A1, '+')
 1030 FORMAT (6X, '-', 51A1, '-', 8X, '-', 51A1, '-')
 1040 FORMAT (1X, F8.1, 17X, F8.1, 17X, F8.1, 4X, G11.4, 14X, G11.4,
     +   10X, G11.4)
 1050 FORMAT (/13X, 'AUTOCORRELATION FUNCTION OF RESIDUALS',
     +   23X, 36H NORMAL PROBABILITY PLOT OF STD RES )
C1060 FORMAT ('+', F5.2, '+', 51A1, '+', 3X, F5.2, '+', 51A1, '+')
 1070 FORMAT (4X, 5H-1.00, 22X, 3H0.0, 21X, 4H1.00, 5X, 4H-2.5, 23X,
     +   3H0.0, 22X, 3H2.5)
C1080 FORMAT ('+', 6X, 2('+', 9A1), '+', 4A1, 'X', 4A1, 2('+', 9A1),
C    *   '+', 10X, 2('+', 9A1), '+', 4A1, 'X', 4A1, 2('+', 9A1), '+')
 1090 FORMAT (// 1X, 13('*')/ 1X, 13H*  WARNING  */ 1X, 13('*')//
     +   54H THE STANDARDIZED RESIDUAL PLOTS HAVE BEEN SUPPRESSED.,
     +   45H  NONE OF THE STANDARDIZED RESIDUALS COULD BE,
     +   10H COMPUTED,/
     +   50H BECAUSE FOR EACH OBSERVATION EITHER THE WEIGHT OR,
     +   48H THE STANDARD DEVIATION OF THE RESIDUAL IS ZERO.)
 2020 FORMAT (1X, F5.2, '+', 51A1, '+', 3X, F5.2, '+', 51A1, '+')
      END

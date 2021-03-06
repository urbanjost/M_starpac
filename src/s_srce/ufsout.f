*UFSOUT
      SUBROUTINE UFSOUT(XAXIS, YAXIS, ISYM, NPTS, BW, IDF, LAG, LAGLST,
     +   NEWPG, ISPCER, LPCV, XPLTMN, XPLTMX, YPLTMN, YPLTMX, ILOG,
     +   PHAS, FREQ, NF, UNIVAR, NMSUB)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRODUCES THE FOURIER BIVARIATE SPECTRUM OUTPUT.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JUNE 10, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   BW,XPLTMN,XPLTMX,YPLTMN,YPLTMX
      INTEGER
     +   IDF,ILOG,ISPCER,LAG,LAGLST,LPCV,NF,NPTS
      LOGICAL
     +   NEWPG,UNIVAR
C
C  ARRAY ARGUMENTS
      REAL
     +   FREQ(NF),PHAS(NF),XAXIS(LPCV),YAXIS(LPCV)
      INTEGER
     +   ISYM(LPCV)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   PI,XMN,XMX,YMN,YMX
      INTEGER
     +   I,IPRT
      LOGICAL
     +   ERROR
C
C  EXTERNAL SUBROUTINES
      EXTERNAL GETPI,IPRINT,PPLMT,PPMN,VERSP
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL BW
C        THE BANDWIDTH.
C     LOGICAL ERROR
C        AN ERROR FLAG
C     REAL FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER IDF
C        THE EFFECTIVE DEGREES OF FREEDOM.
C     INTEGER IERR
C        THE ERROR FLAG.
C     INTEGER ILOG
C        A CODE USED TO SPECIFY THE TYPE OF PLOT, WHERE IF
C        ILOG = 0 THE PLOT IS LINEAR/LINEAR, IF
C        ILOG = 1 THE PLOT IS LOG/LINEAR, IF
C        ILOG = 2 THE PLOT IS LINEAR/LOG, AND IF
C        ILOG = 3 THE PLOT IS LOG/LOG.
C     INTEGER IPRT
C        THE LOGICAL UNIT NUMBER FOR THE OUTPUT.
C     INTEGER ISPCER
C        A VARIABLE USED TO DESIGNATE AN ERROR IN THE SPECTRUM
C        VALUES.
C     INTEGER ISYM(LPCV)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT.
C     INTEGER LAGLST
C        THE LAST LAG BEFORE MISSING DATA CAUSED THE ACVF OF EITHER
C        SERIES 1 OR 2 NOT TO BE COMPUTED.
C     INTEGER LPCV
C        THE LENGTH OF THE VECTORS USED FOR PLOTTING.
C     LOGICAL NEWPG
C        THE LOGICAL VARIABLE USED TO DETERMINE IF OUTPUT
C        WILL BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE CALLING SUBROUTINE.
C     INTEGER NPTS
C        THE NUMBER OF CO-ORDINATES TO BE PLOTTED.
C     REAL PHAS(NF)
C        THE PHASE COMPONENT OF THE BIVARIATE SPECTRA.
C     REAL PI
C        THE VALUE OF PI.
C     LOGICAL UNIVAR
C        THE LOGICAL VARIABLE USED TO DETERMINE IF THE OUTPUT
C        IS FOR UNIVARIATE (TRUE) OR BIVARIATE (FALSE) SPECTRA.
C     REAL XAXIS(LPCV)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     REAL XMN, XMX
C        *
C     REAL XPLTMN, XPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE X AXIS.
C     REAL YAXIS(LPCV)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C     REAL YMN, YMX
C        *
C     REAL YPLTMN, YPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE Y AXIS.
C
C      SET LOGICAL UNIT NUMBER FOR OUTPUT AND SET OUTPUT WIDTH.
C
      CALL IPRINT(IPRT)
C
      CALL GETPI(PI)
C
      IF (NEWPG) WRITE (IPRT,1010)
      IF (ISPCER.LE.1) GO TO 10
      CALL VERSP(.TRUE.)
      WRITE (IPRT,1060) LAGLST, LAG
      RETURN
C
   10 CONTINUE
      CALL VERSP(.TRUE.)
      IF (.NOT.UNIVAR) WRITE (IPRT,1070)
      IF (UNIVAR) WRITE (IPRT,1080)
      WRITE (IPRT,1020) LAG, BW, IDF
      IF (ISPCER.EQ.0) GO TO 20
      WRITE (IPRT,1050)
      GO TO 30
C
   20 CONTINUE
C
C     PRINT PLOTS
C
C     PLOT SQUARED COHERENCY COMPONENT OF SPECTRUM
C
      CALL PPLMT(YAXIS, YAXIS, XAXIS, XAXIS, NPTS, 1, LPCV, YPLTMN,
     +  YPLTMX, YMN, YMX, XPLTMN, XPLTMX, XMN, XMX, ERROR, NMSUB,
     +  .FALSE.)
      IF (.NOT.ERROR)
     +  CALL PPMN(YAXIS, YAXIS, XAXIS, XAXIS, NPTS, 1, LPCV, 1, ISYM,
     +  LPCV, 0, -1, YMN, YMX, XMN, XMX, .FALSE., ILOG)
      IF (XPLTMN.EQ.0.0E0 .AND. XPLTMX.EQ.0.5E0) WRITE (IPRT, 1030)
C
   30 IF (UNIVAR) RETURN
      DO 40 I=1,NF
         XAXIS(I) = FREQ(I)
         XAXIS(NF+I) = FREQ(I)
         YAXIS(I) = PHAS(I)
         IF (PHAS(I).GT.0.0E0) THEN
            YAXIS(NF+I) = PHAS(I) - 2*PI
         ELSE IF (PHAS(I).LT.0.0E0) THEN
            YAXIS(NF+I) = PHAS(I) + 2*PI
         ELSE
            YAXIS(NF+I) = 0.0E0
         END IF
   40 CONTINUE
C
C     PLOT SMOOTHED PHASE COMPONENT OF SPECTRUM
C
      WRITE (IPRT,1010)
      CALL VERSP(.TRUE.)
      WRITE (IPRT,1000)
      WRITE (IPRT,1020) LAG, BW, IDF
      CALL PPLMT(YAXIS, YAXIS, XAXIS, XAXIS, 2*NF, 1, 2*NF, -2*PI, 2*PI,
     +  YMN, YMX, XPLTMN, XPLTMX, XMN, XMX, ERROR, NMSUB, .FALSE.)
      IF (ERROR) THEN
        IERR = 1
      ELSE
        CALL PPMN(YAXIS, YAXIS, XAXIS, XAXIS,
     +            2*NF, 1, 2*NF, 0, ISYM, LPCV,
     +            0, -1, YMN, YMX, XMN, XMX, .FALSE., ILOG)
        IF (XPLTMN.EQ.0.0E0 .AND. XPLTMX.EQ.0.5E0) WRITE (IPRT, 1030)
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (50H -- SMOOTHED FOURIER SPECTRUM (PHASE COMPONENT) --)
 1010 FORMAT ('1')
 1020 FORMAT (45H    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=, I5, 1X,
     +   5H/ BW=, F6.4, 1X, 6H/ EDF=, I6, ')')
 1030 FORMAT (5H+FREQ/7H PERIOD, 9X, 3HINF, 7X, 3H20., 7X, 3H10., 8X,
     +   6H6.6667, 4X, 2H5., 8X, 2H4., 8X, 6H3.3333, 4X, 6H2.8571, 4X,
     +   3H2.5, 7X, 6H2.2222, 4X, 2H2.)
C1040 FORMAT (5H+FREQ/7H PERIOD, 9X, 3HINF, 7X, 3H10., 7X, 2H5., 8X,
C    *   6H3.3333, 4X, 3H2.5, 7X, 2H2.)
 1050 FORMAT (//39H THE PLOT HAS BEEN SUPRESSED BECAUSE NO/
     +   40H POSITIVE SPECTRUM VALUES WERE COMPUTED.)
 1060 FORMAT (//50H THE LARGEST LAG WINDOW TRUNCATION POINT WHICH CAN/
     +   12H BE USED IS , I5, '.'/34H THE SPECTRUM FOR THE REQUESTED LA,
     +   8HG WINDOW, 10H POINT OF , I5, ','/24H THEREFORE, CANNOT BE CO,
     +   7HMPUTED.)
 1070 FORMAT (48H -- SMOOTHED FOURIER SPECTRUM (SQUARED COHERENCY,
     +   46H COMPONENT) (+), 95 PCT. CONFIDENCE LIMITS (.),
     +   38H AND 95 PCT. SIGNIFICANCE LEVEL (-) --)
 1080 FORMAT (32H -- SMOOTHED FOURIER SPECTRUM --)
      END

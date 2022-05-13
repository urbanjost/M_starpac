*BFSDRV
      SUBROUTINE BFSDRV(Y1, Y2, YMISS1, YMISS2, CCOV, NLPPC, SPCF1,
     +   SPCF2, NF, FMIN, FMAX, FREQ, N, NW, LAGMAX, LAGS, LAGMX1,
     +   WORK, LWORK, DELTA, ISYM, XAXIS, YAXIS, LPCV, ALPHA, NPRT,
     +   WINDOW, ICCOV, JCCOV, M, INDEX1, INDEX2, CSPC2, PHAS, ICSPC2,
     +   IPHAS, CODD, CEVEN, W, LW, NMSUB, LDSMIN, LDSTAK, OPTION,
     +   NFFT, INLPPC, JNLPPC, LY)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE CONTROLING ROUTINE FOR TIME SERIES FOURIER
C     SPECTRUM ANALYSIS .
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   ALPHA,DELTA,FMAX,FMIN,YMISS1,YMISS2
      INTEGER
     +   ICCOV,ICSPC2,INDEX1,INDEX2,INLPPC,IPHAS,JCCOV,JNLPPC,
     +   LAGMAX,LAGMX1,LDSMIN,LDSTAK,LPCV,LW,LWORK,LY,M,N,NF,NFFT,
     +   NPRT,NW
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   CCOV(*),CEVEN(*),CODD(*),CSPC2(*),FREQ(*),PHAS(*),SPCF1(*),
     +   SPCF2(*),W(*),WORK(*),XAXIS(*),Y1(*),Y2(*),YAXIS(*)
      INTEGER
     +   ISYM(*),LAGS(*),NLPPC(*)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL WINDOW
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   BW,DF,FMN,FMX,XPLTMN,XPLTMX,YMEAN1,YMEAN2,YPLTMN,YPLTMX
      INTEGER
     +   I,ILOG,ISPCER,J,K,LAG,LAGLST,NFUSED,NPTS,NWUSED
      LOGICAL
     +   NEWPG,UNIVAR
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   LSTLAG
      EXTERNAL LSTLAG
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACVF,ACVFF,ACVFM,BFSER,BFSLAG,BFSMN,CCVF,CCVFF,CCVFM,
     +   DFBW,DFBWM,SETFRQ,UFSEST,UFSOUT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC INT,MAX,MIN
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C  STATEMENT FUNCTIONS
      INTEGER
     +   I3C,I3N
C
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION BW
C        THE BANDWIDTH.
C     DOUBLE PRECISION CCOV(ICCOV,JCCOV,M)
C        THE COVARIANCES.
C     DOUBLE PRECISION CEVEN(LAGMX1)
C        THE SUMS OF THE COVARIANCES FOR EACH LAG.
C     DOUBLE PRECISION CODD(LAGMX1)
C        THE DIFFERENCES OF THE COVARIANCES FOR EACH LAG.
C     DOUBLE PRECISION CSPC2(ICSPC2,NW)
C        THE SQUARED COHERENCY COMPONENT OF THE BIVARIATE SPECTRA.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION DF
C        THE EFFECTIVE DEGREES OF FREEDOM.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCES AT WHICH THE
C        SPECTRUM IS TO BE COMPUTED.
C     DOUBLE PRECISION FMN, FMX
C        *
C     DOUBLE PRECISION FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER INDEX1, INDEX2
C        THE INDICES OF THE COVARIANCES OF THE TWO SERIES.
C     INTEGER INLPPC
C        THE FIRST DIMENSION OF THE ARRAY NLPPC.
C     INTEGER IPHAS
C        THE FIRST DIMENSION OF THE ARRAY PHAS.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER ILOG
C        A CODE USED TO SPECIFY THE TYPE OF PLOT, WHERE IF
C        ILOG = 0 THE PLOT IS LINEAR/LINEAR, IF
C        ILOG = 1 THE PLOT IS LOG/LINEAR, IF
C     INTEGER ISPCER
C        AN ERROR FLAG USED FOR THE SPECTRUM PLOTS.
C     INTEGER ISYM(LPCV)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER I3C
C        STATEMENT FUNCTION FOR FINDING LOCATIONS WITHIN CCOV.
C     INTEGER I3N
C        STATEMENT FUNCTION FOR FINDING LOCATIONS WITHIN NLPPC.
C     INTEGER JCCOV
C        THE SECOND DIMENSION OF CCOV
C     INTEGER JNLPPC
C        THE SECOND DIMENSION OF NLPPC
C     INTEGER LAG
C        THE LAG WINDWO TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LAGLST
C        THE LAST LAG BEFORE MISSING DATA CAUSED AN ACVF
C        TO BE UNABLE TO BE COMPUTED.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGMX1
C        THE VALUE LAGMAX+1.
C     INTEGER LAGS(NW)
C        THE ARRAY USED TO STORE THE LAG WINDOW TRUCCATION
C        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     INTEGER LPCV
C        THE LENGTH OF THE VECTORS USED FOR PLOTTING.
C     INTEGER LWORK
C        THE LENGTH OF THE VECTOR W.
C     INTEGER LY
C        THE LENGTH OF THE VECTORS Y1 AND Y2.
C     INTEGER M
C        THE NUMBER OF SERIES FOR WHICH THE COVARIANCES WERE COMPUTED
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     LOGICAL NEWPG
C        THE LOGICAL VARIABLE USED TO DETERMINE IF OUTPUT
C        WILL BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     INTEGER NFUSED
C        THE NUMBER OF FREQUENCIES ACTUALLY USED.
C     INTEGER NLPPC(INLPPC,JNLPPC,M)
C         THE ARRAY CONTAINING THE NUMBER OF LAG PRODUCT PAIRS.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NPRT
C        A CODE USED TO SPECIFY THE TYPE OF PLOT.
C        IF NPRT < 0 THE PLOT IS DECIBLES/LINEAR
C        IF NPRT = 0 THE PLOT IS SUPPRESSED.
C        IF NPRT > 0 THE PLOT IS LOG/LINEAR
C     INTEGER NPTS
C        THE NUMBER OF X, Y CO-ORDINATES TO BE PLOTTED.
C     INTEGER NW
C        THE VARIABLE USED TO DETERMINE THE NUMBER OF DIFFERENT
C        BANDWIDTHS TO BE USED.
C     INTEGER NWUSED
C        THE NUMBER OF DIFFERENT BANDWIDTHS ACTUALLY USED.
C     LOGICAL OPTION(4)
C        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
C        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
C        OR NOT (FALSE).
C     DOUBLE PRECISION PHAS(IPHAS,NW)
C        THE PHASE COMPONENT OF THE BIVARIATE SPECTRUM.
C     DOUBLE PRECISION SPCF1(NF), SPCF2(NF)
C        THE ARRAYS IN WHICH THE SPECTRUM IS STORED.
C     LOGICAL UNIVAR
C        THE LOGICAL VARIABLE USED TO DETERMINE IF THE OUTPUT
C        IS FOR UNIVARIATE (TRUE) OR BIVARIATE (FALSE) SPECTRA.
C     DOUBLE PRECISION W(LW)
C        THE VECTOR OF WINDOWS.
C     EXTERNAL WINDOW
C        THE SUBROUTINE USED TO COMPUTE THE WINDOW.
C     DOUBLE PRECISION WORK(LWORK)
C        THE VECTOR OF WORK SPACE.
C     DOUBLE PRECISION XAXIS(LPCV)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION XPLTMN, XPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE X AXIS.
C     DOUBLE PRECISION YAXIS(LPCV)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOTS.
C     DOUBLE PRECISION YMEAN1, YMEAN2
C        THE MEAN OF THE OBSERVED TIME SERIES
C     DOUBLE PRECISION YMISS1, YMISS2
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IN THE SERIES IS MISSING.  IF Y(I) = YMISS,
C        THE VALUE IS ASSUMED MISSING, OTHERWISE IT IS NOT.
C     DOUBLE PRECISION YPLTMN, YPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE Y AXIS.
C     DOUBLE PRECISION Y1(N), Y2(N)
C         THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C
C
C  STATEMENT FUNCTION DEFINITIONS
C
      I3C(I,J,K) = I + (J-1)*ICCOV + (K-1)*JCCOV*ICCOV
      I3N(I,J,K) = I + (J-1)*INLPPC + (K-1)*JNLPPC*INLPPC
C
      NFUSED = NF
      IF (OPTION(4)) THEN
        FMN = MAX(FMIN, 0.0D0)
        FMX = MIN(FMAX, 0.5D0)
        IF (FMN.GE.FMX) THEN
          FMN = 0.0D0
          FMX = 0.5D0
        END IF
      ELSE
C
C       SET VARIOUS VALUES FOR SHORT FORMS OF CALL STATEMENT
C
        NPRT = -1
        FMN = 0.0D0
        FMX = 0.5D0
        LAGMX1 = LAGMAX + 1
      END IF
C
C     CHECK FOR ERRORS
C
      CALL BFSER(NMSUB, N, LAGMAX, ICCOV, JCCOV, INLPPC, JNLPPC, M,
     +           INDEX1, INDEX2, ICSPC2, IPHAS, NF, NW, LAGS,
     +           LDSTAK, LDSMIN, LY, NFFT, OPTION)
C
      IF (IERR.EQ.1) RETURN
C
C     SET VARIOUS PROGRAM PARAMETERS.
C
      ALPHA = 0.95D0
      DELTA = 1.0D0
C
C     COMPUTE COVARIANCES
C
      LAGLST = LAGMAX
      IF (OPTION(1)) THEN
        CALL ACVFF(Y1, N, NFFT, YMEAN1,
     +             CCOV(I3C(1,INDEX1,INDEX1)),
     +             LAGMAX, ICCOV, N, WORK, NFFT)
        CALL ACVFF(Y2, N, NFFT, YMEAN2,
     +             CCOV(I3C(1,INDEX2,INDEX2)),
     +             LAGMAX, ICCOV, N, WORK, NFFT)
        CALL CCVFF(Y1, Y2, N, NFFT, LAGMAX,
     +             CCOV(I3C(1,INDEX1,INDEX2)),
     +             CCOV(I3C(1,INDEX2,INDEX1)), ICCOV, N, WORK, LWORK)
C
      ELSE
        IF (OPTION(3)) THEN
          IF (OPTION(2)) LAGLST = LSTLAG(NLPPC,LAGMAX,INLPPC)
        ELSE
          IF (OPTION(2)) THEN
            CALL ACVFM(Y1, YMISS1, N, YMEAN1,
     +                 CCOV(I3C(1,INDEX1,INDEX1)),
     +                 LAGMAX, LAGLST, NLPPC, ICCOV)
            CALL ACVFM(Y2, YMISS2, N, YMEAN2,
     +                 CCOV(I3C(1,INDEX2,INDEX2)),
     +                 LAGMAX, LAGLST, NLPPC, ICCOV)
            CALL CCVFM(Y1, YMISS1, Y2, YMISS2, N, LAGMAX, YMEAN1,
     +                 YMEAN2, CCOV(I3C(1,INDEX1,INDEX2)),
     +                 CCOV(I3C(1,INDEX2,INDEX1)), ICCOV,
     +                 NLPPC(I3N(1,INDEX1,INDEX2)),
     +                 NLPPC(I3N(1,INDEX2,INDEX1)))
C
          ELSE
            CALL ACVF(Y1, N, YMEAN1, CCOV(I3C(1,INDEX1,INDEX1)), LAGMAX,
     +                ICCOV)
            CALL ACVF(Y2, N, YMEAN2, CCOV(I3C(1,INDEX2,INDEX2)), LAGMAX,
     +                ICCOV)
            CALL CCVF(Y1, Y2, N, LAGMAX, YMEAN1, YMEAN2,
     +                CCOV(I3C(1,INDEX1,INDEX2)),
     +                CCOV(I3C(1,INDEX2,INDEX1)), ICCOV)
          END IF
        END IF
      END IF
C
      IF (LAGLST.LE.0) THEN
C
C     AN ERROR HAS BEEN DETECTED
C
         IERR = 2
         RETURN
      END IF
C
C     COMPUTE THE VECTOR OF LAG WINDOW TRUNCATION POINTS, ORDERED
C     SMALLEST TO LARGEST.
C
      NWUSED = NW
      IF (.NOT.OPTION(4)) CALL BFSLAG(CCOV, LAGLST, LAGS, N, NW, NWUSED,
     +                                ICCOV, JCCOV, INDEX1, INDEX2)
C
C     BEGIN COMPUTING FOURIER SPECTRUM FOR SERIES
C
      UNIVAR = .FALSE.
C
      ILOG = 0
C
      XPLTMN = FMN
      XPLTMX = FMX
C
      YPLTMN = 0.0D0
      YPLTMX = 1.0D0
C
C     SET FREQUENCIES FOR THE SPECTRUM.
C
      CALL SETFRQ(FREQ, NF, 1, FMN, FMX, DELTA)
C
C     COMPUTE AND PLOT SPECTRUM VALUES.
C
      NEWPG = .FALSE.
C
C     COMPUTE THE EVEN AND ODD CCVF ESTIMATES
C
      CEVEN(1) = CCOV(I3C(1,INDEX1,INDEX2))
      CODD(1) = 0.0D0
      DO 30 I=1,LAGLST
         CEVEN(I+1) = 0.5D0*
     +                (CCOV(I3C(I+1,INDEX1,INDEX2))+
     +                 CCOV(I3C(I+1,INDEX2,INDEX1)))
         CODD(I+1) = 0.5D0*
     +               (CCOV(I3C(I+1,INDEX1,INDEX2))-
     +                CCOV(I3C(I+1,INDEX2,INDEX1)))
   30 CONTINUE
C
      DO 60 I=1,NWUSED
         LAG = LAGS(I)
         IF (LAG.GT.LAGLST) THEN
            ISPCER = 2
            DF = 0.0D0
         ELSE
C
            ISPCER = 0
C
C     COMPUTE THE WINDOW, AND EFFECTIVE DEGREES OF FREEDOM AND
C     BANDWIDTH BASED ON THE WINDOW
C
            CALL WINDOW(LAG, W, LW)
            IF (OPTION(2)) THEN
               CALL DFBWM(N, LAG, W, LW, NLPPC(I3N(1,INDEX1,INDEX2)),
     +                    NLPPC(I3N(1,INDEX2,INDEX1)), INLPPC, DF, BW)
            ELSE
               CALL DFBW(N, LAG, W, LW, DF, BW)
            END IF
C
C     COMPUTE THE SPECTRUM FOR EACH INDIVIDUAL SERIES
C
            CALL UFSEST(CCOV(I3C(1,INDEX1,INDEX1)), W, LAG, SPCF1,
     +                  NFUSED, ICCOV, LAGMAX, NF, FREQ, DELTA)
C
            CALL UFSEST(CCOV(I3C(1,INDEX2,INDEX2)), W, LAG, SPCF2,
     +                  NFUSED, ICCOV, LAGMAX, NF, FREQ, DELTA)
C
            CALL BFSMN(SPCF1, SPCF2, CEVEN, CODD, W, LW, LAG, DF, NPRT,
     +                 NF, CSPC2(1+(I-1)*ICSPC2), PHAS(1+(I-1)*IPHAS),
     +                 FREQ, NPTS, XAXIS,
     +                 YAXIS, ISYM, LPCV, ALPHA, LAGMX1, DELTA)
C
            IF (NPRT.EQ.0) GO TO 60
C
         END IF
         CALL UFSOUT(XAXIS, YAXIS, ISYM, NPTS, BW, INT(DF+0.5D0), LAG,
     +               LAGMAX, NEWPG, ISPCER, NFUSED+5, XPLTMN, XPLTMX,
     +               YPLTMN, YPLTMX, ILOG, PHAS(1+(I-1)*IPHAS), FREQ,
     +               NF, UNIVAR, NMSUB)
C
         NEWPG = .TRUE.
C
   60 CONTINUE
C
      RETURN
C
      END

*SPPLTL
      SUBROUTINE SPPLTL (SPCMN, SPCMX, ALOW, AUP, YPLTMN, YPLTMX,
     +   CILOW, CIMID, CIUP)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SETS VARIOUS Y AXIS LIMITS FOR DECIBLE
C     SPECTRUM PLOTS.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   ALOW,AUP,CILOW,CIMID,CIUP,SPCMN,SPCMX,YPLTMN,YPLTMX
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   RNGMN,YMAX,YMIN
C
C  INTRINSIC FUNCTIONS
      INTRINSIC LOG10
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ALOW, AUP
C        FACTORS USED TO COMPUTE THE CONFIDENCE INTERVALS.
C     DOUBLE PRECISION CILOW, CIMID, CIUP
C        THE Y CORDNATES FOR THE LOWER  MID AND UPPER CONFIDENCE
C        INTERVAL POINTS.
C     DOUBLE PRECISION RNGMN
C        THE MINIMUM Y AXIS RANGE FOR THE PLOT.
C     DOUBLE PRECISION SPCMN, SPCMX
C        THE MINIMUM AND MAXIMUM SPECTRAL VALUE TO BE PLOTTED.
C     DOUBLE PRECISION YMAX, YMIN
C        THE MAXIMUM AND MINIMUM ACTUAL SPECTRUM VALUE TO BE PLOTTED.
C     DOUBLE PRECISION YPLTMN, YPLTMX
C        THE MINIMUM AND MAXIMUM VAUES TO BE PLOTTED FOR THE Y AXIS.
C
C     SET CO-ORDINATES FOR DECIBLE PLOTS
C
      YMAX = LOG10(SPCMX)
      YMIN = LOG10(SPCMN)
C
      YPLTMX = SPCMX
      YPLTMN = SPCMN
      RNGMN = 2.0D0 * (LOG10(AUP) - LOG10(ALOW))
      IF (YMAX - YMIN .GE. RNGMN) GO TO 10
C
      YPLTMX = 10.0D0 ** (YMAX + (RNGMN - YMAX + YMIN) * 0.5D0)
      YPLTMN = 10.0D0 ** (YMIN - (RNGMN - YMAX + YMIN) * 0.5D0)
C
   10 CIUP = YPLTMX
      CIMID = CIUP / AUP
      CILOW = CIMID * ALOW
C
      RETURN
      END

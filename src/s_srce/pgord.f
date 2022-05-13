*PGORD
      SUBROUTINE PGORD (PER, NPTS, YAXIS, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRODUCES CO-ORDINATES FOR THE PERIODOGRAM PLOT.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   NPRT,NPTS
C
C  ARRAY ARGUMENTS
      REAL
     +   PER(NPTS),YAXIS(NPTS)
C
C  LOCAL SCALARS
      REAL
     +   FPLM
      INTEGER
     +   I
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  INTRINSIC FUNCTIONS
      INTRINSIC IABS,LOG10
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER NPRT
C        THE VARIABLE CONTROLING PRINTED OUTPUT, WHERE
C        IF NPRT .LE. -2, THE OUTPUT CONSISTS OF A PAGE PLOT OF THE
C                         PERIODOGRAM ON A LOG-LINEAR SCALE,
C        IF NPRT .EQ. -1, THE OUTPUT CONSISTS OF A PAGE PLOT OF THE
C                         PERIODOGRAM IN DECIBELS ON A LINEAR SCALE,
C        IF NPRT .EQ.  0, THE OUTPUT IS SUPPRESSED,
C        IF NPRT .EQ.  1, THE OUTPUT CONSISTS OF A VERTICAL PLOT OF THE
C                         PERIODOGRAM IN DECIBELS ON A LINEAR SCALE.
C        IF NPRT .GE.  2, THE OUTPUT CONSISTS OF A VERTICAL PLOT OF THE
C                         PERIODOGRAM ON A LOG-LINEAR SCALE,
C     INTEGER NPTS
C        THE NUMBER OF FREQUENCIES FOR WHICH THE SPECTRAL ESTIMATES
C        ARE ESTIMATED.
C     REAL PER(NPTS)
C        THE ARRAY CONTAINING THE PERIODOGRAM VALUES.
C     REAL YAXIS(NPTS)
C        THE Y CO-ORDINATES FOR THE PERIODOGRAM PLOTS.
C
      FPLM = R1MACH(2)
C
C     THE FIRST VALUE SHOULD BE ZERO, SO NO ATTEMPT IS MADE TO PLOT IT.
C
      YAXIS(1) = FPLM
C
      DO 10 I = 2, NPTS
         YAXIS(I) = FPLM
         IF (PER(I) .LE. 0.0E0) GO TO 10
            YAXIS(I) = PER(I)
            IF (IABS(NPRT) .EQ. 1) YAXIS(I) = 10.0E0*LOG10(YAXIS(I))
   10 CONTINUE
C
      RETURN
C
      END

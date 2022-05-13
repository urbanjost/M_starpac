*PLTPLX
      SUBROUTINE PLTPLX(POINT, YMN, SCALE, IPOINT, IEND)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE POINT LOCATION IN THE PLOT STRING.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JANUARY 21, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   POINT,SCALE,YMN
      INTEGER
     +   IEND,IPOINT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IEND
C        THE NUMBER OF LOCATIONS IN THE PLOT STRING.
C     INTEGER IPOINT
C        THE LOCATION IN THE PLOT STRING OF THE VALUE BEING PLOTTED.
C     REAL POINT
C        THE VALUE TO BE PLOTTED.
C     REAL SCALE
C        THE SCALE INTERVAL OF THE PLOT.
C     REAL YMN
C        THE GRAPH AXIS LOWER LIMITS ACTUALLY USED.
C
      IPOINT = (POINT-YMN)/SCALE + 2.5
      IF (IPOINT .LT. 2) IPOINT = 1
      IF (IPOINT .GT. IEND) IPOINT = IEND
      RETURN
      END

*NLERR
      SUBROUTINE NLERR (ICNVCD, ISKULL)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SETS THE ERROR FLAG IERR BASED ON THE CONVERGENCE
C     CODE RETURNED BY NL2.
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
      INTEGER
     +   ICNVCD
C
C  ARRAY ARGUMENTS
      INTEGER
     +   ISKULL(10)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER ICNVCD
C        THE CONVERGENCE CODE FROM NL2.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER ISKULL(10)
C        AN ERROR MESSAGE INDICATOR VARIABLE.
C
C     INITIALIZE MESSAGE INDICATOR VARIABLE
C
      DO 5 I = 1, 10
         ISKULL(I) = 0
    5 CONTINUE
C
C     SET ERROR FLAG
C
      GO TO (10, 10, 20, 20, 20, 20, 40, 50, 60, 60, 10, 30, 10, 10,
     +   10), ICNVCD
C
C     BAD VALUE
C
   10 IERR = 1
      RETURN
C
C     ACCEPTABLE STOPPING CONDITION
C
   20 IERR = 0
      RETURN
C
C     INITIAL VARIANCE COMPUTATION OVERFLOWS
C
   30 IERR = 2
      ISKULL(2) = 1
      RETURN
C
C     SINGULAR CONVERGENCE
C
   40 IERR = 3
      ISKULL(3) = 1
      RETURN
C
C     FALSE CONVERGENCE
C
   50 IERR = 5
      ISKULL(5) = 1
      RETURN
C
C     ITERATION OR FUNCTION EVALUATION LIMIT
C
   60 IERR = 6
      ISKULL(6) = 1
      RETURN
C
      END

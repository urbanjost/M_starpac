*FITEXT
      SUBROUTINE FITEXT(RSS, YSS, EXACT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE CHECKS WHETHER THE FIT IS EXACT TO MACHINE
C     PRECISION.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   RSS,YSS
      LOGICAL
     +   EXACT
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FPLRS,RSSTST
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   D1MACH
      EXTERNAL D1MACH
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL EXACT
C        AN INDICATOR VALUE USED TO DESIGNATE WHETHER THE FIT
C        WAS EXACT TO MACHINE PRECISION (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION FPLRS
C        THE FLOATING POINT LARGEST RELATIVE SPACING.
C     DOUBLE PRECISION RSS
C        THE RESIDUAL SUM OF SQUARES.
C     DOUBLE PRECISION RSSTST
C        THE VALUE FOR TESTING WHETHER THE RESIDUAL SUM OF SQUARES
C        IS ZERO (TO WITHIN MACHINE PRECISION).
C     DOUBLE PRECISION YSS
C        THE SUM OF SQUARES OF THE DEPENDENT VARIABLE Y.
C
      FPLRS = D1MACH(4)
C
C     TEST FOR EXACT FIT
C
      EXACT = .FALSE.
      RSSTST = RSS
      IF (YSS.GT.0.0D0) RSSTST = RSSTST / YSS
      RSSTST = SQRT(RSSTST)
      IF (RSSTST.LT.10.0D0*FPLRS) EXACT = .TRUE.
C
      RETURN
C
      END

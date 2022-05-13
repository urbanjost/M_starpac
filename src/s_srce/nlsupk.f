*NLSUPK
      SUBROUTINE NLSUPK(PARE, NPARE, PAR, MASK, NPAR)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE UNPACKS A VECTOR PARE INTO A VECTOR PAR, BY
C     PLACING SUCCEDING ELEMENTS OF PARE INTO ELEMENTS OF PAR
C     WHICH CORRESPOND TO ELEMENTS OF MASK WITH THE VALUE 1.
C     OTHER ELEMENTS OF MASK SHOULD BE 0.  THE NUMBER OF ELEMENTS
C     NPARE IN PARE SHOULD EQUAL THE NUMBER OF ELEMENTS OF
C     MASK WHICH ARE 1.
C
C     WRITTEN BY - JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   NPAR,NPARE
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),PARE(NPAR)
      INTEGER
     +   MASK(NPAR)
C
C  LOCAL SCALARS
      INTEGER
     +   I,JPK
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER JPK
C        AN INDEX VARIABLE.
C     INTEGER MASK(NPAR)
C        INPUT PARAMETER.  THE MASK GOVERNING THE PACKING OF PAR.
C        ELEMENTS OF MASK ARE 1 IF THE CORRESPONDING ELEMENT OF PAR
C        WAS ELIMINATED IN PARE, 0 IF IT WAS INCLUDED.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE OPTIMIZED.
C     REAL PARE(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS BEING OPTIMIZED,
C        NOT INCLUDING THOSE WHOSE VALUES ARE FIXED.
C
C     COMMENCE BODY OF ROUTINE
C
      JPK = 0
      DO 20 I=1,NPAR
         IF (MASK(I).NE.0) GO TO 20
         JPK = JPK + 1
         PAR(I) = PARE(JPK)
   20 CONTINUE
      RETURN
      END

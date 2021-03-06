*NLSPK
      SUBROUTINE NLSPK(PAR, MASK, NPAR, PPAR, NPPAR)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PACKS A VECTOR PAR INTO A VECTOR PPAR, BY
C     OMITTING FROM THE PACKED VERSION THOSE ELEMENTS OF THE
C     UNPACKED VERSION CORRESPONDING TO ELEMENTS OF MASK WHICH
C     HAVE THE VALUE 1.  OTHER ELEMENTS OF MASK SHOULD BE ZERO.
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
     +   NPAR,NPPAR
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),PPAR(NPPAR)
      INTEGER
     +   MASK(NPAR)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPPAR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL PAR(NPAR)
C        INPUT PARAMETER.  THE UNPACKED VECTOR.
C     INTEGER I
C        LOOP PARAMETER.
C     INTEGER IPPAR
C        CURRENT ELEMENT OF PPAR.  RANGES FROM 0 (ON INITIALIZATION)
C        TO NPPAR.
C     INTEGER MASK(NPAR)
C        INPUT PARAMETER.  THE MASK GOVERNING THE PACKING OF PAR.
C        ELEMENTS OF MASK ARE 1 IF THE CORRESPONDING ELEMENT OF PAR
C        IS TO BE ELIMINATED IN PPAR, 0 IF IT IS TO BE INCLUDED.
C     INTEGER NPAR
C        INPUT PARAMETER.  THE LENGTH OF PAR.
C     INTEGER NPPAR
C        INPUT PARAMETER.  THE LENGTH OF PPAR.
C     REAL PPAR(NPPAR)
C        OUTPUT PARAMETER.  THE PACKED VERSION OF PAR.  SEE INITIAL
C        DESCRIPTION.
C
C     COMMENCE BODY OF ROUTINE
C
      IPPAR = 0
      DO 10 I=1,NPAR
         IF (MASK(I).NE.0) GO TO 10
         IPPAR = IPPAR + 1
         PPAR(IPPAR) = PAR(I)
   10 CONTINUE
      RETURN
      END

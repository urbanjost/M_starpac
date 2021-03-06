*SETROW
      SUBROUTINE SETROW (NROW, XM, N, M, IXM, NROWU)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SELECTS THE ROW USED BY THE DERIVATIVE CHECKING
C     PROCEDURE.
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
      INTEGER
     +   IXM,M,N,NROW,NROWU
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   XM(IXM,M)
C
C  LOCAL SCALARS
      INTEGER
     +   I,J
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY XM.
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS OF DATA.
C     INTEGER NROW, NROWU
C        THE USER-SUPPLIED NUMBER OF THE ROW OF THE INDEPENDENT
C        VARIABLE ARRAY AT WHICH THE DERIVATIVE IS TO BE CHECKED,
C        AND THE NUMBER OF THE ROW ACTUALLY USED.
C     DOUBLE PRECISION XM(IXM,M)
C        THE INDEPENDENT VARIABLE MATRIX.
C
      NROWU = NROW
C
      IF ((NROWU.GE.1) .AND. (NROWU.LE.N)) RETURN
C
C     SELECT FIRST ROW OF INDEPENDENT VARIABLES WHICH CONTAINS NO ZEROS
C     IF THERE IS ONE, OTHERWISE FIRST ROW IS USED.
C
      DO 20 I = 1, N
         DO 10 J = 1, M
            IF (XM(I,J) .EQ. 0.0D0) GO TO 20
   10    CONTINUE
         NROWU = I
         RETURN
   20 CONTINUE
C
      NROWU = 1
C
      RETURN
      END

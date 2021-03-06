*SETIV
      SUBROUTINE SETIV(VECTOR, N, VALUE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SETS THE FIRST N ELEMENTS OF AN INTEGER VECTOR
C
C     WRITTEN BY  -  JOHN E. KOONTZ
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C        ADAPTED FROM SETRV, WRITTEN BY LINDA L. MITCHELL
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,VALUE
C
C  ARRAY ARGUMENTS
      INTEGER
     +   VECTOR(N)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        *
C     INTEGER N
C        NUMBER OF ELEMENTS TO SET
C     INTEGER VALUE
C        VALUE TO WHICH THE ELEMENTS ARE TO BE SET
C     INTEGER VECTOR(N)
C        VECTOR WHOSE FIRST N ELEMENTS ARE TO BE SET.
C
      DO 10 I=1,N
         VECTOR(I) = VALUE
   10 CONTINUE
C
      RETURN
C
      END

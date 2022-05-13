*GENI
      SUBROUTINE GENI(IVECT, N, IINIT, ISTP)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     PUT VALUES IINIT STEP ISTP THROUGH IINIT + (N - 1)*ISTP INTO
C     A VECTOR IVECT OF LENGTH N.  NO ERROR CHECKING IS DONE.
C
C     WRITTEN BY - JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING LAB/BOULDER
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IINIT,ISTP,N
C
C  ARRAY ARGUMENTS
      INTEGER
     +   IVECT(N)
C
C  LOCAL SCALARS
      INTEGER
     +   I,J
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        INITIALIZATION VALUE.
C     INTEGER IINIT, ISTP
C        INPUT PARAMETERS.  THE INITIAL VALUE AND THE INCREMENT USED
C        IN CREATING THE INITIALIZATION VALUES.
C     INTEGER IVECT(N)
C        OUTPUT PARAMETER.  THE VECTOR INTO WHICH TO PUT THE VALUES
C        IINIT, IINIT + ISTP, ..., IINIT + (N - 1)*ISTP.
C     INTEGER J
C        LOOP PARAMETER.
C     INTEGER N
C        INPUT PARAMETER.  THE LENGTH OF IVECT.
C
      I = IINIT
      DO 10 J=1,N
         IVECT(J) = I
         I = I + ISTP
   10 CONTINUE
      RETURN
      END

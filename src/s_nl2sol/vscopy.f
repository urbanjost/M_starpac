*VSCOPY
      SUBROUTINE VSCOPY(P, Y, S)
C
C  ***  SET P-VECTOR Y TO SCALAR S  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   S
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      REAL
     +   Y(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C
      DO 10 I = 1, P
 10      Y(I) = S
      RETURN
      END

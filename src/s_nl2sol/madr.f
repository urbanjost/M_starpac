*MADR
      SUBROUTINE MADR(N, P, X, NF, R, UIPARM, URPARM, UFPARM)
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NF,P
C
C  ARRAY ARGUMENTS
      REAL
     +   R(N),URPARM(*),X(P)
      INTEGER
     +   UIPARM(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL UFPARM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC COS,SIN
C
      R(1) = X(1)**2 + X(2)**2 + X(1)*X(2)
      R(2) = SIN(X(1))
      R(3) = COS(X(2))
      RETURN
      END

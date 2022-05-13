*PARZEN
      SUBROUTINE PARZEN (LAG, W, LW)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE COMPUTES AND STORES THE PARZEN LAG WINDOW
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LAG,LW
C
C  ARRAY ARGUMENTS
      REAL
     +   W(LW)
C
C  LOCAL SCALARS
      INTEGER
     +   K,L
C
C  INTRINSIC FUNCTIONS
      INTRINSIC REAL
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER K
C        AN INDEX VARIABLE
C     INTEGER L
C        THE VALUE LAG/2.
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LW
C        THE LENGTH OF THE VECTOR W.
C     REAL W(LW)
C        THE VECTOR OF LAG WINDOWS.
C
      L = LAG/2
      W(1) = 1.0E0
      IF (L.LE.0) GO TO 15
      DO 10 K = 1, L
         W(K+1) = REAL(K) / REAL(LAG)
         W(K+1) = 1.0E0 + 6.0E0 * W(K+1) * W(K+1) * (W(K+1) - 1.0E0)
   10 CONTINUE
C
   15 CONTINUE
      L = L + 1
      DO 20 K = L, LAG
         W(K+1) = 1.0E0 - REAL(K) / REAL(LAG)
         W(K+1) = 2.0E0 * W(K+1) * W(K+1) * W(K+1)
   20 CONTINUE
C
      RETURN
      END

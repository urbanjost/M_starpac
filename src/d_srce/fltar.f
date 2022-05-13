*FLTAR
      SUBROUTINE FLTAR (Y, N, IAR, PHI, YF, NYF)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE FILTERS THE INPUT SERIES Y USING THE IAR TERMS
C     OF THE AUTOREGRESSIVE FILTER PHI, COPYING THE FILTERED SERIES
C     INTO YF.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DEVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IAR,N,NYF
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PHI(*),Y(*),YF(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   TEMP
      INTEGER
     +   I,I1,J,K
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER IAR
C        THE NUMBER OF FILTER TERMS.
C     INTEGER I1, J, K
C        INDEXING VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES Y.
C     INTEGER NYF
C        THE NUMBER OF OBSERVATIONS IN THE FILTERED SERIES YF.
C     DOUBLE PRECISION PHI(IAR)
C        THE ARRAY IN WHICH THE FILTER COEFFICIENTS ARE STORED.
C     DOUBLE PRECISION TEMP
C        A TEMPORARY STORAGE LOCATION.
C     DOUBLE PRECISION Y(N)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES.
C     DOUBLE PRECISION YF(N)
C        THE VECTOR IN WHICH THE FILTERED SERIES IS RETURNED.
C
      DO 10 I = 1, N
         YF(I) = Y(I)
   10 CONTINUE
C
      NYF = N - IAR
C
      DO 30 I = 1, NYF
         K = I + IAR
         TEMP = YF(K)
         DO 20 J = 1, IAR
            K = K - 1
            TEMP = TEMP - PHI(J) * YF(K)
   20    CONTINUE
         YF(I) = TEMP
   30 CONTINUE
C
      I1 = NYF + 1
C
      DO 40 I = I1, N
         YF(I) = 0.0D0
   40 CONTINUE
      RETURN
      END

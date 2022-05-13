*DOTCM
      SUBROUTINE DOTCM (Y, YMEAN, YMISS, NY, X, XMEAN, XMISS,
     +   NX, DOTXY, NDOTXY)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE DOT PRODUCT OF TWO
C     SERIES WITH MISSING DATA, CENTERED ABOUT THEIR RESPECTIVE MEANS.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 21, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   DOTXY,XMEAN,XMISS,YMEAN,YMISS
      INTEGER
     +   NDOTXY,NX,NY
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(NX),Y(NY)
C
C  LOCAL SCALARS
      INTEGER
     +   I,M
C
C  EXTERNAL FUNCTIONS
      LOGICAL
     +   MVCHK
      EXTERNAL MVCHK
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DOTXY
C        THE DOT PRODUCT OF THE SERIES (Y(I) - YMEAN) AND
C        (X(I) - XMEAN).
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER M
C        THE SMALLER OF THE NUMBER OF OBSERVATIONS IN X AND Y
C     INTEGER NDOTXY
C        THE NUMBER OF OBSERVATIONS USED TO COMPUTE DOTXY
C     INTEGER NX, NY
C        THE NUMBER OF OBSERVATIONS IN SERIES X AND Y, RESPECTIVELY.
C     DOUBLE PRECISION X(NX)
C        THE VECTOR CONTAINING THE SECOND SERIES
C     DOUBLE PRECISION XMEAN
C        THE MEAN OF THE SECOND SERIES.
C     DOUBLE PRECISION XMISS
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IN THE SERIES IS MISSING.  IF X(I) = XMISS,
C        THE VALUE IS ASSUMED MISSING, OTHERWISE IT IS NOT.
C     DOUBLE PRECISION Y(NY)
C        THE VECTOR CONTAINING THE FIRST SERIES
C     DOUBLE PRECISION YMEAN
C        THE MEAN OF THE FIRST SERIES.
C     DOUBLE PRECISION YMISS
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IN THE SERIES IS MISSING.  IF Y(I) = YMISS,
C        THE VALUE IS ASSUMED MISSING, OTHERWISE IT IS NOT.
C
      NDOTXY = 0
      DOTXY = 0.0D0
      M = MIN(NY, NX)
      DO 10 I = 1, M
         IF (MVCHK(Y(I), YMISS) .OR. MVCHK(X(I), XMISS)) GO TO 10
         DOTXY = DOTXY + (Y(I) - YMEAN) * (X(I) - XMEAN)
         NDOTXY = NDOTXY + 1
   10 CONTINUE
      RETURN
      END

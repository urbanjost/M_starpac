*CCFMNF
      SUBROUTINE CCFMNF (Y1, Y2, N, NFFT, LAGMAX, NCC, CCOV11, CCOV22,
     +   CCOV12, CCOV21, LCCOV, RHOC, SDRHOC, NPRT, LYFFT, WORK, LWORK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN SUBROUTINE FOR COMPUTING CROSS CORRELATIONS AND
C     THEIR STANDARD ERRORS OF A TIME SERIES USING A FFT.
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
      INTEGER
     +   LAGMAX,LCCOV,LWORK,LYFFT,N,NCC,NFFT,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   CCOV11(LCCOV),CCOV12(LCCOV),CCOV21(LCCOV),CCOV22(LCCOV),
     +   RHOC(NCC),SDRHOC(NCC),WORK(LWORK),Y1(LYFFT),Y2(LYFFT)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FAC
      INTEGER
     +   I,I0,IM,IP
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CCFSD,CCVFF
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CCOV11(LCCOV)
C        THE ARRAY IN WHICH THE AUTOCOVARIANCE FUNCTION ESTIMATES
C        FOR THE FIRST SERIES IS STORED.
C     DOUBLE PRECISION CCOV12(LCCOV), CCOV21(LCCOV)
C        THE ARRAYS IN WHICH THE CROSS COVARIANCE FUNCTION
C        ESTIMATES FOR THE FIRST SERIES LAGGED BEHIND THE SECOND
C        AND VISA VERSA, ARE STORED.
C     DOUBLE PRECISION CCOV22(LCCOV)
C        THE ARRAY IN WHICH THE AUTOCOVARIANCE FUNCTION ESTIMATES
C        FOR THE SECOND SERIES IS STORED.
C     DOUBLE PRECISION FAC
C        THE INVERSE OF THE SQUARE ROOT OF THE PRODUCT OF THE
C        AUTOCOVARIANCES AT LAG ZERO.
C     INTEGER I
C        THE INDEXING VARIABLE FOR THE LAG VALUE.
C     INTEGER  IM, IP, I0
C        THE LOCATIONS IN THE CCF RELATED ARRAYS
C        OF THE LAG -I, I, AND 0, RESPECTIVELY.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LCCOV
C        THE DIMENSION OF THE COVARANCE ARRAYS.
C     INTEGER LWORK
C        THE DIMENSION OF THE WORK ARRAY.
C     INTEGER LYFFT
C        THE DIMENSION OF THE DATA ARRAYS.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NCC
C        THE NUMBER OF CROSS CORRELATIONS COMPUTED.
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO CONTROL COMPUTATIONS NEEDED
C        ONLY FOR PRINTED OUTPUT.
C     DOUBLE PRECISION RHOC(NCC)
C        THE ARRAY IN WHICH THE CROSS CORRELATIONS ARE STORED
C     DOUBLE PRECISION SDRHOC(NCC)
C        THE ARRAY CONTAINING THE STD. ERRORS OF THE CROSS CORRELATIONS.
C        ARE STORED.
C     DOUBLE PRECISION WORK(LWORK)
C        THE WORK ARRAY.
C     DOUBLE PRECISION Y1(LYFFT), Y2(LYFFT)
C        THE VECTORS CONTAINING THE OBSERVED SERIES
C
C     COMPUTE THE CROSS CORRELATIONS.
C
      CALL CCVFF (Y1, Y2, N, NFFT, LAGMAX, CCOV12, CCOV21, LCCOV, LYFFT,
     +   WORK, LWORK)
C
      IF (NPRT .EQ. 0 .OR. CCOV11(1)*CCOV22(1) .EQ. 0.0D0) RETURN
C
      FAC = 1.0D0 / SQRT(CCOV11(1) * CCOV22(1))
C
      I0 = LAGMAX + 1
      RHOC(I0) = CCOV12(1) * FAC
      DO 10 I = 1, LAGMAX
         IP = I0 + I
         RHOC(IP) = CCOV12(I+1) * FAC
C
         IM = I0 - I
         RHOC(IM) = CCOV21(I+1) * FAC
   10 CONTINUE
C
C     COMPUTE STANDARD ERROR OF THE CROSSCORRELATIONS.
C
      CALL CCFSD (CCOV11, CCOV22, SDRHOC, LAGMAX, NCC, N, LCCOV)
C
      RETURN
      END

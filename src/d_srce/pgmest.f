*PGMEST
      SUBROUTINE PGMEST (YFFT, NFFT, NF, CNST, PER, LPER)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE PERIODOGRAM ESTIMATES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   CNST
      INTEGER
     +   LPER,NF,NFFT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PER(LPER),YFFT(NFFT)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FAC
      INTEGER
     +   I,ISN,NFFT2
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FFT,REALTR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CNST
C        THE VARIANCE OF THE OBSERVED TIME SERIES TIMES THE NUMBER OF
C        OBSERVATIONS IN THE SERIES IF CALLED BY IPGM,
C        OR 1.0D0 IF CALLED BY PGM.
C     DOUBLE PRECISION FAC
C        A FACTOR USED FOR COMPUTATIONS OF THE INTEGRATED PERIODOGRAM.
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER ISN
C        A CODE USED FOR THE FFT.
C     INTEGER LPER
C        THE LENGTH OF THE PERIODOGRAM ARRAY.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE PERIODOGRAM IS
C        COMPUTED.
C     INTEGER NFFT
C        THE EFFECTIVE NUMBER OF OBSERVATIONS FOR THE FFT TRANSFORM.
C     INTEGER NFFT2
C        THE EFFECTIVE NUMBER OF COMPLEX OBSERVATIONS FOR THE FFT
C        TRANSFORM.
C     DOUBLE PRECISION PER(LPER)
C        THE PERIODOGRAM.
C     DOUBLE PRECISION YFFT(NFFT)
C        THE CENTERED SERIES.
C
C     COMPUTE THE FOURIER COEFFICIENTS
C
      NFFT2 = (NFFT-2) / 2
      ISN = 2
C
      CALL FFT (YFFT(1), YFFT(2), NFFT2, NFFT2, NFFT2, ISN)
      CALL REALTR (YFFT(1), YFFT(2), NFFT2, ISN)
C
      FAC = 0.5D0 / (CNST * (NFFT-2))
C
      NF = NFFT2 + 1
C
      DO 10 I = 1, NF
         PER(I) = (YFFT(2*I-1)*YFFT(2*I-1) + YFFT(2*I)*YFFT(2*I)) * FAC
   10 CONTINUE
C
      RETURN
      END

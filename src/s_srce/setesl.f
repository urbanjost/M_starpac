*SETESL
      SUBROUTINE SETESL(N, NDIV, NFFT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE SMALLEST VALUE OF NFFT WHICH
C     EQUALS OR EXCEEDS N + 2, SUCH THAT NFFT - 2 IS
C     1. DIVISIBLE BY NDIV,
C     2. HAS NO MORE THAN 11 PRIME FACTORS,
C     3. HAS NO PRIME FACTOR GREATER THAN 23, AND
C     4. THE PRODUCT OF THE SQUARE FREE PRIME FACTORS OF
C        (NFFT-2)/NDIV DO NOT EXCEED 210 IF NDIV = 2, AND
C                                    105 IF NDIV = 4.
C     THE VALUE OF NFFT THUS MEET THE REQUIREMENTS OF
C     THE EXTENDED LENGTH OF THE SERIES REQUIRED FOR ANY ROUTINE
C     USING THE SINGLETON FFT PROVIDING THE PROPER VALUE OF NDIV
C     IS CHOSEN.
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
      INTEGER
     +   N,NDIV,NFFT
C
C  LOCAL SCALARS
      INTEGER
     +   I,NPF,NSFP
C
C  LOCAL ARRAYS
      INTEGER
     +   IPF(50),IPFEXP(50)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FACTOR
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MOD
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C       AN INDEX VARIABLE.
C     INTEGER IPF(50), IPFEXP(50)
C        THE VECTORS OF PRIME FACTORS OF NFFT AND THEIR EXPONENTS,
C        RESPECTIVELY, WHERE THE LENGTH OF THESE VECTORS IS
C        SUFFICIENT TO ACCOMODATE THE PRIME FACTORS OF AN INTEGER
C        UP TO 2 ** 128 (APPROXIMATELY 10 ** 40).
C     INTEGER N
C        THE NUMBER UPON WHICH NFFT IS BASED.
C     INTEGER NDIV
C        A REQUIRED FACTOR OF NFFT - 2.
C     INTEGER NFFT
C        THE RETURNED VALUE WHICH MEETS THE ABOVE DESCRIPTION.
C     INTEGER NPF
C        THE NUMBER OF PRIME FACTORS IN NFFT.
C     INTEGER NSFP
C        THE PRODUCT OF THE NON SQUARE FACTORS.
C
      NFFT = N
      IF (NFFT.LE.0) RETURN
      IF (MOD(NFFT, NDIV) .NE. 0) NFFT = NFFT + NDIV - MOD(NFFT, NDIV)
      NFFT = NFFT - NDIV
   20 NFFT = NFFT + NDIV
      CALL FACTOR(NFFT/NDIV, NPF, IPF, IPFEXP)
      IF ((NPF.GE.11) .OR. (IPF(NPF).GT.23)) GO TO 20
      NSFP = 1
      IF (NDIV.EQ.4) NSFP = 2
      DO 30 I = 1, NPF
         IF (MOD(IPFEXP(I), 2).EQ.1) NSFP = NSFP * IPF(I)
   30 CONTINUE
      IF (NSFP .GE. 210) GO TO 20
C
      NFFT = NFFT + 2
C
      RETURN
C
      END

*OANOVA
      SUBROUTINE OANOVA(YSUM, RED, NPAR, RVAR, NNZW, TEMP, IPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     COMPUTE AND PRINT ANALYSIS OF VARIANCE
C
C     WRITTEN BY DAVID HOGBEN, SEL, NBS.   10/09/69.
C
C     THIS ROUTINE WAS ADAPTED FROM THE OMNITAB ROUTINE OANOVA
C     BY - -
C
C     JANET R. DONALDSON
C     STATISTICAL ENGINEERING DIVISION
C     NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 29, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   RVAR,YSUM
      INTEGER
     +   IPRT,NNZW,NPAR
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   RED(NPAR),TEMP(NPAR)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ASUM,CR,F1,F2,FPLM,PF1,PF2,RESMS,RESSS,SSU,V1F2,VR
      INTEGER
     +   I,K,NSUA
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   CDFF,D1MACH
      EXTERNAL CDFF,D1MACH
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ASUM
C        *
C     DOUBLE PRECISION CR
C        *
C     DOUBLE PRECISION FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     DOUBLE PRECISION F1
C        *
C     DOUBLE PRECISION F2
C        *
C     INTEGER I
C        AN INDEX.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER K
C        *
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NSUA
C        *
C     DOUBLE PRECISION PF1
C        *
C     DOUBLE PRECISION PF2
C        *
C     DOUBLE PRECISION RED(NPAR)
C        THE REDUCTION TO THE SUM OF SQUARES DUE TO EACH PARAMETER.
C     DOUBLE PRECISION RESMS
C        *
C     DOUBLE PRECISION RESSS
C        *
C     DOUBLE PRECISION RVAR
C        THE RESIDUAL VARIANCE.
C     DOUBLE PRECISION SSU
C        *
C     DOUBLE PRECISION TEMP(NPAR)
C        A WORK VECTOR.
C     DOUBLE PRECISION VR
C        *
C     DOUBLE PRECISION V1F2
C        *
C     DOUBLE PRECISION YSUM
C        THE SUM OF THE WEIGHTED DEPENDENT VARIABLES SQUARED.
C
C
      FPLM = D1MACH(2)
C
      RESMS = YSUM/NNZW
      NSUA = NNZW
      WRITE (IPRT,1000)
      ASUM = 0.0D0
      VR = NNZW-NPAR
      RESSS = VR*RVAR
      TEMP(NPAR) = RESSS
      IF (NPAR.EQ.1) GO TO 20
      DO 10 I=2,NPAR
         K = NPAR + 2 - I
         TEMP(K-1) = TEMP(K) + RED(K)
   10 CONTINUE
   20 V1F2 = NPAR+1
      SSU = NNZW
      DO 50 I=1,NPAR
         NSUA = NSUA - 1
         ASUM = ASUM + RED(I)
         SSU = SSU - 1.0D0
         CR = ASUM/I
         RESMS = 0.0D0
         IF (SSU.GT.0.0D0) RESMS = TEMP(I)/SSU
         V1F2 = V1F2 - 1.0D0
C
C     NEVER POOL
C
         IF (RVAR.GT.0.0D0) GO TO 30
         F1 = FPLM
         F2 = FPLM
         PF1 = 0.0D0
         PF2 = 0.0D0
         GO TO 40
   30    F1 = RED(I)/RVAR
         PF1 = 1.0D0 - CDFF(F1,1.0D0,VR)
C
C     TEST HIGHER SUB-HYPOTHESES
C
         F2 = (TEMP(I)+RED(I)-RESSS)/V1F2/RVAR
         PF2 = 1.0D0 - CDFF(F2,V1F2,VR)
   40    CONTINUE
         WRITE (IPRT,1010) I, RED(I), CR, I, RESMS, NSUA, F1, PF1, F2,
     +      PF2
   50 CONTINUE
      WRITE (IPRT,1020) RESSS, NSUA
      WRITE (IPRT,1030) YSUM, NNZW
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (////50X, 20HANALYSIS OF VARIANCE/24X, 16H-DEPENDENT ON OR,
     +   33HDER VARIABLES ARE ENTERED, UNLESS, 21H VECTORS ARE ORTHOGON,
     +   3HAL-//
     +   1X, 5H PAR , 4X, 14HSUM OF SQUARES, 63X,
     +   19H------ PAR=0 ------, 4X, 19H------ PARS=0 -----/
     +   1X, 5HINDEX, 4X, 14HRED DUE TO PAR, 7X, 10HCUM MS RED,
     +   6X, 9HDF(MSRED), 6X, 10HCUM RES MS, 6X, 7HDF(RMS), 5X,
     +   'F', 8X, 7HPROB(F), 7X, 'F', 8X, 7HPROB(F)/)
 1010 FORMAT (1X, I3, 6X, G16.9, 3X, G16.9, 1X, I6, 8X, G16.9, 1X, I5,
     +   4X, G12.6, F7.3, 4X, G12.6, F7.3)
 1020 FORMAT (/1X, 10HRESIDUAL  , 1X, G14.7, 20X, I6)
 1030 FORMAT (1X, 10HTOTAL     , 1X, G14.7, 20X, I6)
      END

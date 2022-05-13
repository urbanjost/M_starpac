*VCVOUT
      SUBROUTINE VCVOUT(NP, VCV, IVCV, EST)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE VARIANCE COVARIANCE MATRIX.
C     IF EST IS TRUE, THE COVARIANCES ARE LISTED ABOVE THE
C     DIAGONAL, THE VARIANCES ON THE DIAGONAL, AND THE CORRELATION
C     COEFFICIENTS BELOW THE DIAGONAL.
C     IF EST IS FALSE, THE STANDARD DEVIATIONS ARE LISTED ON THE
C     DIAGONAL, AND THE CORRELATION COEFFICIENTS ARE BELOW THE
C     DIAGONAL.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IVCV,NP
      LOGICAL
     +   EST
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   VCV(IVCV,NP)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   DEN,SVCVII,SVCVJJ
      INTEGER
     +   I,IPRT,J,K,MODE
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,MATPRT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DEN
C        DENOMINATOR OF (I, J) CORRELATION COEFFICIENT
C     LOGICAL EST
C        AN INDICATOR USED TO DESIGNATE WHETHER THE VCV TO BE PRINTED
C        IS OF THE ESTIMATED PARAMETERS (TRUE) OR NOT (FALSE).
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IVCV
C        THE EXACT FIRST DIMENSION OF THE MATRIX VCV.
C     INTEGER J
C        THE INDEX OF THE PARAMETER BEING EXAMINED.
C     INTEGER K
C        AN INDEX VARIABLE.
C     INTEGER MODE
C        IF MODE IS 1, PRINT FULL MATRIX.
C        IF MODE IS 2, PRINT LOWER TRIANGLE WITH SQUARE ROOTS OF
C                      OF THE DIAGONAL.
C     INTEGER NP
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     DOUBLE PRECISION SVCVII, SVCVJJ
C        SQUARE ROOTS OF VCV(I, I) AND VCV(J, J)
C     DOUBLE PRECISION VCV(IVCV,NP)
C        THE VARIANCE COVARIANCE MATRIX.
C
C     COMMENCE BODY OF ROUTINE
C
      CALL IPRINT(IPRT)
C
C     DETERMINE WHETHER TO ISSUE NEGATIVE VARIANCE WARNING
C
      MODE = 2
      DO 10 I=1,NP
         IF (VCV(I,I).GT.0.0D0) GO TO 10
         WRITE (IPRT,1000)
         IF (EST) WRITE (IPRT,1050)
         WRITE (IPRT,1010)
         MODE = 0
         GO TO 70
   10 CONTINUE
C
      IF (EST) GO TO 20
C
C     PRINT HEADING FOR CORRELATION ROUTINES
C
      WRITE (IPRT,1040)
      WRITE (IPRT,1030)
      MODE = 2
      GO TO 30
C
   20 CONTINUE
C
C     PRINT HEADING FOR ESTIMATION ROUTINES
C
      WRITE (IPRT,1050)
      WRITE (IPRT,1020)
      MODE = 1
C
   30 CONTINUE
C
C     COMPUTE THE CORRELATION COEFFICIENTS AND STORE IN THE BOTTOM HALF
C     OF THE VARIANCE COVARIANCE MATRIX
C
      IF (NP.EQ.1) GO TO 60
      DO 50 J=2,NP
         K = J - 1
         SVCVJJ = 0.0D0
         IF (VCV(J,J).GT.0.0D0) SVCVJJ = SQRT(VCV(J,J))
         DO 40 I=1,K
            SVCVII = 0.0D0
            IF (VCV(I,I).GT.0.0D0) SVCVII = SQRT(VCV(I,I))
            DEN = SVCVII*SVCVJJ
            IF (DEN.LE.0.0D0) VCV(J,I) = 0.0D0
            IF (DEN.GT.0.0D0) VCV(J,I) = VCV(J,I)/DEN
   40    CONTINUE
   50 CONTINUE
C
   60 CONTINUE
C
   70 CALL MATPRT(VCV, VCV, NP, IPRT, MODE, 1, IVCV)
C
C     RESTORE THE VCV MATRIX
C
      IF (NP.EQ.1) RETURN
      DO 90 J=2,NP
         K = J - 1
         DO 80 I=1,K
            VCV(J,I) = VCV(I,J)
   80    CONTINUE
   90 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/18H COVARIANCE MATRIX)
 1010 FORMAT (/39H     NONPOSITIVE VARIANCES ENCOUNTERED./10H     CORRE,
     +   39HLATION COEFFICIENTS CANNOT BE COMPUTED.)
 1020 FORMAT (4X, 36H- COVARIANCES ARE ABOVE THE DIAGONAL/4X, 7H- VARIA,
     +   24HNCES ARE ON THE DIAGONAL/4X, 27H- CORRELATION COEFFICIENTS ,
     +   22HARE BELOW THE DIAGONAL)
 1030 FORMAT (4X, 41H- STANDARD DEVIATIONS ARE ON THE DIAGONAL/4X,
     +   49H- CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL)
 1040 FORMAT (/19H CORRELATION MATRIX)
 1050 FORMAT (/45H VARIANCE-COVARIANCE AND CORRELATION MATRICES,
     +   28H OF THE ESTIMATED PARAMETERS/ 1X, 72('-')/)
      END

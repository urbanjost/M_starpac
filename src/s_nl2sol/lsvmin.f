*LSVMIN
      REAL FUNCTION LSVMIN(P, L, X, Y)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  ESTIMATE SMALLEST SING. VALUE OF PACKED LOWER TRIANG. MATRIX L
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      REAL
     +   L(1),X(P),Y(P)
C
C  LOCAL SCALARS
      REAL
     +   B,HALF,ONE,PSJ,R9973,SMINUS,SPLUS,T,XMINUS,XPLUS,ZERO
      INTEGER
     +   I,II,IX,J,J0,JI,JJ,JJJ,JM1,PPLUS1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   V2NORM
      EXTERNAL V2NORM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MOD
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER P
C     REAL L(1), X(P), Y(P)
C     DIMENSION L(P*(P+1)/2)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  PURPOSE  ***
C
C     THIS FUNCTION RETURNS A GOOD OVER-ESTIMATE OF THE SMALLEST
C     SINGULAR VALUE OF THE PACKED LOWER TRIANGULAR MATRIX L.
C
C  ***  PARAMETER DESCRIPTION  ***
C
C  P (IN)  = THE ORDER OF L.  L IS A  P X P  LOWER TRIANGULAR MATRIX.
C  L (IN)  = ARRAY HOLDING THE ELEMENTS OF  L  IN ROW ORDER, I.E.
C             L(1,1), L(2,1), L(2,2), L(3,1), L(3,2), L(3,3), ETC.
C  X (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN X IS A NORMALIZED
C             APPROXIMATE LEFT SINGULAR VECTOR CORRESPONDING TO THE
C             SMALLEST SINGULAR VALUE.  THIS APPROXIMATION MAY BE VERY
C             CRUDE.  IF LSVMIN RETURNS ZERO, THEN SOME COMPONENTS OF X
C             ARE ZERO AND THE REST RETAIN THEIR INPUT VALUES.
C  Y (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN Y = (L**-1)*X IS AN
C             UNNORMALIZED APPROXIMATE RIGHT SINGULAR VECTOR CORRESPOND-
C             ING TO THE SMALLEST SINGULAR VALUE.  THIS APPROXIMATION
C             MAY BE CRUDE.  IF LSVMIN RETURNS ZERO, THEN Y RETAINS ITS
C             INPUT VALUE.  THE CALLER MAY PASS THE SAME VECTOR FOR X
C             AND Y (NONSTANDARD FORTRAN USAGE), IN WHICH CASE Y OVER-
C             WRITES X (FOR NONZERO LSVMIN RETURNS).
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C     THERE ARE NO USAGE RESTRICTIONS.
C
C  ***  ALGORITHM NOTES  ***
C
C     THE ALGORITHM IS BASED ON (1), WITH THE ADDITIONAL PROVISION THAT
C     LSVMIN = 0 IS RETURNED IF THE SMALLEST DIAGONAL ELEMENT OF L
C     (IN MAGNITUDE) IS NOT MORE THAN THE UNIT ROUNDOFF TIMES THE
C     LARGEST.  THE ALGORITHM USES A RANDOM NUMBER GENERATOR PROPOSED
C     IN (4), WHICH PASSES THE SPECTRAL TEST WITH FLYING COLORS -- SEE
C     (2) AND (3).
C
C  ***  SUBROUTINES AND FUNCTIONS CALLED  ***
C
C        V2NORM - FUNCTION, RETURNS THE 2-NORM OF A VECTOR.
C
C  ***  REFERENCES  ***
C
C     (1) CLINE, A., MOLER, C., STEWART, G., AND WILKINSON, J.H.(1977),
C         AN ESTIMATE FOR THE CONDITION NUMBER OF A MATRIX, REPORT
C         TM-310, APPLIED MATH. DIV., ARGONNE NATIONAL LABORATORY.
C
C     (2) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL
C         RANDOM-NUMBER GENERATORS --  AN EMPIRICAL VIEW,
C         MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV.
C
C     (3) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2
C         (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS.
C
C     (4) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER
C         GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18,
C         PP. 586-593.
C
C  ***  HISTORY  ***
C
C     DESIGNED AND CODED BY DAVID M GAY (WINTER 1977/SUMMER 1978).
C
C  ***  GENERAL  ***
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, II, IX, J, JI, JJ, JJJ, JM1, J0, PPLUS1
C     REAL B, PSJ, SMINUS, SPLUS, T, XMINUS, XPLUS
C
C  ***  CONSTANTS  ***
C
C     REAL HALF, ONE, R9973, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL V2NORM
C     REAL V2NORM
C
      DATA IX/2/
      DATA HALF/0.5E0/, ONE/1.0E0/, R9973/9973.0E0/, ZERO/0.0E0/
C
C  ***  BODY  ***
C
C  ***  FIRST CHECK WHETHER TO RETURN LSVMIN = 0 AND INITIALIZE X  ***
C
      II = 0
      DO 10 I = 1, P
         X(I) = ZERO
         II = II + I
         IF (L(II) .EQ. ZERO) GO TO 300
 10      CONTINUE
      IF (MOD(IX, 9973) .EQ. 0) IX = 2
      PPLUS1 = P + 1
C
C  ***  SOLVE (L**T)*X = B, WHERE THE COMPONENTS OF B HAVE RANDOMLY
C  ***  CHOSEN MAGNITUDES IN (.5,1) WITH SIGNS CHOSEN TO MAKE X LARGE.
C
C     DO J = P TO 1 BY -1...
      DO 100 JJJ = 1, P
         J = PPLUS1 - JJJ
C       ***  DETERMINE X(J) IN THIS ITERATION. NOTE FOR I = 1,2,...,J
C       ***  THAT X(I) HOLDS THE CURRENT PARTIAL SUM FOR ROW I.
         IX = MOD(3432*IX, 9973)
         B = HALF*(ONE + IX/R9973)
         XPLUS = (B - X(J))
         XMINUS = (-B - X(J))
         SPLUS = ABS(XPLUS)
         SMINUS = ABS(XMINUS)
         JM1 = J - 1
         J0 = J*JM1/2
         JJ = J0 + J
         XPLUS = XPLUS/L(JJ)
         XMINUS = XMINUS/L(JJ)
         IF (JM1 .EQ. 0) GO TO 30
         DO 20 I = 1, JM1
              JI = J0 + I
              SPLUS = SPLUS + ABS(X(I) + L(JI)*XPLUS)
              SMINUS = SMINUS + ABS(X(I) + L(JI)*XMINUS)
 20           CONTINUE
 30      IF (SMINUS .GT. SPLUS) XPLUS = XMINUS
         X(J) = XPLUS
C       ***  UPDATE PARTIAL SUMS  ***
         IF (JM1 .EQ. 0) GO TO 100
         DO 40 I = 1, JM1
              JI = J0 + I
              X(I) = X(I) + L(JI)*XPLUS
 40           CONTINUE
 100     CONTINUE
C
C  ***  NORMALIZE X  ***
C
      T = ONE/V2NORM(P, X)
      DO 110 I = 1, P
 110     X(I) = T*X(I)
C
C  ***  SOLVE L*Y = X AND RETURN SVMIN = 1/TWONORM(Y)  ***
C
      DO 200 J = 1, P
         PSJ = ZERO
         JM1 = J - 1
         J0 = J*JM1/2
         IF (JM1 .EQ. 0) GO TO 130
         DO 120 I = 1, JM1
              JI = J0 + I
              PSJ = PSJ + L(JI)*Y(I)
 120          CONTINUE
 130     JJ = J0 + J
         Y(J) = (X(J) - PSJ)/L(JJ)
 200     CONTINUE
C
      LSVMIN = ONE/V2NORM(P, Y)
      GO TO 999
C
 300  LSVMIN = ZERO
 999  RETURN
C  ***  LAST CARD OF LSVMIN FOLLOWS  ***
      END

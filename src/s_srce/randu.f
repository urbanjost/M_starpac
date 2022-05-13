*RANDU
      REAL FUNCTION RANDU(JD)
C***BEGIN PROLOGUE  RANDU  (ORIGINALLY UNI)
C***DATE WRITTEN   810915
C***REVISION DATE  900315
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    BLUE, JAMES
C             KAHANER, DAVID
C             APPLIED AND COMPUTATIONAL MATHEMATICS DIVISTION, NIST
C
C             MARSAGLIA, GEORGE
C             COMPUTER SCIENCE DEPT., WASH STATE UNIV
C
C             MODIFIED BY -
C             DONALDSON, JANET
C             APPLIED AND COMPUTATIONAL MATHEMATICS DIVISTION, NIST
C
C***PURPOSE  THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON
C             (0,1] AND CAN BE USED ON ANY COMPUTER WITH WHICH ALLOWS
C             INTEGERS AT LEAST AS LARGE AS 32767.
C***DESCRIPTION
C
C       THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON THE
C       INTERVAL (0,1].  IT CAN BE USED WITH ANY COMPUTER WHICH ALLOWS
C       INTEGERS AT LEAST AS LARGE AS 32767.
C
C
C   USE
C       FIRST TIME....
C                   Z = RANDU(JD)
C                     HERE JD IS ANY  N O N - Z E R O  INTEGER.
C                     THIS CAUSES INITIALIZATION OF THE PROGRAM
C                     AND THE FIRST RANDOM NUMBER TO BE RETURNED AS Z.
C       SUBSEQUENT TIMES...
C                   Z = RANDU(0)
C                     CAUSES THE NEXT RANDOM NUMBER TO BE RETURNED AS Z.
C
C
C===================================================================
C   NOTE: USERS WHO WISH TO TRANSPORT THIS PROGRAM FROM ONE COMPUTER
C         TO ANOTHER SHOULD READ THE FOLLOWING INFORMATION:
C
C   MACHINE DEPENDENCIES...
C      MDIG = A LOWER BOUND ON THE NUMBER OF BINARY DIGITS AVAILABLE
C              FOR REPRESENTING INTEGERS, INCLUDING THE SIGN BIT.
C              THIS VALUE MUST BE AT LEAST 16, BUT MAY BE INCREASED
C              IN LINE WITH REMARK A BELOW.
C
C   REMARKS...
C     A. THIS PROGRAM CAN BE USED IN TWO WAYS:
C        (1) TO OBTAIN REPEATABLE RESULTS ON DIFFERENT COMPUTERS,
C            SET 'MDIG' TO THE SMALLEST OF ITS VALUES ON EACH, OR,
C        (2) TO ALLOW THE LONGEST SEQUENCE OF RANDOM NUMBERS TO BE
C            GENERATED WITHOUT CYCLING (REPEATING) SET 'MDIG' TO THE
C            LARGEST POSSIBLE VALUE.
C     B. THE SEQUENCE OF NUMBERS GENERATED DEPENDS ON THE INITIAL
C          INPUT 'JD' AS WELL AS THE VALUE OF 'MDIG'.
C          IF MDIG=16 ONE SHOULD FIND THAT
C            THE FIRST EVALUATION
C              Z=RANDU(305) GIVES Z=.027832881...
C            THE SECOND EVALUATION
C              Z=RANDU(0) GIVES   Z=.56102176...
C            THE THIRD EVALUATION
C              Z=RANDU(0) GIVES   Z=.41456343...
C            THE THOUSANDTH EVALUATION
C              Z=RANDU(0) GIVES   Z=.19797357...
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  I1MACH,XERROR
C***END PROLOGUE  RANDU
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   JD
C
C  LOCAL SCALARS
      REAL
     +   ONE,ZERO
      INTEGER
     +   I,J,J0,J1,JSEED,K,K0,K1,M1,M2,MDIG
C
C  LOCAL ARRAYS
      INTEGER
     +   M(17)
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   I1MACH
      EXTERNAL I1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL XERROR
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MIN,MOD,REAL
C
C  SAVE STATEMENT
      SAVE I,J,M,M1,M2
C
C
      DATA M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),M(9),M(10),M(11),
     +     M(12),M(13),M(14),M(15),M(16),M(17)/30788,23052,2053,19346,
     +     10646,19427,23975,19049,10949,19693,29746,26748,2796,23890,
     +     29168,31924,16499/
      DATA M1,M2,I,J/32767,256,5,17/
      DATA ZERO,ONE /0.0E0,1.0E0/
C
C***FIRST EXECUTABLE STATEMENT  RANDU
      IF (JD.NE.0) THEN
C  FILL
          MDIG = I1MACH(8) + 1
C
C  MODIFICATION SO SAME NUMBERS WILL BE GENERATED ON ALL MACHINES
C  WITH I1MACH(8) AT LEAST 31
C
          MDIG = MIN(MDIG,32)
C
C  BE SURE THAT MDIG AT LEAST 16...
          IF (MDIG.LT.16) CALL XERROR('RANDU--MDIG LESS THAN 16',22,1,2)
          M1 = 2** (MDIG-2) + (2** (MDIG-2)-1)
          M2 = 2** (MDIG/2)
          JSEED = MIN(ABS(JD),M1)
          IF (MOD(JSEED,2).EQ.0) JSEED = JSEED - 1
          K0 = MOD(9069,M2)
          K1 = 9069/M2
          J0 = MOD(JSEED,M2)
          J1 = JSEED/M2
          DO 10 I = 1,17
              JSEED = J0*K0
              J1 = MOD(JSEED/M2+J0*K1+J1*K0,M2/2)
              J0 = MOD(JSEED,M2)
              M(I) = J0 + M2*J1
   10     CONTINUE
          I = 5
          J = 17
      END IF
C  BEGIN MAIN LOOP HERE
      K = M(I) - M(J)
      IF (K.LT.0) K = K + M1
      M(J) = K
      I = I - 1
      IF (I.EQ.0) I = 17
      J = J - 1
      IF (J.EQ.0) J = 17
      RANDU = REAL(K)/REAL(M1)
C
C  MODIFICATION SO RANDOM NUMBERS IN (0,1] RATHER THAN [0,1)
C
      IF (RANDU.EQ.ZERO) RANDU = ONE
      END

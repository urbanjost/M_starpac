*GAMIT
      REAL FUNCTION GAMIT (A, X)
C JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
C
C EVALUATE TRICOMI-S INCOMPLETE GAMMA FUNCTION DEFINED BY
C
C GAMIT = X**(-A)/GAMMA(A) * INTEGRAL T = 0 TO X OF EXP(-T) * T**(A-1.)
C
C AND ANALYTIC CONTINUATION FOR A .LE. 0.0.  GAMMA(X) IS THE COMPLETE
C GAMMA FUNCTION OF X.  GAMIT IS EVALUATED FOR ARBITRARY REAL VALUES OF
C A AND FOR NON-NEGATIVE VALUES OF X (EVEN THOUGH GAMIT IS DEFINED FOR
C X .LT. 0.0), EXCEPT THAT FOR X = 0 AND A .LE. 0.0, GAMIT IS INFINITE,
C A FATAL ERROR.
C
C      A SLIGHT DETERIORATION OF 2 OR 3 DIGITS ACCURACY WILL OCCUR WHEN
C GAMIT IS VERY LARGE OR VERY SMALL IN ABSOLUTE VALUE, BECAUSE LOG-
C ARITHMIC VARIABLES ARE USED.  ALSO, IF THE PARAMETER A IS VERY CLOSE
C TO A NEGATIVE INTEGER (BUT NOT A NEGATIVE INTEGER), THERE IS A LOSS
C OF ACCURACY, WHICH IS REPORTED IF THE RESULT IS LESS THAN HALF
C MACHINE PRECISION.
C
C REF. -- W. GAUTSCHI, AN EVALUATION PROCEDURE FOR INCOMPLETE GAMMA
C FUNCTIONS, ACM TRANS. MATH. SOFTWARE.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL A,X
C
C  LOCAL SCALARS
      REAL AEPS,AINTA,ALGAP1,ALNEPS,ALNG,ALX,BOT,H,SGA,SGNGAM,SQEPS,T
C
C  EXTERNAL FUNCTIONS
      REAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
      EXTERNAL ALNGAM,GAMR,R1MACH,R9GMIT,R9LGIC,R9LGIT
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ALGAMS,XERCLR,XERROR
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP,LOG,SIGN,SQRT
C
      DATA ALNEPS, SQEPS, BOT / 3*0.0 /
C
      IF (ALNEPS.NE.0.0) GO TO 10
      ALNEPS = -LOG(R1MACH(3))
      SQEPS = SQRT(R1MACH(4))
      BOT = LOG(R1MACH(1))
C
 10   IF (X.LT.0.0) CALL XERROR ('GAMIT   X IS NEGATIVE', 21, 2, 2)
C
      IF (X.NE.0.0) ALX = LOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      AINTA = AINT (A+0.5*SGA)
      AEPS = A - AINTA
C
      IF (X.GT.0.0) GO TO 20
      GAMIT = 0.0
      IF (AINTA.GT.0.0 .OR. AEPS.NE.0.0) GAMIT = GAMR(A+1.0)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 40
      IF (A.GE.(-0.5) .OR. AEPS.NE.0.0) CALL ALGAMS (A+1.0, ALGAP1,
     1  SGNGAM)
      GAMIT = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
C
 40   IF (A.LT.X) GO TO 50
      T = R9LGIT (A, X, ALNGAM(A+1.0))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = EXP(T)
      RETURN
C
 50   ALNG = R9LGIC (A, X, ALX)
C
C EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
C
      H = 1.0
      IF (AEPS.EQ.0.0 .AND. AINTA.LE.0.0) GO TO 60
      CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      T = LOG(ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGA*SGNGAM*EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 60
      CALL XERCLR
      CALL XERROR ('GAMIT   RESULT LT HALF PRECISION', 32, 1, 1)
C
 60   T = -A*ALX + LOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = SIGN (EXP(T), H)
      RETURN
C
 70   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      GAMIT = -SGA*SGNGAM*EXP(T)
      RETURN
C
      END

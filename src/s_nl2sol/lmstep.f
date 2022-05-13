*LMSTEP
      SUBROUTINE LMSTEP(D, G, IERR, IPIVOT, KA, P, QTR, R, STEP, V, W)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  COMPUTE LEVENBERG-MARQUARDT STEP USING MORE-HEBDEN TECHNIQUE  **
C  ***  NL2SOL VERSION 2.2.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IERR,KA,P
C
C  ARRAY ARGUMENTS
      REAL
     +   D(P),G(P),QTR(P),R(1),STEP(P),V(21),W(1)
      INTEGER
     +   IPIVOT(P)
C
C  LOCAL SCALARS
      REAL
     +   A,ADI,ALPHAK,B,D1,D2,DFAC,DFACSQ,DST,DTOL,EIGHT,HALF,LK,
     +   NEGONE,OLDPHI,ONE,P001,PHI,PHIMAX,PHIMIN,PSIFAC,RAD,SI,SJ,
     +   SQRTAK,T,THREE,TTOL,TWOPSI,UK,WL,ZERO
      INTEGER
     +   DGNORM,DST0,DSTNRM,DSTSAV,EPSLON,GTSTEP,I,I1,IP1,J1,K,
     +   KALIM,L,LK0,NREDUC,PHIPIN,PHMNFC,PHMXFC,PP1O2,PREDUC,RAD0,
     +   RADIUS,RES,RES0,RMAT,RMAT0,STPPAR,UK0
C
C  EXTERNAL FUNCTIONS
      REAL
     +   DOTPRD,V2NORM
      EXTERNAL DOTPRD,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL LITVMU,LIVMUL,VCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX,MIN,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IERR, KA, P
C     INTEGER IPIVOT(P)
C     REAL D(P), G(P), QTR(P), R(1), STEP(P), V(21), W(1)
C     DIMENSION W(P*(P+5)/2 + 4)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  PURPOSE  ***
C
C        GIVEN THE R MATRIX FROM THE QR DECOMPOSITION OF A JACOBIAN
C     MATRIX, J, AS WELL AS Q-TRANSPOSE TIMES THE CORRESPONDING
C     RESIDUAL VECTOR, RESID, THIS SUBROUTINE COMPUTES A LEVENBERG-
C     MARQUARDT STEP OF APPROXIMATE LENGTH V(RADIUS) BY THE MORE-
C     TECHNIQUE.
C
C  ***  PARAMETER DESCRIPTION  ***
C
C      D (IN)  = THE SCALE VECTOR.
C      G (IN)  = THE GRADIENT VECTOR (J**T)*R.
C   IERR (I/O) = RETURN CODE FROM QRFACT OR QRFGS -- 0 MEANS R HAS
C             FULL RANK.
C IPIVOT (I/O) = PERMUTATION ARRAY FROM QRFACT OR QRFGS, WHICH COMPUTE
C             QR DECOMPOSITIONS WITH COLUMN PIVOTING.
C     KA (I/O).  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST CALL ON
C             LMSTEP FOR THE CURRENT R AND QTR.  ON OUTPUT KA CON-
C             TAINS THE NUMBER OF HEBDEN ITERATIONS NEEDED TO DETERMINE
C             STEP.  KA = 0 MEANS A GAUSS-NEWTON STEP.
C      P (IN)  = NUMBER OF PARAMETERS.
C    QTR (IN)  = (Q**T)*RESID = Q-TRANSPOSE TIMES THE RESIDUAL VECTOR.
C      R (IN)  = THE R MATRIX, STORED COMPACTLY BY COLUMNS.
C   STEP (OUT) = THE LEVENBERG-MARQUARDT STEP COMPUTED.
C      V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW.
C      W (I/O) = WORKSPACE OF LENGTH P*(P+5)/2 + 4.
C
C  ***  ENTRIES IN V  ***
C
C V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G.
C V(DSTNRM) (I/O) = 2-NORM OF D*STEP.
C V(DST0)   (I/O) = 2-NORM OF GAUSS-NEWTON STEP (FOR NONSING. J).
C V(EPSLON) (IN) = MAX. REL. ERROR ALLOWED IN TWONORM(R)**2 MINUS
C             TWONORM(R - J*STEP)**2.  (SEE ALGORITHM NOTES BELOW.)
C V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP.
C V(NREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
C             FOR A GAUSS-NEWTON STEP.
C V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP
C             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE
C             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS).
C V(PHMXFC) (IN)  (SEE V(PHMNFC).)
C V(PREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
C             BY THE STEP RETURNED.
C V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION.
C V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL.
C V(STPPAR) (I/O) = MARQUARDT PARAMETER (OR ITS NEGATIVE IF THE SPECIAL
C             CASE MENTIONED BELOW IN THE ALGORITHM NOTES OCCURS).
C
C NOTE -- SEE DATA STATEMENT BELOW FOR VALUES OF ABOVE SUBSCRIPTS.
C
C  ***  USAGE NOTES  ***
C
C     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF
C     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT
C     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS
C     WHY MANY PARAMETERS ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE
C     WITH KA = -1), THE CALLER NEED ONLY HAVE INITIALIZED D, G, KA, P,
C     QTR, R, V(EPSLON), V(PHMNFC), V(PHMXFC), V(RADIUS), AND V(RAD0).
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST-
C     SQUARES) PACKAGE (REF. 1).
C
C  ***  ALGORITHM NOTES  ***
C
C     THIS CODE IMPLEMENTS THE STEP COMPUTATION SCHEME DESCRIBED IN
C     REFS. 2 AND 4.  FAST GIVENS TRANSFORMATIONS (SEE REF. 3, PP. 60-
C     62) ARE USED TO COMPUTE STEP WITH A NONZERO MARQUARDT PARAMETER.
C        A SPECIAL CASE OCCURS IF J IS (NEARLY) SINGULAR AND V(RADIUS)
C     IS SUFFICIENTLY LARGE.  IN THIS CASE THE STEP RETURNED IS SUCH
C     THAT  TWONORM(R)**2 - TWONORM(R - J*STEP)**2  DIFFERS FROM ITS
C     OPTIMAL VALUE BY LESS THAN V(EPSLON) TIMES THIS OPTIMAL VALUE,
C     WHERE J AND R DENOTE THE ORIGINAL JACOBIAN AND RESIDUAL.  (SEE
C     REF. 2 FOR MORE DETAILS.)
C
C  ***  FUNCTIONS AND SUBROUTINES CALLED  ***
C
C DOTPRD - RETURNS INNER PRODUCT OF TWO VECTORS.
C LITVMU - APPLY INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
C LIVMUL - APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
C VCOPY  - COPIES ONE VECTOR TO ANOTHER.
C V2NORM - RETURNS 2-NORM OF A VECTOR.
C
C  ***  REFERENCES  ***
C
C 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1980), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, (SUBMITTED TO ACM
C             TRANS. MATH. SOFTWARE).
C 2.  GAY, D.M. (1979), COMPUTING OPTIMAL ELLIPTICALLY CONSTRAINED
C             STEPS, MRC TECH. SUMMARY REPORT NO. 2013, MATH RESEARCH
C             CENTER, UNIV. OF WISCONSIN-MADISON.
C 3.  LAWSON, C.L., AND HANSON, R.J. (1974), SOLVING LEAST SQUARES
C             PROBLEMS, PRENTICE-HALL, ENGLEWOOD CLIFFS, N.J.
C 4.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN-
C             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES
C             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER-
C             VERLAG, BERLIN AND NEW YORK.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY.
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER DSTSAV, I, IP1, I1, J1, K, KALIM, L, LK0, PHIPIN,
C    1        PP1O2, RES, RES0, RMAT, RMAT0, UK0
C     REAL A, ADI, ALPHAK, B, DFACSQ, DST, DTOL, D1, D2,
C    1                 LK, OLDPHI, PHI, PHIMAX, PHIMIN, PSIFAC, RAD,
C    2                 SI, SJ, SQRTAK, T, TWOPSI, UK, WL
C
C     ***  CONSTANTS  ***
C     REAL DFAC, EIGHT, HALF, NEGONE, ONE, P001, THREE,
C    1                 TTOL, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DOTPRD, LITVMU, LIVMUL, VCOPY, V2NORM
C     REAL DOTPRD, V2NORM
C
C  ***  SUBSCRIPTS FOR V  ***
C
C     INTEGER DGNORM, DSTNRM, DST0, EPSLON, GTSTEP, NREDUC, PHMNFC,
C    1        PHMXFC, PREDUC, RADIUS, RAD0, STPPAR
      DATA DGNORM/1/, DSTNRM/2/, DST0/3/, EPSLON/19/,
     +     GTSTEP/4/, NREDUC/6/, PHMNFC/20/,
     +     PHMXFC/21/, PREDUC/7/, RADIUS/8/,
     +     RAD0/9/, STPPAR/5/
C
      DATA DFAC/256.0E0/, EIGHT/8.0E0/, HALF/0.5E0/, NEGONE/-1.0E0/,
     +     ONE/1.0E0/, P001/1.0E-3/, THREE/3.0E0/, TTOL/2.5E0/,
     +     ZERO/0.0E0/
C
C  ***  BODY  ***
C
C     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK AND UK,
C     ***  THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR NONSING. J)
C     ***  AND THE VALUE RETURNED AS V(DSTNRM) ARE STORED AT W(LK0),
C     ***  W(UK0), W(PHIPIN), AND W(DSTSAV) RESPECTIVELY.
      ALPHAK = 0.0E0
      PSIFAC = 0.0E0
      LK0 = P + 1
      PHIPIN = LK0 + 1
      UK0 = PHIPIN + 1
      DSTSAV = UK0 + 1
      RMAT0 = DSTSAV
C     ***  A COPY OF THE R-MATRIX FROM THE QR DECOMPOSITION OF J IS
C     ***  STORED IN W STARTING AT W(RMAT), AND A COPY OF THE RESIDUAL
C     ***  VECTOR IS STORED IN W STARTING AT W(RES).  THE LOOPS BELOW
C     ***  THAT UPDATE THE QR DECOMP. FOR A NONZERO MARQUARDT PARAMETER
C     ***  WORK ON THESE COPIES.
      RMAT = RMAT0 + 1
      PP1O2 = P * (P + 1) / 2
      RES0 = PP1O2 + RMAT0
      RES = RES0 + 1
      RAD = V(RADIUS)
      IF (RAD .GT. ZERO)
     +   PSIFAC = V(EPSLON)/((EIGHT*(V(PHMNFC) + ONE) + THREE) * RAD**2)
      PHIMAX = V(PHMXFC) * RAD
      PHIMIN = V(PHMNFC) * RAD
C     ***  DTOL, DFAC, AND DFACSQ ARE USED IN RESCALING THE FAST GIVENS
C     ***  REPRESENTATION OF THE UPDATED QR DECOMPOSITION.
      DTOL = ONE/DFAC
      DFACSQ = DFAC*DFAC
C     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF
C     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT.
      OLDPHI = ZERO
      LK = ZERO
      UK = ZERO
      KALIM = KA + 12
C
C  ***  START OR RESTART, DEPENDING ON KA  ***
C
      IF (KA) 10, 20, 370
C
C  ***  FRESH START -- COMPUTE V(NREDUC)  ***
C
 10   KA = 0
      KALIM = 12
      K = P
      IF (IERR .NE. 0) K = ABS(IERR) - 1
      V(NREDUC) = HALF*DOTPRD(K, QTR, QTR)
C
C  ***  SET UP TO TRY INITIAL GAUSS-NEWTON STEP  ***
C
 20   V(DST0) = NEGONE
      IF (IERR .NE. 0) GO TO 90
C
C  ***  COMPUTE GAUSS-NEWTON STEP  ***
C
C     ***  NOTE -- THE R-MATRIX IS STORED COMPACTLY BY COLUMNS IN
C     ***  R(1), R(2), R(3), ...  IT IS THE TRANSPOSE OF A
C     ***  LOWER TRIANGULAR MATRIX STORED COMPACTLY BY ROWS, AND WE
C     ***  TREAT IT AS SUCH WHEN USING LITVMU AND LIVMUL.
      CALL LITVMU(P, W, R, QTR)
C     ***  TEMPORARILY STORE PERMUTED -D*STEP IN STEP.
      DO 60 I = 1, P
         J1 = IPIVOT(I)
         STEP(I) = D(J1)*W(I)
 60      CONTINUE
      DST = V2NORM(P, STEP)
      V(DST0) = DST
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX) GO TO 410
C     ***  IF THIS IS A RESTART, GO TO 110  ***
      IF (KA .GT. 0) GO TO 110
C
C  ***  GAUSS-NEWTON STEP WAS UNACCEPTABLE.  COMPUTE L0  ***
C
      DO 70 I = 1, P
         J1 = IPIVOT(I)
         STEP(I) = D(J1)*(STEP(I)/DST)
 70      CONTINUE
      CALL LIVMUL(P, STEP, R, STEP)
      T = ONE / V2NORM(P, STEP)
      W(PHIPIN) = (T/DST)*T
      LK = PHI*W(PHIPIN)
C
C  ***  COMPUTE U0  ***
C
 90   DO 100 I = 1, P
 100     W(I) = G(I)/D(I)
      V(DGNORM) = V2NORM(P, W)
      UK = V(DGNORM)/RAD
      IF (UK .LE. ZERO) GO TO 390
C
C     ***  ALPHAK WILL BE USED AS THE CURRENT MARQUARDT PARAMETER.  WE
C     ***  USE MORE*S SCHEME FOR INITIALIZING IT.
      ALPHAK = ABS(V(STPPAR)) * V(RAD0)/RAD
C
C
C  ***  TOP OF LOOP -- INCREMENT KA, COPY R TO RMAT, QTR TO RES  ***
C
 110  KA = KA + 1
      CALL VCOPY(PP1O2, W(RMAT), R)
      CALL VCOPY(P, W(RES), QTR)
C
C  ***  SAFEGUARD ALPHAK AND INITIALIZE FAST GIVENS SCALE VECTOR.  ***
C
      IF (ALPHAK .LE. ZERO .OR. ALPHAK .LT. LK .OR. ALPHAK .GE. UK)
     +             ALPHAK = UK * MAX(P001, SQRT(LK/UK))
      SQRTAK = SQRT(ALPHAK)
      DO 120 I = 1, P
 120     W(I) = ONE
C
C  ***  ADD ALPHAK*D AND UPDATE QR DECOMP. USING FAST GIVENS TRANS.  ***
C
      DO 270 I = 1, P
C        ***  GENERATE, APPLY 1ST GIVENS TRANS. FOR ROW I OF ALPHAK*D.
C        ***  (USE STEP TO STORE TEMPORARY ROW)  ***
         L = I*(I+1)/2 + RMAT0
         WL = W(L)
         D2 = ONE
         D1 = W(I)
         J1 = IPIVOT(I)
         ADI = SQRTAK*D(J1)
         IF (ADI .GE. ABS(WL)) GO TO 150
 130     A = ADI/WL
         B = D2*A/D1
         T = A*B + ONE
         IF (T .GT. TTOL) GO TO 150
         W(I) = D1/T
         D2 = D2/T
         W(L) = T*WL
         A = -A
         DO 140 J1 = I, P
              L = L + J1
              STEP(J1) = A*W(L)
 140          CONTINUE
         GO TO 170
C
 150     B = WL/ADI
         A = D1*B/D2
         T = A*B + ONE
         IF (T .GT. TTOL) GO TO 130
         W(I) = D2/T
         D2 = D1/T
         W(L) = T*ADI
         DO 160 J1 = I, P
              L = L + J1
              WL = W(L)
              STEP(J1) = -WL
              W(L) = A*WL
 160          CONTINUE
C
 170     IF (I .EQ. P) GO TO 280
C
C        ***  NOW USE GIVENS TRANS. TO ZERO ELEMENTS OF TEMP. ROW  ***
C
         IP1 = I + 1
         DO 260 I1 = IP1, P
              L = I1*(I1+1)/2 + RMAT0
              WL = W(L)
              SI = STEP(I1-1)
              D1 = W(I1)
C
C             ***  RESCALE ROW I1 IF NECESSARY  ***
C
              IF (D1 .GE. DTOL) GO TO 190
                   D1 = D1*DFACSQ
                   WL = WL/DFAC
                   K = L
                   DO 180 J1 = I1, P
                        K = K + J1
                        W(K) = W(K)/DFAC
 180                    CONTINUE
C
C             ***  USE GIVENS TRANS. TO ZERO NEXT ELEMENT OF TEMP. ROW
C
 190          IF (ABS(SI) .GT. ABS(WL)) GO TO 220
              IF (SI .EQ. ZERO) GO TO 260
 200          A = SI/WL
              B = D2*A/D1
              T = A*B + ONE
              IF (T .GT. TTOL) GO TO 220
              W(L) = T*WL
              W(I1) = D1/T
              D2 = D2/T
              DO 210 J1 = I1, P
                   L = L + J1
                   WL = W(L)
                   SJ = STEP(J1)
                   W(L) = WL + B*SJ
                   STEP(J1) = SJ - A*WL
 210               CONTINUE
              GO TO 240
C
 220          B = WL/SI
              A = D1*B/D2
              T = A*B + ONE
              IF (T .GT. TTOL) GO TO 200
              W(I1) = D2/T
              D2 = D1/T
              W(L) = T*SI
              DO 230 J1 = I1, P
                   L = L + J1
                   WL = W(L)
                   SJ = STEP(J1)
                   W(L) = A*WL + SJ
                   STEP(J1) = B*SJ - WL
 230               CONTINUE
C
C             ***  RESCALE TEMP. ROW IF NECESSARY  ***
C
 240          IF (D2 .GE. DTOL) GO TO 260
                   D2 = D2*DFACSQ
                   DO 250 K = I1, P
 250                    STEP(K) = STEP(K)/DFAC
 260          CONTINUE
 270     CONTINUE
C
C  ***  COMPUTE STEP  ***
C
 280  CALL LITVMU(P, W(RES), W(RMAT), W(RES))
C     ***  RECOVER STEP AND STORE PERMUTED -D*STEP AT W(RES)  ***
      DO 290 I = 1, P
         J1 = IPIVOT(I)
         K = RES0 + I
         T = W(K)
         STEP(J1) = -T
         W(K) = T*D(J1)
 290     CONTINUE
      DST = V2NORM(P, W(RES))
      PHI = DST - RAD
      IF (PHI .LE. PHIMAX .AND. PHI .GE. PHIMIN) GO TO 430
      IF (OLDPHI .EQ. PHI) GO TO 430
      OLDPHI = PHI
C
C  ***  CHECK FOR (AND HANDLE) SPECIAL CASE  ***
C
      IF (PHI .GT. ZERO) GO TO 310
         IF (KA .GE. KALIM) GO TO 430
              TWOPSI = ALPHAK*DST*DST - DOTPRD(P, STEP, G)
              IF (ALPHAK .GE. TWOPSI*PSIFAC) GO TO 310
                   V(STPPAR) = -ALPHAK
                   GO TO 440
C
C  ***  UNACCEPTABLE STEP -- UPDATE LK, UK, ALPHAK, AND TRY AGAIN  ***
C
 300  IF (PHI .LT. ZERO) UK = MIN(UK, ALPHAK)
      GO TO 320
 310  IF (PHI .LT. ZERO) UK = ALPHAK
 320  DO 330 I = 1, P
         J1 = IPIVOT(I)
         K = RES0 + I
         STEP(I) = D(J1) * (W(K)/DST)
 330     CONTINUE
      CALL LIVMUL(P, STEP, W(RMAT), STEP)
      DO 340 I = 1, P
 340     STEP(I) = STEP(I) / SQRT(W(I))
      T = ONE / V2NORM(P, STEP)
      ALPHAK = ALPHAK + T*PHI*T/RAD
      LK = MAX(LK, ALPHAK)
      GO TO 110
C
C  ***  RESTART  ***
C
 370  LK = W(LK0)
      UK = W(UK0)
      IF (V(DST0) .GT. ZERO .AND. V(DST0) - RAD .LE. PHIMAX) GO TO 20
      ALPHAK = ABS(V(STPPAR))
      DST = W(DSTSAV)
      PHI = DST - RAD
      T = V(DGNORM)/RAD
      IF (RAD .GT. V(RAD0)) GO TO 380
C
C        ***  SMALLER RADIUS  ***
         UK = T
         IF (ALPHAK .LE. ZERO) LK = ZERO
         IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
         GO TO 300
C
C     ***  BIGGER RADIUS  ***
 380  IF (ALPHAK .LE. ZERO .OR. UK .GT. T) UK = T
      LK = ZERO
      IF (V(DST0) .GT. ZERO) LK = MAX(LK, (V(DST0)-RAD)*W(PHIPIN))
      GO TO 300
C
C  ***  SPECIAL CASE -- RAD .LE. 0 OR (G = 0 AND J IS SINGULAR)  ***
C
 390  V(STPPAR) = ZERO
      DST = ZERO
      LK = ZERO
      UK = ZERO
      V(GTSTEP) = ZERO
      V(PREDUC) = ZERO
      DO 400 I = 1, P
 400     STEP(I) = ZERO
      GO TO 450
C
C  ***  ACCEPTABLE GAUSS-NEWTON STEP -- RECOVER STEP FROM W  ***
C
 410  ALPHAK = ZERO
      DO 420 I = 1, P
         J1 = IPIVOT(I)
         STEP(J1) = -W(I)
 420     CONTINUE
C
C  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  ***
C
 430  V(STPPAR) = ALPHAK
 440  V(GTSTEP) = DOTPRD(P, STEP, G)
      V(PREDUC) = HALF * (ALPHAK*DST*DST - V(GTSTEP))
 450  V(DSTNRM) = DST
      W(DSTSAV) = DST
      W(LK0) = LK
      W(UK0) = UK
      V(RAD0) = RAD
C
      RETURN
C
C  ***  LAST CARD OF LMSTEP FOLLOWS  ***
      END

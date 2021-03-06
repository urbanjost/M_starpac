*ASSESS
      SUBROUTINE ASSESS (D, IV, P, STEP, STLSTG, V, X, X0)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  ASSESS CANDIDATE STEP (NL2SOL VERSION 2.2)  ***
C
C  ***  PURPOSE  ***
C
C        THIS SUBROUTINE IS CALLED BY AN UNCONSTRAINED MINIMIZATION
C     ROUTINE TO ASSESS THE NEXT CANDIDATE STEP.  IT MAY RECOMMEND ONE
C     OF SEVERAL COURSES OF ACTION, SUCH AS ACCEPTING THE STEP, RECOM-
C     PUTING IT USING THE SAME OR A NEW QUADRATIC MODEL, OR HALTING DUE
C     TO CONVERGENCE OR FALSE CONVERGENCE.  SEE THE RETURN CODE LISTING
C     BELOW.
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
     +   D(P),STEP(P),STLSTG(P),V(*),X(P),X0(P)
      INTEGER
     +   IV(*)
C
C  LOCAL SCALARS
      REAL
     +   EMAX,GTS,HALF,ONE,RELDX1,RFAC1,TEMP,TWO,XMAX,ZERO
      INTEGER
     +   AFCTOL,DECFAC,DST0,DSTNRM,DSTSAV,F,F0,FDIF,FLSTGD,GTSLST,
     +   GTSTEP,I,INCFAC,IRC,LMAX0,MLSTGD,MODEL,NFC,NFCALL,NFGCAL,
     +   NREDUC,PLSTGD,PREDUC,RADFAC,RADINC,RDFCMN,RDFCMX,RELDX,
     +   RESTOR,RFCTOL,STAGE,STGLIM,STPPAR,SWITCH,TOOBIG,TUNER1,
     +   TUNER2,TUNER3,XCTOL,XFTOL,XIRC
      LOGICAL
     +   GOODX
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH,RELDST
      EXTERNAL R1MACH,RELDST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL VCOPY
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C--------------------------  PARAMETER USAGE  --------------------------
C
C     IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
C             BELOW OF IV VALUES REFERENCED.
C      D (IN)  SCALE VECTOR USED IN COMPUTING V(RELDX) -- SEE BELOW.
C      P (IN)  NUMBER OF PARAMETERS BEING OPTIMIZED.
C   STEP (I/O) ON INPUT, STEP IS THE STEP TO BE ASSESSED.  IT IS UN-
C             CHANGED ON OUTPUT UNLESS A PREVIOUS STEP ACHIEVED A
C             BETTER OBJECTIVE FUNCTION REDUCTION, IN WHICH CASE STLSTG
C             WILL HAVE BEEN COPIED TO STEP.
C STLSTG (I/O) WHEN ASSESS RECOMMENDS RECOMPUTING STEP EVEN THOUGH THE
C             CURRENT (OR A PREVIOUS) STEP YIELDS AN OBJECTIVE FUNC-
C             TION DECREASE, IT SAVES IN STLSTG THE STEP THAT GAVE THE
C             BEST FUNCTION REDUCTION SEEN SO FAR (IN THE CURRENT ITERA-
C             TION).  IF THE RECOMPUTED STEP YIELDS A LARGER FUNCTION
C             VALUE, THEN STEP IS RESTORED FROM STLSTG AND
C             X = X0 + STEP IS RECOMPUTED.
C      V (I/O) REAL PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
C             BELOW OF V VALUES REFERENCED.
C      X (I/O) ON INPUT, X = X0 + STEP IS THE POINT AT WHICH THE OBJEC-
C             TIVE FUNCTION HAS JUST BEEN EVALUATED.  IF AN EARLIER
C             STEP YIELDED A BIGGER FUNCTION DECREASE, THEN X IS
C             RESTORED TO THE CORRESPONDING EARLIER VALUE.  OTHERWISE,
C             IF THE CURRENT STEP DOES NOT GIVE ANY FUNCTION DECREASE,
C             THEN X IS RESTORED TO X0.
C     X0 (IN)  INITIAL OBJECTIVE FUNCTION PARAMETER VECTOR (AT THE
C             START OF THE CURRENT ITERATION).
C
C  ***  IV VALUES REFERENCED  ***
C
C    IV(IRC) (I/O) ON INPUT FOR THE FIRST STEP TRIED IN A NEW ITERATION,
C             IV(IRC) SHOULD BE SET TO 3 OR 4 (THE VALUE TO WHICH IT IS
C             SET WHEN STEP IS DEFINITELY TO BE ACCEPTED).  ON INPUT
C             AFTER STEP HAS BEEN RECOMPUTED, IV(IRC) SHOULD BE
C             UNCHANGED SINCE THE PREVIOUS RETURN OF ASSESS.
C                ON OUTPUT, IV(IRC) IS A RETURN CODE HAVING ONE OF THE
C             FOLLOWING VALUES...
C                  1 = SWITCH MODELS OR TRY SMALLER STEP.
C                  2 = SWITCH MODELS OR ACCEPT STEP.
C                  3 = ACCEPT STEP AND DETERMINE V(RADFAC) BY GRADIENT
C                       TESTS.
C                  4 = ACCEPT STEP, V(RADFAC) HAS BEEN DETERMINED.
C                  5 = RECOMPUTE STEP (USING THE SAME MODEL).
C                  6 = RECOMPUTE STEP WITH RADIUS = V(LMAX0) BUT DO NOT
C                       EVAULATE THE OBJECTIVE FUNCTION.
C                  7 = X-CONVERGENCE (SEE V(XCTOL)).
C                  8 = RELATIVE FUNCTION CONVERGENCE (SEE V(RFCTOL)).
C                  9 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE.
C                 10 = ABSOLUTE FUNCTION CONVERGENCE (SEE V(AFCTOL)).
C                 11 = SINGULAR CONVERGENCE (SEE V(LMAX0)).
C                 12 = FALSE CONVERGENCE (SEE V(XFTOL)).
C                 13 = IV(IRC) WAS OUT OF RANGE ON INPUT.
C             RETURN CODE I HAS PRECDENCE OVER I+1 FOR I = 9, 10, 11.
C IV(MLSTGD) (I/O) SAVED VALUE OF IV(MODEL).
C  IV(MODEL) (I/O) ON INPUT, IV(MODEL) SHOULD BE AN INTEGER IDENTIFYING
C             THE CURRENT QUADRATIC MODEL OF THE OBJECTIVE FUNCTION.
C             IF A PREVIOUS STEP YIELDED A BETTER FUNCTION REDUCTION,
C             THEN IV(MODEL) WILL BE SET TO IV(MLSTGD) ON OUTPUT.
C IV(NFCALL) (IN)  INVOCATION COUNT FOR THE OBJECTIVE FUNCTION.
C IV(NFGCAL) (I/O) VALUE OF IV(NFCALL) AT STEP THAT GAVE THE BIGGEST
C             FUNCTION REDUCTION THIS ITERATION.  IV(NFGCAL) REMAINS
C             UNCHANGED UNTIL A FUNCTION REDUCTION IS OBTAINED.
C IV(RADINC) (I/O) THE NUMBER OF RADIUS INCREASES (OR MINUS THE NUMBER
C             OF DECREASES) SO FAR THIS ITERATION.
C IV(RESTOR) (OUT) SET TO 0 UNLESS X AND V(F) HAVE BEEN RESTORED, IN
C             WHICH CASE ASSESS SETS IV(RESTOR) = 1.
C  IV(STAGE) (I/O) COUNT OF THE NUMBER OF MODELS TRIED SO FAR IN THE
C             CURRENT ITERATION.
C IV(STGLIM) (IN)  MAXIMUM NUMBER OF MODELS TO CONSIDER.
C IV(SWITCH) (OUT) SET TO 0 UNLESS A NEW MODEL IS BEING TRIED AND IT
C             GIVES A SMALLER FUNCTION VALUE THAN THE PREVIOUS MODEL,
C             IN WHICH CASE ASSESS SETS IV(SWITCH) = 1.
C IV(TOOBIG) (IN)  IS NONZERO IF STEP WAS TOO BIG (E.G. IF IT CAUSED
C             OVERFLOW).
C   IV(XIRC) (I/O) VALUE THAT IV(IRC) WOULD HAVE IN THE ABSENCE OF
C             CONVERGENCE, FALSE CONVERGENCE, AND OVERSIZED STEPS.
C
C  ***  V VALUES REFERENCED  ***
C
C V(AFCTOL) (IN)  ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.  IF THE
C             ABSOLUTE VALUE OF THE CURRENT FUNCTION VALUE V(F) IS LESS
C             THAN V(AFCTOL), THEN ASSESS RETURNS WITH IV(IRC) = 10.
C V(DECFAC) (IN)  FACTOR BY WHICH TO DECREASE RADIUS WHEN IV(TOOBIG) IS
C             NONZERO.
C V(DSTNRM) (IN)  THE 2-NORM OF D*STEP.
C V(DSTSAV) (I/O) VALUE OF V(DSTNRM) ON SAVED STEP.
C   V(DST0) (IN)  THE 2-NORM OF D TIMES THE NEWTON STEP (WHEN DEFINED,
C             I.E., FOR V(NREDUC) .GE. 0).
C      V(F) (I/O) ON BOTH INPUT AND OUTPUT, V(F) IS THE OBJECTIVE FUNC-
C             TION VALUE AT X.  IF X IS RESTORED TO A PREVIOUS VALUE,
C             THEN V(F) IS RESTORED TO THE CORRESPONDING VALUE.
C   V(FDIF) (OUT) THE FUNCTION REDUCTION V(F0) - V(F) (FOR THE OUTPUT
C             VALUE OF V(F) IF AN EARLIER STEP GAVE A BIGGER FUNCTION
C             DECREASE, AND FOR THE INPUT VALUE OF V(F) OTHERWISE).
C V(FLSTGD) (I/O) SAVED VALUE OF V(F).
C     V(F0) (IN)  OBJECTIVE FUNCTION VALUE AT START OF ITERATION.
C V(GTSLST) (I/O) VALUE OF V(GTSTEP) ON SAVED STEP.
C V(GTSTEP) (IN)  INNER PRODUCT BETWEEN STEP AND GRADIENT.
C V(INCFAC) (IN)  MINIMUM FACTOR BY WHICH TO INCREASE RADIUS.
C  V(LMAX0) (IN)  MAXIMUM REASONABLE STEP SIZE (AND INITIAL STEP BOUND).
C             IF THE ACTUAL FUNCTION DECREASE IS NO MORE THAN TWICE
C             WHAT WAS PREDICTED, IF A RETURN WITH IV(IRC) = 7, 8, 9,
C             OR 10 DOES NOT OCCUR, IF V(DSTNRM) .GT. V(LMAX0), AND IF
C             V(PREDUC) .LE. V(RFCTOL) * ABS(V(F0)), THEN ASSESS RE-
C             TURNS WITH IV(IRC) = 11.  IF SO DOING APPEARS WORTHWHILE,
C             THEN ASSESS REPEATS THIS TEST WITH V(PREDUC) COMPUTED FOR
C             A STEP OF LENGTH V(LMAX0) (BY A RETURN WITH IV(IRC) = 6).
C V(NREDUC) (I/O)  FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
C             NEWTON STEP.  IF ASSESS IS CALLED WITH IV(IRC) = 6, I.E.,
C             IF V(PREDUC) HAS BEEN COMPUTED WITH RADIUS = V(LMAX0) FOR
C             USE IN THE SINGULAR CONVERVENCE TEST, THEN V(NREDUC) IS
C             SET TO -V(PREDUC) BEFORE THE LATTER IS RESTORED.
C V(PLSTGD) (I/O) VALUE OF V(PREDUC) ON SAVED STEP.
C V(PREDUC) (I/O) FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
C             CURRENT STEP.
C V(RADFAC) (OUT) FACTOR TO BE USED IN DETERMINING THE NEW RADIUS,
C             WHICH SHOULD BE V(RADFAC)*DST, WHERE  DST  IS EITHER THE
C             OUTPUT VALUE OF V(DSTNRM) OR THE 2-NORM OF
C             DIAG(NEWD)*STEP  FOR THE OUTPUT VALUE OF STEP AND THE
C             UPDATED VERSION, NEWD, OF THE SCALE VECTOR D.  FOR
C             IV(IRC) = 3, V(RADFAC) = 1.0 IS RETURNED.
C V(RDFCMN) (IN)  MINIMUM VALUE FOR V(RADFAC) IN TERMS OF THE INPUT
C             VALUE OF V(DSTNRM) -- SUGGESTED VALUE = 0.1.
C V(RDFCMX) (IN)  MAXIMUM VALUE FOR V(RADFAC) -- SUGGESTED VALUE = 4.0.
C  V(RELDX) (OUT) SCALED RELATIVE CHANGE IN X CAUSED BY STEP, COMPUTED
C             BY FUNCTION  RELDST  AS
C                 MAX (D(I)*ABS(X(I)-X0(I)), 1 .LE. I .LE. P) /
C                    MAX (D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P).
C             IF AN ACCEPTABLE STEP IS RETURNED, THEN V(RELDX) IS COM-
C             PUTED USING THE OUTPUT (POSSIBLY RESTORED) VALUES OF X
C             AND STEP.  OTHERWISE IT IS COMPUTED USING THE INPUT
C             VALUES.
C V(RFCTOL) (IN)  RELATIVE FUNCTION CONVERGENCE TOLERANCE.  IF THE
C             ACTUAL FUNCTION REDUCTION IS AT MOST TWICE WHAT WAS PRE-
C             DICTED AND  V(NREDUC) .LE. V(RFCTOL)*ABS(V(F0)),  THEN
C             ASSESS RETURNS WITH IV(IRC) = 8 OR 9.  SEE ALSO V(LMAX0).
C V(STPPAR) (IN)  MARQUARDT PARAMETER -- 0 MEANS FULL NEWTON STEP.
C V(TUNER1) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
C             REDUCTION WAS MUCH LESS THAN EXPECTED.  SUGGESTED
C             VALUE = 0.1.
C V(TUNER2) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
C             REDUCTION WAS LARGE ENOUGH TO ACCEPT STEP.  SUGGESTED
C             VALUE = 10**-4.
C V(TUNER3) (IN)  TUNING CONSTANT USED TO DECIDE IF THE RADIUS
C             SHOULD BE INCREASED.  SUGGESTED VALUE = 0.75.
C  V(XCTOL) (IN)  X-CONVERGENCE CRITERION.  IF STEP IS A NEWTON STEP
C             (V(STPPAR) = 0) HAVING V(RELDX) .LE. V(XCTOL) AND GIVING
C             AT MOST TWICE THE PREDICTED FUNCTION DECREASE, THEN
C             ASSESS RETURNS IV(IRC) = 7 OR 9.
C  V(XFTOL) (IN)  FALSE CONVERGENCE TOLERANCE.  IF STEP GAVE NO OR ONLY
C             A SMALL FUNCTION DECREASE AND V(RELDX) .LE. V(XFTOL),
C             THEN ASSESS RETURNS WITH IV(IRC) = 12.
C
C-------------------------------  NOTES  -------------------------------
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C        THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR
C     LEAST-SQUARES) PACKAGE.  IT MAY BE USED IN ANY UNCONSTRAINED
C     MINIMIZATION SOLVER THAT USES DOGLEG, GOLDFELD-QUANDT-TROTTER,
C     OR LEVENBERG-MARQUARDT STEPS.
C
C  ***  ALGORITHM NOTES  ***
C
C        SEE (1) FOR FURTHER DISCUSSION OF THE ASSESSING AND MODEL
C     SWITCHING STRATEGIES.  WHILE NL2SOL CONSIDERS ONLY TWO MODELS,
C     ASSESS IS DESIGNED TO HANDLE ANY NUMBER OF MODELS.
C
C  ***  USAGE NOTES  ***
C
C        ON THE FIRST CALL OF AN ITERATION, ONLY THE I/O VARIABLES
C     STEP, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), AND
C     V(PREDUC) NEED HAVE BEEN INITIALIZED.  BETWEEN CALLS, NO I/O
C     VALUES EXECPT STEP, X, IV(MODEL), V(F) AND THE STOPPING TOLER-
C     ANCES SHOULD BE CHANGED.
C        AFTER A RETURN FOR CONVERGENCE OR FALSE CONVERGENCE, ONE CAN
C     CHANGE THE STOPPING TOLERANCES AND CALL ASSESS AGAIN, IN WHICH
C     CASE THE STOPPING TESTS WILL BE REPEATED.
C
C  ***  REFERENCES  ***
C
C     (1) DENNIS, J.E., JR., GAY, D.M., AND WELSCH, R.E. (1980),
C        AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM,
C        SUBMITTED TO ACM TRANS. MATH. SOFTWARE.
C
C     (2) POWELL, M.J.D. (1970)  A FORTRAN SUBROUTINE FOR SOLVING
C        SYSTEMS OF NONLINEAR ALGEBRAIC EQUATIONS, IN NUMERICAL
C        METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS, EDITED BY
C        P. RABINOWITZ, GORDON AND BREACH, LONDON.
C
C  ***  HISTORY  ***
C
C        JOHN DENNIS DESIGNED MUCH OF THIS ROUTINE, STARTING WITH
C     IDEAS IN (2). ROY WELSCH SUGGESTED THE MODEL SWITCHING STRATEGY.
C        DAVID GAY AND STEPHEN PETERS CAST THIS SUBROUTINE INTO A MORE
C     PORTABLE FORM (WINTER 1977), AND DAVID GAY CAST IT INTO ITS
C     PRESENT FORM (FALL 1978).
C
C  ***  GENERAL  ***
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
C     MCS-7906671.
C
C------------------------  EXTERNAL QUANTITIES  ------------------------
C
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL RELDST, VCOPY
C     REAL R1MACH, RELDST
C
C VCOPY.... COPIES ONE VECTOR TO ANOTHER.
C
C/
C  ***  NO COMMON BLOCKS  ***
C
C--------------------------  LOCAL VARIABLES  --------------------------
C
C     LOGICAL GOODX
C     INTEGER I, NFC
C     REAL EMAX, GTS, HALF, ONE, RELDX1, RFAC1,
C    +                 TEMP, TWO, XMAX, ZERO
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER AFCTOL, DECFAC, DSTNRM, DSTSAV, DST0, F, FDIF, FLSTGD, F0,
C    1        GTSLST, GTSTEP, INCFAC, IRC, LMAX0, MLSTGD, MODEL, NFCALL,
C    2        NFGCAL, NREDUC, PLSTGD, PREDUC, RADFAC, RADINC, RDFCMN,
C    3        RDFCMX, RELDX, RESTOR, RFCTOL, STAGE, STGLIM, STPPAR,
C    4        SWITCH, TOOBIG, TUNER1, TUNER2, TUNER3, XCTOL, XFTOL,
C    5        XIRC
C
C  ***  DATA INITIALIZATIONS  ***
C
      DATA HALF/0.5E0/, ONE/1.0E0/, TWO/2.0E0/, ZERO/0.0E0/
C
      DATA IRC/3/, MLSTGD/4/, MODEL/5/, NFCALL/6/,
     +     NFGCAL/7/, RADINC/8/, RESTOR/9/, STAGE/10/,
     +     STGLIM/11/, SWITCH/12/, TOOBIG/2/, XIRC/13/
      DATA AFCTOL/31/, DECFAC/22/, DSTNRM/2/, DST0/3/,
     +     DSTSAV/18/, F/10/, FDIF/11/, FLSTGD/12/, F0/13/,
     +     GTSLST/14/, GTSTEP/4/, INCFAC/23/,
     +     LMAX0/35/, NREDUC/6/, PLSTGD/15/, PREDUC/7/,
     +     RADFAC/16/, RDFCMN/24/, RDFCMX/25/,
     +     RELDX/17/, RFCTOL/32/, STPPAR/5/, TUNER1/26/,
     +     TUNER2/27/, TUNER3/28/, XCTOL/33/, XFTOL/34/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      NFC = IV(NFCALL)
      IV(SWITCH) = 0
      IV(RESTOR) = 0
      RFAC1 = ONE
      GOODX = .TRUE.
      I = IV(IRC)
      IF (I .GE. 1 .AND. I .LE. 12)
     +             GO TO (20,30,10,10,40,360,290,290,290,290,290,140), I
         IV(IRC) = 13
         GO TO 999
C
C  ***  INITIALIZE FOR NEW ITERATION  ***
C
 10   IV(STAGE) = 1
      IV(RADINC) = 0
      V(FLSTGD) = V(F0)
      IF (IV(TOOBIG) .EQ. 0) GO TO 90
         IV(STAGE) = -1
         IV(XIRC) = I
         GO TO 60
C
C  ***  STEP WAS RECOMPUTED WITH NEW MODEL OR SMALLER RADIUS  ***
C  ***  FIRST DECIDE WHICH  ***
C
 20   IF (IV(MODEL) .NE. IV(MLSTGD)) GO TO 30
C        ***  OLD MODEL RETAINED, SMALLER RADIUS TRIED  ***
C        ***  DO NOT CONSIDER ANY MORE NEW MODELS THIS ITERATION  ***
         IV(STAGE) = IV(STGLIM)
         IV(RADINC) = -1
         GO TO 90
C
C  ***  A NEW MODEL IS BEING TRIED.  DECIDE WHETHER TO KEEP IT.  ***
C
 30   IV(STAGE) = IV(STAGE) + 1
C
C     ***  NOW WE ADD THE POSSIBILTIY THAT STEP WAS RECOMPUTED WITH  ***
C     ***  THE SAME MODEL, PERHAPS BECAUSE OF AN OVERSIZED STEP.     ***
C
 40   IF (IV(STAGE) .GT. 0) GO TO 50
C
C        ***  STEP WAS RECOMPUTED BECAUSE IT WAS TOO BIG.  ***
C
         IF (IV(TOOBIG) .NE. 0) GO TO 60
C
C        ***  RESTORE IV(STAGE) AND PICK UP WHERE WE LEFT OFF.  ***
C
         IV(STAGE) = -IV(STAGE)
         I = IV(XIRC)
         GO TO (20, 30, 90, 90, 70), I
C
 50   IF (IV(TOOBIG) .EQ. 0) GO TO 70
C
C  ***  HANDLE OVERSIZE STEP  ***
C
      IF (IV(RADINC) .GT. 0) GO TO 80
         IV(STAGE) = -IV(STAGE)
         IV(XIRC) = IV(IRC)
C
 60      V(RADFAC) = V(DECFAC)
         IV(RADINC) = IV(RADINC) - 1
         IV(IRC) = 5
         GO TO 999
C
 70   IF (V(F) .LT. V(FLSTGD)) GO TO 90
C
C     *** THE NEW STEP IS A LOSER.  RESTORE OLD MODEL.  ***
C
      IF (IV(MODEL) .EQ. IV(MLSTGD)) GO TO 80
         IV(MODEL) = IV(MLSTGD)
         IV(SWITCH) = 1
C
C     ***  RESTORE STEP, ETC. ONLY IF A PREVIOUS STEP DECREASED V(F).
C
 80   IF (V(FLSTGD) .GE. V(F0)) GO TO 90
         IV(RESTOR) = 1
         V(F) = V(FLSTGD)
         V(PREDUC) = V(PLSTGD)
         V(GTSTEP) = V(GTSLST)
         IF (IV(SWITCH) .EQ. 0) RFAC1 = V(DSTNRM) / V(DSTSAV)
         V(DSTNRM) = V(DSTSAV)
         NFC = IV(NFGCAL)
         GOODX = .FALSE.
C
C
C  ***  COMPUTE RELATIVE CHANGE IN X BY CURRENT STEP  ***
C
 90   RELDX1 = RELDST(P, D, X, X0)
C
C  ***  RESTORE X AND STEP IF NECESSARY  ***
C
      IF (GOODX) GO TO 105
      DO 100 I = 1, P
         STEP(I) = STLSTG(I)
         X(I) = X0(I) + STLSTG(I)
 100     CONTINUE
C
 105  V(FDIF) = V(F0) - V(F)
      TEMP = 0.0
      IF (V(PREDUC).GT.R1MACH(1)/V(TUNER2)) TEMP = V(TUNER2) * V(PREDUC)
      IF (V(FDIF).GT.TEMP) GO TO 120
C
C        ***  NO (OR ONLY A TRIVIAL) FUNCTION DECREASE
C        ***  -- SO TRY NEW MODEL OR SMALLER RADIUS
C
         V(RELDX) = RELDX1
         IF (V(F) .LT. V(F0)) GO TO 110
              IV(MLSTGD) = IV(MODEL)
              V(FLSTGD) = V(F)
              V(F) = V(F0)
              CALL VCOPY(P, X, X0)
              IV(RESTOR) = 1
              GO TO 115
 110     IV(NFGCAL) = NFC
 115     IV(IRC) = 1
         IF (IV(STAGE) .LT. IV(STGLIM)) GO TO 130
              IV(IRC) = 5
              IV(RADINC) = IV(RADINC) - 1
              GO TO 130
C
C  ***  NONTRIVIAL FUNCTION DECREASE ACHIEVED  ***
C
 120  IV(NFGCAL) = NFC
      RFAC1 = ONE
      IF (GOODX) V(RELDX) = RELDX1
      V(DSTSAV) = V(DSTNRM)
      IF (V(FDIF) .GT. V(PREDUC)*V(TUNER1)) GO TO 200
C
C  ***  DECREASE WAS MUCH LESS THAN PREDICTED -- EITHER CHANGE MODELS
C  ***  OR ACCEPT STEP WITH DECREASED RADIUS.
C
      IF (IV(STAGE) .GE. IV(STGLIM)) GO TO 125
C        ***  CONSIDER SWITCHING MODELS  ***
         IV(IRC) = 2
         GO TO 130
C
C     ***  ACCEPT STEP WITH DECREASED RADIUS  ***
C
 125  IV(IRC) = 4
C
C  ***  SET V(RADFAC) TO FLETCHER*S DECREASE FACTOR  ***
C
 130  IV(XIRC) = IV(IRC)
      EMAX = V(GTSTEP) + V(FDIF)
      V(RADFAC) = HALF * RFAC1
      IF (EMAX .LT. V(GTSTEP)) V(RADFAC) = RFAC1 * MAX(V(RDFCMN),
     +                                           HALF * V(GTSTEP)/EMAX)
C
C  ***  DO FALSE CONVERGENCE TEST  ***
C
 140  IF (V(RELDX) .LE. V(XFTOL)) GO TO 160
         IV(IRC) = IV(XIRC)
         IF (V(F) .LT. V(F0)) GO TO 230
              GO TO 300
C
 160  IV(IRC) = 12
      GO TO 310
C
C  ***  HANDLE GOOD FUNCTION DECREASE  ***
C
 200  IF (V(FDIF) .LT. (-V(TUNER3) * V(GTSTEP))) GO TO 260
C
C     ***  INCREASING RADIUS LOOKS WORTHWHILE.  SEE IF WE JUST
C     ***  RECOMPUTED STEP WITH A DECREASED RADIUS OR RESTORED STEP
C     ***  AFTER RECOMPUTING IT WITH A LARGER RADIUS.
C
      IF (IV(RADINC) .LT. 0) GO TO 260
      IF (IV(RESTOR) .EQ. 1) GO TO 260
C
C        ***  WE DID NOT.  TRY A LONGER STEP UNLESS THIS WAS A NEWTON
C        ***  STEP.
C
         V(RADFAC) = V(RDFCMX)
         GTS = V(GTSTEP)
         IF (V(FDIF) .LT. (HALF/V(RADFAC) - ONE) * GTS)
     +            V(RADFAC) = MAX(V(INCFAC), HALF*GTS/(GTS + V(FDIF)))
         IV(IRC) = 4
         IF (V(STPPAR) .EQ. ZERO) GO TO 300
C             ***  STEP WAS NOT A NEWTON STEP.  RECOMPUTE IT WITH
C             ***  A LARGER RADIUS.
              IV(IRC) = 5
              IV(RADINC) = IV(RADINC) + 1
C
C  ***  SAVE VALUES CORRESPONDING TO GOOD STEP  ***
C
 230  V(FLSTGD) = V(F)
      IV(MLSTGD) = IV(MODEL)
      CALL VCOPY(P, STLSTG, STEP)
      V(DSTSAV) = V(DSTNRM)
      IV(NFGCAL) = NFC
      V(PLSTGD) = V(PREDUC)
      V(GTSLST) = V(GTSTEP)
      GO TO 300
C
C  ***  ACCEPT STEP WITH RADIUS UNCHANGED  ***
C
 260  V(RADFAC) = ONE
      IV(IRC) = 3
      GO TO 300
C
C  ***  COME HERE FOR A RESTART AFTER CONVERGENCE  ***
C
 290  IV(IRC) = IV(XIRC)
      IF (V(DSTSAV) .GE. ZERO) GO TO 310
         IV(IRC) = 12
         GO TO 310
C
C  ***  PERFORM CONVERGENCE TESTS  ***
C
 300  IV(XIRC) = IV(IRC)
 310  IF (ABS(V(F)) .LT. V(AFCTOL)) IV(IRC) = 10
      IF (HALF * V(FDIF) .GT. V(PREDUC)) GO TO 999
      EMAX = 0.0
      IF (ABS(V(F0)).GT.R1MACH(1)/V(RFCTOL))
     +   EMAX = V(RFCTOL) * ABS(V(F0))
      IF (V(DSTNRM) .GT. V(LMAX0) .AND. V(PREDUC) .LE. EMAX)
     +                       IV(IRC) = 11
      IF (V(DST0) .LT. ZERO) GO TO 320
      I = 0
      IF ((V(NREDUC) .GT. ZERO .AND. V(NREDUC) .LE. EMAX) .OR.
     +    (V(NREDUC) .EQ. ZERO. AND. V(PREDUC) .EQ. ZERO))  I = 2
      IF (V(STPPAR) .EQ. ZERO .AND. V(RELDX) .LE. V(XCTOL)) I = I + 1
      IF (I .GT. 0) IV(IRC) = I + 6
C
C  ***  CONSIDER RECOMPUTING STEP OF LENGTH V(LMAX0) FOR SINGULAR
C  ***  CONVERGENCE TEST.
C
 320  IF (ABS(IV(IRC)-3) .GT. 1 .AND. IV(IRC) .NE. 12) GO TO 999
      IF (V(DSTNRM) .GT. V(LMAX0)) GO TO 330
         IF (V(PREDUC) .GE. EMAX) GO TO 999
              IF (V(DST0) .LT. ZERO) GO TO 340
                   IF (HALF * V(DST0) .LE. V(LMAX0)) GO TO 999
                        GO TO 340
 330  IF (HALF * V(DSTNRM) .LE. V(LMAX0)) GO TO 999
      XMAX = V(LMAX0) / V(DSTNRM)
      IF (XMAX * (TWO - XMAX) * V(PREDUC) .GE. EMAX) GO TO 999
 340  IF (V(NREDUC) .LT. ZERO) GO TO 370
C
C  ***  RECOMPUTE V(PREDUC) FOR USE IN SINGULAR CONVERGENCE TEST  ***
C
      V(GTSLST) = V(GTSTEP)
      V(DSTSAV) = V(DSTNRM)
      IF (IV(IRC) .EQ. 12) V(DSTSAV) = -V(DSTSAV)
      V(PLSTGD) = V(PREDUC)
      IV(IRC) = 6
      CALL VCOPY(P, STLSTG, STEP)
      GO TO 999
C
C  ***  PERFORM SINGULAR CONVERGENCE TEST WITH RECOMPUTED V(PREDUC)  ***
C
 360  V(GTSTEP) = V(GTSLST)
      V(DSTNRM) = ABS(V(DSTSAV))
      CALL VCOPY(P, STEP, STLSTG)
      IV(IRC) = IV(XIRC)
      IF (V(DSTSAV) .LE. ZERO) IV(IRC) = 12
      V(NREDUC) = -V(PREDUC)
      V(PREDUC) = V(PLSTGD)
 370  IF (-V(NREDUC) .LE. V(RFCTOL) * ABS(V(F0))) IV(IRC) = 11
C
 999  RETURN
C
C  ***  LAST CARD OF ASSESS FOLLOWS  ***
      END

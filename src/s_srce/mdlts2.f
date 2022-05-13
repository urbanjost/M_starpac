*MDLTS2
      SUBROUTINE MDLTS2 (PAR, RESTS, Y, NPAR, N, NFAC, MSPECT, PMU,
     +  PARDF, NPARDF, T, TEMP, PARAR, PARMA, MBO, N1, N2, IFLAG)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MODEL ROUTINE FOR PACKS SPECIFICATION OF
C     BOX-JENKINS MODELS.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JANUARY 4, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   PMU
      INTEGER
     +   IFLAG,MBO,N,N1,N2,NFAC,NPAR,NPARDF
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),PARAR(*),PARDF(*),PARMA(*),RESTS(N1:N2),T(*),
     +   TEMP(*),Y(N)
      INTEGER
     +   MSPECT(NFAC,4)
C
C  LOCAL SCALARS
      REAL
     +   FPLPM,RESMAX,WTEST
      INTEGER
     +   I,IMOD,IMOD1,IPAR,IPQ,ISTART,J,K,L,MAXORD,MBO1,NP,NPARAR,
     +   NPARMA
      LOGICAL
     +   PARLE1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,LOG,MOD,SIGN,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL FPLPM
C        THE FLOATING POINT LARGEST POSITIVE MAGNITUDE.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IFLAG
C        AN INDICATOR VARIABLE DESIGNATING WHETHER THE BACK FORECASTS
C        WERE ESSENTIALLY ZERO (IFLAG=0) OR NOT (IFLAG=1).
C     INTEGER IMOD
C        AN INDEX VARIABLE.
C     INTEGER IPAR
C        AN INDEX VARIABLE.
C     INTEGER IPQ
C        AN INDEX VARIABLE.
C     INTEGER ISTART
C        ***
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER K
C        AN INDEX VARIABLE.
C     INTEGER L
C        AN INDEX VARIABLE.
C     INTEGER MAXORD
C        THE LARGEST BACK ORDER.
C     INTEGER MBO
C        THE MAXIMUM BACK ORDER OPERATOR.
C     INTEGER MBO1
C        THE VALUE MBO+1
C     INTEGER MSPECT(NFAC,4)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NP
C        THE NUMBER OF PARAMETERS IN THE EXPANDED TERM.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARAR
C        THE NUMBER OF AUTOREGRESSIVE PARAMETERS
C     INTEGER NPARDF
C        THE ORDER OF THE EXPANDED DIFFERENCE FILTER.
C     INTEGER NPARMA
C        THE LENGTH OF THE VECTOR PARMA
C     INTEGER N1
C        THE LOWER BOUND FOR RESTS.
C     INTEGER N2
C        THE UPPER BOUND FOR RESTS.
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     REAL PARAR(MBO)
C        THE AUTOREGRESSIVE PARAMETERS
C     REAL PARDF(NPARDF)
C        THE VECTOR CONTAINING THE DIFFERENCE FILTER PARAMETERS.
C     LOGICAL PARLE1
C        A FLAG INDICATING WHETHER ALL OF THE MOVING AVERAGE PARAMETERS
C        ARE LESS THAN OR EQUAL TO 1 (PARLE1 = .TRUE.) OR NOT
C        (PARLE1 = .FALSE.)
C     REAL PARMA(MBO)
C        THE MOVING AVERAGE PARAMETERS
C     REAL PMU
C        THE VALUE OF MU, I.E., THE TREND OR MEAN.
C     REAL RESMAX
C        THE LARGEST POSSIBLE RESIDUAL WHICH WILL STILL AVOID OVERFLOW.
C     REAL RESTS(N1:N2)
C        THE PREDICTED VALUE OF THE FIT.
C     REAL T(2*MBO)
C        A TEMPORARY WORK VECTOR.
C     REAL TEMP(MBO)
C        A TEMPORARY WORK VECTOR
C     REAL WTEST
C        THE TEST VALUE USED TO DETERMINE IF THE DIFFERENCED SERIES
C        BACK FORECAST IS EFFECTIVELY ZERO OR NOT.
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C
      FPLPM = R1MACH(2)
C
C     ZERO THE PARAMETER ARRAYS PARAR AND PARMA
C
      DO 10 I=1,MBO
         T(I) = 0.0E0
         TEMP(I) = 0.0E0
   10 CONTINUE
C
      NP = 0
      IPAR = 0
      NPARAR = 0
      ISTART = 0
C
C     EXPAND THE MODEL AND STORE AUTOREGRESSIVE PARAMETERS IN PARAR
C     AND MOVING AVERAGE PARAMETERS IN PARMA
C
      DO 110 IPQ = 1, 3, 2
         DO 100 L=1,NFAC
            IF (MSPECT(L,IPQ).EQ.0) GO TO 100
            MAXORD = MSPECT(L,IPQ)*MSPECT(L,4)
            DO 90 K = MSPECT(L,4), MAXORD, MSPECT(L,4)
               IPAR = IPAR + 1
               TEMP(K) = TEMP(K) + PAR(IPAR)
               DO 80 I = 1, NP
                  TEMP(K+I) = TEMP(K+I) - T(I)*PAR(IPAR)
   80          CONTINUE
   90       CONTINUE
            NP = NP + MAXORD
            DO 95 K = 1, NP
               T(K) = TEMP(K)
   95       CONTINUE
  100    CONTINUE
          IF (IPQ.NE.3) THEN
            IPAR = IPAR + 1
            PMU = PAR(IPAR)
            NPARAR = NP
            DO 105 K =1, NPARAR
               PARAR(K) = T(K)
               T(K) = 0.0E0
               TEMP(K) = 0.0E0
  105       CONTINUE
            NP = 0
         END IF
  110 CONTINUE
      NPARMA = NP
      PARLE1 = .TRUE.
      DO 115 K =1, NPARMA
         PARMA(K) = T(K)
         IF (ABS(PARMA(K)).GT.1.0E0) PARLE1 = .FALSE.
  115 CONTINUE
C
C     COMPUTE FITTED VALUES AND RESIDUALS FOR MODEL.
C
C     COMPUTE W, THE DIFFERENCED SERIES MINUS ITS MEAN, AND STORE IN
C     RESTS(NPARDF+1) TO RESTS(N2)
C
      DO 140 I = NPARDF+1, N2, 1
         RESTS(I) = Y(I) - PMU
         DO 130 J = 1,NPARDF
            RESTS(I) = RESTS(I) - PARDF(J)*Y(I-J)
  130    CONTINUE
  140 CONTINUE
      WTEST = ABS(RESTS(NPARDF+1))*0.01
C
C     BACK FORECAST THE ERROR, E, FOR I = N-NPARAR TO NPARDF+1, AND
C     THE DIFFERENCED SERIES FOR I = NPARDF TO N1
C
      MBO1 = MBO+1
      IFLAG = 0
      DO 170 I = N2-NPARAR,NPARDF+1,-1
         IMOD = MOD(I+1-N1,MBO1) + 1
         T(IMOD) = RESTS(I)
         DO 150 J = 1,NPARAR
            T(IMOD) = T(IMOD) - PARAR(J)*RESTS(I+J)
  150    CONTINUE
         DO 160 J = 1,NPARMA
            IF ((I+J.GT.NPARDF) .AND. (I+J.LE.N))
     +         T(IMOD) = T(IMOD) + PARMA(J)*T(MOD(I+J+1-N1,MBO1)+1)
  160    CONTINUE
  170 CONTINUE
      DO 175 I = NPARDF,N1,-1
         IMOD = MOD(I+1-N1,MBO1) + 1
         RESTS(I) = 0.0E0
         DO 163 J = 1,NPARAR
            RESTS(I) = RESTS(I) + PARAR(J)*RESTS(I+J)
  163    CONTINUE
         DO 166 J = 1,NPARMA
            IF ((I+J.GT.NPARDF) .AND. (I+J.LE.N))
     +         RESTS(I) = RESTS(I) -
     +                    PARMA(J)*T(MOD(I+J+1-N1,MBO1)+1)
  166    CONTINUE
         ISTART = I
         IF ((ISTART.LE.1) .AND. (ABS(RESTS(I)).LE.WTEST)) GO TO 180
  175 CONTINUE
      IFLAG = 1
C
C     COMPUTE RESIDUALS AND STORE VALUES IN RESTS
C
  180 CONTINUE
      DO 210 I = ISTART,N2,1
         IMOD = MOD(I+1-N1,MBO1) + 1
         T(IMOD) = RESTS(I)
         DO 190 J = 1,NPARAR
            IF (I-J.GE.ISTART) T(IMOD) = T(IMOD) - PARAR(J)*RESTS(I-J)
  190    CONTINUE
C
         IF (PARLE1) THEN
C
C     COMPUTE RESIDUALS WHERE THERE IS NO CHANCE OF OVERFLOW
C
            DO 200 J = 1,NPARMA
               IF (I-J.GE.ISTART)
     +            T(IMOD) = T(IMOD) + PARMA(J)*T(MOD(I-J+1-N1,MBO1)+1)
  200       CONTINUE
         ELSE
C
C     COMPUTE RESIDUALS WHERE THERE IS A CHANCE OF OVERFLOW
C
            DO 205 J = 1,NPARMA
               IF (I-J.GE.ISTART) THEN
                  IMOD1 = MOD(I-J+1-N1,MBO1)+1
                  IF (PARMA(J).NE.0.0E0 .AND. T(IMOD1).NE.0.0E0) THEN
                     IF (LOG(ABS(PARMA(J)))+LOG(ABS(T(IMOD1))).LT.
     +                         LOG(FPLPM)
     +                     .AND.
     +                     (SIGN(1.0E0,T(IMOD)).NE.
     +                         SIGN(1.0E0,PARMA(J)*T(IMOD1))
     +                     .OR.
     +                     LOG(ABS(PARMA(J)))+LOG(ABS(T(IMOD1))).LT.
     +                         LOG(FPLPM-ABS(T(IMOD))))) THEN
                        T(IMOD) = T(IMOD) + PARMA(J)*T(IMOD1)
                     ELSE
                        GO TO 300
                     END IF
                  END IF
               END IF
  205       CONTINUE
         END IF
         IF (I-MBO.GE.ISTART) THEN
            RESTS(I-MBO) = T(MOD(I-MBO+1-N1,MBO1)+1)
         END IF
  210 CONTINUE
      DO 220 I = N-MBO+1,N
        RESTS(I) = T(MOD(I-MBO+2-N1,MBO1)+1)
  220 CONTINUE
C
      DO 230 I = N1, ISTART-1
         RESTS(I) = 0.0E0
  230 CONTINUE
C
      RETURN
C
C     SET RESIDUALS TO LARGEST POSSIBLE VALUE
C
  300 RESMAX = SQRT(FPLPM/(N2-N1+1))
      DO 310 I=N1,N2
         RESTS(I) = RESMAX
  310 CONTINUE
C
      RETURN
C
      END

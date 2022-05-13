*MDLTS1
      SUBROUTINE MDLTS1 (PAR, NPAR, XM, N, M, IXM, RESTS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR ESTIMATING BOX-JENKINS
C     ARIMA MODELS.
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
      INTEGER
     +   IXM,M,N,NPAR
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(NPAR),RESTS(NRESTS),XM(IXM,M)
C
C  SCALARS IN COMMON
      INTEGER
     +   IFLAG,MBO,MBOL,MSPECT,NFACT,NPARAR,NPARDF,NPARMA,NRESTS,
     +   PARAR,PARDF,PARMA,T,TEMP
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   PMU
      INTEGER
     +   I,I1
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   RSTAK(12)
      INTEGER
     +   ISTAK(12)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL MDLTS2
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /MDLTSC/MSPECT,NFACT,PARDF,NPARDF,PARAR,NPARAR,PARMA,
     +   NPARMA,MBO,MBOL,T,TEMP,NRESTS,IFLAG
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (DSTAK(1),RSTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IFLAG
C        AN INDICATOR VARIABLE DESIGNATING WHETHER THE BACK FORECASTS
C        WERE ESSENTIALLY ZERO (IFLAG=0) OR NOT (IFLAG=1).
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IXM
C        THE FIRST DIMENSION OF MATRIX XM.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MBO
C        THE MAXIMUM BACK ORDER OPERATOR.
C     INTEGER MBOL
C        THE MAXIMUM BACK ORDER ON THE LEFT
C     INTEGER MSPECT
C        THE STARTING LOCATION IN THE WORK SPACE FOR
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFACT
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARAR
C        THE NUMBER OF AUTOREGRESSIVE PARAMETERS
C     INTEGER NPARDF
C        THE ORDER OF THE EXPANDED DIFFERENCE FILTER.
C     INTEGER NPARMA
C        THE LENGTH OF THE VECTOR PARMA
C     INTEGER NRESTS
C        THE MAXIMUM NUMBER OF RESIDUALS TO BE COMPUTED.
C     DOUBLE PRECISION PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     INTEGER PARAR
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        THE AUTOREGRESSIVE PARAMETERS
C     INTEGER PARDF
C        THE STARTING LOCATION IN THE WORK SPACE FOR
C        THE VECTOR CONTAINING THE DIFFERENCE FILTER PARAMETERS
C     INTEGER PARMA
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        THE MOVING AVERAGE PARAMETERS
C     DOUBLE PRECISION PMU
C        THE VALUE OF MU, I.E., THE TREND OR MEAN.
C     DOUBLE PRECISION RESTS(NRESTS)
C        THE RESIDUALS FROM THE ARIMA MODEL.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER T
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        A TEMPORARY WORK VECTOR.
C     INTEGER TEMP
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        A TEMPORARY WORK VECTOR
C     DOUBLE PRECISION XM(IXM,M)
C        THE INDEPENDENT VARIABLE.
C
C
C     COMPUTE RESIDUALS
C
      CALL MDLTS2 (PAR, RESTS, XM(1,1), NPAR, N, NFACT, ISTAK(MSPECT),
     +  PMU, RSTAK(PARDF), NPARDF, RSTAK(T), RSTAK(TEMP), RSTAK(PARAR),
     +  RSTAK(PARMA), MBO, N-NRESTS+1, N, IFLAG)
C
C     COMPUTE PREDICTED VALUES
C
      I1=NRESTS-N
      DO 20 I = 1,N
        I1=I1+1
        RESTS(I) = XM(I1,1)-RESTS(I1)
   20 CONTINUE
C
      RETURN
      END

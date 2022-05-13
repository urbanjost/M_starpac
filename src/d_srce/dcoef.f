*DCOEF
      SUBROUTINE DCOEF (NDF, ND, IOD, NPARDF, PARDF, MBO, WORK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE EXPANDS THE DIFFERENCE FILTER SPECIFIED BY NDF,
C     IOD AND ND INTO PARDF.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DEVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   MBO,NDF,NPARDF
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PARDF(*),WORK(*)
      INTEGER
     +   IOD(*),ND(*)
C
C  LOCAL SCALARS
      INTEGER
     +   K,KK,L,NTIMES,NWORK1,NWORK2
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   NCHOSE
      EXTERNAL NCHOSE
C
C  EXTERNAL SUBROUTINES
      EXTERNAL MULTBP
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IOD(NDF)
C        THE ORDER OF EACH OF THE DIFFERENCE FACTORS.
C     INTEGER K
C        AN INDEX VARIABLE.
C     INTEGER KK
C        AN INDEX VARIABLE.
C     INTEGER L
C        AN INDEX VARIABLE.
C     INTEGER MBO
C        THE MAXIMUM BACK ORDER OPERATOR.
C     INTEGER ND(NDF)
C        THE NUMBER OF TIMES EACH DIFFERENCE FACTOR IS TO BE APPLIED.
C     INTEGER NDF
C        THE NUMBER OF DIFFERENCE FACTORS
C     INTEGER NPARDF
C        THE ORDER OF THE EXPANDED DIFFERENCE FILTER.
C     INTEGER NTIMES
C        THE NUMBER OF TIMES A GIVEN DIFFERENCE FACTOR IS TO BE APPLIED.
C     INTEGER NWORK1
C        THE NUMBER OF TERMS IN THE FIRST COLUMN OF WORK.
C     INTEGER NWORK2
C        THE NUMBER OF TERMS IN THE SECOND COLUMN OF WORK
C     DOUBLE PRECISION PARDF(MBO)
C        THE VECTOR CONTAINING THE DIFFERENCE FILTER PARAMETERS.
C     DOUBLE PRECISION WORK(MBO,2)
C        A WORK ARRAY NECESSARY TO EXPAND THE DIFFERENCE FILTER.
C
      NPARDF = 0
C
      DO 30 L = 1, NDF
         IF (ND(L).EQ.0) GO TO 30
         NTIMES = ND(L)
         NWORK1 = IOD(L) * ND(L)
         DO 10 K = 1, NWORK1
            WORK(K) = 0.0D0
   10    CONTINUE
         DO 20 K = 1, NTIMES
            KK = K * IOD(L)
            WORK(KK) = ((-1)**(K+1)) * NCHOSE(NTIMES, K)
   20    CONTINUE
         NWORK2 = NWORK1 + NPARDF
         CALL MULTBP (WORK(1), NWORK1, PARDF, NPARDF, WORK(MBO+1),
     +      NWORK2, MBO)
   30 CONTINUE
      RETURN
      END

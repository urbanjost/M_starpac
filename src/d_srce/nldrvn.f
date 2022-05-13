*NLDRVN
      SUBROUTINE NLDRVN (MDL, DRV, DONE, IFIXD, PAR, NPAR, XM, N, M,
     +   IXM, PVT, D, WEIGHT, WT, LWT, STPT, LSTPT, SCL, LSCL)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE NUMERICAL APPROXIMATIONS TO THE
C     DERIVATIVE MATRIX (JACOBIAN).
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IXM,LSCL,LSTPT,LWT,M,N,NPAR
      LOGICAL
     +   DONE,WEIGHT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   D(N,NPAR),PAR(NPAR),PVT(N),SCL(LSCL),STPT(LSTPT),WT(LWT),
     +   XM(IXM,M)
      INTEGER
     +   IFIXD(NPAR)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL DRV,MDL
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   PJ,STPJ,WTSQRT
      INTEGER
     +   I,J,JPK
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX,SIGN,SQRT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION D(N,NPAR)
C        THE FIRST DERIVATIVE OF THE MODEL (JACOBIAN).
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        DERIVATIVE (JACOBIAN) MATRIX OF THE MODEL.
C     LOGICAL DONE
C        THE VARIABLE USED TO INDICATE WHETHER THIS IS THE FINAL
C        COMPUTATION OF THE JACOBIAN OR NOT.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IXM
C        THE FIRST DIMENSION OF MATRIX XM.
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER JPK
C        AN INDEX VARIABLE.
C     INTEGER LSCL
C        THE DIMENSION OF VECTOR SCL.
C     INTEGER LSTPT
C        THE DIMENSION OF VECTOR STPT.
C     INTEGER LWT
C        THE DIMENSION OF VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     DOUBLE PRECISION PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     DOUBLE PRECISION PJ
C        A TEMPORARY LOCATION FOR STORAGE OF THE JTH PARAMETER.
C     DOUBLE PRECISION PVT(N)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES.
C     DOUBLE PRECISION SCL(LSCL)
C        THE SCALE VALUES.
C     DOUBLE PRECISION STPT(LSTPT)
C        THE STEP SIZE ARRAY.
C     DOUBLE PRECISION STPJ
C        THE JTH STEP SIZE.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION WT(LWT)
C        THE USER SUPPLIED WEIGHTS.
C     DOUBLE PRECISION WTSQRT
C        THE SQUARE ROOT OF THE USER SUPPLIED WEIGHTS.
C     DOUBLE PRECISION XM(IXM,M)
C        THE INDEPENDENT VARIABLE.
C
C     COMPUTE FINITE-DIFFERENCE JACOBIAN OF THE OPTIMIZED PARAMETERS
C
      JPK = 0
C
      DO 20 J=1,NPAR
         IF (IFIXD(J).EQ.0) THEN
            JPK = JPK + 1
            PJ = PAR(J)
            IF (SCL(JPK).EQ.0.0D0) THEN
               IF (PAR(J).NE.0.0D0) THEN
                  STPJ = STPT(J)*SIGN(1.0D0,PAR(J))*ABS(PAR(J))
               ELSE
                  STPJ = STPT(J)
               END IF
            ELSE
               STPJ = STPT(J)*
     +                SIGN(1.0D0,PAR(J))*MAX(ABS(PAR(J)),1.0D0/
     +                ABS(SCL(JPK)))
            END IF
C
            STPJ = STPJ + PAR(J)
            STPJ = STPJ - PAR(J)
C
            PAR(J) = PJ + STPJ
            CALL MDL(PAR, NPAR, XM, N, M, IXM, D(1,J))
C
            DO 10 I=1,N
               WTSQRT = 1.0D0
               IF (WEIGHT .AND. (.NOT.DONE)) WTSQRT = SQRT(WT(I))
               D(I,JPK) = WTSQRT*(PVT(I)-D(I,J))/STPJ
   10       CONTINUE
C
            PAR(J) = PJ
         END IF
   20 CONTINUE
C
      RETURN
C
      END

*DRV3
      SUBROUTINE DRV3(PAR, NPAR, XM, N, M, IXM, D)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     DERIVATIVE FUNCTION FOR NLS FAMILY EXERCISER SUBROUTINE MDL3
C
C     WRITTEN BY  -  LINDA L. MITCHELL
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IXM,M,N,NPAR
C
C  ARRAY ARGUMENTS
      REAL
     +   D(N,NPAR),PAR(NPAR),XM(IXM,M)
C
C  LOCAL SCALARS
      INTEGER
     +   I,J
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL PAR(NPAR)
C        MODEL PARAMETERS
C     REAL D(N,NPAR)
C        THE FIRST DERIVATIVE WITH RESPECT TO THE ITH PARAMETER
C     INTEGER I
C        ROW MARKER
C     INTEGER IXM
C        ACTUAL FIRST DIMENSION OF XM
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLESC
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS
C     REAL XM(IXM,M)
C        INDEPENDENT VARIABLE
C
      DO 20 I=1,N
         DO 10 J=1,NPAR
            D(I,J) = XM(I,J)
   10    CONTINUE
   20 CONTINUE
C
      RETURN
C
      END

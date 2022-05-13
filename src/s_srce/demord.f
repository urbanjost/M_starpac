*DEMORD
      SUBROUTINE DEMORD (PHAS1, PHAS2, NDEM, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE SETS UP THE DATA FOR THE PHASE PLOTS.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NDEM
C
C  ARRAY ARGUMENTS
      REAL
     +   PHAS1(N),PHAS2(N)
C
C  LOCAL SCALARS
      REAL
     +   PI
      INTEGER
     +   I
C
C  EXTERNAL SUBROUTINES
      EXTERNAL GETPI
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES BEING DEMODULATED.
C     INTEGER NDEM
C        THE NUMBER OF VALUES IN THE DEMODULATED SERIES, I.E., IN
C        THE AMPLITUDE AND PHASE ARRAYS.
C     REAL PHAS1(N), PHAS2(N)
C        THE ARRAYS CONTAINING THE PRIMARY AND SECONDARY PHASE
C        ESTIMATES, RESPECTIVELY.
C     REAL PI
C        THE VALUE OF PI.
C
      CALL GETPI(PI)
C
      DO 10 I = 1, NDEM
         PHAS2(I) = 0.0E0
         IF (PHAS1(I) .GT. 0.0E0) PHAS2(I) = PHAS1(I) - 2.0E0*PI
         IF (PHAS1(I) .LT. 0.0E0) PHAS2(I) = PHAS1(I) + 2.0E0*PI
   10 CONTINUE
C
      RETURN
      END

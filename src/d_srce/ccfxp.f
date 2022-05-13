*CCFXP
      SUBROUTINE CCFXP (STORE, LAGMAX, M, CCOV, ICCOV, JCCOV, MISS,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO LIST THE COMPUTED RESULTS FROM THE TIME SERIES
C     CROSS CORRELATION ROUTINES.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   CMISS
      INTEGER
     +   ICCOV,INLPPC,JCCOV,JNLPPC,LAGMAX,M
      LOGICAL
     +   MISS,STORE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   CCOV(ICCOV,JCCOV,*)
      INTEGER
     +   NLPPC(INLPPC,JNLPPC,*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT,J,K,L,L1,LAG
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   CCF(16)
C
C  EXTERNAL FUNCTIONS
      LOGICAL
     +   MVCHK
      EXTERNAL MVCHK
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CCF(16)
C        AN ARRAY USED FOR PRINTING THE CCF.
C     DOUBLE PRECISION CCOV(ICCOV,JCCOV,M)
C        THE CROSS COVARIANCE ARRAY.
C     DOUBLE PRECISION CMISS
C        THE MISSING VALUE CODE FOR THE RETURNED CCVF ESTIMATES
C        (VECTOR CCOV).
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER ICCOV
C        THE FIRST DIMENSION OF THE ARRAY CCOV.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C     INTEGER INLPPC
C        THE FIRST DIMENSION OF THE ARRAY NLPPC.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER J
C        AN INDEXING VARIABLE.
C     INTEGER JCCOV, JNLPPC
C        THE SECOND DIMENSIONS OF THE ARRAYS CCOV AND NLPPC,
C        RESPECTIVELY.
C     INTEGER K
C        AN INDEXING VARIABLE.
C     INTEGER L1
C        AN INDEX VARIABLE.
C     INTEGER LAG
C        THE LAG VALUE AT WHICH THE DATA IS BEING PRINTED.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE REQUESTED.
C     INTEGER M
C        THE NUMBER OF SERIES IN THE MULTIVARIATE TIME SERIES YM.
C     LOGICAL MISS
C        THE VALUE INDICATING WHETHER THE ANALYSIS INCLUDED MISSING
C        DATA (TRUE) OR NOT (FALSE).
C     INTEGER NLPPC(INLPPC,JNLPPC,M)
C        THE ARRAY CONTAINING THE NUMBER OF LAGGED PRODUCT PAIRS
C        USED TO COMPUTE EACH CCVF ESTIMATE.
C     LOGICAL STORE
C        THE VALUE INDICATING WHETHER THE RESULTS WERE RETURNED
C        TO THE USER (TRUE) OR NOT (FALSE).
C
      CALL IPRINT(IPRT)
C
C     PRINT IERR
C
      WRITE (IPRT, 1000) IERR
C
      IF (IERR.NE.0) RETURN
C
C     CHECK FOR STORED RESULTS
C
      IF (.NOT.STORE) RETURN
C
C     PRINT HEADING FOR CCVF
C
      WRITE (IPRT, 1010)
      WRITE (IPRT, 1040) ((J,K, K=1,M), J=1,M)
C
C     PRINT CROSS COVARIANCES
C
      LAG = 0
      WRITE (IPRT, 1060) LAG, ((CCOV(1,J,K), K=1,M), J=1,M)
      DO 10 LAG = 1, LAGMAX
         WRITE (IPRT, 1060) LAG, ((CCOV(LAG+1,J,K), K=1,M), J=1,M)
   10 CONTINUE
C
C     PRINT HEADING FOR CCF
C
      WRITE (IPRT, 1020)
      WRITE (IPRT, 1040) ((J,K, K=1,M), J=1,M)
C
C     PRINT CROSS CORRELATIONS
C
      LAG = 0
      I = 0
      DO 30 J = 1, M
         DO 20 K = 1, M
            I = I + 1
            CCF(I) = CCOV(1,J,K) / SQRT(CCOV(1,J,J)*CCOV(1,K,K))
   20    CONTINUE
   30 CONTINUE
      WRITE (IPRT, 1060) LAG, (CCF(L), L=1,I)
C
      DO 60 LAG = 1, LAGMAX
         I = 0
         DO 50 J = 1, M
            DO 40 K = 1, M
               I = I + 1
               IF (.NOT.MISS) GO TO 35
               CCF(I) = CMISS
               IF (MVCHK(CCOV(LAG+1,J,K),CMISS)) GO TO 40
   35          CCF(I) = CCOV(LAG+1,J,K) / SQRT(CCOV(1,J,J)*CCOV(1,K,K))
   40       CONTINUE
   50    CONTINUE
         WRITE (IPRT, 1060) LAG, (CCF(L1), L1=1,I)
   60 CONTINUE
C
C     CHECK FOR MISSING VALUES
C
      IF (.NOT.MISS) RETURN
C
C     PRINT HEADING FOR NUMBERS OF LAGGED PRODUCT PAIRS
C
      WRITE (IPRT, 1030)
      WRITE (IPRT, 1040) ((J,K, K=1,M), J=1,M)
C
C     PRINT NUMBERS OF LAGGED PRODUCT PAIRS FOR EACH CCVF
C
      LAG = 0
      WRITE (IPRT, 1070) LAG, ((NLPPC(1,J,K), K=1,M), J=1,M)
      DO 70 LAG = 1, LAGMAX
         WRITE (IPRT, 1070) LAG, ((NLPPC(LAG+1,J,K), K=1,M), J=1,M)
   70 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (//8H IERR = , I5)
 1010 FORMAT (// 6X, 6H  CCVF)
 1020 FORMAT (// 6X, 6H   CCF)
 1030 FORMAT (// 6X, 6H NLPPC)
 1040 FORMAT (1X, 3HLAG, 16(5X, I1, ',', I1))
 1060 FORMAT (1X, I3, 16F8.4)
 1070 FORMAT (1X, I3, 16I8)
      END

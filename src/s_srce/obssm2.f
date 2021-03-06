*OBSSM2
      SUBROUTINE OBSSM2(N, Y, PVT, SDPVT, RES, SDREST, IFIRST, ILAST)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBOUTINE LISTS THE DATA SUMMARY FOR THE ARIMA ESTIMATION
C     SUBROUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IFIRST,ILAST,N
C
C  ARRAY ARGUMENTS
      REAL
     +   PVT(N),RES(N),SDPVT(N),SDREST(N),Y(N)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   FPLM
      INTEGER
     +   I,IPRT
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIRST, ILAST
C        THE FIRST AND LAST INDICES TO BE LISTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     REAL PVT(N)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES.
C     REAL RES(N)
C        THE RESIDUALS FROM THE FIT.
C     REAL SDPVT(N)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     REAL SDREST(N)
C        THE STANDARDIZED RESIDUALS.
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C
      FPLM = R1MACH(2)
C
      CALL IPRINT(IPRT)
C
      DO 140 I=IFIRST,ILAST
C
C     PRINT DATA SUMMARY.
C
         IF ((SDPVT(I).NE.FPLM) .AND. (SDREST(I).NE.FPLM))
     +      WRITE (IPRT, 1060) I, Y(I), PVT(I), SDPVT(I), RES(I),
     +      SDREST(I)
         IF ((SDPVT(I).NE.FPLM) .AND. (SDREST(I).EQ.FPLM))
     +      WRITE (IPRT, 1050) I, Y(I), PVT(I), SDPVT(I), RES(I)
         IF ((SDPVT(I).EQ.FPLM) .AND. (SDREST(I).EQ.FPLM))
     +      WRITE (IPRT, 1080) I, Y(I), PVT(I), RES(I)
C
  140 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1050 FORMAT (1X, I4, 4E16.8, 4X, 4HNC *, 1X, E9.3)
 1060 FORMAT (1X, I4, 4E16.8, 1X, F7.2, 1X, E9.3)
 1080 FORMAT (1X, I4, 2E16.8, 8X, 4HNC *, 4X, E16.8, 4X, 4HNC *,
     +   1X, E9.3)
      END

*AIMFS
      SUBROUTINE AIMFS(Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK,
     +   NFCST, NFCSTO, IFCSTO, NPRT, FCST, IFCST, FCSTSD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE SUBROUTINE FOR ARIMA ESTIMATION
C     (CONTROL CALL).
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
     +   IFCST,LDSTAK,N,NFAC,NFCST,NFCSTO,NPAR,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   FCST(*),FCSTSD(*),PAR(*),Y(*)
      INTEGER
     +   IFCSTO(*),MSPEC(4,*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT,NFCSTU
      LOGICAL
     +   SAVE
C
C  LOCAL ARRAYS
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMFCNT,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL FCST(IFCST,NFCSTO)
C        THE STORAGE ARRAY FOR THE FORECASTS.
C     REAL FCSTSD(NFCST)
C        THE STORAGE ARRAY FOR THE STANDARD DEVIATIONS OF THE FORECASTS.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFCST
C        THE FIRST DIMENSION OF THE ARRAY FCST.
C     INTEGER IFCSTO(NFCSTO)
C        THE INDICES OF THE ORIGINS FOR THE FORECASTS.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER MSPEC(4,NFAC)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NFCST
C        THE NUMBER OF FORECASTS.
C     INTEGER NFCSTO
C        THE NUMBER OF THE ORIGINS.
C     INTEGER NFCSTU
C        THE NUMBER OF FORCASTES ACTUALLY USED.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE ROUTINE CALLING THE ERROR CHECKING ROUTINE
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
C        (FALSE).
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'A','I','M','F','S',' '/
C
C     SET VARIOUS PROGRAM PARAMETERS
C
      SAVE = .TRUE.
C
      IF ((NFCST.GE.1) .AND. (NFCST.LE.N)) THEN
         NFCSTU = NFCST
      ELSE
         NFCSTU = (N/10)+1
      END IF
C
      CALL AMFCNT(Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK, NFCSTU,
     +   MAX(1,NFCSTO), IFCSTO, NPRT, FCST, IFCST, FCSTSD, NMSUB, SAVE)
C
      IF (IERR.NE.1) RETURN
C
C     PRINT PROPER CALL SEQUENCE
C
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/' THE CORRECT FORM OF THE CALL STATEMENT IS'//
     +  '       CALL AIMFS (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK,'/
     +  '      +            NFCST, NFCSTO, IFCSTO, NPRT, FCST, IFCST,',
     +  ' FCSTSD)')
      END

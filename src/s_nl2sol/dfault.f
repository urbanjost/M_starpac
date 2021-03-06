*DFAULT
      SUBROUTINE DFAULT(IV, V)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C
C  VARIABLE DECLARATIONS
C
C  ARRAY ARGUMENTS
      REAL
     +   V(45)
      INTEGER
     +   IV(25)
C
C  LOCAL SCALARS
      REAL
     +   MACHEP,MEPCRT,ONE,SQTEPS,THREE
      INTEGER
     +   AFCTOL,COSMIN,COVPRT,COVREQ,D0INIT,DECFAC,DELTA0,DFAC,
     +   DINIT,DLTFDC,DLTFDJ,DTYPE,EPSLON,FUZZ,INCFAC,INITS,JTINIT,
     +   LMAX0,MXFCAL,MXITER,OUTLEV,PARPRT,PHMNFC,PHMXFC,PRUNIT,
     +   RDFCMN,RDFCMX,RFCTOL,RLIMIT,SOLPRT,STATPR,TUNER1,TUNER2,
     +   TUNER3,TUNER4,TUNER5,X0PRT,XCTOL,XFTOL
C
C  EXTERNAL FUNCTIONS
      REAL
     +   RMDCON
      INTEGER
     +   IMDCON
      EXTERNAL RMDCON,IMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  ***  SUPPLY NL2SOL (VERSION 2.2) DEFAULT VALUES TO IV AND V  ***
C
C     INTEGER IV(25)
C     REAL V(45)
C/+
C     REAL MAX
C/
C     EXTERNAL IMDCON, RMDCON
C     INTEGER IMDCON
C     REAL RMDCON
C
C     REAL MACHEP, MEPCRT, ONE, SQTEPS, THREE
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER AFCTOL, COSMIN, COVPRT, COVREQ, DECFAC, DELTA0, DFAC,
C    1        DINIT, DLTFDC, DLTFDJ, DTYPE, D0INIT, EPSLON, FUZZ,
C    2        INCFAC, INITS, JTINIT, LMAX0, MXFCAL, MXITER, OUTLEV,
C    3        PARPRT, PHMNFC, PHMXFC, PRUNIT, RDFCMN, RDFCMX,
C    4        RFCTOL, RLIMIT, SOLPRT, STATPR, TUNER1, TUNER2, TUNER3,
C    5        TUNER4, TUNER5, XCTOL, XFTOL, X0PRT
C
      DATA ONE/1.0E0/, THREE/3.0E0/
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA COVPRT/14/, COVREQ/15/, DTYPE/16/, INITS/25/,
     +     MXFCAL/17/, MXITER/18/, OUTLEV/19/,
     +     PARPRT/20/, PRUNIT/21/, SOLPRT/22/,
     +     STATPR/23/, X0PRT/24/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA AFCTOL/31/, COSMIN/43/, DECFAC/22/,
     +     DELTA0/44/, DFAC/41/, DINIT/38/, DLTFDC/40/,
     +     DLTFDJ/36/, D0INIT/37/, EPSLON/19/, FUZZ/45/,
     +     INCFAC/23/, JTINIT/39/, LMAX0/35/, PHMNFC/20/,
     +     PHMXFC/21/, RDFCMN/24/, RDFCMX/25/,
     +     RFCTOL/32/, RLIMIT/42/, TUNER1/26/,
     +     TUNER2/27/, TUNER3/28/, TUNER4/29/,
     +     TUNER5/30/, XCTOL/33/, XFTOL/34/
C
C-----------------------------------------------------------------------
C
      IV(1) = 12
      IV(COVPRT) = 1
      IV(COVREQ) = 1
      IV(DTYPE) = 1
      IV(INITS) = 0
      IV(MXFCAL) = 200
      IV(MXITER) = 150
      IV(OUTLEV) = 1
      IV(PARPRT) = 1
      IV(PRUNIT) = IMDCON(1)
      IV(SOLPRT) = 1
      IV(STATPR) = 1
      IV(X0PRT) = 1
C
      MACHEP = RMDCON(3)
      V(AFCTOL) = 1.0E-20
      IF (MACHEP .GT. 1.0E-10) V(AFCTOL) = MACHEP**2
      V(COSMIN) = MAX(1.0E-6, 1.0E2 * MACHEP)
      V(DECFAC) = 0.5E0
      SQTEPS = RMDCON(4)
      V(DELTA0) = SQTEPS
      V(DFAC) = 0.6E0
      V(DINIT) = 0.0E0
      MEPCRT = MACHEP ** (ONE/THREE)
      V(DLTFDC) = MEPCRT
      V(DLTFDJ) = SQTEPS
      V(D0INIT) = 1.0E0
      V(EPSLON) = 0.1E0
      V(FUZZ) = 1.5E0
      V(INCFAC) = 2.0E0
      V(JTINIT) = 1.0E-6
      V(LMAX0) = 100.0E0
      V(PHMNFC) = -0.1E0
      V(PHMXFC) = 0.1E0
      V(RDFCMN) = 0.1E0
      V(RDFCMX) = 4.0E0
      V(RFCTOL) = MAX(1.0E-10, MEPCRT**2)
      V(RLIMIT) = RMDCON(5)
      V(TUNER1) = 0.1E0
      V(TUNER2) = 1.0E-4
      V(TUNER3) = 0.75E0
      V(TUNER4) = 0.5E0
      V(TUNER5) = 0.75E0
      V(XCTOL) = SQTEPS
      V(XFTOL) = 1.0E2 * MACHEP
C
      RETURN
C  ***  LAST CARD OF DFAULT FOLLOWS  ***
      END

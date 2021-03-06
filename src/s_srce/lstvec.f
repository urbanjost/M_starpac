*LSTVEC
      SUBROUTINE LSTVEC(N, VEC)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE INDICES AND ELEMENT VALUES
C     OF THE VECTOR VEC.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      REAL
     +   VEC(N)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IMAX,IMIN,INDEX,IPRT,NPERL
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   INPERL
      EXTERNAL INPERL
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER IMAX, IMIN
C        THE LARGEST AND SMALLEST INDEX VALUE TO BE PRINTED ON EACH
C        LINE.
C     INTEGER INDEX
C        THE INDEX VALUE TO BE PRINTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER N
C        THE NUMBER OF VALUES TO BE PRINTED IN THE INPUT VECTOR.
C     INTEGER NPERL
C        THE NUMBER OF VALUES TO BE PRINTED PER LINE.
C     REAL VEC(N)
C        THE VECTOR OF VALUES TO BE PRINTED.
C
      CALL IPRINT(IPRT)
C
      NPERL = INPERL(0)
C
      DO 10 I = 1, N, NPERL
         IMIN = I
         IMAX = MIN(I+NPERL-1, N)
         WRITE(IPRT, 1010) (INDEX, INDEX = IMIN, IMAX)
         WRITE(IPRT, 1020) (VEC(INDEX), INDEX = IMIN, IMAX)
   10 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1010 FORMAT(10X, 5HINDEX, I5, 6I15)
 1020 FORMAT(10X, 5HVALUE, 7(1X, G14.7)/)
C
      END

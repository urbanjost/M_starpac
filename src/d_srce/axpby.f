*AXPBY
      SUBROUTINE AXPBY(N,SA,SX,INCX,SB,SY,INCY,SZ,INCZ)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE IS ADAPTED FROM BLAS SUBROUTINE DAXPY.
C
C     OVERWRITE DOUBLE PRECISION SZ WITH DOUBLE PRECISION SA*SX + SB*SY.
C     FOR I = 0 TO N-1, REPLACE  SZ(LZ+I*INCZ) WITH SA*SX(LX+I*INCX) +
C     SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C     AND LY AND LZ ARE DEFINED IN A SIMILAR WAY USING INCY AND INCZ,
C     RESPECTIVELY.
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
      DOUBLE PRECISION
     +   SA,SB
      INTEGER
     +   INCX,INCY,INCZ,N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   SX(*),SY(*),SZ(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IX,IY,IZ,M,MP1,NS
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MOD
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION SX(N), SY(N), SZ(N)
C
      IF(N.LE.0) RETURN
      IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1) .AND. (INCZ .EQ. 1))
     +   GO TO 20
      IF ((INCX .GE. 2) .AND. (INCX .EQ. INCY) .AND. (INCX .EQ. INCZ))
     +   GO TO 60
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IZ = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      IF(INCZ.LT.0)IZ = (-N+1)*INCZ + 1
      DO 10 I = 1,N
        SZ(IZ) = SA*SX(IX) + SB*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
        IZ = IZ + INCZ
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SZ(I) = SA*SX(I) + SB*SY(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SZ(I) = SA*SX(I) + SB*SY(I)
        SZ(I+1) = SA*SX(I+1) + SB*SY(I+1)
        SZ(I+2) = SA*SX(I+2) + SB*SY(I+2)
        SZ(I+3) = SA*SX(I+3) + SB*SY(I+3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
      DO 70 I=1,NS,INCX
        SZ(I) = SA*SX(I) + SB*SY(I)
   70 CONTINUE
      RETURN
      END

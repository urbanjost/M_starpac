*XERABT
      SUBROUTINE XERABT(MESSG,NMESSG)
C
C     LATEST REVISION  -  JANUARY 24, 1990 (JRD)
C
C     ABSTRACT
C        ***NOTE*** MACHINE DEPENDENT ROUTINE
C        XERABT ABORTS THE EXECUTION OF THE PROGRAM.
C        THE ERROR MESSAGE CAUSING THE ABORT IS GIVEN IN THE CALLING
C        SEQUENCE IN CASE ONE NEEDS IT FOR PRINTING ON A DAYFILE,
C        FOR EXAMPLE.
C
C     DESCRIPTION OF PARAMETERS
C        MESSG AND NMESSG ARE AS IN XERROR, EXCEPT THAT NMESSG MAY
C        BE ZERO, IN WHICH CASE NO MESSAGE IS BEING SUPPLIED.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  7 JUNE 1978
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER NMESSG
C
C  ARRAY ARGUMENTS
      CHARACTER MESSG(1)*4
C
      STOP
      END

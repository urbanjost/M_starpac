# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
# WORK IN PROGRESS
* This is just a collection of files to use as seeds for a starpac
library callable from fpm(Fortran Package Manager). *

    https://urbanjost.github.io/M_starpac/

[details="User's Guide"]
```text
                                  USER'S GUIDE

                                       TO

                                    STARPAC

                THE STANDARDS TIME SERIES AND REGRESSION PACKAGE

                ------------------------------------------------

                 National Institute of Standards and Technology
                  (formerly the National Bureau of Standards)
                         Internal Report NBSIR 86-3448
                            (text revised  3/15/90)
                   (STARPAC sample output revised 03/15/90)



                               Janet R. Donaldson
                                 Peter V. Tryon








                          U.S. Department of Commerce
                  Center for Computing and Applied Mathematics
                 National Institute of Standards and Technology
                          Boulder, Colorado 80303-3328
1











                                  In memory of
                                 Peter V. Tryon
                                  1941 - 1982












































                                      <ii>
1











                                   Disclaimer

                  No warranties, express or implied,  are made
                  by  the  distributors  or   developers  that
                  STARPAC or its constituent parts are free of
                  error.   They should  not be  relied upon as
                  the  sole basis for solving a  problem whose
                  incorrect solution could result in injury to
                  person or  property.   If  the  programs are
                  employed  in  such a  manner, it  is  at the
                  user's own  risk and  the  distributors  and
                  developers  disclaim all liability  for such
                  misuse.

                  Computers  have  been   identified  in  this
                  paper  in  order  to  adequately specify the
                  sample  programs  and  test  results.   Such
                  identification does  not  imply  recommendation
                  or   endorsement   by   the  National
                  Institute  of  Standards and Technology, nor
                  does it imply  that the equipment identified
                  is necessarily  the  best  available for the
                  purpose.
























                                      <iii>
1                                   Preface


      STARPAC, the Standards Time  Series and Regression Package, is  a library
 of  Fortran  subroutines  for  statistical  data  analysis  developed  by  the
 Statistical Engineering Division (SED) of the National Institute  of Standards
 and  Technology  (NIST),  formerly  the  National  Bureau  of Standards (NBS),
 Boulder, Colorado.   Earlier versions of this  library were distributed by the
 SED  under  the  name  STATLIB  [Tryon  and  Donaldson, 1978].   Chapter 1 and
 chapter 9 of this document were previously distributed as  NBS Technical Notes
 1068-1  and  1068-2,   respectively   [Donaldson and Tryon,  1983a and 1983b].
 STARPAC   incorporates   many   changes  to   STATLIB,  including   additional
 statistical techniques, improved algorithms and enhanced portability.

      STARPAC consists of  families of subroutines for nonlinear  least squares
 regression, time series analysis  (in both  time and frequency  domains), line
 printer  graphics,  basic  statistical  analysis,  and  linear  least  squares
 regression.  These subroutines feature:

      * ease of use, alone and with other Fortran subroutine libraries;

      * extensive error handling facilities;

      * comprehensive printed reports;

      * no problem size restrictions other than effective machine size; and

      * portability.

      Notation,  format  and naming  conventions are  constant  throughout  the
 STARPAC  documentation,   allowing  the  documentation  for  each   family  of
 subroutines to  be used alone  or in  conjunction with  the  documentation for
 another.

      STARPAC  is  written  in  ANSI Fortran   77  [American National Standards
 Institute, 1977].   Workspace  and  machine-dependent constants  are  supplied
 using subroutines based  on  the  Bell Laboratories  "Framework for a Portable
 Library"  [Fox et al., 1978a].   We have  also used  subroutines  from LINPACK
 [Dongarra et al., 1979],  from  the  "Basic  Linear  Algebra  Subprograms  for
 Fortran Usage" [Lawson et al., 1979],  from  DATAPAC [Filliben, 1977] and from
 the portable special function subroutines  of Fullerton [1977].   The analyses
 generated  by  several  of  the  subroutine  families  have  been adapted from
 OMNITAB II [Hogben et al., 1971]; users are  directed  to  Peavy et al. [1985]
 for information about OMNITAB 80, the current version of OMNITAB.

      Computer facilities for the STARPAC project have been provided in part by
 the  National Oceanic  and Atmospheric Administration  Mountain Administrative
 Support  Center  Computer  Division,  Boulder,  Colorado,  and  we  gratefully
 acknowledge their support.   The STARPAC  subroutine library is the result  of
 the programming  efforts of  Janet  R.  Donaldson  and John  E.  Koontz,  with
 assistance from Ginger A. Caldwell,  Steven M. Keefer, and Linda  L. Mitchell.
 Valuable  contributions have  also been made  by each  of the  members  of the
 Statistical Engineering Division in Boulder, and from many within  the STARPAC
 user community.   We are grateful for the many valuable  comments that we have
 received on early drafts of  the STARPAC documentation; we wish  especially to




                                      <iv>
1thank Paul T. Boggs, Ginger A. Caldwell, Sally E. Howe, John E.  Koontz, James
 T.  Ringland, Ralph J.  Slutz, and  Dominic F. Vecchia.   Finally, we  wish to
 thank Lorna Buhse for excellent manuscript support.


                                                      Janet R. Donaldson
                                                      Peter V. Tryon (deceased)
                                                      October 1985
                                                      Revised September 1987
                                                      Revised February 1990

















































                                      <v>
1                                    Contents


 Preface
 1.   INTRODUCTION TO USING STARPAC
      A.  Overview of STARPAC and Its Contents
      B.  Documentation Conventions
      C.  A Sample Program
      D.  Using STARPAC
          D.1  The PROGRAM Statement
          D.2  The Dimension Statements
          D.3  The CALL Statements
          D.4  STARPAC Output
          D.5  STARPAC Error Handling
          D.6  Common Programming Errors When Using STARPAC
 2.   LINE PRINTER GRAPHICS
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Details
      F.  Examples
 3.   NORMAL RANDOM NUMBER GENERATION
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments
      E.  Computational Methods
      F.  Example
      G.  Acknowledgments
 4.   HISTOGRAMS
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Printed Output
      F.  Example
      G.  Acknowledgments
 5.   STATISTICAL ANALYSIS OF A UNIVARIATE SAMPLE
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Printed Output
      F.  Example
      G.  Acknowledgments









                                      <vi>
16.   ONE-WAY ANALYSIS OF VARIANCE
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Printed Output
      F.  Example
      G.  Acknowledgments
 7.   CORRELATION ANALYSIS
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Print
      F.  Example
      G.  Acknowledgments
 8.   LINEAR LEAST SQUARES
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  The Linear Least Squares Algorithm
          E.2  Computed Results and Printed Output
      F.  Examples
      G.  Acknowledgments
 9.   NONLINEAR LEAST SQUARES
      A.  Introduction
      B.  Subroutine Descriptions
          B.1  Nonlinear Least Squares Estimation Subroutines
          B.2  Derivative Step Size Selection Subroutines
          B.3  Derivative Checking Subroutines
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
               E.1.a  Nonlinear Least Squares Estimation
               E.1.b  Derivative Step Size Selection
               E.1.c  Derivative Checking
          E.2  Computed Results and Printed Output
               E.2.a  The Nonlinear Least Squares Estimation
                      Subroutines
               E.2.b  The Derivative Step Size Selection Subroutines
               E.2.c  The Derivative Checking Subroutines
      F.  Examples
      G.  Acknowledgments









                                      <vii>
110.  DIGITAL FILTERING
      A.  Introduction
      B.  Subroutine Descriptions
          B.1  Symmetric Linear Filter Subroutines
          B.2  Autoregressive or Difference Linear Filter
               Subroutines
          B.3  Gain and Phase Function Subroutines
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Printed Output
      F.  Examples
      G.  Acknowledgments
 11.  COMPLEX DEMODULATION
      A.  Introduction
      B.  Subroutine Descriptions
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
          E.2  Computed Results and Printed Output
      F.  Example
      G.  Acknowledgments
 12.  CORRELATION AND SPECTRUM ANALYSIS
      A.  Introduction
      B.  Subroutine Descriptions
          B.1  Correlation Analysis
               B.1.a  Univariate Series
               B.1.b  Bivariate Series
          B.2  Spectrum Estimation
               B.2.a  Univariate Series
               B.2.b  Bivariate Series
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
               E.1.a  Correlation Analysis
               E.1.b  Spectrum Analysis
          E.2  Computed Results and Printed Output
               E.2.a  Correlation Analysis
               E.2.b  Spectrum Analysis
      F.  Examples
      G.  Acknowledgments















                                      <viii>
113.  ARIMA MODELING
      A.  Introduction
      B.  Subroutine Descriptions
          B.1  ARIMA Estimation Subroutines
          B.2  ARIMA Forecasting Subroutines
      C.  Subroutine Declaration and CALL Statements
      D.  Dictionary of Subroutine Arguments and COMMON Variables
      E.  Computational Methods
          E.1  Algorithms
               E.1.a  ARIMA Estimation
               E.1.b  ARIMA Forecasting
          E.2  Computed Results and Printed Output
               E.2.a  The ARIMA Estimation Subroutines
               E.2.b  The ARIMA Forecasting Subroutines
      F.  Examples
      G.  Acknowledgments
 Appendix A.  CONTINUITY OF VERTICAL PLOTS ON THE CDC CYBER 840 AND 855
 Appendix B.  WEIGHTED LEAST SQUARES
 Appendix C.  ESTIMATING THE NUMBER OF RELIABLE DIGITS IN THE RESULTS OF A
              FUNCTION
 Appendix D.  LIST OF STARPAC SUBPROGRAM NAMES
              D.1  Subprograms specifically written for STARPAC
              D.2  Subprograms from NL2SOL
              D.3  Subprograms from miscellaneous public domain sources
              D.4  Subprograms from LINPACK and BLAS
              D.5  Subprograms specifying machine dependent constants
 Appendix E.  LIST OF STARPAC LABELED COMMON NAMES
 References.
























                                      <ix>
1                                   STARPAC

                The Standards Time Series and Regression Package


                     Janet R. Donaldson and Peter V. Tryon

                          U.S. Department of Commerce
                  Center for Computing and Applied Mathematics
                 National Institute of Standards and Technology
                          Boulder, Colorado 80303-3328

          STARPAC, the Standards Time Series and Regression Package,  is a
      library  of  Fortran  subroutines  for  statistical   data  analysis
      developed by  the Statistical  Engineering Division of  the National
      Institute  of   Standards  and  Technology  (formerly  the  National
      Bureau of Standards), Boulder, Colorado.   Earlier versions  of this
      library were distributed by  the SED  under the name  STATLIB [Tryon
      and Donaldson, 1978].  STARPAC incorporates many changes to STATLIB,
      including additional statistical techniques, improved algorithms and
      enhanced portability.

          STARPAC emphasizes  the statistical  interpretation  of results,
      and, for  this reason,  comprehensive printed  reports  of auxiliary
      statistical information, often in graphical form,  are automatically
      provided to augment the basic statistical computations  performed by
      each  user-callable STARPAC subroutine.   STARPAC thus  provides the
      best  features of  many stand-alone  statistical  software  programs
      within the flexible environment of a subroutine library.


      Key   words:   data  analysis;  nonlinear  least  squares;  STARPAC;
      statistical computing;  statistical subroutine  library; statistics;
      STATLIB; time series analysis.

                                      <x>
1-----                             CHAPTER 1                              -----

                         INTRODUCTION TO USING STARPAC


 A.  Overview of STARPAC and Its Contents

      STARPAC is a portable library of approximately 150 user-callable ANSI  77
 Fortran  subroutines for  statistical data analysis.   Designed primarily  for
 time series analysis and for nonlinear least squares regression,  STARPAC also
 includes subroutines for normal random number generation, line  printer plots,
 basic statistical analyses and linear least squares.  Emphasis has been placed
 on facilitating  the  interpretation of  statistical analyses,  and,  for this
 reason,  comprehensive printed  reports of auxiliary  statistical information,
 often  in graphical  form,  are automatically  provided to  augment  the basic
 statistical  computations performed by each user-callable  STARPAC subroutine.
 STARPAC thus  provides  the  best  features of  many  stand-alone  statistical
 software programs within the flexible environment of a subroutine library.

      STARPAC is  designed to be easy to  use; in  many situations, only  a few
 lines of elementary Fortran code are required for the users' main programs.  A
 fundamental STARPAC  philosophy  is  to  provide  two  or  more  user-callable
 subroutines  for each method of analysis:   one which minimizes the complexity
 of the CALL statement, automatically producing a comprehensive  printed report
 of the  results; and  one or more  others which  provide user  control  of the
 computations, allow suppression of all or part of the printed  reports, and/or
 provide storage of computed results for further analyses.

      STARPAC was developed and is  maintained by the  Center for Computing and
 Applied  Mathematics  of  the  National Institute of  Standards and Technology
 (NIST), Boulder, Colorado.  Users' comments and suggestions,  which  have  had
 significant impact already,  are highly valued and always welcomed.   Comments
 and suggestions should be directed to:

                          Janet R. Donaldson
                          NIST Center for Computing and Applied Mathematics
                          Mail Code 719
                          325 Broadway
                          Boulder, CO 80303-3328.


 B.  Documentation Conventions

      The documentation  for  the various  STARPAC subroutine  families  uses a
 standard  format description  of  the  information  needed to  use  a  STARPAC
 subroutine, including one or more examples.

      References to chapter sections within the STARPAC documentation refer  to
 the  identified  section  within  the current chapter unless explicitly stated
 otherwise.  Figures are identified by the section in  which  they  occur.  For
 example,  figure  B-1  refers to the first figure in section B of this chapter
 (chapter 1).

      Names of INTEGER  and REAL  STARPAC subroutine  arguments  are consistent
 with the implicit Fortran  convention for specifying variable type.   That is,
 variable  names beginning with I through  N are type INTEGER while  all others
 are type REAL unless  otherwise explicitly typed DOUBLE PRECISION  or COMPLEX.
 All dimensioned variables are explicitly declared in STARPAC  documentation by

                                     <1-1>
1means  of  INTEGER,   REAL,  DOUBLE  PRECISION,  or  COMPLEX   statements,  as
 appropriate.   The  convention  used to  specify the  dimension  statements is
 discussed below in section D.2.

      The precision of the STARPAC library is indicated in the printed  reports
 generated by STARPAC:  an S following the STARPAC version number in the output
 heading indicates the single precision  version  is  being  used,  while  a  D
 indicates the double precision version.  The STARPAC documentation is designed
 for  use with both single and double precision versions.  Subroutine arguments
 which are double precision in both versions  are  declared  DOUBLE  PRECISION;
 similarly,  arguments which are single precision in both versions are declared
 REAL.  Arguments whose precision is  dependent  upon  whether  the  single  or
 double precision version of STARPAC is being used are declared <real>.  If the
 double  precision version of the STARPAC library is being used,  then the user
 should substitute DOUBLE PRECISION for <real>; if the single precision version
 is being used,  then  the  user  should  substitute  REAL  for  <real>.  Other
 precision-dependent features are discussed as they occur.


 C.  A Sample Program

      The  sample  program  shown  below  illustrates   the  use   of   STARPAC
 subroutines.   The code shown  is  portable ANSI 77 Fortran.   Section D below
 uses this example to discuss Fortran programming as it relates to STARPAC.

      The data used in this example are 84 relative humidity measurements taken
 at  Pikes  Peak,  Colorado.  STARPAC  subroutine PP,  documented in chapter 2,
 plots the data versus time-order and STARPAC subroutine  STAT,  documented  in
 chapter 5, prints a comprehensive statistical analysis of the data.






























                                     <1-2>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE STAT AND PP USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y AND X MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(100), X(100)
       DOUBLE PRECISION DSTAK(100)
 C
       COMMON /CSTAK/ DSTAK
       COMMON /ERRCHK/ IERR
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     DEFINE LDSTAK, THE LENGTH OF DSTAK
 C
       LDSTAK = 100
 C
 C     READ NUMBER OF OBSERVATIONS INTO N AND
 C          DATA INTO VECTOR Y
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     CREATE A VECTOR OF ORDER INDICES IN X
 C
       DO 10 I=1,N
         X(I) = I
    10 CONTINUE
 C
 C     PRINT TITLE, PLOT OF DATA AND ERROR INDICATOR
 C
       WRITE (IPRT,102)
       CALL PP (Y, X, N)
       WRITE (IPRT,103) IERR
 C
 C     PRINT TITLE, STATISTICAL ANALYSIS OF DATA AND ERROR INDICATOR
 C
       WRITE (IPRT,102)
       CALL STAT (Y, N, LDSTAK)
       WRITE (IPRT,103) IERR
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (11F7.4)
   102 FORMAT ('1DAVIS-HARRISON PIKES PEAK RELATIVE HUMIDITY DATA')
   103 FORMAT (' IERR = ', I1)
       END


                                     <1-3>
1Data:

    84
 0.6067 0.6087 0.6086 0.6134 0.6108 0.6138 0.6125 0.6122 0.6110 0.6104 0.7213
 0.7078 0.7021 0.7004 0.6981 0.7242 0.7268 0.7418 0.7407 0.7199 0.6225 0.6254
 0.6252 0.6267 0.6218 0.6178 0.6216 0.6192 0.6191 0.6250 0.6188 0.6233 0.6225
 0.6204 0.6207 0.6168 0.6141 0.6291 0.6231 0.6222 0.6252 0.6308 0.6376 0.6330
 0.6303 0.6301 0.6390 0.6423 0.6300 0.6260 0.6292 0.6298 0.6290 0.6262 0.5952
 0.5951 0.6314 0.6440 0.6439 0.6326 0.6392 0.6417 0.6412 0.6530 0.6411 0.6355
 0.6344 0.6623 0.6276 0.6307 0.6354 0.6197 0.6153 0.6340 0.6338 0.6284 0.6162
 0.6252 0.6349 0.6344 0.6361 0.6373 0.6337 0.6383
















































                                     <1-4>
1DAVIS-HARRISON PIKES PEAK RELATIVE HUMIDITY DATA
                                                                                                         STARPAC 2.08S (03/15/90)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .7418 -                     + +                                                                               -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .7271 -                    +                                                                                  -
                I                   +                                                                                   I
                I             +          +                                                                              I
                I                                                                                                       I
                I                                                                                                       I
          .7125 -                                                                                                       -
                I                                                                                                       I
                I              +                                                                                        I
                I                                                                                                       I
                I               + +                                                                                     I
          .6978 -                  +                                                                                    -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .6831 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .6685 -                                                                                                       -
                I                                                                                                       I
                I                                                                                  +                    I
                I                                                                                                       I
                I                                                                                                       I
          .6538 -                                                                             +                         -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                      ++                               I
                I                                                          +               + + +                        I
          .6391 -                                                        +                +                           + -
                I                                                    +                          +     +         + + +   I
                I                                                     +                  +        +       ++     +   +  I
                I                                              +   +   ++   + ++ +   +               +                  I
                I                             +                              +    +                 +       +           I
          .6244 -                          + +       + +        + +                                            +        -
                I                         +    + +       +++     +                                                      I
                I                               +  ++ +                                                 +               I
                I                                           +                                            +    +         I
                I     + +++                                  +                                                          I
          .6098 -  ++  +    ++                                                                                          -
                I +                                                                                                     I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .5951 -                                                                  ++                                   -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                 1.0000    9.3000   17.6000   25.9000   34.2000   42.5000   50.8000   59.1000   67.4000   75.7000   84.0000
 IERR = 0


                                                              <1-5>
1DAVIS-HARRISON PIKES PEAK RELATIVE HUMIDITY DATA
                                                                                                         STARPAC 2.08S (03/15/90)
+STATISTICAL ANALYSIS


     N =    84


     FREQUENCY DISTRIBUTION (1-6)            5    25    35     8     1     0     0     4     4     2


     MEASURES OF LOCATION (2-2)                                  MEASURES OF DISPERSION (2-6)

          UNWEIGHTED MEAN          =  6.3734048E-01                    WTD STANDARD DEVIATION   =  3.2405213E-02
          WEIGHTED MEAN            =  6.3734048E-01                    WEIGHTED S.D. OF MEAN    =  3.5356987E-03
          MEDIAN                   =  6.2915000E-01                    RANGE                    =  1.4670000E-01
          MID-RANGE                =  6.6845000E-01                    MEAN DEVIATION           =  2.1076417E-02
          25 PCT UNWTD TRIMMED MEAN=  6.2885952E-01                    VARIANCE                 =  1.0500979E-03
          25 PCT WTD TRIMMED MEAN  =  6.2885952E-01                    COEF. OF. VAR. (PERCENT) =  5.0844430E+00



                    A TWO-SIDED 95 PCT CONFIDENCE INTERVAL FOR MEAN IS 6.3030811E-01 TO  6.4437284E-01 (2-2)
                    A TWO-SIDED 95 PCT CONFIDENCE INTERVAL FOR S.D. IS 2.8137110E-02 TO  3.8211736E-02 (2-7)



     LINEAR TREND STATISTICS (5-1)                               OTHER STATISTICS

          SLOPE                    = -2.4736661E-04                    MINIMUM                  =  5.9510000E-01
          S.D. OF SLOPE            =  1.4414086E-04                    MAXIMUM                  =  7.4180000E-01
          SLOPE/S.D. OF SLOPE = T  = -1.7161450E+00                    BETA ONE                 =  3.7288258E+00
          PROB EXCEEDING ABS VALUE OF OBS T =  .090                    BETA TWO                 =  5.9283926E+00
                                                                       WTD SUM OF VALUES        =  5.3536600E+01
                                                                       WTD SUM OF SQUARES       =  3.4208200E+01
     TESTS FOR NON-RANDOMNESS                                          WTD SUM OF DEV SQUARED   =  8.7158122E-02
                                                                       STUDENTS T               =  1.8025871E+02
          NO. OF RUNS UP AND DOWN  =   47                              WTD SUM ABSOLUTE VALUES  =  5.3536600E+01
          EXPECTED NO. OF RUNS     =   55.7                            WTD AVE ABSOLUTE VALUES  =  6.3734048E-01
          S.D. OF NO. OF RUNS      =    3.82
          MEAN SQ SUCCESSIVE DIFF  =    3.6382337E-04
          MEAN SQ SUCC DIFF/VAR    =     .346


          DEVIATIONS FROM WTD MEAN

               NO. OF + SIGNS      =   22
               NO. OF - SIGNS      =   62
               NO. OF RUNS         =   14
               EXPECTED NO. OF RUNS=   33.5
               S.D. OF RUNS        =    3.51
               DIFF./S.D. OF RUNS  =   -5.550



 NOTE - ITEMS IN PARENTHESES REFER TO PAGE NUMBER IN NBS HANDBOOK 91 (NATRELLA, 1966)
 IERR = 0


                                                              <1-6>
1D.  Using STARPAC

      The following  subsections provide general information needed  when using
 STARPAC,  including a  discussion  of Fortran  programming as  it  relates  to
 STARPAC usage.   Although only elementary knowledge of Fortran is  required to
 use STARPAC, users may still have to consult with a Fortran text  and/or their
 Computing Center staff when questions arise.


 D.1  The PROGRAM Statement

      The PROGRAM statement is used to name the user's main program.   The name
 EXAMPL is  assigned to the  main program  in this example.   The program  name
 cannot  be  the name  of any  variable  in the  user's main  program  and,  in
 addition, cannot be the name of any other subroutine or function called during
 execution of  the user's code.   Specifically, it  cannot be  the name of  any
 subroutine within STARPAC.  To ensure that the name of a STARPAC subroutine is
 not inadvertently  chosen for  the  name of  the main  program,  users  should
 consult with  the local installer of STARPAC  to obtain a list of  the STARPAC
 subroutine names.


 D.2  The Dimension Statements

      The user's program must include dimension statements to define  the sizes
 and types of the  vectors, matrices  and three-dimensional arrays  required by
 each STARPAC subroutine used; STARPAC  itself has  no inherent upper  limit on
 problem size.

   Within the STARPAC documentation for the subroutine declaration and CALL
 statements, lowercase identifiers in the dimension statements
 represent integer constants which must equal or exceed the value
 of the identically-spelled uppercase argument. For example,
 if the documentation specifies the minimum dimension of a variable
 as <real> XM(n,m), and if the number of observations N is 15,
 and the number of columns of data M is 3, then (assuming the single
 precision version of STARPAC is being used) the minimum array size is
 given by the dimension statement REAL XM(15,3).

      The  exact dimensions  assigned  to some  vectors and  matrices  must  be
 supplied in the CALL statements to some STARPAC subroutines.  For example, the
 argument IXM  is defined  as "the exact  value of  the first dimension  of the
 matrix XM as  declared in the calling program."   Continuing the  example from
 the preceding paragraph, if the  statement REAL XM(20,5) is used  to dimension
 the matrix XM for a particular subroutine, and IXM is an argument in  the CALL
 statement, then IXM must have the value 20 regardless of the value assigned to
 the variable N.

      Many STARPAC subroutines  require a work area for  internal computations.
 This  work area is provided by  the DOUBLE PRECISION vector DSTAK.   The rules
 for defining DSTAK are as follows.

      1. Programs which call  subroutines requiring the work vector  DSTAK must
         include  the  statements

                        DOUBLE PRECISION DSTAK (ldstak)
                        COMMON /CSTAK/ DSTAK

         where ldstak indicates the integer constant used to dimension DSTAK.

                                     <1-7>
1
      2. Since all STARPAC subroutines use the same work vector, the  length of
         DSTAK must equal or exceed the  longest length required by any  of the
         individual STARPAC subroutines called by the user's program.

      3. The length, LDSTAK, of the work vector DSTAK must be specified  in the
         CALL statement of any STARPAC subroutine using DSTAK to enable STARPAC
         to verify that there will be sufficient work area for the problem.

      It is recommended that a variable LDSTAK  be set to the length  of DSTAK,
 and that this variable be used in each CALL statement requiring the  length of
 DSTAK to be specified.   Then, if a future modification to the  user's program
 requires the length of DSTAK to  be changed, the only alterations  required in
 the existing code would be to the DOUBLE PRECISION dimension statement  and to
 the statement which assigns the length of DSTAK to LDSTAK.

      STARPAC  manages  its  work area using subroutines modeled after those in
 ACM Algorithm 528:  Framework for a  Portable  Library  [Fox  et  al.  1978a].
 Although STARPAC and the Framework share the same COMMON for their work areas,
 there  are differences between the STARPAC management subroutines and those of
 the  Framework.   In  particular,   the  STARPAC  management  subroutines 
 reinitialize  DSTAK each  time the user invokes a STARPAC subroutine requiring
 work area, destroying all data previously stored in DSTAK;  the Framework only
 initializes  DSTAK  the  first  time  any  of  its  management subroutines are
 invoked,  preserving work area allocations still in use.  Thus,  users must be
 cautious  when  utilizing  STARPAC  with  other  libraries  which  employ  the
 Framework, such as PORT [Fox et al., 1978b].

      The sample program shown in figure C-1 provides an example of the use  of
 dimensioned  variables  with  STARPAC.   The  REAL  vector  Y,  used  by  both
 subroutines PP and STAT,  contains the 84 relative humidity measurements;  its
 minimum length, N (the number of observations), is 84.  The REAL vector X used
 by  subroutine  PP  contains the corresponding time order indices of the data;
 its minimum length is also 84.  The DOUBLE PRECISION vector DSTAK contains the
 work area needed by STAT for intermediate computations; its minimum length, 49
 in this case,  is defined in section D of chapter  4.  In  this  example,  the
 dimensions of Y,  X,  and DSTAK,  are each 100, exceeding the required minimum
 values.


 D.3  The CALL Statements

      The STARPAC CALL statement arguments provide the interface for specifying
 the data to be used, controlling the computations, and providing space for any
 returned results.   The CALL  statements used  in the example  (fig. C-1)  are
 CALL PP(Y, X, N) and CALL STAT(Y, N, LDSTAK).   Note that scalar arguments may
 be specified either by a variable preset to the desired value, as was  done in
 the  example,  or   by  the  actual  numerical  values.    For  example,  CALL
 PP(Y, X, 84)  and CALL STAT(Y, 84, 100)  could have  been used instead  of the
 forms shown.   We recommend  using variables rather than the  actual numerical
 values in order to simplify future changes in the program.  When variables are
 used, changes need to be made in only  one place; numerical values have  to be
 changed every  place they occur.   The use  of variables can also clarify  the
 meaning of the program.





                                     <1-8>
1D.4  STARPAC Output

      Most STARPAC  subroutines produce extensive printed reports,  freeing the
 user from formatting and  printing all  statistics of interest.   The standard
 output device is used for these reports.   The user has the options of titling
 the reports and changing the output device.

      The first page of the report from each STARPAC subroutine does  not start
 on a new page.  This allows the user to supply titles.  For example,

         WRITE (6, 100)
     100 FORMAT ('1DAVIS-HARRISON PIKES PEAK RELATIVE HUMIDITY DATA')
         CALL PP (Y, X, N)

 will print the title DAVIS-HARRISON PIKES PEAK RELATIVE HUMIDITY DATA  on  the
 top line of a new page,  immediately preceding the output produced by the call
 to subroutine PP.   Users should note that titles more than one line in length
 can cause a  printed report  designed for one page to extend beyond the bottom
 of the page.

      The unit number, IPRT, of the  output device used by STARPAC  is returned
 by STARPAC subroutine IPRINT.   Users can change the output device unit number
 by including with their program  a subroutine IPRINT which will  supersede the
 STARPAC subroutine of the same name.  The subroutine must have the form

         SUBROUTINE IPRINT(IPRT)
         IPRT = u
         RETURN
         END

 where u  is an integer value specifying  the output unit to which  all STARPAC
 output will be written.


 D.5  STARPAC Error Handling

      STARPAC  provides extensive error-checking facilities which  include both
 printed  reports and  a  program-accessible  error flag  variable.   There are
 essentially two types of errors STARPAC can detect.

      The first type  of error involves incorrect problem  specification, i.e.,
 one or more of the input arguments in the subroutine statement has an improper
 value.   For example,  the number of observations, N,  might have an obviously
 meaningless non-positive value.  In the case of improper problem specification
 STARPAC generates a  printed report  identifying the subroutine  involved, the
 error  detected, and the proper form  of the  subroutine CALL statement.   The
 latter  is  provided  because  improper  input  is  often  the  result  of  an
 incorrectly specified subroutine argument list.

      A second type of error can be thought of as a computation error:   either
 the initiated calculation cannot be  completed or the results from  the called
 subroutine are questionable.   For example,  when the least squares model  and
 data are found to be  singular, the desired computations cannot  be completed;
 when one or more of the standardized residuals from a least squares fit cannot
 be computed  because the  standard  deviation of  the residual  is  zero,  the
 results of  the error  estimates  from the  least squares  regression  may  be
 questionable.   If a computation error is detected, STARPAC generates a report
 which identifies the error, and, to aid  the user in determining the  cause of

                                     <1-9>
1the error, summarizes the completed results in a printed report.

      STARPAC error reports cannot be suppressed,  even when the normal  output
 from  the  STARPAC  subroutine  has  been suppressed.  (STARPAC output must be
 directed to a separate output device [see section D.4] when users do not  want
 any STARPAC reports displayed under any conditions.) Because  of  this,  users
 seldom have to consciously handle STARPAC error conditions in their code.

      When proper execution of the user's program depends on knowing whether or
 not an error has been detected, the error flag can be examined from within the
 user's code.  When access to the error flag is desired, the statement

                              COMMON /ERRCHK/ IERR

 must be placed with the Fortran declaration statements in the user's  program.
 Following the execution of a STARPAC subroutine, the variable IERR will be set
 to  zero  if  no errors were detected,  and to a nonzero value otherwise;  the
 value of IERR may indicate the type of error [e.g., see chapter 9,  section D,
 argument  IERR].  If  the  CALL  statement is followed with a statement of the
 form

                             IF (IERR .NE. 0) STOP

 then the program will stop  when an  error is detected.   (In figure  C-1, the
 value  of IERR  is printed following  each CALL  statement to  show  the value
 returned.)


 D.6  Common Programming Errors When Using STARPAC

      STARPAC error-checking procedures catch many programming errors and print
 informative diagnostics  when such  errors are detected.   However, there  are
 some  errors  which  STARPAC cannot  detect.   The  more  common of  these are
 discussed below.

      1. The most common error  involves array dimensions which are  too small.
         Although certain arguments are checked by STARPAC to verify that array
         dimensions  are  adequate, if  incorrect information  is  supplied  to
         STARPAC, or if the dimension of an  array which is not checked  is too
         small, the program  will produce  erroneous results  and/or  will stop
         prematurely.   Users should  check the  dimension statements in  their
         program whenever difficulties are encountered in using STARPAC.

      2. The second most common error involves incorrect CALL  statements, that
         is,  CALL  statements   in  which  the  STARPAC  subroutine   name  is
         misspelled,  the  arguments  are  incorrectly  ordered,  one  or  more
         arguments are omitted,  or the  argument types (INTEGER,  REAL, DOUBLE
         PRECISION, and  COMPLEX) are incorrect.   Users having  problems using
         STARPAC should  carefully check their declaration and  CALL statements
         to verify that they agree with the documentation.

      3. The third most  common error  involves incorrect specification  of the
         work vector DSTAK.   Programs which call STARPAC subroutines requiring
         work area must  include both the DOUBLE PRECISION  statement dimension
         DSTAK and the COMMON /CSTAK/ DSTAK statement.

      4. The final  common error involves user-supplied subroutines  which have
         the same name as a  subroutine in  the STARPAC library.   Users should

                                     <1-10>
1        consult  with the local installer of  STARPAC to obtain a list  of all
         STARPAC subroutine names.  This list can then be used to ensure that a
         STARPAC subroutine name has not been duplicated.

      Users  who have  not found  the  cause of  a problem  after  checking the
 possibilities  mentioned  above should  consult with  their  Computing  Center
 advisers.




















































                                     <1-11>
1-----                             CHAPTER 2                              -----

                             LINE PRINTER GRAPHICS


 A.  Introduction

      STARPAC  contains 36 subroutines  for producing  2 basic  styles  of line
 printer plots.

      The first,  called a page plot, uses  a single 11 x 14 inch page  of line
 printer paper.

      The second, called a vertical plot, is designed for plotting time series.
 The user specifies only the y-axis values since the x-axis values (independent
 variable)  are assumed to be  equally spaced  and ordered consecutively.   The
 independent variable in the resulting plot is oriented vertically and the plot
 will continue over as many pages as necessary to plot one point per line.

      Within  these two  basic  styles the  user has  many  options,  including
 controlling the plot symbol, plotting multivariate values, designating missing
 data, using log scales and specifying plot limits and plot size.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  actual declaration and CALL statements are given in section
 C and the subroutine arguments are defined in section D.  The algorithms  used
 and  output  produced by these subroutines are discussed in section E.  Sample
 programs are shown in section F.


 B.  Subroutine Descriptions

      PP (Page Plot) and VP (Vertical Plot) are the  simplest  of  the  STARPAC
 line printer plot subroutines.  For both, plot limits are automatically set by
 the  range  of  the  data  and other control parameters are set to the default
 values given in section D.  The remaining plotting subroutines  are  named  by
 adding letters to the beginning and/or end of PP and VP.

      * Prefix S (e.g., SPP) indicates  the user  controls the plot  symbol for
        each point.

      * Prefix  M  indicates the  subroutine will  accept  multivariate  y-axis
        values (e.g., MPP).

      * Suffix M subroutines allow data with missing observations to be plotted
        (e.g., VPM).

      * Suffix L indicates log scales can be used (e.g., PPL).

      * Suffix C subroutines allow control of the parameters which  specify the
        plot limits, size, scale, etc. (e.g., VPC).


      The  following  table,  which  indicates  the capabilities of each of the
 STARPAC plotting subroutines,  can  be  used  to  select  from  the  available
 subroutines.  Subroutine  declaration and CALL statements are given in section
 C, listed in the same order as in the table.


                                     <2-1>
1
                          STARPAC Plotting Subroutines


             Plot
            Symbol   Multiple   Page   Vertical   Missing    Log      Control
            Control   Y-Axis    Plot     Plot       Data    Scale   Parameters
   Name       (S)       (M)     (PP)     (VP)       (M)      (L)        (C)

   PP                            *
   PPL                           *                            *
   PPC                           *                            *          *
   PPM                           *                   *
   PPML                          *                   *        *
   PPMC                          *                   *        *          *

  SPP          *                 *
  SPPL         *                 *                            *
  SPPC         *                 *                            *          *
  SPPM         *                 *                   *
  SPPML        *                 *                   *        *
  SPPMC        *                 *                   *        *          *

  MPP                   *        *
  MPPL                  *        *                            *
  MPPC                  *        *                            *          *
  MPPM                  *        *                   *
  MPPML                 *        *                   *        *
  MPPMC                 *        *                   *        *          *

   VP                                     *
   VPL                                    *                   *
   VPC                                    *                   *          *
   VPM                                    *          *
   VPML                                   *          *        *
   VPMC                                   *          *        *          *

  SVP          *                          *
  SVPL         *                          *                   *
  SVPC         *                          *                   *          *
  SVPM         *                          *          *
  SVPML        *                          *          *        *
  SVPMC        *                          *          *        *          *

  MVP                   *                 *
  MVPL                  *                 *                   *
  MVPC                  *                 *                   *          *
  MVPM                  *                 *          *
  MVPML                 *                 *          *        *
  MVPMC                 *                 *          *        *          *









                                     <2-2>
1C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given  in  section  D  and
 section  F,  respectively.  The  conventions  used  to  present  the following
 declaration and CALL statements are given in chapter 1, sections B and D.


                                   Page Plots

 PP:     Print Y versus X  scatterplot; linear axes; default control values and
         axis limits; no missing values allowed

         <real> Y(n), X(n)
         :
         :
         CALL PP (Y, X, N)

                                      ===

 PPL:    Print  Y versus  X scatterplot; log  or linear  axes; default  control
         values and axis limits; no missing values allowed

         <real> Y(n), X(n)
         :
         :
         CALL PPL (Y, X, N, ILOG)

                                      ===

 PPC:    Print  Y  versus X  scatterplot; log  or  linear  axes;  user-supplied
         control values and axis limits; no missing values allowed

         <real> Y(n), X(n), YLB, YUB, XLB, XUB
         :
         :
         CALL PPC (Y, X, N, ILOG,
        +          ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===

 PPM:    Print  Y versus X scatterplot; linear axes; default control values and
         axis limits; missing values allowed

         <real> Y(n), YMISS, X(n), XMISS
         :
         :
         CALL PPM (Y, YMISS, X, XMISS, N)

                                      ===










                                     <2-3>
1PPML:   Print  Y versus  X scatterplot; log  or linear  axes; default  control
         values and axis limits; missing values allowed

         <real> Y(n), YMISS, X(n), XMISS
         :
         :
         CALL PPML (Y, YMISS, X, XMISS, N, ILOG)

                                      ===

 PPMC:   Print  Y  versus  X  scatterplot;  log or  linear  axes; user-supplied
         control values and axis limits; missing values allowed

         <real> Y(n), YMISS, X(n), XMISS, YLB, YUB, XLB, XUB
         :
         :
         CALL PPMC (Y, YMISS, X, XMISS, N, ILOG,
        +           ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===

 SPP:    Print  Y versus X  scatterplot with individual plot  symbols specified
         by  user; linear  axes;  default control  values and  axis  limits; no
         missing values allowed

         <real> Y(n), X(n)
         INTEGER ISYM(n)
         :
         :
         CALL SPP (Y, X, N, ISYM)

                                      ===

 SPPL:   Print  Y versus X  scatterplot with individual plot  symbols specified
         by user; log or linear  axes; default control values and  axis limits;
         no missing values allowed

         <real> Y(n), X(n)
         INTEGER ISYM(n)
         :
         :
         CALL SPPL (Y, X, N, ISYM, ILOG)

                                      ===

 SPPC:   Print  Y versus X  scatterplot with individual plot  symbols specified
         by  user; log or  linear axes;  user-supplied control values  and axis
         limits; no missing values allowed

         <real> Y(n), X(n), YLB, YUB, XLB, XUB
         INTEGER ISYM(n)
         :
         :
         CALL SPPC (Y, X, N, ISYM, ILOG,
        +           ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===


                                     <2-4>
1
 SPPM:   Print  Y versus X  scatterplot with individual plot  symbols specified
         by user; linear axes; default control values and axis  limits; missing
         values allowed

         <real> Y(n), YMISS, X(n), XMISS
         INTEGER ISYM(n)
         :
         :
         CALL SPPM (Y, YMISS, X, XMISS, N, ISYM)

                                      ===

 SPPML:  Print  Y versus X  scatterplot with individual plot  symbols specified
         by user; log or linear  axis; default control values and  axis limits;
         missing values allowed

         <real> Y(n), YMISS, X(n), XMISS
         INTEGER ISYM(n)
         :
         :
         CALL SPPML (Y, YMISS, X, XMISS, N, ISYM, ILOG)

                                      ===

 SPPMC:  Print  Y versus X  scatterplot with individual plot  symbols specified
         by  user; log or  linear axes;  user-supplied control values  and axis
         limits; missing values allowed

         <real> Y(n), YMISS, X(n), XMISS, YLB, YUB, XLB, XUB
         INTEGER ISYM(n)
         :
         :
         CALL SPPMC (Y, YMISS, X, XMISS, N, ISYM, ILOG,
        +            ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===

 MPP:    Print   plot  of multiple Y vectors versus  a common  X vector; linear
         axes; default  control  values  and  axis limits;  no  missing  values
         allowed

         <real> YM(n,m), X(n)
         :
         :
         CALL MPP (YM, X, N, M, IYM)

                                      ===











                                     <2-5>
1MPPL:   Print  plot of  multiple Y  vectors versus a  common X  vector; log or
         linear axes; default control values and axis limits; no missing values
         allowed

         <real> YM(n,m), X(n)
         :
         :
         CALL MPPL (YM, X, N, M, IYM, ILOG)

                                      ===

 MPPC:   Print  plot of  multiple Y  vectors versus a  common X  vector; log or
         linear axes; user-supplied control values and axis limits;  no missing
         values allowed

         <real> YM(n,m), X(n), YLB, YUB, XLB, XUB
         :
         :
         CALL MPPC (YM, X, N, M, IYM, ILOG,
        +           ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===

 MPPM:   Print  plot of  multiple Y  vectors versus a  common X  vector; linear
         axes; default control values and axis limits; missing values allowed

         <real> YM(n,m), YMMISS(m), X(n), XMISS
         :
         :
         CALL MPPM (YM, YMMISS, X, XMISS, N, M, IYM)

                                      ===

 MPPML:  Print  plot of  multiple Y  vectors versus a  common X  vector; log or
         linear axes; default  control values  and axis limits;  missing values
         allowed

         <real> YM(n,m), YMMISS(m), X(n), XMISS
         :
         :
         CALL MPPML (YM, YMMISS, X, XMISS, N, M, IYM, ILOG)

                                      ===

 MPPMC:  Print  plot of  multiple Y  vectors versus a  common X  vector; log or
         linear axes;  user-supplied control  values and  axis  limits; missing
         values allowed

         <real> YM(n,m), YMMISS(m), X(n), XMISS, YLB, YUB, XLB, XUB
         :
         :
         CALL MPPMC (YM, YMMISS, X, XMISS, N, M, IYM, ILOG,
        +            ISIZE, NOUT, YLB, YUB, XLB, XUB)

                                      ===




                                     <2-6>
1                                Vertical Plots

 VP:     Print  vertical plot  of Y versus input order;  linear  axes;  default
         control values and axis limits; no missing values allowed

         <real> Y(n)
         :
         :
         CALL VP (Y, N, NS)

                                      ===

 VPL:    Print  vertical plot of Y versus input order; log or linear horizontal
         (Y) axis; default control  values and  axis limits; no  missing values
         allowed

         <real> Y(n)
         :
         :
         CALL VPL (Y, N, NS, ILOG)

                                      ===

 VPC:    Print  vertical plot of Y versus input order; log or linear horizontal
         (Y) axis; user-supplied  control values  and axis  limits;  no missing
         values allowed

         <real> Y(n), YLB, YUB, XLB, XINC
         :
         :
         CALL VPC (Y, N, NS, ILOG,
        +          ISIZE, IRLIN, IBAR, YLB, YUB, XLB, XINC)

                                      ===

 VPM:    Print  vertical plot  of Y  versus input order;  linear axis;  default
         control values and axis limits; missing values allowed

         <real> Y(n), YMISS
         :
         :
         CALL VPM (Y, YMISS, N, NS)

                                      ===

 VPML:   Print  vertical plot of Y versus input order; log or linear horizontal
         (Y) axis;  default  control values  and axis  limits;  missing  values
         allowed

         <real> Y(n), YMISS
         :
         :
         CALL VPML (Y, YMISS, N, NS, ILOG)

                                      ===




                                     <2-7>
1VPMC:   Print  vertical plot of Y versus input order; log or linear horizontal
         (Y) axis; user-supplied control values and axis limits; missing values
         allowed

         <real> Y(n), YMISS, YLB, YUB, XLB, XINC
         :
         :
         CALL VPMC (Y, YMISS, N, NS, ILOG,
        +           ISIZE, IRLIN, IBAR, YLB, YUB, XLB, XINC)

                                      ===

 SVP:    Print  vertical  plot  of Y  versus input order  with individual  plot
         symbols specified by  user; linear  axis; default  control  values and
         axis limits; no missing values allowed

         <real> Y(n)
         INTEGER ISYM(n)
         :
         :
         CALL SVP (Y, N, NS, ISYM)

                                      ===

 SVPL:   Print  vertical  plot  of Y  versus input order  with individual  plot
         symbols specified by user; log or linear horizontal (Y)  axis; default
         control values and axis limits; no missing values allowed

         <real> Y(n)
         INTEGER ISYM(n)
         :
         :
         CALL SVPL (Y, N, NS, ISYM, ILOG)

                                      ===

 SVPC:   Print  vertical  plot  of Y  versus input order  with individual  plot
         symbols  specified  by  user;  log  or  linear  horizontal  (Y)  axis;
         user-supplied  control  values  and  axis limits;  no  missing  values
         allowed

         <real> Y(n), YLB, YUB, XLB, XINC
         INTEGER ISYM(n)
         :
         :
         CALL SVPC (Y, N, NS, ISYM, ILOG,
        +           ISIZE, IREFLN, IBAR, YLB, YUB, XLB, XINC)

                                      ===










                                     <2-8>
1SVPM:   Print  vertical plot  of Y  versus  input order  with individual  plot
         symbols specified by  user; linear  axis; default  control  values and
         axis limits; missing values allowed

         <real> Y(n), YMISS
         INTEGER ISYM(n)
         :
         :
         CALL SVPM (Y, YMISS, N, NS, ISYM)

                                      ===

 SVPML:  Print  vertical plot  of Y  versus  input order  with individual  plot
         symbols specified by user; log or linear horizontal (Y)  axis; default
         control values and axis limits; missing values allowed

         <real> Y(n), YMISS
         INTEGER ISYM(n)
         :
         :
         CALL SVPML (Y, YMISS, N, NS, ISYM, ILOG)

                                      ===

 SVPMC:  Print  vertical plot  of Y  versus  input order  with individual  plot
         symbols  specified  by  user;  log  or  linear  horizontal  (Y)  axis;
         user-supplied control values and axis limits; missing values allowed

         <real> Y(n), YMISS, YLB, YUB, XLB, XINC
         INTEGER ISYM(n)
         :
         :
         CALL SVPMC (Y, YMISS, N, NS, ISYM, ILOG,
        +            ISIZE, IRLIN, IBAR, YLB, YUB, XLB, XINC)

                                      ===

 MVP:    Print  vertical plot  of multiple Y vectors versus input order; linear
         axis; default  control  values  and  horizontal (Y)  axis  limits;  no
         missing values allowed

         <real> YM(n,m)
         :
         :
         CALL MVP (YM, N, M, IYM, NS)

                                      ===












                                     <2-9>
1MVPL:   Print  vertical  plot of multiple Y vectors versus input order; log or
         linear horizontal (Y) axis; default control values and axis limits; no
         missing values allowed

         <real> YM(n,m)
         :
         :
         CALL MVPL (YM, N, M, IYM, NS, ILOG)

                                      ===

 MVPC:   Print  vertical  plot of multiple Y vectors versus input order; log or
         linear  horizontal  (Y) axis;  user-supplied control  values  and axis
         limits; no missing values allowed

         <real> YM(n,m), YLB, YUB, XLB, XINC
         :
         :
         CALL MVPC (YM, N, M, IYM, NS, ILOG,
        +           ISIZE, YLB, YUB, XLB, XINC)

                                      ===

 MVPM:   Print  vertical  plot of multiple Y vectors versus input order; linear
         axis; default control values and axis limits; missing values allowed

         <real> YM(n,m), YMMISS(m)
         :
         :
         CALL MVPM (YM, YMMISS, N, M, IYM, NS)

                                      ===

 MVPML:  Print  vertical  plot of multiple Y vectors versus input order; log or
         linear horizontal (Y)  axis; default  control values and  axis limits;
         missing values allowed

         <real> YM(n,m), YMMISS(m)
         :
         :
         CALL MVPML (YM, YMMISS, N, M, IYM, NS, ILOG)

                                      ===

 MVPMC:  Print  vertical  plot of multiple Y vectors versus input order; log or
         linear  horizontal  (Y) axis;  user-supplied control  values  and axis
         limits; missing values allowed

         <real> YM(n,m), YMMISS(m), YLB, YUB, XLB, XINC
         :
         :
         CALL MVPMC (YM, YMMISS, N, M, IYM, NS, ILOG,
        +            ISIZE, YLB, YUB, XLB, XINC)

                                      ===




                                     <2-10>
1D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates that the  argument is input  to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 IBAR    --> The indicator variable used to designate  whether a  vertical plot
             is to be a bar plot or not.  Bar plots connect the plotted  points
             to  the  reference line [see argument IRLIN] with a string of plot
             symbols,  as is done for example in the  correlation  plots.  [See
             chapter 12.] If IBAR >= 1,  the plot is a bar plot.  If IBAR <= 0,
             it is not.  The default value is IBAR = 0.  When IBAR  is  not  an
             argument  in  the  subroutine  CALL statement the default value is
             used.

 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.] Note that using (or not using) the error  flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected and that  the plot
                      was completed satisfactorily.

             IERR = 1 indicates that  improper input was detected or  that some
                      error prevented the plot from being completed.

 ILOG    --> The indicator variable used  to designate whether the axes  are to
             be on  a log  or linear scale.   ILOG is  a two-digit integer, pq,
             where the value of p is used to designate the scale of  the x-axis
             and the value of q is used  to designate the scale of  the y-axis.
             If p = 0  (q = 0) the x-axis  (y-axis) is  on a  linear  scale; if
             p <> 0 (q <> 0) the  x-axis  (y-axis)  is  on  a  log  scale.  For
             vertical plots, the value of q is used to specify the scale on the
             horizontal-axis  and the value of p is ignored.  The default value
             is ILOG = 0, corresponding to linear scale for both the x-axis and
             the y-axis.  When ILOG is not an argument in the  subroutine  CALL
             statement the default value is used.

 IRLIN   --> The indicator  variable  used to  designate  whether  zero  or the
             series mean is to be plotted as a reference line on  the  vertical
             plots  or  whether no  reference line should be used.  If IRLIN <=
             -1,  no reference line is plotted.  If IRLIN = 0, a reference line
             is plotted showing the location of zero on the plot.  If IRLIN  >=
             1,  a  reference  line  is  plotted  showing the series mean.  The
             default value is IRLIN = -1.  When IRLIN is not an argument in the
             subroutine CALL statement the default value is used.

 ISIZE   --> The indicator variable used to designate the size of a  page plot.
             ISIZE is a two-digit integer, pq, where the value of p is  used to
             designate the  size of the x-axis and  the value  of q is  used to


                                     <2-11>
1            designate  the size  of the y-axis.   If p = 0  (q = 0) the x-axis
             (y-axis) is the maximum possible, 101 columns (51 rows), i.e., 101
             (51) plot positions.  If p <> 0 (q <> 0) the  x-axis  (y-axis)  is
             half  the  maximum,  or 51 columns (26 rows).  For vertical plots,
             the value of q is used to specify the size of the  horizontal-axis
             and  the  value  of p is ignored.  The default value is ISIZE = 0,
             corresponding to a plot of 51 rows by 101 columns.  When ISIZE  is
             not an argument in the subroutine CALL statement the default value
             is used.

 ISYM    --> The  vector of  dimension  at least  N that  contains  the  values
             designating the plotting symbol to  be used  for each point.   The
             plot symbols designated  by each possible integer value  are given
             below.

             ISYM()<= 1 -> +  ISYM() = 9 -> E  ISYM() =16 -> L  ISYM() =23 -> S
             ISYM() = 2 -> .  ISYM() =10 -> F  ISYM() =17 -> M  ISYM() =24 -> T
             ISYM() = 3 -> *  ISYM() =11 -> G  ISYM() =18 -> N  ISYM() =25 -> U
             ISYM() = 4 -> -  ISYM() =12 -> H  ISYM() =19 -> 0  ISYM() =26 -> V
             ISYM() = 5 -> A  ISYM() =13 -> I  ISYM() =20 -> P  ISYM() =27 -> W
             ISYM() = 6 -> B  ISYM() =14 -> J  ISYM() =21 -> Q  ISYM() =28 -> Y
             ISYM() = 7 -> C  ISYM() =15 -> K  ISYM() =22 -> R  ISYM()>=29 -> Z
             ISYM() = 8 -> D

 IYM     --> The exact  value of  the  first  dimension  of the  matrix  YM  as
             specified in the calling program.

 M       --> The number of columns of data in YM.

 N       --> The number of observations.

 NOUT    --> The number of points falling  outside the plot limits that  are to
             be listed following the plot.  If NOUT >= 1,  a message giving the
             total number of points falling outside the plot limits and a  list
             of the coordinates of these points (up to a maximum of NOUT or 50,
             whichever  is  smaller)  is printed.  If NOUT = 0,  only a message
             listing the  number  of  points  falling  outside  the  limits  is
             printed.  If  NOUT  <  0,  no  points are listed and no message is
             given.  The default value is  NOUT  =  0.  When  NOUT  is  not  an
             argument  in  the  subroutine  CALL statement the default value is
             used.

 NS      --> The sampling frequency of  the points plotted on a  vertical plot.
             If NS = 1, every point  is plotted; if NS = 2, every  second point
             is plotted;  if NS = 3,  every third point is  plotted,  etc.  The
             default  value is NS = 1.  When NOUT <= 0 or NS is not an argument
             in the subroutine CALL statement the default value is used.

 X       --> The  vector of  dimension  at least  N that  contains  the  x-axis
             values.

 XINC    --> The  increment to  be  used for  labeling the  x-axis  (i.e.,  the
             vertical-axis)  on vertical  plots.   The x-axis  labels are  XLB,
             XLB + NS*XINC,  XLB + 2*NS*XINC,  etc.    The  default   value  is
             XINC = 1.0.   When XINC is not an argument in  the subroutine CALL
             statement the default value is used.



                                     <2-12>
1XLB     --> The lower bound for the x-axis.

             For page plots:

               The default value is the smallest x-axis value within  the range
               of the y-axis values to be plotted.  If XLB >= XUB,  the default
               value is used.

             For vertical plots:

               The default value is 1.0.

             For both page and vertical plots,  when XLB is not an  argument in
             the subroutine CALL statements  the default  value is used.   (The
             plot limits may be adjusted slightly from the user-supplied values
             when the plotting subroutine uses a log scale.)

 XMISS   --> The missing value code used within the vector X to indicate that a
             value is missing.   The user must indicate missing observations by
             putting  the  missing   value  code  in  place  of   each  missing
             observation.   Missing data are not indicated on page plots in any
             way.

 XUB     --> The  upper bound for the x-axis.  The default value is the largest
             x-axis value within the range of the y-axis values to  be plotted.
             If  XLB  >=  XUB  or XUB is not an argument in the subroutine CALL
             statement the default value is  used.  (The  plot  limits  may  be
             adjusted  slightly  from the user-supplied value when the plotting
             subroutines use a log scale.)

 Y       --> The  vector of  dimension  at least  N that  contains  the  y-axis
             values.

 YLB     --> The lower bound for the y-axis.  The default value is the smallest
             y-axis value within the range of the x-axis values to  be plotted.
             If YLB >= YUB or YLB is not an argument  in  the  subroutine  CALL
             statement the default value is  used.  (The  plot  limits  may  be
             adjusted  slightly  from the user-supplied value when the plotting
             subroutines use a log scale.)

 YM      --> The matrix of dimension at least N by M whose columns each contain
             one of the M sets of N observations to be plotted against a common
             X vector.

 YMISS   --> The missing  value code used  within the  vector Y  to  indicate a
             value is missing.   The user must indicate missing observations by
             putting  the  missing   value  code  in  place  of   each  missing
             observation.   Missing data are indicated on vertical plots by the
             word  "MISSING" next to  the right  axis of the  appropriate line.
             Missing data are not indicated on page plots in any way.

 YMMISS  --> The vector of dimension at  least M  that contains the  codes used
             within each of the M columns of YM to indicate a value is missing,
             where the first element of  YMMISS is  the missing value  code for
             the  first column  of YM,  etc.   The user  must indicate  missing
             observations  by  putting the  appropriate missing  value  code in
             place of each missing observation.   Missing data are indicated on


                                     <2-13>
1            vertical plots by the word "MISSING" next to the right axis of the
             appropriate line.  Missing data are not indicated on page plots in
             any way.

 YUB     --> The upper bound for the y-axis.   The default value is the largest
             y-axis value within the range of the x-axis values to  be plotted.
             If  YLB  >=  YUB  or YUB is not an argument in the subroutine CALL
             statement the default value is  used.  (The  plot  limits  may  be
             adjusted  slightly  from the user-supplied value when the plotting
             subroutines use a log scale.)


 E.  Computational Details

      Plotting Symbols.   The plotting symbol used depends on the  type of plot
 and whether or not more than one point falls on a given plot position.  If two
 to  nine points fall on a  single plot position, the integer  corresponding to
 the number of points is  used as the plotting symbol.   When 10 or more values
 fall on a single position the plotting symbol X is used.  This is the only way
 that integers or X are used as plot symbols.

      Subroutines  without an  S or  M  prefix use  the plotting  symbol  +  to
 indicate one point on a single printer position.

      For subroutines with an S prefix the user-supplied vector ISYM of integer
 values is  used to  specify the  plotting  symbol for  each data  point.   The
 Fourier spectrum plot shown in chapter 12 is an example of this option.

      Subroutines  with an M prefix use  a different letter as the  plot symbol
 for each column of the matrix of the dependent variables (y-axis):   A for the
 first, B for the second, ..., Z for columns  25 and beyond, with X  still only
 used to indicate that 10 or more points fell on a single plot position.

      Continuity of Vertical Plots.     Normally,   a    line   printer    will
 automatically  provide margins at the top  and bottom of each page,  causing a
 break in the continuity of a vertical plot or any other output continuing over
 two  or more pages.  However, these automatic page-ejects can be suppressed by
 the user on many systems.   Appendix A gives the control sequence necessary to
 suppress automatic  page-ejects on a Cyber computer.   Users of  other systems
 should  consult  their  Computer  Center  staff  for  any   equivalent  method
 available.


 F.  Examples

      A sample program,  its data and the  resulting  output  for  the  STARPAC
 plotting  routine  MPP is listed at the end of this section.  The program uses
 MPP to display the 12 years of monthly airline data listed on page 531 of  Box
 and Jenkins [1976] versus month.  The year is indicated by the plotting symbol
 (A = 1949, B = 1950, etc.).

      Other examples of STARPAC plots can be found in the output of many of the
 subroutines  discussed  elsewhere.  The  output  from the complex demodulation
 subroutines includes a sample of the simple vertical plot style  (VP)  and  of
 the  vertical plot of multivariate data style (MVPC) (chapter 11);  the output
 from the autocorrelation and cross correlation subroutines  includes  vertical
 plots using the bar plot option style (VPC) (chapter 12);  and the output from


                                     <2-14>
1the Fourier spectrum subroutines (chapter 12) is  produced  using  the  symbol
 plot style (SPPC).

 Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE MPP USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF X AND YM MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL X(20), YM(20,20)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       IYM = 20
 C
 C     READ NUMBER OF OBSERVATIONS AND NUMBER OF COLUMNS OF DATA
 C          X-AXIS VALUES
 C          Y-AXIS VALUES
 C
       READ (5,100) N, M
       READ (5,101) (X(I), I=1,N)
       READ (5,101) ((YM(I,J), I=1,N), J=1,M)
 C
 C     PRINT TITLE AND CALL MPP FOR PLOT
 C
       WRITE (IPRT,102)
       CALL MPP (YM, X, N, M, IYM)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (12F6.1)
   102 FORMAT ('1RESULTS OF STARPAC PLOT SUBROUTINE MPP')
       END















                                     <2-15>
1Data:

    12   12
    1.0   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0
  112.0 118.0 132.0 129.0 121.0 135.0 148.0 148.0 136.0 119.0 104.0 118.0
  115.0 126.0 141.0 135.0 125.0 149.0 170.0 170.0 158.0 133.0 114.0 140.0
  145.0 150.0 178.0 163.0 172.0 178.0 199.0 199.0 184.0 162.0 146.0 166.0
  171.0 180.0 193.0 181.0 183.0 218.0 230.0 242.0 209.0 191.0 172.0 194.0
  196.0 196.0 236.0 235.0 229.0 243.0 264.0 272.0 237.0 211.0 180.0 201.0
  204.0 188.0 235.0 227.0 234.0 264.0 302.0 293.0 259.0 229.0 203.0 229.0
  242.0 233.0 267.0 269.0 270.0 315.0 364.0 347.0 312.0 274.0 237.0 278.0
  284.0 277.0 317.0 313.0 318.0 374.0 413.0 405.0 355.0 306.0 271.0 306.0
  315.0 301.0 356.0 348.0 355.0 422.0 465.0 467.0 404.0 347.0 305.0 336.0
  340.0 318.0 362.0 348.0 363.0 435.0 491.0 505.0 404.0 359.0 310.0 337.0
  360.0 342.0 406.0 396.0 420.0 472.0 548.0 559.0 463.0 407.0 362.0 405.0
  417.0 391.0 419.0 461.0 472.0 535.0 622.0 606.0 508.0 461.0 390.0 432.0











































                                     <2-16>
1RESULTS OF STARPAC PLOT SUBROUTINE MPP
                                                                                                         STARPAC 2.08S (03/15/90)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
       622.0000 -                                                        L                                              -
                I                                                                                                       I
                I                                                                 L                                     I
                I                                                                                                       I
                I                                                                                                       I
       570.2000 -                                                                                                       -
                I                                                                 K                                     I
                I                                                        K                                              I
                I                                              L                                                        I
                I                                                                                                       I
       518.4000 -                                                                                                       -
                I                                                                 J        L                            I
                I                                                                                                       I
                I                                                        J                                              I
                I                                     L        K                                                        I
       466.6000 -                                                        I        I        K                            -
                I                            L                                                      L                   I
                I                                                                                                       I
                I                                              J                                                      L I
                I                                     K        I                                                        I
       414.8000 - L                 L                                    H                                              -
                I                   K                                             H        2        K                 K I
                I          L                 K                                                               L          I
                I                                                                                                       I
                I                                              H                                                        I
       363.0000 - K                 J                 J                  G                          J        K          -
                I                   I        2        I                                    H                            I
                I J        K                                                      G                 I                   I
                I                                                                                                     2 I
                I          J        H                 H                                                                 I
       311.2000 - I                          H                 G                           G                 J          -
                I          I                                             F                          H        I        H I
                I                                                                 F                                     I
                I H        H                                                                                          G I
                I                   G        G        G                           E                 G        H          I
       259.4000 -                                              F         E                 F                            -
                I                                                                                                       I
                I G                 2        E        F        E                  D        E                 G          I
                I          G                 F        E                  D                          F                 F I
                I                                              D                                                        I
       207.6000 - F                                                                        D        E        F          -
                I E        E        D                                    C        C                                   2 I
                I          F                          D                                    C        D                   I
                I          D        C        D        C        C                                             2          I
                I D                          C                           B        B                 C                 C I
       155.8000 -                                                                          B                            -
                I C        C        B                          B         A        A                          C          I
                I                   A        B                 A                           A        B                 B I
                I          B                 A        2                                                                 I
                I 2        A                                                                        A        B        A I
       104.0000 -                                                                                            A          -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                 1.0000    2.1000    3.2000    4.3000    5.4000    6.5000    7.6000    8.7000    9.8000   10.9000   12.0000



                                                              <2-17>
1-----                             CHAPTER 3                              -----

                        NORMAL RANDOM NUMBER GENERATION


 A.  Introduction

      STARPAC  contains  two  subroutines  for generating pseudo-random numbers
 (noise) which obey  a  normal  probability  law  with  mean  mu  and  standard
 deviation  sigma.  Such  random  numbers  are often useful for evaluating data
 analysis procedures or computer programs.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in  section  D.  The  algorithm  used  by
 these subroutines is discussed in section E.  A sample program showing the use
 of these subroutines is given in section F.


 B.  Subroutine Descriptions

      STARPAC subroutine NRAND generates  a vector  of standard (zero  mean and
 unit  standard  deviation) normal  (Gaussian) random  numbers.   There  is  no
 printed output from this subroutine.

      STARPAC subroutine  NRANDC  generates  Gaussian  noise with  mean  mu and
 standard deviation sigma using the transformation

                                z = sigma*y + mu

 where

 y       is a standard normal pseudo-random number;

 mu      is the desired mean (see section D, argument YMEAN); and

 sigma   is the desired standard deviation (see section D, argument SIGMA).

 There is no printed output from NRANDC.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument  definitions  and  a sample program are given in section D and
 section F,  respectively.  The  conventions  used  to  present  the  following
 declaration and CALL statements are given in chapter 1, sections B and D.


 NRAND:   Generate  a vector of normal pseudo-random numbers with zero mean and
          unit standard deviation

          <real> Y(n)
          :
          :
          CALL NRAND (Y, N, ISEED)

                                      ===


                                     <3-1>
1NRANDC:  Generate a vector of normal pseudo-random numbers with mean YMEAN and
          standard deviation SIGMA

          <real> Y(n)
          :
          :
          CALL NRANDC (Y, N, ISEED, YMEAN, SIGMA)

                                      ===


 D.  Dictionary of Subroutine Arguments

 NOTE:   --> indicates that the  argument is input to the  subroutine and  that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section  D.5.]  Note that using (or not using) the error flag will
             not affect the  printed  error  messages  that  are  automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 ISEED   --> The basis for the pseudo-random number generation.  ISEED must lie
             between  0 and 2**31 - 1, inclusive.   If ISEED is not equal to 0,
             ISEED must be odd.

             For ISEED > 0,  use  of  the  same  value  of  ISEED  will  always
                             generate the same data set.

             For ISEED = 0,  a different data set will be produced by each call
                             to NRAND or NRANDC in the user's program, although
                             the numbers generated will not differ from run  to
                             run.

 N       --> The number of random numbers to be generated.

 SIGMA   --> The standard deviation of the generated random numbers.

 Y       <-- The  vector of dimension  at least  N that contains  the generated
             normal pseudo-random numbers.

 YMEAN   --> The mean of the generated random numbers.


 E.  Computational Methods

      The normal pseudo-random number generation procedure is that of Marsaglia

                                     <3-2>
1and Tsang [1984].   The same pseudo-random numbers (to within final  round-off
 error) will  be generated  on all  computers with  at least  32 binary  digits
 available for  representing integers.   The code was written  by Boisvert  and
























































                                     <3-3>
1Kahanar of the NIST Applied and Computational Mathematics Division.


 F.  Example

      The sample program shown below illustrates the  use  of  both  NRAND  and
 NRANDC.  NRAND  is  used to generate a standard normal pseudo-random sample of
 size 50 from a normal population with zero mean and unit  standard  deviation.
 NRANDC  is  then  used  to generate a sequence of normal pseudo-random numbers
 with a mean of 4 and standard deviation 0.5.  The same seed is used  for  both
 NRAND  and  NRANDC.  Therefore,  the values generated by NRANDC are YMEAN plus
 SIGMA times the values generated by NRAND, i.e.,

            YMEAN(I,2) = YMEAN + SIGMA*YMEAN(I,1) for I = 1, ..., N.

 The generated random numbers are displayed using STARPAC plot  subroutine MVP.
 There is no output from NRAND and NRANDC.










































                                     <3-4>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE NRAND AND NRANDC AND DISPLAY RESULTS WITH MVP USING
 C     SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF YM MUST BE CHANGED TO DOUBLE PRECISION IF
 C          DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL YM(100,2)
 C
 C     SET UP OUTPUT FILE
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       IYM = 100
 C
 C     SET THE SEED
 C         THE NUMBER OF VALUES TO BE GENERATED
 C         THE NUMBER OF SETS OF DATA TO BE GENERATED
 C
       ISEED = 531
       N = 50
       M = 2
 C
 C     GENERATE STANDARD NORMAL PSEUDO-RANDOM NUMBERS INTO COLUMN 1 OF YM
 C     AND      NORMAL PSEUDO-RANDOM NUMBERS WITH MEAN 4.0 AND
 C              STANDARD DEVIATION 0.5 INTO COLUMN 2 OF YM
 C
       CALL NRAND (YM(1,1), N, ISEED)
 C
       YMEAN = 4.0
       SIGMA = 0.5
       CALL NRANDC (YM(1,2), N, ISEED, YMEAN, SIGMA)
 C
 C     PRINT TITLE AND CALL MVP TO PLOT RESULTS,
 C     SAMPLING EVERY OBSERVATION
 C
       WRITE (IPRT,100)
       CALL MVP (YM, N, M, IYM, 1)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT ('1RESULTS OF STARPAC NORMAL PSEUDO-RANDOM NUMBER',
      1  ' GENERATION SUBROUTINES',
      2  ' DISPLAYED WITH STARPAC PLOT SUBROUTINE MVP')
       END


 Data:

 NO DATA NEEDED FOR THIS EXAMPLE


                                     <3-5>
1RESULTS OF STARPAC NORMAL PSEUDO-RANDOM NUMBER GENERATION SUBROUTINES DISPLAYED WITH STARPAC PLOT SUBROUTINE MVP
                                                                                                         STARPAC 2.08S (03/15/90)
             -2.4249   -1.6383    -.8518    -.0652     .7213    1.5078    2.2944    3.0809    3.8675    4.6540    5.4406
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                         A                                                     B                     I
  2.0000      I                                                A                                         B          I
  3.0000      I                                                  A                                        B         I
  4.0000      I                     A                                                       B                       I
  5.0000      I                                  A                                                B                 I
  6.0000      I                 A                                                         B                         I
  7.0000      I                                       A                                              B              I
  8.0000      I                                                          A                                    B     I
  9.0000      I                                           A                                            B            I
  10.000      I                                                      A                                      B       I
  11.000      I                  A                                                         B                        I
  12.000      I                             A                                                   B                   I
  13.000      I                        A                                                      B                     I
  14.000      I                                      A                                              B               I
  15.000      I                                           A                                            B            I
  16.000      I                                       A                                              B              I
  17.000      I                                                                   A                                BI
  18.000      I                               A                                                  B                  I
  19.000      I               A                                                          B                          I
  20.000      I                                                      A                                      B       I
  21.000      I                           A                                                    B                    I
  22.000      I                                   A                                                B                I
  23.000      I                        A                                                     B                      I
  24.000      I                                              A                                          B           I
  25.000      IA                                                                 B                                  I
  26.000      I             A                                                           B                           I
  27.000      I                                          A                                            B             I
  28.000      I                                                         A                                     B     I
  29.000      I               A                                                          B                          I
  30.000      I                                      A                                              B               I
  31.000      I                                    A                                               B                I
  32.000      I                                  A                                                B                 I
  33.000      I                          A                                                    B                     I
  34.000      I                                          A                                             B            I
  35.000      I                                    A                                               B                I
  36.000      I                A                                                         B                          I
  37.000      I                                A                                                 B                  I
  38.000      I                                             A                                           B           I
  39.000      I                            A                                                   B                    I
  40.000      I                       A                                                      B                      I
  41.000      I                                          A                                            B             I
  42.000      I                                  A                                                B                 I
  43.000      I                           A                                                    B                    I
  44.000      I                 A                                                         B                         I
  45.000      I                 A                                                         B                         I
  46.000      I                                      A                                              B               I
  47.000      I                            A                                                   B                    I
  48.000      I                                            A                                           B            I
  49.000      I                           A                                                    B                    I
  50.000      I                 A                                                         B                         I





                                                              <3-6>
1G.  Acknowledgments

      The  code used  to  generate the  pseudo-random numbers  was  written  by
 Boisvert  and  Kahanar  of  the  NIST  Applied and  Computational  Mathematics
 Division.






















































                                     <3-7>
1-----                             CHAPTER 4                              -----

                                   HISTOGRAMS

 A.  Introduction

      STARPAC   contains   two  subroutines  for  producing  histograms.   Both
 subroutines  produce a one-page  printout which  includes, in addition  to the
 histogram, a number  of summary statistics (mean, median,  standard deviation,
 cell fractions, etc.) and several tests for normality.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      HIST provides  the  analysis  described  in  section  A  using  a  preset
 procedure for choosing the number of cells.  The lower and upper bounds of the
 histogram are chosen from the range of the observations.

      HISTC provides the same analysis as  HIST but allows the user  to specify
 the number of cells and the  upper and lower cell boundaries.   Statistics are
 based only on the data within the user-supplied bounds.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The conventions used to present the following  declaration  and
 CALL statements are given in chapter 1, sections B and D.


 HIST:   Compute  and print a histogram  and summary statistics, with automatic
         selection of number of cells and cell boundaries

         <real> Y(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL HIST (Y, N, LDSTAK)

                                      ===












                                     <4-1>
1HISTC:  Compute  and  print a  histogram  and  summary  statistics  with  user
         control of number of cells and cell boundaries

         <real> Y(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL HISTC (Y, N, NCELL, YLB, YUB, LDSTAK)

                                      ===


 D.  Dictionary of Subroutine Arguments

 NOTE:   --> indicates that the  argument is input to  the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.]  Note that using (or not using) the error  flag will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision version of  STARPAC is  being  used  P = 0.5,
             otherwise P = 1.0.  [See  chapter 1, section B.]

             For HIST:   LDSTAK >= (17+N)/2 + 26*P

             For HISTC:  LDSTAK >= (17+N)/2 + max(NCELL, 26)*P

 N       --> The number of observations.

 NCELL   --> The number  of cells in the histogram.   If NCELL <= 0 or NCELL is
             not an argument  in the  subroutine CALL statement  the subroutine
             will choose the number of cells.

 Y       --> The  vector of dimension  at least  N that  contains  the observed
             data.


                                     <4-2>
1YLB     --> The  lower  bound for  constructing the  histogram.   The interval
             [YLB, YUB]  is divided into NCELL increments.  If YLB >= YUB,  the
             lower and upper bounds for constructing the histogram will  be the
             minimum and maximum observations.

 YUB     --> The  upper  bound for  constructing the  histogram.   The interval
             [YLB, YUB]  is divided into NCELL increments.  If YLB >= YUB,  the
             lower and upper bounds for constructing the histogram will  be the
             minimum and maximum observations.


 E.  Computational Methods

 E.1  Algorithms

      The code and output for the histogram subroutines are modeled after early
 versions of MINITAB [Ryan et al., 1974].


 E.2  Computed Results and Printed Output

      The output from the histogram family of subroutines includes a summary of
 the input data in addition to the actual histogram.  This summary includes the
 following  information,  where  the  actual  output  headings are given by the
 uppercase  phrases enclosed in angle braces (<...>).  Results which correspond
 to subroutine CALL statements arguments are identified by the argument name in
 uppercase.  In the  formulas,  x(1),  x(2),  ...,  x(k)  denotes  the  ordered
 observations of Y such that YLB <= Y(i) <= YUB,  i = 1,  ..., N, i.e., x(1) is
 the smallest observation of Y such that YLB  <=  Y(i),  x(k)  is  the  largest
 observation  of  Y  such  that  Y(i)  <=  YUB,  etc.  The value of expressions
 enclosed in square brackets,  e.g.,  [(k/2) + 1],  is the largest integer less
 than or equal to the value of the expression.

      * <NUMBER OF OBSERVATIONS>, N

      * <MINIMUM OBSERVATION>

      * <MAXIMUM OBSERVATION>

      * <HISTOGRAM LOWER BOUND>, YLB

      * <HISTOGRAM UPPER BOUND>, YUB

      * <NUMBER OF CELLS>, NCELL

      * <OBSERVATIONS USED>, k, where

          k = the number of observations for which
              YLB <= Y(i) <= YUB, i = 1, ..., N

      * <MIN. OBSERVATION USED>, x(1), where

          x(1) = the smallest observation such that YLB <= Y(i), i = 1, ..., N

      * <MAX. OBSERVATION USED>, x(k), where

          x(k) = the largest observation such that Y(i) <= YUB, i = 1, ..., N


                                     <4-3>
1     * <MEAN VALUE>, xmean, where

                       k
          xmean = 1/k SUM x(i)
                      i=1

      * <MEDIAN VALUE>, xmedian, where

          xmedian = x([(k+1)/2]) if k is odd

          xmedian = 0.5*( x([k/2]) + x([(k/2)+1]) ) if k is even

      * <25 PCT TRIMMED MEAN>, xtrim, where

                                   k-[k/4]
          xtrim = 1/(k-(2*[k/4]))    SUM    x(i)
                                  i=1+[k/4]

      * <STANDARD DEVIATION>, s, where

                             k
          s = sqrt[ 1/(k-1) SUM (x(i)-xmean)**2 ]
                            i=1

      * <MEAN DEV./STD. DEV.>, r, where

                       k
          r = 1/(s*k) SUM |x(i)-xmean|
                      i=1

      * <sqrt(BETA ONE)>, sqrt(beta1), where

                                        k
          beta1 = k/((k-1)**3 * s**6) (SUM (x(i)-xmean)**3)**2
                                       i=1

      * <BETA TWO>, beta2, where,

                                       k
          beta2 = k/((k-1)**2 * s**4) SUM (x(i)-xmean)**4
                                      i=1

      Information provided  for each cell, u = 1, ..., NCELL, of  the histogram
 includes the following.

      * <INTERVAL MID POINT>, cmid, where

          cmid = YLB + (YUB-YLB)/(2*NCELL)

      * <CUM. FRACT.>, the cumulative fraction of the observations which are in
          cells 1 through u

      * <1-CUM.  FRACT.>, the cumulative fraction of the observations which are
          in cells u through NCELL

      * <CELL FRACT.>, the fraction of the observations which are in cell u

      * <NO. OBS.>, the actual number of observations which are in cell u.

                                     <4-4>
1
      The histogram itself displays  the actual number of observations  in each
 cell  when the  largest number of  observations per  cell does not  exceed 50.
 When  the largest  number of observations  per cell  does exceed  50  then the
 histogram displays the cell fraction.


 F.  Example

      The example program below uses HIST to analyze the 39 measurements of the
 velocity of light shown on page 81 of Mandel [1964].
















































                                     <4-5>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE HIST USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION IF
 C          DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(200)
       DOUBLE PRECISION DSTAK(200)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED DATA
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL HIST TO ANALYZE RESULTS
 C
       WRITE (IPRT,102)
       CALL HIST (Y, N, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (13F5.1)
   102 FORMAT ('1RESULTS OF STARPAC HISTOGRAM SUBROUTINE HIST')
       END


 Data:

    39
   0.4  0.6  1.0  1.0  1.0  0.5  0.6  0.7  1.0  0.6  0.2  1.9  0.2
   0.4  0.0 -0.4 -0.3  0.0 -0.4 -0.3  0.1 -0.1  0.2 -0.5  0.3 -0.1
   0.2 -0.2  0.8  0.5  0.6  0.8  0.7  0.7  0.2  0.5  0.7  0.8  1.1









                                     <4-6>
1RESULTS OF STARPAC HISTOGRAM SUBROUTINE HIST
                                                                                                         STARPAC 2.08S (03/15/90)
 HISTOGRAM

 NUMBER OF OBSERVATIONS =              39
 MINIMUM OBSERVATION    = -5.00000000E-01
 MAXIMUM OBSERVATION    =  1.90000000E+00

 HISTOGRAM LOWER BOUND  = -5.00000000E-01
 HISTOGRAM UPPER BOUND  =  1.90000000E+00

 NUMBER OF CELLS        =               9
 OBSERVATIONS USED      =              39           25 PCT TRIMMED MEAN =  4.23809524E-01
 MIN. OBSERVATION USED  = -5.00000000E-01           STANDARD DEVIATION  =  5.06689395E-01
 MAX. OBSERVATION USED  =  1.90000000E+00           MEAN DEV./STD. DEV. =  7.99040249E-01
 MEAN VALUE             =  4.10256410E-01           SQRT(BETA ONE)      =  3.10646614E-01
 MEDIAN VALUE           =  5.00000000E-01           BETA TWO            =  3.33793260E+00


 FOR A NORMAL DISTRIBUTION, THE VALUES (MEAN DEVIATION/STANDARD DEVIATION), SQRT(BETA ONE), AND BETA TWO ARE APPROXIMATELY
 0.8, 0.0 AND 3.0, RESPECTIVELY.  TO TEST THE NULL HYPOTHESIS OF NORMALITY, SEE TABLES OF CRITICAL VALUES PP. 207-208,
 BIOMETRIKA TABLES FOR STATISTICIANS, VOL. 1.  SEE PP. 67-68 FOR A DISCUSSION OF THESE TESTS.



     INTERVAL     CUM.   1-CUM.   CELL   NO.                   NUMBER OF OBSERVATIONS
     MID POINT   FRACT.  FRACT.  FRACT.  OBS.
+                                              0        10        20        30        40        50
    ------------------------------------------  +---------+---------+---------+---------+---------+
   -3.666667E-01   .128   1.000    .128     5    +++++
   -1.000000E-01   .256    .872    .128     5    +++++
    1.666667E-01   .410    .744    .154     6    ++++++
    4.333333E-01   .564    .590    .154     6    ++++++
    7.000000E-01   .846    .436    .282    11    +++++++++++
    9.666667E-01   .949    .154    .103     4    ++++
    1.233333E+00   .974    .051    .026     1    +
    1.500000E+00   .974    .026    .000     0
    1.766667E+00  1.000    .026    .026     1    +





















                                                              <4-7>
1G.  Acknowledgments

      The code and output for the  histogram subroutines is modeled on  that in
 early versions of MINITAB [Ryan et al., 1974].























































                                     <4-8>
1-----                            CHAPTER 5                               -----

                  STATISTICAL ANALYSIS OF A UNIVARIATE SAMPLE


 A.  Introduction

      STARPAC contains 4 subroutines for performing a comprehensive statistical
 analysis of  a univariate sample.   They each  compute 53 different statistics
 which summarize the sample  through measures of location (mean,  median, etc),
 measures  of  dispersion  (standard  deviation,  mean  deviation,   etc.)  and
 diagnostic features  such  as tests  for outliers,  non-normality,  trends and
 non-randomness (assuming  the input order  of the  data is  a  meaningful time
 sequence).  Common statistics such as Student's t and confidence intervals for
 the mean and standard deviation are also included.   NBS Technical Note 756, A
 User's Guide to the OMNITAB Command "STATISTICAL ANALYSIS," by H. H. Ku [1973]
 provides a complete discussion of  the output  of these subroutines,  which is
 the same output as that provided by the OMNITAB II Command STATISTICAL [Hogben
 et al., 1971].

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      STAT computes and prints  the  53  descriptive  statistics  described  in
 section A.

      STATS  provides the same analysis as STAT but allows the user to suppress
 the printed output and store the computed statistics for further use.

      STATW and STATWS perform a weighted analysis and are  otherwise identical
 to STAT and STATS, respectively.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.


 STAT:   Compute and print 53 statistics describing the input data

         <real> Y(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL STAT (Y, N, LDSTAK)

                                      ===



                                     <5-1>
1STATS:  Compute and optionally print 53 statistics describing input data;
         return statistics

         <real> Y(n), STS(53)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL STATS (Y, N, LDSTAK, STS, NPRT)

                                      ===

 STATW:  Compute and print 53 statistics describing weighted input data

         <real> Y(n), WT(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL STATW (Y, WT, N, LDSTAK)

                                      ===

 STATWS: Compute and optionally print 53 statistics describing weighted input
         data; return statistics

         <real> Y(n), WT(n), STS(53)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL STATWS (Y, WT, N, LDSTAK, STS, NPRT)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that  the argument is input to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.






                                     <5-2>
1IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.]  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed (N/2) + 7.

 N       --> The number of observations.

 NPRT    --> The argument controlling printed output.

             If NPRT  = 0 the printed output is suppressed.

             If NPRT <> 0 the printed output is provided.

 STS     --> The  vector of dimension  at least  53 that contains  the computed
             statistics.   The contents  of STS  are listed  below,  along with
             applicable references;  the number in parenthesis is  the location
             within STS that the given statistic is stored.  In the formulas, x
             denotes the ordered observations of Y  for  which  the  weight  is
             nonzero,  i.e.,  x(1)  is  the  smallest  observation  of Y with a
             nonzero weight,  x(k) is the  largest  observation  of  Y  with  a
             nonzero  weight,  etc.  The weight associated with x(i) is denoted
             w(i).   Zero  weighted  observations  are  not  included  in   the
             analysis.  The  value  of expressions enclosed in square brackets,
             e.g.,  [(k/2) + 1],  is the largest integer less than or equal  to
             the value of the expression.

             (1)  NUMBER OF OBSERVATIONS, N

             (2)  NUMBER OF NONZERO WEIGHTS, k

             (3)  UNWEIGHTED MEAN [Dixon and Massey, 1957, p. 14],

                               k
                  xmean = 1/k SUM x(i)
                              i=1

             (4)  WEIGHTED MEAN [Brownlee, 1965, pp. 95-97],

                              k
                             SUM w(i)*x(i)
                             i=1
                  xwtdmean = ---------------
                                 k
                                SUM w(i)
                                i=1

             (5)  MEDIAN [Dixon and Massey, 1957, p. 70],

                  xmedian = x([(k+1)/2]) if k is odd

                  xmedian = 0.5*( x([k/2]) + x([(k/2)+1]) ) if k is even

                                     <5-3>
1
             (6)  MID-RANGE [Dixon and Massey, 1957, p. 71],

                  xmid = 0.5*(x(1) + x(k))

             (7)  25 PCT UNWTD TRIMMED MEAN [Crow and Siddiqui, 1967],

                                         k-[k/4]
                  xtrim = 1/[k-2*[k/4]]   SUM    x(i)
                                        i=1+[k/4]

             (8)  25 PCT WTD TRIMMED MEAN,

                              k-[k/4]
                                SUM    w(i)*x(i)
                             i=1+[k/4]
                  xwtdtrim = -------------------
                                k-[k/4]
                                  SUM    w(i)
                               i=1+[k/4]

             (9)  WEIGHTED STANDARD DEVIATION [Snedecor and Cochran, 1967,
                  p. 44],

                                     k
                  s = sqrt( 1/(k-1) SUM  w(i)*(x(i)-xwtdmean)**2 )
                                    i=1

             (10) WTD S.D. OF MEAN [Brownlee, 1965, p. 80],

                                     k
                  smean = s / sqrt( SUM w(i) )
                                    i=1

             (11) RANGE [Snedecor and Cochran, 1967, p. 39],

                  xrange = x(k) - x(1)

             (12) MEAN DEVIATION [Duncan, 1965, p. 50],

                                  k
                  xmeandev = 1/k SUM |x(i)-xwtdmean|
                                 i=1

             (13) VARIANCE [Snedecor and Cochran, 1967, p. 44],

                                k
                  s2 = 1/(k-1) SUM  w(i)*(x(i)-xwtdmean)**2
                               i=1

             (14) COEF. OF VAR. (PERCENT) [Snedecor and Cochran, 1967, p. 62],

                  cvar = |100*s/xwtdmean|






                                     <5-4>
1            (15) LOWER CONFIDENCE LIMIT, MEAN [Natrella, 1966, pp. 2-2, 2-3],

                  xmean - t(0.025)*smean

                  where  t(0.025) is the  appropriate  t-statistic  with  (k-1)
                  degrees of freedom

             (16) UPPER CONFIDENCE LIMIT, MEAN [Natrella, 1966, pp. 2-2, 2-3],

                  xmean + t(0.025)*smean

                  where  t(0.025) is the  appropriate  t-statistic  with  (k-1)
                  degrees of freedom

             (17) LOWER CONFIDENCE LIMIT, S.D. [Natrella, 1966, p.  2-7],

                  s*sqrt((k-1)/chi2(0.975))

                  where  chi2(0.975)  is  the  appropriate chi-square statistic
                  with (k-1) degrees of freedom

             (18) UPPER CONFIDENCE LIMIT, S.D. [Natrella, 1966, p.  2-7],

                  s*sqrt((k-1)/chi2(0.025))

                  where chi2(0.025) is  the  appropriate  chi-square  statistic
                  with (k-1) degrees of freedom

             (19) SLOPE [Fisher, 1950, p. 136],

                                       k
                  B = 12/(k*(k**2-1)) SUM i*(x(i) - xwtdmean)
                                      i=1

             (20) S.D. OF SLOPE,

                                                    k
                       sqrt( -B**2*k*(k**2-1) + 12 SUM (x(i)-xwtdmean)**2 )
                                                   i=1
                  sB = ----------------------------------------------------
                                  sqrt( k*(k**2-1)*(k-2) )

             (21) SLOPE/S.D. OF SLOPE = T,

                  t0 = B / sB with (k-2) degrees of freedom

             (22) PROB EXCEEDING ABS VALUE OF OBS T [Brownlee, 1965, p. 344],

                  Prob (t < -|t0| and t > +|t0| )

             (23) NO. OF RUNS UP AND DOWN [Brownlee, 1965, p. 223], r

             (24) EXPECTED NO. OF RUNS [Bradley, 1965, p. 279],

                  E(r) = (2k-1)/3




                                     <5-5>
1            (25) S.D. OF NO. OF RUNS [Bradley, 1965, p. 279],

                  rsd = sqrt( [(16*k-29)/90] )

             (26) MEAN SQ SUCCESSIVE DIFF [Brownlee, 1965, p. 222],

                               k-1
                  D2 = 1/(k-1) SUM (x(i+1)-x(i))**2
                               i=1

             (27) MEAN SQ SUCC DIFF/VAR [Brownlee, 1965, p. 222],

                  D2/s2

             (28) NO. OF + SIGNS,

                  u = number of times sign of (x(i)-xwtdmean) is positive

             (29) NO. OF - SIGNS,

                  v = number of times sign of (x(i)-xwtdmean) is negative

             (30) NO. OF RUNS [Brownlee, 1965, p. 224],

                  RUNS = 1 + number of changes in sign of (x(i)-xwtdmean)

             (31) EXPECTED NO. OF RUNS [Brownlee, 1965, p. 227],

                  E(RUNS) = 1 + (2*u*v)/k

             (32) S.D. OF RUNS [Brownlee, 1965, p. 230],

                  RUNSsd =  sqrt( 2*u*v*(2*u*v-u-v)/((u+v)**2*(k-1)) )

             (33) DIFF./S.D. OF RUNS [Brownlee, 1965, p. 230],

                  [RUNS - E(RUNS)] / RUNSsd

             (34) MINIMUM [Natrella, 1966, p. 19-1],

                  x(1) = smallest value with nonzero weight

             (35) MAXIMUM [Natrella, 1966, p. 19-3],

                  x(k) = largest value with nonzero weight

             (36) BETA ONE [Snedecor and Cochran, 1967, p. 86],

                                                 k
                  beta1 = k/((k-1)**3*s**6) * ( SUM (x(i)-xwtdmean)**3 )**2
                                                i=1

             (37) BETA TWO [Snedecor and Cochran, 1967, p. 87],

                                             k
                  beta2 = k/((k-1)**2*s**4) SUM (x(i)-xwtdmean)**4
                                            i=1


                                     <5-6>
1            (38) WTD SUM OF VALUES,

                   k
                  SUM w(i)*x(i)
                  i=1

             (39) WTD SUM OF SQUARES,

                   k
                  SUM w(i)*x(i)**2
                  i=1

             (40) WTD SUM OF DEVS SQUARED,

                   k
                  SUM w(i)*(x(i)-xwtdmean)**2
                  i=1

             (41) STUDENT'S T [Brownlee, 1965, p. 296],

                             k
                  t = sqrt( SUM w(i) )*xwtdmean/s with (k-1) degrees of freedom
                            i=1

             (42) WTD SUM ABSOLUTE VALUES,

                   k
                  SUM w(i)*|xi|
                  i=1

             (43) WTD AVE ABSOLUTE VALUES,

                   k
                  SUM w(i)*|xi|
                  i=1
                  -------------
                      k
                     SUM w(i)
                     i=1

             (44-53) FREQUENCY DISTRIBUTION [Freund and Williams, 1958, p. 17].

 WT      <-- The vector of dimension at least  N that contains the weights.   A
             zero  weight  excludes  the  corresponding  observation  from  the
             analysis.   If the  weights are  all equal to  1.0, the  resulting
             analysis is equivalent to an unweighted analysis.

 Y       --> The vector of dimension at least N that contains the observations.
             The tests for trend  and randomness will not be  meaningful unless
             the sample is ordered with respect to time or some  other relevant
             variable.








                                     <5-7>
1E.  Computational Methods

 E.1.  Algorithms

      Formulas  for  the  computed  statistics  are  given  in  section D under
 argument STS.  The code for the statistical analysis  subroutines  is  adapted
 from OMNITAB II [Hogben et al., 1971].


 E.2.  Computed Results and Printed Output

      The output consists of a one-page display of the 53 statistics listed for
 argument STS in section D.  The argument, NPRT, controlling the printed output
 is discussed in section D.


 F.  Example

      The example program below uses STAT to analyze the 39 measurements of the
 velocity of light shown on page 81 of Mandel [1964].







































                                     <5-8>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE STAT USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION IF
 C          DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(200)
       DOUBLE PRECISION DSTAK(200)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED DATA
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL STAT TO PERFORM STATISTICAL ANALYSIS
 C
       WRITE (IPRT,102)
       CALL STAT (Y, N, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (13F5.1)
   102 FORMAT ('1RESULTS OF STARPAC STATISTICAL ANALYSIS',
      *  ' SUBROUTINE STAT')
       END


 Data:

    39
   0.4  0.6  1.0  1.0  1.0  0.5  0.6  0.7  1.0  0.6  0.2  1.9  0.2
   0.4  0.0 -0.4 -0.3  0.0 -0.4 -0.3  0.1 -0.1  0.2 -0.5  0.3 -0.1
   0.2 -0.2  0.8  0.5  0.6  0.8  0.7  0.7  0.2  0.5  0.7  0.8  1.1








                                     <5-9>
1RESULTS OF STARPAC STATISTICAL ANALYSIS SUBROUTINE STAT
                                                                                                         STARPAC 2.08S (03/15/90)
+STATISTICAL ANALYSIS


     N =    39


     FREQUENCY DISTRIBUTION (1-6)            5     3     8     3     7     7     5     0     0     1


     MEASURES OF LOCATION (2-2)                                  MEASURES OF DISPERSION (2-6)

          UNWEIGHTED MEAN          =  4.1025641E-01                    WTD STANDARD DEVIATION   =  5.0668940E-01
          WEIGHTED MEAN            =  4.1025641E-01                    WEIGHTED S.D. OF MEAN    =  8.1135237E-02
          MEDIAN                   =  5.0000000E-01                    RANGE                    =  2.4000000E+00
          MID-RANGE                =  7.0000000E-01                    MEAN DEVIATION           =  4.0486522E-01
          25 PCT UNWTD TRIMMED MEAN=  4.2380952E-01                    VARIANCE                 =  2.5673414E-01
          25 PCT WTD TRIMMED MEAN  =  4.2380952E-01                    COEF. OF. VAR. (PERCENT) =  1.2350554E+02



                    A TWO-SIDED 95 PCT CONFIDENCE INTERVAL FOR MEAN IS 2.4600615E-01 TO  5.7450667E-01 (2-2)
                    A TWO-SIDED 95 PCT CONFIDENCE INTERVAL FOR S.D. IS 4.1408928E-01 TO  6.5301023E-01 (2-7)



     LINEAR TREND STATISTICS (5-1)                               OTHER STATISTICS

          SLOPE                    = -4.0080972E-03                    MINIMUM                  = -5.0000000E-01
          S.D. OF SLOPE            =  7.2760495E-03                    MAXIMUM                  =  1.9000000E+00
          SLOPE/S.D. OF SLOPE = T  = -5.5086172E-01                    BETA ONE                 =  9.6501319E-02
          PROB EXCEEDING ABS VALUE OF OBS T =  .585                    BETA TWO                 =  3.3379326E+00
                                                                       WTD SUM OF VALUES        =  1.6000000E+01
                                                                       WTD SUM OF SQUARES       =  1.6320000E+01
     TESTS FOR NON-RANDOMNESS                                          WTD SUM OF DEV SQUARED   =  9.7558974E+00
                                                                       STUDENTS T               =  5.0564517E+00
          NO. OF RUNS UP AND DOWN  =   23                              WTD SUM ABSOLUTE VALUES  =  2.0600000E+01
          EXPECTED NO. OF RUNS     =   25.7                            WTD AVE ABSOLUTE VALUES  =  5.2820513E-01
          S.D. OF NO. OF RUNS      =    2.57
          MEAN SQ SUCCESSIVE DIFF  =    2.8289474E-01
          MEAN SQ SUCC DIFF/VAR    =    1.102


          DEVIATIONS FROM WTD MEAN

               NO. OF + SIGNS      =   20
               NO. OF - SIGNS      =   19
               NO. OF RUNS         =    8
               EXPECTED NO. OF RUNS=   20.5
               S.D. OF RUNS        =    3.08
               DIFF./S.D. OF RUNS  =   -4.056



 NOTE - ITEMS IN PARENTHESES REFER TO PAGE NUMBER IN NBS HANDBOOK 91 (NATRELLA, 1966)



                                                              <5-10>
1G.  Acknowledgments

      The code and output  for the statistical analysis subroutines  is adapted
 from OMNITAB II [Hogben et al., 1971].























































                                     <5-11>
1-----                            CHAPTER 6                               -----

                          ONE-WAY ANALYSIS OF VARIANCE


 A.  Introduction

      STARPAC contains two subroutines  for one-way analysis of variance.   The
 output from these subroutines  includes the  usual analysis of  variance table
 plus  the robust Kruskal-Wallis rank test.   Comprehensive summary  statistics
 are also given,  including means, standard deviations, standard  deviations of
 the  mean and confidence intervals for  the mean of each group.   Within group
 standard deviations, standard deviations of the mean and 95-percent confidence
 intervals for  the mean are  given assuming  fixed effects and  random effects
 models  and also assuming ungrouped data.   The output also includes pair-wise
 multiple  comparisons  using  the  Newman-Keuls and Scheffe'  techniques;  the
 Cochran's C, the Bartlett-Box F  and variance  ratio tests for  homogeneity of
 variances; and the random effects model components of variance estimate.   The
 analysis performed by these subroutines is  the same as that performed  by the
 OMNITAB II command  ONEWAY  [Hogben  et al., 1971].   A reference  for one-way
 analysis of variance is Brownlee [1965], chapter 10.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      Subroutine AOV1  computes  and prints  the one-way  analysis  of variance
 described in section A.

      Subroutine AOV1S provides the same  analysis as AOV1 but allows  the user
 to suppress the printed output and to store the number of observations in each
 group, the group means and the group standard deviations.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.


 AOV1:   Compute and print a one-way analysis of variance of the input data

         <real> Y(n), TAG(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL AOV1 (Y, TAG, N, LDSTAK)

                                      ===



                                     <6-1>
1AOV1S:  Compute and  optionally  print a  one-way  analysis of variance of the
         input data; return tag value of each group,  number of observations in
         each group, group averages, and group standard deviations

         <real> Y(n), TAG(n), GSTAT(igstat,4)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL AOV1S (Y, TAG, N, LDSTAK, NPRT, GSTAT, IGSTAT, NG)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates that the  argument is  input to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 GSTAT   <-- The matrix of dimension at least NG by four whose columns contain,
             in order, the tag value  of the  group, number of  observations in
             the group, group average and group standard deviation.  The groups
             are in order of ascending tag values.

 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.]  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 IGSTAT  --> The exact  value of the  first dimension  of the  matrix  GSTAT as
             specified in the calling program.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single  precision version of  STARPAC is  being used P  = 0.5,
             otherwise P = 1.0 [see chapter 1, section B].

             For AOV1:   LDSTAK >= 22 + N + (8*NG+N)*P

             For AOV1S:  LDSTAK >= 22 + N + (4*NG+N)*P


                                     <6-2>
1N       --> The total number of observations (all of the groups combined).

 NG      <-- The number of distinct  groups, that  is, the number  of different
             positive tag values.

 NPRT    --> The parameter controlling printed output.

             If NPRT  = 0 the printed output is suppressed.

             If NPRT <> 0 the printed output is provided.

 TAG     --> The vector of dimension at least N that contains the tag value for
             each observation.  The  tag  values  may  be  any  <real>  number.
             Groups  are  formed from observations having the same positive tag
             value.  Observations having a zero or negative tag value will  not
             be included in the analysis.

 Y       --> The  vector of dimension  at least  N that  contains  the observed
             data.   The order  of the  observations in  Y  is arbitrary  since
             groups are specified  by the values in the  corresponding elements
             of vector TAG.


 E.  Computational Methods

 E.1  Algorithms

      The code and output for the one-way analysis of variance  subroutines are
 adapted from OMNITAB II [Hogben et al., 1971].  The computations performed are
 discussed below.


 E.2  Computed Results and Printed Output

      The  argument  controlling  the  printed  output,  NPRT,  is discussed in
 section D.

      Each of the five sections of automatic printing is described below  under
 the  headings  which appear on the printed page as shown in example F-1b.  The
 discussion is taken from the OMNITAB II User's Reference Manual [Hogben et al.
 1971], pages 122 to 124.

      ANALYSIS OF VARIANCE.  The traditional analysis of variance for a one-way
 classification is printed.   This shows the sources of variation, the  degrees
 of freedom, the sums of squares, the mean squares, the F-ratio for testing for
 differences  between group means  and the  significance level of  the F-ratio.
 The usual  assumptions  of normality,  independence and  constant  variance of
 measurement errors are made.   A discussion of the statistical treatment of  a
 one-way classification can be found in section 10.2 of Brownlee [1965].

      If the significance level of the F ratio (F PROB.) for

         F = (Between Groups Mean Square) / (Within Groups Mean Square)

 is less  than 0.10 and the number  of groups  exceeds two, the  between groups
 (means) sum of squares is separated into two components:  the first associated
 with the slope (one degree of freedom) and the second  representing deviations
 about  the straight line regression of  group averages on group number.   This

                                     <6-3>
1information, which does not appear in a traditional analysis of  variance, can
 be used  to  examine  the  effect  of  time.  A  discussion  of  some  of  the
 statistical  aspects  of this procedure are found in section 11.12 of Brownlee
 [1965].

      Following the above mentioned  analysis of variance, the results  for the
 Kruskal-Wallis non-parametric H-test for testing for differences between group
 means  (averages) are  printed.   The value  of H  is printed  along  with its
 significance level (F PROB.).   The H-test uses the ranks of the  measurements
 and avoids  any  assumption  about  the distribution  of  measurement  errors.
 Details of this test may be found in section 7.7 of Brownlee [1965].


      ESTIMATES.   The following  items are printed for each  group:
        (1) group number,
        (2) number of observations in the group,
        (3) mean,
        (4) within standard deviation,
        (5) standard deviation of the mean,
        (6) minimum (i.e., smallest) observation,
        (7) maximum (i.e., largest) observation,
        (8) the sum of the ranks of the observations, and
        (9) a 95-percent confidence interval for the mean.

 The  results are  printed  with  the  group  numbers  (tags)  in  consecutive,
 increasing order regardless of the order in which the numbers were entered.

      In printing  the  means  and  standard  deviations  of  the  groups,  the
 characters  + and - are put immediately after the largest and smallest values.
 If two or more values are tied for the largest value,  the character + is  put
 immediately  after  all  of the tied values.  Ties for the smallest values are
 handled in an analogous manner  using  the  character  -.  If  the  number  of
 observations  in  a group equals one,  ESTIMATE NOT AVAILABLE is printed under
 WITHINS S.D.  and S.D.  OF MEAN.  Also, ********** TO ********** appears under
 95PCT CONF INT FOR MEAN.

      The total number  of observations, mean, minimum observation  and maximum
 observation are also printed for the whole dataset combined.  In addition, the
 within standard  deviation,  standard deviation  of the  mean  and  95-percent
 confidence  interval for the mean are  printed for three different models:   a
 fixed effects model (Model I), a random  effects model (MODEL II) and  a model
 which assumes that all observations  are from a single group.   The confidence
 limits are formed by taking the grand mean plus (and minus) the product of the
 percentage point of Student's t distribution and the standard deviation of the
 mean.  Let k be the number of groups and n be the total number of observations
 with positive tag.   Then, the  standard deviation of the mean  is the square
 root of the variance of the mean formed as follows:

   Model                    Variance                 Variance of     Degrees of
                                                        Mean           Freedom

     I           VI  = Within groups mean square        VI/n             n-k

                        k
    II           VII = SUM (Y(i)-mean(Y))**2/(k-1)      VII/k            k-1
                       i=1

 Ungrouped       Vu  = Total mean square                Vu/n             n-1

                                     <6-4>
1

      PAIRWISE MULTIPLE COMPARISON OF MEANS.   This section only appears if the
 significance  level (value of F PROB.)  of the between groups F-ratio  is less
 than 0.10.   The Newman-Keuls-Hartley procedure is not performed if the number
 of  measurements with  positive tag  is  less than  four plus  the  number  of
 groups.

      The purpose of this comparison is to divide the groups in such a way that
 the group means within a division are not significantly different at  the 0.05
 significance level,  whereas  the  group  means  in  different  divisions  are
 significantly different at the 0.05 level.  Two different procedures are used:
 the  Newman-Keuls-Hartley method and the Scheffe' method.  The two methods are
 similar but not identical and frequently give slightly different results.  The
 Newman-Keuls-Hartley method is described in section 10.6  of  Snedecor  [1956]
 and  section  10.8  of  Snedecor  and  Cochran [1967].  The Scheffe' method is
 discussed in section 10.3 of  Brownlee  [1965].  Groups  are  separated  by  a
 string of five asterisks.  If two divisions have no group means in common, the
 two divisions are separated by two strings of five asterisks.

      Both the  Newman-Keuls-Hartley method  and  the  Scheffe'  method require
 percentage  points of the  studentized range:  an  approximation developed  by
 Mandel  is used  to compute them.   Since the  Newman-Keuls-Hartley method  is
 designed for  use when the number of  observations in each group is  the same,
 the number  of observations in each of  the two  groups is approximated  by m,
 where

                           1/m = (1/2)*(1/mi +  1/mj)

 and mi  and  mj are  the actual  number  of measurements  in each  of  the two
 groups.


      TEST FOR HOMOGENEITY OF VARIANCES.   The usual analysis of variance for a
 one-way classification assumes that the variances of all groups are  the same.
 This section of output provides information for assessing the validity of this
 assumption.   Small  values  of the  significance level  P  indicate  lack  of
 homogeneity of variance.   The Cochran's  C statistic printed is discussed  on
 page 180 of  Dixon and  Massey [1957]  and  in more  detail in  chapter  15 of
 Eisenhart  et al.  [1947].   The  Bartlett-Box  F-test  is  a modification  of
 Bartlett's test  which  uses the  F-distribution rather  than  the chi-squared
 distribution and is less sensitive to non-normality.  It is discussed on pages
 179 and  180  of Dixon  and Massey  [1957].   A  table  of critical  values of
 (maximum variance)/(minimum variance) for equal sample sizes is given on pages
 100 and 101 of Owen [1962].

      If either P value is less than or equal to 0.10, the  approximate between
 mean  F-test in the  presence of  heterogeneous variance and  its significance
 level P are also printed.  This approximate F-test for testing for differences
 between  means  is described  on  pages  287-289  of  Snedecor  [1956].   This
 information does not appear in figure F-1b because both P values (significance
 levels) exceed 0.10.


      MODEL  II  -  COMPONENTS  OF  VARIANCE.  This  is  the  usual analysis of
 variance estimate for the between component in a random effects  model  (Model
 II).  For a discussion of this analysis,  see section 10.6 and section 10.7 of
 Brownlee [1965].

                                     <6-5>
1

 F.  Example

      The example program below uses AOV1 to analyze 16 determinations  of  the
 gravitational  constant,  grouped  according  to the material used to make the
 measurements.  A discussion of this example can be found on pages  314-316  of
 Brownlee [1965].



















































                                     <6-6>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE AOV1 USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y AND TAG MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(20), TAG(20)
       DOUBLE PRECISION DSTAK(200)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED DATA
 C          TAG VALUES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
       READ (5,101) (TAG(I), I=1,N)
 C
 C     PRINT TITLE AND CALL AOV1 FOR ANALYSIS OF VARIANCE
 C
       WRITE (IPRT,102)
       CALL AOV1 (Y, TAG, N, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (20F5.1)
   102 FORMAT ('1RESULTS OF STARPAC ONE-WAY ANALYSIS OF VARIANCE',
      *  ' SUBROUTINE AOV1')
       END


 Data:

    16
  83.0 81.0 76.0 78.0 79.0 72.0 61.0 61.0 67.0 67.0 64.0 78.0 71.0 75.0 72.0 74.0
   1.0  1.0  1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0  2.0  3.0  3.0  3.0  3.0  3.0







                                     <6-7>
1RESULTS OF STARPAC ONE-WAY ANALYSIS OF VARIANCE SUBROUTINE AOV1
                                                                                                         STARPAC 2.08S (03/15/90)



                                                ANALYSIS OF VARIANCE


 *GROUP NUMBERS HAVE BEEN ASSIGNED ACCORDING TO TAG VALUES GIVEN, WHERE THE SMALLEST TAG GREATER THAN ZERO HAS BEEN ASSIGNED  *
 *GROUP NUMBER 1, THE NEXT SMALLEST, GROUP NUMBER 2, ETC.  TAGS LESS THAN OR EQUAL TO ZERO HAVE NOT BEEN INCLUDED IN ANALYSIS.*
 *NUMBER OF VALUES EXCLUDED FROM ANALYSIS IS    0                                                                             *

                 SOURCE              D.F.    SUM OF SQUARES     MEAN SQUARES         F RATIO    F PROB.

                 BETWEEN GROUPS        2      5.651042E+02      2.825521E+02       .261E+02      .000
                    SLOPE                 1      6.450893E+01      6.450893E+01      .141E+01      .257
                    DEVS. ABOUT LINE      1      5.005952E+02      5.005952E+02      .462E+02      .000
                 WITHIN GROUPS        13      1.408333E+02      1.083333E+01
                 TOTAL                15      7.059375E+02


           KRUSKAL-WALLIS RANK TEST FOR DIFFERENCE BETWEEN GROUP MEANS * H =   .114E+02, F PROB =  .000 (APPROX.)

                                                       ESTIMATES
                                                                                                SUM OF
      TAG           NO.      MEAN       WITHIN S.D.  S.D. OF MEAN     MINIMUM       MAXIMUM      RANKS   95PCT CONF INT FOR MEAN

   1.000000E+00       6   7.81667E+01+  3.86868E+00+  1.57938E+00   7.20000E+01   8.30000E+01     76.0  7.41067E+01 TO 8.22266E+01
   2.000000E+00       5   6.40000E+01-  3.00000E+00   1.34164E+00   6.10000E+01   6.70000E+01     15.0  6.02750E+01 TO 6.77250E+01
   3.000000E+00       5   7.40000E+01   2.73861E+00-  1.22474E+00   7.10000E+01   7.80000E+01     45.0  7.05996E+01 TO 7.74004E+01

           TOTAL     16   7.24375E+01                               6.10000E+01   8.30000E+01

                 FIXED EFFECTS MODEL    3.29140E+00   8.22851E-01                                       7.06594E+01 TO 7.42156E+01
                 RANDOM EFFECTS MODEL   7.29576E+00   4.21221E+00                                       5.43138E+01 TO 9.05612E+01
                 UNGROUPED DATA         6.86021E+00   1.71505E+00                                       6.87815E+01 TO 7.60935E+01

 PAIRWISE MULTIPLE COMPARISON OF MEANS.  THE MEANS ARE PUT IN INCREASING ORDER IN GROUPS SEPARATED BY *****.  A MEAN IS
 ADJUDGED NON-SIGNIFICANTLY DIFFERENT FROM ANY MEAN IN THE SAME GROUP AND SIGNIFICANTLY DIFFERENT AT THE .05 LEVEL FROM
 ANY MEAN IN ANOTHER GROUP.  ***** ***** INDICATES ADJACENT GROUPS HAVE NO COMMON MEAN.

   NEWMAN-KEULS TECHNIQUE, HARTLEY MODIFICATION. (APPROXIMATE IF GROUP NUMBERS ARE UNEQUAL.)
    6.40000E+01,
   ***** *****
    7.40000E+01, 7.81667E+01,

   SCHEFFE TECHNIQUE.
    6.40000E+01,
   ***** *****
    7.40000E+01, 7.81667E+01,

 TESTS FOR HOMOGENEITY OF VARIANCES.
       COCHRANS C = MAX. VARIANCE/SUM(VARIANCES) =   .4756, P =   .395 (APPROX.)
       BARTLETT-BOX F =      .269, P =   .764
       MAXIMUM VARIANCE / MINIMUM VARIANCE =         1.9956

 MODEL II - COMPONENTS OF VARIANCE.
       ESTIMATE OF BETWEEN COMPONENT   5.114706E+01

                                                              <6-8>
1G.  Acknowledgments

      The code and output for the one-way analysis of variance  subroutines  is
 adapted  from OMNITAB II [Hogben et al., 1971].  The discussion of the results
 (section E.2) is taken from the OMNITAB II User's Reference Manual [Hogben  et
 al., 1971].





















































                                     <6-9>
1-----                             CHAPTER 7                              -----

                              CORRELATION ANALYSIS


 A.  Introduction

 STARPAC contains two subroutines for correlation analysis of a multivariate 
 random sample.   The analysis  provided by  these subroutines consists  of
 seven tables,  which, when used  together, aid  the user in  using correlation
 techniques effectively for prediction and model building.  The analysis is the
 same as that provided by  the OMNITAB II  command  CORRELATION [Hogben  et al.
 1971].  For further information on correlation techniques users should consult
 Kendall and Stuart [1973], Brownlee [1965] and Anderson [1958].

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      Subroutine CORR computes and prints 1) the simple correlation  matrix and
 2)  the significance levels  of the  simple correlation  coefficients;  3) the
 partial correlation  coefficients  and 4)  their significance  levels;  5) the
 Spearman rank correlation coefficients; 6) a test for a quadratic relationship
 among the variables; and 7) 95-percent and 99-percent confidence intervals for
 the simple correlation coefficients.

      Subroutine  CORRS provides  the  same analysis  as CORR  but  returns the
 variance-covariance matrix used to compute the correlation coefficients.   The
 user can also optionally suppress the printed output.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.


 CORR:   Compute  and  print a  correlation analysis of  a multivariate  random
         sample

         <real> YM(n,m)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL CORR (YM, N, M, IYM, LDSTAK)

                                      ===





                                     <7-1>
1CORRS:  Compute and optionally print a  correlation analysis of a multivariate
         random sample; return variance-covariance matrix

         <real> YM(n,m), VCV(m,m)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL CORRS (YM, N, M, IYM, LDSTAK, NPRT, VCV, IVCV)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the  argument is input to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.]  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 IVCV    --> The exact  value of  the  first dimension  of the  matrix  VCV  as
             specified in the calling program.

 IYM     --> The exact  value of  the  first  dimension  of the  matrix  YM  as
             specified in the calling  program.

 LDSTAK  --> The length of the double precision workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision version of  STARPAC is  being  used  P = 0.5,
             otherwise P = 1.0.  [See chapter 1, section B.]

             For CORR:

               LDSTAK >= (47+max(N,M))/2 + (max(N,M)+N*(M+3)+M+7*M**2)*P





                                     <7-2>
1            For CORRS:

               LDSTAK >= (47+IO*max(N,M))/2 + IO*(max(N,M)+N*(M+3)+M+6*M**2)*P

               where IO = 0 if NPRT = 0 and IO = 1 if NPRT <> 0.

 M       --> The number of  variables measured  for each observation,  that is,
             the number of columns of data in YM.

 N       --> The number of observations, that is, the number of rows of data in
             YM.

 NPRT    --> The argument controlling printed output.

             If NPRT = 0 the printed output is suppressed.

             If NPRT ] 0 the printed output is provided.

 VCV     <-- The  matrix  of dimension  at  least  M  by M  that  contains  the
             variance-covariance matrix,

                          1      N
             VCV(j,k) = ----- * SUM (YM(i,j)-YjMEAN)*(YM(i,k)-YkMEAN)
                        (N-1)   i=1

                            1    N                       1    N
             where YjMEAN = - * SUM YM(i,j) and YkMEAN = - * SUM YM(i,k).
                            N   i=1                      N   i=1

 YM      --> The matrix of dimension at least N by M that contains the observed
             multivariate data.   The element in the ith row and  jth column is
             the ith observation on the jth variable.


 E.  Computational Methods

 E.1  Algorithms

      Formulas for the computed tables are given in section E.2.  The code  and
 output  for  the  correlation analysis subroutines are adapted from OMNITAB II
 [Hogben et al., 1971].


 E.2  Computed Results and Printed Output

      The argument controlling  the  printed  output,  NPRT,  is  discussed  in
 section D.

      STARPAC correlation  analysis subroutines compute and print  a seven-part
 correlation  analysis.  The output  for each  part  is discussed  individually
 below;  the text  for this  discussion  is taken  from the  OMNITAB  II User's
 Reference Manual [Hogben et al., 1971].

      The  Simple  Correlation  Coefficient  Matrix.  The (j,k)th entry of this
 matrix is the simple (product moment) correlation  coefficient,  rjk,  between
 the data in columns j and k defined by

                    rjk = VCV(j,k)/sqrt( VCV(j,j)*VCV(k,k) )

                                     <7-3>
1
 Note that when more than two variables are under study, the simple correlation
 coefficient can be seriously distorted by the effect of other variables.   The
 partial  correlation coefficient  (see  below) can  be used  to  identify such
 distortion.

      The  Significance  Levels  of  the  Simple Correlation Coefficients.  The
 (j,k)th entry of  this  table  is  the  significance  level,  Sr(j,k)  of  the
 corresponding partial correlation coefficient, rjk,

                 Sr(j,k) = probability of exceeding F0(1,N-2)

 where F0(1,N-2) is an F-statistic with 1 and N-2 degrees of freedom,

                                        (N-2)*rjk**2
                            F0(1,N-2) = ------------
                                          1-rjk**2

 If  the "true" correlation coefficient is  equal to zero, then Sr(j,k)  is the
 probability that in a random sample (of the same size) the absolute value of a
 sample correlation coefficient will exceed the absolute value of  the observed
 correlation coefficient, rjk.

      The   Partial   Correlation   Coefficients.   The   partial   correlation
 coefficient,  pjk,  between the data in  columns  j  and  k,  j<>k,  with  the
 remaining variables fixed, i.e., held constant, is given by

                            pjk = -cjk/sqrt(cjj*ckk)

 where cjk denotes the (j,k)th element of the inverse of the simple correlation
 matrix.   Because the partial correlation coefficient measures the correlation
 between two variables after eliminating the effect of the remaining variables,
 it  may avoid the  distortion suffered  by the simple  correlation coefficient
 when more  than  two variables  are under  study.   The user  should therefore
 compare  the  simple correlation  coefficients with  the  partial  correlation
 coefficients.   Any "large"  discrepancy indicates  that one  or  more of  the
 remaining variables is having an  important effect on the relationship.   [See
 Kendall and Stuart, 1961, section 27.5, page 318.]

      The  Significance  Levels  of  the Partial Correlation Coefficients.  The
 (j,k)th entry of  this  table  is  the  significance  level,  Sp(j,k)  of  the
 corresponding partial correlation coefficient, pjk,

                 Sp(j,k) = probability of exceeding F0(1,N-M)

 where F0(1,N-M) is an F-statistic with 1 and N-M degrees of freedom,

                                        (N-M)*pjk**2
                            F0(1,N-M) = ------------
                                          1-pjk**2

 If the "true" partial correlation  coefficient is equal to zero,  then Sp(j,k)
 is the  probability that in a random  sample (of  the same size)  the absolute
 value of a partial correlation  coefficient will exceed the absolute  value of
 the observed partial correlation coefficient, pjk.

      Spearman Rank Correlation Coefficient.   The rank correlation coefficient
 does  not  require the  assumption  that  the  data  have a  bivariate  normal

                                     <7-4>
1distribution.  The Spearman rank correlation coefficient, sjk, for the data in
 columns j and k is computed from

                               A - Djk**2 - Tj - Tk
                     sjk = ----------------------------- ,
                           sqrt( (A - 2*Tj)*(A - 2*Tk) )

 where

                          N
                Djk**2 = SUM (rank(YM(i,j)) - rank(YM(i,k)))**2,
                         i=1

                              A = (N-1)*N*(N+1)/6,

                       Tj = (1/12) SUM (tj-1)*tj*(tj+1),
                                    j

                       Tk = (1/12) SUM (tk-1)*tk*(tk+1),
                                    k

       tj = number of ties in a set of tied values in column j of YM, and

       tk = number of ties in a set of tied values in column k of YM.

      The quantities Tj and Tk adjust  for ties in the ranks.   If there are no
 ties, Tj and Tk equal zero and

                             sjk = 1 - (Djk**2 / A).

 A comparison should be made between the rank correlation coefficients  and the
 corresponding  simple and partial correlation coefficients.   Again, a "large"
 discrepancy between  two  comparable  coefficients  is an  indicator  of  some
 abnormality in the data.  See Kendall [1948] for further details.

      Significance Level of Quadratic Fit over Linear Fit.   Underlying the use
 of a  correlation coefficient is  the assumption  that the  two  variables are
 linearly  related.   The results  in this  part  are useful  in assessing  the
 validity of this assumption of linearity.  The variables are all assumed to be
 normally distributed.   The numbers printed are the significance levels  for a
 F-test of the hypothesis that the quadratic term in a quadratic model is zero.
 The F-statistic used is

                          RSS(linear model) - RSS(quadratic model)
              F0(1,N-3) = ----------------------------------------
                                 RSS(quadratic model)/(N-3)

 with 1 and (N-3) degrees of freedom, where RSS is the residual sum  of squares
 function.   The  values of  1 and  (N-3)  are  printed  in the  heading.   The
 significance level, SQ, is then computed as

                               SQ = Prob(F > F0).

 Small values of the significance level (less than 0.05, for  example) indicate
 lack of linearity.  The test results differ depending upon which variable of a
 pair is  considered the dependent  variable and  which one  is  considered the
 independent  (or predictor)  variable.   Hence, the  entire table  is printed,
 rather than just the lower triangle.  The diagonal entries are always equal to

                                     <7-5>
1one  and  have  no  particular  relevance.   Tests  of  hypotheses  in  linear
 regression are discussed in section 13.8 of Brownlee [1965].

     Confidence Intervals For Simple Correlation Coefficients.  Both 95-percent
 and 99-percent  confidence intervals  for the simple  correlation coefficients
 are printed in this two-way table.   There are two entries in each cell of the
 table.   The values  0.95 and 0.99 are printed  along the upper left to  lower
 right  diagonal.   The  95-percent  confidence  limits are  printed  below the
 diagonal and the 99-percent confidence limits are printed above  the diagonal.
 The number in  the lower left of each  cell is the lower confidence  limit and
 the number in the upper right is the upper confidence limit.

      The  confidence intervals  are  based on  a normal  approximation.   [See
 Morrison, 1967, chapter 3, page 101.]  They are computed as follows:

                    Lower confidence limit:  tanh( z - u/sqrt(N-3) )

                    Upper confidence limit:  tanh( z + u/sqrt(N-3) )

 where

                       inv(tanh)(rjk)
                 z = <
                       0.5 * ln( (1+rjk)/(1-rjk) )
 and

                       1.96 for 95-percent confidence interval
                 u = <
                       2.58 for 99-percent confidence interval.


 F.  Example

      The sample program shown below illustrates the use of CORR.  The data are
 taken from Draper and Smith [1968], page 216.  The data are part of a study to
 determine the effect of relative urbanization, educational level, and relative
 income (columns 1, 2, and 3, respectively) on product usage (column 4).






















                                     <7-6>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE CORR USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF YM MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL YM(20,6)
       DOUBLE PRECISION DSTAK(200)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
       IYM = 20
 C
 C     READ NUMBER OF OBSERVATIONS AND NUMBER OF VARIABLES
 C          OBSERVED MULTIVARIATE DATA
 C
       READ (5,100) N, M
       READ (5,101) ((YM(I,J), J=1,M), I=1,N)
 C
 C     PRINT TITLE AND CALL CORR TO PERFORM CORRELATION ANALYSIS
 C
       WRITE (IPRT,102)
       CALL CORR (YM, N, M, IYM, LDSTAK)
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (4F6.1)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' CORRELATION ANALYSIS SUBROUTINE CORR')
       END
















                                     <7-7>
1Data:

     9    4
  42.2  11.2  31.9 167.1
  48.6  10.6  13.2 174.4
  42.6  10.6  28.7 160.8
  39.0  10.4  26.1 162.0
  34.7   9.3  30.1 140.8
  44.5  10.8   8.5 174.6
  39.1  10.7  24.3 163.7
  40.1  10.0  18.6 174.5
  45.9  12.0  20.4 185.7















































                                     <7-8>
1RESULTS OF STARPAC CORRELATION ANALYSIS SUBROUTINE CORR
                                                                                                         STARPAC 2.08S (03/15/90)

 CORRELATION ANALYSIS FOR  4 VARIABLES WITH    9 OBSERVATIONS


 CORRELATION MATRIX
    - STANDARD DEVIATIONS ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

    COLUMN         1                2                3                4

         1     4.1764552
         2     .68374212        .74628711
         3    -.61596989       -.17249295        7.9279218
         4     .80175221        .76795025       -.62874595        12.645157


 SIGNIFICANCE LEVELS OF SIMPLE CORRELATION COEFFICIENTS (ASSUMING NORMALITY)

    COLUMN         1                2                3                4

         1    0.
         2     .42269734E-01   0.
         3     .77358710E-01    .65719690       0.
         4     .93530479E-02    .15660429E-01    .69711268E-01   0.


 PARTIAL CORRELATION COEFFICIENTS WITH  2 REMAINING VARIABLES FIXED

    COLUMN         1                2                3                4

         1     1.0000000
         2     .43170939        1.0000000
         3    -.45663591        .69717000        1.0000000
         4     .10539042        .72682008       -.64778927        1.0000000


 SIGNIFICANCE LEVELS OF PARTIAL CORRELATION COEFFICIENTS (ASSUMING NORMALITY)

    COLUMN         1                2                3                4

         1    0.
         2     .33344895       0.
         3     .30301787        .81676691E-01   0.
         4     .82207563        .64249351E-01    .11565980       0.













                                                              <7-9>
1SPEARMAN RANK CORRELATION COEFFICIENTS (ADJUSTED FOR TIES)

    COLUMN         1                2                3                4

         1     1.0000000
         2     .61088401        1.0000000
         3    -.56666667       -.12552411        1.0000000
         4     .68333333        .60251573       -.71666667        1.0000000


 SIGNIFICANCE LEVEL OF QUADRATIC FIT OVER LINEAR FIT BASED ON F RATIO WITH 1 AND    6 DEGREES OF FREEDOM
 (FOR EXAMPLE,  .1703 IS THE SIGNIFICANCE LEVEL OF THE QUADRATIC TERM WHEN COLUMN  2 IS FITTED TO COLUMN  1)

    COLUMN         1                2                3                4

         1     1.0000000        .40442827        .94936019        .85222263
         2     .17034675        1.0000000        .80987340        .93773837
         3     .71654009        .56763094        1.0000000        .84988526
         4     .15652586        .59973987        .36810270        1.0000000


 CONFIDENCE INTERVALS FOR SIMPLE CORRELATION COEFFICIENTS (USING FISHER TRANSFORMATION)
 95 PER CENT LIMITS BELOW DIAGONAL, 99 PER CENT LIMITS ABOVE DIAGONAL

    COLUMN         1                2                3                4

         1     99.000000        .95517075        .32129731        .97349304
             95.000000       -.21219610       -.94361630        .51874096E-01

         2     .92694792        99.000000        .70508573        .96846093
             .35940598E-01    95.000000       -.84136050       -.36249864E-01

         3     .81486051E-01    .55523436        99.000000        .30247207
            -.90845980       -.75062574        95.000000       -.94585730

         4     .95654889        .94838422        .60737579E-01    99.000000
             .29437227        .21190035       -.91203488        95.000000






















                                                              <7-10>
1G.  Acknowledgments

      The  code  and  output  for  the  correlation subroutines is adapted from
 OMNITAB II [Hogben et al., 1971].  The discussion of the results (section E.2)
 is taken from the OMNITAB II Users' Reference Manual [Hogben et al., 1971].






















































                                     <7-11>
1-----                             CHAPTER 8                              -----

                              LINEAR LEAST SQUARES


 A.  Introduction

      STARPAC contains  eight subroutines  for linear  least  squares analysis.
 For four of these, the user specifies the model by supplying the design matrix
 (the matrix whose columns are the independent variables plus a column  of ones
 if a  constant term  is being  estimated).   The other  four perform  the same
 analysis for  the special case of a  polynomial model, where the need  for the
 user to explicitly create the design matrix is eliminated.

      Each of the subroutines  described  in  this  chapter  assumes  that  the
 observations of the dependent variable,  Y(i),  which are measured with error,
 are modeled by

                     NPAR
              Y(i) = SUM  PAR(j)*x(i,j) + e(i)  for i = 1, ..., N,
                     j=1

 where

 N       is the number of observations;

 NPAR    is the number of parameters;

 x(i,j)  is  the jth element of the ith row of the design matrix (for the 
         user-specified model, x(i,j) = XM(i,j) for i = 1,  ...,  N and j = 1,  ...,
         NPAR,  while  for  the polynomial model,  x(i,j) = X(i)**(j-1) for i =
         1,..., N and j = 1, ..., NPAR);

 PAR     is the vector of the NPAR unknown parameters (coefficients); and

 e(i)    is the random error in the ith observation.

 The least squares solution,  PARhat, is that which minimizes (with respect to
 PAR) the residual sum of squares function

                         N
             RSS(PAR) = SUM e(i)**2*wt(i)
                        i=1

                         N                 NPAR
                      = SUM wt(i)*( Y(i) - SUM PAR(j)*x(i,j) )**2
                        i=1                j=1

 where ''hat'' (e.g., PARhat, Yhat, etc.) denotes the estimated quantity, and

 wt(i)   is the  weight  assigned  to the  ith  observation  (wt(i)  =  1.0  in
         the  ''unweighted''  case).   Appendix  B   discusses   several common
         applications for weighted least squares.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The declaration and CALL statements are given in section  C  and
 the  subroutine  arguments  are defined in section D.  The algorithms used and


                                     <8-1>
1output produced by these  subroutines  are  discussed  in  section  E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      The linear least squares estimation subroutines permit both  weighted and
 unweighted analysis.  The user has two levels of control over the computations
 and printed output.

         * In level one,  a four-part printed report is  automatically provided
           and  the residuals  are  returned to  the user  via  the  subroutine
           argument list.

         * Level two also returns the residuals.  It allows the user to specify
           the  amount  of  printed  output,  and,  in  addition,  returns  the
           following estimated values via the argument list:
               - the estimated parameters;
               - the residual standard deviation;
               - the predicted values;
               - the standard deviations of the predicted values;
               - the standardized residuals; and
               - the variance-covariance matrix of the estimated parameters.

      The simplest  of the linear  least squares  subroutines are  LLS  for the
 user-specified  model  and  LLSP  for  the  polynomial  model.   They  perform
 unweighted analyses,  provide  a  four-part  printed  report  and  return  the
 residuals  via  the  argument  list  (level  one  control).    The  other  six
 subroutines provide greater  flexibility by adding the weighting  and/or level
 two control features.  These features are each indicated by a different suffix
 letter on the subroutine name (e.g., LLSS and LLSPWS).

         * Suffix W indicates user-supplied weights.

         * Suffix S indicates level two control of the computations and output.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument  definitions  and  sample programs are given in sections D and
 section F,  respectively.  The  conventions  used  to  present  the  following
 declaration and CALL statements are given in chapter 1, sections B and D.


 LLS:    Compute  and  print  a  four-part   unweighted  linear  least  squares
         analysis with user-specified model (design matrix); return residuals

         <real> Y(n), XM(n,npar), RES(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLS (Y, XM, N, IXM, NPAR, RES, LDSTAK)

                                      ===




                                     <8-2>
1LLSS:   Compute  and  optionally  print a  four-part unweighted  linear  least
         squares  analysis with  user-specified model  (design  matrix); return
         residuals, parameter estimates, residual standard deviation, predicted
         values,  standard deviations  of the  predicted  values,  standardized
         residuals, and variance-covariance matrix of parameters

         <real> Y(n), XM(n,npar), RES(n)
         <real> PAR(npar), PV(n), SDPV(n), SDRES(n), VCV(npar,npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSS (Y, XM, N, IXM, NPAR, RES, LDSTAK,
        +           NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 LLSW:   Compute  and print a four-part  weighted linear least squares analysis
         with user-specified model (design matrix); return residuals

         <real> Y(n), XM(n,npar), RES(n)
         <real> WT(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSW (Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)

                                      ===

 LLSWS:  Compute  and  optionally  print  a  four-part  weighted  linear  least
         squares  analysis with  user-specified model  (design  matrix); return
         residuals, parameter estimates, residual standard deviation, predicted
         values,  standard deviations  of the  predicted  values,  standardized
         residuals, and variance-covariance matrix of parameters

         <real> Y(n), XM(n,npar), RES(n)
         <real> WT(n), PAR(npar), PV(n), SDPV(n), SDRES(n), VCV(npar,npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSWS (Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
        +            NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 LLSP:   Compute  and   print  a  four-part  unweighted  linear  least  squares
         analysis with polynomial model (design matrix); return residuals

         <real> Y(n), X(n), RES(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSP (Y, X, N, NDEG, RES, LDSTAK)

                                      ===

                                     <8-3>
1
 LLSPS:  Compute  and  optionally  print a  four-part unweighted  linear  least
         squares  analysis  with  polynomial  model  (design   matrix);  return
         residuals, parameter estimates, residual standard deviation, predicted
         values,  standard deviations  of the  predicted  values,  standardized
         residuals, and variance-covariance matrix of parameters

         <real> Y(n), X(n), RES(n)
         <real> PAR(npar), PV(n), SDPV(n), SDRES(n), VCV(npar,npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSPS (Y, X, N, NDEG, RES, LDSTAK,
        +            NPRT, LPAR, PAR, NPAR, RSD,
        +            PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 LLSPW:  Compute  and  print a four-part weighted linear least squares analysis
         with polynomial model (design matrix); return residuals

         <real> Y(n), X(n), RES(n)
         <real> WT(n)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSPW (Y, WT, X, N, NDEG, RES, LDSTAK)

                                      ===

 LLSPWS: Compute  and  optionally  print  a  four-part  weighted  linear  least
         squares  analysis  with  polynomial  model  (design   matrix);  return
         residuals, parameter estimates, residual standard deviation, predicted
         values,  standard deviations  of the  predicted  values,  standardized
         residuals, and variance-covariance matrix of parameters

         <real> Y(n), X(n), RES(n)
         <real> WT(n), PAR(npar), PV(n), SDPV(n), SDRES(n), VCV(npar,npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         :
         :
         CALL LLSPWS (Y, WT, X, N, NDEG, RES, LDSTAK,
        +             NPRT, LPAR, PAR, NPAR, RSD,
        +             PV, SDPV, SDRES, VCV, IVCV)

                                      ===










                                     <8-4>
1D.    Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates that the  argument is input  to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 IERR    ... An  error  flag  returned  in  COMMON  /ERRCHK/  [see  chapter  1,
             section D.5].  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates  that  no errors  were detected  and  that  the
                      least squares solution was found.

             IERR = 1 indicates that improper input was detected.

             IERR = 2 indicates that  the model  is  computationally  singular,
                      which may mean the  model has  too many parameters.   The
                      user should examine  the model and data to  determine and
                      remove the cause of the singularity.

             IERR = 3 indicates that at least one of the standardized residuals
                      could not be computed because its standard  deviation was
                      zero.   The validity of the variance-covariance matrix is
                      questionable.

 IVCV    --> The exact  value of  the  first dimension  of the  matrix  VCV  as
             specified in the calling program.

 IXM     --> The exact  value of  the  first  dimension  of the  matrix  XM  as
             specified in the calling program.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed the  value given  below, where if  the single
             precision  version  of STARPAC  is being  used  P = 0.5, otherwise
             P = 1.0 [see chapter 1, section B].

             LDSTAK >= 28 + [6*N+1+5*NPAR+NPAR*N+2*NPAR**2]*P

 LPAR    --> The actual length of the  vector PAR  as specified in  the calling
             program.

 N       --> The number of observations.

 NDEG    --> The  degree of  the polynomial  model.   The  number  of estimated
             parameters is NPAR = NDEG + 1.

                                     <8-5>
1
 NPAR    --- The number  of parameters to be estimated.   NPAR is  input to the
             subroutines with a user-specified model; for the  subroutines with
             a polynomial model, NPAR = NDEG + 1 is returned.

 NPRT    --> The argument  controlling printed  output.   NPRT is  a four-digit
             integer, where the value of  the Ith digit (counting from  left to
             right) controls the Ith section of the output.

             If the Ith digit =  0, the output from the Ith section is
                                    suppressed;
                              =  1, the brief form of the Ith section is given;
                              >= 2, the full form of the Ith section is given.

             The default value for NPRT is 1112.  If the user-supplied value of
             NPRT  is  less  than zero  or  NPRT  is  not an  argument  in  the
             subroutine CALL statement the default value will be used.

             A  full  discussion  of the printed output is given in section E.2
             and is summarized as follows.

             Section 1 provides  information for each observation based  on the
                       solution.   Brief  output includes  information for  the
                       first  40 observations,  while full output  provides the
                       information for all of the data.

             Section 2 is a set of four residual plots.   Brief output and full
                       output are the same for this section.

             Section 3 is  an  analysis of  variance.   Brief  output  and full
                       output are the same for this section.

             Section 4 is the final summary of the estimated parameters.  Brief
                       output does not include printing the variance-covariance
                       matrix while full output does.

 PAR     <-- The vector of dimension at least NPAR that contains  the estimated
             parameter values.

 PV      <-- The  vector of dimension  at least  N that contains  the predicted
             values of the dependent variable,

                     NPAR
             PV(i) = SUM  PARhat(j)*x(i,j) = Yhat(i) for i = 1, ..., N.
                     j=1

 RES     <-- The vector of dimension at least N that contains the  residuals at
             the solution,

                             NPAR
             RES(i) = Y(i) - SUM  PARhat(j) x(i,j)
                             j=1

                    = Y(i) - Yhat(i) = ehat(i)

             for i = 1, ..., N.



                                     <8-6>
1RSD     <-- The residual standard deviation at the solution,

             RSD = sqrt( RSS(PARhat)/(Nnzw-NPAR) )

             where Nnzw is the number of observations with nonzero weights.

 SDPV    <-- The  vector of dimension  at least  N that  contains  the standard
             deviation of each predicted value at the solution,

             SDPV(i) = the ith diagonal element of sqrt[D*VCV*trans(D)]

             for i = 1, ..., N,  where D is the design matrix,  D(i,j) = x(i,j)
             with x(i,j) defined in section A, and trans(D) is the transpose of
             D.

 SDRES   <-- The vector of dimension at least N that contains  the standardized
             residuals at the solution,  i.e., the ith residual divided  by its
             estimated standard deviation,

             SDRES(i) = RES(i)/sqrt( RSD**2/wt(i) - SDPV(i)**2 )

             for i = 1, ..., N.

 VCV     <-- The matrix of dimension at  least NPAR  by NPAR that  contains the
             variance-covariance  matrix of  the estimated  parameters  at  the
             solution,

             VCV = RSD**2*inv(trans(D)*W*D)

             where W is the N by N diagonal matrix of weights,

             W = diag( wt(i), i=1, ..., N),

             and D is the design matrix, D(i,j) = x(i,j) with x(i,j) defined in
             section A, trans(D) is the transpose  of  D,  and  inv(.)  is  the
             inverse of the designated matrix.

 WT      --> The  vector of dimension  at least  N that  contains  the weights.
             Negative weights are not allowed and the number of nonzero weights
             must equal or exceed the number of parameters being estimated.   A
             zero  weight eliminates  the corresponding  observation  from  the
             analysis,  although  the residual,  the predicted  value  and  the
             standard  deviation  of the  predicted value  of  a  zero-weighted
             observation are still computed [see Appendix B].

 X       --> The vector of dimension  at least N that contains  the independent
             variable used to  construct the  design matrix for  the polynomial
             model.

 XM      --> The matrix  of dimension  at  least N  by NPAR  that  contains the
             design matrix, i.e., the matrix whose columns are  the independent
             variables plus  a column  of  ones if  a constant  term  is  being
             estimated.

 Y       --> The  vector of dimension  at least  N that contains  the dependent
             variable.



                                     <8-7>
1E.  Computational Methods

 E.1.  The Linear Least Squares Algorithm

      The  linear   least  squares   estimation  subroutines  use   a  modified
 Gram-Schmidt algorithm [Davis, 1962; Walsh, 1962].  The printed output for the
 linear least squares subroutines has been modeled on the linear  least squares
 output used by OMNITAB II [Hogben et al., 1971].


 E.2  Computed Results and Printed Output

      The argument controlling  the  printed  output,  NPRT,  is  discussed  in
 section D.

      The  output from the linear least squares estimation subroutines consists
 of four sections, several of which include tables summarizing the results.  In
 the following descriptions,  the  actual  table  headings  are  given  by  the
 uppercase   phrases   enclosed  in  angle  brackets  (<...>).   Results  which
 correspond to input  or  returned  subroutine  CALL  statement  arguments  are
 identified by the argument name in uppercase (not enclosed in angle braces).

 Section 1 provides the following information  for each observation, i, i  = 1,
           ..., N, based on the solution.

       * <ROW>:  the row number of the observation.

       * <PREDICTOR  VALUES>:  the  values for up to the first three columns of
         the independent variable  (design  matrix).  For  subroutines  with  a
         user-supplied  model,  this  is  up  to the first three columns of the
         matrix XM (excluding the first column if it is all ones,  indicating a
         constant  term);  for  the  polynomial model subroutines,  this is the
         variable X.

       * <DEPENDENT VARIABLE>:  the values of the dependent variable, Y.

       * <PREDICTED VALUE>:  the predicted values, PV, from the fit.

       * <STD DEV OF PRED VALUE>:  the  standard  deviations of  the  predicted
         values, SDPV.

       * <RESIDUAL>:  the error estimate, RES.

       * <STD RES>:  the standardized residual, SDRES.

       * <WEIGHT>:  the user-supplied  weights, WT, printed only when  weighted
         analysis is performed.


 Section 2 displays the following plots of the standardized residuals.

       * The standardized residuals versus row numbers.

       * The standardized residuals versus predicted values.

       * The autocorrelation function of the (non-standardized) residuals.

       * The normal probability plot of the standardized residuals.

                                     <8-8>
1

 Section 3 provides an  analysis of  variance.   The results  of this  analysis
           depend upon the order of the columns of the design matrix unless the
           columns  are  orthogonal.    The  analysis  includes  the  following
           information.

       * <PAR INDEX>:  the index, j, of the parameter being examined, PAR(j).

       * <SUM OF SQUARES RED DUE TO PAR>:  SSj,  the reduction in  the  sum  of
         squares due to fitting PAR(j)  after  having  fit  parameters  PAR(1),
         PAR(2),  ...,  PAR(j-1).  SSj  depends  on the order of the parameters
         unless  the  design  matrix  has  orthogonal  columns.   This   is   a
         decomposition of the total sum of squares, TSS, into NPAR + 1 parts

                             NPAR
         TSS = RSS(PARhat) + SUM SSj .
                             j=1

         The residual sum of squares and total sum of squares is also listed in
         this column.

       * <CUM MS RED>:  the cumulative mean square reduction,

                   j
         MSREDj = SUM SSk/j for j = 1, ..., NPAR.
                  k=1

       * <DF(MSRED)>:  the degrees of freedom associated  with  the  cumulative
         mean  square  reduction for each parameter,  DF(MSREDj) = j for j = 1,
         ..., NPAR.  The degrees of freedom for the residuals,  Nnzw-NPAR,  and
         the  total  degrees  of  freedom,  Nnzw,  where  Nnzw is the number of
         observations with nonzero weights, are also listed in this column.

       * <CUM RES MS>:  the cumulative residual mean square,

                                NPAR
         RMSj = (RSS(PARhat) +  SUM  SSk) / (Nnzw-j) for j = 1, ..., NPAR.
                               k=j+1

       * <DF(RMS)>:  the degrees of  freedom  associated  with  the  cumulative
         residual mean square for each parameter,  DF(RMSj) = Nnzw-j for j = 1,
         ..., NPAR.

       * <PAR=0>,  <F> and <PROB(F)>:  the F-ratio and its  significance  level
         under the null hypotheses that PAR(j) is zero after allowance has been
         made for parameters PAR(1), PAR(2), ..., PAR(j-1).  This F-ratio is

         Fj = SSj/RMSnpar for j = 1, ..., NPAR,

         with   1  and Nnzw-NPAR  degrees of  freedom.  The significance  level
         listed is the  probability of  exceeding the calculated  F-ratio under
         the null hypothesis that the corresponding parameter, PAR(j), is zero.

       * <PARS=0>,  <F> and <PROB(F)>:  the F-ratio and its significance  level
         under  the  null  hypothesis  that parameters PAR(j),  PAR(j+1),  ...,
         PAR(NPAR),  are zero after allowance  has  been  made  for  parameters
         PAR(1), PAR(2), ..., PAR(j-1).  This F-ratio is

                                     <8-9>
1
               NPAR
         Fj = (SUM SSk/(NPAR-j+1)) / RMSnpar for j = 1, ..., NPAR,
               k=j

         with  NPAR-j+1  and Nnzw-NPAR  degrees of  freedom.   The significance
         level listed is  the probability  of exceeding the  calculated F-ratio
         under the null hypothesis that all of the parameters PAR(j), PAR(j+1),
         ..., PAR(NPAR) are zero.

         The  numerator of this ratio is the extra sum of squares accounted for
         by  inclusion  of  the  terms  PAR(j)*x(j)  + PAR(j+1)*x(j+1) + ...  +
         PAR(npar)*x(npar) in the model, divided by its degrees of freedom; the
         denominator is the residual mean square of the full model.  This ratio
         is a means of comparing the extra sum of squares to its expected value
         as estimated by the residual mean square.  When the terms of the model
         have a logical order of entry,  this series of F-tests can be used  to
         judge  how  many terms should be included in the model [see Draper and
         Smith, 1981, pages 97 and 98].


 Section 4 summarizes  the  following  information about  the  final  parameter
           estimates and their variances.

       * The variance-covariance matrix, VCV, of the estimated  parameters, and
         the corresponding correlation matrix,

         rjk = VCV(j,k)/sqrt(VCV(j,j)*VCV(k,k)) for j = 1, ..., NPAR
                                                and k = 1, ..., NPAR.

       * <ESTIMATES FROM FIT>:

         - <ESTIMATED PARAMETER>: the final estimate for each parameter, PAR(j)
           for j = 1, ..., NPAR.

         - <SD OF PAR>:  the standard deviation  of  the  estimated  parameter,
           sqrt(VCV(j,j)) for j = 1, ..., NPAR.

         - <T(PAR=0)>: the Student's t statistic under the null hypothesis that
           PAR(j) is actually zero,

           T(PAR=0)j = PAR(j)/sqrt(VCV(j,j)) for j = 1, ..., NPAR.

         - <PROB(T)>:  the two-sided  significance level of T(PAR=0)j.  This is
           the probability of exceeding  the  given  t  value  under  the  null
           hypothesis that the parameter PAR(j) is actually zero.

         - <ACC DIG>:  an estimate of the number  of  reliable  digits  in  the
           parameter  estimates,  i.e.,  an  indication  of  the  computational
           accuracy of the solution.  A computationally accurate solution  will
           produce  values  between  DIGITS-2  and DIGITS,  where DIGITS is the
           number of decimal digits carried by the user's computer for a single
           precision value when the single  precision  version  of  STARPAC  is
           being  used  and  is the number carried for a double precision value
           otherwise.   Values   less   than   DIGITS-4   may   indicate   some
           computational difficulty such as poor scaling or near singularity.

       * <ESTIMATES  FROM  FIT  OMITTING  LAST  PREDICTOR  VALUE>:   <ESTIMATED

                                     <8-10>
1        PARAMETER>,  <SD  OF  PAR>,  <T(PAR=0)> and <PROB(T)> values for a fit
         omitting the last column of the design matrix and  thus  omitting  the
         last parameter from the model.

       * The residual standard deviation, RSD.

       * The residual degrees of  freedom, Nnzw-NPAR, where Nnzw is  the number
         of observations with nonzero weights.

       * The squared multiple correlation coefficient,

                                                               N
                                                              SUM wt(i)*Y(i)
                             RSS(PARhat)                      i=1
         R2 = 1.0  -  ------------------------     where Yw = -------------- .
                       N                                          N
                      SUM wt(i) * (Y(i)-Yw)**2                   SUM wt(i)
                      i=1                                        i=1

         R2 is a measure of how well the fitted equation accounts for the total
         variation of the dependent variable, Y.   It is only computed when the
         first parameter of the model is a constant, i.e., when the elements of
         the first column of the design matrix are all equal.

       * an approximation to the condition number of the design matrix, D(i,j)=
         x(i,j) with x(i,j) defined in section A, under the assumption that the
         absolute error in each column of D is roughly equal. The approximation 
	 will be meaningless if this assumption is not valid; otherwise it
         will usually underestimate the  actual condition number by a factor of
         from 2 to 10  [see Dongarra et al., 1979, p. 1.8 to 1.12].  (Note that
         the condition number returned by the linear least squares  subroutines
         is not  exactly the  same as  that  returned  by the  nonlinear  least
         subroutines  because of differences  in the  computational  procedures
         used by the two families of subroutines.)


 F.  Examples

      User-Specified Model (Design Matrix).  In the first example below, LLS is
 used to compute the least squares solution for the example given on pages  61-
 65 of Daniel and Wood [1971].  The results for this problem are also discussed
 in Draper and Smith [1981], pages 372 and 373.


      Polynomial Model (Design Matrix).  In the second example, LLSP is used to
 compute the least squares solution for the example given on page 311 of Miller
 and Freund [1977].












                                     <8-11>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE LLS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, XM AND RES MUST BE CHANGED TO
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(30), XM(30,5), RES(30)
       DOUBLE PRECISION DSTAK(500)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 500
       IXM = 30
 C
 C     READ NUMBER OF OBSERVATIONS AND NUMBER OF UNKNOWN PARAMETERS
 C          INDEPENDENT VARIABLES
 C          DEPENDENT VARIABLES
 C
       READ (5,100) N, NPAR
       READ (5,101) ((XM(I,J), I=1,N), J=1,NPAR)
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL LLS TO PERFORM LINEAR LEAST SQUARES ANALYSIS
 C
       WRITE (IPRT,102)
       CALL LLS (Y, XM, N, IXM, NPAR, RES, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (21F3.0)
   102 FORMAT ('1RESULTS FROM STARPAC',
      *  ' LINEAR LEAST SQUARES SUBROUTINE LLS')
       END


 Data:

    21    4
   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
  80 80 75 62 62 62 62 62 58 58 58 58 58 58 50 50 50 50 50 56 70
  27 27 25 24 22 23 24 24 23 18 18 17 18 19 18 18 19 19 20 20 20
  89 88 90 87 87 87 93 93 87 80 89 88 82 93 89 86 72 79 80 82 91
  42 37 37 28 18 18 19 20 15 14 14 13 11 12  8  7  8  8  9 15 15



                                     <8-12>
1RESULTS FROM STARPAC LINEAR LEAST SQUARES SUBROUTINE LLS
                                                                                                         STARPAC 2.08S (03/15/90)
+***************************************************************
 *  LINEAR LEAST SQUARES ESTIMATION WITH USER-SPECIFIED MODEL  *
 ***************************************************************


 RESULTS FROM LEAST SQUARES FIT
 -------------------------------

                                                     DEPENDENT       PREDICTED      STD DEV OF                         STD
  ROW             PREDICTOR VALUES                    VARIABLE         VALUE        PRED VALUE        RESIDUAL         RES

    1  80.000000      27.000000      89.000000       42.000000       38.765363       1.7810630       3.2346372        1.19
    2  80.000000      27.000000      88.000000       37.000000       38.917485       1.8285238      -1.9174853        -.72
    3  75.000000      25.000000      90.000000       37.000000       32.444467       1.3553032       4.5555330        1.55
    4  62.000000      24.000000      87.000000       28.000000       22.302226       1.1626690       5.6977742        1.88
    5  62.000000      22.000000      87.000000       18.000000       19.711654       .74116599      -1.7116536        -.54
    6  62.000000      23.000000      87.000000       18.000000       21.006940       .90284068      -3.0069397        -.97
    7  62.000000      24.000000      93.000000       19.000000       21.389491       1.5186314      -2.3894907        -.83
    8  62.000000      24.000000      93.000000       20.000000       21.389491       1.5186314      -1.3894907        -.48
    9  58.000000      23.000000      87.000000       15.000000       18.144379       1.2143513      -3.1443789       -1.05
   10  58.000000      18.000000      80.000000       14.000000       12.732806       1.4506369       1.2671941         .44
   11  58.000000      18.000000      89.000000       14.000000       11.363703       1.2770503       2.6362968         .88
   12  58.000000      17.000000      88.000000       13.000000       10.220540       1.5114771       2.7794604         .97
   13  58.000000      18.000000      82.000000       11.000000       12.428561       1.2872987      -1.4285609        -.48
   14  58.000000      19.000000      93.000000       12.000000       12.050499       1.4714398      -.50499291E-01    -.02
   15  50.000000      18.000000      89.000000       8.0000000       5.6385816       1.4154780       2.3614184         .81
   16  50.000000      18.000000      86.000000       7.0000000       6.0949492       1.1742308       .90505080         .30
   17  50.000000      19.000000      72.000000       8.0000000       9.5199506       2.0821373      -1.5199506        -.61
   18  50.000000      19.000000      79.000000       8.0000000       8.4550930       1.2997464      -.45509295        -.15
   19  50.000000      20.000000      80.000000       9.0000000       9.5982566       1.3549990      -.59825656        -.20
   20  56.000000      20.000000      82.000000       15.000000       13.587853       .91842680       1.4121473         .45
   21  70.000000      20.000000      91.000000       15.000000       22.237713       1.7300647      -7.2377129       -2.64


























                                                              <8-13>
1                                                                                                        STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH USER-SPECIFIED MODEL, CONTINUED

                     STD RES VS ROW NUMBER                                     STD RES VS PREDICTED VALUES
  3.75++---------+---------+----+----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  2.25+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -        *                                          -        -                         *                         -
      -     *                                             -        -                                        *          -
      -*                                                  -        -                                                  *-
   .75+                         *  *      *               +     .75+*      * *                                         +
      -                                                *  -        -            *                                      -
      -                       *              *            -        - *         *                                       -
      -                                 *                 -        -          *                                        -
      -                                           * *     -        -    * *                                            -
  -.75+   *      *       *           *         *          +    -.75+      *   *          *  *                         *+
      -             * *    *                              -        -                   *   **                          -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
 -2.25+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                  *-        -                         *                         -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
 -3.75++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
      1.0                     11.0                     21.0      5.639                    22.28                38.92

             AUTOCORRELATION FUNCTION OF RESIDUALS                        NORMAL PROBABILITY PLOT OF STD RES
     1++---------+---------+----***--+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                       ***                         -        -                                                   -
      -                         *                         -        -                                                   -
      -                         *                         -        -                                                   -
      -                   *******                         -        -                                                   -
     6+                   *******                         +    2.25+                                                   +
      -                         **                        -        -                                                   -
      -                         ****                      -        -                                            *      -
      -                       ***                         -        -                                       *           -
      -                      ****                         -        -                                    *              -
    11+                         ***                       +     .75+                               ** *                +
      -                         ******                    -        -                             *                     -
      -                         *                         -        -                           **                      -
      -                         **                        -        -                          *                        -
      -                         ***                       -        -                        **                         -
    16+                         **                        +    -.75+                  ** ***                           +
      -                     *****                         -        -           *  * *                                  -
      -                    ******                         -        -                                                   -
      -                         ***                       -        -                                                   -
      -                      ****                         -        -                                                   -
    21+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -      *                                            -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    26++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
    -1.00                      0.0                     1.00     -2.5                       0.0                      2.5
                                                              <8-14>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH USER-SPECIFIED MODEL, CONTINUED




                                                  ANALYSIS OF VARIANCE
                        -DEPENDENT ON ORDER VARIABLES ARE ENTERED, UNLESS VECTORS ARE ORTHOGONAL-

  PAR     SUM OF SQUARES                                                               ------ PAR=0 ------    ------ PARS=0 -----
 INDEX    RED DUE TO PAR       CUM MS RED      DF(MSRED)      CUM RES MS      DF(RMS)     F        PROB(F)       F        PROB(F)

   1        6448.76190         6448.76190          1          103.461905        20     613.035       .000     198.185       .000
   2        1750.12199         4099.44195          2          16.7955845        19     166.371       .000     59.9022       .000
   3        130.320772         2776.40156          3          10.4886297        18     12.3886       .003     6.66797       .007
   4        9.96537226         2084.79251          4          10.5194095        17     .947332       .344     .947332       .344

 RESIDUAL     178.8300                            17
 TOTAL        8518.000                            21








































                                                              <8-15>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH USER-SPECIFIED MODEL, CONTINUED

 VARIANCE-COVARIANCE AND CORRELATION MATRICES OF THE ESTIMATED PARAMETERS
 ------------------------------------------------------------------------

    - COVARIANCES ARE ABOVE THE DIAGONAL
    - VARIANCES ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

    COLUMN         1                2                3                4

         1     141.51474        .28758711       -.65179437       -1.6763208
         2     .17926325        .18186730E-01   -.36510675E-01   -.71435215E-02
         3    -.14887895       -.73564128        .13544186        .10476827E-04
         4    -.90159992       -.33891642        .18214234E-03    .24427828E-01




 ------------------------- ESTIMATES FROM FIT ------------------------
+                                                                       ---- ESTIMATES FROM FIT OMITTING LAST PREDICTOR VALUE ----

  ESTIMATED PARAMETER       SD OF PAR     T(PAR=0)   PROB(T)  ACC DIG*
+                                                                       ESTIMATED PARAMETER       SD OF PAR     T(PAR=0)   PROB(T)

   1   -39.9196744         11.8959969      -3.356      .004     14.1          -50.3588401         5.13832806      -9.801      .000
   2    .715640200         .134858185       5.307      .000     14.1           .671154441         .126691047       5.298      .000
   3    1.29528612         .368024265       3.520      .003     13.8           1.29535137         .367485444       3.525      .002
   4   -.152122519         .156294043      -.9733      .344     14.1


 RESIDUAL STANDARD DEVIATION               3.243364                                                               3.238615
 BASED ON DEGREES OF FREEDOM         21 -  4 =   17
+                                                                                                          21 -  3 =   18

 MULTIPLE CORRELATION COEFFICIENT SQUARED     .9136

 APPROXIMATE CONDITION NUMBER             1047.370


 * THE NUMBER OF CORRECTLY COMPUTED DIGITS IN EACH PARAMETER USUALLY DIFFERS BY LESS THAN 1 FROM THE VALUE GIVEN HERE.

















                                                              <8-16>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE LLSP USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, X AND RES MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(30), X(30), RES(30)
       DOUBLE PRECISION DSTAK(500)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 500
 C
 C     READ NUMBER OF OBSERVATIONS AND DEGREE OF THE POLYNOMIAL TO BE FIT
 C          INDEPENDENT AND DEPENDENT VARIABLES
 C
       READ (5,100) N, NDEG
       READ (5,101) (X(I), I=1,N)
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL LLSP TO PERFORM LINEAR LEAST SQUARES ANALYSIS
 C
       WRITE (IPRT,102)
       CALL LLSP (Y, X, N, NDEG, RES, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (10F5.1)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' LINEAR LEAST SQUARES SUBROUTINE LLSP')
       END


 Data:

     9    2
   0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0
  12.0 10.5 10.0  8.0  7.0  8.0  7.5  8.5  9.0








                                     <8-17>
1RESULTS OF STARPAC LINEAR LEAST SQUARES SUBROUTINE LLSP
                                                                                                         STARPAC 2.08S (03/15/90)
+***********************************************************
 *  LINEAR LEAST SQUARES ESTIMATION WITH POLYNOMIAL MODEL  *
 ***********************************************************


 RESULTS FROM LEAST SQUARES FIT
 -------------------------------

                                                     DEPENDENT       PREDICTED      STD DEV OF                         STD
  ROW             PREDICTOR VALUES                    VARIABLE         VALUE        PRED VALUE        RESIDUAL         RES

    1                0.                              12.000000       12.184848       .41999992      -.18484848        -.61
    2                 1.0000000                      10.500000       10.521212       .27284429      -.21212121E-01    -.05
    3                 2.0000000                      10.000000       9.2233766       .23159593       .77662338        1.68
    4                 3.0000000                      8.0000000       8.2913420       .24891687      -.29134199        -.64
    5                 4.0000000                      7.0000000       7.7251082       .26115476      -.72510823       -1.63
    6                 5.0000000                      8.0000000       7.5246753       .24891687       .47532468        1.05
    7                 6.0000000                      7.5000000       7.6900433       .23159593      -.19004329        -.41
    8                 7.0000000                      8.5000000       8.2212121       .27284429       .27878788         .64
    9                 8.0000000                      9.0000000       9.1181818       .41999992      -.11818182        -.39






































                                                              <8-18>
1                                                                                                        STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH POLYNOMIAL MODEL, CONTINUED

                     STD RES VS ROW NUMBER                                     STD RES VS PREDICTED VALUES
  3.75++---------+---------+----+----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  2.25+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -             *                                     -        -                  *                                -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
   .75+                               *                   +     .75+*                                                  +
      -                                            *      -        -       *                                           -
      -                                                   -        -                                                   -
      -      *                                            -        -                                *                  -
      -                                      *           *-        -  *              *                                 -
  -.75+*                  *                               +    -.75+        *                                         *+
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                         *                         -        -  *                                                -
      -                                                   -        -                                                   -
 -2.25+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
 -3.75++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
      1.0                      5.0                      9.0      7.525                    9.855                12.18

             AUTOCORRELATION FUNCTION OF RESIDUALS                        NORMAL PROBABILITY PLOT OF STD RES
     1++---------+------*********----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                 *********                         -        -                                                   -
      -                         ****                      -        -                                                   -
      -                         *                         -        -                                                   -
      -                         ***                       -        -                                                   -
     6+                         *                         +    2.25+                                                   +
      -                         *                         -        -                                                   -
      -                         *                         -        -                                        *          -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    11+                                                   +     .75+                                  *                +
      -                                                   -        -                              *                    -
      -                                                   -        -                                                   -
      -                                                   -        -                           *                       -
      -                                                   -        -                       * *                         -
    16+                                                   +    -.75+                *   *                              +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -           *                                       -
      -                                                   -        -                                                   -
    21+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    26++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
    -1.00                      0.0                     1.00     -2.5                       0.0                      2.5
                                                              <8-19>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH POLYNOMIAL MODEL, CONTINUED




                                                  ANALYSIS OF VARIANCE
                        -DEPENDENT ON ORDER VARIABLES ARE ENTERED, UNLESS VECTORS ARE ORTHOGONAL-

  PAR     SUM OF SQUARES                                                               ------ PAR=0 ------    ------ PARS=0 -----
 INDEX    RED DUE TO PAR       CUM MS RED      DF(MSRED)      CUM RES MS      DF(RMS)     F        PROB(F)       F        PROB(F)

   1        720.027778         720.027778          1          2.59027778         8     2696.46       .000     922.687       .000
   2        8.81666667         364.422222          2          1.70079365         7     33.0178       .001     35.8017       .000
   3        10.3033911         246.382612          3          .267027417         6     38.5855       .001     38.5855       .001

 RESIDUAL     1.602165                             6
 TOTAL        740.7500                             9









































                                                              <8-20>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+LINEAR LEAST SQUARES ESTIMATION WITH POLYNOMIAL MODEL, CONTINUED

 VARIANCE-COVARIANCE AND CORRELATION MATRICES OF THE ESTIMATED PARAMETERS
 ------------------------------------------------------------------------

    - COVARIANCES ARE ABOVE THE DIAGONAL
    - VARIANCES ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

    COLUMN         1                2                3

         1     .17639993       -.82535747E-01    .80917399E-02
         2    -.80268762        .59936673E-01   -.69357771E-02
         3     .65431992       -.96215765        .86697213E-03




 ------------------------- ESTIMATES FROM FIT ------------------------
+                                                                       ---- ESTIMATES FROM FIT OMITTING LAST PREDICTOR VALUE ----

  ESTIMATED PARAMETER       SD OF PAR     T(PAR=0)   PROB(T)  ACC DIG*
+                                                                       ESTIMATED PARAMETER       SD OF PAR     T(PAR=0)   PROB(T)

   1    12.1848485         .419999917       29.01      .000     14.1           10.4777778         .801574729       13.07      .000
   2   -1.84653680         .244819675      -7.542      .000     14.1          -.383333333         .168364369      -2.277      .057
   3    .182900433         .294443905E-01   6.212      .001     14.0


 RESIDUAL STANDARD DEVIATION               .5167470                                                               1.304145
 BASED ON DEGREES OF FREEDOM          9 -  3 =    6
+                                                                                                           9 -  2 =    7

 MULTIPLE CORRELATION COEFFICIENT SQUARED     .9227

 APPROXIMATE CONDITION NUMBER             118.0340


 * THE NUMBER OF CORRECTLY COMPUTED DIGITS IN EACH PARAMETER USUALLY DIFFERS BY LESS THAN 1 FROM THE VALUE GIVEN HERE.



















                                                              <8-21>
1G.  Acknowledgments

      The code and printed output for the linear least squares  subroutines has
 been  modeled on the linear least  squares code and output used  by OMNITAB II
 [Hogben et al., 1971].






















































                                     <8-22>
1-----                             CHAPTER 9                              -----

                            NONLINEAR LEAST SQUARES


 A.  Introduction

      STARPAC contains 16 user-callable subroutines for nonlinear least squares
 regression.  Twelve of these are estimation subroutines that compute the least
 squares solution as described below, performing either weighted  or unweighted
 regression  with either  numerically approximated or  user-supplied (analytic)
 derivatives.   The estimation subroutines allow three levels of control of the
 computations and printed output, and allow the user to specify a subset of the
 parameters to be treated as constants,  with their values held fixed  at their
 input  values.   This last  feature allows  the user  to  examine the  results
 obtained  estimating various  subsets  of the  parameters of  a  general model
 without rewriting  the  model subroutine  for each  subset.   The  other  four
 subroutines described  in  this chapter  are utility  procedures  which choose
 optimum step  sizes  for numerically  approximating the  derivative  and which
 verify the correctness of user-supplied (analytic) derivatives.

      Each  of the  subroutines  described in  this chapter  assumes  that  the
 observations of the dependent variable, y(i), are modeled by

                 y(i) = f(x(i),PAR) + e(i)  for i = 1, ..., N,

 where

 N       is the number of observations;

 f       is  the function (nonlinear  in its  parameters) that  models  the ith
         observation;

 x(i)    is the vector of the M independent variables at the ith observation;

 PAR     is the vector of the NPAR model parameters; and

 e(i)    is  the unobservable random  error in  the ith  observation,  which is
         estimated by the ith residual.

      The least squares estimates of the parameters, PAR, are obtained using an
 iterative procedure that requires the matrix of  partial  derivatives  of  the
 model with respect to each parameter,

                    D(i,k) = partial [ f(x(i),PAR) wrt PAR(k) ]

 for i = 1, ..., N and k = 1, ..., NPAR.

 The  derivative   matrix  may   be  supplied   analytically   or  approximated
 numerically.

       The least squares solution is that which minimizes (with respect to PAR)
 the residual sum of squares function,

                    N                                   N
        RSS(PAR) = SUM wt(i)*(y(i) - f(x(i),PAR))**2 = SUM wt(i)*e(i)**2
                   i=1                                 i=1


                                     <9-1>
1here

 wt(i)   is the weight assigned to the ith observation  (wt(i)  =  1.0  in  the
         ''unweighted'' case). Appendix B discusses several common applications
         for weighted least squares.

      The user must supply both initial  values  for  the  parameters  and  the
 subroutine NLSMDL (described in section D) used to compute f(x(i),PAR), i = 1,
 ...,  N,  i.e.,  the predicted values of  the  dependent  variable  given  the
 independent  variables  and the parameter values from each iteration.  Initial
 parameter  values  should  be  chosen  with  care,   since  good  values   can
 significantly reduce computing time.

      STARPAC  provides  a variety of subroutines to accommodate many levels of
 user sophistication and problem difficulty.  Users are directed to  section  B
 for  a  brief  description  of  the  subroutines.  The  declaration  and  CALL
 statements are given in section C, and the subroutine arguments are defined in
 section D.  The algorithms used and the output produced by  these  subroutines
 are  discussed  in  section  E.  Sample programs and their output are shown in
 section F.


 B.  Subroutine Descriptions

 B.1  Nonlinear Least Squares Estimation Subroutines

      The simplest of  the 12  nonlinear least squares  estimation subroutines,
 NLS, requires neither  user-supplied weights  nor analytic  derivatives.   The
 estimated results and a variety of statistics are automatically  summarized in
 a  five-part printed report,  and the  estimated parameters and  residuals are
 returned  to the user  via the  subroutine argument  list  (level one control,
 described  below).  Most nonlinear least squares problems  can be solved using
 NLS.

      The other 11  estimation subroutines  add the  weighting,  derivative and
 level two and three control features both singly and in combination, providing
 greater  flexibility to  the user  at  the price  of less  simplicity.   These
 features are indicated by the  suffix letter(s) on the subroutine  name (e.g.,
 NLSS and NLSWDC).

         * Suffix W indicates user-supplied weights.

         * Suffix D indicates user-supplied (analytic) derivatives.

         * Suffix C indicates level two control of the computations.

         * Suffix S indicates level three control of the computations.

      The  three levels  of  computation  and  printed output  control  are  as
 follows.

         * In level one,  a five-part printed report,  discussed in  detail  in
           section  E.2.a,  is  automatically  provided and the estimated model
           parameters and residuals are returned to the user via  the  argument
           list.

         * Level two also returns the estimated parameters and residuals,  and,
           in addition, allows the user to supply arguments to indicate

                                     <9-2>
1          - a subset of the model parameters to be treated as constants,  with
             their values held fixed at their input values;
           - either the step sizes used to compute the numerical approximations
             to the derivative, or, when user-supplied analytic derivatives are
             used, whether they will be checked;
           - the maximum number of iterations allowed;
           - the convergence criteria;
           - the scale (i.e., the typical size) of each parameter;
           - the  maximum  change  allowed  in  the  parameters  at  the  first
             iteration;
           - how the variance-covariance matrix is to be approximated; and
           - the amount of printed output desired.

         * Level three  has all the  features of  level two,  and,  in addition
           returns the following estimated values via the argument list:
           - the  number of nonzero weighted observations (only when a weighted
             analysis is performed);
           - the number of parameters actually estimated;
           - the residual standard deviation;
           - the predicted values;
           - the standard deviations of the predicted values;
           - the standardized residuals; and
           - the variance-covariance matrix of the estimated parameters.


 B.2  Derivative Step Size Selection Subroutines

      When the partial derivatives used in the nonlinear least squares solution
 are not available analytically,  STARPAC subroutines approximate them  
 numerically. In this case, the subroutines can  select  optimum step sizes for
 approximating the derivatives [see section  E.1.b].  The  user  also  has  the
 option  of  computing these step sizes independently of the estimation process
 by calling either of the two step size  selection  subroutines  directly.  For
 example,  when  planning  to  use  the  parameter  fixing capability [argument
 IFIXED] to examine several subsets of  the  parameters  of  a  general  model,
 computing  the  step sizes first and passing them to the estimation subroutine
 is more efficient than recomputing them each time the estimation subroutine is
 called.

      The simplest of  the two  user-callable step size  selection subroutines,
 STPLS, summarizes the step size selection information for each parameter  in a
 printed  report and  returns the step  sizes to  the user  via  the subroutine
 argument list.

      The second step  size selection  subroutine, STPLSC,  differs  from STPLS
 only in that it enables the user to supply arguments to indicate
          - the number of reliable digits in the model results;
          - the number  of  exemptions  allowed  by  the  acceptance  criteria,
            specified  as a proportion of the total number of observations (see
            section E.1.b);
          - the scale (i.e., the typical size) of each parameter; and
          - the amount of printed output desired.


 B.3  Derivative Checking Subroutines

      When the partial derivatives used in the nonlinear least squares solution
 are available analytically,  the user can code them for use by the  estimation

                                     <9-3>
1subroutines  [see  section  D,  argument NLSDRV].  Because coding errors are a
 common  problem  with  user-supplied  derivatives,   the  STARPAC   estimation
 subroutines  automatically  check the validity of the user-supplied derivative
 code by comparing its results  to  numerically  approximated  values  for  the
 derivative.   When  the  results  are  questionable,  the  checking  procedure
 attempts to determine whether the problem lies with the user's  code  or  with
 the accuracy of the numerical approximation [see section E.1.c].  Although the
 checking  procedure  is  automatically available to the estimation subroutines
 which accept user-supplied  derivatives,  the  user  may  want  to  check  the
 derivative code independently of the estimation process.  In these cases,  the
 user can call either of the two derivative checking subroutines directly,  and
 suppress checking by the  estimation  subroutines  [see  section  D,  argument
 IDRVCK].

      The  simplest   of  the  two  derivative  checking   subroutines,  DCKLS,
 summarizes the results of the check in a printed report.

      The second of  the derivative  checking subroutine, DCKLSC,  differs from
 DCKLS only in that it enables the user to supply arguments to indicate
          - the number of reliable digits in the model results;
          - the agreement tolerance;
          - the scale (i.e., the typical size) of each parameter;
          - the row at which the derivative is to be checked; and
          - the amount of printed output desired.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.


                 Nonlinear Least Squares Estimation Subroutines

      The <basic declaration block> identifies declaration statements  that are
 needed by all of the nonlinear least squares estimation subroutines.  The user
 should substitute the following four statements for each occurrence  of <basic
 declaration block> given below.

                    <real> Y(n), XM(n,m), PAR(npar), RES(n)
                    DOUBLE PRECISION DSTAK(ldstak)
                    COMMON /CSTAK/ DSTAK
                    EXTERNAL NLSMDL

                                      ===

 NLS:    Compute  and  print  a  five-part   weighted nonlinear  least  squares
         analysis with  numerically approximated derivatives;  return parameter
         estimates and residuals

         <basic declaration block>
         :
         :
         CALL NLS (Y, XM, N, M, IXM, NLSMDL,
        +          PAR, NPAR, RES, LDSTAK)

                                      ===

                                     <9-4>
1
 NLSC:   Compute  and  optionally print a five-part  unweighted nonlinear least
         squares  analysis  with  numerically  approximated  derivatives  using
         user-supplied  control   values;  return   parameter   estimates   and
         residuals

         <basic declaration block>
         INTEGER IFIXED(npar)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         :
         :
         CALL NLSC (Y, XM, N, M, IXM, NLSMDL,
        +           PAR, NPAR, RES, LDSTAK,
        +           IFIXED, STP, MIT, STOPSS, STOPP,
        +           SCALE, DELTA, IVAPRX, NPRT)

                                      ===

 NLSS:   Compute  and optionally  print a five-part  unweighted nonlinear least
         squares  analysis  with  numerically  approximated  derivatives  using
         user-supplied control  values; return parameter  estimates, residuals,
         number of  nonzero  weights, number of parameters  estimated, residual
         standard  deviation,  predicted  values, standard  deviations  of  the
         predicted  values  and  variance-covariance matrix  of  the  estimated
         parameters

         <basic declaration block>
         INTEGER IFIXED(npar)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         <real> RSD, PV(n), SDPV(n), SDRES(n), VCV(npare,npare)
         :
         :
         CALL NLSS (Y, XM, N, M, IXM, NLSMDL,
        +           PAR, NPAR, RES, LDSTAK,
        +           IFIXED, STP, MIT, STOPSS, STOPP,
        +           SCALE, DELTA, IVAPRX, NPRT,
        +           NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 NLSW:   Compute  and  print  a  five-part   weighted nonlinear  least  squares
         analysis with  numerically approximated derivatives;  return parameter
         estimates and residuals

         <basic declaration block>
         <real> WT(n)
         :
         :
         CALL NLSW (Y, WT, XM, N, M, IXM, NLSMDL,
        +           PAR, NPAR, RES, LDSTAK)

                                      ===







                                     <9-5>
1NLSWC:  Compute  and  optionally  print a five-part  weighted nonlinear  least
         squares  analysis  with  numerically  approximated  derivatives  using
         user-supplied control values; return parameter estimates and residuals

         <basic declaration block>
         INTEGER IFIXED(npar)
         <real> WT(n)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         :
         :
         CALL NLSWC (Y, WT, XM, N, M, IXM, NLSMDL,
        +            PAR, NPAR, RES, LDSTAK,
        +            IFIXED, STP, MIT, STOPSS, STOPP,
        +            SCALE, DELTA, IVAPRX, NPRT)

                                      ===

 NLSWS:  Compute  and  optionally  print a five-part  weighted nonlinear  least
         squares  analysis  with  numerically  approximated  derivatives  using
         user-supplied control  values; return parameter  estimates, residuals,
         number of  nonzero weights,  number of parameters  estimated, residual
         standard  deviation,  predicted  values, standard  deviations  of  the
         predicted  values  and  variance-covariance matrix  of  the  estimated
         parameters

         <basic declaration block>
         INTEGER IFIXED(npar)
         <real> WT(n)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         <real> RSD, PV(n), SDPV(n), SDRES(n), VCV(npare,npare)
         :
         :
         CALL NLSWS (Y, WT, XM, N, M, IXM, NLSMDL,
        +            PAR, NPAR, RES, LDSTAK,
        +            IFIXED, STP, MIT, STOPSS, STOPP,
        +            SCALE, DELTA, IVAPRX, NPRT,
        +            NNZW, NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 NLSD:   Compute  and  print  a five-part   unweighted nonlinear least  squares
         analysis with  user-supplied derivatives;  return  parameter estimates
         and residuals

         <basic declaration block>
         EXTERNAL NLSDRV
         :
         :
         CALL NLSD (Y, XM, N, M, IXM, NLSMDL, NLSDRV,
        +           PAR, NPAR, RES, LDSTAK)

                                      ===







                                     <9-6>
1NLSDC:  Compute  and  optionally print a five-part  unweighted nonlinear least
         squares analysis  with user-supplied  derivatives  using user-supplied
         control values; return parameter estimates and residuals

         <basic declaration block>
         EXTERNAL NLSDRV
         INTEGER IFIXED(npar)
         <real> STOPSS, STOPP, SCALE(npar), DELTA
         :
         :
         CALL NLSDC (Y, XM, N, M, IXM, NLSMDL, NLSDRV,
        +            PAR, NPAR, RES, LDSTAK,
        +            IFIXED, IDRVCK, MIT, STOPSS, STOPP,
        +            SCALE, DELTA, IVAPRX, NPRT)

                                      ===

 NLSDS:  Compute  and  optionally print a five-part  unweighted nonlinear least
         squares analysis  with user-supplied  derivatives  using user-supplied
         control  values;  return  parameter estimates,  residuals,  number  of
         parameters estimated,  residual standard deviation,  predicted values,
         standard  deviations of  the predicted values  and variance-covariance
         matrix of the estimated parameters

         <basic declaration block>
         EXTERNAL NLSDRV
         INTEGER IFIXED(npar)
         <real> STOPSS, STOPP, SCALE(npar), DELTA
         <real> RSD, PV(n), SDPV(n), SDRES(n), VCV(npare,npare)
         :
         :
         CALL NLSDS (Y, XM, N, M, IXM, NLSMDL, NLSDRV,
        +            PAR, NPAR, RES, LDSTAK,
        +            IFIXED, IDRVCK, MIT, STOPSS, STOPP,
        +            SCALE, DELTA, IVAPRX, NPRT,
        +            NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===

 NLSWD:  Compute  and  print  a  five-part   weighted nonlinear  least  squares
         analysis with  user-supplied derivatives;  return  parameter estimates
         and residuals

         <basic declaration block>
         EXTERNAL NLSDRV
         <real> WT(n)
         :
         :
         CALL NLSWD (Y, WT, XM, N, M, IXM, NLSMDL, NLSDRV,
        +            PAR, NPAR, RES, LDSTAK)

                                      ===







                                     <9-7>
1NLSWDC: Compute  and  optionally  print a five-part  weighted nonlinear  least
         squares analysis  with user-supplied  derivatives  using user-supplied
         control values; return parameter estimates and residuals

         <basic declaration block>
         EXTERNAL NLSDRV
         INTEGER IFIXED(npar)
         <real> WT(n)
         <real> STOPSS, STOPP, SCALE(npar), DELTA
         :
         :
         CALL NLSWDC (Y, WT, XM, N, M, IXM, NLSMDL, NLSDRV,
        +             PAR, NPAR, RES, LDSTAK,
        +             IFIXED, IDRVCK, MIT, STOPSS, STOPP,
        +             SCALE, DELTA, IVAPRX, NPRT)

                                      ===

 NLSWDS: Compute  and  optionally  print a five-part  weighted nonlinear  least
         squares analysis  with user-supplied  derivatives  using user-supplied
         control  values;  return  parameter estimates,  residuals,  number  of
         nonzero  weights,  number  of parameters estimated,  residual standard
         deviation,  predicted values,  standard deviations  of  the  predicted
         values and variance-covariance matrix of the estimated parameters

         <basic declaration block>
         EXTERNAL NLSDRV
         INTEGER IFIXED(npar)
         <real> WT(n)
         <real> STOPSS, STOPP, SCALE(npar), DELTA
         <real> RSD, PV(n), SDPV(n), SDRES(n), VCV(npare,npare)
         :
         :
         CALL NLSWDS (Y, WT, XM, N, M, IXM, NLSMDL, NLSDRV,
        +             PAR, NPAR, RES, LDSTAK,
        +             IFIXED, IDRVCK, MIT, STOPSS, STOPP,
        +             SCALE, DELTA, IVAPRX, NPRT,
        +             NNZW, NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===


                        Step Size Selection Subroutines

 STPLS:  Compute  and  print  optimum step sizes for numerically  approximating
         derivatives; return selected step sizes

         <real> XM(n,m), PAR(npar), STP(npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         EXTERNAL NLSMDL
         :
         :
         CALL STPLS (XM, N, M, IXM, NLSMDL, PAR, NPAR, LDSTAK, STP)

                                      ===



                                     <9-8>
1STPLSC: Compute  and  optionally  print  optimum  step sizes  for  numerically
         approximating derivatives  using user-supplied control  values; return
         selected step sizes

         <real> XM(n,m), PAR(npar), STP(npar)
         <real> EXMPT, SCALE(npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         EXTERNAL NLSMDL
         :
         :
         CALL STPLSC (XM, N, M, IXM, NLSMDL, PAR, NPAR, LDSTAK, STP,
        +             NETA, EXMPT, SCALE, NPRT)

                                      ===


                        Derivative Checking Subroutines

 DCKLS:  Perform and print derivative checking analysis; return error code

         <real> XM(n,m), PAR(npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         EXTERNAL NLSMDL, NLSDRV
         :
         :
         CALL DCKLS (XM, N, M, IXM, NLSMDL, NLSDRV, PAR, NPAR, LDSTAK)

                                      ===

 DCKLSC: Perform and optionally print derivative checking analysis using
         user-supplied control values; return error code

         <real> XM(n,m), PAR(npar)
         <real> SCALE(npar)
         DOUBLE PRECISION DSTAK(ldstak)
         COMMON /CSTAK/ DSTAK
         EXTERNAL NLSMDL, NLSDRV
         :
         :
         CALL DCKLSC (XM, N, M, IXM, NLSMDL, NLSDRV, PAR, NPAR, LDSTAK,
        +             NETA, NTAU, SCALE, NROW, NPRT)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the argument is input to  the subroutine  and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.

                                     <9-9>
1

 D       <-- The matrix of exact dimension N by NPAR that contains  the partial
             derivatives  of   the  model  with  respect  to   each  parameter,
             PAR(k), k = 1, ..., NPAR.  This argument is used within derivative
             subroutine NLSDRV [see argument NLSDRV below].

 DELTA   --> The maximum scaled change  allowed in the parameters at  the first
             iteration [see section E.1.a].  The default value is  100.0.  When
             DELTA  <=  0.0  or when DELTA is not an argument of the subroutine
             CALL statement the default value  is  used.  A  smaller  value  of
             DELTA  may  be  appropriate  if,   at  the  first  iteration,  the
             computation of the predicted values from the user's model subroutine
	     produces an arithmetric overflow or the parameters leave the
             region of interest in parameter space.  A  reasonable  alternative
             to  the  default  value  of  DELTA is an upper bound to the scaled
             change that the estimated parameters should be allowed to make  on
             the first iteration,

             DELTA = min(|del(PAR(k))|/SCALE(k), for k = 1, ..., NPAR)

             where del(PAR(k)) is  the  maximum  change  allowed  for  the  kth
             parameter at the first iteration.

 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 EXMPT   --> The  proportion  used  to  compute  the  number  of  observations,
             a = EXMPT*N, for which the forward difference  quotient derivative
             with  respect to a  given parameter  is exempted from  meeting the
             acceptance  criteria  for step size selection [see section E.1.b].
             The default value for EXMPT is 0.1 (10 percent).  When  the  
	     user-supplied value is outside the range [0.0, 1.0], or when EXMPT is
             not an argument of the  subroutine  CALL  statement,  the  default
             value is used.

 IDRVCK  --> The indicator  variable  used  to  designate whether  or  not  the
             user-supplied  derivative  subroutine  is  to  be  checked.   When
             IDRVCK  <> 0 the derivative is checked,  and when IDRVCK = 0 it is
             not.  The default value is IDRVCK <> 0.  When  IDRVCK  is  not  an
             argument  of  the  subroutine  CALL statement the default value is
             used.

 IERR    ... An  error  flag  returned  in  COMMON  /ERRCHK/  [see  chapter  1,
             section D.5].  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             For the estimation subroutines:

               IERR = 0 indicates that  no errors  were detected, and  that the
                        iterations converged satisfactorily.

               IERR = 1 indicates  that  improper  input was  detected.



                                     <9-10>
1              IERR = 2 indicates that  the computation of the residual  sum of
                        squares using the initial parameter values  produced an
                        arithmetic overflow.   The user  should reduce the size
                        of DELTA or should supply new starting values.

               IERR = 3 indicates that  the model is  computationally singular,
                        which means the model has too many parameters  near the
                        solution.   The user should examine the model and  data
                        to determine and remove the cause of the singularity.

               IERR = 4 indicates  that  at  least  one  of   the  standardized
                        residuals  could not  be computed because  its standard
                        deviation was  zero.   The validity  of the  covariance
                        matrix is questionable.

               IERR = 5 indicates false convergence [see section E.1.a].

               IERR = 6 indicates  that  convergence  was not  reached  in  the
                        allowed number of iterations or model subroutine  calls
                        [see argument MIT].

               IERR = 7 indicates that the variance-covariance matrix could not
                        be computed.

             For the step size selection subroutines:

               IERR = 0 indicates that  no errors  were detected, and  that all
                        the step sizes satisfied the selection criteria.

               IERR = 1 indicates that improper input was detected.

               IERR = 2 indicates that one or  more of  the step sizes  did not
                        satisfy the selection criteria.

             For the derivative checking subroutines:

               IERR = 0 indicates that  no errors  were detected, and  that the
                        user-supplied derivative code appears to be correct.

               IERR = 1 indicates that improper input was detected.

               IERR = 2 indicates that  the user-supplied  derivative  code and
                        numerical derivatives  do not  agree for  at  least one
                        parameter, but  that in  each case of  disagreement the
                        accuracy of the numerical derivatives  is questionable.
                        Further testing is suggested.

               IERR = 3 indicates that  the user-supplied  derivative  code and
                        numerical derivatives  do not  agree for  at  least one
                        parameter, and in at least one instance of disagreement
                        there is no reason to doubt the numerical derivatives.

 IFIXED  --> The vector of dimension at least NPAR that contains values used to
             indicate whether  the  corresponding parameter  in PAR  is  to  be
             treated  as  a  fixed  constant   or  is  to  be  estimated.    If
             IFIXED(I) > 0, PAR(I) will be  held fixed  at its input  value; if
             IFIXED(I) = 0, PAR(I)  will be  estimated using the  least squares
             procedure  described  in  section  A.   The  default  values   are

                                     <9-11>
1            IFIXED(I)  =  0,  I  =  1,  ...,  NPAR,  i.e.,  all parameters are
             estimated.  When IFIXED(1)  <=  -1,  or  when  IFIXED  is  not  an
             argument of the subroutine CALL statement,  the default value will
             be used.

 IVCV    --> The exact  value of  the  first dimension  of the  matrix  VCV  as
             specified in the calling program.

 IVAPRX  --> The indicator variable used to specify how the variance-covariance
             matrix,  VCV, is  to be  approximated.   Three approximations  are
             available:

             (1) VCV = RSD**2 * inv(trans(Dhat)*W*Dhat)

             (2) VCV = RSD**2 * inv(Hhat)

             (3) VCV = RSD**2 * inv(Hhat) * (trans(Dhat)*W*Dhat) * inv(Hhat)

             where

             trans(.) indicates the transpose of the designated matrix;

             inv(.) indicates the inverse of the designated matrix;

             Hhat is the matrix of second partial derivatives of the model with
               respect to each parameter (the Hessian matrix), evaluated at the
               solution

               = trans(Dhat)*W*Dhat +
                   N
                 (SUM  e(i)*wt(i)*(second partial of e(i) wrt PAR(j) & PAR(k)
                  i=1
                                   for j = 1, ..., NPAR & k = 1, ..., NPAR));

             W is an N by N diagonal matrix of weights,

               W = diag(wt(i), i = 1, ..., N),

             when a weighted analysis is performed, and  is the identity matrix
             otherwise,  and

             Dhat is the matrix that contains the partial  derivatives  of  the
               model  with  respect  to  each  parameter (the Jacobian matrix),
               evaluated at the solution.

             Approximation  (1)  is  based  on  the  assumption   that   H   is
             approximately  equal  to  trans(D)*W*D  because  the residuals are
             sufficiently small at the solution;  approximation (2) is based on
             the  assumption  that  the  necessary  conditions  for  asymptotic
             maximum likelihood theory have been met;  and approximation (3) is
             based   on  the  assumption  that  the  necessary  conditions  for
             asymptotic maximum likelihood theory may be violated.  The results
             of  a  study  by  Donaldson  and  Schnabel  [1987]  indicate  that
             approximation  (1)  is  preferable  because  it  is  simple,  less
             expensive,  more numerically stable and at least  as  accurate  as
             approximations  (2)  and (3).  However,  all approximations to the
             variance-covariance  matrix  are  subject  to  sampling  variation
             because  they  are  computed using the estimated parameter values.

                                     <9-12>
1            The  variance-covariance  matrix  computed  for   any   particular
             nonlinear least squares solution should thus be regarded as only a
             rough estimate [Bard, 1974; Donaldson and Schnabel, 1987].

             If IVAPRX = 1 or 4 then approximation (1) is used;
                       = 2 or 5 then approximation (2) is used; and
                       = 3 or 6 then approximation (3) is used.

             If IVAPRX = 1, 2, or 3,   then,   when    user-supplied   analytic
             derivatives are available [see argument NLSDRV], they are  used to
             compute  VCV; if  IVAPRX = 4,  5, or  6, then  only  the predicted
             values from the model subroutine  are used  to compute VCV.   When
             analytic derivatives  are  available, options  1, 2,  or  3,  will
             generally result in a faster, more accurate computation of VCV.

             The  default  value for  IVAPRX is 1.   When  argument  IVAPRX  is
             outside the range [1, 6], or when IVAPRX is not an argument of the
             subroutine CALL statement, then the default value will be used.

 IXM     --> The exact  value of  the  first  dimension  of the  matrix  XM  as
             specified in the calling program.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision  version of  STARPAC is  being  used P = 0.5,
             otherwise P = 1.0 [see chapter 1, section B].

             For NLS, NLSC, NLSS, NLSW, NLSWC and NLSWS:

               LDSTAK >= 27 + max(IS*(N+NPAR), 30+NPARE) +

                         max(IS*10*N, 94+N*(3+NPAR)+(3*NPARE**2+37*NPARE)/2)*P

               with IS = 1 if default values are used for the  derivative step
               sizes, and IS = 0 otherwise.

             For NLSD, NLSDC, NLSDS, NLSWD, NLSWDC and NLSWDS:

               LDSTAK >= 45 + NPAR + (94+N*(3+NPAR)+(3*NPARE**2+35*NPARE)/2)*P

             For STPLS and STPLSC:

               LDSTAK >= 27 + (N+NPAR) + 10*N*P

             For DCKLS and DCKLSC:

               LDSTAK >= 14 + NPAR + (N*NPAR+N+NPAR)*P

 M       --> The number of  independent variables, i.e., the number  of columns
             of data in XM.

 MIT     --> The  maximum number of iterations allowed.   This argument is also
             used to  compute  the maximum  number of  model  subroutine calls,
             (2*MIT).   The iterations  will stop  if either limit  is reached,
             although,  as a rule,  the maximum  number of  iterations  will be
             reached  first.   The  default  value for  the maximum  number  of
             iterations is 21.  When MIT <= 0 or when MIT is not an argument of
             the subroutine CALL statement the default value will be used.

                                     <9-13>
1
 N       --> The number of observations.

 NETA    --> The number of reliable decimal digits in the predicted values (PV)
             computed  by the user's model subroutine.   The default  value for
             NETA is  experimentally determined  by the procedure  described in
             Appendix C.   The default  value will be used when  NETA is not an
             argument   in   the  subroutine   CALL  statement,  or   when  the
             user-supplied  value  of NETA  is outside  the  range [1, DIGITS],
             where DIGITS is the number of decimal digits carried by the user's
             computer for a  single precision  value when the  single precision
             version of STARPAC is being used  and is the number carried  for a
             double precision value otherwise.

 NLSDRV  *** The name of the user-supplied subroutine that computes the partial
             derivative matrix (Jacobian).   This argument must be listed in an
             EXTERNAL  statement  in  the  program  which  calls   the  STARPAC
             estimation  or  derivative checking subroutine.   The form  of the
             derivative subroutine  argument list  and  dimensioning statements
             must be  exactly as shown  below, although  if there  is  only one
             independent variable (M = 1), XM  may be  declared to be  a vector
             with dimension IXM.

               SUBROUTINE NLSDRV (PAR, NPAR, XM, N, M, IXM, D)
               <real> PAR(NPAR), XM(IXM,M), D(N,NPAR)

               < Computations for D(I,J), I = 1, ..., N and J = 1, ..., NPAR >

               RETURN
               END

 NLSMDL  *** The  name  of  the  user-supplied  subroutine  that  computes  the
             predicted value  of the  dependent variable given  the independent
             variables and the current  values of  the model parameters.   This
             argument must be listed  in an  EXTERNAL statement in  the program
             which calls  the STARPAC  estimation, step size  selection, and/or
             derivative checking subroutines.  The form of the model subroutine
             argument list and dimensioning statements must be exactly as shown
             below, although if there is only one independent variable (M = 1),
             XM may be declared to be a vector with dimension IXM.

               SUBROUTINE NLSMDL (PAR, NPAR, XM, N, M, IXM, PV)
               <real> PAR(NPAR), XM(IXM,M), PV(N)

               < Computations for PV(I), I = 1, ..., N >

               RETURN
               END

 NNZW    <-- The number of observations with nonzero weights.  N.B., this value
             is returned by the estimation subroutines.

 NPAR    --> The number of parameters  in the model, including both  those held
             fixed  at  their  starting  values  and  those  which  are  to  be
             estimated.




                                     <9-14>
1NPARE   <-- The number of  parameters actually estimated, i.e., the  number of
             zero  elements  in  IFIXED.  N.B.,  this  value is returned by the
             estimation subroutines.

 NPRT    --> The argument controlling printed output.

             For the estimation subroutines:

               NPRT is  a five-digit integer,  in which  the value  of  the Ith
               digit (counting from left to  right) is used to control  the Ith
               section of the output.

               If the Ith digit = 0 the output from the Ith section is
                                    suppressed;
                                = 1 the brief form of the Ith section is given;
                                >=2 the full form of the Ith section is given.

               The  default value for NPRT is 11112.  When NPRT <= -1,  or when
               NPRT is not an argument in the subroutine  CALL  statement,  the
               default value will be used.  If the convergence criteria are not
               satisfied the subroutine gives a suitable warning and provides a
               printed  report  even  if  NPRT  =  0.  A full discussion of the
               printed output is given in section E.2.a and  is  summarized  as
               follows.

               Section 1 lists  the  starting  estimates  and  control  values.
                         Brief output and  full output  are the  same  for this
                         section.

               Section 2 reports  the results of the iterations.   Brief output
                         includes  information only  about the  first  and last
                         iteration while full output includes information about
                         all of the iterations.

               Section 3 provides information for each observation based on the
                         final solution.  Brief output includes information for
                         the  first 40 observations while full  output provides
                         the information for all of the data.

               Section 4 is a  set of  four residual plots.   Brief output  and
                         full output are the same for this section.

               Section 5 is  the final  summary of  the  estimated  parameters.
                         Brief  output  does not include printing the 
			 variance-covariance matrix while full output does.

             For the step size selection and derivative checking subroutines:

               If NPRT  = 0 the printed output is suppressed.

               If NPRT <> 0 the printed output is provided.

               When  the acceptance criteria  are not  met a printed  report is
               provided even if NPRT = 0.

 NROW    --> The  row  of   the  independent  variable  matrix  at   which  the
             user-supplied derivative code is to be checked.  The default value
             is the first row with no independent variables equal to zero; when

                                     <9-15>
1            all rows have one or more independent variables equal to zero, row
             one will be  used for the default value.   When the  user-supplied
             value is outside the range [1, N] or when NROW is not  an argument
             of the subroutine CALL statement the default value will be used.

 NTAU    --> The agreement tolerance,  i.e., the number of digits  of agreement
             required between the user-supplied derivatives and the derivatives
             numerically approximated  by the  derivative  checking subroutine.
             The default value is NETA/4.  When the user-supplied value of NTAU
             is outside the range [1, NETA/2]  or when NTAU is not  an argument
             of the subroutine CALL statement the default value will be used.

 PAR     --- The vector of dimension at least NPAR that contains  the parameter
             values.   For all  estimation subroutines it must contain  initial
             values  for the parameters  on input  and will  contain  the final
             values on  return.   For  the  step size  and derivative  checking
             subroutines it  must  contain the  parameter values  at  which the
             operations are to be performed.

 PV      <-- The  vector of dimension  at least  N that contains  the predicted
             values of the dependent variable at the solution,

             PV(i) = f(x(i),PAR)   for i = 1, .., N.

 RES     <-- The vector of dimension at least N that contains the  residuals at
             the solution,

             RES(i) = y(i) - f(x(i),PAR) = e(i) for i = 1, ..., N.

 RSD     <-- The residual standard deviation at the solution,

             RSD = sqrt(RSS(PAR)/(NNZW-NPARE)).

 SCALE   --> The vector of dimension at least NPAR that contains the  scale, or
             typical size,  of each  parameter.   The vector  SCALE is  used to
             normalize the size of each parameter so that

             |PAR(j)/SCALE(j)| approximates |PAR(k)/SCALE(k)|

             for k = 1, ..., NPAR and j = 1, ..., NPAR.

             Values of |SCALE(k)| > |PAR(k)| can be used to increase  the  step
             size  in cases where the model function is known to be insensitive
             to small changes in the value PAR(k).

             For the estimation subroutines:

               The  default  values  for  SCALE  are  selected  by  the  NL2SOL
               algorithm [Dennis et al.,  1981a,b]  and  are  updated  at  each
               iteration.  When SCALE is not an argument in the subroutine CALL
               statement  or when the user-supplied value for SCALE(1) <= 0 the
               default procedure will be used  to  select  scale  values.  When
               SCALE(1) > 0,  values of SCALE(k) <= 0 for k = 2, ..., NPAR will
               be interpreted as an input  error.  User-supplied  scale  values
               may  be either a vector of the typical size of each parameter or
               a vector of ones if the typical  sizes  of  the  parameters  are
               roughly  equal;  user-supplied scale values can sometimes result


                                     <9-16>
1              in reduced computing time since these values are not updated  at
               each iteration.

             For the derivative checking and step size selection subroutines:

               The default values of SCALE are defined for k = 1, ..., NPAR as:

               SCALE(k) = 1.0      if PAR(k) = 0.0

               SCALE(k) = |PAR(k)| otherwise

               where PAR(k) is the input value of the k-th parameter.

               When SCALE is not an argument in the subroutine  CALL  statement
               or  when  the  user-supplied value of |SCALE(k)| <= |PAR(k)| the
               default value for SCALE(k) is used.  When  SCALE(1)  <=  0,  the
               default  values  will be used for SCALE(k),  k = 1,  ...,  NPAR.
               When SCALE(1) > 0, values of SCALE(k) <= 0 for k = 2, ...,  NPAR
               will be interpreted as an input error.

 SDPV    <-- The vector of dimension at least N that contains  an approximation
             to  the  standard  deviation  of  each  predicted  value   at  the
             solution,

             SDPV(i) = the ith diagonal element of sqrt(Dhat*VCV*trans(Dhat))

             for i = 1, ..., N, where

             Dhat(i,j) = partial [ f(x(i),PAR) wrt PAR(j) ]

             for i = 1, ..., N and j = 1, ..., NPAR, evaluated at the solution,
             and trans(Dhat) is the transpose of Dhat.

             This  approximation  is based  on a  linearization of the model in
             the   neighborhood   of   the   solution;   the  validity  of  the
             approximation depends on the  nonlinearity  of  the  model.   This
             approximation  may  be  extremely  inaccurate for a problem with a
             highly nonlinear model.

 SDRES   <-- The vector of dimension at least N that contains  an approximation
             to the standardized residuals at the solution,

             SDRES(i) = RES(i)/sqrt[(RSD**2/WT(i)) - SDPV(i)**2]

             for i = 1, ..., N,  which  is  the  ith residual  divided  by  its
             individual  estimated standard deviation.   This approximation  is
             based on a linearization of  the model in the neighborhood  of the
             solution;  the  validity  of  the  approximation  depends  on  the
             nonlinearity  of the model.   This approximation  may be extremely
             inaccurate for a problem with a highly nonlinear model.

 STOPP   --> The stopping value for  the convergence test based on  the maximum
             scaled relative  change  in  the  parameters at  the  most  recent
             iteration.   The convergence criterion is satisfied if the current
             step is a Newton step and




                                     <9-17>
1              max[|PARc(k)-PARp(k)|/SCALE(k) for k = 1, ..., NPAR]
             --------------------------------------------------------  < STOPP.
             max[(|PARc(k)|+|PARp(k)|)/SCALE(k) for k = 1, ..., NPAR]

             where PARc(k) and PARp(k) indicate the current value and the value
             from the previous iteration,  respectively,  of the kth  parameter
             [see Dennis et  al.  1981a].  This  convergence  test  is  roughly
             equivalent  to  the  test  based on the maximum relative change in
             each parameter as measured by

             max(|PARc(k)-PARp(k)|/|PARp(k)| for k = 1,  ..., NPAR).

             STOPP  is  not  a scale-dependent value;  if its value is 10**(-4)
             then this criteria will be met when the first four digits of  each
             parameter  are the same at two successive iterations regardless of
             the size of the parameter values.

             The default value is approximately 10**(-DIGITS/2),  where  DIGITS
             is the number of decimal digits carried by the user's computer for
             a  single  precision  value  when  the single precision version of
             STARPAC is being used and is  the  number  carried  for  a  double
             precision value otherwise.  When the user-supplied value for STOPP
             is  outside  the  interval  [0.0,  1.0]  or  when  STOPP is not an
             argument of the subroutine CALL statement the default  value  will
             be used.

 STOPSS  --> The stopping value for the convergence test based on the  ratio of
             the  forecasted  change  in   the   residual   sum   of   squares,
             fcst(RSS(PAR)),  to the residual sum of squares from the  previous
             iteration.  The  convergence  criterion  is  satisfied  if certain
             conditions are met and

             fcst(RSS(PAR))/RSS(PARp) < STOPSS,

             where the notation is described in  the  description  of  argument
             STOPP  [see  Dennis  et  al.,  1981a].  This  convergence  test is
             roughly equivalent to the test based on the relative change in the
             residual standard deviation between two iterations as measured  by
             (RSDc - RSDl)/RSDc.  STOPSS is not a scale-dependent value; if its
             value  is  10**(-5)  this criteria will be met when the first five
             digits of the  residual  sum  of  squares  are  the  same  at  two
             successive  iterations  regardless of the size of the residual sum
             of squares.

             The default value is approximately the maximum  of  10**(-10)  and
             10**(-2*DIGITS/3),  where DIGITS is the number of  decimal  digits
             carried  by  the user's computer for a single precision value when
             the single precision version of STARPAC is being used and  is  the
             number  carried  for a double precision value otherwise.  When the
             user-supplied   value   for   STOPSS   is   outside  the  interval
             [10**(-DIGITS),  0.1] or when STOPSS is not  an  argument  of  the
             subroutine CALL statement the default value will be used.

 STP     --- The vector of dimension  at least NPAR that contains  the relative
             step sizes used to approximate the derivative matrix  numerically.
             It  is  input  to the estimation subroutines and returned from the
             step size selection subroutines.  The procedure used to select the
             default values is described in section E.1.b.  For the  estimation

                                     <9-18>
1            subroutines,  when  STP  is not an argument of the subroutine CALL
             statement or when STP(1) <= 0 the default values will be used  for
             all  of the step sizes,  and when STP(1) > 0 values of STP(k) <= 0
             for k = 2, ..., NPAR will be interpreted as an input error.

 VCV     <-- The matrix of dimension at least NPARE by NPARE that  contains the
             variance-covariance matrix of the estimated parameters, 
	     approximated as designated by argument IVAPRX. The parameters which are
             held fixed [see argument IFIXED] are not included in the 
	     variance-covariance matrix.

             The approximation of the variance-covariance matrix is based  on a
             linearization of the  model in  the neighborhood of  the solution;
             the validity of  the approximation depends on the  nonlinearity of
             the model.   This approximation may be extremely inaccurate for  a
             problem with a highly nonlinear model.

 WT      --> The  vector of dimension  at least  N that  contains  the weights.
             Negative weights are not allowed and the number of nonzero weights
             must equal or exceed the number of parameters being estimated.   A
             zero  weight eliminates  the corresponding  observation  from  the
             analysis,  although  the residual,  the predicted  value  and  the
             standard  deviation  of  the  predicted value  are  still computed
             [see Appendix B].

 XM      --> The matrix of dimension at least N by M whose jth  column contains
             the N values of the jth independent variable, j = 1, ..., M.

 Y       --> The  vector of dimension  at least  N that contains  the dependent
             variable.


 E.  Computational Methods

 E.1  Algorithms

 E.1.a  Nonlinear Least Squares Estimation

      The  nonlinear  least  squares  estimation  subroutines  use  the  NL2SOL
 software package written by Dennis et al., [1981a,b].  The observations of the
 dependent  variable,  which are measured with error,  are iteratively fit to a
 nonlinear model by minimizing the sums of squares of the errors  as  described
 in section A.  The iterations continue until the convergence criteria based on
 the  change  in  the  parameter  values  or in the residual sum of squares are
 satisfied [see section D,  arguments STOPP and STOPSS],  the maximum number of
 iterations  (or  model  subroutine calls) is reached [see section D,  argument
 MIT],  or the iterations are terminated due to singularity  in  the  model  or
 false  convergence.  All  but  the  first  of  these  stopping  conditions may
 indicate  computational problems and will produce an error report [see chapter
 1, section D.5].

      Singular convergence means  that the model contains too  many parameters,
 at least near the solution,  while false convergence can indicate  that either
 STOPSS or STOPP is set too small for the  accuracy to which the model  and its
 derivatives are being computed or that  there is an error or  discontinuity in
 the derivative.   Users should  examine their models to determine  and correct
 the underlying cause of singular or false convergence.


                                     <9-19>
1     Iterative procedures  for solving  nonlinear least  squares  problems are
 discussed in Dennis and Schnabel  [1983], Draper and Smith [1981]  and Kennedy
 and Gentle [1980].   The specific procedure used in STARPAC is as follows.  At
 the current iteration the values of the parameter vector PARc are given by

       PARc = PARp - inv(trans(Dp)*W*Dp + Sp + Gp)*trans(Dp)*W*trans(ep)

 subject to the restriction that

                   NPAR
            sqrt ( SUM   [(PARc(k) - PARp(k))/SCALE(k)]**2 ) <= dp,
                   k=1

 where

 trans(.) is the transpose of the designated matrix.

 PARp    is the vector of the NPAR estimated parameter values from the previous
         iteration.

 Dp      is the N by NPAR matrix of the partial derivatives  evaluated at PARp,

         D(i,k) = partial[ f(x(i),PAR) wrt PAR(k) ]

         for i = 1, ..., N and k = 1, ..., NPAR.

 W       is an N by N diagonal matrix of user-supplied weights,

         W = diag(wt(i), i = 1, ..., N)

         when  a weighted  analysis  is performed  and is  the  identity matrix
         otherwise.

 Sp      is  an APPROXIMATION  to the exact  term Sp*  in the matrix  of second
         order terms (Hessian) of  the Taylor series expansion of  the residual
         sum of squares function,

                     N
         Sp*(j,k) = SUM [ep(i)*wt(i)*
                    i=1

                         (second partial ep(i) wrt PARp(j) & PARp(k))],

         for j = 1, ..., NPAR and k = 1, ..., NPAR.

 ep      is the vector of the N residuals from the previous iteration.

 dp      is  the adaptively  chosen  diameter of  the trust  region,  i.e., the
         region in which the  local approximation to the user's  model function
         is reliable.  At each iteration,  dp is computed based on  information
         from  the  previous  iteration.  At  the first iteration,  the initial
         value, d0,  is supplied by argument DELTA which can be used to control
         the change in the parameters permitted at the first iteration.

 Gp      is an NPAR by NPAR diagonal matrix,

         Gp = diag(gp/SCALE(k), k = 1, ..., NPAR),


                                     <9-20>
1        where gp is  chosen to  approximate the  smallest  non-negative number
         such that the restriction given above on the size of the change in the
         parameters is satisfied.

      The second order term Sp*,  which is expensive and difficult  to  compute
 accurately,   is   important  only  if  it  is  large  compared  to  the  term
 trans(Dp)*W*Dp,  that is,  when the residuals are large or the model is highly
 nonlinear.  When  Sp*  is  large compared to trans(Dp)*W*Dp,  algorithms which
 ignore it,  such as Levenberg-Marquardt or Gauss-Newton,  may converge slowly.
 The  NL2SOL  algorithm  used  by  STARPAC,  however,  adaptively  decides when
 inclusion of  this  term  is  necessary  for  reliable  results  and  uses  an
 inexpensive approximation to Sp* in those cases.

      The matrix,  D,  of partial derivatives of the model with respect to each
 parameter is either computed analytically using  a  user-supplied  subroutine,
 NLSDRV,  or  is numerically approximated using forward difference quotients as
 described  in  section  E.1.b.   When   the   derivatives   are   approximated
 numerically,  the  least squares solution,  especially the variance-covariance
 matrix,  can be sensitive to the step sizes used for  the  approximation.  The
 user may want to use STARPAC subroutines STPLS or STPLSC to recompute the step
 sizes  at  the  solution provided by the estimation subroutines to assure that
 the  step  sizes  which  were  used  are  still  acceptable.  If  there  is  a
 significant  change  in  the  step  size  the least squares solution should be
 recomputed with the new step sizes from the current  point.  In  addition,  if
 the  estimation  subroutine  has  convergence  problems  the  user may want to
 recompute the step sizes with the most recent parameter values  to  see  if  a
 change  in the curvature of the model,  which will be reflected as a change in
 the optimum step sizes, is causing the problem.

      Dennis et al.   [1981a] provides a detailed description of the  algorithm
 used  in STARPAC.  STARPAC also includes the subroutines NL2SOL,  NL2SNO,  and
 NL2ITR, which they reference, and which can be used as documented by them [see
 Dennis et al., 1981b].


 E.1.b  Derivative Step Size Selection

      The STARPAC step size selection subroutines use an algorithm developed by
 Schnabel [1982] to compute  optimum step  sizes for approximating  the partial
 derivatives  of the  model  with  respect  to each  parameter.   Briefly,  the
 relative  step sizes selected  by these  subroutines are  those  which produce
 forward difference quotient approximations to the derivative, Dfd,  that agree
 reasonably well with the central difference quotient approximations, Dcd.  The
 central difference  quotient  approximations are  twice as  accurate  but also
 twice as expensive to compute.   Since the  additional accuracy is not usually
 needed, central  difference  quotient  approximations  are  not  used  by  the
 estimation subroutines.

      The number of reliable digits in  these derivatives is a function  of the
 step sizes used to compute them.  Given properly chosen step sizes, the number
 of  reliable  digits  in  Dfd  and  Dcd  will  be  approximately  h/2  and  h,
 respectively, where  h is  the  number of  reliable digits  in  the  predicted
 values, PV,  from the user's model subroutine.   For example, if the predicted
 values are  computed using an  iterative procedure  (such as  quadrature  or a
 solution of partial differential equations) which is expected to  provide five
 good digits, then h would be five; if the predicted values are calculated from
 a simple algebraic expression  translated directly  into Fortran code,  then h


                                     <9-21>
1would (usually) be the number of decimal digits carried by the user's computer
 for the results.

      The relative step size for PAR(k), k = 1, ..., NPAR, is initially

            STP(k) = 2*sqrt*(10**(-NETA)/q)   for k = 1, ..., NPAR,

 where

 q       is  the average curvature  (estimated by  STARPAC) of  the  model with
         respect to PAR(k).

 The forward difference quotient approximations with respect to PAR(k),  k = 1,
 ..., NPAR are then

                       f(x(i),PARk) - f(x(i),PAR)
           Dfd(i,k) = ----------------------------   for i = 1, ..., N,
                      STP(k)*SCALE(k)*SIGN(PAR(k))

 where

 f       is the function which models the ith observation.

 x(i)    is the vector of the values of the M independent variables at  the ith
         observation.

 PAR     is the vector of the NPAR parameter values.

 PARk    is a vector which has the same values  as  PAR  except  that  the  kth
         parameter is equal to

         PAR(k) + STP(k)*SCALE(k)*SIGN(PAR(k)).

 SIGN    is a function which returns the sign of its argument.

      The  central  difference  approximations  to the  model  derivative  with
 respect to PAR(k), k = 1, ..., NPAR, are

                      f(x(i),PARpk) - f(x(i),PARmk)
   Dcd(i,k) = --------------------------------------------   for i = 1, ..., N,
              (3*10**(-NETA))**(1/3)*SCALE(k)*SIGN(PAR(k))

 where

 PARpk   is  a  vector  which  has  the  same values as PAR except that the kth
         parameter is equal to

         PAR(k) + (3*10**(-NETA))**(1/3)*SCALE(k)*SIGN(PAR(k)).

 PARmk   is a vector which  has the  same values as  PAR  except  that  the kth
         parameter is equal to

         PAR(k) - (3*10**(-NETA))**(1/3)*SCALE(k)*SIGN(PAR(k)).

      The relative  step size is  considered acceptable  if, for  at  least N-a
 observations,

      |Dfd(i,k) - Dcd(i,k)| <= min(10**(-NETA/4), .02)  for i = 1, ..., N,

                                     <9-22>
1
 where a is  the  number  of  observations  exempted  from  meeting  the  above
 acceptance criterion [see section D, argument EXMPT].  If the step size is not
 acceptable,  it is adjusted by factors of 10 until the  condition  is  met  or
 until  no further decrease in the number of failures can be made,  although in
 no case will the selected relative step size be greater than 1.0 or less  than
 10**(-NETA).

      Note that the step  size selection  subroutines will return  the selected
 step  sizes even when the number  of failures exceeds the allowed  value; this
 condition will be noted  by the  value of IERR.   The detailed  printed output
 should always be examined for  problems discovered by the step  size selection
 subroutines.


 E.1.c  Derivative Checking

      The STARPAC derivative checking subroutines use an algorithm developed by
 Schnabel [1982]  to  determine the  validity of  the  user-supplied derivative
 subroutine.  The user-supplied derivative subroutine is considered correct for
 a given row i, i = 1, ..., N, and parameter PAR(k), k = 1, ..., NPAR, if

                    |Dfd(i,k) - D(i,k)| <= 10**(-t)*|D(i,k)|
 where

 D       is the derivative computed by the user's subroutine.

 Dfd     is the  forward difference  quotient approximation  to  the derivative
         described in section E.1.b.

 t       is the  agreement  tolerance,  i.e.,  number of  digits  of  agreement
         required between  D and Dfd, which must  be less than or equal  to the
         number of good digits in Dfd [see section D, argument NTAU].

      When the agreement  tolerance is  not satisfied  the  checking subroutine
 attempts  to determine  whether the disagreement  is due  to an  error  in the
 user's code or is due to the inaccuracy of the difference  quotient 
 approximation,  caused either by high curvature  in the user's model or  by 
 significant roundoff error.

      The  derivative  checking  subroutines  each  check  only  one row of the
 derivative matrix.  The user should examine the row at which  the  derivatives
 were  checked  to  ensure  that  some  relation  between  the  parameters  and
 independent variables,  such as a zero-valued independent variable or a factor
 (x(i) - PAR(k)) when x(i) = PAR(k), is not hiding the effect of an incorrectly
 computed derivative.  Checking only one row is appropriate since the same code
 is frequently used to compute the model function and derivatives at each row i
 = 1, ...,  N,  as is the case in the examples shown in section F.  If the code
 used  to  express  the model function and derivatives is not the same for each
 row,  then each distinct section of the  code  should  be  checked  by  making
 multiple calls to DCKLSC with argument NROW set to a row within each section.








                                     <9-23>
1E.2  Computed Results and Printed Output

 E.2.a  The Nonlinear Least Squares Estimation Subroutines

      The argument controlling  the  printed  output,  NPRT,  is  discussed  in
 section D.

      The  output  from  the  nonlinear  least  squares  estimation subroutines
 consists of five sections,  several of which include  tables  summarizing  the
 results.  In  the following descriptions,  the actual table headings are given
 by the uppercase phrases enclosed  in  angle  braces  (<...>).  Results  which
 correspond  to  input  or  returned  subroutine  CALL  statement arguments are
 identified by the argument name in uppercase (not enclosed in angle braces).


 Section 1 provides a summary of the initial estimates and control values.   It
           lists the following information.

       * The initial values of the parameters, PAR, and whether they are  to be
         held fixed or not as specified by IFIXED.

       * The scale values, SCALE.

       * Either the step sizes used to approximate the derivatives numerically,
         or, when user-supplied (analytic) derivatives are used, the results of
         the checking procedure; and the control values used in  these 
	 computations as applicable [see section E.1.b and section E.1.c].

       * The number of observations, N.

       * The number of observations with nonzero weights, NNZW.

       * The number of independent variables, M.

       * The maximum number of iterations allowed, MIT.

       * The maximum number of model subroutine calls allowed.

       * The two convergence criteria, STOPSS and STOPP.

       * The maximum change in  the parameters allowed at the  first iteration,
         DELTA.

       * The residual  sum  of squares  computed using  the  starting parameter
         values.

       * The  residual  standard deviation,  RSD, computed  using  the starting
         parameter values.


 Section 2 lists selected  information about  each iteration  and  includes the
           reason the iterations were terminated.  The information provided for
           each iteration includes the following.

       * The iteration number.




                                     <9-24>
1      * <MODEL  CALLS>:   the total number of times since execution began that
         the user's model subroutine  has  been  called,  not  including  calls
         required to approximate the derivatives numerically.

       * <RSD>:   the  residual standard deviation computed using the parameter
         values from the current iteration.

       * <RSS>: the residual sum of squares computed using the parameter values
         from the current iteration.

       * <REL CHNG RSS>:  the relative change in the residual  sum  of  squares
         caused by the current iteration.

       * <FORECASTED  REL  CHNG  RSS>:   the  forecasted relative change in the
         residual sum of squares at the current  iteration,  and  whether  this
         value was checked against STOPSS (<CHKD> = Y) or not (<CHKD> = N).

       * <REL CHNG PAR>:   the maximum scaled relative change in the parameters
         at the current iteration,  and whether this value was checked  against
         STOPP (<CHKD> = Y) or not (<CHKD> = N).

       * <CURRENT PARAMETER VALUES>:   the estimated parameter values resulting
         from the current iteration.


 Section 3 provides the following information for each observation, i = 1, ...,
           N, based on the final solution.

       * <ROW>:  the row number of the observations.

       * <PREDICTOR VALUES>:  the values for up to the first three  columns  of
         the independent variable matrix, XM, not including the first column if
         it is constant.

       * <DEPENDENT VARIABLE>:  the values of the dependent variable, Y.

       * <PREDICTED VALUE>:  the estimated predicted values, PV, from the fit.

       * <STD  DEV  OF  PRED VALUE>:   the standard deviations of the predicted
         values, SDPV.

       * <RESIDUAL>:  the error estimates, RES.

       * <STD RES>:  the standardized residuals, SDRES.

       * <WEIGHT>:  the user-supplied weights,  WT,  printed only when weighted
         analysis is performed.


 Section 4 displays the following plots of the standardized residuals.

       * The standardized residuals versus row numbers.

       * The standardized residuals versus predicted values.

       * The autocorrelation function of the (non-standardized) residuals.

       * The normal probability plot of the standardized residuals.

                                     <9-25>
1

 Section 5 summarizes  the  following  information about  the  final  parameter
          estimates and their variances.

       * The  variance-covariance  matrix,  VCV,  of  the  estimated  (unfixed)
         parameters, and the corresponding correlation matrix,

         rjk = VCV(j,k) / sqrt(VCV(j,j)*VCV(k,k))  for j = 1, ..., NPARE
                                                   and k = 1, ..., NPARE.

       * <PARAMETER>:  the final value of each parameter,  PAR(k),  k = 1, ...,
         NPAR.

       * <SD OF PAR>:  the standard deviation of each estimated parameter,

         sqrt(VCV(k,k))   for k = 1, ..., NPAR.

       * <RATIO>:  the ratio  of  each  estimated  parameter  to  its  standard
         deviation,

         RATIO(k) = PAR(k) / sqrt(VCV(k,k))  for k = 1, ..., NPAR.

       * <APPROXIMATE 95-PERCENT CONFIDENCE LIMITS>:    the  lower   and  upper
         95-percent confidence  limits for  each parameter, computed  using the
         appropriate  value of  the Student's t  distribution  with  NNZW-NPARE
         degrees of freedom.

       * the residual sum of squares, RSS(PAR).

       * the residual standard deviation at the solution, RSD.

       * the residual degrees of freedom, NNZW-NPARE.

       * an approximation to the  condition number of the derivative  matrix, D
         (the Jacobian), under the  assumption that the absolute error  in each
         column  of D is roughly equal.   The approximation will be meaningless
         if  this  assumption   is  not   valid;  otherwise  it   will  usually
         underestimate the actual condition number by a factor of from 2  to 10
         [see Dongarra et al., 1979, p.  9.5].  (Note that the condition number
         returned by the nonlinear least squares subroutines is not exactly the
         same as  that  returned  by the  linear least  subroutines  because of
         differences in the  computational procedures  used by the two families
         of subroutines.)

 NOTE:   The  standard  deviation of  the predicted  values,  the  standardized
 residuals,  the  variance-covariance matrix,  the standard  deviations  of the
 parameters and  the 95-percent  confidence limits  on the  parameters  are all
 based  on a  linear approximation  to  the  model  in a  neighborhood  of  the
 solution; the validity of  this approximation  depends on the  nonlinearity of
 the  model.   The  statistics  based on  this approximation  may  be extremely
 inaccurate for a problem with a highly nonlinear model.


 E.2.b  The Derivative Step Size Selection Subroutines

      The  argument  controlling  the  printed  output,  NPRT,  is discussed in
 section D.

                                     <9-26>
1
      The output from the step size selection subroutines consists of a summary
 of the input and control values and, for each parameter, the selected relative
























































                                     <9-27>
1step size, the number of observations at which this step size failed  the step
 size selection criteria and the row numbers at which the failures occurred.


 E.2.c  The Derivative Checking Subroutines

      The argument controlling  the  printed  output,  NPRT,  is  discussed  in
 section D.

      The output for the derivative checking subroutines consists of  a summary
 of  the input and control values  and the  results of the  derivative checking
 test with respect to each of the model parameters,  PAR(k),  k = 1, ..., NPAR.
 The possible test results are:

 OK -

       * The user-supplied derivative and the numerical derivative agree to the
         required number of digits.

 QUESTIONABLE -

       * The user-supplied derivative and the approximated derivative  agree to
         the  required number of digits but  both are equal to zero.   The user
         should recheck the derivative at another row.

       * The user-supplied  derivative and  the approximated derivative  do not
         agree to the required  number of digits but the  user-supplied 
	 derivative is identically  zero and  the approximated  derivative
	 is nearly zero.  The user should recheck the derivative at another row.

       * The user-supplied derivative and the approximated  derivative disagree
         but the user-supplied derivative is identically zero.  The user should
         recheck the derivative at another row.

       * The user-supplied derivative and the approximated  derivative disagree
         but  the  validity  of  the approximated  derivative  is  questionable
         because either the ratio of the relative curvature of the model to the
         slope of the model is too high or SCALE(k) is wrong.

       * The user-supplied derivative and the approximated  derivative disagree
         but the validity  of the estimated derivative is  questionable because
         the ratio of the relative curvature of  the model to the slope  of the
         model is too high.

 INCORRECT -

       * The user-supplied derivative and the approximated derivative disagree,
         and there is no reason  to question  the accuracy of  the approximated
         derivative.


 F.  Examples

      The  sample programs  of this section  use the  model and  data  given in
 example one, pages 428 to 441 of Daniel and Wood [1980]; the model is

             f(x(i),b) = PAR(1)*x(i,1)**PAR(2)  for i = 1, ..., N.


                                     <9-28>
1
      Nonlinear Least Squares Estimation.  In the first example program  below,
 NLS   is  used  to  compute  the  least  squares  solution  using  numerically
 approximated derivatives.  In the second example  program,  NLSD  is  used  to
 compute the least squares solution given analytic derivatives.


      Derivative  Step Size Selection.  In the third example program,  STPLS is
 used to compute the optimum  step  sizes  for  numerically  approximating  the
 derivatives with respect to each of the parameters, PAR(k), k = 1, 2.


      Derivative Checking.  In the fourth example program below,  DCKLS is used
 to  check  the  validity  of  a  user-supplied  derivative   subroutine.   The
 derivative  subroutine  has  been  intentionally coded incorrectly in order to
 display the report obtained when the derivative checking subroutine determines
 the derivatives are incorrect,  and the starting parameter  values  have  been
 chosen  in  order  to  display  the  report obtained when the test results are
 questionable.








































                                     <9-29>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE NLS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, XM, PAR AND RES MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(10), XM(10,5), PAR(5), RES(10)
       DOUBLE PRECISION DSTAK(200)
 C
       EXTERNAL NLSMDL
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
       IXM = 10
 C
 C     READ NUMBER OF PARAMETERS
 C          STARTING VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS AND NUMBER OF INDEPENDENT VARIABLES
 C          INDEPENDENT AND DEPENDENT VARIABLES
 C
       READ (5,100) NPAR
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,100) N, M
       READ (5,101) ((XM(I,J), I=1,N), J=1,M), (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL NLS TO PERFORM NONLINEAR REGRESSION
 C     WITH NUMERICALLY APPROXIMATED DERIVATIVES
 C
       WRITE (IPRT,102)
       CALL NLS (Y, XM, N, M, IXM, NLSMDL, PAR, NPAR, RES, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (6F6.3)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' NONLINEAR LEAST SQUARES SUBROUTINE NLS')
       END
       SUBROUTINE NLSMDL (PAR, NPAR, XM, N, M, IXM, PV)
 C
 C     SUBROUTINE TO COMPUTE PREDICTED VALUES OF DEPENDENT VARIABLE
 C
 C     N.B. DECLARATION OF PAR, XM AND PV MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C

                                     <9-30>
1      REAL PAR(NPAR), XM(IXM,M), PV(N)
 C
       DO 10 I = 1, N
         PV(I) = PAR(1) * XM(I, 1) ** PAR(2)
    10 CONTINUE
 C
       RETURN
       END


 Data:

     2
  0.725 4.000
     6    1
  1.309 1.471 1.490 1.565 1.611 1.680
  2.138 3.421 3.597 4.340 4.882 5.660










































                                     <9-31>
1RESULTS OF STARPAC NONLINEAR LEAST SQUARES SUBROUTINE NLS
                                                                                                         STARPAC 2.08S (03/15/90)
+**********************************************************************************
 *  NONLINEAR LEAST SQUARES ESTIMATION WITH NUMERICALLY APPROXIMATED DERIVATIVES  *
 **********************************************************************************


 SUMMARY OF INITIAL CONDITIONS
 ------------------------------


                                                  STEP SIZE FOR    OBSERVATIONS FAILING STEP SIZE SELECTION CRITERIA
                                                  APPROXIMATING                     *
       PARAMETER STARTING VALUE      SCALE          DERIVATIVE       COUNT     NOTES     ROW NUMBER
 INDEX  FIXED      (PAR)            (SCALE)           (STP)                     F C

   1      NO    .72500000           DEFAULT       .46415888E-04         0
   2      NO    4.0000000           DEFAULT       .38782913E-06         0


 *  NOTES.  A PLUS (+) IN THE COLUMNS HEADED F OR C HAS THE FOLLOWING MEANING.

    F - NUMBER OF OBSERVATIONS FAILING STEP SIZE SELECTION CRITERIA EXCEEDS
        NUMBER OF EXEMPTIONS ALLOWED.

    C - HIGH CURVATURE IN THE MODEL IS SUSPECTED AS THE CAUSE OF
        ALL FAILURES NOTED.

 NUMBER OF RELIABLE DIGITS IN MODEL RESULTS                         (NETA)    13

 PROPORTION OF OBSERVATIONS EXEMPTED FROM SELECTION CRITERIA       (EXMPT)   .1000

 NUMBER OF OBSERVATIONS EXEMPTED FROM SELECTION CRITERIA                       1

 NUMBER OF OBSERVATIONS                                                (N)     6

 NUMBER OF INDEPENDENT VARIABLES                                       (M)     1

 MAXIMUM NUMBER OF ITERATIONS ALLOWED                                (MIT)    21

 MAXIMUM NUMBER OF MODEL SUBROUTINE CALLS ALLOWED                             42

 CONVERGENCE CRITERION FOR TEST BASED ON THE

      FORECASTED RELATIVE CHANGE IN RESIDUAL SUM OF SQUARES       (STOPSS)   .3696E-09
      MAXIMUM SCALED RELATIVE CHANGE IN THE PARAMETERS             (STOPP)   .8425E-07


 MAXIMUM CHANGE ALLOWED IN THE PARAMETERS AT THE FIRST ITERATION   (DELTA)   100.0

 RESIDUAL SUM OF SQUARES FOR INPUT PARAMETER VALUES                          .1472E-01

 RESIDUAL STANDARD DEVIATION FOR INPUT PARAMETER VALUES              (RSD)   .6067E-01







                                                              <9-32>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH NUMERICALLY APPROXIMATED DERIVATIVES, CONTINUED


 ITERATION NUMBER    1
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         2       .3390E-01       .4597E-02       .6877       .7109       Y   .1790E-01   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2
          VALUE   .7679852       3.859309


 ITERATION NUMBER    4
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         5       .3285E-01       .4317E-02       .8936E-12   .7028E-12   Y   .1123E-07   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2
          VALUE   .7688623       3.860406

 ***** PARAMETER AND RESIDUAL SUM OF SQUARES CONVERGENCE *****































                                                              <9-33>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH NUMERICALLY APPROXIMATED DERIVATIVES, CONTINUED


 RESULTS FROM LEAST SQUARES FIT
 -------------------------------

                                                     DEPENDENT       PREDICTED      STD DEV OF                         STD
  ROW             PREDICTOR VALUES                    VARIABLE         VALUE        PRED VALUE        RESIDUAL         RES

    1                 1.3090000                      2.1380000       2.1741175       .22079043E-01  -.36117492E-01   -1.48
    2                 1.4710000                      3.4210000       3.4111549       .16469586E-01   .98450831E-02     .35
    3                 1.4900000                      3.5970000       3.5844108       .15615321E-01   .12589151E-01     .44
    4                 1.5650000                      4.3400000       4.3326419       .14065814E-01   .73580833E-02     .25
    5                 1.6110000                      4.8820000       4.8453073       .16512112E-01   .36692701E-01    1.29
    6                 1.6800000                      5.6600000       5.6968365       .26183728E-01  -.36836492E-01   -1.86











































                                                              <9-34>
1                                                                                                        STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH NUMERICALLY APPROXIMATED DERIVATIVES, CONTINUED

                     STD RES VS ROW NUMBER                                     STD RES VS PREDICTED VALUES
  3.75++---------+---------+----+----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  2.25+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                        *          -        -                                      *            -
   .75+                                                   +     .75+                                                   +
      -                                                   -        -                                                   -
      -          *         *         *                    -        -                  * *          *                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  -.75+                                                   +    -.75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -*                                                  -        -*                                                  -
      -                                                  *-        -                                                  *-
 -2.25+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
 -3.75++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
      1.0                      3.5                      6.0      2.174                    3.935                5.697

             AUTOCORRELATION FUNCTION OF RESIDUALS                        NORMAL PROBABILITY PLOT OF STD RES
     1++---------+-------********----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                        **                         -        -                                                   -
      -                       ***                         -        -                                                   -
      -                **********                         -        -                                                   -
      -                         ********                  -        -                                                   -
     6+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                     *             -
    11+                                                   +     .75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                       *   *   *                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    16+                                                   +    -.75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                   *                               -
      -                                                   -        -             *                                     -
    21+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    26++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
    -1.00                      0.0                     1.00     -2.5                       0.0                      2.5
                                                              <9-35>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH NUMERICALLY APPROXIMATED DERIVATIVES, CONTINUED



 VARIANCE-COVARIANCE AND CORRELATION MATRICES OF THE ESTIMATED (UNFIXED) PARAMETERS
 ----------------------------------------------------------------------------------

    - APPROXIMATION BASED ON ASSUMPTION THAT RESIDUALS ARE SMALL
    - COVARIANCES ARE ABOVE THE DIAGONAL
    - VARIANCES ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

 COLUMN         1                2

      1      .3342304E-03    -.9369370E-03
      2     -.9907719         .2675639E-02



 ESTIMATES FROM LEAST SQUARES FIT
 ---------------------------------


                                                                     APPROXIMATE
                                                             95 PERCENT CONFIDENCE LIMITS
 INDEX  FIXED   PARAMETER        SD OF PAR       RATIO            LOWER            UPPER

   1      NO    .76886226        .18281968E-01   42.06        .71810338        .81962114
   2      NO    3.8604056        .51726577E-01   74.63        3.7167896        4.0040216


 RESIDUAL SUM OF SQUARES                  .4317308E-02

 RESIDUAL STANDARD DEVIATION              .3285311E-01
 BASED ON DEGREES OF FREEDOM        6 -   2 =    4

 APPROXIMATE CONDITION NUMBER             20.87491





















                                                              <9-36>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE NLSD USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, XM, PAR AND RES MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(10), XM(10,5), PAR(5), RES(10)
       DOUBLE PRECISION DSTAK(200)
 C
       EXTERNAL NLSMDL, NLSDRV
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
       IXM = 10
 C
 C     READ NUMBER OF PARAMETERS
 C          STARTING VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS AND NUMBER OF INDEPENDENT VARIABLES
 C          INDEPENDENT AND DEPENDENT VARIABLES
 C
       READ (5,100) NPAR
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,100) N, M
       READ (5,101) ((XM(I,J), I=1,N), J=1,M), (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL NLSD TO PERFORM NONLINEAR REGRESSION
 C     WITH USER-SUPPLIED DERIVATIVES
 C
       WRITE (IPRT,102)
       CALL NLSD (Y, XM, N, M, IXM, NLSMDL, NLSDRV, PAR, NPAR, RES,
      *  LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (6F6.3)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' NONLINEAR LEAST SQUARES SUBROUTINE NLSD')
       END
       SUBROUTINE NLSMDL (PAR, NPAR, XM, N, M, IXM, PV)
 C
 C     SUBROUTINE TO COMPUTE PREDICTED VALUES OF DEPENDENT VARIABLE
 C
 C     N.B. DECLARATION OF PAR, XM AND PV MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.

                                     <9-37>
1C
       REAL PAR(NPAR), XM(IXM,M), PV(N)
 C
       DO 10 I = 1, N
         PV(I) = PAR(1) * XM(I, 1) ** PAR(2)
    10 CONTINUE
 C
       RETURN
       END
       SUBROUTINE NLSDRV (PAR, NPAR, XM, N, M, IXM, D)
 C
 C     SUBROUTINE TO COMPUTE THE PARTIAL DERIVATIVE (JACOBIAN) MATRIX
 C
 C     N.B. DECLARATION OF PAR, XM AND D MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL PAR(NPAR), XM(IXM,M), D(N,NPAR)
 C
       DO 10 I = 1, N
         D(I,1) = XM(I,1) ** PAR(2)
         D(I,2) = PAR(1) * XM(I,1) ** PAR(2) * ALOG(XM(I,1))
    10 CONTINUE
 C
       RETURN
       END


 Data:

     2
  0.725 4.000
     6    1
  1.309 1.471 1.490 1.565 1.611 1.680
  2.138 3.421 3.597 4.340 4.882 5.660

























                                     <9-38>
1RESULTS OF STARPAC NONLINEAR LEAST SQUARES SUBROUTINE NLSD
                                                                                                         STARPAC 2.08S (03/15/90)
+***********************************************************************
 *  NONLINEAR LEAST SQUARES ESTIMATION WITH USER-SUPPLIED DERIVATIVES  *
 ***********************************************************************


 SUMMARY OF INITIAL CONDITIONS
 ------------------------------



                                                    DERIVATIVE
       PARAMETER STARTING VALUE      SCALE          ASSESSMENT
 INDEX  FIXED      (PAR)            (SCALE)

   1      NO    .72500000           DEFAULT             OK
   2      NO    4.0000000           DEFAULT             OK

 NUMBER OF RELIABLE DIGITS IN MODEL RESULTS                         (NETA)    13

 NUMBER OF DIGITS IN DERIVATIVE CHECKING AGREEMENT TOLERANCE        (NTAU)     4

 ROW NUMBER AT WHICH DERIVATIVES WERE CHECKED                       (NROW)     1
   -VALUES OF THE INDEPENDENT VARIABLES AT THIS ROW
          INDEX    1
          VALUE   1.309000

 NUMBER OF OBSERVATIONS                                                (N)     6

 NUMBER OF INDEPENDENT VARIABLES                                       (M)     1

 MAXIMUM NUMBER OF ITERATIONS ALLOWED                                (MIT)    21

 MAXIMUM NUMBER OF MODEL SUBROUTINE CALLS ALLOWED                             42

 CONVERGENCE CRITERION FOR TEST BASED ON THE

      FORECASTED RELATIVE CHANGE IN RESIDUAL SUM OF SQUARES       (STOPSS)   .3696E-09
      MAXIMUM SCALED RELATIVE CHANGE IN THE PARAMETERS             (STOPP)   .8425E-07


 MAXIMUM CHANGE ALLOWED IN THE PARAMETERS AT THE FIRST ITERATION   (DELTA)   100.0

 RESIDUAL SUM OF SQUARES FOR INPUT PARAMETER VALUES                          .1472E-01

 RESIDUAL STANDARD DEVIATION FOR INPUT PARAMETER VALUES              (RSD)   .6067E-01













                                                              <9-39>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH USER-SUPPLIED DERIVATIVES, CONTINUED


 ITERATION NUMBER    1
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         2       .3390E-01       .4597E-02       .6877       .7109       Y   .1790E-01   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2
          VALUE   .7679852       3.859309


 ITERATION NUMBER    4
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         5       .3285E-01       .4317E-02      -.3214E-13   .6352E-12   Y   .1068E-07   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2
          VALUE   .7688623       3.860406

 ***** PARAMETER AND RESIDUAL SUM OF SQUARES CONVERGENCE *****































                                                              <9-40>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH USER-SUPPLIED DERIVATIVES, CONTINUED


 RESULTS FROM LEAST SQUARES FIT
 -------------------------------

                                                     DEPENDENT       PREDICTED      STD DEV OF                         STD
  ROW             PREDICTOR VALUES                    VARIABLE         VALUE        PRED VALUE        RESIDUAL         RES

    1                 1.3090000                      2.1380000       2.1741175       .22079044E-01  -.36117523E-01   -1.48
    2                 1.4710000                      3.4210000       3.4111549       .16469585E-01   .98450648E-02     .35
    3                 1.4900000                      3.5970000       3.5844109       .15615321E-01   .12589135E-01     .44
    4                 1.5650000                      4.3400000       4.3326419       .14065814E-01   .73580808E-02     .25
    5                 1.6110000                      4.8820000       4.8453073       .16512112E-01   .36692709E-01    1.29
    6                 1.6800000                      5.6600000       5.6968365       .26183727E-01  -.36836464E-01   -1.86











































                                                              <9-41>
1                                                                                                        STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH USER-SUPPLIED DERIVATIVES, CONTINUED

                     STD RES VS ROW NUMBER                                     STD RES VS PREDICTED VALUES
  3.75++---------+---------+----+----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  2.25+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                        *          -        -                                      *            -
   .75+                                                   +     .75+                                                   +
      -                                                   -        -                                                   -
      -          *         *         *                    -        -                  * *          *                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
  -.75+                                                   +    -.75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -*                                                  -        -*                                                  -
      -                                                  *-        -                                                  *-
 -2.25+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
 -3.75++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
      1.0                      3.5                      6.0      2.174                    3.935                5.697

             AUTOCORRELATION FUNCTION OF RESIDUALS                        NORMAL PROBABILITY PLOT OF STD RES
     1++---------+-------********----+---------+---------++    3.75++---------+---------+----+----+---------+---------++
      -                        **                         -        -                                                   -
      -                       ***                         -        -                                                   -
      -                **********                         -        -                                                   -
      -                         ********                  -        -                                                   -
     6+                                                   +    2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                     *             -
    11+                                                   +     .75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                       *   *   *                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    16+                                                   +    -.75+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                   *                               -
      -                                                   -        -             *                                     -
    21+                                                   +   -2.25+                                                   +
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
      -                                                   -        -                                                   -
    26++---------+---------+----+----+---------+---------++   -3.75++---------+---------+----+----+---------+---------++
    -1.00                      0.0                     1.00     -2.5                       0.0                      2.5
                                                              <9-42>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION WITH USER-SUPPLIED DERIVATIVES, CONTINUED



 VARIANCE-COVARIANCE AND CORRELATION MATRICES OF THE ESTIMATED (UNFIXED) PARAMETERS
 ----------------------------------------------------------------------------------

    - APPROXIMATION BASED ON ASSUMPTION THAT RESIDUALS ARE SMALL
    - COVARIANCES ARE ABOVE THE DIAGONAL
    - VARIANCES ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

 COLUMN         1                2

      1      .3342306E-03    -.9369379E-03
      2     -.9907719         .2675642E-02



 ESTIMATES FROM LEAST SQUARES FIT
 ---------------------------------


                                                                     APPROXIMATE
                                                             95 PERCENT CONFIDENCE LIMITS
 INDEX  FIXED   PARAMETER        SD OF PAR       RATIO            LOWER            UPPER

   1      NO    .76886229        .18281974E-01   42.06        .71810339        .81962119
   2      NO    3.8604055        .51726611E-01   74.63        3.7167894        4.0040216


 RESIDUAL SUM OF SQUARES                  .4317308E-02

 RESIDUAL STANDARD DEVIATION              .3285311E-01
 BASED ON DEGREES OF FREEDOM        6 -   2 =    4

 APPROXIMATE CONDITION NUMBER             20.87492





















                                                              <9-43>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE STPLS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF XM, PAR AND STP MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL XM(10,5), PAR(5), STP(5)
       DOUBLE PRECISION DSTAK(200)
 C
       EXTERNAL NLSMDL, DERIV
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
       IXM = 10
 C
 C     READ NUMBER OF PARAMETERS
 C          STARTING VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS AND NUMBER OF INDEPENDENT VARIABLES
 C          INDEPENDENT VARIABLES
 C
       READ (5,100) NPAR
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,100) N, M
       READ (5,101) ((XM(I,J), I=1,N), J=1,M)
 C
 C     PRINT TITLE AND CALL STPLS TO SELECT STEP SIZES FOR
 C     APPROXIMATING DERIVATIVES
 C
       WRITE (IPRT,102)
       CALL STPLS (XM, N, M, IXM, NLSMDL, PAR, NPAR, LDSTAK, STP)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (6F6.3)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' DERIVATIVE STEP SIZE SELECTION SUBROUTINE STPLS')
       END
       SUBROUTINE NLSMDL (PAR, NPAR, XM, N, M, IXM, PV)
 C
 C     SUBROUTINE TO COMPUTE PREDICTED VALUES OF DEPENDENT VARIABLE
 C
 C     N.B. DECLARATION OF PAR, XM AND PV MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C

                                     <9-44>
1      REAL PAR(NPAR), XM(IXM,M), PV(N)
 C
       DO 10 I = 1, N
         PV(I) = PAR(1) * XM(I, 1) ** PAR(2)
    10 CONTINUE
 C
       RETURN
       END


 Data:

     2
  0.725 4.000
     6    1
  1.309 1.471 1.490 1.565 1.611 1.680











































                                     <9-45>
1RESULTS OF STARPAC DERIVATIVE STEP SIZE SELECTION SUBROUTINE STPLS
                                                                                                         STARPAC 2.08S (03/15/90)
+**********************************
 * DERIVATIVE STEP SIZE SELECTION *
 **********************************


                                                  STEP SIZE FOR    OBSERVATIONS FAILING STEP SIZE SELECTION CRITERIA
                  PARAMETER                       APPROXIMATING                     *
                STARTING VALUE       SCALE          DERIVATIVE       COUNT     NOTES     ROW NUMBER(S)
 INDEX             (PAR)            (SCALE)           (STP)                     F C

   1            .72500000           DEFAULT       .46415888E-04         0
   2            4.0000000           DEFAULT       .38782913E-06         0


 *  NOTES.  A PLUS (+) IN THE COLUMNS HEADED F OR C HAS THE FOLLOWING MEANING.

    F - NUMBER OF OBSERVATIONS FAILING STEP SIZE SELECTION CRITERIA EXCEEDS
        NUMBER OF EXEMPTIONS ALLOWED.

    C - HIGH CURVATURE IN THE MODEL IS SUSPECTED AS THE CAUSE OF
        ALL FAILURES NOTED.

 NUMBER OF RELIABLE DIGITS IN MODEL RESULTS                         (NETA)    13

 PROPORTION OF OBSERVATIONS EXEMPTED FROM SELECTION CRITERIA       (EXMPT)   .1000

 NUMBER OF OBSERVATIONS EXEMPTED FROM SELECTION CRITERIA                       1

 NUMBER OF OBSERVATIONS                                                (N)     6





























                                                              <9-46>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE DCKLS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF XM AND PAR MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL XM(10,5), PAR(5)
       DOUBLE PRECISION DSTAK(200)
 C
       EXTERNAL NLSMDL, NLSDRV
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 200
       IXM = 10
 C
 C     READ NUMBER OF PARAMETERS
 C          STARTING VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS AND NUMBER OF INDEPENDENT VARIABLES
 C          INDEPENDENT VARIABLES
 C
       READ (5,100) NPAR
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,100) N, M
       READ (5,101) ((XM(I,J), I=1,N), J=1,M)
 C
 C     PRINT TITLE AND CALL DCKLS TO PERFORM DERIVATIVE CHECKING
 C
       WRITE (IPRT,102)
       CALL DCKLS (XM, N, M, IXM, NLSMDL, NLSDRV, PAR, NPAR, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (2I5)
   101 FORMAT (6F6.3)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' DERIVATIVE CHECKING SUBROUTINE DCKLS')
       END
       SUBROUTINE NLSMDL (PAR, NPAR, XM, N, M, IXM, PV)
 C
 C     SUBROUTINE TO COMPUTE PREDICTED VALUES OF DEPENDENT VARIABLE
 C
 C     N.B. DECLARATION OF PAR, XM AND PV MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL PAR(NPAR), XM(IXM,M), PV(N)

                                     <9-47>
1C
       DO 10 I = 1, N
         PV(I) = PAR(1) * XM(I, 1) ** PAR(2)
    10 CONTINUE
 C
       RETURN
       END
       SUBROUTINE NLSDRV (PAR, NPAR, XM, N, M, IXM, D)
 C
 C     SUBROUTINE TO COMPUTE THE PARTIAL DERIVATIVE (JACOBIAN) MATRIX
 C
 C     DERIVATIVE WITH RESPECT TO FIRST PARAMETER HAS BEEN CODED
 C     INCORRECTLY TO DEMONSTRATE ERROR DETECTION CAPABILITIES
 C
 C     N.B. DECLARATION OF PAR, XM AND D MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL PAR(NPAR), XM(IXM,M), D(N,NPAR)
 C
       DO 10 I = 1, N
         D(I,1) = XM(I,1) * PAR(2)
         D(I,2) = PAR(1) * XM(I,1) ** PAR(1) * ALOG(XM(I,1))
    10 CONTINUE
 C
       RETURN
       END


 Data:

     2
  0.000 4.000
     6    1
  1.309 1.471 1.490 1.565 1.611 1.680

























                                     <9-48>
1RESULTS OF STARPAC DERIVATIVE CHECKING SUBROUTINE DCKLS
                                                                                                         STARPAC 2.08S (03/15/90)
+***********************
 * DERIVATIVE CHECKING *
 ***********************



                                                              *
                 PARAMETER                          DERIVATIVE
               STARTING VALUE        SCALE          ASSESSMENT
 INDEX             (PAR)            (SCALE)

   1           0.                   DEFAULT          INCORRECT
   2            4.0000000           DEFAULT        QUESTIONABLE (1)

 * NUMBERS IN PARENTHESES REFER TO THE FOLLOWING NOTES.

  (1) USER-SUPPLIED AND APPROXIMATED DERIVATIVES AGREE, BUT
      BOTH ARE ZERO.  RECHECK AT ANOTHER ROW.

 NUMBER OF RELIABLE DIGITS IN MODEL RESULTS                         (NETA)    14

 NUMBER OF DIGITS IN DERIVATIVE CHECKING AGREEMENT TOLERANCE        (NTAU)     4

 ROW NUMBER AT WHICH DERIVATIVES WERE CHECKED                       (NROW)     1
   -VALUES OF THE INDEPENDENT VARIABLES AT THIS ROW
          INDEX    1
          VALUE   1.309000

 NUMBER OF OBSERVATIONS                                                (N)     6





























                                                              <9-49>
1G.  Acknowledgments

      The subroutines used to compute the nonlinear least squares  solution are
 those  referenced  in  Dennis  et al.  [1981].   The algorithms used to select
 optimum  step  sizes  for  numerical  derivatives,   and  to  check   analytic
 derivatives  were  developed  by  Schnabel [1982].  The printed output for the
 nonlinear least squares subroutines has  been  modeled  on  the  linear  least
 squares output used by OMNITAB II [Hogben et al., 1971].



















































                                     <9-50>
1-----                            CHAPTER 10                              -----

                              DIGITAL FILTERING


 A.  Introduction

      STARPAC contains 16 subroutines for digital filtering time series.  These
 include subroutines which:   compute a least squares approximation to an ideal
 low-pass filter;  perform symmetric  linear filter  operations;  sample values
 from a series;  perform autoregressive (or difference) filter  operations; and
 compute  the gain  function of any  symmetric linear  filter and the  gain and
 phase functions of any autoregressive (or difference) filter.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 the  output  produced by these subroutines are discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

 B.1  Symmetric Linear Filter Subroutines

      Subroutine LPCOEF computes symmetric linear low-pass  filter coefficients
 using a  least squares  approximation  to an  ideal low-pass  filter  that has
 convergence  factors  which reduce  overshoot and  ripple  [Bloomfield, 1976].
 This low-pass filter has a transfer function which changes  from approximately
 one to zero in a transition band about the ideal cutoff frequency FC,  that is
 from (FC - 1/K) to (FC + 1/K),  as discussed  in  section  6.4  of  Bloomfield
 [1976].  The  user  must  specify  the  cutoff  frequency in cycles per sample
 interval and the number of filter coefficients,  which must be odd.  The  user
 must also choose the number of filter terms, K, so that (FC - 1/K) > 0 and (FC
 + 1/K) < 0.5.  In addition, K must be chosen as a compromise between:

      1)  A sharp cutoff, that is, 1/K small; and

      2)  Minimizing the number of data points lost by the filtering operations
          ((K-1)/2 data points will be lost from each end of the series).

 The subroutine returns the normalized low-pass filter coefficients.   There is
 no printed output.

      For any  low-pass  filter  there  is  a  corresponding  high-pass  filter
 equivalent to  subtracting  the low-pass  filtered series  from  the  original
 series.    Subroutine  HPCOEF   returns  symmetric   linear  high-pass  filter
 coefficients  computed  from user  supplied low-pass  symmetric  linear filter











                                    <10-1>
1coefficients.   The number  of filter coefficients must be  odd.   There is no
 printed output.

      Subroutine MAFLT performs a simple moving average filter operation on the
 input series using the simple moving average filter defined by

                        HMA(J) = 1/K for J = 1, ..., K.

 The user must specify the number of filter coefficients, K, which must be odd;
 the subroutine returns the filtered  series and the number of  observations in
 the filtered series.  There is no printed output.

      Subroutine SLFLT performs  a symmetric linear filter operation  with user
 supplied coefficients and returns the filtered series to the user.  The filter
 coefficients must be  normalized on input to SLFLT.   The filtered  series and
 the number of observations in  the filtered series are returned.   There is no
 printed output.

      Subroutine SAMPLE samples  every NSth  observation from an  input series.
 If  the input  series was  obtained  using an  NS term  low-pass  filter, this
 sampling  rate  removes  the  autocorrelation  introduced  by   the  filtering
 operation.  This subroutine returns the series of sampled observations and the
 number of observations in the series.  There is no printed output.

      Subroutine LOPASS computes low-pass filter coefficients as  described for
 subroutine LPCOEF  and  then performs  the filtering  operation  described for
 subroutine  SLFLT.  The user must specify the  cutoff frequency  in cycles per
 sample interval  and the  number of  filter  terms, which  must be  odd.   The
 subroutine returns the normalized filter coefficients, the filtered series and
 the number  of  observations in  the filtered  series.   There  is  no printed
 output.

      Subroutine  HIPASS computes the high-pass filter  coefficients equivalent
 to using HPCOEF with the input low-pass filter coefficients supplied by LPCOEF
 and performs the filtering operation described for subroutine SLFLT.  The user
 must specify the cutoff frequency in cycles per sample interval and the number
 of  filter  terms, which  must be  odd.   The  subroutine  returns the  filter
 coefficients,  the filtered  series  and the  number of  observations  in  the
 filtered series.  There is no printed output.


 B.2  Autoregressive or Difference Linear Filter Subroutines

      Subroutine ARFLT  subtracts  the series  mean from  each  observation and
 performs an autoregressive  linear filter operation with user  supplied filter
 coefficients.   This subroutine returns the filtered series and  the number of
 observations in the filtered series.  There is no printed output.

      Subroutine DIF performs a first difference filter operation on  the input
 series.   It returns the differenced series and the  number of observations in
 the differenced series.  There is no printed output.

      Subroutine DIFC performs  a user  controlled differencing operation.   It
 returns the order of the difference filter specified, the coefficients  of the
 difference filter, the differenced  series and  the number of  observations in
 the differenced series.   This subroutine can be used as a high-pass filter or
 for differencing series in the style  of Box and Jenkins [1976].   There is no
 printed output.

                                    <10-2>
1
      Subroutines DIFM and DIFMC are  the same  as DIF and  DIFC, respectively,
 except that the input  series may contain missing data.   A missing value code
 must  be used  within the  input  series to  specify time  points  without  an
 observed value.  The difference between a missing value and an observed value,
 or  between  two  missing values,  will  result  in  a missing  value  in  the
 differenced series; the missing value  code used in the differenced  series is
 also returned  to the  user.   Users should  note that  the number of  missing
 values may be significantly increased by the differencing operation.


 B.3  Gain and Phase Function Subroutines

      Subroutine GFSLF computes the gain function of a symmetric linear filter.
 The printed output consists of a plot of the gain function versus frequency.

      Subroutine  GFSLFS is the same as  GFSLF but  allows the user  to specify
 various options which are preset in GFSLF, including the frequency  range, the
 number of frequencies for which the gain  function is to be computed  and  the
 type  of plot to be used.   In addition, the gain function is returned  to the
 user, permitting the use of other methods of displaying the results.

      Subroutine  GFARF  computes  the  gain  and  phase  functions  of  either
 autoregressive or difference filters.   The output  consists of a plot of  the
 gain and phase functions versus frequency.

      Subroutine  GFARFS is the same as  GFARF but  provides the user  with the
 same options as are available for subroutine GFSLFS.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.


           Subroutines Supporting Symmetric Linear Filter Operations

 LPCOEF:  Compute  symmetric linear low-pass filter coefficients; return filter
          coefficients (no printed output)

          <real> HLP(k)
          :
          :
          CALL LPCOEF (FC, K, HLP)

                                      ===

 HPCOEF:  Compute  symmetric  linear  high-pass   filter  coefficients;  return
          coefficients (no printed output)

          <real> HLP(k), HHP(k)
          :
          :
          CALL HPCOEF (HLP, K, HHP)

                                      ===

                                    <10-3>
1
 MAFLT:   Perform  simple  moving  average; return filtered series (no  printed
          output)

          <real> Y(n), YF(n)
          :
          :
          CALL MAFLT (Y, N, K, YF, NYF)

                                      ===

 SLFLT:   Perform  symmetric linear  filter operation with user-supplied filter
          coefficients; return filtered series (no printed output)

          <real> Y(n), H(k), YF(n)
          :
          :
          CALL SLFLT (Y, N, K, H, YF, NYF)

                                      ===

 SAMPLE:  Sample  (extract)  every  NSth  observation  from  a  series;  return
          sampled series (no printed output)

          <real> Y(n), YS(n)
          :
          :
          CALL SAMPLE (Y, N, NS, YS, NYS)

                                      ===

 LOPASS:  Filter series  with symmetric linear low-pass filter; return filtered
          series (no printed output)

          <real> Y(n), HLP(k), YF(n)
          :
          :
          CALL LOPASS (Y, N, FC, K, HLP, YF, NYF)

                                      ===

 HIPASS:  Filter   series  with  symmetric  linear  high-pass   filter;  return
          filtered series (no printed output)

          <real> Y(n), HHP(k), YF(n)
          :
          :
          CALL HIPASS (Y, N, FC, K, HHP, YF, NYF)

                                      ===









                                    <10-4>
1         Subroutines for Autoregressive or Difference Linear Filters

 ARFLT:   Perform  autoregressive  filter  operation with  user-supplied filter
          coefficients; return filtered series (no printed output)

          <real> Y(n), PHI(iar), YF(n)
          :
          :
          CALL ARFLT (Y, N, IAR, PHI, YF, NYF)

                                      ===

 DIF:     Perform first-difference  filter operation; return differenced series
          (no printed output)

          <real> Y(n), YF(n)
          :
          :
          CALL DIF (Y, N, YF, NYF)

                                      ===

 DIFC:    Perform  user-specified    difference    filter   operation;   return
          differenced series (no printed output)

          INTEGER ND(nfac), IOD(nfac)
          <real> Y(n), YF(n), PHI(iar)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL DIFC (Y, N,
         +           NFAC, ND, IOD, IAR, PHI, LPHI,
         +           YF, NYF, LDSTAK)

                                      ===

 DIFM:    Perform  first-difference  filter  operation on  series with  missing
          data; return differenced series (no printed output)

          <real> Y(n), YF(n)
          :
          :
          CALL DIFM (Y, YMISS, N, YF, YFMISS, NYF)

                                      ===













                                    <10-5>
1DIFMC:   Perform  user-specified  difference  filter operation  on series with
          missing data; return differenced series (no printed output)

          INTEGER ND(nfac), IOD(nfac)
          <real> Y(n), YF(n), PHI(iar)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL DIFMC (Y, YMISS, N,
         +            NFAC, ND, IOD, IAR, PHI, LPHI,
         +            YF, YFMISS, NYF, LDSTAK)

                                      ===


             Subroutines for Computing the Gain and Phase Functions

 GFSLF:   Compute and plot gain function of symmetric linear filter

          <real> H(k)
          :
          :
          CALL GFSLF (H, K)

                                      ===

 GFSLFS:  Compute and optionally  plot gain function of symmetric linear filter
          with user-supplied  control values;  return gain function  values and
          corresponding frequency values

          <real> H(k), GAIN(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL GFSLFS (H, K,
         +             NF, FMIN, FMAX, GAIN, FREQ, NPRT, LDSTAK)

                                      ===

 GFARF:   Compute   and  plot  gain and  phase functions of  autoregressive  or
          difference filter

          <real> PHI(iar)
          :
          :
          CALL GFARF (PHI, IAR)

                                      ===









                                    <10-6>
1GFARFS:  Compute   and  optionally   plot  gain   and   phase   functions   of
          autoregressive  or  difference  filter;  with  user-supplied  control
          values;  return  gain and  phase function  values  and  corresponding
          frequency values

          <real> PHI(iar), GAIN(nf), PHAS(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL GFARFS (PHI, IAR,
         +             NF, FMIN, FMAX, GAIN, PHAS, FREQ, NPRT, LDSTAK)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the argument is  input to the  subroutine and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DSTAK   ... The DOUBLE  PRECISION  vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first LDSTAK locations of DSTAK are overwritten  during subroutine
             execution.

 FC      --> The  cutoff  frequency  for  the  filter,  in  cycles  per  sample
             interval.  FC must lie between 0.0 and 0.5.

 FMAX    --> The maximum frequency, in cycles per sample interval, at which the
             gain and phase functions are computed (0.0 <= MIN < FMAX <=  0.5).
             The  default  value  is 0.5.  If FMAX is outside the range FMIN to
             0.5 or is not an argument in the CALL statement the default  value
             is used.

 FMIN    --> The minimum frequency, in cycles per sample interval, at which the
             gain and phase functions are computed (0.0 <= FMIN < FMAX <= 0.5).
             The  default  value  is  0.5.  If FMIN is outside the range 0.0 to
             FMAX or is not an argument in the CALL statement the default value
             is used.

 FREQ    <-- The vector of dimension at least NF containing the  NF frequencies
             at which the gain and phase functions are computed.








                                    <10-7>
1GAIN    <-- The  vector of  dimension  at  least  NF containing  the  NF  gain
             function values over the range FMIN to FMAX.

             For symmetric linear filters:

               The gain function of a symmetric linear filter is

                            K
               GAIN(I) = | SUM H(J) * cos[2*pi*|Km-J|*(FMIN+DELi)] |
                           J=1

               for I = 1, ..., NF, where

               Km is the midpoint of the symmetric filter, Km = (K+1)/2, and

               DELi is the frequency increment, defined as

                 DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1)

               There is no phase change in a symmetric linear filter.

             For autoregressive (or difference) filters:

               The  gain   and  phase   functions  of  an   autoregressive  (or
               difference) filter are

               GAIN(I) =

                        IAR
                   |1 - SUM PHI(J) * cos[2*pi*J*(FMIN+DELi)]
                        J=1

                            IAR
                      - i * SUM PHI(J) * sin[2*pi*J*(FMIN+DELi)] |
                            J=1

               and

                                  IAR
                                  SUM PHI(J) * sin[i*pi*J*(FMIN+DELi)]
                                  J=1
               PHAS(I) = Arctan ----------------------------------------
                                    IAR
                                1 - SUM PHI(J) * cos[i*pi*J*(FMIN+DELi)]
                                    J=1

               for I = 1, 2, ..., NF, where

               i  is the complex value sqrt(-1); and

               DELi is the frequency increment, defined as

                 DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1).

 H       --> The vector of dimension at least K containing filter coefficients,
             which must be symmetric about H[(K+1)/2].



                                    <10-8>
1HHP     <-- The  vector of dimension  at least  K containing  the  K high-pass
             filter  coefficients, which are symmetric about HHP[(K+1)/2].  The
             high-pass  filter  coefficients  are computed  from  the  low-pass
             coefficients by

                            HLP(J)
             HHP(J) = 1 - ----------   for J = Km,
                           K
                          SUM HLP(J)
                          L=1

                            HLP(J)
             HHP(J) =   - ----------   for J = 1, ..., Km-1, Km+1, ..., K.
                           K
                          SUM HLP(J)
                          L=1

 HLP     --- The  vector of  dimension  at least  K containing  the  K low-pass
             filter  coefficients, which must be symmetric  about HLP[(K+1)/2].
             HLP must be input to HPCOEF; it is returned by LPCOEF  and LOPASS.

             For LPCOEF and LOPASS, HLP is defined by

                            K
             HLP(J) = hJ / SUM hI for J = 1, ..., K,
                           I=1

             where

             hJ is computed by

               hJ = 2*FC                                        for J = Km

                     sin[2*pi*|Km-J|*FC]   sin[2*pi*|Km-J|/K]   for j = 1, ...,
               hJ =  ------------------- * ------------------        Km-1,Km+1,
                         2*pi*|Km-J|          2*pi*|Km-J|/K          ..., K,

             with Km the midpoint of the filter, Km = (K+1)/2.

             This low-pass filter  has a  transfer function which  changes from
             approximately  one to zero  in a  transition band about  the ideal
             cutoff frequency FC,  that is  from (FC - 1/K)  to  (FC + 1/K), as
             shown in figure 6.7 of Bloomfield [1976].

 IAR     --- The number  of coefficients in the autoregressive  (or difference)
             filter,  including zero  coefficients.   Equivalently, IAR  is the
             maximum order  of the backward shift operator.   IAR must be input
             to ARFLT, GFARF and GFARFS; it is returned by DIFC and DIFMC.  IAR
             is defined by

                   NDF
             IAR = SUM IOD(J) * ND(J).
                   J=1






                                    <10-9>
1IERR    ... An  error  flag  returned  in  COMMON  /ERRCHK/   [see  chapter 1,
             section D.5].  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 IOD     --> The vector of dimension  at least NFAC containing the  NFAC values
             designating the order of each difference factor.

 K       --> The number of coefficients in the symmetric linear filter.  K must
             be odd.   For LPCOEF, LOPASS, and HIPASS, the user must choose the
             number   of   filter  terms,   K,  so   that   (FC - 1/K) > 0  and
             (FC + 1/K) < 0.5.   In addition, K must be chosen as  a compromise
             between:

             1)  A sharp cutoff, that is, 1/K small; and

             2)  Minimizing the number  of data  points lost  by  the filtering
                 operations. ((K-1)/2 data points will be lost from each end of
                 the series.)

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision  version of  STARPAC is  being  used P = 0.5,
             otherwise P = 1.0 [see chapter 1, section B].

             For DIFC and DIFMC:

                                 NFAC
               LDSTAK >= 7 + 2 * SUM ND(J)*IOD(J)*P
                                 J=1

             For GFSLFS and GFARFS:

               LDSTAK >= (11+IO*(9+NF))/2 + IO*2*NF*P

               where IO = 0 if NPRT = 0, and IO = 1 if NPRT <> 0.

 LPHI    --> The length of the vector PHI.  LPHI must equal or exceed IAR.

 N       --> The number of observations in the time series.  The minimum number
             of observations is three.

 ND      --> The vector of dimension  at least NFAC containing the  NFAC values
             designating the number of  times each  difference factor is  to be
             applied.

 NF      --> The number of frequencies  at which  the gain and  phase functions
             are  to be computed.  The default  value is 101.   If NF is not an
             argument of the  subroutine CALL  statement the  default  value is
             used.

 NFAC    --> The number of difference factors.



                                    <10-10>
1NPRT    --> The argument controlling printed output.

             If NPRT < 0, the output consists  of a  plot of the  gain function
                          versus  frequency,   where  the   gain   function  is
                          expressed in  decibels and  is adjusted  so  that the
                          peak  is at zero.   For GFARFS  only  the output also
                          includes   a  plot  of  the  phase   function  versus
                          frequency.

             If NPRT = 0, the automatic printout is suppressed.

             If NPRT > 0, the output consists of a log-linear plot of  the gain
                          function versus  frequency.   For  GFARFS  only,  the
                          output also  includes  a  plot of the  phase function
                          versus frequency.

             The  default value  is -1.  If  NPRT is  not an  argument  of  the
             subroutine CALL statement the default value is used.

 NS      --> The sample rate, 1 <= NS <= N.

 NYF     <-- The number of observations in YF.

 NYS     <-- The number of observations in YS.

 PHAS    <-- The  vector of  dimension  at least  NF containing  the  NF  phase
             function values over the range FMIN to FMAX.

 PHI     --- The vector of dimension at least NF containing  IAR autoregressive
             or  difference filter coefficients.   PHI must  be input to ARFLT,
             GFARF, and GFARFS; it is returned by DIFC and DIFMC.

             For DIFC and DIFMC the difference filter coefficients are obtained
             by expanding the difference operator

               NFAC
             PRODUCT (1-B[IOD(J)])**ND(J) =
               J=1

                      1 - PHI(1)*B[1] - PHI(2)*B[2] - ... - PHI(IAR)*B[IAR]

             where

             B[i]   is the backward shift operator, defined by

                    B[i]*y(t) = y(t-i);

             PHI(i) is the ith difference filter coefficient,  which will be  a
                    positive  or  negative  integer  if  the ith order backward
                    shift operator B[i] is used,  and zero  if  the  ith  order
                    backward shift operator is unused.

 Y       --> The vector of dimension  at least N containing the  N observations
             of the time series.





                                    <10-11>
1YF      <-- The vector of dimension at  least N  containing the NYF  values of
             the  filtered  series.    The  filtered   series  will   start  in
             YF(1); YF(NYF+1) through YF(N) will be set to zero.

             For symmetric linear filters:

               The filtered series obtained by applying a moving average filter
               is defined by

                        K
               YF(I) = SUM H(J) * Y(I+K-J) for I = 1, ..., NYF,
                       J=1

               where

               NYF is   the   number  of   values  in   the   filtered  series,
                   NYF = N - (K-1), reflecting  the (K-1)/2  data  points  lost
                   from each  end  of  the  original series  by  the  filtering
                   operation.

             For autoregressive or difference filters:

               The  filtered  series  obtained  using  an   autoregressive  (or
               difference) filter is computed by

                                  IAR
               YF(I) = Z(I+IAR) - SUM PHI(J)*Z(I+IAR-J)  for I = 1, ..., NYF,
                                  J=1

               where

               Z   is the N  observation time series being filtered  which, for
                   ARFLT is  the input series  Y minus  its mean and,  for DIF,
                   DIFC, DIFM and DIFMC is the input series Y;

               NYF is  the  number  of  observations in  the  filtered  series,
                   NYF = N-IAR, reflecting  the IAR  data points lost  from the
                   beginning of the original series by the filtering operation.

 YFMISS  <-- The  missing value  code  used within  the filtered  series  YF to
             indicate that a value could not be computed due to missing data.

 YMISS   --> The missing value code used within the input series Y  to indicate
             that an observation is missing.

 YS      <-- The vector of dimension at least N containing the series formed by
             sampling every NSth element of Y,

             YS(J) = Y((J-1)*NS + 1) for J = 1, ..., NYS,

             where NYS, the number  of observations  in the sampled  series, is
             returned  by subroutine SAMPLE.   The series  will start in YS(1);
             YS(NYS+1) through YS(N) will be set to zero.






                                    <10-12>
1E.  Computational Methods

 E.1  Algorithms

      The  code for  computing  the low-pass  filter coefficients  is  based on
 subroutine LOPASS, listed on  page 149  of Bloomfield [1976].   The transforms
 used to compute the gain function of symmetric filters and the gain  and phase
 functions autoregressive (or  difference) filters are based on  the algorithms
 shown on pages 311 and 420, respectively, of Jenkins and Watts [1968].


 E.2  Computed Results and Printed Output

      Except  for  the  gain  and  phase function subroutines,  STARPAC digital
 filtering subroutines do not produce printed output.  For the gain  and  phase
 function subroutines, the argument controlling the printed results is NPRT and
 is  discussed  in  section  D;  the  output  from  the gain and phase function
 subroutines consists of line printer plots of the gain and phase  function  of
 the input filter.


 F.  Examples

      In  the example program below,  DIF is used to filter the input series Y;
 VP (documented in chapter 2) is used to display the log of the original series
 and the differenced series;  and GFARF is used to  plot  the  gain  and  phase
 functions  of  the  first  difference  filter.  The  data used are the natural
 logarithm of Series G, the airline data, listed on page 531 of Box and Jenkins
 [1976].  The formulas for the gain and phase functions of a  first  difference
 filter  and  a  plot of the corresponding gain function are shown on pages 296
 and 9, respectively, of Jenkins and Watts [1968].




























                                    <10-13>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE DIF AND GFARF USING SINGLE PRECISION VERSION OF
 C     STARPAC
 C
 C     N.B. DECLARATION OF Y, YF AND PHI MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(200), YF(200), PHI(5)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     COMPUTE LOG OF DATA
 C
       DO 10 I = 1, N
         Y(I) = ALOG(Y(I))
    10 CONTINUE
 C
 C     CALL DIF TO PERFORM DIFFERENCING OPERATION
 C
       CALL DIF (Y, N, YF, NYF)
 C
 C     PRINT TITLE AND CALL VP TO DISPLAY LOG OF ORIGINAL SERIES
 C
       WRITE (IPRT,102)
       CALL VP (Y, N, 1)
 C
 C     PRINT TITLE AND CALL VP TO DISPLAY DIFFERENCED SERIES
 C
       WRITE (IPRT,103)
       CALL VP (YF, NYF, 1)
 C
 C     SET PARAMETERS FOR FIRST DIFFERENCE FILTER
 C
       PHI(1) = 1.0
       IAR = 1
 C
 C     PRINT TITLE AND CALL GFARF TO COMPUTE GAIN AND PHASE OF
 C     FIRST DIFFERENCE FILTER
 C
       WRITE (IPRT,104)
       CALL GFARF (PHI, IAR)
 C
 C     FORMAT STATEMENTS
 C

                                    <10-14>
1  100 FORMAT (4I5)
   101 FORMAT (12F6.1)
   102 FORMAT ('1LOG OF ORIGINAL SERIES DISPLAYED WITH STARPAC PLOT',
      1  ' SUBROUTINE VP')
   103 FORMAT ('1RESULTS OF STARPAC FIRST DIFFERENCE DIGITAL FILTERING',
      1  ' SUBROUTINE DIF DISPLAYED WITH STARPAC PLOT SUBROUTINE VP')
   104 FORMAT ('1RESULTS OF STARPAC',
      *  ' GAIN AND PHASE FUNCTION SUBROUTINE GFARF')
       END

 Data:

   144
  112.0 118.0 132.0 129.0 121.0 135.0 148.0 148.0 136.0 119.0 104.0 118.0
  115.0 126.0 141.0 135.0 125.0 149.0 170.0 170.0 158.0 133.0 114.0 140.0
  145.0 150.0 178.0 163.0 172.0 178.0 199.0 199.0 184.0 162.0 146.0 166.0
  171.0 180.0 193.0 181.0 183.0 218.0 230.0 242.0 209.0 191.0 172.0 194.0
  196.0 196.0 236.0 235.0 229.0 243.0 264.0 272.0 237.0 211.0 180.0 201.0
  204.0 188.0 235.0 227.0 234.0 264.0 302.0 293.0 259.0 229.0 203.0 229.0
  242.0 233.0 267.0 269.0 270.0 315.0 364.0 347.0 312.0 274.0 237.0 278.0
  284.0 277.0 317.0 313.0 318.0 374.0 413.0 405.0 355.0 306.0 271.0 306.0
  315.0 301.0 356.0 348.0 355.0 422.0 465.0 467.0 404.0 347.0 305.0 336.0
  340.0 318.0 362.0 348.0 363.0 435.0 491.0 505.0 404.0 359.0 310.0 337.0
  360.0 342.0 406.0 396.0 420.0 472.0 548.0 559.0 463.0 407.0 362.0 405.0
  417.0 391.0 419.0 461.0 472.0 535.0 622.0 606.0 508.0 461.0 390.0 432.0


































                                    <10-15>
1LOG OF ORIGINAL SERIES DISPLAYED WITH STARPAC PLOT SUBROUTINE VP
                                                                                                         STARPAC 2.08S (03/15/90)
              4.6444    4.8232    5.0021    5.1810    5.3598    5.5387    5.7175    5.8964    6.0752    6.2541    6.4329
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I    +                                                                                                I   4.7185
  2.0000      I       +                                                                                             I   4.7707
  3.0000      I             +                                                                                       I   4.8828
  4.0000      I            +                                                                                        I   4.8598
  5.0000      I        +                                                                                            I   4.7958
  6.0000      I               +                                                                                     I   4.9053
  7.0000      I                    +                                                                                I   4.9972
  8.0000      I                    +                                                                                I   4.9972
  9.0000      I               +                                                                                     I   4.9127
  10.000      I        +                                                                                            I   4.7791
  11.000      I+                                                                                                    I   4.6444
  12.000      I       +                                                                                             I   4.7707
  13.000      I      +                                                                                              I   4.7449
  14.000      I           +                                                                                         I   4.8363
  15.000      I                 +                                                                                   I   4.9488
  16.000      I               +                                                                                     I   4.9053
  17.000      I          +                                                                                          I   4.8283
  18.000      I                    +                                                                                I   5.0039
  19.000      I                           +                                                                         I   5.1358
  20.000      I                           +                                                                         I   5.1358
  21.000      I                       +                                                                             I   5.0626
  22.000      I              +                                                                                      I   4.8903
  23.000      I     +                                                                                               I   4.7362
  24.000      I                 +                                                                                   I   4.9416
  25.000      I                   +                                                                                 I   4.9767
  26.000      I                    +                                                                                I   5.0106
  27.000      I                              +                                                                      I   5.1818
  28.000      I                         +                                                                           I   5.0938
  29.000      I                            +                                                                        I   5.1475
  30.000      I                              +                                                                      I   5.1818
  31.000      I                                    +                                                                I   5.2933
  32.000      I                                    +                                                                I   5.2933
  33.000      I                                +                                                                    I   5.2149
  34.000      I                         +                                                                           I   5.0876
  35.000      I                   +                                                                                 I   4.9836
  36.000      I                          +                                                                          I   5.1120
  37.000      I                            +                                                                        I   5.1417
  38.000      I                               +                                                                     I   5.1930
  39.000      I                                   +                                                                 I   5.2627
  40.000      I                               +                                                                     I   5.1985
  41.000      I                                +                                                                    I   5.2095
  42.000      I                                         +                                                           I   5.3845
  43.000      I                                            +                                                        I   5.4381
  44.000      I                                               +                                                     I   5.4889
  45.000      I                                       +                                                             I   5.3423
  46.000      I                                  +                                                                  I   5.2523
  47.000      I                            +                                                                        I   5.1475
  48.000      I                                   +                                                                 I   5.2679
  49.000      I                                   +                                                                 I   5.2781
  50.000      I                                   +                                                                 I   5.2781
  51.000      I                                              +                                                      I   5.4638
  52.000      I                                              +                                                      I   5.4596
  53.000      I                                            +                                                        I   5.4337
  54.000      I                                               +                                                     I   5.4931
  55.000      I                                                    +                                                I   5.5759
                                                             <10-16>
1 56.000      I                                                      +                                              I   5.6058
  57.000      I                                              +                                                      I   5.4681
  58.000      I                                        +                                                            I   5.3519
  59.000      I                               +                                                                     I   5.1930
  60.000      I                                     +                                                               I   5.3033
  61.000      I                                      +                                                              I   5.3181
  62.000      I                                 +                                                                   I   5.2364
  63.000      I                                              +                                                      I   5.4596
  64.000      I                                            +                                                        I   5.4250
  65.000      I                                             +                                                       I   5.4553
  66.000      I                                                    +                                                I   5.5759
  67.000      I                                                            +                                        I   5.7104
  68.000      I                                                          +                                          I   5.6802
  69.000      I                                                   +                                                 I   5.5568
  70.000      I                                            +                                                        I   5.4337
  71.000      I                                     +                                                               I   5.3132
  72.000      I                                            +                                                        I   5.4337
  73.000      I                                               +                                                     I   5.4889
  74.000      I                                             +                                                       I   5.4510
  75.000      I                                                     +                                               I   5.5872
  76.000      I                                                     +                                               I   5.5947
  77.000      I                                                     +                                               I   5.5984
  78.000      I                                                              +                                      I   5.7526
  79.000      I                                                                      +                              I   5.8972
  80.000      I                                                                   +                                 I   5.8493
  81.000      I                                                             +                                       I   5.7430
  82.000      I                                                      +                                              I   5.6131
  83.000      I                                              +                                                      I   5.4681
  84.000      I                                                       +                                             I   5.6276
  85.000      I                                                        +                                            I   5.6490
  86.000      I                                                       +                                             I   5.6240
  87.000      I                                                              +                                      I   5.7589
  88.000      I                                                              +                                      I   5.7462
  89.000      I                                                              +                                      I   5.7621
  90.000      I                                                                        +                            I   5.9243
  91.000      I                                                                             +                       I   6.0234
  92.000      I                                                                            +                        I   6.0039
  93.000      I                                                                     +                               I   5.8721
  94.000      I                                                            +                                        I   5.7236
  95.000      I                                                      +                                              I   5.6021
  96.000      I                                                            +                                        I   5.7236
  97.000      I                                                              +                                      I   5.7526
  98.000      I                                                           +                                         I   5.7071
  99.000      I                                                                     +                               I   5.8749
  100.00      I                                                                    +                                I   5.8522
  101.00      I                                                                     +                               I   5.8721
  102.00      I                                                                              +                      I   6.0450
  103.00      I                                                                                    +                I   6.1420
  104.00      I                                                                                    +                I   6.1463
  105.00      I                                                                            +                        I   6.0014
  106.00      I                                                                   +                                 I   5.8493
  107.00      I                                                            +                                        I   5.7203
  108.00      I                                                                  +                                  I   5.8171
  109.00      I                                                                  +                                  I   5.8289
  110.00      I                                                              +                                      I   5.7621
  111.00      I                                                                      +                              I   5.8916
  112.00      I                                                                    +                                I   5.8522
  113.00      I                                                                      +                              I   5.8944
  114.00      I                                                                                +                    I   6.0753
                                                             <10-17>
1 115.00      I                                                                                       +             I   6.1964
  116.00      I                                                                                        +            I   6.2246
  117.00      I                                                                            +                        I   6.0014
  118.00      I                                                                     +                               I   5.8833
  119.00      I                                                             +                                       I   5.7366
  120.00      I                                                                  +                                  I   5.8201
  121.00      I                                                                     +                               I   5.8861
  122.00      I                                                                   +                                 I   5.8348
  123.00      I                                                                            +                        I   6.0064
  124.00      I                                                                           +                         I   5.9814
  125.00      I                                                                              +                      I   6.0403
  126.00      I                                                                                     +               I   6.1570
  127.00      I                                                                                             +       I   6.3063
  128.00      I                                                                                              +      I   6.3261
  129.00      I                                                                                   +                 I   6.1377
  130.00      I                                                                            +                        I   6.0088
  131.00      I                                                                      +                              I   5.8916
  132.00      I                                                                            +                        I   6.0039
  133.00      I                                                                              +                      I   6.0331
  134.00      I                                                                          +                          I   5.9687
  135.00      I                                                                              +                      I   6.0379
  136.00      I                                                                                   +                 I   6.1334
  137.00      I                                                                                     +               I   6.1570
  138.00      I                                                                                            +        I   6.2823
  139.00      I                                                                                                    +I   6.4329
  140.00      I                                                                                                   + I   6.4069
  141.00      I                                                                                         +           I   6.2305
  142.00      I                                                                                   +                 I   6.1334
  143.00      I                                                                          +                          I   5.9661
  144.00      I                                                                                +                    I   6.0684





























                                                             <10-18>
1RESULTS OF STARPAC FIRST DIFFERENCE DIGITAL FILTERING SUBROUTINE DIF DISPLAYED WITH STARPAC PLOT SUBROUTINE VP
                                                                                                         STARPAC 2.08S (03/15/90)
              -.2231    -.1785    -.1339    -.0893    -.0446     .0000     .0446     .0893     .1339     .1785     .2231
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                                                              +                                      I   .52186E-01
  2.0000      I                                                                           +                         I   .11212
  3.0000      I                                             +                                                       I  -.22990E-01
  4.0000      I                                    +                                                                I  -.64022E-01
  5.0000      I                                                                           +                         I   .10948
  6.0000      I                                                                       +                             I   .91937E-01
  7.0000      I                                                  +                                                  I  0.
  8.0000      I                               +                                                                     I  -.84557E-01
  9.0000      I                    +                                                                                I  -.13353
  10.000      I                    +                                                                                I  -.13473
  11.000      I                                                                              +                      I   .12629
  12.000      I                                            +                                                        I  -.25752E-01
  13.000      I                                                                      +                              I   .91350E-01
  14.000      I                                                                           +                         I   .11248
  15.000      I                                        +                                                            I  -.43485E-01
  16.000      I                                 +                                                                   I  -.76961E-01
  17.000      I                                                                                         +           I   .17563
  18.000      I                                                                                +                    I   .13185
  19.000      I                                                  +                                                  I  0.
  20.000      I                                  +                                                                  I  -.73203E-01
  21.000      I           +                                                                                         I  -.17225
  22.000      I               +                                                                                     I  -.15415
  23.000      I                                                                                                +    I   .20544
  24.000      I                                                          +                                          I   .35091E-01
  25.000      I                                                          +                                          I   .33902E-01
  26.000      I                                                                                        +            I   .17115
  27.000      I                              +                                                                      I  -.88033E-01
  28.000      I                                                              +                                      I   .53744E-01
  29.000      I                                                          +                                          I   .34289E-01
  30.000      I                                                                           +                         I   .11152
  31.000      I                                                  +                                                  I  0.
  32.000      I                                +                                                                    I  -.78369E-01
  33.000      I                     +                                                                               I  -.12734
  34.000      I                           +                                                                         I  -.10399
  35.000      I                                                                               +                     I   .12838
  36.000      I                                                         +                                           I   .29676E-01
  37.000      I                                                             +                                       I   .51293E-01
  38.000      I                                                                  +                                  I   .69733E-01
  39.000      I                                    +                                                                I  -.64193E-01
  40.000      I                                                    +                                                I   .10989E-01
  41.000      I                                                                                         +           I   .17501
  42.000      I                                                              +                                      I   .53584E-01
  43.000      I                                                             +                                       I   .50858E-01
  44.000      I                 +                                                                                   I  -.14660
  45.000      I                              +                                                                      I  -.90061E-01
  46.000      I                           +                                                                         I  -.10478
  47.000      I                                                                             +                       I   .12036
  48.000      I                                                    +                                                I   .10257E-01
  49.000      I                                                  +                                                  I  0.
  50.000      I                                                                                            +        I   .18572
  51.000      I                                                 +                                                   I  -.42463E-02
  52.000      I                                            +                                                        I  -.25864E-01
  53.000      I                                                               +                                     I   .59339E-01
  54.000      I                                                                     +                               I   .82888E-01
  55.000      I                                                         +                                           I   .29853E-01
                                                             <10-19>
1 56.000      I                   +                                                                                 I  -.13774
  57.000      I                        +                                                                            I  -.11620
  58.000      I              +                                                                                      I  -.15890
  59.000      I                                                                           +                         I   .11035
  60.000      I                                                     +                                               I   .14815E-01
  61.000      I                                +                                                                    I  -.81678E-01
  62.000      I                                                                                                    +I   .22314
  63.000      I                                          +                                                          I  -.34635E-01
  64.000      I                                                         +                                           I   .30371E-01
  65.000      I                                                                             +                       I   .12063
  66.000      I                                                                                +                    I   .13448
  67.000      I                                           +                                                         I  -.30254E-01
  68.000      I                      +                                                                              I  -.12334
  69.000      I                      +                                                                              I  -.12311
  70.000      I                       +                                                                             I  -.12052
  71.000      I                                                                             +                       I   .12052
  72.000      I                                                              +                                      I   .55216E-01
  73.000      I                                          +                                                          I  -.37899E-01
  74.000      I                                                                                 +                   I   .13621
  75.000      I                                                    +                                                I   .74627E-02
  76.000      I                                                   +                                                 I   .37106E-02
  77.000      I                                                                                     +               I   .15415
  78.000      I                                                                                  +                  I   .14458
  79.000      I                                       +                                                             I  -.47829E-01
  80.000      I                          +                                                                          I  -.10632
  81.000      I                     +                                                                               I  -.12988
  82.000      I                 +                                                                                   I  -.14507
  83.000      I                                                                                      +              I   .15956
  84.000      I                                                       +                                             I   .21353E-01
  85.000      I                                            +                                                        I  -.24957E-01
  86.000      I                                                                                +                    I   .13488
  87.000      I                                               +                                                     I  -.12699E-01
  88.000      I                                                      +                                              I   .15848E-01
  89.000      I                                                                                      +              I   .16220
  90.000      I                                                                        +                            I   .99192E-01
  91.000      I                                              +                                                      I  -.19561E-01
  92.000      I                    +                                                                                I  -.13177
  93.000      I                 +                                                                                   I  -.14853
  94.000      I                       +                                                                             I  -.12147
  95.000      I                                                                             +                       I   .12147
  96.000      I                                                        +                                            I   .28988E-01
  97.000      I                                        +                                                            I  -.45462E-01
  98.000      I                                                                                        +            I   .16782
  99.000      I                                             +                                                       I  -.22728E-01
  100.00      I                                                      +                                              I   .19915E-01
  101.00      I                                                                                         +           I   .17289
  102.00      I                                                                        +                            I   .97032E-01
  103.00      I                                                   +                                                 I   .42919E-02
  104.00      I                  +                                                                                  I  -.14491
  105.00      I                +                                                                                    I  -.15209
  106.00      I                     +                                                                               I  -.12901
  107.00      I                                                                        +                            I   .96799E-01
  108.00      I                                                     +                                               I   .11834E-01
  109.00      I                                   +                                                                 I  -.66894E-01
  110.00      I                                                                               +                     I   .12959
  111.00      I                                         +                                                           I  -.39442E-01
  112.00      I                                                           +                                         I   .42200E-01
  113.00      I                                                                                           +         I   .18094
  114.00      I                                                                             +                       I   .12110
                                                             <10-20>
1 115.00      I                                                        +                                            I   .28114E-01
  116.00      I+                                                                                                    I  -.22314
  117.00      I                        +                                                                            I  -.11809
  118.00      I                 +                                                                                   I  -.14675
  119.00      I                                                                     +                               I   .83511E-01
  120.00      I                                                                 +                                   I   .66021E-01
  121.00      I                                       +                                                             I  -.51293E-01
  122.00      I                                                                                        +            I   .17154
  123.00      I                                            +                                                        I  -.24939E-01
  124.00      I                                                               +                                     I   .58841E-01
  125.00      I                                                                            +                        I   .11672
  126.00      I                                                                                   +                 I   .14930
  127.00      I                                                      +                                              I   .19874E-01
  128.00      I        +                                                                                            I  -.18842
  129.00      I                     +                                                                               I  -.12891
  130.00      I                        +                                                                            I  -.11717
  131.00      I                                                                           +                         I   .11224
  132.00      I                                                         +                                           I   .29199E-01
  133.00      I                                    +                                                                I  -.64379E-01
  134.00      I                                                                 +                                   I   .69163E-01
  135.00      I                                                                       +                             I   .95527E-01
  136.00      I                                                       +                                             I   .23581E-01
  137.00      I                                                                              +                      I   .12529
  138.00      I                                                                                    +                I   .15067
  139.00      I                                            +                                                        I  -.26060E-01
  140.00      I          +                                                                                          I  -.17640
  141.00      I                            +                                                                        I  -.97083E-01
  142.00      I             +                                                                                       I  -.16725
  143.00      I                                                                         +                           I   .10228






























                                                             <10-21>
1RESULTS OF STARPAC GAIN AND PHASE FUNCTION SUBROUTINE GFARF
                                                                                                         STARPAC 2.08S (03/15/90)
 GAIN FUNCTION OF   1 TERM AUTOREGRESSIVE, OR DIFFERENCE, FILTER
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .0000 -                                                                                   +++++++++++++++++++ -
                I                                                                      +++++++++++++                    I
                I                                                              ++++++++                                 I
                I                                                       +++++++                                         I
                I                                                  +++++                                                I
        -1.8039 -                                             +++++                                                     -
                I                                         ++++                                                          I
                I                                      +++                                                              I
                I                                  ++++                                                                 I
                I                                ++                                                                     I
        -3.6078 -                             +++                                                                       -
                I                           ++                                                                          I
                I                         ++                                                                            I
                I                       ++                                                                              I
                I                     ++                                                                                I
        -5.4117 -                   ++                                                                                  -
                I                  +                                                                                    I
                I                 +                                                                                     I
                I               ++                                                                                      I
                I              +                                                                                        I
        -7.2156 -             +                                                                                         -
                I            +                                                                                          I
                I           +                                                                                           I
                I                                                                                                       I
                I          +                                                                                            I
        -9.0195 -         +                                                                                             -
                I                                                                                                       I
                I        +                                                                                              I
                I       +                                                                                               I
                I                                                                                                       I
       -10.8234 -                                                                                                       -
                I      +                                                                                                I
                I                                                                                                       I
                I     +                                                                                                 I
                I                                                                                                       I
       -12.6273 -                                                                                                       -
                I                                                                                                       I
                I    +                                                                                                  I
                I                                                                                                       I
                I                                                                                                       I
       -14.4312 -                                                                                                       -
                I                                                                                                       I
                I   +                                                                                                   I
                I                                                                                                       I
                I                                                                                                       I
       -16.2351 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
       -18.0390 -  +                                                                                                    -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <10-22>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 PHASE FUNCTION OF   1 TERM AUTOREGRESSIVE, OR DIFFERENCE, FILTER
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         3.1416 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         2.5133 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         1.8850 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         1.2566 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .6283 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .0000 -                                                                                                 ++++  -
                I                                                                                         ++++++++      I
                I                                                                                 ++++++++              I
                I                                                                         ++++++++                      I
                I                                                                 ++++++++                              I
         -.6283 -                                                         ++++++++                                      -
                I                                                  +++++++                                              I
                I                                          ++++++++                                                     I
                I                                  ++++++++                                                             I
                I                          ++++++++                                                                     I
        -1.2566 -                  ++++++++                                                                             -
                I          ++++++++                                                                                     I
                I  ++++++++                                                                                             I
                I +                                                                                                     I
                I                                                                                                       I
        -1.8850 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -2.5133 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -3.1416 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <10-23>
1G.  Acknowledgments

      The  code for  computing  the low-pass  filter coefficients  is  based on
 subroutine LOPASS, listed on  page 149  of Bloomfield [1976].   The transforms
 used to compute the gain function of symmetric filters and the gain  and phase
 functions  of  autoregressive  (or  difference)  filters  are  based   on  the
 algorithms shown  on pages 311  and 420,  respectively, of  Jenkins  and Watts
 [1968].



















































                                    <10-24>
1-----                             CHAPTER 11                             -----

                              COMPLEX DEMODULATION


 A.  Introduction

      STARPAC contains  two  subroutines which  find the  amplitude  and  phase
 functions of a demodulated series as described in Bloomfield [1976].

      The demodulated series w(t) is formed by multiplying the  observed series
 by a complex sinusoid at the demodulation frequency.  If the observed series Y
 is a sinusoid of the nominal demodulation frequency FD with  varying amplitude
 and phase plus noise, that is,

            Y(t) = R(t) * cos[2*pi*FD*t + phi(t)] + a(t)

                 = 0.5 * R(t) * (exp[i*(2*pi*FD*t+phi(t))] +

                                 exp[-i*(2*pi*FD*t+phi(t))]) + a(t)

 then the demodulated series may be represented by

    w(t) = exp[-i*2*pi*FD] * Y(t)

         = 0.5 * R(t) * exp[i*phi(t)] +

           0.5 * R(t) * exp[-i*(4*pi*FD*t+phi(t))] + exp[-i*2*pi*FD] * a(t)

 for t = 1, ..., N, where

 N       is the number of time points in the series;

 i       is the complex value sqrt(-1);

 FD      is the demodulation frequency in cycles per sample interval;

 R(t)    is the amplitude component of the observed series at time t;

 phi(t)  is the phase component of the observed series at time t;

 a(t)    is the noise at time t.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      Subroutine  DEMOD computes the smoothed amplitude and phase components of
 the demodulated series.  The user must  specify  the  demodulation  frequency,
 along  with  the  number of filter terms and the cutoff frequency which define
 the low-pass filter utilized to smooth the  demodulated  series.  Output  from
 DEMOD  consists  of  plots  of  the  amplitude and phase functions.  The phase
 function plot reduces discontinuities using the method suggested by Bloomfield

                                    <11-1>
1[1976].  As shown  in  the example provided in section F, this method displays
 both the principle phase value,  which is defined to lie in the range  -pi  to
 pi, and the principle phase value plus or minus 2*pi, where the sign is chosen
 such that the second value lies in the range -2*pi to 2*pi.

      Subroutine DEMODS is the same as DEMOD except that the computed amplitude
 and phase functions are returned to the user and the printed  output described
 for DEMOD is optional.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The conventions used to present the following  declaration  and
 CALL statements are given in chapter 1, sections B and D.

 DEMOD:   Compute  and plot the results of a complex  demodulation of the input
          series

          <real> Y(n)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL DEMOD (Y, N, FD, FC, K, LDSTAK)

                                      ===

 DEMODS:  Compute and optionally  plot the results of a complex demodulation of
          the input series; return amplitude and phase functions of demodulated
          series

          <real> Y(n), AMPL(n), PHAS(n)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL DEMODS (Y, N, FD, FC, K,
         +             AMPL, PHAS, NDEM, NPRT, LDSTAK)

                                      ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the argument  is input to the subroutine and  that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 AMPL    <-- The vector  of  dimension  at least N-(K-1) that contains the NDEM
             values of the smoothed amplitude function of the observed  series,

                                    <11-2>
1            AMPL(I) = R(t), where R(t) is defined in section A and the index I
             is  computed as I = t - (K-1)/2 for t = (K+1)/2 to N-(K-1)/2.  The
             stored values of the amplitude function  will  start  in  AMPL(1);
             AMPL(NDEM+1) to AMPL(N) will be set to zero.

 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 FC      --> The cutoff frequency for the low-pass filter in cycles  per sample
             interval.  FC must lie between 1/K and FD.

 FD      --> The demodulation frequency in cycles per sample interval.  FD must
             lie between 0.0 and 0.5.

 IERR    ... An  error  flag  returned  in  COMMON  /ERRCHK/   [see  chapter 1,
             section D.5].  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 K       --> The number  of terms in  the low-pass  filter used to  extract the
             amplitude and phase  functions.   K  must be odd.   The user  must
             choose the number of  filter terms, K, so that  (FC - 1/K) > 0 and
             (FC + 1/K) < 0.5.   In addition, K must be chosen as  a compromise
             between:

             1)  A sharp cutoff, that is, 1/K small; and

             2)  Minimizing the number  of data  points lost  by  the filtering
                 operations.  ((K - 1)/2 data points will be lost from each end
                 of the series.)

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision  version of  STARPAC is  being  used P = 0.5,
             otherwise P = 1.0 [see chapter 1, section B].

             For DEMOD:

               LDSTAK >= 10 + (3*N+K)*P

             For DEMODS:

               LDSTAK >= 9 + (IO*2*N+K)*P

               where IO = 0 if NPRT = 0 and IO = 1 if NPRT <> 0.

 N       --> The number of observations, which must equal or exceed 17.

 NDEM    <-- The number of observations in AMPL and PHAS, NDEM = N - (K-1).




                                    <11-3>
1NPRT    --> The variable controlling printed output.

             If NPRT = 0, the automatic printout is suppressed.

             If NPRT <> 0, the automatic printout is provided.

             The default  value is  1.   If  NPRT  is not  an argument  of  the
             subroutine CALL statement the default value is used.

 PHAS    <-- The vector of dimension at least NDEM =  N-(K-1) that contains the
             NDEM primary values of the smoothed phase function of the observed
             series, PHAS(I) = phi(t), where phi(t) is defined in section A and
             the  index  I is computed as I = t - (K-1)/2 for t = (K+1)/2 to N-
             (K-1)/2.  The stored values of the phase function  will  start  in
             PHAS(1); PHAS(NDEM+1) to PHAS(N) will be set to zero.

 Y       --> The vector of dimension at least N that contains  the observations
             of the time series.


 E.  Computational Methods

 E.1  Algorithms

      The STARPAC code for performing complex demodulation was adapted from the
 subroutines  given  on  pages  147  to  150 of Bloomfield [1976].  As noted in
 Bloomfield,  the first term of the demodulated series defined in section A  is
 centered  about  zero  frequency while the remaining two terms are centered at
 frequencies FD and 2*FD.  Thus,  the first term  can  be  separated  from  the
 others   using   the   low-pass  filter  described  in  chapter  10  (with  FC
 approximately FD/2), resulting in the complex filtered series

                                   K
                          YF(t) = SUM HLP(J)*w(t+Km-J)
                                  J=1

                                = alpha(t) + i*beta(t)

 which is approximately

           0.5*R(t)*exp[i*phi(t)]   for t = Km, Km+1, ..., (N-Km+1),

 where

 K         is the number of filter terms;

 Km        is the midpoint of the filter, Km = (K+1)/2;

 HLP(J)    is  the  Jth  low-pass  filter  coefficient,  defined in chapter 10,
           section D;

 alpha(t)  is the real part of the filtered series;

 beta(t)   is the imaginary part of the filtered series.

 The smoothed estimates of the amplitude, Rhat,  and phase,  phihat,  functions
 can then be extracted from the filtered series using


                                    <11-4>
1                  Rhat(t) = 2*sqrt[alpha(t)**2 + beta(t)**2]

 and

                     phihat(t) = arctan[alpha(t)/beta(t)].

 Note  that (K-1)/2  points have been  lost from  each end  of  the demodulated
 series by the filtering operation.


 E.2  Computed Results and Printed Output

      The  argument  controlling  the  printed  output,  NPRT,  is discussed in
 section D.  The output consists of plots of the smoothed amplitude  and  phase
 functions,  and  a  list  of the demodulation frequency,  cutoff frequency and
 number of terms in the low-pass filter used to smooth the demodulated series.


 F.  Example

      In the example program below, DEMOD is used to estimate the amplitude and
 phase function corresponding to the input series Y.  The  data  used  are  the
 Wolf  sunspot  numbers  for  the  years 1700 to 1960 as tabulated by Waldmeier
 [l961].  Further discussion of this example can be found on pages 137  to  141
 of Bloomfield [1976].


































                                    <11-5>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE DEMOD USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(300)
       DOUBLE PRECISION DSTAK(500)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 500
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     SET DEMODULATION FREQUENCY
 C         CUTOFF FREQUENCY
 C         NUMBER OF TERMS IN THE LOW-PASS FILTER
 C
       FD = 1.0 / 11.0
       FC = 1.0 / 22.0
       K  = 41
 C
 C     PRINT TITLE AND CALL DEMOD FOR COMPLEX DEMODULATION ANALYSIS
 C
       WRITE (IPRT,102)
       CALL DEMOD (Y, N, FD, FC, K, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (10F7.2)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' COMPLEX DEMODULATION SUBROUTINE DEMOD')
       END








                                    <11-6>
1Data:

   261
    5.00  11.00  16.00  23.00  36.00  58.00  29.00  20.00  10.00   8.00
    3.00   0.00   0.00   2.00  11.00  27.00  47.00  63.00  60.00  39.00
   28.00  26.00  22.00  11.00  21.00  40.00  78.00 122.00 103.00  73.00
   47.00  35.00  11.00   5.00  16.00  34.00  70.00  81.00 111.00 101.00
   73.00  40.00  20.00  16.00   5.00  11.00  22.00  40.00  60.00  80.90
   83.40  47.70  47.80  30.70  12.20   9.60  10.20  32.40  47.60  54.00
   62.90  85.90  61.20  45.10  36.40  20.90  11.40  37.80  69.80 106.10
  100.80  81.60  66.50  34.80  30.60   7.00  19.80  92.50 154.40 125.90
   84.80  68.10  38.50  22.80  10.20  24.10  82.90 132.00 130.90 118.10
   89.90  66.60  60.00  46.90  41.00  21.30  16.00   6.40   4.10   6.80
   14.50  34.00  45.00  43.10  47.50  42.20  28.10  10.10   8.10   2.50
    0.00   1.40   5.00  12.20  13.90  35.40  45.80  41.10  30.10  23.90
   15.60   6.60   4.00   1.80   8.50  16.60  36.30  49.60  64.20  67.00
   70.90  47.80  27.50   8.50  13.20  56.90 121.50 138.30 103.20  85.70
   64.60  36.70  24.20  10.70  15.00  40.10  61.50  98.50 124.70  96.30
   66.60  64.50  54.10  39.00  20.60   6.70   4.30  22.70  54.80  93.80
   95.80  77.20  59.10  44.00  47.00  30.50  16.30   7.30  37.60  74.00
  139.00 111.20 101.60  66.20  44.70  17.00  11.30  12.40   3.40   6.00
   32.30  54.30  59.70  63.70  63.50  52.20  25.40  13.10   6.80   6.30
    7.10  35.60  73.00  85.10  78.00  64.00  41.80  26.20  26.70  12.10
    9.50   2.70   5.00  24.40  42.00  63.50  53.80  62.00  48.50  43.90
   18.60   5.70   3.60   1.40   9.60  47.40  57.10 103.90  80.60  63.60
   37.60  26.10  14.20   5.80  16.70  44.30  63.90  69.00  77.80  64.90
   35.70  21.20  11.10   5.70   8.70  36.10  79.70 114.40 109.60  88.80
   67.80  47.50  30.60  16.30   9.60  33.20  92.60 151.60 136.30 134.70
   83.90  69.40  31.50  13.90   4.40  38.00 141.70 190.20 184.80 159.00
  112.30





























                                    <11-7>
1RESULTS OF STARPAC COMPLEX DEMODULATION SUBROUTINE DEMOD
                                                                                                         STARPAC 2.08S (03/15/90)

 TIME SERIES DEMODULATION

 DEMODULATION FREQUENCY IS .09090909
 CUTOFF FREQUENCY IS       .04545455
 THE NUMBER OF TERMS IN THE FILTER IS    41



 PLOT OF AMPLITUDE OF SMOOTHED DEMODULATED SERIES

 LOCATION OF MEAN IS GIVEN BY PLOT CHARACTER M
             19.7732   24.0727   28.3722   32.6717   36.9711   41.2706   45.5701   49.8696   54.1690   58.4685   62.7680
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                         +                   M                                                       I   30.652
  2.0000      I                              +              M                                                       I   32.848
  3.0000      I                                    +        M                                                       I   35.191
  4.0000      I                                         +   M                                                       I   37.582
  5.0000      I                                             M +                                                     I   39.925
  6.0000      I                                             M      +                                                I   42.146
  7.0000      I                                             M           +                                           I   44.216
  8.0000      I                                             M               +                                       I   46.153
  9.0000      I                                             M                    +                                  I   47.948
  10.000      I                                             M                       +                               I   49.529
  11.000      I                                             M                          +                            I   50.786
  12.000      I                                             M                            +                          I   51.633
  13.000      I                                             M                             +                         I   52.052
  14.000      I                                             M                             +                         I   52.071
  15.000      I                                             M                            +                          I   51.754
  16.000      I                                             M                           +                           I   51.178
  17.000      I                                             M                         +                             I   50.429
  18.000      I                                             M                       +                               I   49.599
  19.000      I                                             M                     +                                 I   48.776
  20.000      I                                             M                    +                                  I   48.011
  21.000      I                                             M                  +                                    I   47.287
  22.000      I                                             M                +                                      I   46.515
  23.000      I                                             M              +                                        I   45.587
  24.000      I                                             M           +                                           I   44.441
  25.000      I                                             M        +                                              I   43.055
  26.000      I                                             M    +                                                  I   41.462
  27.000      I                                             M+                                                      I   39.721
  28.000      I                                          +  M                                                       I   37.931
  29.000      I                                      +      M                                                       I   36.245
  30.000      I                                   +         M                                                       I   34.853
  31.000      I                                 +           M                                                       I   33.877
  32.000      I                               +             M                                                       I   33.306
  33.000      I                               +             M                                                       I   33.020
  34.000      I                              +              M                                                       I   32.861
  35.000      I                              +              M                                                       I   32.690
  36.000      I                             +               M                                                       I   32.398
  37.000      I                            +                M                                                       I   31.913
  38.000      I                           +                 M                                                       I   31.235
  39.000      I                         +                   M                                                       I   30.500
  40.000      I                        +                    M                                                       I   29.924
  41.000      I                       +                     M                                                       I   29.721
  42.000      I                        +                    M                                                       I   30.052
  43.000      I                          +                  M                                                       I   30.936
                                                             <11-8>
1 44.000      I                             +               M                                                       I   32.234
  45.000      I                                +            M                                                       I   33.719
  46.000      I                                    +        M                                                       I   35.148
  47.000      I                                       +     M                                                       I   36.375
  48.000      I                                         +   M                                                       I   37.438
  49.000      I                                            +M                                                       I   38.520
  50.000      I                                             M +                                                     I   39.817
  51.000      I                                             M    +                                                  I   41.445
  52.000      I                                             M         +                                             I   43.391
  53.000      I                                             M              +                                        I   45.493
  54.000      I                                             M                   +                                   I   47.509
  55.000      I                                             M                       +                               I   49.231
  56.000      I                                             M                          +                            I   50.580
  57.000      I                                             M                            +                          I   51.625
  58.000      I                                             M                              +                        I   52.548
  59.000      I                                             M                                 +                     I   53.573
  60.000      I                                             M                                    +                  I   54.903
  61.000      I                                             M                                        +              I   56.621
  62.000      I                                             M                                            +          I   58.605
  63.000      I                                             M                                                 +     I   60.549
  64.000      I                                             M                                                    +  I   62.053
  65.000      I                                             M                                                      +I   62.768
  66.000      I                                             M                                                     + I   62.470
  67.000      I                                             M                                                  +    I   61.053
  68.000      I                                             M                                            +          I   58.498
  69.000      I                                             M                                    +                  I   54.856
  70.000      I                                             M                         +                             I   50.271
  71.000      I                                             M             +                                         I   44.992
  72.000      I                                             M+                                                      I   39.363
  73.000      I                                 +           M                                                       I   33.810
  74.000      I                     +                       M                                                       I   28.878
  75.000      I             +                               M                                                       I   25.219
  76.000      I        +                                    M                                                       I   23.401
  77.000      I         +                                   M                                                       I   23.440
  78.000      I           +                                 M                                                       I   24.701
  79.000      I               +                             M                                                       I   26.348
  80.000      I                  +                          M                                                       I   27.725
  81.000      I                    +                        M                                                       I   28.463
  82.000      I                    +                        M                                                       I   28.457
  83.000      I                   +                         M                                                       I   27.798
  84.000      I                +                            M                                                       I   26.694
  85.000      I             +                               M                                                       I   25.380
  86.000      I          +                                  M                                                       I   24.053
  87.000      I       +                                     M                                                       I   22.852
  88.000      I     +                                       M                                                       I   21.861
  89.000      I   +                                         M                                                       I   21.114
  90.000      I  +                                          M                                                       I   20.585
  91.000      I +                                           M                                                       I   20.228
  92.000      I +                                           M                                                       I   19.997
  93.000      I+                                            M                                                       I   19.856
  94.000      I+                                            M                                                       I   19.780
  95.000      I+                                            M                                                       I   19.773
  96.000      I+                                            M                                                       I   19.881
  97.000      I +                                           M                                                       I   20.172
  98.000      I  +                                          M                                                       I   20.663
  99.000      I    +                                        M                                                       I   21.322
  100.00      I     +                                       M                                                       I   22.121
  101.00      I        +                                    M                                                       I   23.043
  102.00      I          +                                  M                                                       I   24.074
                                                             <11-9>
1 103.00      I             +                               M                                                       I   25.170
  104.00      I               +                             M                                                       I   26.225
  105.00      I                 +                           M                                                       I   27.089
  106.00      I                  +                          M                                                       I   27.651
  107.00      I                   +                         M                                                       I   27.926
  108.00      I                   +                         M                                                       I   28.094
  109.00      I                    +                        M                                                       I   28.456
  110.00      I                      +                      M                                                       I   29.302
  111.00      I                          +                  M                                                       I   30.831
  112.00      I                               +             M                                                       I   33.091
  113.00      I                                      +      M                                                       I   35.983
  114.00      I                                             M                                                       I   39.318
  115.00      I                                             M        +                                              I   42.866
  116.00      I                                             M                +                                      I   46.383
  117.00      I                                             M                        +                              I   49.656
  118.00      I                                             M                              +                        I   52.532
  119.00      I                                             M                                    +                  I   54.905
  120.00      I                                             M                                        +              I   56.647
  121.00      I                                             M                                          +            I   57.612
  122.00      I                                             M                                          +            I   57.717
  123.00      I                                             M                                         +             I   56.979
  124.00      I                                             M                                     +                 I   55.517
  125.00      I                                             M                                 +                     I   53.529
  126.00      I                                             M                           +                           I   51.257
  127.00      I                                             M                      +                                I   48.939
  128.00      I                                             M                 +                                     I   46.789
  129.00      I                                             M             +                                         I   44.991
  130.00      I                                             M          +                                            I   43.671
  131.00      I                                             M        +                                              I   42.836
  132.00      I                                             M       +                                               I   42.376
  133.00      I                                             M      +                                                I   42.143
  134.00      I                                             M      +                                                I   41.978
  135.00      I                                             M     +                                                 I   41.754
  136.00      I                                             M    +                                                  I   41.392
  137.00      I                                             M   +                                                   I   40.857
  138.00      I                                             M +                                                     I   40.153
  139.00      I                                             M+                                                      I   39.348
  140.00      I                                            +M                                                       I   38.625
  141.00      I                                           + M                                                       I   38.249
  142.00      I                                           + M                                                       I   38.470
  143.00      I                                             M+                                                      I   39.409
  144.00      I                                             M   +                                                   I   41.011
  145.00      I                                             M        +                                              I   43.087
  146.00      I                                             M              +                                        I   45.406
  147.00      I                                             M                   +                                   I   47.742
  148.00      I                                             M                        +                              I   49.872
  149.00      I                                             M                            +                          I   51.603
  150.00      I                                             M                               +                       I   52.815
  151.00      I                                             M                                +                      I   53.486
  152.00      I                                             M                                 +                     I   53.642
  153.00      I                                             M                                +                      I   53.276
  154.00      I                                             M                              +                        I   52.344
  155.00      I                                             M                          +                            I   50.815
  156.00      I                                             M                     +                                 I   48.720
  157.00      I                                             M               +                                       I   46.180
  158.00      I                                             M         +                                             I   43.393
  159.00      I                                             M  +                                                    I   40.591
  160.00      I                                          +  M                                                       I   37.983
  161.00      I                                     +       M                                                       I   35.722
                                                             <11-10>
1 162.00      I                                 +           M                                                       I   33.921
  163.00      I                              +              M                                                       I   32.684
  164.00      I                             +               M                                                       I   32.106
  165.00      I                             +               M                                                       I   32.232
  166.00      I                               +             M                                                       I   32.994
  167.00      I                                  +          M                                                       I   34.207
  168.00      I                                     +       M                                                       I   35.618
  169.00      I                                        +    M                                                       I   36.974
  170.00      I                                           + M                                                       I   38.080
  171.00      I                                            +M                                                       I   38.791
  172.00      I                                             M                                                       I   39.002
  173.00      I                                            +M                                                       I   38.724
  174.00      I                                           + M                                                       I   38.055
  175.00      I                                        +    M                                                       I   37.128
  176.00      I                                      +      M                                                       I   36.047
  177.00      I                                   +         M                                                       I   34.859
  178.00      I                                +            M                                                       I   33.594
  179.00      I                             +               M                                                       I   32.295
  180.00      I                          +                  M                                                       I   31.060
  181.00      I                        +                    M                                                       I   29.998
  182.00      I                      +                      M                                                       I   29.205
  183.00      I                     +                       M                                                       I   28.723
  184.00      I                    +                        M                                                       I   28.547
  185.00      I                     +                       M                                                       I   28.668
  186.00      I                      +                      M                                                       I   29.109
  187.00      I                        +                    M                                                       I   29.919
  188.00      I                          +                  M                                                       I   31.119
  189.00      I                              +              M                                                       I   32.649
  190.00      I                                  +          M                                                       I   34.381
  191.00      I                                      +      M                                                       I   36.171
  192.00      I                                          +  M                                                       I   37.886
  193.00      I                                             M+                                                      I   39.408
  194.00      I                                             M  +                                                    I   40.621
  195.00      I                                             M    +                                                  I   41.430
  196.00      I                                             M     +                                                 I   41.815
  197.00      I                                             M     +                                                 I   41.858
  198.00      I                                             M     +                                                 I   41.681
  199.00      I                                             M    +                                                  I   41.339
  200.00      I                                             M   +                                                   I   40.812
  201.00      I                                             M +                                                     I   40.056
  202.00      I                                             M                                                       I   39.052
  203.00      I                                          +  M                                                       I   37.838
  204.00      I                                       +     M                                                       I   36.522
  205.00      I                                    +        M                                                       I   35.275
  206.00      I                                  +          M                                                       I   34.325
  207.00      I                                 +           M                                                       I   33.954
  208.00      I                                  +          M                                                       I   34.416
  209.00      I                                     +       M                                                       I   35.740
  210.00      I                                          +  M                                                       I   37.716
  211.00      I                                             M +                                                     I   40.017
  212.00      I                                             M      +                                                I   42.334
  213.00      I                                             M           +                                           I   44.438
  214.00      I                                             M               +                                       I   46.206
  215.00      I                                             M                   +                                   I   47.597
  216.00      I                                             M                     +                                 I   48.667
  217.00      I                                             M                       +                               I   49.589
  218.00      I                                             M                          +                            I   50.609
  219.00      I                                             M                             +                         I   51.878
  220.00      I                                             M                                +                      I   53.421
                                                             <11-11>
1 221.00      I                                             M                                    +                  I   55.145


























































                                                             <11-12>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 PLOT OF PHASE OF SMOOTHED DEMODULATED SERIES

 LOCATION OF ZERO IS GIVEN BY PLOT CHARACTER 0
             -6.2832   -5.0265   -3.7699   -2.5133   -1.2566     .0000    1.2566    2.5133    3.7699    5.0265    6.2832
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                       B                          0                      A                           I
  2.0000      I                       B                          0                      A                           I
  3.0000      I                        B                         0                       A                          I
  4.0000      I                        B                         0                       A                          I
  5.0000      I                        B                         0                       A                          I
  6.0000      I                         B                        0                        A                         I
  7.0000      I                         B                        0                        A                         I
  8.0000      I                         A                        0                        B                         I
  9.0000      I                         A                        0                        B                         I
  10.000      I                         A                        0                        B                         I
  11.000      I                          A                       0                         B                        I
  12.000      I                          A                       0                         B                        I
  13.000      I                          A                       0                         B                        I
  14.000      I                          A                       0                         B                        I
  15.000      I                          A                       0                         B                        I
  16.000      I                          A                       0                         B                        I
  17.000      I                          A                       0                         B                        I
  18.000      I                          A                       0                         B                        I
  19.000      I                          A                       0                         B                        I
  20.000      I                          A                       0                         B                        I
  21.000      I                          A                       0                         B                        I
  22.000      I                          A                       0                         B                        I
  23.000      I                          A                       0                         B                        I
  24.000      I                          A                       0                         B                        I
  25.000      I                         A                        0                        B                         I
  26.000      I                         A                        0                        B                         I
  27.000      I                         B                        0                        A                         I
  28.000      I                         B                        0                        A                         I
  29.000      I                         B                        0                        A                         I
  30.000      I                         B                        0                        A                         I
  31.000      I                         B                        0                        A                         I
  32.000      I                         B                        0                        A                         I
  33.000      I                         B                        0                        A                         I
  34.000      I                         B                        0                        A                         I
  35.000      I                        B                         0                       A                          I
  36.000      I                        B                         0                       A                          I
  37.000      I                        B                         0                       A                          I
  38.000      I                        B                         0                       A                          I
  39.000      I                        B                         0                       A                          I
  40.000      I                         B                        0                        A                         I
  41.000      I                         A                        0                        B                         I
  42.000      I                          A                       0                         B                        I
  43.000      I                           A                      0                          B                       I
  44.000      I                           A                      0                          B                       I
  45.000      I                            A                     0                           B                      I
  46.000      I                             A                    0                            B                     I
  47.000      I                             A                    0                            B                     I
  48.000      I                              A                   0                             B                    I
  49.000      I                               A                  0                              B                   I
  50.000      I                                A                 0                               B                  I
  51.000      I                                 A                0                                B                 I
  52.000      I                                   A              0                                  B               I
                                                             <11-13>
1 53.000      I                                    A             0                                   B              I
  54.000      I                                     A            0                                    B             I
  55.000      I                                      A           0                                     B            I
  56.000      I                                      A           0                                     B            I
  57.000      I                                       A          0                                      B           I
  58.000      I                                        A         0                                       B          I
  59.000      I                                          A       0                                         B        I
  60.000      I                                           A      0                                          B       I
  61.000      I                                            A     0                                           B      I
  62.000      I                                             A    0                                            B     I
  63.000      I                                             A    0                                            B     I
  64.000      I                                              A   0                                             B    I
  65.000      I                                              A   0                                             B    I
  66.000      I                                               A  0                                              B   I
  67.000      I                                               A  0                                              B   I
  68.000      I                                               A  0                                              B   I
  69.000      I                                               A  0                                              B   I
  70.000      I                                               A  0                                              B   I
  71.000      I                                              A   0                                             B    I
  72.000      I                                              A   0                                             B    I
  73.000      I                                             A    0                                            B     I
  74.000      I                                           A      0                                          B       I
  75.000      I                                         A        0                                        B         I
  76.000      I                                       A          0                                      B           I
  77.000      I                                    A             0                                   B              I
  78.000      I                                  A               0                                 B                I
  79.000      I                                 A                0                                B                 I
  80.000      I                               A                  0                              B                   I
  81.000      I                               A                  0                              B                   I
  82.000      I                              A                   0                             B                    I
  83.000      I                             A                    0                            B                     I
  84.000      I                             A                    0                            B                     I
  85.000      I                             A                    0                            B                     I
  86.000      I                            A                     0                           B                      I
  87.000      I                            A                     0                           B                      I
  88.000      I                           A                      0                          B                       I
  89.000      I                          A                       0                         B                        I
  90.000      I                          A                       0                         B                        I
  91.000      I                         B                        0                        A                         I
  92.000      I                        B                         0                       A                          I
  93.000      I                       B                          0                      A                           I
  94.000      I                      B                           0                     A                            I
  95.000      I                      B                           0                     A                            I
  96.000      I                     B                            0                    A                             I
  97.000      I                     B                            0                    A                             I
  98.000      I                    B                             0                   A                              I
  99.000      I                    B                             0                   A                              I
  100.00      I                   B                              0                  A                               I
  101.00      I                  B                               0                 A                                I
  102.00      I                  B                               0                 A                                I
  103.00      I                 B                                0                A                                 I
  104.00      I                B                                 0               A                                  I
  105.00      I                B                                 0               A                                  I
  106.00      I                B                                 0               A                                  I
  107.00      I                B                                 0               A                                  I
  108.00      I                B                                 0               A                                  I
  109.00      I                 B                                0                A                                 I
  110.00      I                  B                               0                 A                                I
  111.00      I                   B                              0                  A                               I
                                                             <11-14>
1 112.00      I                     B                            0                    A                             I
  113.00      I                      B                           0                     A                            I
  114.00      I                      B                           0                     A                            I
  115.00      I                       B                          0                      A                           I
  116.00      I                        B                         0                       A                          I
  117.00      I                        B                         0                       A                          I
  118.00      I                         B                        0                        A                         I
  119.00      I                         A                        0                        B                         I
  120.00      I                          A                       0                         B                        I
  121.00      I                          A                       0                         B                        I
  122.00      I                          A                       0                         B                        I
  123.00      I                          A                       0                         B                        I
  124.00      I                          A                       0                         B                        I
  125.00      I                          A                       0                         B                        I
  126.00      I                          A                       0                         B                        I
  127.00      I                         A                        0                        B                         I
  128.00      I                         A                        0                        B                         I
  129.00      I                         B                        0                        A                         I
  130.00      I                        B                         0                       A                          I
  131.00      I                        B                         0                       A                          I
  132.00      I                        B                         0                       A                          I
  133.00      I                       B                          0                      A                           I
  134.00      I                      B                           0                     A                            I
  135.00      I                      B                           0                     A                            I
  136.00      I                     B                            0                    A                             I
  137.00      I                     B                            0                    A                             I
  138.00      I                     B                            0                    A                             I
  139.00      I                    B                             0                   A                              I
  140.00      I                    B                             0                   A                              I
  141.00      I                     B                            0                    A                             I
  142.00      I                     B                            0                    A                             I
  143.00      I                     B                            0                    A                             I
  144.00      I                      B                           0                     A                            I
  145.00      I                      B                           0                     A                            I
  146.00      I                      B                           0                     A                            I
  147.00      I                       B                          0                      A                           I
  148.00      I                       B                          0                      A                           I
  149.00      I                       B                          0                      A                           I
  150.00      I                       B                          0                      A                           I
  151.00      I                       B                          0                      A                           I
  152.00      I                       B                          0                      A                           I
  153.00      I                      B                           0                     A                            I
  154.00      I                      B                           0                     A                            I
  155.00      I                      B                           0                     A                            I
  156.00      I                     B                            0                    A                             I
  157.00      I                     B                            0                    A                             I
  158.00      I                    B                             0                   A                              I
  159.00      I                    B                             0                   A                              I
  160.00      I                   B                              0                  A                               I
  161.00      I                   B                              0                  A                               I
  162.00      I                  B                               0                 A                                I
  163.00      I                  B                               0                 A                                I
  164.00      I                  B                               0                 A                                I
  165.00      I                  B                               0                 A                                I
  166.00      I                   B                              0                  A                               I
  167.00      I                   B                              0                  A                               I
  168.00      I                   B                              0                  A                               I
  169.00      I                   B                              0                  A                               I
  170.00      I                   B                              0                  A                               I
                                                             <11-15>
1 171.00      I                   B                              0                  A                               I
  172.00      I                   B                              0                  A                               I
  173.00      I                   B                              0                  A                               I
  174.00      I                   B                              0                  A                               I
  175.00      I                   B                              0                  A                               I
  176.00      I                   B                              0                  A                               I
  177.00      I                  B                               0                 A                                I
  178.00      I                  B                               0                 A                                I
  179.00      I                 B                                0                A                                 I
  180.00      I                B                                 0               A                                  I
  181.00      I                B                                 0               A                                  I
  182.00      I               B                                  0              A                                   I
  183.00      I              B                                   0             A                                    I
  184.00      I             B                                    0            A                                     I
  185.00      I             B                                    0            A                                     I
  186.00      I            B                                     0           A                                      I
  187.00      I            B                                     0           A                                      I
  188.00      I            B                                     0           A                                      I
  189.00      I           B                                      0          A                                       I
  190.00      I           B                                      0          A                                       I
  191.00      I           B                                      0          A                                       I
  192.00      I           B                                      0          A                                       I
  193.00      I           B                                      0          A                                       I
  194.00      I           B                                      0          A                                       I
  195.00      I           B                                      0          A                                       I
  196.00      I           B                                      0          A                                       I
  197.00      I           B                                      0          A                                       I
  198.00      I           B                                      0          A                                       I
  199.00      I            B                                     0           A                                      I
  200.00      I            B                                     0           A                                      I
  201.00      I            B                                     0           A                                      I
  202.00      I            B                                     0           A                                      I
  203.00      I             B                                    0            A                                     I
  204.00      I             B                                    0            A                                     I
  205.00      I              B                                   0             A                                    I
  206.00      I              B                                   0             A                                    I
  207.00      I               B                                  0              A                                   I
  208.00      I                B                                 0               A                                  I
  209.00      I                B                                 0               A                                  I
  210.00      I                 B                                0                A                                 I
  211.00      I                 B                                0                A                                 I
  212.00      I                  B                               0                 A                                I
  213.00      I                  B                               0                 A                                I
  214.00      I                  B                               0                 A                                I
  215.00      I                  B                               0                 A                                I
  216.00      I                  B                               0                 A                                I
  217.00      I                  B                               0                 A                                I
  218.00      I                   B                              0                  A                               I
  219.00      I                   B                              0                  A                               I
  220.00      I                   B                              0                  A                               I
  221.00      I                    B                             0                   A                              I








                                                             <11-16>
1G.  Acknowledgments

      The code for performing the complex demodulation was adapted from the
 subroutines given on pages 147 to 150 of Bloomfield [1976].























































                                    <11-17>
1-----                            CHAPTER 12                              -----

                      CORRELATION AND SPECTRUM ANALYSIS

 A.  Introduction

      STARPAC contains 50 subroutines for time series correlation  and spectrum
 estimation.   Both univariate and bivariate series can be analyzed.   Included
 are subroutines that compute  the correlation function using the  fast Fourier
 transform and that accept time series with missing observations.  The user may
 choose from  spectrum analysis subroutines implementing the  classical Fourier
 transformed  covariance  function techniques  presented in  Jenkins  and Watts
 [1968], the autoregressive or rational spectrum techniques described  by Jones
 [1971] or the  direct Fourier transform (periodogram) techniques  discussed in
 Bloomfield [1976].

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The  declaration  and CALL statements are given in section C and
 the subroutine arguments are defined in section D.  The  algorithms  used  and
 output  produced  by  these  subroutines  are  discussed in section E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

      STARPAC correlation  and spectrum  analysis subroutines are  divided into
 seven families.   For correlation analysis of univariate and bivariate  series
 there are two families of subroutines supporting
      1.  Autocorrelation Analysis and
      2.  Cross Correlation Analysis.
 For spectrum estimation there are four families of subroutines  for univariate
 series and one family for bivariate series supporting
      3.  Univariate Spectrum  Estimation Using  the Fourier  Transform  of the
          Autocorrelation  Function,
      4.  Univariate Spectrum Estimation Using Autoregressive Models,
      5.  Univariate  Spectrum Estimation  Using the  Direct Fourier Transform,
      6.  Univariate Series Utilities and
      7.  Bivariate  Spectrum  Estimation Using  the Fourier  Transform  of the
          Cross Correlation Function.

      In  general, each family  of subroutines  has one basic  subroutine which
 performs the  desired computations with  a minimum  of user input.   The other
 subroutines in  each family provide  greater flexibility  to the  user  at the
 price of more input.   The features of these subroutines are indicated  by the
 suffix letter(s) on the  subroutine (e.g., ACFM and BFSFS).   Not all features
 are available  for each family.   Features which  are common to more than  one
 family are described here.  Features which are unique to a specific family are
 described in the subsections below.

      * Suffix S indicates that the user is allowed to specify  various options
        which are  preset in the  simplest call  and that  certain  results are
        returned to the user  via the subroutine CALL statement.   In the  
	subsections  that  follow,  the  specific  details  of this feature are
        discussed individually for each family of subroutines.

      * Suffix M indicates that the series  contains  missing data.   A missing
        value code  must be  used  within the  series to  specify  time  points
        without an observed value.  There is no limit on the percentage of data

                                    <12-1>
1       that can be missing.   However, because the correlation matrix computed
        from a series with missing values is not necessarily positive definite,
        the partial autocorrelation function estimates and autoregressive order
        selection  statistics are  not  computed and  caution must  be  used in
        interpreting those results which are provided.  Analysis of time series
        with missing values is discussed in Jones [1971].

      * Suffix  F indicates  that  the  covariances  are  computed   using  the
        Singleton  [1969] fast  Fourier transform  (FFT).   When the  number of
        observations in the series is large this method of computation  is more
        efficient  than  the  direct  computation  normally  used  by  STARPAC.
        Subroutines with an F suffix  reduce the amount of workspace  needed by
        using the vector originally containing the data as workspace;  the data
        must be copied into  another vector prior to calling  these subroutines
        if  the  data are  to be  preserved.   These subroutines  automatically
        extend the length of the input series by appending enough zeros to meet
        the  requirements of this FFT code;  the length  of the vector  used to
        pass the data to these  subroutines must therefore equal or  exceed the
        extended series length, NFFT, as discussed in section D.

      * Suffix V indicates that the user inputs the covariances rather than the
        original series,  thus avoiding  a time-consuming recomputation  of the
        covariance function  if  it is  already available,  for  example,  from
        subroutines ACFS, ACFFS, ACFMS, CCFS, CCFFS or CCFMS.


 B.1  Correlation Analysis

 B.1.a  Univariate Series

      Autocorrelation Analysis.    STARPAC's  autocorrelation   function  (acf)
 subroutines compute  and plot the autocorrelation function  estimates; compute
 the large lag standard error of  the estimates; perform a chi-squared  test of
 the null  hypothesis that  the  series is  white noise;  compute  the  partial
 autocorrelation  function  coefficients  estimates; and,  using  the  modified
 Akaike  information  criterion   [Akaike,  1974],  select  the  order   of  an
 autoregressive process which models the series and estimate the  parameters of
 this autoregressive model.   The user should note that a purely autoregressive
 model may approximate the true  structure of  the model with  an unnecessarily
 large  number of  terms.   Such  an  autoregressive model  must be  used  with
 discretion since the true  structure might actually be more  complex including
 moving average components, harmonic terms or some mixture of deterministic and
 stochastic elements.  For some purposes, a purely autoregressive approximation
 may be useful.   In other cases, careful model identification can  lead to the
 discovery  of more detailed structure of  the data  or to a  more parsimonious
 model.

      The simplest of the Auto Correlation Function subroutines is  ACF,  which
 performs  the  basic analysis described in the preceding paragraph.  The other
 autocorrelation analysis subroutines provide the same basic  analysis  as  ACF
 while adding the features indicated above by suffixes S, M, F, MS and FS.

      For the ACF family of subroutines,  the suffix S feature allows  the user
 to indicate:
      1)  the maximum  lag value for  which the  correlation function is  to be
          computed; and
      2)  the amount of printed output.
 The  acf subroutines with  suffix S  also return  the  autocovariance function

                                    <12-2>
1estimates  and the coefficients  of the  selected autoregressive model  to the
 user via the subroutine CALL statement.

      The ACF family of  subroutines also  includes subroutine ACFD,  where the
 suffix D indicates that the  autocorrelation analysis will be performed  for a
 sequence  of differenced series.   The difference  factors are provided by the
 user.   If  the  number of  difference factors,  NFAC,  is  greater  than one,
 difference  factors beyond the first are  applied to the input series  Y(t) to
 yield a series Z(t) given by

                           NFAC
                Z(t) = [ PRODUCT (1 - B[IOD(J)])**ND(J) ] * Y(t)
                           J=2

 where the B[k] indicates the backward shift operator defined by

                               B[k]*Y(t) = Y(t-k)

 and  IOD and ND are defined in section D.  If the number of difference factors
 is equal to one, Z(t) = Y(t).  In either case, the autocorrelation analysis is
 performed first on the series Z and then  on  series  Z  with  the  difference
 factor (1-B[IOD(1)]) applied 1 to ND(1) times.  This produces ND(1) + 1 passes
 of the basic ACF analysis.


 B.1.b  Bivariate Series

      Cross   Correlation   Analysis.   STARPAC's  cross  correlation  analysis
 subroutines compute and plot the cross correlation function  coefficients  and
 provide  the  large  lag standard error of these estimates.  Subroutine CCF is
 the simplest of the Cross Correlation Function  subroutines.  The  other  five
 cross  correlation analysis subroutines provide the same basic analysis as CCF
 while adding the features indicated above by suffixes S, M, F, MS and FS.

      For the CCF family of  subroutines, suffix S indicates that  the analysis
 is provided for each pair of series of a multivariate time series.  The suffix
 S feature also allows the user to indicate:
      1)  the maximum  lag value for  which the  correlation function is  to be
          computed; and
      2)  the amount of printed output.
 In  addition, the  cross  covariance function  estimates are  returned  to the
 user.


 B.2  Spectrum Estimation

 B.2.a  Univariate Series

      Univariate  Spectrum  Estimation  Using  the  Fourier  Transform  of  the
 Autocorrelation  Function.  The  UFS  (Univariate Fourier Spectrum) subroutine
 family computes the estimated spectrum  from  the  Fourier  transform  of  the
 autocovariance  function (acvf) as discussed in Jenkins and Watts [1968].  The
 spectrum is smoothed using Parzen windows with the bandwidth controlled either
 by the user or a window-closing algorithm.  The principal output from each  of
 these subroutines consists of plots of the estimated spectrum.

      Subroutine  UFS has  the  simplest  CALL  statement  of  this  family  of
 subroutines.   The  printed  output  consists  of  four  spectrum  plots  with

                                    <12-3>
1successively narrower bandwidths.   Each spectrum is displayed in decibels (10
 times the base 10 logarithm of the power spectrum) scaled so that  the maximum
 value  plotted is  zero.   The  length  of  the  upper  and  lower  95-percent
 confidence intervals  and the bandwidth  for each  spectrum are  shown  on the
 plots.

      The other nine univariate Fourier spectrum estimation subroutines provide
 the  same basic analysis as UFS  while adding the features indicated  above by
 suffixes S, M, F, V, FS, MS, VS, MV and MVS.

      For the UFS family of subroutines,  the suffix S feature allows  the user
 to indicate:
      1)  the number  of different window  bandwidths to  be used  and  the lag
          window truncation point for each;
      2)  the  frequency range  and  the number  of frequencies  for  which the
          spectrum is to be computed; and
      3)  whether the plot is to be in decibels or on a logarithmic scale.
 In addition, the spectrum values are returned to the user.


      Univariate  Spectrum  Estimation  Using  Autoregressive  Models.  STARPAC
 Univariate  Autoregressive  Spectrum  estimation  subroutines   (UAS   family)
 approximate  an  input  series  with  an  autoregressive model and compute the
 corresponding theoretical spectrum for that model.  For comparative  purposes,
 the  plot  of  the  autoregressive  spectrum is superimposed against a Fourier
 spectrum plot.

      Subroutine UAS is the simplest of the autoregressive  spectrum estimation
 subroutines.  It uses the modified Akaike information criterion [Akaike, 1974]
 to select the order of the autoregressive model to be used. The 
 autoregressive coefficients are then computed from the autocovariance function
 using the Levinson-Durbin recursive method [Box and Jenkins, 1976] for solving
 the Yule-Walker equations. The lag window truncation point used for the Fourier
 spectrum is half the maximum  truncation point which would have  been selected
 by  subroutine  UFS  [see  section D,  argument LAGS].  The output consists of
 several  autoregressive  order  selection  statistics  and  a  plot   of   the
 autoregressive and Fourier spectra in decibels (10 times the base 10 logarithm
 of  the  power  spectrum)  scaled such that the maximum value plotted is zero.
 The bandwidth and the length of the 95-percent  confidence  interval  for  the
 Fourier  spectrum  are  shown  on  the  plot.  This  Fourier  spectrum and its
 confidence  intervals  should  be  used  in  interpreting  the  autoregressive
 spectrum  since  confidence  intervals  are  not  computed  by STARPAC for the
 autoregressive spectrum.  (The bandwidth is not relevant to the autoregressive
 spectrum.)

      The other five autoregressive spectrum subroutines provide the same basic
 analysis  as subroutine  UAS  while adding  the features  indicated  above  by
 suffixes S, F, V, FS and VS.

      For the  autoregressive  spectrum family  of subroutines,  the  suffix  S
 feature allows the user to indicate:
      1)  the order of the autoregressive model to be used for  the 
          autoregressive spectrum;
      2)  the lag window truncation point to be used for the Fourier spectrum;
      3)  the frequency range and  the number of frequencies within  this range
          at which the spectrum is to be computed; and
      4)  whether the plot is to be in decibels or on a logarithmic scale.
 In addition, the autoregressive and Fourier spectra are returned to  the user.

                                    <12-4>
1The user  should be cautious  about using  high order models  without checking
 order selection statistics since such models can produce spurious peaks in the
 spectrum.


      Univariate Spectrum Estimation Using the Direct  Fourier  Transform.  The
 STARPAC  direct  Fourier  transform  subroutines  (PGM  family)  implement the
 PeriodoGraM approach to time series analysis discussed in  Bloomfield  [1976].
 Subroutines  are  included for computing the raw periodogram and for computing
 and plotting the integrated periodogram (or cumulative spectrum).

      Subroutine PGM computes the periodogram of the series  as  described  for
 argument  PER  in  section  D, using  zeros  to extend the length of the input
 series.  Output consists of a plot of the computed periodogram in decibels (10
 times the base 10 logarithm of the periodogram estimates) scaled so  that  the
 maximum  value  plotted  is zero.  The input series must be either centered or
 tapered by the user before PGM is called.

      PGMS  provides the  same basic analysis  as PGM  but allows  the  user to
 indicate:
      1)  whether zeros or the series mean is used to extend the series;
      2)  the length of the extended series; and
      3)  the amount of printed output.
 In addition, the periodogram values and their frequencies are returned  to the
 user via the subroutine CALL statement.

      The integrated  periodogram subroutine  IPGM first  subtracts  the series
 mean  from the series, then extends  the input  series with zeros  and finally
 computes the normalized integrated periodogram.  Output consists of a one-page
 plot  of the integrated  periodogram, accompanied  by 95-percent  contours for
 testing the  null hypothesis  of white noise.   The integrated  periodogram is
 discussed in chapter 8 of Box and Jenkins [1976].

      The  other  three integrated  periodogram subroutines  add  the  features
 indicated by  suffixes S, P and PS.   The suffix  S option allows the user  to
 control the amount of  printed output;  the integrated periodogram  values and
 their  corresponding frequencies  are  also  returned  to  the  user  via  the
 subroutine CALL statement.  The suffix P option indicates that the user inputs
 the  periodogram   rather   than  the   original  series,   thus   avoiding  a
 time-consuming recomputation of  the periodogram  if it is  already available,
 for example, from subroutine PGMS.


      Utilities.   STARPAC   includes   utility   subroutines   for   centering
 (subtracting the mean) and tapering the observed series, periodogram smoothing
 and  for computing the Fourier coefficients of the series.  These routines are
 particularly useful when using direct Fourier techniques such as PGM and IPGM.

      Subroutine CENTER subtracts the  series mean from the series  and returns
 the centered series.  There is no printed output.

      Subroutine TAPER centers the input series and applies  the  
 split-cosine-bell taper described for argument YT in section D. The user 
 specifies the total proportion of the series to be tapered. The centered 
 tapered series is returned to the user and can be used as input to subroutine
 PGM or PGMS.  There is no printed output.



                                    <12-5>
1     Subroutine FFTLEN computes for an observed series length, N,  the minimum
 extended  series length,  NFFT,  which  will  meet  the  requirements  of  the
 Singleton FFT code.   The value  of the extended series length  is returned to
 the user.  There is no printed output.

      Subroutine MDFLT smooths the input periodogram by applying a  sequence of
 modified Daniell filters as discussed in chapter 7 of Bloomfield [1976].   The
 filtered  series  is  returned  to the  user.   There  is  no printed  output.
 Subroutine MDFLT takes advantage of  the symmetry of the periodogram  to avoid
 losing values from  the ends of the series.   It should  therefore not be used
 for input series that are not symmetric about their end values.  Other digital
 filtering  routines, such as those described  in chapter 10, may also  be used
 for periodogram smoothing but end effect losses will be incurred.

      Subroutine FFTR computes the  Fourier coefficients of an input  series of
 single precision observations.  There is no printed output.


 B.2.b  Bivariate Series

      Bivariate Spectrum Estimation Using the Fourier Transform  of  the  Cross
 Correlation  Function.  The BFS (Bivariate Fourier Spectrum) subroutine family
 computes the estimated spectrum from the Fourier transform of  the  covariance
 function  as  discussed in Jenkins and Watts [1968].  The spectrum is smoothed
 using Parzen windows with the bandwidth controlled either by  the  user  or  a
 window-closing algorithm.  The principal output from each of these subroutines
 consists  of  plots of the squared coherency and phase components of the cross
 spectrum.  The phase function plots reduce discontinuities  using  the  method
 suggested  by  Bloomfield [1976].  As shown in the example in section F,  this
 method displays both the principle phase value, which is defined to lie in the
 range -pi to pi,  and the principle phase value plus or minus 2*pi,  where the
 sign is chosen such that the second value lies in the range -2*pi to 2*pi.

      Subroutine  BFS  provides the basic analysis with a brief CALL statement.
 The printed output consists of four spectrum plot pairs (a  squared  coherency
 plot  and  a  phase  plot) with successively narrower bandwidths chosen by the
 window-closing algorithm.  The upper and lower 95-percent confidence intervals
 and the 95-percent significance levels are shown on the coherency plots.

      The other nine bivariate Fourier spectrum estimation  subroutines provide
 the basic analysis described for BFS while adding the features indicated above
 by suffixes S, M, F, V, FS, MS, VS, MV and MVS.

      For the BFS family of subroutines,  the suffix S feature allows  the user
 to indicate:
      1)  the number  of different window  bandwidths to  be used  and  the lag
          window truncation point for each;
      2)  the  frequency range  and  the number  of frequencies  for  which the
          spectrum is to be computed; and
      3)  whether the plot is to be in decibels or on a logarithmic scale.
 In addition, the squared coherency and phase values are returned to the user.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.

                                    <12-6>
1

                   Subroutines for Autocorrelation Analysis

 ACF:     Compute and print a two-part auto and  partial  correlation  analysis
          of a series,  select the  order of an  autoregressive  process  which
          models the series, and estimate the parameters of this model

          <real> Y(n)
          :
          :
          CALL ACF (Y, N)

                                      ===

 ACFS:    Compute  and optionally print a two-part auto and partial correlation
          analysis of a series,  select the order  of an autoregressive process
          which models the series,  and estimate the  parameters of this model;
          use  user-supplied  control values;  return  autocovariance function,
          and order and parameter estimates of selected autoregressive model

          <real> Y(n), ACOV(lagmax+1), PHI(lagmax)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL ACFS (Y, N,
         +           LAGMAX, LACOV, ACOV, IAR, PHI, NPRT, LDSTAK)

                                      ===

 ACFM:    Compute  and  print  a two-part auto and partial correlation analysis
          of a series with missing observations

          <real> Y(n)
          :
          :
          CALL ACFM (Y, YMISS, N)

                                      ===

 ACFMS:   Compute  and optionally print a two-part auto and partial correlation
          analysis  of  a  series  with missing observations; use user-supplied
          control values; return autocovariance function

          INTEGER NLPPA(lagmax+1)
          <real> Y(n), ACOV(lagmax+1)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL ACFMS (Y, YMISS, N,
         +            LAGMAX, LACOV, ACOV, AMISS, NLPPA, NPRT, LDSTAK)

                                      ===




                                    <12-7>
1ACFF:    Compute  and  print  a two-part auto and partial correlation analysis
          of  a  series,  select  the  order of an autoregressive process which
          models the series, and estimate  the  parameters  of this model;  use
          FFT for computations

          <real> YFFT(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL ACFF (YFFT, N, LYFFT, LDSTAK)

                                      ===

 ACFFS:   Compute  and optionally print a two-part auto and partial correlation
          analysis of a series,  select the order  of an autoregressive process
          which models the series,  and estimate the  parameters of this model;
          use  FFT for computations;  use user-supplied control values;  return
          autocovariance  function,   and  order  and  parameter  estimates  of
          selected autoregressive model

          <real> YFFT(nfft), ACOV(lagmax+1), PHI(lagmax)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL ACFFS (YFFT, N, LYFFT, LDSTAK,
         +            LAGMAX, LACOV, ACOV, IAR, PHI, NPRT)

                                      ===

 ACFD:    Compute  and  print  a two-part auto and partial correlation analysis
          of  a  sequence  of  differenced  series,  select  the  order  of  an
          autoregressive process which models each  series,  and  estimate  the
          parameters of these models

          INTEGER ND(nfac), IOD(nfac)
          <real> Y(n)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL ACFD (Y, N, LAGMAX, NFAC, ND, IOD, LDSTAK)

                                      ===














                                    <12-8>
1                  Subroutines for Cross Correlation Analysis

 CCF:     Compute  and print a two-part cross correlation analysis of a pair of
          series

          <real> Y1(n), Y2(n)
          :
          :
          CALL CCF (Y1, Y2, N)

                                      ===

 CCFS:    Compute  and  optionally print a two-part cross correlation  analysis
          of a  multivariate series using user-supplied control  values; return
          cross covariance function

          <real> YM(n,m), CCOV(lagmax+1,m,m)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL CCFS (YM, N, M, IYM,
         +           LAGMAX, CCOV, ICCOV, JCCOV, NPRT, LDSTAK)

                                      ===

 CCFM:    Compute  and print a two-part cross correlation analysis of a pair of
          series with missing observations

          <real> Y1(n), Y2(n)
          :
          :
          CALL CCFM (Y1, YMISS1, Y2, YMISS2, N)

                                      ===

 CCFMS:   Compute  and  optionally print a two-part cross correlation  analysis
          of  a  multivariate  series  with missing  observations  using  
	  user-supplied control values; return cross covariance function

          INTEGER NLPPC(lagmax+1,m,m)
          <real> YM(n,m), YMMISS(m), CCOV(lagmax+1,m,m)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL CCFMS (YM, YMMISS, N, M, IYM,
         +            LAGMAX, CCOV, CMISS, ICCOV, JCCOV,
         +            NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)

                                      ===








                                    <12-9>
1CCFF:    Compute  and print a two-part cross correlation analysis of a pair of
          series; use FFT for computations

          <real> YFFT1(nfft), YFFT2(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL CCFF(YFFT1, YFFT2, N, LYFFT, LDSTAK)

                                      ===

 CCFFS:   Compute  and  optionally print a two-part cross correlation  analysis
          of a multivariate series using user-supplied control values;  use FFT
          for computations; return cross covariance function

          <real> YMFFT(nfft,m), CCOV(lagmax+1,m,m)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL CCFFS (YMFFT, N, M, IYMFFT,
         +            LAGMAX, CCOV, ICCOV, JCCOV, NPRT, LDSTAK)

                                      ===


                 Subroutines for Univariate Spectrum Estimation
          Using the Fourier Transform of the Autocorrelation Function

 UFS:     Compute and print a univariate Fourier spectrum analysis of a series

          <real> Y(n)
          :
          :
          CALL UFS (Y, N)

                                      ===

 UFSS:    Compute and optionally  print a univariate Fourier  spectrum analysis
          of  a  series  using  user-supplied control  values;  return  Fourier
          spectrum and corresponding frequencies

          INTEGER LAGS(nw)
          <real> Y(n), SPCF(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSS (Y, N,
         +           NW, LAGS, NF, FMIN, FMAX, NPRT,
         +           SPCF, ISPCF, FREQ, LDSTAK)

                                      ===





                                    <12-10>
1UFSF:    Compute and print a univariate Fourier spectrum analysis of a series;
          use FFT for computations

          <real> YFFT(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSF (YFFT, N, LYFFT, LDSTAK)

                                      ===

 UFSFS:   Compute and optionally  print a univariate Fourier  spectrum analysis
          of a series using user-supplied control values; use FFT  for 
	  computations; return Fourier spectrum and corresponding frequencies

          INTEGER LAGS(nw)
          <real> YFFT(nfft), SPCF(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSFS (YFFT, N, LYFFT, LDSTAK,
         +           NW, LAGS, NF, FMIN, FMAX, NPRT,
         +           SPCF, ISPCF, FREQ)

                                      ===

 UFSM:    Compute and print a univariate Fourier  spectrum analysis of a series
          with missing observations

          <real> Y(n)
          :
          :
          CALL UFSM (Y, YMISS, N)

                                      ===

 UFSMS:   Compute and optionally  print a univariate Fourier  spectrum analysis
          of a  series with  missing observations  using  user-supplied control
          values; return Fourier spectrum and corresponding frequencies

          INTEGER LAGS(nw)
          <real> Y(n), SPCF(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSMS (Y, YMISS, N,
         +            NW, LAGS, NF, FMIN, FMAX, NPRT,
         +            SPCF, ISPCF, FREQ, LDSTAK)

                                      ===






                                    <12-11>
1UFSV:    Compute and print a univariate Fourier spectrum analysis of a series;
          input covariances rather than original series

          <real> ACOV(lagmax+1)
          :
          :
          CALL UFSV (ACOV, LAGMAX, N)

                                      ===

 UFSVS:   Compute and optionally  print a univariate Fourier  spectrum analysis
          of a  series using  user-supplied control  values;  input covariances
          rather than original series; return Fourier spectrum  and 
	  corresponding frequencies

          INTEGER LAGS(nw)
          <real> ACOV(lagmax+1), SPCF(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSVS (ACOV, LAGMAX, N,
         +            NW, LAGS, NF, FMIN, FMAX, NPRT,
         +            SPCF, ISPCF, FREQ, LDSTAK)

                                      ===

 UFSMV:   Compute and print a univariate Fourier spectrum  analysis of a series
          with  missing observations;  input covariances  rather  than original
          series

          INTEGER NLPPA(lagmax+1)
          <real> ACOV(lagmax+1)
          :
          :
          CALL UFSMV(ACOV, NLPPA, LAGMAX, N)

                                      ===

 UFSMVS:  Compute and optionally  print a univariate Fourier  spectrum analysis
          of a  series with  missing observations  using  user-supplied control
          values; input covariances rather than original series; return Fourier
          spectrum and corresponding frequencies

          INTEGER NLPPA(lagmax+1)
          <real> ACOV(lagmax+1), SPCF(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UFSMVS (ACOV, NLPPA, LAGMAX, N,
         +             NW, LAGS, NF, FMIN, FMAX, NPRT,
         +             SPCF, ISPCF, FREQ, LDSTAK)

                                      ===




                                    <12-12>
1  Subroutines for Univariate Spectrum Estimation Using Autoregressive Models

 UAS:     Compute and print a  univariate autoregressive spectrum analysis of a
          series

          <real> Y(n)
          :
          :
          CALL UAS (Y, N)

                                      ===

 UASS:    Compute and  optionally  print a univariate  autoregressive  spectrum
          analysis  of  a series  using user-supplied  control  values;  return
          autoregressive and Fourier spectrum and corresponding frequencies

          <real> Y(n), PHI(lagmax), SPCF(nf), SPCA(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UASS (Y, N,
         +           IAR, PHI, LAGMAX, LAG, NF, FMIN, FMAX, NPRT,
         +           SPCA, SPCF, FREQ, LDSTAK)

                                      ===

 UASF:    Compute  and print a univariate autoregressive spectrum analysis of a
          series; use FFT for computations

          <real> YFFT(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UASF(YFFT, N, LYFFT, LDSTAK)

                                      ===

 UASFS:   Compute  and  optionally  print a  univariate autoregressive spectrum
          analysis of a series using user-supplied control values; use  FFT for
          computations;  return   autoregressive  and   Fourier   spectrum  and
          corresponding  frequencies

          <real> YFFT(nfft),  PHI(lagmax), SPCA(nf), SPCF(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UASFS (YFFT, N, LYFFT, LDSTAK,
         +            IAR, PHI, LAGMAX, LAG, NF, FMIN, FMAX, NPRT,
         +            SPCA, SPCF, FREQ)

                                      ===





                                    <12-13>
1UASV:    Compute  and print a univariate autoregressive spectrum analysis of a
          series; input covariances rather than original series

          <real> ACOV (lagmax+1)
          :
          :
          CALL UASV (ACOV, LAGMAX, N)

                                      ===

 UASVS:   Compute and  optionally  print a univariate  autoregressive  spectrum
          analysis  of  a  series  using user-supplied  control  values;  input
          covariances  rather than  original series; return  autoregressive and
          Fourier spectrum and corresponding frequencies

          <real> ACOV(lagmax+1), Y(n), PHI(lagmax)
          <real> SPCA(nf), SPCF(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL UASVS (ACOV, LAGMAX, Y, N,
         +            IAR, PHI, LAG, NF, FMIN, FMAX, NPRT,
         +            SPCA, SPCF, FREQ, LDSTAK)

                                      ===


            Subroutines for Univariate Spectrum Estimation Using the
                            Direct Fourier Transform

 PGM:     Compute and print a  periodogram  analysis  of a series;  use FFT for
          computations

          <real> YFFT(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL PGM (YFFT, N, LYFFT, LDSTAK)

                                      ===

 PGMS:    Compute and optionally print a periodogram analysis of a series;  use
          FFT  for computations; return periodogram and  corresponding 
	  frequencies

          <real> YFFT(nfft), PER(nf), FREQ(nf)
          :
          :
          CALL PGMS (YFFT, N, NFFT, LYFFT,
         +           IEXTND, NF, PER, LPER, FREQ, LFREQ, NPRT)

                                      ===





                                    <12-14>
1IPGM:    Compute  and print  an integrated  periodogram analysis of  a series;
          use FFT for computations

          <real> YFFT(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL IPGM (YFFT, N, LYFFT, LDSTAK)

                                      ===

 IPGMS:   Compute  and optionally print an integrated periodogram analysis of a
          series; use  FFT for computations; return integrated  periodogram and
          corresponding frequencies

          <real> YFFT(nfft), PERI(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL IPGMS (YFFT, N, LYFFT, LDSTAK,
         +            NF, PERI, LPERI, FREQ, LFREQ, NPRT)

                                      ===

 IPGMP:   Compute  and print  an integrated  periodogram analysis of  a series;
          input periodogram rather than original series

          <real> PER(nf), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL IPGMP (PER, FREQ, NF, N, LDSTAK)

                                      ===

 IPGMPS:  Compute  and optionally print an integrated periodogram analysis of a
          series; input  periodogram rather than original series;  return 
	  integrated periodogram and corresponding frequencies

          <real> PER(nf), FREQ(nf), PERI(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL IPGMPS (PER, FREQ, NF, N, LDSTAK, PERI, NPRT)

                                      ===









                                    <12-15>
1                             Utility Subroutines

 CENTER:  Subtract  the  series mean from each observation of  a series; return
          the centered series (no printed output)

          <real> Y(n), YC(n)
          :
          :
          CALL CENTER (Y, N, YC)

                                      ===

 TAPER:   Center  a series  about its mean and apply a split-cosine-bell taper;
          return the tapered series (no printed output)

          <real> Y(n), YT(n)
          :
          :
          CALL TAPER (Y, N, TAPERP, YT)

                                      ===

 FFTLEN:  Compute the  minimum extended  series length for  using the Singleton
          FFT; return the extended series length (no printed output)

          CALL FFTLEN (N, NDIV, NFFT)

                                      ===

 MDFLT:   Smooth a periodogram  by applying  a  sequence  of  modified  Daniell
          filters; return the smoothed periodogram (no printed output)

          INTEGER KMD(nk)
          <real> PER(nf), PERF(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL MDFLT (PER, NF, NK, KMD, PERF, LDSTAK)

                                      ===

 FFTR:    Compute the  Fourier  coefficients  of  an  input  series  of  <real>
          observations; return the Fourier coefficients (no printed output)

          <real> YFFT(n), AB(nfft)
          :
          :
          CALL FFTR (YFFT, N, NFFT, IEXTND, NF, AB, LAB)

                                      ===








                                    <12-16>
1           Subroutines for Bivariate Spectrum Estimation Using the
              Fourier Transform of the Cross Correlation Function

 BFS:     Compute and print a bivariate Fourier  spectrum analysis of a pair of
          series

          <real> Y1(n), Y2(n)
          :
          :
          CALL BFS (Y1, Y2, N)

                                      ===

 BFSS:    Compute and optionally  print a bivariate  Fourier spectrum  analysis
          of a  pair  of  series  using user-supplied  control  values;  return
          squared coherency and phase components of the cross spectrum  and the
          corresponding frequencies

          INTEGER LAGS(nw)
          <real> Y1(n), Y2(n), CSPC2(nf,nw), PHAS(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSS (Y1, Y2, N,
         +           NW, LAGS, NF, FMIN, FMAX, NPRT,
         +           CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)

                                      ===

 BFSF:    Compute and  print a bivariate Fourier spectrum analysis of a pair of
          series; use FFT for computations

          <real> YFFT1(nfft), YFFT2(nfft)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSF (YFFT1, YFFT2, N, LYFFT, LDSTAK)

                                      ===


















                                    <12-17>
1BFSFS:   Compute and optionally  print a bivariate  Fourier spectrum  analysis
          of a pair of series  using user-supplied control values; use  FFT for
          computations; return  squared coherency  and phase components  of the
          cross spectrum and the corresponding frequencies

          INTEGER LAGS(nw)
          <real> YFFT1(nfft), YFFT2(nfft), CSPC2(nf,nw), PHAS(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSFS (YFFT1, YFFT2, N, LYFFT, LDSTAK,
         +            NW, LAGS, NF, FMIN, FMAX, NPRT,
         +            CSPC2, ICSPC2, PHAS, IPHAS, FREQ)

                                      ===

 BFSM:    Compute and print a  bivariate Fourier spectrum analysis of a pair of
          series with missing observations

          <real> Y1(n), Y2(n)
          :
          :
          CALL BFSM (Y1, YMISS1, Y2, YMISS2, N)

                                      ===

 BFSMS:   Compute  and optionally  print a bivariate Fourier spectrum  analysis
          of a pair  of series  with missing  observations  using user-supplied
          control values; return squared coherency and phase components  of the
          cross spectrum and the corresponding frequencies

          INTEGER LAGS(nw)
          <real> Y1(n), Y2(n), CSPC2(nf,nw), PHAS(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSMS (Y1, YMISS1, Y2, YMISS2, N,
         +            NW, LAGS, NF, FMIN, FMAX, NPRT,
         +            CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)

                                      ===

 BFSV:    Compute  and print a bivariate Fourier spectrum analysis of a pair of
          series; input covariances rather than original series

          <real> CCOV(lagmax+1,m,m)
          :
          :
          CALL BFSV (CCOV, INDEX1, INDEX2, N, LAGMAX, ICCOV, JCCOV)

                                      ===






                                    <12-18>
1BFSVS:   Compute and optionally  print a bivariate  Fourier spectrum  analysis
          of a  pair  of  series  using  user-supplied  control  values;  input
          covariances rather than original series; return squared coherency and
          phase components of the cross spectrum and the corresponding 
	  frequencies

          INTEGER LAGS(nw)
          <real> CCOV(lagmax+1,m,m), CSPC2(nf,nw), PHAS(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSVS (CCOV, INDEX1, INDEX2, N, ICCOV, JCCOV,
         +            NW, LAG, NF, FMIN, FMAX, NPRT,
         +            CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)

                                      ===

 BFSMV:   Compute  and print a bivariate Fourier spectrum analysis of a pair of
          series  with  missing  observations; input  covariances  rather  than
          original series

          INTEGER NLPPC(lagmax+1,m,m)
          <real> CCOV(lagmax+1,m,m)
          :
          :
          CALL BFSMV (CCOV, NLPPC, INDEX1, INDEX2, N, LAGMAX,
         +            ICCOV, JCCOV, INLPPC, JNLPPC)

                                      ===

 BFSMVS:  Compute and optionally  print a bivariate  Fourier spectrum  analysis
          of a pair  of series  with missing  observations  using user-supplied
          control values; input covariances rather than original series; return
          squared coherency and phase components of the cross spectrum  and the
          corresponding frequencies

          INTEGER NLPPC(lagmax+1,m,m), LAGS(nw)
          <real> CCOV(n,m,m), CSPC2(nf,nw), PHAS(nf,nw), FREQ(nf)
          DOUBLE PRECISION DSTAK(ldstak)
          COMMON /CSTAK/ DSTAK
          :
          :
          CALL BFSMVS (CCOV, NLPPC, INDEX1, INDEX2, N,
         +             ICCOV, JCCOV, INLPPC, JNLPPC,
         +             NW, LAGS, NF, FMIN, FMAX, NPRT,
         +             CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)

                                      ===










                                    <12-19>
1D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the  argument is input to the subroutine  and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 AB      <-- The vector of dimension at least 2*NF  that contains  the  NF real
             and the NF imaginary components of the Fourier coefficients of the
             input  series.  The real and imaginary components  of the  Fourier
             coefficients  are  returned  in  the array AB stored such that the
             real and imaginary parts of the complex Fourier coefficients are

             AB(2*I-1) + i*AB(2*I) for I = 1, ..., NF,

             where

             i is the complex value sqrt(-1); and

             NF is the number of harmonic frequencies, NF = NFFT/2.

             The vector YFFT used to input the observed series can also be used
             in place of the vector AB to conserve storage space.

 ACOV    --- The vector  of  dimension at  least LAGMAX + 1  that  contains the
             LAGMAX+1 autocovariance function (acvf) estimates for  lags i = 0,
             ..., LAGMAX.   ACOV is  returned from the autocorrelation analysis
             subroutines; it is input to the spectrum analysis subroutines.

             The acvf estimate for lag i, c(i), is defined by

                                    N-i       _          _
                            (N-i) * SUM (Y(t)-Y)*(Y(t+i)-Y)*u(t)*u(t+i)|
                                    t=1
             c(i) = c(-i) = ---------------------------------------------
                                          N-i
                                        N SUM u(t)*u(t+i)
                                          t=1

             for lags i = 0 to LAGMAX, where
             _
             Y is the average of all observed values; and

             u(t) is an indicator variable, defined by

               u(t) = 1 if Y(t) is observed (Y(t) <> YMISS), and

               u(t) = 0 if Y(t) is missing (Y(t) = YMISS).

             The acvf are stored in ACOV such that ACOV(I) = c(I-1) for I  = 1,
             ..., LAGMAX+1.   When there are no missing observations the  above


                                    <12-20>
1            formula  for the  acvf  reduces  to  that of  the  usual  positive
             definite acvf estimator.

 AMISS   <-- The  missing value  code  used within  ACOV to  indicate  that the
             autocovariance function  at  a given  lag could  not  be  computed
             because of missing data.

 CCOV    --- The three-dimensional  array  of dimension at least LAGMAX+1  by M
             by M that contains the cross covariance function (ccvf) estimates.
             CCOV is returned by the cross correlation analysis subroutines; it
             is input to the spectrum analysis subroutines.

             The ccvf estimate, cjk(i), for lag i between the series  stored in
             the jth and kth column of YM is defined by

             cjk(i) = ckj(-i)

                              N-i          __              __
                      (N-i) * SUM (YM(t,j)-YMj)*(YM(t+i,k)-YMk)*uj(t)*uk(t+i)
                              t=1
                    = -------------------------------------------------------
                                           N-i
                                       N * SUM uj(t)*uk(t+i)
                                           t=1

             for lags i = 0, ..., LAGMAX and for  series j = 1, ..., M and  k =
             1, ..., M, where
             __
             YMj is the average of all  the observed values for the  jth column
             of YM;

             uj(t) is an indicator variable, defined by

               uj(t) = 1 if YM(t,j) is observed (YM(t,j) <> YMMISS(J))

               uj(t) = 0 if YM(t,j) is missing (YM(t,j) = YMMISS(J)).

             When there are no  missing observations the above formula  for the
             ccvf reduces to that of the usual positive definite ccvf estimator
             and  when j = k the  above formula  is that of  the autocovariance
             function for the jth series.

             The ccvf are stored in CCOV such that

             CCOV(I,J,K) = cjk(i-1) = ckj(-i+1) for I = 1, ..., LAGMAX+1.

             The appropriate formulas  for the  ccvf estimates computed  by CCF
             and  CCFM can be  obtained by  letting M = 2  and  by substituting
             Y1(t) for YM(t,1) and Y2(t) for YM(t,2) in the above.

 CMISS   <-- The missing  value code used  in CCOV  to indicate that  the cross
             covariance function at a  given lag could not be  computed because
             of missing data.

 CSPC2   <-- The matrix  of dimension  at  least NF  by NW  that  contains  the
             squared  coherency  component of  the cross  spectra  between  two
             series, designated j and k, respectively.  CSPC2(I,L) contains the


                                    <12-21>
1            smoothed squared  coherency value  for the Ith  frequency computed
             using the Lth  lag point specified in LAGS.   The returned  values
             are expressed  as  shown  below  regardless  of  the  plot  option
             selected.

             The estimated  smoothed squared  coherency component of  the cross
             spectrum is

                          COSPEC(I,L)**2 + QSPEC(I,L)**2
             CSPC2(I,L) = ------------------------------
                              SPCFj(I,L)*SPCFk(I,L)

             for I = 1, ..., NF and L = 1, ..., NW, where

             COSPEC(I,L) is the smoothed co-spectra for the Lth lag point,

               COSPEC(I,L) = EV[0] +

                                 LAGS(L)-1
                             2 *    SUM    EV(P)*WL(P)*cos[2*pi*P*(FMIN+DELi)];
                                    P=1

             QSPEC(I,L) is the smoothed quadrature  spectra  for  the  Lth  lag
               point,

                                 LAGS(L)-1
               QSPEC(I,L) =  2 *    SUM    OD(P)*WL(P)*sin[2*pi*P*(FMIN+DELi)];
                                    P=1

             WL(P)  is the  Parzen  lag window  function for  the  Lth  window,
               defined by

               WL(P) = 1                                        P  =  0

               WL(P) = 1 - 6*[|P|/LAGS(L)]**2 +

                       6*[|P|/LAGS(L)]***3                1 <= |P| <= LAGS(L)/2

               WL(P) = 2*[1-(|P|/LAGS(L))]**3     LAGS(L)/2  < |P| <= LAGS(L)

               WL(P) = 0                            LAGS(L)  < |P| ;

             DELi is the frequency increment,

               DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1);

             EV(P) and OD(P) are functions of CCOV defined by

               EV(P) = CCOV(P+1,j,k) + CCOV(P+1,k,j) for P = 0, ..., LAGS(L)

               OD(P) = CCOV(P+1,j,k) - CCOV(P+1,k,j) for P = 0, ..., LAGS(L);

             SPCFj(I,L)  and SPCFk(I,L)  are the  univariate  Fourier  spectrum
               estimates for series  j and k, respectively, computed  using the
               Lth lag window truncation point [see argument SPCF].




                                    <12-22>
1            Note that the modifications necessary for series with missing data
             are included in the computation of CCOV.

 DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 FMAX    --> The maximum frequency, in cycles per sample interval, at which the
             spectrum  is  to  be  computed  (0.0  <=  FMIN < FMAX F 0.5).  The
             default value is 0.5.  If FMAX is outside the range FMIN to 0.5 or
             is not an argument in the CALL  statement  the  default  value  is
             used.

 FMIN    --> The minimum frequency, in cycles per sample interval, at which the
             spectrum  is  to  be  computed  (0.0  <= FMIN < FMAX <= 0.5).  The
             default value is 0.0.  If FMIN is outside the range 0.0 to FMAX or
             is not an argument in the CALL  statement  the  default  value  is
             used.

 FREQ    --- The vector of dimension at least NF that contains the NF frequency
             values   at  which  the  spectrum is  computed.  FREQ is  an input
             argument to subroutines IPGMP and IPGMPS.   The values of FREQ are
             returned by  all  other subroutines  including it  in  their  CALL
             statements.

 IAR     --- The order of  the autoregressive  model chosen to  approximate the
             series and  the order of  the autoregressive  model to be  used in
             computing  the autoregressive spectrum.   IAR is  returned by  the
             autocorrelation  subroutines; it  is input  to  the autoregressive
             spectrum analysis subroutines.

             For the autoregressive  spectrum subroutines:

               The order  IAR  of  the  autoregressive model  must  not  exceed
               LAGMAX.  In no case may IAR exceed N/2.

               If IAR>0, the user  must supply  the coefficients for  the order
                         IAR autoregressive  model in  the vector PHI.   (These
                         coefficients are available, for example,  from 
			 subroutine ACFS  and  ACFFS or  perhaps from the ARIMA model
                         fitting procedure discussed in chapter 13.)

               If IAR<0, the coefficients  for the  order  |IAR| autoregressive
                         model will be computed by the STARPAC subroutine using
                         Durbin's recursive method.

               If IAR=0, the order and model coefficients will be  chosen using
                         the model selection criteria discussed below in {E.2.b
                         and THE INPUT VALUE OF IAR WILL BE OVERWRITTEN  BY THE
                         SELECTED  VALUE.   If the  IAR = 0 option  is used the
                         user MUST specify the order with a variable, i.e.,

                           IAR = 0
                           CALL ASPS(..., IAR, ...)

                         NOT


                                    <12-23>
1
                           CALL ASPS(..., 0, ...).

                         The  latter  will  cause  the  value  of  zero  to  be
                         redefined on  many computers  including  the CYBER 840
                         and 855.

 ICCOV   --> The exact value of the first dimension of CCOV as specified in the
             calling program.

 ICSPC2  --> The exact value of the  first dimension  of CSPC2 as  specified in
             the calling program.

 IERR    ... An  error  flag  returned  in  COMMON  /ERRCHK/   [see  chapter 1,
             section D.5].  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             IERR = 0 indicates that no errors were detected.

             IERR = 1 indicates that improper input was detected.

 IEXTND  --> The indicator  variable  used to  designate whether  zero  or  the
             series mean is to be  used to  extend the series.   If IEXTND = 0,
             zero will be used.  If IEXTND <> 0, the series mean will be used.

 INDEX1  --> The  index of  the  first  series  used  in  computing  the  cross
             covariance function, corresponding to the subscript j used  in the
             definition of CCOV.

 INDEX2  --> The  index of  the  second  series  used in  computing  the  cross
             covariance function, corresponding to the subscript k used  in the
             definition of CCOV.

 INLPPC  --> The exact value of the  first dimension  of NLPPC as  specified in
             the calling program.

 IOD     --> The  vector of  dimension  at least  NFAC that  contains  the NFAC
             values designating the order of each difference factor.

 IPHAS   --> The exact value of the first dimension of PHAS as specified in the
             calling program.

 ISPCF   --> The exact value of the first dimension of SPCF as specified in the
             calling program.

 IYM     --> The exact  value of  the  first  dimension  of the  matrix  YM  as
             specified in the calling program.

 IYMFFT  --> The exact  value of the  first dimension  of the  matrix  YMFFT as
             specified in the calling program.

 JCCOV   --> The exact value of the  second dimension  of CCOV as  specified in
             the calling program.

 JNLPPC  --> The exact value of the  second dimension of NLPPC as  specified in
             the calling program.


                                    <12-24>
1
 KMD     --> The vector of dimension at least NK that contains the  NK modified
             Daniell filter lengths.  All values in KMD must be even.

 LAB     --> The length of the vector AB.  LAB must equal or exceed NFFT.

 LACOV   --> The length of  the vectors  ACOV and NLPPA.   LACOV must  equal or
             exceed LAGMAX + 1.

 LAG     <-> The  lag window  truncation  point to  be used  for  computing the
             Fourier  spectrum.    The  default   value  is  half  the  maximum
             truncation  point  selected by  the algorithm  described  for  the
             default  values of LAGS.   If LAG  is not an argument of  the CALL
             statement or  if LAG is  outside the  range [1,  N-1]  the default
             value  will be used.  If the  user supplied value for LAG  is less
             than or equal to zero, THE INPUT VALUE WILL BE OVERWRITTEN  BY THE
             SELECTED VALUE.  The user MUST therefore  specify the  lag  window
             truncation point  with a variable  in this case, i.e.,

               LAG = 0
               CALL ASPS(..., LAG, ...)

             NOT

               CALL ASPS(..., 0, ...).

             The  latter will cause the value  of zero to be redefined  on many
             computers including the CYBER 840 and 855.

 LAGMAX  --> The maximum lag  value for which the correlation  coefficients are
             computed.   The default value of LAGMAX is selected  by STARPAC as
             follows.

               LAGMAX =  100      if   301 <= N.

               LAGMAX =  N/3      if    96 <= N <= 300.

               LAGMAX =   32      if    33 <= N <=  95.

               LAGMAX =  N-1      if          N <=  32.

             If LAGMAX is less than or equal  to zero or if neither  LAGMAX nor
             LAGS is an argument in  the CALL  statement the  default  value is
             used.   When LAGS is an argument in the  CALL statement, LAGMAX is
             the largest value in LAGS.

 LAGS    --> The vector  of dimension  at  least NW  that contains  the  NW lag
             window truncation points,  1 <= LAGS(i) <= N-1 for i = 1, ...  NW.
             By default,  four lag  window  truncation  points  are  used.  The
             smallest  lag  window truncation point T1 is selected by examining
             the covariance function.  It is chosen to be 3/16  times  the  lag
             value  beyond  which the covariance function remains less than the
             95-percent confidence limits for white noise,  although in no case
             is  T1  greater  than  LAGMAX/8.  The  values of the remaining lag
             window truncation points are then specified by T2  =  2*T1,  T3  =
             4*T1   and   T4  =  8*T1,   resulting  in  progressively  narrower
             bandwidths.   The  procedure  of  using   progressively   narrower


                                    <12-25>
1            bandwidths  to  compute  the  spectrum  of  a given time series is
             called window closing.  It  is  discussed  in  Jenkins  and  Watts
             [1968].  If LAGS is not an argument in the CALL statement the four
             default values are used.

 LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must  equal  or  exceed the appropriate value given below.  In the
             following specifications of the value of  LDSTAK,  if  the  single
             precision version of STARPAC is being used P = 0.5,  otherwise P =
             1.0 [see chapter 1, section B].  Also, IO = 0 if NPRT = 0 and IO =
             1 if NPRT <> 0.

             For autocorrelation subroutine

               ACFS:    LDSTAK >= 12 + [5*LAGMAX + 1]*P

               ACFMS:   LDSTAK >= IO*(13 + [6*LAGMAX + 1]*P)

               ACFF:    LDSTAK >= 6 + NFFT*P

               ACFFS:   LDSTAK >= 12 + [4*LAGMAX + NFFT + 1]*P

               ACFD:    LDSTAK >= 16 + [7*LAGMAX + N + 2]*P

             For cross correlation subroutine

               CCFS:    LDSTAK >= 12 + [2*M + IO*(4*LAGMAX + 2)]*P

               CCFMS:   LDSTAK >= (26+M)/2 + [2*M + IO*(4*LAGMAX + 2)]*P

               CCFF:    LDSTAK >= 6 + NFFT*P

               CCFFS:   LDSTAK >= 13 + [NFFT + 2*M + IO*(4*LAGMAX + 2)]*P

             For univariate Fourier spectrum subroutine

               UFSS:    LDSTAK >= [26 + IO*(NF+5)]/2 +
                                  (2*LAGMAX + 2 + IO*(2*NF+10)|*P

               UFSF:    LDSTAK >= 7 + NFFT*P

               UFSFS:   LDSTAK >= [26 + IO*(NF+5)]/2 +
                                  (LAGMAX + 1 + NFFT + IO*(2*NF+10)|*P

               UFSMS:   LDSTAK >= [29 + LAGMAX + 1 + IO*(NF+5)]/2 +
                                  (2*LAGSMX + 2 + IO*(2*NF+10)|*P

               UFSVS:   LDSTAK >= [23 + IO*(NF+5)]/2 +
                                  (LAGMAX + 1 + IO*(2*NF+10)|*P

               UFSMVS:  LDSTAK >= [23 + IO*(NF+5)]/2 +
                                  (LAGMAX + 1 + IO*(2*NF+10)|*P







                                    <12-26>
1            For autoregressive spectrum subroutine

               UASS:    LDSTAK >= [32 + IO*(2*NF+5)]/2 +
                                  [2*LAGMAX + 2 + IA*(3*LAGMAX+1) +
                                                      IO*(4*NF+10)]*P

               UASF:    LDSTAK >= 7 + NFFT*P

               UASFS:   LDSTAK >= [32 + IO*(2*NF+5)]/2 +
                                  [2*LAGMAX + 2 + NFFT + IA*(3*LAGMAX+1) +
                                                             IO*(4*NF+10)]*P

               UASVS:   LDSTAK >= [29 + IO*(2*NF+5)]/2 +
                                  [LAGMAX + 1 + IA*(3*LAGMAX+1) +
                                                    IO*(4*NF+10)]*P

               where IA = 0 if IAR <> 0 and IA = 1 if IAR = 0.

             For direct Fourier spectrum subroutine

               PGM:     LDSTAK >= 9 + NFFT*P

               IPGM:    LDSTAK >= [123 + NFFT/2]/2 + [2*NFFT + 206]*P

               IPGMS:   LDSTAK >= IO*([123  + NFFT/2]/2 + [2*NFFT + 206]*P)

               IPGMP:   LDSTAK >= [126 + NF]/2 + [3*NF + 206]*P

               IPGMPS:  LDSTAK >= IO*([123 + NF]/2 + [2*NF + 206]*P)

             For utility subroutine

               MDFLT:   LDSTAK >= 7 + NF*P

             For bivariate Fourier spectrum subroutine

               BFSS:    LDSTAK >= [38+IO*4*NF]/2 + [7*LAGMAX+7+IO*8*NF]*P

               BFSF:    LDSTAK >= 7 + NFFT*P

               BFSFS:   LDSTAK >= [38 + IO*4*NF]/2 +
                                  [6*LAGMAX + 6 + NFFT + IO*8*NF]*P

               BFSMS:   LDSTAK >= [45 + 4*LAGMAX + IO*4*NF]/2 +
                                  [7*LAGMAX + 7 + 2*NF + IO*8*NF]*P

               BFSVS:   LDSTAK >= [35 + IO*4*NF]/2 +
                                  [3*LAGMAX + 3 + 2*NF + IO*8*NF]*P

               BFSMVS:  LDSTAK >= [35 + IO*4*NF]/2 +
                                  [3*LAGMAX + 3 + 2*NF + IO*8*NF]*P

 LFREQ   --> The length of the vector FREQ.  LFREQ must equal or exceed NF.

 LPER    --> The length of the vector PER.  LPER must equal or exceed NF.

 LPERI   --> The length of the vector PERI.  LPERI must equal or exceed NF.


                                    <12-27>
1
 LYFFT   --> The maximum  length  allowed for  the extended  series  created in
             YFFT.   LYFFT  must  equal or  exceed the  actual  extended series
             length, NFFT.   If too  small a  value of LYFFT  is used  an error
             message giving the correct value is generated.

 M       --> The number of columns of data in YM.

 N       --> The  number of time  points in  each series.  For  the correlation
             analysis subroutines the minimum  number of time points is  3; for
             the  spectrum  analysis subroutines  the minimum  number  of  time
             points is 17.

 ND      --> The  vector of dimension  at least  NFAC that contains  the values
             designating the number of  times each  difference factor is  to be
             applied.

 NDIV    --> A required  constant  used by  FFTLEN to  determine  the  extended
             series length [see NFFT].  NDIV must be two for a simple  FFT.  It
             must be four when the covariance function is being computed.

 NF      <-> The number of frequencies at which the spectrum is computed.

             For the  Fourier and  autoregressive spectrum subroutines  and the
             utility subroutines:

               NF is an input argument.  The default value of NF is 101.  If NF
               is not an argument of the CALL statement the default  value will
               be used.

             For the direct Fourier spectrum subroutines:

               NF  is returned by the subroutine.   The returned  value is NF =
               NFFT/2.

 NFAC    --> The number of difference factors.

 NFFT    --- The minimum extended series length that meets the  requirements of
             the Singleton FFT code.   NFFT is  returned by FFTLEN; it must  be
             input  to PGMS and FFTR; it  is computed internally by all  of the
             other subroutines using the FFT.

             The extended length of the series is

             NFFT = N + LAGMAX + K

             where N is the number of observations in the series.  K  >=  2  is
             chosen  so that NFFT is as small as possible,  NFFT-2 is divisible
             by 4 and has no prime factors greater than 23,  and the product of
             the  square-free  prime factors of NFFT-2 does not exceed 209.  In
             general,  NFFT will be less than N+LAGMAX+100 for the correlation,
             Fourier spectrum and autoregressive spectrum subroutines.  STARPAC
             subroutine FFTLEN, when called using the sequence

             CALL FFTLEN(N+LAGMAX, NDIV, NFFT)

             with NDIV = 2 when a simple FFT is to be computed and


                                    <12-28>
1
                  NDIV = 4 when the covariance function is to be computed,

             returns  the  value  of  NFFT  that  is  actually  used  by  these
             subroutines.

 NK      --> The number of modified Daniell filters to be applied.

 NLPPA   <-- The vector  of  dimension  at  least LAGMAX+1  that  contains  the
             LAGMAX+1 values  designating the  number of  lagged  product pairs
             used to compute the autocovariance function at each lag,

                        N-i
             NLPPA(I) = SUM u(t)*u(t+i)  for i = 0, ..., LAGMAX and I = i + 1,
                        t=1

             where u(t) = 1 if Y(t) <> YMISS, and u(t) = 0 if Y(t) = YMISS.

 NLPPC   <-- The three-dimensional array of dimension at least LAGMAX+1 by M by
             M that contains the number of lagged product pairs used to compute
             the  cross covariance function  for each  pair of  series  at each
             lag,

                            N-i
             NLPPC(I,J,K) = SUM uj(t)*uk(t+i) for i = 0, ..., LAGMAX
                            t=0               and I = i + 1,

             where the indices J and K correspond directly to the subscripts  j
             and k,  respectively,  for j = 1,  ...,  M and k = 1, ..., M;  and
             uj(t) = 1 if YM(t,j) <> YMMISS(j) and  uj(t)  =  0  if  YM(t,j)  =
             YMMISS(j).

 NPRT    --> The argument controlling printed output.   In each case, when NPRT
             is not an argument  in the  subroutine CALL statement  the default
             value is used.

             For the correlation analysis subroutines:

               If NPRT = 0, the printed output is suppressed.

               If NPRT <> 0, the printed output is provided.

               The default value is NPRT <> 0.

             For the Fourier and autoregressive spectrum analysis subroutines:

               If NPRT <= -1, the spectra are plotted  in decibels on  a linear
                              scale adjusted so that the peak is at zero.

               If NPRT  =  0, the printed output is suppressed.

               If NPRT >=  1, the spectra are plotted on a log-linear scale.

               The default value is NPRT = -1.





                                    <12-29>
1            For the periodogram subroutines:

               If NPRT <= -2, the printed output consists of a page plot of the
                              periodogram  versus  frequency  on  a  log-linear
                              scale.

               If NPRT  = -1, the printed output consists of a page plot of the
                              periodogram in decibels on a linear scale.

               If NPRT  =  0, the printed output is suppressed.

               If NPRT  =  1, the printed output consists of a vertical plot of
                              the periodogram in  decibels on  a  linear scale.
                              The vertical axis shows the frequency in the left
                              margin.

               If NPRT >=  2, the printed output consists of a vertical plot of
                              the periodogram on a log scale. The vertical axis
                              shows the frequency in the left margin.

               The default value is NPRT = -1.

             For the integrated periodogram subroutines:

               If NPRT = 0, the printed output is suppressed.

               If NPRT <> 0, the printed output is provided.

               The default value is NPRT <> 0.

 NW      --> The number of different window bandwidths to be used.  The default
             value is four for the Fourier spectrum subroutines and one for the
             autoregressive spectrum subroutines.   When NW  is not an argument
             of the subroutine CALL statement the default value is used.

 PER     --- The  vector  of  dimension  at  least  NF  that  contains  the  NF
             periodogram values.   PER is returned by subroutine PGMS; it  must
             be input to subroutines IPGMP, IPGMPS and MDFLT.   The vector YFFT
             used  to input  the observed series  to PGMS  can also be  used in
             place of the vector PER to conserve storage space.  (The values in
             YFFT  will be overwritten even when  YFFT is not used in  place of
             PER.)

             The series periodogram is computed at the harmonic frequencies

             fk = k/(NFFT-2) for k = 0, ..., (NFFT-2)/2

             by a direct Fourier transformation of the series using

             PER(I) = A(I)**2 + B(I)**2 for I = 1, ..., NF,

             where

             NF is the number of harmonic frequencies at which  the periodogram
               is computed, NF = NFFT/2;

             A(I) is the real component of the Fourier coefficient,


                                    <12-30>
1
                         2     NFFT
               A(I) = ------ * SUM YCT(t)*cos[2*pi*t*(I-1)/(NFFT-2)];
                      NFFT-2   t=1

             B(I) is the imaginary component of the Fourier coefficient,

                         2     NFFT
               B(I) = ------ * SUM YCT(t)*sin[2*pi*t*(I-1)/(NFFT-2)]; and
                      NFFT-2   t=1

             YCT(t) is the value of  the centered (or tapered) input  series at
               time t [see arguments YC and YT].

 PERF    <-- The vector of dimension at least NF that contains the NF values of
             the  periodogram  smoothed  by  applying a  sequence  of  modified
             Daniell filters to the raw periodogram.   The sequence of filtered
             series is

                         KMD(p)
             PERF(p,I) =  SUM  hP(j)*PERF(p-1,tI)  for I = 1, ..., NF
                          P=1                      and p = 1, ..., NK,

             where

             KMD(p) is the number of terms  in the  pth filter (KMD(p)  must be
               even);

             hP(j) is the jth filter coefficient of the pth filter, defined by

               hP(j) = 1/(2*KMD(p))  for j = 1 and KMD(p)

               hP(j) = 1/KMD(p)      for j = 2, ..., KMD(p)-1;

             PERF(p,.) is the periodogram smoothed using the first p filters
               (PERF(0,.) is the input periodogram and PERF(NK,.) is the series
               actually returned to the user);

             tI is index of PERF(p-1,.) defined by

               tI =  2 - [I-1+j+KMD(P)/2]                  1 <= I <=KMD(P)/2

               tI =       I-1+j-KMD(P)/2          KMD(P)/2+1 <= I <= N-KMD(P)/2

               tI = 2N - [I+1+j+KMD(P)/2]       N-KMD(P)/2+1 <= I <= N

 PERI    <-- The vector of dimension  at least NF that contains  the integrated
             periodogram values.   The vector  YFFT used to input the  observed
             series  to IPGMS  may be  used  in place  of the  vector  PERI  to
             conserve storage space.

             The integrated periodogram is defined by

                          1      I
             PERI(I) = ------ * SUM PER(K) for I = 1, ..., NF
                       N*YVAR   K=1



                                    <12-31>
1            where

                     1     N        _            _   1    N
             YVAR = --- * SUM (Y(t)-Y)**2  with  Y = - * SUM Y(t).
                    N-1   t=1                        N   t=1

 PHAS    <-- The matrix  of dimension  at  least NF  by NW  that  contains  the
             smoothed phase component of the cross spectra between  two series,
             designated j  and k, respectively.   PHAS(I,L) contains  the phase
             value  for the Ith  frequency computed  using the  Lth  lag window
             truncation point specified in LAGS.

             The estimated phase component of the cross spectrum is

                                QSPEC(I,L)**2
             PHAS(I,L) = arctan --------------
                                COSPEC(I,L)**2

             for I = 1, ..., NF and L = 1, ..., NW, where

             COSPEC(I,L) is the smoothed co-spectra for the Lth lag point,

               COSPEC(I,L) = EV[0] +

                                 LAGS(L)-1
                             2 *    SUM    EV(P)*WL(P)*cos[2*pi*P*(FMIN+DELi)];
                                    P=1

             QSPEC(I,L) is the smoothed quadrature  spectra  for  the  Lth  lag
               point,

                                 LAGS(L)-1
               QSPEC(I,L) =  2 *    SUM    OD(P)*WL(P)*sin[2*pi*P*(FMIN+DELi)];
                                    P=1

             WL(P)  is the  Parzen  lag window  function for  the  Lth  window,
               defined by

               WL(P) = 1                                        P  =  0

               WL(P) = 1 - 6*[|P|/LAGS(L)]**2 +

                       6*[|P|/LAGS(L)]***3                1 <= |P| <= LAGS(L)/2

               WL(P) = 2*[1-(|P|/LAGS(L))]**3     LAGS(L)/2  < |P| <= LAGS(L)

               WL(P) = 0                            LAGS(L)  < |P| ;

             DELi is the frequency increment,

               DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1);

             EV(P) and OD(P) are functions of CCOV defined by

               EV(P) = CCOV(P+1,j,k) + CCOV(P+1,k,j) for P = 0, ..., LAGS(L)

               OD(P) = CCOV(P+1,j,k) - CCOV(P+1,k,j) for P = 0, ..., LAGS(L);


                                    <12-32>
1
 PHI     --- The  vector of dimension  at least  LAGMAX that  contains  the IAR
             coefficients  of  the  order  IAR autoregressive  model.   PHI  is
             returned  by  the  autocorrelation  subroutines; it  is  input  or
             returned by  the autoregressive  spectrum  estimation  subroutines
             depending on the value of IAR.

 SPCA    <-- The vector of dimension at least NF that contains the NF values of
             the autoregressive  spectrum.   The autoregressive  spectrum 
	     estimates are defined by

                                  IAR
             SPCA(I) = SI2 * [1 - SUM PHI(P)*exp(i*2*pi*P*(FMIN+DELi))]**(-2)
                                  P=1

             for I = 1, ..., NF, where

             SI2 is the residual or one step prediction variance of  the  order
               IAR autoregressive model,

                        N                 IAR
               SI2 = ------- * [ACOV(1) - SUM PHI(P)*ACOV(P+1)],
                     N-IAR-1              P=1

             i is the complex value sqrt(-1); and

             DELi is the frequency increment,

               DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1).

 SPCF    <-- The  array of  dimension at least  NF by  NW that contains  the NF
             values of  the Fourier spectrum  for each  of the NW  lag windows.
             For the autoregressive  spectrum analysis subroutines, NW =  1 and
             SPCF may  be dimensioned  as  a  vector  of length  at  least  NF.
             SPCF(I,L)  contains  the  spectrum  value for  the  Ith  frequency
             computed  using the Lth  lag point  (see arguments LAG  and LAGS).
             The  returned  spectrum  values  are  expressed  as   shown  below
             regardless of the plot option selected.

             The estimated Fourier spectrum values are

             SPCF(I,L) = ACOV(1) +

                             LAGMAX
                         2 *  SUM  ACOV(P+1)*WL(P)*cos[2*pi*P*(FMIN+DELi)]
                              P=1

             for I = 1, ..., NF and L = 1, ..., NW, where











                                    <12-33>
1            WL(P)  is the  Parzen  lag window  function for  the  Lth  window,
               defined by

               WL(P) = 1                                        P  =  0

               WL(P) = 1 - 6*[|P|/LAGS(L)]**2 +

                       6*[|P|/LAGS(L)]***3                1 <= |P| <= LAGS(L)/2

               WL(P) = 2*[1-(|P|/LAGS(L))]**3     LAGS(L)/2  < |P| <= LAGS(L)

               WL(P) = 0                            LAGS(L)  < |P| ;

             DELi is the frequency increment,

               DELi = 2*(I-1)*(FMAX-FMIN)/(NF-1).

             Note that the modifications necessary for series with missing data
             are included in the computation of ACOV.

 TAPERP  --> The total  proportion of  the  input data  to be  tapered  using a
             split-cosine-bell taper [see argument YT].   If TAPER <= 0.0,  the
             tapered  series is  identical  to  the  centered series,  YC.   If
             TAPERP >= 1.0, a 100-percent taper is applied to the series.

 Y       --> The vector of dimension at  least N  that contains the  N 
             observations of a time series.

 YC      <-- The vector of dimension at least N that contains the centered time
             series.  The centered series is
                            _
             YC(t) = Y(t) - Y for t = 1, ..., N
                   _
             where Y is the mean of the input series

             _          N
             Y = (1/N) SUM Y(t).
                       t=1

 YFFT    <-> The  vector  of   dimension  at   least  NFFT  containing   the  N
             observations of a time series to  be analyzed using the FFT.   The
             length  of the  vector YFFT must  be greater  than N to  allow for
             extending the series.   Note that  the input series is overwritten
             by the FFT computations.

 YFFT1   <-> The  vector  of  dimension  at   least  NFFT   containing   the  N
             observations of  the first of  a pair  of series  selected  from a
             multivariate time series to be analyzed using the FFT.  The length
             of the vector YFFT must be  greater than N to allow  for extending
             the series.   Note that the input series is overwritten by the FFT
             computations.

 YFFT2   <-> The  vector  of  dimension  at   least  NFFT   containing   the  N
             observations of  the second of  a pair  of series selected  from a
             multivariate time series to be analyzed using the FFT.  The length
             of the vector YFFT must be  greater than N to allow  for extending



                                    <12-34>
1            the series.   Note that the input series is overwritten by the FFT
             computations.

 YM      --> The matrix  of dimension at least N  by M each of whose  M columns
             contains the N observations of a multivariate time series.

 YMFFT   <-> The matrix of dimension at least NFFT by M each of whose M columns
             contains the N observations  of a  multivariate time series  to be
             analyzed  using the FFT.   The first  dimension of the array YMFFT
             must  be greater than N to  allow for extending the series.   Note
             that the input series are overwritten by the FFT computations.

 YMISS   --> The missing value code used within the input series Y  to indicate
             that an observation is missing.

 YMISS1  --> The missing value code used within the input series Y1 to indicate
             that an observation is missing.

 YMISS2  --> The missing value code used within the input series Y2 to indicate
             that an observation is missing.

 YMMISS  --> The vector  of dimension at  least M  that contains the  M missing
             value codes (one for each column) used within each of the M series
             contained in YM to indicate that an observation is missing.

 YT      <-- The vector of dimension at least N that contains the  tapered time
             series.  The tapered series is

             YT(t) = W(t)*YC(t) for t = 1, ..., N,

             where

             W(t) is the split-cosine-bell taper weight used at time t,

               W(t) = 0.5*(1-cos[pi*(t-0.5)/m])                   1 <= t <= m

               W(t) =   1                                     m + 1 <= t <= N-m

               W(t) = 0.5*(1-cos[pi*(N-t+0.5)/m])           N-m + 1 <= t <= N

             with m computed such that 2m/N = TAPERP is the proportion  of data
             to be tapered.  The importance of tapering is discussed in chapter
             5 of Bloomfield [1976].

 Y1      --> The vector of dimension at  least N  that contains the  N 
             observations of the first series of a multivariate time series pair.

 Y2      --> The vector of dimension at  least N  that contains the  N 
             observations of the second series of a multivariate time series pair.


 E.  Computational Methods

 E.1  Algorithms

 E.1.a  Correlation Analysis



                                    <12-35>
1     The code for computing the autocovariance and cross  covariance functions
 using a  fast Fourier  transform  are based  on subroutines  written  by Jones
 [1971].   The subroutines  which compute the fast Fourier  transform and which
 compute the  Fourier transform of  real data  are those  written  by Singleton
 [1969] and  the subroutine  to  compute the  cosine transform  was  written by
 Richard H. Jones.   A discussion of the use of the  fast Fourier transform for
 computing  the correlation  function can  be  found on  pages 165  to  167  of
 Bloomfield [1976].   Zeros are automatically appended to the end of the series
 to meet the requirements of  computing the correlation function using  an FFT.
 Correlation function subroutines  not using  the FFT  compute  the correlation
 function directly using the formulas given in section D for arguments ACOV and
 CCOV.

      The autoregressive  model coefficients  given are  the  Yule-Walker 
 estimates, computed as described in Appendix A3.2 of Box and Jenkins [1976].  This
 is NOT  a recommended estimation  procedure but  rather a  means  of providing
 starting values for a  maximum likelihood  or least squares  estimation
 procedure.  Chapter 13 provides subroutines for estimating the least squares values
 of the parameters of an autoregressive model using these starting values.


 E.1.b  Spectrum Analysis

      The transformation used to compute the  univariate  Fourier  spectrum  is
 performed  using  the algorithm shown on page 311 of Jenkins and Watts [1968],
 modified to correspond to the definition of the spectrum  given  for  argument
 SPCF  in  section D.  (The computed spectrum is half that shown in Jenkins and
 Watts.) The bivariate or cross Fourier spectrum is discussed in chapters 8 and
 9 of Jenkins and Watts.  Smoothed squared coherency and phase  estimation  are
 discussed in detail on pages 377-380 of Jenkins and Watts.  The algorithm used
 is  described on pages 418-420.  The confidence interval and significance test
 for coherency  is  discussed  in  Bloomfield  [1976]  on  pages  224-228.  The
 modifications  necessary  for  series  with  missing  data are included in the
 computation of the covariance function described in section D under  arguments
 ACOV  and  CCOV.  Covariance  function  values  computed  using a fast Fourier
 transform use the FFT code written by Singleton [1968];  a discussion of  this
 technique  for  computing the covariance function can be found on pages 165 to
 167 of Bloomfield [1976].  Subroutines using the  FFT  automatically  subtract
 the  series  mean  and  append  zeros  to  the  end  of the series to meet the
 requirements of computing the covariance function using Singleton's code.

      The code for computing the autoregressive order selection  statistics and
 autoregressive spectrum estimates  is based  on subroutines  written  by Jones
 [1971].   The coefficients PHI(k), k = 1, ..., IAR of the autoregressive model
 are  computed  from  the  autocovariance function  using  the  Levinson-Durbin
 recursive method for  solving the Yule-Walker equations discussed  in Appendix
 A3.2  of  Box  and Jenkins  [1976].   The  order  of the  autoregressive model
 selected by STARPAC is that  having the  lowest Akaike final  prediction error
 [Akaike, 1974].

      Subroutines PGM, IPGM, and IPGMS automatically append zeros to the end of
 the series to meet the length requirements of the Singleton FFT code [1969].

      Subroutine PGMS and FFTR do not automatically center or extend  the input
 series.   Centering and extending the series by appending  either zeros or the
 series mean must be done by the user.  The resulting computations are wrong if
 the extension does not comply with the code requirements.


                                    <12-36>
1
      The subroutines for the split-cosine-bell taper and the  modified Daniell
 filter operation were adapted from subroutines TAPER and MODDAN given on pages
 116 and 178 of Bloomfield [1976].


 E.2  Computed Results and Printed Output

 E.2.a  Correlation Analysis

      The Autocorrelation Subroutines.  The argument  controlling  the  printed
 output,  NPRT, is discussed in section D.  The output from the autocorrelation
 analysis subroutines includes summary statistics describing the input  series;
 the  autocorrelation  function  estimates and their large lag standard errors;
 the chi-squared test statistic for white noise and its significance level; the
 partial autocorrelation function estimates;  autoregressive model order 
 selection statistics; and parameter estimates for the selected model.  In addition,
 plots are provided of the autocorrelation and partial autocorrelation function
 estimates.

      The estimated  coefficients of the autocorrelation function  (acf), r(i),
 are computed with the formula

              r(i) = ACOV(i+1)/ACOV(1) for lag i = 1, ..., LAGMAX.

      The estimated standard errors (SE) of the acf are the large  lag standard
 errors discussed  on pages 34 to 36  of Box and Jenkins [1976].   The standard
 error listed at  lag i  is computed assuming  that the  acf is  zero  for lags
 greater than i.  The formula is given by

                                                      q
                            (N-i) * sqrt[(N-i) + 2 * SUM  (N-i-v)*r(v)**2]
                                                     v=1
     SE[r(i)] = SE[r(-i)] = ----------------------------------------------
                                            N*NLPPA(i+1)

 for i = 1, ..., LAGMAX, where q is the  minimum of n - i and i - 1.   The user
 must  select the "correct" large lag  standard error based on whether  the acf
 estimates  are compatible  with the assumption  that the  acf is zero  for the
 selected  lag  value and  beyond.   This  technique  is intended  primarily to
 provide a  means of  selecting the  largest  lag value  for which  the  acf is
 significantly different from zero rather than to provide actual standard error
 estimates for each acf.   The large lag standard error at the  selected lag is
 valid for all acf estimates at lags greater than or equal to the selected lag;
 the  standard error estimates for lag  values less  than the selected  lag are
 meaningless.

     As discussed  on pages  64-66  of Box  and Jenkins  [1976],  the estimated
 partial  autocorrelation function (pacf)  coefficients, p(i)  for i =  1, ...,
 LAGMAX, are estimates of the  ith coefficient in an autoregressive  process of
 order  i.   Use of  the pacf,  in conjunction  with  the acf,  to identify  an
 autoregressive moving average process is discussed on pages 174-193 of Box and
 Jenkins.   It should  be noted, however, that the  method used  to compute the
 pacf  [discussed in  Appendix A3.2 of  Box and  Jenkins] is very  sensitive to
 rounding errors.   The pacf estimates may not be reliable if the values of the
 parameters are close to the nonstationarity boundaries.



                                    <12-37>
1     The order IAR of the autoregressive model selected for the series is that
 having the lowest Akaike final prediction error [Akaike, 1974].   The modified
 Akaike information criteria (AIC) is then computed as

                 AIC(J) = N*Ln(FPEJ) for J = 1, ..., LAGMAX,

 where

 LAGMAX  is the maximum lag value for which the acvf has been estimated  and is
         therefore the  maximum  order of  autoregressive model  which  can  be
         selected;

 FPEJ    is the Akaike final prediction error, normalized so the  minimum final
         prediction error is one, that is,

                  S(J)**2*(N+J+1)
         FPEJ = ------------------- ;
                S(IAR)**2*(N+IAR+1)

 S(J)**2 is the  residual or  one  step  prediction  variance of  the  order  J
         autoregressive model, computed as

                     N                 J
         S(J)**2 = ----- * [ACOV(1) - SUM PHI(P)*ACOV(P+1)];
                   N-J-1              P=1

 S(IAR)**2 is the residual or one step  prediction  variance  of  the  selected
         model;

 IAR     is the order of the AR model which produced the minimum AIC;

 N       is the number of observations in the series.

 The AIC of the selected model order will always be zero.

      The partial F-test of the change in the residual sum of squares is a test
 of  the null hypothesis that the  order J pacf is zero.   The F-ratio is given
 by

                           S(J-1)**2 - S(J)**2   (N-1-J)*PHI(J)**2
              F(1,N-J-1) = ------------------- = -----------------
                                 S(J)**2            1-PHI(J)**2

 where  (S(J-1)**2 - S(J)**2) is the change in the residual variance due to the
 hypothesis that PHI(J) = 0 (1 degree of freedom);  and S(J)**2 is the residual
 variance of the order J autoregressive model (N-J-1 degrees of freedom).

 The significance level of each F-ratio is also listed.


      The  Cross Correlation Subroutines.  The argument controlling the printed
 output,  NPRT,  is  discussed  in  section  D.   The  output  from  the  cross
 correlation  analysis  subroutines  includes summary statistics describing the
 two input series;  the  cross  correlation  coefficient  estimates  and  their
 standard errors; and a plot of the cross correlation estimates.




                                    <12-38>
1     The  cross correlation function  (ccf) estimate  for the  pair  of series
 stored in the jth and kth columns of YM is

                                      CCOV(i+1,j,k)
                     rjk(i) = -----------------------------
                              sqrt(CCOV(1,j,j)*CCOV(1,k,k))

 for lag i = 0, ..., LAGMAX, j = 1, ..., M, and k = 1, ..., M, where

                               rjk(i) = rkj(-i).

      The estimated standard errors (SE) of the ccf are computed as

         SE[rjk(i)] = SE[rkj(-i)]

                                                q
                      (N-i) * sqrt[(N-i) + 2 * SUM (N-i-v)*rjj(v)*rkk(v)]
                                               v=1
                    = ---------------------------------------------------
                                       N*NLPPC(i+1,j,k)

 for lag i = 0,  ..., LAGMAX and for series  j = 1, ...,  M and k = 1,  ..., M,
 where  q is  the  smaller of  N - i and  i - 1.   This standard  error formula
 assumes that  there are no  cross correlations  between the  two  series being
 compared.   It does  not assume the two series  are white noise, although when
 both series  have been  prewhitened  the summed  portion of  the  equation  is
 approximately zero and the  results will  approximately equal N-1-2,  which is
 the  standard error assuming white noise.   However, the importance of 
 prewhitening before performing cross  correlation analysis  is clear from  the above
 equation.   The presence  of autocorrelation in either one  or both series can
 cause  a significant  increase  in  the  variance  of  the  cross  correlation
 estimates, and,  as shown  in the example  on page  338 of  Jenkins  and Watts
 [1968], the cross correlation function estimates can become  quite unreliable.
 Autocorrelation in the input series should be suspected if the  standard error
 estimates computed by the subroutine are not approximately N-1-2.   It is best
 to use ACF to routinely check  for autocorrelation in the input  series before
 performing cross  correlation  analyses  and  to  prewhiten  the  series  when
 necessary.  The digital filtering subroutines [chapter 10] or the 
 autoregressive model  subroutines can be  used for  prewhitening the input  series before
 cross correlation analysis is performed.


 E.2.b  Spectrum Analysis

      The  Univariate  Fourier  Spectrum  Analysis  Subroutines.  The  argument
 controlling the printed output,  NPRT,  is discussed in section D.  The output
 from  the  univariate  Fourier  spectrum analysis subroutines consists of four
 spectrum  plots  with  successively  narrower  bandwidths.  Each  spectrum  is
 plotted  either  in  decibels  (10  times  the  base 10 logarithm of the power
 spectrum),  scaled so that  the  maximum  value  plotted  is  zero,  or  on  a
 logarithmic  scale.  The  95-percent confidence interval and the bandwidth for
 each spectrum are shown on the plots.

      The bandwidth  of  a  Parzen  lag  window  with  truncation  point  T  is
 approximately  1.85/T  for  T<<N  and a low percentage of missing values.  The
 actual bandwidth is



                                    <12-39>
1                                        1
                          BW = -----------------------
                                T   W(P)**2*(N-|P|)**2
                               SUM  ------------------
                               P=-T     N*NLPPA(P)

 where W(P) is the Parzen lag  window  function  defined  in  section  D  under
 argument  SPCF.  The  bandwidth  is  indicated by the displayed "B * W" in the
 upper-right portion of the plot.

      The effective  degrees of freedom  (edf) in  a spectrum estimate  using a
 Parzen lag window is

                                 edf = 2*BW*N.

 For large values  of N  and a low  percentage of  missing values  this  can be
 approximated by 3.71*N/T.

      A 95-percent confidence interval for the estimated spectrum is

                  SPCF(I)*edf       SPCF(I)*edf
                --------------- , ---------------   for I = 1, ..., NF,
                CHI2(edf,0.975)   CHI2(edf,0.025)

 where CHI2(edf,alpha) is the alpha-percent point  function  for  a  chi-square
 distribution with edf degrees of freedom.  Note that when the logarithm of the
 spectrum is plotted this confidence interval has constant width over frequency
 although  it  is  not symmetric.  In this case,  the lower confidence interval
 limit is

                  SPCF(I)*edf                             edf
            log --------------- = log[SPCF(I)] + log ---------------
                CHI2(edf,0.975)                      CHI2(edf,0.975)

 and the upper confidence interval limit is

                  SPCF(I)*edf                             edf
            log --------------- = log[SPCF(I)] + log ---------------
                CHI2(edf,0.025)                      CHI2(edf,0.025)

 for I = 1, ..., NF.

      The  width of  the  confidence interval  is indicated  by  the  displayed
 "C * I" aligned  vertically in  the  upper-right  portion  of the  plot.   The
 asterisk (*)  separates the upper  limit from  the lower limit.   The user  is
 reminded that the confidence interval applies to individual  frequency points,
 not to the spectrum as a whole.


      The Autoregressive Spectrum Subroutines.  The  argument  controlling  the
 printed  output,  NPRT,  is  discussed  in  section  D.  The  output  from the
 autoregressive spectrum  analysis  subroutines  consists  of  a  plot  of  the
 autoregressive  and  Fourier  spectra either in decibels (10 times the base 10
 logarithm of the power spectrum) scaled such that the maximum value plotted is
 zero or on a logarithmic scale.  The output also includes  the  autoregressive
 order  selection  statistics  when  the  user does not supply the value of the
 order via argument IAR.  The  bandwidth  and  the  length  of  the  95-percent


                                    <12-40>
1confidence  interval  for  the  Fourier  spectrum  are shown on the plot.  The
 bandwidth  is  not  relevant  to  the  autoregressive  spectrum.  The  Fourier
 spectrum  and  its  confidence  intervals  should  be used in interpreting the
 autoregressive spectrum since the confidence intervals for the  autoregressive
 spectrum are not computed.

      The autoregressive order selection statistics include the modified Akaike
 information  criteria (AIC);  the results of a partial F-test of the change in
 the  residual  sum  of  squares  with  the  addition  of   the   most   recent
 autoregressive parameter;  the order of the autoregressive model selected; and
 the Yule-Walker estimates of the parameters of the selected model.

      The modified Akaike information criteria (AIC) is

                  AIC(J) = N*ln(FPEJ) for J = 0, ..., LAGMAX,

 where

 LAGMAX  is  the maximum  lag  value  for  which the  autocovariance  has  been
         computed, and is therefore also the maximum order autoregressive model
         which can be selected;

 FPEJ    is the Akaike final  prediction error, normalized so that  the minimum
         final prediction variance S(IAR)**2 is one, that is,

                  (N+J+1)*S(J)**2
         FPEJ = ------------------- ;
                (N+IAR+1)*S(IAR)**2

 S(J)**2 is the residual or one step prediction variance of the order J model;

 S(IAR)**2 is the residual or one step  prediction  variance  of  the  selected
                model;

 IAR     is  the order of the model  which produced the minimum AIC,  i.e., the
         order  having  the lowest  Akaike   final  prediction  error  [Akaike,
         1974].

 The AIC of the selected model order will always be zero.

      The partial F-test of the change in the residual sum of squares is a test
 of the null hypothesis that the Jth autoregressive coefficient of the  order J
 model is zero.  The F-ratio is given by

                                     S(J-1)**2 - S(J)**2
                        F(1,N-J-1) = -------------------
                                           S(J)**2

 where  (S(J-1)**2 - S(J)**2) is the change in the residual variance due to the
 hypothesis that PHI(J) = 0  (one  degree  of  freedom);  and  S(J)**2  is  the
 residual  variance  of  the  order  J  autoregressive  model (N-J-1 degrees of
 freedom).  The significance level of each F-ratio is also listed.

      The estimated  Fourier spectrum and its bandwidth and confidence interval
 are computed as described above.




                                    <12-41>
1     The Direct Fourier Transform Subroutines.  The argument  controlling  the
 printed  output,  NPRT,  is  discussed  in  section  D.  The  output  from the
 periodogram subroutines  consists  of  a  plot  of  the  computed  periodogram
 displayed  either  in  decibels  (10  times  the  base  10  logarithm  of  the
 periodogram estimates),  scaled so that the maximum value plotted is zero,  or
 on   a   logarithmic  scale.   The  output  from  the  integrated  periodogram
 subroutines  consists  of  a  one-page  plot  of  the  integrated  periodogram
 accompanied  by 95-percent confidence bands for testing the null hypothesis of
 white noise.  These  bands  are  the  Kolmogoroff-Smirnov  probability  limits
 applicable to cumulative distribution functions [Hald, 1952].


      The  Bivariate  Fourier  Spectrum.  The  argument controlling the printed
 output,  NPRT,  is discussed in section  D.  The  output  from  the  bivariate
 Fourier  spectrum  subroutines consists of four spectrum plot pairs (a squared
 coherency plot and a phase plot) with successively narrower bandwidths  chosen
 by  the window-closing algorithm.  The 95-percent confidence intervals and the
 95-percent significance levels are shown on the coherency plots.

      The bandwidth  of  a  Parzen  lag  window  with  truncation  point  T  is
 approximately  3.71*N/T  for T<<N and a low percentage of missing values.  The
 actual bandwidth is

                                         1
                          BW = -----------------------
                                T   W(P)**2*(N-|P|)**2
                               SUM  ------------------
                               P=-T     N*NLPPA(P)

 NLPPC(k) is the number of lagged product pairs used to compute the covariance
 function for lag k, that is

                                   N-|k|
                        NLPPC(k) =  SUM  u1(t)*u2(t-|k|)
                                    t=1

 where u(t) is an indicator variable:   u(t) = 1 if Y(t) has been observed, and
 u(t) = 0 if Y(t) is missing.


 F.  Examples

      Autocorrelation Analysis.  In the first example  program  below,  ACF  is
 used to compute and plot the autocorrelations,  partial autocorrelations,  and
 autoregressive order selection statistics of the input  series  Y,  where  the
 data used is the series X1,  a simulated first order autoregressive model with
 PHI(1) = 0.6, listed on page 362 in Jenkins and Watts [1968].  The correlation
 estimates for this series are given on page 421  of  Jenkins  and  Watts.  The
 first  two pages of output from this example show the autocorrelations and the
 second two  pages  show  the  partial  autocorrelations.  Note  that  for  all
 correlation plots the value of the lag,  i, is is shown on the left margin and
 the actual correlation coefficients are shown on the right margin of  each  of
 the plots.


      Cross  Correlation Analysis.  In the second example program,  CCF is used
 to compute and plot the cross correlations between the first 50 values of  the


                                    <12-42>
1series  X1  and  X2  listed  on  page  361  of  Jenkins and Watts [1968].  The
 correlation estimates for these series are given on page 420  of  Jenkins  and
 Watts.  Note that for all correlation plots the value of the lag,  i, is shown
 on the left margin and the actual correlation coefficients are  shown  on  the
 right margin of each of the plots.


      Univariate Fourier Spectrum Analysis.  In the third example program,  UFS
 is  used to compute and plot the univariate Fourier spectrum estimates for the
 first 50 values of the series listed on page 318 of Jenkins and Watts  [1968];
 the  spectrum  of  this  series is shown for various bandwidths on page 270 of
 Jenkins and Watts [1968].  The four plots  produced  by  UFS  show  the  power
 spectrum  (in  decibels  versus frequency).  Each spectrum plotted is computed
 using a narrower bandwidth than the preceding one.  The bandwidth,  lag window
 truncation  point and effective degrees of freedom for the spectrum are listed
 at the top of each page.  The bandwidth and the 95-percent confidence interval
 of the individual spectrum values are shown graphically on each plot.


      Autoregressive Spectrum Analysis.  In the fourth example program,  UAS is
 used  to  compute  and plot the univariate autoregressive and Fourier spectrum
 estimates for the first 50 values listed on page  318  of  Jenkins  and  Watts
 [1968].  The  theoretical and Fourier spectra of this series are shown on page
 270 of Jenkins and Watts [1968].  The first page of output from  this  example
 lists  the  autoregressive order selection statistics while the second gives a
 plot of the autoregressive and Fourier  spectra.  The  bandwidth,  lag  window
 truncation  point and effective degrees of freedom of the Fourier spectrum are
 listed at the top of the plot.  The 95-percent  confidence  interval  and  the
 bandwidth  of  the Fourier spectrum are displayed on the plot.  The concept of
 bandwidth does  not  apply  to  the  autoregressive  spectrum  and  confidence
 intervals for the autoregressive spectrum are not computed.


      Direct  Fourier Transform of a Univariate Series and Utility Subroutines.
 In the fifth example program below,  TAPER is used to center the input  series
 and  to  taper  10 percent of the data at the ends of the series (5 percent at
 each end) and PGMS is used to compute  the  raw  periodogram  of  the  tapered
 series.  MDFLT  is  then  used  to smooth the raw periodogram returned by PGMS
 with three applications of a modified Daniell  filter  of  length  eight.  PPL
 (chapter  3)  is  used  to  produce a log plot the smoothed periodogram versus
 frequency.  (VPL could also have been used to display the results.)  The  data
 used  are the Wolf sunspot numbers from 1700 to 1960 as tabulated by Waldmeier
 [1961].  The raw and smoothed periodograms of the tapered series are shown  on
 pages 95 and 176,  respectively,  of Bloomfield [1976].  Note that there is no
 printed output from TAPER or MDFLT.


      Bivariate Fourier Spectrum Analysis.  In the last example program of this
 chapter,  BFS is used to compute and plot the cross spectrum for the series X1
 and X2 listed on page 361 of Jenkins and Watts [1968].  The squared  coherency
 and  phase  spectra for these series are shown on pages 387 and 388 of Jenkins
 and Watts.







                                    <12-43>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE ACF USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(300)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL ACF FOR AUTOCORRELATION ANALYSIS
 C
       WRITE (IPRT,102)
       CALL ACF (Y, N)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (10F6.2)
   102 FORMAT ('1RESULTS OF STARPAC AUTOCORRELATION SUBROUTINE ACF')
       END

 Data:

   100
  -2.07 -1.15  0.69 -0.46 -1.49 -0.70 -1.07 -0.69 -0.68  1.27
  -1.05 -0.05 -0.84 -0.62 -0.49 -1.29 -0.49 -1.06 -0.38 -0.52
  -0.13  1.30 -1.51 -0.43 -1.33 -0.78  0.31 -0.95 -0.90 -0.30
  -1.02 -0.53  0.15  1.40  1.22  0.59  0.70  1.70  2.78  1.98
   1.39  1.85  2.60  0.51  2.77  1.16  1.07 -0.48 -0.52  0.37
   0.00 -1.99 -1.75  0.70  0.73  1.16  0.06 -0.02  1.10 -0.35
  -1.67 -1.57  1.16  1.84  3.35  0.40  0.45  1.30  0.93  1.17
  -1.74 -1.28 -0.07  1.50  0.53  0.20 -0.42  1.18  0.82  1.50
   2.92  1.18  1.23  3.16  0.79  0.68  1.14  1.02  1.02 -0.71
  -0.17 -1.50 -0.26 -0.38  0.93 -0.33 -1.12 -2.95 -2.09 -1.11










                                    <12-44>
1RESULTS OF STARPAC AUTOCORRELATION SUBROUTINE ACF
                                                                                                         STARPAC 2.08S (03/15/90)
 AUTOCORRELATION ANALYSIS

 AVERAGE OF THE SERIES                 =   .1450000
 STANDARD DEVIATION OF THE SERIES      =   1.277758
 NUMBER OF TIME POINTS                 =        100
 LARGEST LAG VALUE USED                =         33



 AUTOCORRELATION FUNCTION ESTIMATE (ACF)


 LAG                     1      2      3      4      5      6      7      8      9     10     11     12
 ACF                   .54    .34    .22    .26    .27    .15    .06    .03    .06    .02   -.04   -.11
 STANDARD ERROR        .10    .12    .13    .13    .14    .14    .14    .14    .14    .14    .14    .14

 LAG                    13     14     15     16     17     18     19     20     21     22     23     24
 ACF                  -.14   -.10   -.10   -.11   -.12   -.07    .03    .04    .02   -.03    .02    .06
 STANDARD ERROR        .14    .14    .14    .14    .14    .14    .14    .14    .14    .14    .13    .13

 LAG                    25     26     27     28     29     30     31     32     33
 ACF                   .05   -.03   -.04   -.06   -.06   -.03   -.09   -.16   -.17
 STANDARD ERROR        .13    .13    .13    .13    .13    .13    .13    .13    .13



 THE CHI SQUARE TEST STATISTIC OF
 THE NULL HYPOTHESIS OF WHITE NOISE    =            79.54
 DEGREES OF FREEDOM                    =               33
 OBSERVED SIGNIFICANCE LEVEL           =            .0000



























                                                             <12-45>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 AUTOCORRELATION FUNCTION ESTIMATE (ACF)

             -1.0000    -.8000    -.6000    -.4000    -.2000     .0000     .2000     .4000     .6000     .8000    1.0000
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                                                  ++++++++++++++++++++++++++++                       I   .53974
  2.0000      I                                                  ++++++++++++++++++                                 I   .33598
  3.0000      I                                                  ++++++++++++                                       I   .21955
  4.0000      I                                                  ++++++++++++++                                     I   .26190
  5.0000      I                                                  +++++++++++++++                                    I   .27456
  6.0000      I                                                  +++++++++                                          I   .15137
  7.0000      I                                                  ++++                                               I   .56358E-01
  8.0000      I                                                  +++                                                I   .30065E-01
  9.0000      I                                                  ++++                                               I   .58569E-01
  10.000      I                                                  ++                                                 I   .17475E-01
  11.000      I                                                +++                                                  I  -.35159E-01
  12.000      I                                            +++++++                                                  I  -.11314
  13.000      I                                           ++++++++                                                  I  -.13858
  14.000      I                                             ++++++                                                  I  -.10145
  15.000      I                                             ++++++                                                  I  -.10083
  16.000      I                                             ++++++                                                  I  -.10666
  17.000      I                                            +++++++                                                  I  -.12337
  18.000      I                                               ++++                                                  I  -.65475E-01
  19.000      I                                                  ++                                                 I   .26352E-01
  20.000      I                                                  +++                                                I   .37046E-01
  21.000      I                                                  ++                                                 I   .24165E-01
  22.000      I                                                 ++                                                  I  -.27097E-01
  23.000      I                                                  ++                                                 I   .15920E-01
  24.000      I                                                  ++++                                               I   .58745E-01
  25.000      I                                                  ++++                                               I   .52299E-01
  26.000      I                                                +++                                                  I  -.30013E-01
  27.000      I                                                +++                                                  I  -.42630E-01
  28.000      I                                               ++++                                                  I  -.61556E-01
  29.000      I                                               ++++                                                  I  -.55425E-01
  30.000      I                                                 ++                                                  I  -.28172E-01
  31.000      I                                             ++++++                                                  I  -.91703E-01
  32.000      I                                          +++++++++                                                  I  -.15668
  33.000      I                                         ++++++++++                                                  I  -.17497




















                                                             <12-46>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 PARTIAL AUTOCORRELATION FUNCTION ESTIMATE (PACF)

 AND AUTOREGRESSIVE ORDER SELECTION STATISTICS




 LAG                     1      2      3      4      5      6      7      8      9     10     11     12
 PACF                  .54    .06    .02    .17    .09   -.11   -.05    .00    .02   -.06   -.03   -.08
 AIC                   .00   1.60   3.56   2.53   3.70   4.59   6.32   8.34  10.32  11.94  13.85  15.18
 F RATIO             40.28    .39    .05   2.92    .80   1.04    .26    .00    .03    .36    .10    .62
 F PROBABILITY         .00    .54    .83    .09    .37    .31    .61    .99    .87    .55    .75    .43

 LAG                    13     14     15     16     17     18     19     20     21     22     23     24
 PACF                 -.06    .01   -.02   -.01   -.01    .07    .10    .00    .01   -.05    .02    .01
 AIC                 16.81  18.83  20.83  22.88  24.94  26.56  27.63  29.72  31.80  33.65  35.72  37.83
 F RATIO               .35    .02    .04    .01    .00    .36    .81    .00    .01    .20    .04    .01
 F PROBABILITY         .55    .89    .84    .93    .94    .55    .37    .98    .91    .66    .84    .90

 LAG                    25     26     27     28     29     30     31     32     33
 PACF                 -.03   -.10    .00   -.08   -.05    .04   -.06   -.10   -.03
 AIC                 39.87  41.10  43.26  44.82  46.78  48.80  50.70  52.03  54.21
 F RATIO               .07    .68    .00    .44    .16    .13    .22    .61    .05
 F PROBABILITY         .80    .41    .97    .51    .69    .72    .64    .44    .83


 ORDER AUTOREGRESSIVE PROCESS SELECTED =              1
 ONE STEP PREDICTION VARIANCE OF PROCESS SELECTED =  1.1688519

 YULE-WALKER ESTIMATES OF THE COEFFICIENTS OF THE AUTOREGRESSIVE PROCESS SELECTED

 COEFFICIENT NUMBER      1
 COEFFICIENT VALUE   .5397
























                                                             <12-47>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 PARTIAL AUTOCORRELATION FUNCTION ESTIMATE (PACF)

             -1.0000    -.8000    -.6000    -.4000    -.2000     .0000     .2000     .4000     .6000     .8000    1.0000
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
  1.0000      I                                                  ++++++++++++++++++++++++++++                       I   .53974
  2.0000      I                                                  ++++                                               I   .63018E-01
  3.0000      I                                                  ++                                                 I   .22135E-01
  4.0000      I                                                  ++++++++++                                         I   .17266
  5.0000      I                                                  ++++++                                             I   .91596E-01
  6.0000      I                                             ++++++                                                  I  -.10523
  7.0000      I                                               ++++                                                  I  -.53166E-01
  8.0000      I                                                  +                                                  I   .93003E-03
  9.0000      I                                                  ++                                                 I   .17874E-01
  10.000      I                                               ++++                                                  I  -.63464E-01
  11.000      I                                                +++                                                  I  -.33574E-01
  12.000      I                                              +++++                                                  I  -.84131E-01
  13.000      I                                               ++++                                                  I  -.63805E-01
  14.000      I                                                  ++                                                 I   .14446E-01
  15.000      I                                                 ++                                                  I  -.21400E-01
  16.000      I                                                  +                                                  I  -.92888E-02
  17.000      I                                                  +                                                  I  -.76472E-02
  18.000      I                                                  ++++                                               I   .66474E-01
  19.000      I                                                  ++++++                                             I   .10021
  20.000      I                                                  +                                                  I   .30543E-02
  21.000      I                                                  ++                                                 I   .13109E-01
  22.000      I                                               ++++                                                  I  -.50444E-01
  23.000      I                                                  ++                                                 I   .22769E-01
  24.000      I                                                  ++                                                 I   .13963E-01
  25.000      I                                                +++                                                  I  -.30205E-01
  26.000      I                                             ++++++                                                  I  -.96142E-01
  27.000      I                                                  +                                                  I  -.47261E-02
  28.000      I                                              +++++                                                  I  -.78259E-01
  29.000      I                                                +++                                                  I  -.48012E-01
  30.000      I                                                  +++                                                I   .42976E-01
  31.000      I                                               ++++                                                  I  -.56478E-01
  32.000      I                                             ++++++                                                  I  -.95175E-01
  33.000      I                                                 ++                                                  I  -.27238E-01




















                                                             <12-48>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE CCF USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y1 AND Y2 MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y1(100), Y2(100)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          SERIES 1
 C          SERIES 2
 C
       READ (5,100) N
       READ (5,101) (Y1(I), I=1,N)
       READ (5,101) (Y2(I), I=1,N)
 C
 C     PRINT TITLE AND CALL CCF FOR CROSS CORRELATION ANALYSIS
 C
       WRITE (IPRT,102)
       CALL CCF (Y1, Y2, N)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (12F6.2)
   102 FORMAT ('1RESULTS OF STARPAC CROSS CORRELATION SUBROUTINE CCF')
       END

 Data:

    50
  -0.88 -0.16 -1.87 -1.12  1.38  2.13  2.76  0.56 -0.69 -1.79 -3.82 -2.38
   1.00  0.70 -0.15  0.98  0.11 -0.35 -0.73  0.89 -1.63 -0.44 -1.37 -1.71
  -1.22 -2.00 -0.22  0.38  1.31  0.71  0.32  0.48 -1.88 -0.94 -1.54 -0.13
   1.02  0.02 -0.77  0.11 -0.60 -0.52 -0.09  1.23  1.46  0.61  0.42  2.16
   3.18  2.10
   0.79  1.12 -1.10 -2.39 -1.75 -0.82 -0.36  1.27  1.75  2.44  0.36 -2.10
  -1.93 -1.30 -1.75 -0.34  0.74  0.49  0.70  0.71  0.09  0.59  1.54  0.14
   0.55 -1.40 -2.55 -1.66 -0.43  0.58  2.18 -0.24  0.58 -0.18 -1.55 -0.64
  -1.09  0.90 -0.66 -0.35  0.48  0.50  0.05 -0.68  0.24  0.58 -1.26 -0.25
   0.25  2.18








                                    <12-49>
1RESULTS OF STARPAC CROSS CORRELATION SUBROUTINE CCF
                                                                                                         STARPAC 2.08S (03/15/90)
 CROSS CORRELATION ANALYSIS

                                          SERIES  1     SERIES  2

 AVERAGE OF THE SERIES                 =  -.5960000E-01 -.9960000E-01
 STANDARD DEVIATION OF THE SERIES      =   1.399513      1.202491
 NUMBER OF TIME POINTS                 =         50            50

 LARGEST LAG VALUE TO BE USED          =         32



 CROSS CORRELATION FUNCTION ESTIMATE (CCF)

 CCF CORRELATES SERIES 1 AT TIME T WITH SERIES 2 AT TIME T + K.
     (IF PEAK CORRELATION OCCURES AT POSITIVE (NEGATIVE) LAG
        THEN SERIES 1 LEADS (LAGS) SERIES 2)

 LAG                   -25    -26    -27    -28    -29    -30    -31    -32
 CCF                  -.10    .06    .16    .16    .15    .11    .00   -.10
 STANDARD ERROR        .15    .14    .14    .14    .13    .13    .13    .12

 LAG                   -13    -14    -15    -16    -17    -18    -19    -20    -21    -22    -23    -24
 CCF                  -.07   -.05   -.11   -.19   -.09    .08    .18    .17    .09   -.06   -.12   -.14
 STANDARD ERROR        .18    .18    .18    .18    .17    .17    .17    .16    .16    .16    .15    .15

 LAG                    -1     -2     -3     -4     -5     -6     -7     -8     -9    -10    -11    -12
 CCF                  -.52   -.60   -.33   -.05    .15    .37    .35    .13   -.05   -.13    .00    .04
 STANDARD ERROR        .21    .21    .21    .21    .20    .20    .20    .20    .19    .19    .19    .19

 LAG                     0
 CCF                  -.02
 STANDARD ERROR        .14

 LAG                     1      2      3      4      5      6      7      8      9     10     11     12
 CCF                   .41    .51    .48    .33    .01   -.18   -.36   -.26   -.04   -.02    .06    .02
 STANDARD ERROR        .21    .21    .21    .21    .20    .20    .20    .20    .19    .19    .19    .19

 LAG                    13     14     15     16     17     18     19     20     21     22     23     24
 CCF                   .04    .05    .08    .20    .14    .12   -.19   -.27   -.20   -.14    .07    .17
 STANDARD ERROR        .18    .18    .18    .18    .17    .17    .17    .16    .16    .16    .15    .15

 LAG                    25     26     27     28     29     30     31     32
 CCF                   .19    .10   -.08   -.08   -.10   -.11   -.07    .02
 STANDARD ERROR        .15    .14    .14    .14    .13    .13    .13    .12












                                                             <12-50>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 CROSS CORRELATION FUNCTION ESTIMATE (CCF)

 CCF CORRELATES SERIES 1 AT TIME T WITH SERIES 2 AT TIME T + K.
     (IF PEAK CORRELATION OCCURES AT POSITIVE (NEGATIVE) LAG
        THEN SERIES 1 LEADS (LAGS) SERIES 2)
             -1.0000    -.8000    -.6000    -.4000    -.2000     .0000     .2000     .4000     .6000     .8000    1.0000
              -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
 -32.000      I                                             ++++++                                                  I  -.98862E-01
 -31.000      I                                                  +                                                  I  -.28858E-02
 -30.000      I                                                  ++++++                                             I   .10712
 -29.000      I                                                  ++++++++                                           I   .14956
 -28.000      I                                                  +++++++++                                          I   .16054
 -27.000      I                                                  +++++++++                                          I   .16150
 -26.000      I                                                  ++++                                               I   .60338E-01
 -25.000      I                                             ++++++                                                  I  -.10388
 -24.000      I                                           ++++++++                                                  I  -.14362
 -23.000      I                                            +++++++                                                  I  -.12479
 -22.000      I                                               ++++                                                  I  -.55619E-01
 -21.000      I                                                  +++++                                              I   .88452E-01
 -20.000      I                                                  +++++++++                                          I   .16802
 -19.000      I                                                  ++++++++++                                         I   .18040
 -18.000      I                                                  +++++                                              I   .77358E-01
 -17.000      I                                             ++++++                                                  I  -.92188E-01
 -16.000      I                                        +++++++++++                                                  I  -.19470
 -15.000      I                                            +++++++                                                  I  -.11176
 -14.000      I                                               ++++                                                  I  -.50033E-01
 -13.000      I                                               ++++                                                  I  -.69442E-01
 -12.000      I                                                  +++                                                I   .37527E-01
 -11.000      I                                                  +                                                  I   .43801E-02
 -10.000      I                                            +++++++                                                  I  -.12992
 -9.0000      I                                               ++++                                                  I  -.53960E-01
 -8.0000      I                                                  ++++++++                                           I   .13037
 -7.0000      I                                                  ++++++++++++++++++                                 I   .34581
 -6.0000      I                                                  +++++++++++++++++++                                I   .36962
 -5.0000      I                                                  ++++++++                                           I   .14975
 -4.0000      I                                               ++++                                                  I  -.52573E-01
 -3.0000      I                                  +++++++++++++++++                                                  I  -.32785
 -2.0000      I                    +++++++++++++++++++++++++++++++                                                  I  -.59800
 -1.0000      I                        +++++++++++++++++++++++++++                                                  I  -.51544
 0.           I                                                 ++                                                  I  -.17402E-01
  1.0000      I                                                  ++++++++++++++++++++++                             I   .41146
  2.0000      I                                                  +++++++++++++++++++++++++++                        I   .51130
  3.0000      I                                                  +++++++++++++++++++++++++                          I   .48166
  4.0000      I                                                  ++++++++++++++++++                                 I   .33171
  5.0000      I                                                  ++                                                 I   .12280E-01
  6.0000      I                                         ++++++++++                                                  I  -.18197
  7.0000      I                                +++++++++++++++++++                                                  I  -.35999
  8.0000      I                                     ++++++++++++++                                                  I  -.25903
  9.0000      I                                                +++                                                  I  -.41510E-01
  10.000      I                                                 ++                                                  I  -.16659E-01
  11.000      I                                                  ++++                                               I   .62016E-01
  12.000      I                                                  ++                                                 I   .18025E-01
  13.000      I                                                  +++                                                I   .40650E-01
  14.000      I                                                  ++++                                               I   .52799E-01
  15.000      I                                                  +++++                                              I   .79604E-01
  16.000      I                                                  +++++++++++                                        I   .20153
  17.000      I                                                  ++++++++                                           I   .13843
                                                             <12-51>
1 18.000      I                                                  +++++++                                            I   .11858
  19.000      I                                        +++++++++++                                                  I  -.19467
  20.000      I                                    +++++++++++++++                                                  I  -.27015
  21.000      I                                        +++++++++++                                                  I  -.19730
  22.000      I                                           ++++++++                                                  I  -.14359
  23.000      I                                                  +++++                                              I   .73960E-01
  24.000      I                                                  ++++++++++                                         I   .17440
  25.000      I                                                  +++++++++++                                        I   .19101
  26.000      I                                                  ++++++                                             I   .99404E-01
  27.000      I                                              +++++                                                  I  -.84161E-01
  28.000      I                                              +++++                                                  I  -.76186E-01
  29.000      I                                             ++++++                                                  I  -.10048
  30.000      I                                             ++++++                                                  I  -.10986
  31.000      I                                               ++++                                                  I  -.66356E-01
  32.000      I                                                  ++                                                 I   .15525E-01












































                                                             <12-52>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE UFS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(300)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL UFS FOR UNIVARIATE FOURIER SPECTRUM ANALYSIS
 C
       WRITE (IPRT,102)
       CALL UFS (Y, N)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (10F6.2)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' UNIVARIATE FOURIER SPECTRUM ANALYSIS SUBROUTINE UFS')
       END

 Data:

    50
  -0.88 -0.12 -0.89 -1.38 -0.07  1.03  2.14  0.35 -1.10 -1.78
  -2.76 -1.77  0.98  1.00 -0.70 -1.01 -1.30 -0.85 -0.46  1.63
   0.06 -0.17 -1.01 -1.04 -0.66 -1.12 -0.51 -0.71 -0.20 -0.13
   0.14  1.59 -0.76 -1.08 -1.77 -1.20  0.45 -0.07 -0.63 -0.35
  -0.87 -0.62  0.28  1.90  2.14  1.05  0.31  1.07  2.67  2.44














                                    <12-53>
1RESULTS OF STARPAC UNIVARIATE FOURIER SPECTRUM ANALYSIS SUBROUTINE UFS
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    4 / BW= .4685 / EDF=    47)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .0000 - ++++++++++                                           C                                                -
                I           +++++++                                                                                     I
                I                  ++++                                                                                 I
                I                      ++++                                                                             I
                I                          +++                                                                          I
         -.9096 -                             +++                                                                       -
                I                                ++                                                                     I
                I                                  +++                                                                  I
                I                                     ++                                                                I
                I                                       ++                                                              I
        -1.8193 -                                         ++                                                            -
                I       B                                   ++         *                                              W I
                I                                             +                                                         I
                I                                              ++                                                       I
                I                                                ++                                                     I
        -2.7289 -                                                  +                                                    -
                I                                                   ++                                                  I
                I                                                     +                                                 I
                I                                                      +                                                I
                I                                                       ++                                              I
        -3.6385 -                                                      I  +                                             -
                I                                                          +                                            I
                I                                                           ++                                          I
                I                                                             +                                         I
                I                                                              +                                        I
        -4.5482 -                                                               +                                       -
                I                                                                +                                      I
                I                                                                 ++                                    I
                I                                                                   +                                   I
                I                                                                    +                                  I
        -5.4578 -                                                                     +                                 -
                I                                                                      +                                I
                I                                                                       +                               I
                I                                                                        ++                             I
                I                                                                          +                            I
        -6.3674 -                                                                           +                           -
                I                                                                            +                          I
                I                                                                             +                         I
                I                                                                              +                        I
                I                                                                               ++                      I
        -7.2771 -                                                                                 +                     -
                I                                                                                  +                    I
                I                                                                                   +                   I
                I                                                                                    ++                 I
                I                                                                                      +                I
        -8.1867 -                                                                                       ++              -
                I                                                                                         +             I
                I                                                                                          ++           I
                I                                                                                            ++         I
                I                                                                                              +++      I
        -9.0963 -                                                                                                 +++++ -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-54>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    8 / BW= .2381 / EDF=    24)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .0000 - +++++++++++++++++++++++                                                     C                         -
                I                        ++++++                                                                         I
                I                              +++                                                                      I
                I                                 ++                                                                    I
                I                                   ++                                                                  I
        -1.1087 -                                     ++                                                                -
                I                                       +                                                               I
                I                                        ++                                                             I
                I                                          +                                                            I
                I                                           +                                                           I
        -2.2173 -                                            +                                                          -
                I                                             +                                                         I
                I                                              +                                                        I
                I                                               +     B                       *                       W I
                I                                                +                                                      I
        -3.3260 -                                                 +                                                     -
                I                                                                                                       I
                I                                                  +                                                    I
                I                                                   +                                                   I
                I                                                    +                                                  I
        -4.4346 -                                                     +                                                 -
                I                                                                                                       I
                I                                                      +                                                I
                I                                                       +                     I                         I
                I                                                        +                                              I
        -5.5433 -                                                                                                       -
                I                                                         +                                             I
                I                                                          +                                            I
                I                                                           +                                           I
                I                                                                                                       I
        -6.6519 -                                                            +                                          -
                I                                                             +                                         I
                I                                                              +                                        I
                I                                                               +                                       I
                I                                                                +                                      I
        -7.7606 -                                                                 +                                     -
                I                                                                  +                                    I
                I                                                                   +                                   I
                I                                                                    +                                  I
                I                                                                     +                                 I
        -8.8692 -                                                                      ++                               -
                I                                                                        ++                             I
                I                                                                          +                            I
                I                                                                           +++                         I
                I                                                                              ++                       I
        -9.9779 -                                                                                +++                    -
                I                                                                                   ++                  I
                I                                                                                     +++               I
                I                                                                                        ++++           I
                I                                                                                            ++++       I
       -11.0865 -                                                                                                ++++++ -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-55>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   16 / BW= .1225 / EDF=    12)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .2541 -                                                                                         C             -
                I                             +++++                                                                     I
                I +++                       ++     ++                                                                   I
                I    ++++                 ++         +                                                                  I
                I        +++             +            +                                                                 I
        -1.1943 -           ++         ++              +                                                                -
                I             +++    ++                                                                                 I
                I                ++++                   +                                                               I
                I                                        +                                                              I
                I                                                                                                       I
        -2.6428 -                                         +                                                             -
                I                                                                                                       I
                I                                          +                                                            I
                I                                           +                                                           I
                I                                                                                                       I
        -4.0912 -                                            +                               B            *           W -
                I                                                                                                       I
                I                                             +                                                         I
                I                                                                                                       I
                I                                              +                                                        I
        -5.5397 -                                               +                                                       -
                I                                                +                                                      I
                I                                                                                                       I
                I                                                 +                                                     I
                I                                                  ++                                                   I
        -6.9881 -                                                    +                                    I             -
                I                                                     +                                                 I
                I                                                      ++                                               I
                I                                                        +                                              I
                I                                                         +                                             I
        -8.4366 -                                                          +                                            -
                I                                                           +                                           I
                I                                                            +                                          I
                I                                                             ++                                        I
                I                                                               +                                       I
        -9.8850 -                                                                +                                      -
                I                                                                 +                                     I
                I                                                                  ++           +++++++                 I
                I                                                                    ++      +++       ++               I
                I                                                                      ++++++            +              I
       -11.3335 -                                                                                         +             -
                I                                                                                          +            I
                I                                                                                           +           I
                I                                                                                            +          I
                I                                                                                             +         I
       -12.7819 -                                                                                              +        -
                I                                                                                               +       I
                I                                                                                                +      I
                I                                                                                                 ++    I
                I                                                                                                   +++ I
       -14.2304 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-56>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   32 / BW= .0650 / EDF=     7)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         1.5084 -                                                                                              C        -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                               ++++                                                                    I
         -.4449 -                              +    +                                                                   -
                I                            ++                                                                         I
                I                           +        +                                                                  I
                I +++++++                  +                                                                            I
                I        ++                           +                                                                 I
        -2.3981 -          ++             +                                                                             -
                I            +++         +             +                                                                I
                I               +                                                                                       I
                I                +      +                                                                               I
                I                 +    +                +                                                               I
        -4.3514 -                  ++ +                                                                  B     *      W -
                I                    +                                                                                  I
                I                                                                                                       I
                I                                        +                                                              I
                I                                                                                                       I
        -6.3046 -                                                                                                       -
                I                                         +                                                             I
                I                                                                                                       I
                I                                                +++                                                    I
                I                                          +   ++   ++                                                  I
        -8.2579 -                                           +++       +                                        I        -
                I                                                      +                                                I
                I                                                       +                                               I
                I                                                        ++                                             I
                I                                                          +                                            I
       -10.2111 -                                                           +                                           -
                I                                                                                 +++                   I
                I                                                            +                   +   ++                 I
                I                                                             +                 +      +                I
                I                                                              ++              +        +               I
       -12.1644 -                                                                +++++                   +              -
                I                                                                     +       +           +             I
                I                                                                      +                                I
                I                                                                       +    +             +            I
                I                                                                        +                  +           I
       -14.1176 -                                                                         +++                           -
                I                                                                                            +          I
                I                                                                                                     + I
                I                                                                                             +      +  I
                I                                                                                              +    +   I
       -16.0709 -                                                                                                  +    -
                I                                                                                               +++     I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
       -18.0241 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-57>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE UAS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y(300)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     PRINT TITLE AND CALL UAS FOR
 C     UNIVARIATE AUTOREGRESSIVE SPECTRUM ANALYSIS
 C
       WRITE (IPRT,102)
       CALL UAS (Y, N)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (10F6.2)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' UNIVARIATE AUTOREGRESSIVE SPECTRUM ANALYSIS SUBROUTINE UAS')
       END

 Data:

    50
  -0.88 -0.12 -0.89 -1.38 -0.07  1.03  2.14  0.35 -1.10 -1.78
  -2.76 -1.77  0.98  1.00 -0.70 -1.01 -1.30 -0.85 -0.46  1.63
   0.06 -0.17 -1.01 -1.04 -0.66 -1.12 -0.51 -0.71 -0.20 -0.13
   0.14  1.59 -0.76 -1.08 -1.77 -1.20  0.45 -0.07 -0.63 -0.35
  -0.87 -0.62  0.28  1.90  2.14  1.05  0.31  1.07  2.67  2.44













                                    <12-58>
1RESULTS OF STARPAC UNIVARIATE AUTOREGRESSIVE SPECTRUM ANALYSIS SUBROUTINE UAS
                                                                                                         STARPAC 2.08S (03/15/90)

 AUTOREGRESSIVE ORDER SELECTION STATISTICS


 LAG                     1      2      3      4      5      6      7      8      9     10     11     12
 AIC                  5.08    .00   1.91   3.54   3.90   5.93   7.30   7.35   9.30  11.37  13.46  15.40
 F RATIO             23.53   7.15    .09    .35   1.48    .01    .57   1.68    .10    .01    .02    .14
 F PROBABILITY         .00    .01    .76    .55    .23    .93    .45    .20    .75    .91    .89    .71

 LAG                    13     14     15     16     17     18     19     20     21     22     23     24
 AIC                 16.80  18.98  20.84  22.94  25.12  27.38  25.94  27.94  29.59  31.88  34.44  36.43
 F RATIO               .55    .00    .24    .10    .06    .04   2.37    .24    .45    .12    .00    .32
 F PROBABILITY         .46    .96    .63    .75    .81    .84    .13    .63    .51    .74    .94    .58

 LAG                    25     26     27     28     29     30     31     32
 AIC                 38.98  41.68  44.55  46.93  49.91  53.00  55.64  59.10
 F RATIO               .08    .04    .00    .25    .04    .04    .25    .00
 F PROBABILITY         .79    .85    .99    .62    .85    .85    .63   1.00


 ORDER AUTOREGRESSIVE PROCESS SELECTED =              2
 ONE STEP PREDICTION VARIANCE OF PROCESS SELECTED =  .88356885

 YULE-WALKER ESTIMATES OF THE COEFFICIENTS OF THE AUTOREGRESSIVE PROCESS SELECTED

 COEFFICIENT NUMBER      1      2
 COEFFICIENT VALUE   .7819 -.3634






























                                                             <12-59>
1                                                                                                        STARPAC 2.08S (03/15/90)
 FOURIER SPECTRUM (+) (LAG WIND. TRUNC. PT.=   16 / BW= .1225 / EDF=    12)
 AND ORDER  2 AUTOREGRESSIVE SPECTRUM (.)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
          .0981 -                       ....                                                              C             -
                I                   ....    ...++                                                                       I
                I                ...         ++..+++                                                                    I
                I +++++        ..          ++    .  +                                                                   I
                I      +++ ....          ++       .. +                                                                  I
        -1.3503 -      ...2++           +           . +                                                                 -
                I .....      ++       ++             . +                                                                I
                I              +++++++                . +                                                               I
                I                                      .                                                                I
                I                                       .+                                                              I
        -2.7988 -                                        .+                                                             -
                I                                                                                                       I
                I                                         .+                                                            I
                I                                          .                                                            I
                I                                           2                                                           I
        -4.2472 -                                            .                               B            *           W -
                I                                            +.                                                         I
                I                                             +.                                                        I
                I                                               .                                                       I
                I                                              + .                                                      I
        -5.6957 -                                                 .                                                     -
                I                                               +  .                                                    I
                I                                                +  .                                                   I
                I                                                 +  .                                                  I
                I                                                  +  .                                                 I
        -7.1441 -                                                   +  .                                  I             -
                I                                                    +  .                                               I
                I                                                     ++ .                                              I
                I                                                       + .                                             I
                I                                                        ++..                                           I
        -8.5926 -                                                          + .                                          -
                I                                                           + .                                         I
                I                                                            + ..                                       I
                I                                                             +  .                                      I
                I                                                              +  ..                                    I
       -10.0410 -                                                               +   .                                   -
                I                                                                +   ..                                 I
                I                                                                 ++   ..         +++                   I
                I                                                                   +    ..   ++++   +++                I
                I                                                                    +++   22+          +               I
       -11.4895 -                                                                       +++  ..          ++             -
                I                                                                              ...         +            I
                I                                                                                 ...       +           I
                I                                                                                    ....    +          I
                I                                                                                        .....2         I
       -12.9379 -                                                                                              2....... -
                I                                                                                               +       I
                I                                                                                                +      I
                I                                                                                                 +     I
                I                                                                                                  ++   I
       -14.3864 -                                                                                                    ++ -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-60>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE DIRECT FOURIER TRANSFORM SUBROUTINES USING SINGLE
 C     PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, YFFT, PER, FREQ AND PERF MUST BE CHANGED TO
 C          DOUBLE PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS
 C          USED.
 C
       REAL Y(600), YFFT(600), PER(300), FREQ(300), PERF(300)
       INTEGER KMD(10)
       DOUBLE PRECISION DSTAK(1000)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 1000
       LPER = 300
       LFREQ = 300
       LYFFT = 600
 C
 C     READ NUMBER OF MODIFIED DANIEL FILTERS TO BE APPLIED
 C          FILTER LENGTHS
 C          NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) NK
       READ (5,100) (KMD(I), I = 1, NK)
       READ (5,100) N
       READ (5,101) (Y(I), I=1,N)
 C
 C     CENTER THE SERIES AND
 C     APPLY A TEN PERCENT TAPER TO THE ENDS OF THE DATA
 C
       TAPERP = 0.10
       CALL TAPER (Y, N, TAPERP, YFFT)
 C
 C     PRINT TITLE AND CALL PGMS TO COMPUTE PERIODOGRAM OF TAPERED AND
 C     EXTENDED SERIES.
 C
       WRITE (IPRT,102)
       NFFT = 514
       IEXTND = 0
       NPRT = (-1)
       CALL PGMS (YFFT, N, NFFT, LYFFT, IEXTND, NF, PER, LPER, FREQ,
      1   LFREQ, NPRT)
 C
 C     APPLY MODIFIED DANIELL FILTERS TO SMOOTH THE PERIODOGRAM

                                    <12-61>
1C
       CALL MDFLT (PER, NF, NK, KMD, PERF, LDSTAK)
 C
 C     PRINT TITLE AND CALL PPL TO DISPLAY SMOOTHED PERIODOGRAM
 C
       WRITE (IPRT,103)
       ILOG = 1
       CALL PPL (PERF, FREQ, NF, ILOG)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (10I5)
   101 FORMAT (10F7.2)
   102 FORMAT ('1RESULTS OF STARPAC PERIODOGRAM SUBROUTINE PGMS')
   103 FORMAT ('1RESULTS OF STARPAC FILTER SUBROUTINE MDFLT',
      *  ' DISPLAYED USING STARPAC PLOT SUBROUTINE PPL')
       END

 Data:

     3
     8    8    8
   261
    5.00  11.00  16.00  23.00  36.00  58.00  29.00  20.00  10.00   8.00
    3.00   0.00   0.00   2.00  11.00  27.00  47.00  63.00  60.00  39.00
   28.00  26.00  22.00  11.00  21.00  40.00  78.00 122.00 103.00  73.00
   47.00  35.00  11.00   5.00  16.00  34.00  70.00  81.00 111.00 101.00
   73.00  40.00  20.00  16.00   5.00  11.00  22.00  40.00  60.00  80.90
   83.40  47.70  47.80  30.70  12.20   9.60  10.20  32.40  47.60  54.00
   62.90  85.90  61.20  45.10  36.40  20.90  11.40  37.80  69.80 106.10
  100.80  81.60  66.50  34.80  30.60   7.00  19.80  92.50 154.40 125.90
   84.80  68.10  38.50  22.80  10.20  24.10  82.90 132.00 130.90 118.10
   89.90  66.60  60.00  46.90  41.00  21.30  16.00   6.40   4.10   6.80
   14.50  34.00  45.00  43.10  47.50  42.20  28.10  10.10   8.10   2.50
    0.00   1.40   5.00  12.20  13.90  35.40  45.80  41.10  30.10  23.90
   15.60   6.60   4.00   1.80   8.50  16.60  36.30  49.60  64.20  67.00
   70.90  47.80  27.50   8.50  13.20  56.90 121.50 138.30 103.20  85.70
   64.60  36.70  24.20  10.70  15.00  40.10  61.50  98.50 124.70  96.30
   66.60  64.50  54.10  39.00  20.60   6.70   4.30  22.70  54.80  93.80
   95.80  77.20  59.10  44.00  47.00  30.50  16.30   7.30  37.60  74.00
  139.00 111.20 101.60  66.20  44.70  17.00  11.30  12.40   3.40   6.00
   32.30  54.30  59.70  63.70  63.50  52.20  25.40  13.10   6.80   6.30
    7.10  35.60  73.00  85.10  78.00  64.00  41.80  26.20  26.70  12.10
    9.50   2.70   5.00  24.40  42.00  63.50  53.80  62.00  48.50  43.90
   18.60   5.70   3.60   1.40   9.60  47.40  57.10 103.90  80.60  63.60
   37.60  26.10  14.20   5.80  16.70  44.30  63.90  69.00  77.80  64.90
   35.70  21.20  11.10   5.70   8.70  36.10  79.70 114.40 109.60  88.80
   67.80  47.50  30.60  16.30   9.60  33.20  92.60 151.60 136.30 134.70
   83.90  69.40  31.50  13.90   4.40  38.00 141.70 190.20 184.80 159.00
  112.30









                                    <12-62>
1RESULTS OF STARPAC PERIODOGRAM SUBROUTINE PGMS
                                                                                                         STARPAC 2.08S (03/15/90)
 SAMPLE PERIODOGRAM (IN DECIBELS)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
        46.7879 -                   +                                                                                   -
                I                                                                                                       I
                I                                                                                                       I
                I                   2 +                                                                                 I
                I                     +                                                                                 I
        41.0117 -   2                +                                                                                  -
                I     +                                                                                                 I
                I     +            +                                                                                    I
                I                 ++      +                                                                             I
                I    +           +     2 +                                                                              I
        35.2356 -        +           ++  +                                                                              -
                I     ++ +      +         ++                                                                            I
                I  +          + + 2                                                                                     I
                I      +                +  ++         +                                                                 I
                I    +         + +       +  +        +                                                                  I
        29.4594 -         +2+                        +     +                                                            -
                I +            +        +    2      + ++                                                                I
                I       + + +++               ++    2    ++ +                                                           I
                I   +                        +     +    2       +                                                       I
                I  +    +                      +          ++  +          +                                              I
        23.6833 -      +                   +          +      +  ++  +                                                   -
                I              ++             +              2 +   +     +               +                              I
                I            +                  + 2               2   +                                         +       I
                I                              +++       +     ++ ++++  2   +    ++                    ++      +        I
                I            +         +          +                   +           ++   2 ++                +         +  I
        17.9071 -        +                      +             +  +     +   + +   +                 +      ++    +   ++  -
                I                                +     +             +     ++     +                +  + +         +     I
                I          +                               +                            2  +      + ++   +   2++ + ++   I
                I                                  +    ++  +       +              +  2    +                  +   ++  + I
                I                                                              +           +      + +  ++ +         +   I
        12.1310 -                                              +        +++  +2 +   + + ++   3+   +   +  +  +           -
                I                                                                +          + +             +    ++     I
                I                                                     +      + +     +         + 2        +           + I
                I                                                              +    +       +   +                +      I
                I                                                      +  +                      +  +          +        I
         6.3548 -                                                                                                       -
                I                                                               +                    ++                 I
                I                                                                              +           + +          I
                I                                                                    +          +                       I
                I                                                                                                       I
          .5787 -                                                                         +    +                        -
                I                                                                                                       I
                I                                                                   +                                   I
                I                                                                                                       I
                I                                                                                                       I
        -5.1975 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
       -10.9736 -                                                          +                                            -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-63>
1RESULTS OF STARPAC FILTER SUBROUTINE MDFLT DISPLAYED USING STARPAC PLOT SUBROUTINE PPL
                                                                                                         STARPAC 2.08S (03/15/90)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
      8545.2427 -                   322                                                                                 -
      8000.0000 -                 +2  +2                                                                                -
                I                 2    +                                                                                I
      6000.0000 -                +      2                                                                               -
                I                +       +                                                                              I
                I               2        +                                                                              I
      4000.0000 -               +        ++                                                                             -
                I 223          +          +                                                                             I
                I    22        +           +                                                                            I
                I     ++       +           +                                                                            I
                I      2                                                                                                I
                I       2     +            +                                                                            I
      2000.0000 -        2    +             +                                                                           -
                I        +   2              +                                                                           I
                I         2  +               +                                                                          I
                I          32                +                                                                          I
                I                            +                                                                          I
                I                                                                                                       I
      1000.0000 -                             +                                                                         -
                I                             +                                                                         I
       800.0000 -                              +                                                                        -
                I                              +                                                                        I
       600.0000 -                              +                                                                        -
                I                               2                                                                       I
                I                               +                                                                       I
       400.0000 -                                22+32323                                                               -
                I                                 ++     3+                                                             I
                I                                         +2                                                            I
                I                                          +2                                                           I
                I                                            3                                                          I
                I                                             2+                                                        I
       200.0000 -                                              2                                                        -
                I                                               3                                                       I
                I                                                22                                                     I
                I                                                 +2                                                    I
                I                                                   3                                                   I
                I                                                    2+                                                 I
       100.0000 -                                                     22                                                -
        80.0000 -                                                       3                                               -
                I                                                        3                                              I
                I                                                         2                                             I
        60.0000 -                                                          3                                            -
                I                                                           2                                           I
                I                                                            3+                                         I
        40.0000 -                                                             +32 +2323233                   2232332322 -
                I                                                                32       2+            32332+          I
                I                                                                          2+        +32                I
                I                                                                           ++      2+                  I
                I                                                                            2     2+                   I
                I                                                                             2+  3                     I
        20.0000 -                                                                              223                      -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000



                                                             <12-64>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE BFS USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y1 AND Y2 MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       REAL Y1(100), Y2(100)
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     READ NUMBER OF OBSERVATIONS
 C          SERIES 1
 C          SERIES 2
 C
       READ (5,100) N
       READ (5,101) (Y1(I), I=1,N)
       READ (5,101) (Y2(I), I=1,N)
 C
 C     PRINT TITLE AND CALL BFS FOR BIVARIATE FOURIER SPECTRUM ANALYSIS
 C
       WRITE (IPRT,102)
       CALL BFS (Y1, Y2, N)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (I5)
   101 FORMAT (12F6.2)
   102 FORMAT ('1RESULTS OF STARPAC',
      *  ' BIVARIATE FOURIER SPECTRUM ANALYSIS SUBROUTINE BFS')
       END





















                                    <12-65>
1Data:

   100
  -0.88 -0.16 -1.87 -1.12  1.38  2.13  2.76  0.56 -0.69 -1.79 -3.82 -2.38
   1.00  0.70 -0.15  0.98  0.11 -0.35 -0.73  0.89 -1.63 -0.44 -1.37 -1.71
  -1.22 -2.00 -0.22  0.38  1.31  0.71  0.32  0.48 -1.88 -0.94 -1.54 -0.13
   1.02  0.02 -0.77  0.11 -0.60 -0.52 -0.09  1.23  1.46  0.61  0.42  2.16
   3.18  2.10  0.37 -0.24  0.57 -0.53  2.44  1.02 -0.53 -2.49 -2.12 -1.04
  -0.12 -1.88 -1.50  1.54  3.33  3.08  1.71  0.79  1.55  0.89 -0.89 -1.18
   0.89  1.71  3.05  0.15 -1.04  0.12  0.08  0.11 -2.62 -1.28  1.07  3.20
   1.92  0.53 -1.08  0.49 -0.58  0.17  1.15 -0.97 -1.63  1.14 -0.67 -0.88
  -0.07  0.24  0.55 -2.16
   0.79  1.12 -1.10 -2.39 -1.75 -0.82 -0.36  1.27  1.75  2.44  0.36 -2.10
  -1.93 -1.30 -1.75 -0.34  0.74  0.49  0.70  0.71  0.09  0.59  1.54  0.14
   0.55 -1.40 -2.55 -1.66 -0.43  0.58  2.18 -0.24  0.58 -0.18 -1.55 -0.64
  -1.09  0.90 -0.66 -0.35  0.48  0.50  0.05 -0.68  0.24  0.58 -1.26 -0.25
   0.25  2.18  2.96  1.56 -0.36 -0.59 -0.12  3.03  2.11  0.78  0.89 -1.45
  -0.36 -0.37 -1.39 -4.19 -0.73 -0.98  0.36  0.06 -1.94 -0.08  0.17  1.00
  -0.05  0.43  0.15  2.69  0.57  0.29  1.10  0.48 -1.06 -2.28 -2.03 -0.75
   1.00  1.71  0.58  1.97  0.99  1.94  2.18  3.14  0.60  0.51  1.35  0.56
   0.11  0.00  2.34  1.88






































                                    <12-66>
1RESULTS OF STARPAC BIVARIATE FOURIER SPECTRUM ANALYSIS SUBROUTINE BFS
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (SQUARED COHERENCY COMPONENT) (+), 95 PCT. CONFIDENCE LIMITS (.) AND 95 PCT. SIGNIFICANCE LEVEL (-) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    4 / BW= .4657 / EDF=    93)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         1.0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .9000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .8000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .7000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .6000 -                                                                                                       -
                I                                                ............                                           I
                I                                             ...            ....                                       I
                I                                          ...                   ..                                     I
                I                                        ..                        ..                                   I
          .5000 -                                      ..                            ..                                 -
                I                                    ..                                ..                               I
                I                                  ..            ++++++++++++            .                              I
                I                                ..           +++            +++          .                             I
                I                               .          +++                  +++        ..                           I
          .4000 -                             ..         ++                        ++        .                          -
                I                            .         ++                            +        .                         I
                I                           .         +                               ++       ..                       I
                I                         ..        ++             ........             +        .                      I
                I                        .        ++           ....        ....          ++       .                     I
          .3000 -                      ..        +           ..                ...         +       .                    -
                I                     .        ++         ...                     ..        +       .                   I
                I                    .        +         ..                          .        ++      .                  I
                I                  ..       ++        ..                             ..        +      ..                I
                I                 .        +         .                                 ..       +       .               I
          .2000 -               ..       ++        ..                                    .       ++      .              -
                I              .        +        ..                                       ..       +      .             I
                I             .       ++       ..                                           .       +      .            I
                I                    +        .                                              ..      ++                 I
                I                  ++       ..                                                 .       +                I
          .1000 -                ++       ..                                                    ..      ++              -
                I              ++       ..                                                        ..      +             I
                I -----------22-------22------------------------------------------------------------22-----22---------- I
                I         +++      ...                                                                ..     +++        I
                I    +++++     ....                                                                     ...     +++     I
          .0000 - +++         .                                                                            .       ++++ -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-67>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (PHASE COMPONENT) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    4 / BW= .4657 / EDF=    93)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         6.2832 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         5.0265 -                                                                                                       -
                I                                                                                                       I
                I              +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          I
                I        ++++++                                                                               ++++      I
                I     +++                                                                                         ++    I
         3.7699 -    +                                                                                              +   -
                I   +                                                                                                +  I
                I ++                                                                                                  + I
                I                                                                                                       I
                I                                                                                                       I
         2.5133 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         1.2566 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -1.2566 -                                                                                                       -
                I                                                                                                       I
                I              +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          I
                I        ++++++                                                                               ++++      I
                I     +++                                                                                         ++    I
        -2.5133 -    +                                                                                              +   -
                I   +                                                                                                +  I
                I ++                                                                                                  + I
                I                                                                                                       I
                I                                                                                                       I
        -3.7699 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -5.0265 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -6.2832 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-68>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (SQUARED COHERENCY COMPONENT) (+), 95 PCT. CONFIDENCE LIMITS (.) AND 95 PCT. SIGNIFICANCE LEVEL (-) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    8 / BW= .2349 / EDF=    47)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         1.0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .9000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .8000 -                                                                                                       -
                I                            ..............                                                             I
                I                         ...              ...                                                          I
                I                       ..                    ...                                                       I
                I                      .                         ...                                                    I
          .7000 -                     .                             ..                                                  -
                I                    .         +++++++++              ..                                                I
                I                   .       +++         +++             .                                               I
                I                  .      ++               +++           ..                                             I
                I                 .      +                    ++           .                                            I
          .6000 -                .      +                       ++          .                                           -
                I                      +                          ++         ..                                         I
                I               .     +                             ++         .                                        I
                I              .     +           .....                +         .                                       I
                I                   +        ....     ...              +         .                                      I
          .5000 -             .             .            ...            ++        .                                     -
                I                  +      ..                ..            +        .                                    I
                I            .    +      .                    ..           +        .                                   I
                I                       .                       ..          +        .                                  I
                I           .    +     .                          .          +        .                                 I
          .4000 -               +     .                            ..         +        .                                -
                I          .                                         .         +        .                               I
                I         .    +     .                                .         +        ..                             I
                I                   .                                  .         +         .                            I
                I        .    +                                         ..        +         .                           I
          .3000 -                  .                                      .        +         .                          -
                I            +    .                                        .        +                                   I
                I                                                           .        +                                  I
                I           +    .                                           .        +                                 I
                I               .                                             .        +                                I
          .2000 -          +                                                   .        +                               -
                I         +    .                                                .        ++                             I
                I                                                                .         +                            I
                I        +    .                                                   ..        +                           I
                I -----------2------------------------------------------------------2--------22------------------------ I
          .1000 -       +                                                            .         ++                       -
                I      +    .                                                         ..         ++                     I
                I     +    .                                                            .          +++                  I
                I    +    .                                                              ...          +++               I
                I   +    .                                                                  ..           ++++++         I
          .0000 - ++                                                                                           ++++++++ -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-69>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (PHASE COMPONENT) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=    8 / BW= .2349 / EDF=    47)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         6.2832 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         5.0265 -                                                                                                       -
                I                                                    +++++++++++++++++                                  I
                I     +++++++++++++++++++++++++++++++++++++++++++++++                 +++++++++++                       I
                I   ++                                                                           ++++++++               I
                I                                                                                        ++++++         I
         3.7699 -  +                                                                                           ++++     -
                I                                                                                                  ++   I
                I +                                                                                                  ++ I
                I                                                                                                       I
                I                                                                                                       I
         2.5133 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         1.2566 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -1.2566 -                                                                                                       -
                I                                                    +++++++++++++++++                                  I
                I     +++++++++++++++++++++++++++++++++++++++++++++++                 +++++++++++                       I
                I   ++                                                                           ++++++++               I
                I                                                                                        ++++++         I
        -2.5133 -  +                                                                                           ++++     -
                I                                                                                                  ++   I
                I +                                                                                                  ++ I
                I                                                                                                       I
                I                                                                                                       I
        -3.7699 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -5.0265 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -6.2832 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-70>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (SQUARED COHERENCY COMPONENT) (+), 95 PCT. CONFIDENCE LIMITS (.) AND 95 PCT. SIGNIFICANCE LEVEL (-) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   16 / BW= .1191 / EDF=    24)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         1.0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
          .9000 -                                                                                                       -
                I                      ..                                                                               I
                I                  ....  .....                                                                          I
                I                ..           ..                                                                        I
                I               .               ...                                                                     I
          .8000 -              .                   ...............                                                      -
                I             .                                   ..                                                    I
                I            .       ++++++                         .                                                   I
                I                  ++      ++                        .                                                  I
                I           .     +          ++                       .                                                 I
          .7000 -                +             ++                      .                                                -
                I          .    +                +                      .                                               I
                I                                 +++     ++++++         .                                              I
                I         .    +                     +++++      ++        .                                             I
                I             +                                   +        ..                                           I
          .6000 -                                                  +         .                                          -
                I        .           ......                         +         ..                                        I
                I            +      .      ..                                   .                                       I
                I                  .         .                       +           .                                      I
                I       .   +     .           .                       +           ..                                    I
          .5000 -                .             .                       +            .                                   -
                I                               .                                                                       I
                I          +    .                ..                     +                                               I
                I                                  .       .....         +                                              I
                I              .                    .......     .         +                                             I
          .4000 -         +                                      .         +                                            -
                I             .                                   .         +                                           I
                I                                                  .         +                                          I
                I        +                                          .         +                                         I
                I            .                                                 +                                        I
          .3000 -                                                    .          ++                                      -
                I                                                     .           +                                     I
                I       +   .                                                      +                                    I
                I -----------------------------------------------------2------------22--------------------------------- I
                I          .                                            .             ++                                I
          .2000 -      +                                                                +++++                           -
                I                                                        .                   ++++                       I
                I         .                                               .                      ++                     I
                I     +                                                    ..                      +                    I
                I                                                            .                      +                   I
          .1000 -        .                                                    ..                     +                  -
                I    +                                                          .                     ++                I
                I       .                                                        ..                     +               I
                I   +                                                              ..                    ++          ++ I
                I  +                                                                                       +++   ++++   I
          .0000 - +                                                                                           +++       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-71>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (PHASE COMPONENT) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   16 / BW= .1191 / EDF=    24)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         6.2832 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
         5.0265 -  +                                                         ++++++++                                   -
                I   ++++                         +++++++              +++++++        +++++                              I
                I       +++++++++++++++++++++++++       ++++++++++++++                    +++                           I
                I                                                                            +++                        I
                I                                                                               +++                     I
         3.7699 -                                                                                  +++                  -
                I                                                                                     ++                I
                I                                                                                       ++              I
                I                                                                                         +             I
                I                                                                                          ++           I
         2.5133 -                                                                                            +          -
                I                                                                                                       I
                I                                                                                             +         I
                I                                                                                                       I
                I                                                                                              +        I
         1.2566 -                                                                                               +       -
                I                                                                                                       I
                I                                                                                                +      I
                I                                                                                                 ++    I
                I                                                                                                   +   I
          .0000 - 2                                                                                                  +2 -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
        -1.2566 -  +                                                         ++++++++                                   -
                I   ++++                         +++++++              +++++++        +++++                              I
                I       +++++++++++++++++++++++++       ++++++++++++++                    +++                           I
                I                                                                            +++                        I
                I                                                                               +++                     I
        -2.5133 -                                                                                  +++                  -
                I                                                                                     ++                I
                I                                                                                       ++              I
                I                                                                                         +             I
                I                                                                                          ++           I
        -3.7699 -                                                                                            +          -
                I                                                                                                       I
                I                                                                                             +         I
                I                                                                                                       I
                I                                                                                              +        I
        -5.0265 -                                                                                               +       -
                I                                                                                                       I
                I                                                                                                +      I
                I                                                                                                 ++    I
                I                                                                                                   +   I
        -6.2832 -                                                                                                    +  -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-72>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (SQUARED COHERENCY COMPONENT) (+), 95 PCT. CONFIDENCE LIMITS (.) AND 95 PCT. SIGNIFICANCE LEVEL (-) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   33 / BW= .0595 / EDF=    12)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         1.0000 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I                                                                                                       I
                I                    ........              ....                                                         I
          .9000 -                  ..        ..           .    .                                                        -
                I           .......            .         .                                                              I
                I          .                    .               .                                                       I
                I                                       .        .                                                      I
                I         .           +++++      .          ++                                                          I
          .8000 -                    +     ++     ..   .   +      .                                                     -
                I        .          +        +      ...       +    ...                                                  I
                I                  +                      +           ..                                                I
                I                 +           +                +        .                                               I
                I            +++++             +                                                                        I
          .7000 -           +                            +                                                              -
                I                                               +                                                       I
                I          +                    +                                                                       I
                I                                       +                                                               I
                I                                +               +                                                      I
          .6000 -         +            ...                                                                              -
                I                     .   .       +    +    ..    +                                                     I
                I                    .     .               .                                                            I
                I                           .      +  +       .    +                                                    I
                I                   .        .      ++              +                                                   I
          .5000 -        +                                           ++                                                 -
                I                  .                      .            +                                                I
                I ----------------------------2----------------2--------2---------------------------------------------- I
                I                 .                                                                                     I
                I            ... .                                       +        +                                     I
          .4000 -       +       .              .         .                       +                                      -
                I           .                                                      +            +                       I
                I                                               .         +     +              + +                      I
                I          .                    .                                   +         +   +                     I
                I                                                          +   +                                        I
          .3000 -                                       .                                    +                          -
                I      +                         .               .          +++      +             +                    I
                I                                                                                                       I
                I         .                                                                 +       +                   I
                I                                 .    .          .                   +                                 I
          .2000 -                                  .                                       +         +                  -
                I     +                             . .            .                   +                                I
                I        .                           .              ..                    +           +                 I
                I                                                     ..                ++             ++          ++++ I
                I                                                                                        ++             I
          .1000 -    +                                                  .                                  +      +     -
                I                                                                                                +      I
                I   +                                                                                       +           I
                I                                                                                            +  +       I
                I ++                                                                                          ++        I
          .0000 -                                                                                                       -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-73>
1
                                                                                                         STARPAC 2.08S (03/15/90)
 -- SMOOTHED FOURIER SPECTRUM (PHASE COMPONENT) --
    (PARZEN WINDOW WITH LAG WIND. TRUNC. PT.=   33 / BW= .0595 / EDF=    12)
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
         6.2832 -                                                                                                       -
                I                                                                                                       I
                I                                                                                                       I
                I  +                                                                                                    I
                I   +                                                              ++++                                 I
         5.0265 -    +                             ++                   +++++++++++    +                                -
                I     ++++                      +++  +++             +++                ++                              I
                I         ++++++++++++++++++++++        +++++++++++++                     +                             I
                I                                                                          +++++                        I
                I                                                                               +++                     I
         3.7699 -                                                                                  ++                   -
                I                                                                                    +                  I
                I                                                                                     +                 I
                I                                                                                      +                I
                I                                                                                       +++             I
         2.5133 -                                                                                          ++           -
                I                                                                                            +          I
                I                                                                                                       I
                I                                                                                             +         I
                I                                                                                                       I
         1.2566 -                                                                                                       -
                I                                                                                              +        I
                I                                                                                                       I
                I                                                                                               ++      I
                I                                                                                                 +++   I
          .0000 - 2                                                                                                  +2 -
                I                                                                                                       I
                I                                                                                                       I
                I  +                                                                                                    I
                I   +                                                              ++++                                 I
        -1.2566 -    +                             ++                   +++++++++++    +                                -
                I     ++++                      +++  +++             +++                ++                              I
                I         ++++++++++++++++++++++        +++++++++++++                     +                             I
                I                                                                          +++++                        I
                I                                                                               +++                     I
        -2.5133 -                                                                                  ++                   -
                I                                                                                    +                  I
                I                                                                                     +                 I
                I                                                                                      +                I
                I                                                                                       +++             I
        -3.7699 -                                                                                          ++           -
                I                                                                                            +          I
                I                                                                                                       I
                I                                                                                             +         I
                I                                                                                                       I
        -5.0265 -                                                                                                       -
                I                                                                                              +        I
                I                                                                                                       I
                I                                                                                               ++      I
                I                                                                                                 +++   I
        -6.2832 -                                                                                                    +  -
                 -I---------I---------I---------I---------I---------I---------I---------I---------I---------I---------I-
                  .0000     .0500     .1000     .1500     .2000     .2500     .3000     .3500     .4000     .4500     .5000
+FREQ
 PERIOD         INF       20.       10.        6.6667    5.        4.        3.3333    2.8571    2.5       2.2222    2.
                                                             <12-74>
1G.  Acknowledgments

      The code for computing the autocovariance and cross  covariance functions
 using a  fast Fourier  transform  are based  on subroutines  written  by Jones
 [1971].   The subroutines  which compute the fast Fourier  transform and which
 compute the  Fourier transform of  real data  are those  written  by Singleton
 [1969] and  the subroutine  to  compute the  cosine transform  was  written by
 Jones.

      The transformation used  to compute  the univariate  Fourier  spectrum is
 performed using the algorithm shown on  page 311 of Jenkins and  Watts [1968].
 The algorithm used is described on pages 418-420.

      The code for computing the autoregressive order selection  statistics and
 autoregressive spectrum estimates  is based  on subroutines  written  by Jones
 [1971].   The coefficients PHI(k), k = 1, ..., IAR of the autoregressive model
 used  by  the  autoregressive  spectrum  subroutines  are  computed  from  the
 autocovariance function using the Levinson-Durbin recursive method for solving
 the  Yule-Walker equations  discussed  in Appendix  A3.2 of  Box  and  Jenkins
 [1976].

      The subroutines for the split-cosine-bell taper and the  modified Daniell
 filter operation were adapted from subroutines TAPER and MODDAN given on pages
 116 and 178 of Bloomfield [1976].



































                                    <12-75>
1-----                             CHAPTER 13                             -----

                                 ARIMA MODELING


 A.  Introduction

      STARPAC contains five user-callable subroutines for  AutoRegressive 
 Integrated Moving  Average (ARIMA) modeling.   Three are  for computing  the least
 squares  estimates of  the parameters  of  an  ARIMA  model and  two  are  for
 computing the minimum mean square error forecasts using these estimates.  Both
 the estimation and forecasting subroutines allow several levels of  control of
 the computations  and printed output.   The estimation  subroutines also allow
 the user to specify a subset of the parameters to be treated as constants with
 their values held  fixed at their input values.   This last feature allows the
 user  to examine  the  results  obtained  estimating various  subsets  of  the
 parameters of a general model without respecifying the model for each subset.

      Each of  the subroutines  in  this chapter  models the input  series, yi,
 i = 1, ..., N, with a user-specified general multiplicative ARIMA  model using
 the techniques discussed  in Box and Jenkins [1976].   Briefly, this  model is
 defined by

         NFAC                           NFAC
      [PRODUCT PHI(p(j),B[s(j)])] * ([PRODUCT BDO(s(j),d(j))]*Y(i) - mu) =
         j=1                            j=1

                          NFAC
                       [PRODUCT  THETA(q(j),B[s(j)]]*a(i)
                          j=1

 for i = 1, ..., N, where

 N is the number of observations in the series;

 NFAC is the number of factors in the model;

 B[s] is the backward shift operator of order s, i.e.,

      B[s]*Y(i) = Y(i-s);

 mu is the expected value of the differenced series, i.e.,

               NFAC
   mu = E ( [PRODUCT BDO(s(j),d(j))]*yi),
               j=1

   which can be used to allow for a deterministic polynomial trend;

 PHI(p(j),B[s(j)]) is the polynomial in B[s(j)] of order p(j), i.e.,

   PHI(p(j),B[s(j)]) = 1 - phi(1,j)*B[s(j)] -

                       phi(2,j)B[2*s(j)] - ... phi(p(j),j)*B[p(j)*s(j)],

   which represents the jth autoregressive factor in the model;

 BDO(s(j),d(j)) is the backward difference operator, i.e.,

                                    <13-1>
1
   BDO(s(j),d(j)) = (1 - B[s(j)])**d(j)

                  = 1 - d(j)*B[s(j)] + ... B[d(j)*s(j)],

   which represents the jth difference factor in the model;

 THETA(q(j),B[s(j)]) is the polynomial in B[s(j)] of order q(j), i.e.,

   THETA(q(j),B[s(j)]) = 1 - theta(1,j)*B[s(j)] -

                         theta(2,j)*B[2*s(j)] - ... theta(q(j),j)*B[q(j)*s(j)],

   which represents the jth moving average factor in the model; and

 a(i) is the unobservable random noise component at the ith observation.

      The least squares estimates of the parameters,

                     phi(1,1), phi(2,1), ..., phi(p(1),1),
                   phi(1,2), phi(2,2), ..., phi(p(2),2), ...,
                phi(1,NFAC), phi(2,NFAC), ..., phi(p(NFAC),NFAC)
                                      mu,
                  theta(1,1), theta(2,1), ..., theta(q(1),1),
                theta(1,2), theta(2,2), ..., theta(q(2),2), ...,
             theta(1,NFAC), theta(2,NFAC), ..., theta(q(NFAC),NFAC)

 are obtained using back forecasts at each iteration as  discussed  in  section
 E.1.a.  The  least  squares  solution is that which minimizes (with respect to
 the parameters) the sum of the squares of the  random  noise  components,  ai,
 i.e.,

                                   N
                                  SUM     a(i)**2
                              i=-infinity

 The iterative procedure used is documented in chapter 9.

      The user must supply both initial values for the parameters and  an array
 specifying the orders p(j), d(j),  q(j) and  s(j) for each  factor j = 1, ...,
 NFAC  in the model.   Initial parameter  values for the estimation subroutines
 should be chosen with care since good values can significantly reduce computer
 time.

      Users  are  directed  to  section  B  for  a  brief  description  of  the
 subroutines.  The declaration and CALL statements are given in section  C  and
 the  subroutine  arguments  are defined in section D.  The algorithms used and
 output produced by these  subroutines  are  discussed  in  section  E.  Sample
 programs and their output are shown in section F.


 B.  Subroutine Descriptions

 B.1  ARIMA Estimation Subroutines

      The   simplest   of   the  three  ARIMA  estimation  subroutines,   AIME,
 automatically summarizes the estimated results and a variety of statistics  in
 a  five-part  printed  report  described  in  section  E.2.a,  and returns the

                                    <13-2>
1estimated parameters and residuals to the user  via  the  subroutine  argument
 list  (level one control).  Most ARIMA estimation problems can be solved using
 AIME.

      The other two estimation subroutines,  AIMEC and AIMES,  provide  greater
 flexibility to the user at the price of more input.

      AIMEC,  like  AIME,  also returns estimated parameters and residuals from
 the fit.  In addition, it allows the user to supply arguments to indicate
               - a subset of the  model parameters to be treated  as constants,
                 with their values held fixed at their input values;
               - the step sizes used to compute the numerical approximations to
                 the derivative as described in chapter 9, section E.1.b;
               - the maximum number of iterations allowed;
               - the convergence criteria;
               - the scale (i.e., the typical size) of each parameter;
               - the maximum  change  allowed in  the parameters at  the  first
                 iteration;
               - how the variance-covariance matrix is to be approximated; and
               - what sections of the five-part printed report are wanted.

      AIMES has all the  features  of  AIMEC  and,  in  addition,  returns  the
 following estimated values via the argument list:
               - the number of parameters actually estimated;
               - the residual standard deviation;
               - the predicted values;
               - the standard deviations of the predicted values;
               - the standardized residuals; and
               - the variance-covariance matrix of the estimated parameters.


 B.2  ARIMA Forecasting Subroutines

      The  simplest   of  the   three  ARIMA  forecasting   subroutines,  AIMF,
 automatically summarizes the  estimated results in a two-part  printed report.
 Forecasts are made using N as the origin  and extending [N/10] + 1  steps into
 the   future,    i.e.,   observations   are   forecast   for    indices   N+1,
 N+2, ..., N+[N/10]+1.  Many forecasting problems can be solved using AIMF.

      The second  forecasting  subroutine, AIMFS,  allows the  user  to  supply
 arguments to  indicate the number of forecasts  to be made, the origins  to be
 used  and  the amount  of printed  output.  This subroutine  also returns  the
 forecasts and their standard deviations via the argument list.


 C.  Subroutine Declaration and CALL Statements

 NOTE:  Argument definitions and sample programs are given in sections D and F,
 respectively.  The  conventions  used to present the following declaration and
 CALL statements are given in chapter 1, sections B and D.

      The <basic declaration block> identifies declaration statements  that are
 needed by all of the  ARIMA estimation and forecasting subroutines.   The user
 should substitute the following four statements for each occurrence  of <basic
 declaration block> given below.




                                    <13-3>
1                        <real> Y(n), PAR(npar)
                         INTEGER MSPEC(4,nfac)
                         DOUBLE PRECISION DSTAK(ldstak)
                         COMMON /CSTAK/ DSTAK

                                      ===


                        Subroutines for ARIMA Estimation

 AIME:   Compute and print  a five-part least squares analysis of the parameter
         estimates of an ARIMA model; return parameter estimates and residuals

         <basic declaration block>
         <real> RES(n)
         :
         :
         CALL AIME (Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK)

                                      ===

 AIMEC:  Compute  and  optionally print  a five-part least  squares analysis of
         the parameter estimates of an ARIMA model using  user-supplied control
         values; return parameter estimates and residuals

         <basic declaration block>
         <real> RES(n)
         INTEGER IFIXED(npar)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         :
         :
         CALL AIMEC  (Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK,
        +             IFIXED, STP, MIT, STOPSS, STOPP,
        +             SCALE, DELTA, IVAPRX, NPRT)

                                      ===

 AIMES:  Compute and  optionally  print a  five-part least squares  analysis of
         the parameter estimates of an ARIMA model using  user-supplied control
         values;  return parameter  estimates, residuals, number  of parameters
         estimated, residual  standard deviation,  predicted  values,  standard
         deviations of  the predicted values and variance-covariance  matrix of
         the estimated parameters

         <basic declaration block>
         <real> RES(n)
         INTEGER IFIXED(npar)
         <real> STP(npar), STOPSS, STOPP, SCALE(npar), DELTA
         <real> RSD, PV(n), SDPV(n), SDRES(n), VCV(npare,npare)
         :
         :
         CALL AIMES (Y, N, M, MSPEC, NFAR, PAR, NPAR, RES, LDSTAK,
        +            IFIXED, STP, MIT, STOPSS, STOPP,
        +            SCALE, DELTA, IVAPRX, NPRT,
        +            NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)

                                      ===


                                    <13-4>
1
                       Subroutines for ARIMA Forecasting

 AIMF:   Compute and  print  the minimum  mean square error  forecasts obtained
         using an ARIMA model

         <basic declaration block>
         :
         :
         CALL AIMF (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK)

                                     ===

 AIMFS:  Compute  and  optionally print the minimum mean square error forecasts
         obtained using an  ARIMA model;  return forecasts  and  their standard
         errors

         <basic declaration block>
         <real> FCST(nfcst,nfcsto), SDFCST(nfcst)
         INTEGER IFCSTO(nfesto)
         :
         :
         CALL AIMFS (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK,
        +            NFCST, NFCSTO, IFCSTO, NPRT, FCST, IFCST, FCSTSD)

                                     ===


 D.  Dictionary of Subroutine Arguments and COMMON Variables

 NOTE:   --> indicates  that the argument is input to  the subroutine  and that
             the input value is preserved;
         <-- indicates that the argument is returned by the subroutine;
         <-> indicates that the argument  is input  to the subroutine  and that
             the input value is overwritten by the subroutine;
         --- indicates that the argument  is input  to some subroutines  and is
             returned by others;
         *** indicates that the argument is a subroutine name;
         ... indicates that the variable is passed via COMMON.


 DELTA   --> The maximum scaled change  allowed in the parameters at  the first
             iteration.  [See section E.1.a.] The default value is 100.0.  When
             DELTA <= 0.0 or when DELTA is not an argument  of  the  subroutine
             CALL  statement  the  default  value  is used.  A smaller value of
             DELTA  may  be  appropriate  if,  at  the  first  iteration,   the
             computation  of the predicted values from the user's model 
	     subroutine produces an arithmetric overflow or the parameters leave  the
             region  of  interest in parameter space.  A reasonable alternative
             to the default value of DELTA is an  upper  bound  to  the  scaled
             change  that the estimated parameters should be allowed to make on
             the first iteration,

             DELTA = min(|del(PAR(k))|/SCALE(k), for k = 1, ..., NPAR)

             where  del(PAR(k))  is  the  maximum  change  allowed  for  the  kth
             parameter at the first iteration.


                                    <13-5>
1DSTAK   ... The  DOUBLE  PRECISION vector  in COMMON /CSTAK/  of  dimension at
             least LDSTAK.  DSTAK provides workspace for the computations.  The
             first  LDSTAK  locations  of  DSTAK  will  be  overwritten  during
             subroutine execution.

 FCST    <-- The array of dimension at least NFCST by NFCSTO that  contains the
             NFCST forecasts computed from each of the NFCSTO origins.

 FCSTSD  <-- The vector of dimension at least NFCST that contains the standard
             deviation of each of the forecasts,

                                     K-1
             FCSTSD(K) = RSD * sqrt[ SUM C(j)**2] for K = 1, ..., NFCST,
                                     j=0

             where

             C(0) = 1 and

             C(j) = c(1)*C(j-1) + ... + c(P+D)*C(j-P-D) - theta(j)

                      NFAC           NFAC
             with P = SUM  p(j); D = SUM  d(j); and
                      j=1             j=1

                  c(i) = the coefficient of Bi in the polynomial defined by

                           NFAC
                         PRODUCT [PHI(p(j),Bs[(j)]) * BDO(s(j),d(j))]
                           j=1

 IERR    ... An  error  flag  returned  in  COMMON /ERRCHK/.   [See  chapter 1,
             section D.5.]  Note that using (or not using) the error flag  will
             not  affect  the  printed  error  messages  that are automatically
             provided.

             For AIME, AIMEC, and AIMES:

               IERR = 0 indicates that  no errors  were detected  and  that the
                        iterations converged satisfactorily.

               IERR = 1 indicates that improper input was detected.

               IERR = 2 indicates that  the computation of the residual  sum of
                        squares using the initial parameter values  produced an
                        arithmetic overflow.   The user  should reduce the size
                        of DELTA or should supply new starting values.

               IERR = 3 indicates that  the model is  computationally singular,
                        which means the model has too many parameters  near the
                        solution.   The user should examine the model and  data
                        to identify and  remove the cause of the singularity.

               IERR = 4 indicates that at least one of the standardized 
	                residuals   could  not  be  computed  because   its  standard
                        deviation was  zero.   The validity  of the   
			variance-covariance matrix is questionable.


                                    <13-6>
1              IERR = 5 indicates false convergence  [see  chapter  9,  section
                        E.1.a].

               IERR = 6 indicates  that  convergence  was not  reached  in  the
                        allowed  number of iterations or model subroutine calls
                        [see argument MIT].

               IERR = 7 indicates that the variance-covariance matrix could not
                        be computed because of computational difficulties.

             For AIMF and AIMFS:

               IERR = 0 indicates that no errors were detected and that all the
                        forecasts were computed.

               IERR = 1 indicates that improper input was detected.

 IFCST   --> The exact  value of  the  first dimension  of the  matrix  FSCT as
             specified in the calling subroutine.

 IFCSTO  --> The vector of dimension  at least NFCSTO that contains  the NFCSTO
             indices to be used as origins,  where 1 <= IFCSTO(J) <= N for J  =
             1,  ...,  NFCSTO.  The default value for each element of IFCSTO is
             N.  When IFCSTO(J) is outside the range [1, N] or IFCSTO is not an
             argument of the subroutine CALL statement  the  default  value  is
             used.

 IFIXED  --> The vector of dimension at least NPAR that contains values used to
             indicate whether  the  corresponding parameter  in PAR  is  to  be
             treated  as  a  fixed  constant   or  is  to  be  estimated.    If
             IFIXED(I) > 0, PAR(I) will be  held fixed  at its input  value; if
             IFIXED(I) = 0, PAR(I)  will be  estimated using the  least squares
             procedure  described  in  section  A.   The  default  values   are
             IFIXED(I)  =  0,  I  =  1,  ...,  NPAR,  i.e.,  all parameters are
             estimated.  When IFIXED(1)  <=  -1,  or  when  IFIXED  is  not  an
             argument of the subroutine CALL statement,  the default value will
             be used.

 IVAPRX  --> The indicator variable used to specify how the variance-covariance
             matrix,  VCV, is  to be  approximated.   Three approximations  are
             available:

             (1) VCV = RSD**2 * inv(trans(Dhat)*W*Dhat)

             (2) VCV = RSD**2 * inv(Hhat)

             (3) VCV = RSD**2 * inv(Hhat) * (trans(Dhat)*W*Dhat) * inv(Hhat)

             where

             trans(.) indicates the transpose of the designated matrix;

             inv(.) indicates the inverse of the designated matrix;

             Hhat is the matrix of second partial derivatives of the model with
               respect to each parameter (the Hessian matrix), evaluated at the
               solution.


                                    <13-7>
1              = trans(Dhat)*W*Dhat +

                   N
                 (SUM  e(i)*wt(i)*(second partial of e(i) wrt PAR(j) & PAR(k)
                  i=1
                                   for j = 1, ..., NPAR & k = 1, ..., NPAR));

             W is an N by N diagonal matrix of weights,

               W = diag(wt(i), i = 1, ..., N),

             when a weighted analysis is performed, and  is the identity matrix
             otherwise,  and

             Dhat is the matrix that contains the partial  derivatives  of  the
               model with respect  to  each  parameter  (the  Jacobian matrix),
               evaluated at the solution.

             Approximation  (1)  is  based  on  the  assumption   that   H   is
             approximately  equal  to  trans(D)*W*D  because  the residuals are
             sufficiently small at the solution;  approximation (2) is based on
             the  assumption  that  the  necessary  conditions  for  asymptotic
             maximum likelihood theory have been met;  and approximation (3) is
             based   on  the  assumption  that  the  necessary  conditions  for
             asymptotic maximum likelihood theory may be violated.  The results
             of  a  study  by  Donaldson  and  Schnabel  [1987]  indicate  that
             approximation  (1)  is  preferable  because  it  is  simple,  less
             expensive,  more numerically stable and at least  as  accurate  as
             approximations  (2)  and (3).  However,  all approximations to the
             variance-covariance  matrix  are  subject  to  sampling  variation
             because  they  are  computed using the estimated parameter values.
             The  variance-covariance  matrix  computed  for   any   particular
             nonlinear least squares solution should thus be regarded as only a
             rough estimate [Bard, 1974; Donaldson and Schnabel, 1987].

             If IVAPRX = 1 or 4 then approximation (1) is used;
                       = 2 or 5 then approximation (2) is used; and
                       = 3 or 6 then approximation (3) is used.

             If IVAPRX = 1, 2, or 3,   then,   when    user-supplied   analytic
             derivatives are available [see argument NLSDRV], they are  used to
             compute  VCV; if  IVAPRX = 4,  5, or  6, then  only  the predicted
             values from the model subroutine  are used  to compute VCV.   When
             analytic derivatives  are  available, options  1, 2,  or  3,  will
             generally result in a faster, more accurate computation of VCV.

             The  default  value for  IVAPRX is 1.   When  argument  IVAPRX  is
             outside the range [1, 6], or when IVAPRX is not an argument of the
             subroutine CALL statement, then the default value will be used.

 IVCV    --> The exact  value of  the  first dimension  of the  matrix  VCV  as
             specified in the calling program.







                                    <13-8>
1LDSTAK  --> The length of the DOUBLE PRECISION workspace vector DSTAK.  LDSTAK
             must equal or exceed  the appropriate value given below,  where if
             the single precision  version of  STARPAC is  being  used P = 0.5,
             otherwise P = 1.0 [see chapter 1, section B].  Also,

                             NFAC                   NFAC
             MBO = maximum [ SUM  (p(j)+d(j))*s(j), SUM  q(j)*s(j) ].
                             j=1                    j=1

             For AIME, AIMEC and AIMES:

               LDSTAK >= 43 + max[IS*(N+NPAR), 30+NPARE] + 2*NFAC +

                         max [ IS*10*N + 6*MBO+606 ,
                               94 + 4*(N+MBO+101) + 5*MBO +
                                 (3*NPARE**2+35*NPARE)/2 ]*P

               where  IS = 1 if default values are used for the derivative step
               sizes, and IS = 0 otherwise.

             For AIMF and AIMFS:

               LDSTAK >= 18 + 2*NFAC + [5*MBO+2*(N+MBO+101)]*P

 MIT     --> The  maximum number of iterations allowed.   This argument is also
             used to  compute  the maximum  number of  model  subroutine calls,
             (2*MIT).   The iterations  will stop  if either limit  is reached,
             although,  as a rule,  the maximum  number of  iterations  will be
             reached  first.   The  default  value for  the maximum  number  of
             iterations is 21.  When MIT <= 0 or when MIT is not an argument of
             the subroutine CALL statement the default value will be used.

 MSPEC   --> The array  of dimension exactly  4 rows  by at least  NFAC columns
             that contains  the orders  p, d, q  and s  for each factor  in the
             model where p(j), j = 1, ..., NFAC must be in row 1,
                         d(j), j = 1, ..., NFAC must be in row 2,
                         q(j), j = 1, ..., NFAC must be in row 3, and
                         s(j), j = 1, ..., NFAC must be in row 4.
             Values of p(j), d(j), q(j) and  s(j), j = 1, ..., NFAC,  must each
             equal or exceed zero.

 N       --> The number of observations.

 NFAC    --> The number of factors in the model.

 NFCST   --> The number  of forecasts  to be  computed.   The default  value is
             [N/10]+1.   When NFCSTF0  or when NFCST is not  an argument of the
             subroutine CALL statement the default value is used.

 NFCSTO  --> The number of forecast origins supplied.   The default value is 1.
             When NFCSTOF0 or when NFCSTO is not an argument of  the subroutine
             CALL statement the default value is used.







                                    <13-9>
1NPAR    --> The number of parameters  in the  model including both  those held
             fixed  at  their  starting  values  and  those  which  are  to  be
             estimated,

                        NFAC
             NPAR = 1 + SUM  [p(j) + q(j)] .
                        j=1

 NPARE   <-- The number of  parameters actually estimated, i.e., the  number of
             zero elements  in IFIXED.   N.B.   This value  is returned  by the
             estimation subroutines.

 NPRT    --> The argument controlling printed output.

             For AIME, AIMEC and AIMES:

               NPRT is  a five-digit integer,  in which  the value  of  the Ith
               digit (counting from left to  right) is used to control  the Ith
               section of the output.

               If the Ith digit = 0 the output from the Ith section is
                                    suppressed;
                                = 1 the brief form of the Ith section is given;
                                >=2 the full form of the Ith section is given.

               The  default value for NPRT is 11112.  When NPRT <= -1,  or when
               NPRT is not an argument in the subroutine  CALL  statement,  the
               default value will be used.  If the convergence criteria are not
               satisfied the subroutine gives a suitable warning and provides a
               printed  report  even  if  NPRT  =  0.  A full discussion of the
               printed output is given in section E.2.a and  is  summarized  as
               follows.

               Section 1 lists  the  starting  estimates  and  control  values.
                         Brief output and  full output  are the  same  for this
                         section.

               Section 2 reports  the results of the iterations.   Brief output
                         includes  information only  about the  first  and last
                         iteration while full output includes information about
                         all of the iterations.

               Section 3 provides information for each observation based on the
                         final solution.  Brief output includes information for
                         the  first 40 observations while full  output provides
                         the information for all of the data.

               Section 4 is a  set of  four residual plots.   Brief output  and
                         full output are the same for this section.

               Section 5 is  the final  summary of  the  estimated  parameters.
                         Brief  output  does not include printing the 
			 variance-covariance matrix while full output does.






                                    <13-10>
1            For AIMF and AIMFS:

               If NPRT = 0 the printed output is suppressed.

               If NPRT <> 0 the printed output is provided.

               The default value for NPRT is 1.  When NPRT <> -1 or  when  NPRT
               is  not an argument in the subroutine CALL statement the default
               value will be used.

 PAR     --- The vector of dimension at least NPAR that contains  the parameter
             values.   For both the estimation and the forecasting subroutines,
             parameter values must be ordered

             phi(1,1), phi(2,1), ..., phi(p(1),1),
             phi(1,2), phi(2,2), ..., phi(p(2),2), ...,
             phi(1,NFAC), phi(2,NFAC), ..., phi(p(NFAC),2)
             mu,
             theta(1,1), theta(2,1), ..., theta(q(1),1),
             theta(1,2), theta(2,2), ..., theta(q(2),2), ...,
             theta(1,NFAC), theta(2,NFAC), ..., theta(q(NFAC),2)

             i.e., the  parameter values  from the  autoregressive  factors are
             first,  followed by mu,  followed by the parameter values from the
             moving  average  factors.  For all estimation subroutines PAR must
             contain initial values  for  the  parameters  on  input  and  will
             contain  the  final values on return.  For the forecasting 
	     subroutines  PAR  must  contain  the  parameter  values  for  which  the
             forecasts are to be computed.

 PV      <-- The  vector of dimension  at least  N that contains  the predicted
             values of the dependent variable at the solution,

             PV(i) = Y(i) - a(i)     for i = 1, .., N.

 RES     <-- The vector of dimension at least N that contains the  residuals at
             the solution,

             RES(i) = a(i)    for i = 1, ..., N.

 RSD     <-- The residual standard deviation at the solution,

                          N
             RSD = sqrt[ SUM RES(i)**2/(N-NPARDF-NPARE) ]
                         i=1

                            NFAC
             where NPARDF = SUM  s(j)*d(j) .
                            j=1

 SCALE   --> The vector of dimension at least NPAR that contains the  scale, or
             typical size,  of each  parameter.   The vector  SCALE is  used to
             normalize the size of each parameter so that at each iteration,

             |PAR(j)/SCALE(j)| approximates |PAR(k)/SCALE(k)|

             for j = 1, ..., NPAR and k = 1, ..., NPAR.


                                    <13-11>
1            Values of |SCALE(k)| > |PAR(k)| can be used to increase  the  step
             size  in cases where the model function is known to be insensitive
             to small changes  in  the  value  PAR(k),  although  normally  the
             default values should be used.

             The  default values for SCALE are selected by the NL2SOL algorithm
             [Dennis et al.,  1981a,b] and are updated at each iteration.  When
             SCALE is not an argument in the subroutine CALL statement or  when
             the  user-supplied  value  for SCALE(1) <= 0 the default procedure
             will be used to select scale values.  When SCALE(1) > 0, values of
             SCALE(k) <= 0 for k = 2, ..., NPAR will be interpreted as an input
             error.  User-supplied scale values may be either a vector  of  the
             typical  size of each parameter or a vector of ones if the typical
             sizes of the parameters are  roughly  equal;  user-supplied  scale
             values can sometimes result in reduced com-puting time since these
             values are not updated at each iteration.

 SDPV    <-- The vector of dimension at least N that contains  an approximation
             to  the  standard  deviation  of  each  predicted  value   at  the
             solution,

             SDPV(i) = the ith diagonal element of sqrt(Dhat*VCV*trans(Dhat))

             for i = 1, ..., N, where

             Dhat(i,j) = partial [ a(i) wrt PAR(j) ]

             for i = 1, ..., N and j = 1, ..., NPAR, evaluated at the solution,
             and trans(Dhat) is the transpose of Dhat.

             This approximation is based on a linearization of the model in the
             neighborhood of  the solution;  the validity of  the approximation
             depends on the nonlinearity of the model.   This approximation may
             be extremely  inaccurate  for a  problem with  a  highly nonlinear
             model.

 SDRES   <-- The vector of dimension at least N that contains  an approximation
             to the standardized residuals at the solution,

             SDRES(i) = RES(i) / sqrt[RSD**2-SDPV(i)**2]

             for i = 1, ..., N,  which  is  the  ith residual  divided  by  its
             individual  estimated standard deviation.   This approximation  is
             based on a linearization of  the model in the neighborhood  of the
             solution;  the  validity  of  the  approximation  depends  on  the
             nonlinearity  of the model.   This approximation  may be extremely
             inaccurate for a problem with a highly nonlinear model.

 STOPP   --> The stopping value for  the convergence test based on  the maximum
             scaled relative  change  in  the  parameters at  the  most  recent
             iteration.   The convergence criterion is satisfied if the current
             step is a Newton step and

               max[|PARc(k)-PARp(k)|/SCALE(k) for k = 1, ..., NPAR]
             --------------------------------------------------------  < STOPP.
             max[(|PARc(k)|+|PARp(k)|)/SCALE(k) for k = 1, ..., NPAR]

             where PARc(k) and PARp(k) indicate the current value and the value

                                    <13-12>
1            from the previous iteration,  respectively,  of the kth  parameter
             [see  Dennis  et  al.  1981a].  This  convergence  test is roughly
             equivalent to the test based on the  maximum  relative  change  in
             each parameter as measured by

             max(|PARc(k)-PARp(k)|/|PARp(k)| for k = 1,  ..., NPAR).

             STOPP  is  not  a scale-dependent value;  if its value is 10**(-4)
             then this criteria will be met when the first four digits of  each
             parameter  are the same at two successive iterations regardless of
             the size of the parameter values.

             The default value is approximately 10**(-DIGITS/2),  where  DIGITS
             is the number of decimal digits carried by the user's computer for
             a  single  precision  value  when  the single precision version of
             STARPAC is being used and is  the  number  carried  for  a  double
             precision value otherwise.  When the user-supplied value for STOPP
             is  outside  the  interval  [0.0,  1.0]  or  when  STOPP is not an
             argument of the subroutine CALL statement the default  value  will
             be used.

 STOPSS  --> The stopping value for the convergence test based on the  ratio of
             the   forecasted   change   in   the   residual  sum  of  squares,
             fcst(RSS(PAR)),  to  the residual sum of squares from the previous
             iteration, RSS(PARp), where

                             N
             RSS(PAR) =     SUM     a(i)**2
                        i=-infinity

             and PARp are the parameter values from the previous iteration.

             The convergence criterion is satisfied if certain  conditions  are
             met and

             fcst(RSS(PAR))/RSS(PARp) < STOPSS

             [see  Dennis  et  al.  1981a].  This  convergence  test is roughly
             equivalent to the  test  based  on  the  relative  change  in  the
             residual  standard  deviation  between  the  current  and previous
             iterations as measured by (RSDc -  RSDl)/RSDc.  STOPSS  is  not  a
             scale-dependent value; if its value is 10**(-5) this criteria will
             be  met  when the first five digits of the residual sum of squares
             are the same at two successive iterations regardless of  the  size
             of the residual sum of squares.

             The default value is approximately the maximum  of  10**(-10)  and
             10**(-2*DIGITS/3),  where DIGITS is the number of  decimal  digits
             carried  by  the user's computer for a single precision value when
             the single precision version of STARPAC is being used and  is  the
             number  carried  for a double precision value otherwise.  When the
             user-supplied   value   for   STOPSS   is   outside  the  interval
             [10**(-DIGITS),  0.1] or when STOPSS is not  an  argument  of  the
             subroutine CALL statement the default value will be used.

 STP     --- The vector of dimension  at least NPAR that contains  the relative
             step sizes used by the estimation subroutines to  approximate  the
             derivative  matrix  numerically.  The procedure used to select the

                                    <13-13>
1            default values is described in chapter 9, section E.1.a.  When STP
             is not an argument of the subroutine CALL statement or when STP(1)
             <= 0 the default values will be used for all of  the  step  sizes;
             when STP(1) > 0,  values of STP(k) <= 0 for k = 2,  ..., NPAR will
             be interpreted as an input error.

 VCV     <-- The matrix of dimension at least NPARE by NPARE that  contains the
             variance-covariance matrix  of the estimated  parameters, 
	     approximated as designated by argument IVAPRX.   The parameters which are
             held fixed (as specified  by argument IFIXED) are not  included in
             the variance-covariance matrix.

             The approximation of the variance-covariance matrix is based  on a
             linearization of the  model in  the neighborhood of  the solution;
             the validity of  the approximation depends on the  nonlinearity of
             the model.   This approximation may be extremely inaccurate for  a
             problem with a highly nonlinear model.

 Y       --> The vector of dimension at least N that contains the  series being
             modeled.


 E.  Computational Methods

 E.1  Algorithms

 E.1.a  ARIMA Estimation

      The ARIMA estimation subroutines use the NL2SOL software package  written
 by Dennis et al.  [1981a,b].   The  observations  of  the  series,  which  are
 measured with error,  are iteratively fit to the ARIMA model by minimizing the
 sums of squares of the  estimated  random  noise  component  as  described  in
 section  A.  The  back forecasting technique discussed on pages 215-220 of Box
 and Jenkins [1976] is used to compute the random noise component.  Up  to  101
 back  forecasts are computed.  The back forecasts are assumed to be negligible
 when their  magnitude  is  less  than  0.01  times  the  first  value  of  the
 differenced  series  (centered  about  its  mean)  obtained  entirely from the
 observed data.  If,  at the last iteration,  the 101st back  forecast  is  not
 negligible a warning message is printed.

     The iterations continue until the convergence criteria based on the change
 in  the  parameter  values  or  in  the  residual sum of squares are satisfied
 [arguments STOPP and STOPSS],  the maximum  number  of  iterations  (or  model
 subroutine calls) is reached [argument MIT],  or the iterations are terminated
 due to singularity in the model or false convergence.  All but  the  first  of
 these stopping conditions may indicate computational problems and will produce
 an error report [see chapter 1, section D.5].  Singular convergence means that
 the  model  contains  too many parameters,  at least near the solution,  while
 false convergence can indicate that either STOPSS or STOPP is  set  too  small
 for  the accuracy to which the model and its derivatives are being computed or
 that there is a discontinuity in  the  derivative.  Iterative  procedures  for
 solving  nonlinear least squares problems are discussed in Dennis and Schnabel
 [1983], Draper and Smith [1981],  and Kennedy and Gentle [1980].  The specific
 procedure  used  in  STARPAC  is  discussed  in  chapter  9 and  Dennis et al.
 [1981a].




                                    <13-14>
1E.1.b ARIMA Forecasting

      The ARIMA forecasting subroutines use the techniques discussed in Box and
 Jenkins [1976],  chapter 5.  The back forecasting technique discussed on pages
 215-220  of  Box  and  Jenkins  [1976]  is  used  to  compute the random noise
 component needed for the forecasts.  Values of a(i) for  i  greater  than  the
 forecast  origin  are  assumed  to  be  zero.  Up  to  101  back forecasts are
 computed.  The  back  forecasts  are  assumed  to  be  negligible  when  their
 magnitude  is  less  than 0.01 times the first value of the differenced series
 (centered about its mean) obtained entirely from the  observed  data.  If  the
 101st  back  forecast  is  not  negligible a warning message is printed.


 E.2  Computed Results and Printed Output

 E.2.a  The ARIMA Estimation Subroutines

      The argument controlling  the  printed  output,  NPRT,  is  discussed  in
 section D.

      The  output  from  the  ARIMA  estimation  subroutines  consists  of five
 sections,  several of which include tables summarizing  the  results.  In  the
 following descriptions,  the actual table headings are given by the  uppercase
 phrases  enclosed  in angle braces (<...>).  Results which correspond to input
 or returned subroutine CALL statement arguments are identified by the argument
 name in uppercase (not enclosed in angle braces).


 Section 1 provides a summary of the initial estimates and control values.   It
           lists the following information.

       * The initial values of the parameters, PAR, and whether they are  to be
         held fixed or not as indicated by argument IFIXED.

       * The scale values, SCALE.

       * The step sizes used to approximate the derivatives numerically, STP.

       * The number of observations, N.

       * The maximum number of iterations allowed, MIT.

       * The maximum number of model subroutine calls allowed.

       * The two convergence criteria, STOPSS and STOPP.

       * The maximum change in  the parameters allowed at the  first iteration,
         DELTA.

       * The residual  sum  of squares  computed using  the  starting parameter
         values.

       * The residual standard deviation computed using the  starting parameter
         values, RSD.





                                    <13-15>
1Section 2 lists selected  information about  each iteration  and  includes the
           reason the iterations were terminated.  The information provided for
           each iteration includes the following.

       * The iteration number.

       * <MODEL  CALLS>:   the total number of times since execution began that
         the user's model subroutine  has  been  called,  not  including  calls
         required to approximate the derivatives numerically.

       * <RSD>:   the  residual standard deviation computed using the parameter
         values from the current iteration.

       * <RSS>: the residual sum of squares computed using the parameter values
         from the current iteration.

       * <REL CHNG RSS>:  the relative change in the residual  sum  of  squares
         caused by the current iteration.

       * <FORECASTED  REL  CHNG  RSS>:   the  forecasted relative change in the
         residual sum of squares at the current  iteration,  and  whether  this
         value was checked against STOPSS (<CHKD> = Y) or not (<CHKD> = N).

       * <REL CHNG PAR>:   the maximum scaled relative change in the parameters
         at the current iteration,  and whether this value was checked  against
         STOPP (<CHKD> = Y) or not (<CHKD> = N).

       * <CURRENT PARAMETER VALUES>:   the estimated parameter values resulting
         from the current iteration.


 Section 3 provides the following information for each observation, i = 1, ...,
           N, based on the final solution.

       * <ROW>:  the row number of the observations.

       * <SERIES>:  the value of the dependent variable, Y.

       * <PREDICTED VALUE>:  the estimated predicted value, PV, from the fit.

       * <STD DEV OF PRED VALUE>:  the standard   deviation  of  the  predicted
         value, SDPV.

       * <RESIDUAL>:  the error estimate, RES.

       * <STD RES>:  the standardized residual, SDRES.


 Section 4 displays the following plots.

       * The standardized residuals versus row numbers.

       * The autocorrelation function of the (non-standardized) residuals.

       * The normal probability plot of the standardized residuals.




                                    <13-16>
1Section 5 summarizes  the  following  information about  the  final  parameter
           estimates and their variances.

       * The  variance-covariance  matrix,  VCV,  of  the  estimated  (unfixed)
         parameters and the corresponding correlation matrix,

         rjk = VCV(j,k)/sqrt[VCV(j,j)*VCV(k,k)]     for j = 1, ..., NPARE
                                                    and k = 1, ..., NPARE.

       * <PARAMETER ESTIMATES  (PAR)>:  the  final  value  of  each  parameter,
         PAR(k), k = 1, ..., NPAR.

       * <STD  DEV  OF  PARAMETER  ESTIMATES>:   the standard deviation of each
         estimated parameter,

         sqrt[VCV(k,k)]     for k = 1, ..., NPAR.

       * <RATIO PAR/SD OF PAR>:  the ratio of each  estimated parameter to  its
         standard deviation,

         RATIO(k) = PAR(k) / sqrt[VCV(k,k)]     for k = 1, ..., NPAR.

       * <APPROXIMATE 95-PERCENTCONFIDENCE LIMITS>:  the lower  and  upper  95-
         percent  confidence  limits  for  each  parameter,  computed using the
         appropriate value of the Student's t distribution with

                NFAC
         N -  [ SUM  s(j)*d(j)] - NPARE
                i=1

         degrees of freedom.

       * The residual sum of squares at the solution.

       * The residual standard deviation at the solution, RSD.

                                                 NFAC
       * The residual degrees of freedom,  N - [ SUM  s(j)*d(j)] - NPARE.
                                                 i=1

       * An approximation to the condition number of the derivative matrix,
           D(i,j) = partial [a(i) wrt PAR(j)]

         for i = 1, ..., N and j = 1, ..., NPAR, evaluated at the solution,

         (the Jacobian), under the  assumption that the absolute error  in each
         column  of D is roughly equal.   The approximation will be meaningless
         if  this  assumption   is  not   valid;  otherwise  it   will  usually
         underestimate the actual condition number by a factor of from 2  to 10
         [see Dongarra et al., 1979, page 9.5].


 NOTE:  The standard  deviation  of  the  predicted  values,  the  standardized
 residuals,  the  variance-covariance  matrix,  the  standard deviations of the
 parameters and the 95-percent confidence limits  on  the  parameters  are  all
 based  on  a  linear  approximation  to  the  model  in  a neighborhood of the
 solution;  the validity of this approximation depends on the  nonlinearity  of


                                    <13-17>
1the  model.  The  statistics  based  on  this  approximation  may be extremely
 inaccurate for a problem with a highly  nonlinear  model.


 E.2.b The ARIMA Forecasting Subroutines

      The  argument  controlling  the  printed  output,  NPRT,  is discussed in
 section D.

      The output from the ARIMA forecasting subroutines consists of  a  summary
 of the model used to produce the forecasts and,  for each origin, a plot and a
 list of the computed forecasts and a 95-percent confidence interval about  the
 forecasts along with the actual series value when known.


 F.  Examples

      The  sample programs  of this section  use the  model and  data  given in
 table 9.1 of Box and Jenkins [1976]; the model is

 BDO(1,1)*BDO(12,1)*yi - mu = (1-theta(1,1)*B[1]) * (1-theta(1,2)*B[12]) * a(i)

 for i = 1, ..., N.


      ARIMA  Estimation.  In  the  first sample program below,  AIME is used to
 compute the least squares estimates of the parameters.


      ARIMA Forecasting.  In the second sample program, AIMF is used to compute
 the minimum mean square error forecasts  using  the  least  squares  estimates
 obtained in example 1.



























                                    <13-18>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE AIME USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y, PAR AND RES MUST BE CHANGED TO DOUBLE
 C          PRECISION IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       INTEGER MSPEC(4,5)
       REAL Y(200), PAR(5), RES(200)
       DOUBLE PRECISION DSTAK(5000)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 5000
 C
 C     READ NUMBER OF FACTORS IN MODEL
 C          VALUES OF P, D, Q AND S FOR EACH FACTOR
 C          NUMBER OF PARAMETERS
 C          STARTING VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) NFAC
       READ (5,100) ((MSPEC(I,J), I=1,4), J=1,NFAC)
       READ (5,100) NPAR, N
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,101) (Y(I), I=1,N)
 C
 C     COMPUTE LOG OF DATA
 C
       DO 10 I = 1, N
         Y(I) = ALOG(Y(I))
    10 CONTINUE
 C
 C     PRINT TITLE AND CALL AIME TO PERFORM ARIMA ESTIMATION ANALYSIS
 C
       WRITE (IPRT,102)
       CALL AIME (Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (4I5)
   101 FORMAT (12F6.1)
   102 FORMAT ('1RESULTS OF STARPAC ARIMA ESTIMATION SUBROUTINE AIME')
       END



                                    <13-19>
1Data:

     2
     0    1    1    1
     0    1    1   12
     3  144
    0.0   0.4   0.6
  112.0 118.0 132.0 129.0 121.0 135.0 148.0 148.0 136.0 119.0 104.0 118.0
  115.0 126.0 141.0 135.0 125.0 149.0 170.0 170.0 158.0 133.0 114.0 140.0
  145.0 150.0 178.0 163.0 172.0 178.0 199.0 199.0 184.0 162.0 146.0 166.0
  171.0 180.0 193.0 181.0 183.0 218.0 230.0 242.0 209.0 191.0 172.0 194.0
  196.0 196.0 236.0 235.0 229.0 243.0 264.0 272.0 237.0 211.0 180.0 201.0
  204.0 188.0 235.0 227.0 234.0 264.0 302.0 293.0 259.0 229.0 203.0 229.0
  242.0 233.0 267.0 269.0 270.0 315.0 364.0 347.0 312.0 274.0 237.0 278.0
  284.0 277.0 317.0 313.0 318.0 374.0 413.0 405.0 355.0 306.0 271.0 306.0
  315.0 301.0 356.0 348.0 355.0 422.0 465.0 467.0 404.0 347.0 305.0 336.0
  340.0 318.0 362.0 348.0 363.0 435.0 491.0 505.0 404.0 359.0 310.0 337.0
  360.0 342.0 406.0 396.0 420.0 472.0 548.0 559.0 463.0 407.0 362.0 405.0
  417.0 391.0 419.0 461.0 472.0 535.0 622.0 606.0 508.0 461.0 390.0 432.0








































                                    <13-20>
1RESULTS OF STARPAC ARIMA ESTIMATION SUBROUTINE AIME
                                                                                                         STARPAC 2.08S (03/15/90)
+*****************************************************************************
 *  NONLINEAR LEAST SQUARES ESTIMATION FOR THE PARAMETERS OF AN ARIMA MODEL  *
 *                             USING BACKFORECASTS                           *
 *****************************************************************************


 SUMMARY OF INITIAL CONDITIONS
 ------------------------------


    MODEL SPECIFICATION

       FACTOR          (P     D     Q)    S

            1           0     1     1     1
            2           0     1     1    12



                                                                           --STEP SIZE FOR
                                         ------PARAMETER                   --APPROXIMATING
 -----------------PARAMETER DESCRIPTION  STARTING VALUES  ----------SCALE  -----DERIVATIVE
 INDEX  ---------TYPE  --ORDER  --FIXED  ----------(PAR)  --------(SCALE)  ----------(STP)

     1             MU      ---       NO    .00000000E+00
+                                                                DEFAULT    .46415888E-09
     2  MA (FACTOR 1)        1       NO    .40000000E+00
+                                                                DEFAULT    .37387871E-08
     3  MA (FACTOR 2)       12       NO    .60000000E+00
+                                                                DEFAULT    .31405989E-08


 NUMBER OF OBSERVATIONS                                                (N)   144

 MAXIMUM NUMBER OF ITERATIONS ALLOWED                                (MIT)    21

 MAXIMUM NUMBER OF MODEL SUBROUTINE CALLS ALLOWED                             42

 CONVERGENCE CRITERION FOR TEST BASED ON THE

      FORECASTED RELATIVE CHANGE IN RESIDUAL SUM OF SQUARES       (STOPSS)   .3696E-09
      MAXIMUM SCALED RELATIVE CHANGE IN THE PARAMETERS             (STOPP)   .8425E-07


 MAXIMUM CHANGE ALLOWED IN THE PARAMETERS AT THE FIRST ITERATION   (DELTA)   100.0

 RESIDUAL SUM OF SQUARES FOR INPUT PARAMETER VALUES                          .1762      (BACKFORECASTS INCLUDED)

 RESIDUAL STANDARD DEVIATION FOR INPUT PARAMETER VALUES              (RSD)   .3710E-01

 BASED ON DEGREES OF FREEDOM  144 -  13 -   3 =  128







                                                             <13-21>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED


 ITERATION NUMBER    1
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         2       .3709E-01       .1761           .5189E-03   .4714E-03   Y   .1195E-01   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2              3
          VALUE  -.1333249E-03   .3979228       .6145143


 ITERATION NUMBER    5
 ----------------------
     MODEL                                                     FORECASTED
     CALLS         RSD             RSS        REL CHNG RSS    REL CHNG RSS    REL CHNG PAR
                                                              VALUE   CHKD    VALUE   CHKD
         6       .3709E-01       .1761           .8255E-11   .1213E-10   Y   .1690E-05   Y

      CURRENT PARAMETER VALUES
          INDEX    1              2              3
          VALUE  -.1405415E-03   .3954493       .6162809

 ***** RESIDUAL SUM OF SQUARES CONVERGENCE *****































                                                             <13-22>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED


 RESULTS FROM LEAST SQUARES FIT
 -------------------------------

                       -----PREDICTED  ----STD DEV OF                  ---STD
  ROW  --------SERIES  ---------VALUE  ----PRED VALUE  ------RESIDUAL  ---RES

    1   .47184989E+01   .47186578E+01   .23375651E-02  -.15897742E-03     .00
    2   .47706846E+01   .47648934E+01   .61625875E-02   .57912016E-02     .16
    3   .48828019E+01   .48947824E+01   .20921259E-02  -.11980469E-01    -.32
    4   .48598124E+01   .48486540E+01   .24469160E-02   .11158358E-01     .30
    5   .47957905E+01   .48222257E+01   .78395571E-02  -.26435108E-01    -.73
    6   .49052748E+01   .49267056E+01   .45016136E-02  -.21430783E-01    -.58
    7   .49972123E+01   .50175551E+01   .39518870E-02  -.20342841E-01    -.55
    8   .49972123E+01   .50095167E+01   .36735960E-02  -.12304444E-01    -.33
    9   .49126549E+01   .49084972E+01   .53741438E-02   .41577016E-02     .11
   10   .47791235E+01   .47746722E+01   .24312596E-02   .44513296E-02     .12
   11   .46443909E+01   .46444172E+01   .24998555E-02  -.26258075E-04     .00
   12   .47706846E+01   .47883601E+01   .28293434E-02  -.17675496E-01    -.48
   13   .47449321E+01   .47855791E+01   .40174359E-02  -.40647013E-01   -1.10
   14   .48362819E+01   .48094434E+01   .62145599E-02   .26838485E-01     .73
   15   .49487599E+01   .49464401E+01   .14990013E-02   .23197915E-02     .06
   16   .49052748E+01   .49149161E+01   .50785271E-03  -.96412749E-02    -.26
   17   .48283137E+01   .48639358E+01   .64208047E-02  -.35622108E-01    -.98
   18   .50039463E+01   .49585091E+01   .46579865E-02   .45437200E-01    1.23
   19   .51357984E+01   .50850892E+01   .37490159E-02   .50709232E-01    1.37
   20   .51357984E+01   .51182302E+01   .54067261E-02   .17568189E-01     .48
   21   .50625950E+01   .50385922E+01   .43893513E-02   .24002855E-01     .65
   22   .48903491E+01   .49177012E+01   .33191572E-02  -.27352053E-01    -.74
   23   .47361984E+01   .47673933E+01   .10008824E-02  -.31194901E-01    -.84
   24   .49416424E+01   .48855743E+01   .29056588E-02   .56068117E-01    1.52
   25   .49767337E+01   .49143196E+01   .47883251E-02   .62414135E-01    1.70
   26   .50106353E+01   .50168153E+01   .75564688E-02  -.61800171E-02    -.17
   27   .51817836E+01   .51305277E+01   .24220978E-02   .51255824E-01    1.38
   28   .50937502E+01   .51243959E+01   .52120786E-02  -.30645701E-01    -.83
   29   .51474945E+01   .50483710E+01   .66029905E-02   .99123467E-01    2.72
   30   .51817836E+01   .52471047E+01   .80246117E-02  -.65321181E-01   -1.80
   31   .52933048E+01   .53191486E+01   .33433028E-02  -.25843802E-01    -.70
   32   .52933048E+01   .53049155E+01   .33631955E-02  -.11610671E-01    -.31
   33   .52149358E+01   .52140413E+01   .44372459E-02   .89444210E-03     .02
   34   .50875963E+01   .50649018E+01   .24239797E-02   .22694498E-01     .61
   35   .49836066E+01   .49368895E+01   .33405697E-02   .46717120E-01    1.26
   36   .51119878E+01   .51282796E+01   .58403723E-02  -.16291861E-01    -.44
   37   .51416636E+01   .51285808E+01   .24975828E-02   .13082783E-01     .35
   38   .51929569E+01   .51892704E+01   .31852286E-02   .36864202E-02     .10
   39   .52626902E+01   .53294127E+01   .39440625E-02  -.66722480E-01   -1.81
   40   .51984970E+01   .52322795E+01   .60354329E-02  -.33782432E-01    -.92
    .               .               .               .               .       .
    .               .               .               .               .       .
    .               .               .               .               .       .
  144   .60684256E+01   .60533746E+01   .64437477E-02   .15051031E-01     .41





                                                             <13-23>
1                                                                                                        STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED

                                                    STD RES VS ROW NUMBER
  3.75++---------+---------+---------+---------+---------+----+----+---------+---------+---------+---------+---------++
      -                                                                                                               -
      -                                                                                                               -
      -                                                                                                               -
      -                      *                                                                                        -
  2.25+                                                                                                      *        +
      -                                *      *                                                                       -
      -                  *                                *                                                           -
      -              *   * *                 *          *         **                                                  -
      -             *            *      *                *    *            *                   *      *               -
   .75+                                   *            *         *                   *        **   * *            *   +
      -          *    *         *                      *               *              *                 ** *  * *   * -
      - **                         *  *        *             *   *         *      * *        *        *              *-
      -*     ***  *             *  *             *            *      *  * * ** * * * *              *     *           -
      -  *  *      *      *    *  *        **      *       **  *     *          *      * *       *       *     * * *  -
  -.75+   ***  *       *      *        *             *              * *  *    *   *     *         *         *      *  +
      -            *    *   *        *           **        *                   *         ** *          *        *     -
      -         *                        *          *           *                          **      *                  -
      -                                      *                                                                        -
      -                      *      *           *                                                                     -
 -2.25+                                                                                         *                     +
      -                                                                                                      *        -
      -                                                                                                               -
      -                                                                                                               -
      -                                               *                                                               -
 -3.75++---------+---------+---------+---------+---------+----+----+---------+---------+---------+---------+---------++
      1.0                                                   72.5                                                  144.0

             AUTOCORRELATION FUNCTION OF RESIDUALS                        NORMAL PROBABILITY PLOT OF STD RES
     1++---------+---------+----*----+---------+---------++   3.75++---------+---------+----+----+---------+---------++
      -                         *                         -       -                                                   -
      -                       ***                         -       -                                                   -
      -                      ****                         -       -                                                   -
      -                         *                         -       -                                                  *-
     6+                         **                        +   2.25+                                               *   +
      -                         *                         -       -                                            **     -
      -                       ***                         -       -                                          **       -
      -                         ***                       -       -                                      ****         -
      -                         *                         -       -                                   ****            -
    11+                         *                         +    .75+                                ****               +
      -                         *                         -       -                              ***                  -
      -                         **                        -       -                           ***                     -
      -                         *                         -       -                         ***                       -
      -                         *                         -       -                     ****                          -
    16+                      ****                         +   -.75+                 *****                             +
      -                        **                         -       -             *****                                 -
      -                         *                         -       -         *****                                     -
      -                         *                         -       -         *                                         -
      -                       ***                         -       -      ***                                          -
    21+                         *                         +  -2.25+     *                                             +
      -                         *                         -       -   *                                               -
      -                         ****                      -       -                                                   -
      -                         *                         -       -                                                   -
      -                         *                         -       -*                                                  -
    26++---------+---------+----*----+---------+---------++  -3.75++---------+---------+----+----+---------+---------++
    -1.00                      0.0                     1.00     -2.5                       0.0                      2.5
                                                             <13-24>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+NONLINEAR LEAST SQUARES ESTIMATION FOR THE PARAMETERS OF AN ARIMA MODEL, CONTINUED


    MODEL SPECIFICATION

       FACTOR          (P     D     Q)    S

            1           0     1     1     1
            2           0     1     1    12



 VARIANCE-COVARIANCE AND CORRELATION MATRICES OF THE ESTIMATED (UNFIXED) PARAMETERS
 ----------------------------------------------------------------------------------

    - APPROXIMATION BASED ON ASSUMPTION THAT RESIDUALS ARE SMALL
    - COVARIANCES ARE ABOVE THE DIAGONAL
    - VARIANCES ARE ON THE DIAGONAL
    - CORRELATION COEFFICIENTS ARE BELOW THE DIAGONAL

 COLUMN         1                2                3

      1      .8126291E-06    -.3384623E-06    -.1984874E-06
      2     -.4601934E-02     .6656526E-02    -.3227447E-03
      3     -.3143769E-02    -.5648059E-01     .4905375E-02



 ESTIMATES FROM LEAST SQUARES FIT
 ---------------------------------

                                         ------PARAMETER  -----STD DEV OF                   ---------------------APPROXIMATE
 -----------------PARAMETER DESCRIPTION  ------ESTIMATES  ------PARAMETER  ----------RATIO  ----95 PERCENT CONFIDENCE LIMITS
 INDEX  ---------TYPE  --ORDER  --FIXED  ----------(PAR)  ------ESTIMATES  PAR/(SD OF PAR)  ----------LOWER  ----------UPPER

     1             MU      ---       NO   -.14054146E-03
+                                                          .90145941E-03   -.15590437E+00   -.16331169E-02    .13520340E-02
     2  MA (FACTOR 1)        1       NO    .39544931E+00
+                                                          .81587537E-01    .48469328E+01    .26036219E+00    .53053643E+00
     3  MA (FACTOR 2)       12       NO    .61628090E+00
+                                                          .70038379E-01    .87991885E+01    .50031609E+00    .73224571E+00


 NUMBER OF OBSERVATIONS                                                (N)   144


 RESIDUAL SUM OF SQUARES                  .1761309      (BACKFORECASTS INCLUDED)

 RESIDUAL STANDARD DEVIATION              .3709478E-01
 BASED ON DEGREES OF FREEDOM  144 -  13 -   3 =  128

 APPROXIMATE CONDITION NUMBER             90.50754






                                                             <13-25>
1Program:

       PROGRAM EXAMPL
 C
 C     DEMONSTRATE AIMF USING SINGLE PRECISION VERSION OF STARPAC
 C
 C     N.B. DECLARATION OF Y AND PAR MUST BE CHANGED TO DOUBLE PRECISION
 C          IF DOUBLE PRECISION VERSION OF STARPAC IS USED.
 C
       INTEGER MSPEC(4,5)
       REAL Y(200), PAR(5)
       DOUBLE PRECISION DSTAK(5000)
 C
       COMMON /CSTAK/ DSTAK
 C
 C     SET UP INPUT AND OUTPUT FILES
 C     [CHAPTER 1, SECTION D.4, DESCRIBES HOW TO CHANGE OUTPUT UNIT.]
 C
       CALL IPRINT(IPRT)
       OPEN (UNIT=IPRT, FILE='FILENM')
       OPEN (UNIT=5, FILE='DATA')
 C
 C     SPECIFY NECESSARY DIMENSIONS
 C
       LDSTAK = 5000
 C
 C     READ NUMBER OF FACTORS IN MODEL
 C          VALUES OF P, D, Q AND S FOR EACH FACTOR
 C          NUMBER OF PARAMETERS
 C          VALUES FOR PARAMETERS
 C          NUMBER OF OBSERVATIONS
 C          OBSERVED SERIES
 C
       READ (5,100) NFAC
       READ (5,100) ((MSPEC(I,J), I=1,4), J=1,NFAC)
       READ (5,100) NPAR, N
       READ (5,101) (PAR(I), I=1,NPAR)
       READ (5,102) (Y(I), I=1,N)
 C
 C     COMPUTE LOG OF DATA
 C
       DO 10 I = 1, N
         Y(I) = ALOG(Y(I))
    10 CONTINUE
 C
 C     PRINT TITLE AND CALL AIMF TO PERFORM ARIMA FORECASTING ANALYSIS
 C
       WRITE (IPRT,103)
       CALL AIMF (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK)
 C
 C     FORMAT STATEMENTS
 C
   100 FORMAT (4I5)
   101 FORMAT (12F6.1)
   102 FORMAT (12F6.3)
   103 FORMAT ('1RESULTS OF STARPAC ARIMA FORECASTING SUBROUTINE AIME')
       END


                                    <13-26>
1Data:

     2
     0    1    1    1
     0    1    1   12
     3  144
  0.000 0.395 0.615
  112.0 118.0 132.0 129.0 121.0 135.0 148.0 148.0 136.0 119.0 104.0 118.0
  115.0 126.0 141.0 135.0 125.0 149.0 170.0 170.0 158.0 133.0 114.0 140.0
  145.0 150.0 178.0 163.0 172.0 178.0 199.0 199.0 184.0 162.0 146.0 166.0
  171.0 180.0 193.0 181.0 183.0 218.0 230.0 242.0 209.0 191.0 172.0 194.0
  196.0 196.0 236.0 235.0 229.0 243.0 264.0 272.0 237.0 211.0 180.0 201.0
  204.0 188.0 235.0 227.0 234.0 264.0 302.0 293.0 259.0 229.0 203.0 229.0
  242.0 233.0 267.0 269.0 270.0 315.0 364.0 347.0 312.0 274.0 237.0 278.0
  284.0 277.0 317.0 313.0 318.0 374.0 413.0 405.0 355.0 306.0 271.0 306.0
  315.0 301.0 356.0 348.0 355.0 422.0 465.0 467.0 404.0 347.0 305.0 336.0
  340.0 318.0 362.0 348.0 363.0 435.0 491.0 505.0 404.0 359.0 310.0 337.0
  360.0 342.0 406.0 396.0 420.0 472.0 548.0 559.0 463.0 407.0 362.0 405.0
  417.0 391.0 419.0 461.0 472.0 535.0 622.0 606.0 508.0 461.0 390.0 432.0








































                                    <13-27>
1RESULTS OF STARPAC ARIMA FORECASTING SUBROUTINE AIME
                                                                                                         STARPAC 2.08S (03/15/90)
+***********************
 *  ARIMA FORECASTING  *
 ***********************


 MODEL SUMMARY
 -------------


    MODEL SPECIFICATION

       FACTOR          (P     D     Q)    S

            1           0     1     1     1
            2           0     1     1    12

                                ------PARAMETER
 --------PARAMETER DESCRIPTION  ------ESTIMATES
 INDEX  ---------TYPE  --ORDER  ----------(PAR)

     1             MU      ---    .00000000E+00
     2  MA (FACTOR 1)        1    .39500000E+00
     3  MA (FACTOR 2)       12    .61500000E+00


 NUMBER OF OBSERVATIONS                                                (N)   144


 RESIDUAL SUM OF SQUARES                  .1761643      (BACKFORECASTS INCLUDED)

 RESIDUAL STANDARD DEVIATION              .3709830E-01
 BASED ON DEGREES OF FREEDOM  144 -  13 -   3 =  128


























                                                             <13-28>
1
                                                                                                         STARPAC 2.08S (03/15/90)
+ARIMA FORECASTING, CONTINUED


 FORECASTS FOR ORIGIN  1


                                                                                   --------------------95  PERCENT
   5.8659397           6.1720075           6.4780753                               --------------CONFIDENCE LIMITS ---------ACTUAL
             6.0189736           6.4780753           6.6311092     ------FORECASTS ----------LOWER ----------UPPER -------IF KNOWN
         I---------I---------I---------I---------I---------I       ------------[X] ------------[(] ------------[)] ------------[*]
    140 I                                   *               I  140                                                   6.4068800
    141 I                        *                          I  141                                                   6.2304814
    142 I                 *                                 I  142                                                   6.1333980
    143 I       *                                           I  143                                                   5.9661467
    144 I.............*.....................................I  144                                                   6.0684256
    145 I             (----X----)                           I  145   6.1452231       6.0718823       6.2185638
    146 I    (----X-----)                                   I  146   6.0059704       5.9202519       6.0916890
    147 I        (-----X------)                             I  147   6.0836941       5.9871722       6.1802160
    148 I              (------X------)                      I  148   6.1916867       6.0854545       6.2979189
    149 I              (-------X------)                     I  149   6.2003273       6.0852010       6.3154537
    150 I                       (--------X-------)          I  150   6.3482934       6.2249123       6.4716744
    151 I                                 (-------X--------)I  151   6.4999921       6.3688751       6.6311092
    152 I                             (--------X--------)   I  152   6.4531809       6.3147595       6.5916023
    153 I                   (--------X---------)            I  153   6.3004010       6.1550419       6.4457602
    154 I            (---------X---------)                  I  154   6.2066753       6.0546948       6.3586559
    155 I(---------X----------)                             I  155   6.0242650       5.8659397       6.1825903
    156 I      (----------X---------)                       I  156   6.1220395       5.9576142       6.2864649
    157 I          (-----------X-----------)                I  157   6.2023681       6.0226252       6.3821111
    158 I(------------X-----------)                         I  158   6.0631155       5.8731573       6.2530738
    159 I     (------------X------------)                   I  159   6.1408392       5.9411876       6.3404907




























                                                             <13-29>
1G.  Acknowledgments

      The  subroutines  used  to  compute  the least squares solution are those
 referenced in Dennis  et al.  [1981].  The algorithm used  to  select  optimum
 step  sizes  for  numerical derivatives was developed by Schnabel [1982].  The
 printed output for the ARIMA estimation subroutines has been modeled  in  part
 on the linear least squares output used by OMNITAB II [Hogben et al., 1971].




















































                                    <13-30>
1-----                             Appendix A                             -----

           CONTINUITY OF VERTICAL PLOTS ON THE CDC CYBER 840 AND 855


      Normally, a line printer  will automatically  provide margins at  the top
 and bottom of each page causing a break  in the continuity of a  vertical plot
 extending  over  two  or more  pages.   However,  these  automatic page-ejects
 within a vertical plot can be suppressed by the user on many systems.   On the
 CDC Cyber 840  and 855  model machines  this  would be  done by  printing  a Q
 carriage control  in column one  immediately before  the call to  the vertical
 plot routine.   Printing an R carriage  control will cancel this effect.   For
 example, the sequence

                      WRITE (6, 100)
                 100  FORMAT (1H1, 'TITLE FOR VERTICAL PLOT')
                      WRITE (6, 101)
                 101  FORMAT ('Q')
                      CALL VP (Y, N)
                      WRITE (6, 102)
                 102  FORMAT ('R')

 will produce  a vertical plot beginning on  a new page, without any  breaks in
 continuity and without affecting the automatic page-ejects in the rest  of the
 output.  Users of other systems should consult their Computer Center staff for
 any equivalent method available.

































                                     <A-1>
1-----                             Appendix B                             -----

                             WEIGHTED LEAST SQUARES


      Weighted  least squares can  be used  to eliminate observations  from the
 analysis and to compensate for unequal variances in the observational errors.

      Observations can be eliminated  from the analysis by using  weight values
 consisting  only of  zeros and ones.   This will  produce the  same results as
 performing an unweighted analysis with the zero-weighted values removed except
 that the predicted values,  the standard  deviations of the  predicted values,
 and the residuals of the zero-weighted data are computed.   There are two main
 reasons for weighting observations zero.  The first is to obtain the predicted
 values and their standard deviations  for a  set of independent  variables not
 included in the observed data.  (This is done by assigning any arbitrary value
 to the  dependent variable of  the desired  set of independent  variables, and
 then  weighting  these values  zero.)   The  second  reason is  to allow  easy
 examination of the effect of  outliers and influential data points.   Outliers
 often appear as large values in residual plots.   Careful checking of the data
 often leads  to confirmation that the data  are in  error, and sometimes  to a
 correction.   When a  cause for  suspicious data cannot  be found,  it may  be
 advisable  to compare the  analysis with  and without  the  questionable data.
 Caution is in order if the estimates or conclusions are highly sensitive  to a
 small amount of suspicious data.   Data that  have a very high influence  on a
 fitted curve may not result in large  residuals, however, even if they  are in
 error.  In fact, extremely influential observations may force the fitted curve
 to be very close, leading to very small residuals.   It is therefore desirable
 to identify influential observations and to compare the results  obtained with
 and  without  these  points.    Several  methods   for  detecting  influential
 observations are discussed in Bement and Williams [1969], Cook [1977], Hoaglin
 and Welsch [1978], and Belsley et al. [1980].

      Using weights to compensate for unequal observational error  variances is
 not as straightforward as  using zero  weights to eliminate  observations from
 the analysis.   When the variances of the observational errors, e(i), are  not
 equal, the unweighted least squares estimates remain unbiased but do  not have
 minimum variance.   Minimum variance  estimates are obtained by using  weights
 wti = 1/Variance[e(i)] when the error variances are known.  If weights must be
 estimated, they should be based on at least 10 degrees of freedom  [see Bement
 and Williams, 1969].   In practice, however, weights are derived from  theory,
 or obtained from the data being fit, and  either of these methods can  do more
 harm than  good.   When  the  need for  weights is  suspected  and  the  error
 variances are not known, first  fit the  data using unweighted  least squares;
 analysis  of the  residuals may confirm  the need  for weighting and  may also
 provide  estimates for  the weights themselves.   If the  need for  weights is
 confirmed, then a statistician should be consulted to assist in  selecting the
 weights and in interpreting the results.











                                     <B-1>
1-----                             Appendix C                             -----

                    ESTIMATING THE NUMBER OF RELIABLE DIGITS
                         IN THE RESULTS OF A FUNCTION


      The  number of  reliable digits,  h,  in the  results of  a  real  valued
 function, g(b), can be estimated in most cases by evaluating

                                         (|g(bj) - [a+j*b]|)
                    h = -log10(    max   -------------------
                               j=-2,...,2        |g(b)|

 where

 bj      is the vector of the NPAR parameters of the function given by,

         bj(k) = b(k) + b(k)*j*10**(-DIGITS/2)     for j = -2, ..., 2,
                                                   and k = 1, ..., NPAR,

         where

         DIGITS is the number of decimal digits carried by the  user's computer
         for  a single precision  value when  the single  precision  version of
         STARPAC is being used and is the number carried for a double precision
         value otherwise.

                     2
 a       = (0.20) * SUM g(bj).
                   j=-2

                     2
 b       = (0.10) * SUM j*g(bj).
                   j=-2

      This procedure may underestimate the number of reliable digits if g(b) is
 extremely nonlinear.   A more elaborate and more robust procedure is described
 in Gill et al. [1981].





















                                     <C-1>
1-----                             Appendix D                             -----

                        LIST OF STARPAC SUBPROGRAM NAMES


 D.1  Subprograms specifically written for STARPAC

 ABSCOM   ACCDIG   ACF      ACFD     ACFDTL   ACFER    ACFF     ACFFS    ACFLST
 ACFM     ACFMN    ACFMNF   ACFMNM   ACFMS    ACFOUT   ACFS     ACFSD    ACFSDM
 ACVF     ACVFF    ACVFM    ADJLMT   AIME     AIMEC    AIMES    AIMF     AIMFS
 AIMX1    AMDRV    AMEAN    AMEANM   AMECNT   AMEDRV   AMEER    AMEFIN   AMEHDR
 AMEISM   AMEMN    AMEOUT   AMEPT1   AMEPT2   AMESTP   AMFCNT   AMFER    AMFHDR
 AMFMN    AMFOUT   AMLST    AMLST1   AOS      AOSLST   AOV1     AOV1ER   AOV1HD
 AOV1MN   AOV1S    AOV1XP   ARCOEF   ARFLT    AXPBY    BACKOP   BFS      BFSDRV
 BFSER    BFSF     BFSFS    BFSLAG   BFSM     BFSMN    BFSMS    BFSMV    BFSMVS
 BFSS     BFSV     BFSVS    CCF      CCFER    CCFF     CCFFS    CCFLST   CCFM
 CCFMN    CCFMNF   CCFMNM   CCFMS    CCFOUT   CCFS     CCFSD    CCFSDM   CCFXP
 CCVF     CCVFF    CCVFM    CDFCHI   CDFF     CDFNML   CDFT     CENTER   CHIRHO
 CMPFD    CNTR     CORR    !CORRER   CORRHD   CORRMN   CORRS    CORRXP   CPYASF
 CPYMSS   CPYVII   DCKCNT   DCKCRV   DCKDRV   DCKER    DCKFPA   DCKHDR   DCKLS
 DCKLSC   DCKLS1   DCKMN    DCKOUT   DCKZRO   DCOEF    DEMDRV   DEMOD    DEMODS
 DEMODU   DEMORD   DEMOUT   DFBW     DFBWM    DIF      DIFC     DIFM     DIFMC
 DIFSER   DOTC     DOTCM    DRV      DRV1A    DRV1B    DRV2     DRV3     DRV4A
 DRV4B    ECVF     EHDR     EIAGE    EIAGEP   EISEQ    EISGE    EISII    EISLE
 EISRNG   EIVEO    EIVEQ    EIVII    ENFFT    ERAGT    ERAGTM   ERAGTP   ERDF
 ERIODD   ERSEI    ERSGE    ERSGT    ERSIE    ERSII    ERSLF    ERSLFS   ERVGT
 ERVGTM   ERVGTP   ERVII    ERVWT    ETAMDL   EXTEND   FACTOR   FFT      FFTCT
 FFTLEN   FFTR     FITEXT   FITPT1   FITPT2   FITSXP   FITXSP   FIXPRT   FLTAR
 FLTARM   FLTMA    FLTMD    FLTSL    GENI     GENR     GETPI    GFAEST   GFARF
 GFARFS   GFORD    GFOUT    GFSEST   GFSLF    GFSLFS   GMEAN    HIPASS   HIST
 HISTC    HPCOEF   HPFLT    HSTER    HSTMN    ICNTI    ICOPY    INPERL   IPGDV
 IPGM     IPGMN    IPGMP    IPGMPS   IPGMS    IPGORD   IPGOUT   IPRINT   LDSCMP
 LLCNT    LLCNTG   LLCNTP   LLER     LLHDRG   LLHDRP   LLS      LLSMN    LLSP
 LLSPS    LLSPW    LLSPWS   LLSS     LLSW     LLSWS    LOGLMT   LOPASS   LPCOEF
 LPFLT    LSTLAG   LSTVCF   LSTVEC   MAFLT    MATPRF   MATPRT   MDFLT    MDLTS1
 MDLTS2   MDLTS3   MDL1     MDL2     MDL3     MDL4     MGS      MODSUM   MPP
 MPPC     MPPL     MPPM     MPPMC    MPPML    MSGX     MULTBP   MVCHK    MVP
 MVPC     MVPL     MVPM     MVPMC    MVPML    NCHOSE   NLCMP    NLCNT    NLCNTA
 NLCNTN   NLDRVA   NLDRVN   NLER     NLERR    NLFIN    NLHDRA   NLHDRN   NLINIT
 NLISM    NLITRP   NLMN     NLOUT    NLS      NLSC     NLSD     NLSDC    NLSDS
 NLSKL    NLSPK    NLSS     NLSUPK   NLSW     NLSWC    NLSWD    NLSWDC   NLSWDS
 NLSWS    NLSX1    NLSX2    NRAND    NRANDC   OANOVA   OBSSM2   OBSSUM   PARZEN
 PGM      PGMEST   PGMMN    PGMS     PGORD    PGOUT    PLINE    PLTCHK   PLTPLX
 PLTSYM   POLAR    PP       PPC      PPCNT    PPFCHS   PPFF     PPFNML   PPFT
 PPL      PPLMT    PPM      PPMC     PPML     PPMN     PRTCNT   RANDN   !RANDU
 RANKO    REALTR   RELCOM   REPCK    SAMPLE   SETESL   SETFRQ   SETIV    SETLAG
 SETRA    SETROW   SETRV    SLFLT    SMPLY    SPCCK    SPP      SPPC     SPPL
 SPPLTC   SPPLTD   SPPLTL   SPPM     SPPMC    SPPML    SRTIR    SRTIRR   SRTRI
 SRTRRI   STAT     STATER   STATS    STATW    STATWS   STAT1    STAT1W   STAT2
 STAT2W   STKCLR   STKGET   STKREL   STKSET   STKST    STPADJ   STPAMO   STPCNT
 STPDRV   STPER    STPHDR   STPLS    STPLSC   STPLS1   STPLS2   STPMN    STPOUT
 STPSEL   SUMBS    SUMDS    SUMID    SUMIDW   SUMOT    SUMSS    SUMTS    SUMWDS
 SUMWSS   SUMWTS   SVP      SVPC     SVPL     SVPM     SVPMC    SVPML    TAPER
 UAS      UASCFT   UASDV    UASER    UASEST   UASF     UASFS    UASORD   UASOUT
 UASS     UASV     UASVAR   UASVS    UFS      UFSDRV   UFSER    UFSEST   UFSF
 UFSFS    UFSLAG   UFSM     UFSMN    UFSMS    UFSMV    UFSMVS   UFSOUT   UFSPCV
 UFSS     UFSV     UFSVS    VCVOTF   VCVOUT   VERSP    VP       VPC      VPCNT
 VPHEAD   VPL      VPLMT    VPM      VPMC     VPML     VPMN     XACF     XAIMD

                                     <D-1>
1XAIMT    XAOV1    XBFS     XCCF     XCORR    XDCKLD   XDCKLE   XDCKLT   XDEMOD
 XDFLT    XHIST    XLLS     XNLSD    XNLSE    XNLST    XNRAND   XPGM     XPP
 XSTAT    XSTPLD   XSTPLE   XSTPLT   XUAS     XUFS     XVP      XXCH1    XXCH2
 XXCH3    XXCH4    XXCH5    XXCH6    XXCH7    XXCH8    XXCH9    XXCH10   XXCH11
 XXCH12   XXCH13


 D.2  Subprograms from NL2SOL

 ASSESS   COVCLC   DFAULT   DOTPRD   DUPDAT   GQTSTP   IMDCON   ITSMRY   LINVRT
 LITVMU   LIVMUL   LMSTEP   LSQRT    LSVMIN   LTSQAR   MADJ     MADR     NL2ITR
 NL2SNO   NL2SOL   NL2X     PARCHK   QAPPLY   QRFACT   RELDST   RMDCON   RPTMUL
 SLUPDT   SLVMUL   STOPX    UFPARM   VAXPY    VCOPY    VSCOPY   V2NORM


 D.3  Subprograms from miscellaneous public domain sources

 ALBETA   ALGAMS   ALNGAM   ALNREL   BETAI    CSEVL    DBETAI   DCSEVL   DERF
 DERFC    DGAMI    DGAMIT   DGAMLM   DGAMMA   DGAMR    DLBETA   DLGAMS   DLNGAM
 DLNREL   D9GMIT   D9LGIC   D9LGIT   D9LGMC   EPRINT   ERF      ERFC     E9RINT
 FDUMP    GAMI     GAMIT    GAMLIM   GAMMA    GAMR     INITDS   INITS    I8SAVE
 J4SAVE   R9GMIT   R9LGIC   R9LGIT   R9LGMC   SETERR   S88FMT   XERABT   XERCLR
 XERCTL   XERPRT   XERROR   XERRWV   XERSAV   XGETF    XGETUA   XSETF


 D.4  Subprograms from LINPACK and BLAS

 DASUM    DAXPY    DCOPY    DDOT     DNRM2    DSCAL    DSIDI    DSIFA    DSWAP
 DTRCO    DTRDI    IDAMAX   ISAMAX   SASUM    SAXPY    SCOPY    SDOT     SNRM2
 SSCAL    SSIDI    SSIFA    SSWAP    STRCO    STRDI


 D.5  Subprograms specifying machine dependent constants

 D1MACH   I1MACH   R1MACH
























                                     <D-2>
1-----                             Appendix E                             -----

                      LIST OF STARPAC LABELED COMMON NAMES


                            CSTAK   ERRCHK   NOTOPT





















































                                     <E-1>
1-----                             REFERENCES                             -----


 Abramowitz, M.; Stegun, I.  (1964).  Handbook of mathematical functions.  Nat.
 Bur. Stand. (U.S.) Appl. Math. Ser. 55.

 Akaike, H.  (1974).   A new  look at  statistical model identification.   IEEE
 Transactions on Automatic Control, AC-19.

 American National  Standards Institute (1977).   ANS FORTRAN  X3.9-1977.   New
 York, NY:  American National Standards Institute.

 Anderson, T. W. (1958).  An introduction to multivariate statistical analysis.
 New York, NY:  John Wiley and Sons.

 Anscombe,  F.  J.; Tukey,  J. W.  (1963).   The  examination  and analysis  of
 residuals.  Technometrics, 5:  141-160.

 Bard,  Y. (1974).   Nonlinear parameter estimation.   New York, NY:   Academic
 Press.

 Belsley, D. A.; Kuh, E.; Welsch, R. E. (1980).   Regression diagnostics.   New
 York, NY:  John Wiley and Sons.

 Bement,  T.  R.; Williams,  J. S.  (1969).   Variance  of  weighted regression
 estimators when sampling errors are independent and heteroscedastic.  J. Amer.
 Statists. Assoc.  64:  1369-1382.

 Bloomfield, P. (1976).   Fourier analysis  of time  series:   an introduction.
 New York, NY:  John Wiley and Sons.

 Box, G. E. P.;  Jenkins, G. M.  (1976).   Time series analysis forecasting and
 control.  San Francisco, CA:  Holden-Day.

 Bradley,  J. V.   (1968).   Distribution-free  statistical  tests.   Englewood
 Cliffs, NJ:  Prentice-Hall.

 Brownlee, K. A.   (1965).   Statistical theory and methodology in science  and
 technology.  Second edition.  New York, NY:  John Wiley and Sons.

 Cook,  D.  R.   (1977).   Detection  of  influential  observations  in  linear
 regression.  Technometrics, 19:  15.

 Crow,  E. L.; Siddiqui, M. M.   (1967).   Robust estimation of locations.   J.
 Amer. Stat. Assoc. 62:  353-389.

 Daniel, C.; Wood, F. S.  (1980).  Fitting equations to data.   Second edition.
 New York, NY:  John Wiley and Sons.

 Davis,  P. J. (1962).   Orthonormalization codes  in numerical  analysis.   In
 Survey of Numerical Analysis, J. Todd, ed.  New York, NY:  McGraw-Hill.

 Dennis, J. E. Jr.; Gay, D. M.; Welsch, R. E.  (1981a).   An adaptive nonlinear
 least-squares algorithm.  ACM Trans. Math. Software, 7(3):  348-368.

 Dennis, J.  E. Jr.;  Gay,  D. M.;  Welsch, R.  E.   (1981b).   Algorithm  573:
 NL2SOL - an adaptive nonlinear least squares algorithm [E4].  ACM Trans. Math.
 Software, 7(3):  369-383.

                                     <REF-1>
1
 Dennis, J. E.  Jr.; Schnabel,  R. B.   (1983).   Numerical methods  for 
 unconstrained  optimization  and  nonlinear  equations.    Englewood  Cliffs,   NJ:
 Prentice-Hall.

 Dixon,  W.  J.;  Massey,  F. J.  Jr.   (1957).   Introduction  to  statistical
 analysis.  Second edition.  New York, NY:  McGraw-Hill.

 Donaldson,  J.  R.; Schnabel,  R. B.  (1987).   Computational experience  with
 confidence  regions  and confidence  intervals for  nonlinear  least  squares.
 Technometrics, 29:  67-82.

 Donaldson,  J.  R.;  Tryon,  P. V.  (1983a).   Introduction  to  STARPAC,  the
 Standards  time series and regression package.   Nat. Bur. Stand. (U.S.) Tech.
 Note 1068-1.

 Donaldson, J.  R.; Tryon, P. V. (1983b).   Nonlinear Least  Squares Regression
 Using STARPAC, the Standards  time series  and regression package.   Nat. Bur.
 Stand. (U.S.) Tech. Note 1068-2.

 Dongarra, J. J.; Moler, C. B.; Bunch, J.  R.; Stewart, G. W. (1979).   LINPACK
 users' guide.  Philadelphia, PA:  SIAM.

 Draper,  N.  R.;  Smith, H.  (1981).   Applied  regression  analysis.   Second
 edition.  New York, NY:  John Wiley and Sons.

 Duncan,  A. J.  (1965).   Quality control  and industrial  statistics.   Third
 edition.  Richard D. Irwin.

 Eisenhart,   C.  (1947).  Significance  of the  largest  of  a  set of  sample
 estimates of variance.  In Techniques of Statistical Analysis.   C. Eisenhart,
 M. W. Hastay and W. A. Wallis, eds.  New York, NY:  McGraw-Hill.

 Filliben, J. J.   (1977).   User's guide to the DATAPAC data analysis package.
 (Unpublished  -  available   from  NIST  Statistical   Engineering   Division/
 Gaithersburg.)

 Fisher, R. A.   (1950).   Statistical methods for research workers.   Eleventh
 edition.  Edenburg:  Oliver and Boyd.

 Fox, P. A.; Hall, A. D.; Schryer, N. L.  (1978a).   Algorithm 528:   framework
 for a portable library [Z].  ACM Trans. Math. Software, 4(2):  177-188.

 Fox,  P. A.; Hall,  A. D.;  Schryer, N. L.   (1978b).   The  PORT mathematical
 subroutine library.  ACM Trans. Math. Software, 4(2):  104-126.

 Freund, J. E.; Williams, F. J. (1958).  Modern business statistics.  Englewood
 Cliffs, NJ:  Prentice-Hall.

 Fullerton, L. W.  (1977).  Portable special function routines.  In Portability
 of Numerical Software: Proceedings.   Lecture Notes in Computer Science:  Vol.
 57, W. Crowell, editor.  Oak Brook, IL:  Springer-Verlag.

 Gill, P. E.; Murray, W.; Wright, M. H.  (1981).   Practical optimization.  New
 York, NY:  Academic Press.

 Hald, A. (1952).  Statistical theory with engineering applications.  New York,
 NY:  John Wiley and Sons.

                                     <REF-2>
1
 Hoaglin, D. C.; Welsch, R. E. (1978).  The hat matrix in regression and ANOVA.
 American Statistician, 32:  17.

 Hogben, D.; Peavy, S. T.; Varner, R. N.   (1971).  OMNITAB II user's reference
 manual.  Nat. Bur. Stand. (U.S.) Tech. Note 552.

 Jenkins, G. M.; Watts, D. C.  (1968).  Spectral analysis and its applications.
 San Francisco, CA:  Holden-Day.

 Jones, R. H. (1971).   Spectrum estimation with missing observations.   Annals
 of the Institute of Statistical Mathematics, 23:  3.

 Kendall, M. G. (1948).   Rank correlation methods.   London:   Charles Griffen
 and Co.

 Kendall,  M.  G.; Stuart,  A.   (1973).   The  advanced theory  of statistics.
 Inference and Relationship, Vol. 2:  London:  Charles Griffin and Co.

 Kennedy, W. J. Jr.; Gentle, J. E.  (1980).  Statistical computing.   New York,
 NY:  Marcel-Dekker Inc.

 Ku,  H. H.  (1973).   A  user's  guide  to  the OMNITAB  command  "STATISTICAL
 ANALYSIS."  Nat. Bur. Stand. (U.S.) Tech. Note 756.

 Lawson, C.; Hanson, R.; Kincaid, D.; Krogh, F.  (1979).   Basic linear algebra
 subprograms for FORTRAN usage.  ACM Trans. Math. Software, 5(3):  308-371.

 Mandel, J. (1964).  The statistical analysis of experimental data.   New York,
 NY:  Interscience.

 Marsaglia,  G.; Tsang,  W. W. (1984).   A fast,  easily implemented method for
 sampling from decreasing  or symmetric  unimodal density  functions.   SIAM J.
 Sci. Stat. Comput. 5 (2):  349-358.

 Miller,  I.; Freund, J. E. (1977).   Probability and statistics for engineers.
 Second edition.  Englewood Cliffs, NJ:  Prentice-Hall.

 Morrison,  D. F. (1967).   Multivariate statistical  methods.   New  York, NY:
 McGraw-Hill.

 Natrella,  M. G. (1966).   Experimental statistics.   Nat. Bur.  Stand. (U.S.)
 Handb. 91.

 Owen, D. B. (1962).   Handbook of statistical tables.   Reading, MA:  
 Addison-Wesley Publishing Co.

 Ryan, T.  A.; Joiner, B.  L.; Ryan,  B. F. (1974).   Student handbook  for the
 MINITAB statistical computing system.  The Pennsylvania State University.

 Schnabel, R. B.  (1982).  Finite difference derivatives - theory and practice.
 (Unpublished - available from NIST Statistical Engineering Division/Boulder.)

 Singleton,  R. C. (1969).   An algorithm  for computing  the mixed radix  fast
 Fourier transform.  IEEE Transactions on Audio and Electroacoustics, 17:  2.

 Snedecor, G. W. (1956).  Statistical methods.  Fifth edition.  Ames, IA:  Iowa
 State University Press.

                                     <REF-3>
1
 Snedecor, G. W.; Cochran, W. G. (1967).  Statistical methods.   Sixth edition.
 Ames, IA:  Iowa State University Press.

 Tryon, P.  V.;  Donaldson,  J.  R. (1978).   STATLIB:   A  library  of Fortran
 subroutines for  statistical analysis  of experimental data.    (Unpublished -
 available from NIST Statistical Engineering Division/Boulder.)

 Waldmeier, M. (1961).   The Sunspot activity in the years 1610-1960.   Zurich:
 Schulthess.

 Walsh,  P. J. (1962).   Algorithm 127 - ORTHO.   Comm.  Assoc.  Comp. Mach. 5:
 511-513.

                                     <REF-4>
```text
[/details]
--
<!--
install.doc
1                            Installation Document
                                      for
                                 STARPAC 2.08


                     John E. Koontz and Janet R. Donaldson
                  Center for Computing and Applied Mathematics
                 National Institute of Standards and Technology
                  (Formerly the National Bureau of Standards)



                      Please direct STARPAC enquiries to:

                               Janet R. Donaldson
                 National Institute of Standards and Technology
                                 Mail Code 719
                                  325 Broadway
                             Boulder, CO 80303-3328

                          303-497-5114 (FTS 320-5114)
               Internet:  jrd@cam.nist.gov;  Bitnet:  jrd@nistcs2


                                   Disclaimer

 Certain computer equipment is identified in this document in order to describe
 the installation procedure adequately.  Such identification does not imply
 recommendation or endorsement by the National Institute of Standards and
 Technology, nor does it imply that the equipment identified is necessarily the
 best available for the purpose.

 No warranties, express or implied, are made by the distributors or developers
 that STARPAC or its constituent parts are free of error.  They should not be
 relied upon as the sole basis for solving a problem whose incorrect solution
 could result in injury to person or property.  If the programs are employed in
 such a manner, it is at the user's own risk and the distributors and
 developers disclaim all liability for such misuse.


                                    Contents

 1.  Preliminaries
 2.  Installation Procedure
     Step 1:  Prepare for Updates and Modifications
     Step 2:  Select Single or Double Precision Version
     Step 3:  Remove STARPAC Code Already Available on Target System
     Step 4:  Set Necessary Machine Dependent Values
     Step 5:  Compile STARPAC and Create Object-Code Library
     Step 6:  Test STARPAC
     Step 7:  Set Up Public Interface
     Step 8:  Set Up Distributor Interface
 3.  Updates
 4.  STARPAC User's Guide
 5.  Information for Advisers
 6.  STARPAC Contacts
 7.  Acknowledgments
 Appendix A:  STARPAC Portability
      A.1  Language Compatibility
      A.2  Handling of Floating Point Underflow
      A.3  Operating System Compatibility
      A.4  STARPAC Output Characteristics
      A.5  STARPAC Compatibility with Other Software
 Appendix B:  Comparing STARPAC Test Results
 Appendix C:  Systems Running STARPAC
 References
 STARPAC Release Installation Registration Form


                               1.  Preliminaries

 STARPAC, an improved descendent of STATLIB 2 [Tryon and Donaldson, 1978], is a
 library of portable Fortran subprograms for conducting statistical analysis.
 STARPAC consists of groups of subprograms, each group dealing with a
 particular analysis task, for example, nonlinear regression.  Within each
 group there are some subprograms with brief argument lists, and others with
 more extended lists.  The subprograms with brief argument lists require little
 more than the data as input arguments and return only the most basic results
 as output arguments, with default control values determining the details of
 the computation of the results and the printing of the report.  The
 subprograms with more extended argument lists include the control values as
 input arguments and return additional results.

 STARPAC has the following features.

 1.  It is portable.  It should be possible to install it on any system with a
     Fortran 77 compiler, an object-code library facility (see section 2, step
     5), adequate memory, and a printer with 132 characters per line (carriage
     control included).

 2.  It is easy to use, performing extensive error checking and providing
     comprehensive printed reports.

 3.  It has good user documentation, which is available on-line.

 4.  It uses reliable numerical algorithms.

 5.  It has no restrictions on problem size other than effective machine size.

 6.  It can be used with most other Fortran-callable subprogram libraries
     including IMSL [IMSL, 1982], LINPACK [Dongarra et al., 1980], and others.

 7.  It can be used with the user's own Fortran code, for example, special
     data transformations or model functions.

 The current STARPAC release is STARPAC 2.08, which contains subprogram
 groups for performing time series analysis and nonlinear least squares
 regression, as well as normal random number generation, line printer plots,
 basic statistical analyses, and linear least squares.  This version
 supersedes all earlier versions of STARPAC, including STARPAC 1, which
 contained only the nonlinear least squares portion of STARPAC 2.08.

 This document explains the layout of the STARPAC release, and provides
 installation instructions as well as supplemental material which the
 installer may find useful.  Section 5 may also be of use to STARPAC users and
 could be included in the installer's local documentation.

 The release tape includes all materials and code necessary to install
 STARPAC.


 *** PHYSICAL CHARACTERISTICS OF TAPE

      A.  ASCII character set.
      B.  1600 cpi density.
      C.  Unlabeled.
      D.  11 files, each terminated by tapemarks.
      E.  Additional tapemark follows tapemark of last file.
      F.  Files consist of 1 or more blocks (physical records).
      G.  Files 1 to 9 have blocks of
          45 line images (logical records) of 80 characters each,
          i.e., 3600 characters;
      H.  Files 10 and 11 have blocks of
          20 line images (logical records) of 132 characters each,
          i.e., 2640 characters;
      I.  Last block of a file may contain fewer than the specified number of
          line images, in which case it is short, not blank filled.



 *** TAPE CONTENTS

 File No.  File Id.       Description
 --------  -----------    -----------

        1  TOC.DOC      - tape characteristics, file structure and
                          table of contents
                          (line image length = 80, block size = 3600)

        2  INSTALL.DOC  - this installation manual
                          (line image length = 80, block size = 3600)

        3  S_SRCE.FOR   - single precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (line image length = 80, block size = 3600)

        4  S_NL2SOL.FOR - single precision STARPAC 2.08 source code from NL2SOL
                          (line image length = 80, block size = 3600)

        5  D_SRCE.FOR   - double precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (line image length = 80, block size = 3600)

        6  D_NL2SOL.FOR - double precision STARPAC 2.08 source code from NL2SOL
                          (line image length = 80, block size = 3600)

        7  MISC.FOR     - single and double precision subprograms from
                          miscellaneous public domain sources, but excluding
                          code from NL2SOL, LINPACK and BLAS and excluding
                          subprograms supplying machine dependent constants
                          (line image length = 80, block size = 3600)

        8  LPK_BLS.FOR  - single and double precision subprograms
                          from LINPACK and BLAS
                          (line image length = 80, block size = 3600)

        9  MDC.FOR      - single and double precision subprograms supplying
                          machine dependent constants
                          (line image length = 80, block size = 3600)

       10  TEST.OUT     - results of test subprograms
                          (line image length = 132, block size = 2640)

       11  GUIDE.DOC    - STARPAC 2.08 user documentation,
                          in line printer format
                          (line image length = 132, block size = 2640)


 In the remaining sections, the information necessary to install and support
 STARPAC 2.08 is supplied.  Section 2 provides a step by step description of
 the installation procedure.  Section 3 describes the update procedure that
 will be used for future modifications and additions, and section 4 describes
 the STARPAC user documentation.  Section 5 lists information we believe will
 be useful for the local STARPAC advisors, and section 6 lists the names and
 addresses of STARPAC contacts at the National Institute of Standards and
 Technology.  The public domain software used by STARPAC is acknowledged in
 section 7.  In appendix A, we summarize STARPAC characteristics which may
 affect portability, and in appendix B we provide information which should
 help installers determine whether differences they observe between their
 STARPAC output and that obtained by us are significant.  Finally, in appendix
 C we list systems running STARPAC.


                           2.  Installation Procedure

 STARPAC and its constituent parts are written in Fortran as defined in the
 1978 standard (ANSI X3.9-1978) [ANSI, 1978], commonly called Fortran 77.  (The
 terms "Fortran V" or "Fortran 5" are occasionally used as synonyms, but may
 refer instead to "enhanced Fortran IV.") We believe that STARPAC is compatible
 with the ANSI X3.9-1978 full language standard.  However, before carrying out
 the installation procedure, it may be useful to glance at appendix A, which
 surveys some potential installation problems.


 Step 1:  Prepare for Updates and Modifications

 It is an unfortunate fact of computing that code always requires some
 modification during its lifetime.  Since STARPAC requires that a few
 modifications be made during the initial installation, we recommend that a
 maintenance system be organized as the first step of the STARPAC installation
 process.  If possible, the installers should set up such a system for the
 STARPAC source code and object code to make it easy to retrieve and modify
 individual subprograms.  The details of this will depend on the target
 system.  On Cyber computers running NOS, for example, the source code could
 be maintained in UPDATE or MODIFY format, and the object code could be
 handled with LIBEDIT.

 Note that the subprograms in each source code file are ordered alphabetically,
 and that each STARPAC subprogram begins with a comment line consisting of an
 asterisk followed immediately by the subprogram name, i.e., *name.  (This is
 the only use of an asterisk in column 1 of the STARPAC source code.) This is
 done to facilitate the creation of the maintenance system and to aid the
 installer in separating the subprograms.  For example, the following Fortran
 program could be used to split source code file D_SRCE.FOR.

      PROGRAM FSPLIT
      CHARACTER TEST1*1,TEST2*6,LINE*73
      INTEGER IUNIT1,IUNIT2,LAST

      IUNIT1 = 30
      IUNIT2 = 31

      OPEN (UNIT=IUNIT1,FILE='D_SRCE.FOR')

      READ (UNIT=IUNIT1,FMT=1000,END=999) TEST1,TEST2,LINE
      LAST = INDEX(TEST2,' ')-1
      IF (LAST.LE.0) LAST = 6
      OPEN (UNIT=IUNIT2,FILE=TEST2(1:LAST)//'.FOR')
      WRITE (UNIT=IUNIT2,FMT=1000) TEST1,TEST2,LINE

  100 READ (UNIT=IUNIT1,FMT=1000,END=999) TEST1,TEST2,LINE
      IF (TEST1.EQ.'*') THEN
         CLOSE (UNIT=IUNIT2)
         LAST = INDEX(TEST2,' ')-1
         IF (LAST.LE.0) LAST = 6
         OPEN (UNIT=IUNIT2,FILE=TEST2(1:LAST)//'.FOR')
      END IF
      WRITE (UNIT=IUNIT2,FMT=1000) TEST1,TEST2,LINE
      GO TO 100

  999 CONTINUE
      CLOSE (UNIT=IUNIT2)

 1000 FORMAT (A1,A6,A73)
      END


 Step 2:  Select Single or Double Precision Version

 STARPAC is supplied in both single and double precision versions.  Both
 versions are complete as they stand, and except for precision, the two
 versions are identical.  Installers can use either precision version.  The
 two versions cannot be used together, however, because the names of most of
 the subprograms are the same in both versions.  Both versions are documented
 in the STARPAC User's Guide.

 Portions of STARPAC, for example the nonlinear least squares procedures, are
 sensitive to the machine precision and require approximately 14 decimal
 places for rational arithmetic.  Somewhat fewer places should still work, but
 six or seven decimal places are definitely too few for general use.  For
 these numerically sensitive procedures, only the simplest problems could be
 solved correctly at such reduced precisions.

 The installers must choose which version of STARPAC to install based upon
 which version supplies adequate precision on the target machine.  As far as
 we know, at present only CDC and Cray Fortrans offer sufficient precision to
 permit general use of the single precision version of STARPAC.  For other
 machines, we recommend general use of the double precision version only.

 If both versions of STARPAC have sufficient precision on the installers'
 machine, then both may be installed; however, since the bulk of the STARPAC
 subprograms have identical names in both the single precision and double
 precision versions, it will be necessary to keep the two versions of STARPAC
 in separate libraries and to use them separately.  Note that when both the
 single and double precision versions are available, there are likely to be
 trade-offs between them.  The double precision version will offer greater
 accuracy in results, while the single precision version will require less
 storage and machine time.


 Step 3:  Remove STARPAC Code Already Available on Target System

 STARPAC 2.08 contains code written specifically for STARPAC, as well as a
 variety of subprograms taken from public domain sources.  Files S_SRCE.FOR
 and D_SRCE.FOR contain those subprograms written specifically for STARPAC;
 files S_NL2SOL.FOR, D_NL2SOL.FOR, MISC.FOR, LPK_BLS.FOR and MDC.FOR contain
 subprograms extracted from public domain sources.  Subprograms within the
 first group must be used as supplied.  The installer may substitute already
 available versions of subprograms in the second group, if desired, subject to
 the limitations listed below.  We strongly advise, however, that STARPAC be
 installed and tested first with the code supplied, so that the installer can
 determine whether substitution results in any changes in STARPAC's behavior.
 The subprograms taken from public domain sources are described in greater
 detail below.

 Files S_NL2SOL.FOR and D_NL2SOL.FOR contain single and double precision
 versions, respectively, of the NL2SOL Adaptive Nonlinear Least Squares
 Algorithm, developed by John E.  Dennis, Jr., David M.  Gay, and Roy E.
 Welsch.  These subprograms supply the core of the STARPAC nonlinear
 regression capabilities.  We use Version 2.2, obtained directly from David M.
 Gay, with some minor changes made in the interest of increased portability.
 We supply with this release of STARPAC the whole of NL2SOL Version 2.2, as
 modified.  Users who wish to use it directly can do so as documented in
 Dennis et al.  [1981a,b].  The NL2SOL developers have continued their work
 since their release of version 2.2.  The most recent versions of NL2SOL
 CANNOT be substituted directly for the version used in STARPAC, and we cannot
 guarantee the compatibility of STARPAC with versions of NL2SOL other than
 version 2.2.

 File MISC.FOR contains single and double precision versions of a variety of
 public domain subprograms, including subprograms from ACM Algorithm 528,
 Framework for a Portable Library [Fox et al., 1978a,b] (excluding subprograms
 D1MACH, I1MACH and R1MACH), and from the portable special function subprogram
 library of Fullerton [1977].  These subprograms may be replaced with locally
 available versions, if desired.

 File LPK_BLS.FOR contains the subset of LINPACK [Dongarra et al., 1979] and
 BLAS [Lawton et al., 1979a,b] subprograms specifically required by STARPAC.
 If LINPACK and the BLAS have already been installed on the local system, the
 previously installed subprograms may be freely substituted for the
 subprograms supplied on the STARPAC tape.  This may be especially desirable
 on machines for which optimized versions of LINPACK and BLAS exist.

 File MDC.FOR contains the 3 subprograms, I1MACH, R1MACH and D1MACH, that
 supply the machine-dependent constants necessary for the single and double
 precision versions of STARPAC to run.  I1MACH, R1MACH and D1MACH are all part
 of ACM Algorithm 528, Framework for a Portable Library [Fox et al., 1978a,b],
 developed at Bell Laboratories.  If this code is already available on the
 target machine, it should be possible to substitute the installed versions of
 these subprograms instead of modifying the STARPAC versions as described in
 the next step.  Installers must be careful, however, to verify that the
 previously installed versions have not been extended so that they are
 incompatible with the STARPAC versions.


 Step 4:  Set Necessary Machine Dependent Values

 STARPAC 2.08 has been designed so that the number of machine dependencies
 that the installer must address are minimized.  For most sites, only the
 three subprograms in file MDC.FOR must be modified for installation.
 However, subprograms STKSET and GETPI in files S_SRCE.FOR and D_SRCE.FOR and
 subprograms SNRM2 and DNRM2 in file LPK_BLS.FOR may also require modification.
 These necessary and possible changes are described in detail below.

 File MDC.FOR contains the 3 subprograms, I1MACH, R1MACH and D1MACH, that
 supply the machine-dependent constants necessary for STARPAC to run.  As
 stated in the previous step, these 3 subprograms are from ACM Algorithm 528,
 developed at Bell Laboratories.  I1MACH supplies integer valued machine
 dependent constants.  R1MACH, which is actually called only in the single
 precision version of STARPAC, supplies single precision machine dependent
 constants.  D1MACH, which is actually called only in the double precision
 version of STARPAC, is the double precision equivalent of R1MACH.

 As noted in the previous step, it may be possible to substitute a previously
 installed version of these subprograms.  If these subprograms are not already
 available on the target machine, the installer will have to update the
 supplied versions so that they return the necessary values.  I1MACH, R1MACH
 and D1MACH as supplied on the tape will return undefined values if they are
 not updated.  All three subprograms are commented to describe the necessary
 changes.  The required constants for a number of common machines are listed
 in DATA statements in these subprograms.  If the target machine's constants
 are available in this fashion, then all that is necessary is to uncomment the
 relevant DATA statements.  If the target machine's constants are not already
 listed in the DATA statements of these subprograms, then the installer will
 need to obtain the appropriate values from local machine documentation, the
 local systems group, or, if necessary, the machine's manufacturer.  WE CANNOT
 HELP YOU DETERMINE THE NECESSARY VALUES.

 Subprogram STKSET in files S_SRCE.FOR and D_SRCE.FOR supplies information on
 the length in storage units of integer, logical, real, double precision, and
 complex values.  STKSET assumes that integers, logicals, and reals are all
 one storage unit long, and double precisions and complexes are two storage
 units long.  Changes are only needed if these assumptions do not hold, which
 is the case, for example, on some mini- and micro-computers.  STARPAC's
 STKSET is based on the Framework subprograms ISTKIN and I0TK00 [Fox et al.,
 1978a,b].  The part that may require modification is from I0TK00.  The basis
 for the assumptions in the original I0TK00 subprogram is explained in Fox et
 al.  [1978a:  123-4].  An installed version of I0TK00 can supply a guideline
 for changes to STKSET; otherwise rely on the comments in STKSET for help.
 (Caution:  Do not simply insert a call to an available I0TK00.  Unlike
 STKSET, I0TK00 does not reinitialize the stack when a second STARPAC
 subprogram is called.)

 Subprogram GETPI, also in files S_SRCE.FOR and D_SRCE.FOR, defines the
 constant pi to 31 digits.  If the precision for STARPAC as installed on the
 target machine exceeds 31 digits, GETPI should be altered appropriately.

 Subprograms SNRM2 and DNRM2 in file LPK_BLS.FOR may also require
 modification.  These two subprograms compute the L2 norm of a vector.  SNRM2
 is called only in the single precision version of STARPAC, while DNRM2 is
 called only in the double precision version.  They contain two constants,
 CUTLO and CUTHI, which are defined for worst cases over a large and current
 range of machines.

      CUTLO = max(SQRT(U/EPS))
            = 4.441E-16 (single precision)
            = 8.232D-11 (double precision)
      CUTHI = min(SQRT(V))
            = 1.304E19 (single precision)
            = 1.304D19 (double precision)

      where

              EPS = smallest number such that EPS + 1.0 > 1.0,
              U   = smallest positive number (underflow limit),
              V   = largest number (overflow limit).

 It is not likely that the installers' machine will constitute a still worse
 case, but, if it does, CUTLO and CUTHI should be redefined.  Note, however,
 that SNRM2 and DNRM2 are both part of the BLAS [Lawton et al., 1979a,b].
 Before considering changing SNRM2 or DMRM2, therefore, the installers should
 determine whether or not the BLAS are installed on their machine already.  If
 they are, then it should be possible to make use of the versions already
 installed.  This is particularly useful if the installers' machine is one of
 those for which specially optimized versions of the BLAS exist.  [See step 3
 above.]

 Note that a few other constants in the STARPAC code are carried only to a
 limited number of places; consequently, results from calculations dependent
 on these constants may not be accurate to more than approximately this number
 of places.  At present, however, we know of no instances in which the
 resultant accuracy is not adequate for STARPAC's purposes.


 Step 5:  Compile STARPAC and Create Object-Code Library

 The single and double precision versions of the STARPAC 2.08 source code are
 each contained in 5 files.  For the single precision version, files
 S_SRCE.FOR, S_NL2SOL.FOR, MISC.FOR, LPK_BLS.FOR and MDC.FOR should be used.
 For the double precision version, files D_SRCE.FOR, D_NL2SOL.FOR, MISC.FOR,
 LPK_BLS.FOR and MDC.FOR should be used.  In each case, the code must be
 modified as described in steps 3 and 4 above.

 In compiling STARPAC, we recommend that the following compiler options be
 used if available.

 1.  Rounded arithmetic (optional with at least CDC compilers).

 2.  No subscript checking.  STARPAC contains correct subscript references
     that access elements beyond the nominal upper limit of certain actual
     parameter arrays.  Run-time subscript checking would produce a deluge of
     incorrect run-time error messages.

 3.  The most extensive set of error messages possible.  In our experience it
     is well worth hearing everything that a compiler has to say about imported
     code.

 Let us emphasize that we do not expect any problems with compilation.
 Modification of some subprograms may be required, however, by features in the
 installer's Fortran compiler or operating system.  Appendix A lists some
 potential problems which the installer may have to address.

 After the STARPAC code has been successfully compiled, a Fortran object-code
 library must be created.  The term object-code library here refers to whatever
 facility a system has for satisfying external references automatically at
 load-time, using a collection of previously compiled subprograms.  The
 existence of such a facility is essential to the convenient use of STARPAC,
 since each user-callable STARPAC subprogram invokes a whole series of
 additional subprograms.  A user would find it very awkward to deal with these
 additional subprograms explicitly.

 Note that the 5 files containing the STARPAC source code should be compiled
 into a single library, not 5 separate libraries.  Creating 5 separate
 libraries and linking them sequentially could result in unsatisfied externals
 since in some instances subprograms in one file reference subprograms from
 other files.


 Step 6:  Test STARPAC

 To make it possible to verify that STARPAC functions reasonably on the
 installers' machine, we have included two sets of test subprograms in files
 S_SRCE.FOR and D_SRCE.FOR.  The first set, which includes all subprograms in
 these two files that begin with XX, e.g., subprogram XXCH1, provides a
 minimal set of tests for determining that STARPAC is running properly.  The
 second set, which includes all subprograms that begin with a single X, e.g.,
 XACF, provides more extensive tests.  Only the first set is documented here.
 Instructions for using the more extensive tests found in the second set are
 available from us if the results from the first set are inconclusive or
 suggest special problems.

 Each test subprogram beginning with XX is associated with the subprograms
 documented in a given chapter of the STARPAC User's Guide, one test
 subprogram per chapter.  The output from these tests is provided in file
 TEST.OUT, which contains Cyber 180/855 output from the single precision
 version of STARPAC (machine precision c. 14 decimal places).

 The STARPAC test subprograms should be invoked using a program of the form

      PARAMETER (LDSTAK=3000)
      DOUBLE PRECISION DSTAK(LDSTAK)
      COMMON /CSTAK/ DSTAK
      COMMON /ERRCHK/ IERR

      CALL IPRINT(IPRT)
      OPEN (UNIT=IPRT,FILE='TESTFL')

      CALL XXCH<n> (LDSTAK)

      END

 where <n> is an integer between 1 and 13, inclusive.

 Each test subprogram tests the user callable STARPAC subprograms documented
 in a given chapter of the STARPAC User's Guide.  Test subprogram names are
 constructed by appending the chapter number to the characters XXCH.  That is,
 subprogram XXCH2 is the test subprogram associated with the STARPAC user
 callable subprograms documented in chapter 2 of the STARPAC User's Guide;
 XXCH3 is the test subprogram associated with subprograms documented in
 chapter 3; etc.  Each test subprogram essentially duplicates the examples
 shown in the given chapter of the STARPAC User's Guide.  (A test subprogram
 for chapter 1, XXCH1, is also available.  It duplicates the test subprograms
 for chapters 2 and 5, however, and it is therefore not absolutely necessary
 to run it.) The single argument of these test subprograms, e.g., LDSTAK in
 the above example, is the amount of work space provided to STARPAC in common
 /CSTAK/ [see chapter 1, section D.2, of the STARPAC User's Guide].  Note that
 the output for each test is preceeded by a line containing the string 1*CH<n>.

 It is possible and sometimes convenient to combine several test runs into one
 by calling one test subprogram after another in a single test program.  There
 could be a local limit to this, however, since at some point the memory
 required by the test's object code and run time storage may exceed the memory
 available on the machine.

 To verify that the installed version of STARPAC is performing correctly, each
 of the test subprograms should be run and their results compared with the
 appropriate test outputs from the release tape.  There will be at least some
 differences between the two sets of outputs.  Appendix B provides comments
 that are intended to help the installer decide if the observed differences
 are significant.  Differences described as insignificant can be ignored.  Any
 others differences are potential problems.  If a problem proves intractable,
 the installers may feel free to contact the distributors [see section 6],
 though we cannot guarantee that we will be able to respond to the problem
 immediately, or to solve it when we can respond.

 In comparing the tape and installed versions of STARPAC test outputs, please
 note that the comparison should be made BEFORE any updates are inserted that
 were included separately with the release or supplied later.  Such updates
 will not be reflected in the output on the tape, and might easily produce
 additional differences between the tape test output and the test output from
 the installed version of STARPAC.  After the tape and installed test outputs
 compare satisfactorily, the additional updates can be inserted and the
 pre-update and post-update test outputs can be compared usefully.

 The optimal way of comparing two test outputs is to use a text comparison
 utility (comparator).  It may be useful to pre-edit one or both of the test
 files before using such a test comparator, however, in order to remove the
 pervasive but insignificant disagreements due to differences in format
 conversion [see appendix B].  It is best to use a comparator that realigns
 the files after any mismatch, reports the extraneous material, and then keeps
 going.  A literal-minded, line-for-line comparator will require a lot of
 human intervention to work at all, although it might still be preferable to
 visual comparison.  The installers must also recognize that some of the the
 test files consist essentially of a series of nearly identical reports.  A
 realigning comparator may get out of synchronization and begin comparing the
 ith report in a series in one file with the (i+1)th report in the other file.

 Note that two significant formatting differences can occur which are not
 easily detected by a comparison utility.  These involve deviations from
 STARPAC's two assumptions about output devices:  (1) that the printer page
 has a width of 132 characters (including the carriage control character), and
 (2) that the page length is at least 60 lines [see appendix A.3].  The
 installer should look for signs of clipped or awkwardly continued lines and
 for system-forced new pages followed closely by STARPAC-forced new pages.
 The best way to find both sorts of problems quickly is to look at an output
 page containing a scatter plot (test XXCH2).


 Step 7:  Set Up Public Interface

 Even when the package is running, some important steps still remain.  One
 such step is to set up the public, or user, interface.  This activity is
 necessary to some extent even if the installers are the only prospective
 users; it is particularly important if a wider circle of users is
 anticipated.

 The public interface consists essentially of three items:  (1) an
 announcement to interested parties that the package is available; (2) the
 local user document, perhaps incorporated in the announcement; and (3) a
 distribution scheme for the STARPAC User's Guide [see section 4].  The local
 user document is necessary because the STARPAC User's Guide says nothing
 about accessing STARPAC on any particular system.  The local user document
 should indicate at least (1) the current version of STARPAC, (2) how to
 access and use the object code library file, (3) how to obtain the STARPAC
 User's Guide, and (4) how to contact the installers or other designated
 advisers in the event of problems.  See section 5 for other information which
 might usefully be included in a local document.

 It is the installers' responsibility to arrange a local distribution scheme
 for the STARPAC User's Guide.


 Step 8:  Set Up Distributor Interface

 Another important but often forgotten step in installation is setting up the
 distributor interface, i.e., reporting back to us.  The principal reason for
 this is to register the installation for receipt of source code updates and
 STARPAC User's Guide revision sets.  To register, fill out the registration
 form at the end of this document and return it to us.  IF WE DO NOT RECEIVE A
 REGISTRATION FORM, WE WILL ASSUME THAT STARPAC WAS NOT INSTALLED AND NO
 FUTURE UPDATES WILL BE NEEDED.

 In addition, we would be pleased if the installers would let us know their
 experiences with installing STARPAC, particularly any problems arising from
 the nature of the code or the installation procedure.  Information of this
 nature will enable us to improve future releases.

 Finally, we would be interested in seeing any local STARPAC documentation
 that is prepared.  This will give us some idea what STARPAC looks like on
 different systems, and may suggest improvements in the STARPAC User's Guide.


                                  3.  Updates

 As indicated in section 2, step 1, it may be necessary to install an update
 to the source code at some later date.  Updates will always change the
 STARPAC version number.  Versions of STARPAC are identified with a system of
 decimal numbers:  1.01, 1.02, ..., 2.01, etc., in which the whole number part
 identifies a major release (resulting from addition of new capabilities),
 while the fractional part reflects addition of a minor modification or
 correction.  The current version of STARPAC is 2.08.

 Existing updates should already be included in any release tape and will be
 reflected in the test output of the tape.  Installers who have a STARPAC
 release and want to obtain future updates can do so by sending in the STARPAC
 Installation Registration form supplied at the end of this document.

 Updates will identify the program units modified and, for each program unit,
 show the necessary changes.  Updates will also indicate what test subprograms
 [see section 2, step 6] should be run to verify that the update functions
 properly.  Since the test subprograms may not be structured to reveal the
 error being corrected, the results of running the test may simply be to
 confirm that nothing has been accidentally corrupted during insertion of the
 update.  Any changes that the update should produce will be noted.

 STARPAC 2.08 - corrects an error in STARPAC's ARIMA modeling and forecasting
                subroutines.  This error causes the estimated model parameter
                values and the computed forecasts to be incorrect when the
                series mean is not zero.  Errors in the reported covariance
                matrix and the standard errors of the estimated ARIMA model
                parameters are also corrected.

                In addition, STARPAC 2.08 includes a modification to the
                uniform random number generator, used by NRAND and NRANDC,
                that will result in the generation of different random numbers
                than those produced by STARPAC 2.07.  It also includes a minor
                modification to the subroutine used to compute the number of
                good digits in the results of a user-supplied nonlinear model,
                used by the NLS family of routines, and changes in the
                declaration blocks for each of the STARPAC subroutines that
                will result in improved portability.

                No changes were made to the user interface to STARPAC.



                            4.  STARPAC User's Guide

 The STARPAC User's Guide is supplied in line printer format in file
 GUIDE.DOC.  This document should be referenced as National Institute of
 Standards and Technology Internal Report NBSIR 86-3448 [Donaldson and Tryon,
 1987].  Future modifications and corrections to the on-line documentation will
 be provided to registered STARPAC installations.  It is the installer's
 responsibility to make this documentation and future modifications to it
 available to local users.

 The on-line documentation uses the standard Fortran conventions for printing
 formatted records.  In particular, the first character of each line is used
 to determine vertical spacing.  For correct printing, if the first character
 of a line is a blank then the printer should vertical space one line, if the
 first character is a + then there should be no vertical advance before
 printing the line, and if the first character is a 1 then the printer should
 advance to the beginning of a new page.  The maximum page size is 132 columns
 (printer control included) by 60 lines.  (Note that the on-line documentation
 actually consists of two types of pages:  descriptive text, and sample output
 from the STARPAC subprograms.  The descriptive text has a maximum width of
 only 80 columns, while the sample output has a maximum width of 132 columns.)

 The installer may modify the on-line documentation as necessary if its format
 is not convenient.  Each section of the documentation begins with a line
 beginning with '1-----'.  This is done in part to aid the installer in
 separating file GUIDE.DOC into smaller units that may be more convenient to
 modify.  For example, the following Fortran program could be used to subdivide
 file GUIDE.DOC.

      PROGRAM CSPLIT
      INTEGER IUNIT1,IUNIT2
      CHARACTER ID*6,TEST*6,LINE*74
      CHARACTER*6 FILEID(31:50)

      DATA ID        /'1-----'/
      DATA FILEID(31)/'PREFAC'/
      DATA FILEID(32)/'CHAP01'/
      DATA FILEID(33)/'CHAP02'/
      DATA FILEID(34)/'CHAP03'/
      DATA FILEID(35)/'CHAP04'/
      DATA FILEID(36)/'CHAP05'/
      DATA FILEID(37)/'CHAP06'/
      DATA FILEID(38)/'CHAP07'/
      DATA FILEID(39)/'CHAP08'/
      DATA FILEID(40)/'CHAP09'/
      DATA FILEID(41)/'CHAP10'/
      DATA FILEID(42)/'CHAP11'/
      DATA FILEID(43)/'CHAP12'/
      DATA FILEID(44)/'CHAP13'/
      DATA FILEID(45)/'APPENA'/
      DATA FILEID(46)/'APPENB'/
      DATA FILEID(47)/'APPENC'/
      DATA FILEID(48)/'APPEND'/
      DATA FILEID(49)/'APPENE'/
      DATA FILEID(50)/'REFERS'/

      IUNIT1 = 30
      OPEN (UNIT=IUNIT1,FILE='GUIDE.DOC')

      IUNIT2 = 31
      OPEN  (UNIT=IUNIT2,FILE=FILEID(IUNIT2))

      READ  (UNIT=IUNIT1,FMT=1000,END=999) TEST,LINE
      WRITE (UNIT=IUNIT2,FMT=1000) TEST,LINE
  100 READ  (UNIT=IUNIT1,FMT=1000,END=999) TEST,LINE
      IF (TEST.EQ.ID) THEN
         CLOSE (UNIT=IUNIT2)
         IUNIT2 = IUNIT2+1
         OPEN (UNIT=IUNIT2,FILE=FILEID(IUNIT2))
         PRINT *, IUNIT2
      END IF
      WRITE (UNIT=IUNIT2,FMT=1000) TEST,LINE
      GO TO 100
  999 CONTINUE
      CLOSE (UNIT=IUNIT2)

 1000 FORMAT (A6,A74)
      END


                          5.  Information for Advisers

 In general, advisers should be familiar with the contents of chapter 1 of the
 STARPAC User's Guide.  Most of it need only be skimmed by an experienced
 Fortran user; however, STARPAC advisors should particularly note the
 following.

 1.  The aspect of STARPAC which most often confuses users is the working
     storage stack DSTAK.  This is described in chapter 1, section D.2.  The
     basic system, taken from the Framework package, is discussed in Fox et
     al.  [1978a:  116-121; 1978b:  184-188].

 2.  The unit number used for STARPAC output can be retrieved with subprogram
     IPRINT:  CALL IPRINT(IPRT) results in the unit number for output being
     returned in the INTEGER variable IPRT.  Instructions for altering the
     unit number to which STARPAC directs its output are found in chapter 1,
     section D.4.

 3.  Chapter 1, section D.6, contains a list of the most common programming
     errors made when using STARPAC.  A very high proportion of all errors
     fall into the classes listed there; the use of variably-dimensioned
     arrays in passing data and returning results seems to be especially
     troublesome.  Many other errors result from undefined variables.

 Most STARPAC error messages are self-explanatory.  There are a few, however,
 which are cryptic.  They are:

      ERROR 1 IN D1MACH - I OUT OF BOUNDS

      DSTAK IS TOO SHORT

      ***** ERROR *****
      THE POINTER FOR ALLOCATION NUMBER  <n> HAS BEEN OVERWRITTEN

      ERROR IN STARPAC. LSTVEC TRIES TO ACCESS MORE ELEMENTS THAN EXIST IN MASK.

 These are all "second string" error messages which can only occur due to
 either internal errors in the STARPAC code, or overwriting of code or data
 storage.  That is, the user should not encounter them unless there is a
 serious problem, either with STARPAC or with the user's own code.  If the
 error is actually in STARPAC, we would appreciate hearing about it.


                                  6.  Contacts

 The installers may contact the distributors as follows.

 Principal Contact:

       Janet R. Donaldson
       National Institute of Standards and Technology
       Mail Code 719
       325 Broadway
       Boulder, CO  80303-3328

       303-497-5114 or FTS 320-5114
       Internet:  jrd@cam.nist.gov;  Bitnet:  jrd@nistcs2


 Backup:

       John E. Koontz
       National Institute of Standards and Technology
       Mail Code 715
       325 Broadway
       Boulder, CO  80303-3328

       303-497-5180 or FTS 320-5180
       Internet:  jek@cam.nist.gov;  Bitnet:  jek@nistcs2


                              7.  Acknowledgments

 STARPAC uses various public-domain code; we would like to acknowledge our use
 of code from the following sources.

 1.  NL2SOL [Dennis et al., 1981a,b; see also section 2, step 3]

 2.  Framework for a Portable Library [Fox et al., 1978b]

 3.  LINPACK [Dongarra et al., 1979]

 4.  The BLAS [Lawson et al., 1979a,b]

 5.  Fullerton's Portable Special Function Subprograms [Fullerton, 1977]

 6.  OMNITAB [Hogben et al., 1971].  (STARPAC derives historically from
     STATLIB 2 [Tryon and Donaldson, 1978], and STATLIB 2 was in part a library
     variant of portions of OMNITAB.)

 7.  DATAPAC [Filliben, 1977].  (DATAPAC supplies certain subprograms used for
     probability distribution computations.)

 The STARPAC release has benefited from suggestions by Elsie Clark, Richard
 Freemire, David Schrader, and Wendell Slattery.

 STARPAC has been granted development and testing facilities on systems
 operated by the MASC Computer Services Division, NOAA (U.S.  Department of
 Commerce); the CAM Scientific Computing Division, NBS (U.S.  Department of
 Commerce); the NCAR Scientific Computing Division; and the USFS Forest Fire
 Laboratory (U.S.  Department of the Interior).  We would like to thank these
 groups, and particularly Ginger Caldwell (NCAR), Elsie Clark (CAM), Richard
 Freemire (CAM), Francis Fujioka (USFS), Sally Howe (CAM), Ken Walton (NCAR),
 and Wesley Wilson (NCAR), who participated in test installations.

 STARPAC's conception and inception owe much to Peter V. Tryon (1941-1982).


                        Appendix A:  STARPAC Portability

 This section is designed to outline and summarize considerations affecting
 STARPAC portability that are not covered in the step by step installation
 instructions of section 2.

 A.1  Language Compatibility

 STARPAC is written in standard Fortran 77 as defined by the ANSI X3.9-1978
 document.  We therefore expect a minimum number of problems arising from
 language incompatibility.  There are, however, a few potential problems.

 First, some compilers restrict the number of arguments which a subprogram may
 have.  When this problem occurs, it is most easily solved by deleting from
 the formal and actual parameter lists a sufficient number of scalar arguments
 of a single type and either (a) passing them via a labeled common or (b)
 packing them into a single array and passing them via the argument list.  (A
 list of labeled common names already used by STARPAC and therefore not to be
 duplicated for this purpose is provided in the STARPAC User's Guide in
 Appendix E.)  Details of this are left for the installer.

 Second, some compilers may flag some STARPAC code as potentially dangerous.

 1.  STARPAC uses some subscripts that appear to a range checker to be out of
     range but, in fact, are not.  Subscript range checking options should be
     turned off when STARPAC is compiled.

 2.  STARPAC equates actual expression parameters to formal parameters which
     could be defined in the subprogram.  Redefinition of the formal parameter
     does not occur in these cases.

 3.  STARPAC equates two formal parameters to one actual parameter during
     subprogram invocations.  If one formal parameter were redefined, the
     value of the other would be undefined under the Standard.  We know of no
     instances in STARPAC where the potential hazards of this operation are
     realized.

 Finally, STARPAC makes two nonstandard assumptions that a compiler may flag.
 (Fox et al., [1978a:  122-123] provide the justification for these
 assumptions.)

 1.  It uses noncorner elements in EQUIVALENCEs of mixed data types.  This
     occurs in Framework subprograms R1MACH and D1MACH, and has been done in
     such a fashion that it should not cause problems.

 2.  It assumes that if a variable is set in a DATA statement and then reset
     by an assignment statement, leaving the subprogram and reentering it will
     not cause the DATA statement value to be restored.  That is, it assumes
     that DATA initialization occurs at initial load time, not whenever a
     subprogram is (re)invoked.


 A.2  Handling of Floating Point Underflow

 For some data sets, STARPAC may result in a floating point underflow.  We
 assume that the result of such an underflow will be set to zero by the target
 system.  For all cases of underflow that we have observed in running STARPAC,
 this produces the desired result.

 We expect that STARPAC users will observe floating point underflows only
 infrequently.  If this is not the case on the target system, please notify
 us.


 A.3  Operating System Compatibility

 STARPAC requires very little of the operating system.  The only direct
 interaction is a requirement that STARPAC know the standard unit number for
 output.  The installers must set this in subprogram I1MACH in the Framework
 portion of STARPAC [see section 2, step 4].  Note that I1MACH returns
 undefined values until modified by the installer.

 Less directly, STARPAC requires some sort of object code library facility if
 STARPAC is to be used conveniently [see section 2, step 5].

 Finally, some STARPAC subprograms require a fair amount of space to run.  The
 installers may need to do some sort of segmenting if their system is small.
 The extent to which this is actually necessary can be determined by running
 the various test subprograms to see if they load.


 A.4  STARPAC Output Characteristics

 STARPAC expects output devices with pages 132 characters wide (1 carriage
 control character plus 131 printing characters), and 60 lines long.  (Note:
 Some STARPAC reports are not paginated, and run for an arbitrary number of
 lines.) Because of the paginated reports, it is advisable to alter or
 suppress, if possible, any system-defined page-length of less than 60 lines.
 If a subprogram call is needed to accomplish this, it may be convenient to
 place it in subprogram VERSP, which is executed before each page of STARPAC
 output is produced.  This problem is known to affect at least some Sperry
 systems.


 A.5  STARPAC Compatibility with Other Software

 There are only two known sources of incompatibility with other software.
 First, the combination of STARPAC with another code might be too large to run
 on a particular system.  Second, an overlap of the names of STARPAC
 subprograms or labeled commons with the names of such entities in other code
 might prevent STARPAC being run in conjunction with the other code [see
 STARPAC User's Guide, appendices D and E, for a list of STARPAC 2.08
 subprogram and common names].  Some systems may have loaders which eliminate
 the second difficulty in one way or another.  In the case of overlapping
 commons, no problem results if the commons' dimensions are compatible and the
 other code does not expect the common to be intact during or after a STARPAC
 call.

 A particular instance of the overlapping common names problem occurs because
 STARPAC uses the Framework for a Portable Library [Fox et al., 1978a,b] in a
 modified version that reinitializes the Framework's CSTAK common each time a
 STARPAC user-callable subprogram is called.  The published version of the
 Framework initializes CSTAK once only - the first time that a Framework
 subprogram is called.  The requirement for what, in effect, are two different
 CSTAK commons may prevent the use of STARPAC with other Framework-using
 software, such as the PORT library [Fox et al., 1978a], since calling STARPAC
 would destroy the contents of PORT's CSTAK.  See Donaldson and Tryon [1987:
 Chapter 1, Section D.2] for further discussion.


                      Appendix B:  Comparing STARPAC Test Results

 Differences between STARPAC outputs run on different machines can occur for a
 variety of reasons.  The most pervasive of the differences between the
 installers' outputs and those on the tape are apt to be those resulting from
 differences in the implementation of the F and G conversion specifications in
 different Fortrans.  Some Fortrans, for example, print '0.xxxx' when others
 print ' .xxxx', and so on.  This kind of difference is obviously
 insignificant.  Also, there will be differences observed in the format of
 many of the values printed in exponential notation if the double precision
 version of STARPAC is installed, since file TEST.OUT contains output from a
 single precision version of STARPAC.  This difference in precision will also
 be reflected in the version number printed at the top of each page of STARPAC
 output:  the single precision version has a version number followed by the
 letter S (e.g., 2.08S) while the double precision version has a version
 number followed by the letter D (e.g., 2.08D).

 Next, bear in mind that either of the files could legitimately contain
 material different from that in the other.  To begin with, certain values
 included in the output are machine dependent.  These include the code printed
 to indicate an uncomputable result, and the default stopping criteria used in
 such iterative processes as nonlinear least squares.  The uncomputed value
 code is the largest legal rational magnitude; in the outputs on the tape,
 which were generated on CDC equipment, this code can be recognized as a value
 of approximately 1.2E+322.  Default stopping criteria are multiples of some
 root of the machine precision (largest relative spacing), and are labeled as
 stopping criteria in the output.

 Differences in machine precision will be reflected in the outputs in a number
 of additional ways.  Trivially, since the computed values are only
 approximately correct, performing the calculations in two different
 precisions can result in approximations which differ in the last few decimal
 places.  In the case of values which are actually approximations to zero, the
 values may differ in all places, and even in scale (approximations to zero
 values can be recognized as very small fractions which may differ between the
 two test outputs).  Results from the random number generator will also differ
 in the last few places, causing corresponding differences in the results from
 STARPAC procedures such as HIST, which are used to analyze and display the
 generated values.  Differences in machine precision can also cause values to
 fall into different intervals in the plotting routines, and the frequency
 distribution included in the STAT family of procedures.

 Installers should expect differences in machine precision and consequent
 differences in stopping criteria to cause iterative processes in different
 versions of STARPAC to take different numbers of iterations to solve the same
 problem, changing the value of the reported number of iterations.

 Another source of different answers is the nonassociativity of machine
 arithmetic.  Different orders of evaluation of the same expression on the
 same machine can produce different answers, with symptoms similar to those
 produced by a slight change in precision.  Though Fortran compilers are
 deterministic, there are a number of ways in which alternate orders of
 evaluation can be selected on the same machine.  It can happen on the same
 machine, but under different compilers, or under the same compiler, but at
 different update or optimization levels.  Installers with CDC equipment who
 find that their results differ from those on the tape should bear this in
 mind, as should all installers who discover that the answers change in
 response to a simple change in the choice of compiler or compiler option, or
 to recompilation with a new version of a compiler.  These differences can
 seem alarming if they are unexpected, but are actually trivial.

 Some test case problems which are relatively more sensitive to differences in
 precision may become unsolvable when precision is reduced.  In such cases, an
 error report may be substituted for a more lengthy solution report.  A
 certain amount of this sort of change is to be expected; however, if a large
 proportion of the test problems become unsolvable under the installed version
 of STARPAC, there is a serious problem.  The most likely scenario for this is
 an injudicious attempt to install STARPAC in single precision on a machine
 that has only six or seven decimal digits at that precision.  In this case,
 the single precison version should be abandoned in favor of the double
 precision version.

 One final source of differences in test outputs is that runtime error or
 status messages may occur in one of the test outputs but not the other, due
 to differences in the systems, compilers, or hardware.


                      Appendix C: Systems Running STARPAC

 Equipment             Operating        Compiler                    Precision
                       System

 CDC Cyber 180/855     NOS 2.5.2        FTN 5.1+670 (opt = 0 & 1)   single

 CDC Cyber 200/205     VSOS 2.3         FTN200                      single

 Concurrent Computer
 Corportation 3230     OS32MT           FORTRAN-VIIZ                double

 VAX 11/780            VMS 2.2          FORTRAN 4.7                 double

 SUN 3/180             BSD UNIX 4.3     F77                         double


                                   References

 American National Standards Institute. (1966).  U.S.A. standard Fortran.  New
 York: United States of America Standards Institute.

 American National Standards Institute. (1978).  American national standard
 programming language FORTRAN - approved April 3, 1978.  New York: American
 National Standards Institute.

 Dennis, J. E. Jr.; Gay, D. M.; Welsch, R. E. (1981a).  An adaptive nonlinear
 least-squares algorithm.  ACM Trans. on Math. Software, 7(3): 348-368.

 Dennis, J.  E. Jr.; Gay, D. M.; Welsch, R.  E. (1981b).  Algorithm 573: NL2SOL
 - an adaptive nonlinear least-squares algorithm.  ACM Trans. on Math.
 Software, 7(3): 369-383.

 Donaldson, J.  R.; Tryon, P. V. (1983a).  Introduction to STARPAC, the
 Standards time series and regression package.  National Institute of Standards
 and Technology Technical Note NBSTN 1068-1.

 Donaldson, J. R.; Tryon, P. V. (1983b).  Nonlinear least squares regression
 using STARPAC, the Standards time series and regression package.  National
 Institute of Standards and Technology Technical Note NBSTN 1068-2.

 Donaldson, J.  R.; Tryon, P.  V.  (1987).  STARPAC, the Standards Time Series
 and Regression Package.  National Institute of Standards and Technology
 Internal Report NBSIR 86-3448.

 Dongarra, J. J.; Moler, C. B.; Bunch, J. R.; Stewart, G. W. (1980).  LINPACK
 users' guide.  Philadelphia: SIAM.

 Dorrenbacher, J.; Paddock, D.; Wisneski, D.; Fosdick, L. D. (1976).  POLISH, a
 Fortran program to edit Fortran programs.  Dept. of Computer Science,
 University of Colorado, Boulder, Report CU-CS-050-76 (revised) May 1976.

 Filliben, J. J. (1977).  User's guide to the DATAPAC data analysis package,
 Version 77.5.  Draft.

 Fox, P. A.; Hall, A. D.; and Schryer, N. L. (1978a).  The PORT mathematical
 subroutine library.  ACM Trans. on Math.  Software, 4(2): 104-126.

 Fox, P.  A.; Hall, A.  D.; and Schryer, N.  L. (1978b).  Algorithm 528:
 framework for a portable library.  ACM Trans. on Math.  Software, 4(2): 177-
 188.

 Fullerton, L.  W.  (1977).  Portable special function routines.  In
 Portability of Numerical Software:  Proceedings.  W.  Crowell, editor.
 (Lecture Notes in Computer Science:  Vol.  57).  Oak Brook, IL:  Springer-Verlag.

 Hogben, D.; Peavy, S. T.; Varner, R. N. (1971).  OMNITAB II user's reference
 manual.  National Institute of Standards and Technology Technical Note NBSTN
 552.

 IMSL. (1982).  IMSL library manual.  Houston: IMSL, Inc.

 Lawson, F. L.; Hanson, R. J.; Kincaid, D.  R.; Krogh, F. T. (1979a).  Basic
 linear algebra subprograms for Fortran usage.  ACM Trans. on Math. Software,
 5(3): 308-323.

 Lawson, F.  L.; Hanson, R. J.; Kincaid, D.  R.; Krogh, F.  T. (1979b).
 Algorithm 539: basic linear algebra subprograms for Fortran usage.  ACM Trans.
 on Math. Software, 5(3): 324-325.

 Ryder, G. G.; Hall, A. D. (1973).  The PFORT verifier.  Computing Science
 Technical Report No. 12, (Revised April 1979.) Murray Hill: Bell Laboratories.

 Siemieniuch, J. L. (n.d.).  APT.  Internal document and source code developed
 by the Numerical Algorithms Group Ltd.

 Tryon, P.  V.; Donaldson, J. R. (1978).  STATLIB - a library of FORTRAN
 subroutines for statistical analysis of experimental data.  Revised June 1,
 1978 and October 2, 1978.  Unpublished.

1
                                                                     (03/16/87)
              STARPAC 2.08 Release Installation Registration Form

 Filling out and returning this form will ensure receipt of all future updates
 to the source code and manual, and of any announcements of new versions.

 We welcome your comments on STARPAC, the STARPAC release, and STARPAC
 portability.  Please feel free to include these on separate sheets.  We are
 particularly interested in any unanticipated changes that you may have had to
 make to the code in order to install STARPAC.


 1.   Contact

         Name:  .........................................................
      Address:  .........................................................
                .........................................................
                .........................................................
                .........................................................
                                                          Zip:  .........
      Telephone:  ....................

 2.   Installer (optional, may report if different from contact)

         Name:  .........................................................
      Address:  .........................................................
                .........................................................
                .........................................................
                .........................................................
                                                          Zip:  .........

 3.   Name of Installation Site (if not clear from address)

      ...................................................................

 4.   STARPAC Operating Environment (use additional sheets as necessary)

                                   Machine 1               Machine 2

             Computer Model:  .....................   ...................
           Operating System:  .....................   ...................
           Fortran Compiler:  .....................   ...................



 Please return to:

       Janet R. Donaldson
       National Institute of Standards and Technology
       Mail Code 719
       325 Broadway
       Boulder, CO  80303-3328

       303-497-5114 or FTS 320-5114
       Internet:  jrd@cam.nist.gov;  Bitnet:  jrd@nistcs2


-->
<!--
toc.doc
 STARPAC 2.08 -- The Standards Time Series and Regression Package


 Direct questions to

      Janet R. Donaldson
      Optimization Group/Applied and Computational Mathematics Division (719)
      National Institute of Standards and Technology
      325 Broadway
      Boulder, CO 80303-3328
      (303) 497-5114
      e-mail:  internet -- jrd@cam.nist.gov
               bitnet   -- jrd@nistcs2


 *** PHYSICAL CHARACTERISTICS OF TAPE

      A.  ASCII character set.
      B.  1600 cpi density.
      C.  Unlabeled.
      D.  11 files, each terminated by tapemarks.
      E.  Additional tapemark follows tapemark of last file.
      F.  Files consist of 1 or more blocks (physical records).
      G.  Files 1 to 9 have blocks of
          45 line images (logical records) of 80 characters each,
          i.e., 3600 characters;
      H.  Files 10 and 11 have blocks of
          20 line images (logical records) of 132 characters each,
          i.e., 2640 characters;
      I.  Last block of a file may contain fewer than the specified number of
          line images, in which case it is short, not blank filled.



 *** TAPE CONTENTS

 File No.  File Id.       Description
 --------  -----------    -----------

        1  TOC.DOC      - tape characteristics, file structure and
                          table of contents
                          (line image length = 80, block size = 3600)

        2  INSTALL.DOC  - this installation manual
                          (line image length = 80, block size = 3600)

        3  S_SRCE.FOR   - single precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (line image length = 80, block size = 3600)

        4  S_NL2SOL.FOR - single precision STARPAC 2.08 source code from NL2SOL
                          (line image length = 80, block size = 3600)

        5  D_SRCE.FOR   - double precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (line image length = 80, block size = 3600)

        6  D_NL2SOL.FOR - double precision STARPAC 2.08 source code from NL2SOL
                          (line image length = 80, block size = 3600)

        7  MISC.FOR     - single and double precision subprograms from
                          miscellaneous public domain sources, but excluding
                          code from NL2SOL, LINPACK and BLAS and excluding
                          subprograms supplying machine dependent constants
                          (line image length = 80, block size = 3600)

        8  LPK_BLS.FOR  - single and double precision subprograms
                          from LINPACK and BLAS
                          (line image length = 80, block size = 3600)

        9  MDC.FOR      - single and double precision subprograms supplying
                          machine dependent constants
                          (line image length = 80, block size = 3600)

       10  TEST.OUT     - results of test subprograms
                          (line image length = 132, block size = 2640)

       11  GUIDE.DOC    - STARPAC 2.08 user documentation,
                          in line printer format
                          (line image length = 132, block size = 2640)

-->
<!--
                               STARPAC 2.08
              The Standards Time Series and Regression Package


 STARPAC 2.08 is a library of portable Fortran subprograms for conducting
 statistical analysis.  It was developed and is maintained by the Center for
 Computing and Applied Mathematics of the National Institute of Standards and
 Technology (formerly the National Bureau of Standards).  STARPAC contains
 subprogram groups for performing time series analysis and nonlinear least
 squares regression, as well as normal random number generation, line printer
 plots, basic statistical analyses, and linear least squares.

 Currently, the primary mode of distribution for STARPAC is magnetic tape.
 The tape release package consists of the 11 files contained in this directory:

        1  TOC.DOC      - tape characteristics, file structure and
                          table of contents 
                          (3318 bytes)

        2  INSTALL.DOC  - the installation manual
                          (59735 bytes)

        3  S_SRCE.FOR   - single precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (2575385 bytes)

        4  S_NL2SOL.FOR - single precision STARPAC 2.08 source code from NL2SOL
                          (194141 bytes)

        5  D_SRCE.FOR   - double precision STARPAC 2.08 source code, excluding
                          code from various public domain sources
                          (2628509 bytes)

        6  D_NL2SOL.FOR - double precision STARPAC 2.08 source code from NL2SOL
                          (195881 bytes)

        7  MISC.FOR     - single and double precision subprograms from
                          miscellaneous public domain sources, but excluding
                          code from NL2SOL, LINPACK and BLAS and excluding
                          subprograms supplying machine dependent constants
                          (118434 bytes)

        8  LPK_BLS.FOR  - single and double precision subprograms
                          from LINPACK and BLAS
                          (76661 bytes)

        9  MDC.FOR      - single and double precision subprograms supplying
                          machine dependent constants
                          (47665 bytes)

       10  TEST.OUT     - results of test subprograms
                          (349444 bytes)

       11  GUIDE.DOC    - STARPAC 2.08 user documentation, 
                          in line printer format
                          (847693 bytes)

 File TOC.DOC is obviously not needed when the files are acquired via FTP.
 File INSTALL.DOC contains a more detailed description of the STARPAC files and
 of the installation process; installers may wish to read it before
 transmitting the remaining files.

 Before proceeding to acquire the STARPAC files, installers should note that
 although STARPAC is supplied in both single and double precision versions, the
 two versions cannot be used together.  This is because the names of most of
 the subprograms are the same in the two versions.  As is documented in the
 Installation Guide, most sites will want to install only one of the two
 versions.  Transmission time can therefore be reduced if the source files for
 only one version are transmitted.

   *  For the single precision version, the source files
          S_SRCE.FOR, S_NL2SOL.FOR, MISC.FOR, LPK_BLS.FOR and MDC.FOR 
      should be used.

   *  For the double precision version, the source files 
          D_SRCE.FOR, D_NL2SOL.FOR, MISC.FOR, LPK_BLS.FOR and MDC.FOR 
      should be used.  

 Files INSTALL.DOC, GUIDE.DOC and TEST.OUT will be needed by all sites.

 It is particularly important that the STARPAC developers be notified that you
 have acquired STARPAC when the package is obtained via FTP.  Without such
 notification, we will not know that your site should receive future notices
 and updates.  Notification can be accomplished either by filling out and
 returning the registration form at the end of the Installation Guide, or by
 sending an e-mail message containing the recipient's name, address, phone
 number and e-mail address to

       Janet R. Donaldson

          Internet:  jrogers@boulder.nist.gov
   
          Address:   National Institute of Standards and Technology
                     Applied and Computational Mathematics Division
                     Mail Code 719
                     325 Broadway
                     Boulder, Colorado  80303-3328

          Phone:     303-497-5114 or FTS 320-5114
   

 We welcome your comments on STARPAC, the STARPAC release, and STARPAC
 portability.
-->
 !1000 FORMAT (105X, 'STARPAC 2.08 (03/15/90)')
 !1010 FORMAT (54X, 'STARPAC 2.08 (03/15/90)')

### MISMATCHES
Calls to PPLMT or PPMN many times have fourth argument as an array when it is a scalar in code.
For now, pass first element of arrays.
```fortran
      CALL PPLMT(YAXIS,YAXIS,XAXIS,XAXIS(1),2*NF,1,2*NF,-2*PI,2*PI,YMN,YMX,XPLTMN,XPLTMX,XMN,XMX,ERROR,NMSUB,.FALSE.)
      CALL PPMN(YAXIS, YAXIS, XAXIS, XAXIS(1),
```
The call to PPMN:
```text
srce/ppmn.ffinc:
!!SUBROUTINE PPMN (YM(IYM,M),YMMISS(M),X(N),xmiss,n,m,iym,ischck,ISYM(LISYM),lisym,isize,nout,ymn,ymx,xmn,xmx,_miss,ilog)
```

No clarity or code consistency and not clearly documented if YMISS is an array or a single value in a lot of the plotting
procedures, and in which routines the Y data is a vector or a matrix. Might need to make some generics.

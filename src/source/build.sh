#!/bin/bash
cd $(dirname $0)
export GITHUB=TRUE
export DEMO_OUTDIR=../../example
export DEMO_SUBDIR=FALSE
GPF_build_module M_starpac M_starpac_s M_starpac_d
cp ../../docs/man3.html ../../docs/index.html
#cp ../../docs/BOOK_M_starpac.html ../../docs/index.html
#ccall ../../test/test_suite_M_starpac.[fF]90
changen -cmd 'mv' .f90 .f -- ../*.f90
exit

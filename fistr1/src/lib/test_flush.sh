#!/bin/sh
#
# File        : test_flush.sh
# Author      : K. Goto (PExProCS, LLC.)
# Last Update : Feb. 14, 2014
#
# Description :
#  Test availability of subroutine flush or flush statement with given
#  F90 compiler. If needed, generate subroutine flush either with flush
#  statement or empty statement, and add the object to target fortran
#  library.
#
#  Following environment variables has to be pre-defined
#  (if not defined, default value in () is used).
#
#  F90       : fortran compiler command ("mpif90")
#  F90FLAGS  : fortran compiler options ("")
#  F90TARGET : fortran library ("../../lib/libffistr.a")
#  AR        : archive command for building static library ("ar ruv")
#  RM        : remove command ("rm -f")
#
logfile=test_flush.log

if [ "${F90}" = "" ]; then
    echo "WARNING: F90 not defined" >> $logfile
    F90="mpif90"
fi
if [ "${F90FLAGS}" = "" ]; then
    echo "WARNING: F90FLAGS not defined" >> $logfile
    F90FLAGS=""
fi
if [ "${F90TARGET}" = "" ]; then
    echo "WARNING: F90TARGET not defined" >> $logfile
    F90TARGET="../../lib/libffistr.a"
fi
if [ "${AR}" = "" ]; then
    echo "WARNING: AR not defined" >> $logfile
    AR="ar ruv"
fi
if [ "${RM}" = "" ]; then
    echo "WARNING: RM not defined" >> $logfile
    RM="rm -f"
fi

${RM} ${logfile}

HAVE_FLUSH_SUB=0
HAVE_FLUSH_STAT=0

run_command() {
    cmd=$*
    echo + ${cmd} >> ${logfile}
    ${cmd} 1>> ${logfile} 2>&1
}

test_f90_flush() {
    ${RM} $1
    run_command ${F90} -o $1 $1.f90
    if [ -f $1 ]; then
        # run_command ./$1
        # out=`./$1`
        ${RM} $1
        # if [ "${out}" = "test_flush" ]; then
        return 0
        # fi
    fi
    return 1
}

cat > test_flush_sub.f90 <<EOF
program test_flush
  implicit none
  write(6,'(A)') 'test_flush'
  call flush(6)
end program test_flush
EOF
run_command cat test_flush_sub.f90
test_f90_flush test_flush_sub && HAVE_FLUSH_SUB=1
echo + HAVE_FLUSH_SUB=$HAVE_FLUSH_SUB >> ${logfile}
${RM} test_flush_sub.f90

if [ $HAVE_FLUSH_SUB -eq 1 ]; then
    exit 0
fi

cat > test_flush_stat.f90 <<EOF
program test_flush
  implicit none
  write(6,'(A)') 'test_flush'
  flush(6)
end program test_flush
EOF
run_command cat test_flush_stat.f90
test_f90_flush test_flush_stat && HAVE_FLUSH_STAT=1
echo + HAVE_FLUSH_STAT=$HAVE_FLUSH_STAT >> ${logfile}
${RM} test_flush_stat.f90

if [ $HAVE_FLUSH_STAT -eq 1 ]; then
    cat > flush_stat.f90 <<EOF
subroutine flush(n)
  implicit none
  integer :: n
  flush(n)
end subroutine flush
EOF
    run_command cat flush_stat.f90
    run_command ${F90} -c flush_stat.f90
    run_command ${AR} ${F90TARGET} flush_stat.o
    ${RM} flush_stat.*
else
    cat > flush_empty.f90 <<EOF
subroutine flush(n)
  implicit none
  integer :: n
end subroutine flush
EOF
    run_command cat flush_empty.f90
    run_command ${F90} -c flush_empty.f90
    run_command ${AR} ${F90TARGET} flush_empty.o
    ${RM} flush_empty.*
fi

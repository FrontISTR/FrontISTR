###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(Fistr_Test)

find_program(BASH_PROGRAM bash)

function(GetProcThreads type np nt)
  if(type STREQUAL "serial")
    set(np 1 PARENT_SCOPE)
    set(nt 1 PARENT_SCOPE)
  elseif(type STREQUAL "openmp")
    set(np 1 PARENT_SCOPE)
    set(nt 2 PARENT_SCOPE)
  elseif(type STREQUAL "mpi")
    set(np 2 PARENT_SCOPE)
    set(nt 1 PARENT_SCOPE)
  elseif(type STREQUAL "hybrid")
    set(np 2 PARENT_SCOPE)
    set(nt 2 PARENT_SCOPE)
  endif()
endfunction(GetProcThreads)
set(parallels serial)
IF(WITH_MPI)
  list(APPEND parallels mpi)
endif()
IF(WITH_OPENMP)
  list(APPEND parallels openmp)
endif()
IF(WITH_OPENMP AND WITH_MPI)
  list(APPEND parallels hybrid)
endif()
foreach(parallel IN LISTS parallels)
  GetProcThreads(${parallel} np nt)

  foreach(target analysis/static analysis/eigen analysis/heat analysis/dynamic solver/CG)
    add_test(NAME Test-${parallel}-${target}
     COMMAND ./test.sh -f ${FrontISTR_BINARY_DIR}/fistr1/fistr1 -e ${hecmw_BINARY_DIR}/tools/hecmw_part1 -r ${hecmw_BINARY_DIR}/tools/rmerge -d ${target} -p ${np} -t ${nt} -vv
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test-${parallel}-${target} PROPERTIES LABELS ${target},${parallel})
  endforeach()
  if(WITH_MUMPS)
    set(target solver/mumps)
    add_test(NAME Test-${parallel}-${target}
     COMMAND ./test.sh -f ${FrontISTR_BINARY_DIR}/fistr1/fistr1 -e ${hecmw_BINARY_DIR}/tools/hecmw_part1 -r ${hecmw_BINARY_DIR}/tools/rmerge -d ${target} -p ${np} -t ${nt} -vv
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test-${parallel}-${target} PROPERTIES LABELS ${target},${parallel})
  endif()
  if(WITH_MKL)
    set(target solver/mkl)
    add_test(NAME Test-${parallel}-${target}
     COMMAND ./test.sh -f ${FrontISTR_BINARY_DIR}/fistr1/fistr1 -e ${hecmw_BINARY_DIR}/tools/hecmw_part1 -r ${hecmw_BINARY_DIR}/tools/rmerge -d ${target} -p ${np} -t ${nt} -vv
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test-${parallel}-${target} PROPERTIES LABELS ${target},${parallel})
  endif()
endforeach()




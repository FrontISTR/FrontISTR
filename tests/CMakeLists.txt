###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

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

# Get the list of tests in specifed directory
function(AppendListTests rootdir list_in TEST_TARGETS)

  set(TMPLIST "${list_in}")
  file(GLOB_RECURSE TEST_ANALYSIS_MESH_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${rootdir}/*.msh")
  foreach(MESH_FILE ${TEST_ANALYSIS_MESH_FILES})
    get_filename_component(DIR_NAME ${MESH_FILE} DIRECTORY)
    list(APPEND TMPLIST ${DIR_NAME})
  endforeach()
  list(REMOVE_DUPLICATES TMPLIST)

  set(TEST_TARGETS "${TMPLIST}" PARENT_SCOPE)
endfunction(AppendListTests)

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

  # Get the default list of tests
  set(TEST_TARGETS "")
  AppendListTests("analysis/*" "${TEST_TARGETS}" TEST_TARGETS)
  AppendListTests("lib/*" "${TEST_TARGETS}" TEST_TARGETS)
  AppendListTests("solver/*" "${TEST_TARGETS}" TEST_TARGETS)

  # Get the optional list of tests
  if(WITH_MUMPS)
    AppendListTests("with_mumps/*" "${TEST_TARGETS}" TEST_TARGETS)
  endif()
  if(WITH_MKL AND ( (parallel STREQUAL "serial") OR (parallel STREQUAL "openmp") ) )
    AppendListTests("with_mkl/*" "${TEST_TARGETS}" TEST_TARGETS)
  endif()
  if(WITH_ML)
    AppendListTests("with_ml/*" "${TEST_TARGETS}" TEST_TARGETS)
  endif()

  foreach(target ${TEST_TARGETS})
    get_filename_component(LAST_TEST_NAME ${target} NAME)
    if(NOT (${LAST_TEST_NAME} STREQUAL WITH_MPI))
      set(TEST_NAME "test_${parallel}_${target}")
      add_test(NAME ${TEST_NAME}
        COMMAND
          ./test.sh
            -f ${FrontISTR_BINARY_DIR}/fistr1/fistr1
            -e ${hecmw_BINARY_DIR}/tools/hecmw_part1
            -r ${hecmw_BINARY_DIR}/tools/rmerge
            -d ${target}
            -p ${np}
            -t ${nt}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
      set_tests_properties(
        ${TEST_NAME}
        PROPERTIES
          LABELS ${target},${parallel}
          ENVIRONMENT CTEST_TEST_NAME=${TEST_NAME}
      )
      unset(TEST_NAME)
    endif()
  endforeach()
endforeach()

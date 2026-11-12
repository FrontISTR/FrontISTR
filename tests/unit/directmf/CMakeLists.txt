###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# The tests use the modules of hecmw, which are written to its binary directory.
set(directmf_unit_names symbolic numeric pivot lu blr)
if(WITH_MPI AND MPIEXEC_EXECUTABLE)
  list(APPEND directmf_unit_names mpi)
endif()
foreach(name IN LISTS directmf_unit_names)
  set(exe directmf_unit_${name})
  add_executable(${exe} test_mf_${name}.f90)
  set_target_properties(${exe} PROPERTIES LINKER_LANGUAGE "CXX")
  target_link_libraries(${exe} PRIVATE hecmw)
  target_include_directories(${exe} PRIVATE ${FrontISTR_INCLUDE_DIRS} ${hecmw_BINARY_DIR})
  target_compile_definitions(${exe} PRIVATE ${FrontISTR_DEFINITIONS})
endforeach()

# The single process tests are registered in the modes without MPI.
set(directmf_unit_parallels serial)
if(WITH_OPENMP)
  list(APPEND directmf_unit_parallels openmp)
endif()
foreach(name symbolic numeric pivot lu blr)
  foreach(parallel IN LISTS directmf_unit_parallels)
    GetProcThreads(${parallel} np nt)
    set(TEST_NAME "test_${parallel}_unit/directmf/${name}")
    add_test(NAME ${TEST_NAME} COMMAND directmf_unit_${name})
    set_tests_properties(${TEST_NAME} PROPERTIES
      LABELS unit/directmf,${parallel}
      ENVIRONMENT OMP_NUM_THREADS=${nt})
  endforeach()
endforeach()

if(NOT "mpi" IN_LIST directmf_unit_names)
  return()
endif()

# test_mf_mpi is registered in every mode and started as test.sh starts fistr1: directly
# on a single process, otherwise through mpiexec with the options test.sh gives Open MPI.
set(directmf_unit_mpiexec_flags ${MPIEXEC_PREFLAGS})
set(directmf_unit_mpi_env)
if("${MPI_C_LIBRARY_VERSION_STRING};${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI")
  list(APPEND directmf_unit_mpiexec_flags --oversubscribe --allow-run-as-root)
  list(APPEND directmf_unit_mpi_env OMPI_MCA_btl_sm_backing_directory=/tmp)
endif()
foreach(parallel IN LISTS parallels)
  GetProcThreads(${parallel} np nt)
  set(TEST_NAME "test_${parallel}_unit/directmf/mpi")
  if(np EQUAL 1)
    add_test(NAME ${TEST_NAME} COMMAND directmf_unit_mpi)
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
      ${directmf_unit_mpiexec_flags} $<TARGET_FILE:directmf_unit_mpi> ${MPIEXEC_POSTFLAGS})
  endif()
  set(env OMP_NUM_THREADS=${nt} ${directmf_unit_mpi_env})
  set_tests_properties(${TEST_NAME} PROPERTIES
    LABELS unit/directmf,${parallel}
    ENVIRONMENT "${env}")
endforeach()

# Only 4 or more processes put the fully summed part of a front on a process grid with two
# columns. PROCESSORS keeps ctest -j from running this test on the cores of other tests.
set(TEST_NAME "test_mpi_unit/directmf/mpi_np4")
add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
  ${directmf_unit_mpiexec_flags} $<TARGET_FILE:directmf_unit_mpi> ${MPIEXEC_POSTFLAGS})
set(env OMP_NUM_THREADS=1 ${directmf_unit_mpi_env})
set_tests_properties(${TEST_NAME} PROPERTIES
  LABELS unit/directmf,mpi
  ENVIRONMENT "${env}"
  PROCESSORS 4)

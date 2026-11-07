###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# The SA-AMG core modules are compiled a second time here without HEC-MW
# (HECMW_SERIAL, HECMW_SAAMG_STANDALONE), as an object library with its own
# module directory so that its modules do not clash with those of hecmw.
set(SAAMG_DIR ${FrontISTR_SOURCE_DIR}/hecmw1/src/solver/precond/saamg)

add_library(saamg_unit_core OBJECT
  ${SAAMG_DIR}/hecmw_precond_SAAMG_util.f90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_param.f90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_comm.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_matrix.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_smoother.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_aggregate.f90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_prolongation.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_coarse.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_coarse_mumps.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_coarse_mf.F90
  ${SAAMG_DIR}/hecmw_precond_SAAMG_core.f90
)
target_compile_definitions(saamg_unit_core PRIVATE HECMW_SERIAL HECMW_WITH_LAPACK HECMW_SAAMG_STANDALONE)
set_target_properties(saamg_unit_core PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

# Registered in the single-process modes only: the core is built here without
# MPI, so the mpi and hybrid modes would check nothing more. Each mode runs in
# its own directory because test_step8 writes a file and reads it back.
set(saamg_unit_parallels serial)
if(WITH_OPENMP)
  list(APPEND saamg_unit_parallels openmp)
endif()
foreach(parallel IN LISTS saamg_unit_parallels)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${parallel})
endforeach()

foreach(n RANGE 1 10)
  set(exe saamg_unit_step${n})
  add_executable(${exe} test_step${n}.f90 $<TARGET_OBJECTS:saamg_unit_core>)
  target_include_directories(${exe} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mod)
  target_link_libraries(${exe} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  add_dependencies(${exe} saamg_unit_core)
  foreach(parallel IN LISTS saamg_unit_parallels)
    GetProcThreads(${parallel} np nt)
    set(TEST_NAME "test_${parallel}_unit/saamg/step${n}")
    add_test(NAME ${TEST_NAME} COMMAND ${exe}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${parallel})
    set_tests_properties(${TEST_NAME} PROPERTIES
      LABELS unit/saamg,${parallel}
      ENVIRONMENT OMP_NUM_THREADS=${nt})
  endforeach()
endforeach()

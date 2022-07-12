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

  foreach(case exA exB exC exD exE exF exG exI spring_boundary)
    set(type static)
    add_test(NAME Test_${type}_${case}_${parallel}
     COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} ${type}/${case} ${np} ${nt}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test_${type}_${case}_${parallel} PROPERTIES LABELS ${type},${case},${parallel})
  endforeach()  

  foreach(case exJ exK)
    set(type eigen)
    add_test(NAME Test_${type}_${case}_${parallel}
     COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} ${type}/${case} ${np} ${nt}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test_${type}_${case}_${parallel} PROPERTIES LABELS ${type},${case},${parallel})
  endforeach()  

  foreach(case exM exN exO exP exQ exR exS exT exU exU2 exV)
    set(type heat)
    add_test(NAME Test_${type}_${case}_${parallel}
     COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} ${type}/${case} ${np} ${nt}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test_${type}_${case}_${parallel} PROPERTIES LABELS ${type},${case},${parallel})
  endforeach()  

  foreach(case exW exX)
    set(type dynamic)
    add_test(NAME Test_${type}_${case}_${parallel}
     COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} ${type}/${case} ${np} ${nt}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(Test_${type}_${case}_${parallel} PROPERTIES LABELS ${type},${case},${parallel})
  endforeach()  

endforeach()




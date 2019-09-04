###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(Fistr_Test)

find_program(BASH_PROGRAM bash)

foreach(case exA exB exC exD exE exF exG exI spring_boundary)
  add_test(
    NAME Static_${case}_Test
   COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} static/${case}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_tests_properties(Static_${case}_Test PROPERTIES LABELS static)
endforeach()

foreach(case exJ exK)
  add_test(
    NAME Eigen_${case}_Test
   COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} eigen/${case}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_tests_properties(Eigen_${case}_Test PROPERTIES LABELS eigen)
endforeach()

foreach(case exM exN exO exP exQ exR exS exT exU exV)
  add_test(
   NAME Heat_${case}_Test
   COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} heat/${case}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_tests_properties(Heat_${case}_Test PROPERTIES LABELS heat)
endforeach()

foreach(case exW exX)
  add_test(
   NAME Dynamic_${case}_Test
   COMMAND ruby ./test_FrontISTR.rb ${FrontISTR_BINARY_DIR} dynamic/${case}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_tests_properties(Dynamic_${case}_Test PROPERTIES LABELS dynamic)
endforeach()

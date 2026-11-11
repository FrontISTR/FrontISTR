###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_executable(ctrl_data_line_test ctrl_data_line_test.c)
target_include_directories(ctrl_data_line_test PRIVATE
  ${FrontISTR_SOURCE_DIR}/fistr1/src/lib/utilities)
target_link_libraries(ctrl_data_line_test PRIVATE fistr)
add_test(NAME test_serial_unit/input/data_line COMMAND ctrl_data_line_test)
set_tests_properties(test_serial_unit/input/data_line PROPERTIES
  LABELS unit/input,serial)

add_executable(solver_duplicate_test solver_duplicate_test.f90)
target_include_directories(solver_duplicate_test PRIVATE
  ${FrontISTR_BINARY_DIR}/fistr1)
target_link_libraries(solver_duplicate_test PRIVATE fistr)
foreach(format cnt inp)
  add_test(NAME test_serial_unit/input/solver_${format}
    COMMAND ${CMAKE_COMMAND}
      -DPROGRAM=$<TARGET_FILE:solver_duplicate_test>
      -DFORMAT=${format}
      -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/check_solver.cmake)
  set_tests_properties(test_serial_unit/input/solver_${format} PROPERTIES
    LABELS unit/input,serial)
endforeach()

add_executable(ctrl_file_query_test ctrl_file_query_test.c)
target_link_libraries(ctrl_file_query_test PRIVATE hecmw)
add_test(NAME test_serial_unit/input/control_file
  COMMAND ctrl_file_query_test ${CMAKE_CURRENT_BINARY_DIR}/ctrl_query.inp)
set_tests_properties(test_serial_unit/input/control_file PROPERTIES
  LABELS unit/input,serial)

if(WITH_METIS AND METIS_FOUND)
  add_test(NAME test_serial_unit/input/partition_messages
  COMMAND ${CMAKE_COMMAND}
    -DPROGRAM=$<TARGET_FILE:hecmw_part1>
    -DINPUT=${FrontISTR_SOURCE_DIR}/tests/analysis/inp/standalone_static/standalone_static.inp
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/partition_messages
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_partition.cmake)
  set_tests_properties(test_serial_unit/input/partition_messages PROPERTIES
    LABELS unit/input,serial)
endif()

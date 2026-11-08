###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

project(FrontISTR_API_Test LANGUAGES C)

enable_testing()
include_directories(../../fistr1/src ../../hecmw1/src)

add_executable(hecmw_mesh_api_test hecmw_mesh_api_test.c)
target_link_libraries(hecmw_mesh_api_test PRIVATE hecmw)
add_test(NAME test_serial_api/mesh COMMAND hecmw_mesh_api_test)
set_tests_properties(test_serial_api/mesh
  PROPERTIES
    LABELS api/mesh,serial
    ENVIRONMENT CTEST_TEST_NAME=test_serial_api/mesh
)

add_executable(hecmw_result_api_test hecmw_result_api_test.c)
target_link_libraries(hecmw_result_api_test PRIVATE hecmw)
add_test(NAME test_serial_api/result COMMAND hecmw_result_api_test)
set_tests_properties(test_serial_api/result
  PROPERTIES
    LABELS api/result,serial
    ENVIRONMENT CTEST_TEST_NAME=test_serial_api/result
)

add_executable(fistr_param_api_test fistr_param_api_test.c)
target_link_libraries(fistr_param_api_test PRIVATE fistr hecmw)
add_test(NAME test_serial_api/param COMMAND fistr_param_api_test)
set_tests_properties(test_serial_api/param
  PROPERTIES
    LABELS api/param,serial
    ENVIRONMENT CTEST_TEST_NAME=test_serial_api/param
)


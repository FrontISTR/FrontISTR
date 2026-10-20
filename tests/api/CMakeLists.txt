###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

project(FrontISTR_API_Test LANGUAGES C)

enable_testing()
include_directories(../../fistr1/src ../../hecmw1/src)

add_executable(fistr_param_api_test fistr_param_api_test.c)
target_link_libraries(fistr_param_api_test PRIVATE fistr hecmw)
add_test(NAME test_param_api COMMAND fistr_param_api_test)
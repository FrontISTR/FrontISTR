###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

# neu2fstr
set(neu2fstr_SOURCES)
set(neu2fstr_INCLUDE_DIRS)
include(neu2fstr/CMakeLists.txt)
add_executable(neu2fstr ${neu2fstr_SOURCES})
target_link_libraries(neu2fstr hecmw fistr ${fistr_LIBRARIES})
target_include_directories(neu2fstr PUBLIC ${neu2fstr_INCLUDE_DIRS})

###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

project(hecmw C CXX Fortran)

#
# Create the hecmw library target first
#
add_library(hecmw STATIC)

#
# Inherit definitions from FrontISTR options
#
target_include_directories(hecmw PUBLIC ${FrontISTR_INCLUDE_DIRS})
target_compile_definitions(hecmw PUBLIC ${FrontISTR_DEFINITIONS})
target_link_libraries(hecmw PRIVATE ${FrontISTR_LIBRARIES})

#
# Add subdirectories that will contribute sources to the hecmw target
#
add_subdirectory(src)

#
# Build tools
#
if(WITH_TOOLS)
  add_subdirectory(tools)
endif()

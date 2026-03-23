###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

project(fistr CXX Fortran C)

#
# Create the fistr library target first
#
add_library(fistr STATIC)

#
# Create the fistr1 executable target
#
add_executable(fistr1)

#
# Inherit definitions from FrontISTR options and link with hecmw
#
target_include_directories(fistr PUBLIC ${FrontISTR_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/hecmw1)
target_compile_definitions(fistr PUBLIC ${FrontISTR_DEFINITIONS})
target_link_libraries(fistr PUBLIC hecmw)

#
# Add subdirectories that will contribute sources to the fistr targets
#
add_subdirectory(src)

#
# Link fistr1 with fistr and configure additional properties
#
target_link_libraries(fistr1 PRIVATE 
  fistr
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)

set_target_properties(fistr1 PROPERTIES LINKER_LANGUAGE "CXX")

target_include_directories(fistr1 PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

#
# Set installing files
#
install(TARGETS fistr1
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)

#
# Build tools
#
if(WITH_TOOLS)
  add_subdirectory(tools)
endif()

###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

project(libfemmpc Fortran)
set(contact_SOURCES)
set(element_SOURCES)
set(util_SOURCES)

include(${CMAKE_CURRENT_LIST_DIR}/../src/lib/contact/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/../src/lib/element/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/../src/lib/utilities/CMakeLists.txt)

list(APPEND contact_SOURCES ${element_SOURCES})
list(APPEND contact_SOURCES ${util_SOURCES})

add_library(femmpc STATIC ${contact_SOURCES})
install(TARGETS femmpc
  DESTINATION lib
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)


###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

#
# Add subdirectories that will contribute sources to the neu2fstr target
#
add_subdirectory(HECD)
add_subdirectory(NFD)
add_subdirectory(converter)

target_sources(neu2fstr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/neu2fstr.cpp
)


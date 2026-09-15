###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

#
# Create the neu2fstr executable target first
#
add_executable(neu2fstr)

# Converter only: do not link solver libraries such as MKL or Trilinos.
#
# Add subdirectories that will contribute sources to the neu2fstr target
# Use second parameter to avoid directory-target name conflict
#
add_subdirectory(neu2fstr neu2fstr_dir)

#
# Set target properties
#
set_target_properties(neu2fstr PROPERTIES LINKER_LANGUAGE "CXX")

#
# Set installing files
#
install(TARGETS neu2fstr
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)

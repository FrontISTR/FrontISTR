###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_subdirectory(bilu)
add_subdirectory(diag)
add_subdirectory(ml)
add_subdirectory(rif)
add_subdirectory(sainv)
add_subdirectory(ssor)

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond.f90
)

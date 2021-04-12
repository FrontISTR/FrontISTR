###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_subdirectory(11)
add_subdirectory(22)
add_subdirectory(33)
add_subdirectory(44)
add_subdirectory(66)
add_subdirectory(nn)

target_sources(hecmw
	PRIVATE
    hecmw_precond.f90
)

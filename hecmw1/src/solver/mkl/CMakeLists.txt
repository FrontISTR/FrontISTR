###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
	PRIVATE
    hecmw_MKL_wrapper.F90
    hecmw_solver_direct_MKL.f90
)

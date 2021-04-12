###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
	PRIVATE
    hecmw_ordering_qmd.f90
    hecmw_ordering_metis.F90
    hecmw_ordering_rcm.f90
    hecmw_ordering.F90
    hecmw_solver_direct.f90
)

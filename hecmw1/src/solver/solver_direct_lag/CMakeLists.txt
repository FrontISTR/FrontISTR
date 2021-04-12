###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
	PRIVATE
    my_hecmw_util_f_lag.f90
    m_irjc_matrix_lag.f90
    m_child_matrix_lag.f90
    m_crs_matrix_lag.f90
    hecmw_solver_direct_serial_lag.f90
)

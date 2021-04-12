###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
	PRIVATE
    hecmw_ML_helper_nn_c.c
    hecmw_ML_wrapper.c
    hecmw_precond_BILU_nn.f90
    hecmw_precond_DIAG_nn.f90
    hecmw_precond_SSOR_nn.f90
    hecmw_precond_SAINV_nn.f90
    hecmw_precond_RIF_nn.f90
    hecmw_ML_helper.f90
    hecmw_ML_helper_nn_f.f90
    hecmw_precond_ML_nn.f90
    hecmw_precond_nn.f90
  PUBLIC
    hecmw_ML_helper.h
    hecmw_ML_helper_nn.h
)


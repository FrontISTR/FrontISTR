###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
  PRIVATE
    hecmw_precond_BILU_44.f90
    hecmw_precond_DIAG_44.f90
    hecmw_precond_SSOR_44.f90
    hecmw_precond_44.f90
)


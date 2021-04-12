###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
  PRIVATE
    hecmw_solver_CG.f90
    hecmw_solver_BiCGSTAB.f90
    hecmw_solver_GMRES.f90
    hecmw_solver_GPBiCG.f90
    hecmw_solver_Iterative.f90
)


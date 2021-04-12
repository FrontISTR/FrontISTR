###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr
	PRIVATE
  fstr_EIG_lanczos_util.f90
  fstr_EIG_setMASS.f90
  fstr_EIG_tridiag.f90
  fstr_EIG_output.f90
  fstr_EIG_lanczos.f90
  fstr_solve_eigen.f90
)

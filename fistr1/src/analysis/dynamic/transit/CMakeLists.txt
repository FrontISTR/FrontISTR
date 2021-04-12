###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr
	PRIVATE
  table_dyn.f90
  dynamic_mat_ass_load.f90
  dynamic_mat_ass_bc.f90
  dynamic_mat_ass_bc_vl.f90
  dynamic_mat_ass_bc_ac.f90
  dynamic_mat_ass_couple.f90
  dynamic_var_init.f90
  dynamic_output.f90
  fstr_dynamic_nlexplicit.f90
  fstr_dynamic_nlimplicit.f90
  fstr_solve_dynamic.f90
)

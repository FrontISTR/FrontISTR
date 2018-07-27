!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw
  use hecmw_allocate
  use hecmw_control
  use hecmw_dist_copy_c2f_f
  use hecmw_dist_copy_f2c_f
  use hecmw_dist_free_f
  use hecmw_dist_print_f
  use hecmw_etype
  use hecmw_comm_group
  use hecmw_io
  use hecmw_local_matrix
  use hecmw_logging
  use hecmw_matrix_ass
  use hecmw_matrix_con
  use hecmw_matrix_contact
  use hecmw_matrix_misc
  use hecmw_matrix_dump
  use hecmw_mpc_prepost
  use hecmw_msg
  use hecmw_msgno
  use hecmw_restart
  use hecmw_result
  !use hecmw_solve_sai_make_33
  use hecmw_solver
  use hecmw_solver_las
  !use hecmw_solver_sai_bicgstab_33
  !use hecmw_solver_sai_gmres_33
  !use hecmw_solver_sai_gpbicg_33
  use hecmw_solver_misc
  use hecmw_solver_sr_11
  use hecmw_solver_sr_11i
  use hecmw_solver_sr_22
  use hecmw_solver_sr_22i
  use hecmw_solver_sr_33
  use hecmw_solver_sr_33i
  use hecmw_solver_sr_mm
  use hecmw_solver_sr_mmi
  use hecmw_util
  use hecmw_visualizer
  use m_hecmw_comm_f
  use m_hecmw_solve_error
  use m_hecmw_solve_init
end module hecmw

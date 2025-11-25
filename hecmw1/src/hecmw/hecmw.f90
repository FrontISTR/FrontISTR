!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw
  use hecmw_allocate
  use hecmw_amplitude
  use hecmw_control
  use hecmw_dist_copy_c2f_f
  use hecmw_dist_copy_f2c_f
  use hecmw_dist_free_f
  use hecmw_dist_print_f
  use hecmw_etype
  use hecmw_comm_group
  use hecmw_varray_int
  use hecmw_io
  use hecmw_local_matrix
  use hecmw_logging
  use hecmw_matrix_ass
  use hecmw_matrix_con
  use hecmw_matrix_contact_lagrange
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
  use hecmw_solver_misc
  use hecmw_solver_sr
  use hecmw_solver_sr_i
  use hecmw_util
  use hecmw_visualizer
  use m_hecmw_comm_f
  use m_hecmw_solve_error
  use m_hecmw_solve_init
  use hecmw_es_mesh_connectivity
end module hecmw

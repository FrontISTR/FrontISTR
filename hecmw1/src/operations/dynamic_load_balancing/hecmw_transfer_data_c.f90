!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Dynamic Load Balancing

subroutine hecmw_transfer_data_f2c (hecMESHnew, adapRES)

  use  hecmw_util
  use  hecmw_io
  use  hecmw_dist_copy_f2c_f
  use  hecmw_result
  type (hecmwST_local_mesh) :: hecMESHnew
  type (hecmwST_result_data):: adapRES

  call hecmw_dlb_f2c_init
  call hecmw_dist_copy_f2c(hecMESHnew, ierr)
  call hecmw_dlb_f2c_finalize
  call test_mesh
  call hecmw_dist_result_copy_f2c(adapRES)
end subroutine hecmw_transfer_data_f2c


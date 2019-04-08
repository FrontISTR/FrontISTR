!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Dynamic Load Balancing

subroutine hecmw_transfer_data_c2f (hecMESHnew, adapRES)

  use  hecmw_util
  use  hecmw_io
  use  hecmw_dist_copy_c2f_f
  use  hecmw_result
  type (hecmwST_local_mesh) :: hecMESHnew
  type (hecmwST_result_data):: adapRES
  write (*,'(/,a)') '#start transfer from c to fortran'

  call hecmw_dist_free(hecMESHnew)
  write (*,'(/,a)') '#ok to deallocate'
  call hecmw_dlb_c2f_init()
  write(*,*) 'my_rank=',hecMESHnew%my_rank
  if(hecMESHnew%my_rank .eq. 0) write (*,'(/,a)') '#ok to initialize'
  call hecmw_dist_copy_c2f(hecMESHnew, ierr)
  write (*,'(/,a)') '#ok to copy'
  print *, hecMESHnew%n_node
  call hecmw_dlb_c2f_finalize
  call hecmw_dist_result_copy_c2f(hecMESHnew, adapRES)
end subroutine hecmw_transfer_data_c2f


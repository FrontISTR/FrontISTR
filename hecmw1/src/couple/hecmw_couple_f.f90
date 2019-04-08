!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

module hecmw_couple_f

  use hecmw_util
  use hecmw_couple_struct_f
  use hecmw_couple_copy_c2f_f
  use hecmw_couple_copy_f2c_f

  implicit none
  private
  public :: hecmw_couple

contains

subroutine hecmw_couple(boundary_id, couple_value)

  character(len=HECMW_NAME_LEN), intent(in)    :: boundary_id
  type(hecmw_couple_value),      intent(inout) :: couple_value
  integer(kind=kint)                           :: ierr

  call hecmw_couple_exec_init_if(ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_copy_f2c(couple_value, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_if(boundary_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_copy_c2f(couple_value, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_exec_finalize_if(ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

end subroutine hecmw_couple

end module hecmw_couple_f

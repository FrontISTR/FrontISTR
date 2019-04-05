!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

module hecmw_couple_finalize_f

  use hecmw_util

  implicit none
  private
  public :: hecmw_couple_finalize

contains

subroutine hecmw_couple_finalize(boundary_id)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: ierr

  call hecmw_couple_finalize_if(boundary_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

end subroutine hecmw_couple_finalize

end module hecmw_couple_finalize_f

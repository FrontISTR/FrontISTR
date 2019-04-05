!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!> \brief I/O and Utility

module hecmw_control
  use hecmw_util
  implicit none

contains

  subroutine hecmw_ctrl_get_control_file(name_ID, filename)
    character(len=HECMW_NAME_LEN) :: name_ID
    character(len=HECMW_FILENAME_LEN) :: filename
    integer(kind=kint) :: ierr

    call hecmw_ctrl_get_control_file_if(name_ID, filename, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_ctrl_get_control_file

end module hecmw_control


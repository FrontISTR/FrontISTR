!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module hecmw_msg
  !      use hecmw_util
  use hecmw_msgno

contains

  subroutine hecmw_strmsg(msgno, msg)
    integer(kind=kint) :: msgno
    character(len=HECMW_MSG_LEN) :: msg

    call hecmw_strmsg_if(msgno, msg)
  end subroutine hecmw_strmsg
end module hecmw_msg

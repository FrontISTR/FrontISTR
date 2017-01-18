!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to deal with time and increment of stress analysis

module m_fstr_TimeInc
  use m_fstr

  implicit none

  real(kind=kreal), private :: current_time  = 0.d0  !< current time
  real(kind=kreal), private :: time_inc = 0.d0  !< current increment of time

  contains

  !! basic functions to get time/timeinc
  real(kind=kreal) function fstr_get_time()
    fstr_get_time = current_time
  end function

  real(kind=kreal) function fstr_get_timeinc()
    fstr_get_timeinc = time_inc
  end function

  real(kind=kreal) function fstr_get_timeinc_base()
    fstr_get_timeinc_base = time_inc
  end function

  !! basic functions to set/update time/timeinc
  subroutine fstr_set_time( time )
    real(kind=kreal), intent(in) :: time  !< input time
    current_time = time
  end subroutine

  subroutine fstr_set_timeinc( dtime )
    real(kind=kreal), intent(in) :: dtime  !< input increment of time
    time_inc = dtime
  end subroutine

  subroutine fstr_set_timeinc_base( dtime_base )
    real(kind=kreal), intent(in) :: dtime_base  !< input increment of time
    time_inc = dtime_base
  end subroutine

  subroutine fstr_proceed_time()
    call fstr_set_time( current_time+time_inc )
  end subroutine

end module m_fstr_TimeInc

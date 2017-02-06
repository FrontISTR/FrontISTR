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

  !! return if step finished
  logical function fstr_TimeInc_isStepFinished( stepinfo )
    type(step_info), intent(in)  :: stepinfo

    real(kind=kreal) :: endtime, remain_time

    endtime = stepinfo%starttime + stepinfo%elapsetime
    remain_time = (endtime-current_time)/stepinfo%elapsetime
    fstr_TimeInc_isStepFinished = (remain_time < 1.d-12)
  end function

  !! set time_inc from stepinfo and Newton Raphson iteration status
  subroutine fstr_TimeInc_SetTimeIncrement( stepinfo, fstrPARAM, substep, NRstatI, NRstatR, AutoINC_stat, Cutback_stat )
    type(step_info), intent(in)       ::  stepinfo
    type(fstr_param), intent(in)      ::  fstrPARAM
    integer(kind=kint), intent(in)    ::  substep
    integer(kind=kint), intent(in)    ::  NRstatI(:)  !previous Newton Raphson Iteration Log
    real(kind=kreal), intent(in)      ::  NRstatR(:)  !previous Newton Raphson Iteration Log
    integer(kind=kint), intent(inout) ::  AutoINC_stat
    integer(kind=kint), intent(in)    ::  Cutback_stat

    real(kind=kreal) :: timeinc0
    real(kind=kreal) :: endtime
    logical          :: to_be_decreased, to_be_increased

    if( substep == 1 .or. stepinfo%inc_type == stepFixedInc ) then
      timeinc0 = stepinfo%initdt
    else ! INCTYPE==AUTO and substep > 2
      timeinc0 = time_inc

      if( Cutback_stat > 0 ) then
        timeinc0 = fstrPARAM%ainc_Rc*time_inc
        write(*,'(2(A,E10.3))') 'time increment is decreased from ', time_inc, ' to ', timeinc0
        AutoINC_stat = -1
      else
        !decrease condition
        to_be_decreased = .false.
        if( NRstatI(knstMAXIT) > fstrPARAM%NRbound_s(knstMAXIT) ) to_be_decreased = .true.
        if( NRstatI(knstSUMIT) > fstrPARAM%NRbound_s(knstSUMIT) ) to_be_decreased = .true.
        if( NRstatI(knstCITER) > fstrPARAM%NRbound_s(knstCITER) ) to_be_decreased = .true.

        !increase condition
        to_be_increased = .true.
        if( NRstatI(knstMAXIT) > fstrPARAM%NRbound_l(knstMAXIT) ) to_be_increased = .false.
        if( NRstatI(knstSUMIT) > fstrPARAM%NRbound_l(knstSUMIT) ) to_be_increased = .false.
        if( NRstatI(knstCITER) > fstrPARAM%NRbound_l(knstCITER) ) to_be_increased = .false.

        ! count # of times that increase/decrease condition has been satisfied
        ! AutoINC_stat < 0 ... decrease condition has been satisfied -AutoINC_stat times
        ! AutoINC_stat > 0 ... increase condition has been satisfied AutoINC_stat times
        if(to_be_decreased) then
          AutoINC_stat = min(AutoINC_stat,0) - 1
        else if(to_be_increased) then
          AutoINC_stat = max(AutoINC_stat,0) + 1
        else
          AutoINC_stat = 0
        end if

        if( AutoINC_stat <= -fstrPARAM%NRtimes_s ) then
          timeinc0 = fstrPARAM%ainc_Rs*time_inc
          write(*,'(2(A,E10.3))') 'time increment is decreased from ', time_inc, ' to ', timeinc0
        else if( AutoINC_stat >= fstrPARAM%NRtimes_l ) then
          timeinc0 = dmin1(fstrPARAM%ainc_Rl*time_inc,stepinfo%maxdt)
          write(*,'(2(A,E10.3))') 'time increment is increased from ', time_inc, ' to ', timeinc0
        end if
      end if
    end if

    if( timeinc0 < stepinfo%mindt ) then
      write(*,'(2(A,E10.3))') 'Error: current time increment ',timeinc0,' is smaller than lower bound',stepinfo%mindt
      stop
    end if

    endtime = stepinfo%starttime + stepinfo%elapsetime
    time_inc = dmin1(timeinc0,endtime-current_time)

  end subroutine


end module m_fstr_TimeInc

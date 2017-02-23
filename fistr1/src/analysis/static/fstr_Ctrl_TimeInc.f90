!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to deal with time and increment of stress analysis

module m_fstr_TimeInc
  use m_fstr
  use m_timepoint

  implicit none

  real(kind=kreal), private :: current_time  = 0.d0  !< current time
  real(kind=kreal), private :: time_inc      = 0.d0  !< current increment of time
  real(kind=kreal), private :: time_inc_base = 0.d0  !< base increment of time

  contains

  !! basic functions to get time/timeinc
  real(kind=kreal) function fstr_get_time()
    fstr_get_time = current_time
  end function

  real(kind=kreal) function fstr_get_timeinc()
    fstr_get_timeinc = time_inc
  end function

  real(kind=kreal) function fstr_get_timeinc_base()
    fstr_get_timeinc_base = time_inc_base
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
    time_inc_base = dtime_base
  end subroutine

  subroutine fstr_proceed_time()
    current_time = current_time + time_inc
  end subroutine

  !! functions to print status
  subroutine fstr_TimeInc_PrintSTATUS_init
    write(ISTA,'(A10,"-+-",A60,"-+-",A40)') REPEAT("-",10),REPEAT("-",60),REPEAT("-",40)
    write(ISTA,'(2A5," | ",2A5,2A7,3A12," | ",A40)') '     ', '     ', ' ', ' # of',  'MAX #', 'TOT #',  '','', '', ''
    write(ISTA,'(2A5," | ",2A5,2A7,3A12," | ",A7,A33)') &
      & 'STEP', 'SUB', 'STAT', ' CONT', 'NEWTON', 'NEWTON', 'START', 'TIME','END', 'MESSAGE',''
    write(ISTA,'(2A5," | ",2A5,2A7,3A12," | ",A40)') &
      & '', 'STEP', '  ', 'ITER',  'ITER', 'ITER',  'TIME', 'INC', 'TIME',''
    write(ISTA,'(A10,"-+-",A60,"-+-",A40)') REPEAT("-",10),REPEAT("-",60),REPEAT("-",40)
  end subroutine

  subroutine fstr_TimeInc_PrintSTATUS( stepinfo, fstrPARAM, totstep, substep, NRstatI, NRstatR, AutoINC_stat, Cutback_stat )
    type(step_info), intent(in)       ::  stepinfo
    type(fstr_param), intent(in)      ::  fstrPARAM
    integer(kind=kint), intent(in)    ::  totstep
    integer(kind=kint), intent(in)    ::  substep
    integer(kind=kint), intent(in)    ::  NRstatI(:)  !previous Newton Raphson Iteration Log
    real(kind=kreal), intent(in)      ::  NRstatR(:)  !previous Newton Raphson Iteration Log
    integer(kind=kint), intent(inout) ::  AutoINC_stat
    integer(kind=kint), intent(in)    ::  Cutback_stat

    character(len=5) :: cstep, csstep, cstate, ccont
    character(len=7) :: cmaxn, ctotn
    character(len=12) :: ctime, cdtime, etime
    character(len=40) :: message

    write(cstep,'(I5)') totstep
    write(csstep,'(I5)') substep
    write(ccont,'(I5)') NRstatI(knstCITER)
    write(cmaxn,'(I7)') NRstatI(knstMAXIT)
    write(ctotn,'(I7)') NRstatI(knstSUMIT)
    write(ctime,'(1pE12.4)') current_time
    write(cdtime,'(1pE12.4)') time_inc

    if( Cutback_stat > 0 ) then
      write(etime,'(1pE12.4)') current_time
      write(cstate,'(I4,A)') Cutback_stat,'F'
      if( NRstatI(knstDRESN) == 1 ) write(message,'(A)') 'Failed to converge due to MAXITER.'
      if( NRstatI(knstDRESN) == 2 ) write(message,'(A)') 'Failed to converge due to MAXRES.'
      if( NRstatI(knstDRESN) == 3 ) write(message,'(A)') 'Failed to converge due to MAXCONTITER.'
      if( Cutback_stat == fstrPARAM%CBbound ) write(message,'(A)') '# of successive cutback reached max.'
    else
      write(etime,'(1pE12.4)') current_time+time_inc
      write(cstate,'(A5)') 'S'
      write(message,'(A)') ''
    endif

    write(ISTA,'(2A5," | ",2A5,2A7,3A12," | ",A40)') cstep, csstep, cstate, ccont, &
      &  cmaxn, ctotn, ctime, cdtime, etime, message

  end subroutine

  subroutine fstr_TimeInc_PrintSTATUS_final(success_flag)
    logical, intent(in) :: success_flag
    write(ISTA,'(A10,"-+-",A60,"-+-",A40)') REPEAT("-",10),REPEAT("-",60),REPEAT("-",40)
    if(success_flag) then
      write(ISTA,'(A)') 'FSTR_SOLVE_NLGEOM HAS COMPLETED SUCCESSFULLY'
    else
      write(ISTA,'(A)') 'FSTR_SOLVE_NLGEOM HAS NOT COMPLETED SUCCESSFULLY'
    end if
  end subroutine

  !! true if step finished
  logical function fstr_TimeInc_isStepFinished( stepinfo )
    type(step_info), intent(in)  :: stepinfo

    real(kind=kreal) :: endtime, remain_time

    endtime = stepinfo%starttime + stepinfo%elapsetime
    remain_time = (endtime-current_time)/stepinfo%elapsetime
    fstr_TimeInc_isStepFinished = (remain_time < 1.d-12)
  end function

  !! true if current time is time point
  logical function fstr_TimeInc_isTimePoint( stepinfo, fstrPARAM )
    type(step_info), intent(in)    :: stepinfo
    type(fstr_param), intent(in)   :: fstrPARAM

    fstr_TimeInc_isTimePoint = .false.
    if( stepinfo%inc_type == stepFixedInc ) return
    if( stepinfo%timepoint_id == 0 ) return

    fstr_TimeInc_isTimePoint = is_at_timepoints(current_time,stepinfo%starttime, &
      &  fstrPARAM%timepoints(stepinfo%timepoint_id))
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
    real(kind=kreal) :: endtime, remain
    logical          :: to_be_decreased, to_be_increased

    if( stepinfo%inc_type == stepFixedInc ) then
      timeinc0 = stepinfo%initdt
    else ! INCTYPE==AUTO
      if( Cutback_stat > 0 ) then
        timeinc0 = fstrPARAM%ainc_Rc*time_inc_base
        if(myrank == 0) write(*,'(2(A,E10.3))') 'time increment is decreased from ', time_inc_base, ' to ', timeinc0
        AutoINC_stat = -1
      else
        if( substep == 1 ) then
          timeinc0 = stepinfo%initdt
          AutoINC_stat = 0
        else
          timeinc0 = time_inc_base

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
            timeinc0 = fstrPARAM%ainc_Rs*timeinc0
            if(myrank == 0) write(*,'(2(A,E10.3))') 'time increment is decreased from ', time_inc_base, ' to ', timeinc0
          else if( AutoINC_stat >= fstrPARAM%NRtimes_l ) then
            timeinc0 = dmin1(fstrPARAM%ainc_Rl*timeinc0,stepinfo%maxdt)
            if(myrank == 0) write(*,'(2(A,E10.3))') 'time increment is increased from ', time_inc_base, ' to ', timeinc0
          end if
        end if
      end if
    end if

    if( timeinc0 < stepinfo%mindt ) then
      if(myrank == 0) then
        write(*,'(2(A,E10.3))') 'Error: current time increment ',timeinc0,' is smaller than lower bound',stepinfo%mindt
        call fstr_TimeInc_PrintSTATUS_final(.false.)
      endif
      stop
    end if

    endtime = stepinfo%starttime + stepinfo%elapsetime
    time_inc_base = dmin1(timeinc0,endtime-current_time)
    if( stepinfo%timepoint_id > 0 ) then
      remain = get_remain_to_next_timepoints(current_time,stepinfo%starttime,fstrPARAM%timepoints(stepinfo%timepoint_id))
      time_inc = dmin1(time_inc_base,remain)
    else
      time_inc = time_inc_base
    end if

  end subroutine

end module m_fstr_TimeInc

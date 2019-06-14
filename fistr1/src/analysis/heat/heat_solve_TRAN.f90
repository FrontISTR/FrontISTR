!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to control nonsteady heat analysis
module m_heat_solve_TRAN
contains

  subroutine heat_solve_TRAN ( hecMESH,hecMAT,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP,total_step )
    use m_fstr
    use m_heat_mat_ass_conductivity
    use m_heat_mat_ass_capacity
    use m_heat_mat_ass_boundary
    use m_heat_init
    use m_heat_solve_main
    use m_solve_lineq
    use m_heat_io
    implicit none
    integer(kind=kint) :: ISTEP, ITM, incr, iterALL, i, inod, mnod, max_step, tstep, interval, bup_n_dof
    integer(kind=kint) :: total_step
    real(kind=kreal)   :: start_time, delta_time, current_time, next_time, total_time, end_time, DELMAX, DELMIN
    real(kind=kreal)   :: val, CHK, tmpmax, dltmp, tmpmax_myrank
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_heat)           :: fstrHEAT
    type(hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint), parameter :: miniter = 4
    logical :: is_end, outflag

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMATmpc)

    if(ISTEP == 1)then
      start_time = 0.0d0
    else
      start_time = 0.0d0
      do i = 1, ISTEP - 1
        start_time = start_time + fstrHEAT%STEP_EETIME(i)
      enddo
    endif

    delta_time = fstrHEAT%STEP_DLTIME(ISTEP)
    end_time = fstrHEAT%STEP_EETIME(ISTEP)
    DELMIN = fstrHEAT%STEP_DELMIN(ISTEP)
    DELMAX = fstrHEAT%STEP_DELMAX(ISTEP)

    is_end = .false.
    hecMAT%NDOF = 1
    hecMAT%Iarray(98) = 1 !Assmebly complete
    hecMAT%X = 0.0d0
    current_time = 0.0d0
    next_time = 0.0d0

    if(fstrHEAT%is_steady == 1)then
      fstrHEAT%beta = 1.0d0
    else
      fstrHEAT%beta = 0.5d0
    endif

    call heat_input_restart(fstrHEAT, hecMESH, total_step, start_time)

    !C--------------------   START TRANSIET LOOP   ------------------------
    tr_loop: do

      if(end_time <= current_time + delta_time + delta_time*1.0d-6) then
        delta_time = end_time - current_time
        next_time = end_time
        is_end = .true.
      else
        next_time = current_time + delta_time
        is_end = .false.
      endif
      if( fstrHEAT%is_steady == 1 ) is_end = .true.
      total_time = start_time + next_time

      if( hecMESH%my_rank.eq.0 ) then
        write(IMSG,"(a,i8,a,1pe12.5,a,1pe12.5)") " ** Increment No. :", total_step, ", current time: ", &
        & current_time, ", delta t: ", delta_time
        write(*,   "(a,i8,a,1pe12.5,a,1pe12.5)") " ** Increment No. :", total_step, ", current time: ", &
        & current_time, ", delta t: ", delta_time
      endif

      if(delta_time < DELMIN)then
        if(hecMESH%my_rank == 0) write(IMSG,*) ' !!! DELTA TIME EXCEEDED TOLERANCE OF TIME INCREMENT'
        call hecmw_abort(hecmw_comm_get_comm())
      endif

      call heat_solve_main(hecMESH, hecMAT, hecMATmpc, fstrPARAM, fstrHEAT, ISTEP, iterALL, total_time, delta_time)

      if(0.0d0 < DELMIN)then
        tmpmax = 0.0d0
        do i = 1, hecMESH%nn_internal
          inod = fstrPARAM%global_local_id(1,i)
          dltmp = fstrHEAT%TEMP0(i) - fstrHEAT%TEMP(i)
          if(tmpmax < dabs(dltmp)) then
            mnod = inod
            tmpmax = dabs(dltmp)
          endif
        enddo
        tmpmax_myrank = tmpmax
        call hecmw_allREDUCE_R1(hecMESH, tmpmax, hecmw_max)

        if(tmpmax_myrank < tmpmax) mnod = -1
        call hecmw_allREDUCE_I1(hecMESH, mnod, hecmw_max)

        if(DELMAX < tmpmax .or. fstrHEAT%is_iter_max_limit)then
          if(hecMESH%my_rank == 0)then
            write(*,*) ' *** EXCEEDED TOLERANCE OF VARIATION IN TEMPERATUTE.'
            write(*,*) ' : NODE NUMBER  = ', mnod, ' : DELTA TEMP   = ', tmpmax
          endif
          delta_time = 0.5d0*delta_time
          cycle tr_loop
        endif

        if(iterALL <= miniter) delta_time = delta_time*1.5d0
      else
        if(fstrHEAT%is_iter_max_limit) call hecmw_abort( hecmw_comm_get_comm() )
      endif

      do i = 1, hecMESH%n_node
        fstrHEAT%TEMP0(i) = fstrHEAT%TEMP(i)
      enddo

      call heat_output_log(hecMESH, fstrPARAM, fstrHEAT, total_step, next_time)
      call heat_output_result(hecMESH, fstrHEAT, total_step, next_time, is_end)
      call heat_output_visual(hecMESH, fstrRESULT, fstrHEAT, total_step, next_time, is_end)
      call heat_output_restart(hecMESH, fstrHEAT, total_step, is_end, next_time)

      total_step = total_step + 1
      current_time = next_time
      if( is_end ) exit
    enddo tr_loop
    !C--------------------   END TRANSIET LOOP   ------------------------

    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMATmpc)

  end subroutine heat_solve_TRAN
end module m_heat_solve_TRAN

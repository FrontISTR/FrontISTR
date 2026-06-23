!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provide a function to control nonsteady heat analysis
module m_heat_solve_TRAN
contains

  subroutine heat_solve_TRAN ( hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP,total_step,current_time )
    use m_fstr
    use m_heat_mat_ass_conductivity
    use m_heat_mat_ass_capacity
    use m_heat_mat_ass_boundary
    use m_heat_init
    use m_heat_solve_main
    use m_solve_lineq
    use m_heat_io
    implicit none
    integer(kind=kint) :: ISTEP, iterALL, i, inod, mnod
    integer(kind=kint) :: total_step
    real(kind=kreal)   :: start_time, delta_time_base, delta_time, current_time, next_time, total_time, end_time, DELMAX, DELMIN
    real(kind=kreal)   :: tmpmax, dltmp, tmpmax_myrank, remain_time
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(fstr_solid)          :: fstrSOLID
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_heat)           :: fstrHEAT
    type(hecmwST_local_mesh), pointer :: hecMESHmpc
    type(hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint), parameter :: miniter = 4
    logical :: is_end, outflag

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    start_time = 0.0d0
    do i = 1, ISTEP - 1
      start_time = start_time + fstrHEAT%STEP_EETIME(i)
    enddo

    total_time = start_time + current_time

    delta_time_base = fstrHEAT%STEP_DLTIME(ISTEP)
    end_time = fstrHEAT%STEP_EETIME(ISTEP)
    DELMIN = fstrHEAT%STEP_DELMIN(ISTEP)
    DELMAX = fstrHEAT%STEP_DELMAX(ISTEP)

    is_end = .false.
    hecMAT%NDOF = 1
    hecMAT%Iarray(98) = 1 !Assembly complete
    hecMAT%X = 0.0d0

    if(fstrHEAT%beta == -1.0d0)then
      if(fstrHEAT%is_steady == 1)then
        fstrHEAT%beta = 1.0d0
      else
        fstrHEAT%beta = 0.5d0
      endif
    endif

    if(fstrHEAT%is_steady /= 1 .and. total_step == 1) then
      call heat_output_result(hecMESH, fstrHEAT, fstrSOLID, 0, total_time, .true.)
      call heat_output_visual(hecMESH, fstrRESULT, fstrHEAT, fstrSOLID, 0, total_time, .true.)
    endif

    !C--------------------   START TRANSIENT LOOP   ------------------------
    tr_loop: do

      if(end_time <= current_time + delta_time_base + delta_time_base*1.0d-6) then
        delta_time_base = end_time - current_time
      endif
      if( 0.0d0 < DELMIN .and. fstrHEAT%timepoint_id > 0 ) then
        remain_time = get_remain_to_next_timepoints(total_time, 0.0d0, fstrPARAM%timepoints(fstrHEAT%timepoint_id))
        delta_time = dmin1(delta_time_base, remain_time)
      else
        delta_time = delta_time_base
      endif
      next_time = current_time + delta_time
      total_time = start_time + next_time

      if( fstrHEAT%is_steady == 1 ) then
        is_end = .true.
      else
        if( (end_time - next_time) / end_time < 1.d-12 ) is_end = .true.
      endif

      if( 0.0d0 < DELMIN .and. fstrHEAT%timepoint_id > 0 ) then
        outflag = is_end .or. is_at_timepoints(total_time, 0.0d0, fstrPARAM%timepoints(fstrHEAT%timepoint_id))
      else
        outflag = is_end
      endif

      if( hecMESH%my_rank.eq.0 ) then
        write(IMSG,"(a,i8,a,1pe12.5,a,1pe12.5)") " ** Increment No. :", total_step, ", total time: ", &
        & total_time, ", delta t: ", delta_time
        write(*,   "(a,i8,a,1pe12.5,a,1pe12.5)") " ** Increment No. :", total_step, ", total time: ", &
        & total_time, ", delta t: ", delta_time
      endif

      if(delta_time_base < DELMIN .and. (.not. is_end))then
        if(hecMESH%my_rank == 0) write(IMSG,*) ' !!! DELTA TIME EXCEEDED TOLERANCE OF TIME INCREMENT'
        call hecmw_abort(hecmw_comm_get_comm())
      endif

      call heat_solve_main(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, &
       & fstrSOLID, fstrPARAM, fstrHEAT, ISTEP, iterALL, total_time, delta_time)

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
          delta_time_base = 0.5d0*delta_time_base
          cycle tr_loop
        endif

        if(iterALL <= miniter) delta_time_base = delta_time_base*1.5d0
      else
        if(fstrHEAT%is_iter_max_limit) call hecmw_abort( hecmw_comm_get_comm() )
      endif

      do i = 1, hecMESH%n_node
        fstrHEAT%TEMP0(i) = fstrHEAT%TEMP(i)
      enddo

      call heat_output_log(hecMESH, fstrPARAM, fstrHEAT, total_step, total_time)
      call heat_output_result(hecMESH, fstrHEAT, fstrSOLID, total_step, total_time, outflag)
      call heat_output_visual(hecMESH, fstrRESULT, fstrHEAT, fstrSOLID, total_step, total_time, outflag)
      call heat_output_restart(hecMESH, fstrHEAT, ISTEP, total_step, next_time, outflag)

      total_step = total_step + 1
      current_time = next_time
      if( is_end ) exit
    enddo tr_loop
    !C--------------------   END TRANSIENT LOOP   ------------------------

    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

  end subroutine heat_solve_TRAN
end module m_heat_solve_TRAN

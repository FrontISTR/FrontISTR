!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides main suboruitne for nonliear calculation.

module m_fstr_solve_NLGEOM
  use m_fstr
  use m_static_lib
  use m_static_output
  use m_fstr_NonLinearMethod
  use m_fstr_Restart
  use fstr_matrix_con_contact
  use m_fstr_TimeInc
  use m_fstr_Cutback
  use mContact
  use m_solve_LINEQ_contact

  implicit none

contains

  !======================================================================!
  !> \brief This module provides main subroutine for nonlinear calculation.
  !>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
  !>              Z. Sun(ASTOM)/2010/11
  !>  \date       2009/08/31
  !>  \version    0.00
  subroutine FSTR_SOLVE_NLGEOM(hecMESH,hecMAT,fstrSOLID,hecLagMAT,fstrPARAM,conMAT)
    type (hecmwST_local_mesh)              :: hecMESH      !< mesh information
    type (hecmwST_matrix    )              :: hecMAT       !< linear equation, its right side modified here
    type (fstr_param       )               :: fstrPARAM    !< analysis control parameters
    type (fstr_solid       )               :: fstrSOLID    !< we need boundary conditions of curr step
    type (hecmwST_matrix_lagrange)         :: hecLagMAT    !< type hecmwST_matrix_lagrange
    type (fstr_info_contactChange)         :: infoCTChange, infoCTChange_bak !< type fstr_info_contactChange
    type (hecmwST_matrix    )              :: conMAT

    integer(kind=kint) :: ndof, nn
    integer(kind=kint) :: j, i, tot_step, step_count, tot_step_print, CBbound
    integer(kind=kint) :: sub_step
    integer(kind=kint) :: restart_step_num, restart_substep_num
    real(kind=kreal)   :: ctime, dtime, endtime, factor
    real(kind=kreal)   :: time_1, time_2
    logical            :: ctchanged, is_OutPoint, is_interaction_active

    if(hecMESH%my_rank==0) call fstr_TimeInc_PrintSTATUS_init

    hecMAT%NDOF = hecMESH%n_dof

    ndof = hecMAT%NDOF
    nn = ndof*ndof

    is_interaction_active = ( associated( fstrSOLID%contacts ) .or. associated( fstrSOLID%embeds ) )

    if( fstrSOLID%TEMP_ngrp_tot>0 .and. hecMESH%hecmw_flag_initcon==1 ) then
      fstrSOLID%last_temp = 0.0d0
      fstrSOLID%temperature = 0.0d0
      do j=1, size(hecMESH%node_init_val_item)
        i = hecMESH%node_init_val_index(j)
        fstrSOLID%last_temp(j) = hecMESH%node_init_val_item(i)
        fstrSOLID%temperature(j) = hecMESH%node_init_val_item(i)
      end do
    endif
    if( fstrSOLID%TEMP_ngrp_tot>0 .and. associated(g_InitialCnd) ) then
      fstrSOLID%last_temp = 0.0d0
      fstrSOLID%temperature = 0.0d0
	  do j=1,size(g_InitialCnd)
          if( g_InitialCnd(j)%cond_name=="temperature" ) then
            if( .not. associated(fstrSOLID%temperature) ) then
                allocate( fstrSOLID%temperature( hecMESH%n_node ) )
                allocate( fstrSOLID%temp_bak( hecMESH%n_node ) )
                allocate( fstrSOLID%last_temp( hecMESH%n_node ) )
            endif
            do i= 1, hecMESH%n_node
              fstrSOLID%last_temp(i) = g_InitialCnd(j)%realval(i)
              fstrSOLID%temperature(i) = fstrSOLID%last_temp(i)
            enddo
          endif
      end do
    endif

    if( associated( fstrSOLID%contacts ) ) then
      call initialize_contact_output_vectors(fstrSOLID,hecMAT)
      call setup_contact_elesurf_for_area( 1, hecMESH, fstrSOLID )
    endif
    if( fstrSOLID%n_embeds > 0 ) call initialize_embed_vectors(fstrSOLID,hecMAT)

    restart_step_num    = 1
    restart_substep_num = 1
    fstrSOLID%unode = 0.0d0
    step_count = 0 !**
    infoCTChange%contactNode_previous = 0
    infoCTChange%contactNode_current  = 0
    if( fstrSOLID%restart_nout < 0 ) then
      call fstr_read_restart(restart_step_num,restart_substep_num,step_count,ctime,dtime,hecMESH,fstrSOLID, &
        fstrPARAM,infoCTChange%contactNode_previous)
      hecMAT%Iarray(98) = 1
      call fstr_set_time( ctime )
      call fstr_set_timeinc_base( dtime )
      fstrSOLID%restart_nout = - fstrSOLID%restart_nout
    else
      call fstr_static_Output( 1, 0, 0.d0, hecMESH, fstrSOLID, fstrPARAM, fstrPR%solution_type, .true. )
    endif

    fstrSOLID%FACTOR = 0.0d0
    call fstr_cutback_init( hecMESH, fstrSOLID, fstrPARAM )
    call fstr_cutback_save( fstrSOLID, infoCTChange, infoCTChange_bak )

    do tot_step=1, fstrSOLID%nstep_tot
      tot_step_print = tot_step+restart_step_num-1
      if(hecMESH%my_rank==0) write(*,*) ''
      if(hecMESH%my_rank==0) write(*,'(a,i5)') ' loading step=',tot_step_print

      if( fstrSOLID%TEMP_ngrp_tot>0 ) then
        do j=1, hecMESH%n_node
          fstrSOLID%temp_bak(j) = fstrSOLID%temperature(j)
        end do
      endif
      call fstr_UpdateState( hecMESH, fstrSOLID, 0.0d0 )

      fstrSOLID%unode_bak(:) = fstrSOLID%unode(:)

      ! -------------------------------------------------------------------------
      !      STEP LOOP
      ! -------------------------------------------------------------------------
      sub_step = restart_substep_num
      do while(.true.)

        ! ----- time history of factor
        call fstr_TimeInc_SetTimeIncrement( fstrSOLID%step_ctrl(tot_step), fstrPARAM, sub_step, &
          &  fstrSOLID%NRstat_i, fstrSOLID%NRstat_r, fstrSOLID%AutoINC_stat, fstrSOLID%CutBack_stat )
        if( fstrSOLID%TEMP_irres > 0 ) then
          fstrSOLID%FACTOR(1) = 0.d0
          fstrSOLID%FACTOR(2) = 1.d0
          call table_nlsta(hecMESH,fstrSOLID,tot_step,fstr_get_time()+fstr_get_timeinc(), factor)
          fstrSOLID%TEMP_FACTOR = factor
        else
          call table_nlsta(hecMESH,fstrSOLID,tot_step,fstr_get_time(),factor)
          fstrSOLID%FACTOR(1) = factor
          call table_nlsta(hecMESH,fstrSOLID,tot_step,fstr_get_time()+fstr_get_timeinc(), factor)
          fstrSOLID%FACTOR(2) = factor
        endif

        if(hecMESH%my_rank==0) then
          write(*,'(A,I0,2(A,E12.4))') ' sub_step= ',sub_step,', &
            &  current_time=',fstr_get_time(), ', time_inc=',fstr_get_timeinc()
          write(*,'(A,2f12.7)') ' loading_factor= ', fstrSOLID%FACTOR
          if( fstrSOLID%TEMP_irres > 0 ) write(*,'(A,2f12.7)') ' readtemp_factor= ', fstrSOLID%TEMP_FACTOR
        endif

        time_1 = hecmw_Wtime()

        ! analysis algorithm ( Newton-Rapshon Method )
        if( .not. is_interaction_active ) then
          call fstr_Newton( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM,   &
            restart_step_num, sub_step, fstr_get_time(), fstr_get_timeinc() )
        else
          if( fstrPARAM%contact_algo == kcaSLagrange ) then
            call fstr_Newton_contactSLag( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT,  &
              restart_step_num, restart_substep_num, sub_step, fstr_get_time(), fstr_get_timeinc(), infoCTChange, conMAT )
          else if( fstrPARAM%contact_algo == kcaALagrange ) then
            call fstr_Newton_contactALag( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM,            &
              restart_step_num, restart_substep_num, sub_step, fstr_get_time(), fstr_get_timeinc(), infoCTChange, conMAT )
          endif
        endif

        ! Time Increment
        if( hecMESH%my_rank == 0 ) call fstr_TimeInc_PrintSTATUS( fstrSOLID%step_ctrl(tot_step), fstrPARAM, &
          &  tot_step_print, sub_step, fstrSOLID%NRstat_i, fstrSOLID%NRstat_r,   &
          &  fstrSOLID%AutoINC_stat, fstrSOLID%CutBack_stat )
        if( fstr_cutback_active() ) then

          if( fstrSOLID%CutBack_stat == 0 ) then ! converged
            call fstr_cutback_save( fstrSOLID, infoCTChange, infoCTChange_bak )  ! save analysis state
            call fstr_proceed_time()             ! current time += time increment

          else                                   ! not converged
            CBbound = fstrPARAM%ainc(fstrSOLID%step_ctrl(tot_step)%AincParam_id)%CBbound
            if( fstrSOLID%CutBack_stat == CBbound ) then
              if( hecMESH%my_rank == 0 ) then
                write(*,*) 'Number of successive cutback reached max number: ',CBbound
                call fstr_TimeInc_PrintSTATUS_final(.false.)
              endif
              call hecmw_abort( hecmw_comm_get_comm() )
            endif
            call fstr_cutback_load( fstrSOLID, infoCTChange, infoCTChange_bak )  ! load analysis state
            call fstr_set_contact_active( infoCTChange%contactNode_current > 0 )

            ! restore matrix structure for slagrange contact analysis
            if( is_interaction_active .and. fstrPARAM%contact_algo == kcaSLagrange ) then
              call fstr_mat_con_contact( tot_step, fstrPARAM%contact_algo, hecMAT, fstrSOLID, hecLagMAT, &
                &  infoCTChange, conMAT, fstr_is_contact_active())
              conMAT%B(:) = 0.0d0
              call solve_LINEQ_contact_init(hecMESH, hecMAT, hecLagMAT, fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH))
            endif
            if( hecMESH%my_rank == 0 ) write(*,*) '### State has been restored at time =',fstr_get_time()

            !stop if # of substeps reached upper bound.
            if( sub_step == fstrSOLID%step_ctrl(tot_step)%num_substep ) then
              if( hecMESH%my_rank == 0 ) then
                write(*,'(a,i5,a,f6.3)') '### Number of substeps reached max number: at total_step=', &
                  & tot_step_print, '  time=', fstr_get_time()
              endif
              call hecmw_abort( hecmw_comm_get_comm())
            endif

            ! output time
            time_2 = hecmw_Wtime()
            if( hecMESH%my_rank==0) write(IMSG,'(a,",",2(I8,","),f10.2)') &
              &  'step, substep, solve (sec) :', tot_step_print, sub_step, time_2 - time_1
            cycle
          endif
        else
          if( fstrSOLID%CutBack_stat > 0 ) stop
          call fstr_proceed_time() ! current time += time increment
        endif

        step_count = step_count + 1

        ! ----- Restart
        if( fstrSOLID%restart_nout > 0) then
          if( mod(step_count,fstrSOLID%restart_nout) == 0 ) then
            call fstr_write_restart(tot_step,tot_step_print,sub_step,step_count,fstr_get_time(),  &
              & fstr_get_timeinc_base(), hecMESH,fstrSOLID,fstrPARAM,.false.,infoCTChange%contactNode_current)
          endif
        endif

        ! ----- Result output (include visualize output)
        is_OutPoint = fstr_TimeInc_isTimePoint( fstrSOLID%step_ctrl(tot_step), fstrPARAM ) &
          & .or. fstr_TimeInc_isStepFinished( fstrSOLID%step_ctrl(tot_step) )
        call fstr_static_Output( tot_step, step_count, fstr_get_time(), hecMESH, fstrSOLID, fstrPARAM, &
          &                      fstrPR%solution_type, is_OutPoint )

        time_2 = hecmw_Wtime()
        if( hecMESH%my_rank==0 ) then
          write(IMSG,'(A,",",2(I8,","),f10.2)') 'step, substep, solve (sec) :', tot_step_print, sub_step, time_2 - time_1
          write(IMSG,'(A,I0,",",1pE15.8)') '### stepcount (for output), time :', step_count, fstr_get_time()
        endif

        !if time reached the end time of step, exit loop.
        if( fstr_TimeInc_isStepFinished( fstrSOLID%step_ctrl(tot_step) ) ) exit

        if( sub_step == fstrSOLID%step_ctrl(tot_step)%num_substep ) then
          if( hecMESH%my_rank == 0 ) then
            write(*,'(a,i5,a,f6.3)') '### Number of substeps reached max number: at total_step=', &
              & tot_step_print, '  time=', fstr_get_time()
          endif
          if( hecMESH%my_rank == 0 ) call fstr_TimeInc_PrintSTATUS_final(.false.)
          stop !stop if # of substeps reached upper bound.
        endif

        sub_step = sub_step + 1
      enddo    !--- end of substep  loop

      ! ----- Restart at the end of step
      if( fstrSOLID%restart_nout > 0 ) then
        call fstr_write_restart(tot_step,tot_step_print,sub_step,step_count,fstr_get_time(),fstr_get_timeinc_base(), &
          &  hecMESH,fstrSOLID,fstrPARAM,.true.,infoCTChange%contactNode_current)
      endif
      restart_substep_num = 1
      if( fstrSOLID%TEMP_irres > 0 ) exit
    enddo      !--- end of tot_step loop

    call fstr_cutback_finalize( fstrSOLID )

    !  message
    if(myrank == 0)then
      call fstr_TimeInc_PrintSTATUS_final(.true.)
      write(IMSG,'("### FSTR_SOLVE_NLGEOM FINISHED!")')
      write(*,'("### FSTR_SOLVE_NLGEOM FINISHED!")')
    endif

  end subroutine FSTR_SOLVE_NLGEOM

  !C================================================================C
  !> \brief This subroutine decide the loading increment considering
  !>        the amplitude definition
  !C================================================================C
  subroutine table_nlsta(hecMESH, fstrSOLID, cstep, time, f_t)
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH    !< hecmw mesh
    type ( fstr_solid         ), intent(in) :: fstrSOLID  !< fstr_solid
    integer(kind=kint), intent(in)          :: cstep      !< curr loading step
    real(kind=kreal), intent(in)            :: time       !< loading time(total time)
    real(kind=kreal), intent(out)           :: f_t        !< loading factor

    integer(kind=kint) :: i
    integer(kind=kint) :: jj_n_amp, jj1, jj2
    integer(kind=kint) :: s1, s2, flag
    real(kind=kreal) :: t_1, t_2, t_t, f_1, f_2, tincre

    s1 = 0; s2 = 0
    jj_n_amp = fstrSOLID%step_ctrl( cstep )%amp_id

    if( jj_n_amp <= 0 ) then  ! Amplitude not defined
      f_t = (time-fstrSOLID%step_ctrl(cstep)%starttime)/fstrSOLID%step_ctrl(cstep)%elapsetime
      if( f_t>1.d0 ) f_t=1.d0
    else
      tincre = fstrSOLID%step_ctrl( cstep )%initdt
      jj1 = hecMESH%amp%amp_index(jj_n_amp - 1)
      jj2 = hecMESH%amp%amp_index(jj_n_amp)

      jj1 = jj1 + 2
      t_t = time-fstrSOLID%step_ctrl(cstep)%starttime

      !      if(jj2 .eq. 0) then
      !         f_t = 1.0
      if(t_t .gt. hecMESH%amp%amp_table(jj2)) then
        f_t = hecMESH%amp%amp_val(jj2)
      else if(t_t .le. hecMESH%amp%amp_table(jj2)) then
        flag=0
        do i = jj1, jj2
          if(t_t .le. hecMESH%amp%amp_table(i)) then
            s2 = i
            s1 = i - 1
            flag = 1
          endif
          if( flag == 1 ) exit
        end do

        t_2 = hecMESH%amp%amp_table(s2)
        t_1 = hecMESH%amp%amp_table(s1)
        f_2 = hecMESH%amp%amp_val(s2)
        f_1 = hecMESH%amp%amp_val(s1)
        if( t_2-t_1 .lt. 1.0e-20) then
          if(myrank == 0) then
            write(imsg,*) 'stop due to t_2-t_1 <= 0'
          endif
          call hecmw_abort( hecmw_comm_get_comm())
        endif
        f_t = ((t_2*f_1 - t_1*f_2) + (f_2 - f_1)*t_t) / (t_2 - t_1)
      endif

    endif

  end subroutine table_nlsta

end module m_fstr_solve_NLGEOM

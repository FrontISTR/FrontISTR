!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains subroutines for nonlinear explicit dynamic analysis

module fstr_dynamic_nlexplicit
  use m_fstr
  use m_static_lib
  use m_dynamic_output
  use m_fstr_EIG_setMASS
  use m_dynamic_mat_ass_bc_ac
  use m_dynamic_mat_ass_bc
  use m_dynamic_mat_ass_bc_vl
  use m_dynamic_mat_ass_load
  use m_fstr_Update
  use m_fstr_Residual
  use m_fstr_Restart
  use m_dynamic_mat_ass_couple
  use m_fstr_rcap_io
  use mContact
  use m_fstr_TimeInc

contains

  !C================================================================C
  !C-- subroutine  fstr_solve_LINEAR_DYNAMIC
  !C================================================================C
  subroutine fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
      ,fstrDYN,fstrRESULT,fstrPARAM,infoCTChange &
      ,fstrCPL, restrt_step_num, restrt_step_count )
    implicit none
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_result_data)            :: fstrRESULT
    type(fstr_param)                     :: fstrPARAM
    type(fstr_dynamic)                   :: fstrDYN
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< fstr_info_contactChange
    type(fstr_couple)                    :: fstrCPL !for COUPLE
    type(hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint), allocatable :: mark(:)
    integer(kind=kint) :: nnod, ndof, nn, numnp
    integer(kind=kint) :: i, j, ids, ide, kk
    integer(kind=kint) :: kkk0, kkk1
    integer(kind=kint) :: ierror
    integer(kind=kint) :: iiii5, iexit
    integer(kind=kint) :: revocap_flag
    real(kind=kreal), allocatable :: prevB(:)
    real(kind=kreal) :: bsize, res
    real(kind=kreal) :: time_1, time_2
    integer(kind=kint) :: restrt_step_num
    integer(kind=kint) :: restrt_step_count
    real(kind=kreal), parameter :: PI = 3.14159265358979323846D0
    integer(kind=kint) :: tot_step, sub_step, step_count
    logical :: is_OutPoint

    call hecmw_mpc_mat_init_explicit(hecMESH, hecMAT, hecMATmpc)

    call hecmw_mat_set_NDOF(hecMAT, hecMESH%n_dof)
    nnod=hecMESH%n_node
    ndof=hecmw_mat_get_NDOF(hecMAT)
    nn=ndof*ndof

    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        allocate( prevB(hecmw_mat_get_NP(hecMAT)*ndof)      ,stat=ierror )
        prevB = 0.0d0
        if( ierror /= 0 ) then
          write(idbg,*) 'stop due to allocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
          write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
        endif
      endif
    endif

    fstrSOLID%dunode(:) =0.d0

    call fstr_prepare_dynamic_explicit( hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYN, &
      & ndof, nnod, restrt_step_count )

    if( restrt_step_count == 0 ) then
      call fstr_dynamic_Output(1, 0, 0.d0, hecMESH, fstrSOLID, fstrDYN, fstrPARAM, .true.)
      call dynamic_output_monit(1, 0, 0.d0, hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)
    end if

    if( associated( fstrSOLID%contacts ) )  then
      call initialize_contact_output_vectors(fstrSOLID,hecMAT)
      call forward_increment_Lagrange(1,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
        & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
    endif

    step_count = restrt_step_count
    do tot_step = 1, fstrSOLID%nstep_tot
      if(hecMESH%my_rank==0) write(*,'(a,i5)') ' loading step=',tot_step

      sub_step = restrt_step_num
      do while(.true.)
        call fstr_TimeInc_SetTimeIncrement( fstrSOLID%step_ctrl(tot_step), fstrPARAM, sub_step, &
          & fstrSOLID%NRstat_i, fstrSOLID%NRstat_r, fstrSOLID%AutoINC_stat, fstrSOLID%CutBack_stat )

        fstrDYN%t_curr = fstr_get_time()
        fstrDYN%t_delta = fstr_get_timeinc()

        step_count = step_count + 1

        call fstr_advance_dynamic_explicit( tot_step, sub_step, &
            hecMESH, hecMAT, hecMATmpc, fstrSOLID, fstrEIG, fstrDYN, fstrPARAM, &
            fstrCPL, infoCTChange, &
            restrt_step_num, ndof, nnod, prevB, &
            fstr_TimeInc_isStepFinished( fstrSOLID%step_ctrl(tot_step) ) )

        ! ----- Result output (include visualize output)
        ! Evaluate isTimePoint before time advance, then OR with isStepFinished after
        is_OutPoint = fstr_TimeInc_isTimePoint( fstrSOLID%step_ctrl(tot_step), fstrPARAM )

          !C-- output result of monitoring node
        call dynamic_output_monit(tot_step, sub_step, fstrDYN%t_curr, hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)

        if( fstrDYN%restart_nout > 0 ) then
          if ( mod(step_count,fstrDYN%restart_nout).eq.0 ) then
            call fstr_write_restart_dyna_nl(tot_step,sub_step,hecMESH,fstrSOLID,fstrDYN,fstrPARAM,&
              .false.,infoCTChange%contactNode_current,step_count)
          end if
        end if

        call fstr_proceed_time()
        fstrDYN%t_curr = fstr_get_time()

        ! isStepFinished must be evaluated after fstr_proceed_time
        is_OutPoint = is_OutPoint .or. fstr_TimeInc_isStepFinished( fstrSOLID%step_ctrl(tot_step) )

          !C-- output new displacement, velocity and acceleration
        call fstr_dynamic_Output(tot_step, step_count, fstrDYN%t_curr, hecMESH, fstrSOLID, fstrDYN, fstrPARAM, is_OutPoint)

        if( fstr_TimeInc_isStepFinished( fstrSOLID%step_ctrl(tot_step) ) ) exit

        if( sub_step == fstrSOLID%step_ctrl(tot_step)%num_substep ) then
          if( hecMESH%my_rank == 0 ) then
            write(*,'(a,i5,a,f6.3)') '### Number of substeps reached max number: at total_step=', &
              & tot_step, '  time=', fstr_get_time()
          endif
          call hecmw_abort( hecmw_comm_get_comm())
        endif

        sub_step = sub_step + 1
      enddo

      if( fstrDYN%restart_nout > 0 ) then
        call fstr_write_restart_dyna_nl(tot_step,sub_step,hecMESH,fstrSOLID,fstrDYN,fstrPARAM,&
          .true.,infoCTChange%contactNode_current,step_count)
      end if
      restrt_step_num = 1
    enddo

    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        deallocate( prevB      ,stat=ierror )
        if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
          write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
        endif
      endif
    endif

    call hecmw_mpc_mat_finalize_explicit(hecMESH, hecMAT, hecMATmpc)

  end subroutine fstr_solve_dynamic_nlexplicit

  !> \brief Prepare initial state for explicit dynamic analysis (central difference).
  !!
  !! Computes the diagonal effective mass VEC1 = (1/dt^2 + ray_m/(2*dt)) * mass
  !! and the virtual past displacements DISP(:,3) = u(-dt), DISP(:,2) = u(-2*dt)
  !! used to start the central-difference scheme. All scheme-specific
  !! coefficients are local to this subroutine.
  subroutine fstr_prepare_dynamic_explicit( hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYN, &
      ndof, nnod, restrt_step_count )
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(inout)     :: hecMAT
    type(fstr_solid), intent(inout)         :: fstrSOLID
    type(fstr_eigen), intent(inout)         :: fstrEIG
    type(fstr_dynamic), intent(inout)       :: fstrDYN
    integer(kind=kint), intent(in)          :: ndof
    integer(kind=kint), intent(in)          :: nnod
    integer(kind=kint), intent(in)          :: restrt_step_count

    integer(kind=kint), allocatable :: mark(:)
    integer(kind=kint) :: j
    real(kind=kreal)   :: a1, a2

    a1 = 1.d0/fstrDYN%t_delta**2
    a2 = 1.d0/(2.d0*fstrDYN%t_delta)

    call setMASS(fstrSOLID,hecMESH,hecMAT,fstrEIG)
    call hecmw_mpc_trans_mass(hecMESH, hecMAT, fstrEIG%mass)

    allocate(mark(hecmw_mat_get_NP(hecMAT) * hecmw_mat_get_NDOF(hecMAT)))
    call hecmw_mpc_mark_slave(hecMESH, hecMAT, mark)

    do j = 1 ,ndof*nnod
      fstrDYN%VEC1(j) = (a1 + a2 *fstrDYN%ray_m) * fstrEIG%mass(j)
      if(mark(j) == 1) fstrDYN%VEC1(j) = 1.d0
      if(dabs(fstrDYN%VEC1(j)) < 1.0e-20) then
        if( hecMESH%my_rank == 0 ) then
          write(*,*) 'stop due to fstrDYN%VEC(j) = 0 ,  j = ', j
          write(imsg,*) 'stop due to fstrDYN%VEC(j) = 0 ,  j = ', j
        end if
        call hecmw_abort( hecmw_comm_get_comm())
      endif
    end do

    deallocate(mark)

    !C-- virtual past displacements for central-difference startup
    if( restrt_step_count == 0 ) then
      do j = 1 ,ndof*nnod
        fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1) - fstrDYN%VEL (j,1)/(2.d0*a2)  + fstrDYN%ACC (j,1)/ (2.d0*a1)
        fstrDYN%DISP(j,2) = fstrDYN%DISP(j,1) - fstrDYN%VEL (j,1)/ a2 + fstrDYN%ACC (j,1)/ (2.d0*a1) * 4.d0
      end do
    endif

  end subroutine fstr_prepare_dynamic_explicit

  !> \brief Advance one time step of explicit dynamic analysis
  subroutine fstr_advance_dynamic_explicit( cstep, istep, &
      hecMESH, hecMAT, hecMATmpc, fstrSOLID, fstrEIG, fstrDYN, fstrPARAM, &
      fstrCPL, infoCTChange, &
      restrt_step_num, ndof, nnod, prevB, &
      is_last_step )
    implicit none
    integer(kind=kint), intent(in)               :: cstep           !< current loading step
    integer(kind=kint), intent(in)               :: istep           !< current sub-step (time step index)
    type(hecmwST_local_mesh), intent(inout)      :: hecMESH
    type(hecmwST_matrix), intent(inout)          :: hecMAT
    type(hecmwST_matrix), pointer, intent(inout) :: hecMATmpc
    type(fstr_solid), intent(inout)              :: fstrSOLID
    type(fstr_eigen), intent(inout)              :: fstrEIG
    type(fstr_dynamic), intent(inout)            :: fstrDYN
    type(fstr_param), intent(inout)              :: fstrPARAM
    type(fstr_couple), intent(inout)             :: fstrCPL
    type(fstr_info_contactChange), intent(inout) :: infoCTChange
    integer(kind=kint), intent(in)               :: restrt_step_num
    integer(kind=kint), intent(in)               :: ndof
    integer(kind=kint), intent(in)               :: nnod
    real(kind=kreal), allocatable, intent(inout) :: prevB(:)
    logical, intent(in)                          :: is_last_step

    integer(kind=kint) :: j, kk, kkk0, kkk1
    integer(kind=kint) :: revocap_flag
    real(kind=kreal)   :: bsize
    real(kind=kreal), parameter :: PI = 3.14159265358979323846D0
    real(kind=kreal)   :: a1, a2
    real(kind=kreal)   :: b1, b2, b3, a3
    real(kind=kreal), pointer :: pB(:), pX(:)

    !C-- central-difference coefficients (depend only on dt)
    a1 = 1.d0/fstrDYN%t_delta**2
    a2 = 1.d0/(2.d0*fstrDYN%t_delta)
    !C-- coupling-only coefficients (kept zero for central difference)
    a3 = 0.d0
    b1 = 0.d0; b2 = 0.d0; b3 = 0.d0

    !C-- mechanical boundary condition
    call dynamic_mat_ass_load (cstep, fstrDYN%t_curr, hecMESH, hecMAT, fstrSOLID, fstrDYN, fstrPARAM)
    pB => hecmw_mat_get_B(hecMAT)
    do j=1, hecMESH%n_node*  hecMESH%n_dof
      pB(j) = pB(j) - fstrSOLID%QFORCE(j)
    end do

    !C ********************************************************************************
    !C for couple analysis
    if( fstrPARAM%fg_couple == 1 ) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        pB => hecmw_mat_get_B(hecMAT)
        do j = 1, hecmw_mat_get_NP(hecMAT) * ndof
          prevB(j) = pB(j)
        enddo
      endif
    endif
    do
      if( fstrPARAM%fg_couple == 1 ) then
        if( fstrPARAM%fg_couple_type==1 .or. &
          fstrPARAM%fg_couple_type==3 .or. &
          fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
        if( fstrPARAM%fg_couple_first /= 0 ) then
          bsize = dfloat( istep ) / dfloat( fstrPARAM%fg_couple_first )
          if( bsize > 1.0 ) bsize = 1.0
          do kkk0 = 1, fstrCPL%coupled_node_n
            kkk1 = 3 * kkk0
            fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
            fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
            fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
          enddo
        endif
        if( fstrPARAM%fg_couple_window > 0 ) then
          j = istep - restrt_step_num + 1
          kk = fstrDYN%n_step - restrt_step_num + 1
          bsize = 0.5*(1.0-cos(2.0*PI*dfloat(j)/dfloat(kk)))
          do kkk0 = 1, fstrCPL%coupled_node_n
            kkk1 = 3 * kkk0
            fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
            fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
            fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
          enddo
        endif
        call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
      endif
      !C ********************************************************************************

      call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)

      pB => hecmw_mat_get_B(hecMATmpc)
      do j = 1 ,ndof*nnod
        pB(j) = pB(j) + 2.d0*a1* fstrEIG%mass(j) * fstrDYN%DISP(j,1)  &
          + (- a1 + a2 * fstrDYN%ray_m) * fstrEIG%mass(j) * fstrDYN%DISP(j,3)
      end do

      !C
      !C-- geometrical boundary condition

      call dynamic_explicit_ass_bc(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrDYN%t_curr)
      call dynamic_explicit_ass_vl(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrDYN%t_curr)
      call dynamic_explicit_ass_ac(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrDYN%t_curr)

      ! Finish the calculation
      pB => hecmw_mat_get_B(hecMATmpc)
      pX => hecmw_mat_get_X(hecMATmpc)
      do j = 1 ,ndof*nnod
        pX(j) = pB(j) / fstrDYN%VEC1(j)
        if(dabs(pX(j)) > 1.0d+5) then
          if( hecMESH%my_rank == 0 ) then
            print *, 'Displacement increment too large, please adjust your step size!',istep,pX(j)
            write(imsg,*) 'Displacement increment too large, please adjust your step size!', &
              & istep,pB(j),fstrDYN%VEC1(j)
          end if
          call hecmw_abort( hecmw_comm_get_comm())
        end if
      end do
      call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)

      !C *****************************************************
      !C for couple analysis
      if( fstrPARAM%fg_couple == 1 ) then
        if( fstrPARAM%fg_couple_type>1 ) then
          pX => hecmw_mat_get_X(hecMAT)
          do j=1, fstrCPL%coupled_node_n
            if( fstrCPL%dof == 3 ) then
              kkk0 = j*3
              kkk1 = fstrCPL%coupled_node(j)*3

              fstrCPL%disp (kkk0-2) = pX(kkk1-2)
              fstrCPL%disp (kkk0-1) = pX(kkk1-1)
              fstrCPL%disp (kkk0  ) = pX(kkk1)

              fstrCPL%velo (kkk0-2) = -b1*fstrDYN%ACC(kkk1-2,1) - b2*fstrDYN%VEL(kkk1-2,1) + &
                b3*( pX(kkk1-2) - fstrDYN%DISP(kkk1-2,1) )
              fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
                b3*( pX(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
              fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
                b3*( pX(kkk1) - fstrDYN%DISP(kkk1,1) )
              fstrCPL%accel(kkk0-2) = -a1*fstrDYN%ACC(kkk1-2,1) - a2*fstrDYN%VEL(kkk1-2,1) + &
                a3*( pX(kkk1-2) - fstrDYN%DISP(kkk1-2,1) )
              fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
                a3*( pX(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
              fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
                a3*( pX(kkk1) - fstrDYN%DISP(kkk1,1) )
            else
              kkk0 = j*2
              kkk1 = fstrCPL%coupled_node(j)*2

              fstrCPL%disp (kkk0-1) = pX(kkk1-1)
              fstrCPL%disp (kkk0  ) = pX(kkk1)

              fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
                b3*( pX(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
              fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
                b3*( pX(kkk1) - fstrDYN%DISP(kkk1,1) )
              fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
                a3*( pX(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
              fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
                a3*( pX(kkk1) - fstrDYN%DISP(kkk1,1) )
            endif
          end do
          call fstr_rcap_send( fstrCPL )
        endif

        select case ( fstrPARAM%fg_couple_type )
          case (4)
            call fstr_rcap_get( fstrCPL )
          case (5)
            call fstr_get_convergence( revocap_flag )
            if( revocap_flag==0 ) then
              pB => hecmw_mat_get_B(hecMAT)
              do j = 1, hecmw_mat_get_NP(hecMAT) * ndof
                pB(j) = prevB(j)
              enddo
              cycle
            endif
          case (6)
            call fstr_get_convergence( revocap_flag )
            if( revocap_flag==0 ) then
              pB => hecmw_mat_get_B(hecMAT)
              do j = 1, hecmw_mat_get_NP(hecMAT) * ndof
                pB(j) = prevB(j)
              enddo
              call fstr_rcap_get( fstrCPL )
              cycle
            else
              if( .not. is_last_step ) call fstr_rcap_get( fstrCPL )
            endif
        end select
      endif
      exit
    enddo

    !C *****************************************************
    !C-- contact corrector
    !C
    pX => hecmw_mat_get_X(hecMAT)
    do j = 1 ,ndof*nnod
      fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
      fstrSOLID%dunode(j)  = pX(j)-fstrDYN%DISP(j,1)
    enddo
    if( associated( fstrSOLID%contacts ) )  then
      !call fstr_scan_contact_state( cstep, fstrDYN%t_delta, kcaSLAGRANGE, hecMESH, fstrSOLID, infoCTChange )
      call forward_increment_Lagrange(cstep,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
        & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
      pX => hecmw_mat_get_X(hecMAT)
      do j = 1 ,ndof*nnod
        pX(j) = pX(j) + fstrSOLID%ddunode(j)
      enddo
    endif

    !C-- new displacement, velocity and acceleration
    pX => hecmw_mat_get_X(hecMAT)
    do j = 1 ,ndof*nnod
      fstrDYN%ACC (j,1) = a1*(pX(j) - 2.d0*fstrDYN%DISP(j,1) + fstrDYN%DISP(j,3))
      fstrDYN%VEL (j,1) = a2*(pX(j) - fstrDYN%DISP(j,3))
      fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
      fstrSOLID%dunode(j)  = pX(j)-fstrDYN%DISP(j,1)
      fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1)
      fstrDYN%DISP(j,1) = pX(j)
      pX(j) = fstrSOLID%dunode(j)
      fstrDYN%kineticEnergy = fstrDYN%kineticEnergy + 0.5d0*fstrEIG%mass(j)*fstrDYN%VEL(j,1)*fstrDYN%VEL(j,1)
    end do

    ! ----- update strain, stress, and internal force
    call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID, fstrDYN%t_curr, fstrDYN%t_delta, 0, fstrDYN%strainEnergy )

    ! ----- update reaction force at constrained DOFs using converged QFORCE
    call fstr_Update_REACTION_SPC( cstep, hecMESH, fstrSOLID )

    do j = 1 ,ndof*nnod
      fstrSOLID%unode(j) = fstrSOLID%unode(j) + fstrSOLID%dunode(j)
    end do
    call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYN%t_delta )

  end subroutine fstr_advance_dynamic_explicit

  !< This subroutine implements Forward increment Lagrange multiplier method( NJ Carpenter et al. Int.J.Num.Meth.Eng.,32(1991),103-128 )
  subroutine forward_increment_Lagrange(cstep,ndof,mmat,hecMESH,fstrSOLID,infoCTChange,wkarray,uc)
    integer, intent(in)                    :: cstep
    integer, intent(in)                    :: ndof
    real(kind=kreal), intent(in)           :: mmat(:)
    type( hecmwST_local_mesh ), intent(in) :: hecMESH       !< type mesh
    type(fstr_solid), intent(inout)        :: fstrSOLID
    type(fstr_info_contactChange)          :: infoCTChange
    real(kind=kreal), intent(out)          :: wkarray(:)
    real(kind=kreal), intent(out)          :: uc(:)
    integer :: i, j, k, m, grpid, slave, nn, iSS, sid, etype, iter
    real(kind=kreal) :: fdum, conv, dlambda, shapefunc(l_max_surface_node), lambda(3)

    call fstr_scan_contact_state_exp( cstep, hecMESH, fstrSOLID, infoCTChange )
    if( .not. infoCTChange%active ) return

    uc = 0.0d0

    iter = 0
    do
      wkarray = 0.0d0
      do i=1,fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          if( .not. is_contact_active(fstrSOLID%contacts(i)%states(j)%state) ) cycle
          if( fstrSOLID%contacts(i)%states(j)%distance>epsilon(1.d0) ) then
            fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
            cycle
          endif
          if( iter==0 ) then
            fstrSOLID%contacts(i)%states(j)%multiplier(:) =0.d0
            fstrSOLID%contacts(i)%states(j)%wkdist =0.d0
            cycle
          endif
          slave = fstrSOLID%contacts(i)%slave(j)

          sid = fstrSOLID%contacts(i)%states(j)%surface
          nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
          etype = fstrSOLID%contacts(i)%master(sid)%etype
          call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
          wkarray( slave ) = -fstrSOLID%contacts(i)%states(j)%multiplier(1)
          do k=1,nn
            iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
            wkarray( iSS ) = wkarray( iSS ) + shapefunc(k) * fstrSOLID%contacts(i)%states(j)%multiplier(1)
          enddo
        enddo
      enddo

      if(iter > 0)then
        do i=1,fstrSOLID%n_contacts
          do j= 1, size(fstrSOLID%contacts(i)%slave)
            if( .not. is_contact_active(fstrSOLID%contacts(i)%states(j)%state) ) cycle
            slave = fstrSOLID%contacts(i)%slave(j)
            sid = fstrSOLID%contacts(i)%states(j)%surface
            nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
            etype = fstrSOLID%contacts(i)%master(sid)%etype
            call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
            fstrSOLID%contacts(i)%states(j)%wkdist = -wkarray( slave )/mmat( (slave-1)*ndof+1 )
            do k=1,nn
              iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
              fstrSOLID%contacts(i)%states(j)%wkdist = fstrSOLID%contacts(i)%states(j)%wkdist  &
                   + shapefunc(k) * wkarray(iSS) / mmat( (iSS-1)*ndof+1 )
            enddo
          enddo
        enddo
      endif

      conv = 0.d0
      wkarray = 0.d0
      do i=1,fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          if( .not. is_contact_active(fstrSOLID%contacts(i)%states(j)%state) ) cycle
          slave = fstrSOLID%contacts(i)%slave(j)
          sid = fstrSOLID%contacts(i)%states(j)%surface
          nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
          etype = fstrSOLID%contacts(i)%master(sid)%etype
          call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
          fdum = 1.d0/mmat( (slave-1)*ndof+1 )
          do k=1,nn
            iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
            fdum = fdum + shapefunc(k)*shapefunc(k)/mmat( (iSS-1)*ndof+1 )
          enddo
          dlambda= (fstrSOLID%contacts(i)%states(j)%distance-fstrSOLID%contacts(i)%states(j)%wkdist) /fdum
          conv = conv + dlambda*dlambda;
          fstrSOLID%contacts(i)%states(j)%multiplier(1) = fstrSOLID%contacts(i)%states(j)%multiplier(1) + dlambda
          if( fstrSOLID%contacts(i)%fcoeff>0.d0 ) then
            if( fstrSOLID%contacts(i)%states(j)%state == CONTACTSLIP ) then
              fstrSOLID%contacts(i)%states(j)%multiplier(2) =             &
              fstrSOLID%contacts(i)%fcoeff * fstrSOLID%contacts(i)%states(j)%multiplier(1)
            else    ! stick
              !      fstrSOLID%contacts(i)%states(j)%multiplier(2) =
            endif
          endif
          lambda = fstrSOLID%contacts(i)%states(j)%multiplier(1)* fstrSOLID%contacts(i)%states(j)%direction
          wkarray((slave-1)*ndof+1:(slave-1)*ndof+3) = lambda(:)
          do k=1,nn
            iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
            wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) = wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) -lambda(:)*shapefunc(k)
          enddo
        enddo
      enddo
      if( dsqrt(conv)<1.d-8 ) exit
      iter = iter+1
    enddo

    do i=1,hecMESH%n_node*ndof
      uc(i) = wkarray(i)/mmat(i)
    enddo
  end subroutine forward_increment_Lagrange

end module fstr_dynamic_nlexplicit

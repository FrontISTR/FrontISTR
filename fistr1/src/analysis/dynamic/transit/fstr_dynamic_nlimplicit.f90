!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains subroutines for nonlinear implicit dynamic analysis

module fstr_dynamic_nlimplicit
  use m_fstr
  use m_dynamic_output
  use m_fstr_EIG_setMASS
  use m_dynamic_mat_ass_bc_ac
  use m_dynamic_mat_ass_bc
  use m_dynamic_mat_ass_bc_vl
  use m_dynamic_mat_ass_load
  use m_fstr_StiffMatrix
  use m_fstr_Update
  use m_fstr_Restart
  use fstr_matrix_con_contact
  use m_fstr_Residual
  use mContact
  use m_addContactStiffness
  use m_solve_LINEQ_contact
  use m_dynamic_init_variables
  use m_fstr_TimeInc
  use m_fstr_Cutback

  !-------- for couple -------
  use m_dynamic_mat_ass_couple
  use m_fstr_rcap_io

  real(kind=kreal), parameter :: PI = 3.14159265358979323846D0

contains

  !> \brief This subroutine provides function of nonlinear implicit dynamic analysis using the Newmark method.
  !> Standard Lagrange multiplier algorithm for contact analysis is included in this subroutine.
  subroutine FSTR_SOLVE_NLGEOM_DYNAMIC_IMPLICIT_CONTACTSLAG(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
      ,fstrCPL,hecLagMAT,restart_step_num,restart_substep_num,infoCTChange  &
      ,conMAT )
    implicit none
    !C-- global variable
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(hecmwST_matrix), pointer        :: hecMAT0
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_result_data)            :: fstrRESULT ! unused - kept for interface compatibility
    type(fstr_param)                     :: fstrPARAM
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_couple)                    :: fstrCPL !for COUPLE
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< fstr_info_contactChange
    type(hecmwST_matrix)                 :: conMAT

    !C-- local variable
    integer(kind=kint) :: nnod, ndof, nn
    integer(kind=kint) :: i, tot_step_print, CBbound
    real(kind=kreal) :: time_1, time_2, factor
    integer(kind=kint) :: sub_step

    integer(kind=kint) :: restart_step_num, restart_substep_num, tot_step, step_count
    integer(kind=kint) :: ctAlgo
    integer(kind=kint) :: max_iter_contact
    real(kind=kreal) :: converg_dlag
    type(fstr_info_contactChange)        :: infoCTChange_bak

    integer(kind=kint) :: n_node_global
    logical :: is_mat_symmetric, is_interaction_active

    if(hecMESH%my_rank==0) call fstr_TimeInc_PrintSTATUS_init

    is_interaction_active = ( associated( fstrSOLID%contacts ) .or. associated( fstrSOLID%embeds ) )

    nullify(hecMAT0)

    ! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

    ctAlgo = fstrPARAM%contact_algo

    if( hecMAT%Iarray(99)==4 .and. .not.fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH) ) then
      write(*,*) ' This type of direct solver is not yet available in such case ! '
      write(*,*) ' Please use intel MKL direct solver !'
      call  hecmw_abort(hecmw_comm_get_comm())
    endif

    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof

    if( associated( fstrSOLID%contacts ) ) call initialize_contact_output_vectors(fstrSOLID,hecMAT)

    !!-- initial value
    time_1 = hecmw_Wtime()

    !C-- check parameters
    if(dabs(fstrDYNAMIC%beta) < 1.0e-20) then
      if( hecMESH%my_rank == 0 ) then
        write(imsg,*) 'stop due to Newmark-beta = 0'
      endif
      call hecmw_abort( hecmw_comm_get_comm())
    endif

    !C-- matrix [M] lumped mass matrix
    if(fstrDYNAMIC%idx_mas == 1) then
      call setMASS(fstrSOLID,hecMESH,hecMAT,fstrEIG)

    !C-- consistent mass matrix
    else if(fstrDYNAMIC%idx_mas == 2) then
      if( hecMESH%my_rank .eq. 0 ) then
        write(imsg,*) 'stop: consistent mass matrix is not yet available !'
      endif
      call hecmw_abort( hecmw_comm_get_comm())
    endif

    hecMAT%Iarray(98) = 1   !Assembly complete
    hecMAT%Iarray(97) = 1   !Need numerical factorization

    !C-- initialize variables
    if( restart_step_num == 1 .and. fstrDYNAMIC%VarInitialize .and. abs(fstrDYNAMIC%ray_m) > 1.0d-15 ) &
      call dynamic_init_varibles( hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYNAMIC, fstrPARAM )

    !C-- output of initial state
    if( restart_step_num == 1 ) then
      call fstr_dynamic_Output(1, 0, 0, 0.d0, hecMESH, fstrSOLID, fstrDYNAMIC, fstrPARAM)
      call dynamic_output_monit(1, 0, 0.d0, hecMESH, fstrPARAM, fstrDYNAMIC, fstrEIG, fstrSOLID)
    endif

    fstrDYNAMIC%VEC3(:) =0.d0
    hecMAT%X(:) =0.d0

    call fstr_save_originalMatrixStructure(hecMAT)
    call fstr_scan_contact_state(restart_step_num, restart_step_num, 0, fstrDYNAMIC%t_delta, &
          ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B)

    call hecmw_mat_copy_profile( hecMAT, conMAT )

    if ( fstr_is_contact_active() ) then
      call fstr_mat_con_contact( restart_step_num, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, &
          infoCTChange, conMAT, fstr_is_contact_active())
    elseif( hecMAT%Iarray(99)==4 ) then
      write(*,*) ' This type of direct solver is not yet available in such case ! '
      write(*,*) ' Please change solver type to intel MKL direct solver !'
      call  hecmw_abort(hecmw_comm_get_comm())
    endif
    is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)
    call solve_LINEQ_contact_init(hecMESH,hecMAT,hecLagMAT,is_mat_symmetric)

    fstrSOLID%FACTOR = 0.0d0
    call fstr_cutback_init( hecMESH, fstrSOLID, fstrPARAM )
    call fstr_cutback_save( fstrSOLID, infoCTChange, infoCTChange_bak )
    
    step_count = 0
    do tot_step=restart_step_num, fstrSOLID%nstep_tot
      tot_step_print = tot_step+restart_step_num-1
      if(hecMESH%my_rank==0) write(*,*) ''
      if(hecMESH%my_rank==0) write(*,'(a,i5)') ' loading step=',tot_step_print

      sub_step = restart_substep_num
      do while(.true.)
        ! ----- time history of factor
        call fstr_TimeInc_SetTimeIncrement( fstrSOLID%step_ctrl(tot_step), fstrPARAM, sub_step, &
          &  fstrSOLID%NRstat_i, fstrSOLID%NRstat_r, fstrSOLID%AutoINC_stat, fstrSOLID%CutBack_stat )

        fstrDYNAMIC%t_curr = fstr_get_time()
        fstrDYNAMIC%t_delta = fstr_get_timeinc()

        call fstr_Newton_dynamic_contactSLag(tot_step, hecMESH, hecMAT, fstrSOLID, fstrEIG, &
            fstrDYNAMIC, fstrPARAM, fstrCPL, hecLagMAT, infoCTChange, conMAT, &
            restart_step_num, hecMAT0, sub_step, fstrDYNAMIC%t_curr, fstrDYNAMIC%t_delta)

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

        !C-- output new displacement, velocity and acceleration
        call fstr_dynamic_Output(tot_step, i, step_count, fstrDYNAMIC%t_curr, hecMESH, fstrSOLID, fstrDYNAMIC, fstrPARAM)

        !C-- output result of monitoring node
        call dynamic_output_monit(tot_step, i, fstrDYNAMIC%t_curr, hecMESH, fstrPARAM, fstrDYNAMIC, fstrEIG, fstrSOLID)

        !---  Restart info
        if( fstrDYNAMIC%restart_nout > 0 ) then
          if( mod(i,fstrDYNAMIC%restart_nout).eq.0 .or. i.eq.fstrDYNAMIC%n_step ) then
            call fstr_write_restart_dyna_nl(tot_step,i,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
              infoCTChange%contactNode_current)
          endif
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
      enddo
      restart_substep_num = 1
      !C-- end of time step loop
    enddo

    if (associated(hecMAT0)) then
      call hecmw_mat_finalize(hecMAT0)
      deallocate(hecMAT0)
    endif

    time_2 = hecmw_Wtime()
    if( hecMESH%my_rank == 0 ) then
      write(ISTA,'(a,f10.2,a)') '         solve (sec) :', time_2 - time_1, 's'
    endif

  end subroutine FSTR_SOLVE_NLGEOM_DYNAMIC_IMPLICIT_CONTACTSLAG

  subroutine fstr_Newton_dynamic_contactSLag(cstep, hecMESH, hecMAT, fstrSOLID, fstrEIG, &
      fstrDYNAMIC, fstrPARAM, fstrCPL, hecLagMAT, infoCTChange, conMAT, &
      restart_step_num, hecMAT0, istep, t_curr, t_delta)
    implicit none
    !C-- arguments
    integer(kind=kint), intent(in)       :: cstep, restart_step_num, istep
    real(kind=kreal), intent(in)         :: t_curr, t_delta
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(hecmwST_matrix), pointer        :: hecMAT0
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_param)                     :: fstrPARAM
    type(fstr_couple)                    :: fstrCPL
    type(hecmwST_matrix_lagrange)        :: hecLagMAT
    type(fstr_info_contactChange)        :: infoCTChange
    type(hecmwST_matrix)                 :: conMAT

    !C-- local variables
    integer(kind=kint) :: j, kk, idm, imm
    integer(kind=kint) :: iter
    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: res, res1, res0, relres
    integer(kind=kint) :: count_step, stepcnt
    real(kind=kreal) :: maxDLag
    integer(kind=kint) :: contact_changed_global
    logical :: is_mat_symmetric
    integer :: istat
    logical :: is_cycle
    integer(kind=kint) :: ctAlgo, max_iter_contact
    integer(kind=kint) :: nnod, ndof, nn
    real(kind=kreal) :: converg_dlag
    real(kind=kreal), allocatable :: coord(:)

    fstrSOLID%NRstat_i(:) = 0 ! logging newton iteration(init)

    !C-- initialize local variables
    ctAlgo = fstrPARAM%contact_algo
    max_iter_contact = fstrSOLID%step_ctrl(cstep)%max_contiter
    converg_dlag = fstrSOLID%step_ctrl(cstep)%converg_lag
    nnod = hecMESH%n_node
    ndof = hecMAT%NDOF
    nn = ndof*ndof
    allocate(coord(hecMESH%n_node*ndof))

    a1 = .5d0/fstrDYNAMIC%beta - 1.d0
    a2 = 1.d0/(fstrDYNAMIC%beta*t_delta)
    a3 = 1.d0/(fstrDYNAMIC%beta*t_delta*t_delta)
    b1 = ( .5d0*fstrDYNAMIC%gamma/fstrDYNAMIC%beta - 1.d0 )*t_delta
    b2 = fstrDYNAMIC%gamma/fstrDYNAMIC%beta - 1.d0
    b3 = fstrDYNAMIC%gamma/(fstrDYNAMIC%beta*t_delta)
    c1 = 1.d0 + fstrDYNAMIC%ray_k*b3
    c2 = a3 + fstrDYNAMIC%ray_m*b3

    if(hecMESH%my_rank==0) then
      write(*,'(A)')'-------------------------------------------------'
      write(*,'('' time step='',i10,'' time='',1pe13.4e3)') istep,t_curr
    endif

    do j = 1 ,ndof*nnod
      fstrDYNAMIC%VEC1(j) = a1*fstrDYNAMIC%ACC(j,1) + a2*fstrDYNAMIC%VEL(j,1)
      fstrDYNAMIC%VEC2(j) = b1*fstrDYNAMIC%ACC(j,1) + b2*fstrDYNAMIC%VEL(j,1)
    enddo

    count_step = 0
    stepcnt = 0

    !C for couple analysis
    do
      fstrSOLID%dunode(:) =0.d0
      ! call fstr_UpdateEPState( hecMESH, fstrSOLID )
      call fstr_solve_dynamic_nlimplicit_couple_init(fstrPARAM, fstrCPL)

    loopFORcontactAnalysis: do while( .TRUE. )
    count_step = count_step + 1

      ! ----- Inner Iteration
      res0   = 0.d0
      res1   = 0.d0
      relres = 1.d0

      do iter = 1, fstrSOLID%step_ctrl(cstep)%max_iter
        stepcnt=stepcnt+1
        if (fstrPARAM%nlgeom) then
          call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, t_curr, t_delta )
        else
          if (.not. associated(hecMAT0)) then
            call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, t_curr, t_delta )
            allocate(hecMAT0)
            call hecmw_mat_init(hecMAT0)
            call hecmw_mat_copy_profile(hecMAT, hecMAT0)
            call hecmw_mat_copy_val(hecMAT, hecMAT0)
          else
            call hecmw_mat_copy_val(hecMAT0, hecMAT)
          endif
        endif

        if( abs(fstrDYNAMIC%ray_k) > 1.0d-15 .or. abs(fstrDYNAMIC%ray_m) > 1.0d-15 ) then
          do j = 1 ,ndof*nnod
            hecMAT%X(j) = fstrDYNAMIC%VEC2(j) - b3*fstrSOLID%dunode(j)
          enddo
        endif
        if( abs(fstrDYNAMIC%ray_k) > 1.0d-15 ) then
          if( hecMESH%n_dof == 3 ) then
            call hecmw_matvec (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3)
          else if( hecMESH%n_dof == 2 ) then
            call hecmw_matvec (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3)
          else if( hecMESH%n_dof == 6 ) then
            call matvec(fstrDYNAMIC%VEC3, hecMAT%X, hecMAT, ndof, hecMAT%D, hecMAT%AU, hecMAT%AL)
          endif
        endif

        !C-- mechanical boundary condition
        call dynamic_mat_ass_load (cstep, t_curr, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, iter )
        do j=1, hecMESH%n_node*  hecMESH%n_dof
          hecMAT%B(j)=hecMAT%B(j)- fstrSOLID%QFORCE(j) + fstrEIG%mass(j)*( fstrDYNAMIC%VEC1(j)-a3*fstrSOLID%dunode(j)   &
            + fstrDYNAMIC%ray_m* hecMAT%X(j) ) + fstrDYNAMIC%ray_k*fstrDYNAMIC%VEC3(j)
        enddo

        !C for couple analysis
        call fstr_solve_dynamic_nlimplicit_couple_pre(hecMESH, hecMAT, fstrSOLID, &
          & fstrPARAM, fstrDYNAMIC, fstrCPL, restart_step_num, istep)

        do j = 1 ,nn*hecMAT%NP
          hecMAT%D(j)  = c1* hecMAT%D(j)
        enddo
        do j = 1 ,nn*hecMAT%NPU
          hecMAT%AU(j) = c1* hecMAT%AU(j)
        enddo
        do j = 1 ,nn*hecMAT%NPL
          hecMAT%AL(j) = c1*hecMAT%AL(j)
        enddo
        do j=1,nnod
          do kk=1,ndof
            idm = nn*(j-1)+1 + (ndof+1)*(kk-1)
            imm = ndof*(j-1) + kk
            hecMAT%D(idm) = hecMAT%D(idm) + c2*fstrEIG%mass(imm)
          enddo
        enddo

        call hecmw_mat_clear( conMAT )
        call hecmw_mat_clear_b( conMAT )
        conMAT%X = 0.0d0

        if( fstr_is_contact_active() ) then
          call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
          call fstr_AddContactStiffness(cstep,iter,conMAT,hecLagMAT,fstrSOLID)
        endif

        !C-- geometrical boundary condition
        call dynamic_mat_ass_bc   (cstep, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, &
          &  fstrPARAM, hecLagMAT, t_curr, stepcnt, conMAT=conMAT)
        call dynamic_mat_ass_bc_vl(cstep, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, &
          &  fstrPARAM, hecLagMAT, t_curr, stepcnt, conMAT=conMAT)
        call dynamic_mat_ass_bc_ac(cstep, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, &
          &  fstrPARAM, hecLagMAT, t_curr, stepcnt, conMAT=conMAT)

        ! ----- check convergence
        res = fstr_get_norm_para_contact(hecMAT,hecLagMAT,conMAT,hecMESH)

        if(iter == 1)then
          res0 = res
        endif

        ! ----- check convergence
        if( .not.fstr_is_contact_active() ) then
          maxDLag = 0.0d0
        elseif( abs(maxDLag) < 1.0d-15) then
          maxDLag = 1.0D0
        endif
        call hecmw_allreduce_R1(hecMESH, maxDlag, HECMW_MAX)

        res = dsqrt(res/res0)
        if( hecMESH%my_rank==0 ) then
            if(hecMESH%my_rank==0) write(*,'(a,i5,a,1pe12.4)')"iter: ",iter,", res: ",res
          write(*,'(a,1e15.7)') ' - MaxDLag =',maxDLag
        endif
        if( res<fstrSOLID%step_ctrl(cstep)%converg .and. maxDLag < converg_dlag ) exit

        !   ----  For Parallel Contact with Multi-Partition Domains
        hecMAT%X = 0.0d0
        call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
        call solve_LINEQ_contact(hecMESH,hecMAT,hecLagMAT,conMAT,istat,1.0D0,fstr_is_contact_active())
        call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)
        ! ----- check matrix solver error
        if( istat /= 0 ) then
          if( hecMESH%my_rank == 0) then
            write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', istep
          end if
          fstrSOLID%NRstat_i(knstDRESN) = 4
          fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
          return
        end if

        ! ----- update external nodal displacement increments
        call hecmw_update_R (hecMESH, hecMAT%X, hecMAT%NP, hecMAT%NDOF)

        ! ----- update the strain, stress, and internal force
        do j=1,hecMESH%n_node*ndof
          fstrSOLID%dunode(j)  = fstrSOLID%dunode(j)+hecMAT%X(j)
        enddo
        call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID, t_curr, &
          &   t_delta,iter, fstrDYNAMIC%strainEnergy )

        if(.not. fstrPARAM%nlgeom) exit

        ! ----- update the Lagrange multipliers
        if( fstr_is_contact_active() ) then
          maxDLag = 0.0d0
          do j=1,hecLagMAT%num_lagrange
            hecLagMAT%lagrange(j) = hecLagMAT%lagrange(j) + hecMAT%X(hecMESH%n_node*ndof+j)
            if(dabs(hecMAT%X(hecMESH%n_node*ndof+j))>maxDLag) maxDLag=dabs(hecMAT%X(hecMESH%n_node*ndof+j))
            !              write(*,*)'Lagrange:', j,hecLagMAT%lagrange(j),hecMAT%X(hecMESH%n_node*ndof+j)
          enddo
        endif
      enddo

      fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
      fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

      ! -----  not convergence
      if( iter>fstrSOLID%step_ctrl(cstep)%max_iter .or. res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) then
        if( hecMESH%my_rank == 0) then
          write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', istep
        end if
        fstrSOLID%NRstat_i(knstCITER) = count_step                              ! logging contact iteration
        fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
        if( iter == fstrSOLID%step_ctrl(cstep)%max_iter ) fstrSOLID%NRstat_i(knstDRESN) = 1
        if( res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) fstrSOLID%NRstat_i(knstDRESN) = 2
        return
      end if

      call fstr_scan_contact_state(cstep, istep, count_step, t_delta, ctAlgo, &
                                   hecMESH, fstrSOLID, infoCTChange, hecMAT%B)

      if( hecMAT%Iarray(99)==4 .and. .not. fstr_is_contact_active() ) then
        write(*,*) ' This type of direct solver is not yet available in such case ! '
        write(*,*) ' Please use intel MKL direct solver !'
        call  hecmw_abort(hecmw_comm_get_comm())
      endif

      is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)
      contact_changed_global=0
      if( fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH) ) then
        exit loopFORcontactAnalysis
      elseif( fstr_is_matrixStructure_changed(infoCTChange) ) then
        call fstr_mat_con_contact( cstep, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, infoCTChange, conMAT, fstr_is_contact_active())
        contact_changed_global=1
      endif
      call hecmw_allreduce_I1(hecMESH,contact_changed_global,HECMW_MAX)
      if (contact_changed_global > 0) then
        call hecmw_mat_clear_b( hecMAT )
        call hecmw_mat_clear_b( conMAT )
        call solve_LINEQ_contact_init(hecMESH,hecMAT,hecLagMAT,is_mat_symmetric)
      endif

      if( count_step > max_iter_contact ) exit loopFORcontactAnalysis

    enddo loopFORcontactAnalysis

    !C for couple analysis
      call fstr_solve_dynamic_nlimplicit_couple_post(hecMESH, hecMAT, fstrSOLID, &
        & fstrPARAM, fstrDYNAMIC, fstrCPL, a1, a2, a3, b1, b2, b3, istep, is_cycle)
      if(is_cycle) cycle
      exit
    enddo

    !C-- new displacement, velocity and acceleration
    fstrDYNAMIC%kineticEnergy = 0.0d0
    do j = 1 ,ndof*nnod
      fstrDYNAMIC%ACC (j,2) = -a1*fstrDYNAMIC%ACC(j,1) - a2*fstrDYNAMIC%VEL(j,1) + &
        a3*fstrSOLID%dunode(j)
      fstrDYNAMIC%VEL (j,2) = -b1*fstrDYNAMIC%ACC(j,1) - b2*fstrDYNAMIC%VEL(j,1) + &
        b3*fstrSOLID%dunode(j)
      fstrDYNAMIC%ACC (j,1) = fstrDYNAMIC%ACC (j,2)
      fstrDYNAMIC%VEL (j,1) = fstrDYNAMIC%VEL (j,2)

      fstrSOLID%unode(j)  = fstrSOLID%unode(j)+fstrSOLID%dunode(j)
      fstrDYNAMIC%DISP(j,2) = fstrSOLID%unode(j)

      fstrDYNAMIC%kineticEnergy = fstrDYNAMIC%kineticEnergy + &
        0.5d0*fstrEIG%mass(j)*fstrDYNAMIC%VEL(j,2)*fstrDYNAMIC%VEL(j,2)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, t_delta )

    deallocate(coord)
  end subroutine fstr_Newton_dynamic_contactSLag

  subroutine fstr_solve_dynamic_nlimplicit_couple_init(fstrPARAM, fstrCPL)
    implicit none
    type(fstr_param)         :: fstrPARAM
    type(fstr_couple)        :: fstrCPL
    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==1 .or. &
          fstrPARAM%fg_couple_type==3 .or. &
          fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
    endif
  end subroutine fstr_solve_dynamic_nlimplicit_couple_init

  subroutine fstr_solve_dynamic_nlimplicit_couple_pre(hecMESH, hecMAT, fstrSOLID, &
    & fstrPARAM, fstrDYNAMIC, fstrCPL, restart_step_num, i)
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_solid)         :: fstrSOLID
    type(fstr_param)         :: fstrPARAM
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_couple)        :: fstrCPL
    integer(kint) :: kkk0, kkk1, j, kk, i, restart_step_num
    real(kreal) :: bsize

    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_first /= 0 ) then
        bsize = dfloat( i ) / dfloat( fstrPARAM%fg_couple_first )
        if( bsize > 1.0 ) bsize = 1.0
        do kkk0 = 1, fstrCPL%coupled_node_n
          kkk1 = 3 * kkk0
          fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
          fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
          fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
        enddo
      endif
      if( fstrPARAM%fg_couple_window > 0 ) then
        j = i - restart_step_num + 1
        kk = fstrDYNAMIC%n_step - restart_step_num + 1
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
  end subroutine fstr_solve_dynamic_nlimplicit_couple_pre

  subroutine fstr_solve_dynamic_nlimplicit_couple_post(hecMESH, hecMAT, fstrSOLID, &
    & fstrPARAM, fstrDYNAMIC, fstrCPL, a1, a2, a3, b1, b2, b3, i, is_cycle)
    implicit none
    type(hecmwST_local_mesh) :: hecMESH ! unused - kept for interface compatibility
    type(hecmwST_matrix)     :: hecMAT ! unused - kept for interface compatibility
    type(fstr_solid)         :: fstrSOLID
    type(fstr_param)         :: fstrPARAM
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_couple)        :: fstrCPL
    integer(kint) :: kkk0, kkk1, j, i, revocap_flag
    real(kreal) :: a1, a2, a3, b1, b2, b3
    logical :: is_cycle

    is_cycle = .false.

    if( fstrPARAM%fg_couple == 1 ) then
      if( fstrPARAM%fg_couple_type>1 ) then
        do j=1, fstrCPL%coupled_node_n
          if( fstrCPL%dof == 3 ) then
            kkk0 = j*3
            kkk1 = fstrCPL%coupled_node(j)*3

            fstrCPL%disp (kkk0-2) = fstrSOLID%unode(kkk1-2) + fstrSOLID%dunode(kkk1-2)
            fstrCPL%disp (kkk0-1) = fstrSOLID%unode(kkk1-1) + fstrSOLID%dunode(kkk1-1)
            fstrCPL%disp (kkk0  ) = fstrSOLID%unode(kkk1  ) + fstrSOLID%dunode(kkk1  )

            fstrCPL%velo (kkk0-2) = -b1*fstrDYNAMIC%ACC(kkk1-2,1) - b2*fstrDYNAMIC%VEL(kkk1-2,1) + &
              b3*fstrSOLID%dunode(kkk1-2)
            fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
              b3*fstrSOLID%dunode(kkk1-1)
            fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
              b3*fstrSOLID%dunode(kkk1)
            fstrCPL%accel(kkk0-2) = -a1*fstrDYNAMIC%ACC(kkk1-2,1) - a2*fstrDYNAMIC%VEL(kkk1-2,1) + &
              a3*fstrSOLID%dunode(kkk1-2)
            fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
              a3*fstrSOLID%dunode(kkk1-1)
            fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
              a3*fstrSOLID%dunode(kkk1)
          else
            kkk0 = j*2
            kkk1 = fstrCPL%coupled_node(j)*2

            fstrCPL%disp (kkk0-1) = fstrSOLID%unode(kkk1-1) + fstrSOLID%dunode(kkk1-1)
            fstrCPL%disp (kkk0  ) = fstrSOLID%unode(kkk1  ) + fstrSOLID%dunode(kkk1  )

            fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
              b3*fstrSOLID%dunode(kkk1-1)
            fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
              b3*fstrSOLID%dunode(kkk1)
            fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
              a3*fstrSOLID%dunode(kkk1-1)
            fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
              a3*fstrSOLID%dunode(kkk1)
          endif
        enddo
        call fstr_rcap_send( fstrCPL )
      endif

      select case ( fstrPARAM%fg_couple_type )
        case (4)
          call fstr_rcap_get( fstrCPL )
        case (5)
          call fstr_get_convergence( revocap_flag )
          if( revocap_flag==0 ) is_cycle = .true.
        case (6)
          call fstr_get_convergence( revocap_flag )
          if( revocap_flag==0 ) then
            call fstr_rcap_get( fstrCPL )
            is_cycle = .true.
          else
            if( i /= fstrDYNAMIC%n_step ) call fstr_rcap_get( fstrCPL )
          endif
      end select
    endif
  end subroutine fstr_solve_dynamic_nlimplicit_couple_post

end module fstr_dynamic_nlimplicit

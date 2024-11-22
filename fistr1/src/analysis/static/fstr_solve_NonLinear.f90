!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions on nonlinear analysis

module m_fstr_NonLinearMethod

  use m_fstr
  use m_static_lib
  use m_static_output

  use m_fstr_spring
  use m_fstr_StiffMatrix
  use m_fstr_Update
  use m_fstr_ass_load
  use m_fstr_AddBC
  use m_fstr_Residual
  use m_fstr_Restart

  implicit none
  ! parameters for line search
  real(kind=kreal), parameter :: C_line_search=2.0, delta=0.2, sigma=0.9, eps_wolfe=1.0d-6
  real(kind=kreal), parameter :: omega_wolfe=0.001, Delta_approx_wolfe=0.7
  real(kind=kreal) :: C_wolfe, Q_Wolfe
  integer(kind=kint), parameter :: maxiter_ls = 10

contains


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method
  subroutine fstr_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, &
      restrt_step_num, sub_step, ctime, dtime )

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange

    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: stepcnt
    integer(kind=kint) :: restrt_step_num
    real(kind=kreal)   :: tt0, tt, res, qnrm, rres, tincr, xnrm, dunrm, rxnrm, pot(3), pot0(3)
    real(kind=kreal), allocatable :: coord(:), P(:)
    logical :: isLinear = .false.
    integer(kind=kint) :: iterStatus

    call fstr_Quasi_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, restrt_step_num, sub_step, ctime, dtime)
    return


    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    if(.not. fstrPR%nlgeom)then
      isLinear = .true.
    endif

    hecMAT%NDOF = hecMESH%n_dof
    NDOF = hecMAT%NDOF

    allocate(P(hecMESH%n_node*NDOF))
    allocate(coord(hecMESH%n_node*ndof))

    tincr = dtime
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0

    P = 0.0d0
    stepcnt = 0
    call fstr_init_Newton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, hecLagMAT, pot, pot0, ndof)

    ! ----- Inner Iteration, lagrange multiplier constant
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt = stepcnt+1

      call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, ctime, tincr )
      call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

      ! ----- Set Boundary condition
      call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, stepcnt)
      call hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
      call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)

      !----- SOLVE [Kt]{du}={R}
      if( sub_step == restrt_step_num .and. iter == 1 ) hecMATmpc%Iarray(98) = 1
      if( iter == 1 ) then
        hecMATmpc%Iarray(97) = 2   !Force numerical factorization
      else
        hecMATmpc%Iarray(97) = 1   !Need numerical factorization
      endif
      hecMATmpc%X = 0.0d0
      call fstr_set_current_config_to_mesh(hecMESHmpc,fstrSOLID,coord)
      call solve_LINEQ(hecMESHmpc,hecMATmpc)
      call fstr_recover_initial_config_to_mesh(hecMESHmpc,fstrSOLID,coord)
      call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)

      ! ----- update the small displacement and the displacement for 1step
      !       \delta u^k => solver's solution
      !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i = 1, hecMESH%n_node*ndof
        fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

      call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)

      ! do i = 1, hecMESH%n_node*ndof
      !   write(6,*) fstrSOLID%dunode(i), fstrSOLID%QFORCE(i), fstrSOLID%GL(i)
      ! enddo

      if( isLinear ) exit

      ! ----- check convergence
      iterStatus = fstr_check_iteration_converged(hecMESH, hecMAT, fstrSOLID, ndof, iter, sub_step, cstep, pot)
      if (iterStatus == kitrConverged) exit
      if (iterStatus == kitrDiverged .or. iterStatus==kitrFloatingError) return
    enddo
    ! ----- end of inner loop

    fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
    fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

    ! ----- update the total displacement
    ! u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    fstrSOLID%CutBack_stat = 0
    deallocate(coord)
    deallocate(P)
    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
  end subroutine fstr_Newton


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method combined with Nested iteration of augmentation calculation as suggested
  !> by Simo & Laursen (Compu & Struct, Vol42, pp97-116, 1992 )
  subroutine fstr_Newton_contactALag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM,                   &
      restart_step_num, restart_substep_num, sub_step, ctime, dtime, infoCTChange, conMAT )
    use mContact
    use m_addContactStiffness
    use m_solve_LINEQ_contact

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (fstr_info_contactChange)        :: infoCTChange  !< fstr_info_contactChange
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    type (hecmwST_matrix)                 :: conMAT

    integer(kind=kint) :: ndof
    integer(kind=kint) :: ctAlgo
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: al_step, n_al_step, stepcnt
    real(kind=kreal)   :: tt0, tt, res, res0, res1, maxv, relres, tincr
    integer(kind=kint) :: restart_step_num, restart_substep_num
    logical            :: convg, ctchange
    integer(kind=kint) :: n_node_global
    integer(kind=kint) :: contact_changed_global
    real(kind=kreal), allocatable :: coord(:)
    integer(kind=kint)  :: istat


    ! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

    ctAlgo = fstrPARAM%contact_algo

    hecMAT%NDOF = hecMESH%n_dof
    ndof = hecMAT%NDOF

    fstrSOLID%NRstat_i(:) = 0 ! logging newton iteration(init)

    allocate(coord(hecMESH%n_node*ndof))

    tincr = dtime
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.0d0

    fstrSOLID%dunode(:) = 0.0d0

    if( cstep == 1 .and. sub_step == restart_substep_num ) then
      call fstr_save_originalMatrixStructure(hecMAT)
      if(hecMESH%my_rank==0) write(*,*) "---Scanning initial contact state---"
      call fstr_scan_contact_state( cstep, sub_step, 0, dtime, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )
      call hecmw_mat_copy_profile( hecMAT, conMAT )
      if ( fstr_is_contact_active() ) then
        call fstr_mat_con_contact(cstep, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, infoCTChange, conMAT, fstr_is_contact_active())
      elseif( hecMAT%Iarray(99)==4 ) then
        write(*, *) ' This type of direct solver is not yet available in such case ! '
        write(*, *) ' Please change the solver type to intel MKL direct solver !'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      call solve_LINEQ_contact_init(hecMESH, hecMAT, hecLagMAT, .true.)
    endif

    hecMAT%X = 0.0d0

    stepcnt = 0

    call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    call hecmw_mat_clear_b(conMAT)

    if( fstr_is_contact_active() ) call fstr_ass_load_contactAlag( hecMESH, fstrSOLID, conMAT%B )

    ! ----- Augmentation loop. In case of no contact, it is inactive
    n_al_step = fstrSOLID%step_ctrl(cstep)%max_contiter
    if( .not. fstr_is_contact_active() ) n_al_step = 1

    do al_step = 1, n_al_step

      if( hecMESH%my_rank == 0) then
        if( n_al_step > 1 ) then
          print *, "Augmentation step:", al_step
          write(IMSG, *) "Augmentation step:", al_step
        endif
      end if

      ! ----- Inner Iteration, lagrange multiplier constant
      res0   = 0.0d0
      res1   = 0.0d0
      relres = 1.0d0

      do iter = 1,fstrSOLID%step_ctrl(cstep)%max_iter
        stepcnt = stepcnt+1

        call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, ctime, tincr )
        call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        call hecmw_mat_clear( conMAT )
        conMAT%X = 0.0d0

        ! ----- Contact
        if( al_step == 1 .and. stepcnt == 1 ) then
          maxv = hecmw_mat_diag_max( hecMAT, hecMESH )
          call fstr_set_contact_penalty( maxv )
        endif
        if( fstr_is_contact_active() )  then
          call fstr_contactBC( cstep, iter, hecMESH, conMAT, fstrSOLID )
        endif

        ! ----- Set Boundary condition
        call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, stepcnt, conMAT)

        !----- SOLVE [Kt]{du}={R}
        ! ----  For Parallel Contact with Multi-Partition Domains
        hecMAT%X = 0.0d0
        call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
        call solve_LINEQ_contact(hecMESH, hecMAT, hecLagMAT, conMAT, istat, 1.0D0, fstr_is_contact_active())
        call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)

        call hecmw_update_R (hecMESH, hecMAT%X, hecMAT%NP, hecMESH%n_dof)

        ! ----- update the small displacement and the displacement for 1step
        !       \delta u^k => solver's solution
        !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
        do i = 1, hecMESH%n_node*ndof
          fstrSOLID%dunode(i) = fstrSOLID%dunode(i)+hecMAT%X(i)
        enddo

        ! ----- update the strain, stress, and internal force
        call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter)

        ! ----- Set residual
        if( fstrSOLID%DLOAD_follow /= 0 .or. fstrSOLID%CLOAD_ngrp_rot /= 0 ) &
          call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID, conMAT)

        if( fstr_is_contact_active() ) then
          call hecmw_mat_clear_b( conMAT )
          call fstr_update_contact0(hecMESH, fstrSOLID, conMAT%B)
        endif
        !    Consider SPC condition
        call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, hecMAT%B)
        call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, conMAT%B)
    
        !res = fstr_get_residual(hecMAT%B, hecMESH)
        res = fstr_get_norm_para_contact(hecMAT,hecLagMAT,conMAT,hecMESH)
        ! ----- Gather global residual
        res = sqrt(res)/n_node_global
        if( iter == 1 ) res0 = res
        if( res0 == 0.0d0 ) then
          res0 = 1.0d0
        else
          relres = dabs( res1-res )/res0
        endif

        if( hecMESH%my_rank == 0 ) then
          write(*, '(a,i3,a,2e15.7)') ' - Residual(',iter,') =', res, relres
        endif

        ! ----- check convergence
        if( res < fstrSOLID%step_ctrl(cstep)%converg  .or.     &
          relres < fstrSOLID%step_ctrl(cstep)%converg_ddisp ) exit
        res1 = res

        ! ----- check divergence and NaN
        if( iter == fstrSOLID%step_ctrl(cstep)%max_iter .or. res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) then
          if( hecMESH%my_rank == 0) then
            write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          end if
          fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
          fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sumofiter)
          fstrSOLID%NRstat_i(knstCITER) = al_step                                 ! logging contact iteration
          fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
          if( iter == fstrSOLID%step_ctrl(cstep)%max_iter ) fstrSOLID%NRstat_i(knstDRESN) = 1
          if( res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) fstrSOLID%NRstat_i(knstDRESN) = 2
          return
        end if

      enddo
      ! ----- end of inner loop

      fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
      fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

      ! ----- deal with contact boundary
      call fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchange )
      call fstr_scan_contact_state( cstep, sub_step, al_step, dtime, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )

      contact_changed_global = 0
      if( fstr_is_matrixStructure_changed(infoCTChange) ) then
        call fstr_mat_con_contact( cstep, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, infoCTChange, conMAT, fstr_is_contact_active())
        contact_changed_global = 1
      endif
      call hecmw_allreduce_I1(hecMESH, contact_changed_global, HECMW_MAX)
      if (contact_changed_global > 0) then
        call hecmw_mat_clear_b( hecMAT )
        call hecmw_mat_clear_b( conMAT )
        call solve_LINEQ_contact_init(hecMESH, hecMAT, hecLagMAT, .true.)
      endif

      if( fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH) .and. .not. ctchange ) exit

      ! ----- check divergence
      if( al_step >= fstrSOLID%step_ctrl(cstep)%max_contiter ) then
        if( hecMESH%my_rank == 0) then
          write(   *,'(a,i5,a,i5)') '     ### Contact failed to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        end if
        fstrSOLID%NRstat_i(knstCITER) = al_step                              ! logging contact iteration
        fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
        fstrSOLID%NRstat_i(knstDRESN) = 3
        return
      end if

      ! ----- Set residual for next newton iteration
      call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT )

      if( fstr_is_contact_active() )  then
        call hecmw_mat_clear_b( conMAT )
        call fstr_update_contact0(hecMESH, fstrSOLID, conMAT%B)
      endif

    enddo
    ! ----- end of augmentation loop

    ! ----- update the total displacement
    ! u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i)+fstrSOLID%dunode(i)
    enddo

    fstrSOLID%NRstat_i(knstCITER) = al_step-1 ! logging contact iteration

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    deallocate(coord)
    fstrSOLID%CutBack_stat = 0
  end subroutine fstr_Newton_contactALag


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson method.
  !> Standard Lagrange multiplier algorithm for contact analysis is incoluded in this subroutine.
  subroutine fstr_Newton_contactSLag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT,                  &
      restart_step_num, restart_substep_num, sub_step, ctime, dtime, infoCTChange, conMAT )

    use mContact
    use m_addContactStiffness
    use m_solve_LINEQ_contact

    integer, intent(in)                    :: cstep        !< current loading step
    type (hecmwST_local_mesh)              :: hecMESH      !< hecmw mesh
    type (hecmwST_matrix)                  :: hecMAT       !< hecmw matrix
    type (fstr_solid)                      :: fstrSOLID    !< fstr_solid
    integer, intent(in)                    :: sub_step     !< substep number of current loading step
    real(kind=kreal), intent(in)           :: ctime     !< current time
    real(kind=kreal), intent(in)           :: dtime     !< time increment
    type (fstr_param)                      :: fstrPARAM    !< type fstr_param
    type (fstr_info_contactChange)         :: infoCTChange !< fstr_info_contactChange
    type (hecmwST_matrix_lagrange)         :: hecLagMAT      !< type hecmwST_matrix_lagrange
    type (hecmwST_matrix)                  :: conMAT

    integer(kind=kint) :: ndof
    integer(kind=kint) :: ctAlgo
    integer(kind=kint) :: i, iter, max_iter_contact
    integer(kind=kint) :: stepcnt, count_step
    real(kind=kreal)   :: tt0, tt, res, res0, res1, relres, tincr, resX
    integer(kind=kint) :: restart_step_num, restart_substep_num
    logical            :: is_mat_symmetric
    integer(kind=kint) :: n_node_global
    integer(kind=kint) :: contact_changed_global
    integer(kint)      :: nndof
    real(kreal)        :: q_residual,x_residual
    real(kind=kreal), allocatable :: coord(:)
    integer(kind=kint)  :: istat


    ! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

    if( hecMAT%Iarray(99) == 4 .and. .not. fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH) ) then
      write(*, *) ' This type of direct solver is not yet available in such case ! '
      write(*, *) ' Please use intel MKL direct solver !'
      call  hecmw_abort( hecmw_comm_get_comm() )
    endif

    ctAlgo = fstrPARAM%contact_algo

    hecMAT%NDOF = hecMESH%n_dof
    ndof = hecMAT%NDOF

    fstrSOLID%NRstat_i(:) = 0 ! logging newton iteration(init)

    allocate(coord(hecMESH%n_node*ndof))

    tincr = dtime
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.0d0

    fstrSOLID%dunode(:)  = 0.0d0

    if( cstep==1 .and. sub_step==restart_substep_num  ) then
      call fstr_save_originalMatrixStructure(hecMAT)
      call fstr_scan_contact_state( cstep, sub_step, 0, dtime, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )
      call hecmw_mat_copy_profile( hecMAT, conMAT )
      if ( fstr_is_contact_active() ) then
        call fstr_mat_con_contact(cstep, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, infoCTChange, conMAT, fstr_is_contact_active())
      elseif( hecMAT%Iarray(99)==4 ) then
        write(*, *) ' This type of direct solver is not yet available in such case ! '
        write(*, *) ' Please change the solver type to intel MKL direct solver !'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH)
      call solve_LINEQ_contact_init(hecMESH, hecMAT, hecLagMAT, is_mat_symmetric)
    endif

    stepcnt = 0

    call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    call hecmw_mat_clear_b(conMAT)

    if( fstr_is_contact_active() )  then
      call fstr_ass_load_contact(cstep, hecMESH, conMAT, fstrSOLID, hecLagMAT)
    endif

    fstrSOLID%dunode(:) = 0.0d0

    count_step = 0

    loopFORcontactAnalysis: do while( .TRUE. )
      count_step = count_step+1

      ! ----- Inner Iteration
      res0   = 0.d0
      res1   = 0.d0
      relres = 1.d0

      do iter = 1, fstrSOLID%step_ctrl(cstep)%max_iter
        call hecmw_BARRIER(hecMESH)
        if( myrank == 0 ) print *,'-------------------------------------------------'
        call hecmw_BARRIER(hecMESH)
        stepcnt = stepcnt+1

        call fstr_StiffMatrix(hecMESH, hecMAT, fstrSOLID, ctime, tincr)
        call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        call hecmw_mat_clear( conMAT )
        conMAT%X = 0.0d0

        if( fstr_is_contact_active() ) then
          call fstr_AddContactStiffness(cstep,iter,conMAT,hecLagMAT,fstrSOLID)
        endif

        ! ----- Set Boundary condition
        call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, stepcnt, conMAT)

        nndof = hecMAT%N*hecMAT%ndof

        !----- SOLVE [Kt]{du}={R}
        ! ----  For Parallel Contact with Multi-Partition Domains
        hecMAT%X = 0.0d0
        call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
        q_residual = fstr_get_norm_para_contact(hecMAT,hecLagMAT,conMAT,hecMESH)
        call solve_LINEQ_contact(hecMESH, hecMAT, hecLagMAT, conMAT, istat, 1.0D0, fstr_is_contact_active())
        call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)
        ! ----- check matrix solver error
        if( istat /= 0 ) then
          if( hecMESH%my_rank == 0) then
            write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          end if
          fstrSOLID%NRstat_i(knstDRESN) = 4
          fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
          return
        end if

        x_residual = fstr_get_x_norm_contact(hecMAT,hecLagMAT,hecMESH)

        call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%X,hecMAT%X,resX)
        resX = sqrt(resX)/n_node_global

        if( hecMESH%my_rank==0 ) then
          write(*,'(a,i3,a,e15.7)') ' - ResidualX    (',iter,') =',resX
          write(*,'(a,i3,a,e15.7)') ' - ResidualX+LAG(',iter,') =',sqrt(x_residual)/n_node_global
          write(*,'(a,i3,a,e15.7)') ' - ResidualQ    (',iter,') =',sqrt(q_residual)/n_node_global
        endif

        ! ----- update the small displacement and the displacement for 1step
        do i = 1, hecMESH%n_node*ndof
          fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
        enddo

        ! ----- update the Lagrange multipliers
        if( fstr_is_contact_active() ) then
          do i = 1, hecLagMAT%num_lagrange
            hecLagMAT%lagrange(i) = hecLagMAT%lagrange(i)+hecMAT%X(hecMESH%n_node*ndof+i)
          enddo
        endif

        ! ----- update the strain, stress, and internal force (only QFORCE)
        call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter)

        ! ----- Set residual
        if( fstrSOLID%DLOAD_follow /= 0 .or. fstrSOLID%CLOAD_ngrp_rot /= 0 ) &
          call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT )

        if( fstr_is_contact_active() )  then
          call hecmw_mat_clear_b( conMAT )
          call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
        endif

        res = fstr_get_norm_para_contact(hecMAT,hecLagMAT,conMAT,hecMESH)

        res = sqrt(res)/n_node_global
        if( iter == 1 ) res0 = res
        if( res0 == 0.0d0 ) then
          res0 =1.0d0
        else
          relres = dabs( res1-res )/res0
        endif
        if( hecMESH%my_rank == 0 ) then
          write(*, '(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
        endif

        ! ----- check convergence
        if( res < fstrSOLID%step_ctrl(cstep)%converg  .or.     &
          relres < fstrSOLID%step_ctrl(cstep)%converg_ddisp ) exit
        res1 = res

        ! ----- check divergence and NaN
        if( iter == fstrSOLID%step_ctrl(cstep)%max_iter .or. res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) then
          if( hecMESH%my_rank == 0) then
            write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          end if
          fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
          fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sumofiter)
          fstrSOLID%NRstat_i(knstCITER) = count_step                              ! logging contact iteration
          fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
          if( iter == fstrSOLID%step_ctrl(cstep)%max_iter ) fstrSOLID%NRstat_i(knstDRESN) = 1
          if( res > fstrSOLID%step_ctrl(cstep)%maxres .or. res /= res ) fstrSOLID%NRstat_i(knstDRESN) = 2
          return
        end if

      enddo
      ! ----- end of inner loop

      fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
      fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

      call fstr_scan_contact_state( cstep, sub_step, count_step, dtime, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )

      if( hecMAT%Iarray(99) == 4 .and. .not. fstr_is_contact_active() ) then
        write(*, *) ' This type of direct solver is not yet available in such case ! '
        write(*, *) ' Please use intel MKL direct solver !'
        call  hecmw_abort( hecmw_comm_get_comm() )
      endif

      is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH)
      contact_changed_global = 0
      if( fstr_is_matrixStructure_changed(infoCTChange) ) then
        call fstr_mat_con_contact( cstep, ctAlgo, hecMAT, fstrSOLID, hecLagMAT, infoCTChange, conMAT, fstr_is_contact_active())
        contact_changed_global = 1
      endif

      if( fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH) ) exit loopFORcontactAnalysis

      call hecmw_allreduce_I1(hecMESH, contact_changed_global, HECMW_MAX)
      if (contact_changed_global > 0) then
        call hecmw_mat_clear_b( hecMAT )
        call hecmw_mat_clear_b( conMAT )
        call solve_LINEQ_contact_init(hecMESH, hecMAT, hecLagMAT, is_mat_symmetric)
      endif

      ! ----- check divergence
      if( count_step >= fstrSOLID%step_ctrl(cstep)%max_contiter ) then
        if( hecMESH%my_rank == 0) then
          write(   *,'(a,i5,a,i5)') '     ### Contact failed to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        end if
        fstrSOLID%NRstat_i(knstCITER) = count_step                              ! logging contact iteration
        fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
        fstrSOLID%NRstat_i(knstDRESN) = 3
        return
      end if

      ! ----- Set residual for next newton iteration
      call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT )

      if( fstr_is_contact_active() )  then
        call hecmw_mat_clear_b( conMAT )
        call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
      endif

    enddo loopFORcontactAnalysis

    fstrSOLID%NRstat_i(knstCITER) = count_step ! logging contact iteration

    ! ----- update the total displacement
    !       u_{n+1} = u_{n} + \Delta u_{n+1}
    do i = 1, hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i)+fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState(hecMESH, fstrSOLID, tincr)
    call fstr_update_contact_TangentForce( fstrSOLID )
    if( fstrSOLID%n_embeds > 0 .and. paraContactFlag ) then
      call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%EMBED_NFORCE,1)
      call fstr_Update_NDForce_SPC( cstep, hecMESH, fstrSOLID, hecMAT%B )
    endif

    deallocate(coord)
    fstrSOLID%CutBack_stat = 0
  end subroutine fstr_Newton_contactSLag

  subroutine fstr_init_Newton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, hecLagMAT, pot, pot0, ndof)
    use m_fstr_Update
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    real(kind=kreal), intent(inout) :: pot(3), pot0(3)
    integer(kind=kint), intent(in) :: ndof

    real(kind=kreal)   :: res

    fstrSOLID%dunode(:) = 0.0d0
    fstrSOLID%NRstat_i(:) = 0 ! logging newton iteration(init)

    call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    !calc initial potential
    !! initialize du for non-zero Dirichlet condition
    call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 1, RHSvector=fstrSOLID%dunode)
    !! update residual vector
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)

    call hecmw_InnerProduct_R(hecMESH, ndof, hecMAT%B, hecMAT%B, res)
    res = sqrt(res)
    pot(1) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,1)
    pot(2) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,2)
    pot(3) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,3)
    pot0(1:3) = pot(1:3)
    if( hecMESH%my_rank == 0 ) then
      write(IMSG,'(A4,7(",",A14))') "iter","residual","p1","p2","p3","diffp1","diffp2","diffp3"
      write(IMSG,'(I4,7(",",1pE14.7))') 0,res,pot(1:3),pot(1:3)-pot0(1:3)
    endif

    !! reset du and stress and strain 
    fstrSOLID%dunode(:) = 0.0d0
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    !end calc initial potential
  end subroutine fstr_init_Newton


  !> \breaf This subroutine calculate residual vector
  subroutine fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    use m_fstr_Update
    implicit none
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter

    ! ----- update the strain, stress, and internal force
    call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter)

    ! ----- Set residual
    if( fstrSOLID%DLOAD_follow /= 0 .or. fstrSOLID%CLOAD_ngrp_rot /= 0 ) &
      & call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM )

    call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID)
  end subroutine fstr_calc_residual_vector

  !> \breaf This subroutine calculate residual vector
  subroutine fstr_calc_residual_vector_with_X(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    use m_fstr_Update
    implicit none
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter

    integer :: i
    real(kind=kreal) :: dunode_bak(hecMAT%ndof*hecMESH%n_node)

    do i=1, hecMAT%ndof*hecMESH%n_node
      dunode_bak(i) = fstrSOLID%dunode(i)
      fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
    enddo
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    do i=1, hecMAT%ndof*hecMESH%n_node
      fstrSOLID%dunode(i) = dunode_bak(i)
    enddo
  end subroutine fstr_calc_residual_vector_with_X


  !> \breaf This function check iteration status
  function fstr_check_iteration_converged(hecMESH, hecMAT, fstrSOLID, ndof, iter, sub_step, cstep, pot) result(iterStatus)
    implicit none
    integer(kind=kint) :: iterStatus

    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer(kind=kint), intent(in) :: ndof
    integer(kind=kint), intent(in) :: iter
    integer(kind=kint), intent(in) :: sub_step, cstep
    real(kind=kreal), intent(inout) :: pot(3)

    real(kind=kreal)   :: res, qnrm, rres, xnrm, dunrm, rxnrm, pot0(3)

    iterStatus = kitrContinue

    call hecmw_InnerProduct_R(hecMESH, ndof, hecMAT%B, hecMAT%B, res)
    res = sqrt(res)
    call hecmw_InnerProduct_R(hecMESH, ndof, hecMAT%X, hecMAT%X, xnrm)
    xnrm = sqrt(xnrm)
    call hecmw_innerProduct_R(hecMESH, ndof, fstrSOLID%QFORCE, fstrSOLID%QFORCE, qnrm)
    qnrm = sqrt(qnrm)
    if (qnrm < 1.0d-8) qnrm = 1.0d0
    if( iter == 1 ) then
      dunrm = xnrm
    else
      call hecmw_InnerProduct_R(hecMESH, ndof, fstrSOLID%dunode, fstrSOLID%dunode, dunrm)
      dunrm = sqrt(dunrm)
    endif
    rres = res/qnrm
    rxnrm = xnrm/dunrm

    pot0(1:3) = pot(1:3)
    pot(1) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,1)
    pot(2) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,2)
    pot(3) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,3)

    if( hecMESH%my_rank == 0 ) then
      if (qnrm == 1.0d0) then
        write(*,"(a,i8,a,1pe11.4,a,1pe11.4)")" iter:", iter, ", residual(abs):", rres, ", disp.corr.:", rxnrm
      else
        write(*,"(a,i8,a,1pe11.4,a,1pe11.4)")" iter:", iter, ", residual:", rres, ", disp.corr.:", rxnrm
      endif
      write(IMSG,'(I4,7(",",1pE14.7))') iter,res,pot(1:3),pot(1:3)-pot0(1:3)
    endif
    if( hecmw_mat_get_flag_diverged(hecMAT) == kNO ) then
      if( rres < fstrSOLID%step_ctrl(cstep)%converg .or. &
          rxnrm < fstrSOLID%step_ctrl(cstep)%converg_ddisp ) then
          iterStatus=kitrConverged
          return
      endif
    endif

    ! ----- check divergence and NaN
    if ( iter == fstrSOLID%step_ctrl(cstep)%max_iter .or. rres > fstrSOLID%step_ctrl(cstep)%maxres) &
      iterStatus = kitrDiverged
    if (rres /= rres ) iterStatus = kitrFloatingError
    if (iterStatus /= kitrContinue) then
      if( hecMESH%my_rank == 0) then
        write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
      end if
      fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
      fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sumofiter)
      fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
      if( iter == fstrSOLID%step_ctrl(cstep)%max_iter ) fstrSOLID%NRstat_i(knstDRESN) = 1
      if( rres > fstrSOLID%step_ctrl(cstep)%maxres .or. rres /= rres ) fstrSOLID%NRstat_i(knstDRESN) = 2
      return
    end if
  end function fstr_check_iteration_converged

  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method
  subroutine fstr_Quasi_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, &
    restrt_step_num, sub_step, ctime, dtime )
    implicit none

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange

    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: stepcnt
    integer(kind=kint) :: restrt_step_num
    real(kind=kreal)   :: tt0, tt, res, qnrm, rres, tincr, xnrm, dunrm, rxnrm, pot(3), pot0(3)
    real(kind=kreal), allocatable :: coord(:), P(:)
    logical :: isLinear = .false.
    integer(kind=kint) :: iterStatus

    real(kind=kreal), allocatable :: z_k(:), s_k(:,:), y_k(:,:), g_prev(:), rho_k(:)
    real(kind=kreal) :: sdoty
    integer, parameter :: n_mem=10
    integer :: len_vector
    integer(kind=kint) :: k

    integer :: u_debug
    integer :: max_iter_bak
    logical :: flag_approx_Wolfe

    max_iter_bak = fstrSOLID%step_ctrl(cstep)%max_iter
    fstrSOLID%step_ctrl(cstep)%max_iter = 100*fstrSOLID%step_ctrl(cstep)%max_iter

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    if(.not. fstrPR%nlgeom)then
      isLinear = .true.
    endif

    hecMAT%NDOF = hecMESH%n_dof
    NDOF = hecMAT%NDOF

    allocate(P(hecMESH%n_node*NDOF))
    allocate(coord(hecMESH%n_node*ndof))
    P = 0.0d0
    stepcnt = 0

    tincr = dtime
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0
    call fstr_init_Newton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, hecLagMAT, pot, pot0, ndof)

    !! initialize du for non-zero Dirichlet condition
    call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 1, RHSvector=fstrSOLID%dunode)
    !! update residual vector
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)

    len_vector = hecMESH%n_node*ndof
    allocate(z_k(len_vector))
    allocate(s_k(len_vector, n_mem))
    allocate(y_k(len_vector, n_mem))
    allocate(g_prev(len_vector))
    allocate(rho_k(n_mem))
    z_k(:) = 0.0d0
    s_k(:,:) = 0.0d0
    y_k(:,:) = 0.0d0
    g_prev(:) = 0.0d0
    rho_k(:) = 0.0d0
    do i=1,hecMESH%n_node*ndof
      y_k(i,1) = -hecMAT%B(i)
    enddo

    ! parameter to judge Wolfe/approx Wolfe selection
    C_wolfe = 0.0d0
    Q_wolfe = 0.0d0
    flag_approx_Wolfe = .false.
    ! ----- Inner Iteration, lagrange multiplier constant
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt = stepcnt+1

      ! ----- calculate search direction by limited BFGS method
      do i=1,hecMESH%n_node*ndof
          g_prev(i) = -hecMAT%B(i)
      enddo
      call fstr_calc_direction_LBFGS(hecMesh, g_prev, s_k, y_k, rho_k, z_k, n_mem)

      ! ! ----- Set Boundary condition
      call fstr_AddBC_to_direction_vector(z_k, hecMESH,fstrSOLID, cstep)
      

      !----- line search of step length
      call fstr_line_search_along_direction(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k)

      ! ----- update the small displacement and the displacement for 1step
      !       \delta u^k => solver's solution
      !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i = 1, hecMESH%n_node*ndof
        fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

      !! set du for non-zero Dirichlet condition
      ! call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 1, RHSvector=fstrSOLID%dunode)
      call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
 
      if( isLinear ) exit

      ! ----- check convergence
      iterStatus = fstr_check_iteration_converged(hecMESH, hecMAT, fstrSOLID, ndof, iter, sub_step, cstep, pot)
      if (iterStatus == kitrConverged) exit
      if (iterStatus == kitrDiverged .or. iterStatus==kitrFloatingError) return
      ! if (iterStatus == kitrDiverged) exit
      ! if (iterStatus==kitrFloatingError) return

      do k=n_mem, 2, -1
        s_k(:,k) = s_k(:,k-1)
        y_k(:,k) = y_k(:,k-1)
        rho_k(k) = rho_k(k-1)
      enddo

      do i = 1, hecMESH%n_node*ndof
        s_k(i,1) = hecMAT%X(i)
        y_k(i,1) = -hecMAT%B(i) - g_prev(i)
      enddo
      call hecmw_innerProduct_R(hecMESH,ndof,s_k(:,1), y_k(:,1), sdoty)
      if (abs(sdoty) < 1.0d-10) then
        rho_k(1) = 0.0d0
      else
        rho_k(1) = 1.0d0/sdoty
      endif
    enddo
    ! ----- end of inner loop

    fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
    fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

    ! ----- update the total displacement
    ! u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    fstrSOLID%CutBack_stat = 0
    deallocate(coord)
    deallocate(P)
    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    fstrSOLID%step_ctrl(cstep)%max_iter = max_iter_bak
  end subroutine fstr_Quasi_Newton

  subroutine fstr_calc_direction_LBFGS(hecMESH, g_prev, s_k, y_k, rho_k, z_k, n_mem)
    implicit none

    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    real(kind=kreal) :: z_k(:), s_k(:,:), y_k(:,:), g_prev(:), rho_k(:)
    integer :: n_mem

    real(kind=kreal), allocatable :: q(:)
    real(kind=kreal) :: alpha(n_mem), beta
    real(kind=kreal) :: sdotq, ysq, gamma, ydotz
    
    integer :: len_vector, ndof
    integer(kind=kint) :: k,i

    ndof = hecMesh%n_dof
    len_vector = hecMESH%n_node*hecMesh%n_dof
    allocate(q(len_vector))

    q(1:len_vector) = g_prev(1:len_vector)
    ! do i=1, len_vector
    !   write(6,*) 'q ', q(i)
    ! enddo

    do k=1, n_mem
      call hecmw_innerProduct_R(hecMESH,ndof,s_k(:,k), q, sdotq)
      alpha(k) = rho_k(k) * sdotq
      do i=1, len_vector
        q(i) = q(i) - alpha(k)*y_k(i,k)
      enddo
    enddo
    call hecmw_innerProduct_R(hecMESH,ndof,y_k(:,1), y_k(:,1), ysq)
    if (abs(rho_k(1)) < 1.0d-10) then
      ! if (ysq < 1.0d-10) then
        gamma = 1.0d0
      ! else
      !   gamma = 1.0d0/ysq
      ! endif
    else
      gamma = 1.0d0/(rho_k(1)*ysq)
    endif

    do i=1, len_vector
      z_k(i) = gamma*q(i)
    enddo

    do k=n_mem, 1, -1
      call hecmw_innerProduct_R(hecMESH,ndof,y_k(:,k), z_k, ydotz)
      beta = rho_k(k)*ydotz
      do i=1, len_vector
        z_k(i)=z_k(i) + s_k(i,k)*(alpha(k)-beta)
      enddo
    enddo

    ! do i=1, len_vector
    !   z_k(i) = -z_k(i)
    ! enddo
    deallocate(q)
  end subroutine fstr_calc_direction_LBFGS

  subroutine fstr_AddBC_to_direction_vector(z_k, hecMESH,fstrSOLID, cstep)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer(kind=kint) :: cstep
    real(kind=kreal) :: z_k(:)

    integer(kind=kint) :: ig0, grpid, ig, ityp, idofS, idofE, iS0, iE0, ik, in, idof, ndof

    ndof = hecMesh%n_dof
    !   ----- Prescibed displacement Boundary Conditions
    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      idofS = ityp/10
      idofE = ityp - idofS*10
      !
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      !
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        !
        do idof = idofS, idofE
            z_k(ndof*(in-1)+idof) = 0.0d0
        enddo
      enddo
    enddo
  end subroutine fstr_AddBC_to_direction_vector

  subroutine fstr_apply_alpha0(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, h_prime, pot)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal) :: h_prime, pot

    hecMat%X(:) = 0.0d0
    call hecmw_innerProduct_R(hecMESH, hecMAT%NDOF, hecMat%B, z_k, h_prime)
    pot = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,1)
  end subroutine fstr_apply_alpha0

  subroutine fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha, h_prime, pot)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal), intent(in) :: alpha
    real(kind=kreal) :: h_prime, pot

    hecMat%X(:) = -alpha*z_k(:)
    call fstr_calc_residual_vector_with_X(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    call hecmw_innerProduct_R(hecMESH, hecMAT%NDOF, hecMat%B, z_k, h_prime)
    pot = fstr_get_potential_with_X(cstep,hecMESH,hecMAT,fstrSOLID,1)
  end subroutine fstr_apply_alpha

  subroutine fstr_line_search_along_direction(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal) :: z_k(:)

    real(kind=kreal) :: alpha_S, alpha_E
    real(kind=kreal) :: h_prime_0, h_prime_S, h_prime_E, h_prime_c
    real(kind=kreal) :: pot_0, pot_S, pot_E, pot_c
    real(kind=kreal) :: alpha_c
    logical :: flag_converged
    integer :: ndof, len_vector
    real(kind=kreal) :: res

    integer(kind=kint) :: i, ierr, iter_ls
    real(kind=kreal) :: z_max
    integer :: dummy

    ndof = hecMAT%NDOF

    call fstr_apply_alpha0(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, h_prime_0, pot_0)
    if (h_prime_0 > 0.0d0) then
      write(6,*) 'residual vector is not directed to potential decretion.', h_prime_0
      stop
    endif

    alpha_S = 0.0d0
    h_prime_S = h_prime_0
    pot_S = pot_0

    len_vector = hecMESH%n_node*hecMesh%n_dof
    if (iter==1) then
      z_max = 0.0d0
      do i=1, len_vector
        z_max = max(z_max, abs(z_k(i)))
      end do
      call MPI_Allreduce                                              &
      &       (MPI_IN_PLACE, z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX,              &
      &        hecMESH%MPI_COMM, ierr)
      if (z_max==0.0d0) then
        write(6,*) 'direction for line search is zero-vector'
        stop
      endif
    else
     z_max = 1.0d0
    endif
    alpha_E = 1.0d0 / z_max /C_line_search

    h_prime_E = h_prime_0
    do while (h_prime_E < 0.0d0)
      alpha_E = alpha_E * C_line_search
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha_E, h_prime_E, pot_E)
      if( hecMESH%my_rank == 0 ) then
        write(6,'(a, 3es27.16e3)') 'alpha_E, h_prime_E, pot_E', alpha_E, h_prime_E, pot_E
      endif
    enddo

		Q_Wolfe = 1 + Delta_approx_wolfe * Q_Wolfe
		C_Wolfe = C_Wolfe + (abs(pot_0)-C_Wolfe) / Q_Wolfe

    flag_converged = .false.
    iter_ls=1
    do while (iter_ls<maxiter_ls)
      call fstr_set_secant(alpha_S, h_prime_S, alpha_E, h_prime_E, alpha_c)
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha_c, h_prime_c, pot_c)

      ! write(6,'(a, 5es27.16e3)') 'alpha_S, alpha_E, alpha_c, h_prime_c, pot_c', alpha_S, alpha_E, alpha_c, h_prime_c, pot_c
      ! call hecmw_InnerProduct_R(hecMESH, ndof, hecMAT%B, hecMAT%B, res)
      ! res = sqrt(res)
      ! write(6,*) 'res: ', res

      if (abs( pot_c - pot_0) <= (omega_Wolfe * C_Wolfe) ) then
        flag_converged = fstr_approx_wolfe_condition(alpha_S, alpha_E, alpha_c, &
          &  h_prime_0, h_prime_E, h_prime_c, pot_0, pot_E, pot_c)
      else
        flag_converged = fstr_wolfe_condition(alpha_S, alpha_E, alpha_c, &
          &  h_prime_0, h_prime_E, h_prime_c, pot_0, pot_E, pot_c)
      endif

      if (flag_converged) exit

      if (h_prime_c < 0.0d0) then
        alpha_S = alpha_c
        h_prime_S = h_prime_c
        pot_S = pot_c
      else
        alpha_E = alpha_c
        h_prime_E = h_prime_c
        pot_E = pot_c
      endif
      iter_ls = iter_ls +1
    enddo
    if( hecMESH%my_rank == 0 ) then
      write(6,'(a, 5es27.16e3)') 'alpha_S, alpha_E, alpha_c, h_prime_c, pot_c', &
          &  alpha_S, alpha_E, alpha_c, h_prime_c, pot_c
    endif    
  end subroutine fstr_line_search_along_direction

  subroutine fstr_set_secant(a, Fa, b, Fb, c)
    implicit none
    real(kind=kreal), intent(in) :: a,b, Fa,Fb
    real(kind=kreal), intent(out) :: c
    if(Fb /= Fa) then
      c = (a*Fb  - b*Fa) / ( Fb - Fa )
    else
      c = 0.5*a + 0.5*b
    endif
  end subroutine fstr_set_secant

  function fstr_wolfe_condition(alpha_0, alpha_a, alpha_c, h_prime_0, &
    &  h_prime_a, h_prime_c, pot_0, pot_a, pot_c) result(flag_converged)
    implicit none
    logical :: flag_converged

    real(kind=kreal), intent(in) :: alpha_0, alpha_a, alpha_c, h_prime_0, h_prime_a, h_prime_c, pot_0, pot_a, pot_c

    real(kind=kreal) :: wolfe1_left, wolfe1_right, wolfe2_left, wolfe2_right

    wolfe1_left = pot_c - pot_0
    wolfe1_right = delta*h_prime_c*alpha_c

    wolfe2_left = h_prime_c
    wolfe2_right = sigma*h_prime_0

    flag_converged = (wolfe1_left<=wolfe1_right) .and. (wolfe2_left >= wolfe2_right)
    ! flag_converged = (abs(h_prime_c) < abs(sigma*h_prime_0))
    ! flag_converged = h_prime_c>=(sigma*h_prime_0)
  end function fstr_wolfe_condition

  function fstr_approx_wolfe_condition(alpha_0, alpha_a, alpha_c, h_prime_0, h_prime_a, &
    &  h_prime_c, pot_0, pot_a, pot_c) result(flag_converged)
    implicit none
    logical :: flag_converged

    real(kind=kreal), intent(in) :: alpha_0, alpha_a, alpha_c
    real(kind=kreal), intent(in) :: h_prime_0, h_prime_a, h_prime_c, pot_0, pot_a, pot_c

    real(kind=kreal) :: wolfe1_left, wolfe1_right, wolfe2_left, wolfe2_right
    real(kind=kreal) :: pot_Eps

    wolfe1_left =  ( 2.0 * delta - 1.0d0 ) * h_prime_0
    wolfe1_right = h_prime_c

    wolfe2_left = h_prime_c
    wolfe2_right = sigma*h_prime_0

    flag_converged = &
      (wolfe1_left>=wolfe1_right) &
      .and. (wolfe2_left >= wolfe2_right) &
      .and. (pot_c <= pot_0 + (eps_wolfe*C_Wolfe))
  end function fstr_approx_wolfe_condition
end module m_fstr_NonLinearMethod

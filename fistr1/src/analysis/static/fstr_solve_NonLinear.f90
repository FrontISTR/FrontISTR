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
    real(kind=kreal)   :: tt0, tt, res, qnrm, rres, tincr, xnrm, dunrm, rxnrm
    real(kind=kreal), allocatable :: coord(:), P(:)
    logical :: isLinear = .false.

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
    fstrSOLID%dunode(:) = 0.0d0
    fstrSOLID%NRstat_i(:) = 0 ! logging newton iteration(init)

    call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

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

      ! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter)

      ! ----- Set residual
      if( fstrSOLID%DLOAD_follow /= 0 .or. fstrSOLID%CLOAD_ngrp_rot /= 0 ) &
        & call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM )

      call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID)

      if( isLinear ) exit

      ! ----- check convergence
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
      if( hecMESH%my_rank == 0 ) then
        if (qnrm == 1.0d0) then
          write(*,"(a,i8,a,1pe11.4,a,1pe11.4)")" iter:", iter, ", residual(abs):", rres, ", disp.corr.:", rxnrm
        else
          write(*,"(a,i8,a,1pe11.4,a,1pe11.4)")" iter:", iter, ", residual:", rres, ", disp.corr.:", rxnrm
        endif
      endif
      if( hecmw_mat_get_flag_diverged(hecMAT) == kNO ) then
        if( rres < fstrSOLID%step_ctrl(cstep)%converg ) exit
        if( rxnrm < fstrSOLID%step_ctrl(cstep)%converg ) exit
      endif

      ! ----- check divergence and NaN
      if( iter == fstrSOLID%step_ctrl(cstep)%max_iter .or. rres > fstrSOLID%step_ctrl(cstep)%maxres .or. rres /= rres ) then
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
          relres < fstrSOLID%step_ctrl(cstep)%converg ) exit
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
          relres < fstrSOLID%step_ctrl(cstep)%converg ) exit
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


end module m_fstr_NonLinearMethod

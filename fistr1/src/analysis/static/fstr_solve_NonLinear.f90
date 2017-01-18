!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions on nonlinear analysis

module m_fstr_NonLinearMethod

  use m_fstr
  use m_fstr_para_contact
  use m_static_lib
  use m_static_output

  use m_fstr_spring
  use m_fstr_StiffMatrix
  use m_fstr_Update
  use m_fstr_ass_load
  use m_fstr_AddBC
  use m_fstr_Residual
  use m_fstr_Restart
  use fstr_matrix_con_contact

  implicit none

  contains


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method
  subroutine fstr_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, &
    restrt_step_num, sub_step )

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (fstrST_matrix_contact_lagrange) :: fstrMAT   !< type fstrST_matrix_contact_lagrange

    integer(kind=kint) :: ndof
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: stepcnt
    real(kind=kreal)   :: tt0, tt, res, res0, res1, relres, tincr
    integer(kind=kint) :: restrt_step_num
    integer(kind=kint) :: n_node_global
    real(kind=kreal), pointer :: coord(:)

    ! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH, n_node_global, HECMW_SUM)

    hecMAT%NDOF = hecMESH%n_dof
    ndof = hecMAT%NDOF

    allocate(coord(hecMESH%n_node*ndof))

    tincr = fstrSOLID%step_ctrl(cstep)%initdt
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0

    call cpu_time(tt0)

    hecMAT%X = 0.0d0

    stepcnt = 0

    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    fstrSOLID%dunode(:) = 0.0d0

    ! ----- Inner Iteration, lagrange multiplier constant
    res1   = 0.0d0
    relres = 1.0d0

    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt = stepcnt+1

      call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )
      call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

      ! ----- Set Boundary condition
      call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, stepcnt)

      !----- SOLVE [Kt]{du}={R}
      if( sub_step == restrt_step_num .and. iter == 1 ) hecMAT%Iarray(98) = 1
      if( iter == 1 ) then
        hecMAT%Iarray(97) = 2   !Force numerical factorization
      else
        hecMAT%Iarray(97) = 1   !Need numerical factorization
      endif
      hecMAT%X = 0.0d0
      call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
      CALL solve_LINEQ(hecMESH,hecMAT,imsg)
      call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)

      if( hecMESH%n_dof == 3 ) then
        call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
      else if( hecMESH%n_dof == 2 ) then
        call hecmw_update_2_R (hecMESH, hecMAT%X, hecMAT%NP)
      endif

      ! ----- update the small displacement and the displacement for 1step
      !       \delta u^k => solver's solution
      !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i = 1, hecMESH%n_node*ndof
        fstrSOLID%dunode(i) = fstrSOLID%dunode(i)+hecMAT%X(i)
      enddo

      ! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, tincr, iter)

      ! ----- Set residual
      if( fstrSOLID%DLOAD_follow /= 0 )                                  &
        call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )

      call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID)
      res = fstr_get_residual(hecMAT%B, hecMESH)

      ! ----- Gather global residual
      res = dsqrt( res )/n_node_global
      if( iter == 1 ) res0 = res
      if( res0 == 0.0d0 ) then
        res0 = 1.0d0
      else
        relres = dabs( res1-res )/res0
      endif
      if( hecMESH%my_rank == 0 ) then
        write(*, '(a, i3, a, 2e15.7)') ' - Residual(', iter, ') =', res, relres
        write(ista, '(''iter='', i5, ''res/res0='', 2e15.7)') iter, res, relres
      endif

      ! ----- check convergence
      if( res < fstrSOLID%step_ctrl(cstep)%converg  .or.     &
        relres < fstrSOLID%step_ctrl(cstep)%converg ) exit
      res1 = res

    enddo
    ! ----- end of inner loop

    ! -----  not convergence
    if( iter > fstrSOLID%step_ctrl(cstep)%max_iter ) then
      if( hecMESH%my_rank == 0 ) then
        write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
      end if
      stop
    end if

    ! ----- update the total displacement
    ! u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i)+fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    call cpu_time(tt)

    if( hecMESH%my_rank == 0 ) then
      write(ISTA,'("### Converged in NR ietration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
    endif

    deallocate(coord)
  end subroutine fstr_Newton


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method combined with Nested iteration of augmentation calculation as suggested
  !> by Simo & Laursen (Compu & Struct, Vol42, pp97-116, 1992 )
  subroutine fstr_Newton_contactALag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM,                   &
    restart_step_num, restart_substep_num, sub_step, infoCTChange )
  use mContact

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (fstr_info_contactChange)        :: infoCTChange  !< fstr_info_contactChange
    type (fstrST_matrix_contact_lagrange) :: fstrMAT !< type fstrST_matrix_contact_lagrange

    integer(kind=kint) :: ndof
    integer(kind=kint) :: ctAlgo
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: al_step, n_al_step, stepcnt
    real(kind=kreal)   :: tt0, tt, res, res0, res1, maxv, relres, tincr
    integer(kind=kint) :: restart_step_num, restart_substep_num
    logical            :: convg, ctchange
    integer(kind=kint) :: n_node_global
    real(kind=kreal), pointer :: coord(:)

    ! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

    ctAlgo = fstrPARAM%contact_algo

    hecMAT%NDOF = hecMESH%n_dof
    ndof = hecMAT%NDOF

    allocate(coord(hecMESH%n_node*ndof))

    tincr = fstrSOLID%step_ctrl(cstep)%initdt
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.0d0

    call cpu_time(tt0)

    if( cstep == restart_step_num .and. sub_step == restart_substep_num ) then
      if(hecMESH%my_rank==0) write(*,*) "---Scanning initial contact state---"
      call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange )
      if(hecMESH%my_rank==0) write(*,*)
    endif

    hecMAT%X = 0.0d0

    stepcnt = 0

    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    ! ----- Augmentation loop. In case of no contact, it is inactive
    n_al_step = 10
    if( .not. fstr_is_contact_active() ) n_al_step = 1

    do al_step = 1, n_al_step

      if( hecMESH%my_rank == 0) then
        if( n_al_step > 1 ) then
          print *, "Augmentation step:", al_step
          write(IMSG, *) "Augmentation step:", al_step
        endif
      end if

      fstrSOLID%dunode(:) = 0.0d0

      ! ----- Inner Iteration, lagrange multiplier constant
      res1   = 0.0d0
      relres = 1.0d0

      do iter = 1,fstrSOLID%step_ctrl(cstep)%max_iter
        stepcnt = stepcnt+1

        call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )
        call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        ! ----- Contact
        if( fstr_is_contact_active() .and. stepcnt==1 )  then
          maxv = hecmw_mat_diag_max( hecMAT, hecMESH )
          call fstr_set_contact_penalty( maxv )
        endif
        if( fstr_is_contact_active() ) then
          call fstr_contactBC( iter, hecMESH, hecMAT, fstrSOLID )
        endif

        ! ----- Set Boundary condition
        call fstr_AddBC(cstep, hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrMAT,stepcnt)

        !----- SOLVE [Kt]{du}={R}
        if( sub_step == restart_step_num .and. iter == 1 ) hecMAT%Iarray(98) = 1
        if( iter == 1 ) then
          hecMAT%Iarray(97) = 2   !Force numerical factorization
        else
          hecMAT%Iarray(97) = 1   !Need numerical factorization
        endif
        hecMAT%X = 0.0d0
        call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
        CALL solve_LINEQ(hecMESH,hecMAT,imsg)
        call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)

        if( hecMESH%n_dof == 3 ) then
          call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
          if( hecMESH%my_rank == 0 ) then
            write(IMSG, *) 'hecmw_update_3_R: OK'
          endif
        else if( hecMESH%n_dof == 2 ) then
          call hecmw_update_2_R (hecMESH, hecMAT%X, hecMAT%NP)
          if( hecMESH%my_rank == 0 ) then
            write(IMSG, *) 'hecmw_update_2_R: OK'
          endif
        endif

        ! ----- update the small displacement and the displacement for 1step
        !       \delta u^k => solver's solution
        !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
        do i = 1, hecMESH%n_node*ndof
          fstrSOLID%dunode(i) = fstrSOLID%dunode(i)+hecMAT%X(i)
        enddo

        ! ----- update the strain, stress, and internal force
        call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, tincr, iter)

        ! ----- Set residual
        if( fstrSOLID%DLOAD_follow /= 0 )                                 &
          call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID)

        if( fstr_is_contact_active() ) then
          call fstr_update_contact0(hecMESH, fstrSOLID, hecMAT%B)
        endif

        res = fstr_get_residual(hecMAT%B, hecMESH)

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
          write(ISTA, '(''iter='',I5,''res/res0='',2E15.7)') iter, res, relres
        endif

        ! ----- check convergence
        if( res < fstrSOLID%step_ctrl(cstep)%converg  .or.     &
          relres < fstrSOLID%step_ctrl(cstep)%converg ) exit
        res1 = res

      enddo
      ! ----- end of inner loop

      ! -----  not convergence
      if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
        if( hecMESH%my_rank==0) then
          write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        end if
        stop
      end if

      ! ----- deal with contact boundary
      convg = .true.
      ctchange = .false.
      if( associated(fstrSOLID%contacts) ) then
        call fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchange )
        call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )
        if( infoCTChange%contact2free+infoCTChange%contact2neighbor+infoCTChange%free2contact > 0 ) &
          ctchange = .true.
      endif
      if( fstr_is_contact_active() ) then
        gnt(2)=gnt(2)/iter
        convg = fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH)
      endif

      ! ----- update the total displacement
      ! u_{n+1} = u_{n} + \Delta u_{n+1}
      do i=1,hecMESH%n_node*ndof
        fstrSOLID%unode(i) = fstrSOLID%unode(i)+fstrSOLID%dunode(i)
      enddo

      if( convg .and. (.not.ctchange) ) exit

    enddo
    ! ----- end of augmentation loop

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    call cpu_time(tt)

    if( hecMESH%my_rank == 0 ) then
      write(ISTA,'("### Converged in NR ietration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
    endif

    deallocate(coord)
  end subroutine fstr_Newton_contactALag


  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson method.
  !> Standard Lagrange multiplier algorithm for contact analysis is incoluded in this subroutine.
  subroutine fstr_Newton_contactSLag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT,                  &
    restart_step_num, restart_substep_num, sub_step, infoCTChange, conMAT )

    use mContact
    use m_addContactStiffness
    use m_solve_LINEQ_contact

    integer, intent(in)                    :: cstep        !< current loading step
    type (hecmwST_local_mesh)              :: hecMESH      !< hecmw mesh
    type (hecmwST_matrix)                  :: hecMAT       !< hecmw matrix
    type (fstr_solid)                      :: fstrSOLID    !< fstr_solid
    integer, intent(in)                    :: sub_step     !< substep number of current loading step
    type (fstr_param)                      :: fstrPARAM    !< type fstr_param
    type (fstr_info_contactChange)         :: infoCTChange !< fstr_info_contactChange
    type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange
    type (hecmwST_matrix), optional        :: conMAT

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
    real(kind=kreal), pointer :: coord(:)

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

    allocate(coord(hecMESH%n_node*ndof))

    tincr = fstrSOLID%step_ctrl(cstep)%initdt
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.0d0

    call cpu_time(tt0)

    fstrSOLID%dunode(:)  = 0.0d0

    if( cstep==restart_step_num.and.sub_step==restart_substep_num  ) then
      call fstr_save_originalMatrixStructure(hecMAT)
      call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )
      if(paraContactFlag.and.present(conMAT)) then
        call copyClearMatrix(hecMAT,conMAT)
      endif
      if ( fstr_is_contact_active() ) then
        ! ----  For Parallel Contact with Multi-Partition Domains
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_mat_con_contact(cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange, conMAT)
        else
          call fstr_mat_con_contact(cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)
        endif
      elseif( hecMAT%Iarray(99)==4 ) then
        write(*, *) ' This type of direct solver is not yet available in such case ! '
        write(*, *) ' Please change the solver type to intel MKL direct solver !'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH)
      call solve_LINEQ_contact_init(hecMESH, hecMAT, fstrMAT, is_mat_symmetric)
    endif

    stepcnt = 0

    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    if( paraContactFlag.and.present(conMAT) ) then
      conMAT%B(:) = 0.0D0
    endif

    if( fstr_is_contact_active() )  then
      ! ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        call fstr_ass_load_contact(cstep, hecMESH, conMAT, fstrSOLID, fstrMAT)
      else
        call fstr_ass_load_contact(cstep, hecMESH, hecMAT, fstrSOLID, fstrMAT)
      endif
    endif

    fstrSOLID%dunode(:) = 0.0d0

    max_iter_contact = 10

    count_step = 0

    loopFORcontactAnalysis: DO WHILE( .TRUE. )
      count_step = count_step+1

      ! ----- Inner Iteration
      res1   = 0.d0
      relres = 1.d0

      do iter = 1, fstrSOLID%step_ctrl(cstep)%max_iter
        call hecmw_BARRIER(hecMESH)
        if( myrank == 0 ) print *,'-------------------------------------------------'
        call hecmw_BARRIER(hecMESH)
        stepcnt = stepcnt+1

        call fstr_StiffMatrix(hecMESH, hecMAT, fstrSOLID, tincr)
        call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

        if( paraContactFlag .and. present(conMAT) ) then
          call hecMAT_clear( conMAT )
        endif

        if( fstr_is_contact_active() ) then
          ! ---- For Parallel Contact with Multi-Partition Domains
          if( paraContactFlag .and. present(conMAT) ) then
            call fstr_AddContactStiffness(cstep,iter,conMAT,fstrMAT,fstrSOLID)
          else
            call fstr_AddContactStiffness(cstep,iter,hecMAT,fstrMAT,fstrSOLID)
          endif
        endif

        ! ----- Set Boundary condition
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, stepcnt, conMAT)
        else
          call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, stepcnt)
        endif

        nndof = hecMAT%N*hecMAT%ndof

        !----- SOLVE [Kt]{du}={R}
        ! ----  For Parallel Contact with Multi-Partition Domains
        hecMAT%X = 0.0d0
        call fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
        if(paraContactFlag.and.present(conMAT)) then
          q_residual = fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH)
          call solve_LINEQ_contact(hecMESH, hecMAT, fstrMAT, 1.0D0, conMAT)
        else
          q_residual = fstr_get_norm_contact('residualForce',hecMESH,hecMAT,fstrSOLID,fstrMAT)
          call solve_LINEQ_contact(hecMESH, hecMAT, fstrMAT)
        endif
        call fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)

        x_residual = fstr_get_x_norm_contact(hecMAT,fstrMAT,hecMESH)

        ! ----- update external nodal displacement increments
        if(paraContactFlag.and.present(conMAT)) then
          call paraContact_update_3_R(hecMESH,hecMAT%X)
        else
          call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
        endif
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
          do i = 1, fstrMAT%num_lagrange
            fstrMAT%lagrange(i) = fstrMAT%lagrange(i)+hecMAT%X(hecMESH%n_node*ndof+i)
          enddo
        endif

        ! ----- update the strain, stress, and internal force (only QFORCE)
        call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID,tincr,iter)

        ! ----- Set residual
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT )
        else
          call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID)
        endif

        if( fstr_is_contact_active() )  then
          if(paraContactFlag.and.present(conMAT)) then
            conMAT%B(1:3*hecMAT%NP) = 0.0D0
            call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID,conMAT)
          else
            call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID)
          endif
        endif

        if( paraContactFlag .and. present(conMAT) ) then
          res = fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH)
        else
          res = fstr_get_norm_contact('residualForce',hecMESH,hecMAT,fstrSOLID,fstrMAT)
        endif

        res = sqrt(res)/n_node_global
        if( iter == 1 ) res0 = res
        if( res0 == 0.0d0 ) then
          res0 =1.0d0
        else
          relres = dabs( res1-res )/res0
        endif
        if( hecMESH%my_rank == 0 ) then
          write(*, '(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
          write(ISTA, '(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
        endif

        ! ----- check convergence
        if( res < fstrSOLID%step_ctrl(cstep)%converg  .or.     &
          relres < fstrSOLID%step_ctrl(cstep)%converg ) exit
        res1 = res

      enddo
      ! ----- end of inner loop

      ! -----  not convergence
      if( iter > fstrSOLID%step_ctrl(cstep)%max_iter ) then
        if( hecMESH%my_rank == 0) then
          write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
          write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
        end if
        stop
      end if

      call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )

      if( hecMAT%Iarray(99) == 4 .and. .not. fstr_is_contact_active() ) then
        write(*, *) ' This type of direct solver is not yet available in such case ! '
        write(*, *) ' Please use intel MKL direct solver !'
        call  hecmw_abort( hecmw_comm_get_comm() )
      endif

      is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID, hecMESH)
      contact_changed_global = 0

      if( fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH) ) then
        exit loopFORcontactAnalysis
      elseif( fstr_is_matrixStructure_changed(infoCTChange) ) then
        ! ----  For Parallel Contact with Multi-Partition Domains
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange, conMAT)
        else
          call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)
        endif
        contact_changed_global = 1
      else

      endif

      call hecmw_allreduce_I1(hecMESH, contact_changed_global, HECMW_MAX)

      if (contact_changed_global > 0) then
        hecMAT%B(:) = 0.0D0
        if( paraContactFlag .and. present(conMAT) ) then
          conMAT%B(:) = 0.0D0
        endif
        call solve_LINEQ_contact_init(hecMESH, hecMAT, fstrMAT, is_mat_symmetric)
      endif
      if( count_step > max_iter_contact ) exit loopFORcontactAnalysis

    ENDDO loopFORcontactAnalysis

    ! ----- update the total displacement
    !       u_{n+1} = u_{n} + \Delta u_{n+1}
    do i = 1, hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i)+fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState(hecMESH, fstrSOLID, tincr)

    call cpu_time(tt)
    if( hecMESH%my_rank == 0) then
      write(ISTA,'("### Converged in contact iteration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
    endif

    deallocate(coord)
  end subroutine fstr_Newton_contactSLag


end module m_fstr_NonLinearMethod

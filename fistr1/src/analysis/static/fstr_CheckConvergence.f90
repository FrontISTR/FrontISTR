!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a unified convergence check for Newton iteration.
!>
!> The subroutine fstr_check_convergence computes residual norms from the
!> assembled residual vector and checks convergence/divergence criteria.
!> It replaces the former fstr_check_iteration_converged and various inline
!> convergence checks in contact and dynamic analysis routines.

module m_fstr_IterationControl
  use m_fstr
  implicit none

  private
  public :: fstr_check_convergence
  public :: fstr_check_convergence_main

contains

  !> \brief Wrapper that calls fstr_check_convergence_main and applies the
  !>        common divergence/NaN handling (status classification, failure logging,
  !>        fstrSOLID stats update).
  !>
  !> The body fstr_check_convergence_main is kept separate so that customized
  !> convergence criteria can be implemented by swapping/extending it without
  !> touching the failure-handling boilerplate here.
  subroutine fstr_check_convergence( &
      hecMESH, hecMAT, fstrSOLID, fstrPR, &
      ndof, iter, sub_step, cstep, &
      residual_vec, nresid, &
      resb, res_prev, &
      n_node_global, &
      iterStatus, &
      maxDLag, converg_dlag )
    implicit none

    type(hecmwST_local_mesh), intent(in)    :: hecMESH
    type(hecmwST_matrix), intent(in)        :: hecMAT
    type(fstr_solid), intent(inout)         :: fstrSOLID
    type(fstr_param), intent(in)            :: fstrPR
    integer(kind=kint), intent(in)          :: ndof
    integer(kind=kint), intent(in)          :: iter
    integer(kind=kint), intent(in)          :: sub_step
    integer(kind=kint), intent(in)          :: cstep
    real(kind=kreal), intent(in)            :: residual_vec(:)
    integer(kind=kint), intent(in)          :: nresid
    real(kind=kreal), intent(inout)         :: resb
    real(kind=kreal), intent(inout)         :: res_prev
    integer(kind=kint), intent(in)          :: n_node_global
    integer(kind=kint), intent(out)         :: iterStatus
    real(kind=kreal), intent(in), optional  :: maxDLag
    real(kind=kreal), intent(in), optional  :: converg_dlag

    real(kind=kreal)    :: res_for_check
    logical             :: do_failure_check, is_dynamic, is_contact

    is_dynamic = (fstrPR%solution_type == kstDYNAMIC)
    is_contact = (n_node_global > 0)

    ! --- core convergence check (customizable) ---
    call fstr_check_convergence_main( &
        hecMESH, hecMAT, fstrSOLID, fstrPR, &
        ndof, iter, cstep, &
        residual_vec, nresid, &
        resb, res_prev, &
        n_node_global, &
        iterStatus, &
        do_failure_check, res_for_check, &
        maxDLag, converg_dlag )

    if( iterStatus == kitrConverged ) return
    if( .not. do_failure_check ) return

    ! --- common divergence / NaN classification ---
    if( res_for_check /= res_for_check ) then
      iterStatus = kitrFloatingError
    else if( iter == fstrSOLID%step_ctrl(cstep)%max_iter .or. &
             res_for_check > fstrSOLID%step_ctrl(cstep)%maxres ) then
      iterStatus = kitrDiverged
    endif

    if( iterStatus == kitrContinue ) return

    ! --- common failure handling: log + stats update ---
    if( hecMESH%my_rank == 0 ) then
      ! Static non-contact path (legacy fstr_check_iteration_converged) also
      ! wrote to ILOG; preserve that behaviour for backward compatibility.
      if( .not. is_dynamic .and. .not. is_contact ) then
        write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
      endif
      write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
    endif
    fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT), iter)
    fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter
    fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1
    if( iterStatus == kitrDiverged .and. &
        iter == fstrSOLID%step_ctrl(cstep)%max_iter ) then
      fstrSOLID%NRstat_i(knstDRESN) = 1
    else
      ! kitrDiverged due to maxres, or kitrFloatingError due to NaN
      fstrSOLID%NRstat_i(knstDRESN) = 2
    endif

  end subroutine fstr_check_convergence

  !> \brief Core convergence check: computes residual norm, applies the per-path
  !>        convergence criterion, and (when not converged) returns the residual
  !>        value used to decide divergence in the caller.
  !>
  !> Branching is by fstrPR%solution_type (static/dynamic) and n_node_global (>0 = contact path).
  !> Currently used paths:
  !>   - Static, non-contact  : fstr_Newton, fstr_QuasiNewton
  !>   - Static, contact      : fstr_Newton_contactALag, fstr_Newton_contactSLag
  !>   - Dynamic, contact     : fstr_Newton_dynamic_contactSLag
  !> (Dynamic non-contact path is currently not used; all dynamic Newton drivers
  !>  go through the contact-aware routine.)
  !>
  !> \param[in]    hecMESH       mesh
  !> \param[in]    hecMAT        matrix (B=residual, X=solution increment)
  !> \param[inout] fstrSOLID     solid data (QFORCE, dunode)
  !> \param[in]    fstrPR        global parameters (solution_type, nlgeom)
  !> \param[in]    ndof          degrees of freedom per node
  !> \param[in]    iter          current Newton iteration number
  !> \param[in]    cstep         current loading step number
  !> \param[in]    residual_vec  assembled residual vector (used by contact paths)
  !> \param[in]    nresid        size of residual_vec
  !> \param[inout] resb          reference residual for self-normalization (set at iter=1)
  !> \param[inout] res_prev      previous normalized residual (for relres computation)
  !> \param[in]    n_node_global total internal node count (>0 for contact path, 0 for non-contact)
  !> \param[out]   iterStatus    kitrConverged or kitrContinue (failure paths handled by caller)
  !> \param[out]   do_failure_check  true if caller should run divergence/NaN check
  !> \param[out]   res_for_check     residual value the caller should use for divergence check
  !> \param[in]    maxDLag       (optional) max change of Lagrange multiplier
  !> \param[in]    converg_dlag  (optional) threshold for maxDLag
  subroutine fstr_check_convergence_main( &
      hecMESH, hecMAT, fstrSOLID, fstrPR, &
      ndof, iter, cstep, &
      residual_vec, nresid, &
      resb, res_prev, &
      n_node_global, &
      iterStatus, &
      do_failure_check, res_for_check, &
      maxDLag, converg_dlag )
    implicit none

    type(hecmwST_local_mesh), intent(in)    :: hecMESH
    type(hecmwST_matrix), intent(in)        :: hecMAT
    type(fstr_solid), intent(inout)         :: fstrSOLID
    type(fstr_param), intent(in)            :: fstrPR
    integer(kind=kint), intent(in)          :: ndof
    integer(kind=kint), intent(in)          :: iter
    integer(kind=kint), intent(in)          :: cstep
    real(kind=kreal), intent(in)            :: residual_vec(:)
    integer(kind=kint), intent(in)          :: nresid
    real(kind=kreal), intent(inout)         :: resb
    real(kind=kreal), intent(inout)         :: res_prev
    integer(kind=kint), intent(in)          :: n_node_global
    integer(kind=kint), intent(out)         :: iterStatus
    logical, intent(out)                    :: do_failure_check
    real(kind=kreal), intent(out)           :: res_for_check
    real(kind=kreal), intent(in), optional  :: maxDLag
    real(kind=kreal), intent(in), optional  :: converg_dlag

    real(kind=kreal) :: res_sq, res_nrm, qnrm, rres, xnrm, dunrm, rxnrm
    real(kind=kreal) :: relres, res_normalized
    logical :: is_dynamic, is_contact

    iterStatus = kitrContinue
    do_failure_check = .false.
    res_for_check = 0.0d0

    is_dynamic = (fstrPR%solution_type == kstDYNAMIC)
    is_contact = (n_node_global > 0)

    ! =========================================================================
    ! PATH 1: Static, non-contact  (QFORCE normalization + disp correction)
    ! =========================================================================
    if( .not. is_dynamic .and. .not. is_contact ) then
      call hecmw_InnerProduct_R(hecMESH, ndof, hecmw_mat_get_B(hecMAT), hecmw_mat_get_B(hecMAT), res_sq)
      res_nrm = sqrt(res_sq)

      call hecmw_InnerProduct_R(hecMESH, ndof, hecmw_mat_get_X(hecMAT), hecmw_mat_get_X(hecMAT), xnrm)
      xnrm = sqrt(xnrm)

      call hecmw_InnerProduct_R(hecMESH, ndof, fstrSOLID%QFORCE, fstrSOLID%QFORCE, qnrm)
      qnrm = sqrt(qnrm)
      if( qnrm < 1.0d-8 ) qnrm = 1.0d0

      if( iter == 1 ) then
        dunrm = xnrm
      else
        call hecmw_InnerProduct_R(hecMESH, ndof, fstrSOLID%dunode, fstrSOLID%dunode, dunrm)
        dunrm = sqrt(dunrm)
      endif

      rres = res_nrm / qnrm
      rxnrm = xnrm / dunrm

      if( hecMESH%my_rank == 0 ) then
        if( qnrm == 1.0d0 ) then
          write(*,"(a,i8,a,1pe11.4,a,1pe11.4)") " iter:", iter, ", residual(abs):", rres, ", disp.corr.:", rxnrm
        else
          write(*,"(a,i8,a,1pe11.4,a,1pe11.4)") " iter:", iter, ", residual:", rres, ", disp.corr.:", rxnrm
        endif
      endif

      ! Convergence (skip if linear solver flagged divergence)
      if( hecmw_mat_get_flag_diverged(hecMAT) == kNO ) then
        if( rres < fstrSOLID%step_ctrl(cstep)%converg .or. &
            rxnrm < fstrSOLID%step_ctrl(cstep)%converg_ddisp ) then
          iterStatus = kitrConverged
          return
        endif
      endif

      do_failure_check = .true.
      res_for_check = rres

    ! =========================================================================
    ! PATH 2: Static, contact  (self-normalization by n_node_global)
    ! =========================================================================
    else if( .not. is_dynamic .and. is_contact ) then
      res_sq = dot_product(residual_vec(1:nresid), residual_vec(1:nresid))
      call hecmw_allreduce_R1(hecMESH, res_sq, hecmw_sum)
      res_normalized = sqrt(res_sq) / dble(n_node_global)

      if( iter == 1 ) resb = res_normalized
      if( resb == 0.0d0 ) then
        resb = 1.0d0
      else
        relres = dabs(res_prev - res_normalized) / resb
      endif

      if( hecMESH%my_rank == 0 ) then
        if( iter == 1 ) then
          write(*, '(a,i3,a,2e15.7)') ' - Residual(', iter, ') =', res_normalized, 0.0d0
        else
          write(*, '(a,i3,a,2e15.7)') ' - Residual(', iter, ') =', res_normalized, relres
        endif
      endif

      ! Convergence (iter==1: absolute only; iter>1: also relative)
      if( iter > 1 ) then
        if( res_normalized < fstrSOLID%step_ctrl(cstep)%converg .or. &
            relres < fstrSOLID%step_ctrl(cstep)%converg_ddisp ) then
          iterStatus = kitrConverged
          res_prev = res_normalized
          return
        endif
      else
        if( res_normalized < fstrSOLID%step_ctrl(cstep)%converg ) then
          iterStatus = kitrConverged
          res_prev = res_normalized
          return
        endif
      endif

      res_prev = res_normalized

      do_failure_check = .true.
      res_for_check = res_normalized

    ! =========================================================================
    ! PATH 3: Dynamic, contact  (self-normalization + maxDLag condition)
    ! =========================================================================
    else if( is_dynamic .and. is_contact ) then
      res_sq = dot_product(residual_vec(1:nresid), residual_vec(1:nresid))
      call hecmw_allreduce_R1(hecMESH, res_sq, hecmw_sum)

      if( iter == 1 ) resb = res_sq
      res_normalized = dsqrt(res_sq / resb)

      ! Only check when nlgeom and ndof /= 4
      if( fstrPR%nlgeom .and. ndof /= 4 ) then
        if( hecMESH%my_rank == 0 ) then
          write(*,'(a,i5,a,1pe12.4)') "iter: ", iter, ", res: ", res_normalized
          write(ISTA,'(''iter='',I5,''- Residual'',E15.7)') iter, res_normalized
          if( present(maxDLag) ) then
            write(*,'(a,1e15.7)') ' - MaxDLag =', maxDLag
            write(ISTA,'(a,1e15.7)') ' - MaxDLag =', maxDLag
          endif
        endif
        if( present(maxDLag) .and. present(converg_dlag) ) then
          if( res_normalized < fstrSOLID%step_ctrl(cstep)%converg .and. &
              maxDLag < converg_dlag ) then
            iterStatus = kitrConverged
            return
          endif
        else
          if( res_normalized < fstrSOLID%step_ctrl(cstep)%converg ) then
            iterStatus = kitrConverged
            return
          endif
        endif

        do_failure_check = .true.
        res_for_check = res_normalized
      endif

    endif

  end subroutine fstr_check_convergence_main

end module m_fstr_IterationControl

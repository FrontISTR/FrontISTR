!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a unified convergence check for Newton iteration.
!>
!> fstr_check_convergence measures the residual against a reference force assembled from the internal force and, in
!> dynamic analysis, the inertia and viscous force. Quantities of different physical dimension (translational force,
!> moment, Lagrange multiplier) are evaluated separately instead of being accumulated into a single norm.
!> fstr_check_linear_solver does the same for the linear solver result, replacing the per-driver inline checks that
!> followed solve_LINEQ / solve_LINEQ_contact.

module m_fstr_IterationControl
  use m_fstr
  implicit none

  private
  public :: fstr_convergence_state
  public :: fstr_check_convergence
  public :: fstr_check_convergence_main
  public :: fstr_check_linear_solver

  !> Norms evaluated by the latest convergence check.
  !>
  !> The variable is held by the Newton driver so that the quantities the decision was made from remain available to
  !> the caller after the check returns.
  type fstr_convergence_state
    real(kind=kreal)   :: fref_t   !< reference force of translational DOFs
    real(kind=kreal)   :: fref_r   !< reference moment of rotational DOFs
    real(kind=kreal)   :: res_t    !< translational residual norm / fref_t
    real(kind=kreal)   :: res_r    !< rotational residual norm / fref_r
    real(kind=kreal)   :: dx_t     !< translational correction norm / increment norm
    real(kind=kreal)   :: dx_r     !< rotational correction norm / increment norm
    real(kind=kreal)   :: dx_l     !< Lagrange correction norm / multiplier norm
    logical            :: has_rot  !< rotational DOFs are evaluated
    logical            :: has_lag  !< Lagrange rows exist in the model
    logical            :: has_dx   !< correction norms are evaluated
    logical            :: abs_t    !< fref_t fell back to the floor, res_t is an absolute value
    logical            :: abs_r    !< fref_r fell back to the floor, res_r is an absolute value
  end type fstr_convergence_state

  !> Lower bound of the reference force. Below it the residual criterion degenerates into an absolute one.
  real(kind=kreal), parameter :: FREF_FLOOR = 1.0d-8

contains

  !> \brief Wrapper that calls fstr_check_convergence_main and applies the common divergence/NaN handling
  !>        (status classification, failure logging, fstrSOLID stats update).
  !>
  !> The body is kept separate so that customized convergence criteria can be implemented by swapping or extending
  !> fstr_check_convergence_main without touching the failure-handling boilerplate here.
  subroutine fstr_check_convergence( hecMESH, hecMAT, fstrSOLID, fstrPR, ndof, iter, sub_step, cstep, &
      residual_vec, cnvstat, iterStatus, hecLagMAT )
    implicit none

    type(hecmwST_local_mesh), intent(in)      :: hecMESH
    type(hecmwST_matrix), intent(in)          :: hecMAT
    type(fstr_solid), intent(inout)           :: fstrSOLID
    type(fstr_param), intent(in)              :: fstrPR
    integer(kind=kint), intent(in)            :: ndof
    integer(kind=kint), intent(in)            :: iter
    integer(kind=kint), intent(in)            :: sub_step
    integer(kind=kint), intent(in)            :: cstep
    real(kind=kreal), intent(in)              :: residual_vec(:)
    type(fstr_convergence_state), intent(out) :: cnvstat
    integer(kind=kint), intent(out)           :: iterStatus
    type(hecmwST_matrix_lagrange), intent(in), optional :: hecLagMAT

    real(kind=kreal)    :: res_for_check
    logical             :: do_failure_check

    ! --- core convergence check (customizable) ---
    call fstr_check_convergence_main( hecMESH, hecMAT, fstrSOLID, fstrPR, ndof, iter, cstep, &
        residual_vec, cnvstat, iterStatus, do_failure_check, res_for_check, hecLagMAT )

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
      write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
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

  !> \brief Core convergence check.
  !>
  !> Residual criterion (evaluated per dimension):
  !>   res_t = ||R_t|| / max(||Q_t||, ||D_t||)  <  CONVERG
  !>   res_r = ||R_r|| / max(||Q_r||, ||D_r||)  <  CONVERG        (ndof=6 only)
  !> Correction criterion:
  !>   dx_t  = ||du_t|| / ||Du_t||              <  CONVERG_DDISP
  !>   dx_r  = ||du_r|| / ||Du_r||              <  CONVERG_DDISP  (ndof=6 only)
  !>   dx_l  = ||dlambda|| / ||lambda||         <  CONVERG_LAG    (Lagrange rows only)
  !> and the decision is
  !>   converged = dx_l .and. ( all dx .or. all res )
  !>
  !> Q is the internal force and D the inertia plus viscous force, so that the reference force stays finite for a body
  !> in free vibration, where the vector sum Q+D vanishes at equilibrium while neither norm does. Pressure DOFs of
  !> ndof=4 and the Lagrange rows of the residual carry a dimension of their own and are represented by dx_l instead of
  !> entering the residual norms.
  !>
  !> In dynamic analysis the check precedes the linear solve, so at iter=1 the correction vector belongs to no
  !> iteration yet and the state is the one before any correction has been applied; nothing is decided there.
  !>
  !> \param[in]    hecMESH       mesh
  !> \param[in]    hecMAT        matrix (X=solution increment)
  !> \param[inout] fstrSOLID     solid data (QFORCE, DFORCE, dunode)
  !> \param[in]    fstrPR        global parameters (solution_type)
  !> \param[in]    ndof          degrees of freedom per node
  !> \param[in]    iter          current Newton iteration number
  !> \param[in]    cstep         current loading step number
  !> \param[in]    residual_vec  assembled residual vector
  !> \param[out]   cnvstat       norms the decision was made from
  !> \param[out]   iterStatus    kitrConverged or kitrContinue (failure paths handled by caller)
  !> \param[out]   do_failure_check  true if caller should run divergence/NaN check
  !> \param[out]   res_for_check     residual value the caller should use for divergence check
  !> \param[in]    hecLagMAT     (optional) Lagrange multipliers of contact analysis
  subroutine fstr_check_convergence_main( hecMESH, hecMAT, fstrSOLID, fstrPR, ndof, iter, cstep, &
      residual_vec, cnvstat, iterStatus, do_failure_check, res_for_check, hecLagMAT )
    implicit none

    type(hecmwST_local_mesh), intent(in)      :: hecMESH
    type(hecmwST_matrix), intent(in)          :: hecMAT           !< X=solution increment
    type(fstr_solid), intent(inout)           :: fstrSOLID
    type(fstr_param), intent(in)              :: fstrPR
    integer(kind=kint), intent(in)            :: ndof
    integer(kind=kint), intent(in)            :: iter
    integer(kind=kint), intent(in)            :: cstep
    real(kind=kreal), intent(in)              :: residual_vec(:)
    type(fstr_convergence_state), intent(out) :: cnvstat
    integer(kind=kint), intent(out)           :: iterStatus
    logical, intent(out)                      :: do_failure_check
    real(kind=kreal), intent(out)             :: res_for_check
    type(hecmwST_matrix_lagrange), intent(in), optional :: hecLagMAT

    integer(kind=kint), parameter :: NSUM = 13
    real(kind=kreal)   :: sq(NSUM)
    integer(kind=kint) :: i, npndof, num_lagrange
    real(kind=kreal)   :: converg, converg_ddisp, converg_lag
    logical            :: is_dynamic, ok_res, ok_dx, ok_lag

    iterStatus = kitrContinue
    do_failure_check = .false.
    res_for_check = 0.0d0

    is_dynamic = (fstrPR%solution_type == kstDYNAMIC)
    converg = fstrSOLID%step_ctrl(cstep)%converg
    converg_ddisp = fstrSOLID%step_ctrl(cstep)%converg_ddisp
    converg_lag = fstrSOLID%step_ctrl(cstep)%converg_lag

    num_lagrange = 0
    if( present(hecLagMAT) ) num_lagrange = hecLagMAT%num_lagrange
    npndof = hecMAT%NP*ndof

    cnvstat%has_rot = ( ndof == 6 )
    cnvstat%has_dx = .not. ( is_dynamic .and. iter == 1 )

    ! --- squares of the norms, reduced in a single collective ---
    sq(:) = 0.0d0
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, residual_vec,      sq(1), sq(2) )
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, fstrSOLID%QFORCE,  sq(3), sq(4) )
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, fstrSOLID%DFORCE,  sq(5), sq(6) )
    if( cnvstat%has_dx ) then
      call fstr_get_sqnorm_dofgroup( hecMESH, ndof, hecMAT%X,         sq(7), sq(8) )
      call fstr_get_sqnorm_dofgroup( hecMESH, ndof, fstrSOLID%dunode, sq(9), sq(10) )
      do i = 1, num_lagrange
        sq(11) = sq(11) + hecMAT%X(npndof+i)*hecMAT%X(npndof+i)
        sq(12) = sq(12) + hecLagMAT%Lagrange(i)*hecLagMAT%Lagrange(i)
      enddo
    endif
    sq(13) = dble(num_lagrange)
    call hecmw_allreduce_R( hecMESH, sq, NSUM, hecmw_sum )

    ! Lagrange rows may be absent from a subdomain while present in the model, so the criteria must be selected from
    ! the reduced count.
    cnvstat%has_lag = ( sq(13) > 0.5d0 )

    ! --- residual relative to the reference force ---
    cnvstat%fref_t = max( sqrt(sq(3)), sqrt(sq(5)) )
    cnvstat%abs_t = ( cnvstat%fref_t < FREF_FLOOR )
    if( cnvstat%abs_t ) cnvstat%fref_t = 1.0d0
    cnvstat%res_t = sqrt(sq(1)) / cnvstat%fref_t

    cnvstat%fref_r = 1.0d0
    cnvstat%abs_r = .false.
    cnvstat%res_r = 0.0d0
    if( cnvstat%has_rot ) then
      cnvstat%fref_r = max( sqrt(sq(4)), sqrt(sq(6)) )
      cnvstat%abs_r = ( cnvstat%fref_r < FREF_FLOOR )
      if( cnvstat%abs_r ) cnvstat%fref_r = 1.0d0
      cnvstat%res_r = sqrt(sq(2)) / cnvstat%fref_r
    endif

    ! --- correction relative to the accumulated increment ---
    ! A vanishing denominator leaves the ratio at 1, i.e. not converged.
    cnvstat%dx_t = 1.0d0
    cnvstat%dx_r = 0.0d0
    cnvstat%dx_l = 0.0d0
    if( cnvstat%has_dx ) then
      if( sq(9) > 0.0d0 ) cnvstat%dx_t = sqrt( sq(7)/sq(9) )
      if( cnvstat%has_rot ) then
        cnvstat%dx_r = 1.0d0
        if( sq(10) > 0.0d0 ) cnvstat%dx_r = sqrt( sq(8)/sq(10) )
      endif
      if( cnvstat%has_lag ) then
        if( sq(12) > 0.0d0 ) then
          cnvstat%dx_l = sqrt( sq(11)/sq(12) )
        else if( sq(11) > 0.0d0 ) then
          cnvstat%dx_l = 1.0d0
        endif
      endif
    endif

    ! A NaN component fails every comparison above: it never passes the test on a denominator, so a ratio would keep
    ! its default. The sums are the only place it survives, and every norm takes it over.
    do i = 1, NSUM-1
      if( sq(i) /= sq(i) ) then
        cnvstat%res_t = sq(i)
        cnvstat%res_r = sq(i)
        cnvstat%dx_t = sq(i)
        cnvstat%dx_r = sq(i)
        cnvstat%dx_l = sq(i)
      endif
    enddo

    if( hecMESH%my_rank == 0 ) call fstr_print_convergence_state( iter, cnvstat )

    if( .not. cnvstat%has_dx ) return

    ok_res = ( cnvstat%res_t < converg )
    ok_dx = ( cnvstat%dx_t < converg_ddisp )
    if( cnvstat%has_rot ) then
      ok_res = ok_res .and. ( cnvstat%res_r < converg )
      ok_dx = ok_dx .and. ( cnvstat%dx_r < converg_ddisp )
    endif
    ok_lag = .true.
    if( cnvstat%has_lag ) ok_lag = ( cnvstat%dx_l < converg_lag )

    if( ok_lag .and. ( ok_dx .or. ok_res ) ) then
      iterStatus = kitrConverged
      return
    endif

    do_failure_check = .true.
    res_for_check = max( cnvstat%res_t, cnvstat%res_r )

  end subroutine fstr_check_convergence_main

  !> \brief Sum of squares of a nodal vector over the internal nodes, split into translational (DOF 1-3) and
  !>        rotational (DOF 4-6) components.
  !>
  !> The result is not reduced across subdomains; the caller reduces all norms at once. DOF 4 of ndof=4, the pressure,
  !> belongs to neither group.
  subroutine fstr_get_sqnorm_dofgroup( hecMESH, ndof, vec, sq_t, sq_r )
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in)       :: ndof
    real(kind=kreal), intent(in)         :: vec(:)
    real(kind=kreal), intent(out)        :: sq_t
    real(kind=kreal), intent(out)        :: sq_r

    integer(kind=kint) :: i, idof, idx

    sq_t = 0.0d0
    sq_r = 0.0d0
    do i = 1, hecMESH%nn_internal
      idx = ndof*(i-1)
      do idof = 1, min(ndof,3)
        sq_t = sq_t + vec(idx+idof)*vec(idx+idof)
      enddo
      if( ndof /= 6 ) cycle
      do idof = 4, 6
        sq_r = sq_r + vec(idx+idof)*vec(idx+idof)
      enddo
    enddo

  end subroutine fstr_get_sqnorm_dofgroup

  !> \brief Print the criteria of one Newton iteration, omitting those that do not apply to the analysis.
  subroutine fstr_print_convergence_state( iter, cnvstat )
    implicit none
    integer(kind=kint), intent(in)           :: iter
    type(fstr_convergence_state), intent(in) :: cnvstat

    character(len=256) :: line

    write(line,'(a,i8)') " iter:", iter
    if( cnvstat%abs_t ) then
      write(line(len_trim(line)+1:),'(a,1pe11.4)') ", res(force,abs):", cnvstat%res_t
    else
      write(line(len_trim(line)+1:),'(a,1pe11.4)') ", res(force):", cnvstat%res_t
    endif
    if( cnvstat%has_rot ) then
      if( cnvstat%abs_r ) then
        write(line(len_trim(line)+1:),'(a,1pe11.4)') ", res(mom,abs):", cnvstat%res_r
      else
        write(line(len_trim(line)+1:),'(a,1pe11.4)') ", res(mom):", cnvstat%res_r
      endif
    endif
    if( cnvstat%has_dx ) then
      write(line(len_trim(line)+1:),'(a,1pe11.4)') ", disp.corr.:", cnvstat%dx_t
      if( cnvstat%has_rot ) then
        write(line(len_trim(line)+1:),'(a,1pe11.4)') ", rot.corr.:", cnvstat%dx_r
      endif
      if( cnvstat%has_lag ) then
        write(line(len_trim(line)+1:),'(a,1pe11.4)') ", lag.corr.:", cnvstat%dx_l
      endif
    endif
    write(*,'(a)') trim(line)

  end subroutine fstr_print_convergence_state

  !> \brief Classify the linear solver result and record why the solve failed.
  !>
  !> Called by every Newton driver right after the linear solve, so the drivers share one definition of a usable
  !> solution.
  subroutine fstr_check_linear_solver( hecMESH, hecMAT, fstrSOLID, cstep, sub_step, iterStatus, istat )
    implicit none

    type(hecmwST_local_mesh), intent(in)     :: hecMESH
    type(hecmwST_matrix), intent(in)         :: hecMAT      !< the matrix handed to the solver
    type(fstr_solid), intent(inout)          :: fstrSOLID
    integer(kind=kint), intent(in)           :: cstep
    integer(kind=kint), intent(in)           :: sub_step
    integer(kind=kint), intent(out)          :: iterStatus  !< kitrContinue if the solution is usable
    integer(kind=kint), intent(in), optional :: istat       !< status returned by solve_LINEQ_contact

    iterStatus = kitrContinue

    ! direct solvers report a genuine failure only through istat
    if( present(istat) ) then
      if( istat /= 0 ) then
        iterStatus = kitrDiverged
        fstrSOLID%NRstat_i(knstDRESN) = 4
      endif
    endif

    if( iterStatus == kitrContinue ) then
      if( hecmw_mat_get_flag_diverged(hecMAT) /= kNO ) then
        ! broke down: indefinite preconditioner/matrix, or NaN
        iterStatus = kitrDiverged
        fstrSOLID%NRstat_i(knstDRESN) = 6
      else if( hecmw_mat_get_solver_type(hecMAT) == 1 .and. &
               hecmw_mat_get_flag_converged(hecMAT) == kNO ) then
        ! ran out of iterations. Restricted to iterative solvers because a direct solver leaves the same flag unset
        ! whenever its residual exceeds the tolerance, and hecmw_solve only warns and keeps going in that case
        iterStatus = kitrDiverged
        fstrSOLID%NRstat_i(knstDRESN) = 5
      endif
    endif

    if( iterStatus == kitrContinue ) return

    if( hecMESH%my_rank == 0 ) then
      write(*,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
    endif
    fstrSOLID%CutBack_stat = fstrSOLID%CutBack_stat + 1

  end subroutine fstr_check_linear_solver

end module m_fstr_IterationControl

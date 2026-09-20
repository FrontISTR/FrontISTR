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
!>
!> Whether an iteration has converged is decided in fstr_decide_convergence and nowhere else. The other routines only
!> prepare the table of criteria it reads: setup_stepInfo_converg of m_step fixes, when the step is read, which
!> criteria are checked against which thresholds, and fstr_evaluate_convergence fills in their values and unchecks
!> the criteria that do not apply to the analysis or to the iteration.

module m_fstr_IterationControl
  use m_fstr
  implicit none

  private
  public :: fstr_convergence_measure
  public :: fstr_convergence_state
  public :: fstr_check_convergence
  public :: fstr_check_convergence_main
  public :: fstr_check_linear_solver

  !> One convergence criterion: a normalized value compared with its threshold.
  type fstr_convergence_measure
    real(kind=kreal)   :: value = 0.0d0      !< normalized value
    real(kind=kreal)   :: tol   = 0.0d0      !< threshold
    logical            :: check = .false.    !< the criterion takes part in the decision (CHECK) or not (SKIP)
    logical            :: ok    = .false.    !< check .and. value < tol
    logical            :: absolute = .false. !< the reference fell back to the floor, value is not normalized
    integer(kind=kint) :: node  = 0          !< global node ID where a max norm occurs, 0 if not located
    integer(kind=kint) :: dof   = 0          !< DOF where a max norm occurs
  end type fstr_convergence_measure

  !> Criteria evaluated by the latest convergence check, indexed by (quantity, DOF group, norm) with the kcnv*
  !> constants. A criterion that does not apply to the analysis or to the iteration is left unchecked.
  !>
  !> The variable is held by the Newton driver so that the quantities the decision was made from remain available to
  !> the caller after the check returns.
  type fstr_convergence_state
    type(fstr_convergence_measure) :: m(2,3,2)
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
  !> The criteria, indexed by (quantity, DOF group, norm), are
  !>   residual,   translation, L2:  ||R_t|| / max(||Q_t||, ||D_t||)
  !>   residual,   rotation,    L2:  ||R_r|| / max(||Q_r||, ||D_r||)   (ndof=6 only)
  !>   correction, translation, L2:  ||du_t|| / ||Du_t||
  !>   correction, rotation,    L2:  ||du_r|| / ||Du_r||               (ndof=6 only)
  !>   correction, Lagrange,    L2:  ||dlambda|| / max(||lambda||, Fref_t)   (Lagrange rows only)
  !> and the same ratios in the max norm, the largest absolute DOF component over the nodes. Each is compared with the
  !> threshold fixed for the step, and fstr_decide_convergence combines them.
  !>
  !> Q is the internal force and D the inertia plus viscous force, so that the reference force stays finite for a body
  !> in free vibration, where the vector sum Q+D vanishes at equilibrium while neither norm does. Pressure DOFs of
  !> ndof=4 and the Lagrange rows of the residual carry a dimension of their own and are represented by the Lagrange
  !> correction instead of entering the residual norms.
  !>
  !> The multiplier of a contact constraint is a contact force, so the translational reference force Fref_t bounds
  !> ||lambda|| from below in its criterion. Normalizing by ||lambda|| alone makes the criterion arbitrarily strict
  !> where the contact carries little of the load, and unsatisfiable at a node entering contact, where lambda is
  !> still zero while dlambda is not.
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
  !> \param[out]   cnvstat       criteria the decision was made from
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

    logical :: has_dx

    iterStatus = kitrContinue
    do_failure_check = .false.
    res_for_check = 0.0d0

    has_dx = .not. ( fstrPR%solution_type == kstDYNAMIC .and. iter == 1 )

    call fstr_evaluate_convergence( hecMESH, hecMAT, fstrSOLID, ndof, cstep, has_dx, residual_vec, cnvstat, hecLagMAT )

    if( hecMESH%my_rank == 0 ) call fstr_print_convergence_state( iter, cnvstat )

    if( .not. has_dx ) return

    if( fstr_decide_convergence( cnvstat ) ) then
      iterStatus = kitrConverged
      return
    endif

    do_failure_check = .true.
    res_for_check = max( cnvstat%m(kcnvResidual, kcnvTranslation, kcnvL2)%value, &
                         cnvstat%m(kcnvResidual, kcnvRotation, kcnvL2)%value )

  end subroutine fstr_check_convergence_main

  !> \brief Combine the criteria into the convergence decision.
  !>
  !> This is the only place where the criteria are combined. A group of criteria is judged by
  !>   cnv_all_ok(ms)     every checked criterion is satisfied; true when none is checked
  !>   cnv_any_check(ms)  at least one criterion is checked
  !> which follow ALL and ANY of an empty group. A group joined by .or. has to be guarded by cnv_any_check, since a
  !> group without a checked criterion would pass otherwise. The norms of one quantity are all required.
  logical function fstr_decide_convergence( cnvstat )
    implicit none
    type(fstr_convergence_state), intent(in) :: cnvstat

    type(fstr_convergence_measure) :: lag(2), dx(4), res(4)

    lag = [ cnvstat%m(kcnvCorrection, kcnvLagrange, :) ]
    dx  = [ cnvstat%m(kcnvCorrection, kcnvTranslation:kcnvRotation, :) ]
    res = [ cnvstat%m(kcnvResidual,   kcnvTranslation:kcnvRotation, :) ]

    fstr_decide_convergence = cnv_all_ok(lag) .and. &
      ( ( cnv_any_check(dx)  .and. cnv_all_ok(dx)  ) .or. &
        ( cnv_any_check(res) .and. cnv_all_ok(res) ) )

  end function fstr_decide_convergence

  !> \brief Every checked criterion of the group is satisfied. True when none is checked.
  logical function cnv_all_ok( ms )
    implicit none
    type(fstr_convergence_measure), intent(in) :: ms(:)

    cnv_all_ok = all( ms%ok .or. .not. ms%check )

  end function cnv_all_ok

  !> \brief At least one criterion of the group is checked.
  logical function cnv_any_check( ms )
    implicit none
    type(fstr_convergence_measure), intent(in) :: ms(:)

    cnv_any_check = any( ms%check )

  end function cnv_any_check

  !> \brief Evaluate the criteria of one Newton iteration and select those that take part in the decision.
  !>
  !> A criterion is checked when the step asks for it and it applies to the analysis and to the iteration: rotation
  !> for ndof=6, Lagrange for a model with Lagrange rows, correction when has_dx.
  subroutine fstr_evaluate_convergence( hecMESH, hecMAT, fstrSOLID, ndof, cstep, has_dx, residual_vec, cnvstat, &
      hecLagMAT )
    implicit none

    type(hecmwST_local_mesh), intent(in)      :: hecMESH
    type(hecmwST_matrix), intent(in)          :: hecMAT           !< X=solution increment
    type(fstr_solid), intent(in)              :: fstrSOLID
    integer(kind=kint), intent(in)            :: ndof
    integer(kind=kint), intent(in)            :: cstep
    logical, intent(in)                       :: has_dx           !< the correction belongs to this iteration
    real(kind=kreal), intent(in)              :: residual_vec(:)
    type(fstr_convergence_state), intent(out) :: cnvstat
    type(hecmwST_matrix_lagrange), intent(in), optional :: hecLagMAT

    integer(kind=kint), parameter :: NSUM = 13
    real(kind=kreal)   :: sq(NSUM)
    integer(kind=kint) :: i, npndof, num_lagrange
    real(kind=kreal)   :: dx_t, dx_r, fref_t
    logical            :: has_rot, has_lag

    num_lagrange = 0
    if( present(hecLagMAT) ) num_lagrange = hecLagMAT%num_lagrange
    npndof = hecMAT%NP*ndof

    has_rot = ( ndof == 6 )

    ! --- squares of the norms, reduced in a single collective ---
    sq(:) = 0.0d0
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, residual_vec,      sq(1), sq(2) )
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, fstrSOLID%QFORCE,  sq(3), sq(4) )
    call fstr_get_sqnorm_dofgroup( hecMESH, ndof, fstrSOLID%DFORCE,  sq(5), sq(6) )
    if( has_dx ) then
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
    has_lag = ( sq(13) > 0.5d0 )

    cnvstat%m(:,:,:)%tol = fstrSOLID%step_ctrl(cstep)%cnv_tol(:,:,:)
    cnvstat%m(:,:,:)%check = fstrSOLID%step_ctrl(cstep)%cnv_check(:,:,:)
    if( .not. has_rot ) cnvstat%m(:,kcnvRotation,:)%check = .false.
    if( .not. has_lag ) cnvstat%m(:,kcnvLagrange,:)%check = .false.
    if( .not. has_dx ) cnvstat%m(kcnvCorrection,:,:)%check = .false.

    ! --- residual relative to the reference force ---
    fref_t = max( sqrt(sq(3)), sqrt(sq(5)) )
    call fstr_set_relative_measure( cnvstat%m(kcnvResidual,kcnvTranslation,kcnvL2), &
        sqrt(sq(1)), fref_t )
    if( has_rot ) then
      call fstr_set_relative_measure( cnvstat%m(kcnvResidual,kcnvRotation,kcnvL2), &
          sqrt(sq(2)), max( sqrt(sq(4)), sqrt(sq(6)) ) )
    endif

    ! --- correction relative to the accumulated increment ---
    ! A vanishing denominator leaves the ratio at 1, i.e. not converged.
    if( has_dx ) then
      dx_t = 1.0d0
      if( sq(9) > 0.0d0 ) dx_t = sqrt( sq(7)/sq(9) )
      dx_r = 1.0d0
      if( sq(10) > 0.0d0 ) dx_r = sqrt( sq(8)/sq(10) )
      cnvstat%m(kcnvCorrection,kcnvTranslation,kcnvL2)%value = dx_t
      if( has_rot ) cnvstat%m(kcnvCorrection,kcnvRotation,kcnvL2)%value = dx_r
      ! The multiplier is a contact force, so the reference force bounds its norm from below.
      if( has_lag ) then
        call fstr_set_relative_measure( cnvstat%m(kcnvCorrection,kcnvLagrange,kcnvL2), &
            sqrt(sq(11)), max( sqrt(sq(12)), fref_t ) )
      endif
    endif

    ! The flags are the same on every subdomain, so the collectives of the max norms are skipped consistently.
    if( any( cnvstat%m(:,:,kcnvMax)%check ) ) then
      call fstr_evaluate_convergence_max( hecMESH, hecMAT, fstrSOLID, ndof, has_dx, residual_vec, cnvstat, &
          hecLagMAT )
    endif

    ! A NaN component fails every comparison above: it never wins a maximum and never passes the test on a denominator,
    ! so a ratio would keep its default. The sums are the only place it survives, and every criterion takes it over.
    do i = 1, NSUM-1
      if( sq(i) /= sq(i) ) cnvstat%m(:,:,:)%value = sq(i)
    enddo

    cnvstat%m(:,:,:)%ok = cnvstat%m(:,:,:)%check .and. ( cnvstat%m(:,:,:)%value < cnvstat%m(:,:,:)%tol )

  end subroutine fstr_evaluate_convergence

  !> \brief Evaluate the max-norm criteria and locate the largest residual and correction components.
  !>
  !> The max norm of a DOF group is the largest absolute DOF component over the internal nodes, and each criterion is
  !> the same ratio as its L2 counterpart with the norms replaced. Where several nodes hold the same maximum, the
  !> largest node ID is reported, with the DOF taken from the subdomain owning that node.
  subroutine fstr_evaluate_convergence_max( hecMESH, hecMAT, fstrSOLID, ndof, has_dx, residual_vec, cnvstat, &
      hecLagMAT )
    implicit none

    type(hecmwST_local_mesh), intent(in)        :: hecMESH
    type(hecmwST_matrix), intent(in)            :: hecMAT           !< X=solution increment
    type(fstr_solid), intent(in)                :: fstrSOLID
    integer(kind=kint), intent(in)              :: ndof
    logical, intent(in)                         :: has_dx           !< the correction belongs to this iteration
    real(kind=kreal), intent(in)                :: residual_vec(:)
    type(fstr_convergence_state), intent(inout) :: cnvstat
    type(hecmwST_matrix_lagrange), intent(in), optional :: hecLagMAT

    integer(kind=kint), parameter :: NMAX = 12
    ! entries of vmax located, and the criteria they belong to
    integer(kind=kint), parameter :: NLOC = 4
    integer(kind=kint), parameter :: IMAX_LOC(NLOC) = [ 1, 2, 7, 8 ]
    integer(kind=kint), parameter :: IQ_LOC(NLOC) = [ kcnvResidual, kcnvResidual, kcnvCorrection, kcnvCorrection ]
    integer(kind=kint), parameter :: IG_LOC(NLOC) = [ kcnvTranslation, kcnvRotation, kcnvTranslation, kcnvRotation ]
    real(kind=kreal)   :: vmax(NMAX), gmax(NMAX)
    integer(kind=kint) :: loc(2,NLOC), gnode(NLOC), gdof(NLOC)
    integer(kind=kint) :: i, k, npndof, num_lagrange
    real(kind=kreal)   :: dx_t, dx_r, fref_t

    num_lagrange = 0
    if( present(hecLagMAT) ) num_lagrange = hecLagMAT%num_lagrange
    npndof = hecMAT%NP*ndof

    vmax(:) = 0.0d0
    loc(:,:) = 0
    call fstr_get_maxabs_dofgroup( hecMESH, ndof, residual_vec,      vmax(1), vmax(2), loc(:,1), loc(:,2) )
    call fstr_get_maxabs_dofgroup( hecMESH, ndof, fstrSOLID%QFORCE,  vmax(3), vmax(4) )
    call fstr_get_maxabs_dofgroup( hecMESH, ndof, fstrSOLID%DFORCE,  vmax(5), vmax(6) )
    if( has_dx ) then
      call fstr_get_maxabs_dofgroup( hecMESH, ndof, hecMAT%X,         vmax(7), vmax(8), loc(:,3), loc(:,4) )
      call fstr_get_maxabs_dofgroup( hecMESH, ndof, fstrSOLID%dunode, vmax(9), vmax(10) )
      do i = 1, num_lagrange
        vmax(11) = max( vmax(11), abs(hecMAT%X(npndof+i)) )
        vmax(12) = max( vmax(12), abs(hecLagMAT%Lagrange(i)) )
      enddo
    endif
    gmax(:) = vmax(:)
    call hecmw_allreduce_R( hecMESH, gmax, NMAX, hecmw_max )

    do k = 1, NLOC
      if( gmax(IMAX_LOC(k)) > vmax(IMAX_LOC(k)) ) loc(1,k) = 0
    enddo
    gnode(:) = loc(1,:)
    call hecmw_allreduce_I( hecMESH, gnode, NLOC, hecmw_max )
    do k = 1, NLOC
      gdof(k) = 0
      if( loc(1,k) == gnode(k) ) gdof(k) = loc(2,k)
    enddo
    call hecmw_allreduce_I( hecMESH, gdof, NLOC, hecmw_max )

    ! --- residual relative to the reference force ---
    fref_t = max( gmax(3), gmax(5) )
    call fstr_set_relative_measure( cnvstat%m(kcnvResidual,kcnvTranslation,kcnvMax), &
        gmax(1), fref_t )
    if( ndof == 6 ) then
      call fstr_set_relative_measure( cnvstat%m(kcnvResidual,kcnvRotation,kcnvMax), &
          gmax(2), max( gmax(4), gmax(6) ) )
    endif

    ! --- correction relative to the accumulated increment, with the same rules as the L2 norms ---
    if( has_dx ) then
      dx_t = 1.0d0
      if( gmax(9) > 0.0d0 ) dx_t = gmax(7) / gmax(9)
      dx_r = 1.0d0
      if( gmax(10) > 0.0d0 ) dx_r = gmax(8) / gmax(10)
      cnvstat%m(kcnvCorrection,kcnvTranslation,kcnvMax)%value = dx_t
      if( ndof == 6 ) cnvstat%m(kcnvCorrection,kcnvRotation,kcnvMax)%value = dx_r
      call fstr_set_relative_measure( cnvstat%m(kcnvCorrection,kcnvLagrange,kcnvMax), &
          gmax(11), max( gmax(12), fref_t ) )
    endif

    do k = 1, NLOC
      cnvstat%m(IQ_LOC(k),IG_LOC(k),kcnvMax)%node = gnode(k)
      cnvstat%m(IQ_LOC(k),IG_LOC(k),kcnvMax)%dof = gdof(k)
    enddo

  end subroutine fstr_evaluate_convergence_max

  !> \brief A norm relative to its reference, or the norm itself when the reference is below FREF_FLOOR.
  subroutine fstr_set_relative_measure( ms, val, fref )
    implicit none
    type(fstr_convergence_measure), intent(inout) :: ms
    real(kind=kreal), intent(in)                  :: val   !< norm to measure
    real(kind=kreal), intent(in)                  :: fref  !< reference norm

    ms%absolute = ( fref < FREF_FLOOR )
    if( ms%absolute ) then
      ms%value = val
    else
      ms%value = val / fref
    endif

  end subroutine fstr_set_relative_measure

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

  !> \brief Largest absolute component of a nodal vector over the internal nodes, split into translational (DOF 1-3)
  !>        and rotational (DOF 4-6) components, with the global node ID and the DOF where it occurs.
  !>
  !> As in fstr_get_sqnorm_dofgroup, the result is not reduced across subdomains.
  subroutine fstr_get_maxabs_dofgroup( hecMESH, ndof, vec, vmax_t, vmax_r, loc_t, loc_r )
    implicit none
    type(hecmwST_local_mesh), intent(in)      :: hecMESH
    integer(kind=kint), intent(in)            :: ndof
    real(kind=kreal), intent(in)              :: vec(:)
    real(kind=kreal), intent(out)             :: vmax_t
    real(kind=kreal), intent(out)             :: vmax_r
    integer(kind=kint), intent(out), optional :: loc_t(2)  !< (global node ID, DOF) of vmax_t, 0 if vec vanishes
    integer(kind=kint), intent(out), optional :: loc_r(2)  !< (global node ID, DOF) of vmax_r, 0 if vec vanishes

    integer(kind=kint) :: i, idof, idx, gid
    integer(kind=kint) :: lt(2), lr(2)
    real(kind=kreal)   :: a

    vmax_t = 0.0d0
    vmax_r = 0.0d0
    lt(:) = 0
    lr(:) = 0
    ! Ties between nodes go to the largest node ID, as in the reduction across subdomains, so that the location does
    ! not depend on the partitioning.
    do i = 1, hecMESH%nn_internal
      idx = ndof*(i-1)
      gid = hecMESH%global_node_ID(i)
      do idof = 1, min(ndof,3)
        a = abs(vec(idx+idof))
        if( a > vmax_t .or. ( a >= vmax_t .and. a > 0.0d0 .and. gid > lt(1) ) ) then
          vmax_t = a
          lt(1) = gid
          lt(2) = idof
        endif
      enddo
      if( ndof /= 6 ) cycle
      do idof = 4, 6
        a = abs(vec(idx+idof))
        if( a > vmax_r .or. ( a >= vmax_r .and. a > 0.0d0 .and. gid > lr(1) ) ) then
          vmax_r = a
          lr(1) = gid
          lr(2) = idof
        endif
      enddo
    enddo
    if( present(loc_t) ) loc_t(:) = lt(:)
    if( present(loc_r) ) loc_r(:) = lr(:)

  end subroutine fstr_get_maxabs_dofgroup

  !> \brief Print the checked criteria of one Newton iteration.
  !>
  !> The L2 norms make up the line of the iteration; the max norms, if any is checked, follow on a line of their own
  !> aligned with the first criterion.
  subroutine fstr_print_convergence_state( iter, cnvstat )
    implicit none
    integer(kind=kint), intent(in)           :: iter
    type(fstr_convergence_state), intent(in) :: cnvstat

    character(len=512) :: line
    integer(kind=kint) :: iq, ig

    write(line,'(a,i8)') " iter:", iter
    do iq = kcnvResidual, kcnvCorrection
      do ig = kcnvTranslation, kcnvLagrange
        if( .not. cnvstat%m(iq,ig,kcnvL2)%check ) cycle
        write(line(len_trim(line)+1:),'(a,a,a,1pe11.4)') ", ", &
            trim(fstr_convergence_label( iq, ig, kcnvL2, cnvstat%m(iq,ig,kcnvL2)%absolute )), ":", &
            cnvstat%m(iq,ig,kcnvL2)%value
      enddo
    enddo
    write(*,'(a)') trim(line)

    line = ""
    do iq = kcnvResidual, kcnvCorrection
      do ig = kcnvTranslation, kcnvLagrange
        if( .not. cnvstat%m(iq,ig,kcnvMax)%check ) cycle
        if( len_trim(line) > 0 ) line = trim(line) // ","
        write(line(len_trim(line)+1:),'(a,a,a,1pe11.4)') " ", &
            trim(fstr_convergence_label( iq, ig, kcnvMax, cnvstat%m(iq,ig,kcnvMax)%absolute )), ":", &
            cnvstat%m(iq,ig,kcnvMax)%value
        if( cnvstat%m(iq,ig,kcnvMax)%node > 0 ) then
          write(line(len_trim(line)+1:),'(a,i0,a,i0,a)') " (node ", cnvstat%m(iq,ig,kcnvMax)%node, &
              ", dof ", cnvstat%m(iq,ig,kcnvMax)%dof, ")"
        endif
      enddo
    enddo
    if( len_trim(line) > 0 ) write(*,'(a,a)') repeat(" ", 15), trim(line)

  end subroutine fstr_print_convergence_state

  !> \brief Label of a criterion in the iteration log.
  function fstr_convergence_label( iq, ig, inorm, absolute ) result( label )
    implicit none
    integer(kind=kint), intent(in) :: iq        !< quantity
    integer(kind=kint), intent(in) :: ig        !< DOF group
    integer(kind=kint), intent(in) :: inorm     !< norm
    logical, intent(in)            :: absolute  !< the value is not normalized
    character(len=32)              :: label

    character(len=5), parameter  :: RES_NAME(3)  = [ 'force', 'mom  ', 'lag  ' ]
    character(len=10), parameter :: CORR_NAME(3) = [ 'disp.corr.', 'rot.corr. ', 'lag.corr. ' ]

    if( iq == kcnvResidual ) then
      label = 'res(' // trim(RES_NAME(ig))
      if( inorm == kcnvMax ) label = trim(label) // ',max'
      if( absolute ) label = trim(label) // ',abs'
      label = trim(label) // ')'
    else
      label = CORR_NAME(ig)
      if( inorm == kcnvMax ) label = trim(label) // '(max)'
    endif

  end function fstr_convergence_label

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

!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C***
!C*** module hecmw_solver_CR
!C***
!
module hecmw_solver_CR

  public :: hecmw_solve_CR

contains
  !C
  !C*** hecmw_solve_CR
  !C
  subroutine hecmw_solve_CR(hecMESH, hecMAT, ITER, RESID, error, Tset, Tsol, Tcomm)
    use hecmw_util
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_matrix_misc
    use hecmw_solver_misc
    use hecmw_solver_las
    use hecmw_solver_scaling
    use hecmw_precond
    use hecmw_jad_type

    implicit none

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint ), intent(inout):: ITER, error
    real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

    integer(kind=kint ) :: N, NP, NDOF, NNDOF
    integer(kind=kint ) :: my_rank
    integer(kind=kint ) :: ITERlog, TIMElog
    real(kind=kreal), pointer :: B(:), X(:)

    real(kind=kreal), dimension(:,:), allocatable :: WW

    integer(kind=kint ) :: MAXIT

    ! local variables
    real   (kind=kreal):: TOL
    integer(kind=kint )::i
    real   (kind=kreal)::S_TIME, S1_TIME, E_TIME, E1_TIME, START_TIME, END_TIME
    real   (kind=kreal)::BNRM2, rTAr, rTAr_old, ApTAp
    real   (kind=kreal)::ALPHA, BETA, DNRM2, DNRM2_TRUE
    real   (kind=kreal)::RESID_REC, RESID_TRUE
    real   (kind=kreal)::t_max, t_min, t_avg, t_sd
    logical :: true_checked

    ! PCR vector definitions
    integer(kind=kint), parameter :: R = 1    ! {z} = [Minv]{rt}
    integer(kind=kint), parameter :: P = 2    ! search direction {p}
    integer(kind=kint), parameter :: AP = 3   ! [A]{p}
    integer(kind=kint), parameter :: AZ = 4   ! [A]{z}
    integer(kind=kint), parameter :: WK = 4   ! preconditioner / true-residual workspace; shares AZ column (lifetimes disjoint)
    integer(kind=kint), parameter :: MAP = 5  ! [Minv][A]{p}
    integer(kind=kint), parameter :: RT = 6   ! recursive residual {rt} = {b} - [A]{x}

    integer(kind=kint), parameter :: N_ITER_CHECK_TRUE_R = 100

    call hecmw_barrier(hecMESH)
    S_TIME = HECMW_WTIME()

    !C===
    !C +-------+
    !C | INIT. |
    !C +-------+
    !C===
    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
    NNDOF = N * NDOF
    my_rank = hecMESH%my_rank
    X => hecMAT%X
    B => hecMAT%B

    ITERlog = hecmw_mat_get_iterlog( hecMAT )
    TIMElog = hecmw_mat_get_timelog( hecMAT )
    MAXIT   = hecmw_mat_get_iter( hecMAT )
    TOL     = hecmw_mat_get_resid( hecMAT )

    error = 0
    rTAr_old = 0.0d0

    allocate (WW(NDOF*NP, 6))
    WW = 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_INIT(hecMAT)
    endif

    !C===
    !C +----------------------+
    !C | SETUP PRECONDITIONER |
    !C +----------------------+
    !C===
    call hecmw_precond_setup(hecMAT, hecMESH, 1)

    !C===
    !C +-----------------------+
    !C | {rt}= {b} - [A]{x}    |
    !C +-----------------------+
    !C===
    call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,RT), Tcomm)

    !C-- compute ||{b}||
    call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
    if (BNRM2.eq.0.d0) then
      ITER = 0
      MAXIT = 0
      RESID = 0.d0
      X = 0.d0
    else
      !C-- Check the true initial residual before constructing PCR vectors
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,RT), WW(:,RT), DNRM2, Tcomm)
      RESID = dsqrt(DNRM2/BNRM2)
      if (RESID.le.TOL) then
        ITER = 0
        MAXIT = 0
      else
        !C-- {z}= [Minv]{rt}
        call hecmw_precond_apply(hecMESH, hecMAT, WW(:,RT), WW(:,R), WW(:,WK), Tcomm)

        !C-- {p}= {z}
        call hecmw_copy_R(NNDOF, WW(:,R), WW(:,P))

        !C-- {az}= [A]{z}, {ap}= {az}
        call hecmw_matvec(hecMESH, hecMAT, WW(:,R), WW(:,AZ), Tcomm)
        call hecmw_copy_R(NNDOF, WW(:,AZ), WW(:,AP))

        !C-- rTAr_old= {z}^T[A]{z}
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,AZ), rTAr_old, Tcomm)
        if (rTAr_old /= rTAr_old) then
          error = HECMW_SOLVER_ERROR_DIVERGE_NAN
          ITER = 0
          MAXIT = 0
        elseif (rTAr_old <= 0.0d0) then
          error = HECMW_SOLVER_ERROR_DIVERGE_MAT
          ITER = 0
          MAXIT = 0
        endif
      endif
    endif

    E_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E_TIME - S_TIME, t_max, t_min, t_avg, t_sd)
      if (hecMESH%my_rank.eq.0) then
        write(*,*) 'Time solver setup'
        write(*,*) '  Max     :',t_max
        write(*,*) '  Min     :',t_min
        write(*,*) '  Avg     :',t_avg
        write(*,*) '  Std Dev :',t_sd
      endif
      Tset = t_max
    else
      Tset = E_TIME - S_TIME
    endif

    Tcomm = 0.d0
    call hecmw_barrier(hecMESH)
    S1_TIME = HECMW_WTIME()

    !C
    !C*************************************************************** Iteration begins
    !C
    do iter = 1, MAXIT

      !C===
      !C +-------------------------+
      !C | {map}= [Minv]{ap}    |
      !C +-------------------------+
      !C===
      ! AZ and WK share one column because [A]{z} is dead while the preconditioner workspace is needed.
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,AP), WW(:,MAP), WW(:,WK), Tcomm)

      !C===
      !C +--------------------------------+
      !C | ApTAp= {ap}^T[Minv]{ap}    |
      !C +--------------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,AP), WW(:,MAP), ApTAp, Tcomm)

      if (ApTAp /= ApTAp) then
        error = HECMW_SOLVER_ERROR_DIVERGE_NAN
        exit
      elseif (ApTAp <= 0.0d0) then
        ! A symmetric positive definite preconditioner must give {ap}^T[Minv]{ap} > 0 for nonzero {ap}.
        error = HECMW_SOLVER_ERROR_DIVERGE_PC
        exit
      endif
      ALPHA = rTAr_old / ApTAp

      !C===
      !C +--------------------------------+
      !C | {x}= {x} + ALPHA*{p}        |
      !C | {rt}= {rt} - ALPHA*{ap}     |
      !C | {z}= {z} - ALPHA*{map}      |
      !C +--------------------------------+
      !C===
      call hecmw_axpy_R(NNDOF,  ALPHA, WW(:,P), X)
      call hecmw_axpy_R(NNDOF, -ALPHA, WW(:,AP), WW(:,RT))
      call hecmw_axpy_R(NNDOF, -ALPHA, WW(:,MAP), WW(:,R))

      !C===
      !C +------------------------------------+
      !C | RESID_REC= ||{rt}|| / ||{b}||     |
      !C +------------------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,RT), WW(:,RT), DNRM2, Tcomm)
      RESID_REC = dsqrt(DNRM2/BNRM2)

      !C===
      !C +------------------------------------------------+
      !C | Explicit residual RESID_TRUE for monitoring   |
      !C +------------------------------------------------+
      !C===
      ! Compute the explicit residual only for convergence monitoring.
      ! Do not replace the recursive PCR residual or restart the search directions.
      ! It is checked when the recursive residual reaches the tolerance, on a
      ! periodic monitoring iteration, or at the maximum iteration count.
      true_checked = .false.
      if (RESID_REC.le.TOL .or. mod(ITER,N_ITER_CHECK_TRUE_R)==0 .or. ITER.eq.MAXIT) then
        ! WK reuses the AZ column, which is free until [A]{z} is formed below.
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,WK), Tcomm)
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,WK), WW(:,WK), DNRM2_TRUE, Tcomm)
        RESID_TRUE = dsqrt(DNRM2_TRUE/BNRM2)
        true_checked = .true.
      endif

      ! Return/display value: the explicit residual on a monitoring iteration,
      ! the recursive residual otherwise.
      if (true_checked) then
        RESID = RESID_TRUE
      else
        RESID = RESID_REC
      endif

      !C##### Iteration history output
      if (my_rank.eq.0.and.ITERlog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID
      !C#####

      if (true_checked) then
        if (RESID_TRUE.le.TOL) exit
        if (ITER.eq.MAXIT) then
          error = HECMW_SOLVER_ERROR_NOCONV_MAXIT
          exit
        endif
      endif

      !C===
      !C +--------------------+
      !C | {az}= [A]{z}      |
      !C +--------------------+
      !C===
      call hecmw_matvec(hecMESH, hecMAT, WW(:,R), WW(:,AZ), Tcomm)

      !C===
      !C +-----------------------------+
      !C | rTAr= {z}^T[A]{z}          |
      !C | BETA= rTAr / rTAr_old      |
      !C +-----------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,AZ), rTAr, Tcomm)

      if (rTAr /= rTAr) then
        error = HECMW_SOLVER_ERROR_DIVERGE_NAN
        exit
      elseif (rTAr <= 0.0d0) then
        error = HECMW_SOLVER_ERROR_DIVERGE_MAT
        exit
      endif

      BETA = rTAr / rTAr_old
      rTAr_old = rTAr

      !C===
      !C +------------------------------+
      !C | {p}= {z} + BETA*{p}        |
      !C | {ap}= {az} + BETA*{ap}     |
      !C +------------------------------+
      !C===
      call hecmw_xpay_R(NNDOF, BETA, WW(:,R), WW(:,P))
      call hecmw_xpay_R(NNDOF, BETA, WW(:,AZ), WW(:,AP))

    enddo
    if (MAXIT.eq.0) ITER = 0
    !C
    !C*************************************************************** Iteration ends
    !C

    call hecmw_solver_scaling_bk(hecMAT)
    !C
    !C-- INTERFACE data EXCHANGE
    !C
    START_TIME = HECMW_WTIME()
    call hecmw_update_R(hecMESH, X, hecMAT%NP, hecMAT%NDOF)
    END_TIME = HECMW_WTIME()
    Tcomm = Tcomm + END_TIME - START_TIME

    deallocate(WW)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_FINALIZE(hecMAT)
    endif

    E1_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E1_TIME - S1_TIME, t_max, t_min, t_avg, t_sd)
      if (hecMESH%my_rank.eq.0) then
        write(*,*) 'Time solver iterations'
        write(*,*) '  Max     :',t_max
        write(*,*) '  Min     :',t_min
        write(*,*) '  Avg     :',t_avg
        write(*,*) '  Std Dev :',t_sd
      endif
      Tsol = t_max
    else
      Tsol = E1_TIME - S1_TIME
    endif

  end subroutine hecmw_solve_CR
end module hecmw_solver_CR

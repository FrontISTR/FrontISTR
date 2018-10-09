!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_CG
!C***
!C
module hecmw_solver_CG

  public :: hecmw_solve_CG

contains
  !C
  !C*** CG_nn
  !C
  subroutine hecmw_solve_CG( hecMESH,  hecMAT, ITER, RESID, error, &
      &                              Tset, Tsol, Tcomm )

    use hecmw_util
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_matrix_misc
    use hecmw_solver_misc
    use hecmw_solver_las
    use hecmw_solver_scaling
    use hecmw_precond
    use hecmw_jad_type
    use hecmw_estimate_condition

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

    integer(kind=kint), parameter ::  R= 1
    integer(kind=kint), parameter ::  Z= 2
    integer(kind=kint), parameter ::  Q= 2
    integer(kind=kint), parameter ::  P= 3
    integer(kind=kint), parameter :: WK= 4

    integer(kind=kint ) :: MAXIT

    ! local variables
    real   (kind=kreal) :: TOL
    integer(kind=kint )::i
    real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME, START_TIME, END_TIME
    real   (kind=kreal)::BNRM2
    real   (kind=kreal)::RHO,RHO1,BETA,C1,ALPHA,DNRM2
    real   (kind=kreal)::t_max,t_min,t_avg,t_sd
    integer(kind=kint) ::ESTCOND
    real   (kind=kreal), allocatable::D(:),E(:)
    real   (kind=kreal)::ALPHA1
    integer(kind=kint) :: n_indef_precond

    integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 50

    call hecmw_barrier(hecMESH)
    S_TIME= HECMW_WTIME()

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
    MAXIT  = hecmw_mat_get_iter( hecMAT )
    TOL   = hecmw_mat_get_resid( hecMAT )
    ESTCOND = hecmw_mat_get_estcond( hecMAT )

    error = 0
    n_indef_precond = 0
    RHO1 = 0.0d0
    ALPHA1 = 0.0d0
    BETA = 0.0d0

    allocate (WW(NDOF*NP, 4))
    WW = 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_INIT(hecMAT)
    endif

    if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
      allocate(D(MAXIT),E(MAXIT-1))
    endif
    !C===
    !C +----------------------+
    !C | SETUP PRECONDITIONER |
    !C +----------------------+
    !C===
    call hecmw_precond_setup(hecMAT, hecMESH, 1)

    !C===
    !C +---------------------+
    !C | {r0}= {b} - [A]{x0} |
    !C +---------------------+
    !C===
    call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)

    !C-- compute ||{b}||
    call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
    if (BNRM2.eq.0.d0) then
      iter = 0
      MAXIT = 0
      RESID = 0.d0
      X = 0.d0
    endif

    E_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E_TIME - S_TIME, &
        t_max, t_min, t_avg, t_sd)
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
    !C************************************************* Conjugate Gradient Iteration start
    !C
    do iter = 1, MAXIT

      !C===
      !C +----------------+
      !C | {z}= [Minv]{r} |
      !C +----------------+
      !C===
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,R), WW(:,Z), WW(:,WK), Tcomm)


      !C===
      !C +---------------+
      !C | {RHO}= {r}{z} |
      !C +---------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,Z), RHO, Tcomm)
      ! if RHO is NaN or Inf then no converge
      if (RHO == 0.d0) then
        ! converged due to RHO==0
        exit
      elseif (iter > 1 .and. RHO*RHO1 <= 0) then
        n_indef_precond = n_indef_precond + 1
        if (n_indef_precond >= 3) then
          ! diverged due to indefinite preconditioner
          error = HECMW_SOLVER_ERROR_DIVERGE_PC
          exit
        endif
      endif

      !C===
      !C +-----------------------------+
      !C | {p} = {z} if      ITER=1    |
      !C | BETA= RHO / RHO1  otherwise |
      !C +-----------------------------+
      !C===
      if ( ITER.eq.1 ) then
        do i = 1, NNDOF
          WW(i,P) = WW(i,Z)
        enddo
      else
        BETA = RHO / RHO1
        do i = 1, NNDOF
          WW(i,P) = WW(i,Z) + BETA*WW(i,P)
        enddo
      endif

      !C===
      !C +--------------+
      !C | {q}= [A] {p} |
      !C +--------------+
      !C===
      call hecmw_matvec(hecMESH, hecMAT, WW(:,P), WW(:,Q), Tcomm)

      !C===
      !C +---------------------+
      !C | ALPHA= RHO / {p}{q} |
      !C +---------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,P), WW(:,Q), C1, Tcomm)
      ! if C1 is NaN or Inf, no converge
      if (C1 <= 0) then
        ! diverged due to indefinite or negative definite matrix
        error = HECMW_SOLVER_ERROR_DIVERGE_MAT
        exit
      endif

      ALPHA= RHO / C1

      !C===
      !C +----------------------+
      !C | {x}= {x} + ALPHA*{p} |
      !C | {r}= {r} - ALPHA*{q} |
      !C +----------------------+
      !C===
      do i = 1, NNDOF
        X(i)  = X(i)    + ALPHA * WW(i,P)
        ! WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo

      if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
      else
        do i = 1, NNDOF
          WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
        enddo
      endif

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)

      RESID= dsqrt(DNRM2/BNRM2)

      !C##### ITERATION HISTORY
      if (my_rank.eq.0.and.ITERLog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID
      !C#####

      if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
        if (ITER == 1) then
          D(1) = 1.d0 / ALPHA
        else
          D(ITER) = 1.d0 / ALPHA + BETA / ALPHA1
          E(ITER-1) = sqrt(BETA) / ALPHA1
        endif
        ALPHA1 = ALPHA
        if (mod(ITER,ESTCOND) == 0) call hecmw_estimate_condition_CG(ITER, D, E)
      endif

      if ( RESID.le.TOL   ) then
        if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) exit
        !C----- recompute R to make sure it is really converged
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
        RESID= dsqrt(DNRM2/BNRM2)
        if ( RESID.le.TOL ) exit
      endif
      if ( ITER .eq.MAXIT ) error = HECMW_SOLVER_ERROR_NOCONV_MAXIT

      RHO1 = RHO

    enddo
    !C
    !C************************************************* Conjugate Gradient Iteration end
    !C
    call hecmw_solver_scaling_bk(hecMAT)
    !C
    !C-- INTERFACE data EXCHANGE
    !C
    START_TIME= HECMW_WTIME()
    call hecmw_update_m_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
    END_TIME = HECMW_WTIME()
    Tcomm = Tcomm + END_TIME - START_TIME

    deallocate (WW)
    !call hecmw_precond_clear(hecMAT)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_FINALIZE(hecMAT)
    endif

    if (ESTCOND /= 0 .and. error == 0 .and. hecMESH%my_rank == 0) then
      call hecmw_estimate_condition_CG(ITER, D, E)
      deallocate(D, E)
    endif

    E1_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E1_TIME - S1_TIME, &
        t_max, t_min, t_avg, t_sd)
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

  end subroutine hecmw_solve_CG

end module     hecmw_solver_CG

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
  subroutine hecmw_solve_CR( hecMESH,  hecMAT, ITER, RESID, error, &
      &                                    Tset, Tsol, Tcomm )
    use hecmw_util
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_matrix_misc
    use hecmw_solver_misc
    use hecmw_solver_las
    use hecmw_solver_scaling
    use hecmw_precond

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
    real(kind=kreal), dimension(2) :: CG

    integer(kind=kint ) :: MAXIT

    ! local variables
    real   (kind=kreal):: TOL
    integer(kind=kint )::i
    real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME, START_TIME, END_TIME
    real   (kind=kreal)::BNRM2, rTAr, rTAr_old, ApTAp
    real   (kind=kreal)::ALPHA, BETA, DNRM2
    real   (kind=kreal)::t_max,t_min,t_avg,t_sd
    integer(kind=kint) ::n_indef_precond

    ! PCR vector definitions (corrected version)
    integer(kind=kint), parameter :: R = 1     ! Residual vector r = M^{-1}(b - Ax)
    integer(kind=kint), parameter :: P = 2     ! Search direction vector p
    integer(kind=kint), parameter :: AP = 3    ! A*p
    integer(kind=kint), parameter :: AR = 4    ! A*r
    integer(kind=kint), parameter :: MAp = 5   ! M^{-1}*A*p
    integer(kind=kint), parameter :: WK = 6    ! Work vector

    integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 100

    call hecmw_barrier(hecMESH)
    S_time= HECMW_WTIME()

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
    n_indef_precond = 0

    allocate (WW(NDOF*NP, 6))
    WW = 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    !C===
    !C +----------------------+
    !C | SETUP PRECONDITIONER |
    !C +----------------------+
    !C===
    call hecmw_precond_setup(hecMAT, hecMESH, 1)

    !C===
    !C +---------------------+
    !C | {r0}= M^{-1}({b} - [A]{x0}) |
    !C +---------------------+
    !C===
    
    ! Use WK as temporary vector
    call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,WK), Tcomm)
    call hecmw_precond_apply(hecMESH, hecMAT, WW(:,WK), WW(:,R), WW(:,MAp), Tcomm)

    !C-- compute ||{b}||
    call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
    if (BNRM2.eq.0.d0) then
      iter = 0
      MAXIT = 0
      RESID = 0.d0
      X = 0.d0
    endif

    !C-- Initialize: p_0 = r_0
    do i=1, NNDOF
      WW(i,P) = WW(i,R)
    enddo

    !C-- Compute A*r_0 (AR = A*r)
    call hecmw_matvec(hecMESH, hecMAT, WW(:,R), WW(:,AR), Tcomm)

    !C-- Compute A*p_0 (= A*r_0)
    do i=1, NNDOF
      WW(i,AP) = WW(i,AR)
    enddo

    !C-- Initial rTAr calculation
    call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,AR), rTAr_old, Tcomm)
    if (rTAr_old <= 0.0d0) then
      ! Initial residual is 0 or indefinite matrix
      if (my_rank.eq.0) then
        write(*,*) 'PCR warning: Initial rTAr <= 0, possible indefinite matrix'
      endif
      ! Set small positive value to continue
      rTAr_old = 1.0d-10
    endif

    E_time = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E_time - S_time, &
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
      Tset = E_time - S_time
    endif

    Tcomm = 0.d0
    call hecmw_barrier(hecMESH)
    S1_time = HECMW_WTIME()
    
    !C
    !C*************************************************************** Iteration begins
    !C
    do iter = 1, MAXIT

      !C===
      !C +--------------------------------+
      !C | rTAr = r_k^T * (A*r_k)        |
      !C +--------------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,AR), rTAr, Tcomm)

      !C===
      !C +------------------------------------+
      !C | Compute M^{-1}*A*p_k               |
      !C +------------------------------------+
      !C===
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,AP), WW(:,MAp), WW(:,WK), Tcomm)

      !C===
      !C +-------------------------------------------+
      !C | ApTAp = (A*p_k)^T * M^{-1} * (A*p_k)     |
      !C +-------------------------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,AP), WW(:,MAp), ApTAp, Tcomm)

      !C===
      !C +------------------------+
      !C | ALPHA = rTAr / ApTAp   |
      !C +------------------------+
      !C===
      if (ApTAp <= 0.0d0) then
        ! Divergence or invalid matrix
        error = HECMW_SOLVER_ERROR_DIVERGE_MAT
        exit
      endif
      ALPHA = rTAr / ApTAp

      !C===
      !C +--------------------+
      !C | x_{k+1} = x_k + ALPHA*p_k |
      !C +--------------------+
      !C===
      do i = 1, NNDOF
        X(i) = X(i) + ALPHA * WW(i,P)
      enddo

      !C===
      !C +--------------------------------------------+
      !C | r_{k+1} = r_k - ALPHA * M^{-1} * A * p_k  |
      !C +--------------------------------------------+
      !C===
      do i = 1, NNDOF
        WW(i,R) = WW(i,R) - ALPHA * WW(i,MAp)
      enddo

      !C===
      !C +-------------------------+
      !C | Periodically recompute residual |
      !C +-------------------------+
      !C===
      if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
        ! r = M^{-1}(b - Ax)
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,WK), Tcomm)
        call hecmw_precond_apply(hecMESH, hecMAT, WW(:,WK), WW(:,R), WW(:,MAp), Tcomm)
      endif

      !C===
      !C +-------------------------+
      !C | Compute A*r_{k+1} |
      !C +-------------------------+
      !C===
      call hecmw_matvec(hecMESH, hecMAT, WW(:,R), WW(:,AR), Tcomm)

      !C===
      !C +------------------------------------------+
      !C | BETA = (r_{k+1}^T * A*r_{k+1}) / (r_k^T * A*r_k) |
      !C +------------------------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,AR), rTAr, Tcomm)
      
      if (rTAr <= 0.0d0) then
        ! rTAr is non-positive, possible numerical issue
        n_indef_precond = n_indef_precond + 1
        if (n_indef_precond >= 3) then
          error = HECMW_SOLVER_ERROR_DIVERGE_PC
          exit
        endif
        ! Set small positive value and try to continue
        rTAr = 1.0d-10
      endif
      
      BETA = rTAr / rTAr_old
      rTAr_old = rTAr

      !C===
      !C +--------------------------------+
      !C | p_{k+1} = r_{k+1} + BETA * p_k |
      !C +--------------------------------+
      !C===
      do i = 1, NNDOF
        WW(i,P) = WW(i,R) + BETA * WW(i,P)
      enddo

      !C===
      !C +--------------------------------+
      !C | A*p_{k+1} = A*r_{k+1} + BETA * A*p_k |
      !C +--------------------------------+
      !C===
      do i = 1, NNDOF
        WW(i,AP) = WW(i,AR) + BETA * WW(i,AP)
      enddo

      !C===
      !C +-------------------------+
      !C | Convergence check |
      !C +-------------------------+
      !C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)

      RESID = dsqrt(DNRM2/BNRM2)

      !C##### Iteration history output
      if (my_rank.eq.0.and.ITERlog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID
      !C#####

      if (RESID.le.TOL) then
        if (mod(ITER,N_ITER_RECOMPUTE_R)==0) exit
        !C----- Recompute residual to confirm convergence
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,WK), Tcomm)
        call hecmw_precond_apply(hecMESH, hecMAT, WW(:,WK), WW(:,R), WW(:,MAp), Tcomm)
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
        RESID = dsqrt(DNRM2/BNRM2)
        if (RESID.le.TOL) exit
      endif
      if (ITER.eq.MAXIT) error = HECMW_SOLVER_ERROR_NOCONV_MAXIT

    enddo
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

    E1_time = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E1_time - S1_time, &
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
      Tsol = E1_time - S1_time
    endif

  end subroutine hecmw_solve_CR
end module hecmw_solver_CR
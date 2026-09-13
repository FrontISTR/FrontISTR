!C
!C*** module hecmw_solver_GroppCG
!C
module hecmw_solver_GroppCG

  public :: hecmw_solve_GroppCG

contains
  !C
  !C*** GroppCG
  !C
  subroutine hecmw_solve_GroppCG( hecMESH, hecMAT, ITER, RESID, error, &
      &                          Tset, Tsol, Tcomm )

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
    integer(kind=kint), intent(inout) :: ITER, error
    real(kind=kreal), intent(inout) :: RESID, Tset, Tsol, Tcomm

    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: my_rank
    integer(kind=kint) :: ITERlog, TIMElog
    real(kind=kreal), pointer :: B(:), X(:)

    real(kind=kreal), dimension(:,:), allocatable :: WW

    integer(kind=kint), parameter :: R  = 1
    integer(kind=kint), parameter :: U  = 2
    integer(kind=kint), parameter :: V  = 3
    integer(kind=kint), parameter :: Q  = 4
    integer(kind=kint), parameter :: P  = 5
    integer(kind=kint), parameter :: S  = 6
    integer(kind=kint), parameter :: WK = 7

    integer(kind=kint) :: MAXIT

    real(kind=kreal) :: TOL
    integer(kind=kint) :: i
    real(kind=kreal) :: S_TIME, S1_TIME, E_TIME, E1_TIME
    real(kind=kreal) :: START_TIME, END_TIME
    real(kind=kreal) :: BNRM2, DNRM2
    real(kind=kreal) :: ALPHA, ALPHA1, BETA
    real(kind=kreal) :: GAMMA, GAMMA1, DELTA
    real(kind=kreal) :: CG(2)
    real(kind=kreal) :: t_max, t_min, t_avg, t_sd
    integer(kind=kint) :: ESTCOND
    real(kind=kreal), allocatable :: D(:), E(:)
    integer(kind=kint) :: n_indef_precond

    integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R = 50

    call hecmw_barrier(hecMESH)
    S_TIME = HECMW_WTIME()

    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
    NNDOF = N * NDOF
    my_rank = hecMESH%my_rank
    X => hecMAT%X
    B => hecMAT%B

    ITERlog = hecmw_mat_get_iterlog(hecMAT)
    TIMElog = hecmw_mat_get_timelog(hecMAT)
    MAXIT = hecmw_mat_get_iter(hecMAT)
    TOL = hecmw_mat_get_resid(hecMAT)
    ESTCOND = hecmw_mat_get_estcond(hecMAT)

    error = 0
    ITER = 0
    RESID = 0.0d0
    n_indef_precond = 0
    ALPHA1 = 0.0d0
    BETA = 0.0d0

    allocate(WW(NDOF*NP, 7))
    WW = 0.0d0

    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)
    call hecmw_mat_integrate(hecMAT)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_INIT(hecMAT)
    endif

    if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
      allocate(D(MAXIT), E(MAXIT-1))
    endif

    call hecmw_precond_setup(hecMAT, hecMESH, 1)

    call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
    call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
    if (BNRM2.eq.0.0d0) then
      MAXIT = 0
      X = 0.0d0
    else
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,R), WW(:,U), WW(:,WK), Tcomm)
      call hecmw_copy_R(NNDOF, WW(:,U), WW(:,P))
      call hecmw_matvec(hecMESH, hecMAT, WW(:,P), WW(:,S), Tcomm)

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,U), CG(1))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,R), CG(2))
      START_TIME = HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 2, HECMW_SUM)
      END_TIME = HECMW_WTIME()
      Tcomm = Tcomm + END_TIME-START_TIME
      GAMMA = CG(1)
      DNRM2 = CG(2)
      RESID = dsqrt(DNRM2/BNRM2)
      if (RESID.le.TOL) MAXIT = 0
    endif

    E_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E_TIME-S_TIME, t_max, t_min, t_avg, t_sd)
      if (hecMESH%my_rank.eq.0) then
        write(*,*) 'Time solver setup'
        write(*,*) '  Max     :', t_max
        write(*,*) '  Min     :', t_min
        write(*,*) '  Avg     :', t_avg
        write(*,*) '  Std Dev :', t_sd
      endif
      Tset = t_max
    else
      Tset = E_TIME-S_TIME
    endif

    Tcomm = 0.0d0
    call hecmw_barrier(hecMESH)
    S1_TIME = HECMW_WTIME()

    do i = 1, MAXIT
      if (GAMMA.eq.0.0d0) then
        ITER = i
        error = HECMW_SOLVER_ERROR_DIVERGE_PC
        exit
      elseif (GAMMA.ne.GAMMA) then
        ITER = i
        error = HECMW_SOLVER_ERROR_DIVERGE_NAN
        exit
      endif

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,P), WW(:,S), DELTA, Tcomm)
      if (DELTA.le.0.0d0) then
        ITER = i
        error = HECMW_SOLVER_ERROR_DIVERGE_MAT
        exit
      elseif (DELTA.ne.DELTA) then
        ITER = i
        error = HECMW_SOLVER_ERROR_DIVERGE_NAN
        exit
      endif

      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,S), WW(:,Q), WW(:,WK), Tcomm)

      ALPHA = GAMMA/DELTA

      if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
        if (i.eq.1) then
          D(1) = 1.0d0/ALPHA
        else
          D(i) = 1.0d0/ALPHA+BETA/ALPHA1
          E(i-1) = dsqrt(BETA)/ALPHA1
        endif
        if (mod(i,ESTCOND).eq.0) call hecmw_estimate_condition_CG(i, D, E)
      endif

      call hecmw_axpy_R(NNDOF, ALPHA, WW(:,P), X)

      if (mod(i,N_ITER_RECOMPUTE_R).eq.0) then
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
        call hecmw_precond_apply(hecMESH, hecMAT, WW(:,R), WW(:,U), WW(:,WK), Tcomm)
      else
        call hecmw_axpy_R(NNDOF, -ALPHA, WW(:,S), WW(:,R))
        call hecmw_axpy_R(NNDOF, -ALPHA, WW(:,Q), WW(:,U))
      endif

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,U), CG(1))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,R), CG(2))
      START_TIME = HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 2, HECMW_SUM)
      END_TIME = HECMW_WTIME()
      Tcomm = Tcomm + END_TIME-START_TIME
      GAMMA1 = CG(1)
      DNRM2 = CG(2)
      RESID = dsqrt(DNRM2/BNRM2)
      ITER = i

      if (RESID.le.TOL) then
        if (mod(i,N_ITER_RECOMPUTE_R).ne.0) then
          call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
          call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
          RESID = dsqrt(DNRM2/BNRM2)
        endif
        if (my_rank.eq.0 .and. ITERlog.eq.1) write(*,'(i7, 1pe16.6)') ITER, RESID
        if (RESID.le.TOL) exit

        call hecmw_precond_apply(hecMESH, hecMAT, WW(:,R), WW(:,U), WW(:,WK), Tcomm)
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,U), GAMMA1, Tcomm)
      else
        if (my_rank.eq.0 .and. ITERlog.eq.1) write(*,'(i7, 1pe16.6)') ITER, RESID
      endif

      if (i.eq.MAXIT) then
        if (mod(i,N_ITER_RECOMPUTE_R).ne.0) then
          call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
          call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
          RESID = dsqrt(DNRM2/BNRM2)
        endif
        if (RESID.gt.TOL) error = HECMW_SOLVER_ERROR_NOCONV_MAXIT
        exit
      endif

      if (GAMMA1.eq.0.0d0) then
        error = HECMW_SOLVER_ERROR_DIVERGE_PC
        exit
      elseif (GAMMA1.ne.GAMMA1) then
        error = HECMW_SOLVER_ERROR_DIVERGE_NAN
        exit
      elseif (GAMMA1*GAMMA.le.0.0d0) then
        n_indef_precond = n_indef_precond + 1
        if (n_indef_precond.ge.3) then
          error = HECMW_SOLVER_ERROR_DIVERGE_PC
          exit
        endif
      endif

      call hecmw_matvec(hecMESH, hecMAT, WW(:,U), WW(:,V), Tcomm)

      BETA = GAMMA1/GAMMA
      call hecmw_xpay_R(NNDOF, BETA, WW(:,U), WW(:,P))
      call hecmw_xpay_R(NNDOF, BETA, WW(:,V), WW(:,S))

      GAMMA = GAMMA1
      ALPHA1 = ALPHA
    enddo

    call hecmw_solver_scaling_bk(hecMAT)

    START_TIME = HECMW_WTIME()
    call hecmw_update_R(hecMESH, X, hecMAT%NP, hecMAT%NDOF)
    END_TIME = HECMW_WTIME()
    Tcomm = Tcomm + END_TIME-START_TIME

    deallocate(WW)

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
      call hecmw_JAD_FINALIZE(hecMAT)
    endif

    if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
      if (error.eq.0 .and. ITER.gt.0) call hecmw_estimate_condition_CG(ITER, D, E)
      deallocate(D, E)
    endif

    E1_TIME = HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E1_TIME-S1_TIME, t_max, t_min, t_avg, t_sd)
      if (hecMESH%my_rank.eq.0) then
        write(*,*) 'Time solver iterations'
        write(*,*) '  Max     :', t_max
        write(*,*) '  Min     :', t_min
        write(*,*) '  Avg     :', t_avg
        write(*,*) '  Std Dev :', t_sd
      endif
      Tsol = t_max
    else
      Tsol = E1_TIME-S1_TIME
    endif

  end subroutine hecmw_solve_GroppCG

end module hecmw_solver_GroppCG
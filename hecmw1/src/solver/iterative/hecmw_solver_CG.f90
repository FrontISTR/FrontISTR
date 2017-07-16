!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_CG
public :: hecmw_solve_CG

contains
subroutine hecmw_solve_CG( hecMESH, hecMAT, ITER, RESID, ERROR, Tset, Tsol, Tcomm )

  use hecmw_util
  use m_hecmw_solve_error
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
  integer(kind=kint ), intent(inout):: ITER, ERROR
  real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm
  integer(kind=kint ) :: N, NP, NDOF, NNDOF, my_rank, ITERlog, TIMElog, MAXIT, ESTCOND
  real   (kind=kreal) :: TOL, S1_TIME, S2_TIME
  real   (kind=kreal), pointer :: B(:), X(:)
  
  integer(kind=kint), parameter ::  iX=0, iR=1, iZ=2, iQ=2, iP=3, iWK=4
  integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 50
  real   (kind=kreal), allocatable,target :: WW(:,:)
  real   (kind=kreal), pointer  :: R(:), Z(:), Q(:), P(:), WK(:)
  real   (kind=kreal), allocatable :: D(:),E(:)
  integer(kind=kint ):: i, n_indef_precond
  real   (kind=kreal):: RHO, RHO1, ALPHA1, BETA, C1, ALPHA, DNRM2, BNRM2

  call hecmw_barrier(hecMESH)
  S1_TIME= HECMW_WTIME()

  N = hecMAT%N
  NP = hecMAT%NP
  NDOF = hecMAT%NDOF
  NNDOF = N * NDOF
  my_rank = hecMESH%my_rank

  ITERlog = hecmw_mat_get_iterlog( hecMAT )
  TIMElog = hecmw_mat_get_timelog( hecMAT )
  MAXIT   = hecmw_mat_get_iter   ( hecMAT )
  TOL     = hecmw_mat_get_resid  ( hecMAT )
  ESTCOND = hecmw_mat_get_estcond( hecMAT )

  ERROR = 0
  n_indef_precond = 0

  allocate (WW (NDOF*NP,0:4))
  WW=0.0d0
  B=>hecMAT%B
  X => WW(:,iX)
  R => WW(:,iR)
  Z => WW(:,iZ)
  Q => WW(:,iQ)
  P => WW(:,iP)
  WK=> WW(:,iWK)
  call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

  IF (hecmw_mat_get_usejad(hecMAT).ne.0) call hecmw_JAD_INIT(hecMAT)
  IF (ESTCOND /= 0 .and. hecMESH%my_rank == 0) allocate(D(MAXIT),E(MAXIT-1))

  call hecmw_precond_setup(hecMAT, hecMESH, 1)

  !C | {r0}= {b} - [A]{x0} |
  call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)

  !C-- compute ||{b}||
  call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
  if (BNRM2.eq.0.d0) then
    iter = 0
    MAXIT = 0
    RESID = 0.d0
    X = 0.d0
  endif

  Tset = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver setup'
    call hecmw_print_time_statistics(hecMESH, Tset)
  endif

  Tcomm = 0.d0
  call hecmw_barrier(hecMESH)
  S1_TIME = HECMW_WTIME()

  !C************************************************* Conjugate Gradient Iteration start
  do iter = 1, MAXIT

    !C | {z}= [Minv]{r} |
    call hecmw_precond_apply(hecMESH, hecMAT, R, Z, WK, Tcomm)

    !C | {RHO}= {r}{z} |
    call hecmw_InnerProduct_R(hecMESH, NDOF, R, Z, RHO, Tcomm)
    ! if RHO is NaN or Inf then no converge
    if (RHO == 0.d0) then
      ! converged due to RHO==0
      exit
    elseif (iter > 1 .and. RHO*RHO1 <= 0) then
      n_indef_precond = n_indef_precond + 1
      if (n_indef_precond >= 3) then
        ! diverged due to indefinite preconditioner
        ERROR = HECMW_SOLVER_ERROR_DIVERGE_PC
        exit
      endif
    endif

    !C | {p} = {z} if      ITER=1    |
    !C | BETA= RHO / RHO1  otherwise |
    if ( ITER.eq.1 ) then
      do i = 1, NNDOF
        WW(i,iP) = WW(i,iZ)
      enddo
     else
       BETA = RHO / RHO1
       do i = 1, NNDOF
         WW(i,iP) = WW(i,iZ) + BETA*WW(i,iP)
       enddo
    endif

    !C | {q}= [A] {p} |
    call hecmw_matvec(hecMESH, hecMAT, P, Q, Tcomm)

    !C | ALPHA= RHO / {p}{q} |
    call hecmw_InnerProduct_R(hecMESH, NDOF, P, Q, C1, Tcomm)
    ! if C1 is NaN or Inf, no converge
    if (C1 <= 0) then
      ! diverged due to indefinite or negative definite matrix
      ERROR = HECMW_SOLVER_ERROR_DIVERGE_MAT
      exit
    endif

    ALPHA= RHO / C1

    !C | {x}= {x} + ALPHA*{p} |
    !C | {r}= {r} - ALPHA*{q} |
    do i = 1, NNDOF
       WW(i,iX)  = WW(i,iX)    + ALPHA * WW(i,iP)
    enddo

    if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
      call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
    else
      do i = 1, NNDOF
        WW(i,iR)= WW(i,iR) - ALPHA * WW(i,iQ)
      enddo
    endif

    call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)

    RESID= dsqrt(DNRM2/BNRM2)

    !C##### ITERATION HISTORY
    if (my_rank.eq.0.and.ITERLog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID

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
      call hecmw_matresid(hecMESH, hecMAT, X, B, P, Tcomm)
      call hecmw_InnerProduct_R(hecMESH, NDOF, P, P, DNRM2, Tcomm)
      RESID= dsqrt(DNRM2/BNRM2)
      if ( RESID.le.TOL ) exit
    endif
    if ( ITER .eq.MAXIT ) ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT
    RHO1 = RHO

  enddo
  !C************************************************* Conjugate Gradient Iteration end

  call hecmw_solver_scaling_bk(hecMAT)

  !C-- INTERFACE data EXCHANGE
  S2_TIME= HECMW_WTIME() 
  call hecmw_update_m_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
  Tcomm = Tcomm + HECMW_WTIME() - S2_TIME
  
  do i = 1, NDOF*NP
     hecMAT%X(i)=WW(i,iX)
  end do 

  deallocate (WW)
  !call hecmw_precond_clear(hecMAT)

  IF (hecmw_mat_get_usejad(hecMAT).ne.0) THEN
    call hecmw_JAD_FINALIZE(hecMAT)
  ENDIF

  if (ESTCOND /= 0 .and. ERROR == 0 .and. hecMESH%my_rank == 0) then
    call hecmw_estimate_condition_CG(ITER, D, E)
    deallocate(D, E)
  endif

  Tsol = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver iterations'
    call hecmw_print_time_statistics(hecMESH, Tsol)
  endif
  
end subroutine hecmw_solve_CG
end module     hecmw_solver_CG
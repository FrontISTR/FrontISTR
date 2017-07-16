!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_BiCGSTAB
public :: hecmw_solve_BiCGSTAB

contains
subroutine hecmw_solve_BiCGSTAB( hecMESH,  hecMAT, ITER, RESID, ERROR, Tset, Tsol, Tcomm )

  use hecmw_util
  use m_hecmw_solve_error
  use hecmw_matrix_misc
  use hecmw_solver_misc
  use hecmw_solver_las
  use hecmw_solver_scaling
  use hecmw_precond
  use hecmw_estimate_condition


  implicit none
  type(hecmwST_local_mesh) :: hecMESH
  type(hecmwST_matrix) :: hecMAT
  integer(kind=kint ), intent(inout):: ITER, ERROR
  real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm
  integer(kind=kint ) :: N, NP, NDOF, NNDOF, my_rank, ITERlog, TIMElog, MAXIT, ESTCOND
  real   (kind=kreal) :: TOL, S1_TIME, S2_TIME
  real   (kind=kreal), pointer :: B(:), X(:)

  integer(kind=kint), parameter ::  iR=1, iRT=2, iP=3, iPT=4, iS=5, iST=1, iT=6, iV=7, iWK=8, iX=9
  integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 100
  real   (kind=kreal), allocatable,target :: WW(:,:)
  real   (kind=kreal), pointer  :: R(:), RT(:), P(:), PT(:), S(:), ST(:), T(:), V(:), WK(:) 
  real   (kind=kreal), dimension(2) :: CG
  integer(kind=kint )::i
  real   (kind=kreal)::RHO, RHO1, BETA, ALPHA, DNRM2, BNRM2, C2, OMEGA

  
  call hecmw_barrier(hecMESH)
  S1_TIME= HECMW_WTIME()

  N = hecMAT%N
  NP = hecMAT%NP
  NDOF = hecMAT%NDOF
  NNDOF = N * NDOF
  my_rank = hecMESH%my_rank
  X => hecMAT%X
  B => hecMAT%B

  ITERlog = hecmw_mat_get_iterlog( hecMAT )
  TIMElog = hecmw_mat_get_timelog( hecMAT )
  MAXIT   = hecmw_mat_get_iter   ( hecMAT )
  TOL     = hecmw_mat_get_resid  ( hecMAT )
  ESTCOND = hecmw_mat_get_estcond( hecMAT )

  ERROR = 0
  allocate (WW (NDOF*NP,8))
  WW=0.0d0
  R => WW(:,iR )
  RT=> WW(:,iRT)
  P => WW(:,iP )
  PT=> WW(:,iPT)
  S => WW(:,iS )
  ST=> WW(:,iST)  
  T => WW(:,iT )  
  V => WW(:,iV )  
  WK=> WW(:,iWK)  
  call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)
  call hecmw_precond_setup(hecMAT, hecMESH, 0)

  !C | {r0}= {b} - [A]{x0} |
  call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)

  !C-- set arbitrary {r_tld}
  do i=1, NNDOF
    WW(i,iRT) = WW(i,iR)
  enddo

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

  !C*************************************************************** iterative procedures start
  do iter = 1, MAXIT

    !C | RHO= {r}{r_tld} |
    call hecmw_InnerProduct_R(hecMESH, NDOF, R, RT, RHO, Tcomm)

    !C | BETA= (RHO/RHO1) * (ALPHA/OMEGA)       |
    !C | {p} = {r} + BETA * ( {p} - OMEGA*{v} ) |
    if ( iter.gt.1 ) then
      BETA = (RHO/RHO1) * (ALPHA/OMEGA)
      do i = 1, NNDOF
        WW(i,iP) = WW(i,iR) + BETA * (WW(i,iP) - OMEGA * WW(i,iV))
      enddo
     else
      do i = 1, NNDOF
        WW(i,iP) = WW(i,iR)
      enddo
    endif


    !C | {p_tld}= [Minv]{p} |
    call hecmw_precond_apply(hecMESH, hecMAT, P, PT, WK, Tcomm)

    !C | {v}= [A] {p_tld} |
    call hecmw_matvec(hecMESH, hecMAT, PT, V, Tcomm)

    !C-- calc. ALPHA
    call hecmw_InnerProduct_R(hecMESH, NDOF, RT, V, C2, Tcomm)
    ALPHA = RHO / C2

    !C-- {s}= {r} - ALPHA*{V}
    do i = 1, NNDOF 
      WW(i,iS) = WW(i,iR) - ALPHA * WW(i,iV)
    enddo

    !C | {s_tld}= [Minv]{s} |
    call hecmw_precond_apply(hecMESH, hecMAT, S, ST, WK, Tcomm)

    !C | {t}= [A] {s_tld} |
    call hecmw_matvec(hecMESH, hecMAT, ST, T, Tcomm)

    !C | OMEGA= ({t}{s}) / ({t}{t}) |
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, S, CG(1))
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, T, CG(2))
    S2_TIME = HECMW_WTIME()
    call hecmw_allreduce_R(hecMESH, CG, 2, HECMW_SUM)
    Tcomm = Tcomm + HECMW_WTIME() - S2_TIME

    OMEGA = CG(1) / CG(2)


    !C | update {x},{r} |
    do i = 1, NNDOF
      X(i) = X(i) + ALPHA * WW(i,iPT) + OMEGA * WW(i,iST)
    enddo

    !C--- recompute R sometimes
    if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
      call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
    else
      do i = 1, NNDOF
        WW(i,iR) = WW(i,iS) - OMEGA * WW(i,iT)
      enddo
    endif

    call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)

    RESID= dsqrt(DNRM2/BNRM2)

    !C##### ITERATION HISTORY
    if (my_rank.eq.0.and.ITERlog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID

    if ( RESID.le.TOL   ) then
      if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) exit
    !C----- recompute R to make sure it is really converged
      call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
      call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)
      RESID= dsqrt(DNRM2/BNRM2)
      if ( RESID.le.TOL ) exit
    endif
    if ( ITER .eq.MAXIT ) ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT

    RHO1 = RHO

  enddo
  !C*************************************************************** iterative procedures end

  call hecmw_solver_scaling_bk(hecMAT)

  !C-- INTERFACE data EXCHANGE
  S2_TIME = HECMW_WTIME()
  call hecmw_update_m_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
  Tcomm = Tcomm + HECMW_WTIME() - S2_TIME

  deallocate (WW)
  !call hecmw_precond_clear(hecMAT)

  Tsol = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver iterations'
    call hecmw_print_time_statistics(hecMESH, Tsol)
  endif

end subroutine hecmw_solve_BiCGSTAB
end module     hecmw_solver_BiCGSTAB

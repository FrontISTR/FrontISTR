!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_GPBiCG
public :: hecmw_solve_GPBiCG

contains
subroutine hecmw_solve_GPBiCG( hecMESH,  hecMAT, ITER, RESID, ERROR, Tset, Tsol, Tcomm )

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
  
  integer(kind=kint), parameter :: iR=1, iRT=2, iT=3, iTT=4, iT0=5, iP=6, iPT=7, iU=8, iW1=9 
  integer(kind=kint), parameter :: iY=10, iZ=11, iWK=12, iW2=13, iZQ=14
  integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 20
  real(kind=kreal), allocatable,target :: WW(:,:)
  real   (kind=kreal), pointer :: R(:), RT(:), T(:), TT(:), T0(:), P(:), PT(:), U(:), W1(:) 
  real   (kind=kreal), pointer :: Y(:), Z(:), WK(:), W2(:), ZQ(:)
  integer(kind=kint ) :: i,j
  real   (kind=kreal) :: RHO, RHO1, BETA, ALPHA, DNRM2, BNRM2, QSI, ETA, COEF1
  real   (kind=kreal), dimension(2) :: RR

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
  
  ERROR= 0
  allocate (WW(NDOF*NP,14));
  do j=1,14
    do i=1,NDOF*NP
      WW(i,j) = 0.0d0
    end do 
  end do 
   R => WW(:, iR);  RT => WW(:,iRT);
   T => WW(:, iT);  TT => WW(:,iTT);  T0 => WW(:,iT0);
   P => WW(:, iP);  PT => WW(:,iPT);   U => WW(:, iU);
   Y => WW(:, iY);   Z => WW(:, iZ);  ZQ => WW(:,iZQ);
  W1 => WW(:,iW1);  W2 => WW(:,iW2);  WK => WW(:,iWK);


  call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)
  call hecmw_precond_setup(hecMAT, hecMESH, 0)

  !C | {r}= {b} - [A]{xini} |
  call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)

  do i= 1, NNDOF
    WW(i,iRT)= WW(i,iR)
  enddo
  !C==
  call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)

  if (BNRM2.eq.0.d0) then
    iter = 0
    MAXIT = 0
    RESID = 0.d0
    X = 0.d0
  endif

  call hecmw_InnerProduct_R(hecMESH, NDOF, RT, R, RHO, Tcomm)

  Tset = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver setup'
    call hecmw_print_time_statistics(hecMESH, Tset)
  endif

  !C*************************************************************** ITERATIVE PROC.
  call hecmw_barrier(hecMESH)
  S1_TIME= HECMW_WTIME()
  do iter= 1, MAXIT
    !C | {r}= [Minv]{r} |
    do j= 1, NNDOF
      WW(j,iWK)=  WW(j,iR)
    enddo

    call hecmw_precond_apply(hecMESH, hecMAT, WK, R, ZQ, Tcomm)

    !C | {p} = {r} + BETA * ( {p} - {u} ) |
    if (iter.gt.1) then
      do j= 1, NNDOF
        WW(j,iP)= WW(j,iR) + BETA*( WW(j,iP)-WW(j,iU))
      enddo
     else
      do j= 1, NNDOF
        WW(j,iP)= WW(j,iR)
      enddo
    endif

    !C | ALPHA= {r_tld}{r}/{r_tld} A{p} |
    !C-- calc. {p_tld}= A{p}
    call hecmw_matvec(hecMESH, hecMAT, P, PT, Tcomm)
    !C-- calc. ALPHA
    call hecmw_InnerProduct_R(hecMESH, NDOF, RT, PT, RHO1, Tcomm)

    ALPHA= RHO / RHO1

    !C | {y}= {t} - {r} - ALPHA{w} + ALPHA{p_tld} |
    !C | {t}= {r}                  - ALPHA{p_tld} |
    do j= 1, NNDOF
      WW(j,iY)= WW(j,iT) - WW(j,iWK) + ALPHA*(-WW(j,iW1) + WW(j,iPT))
      WW(j,iT)= WW(j,iWK) - ALPHA*WW(j,iPT)
    enddo
    !C | {t_tld}= [A][Minv]{t} |
    !C-- calc. {t_tld} and {t0} by [M] inversion
    !C         {W2}   = [Minv]{p_tld}
    call hecmw_precond_apply(hecMESH, hecMAT, T, TT, ZQ, Tcomm)
    call hecmw_precond_apply(hecMESH, hecMAT, T0, W2, ZQ, Tcomm)
    do i= 1, NNDOF
      WW(i,iT0)= WW(i,iW2)
    enddo
    call hecmw_precond_apply(hecMESH, hecMAT, PT, W2, ZQ, Tcomm)

    !C-- calc. [A]{t_tld}
    call hecmw_matvec(hecMESH, hecMAT, TT, WK, Tcomm)

    do i= 1, NNDOF
      WW(i,iTT)= WW(i,iWK)
    enddo

    !C | calc. QSI and ETA |
    !call pol_coef(iter, T, TT, Y, QSI, ETA)
    !call pol_coef_vanilla(iter, T, TT, Y, QSI, ETA)
    call pol_coef_vanilla2(iter,NDOF, T, TT, Y, QSI, ETA)


    !C | {u} = QSI [Minv]{pt} + ETA([Minv]{t0}-[Minv]{r}+BETA*{u} |
    !C | {z} = QSI [Minv]{r}  + ETA*{z} - ALPHA*{u}               |
    !C-- compute. {u},{z}
    if (iter.gt.1) then
      do j= 1, NNDOF
        WW(j,iU)= QSI* W2(j) + ETA*(WW(j,iT0) - WW(j,iR) + BETA*WW(j,iU))
        WW(j,iZ)= QSI* WW(j,iR) + ETA* WW(j,iZ) -          ALPHA*WW(j,iU)
      enddo
    else
      do j= 1, NNDOF
        WW(j,iU)= QSI* W2(j) + ETA*(WW(j,iT0) - WW(j,iR))
        WW(j,iZ)= QSI* WW(j,iR) + ETA* WW(j,iZ) -          ALPHA*WW(j,iU)
      enddo
    endif

    !C | update {x},{r},{w} |
    do j= 1, NNDOF
       X(j)= X(j) + ALPHA*WW(j,iP) + WW(j,iZ)
      WW(j,iT0)= WW(j,iT)
    enddo

    !C--- recompute R sometimes
    if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
      call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
    else
      do j= 1, NNDOF
        WW(j,iR)= WW(j,iT) - ETA*WW(j,iY) - QSI*WW(j,iTT)
      enddo
    endif

    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, R,  R, RR(1))
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, R, RT, RR(2))
    S2_TIME= HECMW_WTIME()
    call hecmw_allreduce_R(hecMESH, RR, 2, HECMW_SUM)
    Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME
    DNRM2 = RR(1)
    COEF1 = RR(2)

    BETA = ALPHA*COEF1 / (QSI*RHO)
    do j= 1, NNDOF
      WW(j,iW1)= WW(j,iTT) + BETA*WW(j,iPT)
    enddo

    RESID= dsqrt(DNRM2/BNRM2)
    RHO  = COEF1

    !C##### ITERATION HISTORY
    if (my_rank.eq.0 .and. ITERlog.eq.1) write (*, '(i5, 1pe16.6)') ITER, RESID

    if (RESID.le.TOL   ) then
      if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) exit
    !C----- recompute R to make sure it is really converged
      call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
      call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)
      RESID= dsqrt(DNRM2/BNRM2)
      if (RESID.le.TOL) exit
    endif
    if ( ITER.eq.MAXIT ) ERROR= HECMW_SOLVER_ERROR_NOCONV_MAXIT
    !C===
  enddo

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
 
  contains
  subroutine pol_coef(iter, T, TT, Y, QSI, ETA)
  !C*** pol_coef : computes QSI and ETA in original GPBiCG way
    implicit none
    integer(kind=kint), intent(in) :: iter
    real(kind=kreal), intent(in) :: T(:), TT(:), Y(:)
    real(kind=kreal), intent(out) :: QSI, ETA

    real(kind=kreal), dimension(5) :: CG
    real(kind=kreal), dimension(2) :: EQ
    real(kind=kreal) :: delta

    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, Y, CG(1)) ! myu1
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, TT, T, CG(2)) ! omega2
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, T, CG(3)) ! omega1
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, TT, Y, CG(4)) ! nyu
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, TT, TT, CG(5)) ! myu2
    S2_TIME= HECMW_WTIME()
    call hecmw_allreduce_R(hecMESH, CG, 5, HECMW_SUM)
    Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME

    if (iter.eq.1) then
      EQ(1)= CG(2)/CG(5) ! omega2 / myu2
      EQ(2)= 0.d0
    else
      delta= (CG(5)*CG(1)-CG(4)*CG(4))         ! myu1*myu2 - nyu^2
      EQ(1)= (CG(1)*CG(2)-CG(3)*CG(4)) / delta ! (myu1*omega2-nyu*omega1)/delta
      EQ(2)= (CG(5)*CG(3)-CG(4)*CG(2)) / delta ! (myu2*omega1-nyu*omega2)/delta
    endif

    QSI= EQ(1)
    ETA= EQ(2)
  end subroutine pol_coef

  subroutine pol_coef_vanilla(iter, T, TT, Y, QSI, ETA)
  !C*** pol_coef_vanilla : computes QSI and ETA with vanilla strategy
  !C      see Fujino, Abe, Sugihara and Nakashima(2013) ISBN978-4-621-08741-1
    implicit none
    integer(kind=kint), intent(in) :: iter
    real(kind=kreal), intent(inout) :: T(:), TT(:), Y(:)
    real(kind=kreal), intent(out) :: QSI, ETA

    real(kind=kreal), dimension(3) :: CG
    real(kind=kreal) :: gamma1, gamma2
    real(kind=kreal) :: c, c_abs

    real(kind=kreal), parameter :: OMEGA = 0.707106781d0

    if (iter.eq.1) then
      gamma1 = 0.d0
      gamma2 = 0.d0
    else
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, Y, CG(1)) ! myu
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, TT, CG(2)) ! nyu
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, T, CG(3)) ! omega
      S2_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
      Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME
      gamma1 = CG(3)/CG(1) ! omega / myu
      gamma2 = CG(2)/CG(1) ! nyu / myu
    !!! COMMENTED OUT because no convergence obtained with following updates.
    !         do j= 1, NNDOF
    !           WW(j,iT)= WW(j,iT) - gamma1*WW(j,iY)
    !           WW(j,iTT)= WW(j,iTT) - gamma2*WW(j,iY)
    !         enddo
    endif

    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, T, CG(1)) ! |r|^2
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, TT, TT, CG(2)) ! |s|^2
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, TT, CG(3)) ! r.s
    S2_TIME= HECMW_WTIME()
    call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
    Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME

    c = CG(3) / dsqrt(CG(1)*CG(2))
    c_abs = dabs(c)
    if (c_abs > OMEGA) then
      QSI = c * dsqrt(CG(1)/CG(2))
    else
      ! QSI = (c / c_abs) * OMEGA * dsqrt(CG(1)/CG(2))
      if (c >= 0.d0) then
        QSI = OMEGA * dsqrt(CG(1)/CG(2))
      else
        QSI = -OMEGA * dsqrt(CG(1)/CG(2))
      endif
    endif
    ETA = gamma1 - QSI*gamma2
  end subroutine pol_coef_vanilla

  subroutine pol_coef_vanilla2(iter,NDOF, T, TT, Y, QSI, ETA)
  !*** pol_coef_vanilla2 : optimized version of pol_coef_vanilla
    implicit none
    integer(kind=kint), intent(in) :: iter,NDOF
    real(kind=kreal), intent(inout) :: T(:), TT(:), Y(:)
    real(kind=kreal), intent(out) :: QSI, ETA

    real(kind=kreal), dimension(6) :: CG
    real(kind=kreal) :: gamma1, gamma2
    real(kind=kreal) :: c, c_abs

    real(kind=kreal), parameter :: OMEGA = 0.707106781d0

    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, T, CG(1)) ! |r|^2
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, TT, TT, CG(2)) ! |s|^2
    call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, T, TT, CG(3)) ! r.s

    if (iter.gt.1) then
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, Y, CG(4)) ! myu
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, TT, CG(5)) ! nyu
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, Y, T, CG(6)) ! omega
      S2_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 6, HECMW_SUM)
      Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME
      gamma1 = CG(6)/CG(4) ! omega / myu
      gamma2 = CG(5)/CG(4) ! nyu / myu
    else
      S2_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
      Tcomm =  Tcomm + HECMW_WTIME() - S2_TIME
      gamma1 = 0.d0
      gamma2 = 0.d0
    endif

    c = CG(3) / dsqrt(CG(1)*CG(2))
    c_abs = dabs(c)
    if (c_abs > OMEGA) then
      QSI = c * dsqrt(CG(1)/CG(2))
    else
      if (c >= 0.d0) then
        QSI = OMEGA * dsqrt(CG(1)/CG(2))
      else
        QSI = -OMEGA * dsqrt(CG(1)/CG(2))
      endif
    endif
    ETA = gamma1 - QSI*gamma2
  end subroutine pol_coef_vanilla2

end subroutine  hecmw_solve_GPBiCG
end module     hecmw_solver_GPBiCG

!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C***
!C*** module hecmw_solver_GPBiCG
!C***
!
module hecmw_solver_GPBiCG

  private
  public :: hecmw_solve_GPBiCG

contains
  !C
  !C*** hecmw_solve_GPBiCG
  !C
  subroutine hecmw_solve_GPBiCG( hecMESH,  hecMAT, ITER, RESID, error, &
      &                              Tset, Tsol, Tcomm )
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

    integer(kind=kint ) ::  N, NP, NDOF, NNDOF
    integer(kind=kint ) ::  my_rank

    real   (kind=kreal), pointer ::  B(:)
    real   (kind=kreal), pointer ::  X(:)

    integer(kind=kint ) :: ITERlog, TIMElog

    real(kind=kreal), dimension(:,:),  allocatable       ::  WW

    real(kind=kreal), dimension(2) :: RR

    integer(kind=kint ) :: MAXIT
    real   (kind=kreal) :: TOL
    integer(kind=kint ) :: i,j
    real   (kind=kreal) :: S_TIME,S1_TIME,E_TIME,E1_TIME
    real   (kind=kreal) :: BNRM2
    real   (kind=kreal) :: RHO,RHO1,BETA,ALPHA,DNRM2
    real   (kind=kreal) :: QSI,ETA,COEF1
    real   (kind=kreal) :: t_max,t_min,t_avg,t_sd

    integer(kind=kint), parameter :: R= 1
    integer(kind=kint), parameter ::RT= 2
    integer(kind=kint), parameter :: T= 3
    integer(kind=kint), parameter ::TT= 4
    integer(kind=kint), parameter ::T0= 5
    integer(kind=kint), parameter :: P= 6
    integer(kind=kint), parameter ::PT= 7
    integer(kind=kint), parameter :: U= 8
    integer(kind=kint), parameter ::W1= 9
    integer(kind=kint), parameter :: Y=10
    integer(kind=kint), parameter :: Z=11
    integer(kind=kint), parameter ::WK=12
    integer(kind=kint), parameter ::W2=13
    integer(kind=kint), parameter ::ZQ=14

    integer(kind=kint), parameter :: N_ITER_RECOMPUTE_R= 20

    !! matrix integration for OpenACC
#ifdef _OPENACC
    integer(kind=kint), allocatable :: indexA(:), itemA(:)
    real(kind=kreal), allocatable   :: A(:)
    integer(kind=kint) :: nn, NPA, pre, pp, jS, jE, k
#endif

    call hecmw_barrier(hecMESH)
    S_TIME= HECMW_WTIME()
    !C
    !C-- INIT.
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

    error= 0
    BETA = 0.0d0

    allocate (WW(NDOF*NP,14))
    WW= 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    !C
    !C-- matrix integration for OpenACC
#ifdef _OPENACC
    nn = NDOF * NDOF
    NPA = NP + hecMAT%NPL + hecMAT%NPU
    allocate (indexA(0:NP), itemA(NPA), A(nn * NPA))
    indexA(0) = 0

    pre = 0
    pp = 0
    !$acc parallel loop private(i, j, k, pre, pp, jS, jE)
    do i = 1, NP
      indexA(i) = i + hecMAT%indexL(i) + hecMAT%indexU(i)

      pre = i - 1 + hecMAT%indexU(i - 1)
      jS= hecMAT%indexL(i - 1) + 1
      jE= hecMAT%indexL(i  )
      do j = jS, jE
        pp = pre + j
        itemA(pp) = hecMAT%itemL(j)
        do k = -nn+1, 0
          A(nn * pp + k) = hecMAT%AL(nn * j + k)
        enddo
      enddo

      pp = i + hecMAT%indexU(i - 1) + hecMAT%indexL(i)
      itemA(pp) = i
      do k = -nn+1, 0
        A(nn * pp + k) = hecMAT%D(nn * i + k)
      enddo

      pre = i + hecMAT%indexL(i)
      jS= hecMAT%indexU(i - 1) + 1
      jE= hecMAT%indexU(i  )
      do j = jS, jE
        pp = pre + j
        itemA(pp) = hecMAT%itemU(j)
        do k = -nn+1, 0
          A(nn * pp + k) = hecMAT%AU(nn * j + k)
        enddo
      enddo
    enddo
    !$acc end parallel
#endif
    !C===
    !C +----------------------+
    !C | SETUP PRECONDITIONER |
    !C +----------------------+
    !C===
    call hecmw_precond_setup(hecMAT, hecMESH, 0)

    !C
    !C +----------------------+
    !C | {r}= {b} - [A]{xini} |
    !C +----------------------+
    !C===
#ifdef _OPENACC
    call hecmw_matresid_A(hecMESH, hecMAT, indexA, itemA, A, X, B, WW(:,R), Tcomm)
#else
    call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
#endif

    !$acc kernels
    !$acc loop independent
    do i= 1, NNDOF
      WW(i,RT)= WW(i,R)
    enddo
    !$acc end kernels
    !C==
    call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
    if (BNRM2.eq.0.d0) then
      iter = 0
      MAXIT = 0
      RESID = 0.d0
      X = 0.d0
    endif

    call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,RT), WW(:,R), RHO, Tcomm)

    E_TIME= HECMW_WTIME()
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
    !C===

    !C
    !C*************************************************************** ITERATIVE PROC.
    !C
    call hecmw_barrier(hecMESH)
    S1_TIME= HECMW_WTIME()
    do iter= 1, MAXIT
      !C
      !C +----------------+
      !C | {r}= [Minv]{r} |
      !C +----------------+
      !C===
      !$acc kernels
      !$acc loop independent
      do j= 1, NNDOF
        WW(j,WK)=  WW(j,R)
      enddo
      !$acc end kernels

#ifdef _OPENACC
      call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,WK), WW(:,R), WW(:,ZQ), Tcomm)
#else
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,WK), WW(:,R), WW(:,ZQ), Tcomm)
#endif
      !C===

      !C
      !C +----------------------------------+
      !C | {p} = {r} + BETA * ( {p} - {u} ) |
      !C +----------------------------------+
      !C===
      if (iter.gt.1) then
        !$acc kernels
        !$acc loop independent
        do j= 1, NNDOF
          WW(j,P)= WW(j,R) + BETA*( WW(j,P)-WW(j,U))
        enddo
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j= 1, NNDOF
          WW(j,P)= WW(j,R)
        enddo
        !$acc end kernels
      endif
      !C===

      !C
      !C +--------------------------------+
      !C | ALPHA= {r_tld}{r}/{r_tld} A{p} |
      !C +--------------------------------+
      !C===

      !C
      !C-- calc. {p_tld}= A{p}
#ifdef _OPENACC
      call hecmw_matvec_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,P), WW(:,PT), Tcomm)
#else
      call hecmw_matvec(hecMESH, hecMAT, WW(:,P), WW(:,PT), Tcomm)
#endif

      !C
      !C-- calc. ALPHA
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,RT), WW(:,PT), RHO1, Tcomm)

      ALPHA= RHO / RHO1
      !C===

      !C
      !C +------------------------------------------+
      !C | {y}= {t} - {r} - ALPHA{w} + ALPHA{p_tld} |
      !C | {t}= {r}                  - ALPHA{p_tld} |
      !C +------------------------------------------+
      !C===
      !$acc kernels
      !$acc loop independent
      do j= 1, NNDOF
        WW(j,Y)= WW(j,T) - WW(j,WK) + ALPHA*(-WW(j,W1) + WW(j,PT))
        WW(j,T)= WW(j,WK) - ALPHA*WW(j,PT)
      enddo
      !$acc end kernels
      !C===

      !C
      !C +-----------------------+
      !C | {t_tld}= [A][Minv]{t} |
      !C +-----------------------+
      !C===

      !C
      !C-- calc. {t_tld} and {t0} by [M] inversion
      !C         {W2}   = [Minv]{p_tld}
      !C
#ifdef _OPENACC
      call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,T), WW(:,TT), WW(:,ZQ), Tcomm)
#else
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,T), WW(:,TT), WW(:,ZQ), Tcomm)
#endif
#ifdef _OPENACC
      call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,T0), WW(:,W2), WW(:,ZQ), Tcomm)
#else
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,T0), WW(:,W2), WW(:,ZQ), Tcomm)
#endif
      !$acc kernels
      !$acc loop independent
      do i= 1, NNDOF
        WW(i,T0)= WW(i,W2)
      enddo
      !$acc end kernels
#ifdef _OPENACC
      call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,PT), WW(:,W2), WW(:,ZQ), Tcomm)
#else
      call hecmw_precond_apply(hecMESH, hecMAT, WW(:,PT), WW(:,W2), WW(:,ZQ), Tcomm)
#endif

      !C===

      !C
      !C-- calc. [A]{t_tld}
#ifdef _OPENACC
      call hecmw_matvec_A(hecMESH, hecMAT, indexA, itemA, A, WW(:,TT), WW(:,WK), Tcomm)
#else
      call hecmw_matvec(hecMESH, hecMAT, WW(:,TT), WW(:,WK), Tcomm)
#endif

      !$acc kernels
      !$acc loop independent
      do i= 1, NNDOF
        WW(i,TT)= WW(i,WK)
      enddo
      !$acc end kernels
      !C===

      !C
      !C +-------------------+
      !C | calc. QSI and ETA |
      !C +-------------------+
      !C===
      !call pol_coef(iter, WW, T, TT, Y, QSI, ETA)
      !call pol_coef_vanilla(iter, WW, T, TT, Y, QSI, ETA)
      call pol_coef_vanilla2(iter, WW, T, TT, Y, QSI, ETA)
      !C===

      !C
      !C +----------------------------------------------------------+
      !C | {u} = QSI [Minv]{pt} + ETA([Minv]{t0}-[Minv]{r}+BETA*{u} |
      !C | {z} = QSI [Minv]{r}  + ETA*{z} - ALPHA*{u}               |
      !C +----------------------------------------------------------+
      !C===

      !C
      !C-- compute. {u},{z}

      if (iter.gt.1) then
        !$acc kernels
        !$acc loop independent
        do j= 1, NNDOF
          WW(j,U)= QSI* WW(j,W2) + ETA*(WW(j,T0) - WW(j,R) + BETA*WW(j,U))
          WW(j,Z)= QSI* WW(j, R) + ETA* WW(j, Z) -          ALPHA*WW(j,U)
        enddo
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j= 1, NNDOF
          WW(j,U)= QSI* WW(j,W2) + ETA*(WW(j,T0) - WW(j,R))
          WW(j,Z)= QSI* WW(j, R) + ETA* WW(j, Z) -          ALPHA*WW(j,U)
        enddo
        !$acc end kernels
      endif
      !C===

      !C
      !C +--------------------+
      !C | update {x},{r},{w} |
      !C +--------------------+
      !C===
      !$acc kernels
      !$acc loop independent
      do j= 1, NNDOF
        X (j)= X(j) + ALPHA*WW(j,P) + WW(j,Z)
        ! WW(j,R)= WW(j,T) - ETA*WW(j,Y) - QSI*WW(j,TT)
        WW(j,T0)= WW(j,T)
      enddo
      !$acc end kernels
      !C
      !C--- recompute R sometimes
      if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) then
#ifdef _OPENACC
        call hecmw_matresid_A(hecMESH, hecMAT, indexA, itemA, A, X, B, WW(:,R), Tcomm)
#else
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
#endif
      else
        !$acc kernels
        !$acc loop independent
        do j= 1, NNDOF
          WW(j,R)= WW(j,T) - ETA*WW(j,Y) - QSI*WW(j,TT)
        enddo
        !$acc end kernels
      endif

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,R), RR(1))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,R), WW(:,RT), RR(2))
      S_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, RR, 2, HECMW_SUM)
      E_TIME= HECMW_WTIME()
      Tcomm =  Tcomm + E_TIME - S_TIME
      DNRM2 = RR(1)
      COEF1 = RR(2)

      BETA = ALPHA*COEF1 / (QSI*RHO)
      !$acc kernels
      !$acc loop independent
      do j= 1, NNDOF
        WW(j,W1)= WW(j,TT) + BETA*WW(j,PT)
      enddo
      !$acc end kernels

      RESID= dsqrt(DNRM2/BNRM2)
      RHO  = COEF1

      !C##### ITERATION HISTORY
      if (my_rank.eq.0 .and. ITERlog.eq.1)                              &
        &    write (*, 1000) ITER, RESID
      1000   format (i5, 1pe16.6)
      !C#####

      if (RESID.le.TOL   ) then
        if ( mod(ITER,N_ITER_RECOMPUTE_R)==0 ) exit
        !C----- recompute R to make sure it is really converged
#ifdef _OPENACC
        call hecmw_matresid_A(hecMESH, hecMAT, indexA, itemA, A, X, B, WW(:,R), Tcomm)
#else
        call hecmw_matresid(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
#endif
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
        RESID= dsqrt(DNRM2/BNRM2)
        if (RESID.le.TOL) exit
      endif
      if ( ITER.eq.MAXIT ) error= HECMW_SOLVER_ERROR_NOCONV_MAXIT
      !C===
    enddo

    call hecmw_solver_scaling_bk(hecMAT)
    !C
    !C-- INTERFACE data EXCHANGE

    S_TIME = HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
    E_TIME = HECMW_WTIME()
    Tcomm = Tcomm + E_TIME - S_TIME

    deallocate (WW)
    !call hecmw_precond_clear(hecMAT)

    E1_TIME= HECMW_WTIME()
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

  contains

    !C
    !C*** pol_coef : computes QSI and ETA in original GPBiCG way
    !C
    subroutine pol_coef(iter, WW, T, TT, Y, QSI, ETA)
      implicit none
      integer(kind=kint), intent(in) :: iter
      real(kind=kreal), intent(in) :: WW(:,:)
      integer(kind=kint), intent(in) :: T, TT, Y
      real(kind=kreal), intent(out) :: QSI, ETA

      real(kind=kreal), dimension(5) :: CG
      real(kind=kreal), dimension(2) :: EQ
      real(kind=kreal) :: delta

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,Y ), CG(1)) ! myu1
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,T ), CG(2)) ! omega2
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,T ), CG(3)) ! omega1
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,Y ), CG(4)) ! nyu
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,TT), CG(5)) ! myu2
      S_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 5, HECMW_SUM)
      E_TIME= HECMW_WTIME()
      Tcomm =  Tcomm + E_TIME - S_TIME

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

    !C
    !C*** pol_coef_vanilla : computes QSI and ETA with vanilla strategy
    !C      see Fujino, Abe, Sugihara and Nakashima(2013) ISBN978-4-621-08741-1
    !C
    subroutine pol_coef_vanilla(iter, WW, T, TT, Y, QSI, ETA)
      implicit none
      integer(kind=kint), intent(in) :: iter
      real(kind=kreal), intent(inout) :: WW(:,:)
      integer(kind=kint), intent(in) :: T, TT, Y
      real(kind=kreal), intent(out) :: QSI, ETA

      real(kind=kreal), dimension(3) :: CG
      real(kind=kreal) :: gamma1, gamma2
      real(kind=kreal) :: c, c_abs

      real(kind=kreal), parameter :: OMEGA = 0.707106781d0

      if (iter.eq.1) then
        gamma1 = 0.d0
        gamma2 = 0.d0
      else
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,Y ), CG(1)) ! myu
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,TT), CG(2)) ! nyu
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,T ), CG(3)) ! omega
        S_TIME= HECMW_WTIME()
        call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
        E_TIME= HECMW_WTIME()
        Tcomm =  Tcomm + E_TIME - S_TIME
        gamma1 = CG(3)/CG(1) ! omega / myu
        gamma2 = CG(2)/CG(1) ! nyu / myu
        !!! COMMENTED OUT because no convergence obtained with following updates.
        !         do j= 1, NNDOF
        !           WW(j,T )= WW(j,T ) - gamma1*WW(j,Y)
        !           WW(j,TT)= WW(j,TT) - gamma2*WW(j,Y)
        !         enddo
      endif

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,T ), WW(:,T ), CG(1)) ! |r|^2
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,TT), CG(2)) ! |s|^2
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,T ), WW(:,TT), CG(3)) ! r.s
      S_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
      E_TIME= HECMW_WTIME()
      Tcomm =  Tcomm + E_TIME - S_TIME

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

    !C
    !C*** pol_coef_vanilla2 : optimized version of pol_coef_vanilla
    !C
    subroutine pol_coef_vanilla2(iter, WW, T, TT, Y, QSI, ETA)
      implicit none
      integer(kind=kint), intent(in) :: iter
      real(kind=kreal), intent(inout) :: WW(:,:)
      integer(kind=kint), intent(in) :: T, TT, Y
      real(kind=kreal), intent(out) :: QSI, ETA

      real(kind=kreal), dimension(6) :: CG
      real(kind=kreal) :: gamma1, gamma2
      real(kind=kreal) :: c, c_abs

      real(kind=kreal), parameter :: OMEGA = 0.707106781d0

      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,T ), WW(:,T ), CG(1)) ! |r|^2
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,TT), CG(2)) ! |s|^2
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,T ), WW(:,TT), CG(3)) ! r.s

      if (iter.gt.1) then
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,Y ), CG(4)) ! myu
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,TT), CG(5)) ! nyu
        call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,T ), CG(6)) ! omega
        S_TIME= HECMW_WTIME()
        call hecmw_allreduce_R(hecMESH, CG, 6, HECMW_SUM)
        E_TIME= HECMW_WTIME()
        Tcomm =  Tcomm + E_TIME - S_TIME
        gamma1 = CG(6)/CG(4) ! omega / myu
        gamma2 = CG(5)/CG(4) ! nyu / myu
      else
        S_TIME= HECMW_WTIME()
        call hecmw_allreduce_R(hecMESH, CG, 3, HECMW_SUM)
        E_TIME= HECMW_WTIME()
        Tcomm =  Tcomm + E_TIME - S_TIME
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

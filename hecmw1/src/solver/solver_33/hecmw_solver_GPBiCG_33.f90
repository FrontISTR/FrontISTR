!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                       Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!C*** 
!C*** module hecmw_solver_GPBiCG_33
!C***
!
      module hecmw_solver_GPBiCG_33
      contains
!C
!C*** hecmw_solve_GPBiCG_33
!C
      subroutine hecmw_solve_GPBiCG_33( hecMESH,  hecMAT, ITER, RESID, ERROR, &
     &                              Tset, Tsol, Tcomm )
      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_matrix_misc
      use hecmw_solver_misc
      use hecmw_solver_las_33
      use hecmw_solver_scaling_33
      use hecmw_precond_33

      implicit none

      type(hecmwST_local_mesh) :: hecMESH
      type(hecmwST_matrix) :: hecMAT
      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      integer(kind=kint ) ::  N, NP, NDOF, NNDOF
      integer(kind=kint ) ::  my_rank

      real   (kind=kreal), pointer ::  B(:)
      real   (kind=kreal), pointer ::  X(:)

      integer(kind=kint ) :: ITERlog, TIMElog

      real(kind=kreal), dimension(:,:),  allocatable       ::  WW

      real(kind=kreal), dimension(5) :: CG
      real(kind=kreal), dimension(2) :: EQ

      integer(kind=kint ) :: MAXIT
      real   (kind=kreal) :: TOL
      integer(kind=kint ) :: i,j
      real   (kind=kreal) :: S_TIME,S1_TIME,E_TIME,E1_TIME
      real   (kind=kreal) :: BNRM2
      real   (kind=kreal) :: RHO,RHO1,BETA,ALPHA,DNRM2
      real   (kind=kreal) :: QSI,ETA,COEF1

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

      ERROR= 0

      allocate (WW(NDOF*NP,14))
      WW= 0.d0

!C
!C-- SCALING
      call hecmw_solver_scaling_fw_33(hecMESH, hecMAT, Tcomm)

!C===
!C +----------------------+
!C | SETUP PRECONDITIONER |
!C +----------------------+
!C===
      call hecmw_precond_33_setup(hecMAT)

!C
!C +----------------------+
!C | {r}= {b} - [A]{xini} |
!C +----------------------+
!C===
      call hecmw_matresid_33(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)

      do i= 1, NNDOF
        WW(i,RT)= WW(i,R)
      enddo
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
      Tset= Tset + E_TIME - S_TIME
!C===

!C
!C*************************************************************** ITERATIVE PROC.
!C
      S1_TIME= HECMW_WTIME()
      do iter= 1, MAXIT
!C
!C +----------------+
!C | {r}= [Minv]{r} |
!C +----------------+
!C===
      do j= 1, NNDOF
        WW(j,WK)=  WW(j,R)
      enddo

      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,WK), WW(:,R), WW(:,ZQ), Tcomm)
!C===

!C
!C +----------------------------------+
!C | {p} = {r} + BETA * ( {p} - {u} ) |
!C +----------------------------------+
!C===
      if (iter.gt.1) then
        do j= 1, NNDOF
          WW(j,P)= WW(j,R) + BETA*( WW(j,P)-WW(j,U))
        enddo
       else
        do j= 1, NNDOF
          WW(j,P)= WW(j,R)
        enddo
      endif
!C===

!C
!C +--------------------------------+
!C | ALPHA= {r_tld}{r}/{r_tld} A{p} |
!C +--------------------------------+
!C===

!C
!C-- calc. {p_tld}= A{p}
      call hecmw_matvec_33(hecMESH, hecMAT, WW(:,P), WW(:,PT), Tcomm)

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
      do j= 1, NNDOF
        WW(j,Y)= WW(j,T) - WW(j,WK) + ALPHA*(-WW(j,W1) + WW(j,PT))
        WW(j,T)= WW(j,WK) - ALPHA*WW(j,PT)
      enddo
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
      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,T), WW(:,TT), WW(:,ZQ), Tcomm)
      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,T0), WW(:,W2), WW(:,ZQ), Tcomm)
      do i= 1, NNDOF
        WW(i,T0)= WW(i,W2)
      enddo
      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,PT), WW(:,W2), WW(:,ZQ), Tcomm)

!C===

!C
!C-- calc. [A]{t_tld}
      call hecmw_matvec_33(hecMESH, hecMAT, WW(:,TT), WW(:,WK), Tcomm)

      do i= 1, NNDOF
        WW(i,TT)= WW(i,WK)
      enddo
!C===

!C
!C +-------------------+
!C | calc. QSI and ETA |
!C +-------------------+
!C===
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,Y ), CG(1))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,T ), CG(2))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,Y ), WW(:,T ), CG(NDOF))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,Y ), CG(4))
      call hecmw_InnerProduct_R_nocomm(hecMESH, NDOF, WW(:,TT), WW(:,TT), CG(5))
      S_TIME= HECMW_WTIME()
      call hecmw_allreduce_R(hecMESH, CG, 5, HECMW_SUM)
      E_TIME= HECMW_WTIME()
      Tcomm =  Tcomm + E_TIME - S_TIME

      if (iter.eq.1) then
        EQ(1)= CG(2)/CG(5)
        EQ(2)= 0.d0
       else
        EQ(1)= (CG(1)*CG(2)-CG(3)*CG(4)) / (CG(5)*CG(1)-CG(4)*CG(4))
        EQ(2)= (CG(5)*CG(3)-CG(4)*CG(2)) / (CG(5)*CG(1)-CG(4)*CG(4))
      endif

      QSI= EQ(1)
      ETA= EQ(2)
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
        do j= 1, NNDOF
          WW(j,U)= QSI* WW(j,W2) + ETA*(WW(j,T0) - WW(j,R) + BETA*WW(j,U))
          WW(j,Z)= QSI* WW(j, R) + ETA* WW(j, Z) -          ALPHA*WW(j,U)
        enddo
      else
        do j= 1, NNDOF
          WW(j,U)= QSI* WW(j,W2) + ETA*(WW(j,T0) - WW(j,R))
          WW(j,Z)= QSI* WW(j, R) + ETA* WW(j, Z) -          ALPHA*WW(j,U)
        enddo
      endif
!C===

!C
!C +--------------------+
!C | update {x},{r},{w} |
!C +--------------------+
!C===
      do j= 1, NNDOF
        X (j)= X(j) + ALPHA*WW(j,P) + WW(j,Z)
        WW(j,R)= WW(j,T) - ETA*WW(j,Y) - QSI*WW(j,TT)
        WW(j,T0)= WW(j,T)
      enddo

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,RT), COEF1, Tcomm)


      BETA = ALPHA*COEF1 / (QSI*RHO)
      do j= 1, NNDOF
        WW(j,W1)= WW(j,TT) + BETA*WW(j,PT)
      enddo

      RESID= dsqrt(DNRM2/BNRM2)
      RHO  = COEF1

!C##### ITERATION HISTORY
      if (my_rank.eq.0 .and. ITERlog.eq.1)                              &
     &    write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
!C#####

      if (RESID.le.TOL   ) exit
      if ( ITER.eq.MAXIT ) ERROR= -300
!C===
      enddo

      call hecmw_solver_scaling_bk_33(hecMAT)
!C
!C-- INTERFACE data EXCHANGE

      S_TIME = HECMW_WTIME()
      call hecmw_update_3_R (hecMESH, X, hecMAT%NP)
      E_TIME = HECMW_WTIME()
      Tcomm = Tcomm + E_TIME - S_TIME

      E1_TIME= HECMW_WTIME()
      Tsol = E1_TIME - S1_TIME

      deallocate (WW)
      call hecmw_precond_33_clear(hecMAT)

      end subroutine  hecmw_solve_GPBiCG_33
      end module     hecmw_solver_GPBiCG_33

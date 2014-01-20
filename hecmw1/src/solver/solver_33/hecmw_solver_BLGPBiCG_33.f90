!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.5                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!
!C*** 
!C*** module hecmw_solver_BLGPBiCG_33
!C***
!
      module hecmw_solver_BLGPBiCG_33
      contains
!C
!C*** hecmw_solve_BLGPBiCG_33
!C
!C      (1) Block ILU(0) 
!C      (2) Block ILU(1) 
!C      (3) Block ILU(2) 
!C

      subroutine hecmw_solve_BLGPBiCG_33                                &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, RESID, SIGMA, SIGMA_DIAG, ITER, ERROR, my_rank,         &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM , PRECOND, iterPREmax,         &
     &                       Tset, Tsol, Tcomm, ITERlog)

      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_solver_SR_33
      use hecmw_precond_BILU_33
      implicit none

      integer(kind=kint ), intent(in):: N, NP, NPU, NPL, my_rank
      integer(kind=kint ), intent(in):: ITERlog
      integer(kind=kint ), intent(in):: NEIBPETOT, iterPREmax, PRECOND
      integer(kind=kint ), intent(in):: SOLVER_COMM
      real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      real(kind=kreal), dimension(3*NP)   , intent(inout):: B, X
      real(kind=kreal), dimension(9*NPL), intent(inout):: AL
      real(kind=kreal), dimension(9*NPU), intent(inout):: AU
      real(kind=kreal), dimension(9*NP ), intent(inout):: D

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)

      real(kind=kreal), dimension(:),    allocatable       :: WS, WR
      real(kind=kreal), dimension(:,:),  allocatable       :: WW

      real(kind=kreal), dimension(5) :: C0, CG
      real(kind=kreal), dimension(2) :: EQ

      integer(kind=kint ) :: R, RT, T, T0, TT, P, PT
      integer(kind=kint ) :: U, W1, Y, Z, WK, W2, MAXIT
      real   (kind=kreal) :: TOL
      integer(kind=kint ) :: i,j,k
      integer(kind=kint ) :: ns, nr
      real   (kind=kreal) :: S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal) :: BNRM20,BNRM2,X1,X2,X3
      real   (kind=kreal) :: RHO,RHO0,RHO1,BETA,ALPHA,DNRM20,DNRM2,RHO10
      real   (kind=kreal) :: WVAL1,WVAL2,WVAL3
      real   (kind=kreal) :: COMMtime,QSI,ETA,COEF10,COEF1

      START_TIME= HECMW_WTIME()
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      if( NEIBPETOT > 0 ) then
        ns = STACK_EXPORT(NEIBPETOT)
        nr = STACK_IMPORT(NEIBPETOT)
      else
        ns = 0
        nr = 0
      end if

      allocate (WW(3*NP,13))
      allocate (WS(3*ns), WR(3*nr))

      WW= 0.d0
      ERROR= 0

      MAXIT = ITER
      TOL   = RESID

      ERROR= 0

       R= 1
      RT= 2
       T= 3
      TT= 4
      T0= 5
       P= 6
      PT= 7
       U= 8
      W1= 9
       Y=10
       Z=11
      WK=12
      W2=13

      COMMtime= 0.d0

!C
!C-- exchanging DIAGONAL components
        WW= 0.d0
        do i= 1, N
          WW(3*i-2,1)= D(9*i-8)
          WW(3*i-2,2)= D(9*i-7)
          WW(3*i-2,3)= D(9*i-6)
          WW(3*i-1,1)= D(9*i-5)
          WW(3*i-1,2)= D(9*i-4)
          WW(3*i-1,3)= D(9*i-3)
          WW(3*i  ,1)= D(9*i-2)
          WW(3*i  ,2)= D(9*i-1)
          WW(3*i  ,3)= D(9*i  )
        enddo

        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,1),                 &
     &       SOLVER_COMM,my_rank )
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,2),                 &
     &       SOLVER_COMM,my_rank )
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,3),                 &
     &       SOLVER_COMM,my_rank )

        do i= N+1, NP
          D(9*i-8)= WW(3*i-2,1)
          D(9*i-7)= WW(3*i-2,2)
          D(9*i-6)= WW(3*i-2,3)
          D(9*i-5)= WW(3*i-1,1)
          D(9*i-4)= WW(3*i-1,2)
          D(9*i-3)= WW(3*i-1,3)
          D(9*i-2)= WW(3*i  ,1)
          D(9*i-1)= WW(3*i  ,2)
          D(9*i  )= WW(3*i  ,3)
        enddo

!C===

!C
!C +-------------------+
!C | ILU decomposition |
!C +-------------------+
!C===
      call hecmw_precond_BILU_33_setup                    &
           &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
           &    SIGMA, SIGMA_DIAG, my_rank,                     &
           &    PRECOND)

!C===

      WW= 0.d0
!C
!C +----------------------+
!C | {r}= {b} - [A]{xini} |
!C +----------------------+
!C===
!C
!C-- INTERFACE data EXCHANGE
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

!C
!C-- BEGIN calculation

      do j= 1, N
           X1= X(3*j-2)
           X2= X(3*j-1)
           X3= X(3*j  )
        WVAL1= B(3*j-2) - D(9*j-8)*X1 - D(9*j-7)*X2 - D(9*j-6)*X3
        WVAL2= B(3*j-1) - D(9*j-5)*X1 - D(9*j-4)*X2 - D(9*j-3)*X3
        WVAL3= B(3*j  ) - D(9*j-2)*X1 - D(9*j-1)*X2 - D(9*j  )*X3

        do k= INL(j-1)+1, INL(j)
          i= IAL(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
          WVAL2= WVAL2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
          WVAL3= WVAL3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
        enddo

        do k= INU(j-1)+1, INU(j)
          i= IAU(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
          WVAL2= WVAL2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
          WVAL3= WVAL3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
        enddo

        WW(3*j-2,R)= WVAL1
        WW(3*j-1,R)= WVAL2
        WW(3*j  ,R)= WVAL3
        WW(3*j-2,RT)= WVAL1
        WW(3*j-1,RT)= WVAL2
        WW(3*j  ,RT)= WVAL3
      enddo

      BNRM20= 0.d0
      RHO0  = 0.0d0

      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
        RHO0  = RHO0 + WW(3*i-2,RT)*WW(3*i-2,R)+WW(3*i-1,RT)*WW(3*i-1,R)&
     &                                         +WW(3*i  ,RT)*WW(3*i  ,R)
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0  , RHO,   1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) then
        iter = 0
        MAXIT = 0
        RESID = 0.d0
        X = 0.d0
      endif

       END_TIME= HECMW_WTIME()
      SETupTIME= END_TIME - START_TIME
!C===

      START_TIME= HECMW_WTIME()
!C
!C*************************************************************** ITERATIVE PROC.
!C
      do iter= 1, MAXIT
!C
!C-- INIT.

      do j= 1, N
        WW(3*j-2,WK)=  WW(3*j-2,R)
        WW(3*j-1,WK)=  WW(3*j-1,R)
        WW(3*j  ,WK)=  WW(3*j  ,R)
      enddo

!C
!C +----------------+
!C | {r}= [Minv]{r} |
!C +----------------+
!C===

!C
!C== Block ILU
      do i= 1+N, NP
        WW(3*i-2,R)= 0.d0
        WW(3*i-1,R)= 0.d0
        WW(3*i  ,R)= 0.d0
      enddo

      call hecmw_precond_BILU_33_apply(N, WW(:,R))

!C==

!C===

!C
!C +----------------------------------+
!C | {p} = {r} + BETA * ( {p} - {u} ) |
!C +----------------------------------+
!C===
      if (iter.gt.1) then
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R) + BETA*( WW(3*j-2,P)-WW(3*j-2,U))
          WW(3*j-1,P)= WW(3*j-1,R) + BETA*( WW(3*j-1,P)-WW(3*j-1,U))
          WW(3*j  ,P)= WW(3*j  ,R) + BETA*( WW(3*j  ,P)-WW(3*j  ,U))
        enddo
       else
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R)
          WW(3*j-1,P)= WW(3*j-1,R)
          WW(3*j  ,P)= WW(3*j  ,R)
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

      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,P) , SOLVER_COMM,     &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      do j= 1, N
           X1= WW(3*j-2,P)
           X2= WW(3*j-1,P)
           X3= WW(3*j  ,P)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3

        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo

        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,PT)= WVAL1
        WW(3*j-1,PT)= WVAL2
        WW(3*j  ,PT)= WVAL3
      enddo

!C
!C-- calc. ALPHA

      RHO10= 0.d0
      do j= 1, N
        RHO10= RHO10+WW(3*j-2,RT)*WW(3*j-2,PT)+WW(3*j-1,RT)*WW(3*j-1,PT)&
     &                                        +WW(3*j  ,RT)*WW(3*j  ,PT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO10, RHO1, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO10, RHO1, 1, MPI_DOUBLE_PRECISION,         &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      ALPHA= RHO / RHO1
!C===

!C
!C +------------------------------------------+
!C | {y}= {t} - {r} - ALPHA{w} + ALPHA{p_tld} |
!C | {t}= {r}                  - ALPHA{p_tld} |
!C +------------------------------------------+
!C===
      do j= 1, N
        WW(3*j-2,Y)= WW(3*j-2,T) - WW(3*j-2,WK)                         &
     &                   + ALPHA*(-WW(3*j-2,W1) + WW(3*j-2,PT))
        WW(3*j-1,Y)= WW(3*j-1,T) - WW(3*j-1,WK)                         &
     &                   + ALPHA*(-WW(3*j-1,W1) + WW(3*j-1,PT))
        WW(3*j  ,Y)= WW(3*j  ,T) - WW(3*j  ,WK)                         &
     &                   + ALPHA*(-WW(3*j  ,W1) + WW(3*j  ,PT))

        WW(3*j-2,T)= WW(3*j-2,WK) - ALPHA*WW(3*j-2,PT)
        WW(3*j-1,T)= WW(3*j-1,WK) - ALPHA*WW(3*j-1,PT)
        WW(3*j  ,T)= WW(3*j  ,WK) - ALPHA*WW(3*j  ,PT)
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

!C
!C== Block ILU
      do i= 1, N
        WW(3*i-2,TT)= WW(3*i-2,T )
        WW(3*i-1,TT)= WW(3*i-1,T )
        WW(3*i  ,TT)= WW(3*i  ,T )
        WW(3*i-2,W2)= WW(3*i-2,PT)
        WW(3*i-1,W2)= WW(3*i-1,PT)
        WW(3*i  ,W2)= WW(3*i  ,PT)
      enddo

      do i= N+1, NP
        WW(3*i-2,TT)= 0.d0
        WW(3*i-1,TT)= 0.d0
        WW(3*i  ,TT)= 0.d0
        WW(3*i-2,W2)= 0.d0
        WW(3*i-1,W2)= 0.d0
        WW(3*i  ,W2)= 0.d0
        WW(3*i-2,T0)= 0.d0
        WW(3*i-1,T0)= 0.d0
        WW(3*i  ,T0)= 0.d0
      enddo

      call hecmw_precond_BILU_33_apply(N, WW(:,TT))
      call hecmw_precond_BILU_33_apply(N, WW(:,W2))
      call hecmw_precond_BILU_33_apply(N, WW(:,T0))

!C==
!C===

!C
!C-- calc. [A]{t_tld}
      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,TT) , SOLVER_COMM,    &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      do j= 1, N
           X1= WW(3*j-2,TT)
           X2= WW(3*j-1,TT)
           X3= WW(3*j  ,TT)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3

        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,TT)
          X2= WW(3*i-1,TT)
          X3= WW(3*i  ,TT)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo

        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,TT)
          X2= WW(3*i-1,TT)
          X3= WW(3*i  ,TT)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,WK)= WVAL1
        WW(3*j-1,WK)= WVAL2
        WW(3*j  ,WK)= WVAL3
      enddo

      do j= 1, N
        WW(3*j-2,TT)= WW(3*j-2,WK)
        WW(3*j-1,TT)= WW(3*j-1,WK)
        WW(3*j  ,TT)= WW(3*j  ,WK)
      enddo
!C===

!C
!C +-------------------+
!C | calc. QSI and ETA |
!C +-------------------+
!C===
      C0(1)= 0.0d0
      C0(2)= 0.0d0
      C0(3)= 0.0d0
      C0(4)= 0.0d0
      C0(5)= 0.0d0
      CG(1)= 0.0d0
      CG(2)= 0.0d0
      CG(3)= 0.0d0
      CG(4)= 0.0d0
      CG(5)= 0.0d0

      do j= 1, N
        C0(1)= C0(1)+WW(3*j-2, Y)*WW(3*j-2, Y)+WW(3*j-1, Y)*WW(3*j-1, Y)&
     &                                        +WW(3*j  , Y)*WW(3*j  , Y)
        C0(2)= C0(2)+WW(3*j-2,TT)*WW(3*j-2, T)+WW(3*j-1,TT)*WW(3*j-1, T)&
     &                                        +WW(3*j  ,TT)*WW(3*j  , T)
        C0(3)= C0(3)+WW(3*j-2, Y)*WW(3*j-2, T)+WW(3*j-1, Y)*WW(3*j-1, T)&
     &                                        +WW(3*j  , Y)*WW(3*j  , T)
        C0(4)= C0(4)+WW(3*j-2,TT)*WW(3*j-2, Y)+WW(3*j-1,TT)*WW(3*j-1, Y)&
     &                                        +WW(3*j  ,TT)*WW(3*j  , Y)
        C0(5)= C0(5)+WW(3*j-2,TT)*WW(3*j-2,TT)+WW(3*j-1,TT)*WW(3*j-1,TT)&
     &                                        +WW(3*j  ,TT)*WW(3*j  ,TT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP (C0, CG, 5, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C0, CG,  5, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

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
        do j= 1, N
          WW(3*j-2,U)= QSI* WW(3*j-2,W2) +                                &
     &                 ETA*(WW(3*j-2,T0) - WW(3*j-2,R) + BETA*WW(3*j-2,U))
          WW(3*j-2,Z)= QSI* WW(3*j-2, R) +                                &
     &                 ETA* WW(3*j-2, Z) -              ALPHA*WW(3*j-2,U)
          WW(3*j-1,U)= QSI* WW(3*j-1,W2) +                                &
     &                 ETA*(WW(3*j-1,T0) - WW(3*j-1,R) + BETA*WW(3*j-1,U))
          WW(3*j-1,Z)= QSI* WW(3*j-1, R) +                                &
     &                 ETA* WW(3*j-1, Z) -              ALPHA*WW(3*j-1,U)
          WW(3*j  ,U)= QSI* WW(3*j  ,W2) +                                &
     &                 ETA*(WW(3*j  ,T0) - WW(3*j  ,R) + BETA*WW(3*j  ,U))
          WW(3*j  ,Z)= QSI* WW(3*j  , R) +                                &
     &                 ETA* WW(3*j  , Z) -              ALPHA*WW(3*j  ,U)
        enddo
      else
        do j= 1, N
          WW(3*j-2,U)= QSI* WW(3*j-2,W2) +                                &
     &                 ETA*(WW(3*j-2,T0) - WW(3*j-2,R))
          WW(3*j-2,Z)= QSI* WW(3*j-2, R) +                                &
     &                 ETA* WW(3*j-2, Z) -              ALPHA*WW(3*j-2,U)
          WW(3*j-1,U)= QSI* WW(3*j-1,W2) +                                &
     &                 ETA*(WW(3*j-1,T0) - WW(3*j-1,R))
          WW(3*j-1,Z)= QSI* WW(3*j-1, R) +                                &
     &                 ETA* WW(3*j-1, Z) -              ALPHA*WW(3*j-1,U)
          WW(3*j  ,U)= QSI* WW(3*j  ,W2) +                                &
     &                 ETA*(WW(3*j  ,T0) - WW(3*j  ,R))
          WW(3*j  ,Z)= QSI* WW(3*j  , R) +                                &
     &                 ETA* WW(3*j  , Z) -              ALPHA*WW(3*j  ,U)
        enddo
      endif
!C===

!C
!C +--------------------+
!C | update {x},{r},{w} |
!C +--------------------+
!C===
      DNRM20= 0.d0
      COEF10= 0.d0

      do j= 1, N
        X (3*j-2)= X(3*j-2) + ALPHA*WW(3*j-2,P) + WW(3*j-2,Z)
        X (3*j-1)= X(3*j-1) + ALPHA*WW(3*j-1,P) + WW(3*j-1,Z)
        X (3*j  )= X(3*j  ) + ALPHA*WW(3*j  ,P) + WW(3*j  ,Z)

        WW(3*j-2,R)= WW(3*j-2,T) - ETA*WW(3*j-2,Y) - QSI*WW(3*j-2,TT)
        WW(3*j-1,R)= WW(3*j-1,T) - ETA*WW(3*j-1,Y) - QSI*WW(3*j-1,TT)
        WW(3*j  ,R)= WW(3*j  ,T) - ETA*WW(3*j  ,Y) - QSI*WW(3*j  ,TT)

        WW(3*j-2,T0)= WW(3*j-2,T)
        WW(3*j-1,T0)= WW(3*j-1,T)
        WW(3*j  ,T0)= WW(3*j  ,T)

        DNRM20= DNRM20 + WW(3*j-2,R)**2+ WW(3*j-1,R)**2+ WW(3*j,R)**2
        COEF10= COEF10 + WW(3*j-2,R)*WW(3*j-2,RT)                       &
     &                 + WW(3*j-1,R)*WW(3*j-1,RT) + WW(3*j,R)*WW(3*j,RT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      call HECMW_allREDUCE_DP1 (COEF10, COEF1, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (COEF10, COEF1, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      BETA = ALPHA*COEF1 / (QSI*RHO)
      do j= 1, N
        WW(3*j-2,W1)= WW(3*j-2,TT) + BETA*WW(3*j-2,PT)
        WW(3*j-1,W1)= WW(3*j-1,TT) + BETA*WW(3*j-1,PT)
        WW(3*j  ,W1)= WW(3*j  ,TT) + BETA*WW(3*j  ,PT)
      enddo

      RESID= dsqrt(DNRM2/BNRM2)
      RHO  = COEF1

!C##### ITERATION HISTORY
        if (my_rank.eq.0.and.ITERlog.eq.1) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
!C#####

      if (RESID.le.TOL   ) exit
      if ( ITER.eq.MAXIT ) ERROR= -300
!C===
      enddo
      END_TIME= HECMW_WTIME()

!C
!C-- INTERFACE data EXCHANGE
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW, WR, WS)
      call hecmw_precond_BILU_33_clear(PRECOND)

      Tset = Tset + SETupTIME
      Tsol = Tsol + END_TIME - START_TIME
      Tcomm= COMMtime

      end subroutine  hecmw_solve_BLGPBiCG_33
      end module     hecmw_solver_BLGPBiCG_33

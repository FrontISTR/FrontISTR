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
      real(kind=kreal), dimension(:)  ,  allocatable       :: ALU

      real(kind=kreal), dimension(:), allocatable :: AUlu0, ALlu0
      real(kind=kreal), dimension(:), allocatable :: Dlu0

      integer(kind=kint), dimension(:), allocatable :: inumFI1L, FI1L
      integer(kind=kint), dimension(:), allocatable :: inumFI1U, FI1U

      real(kind=kreal), dimension(5) :: C0, CG
      real(kind=kreal), dimension(2) :: EQ

      integer(kind=kint ) :: R, RT, T, T0, TT, P, PT
      integer(kind=kint ) :: U, W1, Y, Z, WK, W2, MAXIT
      real   (kind=kreal) :: TOL
      integer(kind=kint ) :: i,j,k,isU,ieU,isL,ieL
      integer(kind=kint ) :: ip
      real   (kind=kreal) :: S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal) :: BNRM20,BNRM2,X1,X2,X3,D11,D22,D33, ALUtmp(3,3)
      real   (kind=kreal) :: RHO,RHO0,RHO1,BETA,ALPHA,DNRM20,DNRM2,RHO10
      real   (kind=kreal) :: WVAL1,WVAL2,WVAL3
      real   (kind=kreal) :: SW1,SW2,SW3,aSW1,aSW2,aSW3,bSW1,bSW2,bSW3,cSW1,cSW2,cSW3
      real   (kind=kreal) :: COMMtime,QSI,ETA,COEF10,COEF1

      START_TIME= HECMW_WTIME()
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      allocate (WW(3*NP,13))
      allocate (WS(3*NP), WR(3*NP))
      allocate (ALU(9*NP))

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
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,1),                 &
     &       SOLVER_COMM,my_rank )
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,2),                 &
     &       SOLVER_COMM,my_rank )
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,3),                 &
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
      if (PRECOND.eq.11) call FORM_ILU1_33
      if (PRECOND.eq.12) call FORM_ILU2_33

      if (PRECOND.eq.10) then
        do ip= 1, NP
          D11= D(9*ip-8) * SIGMA_DIAG
          D22= D(9*ip-4) * SIGMA_DIAG
          D33= D(9*ip  ) * SIGMA_DIAG
          call ILU1a33 (ALUtmp,                                    &
     &                  D11      , D(9*ip-7), D(9*ip-6),                &
     &                  D(9*ip-5), D22      , D(9*ip-3),                &
     &                  D(9*ip-2), D(9*ip-1), D33)
          ALU(9*ip-8)= ALUtmp(1,1)
          ALU(9*ip-7)= ALUtmp(1,2)
          ALU(9*ip-6)= ALUtmp(1,3)
          ALU(9*ip-5)= ALUtmp(2,1)
          ALU(9*ip-4)= ALUtmp(2,2)
          ALU(9*ip-3)= ALUtmp(2,3)
          ALU(9*ip-2)= ALUtmp(3,1)
          ALU(9*ip-1)= ALUtmp(3,2)
          ALU(9*ip  )= ALUtmp(3,3)
        enddo
      endif

      if (PRECOND.eq.11.or.PRECOND.eq.12) then
        do ip= 1, NP
          call ILU1a33 (ALUtmp,                                    &
     &                  Dlu0(9*ip-8), Dlu0(9*ip-7), Dlu0(9*ip-6),       &
     &                  Dlu0(9*ip-5), Dlu0(9*ip-4), Dlu0(9*ip-3),       &
     &                  Dlu0(9*ip-2), Dlu0(9*ip-1), Dlu0(9*ip  ))
          ALU(9*ip-8)= ALUtmp(1,1)
          ALU(9*ip-7)= ALUtmp(1,2)
          ALU(9*ip-6)= ALUtmp(1,3)
          ALU(9*ip-5)= ALUtmp(2,1)
          ALU(9*ip-4)= ALUtmp(2,2)
          ALU(9*ip-3)= ALUtmp(2,3)
          ALU(9*ip-2)= ALUtmp(3,1)
          ALU(9*ip-1)= ALUtmp(3,2)
          ALU(9*ip  )= ALUtmp(3,3)
        enddo
      endif
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
!C== Block SSOR
      if (PRECOND.eq.10) then

        do i= 1+N, NP
          WW(3*i-2,R)= 0.d0
          WW(3*i-1,R)= 0.d0
          WW(3*i  ,R)= 0.d0
        enddo
!C
!C-- FORWARD
        do i= 1, N
          SW1= WW(3*i-2,R)
          SW2= WW(3*i-1,R)
          SW3= WW(3*i  ,R)
          isL= INL(i-1)+1
          ieL= INL(i)
          do j= isL, ieL
              k= IAL(j)
             X1= WW(3*k-2,R)
             X2= WW(3*k-1,R)
             X3= WW(3*k  ,R)
            SW1= SW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            SW2= SW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            SW3= SW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
          enddo
          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,R)= X1
          WW(3*i-1,R)= X2
          WW(3*i  ,R)= X3
        enddo
!C
!C-- BACKWARD
        do i= N, 1, -1
          isU= INU(i-1) + 1
          ieU= INU(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= ieU, isU, -1
              k= IAU(j)
             X1= WW(3*k-2,R)
             X2= WW(3*k-1,R)
             X3= WW(3*k  ,R)
            SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
          enddo
          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,R)=  WW(3*i-2,R) - X1
          WW(3*i-1,R)=  WW(3*i-1,R) - X2
          WW(3*i  ,R)=  WW(3*i  ,R) - X3
        enddo

      endif

!C
!C== Block ILU
      if (PRECOND.eq.11.or.PRECOND.eq.12) then
        do i= 1+N, NP
          WW(3*i-2,R)= 0.d0
          WW(3*i-1,R)= 0.d0
          WW(3*i  ,R)= 0.d0
        enddo

!C
!C-- FORWARD
        do i= 1, N
          SW1= WW(3*i-2,R)
          SW2= WW(3*i-1,R)
          SW3= WW(3*i  ,R)
          isL= inumFI1L(i-1)+1
          ieL= inumFI1L(i)
          do j= isL, ieL
              k= FI1L(j)
             X1= WW(3*k-2,R)
             X2= WW(3*k-1,R)
             X3= WW(3*k  ,R)
            SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
            SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
            SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
          enddo

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,R)=  X1
          WW(3*i-1,R)=  X2
          WW(3*i  ,R)=  X3
        enddo
      
!C
!C-- BACKWARD
        do i= N, 1, -1
          isU= inumFI1U(i-1) + 1
          ieU= inumFI1U(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= ieU, isU, -1
              k= FI1U(j)
             X1= WW(3*k-2,R)
             X2= WW(3*k-1,R)
             X3= WW(3*k  ,R)
            SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
            SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
            SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          enddo

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,R)=  WW(3*i-2,R) - X1
          WW(3*i-1,R)=  WW(3*i-1,R) - X2
          WW(3*i  ,R)=  WW(3*i  ,R) - X3
        enddo
      endif
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
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,P) , SOLVER_COMM,     &
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
!C== Block SSOR
      if (PRECOND.eq.10) then
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

!C
!C-- FORWARD
        do i= 1, N
          aSW1= WW(3*i-2,TT)
          aSW2= WW(3*i-1,TT)
          aSW3= WW(3*i  ,TT)
          bSW1= WW(3*i-2,W2)
          bSW2= WW(3*i-1,W2)
          bSW3= WW(3*i  ,W2)
          cSW1= WW(3*i-2,T0)
          cSW2= WW(3*i-1,T0)
          cSW3= WW(3*i  ,T0)
          isL= INL(i-1)+1
          ieL= INL(i)
          do j= isL, ieL
               k= IAL(j)
              X1= WW(3*k-2,TT)
              X2= WW(3*k-1,TT)
              X3= WW(3*k  ,TT)
            aSW1= aSW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            aSW2= aSW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            aSW3= aSW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
              X1= WW(3*k-2,W2)
              X2= WW(3*k-1,W2)
              X3= WW(3*k  ,W2)
            bSW1= bSW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            bSW2= bSW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            bSW3= bSW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
              X1= WW(3*k-2,T0)
              X2= WW(3*k-1,T0)
              X3= WW(3*k  ,T0)
            cSW1= cSW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            cSW2= cSW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            cSW3= cSW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
          enddo
          X1= aSW1
          X2= aSW2
          X3= aSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,TT)= X1
          WW(3*i-1,TT)= X2
          WW(3*i  ,TT)= X3

          X1= bSW1
          X2= bSW2
          X3= bSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,W2)= X1
          WW(3*i-1,W2)= X2
          WW(3*i  ,W2)= X3

          X1= cSW1
          X2= cSW2
          X3= cSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,T0)= X1
          WW(3*i-1,T0)= X2
          WW(3*i  ,T0)= X3
        enddo
!C
!C-- BACKWARD
        do i= N, 1, -1
           isU= INU(i-1) + 1
           ieU= INU(i) 
          aSW1= 0.d0
          aSW2= 0.d0
          aSW3= 0.d0
          bSW1= 0.d0
          bSW2= 0.d0
          bSW3= 0.d0
          cSW1= 0.d0
          cSW2= 0.d0
          cSW3= 0.d0
          do j= ieU, isU, -1
               k= IAU(j)
              X1= WW(3*k-2,TT)
              X2= WW(3*k-1,TT)
              X3= WW(3*k  ,TT)
            aSW1= aSW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            aSW2= aSW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            aSW3= aSW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
              X1= WW(3*k-2,W2)
              X2= WW(3*k-1,W2)
              X3= WW(3*k  ,W2)
            bSW1= bSW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            bSW2= bSW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            bSW3= bSW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
              X1= WW(3*k-2,T0)
              X2= WW(3*k-1,T0)
              X3= WW(3*k  ,T0)
            cSW1= cSW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            cSW2= cSW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            cSW3= cSW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
          enddo

          X1= aSW1
          X2= aSW2
          X3= aSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,TT)=  WW(3*i-2,TT) - X1
          WW(3*i-1,TT)=  WW(3*i-1,TT) - X2
          WW(3*i  ,TT)=  WW(3*i  ,TT) - X3

          X1= bSW1
          X2= bSW2
          X3= bSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,W2)=  WW(3*i-2,W2) - X1
          WW(3*i-1,W2)=  WW(3*i-1,W2) - X2
          WW(3*i  ,W2)=  WW(3*i  ,W2) - X3

          X1= cSW1
          X2= cSW2
          X3= cSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,T0)=  WW(3*i-2,T0) - X1
          WW(3*i-1,T0)=  WW(3*i-1,T0) - X2
          WW(3*i  ,T0)=  WW(3*i  ,T0) - X3
        enddo

      endif

!C
!C== Block ILU
      if (PRECOND.eq.11.or.PRECOND.eq.12) then
      do i= 1, N
        WW(3*i-2,TT)= WW(3*i-2,T )
        WW(3*i-1,TT)= WW(3*i-1,T )
        WW(3*i  ,TT)= WW(3*i  ,T )
        WW(3*i-2,W2)= WW(3*i-2,PT)
        WW(3*i-1,W2)= WW(3*i-1,PT)
        WW(3*i  ,W2)= WW(3*i  ,PT)
      enddo

      do i= 1+N, NP
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

!C
!C-- FORWARD
        do i= 1, N
          aSW1= WW(3*i-2,TT)
          aSW2= WW(3*i-1,TT)
          aSW3= WW(3*i  ,TT)
          bSW1= WW(3*i-2,W2)
          bSW2= WW(3*i-1,W2)
          bSW3= WW(3*i  ,W2)
          cSW1= WW(3*i-2,T0)
          cSW2= WW(3*i-1,T0)
          cSW3= WW(3*i  ,T0)
          isL= inumFI1L(i-1)+1
          ieL= inumFI1L(i)

          do j= isL, ieL
               k= FI1L(j)
              X1= WW(3*k-2,TT)
              X2= WW(3*k-1,TT)
              X3= WW(3*k  ,TT)
            aSW1= aSW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
            aSW2= aSW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
            aSW3= aSW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
              X1= WW(3*k-2,W2)
              X2= WW(3*k-1,W2)
              X3= WW(3*k  ,W2)
            bSW1= bSW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
            bSW2= bSW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
            bSW3= bSW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
              X1= WW(3*k-2,T0)
              X2= WW(3*k-1,T0)
              X3= WW(3*k  ,T0)
            cSW1= cSW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
            cSW2= cSW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
            cSW3= cSW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
          enddo

          X1= aSW1
          X2= aSW2
          X3= aSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,TT)=  X1
          WW(3*i-1,TT)=  X2
          WW(3*i  ,TT)=  X3

          X1= bSW1
          X2= bSW2
          X3= bSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,W2)=  X1
          WW(3*i-1,W2)=  X2
          WW(3*i  ,W2)=  X3

          X1= cSW1
          X2= cSW2
          X3= cSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,T0)=  X1
          WW(3*i-1,T0)=  X2
          WW(3*i  ,T0)=  X3
        enddo
      
!C
!C-- BACKWARD
        do i= N, 1, -1
           isU= inumFI1U(i-1) + 1
           ieU= inumFI1U(i) 
          aSW1= 0.d0
          aSW2= 0.d0
          aSW3= 0.d0
          bSW1= 0.d0
          bSW2= 0.d0
          bSW3= 0.d0
          cSW1= 0.d0
          cSW2= 0.d0
          cSW3= 0.d0

          do j= ieU, isU, -1
               k= FI1U(j)
              X1= WW(3*k-2,TT)
              X2= WW(3*k-1,TT)
              X3= WW(3*k  ,TT)
            aSW1= aSW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
            aSW2= aSW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
            aSW3= aSW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
              X1= WW(3*k-2,W2)
              X2= WW(3*k-1,W2)
              X3= WW(3*k  ,W2)
            bSW1= bSW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
            bSW2= bSW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
            bSW3= bSW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
              X1= WW(3*k-2,T0)
              X2= WW(3*k-1,T0)
              X3= WW(3*k  ,T0)
            cSW1= cSW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
            cSW2= cSW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
            cSW3= cSW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          enddo

          X1= aSW1
          X2= aSW2
          X3= aSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,TT)=  WW(3*i-2,TT) - X1
          WW(3*i-1,TT)=  WW(3*i-1,TT) - X2
          WW(3*i  ,TT)=  WW(3*i  ,TT) - X3

          X1= bSW1
          X2= bSW2
          X3= bSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,W2)=  WW(3*i-2,W2) - X1
          WW(3*i-1,W2)=  WW(3*i-1,W2) - X2
          WW(3*i  ,W2)=  WW(3*i  ,W2) - X3

          X1= cSW1
          X2= cSW2
          X3= cSW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,T0)=  WW(3*i-2,T0) - X1
          WW(3*i-1,T0)=  WW(3*i-1,T0) - X2
          WW(3*i  ,T0)=  WW(3*i  ,T0) - X3
        enddo
      endif
!C==
!C===

!C
!C-- calc. [A]{t_tld}
      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,TT) , SOLVER_COMM,    &
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

      deallocate (WW, WR, WS, ALU)
      if (PRECOND.ne.10) then
        deallocate (AUlu0, ALlu0, Dlu0)
        deallocate (inumFI1L, inumFI1U, FI1L, FI1U)
      endif

      Tset = Tset + SETupTIME
      Tsol = Tsol + END_TIME - START_TIME
      Tcomm= COMMtime

      contains


!C
!C***
!C*** FORM_ILU1_33
!C***
!C
!C    form ILU(1) matrix
!C
      subroutine FORM_ILU1_33
      implicit none
      integer(kind=kint), dimension(:), allocatable :: IW1 , IW2
      integer(kind=kint), dimension(:), allocatable :: IWsL, IWsU
      real (kind=kreal),  dimension(3,3) :: RHS_Aij, DkINV, Aik, Akj
      real (kind=kreal)  :: D12,D13,D21,D23,D31,D32
      integer(kind=kint) :: NPLf1,NPLf2,NPUf1
      integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
      integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3

!C
!C +--------------+
!C | find fill-in |
!C +--------------+
!C===

!C
!C-- count fill-in
      allocate (IW1(NP) , IW2(NP))
      allocate (inumFI1L(0:NP), inumFI1U(0:NP))

      inumFI1L= 0
      inumFI1U= 0

      NPLf1= 0
      NPUf1= 0
      do i= 2, NP
        icou= 0
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 1
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              inumFI1L(i)= inumFI1L(i)+1
                  IW1(jj)= 1
            endif  
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              inumFI1U(i)= inumFI1U(i)+1
                  IW1(jj)= 1
            endif  
          enddo
        enddo
        NPLf1= NPLf1 + inumFI1L(i)
        NPUf1= NPUf1 + inumFI1U(i)
      enddo

!C
!C-- specify fill-in
      allocate (IWsL(0:NP), IWsU(0:NP))
      allocate (FI1L (NPL+NPLf1), FI1U (NPU+NPUf1))
      allocate (ALlu0(9*(NPL+NPLf1)), AUlu0(9*(NPU+NPUf1)))

      FI1L= 0
      FI1U= 0

      IWsL= 0
      IWsU= 0
      do i= 1, NP
        IWsL(i)= INL(i)-INL(i-1) + inumFI1L(i) + IWsL(i-1)
        IWsU(i)= INU(i)-INU(i-1) + inumFI1U(i) + IWsU(i-1)
      enddo

      do i= 2, NP
        icouL= 0
        icouU= 0
        inumFI1L(i)= inumFI1L(i-1) + inumFI1L(i)
        inumFI1U(i)= inumFI1U(i-1) + inumFI1U(i)
        icou= 0
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 1
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
                   icouL           = icouL + 1
              FI1L(icouL+IWsL(i-1)+INL(i)-INL(i-1))= jj
                  IW1(jj)          = 1
            endif
            if (IW1(jj).eq.0 .and. jj.gt.i) then
                   icouU           = icouU + 1
              FI1U(icouU+IWsU(i-1)+INU(i)-INU(i-1))= jj
                  IW1(jj)          = 1
            endif
          enddo
        enddo
      enddo
!C=== 

!C
!C +-------------------------------------------------+
!C | SORT and RECONSTRUCT matrix considering fill-in |
!C +-------------------------------------------------+
!C===
      ALlu0= 0.d0
      AUlu0= 0.d0
      iSL  = 0
      iSU  = 0
      do i= 1, NP
        icouL1=      INL(i) -      INL(i-1)
        icouL2= inumFI1L(i) - inumFI1L(i-1)
        icouL3= icouL1 + icouL2
        icouU1=      INU(i) -      INU(i-1)
        icouU2= inumFI1U(i) - inumFI1U(i-1)
        icouU3= icouU1 + icouU2
!C
!C-- LOWER part
        icou0= 0
        do k= INL(i-1)+1, INL(i)
              icou0 = icou0 + 1
          IW1(icou0)= IAL(k)
        enddo

        do k= inumFI1L(i-1)+1, inumFI1L(i)
              icou0 = icou0 + 1
          IW1(icou0)= FI1L(icou0+IWsL(i-1))
        enddo

        do k= 1, icouL3
          IW2(k)= k
        enddo
        call fill_in_S33_SORT (IW1, IW2, icouL3, NP)

        do k= 1, icouL3
          FI1L (k+isL)= IW1(k)
          ik= IW2(k)
          if (ik.le.INL(i)-INL(i-1)) then
            kk1= 9*( k+isL)
            kk2= 9*(ik+INL(i-1))
            ALlu0(kk1-8)= AL(kk2-8)
            ALlu0(kk1-7)= AL(kk2-7)
            ALlu0(kk1-6)= AL(kk2-6)
            ALlu0(kk1-5)= AL(kk2-5)
            ALlu0(kk1-4)= AL(kk2-4)
            ALlu0(kk1-3)= AL(kk2-3)
            ALlu0(kk1-2)= AL(kk2-2)
            ALlu0(kk1-1)= AL(kk2-1)
            ALlu0(kk1  )= AL(kk2  )
          endif
        enddo
!C
!C-- UPPER part
        icou0= 0
        do k= INU(i-1)+1, INU(i)
              icou0 = icou0 + 1
          IW1(icou0)= IAU(k)
        enddo

        do k= inumFI1U(i-1)+1, inumFI1U(i)
              icou0 = icou0 + 1
          IW1(icou0)= FI1U(icou0+IWsU(i-1))
        enddo

        do k= 1, icouU3
          IW2(k)= k
        enddo
        call fill_in_S33_SORT (IW1, IW2, icouU3, NP)

        do k= 1, icouU3
          FI1U (k+isU)= IW1(k)
          ik= IW2(k)
          if (ik.le.INU(i)-INU(i-1)) then
            kk1= 9*( k+isU)
            kk2= 9*(ik+INU(i-1))
            AUlu0(kk1-8)= AU(kk2-8)
            AUlu0(kk1-7)= AU(kk2-7)
            AUlu0(kk1-6)= AU(kk2-6)
            AUlu0(kk1-5)= AU(kk2-5)
            AUlu0(kk1-4)= AU(kk2-4)
            AUlu0(kk1-3)= AU(kk2-3)
            AUlu0(kk1-2)= AU(kk2-2)
            AUlu0(kk1-1)= AU(kk2-1)
            AUlu0(kk1  )= AU(kk2  )
          endif
        enddo

        iSL= iSL + icouL3
        iSU= iSU + icouU3
      enddo

!C===
      do i= 1, NP
        inumFI1L(i)= IWsL(i)
        inumFI1U(i)= IWsU(i)
      enddo
      deallocate (IWsL, IWsU)

!C
!C +----------------------+
!C | ILU(1) factorization |
!C +----------------------+
!C===
      allocate (Dlu0(9*NP))
      Dlu0= D
      do i=1,NP
        Dlu0(9*k-8)=Dlu0(9*k-8)*SIGMA_DIAG
        Dlu0(9*k-4)=Dlu0(9*k-4)*SIGMA_DIAG
        Dlu0(9*k  )=Dlu0(9*k  )*SIGMA_DIAG
      enddo

      do i= 2, NP
        IW1= 0
        IW2= 0

        do k= inumFI1L(i-1)+1, inumFI1L(i)
          IW1(FI1L(k))= k
        enddo

        do k= inumFI1U(i-1)+1, inumFI1U(i)
          IW2(FI1U(k))= k
        enddo

        do kk= INL(i-1)+1, INL(i)
          k= IAL(kk)
          D11= Dlu0(9*k-8)
          D12= Dlu0(9*k-7)
          D13= Dlu0(9*k-6)
          D21= Dlu0(9*k-5)
          D22= Dlu0(9*k-4)
          D23= Dlu0(9*k-3)
          D31= Dlu0(9*k-2)
          D32= Dlu0(9*k-1)
          D33= Dlu0(9*k  )

          call ILU1a33 (DkINV, D11,D12,D13,D21,D22,D23,D31,D32,D33)

          do kk1= inumFI1L(i-1)+1, inumFI1L(i)
            if (k.eq.FI1L(kk1)) then
              Aik(1,1)= ALlu0(9*kk1-8)
              Aik(1,2)= ALlu0(9*kk1-7)
              Aik(1,3)= ALlu0(9*kk1-6)
              Aik(2,1)= ALlu0(9*kk1-5)
              Aik(2,2)= ALlu0(9*kk1-4)
              Aik(2,3)= ALlu0(9*kk1-3)
              Aik(3,1)= ALlu0(9*kk1-2)
              Aik(3,2)= ALlu0(9*kk1-1)
              Aik(3,3)= ALlu0(9*kk1  )
              exit
            endif
          enddo
  
          do jj= INU(k-1)+1, INU(k)
            j= IAU(jj)
            do jj1= inumFI1U(k-1)+1, inumFI1U(k)
              if (j.eq.FI1U(jj1)) then
                Akj(1,1)= AUlu0(9*jj1-8)
                Akj(1,2)= AUlu0(9*jj1-7)
                Akj(1,3)= AUlu0(9*jj1-6)
                Akj(2,1)= AUlu0(9*jj1-5)
                Akj(2,2)= AUlu0(9*jj1-4)
                Akj(2,3)= AUlu0(9*jj1-3)
                Akj(3,1)= AUlu0(9*jj1-2)
                Akj(3,2)= AUlu0(9*jj1-1)
                Akj(3,3)= AUlu0(9*jj1  )
                exit
              endif
            enddo

            call ILU1b33 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(9*i-8)= Dlu0(9*i-8) - RHS_Aij(1,1)
              Dlu0(9*i-7)= Dlu0(9*i-7) - RHS_Aij(1,2)
              Dlu0(9*i-6)= Dlu0(9*i-6) - RHS_Aij(1,3)
              Dlu0(9*i-5)= Dlu0(9*i-5) - RHS_Aij(2,1)
              Dlu0(9*i-4)= Dlu0(9*i-4) - RHS_Aij(2,2)
              Dlu0(9*i-3)= Dlu0(9*i-3) - RHS_Aij(2,3)
              Dlu0(9*i-2)= Dlu0(9*i-2) - RHS_Aij(3,1)
              Dlu0(9*i-1)= Dlu0(9*i-1) - RHS_Aij(3,2)
              Dlu0(9*i  )= Dlu0(9*i  ) - RHS_Aij(3,3)
            endif   

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(9*ij0-8)= ALlu0(9*ij0-8) - RHS_Aij(1,1)
              ALlu0(9*ij0-7)= ALlu0(9*ij0-7) - RHS_Aij(1,2)
              ALlu0(9*ij0-6)= ALlu0(9*ij0-6) - RHS_Aij(1,3)
              ALlu0(9*ij0-5)= ALlu0(9*ij0-5) - RHS_Aij(2,1)
              ALlu0(9*ij0-4)= ALlu0(9*ij0-4) - RHS_Aij(2,2)
              ALlu0(9*ij0-3)= ALlu0(9*ij0-3) - RHS_Aij(2,3)
              ALlu0(9*ij0-2)= ALlu0(9*ij0-2) - RHS_Aij(3,1)
              ALlu0(9*ij0-1)= ALlu0(9*ij0-1) - RHS_Aij(3,2)
              ALlu0(9*ij0  )= ALlu0(9*ij0  ) - RHS_Aij(3,3)
            endif   

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(9*ij0-8)= AUlu0(9*ij0-8) - RHS_Aij(1,1)
              AUlu0(9*ij0-7)= AUlu0(9*ij0-7) - RHS_Aij(1,2)
              AUlu0(9*ij0-6)= AUlu0(9*ij0-6) - RHS_Aij(1,3)
              AUlu0(9*ij0-5)= AUlu0(9*ij0-5) - RHS_Aij(2,1)
              AUlu0(9*ij0-4)= AUlu0(9*ij0-4) - RHS_Aij(2,2)
              AUlu0(9*ij0-3)= AUlu0(9*ij0-3) - RHS_Aij(2,3)
              AUlu0(9*ij0-2)= AUlu0(9*ij0-2) - RHS_Aij(3,1)
              AUlu0(9*ij0-1)= AUlu0(9*ij0-1) - RHS_Aij(3,2)
              AUlu0(9*ij0  )= AUlu0(9*ij0  ) - RHS_Aij(3,3)
            endif   

          enddo
	enddo
      enddo

      deallocate (IW1, IW2)
!C===
      end subroutine FORM_ILU1_33

!C
!C***
!C*** FORM_ILU2_33
!C***
!C
!C    form ILU(2) matrix
!C
      subroutine FORM_ILU2_33
      implicit none
      integer(kind=kint), dimension(:), allocatable:: IW1 , IW2
      integer(kind=kint), dimension(:), allocatable:: IWsL, IWsU
      integer(kind=kint), dimension(:), allocatable:: iconFI1L, iconFI1U
      integer(kind=kint), dimension(:), allocatable:: inumFI2L, inumFI2U
      integer(kind=kint), dimension(:), allocatable::     FI2L,     FI2U
      real (kind=kreal), dimension(3,3) :: RHS_Aij, DkINV, Aik, Akj
      real (kind=kreal)  :: D12,D13,D21,D23,D31,D32
      integer(kind=kint) :: NPLf1,NPLf2,NPUf1,NPUF2
      integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj,iAS
      integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3
      integer(kind=kint) :: icouL,icouL1,icouL2,icouL3,iconIK,iconKJ

!C
!C +------------------+  
!C | find fill-in (1) |
!C +------------------+
!C===

!C
!C-- count fill-in
      allocate (IW1(NP) , IW2(NP))
      allocate (inumFI2L(0:NP), inumFI2U(0:NP))

      inumFI2L= 0
      inumFI2U= 0

      NPLf1= 0
      NPUf1= 0
      do i= 2, NP
        icou= 0
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 1
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              inumFI2L(i)= inumFI2L(i)+1
                  IW1(jj)= 1
            endif  
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              inumFI2U(i)= inumFI2U(i)+1
                  IW1(jj)= 1
            endif  
          enddo
        enddo
        NPLf1= NPLf1 + inumFI2L(i)
        NPUf1= NPUf1 + inumFI2U(i)
      enddo

!C
!C-- specify fill-in
      allocate (IWsL(0:NP), IWsU(0:NP))
      allocate (FI2L (NPLf1), FI2U (NPUf1))

      FI2L= 0
      FI2U= 0

      do i= 2, NP
        icouL= 0
        icouU= 0
        inumFI2L(i)= inumFI2L(i-1) + inumFI2L(i)
        inumFI2U(i)= inumFI2U(i-1) + inumFI2U(i)
        icou= 0
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 1
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
                   icouL = icouL + 1
              FI2L(icouL+inumFI2L(i-1))= jj
                  IW1(jj)= 1
            endif
            if (IW1(jj).eq.0 .and. jj.gt.i) then
                   icouU = icouU + 1
              FI2U(icouU+inumFI2U(i-1))= jj
                  IW1(jj)= 1
            endif
          enddo
        enddo
      enddo
!C=== 

!C
!C +------------------+  
!C | find fill-in (2) |
!C +------------------+
!C===
      allocate (inumFI1L(0:NP), inumFI1U(0:NP))

      NPLf2= 0
      NPUf2= 0
      inumFI1L= 0
      inumFI1U= 0
!C
!C-- count fill-in
      do i= 2, NP
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 2
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 2
        enddo

        do L= inumFI2L(i-1)+1, inumFI2L(i)
          IW1(FI2L(L))= 1
        enddo

        do L= inumFI2U(i-1)+1, inumFI2U(i)
          IW1(FI2U(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= inumFI2U(kk-1) + 1
          iEj= inumFI2U(kk)
          do j= iSj, iEj
            jj= FI2U(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              inumFI1L(i)= inumFI1L(i) + 1
                  IW1(jj)= 1
            endif  
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              inumFI1U(i)= inumFI1U(i) + 1
                  IW1(jj)= 1
            endif  
          enddo
        enddo

        iSk= inumFI2L(i-1)+1
        iEk= inumFI2L(i)
        do k= iSk, iEk
          kk= FI2L(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              inumFI1L(i)= inumFI1L(i) + 1
                  IW1(jj)= 1
            endif  
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              inumFI1U(i)= inumFI1U(i) + 1
                  IW1(jj)= 1
            endif  
          enddo
        enddo
        NPLf2= NPLf2 + inumFI1L(i)
        NPUf2= NPUf2 + inumFI1U(i)
      enddo

!C
!C-- specify fill-in
      allocate (FI1L(NPL+NPLf1+NPLf2))
      allocate (FI1U(NPU+NPUf1+NPUf2))

      allocate (iconFI1L(NPL+NPLf1+NPLf2))
      allocate (iconFI1U(NPU+NPUf1+NPUf2))

      IWsL= 0
      IWsU= 0
      do i= 1, NP
        IWsL(i)= INL(i)-INL(i-1) + inumFI2L(i)-inumFI2L(i-1) +          &
     &                             inumFI1L(i) + IWsL(i-1)
        IWsU(i)= INU(i)-INU(i-1) + inumFI2U(i)-inumFI2U(i-1) +          &
     &                             inumFI1U(i) + IWsU(i-1)
      enddo

      do i= 2, NP
        icouL= 0
        icouU= 0
        inumFI1L(i)= inumFI1L(i-1) + inumFI1L(i)
        inumFI1U(i)= inumFI1U(i-1) + inumFI1U(i)
        icou= 0
        IW1= 0
        IW1(i)= 1
        do L= INL(i-1)+1, INL(i)
          IW1(IAL(L))= 1
        enddo
        do L= INU(i-1)+1, INU(i)
          IW1(IAU(L))= 1
        enddo

        do L= inumFI2L(i-1)+1, inumFI2L(i)
          IW1(FI2L(L))= 1
        enddo

        do L= inumFI2U(i-1)+1, inumFI2U(i)
          IW1(FI2U(L))= 1
        enddo

        iSk= INL(i-1) + 1
        iEk= INL(i)
        do k= iSk, iEk
          kk= IAL(k)
          iSj= inumFI2U(kk-1) + 1
          iEj= inumFI2U(kk  )
          do j= iSj, iEj
            jj= FI2U(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              iAS= INL(i)-INL(i-1)+inumFI2L(i)-inumFI2L(i-1)+IWsL(i-1)
                   icouL     = icouL + 1
              FI1L(icouL+iAS)= jj
                  IW1(jj)    = 1
            endif
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              iAS= INU(i)-INU(i-1)+inumFI2U(i)-inumFI2U(i-1)+IWsU(i-1)
                   icouU     = icouU + 1
              FI1U(icouU+iAS)= jj
                  IW1(jj)    = 1
            endif
          enddo
        enddo

        iSk= inumFI2L(i-1) + 1
        iEk= inumFI2L(i)
        do k= iSk, iEk
          kk= FI2L(k)
          iSj= INU(kk-1) + 1
          iEj= INU(kk  )
          do j= iSj, iEj
            jj= IAU(j)
            if (IW1(jj).eq.0 .and. jj.lt.i) then
              iAS= INL(i)-INL(i-1)+inumFI2L(i)-inumFI2L(i-1)+IWsL(i-1)
                   icouL     = icouL + 1
              FI1L(icouL+iAS)= jj
                  IW1(jj)    = 1
            endif
            if (IW1(jj).eq.0 .and. jj.gt.i) then
              iAS= INU(i)-INU(i-1)+inumFI2U(i)-inumFI2U(i-1)+IWsU(i-1)
                   icouU     = icouU + 1
              FI1U(icouU+iAS)= jj
                  IW1(jj)    = 1
            endif
          enddo
        enddo
      enddo
!C=== 

!C
!C +-------------------------------------------------+
!C | SORT and RECONSTRUCT matrix considering fill-in |
!C +-------------------------------------------------+
!C===
      allocate (ALlu0(9*(NPL+NPLf1+NPLf2)))
      allocate (AUlu0(9*(NPU+NPUf1+NPUf2)))

      ALlu0= 0.d0
      AUlu0= 0.d0
      iSL  = 0
      iSU  = 0

      iconFI1L= 0
      iconFI1U= 0

      do i= 1, NP
        icouL1=      INL(i) -      INL(i-1)
        icouL2= inumFI2L(i) - inumFI2L(i-1) + icouL1
        icouL3= inumFI1L(i) - inumFI1L(i-1) + icouL2

        icouU1=      INU(i) -      INU(i-1)
        icouU2= inumFI2U(i) - inumFI2U(i-1) + icouU1
        icouU3= inumFI1U(i) - inumFI1U(i-1) + icouU2

!C
!C-- LOWER part
        icou= 0
        do k= INL(i-1)+1, INL(i)
              icou = icou + 1
          IW1(icou)= IAL(k)
        enddo

        icou= 0
        do k= inumFI2L(i-1)+1, inumFI2L(i)
              icou        = icou + 1
          IW1(icou+icouL1)= FI2L(k)
        enddo

        icou= 0
        do k= inumFI1L(i-1)+1, inumFI1L(i)
               icou        = icou + 1
          IW1(icou+icouL2)= FI1L(icou+icouL2+iSL)
        enddo

        do k= 1, icouL3
          IW2(k)= k
        enddo

        call fill_in_S33_SORT (IW1, IW2, icouL3, NP)

        do k= 1, icouL3
          FI1L (k+isL)= IW1(k)
          ik= IW2(k)
          if (ik.le.INL(i)-INL(i-1)) then
            kk1= 9*( k+isL)
            kk2= 9*(ik+INL(i-1))
            ALlu0(kk1-8)= AL(kk2-8)
            ALlu0(kk1-7)= AL(kk2-7)
            ALlu0(kk1-6)= AL(kk2-6)
            ALlu0(kk1-5)= AL(kk2-5)
            ALlu0(kk1-4)= AL(kk2-4)
            ALlu0(kk1-3)= AL(kk2-3)
            ALlu0(kk1-2)= AL(kk2-2)
            ALlu0(kk1-1)= AL(kk2-1)
            ALlu0(kk1  )= AL(kk2  )
          endif
        enddo

        icou= 0
        do k= INL(i-1)+1, INL(i)
              icou = icou + 1
          IW1(icou)= 0
        enddo

        icou= 0
        do k= inumFI2L(i-1)+1, inumFI2L(i)
              icou        = icou + 1
          IW1(icou+icouL1)= 1
        enddo

        icou= 0
        do k= inumFI1L(i-1)+1, inumFI1L(i)
               icou        = icou + 1
          IW1(icou+icouL2)= 2
        enddo

        do k= 1, icouL3
          iconFI1L(k+iSL)= IW1(IW2(k))
        enddo
!C
!C-- UPPER part
        icou= 0
        do k= INU(i-1)+1, INU(i)
              icou = icou + 1
          IW1(icou)= IAU(k)
        enddo

        icou= 0
        do k= inumFI2U(i-1)+1, inumFI2U(i)
              icou        = icou + 1
          IW1(icou+icouU1)= FI2U(k)
        enddo

        icou= 0
        do k= inumFI1U(i-1)+1, inumFI1U(i)
               icou        = icou + 1
          IW1(icou+icouU2)= FI1U(icou+icouU2+iSU)
        enddo

        do k= 1, icouU3
          IW2(k)= k
        enddo
        call fill_in_S33_SORT (IW1, IW2, icouU3, NP)

        do k= 1, icouU3
          FI1U (k+isU)= IW1(k)
          ik= IW2(k)
          if (ik.le.INU(i)-INU(i-1)) then
            kk1= 9*( k+isU)
            kk2= 9*(ik+INU(i-1))
            AUlu0(kk1-8)= AU(kk2-8)
            AUlu0(kk1-7)= AU(kk2-7)
            AUlu0(kk1-6)= AU(kk2-6)
            AUlu0(kk1-5)= AU(kk2-5)
            AUlu0(kk1-4)= AU(kk2-4)
            AUlu0(kk1-3)= AU(kk2-3)
            AUlu0(kk1-2)= AU(kk2-2)
            AUlu0(kk1-1)= AU(kk2-1)
            AUlu0(kk1  )= AU(kk2  )
          endif
        enddo

        icou= 0
        do k= INU(i-1)+1, INU(i)
              icou = icou + 1
          IW1(icou)= 0
        enddo

        icou= 0
        do k= inumFI2U(i-1)+1, inumFI2U(i)
              icou        = icou + 1
          IW1(icou+icouU1)= 1
        enddo

        icou= 0
        do k= inumFI1U(i-1)+1, inumFI1U(i)
               icou        = icou + 1
          IW1(icou+icouU2)= 2
        enddo

        do k= 1, icouU3
          iconFI1U(k+iSU)= IW1(IW2(k))
        enddo

        iSL= iSL + icouL3
        iSU= iSU + icouU3
      enddo
!C===
      do i= 1, NP
        inumFI1L(i)= IWsL(i)
        inumFI1U(i)= IWsU(i)
      enddo

      deallocate (IWsL, IWsU)
      deallocate (inumFI2L, inumFI2U)
      deallocate (    FI2L,     FI2U)

!C
!C +----------------------+
!C | ILU(2) factorization |
!C +----------------------+
!C===
      allocate (Dlu0(9*NP))
      Dlu0= D
      do i=1,NP
        Dlu0(9*k-8)=Dlu0(9*k-8)*SIGMA_DIAG
        Dlu0(9*k-4)=Dlu0(9*k-4)*SIGMA_DIAG
        Dlu0(9*k  )=Dlu0(9*k  )*SIGMA_DIAG
      enddo

      do i= 2, NP
        IW1= 0
        IW2= 0

        do k= inumFI1L(i-1)+1, inumFI1L(i)
          IW1(FI1L(k))= k
        enddo

        do k= inumFI1U(i-1)+1, inumFI1U(i)
          IW2(FI1U(k))= k
        enddo

        do kk= inumFI1L(i-1)+1, inumFI1L(i)
          k= FI1L(kk)
          iconIK= iconFI1L(kk)

          D11= Dlu0(9*k-8)
          D12= Dlu0(9*k-7)
          D13= Dlu0(9*k-6)
          D21= Dlu0(9*k-5)
          D22= Dlu0(9*k-4)
          D23= Dlu0(9*k-3)
          D31= Dlu0(9*k-2)
          D32= Dlu0(9*k-1)
          D33= Dlu0(9*k  )

          call ILU1a33 (DkINV, D11,D12,D13,D21,D22,D23,D31,D32,D33)

          Aik(1,1)= ALlu0(9*kk-8)
          Aik(1,2)= ALlu0(9*kk-7)
          Aik(1,3)= ALlu0(9*kk-6)
          Aik(2,1)= ALlu0(9*kk-5)
          Aik(2,2)= ALlu0(9*kk-4)
          Aik(2,3)= ALlu0(9*kk-3)
          Aik(3,1)= ALlu0(9*kk-2)
          Aik(3,2)= ALlu0(9*kk-1)
          Aik(3,3)= ALlu0(9*kk  )
  
          do jj= inumFI1U(k-1)+1, inumFI1U(k)
            j= FI1U(jj)
            iconKJ= iconFI1U(jj)

            if ((iconIK+iconKJ).lt.2) then
            Akj(1,1)= AUlu0(9*jj-8)
            Akj(1,2)= AUlu0(9*jj-7)
            Akj(1,3)= AUlu0(9*jj-6)
            Akj(2,1)= AUlu0(9*jj-5)
            Akj(2,2)= AUlu0(9*jj-4)
            Akj(2,3)= AUlu0(9*jj-3)
            Akj(3,1)= AUlu0(9*jj-2)
            Akj(3,2)= AUlu0(9*jj-1)
            Akj(3,3)= AUlu0(9*jj  )

            call ILU1b33 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(9*i-8)= Dlu0(9*i-8) - RHS_Aij(1,1)
              Dlu0(9*i-7)= Dlu0(9*i-7) - RHS_Aij(1,2)
              Dlu0(9*i-6)= Dlu0(9*i-6) - RHS_Aij(1,3)
              Dlu0(9*i-5)= Dlu0(9*i-5) - RHS_Aij(2,1)
              Dlu0(9*i-4)= Dlu0(9*i-4) - RHS_Aij(2,2)
              Dlu0(9*i-3)= Dlu0(9*i-3) - RHS_Aij(2,3)
              Dlu0(9*i-2)= Dlu0(9*i-2) - RHS_Aij(3,1)
              Dlu0(9*i-1)= Dlu0(9*i-1) - RHS_Aij(3,2)
              Dlu0(9*i  )= Dlu0(9*i  ) - RHS_Aij(3,3)
            endif   

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(9*ij0-8)= ALlu0(9*ij0-8) - RHS_Aij(1,1)
              ALlu0(9*ij0-7)= ALlu0(9*ij0-7) - RHS_Aij(1,2)
              ALlu0(9*ij0-6)= ALlu0(9*ij0-6) - RHS_Aij(1,3)
              ALlu0(9*ij0-5)= ALlu0(9*ij0-5) - RHS_Aij(2,1)
              ALlu0(9*ij0-4)= ALlu0(9*ij0-4) - RHS_Aij(2,2)
              ALlu0(9*ij0-3)= ALlu0(9*ij0-3) - RHS_Aij(2,3)
              ALlu0(9*ij0-2)= ALlu0(9*ij0-2) - RHS_Aij(3,1)
              ALlu0(9*ij0-1)= ALlu0(9*ij0-1) - RHS_Aij(3,2)
              ALlu0(9*ij0  )= ALlu0(9*ij0  ) - RHS_Aij(3,3)
            endif   

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(9*ij0-8)= AUlu0(9*ij0-8) - RHS_Aij(1,1)
              AUlu0(9*ij0-7)= AUlu0(9*ij0-7) - RHS_Aij(1,2)
              AUlu0(9*ij0-6)= AUlu0(9*ij0-6) - RHS_Aij(1,3)
              AUlu0(9*ij0-5)= AUlu0(9*ij0-5) - RHS_Aij(2,1)
              AUlu0(9*ij0-4)= AUlu0(9*ij0-4) - RHS_Aij(2,2)
              AUlu0(9*ij0-3)= AUlu0(9*ij0-3) - RHS_Aij(2,3)
              AUlu0(9*ij0-2)= AUlu0(9*ij0-2) - RHS_Aij(3,1)
              AUlu0(9*ij0-1)= AUlu0(9*ij0-1) - RHS_Aij(3,2)
              AUlu0(9*ij0  )= AUlu0(9*ij0  ) - RHS_Aij(3,3)
            endif   
          endif
          enddo
        enddo
      enddo

      deallocate (IW1, IW2)
      deallocate (iconFI1L, iconFI1U)
!C===      
      end subroutine FORM_ILU2_33

      end subroutine  hecmw_solve_BLGPBiCG_33
      end module     hecmw_solver_BLGPBiCG_33

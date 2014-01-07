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

!C
!C*** 
!C*** module hecmw_solver_BLCG_33
!C***
!C
      module hecmw_solver_BLCG_33
      contains
!C
!C*** BLCG_33
!C
!C    BLCG_33 solves the linear system Ax = b with 3*3 block matrix 
!C    using the Conjugate Gradient iterative method with the following
!C    FULL-BLOCK TYPE preconditioners :
!C
!C      (1) Block IC(0) with Additive Schwartz Domain Decomposition
!C      (2) Block IC(1) with Additive Schwartz Domain Decomposition
!C      (3) Block IC(2) with Additive Schwartz Domain Decomposition
!C
      subroutine hecmw_solve_BLCG_33                                    &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, RESID, SIGMA, SIGMA_DIAG, ITER, ERROR, my_rank,         &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM , PRECOND, iterPREmax,         &
     &                       Tset, Tsol, Tcomm, ITERlog)

      use  hecmw_util
      use  hecmw_solver_SR_33
      use m_hecmw_solve_error
      use m_hecmw_comm_f

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

      real(kind=kreal), dimension(:),    allocatable :: WS, WR
      real(kind=kreal), dimension(:),    allocatable :: SCALE
      real(kind=kreal), dimension(:,:),  allocatable :: WW
      real(kind=kreal), dimension(:),    allocatable :: ALU

      real(kind=kreal), dimension(:), allocatable :: AUlu0, ALlu0
      real(kind=kreal), dimension(:), allocatable :: Dlu0

      integer(kind=kint), dimension(:), allocatable :: inumFI1L, FI1L
      integer(kind=kint), dimension(:), allocatable :: inumFI1U, FI1U

      integer(kind=kint ) :: P, Q, R, Z, ZP, MAXIT

      ! local variables
      integer(kind=kint )::i,j,k,isU,ieU,isL,ieL,iterPRE,indexA,indexB,indexC
      integer(kind=kint )::ip,ip1,ip2,ip3,inod,iq1,iq2,iq3
      integer(kind=kint ):: ns, nr
      real   (kind=kreal):: TOL, D11,D22,D33, ALUtmp(3,3)
      real   (kind=kreal)::S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2,X3
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,C10,C1,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::COMMtime
      real   (kind=kreal)::WVAL1,WVAL2,WVAL3,SW1,SW2,SW3,WV1,WV2,WV3

      START_TIME= HECMW_WTIME()
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0

      if( NEIBPETOT > 0 ) then
        ns = STACK_EXPORT(NEIBPETOT)
        nr = STACK_IMPORT(NEIBPETOT)
      else
        ns = 0
        nr = 0
      end if

      allocate (WW(3*NP,4))
      allocate (WS(3*ns), WR(3*nr))
      allocate (SCALE(3*NP))
      allocate (ALU(9*N))

      COMMtime= 0.d0

      SCALE= 1.d0

      R = 1
      Z = 2
      Q = 2
      P = 3
      ZP= 4
      
      MAXIT  = ITER
       TOL   = RESID           

!C
!C-- exchanging DIAGONAL components
        WW= 0.d0
!*voption indep (WW,D)
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

!*voption indep (D,WW)
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

!C
!C-- SCALING
        do i= 1, NP
          SCALE (3*i-2)= 1.d0/dsqrt(dabs(D(9*i-8)))
          SCALE (3*i-1)= 1.d0/dsqrt(dabs(D(9*i-4)))
          SCALE (3*i  )= 1.d0/dsqrt(dabs(D(9*i  )))
        enddo

      call hecmw_solve_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, SCALE, SOLVER_COMM,my_rank)

      do i= 1, NP
        ip1= 3*i-2
        ip2= 3*i-1
        ip3= 3*i 
        D(9*i-8)= D(9*i-8)*SCALE(ip1)*SCALE(ip1)
        D(9*i-7)= D(9*i-7)*SCALE(ip1)*SCALE(ip2)
        D(9*i-6)= D(9*i-6)*SCALE(ip1)*SCALE(ip3)
        D(9*i-5)= D(9*i-5)*SCALE(ip2)*SCALE(ip1)
        D(9*i-4)= D(9*i-4)*SCALE(ip2)*SCALE(ip2)
        D(9*i-3)= D(9*i-3)*SCALE(ip2)*SCALE(ip3)
        D(9*i-2)= D(9*i-2)*SCALE(ip3)*SCALE(ip1)
        D(9*i-1)= D(9*i-1)*SCALE(ip3)*SCALE(ip2)
        D(9*i  )= D(9*i  )*SCALE(ip3)*SCALE(ip3)

          isL= INL(i-1) + 1
          ieL= INL(i  ) 
!*voption indep (IAL,AL,SCALE)
          do k= isL, ieL
             inod= IAL(k)
             iq1= 3*inod - 2
             iq2= 3*inod - 1
             iq3= 3*inod 
            AL(9*k-8)= AL(9*k-8)*SCALE(ip1)*SCALE(iq1)
            AL(9*k-7)= AL(9*k-7)*SCALE(ip1)*SCALE(iq2)
            AL(9*k-6)= AL(9*k-6)*SCALE(ip1)*SCALE(iq3)
            AL(9*k-5)= AL(9*k-5)*SCALE(ip2)*SCALE(iq1)
            AL(9*k-4)= AL(9*k-4)*SCALE(ip2)*SCALE(iq2)
            AL(9*k-3)= AL(9*k-3)*SCALE(ip2)*SCALE(iq3)
            AL(9*k-2)= AL(9*k-2)*SCALE(ip3)*SCALE(iq1)
            AL(9*k-1)= AL(9*k-1)*SCALE(ip3)*SCALE(iq2)
            AL(9*k  )= AL(9*k  )*SCALE(ip3)*SCALE(iq3)
          enddo

          isU= INU(i-1) + 1
          ieU= INU(i  ) 
!*voption indep (IAU,AU,SCALE)
          do k= isU, ieU
             inod= IAU(k)
             iq1= 3*inod - 2
             iq2= 3*inod - 1
             iq3= 3*inod 
            AU(9*k-8)= AU(9*k-8)*SCALE(ip1)*SCALE(iq1)
            AU(9*k-7)= AU(9*k-7)*SCALE(ip1)*SCALE(iq2)
            AU(9*k-6)= AU(9*k-6)*SCALE(ip1)*SCALE(iq3)
            AU(9*k-5)= AU(9*k-5)*SCALE(ip2)*SCALE(iq1)
            AU(9*k-4)= AU(9*k-4)*SCALE(ip2)*SCALE(iq2)
            AU(9*k-3)= AU(9*k-3)*SCALE(ip2)*SCALE(iq3)
            AU(9*k-2)= AU(9*k-2)*SCALE(ip3)*SCALE(iq1)
            AU(9*k-1)= AU(9*k-1)*SCALE(ip3)*SCALE(iq2)
            AU(9*k  )= AU(9*k  )*SCALE(ip3)*SCALE(iq3)
          enddo
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
        do ip= 1, N
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
        do ip= 1, N
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

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===
!*voption indep (B,SCALE)
      do i= 1, N
        B(3*i-2)= B(3*i-2) * SCALE(3*i-2)
        B(3*i-1)= B(3*i-1) * SCALE(3*i-1)
        B(3*i  )= B(3*i  ) * SCALE(3*i  )
      enddo
!C
!C-- INTERFACE data EXCHANGE

      call hecmw_solve_SEND_RECV_33                                     &
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
      enddo

      BNRM20= 0.d0
!*voption indep (B)
      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
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
      do iter= 1, MAXIT
!      call MPI_BARRIER(SOLVER_COMM,ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

      indexA= ZP
      indexB= Z
      indexC= R
!C
!C== Block SSOR
      if (PRECOND.eq.10) then

      do i= 1, N
        WW(3*i-2,indexA)= WW(3*i-2,indexC)
        WW(3*i-1,indexA)= WW(3*i-1,indexC)
        WW(3*i  ,indexA)= WW(3*i  ,indexC)
      enddo
      do i= 1+N, NP
        WW(3*i-2,indexA)= 0.d0
        WW(3*i-1,indexA)= 0.d0
        WW(3*i  ,indexA)= 0.d0
      enddo

      do i= 1, NP
        WW(3*i-2,indexB )= 0.d0
        WW(3*i-1,indexB )= 0.d0
        WW(3*i  ,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax

!C
!C-- FORWARD
        do i= 1, N
          SW1= WW(3*i-2,indexA)
          SW2= WW(3*i-1,indexA)
          SW3= WW(3*i  ,indexA)
          isL= INL(i-1)+1
          ieL= INL(i)
         do j= isL, ieL
              k= IAL(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
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
          WW(3*i-2,indexA)= X1
          WW(3*i-1,indexA)= X2
          WW(3*i  ,indexA)= X3
        enddo

!C
!C-- BACKWARD
        do i= N, 1, -1
          isU= INU(i-1) + 1
          ieU= INU(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= isU, ieU
              k= IAU(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
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
          WW(3*i-2,indexA)=  WW(3*i-2,indexA) - X1
          WW(3*i-1,indexA)=  WW(3*i-1,indexA) - X2
          WW(3*i  ,indexA)=  WW(3*i  ,indexA) - X3
        enddo

!C
!C-- additive Schwartz

        do i= 1, N
          WW(3*i-2,indexB)= WW(3*i-2,indexB) + WW(3*i-2,indexA)
          WW(3*i-1,indexB)= WW(3*i-1,indexB) + WW(3*i-1,indexA)
          WW(3*i  ,indexB)= WW(3*i  ,indexB) + WW(3*i  ,indexA)
        enddo

        if (iterPRE.eq.iterPREmax) exit

        S_TIME= HECMW_WTIME()
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexB),            &
     &       SOLVER_COMM,  my_rank)
        E_TIME= HECMW_WTIME()
        COMMtime = COMMtime + E_TIME - S_TIME

        do j= 1, N
           X1= WW(3*j-2,indexB)
           X2= WW(3*j-1,indexB)
           X3= WW(3*j  ,indexB)
          WV1= WW(3*j-2,indexC) - D(9*j-8)*X1-D(9*j-7)*X2-D(9*j-6)*X3
          WV2= WW(3*j-1,indexC) - D(9*j-5)*X1-D(9*j-4)*X2-D(9*j-3)*X3
          WV3= WW(3*j  ,indexC) - D(9*j-2)*X1-D(9*j-1)*X2-D(9*j  )*X3
          do k= INL(j-1)+1, INL(j)
              i= IAL(k)
             X1= WW(3*i-2,indexB)
             X2= WW(3*i-1,indexB)
             X3= WW(3*i  ,indexB)
            WV1= WV1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
            WV2= WV2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
            WV3= WV3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
          enddo
          do k= INU(j-1)+1, INU(j)
              i= IAU(k)
             X1= WW(3*i-2,indexB)
             X2= WW(3*i-1,indexB)
             X3= WW(3*i  ,indexB)
            WV1= WV1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
            WV2= WV2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
            WV3= WV3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
          enddo

          WW(3*j-2,indexA)= WV1
          WW(3*j-1,indexA)= WV2
          WW(3*j  ,indexA)= WV3
        enddo

      enddo
      endif

!C
!C== Block ILU
      if (PRECOND.eq.11.or.PRECOND.eq.12) then
      do i= 1, N
        WW(3*i-2,indexA)= WW(3*i-2,indexC)
        WW(3*i-1,indexA)= WW(3*i-1,indexC)
        WW(3*i  ,indexA)= WW(3*i  ,indexC)
      enddo
      do i= 1+N, NP
        WW(3*i-2,indexA)= 0.d0
        WW(3*i-1,indexA)= 0.d0
        WW(3*i  ,indexA)= 0.d0
      enddo

      do i= 1, NP
        WW(3*i-2,indexB )= 0.d0
        WW(3*i-1,indexB )= 0.d0
        WW(3*i  ,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax

!C
!C-- FORWARD

        do i= 1, N
          SW1= WW(3*i-2,indexA)
          SW2= WW(3*i-1,indexA)
          SW3= WW(3*i  ,indexA)
          isL= inumFI1L(i-1)+1
          ieL= inumFI1L(i)
         do j= isL, ieL
              k= FI1L(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
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
          WW(3*i-2,indexA)= X1
          WW(3*i-1,indexA)= X2
          WW(3*i  ,indexA)= X3
        enddo

!C
!C-- BACKWARD

        do i= N, 1, -1
          isU= inumFI1U(i-1) + 1
          ieU= inumFI1U(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= isU, ieU
              k= FI1U(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
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
          WW(3*i-2,indexA)=  WW(3*i-2,indexA) - X1
          WW(3*i-1,indexA)=  WW(3*i-1,indexA) - X2
          WW(3*i  ,indexA)=  WW(3*i  ,indexA) - X3
        enddo

!C
!C-- additive Schwartz

      S_TIME= HECMW_WTIME()
        call hecmw_solve_SEND_RECV_33                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexA) ,           &
     &       SOLVER_COMM, my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

        do i= 1, NP
          WW(3*i-2,indexB)= WW(3*i-2,indexB) + WW(3*i-2,indexA)         
          WW(3*i-1,indexB)= WW(3*i-1,indexB) + WW(3*i-1,indexA)         
          WW(3*i  ,indexB)= WW(3*i  ,indexB) + WW(3*i  ,indexA)         
        enddo

        if (iterPRE.eq.iterPREmax) exit

        do j= 1, N
           X1= WW(3*j-2,indexB)
           X2= WW(3*j-1,indexB)
           X3= WW(3*j  ,indexB)
          WV1= WW(3*j-2,indexC) - D(9*j-8)*X1-D(9*j-7)*X2-D(9*j-6)*X3
          WV2= WW(3*j-1,indexC) - D(9*j-5)*X1-D(9*j-4)*X2-D(9*j-3)*X3
          WV3= WW(3*j  ,indexC) - D(9*j-2)*X1-D(9*j-1)*X2-D(9*j  )*X3
          do k= INL(j-1)+1, INL(j)
              i= IAL(k)
             X1= WW(3*i-2,indexB)
             X2= WW(3*i-1,indexB)
             X3= WW(3*i  ,indexB)
            WV1= WV1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
            WV2= WV2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
            WV3= WV3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
          enddo
          do k= INU(j-1)+1, INU(j)
              i= IAU(k)
             X1= WW(3*i-2,indexB)
             X2= WW(3*i-1,indexB)
             X3= WW(3*i  ,indexB)
            WV1= WV1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
            WV2= WV2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
            WV3= WV3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
          enddo
          WW(3*j-2,indexA)= WV1
          WW(3*j-1,indexA)= WV2
          WW(3*j  ,indexA)= WV3
        enddo

      enddo
      endif
!C==
!C===
      
!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0
!*voption indep (WW)
      do i= 1, N
        RHO0= RHO0 + WW(3*i-2,R)*WW(3*i-2,Z) + WW(3*i-1,R)*WW(3*i-1,Z)  &
     &             + WW(3*i  ,R)*WW(3*i  ,Z)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1(RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then
!*voption indep (WW)
        do i= 1, N
          WW(3*i-2,P)= WW(3*i-2,Z)
          WW(3*i-1,P)= WW(3*i-1,Z)
          WW(3*i  ,P)= WW(3*i  ,Z)
        enddo
       else
         BETA= RHO / RHO1
!*voption indep (WW)
         do i= 1, N
           WW(3*i-2,P)= WW(3*i-2,Z) + BETA*WW(3*i-2,P)
           WW(3*i-1,P)= WW(3*i-1,Z) + BETA*WW(3*i-1,P)
           WW(3*i  ,P)= WW(3*i  ,Z) + BETA*WW(3*i  ,P)
         enddo
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        

!C
!C-- INTERFACE data EXCHANGE

      S_TIME= HECMW_WTIME()
      call hecmw_solve_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,P) , SOLVER_COMM,     &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

!C
!C-- BEGIN calculation
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

        WW(3*j-2,Q)= WVAL1
        WW(3*j-1,Q)= WVAL2
        WW(3*j  ,Q)= WVAL3

      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
!*voption indep (WW)
      do i= 1, N
        C10= C10 + WW(3*i-2,P)*WW(3*i-2,Q) + WW(3*i-1,P)*WW(3*i-1,Q)    &
     &           + WW(3*i  ,P)*WW(3*i  ,Q)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C10, C1,  HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===
!*voption indep (X,WW)
      do i= 1, N
         X(3*i-2)  = X (3*i-2)   + ALPHA * WW(3*i-2,P)
         X(3*i-1)  = X (3*i-1)   + ALPHA * WW(3*i-1,P)
         X(3*i  )  = X (3*i  )   + ALPHA * WW(3*i  ,P)
        WW(3*i-2,R)= WW(3*i-2,R) - ALPHA * WW(3*i-2,Q)
        WW(3*i-1,R)= WW(3*i-1,R) - ALPHA * WW(3*i-1,Q)
        WW(3*i  ,R)= WW(3*i  ,R) - ALPHA * WW(3*i  ,Q)
      enddo

      DNRM20= 0.d0
!*voption indep (WW)
      do i= 1, N
        DNRM20= DNRM20 + WW(3*i-2,R)**2 + WW(3*i-1,R)**2                &
     &                                  + WW(3*i  ,R)**2
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

        RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
        if (my_rank.eq.0.and.ITERlog.eq.1) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
! 1010   format (1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300

        RHO1 = RHO
      enddo
!C===

!C
!C-- INTERFACE data EXCHANGE
   30 continue
      END_TIME= HECMW_WTIME()

!*voption indep (X,B,SCALE)
      do i= 1, N
        X(3*i-2)= X(3*i-2) * SCALE(3*i-2)
        X(3*i-1)= X(3*i-1) * SCALE(3*i-1)
        X(3*i  )= X(3*i  ) * SCALE(3*i  )
        B(3*i-2)= B(3*i-2) / SCALE(3*i-2)
        B(3*i-1)= B(3*i-1) / SCALE(3*i-1)
        B(3*i  )= B(3*i  ) / SCALE(3*i  )
      enddo

      do i= 1, NP
        ip1= 3*i-2
        ip2= 3*i-1
        ip3= 3*i 
        D(9*i-8)= D(9*i-8)/(SCALE(ip1)*SCALE(ip1))
        D(9*i-7)= D(9*i-7)/(SCALE(ip1)*SCALE(ip2))
        D(9*i-6)= D(9*i-6)/(SCALE(ip1)*SCALE(ip3))
        D(9*i-5)= D(9*i-5)/(SCALE(ip2)*SCALE(ip1))
        D(9*i-4)= D(9*i-4)/(SCALE(ip2)*SCALE(ip2))
        D(9*i-3)= D(9*i-3)/(SCALE(ip2)*SCALE(ip3))
        D(9*i-2)= D(9*i-2)/(SCALE(ip3)*SCALE(ip1))
        D(9*i-1)= D(9*i-1)/(SCALE(ip3)*SCALE(ip2))
        D(9*i  )= D(9*i  )/(SCALE(ip3)*SCALE(ip3))

          isL= INL(i-1) + 1
          ieL= INL(i  ) 
!*voption indep (IAL,AL,SCALE)
          do k= isL, ieL
             inod= IAL(k)
             iq1= 3*inod - 2
             iq2= 3*inod - 1
             iq3= 3*inod 
            AL(9*k-8)= AL(9*k-8)/(SCALE(ip1)*SCALE(iq1))
            AL(9*k-7)= AL(9*k-7)/(SCALE(ip1)*SCALE(iq2))
            AL(9*k-6)= AL(9*k-6)/(SCALE(ip1)*SCALE(iq3))
            AL(9*k-5)= AL(9*k-5)/(SCALE(ip2)*SCALE(iq1))
            AL(9*k-4)= AL(9*k-4)/(SCALE(ip2)*SCALE(iq2))
            AL(9*k-3)= AL(9*k-3)/(SCALE(ip2)*SCALE(iq3))
            AL(9*k-2)= AL(9*k-2)/(SCALE(ip3)*SCALE(iq1))
            AL(9*k-1)= AL(9*k-1)/(SCALE(ip3)*SCALE(iq2))
            AL(9*k  )= AL(9*k  )/(SCALE(ip3)*SCALE(iq3))
          enddo

          isU= INU(i-1) + 1
          ieU= INU(i  ) 
!*voption indep (IAU,AU,SCALE)
          do k= isU, ieU
             inod= IAU(k)
             iq1= 3*inod - 2
             iq2= 3*inod - 1
             iq3= 3*inod 
            AU(9*k-8)= AU(9*k-8)/(SCALE(ip1)*SCALE(iq1))
            AU(9*k-7)= AU(9*k-7)/(SCALE(ip1)*SCALE(iq2))
            AU(9*k-6)= AU(9*k-6)/(SCALE(ip1)*SCALE(iq3))
            AU(9*k-5)= AU(9*k-5)/(SCALE(ip2)*SCALE(iq1))
            AU(9*k-4)= AU(9*k-4)/(SCALE(ip2)*SCALE(iq2))
            AU(9*k-3)= AU(9*k-3)/(SCALE(ip2)*SCALE(iq3))
            AU(9*k-2)= AU(9*k-2)/(SCALE(ip3)*SCALE(iq1))
            AU(9*k-1)= AU(9*k-1)/(SCALE(ip3)*SCALE(iq2))
            AU(9*k  )= AU(9*k  )/(SCALE(ip3)*SCALE(iq3))
          enddo
        enddo

      call hecmw_solve_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW, WR, WS, SCALE, ALU)
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
      integer(kind=kint), dimension(:), allocatable :: IW1, IW2
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
      integer(kind=kint) :: NPLf1,NPLf2,NPUf1,NPUf2,iAS,iconIK,iconKJ
      integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
      integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3

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

      end subroutine  hecmw_solve_BLCG_33
      end module     hecmw_solver_BLCG_33

!C
!C***
!C*** fill_in_S33_SORT
!C*** 
!C
      subroutine fill_in_S33_SORT (STEM, INUM, N, NP)
      use hecmw_util
      implicit none
      integer(kind=kint) :: N, NP
      integer(kind=kint), dimension(NP) :: STEM
      integer(kind=kint), dimension(NP) :: INUM
      integer(kind=kint), dimension(:), allocatable :: ISTACK
      integer(kind=kint) :: M,NSTACK,jstack,l,ir,ip,i,j,k,ss,ii,temp,it

      allocate (ISTACK(-NP:+NP)) 

      M     = 100
      NSTACK= NP

      jstack= 0
      l     = 1
      ir    = N

      ip= 0
 1    continue
      ip= ip + 1

      if (ir-l.lt.M) then
        do 12 j= l+1, ir
          ss= STEM(j)
          ii= INUM(j)

          do 11 i= j-1,1,-1
            if (STEM(i).le.ss) goto 2
            STEM(i+1)= STEM(i)
            INUM(i+1)= INUM(i)
 11       continue
          i= 0

 2        continue
            STEM(i+1)= ss
            INUM(i+1)= ii
 12     continue

        if (jstack.eq.0) then
          deallocate (ISTACK)
          return
        endif

        ir = ISTACK(jstack)
         l = ISTACK(jstack-1)
        jstack= jstack - 2
       else

        k= (l+ir) / 2
            temp = STEM(k)
        STEM(k)  = STEM(l+1)
        STEM(l+1)= temp

              it = INUM(k)
        INUM(k)  = INUM(l+1)
        INUM(l+1)= it

        if (STEM(l+1).gt.STEM(ir)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(ir)
          STEM(ir )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(ir)
          INUM(ir )= it
        endif

        if (STEM(l).gt.STEM(ir)) then
             temp = STEM(l)
          STEM(l )= STEM(ir)
          STEM(ir)= temp
               it = INUM(l)
          INUM(l )= INUM(ir)
          INUM(ir)= it
        endif

        if (STEM(l+1).gt.STEM(l)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(l)
          STEM(l  )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(l)
          INUM(l  )= it
        endif

        i= l + 1
        j= ir

        ss= STEM(l)
        ii= INUM(l)

 3      continue
          i= i + 1
          if (STEM(i).lt.ss) goto 3

 4      continue
          j= j - 1
          if (STEM(j).gt.ss) goto 4

        if (j.lt.i)        goto 5

        temp   = STEM(i)
        STEM(i)= STEM(j)
        STEM(j)= temp

        it     = INUM(i)
        INUM(i)= INUM(j)
        INUM(j)= it

        goto 3

 5      continue

        STEM(l)= STEM(j)
        STEM(j)= ss
        INUM(l)= INUM(j)
        INUM(j)= ii

        jstack= jstack + 2

        if (jstack.gt.NSTACK) then
          write (*,*) 'NSTACK overflow'
          stop
        endif

        if (ir-i+1.ge.j-1) then
          ISTACK(jstack  )= ir
          ISTACK(jstack-1)= i
          ir= j-1
         else
          ISTACK(jstack  )= j-1
          ISTACK(jstack-1)= l
          l= i
        endif

      endif

      goto 1

      end

!C
!C***
!C*** ILU1a33
!C***
!C
!C    computes LU factorization of 3*3 Diagonal Block
!C
      subroutine ILU1a33 (ALU, D11,D12,D13,D21,D22,D23,D31,D32,D33)
      use hecmw_util
      implicit none
      real(kind=kreal) :: ALU(3,3), PW(3)
      real(kind=kreal) :: D11,D12,D13,D21,D22,D23,D31,D32,D33
      integer(kind=kint) :: i,j,k,L
 
      ALU(1,1)= D11
      ALU(1,2)= D12
      ALU(1,3)= D13
      ALU(2,1)= D21
      ALU(2,2)= D22
      ALU(2,3)= D23
      ALU(3,1)= D31
      ALU(3,2)= D32
      ALU(3,3)= D33

      do k= 1, 3
        ALU(k,k)= 1.d0/ALU(k,k)
        do i= k+1, 3
          ALU(i,k)= ALU(i,k) * ALU(k,k)
          do j= k+1, 3
            PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
          enddo
          do j= k+1, 3
            ALU(i,j)= PW(j)
          enddo
        enddo
      enddo

      return
      end
      
!C
!C***
!C*** ILU1b33
!C***
!C
!C    computes L_ik * D_k_INV * U_kj at ILU factorization
!C    for 3*3 Block Type Matrix
!C
      subroutine ILU1b33 (RHS_Aij, DkINV, Aik, Akj)
      use hecmw_util
      implicit none
      real(kind=kreal) :: RHS_Aij(3,3), DkINV(3,3), Aik(3,3), Akj(3,3)
      real(kind=kreal) :: X1,X2,X3

!C
!C-- 1st Col.
      X1= Akj(1,1)
      X2= Akj(2,1)
      X3= Akj(3,1)

        X2= X2 - DkINV(2,1)*X1
        X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

        X3= DkINV(3,3)*  X3
        X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
        X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

        RHS_Aij(1,1)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
        RHS_Aij(2,1)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
        RHS_Aij(3,1)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

!C
!C-- 2nd Col.
      X1= Akj(1,2)
      X2= Akj(2,2)
      X3= Akj(3,2)

        X2= X2 - DkINV(2,1)*X1
        X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

        X3= DkINV(3,3)*  X3
        X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
        X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

        RHS_Aij(1,2)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
        RHS_Aij(2,2)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
        RHS_Aij(3,2)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

!C
!C-- 3rd Col.
      X1= Akj(1,3)
      X2= Akj(2,3)
      X3= Akj(3,3)

        X2= X2 - DkINV(2,1)*X1
        X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

        X3= DkINV(3,3)*  X3
        X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
        X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

        RHS_Aij(1,3)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
        RHS_Aij(2,3)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
        RHS_Aij(3,3)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

      return
      end

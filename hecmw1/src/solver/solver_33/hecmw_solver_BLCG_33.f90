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
      use  hecmw_precond_BILU_33
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

      integer(kind=kint ) :: P, Q, R, Z, ZP, MAXIT

      ! local variables
      integer(kind=kint )::i,j,k,isU,ieU,isL,ieL,iterPRE,indexA,indexB,indexC
      integer(kind=kint )::ip1,ip2,ip3,inod,iq1,iq2,iq3
      integer(kind=kint ):: ns, nr
      real   (kind=kreal):: TOL
      real   (kind=kreal)::S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2,X3
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,C10,C1,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::COMMtime
      real   (kind=kreal)::WVAL1,WVAL2,WVAL3,WV1,WV2,WV3

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
      call hecmw_precond_BILU_33_setup                    &
           &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
           &    SIGMA, SIGMA_DIAG, my_rank,                     &
           &    PRECOND)

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
!C== Block ILU
!C==
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

        call hecmw_precond_BILU_33_apply(N, WW(:,indexA))

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

      deallocate (WW, WR, WS, SCALE)
      call hecmw_precond_BILU_33_clear(PRECOND)

      Tset = Tset + SETupTIME
      Tsol = Tsol + END_TIME - START_TIME
      Tcomm= COMMtime

      end subroutine  hecmw_solve_BLCG_33
      end module     hecmw_solver_BLCG_33

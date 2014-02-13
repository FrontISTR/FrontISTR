!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
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
!C*** module hecmw_solver_BLCG_22
!C***
!C
      module hecmw_solver_BLCG_22
      contains
!C
!C*** BLCG_22
!C
!C    BLCG_22 solves the linear system Ax = b with 2*2 block matrix
!C    using the Conjugate Gradient iterative method with the following
!C    FULL-BLOCK TYPE preconditioners :
!C
!C      (1) Block IC(0) with Additive Shcwartz Domain Decomposition
!C      (2) Block IC(1) with Additive Shcwartz Domain Decomposition
!C      (3) Block IC(2) with Additive Shcwartz Domain Decomposition
!C
!C
      subroutine hecmw_solve_BLCG_22                                    &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, RESID, SIGMA, SIGMA_DIAG, ITER, ERROR, my_rank,         &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM , PRECOND, iterPREmax,         &
     &                       Tset, Tsol, Tcomm, ITERlog)

      use  hecmw_util
      use  hecmw_solver_SR_22
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

      real(kind=kreal), dimension(2*NP)   , intent(inout):: B, X
      real(kind=kreal), dimension(4*NPL), intent(inout):: AL
      real(kind=kreal), dimension(4*NPU), intent(inout):: AU
      real(kind=kreal), dimension(4*NP ), intent(inout):: D

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)

      real(kind=kreal), dimension(:),    allocatable       :: WS, WR
      real(kind=kreal), dimension(:),    allocatable, save :: SCALE
      real(kind=kreal), dimension(:,:),  allocatable       :: WW
      real(kind=kreal), dimension(:),    allocatable, save :: ALU

      real(kind=kreal), dimension(:), allocatable :: AUlu0, ALlu0
      real(kind=kreal), dimension(:), allocatable :: Dlu0

      integer(kind=kint), dimension(:), allocatable :: inumFI1L, FI1L
      integer(kind=kint), dimension(:), allocatable :: inumFI1U, FI1U

      integer(kind=kint ) :: P, Q, R, Z, ZP, MAXIT
      real   (kind=kreal) :: TOL, COMMtime

      ! local variables
      integer(kind=kint )::i,j,jj,jj1,ij0,k,kk,ik,kk1,kk2,iSj,iEj,iSk,iEk,isU,ieU,isL,ieL
      integer(kind=kint )::iterPRE,ip,ip1,ip2,iq1,iq2,iAS
      integer(kind=kint )::indexA,indexB,indexC,inod,NPLf1,NPLf2,NPUf1,NPUf2
      integer(kind=kint )::icou,icou0,icouL,L,icouL1,icouL2,icouL3
      integer(kind=kint )::icouU,icouU1,icouU2,icouU3,iconIK,iconKJ
      integer(kind=kint )::ns,nr
      real   (kind=kreal)::S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,C10,C1,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::WVAL1,WVAL2,SW1,SW2,WV1,WV2
      real   (kind=kreal)::D11,D22,D12,D21,ALUtmp(2,2)

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

      allocate (WW(2*NP,4))
      allocate (WS(2*ns), WR(2*nr))
      allocate (SCALE(2*NP))
      allocate (ALU(4*N))

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
!C-- SCALING
        do i= 1, NP
          SCALE (2*i-1)= 1.d0/dsqrt(dabs(D(4*i-3)))
          SCALE (2*i  )= 1.d0/dsqrt(dabs(D(4*i  )))
        enddo

      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, SCALE, SOLVER_COMM,my_rank)

      do i= 1, NP
        ip1= 2*i-1
        ip2= 2*i
        D(4*i-3)= D(4*i-3)*SCALE(ip1)*SCALE(ip1)
        D(4*i-2)= D(4*i-2)*SCALE(ip1)*SCALE(ip2)
        D(4*i-1)= D(4*i-1)*SCALE(ip2)*SCALE(ip1)
        D(4*i  )= D(4*i  )*SCALE(ip2)*SCALE(ip2)

          isL= INL(i-1) + 1
          ieL= INL(i  )
!*voption indep (IAL,AL,SCALE)
          do k= isL, ieL
             inod= IAL(k)
             iq1= 2*inod - 1
             iq2= 2*inod
            AL(4*k-3)= AL(4*k-3)*SCALE(ip1)*SCALE(iq1)
            AL(4*k-2)= AL(4*k-2)*SCALE(ip1)*SCALE(iq2)
            AL(4*k-1)= AL(4*k-1)*SCALE(ip2)*SCALE(iq1)
            AL(4*k  )= AL(4*k  )*SCALE(ip2)*SCALE(iq2)
          enddo

          isU= INU(i-1) + 1
          ieU= INU(i  )
!*voption indep (IAU,AU,SCALE)
          do k= isU, ieU
             inod= IAU(k)
             iq1= 2*inod - 1
             iq2= 2*inod
            AU(4*k-3)= AU(4*k-3)*SCALE(ip1)*SCALE(iq1)
            AU(4*k-2)= AU(4*k-2)*SCALE(ip1)*SCALE(iq2)
            AU(4*k-1)= AU(4*k-1)*SCALE(ip2)*SCALE(iq1)
            AU(4*k  )= AU(4*k  )*SCALE(ip2)*SCALE(iq2)
          enddo
        enddo
!C===

!C
!C +-------------------+
!C | ILU decomposition |
!C +-------------------+
!C===
      if (PRECOND.eq.11) call FORM_ILU1_22
      if (PRECOND.eq.12) call FORM_ILU2_22

      if (PRECOND.eq.10) then
        do ip= 1, N
          D11= D(4*ip-3) * SIGMA_DIAG
          D22= D(4*ip  ) * SIGMA_DIAG
          call ILU1a22 (ALUtmp,                                    &
     &                  D11, D(4*ip-2), D(4*ip-1), D22)
          ALU(4*ip-3)= ALUtmp(1,1)
          ALU(4*ip-2)= ALUtmp(1,2)
          ALU(4*ip-1)= ALUtmp(2,1)
          ALU(4*ip  )= ALUtmp(2,2)
        enddo
      endif

      if (PRECOND.eq.11.or.PRECOND.eq.12) then
        do ip= 1, N
          call ILU1a22 (ALUtmp,                                    &
     &                  Dlu0(4*ip-3), Dlu0(4*ip-2),                     &
     &                  Dlu0(4*ip-1), Dlu0(4*ip))
          ALU(4*ip-3)= ALUtmp(1,1)
          ALU(4*ip-2)= ALUtmp(1,2)
          ALU(4*ip-1)= ALUtmp(2,1)
          ALU(4*ip  )= ALUtmp(2,2)

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
        B(2*i-1)= B(2*i-1) * SCALE(2*i-1)
        B(2*i  )= B(2*i  ) * SCALE(2*i  )
      enddo

!C
!C-- INTERFACE data EXCHANGE

      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      do j= 1, N
        X1= X(2*j-1)
        X2= X(2*j  )
        WVAL1= B(2*j-1) - D(4*j-3)*X1 - D(4*j-2)*X2
        WVAL2= B(2*j  ) - D(4*j-1)*X1 - D(4*j  )*X2
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= X(2*i-1)
          X2= X(2*i  )
          WVAL1= WVAL1 - AL(4*k-3)*X1 - AL(4*k-2)*X2
          WVAL2= WVAL2 - AL(4*k-1)*X1 - AL(4*k  )*X2
        enddo

        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= X(2*i-1)
          X2= X(2*i  )
          WVAL1= WVAL1 - AU(4*k-3)*X1 - AU(4*k-2)*X2
          WVAL2= WVAL2 - AU(4*k-1)*X1 - AU(4*k  )*X2
        enddo
        WW(2*j-1,R)= WVAL1
        WW(2*j  ,R)= WVAL2
      enddo

      BNRM20= 0.d0
      do i= 1, N
        BNRM20= BNRM20 + B(2*i-1)**2 + B(2*i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM,SOLVER_COMM)
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0
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
        WW(2*i-1,indexA)= WW(2*i-1,indexC)
        WW(2*i  ,indexA)= WW(2*i  ,indexC)
      enddo

      do i= 1, NP
        WW(2*i-1,indexB )= 0.d0
        WW(2*i  ,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax

!C
!C-- FORWARD
        do i= 1, N
          SW1= WW(2*i-1,indexA)
          SW2= WW(2*i  ,indexA)
          isL= INL(i-1)+1
          ieL= INL(i)
         do j= isL, ieL
              k= IAL(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 - AL(4*j-3)*X1 - AL(4*j-2)*X2
            SW2= SW2 - AL(4*j-1)*X1 - AL(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= X1
          WW(2*i  ,indexA)= X2
        enddo

!C
!C-- BACKWARD

        do i= N, 1, -1
          isU= INU(i-1) + 1
          ieU= INU(i)
          SW1= 0.d0
          SW2= 0.d0
          do j= isU, ieU
              k= IAU(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 + AU(4*j-3)*X1 + AU(4*j-2)*X2
            SW2= SW2 + AU(4*j-1)*X1 + AU(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= WW(2*i-1,indexA) - X1
          WW(2*i  ,indexA)= WW(2*i  ,indexA) - X2
        enddo

!C
!C-- additive Schwartz

        do i= 1, N
          WW(2*i-1,indexB)= WW(2*i-1,indexB) + WW(2*i-1,indexA)
          WW(2*i  ,indexB)= WW(2*i  ,indexB) + WW(2*i  ,indexA)
        enddo


        if (iterPRE.eq.iterPREmax) goto 750

      S_TIME= HECMW_WTIME()
        call hecmw_solve_SEND_RECV_22                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexB) ,           &
     &       SOLVER_COMM, my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

        do j= 1, N
           X1= WW(2*j-1,indexB)
           X2= WW(2*j  ,indexB)
          WV1= WW(2*j-1,indexC) - D(4*j-3)*X1-D(4*j-2)*X2
          WV2= WW(2*j  ,indexC) - D(4*j-1)*X1-D(4*j  )*X2
          do k= INL(j-1)+1, INL(j)
              i= IAL(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AL(4*k-3)*X1 - AL(4*k-2)*X2
            WV2= WV2 - AL(4*k-1)*X1 - AL(4*k  )*X2
          enddo
          do k= INU(j-1)+1, INU(j)
              i= IAU(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AU(4*k-3)*X1 - AU(4*k-2)*X2
            WV2= WV2 - AU(4*k-1)*X1 - AU(4*k  )*X2
          enddo
          WW(2*j-1,indexA)= WV1
          WW(2*j  ,indexA)= WV2
        enddo

 750    continue

      enddo
      endif

!C
!C== Block ILU
      if (PRECOND.eq.11.or.PRECOND.eq.12) then
      do i= 1, N
        WW(2*i-1,indexA)= WW(2*i-1,indexC)
        WW(2*i  ,indexA)= WW(2*i  ,indexC)
      enddo

      do i= 1, NP
        WW(2*i-1,indexB )= 0.d0
        WW(2*i  ,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax
        do i= 1+N, NP
          WW(2*i-1,indexA)= 0.d0
          WW(2*i  ,indexA)= 0.d0
        enddo

!C
!C-- FORWARD
        do i= 1, N
          SW1= WW(2*i-1,indexA)
          SW2= WW(2*i  ,indexA)
          isL= inumFI1L(i-1)+1
          ieL= inumFI1L(i)
         do j= isL, ieL
              k= FI1L(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 - ALlu0(4*j-3)*X1 - ALlu0(4*j-2)*X2
            SW2= SW2 - ALlu0(4*j-1)*X1 - ALlu0(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= X1
          WW(2*i  ,indexA)= X2
        enddo

!C
!C-- BACKWARD

        do i= N, 1, -1
          isU= inumFI1U(i-1) + 1
          ieU= inumFI1U(i)
          SW1= 0.d0
          SW2= 0.d0
          do j= isU, ieU
              k= FI1U(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 + AUlu0(4*j-3)*X1 + AUlu0(4*j-2)*X2
            SW2= SW2 + AUlu0(4*j-1)*X1 + AUlu0(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= WW(2*i-1,indexA) - X1
          WW(2*i  ,indexA)= WW(2*i  ,indexA) - X2
        enddo

!C
!C-- additive Schwartz

      S_TIME= HECMW_WTIME()
        call hecmw_solve_SEND_RECV_22                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexA) ,           &
     &       SOLVER_COMM, my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

        do i= 1, NP
          WW(2*i-1,indexB)= WW(2*i-1,indexB) + WW(2*i-1,indexA)
          WW(2*i  ,indexB)= WW(2*i  ,indexB) + WW(2*i  ,indexA)
        enddo


        if (iterPRE.eq.iterPREmax) goto 850

        do j= 1, N
           X1= WW(2*j-1,indexB)
           X2= WW(2*j  ,indexB)
          WV1= WW(2*j-1,indexC) - D(4*j-3)*X1-D(4*j-2)*X2
          WV2= WW(2*j  ,indexC) - D(4*j-1)*X1-D(4*j  )*X2
          do k= INL(j-1)+1, INL(j)
              i= IAL(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AL(4*k-3)*X1 - AL(4*k-2)*X2
            WV2= WV2 - AL(4*k-1)*X1 - AL(4*k  )*X2
          enddo
          do k= INU(j-1)+1, INU(j)
              i= IAU(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AU(4*k-3)*X1 - AU(4*k-2)*X2
            WV2= WV2 - AU(4*k-1)*X1 - AU(4*k  )*X2
          enddo
          WW(2*j-1,indexA)= WV1
          WW(2*j  ,indexA)= WV2
        enddo

 850    continue

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
        RHO0= RHO0 + WW(2*i-1,R)*WW(2*i-1,Z) + WW(2*i  ,R)*WW(2*i  ,Z)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
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
          WW(2*i-1,P)= WW(2*i-1,Z)
          WW(2*i  ,P)= WW(2*i  ,Z)
        enddo
       else
         BETA= RHO / RHO1
!*voption indep (WW)
         do i= 1, N
           WW(2*i-1,P)= WW(2*i-1,Z) + BETA*WW(2*i-1,P)
           WW(2*i  ,P)= WW(2*i  ,Z) + BETA*WW(2*i  ,P)
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
      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,P) , SOLVER_COMM,     &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      indexA= P
      indexB= Q
!C
!C-- BEGIN calculation

      do j= 1, N
           X1= WW(2*j-1,indexA)
           X2= WW(2*j  ,indexA)
        WVAL1= D(4*j-3)*X1 + D(4*j-2)*X2
        WVAL2= D(4*j-1)*X1 + D(4*j  )*X2
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(2*i-1,indexA)
          X2= WW(2*i  ,indexA)
          WVAL1= WVAL1 + AL(4*k-3)*X1 + AL(4*k-2)*X2
          WVAL2= WVAL2 + AL(4*k-1)*X1 + AL(4*k  )*X2
        enddo
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(2*i-1,indexA)
          X2= WW(2*i  ,indexA)
          WVAL1= WVAL1 + AU(4*k-3)*X1 + AU(4*k-2)*X2
          WVAL2= WVAL2 + AU(4*k-1)*X1 + AU(4*k  )*X2
        enddo

        WW(2*j-1,indexB)= WVAL1
        WW(2*j  ,indexB)= WVAL2
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
        C10= C10 + WW(2*i-1,P)*WW(2*i-1,Q) + WW(2*i,P)*WW(2*i,Q)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C10, C1, HECMW_SUM, SOLVER_COMM )
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
         X(2*i-1)  = X (2*i-1)   + ALPHA * WW(2*i-1,P)
         X(2*i  )  = X (2*i  )   + ALPHA * WW(2*i  ,P)
        WW(2*i-1,R)= WW(2*i-1,R) - ALPHA * WW(2*i-1,Q)
        WW(2*i  ,R)= WW(2*i  ,R) - ALPHA * WW(2*i  ,Q)
      enddo

      DNRM20= 0.d0
!*voption indep (WW)
      do i= 1, N
        DNRM20= DNRM20 + WW(2*i-1,R)**2 + WW(2*i,R)**2
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2,HECMW_SUM, SOLVER_COMM )
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
        X(2*i-1)= X(2*i-1) * SCALE(2*i-1)
        X(2*i  )= X(2*i  ) * SCALE(2*i  )
        B(2*i-1)= B(2*i-1) / SCALE(2*i-1)
        B(2*i  )= B(2*i  ) / SCALE(2*i  )
      enddo

      do i= 1, NP
        ip1= 2*i-1
        ip2= 2*i
        D(4*i-3)= D(4*i-3)/(SCALE(ip1)*SCALE(ip1))
        D(4*i-2)= D(4*i-2)/(SCALE(ip1)*SCALE(ip2))
        D(4*i-1)= D(4*i-1)/(SCALE(ip2)*SCALE(ip1))
        D(4*i  )= D(4*i  )/(SCALE(ip2)*SCALE(ip2))

          isL= INL(i-1) + 1
          ieL= INL(i  )
!*voption indep (IAL,AL,SCALE)
          do k= isL, ieL
             inod= IAL(k)
             iq1= 2*inod - 1
             iq2= 2*inod
            AL(4*k-3)= AL(4*k-3)/(SCALE(ip1)*SCALE(iq1))
            AL(4*k-2)= AL(4*k-2)/(SCALE(ip1)*SCALE(iq2))
            AL(4*k-1)= AL(4*k-1)/(SCALE(ip2)*SCALE(iq1))
            AL(4*k  )= AL(4*k  )/(SCALE(ip2)*SCALE(iq2))
          enddo

          isU= INU(i-1) + 1
          ieU= INU(i  )
!*voption indep (IAU,AU,SCALE)
          do k= isU, ieU
             inod= IAU(k)
             iq1= 2*inod - 1
             iq2= 2*inod
            AU(4*k-3)= AU(4*k-3)/(SCALE(ip1)*SCALE(iq1))
            AU(4*k-2)= AU(4*k-2)/(SCALE(ip1)*SCALE(iq2))
            AU(4*k-1)= AU(4*k-1)/(SCALE(ip2)*SCALE(iq1))
            AU(4*k  )= AU(4*k  )/(SCALE(ip2)*SCALE(iq2))
          enddo
        enddo

      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW, WR, WS, SCALE)
      if (PRECOND.ne.10) then
        deallocate (AUlu0, ALlu0, Dlu0)
        deallocate (inumFI1L, inumFI1U, FI1L, FI1U, ALU)
      endif

      Tset = Tset + SETupTIME
      Tsol = Tsol + END_TIME - START_TIME
      Tcomm= COMMtime
      contains

!C
!C***
!C*** FORM_ILU1_22
!C***
!C
!C    form ILU(1) matrix
!C
      subroutine FORM_ILU1_22
      integer(kind=kint), dimension(:), allocatable :: IW1 , IW2
      integer(kind=kint), dimension(:), allocatable :: IWsL, IWsU
      real (kind=kreal),  dimension(2,2) :: RHS_Aij, DkINV, Aik, Akj

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
      allocate (ALlu0(4*(NPL+NPLf1)), AUlu0(4*(NPU+NPUf1)))

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
        call fill_in_S22_SORT (IW1, IW2, icouL3, NP)

        do k= 1, icouL3
          FI1L (k+isL)= IW1(k)
          ik= IW2(k)
          if (ik.le.INL(i)-INL(i-1)) then
            kk1= 4*( k+isL)
            kk2= 4*(ik+INL(i-1))
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
        call fill_in_S22_SORT (IW1, IW2, icouU3, NP)

        do k= 1, icouU3
          FI1U (k+isU)= IW1(k)
          ik= IW2(k)
          if (ik.le.INU(i)-INU(i-1)) then
            kk1= 4*( k+isU)
            kk2= 4*(ik+INU(i-1))
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
      allocate (Dlu0(4*NP))
      Dlu0= D

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
          D11= Dlu0(4*k-3)
          D12= Dlu0(4*k-2)
          D21= Dlu0(4*k-1)
          D22= Dlu0(4*k  )

          call ILU1a22 (DkINV, D11,D12, D21,D22)

          do kk1= inumFI1L(i-1)+1, inumFI1L(i)
            if (k.eq.FI1L(kk1)) then
              Aik(1,1)= ALlu0(4*kk1-3)
              Aik(1,2)= ALlu0(4*kk1-2)
              Aik(2,1)= ALlu0(4*kk1-1)
              Aik(2,2)= ALlu0(4*kk1  )
              exit
            endif
          enddo

          do jj= INU(k-1)+1, INU(k)
            j= IAU(jj)
            do jj1= inumFI1U(k-1)+1, inumFI1U(k)
              if (j.eq.FI1U(jj1)) then
                Akj(1,1)= AUlu0(4*jj1-3)
                Akj(1,2)= AUlu0(4*jj1-2)
                Akj(2,1)= AUlu0(4*jj1-1)
                Akj(2,2)= AUlu0(4*jj1  )
                exit
              endif
            enddo

            call ILU1b22 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(4*i-3)= Dlu0(4*i-3) - RHS_Aij(1,1)
              Dlu0(4*i-2)= Dlu0(4*i-2) - RHS_Aij(1,2)
              Dlu0(4*i-1)= Dlu0(4*i-1) - RHS_Aij(2,1)
              Dlu0(4*i  )= Dlu0(4*i  ) - RHS_Aij(2,2)
            endif

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(4*ij0-3)= ALlu0(4*ij0-3) - RHS_Aij(1,1)
              ALlu0(4*ij0-2)= ALlu0(4*ij0-2) - RHS_Aij(1,2)
              ALlu0(4*ij0-1)= ALlu0(4*ij0-1) - RHS_Aij(2,1)
              ALlu0(4*ij0  )= ALlu0(4*ij0  ) - RHS_Aij(2,2)
            endif

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(4*ij0-3)= AUlu0(4*ij0-3) - RHS_Aij(1,1)
              AUlu0(4*ij0-2)= AUlu0(4*ij0-2) - RHS_Aij(1,2)
              AUlu0(4*ij0-1)= AUlu0(4*ij0-1) - RHS_Aij(2,1)
              AUlu0(4*ij0  )= AUlu0(4*ij0  ) - RHS_Aij(2,2)
            endif

          enddo
        enddo
      enddo

      deallocate (IW1, IW2)
!C===
      end subroutine FORM_ILU1_22

!C
!C***
!C*** FORM_ILU2_22
!C***
!C
!C    form ILU(2) matrix
!C
      subroutine FORM_ILU2_22

      integer(kind=kint), dimension(:), allocatable:: IW1 , IW2
      integer(kind=kint), dimension(:), allocatable:: IWsL, IWsU
      integer(kind=kint), dimension(:), allocatable:: iconFI1L, iconFI1U
      integer(kind=kint), dimension(:), allocatable:: inumFI2L, inumFI2U
      integer(kind=kint), dimension(:), allocatable::     FI2L,     FI2U
      real (kind=kreal), dimension(2,2) :: RHS_Aij, DkINV, Aik, Akj
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
      allocate (ALlu0(4*(NPL+NPLf1+NPLf2)))
      allocate (AUlu0(4*(NPU+NPUf1+NPUf2)))

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

        call fill_in_S22_SORT (IW1, IW2, icouL3, NP)

        do k= 1, icouL3
          FI1L (k+isL)= IW1(k)
          ik= IW2(k)
          if (ik.le.INL(i)-INL(i-1)) then
            kk1= 4*( k+isL)
            kk2= 4*(ik+INL(i-1))
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
        call fill_in_S22_SORT (IW1, IW2, icouU3, NP)

        do k= 1, icouU3
          FI1U (k+isU)= IW1(k)
          ik= IW2(k)
          if (ik.le.INU(i)-INU(i-1)) then
            kk1= 4*( k+isU)
            kk2= 4*(ik+INU(i-1))
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
      allocate (Dlu0(4*NP))
      Dlu0= D

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

          D11= Dlu0(4*k-3)
          D12= Dlu0(4*k-2)
          D21= Dlu0(4*k-1)
          D22= Dlu0(4*k  )

          call ILU1a22 (DkINV, D11,D12,D21,D22)

          Aik(1,1)= ALlu0(4*kk-3)
          Aik(1,2)= ALlu0(4*kk-2)
          Aik(2,1)= ALlu0(4*kk-1)
          Aik(2,2)= ALlu0(4*kk  )

          do jj= inumFI1U(k-1)+1, inumFI1U(k)
            j= FI1U(jj)
            iconKJ= iconFI1U(jj)

            if ((iconIK+iconKJ).lt.2) then
            Akj(1,1)= AUlu0(4*jj-3)
            Akj(1,2)= AUlu0(4*jj-2)
            Akj(2,1)= AUlu0(4*jj-1)
            Akj(2,2)= AUlu0(4*jj  )

            call ILU1b22 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(4*i-3)= Dlu0(4*i-3) - RHS_Aij(1,1)
              Dlu0(4*i-2)= Dlu0(4*i-2) - RHS_Aij(1,2)
              Dlu0(4*i-1)= Dlu0(4*i-1) - RHS_Aij(2,1)
              Dlu0(4*i  )= Dlu0(4*i  ) - RHS_Aij(2,2)
            endif

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(4*ij0-3)= ALlu0(4*ij0-3) - RHS_Aij(1,1)
              ALlu0(4*ij0-2)= ALlu0(4*ij0-2) - RHS_Aij(1,2)
              ALlu0(4*ij0-1)= ALlu0(4*ij0-1) - RHS_Aij(2,1)
              ALlu0(4*ij0  )= ALlu0(4*ij0  ) - RHS_Aij(2,2)
            endif

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(4*ij0-3)= AUlu0(4*ij0-3) - RHS_Aij(1,1)
              AUlu0(4*ij0-2)= AUlu0(4*ij0-2) - RHS_Aij(1,2)
              AUlu0(4*ij0-1)= AUlu0(4*ij0-1) - RHS_Aij(2,1)
              AUlu0(4*ij0  )= AUlu0(4*ij0  ) - RHS_Aij(2,2)
            endif
          endif
         enddo
        enddo
      enddo

      deallocate (IW1, IW2)
      deallocate (iconFI1L, iconFI1U)
!C===
      end subroutine FORM_ILU2_22

      end subroutine  hecmw_solve_BLCG_22
      end module     hecmw_solver_BLCG_22

!C
!C***
!C*** fill_in_S22_SORT
!C***
!C
      subroutine fill_in_S22_SORT (STEM, INUM, N, NP)
      use hecmw_util
      implicit none
      integer(kind=kint) :: N, NP
      integer(kind=kint), dimension(NP) :: STEM
      integer(kind=kint), dimension(NP) :: INUM
      integer(kind=kint), dimension(:), allocatable :: ISTACK
      ! local
      integer(kind=kint)::M,NSTACK,jstack,l,ir,ip,j,ss,ii,i,k,temp,it

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

      end subroutine fill_in_S22_SORT

!C
!C***
!C*** ILU1a22
!C***
!C
!C    computes LU factorization of 2*2 Diagonal Block
!C
      subroutine ILU1a22 (ALU, D11,D12,D21,D22)
      use hecmw_util
      implicit none
      real(kind=kreal)::ALU(2,2), PW(2)
      real(kind=kreal)::D11,D12,D21,D22
      integer(kind=kint)::k,L,i,j

      ALU(1,1)= D11
      ALU(1,2)= D12
      ALU(2,1)= D21
      ALU(2,2)= D22

      do k= 1, 2
        ALU(k,k)= 1.d0/ALU(k,k)
        do i= k+1, 2
          ALU(i,k)= ALU(i,k) * ALU(k,k)
          do j= k+1, 2
            PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
          enddo
          do j= k+1, 2
            ALU(i,j)= PW(j)
          enddo
        enddo
      enddo

      return
      end subroutine ILU1a22

!C
!C***
!C*** ILU1b22
!C***
!C
!C    computes L_ik * D_k_INV * U_kj at ILU factorization
!C    for 2*2 Block Type Matrix
!C
      subroutine ILU1b22 (RHS_Aij, DkINV, Aik, Akj)
      use hecmw_util
      implicit none
      real(kind=kreal)::RHS_Aij(2,2), DkINV(2,2), Aik(2,2), Akj(2,2)
      real(kind=kreal)::X1,X2
!C
!C-- 1st Col.
      X1= Akj(1,1)
      X2= Akj(2,1)

        X2= X2 - DkINV(2,1)*X1

        X2= DkINV(2,2)*  X2
        X1= DkINV(1,1)*( X1 - DkINV(1,2)*X2)

        RHS_Aij(1,1)=  Aik(1,1)*X1 + Aik(1,2)*X2
        RHS_Aij(2,1)=  Aik(2,1)*X1 + Aik(2,2)*X2

!C
!C-- 2nd Col.
      X1= Akj(1,2)
      X2= Akj(2,2)

        X2= X2 - DkINV(2,1)*X1
        X2= DkINV(2,2)*  X2
        X1= DkINV(1,1)*( X1 - DkINV(1,2)*X2)

        RHS_Aij(1,2)=  Aik(1,1)*X1 + Aik(1,2)*X2
        RHS_Aij(2,2)=  Aik(2,1)*X1 + Aik(2,2)*X2

      return
      end subroutine ILU1b22


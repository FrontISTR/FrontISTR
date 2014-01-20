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
!C*** module hecmw_solver_BLBiCGSTAB_33
!C***
!
      module hecmw_solver_BLBiCGSTAB_33
      contains
!C
!C*** hecmw_solve_BLBiCGSTAB_33
!C
!C      (1) Block ILU(0) 
!C      (2) Block ILU(1) 
!C      (3) Block ILU(2) 
!C

      subroutine hecmw_solve_BLBiCGSTAB_33                              &
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

      real   (kind=kreal), dimension(2)                :: C0, CG

      integer(kind=kint ):: R, RT, P, PT, S, ST, T, V, MAXIT
      real   (kind=kreal):: TOL
      integer(kind=kint )::i,j,k
      integer(kind=kint ):: ns, nr
      real   (kind=kreal)::S_TIME,E_TIME,START_TIME,END_TIME,SETupTIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2,X3,C20,C2
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::WVAL1,WVAL2,WVAL3
      real   (kind=kreal)::COMMtime,OMEGA

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

      allocate (WW(3*NP,8))
      allocate (WS(3*ns), WR(3*nr))

      WW   = 0.d0
      ERROR= 0

      MAXIT = ITER
      TOL   = RESID

      ERROR= 0

       R= 1
      RT= 2
       P= 3
      PT= 4
       S= 5
      ST= 6
       T= 7
       V= 8

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

      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2,HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) then
        iter = 0
        MAXIT = 0
        RESID = 0
        X = 0.d0
      endif

       END_TIME= HECMW_WTIME()
      SETupTIME= END_TIME - START_TIME
!C===

      START_TIME= HECMW_WTIME()
!C===

!C
!C*************************************************************** iterative procedures
!C
      do iter= 1, MAXIT
!C
!C +-----------------+
!C | RHO= {r}{r_tld} |
!C +-----------------+
!C===
      RHO0= 0.d0
!*voption indep (WW)
      do j= 1, N
        RHO0= RHO0+WW(3*j-2,RT)*WW(3*j-2,R)+WW(3*j-1,RT)*WW(3*j-1,R)    &
     &                                     +WW(3*j  ,RT)*WW(3*j  ,R)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME
!C===

!C
!C +----------------------------------------+
!C | BETA= (RHO/RHO1) * (ALPHA/OMEGA)       |
!C | {p} = {r} + BETA * ( {p} - OMEGA*{v} ) |
!C +----------------------------------------+
!C===
      if (iter.gt.1) then
        BETA= (RHO/RHO1) * (ALPHA/OMEGA)
!*voption indep (WW)
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R)+BETA*(WW(3*j-2,P)-OMEGA*WW(3*j-2,V))
          WW(3*j-1,P)= WW(3*j-1,R)+BETA*(WW(3*j-1,P)-OMEGA*WW(3*j-1,V))
          WW(3*j  ,P)= WW(3*j  ,R)+BETA*(WW(3*j  ,P)-OMEGA*WW(3*j  ,V))
        enddo
       else
!*voption indep (WW)
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R)
          WW(3*j-1,P)= WW(3*j-1,R)
          WW(3*j  ,P)= WW(3*j  ,R)
        enddo
      endif
!C===

!C
!C +--------------------+
!C | {p_tld}= [Minv]{p} |
!C +--------------------+
!C===

!C
!C== Block ILU
!*voption indep (WW)
      do i= 1, N
        WW(3*i-2,PT)= WW(3*i-2,P)
        WW(3*i-1,PT)= WW(3*i-1,P)
        WW(3*i  ,PT)= WW(3*i  ,P)
      enddo

!*voption indep (WW)
      do i= N+1, NP
        WW(3*i-2,PT)= 0.d0
        WW(3*i-1,PT)= 0.d0
        WW(3*i  ,PT)= 0.d0
      enddo

      call hecmw_precond_BILU_33_apply(N, WW(:,PT))

!C==
!C===

!C
!C +-------------------+
!C | {v} = [A] {p_tld} |
!C +-------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,PT) , SOLVER_COMM,    &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

!C
!C-- BEGIN calculation
      do j= 1, N
           X1= WW(3*j-2,PT)
           X2= WW(3*j-1,PT)
           X3= WW(3*j  ,PT)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3

!*voption indep (IAL,AL,WW)
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,PT)
          X2= WW(3*i-1,PT)
          X3= WW(3*i  ,PT)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo

!*voption indep (IAU,AU,WW)
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,PT)
          X2= WW(3*i-1,PT)
          X3= WW(3*i  ,PT)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,V)= WVAL1
        WW(3*j-1,V)= WVAL2
        WW(3*j  ,V)= WVAL3
      enddo
!C===

!C
!C-- calc. ALPHA

      C20= 0.d0
!*voption indep (WW)
      do j= 1, N
        C20= C20 + WW(3*j-2,RT)*WW(3*j-2,V) + WW(3*j-1,RT)*WW(3*j-1,V)  &
     &                                      + WW(3*j  ,RT)*WW(3*j  ,V)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C20, C2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C20, C2, 1, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr) 
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      ALPHA= RHO / C2

!C
!C-- {s}= {r} - ALPHA*{V}
!*voption indep (WW)
      do j= 1, N
        WW(3*j-2,S)= WW(3*j-2,R) - ALPHA*WW(3*j-2,V)
        WW(3*j-1,S)= WW(3*j-1,R) - ALPHA*WW(3*j-1,V)
        WW(3*j  ,S)= WW(3*j  ,R) - ALPHA*WW(3*j  ,V)
      enddo

!C
!C +--------------------+
!C | {s_tld}= [Minv]{s} |
!C +--------------------+
!C===

!C
!C== Block ILU
!*voption indep (WW)
      do i= 1, N
        WW(3*i-2,ST)= WW(3*i-2,S)
        WW(3*i-1,ST)= WW(3*i-1,S)
        WW(3*i  ,ST)= WW(3*i  ,S)
      enddo

!*voption indep (WW)
      do i= N+1, NP
        WW(3*i-2,ST)= 0.d0
        WW(3*i-1,ST)= 0.d0
        WW(3*i  ,ST)= 0.d0
      enddo

      call hecmw_precond_BILU_33_apply(N, WW(:,ST))

!C==
!C===

!C
!C +------------------+
!C | {t} = [A]{s_tld} |
!C +------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,ST) , SOLVER_COMM,    &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

!C
!C-- BEGIN calculation
      do j= 1, N
           X1= WW(3*j-2,ST)
           X2= WW(3*j-1,ST)
           X3= WW(3*j  ,ST)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3

!*voption indep (IAL,AL,WW)
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,ST)
          X2= WW(3*i-1,ST)
          X3= WW(3*i  ,ST)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo

!*voption indep (IAU,AU,WW)
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,ST)
          X2= WW(3*i-1,ST)
          X3= WW(3*i  ,ST)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,T)= WVAL1
        WW(3*j-1,T)= WVAL2
        WW(3*j  ,T)= WVAL3
      enddo
!C===

!C
!C +----------------------------+
!C | OMEGA= ({t}{s}) / ({t}{t}) |
!C +----------------------------+
!C===
      C0(1)= 0.d0
      C0(2)= 0.d0
      CG(1)= 0.d0
      CG(2)= 0.d0

!*voption indep (WW)
      do j= 1, N
        C0(1)= C0(1) + WW(3*j-2,T)*WW(3*j-2,S) + WW(3*j-1,T)*WW(3*j-1,S)&
     &                                         + WW(3*j  ,T)*WW(3*j  ,S)
        C0(2)= C0(2) + WW(3*j-2,T)*WW(3*j-2,T) + WW(3*j-1,T)*WW(3*j-1,T)&
     &                                         + WW(3*j  ,T)*WW(3*j  ,T)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP (C0, CG, 2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C0, CG, 2, MPI_DOUBLE_PRECISION,              &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      OMEGA= CG(1) / CG(2)
!C===

!C
!C +----------------+
!C | update {x},{r} |
!C +----------------+
!C===
      DNRM20= 0.d0
!*voption indep (X,WW)
      do j= 1, N
        X (3*j-2)= X(3*j-2) + ALPHA*WW(3*j-2,PT) + OMEGA*WW(3*j-2,ST)
        X (3*j-1)= X(3*j-1) + ALPHA*WW(3*j-1,PT) + OMEGA*WW(3*j-1,ST)
        X (3*j  )= X(3*j  ) + ALPHA*WW(3*j  ,PT) + OMEGA*WW(3*j  ,ST)
        WW(3*j-2,R)= WW(3*j-2,S) - OMEGA*WW(3*j-2,T)
        WW(3*j-1,R)= WW(3*j-1,S) - OMEGA*WW(3*j-1,T)
        WW(3*j  ,R)= WW(3*j  ,S) - OMEGA*WW(3*j  ,T)
        DNRM20= DNRM20+WW(3*j-2,S)**2+WW(3*j-1,S)**2+WW(3*j,S)**2
      enddo

      RHO1= RHO


      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      RESID= dsqrt(DNRM2/BNRM2)

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

      end subroutine  hecmw_solve_BLBiCGSTAB_33
      end module     hecmw_solver_BLBiCGSTAB_33

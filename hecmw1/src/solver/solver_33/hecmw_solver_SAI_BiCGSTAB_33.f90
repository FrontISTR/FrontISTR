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
!C*** module hecmw_solver_SAI_BiCGSTAB_33
!C***
!C
      module hecmw_solver_SAI_BiCGSTAB_33
      contains
!C
!C*** hecmw_solve_SAI_BiCGSTAB_33
!C
      subroutine hecmw_solve_SAI_BiCGSTAB_33                            &
     &   (N, NP, NPL, AMAT, AMATsai, INDEX, ITEM, INDEXsai, ITEMsai,    &
     &    B, X, RESID, ITER, ERROR, my_rank,                            &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM, Tset, Tsol, Tcomm, ITERlog)

      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_solver_SR_33

      implicit none

      integer(kind=kint ), intent(in):: N, NP, NPL, my_rank
      integer(kind=kint ), intent(in):: NEIBPETOT
      integer(kind=kint ), intent(in):: SOLVER_COMM, ITERlog

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      real(kind=kreal), dimension(3*NP)   , intent(inout):: B, X
      real(kind=kreal), dimension(9*NPL), intent(inout):: AMAT, AMATsai

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INDEX,INDEXsai
      integer(kind=kint ), dimension(  NPL),intent(in) ::  ITEM, ITEMsai

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)

      real(kind=kreal), dimension(:),    allocatable       :: WS, WR
      real(kind=kreal), dimension(:,:),  allocatable       :: WW

      integer(kind=kint ) :: MAXIT

      ! local variables
      real   (kind=kreal):: TOL
      integer(kind=kint )::i,j,k
      real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2,X3,C20,C2
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::WVAL1,WVAL2,WVAL3
      real   (kind=kreal)::OMEGA

      real   (kind=kreal), dimension(2)                :: C0, CG

      integer(kind=kint), parameter :: R = 1
      integer(kind=kint), parameter :: RT= 2
      integer(kind=kint), parameter :: P = 3
      integer(kind=kint), parameter :: PT= 4
      integer(kind=kint), parameter :: S = 5
      integer(kind=kint), parameter :: ST= 6
      integer(kind=kint), parameter :: T = 7
      integer(kind=kint), parameter :: V = 8
      integer(kind=kint), parameter :: PQ= 9

      S_time= HECMW_WTIME()
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0
      allocate (WW(3*NP,9), WS(3*NP), WR(3*NP))

      MAXIT= ITER
      TOL  = RESID
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

!C
!C-- BEGIN calculation

      do j= 1, N
        WVAL1= B(3*j-2)
        WVAL2= B(3*j-1)
        WVAL3= B(3*j  )

        do k= INDEX(j-1)+1, INDEX(j)
          i= ITEM(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AMAT(9*k-8)*X1-AMAT(9*k-7)*X2-AMAT(9*k-6)*X3
          WVAL2= WVAL2 - AMAT(9*k-5)*X1-AMAT(9*k-4)*X2-AMAT(9*k-3)*X3
          WVAL3= WVAL3 - AMAT(9*k-2)*X1-AMAT(9*k-1)*X2-AMAT(9*k  )*X3
        enddo

        WW(3*j-2,R)= WVAL1
        WW(3*j-1,R)= WVAL2
        WW(3*j  ,R)= WVAL3

        WW(3*j-2,RT)= WVAL1
        WW(3*j-1,RT)= WVAL2
        WW(3*j  ,RT)= WVAL3

        WW(3*j-2,PT)= 0.d0
        WW(3*j-1,PT)= 0.d0
        WW(3*j  ,PT)= 0.d0
        WW(3*j-2,ST)= 0.d0
        WW(3*j-1,ST)= 0.d0
        WW(3*j  ,ST)= 0.d0
        WW(3*j-2,PQ)= 0.d0
        WW(3*j-1,PQ)= 0.d0
        WW(3*j  ,PQ)= 0.d0
      enddo

      BNRM20= 0.d0
      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      E_time= HECMW_WTIME()

      Tset= Tset + E_time - S_time

      iter= 0
!C===

      S1_time= HECMW_WTIME()
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
      do j= 1, N
        RHO0= RHO0+WW(3*j-2,RT)*WW(3*j-2,R)+WW(3*j-1,RT)*WW(3*j-1,R)    &
     &                                     +WW(3*j  ,RT)*WW(3*j  ,R)
      enddo

      S_time= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time
!C===

!C
!C +----------------------------------------+
!C | BETA= (RHO/RHO1) * (ALPHA/OMEGA)       |
!C | {p} = {r} + BETA * ( {p} - OMEGA*{v} ) |
!C +----------------------------------------+
!C===
      if (iter.gt.1) then
        BETA= (RHO/RHO1) * (ALPHA/OMEGA)
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R)+BETA*(WW(3*j-2,P)-OMEGA*WW(3*j-2,V))
          WW(3*j-1,P)= WW(3*j-1,R)+BETA*(WW(3*j-1,P)-OMEGA*WW(3*j-1,V))
          WW(3*j  ,P)= WW(3*j  ,R)+BETA*(WW(3*j  ,P)-OMEGA*WW(3*j  ,V))
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
!C +--------------------+
!C | {p_tld}= [Minv]{p} |
!C +--------------------+
!C===

      S_time= HECMW_WTIME()
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,P) , SOLVER_COMM,     &
     &     my_rank)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      do j= 1, N
        WVAL1= 0.d0
        WVAL2= 0.d0
        WVAL3= 0.d0
        do k= INDEXsai(j-1)+1, INDEXsai(j)
           i= ITEMsai(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AMATsai(9*k-8)*X1 + AMATsai(9*k-7)*X2          &
     &                                     + AMATsai(9*k-6)*X3
          WVAL2= WVAL2 + AMATsai(9*k-5)*X1 + AMATsai(9*k-4)*X2          &
     &                                     + AMATsai(9*k-3)*X3
          WVAL3= WVAL3 + AMATsai(9*k-2)*X1 + AMATsai(9*k-1)*X2          &
     &                                     + AMATsai(9*k  )*X3
        enddo

        WW(3*j-2,PT)= WVAL1
        WW(3*j-1,PT)= WVAL2
        WW(3*j  ,PT)= WVAL3
      enddo
!C===

!C
!C +-------------------+
!C | {v} = [A] {p_tld} |
!C +-------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

      S_time= HECMW_WTIME()
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,PT) , SOLVER_COMM,    &
     &     my_rank)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      do j= 1, N
        WVAL1= 0.d0
        WVAL2= 0.d0
        WVAL3= 0.d0
        do k= INDEX(j-1)+1, INDEX(j)
           i= ITEM(k)
          X1= WW(3*i-2,PT)
          X2= WW(3*i-1,PT)
          X3= WW(3*i  ,PT)
          WVAL1= WVAL1 + AMAT(9*k-8)*X1+AMAT(9*k-7)*X2+AMAT(9*k-6)*X3
          WVAL2= WVAL2 + AMAT(9*k-5)*X1+AMAT(9*k-4)*X2+AMAT(9*k-3)*X3
          WVAL3= WVAL3 + AMAT(9*k-2)*X1+AMAT(9*k-1)*X2+AMAT(9*k  )*X3
        enddo

        WW(3*j-2,V)= WVAL1
        WW(3*j-1,V)= WVAL2
        WW(3*j  ,V)= WVAL3
      enddo
!C===

!C
!C-- calc. ALPHA

      C20= 0.d0
      do j= 1, N
        C20= C20 + WW(3*j-2,RT)*WW(3*j-2,V) + WW(3*j-1,RT)*WW(3*j-1,V)  &
     &                                      + WW(3*j  ,RT)*WW(3*j  ,V)
      enddo

      S_time= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C20, C2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C20, C2, 1, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      ALPHA= RHO / C2

!C
!C-- {s}= {r} - ALPHA*{V}
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

      S_time= HECMW_WTIME()
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,S) , SOLVER_COMM,     &
     &     my_rank)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      do j= 1, N
        WVAL1= 0.d0
        WVAL2= 0.d0
        WVAL3= 0.d0
        do k= INDEXsai(j-1)+1, INDEXsai(j)
           i= ITEMsai(k)
          X1= WW(3*i-2,S)
          X2= WW(3*i-1,S)
          X3= WW(3*i  ,S)
          WVAL1= WVAL1 + AMATsai(9*k-8)*X1 + AMATsai(9*k-7)*X2          &
     &                                     + AMATsai(9*k-6)*X3
          WVAL2= WVAL2 + AMATsai(9*k-5)*X1 + AMATsai(9*k-4)*X2          &
     &                                     + AMATsai(9*k-3)*X3
          WVAL3= WVAL3 + AMATsai(9*k-2)*X1 + AMATsai(9*k-1)*X2          &
     &                                     + AMATsai(9*k  )*X3
        enddo

        WW(3*j-2,ST)= WVAL1
        WW(3*j-1,ST)= WVAL2
        WW(3*j  ,ST)= WVAL3
      enddo
!C===

!C
!C +------------------+
!C | {t} = [A]{s_tld} |
!C +------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      S_time= HECMW_WTIME()
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,ST) , SOLVER_COMM,    &
     &     my_rank)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      do j= 1, N
        WVAL1= 0.d0
        WVAL2= 0.d0
        WVAL3= 0.d0
        do k= INDEX(j-1)+1, INDEX(j)
           i= ITEM(k)
          X1= WW(3*i-2,ST)
          X2= WW(3*i-1,ST)
          X3= WW(3*i  ,ST)
          WVAL1= WVAL1 + AMAT(9*k-8)*X1+AMAT(9*k-7)*X2+AMAT(9*k-6)*X3
          WVAL2= WVAL2 + AMAT(9*k-5)*X1+AMAT(9*k-4)*X2+AMAT(9*k-3)*X3
          WVAL3= WVAL3 + AMAT(9*k-2)*X1+AMAT(9*k-1)*X2+AMAT(9*k  )*X3
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

      do j= 1, N
        C0(1)= C0(1) + WW(3*j-2,T)*WW(3*j-2,S) + WW(3*j-1,T)*WW(3*j-1,S)&
     &                                         + WW(3*j  ,T)*WW(3*j  ,S)
        C0(2)= C0(2) + WW(3*j-2,T)*WW(3*j-2,T) + WW(3*j-1,T)*WW(3*j-1,T)&
     &                                         + WW(3*j  ,T)*WW(3*j  ,T)
      enddo

      S_time= HECMW_WTIME()
      call HECMW_allREDUCE_DP (C0, CG, 2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C0, CG, 2, MPI_DOUBLE_PRECISION,              &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      OMEGA= CG(1) / CG(2)
!C===

!C
!C +----------------+
!C | update {x},{r} |
!C +----------------+
!C===
      DNRM20= 0.d0
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

      S_time= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      E_time= HECMW_WTIME()
      Tcomm= Tcomm + E_time - S_time

      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
        if (my_rank.eq.0.and.ITERlog.eq.1) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
!C#####

      if (RESID.le.TOL   ) exit
      if ( ITER.eq.MAXIT ) ERROR= -300
!C===

      enddo
      E1_time= HECMW_WTIME()

!C
!C-- INTERFACE data EXCHANGE
      call hecmw_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM, my_rank)

      deallocate (WW)
      deallocate (WR)
      deallocate (WS)

      Tsol= E1_time - S1_time

      end subroutine hecmw_solve_SAI_BiCGSTAB_33
      end module     hecmw_solver_SAI_BiCGSTAB_33
